#!/usr/bin/env python3
"""Measure how often each signal diagnostic fires on known planets, inverted
(noise-only) controls and real eclipsing binaries.

Every case is searched with production Profile J and five passes, then
vetting/signal_diagnostics.py runs on the output. The diagnostics are
report-only, so this measures flag rates, not detection performance; no
threshold is changed and no candidate is removed.

For planets the flags are read on the candidate that matches the catalogue
ephemeris (if any); for controls and binaries on the highest-jackknife
transit candidate. Only the calibration split is accepted, so the blind split
stays untouched. Outputs go to --output only.
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import evaluate_benchmark as ev
from vetting import signal_diagnostics as sd
from vetting.io import read_table
from vetting.tier import jackknife_score

csv.field_size_limit(2**31 - 1)
PROFILE_J = ["--search-profile", "fast", "--transit-search-passes", "5"]
FLAGS = ("detection_tier", "depth_inconsistent", "masking_dependent", "modulation_explained",
         "interleaved_secondary", "one_parity_events")


def search(binary: str, curve: Path, out: Path, threads: int, timeout: float) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    candidates = out / "candidates.csv"
    if not candidates.exists():
        # Write beside the target and rename on success, so an interrupted
        # run never leaves a partial file that a resumed run would reuse.
        partial = out / "candidates.partial.csv"
        partial.unlink(missing_ok=True)
        subprocess.run([binary, "--input", str(curve), "--output", str(partial), "--threads", str(threads),
                        *PROFILE_J], check=True, capture_output=True, text=True, timeout=timeout)
        partial.replace(candidates)
    return candidates


def flags_for(report: dict, label: str, tier: str) -> dict:
    out = {flag: False for flag in FLAGS}
    out["detection_tier"] = tier == "detection"
    for group in report["event_groups"]:
        members = {m["detection"]: m for m in group["members"]}
        if label not in members:
            continue
        issues = group["depth_check"]["issues"]
        out["depth_inconsistent"] = any(i.startswith(f"{label}:") and "box depth" in i and "local depth" in i
                                        for i in issues) or any(i.startswith("reported depths") for i in issues)
        out["masking_dependent"] = any(i.startswith(f"{label}:") and "unmasked light curve" in i for i in issues)
        out["one_parity_events"] = members[label]["parity_coverage"].endswith("(effective period 2P)")
    for conflict in report["variability_conflicts"]:
        if conflict["detection"] == label:
            verdict = conflict.get("modulation_test", {}).get("verdict", "")
            out["modulation_explained"] = verdict.startswith("competing explanation")
    out["interleaved_secondary"] = any(label in pair["detections"] and abs(pair["phase_of_second_events"] - 0.5)
                                       <= sd.SECONDARY_PHASE_WINDOW for pair in report["interleaved_pairs"])
    return out


def one_case(job: tuple) -> dict:
    kind, case_id, curve, truth, binary, out_root, threads, timeout = job
    out = Path(out_root) / kind / case_id
    result = {"class": kind, "case_id": case_id}
    try:
        candidates = search(binary, Path(curve), out, threads, timeout)
    except subprocess.TimeoutExpired:
        return {**result, "status": "search timeout"}
    except subprocess.CalledProcessError as error:
        return {**result, "status": f"search failed: {error.stderr.strip()[-200:]}"}
    rows = ev.read_csv(candidates)
    transits = [r for r in rows if r["class"] in ("periodic_dip", "duo_transit")]
    if truth:
        span = 30.0
        chosen = next((r for r in transits
                       if any(ev.close_ephemeris(e, truth, span, 1.0) for e in ev.ephemerides(r))), None)
    else:
        chosen = max(transits, key=jackknife_score, default=None)
    if chosen is None:
        return {**result, "status": "no matching transit candidate" if truth else "no transit candidate"}
    data = read_table(Path(curve))
    report = sd.diagnose(rows, None, data[:, 0], data[:, 1])
    (out / "diagnostics.json").write_text(json.dumps(report, indent=1) + "\n")
    label = f"astrosearch:{chosen['rank']}"
    return {**result, "status": "evaluated", "rank": chosen["rank"], "period_days": chosen["period_days"],
            "jackknife_score": round(jackknife_score(chosen), 3), **flags_for(report, label, chosen.get("tier", ""))}


def summarize(rows: list[dict]) -> dict:
    out = {}
    for kind in sorted({r["class"] for r in rows}):
        subset = [r for r in rows if r["class"] == kind]
        evaluated = [r for r in subset if r["status"] == "evaluated"]
        status = {}
        for r in subset:
            status[r["status"]] = status.get(r["status"], 0) + 1
        out[kind] = {"cases": len(subset), "evaluated": len(evaluated), "status": status}
        for flag in FLAGS:
            hits = sum(bool(r[flag]) for r in evaluated)
            out[kind][flag] = [hits, len(evaluated), ev.wilson(hits, len(evaluated))]
            tiered = [r for r in evaluated if r["detection_tier"]]
            if flag != "detection_tier":
                hits_t = sum(bool(r[flag]) for r in tiered)
                out[kind][f"{flag}_among_detection_tier"] = [hits_t, len(tiered), ev.wilson(hits_t, len(tiered))]
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--known-answer", type=Path, required=True, help="known-answer truth-manifest.csv")
    p.add_argument("--split", default="calibration", choices=("calibration",),
                   help="only the calibration split is allowed; the blind split stays untouched")
    p.add_argument("--eb-survey", type=Path, help="hard-negative survey folder with */combined.csv")
    p.add_argument("--binary", required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--timeout", type=float, default=3600)
    p.add_argument("--limit", type=int, default=0, help="cases per class (0: all)")
    a = p.parse_args(argv)
    root = a.known_answer.resolve().parent
    jobs = []
    cases = [c for c in ev.read_csv(a.known_answer) if c["split"] == a.split]
    for kind, truth_class in (("planet", "confirmed_planet"), ("control", "inverted_null")):
        chosen = [c for c in cases if c["truth_class"] == truth_class]
        chosen = chosen[:a.limit] if a.limit else chosen
        jobs += [(kind, c["case_id"], str(root / "lightcurves" / f"{c['case_id']}.csv"),
                  c if kind == "planet" else None, a.binary, str(a.output), a.threads, a.timeout) for c in chosen]
    if a.eb_survey:
        curves = sorted(x for x in a.eb_survey.rglob("combined.csv") if "mast-cache" not in x.parts)
        curves = curves[:a.limit] if a.limit else curves
        jobs += [("eclipsing_binary", f"{c.parent.parent.name}-{c.parent.name}", str(c), None, a.binary,
                  str(a.output), a.threads, a.timeout) for c in curves]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        rows = list(pool.map(one_case, jobs, chunksize=1))
    a.output.mkdir(parents=True, exist_ok=True)
    fields = ["class", "case_id", "status", "rank", "period_days", "jackknife_score", *FLAGS]
    with (a.output / "cases.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    summary = {"search": PROFILE_J, "split": a.split, "binary": a.binary, "classes": summarize(rows)}
    (a.output / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary["classes"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
