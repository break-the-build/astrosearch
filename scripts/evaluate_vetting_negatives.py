#!/usr/bin/env python3
"""How often do astrophysical false positives survive light-curve vetting?

The search finds eclipses by design (every cataloged eclipsing binary in the
partial hard-negative run exceeded the detection threshold), so the useful
question is whether the gate's light-curve vetting rejects them while
keeping real planets. Each case's candidate goes through the same checks the
gate applies before catalog, centroid and FPP stages:

  * secondary-eclipse scan (vetting.diagnostics.secondary_scan)
  * odd/even depth consistency
  * grazing / V-shape from the gate's batman fit (ctoi_gate.physical_fit,
    local covariance)
  * calibrated detection tier (jackknife score)

Classes:
  planet      known-answer confirmed planets, the candidate matching the truth
  eb          TESS Ten Thousand Catalog eclipsing binaries (Kostov et al. 2025),
              the top transit candidate of the survey search
  systematic  known-answer inverted-flux controls, top transit candidate
  blend       synthetic: a background eclipsing binary diluted to planet-like
              depth (secondary eclipse, V-shaped eclipses) injected into real
              light curves with their planet removed; vetted at the injected
              ephemeris

Light-curve vetting cannot reject a blend whose shape is planet-like and
whose secondary is undetectable; centroid and neighbor checks (pixel data,
TIC) exist for those and are not exercised here.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import evaluate_benchmark as ev
from ctoi_gate import physical_fit
from vetting import diagnostics
from vetting.io import Candidate, phase_offset, read_lightcurve, read_table
from vetting.tier import jackknife_score

csv.field_size_limit(2**31 - 1)
TIER = json.loads((Path(__file__).resolve().parents[1] / "data/calibration/detection_threshold.json")
                  .read_text())["threshold"]


def to_candidate(row: dict) -> Candidate:
    return Candidate(rank=int(row["rank"]), kind=row["class"], score=float(row["score"]),
                     period=float(row["period_days"]), epoch=float(row["epoch"]),
                     duration=float(row["duration_days"]), depth_ppm=float(row["depth_ppm"]),
                     n_events=int(row["n_events"]), event_epochs=(), row=row)


def vet(lightcurve: Path, candidates: Path, row: dict, work: Path) -> dict:
    lc = read_lightcurve(lightcurve)
    cand = to_candidate(row)
    fit = physical_fit(lightcurve, candidates, cand, work, mcmc_steps=0)
    period = fit.get("period") or cand.period
    epoch = fit.get("epoch") or cand.epoch
    duration = fit.get("duration") or cand.duration
    events = diagnostics.event_depths(lc, epoch, period, duration)
    checks = {"secondary": diagnostics.secondary_scan(lc, epoch, period, duration)["status"],
              "odd_even": diagnostics.odd_even(events)["status"],
              "grazing": diagnostics.grazing(fit)["status"] if "error" not in fit else "NOT_EVALUATED"}
    rejected_by = [k for k, v in checks.items() if v == "FAIL"]
    return {**checks, "rejected_by": rejected_by, "lc_vetting_pass": not rejected_by,
            "detection_tier": jackknife_score(row) >= TIER, "period": period,
            "impact": fit.get("impact"), "rp_rs": fit.get("rp_rs")}


def transit_rows(path: Path) -> list[dict]:
    with path.open(newline="") as stream:
        return [r for r in csv.DictReader(stream) if r.get("class") in ("periodic_dip", "duo_transit")]


def job_planet_or_null(args) -> dict:
    case, root = args
    path = Path(root) / case["candidates_csv"]
    rows = transit_rows(path)
    if case["truth_class"] == "confirmed_planet":
        span = float(case.get("baseline_days") or 0) or 1.0
        row = next((r for r in rows if any(ev.close_ephemeris(e, case, span, 1.0)
                                           for e in ev.ephemerides(r))), None)
        label = "planet"
    else:
        row = max(rows, key=lambda r: float(r["score"]), default=None)
        label = "systematic"
    base = {"case": case["case_id"], "class": label}
    if row is None:
        return {**base, "status": "no candidate"}
    with tempfile.TemporaryDirectory() as tmp:
        try:
            return {**base, "status": "vetted",
                    **vet(Path(root) / "lightcurves" / f"{case['case_id']}.csv", path, row, Path(tmp))}
        except Exception as error:
            return {**base, "status": f"error: {type(error).__name__}: {error}"}


def job_eb(args) -> dict:
    case_dir, period_catalog = args
    case_dir = Path(case_dir)
    base = {"case": case_dir.parent.name, "class": "eb", "catalog_period": period_catalog}
    rows = transit_rows(case_dir / "candidates.csv")
    if not rows:
        return {**base, "status": "no candidate"}
    row = max(rows, key=lambda r: float(r["score"]))
    with tempfile.TemporaryDirectory() as tmp:
        try:
            result = vet(case_dir / "combined.csv", case_dir / "candidates.csv", row, Path(tmp))
        except Exception as error:
            return {**base, "status": f"error: {type(error).__name__}: {error}"}
    ratio = result["period"] / period_catalog if period_catalog else None
    return {**base, "status": "vetted", **result, "period_ratio_to_catalog": ratio}


def trapezoid(time, period, epoch, duration, depth, ingress_fraction):
    """Symmetric trapezoid; ingress_fraction = T12 / T14 (0.5 is a V)."""
    x = np.abs(phase_offset(time, epoch, period)) / (duration / 2)
    flat = 1 - 2 * ingress_fraction
    shape = np.clip((1 - x) / max(1 - flat, 1e-9), 0, 1)
    return np.where(x <= flat, depth, depth * shape)


def job_blend(args) -> dict:
    case, root, seed = args
    rng = random.Random(seed)
    data = read_table(Path(root) / "lightcurves" / f"{case['case_id']}.csv")
    known = np.abs(phase_offset(data[:, 0], float(case["epoch_days"]), float(case["period_days"])))
    data = data[known > float(case["duration_days"])]
    time, flux, error = data[:, 0], data[:, 1], data[:, 2]
    period = math.exp(rng.uniform(math.log(0.8), math.log(10)))
    epoch = float(time[0] + rng.uniform(0, period))
    duration = min(0.3, 0.06 * period ** (1 / 3) * rng.uniform(1.0, 2.0))
    inside = np.abs(phase_offset(time, epoch, period)) <= duration / 2
    scatter = 1.4826 * float(np.median(np.abs(np.diff(flux)))) / math.sqrt(2)
    depth = rng.uniform(15, 40) * scatter / math.sqrt(max(inside.sum(), 1))
    secondary_ratio = rng.uniform(0.1, 0.6)
    ingress = rng.uniform(0.3, 0.5)
    model = (trapezoid(time, period, epoch, duration, depth, ingress)
             + trapezoid(time, period, epoch + period / 2, duration, depth * secondary_ratio, ingress))
    injected = flux * (1 - model)
    base = {"case": case["case_id"], "class": "blend", "depth_ppm": depth * 1e6,
            "secondary_ratio": secondary_ratio, "ingress_fraction": ingress}
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        curve = tmp / "lc.csv"
        np.savetxt(curve, np.column_stack([time, injected, error]), delimiter=",",
                   header="time,flux,flux_err", comments="", fmt="%.8f")
        row = {"rank": "1", "class": "periodic_dip", "score": "50", "period_days": f"{period:.8f}",
               "epoch": f"{epoch:.8f}", "duration_days": f"{duration:.8f}",
               "depth_ppm": f"{depth * 1e6:.3f}", "n_events": "10", "class_rank": "1",
               "event_evidence": "", "alias_ephemerides": ""}
        candidates = tmp / "candidates.csv"
        with candidates.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(row))
            writer.writeheader()
            writer.writerow(row)
        try:
            result = vet(curve, candidates, row, tmp)
        except Exception as error:
            return {**base, "status": f"error: {type(error).__name__}: {error}"}
    result.pop("detection_tier")  # not a search output; the tier is not meaningful here
    return {**base, "status": "vetted", **result}


def summarize(rows: list[dict]) -> dict:
    out = {}
    for label in ("planet", "eb", "systematic", "blend"):
        sub = [r for r in rows if r["class"] == label]
        vetted = [r for r in sub if r.get("status") == "vetted"]
        passed = [r for r in vetted if r["lc_vetting_pass"]]
        block = {"cases": len(sub), "vetted": len(vetted),
                 "lc_vetting_pass": [len(passed), len(vetted), ev.wilson(len(passed), len(vetted))],
                 "rejected_by": {k: sum(k in r["rejected_by"] for r in vetted)
                                 for k in ("secondary", "odd_even", "grazing")}}
        if label != "blend":
            both = [r for r in passed if r.get("detection_tier")]
            block["tier_and_lc_vetting_pass"] = [len(both), len(vetted), ev.wilson(len(both), len(vetted))]
        out[label] = block
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--known-answer", type=Path, required=True, help="known-answer truth manifest")
    p.add_argument("--eb-survey", type=Path, help="hard-negative survey output folder (omit to skip EBs)")
    p.add_argument("--eb-targets", type=Path, help="probed EB targets (catalog periods)")
    p.add_argument("--split", default="blind")
    p.add_argument("--planets", type=int, default=120)
    p.add_argument("--blends", type=int, default=80)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--seed", type=int, default=20261010)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    root = a.known_answer.resolve().parent
    cases = [c for c in ev.read_csv(a.known_answer)
             if c["split"] == a.split and (root / c["candidates_csv"]).exists()]
    rng = random.Random(a.seed)
    planets = [c for c in cases if c["truth_class"] == "confirmed_planet"]
    nulls = [c for c in cases if c["truth_class"] == "inverted_null"]
    rng.shuffle(planets)
    rng.shuffle(nulls)
    planets, nulls = planets[:a.planets], nulls[:a.planets]
    hosts = rng.sample(planets, min(a.blends, len(planets)))
    periods = ({r["tic_id"]: float(r["period_days"]) for r in ev.read_csv(a.eb_targets) if r.get("period_days")}
               if a.eb_targets else {})
    eb_dirs = sorted(p.parent for p in a.eb_survey.glob("*/TIC*/*/candidates.csv")) if a.eb_survey else []
    jobs = ([(job_planet_or_null, (c, str(root))) for c in planets + nulls]
            + [(job_eb, (str(d), periods.get(d.parent.name[3:]))) for d in eb_dirs]
            + [(job_blend, (c, str(root), a.seed + i)) for i, c in enumerate(hosts)])
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futures = [pool.submit(fn, arg) for fn, arg in jobs]
        rows = [f.result() for f in futures]
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with a.output.with_suffix(".csv").open("w", newline="") as stream:
        fields = sorted({k for r in rows for k in r})
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    summary = {"tier_threshold": TIER, "split": a.split, **summarize(rows)}
    a.output.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
