#!/usr/bin/env python3
"""Grid-cap trade-off study: recall versus runtime for the period-trial cap.

The default cap (250,000 trials) binds for almost every multi-year
known-answer benchmark light curve, so the grid is coarser than its own
drift rule. Three variants are compared on the same calibration-split hosts
(each host contributes its confirmed planet and its inverted-flux null):

  A  default cap (250k)
  B  cap raised to 1,000,000
  C  TIC density with the conservative duration floor (b <= 0.95, e <= 0.9)
  D  A + 10-minute binning
  E  TIC density with the physical floor (half the central duration)
  F  D + E
  G  F with the cap raised to 1,000,000 (resolution-complete under the physical prior)
  H  F, cap 2,000,000, two-stage search with factor 7 (coarse screen + full refinement)
  I  H with the coarse stage ranked by the fast folding algorithm
  J  candidate default: I plus the physical minimum orbital period

  prepare   sample hosts, fetch TIC densities, write one run manifest per variant
  evaluate  recall at matched null trigger rates (thresholds from each
            variant's own nulls), paired McNemar tests, and CPU per case

Only the calibration split is used, so the blind split stays untouched for
evaluating whatever default is adopted.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
from pathlib import Path

import evaluate_benchmark as ev
from vetting.tier import statistic

csv.field_size_limit(2**31 - 1)
VARIANTS = ("A_cap250k", "B_cap1M", "C_density", "D_bin10", "E_physical", "F_bin10_physical",
            "G_bin10_physical_cap1M", "H_two_stage", "I_ffa_two_stage", "J_default_candidate")
#: Variants that pass the host TIC density to the search.
DENSITY_VARIANTS = ("C_density", "E_physical", "F_bin10_physical", "G_bin10_physical_cap1M",
                    "H_two_stage", "I_ffa_two_stage", "J_default_candidate")
RHO_SUN_G_CM3 = 1.408


def tic_densities(tic_ids: list[str]) -> dict[str, float]:
    """TIC 8.2 stellar density (g/cm^3), from `rho` or mass/radius."""
    from astroquery.mast import Catalogs

    from vetting import network

    network.install_default_timeouts()
    out: dict[str, float] = {}
    for start in range(0, len(tic_ids), 100):
        chunk = [int(t) for t in tic_ids[start:start + 100]]
        table = Catalogs.query_criteria(catalog="Tic", ID=chunk)
        for row in table:
            rho, mass, rad = (row[k] for k in ("rho", "mass", "rad"))
            value = None
            if rho is not None and not hasattr(rho, "mask") and math.isfinite(float(rho)) and rho > 0:
                value = float(rho)
            elif all(v is not None and not hasattr(v, "mask") for v in (mass, rad)) and rad > 0:
                value = float(mass) / float(rad) ** 3
            if value and math.isfinite(value):
                out[str(int(row["ID"]))] = value * RHO_SUN_G_CM3
    return out


def prepare(a) -> None:
    root = a.manifest.resolve().parent
    cases = [c for c in ev.read_csv(a.manifest) if c["split"] == "calibration"]
    by_host: dict[str, dict] = {}
    for case in cases:
        by_host.setdefault(case["host_id"], {})[case["truth_class"]] = case
    hosts = sorted(h for h, v in by_host.items() if {"confirmed_planet", "inverted_null"} <= set(v))
    random.Random(a.seed).shuffle(hosts)
    hosts = sorted(hosts[:a.hosts])
    densities = tic_densities(hosts)
    a.out.mkdir(parents=True, exist_ok=True)
    selected = [by_host[h][k] for h in hosts for k in ("confirmed_planet", "inverted_null")]
    for variant in VARIANTS:
        with (a.out / f"{variant}.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=["case_id", "data_csv", "candidates_csv",
                                                        "stellar_density_g_cm3"])
            writer.writeheader()
            for case in selected:
                density = densities.get(case["host_id"]) if variant in DENSITY_VARIANTS else None
                writer.writerow({"case_id": case["case_id"],
                                 "data_csv": str(root / "lightcurves" / f"{case['case_id']}.csv"),
                                 "candidates_csv": f"{variant}/{case['case_id']}/candidates.csv",
                                 "stellar_density_g_cm3": f"{density:.6g}" if density else ""})
    truth = a.out / "truth.csv"
    with truth.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(selected[0]))
        writer.writeheader()
        writer.writerows(selected)
    print(json.dumps({"hosts": len(hosts), "cases": len(selected),
                      "hosts_with_density": sum(h in densities for h in hosts)}))


def case_result(case: dict, candidates: Path, name: str) -> dict:
    with candidates.open(newline="") as stream:
        rows = [r for r in csv.DictReader(stream) if r.get("class") in ("periodic_dip", "duo_transit")]
    span = float(case.get("baseline_days") or 0) or 1.0
    matched = None
    if case["truth_class"] == "confirmed_planet":
        matched = next((statistic(r, name) for r in rows
                        if any(ev.close_ephemeris(e, case, span, 1.0) for e in ev.ephemerides(r))), None)
    run = json.loads(candidates.with_suffix(".run.json").read_text())
    metrics = run.get("metrics", {})
    return {"best": max((statistic(r, name) for r in rows), default=-math.inf), "matched": matched,
            "cpu": metrics.get("cpu_seconds_sampled"), "wall": metrics.get("wall_seconds"),
            "trials": metrics.get("transit_trials")}


def evaluate(a) -> None:
    from scipy import stats

    truth = {c["case_id"]: c for c in ev.read_csv(a.out / "truth.csv")}
    variants = a.variants or list(VARIANTS)
    unknown = sorted(set(variants) - set(VARIANTS))
    if unknown or "A_cap250k" not in variants:
        raise SystemExit(f"--variants must name known variants including A_cap250k; unknown: {unknown}")
    results = {}
    for variant in variants:
        rows = {}
        for case_id, case in truth.items():
            path = a.out / variant / case_id / "candidates.csv"
            if path.exists():
                rows[case_id] = {**case_result(case, path, a.statistic), "truth": case["truth_class"]}
        results[variant] = rows
    complete = set.intersection(*(set(r) for r in results.values()))
    out = {"statistic": a.statistic, "complete_cases": len(complete), "variants": {}}
    suffix = "" if variants == list(VARIANTS) else "_" + "-".join(v.split("_")[0] for v in variants)
    hits = {}
    for variant, rows in results.items():
        rows = {k: v for k, v in rows.items() if k in complete}
        nulls = sorted((v["best"] for v in rows.values() if v["truth"] == "inverted_null"), reverse=True)
        planets = {k: v for k, v in rows.items() if v["truth"] == "confirmed_planet"}
        summary = {"planets": len(planets), "nulls": len(nulls),
                   "in_candidate_list": sum(v["matched"] is not None for v in planets.values()),
                   "cpu_seconds_total": round(sum(v["cpu"] or 0 for v in rows.values()), 1),
                   "cpu_seconds_median": sorted(v["cpu"] or 0 for v in rows.values())[len(rows) // 2],
                   "trials_median": sorted(v["trials"] or 0 for v in rows.values())[len(rows) // 2]}
        for far in (0.05, 0.10):
            threshold = nulls[min(int(far * len(nulls)), len(nulls) - 1)]
            found = {k for k, v in planets.items() if v["matched"] is not None and v["matched"] > threshold}
            hits[(variant, far)] = found
            summary[f"recall_at_far_{far:.2f}"] = [len(found), len(planets), ev.wilson(len(found), len(planets))]
        out["variants"][variant] = summary
    for variant in (v for v in variants if v != "A_cap250k"):
        for far in (0.05, 0.10):
            gain = len(hits[(variant, far)] - hits[("A_cap250k", far)])
            loss = len(hits[("A_cap250k", far)] - hits[(variant, far)])
            p = float(stats.binomtest(gain, gain + loss, 0.5).pvalue) if gain + loss else 1.0
            out.setdefault("paired_vs_A", {})[f"{variant}@{far:.2f}"] = {"gained": gain, "lost": loss, "mcnemar_p": p}
    (a.out / f"evaluation_{a.statistic}{suffix}.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    prep.add_argument("manifest", type=Path, help="known-answer truth manifest")
    prep.add_argument("out", type=Path)
    prep.add_argument("--hosts", type=int, default=80)
    prep.add_argument("--seed", type=int, default=20261010)
    evaluation = sub.add_parser("evaluate")
    evaluation.add_argument("out", type=Path)
    evaluation.add_argument("--statistic", default="jackknife_score", choices=("score", "jackknife_score"))
    evaluation.add_argument("--variants", nargs="+",
                            help="evaluate only these variants (must include A_cap250k); default all")
    a = p.parse_args()
    prepare(a) if a.command == "prepare" else evaluate(a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
