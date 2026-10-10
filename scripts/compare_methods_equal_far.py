#!/usr/bin/env python3
"""astrosearch (raw score, jackknife score) vs astropy BLS at equal blind-split
null trigger rates, on the known-answer benchmark (Experiment F follow-up).

Reads the per-case BLS cache written by baseline_bls_benchmark.py. For each
target false-alarm rate, each method's threshold is the blind null quantile,
and recall is the fraction of blind planets matched at or above it. This is
an evaluation of ranking quality (an ROC point), not a deployable threshold.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import baseline_bls_benchmark as bls  # noqa: E402
import evaluate_benchmark as ev  # noqa: E402
from vetting.tier import statistic  # noqa: E402


def astrosearch_rows(cases, root, name):
    rows = []
    for case in cases:
        transit, _ = ev.read_scores(case, root, 0.0)
        span = float(case.get("baseline_days") or 0) or 1.0
        matched = None
        if case["truth_class"] == "confirmed_planet":
            matched = next((statistic(r, name) for r in transit
                            if any(ev.close_ephemeris(e, case, span, 1.0) for e in ev.ephemerides(r))), None)
        rows.append({"case_id": case["case_id"], "truth": case["truth_class"],
                     "best": max((statistic(r, name) for r in transit), default=-math.inf),
                     "matched": matched})
    return rows


def recall_at(rows, far):
    nulls = sorted((r["best"] for r in rows if r["truth"] == "inverted_null"), reverse=True)
    threshold = nulls[int(far * len(nulls))]
    planets = [r for r in rows if r["truth"] == "confirmed_planet"]
    hits = {r["case_id"] for r in planets if r["matched"] is not None and r["matched"] > threshold}
    return hits, len(planets)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("manifest", type=Path)
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    root = a.manifest.parent
    cases = [c for c in ev.read_csv(a.manifest)
             if (root / c["candidates_csv"]).exists() and c["split"] == "blind"
             and (a.cache / f"{c['case_id']}.json").exists()]
    methods = {"astrosearch_score": astrosearch_rows(cases, root, "score"),
               "astrosearch_jackknife": astrosearch_rows(cases, root, "jackknife_score"),
               "bls_sde": [bls.run_case((c, str(root), str(a.cache))) for c in cases]}
    for row in methods["bls_sde"]:
        row["truth"] = row.pop("truth")
    out = {"blind_cases": len(cases), "points": []}
    from scipy import stats
    for far in (0.02, 0.04, 0.05, 0.07, 0.10):
        point = {"far": far}
        hits = {}
        for name, rows in methods.items():
            hits[name], total = recall_at(rows, far)
            point[name] = [len(hits[name]), total]
        only_j = len(hits["astrosearch_jackknife"] - hits["bls_sde"])
        only_b = len(hits["bls_sde"] - hits["astrosearch_jackknife"])
        point["jackknife_vs_bls_mcnemar"] = {
            "jackknife_only": only_j, "bls_only": only_b,
            "p": float(stats.binomtest(only_j, only_j + only_b, 0.5).pvalue) if only_j + only_b else 1.0}
        out["points"].append(point)
        print(json.dumps(point))
    a.output.write_text(json.dumps(out, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
