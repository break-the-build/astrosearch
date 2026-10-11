#!/usr/bin/env python3
"""Check whether target-yield ranks predict blind known-answer recoveries."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import calibrate_detection_threshold as calibration
import evaluate_benchmark as ev


def read_rows(path: Path) -> list[dict]:
    return ev.read_csv(path)


def rankdata(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        mean_rank = (start + 1 + end) / 2
        for position in range(start, end):
            ranks[order[position]] = mean_rank
        start = end
    return ranks


def spearman(a: list[float], b: list[float]) -> float | None:
    if len(a) < 2:
        return None
    ra, rb = rankdata(a), rankdata(b)
    ma, mb = sum(ra) / len(ra), sum(rb) / len(rb)
    cov = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    va = sum((x - ma) ** 2 for x in ra)
    vb = sum((y - mb) ** 2 for y in rb)
    return cov / (va * vb) ** 0.5 if va and vb else None


def wilson(k: int, n: int) -> list[float] | None:
    if not n:
        return None
    z, p = 1.959963984540054, k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / d
    return [c - h, c + h]


def yield_map(path: Path) -> dict[str, float]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return {str(row["tic_id"]): float(row["expected_detectable_planets"])
                for row in csv.DictReader(stream)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True,
                        help="case-metadata.csv containing host_id and TIC stellar fields")
    parser.add_argument("--ranking", type=Path, required=True,
                        help="astrosearch-target-yield output for the same hosts")
    parser.add_argument("--calibration", type=Path, required=True,
                        help="calibrate_detection_threshold.py output for this profile")
    parser.add_argument("--alternative-ranking", type=Path,
                        help="optional ranking under an alternative discovery-prior assumption")
    parser.add_argument("--bins", type=int, default=10)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    threshold_record = json.loads(args.calibration.read_text())
    threshold = float(threshold_record["evaluation"]["threshold_from_calibration_split"])
    rows = read_rows(args.manifest)
    candidates_root = args.manifest.parent
    metadata = {r["case_id"]: r for r in read_rows(args.metadata)}
    yields = yield_map(args.ranking)
    blind = []
    missing = []
    for case in rows:
        if case["split"] != "blind" or case["truth_class"] != "confirmed_planet":
            continue
        meta = metadata.get(case["case_id"])
        if not meta or str(meta.get("host_id")) not in yields:
            missing.append(case["case_id"])
            continue
        _, matched = calibration.case_values(case, candidates_root, "jackknife_score")
        blind.append({"case_id": case["case_id"], "tic_id": str(meta["host_id"]),
                      "yield": yields[str(meta["host_id"])],
                      "recovered": matched is not None and matched >= threshold})
    if not blind:
        parser.error("no blind positive cases joined to target yields")
    ordered = sorted(blind, key=lambda r: (r["yield"], r["tic_id"]))
    bins = min(args.bins, len(ordered))
    grouped = []
    for i in range(bins):
        subset = ordered[i * len(ordered) // bins:(i + 1) * len(ordered) // bins]
        recovered = sum(r["recovered"] for r in subset)
        grouped.append({"bin": i + 1, "n": len(subset),
                        "yield_min": subset[0]["yield"], "yield_max": subset[-1]["yield"],
                        "yield_median": sorted(r["yield"] for r in subset)[len(subset) // 2],
                        "recovered": recovered, "recall": recovered / len(subset),
                        "recall_wilson_95": wilson(recovered, len(subset))})
    recovered_yields = [r["yield"] for r in blind if r["recovered"]]
    missed_yields = [r["yield"] for r in blind if not r["recovered"]]
    auc = (sum(x > y for x in recovered_yields for y in missed_yields)
           + 0.5 * sum(x == y for x in recovered_yields for y in missed_yields)) / (
               len(recovered_yields) * len(missed_yields)) if recovered_yields and missed_yields else None

    report = {"threshold": threshold, "calibration_profile": threshold_record.get("context", {}),
              "blind_positive_cases_joined": len(blind), "missing_case_ids": missing,
              "recovered": sum(r["recovered"] for r in blind),
              "binned_by_expected_detectable_yield": grouped,
              "recovery_yield_auc": auc,
              "spearman_yield_vs_recovery": spearman(
                  [r["yield"] for r in blind], [float(r["recovered"]) for r in blind])}
    if args.alternative_ranking:
        alt = yield_map(args.alternative_ranking)
        ids = sorted(yields.keys() & alt.keys())
        k = max(1, len(ids) // 10)
        base_top = set(sorted(ids, key=lambda tic: yields[tic], reverse=True)[:k])
        alt_top = set(sorted(ids, key=lambda tic: alt[tic], reverse=True)[:k])
        report["prior_sensitivity"] = {
            "common_targets": len(ids), "top_decile_size": k,
            "spearman_rank_correlation": spearman([yields[t] for t in ids], [alt[t] for t in ids]),
            "top_decile_overlap_count": len(base_top & alt_top),
            "top_decile_overlap_fraction_of_baseline": len(base_top & alt_top) / k,
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
