#!/usr/bin/env python3
"""Paired comparison of two frozen blind benchmark reports.

This compares results already evaluated at thresholds selected on calibration
controls. It does not select thresholds from the blind set.
"""
from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path


def exact_mcnemar_p(first_only: int, second_only: int) -> float:
    discordant = first_only + second_only
    if discordant == 0:
        return 1.0
    tail = sum(math.comb(discordant, k)
               for k in range(min(first_only, second_only) + 1)) / (2 ** discordant)
    return min(1.0, 2.0 * tail)


def paired_summary(first: dict[str, bool], second: dict[str, bool],
                   bootstrap_replicates: int = 50_000,
                   seed: int = 20261010) -> dict:
    if first.keys() != second.keys() or not first:
        raise ValueError("paired outcomes must have the same nonempty case IDs")
    ids = sorted(first)
    differences = [int(first[key]) - int(second[key]) for key in ids]
    first_only = sum(a and not second[key] for key, a in first.items())
    second_only = sum(b and not first[key] for key, b in second.items())
    both = sum(first[key] and second[key] for key in ids)
    neither = len(ids) - first_only - second_only - both

    rng = random.Random(seed)
    n = len(ids)
    samples = []
    for _ in range(bootstrap_replicates):
        samples.append(sum(differences[rng.randrange(n)] for _ in range(n)) / n)
    samples.sort()
    lo = samples[int(0.025 * bootstrap_replicates)]
    hi = samples[min(bootstrap_replicates - 1,
                     math.ceil(0.975 * bootstrap_replicates) - 1)]
    return {
        "trials": n,
        "first_rate": sum(first.values()) / n,
        "second_rate": sum(second.values()) / n,
        "difference_first_minus_second": sum(differences) / n,
        "paired_bootstrap_95_percentile_interval": [lo, hi],
        "first_only": first_only,
        "second_only": second_only,
        "both": both,
        "neither": neither,
        "mcnemar_exact_two_sided_p": exact_mcnemar_p(first_only, second_only),
        "bootstrap_replicates": bootstrap_replicates,
        "bootstrap_seed": seed,
    }


def compare_reports(first_report: dict, second_report: dict,
                    bootstrap_replicates: int = 50_000,
                    seed: int = 20261010) -> dict:
    first_rows = {row["case_id"]: row
                  for row in first_report["evaluation"]["case_outcomes"]}
    second_rows = {row["case_id"]: row
                   for row in second_report["evaluation"]["case_outcomes"]}
    if first_rows.keys() != second_rows.keys():
        raise ValueError("reports do not contain the same blind cases")
    by_class = {}
    for truth, metric in (("confirmed_planet", "planet_recovered"),
                          ("inverted_null", "any_trigger"), ("null", "any_trigger")):
        ids = [case_id for case_id, row in first_rows.items()
               if row["truth_class"] == truth]
        if not ids:
            continue
        by_class[truth] = paired_summary(
            {case_id: bool(first_rows[case_id][metric]) for case_id in ids},
            {case_id: bool(second_rows[case_id][metric]) for case_id in ids},
            bootstrap_replicates, seed,
        )
    return {
        "first_statistic": first_report.get("statistic", "unspecified"),
        "second_statistic": second_report.get("statistic", "unspecified"),
        "thresholds": {
            "first": first_report.get("calibrated_threshold"),
            "second": second_report.get("calibrated_threshold"),
        },
        "comparison": "paired host-level bootstrap and exact McNemar test; thresholds come from calibration data",
        "blind_results": by_class,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("first", type=Path)
    parser.add_argument("second", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap-replicates", type=int, default=50_000)
    parser.add_argument("--seed", type=int, default=20261010)
    args = parser.parse_args()
    if args.bootstrap_replicates < 1:
        parser.error("bootstrap replicates must be positive")
    try:
        result = compare_reports(json.loads(args.first.read_text()),
                                 json.loads(args.second.read_text()),
                                 args.bootstrap_replicates, args.seed)
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        parser.error(str(error))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
