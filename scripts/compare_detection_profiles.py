#!/usr/bin/env python3
"""Compare two searches on paired known-answer cases at their calibrated tiers.

Thresholds must come from calibrate_detection_threshold.py outputs. They are
selected on calibration controls; this script reports blind performance and
paired gains/losses without using blind labels to set a threshold.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import calibrate_detection_threshold as calibration
import evaluate_benchmark as ev


def wilson(successes: int, trials: int) -> list[float] | None:
    if not trials:
        return None
    z = 1.959963984540054
    p = successes / trials
    denom = 1 + z * z / trials
    center = (p + z * z / (2 * trials)) / denom
    half = z * ((p * (1 - p) / trials + z * z / (4 * trials * trials)) ** 0.5) / denom
    return [center - half, center + half]


def lookup(manifest: Path, threshold_path: Path, statistic: str) -> tuple[dict, dict]:
    config = json.loads(threshold_path.read_text())
    threshold = float(config["evaluation"]["threshold_from_calibration_split"])
    min_events = int(config.get("event_structure_min_events", 0))
    rows = ev.read_csv(manifest)
    root = manifest.parent
    results = {}
    for case in rows:
        if case["split"] == "calibration" and case["truth_class"] == "confirmed_planet":
            continue
        best, matched = calibration.case_values(case, root, statistic, min_events)
        sidecar = (root / case["candidates_csv"]).with_name("candidates.run.json")
        metrics = json.loads(sidecar.read_text()).get("metrics", {}) if sidecar.is_file() else {}
        results[case["case_id"]] = {
            "split": case["split"], "truth": case["truth_class"],
            "trigger": best >= threshold,
            "recovered": matched is not None and matched >= threshold,
            "best": best, "matched": matched,
            "cpu_seconds": metrics.get("cpu_seconds_sampled"),
            "wall_seconds": metrics.get("wall_seconds"),
        }
    return {"threshold": threshold, "statistic": statistic,
            "min_events": min_events}, results


def paired(a: dict[str, dict], b: dict[str, dict], split: str, truth: str,
           field: str) -> dict:
    keys = sorted(k for k in a.keys() & b.keys()
                  if a[k]["split"] == split and b[k]["split"] == split
                  and a[k]["truth"] == truth and b[k]["truth"] == truth)
    a_only = sum(a[k][field] and not b[k][field] for k in keys)
    b_only = sum(b[k][field] and not a[k][field] for k in keys)
    both = sum(a[k][field] and b[k][field] for k in keys)
    neither = len(keys) - a_only - b_only - both
    return {"n": len(keys), "profile_a_only": a_only, "profile_b_only": b_only,
            "both": both, "neither": neither,
            "paired_difference_b_minus_a": (b_only - a_only) / len(keys) if keys else None,
            "mcnemar_exact_two_sided_p": mcnemar_exact(a_only, b_only),
            "discordant_case_ids_a_only": [k for k in keys if a[k][field] and not b[k][field]],
            "discordant_case_ids_b_only": [k for k in keys if b[k][field] and not a[k][field]]}


def mcnemar_exact(a_only: int, b_only: int) -> float | None:
    discordant = a_only + b_only
    if not discordant:
        return None
    k = min(a_only, b_only)
    tail = sum(math.comb(discordant, i) for i in range(k + 1)) / (2 ** discordant)
    return min(1.0, 2 * tail)


def metrics(result: dict[str, dict], split: str, truth: str, field: str) -> dict:
    rows = [r for r in result.values() if r["split"] == split and r["truth"] == truth]
    successes = sum(r[field] for r in rows)
    return {"successes": successes, "trials": len(rows),
            "rate": successes / len(rows) if rows else None,
            "wilson_95": wilson(successes, len(rows))}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-a-manifest", type=Path, required=True)
    parser.add_argument("--profile-a-calibration", type=Path, required=True)
    parser.add_argument("--profile-a-name", default="standard")
    parser.add_argument("--profile-b-manifest", type=Path, required=True)
    parser.add_argument("--profile-b-calibration", type=Path, required=True)
    parser.add_argument("--profile-b-name", default="profile-j")
    parser.add_argument("--statistic", default="jackknife_score")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    a_context, a = lookup(args.profile_a_manifest, args.profile_a_calibration, args.statistic)
    b_context, b = lookup(args.profile_b_manifest, args.profile_b_calibration, args.statistic)
    common = sorted(a.keys() & b.keys())
    if not common:
        parser.error("the manifests have no case IDs in common")
    resource_pairs = [key for key in common if a[key].get("cpu_seconds") is not None
                      and b[key].get("cpu_seconds") is not None]
    cpu_a = sum(a[key]["cpu_seconds"] for key in resource_pairs)
    cpu_b = sum(b[key]["cpu_seconds"] for key in resource_pairs)
    wall_pairs = [key for key in common if a[key].get("wall_seconds") is not None
                  and b[key].get("wall_seconds") is not None]
    resource_totals = {}
    for name, values in ((args.profile_a_name, a), (args.profile_b_name, b)):
        cpu_rows = [r["cpu_seconds"] for r in values.values() if r.get("cpu_seconds") is not None]
        wall_rows = [r["wall_seconds"] for r in values.values() if r.get("wall_seconds") is not None]
        resource_totals[name] = {
            "cpu_measured_cases": len(cpu_rows), "cpu_seconds_sum": sum(cpu_rows),
            "per_case_wall_measured_cases": len(wall_rows),
            "summed_per_case_wall_seconds": sum(wall_rows)}
    report = {"profiles": {args.profile_a_name: a_context, args.profile_b_name: b_context},
              "paired_common_case_ids": len(common), "statistic": args.statistic,
              "resource_totals_by_profile": resource_totals,
              "paired_resources": {
                  "cpu_measured_cases": len(resource_pairs),
                  "cpu_seconds_sum": {args.profile_a_name: cpu_a, args.profile_b_name: cpu_b},
                  "cpu_speedup_a_over_b": cpu_a / cpu_b if cpu_b else None,
                  "summed_per_case_wall_seconds": {
                      args.profile_a_name: sum(a[k]["wall_seconds"] for k in wall_pairs),
                      args.profile_b_name: sum(b[k]["wall_seconds"] for k in wall_pairs)},
                  "wall_measured_cases": len(wall_pairs)},
              "blind": {
                  "positive_recall": {
                      args.profile_a_name: metrics(a, "blind", "confirmed_planet", "recovered"),
                      args.profile_b_name: metrics(b, "blind", "confirmed_planet", "recovered"),
                      "paired": paired(a, b, "blind", "confirmed_planet", "recovered")},
                  "null_trigger": {
                      args.profile_a_name: metrics(a, "blind", "inverted_null", "trigger"),
                      args.profile_b_name: metrics(b, "blind", "inverted_null", "trigger"),
                      "paired": paired(a, b, "blind", "inverted_null", "trigger")}},
              "calibration_null_trigger": {
                  args.profile_a_name: metrics(a, "calibration", "inverted_null", "trigger"),
                  args.profile_b_name: metrics(b, "calibration", "inverted_null", "trigger")}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
