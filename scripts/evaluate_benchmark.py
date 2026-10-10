#!/usr/bin/env python3
"""Calibrate a threshold on calibration controls and evaluate a blind manifest.

Supported truth classes are confirmed_planet, inverted_null, and
known_false_positive. A case is one unique host/target; paired positive and
inverted-control rows share host_id but are summarized separately. The blind
split is never used to select the threshold.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
csv.field_size_limit(2**31 - 1)  # candidate rows carry long evidence/alias cells


def wilson(successes: int, trials: int, z: float = 1.959963984540054):
    if trials == 0:
        return None
    p = successes / trials
    den = 1 + z * z / trials
    center = (p + z * z / (2 * trials)) / den
    half = z * math.sqrt(p * (1 - p) / trials + z * z / (4 * trials * trials)) / den
    return [max(0.0, center - half), min(1.0, center + half)]


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def positive_number(row: dict, *keys: str) -> float:
    for key in keys:
        try:
            value = float(row[key])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(value):
            return value
    raise ValueError(f"missing finite numeric field; expected one of {keys}")


def ephemerides(row: dict):
    result = []
    try:
        result.append((float(row["period_days"]), float(row.get("epoch", "")),
                       float(row["duration_days"])))
    except (KeyError, TypeError, ValueError):
        pass
    raw = row.get("alias_ephemerides", "")
    for alias in raw.split(";"):
        fields = alias.split(":")
        if len(fields) >= 3:
            try:
                result.append((float(fields[0]), float(fields[1]), float(fields[2])))
            except ValueError:
                continue
    return result


def close_ephemeris(candidate, truth: dict, span: float, tolerance: float) -> bool:
    try:
        period = positive_number(truth, "period_days")
        epoch = positive_number(truth, "epoch_days", "epoch_btjd")
        duration = positive_number(truth, "duration_days")
    except ValueError:
        return False
    cp, ce, cd = candidate
    if period <= 0 or cp <= 0 or duration <= 0 or cd <= 0:
        return False
    # Error in the published ephemeris can be carried in the truth manifest.
    period_error = max(0.0, float(truth.get("period_error_days") or 0))
    epoch_error = max(0.0, float(truth.get("epoch_error_days") or 0))
    sigma_epoch = math.hypot(epoch_error, span * period_error / period)
    allowed = tolerance * duration + 3.0 * sigma_epoch
    drift = abs(cp - period) * span / period
    phase = abs((ce - epoch) % period)
    phase = min(phase, period - phase)
    duration_consistent = 0.5 <= cd / duration <= 2.0
    return drift <= allowed and phase <= allowed and duration_consistent


def read_scores(case: dict, manifest_dir: Path, minimum_score: float):
    candidate_path = (manifest_dir / case["candidates_csv"]).resolve()
    rows = [row for row in read_csv(candidate_path)
            if row.get("class") in ("periodic_dip", "duo_transit")]
    eligible = []
    for row in rows:
        try:
            score = float(row["score"])
            if math.isfinite(score) and score >= minimum_score:
                eligible.append((score, row))
        except (KeyError, ValueError):
            continue
    return rows, eligible


def summarize(cases: list[dict], threshold: float, manifest_dir: Path,
              ephemeris_tolerance: float):
    by_class: dict[str, list[dict]] = defaultdict(list)
    by_stratum: dict[str, list[dict]] = defaultdict(list)
    outcomes = []
    for case in cases:
        truth_class = case["truth_class"].strip().lower()
        if truth_class not in {"confirmed_planet", "inverted_null", "null", "known_false_positive"}:
            raise ValueError(f"{case.get('case_id')}: invalid truth_class {truth_class!r}")
        rows, _ = read_scores(case, manifest_dir, -math.inf)
        eligible = []
        for row in rows:
            try:
                score = float(row["score"])
                if math.isfinite(score) and score >= threshold:
                    eligible.append(row)
            except (KeyError, ValueError):
                continue
        is_planet = truth_class == "confirmed_planet"
        span = float(case.get("baseline_days") or 0)
        recovered = False
        if is_planet:
            if span <= 0:
                raise ValueError(f"{case.get('case_id')}: baseline_days is required")
            recovered = any(
                close_ephemeris(ephemeris, case, span, ephemeris_tolerance)
                for row in eligible for ephemeris in ephemerides(row))
        outcome = {
            "case_id": case.get("case_id", ""),
            "host_id": case.get("host_id", case.get("pair_id", case.get("case_id", ""))),
            "truth_class": truth_class,
            "trigger_count": len(eligible),
            "max_periodic_score": max((float(row["score"]) for row in rows
                                        if _finite_score(row)), default=None),
            "planet_recovered": recovered if is_planet else None,
            "any_trigger": bool(eligible),
            "mission": case.get("mission", "unknown"),
            "period_stratum": case.get("period_stratum", "unknown"),
            "depth_stratum": case.get("depth_stratum", "unknown"),
            "cadence_days": _optional_float(case.get("cadence_days")),
            "depth_ppm_proxy": _optional_float(case.get("depth_ppm_proxy")),
            "stellar_teff_k": _optional_float(case.get("stellar_teff_k")),
            "magnitude_stratum": case.get("magnitude_stratum", "unknown"),
        }
        outcomes.append(outcome)
        by_class[truth_class].append(outcome)
        by_stratum[f"{outcome['mission']}:{outcome['period_stratum']}:{outcome['depth_stratum']}:{truth_class}"].append(outcome)

    def rate(group, key):
        successes = sum(bool(row[key]) for row in group)
        return {"successes": successes, "trials": len(group),
                "fraction": successes / len(group) if group else None,
                "wilson_95": wilson(successes, len(group))}

    positives = by_class["confirmed_planet"]
    negatives = by_class["inverted_null"] + by_class["null"]
    hard_fp = by_class["known_false_positive"]
    return {
        "blind_case_count": len(cases),
        "confirmed_planet_recovery": rate(positives, "planet_recovered"),
        "inverted_or_null_control_trigger_rate": rate(negatives, "any_trigger"),
        "known_false_positive_trigger_rate": rate(hard_fp, "any_trigger"),
        "strata": {
            key: {
                "planet_recovery": rate(group, "planet_recovered") if key.endswith("confirmed_planet") else None,
                "trigger_rate": rate(group, "any_trigger") if not key.endswith("confirmed_planet") else None,
            }
            for key, group in sorted(by_stratum.items())
        },
        "case_outcomes": outcomes,
    }


def _finite_score(row):
    try:
        return math.isfinite(float(row["score"]))
    except (KeyError, TypeError, ValueError):
        return False


def _optional_float(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def choose_threshold(cases: list[dict], manifest_dir: Path, target_fpr: float,
                     minimum_score: float) -> tuple[float, dict]:
    controls = [case for case in cases
                if case["truth_class"].strip().lower() in {"inverted_null", "null"}]
    if not controls:
        raise ValueError("calibration split contains no inverted/null controls")
    maxima = []
    for case in controls:
        rows, _ = read_scores(case, manifest_dir, minimum_score)
        maxima.append(max((float(row["score"]) for row in rows if _finite_score(row)),
                          default=-math.inf))
    finite_scores = sorted({score for score in maxima if math.isfinite(score)})
    if not finite_scores:
        return minimum_score, {
            "control_cases": len(controls), "triggered_controls": 0,
            "empirical_trigger_rate": 0.0, "wilson_95": wilson(0, len(controls)),
            "target_trigger_rate": target_fpr,
            "note": "No finite periodic scores in calibration controls; using the configured score floor.",
        }
    thresholds = [*finite_scores, math.nextafter(finite_scores[-1], math.inf)]
    feasible = []
    for threshold in thresholds:
        alarms = sum(score >= threshold for score in maxima)
        fpr = alarms / len(maxima)
        if fpr <= target_fpr:
            feasible.append((threshold, alarms, fpr))
    threshold, alarms, fpr = min(feasible, key=lambda item: item[0])
    return threshold, {
        "control_cases": len(controls), "triggered_controls": alarms,
        "empirical_trigger_rate": fpr, "wilson_95": wilson(alarms, len(controls)),
        "target_trigger_rate": target_fpr,
    }


def validate_manifest(cases: list[dict], calibration_split: str, blind_split: str) -> None:
    ids = [row.get("case_id", "") for row in cases]
    if not all(ids) or len(ids) != len(set(ids)):
        raise ValueError("case_id values must be present and unique")
    host_splits: dict[str, set[str]] = defaultdict(set)
    host_classes: dict[tuple[str, str], int] = defaultdict(int)
    for row in cases:
        host = row.get("host_id") or row.get("pair_id") or row["case_id"]
        split = row.get("split", "")
        if split not in {calibration_split, blind_split}:
            raise ValueError(f"{row['case_id']}: split must be {calibration_split!r} or {blind_split!r}")
        truth = row.get("truth_class", "").lower()
        host_splits[host].add(split)
        if truth in {"confirmed_planet", "inverted_null", "null"}:
            host_classes[(host, truth)] += 1
    if any(len(splits) > 1 for splits in host_splits.values()):
        raise ValueError("a host appears in both calibration and blind splits")
    if any(count > 1 for count in host_classes.values()):
        raise ValueError("benchmark requires at most one case per host and truth class")
    positives = {host for host, truth in host_classes if truth == "confirmed_planet"}
    controls = {host for host, truth in host_classes if truth in {"inverted_null", "null"}}
    if positives and controls and positives != controls:
        raise ValueError("paired benchmark requires one positive and one inverted/null control per host")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path,
                        help="truth manifest; keep labels/ephemerides private until blind evaluation")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--minimum-score", type=float, default=0.0,
                        help="optional lower floor during calibration threshold search")
    parser.add_argument("--fixed-threshold", type=float,
                        help="evaluate a pre-registered threshold instead of calibrating one")
    parser.add_argument("--target-false-alarm-rate", type=float, default=0.05)
    parser.add_argument("--calibration-split", default="calibration")
    parser.add_argument("--evaluation-split", default="blind")
    parser.add_argument("--ephemeris-tolerance", type=float, default=1.0,
                        help="match tolerance in transit durations, plus 3-sigma catalog timing uncertainty")
    args = parser.parse_args()
    if (not math.isfinite(args.minimum_score)
            or not 0 < args.target_false_alarm_rate < 1
            or not math.isfinite(args.ephemeris_tolerance)
            or args.ephemeris_tolerance <= 0):
        parser.error("invalid score floor, false-alarm target, or ephemeris tolerance")
    manifest = args.manifest.expanduser().resolve()
    cases = read_csv(manifest)
    try:
        validate_manifest(cases, args.calibration_split, args.evaluation_split)
    except ValueError as error:
        parser.error(str(error))
    calibration = [row for row in cases if row.get("split") == args.calibration_split]
    blind = [row for row in cases if row.get("split") == args.evaluation_split]
    if not calibration or not blind:
        parser.error("manifest must contain nonempty calibration and blind splits")
    if args.fixed_threshold is None:
        threshold, calibration_summary = choose_threshold(
            calibration, manifest.parent, args.target_false_alarm_rate, args.minimum_score)
        calibration_summary["calibration_split"] = args.calibration_split
    else:
        if not math.isfinite(args.fixed_threshold):
            parser.error("fixed threshold must be finite")
        threshold = args.fixed_threshold
        calibration_summary = {"calibration_split": args.calibration_split,
                               "threshold_was_pre_registered": True}
    evaluation = summarize(blind, threshold, manifest.parent, args.ephemeris_tolerance)
    n_blind_hosts = len({row.get("host_id", row.get("pair_id", row.get("case_id")))
                         for row in blind if row["truth_class"].lower() == "confirmed_planet"})
    n_per_class = min(evaluation["confirmed_planet_recovery"]["trials"],
                      evaluation["inverted_or_null_control_trigger_rate"]["trials"])
    recovery_ci = evaluation["confirmed_planet_recovery"]["wilson_95"]
    report = {
        "manifest": str(manifest), "minimum_candidate_score_floor": args.minimum_score,
        "calibrated_threshold": threshold, "calibration": calibration_summary,
        "evaluation": evaluation, "blind_independent_positive_hosts": n_blind_hosts,
        "nominal_worst_case_95_percent_half_width": 1.96 * math.sqrt(0.25 / n_per_class) if n_per_class else None,
        "observed_recovery_wilson_95_half_width": (recovery_ci[1] - recovery_ci[0]) / 2 if recovery_ci else None,
        "underpowered_for_approximately_7pp_precision": n_per_class < 200,
        "threshold_selection": "Lowest score threshold whose calibration inverted/null-control trigger rate is at or below the target; blind data are not used to select it." if args.fixed_threshold is None else "Pre-registered fixed threshold.",
        "warning": "Known-answer recovery is conditional on published confirmed planets and selected public light curves. Inverted curves approximate noise/systematic controls and do not model the full astrophysical false-positive population. This is not a survey-wide completeness, reliability, or planet-confirmation estimate.",
    }
    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
