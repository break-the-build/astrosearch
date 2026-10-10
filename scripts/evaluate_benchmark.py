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
import re
from collections import defaultdict
from pathlib import Path

from vetting.tier import STATISTICS, statistic as candidate_statistic
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


def qlp_family(products_json: str) -> tuple[str, int]:
    try:
        products = json.loads(products_json or "[]")
    except (TypeError, json.JSONDecodeError):
        return "unknown", 0
    sectors = set()
    for product in products if isinstance(products, list) else []:
        for key in ("sector", "sequence_number", "productFilename", "obs_id"):
            raw = product.get(key, "")
            try:
                sector = int(float(raw)) if key in {"sector", "sequence_number"} else None
            except (TypeError, ValueError, OverflowError):
                sector = None
            if sector is None:
                match = re.search(r"(?:s|sector[_-]?)(\d{1,4})", str(raw), re.IGNORECASE)
                sector = int(match.group(1)) if match else None
            if sector and 1 <= sector <= 999:
                sectors.add(sector)
                break
    if not sectors:
        return "unknown", 0
    families = {"TGLC-based QLP (S94+)" if sector >= 94 else "legacy QLP"
                for sector in sectors}
    return (next(iter(families)) if len(families) == 1 else "mixed QLP eras",
            len(sectors))


def enrich_product_strata(cases: list[dict], metadata_path: Path | None) -> None:
    """Backfill product-era strata from the benchmark's host metadata snapshot."""
    if metadata_path is None:
        return
    metadata = read_csv(metadata_path)
    by_host = {row.get("host_id", ""): row for row in metadata}
    for case in cases:
        host = case.get("host_id", "")
        info = by_host.get(host)
        if not info:
            continue
        family, sector_count = qlp_family(info.get("mast_products_json", ""))
        case.setdefault("extraction_family", family)
        case.setdefault("sector_count", sector_count)
        case.setdefault("cadence_stratum", info.get("cadence_stratum", "unknown"))


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


def read_scores(case: dict, manifest_dir: Path, minimum_score: float,
                statistic: str = "score"):
    candidate_path = (manifest_dir / case["candidates_csv"]).resolve()
    rows = [row for row in read_csv(candidate_path)
            if row.get("class") in ("periodic_dip", "duo_transit")]
    eligible = []
    for row in rows:
        try:
            score = candidate_statistic(row, statistic)
            if math.isfinite(score) and score >= minimum_score:
                eligible.append((score, row))
        except (KeyError, TypeError, ValueError):
            continue
    return rows, eligible


def summarize(cases: list[dict], threshold: float, manifest_dir: Path,
              ephemeris_tolerance: float, statistic: str = "score"):
    by_class: dict[str, list[dict]] = defaultdict(list)
    by_stratum: dict[str, list[dict]] = defaultdict(list)
    outcomes = []
    for case in cases:
        truth_class = case["truth_class"].strip().lower()
        if truth_class not in {"confirmed_planet", "inverted_null", "null", "known_false_positive"}:
            raise ValueError(f"{case.get('case_id')}: invalid truth_class {truth_class!r}")
        rows, scored_rows = read_scores(case, manifest_dir, -math.inf, statistic)
        eligible = [row for score, row in scored_rows if score >= threshold]
        is_planet = truth_class == "confirmed_planet"
        span = float(case.get("baseline_days") or 0)
        recovered = False
        matching_rows = []
        if is_planet:
            if span <= 0:
                raise ValueError(f"{case.get('case_id')}: baseline_days is required")
            matching_rows = [row for row in rows
                             if any(close_ephemeris(ephemeris, case, span, ephemeris_tolerance)
                                    for ephemeris in ephemerides(row))]
            recovered = any(row in eligible for row in matching_rows)
        match_rank = None
        if matching_rows:
            ranks = []
            for row in matching_rows:
                try:
                    ranks.append(int(row.get("class_rank") or row.get("rank") or 0))
                except (TypeError, ValueError):
                    pass
            match_rank = min((rank for rank in ranks if rank > 0), default=None)
        outcome = {
            "case_id": case.get("case_id", ""),
            "host_id": case.get("host_id", case.get("pair_id", case.get("case_id", ""))),
            "truth_class": truth_class,
            "false_positive_type": case.get("false_positive_type", "unspecified"),
            "trigger_count": len(eligible),
            "max_periodic_score": max((float(row["score"]) for row in rows
                                        if _finite_score(row)), default=None),
            "max_statistic": max((candidate_statistic(row, statistic) for row in rows
                                  if _finite_statistic(row, statistic)), default=None),
            "planet_recovered": recovered if is_planet else None,
            "matched_candidate_rank": match_rank,
            "any_trigger": bool(eligible),
            "mission": case.get("mission", "unknown"),
            "extraction_family": case.get("extraction_family", "unknown"),
            "cadence_stratum": case.get("cadence_stratum", "unknown"),
            "sector_count": _optional_float(case.get("sector_count")),
            "coverage_stratum": ("short" if span < 90 else "medium" if span < 365
                                 else "long"),
            "period_stratum": case.get("period_stratum", "unknown"),
            "depth_stratum": case.get("depth_stratum", "unknown"),
            "cadence_days": _optional_float(case.get("cadence_days")),
            "depth_ppm_proxy": _optional_float(case.get("depth_ppm_proxy")),
            "stellar_teff_k": _optional_float(case.get("stellar_teff_k")),
            "magnitude_stratum": case.get("magnitude_stratum", "unknown"),
        }
        outcomes.append(outcome)
        by_class[truth_class].append(outcome)
        sector_count = outcome["sector_count"]
        sector_band = ("unknown" if sector_count is None else "1-3" if sector_count <= 3
                       else "4-10" if sector_count <= 10 else "11+")
        dimensions = {
            "mission": outcome["mission"],
            "extraction_family": outcome["extraction_family"],
            "cadence": outcome["cadence_stratum"],
            "baseline": outcome["coverage_stratum"],
            "sector_count": sector_band,
            "period": outcome["period_stratum"],
            "depth": outcome["depth_stratum"],
            "magnitude": outcome["magnitude_stratum"],
        }
        for dimension, value in dimensions.items():
            by_stratum[f"{dimension}={value}:{truth_class}"].append(outcome)

    def rate(group, key):
        successes = sum(bool(row[key]) for row in group)
        return {"successes": successes, "trials": len(group),
                "fraction": successes / len(group) if group else None,
                "wilson_95": wilson(successes, len(group))}

    positives = by_class["confirmed_planet"]
    negatives = by_class["inverted_null"] + by_class["null"]
    hard_fp = by_class["known_false_positive"]
    hard_fp_by_type: dict[str, list[dict]] = defaultdict(list)
    for row in hard_fp:
        hard_fp_by_type[row["false_positive_type"]].append(row)
    def rank_rate(max_rank):
        group = [row for row in positives
                 if row["matched_candidate_rank"] is not None
                 and (max_rank is None or row["matched_candidate_rank"] <= max_rank)]
        successes, trials = len(group), len(positives)
        return {"successes": successes, "trials": trials,
                "fraction": successes / trials if trials else None,
                "wilson_95": wilson(successes, trials)}
    return {
        "blind_case_count": len(cases),
        "confirmed_planet_recovery": rate(positives, "planet_recovered"),
        "inverted_or_null_control_trigger_rate": rate(negatives, "any_trigger"),
        "known_false_positive_trigger_rate": rate(hard_fp, "any_trigger"),
        "known_false_positive_trigger_rate_by_type": {
            label: rate(group, "any_trigger")
            for label, group in sorted(hard_fp_by_type.items())
        },
        "planet_recovery_by_candidate_rank": {
            "top_1": rank_rate(1), "top_5": rank_rate(5),
            "top_10": rank_rate(10), "any_rank": rank_rate(None),
        },
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


def _finite_statistic(row, statistic: str) -> bool:
    try:
        return math.isfinite(candidate_statistic(row, statistic))
    except (KeyError, TypeError, ValueError):
        return False


def _optional_float(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def choose_threshold(cases: list[dict], manifest_dir: Path, target_fpr: float,
                     minimum_score: float,
                     statistic: str = "score") -> tuple[float, dict]:
    controls = [case for case in cases
                if case["truth_class"].strip().lower() in {"inverted_null", "null"}]
    if not controls:
        raise ValueError("calibration split contains no inverted/null controls")
    maxima = []
    for case in controls:
        _, eligible = read_scores(case, manifest_dir, minimum_score, statistic)
        maxima.append(max((score for score, _ in eligible), default=-math.inf))
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
    parser.add_argument("--case-metadata", type=Path,
                        help="optional frozen metadata CSV to backfill QLP processing-era strata")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--minimum-score", type=float, default=0.0,
                        help="optional lower floor for the selected statistic")
    parser.add_argument("--statistic", choices=STATISTICS, default="jackknife_score",
                        help="candidate statistic to calibrate and evaluate (default: production detection tier)")
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
    metadata_path = args.case_metadata.expanduser().resolve() if args.case_metadata else None
    if metadata_path is not None and not metadata_path.is_file():
        parser.error("--case-metadata file does not exist")
    enrich_product_strata(cases, metadata_path)
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
            calibration, manifest.parent, args.target_false_alarm_rate,
            args.minimum_score, args.statistic)
        calibration_summary["calibration_split"] = args.calibration_split
    else:
        if not math.isfinite(args.fixed_threshold):
            parser.error("fixed threshold must be finite")
        threshold = args.fixed_threshold
        calibration_summary = {"calibration_split": args.calibration_split,
                               "threshold_was_pre_registered": True}
    evaluation = summarize(blind, threshold, manifest.parent,
                           args.ephemeris_tolerance, args.statistic)
    n_blind_hosts = len({row.get("host_id", row.get("pair_id", row.get("case_id")))
                         for row in blind if row["truth_class"].lower() == "confirmed_planet"})
    n_per_class = min(evaluation["confirmed_planet_recovery"]["trials"],
                      evaluation["inverted_or_null_control_trigger_rate"]["trials"])
    recovery_ci = evaluation["confirmed_planet_recovery"]["wilson_95"]
    report = {
        "manifest": str(manifest), "minimum_candidate_score_floor": args.minimum_score,
        "statistic": args.statistic, "calibrated_threshold": threshold,
        "calibration": calibration_summary,
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
    console_evaluation = {key: value for key, value in evaluation.items()
                          if key != "case_outcomes"}
    print(json.dumps({"statistic": report["statistic"],
                      "calibrated_threshold": threshold,
                      "calibration": calibration_summary,
                      "evaluation": console_evaluation,
                      "underpowered_for_approximately_7pp_precision":
                          report["underpowered_for_approximately_7pp_precision"]},
                     indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
