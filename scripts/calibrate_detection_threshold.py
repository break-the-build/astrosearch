#!/usr/bin/env python3
"""Calibrate the astrosearch detection tier from the known-answer benchmark.

The search reports every transit candidate with score >= 6, but on the
known-answer benchmark 97% of inverted-flux null controls produce one, so 6 is
a screening floor, not a detection threshold.

Two thresholds are written, both at the target null trigger rate:

* `evaluation`: chosen on the calibration split's null controls only; the
  blind split's planet recall and null trigger rate at that value are the
  honest out-of-sample performance;
* `threshold` (deployment): chosen on all null controls, after the
  evaluation, for use by the search (`--detection-threshold`) and the gate.

The statistic is the jackknife score by default (score with the strongest
single event removed; null maxima are dominated by one event). The file
records the statistic and its context (data product, passes, binary).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import evaluate_benchmark as ev
from vetting.tier import STATISTICS, statistic

csv.field_size_limit(2**31 - 1)


def case_values(case: dict, root: Path, name: str) -> tuple[float, float | None]:
    """(maximum statistic over transit candidates, statistic of the candidate
    matching the true ephemeris or None)."""
    rows, _ = ev.read_scores(case, root, 0.0)
    best = max((statistic(r, name) for r in rows), default=0.0)
    matched = None
    if case["truth_class"] == "confirmed_planet":
        span = float(case.get("baseline_days") or 0) or 1.0
        matched_scores = [statistic(r, name) for r in rows
                          if any(ev.close_ephemeris(e, case, span, 1.0)
                                 for e in ev.ephemerides(r))]
        matched = max(matched_scores, default=None)
    return best, matched


def required_cases(rows: list[dict], root: Path) -> tuple[list[dict], list[dict]]:
    """Return every case used for threshold selection/evaluation and intentionally
    unsearched calibration positives. Missing required outputs are an error, not
    silently dropped from a denominator.
    """
    required = [row for row in rows
                if (row["split"] == "calibration"
                    and row["truth_class"] in {"inverted_null", "null"})
                or row["split"] == "blind"]
    skipped = [row for row in rows
               if row["split"] == "calibration" and row["truth_class"] == "confirmed_planet"]
    missing = [row["case_id"] for row in required
               if not (root / row["candidates_csv"]).is_file()]
    if missing:
        preview = ", ".join(missing[:5])
        suffix = " …" if len(missing) > 5 else ""
        raise ValueError(f"missing results for {len(missing)} required benchmark cases: {preview}{suffix}")
    return required, skipped


def threshold_for(maxima: list[float], target: float) -> float:
    """Smallest observed value s with #(maxima >= s) / N <= target."""
    for value in sorted(set(maxima)):
        if sum(m >= value for m in maxima) / len(maxima) <= target:
            return value
    return math.nextafter(max(maxima), math.inf)


def rate(k: int, n: int) -> dict:
    return {"count": k, "total": n, "rate": k / n if n else None, "wilson_95": ev.wilson(k, n)}


def run_context(cases: list[dict], root: Path) -> dict:
    for case in cases:
        run = (root / case["candidates_csv"]).with_name("candidates.run.json")
        if run.exists():
            state = json.loads(run.read_text()).get("state", {})
            config = state.get("config", {})
            return {"binary_sha256": state.get("binary_sha256"),
                    "transit_search_passes": config.get("transit_search_passes"),
                    "min_period": config.get("min_period"), "max_period": config.get("max_period")}
    return {}


def split_metrics(records: list[dict], split: str, threshold: float) -> dict:
    nulls = [r for r in records if r["split"] == split and r["truth"] == "inverted_null"]
    planets = [r for r in records if r["split"] == split and r["truth"] == "confirmed_planet"]
    hit = lambda r: r["matched"] is not None and r["matched"] >= threshold  # noqa: E731
    return {"null_trigger": rate(sum(r["best"] >= threshold for r in nulls), len(nulls)),
            "null_trigger_at_screening_score": rate(sum(r["best_score"] >= 6.0 for r in nulls), len(nulls)),
            "planet_recall": rate(sum(hit(r) for r in planets), len(planets)),
            "planet_in_candidate_list": rate(sum(r["matched"] is not None for r in planets), len(planets)),
            "planet_recall_by_depth": {d: rate(sum(hit(r) for r in planets if r["depth"] == d),
                                               sum(r["depth"] == d for r in planets))
                                       for d in ("shallow", "medium", "deep")}}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("manifest", type=Path, help="truth manifest (rows with missing results are skipped)")
    p.add_argument("--statistic", choices=STATISTICS, default="jackknife_score")
    p.add_argument("--target", type=float, default=0.05, help="null trigger rate")
    p.add_argument("--flux-origin", default="QLP", help="data product the benchmark used")
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    root = a.manifest.parent
    all_rows = ev.read_csv(a.manifest)
    try:
        cases, skipped_calibration_positives = required_cases(all_rows, root)
    except ValueError as error:
        p.error(str(error))
    records = []
    for case in cases:
        best, matched = case_values(case, root, a.statistic)
        best_score, _ = case_values(case, root, "score") if a.statistic != "score" else (best, None)
        records.append({"split": case["split"], "truth": case["truth_class"],
                        "depth": case["depth_stratum"], "best": best, "matched": matched,
                        "best_score": best_score})
    nulls = lambda split=None: [r["best"] for r in records  # noqa: E731
                                if r["truth"] == "inverted_null" and (split is None or r["split"] == split)]
    evaluation = threshold_for(nulls("calibration"), a.target)
    deployment = threshold_for(nulls(), a.target)
    out = {"statistic": a.statistic, "threshold": deployment, "target_false_alarm_rate": a.target,
           "screening_threshold": 6.0,
           "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
           "manifest": str(a.manifest), "cases_used": len(cases),
           "case_coverage": {
               "required_cases": len(cases), "missing_required_results": 0,
               "calibration_positive_cases_not_searched_by_design": len(skipped_calibration_positives),
               "unreferenced_result_sidecars": max(0, len(list(root.glob("results/*/*.run.json"))) - len(cases)),
           },
           "context": {**run_context(cases, root), "flux_origin": a.flux_origin},
           "evaluation": {"threshold_from_calibration_split": evaluation,
                          "calibration": split_metrics(records, "calibration", evaluation),
                          "blind": split_metrics(records, "blind", evaluation)},
           "deployment": {"threshold_from_all_controls": deployment,
                          "all_null_trigger": rate(sum(v >= deployment for v in nulls()), len(nulls()))}}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
    blind = out["evaluation"]["blind"]
    print(json.dumps({"statistic": a.statistic, "evaluation_threshold": evaluation,
                      "deployment_threshold": deployment, "cases_used": len(cases)}))
    print("blind recall", blind["planet_recall"], "blind null", blind["null_trigger"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
