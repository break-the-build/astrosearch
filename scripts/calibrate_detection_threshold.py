#!/usr/bin/env python3
"""Calibrate the astrosearch detection tier from the known-answer benchmark.

The search reports every transit candidate with score >= 6, but on the
known-answer benchmark 97% of inverted-flux null controls produce one, so 6 is
a screening threshold, not a detection threshold. This script picks the
smallest score whose trigger rate on the calibration split's null controls is
at most the target false-alarm rate, then reports blind-split recall and null
trigger rate at that score (the blind split never influences the threshold).

Writes a JSON file read by ctoi_gate.py (criterion 5) and by the README's
detection-tier table. The calibration is specific to its context (data
product, search passes, binary); the file records that context.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import evaluate_benchmark as ev  # noqa: E402

csv.field_size_limit(2**31 - 1)


def case_maximum(case: dict, root: Path) -> float:
    rows, _ = ev.read_scores(case, root, 0.0)
    return max((float(r["score"]) for r in rows), default=0.0)


def matched_score(case: dict, root: Path) -> float | None:
    """Score of the first transit candidate matching the true ephemeris."""
    rows, _ = ev.read_scores(case, root, 0.0)
    span = float(case.get("baseline_days") or 0) or 1.0
    for row in rows:
        if any(ev.close_ephemeris(e, case, span, 1.0) for e in ev.ephemerides(row)):
            return float(row["score"])
    return None


def threshold_for(maxima: list[float], target: float) -> float:
    """Smallest observed score s with #(maxima >= s) / N <= target."""
    for score in sorted(set(maxima)):
        if sum(m >= score for m in maxima) / len(maxima) <= target:
            return score
    return math.nextafter(max(maxima), math.inf)


def rate(k: int, n: int) -> dict:
    return {"count": k, "total": n, "rate": k / n if n else None, "wilson_95": ev.wilson(k, n)}


def run_context(cases: list[dict], root: Path) -> dict:
    for case in cases:
        run = (root / case["candidates_csv"]).with_name("candidates.run.json")
        if run.exists():
            state = json.loads(run.read_text()).get("state", {})
            return {"binary_sha256": state.get("binary_sha256"),
                    "transit_search_passes": state.get("config", {}).get("transit_search_passes"),
                    "min_period": state.get("config", {}).get("min_period"),
                    "max_period": state.get("config", {}).get("max_period")}
    return {}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("manifest", type=Path, help="truth manifest (rows with missing results are skipped)")
    p.add_argument("--target", type=float, default=0.05, help="calibration null trigger rate")
    p.add_argument("--flux-origin", default="QLP", help="data product the benchmark used")
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    root = a.manifest.parent
    cases = [c for c in ev.read_csv(a.manifest) if (root / c["candidates_csv"]).exists()]
    skipped = len(ev.read_csv(a.manifest)) - len(cases)

    def select(split, truth):
        return [c for c in cases if c["split"] == split and c["truth_class"] == truth]

    cal_nulls = [case_maximum(c, root) for c in select("calibration", "inverted_null")]
    threshold = threshold_for(cal_nulls, a.target)
    out = {"threshold": threshold, "target_false_alarm_rate": a.target,
           "screening_threshold": 6.0,
           "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
           "manifest": str(a.manifest), "cases_used": len(cases), "cases_missing_results": skipped,
           "context": {**run_context(cases, root), "flux_origin": a.flux_origin},
           "calibration": {"null_trigger": rate(sum(m >= threshold for m in cal_nulls), len(cal_nulls))}}
    for split in ("calibration", "blind"):
        nulls = [case_maximum(c, root) for c in select(split, "inverted_null")]
        planets = [matched_score(c, root) for c in select(split, "confirmed_planet")]
        block = out.setdefault(split, {})
        block["null_trigger"] = rate(sum(m >= threshold for m in nulls), len(nulls))
        block["null_trigger_at_screening"] = rate(sum(m >= 6.0 for m in nulls), len(nulls))
        block["planet_recall"] = rate(sum(s is not None and s >= threshold for s in planets), len(planets))
        block["planet_in_candidate_list"] = rate(sum(s is not None for s in planets), len(planets))
        by_depth = {}
        for depth in ("shallow", "medium", "deep"):
            sub = [matched_score(c, root) for c in select(split, "confirmed_planet")
                   if c["depth_stratum"] == depth]
            by_depth[depth] = rate(sum(s is not None and s >= threshold for s in sub), len(sub))
        block["planet_recall_by_depth"] = by_depth
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: out[k] for k in ("threshold", "cases_used", "cases_missing_results")}))
    print("blind recall", out["blind"]["planet_recall"], "blind null", out["blind"]["null_trigger"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
