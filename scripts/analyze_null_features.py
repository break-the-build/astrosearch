#!/usr/bin/env python3
"""What distinguishes high-scoring null-control candidates from real planets?

For each benchmark case, takes the top transit candidate of each inverted
null control and the matched candidate of each confirmed planet, and computes
per-candidate features from the search CSV:

* single-event share: largest per-event SNR^2 over the sum (1 = one event
  carries the whole signal);
* jackknife score: score * sqrt(1 - share), i.e. the signal left after
  dropping the strongest event, assuming score^2 ~ sum of event SNR^2;
* n_events, odd_even_sigma, noise_factor, white_score / score.

Reports the separation (ROC AUC, planets vs null maxima) for each feature on
the calibration split only, so feature choice does not peek at the blind split.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

import evaluate_benchmark as ev

csv.field_size_limit(2**31 - 1)


def event_snrs(row: dict) -> list[float]:
    out = []
    for item in (row.get("event_evidence") or "").split(";"):
        fields = item.split(":")
        if len(fields) >= 3:
            try:
                out.append(float(fields[2]))
            except ValueError:
                continue
    return out


def features(row: dict) -> dict:
    score = float(row["score"])
    snr = np.array(event_snrs(row))
    positive = np.clip(snr, 0, None) ** 2
    share = float(positive.max() / positive.sum()) if positive.size and positive.sum() > 0 else 1.0
    number = lambda key: float(row[key]) if row.get(key) not in (None, "") else math.nan  # noqa: E731
    return {"score": score, "single_event_share": share,
            "jackknife_score": score * math.sqrt(max(0.0, 1 - share)),
            "n_events": float(row["n_events"]), "odd_even_sigma": number("odd_even_sigma"),
            "noise_factor": number("noise_factor"),
            "white_over_score": number("white_score") / score if score else math.nan}


def auc(positives: list[float], negatives: list[float]) -> float:
    pos, neg = np.array(positives), np.array(negatives)
    pos, neg = pos[np.isfinite(pos)], neg[np.isfinite(neg)]
    if not len(pos) or not len(neg):
        return math.nan
    greater = (pos[:, None] > neg[None, :]).sum() + 0.5 * (pos[:, None] == neg[None, :]).sum()
    return float(greater / (len(pos) * len(neg)))


def collect(cases: list[dict], root: Path) -> list[dict]:
    rows = []
    for case in cases:
        transit, _ = ev.read_scores(case, root, 0.0)
        if not transit:
            continue
        if case["truth_class"] == "inverted_null":
            chosen = max(transit, key=lambda r: float(r["score"]))
        else:
            span = float(case.get("baseline_days") or 0) or 1.0
            chosen = next((r for r in transit
                           if any(ev.close_ephemeris(e, case, span, 1.0) for e in ev.ephemerides(r))), None)
            if chosen is None:
                continue
        rows.append({"case_id": case["case_id"], "split": case["split"], "truth": case["truth_class"],
                     "depth": case["depth_stratum"], **features(chosen)})
    return rows


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("manifest", type=Path)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    root = a.manifest.parent
    cases = [c for c in ev.read_csv(a.manifest) if (root / c["candidates_csv"]).exists()]
    rows = collect(cases, root)
    keys = ["score", "jackknife_score", "single_event_share", "n_events", "odd_even_sigma",
            "noise_factor", "white_over_score"]
    cal = [r for r in rows if r["split"] == "calibration"]
    planets = [r for r in cal if r["truth"] == "confirmed_planet"]
    nulls = [r for r in cal if r["truth"] == "inverted_null"]
    summary = {k: {"auc_planet_vs_null": auc([r[k] for r in planets], [r[k] for r in nulls]),
                   "planet_median": float(np.nanmedian([r[k] for r in planets])),
                   "null_median": float(np.nanmedian([r[k] for r in nulls]))} for k in keys}
    a.output.write_text(json.dumps({"calibration_split": summary, "rows": rows}, indent=1) + "\n")
    for k, v in summary.items():
        print(f"{k:20s} AUC {v['auc_planet_vs_null']:.3f}  planet {v['planet_median']:.3f}  null {v['null_median']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
