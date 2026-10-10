#!/usr/bin/env python3
"""Box Least Squares baseline on the known-answer benchmark (Experiment F).

Runs astropy BoxLeastSquares on exactly the light curves astrosearch searched
and scores both methods the same way: the per-case statistic is the best
transit-like peak; the detection threshold is the calibration-split null
control quantile at the target trigger rate; recall is the fraction of blind
planets whose true ephemeris matches a peak at or above that threshold
(evaluate_benchmark.close_ephemeris, the rule used for astrosearch). The blind
split never informs the threshold.

BLS setup (documented, not tuned on the blind split): per-segment 1-day
running-median flattening, binning to 10 minutes, durations 1-8 h, periods
from 0.5 d to min(30 d, span/2) on a log grid whose phase drift over the
baseline stays within the shortest duration, likelihood objective, SDE
statistic, the top 20 distinct peaks kept per case.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import evaluate_benchmark as ev
from vetting.fap import segments
from vetting.io import read_table

csv.field_size_limit(2**31 - 1)
DURATIONS_H = (1.0, 2.0, 4.0, 8.0)
BIN_DAYS = 10 / 1440


def flatten(time: np.ndarray, flux: np.ndarray, window: float = 1.0) -> np.ndarray:
    out = np.empty_like(flux)
    for idx in segments(time):
        t, f = time[idx], flux[idx]
        lo = np.searchsorted(t, t - window / 2)
        hi = np.searchsorted(t, t + window / 2, side="right")
        out[idx] = f / np.array([np.median(f[a:b]) for a, b in zip(lo, hi)])
    return out


def bin_series(time, flux, error, width=BIN_DAYS):
    keys = np.floor((time - time[0]) / width).astype(np.int64)
    _, index = np.unique(keys, return_inverse=True)
    w = 1 / error ** 2
    w_sum = np.bincount(index, w)
    return (np.bincount(index, w * time) / w_sum, np.bincount(index, w * flux) / w_sum,
            1 / np.sqrt(w_sum))


def search(path: Path, peaks: int = 20) -> dict:
    from astropy.timeseries import BoxLeastSquares

    data = read_table(path)
    data = data[np.argsort(data[:, 0], kind="stable")]
    time, flux, error = data[:, 0], data[:, 1], data[:, 2]
    good = np.isfinite(flux) & np.isfinite(error) & (error > 0) & (flux > 0)
    time, flux, error = time[good], flux[good], error[good]
    flux = flatten(time, flux)
    time, flux, error = bin_series(time, flux, error)
    span = time[-1] - time[0]
    durations = np.array(DURATIONS_H) / 24
    max_period = min(30.0, span / 2)
    step = durations.min() / span
    grid = np.exp(np.arange(np.log(0.5), np.log(max_period), step))
    result = BoxLeastSquares(time, flux, dy=error).power(grid, durations, objective="likelihood")
    power = np.asarray(result.power)
    sde = (power - power.mean()) / power.std()
    out, used = [], []
    for i in np.argsort(power)[::-1]:
        p = float(result.period[i])
        if any(abs(p - q) / q < 0.01 for q in used):
            continue
        used.append(p)
        if result.depth[i] <= 0:
            continue
        out.append({"period": p, "epoch": float(result.transit_time[i]),
                    "duration": float(result.duration[i]), "sde": float(sde[i])})
        if len(out) == peaks:
            break
    return {"peaks": out, "periods_searched": int(len(grid)), "points_binned": int(len(time))}


def run_case(job: tuple[dict, str, str]) -> dict:
    case, root, cache = job
    cache_path = Path(cache) / f"{case['case_id']}.json"
    if cache_path.exists():
        result = json.loads(cache_path.read_text())
    else:
        start = time.perf_counter()
        result = search(Path(root) / "lightcurves" / f"{case['case_id']}.csv")
        result["seconds"] = time.perf_counter() - start
        cache_path.write_text(json.dumps(result))
    span = float(case.get("baseline_days") or 0) or 1.0
    best = max((p["sde"] for p in result["peaks"]), default=-math.inf)
    matched = None
    if case["truth_class"] == "confirmed_planet":
        matched = next((p["sde"] for p in result["peaks"]
                        if ev.close_ephemeris((p["period"], p["epoch"], p["duration"]), case, span, 1.0)),
                       None)
    return {"case_id": case["case_id"], "split": case["split"], "truth": case["truth_class"],
            "depth": case["depth_stratum"], "best": best, "matched": matched,
            "seconds": result.get("seconds")}


def astrosearch_rows(cases: list[dict], root: Path) -> list[dict]:
    rows = []
    for case in cases:
        transit, _ = ev.read_scores(case, root, 0.0)
        span = float(case.get("baseline_days") or 0) or 1.0
        matched = None
        if case["truth_class"] == "confirmed_planet":
            matched = next((float(r["score"]) for r in transit
                            if any(ev.close_ephemeris(e, case, span, 1.0) for e in ev.ephemerides(r))),
                           None)
        rows.append({"case_id": case["case_id"], "split": case["split"], "truth": case["truth_class"],
                     "depth": case["depth_stratum"],
                     "best": max((float(r["score"]) for r in transit), default=-math.inf),
                     "matched": matched})
    return rows


def matched_far(rows: list[dict], target: float) -> dict:
    cal = sorted(r["best"] for r in rows if r["split"] == "calibration" and r["truth"] == "inverted_null")
    threshold = next((s for s in sorted(set(cal)) if sum(c >= s for c in cal) / len(cal) <= target),
                     math.inf)

    def rate(k, n):
        return {"count": k, "total": n, "rate": k / n if n else None, "wilson_95": ev.wilson(k, n)}

    out = {"threshold": threshold}
    for split in ("calibration", "blind"):
        planets = [r for r in rows if r["split"] == split and r["truth"] == "confirmed_planet"]
        nulls = [r for r in rows if r["split"] == split and r["truth"] == "inverted_null"]
        out[split] = {
            "recall": rate(sum(r["matched"] is not None and r["matched"] >= threshold for r in planets),
                           len(planets)),
            "null_trigger": rate(sum(r["best"] >= threshold for r in nulls), len(nulls)),
            "recall_by_depth": {d: rate(sum(r["matched"] is not None and r["matched"] >= threshold
                                            for r in planets if r["depth"] == d),
                                        sum(r["depth"] == d for r in planets))
                                for d in ("shallow", "medium", "deep")},
        }
    return out


def paired_difference(a: list[dict], b: list[dict], threshold_a: float, threshold_b: float) -> dict:
    """McNemar test on blind planets detected by one method but not the other."""
    from scipy import stats

    by_id = {r["case_id"]: r for r in b}
    only_a = only_b = 0
    for r in a:
        if r["split"] != "blind" or r["truth"] != "confirmed_planet" or r["case_id"] not in by_id:
            continue
        hit_a = r["matched"] is not None and r["matched"] >= threshold_a
        other = by_id[r["case_id"]]
        hit_b = other["matched"] is not None and other["matched"] >= threshold_b
        only_a += hit_a and not hit_b
        only_b += hit_b and not hit_a
    n = only_a + only_b
    p = float(stats.binomtest(only_a, n, 0.5).pvalue) if n else 1.0
    return {"astrosearch_only": only_a, "bls_only": only_b, "mcnemar_exact_p": p}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("manifest", type=Path)
    p.add_argument("--workers", type=int, default=12)
    p.add_argument("--target", type=float, default=0.05)
    p.add_argument("--cache", type=Path, required=True, help="per-case BLS result directory")
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    root = a.manifest.parent
    cases = [c for c in ev.read_csv(a.manifest) if (root / c["candidates_csv"]).exists()]
    a.cache.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        bls = list(pool.map(run_case, [(c, str(root), str(a.cache)) for c in cases], chunksize=2))
    ours = astrosearch_rows(cases, root)
    bls_summary, our_summary = matched_far(bls, a.target), matched_far(ours, a.target)
    out = {"cases": len(cases), "target_false_alarm_rate": a.target,
           "astrosearch": our_summary, "bls": bls_summary,
           "paired_blind": paired_difference(ours, bls, our_summary["threshold"], bls_summary["threshold"]),
           "bls_seconds_median": float(np.median([r["seconds"] for r in bls if r["seconds"]] or [0]))}
    a.output.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: out[k] for k in ("cases", "paired_blind")}))
    for name in ("astrosearch", "bls"):
        print(name, out[name]["threshold"], out[name]["blind"]["recall"], out[name]["blind"]["null_trigger"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
