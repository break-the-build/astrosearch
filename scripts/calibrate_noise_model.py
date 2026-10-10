#!/usr/bin/env python3
"""Calibrate 1-hour photometric noise against TESS magnitude from real light curves.

For each known-answer benchmark light curve (QLP full-frame-image
photometry), the known planet's transits are masked, the flux is binned to
1 hour within observing segments, and the noise is the robust scatter of
adjacent 1-hour bin differences divided by sqrt(2), which ignores slow
trends. The output is the median 1-hour noise (ppm) per 0.5-magnitude bin,
used by target_yield.py to predict transit S/N for stars not yet downloaded.
"""
from __future__ import annotations

import argparse
import csv
import json
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from vetting.io import phase_offset, read_table

csv.field_size_limit(2**31 - 1)
HOUR = 1 / 24


def one_hour_noise(job: tuple[str, dict]) -> tuple[float, float, float] | None:
    root, case = job
    path = Path(root) / "lightcurves" / f"{case['case_id']}.csv"
    if not path.exists():  # excluded cases have no light curve
        return None
    data = read_table(path)
    time, flux = data[:, 0], data[:, 1]
    try:
        period, epoch = float(case["period_days"]), float(case["epoch_days"])
        duration = float(case["duration_days"])
        keep = np.abs(phase_offset(time, epoch, period)) > duration
        time, flux = time[keep], flux[keep]
    except (KeyError, ValueError):
        pass
    flux = flux / np.median(flux)
    bins = np.floor((time - time[0]) / HOUR).astype(np.int64)
    unique, index, counts = np.unique(bins, return_inverse=True, return_counts=True)
    means = np.bincount(index, flux) / counts
    full = counts >= max(2, int(0.5 * np.median(counts)))
    adjacent = np.flatnonzero(np.diff(unique) == 1)
    adjacent = adjacent[full[adjacent] & full[adjacent + 1]]
    if adjacent.size < 50:
        return None
    differences = means[adjacent + 1] - means[adjacent]
    sigma = 1.4826 * np.median(np.abs(differences - np.median(differences))) / np.sqrt(2)
    # Point scatter as injection_recovery.py measures it; beta is how much
    # worse the 1-hour noise is than white noise averaged over an hour.
    point = 1.4826 * np.median(np.abs(np.diff(flux))) / np.sqrt(2)
    per_hour = HOUR / np.median(np.diff(time))
    beta = sigma / (point / np.sqrt(per_hour))
    return float(case["tmag"]), float(sigma * 1e6), float(beta)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("metadata", type=Path, help="known-answer case-metadata.csv")
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--bin-width", type=float, default=0.5, help="magnitude bin width")
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    root = a.metadata.resolve().parent
    with a.metadata.open(newline="") as stream:
        tmag_by_pair = {c["case_id"]: c["tmag"] for c in csv.DictReader(stream) if c.get("tmag")}
    # Light curves are stored under the blinded truth-manifest IDs; pair_id
    # links each planet-bearing curve back to its catalogue metadata.
    with (root / "truth-manifest.csv").open(newline="") as stream:
        cases = [{**row, "tmag": tmag_by_pair[row["pair_id"]]} for row in csv.DictReader(stream)
                 if row["truth_class"] == "confirmed_planet" and row["pair_id"] in tmag_by_pair]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        measured = [m for m in pool.map(one_hour_noise, [(str(root), c) for c in cases], chunksize=4) if m]
    tmag = np.array([m[0] for m in measured])
    sigma = np.array([m[1] for m in measured])
    beta = np.array([m[2] for m in measured])
    edges = np.arange(np.floor(tmag.min() / a.bin_width) * a.bin_width,
                      tmag.max() + a.bin_width, a.bin_width)
    table = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        inside = (tmag >= lo) & (tmag < hi)
        if inside.sum() >= 5:
            table.append({"tmag": round(float(np.median(tmag[inside])), 3),
                          "sigma_1h_ppm": round(float(np.median(sigma[inside])), 1),
                          "n": int(inside.sum()),
                          "sigma_1h_ppm_p16": round(float(np.percentile(sigma[inside], 16)), 1),
                          "sigma_1h_ppm_p84": round(float(np.percentile(sigma[inside], 84)), 1),
                          "beta_1h": round(float(np.median(beta[inside])), 3)})
    out = {"created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
           "source": str(a.metadata), "product": "QLP full-frame-image light curves",
           "method": "masked known planet; 1-h bins; robust std of adjacent-bin differences / sqrt(2)",
           "beta_1h": "sigma_1h / (point scatter / sqrt(points per hour)); point scatter as in injection_recovery.py",
           "light_curves": len(measured), "table": table}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=1) + "\n")
    for row in table:
        print(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
