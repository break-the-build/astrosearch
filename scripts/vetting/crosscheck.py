"""Independent confirmation: a separate search algorithm and separate products.

Agreement between astrosearch and an independent Box Least Squares search
(astropy, on a Savitzky-Golay-flattened light curve) or Transit Least Squares,
and between data products reducing the same pixels differently, is stronger
evidence than any single score.
"""
from __future__ import annotations

import numpy as np

from .diagnostics import combine, event_depths
from .io import LightCurve, phase_offset, read_lightcurve


def flatten(lc: LightCurve, epoch: float, period: float, duration: float,
            window_days: float = 1.0) -> np.ndarray:
    """Per-part running-median flattening with the candidate's transits
    excluded from the trend, so the dips are not absorbed."""
    flat = np.empty_like(lc.flux)
    masked = np.abs(phase_offset(lc.time, epoch, period)) <= duration
    for part in np.unique(lc.part):
        idx = np.flatnonzero(lc.part == part)
        t, f, m = lc.time[idx], lc.flux[idx], masked[idx]
        trend = np.empty(len(idx))
        half = window_days / 2
        lo = hi = 0
        for j, tj in enumerate(t):
            while t[lo] < tj - half:
                lo += 1
            while hi < len(t) and t[hi] <= tj + half:
                hi += 1
            window = f[lo:hi][~m[lo:hi]]
            trend[j] = np.median(window) if window.size else np.median(f[lo:hi])
        flat[idx] = f / trend
    return flat


def bls_recovery(lc: LightCurve, period: float, epoch: float, duration: float,
                 min_period: float = 0.5, max_period: float | None = None,
                 peaks: int = 5) -> dict:
    """Run astropy BLS over the full range; recovered when any of the top
    peaks matches P, P/2, or 2P within 0.5% and phase within a duration."""
    from astropy.timeseries import BoxLeastSquares

    flux = flatten(lc, epoch, period, duration)
    span = lc.time[-1] - lc.time[0]
    max_period = min(max_period or max(30.0, 2.5 * period), span / 2)
    # Bin to a sixth of the duration: the box shape is preserved and N drops.
    width = duration / 6
    keys = np.floor((lc.time - lc.time[0]) / width).astype(np.int64)
    _, index, counts = np.unique(keys, return_inverse=True, return_counts=True)
    weights = 1 / lc.error ** 2
    w_sum = np.bincount(index, weights)
    time = np.bincount(index, weights * lc.time) / w_sum
    binned = np.bincount(index, weights * flux) / w_sum
    error = 1 / np.sqrt(w_sum)
    bls = BoxLeastSquares(time, binned, dy=error)
    durations = np.array([0.67, 1.0, 1.5]) * duration
    durations = durations[durations < min_period / 2]
    # Log-spaced periods whose accumulated phase drift over the baseline stays
    # within the shortest duration (astropy's autoperiod is far denser here).
    step = durations.min() / span
    grid = np.exp(np.arange(np.log(min_period), np.log(max_period), step))
    result = bls.power(grid, durations, objective="likelihood")
    power = np.asarray(result.power)
    sde = (power - power.mean()) / power.std()
    order = np.argsort(power)[::-1]
    distinct, used = [], []
    for i in order:
        p = float(result.period[i])
        if any(abs(p - q) / q < 0.01 for q in used):
            continue
        used.append(p)
        distinct.append({"period": p, "epoch": float(result.transit_time[i]),
                         "depth_ppm": float(result.depth[i] * 1e6), "sde": float(sde[i])})
        if len(distinct) == 50:
            break
    top = distinct[:peaks]
    near = np.zeros(len(grid), dtype=bool)
    for ratio in (1.0, 0.5, 2.0):
        near |= np.abs(np.asarray(result.period) - period * ratio) / (period * ratio) <= 0.005
    candidate_rank = next((n for n, peak in enumerate(distinct, start=1)
                           if any(abs(peak["period"] - period * r) / (period * r) <= 0.005
                                  for r in (1.0, 0.5, 2.0))), None)
    candidate_sde = float(sde[near].max()) if near.any() else None
    match = None
    for rank, peak in enumerate(top, start=1):
        for ratio in (1.0, 0.5, 2.0):
            target = period * ratio
            if abs(peak["period"] - target) / target <= 0.005:
                phase = abs(phase_offset(np.array([peak["epoch"]]), epoch, min(target, period))[0])
                if phase <= duration:
                    match = {"rank": rank, "alias_ratio": ratio, **peak}
                    break
        if match:
            break
    return {"status": "PASS" if match else "FAIL",
            "algorithm": "astropy BoxLeastSquares (likelihood), binned to duration/6",
            "periods_searched": int(len(grid)), "top_peaks": top, "match": match,
            "candidate_peak_rank": candidate_rank, "candidate_sde": candidate_sde,
            "rule": "PASS when P, P/2, or 2P is among the top 5 distinct peaks"}


def product_depths(main: LightCurve, others: list, epoch: float, period: float,
                   duration: float) -> dict:
    """Depth at the ephemeris in each independent product; consistent when the
    chi^2 across products has p >= 0.01."""
    from scipy import stats

    rows = [{"product": "primary", **combine(event_depths(main, epoch, period, duration))}]
    for path in others:
        lc = read_lightcurve(path)
        rows.append({"product": str(path), **combine(event_depths(lc, epoch, period, duration))})
    usable = [r for r in rows if r.get("depth_ppm") is not None]
    if len(usable) < 2:
        return {"products": rows, "status": "NOT_EVALUATED",
                "reason": "needs the signal covered in at least two products"}
    d = np.array([r["depth_ppm"] for r in usable])
    w = 1 / np.array([r["sigma_ppm"] for r in usable]) ** 2
    mean = (w * d).sum() / w.sum()
    chi2 = float((w * (d - mean) ** 2).sum())
    p = float(stats.chi2.sf(chi2, len(usable) - 1))
    detected = all(r["snr"] >= 3 for r in usable)
    return {"products": rows, "chi2": chi2, "p_value": p, "all_detected_3sigma": detected,
            "status": "PASS" if p >= 0.01 and detected else "FAIL"}
