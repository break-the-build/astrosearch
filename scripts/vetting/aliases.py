"""Period aliases of two-transit (duo) candidates and future transit windows.

With two transits at t1 < t2, any period P_k = (t2 - t1) / k is possible.
An alias is ruled out when it predicts a transit inside observed data that
shows no dip: the measured depth there is more than 3 sigma below the duo's
depth. Surviving aliases each get predicted windows.
"""
from __future__ import annotations

import numpy as np

from .diagnostics import combine, event_depths
from .io import LightCurve


def duo_aliases(lc: LightCurve, t1: float, t2: float, duration: float, depth_ppm: float,
                depth_sigma_ppm: float, min_period: float = 0.5, max_aliases: int = 200) -> list[dict]:
    separation = t2 - t1
    out = []
    for k in range(1, max_aliases + 1):
        period = separation / k
        if period < min_period:
            break
        events = event_depths(lc, t1, period, duration)
        others = [e for e in events
                  if abs(e["midpoint"] - t1) > duration and abs(e["midpoint"] - t2) > duration]
        contradicted = [e for e in others
                        if e["depth_ppm"] < depth_ppm - 3 * np.hypot(e["sigma_ppm"], depth_sigma_ppm)]
        summary = combine(others)
        out.append({"k": k, "period_days": period, "covered_other_epochs": len(others),
                    "contradicting_epochs": len(contradicted),
                    "other_epochs_depth_ppm": summary["depth_ppm"],
                    "allowed": not contradicted})
    return out


def predicted_windows(epoch: float, period: float, after: float, through: float,
                      epoch_sigma: float = 0.0, period_sigma: float = 0.0,
                      covariance: float = 0.0, duration: float = 0.0) -> list[dict]:
    """Mid-transit times in (after, through] with 1-sigma timing errors
    propagated from (epoch, period, covariance), plus +/-3 sigma windows that
    include the transit duration."""
    first = int(np.ceil((after - epoch) / period))
    last = int(np.floor((through - epoch) / period))
    windows = []
    for n in range(first, last + 1):
        mid = epoch + n * period
        sigma = float(np.sqrt(max(epoch_sigma ** 2 + n * n * period_sigma ** 2
                                  + 2 * n * covariance, 0.0)))
        half = 3 * sigma + duration / 2
        windows.append({"epoch_number": n, "midpoint_btjd": mid, "sigma_days": sigma,
                        "window_start_btjd": mid - half, "window_end_btjd": mid + half})
    return windows
