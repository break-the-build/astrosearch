"""Light-curve vetting diagnostics at a fixed ephemeris.

Depths use a local baseline per event (1.5-3 durations from mid-event), so
slow trends in the un-detrended, sector-normalized flux do not bias them.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import stats

from .io import LightCurve, phase_offset


def event_depths(lc: LightCurve, epoch: float, period: float, duration: float,
                 min_in: int = 2, min_out: int = 4) -> list[dict]:
    """Weighted in-transit depth of every covered event, in ppm."""
    offset = phase_offset(lc.time, epoch, period)
    cycle = np.rint((lc.time - epoch) / period).astype(int)
    weight = 1.0 / lc.error ** 2
    inside = np.abs(offset) <= duration / 2
    near = (np.abs(offset) >= 1.5 * duration) & (np.abs(offset) <= 3.0 * duration)
    events = []
    for k in np.unique(cycle[inside | near]):
        mine = cycle == k
        i, o = mine & inside, mine & near
        if i.sum() < min_in or o.sum() < min_out:
            continue
        w_in, w_out = weight[i].sum(), weight[o].sum()
        mean_in = (weight[i] * lc.flux[i]).sum() / w_in
        mean_out = (weight[o] * lc.flux[o]).sum() / w_out
        events.append({
            "cycle": int(k), "midpoint": float(epoch + k * period),
            "depth_ppm": float((mean_out - mean_in) / mean_out * 1e6),
            "sigma_ppm": float(math.sqrt(1 / w_in + 1 / w_out) / mean_out * 1e6),
            "points": int(i.sum()), "part": int(np.bincount(lc.part[i]).argmax()),
        })
    return events


def combine(events: list[dict]) -> dict:
    """Inverse-variance mean depth, its error, and the chi^2 consistency."""
    if not events:
        return {"depth_ppm": None, "sigma_ppm": None, "events": 0}
    d = np.array([e["depth_ppm"] for e in events])
    s = np.array([e["sigma_ppm"] for e in events])
    w = 1 / s ** 2
    mean = float((w * d).sum() / w.sum())
    chi2 = float((w * (d - mean) ** 2).sum())
    dof = len(events) - 1
    return {"depth_ppm": mean, "sigma_ppm": float(1 / math.sqrt(w.sum())), "events": len(events),
            "snr": mean * math.sqrt(w.sum()), "chi2": chi2, "dof": dof,
            "p_value": float(stats.chi2.sf(chi2, dof)) if dof > 0 else None}


def odd_even(events: list[dict]) -> dict:
    odd = combine([e for e in events if e["cycle"] % 2])
    even = combine([e for e in events if e["cycle"] % 2 == 0])
    if odd["depth_ppm"] is None or even["depth_ppm"] is None:
        return {"status": "NOT_EVALUATED", "reason": "needs both odd and even events"}
    sigma = math.hypot(odd["sigma_ppm"], even["sigma_ppm"])
    z = abs(odd["depth_ppm"] - even["depth_ppm"]) / sigma
    return {"status": "PASS" if z < 3 else "FAIL", "odd_ppm": odd["depth_ppm"],
            "even_ppm": even["depth_ppm"], "difference_sigma": z}


def part_consistency(lc: LightCurve, events: list[dict]) -> dict:
    """Depth per sector/product part and whether they agree (chi^2 p >= 0.01)."""
    rows = []
    for part in sorted({e["part"] for e in events}):
        summary = combine([e for e in events if e["part"] == part])
        meta = lc.parts[part] if lc.parts else {}
        rows.append({"part": part, "sector": meta.get("sector"), "author": meta.get("author"),
                     "events": summary["events"], "depth_ppm": summary["depth_ppm"],
                     "sigma_ppm": summary["sigma_ppm"]})
    if len(rows) < 2:
        return {"parts": rows, "p_value": None, "status": "NOT_EVALUATED",
                "reason": "events in fewer than two sectors/products"}
    d = np.array([r["depth_ppm"] for r in rows])
    w = 1 / np.array([r["sigma_ppm"] for r in rows]) ** 2
    mean = (w * d).sum() / w.sum()
    chi2 = float((w * (d - mean) ** 2).sum())
    p = float(stats.chi2.sf(chi2, len(rows) - 1))
    return {"parts": rows, "chi2": chi2, "dof": len(rows) - 1, "p_value": p,
            "status": "PASS" if p >= 0.01 else "FAIL"}


def scan_threshold(period: float, duration: float, false_alarm: float = 0.01) -> tuple[float, float]:
    """One-sided Gaussian threshold for the strongest of the scan's
    independent phase trials (one per duration), and that trial count."""
    trials = max((period - 3 * duration) / duration, 1.0)
    return float(stats.norm.isf(false_alarm / trials)), trials


def secondary_scan(lc: LightCurve, epoch: float, period: float, duration: float,
                   step_fraction: float = 1 / 3, false_alarm: float = 0.01,
                   half_sigma: float = 3.0) -> dict:
    """Eclipse-like dips at every phase away from the transit (eccentric
    secondaries need not sit at phase 0.5), modeled on the Kepler Robovetter
    model-shift test (Coughlin et al. 2016, ApJS 224, 12):

    * Per-event errors are white; the robust scatter of the SNR across all
      scanned phases measures how much correlated noise inflates it (the
      red-noise factor, floored at 1). SNRs are divided by it.
    * The strongest phase is compared with a threshold that holds the scan's
      false-alarm probability at `false_alarm` given ~P/T independent phases.
    * Phase 0.5 is a single pre-specified trial, tested at `half_sigma`.
    * A real eclipse recurs, so the strongest phase must stay above threshold
      with its most significant single event removed (jackknife), and must
      exceed the strongest positive excursion (bump) of the scan by the
      Robovetter FA2 margin, sqrt(2) erfcinv(T/P); systematics make bumps
      as readily as dips.
    """
    step = max(duration * step_fraction, 1e-4)
    shifts = np.arange(1.5 * duration, period - 1.5 * duration, step)
    scanned = []
    for shift in shifts:
        events = event_depths(lc, epoch + shift, period, duration)
        summary = combine(events)
        if summary["depth_ppm"] is not None and summary["events"] >= 1:
            scanned.append({"phase": float(shift / period), **summary, "_events": events})
    half_events = event_depths(lc, epoch + period / 2, period, duration)
    half = combine(half_events)
    threshold, trials = scan_threshold(period, duration, false_alarm)
    if not scanned:
        return {"status": "NOT_EVALUATED", "strongest": None, "phase_half": half,
                "phases_tested": int(len(shifts)), "threshold_sigma": threshold}
    snr = np.array([s["snr"] for s in scanned])
    best = max(scanned, key=lambda s: s["snr"])
    best_events = best.pop("_events")
    for s in scanned:
        s.pop("_events", None)
    # Phases whose local baseline (1.5-3 durations away) overlaps the transit
    # or the strongest dip show mirror-image "bumps"; leave them out of the
    # noise and bump statistics.
    shift = np.array([s["phase"] for s in scanned]) * period
    from_transit = np.minimum(shift, period - shift)
    gap = np.abs(shift - best["phase"] * period)
    from_best = np.minimum(gap, period - gap)
    clean = (from_transit >= 3.5 * duration) & (from_best >= 3.5 * duration)
    reference = snr[clean] if clean.sum() >= 5 else snr
    red = max(1.0, float(1.4826 * np.median(np.abs(reference - np.median(reference)))))
    best_sigma = best["snr"] / red
    jackknife_sigma = _jackknife_snr(best_events) / red
    bump_sigma = max(0.0, -float(reference.min())) / red
    margin = float(stats.norm.isf(min(duration / period, 0.5)))
    best_failed = (best_sigma >= threshold and jackknife_sigma >= threshold
                   and best_sigma - bump_sigma >= margin)
    half_snr = half.get("snr")
    half_sigma_value = None if half_snr is None else half_snr / red
    half_jackknife = _jackknife_snr(half_events) / red if half_events else None
    half_failed = (half_sigma_value is not None and half_sigma_value >= half_sigma
                   and half_jackknife is not None and half_jackknife >= half_sigma)
    return {"status": "FAIL" if best_failed or half_failed else "PASS", "strongest": best,
            "phase_half": half, "phases_tested": int(len(shifts)), "independent_trials": trials,
            "red_noise_factor": red, "threshold_sigma": threshold,
            "strongest_sigma_red": best_sigma, "strongest_jackknife_sigma_red": jackknife_sigma,
            "strongest_bump_sigma_red": bump_sigma, "bump_margin_sigma": margin,
            "phase_half_sigma_red": half_sigma_value, "phase_half_jackknife_sigma_red": half_jackknife,
            "half_threshold_sigma": half_sigma, "false_alarm": false_alarm,
            "method": ("model-shift style: red-noise-scaled SNR vs look-elsewhere threshold, "
                       "single-event jackknife, and positive-excursion margin")}


def _jackknife_snr(events: list[dict]) -> float:
    """Smallest combined SNR after dropping any one event (0 with < 2 events)."""
    if len(events) < 2:
        return 0.0
    return float(min(combine(events[:i] + events[i + 1:])["snr"] for i in range(len(events))))


def quality_overlap(lc: LightCurve, events: list[dict], duration: float) -> dict:
    """Events within a duration of a momentum dump or inside scattered light."""
    dumps = [t for meta in lc.parts for t in meta.get("momentum_dump_times", [])]
    spans = [s for meta in lc.parts for s in meta.get("scattered_light_intervals", [])]
    interpreted = any("momentum_dump_times" in meta for meta in lc.parts)
    flagged = []
    for e in events:
        near_dump = any(abs(t - e["midpoint"]) <= duration for t in dumps)
        in_light = any(lo - duration / 2 <= e["midpoint"] <= hi + duration / 2 for lo, hi in spans)
        if near_dump or in_light:
            flagged.append({"cycle": e["cycle"], "midpoint": e["midpoint"],
                            "momentum_dump": near_dump, "scattered_light": in_light})
    fraction = len(flagged) / len(events) if events else 0.0
    return {"status": ("NOT_EVALUATED" if not interpreted else
                       ("FAIL" if fraction > 0.5 else "PASS")),
            "flagged_events": flagged, "fraction": fraction,
            "note": "only SPOC-family quality bits are interpreted"}


def grazing(fit: dict) -> dict:
    """V-shaped (grazing) geometry: impact parameter beyond 1 - Rp/R*."""
    b, k = fit.get("impact"), fit.get("rp_rs")
    if b is None or k is None:
        return {"status": "NOT_EVALUATED"}
    is_grazing = b > 1 - k
    return {"status": "FAIL" if is_grazing else "PASS", "impact": b, "rp_rs": k,
            "note": "grazing geometry is common for eclipsing binaries" if is_grazing else ""}


def sector_coverage(lc: LightCurve, events: list[dict]) -> dict:
    sectors = sorted({str(lc.parts[e["part"]].get("sector")) if lc.parts else str(e["part"])
                      for e in events})
    authors = sorted({str(lc.parts[e["part"]].get("author")) for e in events}) if lc.parts else []
    return {"sectors_with_events": sectors, "authors_with_events": authors,
            "events_per_sector": {s: sum(1 for e in events if (str(lc.parts[e["part"]].get("sector"))
                                                              if lc.parts else str(e["part"])) == s)
                                  for s in sectors}}
