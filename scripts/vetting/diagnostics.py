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


def secondary_scan(lc: LightCurve, epoch: float, period: float, duration: float,
                   step_fraction: float = 1 / 3) -> dict:
    """Eclipse-like dips at every phase away from the transit (eccentric
    secondaries need not sit at phase 0.5). Reports the most significant
    phase and phase 0.5 separately; the scan's own look-elsewhere effect
    means a single-phase 3 sigma is not significant, so the scan uses 4 sigma.
    """
    step = max(duration * step_fraction, 1e-4)
    shifts = np.arange(1.5 * duration, period - 1.5 * duration, step)
    best = None
    for shift in shifts:
        summary = combine(event_depths(lc, epoch + shift, period, duration))
        if summary["depth_ppm"] is None or summary["events"] < 1:
            continue
        if best is None or summary["snr"] > best["snr"]:
            best = {"phase": float(shift / period), **summary}
    half = combine(event_depths(lc, epoch + period / 2, period, duration))
    best_snr = best["snr"] if best else None
    half_snr = half.get("snr")
    failed = (best_snr is not None and best_snr >= 4) or (half_snr is not None and half_snr >= 3)
    return {"status": "FAIL" if failed else ("PASS" if best else "NOT_EVALUATED"),
            "strongest": best, "phase_half": half, "phases_tested": int(len(shifts))}


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
