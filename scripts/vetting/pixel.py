"""Target-pixel difference imaging: where on the detector does the flux drop?

The in-transit minus out-of-transit image of an on-target transit peaks at the
target; a blended eclipsing binary shifts it toward the contaminant. Three
details matter for honest results:

* the target position must be propagated with proper motion to the TPF epoch
  (TIC positions are J2000; LHS 3844 moves ~1.3 arcsec/yr, ~1 pixel by 2019);
* the centroid uses only pixels whose difference is significant, so noise in
  the stamp's empty pixels does not drag it toward the stamp center;
* the offset is judged against its own uncertainty, estimated by
  resampling transit events, not against a fixed one-pixel cut alone.
"""
from __future__ import annotations

import numpy as np

TESS_PIXEL_ARCSEC = 21.0


def event_windows(time: np.ndarray, epoch: float, period: float, duration: float):
    offset = np.abs((time - epoch + period / 2) % period - period / 2)
    cycle = np.rint((time - epoch) / period).astype(int)
    return offset <= duration / 2, offset >= 2.0 * duration, cycle


def difference_image(cube: np.ndarray, inside: np.ndarray, outside: np.ndarray):
    in_image = np.nanmedian(cube[inside], axis=0)
    out_image = np.nanmedian(cube[outside], axis=0)
    difference = out_image - in_image
    scatter = 1.4826 * np.nanmedian(np.abs(cube[outside] - out_image), axis=0)
    uncertainty = scatter * np.sqrt(np.pi / 2) * np.sqrt(1.0 / inside.sum() + 1.0 / outside.sum())
    snr = np.divide(difference, uncertainty, out=np.zeros_like(difference),
                    where=np.isfinite(uncertainty) & (uncertainty > 0))
    return difference, snr


def significant_centroid(difference: np.ndarray, snr: np.ndarray, min_snr: float = 3.0):
    """Flux-weighted centroid (x, y) of pixels with SNR >= min_snr; None if none."""
    mask = np.isfinite(difference) & (snr >= min_snr) & (difference > 0)
    if not mask.any():
        return None
    weights = np.where(mask, difference, 0.0)
    yy, xx = np.indices(difference.shape)
    total = weights.sum()
    return float((xx * weights).sum() / total), float((yy * weights).sum() / total)


def centroid_with_uncertainty(time: np.ndarray, cube: np.ndarray, epoch: float, period: float,
                              duration: float, draws: int = 200, seed: int = 20261010,
                              min_snr: float = 3.0) -> dict:
    """Centroid of the significant difference pixels and its scatter when
    whole transit events are resampled with replacement."""
    inside, outside, cycle = event_windows(time, epoch, period, duration)
    if inside.sum() < 3 or outside.sum() < 10:
        return {"status": "NOT_EVALUATED",
                "reason": f"insufficient cadences: in={int(inside.sum())}, out={int(outside.sum())}"}
    difference, snr = difference_image(cube, inside, outside)
    centroid = significant_centroid(difference, snr, min_snr)
    if centroid is None:
        return {"status": "NOT_EVALUATED", "reason": f"no difference pixel with SNR >= {min_snr}",
                "max_pixel_snr": float(np.nanmax(snr))}
    events = np.unique(cycle[inside])
    rng = np.random.default_rng(seed)
    samples = []
    for _ in range(draws if len(events) >= 2 else 0):
        chosen = rng.choice(events, size=len(events), replace=True)
        resampled = np.concatenate([np.flatnonzero(inside & (cycle == e)) for e in chosen])
        boot_in = np.zeros_like(inside)
        boot_in[resampled] = True
        boot = significant_centroid(*difference_image(cube, boot_in, outside), min_snr)
        if boot is not None:
            samples.append(boot)
    sigma = (float(np.std([s[0] for s in samples], ddof=1)), float(np.std([s[1] for s in samples], ddof=1))) \
        if len(samples) >= 10 else (None, None)
    return {"status": "OK", "centroid_xy": list(centroid), "centroid_sigma_xy": list(sigma),
            "bootstrap_draws": len(samples), "events": int(len(events)),
            "in_transit_cadences": int(inside.sum()), "out_of_transit_cadences": int(outside.sum()),
            "max_pixel_snr": float(np.nanmax(snr)), "difference": difference, "snr": snr}


def offset_significance(centroid_xy, sigma_xy, target_xy) -> dict:
    """Offset (pixels and arcsec) and its significance; sigma floors at 0.05 px."""
    dx, dy = centroid_xy[0] - target_xy[0], centroid_xy[1] - target_xy[1]
    pixels = float(np.hypot(dx, dy))
    if sigma_xy[0] is None:
        return {"offset_pixels": pixels, "offset_arcsec": pixels * TESS_PIXEL_ARCSEC,
                "offset_sigma": None}
    sx, sy = max(sigma_xy[0], 0.05), max(sigma_xy[1], 0.05)
    significance = float(np.hypot(dx / sx, dy / sy))
    return {"offset_pixels": pixels, "offset_arcsec": pixels * TESS_PIXEL_ARCSEC,
            "offset_sigma": significance}


def centroid_status(report: dict, max_arcsec: float = TESS_PIXEL_ARCSEC,
                    min_sigma: float = 3.0) -> str:
    """FAIL only for an offset that is both larger than a pixel and >= 3 sigma."""
    arcsec = report.get("centroid_separation_arcsec")
    sigma = report.get("centroid_offset_sigma")
    if arcsec is None or not np.isfinite(arcsec):
        return "NOT_EVALUATED"
    if arcsec > max_arcsec and (sigma is None or sigma >= min_sigma):
        return "FAIL"
    return "PASS"
