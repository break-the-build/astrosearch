#!/usr/bin/env python3
"""Make a sector-level TESS target-pixel in/out transit difference diagnostic."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tic-id", required=True, type=int)
    parser.add_argument("--sector", required=True, type=int)
    parser.add_argument("--period", required=True, type=float)
    parser.add_argument("--epoch-btjd", required=True, type=float)
    parser.add_argument("--duration-days", required=True, type=float)
    parser.add_argument("--ra-deg", required=True, type=float)
    parser.add_argument("--dec-deg", required=True, type=float)
    parser.add_argument("--author", default="SPOC")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.period <= 0 or args.duration_days <= 0:
        parser.error("period and duration must be positive")
    try:
        import astropy.units as u
        import lightkurve as lk
        from astropy.coordinates import SkyCoord
    except ImportError as error:
        raise SystemExit("Install Python dependencies with `uv sync`") from error

    result = lk.search_targetpixelfile(f"TIC {args.tic_id}", mission="TESS",
        sector=args.sector, author=args.author)
    if len(result) == 0:
        raise SystemExit("No matching TESS target-pixel product was found")
    tpf = result[[0]].download(quality_bitmask="hard")
    if tpf is None:
        raise SystemExit("MAST did not return a downloadable target-pixel product")
    time = np.asarray(tpf.time.value, dtype=float)
    quality = np.asarray(tpf.quality, dtype=int)
    cube = np.ma.filled(np.asarray(tpf.flux.value), np.nan).astype(float)
    valid = (quality == 0) & np.isfinite(time) & np.isfinite(cube).all(axis=(1, 2))
    time, cube = time[valid], cube[valid]
    phase_distance = np.abs((time - args.epoch_btjd + args.period / 2) % args.period - args.period / 2)
    in_transit = phase_distance <= args.duration_days / 2
    out_of_transit = phase_distance >= 2.0 * args.duration_days
    if in_transit.sum() < 3 or out_of_transit.sum() < 10:
        raise SystemExit(f"Insufficient usable cadences: in={in_transit.sum()}, out={out_of_transit.sum()}")
    in_image = np.nanmedian(cube[in_transit], axis=0)
    out_image = np.nanmedian(cube[out_of_transit], axis=0)
    difference = out_image - in_image
    scatter = 1.4826 * np.nanmedian(np.abs(cube[out_of_transit] - out_image), axis=0)
    uncertainty = scatter * np.sqrt(np.pi / 2) * np.sqrt(
        1.0 / in_transit.sum() + 1.0 / out_of_transit.sum())
    snr_image = np.divide(difference, uncertainty, out=np.zeros_like(difference),
                          where=np.isfinite(uncertainty) & (uncertainty > 0))
    weights = np.where(np.isfinite(difference), np.maximum(difference, 0), 0)
    yy, xx = np.indices(weights.shape)
    weight_sum = float(weights.sum())
    if weight_sum <= 0:
        raise SystemExit("Difference image has no positive in-transit flux deficit")
    centroid_x = float((xx * weights).sum() / weight_sum)
    centroid_y = float((yy * weights).sum() / weight_sum)
    target = SkyCoord(args.ra_deg * u.deg, args.dec_deg * u.deg)
    target_x, target_y = tpf.wcs.world_to_pixel(target)
    difference_sky = tpf.wcs.pixel_to_world(centroid_x, centroid_y)
    separation_arcsec = float(target.separation(difference_sky).to_value(u.arcsec))
    args.output.mkdir(parents=True, exist_ok=True)
    np.savetxt(args.output / "difference_image.csv", difference, delimiter=",")
    np.savetxt(args.output / "difference_snr.csv", snr_image, delimiter=",")
    report = {
        "tic_id": args.tic_id, "sector": args.sector, "author": args.author,
        "period_days": args.period, "epoch_btjd": args.epoch_btjd,
        "duration_days": args.duration_days, "usable_cadences": int(len(time)),
        "in_transit_cadences": int(in_transit.sum()),
        "out_of_transit_cadences": int(out_of_transit.sum()),
        "target_pixel_xy": [float(target_x), float(target_y)],
        "positive_difference_centroid_xy": [centroid_x, centroid_y],
        "centroid_separation_arcsec": separation_arcsec,
        "pixel_difference_csv": str(args.output / "difference_image.csv"),
        "pixel_difference_snr_csv": str(args.output / "difference_snr.csv"),
        "warning": "Exploratory single-sector image diagnostic. Pixel S/N assumes independent cadences; systematics, dilution, aperture choice, and limited transit counts can invalidate the centroid estimate. Inspect images and centroid uncertainties manually.",
    }
    (args.output / "pixel_vet.json").write_text(json.dumps(report, indent=2) + "\n")
    print(args.output / "pixel_vet.json")


if __name__ == "__main__":
    main()
