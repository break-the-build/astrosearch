#!/usr/bin/env python3
"""Sector-level TESS target-pixel in/out-of-transit difference diagnostic.

Downloads one target-pixel file, builds the difference image, and reports the
centroid of its significant pixels relative to the target's proper-motion-
corrected position, with an event-bootstrap uncertainty (vetting.pixel).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from vetting import network, pixel, stellar

WARNING = ("Single-sector image diagnostic. Pixel S/N assumes independent cadences; systematics, "
           "crowding and few transits can bias the centroid. Inspect the images before relying on it.")


def run_pixel_check(tic_id: int, sector: int, period: float, epoch_btjd: float, duration_days: float,
                    ra_deg: float, dec_deg: float, output: Path, author: str = "SPOC",
                    pm_ra: float | None = None, pm_dec: float | None = None) -> dict:
    """Run the diagnostic and write pixel_vet.json plus the images to `output`."""
    import astropy.units as u
    import lightkurve as lk
    from astropy.coordinates import SkyCoord

    result = lk.search_targetpixelfile(f"TIC {tic_id}", mission="TESS", sector=sector, author=author)
    if len(result) == 0:
        raise RuntimeError(f"no {author} target-pixel file for TIC {tic_id} sector {sector}")
    tpf = result[[0]].download(quality_bitmask="hard")
    if tpf is None:
        raise RuntimeError("MAST did not return a downloadable target-pixel product")
    time = np.asarray(tpf.time.value, dtype=float)
    quality = np.asarray(tpf.quality, dtype=int)
    cube = np.ma.filled(np.asarray(tpf.flux.value), np.nan).astype(float)
    valid = (quality == 0) & np.isfinite(time) & np.isfinite(cube).all(axis=(1, 2))
    time, cube = time[valid], cube[valid]
    measured = pixel.centroid_with_uncertainty(time, cube, epoch_btjd, period, duration_days)
    year = stellar.btjd_to_year(float(np.median(time)))
    ra_epoch, dec_epoch = stellar.propagate_position(ra_deg, dec_deg, pm_ra, pm_dec,
                                                     stellar.TIC_EPOCH_YEAR, year)
    target = SkyCoord(ra_epoch * u.deg, dec_epoch * u.deg)
    target_xy = [float(v) for v in tpf.wcs.world_to_pixel(target)]
    report = {"tic_id": tic_id, "sector": sector, "author": author, "period_days": period,
              "epoch_btjd": epoch_btjd, "duration_days": duration_days,
              "usable_cadences": int(len(time)), "target_pixel_xy": target_xy,
              "target_epoch_year": year, "proper_motion_applied": pm_ra is not None and pm_dec is not None,
              "warning": WARNING}
    output.mkdir(parents=True, exist_ok=True)
    if measured["status"] != "OK":
        report.update(status=measured["status"], reason=measured.get("reason"))
    else:
        np.savetxt(output / "difference_image.csv", measured.pop("difference"), delimiter=",")
        np.savetxt(output / "difference_snr.csv", measured.pop("snr"), delimiter=",")
        offset = pixel.offset_significance(measured["centroid_xy"], measured["centroid_sigma_xy"], target_xy)
        sky = tpf.wcs.pixel_to_world(*measured["centroid_xy"])
        report.update(measured)
        report.update(centroid_separation_arcsec=float(target.separation(sky).to_value(u.arcsec)),
                      centroid_offset_sigma=offset["offset_sigma"],
                      pixel_difference_csv=str(output / "difference_image.csv"),
                      pixel_difference_snr_csv=str(output / "difference_snr.csv"))
        report["status"] = pixel.centroid_status(report)
    (output / "pixel_vet.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tic-id", required=True, type=int)
    parser.add_argument("--sector", required=True, type=int)
    parser.add_argument("--period", required=True, type=float)
    parser.add_argument("--epoch-btjd", required=True, type=float)
    parser.add_argument("--duration-days", required=True, type=float)
    parser.add_argument("--ra-deg", required=True, type=float)
    parser.add_argument("--dec-deg", required=True, type=float)
    parser.add_argument("--pm-ra", type=float, help="mas/yr including cos(dec) (TIC pmRA)")
    parser.add_argument("--pm-dec", type=float, help="mas/yr (TIC pmDEC)")
    parser.add_argument("--author", default="SPOC")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    network.install_default_timeouts()
    if args.period <= 0 or args.duration_days <= 0:
        parser.error("period and duration must be positive")
    try:
        report = run_pixel_check(args.tic_id, args.sector, args.period, args.epoch_btjd,
                                 args.duration_days, args.ra_deg, args.dec_deg, args.output,
                                 args.author, args.pm_ra, args.pm_dec)
    except (RuntimeError, ImportError) as error:
        raise SystemExit(str(error)) from error
    print(args.output / "pixel_vet.json", report.get("status"))


if __name__ == "__main__":
    main()
