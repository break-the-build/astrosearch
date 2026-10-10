#!/usr/bin/env python3
"""Iterative limb-darkened Transit Least Squares on a normalized light curve."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


def read_curve(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"{path}: no light-curve rows")
    names = set(rows[0])

    def values(*keys, required=True, default=np.nan):
        key = next((name for name in keys if name in names), None)
        if key is None:
            if required:
                raise ValueError(f"{path}: missing column from {keys}")
            return np.full(len(rows), default)
        return np.asarray([float(row[key]) if str(row.get(key, "")).strip() else default
                           for row in rows], dtype=float)

    time = values("time", "btjd", "bjd")
    flux = values("flux", "pdcsap_flux", "sap_flux")
    error = values("flux_err", "error", "pdcsap_flux_err", "sap_flux_err", required=False)
    sector_key = next((name for name in ("sector", "sector_number") if name in names), None)
    sector = (np.asarray([row[sector_key] for row in rows], dtype=str)
              if sector_key else np.full(len(rows), "all"))
    good = np.isfinite(time) & np.isfinite(flux) & (flux > 0)
    if np.isfinite(error).any():
        good &= np.isfinite(error) & (error > 0)
    time, flux, error, sector = (array[good] for array in (time, flux, error, sector))
    order = np.argsort(time, kind="stable")
    time, flux, error, sector = (array[order] for array in (time, flux, error, sector))
    if len(time) < 50:
        raise ValueError("need at least 50 valid observations")
    # Remove sector constants while preserving within-sector transit shapes.
    for label in np.unique(sector):
        in_sector = sector == label
        center = float(np.median(flux[in_sector]))
        flux[in_sector] /= center
        if np.isfinite(error[in_sector]).any():
            error[in_sector] /= center
    invalid = ~np.isfinite(error) | (error <= 0)
    if invalid.any():
        valid_errors = error[~invalid]
        fallback = (float(np.median(valid_errors)) if valid_errors.size else
                    max(1.4826 * float(np.median(np.abs(flux - np.median(flux)))), 1e-6))
        error[invalid] = fallback
    return time, flux, error


def alias_label(period: float, prior: list[dict]) -> str:
    for signal in prior:
        ratio = period / signal["period_days"]
        nearest = round(2 * ratio) / 2
        if nearest > 0 and abs(ratio - nearest) / nearest <= 0.01:
            return "integer" if nearest.is_integer() else "half-integer"
    return "none"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lightcurve", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--min-period", type=float, default=0.15)
    parser.add_argument("--max-period", type=float, default=30.0)
    parser.add_argument("--max-planets", type=int, default=5)
    parser.add_argument("--min-sde", type=float, default=7.0)
    parser.add_argument("--depth-min-ppm", type=float, default=30.0)
    parser.add_argument("--oversampling", type=int, default=3)
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--stellar-radius", type=float,
                        help="stellar radius in solar units; improves the TLS duration grid")
    parser.add_argument("--stellar-mass", type=float,
                        help="stellar mass in solar units; improves the TLS duration grid")
    parser.add_argument("--limb-darkening", default="0.4,0.2",
                        help="quadratic TESS-band u1,u2 coefficients")
    args = parser.parse_args()
    if not (0 < args.min_period < args.max_period and args.max_planets > 0
            and args.min_sde >= 0 and args.depth_min_ppm > 0
            and args.oversampling >= 1 and args.threads >= 0):
        parser.error("invalid search configuration")
    try:
        limb_darkening = tuple(float(value) for value in args.limb_darkening.split(","))
    except ValueError:
        parser.error("--limb-darkening must be two comma-separated numbers")
    if len(limb_darkening) != 2 or not all(math.isfinite(value) for value in limb_darkening):
        parser.error("--limb-darkening must be two finite coefficients")
    if any(value is not None and (not math.isfinite(value) or value <= 0)
           for value in (args.stellar_radius, args.stellar_mass)):
        parser.error("stellar radius and mass must be finite and positive")
    try:
        from transitleastsquares import model as tls_model, transit_mask
    except ImportError as error:
        parser.error(f"install TLS support with uv sync --extra tls: {error}")

    time, flux, error = read_curve(args.lightcurve)
    keep = np.ones(len(time), dtype=bool)
    output_rows = []
    for iteration in range(args.max_planets):
        if keep.sum() < 50:
            break
        search = tls_model(time[keep], flux[keep], dy=error[keep], verbose=False)
        kwargs = {
            "period_min": args.min_period,
            "period_max": min(args.max_period, (time[keep][-1] - time[keep][0]) / 2),
            "n_transits_min": 2,
            "transit_depth_min": args.depth_min_ppm * 1e-6,
            "oversampling_factor": args.oversampling,
            "transit_template": "default",
            "u": list(limb_darkening),
            "limb_dark": "quadratic",
            "show_progress_bar": False,
            "verbose": False,
        }
        if args.stellar_radius is not None:
            kwargs["R_star"] = args.stellar_radius
        if args.stellar_mass is not None:
            kwargs["M_star"] = args.stellar_mass
        if kwargs["period_max"] <= kwargs["period_min"]:
            break
        if args.threads:
            kwargs["use_threads"] = args.threads
        result = search.power(**kwargs)
        sde = float(result.SDE)
        if not math.isfinite(sde) or sde < args.min_sde:
            break
        fap = getattr(result, "FAP", math.nan)
        row = {
            "iteration": iteration + 1, "period_days": float(result.period),
            "epoch_days": float(result.T0), "duration_days": float(result.duration),
            "depth_ppm": float(result.depth) * 1e6, "radius_ratio": float(result.rp_rs),
            "sde": sde, "sde_raw": float(result.SDE_raw),
            "fap_white_noise": float(fap) if fap is not None and math.isfinite(float(fap)) else "",
            "odd_even_mismatch_sigma": float(result.odd_even_mismatch)
            if result.odd_even_mismatch is not None else "",
            "transit_count": int(result.transit_count),
            "distinct_transit_count": int(result.distinct_transit_count),
            "alias_type_vs_previous": alias_label(float(result.period), output_rows),
            "engine": "transitleastsquares-default-limb-darkened-template",
        }
        output_rows.append(row)
        in_transit = transit_mask(time, row["period_days"], row["duration_days"], row["epoch_days"])
        if not np.any(in_transit & keep):
            break
        keep &= ~in_transit

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(output_rows[0]) if output_rows else [
        "iteration", "period_days", "epoch_days", "duration_days", "depth_ppm",
        "radius_ratio", "sde", "sde_raw", "fap_white_noise", "odd_even_mismatch_sigma",
        "transit_count", "distinct_transit_count", "alias_type_vs_previous", "engine"]
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output_rows)
    args.output.with_suffix(".json").write_text(json.dumps({
        "input": str(args.lightcurve), "signals": len(output_rows),
        "remaining_observations": int(keep.sum()), "minimum_sde": args.min_sde,
        "warning": "TLS FAP assumes white noise; use inverted-curve controls to assess systematics.",
    }, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.output}; iterative TLS signals: {len(output_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
