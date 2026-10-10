#!/usr/bin/env python3
"""Calibrate astrosearch on real per-sector light curves.

Runs inverted-flux and per-segment circular-shift null searches, randomized
finite-exposure transit injections, and a leave-one-season-out ephemeris check.
The report is not planet validation or a calibrated survey-wide false-alarm rate.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
from pathlib import Path

import numpy as np
csv.field_size_limit(2**31 - 1)  # candidate rows carry long evidence/alias cells


def load_curve(path: Path) -> np.ndarray:
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=float, encoding=None)
    names = data.dtype.names or ()
    def col(*keys):
        key = next((x for x in keys if x in names), None)
        if key is None:
            raise ValueError(f"{path}: missing one of columns {keys}")
        return np.asarray(data[key], dtype=float)
    t, f = col("time", "btjd", "bjd"), col("flux", "pdcsap_flux", "sap_flux")
    try:
        e = col("flux_err", "error", "pdcsap_flux_err", "sap_flux_err")
    except ValueError:
        e = np.full(len(t), np.nan)
    good = np.isfinite(t) & np.isfinite(f) & (f > 0)
    if np.isfinite(e).any():
        good &= np.isfinite(e) & (e > 0)
    else:
        e = np.full(len(t), np.nan)
    exposure = col("exposure_days") if "exposure_days" in names else None
    out = np.column_stack((t[good], f[good], e[good],
                           exposure[good] if exposure is not None else np.full(int(good.sum()), np.nan)))
    out = out[np.argsort(out[:, 0], kind="stable")]
    if not np.isfinite(out[:, 3]).any():
        out[:, 3] = float(np.median(np.diff(out[:, 0])))
    else:
        invalid = ~np.isfinite(out[:, 3]) | (out[:, 3] <= 0)
        out[invalid, 3] = float(np.nanmedian(out[~invalid, 3]))
    return out


def noise_density(curve: np.ndarray, exposure_s: float) -> float:
    t, f = curve[:, 0], curve[:, 1]
    if len(t) < 30 or exposure_s <= 0:
        return math.nan
    dt, df = np.diff(t), np.diff(f)
    positive = dt[np.isfinite(dt) & (dt > 0)]
    if not len(positive):
        return math.nan
    cadence = float(np.median(positive))
    d = df[np.isfinite(df) & np.isfinite(dt) & (dt <= 3 * cadence)]
    if len(d) < 20:
        return math.nan
    return 1.4826 * float(np.median(np.abs(d - np.median(d)))) / math.sqrt(2) * math.sqrt(exposure_s)


def write_curve(path: Path, curve: np.ndarray) -> None:
    np.savetxt(path, curve[:, :3], delimiter=",", header="time,flux,flux_err", comments="", fmt="%.15g")


def run_search(binary: Path, curve_path: Path, output: Path, args, pmin=None, pmax=None) -> None:
    cmd = [str(binary), "--input", str(curve_path), "--output", str(output),
           "--min-period", str(pmin if pmin is not None else args.min_period),
           "--max-period", str(pmax if pmax is not None else args.max_period),
           "--threads", str(args.threads), "--transit-candidates", str(args.candidates),
           "--trend-method", str(args.trend_method)]
    if args.stellar_density is not None:
        cmd += ["--stellar-density", str(args.stellar_density)]
    if args.period_trials:
        cmd += ["--periods", str(args.period_trials)]
    with (output.with_suffix(".log")).open("w") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)


def read_candidates(path: Path):
    if not path.exists():
        return []
    with path.open(newline="") as f:
        return [r for r in csv.DictReader(f) if r.get("class") in ("periodic_dip", "duo_transit")]


def phase_distance(a: float, b: float, period: float) -> float:
    d = abs((a - b) % period)
    return min(d, period - d)


def matching_candidate(rows, period, epoch, duration, span):
    for rank, row in enumerate(rows, 1):
        try:
            p, e, d = float(row["period_days"]), float(row["epoch"]), float(row["duration_days"])
        except (KeyError, ValueError):
            continue
        if abs(p - period) * span / period <= duration and phase_distance(e, epoch, period) <= duration:
            return {**row, "candidate_rank": rank}
    return None


def quadratic_transit_flux(separation, radius_ratio, u1, u2):
    """Numerically integrate a quadratic-limb-darkened stellar disk occultation."""
    z = np.asarray(separation, dtype=float).reshape(-1)
    nodes, weights = np.polynomial.legendre.leggauss(128)
    radius = (nodes + 1) / 2
    radial_weight = weights / 2
    mu = np.sqrt(np.maximum(0.0, 1 - radius * radius))
    intensity = 1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2
    annulus = radial_weight * radius * intensity
    total = float(annulus.sum())
    z2 = z[:, None]
    r2 = radius[None, :]
    denom = 2 * r2 * z2
    with np.errstate(divide="ignore", invalid="ignore"):
        cosine = (r2 * r2 + z2 * z2 - radius_ratio * radius_ratio) / denom
    fraction = np.where(cosine <= -1, 1.0,
                        np.where(cosine >= 1, 0.0, np.arccos(np.clip(cosine, -1, 1)) / np.pi))
    central = z < 1e-12
    if central.any():
        fraction[central] = (radius[None, :] <= radius_ratio).astype(float)
    blocked = fraction @ annulus
    return 1 - blocked / total


def radius_ratio_for_depth(depth, impact_parameter, u1, u2):
    """Choose Rp/R* so the model's mid-transit depth matches the requested depth."""
    target = depth * 1e-6
    lo, hi = 1e-5, 0.8
    for _ in range(32):
        mid = 0.5 * (lo + hi)
        model_depth = 1 - float(quadratic_transit_flux(
            np.array([impact_parameter]), mid, u1, u2)[0])
        if model_depth < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def trapezoid_shape(offset, duration, ingress_fraction):
    """Unit-depth transit profile; ingress_fraction is each ramp / total duration."""
    half = 0.5 * duration
    ingress = max(0.5 * duration * ingress_fraction, np.finfo(float).eps)
    distance = np.abs(offset)
    return np.clip((half - distance) / ingress, 0.0, 1.0)


def inject(curve, period, phase, duration, depth_ppm, shape="quadratic",
           ingress_fraction=0.15, integration_samples=7, impact_parameter=0.4,
           limb_darkening=(0.4, 0.2), exposure_days=None):
    out = curve.copy()
    t0 = float(curve[0, 0] + phase * period)
    centers = np.arange(t0 - period, curve[-1, 0] + period, period)
    if exposure_days is None:
        exposure = np.full(len(out), float(np.median(np.diff(out[:, 0]))))
    else:
        exposure = np.asarray(exposure_days, dtype=float)
    u1, u2 = limb_darkening
    radius_ratio = (radius_ratio_for_depth(depth_ppm, impact_parameter, u1, u2)
                    if shape == "quadratic" else None)
    offsets = (np.arange(integration_samples) + 0.5) / integration_samples - 0.5
    for center in centers:
        near = np.abs(out[:, 0] - center) <= 0.5 * duration + exposure
        if not np.any(near):
            continue
        if shape == "box":
            profile = (np.abs(out[near, 0] - center) <= 0.5 * duration).astype(float)
            model_flux = 1 - depth_ppm * 1e-6 * profile
        elif shape == "trapezoid":
            sample_times = (out[near, 0, None]
                            + exposure[near, None] * offsets[None, :])
            profile = trapezoid_shape(sample_times - center, duration, ingress_fraction)
            model_flux = 1 - depth_ppm * 1e-6 * profile.mean(axis=1)
        else:
            velocity = math.sqrt(max((1 + radius_ratio) ** 2 - impact_parameter ** 2, 1e-12)) / (duration / 2)
            sample_times = (out[near, 0, None]
                            + exposure[near, None] * offsets[None, :])
            separation = np.sqrt(impact_parameter ** 2
                                 + (velocity * (sample_times - center)) ** 2)
            model_flux = quadratic_transit_flux(separation.ravel(), radius_ratio, u1, u2)
            model_flux = model_flux.reshape(sample_times.shape).mean(axis=1)
        out[near, 1] *= model_flux
    return out, t0


def covered_transits(curve, period, epoch, duration, exposure_days):
    """Count transit epochs with at least two in-transit samples (search eligible)."""
    times = curve[:, 0]
    first = math.ceil((times[0] - epoch) / period)
    last = math.floor((times[-1] - epoch) / period)
    covered = 0
    expected = max(0, last - first + 1)
    for cycle in range(first, last + 1):
        center = epoch + cycle * period
        mask = np.abs(times - center) <= duration / 2 + exposure_days / 2
        covered += int(mask.sum() >= 2)
    return expected, covered


def circular_shift_null(curve, rng):
    """Independently circular-shift each continuous segment, retaining local noise."""
    out = curve.copy()
    t = out[:, 0]
    delta = np.diff(t)
    positive = delta[np.isfinite(delta) & (delta > 0)]
    if not len(positive):
        return out
    cadence = float(np.median(positive))
    breaks = np.flatnonzero(delta > max(5.0 * cadence, 0.25)) + 1
    bounds = np.r_[0, breaks, len(out)]
    for start, stop in zip(bounds[:-1], bounds[1:]):
        n = stop - start
        if n > 2:
            shift = int(rng.integers(1, n))
            out[start:stop, 1] = np.roll(out[start:stop, 1], shift)
            out[start:stop, 2] = np.roll(out[start:stop, 2], shift)
    return out


def wilson_interval(successes, trials, z=1.96):
    if trials == 0:
        return None, None
    p = successes / trials
    den = 1 + z * z / trials
    center = (p + z * z / (2 * trials)) / den
    half = z * math.sqrt(p * (1 - p) / trials + z * z / (4 * trials * trials)) / den
    return max(0.0, center - half), min(1.0, center + half)


def heldout_depth(curve, period, epoch, duration):
    t, f, err = curve[:, :3].T
    estimates = []
    first = math.ceil((t[0] - epoch) / period)
    last = math.floor((t[-1] - epoch) / period)
    for cycle in range(first, last + 1):
        center = epoch + cycle * period
        inside = np.abs(t - center) <= duration / 2
        outside = (np.abs(t - center) >= 1.5 * duration) & (np.abs(t - center) <= 3.0 * duration)
        if inside.sum() < 2 or outside.sum() < 4:
            continue
        if np.isfinite(err).all():
            wi, wo = 1 / err[inside]**2, 1 / err[outside]**2
            transit = np.average(f[inside], weights=wi)
            baseline = np.average(f[outside], weights=wo)
            sigma = math.sqrt(1 / wi.sum() + 1 / wo.sum())
        else:
            transit, baseline = float(np.median(f[inside])), float(np.median(f[outside]))
            sigma = 1.4826 * float(np.median(np.abs(f[outside] - baseline))) / math.sqrt(inside.sum())
        if sigma > 0 and math.isfinite(sigma):
            estimates.append(((baseline - transit), sigma))
    if not estimates:
        return {"transits_covered": 0, "depth": None, "depth_sigma": None}
    weights = np.array([1 / s**2 for _, s in estimates])
    depths = np.array([d for d, _ in estimates])
    return {"transits_covered": len(estimates),
            "depth": float(np.average(depths, weights=weights)),
            "depth_sigma": float(np.sum(depths * weights) / math.sqrt(weights.sum()))}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sector-dir", required=True, type=Path, help="Folder containing sector_NNN.csv light curves")
    p.add_argument("--output", required=True, type=Path, help="Validation output directory")
    p.add_argument("--binary", type=Path, help="astrosearch executable (default target/release/astrosearch)")
    p.add_argument("--min-period", type=float, default=1.0)
    p.add_argument("--max-period", type=float, default=15.0)
    p.add_argument("--stellar-density", type=float,
                   help="optional target mean density in g/cm^3 for the adaptive transit grid")
    p.add_argument("--max-sector-noise-ratio", type=float, default=5.0)
    p.add_argument("--trend-method", choices=("median", "biweight"), default="biweight")
    p.add_argument("--depths-ppm", default="100,200,400,800")
    p.add_argument("--injections-per-depth", type=int, default=5,
                   help="Random realizations per depth; use at least 50 per bin for useful completeness estimates")
    p.add_argument("--duration-hours", default="1,2,4")
    p.add_argument("--injection-shape", choices=("quadratic", "trapezoid", "box"), default="quadratic",
                   help="quadratic limb-darkened transit is default; trapezoid/box are controls")
    p.add_argument("--ingress-fraction", type=float, default=0.15,
                   help="each trapezoid ingress/egress as a fraction of total duration")
    p.add_argument("--max-impact-parameter", type=float, default=0.9,
                   help="uniform upper bound for injected impact parameter")
    p.add_argument("--limb-darkening", default="0.4,0.2",
                   help="quadratic TESS-band u1,u2; replace with stellar-model coefficients")
    p.add_argument("--null-trials", type=int, default=10,
                   help="independent per-segment circular-shift null searches")
    p.add_argument("--period-trials", type=int, default=0)
    p.add_argument("--candidates", type=int, default=256)
    p.add_argument("--threads", type=int, default=0)
    p.add_argument("--seed", type=int, default=20261009)
    p.add_argument("--max-candidates", type=int, default=10)
    p.add_argument("--fold-trials", type=int, default=1200)
    args = p.parse_args()
    if not (0 < args.min_period < args.max_period and args.injections_per_depth > 0):
        p.error("invalid period range or injection count")
    if (args.null_trials < 0 or not (0.01 <= args.ingress_fraction <= 0.49)
            or not (0.0 <= args.max_impact_parameter < 1.0)):
        p.error("null-trials must be nonnegative, ingress-fraction in [0.01, 0.49], and max-impact-parameter in [0, 1)")
    try:
        limb_darkening = tuple(float(value) for value in args.limb_darkening.split(","))
    except ValueError:
        p.error("--limb-darkening must be two comma-separated coefficients")
    if len(limb_darkening) != 2 or not all(math.isfinite(value) for value in limb_darkening):
        p.error("--limb-darkening must be two finite coefficients")
    if args.stellar_density is not None and (not math.isfinite(args.stellar_density) or args.stellar_density <= 0):
        p.error("--stellar-density must be finite and positive")
    if args.max_sector_noise_ratio < 0 or 0 < args.max_sector_noise_ratio < 1:
        p.error("--max-sector-noise-ratio must be 0 (disabled) or at least 1")
    root = Path(__file__).resolve().parents[1]
    binary = args.binary.resolve() if args.binary else root / "target/release/astrosearch"
    if not binary.is_file():
        p.error(f"astrosearch binary not found: {binary}; run cargo build --release")
    args.output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    files = sorted(args.sector_dir.glob("sector_*.csv"))
    if not files:
        p.error("no sector_*.csv files found")
    products = []
    for path in files:
        curve = load_curve(path)
        try:
            sector = int(path.stem.split("_")[-1])
        except ValueError:
            sector = len(products) + 1
        exposure = 86400 * float(np.median(curve[:, 3]))
        products.append({"sector": sector, "path": path, "curve": curve,
                         "noise_density": noise_density(curve, exposure),
                         "exposure_s": exposure})
    valid = [x["noise_density"] for x in products if math.isfinite(x["noise_density"]) and x["noise_density"] > 0]
    reference = float(np.median(valid)) if valid else math.nan
    for x in products:
        x["noise_ratio"] = x["noise_density"] / reference if math.isfinite(reference) and math.isfinite(x["noise_density"]) else None
        x["included"] = x["noise_ratio"] is None or args.max_sector_noise_ratio <= 0 or x["noise_ratio"] <= args.max_sector_noise_ratio
    included = [x for x in products if x["included"]]
    product_arrays = [x["curve"] for x in included]
    with_exposure = np.concatenate(product_arrays)
    with_exposure = with_exposure[np.argsort(with_exposure[:, 0], kind="stable")]
    with_exposure = with_exposure[np.r_[True, np.diff(with_exposure[:, 0]) > 0]]
    exposure_days = with_exposure[:, 3]
    curve = with_exposure[:, :3]
    if len(curve) < 50:
        p.error("fewer than 50 usable observations after sector-noise filtering")
    combined_path = args.output / "combined.csv"
    write_curve(combined_path, curve)
    qc_rows = [{"sector": x["sector"], "observations": len(x["curve"]),
                "exposure_s": x["exposure_s"], "noise_density": x["noise_density"],
                "noise_ratio": x["noise_ratio"], "included": x["included"]} for x in products]
    with (args.output / "sector_quality.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=qc_rows[0].keys()); w.writeheader(); w.writerows(qc_rows)

    base_path = args.output / "base_candidates.csv"
    run_search(binary, combined_path, base_path, args)

    flipped = curve.copy()
    flipped[:, 1] = 2 * float(np.median(flipped[:, 1])) - flipped[:, 1]
    inverted_path = args.output / "inverted.csv"
    write_curve(inverted_path, flipped)
    inverted_candidates = args.output / "inverted_candidates.csv"
    run_search(binary, inverted_path, inverted_candidates, args)
    null_scores = [float(r["score"]) for r in read_candidates(inverted_candidates) if r.get("score")]

    shifted_null_maxima = []
    for trial in range(args.null_trials):
        null_curve = circular_shift_null(curve, rng)
        null_input = args.output / f"shift_null_{trial + 1:03d}.csv"
        null_output = args.output / f"shift_null_{trial + 1:03d}_candidates.csv"
        write_curve(null_input, null_curve)
        run_search(binary, null_input, null_output, args)
        scores = [float(r["score"]) for r in read_candidates(null_output)
                  if r.get("score") and math.isfinite(float(r["score"]))]
        shifted_null_maxima.append(max(scores, default=0.0))

    depths = [float(x) for x in args.depths_ppm.split(",") if x.strip()]
    durations = [float(x) / 24 for x in args.duration_hours.split(",") if x.strip()]
    injections = []
    span = float(curve[-1, 0] - curve[0, 0])
    for depth in depths:
        for trial in range(args.injections_per_depth):
            period = float(np.exp(rng.uniform(math.log(args.min_period), math.log(args.max_period))))
            phase = float(rng.random())
            duration = float(rng.choice(durations))
            impact = float(rng.uniform(0, args.max_impact_parameter))
            injected, epoch = inject(curve, period, phase, duration, depth,
                                     shape=args.injection_shape,
                                     ingress_fraction=args.ingress_fraction,
                                     impact_parameter=impact,
                                     limb_darkening=limb_darkening,
                                     exposure_days=exposure_days)
            expected_events, covered_events = covered_transits(
                curve, period, epoch, duration, exposure_days)
            inp = args.output / "injection_input.csv"
            out = args.output / "injection_candidates.csv"
            write_curve(inp, injected)
            run_search(binary, inp, out, args)
            hit = matching_candidate(read_candidates(out), period, epoch, duration, span)
            injections.append({"depth_ppm": depth, "trial": trial, "period_days": period,
                "phase_fraction": phase, "duration_days": duration,
                "injection_shape": args.injection_shape,
                "impact_parameter": impact,
                "radius_ratio": (radius_ratio_for_depth(depth, impact, *limb_darkening)
                                 if args.injection_shape == "quadratic" else None),
                "expected_transit_epochs": expected_events,
                "covered_transit_epochs": covered_events,
                "search_coverage_eligible": covered_events >= 2,
                "recovered": hit is not None,
                "recovered_rank": hit.get("candidate_rank") if hit else None,
                "recovered_period_days": hit.get("period_days") if hit else None,
                "recovered_score": hit.get("score") if hit else None})
    with (args.output / "injection_recovery.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=injections[0].keys()); w.writeheader(); w.writerows(injections)

    # Discover candidates without one observing year, then test those
    # training-only ephemerides in the withheld year's measurements. This
    # avoids choosing candidates on the same data used for the prediction.
    folds = []
    season = np.floor((curve[:, 0] - curve[0, 0]) / 365.25).astype(int)
    for year in sorted(set(season)):
        held = season == year
        train = curve[~held]
        test = curve[held]
        if len(train) < 50 or len(test) < 10:
            continue
        train_in = args.output / "fold_train.csv"
        train_out = args.output / "fold_candidates.csv"
        write_curve(train_in, train)
        fold_args = argparse.Namespace(**vars(args))
        fold_args.period_trials = args.fold_trials
        run_search(binary, train_in, train_out, fold_args)
        for rank, candidate in enumerate(read_candidates(train_out)[:args.max_candidates], 1):
            try:
                pfit = float(candidate["period_days"])
                efit = float(candidate["epoch"])
                duration = float(candidate["duration_days"])
            except (KeyError, ValueError):
                continue
            result = heldout_depth(test, pfit, efit, duration)
            folds.append({"training_candidate_rank": rank, "training_period_days": pfit,
                "training_epoch": efit, "duration_days": duration,
                "held_out_season": int(year), **result})
    if folds:
        with (args.output / "holdout_seasons.csv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=folds[0].keys()); w.writeheader(); w.writerows(folds)

    by_depth = {}
    for depth in depths:
        group = [x for x in injections if x["depth_ppm"] == depth]
        recovered = sum(x["recovered"] for x in group)
        low, high = wilson_interval(recovered, len(group))
        ranks = [x["recovered_rank"] for x in group if x["recovered_rank"] is not None]
        eligible = [x for x in group if x["search_coverage_eligible"]]
        eligible_recovered = sum(x["recovered"] for x in eligible)
        eligible_low, eligible_high = wilson_interval(eligible_recovered, len(eligible))
        by_depth[str(depth)] = {
            "recovered_any_rank": recovered,
            "trials": len(group),
            "recovery_fraction_any_rank": recovered / len(group),
            "wilson_95_interval_any_rank": [low, high],
            "coverage_eligible_trials_at_least_two_epochs": len(eligible),
            "recovered_among_coverage_eligible": eligible_recovered,
            "recovery_fraction_given_two_epoch_coverage": (
                eligible_recovered / len(eligible) if eligible else None),
            "wilson_95_interval_given_two_epoch_coverage": [eligible_low, eligible_high],
            "trials_with_insufficient_coverage": len(group) - len(eligible),
            "recovered_top_1": sum(rank <= 1 for rank in ranks),
            "recovered_top_5": sum(rank <= 5 for rank in ranks),
            "recovered_top_10": sum(rank <= 10 for rank in ranks),
        }
    real_max = max((float(r["score"]) for r in read_candidates(base_path)
                    if r.get("score") and math.isfinite(float(r["score"]))), default=0.0)
    null_exceedances = sum(score >= real_max for score in shifted_null_maxima)
    empirical_tail = ((null_exceedances + 1) / (len(shifted_null_maxima) + 1)
                      if shifted_null_maxima else None)
    report = {"trend_method": args.trend_method,
        "input_sectors": qc_rows, "included_sectors": [x["sector"] for x in included],
        "excluded_sectors": [x["sector"] for x in products if not x["included"]],
        "sector_noise_threshold": args.max_sector_noise_ratio,
        "inverted_null": {"candidate_count": len(null_scores), "maximum_score": max(null_scores, default=0.0),
                          "scores": null_scores},
        "circular_shift_null": {
            "trials": len(shifted_null_maxima),
            "maximum_scores": shifted_null_maxima,
            "minimum_resolvable_plus_one_tail_fraction": (
                1 / (len(shifted_null_maxima) + 1) if shifted_null_maxima else None),
            "real_maximum_score": real_max,
            "null_maxima_at_least_real_maximum": null_exceedances,
            "plus_one_empirical_tail_fraction": empirical_tail,
            "warning": "Exploratory max-score comparison only; not a calibrated false-alarm probability."},
        "injection_recovery": by_depth,
        "injection_model": {"shape": args.injection_shape,
                            "ingress_fraction": args.ingress_fraction,
                            "limb_darkening_quadratic_u1_u2": limb_darkening,
                            "impact_parameter_distribution": [0.0, args.max_impact_parameter],
                            "exposure_source": "per-sector median cadence from the product timestamps",
                            "exposure_integration_samples": 7,
                            "warning": "Quadratic limb-darkened circular-orbit model with linear projected motion; supply stellar-model TESS-band coefficients and extend the geometry priors for publication-grade completeness."},
        "holdout_season_folds": len(folds),
        "warning": "Exploratory calibration only: few realizations have wide intervals. Results are target-specific and are not survey completeness, a calibrated FAP, or planet validation."}
    (args.output / "validation_report.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
