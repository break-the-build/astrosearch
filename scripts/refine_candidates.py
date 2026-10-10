#!/usr/bin/env python3
"""Refit periodic-dip triggers with a finite-exposure quadratic transit model.

This is a candidate-refinement tool, not an astrophysical validation or FPP
calculation. It fits one candidate at a time with circular-orbit geometry and
fixed quadratic limb-darkening coefficients. Supply a light curve CSV and the
candidate CSV emitted by astrosearch.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.optimize import least_squares
csv.field_size_limit(2**31 - 1)  # candidate rows carry long evidence/alias cells


_GL_X, _GL_W = np.polynomial.legendre.leggauss(64)


def occultation_flux(separation, radius_ratio, u1, u2):
    """Quadratic-limb-darkened occultation with adaptive radial intervals.

    Integration intervals follow the exact full/partial occultation boundaries,
    avoiding the fixed radial grid's loss of accuracy for small planets.
    """
    z = np.asarray(separation, dtype=float).reshape(-1)
    p = float(radius_ratio)
    total_r = (_GL_X + 1) / 2
    total_w = _GL_W / 2
    total_mu = np.sqrt(np.maximum(0.0, 1 - total_r * total_r))
    total_i = 1 - u1 * (1 - total_mu) - u2 * (1 - total_mu) ** 2
    total = float(np.sum(total_w * total_r * total_i))
    blocked = np.zeros(len(z))
    full_hi = np.clip(p - z, 0.0, 1.0)
    full_r = full_hi[:, None] * total_r[None, :]
    full_mu = np.sqrt(np.maximum(0.0, 1 - full_r * full_r))
    full_i = 1 - u1 * (1 - full_mu) - u2 * (1 - full_mu) ** 2
    blocked += full_hi * full_hi * 0.5 * np.sum(
        total_w[None, :] * total_r[None, :] * full_i, axis=1)
    lo = np.abs(p - z)
    hi = np.minimum(1.0, p + z)
    width = np.maximum(0.0, hi - lo)
    r = lo[:, None] + width[:, None] * total_r[None, :]
    mu = np.sqrt(np.maximum(0.0, 1 - r * r))
    intensity = 1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2
    denom = 2 * r * z[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        cosine = (r * r + z[:, None] * z[:, None] - p * p) / denom
    fraction = np.where(width[:, None] > 0,
                        np.arccos(np.clip(cosine, -1, 1)) / np.pi, 0.0)
    partial = width * np.sum(total_w[None, :] * r * intensity * fraction, axis=1)
    blocked += partial
    # At z=0 the partial interval collapses and the full-disc integral above
    # correctly integrates the occulted central disk.
    return 1 - blocked / total


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"{path}: no rows")
    names = set(rows[0])

    def column(*keys, required=True, default=None):
        key = next((item for item in keys if item in names), None)
        if key is None:
            if required:
                raise ValueError(f"{path}: missing column (expected one of {keys})")
            return np.full(len(rows), default, dtype=object)
        return np.array([row.get(key, "") for row in rows], dtype=object)

    def numbers(*keys, required=True, default=np.nan):
        values = column(*keys, required=required, default=default)
        return np.array([float(value) if str(value).strip() else default for value in values], dtype=float)

    time = numbers("time", "btjd", "bjd")
    flux = numbers("flux", "pdcsap_flux", "sap_flux")
    error = numbers("flux_err", "error", "pdcsap_flux_err", "sap_flux_err", required=False)
    exposure = numbers("exposure_days", required=False)
    sector = column("sector", required=False, default="0").astype(str)
    good = np.isfinite(time) & np.isfinite(flux) & (flux > 0)
    time, flux, error, exposure, sector = (v[good] for v in (time, flux, error, exposure, sector))
    order = np.argsort(time, kind="stable")
    time, flux, error, exposure, sector = (v[order] for v in (time, flux, error, exposure, sector))
    unique = np.r_[True, np.diff(time) > 0]
    time, flux, error, exposure, sector = (v[unique] for v in (time, flux, error, exposure, sector))
    if len(time) < 20:
        raise ValueError("need at least 20 unique finite observations")
    normalization = float(np.nanmedian(flux))
    if not math.isfinite(normalization) or normalization <= 0:
        raise ValueError("could not determine a positive flux normalization")
    flux = flux / normalization
    error = error / normalization
    sidecar = path.with_suffix(".meta.json")
    column_given = np.isfinite(exposure).any()
    if not column_given and sidecar.exists():
        # export_tess.py sidecars give each sector's exposure; mixed-cadence
        # fits otherwise smear 30-minute samples as if they were 2-minute.
        meta = json.loads(sidecar.read_text())
        for part in meta.get("parts", [meta]):
            lo, hi, seconds = part.get("time_min"), part.get("time_max"), part.get("exptime_s")
            if lo is not None and hi is not None and seconds:
                exposure[(time >= lo - 1e-9) & (time <= hi + 1e-9)] = float(seconds) / 86400.0
    if np.isfinite(exposure).all():
        exposure_source = ("exposure_days column" if column_given
                           else "per-sector exptime from the .meta.json sidecar")
    elif not np.isfinite(exposure).any():
        exposure[:] = np.median(np.diff(time))
        exposure_source = "estimated from median timestamp spacing (supply exposure_days for mixed cadence)"
    else:
        fallback = np.nanmedian(exposure[np.isfinite(exposure) & (exposure > 0)])
        exposure[~np.isfinite(exposure) | (exposure <= 0)] = fallback
        exposure_source = "exposure_days column"
    if not np.isfinite(error).any():
        error[:] = 1.4826 * np.median(np.abs(flux - np.median(flux)))
        if not np.isfinite(error[0]) or error[0] <= 0:
            error[:] = np.std(flux)
        error_source = "global robust scatter fallback"
    else:
        missing = ~np.isfinite(error) | (error <= 0)
        fallback = np.nanmedian(error[~missing])
        error[missing] = fallback
        error_source = "input uncertainty column with median fill for invalid rows"
    if np.any(~np.isfinite(error)) or np.any(error <= 0) or np.any(exposure <= 0):
        raise ValueError("could not determine positive uncertainties and exposures")
    return time, flux, error, exposure, sector, exposure_source, error_source


def model_factory(time, exposure, u1, u2, samples):
    groups = []
    # Keep per-product integrations separate for mixed-cadence light curves.
    keys = np.round(exposure, 10)
    for value in np.unique(keys):
        indices = np.flatnonzero(keys == value)
        groups.append((indices, float(value)))

    try:
        import batman
    except ImportError:
        batman = None
    if batman is not None:
        params = batman.TransitParams()
        params.t0 = float(time[0])
        params.per = 2.0
        params.rp = 0.05
        params.a = 10.0
        params.inc = 88.0
        params.ecc = 0.0
        params.w = 90.0
        params.u = [u1, u2]
        params.limb_dark = "quadratic"
        models = [(indices, batman.TransitModel(
            params, time[indices], supersample_factor=samples, exp_time=exp_time))
            for indices, exp_time in groups]

        def evaluate_batman(t0, period, radius_ratio, a_rs, impact):
            if not (a_rs > 1 + radius_ratio and impact < 1 + radius_ratio):
                return np.ones(len(time))
            params.t0, params.per, params.rp, params.a = t0, period, radius_ratio, a_rs
            params.inc = math.degrees(math.acos(min(1.0, max(-1.0, impact / a_rs))))
            output = np.ones(len(time))
            for indices, transit_model in models:
                output[indices] = transit_model.light_curve(params)
            return output

        return evaluate_batman, "batman-quadratic-finite-exposure"

    def evaluate(t0, period, radius_ratio, a_rs, impact):
        if not (a_rs > 1 + radius_ratio and impact < 1 + radius_ratio):
            return np.ones(len(time))
        inclination_cosine = min(1.0, max(-1.0, impact / a_rs))
        output = np.ones(len(time))
        offsets = ((np.arange(samples) + 0.5) / samples - 0.5)
        for indices, exp_time in groups:
            sample_times = time[indices, None] + exp_time * offsets[None, :]
            angle = 2 * np.pi * (sample_times - t0) / period
            separation = a_rs * np.sqrt(
                np.sin(angle) ** 2 + inclination_cosine ** 2 * np.cos(angle) ** 2)
            output[indices] = occultation_flux(separation.ravel(), radius_ratio, u1, u2).reshape(
                sample_times.shape).mean(axis=1)
        return output

    return evaluate, "custom-quadratic-radial-integration-finite-exposure"


def vet_diagnostics(time, flux, error, sector, epoch, period, duration):
    """Report inexpensive secondary-eclipse and per-sector depth checks."""
    phase = (time - epoch + period / 2) % period - period / 2
    transit_distance = np.abs(phase)
    secondary_distance = period / 2 - transit_distance
    sector_depths = []
    for label in np.unique(sector):
        in_sector = sector == label
        inside = in_sector & (transit_distance <= duration / 2)
        outside = in_sector & (transit_distance >= 1.5 * duration) & (transit_distance <= 3 * duration)
        if inside.sum() >= 2 and outside.sum() >= 4:
            baseline = float(np.median(flux[outside]))
            measured = float(baseline - np.median(flux[inside])) * 1e6
            sigma = math.sqrt(float(np.median(error[inside]) ** 2 / inside.sum()
                                    + np.median(error[outside]) ** 2 / outside.sum())) * 1e6
            sector_depths.append({"sector": str(label), "depth_ppm": measured,
                                  "depth_sigma_ppm_approx": sigma,
                                  "in_transit_points": int(inside.sum()),
                                  "baseline_points": int(outside.sum())})
    weights = np.array([1 / max(row["depth_sigma_ppm_approx"], 1e-9) ** 2
                        for row in sector_depths])
    depth_values = np.array([row["depth_ppm"] for row in sector_depths])
    if len(sector_depths) >= 2:
        weighted_mean = float(np.average(depth_values, weights=weights))
        consistency_chi2 = float(np.sum(weights * (depth_values - weighted_mean) ** 2))
        consistency_dof = len(sector_depths) - 1
    else:
        consistency_chi2, consistency_dof = None, None

    secondary_supported = 6 * duration < period
    secondary = ((secondary_distance <= duration / 2)
                 & (secondary_distance >= 0))
    secondary_out = ((secondary_distance >= 1.5 * duration)
                     & (secondary_distance <= 3 * duration))
    if secondary_supported and secondary.sum() >= 2 and secondary_out.sum() >= 4:
        baseline = float(np.median(flux[secondary_out]))
        depth = float(baseline - np.median(flux[secondary])) * 1e6
        sigma = math.sqrt(float(np.median(error[secondary]) ** 2 / secondary.sum()
                                + np.median(error[secondary_out]) ** 2 / secondary_out.sum())) * 1e6
        secondary_result = {"depth_ppm": depth, "depth_sigma_ppm_approx": sigma,
                            "significance_approx": depth / sigma if sigma > 0 else None,
                            "in_eclipse_points": int(secondary.sum())}
    else:
        secondary_result = {"depth_ppm": None, "depth_sigma_ppm_approx": None,
                            "significance_approx": None,
                            "in_eclipse_points": int(secondary.sum())}
    return {"secondary_eclipse_phase_0p5": secondary_result,
            "sector_transit_depths": sector_depths,
            "sector_depth_consistency_chi2": consistency_chi2,
            "sector_depth_consistency_dof": consistency_dof}


#: Above this many observed transits the fit uses a sparse finite-difference Jacobian.
SPARSE_JACOBIAN_EVENTS = 40


def event_jacobian_sparsity(event_index, n_events):
    """Sparsity of d(residual)/d(baseline, slope) per event: sample i depends
    only on its own event's two nuisance parameters."""
    rows = np.repeat(np.arange(len(event_index)), 2)
    columns = np.ravel(np.column_stack([2 * event_index, 2 * event_index + 1]))
    return sparse.csr_matrix((np.ones(len(rows)), (rows, columns)),
                             shape=(len(event_index), 2 * n_events))


def fit_candidate(time, flux, error, exposure, sector, row, args):
    period0, epoch0 = float(row["period_days"]), float(row["epoch"])
    duration0 = max(float(row["duration_days"]), np.median(exposure) * 2)
    depth_ppm = max(float(row.get("depth_ppm") or 300), 1.0)
    if period0 <= 0 or not math.isfinite(period0 + epoch0 + duration0):
        raise ValueError("candidate has invalid ephemeris or duration")
    span = float(time[-1] - time[0])
    epoch_radius = max(duration0, args.epoch_window_durations * duration0)
    epoch_number = np.rint((time - epoch0) / period0).astype(int)
    phase_offset = time - (epoch0 + epoch_number * period0)
    keep = np.abs(phase_offset) <= args.fit_window_durations * duration0
    indices = np.flatnonzero(keep)
    if len(indices) < 20:
        raise ValueError("fewer than 20 observations near predicted transit windows")
    t, y, sigma, exp, labels, cycle = (v[indices] for v in (time, flux, error, exposure, sector, epoch_number))
    event_keys = np.array([f"{s}:{n}" for s, n in zip(labels, cycle)])
    unique_events = list(dict.fromkeys(event_keys.tolist()))
    if len(unique_events) < 2:
        raise ValueError("fewer than two observed transit epochs; physical fit is underconstrained")
    position = {key: i for i, key in enumerate(unique_events)}
    event_index = np.array([position[key] for key in event_keys])
    # Initialize a/R* from the small-planet, circular duration approximation.
    b0 = min(0.6, max(0.05, args.initial_impact))
    a0 = max(2.0, period0 / (math.pi * duration0) * math.sqrt(max(1 - b0 * b0, 0.05)))
    radius0 = min(0.3, max(0.001, math.sqrt(depth_ppm * 1e-6)))
    u1, u2 = args.limb_darkening
    transit, model_backend = model_factory(t, exp, u1, u2, args.integration_samples)
    event_center = np.array([epoch0 + int(key.rsplit(":", 1)[1]) * period0 for key in unique_events])
    # Shared parameters: t0, period, Rp/R*, a/R*, impact parameter.
    # Each observed transit gets a nuisance baseline and slope in both null
    # and transit fits so slow local trends do not masquerade as depth.
    n_events = len(unique_events)
    p0 = np.r_[epoch0, period0, radius0, a0, b0, np.tile([1.0, 0.0], n_events)]
    period_half_width = max(duration0 / max(span, period0) * period0, period0 * args.period_window_fraction)
    lower = np.r_[epoch0 - epoch_radius, period0 - period_half_width, 1e-4, 1.1, 0.0,
                  np.tile([0.8, -0.1], n_events)]
    upper = np.r_[epoch0 + epoch_radius, period0 + period_half_width, 0.5, 1000.0, 0.99,
                  np.tile([1.2, 0.1], n_events)]
    p0 = np.minimum(np.maximum(p0, lower + 1e-9), upper - 1e-9)
    centered_time = (t - event_center[event_index]) / max(duration0, 1e-8)

    def residual(parameters, with_transit):
        if with_transit:
            baseline_parameters = parameters[5:]
            profile = transit(*parameters[:5])
        else:
            baseline_parameters = parameters
            profile = np.ones(len(t))
        baseline = (baseline_parameters[0::2][event_index]
                    + baseline_parameters[1::2][event_index] * centered_time)
        return (y - baseline * profile) / sigma

    # Each event's baseline and slope touch only that event's samples, so the
    # Jacobian is sparse; declaring it lets finite differences perturb all
    # baseline columns at once (about 7 model evaluations per Jacobian instead
    # of 5 + 2 x events, which made multi-year fits take tens of minutes).
    # With few events the exact dense solver is fast and reaches a marginally
    # better optimum than the iterative sparse one, so sparsity is used only
    # when there are many events (measured: 203 events, 3482 s -> 181 s, and
    # the dense fit had not converged; 28 events, 59 s either way).
    transit_sparsity = null_sparsity = None
    if n_events > SPARSE_JACOBIAN_EVENTS:
        null_sparsity = event_jacobian_sparsity(event_index, n_events)
        transit_sparsity = sparse.hstack([sparse.csr_matrix(np.ones((len(t), 5))),
                                          null_sparsity]).tocsr()
    result = least_squares(lambda x: residual(x, True), p0, bounds=(lower, upper),
                           loss="soft_l1", f_scale=args.robust_scale, max_nfev=args.max_evaluations,
                           jac_sparsity=transit_sparsity)
    null = least_squares(lambda x: residual(x, False), p0[5:],
                         bounds=(lower[5:], upper[5:]), loss="soft_l1",
                         f_scale=args.robust_scale, max_nfev=args.max_evaluations,
                         jac_sparsity=null_sparsity)
    physical_resid = residual(result.x, True)
    # The null model has only per-event baselines; pass its own parameters.
    null_resid = residual(null.x, False)
    chi2, null_chi2 = float(physical_resid @ physical_resid), float(null_resid @ null_resid)
    delta = max(0.0, null_chi2 - chi2)
    k_transit, k_null, n = len(result.x), len(null.x), len(t)
    bic_transit = chi2 + k_transit * math.log(n)
    bic_null = null_chi2 + k_null * math.log(n)
    errors = [None] * 5
    epoch_period_covariance = None
    try:
        jacobian = result.jac.toarray() if sparse.issparse(result.jac) else result.jac
        covariance = np.linalg.pinv(jacobian.T @ jacobian) * chi2 / max(1, n - k_transit)
        errors = np.sqrt(np.maximum(0, np.diag(covariance)[:5])).tolist()
        epoch_period_covariance = float(covariance[0, 1])
    except np.linalg.LinAlgError:
        pass
    posterior = {}
    if args.mcmc_steps > 0 and result.success:
        try:
            import emcee
        except ImportError as exception:
            raise RuntimeError("install optional inference dependencies with `uv sync --extra inference`") from exception

        def log_probability(parameters):
            if np.any(parameters <= lower) or np.any(parameters >= upper):
                return -np.inf
            values = residual(parameters, True)
            return float(-0.5 * values @ values)

        nwalkers = max(32, 2 * len(result.x) + 2)
        scales = np.r_[max(duration0 * 0.02, 1e-6),
                       max(period_half_width * 0.02, 1e-8),
                       max(result.x[2] * 0.02, 1e-6),
                       max(result.x[3] * 0.02, 1e-5), 0.01,
                       np.full(len(result.x) - 5, 1e-4)]
        rng = np.random.default_rng(args.mcmc_seed)
        initial = result.x + rng.normal(size=(nwalkers, len(result.x))) * scales
        initial = np.minimum(np.maximum(initial, lower + 1e-10), upper - 1e-10)
        sampler = emcee.EnsembleSampler(nwalkers, len(result.x), log_probability)
        sampler.run_mcmc(initial, args.mcmc_steps, progress=False)
        flat = sampler.get_chain(discard=args.mcmc_steps // 2, flat=True)
        quantiles = np.quantile(flat[:, :5], [0.16, 0.5, 0.84], axis=0)
        posterior = {
            "mcmc_median_t0_period_rprs_a_rs_impact": quantiles[1].tolist(),
            "mcmc_16_84_bounds": [[float(v) for v in quantiles[0]],
                                  [float(v) for v in quantiles[2]]],
            "mcmc_steps_per_walker": args.mcmc_steps,
            "mcmc_walkers": nwalkers,
            "mcmc_mean_acceptance_fraction": float(np.mean(sampler.acceptance_fraction)),
            "mcmc_seed": args.mcmc_seed,
            "mcmc_convergence_note": "posterior sample summary only; inspect chain mixing and autocorrelation before reporting uncertainties",
        }
    fitted_t14 = (float(result.x[1]) / math.pi
                  * math.asin(min(1.0, math.sqrt(max((1 + result.x[2]) ** 2 - result.x[4] ** 2, 0.0))
                                                  / max(result.x[3] * math.sqrt(max(1 - (result.x[4] / result.x[3]) ** 2, 1e-12)), 1e-12))))
    vet = vet_diagnostics(time, flux, error, sector, float(result.x[0]),
                          float(result.x[1]), fitted_t14)
    return {
        "candidate_class_rank": row.get("class_rank", ""),
        "input_period_days": period0, "input_epoch_btjd": epoch0,
        "fitted_period_days": float(result.x[1]), "fitted_period_sigma_days_local": errors[1],
        "fitted_epoch_btjd": float(result.x[0]), "fitted_epoch_sigma_days_local": errors[0],
        "fitted_radius_ratio": float(result.x[2]), "fitted_radius_ratio_sigma_local": errors[2],
        "fitted_a_over_rstar": float(result.x[3]), "fitted_impact_parameter": float(result.x[4]),
        "fitted_a_over_rstar_sigma_local": errors[3],
        "fitted_impact_parameter_sigma_local": errors[4],
        "fitted_epoch_period_covariance_local": epoch_period_covariance,
        "fitted_duration_days_circular_approx": fitted_t14,
        "observed_transit_windows": n_events, "fit_observations": n,
        "delta_chi2_vs_local_baseline": delta, "delta_bic_null_minus_transit": bic_null - bic_transit,
        "optimizer_success": bool(result.success), "optimizer_message": result.message,
        **posterior,
        **{key: json.dumps(value, sort_keys=True) if isinstance(value, (list, dict)) else value
           for key, value in vet.items()},
        "model": model_backend,
        "interpretation": "screening fit only; local covariance omits ephemeris aliases, stellar priors, correlated noise, and model uncertainty",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lightcurve", required=True, type=Path)
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--limb-darkening", default="0.4,0.2")
    parser.add_argument("--integration-samples", type=int, default=9)
    parser.add_argument("--initial-impact", type=float, default=0.4)
    parser.add_argument("--epoch-window-durations", type=float, default=1.0)
    parser.add_argument("--fit-window-durations", type=float, default=2.5)
    parser.add_argument("--period-window-fraction", type=float, default=0.01)
    parser.add_argument("--robust-scale", type=float, default=1.0)
    parser.add_argument("--max-evaluations", type=int, default=3000)
    parser.add_argument("--mcmc-steps", type=int, default=0,
                        help="optional emcee steps per walker; 0 keeps deterministic optimization only")
    parser.add_argument("--mcmc-seed", type=int, default=20261009)
    parser.add_argument("--max-candidates", type=int, default=12)
    args = parser.parse_args()
    try:
        args.limb_darkening = tuple(float(x) for x in args.limb_darkening.split(","))
    except ValueError:
        parser.error("--limb-darkening must be u1,u2")
    if (len(args.limb_darkening) != 2 or not all(math.isfinite(x) for x in args.limb_darkening)
            or args.integration_samples < 1 or args.max_candidates < 1
            or args.fit_window_durations <= 1 or args.period_window_fraction <= 0
            or args.mcmc_steps < 0):
        parser.error("invalid limb-darkening coefficients, fit windows, or iteration settings")
    time, flux, error, exposure, sector, exposure_source, error_source = read_csv(args.lightcurve)
    with args.candidates.open(newline="", encoding="utf-8-sig") as stream:
        candidates = [row for row in csv.DictReader(stream) if row.get("class") in ("periodic_dip", "duo_transit")]
    results = []
    for row in candidates[:args.max_candidates]:
        try:
            results.append({"candidate_rank": row.get("rank", ""),
                            "candidate_class": row.get("class", ""),
                            **fit_candidate(time, flux, error, exposure, sector, row, args)})
        except Exception as exception:
            results.append({"candidate_rank": row.get("rank", ""),
                            "candidate_class": row.get("class", ""),
                            "candidate_class_rank": row.get("class_rank", ""),
                            "input_period_days": row.get("period_days", ""),
                            "input_epoch_btjd": row.get("epoch", ""),
                            "fit_error": f"{type(exception).__name__}: {exception}"})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        fields = sorted({key for row in results for key in row})
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows(results)
    metadata = {"lightcurve": str(args.lightcurve), "candidates": str(args.candidates),
                "exposure_source": exposure_source, "uncertainty_source": error_source,
                "limb_darkening_u1_u2": args.limb_darkening,
                "model": sorted({row["model"] for row in results if row.get("model")}),
                "warning": "Screening fit only; not planet validation, FPP, or a final transit solution.",
                "input_candidates": len(candidates), "fits_written": len(results)}
    args.output.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"wrote {args.output} and {args.output.with_suffix('.json')}")


if __name__ == "__main__":
    main()
