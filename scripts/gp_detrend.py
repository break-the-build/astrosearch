#!/usr/bin/env python3
"""Flatten a light curve with a Gaussian process for active (variable) stars.

Fits a celerite2 stochastically driven, damped harmonic oscillator (SHO) GP to
each segment, with a candidate's transits optionally masked so the GP does not
absorb them, and writes flux / GP prediction. Search the result with
`astrosearch --trend-method none`, and compare against the default biweight
with injection-recovery (e.g. scripts/real_curve_validation.py) before
preferring either for a target class.

usage: gp_detrend.py INPUT.csv OUTPUT.csv [--mask-period P --mask-epoch T0 --mask-duration D]
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vetting.fap import segments  # noqa: E402
from vetting.io import phase_offset, read_table  # noqa: E402


def fit_segment(t: np.ndarray, y: np.ndarray, e: np.ndarray, use: np.ndarray) -> np.ndarray:
    import celerite2
    from celerite2 import terms

    mean = float(np.median(y[use]))
    resid = y - mean
    scale = float(np.std(resid[use])) or 1e-4

    def gp_for(log_params):
        log_sigma, log_rho, log_jitter = log_params
        kernel = terms.SHOTerm(sigma=math.exp(log_sigma), rho=math.exp(log_rho), Q=1 / math.sqrt(2))
        gp = celerite2.GaussianProcess(kernel, mean=0.0)
        gp.compute(t[use], diag=e[use] ** 2 + math.exp(2 * log_jitter), quiet=True)
        return gp

    def negative_log_likelihood(log_params):
        try:
            return -gp_for(log_params).log_likelihood(resid[use])
        except Exception:
            return 1e25

    # Variability timescale bounded to >= 1 day so transit-length dips stay in the residuals.
    start = [math.log(scale), math.log(3.0), math.log(float(np.median(e)) * 0.1 + 1e-9)]
    bounds = [(math.log(scale) - 5, math.log(scale) + 5), (math.log(1.0), math.log(50.0)),
              (math.log(1e-7), math.log(scale) + 2)]
    best = minimize(negative_log_likelihood, start, method="L-BFGS-B", bounds=bounds)
    gp = gp_for(best.x)
    return mean + gp.predict(resid[use], t=t, return_cov=False)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--mask-period", type=float)
    p.add_argument("--mask-epoch", type=float)
    p.add_argument("--mask-duration", type=float)
    a = p.parse_args()
    data = read_table(a.input)
    t, y, e = data[:, 0], data[:, 1], data[:, 2]
    use = np.ones(len(t), dtype=bool)
    if a.mask_period and a.mask_epoch is not None and a.mask_duration:
        use = np.abs(phase_offset(t, a.mask_epoch, a.mask_period)) > a.mask_duration
    model = np.empty_like(y)
    for idx in segments(t):
        if use[idx].sum() < 20:
            model[idx] = np.median(y[idx])
            continue
        model[idx] = fit_segment(t[idx], y[idx], e[idx], use[idx])
    rows = np.column_stack([t, y / model, e / model])
    np.savetxt(a.output, rows, delimiter=",", header="time,flux,flux_err", comments="", fmt="%.8f")
    print(f"wrote {a.output}: {len(t)} rows, {len(segments(t))} segments")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
