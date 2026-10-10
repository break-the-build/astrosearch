#!/usr/bin/env python3
"""Measure the secondary-eclipse scan on known-answer benchmark planets.

For each confirmed planet, the scan runs at the catalog ephemeris on the raw
light curve. Two rules are compared on the same scan output:

* baseline: raw white-noise SNR >= 4 at any phase, or >= 3 at phase 0.5
  (the rule before the review-board calibration);
* calibrated: SNR divided by the red-noise factor, against the look-elsewhere
  threshold (vetting.diagnostics.secondary_scan).

Specificity: confirmed planets should pass (a hot-Jupiter secondary is
possible but rare at TESS precision). Power: the same light curves get an
injected phase-0.5 box with half the planet's depth proxy (an eclipsing
binary-like secondary), which should fail.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from vetting import diagnostics  # noqa: E402
from vetting.io import LightCurve, phase_offset, read_table  # noqa: E402

csv.field_size_limit(2**31 - 1)


def load(path: Path) -> LightCurve:
    data = read_table(path)
    data = data[np.argsort(data[:, 0], kind="stable")]
    error = data[:, 2]
    if not np.all(np.isfinite(error) & (error > 0)):
        fallback = np.nanmedian(np.abs(np.diff(data[:, 1]))) / np.sqrt(2) * 1.4826
        error = np.where(np.isfinite(error) & (error > 0), error, fallback)
    return LightCurve(data[:, 0], data[:, 1], error, np.zeros(len(data), dtype=int), ())


def rules(scan: dict) -> dict:
    best = (scan.get("strongest") or {}).get("snr")
    half = (scan.get("phase_half") or {}).get("snr")
    baseline = (best is not None and best >= 4) or (half is not None and half >= 3)
    return {"baseline_fail": bool(baseline), "calibrated_fail": scan["status"] == "FAIL",
            "raw_best_snr": best, "raw_half_snr": half,
            "red_noise_factor": scan.get("red_noise_factor"),
            "threshold_sigma": scan.get("threshold_sigma")}


def evaluate(job: tuple[dict, str]) -> dict:
    case, root = job
    lc = load(Path(root) / "lightcurves" / f"{case['case_id']}.csv")
    period, epoch = float(case["period_days"]), float(case["epoch_days"])
    duration = float(case["duration_days"])
    depth = float(case.get("depth_ppm_proxy") or 0) * 1e-6
    out = {"case_id": case["case_id"], "split": case["split"], "depth_stratum": case["depth_stratum"],
           "period": period, "duration": duration}
    out.update({f"clean_{k}": v for k, v in rules(diagnostics.secondary_scan(lc, epoch, period, duration)).items()})
    flux = lc.flux.copy()
    flux[np.abs(phase_offset(lc.time, epoch + period / 2, period)) <= duration / 2] -= 0.5 * depth
    injected = LightCurve(lc.time, flux, lc.error, lc.part, lc.parts)
    out.update({f"injected_{k}": v for k, v in
                rules(diagnostics.secondary_scan(injected, epoch, period, duration)).items()})
    return out


def summarize(rows: list[dict]) -> dict:
    def frac(key, subset):
        k = sum(r[key] for r in subset)
        return {"count": k, "total": len(subset), "rate": k / len(subset) if subset else None}
    summary = {}
    groups = ([("all", rows)] + [(s, [r for r in rows if r["split"] == s]) for s in ("calibration", "blind")]
              + [(d, [r for r in rows if r["depth_stratum"] == d]) for d in ("shallow", "medium", "deep")])
    for label, subset in groups:
        summary[label] = {
            "planets_failing_baseline": frac("clean_baseline_fail", subset),
            "planets_failing_calibrated": frac("clean_calibrated_fail", subset),
            "injected_secondary_detected_baseline": frac("injected_baseline_fail", subset),
            "injected_secondary_detected_calibrated": frac("injected_calibrated_fail", subset),
        }
    factors = [r["clean_red_noise_factor"] for r in rows if r["clean_red_noise_factor"]]
    summary["red_noise_factor_quantiles"] = {q: float(np.quantile(factors, q)) for q in (0.1, 0.5, 0.9)}
    return summary


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("manifest", type=Path)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--max-period-over-duration", type=float, default=400.0,
                   help="skip very long scans (period/duration above this) for runtime")
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    root = a.manifest.parent
    with a.manifest.open(newline="") as stream:
        cases = [c for c in csv.DictReader(stream) if c["truth_class"] == "confirmed_planet"]
    cases = [c for c in cases
             if float(c["period_days"]) / float(c["duration_days"]) <= a.max_period_over_duration]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        rows = list(pool.map(evaluate, [(c, str(root)) for c in cases], chunksize=4))
    with a.output.with_suffix(".csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = summarize(rows)
    a.output.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary["all"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
