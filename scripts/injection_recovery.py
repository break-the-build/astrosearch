#!/usr/bin/env python3
"""Injection-recovery into real TESS light curves (Experiment B).

Hosts are confirmed-planet light curves from the known-answer benchmark with
the known planet's transits removed (so real systematics and stellar noise
remain). Into each, a limb-darkened transit (batman, quadratic u = 0.40,
0.25, b = 0) is injected at a period drawn log-uniformly from 0.7-12 d, a
duration from a solar-density scaling, and a depth set by a target
white-noise SNR. astrosearch searches the result (one pass, 0.5-15 d); an
injection is recovered when a transit candidate matches its ephemeris
(evaluate_benchmark.close_ephemeris). Recovery is reported at the screening
floor (6) and at the calibrated detection tier, by SNR bin.

The detection tier comes from data/calibration/detection_threshold.json and
is not adjusted here; hosts are drawn with a fixed seed.
"""
from __future__ import annotations

import argparse
import csv
import json
import shlex
import math
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import evaluate_benchmark as ev
from vetting.io import phase_offset, read_table
from vetting.tier import jackknife_score

csv.field_size_limit(2**31 - 1)
SNR_LEVELS = (7.0, 10.0, 15.0, 20.0, 30.0, 40.0, 60.0)
LIMB_DARKENING = (0.40, 0.25)


def solar_duration(period: float) -> float:
    """Central-transit duration (days) around a Sun-like star."""
    return float(np.clip(13 / 24 * (period / 365.25) ** (1 / 3), 0.04, 0.35))


def transit_model(time: np.ndarray, period: float, epoch: float, duration: float,
                  depth: float) -> np.ndarray:
    import batman

    params = batman.TransitParams()
    params.t0, params.per, params.inc, params.ecc, params.w = epoch, period, 90.0, 0.0, 90.0
    params.rp = math.sqrt(depth)
    # T14 ~ P / (pi a/R*) (1 + k) for b = 0.
    params.a = max(period * (1 + params.rp) / (math.pi * duration), 1.5)
    params.limb_dark, params.u = "quadratic", list(LIMB_DARKENING)
    return batman.TransitModel(params, time).light_curve(params)


def point_scatter(flux: np.ndarray) -> float:
    return float(1.4826 * np.median(np.abs(np.diff(flux))) / math.sqrt(2))


def build_injection(case: dict, root: Path, rng: np.random.Generator, snr: float):
    data = read_table(root / "lightcurves" / f"{case['case_id']}.csv")
    data = data[np.argsort(data[:, 0], kind="stable")]
    known = np.abs(phase_offset(data[:, 0], float(case["epoch_days"]), float(case["period_days"])))
    data = data[known > float(case["duration_days"])]
    time, flux, error = data[:, 0], data[:, 1], data[:, 2]
    period = float(np.exp(rng.uniform(np.log(0.7), np.log(12.0))))
    epoch = float(time[0] + rng.uniform(0, period))
    duration = solar_duration(period)
    inside = np.abs(phase_offset(time, epoch, period)) <= duration / 2
    if inside.sum() < 10:
        return None
    # Limb darkening lowers the mean in-transit depth to ~0.85 of the peak.
    depth = snr * point_scatter(flux) / math.sqrt(inside.sum()) / 0.85
    model = transit_model(time, period, epoch, duration, depth)
    truth = {"period_days": period, "epoch_days": epoch, "duration_days": duration,
             "depth_ppm": depth * 1e6, "snr_target": snr, "in_transit_points": int(inside.sum())}
    return np.column_stack([time, flux * model, error]), truth


def run_one(job: tuple) -> dict:
    case, root, work, binary, snr, seed, threshold, extra = job
    rng = np.random.default_rng(seed)
    built = build_injection(case, Path(root), rng, snr)
    if built is None:
        return {"case_id": case["case_id"], "snr_target": snr, "skipped": "too few in-transit points"}
    rows, truth = built
    curve = Path(work) / f"{case['case_id']}_snr{int(snr)}_{seed}.csv"
    out = curve.with_suffix(".candidates.csv")
    np.savetxt(curve, rows, delimiter=",", header="time,flux,flux_err", comments="", fmt="%.8f")
    subprocess.run([binary, "--input", str(curve), "--output", str(out), "--threads", "1",
                    "--transit-search-passes", "1", "--min-period", "0.5", "--max-period", "15",
                    *shlex.split(extra)],
                   check=True, capture_output=True, text=True)
    with out.open(newline="") as stream:
        transit = [r for r in csv.DictReader(stream) if r["class"] in ("periodic_dip", "duo_transit")]
    span = float(rows[-1, 0] - rows[0, 0])
    matched_row = next((r for r in transit
                        if any(ev.close_ephemeris(e, truth, span, 1.0) for e in ev.ephemerides(r))), None)
    matched = float(matched_row["score"]) if matched_row else None
    jackknife = jackknife_score(matched_row) if matched_row else None
    curve.unlink()
    return {"case_id": case["case_id"], "split": case["split"], **truth, "matched_score": matched,
            "recovered_screening": matched is not None and matched >= 6.0,
            "matched_jackknife": jackknife,
            "recovered_detection": jackknife is not None and jackknife >= threshold}


def summarize(rows: list[dict], threshold: float) -> dict:
    done = [r for r in rows if "skipped" not in r]
    out = {"threshold": threshold, "injections": len(done), "skipped": len(rows) - len(done), "by_snr": {}}
    for snr in SNR_LEVELS:
        sub = [r for r in done if r["snr_target"] == snr]
        k1 = sum(r["recovered_screening"] for r in sub)
        k2 = sum(r["recovered_detection"] for r in sub)
        out["by_snr"][str(snr)] = {"n": len(sub),
                                   "screening": [k1, ev.wilson(k1, len(sub))],
                                   "detection": [k2, ev.wilson(k2, len(sub))]}
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("manifest", type=Path)
    p.add_argument("--binary", required=True)
    p.add_argument("--hosts", type=int, default=30)
    p.add_argument("--max-observations", type=int, default=20000)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--seed", type=int, default=20261010)
    p.add_argument("--extra-search-args", default="", help="additional astrosearch options")
    p.add_argument("--calibration", type=Path,
                   default=Path(__file__).resolve().parents[1] / "data/calibration/detection_threshold.json")
    p.add_argument("--work-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    root = a.manifest.parent
    threshold = float(json.loads(a.calibration.read_text())["threshold"])
    hosts = [c for c in ev.read_csv(a.manifest) if c["truth_class"] == "confirmed_planet"
             and float(c["observations"]) <= a.max_observations]
    rng = np.random.default_rng(a.seed)
    chosen = [hosts[i] for i in rng.choice(len(hosts), size=min(a.hosts, len(hosts)), replace=False)]
    a.work_dir.mkdir(parents=True, exist_ok=True)
    jobs = [(case, str(root), str(a.work_dir), a.binary, snr, a.seed + 1000 * i + j, threshold,
             a.extra_search_args)
            for i, case in enumerate(chosen) for j, snr in enumerate(SNR_LEVELS)]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        rows = list(pool.map(run_one, jobs))
    with a.output.with_suffix(".csv").open("w", newline="") as stream:
        fields = sorted({k for r in rows for k in r})
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    summary = summarize(rows, threshold)
    a.output.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary["by_snr"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
