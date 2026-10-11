#!/usr/bin/env python3
"""Injection-recovery into real TESS light curves (Experiment B).

Hosts are confirmed-planet light curves from the known-answer benchmark with
the known planet's transits removed (so real systematics and stellar noise
remain). Into each, a limb-darkened transit (batman, quadratic u = 0.40,
0.25, b = 0) is injected at a period drawn log-uniformly from 0.7-12 d, a
duration from a solar-density scaling, and a depth set by a target
white-noise SNR. astrosearch searches the result with five-pass Profile J by
default (0.5-15 d); an
injection is recovered when a transit candidate matches its ephemeris
(evaluate_benchmark.close_ephemeris). Recovery is reported at the screening
floor (6) and at the calibrated detection tier, by SNR bin.

The detection tier comes from a calibration JSON and is not adjusted here.
``--threshold-source evaluation`` selects the calibration-split tier for an
unbiased paired profile comparison. Hosts are drawn with a fixed seed.
"""
from __future__ import annotations

import argparse
import csv
import json
import shlex
import math
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

try:
    import resource
except ImportError:  # pragma: no cover - unavailable on Windows
    resource = None

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
    # batman 2.5.x imports distutils directly; setuptools provides its
    # compatibility shim on Python 3.12 where distutils was removed.
    import setuptools  # noqa: F401
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
    case, root, work, binary, snr, seed, threshold, extra, timeout, transit_passes, profile = job
    rng = np.random.default_rng(seed)
    built = build_injection(case, Path(root), rng, snr)
    if built is None:
        return {"case_id": case["case_id"], "seed": seed, "snr_target": snr,
                "skipped": "too few in-transit points"}
    rows, truth = built
    curve = Path(work) / f"{case['case_id']}_snr{int(snr)}_{seed}.csv"
    out = curve.with_suffix(".candidates.csv")
    np.savetxt(curve, rows, delimiter=",", header="time,flux,flux_err", comments="", fmt="%.8f")
    child_before = (resource.getrusage(resource.RUSAGE_CHILDREN)
                    if resource is not None else None)
    search_started = time.monotonic()
    error = None
    timed_out = False
    try:
        subprocess.run([binary, "--search-profile", profile,
                        "--input", str(curve), "--output", str(out), "--threads", "1",
                        "--transit-search-passes", str(transit_passes),
                        "--min-period", "0.5", "--max-period", "15",
                        *shlex.split(extra)], check=True, capture_output=True, text=True,
                       timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        error = f"search timed out after {timeout:g}s"
    except subprocess.CalledProcessError as failure:
        detail = (failure.stderr or failure.stdout or "").strip().replace("\n", " ")
        error = f"search failed exit={failure.returncode}" + (f": {detail[:500]}" if detail else "")
    search_wall = time.monotonic() - search_started
    search_cpu = None
    if resource is not None and child_before is not None:
        child_after = resource.getrusage(resource.RUSAGE_CHILDREN)
        search_cpu = ((child_after.ru_utime + child_after.ru_stime)
                      - (child_before.ru_utime + child_before.ru_stime))
    if error:
        return {"case_id": case["case_id"], "seed": seed, "split": case["split"], "snr_target": snr,
                "search_wall_seconds": search_wall, "search_cpu_seconds": search_cpu,
                "timed_out": timed_out, "error": error}
    with out.open(newline="") as stream:
        transit = [r for r in csv.DictReader(stream) if r["class"] in ("periodic_dip", "duo_transit")]
    span = float(rows[-1, 0] - rows[0, 0])
    matched_row = next((r for r in transit
                        if any(ev.close_ephemeris(e, truth, span, 1.0) for e in ev.ephemerides(r))), None)
    matched = float(matched_row["score"]) if matched_row else None
    jackknife = jackknife_score(matched_row) if matched_row else None
    curve.unlink()
    return {"case_id": case["case_id"], "seed": seed, "split": case["split"], **truth,
            "search_wall_seconds": search_wall, "search_cpu_seconds": search_cpu,
            "matched_score": matched,
            "recovered_screening": matched is not None and matched >= 6.0,
            "matched_jackknife": jackknife,
            "recovered_detection": jackknife is not None and jackknife >= threshold,
            "timed_out": False}


def summarize(rows: list[dict], threshold: float) -> dict:
    done = [r for r in rows if "skipped" not in r and "error" not in r]
    out = {"threshold": threshold, "injections": len(done),
           "skipped": sum("skipped" in r for r in rows),
           "failed": sum("error" in r for r in rows),
           "timed_out": sum(bool(r.get("timed_out")) for r in rows), "by_snr": {}}
    searched = [r for r in rows if "skipped" not in r]
    measured_cpu = [r["search_cpu_seconds"] for r in searched
                    if r.get("search_cpu_seconds") is not None]
    out["search_metrics"] = {
        "searches_attempted": len(searched),
        "cpu_seconds_sum": sum(measured_cpu) if measured_cpu else None,
        "cpu_measured": len(measured_cpu),
        "wall_seconds_sum": sum(r.get("search_wall_seconds", 0.0) for r in searched),
    }
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
    p.add_argument("--search-timeout", type=float, default=3600,
                   help="maximum seconds for one injected-curve search")
    p.add_argument("--search-profile", choices=("fast", "standard"), default="fast",
                   help="search configuration (default: Profile J / fast)")
    p.add_argument("--transit-search-passes", type=int,
                   help="number of iterative transit-search passes (defaults: 5 for fast, 1 for standard)")
    p.add_argument("--extra-search-args", default="", help="additional astrosearch options")
    p.add_argument("--calibration", type=Path,
                   help="detection-tier calibration; defaults to one matched to --search-profile")
    p.add_argument("--threshold-source", choices=("deployment", "evaluation"),
                   default="deployment",
                   help="use the all-controls deployment threshold or the calibration-split threshold")
    p.add_argument("--work-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.transit_search_passes is None:
        a.transit_search_passes = 5 if a.search_profile == "fast" else 1
    if a.calibration is None:
        calibration_name = ("profile-j-fivepass-20261010-fa05-verified.json"
                            if a.search_profile == "fast" else "detection_threshold.json")
        a.calibration = Path(__file__).resolve().parents[1] / "data/calibration" / calibration_name
    if a.transit_search_passes < 1:
        p.error("--transit-search-passes must be at least 1")
    root = a.manifest.parent
    calibration = json.loads(a.calibration.read_text())
    threshold = float(calibration["threshold"] if a.threshold_source == "deployment"
                      else calibration["evaluation"]["threshold_from_calibration_split"])
    hosts = [c for c in ev.read_csv(a.manifest) if c["truth_class"] == "confirmed_planet"
             and float(c["observations"]) <= a.max_observations]
    rng = np.random.default_rng(a.seed)
    chosen = [hosts[i] for i in rng.choice(len(hosts), size=min(a.hosts, len(hosts)), replace=False)]
    a.work_dir.mkdir(parents=True, exist_ok=True)
    batch_started = time.monotonic()
    jobs = [(case, str(root), str(a.work_dir), a.binary, snr, a.seed + 1000 * i + j, threshold,
            a.extra_search_args, a.search_timeout, a.transit_search_passes, a.search_profile)
            for i, case in enumerate(chosen) for j, snr in enumerate(SNR_LEVELS)]
    rows = []
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futures = [pool.submit(run_one, job) for job in jobs]
        for index, future in enumerate(as_completed(futures), 1):
            row = future.result()
            rows.append(row)
            status = "failed" if "error" in row else "skipped" if "skipped" in row else "searched"
            print(f"[{index}/{len(jobs)}] {row['case_id']} S/N {row['snr_target']:g}: {status}",
                  flush=True)
    with a.output.with_suffix(".csv").open("w", newline="") as stream:
        fields = sorted({k for r in rows for k in r})
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    summary = summarize(rows, threshold)
    summary["threshold_source"] = a.threshold_source
    summary["transit_search_passes"] = a.transit_search_passes
    summary["search_profile"] = a.search_profile
    summary["calibration_file"] = str(a.calibration.resolve())
    summary["batch_wall_seconds"] = time.monotonic() - batch_started
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary["by_snr"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
