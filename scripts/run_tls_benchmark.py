#!/usr/bin/env python3
"""Run iterative TLS on a frozen benchmark and calibrate on calibration nulls.

The known-answer truth manifest is read here to map each case to its public
light curve; it is never passed to TLS. Threshold selection uses only
calibration inverted/null controls. TLS's built-in white-noise FAP is not used.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import importlib.metadata
import json
import math
import sys
from pathlib import Path

from run_metrics import run_with_metrics
import evaluate_benchmark as evaluation
import compare_detection_statistics as paired_comparison


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_tls(raw: Path, output: Path) -> None:
    with raw.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    fields = ("rank", "class", "score", "jackknife_score", "period_days",
              "epoch", "duration_days", "depth_ppm", "n_events", "note")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for rank, row in enumerate(rows, 1):
            try:
                sde = float(row["sde"])
                period = float(row["period_days"])
                epoch = float(row["epoch_days"])
                duration = float(row["duration_days"])
                depth = float(row["depth_ppm"])
                events = int(row["transit_count"])
            except (KeyError, TypeError, ValueError):
                continue
            if not all(math.isfinite(value) for value in (sde, period, epoch, duration, depth)):
                continue
            writer.writerow({"rank": rank, "class": "periodic_dip", "score": sde,
                             "jackknife_score": sde, "period_days": period,
                             "epoch": epoch, "duration_days": duration,
                             "depth_ppm": depth, "n_events": events,
                             "note": f"TLS SDE={sde:.6g}; raw_output={raw.name}"})


def run_case(case: dict, data: Path, output: Path, tls_script: Path,
             config: dict, args) -> tuple[dict, dict]:
    """Run or resume one independent light curve; safe for outer concurrency."""
    case_id = case["case_id"]
    case_dir = output / case_id
    raw = case_dir / "tls_raw.csv"
    candidates = case_dir / "candidates.csv"
    sidecar = case_dir / "run.json"
    state = {"input_sha256": sha256(data), "config": config}
    try:
        cached = json.loads(sidecar.read_text()) if sidecar.is_file() else {}
    except (OSError, json.JSONDecodeError):
        cached = {}
    def compatible_state(previous: dict) -> bool:
        if previous.get("input_sha256") != state["input_sha256"]:
            return False
        old_config = dict(previous.get("config", {}))
        new_config = dict(state["config"])
        # Worker and backend thread counts affect resource use only; TLS scans
        # each trial independently, so a resumed run can use available CPUs.
        for key in ("workers", "threads"):
            old_config.pop(key, None)
            new_config.pop(key, None)
        return old_config == new_config

    if (not args.rerun and raw.is_file() and compatible_state(cached.get("state", {}))
            and cached.get("raw_sha256") == sha256(raw)):
        metrics = cached.get("metrics", {})
    else:
        case_dir.mkdir(parents=True, exist_ok=True)
        command = [sys.executable, str(tls_script), "--lightcurve", str(data),
                   "--output", str(raw), "--min-period", str(args.min_period),
                   "--max-period", str(args.max_period), "--max-planets",
                   str(args.max_planets), "--min-sde", "0", "--oversampling",
                   str(args.oversampling), "--threads", str(args.threads)]
        metrics = run_with_metrics(command, case_dir / "tls.log", args.timeout)
        if metrics["timed_out"] or metrics["exit_code"] != 0 or not raw.is_file():
            raise RuntimeError(f"{case_id}: TLS failed/timeout; inspect {case_dir / 'tls.log'}")
        temporary = sidecar.with_suffix(".json.tmp")
        temporary.write_text(json.dumps({"state": state, "command": command,
                                         "raw_sha256": sha256(raw), "metrics": metrics},
                                        indent=2) + "\n")
        temporary.replace(sidecar)
    normalize_tls(raw, candidates)
    row = dict(case)
    row["candidates_csv"] = str(candidates)
    return row, {"case_id": case_id, **metrics}


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth-manifest", required=True, type=Path)
    parser.add_argument("--run-manifest", required=True, type=Path,
                        help="label-free manifest mapping case ids to lightcurve CSVs")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--tls-script", type=Path, default=root / "scripts/tls_search.py")
    parser.add_argument("--min-period", type=float, default=0.15)
    parser.add_argument("--max-period", type=float, default=30.0)
    parser.add_argument("--max-planets", type=int, default=12)
    parser.add_argument("--oversampling", type=int, default=3)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--workers", type=int, default=1,
                        help="number of independent TLS light curves to search concurrently")
    parser.add_argument("--timeout", type=float, default=3600)
    parser.add_argument("--target-trigger-rate", type=float, default=0.05)
    parser.add_argument("--reference-report", type=Path,
                        help="astrosearch evaluation.json; reports a paired blind comparison")
    parser.add_argument("--rerun", action="store_true")
    args = parser.parse_args()
    if (not 0 < args.min_period < args.max_period or args.max_planets < 1
            or args.oversampling < 1 or args.threads < 1 or args.workers < 1 or args.timeout <= 0
            or not 0 < args.target_trigger_rate < 1):
        parser.error("invalid period range, planet count, resources, timeout, or trigger rate")
    truth_path = args.truth_manifest.expanduser().resolve()
    run_path = args.run_manifest.expanduser().resolve()
    output = args.output.expanduser().resolve()
    tls_script = args.tls_script.expanduser().resolve()
    if not truth_path.is_file() or not run_path.is_file() or not tls_script.is_file():
        parser.error("truth manifest, run manifest, or TLS script is missing")

    try:
        tls_version = importlib.metadata.version("transitleastsquares")
    except importlib.metadata.PackageNotFoundError:
        parser.error("install optional TLS dependency with uv sync --extra tls")
    truth = evaluation.read_csv(truth_path)
    run_rows = evaluation.read_csv(run_path)
    data_by_id = {row["case_id"]: (run_path.parent / row["data_csv"]).resolve()
                  for row in run_rows}
    required = [row for row in truth if
                (row["split"] == "calibration" and row["truth_class"] in {"inverted_null", "null"})
                or row["split"] == "blind"]
    output.mkdir(parents=True, exist_ok=True)
    config = {"min_period": args.min_period, "max_period": args.max_period,
              "max_planets": args.max_planets, "oversampling": args.oversampling,
              "threads": args.threads, "workers": args.workers, "min_sde": 0.0,
              "tls_version": tls_version,
              "tls_script_sha256": sha256(tls_script)}
    jobs = []
    for case in required:
        case_id = case["case_id"]
        data = data_by_id.get(case_id)
        if data is None or not data.is_file():
            parser.error(f"{case_id}: required light curve missing from run manifest")
        jobs.append((case, data))
    normalized_truth = []
    measurements = []
    case_order = {case["case_id"]: index for index, case in enumerate(required, 1)}
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        future_cases = {executor.submit(run_case, case, data, output, tls_script,
                                        config, args): case
                        for case, data in jobs}
        finished = 0
        for future in as_completed(future_cases):
            case = future_cases[future]
            try:
                row, metrics = future.result()
            except Exception as error:
                for pending in future_cases:
                    pending.cancel()
                parser.error(str(error))
            normalized_truth.append(row)
            measurements.append(metrics)
            finished += 1
            count = len(evaluation.read_csv(Path(row["candidates_csv"])))
            print(f"[{finished}/{len(required)}] {case['case_id']}: {count} TLS candidates",
                  flush=True)
    normalized_truth.sort(key=lambda row: case_order[row["case_id"]])
    measurements.sort(key=lambda row: case_order[row["case_id"]])

    calibration = [row for row in normalized_truth
                   if row["split"] == "calibration"]
    blind = [row for row in normalized_truth if row["split"] == "blind"]
    threshold, calibration_stats = evaluation.choose_threshold(
        calibration, output, args.target_trigger_rate, 0.0, "score")
    blind_result = evaluation.summarize(blind, threshold, output, 1.0, "score")
    summary = {
        "method": "TLS iterative limb-darkened template",
        "statistic": "tls_sde",
        "tls_version": tls_version,
        "config": config, "threshold_statistic": "TLS SDE",
        "calibrated_threshold": threshold,
        "threshold_from_calibration_controls": threshold,
        "calibration": calibration_stats,
        "evaluation": blind_result,
        "blind": blind_result,
        "runtime": {
            "total_wall_seconds": sum(x.get("wall_seconds", 0) for x in measurements),
            "median_wall_seconds": sorted(x.get("wall_seconds", 0) for x in measurements)[len(measurements)//2],
            "max_sampled_peak_rss_bytes": max((x.get("peak_rss_bytes_sampled") or 0 for x in measurements),
                                               default=0),
            "cases": len(measurements),
        },
        "warning": "TLS SDE is a ranking statistic; empirical threshold uses calibration controls. TLS white-noise FAP is not used.",
    }
    if args.reference_report:
        reference = json.loads(args.reference_report.expanduser().read_text())
        summary["paired_reference_comparison"] = paired_comparison.compare_reports(
            summary, reference)
    (output / "evaluation.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (output / "case_metrics.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(measurements[0]))
        writer.writeheader()
        writer.writerows(measurements)
    print(json.dumps({"threshold": threshold,
                      "blind_planet_recovery": blind_result["confirmed_planet_recovery"],
                      "blind_control_trigger_rate": blind_result["inverted_or_null_control_trigger_rate"]},
                     indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
