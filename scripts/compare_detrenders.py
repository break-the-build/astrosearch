#!/usr/bin/env python3
"""Run paired real-curve injection/null calibration for median and biweight trends."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--min-period", type=float, default=1.0)
    parser.add_argument("--max-period", type=float, default=15.0)
    parser.add_argument("--injections-per-depth", type=int, default=50)
    parser.add_argument("--null-trials", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20261009)
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--period-trials", type=int, default=0)
    parser.add_argument("--depths-ppm", default="100,200,400,800")
    args = parser.parse_args()
    if args.injections_per_depth < 1 or args.null_trials < 1:
        parser.error("use positive injection and null trial counts")
    root = Path(__file__).resolve().parents[1]
    validator = root / "scripts/real_curve_validation.py"
    args.output.mkdir(parents=True, exist_ok=True)
    reports = {}
    for method in ("median", "biweight"):
        destination = args.output / method
        command = [sys.executable, str(validator), "--sector-dir", str(args.sector_dir),
                   "--output", str(destination), "--min-period", str(args.min_period),
                   "--max-period", str(args.max_period), "--injections-per-depth",
                   str(args.injections_per_depth), "--null-trials", str(args.null_trials),
                   "--seed", str(args.seed), "--threads", str(args.threads),
                   "--depths-ppm", args.depths_ppm, "--trend-method", method]
        if args.binary:
            command.extend(("--binary", str(args.binary)))
        if args.period_trials:
            command.extend(("--period-trials", str(args.period_trials)))
        subprocess.run(command, check=True)
        reports[method] = json.loads((destination / "validation_report.json").read_text())
    depths = sorted(set(reports["median"]["injection_recovery"]) |
                    set(reports["biweight"]["injection_recovery"]), key=float)
    comparison = {
        "input_sectors": str(args.sector_dir), "same_seed_paired_trials": args.seed,
        "injections_per_depth": args.injections_per_depth, "null_trials": args.null_trials,
        "trend_methods": {name: reports[name]["trend_method"] for name in reports},
        "injection_recovery": {
            depth: {
                method: reports[method]["injection_recovery"].get(depth)
                for method in ("median", "biweight")
            } for depth in depths
        },
        "maximum_null_scores": {
            method: reports[method]["circular_shift_null"]["maximum_scores"]
            for method in ("median", "biweight")
        },
        "caution": "Paired target-specific diagnostics only. Compare full recovery intervals and null tails across many quiet, active, and gapped stars before changing the default detrending method.",
    }
    output = args.output / "detrender_comparison.json"
    output.write_text(json.dumps(comparison, indent=2) + "\n")
    print(output)


if __name__ == "__main__":
    main()
