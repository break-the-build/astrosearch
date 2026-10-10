#!/usr/bin/env python3
"""Run the TIC 4206066 real-signal benchmark as a search-quality release gate.

The published ~3.18-day signal is a candidate, not a confirmed planet. This
check only tests whether the configured search returns a matching periodic dip;
it does not validate the astrophysical interpretation.
"""
import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path
csv.field_size_limit(2**31 - 1)  # candidate rows carry long evidence/alias cells


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "data/tic4206066/combined.csv")
    parser.add_argument("--binary", type=Path, default=root / "target/release/astrosearch")
    parser.add_argument("--output", type=Path, default=root / "data/tic4206066/benchmark_candidates.csv")
    parser.add_argument("--expected-period", type=float, default=3.18)
    parser.add_argument("--reference-epoch", type=float, required=True,
                        help="published reference transit midpoint in the input time system (BTJD)")
    parser.add_argument("--max-timing-error-hours", type=float, default=1.0,
                        help="maximum phase offset and accumulated baseline period drift")
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--trend-method", choices=("median", "biweight"), default="biweight")
    args = parser.parse_args()
    binary = args.binary.expanduser()
    if not binary.is_file():
        found = shutil.which(str(args.binary))
        binary = Path(found) if found else Path(".")
    if not args.input.is_file() or not binary.is_file():
        parser.error("input light curve or astrosearch executable not found")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    command = [str(binary), "--input", str(args.input), "--output", str(args.output),
               "--threads", str(args.threads), "--trend-method", args.trend_method]
    subprocess.run(command, check=True)
    data = []
    with args.input.open(newline="") as stream:
        for row in csv.DictReader(stream):
            try:
                data.append(float(row.get("time", row.get("btjd", row.get("bjd", "nan")))))
            except (TypeError, ValueError):
                pass
    if not data or args.expected_period <= 0 or args.max_timing_error_hours <= 0:
        parser.error("expected period and timing tolerance must be positive and input times must be readable")
    span = max(data) - min(data)
    tolerance_days = args.max_timing_error_hours / 24
    with args.output.open(newline="") as stream:
        candidates = [row for row in csv.DictReader(stream)
                      if row.get("class") in ("periodic_dip", "duo_transit")]
    matches = []
    for rank, row in enumerate(candidates, 1):
        try:
            period = float(row["period_days"])
            events = int(row["n_events"])
        except (KeyError, ValueError):
            continue
        period_drift = abs(period - args.expected_period) * span / args.expected_period
        phase_offset = abs((float(row["epoch"]) - args.reference_epoch) % args.expected_period)
        phase_offset = min(phase_offset, args.expected_period - phase_offset)
        if period_drift <= tolerance_days and phase_offset <= tolerance_days and events >= 2:
            matches.append({"rank_within_output": rank, "period_days": period,
                            "period_drift_days_over_baseline": period_drift,
                            "phase_offset_days": phase_offset,
                            "n_events": events, "score": float(row["score"]),
                            "epoch_btjd": float(row["epoch"]),
                            "duration_days": float(row["duration_days"])})
    result = {
        "benchmark": "TIC 4206066 published ~3.18-day transit-like candidate",
        "input": str(args.input), "command": command,
        "expected_period_days": args.expected_period,
        "reference_epoch_btjd": args.reference_epoch,
        "baseline_days": span,
        "max_timing_error_hours": args.max_timing_error_hours,
        "passed": bool(matches), "matching_candidates": matches,
        "interpretation": "Search-recovery gate only; the signal is not confirmed as a planet.",
    }
    report_path = args.output.with_suffix(".benchmark.json")
    report_path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    print(f"benchmark report: {report_path}", file=sys.stderr)
    return 0 if matches else 1


if __name__ == "__main__":
    raise SystemExit(main())
