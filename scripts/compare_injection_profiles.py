#!/usr/bin/env python3
"""Report paired S/N recovery and resource differences between injections."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_csv(path: Path) -> dict[tuple[str, str, str], dict]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return {(r["case_id"], r["snr_target"], r["seed"]): r for r in csv.DictReader(stream)}


def wilson(k: int, n: int) -> list[float] | None:
    if not n:
        return None
    z, p = 1.959963984540054, k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / d
    return [c - h, c + h]


def flag(row: dict, name: str) -> bool:
    return row.get(name, "").strip().lower() in {"true", "1"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-a-csv", type=Path, required=True)
    parser.add_argument("--profile-a-summary", type=Path, required=True)
    parser.add_argument("--profile-a-name", default="standard")
    parser.add_argument("--profile-b-csv", type=Path, required=True)
    parser.add_argument("--profile-b-summary", type=Path, required=True)
    parser.add_argument("--profile-b-name", default="profile-j")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    a, b = read_csv(args.profile_a_csv), read_csv(args.profile_b_csv)
    common = sorted(a.keys() & b.keys())
    by_snr = {}
    for snr in ("7.0", "10.0", "15.0", "20.0", "30.0", "40.0", "60.0"):
        ids = [key for key in common if key[1] == snr and "error" not in a[key]
               and "error" not in b[key] and "skipped" not in a[key] and "skipped" not in b[key]]
        av = [flag(a[key], "recovered_detection") for key in ids]
        bv = [flag(b[key], "recovered_detection") for key in ids]
        a_recovered, b_recovered = sum(av), sum(bv)
        a_only = sum(x and not y for x, y in zip(av, bv))
        b_only = sum(y and not x for x, y in zip(av, bv))
        by_snr[snr] = {args.profile_a_name: {"recovered": a_recovered, "n": len(ids),
                                             "rate": a_recovered / len(ids) if ids else None,
                                             "wilson_95": wilson(a_recovered, len(ids))},
                       args.profile_b_name: {"recovered": b_recovered, "n": len(ids),
                                             "rate": b_recovered / len(ids) if ids else None,
                                             "wilson_95": wilson(b_recovered, len(ids))},
                       "paired": {"a_only": a_only, "b_only": b_only,
                                  "both": sum(x and y for x, y in zip(av, bv)),
                                  "neither": sum(not x and not y for x, y in zip(av, bv)),
                                  "b_minus_a": (b_only - a_only) / len(ids) if ids else None}}
    a_summary = json.loads(args.profile_a_summary.read_text())
    b_summary = json.loads(args.profile_b_summary.read_text())
    report = {"paired_injections": len(common),
              "valid_pairs": sum(row[args.profile_a_name]["n"] for row in by_snr.values()),
              "paired_key_coverage": len(common) / len(a) if a else None,
              "resource_metrics": {args.profile_a_name: a_summary.get("search_metrics"),
                                   args.profile_b_name: b_summary.get("search_metrics")},
              "batch_wall_seconds": {args.profile_a_name: a_summary.get("batch_wall_seconds"),
                                     args.profile_b_name: b_summary.get("batch_wall_seconds")},
              "failed_or_skipped": {args.profile_a_name: {"failed": a_summary.get("failed"),
                                                           "skipped": a_summary.get("skipped"),
                                                           "timed_out": a_summary.get("timed_out")},
                                    args.profile_b_name: {"failed": b_summary.get("failed"),
                                                           "skipped": b_summary.get("skipped"),
                                                           "timed_out": b_summary.get("timed_out")}},
              "by_snr": by_snr}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
