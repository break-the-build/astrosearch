#!/usr/bin/env python3
"""Run the frozen, label-free input manifest through astrosearch.

This stage reads only run-manifest.csv, which contains no truth labels or
ephemerides. It is resumable when data, executable, and search configuration
hashes match the per-case sidecar.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import json
import sys
from pathlib import Path
from run_metrics import build_identity, run_with_metrics


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_case(index: int, total: int, row: dict, manifest_dir: Path,
             binary: Path, binary_hash: str, config: dict, build: dict,
             args) -> tuple[int, str, str]:
    case_id = row["case_id"]
    data_path = (manifest_dir / row["data_csv"]).resolve()
    candidate_path = (manifest_dir / row["candidates_csv"]).resolve()
    sidecar = candidate_path.with_suffix(".run.json")
    if not data_path.is_file():
        return index, case_id, f"missing data {data_path}"
    state = {"input_sha256": sha256(data_path), "binary_sha256": binary_hash,
             "config": config, "build_identity": build, "metrics_schema": 1}
    try:
        if (not args.rerun and candidate_path.is_file() and sidecar.is_file()
                and json.loads(sidecar.read_text()).get("state") == state):
            return index, case_id, "cached"
        candidate_path.parent.mkdir(parents=True, exist_ok=True)
        stage = candidate_path.with_name(f".{candidate_path.name}.tmp")
        command = [str(binary), "--input", str(data_path), "--output", str(stage),
                   "--threads", str(args.threads), "--min-period", str(args.min_period),
                   "--max-period", str(args.max_period), "--candidate-grouping", args.grouping,
                   "--transit-search-passes", str(args.transit_search_passes)]
        if args.periods is not None:
            command.extend(("--periods", str(args.periods)))
        log_path = candidate_path.with_suffix(".log")
        metrics = run_with_metrics(command, log_path, timeout=args.search_timeout)
        if metrics["timed_out"] or metrics["exit_code"] != 0 or not stage.is_file():
            stage.unlink(missing_ok=True)
            Path(str(stage) + ".provenance.json").unlink(missing_ok=True)
            return index, case_id, ("search timed out" if metrics["timed_out"] else
                                    f"search failed exit={metrics['exit_code']} (see {log_path})")
        # The search writes provenance beside its (staging) output; keep it
        # next to the published candidates instead of under the temporary name.
        provenance = Path(str(stage) + ".provenance.json")
        if provenance.is_file():
            record = json.loads(provenance.read_text())
            metrics.update(observations=record.get("observations"), span_days=record.get("span_days"),
                           transit_trials=record.get("transit_trials"))
            provenance.replace(Path(str(candidate_path) + ".provenance.json"))
        stage.replace(candidate_path)
        temporary = sidecar.with_name(f".{sidecar.name}.tmp")
        temporary.write_text(json.dumps({"state": state, "command": command,
                                         "candidate_sha256": sha256(candidate_path),
                                         "metrics": metrics}, indent=2) + "\n")
        temporary.replace(sidecar)
        return index, case_id, "searched"
    except Exception as error:
        candidate_path.with_name(f".{candidate_path.name}.tmp").unlink(missing_ok=True)
        return index, case_id, f"{type(error).__name__}: {error}"


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=root / "data/benchmarks/generated/run-manifest.csv")
    parser.add_argument("--binary", type=Path, default=root / "target/release/astrosearch")
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--workers", type=int, default=1,
                        help="parallel target searches; set --threads 1 to avoid oversubscribing CPUs")
    parser.add_argument("--min-period", type=float, default=0.15)
    parser.add_argument("--max-period", type=float, default=30)
    parser.add_argument("--periods", type=int)
    parser.add_argument("--grouping", choices=("ephemeris", "events"), default="events")
    parser.add_argument("--transit-search-passes", type=int, default=5)
    parser.add_argument("--search-timeout", type=float, default=3600,
                        help="maximum seconds for one target search")
    parser.add_argument("--rerun", action="store_true", help="ignore matching resume metadata")
    args = parser.parse_args()
    manifest = args.manifest.expanduser().resolve()
    binary = args.binary.expanduser().resolve()
    if not manifest.is_file() or not binary.is_file():
        parser.error("run manifest or astrosearch executable not found")
    if (not 0 < args.min_period < args.max_period or args.workers < 1
            or not 1 <= args.transit_search_passes <= 12
            or args.search_timeout <= 0):
        parser.error("invalid period bounds, worker count, or transit search pass count")
    build = build_identity(root, binary)
    config = {"threads": args.threads, "workers": args.workers,
              "min_period": args.min_period,
              "max_period": args.max_period, "periods": args.periods,
              "grouping": args.grouping,
              "transit_search_passes": args.transit_search_passes,
              "search_timeout_s": args.search_timeout}
    binary_hash = sha256(binary)
    with manifest.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    failed = 0
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(run_case, i, len(rows), row, manifest.parent,
                                   binary, binary_hash, config, build, args)
                   for i, row in enumerate(rows, 1)]
        for future in as_completed(futures):
            index, case_id, status = future.result()
            stream = sys.stderr if status not in {"searched", "cached"} else sys.stdout
            print(f"[{index}/{len(rows)}] {case_id}: {status}", file=stream, flush=True)
            failed += status not in {"searched", "cached"}
    print(f"completed {len(rows) - failed}/{len(rows)} cases; {failed} failed")
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
