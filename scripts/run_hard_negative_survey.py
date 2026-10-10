#!/usr/bin/env python3
"""Run the standard TESS survey over disjoint hard-negative target chunks."""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def write_targets(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", type=Path,
                        default=root / "data/benchmarks/hard-negatives/survey-targets.csv")
    parser.add_argument("--output", type=Path,
                        default=root / "data/benchmarks/hard-negatives/survey-parallel")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--threads", type=int, default=6,
                        help="Rust search threads per independent survey process")
    parser.add_argument("--minimum-exposure", type=float, default=20)
    parser.add_argument("--network-timeout", type=float, default=60)
    parser.add_argument("--pause", type=float, default=0.5)
    parser.add_argument("--rerun", action="store_true")
    args = parser.parse_args()
    if (args.workers < 1 or args.threads < 1 or args.network_timeout <= 0
            or args.pause < 0 or args.minimum_exposure < 0):
        parser.error("invalid worker, thread, timeout, or pause setting")
    targets = args.targets.expanduser().resolve()
    output = args.output.expanduser().resolve()
    with targets.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        fields = list(reader.fieldnames or [])
        all_rows = list(reader)
    rows = ([row for row in all_rows if row.get("probe_status") == "available"]
            if "probe_status" in fields else all_rows)
    if not rows:
        parser.error("target list is empty")
    output.mkdir(parents=True, exist_ok=True)
    count = min(args.workers, len(rows))
    by_author: dict[str, list[dict]] = {}
    for row in rows:
        author = row.get("preferred_author", "").strip() or "QLP"
        by_author.setdefault(author, []).append(row)
    chunks = []
    for author in sorted(by_author):
        for index in range(count):
            chunk = by_author[author][index::count]
            if chunk:
                chunks.append((author, chunk))
    script = root / "scripts/tess_survey.py"
    binary = root / "target/release/astrosearch"
    cache = output / "mast-cache"
    cache.mkdir(parents=True, exist_ok=True)

    def run_chunk(index: int, author: str, chunk: list[dict]) -> tuple[int, int, int, str]:
        safe_author = author.lower().replace("/", "-")
        chunk_dir = output / f"{safe_author}-chunk-{index + 1:02d}"
        chunk_dir.mkdir(parents=True, exist_ok=True)
        list_path = chunk_dir / "targets.csv"
        write_targets(list_path, fields, chunk)
        command = [sys.executable, str(script), "--targets", str(list_path),
                   "--output", str(chunk_dir), "--binary", str(binary),
                   "--authors", author, "--min-sector", "1", "--min-exptime", str(args.minimum_exposure),
                   "--max-sector-noise-ratio", "0", "--min-baseline", "0",
                   "--min-teff", "2700", "--max-teff", "10000",
                   "--allow-missing-distance", "--max-radius", "100", "--max-tmag", "15",
                   "--min-period", "0.15", "--max-period", "30", "--threads", str(args.threads),
                   "--trend-method", "biweight", "--trend-days", "1.05",
                   "--segment-gap-days", "0.5", "--transit-search-passes", "1",
                   "--candidate-grouping", "events", "--search-timeout", "3600",
                   "--network-timeout", str(args.network_timeout), "--pause", str(args.pause),
                   "--ledger", str(chunk_dir / "ledger.sqlite"),
                   "--candidate-registry", str(chunk_dir / "candidate_registry.sqlite"),
                   "--cache-dir", str(cache)]
        if args.rerun:
            command.append("--revisit")
        log = chunk_dir / "survey.log"
        with log.open("w", encoding="utf-8") as stream:
            result = subprocess.run(command, cwd=root, stdout=stream,
                                    stderr=subprocess.STDOUT, check=False)
        manifests = list(chunk_dir.glob("TIC*/manifest.json"))
        searched = 0
        for path in manifests:
            try:
                data = json.loads(path.read_text())
            except (OSError, json.JSONDecodeError):
                continue
            searched += any(item.get("author") == "QLP" and item.get("status") == "searched"
                            for item in data.get("authors", []))
        return index, result.returncode, searched, author

    failures = 0
    with ThreadPoolExecutor(max_workers=count) as executor:
        future_chunks = {executor.submit(run_chunk, index, author, chunk): (index, author)
                         for index, (author, chunk) in enumerate(chunks)}
        for future in as_completed(future_chunks):
            index, code, searched, author = future.result()
            print(f"{author} chunk {index + 1}/{len(chunks)}: searched={searched}, exit={code}",
                  flush=True)
            failures += code != 0
    summary = {"targets_in_input": len(all_rows), "targets_with_supported_products": len(rows),
               "workers": args.workers,
               "threads_per_rust_search": args.threads,
               "tasks": len(chunks), "successful_chunks": len(chunks) - failures,
               "output": str(output),
               "note": "Each worker has an independent survey ledger and output tree; MAST cache is shared."}
    (output / "batch-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
