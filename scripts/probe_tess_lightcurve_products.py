#!/usr/bin/env python3
"""Find eligible public TESS reductions for a list of TIC targets."""
from __future__ import annotations

import argparse
import csv
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

PREFERENCE = {"QLP": 0, "TESS-SPOC": 1, "SPOC": 2}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--minimum-exposure", type=float, default=20,
                        help="accept FFI and short-cadence products (seconds)")
    parser.add_argument("--pause", type=float, default=0.2)
    args = parser.parse_args()
    if args.workers < 1 or args.minimum_exposure < 0 or args.pause < 0:
        parser.error("invalid worker, exposure, or pause setting")
    import lightkurve as lk

    with args.targets.expanduser().open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        parser.error("target list is empty")

    def probe(row: dict) -> dict:
        try:
            result = lk.search_lightcurve(f"TIC {row['tic_id']}", mission="TESS")
            choices = []
            columns = set(result.table.colnames)
            sector_key = "sector" if "sector" in columns else "sequence_number" if "sequence_number" in columns else None
            if sector_key:
                for index in range(len(result)):
                    item = result.table[index]
                    try:
                        author = str(item["author"]).strip()
                        exposure = float(item["exptime"])
                        sector = int(item[sector_key])
                    except (KeyError, TypeError, ValueError):
                        continue
                    if author not in PREFERENCE or exposure < args.minimum_exposure:
                        continue
                    filename = str(item["productFilename"]) if "productFilename" in columns else ""
                    choices.append((PREFERENCE[author], exposure, sector, author, filename))
            if not choices:
                return {**row, "preferred_author": "", "available_sectors": 0,
                        "available_exposures_s": "", "probe_status": "no_eligible_product"}
            author = min(choices)[3]
            selected = [item for item in choices if item[3] == author]
            return {**row, "preferred_author": author,
                    "available_sectors": len({item[2] for item in selected}),
                    "available_exposures_s": ";".join(map(str, sorted({int(x[1]) for x in selected}))),
                    "probe_status": "available"}
        except Exception as error:
            return {**row, "preferred_author": "", "available_sectors": 0,
                    "available_exposures_s": "",
                    "probe_status": f"query_error:{type(error).__name__}:{error}"}
        finally:
            if args.pause:
                time.sleep(args.pause)

    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(probe, row): row for row in rows}
        for future in as_completed(futures):
            results.append(future.result())
            if len(results) % 10 == 0 or len(results) == len(rows):
                available = sum(row["probe_status"] == "available" for row in results)
                print(f"probed {len(results)}/{len(rows)}; eligible={available}", flush=True)
    order = {row["tic_id"]: index for index, row in enumerate(rows)}
    results.sort(key=lambda row: order[row["tic_id"]])
    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in results for key in row))
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)
    summary = {"targets": len(results),
               "available": sum(row["probe_status"] == "available" for row in results),
               "by_author": {},
               "errors": sum(str(row["probe_status"]).startswith("query_error") for row in results),
               "output_csv": str(output),
               "author_preference": ["QLP", "TESS-SPOC", "SPOC"]}
    for author in PREFERENCE:
        summary["by_author"][author] = sum(row["preferred_author"] == author for row in results)
    output.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
