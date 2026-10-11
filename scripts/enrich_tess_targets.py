#!/usr/bin/env python3
"""Attach real MAST TESS coverage metadata to a TIC target list.

This should run before target-yield ranking when actual sectors are known. It
queries Lightkurve's same MAST search endpoint as the survey, caches each TIC
for resumability, and can emit a second CSV limited to targets meeting the
survey's recent-sector and multiyear-coverage requirements.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import math
import os
import random
import tempfile
import time
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter


class RequestTimeoutAdapter(HTTPAdapter):
    def __init__(self, timeout_seconds: float):
        self.request_timeout = timeout_seconds
        super().__init__(max_retries=0)

    def send(self, request, **kwargs):
        if kwargs.get("timeout") is None:
            kwargs["timeout"] = self.request_timeout
        return super().send(request, **kwargs)


def parse_tic(raw: str) -> str:
    value = str(raw).strip().upper().removeprefix("TIC").strip()
    if value.endswith(".0") and value[:-2].isdigit():
        value = value[:-2]
    if not value.isascii() or not value.isdigit() or int(value) <= 0:
        raise ValueError(f"invalid TIC ID: {raw!r}")
    return value


def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = None
    try:
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent,
                prefix=f".{path.name}.", suffix=".tmp", delete=False) as f:
            temp = Path(f.name)
            json.dump(data, f, indent=2, sort_keys=True)
            f.write("\n")
            f.flush()
            os.fsync(f.fileno())
        os.replace(temp, path)
    finally:
        if temp is not None and temp.exists():
            temp.unlink()


def query_target(tic: str, authors: list[str], min_exptime: float) -> dict:
    import lightkurve as lk

    products = []
    errors = []
    for author in authors:
        result = None
        for attempt in range(3):
            try:
                result = lk.search_lightcurve(f"TIC {tic}", mission="TESS", author=author)
                break
            except (requests.exceptions.RequestException, TimeoutError) as error:
                if attempt == 2:
                    errors.append(f"{author}: {type(error).__name__}: {error}")
                else:
                    time.sleep((2 ** attempt) * random.uniform(0.8, 1.2))
            except Exception as error:
                errors.append(f"{author}: {type(error).__name__}: {error}")
                break
        if result is None:
            continue
        try:
            names = set(result.table.colnames)
            sector_col = "sector" if "sector" in names else "sequence_number" if "sequence_number" in names else None
            if sector_col is None or "exptime" not in names:
                errors.append(f"{author}: no sector/exptime metadata")
                continue
            for row in result.table:
                try:
                    sector = int(row[sector_col])
                    exptime = float(row["exptime"])
                    t0 = float(row["t_min"])
                    t1 = float(row["t_max"])
                except (TypeError, ValueError, KeyError):
                    continue
                if not math.isfinite(exptime) or not math.isfinite(t0) or not math.isfinite(t1) or t1 <= t0:
                    continue
                products.append({"author": author, "sector": sector, "exptime_s": exptime,
                                 "t_min": t0, "t_max": t1})
        except Exception as error:  # retain partial author results and exact error text
            errors.append(f"{author}: {type(error).__name__}: {error}")

    by_author = {}
    for product in products:
        if product["exptime_s"] < min_exptime:
            continue
        by_author.setdefault(product["author"], {})
        sector = product["sector"]
        previous = by_author[product["author"]].get(sector)
        if previous is None or product["exptime_s"] < previous["exptime_s"]:
            by_author[product["author"]][sector] = product

    sector_spans = {}
    for product in products:
        sector = product["sector"]
        span = min(27.4, product["t_max"] - product["t_min"])
        sector_spans[sector] = max(sector_spans.get(sector, 0.0), span)
    author_coverage = {}
    for author, rows in by_author.items():
        selected = list(rows.values())
        author_coverage[author] = {
            "sectors": sorted(rows),
            "baseline_days": max(x["t_max"] for x in selected) - min(x["t_min"] for x in selected)
                              if selected else 0.0,
        }
    return {
        "status": "ok" if products else "no_products",
        "sectors": sorted(sector_spans),
        "observed_days": round(sum(sector_spans.values()), 4),
        "spoc_2min": any(x["author"] == "SPOC" and x["exptime_s"] <= 120.0 for x in products),
        "spoc_2min_sectors": sorted({x["sector"] for x in products
                                      if x["author"] == "SPOC" and x["exptime_s"] <= 120.0}),
        "authors": author_coverage,
        "errors": errors,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("targets", type=Path, help="TIC target CSV, preferably yield-ranked")
    parser.add_argument("--output", type=Path, required=True, help="all input rows with MAST coverage columns")
    parser.add_argument("--eligible-output", type=Path, required=True,
                        help="rows with a usable author baseline and a recent sector")
    parser.add_argument("--cache", type=Path, required=True, help="JSON cache directory; safe to reuse")
    parser.add_argument("--limit", type=int, help="query only the first N rows")
    parser.add_argument("--workers", type=int, default=2, help="concurrent MAST searches (default: 2)")
    parser.add_argument("--timeout", type=float, default=45.0, help="HTTP connect/read timeout in seconds")
    parser.add_argument("--min-exptime", type=float, default=200.0,
                        help="minimum exposure used for survey-eligible sector coverage")
    parser.add_argument("--min-sector", type=int, default=90)
    parser.add_argument("--min-baseline", type=float, default=1095.0)
    parser.add_argument("--authors", nargs="+", default=["QLP", "SPOC"])
    args = parser.parse_args(argv)
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    if args.workers < 1 or args.timeout <= 0 or args.min_exptime < 0 or args.min_baseline <= 0:
        parser.error("workers, timeout, exposure and baseline values must be positive")

    from astroquery.mast import Observations
    Observations._session.mount("https://", RequestTimeoutAdapter(args.timeout))
    with args.targets.open(newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        parser.error("target file is empty")
    limit = len(rows) if args.limit is None else min(len(rows), args.limit)
    selected = rows[:limit]
    cache = {}
    pending = []
    for row in selected:
        tic = parse_tic(row.get("tic_id", row.get("tic", row.get("ID", ""))))
        path = args.cache / f"tic{tic}.json"
        if path.exists():
            try:
                cached = json.loads(path.read_text(encoding="utf-8"))
                if "spoc_2min_sectors" in cached:
                    cache[tic] = cached
                    continue
            except (OSError, json.JSONDecodeError):
                pass
        pending.append(tic)

    def run(tic):
        # Bounded retry handles short network disturbances without long silent hangs.
        last = None
        for attempt in range(3):
            try:
                return tic, query_target(tic, args.authors, args.min_exptime)
            except (requests.exceptions.RequestException, TimeoutError) as error:
                last = error
                if attempt < 2:
                    time.sleep((2 ** attempt) * random.uniform(0.8, 1.2))
        return tic, {"status": "error", "sectors": [], "observed_days": 0.0,
                     "spoc_2min": False, "authors": {}, "errors": [f"{type(last).__name__}: {last}"]}

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run, tic): tic for tic in pending}
        for future in concurrent.futures.as_completed(futures):
            tic, result = future.result()
            cache[tic] = result
            write_json(args.cache / f"tic{tic}.json", result)
            print(f"cached TIC {tic}: {result['status']} ({len(result.get('sectors', []))} sectors)", flush=True)

    columns = ["sectors", "observed_days", "spoc_2min", "spoc_2min_sectors", "sector_list",
               "qlp_sectors", "qlp_baseline_days", "spoc_sectors", "spoc_baseline_days",
               "coverage_eligible", "coverage_errors"]
    enriched = []
    eligible = []
    for index, original in enumerate(rows):
        row = dict(original)
        if index < limit:
            tic = parse_tic(row.get("tic_id", row.get("tic", row.get("ID", ""))))
            result = cache[tic]
            row["sectors"] = len(result.get("sectors", []))
            row["observed_days"] = result.get("observed_days", 0.0)
            row["spoc_2min"] = str(bool(result.get("spoc_2min"))).lower()
            row["spoc_2min_sectors"] = ";".join(map(str, result.get("spoc_2min_sectors", [])))
            row["sector_list"] = ";".join(map(str, result.get("sectors", [])))
            for author, prefix in (("QLP", "qlp"), ("SPOC", "spoc")):
                coverage = result.get("authors", {}).get(author, {})
                row[f"{prefix}_sectors"] = ";".join(map(str, coverage.get("sectors", [])))
                row[f"{prefix}_baseline_days"] = round(float(coverage.get("baseline_days", 0.0)), 4)
            row["coverage_errors"] = " | ".join(result.get("errors", []))
            row["coverage_eligible"] = str(any(
                float(result.get("authors", {}).get(a, {}).get("baseline_days", 0.0)) >= args.min_baseline
                and any(s >= args.min_sector for s in result.get("authors", {}).get(a, {}).get("sectors", []))
                for a in args.authors)).lower()
        else:
            for column in columns:
                row.setdefault(column, "")
        enriched.append(row)
        if row.get("coverage_eligible") == "true":
            eligible.append(row)

    fields = list(dict.fromkeys([*rows[0].keys(), *columns]))
    for path, output_rows in ((args.output, enriched), (args.eligible_output, eligible)):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(output_rows)
    print(f"coverage-enriched {limit} of {len(rows)} rows; {len(eligible)} meet >=S{args.min_sector} "
          f"and >= {args.min_baseline:g} d in one requested author")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
