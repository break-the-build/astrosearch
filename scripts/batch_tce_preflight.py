#!/usr/bin/env python3
"""Check a survey's detection-tier signals against SPOC TCE sector files once.

Unlike calling candidate_preflight.py separately for every TIC, this batches
all targets, downloads each relevant MAST TCE CSV once, and records exact file
coverage. The result is a catalog screen, not candidate validation.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import json
import re
import tempfile
import time
from pathlib import Path

import requests

from candidate_preflight import Links, ephemeris_matches, select_tce_files


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent,
            prefix=f".{path.name}.", suffix=".tmp", delete=False) as f:
        temp = Path(f.name)
        json.dump(value, f, indent=2, sort_keys=True)
        f.write("\n")
        f.flush()
    temp.replace(path)


def tic_from_path(path: Path) -> str | None:
    for part in path.parts:
        match = re.fullmatch(r"TIC(\d+)", part)
        if match:
            return match.group(1)
    return None


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True,
                        help="completed survey output with TIC*/QLP/candidates.csv")
    parser.add_argument("--coverage", type=Path, required=True,
                        help="CSV containing tic_id and spoc_2min_sectors")
    parser.add_argument("--output", type=Path, required=True, help="JSON TCE check report")
    parser.add_argument("--min-sector", type=int, default=90,
                        help="only query sectors at or above this value (default: 90)")
    parser.add_argument("--timeout", type=float, default=90.0)
    args = parser.parse_args(argv)
    if args.min_sector < 1 or args.timeout <= 0:
        parser.error("min-sector and timeout must be positive")

    requested: dict[str, set[int]] = {}
    with args.coverage.open(newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            tic = str(row.get("tic_id", "")).strip()
            if not tic.isdigit():
                continue
            sectors = set()
            for value in str(row.get("spoc_2min_sectors", "") or "").replace(",", ";").split(";"):
                try:
                    sector = int(value)
                    if sector >= args.min_sector:
                        sectors.add(sector)
                except ValueError:
                    continue
            requested[tic] = sectors

    candidates: dict[str, list[dict]] = {}
    for path in sorted(args.run_dir.glob("TIC*/QLP/candidates.csv")):
        tic = tic_from_path(path)
        if tic is None:
            continue
        with path.open(newline="", encoding="utf-8-sig") as f:
            rows = list(csv.DictReader(f))
        selected = [row for row in rows if row.get("class") in ("periodic_dip", "duo_transit")
                    and row.get("tier") == "detection"]
        if selected:
            candidates[tic] = selected

    target_tics = sorted(set(candidates) & set(requested), key=int)
    sector_union = sorted(set().union(*(requested[tic] for tic in target_tics))) if target_tics else []
    report = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "run_dir": str(args.run_dir.resolve()),
        "coverage_csv": str(args.coverage.resolve()),
        "min_sector": args.min_sector,
        "scope_note": "Checks only SPOC 2-minute sectors at/above min_sector; absence does not rule out a TCE in other sectors.",
        "target_count": len(target_tics),
        "requested_sector_union": sector_union,
        "files": [],
        "targets": {tic: {"requested_sectors": sorted(requested[tic]), "candidate_rows": len(candidates[tic]),
                           "tce_rows": [], "status": "NOT_EVALUATED" if not requested[tic] else "NO_MATCH"}
                    for tic in target_tics},
        "errors": [],
    }
    index_url = "https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_tce.html"
    session = requests.Session()
    for attempt in range(3):
        try:
            response = session.get(index_url, timeout=args.timeout)
            response.raise_for_status()
            links = Links()
            links.feed(response.text)
            file_map = select_tce_files(links.hrefs, sector_union, index_url)
            break
        except Exception as error:
            if attempt == 2:
                report["errors"].append(f"TCE index: {type(error).__name__}: {error}")
                file_map = {}
            else:
                time.sleep(2 ** attempt)

    for url, (start, stop) in sorted(file_map.items()):
        try:
            response = session.get(url, timeout=args.timeout)
            response.raise_for_status()
            content = response.content
            # MAST DVR-TCE tables have comment metadata before the CSV header.
            rows = list(csv.DictReader(line for line in response.text.splitlines()
                                       if not line.lstrip().startswith("#")))
            file_record = {"url": url, "sector_start": start, "sector_stop": stop,
                           "rows": len(rows), "sha256": hashlib.sha256(content).hexdigest(), "status": "ok"}
            report["files"].append(file_record)
        except Exception as error:
            report["files"].append({"url": url, "sector_start": start, "sector_stop": stop,
                                    "status": "error", "error": f"{type(error).__name__}: {error}"})
            continue
        if not rows:
            continue
        keys = {re.sub(r"[^a-z0-9]", "", key.lower()): key for key in rows[0]}
        id_key = next((keys[k] for k in ("ticid", "tic_id", "targetid", "tessinputcatalogid") if k in keys), None)
        period_key = next((keys[k] for k in ("tceperiod", "period", "plorbper") if k in keys), None)
        epoch_key = next((keys[k] for k in ("tcetime0bt", "tceepochbtjd", "tceepoch", "tceepochbjd",
                                            "tceepochbkjd", "tcetime0") if k in keys), None)
        duration_key = next((keys[k] for k in ("tceduration", "duration", "pltrandurh") if k in keys), None)
        if not id_key:
            report["errors"].append(f"unrecognized TCE schema: {url}")
            continue
        for row in rows:
            try:
                tic = str(int(float(row[id_key])))
            except (TypeError, ValueError):
                continue
            if tic not in candidates or not any(start <= s <= stop for s in requested.get(tic, set())):
                continue
            try:
                period = float(row[period_key]) if period_key else None
            except (TypeError, ValueError):
                period = None
            try:
                epoch = float(row[epoch_key]) if epoch_key else None
            except (TypeError, ValueError):
                epoch = None
            if epoch is not None and epoch > 1_000_000:
                epoch -= 2_457_000
            try:
                duration = float(row[duration_key]) / 24.0 if duration_key else None
            except (TypeError, ValueError):
                duration = None
            ephemerides = ephemeris_matches(candidates[tic], period, epoch, duration)
            period_match = any(
                period is not None and period > 0 and
                abs(period - float(candidate["period_days"])) / float(candidate["period_days"]) <= .01
                for candidate in candidates[tic] if candidate.get("period_days")
            )
            report["targets"][tic]["tce_rows"].append({
                "file": url, "period_days": period, "epoch_btjd": epoch,
                "period_match_within_1pct": period_match,
                "ephemeris_matches": ephemerides, "row": row,
            })

    successful_ranges = [(item["sector_start"], item["sector_stop"])
                         for item in report["files"] if item.get("status") == "ok"]
    failed_ranges = [(item["sector_start"], item["sector_stop"])
                     for item in report["files"] if item.get("status") != "ok"]
    for tic, result in report["targets"].items():
        covered = all(any(start <= sector <= stop for start, stop in successful_ranges)
                      for sector in result["requested_sectors"])
        has_failed_coverage = any(any(start <= sector <= stop for start, stop in failed_ranges)
                                  for sector in result["requested_sectors"])
        result["status"] = "MATCH" if any(
            item["period_match_within_1pct"] or item["ephemeris_matches"] for item in result["tce_rows"]
        ) else "NO_DATA" if not result["requested_sectors"] else (
            "PARTIAL" if report["errors"] or has_failed_coverage or not covered else "NO_MATCH")
    # Do not replace a previous successful snapshot with an empty report when
    # the catalog index is temporarily unreachable.
    previous_success = None
    current_incomplete = bool(report["errors"]) or any(
        item.get("status") != "ok" for item in report["files"])
    if current_incomplete and args.output.exists():
        try:
            old = json.loads(args.output.read_text(encoding="utf-8"))
            if old.get("files") and not old.get("errors"):
                previous_success = old
        except (OSError, json.JSONDecodeError):
            pass
    if previous_success is None:
        write_json(args.output, report)
    displayed_report = previous_success or report
    for tic, result in displayed_report["targets"].items():
        preflight = args.run_dir / f"TIC{tic}" / "preflight.md"
        if not preflight.exists():
            continue
        text = preflight.read_text(encoding="utf-8")
        line = (f"- Batch SPOC TCE check: **{result['status']}** in sectors "
                f"{', '.join(map(str, result['requested_sectors'])) or 'none'}; see "
                f"`{args.output.name}` for checked-file hashes and any matching rows. "
                "This check is limited to recent SPOC 2-minute sectors.")
        stale = "- SPOC TCE files were not queried because --tce-sectors was not supplied."
        if stale in text:
            text = text.replace(stale, line)
        elif "Batch SPOC TCE check:" in text:
            text = re.sub(r"^- Batch SPOC TCE check:.*$", line, text, count=1, flags=re.MULTILINE)
        else:
            text += "\n\n## Batch SPOC TCE check\n\n" + line + "\n"
        preflight.write_text(text, encoding="utf-8")
    file_errors = sum(item.get("status") != "ok" for item in report["files"])
    error_count = len(report["errors"]) + file_errors
    print(f"checked {len(displayed_report['files'])} TCE files once for {len(target_tics)} targets; "
          f"matches={sum(v['status']=='MATCH' for v in displayed_report['targets'].values())}; "
          f"errors={error_count}" + ("; retained previous successful snapshot" if previous_success else ""))
    return 1 if error_count else 0


if __name__ == "__main__":
    raise SystemExit(main())
