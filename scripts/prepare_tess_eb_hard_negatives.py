#!/usr/bin/env python3
"""Select a deterministic, catalog-labeled TESS EB cohort for hard-negative tests.

The source is the published table of 7,936 new, vetted TESS eclipsing binaries
from Kostov et al. (2025), mirrored by VizieR as J/ApJS/279/50/table3. This
script selects target IDs and metadata only; light curves and Rust outputs are
obtained separately with the regular TESS survey pipeline.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

SOURCE_REFERENCE = (
    "Kostov et al. 2025, ApJS 279, 50, doi:10.3847/1538-4365/ade2d8; "
    "VizieR J/ApJS/279/50/table3, doi:10.26093/cds/vizier/bb-8q-1j"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_catalog(path: Path) -> list[dict]:
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    header = next((i for i, line in enumerate(lines)
                   if line.startswith("TIC\tRAJ2000\t")), None)
    if header is None:
        raise ValueError("VizieR TSV header not found")
    columns = lines[header].split("\t")
    result = []
    for line in lines[header + 1:]:
        if not line or line.startswith("-") or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != len(columns):
            continue
        row = dict(zip(columns, fields))
        try:
            tic = str(int(row["TIC"].strip()))
            period = float(row["Per"])
            tmag = float(row["Tmag"])
            teff = float(row["Teff"])
            epoch = float(row["T0-pri"])
            duration_hours = float(row["Dur-pri"])
            sectors = int(row["Sect"])
        except (KeyError, TypeError, ValueError, OverflowError):
            continue
        if (not all(math.isfinite(v) for v in (period, tmag, teff, epoch, duration_hours))
                or not 0.15 <= period <= 30 or duration_hours <= 0 or sectors < 1):
            continue
        result.append({
            "tic_id": tic, "host_id": tic, "period_days": period,
            "epoch_days": epoch, "duration_days": duration_hours / 24,
            "teff_k": teff, "tmag": tmag, "catalog_sectors": sectors,
            "ra_deg": row["RAJ2000"].strip(), "dec_deg": row["DEJ2000"].strip(),
            "depth_primary_ppt": row["Depth-pri"].strip(),
            "depth_secondary_ppt": row["Depth-sec"].strip(),
            "label_reference": SOURCE_REFERENCE,
        })
    return result


def period_band(period: float) -> str:
    return "short" if period < 3 else "medium" if period < 10 else "long"


def mag_band(tmag: float) -> str:
    return "bright" if tmag < 12 else "faint"


def choose(rows: list[dict], excluded_hosts: set[str], limit: int, seed: int) -> list[dict]:
    strata: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        if row["host_id"] in excluded_hosts:
            continue
        # Exclude the hottest stars for a more representative transit-vetting
        # sample while retaining a broad magnitude and period range.
        if not 2700 <= row["teff_k"] <= 10000 or row["tmag"] > 15:
            continue
        band = f"{period_band(row['period_days'])}:{mag_band(row['tmag'])}"
        key = hashlib.sha256(f"{seed}:{row['tic_id']}".encode()).hexdigest()
        strata[band].append({**row, "selection_hash": key, "stratum": band})
    for values in strata.values():
        values.sort(key=lambda row: row["selection_hash"])
    selected = []
    keys = sorted(strata)
    cursor = {key: 0 for key in keys}
    while len(selected) < limit:
        advanced = False
        for key in keys:
            if cursor[key] < len(strata[key]):
                selected.append(strata[key][cursor[key]])
                cursor[key] += 1
                advanced = True
                if len(selected) == limit:
                    break
        if not advanced:
            break
    return selected


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("selection returned no targets")
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path,
                        default=root / "data/benchmarks/hard-negatives/tess-eb-catalog-table3.tsv")
    parser.add_argument("--known-answer-metadata", type=Path,
                        default=root / "data/benchmarks/generated-20261009/case-metadata.csv")
    parser.add_argument("--known-planets", type=Path,
                        default=root / "data/benchmarks/generated-20261009/nasa-pscomppars-tess-transiting.csv")
    parser.add_argument("--output-dir", type=Path,
                        default=root / "data/benchmarks/hard-negatives")
    parser.add_argument("--limit", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20261010)
    args = parser.parse_args()
    if args.limit < 1:
        parser.error("--limit must be positive")
    catalog = args.catalog.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    known_answer_hosts = set()
    if args.known_answer_metadata.is_file():
        with args.known_answer_metadata.open(newline="", encoding="utf-8-sig") as stream:
            known_answer_hosts = {row.get("tic_id", "").strip()
                                  for row in csv.DictReader(stream) if row.get("tic_id")}
    catalog_planet_hosts = set()
    if args.known_planets.is_file():
        with args.known_planets.open(newline="", encoding="utf-8-sig") as stream:
            for row in csv.DictReader(stream):
                tic = row.get("tic_id", "").strip().upper().removeprefix("TIC").strip()
                if tic.isdigit():
                    catalog_planet_hosts.add(tic)
    excluded = known_answer_hosts | catalog_planet_hosts
    rows = read_catalog(catalog)
    selected = choose(rows, excluded, args.limit, args.seed)
    if len(selected) < args.limit:
        parser.error(f"only {len(selected)} eligible non-overlapping catalog targets")
    write_csv(output_dir / "selected-eb-targets.csv", selected)
    target_rows = [{"tic_id": r["tic_id"], "teff_k": r["teff_k"],
                    "tmag": r["tmag"], "period_days": r["period_days"],
                    "epoch_days": r["epoch_days"], "duration_days": r["duration_days"],
                    "radius_solar": "", "mass_solar": ""}
                   for r in selected]
    write_csv(output_dir / "survey-targets.csv", target_rows)
    report = {
        "catalog": "VizieR J/ApJS/279/50/table3",
        "catalog_reference": SOURCE_REFERENCE,
        "catalog_url": "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJS/279/50/table3&-out.max=10000",
        "catalog_sha256": sha256(catalog), "selection_seed": args.seed,
        "known_planets_catalog": str(args.known_planets.resolve()) if args.known_planets.is_file() else None,
        "known_planets_catalog_sha256": sha256(args.known_planets.resolve()) if args.known_planets.is_file() else None,
        "known_answer_metadata_sha256": sha256(args.known_answer_metadata.resolve()) if args.known_answer_metadata.is_file() else None,
        "selected_targets": len(selected), "catalog_rows_eligible": len(rows),
        "known_answer_hosts_excluded": len(known_answer_hosts),
        "NASA_TESS_planet_hosts_excluded": len(catalog_planet_hosts),
        "unique_hosts_excluded": len(excluded),
        "criteria": {"period_days": [0.15, 30], "teff_k": [2700, 10000],
                     "tmag_max": 15, "strata": "period short/medium/long x Tmag bright/faint",
                     "exclude_hosts_in_full_frozen_TESS_planet_catalog": True,
                     "split": "all held out as blind hard negatives; threshold calibrated on independent known-answer inverted controls"},
        "selected_by_stratum": dict(sorted(Counter(r["stratum"] for r in selected).items())),
        "host_ids": [r["host_id"] for r in selected],
    }
    (output_dir / "selection.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "host_ids"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
