#!/usr/bin/env python3
"""Benchmark search recovery and vetting on known planets and known eclipsing binaries.

  build  Query the NASA Exoplanet Archive for literature ephemerides of the
         listed planets and for deep TFOPWG false-positive (FP) TOIs, mostly eclipsing binaries, and write a
         manifest (so ephemerides come from the archive, not from memory).
  run    For each target, export SPOC 2-minute light curves, search them with
         astrosearch, find the rank of each known planet, and run the CTOI gate
         without nulls or FPP. Planets should pass criteria 2 and 3 and fail 4
         (they are known); false positives should fail 2 or 3.

usage: benchmark_known_systems.py build MANIFEST.csv
       benchmark_known_systems.py run MANIFEST.csv OUT_DIR [--binary B] [--sectors N]
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import candidate_preflight as preflight
from vetting import network

ROOT = Path(__file__).resolve().parents[1]
PLANETS = ["LHS 3844 b", "L 98-59 b", "L 98-59 c", "L 98-59 d", "TOI-700 b", "TOI-700 c"]
FALSE_POSITIVE_COUNT = 4


def build(manifest: Path) -> None:
    names = ",".join(f"'{n}'" for n in PLANETS)
    planets = preflight.tap(
        "select pl_name,tic_id,pl_orbper,pl_tranmid,pl_trandur as pl_trandurh,pl_trandep from pscomppars "
        f"where pl_name in ({names})")
    ebs = preflight.tap(
        "select top 40 toi,tid,tfopwg_disp,pl_orbper,pl_tranmid,pl_trandurh,pl_trandep,st_tmag "
        "from toi where tfopwg_disp='FP' and st_tmag < 10.5 and pl_orbper between 1 and 8 "
        "and pl_trandep > 1000 "
        "order by toi asc")
    rows = []
    for p in planets:
        rows.append({"label": p["pl_name"], "kind": "planet",
                     "tic": str(p["tic_id"]).replace("TIC", "").strip(),
                     "period_days": p["pl_orbper"], "epoch_btjd": float(p["pl_tranmid"]) - 2_457_000,
                     "duration_hours": p["pl_trandurh"], "source": "NASA Exoplanet Archive pscomppars"})
    for e in ebs[:FALSE_POSITIVE_COUNT]:
        rows.append({"label": f"TOI-{e['toi']}", "kind": "false_positive", "tic": str(e["tid"]),
                     "period_days": e["pl_orbper"], "epoch_btjd": float(e["pl_tranmid"]) - 2_457_000,
                     "duration_hours": e["pl_trandurh"],
                     "source": "NASA Exoplanet Archive TOI table, tfopwg_disp=FP, depth > 1000 ppm"})
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {manifest}: {len(rows)} entries")


def product_specs(tic: str, author: str, count: int, exptime: int | None = None,
                  sectors: list[int] | None = None) -> list[str]:
    """`sector:AUTHOR:exptime` export specs for up to `count` sectors (all of
    `sectors` when given) with light curves from `author`; [] when none exist."""
    import lightkurve as lk
    query = {"author": author, **({"exptime": exptime} if exptime else {})}
    result = lk.search_lightcurve(f"TIC {tic}", **query)
    if len(result) == 0 or "mission" not in result.table.colnames:
        return []
    found = {}
    for mission, seconds in zip(result.table["mission"], result.table["exptime"]):
        word = str(mission).split()[-1]
        if word.isdigit():
            found.setdefault(int(word), int(float(seconds)))
    chosen = [s for s in sorted(found) if sectors is None or s in sectors]
    if sectors is None:
        chosen = chosen[:count]
    return [f"{s}:{author}:{found[s]}" for s in chosen]


def export_product(python: str, directory: Path, tic: str, specs: list[str]) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    subprocess.run([python, str(ROOT / "scripts/export_tess.py"), str(directory), f"TIC {tic}", *specs],
                   check=True, capture_output=True, text=True)
    return directory / "combined.csv"


def match_rank(candidates: Path, period: float, epoch: float, duration_days: float):
    with candidates.open(newline="") as stream:
        rows = [r for r in csv.DictReader(stream) if r["class"] in ("periodic_dip", "duo_transit")]
    for position, row in enumerate(rows, start=1):
        if preflight.ephemeris_matches([row], period, epoch, duration_days):
            return position, int(row["rank"]), row["class"]
    return None, None, None


def gate_record(report: dict) -> dict:
    record = {key: value["status"] for key, value in report["criteria"].items()}
    checks = report["criteria"]["3_no_false_positive_indicators"].get("checks", {})
    record.update({f"fp_{key}": value for key, value in checks.items()})
    record.update(secondary=report["secondary"]["status"], odd_even=report["odd_even"]["status"],
                  density=report["density"].get("status"), bls=report["bls"]["status"],
                  ruwe=report.get("ruwe"),
                  siblings_masked=len((report.get("siblings") or {}).get("masked", [])),
                  products=(report.get("products") or {}).get("status"))
    return record


def run_target(tic: str, entries: list[dict], out: Path, binary: Path, sector_count: int,
               python: str, second_author: str | None, skip_search: bool) -> list[dict]:
    target_dir = out / f"tic{tic}"
    target_dir.mkdir(parents=True, exist_ok=True)
    candidates = target_dir / "candidates.csv"
    specs = product_specs(tic, "SPOC", sector_count, exptime=120)
    if not specs:
        return [{"label": e["label"], "kind": e["kind"], "tic": tic, "status": "no SPOC 2-min data"}
                for e in entries]
    if not (skip_search and candidates.exists()):
        export_product(python, target_dir, tic, specs)
        subprocess.run([str(binary), "--input", str(target_dir / "combined.csv"),
                        "--output", str(candidates), "--min-period", "0.3", "--max-period", "20"],
                       check=True, capture_output=True, text=True)
    products = []
    if second_author:
        sectors = [int(s.split(":")[0]) for s in specs]
        extra = product_specs(tic, second_author, sector_count, sectors=sectors)
        second = target_dir / second_author.lower() / "combined.csv"
        if extra and not (skip_search and second.exists()):
            export_product(python, second.parent, tic, extra)
        products = [second] if second.exists() else []
    records = []
    for entry in entries:
        period, epoch = float(entry["period_days"]), float(entry["epoch_btjd"])
        duration = float(entry["duration_hours"] or 2) / 24
        position, rank, cls = match_rank(candidates, period, epoch, duration)
        record = {"label": entry["label"], "kind": entry["kind"], "tic": tic,
                  "sectors": " ".join(specs), "transit_rank": position, "class": cls,
                  "second_product": str(products[0]) if products else None}
        if rank is not None:
            gate_dir = target_dir / f"gate_{entry['label'].replace(' ', '_')}"
            done = subprocess.run(
                [python, str(ROOT / "scripts/ctoi_gate.py"), "--tic", tic,
                 "--lightcurve", str(target_dir / "combined.csv"), "--candidates", str(candidates),
                 "--rank", str(rank), "--binary", str(binary), "--null-trials", "0", "--skip-fpp",
                 "--mcmc-steps", "0", "--out-dir", str(gate_dir),
                 *(["--products", *map(str, products)] if products else [])],
                check=False, capture_output=True, text=True)
            report_path = gate_dir / "gate.json"
            if done.returncode == 0 and report_path.exists():
                record.update(gate_record(json.loads(report_path.read_text())))
            else:
                record["gate_error"] = done.stderr.strip().splitlines()[-1:] or ["no gate.json"]
        records.append(record)
    return records


def run(manifest: Path, out: Path, binary: Path, sector_count: int, python: str,
        second_author: str | None = None, skip_search: bool = False) -> None:
    with manifest.open(newline="") as stream:
        entries = list(csv.DictReader(stream))
    results = []
    for tic in dict.fromkeys(e["tic"] for e in entries):
        mine = [e for e in entries if e["tic"] == tic]
        try:
            records = run_target(tic, mine, out, binary, sector_count, python, second_author, skip_search)
        except (subprocess.CalledProcessError, OSError, ValueError) as error:
            detail = getattr(error, "stderr", "") or str(error)
            records = [{"label": e["label"], "kind": e["kind"], "tic": tic,
                        "status": f"error: {type(error).__name__}: {detail.strip()[-300:]}"}
                       for e in mine]
        for record in records:
            print(json.dumps(record))
        results.extend(records)
    summary = out / "summary.csv"
    fields = sorted({k for r in results for k in r})
    with summary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)
    print(f"wrote {summary}")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)
    b = sub.add_parser("build")
    b.add_argument("manifest", type=Path)
    r = sub.add_parser("run")
    r.add_argument("manifest", type=Path)
    r.add_argument("out", type=Path)
    r.add_argument("--binary", type=Path, default=ROOT / "target/release/astrosearch")
    r.add_argument("--sectors", type=int, default=2)
    r.add_argument("--python", default=sys.executable)
    r.add_argument("--second-author", default="QLP",
                   help="independent product exported for the same sectors (empty to skip)")
    r.add_argument("--skip-search", action="store_true",
                   help="reuse existing exports and candidates; rerun only the gate")
    a = p.parse_args()
    network.install_default_timeouts()
    if a.command == "build":
        build(a.manifest)
    else:
        run(a.manifest, a.out, a.binary, a.sectors, a.python, a.second_author or None, a.skip_search)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
