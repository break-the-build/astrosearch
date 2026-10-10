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

sys.path.insert(0, str(Path(__file__).resolve().parent))
import candidate_preflight as preflight  # noqa: E402

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


def spoc_sectors(tic: str, count: int) -> list[str]:
    import lightkurve as lk
    result = lk.search_lightcurve(f"TIC {tic}", author="SPOC", exptime=120)
    sectors = sorted({int(m.split()[-1]) for m in result.table["mission"] if m.split()[-1].isdigit()})
    return [f"{s}:SPOC:120" for s in sectors[:count]]


def match_rank(candidates: Path, period: float, epoch: float, duration_days: float):
    with candidates.open(newline="") as stream:
        rows = [r for r in csv.DictReader(stream) if r["class"] in ("periodic_dip", "duo_transit")]
    for position, row in enumerate(rows, start=1):
        if preflight.ephemeris_matches([row], period, epoch, duration_days):
            return position, int(row["rank"]), row["class"]
    return None, None, None


def run(manifest: Path, out: Path, binary: Path, sector_count: int, python: str) -> None:
    with manifest.open(newline="") as stream:
        entries = list(csv.DictReader(stream))
    results = []
    for tic in dict.fromkeys(e["tic"] for e in entries):
        target_dir = out / f"tic{tic}"
        target_dir.mkdir(parents=True, exist_ok=True)
        specs = spoc_sectors(tic, sector_count)
        if not specs:
            results.extend({"label": e["label"], "status": "no SPOC 2-min data"}
                           for e in entries if e["tic"] == tic)
            continue
        subprocess.run([python, str(ROOT / "scripts/export_tess.py"), str(target_dir), f"TIC {tic}",
                        *specs], check=True, capture_output=True, text=True)
        candidates = target_dir / "candidates.csv"
        subprocess.run([str(binary), "--input", str(target_dir / "combined.csv"),
                        "--output", str(candidates), "--min-period", "0.3", "--max-period", "20"],
                       check=True, capture_output=True, text=True)
        for entry in (e for e in entries if e["tic"] == tic):
            period, epoch = float(entry["period_days"]), float(entry["epoch_btjd"])
            duration = float(entry["duration_hours"] or 2) / 24
            position, rank, cls = match_rank(candidates, period, epoch, duration)
            record = {"label": entry["label"], "kind": entry["kind"], "tic": tic,
                      "sectors": " ".join(specs), "transit_rank": position, "class": cls}
            if rank is not None:
                gate_dir = target_dir / f"gate_{entry['label'].replace(' ', '_')}"
                subprocess.run([python, str(ROOT / "scripts/ctoi_gate.py"), "--tic", tic,
                                "--lightcurve", str(target_dir / "combined.csv"),
                                "--candidates", str(candidates), "--rank", str(rank),
                                "--binary", str(binary), "--null-trials", "0", "--skip-fpp",
                                "--mcmc-steps", "0", "--out-dir", str(gate_dir)],
                               check=False, capture_output=True, text=True)
                report_path = gate_dir / "gate.json"
                if report_path.exists():
                    report = json.loads(report_path.read_text())
                    record.update({key: value["status"] for key, value in report["criteria"].items()})
                    record["secondary"] = report["secondary"]["status"]
                    record["odd_even"] = report["odd_even"]["status"]
                    record["density"] = report["density"].get("status")
                    record["bls"] = report["bls"]["status"]
            results.append(record)
            print(json.dumps(record))
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
    a = p.parse_args()
    if a.command == "build":
        build(a.manifest)
    else:
        run(a.manifest, a.out, a.binary, a.sectors, a.python)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
