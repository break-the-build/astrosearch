"""Download nearby TESS-SPOC FFI light curves for null and injection calibration.

usage: fetch_calibration_curves.py OUT_DIR TARGET SECTOR EXPTIME RADIUS_ARCSEC COUNT

Takes the COUNT nearest stars to TARGET (excluding TARGET) that have a TESS-SPOC
light curve in SECTOR at EXPTIME seconds. Stars share the target's camera and
systematics. Each curve is quality-masked ("hard"), NaN-free, and normalized
to its median, matching export_tess.py. Writes OUT_DIR/manifest.csv.
"""
import csv
import hashlib
import sys
from pathlib import Path

import lightkurve as lk
import numpy as np


def export(lc: lk.LightCurve, path: Path) -> int:
    lc = lc.remove_nans().normalize()
    rows = np.column_stack([
        np.asarray(lc.time.value, float),
        np.ma.filled(lc.flux.value, np.nan),
        np.ma.filled(lc.flux_err.value, np.nan),
    ])
    rows = rows[np.all(np.isfinite(rows), axis=1) & (rows[:, 1] > 0)]
    np.savetxt(path, rows, delimiter=",", header="time,flux,flux_err", comments="", fmt="%.8f")
    return len(rows)


def main() -> None:
    out_dir = Path(sys.argv[1])
    target, sector, exptime = sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    radius, count = float(sys.argv[5]), int(sys.argv[6])
    out_dir.mkdir(parents=True, exist_ok=True)
    result = lk.search_lightcurve(target, radius=radius, sector=sector,
                                  author="TESS-SPOC", exptime=exptime)
    table = result.table
    order = np.argsort(np.asarray(table["distance"], float), kind="stable")
    central = target.replace("TIC", "").strip()
    manifest, seen = [], set()
    for index in order:
        name = str(table["target_name"][index]).strip()
        if name == central or name in seen:
            continue
        seen.add(name)
        path = out_dir / f"tic{name}_s{sector:02d}_{exptime}.csv"
        try:
            lc = result[int(index)].download(quality_bitmask="hard",
                                             download_dir=str(out_dir / "cache"))
            rows = export(lc, path)
        except Exception as error:  # a broken MAST product should not stop the set
            print(f"skip TIC {name}: {error}", file=sys.stderr)
            continue
        if rows < 500:
            path.unlink()
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        manifest.append({"tic": name, "sector": sector, "exptime": exptime, "rows": rows,
                         "distance_arcsec": float(table["distance"][index]),
                         "file": path.name, "sha256": digest})
        print(f"{len(manifest):3d} TIC {name}: {rows} rows")
        if len(manifest) >= count:
            break
    with (out_dir / "manifest.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)


if __name__ == "__main__":
    main()
