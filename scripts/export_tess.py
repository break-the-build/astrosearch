"""Download TESS light curves and export astrosearch CSVs with provenance sidecars.

usage: export_tess.py OUT_DIR TARGET SECTOR:AUTHOR:EXPTIME [...]

Each sector is quality-masked ("hard"), NaN-free, and normalized to its own
median; OUT_DIR/combined.csv is the time-sorted concatenation. Every CSV gets a
`<name>.meta.json` sidecar recording the product, flux column, crowding
keywords (whether the pipeline already corrected dilution), and the times of
momentum dumps and scattered-light flags, which transit vetting checks against.
"""
import json
import sys
from pathlib import Path

import lightkurve as lk

from vetting.network import install_default_timeouts
import numpy as np

# SPOC quality bits (TESS Science Data Products Description Document).
MOMENTUM_DUMP = 32
SCATTERED_LIGHT = 2048 | 4096
SPOC_AUTHORS = {"SPOC", "TESS-SPOC"}
HEADER_KEYS = ["CROWDSAP", "FLFRCSAP", "TICID", "SECTOR", "CAMERA", "CCD", "TESSMAG",
               "TEFF", "RADIUS", "RA_OBJ", "DEC_OBJ"]


def intervals(times: np.ndarray, cadence: float) -> list[list[float]]:
    """Merge flagged cadences closer than two cadences into [start, end] spans."""
    spans: list[list[float]] = []
    for t in np.sort(times):
        if spans and t - spans[-1][1] <= 2 * cadence:
            spans[-1][1] = float(t)
        else:
            spans.append([float(t), float(t)])
    return spans


def quality_events(raw: lk.LightCurve, author: str) -> dict:
    """Momentum dumps and scattered-light spans from the unmasked light curve."""
    if author not in SPOC_AUTHORS:
        return {"quality_convention": f"{author}: not interpreted"}
    time = np.asarray(raw.time.value, float)
    quality = np.asarray(raw.quality.value if hasattr(raw.quality, "value") else raw.quality, int)
    cadence = float(np.median(np.diff(time))) if len(time) > 1 else 0.0
    return {
        "quality_convention": "SPOC bits: 32 momentum dump, 2048|4096 scattered light",
        "momentum_dump_times": [float(t) for t in time[(quality & MOMENTUM_DUMP) != 0]],
        "scattered_light_intervals": intervals(time[(quality & SCATTERED_LIGHT) != 0], cadence),
    }


def metadata(lc: lk.LightCurve, raw: lk.LightCurve, target: str, sector: int, author: str,
             exptime: int, rows: np.ndarray) -> dict:
    header = {key: lc.meta.get(key) for key in HEADER_KEYS if lc.meta.get(key) is not None}
    flux_origin = str(lc.meta.get("FLUX_ORIGIN", "")).lower()
    # SPOC-family PDCSAP flux has the CROWDSAP dilution correction applied;
    # other products are flagged unknown so vetting applies the TIC ratio.
    corrected = (author in SPOC_AUTHORS and flux_origin.startswith("pdcsap")
                 and "CROWDSAP" in header)
    return {
        "target": target, "sector": sector, "author": author, "exptime_s": exptime,
        "flux_origin": flux_origin, "crowding_corrected": corrected if corrected else None,
        "header": {k: (float(v) if isinstance(v, (int, float, np.floating)) else str(v))
                   for k, v in header.items()},
        "rows": int(len(rows)), "time_min": float(rows[0, 0]), "time_max": float(rows[-1, 0]),
        "median_flux_err": float(np.median(rows[:, 2])),
        **quality_events(raw, author),
    }


def export(out_dir: Path, target: str, spec: str):
    sector, author, exptime = spec.split(":")
    result = lk.search_lightcurve(target, sector=int(sector), author=author, exptime=int(exptime))
    if not len(result):
        print("missing", spec)
        return None
    lc = result[0].download(quality_bitmask="hard", download_dir=str(out_dir / "cache"))
    raw = result[0].download(quality_bitmask="none", download_dir=str(out_dir / "cache"))
    lc = lc.remove_nans().normalize()
    lc = lc[np.isfinite(lc.flux_err.value) & (lc.flux.value > 0)]
    rows = np.column_stack([np.asarray(lc.time.value, float), np.ma.filled(lc.flux.value, np.nan),
                            np.ma.filled(lc.flux_err.value, np.nan)])
    rows = rows[np.all(np.isfinite(rows), axis=1)]
    name = out_dir / f"s{int(sector):02d}_{author}_{exptime}.csv"
    np.savetxt(name, rows, delimiter=",", header="time,flux,flux_err", comments="", fmt="%.8f")
    meta = metadata(lc, raw, target, int(sector), author, int(exptime), rows)
    meta["file"] = name.name
    name.with_suffix(".meta.json").write_text(json.dumps(meta, indent=1, sort_keys=True) + "\n")
    print(name.name, len(rows), "rows, median err ppm", round(np.median(rows[:, 2]) * 1e6))
    return rows, meta


def main() -> None:
    install_default_timeouts()
    out_dir, target = Path(sys.argv[1]), sys.argv[2]
    out_dir.mkdir(parents=True, exist_ok=True)
    combined, metas = [], []
    for spec in sys.argv[3:]:
        exported = export(out_dir, target, spec)
        if exported:
            combined.append(exported[0])
            metas.append(exported[1])
    rows = np.vstack(combined)
    rows = rows[np.argsort(rows[:, 0])]
    np.savetxt(out_dir / "combined.csv", rows, delimiter=",", header="time,flux,flux_err",
               comments="", fmt="%.8f")
    (out_dir / "combined.meta.json").write_text(json.dumps(
        {"target": target, "parts": metas, "rows": int(len(rows))}, indent=1, sort_keys=True) + "\n")
    print("combined", len(rows))


if __name__ == "__main__":
    main()
