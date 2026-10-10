"""Is the candidate already known? TOI, CTOI, confirmed planets, SPOC TCEs.

Reuses the NASA Exoplanet Archive and MAST helpers in candidate_preflight.py
and adds the ExoFOP CTOI table. An empty or failed query is reported as such,
never as evidence of novelty.
"""
from __future__ import annotations

import csv
import io as _io

import candidate_preflight as preflight

CTOI_URL = "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv"


def candidate_row(period: float, epoch: float, duration: float) -> dict:
    return {"period_days": period, "epoch": epoch, "duration_days": duration, "class_rank": "1"}


def query_ctoi(tic_id: int) -> list[dict]:
    response = preflight.SESSION.get(CTOI_URL, timeout=90)
    response.raise_for_status()
    rows = list(csv.DictReader(_io.StringIO(response.text)))
    key = next((k for k in rows[0] if k.lower().replace(" ", "") in ("ticid", "tic")), None) if rows else None
    return [r for r in rows if key and str(r.get(key, "")).strip() == str(tic_id)]


def ctoi_matches(rows: list[dict], candidate: dict) -> list[dict]:
    matches = []
    for row in rows:
        lowered = {k.lower(): v for k, v in row.items()}
        try:
            period = float(lowered.get("period (days)") or lowered.get("period") or "nan")
            epoch = float(lowered.get("epoch (bjd)") or lowered.get("epoch") or "nan")
            duration = float(lowered.get("duration (hours)") or 0) / 24
        except ValueError:
            continue
        if epoch > 1_000_000:
            epoch -= 2_457_000
        for match in preflight.ephemeris_matches([candidate], period, epoch, duration):
            matches.append({"ctoi": lowered.get("ctoi"), **match})
    return matches


def known_object_check(tic_id: int, period: float, epoch: float, duration: float,
                       tce_sectors: list[int]) -> dict:
    candidate = candidate_row(period, epoch, duration)
    report, errors = {}, {}

    def attempt(name, fn):
        try:
            report[name] = fn()
        except Exception as error:  # network/schema failures are reported, not hidden
            errors[name] = f"{type(error).__name__}: {error}"

    attempt("toi", lambda: preflight.tap(
        "select toi,ctoi_alias,tid,tfopwg_disp,pl_orbper,pl_tranmid,pl_trandurh,pl_trandep "
        f"from toi where tid={int(tic_id)}"))
    attempt("confirmed", lambda: preflight.tap(
        "select pl_name,tic_id,pl_orbper,pl_tranmid,pl_trandur as pl_trandurh,discoverymethod "
        f"from ps where tic_id='TIC {int(tic_id)}'"))
    attempt("ctoi", lambda: query_ctoi(tic_id))
    if tce_sectors:
        attempt("spoc_tce", lambda: preflight.check_spoc_tces(tic_id, [candidate], tce_sectors))
    matches = {
        "toi": preflight.match_catalog_ephemerides(report.get("toi"), [candidate]),
        "confirmed": preflight.match_catalog_ephemerides(report.get("confirmed"), [candidate]),
        "ctoi": ctoi_matches(report.get("ctoi") or [], candidate),
        "spoc_tce": [m for m in (report.get("spoc_tce") or {}).get("matches", [])
                     if m.get("ephemeris_matches")],
    }
    known = any(matches.values())
    complete = not errors and "toi" in report and "confirmed" in report and "ctoi" in report
    status = "FAIL" if known else ("PASS" if complete else "NOT_EVALUATED")
    return {"status": status, "matches": matches, "errors": errors,
            "rows_returned": {k: len(v) if isinstance(v, list) else None for k, v in report.items()},
            "tce_sectors_checked": tce_sectors,
            "note": "PASS means no ephemeris match in these catalogs on this date, not novelty"}
