"""Persistent, cross-run bookkeeping for TESS survey targets."""
from __future__ import annotations

import hashlib
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _digest(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def metadata_digest(row: dict[str, Any]) -> str:
    aliases = {
        "tic_id": ("tic_id", "tic", "ID"),
        "teff_k": ("teff_k", "teff", "Teff"),
        "distance_pc": ("distance_pc", "distance", "d"),
        "radius_solar": ("radius_solar", "rad", "radius"),
        "mass_solar": ("mass_solar", "mass", "Mstar", "mstar"),
        "tmag": ("tmag", "Tmag"),
        "stellar_density_g_cm3": ("stellar_density_g_cm3", "rho_g_cm3", "rho"),
    }
    canonical = {}
    for name, keys in aliases.items():
        canonical[name] = next((row[key] for key in keys if row.get(key) not in (None, "")), None)
    return _digest(canonical)


def config_digest(args: Any) -> str:
    """Only settings that change the resulting data/search are fingerprinted."""
    profile = getattr(args, "search_profile", "standard")
    # The standard profile adds no key, so existing ledgers keep their digests.
    extra = {} if profile == "standard" else {"search_profile": profile}
    return _digest({**extra,
        "authors": list(args.authors),
        "min_sector": args.min_sector,
        "min_exptime": args.min_exptime,
        "min_baseline": args.min_baseline,
        "min_teff": args.min_teff,
        "max_teff": args.max_teff,
        "max_distance_pc": args.max_distance_pc,
        "max_radius": args.max_radius,
        "max_tmag": args.max_tmag,
        "max_sector_noise_ratio": args.max_sector_noise_ratio,
        "min_period": args.min_period,
        "max_period": args.max_period,
        "periods": args.periods,
        "threads": args.threads,
        "trend_method": args.trend_method,
        "trend_days": getattr(args, "trend_days", 1.05),
        "segment_gap_days": getattr(args, "segment_gap_days", 0.5),
        "transit_search_passes": getattr(args, "transit_search_passes", 5),
        "candidate_grouping": args.candidate_grouping,
        "tls": getattr(args, "tls", False),
        "tls_max_planets": getattr(args, "tls_max_planets", 5),
        "tls_min_sde": getattr(args, "tls_min_sde", 7.0),
    })


def connect(path: Path) -> sqlite3.Connection:
    path = path.expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(path, timeout=30)
    db.row_factory = sqlite3.Row
    db.execute("PRAGMA journal_mode=WAL")
    db.execute("PRAGMA busy_timeout=30000")
    db.executescript("""
        CREATE TABLE IF NOT EXISTS survey_runs (
            run_id TEXT PRIMARY KEY,
            started_utc TEXT NOT NULL,
            finished_utc TEXT,
            targets_path TEXT NOT NULL,
            decision_engine TEXT NOT NULL,
            config_sha256 TEXT NOT NULL,
            status TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS target_state (
            tic_id TEXT PRIMARY KEY,
            status TEXT NOT NULL,
            last_run_id TEXT NOT NULL,
            updated_utc TEXT NOT NULL,
            metadata_sha256 TEXT NOT NULL,
            config_sha256 TEXT NOT NULL,
            output_path TEXT,
            result_json TEXT NOT NULL,
            FOREIGN KEY(last_run_id) REFERENCES survey_runs(run_id)
        );
        CREATE INDEX IF NOT EXISTS target_state_updated_idx
            ON target_state(updated_utc);
    """)
    db.commit()
    return db


def begin_run(db: sqlite3.Connection, run_id: str, targets_path: Path,
              decision_engine: str, config_sha: str) -> None:
    db.execute(
        "INSERT INTO survey_runs(run_id,started_utc,targets_path,decision_engine,config_sha256,status) "
        "VALUES(?,?,?,?,?, 'running')",
        (run_id, datetime.now(timezone.utc).isoformat(), str(targets_path.resolve()),
         decision_engine, config_sha),
    )
    db.commit()


def should_defer(db: sqlite3.Connection, tic_id: str, row: dict[str, Any],
                 config_sha: str, cooldown_days: int, revisit: bool) -> bool:
    if revisit or cooldown_days <= 0:
        return False
    saved = db.execute("SELECT * FROM target_state WHERE tic_id=?", (tic_id,)).fetchone()
    if not saved or saved["status"] not in (
        "processed", "imported", "completed_no_data", "insufficient_multiyear_baseline"
    ):
        return False
    if saved["metadata_sha256"] != metadata_digest(row) or saved["config_sha256"] != config_sha:
        return False
    try:
        updated = datetime.fromisoformat(saved["updated_utc"])
    except ValueError:
        return False
    age = datetime.now(timezone.utc) - updated.astimezone(timezone.utc)
    return age.total_seconds() < cooldown_days * 86400


def claim_target(db: sqlite3.Connection, run_id: str, tic_id: str,
                 row: dict[str, Any], config_sha: str,
                 lease_hours: float = 24.0) -> bool:
    """Atomically claim a target; recover an abandoned claim after its lease."""
    db.execute("BEGIN IMMEDIATE")
    saved = db.execute("SELECT status, updated_utc FROM target_state WHERE tic_id=?",
                       (tic_id,)).fetchone()
    if saved and saved["status"] == "running":
        try:
            updated = datetime.fromisoformat(saved["updated_utc"])
            age = datetime.now(timezone.utc) - updated.astimezone(timezone.utc)
        except (ValueError, TypeError):
            age = None
        if age is not None and age.total_seconds() < lease_hours * 3600:
            db.rollback()
            return False
    now = datetime.now(timezone.utc).isoformat()
    db.execute(
        "INSERT INTO target_state(tic_id,status,last_run_id,updated_utc,metadata_sha256,config_sha256,output_path,result_json) "
        "VALUES(?, 'running', ?, ?, ?, ?, NULL, '{}') "
        "ON CONFLICT(tic_id) DO UPDATE SET status='running', last_run_id=excluded.last_run_id, "
        "updated_utc=excluded.updated_utc, metadata_sha256=excluded.metadata_sha256, "
        "config_sha256=excluded.config_sha256, output_path=NULL, result_json='{}'",
        (tic_id, run_id, now, metadata_digest(row), config_sha),
    )
    db.commit()
    return True


def record_target(db: sqlite3.Connection, run_id: str, tic_id: str,
                  row: dict[str, Any], config_sha: str,
                  result: dict[str, Any], status: str | None = None) -> bool:
    """Write a target's result if this run still holds its claim.

    Returns False, writing nothing, when another run recovered the claim
    after this run's lease expired; the newer run's state must not be
    overwritten by a stale result."""
    now = datetime.now(timezone.utc).isoformat()
    final_status = status or str(result.get("status", "error"))
    values = (final_status, now, metadata_digest(row), config_sha, result.get("output"),
              json.dumps(result, sort_keys=True, default=str))
    db.execute("BEGIN IMMEDIATE")
    updated = db.execute(
        "UPDATE target_state SET status=?, updated_utc=?, metadata_sha256=?, config_sha256=?, "
        "output_path=?, result_json=? WHERE tic_id=? AND last_run_id=?",
        (*values, tic_id, run_id)).rowcount
    if updated == 0:
        exists = db.execute("SELECT 1 FROM target_state WHERE tic_id=?", (tic_id,)).fetchone()
        if exists:
            db.rollback()
            return False
        db.execute(
            "INSERT INTO target_state(status,updated_utc,metadata_sha256,config_sha256,output_path,"
            "result_json,tic_id,last_run_id) VALUES(?,?,?,?,?,?,?,?)", (*values, tic_id, run_id))
    db.commit()
    return True


def finish_run(db: sqlite3.Connection, run_id: str, status: str = "complete") -> None:
    if status in ("failed", "partial", "interrupted"):
        db.execute(
            "UPDATE target_state SET status='error', updated_utc=?, "
            "result_json=? "
            "WHERE last_run_id=? AND status='running'",
            (datetime.now(timezone.utc).isoformat(),
             json.dumps({"status": "error", "reason": "run interrupted before target commit"}),
             run_id),
        )
    db.execute("UPDATE survey_runs SET finished_utc=?, status=? WHERE run_id=?",
               (datetime.now(timezone.utc).isoformat(), status, run_id))
    db.commit()
