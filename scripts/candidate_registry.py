#!/usr/bin/env python3
"""Append-only registry for candidate detections and human-reviewed status changes."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

STATUSES = ("DETECTED_SIGNAL", "VETTED_CANDIDATE",
            "INDEPENDENTLY_VALIDATED_CANDIDATE", "CONFIRMED_DISCOVERY")
NUMERIC_FIELDS = ("class_rank", "score", "jackknife_score", "white_score", "noise_factor",
                  "period_days", "epoch", "duration_days", "depth_ppm",
                  "amplitude_ppm", "n_events", "odd_even_sigma")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def optional_number(value: str | None):
    if value in (None, ""):
        return None
    try:
        result = float(value)
    except ValueError:
        return None
    return result if result == result and abs(result) != float("inf") else None


def initialize(connection: sqlite3.Connection) -> None:
    connection.execute("PRAGMA foreign_keys=ON")
    connection.executescript("""
        CREATE TABLE IF NOT EXISTS runs (
            run_id TEXT PRIMARY KEY,
            target_id TEXT NOT NULL,
            created_utc TEXT NOT NULL,
            source_revision TEXT NOT NULL,
            input_sha256 TEXT NOT NULL,
            candidates_sha256 TEXT NOT NULL,
            provenance_json TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS detections (
            detection_id TEXT PRIMARY KEY,
            run_id TEXT NOT NULL REFERENCES runs(run_id),
            candidate_rank INTEGER,
            candidate_class TEXT NOT NULL,
            class_rank INTEGER,
            score REAL,
            jackknife_score REAL,
            white_score REAL,
            noise_factor REAL,
            period_days REAL,
            epoch REAL,
            duration_days REAL,
            depth_ppm REAL,
            amplitude_ppm REAL,
            n_events INTEGER,
            odd_even_sigma REAL,
            tier TEXT,
            event_evidence TEXT,
            alias_ephemerides TEXT,
            note TEXT,
            raw_candidate_json TEXT NOT NULL,
            status TEXT NOT NULL CHECK(status IN (
                'DETECTED_SIGNAL','VETTED_CANDIDATE',
                'INDEPENDENTLY_VALIDATED_CANDIDATE','CONFIRMED_DISCOVERY')),
            created_utc TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS status_history (
            detection_id TEXT NOT NULL REFERENCES detections(detection_id),
            previous_status TEXT,
            new_status TEXT NOT NULL,
            changed_utc TEXT NOT NULL,
            actor TEXT NOT NULL,
            evidence TEXT NOT NULL,
            PRIMARY KEY(detection_id, new_status)
        );
    """)
    columns = {row[1] for row in connection.execute("PRAGMA table_info(detections)")}
    for name in ("jackknife_score", "white_score", "noise_factor"):
        if name not in columns:
            connection.execute(f"ALTER TABLE detections ADD COLUMN {name} REAL")


def stable_json(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def ingest(database: Path, target_id: str, input_path: Path, candidates_path: Path,
           source_revision: str, provenance_path: Path | None = None) -> dict:
    input_path = input_path.resolve()
    candidates_path = candidates_path.resolve()
    if not input_path.is_file() or not candidates_path.is_file():
        raise ValueError("input light curve and candidates CSV must exist")
    input_hash = sha256(input_path)
    candidates_hash = sha256(candidates_path)
    provenance = {
        "input_path": str(input_path), "input_sha256": input_hash,
        "candidates_path": str(candidates_path), "candidates_sha256": candidates_hash,
        "source_revision": source_revision,
    }
    if provenance_path:
        provenance_path = provenance_path.resolve()
        if not provenance_path.is_file():
            raise ValueError("provenance file does not exist")
        provenance["manifest_path"] = str(provenance_path)
        provenance["manifest_sha256"] = sha256(provenance_path)
        try:
            with provenance_path.open(encoding="utf-8") as stream:
                provenance["manifest"] = json.load(stream)
        except (UnicodeDecodeError, json.JSONDecodeError) as error:
            raise ValueError(f"provenance file must contain valid UTF-8 JSON: {error}") from error
    run_id = hashlib.sha256(stable_json({"target_id": target_id,
                                         "source_revision": source_revision,
                                         "input_sha256": input_hash,
                                         "candidates_sha256": candidates_hash}).encode()).hexdigest()
    created = datetime.now(timezone.utc).isoformat(timespec="seconds")
    database.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(database) as connection:
        initialize(connection)
        if connection.execute("SELECT 1 FROM runs WHERE run_id=?", (run_id,)).fetchone():
            return {"run_id": run_id, "inserted": False, "detection_count": 0,
                    "reason": "identical run already registered"}
        connection.execute("INSERT INTO runs VALUES (?,?,?,?,?,?,?)",
                           (run_id, target_id, created, source_revision, input_hash,
                            candidates_hash, stable_json(provenance)))
        inserted = 0
        with candidates_path.open(newline="", encoding="utf-8-sig") as stream:
            for row in csv.DictReader(stream):
                rank = optional_number(row.get("rank"))
                candidate_key = stable_json({"run_id": run_id, "rank": rank,
                                             "candidate": row})
                detection_id = hashlib.sha256(candidate_key.encode()).hexdigest()
                values = [optional_number(row.get(field)) for field in NUMERIC_FIELDS]
                values[0] = int(values[0]) if values[0] is not None else None
                values[10] = int(values[10]) if values[10] is not None else None
                connection.execute("""INSERT INTO detections (
                    detection_id,run_id,candidate_rank,candidate_class,class_rank,score,
                    jackknife_score,white_score,noise_factor,
                    period_days,epoch,duration_days,depth_ppm,amplitude_ppm,n_events,
                    odd_even_sigma,tier,event_evidence,alias_ephemerides,note,
                    raw_candidate_json,status,created_utc
                    ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                    (detection_id, run_id, int(rank) if rank is not None else None,
                     row.get("class", "unknown"), *values, row.get("tier"),
                     row.get("event_evidence"), row.get("alias_ephemerides"),
                     row.get("note"), stable_json(row), "DETECTED_SIGNAL", created))
                connection.execute("INSERT INTO status_history VALUES (?,?,?,?,?,?)",
                                   (detection_id, None, "DETECTED_SIGNAL", created,
                                    "system", "Imported from hashed candidate output"))
                inserted += 1
    return {"run_id": run_id, "inserted": True, "detection_count": inserted,
            "input_sha256": input_hash, "candidates_sha256": candidates_hash}


def transition(database: Path, detection_id: str, new_status: str,
               actor: str, evidence: str) -> dict:
    if new_status not in STATUSES:
        raise ValueError(f"status must be one of {', '.join(STATUSES)}")
    if not actor.strip() or not evidence.strip():
        raise ValueError("a human actor and evidence reference are required")
    with sqlite3.connect(database) as connection:
        initialize(connection)
        row = connection.execute("SELECT status FROM detections WHERE detection_id=?",
                                 (detection_id,)).fetchone()
        if row is None:
            raise ValueError("unknown detection_id")
        previous = row[0]
        if STATUSES.index(new_status) != STATUSES.index(previous) + 1:
            raise ValueError(f"promotion must advance one stage from {previous}")
        changed = datetime.now(timezone.utc).isoformat(timespec="seconds")
        connection.execute("UPDATE detections SET status=? WHERE detection_id=?",
                           (new_status, detection_id))
        connection.execute("INSERT INTO status_history VALUES (?,?,?,?,?,?)",
                           (detection_id, previous, new_status, changed,
                            actor.strip(), evidence.strip()))
    return {"detection_id": detection_id, "previous_status": previous,
            "status": new_status, "actor": actor.strip(), "evidence": evidence.strip()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    add = sub.add_parser("ingest", help="append one hashed search run")
    add.add_argument("--database", type=Path, required=True)
    add.add_argument("--target-id", required=True)
    add.add_argument("--input", type=Path, required=True)
    add.add_argument("--candidates", type=Path, required=True)
    add.add_argument("--source-revision", required=True)
    add.add_argument("--provenance", type=Path, required=True,
                     help="JSON manifest with observing products and data provenance")
    promote = sub.add_parser("promote", help="record a human-evidenced stage transition")
    promote.add_argument("--database", type=Path, required=True)
    promote.add_argument("--detection-id", required=True)
    promote.add_argument("--status", choices=STATUSES[1:], required=True)
    promote.add_argument("--actor", required=True)
    promote.add_argument("--evidence", required=True)
    listing = sub.add_parser("list", help="list registered detections")
    listing.add_argument("--database", type=Path, required=True)
    listing.add_argument("--target-id")
    args = parser.parse_args()
    try:
        if args.command == "ingest":
            result = ingest(args.database, args.target_id, args.input,
                            args.candidates, args.source_revision, args.provenance)
        elif args.command == "promote":
            result = transition(args.database, args.detection_id, args.status,
                                args.actor, args.evidence)
        else:
            with sqlite3.connect(args.database) as connection:
                initialize(connection)
                query = "SELECT detection_id,run_id,candidate_rank,candidate_class,score,jackknife_score,period_days,epoch,duration_days,depth_ppm,status FROM detections"
                values = ()
                if args.target_id:
                    query += " WHERE run_id IN (SELECT run_id FROM runs WHERE target_id=?)"
                    values = (args.target_id,)
                query += " ORDER BY created_utc,detection_id"
                columns = [column[0] for column in connection.execute(query, values).description]
                result = [dict(zip(columns, row)) for row in connection.execute(query, values)]
    except (OSError, ValueError, sqlite3.Error, csv.Error) as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
