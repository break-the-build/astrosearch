#!/usr/bin/env python3
"""Validate and freeze an externally curated astrophysical false-positive set.

Input CSV columns: case_id, host_id, split, object_type, label_reference,
data_csv, candidates_csv. Each row must represent a documented non-planet
disposition (for example eclipsing binary, blend, variable, or artifact).
Paths are resolved relative to the labels CSV. This tool does not infer labels.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

REQUIRED = ("case_id", "host_id", "split", "object_type", "label_reference",
            "data_csv", "candidates_csv")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def freeze(labels: Path, output: Path) -> dict:
    with labels.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        missing = set(REQUIRED) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"missing required columns: {', '.join(sorted(missing))}")
        rows = list(reader)
    if not rows:
        raise ValueError("labels file contains no cases")
    seen_cases, host_splits, frozen = set(), {}, []
    for number, row in enumerate(rows, 2):
        case = row["case_id"].strip()
        host = row["host_id"].strip()
        split = row["split"].strip().lower()
        object_type = row["object_type"].strip()
        reference = row["label_reference"].strip()
        if not case or not host or not object_type or not reference:
            raise ValueError(f"line {number}: case, host, object_type, and label_reference are required")
        if split not in {"calibration", "blind"}:
            raise ValueError(f"line {number}: split must be calibration or blind")
        if case in seen_cases:
            raise ValueError(f"line {number}: duplicate case_id {case}")
        if host in host_splits and host_splits[host] != split:
            raise ValueError(f"line {number}: host {host} appears in both splits")
        seen_cases.add(case)
        host_splits[host] = split
        data = (labels.parent / row["data_csv"].strip()).resolve()
        candidates = (labels.parent / row["candidates_csv"].strip()).resolve()
        if not data.is_file() or not candidates.is_file():
            raise ValueError(f"line {number}: missing curve or candidate output")
        frozen.append({
            **row, "case_id": case, "host_id": host, "split": split,
            "truth_class": "known_false_positive", "false_positive_type": object_type,
            "label_reference": reference,
            "data_csv": str(data), "candidates_csv": str(candidates),
            "data_sha256": sha256(data), "candidates_sha256": sha256(candidates),
        })
    blind_types = sorted({r["false_positive_type"] for r in frozen if r["split"] == "blind"})
    if not blind_types:
        raise ValueError("need at least one blind known-false-positive case")
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in frozen for key in row))
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(frozen)
    return {
        "case_count": len(frozen),
        "calibration_count": sum(r["split"] == "calibration" for r in frozen),
        "blind_count": sum(r["split"] == "blind" for r in frozen),
        "blind_false_positive_types": blind_types,
        "labels_sha256": sha256(labels), "manifest_sha256": sha256(output),
        "label_source": "externally supplied label_reference; no label inference",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        print(json.dumps(freeze(args.labels.expanduser().resolve(),
                                args.output.expanduser().resolve()), indent=2))
    except (OSError, ValueError, csv.Error) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
