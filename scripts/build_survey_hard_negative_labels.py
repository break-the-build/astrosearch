#!/usr/bin/env python3
"""Create hard-negative labels from completed TESS EB survey manifests."""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

REFERENCE = (
    "Vetted eclipsing binary in the TESS Ten Thousand Catalog: "
    "Kostov et al. 2025, ApJS 279, 50, "
    "https://doi.org/10.3847/1538-4365/ade2d8; "
    "VizieR J/ApJS/279/50/table3, https://doi.org/10.26093/cds/vizier/bb-8q-1j"
)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--survey-output", type=Path, required=True)
    parser.add_argument("--selection", type=Path,
                        default=root / "data/benchmarks/hard-negatives/selected-eb-targets.csv")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    survey = args.survey_output.expanduser().resolve()
    selection = args.selection.expanduser().resolve()
    output = args.output.expanduser().resolve()
    with selection.open(newline="", encoding="utf-8-sig") as stream:
        selected = {row["tic_id"]: row for row in csv.DictReader(stream)}
    fields = ("case_id", "host_id", "split", "object_type", "label_reference",
              "data_csv", "candidates_csv")
    rows = []
    for manifest_path in sorted(survey.rglob("TIC*/manifest.json")):
        try:
            manifest = json.loads(manifest_path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        tic = str(manifest.get("tic_id", ""))
        if tic not in selected:
            continue
        preferred = selected[tic].get("preferred_author", "").strip()
        usable_authors = [preferred] if preferred else ["QLP", "TESS-SPOC", "SPOC"]
        result = next((item for author in usable_authors
                       for item in manifest.get("authors", [])
                       if item.get("author") == author and item.get("status") == "searched"), None)
        if not result:
            continue
        data = Path(result.get("combined_csv", "")).resolve()
        candidates = Path(result.get("candidates_csv", "")).resolve()
        if not data.is_file() or not candidates.is_file():
            continue
        rows.append({
            "case_id": f"eb-tic-{tic}", "host_id": tic, "split": "blind",
            "object_type": "vetted_eclipsing_binary", "label_reference": REFERENCE,
            "data_csv": os.path.relpath(data, output.parent),
            "candidates_csv": os.path.relpath(candidates, output.parent),
        })
    if not rows:
        parser.error("no successfully searched selected EB targets found")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"selected_targets": len(selected), "searched_ebs": len(rows),
                      "blind_hard_negatives": len(rows), "labels_csv": str(output)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
