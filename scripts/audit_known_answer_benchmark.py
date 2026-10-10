#!/usr/bin/env python3
"""Audit the frozen inputs, opaque run manifest, and outputs of a known-answer release.

This checks exact required-case coverage, input and output hashes, and uniform
search configuration before benchmark metrics are interpreted. Calibration
positives are deliberately excluded from required outputs; calibration controls
and every blind case are required.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def audit(run_manifest: Path, truth_manifest: Path, provenance_path: Path) -> dict:
    root = run_manifest.resolve().parent
    run_rows = read_csv(run_manifest)
    truth_rows = read_csv(truth_manifest)
    truth = {row["case_id"]: row for row in truth_rows}
    if len(truth) != len(truth_rows):
        raise ValueError("truth manifest contains duplicate case IDs")

    expected = {
        row["case_id"] for row in truth_rows
        if (row["split"] == "calibration"
            and row["truth_class"] in {"inverted_null", "null"})
        or row["split"] == "blind"
    }
    run_ids = [row["case_id"] for row in run_rows]
    if len(run_ids) != len(set(run_ids)):
        raise ValueError("run manifest contains duplicate case IDs")
    run_set = set(run_ids)
    if run_set != expected:
        missing = sorted(expected - run_set)
        unexpected = sorted(run_set - expected)
        raise ValueError(f"run/truth coverage mismatch: missing={len(missing)}, unexpected={len(unexpected)}")
    if any(key in run_rows[0] for key in ("truth_class", "host_id", "planet_name",
                                           "period_days", "epoch_days")):
        raise ValueError("label-bearing fields leaked into the run manifest")

    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    split_counts = Counter((row["split"], row["truth_class"]) for row in truth_rows)
    missing_products = []
    input_hash_mismatches = []
    candidate_hash_mismatches = []
    configs: set[str] = set()
    binary_hashes: set[str] = set()
    source_hashes: set[str] = set()
    telemetry_count = 0
    schemas = Counter()
    capped = []
    uncapped = 0
    requested_trials = []

    for row in run_rows:
        case_id = row["case_id"]
        data = (root / row["data_csv"]).resolve()
        candidate = (root / row["candidates_csv"]).resolve()
        sidecar = candidate.with_suffix(".run.json")
        log = candidate.with_suffix(".log")
        if not all(path.is_file() for path in (data, candidate, sidecar, log)):
            missing_products.append(case_id)
            continue

        state_file = json.loads(sidecar.read_text(encoding="utf-8"))
        state = state_file.get("state", {})
        if sha256(data) != state.get("input_sha256") or state.get("input_sha256") != truth[case_id].get("data_sha256"):
            input_hash_mismatches.append(case_id)
        if sha256(candidate) != state_file.get("candidate_sha256"):
            candidate_hash_mismatches.append(case_id)
        configs.add(json.dumps(state.get("config", {}), sort_keys=True))
        binary_hashes.add(state.get("binary_sha256", ""))
        build = state.get("build_identity", {})
        if build.get("source_tree_sha256"):
            source_hashes.add(build["source_tree_sha256"])
        telemetry_count += bool(state_file.get("metrics"))
        with candidate.open(newline="", encoding="utf-8-sig") as stream:
            schemas[tuple(next(csv.reader(stream), []))] += 1

        match = re.search(r"needs (\d+) trials.*capped at (\d+)",
                          log.read_text(encoding="utf-8", errors="replace"))
        if match:
            capped.append((int(match.group(1)), int(match.group(2))))
            requested_trials.append(int(match.group(1)))
        else:
            uncapped += 1

    if missing_products or input_hash_mismatches or candidate_hash_mismatches:
        raise ValueError(
            f"output integrity failed: missing={len(missing_products)}, "
            f"input_hash_mismatches={len(input_hash_mismatches)}, "
            f"candidate_hash_mismatches={len(candidate_hash_mismatches)}"
        )
    if len(configs) != 1 or len(binary_hashes) != 1 or "" in binary_hashes:
        raise ValueError("results do not share one recorded search configuration and binary")

    config = json.loads(next(iter(configs)))
    expected_search = provenance.get("benchmark_search", {})
    context_matches = (
        config.get("threads") == expected_search.get("threads_per_process")
        and config.get("transit_search_passes") == expected_search.get("transit_search_passes")
        and config.get("grouping") == expected_search.get("grouping")
        and config.get("min_period") == expected_search.get("min_period_days")
        and config.get("max_period") == expected_search.get("max_period_days")
    )

    sidecars = list((root / "results").glob("*/candidates.run.json"))
    demanded = sorted(requested_trials)
    return {
        "release": root.name,
        "run_case_count": len(run_rows),
        "truth_case_count": len(truth_rows),
        "truth_counts": {f"{split}:{truth_class}": count
                         for (split, truth_class), count in sorted(split_counts.items())},
        "blind_positive_host_count": len({row["host_id"] for row in truth_rows
                                           if row["split"] == "blind"
                                           and row["truth_class"] == "confirmed_planet"}),
        "calibration_positives_not_searched_by_design": sum(
            row["split"] == "calibration" and row["truth_class"] == "confirmed_planet"
            for row in truth_rows),
        "run_manifest_has_no_truth_labels": True,
        "required_outputs_complete": True,
        "input_hashes_match_truth_and_sidecars": True,
        "candidate_hashes_match_sidecars": True,
        "unique_binary_sha256": next(iter(binary_hashes)),
        "source_tree_sha256_values": sorted(source_hashes),
        "search_config": config,
        "search_config_matches_provenance": context_matches,
        "candidate_output_schemas": {";".join(schema): count for schema, count in schemas.items()},
        "period_grid": {
            "capped_cases": len(capped),
            "uncapped_cases": uncapped,
            "cap_fraction": len(capped) / len(run_rows) if run_rows else None,
            "median_requested_trials_capped_cases": demanded[len(demanded) // 2] if demanded else None,
            "maximum_requested_trials": max(demanded) if demanded else None,
            "actual_trial_cap": capped[0][1] if capped else None,
        },
        "unreferenced_result_sidecars": max(0, len(sidecars) - len(run_rows)),
        "runtime_measurement_case_count": telemetry_count,
        "runtime_measurements_recorded": telemetry_count == len(run_rows),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-manifest", type=Path, required=True)
    parser.add_argument("--truth-manifest", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        report = audit(args.run_manifest, args.truth_manifest, args.provenance)
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        parser.error(str(error))
    encoded = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
