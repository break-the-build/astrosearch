"""Detection-tier statistics computed from astrosearch candidate CSV rows.

Mirrors `Candidate::jackknife_score` in the Rust crate so the calibration and
the gate can evaluate CSVs written before the column existed.
"""
from __future__ import annotations

import math

STATISTICS = ("score", "jackknife_score")


def event_snrs(row: dict) -> list[float]:
    """Per-event S/N from `event_evidence` (`midpoint:depth:snr:points;...`)."""
    out = []
    for item in (row.get("event_evidence") or "").split(";"):
        fields = item.split(":")
        if len(fields) >= 3:
            try:
                out.append(float(fields[2]))
            except ValueError:
                continue
    return out


def single_event_share(row: dict) -> float:
    squares = [max(s, 0.0) ** 2 for s in event_snrs(row)]
    total = sum(squares)
    return max(squares) / total if total > 0 else 1.0


def jackknife_score(row: dict) -> float:
    if row.get("jackknife_score") not in (None, ""):
        return float(row["jackknife_score"])
    return float(row["score"]) * math.sqrt(max(0.0, 1.0 - single_event_share(row)))


def statistic(row: dict, name: str) -> float:
    if name == "score":
        return float(row["score"])
    if name == "jackknife_score":
        return jackknife_score(row)
    raise ValueError(f"unknown detection statistic {name!r}; expected one of {STATISTICS}")
