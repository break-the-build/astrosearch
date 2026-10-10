"""Read astrosearch light curves, their provenance sidecars, and candidate CSVs."""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
csv.field_size_limit(2**31 - 1)  # candidate rows carry long evidence/alias cells

TRANSIT_CLASSES = ("periodic_dip", "duo_transit")


@dataclass(frozen=True)
class LightCurve:
    time: np.ndarray
    flux: np.ndarray
    error: np.ndarray
    #: Index into `parts` for each sample (which sector/product it came from).
    part: np.ndarray
    #: Per-part metadata from `*.meta.json` (may be empty dicts if absent).
    parts: tuple = field(default_factory=tuple)

    def sector_of(self, index: int) -> str:
        meta = self.parts[self.part[index]] if self.parts else {}
        return str(meta.get("sector", self.part[index]))


@dataclass(frozen=True)
class Candidate:
    rank: int
    kind: str
    score: float
    period: float
    epoch: float
    duration: float
    depth_ppm: float
    n_events: int
    event_epochs: tuple
    row: dict


#: BTJD = BJD_TDB - 2457000. Times above this are taken as full BJD.
BTJD_OFFSET = 2457000.0
FULL_BJD_MINIMUM = 2_400_000.0
TIME_COLUMNS = ("time", "btjd", "bjd", "t")


def normalize_time(time: np.ndarray) -> tuple[np.ndarray, float]:
    """Express times in BTJD; returns the times and the offset subtracted.
    astrosearch accepts `bjd` columns, but the gate's windows, catalogs and
    CTOI row assume BTJD, so full BJD inputs are shifted here once."""
    if time.size and np.nanmedian(time) > FULL_BJD_MINIMUM:
        return time - BTJD_OFFSET, BTJD_OFFSET
    return time, 0.0


def normalize_epoch(epoch: float) -> float:
    return epoch - BTJD_OFFSET if epoch > FULL_BJD_MINIMUM else epoch


def read_table(path: Path) -> np.ndarray:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        fields = {name.strip().lower(): name for name in reader.fieldnames or []}
        key = next((fields[c] for c in TIME_COLUMNS if c in fields), None)
        if key is None:
            raise ValueError(f"{path}: needs one of the time columns {TIME_COLUMNS}")
        rows = [(float(r[key]), float(r["flux"]), float(r.get("flux_err") or "nan"))
                for r in reader]
    data = np.array(rows, dtype=float).reshape(-1, 3)
    data[:, 0], _ = normalize_time(data[:, 0])
    return data


def read_lightcurve(path: Path) -> LightCurve:
    """A combined CSV with `combined.meta.json` maps samples back to parts by
    time range; a single-sector CSV uses its own `<name>.meta.json`."""
    data = read_table(path)
    order = np.argsort(data[:, 0], kind="stable")
    data = data[order]
    parts: list[dict] = []
    combined = path.with_suffix(".meta.json")
    if combined.exists():
        meta = json.loads(combined.read_text())
        parts = list(meta.get("parts", [meta]))
    part = np.zeros(len(data), dtype=int)
    for index, meta in enumerate(parts):
        lo, hi = meta.get("time_min"), meta.get("time_max")
        if lo is not None and hi is not None:
            part[(data[:, 0] >= lo - 1e-9) & (data[:, 0] <= hi + 1e-9)] = index
    error = data[:, 2]
    if not np.all(np.isfinite(error) & (error > 0)):
        fallback = np.nanmedian(np.abs(np.diff(data[:, 1]))) / np.sqrt(2) * 1.4826
        error = np.where(np.isfinite(error) & (error > 0), error, fallback)
    return LightCurve(data[:, 0], data[:, 1], error, part, tuple(parts))


def parse_event_epochs(text: str) -> tuple:
    """`event_evidence` cells encode `midpoint:depth:snr:points;...`."""
    epochs = []
    for item in (text or "").split(";"):
        head = item.split(":", 1)[0].strip()
        if head:
            try:
                epochs.append(float(head))
            except ValueError:
                continue
    return tuple(epochs)


def read_candidates(path: Path, classes=TRANSIT_CLASSES) -> list[Candidate]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = [r for r in csv.DictReader(stream) if r.get("class") in classes]
    return [Candidate(
        rank=int(r["rank"]), kind=r["class"], score=float(r["score"]),
        period=float(r["period_days"]), epoch=normalize_epoch(float(r["epoch"])),
        duration=float(r["duration_days"]), depth_ppm=float(r["depth_ppm"]),
        n_events=int(r["n_events"]),
        event_epochs=tuple(normalize_epoch(t) for t in parse_event_epochs(r.get("event_evidence", ""))),
        row=r) for r in rows]


def phase_offset(time: np.ndarray, epoch: float, period: float) -> np.ndarray:
    """Signed time from the nearest predicted mid-transit."""
    return (time - epoch + 0.5 * period) % period - 0.5 * period


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=1, sort_keys=True, default=_jsonable) + "\n")


def _jsonable(value):
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"not JSON serializable: {type(value)}")
