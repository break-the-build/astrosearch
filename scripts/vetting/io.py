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
    #: Rows the reader dropped or merged (None for in-memory curves).
    stats: "TableStats | None" = None

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


FLUX_COLUMNS = ("flux", "pdcsap_flux", "sap_flux", "relative_flux")
ERROR_COLUMNS = ("error", "flux_err", "flux_error", "uncertainty")
PAIRED_ERRORS = {"pdcsap_flux": "pdcsap_flux_err", "sap_flux": "sap_flux_err"}
QUALITY_COLUMNS = ("quality", "quality_flag")


@dataclass(frozen=True)
class TableStats:
    """What the reader kept and dropped, mirroring the Rust search's warnings."""
    rows_read: int
    rows_skipped: int
    duplicates_merged: int


def _number(value) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _header(lines: list[str]) -> tuple[int, str]:
    for index, line in enumerate(lines):
        text = line.lstrip("﻿").strip()
        if text and not text.startswith("#"):
            return index, text
    raise ValueError("empty input")


def read_table_with_stats(path: Path) -> tuple[np.ndarray, TableStats]:
    """Parse a light-curve CSV/TSV with the same rules as astrosearch
    (src/input.rs): the first column alias found in preference order; rows
    skipped when quality is not exactly 0, or time/flux are non-finite or flux
    is non-positive; non-finite or non-positive errors become NaN; exact
    duplicate times merged to the median flux and median positive error.
    Times are returned in BTJD."""
    lines = Path(path).read_text(encoding="utf-8-sig").splitlines()
    start, header = _header(lines)
    delimiter = "\t" if "\t" in header else ","
    body = [line for line in lines[start + 1:] if line.strip() and not line.lstrip().startswith("#")]
    reader = csv.reader([header, *body], delimiter=delimiter, skipinitialspace=True)
    names = [name.strip().lower() for name in next(reader)]

    def column(keys):
        return next((names.index(k) for k in keys if k in names), None)

    time_col = column(TIME_COLUMNS)
    flux_col = column(FLUX_COLUMNS)
    if time_col is None:
        raise ValueError(f"{path}: header needs a time/BTJD/BJD column {TIME_COLUMNS}")
    if flux_col is None:
        raise ValueError(f"{path}: header needs a flux column {FLUX_COLUMNS}")
    paired = PAIRED_ERRORS.get(names[flux_col])
    error_col = column((paired,)) if paired else None
    error_col = error_col if error_col is not None else column(ERROR_COLUMNS)
    quality_col = column(QUALITY_COLUMNS)
    rows, skipped, read = [], 0, 0
    for record in reader:
        read += 1
        cell = lambda i: _number(record[i].strip()) if i is not None and i < len(record) else None  # noqa: E731
        if quality_col is not None and cell(quality_col) != 0.0:
            skipped += 1
            continue
        t, f = cell(time_col), cell(flux_col)
        if t is None or f is None or not (np.isfinite(t) and np.isfinite(f) and f > 0):
            skipped += 1
            continue
        e = cell(error_col)
        rows.append((t, f, e if e is not None and np.isfinite(e) and e > 0 else np.nan))
    data = np.array(rows, dtype=float).reshape(-1, 3)
    data = data[np.argsort(data[:, 0], kind="stable")]
    merged, duplicates = _merge_duplicates(data)
    merged[:, 0], _ = normalize_time(merged[:, 0])
    return merged, TableStats(read, skipped, duplicates)


def _merge_duplicates(data: np.ndarray) -> tuple[np.ndarray, int]:
    if len(data) < 2 or np.all(np.diff(data[:, 0]) > 0):
        return data, 0
    times, starts = np.unique(data[:, 0], return_index=True)
    out = np.empty((len(times), 3))
    bounds = [*starts[1:], len(data)]
    for i, (lo, hi) in enumerate(zip(starts, bounds)):
        group = data[lo:hi]
        errors = group[:, 2][np.isfinite(group[:, 2])]
        out[i] = (times[i], np.median(group[:, 1]), np.median(errors) if errors.size else np.nan)
    return out, len(data) - len(times)


def read_table(path: Path) -> np.ndarray:
    return read_table_with_stats(path)[0]


def read_lightcurve(path: Path) -> LightCurve:
    """A combined CSV with `combined.meta.json` maps samples back to parts by
    time range; a single-sector CSV uses its own `<name>.meta.json`."""
    data, stats = read_table_with_stats(path)
    if len(data) < 30:
        raise ValueError(f"{path}: need at least 30 valid observations, found {len(data)} "
                         f"({stats.rows_skipped} of {stats.rows_read} rows skipped)")
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
    return LightCurve(data[:, 0], data[:, 1], error, part, tuple(parts), stats)


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
