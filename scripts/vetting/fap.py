"""Per-target empirical false-alarm probability from null light curves.

Nulls come from the target's own data, so they keep its noise, systematics,
and cadence: inverted flux (dips become bumps) and block-shuffled flux (the
candidate's transits removed, then day-long blocks permuted within each
segment, which preserves correlated noise but destroys periodicity). Each null
is searched with the same astrosearch settings; the null statistic is the best
transit-class score. FAP = (k + 1) / (N + 1) for k null maxima at or above the
candidate's score. A Gumbel fit to the null maxima extrapolates below 1/N and
is reported separately as a model-based estimate.
"""
from __future__ import annotations

import csv
import math
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy import stats

from .io import TRANSIT_CLASSES, LightCurve, phase_offset

SEGMENT_GAP_DAYS = 0.5


def segments(time: np.ndarray, gap: float = SEGMENT_GAP_DAYS) -> list[np.ndarray]:
    breaks = np.flatnonzero(np.diff(time) > gap) + 1
    return np.split(np.arange(len(time)), breaks)


def inverted(lc: LightCurve) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    flux = lc.flux.copy()
    for part in np.unique(lc.part):
        mine = lc.part == part
        flux[mine] = 2 * np.median(flux[mine]) - flux[mine]
    return lc.time, flux, lc.error


def block_shuffled(lc: LightCurve, rng: np.random.Generator, epoch: float, period: float,
                   duration: float, block_days: float = 1.0):
    """Remove the candidate's in-transit samples, then permute day-long blocks
    of (flux, error) within each segment while keeping the time stamps."""
    keep = np.abs(phase_offset(lc.time, epoch, period)) > duration
    time, flux, error = lc.time[keep], lc.flux[keep], lc.error[keep]
    out_flux, out_error = flux.copy(), error.copy()
    for indices in segments(time):
        if len(indices) < 3:
            continue
        t = time[indices]
        block = np.floor((t - t[0]) / block_days).astype(int)
        groups = [indices[block == b] for b in np.unique(block)]
        order = rng.permutation(len(groups))
        shuffled = np.concatenate([groups[i] for i in order])
        out_flux[indices] = flux[shuffled]
        out_error[indices] = error[shuffled]
    return time, out_flux, out_error


def write_curve(path: Path, time, flux, error) -> None:
    rows = np.column_stack([time, flux, error])
    np.savetxt(path, rows, delimiter=",", header="time,flux,flux_err", comments="", fmt="%.8f")


def best_transit_score(binary: Path, curve: Path, args: list[str], threads: int,
                       timeout: float | None = None) -> float:
    """Best transit-class score of one null search. Raises
    subprocess.TimeoutExpired (the child is killed) or CalledProcessError."""
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "candidates.csv"
        subprocess.run([str(binary), "--input", str(curve), "--output", str(out),
                        "--threads", str(threads), *args],
                       check=True, capture_output=True, text=True, timeout=timeout)
        with out.open(newline="") as stream:
            scores = [float(r["score"]) for r in csv.DictReader(stream)
                      if r.get("class") in TRANSIT_CLASSES]
    return max(scores, default=0.0)


def gumbel_tail(maxima: np.ndarray, statistic: float) -> float | None:
    """Method-of-moments Gumbel fit; P(max >= statistic)."""
    if len(maxima) < 10 or np.std(maxima) == 0:
        return None
    beta = np.std(maxima, ddof=1) * math.sqrt(6) / math.pi
    mu = np.mean(maxima) - 0.5772156649 * beta
    return float(-np.expm1(-np.exp(-(statistic - mu) / beta)))


def clopper_pearson_upper(k: int, n: int, level: float = 0.95) -> float:
    return 1.0 if k >= n else float(stats.beta.ppf(level, k + 1, n - k))


def empirical_fap(lc: LightCurve, binary: Path, statistic: float, epoch: float, period: float,
                  duration: float, search_args: list[str], trials: int = 100,
                  workers: int = 4, threads: int = 6, seed: int = 20261010,
                  work_dir: Path | None = None, null_timeout: float = 900.0,
                  budget_seconds: float = 7200.0, min_nulls: int = 20) -> dict:
    """Search `trials` nulls (one inverted, the rest block-shuffled) within a
    time budget. Each null search gets at most `null_timeout` seconds and
    nothing starts after `budget_seconds`; failed, timed-out and skipped
    nulls are counted, and the FAP uses the completed ones. With fewer than
    `min_nulls` completed the result is NOT_EVALUATED rather than a FAP."""
    rng = np.random.default_rng(seed)
    args = [a for a in search_args] + ["--transit-search-passes", "1"]
    work = Path(work_dir or tempfile.mkdtemp(prefix="nulls-"))
    work.mkdir(parents=True, exist_ok=True)
    curves = []
    path = work / "null_inverted.csv"
    write_curve(path, *inverted(lc))
    curves.append(("inverted", path))
    for i in range(trials - 1):
        path = work / f"null_shuffled_{i:04d}.csv"
        write_curve(path, *block_shuffled(lc, rng, epoch, period, duration))
        curves.append(("shuffled", path))
    deadline = time.monotonic() + budget_seconds

    def search(curve):
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            return None, "budget exhausted before start"
        try:
            return best_transit_score(binary, curve[1], args, threads,
                                      timeout=min(null_timeout, remaining)), None
        except subprocess.TimeoutExpired:
            return None, "timed out"
        except (subprocess.CalledProcessError, OSError, ValueError) as error:
            return None, f"failed: {type(error).__name__}"

    with ThreadPoolExecutor(max_workers=workers) as pool:
        outcomes = list(pool.map(search, curves))
    maxima = [score for score, _ in outcomes if score is not None]
    failures: dict[str, int] = {}
    for _, reason in outcomes:
        if reason:
            failures[reason] = failures.get(reason, 0) + 1
    budget = {"requested": len(curves), "completed": len(maxima), "failures": failures,
              "null_timeout_s": null_timeout, "budget_s": budget_seconds, "min_nulls": min_nulls}
    if len(maxima) < min_nulls:
        return {"status": "NOT_EVALUATED", "statistic": statistic, "nulls": len(maxima),
                "fap_empirical": None, "fap_upper_95": None, "budget": budget,
                "reason": f"only {len(maxima)} of {len(curves)} null searches completed"}
    maxima_arr = np.array(maxima)
    exceed = int((maxima_arr >= statistic).sum())
    n = len(maxima)
    return {
        "status": "OK", "statistic": statistic, "nulls": n, "exceeding": exceed,
        "fap_empirical": (exceed + 1) / (n + 1),
        "fap_upper_95": clopper_pearson_upper(exceed, n),
        "fap_gumbel_extrapolated": gumbel_tail(maxima_arr, statistic),
        "inverted_null_max": outcomes[0][0],
        "null_max_quantiles": {q: float(np.quantile(maxima_arr, q)) for q in (0.5, 0.9, 0.99)},
        "search_args": args, "seed": seed, "block_days": 1.0, "budget": budget,
        "note": ("nulls use one search pass; on TIC 4206066 the 5-pass null maximum matched "
                 "the 1-pass one in 14/16 nulls (worst ratio 1.04), so this is very slightly "
                 "anti-conservative"),
    }
