"""Separate grid-resolution losses from ranking losses for a known periodic dip.

usage: diagnose_missed_signal.py BINARY INPUT PERIOD OUT_DIR [--fraction F]

For each scoring statistic (white, correlated):
  full oracle    dense period (20,000 trials) and duration (64) grid within
                 PERIOD * (1 +/- F): the best score the statistic can give the signal
  period oracle  the same dense periods with the default 10 durations
  grid best      best trial within one duration of the signal in the default full
                 search (the score the signal actually receives)
  rank           position among deduplicated periodic dips with no 12-candidate
                 cap, and the scores of stronger peaks

duration-grid loss = full oracle - period oracle
period-grid loss   = period oracle - grid best
ranking loss       = how many distinct peaks outrank the grid best
"""
import argparse
import csv
import re
import subprocess
from pathlib import Path


def run(binary: str, args: list[str], log: Path) -> float:
    result = subprocess.run([binary, *args], capture_output=True, text=True, check=True)
    log.write_text(result.stderr)
    return float(re.search(r"across ([0-9.]+) days", result.stderr).group(1))


def dips(path: Path) -> list[dict]:
    with path.open(newline="") as stream:
        rows = [r for r in csv.DictReader(stream) if r["class"] in ("periodic_dip", "duo_transit")]
    for row in rows:
        for key in ("score", "white_score", "period_days", "epoch", "duration_days", "depth_ppm"):
            row[key] = float(row[key])
    return rows


def near(row: dict, period: float, epoch: float, duration: float, span: float) -> bool:
    drift = abs(row["period_days"] - period) * span / period
    phase = abs((row["epoch"] - epoch + period / 2) % period - period / 2)
    return drift <= duration and phase <= duration


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("binary")
    parser.add_argument("input")
    parser.add_argument("period", type=float)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--fraction", type=float, default=0.002)
    a = parser.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    lo, hi = a.period * (1 - a.fraction), a.period * (1 + a.fraction)
    common = ["--input", a.input, "--min-transit-score", "3", "--transit-candidates", "100000"]
    print("| scoring | full oracle | period oracle | grid best | duration-grid loss | "
          "period-grid loss | uncapped rank | stronger peaks (top 3: P d / score / events) |")
    print("| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |")
    for scoring in ("white", "correlated"):
        tag = a.out_dir / scoring
        oracle = {}
        for name, durations in (("full_oracle", "64"), ("period_oracle", "10")):
            out = Path(f"{tag}_{name}.csv")
            run(a.binary, [*common, "--scoring", scoring, "--min-period", str(lo),
                           "--max-period", str(hi), "--periods", "20000",
                           "--durations", durations, "--output", str(out)],
                Path(f"{tag}_{name}.log"))
            oracle[name] = max(dips(out), key=lambda r: r["score"])
        truth = oracle["full_oracle"]
        period, epoch, duration = truth["period_days"], truth["epoch"], truth["duration_days"]
        out, trials = Path(f"{tag}_full.csv"), Path(f"{tag}_periodogram.csv")
        span = run(a.binary, [*common, "--scoring", scoring, "--output", str(out),
                              "--periodogram", str(trials)], Path(f"{tag}_full.log"))
        close = [r for r in dips(trials) if near(r, period, epoch, duration, span)]
        grid = max((r["score"] for r in close), default=float("nan"))
        ranked = dips(out)
        rank = next((i + 1 for i, r in enumerate(ranked)
                     if near(r, period, epoch, duration, span)), None)
        stronger = ranked[: (rank - 1) if rank else len(ranked)]
        top = "; ".join(f"{r['period_days']:.3f} / {r['score']:.1f} / {r['n_events']}"
                        for r in stronger[:3])
        full, per = truth["score"], oracle["period_oracle"]["score"]
        print(f"| {scoring} | {full:.2f} | {per:.2f} | {grid:.2f} | {full - per:.2f} | "
              f"{per - grid:.2f} | {rank or 'absent'} ({len(stronger)} stronger) | {top} |")
        print(f"<!-- {scoring}: oracle P={period:.5f} epoch={epoch:.4f} "
              f"duration={duration:.4f} depth={truth['depth_ppm']:.0f} "
              f"white_score={truth['white_score']:.2f}; {len(close)} grid trials near -->")


if __name__ == "__main__":
    main()
