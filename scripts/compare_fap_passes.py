#!/usr/bin/env python3
"""Do single-pass null searches understate the FAP of a multi-pass search?

The gate's empirical FAP searched its nulls with one pass, while candidates
come from a five-pass search. The null statistic is the best transit score
over everything the search reports; extra passes can only add candidates, so
the multi-pass null maximum is >= the single-pass one. This script measures
how often and by how much, on block-shuffled and inverted nulls of one target.
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from vetting import fap  # noqa: E402
from vetting.io import read_lightcurve  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lightcurve", type=Path, required=True)
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--epoch", type=float, required=True)
    p.add_argument("--period", type=float, required=True)
    p.add_argument("--duration", type=float, required=True)
    p.add_argument("--trials", type=int, default=16)
    p.add_argument("--passes", type=int, default=5)
    p.add_argument("--workers", type=int, default=6)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--work-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    lc = read_lightcurve(a.lightcurve)
    rng = np.random.default_rng(20261010)
    a.work_dir.mkdir(parents=True, exist_ok=True)
    curves = [a.work_dir / "inverted.csv"]
    fap.write_curve(curves[0], *fap.inverted(lc))
    for i in range(a.trials - 1):
        curves.append(a.work_dir / f"shuffled_{i:03d}.csv")
        fap.write_curve(curves[-1], *fap.block_shuffled(lc, rng, a.epoch, a.period, a.duration))
    jobs = [(c, n) for c in curves for n in (1, a.passes)]

    def score(job):
        curve, passes = job
        return fap.best_transit_score(a.binary, curve, ["--transit-search-passes", str(passes)], a.threads)

    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        scores = list(pool.map(score, jobs))
    single = np.array(scores[0::2])
    multi = np.array(scores[1::2])
    out = {"trials": len(curves), "passes": a.passes,
           "single_pass_max": single.tolist(), "multi_pass_max": multi.tolist(),
           "multi_exceeds_single": int((multi > single + 1e-9).sum()),
           "median_ratio": float(np.median(multi / np.maximum(single, 1e-9))),
           "max_ratio": float(np.max(multi / np.maximum(single, 1e-9)))}
    a.output.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: v for k, v in out.items() if not k.endswith("_max")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
