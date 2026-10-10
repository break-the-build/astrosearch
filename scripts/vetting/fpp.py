"""TRICERATOPS false-positive probability (Giacalone et al. 2021).

Optional: requires `triceratops` (pyproject extra `fpp`) and network access to
MAST and Gaia. Uses SPOC apertures; sectors without one are dropped and named
in the report. FPP < 0.015 and NFPP < 0.001 are the published thresholds for
statistical validation, which also requires high-resolution imaging.
"""
from __future__ import annotations

import numpy as np

from .crosscheck import flatten
from .io import LightCurve, phase_offset


def triceratops_fpp(tic_id: int, lc: LightCurve, sectors: list[int], period: float,
                    epoch: float, duration: float, depth_ppm: float,
                    draws: int = 200_000, window_durations: float = 3.0) -> dict:
    try:
        import triceratops.triceratops as tr
    except ImportError:
        return {"status": "NOT_EVALUATED",
                "reason": "triceratops not installed (uv sync --extra fpp)"}
    try:
        probe = tr.target(ID=int(tic_id), sectors=np.array(sectors))
        apertures = probe.get_spoc_apertures()
        usable = [s for s, ap in zip(sectors, apertures) if ap is not None]
        if not usable:
            return {"status": "NOT_EVALUATED", "reason": "no SPOC apertures for these sectors",
                    "sectors": sectors}
        target = probe if len(usable) == len(sectors) else tr.target(ID=int(tic_id),
                                                                     sectors=np.array(usable))
        aps = [ap for ap in apertures if ap is not None]
        target.calc_depths(tdepth=depth_ppm * 1e-6, all_ap_pixels=aps)
        flux = flatten(lc, epoch, period, duration)
        offset = phase_offset(lc.time, epoch, period)
        near = np.abs(offset) <= window_durations * duration
        target.calc_probs(time=offset[near], flux_0=flux[near],
                          flux_err_0=float(np.median(lc.error[near])), P_orb=period,
                          N=draws, verbose=0)
        fpp, nfpp = float(target.FPP), float(target.NFPP)
    except Exception as error:  # network, catalog, or model failures
        return {"status": "NOT_EVALUATED", "reason": f"{type(error).__name__}: {error}"}
    return {"status": "PASS" if fpp < 0.01 else "FAIL", "fpp": fpp, "nfpp": nfpp,
            "sectors_used": usable, "sectors_dropped": [s for s in sectors if s not in usable],
            "draws": draws,
            "validation_thresholds": "FPP < 0.015 and NFPP < 0.001 plus high-resolution imaging"}
