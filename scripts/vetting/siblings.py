"""Mask other detected transit signals ("siblings") before vetting a candidate.

In a multi-planet system the other planets' transits fall at arbitrary phases
of the candidate's period, so an unmasked light curve fails the secondary
eclipse scan and the odd/even test for every planet but the dominant one
(L 98-59 d in the known-systems benchmark). Siblings are masked only when they
are significant and not commensurate with the candidate: a same-period dip at
another phase is exactly the eclipsing-binary secondary the scan must see.
"""
from __future__ import annotations

import numpy as np

from .io import Candidate, LightCurve, phase_offset

#: Period ratios whose events could be the candidate's own secondary or alias.
COMMENSURATE_RATIOS = (1 / 3, 0.5, 1.0, 2.0, 3.0)


def commensurate(target: Candidate, other: Candidate, baseline: float) -> bool:
    """True when `other`'s period is r * target period closely enough that its
    events stay within a duration of a fixed target phase over the baseline."""
    tolerance = max(target.duration, other.duration)
    for ratio in COMMENSURATE_RATIOS:
        drift = abs(other.period - ratio * target.period) * baseline / max(other.period, 1e-9)
        if drift < tolerance:
            return True
    return False


def in_transit(lc: LightCurve, cand: Candidate, factor: float = 1.0) -> np.ndarray:
    return np.abs(phase_offset(lc.time, cand.epoch, cand.period)) <= factor * cand.duration / 2


def subset(lc: LightCurve, keep: np.ndarray) -> LightCurve:
    return LightCurve(lc.time[keep], lc.flux[keep], lc.error[keep], lc.part[keep], lc.parts)


def mask_siblings(lc: LightCurve, target: Candidate, others: list[Candidate], baseline: float,
                  min_score_fraction: float = 0.5, max_siblings: int = 5,
                  mask_factor: float = 1.5, max_overlap: float = 0.2) -> tuple[LightCurve, dict]:
    """Return a light curve with significant, non-commensurate siblings'
    in-transit samples (mask_factor x duration) removed, and a record of
    what was masked or kept and why."""
    target_window = in_transit(lc, target, mask_factor)
    masked, kept = [], []
    remove = np.zeros(lc.time.size, dtype=bool)
    ranked = sorted((o for o in others if o.rank != target.rank), key=lambda o: -o.score)
    for other in ranked:
        entry = {"rank": other.rank, "period": other.period, "epoch": other.epoch,
                 "duration": other.duration, "score": other.score}
        if other.score < min_score_fraction * target.score:
            kept.append({**entry, "reason": "score below sibling threshold"})
            continue
        if len(masked) >= max_siblings:
            kept.append({**entry, "reason": "sibling limit reached"})
            continue
        if commensurate(target, other, baseline):
            kept.append({**entry, "reason": "commensurate period"})
            continue
        window = in_transit(lc, other, mask_factor)
        overlap = float((window & target_window).sum() / max(window.sum(), 1))
        if overlap > max_overlap:
            kept.append({**entry, "reason": f"overlaps candidate transits ({overlap:.0%})"})
            continue
        remove |= window
        masked.append({**entry, "samples_removed": int(window.sum())})
    record = {"masked": masked, "kept": kept, "samples_removed": int(remove.sum()),
              "min_score_fraction": min_score_fraction, "mask_factor": mask_factor}
    return (subset(lc, ~remove) if remove.any() else lc), record
