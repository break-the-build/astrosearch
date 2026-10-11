"""Evidence-preserving diagnostics for transit-like detections.

Three report-only checks. None of them edits, drops or re-ranks a candidate
row; each returns annotations that cite the original rows.

1. Event groups. Detections (astrosearch rows and TLS iterations) are linked
   when they share observed events: two event midpoints match within
   max(half the summed durations, EVENT_PHASE_TOLERANCE of the shorter period),
   and a pair is linked when at least MIN_EVENT_OVERLAP of the smaller event
   set is shared. Linked rows describe one event pattern; their period ratio
   (same-period, 2:1, 3:2, ...) records the alias relationship.
2. Depth consistency. astrosearch reports either the box-fit depth or a local
   depth remeasured after transit-masked detrending (see src/refine.rs); its
   note keeps the box depth. Within one event pattern every reported depth is
   listed with the estimator that produced it, and large disagreements are
   flagged instead of resolved.
3. Variability conflict. A transit-like detection whose period is related
   (1/2, 1, 2) to a reported periodic_variability signal is tested directly:
   harmonics of the modulation are fitted per segment to out-of-transit data,
   and the in-transit depth is measured before and after removing them. A dip
   that vanishes is explained by the modulation; one that survives is not.
   The verdict is a competing explanation, never a classification.

4. Interleaved pairs. Two detections at the same or 2:1 periods that share
   few events but whose events sit at a consistent phase of the longer
   period are reported with both depths. A phase near 0.5 with unequal depths
   is the primary/secondary signature of an eclipsing binary; a multi-pass
   search often reports the secondary as a "half-period" row whose events
   all fall on one cycle parity, because the primary eclipses were masked.

With the light curve supplied, every detection's depth is also measured on
the unmasked data at its own ephemeris. astrosearch masks earlier winners
before later passes; on a smoothly varying star that raises the remaining
out-of-transit mean, so a later-pass depth can exceed anything in the
unmasked light curve. Such rows are flagged with the pass and masked
fraction astrosearch records in `note`.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

import numpy as np

EVENT_PHASE_TOLERANCE = 0.02
MIN_EVENT_OVERLAP = 0.5
DEPTH_RATIO_LIMIT = 3.0
LOCAL_BOX_WARNING = 0.25
RELATED_RATIOS = (0.5, 1.0, 2.0)
PERIOD_TOLERANCE = 0.002
HARMONICS = 2
SEGMENT_GAP_DAYS = 0.5
EXPLAINED_BELOW, SURVIVES_ABOVE = 0.3, 0.7

BOX_DEPTH = re.compile(r"\(box depth (-?[0-9.]+) ppm\)")
LOCAL_NOT_POSITIVE = re.compile(r"transit-masked local depth (-?[0-9.]+) ppm is not positive")
REMEASURED = "depth remeasured after transit-masked detrending"
PASS_NOTE = re.compile(r"found in search pass (\d+) with ([0-9.]+)% of samples masked")
UNMASKED_WARNING = 0.25
PHASE_SCATTER_LIMIT = 0.02
SECONDARY_PHASE_WINDOW = 0.05
MIN_OTHER_PHASE = 0.2  # smaller offsets are usually adjacent windows on one broad feature


@dataclass(frozen=True)
class Detection:
    source: str
    row_id: str
    period: float
    epoch: float
    duration: float
    reported_depth_ppm: float
    depth_estimator: str
    box_depth_ppm: float | None
    local_depth_ppm: float | None
    score: float | None
    odd_even_sigma: float | None
    events: tuple[float, ...] = field(repr=False)
    search_pass: int | None = None
    masked_percent: float | None = None

    def label(self) -> str:
        return f"{self.source}:{self.row_id}"


def _number(raw) -> float | None:
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def from_astrosearch(rows: list[dict]) -> list[Detection]:
    """Transit-like rows with their depth provenance parsed from `note`."""
    out = []
    for row in rows:
        if row.get("class") not in ("periodic_dip", "duo_transit"):
            continue
        note = row.get("note", "")
        reported = _number(row.get("depth_ppm")) or 0.0
        box = BOX_DEPTH.search(note)
        not_positive = LOCAL_NOT_POSITIVE.search(note)
        if REMEASURED in note:
            estimator, local = "local (transit-masked detrending)", reported
            box_depth = float(box.group(1)) if box else None  # older outputs lack it
        elif not_positive:
            estimator, local, box_depth = "box (local not positive)", float(not_positive.group(1)), reported
        else:
            estimator, local, box_depth = "box", None, reported
        events = tuple(sorted(float(e.split(":")[0]) for e in row.get("event_evidence", "").split(";") if e))
        passed = PASS_NOTE.search(note)
        out.append(Detection("astrosearch", str(row.get("rank", "")), float(row["period_days"]),
                             float(row["epoch"]), float(row["duration_days"]), reported, estimator,
                             box_depth, local, _number(row.get("score")),
                             _number(row.get("odd_even_sigma")), events,
                             int(passed.group(1)) if passed else (1 if REMEASURED in note or box else None),
                             float(passed.group(2)) if passed else None))
    return out


def observed_events(time: np.ndarray, period: float, epoch: float, duration: float) -> tuple[float, ...]:
    """Predicted midpoints with at least one sample within half a duration."""
    if time.size == 0:
        return ()
    first = math.floor((time.min() - epoch) / period)
    last = math.ceil((time.max() - epoch) / period)
    mids = epoch + period * np.arange(first, last + 1)
    sorted_time = np.sort(time)
    lo = np.searchsorted(sorted_time, mids - duration / 2, side="left")
    hi = np.searchsorted(sorted_time, mids + duration / 2, side="right")
    return tuple(float(m) for m in mids[hi > lo])


def from_tls(rows: list[dict], time: np.ndarray) -> list[Detection]:
    out = []
    for row in rows:
        period, epoch, duration = (float(row[k]) for k in ("period_days", "epoch_days", "duration_days"))
        raw = row.get("observed_transit_times", "")
        events = (tuple(sorted(float(x) for x in raw.split(";") if x)) if raw
                  else observed_events(time, period, epoch, duration))
        depth = float(row["depth_ppm"])
        out.append(Detection("tls", str(row.get("iteration", "")), period, epoch, duration, depth,
                             "tls model depth", None, None, _number(row.get("sde")),
                             _number(row.get("odd_even_mismatch_sigma")), events))
    return out


def shared_events(a: Detection, b: Detection) -> int:
    tolerance = max(0.5 * (a.duration + b.duration), EVENT_PHASE_TOLERANCE * min(a.period, b.period))
    i = j = common = 0
    while i < len(a.events) and j < len(b.events):
        difference = a.events[i] - b.events[j]
        if abs(difference) <= tolerance:
            common, i, j = common + 1, i + 1, j + 1
        elif difference < 0:
            i += 1
        else:
            j += 1
    return common


def period_relation(a: float, b: float) -> str:
    """Ratio of the longer period to the shorter, named when it is within 1%
    of an integer or half-integer."""
    ratio = max(a, b) / min(a, b)
    nearest = round(2 * ratio) / 2
    if nearest >= 1 and abs(ratio - nearest) / nearest <= 0.01:
        if nearest == 1:
            return "same-period"
        if nearest.is_integer():
            return f"{int(nearest)}:1 harmonic"
        return f"{int(2 * nearest)}:2 harmonic"
    return f"unrelated (ratio {ratio:.4f})"


def event_groups(detections: list[Detection], unmasked: dict[str, float] | None = None) -> list[dict]:
    parent = list(range(len(detections)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    relations = []
    for i in range(len(detections)):
        for j in range(i + 1, len(detections)):
            a, b = detections[i], detections[j]
            if not a.events or not b.events:
                continue
            common = shared_events(a, b)
            fraction = common / min(len(a.events), len(b.events))
            if fraction >= MIN_EVENT_OVERLAP:
                parent[find(j)] = find(i)
                relations.append({"a": a.label(), "b": b.label(), "shared_events": common,
                                  "overlap_fraction": round(fraction, 3),
                                  "period_relation": period_relation(a.period, b.period)})
    groups: dict[int, list[int]] = {}
    for i in range(len(detections)):
        groups.setdefault(find(i), []).append(i)
    out = []
    for number, members in enumerate(sorted(groups.values(), key=lambda m: min(m))):
        rows = [detections[i] for i in members]
        labels = {d.label() for d in rows}
        unmasked = unmasked or {}
        out.append({"group": number + 1,
                    "members": [{**member_summary(d), "unmasked_depth_ppm": unmasked.get(d.label())} for d in rows],
                    "relations": [r for r in relations if r["a"] in labels],
                    "depth_check": depth_check(rows, unmasked)})
    return out


def parity_coverage(d: Detection) -> str:
    """Whether observed events fall on both even and odd cycles of the
    detection's own ephemeris; one parity means an effective period of 2P."""
    if not d.events:
        return "no events"
    parities = {round((e - d.epoch) / d.period) % 2 for e in d.events}
    return "both parities" if len(parities) == 2 else (
        "even cycles only (effective period 2P)" if parities == {0} else "odd cycles only (effective period 2P)")


def _effective(d: Detection) -> float:
    return 2 * d.period if parity_coverage(d).endswith("(effective period 2P)") else d.period


def _effective_epoch(d: Detection) -> float:
    """An epoch on the cycles that hold the events: a one-parity row's own
    epoch can sit on the other parity, where it has no data."""
    return d.epoch + d.period if parity_coverage(d).startswith("odd") else d.epoch


def interleaved_pairs(detections: list[Detection]) -> list[dict]:
    out = []
    for i, a in enumerate(detections):
        for b in detections[i + 1:]:
            if len(a.events) < 3 or len(b.events) < 3:
                continue
            long_, short = (a, b) if _effective(a) >= _effective(b) else (b, a)
            reference = _effective(long_)
            if not period_relation(reference, _effective(short)).startswith(("same-period", "2:1")):
                continue
            if shared_events(a, b) >= MIN_EVENT_OVERLAP * min(len(a.events), len(b.events)):
                continue
            phases = ((np.array(short.events) - _effective_epoch(long_)) / reference) % 1.0
            angle = np.angle(np.mean(np.exp(2j * np.pi * phases)))
            scatter = float(np.sqrt(max(0.0, -2 * np.log(abs(np.mean(np.exp(2j * np.pi * phases)))))) / (2 * np.pi))
            phase = float((angle / (2 * np.pi)) % 1.0)
            offset = min(phase, 1 - phase)
            if scatter > PHASE_SCATTER_LIMIT or offset * reference <= long_.duration + short.duration:
                continue
            if abs(phase - 0.5) > SECONDARY_PHASE_WINDOW and offset < MIN_OTHER_PHASE:
                continue
            depth_long = long_.box_depth_ppm if long_.box_depth_ppm is not None else long_.reported_depth_ppm
            depth_short = short.box_depth_ppm if short.box_depth_ppm is not None else short.reported_depth_ppm
            secondary_like = abs(phase - 0.5) <= SECONDARY_PHASE_WINDOW
            out.append({
                "detections": [long_.label(), short.label()], "reference_period_days": reference,
                "phase_of_second_events": round(phase, 4), "phase_scatter": round(scatter, 4),
                "box_depths_ppm": [depth_long, depth_short],
                "depth_ratio": round(max(depth_long, depth_short) / min(depth_long, depth_short), 3)
                if min(depth_long, depth_short) > 0 else None,
                "interpretation": ("events interleave at phase ~0.5 of one period: consistent with primary and "
                                   "secondary eclipses of an eclipsing binary (a very shallow second dip could "
                                   "instead be a planetary occultation); also consider the half-period alias"
                                   if secondary_like else
                                   "events interleave at a fixed phase of one period: possible eccentric eclipsing "
                                   "binary, two signals sharing a period, or offset windows on one broad feature"),
            })
    return out


def member_summary(d: Detection) -> dict:
    return {"detection": d.label(), "period_days": d.period, "epoch": d.epoch, "duration_days": d.duration,
            "reported_depth_ppm": d.reported_depth_ppm, "depth_estimator": d.depth_estimator,
            "box_depth_ppm": d.box_depth_ppm, "local_depth_ppm": d.local_depth_ppm,
            "observed_events": len(d.events), "score": d.score, "odd_even_sigma": d.odd_even_sigma,
            "search_pass": d.search_pass, "masked_percent": d.masked_percent,
            "parity_coverage": parity_coverage(d)}


def depth_check(rows: list[Detection], unmasked: dict[str, float] | None = None) -> dict:
    """List every depth with its estimator; flag disagreement, never pick one."""
    issues = []
    for d in rows:
        direct = (unmasked or {}).get(d.label())
        box = d.box_depth_ppm if d.box_depth_ppm is not None else d.reported_depth_ppm
        if direct is not None and box and box > 0 and direct < UNMASKED_WARNING * box:
            where = (f" (found in search pass {d.search_pass} with {d.masked_percent:.1f}% of samples masked)"
                     if d.masked_percent is not None else "")
            issues.append(f"{d.label()}: depth {box:.1f} ppm but only {direct:.1f} ppm at this ephemeris in the "
                          f"unmasked light curve{where}; the dip depends on masking or detrending, not the data alone")
    for d in rows:
        if d.box_depth_ppm and d.local_depth_ppm is not None and d.box_depth_ppm > 0:
            ratio = d.local_depth_ppm / d.box_depth_ppm
            if ratio < LOCAL_BOX_WARNING:
                issues.append(f"{d.label()}: local depth {d.local_depth_ppm:.1f} ppm is {100 * ratio:.1f}% of "
                              f"box depth {d.box_depth_ppm:.1f} ppm (dip no sharper than nearby variation)")
    positive = [d.reported_depth_ppm for d in rows if d.reported_depth_ppm > 0]
    spread = max(positive) / min(positive) if len(positive) >= 2 else 1.0
    if spread > DEPTH_RATIO_LIMIT:
        issues.append(f"reported depths in one event pattern span {min(positive):.1f}-{max(positive):.1f} ppm "
                      f"({spread:.0f}x); they come from different estimators: "
                      + "; ".join(f"{d.label()} {d.reported_depth_ppm:.1f} ppm [{d.depth_estimator}]" for d in rows))
    return {"status": "INCONSISTENT" if issues else "consistent", "issues": issues,
            "note": "all values are kept; no depth is preferred by this check"}


def related_variability(det: Detection, variability: list[dict]) -> list[dict]:
    matches = []
    for row in variability:
        period = _number(row.get("period_days"))
        if not period:
            continue
        for ratio in RELATED_RATIOS:
            if abs(det.period / ratio - period) / period <= PERIOD_TOLERANCE:
                matches.append({"variability_rank": row.get("rank"), "variability_period_days": period,
                                "transit_period_over_variability_period": ratio,
                                "amplitude_ppm": _number(row.get("amplitude_ppm")),
                                "score": _number(row.get("score"))})
    return sorted(matches, key=lambda m: -(m["score"] or 0))


def segments(time: np.ndarray) -> list[np.ndarray]:
    order = np.argsort(time, kind="stable")
    breaks = np.flatnonzero(np.diff(time[order]) > SEGMENT_GAP_DAYS) + 1
    return np.split(order, breaks)


def _box_depth(residual: np.ndarray, inside: np.ndarray) -> float:
    return float(np.mean(residual[~inside]) - np.mean(residual[inside])) * 1e6


def in_transit(time: np.ndarray, det: Detection) -> np.ndarray:
    phase = ((time - det.epoch + det.period / 2) % det.period) - det.period / 2
    return np.abs(phase) <= det.duration / 2


def unmasked_depth(time: np.ndarray, flux: np.ndarray, det: Detection) -> float | None:
    """Box depth at the detection's ephemeris on all samples, after only a
    per-segment constant and linear term fitted out of transit."""
    inside = in_transit(time, det)
    if inside.sum() < 3 or (~inside).sum() < 10:
        return None
    residual = np.full_like(flux, np.nan)
    for index in segments(time):
        t = time[index] - time[index].mean()
        f = flux[index] / np.median(flux[index]) - 1.0
        base = np.column_stack([np.ones_like(t), t])
        rows = ~inside[index]
        if rows.sum() > 4:
            residual[index] = f - base @ np.linalg.lstsq(base[rows], f[rows], rcond=None)[0]
    usable = np.isfinite(residual)
    if (usable & inside).sum() < 3:
        return None
    return round(_box_depth(residual[usable], inside[usable]), 2)


def modulation_test(time: np.ndarray, flux: np.ndarray, det: Detection, variability_period: float,
                    harmonics: int = HARMONICS) -> dict:
    """Depth at the detection's ephemeris before and after removing
    `harmonics` harmonics of the modulation, fitted per segment to
    out-of-transit samples with a constant and linear term."""
    inside = in_transit(time, det)
    if inside.sum() < 3 or (~inside).sum() < 10:
        return {"status": "NOT_EVALUATED", "reason": "too few in- or out-of-transit samples"}
    before = np.empty_like(flux)
    after = np.empty_like(flux)
    fitted = 0
    for index in segments(time):
        t = time[index] - time[index].mean()
        f = flux[index] / np.median(flux[index]) - 1.0
        base = np.column_stack([np.ones_like(t), t])
        waves = [fn(2 * np.pi * k * time[index] / variability_period)
                 for k in range(1, harmonics + 1) for fn in (np.sin, np.cos)]
        full = np.column_stack([base, *waves])
        fit_rows = ~inside[index]
        if fit_rows.sum() <= full.shape[1] + 2:
            before[index] = after[index] = np.nan
            continue
        fitted += 1
        before[index] = f - base @ np.linalg.lstsq(base[fit_rows], f[fit_rows], rcond=None)[0]
        after[index] = f - full @ np.linalg.lstsq(full[fit_rows], f[fit_rows], rcond=None)[0]
    usable = np.isfinite(before) & np.isfinite(after)
    if (usable & inside).sum() < 3:
        return {"status": "NOT_EVALUATED", "reason": "no fitted segment contains in-transit samples"}
    depth_before = _box_depth(before[usable], inside[usable])
    depth_after = _box_depth(after[usable], inside[usable])
    retained = depth_after / depth_before if depth_before > 0 else None
    if retained is None:
        verdict = "no positive dip at this ephemeris before modulation removal"
    elif retained < EXPLAINED_BELOW:
        verdict = ("competing explanation: periodic modulation accounts for most of the dip "
                   "(stellar variability, ellipsoidal or reflection modulation, or an alias of it)")
    elif retained > SURVIVES_ABOVE:
        verdict = "dip survives removal of the modulation; the modulation does not explain it"
    else:
        verdict = "modulation explains part of the dip; both explanations remain open"
    return {"status": "OK", "variability_period_days": variability_period, "harmonics": harmonics,
            "segments_fitted": fitted, "in_transit_samples": int((usable & inside).sum()),
            "depth_before_ppm": round(depth_before, 2), "depth_after_ppm": round(depth_after, 2),
            "retained_fraction": None if retained is None else round(retained, 4), "verdict": verdict,
            "method": "per-segment least squares on out-of-transit samples: constant + linear "
                      f"(before) and + {harmonics} harmonics of the variability period (after)"}


def variability_conflicts(detections: list[Detection], variability: list[dict],
                          time: np.ndarray | None, flux: np.ndarray | None) -> list[dict]:
    out = []
    for det in detections:
        matches = related_variability(det, variability)
        entry = {"detection": det.label(), "period_days": det.period, "related_variability": matches[:3]}
        if not matches:
            entry["status"] = "no related periodic_variability signal reported"
        elif time is None or flux is None:
            entry["status"] = "related modulation reported; light curve not supplied, so not tested"
        else:
            entry["status"] = "tested"
            entry["modulation_test"] = modulation_test(time, flux, det, matches[0]["variability_period_days"])
        out.append(entry)
    return out


def diagnose(astrosearch_rows: list[dict], tls_rows: list[dict] | None = None,
             time: np.ndarray | None = None, flux: np.ndarray | None = None) -> dict:
    detections = from_astrosearch(astrosearch_rows)
    if tls_rows:
        if time is None:
            raise ValueError("TLS rows need the light-curve times to locate observed events")
        detections += from_tls(tls_rows, time)
    variability = [r for r in astrosearch_rows if r.get("class") == "periodic_variability"]
    unmasked = ({d.label(): unmasked_depth(time, flux, d) for d in detections}
                if time is not None and flux is not None else {})
    return {"detections": len(detections),
            "parameters": {"event_phase_tolerance": EVENT_PHASE_TOLERANCE, "min_event_overlap": MIN_EVENT_OVERLAP,
                           "depth_ratio_limit": DEPTH_RATIO_LIMIT, "local_box_warning": LOCAL_BOX_WARNING,
                           "related_ratios": list(RELATED_RATIOS), "period_tolerance": PERIOD_TOLERANCE,
                           "harmonics": HARMONICS, "explained_below": EXPLAINED_BELOW,
                           "survives_above": SURVIVES_ABOVE, "unmasked_warning": UNMASKED_WARNING,
                           "phase_scatter_limit": PHASE_SCATTER_LIMIT,
                           "secondary_phase_window": SECONDARY_PHASE_WINDOW,
                           "min_other_phase": MIN_OTHER_PHASE},
            "event_groups": event_groups(detections, {k: v for k, v in unmasked.items() if v is not None}),
            "interleaved_pairs": interleaved_pairs(detections),
            "variability_conflicts": variability_conflicts(detections, variability, time, flux),
            "provenance": "report-only: input candidate rows are not modified, filtered or re-ranked"}
