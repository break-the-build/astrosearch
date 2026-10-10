"""The six-point path to a credible CTOI, as pure functions over check results.

Each criterion is PASS, FAIL, NOT_EVALUATED (missing data, tool, or network),
or ACTION_REQUIRED (a human step). A candidate is ready for human review only
when every automated criterion passes; the gate never submits anything.
"""
from __future__ import annotations

import math

PASS, FAIL, NE, ACTION = "PASS", "FAIL", "NOT_EVALUATED", "ACTION_REQUIRED"
TESS_PIXEL_ARCSEC = 21.0
RUWE_LIMIT = 1.4


def worst(statuses) -> str:
    statuses = list(statuses)
    if FAIL in statuses:
        return FAIL
    if NE in statuses or not statuses:
        return NE
    return PASS


def criterion_multisector(coverage: dict, parts: dict, products: dict | None) -> dict:
    sectors = coverage.get("sectors_with_events", [])
    authors = coverage.get("authors_with_events", [])
    checks = {
        "two_sectors_with_events": PASS if len(sectors) >= 2 else FAIL,
        "depth_consistent_across_sectors": parts.get("status", NE),
        "two_independent_products": (PASS if len(authors) >= 2 or (products or {}).get("status") == PASS
                                     else (products or {}).get("status", FAIL)),
    }
    return {"status": worst(checks.values()), "checks": checks, "sectors": sectors,
            "authors": authors}


def criterion_physical_fit(fit: dict, density: dict, min_delta_bic: float = 10.0) -> dict:
    delta_bic = fit.get("delta_bic")
    bic = NE if delta_bic is None else (PASS if delta_bic >= min_delta_bic else FAIL)
    density_status = {"CONSISTENT": PASS, "TENSION": FAIL, "INCONSISTENT": FAIL}.get(
        density.get("status"), NE)
    return {"status": worst([bic, density_status]),
            "checks": {"transit_model_preferred_delta_bic": bic,
                       "stellar_density_consistent": density_status},
            "delta_bic": delta_bic, "density": density}


def criterion_false_positive(secondary: dict, odd_even: dict, grazing: dict, quality: dict,
                             blends: list[dict], ruwe: float | None,
                             pixel_reports: list[dict], duo: bool) -> dict:
    centroid = NE
    if pixel_reports:
        offsets = [r.get("centroid_separation_arcsec") for r in pixel_reports]
        offsets = [o for o in offsets if o is not None and math.isfinite(o)]
        centroid = (FAIL if any(o > TESS_PIXEL_ARCSEC for o in offsets) else PASS) if offsets else NE
    checks = {
        "no_secondary_eclipse": secondary.get("status", NE),
        # A duo has one transit of each parity at most, so odd/even is undefined.
        "odd_even_consistent": PASS if duo and odd_even.get("status") == NE else odd_even.get("status", NE),
        "not_grazing": grazing.get("status", NE),
        "not_dominated_by_momentum_dumps_or_scattered_light": (
            PASS if quality.get("status") == NE else quality.get("status", NE)),
        "no_centroid_offset": centroid,
        "no_neighbor_able_to_mimic_depth": FAIL if any(b["can_produce_depth"] for b in blends) else PASS,
        "ruwe_below_1_4": NE if ruwe is None else (PASS if ruwe < RUWE_LIMIT else FAIL),
    }
    return {"status": worst(checks.values()), "checks": checks,
            "blending_neighbors": [b for b in blends if b["can_produce_depth"]]}


def criterion_not_known(known: dict) -> dict:
    return {"status": known.get("status", NE), "details": known}


def criterion_significance(fap: dict | None, fpp: dict | None, limit: float = 0.01) -> dict:
    fap_status = NE
    if fap:
        fap_status = PASS if fap["fap_upper_95"] < limit else FAIL
    fpp_status = NE if not fpp else fpp.get("status", NE)
    return {"status": worst([fap_status, fpp_status]),
            "checks": {"empirical_fap_upper95_below_1pct": fap_status,
                       "triceratops_fpp_below_1pct": fpp_status},
            "fap": fap, "fpp": fpp}


def criterion_ephemeris(fit: dict, windows: list[dict]) -> dict:
    sigmas = [fit.get("epoch_sigma"), fit.get("period_sigma")]
    ok = all(s is not None and math.isfinite(s) and s > 0 for s in sigmas) and bool(windows)
    return {"status": ACTION if ok else FAIL,
            "action": "publish the predicted windows with a timestamp (e.g. Zenodo/OSF), then "
                      "submit the CTOI on ExoFOP; both are human steps",
            "next_windows": windows[:5]}


def verdict(criteria: dict) -> dict:
    automated = {k: v["status"] for k, v in criteria.items() if v["status"] != ACTION}
    blockers = [k for k, s in automated.items() if s != PASS]
    ready = not blockers
    return {"ready_for_human_review": ready, "blockers": blockers,
            "summary": ("All automated checks pass. Ready for human review and the remaining "
                        "manual steps; this is still a candidate, not a planet."
                        if ready else
                        "Not ready for a CTOI: " + ", ".join(blockers) + ".")}


def ctoi_row(tic_id: int, fit: dict, depth: dict, planet: dict, notes: str) -> dict:
    """Fields of the ExoFOP CTOI upload form (BJD_TDB epochs)."""
    radius = planet.get("radius_earth") or {}
    duration_h = (fit.get("duration") or 0) * 24
    return {
        "tic_id": tic_id,
        "period_days": fit.get("period"), "period_unc_days": fit.get("period_sigma"),
        "epoch_bjd": None if fit.get("epoch") is None else fit["epoch"] + 2_457_000,
        "epoch_unc_days": fit.get("epoch_sigma"),
        "depth_ppm": depth.get("depth_ppm"), "depth_unc_ppm": depth.get("sigma_ppm"),
        "duration_hours": duration_h,
        "planet_radius_earth": radius.get("median"),
        "planet_radius_unc_earth": (None if not radius else 0.5 * (radius["minus"] + radius["plus"])),
        "notes": notes,
    }


def render_markdown(tic_id: int, candidate: dict, criteria: dict, decision: dict,
                    planet: dict, ctoi: dict) -> str:
    lines = [f"# CTOI gate: TIC {tic_id}", "",
             f"Candidate: {candidate['kind']} rank {candidate['rank']}, P = {candidate['period']:.6f} d, "
             f"score {candidate['score']:.2f}", "",
             f"**Verdict:** {decision['summary']}", "",
             "> Automated screening only. Passing every check does not make this a planet; "
             "it makes it worth a human vetting review and follow-up.", "",
             "| # | Criterion | Status |", "| --- | --- | --- |"]
    names = {"1_multisector_multiproduct": "Signal in 2+ sectors and 2+ products, consistent depth",
             "2_physical_fit_density": "Physical transit fit preferred; stellar density consistent",
             "3_no_false_positive_indicators": "No secondary, odd/even, grazing, centroid, neighbor, RUWE flags",
             "4_not_already_known": "Not a known TOI, CTOI, TCE, or confirmed planet",
             "5_significance": "Empirical FAP and TRICERATOPS FPP below 1%",
             "6_preregistered_ephemeris": "Ephemeris with uncertainties, predicted windows"}
    for i, (key, result) in enumerate(criteria.items(), start=1):
        lines.append(f"| {i} | {names.get(key, key)} | {result['status']} |")
    lines += ["", "## Sub-checks", ""]
    for key, result in criteria.items():
        for check, status in (result.get("checks") or {}).items():
            lines.append(f"- {key} / {check}: **{status}**")
    lines += ["", "## Derived planet parameters (if on target)", ""]
    for key in ("radius_earth", "semi_major_axis_au", "insolation_earth",
                "equilibrium_temperature_k", "transit_density_g_cm3"):
        value = planet.get(key)
        if value:
            lines.append(f"- {key}: {value['median']:.4g} (-{value['minus']:.2g} / +{value['plus']:.2g})")
    lines += ["", "## CTOI fields (draft, for human review)", ""]
    lines += [f"- {k}: {v}" for k, v in ctoi.items()]
    return "\n".join(lines) + "\n"
