#!/usr/bin/env python3
"""Run the six-point CTOI readiness gate on one astrosearch transit candidate.

  1. Signal in 2+ sectors and 2+ data products with a consistent depth.
  2. A physical transit fit is preferred and its stellar density matches TIC.
  3. No secondary eclipse, odd/even difference, grazing shape, centroid
     offset, neighbor able to mimic the depth, or high Gaia RUWE.
  4. Not already a TOI, CTOI, SPOC TCE, or confirmed planet.
  5. Empirical per-target FAP and TRICERATOPS FPP below 1%.
  6. Ephemeris with uncertainties and predicted windows (publication and the
     ExoFOP submission remain human actions).

Writes gate.json, gate.md, ctoi_row.csv, and ephemeris_windows.csv to
--out-dir. Never submits anything and never calls a candidate a planet.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from vetting import aliases, crosscheck, diagnostics, fap, fpp, gate, provenance, siblings, stellar  # noqa: E402
from vetting.io import (  # noqa: E402
    TRANSIT_CLASSES, normalize_epoch, read_candidates, read_lightcurve, write_json)

ROOT = Path(__file__).resolve().parents[1]


def number(row: dict, key: str):
    try:
        value = float(row.get(key, ""))
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def physical_fit(lightcurve: Path, candidates: Path, candidate, out_dir: Path,
                 mcmc_steps: int) -> dict:
    """Fit the selected candidate with refine_candidates.py."""
    with candidates.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        rows = [r for r in reader if r.get("rank") == str(candidate.rank)]
        fields = reader.fieldnames
    single = out_dir / "candidate.csv"
    with single.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    output = out_dir / "fit.csv"
    subprocess.run([sys.executable, str(ROOT / "scripts/refine_candidates.py"),
                    "--lightcurve", str(lightcurve), "--candidates", str(single),
                    "--output", str(output), "--max-candidates", "1",
                    "--mcmc-steps", str(mcmc_steps)], check=True, capture_output=True, text=True)
    with output.open(newline="") as stream:
        row = next(csv.DictReader(stream))
    if row.get("fit_error"):
        return {"error": row["fit_error"]}
    fit = {
        "period": number(row, "fitted_period_days"),
        "period_sigma": number(row, "fitted_period_sigma_days_local"),
        "epoch": number(row, "fitted_epoch_btjd"),
        "epoch_sigma": number(row, "fitted_epoch_sigma_days_local"),
        "covariance": number(row, "fitted_epoch_period_covariance_local") or 0.0,
        "rp_rs": number(row, "fitted_radius_ratio"),
        "rp_rs_sigma": number(row, "fitted_radius_ratio_sigma_local"),
        "a_rs": number(row, "fitted_a_over_rstar"),
        "a_rs_sigma": number(row, "fitted_a_over_rstar_sigma_local"),
        "impact": number(row, "fitted_impact_parameter"),
        "impact_sigma": number(row, "fitted_impact_parameter_sigma_local"),
        "duration": number(row, "fitted_duration_days_circular_approx"),
        "delta_bic": number(row, "delta_bic_null_minus_transit"),
        "observed_windows": number(row, "observed_transit_windows"),
        "uncertainty_source": "local covariance",
    }
    # Local covariance understates a/R* (degenerate with impact parameter) and
    # ephemeris errors; when an MCMC sample exists, use its 16-84% half-widths.
    if row.get("mcmc_median_t0_period_rprs_a_rs_impact"):
        median = json.loads(row["mcmc_median_t0_period_rprs_a_rs_impact"])
        lower, upper = json.loads(row["mcmc_16_84_bounds"])
        names = ["epoch", "period", "rp_rs", "a_rs", "impact"]
        for i, name in enumerate(names):
            fit[name] = median[i]
            fit[f"{name}_sigma"] = 0.5 * (upper[i] - lower[i])
        fit["uncertainty_source"] = f"emcee posterior ({mcmc_steps} steps per walker)"
    # refine_candidates.py fits in the input's own time system; the gate works in BTJD.
    if fit.get("epoch") is not None:
        fit["epoch"] = normalize_epoch(fit["epoch"])
    return fit


def correct_dilution(lc, events: list[dict], contratio: float | None) -> tuple[list[dict], float]:
    """Scale depths from parts lacking a pipeline crowding correction by
    1 + contratio; returns corrected events and the mean correction factor."""
    if not contratio:
        return events, 1.0
    corrected, factors = [], []
    for e in events:
        meta = lc.parts[e["part"]] if lc.parts else {}
        factor = 1.0 if meta.get("crowding_corrected") else 1.0 + contratio
        factors.append(factor)
        corrected.append({**e, "depth_ppm": e["depth_ppm"] * factor,
                          "sigma_ppm": e["sigma_ppm"] * factor, "dilution_factor": factor})
    return corrected, float(np.mean(factors)) if factors else 1.0


def next_tess_sectors(tic_id: int, ra: float, dec: float) -> dict:
    try:
        from tess_stars2px import tess_stars2px_function_entry
    except ImportError:
        return {"status": gate.NE, "reason": "tess-point not installed"}
    result = tess_stars2px_function_entry(int(tic_id), float(ra), float(dec))
    sectors = sorted({int(s) for s in result[3] if int(s) > 0})
    return {"observing_sectors": sectors, "source": "tess-point (planned pointings)"}



def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tic", type=int, required=True)
    p.add_argument("--lightcurve", type=Path, required=True,
                   help="combined CSV searched by astrosearch (with .meta.json sidecar)")
    p.add_argument("--candidates", type=Path, required=True, help="astrosearch candidate CSV")
    p.add_argument("--rank", type=int, help="global rank of the candidate (default: best transit)")
    p.add_argument("--binary", type=Path, default=ROOT / "target/release/astrosearch")
    p.add_argument("--search-args", default="", help="astrosearch options used for the original search")
    p.add_argument("--products", type=Path, nargs="*", default=[],
                   help="independent product CSVs covering the signal (e.g. QLP and SPOC)")
    p.add_argument("--pixel-reports", type=Path, nargs="*", default=[],
                   help="pixel_vet.json files from pixel_difference_check.py")
    p.add_argument("--tce-sectors", default="")
    p.add_argument("--null-trials", type=int, default=300)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--threads", type=int, default=6)
    p.add_argument("--mcmc-steps", type=int, default=2000,
                   help="emcee steps per walker for posterior uncertainties (0: local covariance)")
    p.add_argument("--predict-days", type=float, default=730.0)
    p.add_argument("--offline", action="store_true", help="skip catalog, Gaia, and FPP queries")
    p.add_argument("--skip-fpp", action="store_true")
    p.add_argument("--detection-calibration", type=Path,
                   default=ROOT / "data/calibration/detection_threshold.json",
                   help="detection-tier calibration from calibrate_detection_threshold.py")
    p.add_argument("--sibling-score-fraction", type=float, default=0.5,
                   help="mask other transit candidates scoring at least this fraction of the "
                        "candidate's score (0 disables sibling masking)")
    p.add_argument("--out-dir", type=Path, required=True)
    a = p.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    search_args = shlex.split(a.search_args)

    lc = read_lightcurve(a.lightcurve)
    candidates = read_candidates(a.candidates, TRANSIT_CLASSES)
    if not candidates:
        raise SystemExit("no transit-class candidates in the candidate CSV")
    candidate = (next(c for c in candidates if c.rank == a.rank) if a.rank
                 else max(candidates, key=lambda c: c.score))
    duo = candidate.kind == "duo_transit"
    report = {"tic": a.tic, "candidate": {k: getattr(candidate, k) for k in
                                          ("rank", "kind", "score", "period", "epoch",
                                           "duration", "depth_ppm", "n_events")}}

    fit = physical_fit(a.lightcurve, a.candidates, candidate, a.out_dir, a.mcmc_steps)
    report["fit"] = fit
    period = fit.get("period") or candidate.period
    epoch = fit.get("epoch") or candidate.epoch
    duration = fit.get("duration") or candidate.duration

    star, neighbors, ruwe = {}, [], None
    observation_year = stellar.btjd_to_year(float(np.median(lc.time)))
    if not a.offline:
        try:
            tic = stellar.query_tic(a.tic, a.out_dir / "cache")
            star = tic["star"]
            neighbors = stellar.neighbors_at_epoch(star, tic["neighbors"], observation_year)
            gaia = stellar.query_gaia_ruwe(star, a.tic, a.out_dir / "cache")
            ruwe = gaia.get("ruwe")
            report["gaia"] = gaia
        except Exception as error:
            report["stellar_error"] = f"{type(error).__name__}: {error}"
    report["star"], report["ruwe"], report["observation_year"] = star, ruwe, observation_year

    # Other significant, non-commensurate transit signals are masked so they
    # cannot masquerade as this candidate's secondary or odd/even difference.
    full_lc = lc
    if a.sibling_score_fraction > 0:
        baseline = float(lc.time[-1] - lc.time[0])
        lc, report["siblings"] = siblings.mask_siblings(
            full_lc, candidate, candidates, baseline, min_score_fraction=a.sibling_score_fraction)

    raw_events = diagnostics.event_depths(lc, epoch, period, duration)
    events, depth_factor = correct_dilution(lc, raw_events, star.get("contratio"))
    depth = diagnostics.combine(events)
    observed_depth = diagnostics.combine(raw_events)
    report["events"], report["depth"], report["observed_depth"] = events, depth, observed_depth
    coverage = diagnostics.sector_coverage(lc, events)
    parts = diagnostics.part_consistency(lc, events)
    oe = diagnostics.odd_even(events)
    secondary = diagnostics.secondary_scan(lc, epoch, period, duration)
    quality = diagnostics.quality_overlap(lc, events, duration)
    shape = diagnostics.grazing(fit)
    planet = stellar.derive_planet(star, fit, depth_correction=depth_factor) if star else {}
    density = stellar.density_consistency(planet.get("transit_density_g_cm3"), star)
    shape_limit = stellar.shape_delta_mag_limit(fit, observed_depth["depth_ppm"] or 0.0)
    blends = stellar.blend_capability(star, neighbors, observed_depth["depth_ppm"] or 0.0,
                                      delta_mag_max=shape_limit["delta_mag_max"])
    report.update(coverage=coverage, parts=parts, odd_even=oe, secondary=secondary,
                  quality=quality, grazing=shape, planet=planet, density=density, blends=blends,
                  blend_shape_limit=shape_limit)

    if duo and len(candidate.event_epochs) >= 2:
        t1, t2 = sorted(candidate.event_epochs)[:2]
        report["duo_aliases"] = aliases.duo_aliases(
            lc, t1, t2, duration, depth["depth_ppm"] or candidate.depth_ppm,
            depth["sigma_ppm"] or candidate.depth_ppm / 3)

    report["fap"] = fap.empirical_fap(
        full_lc, a.binary, candidate.score, epoch, period, duration, search_args,
        trials=a.null_trials, workers=a.workers, threads=a.threads,
        work_dir=a.out_dir / "nulls") if a.null_trials > 0 else None
    report["bls"] = crosscheck.bls_recovery(lc, period, epoch, duration)
    report["products"] = (crosscheck.product_depths(lc, a.products, epoch, period, duration)
                          if a.products else None)

    sectors = sorted({int(m["sector"]) for m in lc.parts if m.get("sector") is not None})
    if a.offline:
        known = {"status": gate.NE, "reason": "offline"}
        fpp_result = {"status": gate.NE, "reason": "offline"}
    else:
        from vetting import catalogs
        tce = [int(s) for s in a.tce_sectors.split(",") if s.strip()] or sectors
        known = catalogs.known_object_check(a.tic, period, epoch, duration, tce)
        fpp_result = ({"status": gate.NE, "reason": "--skip-fpp"} if a.skip_fpp else
                      fpp.triceratops_fpp(a.tic, lc, sectors, period, epoch, duration,
                                          observed_depth["depth_ppm"] or candidate.depth_ppm))
    report["known"], report["fpp"] = known, fpp_result

    end = float(lc.time[-1])
    windows = aliases.predicted_windows(epoch, period, end, end + a.predict_days,
                                        fit.get("epoch_sigma") or 0.0, fit.get("period_sigma") or 0.0,
                                        fit.get("covariance") or 0.0, duration)
    if star.get("ra") is not None:
        report["tess_sectors"] = next_tess_sectors(a.tic, star["ra"], star["dec"])
    pixel_reports = [json.loads(Path(r).read_text()) for r in a.pixel_reports]

    criteria = {
        "1_multisector_multiproduct": gate.criterion_multisector(coverage, parts, report["products"]),
        "2_physical_fit_density": gate.criterion_physical_fit(fit, density),
        "3_no_false_positive_indicators": gate.criterion_false_positive(
            secondary, oe, shape, quality, blends, ruwe, pixel_reports, duo),
        "4_not_already_known": gate.criterion_not_known(known),
        "5_significance": gate.criterion_significance(
            report["fap"], fpp_result,
            tier=gate.detection_tier(candidate.score, gate.load_calibration(a.detection_calibration))),
        "6_preregistered_ephemeris": gate.criterion_ephemeris(fit, windows),
    }
    decision = gate.verdict(criteria)
    notes = (f"astrosearch {candidate.kind}; FAP {report['fap']['fap_empirical']:.3g}"
             if report["fap"] else f"astrosearch {candidate.kind}")
    ctoi = gate.ctoi_row(a.tic, {**fit, "period": period, "epoch": epoch, "duration": duration},
                         depth, planet, notes)
    report.update(criteria=criteria, verdict=decision, ctoi_row=ctoi, ephemeris_windows=windows)
    binary_provenance = Path(str(a.candidates) + ".provenance.json")
    report["provenance"] = provenance.record({
        "lightcurve": a.lightcurve, "candidates": a.candidates,
        "detection_calibration": a.detection_calibration, "binary": a.binary})
    if binary_provenance.exists():
        report["provenance"]["search"] = json.loads(binary_provenance.read_text())

    write_json(a.out_dir / "gate.json", report)
    (a.out_dir / "gate.md").write_text(gate.render_markdown(
        a.tic, report["candidate"], criteria, decision, planet, ctoi))
    with (a.out_dir / "ctoi_row.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(ctoi))
        writer.writeheader()
        writer.writerow(ctoi)
    if windows:
        with (a.out_dir / "ephemeris_windows.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(windows[0]))
            writer.writeheader()
            writer.writerows(windows)
    print(decision["summary"])
    for key, result in criteria.items():
        print(f"  {key}: {result['status']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
