#!/usr/bin/env python3
"""Rank stars by the expected number of new transiting planets this search can find.

For each star, the planet population is integrated over period and radius:

    Y = sum_{P,R} f(P, R | Teff) dlnP dlnR * p_transit * C_ours(SNR) * (1 - C_prior(SNR))

  f          SAG13 occurrence (Kopparapu et al. 2018, ApJ 856, 122) with the
             inner-period fall-off below ~10 d (Mulders et al. 2015, ApJ 798,
             112: P^1.5) and the stellar-type factors of Mulders et al. 2015
             (ApJ 814, 130): about 2x more small planets around M dwarfs than
             G dwarfs, fewer giants.
  p_transit  (R* + Rp) / a, or the coplanar conditional probability when the
             star already has a known transiting planet.
  SNR        predicted from TESS magnitude with 1-hour noise measured on 384
             QLP light curves (data/calibration/noise_vs_tmag.json), or a
             measured cdpp_1h_ppm column, on the same white-noise scale as
             scripts/injection_recovery.py.
  C_ours     measured injection-recovery completeness of the deployed
             detection tier (data/calibration/completeness.json).
  C_prior    assumed completeness of SPOC and QLP searches
             (data/calibration/discovery_priors.json) so already-searched
             parameter space does not count as new.

Y is an expected count; 1 - exp(-Y) is the probability of at least one new
detectable planet under this model. It is a target-selection heuristic built
on population statistics, never evidence that any given star has a planet.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np

CALIBRATION_ENV = "ASTROSEARCH_CALIBRATION_DIR"
DEFAULT_CALIBRATION = Path(__file__).resolve().parents[1] / "data" / "calibration"
ALIASES = {"teff_k": ("teff_k", "teff", "Teff"), "radius_solar": ("radius_solar", "rad", "radius"),
           "tmag": ("tmag", "Tmag"), "mass_solar": ("mass_solar", "mass", "Mstar", "mstar"),
           "stellar_density_g_cm3": ("stellar_density_g_cm3", "rho_g_cm3", "rho"),
           "cdpp_1h_ppm": ("cdpp_1h_ppm", "cdpp1h_ppm", "cdpp_ppm"),
           "contamination_ratio": ("contamination_ratio", "contratio", "contam_ratio"),
           "observed_days": ("observed_days",), "sectors": ("sectors", "n_sectors"),
           "ra": ("ra", "RA", "ra_deg"), "dec": ("dec", "Dec", "DEC", "dec_deg")}
SOLAR_DENSITY = 1.41
R_EARTH_IN_SOLAR = 0.009168
R_SUN_IN_AU = 0.00465047
DAYS_PER_SECTOR = 24.0  # usable days per 27.4-day sector after gaps
# SAG13 (Kopparapu et al. 2018, Table 3): dN/dlnR dlnP = Gamma R^alpha P^beta, P in years.
SAG13 = ((3.4, 0.38, -0.19, 0.26), (math.inf, 0.73, -1.18, 0.59))
INNER_BREAK_DAYS, INNER_SLOPE = 10.0, 1.5
SMALL_PLANET_LIMIT = 3.4
MUTUAL_INCLINATION_DEG = 1.5  # Fabrycky et al. 2014 (ApJ 790, 146): ~1-2 degrees


@dataclass(frozen=True)
class Pipeline:
    name: str
    threshold: float
    width: float
    efficiency: float
    min_tmag: float
    max_tmag: float
    applies_if_column: str | None


@dataclass(frozen=True)
class YieldModel:
    noise_tmag: np.ndarray
    noise_log_sigma: np.ndarray
    beta: np.ndarray
    completeness_snr: np.ndarray
    completeness: np.ndarray
    pipelines: tuple[Pipeline, ...]


@dataclass(frozen=True)
class Settings:
    min_period: float = 0.5
    max_period: float = 30.0
    min_transits: float = 3.0
    assumed_sectors: float = 4.0
    period_points: int = 48
    radius_points: int = 36
    min_radius_earth: float = 0.5
    max_radius_earth: float = 20.0


def value(row: dict, name: str, positive: bool = True) -> float | None:
    for key in ALIASES.get(name, (name,)):
        try:
            x = float(row.get(key, ""))
        except (TypeError, ValueError):
            continue
        if math.isfinite(x) and (x > 0 or (not positive and x >= 0)):
            return x
    return None


def truthy(raw) -> bool:
    return str(raw or "").strip().lower() in ("1", "true", "yes", "y", "t")


def load_model(directory: Path | None = None) -> YieldModel:
    root = Path(directory or os.environ.get(CALIBRATION_ENV) or DEFAULT_CALIBRATION)
    try:
        noise = json.loads((root / "noise_vs_tmag.json").read_text())["table"]
        curve = json.loads((root / "completeness.json").read_text())["curve"]
        priors = json.loads((root / "discovery_priors.json").read_text())["pipelines"]
    except (OSError, KeyError, json.JSONDecodeError) as error:
        raise ValueError(f"target-yield calibration unreadable in {root} ({error}); "
                         f"set {CALIBRATION_ENV} or pass the calibration directory") from error
    if len(noise) < 2 or len(curve) < 2:
        raise ValueError("target-yield calibration needs at least two noise and completeness points")
    noise = sorted(noise, key=lambda r: r["tmag"])
    curve = sorted(curve, key=lambda r: r["snr"])
    pipelines = tuple(Pipeline(p["name"], float(p["threshold"]), float(p.get("width", 1.0)),
                               float(p["efficiency"]), float(p.get("min_tmag", -99)),
                               float(p.get("max_tmag", 99)), p.get("applies_if_column"))
                      for p in priors)
    # Noise cannot fall with magnitude; small dips between bins are scatter.
    sigma = np.maximum.accumulate(np.array([r["sigma_1h_ppm"] for r in noise], float))
    return YieldModel(np.array([r["tmag"] for r in noise], float), np.log(sigma * 1e-6),
                      np.array([r.get("beta_1h", 1.0) for r in noise], float),
                      np.array([r["snr"] for r in curve], float),
                      np.array([r["recovered"] / r["n"] for r in curve], float), pipelines)


def noise_at(model: YieldModel, tmag: float) -> tuple[float, float]:
    """1-hour noise (fractional) and beta, log-linear in Tmag. Brighter than
    the table the noise stays at the brightest bin (a systematics floor);
    fainter, the last segment is extrapolated."""
    x, y = model.noise_tmag, model.noise_log_sigma
    if tmag >= x[-1]:
        log_sigma = y[-1] + (tmag - x[-1]) * (y[-1] - y[-2]) / (x[-1] - x[-2])
    else:
        log_sigma = float(np.interp(tmag, x, y))
    return math.exp(log_sigma), float(np.interp(tmag, x, model.beta))


def occurrence(period: np.ndarray, radius: np.ndarray, teff: float | None) -> np.ndarray:
    """Planets per star per dlnP dlnR."""
    small = radius < SMALL_PLANET_LIMIT
    gamma = np.where(small, SAG13[0][1], SAG13[1][1])
    alpha = np.where(small, SAG13[0][2], SAG13[1][2])
    beta = np.where(small, SAG13[0][3], SAG13[1][3])
    years = np.maximum(period, INNER_BREAK_DAYS) / 365.25
    rate = gamma * radius ** alpha * years ** beta
    rate = rate * np.minimum(1.0, period / INNER_BREAK_DAYS) ** INNER_SLOPE
    if teff is not None:
        small_factor = 2.0 if teff < 3900 else 1.4 if teff < 5300 else 1.0 if teff < 6000 else 0.67
        giant_factor = 0.5 if teff < 3900 else 1.0
        rate = rate * np.where(small, small_factor, giant_factor)
    return rate


def prior_completeness(model: YieldModel, snr_red: np.ndarray, tmag: float, row: dict) -> np.ndarray:
    missed = np.ones_like(snr_red)
    for p in model.pipelines:
        if not p.min_tmag <= tmag < p.max_tmag:
            continue
        if p.applies_if_column and not truthy(row.get(p.applies_if_column)):
            continue
        found = p.efficiency / (1.0 + np.exp(-(snr_red - p.threshold) / p.width))
        missed = missed * (1.0 - found)
    return 1.0 - missed


def known_periods(row: dict) -> list[float]:
    raw = str(row.get("known_planet_periods", "") or "")
    periods = []
    for token in raw.replace(",", ";").split(";"):
        try:
            x = float(token)
        except ValueError:
            continue
        if math.isfinite(x) and x > 0:
            periods.append(x)
    return periods


def expected_new_planets(row: dict, model: YieldModel, settings: Settings = Settings()) -> dict | None:
    """Expected count of new planets detectable around one star, or None if
    the star lacks a radius or any noise estimate."""
    radius = value(row, "radius_solar")
    tmag, cdpp = value(row, "tmag"), value(row, "cdpp_1h_ppm")
    if radius is None or (tmag is None and cdpp is None):
        return None
    density = value(row, "stellar_density_g_cm3")
    mass = (density / SOLAR_DENSITY * radius ** 3 if density else
            value(row, "mass_solar") or radius)  # M ~ R for dwarfs without a mass
    sigma, beta = noise_at(model, tmag) if tmag is not None else (None, float(np.median(model.beta)))
    if cdpp is not None:
        sigma = cdpp * 1e-6
    dilution = 1.0 / (1.0 + (value(row, "contamination_ratio", positive=False) or 0.0))
    days = value(row, "observed_days") or DAYS_PER_SECTOR * (value(row, "sectors") or settings.assumed_sectors)
    tmag_for_priors = tmag if tmag is not None else 99.0

    lnp = np.linspace(math.log(settings.min_period), math.log(settings.max_period), settings.period_points)
    lnr = np.linspace(math.log(settings.min_radius_earth), math.log(settings.max_radius_earth),
                      settings.radius_points)
    period, rp = np.exp(lnp)[:, None], np.exp(lnr)[None, :]
    cell = (lnp[1] - lnp[0]) * (lnr[1] - lnr[0])
    a_over_r = (mass * (period / 365.25) ** 2) ** (1 / 3) / (radius * R_SUN_IN_AU)
    k = rp * R_EARTH_IN_SOLAR / radius
    geometric = np.clip((1 + k) / a_over_r, 0.0, 1.0)
    neighbours = known_periods(row)
    if neighbours:
        sigma_i = math.radians(MUTUAL_INCLINATION_DEG)
        # A known transiting planet means the system is seen nearly edge-on.
        geometric = np.vectorize(math.erf)((1 + k) / a_over_r / (math.sqrt(2) * sigma_i))
        near_known = np.any([np.abs(np.log(period / p)) < 0.02 for p in neighbours], axis=0)
        geometric = np.where(near_known, 0.0, geometric)
    duration_h = 6.0 * period / a_over_r  # mean chord: (pi/4) * central duration
    transits = days / period
    snr_white = 0.85 * k ** 2 * dilution * np.sqrt(transits * duration_h) * beta / sigma
    ours = np.interp(snr_white, model.completeness_snr, model.completeness, left=0.0, right=1.0)
    ours = np.where(transits >= settings.min_transits, ours, 0.0)
    prior = prior_completeness(model, snr_white / beta, tmag_for_priors, row)
    weight = occurrence(period, rp, value(row, "teff_k")) * cell * geometric
    detectable = float((weight * ours).sum())
    new = float((weight * ours * (1.0 - prior)).sum())
    return {"expected_new_planets": new, "probability_new_planet": 1.0 - math.exp(-new),
            "expected_detectable_planets": detectable, "sigma_1h_ppm": sigma * 1e6,
            "noise_source": "cdpp_1h_ppm" if cdpp is not None else "tmag_calibration",
            "observed_days": days, "known_transiting_host": bool(neighbours)}


def sky_cells(rows: list[dict], results: list[dict | None], cell_deg: float) -> list[dict]:
    """Sum expected new planets in roughly equal-area RA/Dec cells."""
    cells: dict[tuple[int, int], dict] = {}
    for row, result in zip(rows, results):
        ra, dec = value(row, "ra", positive=False), row_dec(row)
        if result is None or ra is None or dec is None:
            continue
        band = int(math.floor((dec + 90.0) / cell_deg))
        centre = -90.0 + (band + 0.5) * cell_deg
        width = cell_deg / max(math.cos(math.radians(centre)), cell_deg / 360.0)
        column = int(math.floor((ra % 360.0) / width))
        slot = cells.setdefault((band, column), {
            "dec_center": round(centre, 3), "ra_center": round(min((column + 0.5) * width, 360.0), 3),
            "targets": 0, "expected_new_planets": 0.0, "best_tic_id": "", "best_expected": -1.0})
        slot["targets"] += 1
        slot["expected_new_planets"] += result["expected_new_planets"]
        if result["expected_new_planets"] > slot["best_expected"]:
            slot["best_expected"] = result["expected_new_planets"]
            slot["best_tic_id"] = row.get("tic_id", "")
    return sorted(cells.values(), key=lambda c: -c["expected_new_planets"])


def row_dec(row: dict) -> float | None:
    for key in ALIASES["dec"]:
        try:
            x = float(row.get(key, ""))
        except (TypeError, ValueError):
            continue
        if math.isfinite(x) and -90.0 <= x <= 90.0:
            return x
    return None


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("targets", type=Path, help="CSV with tic_id, radius_solar, tmag (or cdpp_1h_ppm); "
                   "optional teff_k, mass_solar, stellar_density_g_cm3, contamination_ratio, "
                   "observed_days or sectors, spoc_2min, known_planet_periods, ra, dec")
    p.add_argument("--output", type=Path, required=True, help="ranked CSV")
    p.add_argument("--sky-cells", type=Path, help="also write expected yield per sky cell")
    p.add_argument("--cell-degrees", type=float, default=10.0)
    p.add_argument("--calibration-dir", type=Path)
    p.add_argument("--min-period", type=float, default=Settings.min_period)
    p.add_argument("--max-period", type=float, default=Settings.max_period)
    p.add_argument("--assumed-sectors", type=float, default=Settings.assumed_sectors,
                   help="sectors assumed when a row has neither observed_days nor sectors")
    a = p.parse_args(argv)
    if not 0 < a.min_period < a.max_period:
        p.error("need 0 < --min-period < --max-period")
    if not 0 < a.cell_degrees <= 90:
        p.error("--cell-degrees must be in (0, 90]")
    try:
        model = load_model(a.calibration_dir)
    except ValueError as error:
        p.error(str(error))
    settings = Settings(min_period=a.min_period, max_period=a.max_period, assumed_sectors=a.assumed_sectors)
    with a.targets.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    results = [expected_new_planets(row, model, settings) for row in rows]
    ranked = sorted(range(len(rows)), key=lambda i: (results[i] is None,
                                                     -(results[i] or {}).get("expected_new_planets", 0.0), i))
    fields = ["rank", "tic_id", "expected_new_planets", "probability_new_planet",
              "expected_detectable_planets", "sigma_1h_ppm", "noise_source", "observed_days",
              "known_transiting_host"]
    write_csv(a.output, fields, [{"rank": n + 1, "tic_id": rows[i].get("tic_id", ""), **(results[i] or {})}
                                 for n, i in enumerate(ranked)])
    if a.sky_cells:
        write_csv(a.sky_cells, ["dec_center", "ra_center", "targets", "expected_new_planets",
                                "best_tic_id", "best_expected"], sky_cells(rows, results, a.cell_degrees))
    scored = [r for r in results if r]
    print(f"ranked {len(scored)} of {len(rows)} targets; expected new planets in total "
          f"{sum(r['expected_new_planets'] for r in scored):.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
