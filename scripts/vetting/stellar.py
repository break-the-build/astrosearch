"""Stellar context (TIC 8.2, Gaia DR3) and derived planet parameters.

Network queries are cached as JSON. The derivations are pure functions with
Monte Carlo uncertainty propagation, so they are testable offline.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

G = 6.67430e-11
R_SUN_M = 6.957e8
M_SUN_KG = 1.98847e30
AU_M = 1.495978707e11
R_SUN_OVER_R_EARTH = 109.076
RHO_SUN_G_CM3 = 1.408
TEFF_SUN = 5772.0
TESS_PIXEL_ARCSEC = 21.0


# --- queries ---------------------------------------------------------------

def _cached(path: Path | None, fetch):
    if path and path.exists():
        return json.loads(path.read_text())
    value = fetch()
    if path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value, indent=1, sort_keys=True, default=float) + "\n")
    return value


def _row_dict(row, keys) -> dict:
    out = {}
    for key in keys:
        try:
            value = row[key]
        except KeyError:
            continue
        if np.ma.is_masked(value):
            out[key] = None
        elif isinstance(value, (np.integer, int)):
            out[key] = int(value)
        elif isinstance(value, (np.floating, float)):
            out[key] = None if not math.isfinite(float(value)) else float(value)
        else:
            out[key] = str(value)
    return out


TIC_KEYS = ["ID", "ra", "dec", "Tmag", "Teff", "e_Teff", "logg", "e_logg", "rad", "e_rad",
            "mass", "e_mass", "rho", "e_rho", "d", "e_d", "contratio", "GAIA", "objType",
            "disposition", "lumclass"]


def query_tic(tic_id: int, cache_dir: Path | None = None) -> dict:
    """TIC 8.2 row for the target plus every TIC source within 3 TESS pixels."""
    def fetch():
        from astroquery.mast import Catalogs
        target = Catalogs.query_criteria(catalog="Tic", ID=int(tic_id))
        if len(target) == 0:
            raise RuntimeError(f"TIC {tic_id} not found")
        star = _row_dict(target[0], TIC_KEYS)
        radius = 3 * TESS_PIXEL_ARCSEC / 3600.0
        cone = Catalogs.query_region(f"{star['ra']} {star['dec']}", radius=radius, catalog="Tic")
        neighbors = []
        for row in cone:
            entry = _row_dict(row, ["ID", "ra", "dec", "Tmag", "dstArcSec", "disposition"])
            if int(entry.get("ID") or 0) != int(tic_id):
                neighbors.append(entry)
        return {"star": star, "neighbors": neighbors, "source": "MAST TIC 8.2"}
    return _cached(cache_dir / f"tic{tic_id}_tic.json" if cache_dir else None, fetch)


def query_gaia_ruwe(ra: float, dec: float, tic_id: int, cache_dir: Path | None = None) -> dict:
    """Gaia DR3 RUWE of the brightest source within 2 arcsec of the TIC position."""
    def fetch():
        from astroquery.gaia import Gaia
        query = ("SELECT TOP 5 source_id, ruwe, phot_g_mean_mag, "
                 f"DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', {ra}, {dec})) AS sep "
                 "FROM gaiadr3.gaia_source "
                 f"WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {2/3600})) "
                 "ORDER BY phot_g_mean_mag ASC")
        rows = Gaia.launch_job(query).get_results()
        if len(rows) == 0:
            return {"ruwe": None, "source": "Gaia DR3", "note": "no source within 2 arcsec"}
        best = _row_dict(rows[0], ["source_id", "ruwe", "phot_g_mean_mag", "sep"])
        return {**best, "source": "Gaia DR3", "matches_within_2_arcsec": len(rows)}
    return _cached(cache_dir / f"tic{tic_id}_gaia.json" if cache_dir else None, fetch)


# --- derivations -------------------------------------------------------------

def _draw(rng, value, sigma, n, positive=True):
    if value is None or not math.isfinite(value):
        return None
    sigma = 0.0 if sigma is None or not math.isfinite(sigma) else abs(sigma)
    samples = rng.normal(value, sigma, n) if sigma > 0 else np.full(n, float(value))
    return np.where(samples > 0, samples, np.nan) if positive else samples


def _summary(samples) -> dict | None:
    if samples is None:
        return None
    finite = samples[np.isfinite(samples)]
    if finite.size == 0:
        return None
    lo, mid, hi = np.percentile(finite, [16, 50, 84])
    return {"median": float(mid), "minus": float(mid - lo), "plus": float(hi - mid)}


def transit_density_g_cm3(a_over_rstar: float, period_days: float) -> float:
    """Stellar density implied by a circular-orbit transit: 3 pi (a/R*)^3 / (G P^2)."""
    period_s = period_days * 86400.0
    return 3 * math.pi * a_over_rstar ** 3 / (G * period_s ** 2) / 1000.0


def derive_planet(star: dict, fit: dict, depth_correction: float = 1.0, albedo: float = 0.3,
                  draws: int = 20000, seed: int = 20261010) -> dict:
    """Planet radius, orbit, insolation, equilibrium temperature, and the
    transit-implied stellar density, each as median and 16/84 offsets.

    `depth_correction` multiplies the transit depth (1 + contamination ratio
    when the photometry was not dilution corrected); Rp/R* scales by its root.
    """
    rng = np.random.default_rng(seed)
    n = draws
    rp_rs = _draw(rng, fit.get("rp_rs"), fit.get("rp_rs_sigma"), n)
    a_rs = _draw(rng, fit.get("a_rs"), fit.get("a_rs_sigma"), n)
    period = _draw(rng, fit.get("period"), fit.get("period_sigma"), n)
    r_star = _draw(rng, star.get("rad"), star.get("e_rad"), n)
    m_star = _draw(rng, star.get("mass"), star.get("e_mass"), n)
    teff = _draw(rng, star.get("Teff"), star.get("e_Teff"), n)
    out = {"albedo_assumed": albedo, "depth_correction": depth_correction,
           "assumptions": "circular orbit, full heat redistribution, Gaussian TIC errors"}
    if rp_rs is not None:
        rp_rs = rp_rs * math.sqrt(depth_correction)
        out["rp_rs_corrected"] = _summary(rp_rs)
        if r_star is not None:
            out["radius_earth"] = _summary(rp_rs * r_star * R_SUN_OVER_R_EARTH)
    if period is not None and m_star is not None:
        a_m = (G * m_star * M_SUN_KG * (period * 86400.0) ** 2 / (4 * math.pi ** 2)) ** (1 / 3)
        out["semi_major_axis_au"] = _summary(a_m / AU_M)
        if r_star is not None and teff is not None:
            flux = (r_star * R_SUN_M / a_m) ** 2 * (teff / TEFF_SUN) ** 4 * (AU_M / R_SUN_M) ** 2
            out["insolation_earth"] = _summary(flux)
            out["equilibrium_temperature_k"] = _summary(
                teff * np.sqrt(r_star * R_SUN_M / (2 * a_m)) * (1 - albedo) ** 0.25)
    if a_rs is not None and period is not None:
        density = 3 * np.pi * a_rs ** 3 / (G * (period * 86400.0) ** 2) / 1000.0
        out["transit_density_g_cm3"] = _summary(density)
    return out


def density_consistency(transit: dict | None, star: dict) -> dict:
    """Compare transit-implied and catalog stellar density in log space.

    CONSISTENT within 3 sigma; TENSION beyond 3 sigma but within a factor 5
    (eccentric orbits shift the transit density); INCONSISTENT otherwise,
    which points to a blend, a giant host, or a wrong catalog parameter.
    """
    rho = star.get("rho")
    if transit is None or rho is None:
        return {"status": "NOT_EVALUATED", "reason": "missing fit a/R* or TIC density"}
    catalog = rho * RHO_SUN_G_CM3
    catalog_sigma = (star.get("e_rho") or 0.0) * RHO_SUN_G_CM3
    fitted = transit["median"]
    fitted_sigma = 0.5 * (transit["minus"] + transit["plus"])
    sigma_log = math.hypot(fitted_sigma / fitted, catalog_sigma / catalog) / math.log(10)
    log_ratio = math.log10(fitted / catalog)
    z = log_ratio / sigma_log if sigma_log > 0 else math.inf
    ratio = 10 ** log_ratio
    if abs(z) <= 3:
        status = "CONSISTENT"
    elif 0.2 <= ratio <= 5:
        status = "TENSION"
    else:
        status = "INCONSISTENT"
    return {"status": status, "transit_g_cm3": fitted, "catalog_g_cm3": catalog,
            "ratio": ratio, "z_log": z}


def blend_capability(star: dict, neighbors: list[dict], depth_ppm: float,
                     radius_pixels: float = 2.0) -> list[dict]:
    """Neighbors within `radius_pixels` TESS pixels that could produce the
    observed depth if fully eclipsed: F_nb / (F_target + F_nb) >= depth."""
    tmag = star.get("Tmag")
    if tmag is None:
        return []
    out = []
    for nb in neighbors:
        separation, nb_tmag = nb.get("dstArcSec"), nb.get("Tmag")
        if separation is None or nb_tmag is None or separation > radius_pixels * TESS_PIXEL_ARCSEC:
            continue
        ratio = 10 ** (-0.4 * (nb_tmag - tmag))
        max_depth_ppm = ratio / (1 + ratio) * 1e6
        out.append({"tic": nb.get("ID"), "separation_arcsec": separation, "tmag": nb_tmag,
                    "delta_tmag": nb_tmag - tmag, "max_eclipse_depth_ppm": max_depth_ppm,
                    "can_produce_depth": max_depth_ppm >= depth_ppm})
    return sorted(out, key=lambda e: e["separation_arcsec"])
