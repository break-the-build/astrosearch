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


TIC_KEYS = ["ID", "ra", "dec", "pmRA", "pmDEC", "Tmag", "Teff", "e_Teff", "logg", "e_logg",
            "rad", "e_rad", "mass", "e_mass", "rho", "e_rho", "d", "e_d", "contratio", "GAIA",
            "objType", "disposition", "lumclass"]
#: TIC 8 positions are at epoch J2000; Gaia DR3 positions at J2016.0.
TIC_EPOCH_YEAR = 2000.0
GAIA_DR3_EPOCH_YEAR = 2016.0
BTJD_OFFSET = 2457000.0
JD_J2000 = 2451545.0


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
            entry = _row_dict(row, ["ID", "ra", "dec", "pmRA", "pmDEC", "Tmag", "dstArcSec",
                                    "disposition"])
            if int(entry.get("ID") or 0) != int(tic_id):
                neighbors.append(entry)
        return {"star": star, "neighbors": neighbors, "source": "MAST TIC 8.2"}
    # v2 adds proper motions; older caches without them are not reused.
    return _cached(cache_dir / f"tic{tic_id}_tic_v2.json" if cache_dir else None, fetch)


def gaia_query_by_dr2_id(dr2_id) -> str:
    """ADQL for the Gaia DR3 counterpart of a Gaia DR2 source (the TIC stores
    DR2 identifiers), nearest match first."""
    text = str(dr2_id).strip()
    if not text.isdigit():
        raise ValueError(f"Gaia DR2 source_id must be an integer, got {dr2_id!r}")
    return ("SELECT g.source_id, g.ruwe, g.phot_g_mean_mag, n.angular_distance AS sep "
            "FROM gaiadr3.dr2_neighbourhood AS n "
            "JOIN gaiadr3.gaia_source AS g ON g.source_id = n.dr3_source_id "
            f"WHERE n.dr2_source_id = {text} ORDER BY n.angular_distance ASC")


def _cone_query(ra: float, dec: float, radius_arcsec: float) -> str:
    return ("SELECT TOP 5 source_id, ruwe, phot_g_mean_mag, "
            f"DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', {ra}, {dec})) AS sep "
            "FROM gaiadr3.gaia_source "
            f"WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, "
            f"{radius_arcsec / 3600})) ORDER BY phot_g_mean_mag ASC")


def query_gaia_ruwe(star: dict, tic_id: int, cache_dir: Path | None = None) -> dict:
    """Gaia DR3 RUWE of the TIC target.

    Matches by the TIC's Gaia DR2 source_id through `gaiadr3.dr2_neighbourhood`;
    without one, cones 2 arcsec around the TIC position propagated with its
    proper motion to the Gaia epoch (J2016.0). A J2000 cone alone misses
    high-proper-motion stars (LHS 3844 moved ~20 arcsec)."""
    def fetch():
        from astroquery.gaia import Gaia
        if star.get("GAIA"):
            rows = Gaia.launch_job(gaia_query_by_dr2_id(star["GAIA"])).get_results()
            if len(rows):
                best = _row_dict(rows[0], ["source_id", "ruwe", "phot_g_mean_mag", "sep"])
                return {**best, "source": "Gaia DR3", "match": "dr2_neighbourhood",
                        "dr2_source_id": str(star["GAIA"])}
        ra, dec = propagate_position(star["ra"], star["dec"], star.get("pmRA"), star.get("pmDEC"),
                                     TIC_EPOCH_YEAR, GAIA_DR3_EPOCH_YEAR)
        rows = Gaia.launch_job(_cone_query(ra, dec, 2.0)).get_results()
        if len(rows) == 0:
            return {"ruwe": None, "source": "Gaia DR3", "match": "cone (J2016 position)",
                    "note": "no source within 2 arcsec of the J2016 position"}
        best = _row_dict(rows[0], ["source_id", "ruwe", "phot_g_mean_mag", "sep"])
        return {**best, "source": "Gaia DR3", "match": "cone (J2016 position)",
                "matches_within_2_arcsec": len(rows)}
    return _cached(cache_dir / f"tic{tic_id}_gaia_dr3.json" if cache_dir else None, fetch)


# --- astrometry --------------------------------------------------------------

def btjd_to_year(btjd: float) -> float:
    """Julian epoch year of a BTJD time (TDB; time-scale offsets are negligible here)."""
    return 2000.0 + (btjd + BTJD_OFFSET - JD_J2000) / 365.25


def propagate_position(ra: float, dec: float, pm_ra, pm_dec, from_year: float,
                       to_year: float) -> tuple[float, float]:
    """Linear proper-motion propagation; pm_ra includes cos(dec), in mas/yr (TIC, Gaia)."""
    if pm_ra is None or pm_dec is None:
        return ra, dec
    years = to_year - from_year
    dec_new = dec + pm_dec * years / 3.6e6
    ra_new = ra + pm_ra * years / 3.6e6 / max(math.cos(math.radians(dec)), 1e-9)
    return ra_new, dec_new


def _separation_arcsec(ra1, dec1, ra2, dec2) -> float:
    r1, d1, r2, d2 = map(math.radians, (ra1, dec1, ra2, dec2))
    hav = math.sin((d2 - d1) / 2) ** 2 + math.cos(d1) * math.cos(d2) * math.sin((r2 - r1) / 2) ** 2
    return math.degrees(2 * math.asin(min(1.0, math.sqrt(hav)))) * 3600


def neighbors_at_epoch(star: dict, neighbors: list[dict], year: float) -> list[dict]:
    """Neighbor separations at `year`, propagating the target and each neighbor
    from the TIC epoch. Entries without positions keep their J2000 `dstArcSec`."""
    if star.get("ra") is None or star.get("dec") is None:
        return [dict(nb) for nb in neighbors]
    ra0, dec0 = propagate_position(star["ra"], star["dec"], star.get("pmRA"), star.get("pmDEC"),
                                   TIC_EPOCH_YEAR, year)
    out = []
    for nb in neighbors:
        entry = {**nb, "dstArcSec_j2000": nb.get("dstArcSec")}
        if nb.get("ra") is not None and nb.get("dec") is not None:
            ra1, dec1 = propagate_position(nb["ra"], nb["dec"], nb.get("pmRA"), nb.get("pmDEC"),
                                           TIC_EPOCH_YEAR, year)
            entry["dstArcSec"] = _separation_arcsec(ra0, dec0, ra1, dec1)
            entry["epoch_year"] = year
        out.append(entry)
    return out


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


def _contact_ratio(k: float, b: float) -> float:
    """T12/T13 for a circular orbit, with T13 = (T14 + T23) / 2; 1 when grazing."""
    t14 = math.sqrt(max((1 + k) ** 2 - b ** 2, 0.0))
    t23 = math.sqrt(max((1 - k) ** 2 - b ** 2, 0.0)) if b < 1 - k else 0.0
    if t14 <= 0:
        return 1.0
    return min(1.0, (t14 - t23) / (t14 + t23))


def shape_delta_mag_limit(fit: dict, depth_ppm: float, n_sigma: float = 2.0) -> dict:
    """Largest magnitude difference a blended eclipsing source could have and
    still produce this transit shape: dm <= 2.5 log10((T12/T13)^2 / depth)
    (Vanderburg et al. 2019, ApJL 881, L19). Rp/R* and b are pushed n_sigma
    toward a longer ingress, which can only raise (loosen) the limit."""
    k, b = fit.get("rp_rs"), fit.get("impact")
    if k is None or b is None or not depth_ppm or depth_ppm <= 0:
        return {"delta_mag_max": None, "reason": "needs fitted Rp/R*, impact, and a positive depth"}
    k_hi = k + n_sigma * abs(fit.get("rp_rs_sigma") or 0.0)
    b_hi = max(b, 0.0) + n_sigma * abs(fit.get("impact_sigma") or 0.0)
    ratio = _contact_ratio(k_hi, b_hi)
    depth = depth_ppm * 1e-6
    limit = 2.5 * math.log10(ratio ** 2 / depth) if ratio ** 2 > depth else 0.0
    return {"delta_mag_max": limit, "t12_over_t13": ratio, "rp_rs_used": k_hi,
            "impact_used": b_hi, "n_sigma": n_sigma,
            "reference": "Vanderburg et al. 2019, ApJL 881, L19"}


def blend_capability(star: dict, neighbors: list[dict], depth_ppm: float,
                     radius_pixels: float = 2.0, delta_mag_max: float | None = None) -> list[dict]:
    """Neighbors within `radius_pixels` TESS pixels that could produce the
    observed depth: bright enough if fully eclipsed (F_nb / (F_target + F_nb)
    >= depth) and, when a shape limit is given, no fainter than delta_mag_max."""
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
        delta = nb_tmag - tmag
        flux_ok = max_depth_ppm >= depth_ppm
        shape_ok = True if delta_mag_max is None else delta <= delta_mag_max
        out.append({"tic": nb.get("ID"), "separation_arcsec": separation, "tmag": nb_tmag,
                    "delta_tmag": delta, "max_eclipse_depth_ppm": max_depth_ppm,
                    "flux_allows_depth": flux_ok, "shape_allows_dilution": shape_ok,
                    "can_produce_depth": flux_ok and shape_ok})
    return sorted(out, key=lambda e: e["separation_arcsec"])
