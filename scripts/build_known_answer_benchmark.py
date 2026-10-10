#!/usr/bin/env python3
"""Build a frozen TESS known-planet benchmark and paired inverted controls.

The default design uses 200 host stars for score calibration and 200 separate
host stars for blind evaluation. Each light curve produces one positive case
and one inverted-flux null control. Catalog and light-curve products are
snapshotted locally with hashes so later runs do not silently change samples.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import hashlib
import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
import requests
from requests.adapters import HTTPAdapter

TAP_URL = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
TAP_QUERY = """select pl_name,hostname,tic_id,disc_facility,tran_flag,
pl_controv_flag,ttv_flag,pl_orbper,pl_orbpererr1,pl_tranmid,pl_tranmiderr1,
pl_trandur,pl_rade,st_rad,st_teff,sy_tmag,pl_tranmid_systemref,
pl_tranmid_reflink,pl_orbper_reflink
from pscomppars
where tran_flag=1 and tic_id is not null and pl_orbper between 0.15 and 30
and pl_tranmid is not null and pl_trandur is not null
and disc_facility like '%Transiting Exoplanet Survey Satellite%'"""


def digest_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


class RequestTimeoutAdapter(HTTPAdapter):
    """Give MAST streaming downloads finite connect and read timeouts."""
    def __init__(self, timeout_seconds: float):
        self.request_timeout = (min(15.0, timeout_seconds), timeout_seconds)
        super().__init__(max_retries=0)

    def send(self, request, **kwargs):
        if kwargs.get("timeout") is None:
            kwargs["timeout"] = self.request_timeout
        return super().send(request, **kwargs)


def stable_key(seed: int, value: str) -> str:
    return hashlib.sha256(f"{seed}:{value}".encode()).hexdigest()


def finite(row: dict, key: str, default: float | None = None) -> float | None:
    try:
        value = float(row.get(key, ""))
        return value if math.isfinite(value) else default
    except (TypeError, ValueError):
        return default


def eligible_planets(rows: list[dict]) -> list[dict]:
    """Keep one published, non-controversial TESS planet per unique TIC."""
    best: dict[str, dict] = {}
    for row in rows:
        tic = re.sub(r"\.0$", "", str(row.get("tic_id", "")).strip())
        tic = re.sub(r"^TIC\s*", "", tic, flags=re.IGNORECASE)
        period = finite(row, "pl_orbper")
        t0 = finite(row, "pl_tranmid")
        duration_hr = finite(row, "pl_trandur")
        if (not tic.isdigit() or not row.get("pl_name") or period is None
                or not 0.15 <= period <= 30 or t0 is None
                or duration_hr is None or duration_hr <= 0):
            continue
        if str(row.get("tran_flag", "1")) not in {"1", "1.0"}:
            continue
        if str(row.get("pl_controv_flag", "0")) in {"1", "1.0"}:
            continue
        # Large TTVs make a linear catalog ephemeris a poor recovery label.
        if str(row.get("ttv_flag", "0")) in {"1", "1.0"}:
            continue
        # PSCompPars can include TOI candidate labels as well as confirmed
        # planet names. A numeric TOI suffix (for example TOI-123.01) is a
        # candidate identifier, not a confirmed planet designation.
        if re.search(r"\.\d+$", str(row["pl_name"]).strip()):
            continue
        period_error = finite(row, "pl_orbpererr1")
        epoch_error = finite(row, "pl_tranmiderr1")
        if period_error is None or epoch_error is None or period_error < 0 or epoch_error < 0:
            continue
        item = dict(row)
        item["tic_id"] = tic
        item["duration_days"] = duration_hr / 24.0
        item["stratum_period"] = ("short" if period < 3 else "medium" if period < 10 else "long")
        radius, star_radius = finite(row, "pl_rade"), finite(row, "st_rad")
        depth_ppm = (1e6 * (radius / (star_radius * 109.076)) ** 2
                     if radius is not None and radius > 0 and star_radius and star_radius > 0
                     else None)
        item["depth_ppm_proxy"] = depth_ppm
        item["stratum_depth"] = ("shallow" if depth_ppm is not None and depth_ppm < 500
                                  else "medium" if depth_ppm is not None and depth_ppm < 5000
                                  else "deep" if depth_ppm is not None else "unknown")
        # A star, not each planet around it, is the sampling unit.
        old = best.get(tic)
        if old is None or stable_key(0, item["pl_name"]) < stable_key(0, old["pl_name"]):
            best[tic] = item
    return sorted(best.values(), key=lambda r: (r["stratum_period"], r["tic_id"]))


def archive_time_to_btjd(value: float) -> float:
    """Convert an absolute Julian transit midpoint to TESS BTJD."""
    return value - 2457000.0 if value > 2_400_000 else value


def qlp_processing_family(products: list[dict]) -> tuple[str, list[int]]:
    """Classify QLP product-era metadata from MAST product names/sequences."""
    sectors = set()
    for product in products:
        for key in ("sector", "sequence_number", "productFilename", "obs_id"):
            raw = product.get(key, "")
            if key in {"sector", "sequence_number"}:
                try:
                    number = int(float(raw))
                    if 1 <= number <= 999:
                        sectors.add(number)
                        break
                except (TypeError, ValueError, OverflowError):
                    pass
            match = re.search(r"(?:s|sector[_-]?)(\d{1,4})", str(raw), re.IGNORECASE)
            if match:
                number = int(match.group(1))
                if 1 <= number <= 999:
                    sectors.add(number)
                    break
    if not sectors:
        return "unknown", []
    families = {"TGLC-based QLP (S94+)" if sector >= 94 else "legacy QLP"
                for sector in sectors}
    family = next(iter(families)) if len(families) == 1 else "mixed-QLP-eras"
    return family, sorted(sectors)


def split_for_host(seed: int, tic_id: str) -> str:
    return "blind" if int(stable_key(seed, tic_id)[:8], 16) % 2 else "calibration"


def select_rows(rows: list[dict], count_per_split: int, seed: int,
                allow_incomplete: bool = False) -> list[dict]:
    """Balance period strata while selecting a stable, host-level sample."""
    groups: dict[str, list[dict]] = {}
    for row in rows:
        groups.setdefault(row["stratum_period"], []).append(row)
    for stratum, group in groups.items():
        groups[stratum] = sorted(group, key=lambda r: stable_key(seed, r["tic_id"]))
    chosen: list[dict] = []
    cursor = {key: 0 for key in groups}
    split_counts = {"blind": 0, "calibration": 0}
    keys = sorted(groups)
    while min(split_counts.values()) < count_per_split:
        advanced = False
        for key in keys:
            group = groups[key]
            while cursor[key] < len(group):
                row = group[cursor[key]]
                cursor[key] += 1
                split = split_for_host(seed, row["tic_id"])
                if split_counts[split] >= count_per_split:
                    continue
                row = dict(row, split=split)
                chosen.append(row)
                split_counts[split] += 1
                advanced = True
                break
        if not advanced:
            break
    if min(split_counts.values()) < count_per_split and not allow_incomplete:
        raise ValueError(
            f"catalog has too few eligible unique hosts for requested design: "
            f"calibration={split_counts['calibration']}, blind={split_counts['blind']}, "
            f"requested={count_per_split} each")
    return chosen


def query_archive(timeout: float) -> tuple[bytes, list[dict]]:
    response = requests.get(TAP_URL, params={"query": TAP_QUERY, "format": "csv"}, timeout=timeout)
    response.raise_for_status()
    payload = response.content
    rows = list(csv.DictReader(payload.decode("utf-8-sig").splitlines()))
    return payload, rows


def export_lightcurve(lc, output: Path, invert: bool = False) -> tuple[int, float, float, float]:
    import numpy as np

    lc = lc.remove_nans().normalize()
    time = np.asarray(lc.time.value, dtype=float)
    flux = np.ma.filled(lc.flux.value, np.nan).astype(float)
    error = np.ma.filled(lc.flux_err.value, np.nan).astype(float)
    good = np.isfinite(time) & np.isfinite(flux) & (flux > 0)
    good &= np.isfinite(error) & (error > 0)
    time, flux, error = time[good], flux[good], error[good]
    order = np.argsort(time, kind="stable")
    time, flux, error = time[order], flux[order], error[order]
    if len(time) < 500:
        raise ValueError(f"only {len(time)} valid observations")
    if invert:
        flux = 2.0 - flux
        if np.any(flux <= 0):
            raise ValueError("inverted flux is non-positive; reject this control")
    data = np.column_stack((time, flux, error))
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(output, data, delimiter=",", header="time,flux,flux_err", comments="", fmt="%.17g")
    cadence = float(np.median(np.diff(time)))
    return len(time), float(time[0]), float(time[-1]), cadence


def expected_events_observed(times: list[float], t0: float, period: float,
                             duration: float) -> int:
    import numpy as np

    observed = np.asarray(times, dtype=float)
    if not len(observed):
        return 0
    first = math.ceil((observed.min() - t0) / period)
    last = math.floor((observed.max() - t0) / period)
    tolerance = max(duration / 2, 2 * float(np.median(np.diff(observed))))
    return sum(bool(np.any(np.abs(observed - (t0 + n * period)) <= tolerance))
               for n in range(first, last + 1))


def fetch_planet_case(planet: dict, seed: int, out: Path, cache: Path,
                      timeout_seconds: float):
    """Download and freeze one host's positive and paired inverted curves."""
    import numpy as np
    import lightkurve as lk
    from astroquery.mast import Observations

    timeout_adapter = RequestTimeoutAdapter(timeout_seconds)
    Observations._session.mount("https://", timeout_adapter)
    Observations._session.mount("http://", timeout_adapter)

    split = planet["split"]
    tic = planet["tic_id"]
    try:
        search = lk.search_lightcurve(f"TIC {tic}", mission="TESS", author="QLP")
        if len(search) == 0:
            raise ValueError("no public TESS QLP light-curve products")
        collection = search.download_all(quality_bitmask="hard", download_dir=str(cache))
        if collection is None or len(collection) == 0:
            raise ValueError("QLP products were listed but none downloaded")
        mast_products = []
        for product in search.table:
            mast_products.append({
                "obs_id": str(product["obs_id"]) if "obs_id" in search.table.colnames else "",
                "productFilename": str(product["productFilename"])
                if "productFilename" in search.table.colnames else "",
                "sequence_number": str(product["sequence_number"])
                if "sequence_number" in search.table.colnames else "",
                "t_min": str(product["t_min"]) if "t_min" in search.table.colnames else "",
                "t_max": str(product["t_max"]) if "t_max" in search.table.colnames else "",
            })
        processing_family, product_sectors = qlp_processing_family(mast_products)
        lightcurve = collection.stitch()
        normalized = lightcurve.remove_nans().normalize()
        times = np.asarray(normalized.time.value, dtype=float)
        period = float(planet["pl_orbper"])
        epoch_btjd = archive_time_to_btjd(float(planet["pl_tranmid"]))
        duration = float(planet["duration_days"])
        period_error = float(planet["pl_orbpererr1"])
        epoch_error = float(planet["pl_tranmiderr1"])
        ephemeris_sigma = math.hypot(epoch_error,
                                     max(0.0, float(np.nanmax(times)) - epoch_btjd)
                                     * period_error / period)
        if 3 * ephemeris_sigma > duration / 2:
            raise ValueError("published linear ephemeris uncertainty exceeds half a transit duration")
        n_events = expected_events_observed(times.tolist(), epoch_btjd, period, duration)
        cadence = float(np.median(np.diff(np.sort(times))))
        if duration < 2 * cadence or n_events < 2:
            raise ValueError(f"known ephemeris has insufficient coverage (events={n_events}, cadence_days={cadence:.6g})")

        planet_case_id = "case-" + stable_key(seed, f"{tic}:positive")[:16]
        control_case_id = "case-" + stable_key(seed, f"{tic}:control")[:16]
        source = out / "lightcurves" / f"{planet_case_id}.csv"
        inverted = out / "lightcurves" / f"{control_case_id}.csv"
        count, tmin, tmax, cadence = export_lightcurve(lightcurve, source)
        export_lightcurve(lightcurve, inverted, invert=True)
        hashes = {"planet": digest_bytes(source.read_bytes()),
                  "inverted": digest_bytes(inverted.read_bytes())}
        metadata = {
            "case_id": "case-" + stable_key(seed, tic)[:16], "host_id": tic, "split": split,
            "mission": "TESS", "period_stratum": planet["stratum_period"],
            "depth_stratum": planet["stratum_depth"], "tic_id": tic,
            "planet_name": planet["pl_name"], "hostname": planet["hostname"],
            "period_days": period, "epoch_days": epoch_btjd, "duration_days": duration,
            "radius_earth": planet.get("pl_rade", ""),
            "depth_ppm_proxy": planet.get("depth_ppm_proxy") or "",
            "period_error_days": planet.get("pl_orbpererr1", ""),
            "epoch_error_days": planet.get("pl_tranmiderr1", ""),
            "stellar_radius_solar": planet.get("st_rad", ""),
            "stellar_teff_k": planet.get("st_teff", ""), "tmag": planet.get("sy_tmag", ""),
            "baseline_days": tmax - tmin, "observed_catalog_events": n_events,
            "cadence_days": cadence, "observations": count,
            "cadence_stratum": "short" if cadence < 0.01 else "long",
            "stellar_type_proxy": ("cool" if finite(planet, "st_teff") and finite(planet, "st_teff") < 4500
                                   else "solar-hot" if finite(planet, "st_teff") else "unknown"),
            "magnitude_stratum": ("bright" if finite(planet, "sy_tmag") and finite(planet, "sy_tmag") < 10
                                  else "faint" if finite(planet, "sy_tmag") else "unknown"),
            "catalog_reference": planet.get("pl_orbper_reflink", ""),
            "epoch_reference": planet.get("pl_tranmid_reflink", ""),
            "period_reference": planet.get("pl_orbper_reflink", ""),
            "time_reference": planet.get("pl_tranmid_systemref", ""),
            "extraction_family": processing_family,
            "product_sectors": json.dumps(product_sectors),
            "mast_products_json": json.dumps(mast_products, sort_keys=True),
            "lightcurve_sha256": hashes["planet"], "inverted_sha256": hashes["inverted"],
        }
        rows = []
        truth = []
        for row_id, data_file, flux_hash, truth_class in (
            (planet_case_id, source, hashes["planet"], "confirmed_planet"),
            (control_case_id, inverted, hashes["inverted"], "inverted_null"),
        ):
            candidate_path = Path("results") / row_id / "candidates.csv"
            # Threshold calibration uses only null controls. Keep calibration
            # positives and their ephemerides in the truth manifest for
            # provenance, but do not search them or expose their labels to the
            # label-free run stage.
            if truth_class == "inverted_null" or split == "blind":
                rows.append({"case_id": row_id, "split": split, "mission": "TESS",
                             "data_csv": str(data_file.relative_to(out)),
                             "candidates_csv": str(candidate_path)})
            truth.append({
                "case_id": row_id, "pair_id": metadata["case_id"], "host_id": tic,
                "split": split, "truth_class": truth_class,
                "candidates_csv": str(candidate_path), "baseline_days": tmax - tmin,
                "period_days": period if truth_class == "confirmed_planet" else "",
                "epoch_days": epoch_btjd if truth_class == "confirmed_planet" else "",
                "duration_days": duration if truth_class == "confirmed_planet" else "",
                "period_error_days": planet.get("pl_orbpererr1", "") if truth_class == "confirmed_planet" else "",
                "epoch_error_days": planet.get("pl_tranmiderr1", "") if truth_class == "confirmed_planet" else "",
                "mission": "TESS", "period_stratum": planet["stratum_period"],
                "depth_stratum": planet["stratum_depth"],
                "depth_ppm_proxy": planet.get("depth_ppm_proxy") or "",
                "cadence_days": cadence, "observations": count,
                "cadence_stratum": "short" if cadence < 0.01 else "long",
                "stellar_teff_k": planet.get("st_teff", ""),
                "magnitude_stratum": metadata["magnitude_stratum"],
                "extraction_family": processing_family,
                "product_sectors": json.dumps(product_sectors),
                "data_sha256": flux_hash, "planet_name": planet["pl_name"],
            })
        return {"split": split, "tic_id": tic, "planet_name": planet["pl_name"],
                "metadata": metadata, "run_rows": rows, "truth_rows": truth,
                "stratum": f"TESS:{planet['stratum_period']}"}
    except Exception as error:
        return {"split": split, "tic_id": tic, "planet_name": planet["pl_name"],
                "failure": f"{type(error).__name__}: {error}"}


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=root / "data/benchmarks/generated")
    parser.add_argument("--cases-per-split", type=int, default=200,
                        help="independent host stars in each calibration/blind split; 200 gives ~±7 pp max 95%% CI")
    parser.add_argument("--seed", type=int, default=20261009)
    parser.add_argument("--timeout", type=float, default=60)
    parser.add_argument("--workers", type=int, default=4,
                        help="concurrent MAST target downloads (reduce to 1 if rate-limited)")
    parser.add_argument("--pilot", action="store_true",
                        help="allow fewer products than requested; output will be explicitly marked underpowered")
    args = parser.parse_args()
    if args.cases_per_split < 1 or args.timeout <= 0 or args.workers < 1:
        parser.error("cases-per-split, timeout, and workers must be positive")

    try:
        import lightkurve as lk
    except ImportError as error:
        parser.error(f"install project Python dependencies first: {error}")
    out = args.output.expanduser().resolve()
    if (out / "provenance.json").exists():
        parser.error("output already contains a benchmark release; choose a new --output directory to preserve it")
    out.mkdir(parents=True, exist_ok=True)
    try:
        catalog_bytes, catalog_rows = query_archive(args.timeout)
    except Exception as error:
        print(f"NASA Exoplanet Archive TAP query failed: {error}", file=sys.stderr)
        return 2
    catalog_path = out / "nasa-pscomppars-tess-transiting.csv"
    catalog_path.write_bytes(catalog_bytes)
    catalog_hash = digest_bytes(catalog_bytes)
    # Fetch an oversampled deterministic pool so broken/missing MAST products
    # do not immediately make the target sample short.
    pool_per_split = math.ceil(args.cases_per_split * 1.25)
    candidates = select_rows(eligible_planets(catalog_rows), pool_per_split,
                             args.seed, allow_incomplete=True)


    truth_rows: list[dict] = []
    run_rows: list[dict] = []
    metadata_rows: list[dict] = []
    failures: list[dict] = []
    split_success = {"calibration": 0, "blind": 0}
    selected_by_stratum: dict[str, int] = {}
    cache = out / "mast-cache"
    cache.mkdir(exist_ok=True)
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for split in ("calibration", "blind"):
            split_candidates = [planet for planet in candidates if planet["split"] == split]
            for start in range(0, len(split_candidates), args.workers):
                if split_success[split] >= args.cases_per_split:
                    break
                batch = split_candidates[start:start + args.workers]
                futures = [executor.submit(fetch_planet_case, planet, args.seed, out, cache,
                                           args.timeout)
                           for planet in batch]
                # Resolve results in stable candidate order so output files and
                # the exact selected sample do not depend on network timing.
                for future in futures:
                    result = future.result()
                    if "failure" in result:
                        failures.append({"tic_id": result["tic_id"],
                                         "planet_name": result["planet_name"],
                                         "split": result["split"], "reason": result["failure"]})
                        print(f"excluded {result['split']}: {result['planet_name']} "
                              f"(TIC {result['tic_id']}): {result['failure']}",
                              file=sys.stderr, flush=True)
                        continue
                    if split_success[split] >= args.cases_per_split:
                        for row in result["run_rows"]:
                            (out / row["data_csv"]).unlink(missing_ok=True)
                        continue
                    metadata_rows.append(result["metadata"])
                    run_rows.extend(result["run_rows"])
                    truth_rows.extend(result["truth_rows"])
                    selected_by_stratum[result["stratum"]] = selected_by_stratum.get(result["stratum"], 0) + 1
                    split_success[split] += 1
                    print(f"[{sum(split_success.values())}/{2 * args.cases_per_split}] {split}: "
                          f"{result['planet_name']} (TIC {result['tic_id']})", flush=True)
    achieved = min(split_success.values())
    if achieved < args.cases_per_split and not args.pilot:
        print(f"under target: calibration={split_success['calibration']}, blind={split_success['blind']}; "
              "rerun with --pilot only to save an explicitly underpowered sample", file=sys.stderr)
    write_csv(out / "run-manifest.csv", run_rows)
    write_csv(out / "truth-manifest.csv", truth_rows)
    write_csv(out / "case-metadata.csv", metadata_rows)
    write_csv(out / "excluded.csv", failures)
    design = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "catalog_source": TAP_URL, "catalog_table": "pscomppars",
        "catalog_doi": "10.26133/NEA13",
        "catalog_query": TAP_QUERY, "catalog_snapshot": catalog_path.name,
        "catalog_sha256": catalog_hash, "seed": args.seed,
        "requested_per_split": args.cases_per_split,
        "catalog_candidate_pool_per_split": pool_per_split,
        "achieved_per_split": split_success,
        "positive_host_count": sum(split_success.values()),
        "paired_inverted_control_count": sum(split_success.values()),
        "period_strata_counts": selected_by_stratum,
        "lightcurve_author": "QLP", "lightkurve_version": lk.__version__,
        "statistical_precision": "At n=200 blind hosts, worst-case 95% Wilson half-width is about 6.9 percentage points. Calibration and blind hosts are disjoint; one planet per TIC; each positive is paired with one inverted control.",
        "pilot": args.pilot,
        "underpowered": args.pilot or achieved < args.cases_per_split,
        "limitations": [
            "Known-planet targets are discovery-selected and do not measure population completeness.",
            "Inverted controls are empirical null controls, not a complete model of astrophysical false positives.",
            "TESS QLP calibrated light curves test light-curve search recovery, not pixel-level validation.",
            "Catalog ephemerides are frozen at build time; TTV-flagged systems and ephemerides with 3-sigma uncertainty over half a transit duration are excluded.",
        ],
    }
    (out / "provenance.json").write_text(json.dumps(design, indent=2) + "\n")
    print(json.dumps(design, indent=2))
    return 0 if achieved >= args.cases_per_split or args.pilot else 3


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("\n")
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    raise SystemExit(main())
