#!/usr/bin/env python3
"""Batch TESS FFI search using Lightkurve and the astrosearch Rust CLI.

Targets CSV: tic_id and TIC/CTL metadata (teff_k, distance_pc, radius_solar,
tmag recommended; stellar_density_g_cm3 optional for the adaptive transit grid).
SPOC and QLP products remain separate for comparison.
"""
from __future__ import annotations
import argparse, atexit, csv, hashlib, json, math, os, random, shutil, signal, subprocess, sys, tempfile, time, uuid
from datetime import datetime, timezone
from pathlib import Path
from survey_ledger import begin_run, config_digest, connect as connect_ledger
from survey_ledger import claim_target, finish_run, record_target, should_defer
import requests
from requests.adapters import HTTPAdapter

ALIASES = {"teff_k": ("teff_k","teff","Teff"), "distance_pc": ("distance_pc","distance","d"),
           "radius_solar": ("radius_solar","rad","radius"), "tmag": ("tmag","Tmag"),
           "mass_solar": ("mass_solar","mass","Mstar","mstar"),
           "stellar_density_g_cm3": ("stellar_density_g_cm3","rho_g_cm3","rho"),
           "cdpp_1h_ppm": ("cdpp_1h_ppm","cdpp1h_ppm","cdpp_ppm"),
           "contamination_ratio": ("contamination_ratio","contratio","contam_ratio")}
_ACTIVE_RUN = None

def finalize_interrupted_run():
    global _ACTIVE_RUN
    if _ACTIVE_RUN is None: return
    ledger,run_id=_ACTIVE_RUN
    try:
        ledger.rollback()
        finish_run(ledger,run_id,"interrupted")
        ledger.close()
    except Exception as error:
        print(f"could not finalize interrupted ledger run: {error}",file=sys.stderr)
    _ACTIVE_RUN=None

atexit.register(finalize_interrupted_run)

class RequestTimeoutAdapter(HTTPAdapter):
    """Apply a finite connect/read timeout to Astroquery's session requests."""
    def __init__(self, timeout_seconds):
        self.request_timeout=timeout_seconds
        super().__init__(max_retries=0)

    def send(self, request, **kwargs):
        if kwargs.get("timeout") is None:
            kwargs["timeout"]=self.request_timeout
        return super().send(request,**kwargs)

def val(row, name):
    for key in ALIASES.get(name, (name,)):
        try:
            x = float(row.get(key, ""))
            if math.isfinite(x) and x > 0: return x
        except (TypeError, ValueError): pass
    return None

def tic(row):
    raw = str(row.get("tic_id", row.get("tic", row.get("ID", "")))).strip()
    if raw.upper().startswith("TIC"):
        raw = raw[3:].strip()
    if raw.endswith(".0") and raw[:-2].isascii() and raw[:-2].isdigit():
        raw = raw[:-2]
    if not raw.isascii() or not raw.isdigit() or int(raw) <= 0:
        raise ValueError("invalid TIC ID")
    return raw

def read_targets(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if not rows: raise ValueError("target list is empty")
    return rows

def accepted(row, a):
    t,d,r,mag = (val(row,k) for k in ("teff_k","distance_pc","radius_solar","tmag"))
    if t is None or not a.min_teff <= t <= a.max_teff: return False, "missing/out-of-range Teff (K/M cut)"
    if a.require_distance and d is None: return False, "missing distance"
    if d is not None and d > a.max_distance_pc: return False, f"distance {d:g} pc exceeds limit"
    if r is not None and r > a.max_radius: return False, f"radius {r:g} R_sun exceeds limit"
    if mag is not None and mag > a.max_tmag: return False, f"Tmag {mag:g} exceeds limit"
    return True, "passes stellar cuts"

def priority(row, nsectors, baseline):
    # Transparent target triage score, never a candidate significance.
    d,r = val(row,"distance_pc"),val(row,"radius_solar")
    score = (.35*min(1.,baseline/2190.) + .25*min(1.,math.log1p(nsectors)/math.log1p(20))
             + .25*max(0.,1.-(d or 100.)/100.) + .15*max(0.,1.-(r or 1.)/1.2))
    return round(100*score,3)

def local_snr_proxy(row):
    """Relative fixed-radius-planet S/N proxy, preferring measured CDPP."""
    radius,mag=val(row,"radius_solar"),val(row,"tmag")
    if radius is None: return None
    cdpp=val(row,"cdpp_1h_ppm")
    if cdpp is not None:
        proxy=1.0/(radius*radius*cdpp)
    elif mag is not None:
        # Photon-limited fallback; not a measured TESS noise level.
        proxy=10.0**(-0.2*(mag-10.0))/(radius*radius)
    else:
        return None
    for key in ALIASES["contamination_ratio"]:
        try:
            contamination=float(row.get(key,""))
            if math.isfinite(contamination) and contamination>=0:
                proxy /= 1.0+contamination
            break
        except (TypeError,ValueError):
            continue
    return proxy

def sector_noise_density(times, flux, exposure_s):
    """Robust adjacent-difference scatter normalized per sqrt(second)."""
    import numpy as np
    if len(times) < 30 or exposure_s <= 0:
        return None
    dt=np.diff(times)
    df=np.diff(flux)
    cadence=np.median(dt[np.isfinite(dt)&(dt>0)])
    if not math.isfinite(cadence) or cadence <= 0:
        return None
    adjacent=np.isfinite(df)&np.isfinite(dt)&(dt<=3.0*cadence)
    if int(adjacent.sum()) < 20:
        return None
    differences=df[adjacent]
    center=float(np.median(differences))
    scatter=1.4826*float(np.median(np.abs(differences-center)))/math.sqrt(2.0)
    return scatter*math.sqrt(exposure_s)

def dependencies():
    try:
        import lightkurve as lk
        import numpy as np
    except ImportError as e:
        raise RuntimeError("Install Python dependencies with uv sync") from e
    return lk,np

def configure_mast_timeout(timeout_seconds):
    from astroquery.mast import Observations
    adapter=RequestTimeoutAdapter(timeout_seconds)
    Observations._session.mount("https://",adapter)

def retry_network(operation, description, attempts=4):
    """Retry transient HTTP/connection failures with bounded exponential backoff."""
    for attempt in range(attempts):
        try:
            return operation()
        except (requests.exceptions.RequestException, TimeoutError) as error:
            response=getattr(error,"response",None)
            status=getattr(response,"status_code",None)
            retryable=(response is None or status==429 or (status is not None and status>=500))
            if not retryable or attempt+1==attempts:
                raise
            delay=min(30.0,2.0**attempt)*random.uniform(0.8,1.2)
            print(f"retrying {description} after {type(error).__name__} in {delay:.1f}s "
                  f"({attempt+1}/{attempts-1})",file=sys.stderr,flush=True)
            time.sleep(delay)

def atomic_write_text(path, text, encoding="utf-8"):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=None
    try:
        with tempfile.NamedTemporaryFile("w",encoding=encoding,dir=path.parent,
                prefix=f".{path.name}.",suffix=".tmp",delete=False) as stream:
            temp=Path(stream.name)
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp,path)
    finally:
        if temp is not None and temp.exists(): temp.unlink()

def atomic_savetxt(path, values, np, **kwargs):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=None
    try:
        with tempfile.NamedTemporaryFile("w",dir=path.parent,prefix=f".{path.name}.",
                suffix=".tmp",delete=False) as stream:
            temp=Path(stream.name)
            np.savetxt(stream,values,**kwargs)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp,path)
    finally:
        if temp is not None and temp.exists(): temp.unlink()

def atomic_write_csv(path, fieldnames, rows):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=None
    try:
        with tempfile.NamedTemporaryFile("w",newline="",encoding="utf-8",dir=path.parent,
                prefix=f".{path.name}.",suffix=".tmp",delete=False) as stream:
            temp=Path(stream.name)
            writer=csv.DictWriter(stream,fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(rows)
            stream.flush(); os.fsync(stream.fileno())
        os.replace(temp,path)
    finally:
        if temp is not None and temp.exists(): temp.unlink()

def sha256_file(path):
    digest=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024),b""):
            digest.update(block)
    return digest.hexdigest()

def binary_path(value):
    if value:
        result=Path(value).expanduser().resolve()
        if not result.is_file(): raise FileNotFoundError(result)
        return result
    root=Path(__file__).resolve().parents[1]
    result=root/"target/release/astrosearch"
    if not result.exists():
        cargo=shutil.which("cargo")
        if not cargo: raise RuntimeError("Run cargo build --release or pass --binary")
        subprocess.run([cargo,"build","--release"],cwd=root,check=True)
    return result

def one_target(row,a,binary):
    lk,np=dependencies()
    ident=tic(row)
    density=val(row,"stellar_density_g_cm3")
    ok,reason=accepted(row,a)
    if not ok: return {"tic_id":ident,"status":"skipped","reason":reason}
    root=a.output.resolve()/f"TIC{ident}"
    root.mkdir(parents=True,exist_ok=True)
    author_runs=[]
    for author in a.authors:
        search=retry_network(
            lambda: lk.search_lightcurve(f"TIC {ident}",mission="TESS",author=author),
            f"MAST search for TIC {ident} ({author})")
        # One product/sector, shortest available FFI exposure; do not merge
        # SPOC and QLP reductions into one detection statistic.
        selected={}
        columns=set(search.table.colnames)
        # MAST's current HLSP result tables expose sector as
        # `sequence_number`; some Lightkurve/product combinations use
        # `sector`. Handle both without aborting an entire target.
        sector_column="sector" if "sector" in columns else (
            "sequence_number" if "sequence_number" in columns else None)
        if sector_column is None:
            author_runs.append({"author":author,"status":"no_sector_metadata"})
            continue
        for i in range(len(search)):
            # Lightkurve 2.6's SearchResult integer indexing incorrectly
            # routes through Astropy column selection; table row access is
            # stable across supported versions.
            item=search.table[i]
            try: sector,exptime=int(item[sector_column]),float(item["exptime"])
            except (KeyError,TypeError,ValueError): continue
            if not math.isfinite(exptime) or exptime<a.min_exptime: continue
            if sector not in selected or exptime<selected[sector]["exptime"]:
                selected[sector]={"index":i,"sector":sector,"exptime":exptime,
                                  "author":str(item["author"]),"mission":str(item["mission"])}
                for column in ("dataURI", "productFilename", "obs_id"):
                    if column in columns:
                        selected[sector][column]=str(item[column])
        recent_sectors=[s for s in selected if s>=a.min_sector]
        if not recent_sectors:
            author_runs.append({"author":author,"status":"no_recent_ffi_products"})
            continue
        out=root/author; out.mkdir(parents=True,exist_ok=True); products=[]
        for sector,meta in sorted(selected.items()):
            lc=retry_network(lambda: search[[meta["index"]]].download(
                quality_bitmask="hard",download_dir=str(a.cache_dir.resolve())),
                f"MAST download for TIC {ident}, sector {sector} ({author})")
            lc=lc.remove_nans().normalize()
            t=np.asarray(lc.time.value,dtype=float)
            f=np.ma.filled(np.asarray(lc.flux.value),np.nan).astype(float)
            e=np.ma.filled(np.asarray(lc.flux_err.value),np.nan).astype(float)
            good=np.isfinite(t)&np.isfinite(f)&np.isfinite(e)&(f>0)&(e>0)
            if int(good.sum())<30: continue
            order=np.argsort(t[good],kind="stable")
            values=np.column_stack((t[good][order],f[good][order],e[good][order],
                np.full(int(good.sum()),meta["exptime"]/86400.0),
                np.full(int(good.sum()),sector,dtype=float)))
            name=f"sector_{sector:03}.csv"; file=out/name
            atomic_savetxt(file,values,np,delimiter=",",header="time,flux,flux_err,exposure_days,sector",
                           comments="",fmt="%.15g")
            noise_density=sector_noise_density(values[:,0],values[:,1],meta["exptime"])
            products.append({**meta,"rows":len(values),"file":name,
                "time_start":float(values[0,0]),"time_end":float(values[-1,0]),
                "noise_density":noise_density,
                "sha256":hashlib.sha256(file.read_bytes()).hexdigest()})
        if not products:
            author_runs.append({"author":author,"status":"no_usable_products"}); continue
        measured=[x["noise_density"] for x in products
                  if x["noise_density"] is not None and x["noise_density"]>0]
        reference=float(np.median(measured)) if measured else None
        for product in products:
            score=product["noise_density"]
            product["noise_ratio"]=(score/reference if score is not None and reference else None)
            product["excluded_noisy"]=(product["noise_ratio"] is not None and
                a.max_sector_noise_ratio>0 and
                product["noise_ratio"]>a.max_sector_noise_ratio)
        used_products=[x for x in products if not x["excluded_noisy"]]
        if not used_products:
            author_runs.append({"author":author,"status":"no_usable_products",
                                "reason":"all products failed the sector quality screen"})
            continue
        arrays=[np.loadtxt(out/x["file"],delimiter=",",skiprows=1) for x in used_products]
        combined=np.concatenate(arrays); combined=combined[np.argsort(combined[:,0],kind="stable")]
        combined_file=out/"combined.csv"
        atomic_savetxt(combined_file,combined,np,delimiter=",",header="time,flux,flux_err,exposure_days,sector",
                       comments="",fmt="%.15g")
        baseline=float(combined[-1,0]-combined[0,0])
        sectors=sorted(x["sector"] for x in used_products)
        qualifies=baseline>=a.min_baseline and any(s>=a.min_sector for s in sectors)
        if not qualifies:
            author_runs.append({"author":author,"status":"insufficient_multiyear_baseline",
                "sectors":sectors,"products":products,
                "excluded_noisy_sectors":[x["sector"] for x in products if x["excluded_noisy"]],
                "baseline_days":baseline,
                "recent_multiyear":False,"priority_score":priority(row,len(sectors),baseline),
                "combined_csv":str(combined_file)})
            continue
        candidates=out/"candidates.csv"
        candidates_stage=out/f".candidates.{uuid.uuid4().hex}.tmp"
        candidates.unlink(missing_ok=True)
        cmd=[str(binary),"--input",str(combined_file),"--output",str(candidates_stage),
             "--min-period",str(a.min_period),"--max-period",str(a.max_period),
             "--threads",str(a.threads),"--trend-method",a.trend_method,
             "--trend-days",str(a.trend_days),"--segment-gap-days",str(a.segment_gap_days),
             "--transit-search-passes",str(a.transit_search_passes),
             "--candidate-grouping",a.candidate_grouping]
        if density is not None: cmd.extend(("--stellar-density",str(density)))
        if a.periods: cmd.extend(("--periods",str(a.periods)))
        log=out/"search.log"
        proc=None
        with tempfile.NamedTemporaryFile("w",dir=out,prefix=".search.",suffix=".tmp",
                                         delete=False) as stream:
            log_tmp=Path(stream.name)
            try:
                proc=subprocess.run(cmd,stdout=stream,stderr=subprocess.STDOUT,
                                    check=False,timeout=a.search_timeout)
            except subprocess.TimeoutExpired:
                stream.write(f"\nsearch timed out after {a.search_timeout}s\n")
            stream.flush(); os.fsync(stream.fileno())
        os.replace(log_tmp,log)
        if proc is None:
            candidates_stage.unlink(missing_ok=True)
            author_runs.append({"author":author,"status":"search_failed",
                "reason":"search timeout","search_timeout_s":a.search_timeout,
                "sectors":sectors,"products":products,"combined_csv":str(combined_file),
                "search_log":str(log)})
            continue
        search_succeeded=proc.returncode==0 and candidates_stage.is_file()
        if search_succeeded:
            os.replace(candidates_stage,candidates)
        else:
            candidates_stage.unlink(missing_ok=True)
        author_runs.append({"author":author,"status":"searched" if search_succeeded else "search_failed",
            "exit_code":proc.returncode,"sectors":sectors,"products":products,
            "excluded_noisy_sectors":[x["sector"] for x in products if x["excluded_noisy"]],
            "baseline_days":baseline,
            "recent_multiyear":qualifies,"priority_score":priority(row,len(sectors),baseline),
            "combined_csv":str(combined_file),"combined_sha256":sha256_file(combined_file),
            "candidates_csv":str(candidates),
            "candidates_sha256":sha256_file(candidates) if search_succeeded else None,
            "search_log":str(log)})
        if search_succeeded and a.tls:
            tls_output=out/f"tls_candidates_{author.lower()}.csv"
            tls_cmd=[sys.executable,str(Path(__file__).with_name("tls_search.py")),
                     "--lightcurve",str(combined_file),"--output",str(tls_output),
                     "--min-period",str(a.min_period),"--max-period",str(a.max_period),
                     "--max-planets",str(a.tls_max_planets),"--min-sde",str(a.tls_min_sde),
                     "--threads",str(a.threads)]
            stellar_radius=val(row,"radius_solar")
            stellar_mass=val(row,"mass_solar")
            if stellar_radius is not None:
                tls_cmd.extend(("--stellar-radius",str(stellar_radius)))
            if stellar_mass is not None:
                tls_cmd.extend(("--stellar-mass",str(stellar_mass)))
            try:
                tls_proc=subprocess.run(tls_cmd,capture_output=True,text=True,
                                        check=False,timeout=a.search_timeout)
                author_runs[-1]["tls_status"]="searched" if tls_proc.returncode==0 else "failed"
                author_runs[-1]["tls_candidates_csv"]=str(tls_output)
                author_runs[-1]["tls_log"]=tls_proc.stdout+tls_proc.stderr
            except subprocess.TimeoutExpired:
                author_runs[-1]["tls_status"]="timeout"
    manifest={"tic_id":ident,"created_utc":datetime.now(timezone.utc).isoformat(),
        "input_catalog_row":row,"authors":author_runs,
        "lightkurve_version":getattr(lk,"__version__","unknown"),
        "search_config":{"min_sector":a.min_sector,"min_exptime_s":a.min_exptime,
            "min_baseline_days":a.min_baseline,"min_period_days":a.min_period,
            "max_period_days":a.max_period,"period_trials":a.periods,"threads":a.threads,
            "search_timeout_s":a.search_timeout,"network_timeout_s":a.network_timeout,
            "trend_method":a.trend_method,
            "trend_days":a.trend_days,"segment_gap_days":a.segment_gap_days,
            "transit_search_passes":a.transit_search_passes,
            "candidate_grouping":a.candidate_grouping,
            "tls_enabled":a.tls,"tls_max_planets":a.tls_max_planets,
            "tls_min_sde":a.tls_min_sde,
            "stellar_radius_solar":val(row,"radius_solar"),
            "stellar_mass_solar":val(row,"mass_solar"),
            "stellar_density_g_cm3":density,
            "max_targets":a.max_targets,
            "max_sector_noise_ratio":a.max_sector_noise_ratio},
        "warning":"Candidates are unvalidated triggers; scores are not false-alarm probabilities."}
    atomic_write_text(root/"manifest.json",json.dumps(manifest,indent=2,default=str)+"\n")
    # A Rust search failure must remain retryable in the shared ledger. Data
    # availability outcomes (for example, no recent products) are completed
    # checks and can use the normal cooldown before being refreshed.
    retryable={"search_failed","no_sector_metadata","no_usable_products"}
    if any(run.get("status") in retryable for run in author_runs):
        status="error"
    elif any(run.get("status")=="searched" for run in author_runs):
        status="processed"
    elif any(run.get("status")=="insufficient_multiyear_baseline" for run in author_runs):
        status="insufficient_multiyear_baseline"
    else:
        status="completed_no_data"
    return {"tic_id":ident,"status":status,"output":str(root),"authors":author_runs}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--targets",required=True,type=Path,help="CSV: tic_id,teff_k,distance_pc,radius_solar,tmag; optional stellar_density_g_cm3")
    p.add_argument("--output",type=Path,default=Path("data/survey"))
    p.add_argument("--binary",help="astrosearch executable; builds release binary by default")
    p.add_argument("--authors",nargs="+",default=["QLP","SPOC"])
    p.add_argument("--min-sector",type=int,default=90)
    p.add_argument("--min-exptime",type=float,default=200,help="Exclude short-cadence products")
    p.add_argument("--max-sector-noise-ratio",type=float,default=5.0,
                   help="Reject sectors above this robust per-cadence-adjusted noise ratio; 0 disables QC")
    p.add_argument("--min-baseline",type=float,default=1095,help="Require about three years and one recent sector")
    p.add_argument("--min-teff",type=float,default=2700); p.add_argument("--max-teff",type=float,default=5200)
    p.add_argument("--max-distance-pc",type=float,default=100)
    p.add_argument("--allow-missing-distance",action="store_true",
                   help="Process targets with unknown distance; default requires a nearby-star measurement")
    p.add_argument("--max-radius",type=float,default=.9); p.add_argument("--max-tmag",type=float,default=14)
    p.add_argument("--min-period",type=float,default=.15); p.add_argument("--max-period",type=float,default=30)
    p.add_argument("--periods",type=int); p.add_argument("--threads",type=int,default=0)
    p.add_argument("--trend-method",choices=("median","biweight"),default="biweight")
    p.add_argument("--trend-days",type=float,default=1.05)
    p.add_argument("--segment-gap-days",type=float,default=.5)
    p.add_argument("--transit-search-passes",type=int,default=5)
    p.add_argument("--candidate-grouping",choices=("ephemeris","events"),default="events",
                   help="group candidates by shared observed events while retaining harmonic alternatives")
    p.add_argument("--tls",action="store_true",help="also run iterative limb-darkened TLS on each searched light curve")
    p.add_argument("--tls-max-planets",type=int,default=5)
    p.add_argument("--tls-min-sde",type=float,default=7.0)
    p.add_argument("--search-timeout",type=float,default=7200,
                   help="Maximum seconds for one Rust search; timed-out targets are retryable")
    p.add_argument("--network-timeout",type=float,default=300,
                   help="Connect/read timeout in seconds for each MAST HTTP request")
    p.add_argument("--pause",type=float,default=1)
    p.add_argument("--max-targets",type=int,
                   help="Optional cap on eligible TICs searched, after priority ordering")
    p.add_argument("--ledger",type=Path,
                   help="Shared SQLite survey ledger (default: project data/survey/survey-ledger.sqlite)")
    p.add_argument("--cache-dir",type=Path,
                   help="Shared MAST download cache (default: project data/survey/mast-cache)")
    p.add_argument("--revisit-after-days",type=int,default=30,
                   help="Reconsider completed targets after this many days; default 30")
    p.add_argument("--revisit",action="store_true",
                   help="Force reconsideration of completed targets regardless of the ledger cooldown")
    a=p.parse_args()
    signal.signal(signal.SIGTERM,lambda signum,frame: sys.exit(128+signum))
    if a.max_targets is not None and a.max_targets < 1:
        p.error("--max-targets must be positive")
    if a.revisit_after_days < 0:
        p.error("--revisit-after-days cannot be negative")
    if not math.isfinite(a.search_timeout) or a.search_timeout <= 0:
        p.error("--search-timeout must be finite and positive")
    if not math.isfinite(a.network_timeout) or a.network_timeout <= 0:
        p.error("--network-timeout must be finite and positive")
    if (not math.isfinite(a.max_sector_noise_ratio) or a.max_sector_noise_ratio < 0
            or 0 < a.max_sector_noise_ratio < 1):
        p.error("--max-sector-noise-ratio must be 0 (disabled) or at least 1")
    a.output.mkdir(parents=True,exist_ok=True)
    a.require_distance = not a.allow_missing_distance
    rows=read_targets(a.targets)
    project_root=Path(__file__).resolve().parents[1]
    ledger_path=(a.ledger.expanduser() if a.ledger else project_root/"data/survey/survey-ledger.sqlite")
    a.cache_dir=(a.cache_dir.expanduser() if a.cache_dir else project_root/"data/survey/mast-cache")
    a.cache_dir.mkdir(parents=True,exist_ok=True)
    ledger=connect_ledger(ledger_path)
    config_sha=config_digest(a)
    run_id=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid.uuid4().hex[:8]
    begin_run(ledger,run_id,a.targets,"local",config_sha)
    global _ACTIVE_RUN
    _ACTIVE_RUN=(ledger,run_id)
    catalog_eligible=[]
    duplicate_rows=set()
    seen_tics=set()
    for index,row in enumerate(rows):
        try:
            ident=tic(row)
            if accepted(row,a)[0]:
                if ident in seen_tics:
                    duplicate_rows.add(index)
                else:
                    seen_tics.add(ident)
                    catalog_eligible.append(index)
        except ValueError:
            continue
    ledger_deferred=set()
    eligible=[]
    for index in catalog_eligible:
        ident=tic(rows[index])
        if should_defer(ledger,ident,rows[index],config_sha,
                        a.revisit_after_days,a.revisit):
            ledger_deferred.add(index)
        else:
            eligible.append(index)
    eligible.sort(key=lambda index:(
        local_snr_proxy(rows[index]) is None,
        -(local_snr_proxy(rows[index]) or 0.0),index))
    selected=eligible[:a.max_targets] if a.max_targets is not None else eligible
    selected_set=set(selected)
    # Record the local physical proxy for every catalog row so selection is
    # auditable and repeatable without an external decision service.
    decision_fields=["input_row","tic_id","eligible","selected","local_snr_proxy","status"]
    decision_rows=[]
    for index,row in enumerate(rows):
        ident=row.get("tic_id",row.get("tic",row.get("ID","")))
        score=local_snr_proxy(row)
        if index in duplicate_rows: state="duplicate_tic_id"
        elif index in ledger_deferred: state="ledger_deferred"
        elif index in selected_set: state="selected"
        elif index in eligible: state="budget_deferred"
        else: state="ineligible"
        decision_rows.append({"input_row":index+1,"tic_id":ident,
            "eligible":index in catalog_eligible,"selected":index in selected_set,
            "local_snr_proxy":"" if score is None else score,"status":state})
    atomic_write_csv(a.output/"local_decisions.csv",decision_fields,decision_rows)
    # Ineligible rows still receive a cheap local skip record; only selected,
    # eligible rows can trigger downloads and Rust searches.
    work_rows=[(i,rows[i]) for i in selected]
    summary=[]
    for index,row in enumerate(rows):
        if index not in selected_set:
            if index in duplicate_rows:
                summary.append({"tic_id":row.get("tic_id",row.get("tic","")),
                    "status":"duplicate_tic_id","reason":"duplicate TIC row; only first eligible row is processed"})
            elif index in ledger_deferred:
                summary.append({"tic_id":row.get("tic_id",row.get("tic","")),
                    "status":"ledger_deferred",
                    "reason":f"already processed within the {a.revisit_after_days}-day ledger cooldown"})
            elif index in eligible:
                summary.append({"tic_id":row.get("tic_id",""),"status":"budget_deferred",
                                "reason":"eligible target was below --max-targets cutoff"})
            else:
                try:
                    ident=tic(row); ok,reason=accepted(row,a)
                    summary.append({"tic_id":ident,"status":"skipped","reason":reason})
                except ValueError as error:
                    summary.append({"tic_id":row.get("tic_id",""),"status":"error","reason":str(error)})
    summary_by_index={}
    binary=binary_path(a.binary)
    if work_rows:
        configure_mast_timeout(a.network_timeout)
    failures=0
    for i,(index,row) in enumerate(work_rows,1):
        ident=tic(row)
        if not claim_target(ledger,run_id,ident,row,config_sha):
            result={"tic_id":ident,"status":"busy",
                    "reason":"another run holds a live target lease"}
            failures+=1
            summary_by_index[index]=result
            print(f"[{i}/{len(work_rows)}] TIC {ident}: busy (not claimed)",flush=True)
            continue
        try: result=one_target(row,a,binary)
        except Exception as e:
            result={"tic_id":row.get("tic_id",""),"status":"error","reason":f"{type(e).__name__}: {e}"}
        if result.get("status")=="error": failures+=1
        summary_by_index[index]=result
        if not record_target(ledger,run_id,ident,row,config_sha,result):
            # Our lease expired and another run took the target over; keep its state.
            result={**result,"status":"lost_claim",
                    "reason":"lease expired and another run recovered this target; result not recorded"}
            failures+=1
            summary_by_index[index]=result
        print(f"[{i}/{len(work_rows)}] TIC {result.get('tic_id')}: {result['status']}",flush=True)
        time.sleep(max(0.,a.pause))
    for index,result in summary_by_index.items():
        summary.append(result)
    # Preserve input order in the summary while recording the actual selection.
    index_for_result={}
    for index,row in enumerate(rows):
        try: index_for_result.setdefault(tic(row),[]).append(index)
        except ValueError: pass
    summary.sort(key=lambda item:min(index_for_result.get(str(item.get("tic_id","")),[len(rows)])))
    atomic_write_text(a.output/"batch_summary.json",json.dumps(summary,indent=2,default=str)+"\n")
    ranking=[]
    for result in summary:
        for run in result.get("authors",[]):
            if run.get("status") in ("searched", "insufficient_multiyear_baseline"):
                ranking.append({"tic_id":result["tic_id"],"author":run["author"],
                    "recent_multiyear":run.get("recent_multiyear",False),"priority_score":run["priority_score"],
                    "sectors":";".join(map(str,run["sectors"])),"baseline_days":run["baseline_days"],
                    "candidates_csv":run.get("candidates_csv","")})
    ranking.sort(key=lambda x:(not x["recent_multiyear"],-x["priority_score"]))
    atomic_write_csv(a.output/"target_ranking.csv",["tic_id","author","recent_multiyear","priority_score",
            "sectors","baseline_days","candidates_csv"],ranking)
    run_status="failed" if failures and failures==len(work_rows) else (
        "partial" if failures else "complete")
    atomic_write_text(a.output/"run_metadata.json",json.dumps({
        "run_id":run_id,"ledger":str(ledger_path.resolve()),"ranking":"local_snr_proxy",
        "config_sha256":config_sha,"targets_csv_sha256":sha256_file(a.targets),
        "astrosearch_binary":str(binary),"binary_sha256":sha256_file(binary),
        "revisit_after_days":a.revisit_after_days,"forced_revisit":a.revisit,
        "catalog_eligible":len(catalog_eligible),"ledger_deferred":len(ledger_deferred),
        "selected":len(selected_set),"failed_targets":failures,"status":run_status,
        "finished_utc":datetime.now(timezone.utc).isoformat()
    },indent=2)+"\n")
    finish_run(ledger,run_id,run_status)
    ledger.close()
    _ACTIVE_RUN=None
    if failures:
        raise SystemExit(1)

if __name__=="__main__":
    try: main()
    except KeyboardInterrupt: sys.exit(130)
