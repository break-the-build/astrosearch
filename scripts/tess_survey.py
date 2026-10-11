#!/usr/bin/env python3
"""Batch TESS FFI search using Lightkurve and the astrosearch Rust CLI.

Targets CSV: tic_id and TIC/CTL metadata (teff_k, distance_pc, radius_solar,
tmag recommended; stellar_density_g_cm3 optional for the adaptive transit grid).
SPOC and QLP products remain separate for comparison.
"""
from __future__ import annotations
import argparse, atexit, csv, hashlib, json, math, os, random, shutil, signal, subprocess, sys, tempfile, time, uuid
import inspect
import multiprocessing
import queue
from contextlib import contextmanager
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, wait, FIRST_COMPLETED
from datetime import datetime, timezone
from pathlib import Path
from survey_ledger import begin_run, config_digest, connect as connect_ledger
from survey_ledger import claim_target, finish_run, record_target, should_defer
from vetting.network import install_default_timeouts
from run_metrics import build_identity as capture_build_identity, run_with_metrics
from candidate_registry import ingest as register_candidates
import target_yield
from search_contract import record_search, calibration_match
import requests
from requests.adapters import HTTPAdapter

ALIASES = {"teff_k": ("teff_k","teff","Teff"), "distance_pc": ("distance_pc","distance","d"),
           "radius_solar": ("radius_solar","rad","radius"), "tmag": ("tmag","Tmag"),
           "mass_solar": ("mass_solar","mass","Mstar","mstar"),
           "stellar_density_g_cm3": ("stellar_density_g_cm3","rho_g_cm3","rho"),
           "cdpp_1h_ppm": ("cdpp_1h_ppm","cdpp1h_ppm","cdpp_ppm"),
           "contamination_ratio": ("contamination_ratio","contratio","contam_ratio")}
_ACTIVE_RUN = None
# Named search configurations. "fast" is grid-study variant J: 10-minute
# binning, physical duration limits, a 7x coarse grid ranked by fast folding
# and refined around the best cells, and a 2M period-trial cap. On 158
# known-answer cases it used 7x less CPU than "standard" and recovered every
# planet "standard" did (plus 3-4 more) at matched false-alarm rates.
SEARCH_PROFILES = {
    "standard": (),
    "fast": ("--bin-minutes", "10", "--physical-duration-floor", "--two-stage", "7",
             "--coarse-ffa", "--max-period-trials", "2000000"),
}

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

def select_lightcurve_products(search, min_exptime):
    """Pick the shortest cadence that meets the survey limit for each sector."""
    selected={}
    columns=set(search.table.colnames)
    sector_column="sector" if "sector" in columns else (
        "sequence_number" if "sequence_number" in columns else None)
    if sector_column is None:
        return selected, False
    for i in range(len(search)):
        item=search.table[i]
        try: sector,exptime=int(item[sector_column]),float(item["exptime"])
        except (KeyError,TypeError,ValueError): continue
        if not math.isfinite(exptime) or exptime<min_exptime: continue
        if sector not in selected or exptime<selected[sector]["exptime"]:
            selected[sector]={"index":i,"sector":sector,"exptime":exptime,
                              "author":str(item["author"]),"mission":str(item["mission"])}
            for column in ("dataURI", "productFilename", "obs_id"):
                if column in columns:
                    selected[sector][column]=str(item[column])
            for column in ("t_min", "t_max"):
                try:
                    value=float(item[column])
                    if math.isfinite(value): selected[sector][column]=value
                except (KeyError, TypeError, ValueError):
                    pass
            if "cached_path" in columns:
                selected[sector]["cached_path"] = str(item["cached_path"])
    return selected, True

def cached_product_search(lk, ident, author, cache_dir, min_exptime):
    """Build a product listing from already-downloaded FITS files, without MAST."""
    from astropy.io import fits
    from astropy.table import Table
    tic_token=f"{int(ident):016d}"
    rows=[]
    for path in Path(cache_dir).expanduser().rglob("*.fits"):
        name=path.name
        if tic_token not in name:
            continue
        is_qlp=name.lower().startswith("hlsp_qlp_tess_ffi_")
        product_author="QLP" if is_qlp else ("SPOC" if "_s_lc.fits" in name.lower() else "")
        if product_author.upper()!=author.upper():
            continue
        try:
            with fits.open(path,memmap=True) as hdus:
                primary=hdus[0].header
                table=next(h.header for h in hdus[1:] if h.header.get("TIMEDEL") is not None)
                sector=int(primary.get("SECTOR",table.get("SECTOR")))
                # Round to the nominal cadence so floating-point FITS values
                # such as 199.99999999999997 s still meet a 200 s cut.
                exptime=round(float(table["TIMEDEL"])*86400.,3)
                t0=float(table["TSTART"]); t1=float(table["TSTOP"])
                tic_header=str(primary.get("TICID",primary.get("OBJECT","")))
                if str(int(ident)) not in tic_header:
                    continue
        except (OSError,ValueError,TypeError,StopIteration):
            continue
        if not math.isfinite(exptime) or exptime<min_exptime or not t0<t1:
            continue
        rows.append({"sector":sector,"exptime":exptime,"author":author,
                     "mission":f"TESS Sector {sector}","dataURI":"local-cache",
                     "productFilename":name,"obs_id":path.stem,"t_min":t0,"t_max":t1,
                     "cached_path":str(path.resolve())})
    if not rows:
        return None
    # Reuse the ordinary per-sector cadence selector and preflight logic.
    table=Table(rows=rows)
    class CachedSearch:
        def __init__(self,table): self.table=table
        def __len__(self): return len(self.table)
    return CachedSearch(table)

def fresh_product_search(lk, ident, author):
    """Bypass Lightkurve's in-process memoization, without clearing shared caches.

    Astroquery's MAST portal requests already use cache=False. Downloaded FITS
    remain reusable. Unwrap only this callable, so other queries are unaffected.
    """
    return retry_network(lambda: inspect.unwrap(lk.search_lightcurve)(
        f"TIC {ident}", mission="TESS", author=author),
        f"fresh MAST product search for TIC {ident} ({author})")

def product_preflight(selected, min_sector, min_baseline):
    """Reject only impossible coverage using complete selected-product envelopes.

    MAST t_min/t_max are MJD envelopes, not quality-screened exposure windows.
    Incomplete envelopes cannot establish an upper bound; download in that case.
    """
    result={"queried_utc":datetime.now(timezone.utc).isoformat(),
            "source":"fresh author-specific MAST product listing; MJD envelopes",
            "sectors":sorted(selected), "products":list(selected.values()),
            "status":"download_required"}
    if not any(sector>=min_sector for sector in selected):
        result["status"]="no_recent_ffi_products"
        return result
    windows=[]
    for product in selected.values():
        t0,t1=product.get("t_min"),product.get("t_max")
        if t0 is None or t1 is None or not t0<t1:
            result["reason"]="incomplete product envelopes; baseline checked after download"
            return result
        windows.append((t0,t1))
    upper=max(t1 for t0,t1 in windows)-min(t0 for t0,t1 in windows)
    result["baseline_upper_bound_days"]=upper
    if upper+1e-6<min_baseline:
        result.update(status="insufficient_product_baseline",
                      reason=f"product envelope {upper:.1f} d is below {min_baseline:g} d")
    return result

def initialize_target_worker(timeout, ready):
    # Separate HTTP sessions and process groups; workers never share a ledger
    # connection. A whole group can be stopped, including its Rust/TLS child.
    if os.name=="posix": os.setsid()
    ready.put(os.getpid())
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    configure_mast_timeout(timeout)
    install_default_timeouts(timeout)

@contextmanager
def target_pool(workers, timeout):
    context=multiprocessing.get_context("spawn")
    ready=context.Queue()
    pool=ProcessPoolExecutor(max_workers=workers,mp_context=context,
        initializer=initialize_target_worker,initargs=(timeout,ready))
    try:
        yield pool
    except BaseException:
        # Initializers publish their group before doing network/import work.
        # Cancel unsent futures and terminate active target groups on interrupt.
        pids=[]
        while True:
            try: pids.append(ready.get_nowait())
            except queue.Empty: break
        for pid in pids:
            try:
                if os.name=="posix": os.killpg(pid,signal.SIGTERM)
                else: os.kill(pid,signal.SIGTERM)
            except ProcessLookupError: pass
        pool.shutdown(wait=True,cancel_futures=True)
        raise
    else:
        pool.shutdown(wait=True)
    finally:
        ready.close()
        ready.join_thread()

def target_results(work_rows, a, binary, claim, pool=None, worker=None):
    """Bound claims and queued work to the worker count; commit in the caller.

    Completion order may differ; selection order and final summary order do not.
    The pause spaces submissions globally rather than multiplying MAST traffic.
    """
    worker=worker or one_target
    def execute(row):
        try: return worker(row,a,binary)
        except Exception as error:
            return {"tic_id":tic(row),"status":"error",
                    "reason":f"{type(error).__name__}: {error}"}
    if pool is None:
        for index,row in work_rows:
            if claim(row): result=execute(row)
            else: result={"tic_id":tic(row),"status":"busy",
                          "reason":"another run holds a live target lease"}
            yield index,row,result
            time.sleep(a.pause)
        return
    pending={}
    rows=iter(work_rows)
    exhausted=False
    next_start=0.0
    while pending or not exhausted:
        while len(pending)<a.target_workers and not exhausted:
            try: index,row=next(rows)
            except StopIteration:
                exhausted=True
                break
            if not claim(row):
                yield index,row,{"tic_id":tic(row),"status":"busy",
                                 "reason":"another run holds a live target lease"}
                continue
            time.sleep(max(0.,next_start-time.monotonic()))
            future=pool.submit(worker,row,a,binary)
            pending[future]=(index,row)
            next_start=time.monotonic()+a.pause
        if not pending: continue
        done,_=wait(pending,return_when=FIRST_COMPLETED)
        for future in sorted(done,key=lambda f:pending[f][0]):
            index,row=pending.pop(future)
            try: result=future.result()
            except Exception as error:
                result={"tic_id":tic(row),"status":"error",
                        "reason":f"{type(error).__name__}: {error}"}
            yield index,row,result

def listed_sectors(row, author):
    raw=str(row.get(f"{author.lower()}_sectors", "") or "")
    sectors=[]
    for token in raw.replace(",", ";").split(";"):
        try: sectors.append(int(token))
        except ValueError: continue
    return sectors

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
    download_workers=max(1,int(getattr(a,"download_workers",2)))
    density=val(row,"stellar_density_g_cm3")
    ok,reason=accepted(row,a)
    if not ok: return {"tic_id":ident,"status":"skipped","reason":reason}
    root=a.output.resolve()/f"TIC{ident}"
    root.mkdir(parents=True,exist_ok=True)
    author_runs=[]
    preflight_by_author={}
    for author in a.authors:
        # One product/sector, shortest available FFI exposure; do not merge
        # SPOC and QLP reductions into one detection statistic.
        min_exptime=(a.min_exptime_spoc if author.upper()=="SPOC" else a.min_exptime)
        if a.offline_cache_only:
            search=cached_product_search(lk,ident,author,a.cache_dir,min_exptime)
            if search is None:
                author_runs.append({"author":author,"status":"no_cached_products",
                    "reason":"no matching usable FITS products were present in the local MAST cache"})
                continue
        else:
            search=fresh_product_search(lk,ident,author)
        selected,has_sector_metadata=select_lightcurve_products(search,min_exptime)
        if not has_sector_metadata:
            author_runs.append({"author":author,"status":"no_sector_metadata"})
            continue
        recent_sectors=[s for s in selected if s>=a.min_sector]
        expected_recent=[s for s in listed_sectors(row,author) if s>=a.min_sector]
        # MAST's indexed observation list can lag the target-yield coverage
        # prefetch. Refresh once before dropping a target that the prefetch
        # says has recent sectors; this avoids a transient empty/stale listing
        # consuming one of the user's --max-targets slots.
        for refresh in range(0 if a.offline_cache_only else 2):
            if recent_sectors or not expected_recent:
                break
            time.sleep(2.0 * (refresh + 1))
            search=fresh_product_search(lk,ident,author)
            selected,has_sector_metadata=select_lightcurve_products(search,min_exptime)
            if not has_sector_metadata:
                break
            recent_sectors=[s for s in selected if s>=a.min_sector]
        if recent_sectors and expected_recent and not set(expected_recent).intersection(recent_sectors):
            print(f"TIC {ident} {author}: refreshed recent-sector set differs from prefetch; "
                  f"prefetch={expected_recent} search={recent_sectors}",file=sys.stderr,flush=True)
        preflight=product_preflight(selected,a.min_sector,a.min_baseline)
        if a.offline_cache_only:
            preflight["source"]="local MAST cache FITS headers; no network query or download"
        preflight_by_author[author]=preflight
        if not recent_sectors:
            author_runs.append({"author":author,"status":"no_recent_ffi_products",
                "product_preflight":preflight,
                "expected_recent_sectors":expected_recent,
                "reason":"MAST returned no usable product at or after the configured minimum sector"})
            continue
        if preflight["status"]=="insufficient_product_baseline":
            author_runs.append({"author":author,"status":preflight["status"],
                "reason":preflight["reason"],"product_preflight":preflight})
            continue
        out=root/author; out.mkdir(parents=True,exist_ok=True)

        def download_sector(item):
            sector,meta=item
            if a.offline_cache_only:
                lc=lk.read(meta["cached_path"],quality_bitmask="hard")
            else:
                lc=retry_network(lambda: search[[meta["index"]]].download(
                    quality_bitmask="hard",download_dir=str(a.cache_dir.resolve())),
                    f"MAST download for TIC {ident}, sector {sector} ({author})")
            flux_column=str(getattr(lc,"flux_col","") or "").upper()
            lc=lc.remove_nans().normalize()
            t=np.asarray(lc.time.value,dtype=float)
            f=np.ma.filled(np.asarray(lc.flux.value),np.nan).astype(float)
            e=np.ma.filled(np.asarray(lc.flux_err.value),np.nan).astype(float)
            good=np.isfinite(t)&np.isfinite(f)&np.isfinite(e)&(f>0)&(e>0)
            if int(good.sum())<30: return None
            order=np.argsort(t[good],kind="stable")
            values=np.column_stack((t[good][order],f[good][order],e[good][order],
                np.full(int(good.sum()),meta["exptime"]/86400.0),
                np.full(int(good.sum()),sector,dtype=float)))
            name=f"sector_{sector:03}.csv"; file=out/name
            atomic_savetxt(file,values,np,delimiter=",",header="time,flux,flux_err,exposure_days,sector",
                           comments="",fmt="%.15g")
            noise_density=sector_noise_density(values[:,0],values[:,1],meta["exptime"])
            return {**meta,"rows":len(values),"file":name,
                "flux_column":flux_column,
                "crowding_corrected":"PDCSAP" in flux_column,
                "time_start":float(values[0,0]),"time_end":float(values[-1,0]),
                "noise_density":noise_density,
                "sha256":hashlib.sha256(file.read_bytes()).hexdigest()}

        def download_stage(items):
            if not items:
                return []
            with ThreadPoolExecutor(max_workers=min(download_workers,len(items)),
                                    thread_name_prefix="mast-sector") as executor:
                return [product for product in executor.map(download_sector,items)
                        if product is not None]

        # Fetch and validate recent sectors first. Only if at least one is
        # usable do we fetch the older baseline; downloads within each stage
        # overlap bounded MAST I/O while preserving the existing early stop.
        recent_items=sorted(((s,m) for s,m in selected.items() if s>=a.min_sector),
                            key=lambda item:item[0])
        older_items=sorted(((s,m) for s,m in selected.items() if s<a.min_sector),
                           key=lambda item:item[0])
        download_started=time.perf_counter()
        products=download_stage(recent_items)
        recent_products_usable=bool(products)
        if recent_products_usable:
            products.extend(download_stage(older_items))
        download_wall_seconds=time.perf_counter()-download_started
        download_metrics={"workers":min(download_workers,max(1,len(recent_items),
            len(older_items))),"wall_seconds":download_wall_seconds,
            "recent_products_requested":len(recent_items),
            "older_products_requested":len(older_items) if recent_products_usable else 0}
        if not products:
            older=[s for s in selected if s<a.min_sector]
            author_runs.append({"author":author,"status":(
                "missing_recent_sector_after_download" if older else "no_usable_products"),
                "reason":"no recent product has at least 30 usable rows",
                "product_preflight":preflight,"older_downloads_avoided":len(older),
                "download_metrics":download_metrics})
            continue
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
                                "reason":"all products failed the sector quality screen",
                                "download_metrics":download_metrics})
            continue
        arrays=[np.loadtxt(out/x["file"],delimiter=",",skiprows=1) for x in used_products]
        combined=np.concatenate(arrays); combined=combined[np.argsort(combined[:,0],kind="stable")]
        combined_file=out/"combined.csv"
        atomic_savetxt(combined_file,combined,np,delimiter=",",header="time,flux,flux_err,exposure_days,sector",
                       comments="",fmt="%.15g")
        combined_parts=[]
        for product in used_products:
            part={**product,
                  "time_min":product["time_start"],"time_max":product["time_end"],
                  "crowding_corrected":product.get("crowding_corrected", True)}
            combined_parts.append(part)
        atomic_write_text(out/"combined.meta.json",json.dumps(
            {"tic_id":ident,"author":author,"rows":int(len(combined)),
             "parts":combined_parts},indent=2,sort_keys=True)+"\n")
        baseline=float(combined[-1,0]-combined[0,0])
        sectors=sorted(x["sector"] for x in used_products)
        has_recent=any(s>=a.min_sector for s in sectors)
        has_baseline=baseline>=a.min_baseline
        qualifies=has_baseline and has_recent
        if not qualifies:
            status="insufficient_multiyear_baseline" if not has_baseline else "missing_recent_sector_after_download"
            reason=(f"usable baseline {baseline:.1f} d is below {a.min_baseline:g} d" if not has_baseline else
                    f"no downloaded usable sector is >= {a.min_sector}")
            author_runs.append({"author":author,"status":status,"reason":reason,
                "sectors":sectors,"products":products,
                "excluded_noisy_sectors":[x["sector"] for x in products if x["excluded_noisy"]],
                "baseline_days":baseline,
                "recent_multiyear":has_baseline and has_recent,
                "expected_recent_sectors":expected_recent,
                "priority_score":priority(row,len(sectors),baseline),
                "combined_csv":str(combined_file),"download_metrics":download_metrics})
            continue
        candidates=out/"candidates.csv"
        candidates_stage=out/f".candidates.{uuid.uuid4().hex}.tmp"
        candidates.unlink(missing_ok=True)
        cmd=[str(binary),"--input",str(combined_file),"--output",str(candidates_stage),
             "--search-profile",a.search_profile,
             "--min-period",str(a.min_period),"--max-period",str(a.max_period),
             "--threads",str(a.threads),"--trend-method",a.trend_method,
             "--trend-days",str(a.trend_days),"--segment-gap-days",str(a.segment_gap_days),
             "--transit-search-passes",str(a.transit_search_passes),
             "--candidate-grouping",a.candidate_grouping]
        if density is not None: cmd.extend(("--stellar-density",str(density)))
        if a.periods: cmd.extend(("--periods",str(a.periods)))
        if a.detection_threshold is not None:
            cmd.extend(("--detection-threshold",str(a.detection_threshold)))
        if a.stop_passes_below_detection:
            cmd.append("--stop-passes-below-detection")
        cmd.extend(SEARCH_PROFILES[a.search_profile])
        log=out/"search.log"
        with tempfile.NamedTemporaryFile("w",dir=out,prefix=".search.",suffix=".tmp",
                                         delete=False) as stream:
            log_tmp=Path(stream.name)
        search_metrics=run_with_metrics(cmd,log_tmp,timeout=a.search_timeout)
        if search_metrics["timed_out"]:
            with log_tmp.open("a",encoding="utf-8") as stream:
                stream.write(f"\nsearch timed out after {a.search_timeout}s\n")
        with log_tmp.open("rb") as stream:
            os.fsync(stream.fileno())
        os.replace(log_tmp,log)
        if search_metrics["timed_out"]:
            candidates_stage.unlink(missing_ok=True)
            author_runs.append({"author":author,"status":"search_failed",
                "reason":"search timeout","search_timeout_s":a.search_timeout,
                "sectors":sectors,"products":products,"combined_csv":str(combined_file),
                "search_log":str(log),"search_metrics":search_metrics})
            continue
        search_succeeded=search_metrics["exit_code"]==0 and candidates_stage.is_file()
        if search_succeeded:
            provenance_stage = Path(str(candidates_stage) + ".provenance.json")
            provenance_final = Path(str(candidates) + ".provenance.json")
            provenance_final.unlink(missing_ok=True)
            if provenance_stage.exists():
                os.replace(provenance_stage, provenance_final)
            os.replace(candidates_stage,candidates)
        else:
            candidates_stage.unlink(missing_ok=True)
        author_runs.append({"author":author,"status":"searched" if search_succeeded else "search_failed",
            "exit_code":search_metrics["exit_code"],"search_metrics":search_metrics,
            "download_metrics":download_metrics,
            "sectors":sectors,"products":products,
            "excluded_noisy_sectors":[x["sector"] for x in products if x["excluded_noisy"]],
            "baseline_days":baseline,
            "recent_multiyear":qualifies,"priority_score":priority(row,len(sectors),baseline),
            "combined_csv":str(combined_file),"combined_sha256":sha256_file(combined_file),
            "candidates_csv":str(candidates),
            "candidates_sha256":sha256_file(candidates) if search_succeeded else None,
            "search_log":str(log)})
        if search_succeeded:
            try:
                search_record = record_search(candidates, binary, author)
                search_record["status"] = "OK"
                calibration_path = a.detection_calibration or ((a.calibration_dir or
                    Path(__file__).resolve().parents[1] / "data/calibration") / "detection_threshold.json")
                calibration = json.loads(calibration_path.read_text()) if calibration_path.exists() else None
                author_runs[-1]["calibration"] = calibration_match(search_record, calibration)
                author_runs[-1]["search_contract_id"] = search_record["contract_id"]
            except (OSError, ValueError, KeyError) as error:
                author_runs[-1]["calibration"] = {"status":"NOT_EVALUATED", "reason":str(error)}
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
        "input_catalog_row":row,"authors":author_runs,"product_preflights":preflight_by_author,
        "source_revision":a.source_revision,"build_identity":a.build_identity,
        "lightkurve_version":getattr(lk,"__version__","unknown"),
        "search_config":{"min_sector":a.min_sector,"min_exptime_s":a.min_exptime,
            "min_baseline_days":a.min_baseline,"min_period_days":a.min_period,
            "max_period_days":a.max_period,"period_trials":a.periods,"threads":a.threads,
            "search_timeout_s":a.search_timeout,"network_timeout_s":a.network_timeout,
            "trend_method":a.trend_method,
            "trend_days":a.trend_days,"segment_gap_days":a.segment_gap_days,
            "transit_search_passes":a.transit_search_passes,
            "detection_threshold":a.detection_threshold,
            "stop_passes_below_detection":a.stop_passes_below_detection,
            "search_profile":a.search_profile,
            "offline_cache_only":a.offline_cache_only,
            "candidate_grouping":a.candidate_grouping,
            "tls_enabled":a.tls,"tls_max_planets":a.tls_max_planets,
            "tls_min_sde":a.tls_min_sde,
            "stellar_radius_solar":val(row,"radius_solar"),
            "stellar_mass_solar":val(row,"mass_solar"),
            "stellar_density_g_cm3":density,
            "max_targets":a.max_targets,
            "download_workers":download_workers,
            "max_sector_noise_ratio":a.max_sector_noise_ratio},
        "warning":"Candidates are unvalidated triggers; scores are not false-alarm probabilities."}
    atomic_write_text(root/"manifest.json",json.dumps(manifest,indent=2,default=str)+"\n")
    registry_path=(a.candidate_registry.expanduser().resolve() if a.candidate_registry
                   else a.output.resolve()/"candidate_registry.sqlite")
    registry_records=[]
    for run in author_runs:
        if run.get("status") != "searched":
            continue
        try:
            registered=register_candidates(
                registry_path, f"TIC{ident}:{run['author']}",
                Path(run["combined_csv"]), Path(run["candidates_csv"]),
                a.source_revision, root/"manifest.json")
            registry_records.append({"author":run["author"],**registered})
        except Exception as error:
            registry_records.append({"author":run["author"],
                                     "error":f"{type(error).__name__}: {error}"})
    if registry_records:
        atomic_write_text(root/"candidate_registry_ingest.json",
                          json.dumps(registry_records,indent=2)+"\n")
    # A Rust search failure must remain retryable in the shared ledger. Data
    # availability outcomes (for example, no recent products) are completed
    # checks and can use the normal cooldown before being refreshed.
    retryable={"search_failed","no_sector_metadata","no_usable_products"}
    has_searched=any(run.get("status")=="searched" for run in author_runs)
    has_retryable=any(run.get("status") in retryable for run in author_runs)
    if has_searched:
        status="partial" if has_retryable else "processed"
    elif has_retryable:
        status="error"
    elif any(run.get("status") in ("insufficient_multiyear_baseline",
                                  "insufficient_product_baseline") for run in author_runs):
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
    p.add_argument("--min-exptime",type=float,default=200,
                   help="Minimum QLP/FFI exposure in seconds (default: 200)")
    p.add_argument("--min-exptime-spoc",type=float,default=120,
                   help="Minimum SPOC exposure in seconds; 120 includes standard 2-minute SPOC products")
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
    p.add_argument("--target-workers",type=int,default=1,
                   help="Concurrent target processes (1-4; default serial). Requires explicit --threads when >1")
    p.add_argument("--download-workers",type=int,default=2,
                   help="Concurrent sector downloads within a target (1-4; total active MAST downloads are capped at 4)")
    p.add_argument("--cpu-budget",type=int,default=os.cpu_count() or 1,
                   help="Maximum target-workers times per-search threads; reserve capacity for other runs")
    p.add_argument("--detection-calibration", type=Path,
                   help="calibration to verify against the resolved search contract; Profile J uses the five-pass J calibration by default")
    p.add_argument("--trend-method",choices=("median","biweight"),default="biweight")
    p.add_argument("--trend-days",type=float,default=1.05)
    p.add_argument("--segment-gap-days",type=float,default=.5)
    p.add_argument("--transit-search-passes",type=int,
                   help="number of masked search passes (defaults to 5 for Profile J and 1 for standard)")
    p.add_argument("--detection-threshold",type=float,
                   help="profile-calibrated jackknife threshold passed to astrosearch (defaults to the selected calibration's deployment threshold)")
    p.add_argument("--stop-passes-below-detection",action="store_true",
                   help="stop iterative searches when the current best score is below the detection tier")
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
    p.add_argument("--candidate-registry",type=Path,
                   help="Candidate database (default: <output>/candidate_registry.sqlite)")
    p.add_argument("--cache-dir",type=Path,
                   help="Shared MAST download cache (default: project data/survey/mast-cache)")
    p.add_argument("--offline-cache-only",action="store_true",
                   help="Search only matching FITS already in --cache-dir; never query or download from MAST")
    p.add_argument("--revisit-after-days",type=int,default=30,
                   help="Reconsider completed targets after this many days; default 30")
    p.add_argument("--revisit",action="store_true",
                   help="Force reconsideration of completed targets regardless of the ledger cooldown")
    p.add_argument("--search-profile",choices=tuple(SEARCH_PROFILES),default="fast",
                   help="search configuration: fast (Profile J; default, about 7x less CPU in the "
                        "158-case study) or standard")
    p.add_argument("--ranking",choices=("yield","snr"),default="yield",
                   help="target order: expected new detectable planets (default; scripts/target_yield.py) "
                        "or the older fixed-radius S/N proxy")
    p.add_argument("--assumed-sectors",type=float,default=target_yield.Settings.assumed_sectors,
                   help="sectors assumed by the yield ranking when a row has neither observed_days nor sectors")
    p.add_argument("--calibration-dir",type=Path,
                   help="yield-ranking calibration directory; default data/calibration")
    a=p.parse_args()
    project_root=Path(__file__).resolve().parents[1]
    if a.transit_search_passes is None:
        a.transit_search_passes = 5 if a.search_profile == "fast" else 1
    if a.detection_calibration is None:
        if a.search_profile == "fast":
            a.detection_calibration = project_root / "data/calibration/profile-j-fivepass-20261010-fa05-verified.json"
        else:
            calibration_dir = a.calibration_dir or project_root / "data/calibration"
            a.detection_calibration = calibration_dir / "detection_threshold.json"
    try:
        detection_calibration = json.loads(a.detection_calibration.read_text())
        if a.detection_threshold is None:
            a.detection_threshold = float(detection_calibration["threshold"])
    except (OSError, ValueError, KeyError, TypeError) as error:
        p.error(f"cannot load a detection threshold from {a.detection_calibration}: {error}")
    signal.signal(signal.SIGTERM,lambda signum,frame: sys.exit(128+signum))
    if a.max_targets is not None and a.max_targets < 1:
        p.error("--max-targets must be positive")
    if not 1<=a.target_workers<=4:
        p.error("--target-workers must be between 1 and 4")
    if not 1<=a.download_workers<=4:
        p.error("--download-workers must be between 1 and 4")
    if a.target_workers*a.download_workers>4:
        p.error("--target-workers times --download-workers must not exceed the four-download MAST budget")
    if a.cpu_budget<1 or a.threads<0:
        p.error("--cpu-budget must be positive and --threads nonnegative")
    if a.cpu_budget>(os.cpu_count() or 1):
        p.error("--cpu-budget exceeds available logical CPUs")
    if a.target_workers>1 and a.threads==0:
        p.error("parallel targets require a fixed positive --threads budget")
    if a.threads==0 and a.cpu_budget<(os.cpu_count() or 1):
        a.threads=a.cpu_budget
    if a.threads and a.target_workers*a.threads>a.cpu_budget:
        p.error("target-workers times threads exceeds --cpu-budget")
    if not math.isfinite(a.pause) or a.pause<0:
        p.error("--pause must be finite and nonnegative")
    if a.revisit_after_days < 0:
        p.error("--revisit-after-days cannot be negative")
    if not math.isfinite(a.search_timeout) or a.search_timeout <= 0:
        p.error("--search-timeout must be finite and positive")
    if not math.isfinite(a.min_exptime) or a.min_exptime < 0:
        p.error("--min-exptime must be finite and nonnegative")
    if not math.isfinite(a.min_exptime_spoc) or a.min_exptime_spoc < 0:
        p.error("--min-exptime-spoc must be finite and nonnegative")
    if not math.isfinite(a.network_timeout) or a.network_timeout <= 0:
        p.error("--network-timeout must be finite and positive")
    if (not math.isfinite(a.max_sector_noise_ratio) or a.max_sector_noise_ratio < 0
            or 0 < a.max_sector_noise_ratio < 1):
        p.error("--max-sector-noise-ratio must be 0 (disabled) or at least 1")
    if not math.isfinite(a.assumed_sectors) or a.assumed_sectors <= 0:
        p.error("--assumed-sectors must be finite and positive")
    if not math.isfinite(a.detection_threshold) or a.detection_threshold <= 0:
        p.error("--detection-threshold must be finite and positive")
    if a.search_profile!="standard" and a.periods:
        p.error("--periods fixes the period grid; it cannot be combined with --search-profile fast")
    if not math.isfinite(a.max_period) or not 0 < a.min_period < a.max_period:
        p.error("need 0 < --min-period < --max-period")
    yield_model=None
    if a.ranking=="yield":
        try:
            yield_model=target_yield.load_model(a.calibration_dir)
        except ValueError as error:
            p.error(str(error))
    a.output.mkdir(parents=True,exist_ok=True)
    a.require_distance = not a.allow_missing_distance
    rows=read_targets(a.targets)
    if not any(val(row,"teff_k") is not None for row in rows):
        p.error("target catalog has no valid Teff values required by the K/M stellar cuts; "
                "join target-yield rankings to the original TIC metadata before surveying")
    a.source_revision="unknown"
    try:
        head=subprocess.run(["git","rev-parse","HEAD"],cwd=project_root,
                            check=True,capture_output=True,text=True).stdout.strip()
        dirty=subprocess.run(["git","diff","--quiet","HEAD"],cwd=project_root).returncode != 0
        a.source_revision=head+("-dirty" if dirty else "")
    except (OSError,subprocess.CalledProcessError):
        pass
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
    yield_settings=target_yield.Settings(min_period=a.min_period,max_period=a.max_period,
                                         assumed_sectors=a.assumed_sectors)
    expected={}
    def ranking_score(index):
        if index not in expected:
            if yield_model is None:
                expected[index]=local_snr_proxy(rows[index])
            else:
                result=target_yield.expected_new_planets(rows[index],yield_model,yield_settings)
                expected[index]=None if result is None else result["expected_new_planets"]
        return expected[index]
    # Highest expected yield (or S/N proxy) first; unscorable rows last.
    eligible.sort(key=lambda index:(ranking_score(index) is None,-(ranking_score(index) or 0.0),index))
    selected=eligible[:a.max_targets] if a.max_targets is not None else eligible
    selected_set=set(selected)
    # Record the local physical proxy for every catalog row so selection is
    # auditable and repeatable without an external decision service.
    decision_fields=["input_row","tic_id","eligible","selected","expected_new_planets",
                     "local_snr_proxy","status"]
    decision_rows=[]
    for index,row in enumerate(rows):
        ident=row.get("tic_id",row.get("tic",row.get("ID","")))
        score=local_snr_proxy(row)
        if index in duplicate_rows: state="duplicate_tic_id"
        elif index in ledger_deferred: state="ledger_deferred"
        elif index in selected_set: state="selected"
        elif index in eligible: state="budget_deferred"
        else: state="ineligible"
        planets=ranking_score(index) if yield_model is not None and index in catalog_eligible else None
        decision_rows.append({"input_row":index+1,"tic_id":ident,
            "eligible":index in catalog_eligible,"selected":index in selected_set,
            "expected_new_planets":"" if planets is None else planets,
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
    a.build_identity=capture_build_identity(project_root,binary)
    a.source_revision=a.build_identity["source_revision"]+(
        "-dirty" if a.build_identity["working_tree_dirty"] else "")
    if work_rows and a.target_workers==1:
        configure_mast_timeout(a.network_timeout)
        install_default_timeouts(a.network_timeout)
    failures=0
    def claim(row):
        return claim_target(ledger,run_id,tic(row),row,config_sha)
    @contextmanager
    def selected_pool():
        if work_rows and a.target_workers>1:
            with target_pool(a.target_workers,a.network_timeout) as pool: yield pool
        else: yield None
    orchestration_started=time.monotonic()
    with selected_pool() as pool:
        for i,(index,row,result) in enumerate(target_results(work_rows,a,binary,claim,pool),1):
            ident=tic(row)
            if result.get("status") in ("error","partial","busy"): failures+=1
            summary_by_index[index]=result
            if result.get("status")!="busy" and not record_target(ledger,run_id,ident,row,config_sha,result):
                # Another run recovered an expired lease; preserve its state.
                result={**result,"status":"lost_claim",
                        "reason":"lease expired and another run recovered this target; result not recorded"}
                failures+=1
                summary_by_index[index]=result
            print(f"[{i}/{len(work_rows)}] TIC {result.get('tic_id')}: {result['status']}",flush=True)
    orchestration_wall=time.monotonic()-orchestration_started
    for index,result in summary_by_index.items():
        summary.append(result)
    # Preserve input order in the summary while recording the actual selection.
    index_for_result={}
    for index,row in enumerate(rows):
        try: index_for_result.setdefault(tic(row),[]).append(index)
        except ValueError: pass
    summary.sort(key=lambda item:min(index_for_result.get(str(item.get("tic_id","")),[len(rows)])))
    atomic_write_text(a.output/"batch_summary.json",json.dumps(summary,indent=2,default=str)+"\n")
    # Keep transient author failures in an author-specific retry queue. A
    # partial target is cooldown-protected in the shared ledger, so rerunning
    # its failed author does not repeat the pipeline that already succeeded.
    retryable_author_statuses={"search_failed","no_sector_metadata","no_usable_products"}
    input_by_tic={}
    for row in rows:
        try: input_by_tic[tic(row)]=row
        except ValueError: pass
    retry_rows={}
    for result in summary:
        ident=str(result.get("tic_id",""))
        if ident not in input_by_tic: continue
        for author_run in result.get("authors",[]):
            if author_run.get("status") not in retryable_author_statuses: continue
            author=author_run.get("author","unknown")
            retry_rows.setdefault(author,[]).append({
                **input_by_tic[ident],"retry_author":author,
                "retry_status":author_run.get("status"),
                "retry_reason":author_run.get("reason","")})
    retry_queue_paths={}
    for author,queue_rows in retry_rows.items():
        suffix="".join(ch.lower() if ch.isalnum() else "_" for ch in str(author)).strip("_")
        queue_path=a.output/f"author-retry-queue-{suffix}.csv"
        fields=list(rows[0].keys())+["retry_author","retry_status","retry_reason"]
        atomic_write_csv(queue_path,fields,queue_rows)
        retry_queue_paths[author]={"path":str(queue_path.resolve()),"targets":len(queue_rows),
            "next_run_authors":[author]}
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
        "run_id":run_id,"ledger":str(ledger_path.resolve()),"ranking":("expected_new_planets" if yield_model is not None else "local_snr_proxy"),
        "config_sha256":config_sha,"targets_csv_sha256":sha256_file(a.targets),
        "astrosearch_binary":str(binary),"binary_sha256":sha256_file(binary),
        "source_revision":a.source_revision,
        "build_identity":a.build_identity,
        "candidate_registry":str((a.candidate_registry.expanduser() if a.candidate_registry
                                   else a.output.resolve()/"candidate_registry.sqlite")),
        "revisit_after_days":a.revisit_after_days,"forced_revisit":a.revisit,
        "catalog_eligible":len(catalog_eligible),"ledger_deferred":len(ledger_deferred),
        "selected":len(selected_set),"failed_targets":failures,"status":run_status,
        "author_retry_queues":retry_queue_paths,
        "orchestration":{"target_workers":a.target_workers,"threads_per_search":a.threads,
            "download_workers_per_target":a.download_workers,
            "max_concurrent_mast_downloads":a.target_workers*a.download_workers,
            "cpu_budget":a.cpu_budget,"submission_pause_s":a.pause,
            "wall_seconds":orchestration_wall,
            "product_preflight_rejections":sum(run.get("status") in
                ("no_recent_ffi_products","insufficient_product_baseline")
                for result in summary for run in result.get("authors",[]))},
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
