#!/usr/bin/env python3
"""Timestamped discovery preflight from astrosearch candidate CSV output.

Automated: NASA Exoplanet Archive TOIs/planet systems and requested SPOC TCE
files. CTOI/ExoFOP, pixels, Gaia neighbors, depth consistency, and TRICERATOPS
remain explicit human review. An empty API result is never called a discovery.
"""
import argparse, csv, datetime, json, os, re, sys, tempfile
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
csv.field_size_limit(2**31 - 1)  # candidate rows carry long evidence/alias cells

TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
SESSION = requests.Session()
SESSION.headers["User-Agent"] = "astrosearch-research/0.1"
SESSION.mount("https://", HTTPAdapter(max_retries=Retry(
    total=4, connect=4, read=3, status=3, backoff_factor=0.5,
    status_forcelist=(429, 500, 502, 503, 504),
    allowed_methods=frozenset({"GET"}), raise_on_status=False)))

def tap(query):
    r = SESSION.get(TAP, params={"query": query, "format": "json"}, timeout=45)
    r.raise_for_status()
    data = r.json()
    if not isinstance(data, list): raise RuntimeError(f"unexpected TAP response: {data}")
    return data

class Links(HTMLParser):
    def __init__(self):
        super().__init__(); self.hrefs=[]
    def handle_starttag(self, tag, attrs):
        if tag == "a":
            href=dict(attrs).get("href")
            if href: self.hrefs.append(href)

def select_tce_files(hrefs, requested_sectors, index_url):
    requested=set(requested_sectors)
    files={}
    for href in hrefs:
        match=re.search(r"-s(\d{4})-s(\d{4})_dvr-tcestats\.csv$",href)
        if not match: continue
        start,stop=map(int,match.groups())
        if any(start <= sector <= stop for sector in requested):
            files[urljoin(index_url,href)]=(start,stop)
    return files

def ephemeris_matches(candidate_rows, period, epoch_btjd, duration_days=None):
    """Match periods including common 1:2 aliases, then compare transit phase."""
    if period is None or epoch_btjd is None or period <= 0:
        return []
    matches=[]
    for candidate in candidate_rows:
        try:
            cp=float(candidate["period_days"]); ce=float(candidate["epoch"])
            cd=float(candidate["duration_days"])
        except (KeyError,TypeError,ValueError):
            continue
        if cp <= 0: continue
        ratio=cp/period
        alias=min((1.0,2.0,0.5),key=lambda value:abs(ratio-value))
        if abs(ratio-alias)/alias > 0.01: continue
        phase_period=min(cp,period)
        phase=abs((ce-epoch_btjd) % phase_period)
        phase=min(phase,phase_period-phase)
        tolerance=max(cd, (duration_days or 0.0))/2
        if phase <= tolerance:
            matches.append({"candidate_class_rank":candidate.get("class_rank"),
                            "period_alias_ratio":alias,"phase_offset_days":phase,
                            "phase_tolerance_days":tolerance})
    return matches

def match_catalog_ephemerides(rows, candidate_rows, epoch_key="pl_tranmid"):
    matches=[]
    for row in rows or []:
        try:
            period=float(row["pl_orbper"])
            epoch=float(row[epoch_key])
            duration=float(row.get("pl_trandurh") or 0)/24
        except (KeyError,TypeError,ValueError):
            continue
        # Archive transit midpoints are BJD; astrosearch candidates are BTJD.
        if epoch > 1_000_000: epoch -= 2_457_000
        for match in ephemeris_matches(candidate_rows,period,epoch,duration):
            matches.append({"catalog_row":row,**match})
    return matches

def check_spoc_tces(tic_id, candidate_rows, requested_sectors):
    """Read the public MAST TCE index and inspect the requested sector files."""
    index="https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_tce.html"
    response=SESSION.get(index,timeout=45); response.raise_for_status()
    links=Links(); links.feed(response.text)
    files=select_tce_files(links.hrefs,requested_sectors,index)
    checked=[]; matches=[]
    for url,(start,stop) in sorted(files.items()):
        r=SESSION.get(url,timeout=90); r.raise_for_status()
        rows=list(csv.DictReader(r.text.splitlines()))
        if not rows: checked.append({"url":url,"status":"empty"}); continue
        keys={re.sub(r"[^a-z0-9]","",k.lower()):k for k in rows[0]}
        id_key=next((keys[k] for k in ("ticid","tic_id","targetid","tessinputcatalogid") if k in keys),None)
        period_key=next((keys[k] for k in ("tceperiod","period","plorbper") if k in keys),None)
        epoch_key=next((keys[k] for k in ("tceepoch","tceepochbtjd","tceepochbjd","tceepochbkjd") if k in keys),None)
        duration_key=next((keys[k] for k in ("tceduration","duration","pltrandurh") if k in keys),None)
        if not id_key:
            checked.append({"url":url,"status":"unrecognized_schema","columns":list(rows[0])}); continue
        target_rows=[]
        for row in rows:
            try: row_tic=str(int(float(row[id_key])))
            except (TypeError,ValueError): continue
            if row_tic==str(tic_id): target_rows.append(row)
        for row in target_rows:
            pvalue=None
            if period_key:
                try: pvalue=float(row[period_key])
                except (TypeError,ValueError): pass
            try: evalue=float(row[epoch_key]) if epoch_key else None
            except (TypeError,ValueError): evalue=None
            if evalue is not None and evalue > 1_000_000: evalue -= 2_457_000
            try: dvalue=float(row[duration_key])/24 if duration_key else None
            except (TypeError,ValueError): dvalue=None
            matches_ephemeris=ephemeris_matches(candidate_rows,pvalue,evalue,dvalue)
            similar=bool(matches_ephemeris) or any(
                pvalue is not None and pvalue>0 and abs(pvalue-float(candidate["period_days"]))/float(candidate["period_days"])<=0.01
                for candidate in candidate_rows if candidate.get("period_days"))
            matches.append({"file":url,"sector_start":start,"sector_stop":stop,
                            "period_match_within_1pct":similar,
                            "ephemeris_matches":matches_ephemeris,
                            "epoch_column":epoch_key,"duration_column":duration_key,"row":row})
        checked.append({"url":url,"sector_start":start,"sector_stop":stop,
                        "rows":len(rows),"tic_rows":len(target_rows),
                        "schema_status":"ok" if period_key else "no_period_column"})
    return {"requested_sectors":sorted(set(requested_sectors)),"files_checked":checked,"matches":matches}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tic-id",required=True,type=int)
    p.add_argument("--candidates",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--tce-sectors",default="",help="Comma-separated sectors to inspect at MAST")
    a=p.parse_args()
    now=datetime.datetime.now(datetime.timezone.utc).isoformat()
    errors={}
    try:
        toi=tap(f"select toi,ctoi_alias,tid,tfopwg_disp,pl_orbper,pl_tranmid,pl_trandurh,pl_trandep,ra,dec from toi where tid={a.tic_id}")
    except Exception as error:
        toi=None; errors["toi_query"]=f"{type(error).__name__}: {error}"
    try:
        planets=tap(f"select pl_name,tic_id,pl_orbper,pl_tranmid,pl_trandur as pl_trandurh,discoverymethod from ps where tic_id='TIC {a.tic_id}'")
    except Exception as error:
        planets=None; errors["confirmed_planet_query"]=f"{type(error).__name__}: {error}"
    with a.candidates.open(newline="",encoding="utf-8-sig") as f: candidates=list(csv.DictReader(f))
    dips=[r for r in candidates if r.get("class") in ("periodic_dip", "duo_transit")][:12]
    sector_ids=[]
    for item in a.tce_sectors.split(","):
        if item.strip(): sector_ids.append(int(item.strip()))
    tce=None
    if sector_ids:
        try:
            tce=check_spoc_tces(a.tic_id,dips,sector_ids)
        except Exception as error:
            errors["spoc_tce_query"]=f"{type(error).__name__}: {error}"
    sectors=", ".join(map(str,sector_ids)) or "not provided"
    toi_url=f"https://exoplanetarchive.ipac.caltech.edu/overview/TIC%20{a.tic_id}"
    exofop_url=f"https://exofop.ipac.caltech.edu/tess/target.php?id={a.tic_id}"
    tce_url="https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_tce.html"
    toi_ephemeris_matches=match_catalog_ephemerides(toi,dips)
    confirmed_ephemeris_matches=match_catalog_ephemerides(planets,dips)
    ctoi_aliases=[row.get("ctoi_alias") for row in (toi or []) if row.get("ctoi_alias")]
    lines=[
        f"# Candidate preflight: TIC {a.tic_id}","",f"Generated UTC: {now}","",
        "> Screening artifact only. A blank catalog response or a high search score does not establish a new planet.","",
        "## Automated catalog checks","",
        f"- NASA Exoplanet Archive TOI rows: **{len(toi) if toi is not None else 'unavailable'}**. TOI is refreshed on a schedule; check the live [TIC record]({toi_url}).",
        f"- Community TOI aliases present in the Archive response: {json.dumps(ctoi_aliases)}. Verify current CTOI status directly in ExoFOP.",
        f"- NASA confirmed-planet rows: **{len(planets) if planets is not None else 'unavailable'}**.",
        f"- TOI period/epoch aliases matching these triggers: {json.dumps(toi_ephemeris_matches,default=str)}",
        f"- Confirmed-planet period/epoch aliases matching these triggers: {json.dumps(confirmed_ephemeris_matches,default=str)}",
        f"- TOI response: {json.dumps(toi,default=str)}" if toi else ("- TOI response: no matching TIC row." if toi is not None else "- TOI response: query failed; see service errors below."),
        f"- Confirmed planets: {json.dumps(planets,default=str)}" if planets else ("- Confirmed planets: no matching TIC row." if planets is not None else "- Confirmed planets: query failed; see service errors below."),
        f"- ExoFOP target / CTOI: [open TIC {a.tic_id}]({exofop_url}) and check current Community TOIs. Authenticated ExoFOP data are not queried here.",
        f"- SPOC TCE product index: [MAST sector TCE tables]({tce_url}); requested sectors: {sectors}.",
        f"- Automated TCE file check: {json.dumps(tce,default=str)}" if tce else ("- SPOC TCE query failed; see service errors below." if sector_ids else "- SPOC TCE files were not queried because --tce-sectors was not supplied."),
        "- A missing row only means no TIC match in the sector files fetched; also inspect multi-sector TCEs and keep catalog coverage limitations in the report.",
        f"- Service errors: {json.dumps(errors,default=str) if errors else 'none'}",
        "","## Candidate triggers from this run","",
        "| Class rank | Period (d) | Epoch (BTJD) | Duration (d) | Depth (ppm) | Score | Events | Odd/even sigma |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for row in dips:
        lines.append("| "+" | ".join(row.get(k,"") for k in ("class_rank","period_days","epoch","duration_days","depth_ppm","score","n_events","odd_even_sigma"))+" |")
    if not dips: lines.append("| No periodic-dip triggers in this CSV | | | | | | | |")
    lines += [
        "","## Required validation before reporting a candidate","",
        "- [ ] Search TOI, CTOI, confirmed-planet, and SPOC TCE catalogs by TIC and by period/epoch aliases. Record date, identifiers, and dispositions.",
        "- [ ] Repeat independently on SPOC and QLP products; compare ephemerides and fitted depths.",
        "- [ ] Fit transit models by sector. Check odd/even depths, a secondary eclipse at phase 0.5, and depth consistency across sectors.",
        "- [ ] Inspect Gaia neighbors, RUWE, TIC contamination ratio, and target-pixel difference images/centroids. Rule out off-target eclipsing binaries.",
        "- [ ] Run TRICERATOPS with the correct target-pixel aperture, Gaia field, sectors, stellar parameters, and dilution assumptions. Save versions, inputs, scenarios, and FPP.",
        "- [ ] Measure recovery and false alarms for this cadence/target class with injections and correlated-noise null curves. The astrosearch score is not a probability.",
        "- [ ] If compelling, publish a time-stamped ephemeris and uncertainty before future transits, then prepare the ExoFOP CTOI submission with the validation package.",
        "","## Transit ephemeris","",
        "Candidate epochs are fitted window centers in BTJD. Propagate fitted epoch and period uncertainties before publishing: T_n = T_0 + nP and sigma(T_n) = sqrt(sigma(T_0)^2 + n^2 sigma(P)^2); include covariance when available. Screening output does not estimate these uncertainties.",
        "","## Sources","",
        "- [NASA Exoplanet Archive TAP](https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html)",
        "- [MAST TESS TCE bulk downloads](https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_tce.html)",
        "- [ExoFOP-TESS](https://exofop.ipac.caltech.edu/tess/)",
        "- [TRICERATOPS documentation](https://triceratops.readthedocs.io/)",""]
    a.output.parent.mkdir(parents=True,exist_ok=True)
    temp=None
    try:
        with tempfile.NamedTemporaryFile("w",encoding="utf-8",dir=a.output.parent,
                prefix=f".{a.output.name}.",suffix=".tmp",delete=False) as stream:
            temp=Path(stream.name); stream.write("\n".join(lines)); stream.flush(); os.fsync(stream.fileno())
        os.replace(temp,a.output)
    finally:
        if temp is not None and temp.exists(): temp.unlink()
    print(a.output)
    return 1 if errors else 0

if __name__=="__main__":
    try: raise SystemExit(main())
    except Exception as e:
        print(f"preflight: {type(e).__name__}: {e}",file=sys.stderr)
        raise SystemExit(2)
