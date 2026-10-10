#!/usr/bin/env python3
"""Write future BJD_TDB transit windows from a fitted ephemeris and uncertainties."""
import argparse, csv, math
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--period",type=float,required=True,help="period in days")
    p.add_argument("--epoch",type=float,required=True,help="reference mid-transit BTJD = BJD_TDB - 2457000")
    p.add_argument("--period-sigma",type=float,required=True,help="period uncertainty in days")
    p.add_argument("--epoch-sigma",type=float,required=True,help="epoch uncertainty in days")
    p.add_argument("--after-btjd",type=float,required=True,help="only emit events later than this BTJD")
    p.add_argument("--through-btjd",type=float,required=True)
    p.add_argument("--covariance",type=float,default=0.,help="covariance(epoch,period) in days squared")
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    if not all(math.isfinite(x) for x in (a.period,a.epoch,a.period_sigma,a.epoch_sigma,a.after_btjd,a.through_btjd,a.covariance)):
        p.error("ephemeris values must be finite")
    if a.period<=0 or a.period_sigma<0 or a.epoch_sigma<0 or a.through_btjd<=a.after_btjd:
        p.error("period must be positive, errors nonnegative, and through-btjd after after-btjd")
    first=math.floor((a.after_btjd-a.epoch)/a.period)+1
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open("w",newline="") as f:
        w=csv.writer(f); w.writerow(("cycle","mid_transit_btjd_tdb","sigma_days","window_3sigma_start_btjd","window_3sigma_end_btjd"))
        n=first
        while True:
            center=a.epoch+n*a.period
            if center>a.through_btjd: break
            variance=a.epoch_sigma**2+n*n*a.period_sigma**2+2*n*a.covariance
            if variance<0: p.error(f"negative propagated variance at cycle {n}; check covariance")
            sigma=math.sqrt(variance)
            w.writerow((n,format(center,".10f"),format(sigma,".10g"),
                        format(center-3*sigma,".10f"),format(center+3*sigma,".10f")))
            n+=1
    print(a.output)

if __name__=="__main__": main()

