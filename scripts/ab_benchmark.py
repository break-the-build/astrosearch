"""Interleaved A/B benchmark of two astrosearch binaries.

Records wall, user and system time, peak RSS (from /usr/bin/time -l) and the
SHA-256 of the candidate CSV for each run. Interleaving A and B spreads
background-load drift evenly over both; user CPU time is less sensitive to
load than wall time, but on an oversubscribed machine it still inflates.
"""
import hashlib
import json
import os
import re
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

a_bin, b_bin, curve, rounds = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
extra = sys.argv[5:]
results = {"A": [], "B": []}
with tempfile.TemporaryDirectory() as tmp:
    for r in range(rounds):
        for label, binary in (("A", a_bin), ("B", b_bin)) if r % 2 == 0 else (("B", b_bin), ("A", a_bin)):
            out = Path(tmp) / f"{label}{r}.csv"
            done = subprocess.run(["/usr/bin/time", "-l", binary, "--input", curve, "--output", str(out), *extra],
                                  capture_output=True, text=True, check=True)
            text = done.stderr
            real, user, sys_t = (float(x) for x in re.search(r"([\d.]+) real\s+([\d.]+) user\s+([\d.]+) sys", text).groups())
            rss = int(re.search(r"(\d+)\s+maximum resident set size", text).group(1))
            sha = hashlib.sha256(out.read_bytes()).hexdigest()[:16]
            results[label].append({"real": real, "user": user, "sys": sys_t, "max_rss_mb": rss / 2**20, "sha": sha})
            print(label, r, results[label][-1], flush=True)
load = os.getloadavg()
summary = {k: {"user_median": statistics.median(x["user"] for x in v),
               "real_median": statistics.median(x["real"] for x in v),
               "rss_mb_median": statistics.median(x["max_rss_mb"] for x in v),
               "digests": sorted({x["sha"] for x in v})} for k, v in results.items()}
summary["user_speedup_A_over_B"] = summary["A"]["user_median"] / summary["B"]["user_median"]
summary["identical_output"] = summary["A"]["digests"] == summary["B"]["digests"] and len(summary["A"]["digests"]) == 1
summary["load_average_at_end"] = load
summary["args"] = extra
print(json.dumps(summary, indent=1))
