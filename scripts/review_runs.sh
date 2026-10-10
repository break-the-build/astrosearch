#!/bin/sh
# Reproduce the TIC 4206066 review runs, inverted-flux null controls, and
# serial/parallel identity checks.
#
# usage: scripts/review_runs.sh [BINARY] [OUT_DIR]
#   BINARY   astrosearch executable (default: target/release/astrosearch)
#   OUT_DIR  where candidate CSVs and logs go (default: target/review)
#
# Known answer: a ~500 ppm transit at P = 3.18 d (Rabtsevich 2026,
# doi:10.5281/zenodo.22967456). A "hit" is a periodic_dip within 0.005 d of 3.182.
set -eu

root=$(cd "$(dirname "$0")/.." && pwd)
bin=${1:-$root/target/release/astrosearch}
out=${2:-$root/target/review}
data=$root/data/tic4206066
mkdir -p "$out"

summarize() { # NAME: elapsed, segments, best dip, and rank/score of the 3.18 d planet
  name=$1
  elapsed=$(awk '/^real/ {print $2}' "$out/$name.log")
  segments=$(grep -c '^segment' "$out/$name.log" || true)
  best=$(awk -F, '$2=="periodic_dip" {printf "P=%.4f score=%.2f depth=%.0f n=%s", $5, $4, $8, $10; exit}' "$out/$name.csv")
  planet=$(awk -F, '$2=="periodic_dip" && ($5-3.182)^2 < 0.005^2 {printf "class_rank=%s score=%.2f depth=%.0f", $3, $4, $8; exit}' "$out/$name.csv")
  echo "$name | ${elapsed}s | segments=$segments | best dip: ${best:-none} | 3.18 d: ${planet:-not found}"
}

run() { # NAME INPUT [ARGS...]
  name=$1; input=$2; shift 2
  /usr/bin/time -p "$bin" --input "$input" --output "$out/$name.csv" "$@" > "$out/$name.log" 2>&1
  summarize "$name"
}

# Null control: inverting normalized flux turns dips into brightenings, so any
# periodic_dip found afterwards comes from noise or systematics.
for f in s06_TESS-SPOC_1800 s32_SPOC_120 s98_QLP_200 combined; do
  awk -F, 'NR == 1 {print; next} {printf "%s,%.8f,%s\n", $1, 2 - $2, $3}' \
    "$data/$f.csv" > "$out/inverted_$f.csv"
done

echo "== Full default search (review commands)"
for f in s06_TESS-SPOC_1800 s32_SPOC_120 s98_QLP_200 combined; do
  run "full_$f" "$data/$f.csv"
done

echo "== Narrow search near 3.18 d"
for trend in 0.5 1.0; do
  for f in s06_TESS-SPOC_1800 s32_SPOC_120 s98_QLP_200 combined; do
    run "narrow_${trend}_$f" "$data/$f.csv" --min-period 3.0 --max-period 3.4 --trend-days "$trend"
  done
done

echo "== Inverted-flux null controls"
for f in s06_TESS-SPOC_1800 s32_SPOC_120 s98_QLP_200 combined; do
  run "null_full_$f" "$out/inverted_$f.csv"
done
run "null_narrow_combined" "$out/inverted_combined.csv" --min-period 3.0 --max-period 3.4

echo "== Serial versus parallel output identity"
run "serial_full_s06" "$data/s06_TESS-SPOC_1800.csv" --threads 1
cmp "$out/serial_full_s06.csv" "$out/full_s06_TESS-SPOC_1800.csv" && echo "full s06: identical"
run "serial_narrow_combined" "$data/combined.csv" --min-period 3.0 --max-period 3.4 --threads 1
cmp "$out/serial_narrow_combined.csv" "$out/narrow_0.5_combined.csv" && echo "narrow combined: identical"
