#!/bin/sh
# Check that --scoring white reproduces a previous build, and time all modes.
#
# usage: scripts/compare_scoring.sh BASELINE_BINARY [NEW_BINARY] [OUT_DIR]
# Compares the columns both builds share (rank..odd_even_sigma and note) for
# full and narrow searches of the TIC 4206066 inputs.
set -eu

root=$(cd "$(dirname "$0")/.." && pwd)
old=$1
new=${2:-$root/target/release/astrosearch}
out=${3:-$root/target/compare-scoring}
data=$root/data/tic4206066
mkdir -p "$out"

shared() { # CSV -> shared columns (first 11 fields and the quoted note)
  awk -F, '{n = $0; sub(/^([^,]*,){11}/, "", n); sub(/^[^"]*/, "", n); print $1","$2","$3","$4","$5","$6","$7","$8","$9","$10","$11","n}' "$1"
}

timed() { # NAME BINARY ARGS...
  name=$1; bin=$2; shift 2
  /usr/bin/time -p "$bin" --output "$out/$name.csv" "$@" > "$out/$name.log" 2>&1
  awk '/^real/ {print $2}' "$out/$name.log"
}

for case in full_s06:s06_TESS-SPOC_1800: full_s32:s32_SPOC_120: full_combined:combined: \
  narrow_s98:s98_QLP_200:narrow narrow_combined:combined:narrow; do
  name=${case%%:*}; rest=${case#*:}; file=${rest%%:*}; kind=${rest#*:}
  set -- --input "$data/$file.csv"
  if [ "$kind" = narrow ]; then set -- "$@" --min-period 3.0 --max-period 3.4; fi
  t_old=$(timed "${name}_baseline" "$old" "$@")
  t_white=$(timed "${name}_white" "$new" "$@" --scoring white)
  t_red=$(timed "${name}_correlated" "$new" "$@")
  shared "$out/${name}_white.csv" > "$out/${name}_white.csv.shared"
  if shared "$out/${name}_baseline.csv" | cmp -s - "$out/${name}_white.csv.shared"; then
    same=identical
  else
    same=DIFFERENT
  fi
  echo "$name | baseline ${t_old}s | white ${t_white}s ($same) | correlated ${t_red}s"
done
