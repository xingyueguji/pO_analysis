#!/usr/bin/env bash
# correction/run_idiso_sf.sh
#
# Run the ELECTRON ID+ISO SF CROSS-CHECK (2026-09-24) and KEEP THE LOGS.
# The definition lives in idiso_sf_common.h: a W-sample efficiency of
# eleMVAIdWP90 && eleMVAIsoWP90, measured with an electron-only W+Z fit that
# runs in the Combine fork (test/run_pO_idiso_sf.sh), separate from the
# nominal fit stream, and compared with the EGM 2025Prompt wp90iso SF.
#
# Usage:  ./run_idiso_sf.sh [skim|inputs|plots|all] [samples...]
#   skim    the separate skim, one job per sample (default: all 7,
#           Data Wp Wm DY DYtau Wptau Wmtau)       -> rootfile/idiso_sf_ele/skim_<s>.root
#   inputs  templates, QCD, Combine inputs + the composition record
#                                                  -> rootfile/idiso_sf_ele/combine_input_idiso_*.root
#   plots   efficiencies, SF, EGM overlay; needs the fit summaries downloaded
#           from lxplus (FORK_TEST overrides the fork's test/ dir)
#                                                  -> plots/idiso_sf_ele/
#   all     skim + inputs (everything before the lxplus fit)
#
# Logs: logs/idiso_sf_ele_skim_<sample>.log, logs/idiso_sf_ele_inputs.log,
#       logs/idiso_sf_ele_plots.log -- ROOT prints the record to stdout only.
#
# bash-3.2 safe (macOS stock /bin/bash): no associative arrays, no mapfile.
set -euo pipefail

cd "$(dirname "$0")"

STAGE="${1:-all}"; [ $# -gt 0 ] && shift
case "$STAGE" in
  skim|inputs|plots|all) ;;
  *) echo "usage: $0 [skim|inputs|plots|all] [samples...]" >&2; exit 2 ;;
esac
SAMPLES="${*:-Data Wp Wm DY DYtau Wptau Wmtau}"

mkdir -p logs rootfile plots

build() {  # $1 = macro; pre-build once so no two jobs race on the ACLiC artifacts
  echo "[BUILD] compiling $1 ..."
  if ! root -l -b -q -e ".L $1+" > "logs/idiso_sf_ele_build.log" 2>&1; then
    echo "[FAIL] compilation of $1 failed -- see logs/idiso_sf_ele_build.log" >&2
    tail -20 logs/idiso_sf_ele_build.log >&2
    exit 1
  fi
}

run_root() {  # $1 = label, $2 = root command, $3 = log, $4 = headline regex
  printf '[RUN ] %-14s -> %s\n' "$1" "$3"
  SECONDS=0
  if root -l -b -q "$2" > "$3" 2>&1; then status="OK"; else status="FAIL"; fi
  # root exits 0 even when the macro bails out -> also check for a FATAL line
  grep -q '^\[FATAL\]' "$3" && status="FAIL"
  grep -E "$4" "$3" 2>/dev/null | sed 's/^/       /' || true
  printf '[%-4s] %-14s %ds   log: %s\n' "$status" "$1" "$SECONDS" "$3"
  [ "$status" = "OK" ]
}

fail=0
if [ "$STAGE" = "skim" ] || [ "$STAGE" = "all" ]; then
  build idiso_sf_skim.C
  for s in $SAMPLES; do
    run_root "skim $s" "idiso_sf_skim.C+(\"$s\")" "logs/idiso_sf_ele_skim_${s}.log" \
             '^\[INPUT\]|^\[RESULT\]|^\[GENMATCH\]|^\[WARN\]|^\[FATAL\]' || fail=$((fail+1))
  done
fi
if [ "$STAGE" = "inputs" ] || [ "$STAGE" = "all" ]; then
  if [ "$fail" -gt 0 ]; then
    echo "[SKIP] inputs: the skim had failures" >&2
  else
    build idiso_sf_inputs.C
    run_root "inputs" "idiso_sf_inputs.C+()" "logs/idiso_sf_ele_inputs.log" \
             '^\[RESULT\]|^\[EPS\]|^\[ZTNP\]|^\[EGMPRED\]|^\[ABCD\]|^\[OUT\]|^\[WARN\]|^\[FATAL\]' || fail=$((fail+1))
  fi
fi
if [ "$STAGE" = "plots" ]; then
  build idiso_sf_plots.C
  run_root "plots" "idiso_sf_plots.C+()" "logs/idiso_sf_ele_plots.log" \
           '^\[FIT\]|^\[SF\]|^\[EGM\]|^\[PULL\]|^\[OUT\]|^\[WARN\]|^\[FATAL\]' || fail=$((fail+1))
fi

echo "============================================================"
if [ "$fail" -eq 0 ]; then
  echo "[DONE] failures=0   logs in $(pwd)/logs/"
else
  echo "[DONE] failures=${fail}  -- check the logs above" >&2
  exit 1
fi
