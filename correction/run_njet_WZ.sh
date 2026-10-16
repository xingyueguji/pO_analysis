#!/usr/bin/env bash
# correction/run_njet_WZ.sh
#
# Run the jet-multiplicity study of W-/Z-tagged events and KEEP THE LOG.
#
# njet_WZ.C re-runs the skim's W and Z selections on the ntuples (the skim
# outputs carry no jets), so besides the njet distributions its console output
# is the EVENT-IDENTITY cross-check of the selection: "[INFO] RunW selected N
# events" for data must equal the skim cutflow N[8] (skim/logs/skim_W*_Data.log)
# and "RunZ selected N" the skim's Z-peak count -- two independent
# implementations of the same selection. ROOT prints it to stdout only, hence
# this wrapper (same convention as run_charge_flip.sh; added 2026-09-29 when
# the electron selection moved to WP90 ID+iso and the check had to be re-run).
#
# Usage:  ./run_njet_WZ.sh [Wmu|Wel|Zmm|Zee|all]      (default: all)
#
# Outputs: logs/njet_WZ_<chan>.log       full console output (THE record)
#          rootfile/njet_<chan>.root     per-sample njet/jet histograms
#          plots/njet_<chan>/            data vs MC njet, jet pT/eta
#
# bash-3.2 safe (macOS stock /bin/bash): no associative arrays, no mapfile.
set -euo pipefail

cd "$(dirname "$0")"

WHICH="${1:-all}"
case "$WHICH" in
  Wmu|Wel|Zmm|Zee) CHANS="$WHICH" ;;
  all)             CHANS="Wmu Wel Zmm Zee" ;;
  *) echo "usage: $0 [Wmu|Wel|Zmm|Zee|all]" >&2; exit 2 ;;
esac

mkdir -p logs rootfile plots

# Pre-build once so the jobs never race on the ACLiC artifacts.
echo "[BUILD] compiling njet_WZ.C ..."
if ! root -l -b -q -e '.L njet_WZ.C+' > logs/njet_WZ_build.log 2>&1; then
  echo "[FAIL] compilation failed -- see logs/njet_WZ_build.log" >&2
  exit 1
fi

run_one() {  # $1 = Wmu | Wel | Zmm | Zee
  chan="$1"
  log="logs/njet_WZ_${chan}.log"

  printf '[RUN ] %-3s -> %s\n' "$chan" "$log"
  SECONDS=0
  if root -l -b -q "njet_WZ.C+(\"${chan}\")" > "$log" 2>&1; then
    status="OK"
  else
    status="FAIL"
  fi
  elapsed=$SECONDS

  # Echo the headline numbers (the selected-event counts = the identity check).
  grep -E '^\[INPUT\]|selected|^ +data|^\[WARN\]|^\[ERR|^\[FATAL' "$log" 2>/dev/null | sed 's/^/       /' || true
  printf '[%-4s] %-3s  %ds   log: %s\n' "$status" "$chan" "$elapsed" "$log"
  [ "$status" = "OK" ] || return 1
}

fail=0
for c in $CHANS; do run_one "$c" || fail=$((fail+1)); done

echo "============================================================"
if [ "$fail" -eq 0 ]; then
  echo "[DONE] failures=0   logs in $(pwd)/logs/"
else
  echo "[DONE] failures=${fail}  -- check the logs above" >&2
  exit 1
fi
