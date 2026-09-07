#!/usr/bin/env bash
# correction/run_charge_flip.sh
#
# Run the lepton charge-misidentification (charge-flip) study on the W signal
# MC and KEEP THE LOG.
#
# charge_flip.C prints the record -- inclusive rates with Clopper-Pearson
# intervals, the per-y and per-pT tables, the match-quality summary, the
# flip-induced charge-asymmetry bias per analysis bin, and (electrons) the
# eleCharge / eleTrkCharge consistency numbers. ROOT writes all of it to
# stdout only, hence this wrapper: same convention as run_qcd_abcd.sh, one log
# per flavor under correction/logs/.
#
# Usage:  ./run_charge_flip.sh [mu|ele|both]      (default: both)
#
# Outputs: logs/charge_flip_<chan>.log       full console output (THE record)
#          rootfile/charge_flip_<chan>.root  count histos + TEfficiency objects
#          plots/charge_flip_<chan>/         rate plots + match-quality shapes
#
# bash-3.2 safe (macOS stock /bin/bash): no associative arrays, no mapfile.
set -euo pipefail

cd "$(dirname "$0")"

WHICH="${1:-both}"
case "$WHICH" in
  mu|ele|both) ;;
  *) echo "usage: $0 [mu|ele|both]" >&2; exit 2 ;;
esac

mkdir -p logs rootfile plots

# Pre-build once so the two jobs never race on the ACLiC artifacts.
echo "[BUILD] compiling charge_flip.C ..."
if ! root -l -b -q -e '.L charge_flip.C+' > logs/charge_flip_build.log 2>&1; then
  echo "[FAIL] compilation failed -- see logs/charge_flip_build.log" >&2
  exit 1
fi

run_one() {  # $1 = mu | ele
  chan="$1"
  log="logs/charge_flip_${chan}.log"

  printf '[RUN ] %-3s -> %s\n' "$chan" "$log"
  SECONDS=0
  if root -l -b -q "charge_flip.C+(\"${chan}\")" > "$log" 2>&1; then
    status="OK"
  else
    status="FAIL"
  fi
  elapsed=$SECONDS

  # Echo the headline numbers so the terminal still tells you what happened.
  grep -E '^\[CONFIG\]|^\[RESULT\]|^\[INPUT\]|^\[WARN\]|^\[ERR|^\[FATAL' "$log" 2>/dev/null | sed 's/^/       /' || true
  printf '[%-4s] %-3s  %ds   log: %s\n' "$status" "$chan" "$elapsed" "$log"
  [ "$status" = "OK" ] || return 1
}

fail=0
if [ "$WHICH" = "both" ] || [ "$WHICH" = "mu" ];  then run_one mu  || fail=$((fail+1)); fi
if [ "$WHICH" = "both" ] || [ "$WHICH" = "ele" ]; then run_one ele || fail=$((fail+1)); fi

echo "============================================================"
if [ "$fail" -eq 0 ]; then
  echo "[DONE] failures=0   logs in $(pwd)/logs/"
else
  echo "[DONE] failures=${fail}  -- check the logs above" >&2
  exit 1
fi
