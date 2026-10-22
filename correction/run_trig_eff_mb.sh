#!/usr/bin/env bash
# correction/run_trig_eff_mb.sh
#
# Run the MB-denominator single-lepton trigger-efficiency study (turn-on,
# data vs W signal MC) and KEEP THE LOG.
#
# trig_eff_mb.C prints the record -- per-pT-bin and per-y-bin tables
# (den/num/eff for data and MC, the data/MC SF), the inclusive pT > 25
# numbers, the MB-prescale checks and the per-run data table. ROOT writes all
# of it to stdout only, hence this wrapper: same convention as
# run_charge_flip.sh, one log per flavor.
#
# Usage:  ./run_trig_eff_mb.sh [mu|ele|both]      (default: both)
#
# Outputs (since 2026-10-20) under <base> = muon_sf_check/trigger (mu) or
# electron_sf_check/trigger (ele), = trig_eff_mb.C's TrigOutBase:
#          <base>/logs/trig_eff_mb_<chan>.log       full console output (THE record)
#          <base>/rootfile/trig_eff_mb_<chan>.root  count histos, TEfficiency objects, SF graphs
#                                                   (read by skim/muon_sf.h, skim/electron_sf.h)
#          <base>/plots/                            turn-on + control plots, trig_eff_<sel>.csv
# The build log stays in correction/logs/trig_eff_mb_build.log (one macro, both flavours).
#
# bash-3.2 safe (macOS stock /bin/bash): no associative arrays, no mapfile.
set -euo pipefail

cd "$(dirname "$0")"

WHICH="${1:-both}"
case "$WHICH" in
  mu|ele|both) ;;
  *) echo "usage: $0 [mu|ele|both]" >&2; exit 2 ;;
esac

base_of() {  # $1 = mu | ele
  if [ "$1" = "mu" ]; then echo muon_sf_check/trigger; else echo electron_sf_check/trigger; fi
}
mkdir -p logs "$(base_of mu)/logs" "$(base_of ele)/logs"

# Pre-build once so the two jobs never race on the ACLiC artifacts.
echo "[BUILD] compiling trig_eff_mb.C ..."
if ! root -l -b -q -e '.L trig_eff_mb.C+' > logs/trig_eff_mb_build.log 2>&1; then
  echo "[FAIL] compilation failed -- see logs/trig_eff_mb_build.log" >&2
  exit 1
fi

run_one() {  # $1 = mu | ele
  chan="$1"
  log="$(base_of "$chan")/logs/trig_eff_mb_${chan}.log"

  printf '[RUN ] %-3s -> %s\n' "$chan" "$log"
  SECONDS=0
  if root -l -b -q "trig_eff_mb.C+(\"${chan}\")" > "$log" 2>&1; then
    status="OK"
  else
    status="FAIL"
  fi
  elapsed=$SECONDS

  # Echo the headline numbers so the terminal still tells you what happened.
  grep -E '^\[CONFIG\]|^\[RESULT\]|^\[INPUT\]|^\[INFO\]|^\[WARN\]|^\[ERR|^\[FATAL' "$log" 2>/dev/null | sed 's/^/       /' || true
  printf '[%-4s] %-3s  %ds   log: %s\n' "$status" "$chan" "$elapsed" "$log"
  [ "$status" = "OK" ] || return 1
}

fail=0
if [ "$WHICH" = "both" ] || [ "$WHICH" = "mu" ];  then run_one mu  || fail=$((fail+1)); fi
if [ "$WHICH" = "both" ] || [ "$WHICH" = "ele" ]; then run_one ele || fail=$((fail+1)); fi

echo "============================================================"
if [ "$fail" -eq 0 ]; then
  echo "[DONE] failures=0   logs in $(pwd)/{muon,electron}_sf_check/trigger/logs/"
else
  echo "[DONE] failures=${fail}  -- check the logs above" >&2
  exit 1
fi
