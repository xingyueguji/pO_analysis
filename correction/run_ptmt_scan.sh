#!/usr/bin/env bash
# correction/run_ptmt_scan.sh
#
# Run the (lepton-pT x m_T) cut-pair scan and KEEP THE LOG (2026-10-16; the
# scan's console output -- the normalization inputs, the baseline / tentative
# / optimum table and the per-pT-cut optimum rows -- is what
# docs/AN_selection_optimization.tex quotes, and ROOT writes it to stdout
# only). Same convention as run_qcd_abcd.sh: one log per flavor under
# correction/logs/, a pre-build so two jobs never race on ACLiC artifacts.
#
# Usage:  ./run_ptmt_scan.sh [mu|ele|both]      (default: both)
#
# Inputs:  the W skim outputs ../skim/rootfile/*_hist.root (the h_pt_mt*
#          scan planes) + ../skim/rootfile/ngen.root (k_s) -- re-run after a
#          re-skim or a change of the EWK sample list (ttbar since 2026-10-16).
# Outputs: logs/ptmt_scan_<chan>.log         full console output (THE record)
#          rootfile/ptmt_scan_<chan>.root    maps, sweeps, profiles
#          plots/ptmt_scan_<chan>/           S/sqrt(S+B) & S/B maps, ROCs, m_T profiles
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

echo "[BUILD] compiling ptmt_scan.C ..."
if ! root -l -b -q -e '.L ptmt_scan.C+' > logs/ptmt_scan_build.log 2>&1; then
  echo "[FAIL] compilation failed -- see logs/ptmt_scan_build.log" >&2
  exit 1
fi

run_one() {  # $1 = mu | ele
  chan="$1"
  log="logs/ptmt_scan_${chan}.log"
  if [ "$chan" = "ele" ]; then arg="true"; else arg="false"; fi

  printf '[RUN ] %-3s -> %s\n' "$chan" "$log"
  SECONDS=0
  if root -l -b -q "ptmt_scan.C+(${arg})" > "$log" 2>&1; then
    status="OK"
  else
    status="FAIL"
  fi
  elapsed=$SECONDS

  # Echo the headline numbers: the normalization inputs and the summary rows.
  grep -E '^\[INFO\]|baseline \(current\)|tentative |OPTIMUM|^\[WARN\]|^\[ERR|^\[FATAL' "$log" 2>/dev/null | sed 's/^/       /' || true
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
