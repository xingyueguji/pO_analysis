#!/usr/bin/env bash
# correction/run_qcd_abcd.sh
#
# Run the data-driven QCD (ABCD) background estimate and KEEP THE LOG.
#
# qcd_abcd.C prints far more than the templates it writes: the region
# composition and closure, the factorisation tests (T in sub-slices of the low-y
# band), the anti-iso window scan, the plane-to-plane transport with its r-scan,
# the multijet fraction of each fitted selection, and the assembled systematic
# budget -> lnN kappa. Those numbers are quoted in the AN
# (../AN_qcd_background.tex), so they must be reproducible, not scrolled past.
# ROOT writes them to stdout only, hence this wrapper: same convention as
# skim/run_all.sh, one log per job under correction/logs/.
#
# Usage:  ./run_qcd_abcd.sh [mu|ele|both]      (default: both)
#
# Outputs: logs/qcd_abcd_<chan>.log      full console output (THE record)
#          rootfile/qcd_abcd_<chan>.root  QCD templates for Combine
#          plots/qcd_abcd_<chan>/         2D planes, closure, shape checks
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

# Pre-build once so the two jobs never race on the ACLiC artifacts (the same
# trap as running two skim/run_all.sh invocations from one directory).
echo "[BUILD] compiling qcd_abcd.C ..."
if ! root -l -b -q -e '.L qcd_abcd.C+' > logs/qcd_abcd_build.log 2>&1; then
  echo "[FAIL] compilation failed -- see logs/qcd_abcd_build.log" >&2
  exit 1
fi

run_one() {  # $1 = mu | ele
  chan="$1"
  arg="false"; [ "$chan" = "ele" ] && arg="true"
  log="logs/qcd_abcd_${chan}.log"

  printf '[RUN ] %-3s -> %s\n' "$chan" "$log"
  SECONDS=0
  if root -l -q -b "qcd_abcd.C+(${arg})" > "$log" 2>&1; then
    status="OK"
  else
    status="FAIL"
  fi
  elapsed=$SECONDS

  # Echo the headline numbers so the terminal still tells you what happened.
  grep -E '^\[CONFIG\]|Transfer factor|=> lnN kappa|=> in-fit ABCD|A40 = B\*C40/D' "$log" 2>/dev/null | sed 's/^/       /' || true
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
