#!/bin/bash
# skim/run_lhe_updown.sh -- the LHE-weight combination step, with a log.
#
#   ./run_lhe_updown.sh [Wmu|Wel|Zmm|Zee|all]      (default: all)
#
# For every MC skim file of the channel, lhe_updown.py turns the member twins
# written by skim.C into Up/Down templates written back into the same file:
#   <hist>_epps21 -> <hist>_nPDFUp/Down      (LHAPDF PDFSet.uncertainty())
#   <hist>_scale  -> <hist>_qcdScaleUp/Down  (muR/muF envelope over all 9 points, per-bin max/min)
#   <hist>_alphas -> <hist>_alphaSUp/Down    (the 0.119 / 0.117 member templates)
# Must be re-run after every re-skim (the skim RECREATEs its output files).
# Log: logs/lhe_updown_<channel>.log (the record of the integral shifts).
# bash-3.2 safe (macOS /bin/bash).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
# shellcheck source=lhe_env.sh
source ./lhe_env.sh
mkdir -p logs

files_for() {
  local s
  case "$1" in
    Wmu) for s in DY Wp Wm DYtau Wptau Wmtau; do echo "rootfile/WToMuNu_pO_PFMet_${s}_hist.root";   done ;;
    Wel) for s in DY Wp Wm DYtau Wptau Wmtau; do echo "rootfile/WToElecNu_pO_PFMet_${s}_hist.root"; done ;;
    Zmm) for s in DY Wp Wm DYtau Wptau Wmtau; do echo "rootfile/ZToMuMu_pO2025_${s}_MC_hist.root";  done ;;
    Zee) for s in DY Wp Wm DYtau Wptau Wmtau; do echo "rootfile/ZToEE_pO2025_${s}_MC_hist.root";    done ;;
    *)   return 1 ;;
  esac
}

run_channel() {
  local c="$1" log="logs/lhe_updown_$1.log" f present=""
  local list
  list=$(files_for "$c") || { echo "[ERR] unknown channel '$c' (Wmu|Wel|Zmm|Zee|all)"; exit 2; }
  for f in $list; do
    if [ -f "$f" ]; then present="$present $f"; else echo "[WARN] $c: missing $f (not skimmed yet?)"; fi
  done
  [ -n "$present" ] || { echo "[ERR] $c: no MC skim files found in rootfile/"; return 1; }
  echo "=== lhe_updown $c  ->  $log"
  # shellcheck disable=SC2086
  if "$LHE_PYTHON" lhe_updown.py $present >"$log" 2>&1; then
    grep -E "^\[LHAPDF\]|^\[RECIPE\]|^\[SUMMARY\]|^\[WARN\]|^\[ERR\]" "$log" || true
  else
    echo "[ERR] $c: lhe_updown.py failed -- tail of $log:"; tail -20 "$log"; return 1
  fi
}

chan="${1:-all}"
rc=0
if [ "$chan" = "all" ]; then
  for c in Wmu Wel Zmm Zee; do run_channel "$c" || rc=1; done
else
  run_channel "$chan" || rc=1
fi
exit $rc
