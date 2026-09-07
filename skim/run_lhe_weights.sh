#!/usr/bin/env bash
# skim/run_lhe_weights.sh
#
# Decode the LHE systematic weights stored in hiEvtAnalyzer/HiTree::ttbar_w
# (217 per event: nominal, QCD scale variations, alternative PDF centrals, two
# nPDF error-set blocks incl. EPPS21) over ALL events of the W+/W-/DY MC, and
# report the relative change of the inclusive cross section per weight index
# plus per-block summaries (scale envelope, Hessian uncertainties, W+/W- ratio,
# charge asymmetry). See the header of lhe_weights.C for the decoded layout.
#
# Outputs:
#   output/lhe_weights.txt             -- per-index table + block summary
#   output/lhe_weights_structure.png   -- slide-ready figure
#   rootfile/lhe_weights.root          -- h_wrel_/h_wrms_/h_sumw_<label>
#   logs/lhe_weights.log               -- full console output (the headline
#                                         block summaries are echoed here)
#
# Usage:
#   ./run_lhe_weights.sh                 # Wp, Wm, DY -- muon-channel files
#   ./run_lhe_weights.sh ele             # the electron-channel files instead
#   ./run_lhe_weights.sh mu 200000       # quick look on the first 200k events
#   ./run_lhe_weights.sh draw            # redraw the figure from rootfile/lhe_weights.root (seconds)
#
# Input paths come from ResolveMCSample (skim_common.h); ~5 min for the three
# full files (the ttbar_w branch is ~0.6 GB compressed per file).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

flavour="${1:-mu}"
nmax="${2:--1}"

mkdir -p rootfile output logs

if [[ "$flavour" == "draw" ]]; then
  log="logs/lhe_weights_draw.log"
  if root -l -q -b 'lhe_weights.C+("draw")' >"$log" 2>&1; then
    echo "[OK] redrew output/lhe_weights_structure.png + _named.png from rootfile/lhe_weights.root   (log: $log)"
    grep "regions" "$log" || true
    exit 0
  else
    echo "[FAIL] lhe_weights.C(\"draw\") returned non-zero -- see $log"
    tail -20 "$log"
    exit 1
  fi
fi

log="logs/lhe_weights_${flavour}.log"

echo "============================================================"
echo "Decoding HiTree::ttbar_w (LHE weights) -- flavour=${flavour} nmax=${nmax}"
echo "============================================================"

if root -l -q -b "lhe_weights.C+(\"Wp,Wm,DY\", \"${flavour}\", ${nmax})" >"$log" 2>&1; then
  echo "[OK] see output/lhe_weights*.txt, output/lhe_weights_structure*.png, rootfile/lhe_weights.root   (log: $log)"
  echo
  sed -n '/=== ttbar_w headline/,$p' "$log"
  exit 0
else
  echo "[FAIL] lhe_weights.C returned non-zero -- see $log"
  tail -20 "$log"
  exit 1
fi
