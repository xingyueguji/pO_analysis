#!/usr/bin/env bash
# correction/run_tnp_charge_ele.sh
#
# Run the ELECTRON TAG-AND-PROBE SPLIT BY PROBE CHARGE (2026-10-20) and KEEP
# THE LOGS. The definitions live in the header of tnp_charge_ele.C: the reco
# link of reco_sf_ele.C and the Z ID+iso link of idiso_sf_skim.C, with the
# probes split into charge x {|eta_SC| < 1.6, forward eta_SC < -1.6, backward
# eta_SC > +1.6} (y = -eta); the answer is SF(e-)/SF(e+) per |eta| half and per
# side ("[RATIO]" lines). The reco link is run twice on the same probes: with
# the SC-energy probe four-vector (reco, EGM's) and the previous photon-object
# one (recopho); the two reco plots show both, with the pp EGM Reco SF in its
# own (eta_SC, pT) cells, and the ID+iso plots carry the EGM wp90iso cells
# ("[EGM]" lines: which pT cells are drawn, by probe fraction).
#
# Outputs (since 2026-10-20) under correction/electron_sf_check/charge_split/:
# Usage:  ./run_tnp_charge_ele.sh [skim|fit|all]      (default: all)
#   skim  the ntuple loops, Data and DY         -> rootfile/skim_{Data,DY}.root
#   fit   per bin eps (fit + cut-and-count), SF -> plots/fit_<link>_<bin>_{pass,fail},
#         and the e-/e+ ratios                     sf_charge[_signed]_{reco,idiso} (.png/.pdf) + tnp_charge_ele.csv
# The fit stage also checks that the cells add up to the parent studies'
# inclusive numbers ("[CHECK]"): run ./run_reco_sf.sh skim and
# ./run_idiso_sf.sh skim first if their skim files are absent.
# Logs: logs/tnp_charge_ele_skim_{Data,DY}.log, logs/tnp_charge_ele_fit.log --
# ROOT prints the record to stdout only.
#
# bash-3.2 safe (macOS stock /bin/bash): no associative arrays, no mapfile.
set -euo pipefail

cd "$(dirname "$0")"

STAGE="${1:-all}"
case "$STAGE" in
  skim|fit|all) ;;
  *) echo "usage: $0 [skim|fit|all]" >&2; exit 2 ;;
esac

OUT=electron_sf_check/charge_split   # = tnp_charge_ele.C's kOutDir / kPlotDir parent
LOG=$OUT/logs
mkdir -p "$LOG" "$OUT/rootfile" "$OUT/plots"

# pre-build once so no two jobs race on the ACLiC artifacts
echo "[BUILD] compiling tnp_charge_ele.C ..."
if ! root -l -b -q -e ".L tnp_charge_ele.C+" > "$LOG/tnp_charge_ele_build.log" 2>&1; then
  echo "[FAIL] compilation of tnp_charge_ele.C failed -- see $LOG/tnp_charge_ele_build.log" >&2
  tail -20 "$LOG/tnp_charge_ele_build.log" >&2
  exit 1
fi

run_root() {  # $1 = label, $2 = root command, $3 = log, $4 = headline regex
  printf '[RUN ] %-10s -> %s\n' "$1" "$3"
  SECONDS=0
  if root -l -b -q "$2" > "$3" 2>&1; then status="OK"; else status="FAIL"; fi
  # root exits 0 even when the macro bails out -> also check for a FATAL line
  grep -q '^\[FATAL\]' "$3" && status="FAIL"
  grep -E "$4" "$3" 2>/dev/null | sed 's/^/       /' || true
  printf '[%-4s] %-10s %ds   log: %s\n' "$status" "$1" "$SECONDS" "$3"
  [ "$status" = "OK" ]
}

fail=0
if [ "$STAGE" = "skim" ] || [ "$STAGE" = "all" ]; then
  for s in Data DY; do
    run_root "skim $s" "tnp_charge_ele.C+(\"skim\",\"$s\")" "$LOG/tnp_charge_ele_skim_${s}.log" \
             '^\[INPUT\]|^\[RAW\]|^\[WARN\]|^\[FATAL\]' || fail=$((fail+1))
  done
fi
if [ "$STAGE" = "fit" ] || [ "$STAGE" = "all" ]; then
  if [ "$fail" -gt 0 ]; then
    echo "[SKIP] fit: the skim had failures" >&2
  else
    run_root "fit" "tnp_charge_ele.C+(\"fit\")" "$LOG/tnp_charge_ele_fit.log" \
             '^\[CHECK\]|^\[RATIO\]|^\[EGM\] .*-> drawn$|^ *\[WARN\]|^\[OUT\]|^\[FATAL\]' || fail=$((fail+1))
  fi
fi

echo "============================================================"
if [ "$fail" -eq 0 ]; then
  echo "[DONE] failures=0   logs in $(pwd)/$LOG/"
else
  echo "[DONE] failures=${fail}  -- check the logs above" >&2
  exit 1
fi
