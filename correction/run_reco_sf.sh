#!/usr/bin/env bash
# correction/run_reco_sf.sh
#
# Run the ELECTRON RECO-EFFICIENCY TAG-AND-PROBE in pO (2026-09-29) and KEEP
# THE LOGS. The definition lives in the header of reco_sf_ele.C: Z->ee pairs of
# a WP90 ID+iso, EG10-matched tag and a supercluster probe (a photon object, NO
# ID), pass = a forest electron on the probe's SC; eps_data from a pass/fail
# mass fit, eps_MC from the DY->ee MC, SF compared with the pp EGM Reco SF.
#
# Usage:  ./run_reco_sf.sh [skim|fit|all]      (default: all)
#   skim  the ntuple loops, Data and DY          -> rootfile/reco_sf_ele/skim_{Data,DY}.root
#   fit   the pass/fail fits, SF, EGM comparison -> plots/reco_sf_ele/fit_<variant>_<bin>_{pass,fail}.{png,pdf},
#                                                   reco_sf_vs_{eta,et}.{png,pdf} + reco_sf_ele.csv
#                                                   + rootfile/reco_sf_ele/reco_sf_ele_result.root = the
#                                                   pO Reco SF, applied by the skim under
#                                                   kEleRecoSource = kRecoTnP only (2026-09-30; since
#                                                   2026-10-01 the pp EGM tables are applied and this is
#                                                   the cross-check) -- under kRecoTnP re-skim the
#                                                   electron MC after a re-run (skim/run_all.sh Wel, Zee)
# Logs: logs/reco_sf_ele_skim_{Data,DY}.log, logs/reco_sf_ele_fit.log -- ROOT
# prints the record to stdout only.
#
# bash-3.2 safe (macOS stock /bin/bash): no associative arrays, no mapfile.
set -euo pipefail

cd "$(dirname "$0")"

STAGE="${1:-all}"
case "$STAGE" in
  skim|fit|all) ;;
  *) echo "usage: $0 [skim|fit|all]" >&2; exit 2 ;;
esac

mkdir -p logs rootfile/reco_sf_ele plots/reco_sf_ele

# pre-build once so no two jobs race on the ACLiC artifacts
echo "[BUILD] compiling reco_sf_ele.C ..."
if ! root -l -b -q -e ".L reco_sf_ele.C+" > logs/reco_sf_ele_build.log 2>&1; then
  echo "[FAIL] compilation of reco_sf_ele.C failed -- see logs/reco_sf_ele_build.log" >&2
  tail -20 logs/reco_sf_ele_build.log >&2
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
    run_root "skim $s" "reco_sf_ele.C+(\"skim\",\"$s\")" "logs/reco_sf_ele_skim_${s}.log" \
             '^\[INPUT\]|^\[RAW\]|^\[WARN\]|^\[FATAL\]' || fail=$((fail+1))
  done
fi
if [ "$STAGE" = "fit" ] || [ "$STAGE" = "all" ]; then
  if [ "$fail" -gt 0 ]; then
    echo "[SKIP] fit: the skim had failures" >&2
  else
    run_root "fit" "reco_sf_ele.C+(\"fit\")" "logs/reco_sf_ele_fit.log" \
             '^\[RESULT\]|^ *\[CC\]|^ *\[WARN\]|^\[MC-TRUTH\]|^\[OUT\]|^\[FATAL\]' || fail=$((fail+1))
  fi
fi

echo "============================================================"
if [ "$fail" -eq 0 ]; then
  echo "[DONE] failures=0   logs in $(pwd)/logs/"
else
  echo "[DONE] failures=${fail}  -- check the logs above" >&2
  exit 1
fi
