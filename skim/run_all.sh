#!/usr/bin/env bash
# skim/run_all.sh
#
# Usage:
#   ./run_all.sh <channel> [samples]
#
# <channel> :  Zmm | Zee | Wmu | Wel | all
# [samples] :  space-separated subset of:
#                Data DY Wp Wm DYtau Wptau Wmtau
#              (default: all 7)
#
# Examples:
#   ./run_all.sh Zmm                  # all 7 samples for Z->mumu
#   ./run_all.sh Wel "Data DY"        # only data + DY MC for W->enu
#   ./run_all.sh all                  # 4 channels x 7 samples = 28 jobs

set -euo pipefail

declare -A MACRO=(
  [Zmm]="DrawDimuonPeak.C"
  [Zee]="DrawDielectronPeak.C"
  [Wmu]="DrawWToMuNu_PFMet.C"
  [Wel]="DrawWToElecNu_PFMet.C"
)

declare -A SAMPLE_ENUM=(
  [Data]=kData
  [DY]=kDY    [Wp]=kWp       [Wm]=kWm
  [DYtau]=kDYtau [Wptau]=kWptau [Wmtau]=kWmtau
)

ALL_CHANNELS=(Zmm Zee Wmu Wel)
ALL_SAMPLES=(Data DY Wp Wm DYtau Wptau Wmtau)

DATA_FILE="root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2025.root"

CHANNEL_ARG="${1:?usage: $0 <Zmm|Zee|Wmu|Wel|all> [samples]}"
shift || true
SAMPLES=("$@")
[[ ${#SAMPLES[@]} -eq 0 ]] && SAMPLES=("${ALL_SAMPLES[@]}")

if [[ "$CHANNEL_ARG" == "all" ]]; then
  CHANNELS=("${ALL_CHANNELS[@]}")
else
  CHANNELS=("$CHANNEL_ARG")
fi

mkdir -p rootfile output logs

fail=0
for ch in "${CHANNELS[@]}"; do
  macro="${MACRO[$ch]:-}"
  [[ -z "$macro" ]] && { echo "[ERR] unknown channel '$ch'"; exit 2; }
  fn="${macro%.C}"
  for s in "${SAMPLES[@]}"; do
    enum="${SAMPLE_ENUM[$s]:-}"
    [[ -z "$enum" ]] && { echo "[ERR] unknown sample '$s'"; exit 2; }
    log="logs/${fn}_${s}.log"
    printf ">>> [%s / %-6s] -> %s\n" "$ch" "$s" "$log"
    if root -l -q -b "${macro}+(\"${DATA_FILE}\", ${enum})" >"$log" 2>&1; then
      echo "    OK"
    else
      echo "    FAIL  (tail $log)"
      fail=$((fail+1))
    fi
  done
done

echo
echo "Done. failures=$fail"
exit $fail