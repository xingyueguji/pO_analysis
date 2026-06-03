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
#   DATA_FILE=root://.../other.root ./run_all.sh Wmu Data   # override data file
#
# Input files (single source of truth): the data-file path is read from
# kDefaultDataFile in skim_common.h, and the MC-file paths from ResolveMCSample
# in the same header. This script does not hardcode any input path.
#
# Implementation note (post-refactor): all four channels are dispatched
# through the unified skim/skim.C macro:
#   root -l -q -b 'skim.C+(<channel-enum>, "<data-file>", <sample-enum>)'
#
# Per job the script reports the input ROOT file that skim.C actually opened
# (printed by skim.C as a "[INPUT] <path>" line and grepped back out of the
# log) so the data/MC file used can be double-checked. The MC file paths are
# resolved inside skim_common.h (ResolveMCSample), which stays the single
# source of truth -- this script never re-derives them.

set -euo pipefail

declare -A CHANNEL_ENUM=(
  [Zmm]=kZmm
  [Zee]=kZee
  [Wmu]=kWmu
  [Wel]=kWel
)

declare -A SAMPLE_ENUM=(
  [Data]=kData
  [DY]=kDY    [Wp]=kWp       [Wm]=kWm
  [DYtau]=kDYtau [Wptau]=kWptau [Wmtau]=kWmtau
)

ALL_CHANNELS=(Zmm Zee Wmu Wel)
ALL_SAMPLES=(Data DY Wp Wm DYtau Wptau Wmtau)

# Data-file path: the single source of truth is `kDefaultDataFile` in
# skim_common.h (which also holds the MC paths via ResolveMCSample). Derive it
# from there so the two never drift. Override for a one-off run by exporting
# DATA_FILE before invoking this script.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMMON_H="${SCRIPT_DIR}/skim_common.h"
if [[ -z "${DATA_FILE:-}" ]]; then
  # `|| true`: keep a failed grep (e.g. header missing) from tripping `set -e`
  # before the guard below can print a useful message.
  DATA_FILE=$(grep -A2 'kDefaultDataFile' "$COMMON_H" 2>/dev/null \
              | grep -oE 'root://[^"]+' | head -1 || true)
fi
[[ -n "${DATA_FILE:-}" ]] || {
  echo "[ERR] could not read the data file (kDefaultDataFile) from $COMMON_H;"
  echo "      set the DATA_FILE env var to override."
  exit 3
}

CHANNEL_ARG="${1:?usage: $0 <Zmm|Zee|Wmu|Wel|all> [samples]}"
shift || true
# Accept samples either as separate args (Data DY Wp) or as a single quoted,
# space-separated string ("Data DY Wp") -- both expand to the same list.
if (( $# > 0 )); then
  read -ra SAMPLES <<< "$*"
else
  SAMPLES=("${ALL_SAMPLES[@]}")
fi

if [[ "$CHANNEL_ARG" == "all" ]]; then
  CHANNELS=("${ALL_CHANNELS[@]}")
else
  CHANNELS=("$CHANNEL_ARG")
fi

# --- Validate up front so the job total is correct before we start ---------
for ch in "${CHANNELS[@]}"; do
  [[ -n "${CHANNEL_ENUM[$ch]:-}" ]] || { echo "[ERR] unknown channel '$ch'"; exit 2; }
done
for s in "${SAMPLES[@]}"; do
  [[ -n "${SAMPLE_ENUM[$s]:-}" ]] || { echo "[ERR] unknown sample '$s'"; exit 2; }
done

mkdir -p rootfile output logs

# --- Progress-bar helpers --------------------------------------------------
# Use the live in-place bar only on a real terminal; when stdout is piped or
# redirected, fall back to plain scrolling lines (no control characters).
if [[ -t 1 ]]; then USE_BAR=1; else USE_BAR=0; fi
BAR_WIDTH=24

bar_str() { # $1=filled $2=empty -> "######------"
  local f e
  f=$(printf '%*s' "$1" '' | tr ' ' '#')
  e=$(printf '%*s' "$2" '' | tr ' ' '-')
  printf '%s%s' "$f" "$e"
}

draw_bar() { # $1=done $2=total $3=label ; in-place, no trailing newline
  (( USE_BAR )) || return 0
  local done=$1 total=$2 label=${3:-}
  local pct=0 filled=0 empty
  if (( total > 0 )); then
    pct=$(( done * 100 / total ))
    filled=$(( done * BAR_WIDTH / total ))
  fi
  empty=$(( BAR_WIDTH - filled ))
  printf '\r\033[K[%s] %3d%%  %2d/%-2d  %s' \
         "$(bar_str "$filled" "$empty")" "$pct" "$done" "$total" "$label"
}

record() { # print a full line; clears the pinned bar first when it is live
  if (( USE_BAR )); then printf '\r\033[K%s\n' "$*"; else printf '%s\n' "$*"; fi
}

# --- Run -------------------------------------------------------------------
total=$(( ${#CHANNELS[@]} * ${#SAMPLES[@]} ))
completed=0
fail=0

echo "============================================================"
printf "Skim run: %d channel(s) x %d sample(s) = %d job(s)\n" \
       "${#CHANNELS[@]}" "${#SAMPLES[@]}" "$total"
echo "============================================================"

for ch in "${CHANNELS[@]}"; do
  chEnum="${CHANNEL_ENUM[$ch]}"
  for s in "${SAMPLES[@]}"; do
    sEnum="${SAMPLE_ENUM[$s]}"
    log="logs/skim_${ch}_${s}.log"

    draw_bar "$completed" "$total" "running ${ch}/${s} ..."

    SECONDS=0
    if root -l -q -b "skim.C+(${chEnum}, \"${DATA_FILE}\", ${sEnum})" >"$log" 2>&1; then
      status="OK"
    else
      status="FAIL"
      fail=$(( fail + 1 ))
    fi
    elapsed=$SECONDS

    # The file skim.C actually opened (its own "[INPUT] ..." line).
    infile=$(grep -m1 '^\[INPUT\] ' "$log" 2>/dev/null | sed 's/^\[INPUT\] //' || true)
    [[ -n "$infile" ]] || infile="(no [INPUT] line -- see $log)"

    completed=$(( completed + 1 ))
    record "$(printf '[%2d/%2d] %-3s / %-6s  %-4s %4ds   log: %s' \
              "$completed" "$total" "$ch" "$s" "$status" "$elapsed" "$log")"
    record "         input: $infile"

    draw_bar "$completed" "$total" "${status}  ${ch}/${s}"
  done
done

(( USE_BAR )) && printf '\n'

echo
echo "Done. ${completed}/${total} jobs, failures=${fail}"
exit "$fail"
