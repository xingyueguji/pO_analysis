#!/bin/bash
# merge_rootfile/run_all_hadd_Aug_20.sh
#
# Merge the Aug-20 EMBEDDED production (POWHEG + Angantyr): run
# hadd_from_list.sh over every Aug_20_*.txt filelist in sequence.
# Aug-20 twin of run_all_hadd_July_29.sh -- same behavior:
#   - AUTODETECTION: the sample set is discovered by globbing Aug_20_*.txt in
#     this directory (no hardcoded list -- write the lists first with
#     ./make_filelists_Aug_20.sh), and a sample whose output .root already
#     exists (non-empty) on EOS is SKIPPED, so an already-merged sample is
#     never redone. Re-merge everything with FORCE=1.
#   - a failed sample does NOT abort the run: failures are collected and
#     reported at the end (exit 1 if any), so one flaky EOS file cannot kill
#     the whole chain;
#   - per-sample logs under ./logs_hadd/ (one summary line per sample on
#     stdout, full hadd output in the log).
#
# Output naming = <list basename>.root, i.e. exactly the manual command
#   ./hadd_from_list.sh Aug_20_MC_Wp_mu.txt $OUT/Aug_20_MC_Wp_mu.root
# repeated for every discovered list.
#
# NB this production is MC ONLY -- the data stays on July-29.
#
# Usage:
#   ./run_all_hadd_Aug_20.sh [njobs]        # njobs -> hadd -j, default 4
#   FORCE=1 ./run_all_hadd_Aug_20.sh        # re-merge even if output exists
#
# NB bash-3.2 safe per the repo convention, though this script only makes
# sense where /eos is mounted (lxplus).

set -u   # deliberately NOT -e: continue past per-sample failures

OUT=/eos/cms/store/group/phys_heavyions/zheng/pO_2026_Aug_20
NJOBS="${1:-4}"
FORCE="${FORCE:-0}"

cd "$(dirname "$0")"
mkdir -p "$OUT"          # use 'eos mkdir -p "$OUT"' if plain mkdir fails on EOS
mkdir -p logs_hadd

# --- autodetect the sample lists ---
shopt -s nullglob
LISTS=(Aug_20_*.txt)
shopt -u nullglob
if [[ ${#LISTS[@]} -eq 0 ]]; then
    echo "ERROR: no Aug_20_*.txt filelists found in $(pwd)"
    echo "       run ./make_filelists_Aug_20.sh first"
    exit 1
fi

FAILED=""
NOK=0
NSKIP=0
TOTAL=${#LISTS[@]}
I=0
START=$(date +%s)

for list in "${LISTS[@]}"; do
    I=$((I + 1))
    s="${list%.txt}"
    outroot="$OUT/$s.root"
    printf '=== [%2d/%2d] %s ===\n' "$I" "$TOTAL" "$s"

    if [[ "$FORCE" != "1" && -s "$outroot" ]]; then
        size=$(ls -l "$outroot" | awk '{print $5}')
        echo "    already merged: $outroot ($size bytes) -- skipped (FORCE=1 to redo)"
        NSKIP=$((NSKIP + 1))
        continue
    fi

    if ./hadd_from_list.sh "$list" "$outroot" "$NJOBS" > "logs_hadd/$s.log" 2>&1; then
        NOK=$((NOK + 1))
        tail -1 "logs_hadd/$s.log"   # the "Done. Merged N files -> ..." line
    else
        echo "    FAILED (see logs_hadd/$s.log)"
        FAILED="$FAILED $s"
    fi
done

echo "==================================================="
echo "Done: $NOK merged, $NSKIP already there, $TOTAL total in $(( $(date +%s) - START ))s -> $OUT"
if [[ -n "$FAILED" ]]; then
    echo "FAILED:$FAILED"
    exit 1
fi
