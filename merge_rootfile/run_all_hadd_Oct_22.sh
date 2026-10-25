#!/bin/bash
# merge_rootfile/run_all_hadd_Oct_22.sh
#
# Merge the Oct-22 production (HiForest 2026_10_22: boost-fixed POWHEG, new
# MiNNLO W/Z+jets, supercluster branches, pp Run-3 electron MVA, DATA re-made):
# run hadd_from_list.sh over every Oct_22_*.txt filelist in sequence.
# Oct-22 twin of run_all_hadd_Aug_20.sh -- same behavior:
#   - AUTODETECTION: the sample set is discovered by globbing Oct_22_*.txt in
#     this directory (write the lists first with ./make_filelists_Oct_22.sh),
#     and a sample whose output .root already exists (non-empty) on EOS is
#     SKIPPED. Re-merge everything with FORCE=1.
#   - a failed sample does NOT abort the run: failures are collected and
#     reported at the end (exit 1 if any);
#   - per-sample logs under ./logs_hadd/.
# New here: the DATA list is merged LAST (one ~50 GB file, the longest job),
# so the MC -- which the boost check needs first -- is ready early; and
# ONLY="Wp_mu Wm_mu DY_mu_Z" restricts the run to the lists whose name
# contains one of the words.
#
# Output naming = <list basename>.root, i.e. exactly the manual command
#   ./hadd_from_list.sh Oct_22_MC_Wp_mu.txt $OUT/Oct_22_MC_Wp_mu.root
#
# Usage (on lxplus):
#   ./run_all_hadd_Oct_22.sh [njobs]                 # njobs -> hadd -j, default 4
#   ONLY="Wp_mu Wm_mu DY_mu_Z" ./run_all_hadd_Oct_22.sh
#   FORCE=1 ./run_all_hadd_Oct_22.sh                 # re-merge even if output exists
#
# NB bash-3.2 safe per the repo convention, though this script only makes
# sense where /eos is mounted (lxplus).

set -u   # deliberately NOT -e: continue past per-sample failures

OUT=/eos/cms/store/group/phys_heavyions/zheng/pO_2026_Oct_22
NJOBS="${1:-4}"
FORCE="${FORCE:-0}"
ONLY="${ONLY:-}"

cd "$(dirname "$0")"
mkdir -p "$OUT"          # use 'eos mkdir -p "$OUT"' if plain mkdir fails on EOS
mkdir -p logs_hadd

# --- autodetect the sample lists: MC first, DATA last ---
shopt -s nullglob
ALL=(Oct_22_MC_*.txt Oct_22_DATA*.txt)
shopt -u nullglob
LISTS=()
for list in "${ALL[@]}"; do
    if [[ -n "$ONLY" ]]; then
        keep=0
        for w in $ONLY; do [[ "$list" == *"$w"* ]] && keep=1; done
        [[ "$keep" -eq 1 ]] || continue
    fi
    LISTS+=("$list")
done
if [[ ${#LISTS[@]} -eq 0 ]]; then
    echo "ERROR: no Oct_22_*.txt filelists found in $(pwd)${ONLY:+ matching ONLY=\"$ONLY\"}"
    echo "       run ./make_filelists_Oct_22.sh first"
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
echo "Next (on the Mac): merge_rootfile/download_Oct_22.sh"
