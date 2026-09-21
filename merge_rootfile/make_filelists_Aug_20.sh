#!/bin/bash
# merge_rootfile/make_filelists_Aug_20.sh
#
# Build the per-sample filelists for the Aug-20 EMBEDDED production
# (POWHEG + Angantyr underlying event, HiForest 2026_08_20).
#
# This is the per-sample, self-driving twin of make_filelist.sh: instead of
# being pointed at ONE sample directory and one output .txt, it walks the whole
# production directory and writes one list per sample, named so that the merge
# driver (run_all_hadd_Aug_20.sh) and skim/skim_common.h::ResolveMCSample can
# consume them with no further translation:
#
#     Aug_20_MC_DY_{mu,ele,tau}_Z.txt      ->  Aug_20_MC_DY_{mu,ele,tau}_Z.root
#     Aug_20_MC_W{p,m}_{mu,ele,tau}.txt    ->  Aug_20_MC_W{p,m}_{mu,ele,tau}.root
#
# (= the July-29 naming with the date token swapped, so LabelFromFname in
# skim/count_ngen.C only needs "Aug_20_MC_" added to its prefix-strip list.)
#
# The low-mass DY samples (DYTo*_M_10_50) are present on EOS but deliberately
# NOT listed here -- the analysis has never used them (user decision, and the
# July-29 production had none either).
#
# There is no DATA in this production directory; data stays on July-29.
#
# Usage:
#   ./make_filelists_Aug_20.sh                       # default production dir
#   ./make_filelists_Aug_20.sh <production_dir>      # override
#   PREFIX=Emb_Aug_20_MC_ ./make_filelists_Aug_20.sh # different output naming
#
# NB bash-3.2 safe per the repo convention (no associative arrays), though it
# only makes sense where /eos is mounted (lxplus).

set -u

PROD_DIR="${1:-/eos/cms/store/group/phys_heavyions/anstahll/CERN/pO2025/HiForest/2026_08_20/MC/POWHEG_9p62TeV_2025Run3}"
PREFIX="${PREFIX:-Aug_20_MC_}"

cd "$(dirname "$0")"

if [[ ! -d "$PROD_DIR" ]]; then
    echo "ERROR: production directory does not exist: $PROD_DIR"
    exit 1
fi

# --- sample map: "<EOS dir token>|<output basename after PREFIX>" ------------
# The dir token is matched as  HiForest_<token>_*  , which is why the DY
# entries carry the explicit _M_50: "DYToEE_M_50" cannot match the
# "HiForest_DYToEE_M_10_50_..." directory.
SAMPLES="
DYToMuMu_M_50|DY_mu_Z
DYToEE_M_50|DY_ele_Z
DYToTauTau_M_50|DY_tau_Z
WpToMuNu|Wp_mu
WpToENu|Wp_ele
WpToTauNu|Wp_tau
WmToMuNu|Wm_mu
WmToENu|Wm_ele
WmToTauNu|Wm_tau
"

# GNU find (-printf) gives us file sizes for free; BSD find does not.
HAVE_PRINTF=0
if find . -maxdepth 0 -printf '' >/dev/null 2>&1; then HAVE_PRINTF=1; fi

echo "==================================================="
echo "Production : $PROD_DIR"
echo "Prefix     : $PREFIX"
echo "==================================================="

FAILED=""
NTOT=0
NSAMP=0

for pair in $SAMPLES; do
    tok="${pair%%|*}"
    out="${PREFIX}${pair##*|}.txt"
    NSAMP=$((NSAMP + 1))

    # --- locate the sample directory (must be unique) ---
    shopt -s nullglob
    dirs=("$PROD_DIR"/HiForest_"$tok"_*)
    shopt -u nullglob

    if [[ ${#dirs[@]} -eq 0 ]]; then
        echo "!! $tok -> no directory HiForest_${tok}_* under the production dir"
        FAILED="$FAILED $tok"
        continue
    fi
    if [[ ${#dirs[@]} -gt 1 ]]; then
        echo "!! $tok -> AMBIGUOUS, ${#dirs[@]} directories match HiForest_${tok}_*:"
        for d in "${dirs[@]}"; do echo "     $(basename "$d")"; done
        FAILED="$FAILED $tok"
        continue
    fi
    sdir="${dirs[0]}"

    # --- collect the ROOT files (skip CRAB failed/ and log/ trees) ---
    find "$sdir/" -type d \( -name failed -o -name log \) -prune \
         -o -type f -name '*.root' -print 2>/dev/null | sort > "$out"

    n=$(wc -l < "$out" | tr -d ' ')
    if [[ "$n" -eq 0 ]]; then
        echo "!! $tok -> 0 ROOT files found in $(basename "$sdir")"
        FAILED="$FAILED $tok"
        continue
    fi

    # --- duplicate basenames => a resubmitted CRAB task would be DOUBLE COUNTED
    ndup=$(awk -F/ '{print $NF}' "$out" | sort | uniq -d | wc -l | tr -d ' ')

    # --- which CRAB submission subdirectories contributed ---
    subs=$(sed "s|^$sdir/||" "$out" | awk -F/ '{print $1}' | sort -u | tr '\n' ' ')

    if [[ "$HAVE_PRINTF" -eq 1 ]]; then
        size=$(find "$sdir/" -type d \( -name failed -o -name log \) -prune \
                    -o -type f -name '*.root' -printf '%s\n' 2>/dev/null \
               | awk '{s+=$1} END {printf "%.2f GB", s/1024/1024/1024}')
    else
        size="n/a"
    fi

    printf '%-22s -> %-28s %5d files  %10s   [%s]\n' \
           "$tok" "$out" "$n" "$size" "${subs% }"
    if [[ "$ndup" -ne 0 ]]; then
        echo "   WARNING: $ndup duplicated file basename(s) across the subdirectories above."
        echo "            A resubmitted CRAB task would be merged TWICE -- check before hadd:"
        echo "            awk -F/ '{print \$NF}' $out | sort | uniq -d | head"
    fi
    NTOT=$((NTOT + n))
done

echo "==================================================="
echo "Wrote lists for $((NSAMP - $(echo $FAILED | wc -w | tr -d ' '))) / $NSAMP samples, $NTOT ROOT files total."
if [[ -n "$FAILED" ]]; then
    echo "FAILED:$FAILED"
    exit 1
fi
echo "Next: ./run_all_hadd_Aug_20.sh"
