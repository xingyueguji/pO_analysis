#!/bin/bash
# merge_rootfile/make_filelists_Oct_22.sh
#
# Build the per-sample filelists for the Oct-22 production
# (HiForest 2026_10_22, announced 2026-10-25):
#   - POWHEG MC remade with the fixed boost for asymmetric beams (the Aug-20
#     hard process sat in the symmetric NN CM frame -- CLAUDE.md TODO
#     "MC rapidity boost");
#   - NEW MiNNLO W+jets / Z+jets samples;
#   - ggNtuple extended with supercluster information;
#   - electron MVA ID / isolation replaced by the pp Run-3 ones;
#   - DATA re-made too (unlike Aug-20, which was MC only).
#
# Twin of make_filelists_Aug_20.sh with two additions: the DATA list, and a
# report of every HiForest_* directory of the production that no table entry
# claims, with a suggested table line -- the MiNNLO directory names were not
# known when this script was written, so the FIRST run tells you what to add.
#
# Output (one list per sample, read by run_all_hadd_Oct_22.sh; the names are
# the Aug-20 ones with the date token swapped):
#     Oct_22_DATA.txt                      every IonPhysics* PD, ONE list
#     Oct_22_MC_DY_{mu,ele,tau}_Z.txt
#     Oct_22_MC_W{p,m}_{mu,ele,tau}.txt
#     Oct_22_MC_TT.txt
#     Oct_22_MC_MiNNLO_<name>.txt          once its table line is filled in
#
# Low-mass DY (DYTo*_M_10_50) is ignored by design, as in every production.
#
# Usage (on lxplus, where /eos is mounted):
#   ./make_filelists_Oct_22.sh                       # default production dir
#   ./make_filelists_Oct_22.sh <production_root>     # override
#   PREFIX=X_ ./make_filelists_Oct_22.sh             # MC naming (default Oct_22_MC_)
#
# Exit status: 0 = every table entry found AND no unclaimed directory;
#              1 = a table entry missing/ambiguous/empty, or unclaimed
#                  directories reported (lists are still written for the rest).
#
# NB bash-3.2 safe per the repo convention (no associative arrays).

set -u

PROD_ROOT="${1:-/eos/cms/store/group/phys_heavyions/anstahll/CERN/pO2025/HiForest/2026_10_22}"
PREFIX="${PREFIX:-Oct_22_MC_}"
DATA_OUT="Oct_22_DATA.txt"

cd "$(dirname "$0")"

if [[ ! -d "$PROD_ROOT" ]]; then
    echo "ERROR: production directory does not exist: $PROD_ROOT"
    exit 1
fi

# --- MC sample map: "<generator key>|<EOS dir token>|<output name after PREFIX>"
# A directory matches when its basename is  HiForest_<token>_*<key>*  and it
# sits under $PROD_ROOT/MC/. The key keeps a POWHEG token from also matching a
# MiNNLO directory of the same process (or vice versa). The DY tokens carry the
# explicit _M_50 so they cannot match HiForest_DYTo*_M_10_50_*.
SAMPLES="
POWHEG|DYToMuMu_M_50|DY_mu_Z
POWHEG|DYToEE_M_50|DY_ele_Z
POWHEG|DYToTauTau_M_50|DY_tau_Z
POWHEG|WpToMuNu|Wp_mu
POWHEG|WpToENu|Wp_ele
POWHEG|WpToTauNu|Wp_tau
POWHEG|WmToMuNu|Wm_mu
POWHEG|WmToENu|Wm_ele
POWHEG|WmToTauNu|Wm_tau
POWHEG|TT|TT
"
# MiNNLO W+jets / Z+jets: add one line per sample once the first run has
# printed the directory names, e.g.
#   MiNNLO|<token>|MiNNLO_Wp_mu
# (the suggested line is printed for every unclaimed directory below).
SAMPLES_MINNLO="
"

# Directory basenames that are ignored by design (glob patterns).
IGNORE="
HiForest_DYTo*_M_10_50_*
"

# GNU find (-printf) gives file sizes for free; BSD find does not.
HAVE_PRINTF=0
if find . -maxdepth 0 -printf '' >/dev/null 2>&1; then HAVE_PRINTF=1; fi

size_of() {  # <dir> -> total size of its ROOT files (failed/ and log/ pruned)
    if [[ "$HAVE_PRINTF" -eq 1 ]]; then
        find "$1/" -type d \( -name failed -o -name log \) -prune \
             -o -type f -name '*.root' -printf '%s\n' 2>/dev/null \
        | awk '{s+=$1} END {printf "%.2f GB", s/1024/1024/1024}'
    else
        echo "n/a"
    fi
}

list_root_files() {  # <dir> -> sorted ROOT files, CRAB failed/ and log/ pruned
    find "$1/" -type d \( -name failed -o -name log \) -prune \
         -o -type f -name '*.root' -print 2>/dev/null | sort
}

echo "==================================================="
echo "Production : $PROD_ROOT"
echo "MC prefix  : $PREFIX"
echo "==================================================="

FAILED=""
CLAIMED=""     # newline-separated list of the directories used
NTOT=0
NSAMP=0

# ============================================================== DATA
# Layout (as in every pass so far):
#   DATA/IonPhysics<N>/HiForest_IonPhysics<N>_pORun2025_..._<date>/<crab>/<000X>/*.root
# All PDs go into ONE list. File basenames repeat ACROSS PDs (every PD has a
# HiForestMiniAOD_10.root), so the duplicate check is per PD: two CRAB task
# directories of the same PD = a resubmission that would be merged twice.
echo "--- DATA"
shopt -s nullglob
ddirs=("$PROD_ROOT"/DATA/*/HiForest_*)
shopt -u nullglob
if [[ ${#ddirs[@]} -eq 0 ]]; then
    echo "!! DATA -> no directory DATA/*/HiForest_* under the production"
    FAILED="$FAILED DATA"
else
    : > "$DATA_OUT"
    npd=0
    for d in "${ddirs[@]}"; do
        npd=$((npd + 1))
        CLAIMED="$CLAIMED
$d"
        list_root_files "$d" >> "$DATA_OUT"
        ntask=$(find "$d" -mindepth 1 -maxdepth 1 -type d ! -name failed ! -name log | wc -l | tr -d ' ')
        nf=$(list_root_files "$d" | wc -l | tr -d ' ')
        ndup=$(list_root_files "$d" | awk -F/ '{print $NF}' | sort | uniq -d | wc -l | tr -d ' ')
        printf '    %-62s %5d files  %d CRAB task(s)' "$(basename "$d")" "$nf" "$ntask"
        [[ "$ntask" -gt 1 ]] && printf '  <- more than one task'
        [[ "$ndup" -ne 0 ]] && printf '  WARNING: %d duplicated basename(s) = double counting' "$ndup"
        printf '\n'
    done
    sort -o "$DATA_OUT" "$DATA_OUT"
    n=$(wc -l < "$DATA_OUT" | tr -d ' ')
    printf '%-22s -> %-28s %5d files  %10s   [%d PDs]\n' "DATA" "$DATA_OUT" "$n" \
           "$(size_of "$PROD_ROOT/DATA")" "$npd"
    if [[ "$n" -eq 0 ]]; then
        echo "!! DATA -> 0 ROOT files"
        FAILED="$FAILED DATA"
    fi
    NTOT=$((NTOT + n))
fi

# ============================================================== MC
echo "--- MC"
for pair in $SAMPLES $SAMPLES_MINNLO; do
    gen="${pair%%|*}"
    rest="${pair#*|}"
    tok="${rest%%|*}"
    out="${PREFIX}${rest##*|}.txt"
    NSAMP=$((NSAMP + 1))

    # --- locate the sample directory (must be unique) ---
    shopt -s nullglob
    dirs=("$PROD_ROOT"/MC/*/HiForest_"$tok"_*"$gen"*)
    shopt -u nullglob

    if [[ ${#dirs[@]} -eq 0 ]]; then
        echo "!! $gen $tok -> no directory MC/*/HiForest_${tok}_*${gen}*"
        FAILED="$FAILED $tok"
        continue
    fi
    if [[ ${#dirs[@]} -gt 1 ]]; then
        echo "!! $gen $tok -> AMBIGUOUS, ${#dirs[@]} directories match:"
        for d in "${dirs[@]}"; do echo "     ${d#$PROD_ROOT/}"; done
        FAILED="$FAILED $tok"
        continue
    fi
    sdir="${dirs[0]}"
    CLAIMED="$CLAIMED
$sdir"

    list_root_files "$sdir" > "$out"
    n=$(wc -l < "$out" | tr -d ' ')
    if [[ "$n" -eq 0 ]]; then
        echo "!! $gen $tok -> 0 ROOT files found in $(basename "$sdir")"
        FAILED="$FAILED $tok"
        continue
    fi

    # --- duplicate basenames => a resubmitted CRAB task would be DOUBLE COUNTED
    ndup=$(awk -F/ '{print $NF}' "$out" | sort | uniq -d | wc -l | tr -d ' ')
    subs=$(sed "s|^$sdir/||" "$out" | awk -F/ '{print $1}' | sort -u | tr '\n' ' ')

    printf '%-22s -> %-28s %5d files  %10s   [%s]\n' \
           "$tok" "$out" "$n" "$(size_of "$sdir")" "${subs% }"
    if [[ "$ndup" -ne 0 ]]; then
        echo "   WARNING: $ndup duplicated file basename(s) across the subdirectories above."
        echo "            A resubmitted CRAB task would be merged TWICE -- check before hadd:"
        echo "            awk -F/ '{print \$NF}' $out | sort | uniq -d | head"
    fi
    NTOT=$((NTOT + n))
done

# ============================================================== unclaimed
# Every HiForest_* directory of the production that no entry above used and
# that is not ignored by design -- typically the MiNNLO samples on the first
# run. Printed with a suggested table line.
echo "--- directories not claimed by any table entry"
NUNCL=0
while IFS= read -r d; do
    [[ -z "$d" ]] && continue
    printf '%s\n' "$CLAIMED" | grep -qxF "$d" && continue
    b=$(basename "$d")
    skip=0
    for pat in $IGNORE; do
        # shellcheck disable=SC2053
        [[ "$b" == $pat ]] && skip=1
    done
    if [[ "$skip" -eq 1 ]]; then
        echo "    ignored by design: ${d#$PROD_ROOT/}"
        continue
    fi
    NUNCL=$((NUNCL + 1))
    tok=$(echo "$b" | sed -e 's/^HiForest_//' -e 's/_.*//')
    echo "    UNCLAIMED: ${d#$PROD_ROOT/}   ($(size_of "$d"))"
    echo "               suggested line: <generator key>|$tok|<output name>"
done < <(find "$PROD_ROOT" -mindepth 2 -maxdepth 4 -type d -name 'HiForest_*' 2>/dev/null | sort)
[[ "$NUNCL" -eq 0 ]] && echo "    none"

echo "==================================================="
nmcfail=$(echo $FAILED | tr ' ' '\n' | grep -v -x DATA | grep -c . )
dataok=ok
[[ " $FAILED " == *" DATA "* ]] && dataok=MISSING
echo "MC lists for $((NSAMP - nmcfail)) / $NSAMP table entries, DATA $dataok, $NTOT ROOT files total."
status=0
if [[ -n "$FAILED" ]]; then
    echo "FAILED:$FAILED"
    status=1
fi
if [[ "$NUNCL" -ne 0 ]]; then
    echo "$NUNCL unclaimed director(y/ies): add a SAMPLES_MINNLO line for each sample to use, or an IGNORE pattern."
    status=1
fi
[[ "$status" -eq 0 ]] && echo "Next: ./run_all_hadd_Oct_22.sh"
exit "$status"
