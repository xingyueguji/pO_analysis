#!/bin/bash
# correction/run_isolation.sh
#
# Logging wrapper for the isolation / electron-ID working-point studies, in the
# same spirit as run_qcd_abcd.sh and skim/run_all.sh: ROOT prints the scan
# tables (AUC, Youden J, efficiencies at the analysis cut) to stdout ONLY, and
# those numbers are quoted in AN_selection_optimization.tex -- a bare
# `root -l -q 'isolation_ele.C+'` silently throws the record away.
#
# Usage:
#   ./run_isolation.sh [mu|ele|both]     # default: both
#
# Input file: taken from pOSkim::kDefaultDataFile in skim/skim_common.h, the
# SINGLE SOURCE OF TRUTH for the production in use -- do NOT hardcode a path
# here or in the macros.  Override a single run with
#   DATA_FILE=/some/other.root ./run_isolation.sh mu
#
# Outputs:
#   rootfile/IsoStudyOutputs_muon_tight.root   (muon, dbeta + nodbeta scans)
#   rootfile/IsoStudyOutputs_electron.root     (electron, 8 ID working points)
#   logs/isolation_{mu,ele}.log                (the full scan tables)
#   plotsROC/summary_muon.png, plotsROC_ele/summary_<ID>.png  (via plot_iso_summary.C)
#
# NB these macros loop the FULL ntuple (~3.9M events) -- ~10 min per channel.

set -u

CHAN="${1:-both}"
case "$CHAN" in
  mu|ele|both) ;;
  *) echo "usage: $0 [mu|ele|both]" >&2; exit 2 ;;
esac

cd "$(dirname "$0")" || exit 1
mkdir -p logs rootfile plots plotsROC plotsROC_ele

# ---- report the input actually in use (same idea as run_all.sh's [INPUT]) ----
HEADER="../skim/skim_common.h"
if [ -n "${DATA_FILE:-}" ]; then
  IN="$DATA_FILE"
  ARG="(\"$IN\")"
  echo "[INPUT] $IN   (DATA_FILE override)"
else
  IN=$(grep -A2 'kDefaultDataFile *=' "$HEADER" \
        | grep -oE '"[^"]+\.root"' | head -1 | tr -d '"' || true)
  ARG=""
  echo "[INPUT] ${IN:-<unresolved>}   (from $HEADER)"
fi
if [ -n "$IN" ] && [ ! -f "$IN" ] && [ "${IN#root://}" = "$IN" ]; then
  echo "[WARN] input file not found locally: $IN" >&2
fi

# ---- pre-build once so parallel/repeat runs cannot race on ACLiC artifacts ---
echo "[BUILD] compiling macros ..."
BUILD_OK=1
for m in isolation_mu_tight isolation_ele plot_iso_summary; do
  if ! root -l -b -q -e ".L ${m}.C+" > "logs/build_${m}.log" 2>&1; then
    echo "[ERR] build failed: $m (see logs/build_${m}.log)" >&2
    BUILD_OK=0
  fi
done
[ "$BUILD_OK" -eq 1 ] || exit 1

run_one() {
  # $1 = macro name, $2 = log tag
  echo "[RUN] $1 -> logs/isolation_$2.log"
  root -l -b -q "$1.C+$ARG" 2>&1 | tee "logs/isolation_$2.log" > /dev/null
  echo "----- $2: control samples -----"
  grep -E "^Selected:" "logs/isolation_$2.log" || echo "  (no 'Selected:' line -- run failed?)"
  echo "----- $2: continuous relIso scan -----"
  # the table between the header row and the closing '====' rule
  awk '/aucSS/{p=1} p&&/^=====/{p=0} p' "logs/isolation_$2.log"
}

case "$CHAN" in
  mu)   run_one isolation_mu_tight mu ;;
  ele)  run_one isolation_ele      ele ;;
  both) run_one isolation_mu_tight mu
        run_one isolation_ele      ele ;;
esac

# ---- summary plots (cheap: reads the stored graphs, no ntuple loop) ----------
echo "[PLOT] plot_iso_summary.C -> plotsROC[_ele]/summary_*.png"
root -l -b -q 'plot_iso_summary.C+' > logs/plot_iso_summary.log 2>&1 \
  || echo "[WARN] plot_iso_summary failed (see logs/plot_iso_summary.log)" >&2

echo "[DONE] logs in correction/logs/, graphs in correction/rootfile/"
