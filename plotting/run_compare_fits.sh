#!/bin/bash
# plotting/run_compare_fits.sh -- two fits of the same data side by side, with a log (2026-10-19).
#
#   ./run_compare_fits.sh [ref_tag] [alt_tag]     (default: leppt_mt40 leppt_mt40_abcdy)
#
# Runs compare_fits.C: fit status + Asimov closure, the SR multijet per rapidity
# group (prefit/postfit from each fit's fitDiagnostics), every nuisance pull, r
# per bin, sigma W+/W-/W/Z, the e/mu ratios, A_ch and R_FB (fiducial) for the
# grand, mu-only and e-only fits in both binnings. The tags are downstream tags
# of plotting/disc_variants.h; run analysis/run_observables.sh <tag> for both
# first (the sigma CSVs and the A_ch / R_FB graphs come from there).
# logs/compare_fits_<ref>_vs_<alt>.log is THE record; the [HEAD] lines are echoed.
# Outputs: plots/compare_fits/<ref>_vs_<alt>/. FORK_TEST overrides the fork
# location (exported absolute, as in run_observables.sh). bash-3.2 safe.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
mkdir -p logs
ref="${1:-leppt_mt40}"
alt="${2:-leppt_mt40_abcdy}"
FORK_TEST="${FORK_TEST:-../../HiggsAnalysis-CombinedLimit/test}"
if [ -d "$FORK_TEST" ]; then FORK_TEST="$(cd "$FORK_TEST" && pwd)"; fi
export FORK_TEST
# pre-build once (a broken build must not look like an empty comparison)
root -l -b -q -e '.L compare_fits.C+' >/dev/null 2>&1 || { echo "[ERR] compare_fits.C does not compile"; root -l -b -q -e '.L compare_fits.C+' 2>&1 | tail -20; exit 1; }
log="logs/compare_fits_${ref}_vs_${alt}.log"
echo "=== compare_fits($ref, $alt)  ->  $log"
root -l -b -q "compare_fits.C+(\"$ref\",\"$alt\")" >"$log" 2>&1
csv="plots/compare_fits/${ref}_vs_${alt}/compare_fits.csv"
if ! grep -q "^\[DONE\]" "$log" || [ ! -s "$csv" ]; then
  echo "[ERR] compare_fits($ref, $alt) did not finish -- tail of $log:"; tail -20 "$log"; exit 1
fi
grep -E "^\[HEAD\]|^\[WARN\]|^\[ERR" "$log" || true
echo "[OK] $log, $csv"
