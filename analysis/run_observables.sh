#!/bin/bash
# =============================================================================
# run_observables.sh -- Module 5 driver: fitted yields -> final observable
# plots for ONE W-discriminant variant (or all three), carrying the disc tag
# through every filename and output folder so variants never overwrite each
# other (2026-08-03).
#
#   ./run_observables.sh [met|leppt|leppt_mt40|all]        (default: met)
#
# Per discriminant it runs, for each channel (mu, ele):
#   1. charge_asym.C   fitted yields -> ../skim/rootfile/charge_asym_fit_<chan>_<disc>.root
#   2. FBratio.C       fitted yields -> ../skim/rootfile/FBratio_fit_<chan>_<disc>.root
# then, from ../plotting:
#   3. observables(false|true, disc)  -> plots[/Elec]/{charge_asym,FBratio}/<disc>/
#   4. observables_overlay(disc)      -> plots/merged/<disc>/
#   5. xsec_fiducial(disc) + xsec_fiducial_diff(mu|ele, disc)
#                                     -> plots/xsec/<disc>/
#
# Input: the fork's out-tree for that variant must exist (fit already run +
# downloaded):  $FORK_TEST/pO_fit_out<suffix>/<chan>/summary/<chan>_fitted_yields.root
# produced by run_pO_fits.sh both all [--disc <disc>]. With 'all', variants
# whose out-tree is missing are SKIPPED with a warning; a single named variant
# errors out instead.
#
# FORK_TEST env var overrides the fork location (default: sibling checkout).
# bash-3.2 compatible (macOS stock bash): no associative arrays.
# =============================================================================
set -u

HERE="$(cd "$(dirname "$0")" && pwd)"
FORK_TEST="${FORK_TEST:-$HERE/../../HiggsAnalysis-CombinedLimit/test}"

disc_suffix() {
  case "$1" in
    met)        echo "" ;;
    leppt)      echo "_leppt" ;;
    leppt_mt40) echo "_leppt_mt40" ;;
    *)          echo "[ERROR] unknown disc '$1' (met|leppt|leppt_mt40)" >&2; exit 1 ;;
  esac
}

# require_file <path> <what to run if missing>
require_file() {
  if [ ! -f "$1" ]; then
    echo "[ERROR] expected output missing: $1"
    echo "        ($2)"
    exit 1
  fi
}

run_one() {
  disc="$1"; strict="$2"
  dsuf="$(disc_suffix "$disc")"

  # ---- inputs: both channels' fitted yields from the matching out-tree ------
  for c in mu ele; do
    F="$FORK_TEST/pO_fit_out${dsuf}/$c/summary/${c}_fitted_yields.root"
    if [ ! -f "$F" ]; then
      if [ "$strict" = "strict" ]; then
        echo "[ERROR] no fitted yields for disc=$disc: $F"
        echo "        Run the fit first: (fork) ./run_pO_fits.sh both all --disc $disc"
        echo "        (then sync_lxplus.sh download if it ran on lxplus)"
        exit 1
      fi
      echo "[SKIP]  disc=$disc: no fitted yields at $F"
      return 0
    fi
  done

  echo "=================================================================="
  echo "[observables] disc=$disc  (out-tree: pO_fit_out${dsuf})"
  echo "=================================================================="
  mkdir -p "$HERE/../skim/rootfile"

  # ---- 1+2: fitted yields -> disc-tagged charge-asym / FB-ratio graphs -----
  for c in mu ele; do
    F="$FORK_TEST/pO_fit_out${dsuf}/$c/summary/${c}_fitted_yields.root"
    CA="../skim/rootfile/charge_asym_fit_${c}_${disc}.root"
    FB="../skim/rootfile/FBratio_fit_${c}_${disc}.root"
    echo "[observables] $c: charge_asym + FBratio from $F"
    (cd "$HERE" && root -l -b -q "charge_asym.C+(\"$F\",\"$CA\")") || exit 1
    (cd "$HERE" && root -l -b -q "FBratio.C+(\"$F\",\"$FB\")")     || exit 1
    # root returns 0 even when a macro bails out early -> verify the outputs
    require_file "$HERE/$CA" "charge_asym.C failed -- see its [ERROR] above"
    require_file "$HERE/$FB" "FBratio.C failed -- see its [ERROR] above"
  done

  # ---- 3: per-channel observable plots -------------------------------------
  (cd "$HERE/../plotting" && root -l -b -q "observables.C+(false, \"$disc\")") || exit 1
  (cd "$HERE/../plotting" && root -l -b -q "observables.C+(true, \"$disc\")")  || exit 1
  require_file "$HERE/../plotting/plots/charge_asym/$disc/chargeAsym.png"      "observables(false) produced no plot"
  require_file "$HERE/../plotting/plots/Elec/charge_asym/$disc/chargeAsym.png" "observables(true) produced no plot"

  # ---- 4: merged mu+e overlay ----------------------------------------------
  (cd "$HERE/../plotting" && root -l -b -q -e "gROOT->LoadMacro(\"observables.C+\"); observables_overlay(\"$disc\");") || exit 1
  require_file "$HERE/../plotting/plots/merged/$disc/chargeAsym_overlay.png" "observables_overlay produced no plot"

  # ---- 5: fiducial cross sections ------------------------------------------
  (cd "$HERE/../plotting" && root -l -b -q "xsec_fiducial.C+(\"$disc\")") || exit 1
  (cd "$HERE/../plotting" && root -l -b -q -e "gROOT->LoadMacro(\"xsec_fiducial.C+\"); xsec_fiducial_diff(false, \"$disc\"); xsec_fiducial_diff(true, \"$disc\");") || exit 1
  require_file "$HERE/../plotting/plots/xsec/$disc/W_fiducial.png"       "xsec_fiducial produced no plot"
  require_file "$HERE/../plotting/plots/xsec/$disc/W_dsigma_deta_mu.png" "xsec_fiducial_diff produced no plot"

  echo "[observables] disc=$disc DONE. Outputs:"
  echo "    ../skim/rootfile/{charge_asym,FBratio}_fit_{mu,ele}_${disc}.root"
  echo "    ../plotting/plots/{charge_asym,FBratio}/$disc/          (muon)"
  echo "    ../plotting/plots/Elec/{charge_asym,FBratio}/$disc/     (electron)"
  echo "    ../plotting/plots/merged/$disc/                         (mu+e overlay)"
  echo "    ../plotting/plots/xsec/$disc/                           (fiducial sigma)"
}

DISC_ARG="${1:-met}"
case "$DISC_ARG" in
  all)
    for d in met leppt leppt_mt40; do run_one "$d" skip; done
    ;;
  met|leppt|leppt_mt40)
    run_one "$DISC_ARG" strict
    ;;
  *)
    echo "usage: $0 [met|leppt|leppt_mt40|all]"
    exit 1
    ;;
esac
echo "[observables] all requested variants processed."
