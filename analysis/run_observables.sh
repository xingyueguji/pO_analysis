#!/bin/bash
# =============================================================================
# run_observables.sh -- Module 5 driver: fitted yields -> final observable
# plots for ONE W-discriminant variant (or all three), carrying the disc tag
# through every filename and output folder so variants never overwrite each
# other (2026-08-03).
#
#   ./run_observables.sh [met|leppt|leppt_mt40|all]        (default: met)
#
# Per discriminant it runs (each block only when its fit outputs exist):
#
# PRIMARY -- the simfit GRAND SIMULTANEOUS FIT (2026-08-04, "comb" channel):
#   charge_asym.C + FBratio.C on pO_fit_out<suffix>/simfit/summary/
#   comb_fitted_yields.root (mu+e-shared r's; errors use the h_cov_yield[_FB]
#   covariance) -> ../skim/rootfile/{charge_asym,FBratio}_fit_comb_<disc>.root,
#   then observables_comb(disc) -> plots/comb/{charge_asym,FBratio}/<disc>/,
#   and xsec_fiducial_comb(disc) -> plots/comb/xsec/<disc>/ (mu+e-summed
#   sigma_fid, cov-propagated, post-fit vs prefit-MC open-marker overlay).
#
# LEGACY -- the per-flavour per-bin fits (kept for comparison), for mu and ele:
#   1. charge_asym.C   fitted yields -> ../skim/rootfile/charge_asym_fit_<chan>_<disc>.root
#   2. FBratio.C       fitted yields -> ../skim/rootfile/FBratio_fit_<chan>_<disc>.root
#   3. observables(false|true, disc)  -> plots[/Elec]/{charge_asym,FBratio}/<disc>/
#   4. observables_overlay(disc)      -> plots/merged/<disc>/
#   5. xsec_fiducial(disc) + xsec_fiducial_diff(mu|ele, disc)
#                                     -> plots/xsec/<disc>/
#
# Input: the fork's out-tree for that variant (fit already run + downloaded):
#   simfit  : $FORK_TEST/pO_fit_out<suffix>/simfit/summary/comb_fitted_yields.root
#             produced by run_pO_fits.sh [--disc <disc>]   (simfit is the default)
#   legacy  : $FORK_TEST/pO_fit_out<suffix>/<chan>/summary/<chan>_fitted_yields.root
#             produced by run_pO_fits.sh both all [--disc <disc>]
# With 'all', variants with no outputs at all are SKIPPED with a warning; a
# single named variant errors out instead.
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

  # ---- what exists for this variant? ----------------------------------------
  # PRIMARY (2026-08-04): the simfit grand fit -> comb_fitted_yields.root.
  # LEGACY: the per-flavour per-bin fits ("wrong way": each per-bin card refits
  # the same Z data) -- kept runnable for comparison, used when present.
  COMBF="$FORK_TEST/pO_fit_out${dsuf}/simfit/summary/comb_fitted_yields.root"
  HAVE_COMB=0; [ -f "$COMBF" ] && HAVE_COMB=1
  HAVE_LEGACY=1
  for c in mu ele; do
    [ -f "$FORK_TEST/pO_fit_out${dsuf}/$c/summary/${c}_fitted_yields.root" ] || HAVE_LEGACY=0
  done

  if [ "$HAVE_COMB" -eq 0 ] && [ "$HAVE_LEGACY" -eq 0 ]; then
    if [ "$strict" = "strict" ]; then
      echo "[ERROR] no fit outputs for disc=$disc:"
      echo "        simfit : $COMBF"
      echo "        legacy : pO_fit_out${dsuf}/{mu,ele}/summary/*_fitted_yields.root"
      echo "        Run the fit first: (fork) ./run_pO_fits.sh --disc $disc          (simfit, DEFAULT)"
      echo "        or:                (fork) ./run_pO_fits.sh both all --disc $disc (legacy + simfit)"
      echo "        (then sync_lxplus.sh download if it ran on lxplus)"
      exit 1
    fi
    echo "[SKIP]  disc=$disc: no fit outputs (neither simfit nor legacy per-flavour)"
    return 0
  fi

  echo "=================================================================="
  echo "[observables] disc=$disc  (out-tree: pO_fit_out${dsuf};  simfit: $([ "$HAVE_COMB" -eq 1 ] && echo yes || echo NO), legacy: $([ "$HAVE_LEGACY" -eq 1 ] && echo yes || echo NO))"
  echo "=================================================================="
  mkdir -p "$HERE/../skim/rootfile"

  # ==== PRIMARY: simfit (grand simultaneous fit) -> comb observables =========
  if [ "$HAVE_COMB" -eq 1 ]; then
    CA="../skim/rootfile/charge_asym_fit_comb_${disc}.root"
    FB="../skim/rootfile/FBratio_fit_comb_${disc}.root"
    echo "[observables] comb (simfit): charge_asym + FBratio from $COMBF"
    (cd "$HERE" && root -l -b -q "charge_asym.C+(\"$COMBF\",\"$CA\")") || exit 1
    (cd "$HERE" && root -l -b -q "FBratio.C+(\"$COMBF\",\"$FB\")")     || exit 1
    # root returns 0 even when a macro bails out early -> verify the outputs
    require_file "$HERE/$CA" "charge_asym.C failed -- see its [ERROR] above"
    require_file "$HERE/$FB" "FBratio.C failed -- see its [ERROR] above"
    (cd "$HERE/../plotting" && root -l -b -q -e "gROOT->LoadMacro(\"observables.C+\"); observables_comb(\"$disc\");") || exit 1
    require_file "$HERE/../plotting/plots/comb/charge_asym/$disc/chargeAsym.png" "observables_comb produced no plot"
    # comb fiducial sigma = r x sigma_gen-fid (per flavour, covariance-propagated)
    # + the per-flavour extraction diagnostic (gen / reco / r x gen / r x reco)
    (cd "$HERE/../plotting" && root -l -b -q -e "gROOT->LoadMacro(\"xsec_fiducial.C+\"); xsec_fiducial_comb(\"$disc\"); xsec_fiducial_diag(\"$disc\");") || exit 1
    require_file "$HERE/../plotting/plots/comb/xsec/$disc/W_fiducial.png" "xsec_fiducial_comb produced no plot"
    require_file "$HERE/../plotting/plots/comb/xsec/$disc/W_xsec_diag_mu.png" "xsec_fiducial_diag produced no plot"
    # rapidity-INCLUSIVE postfit stacks (per flavour, Wp/Wm/W; exact because the
    # simfit parameters are pure normalizations -- see postfit_incl.C header)
    (cd "$HERE/../plotting" && root -l -b -q "postfit_incl.C+(\"$disc\")") || exit 1
    require_file "$HERE/../plotting/plots/comb/postfit_incl/$disc/postfit_mu_W.png" "postfit_incl produced no plot"
  else
    echo "[note] disc=$disc: no simfit output -> comb (primary) chain skipped"
  fi

  # ==== LEGACY per-flavour chain (per-bin fits; comparison only) =============
  if [ "$HAVE_LEGACY" -eq 1 ]; then
    # ---- 1+2: fitted yields -> disc-tagged charge-asym / FB-ratio graphs ----
    for c in mu ele; do
      F="$FORK_TEST/pO_fit_out${dsuf}/$c/summary/${c}_fitted_yields.root"
      CA="../skim/rootfile/charge_asym_fit_${c}_${disc}.root"
      FB="../skim/rootfile/FBratio_fit_${c}_${disc}.root"
      echo "[observables] $c (legacy): charge_asym + FBratio from $F"
      (cd "$HERE" && root -l -b -q "charge_asym.C+(\"$F\",\"$CA\")") || exit 1
      (cd "$HERE" && root -l -b -q "FBratio.C+(\"$F\",\"$FB\")")     || exit 1
      require_file "$HERE/$CA" "charge_asym.C failed -- see its [ERROR] above"
      require_file "$HERE/$FB" "FBratio.C failed -- see its [ERROR] above"
    done

    # ---- 3: per-channel observable plots ------------------------------------
    (cd "$HERE/../plotting" && root -l -b -q "observables.C+(false, \"$disc\")") || exit 1
    (cd "$HERE/../plotting" && root -l -b -q "observables.C+(true, \"$disc\")")  || exit 1
    require_file "$HERE/../plotting/plots/charge_asym/$disc/chargeAsym.png"      "observables(false) produced no plot"
    require_file "$HERE/../plotting/plots/Elec/charge_asym/$disc/chargeAsym.png" "observables(true) produced no plot"

    # ---- 4: merged mu+e overlay ---------------------------------------------
    (cd "$HERE/../plotting" && root -l -b -q -e "gROOT->LoadMacro(\"observables.C+\"); observables_overlay(\"$disc\");") || exit 1
    require_file "$HERE/../plotting/plots/merged/$disc/chargeAsym_overlay.png" "observables_overlay produced no plot"

    # ---- 5: fiducial cross sections (needs the legacy incl fits' CSVs) ------
    (cd "$HERE/../plotting" && root -l -b -q "xsec_fiducial.C+(\"$disc\")") || exit 1
    (cd "$HERE/../plotting" && root -l -b -q -e "gROOT->LoadMacro(\"xsec_fiducial.C+\"); xsec_fiducial_diff(false, \"$disc\"); xsec_fiducial_diff(true, \"$disc\");") || exit 1
    require_file "$HERE/../plotting/plots/xsec/$disc/W_fiducial.png"       "xsec_fiducial produced no plot"
    require_file "$HERE/../plotting/plots/xsec/$disc/W_dsigma_deta_mu.png" "xsec_fiducial_diff produced no plot"
  else
    echo "[note] disc=$disc: no legacy per-flavour outputs -> mu/ele plots, overlay + xsec skipped"
  fi

  echo "[observables] disc=$disc DONE. Outputs:"
  if [ "$HAVE_COMB" -eq 1 ]; then
    echo "    ../skim/rootfile/{charge_asym,FBratio}_fit_comb_${disc}.root"
    echo "    ../plotting/plots/comb/{charge_asym,FBratio}/$disc/     (PRIMARY: simfit mu+e)"
    echo "    ../plotting/plots/comb/xsec/$disc/                      (PRIMARY: simfit fiducial sigma vs MC)"
    echo "    ../plotting/plots/comb/postfit_incl/$disc/              (PRIMARY: y-inclusive postfit stacks)"
  fi
  if [ "$HAVE_LEGACY" -eq 1 ]; then
    echo "    ../skim/rootfile/{charge_asym,FBratio}_fit_{mu,ele}_${disc}.root"
    echo "    ../plotting/plots/{charge_asym,FBratio}/$disc/          (legacy muon)"
    echo "    ../plotting/plots/Elec/{charge_asym,FBratio}/$disc/     (legacy electron)"
    echo "    ../plotting/plots/merged/$disc/                         (legacy mu+e overlay)"
    echo "    ../plotting/plots/xsec/$disc/                           (legacy fiducial sigma)"
  fi
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
