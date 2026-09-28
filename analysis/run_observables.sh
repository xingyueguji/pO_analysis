#!/bin/bash
# =============================================================================
# run_observables.sh -- Module 5 driver: fitted yields -> final observable
# plots for ONE W-discriminant variant (or both), carrying the disc tag
# through every filename and output folder so variants never overwrite each
# other (2026-08-03).
#
#   ./run_observables.sh [met|leppt_mt40|all]              (default: met)
#
# leppt_mt40 is the PRIMARY discriminant, met the backup. The third variant
# "leppt" (plain W selection, no m_T cut) was dropped 2026-08-16 and retired
# 2026-09-21 -- mtandmet.C no longer writes its input or its stacks.
#
# Per discriminant it runs (each block only when its fit outputs exist):
#
# PRIMARY -- the simfit GRAND SIMULTANEOUS FIT (2026-08-04, "comb" channel):
#   fiducial_yields.C on the fit (r x sigma_gen per bin + the r covariance;
#   mu+e-shared r's) -> ../skim/rootfile/fidyields_comb_<disc>.root,
#   charge_asym.C + FBratio.C on it -> {charge_asym,FBratio}_fid_comb_<disc>.root
#   (FIDUCIAL since 2026-09-22 -- user: every observable from r x sigma_gen,
#   never raw counts), then observables_comb(disc) -> plots/comb/{charge_asym,FBratio}/<disc>/,
#   and xsec_fiducial_comb(disc) -> plots/comb/xsec/<disc>/ (sigma = r x
#   sigma_gen, cov-propagated, stat bars + syst boxes), the (sigma_W, sigma_Z)
#   contour and the y-inclusive postfit stacks.
#
# PER-FLAVOUR -- the mu-only / e-only simultaneous fits (fork mode flavfit,
# 2026-09-22; the same model and treatment as the grand fit, one flavour each):
#   fiducial_yields.C: r_i x sigma_gen-fid,i + the r covariance of each fit
#   (and of the grand one, as the reference) -> ../skim/rootfile/fidyields_<tag>_<disc>.root,
#   charge_asym.C + FBratio.C on those -> {charge_asym,FBratio}_fid_<tag>_<disc>.root
#   (ACCEPTANCE-CORRECTED: the count-based R_FB of mu and e is not comparable,
#   see fiducial_yields.C), then the mu-vs-e OVERLAYS into plots/flavfit/:
#     observables_flav(disc)          {charge_asym,FBratio}/<disc>/   (+ mu-vs-e chi2)
#     xsec_fiducial_flav(disc)        xsec/<disc>/  (sigma W+/W-/W incl. the grand
#                                     fit, d(sigma)/d(eta) per charge, e/mu ratios)
#     xsec_contour_WZ_flav(disc, bn)  xsec/<disc>/xsec_contour_WZ_<bn>  (lab AND fb)
#   and each fit's own plots: xsec_contour_WZ_fit -> xsec/<disc>/..._<bn>_<flav>,
#   postfit_incl_fit -> postfit_incl/<disc>/postfit_<flav>_{Wp,Wm,W}.
# (The LEGACY per-flavour per-bin chain -- observables(isElec), observables_overlay,
#  xsec_fiducial(_diff) on <chan>_fitted_yields.root -- was removed with those
#  fits on 2026-09-22.)
#
# Input: the fork's out-tree for that variant (fit already run + downloaded):
#   simfit  : $FORK_TEST/pO_fit_out<suffix>/simfit/summary/comb_fitted_yields.root
#             produced by run_pO_fits.sh [--disc <disc>]   (simfit is the default)
#   flavfit : $FORK_TEST/pO_fit_out<suffix>/simfit_{mu,ele}/summary/simfit_{mu,ele}_fitted_yields.root
#             produced by run_pO_fits.sh both flavfit [--disc <disc>]  (or mode all)
# With 'all', variants with no outputs at all are SKIPPED with a warning; a
# single named variant errors out instead.
#
# FORK_TEST env var overrides the fork location (default: sibling checkout);
# it is exported ABSOLUTE, so the ROOT macros (plotting/fit_variants.h) read
# exactly the tree this script checked.
# bash-3.2 compatible (macOS stock bash): no associative arrays.
# =============================================================================
set -u

HERE="$(cd "$(dirname "$0")" && pwd)"
FORK_TEST="${FORK_TEST:-$HERE/../../HiggsAnalysis-CombinedLimit/test}"
# absolute + exported: the macros run from plotting/ and analysis/ (different
# CWDs), and plotting/fit_variants.h reads $FORK_TEST -- one tree for everyone
if [ -d "$FORK_TEST" ]; then FORK_TEST="$(cd "$FORK_TEST" && pwd)"; fi
export FORK_TEST

disc_suffix() {
  case "$1" in
    met)        echo "" ;;
    leppt_mt40) echo "_leppt_mt40" ;;
    leppt)      echo "[ERROR] disc 'leppt' (plain W selection) was retired 2026-09-21" \
                     "-- use leppt_mt40 (primary) or met (backup)" >&2; exit 1 ;;
    *)          echo "[ERROR] unknown disc '$1' (met|leppt_mt40)" >&2; exit 1 ;;
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

# fid_observables <fit tag> <fork work dir> <disc>: the FIDUCIAL yields
# r x sigma_gen of one fit (fiducial_yields.C, with the fit's r covariance)
# and the charge asymmetry + F/B ratio built from them -- the ONLY inputs of
# every observable since 2026-09-22 (user: "all the observables should be
# using r * gen xsec, never raw counts" -- there is no dedicated efficiency or
# acceptance correction; r x sigma_gen applies it from MC). The raw fitted
# yields r x S stay the fit's record (and feed the postfit stacks).
fid_observables() {
  T="$1"; SUMD="$2/summary"; d="$3"
  FY="../skim/rootfile/fidyields_${T}_${d}.root"
  CA="../skim/rootfile/charge_asym_fid_${T}_${d}.root"
  FB="../skim/rootfile/FBratio_fid_${T}_${d}.root"
  echo "[observables] $T: fiducial yields (r x sigma_gen) -> charge_asym + FBratio"
  (cd "$HERE" && root -l -b -q "fiducial_yields.C+(\"$SUMD/${T}_W_yields.csv\",\"$SUMD/${T}_fitted_yields.root\",\"$FY\")") || exit 1
  require_file "$HERE/$FY" "fiducial_yields.C failed -- see its [ERROR] above"
  (cd "$HERE" && root -l -b -q "charge_asym.C+(\"$FY\",\"$CA\")") || exit 1
  (cd "$HERE" && root -l -b -q "FBratio.C+(\"$FY\",\"$FB\")")     || exit 1
  # root returns 0 even when a macro bails out early -> verify the outputs
  require_file "$HERE/$CA" "charge_asym.C failed -- see its [ERROR] above"
  require_file "$HERE/$FB" "FBratio.C failed -- see its [ERROR] above"
}

run_one() {
  disc="$1"; strict="$2"
  # NB disc_suffix's `exit 1` fires inside this command substitution's subshell
  # and would NOT stop the script -- an unknown tag would silently get met's
  # empty suffix and land its plots in a mislabeled folder. Check the status.
  dsuf="$(disc_suffix "$disc")" || exit 1

  # ---- what exists for this variant? ----------------------------------------
  # PRIMARY (2026-08-04): the simfit grand fit -> comb_fitted_yields.root.
  # PER-FLAVOUR (2026-09-22): the flavfit trees simfit_{mu,ele}/ -> the mu-vs-e
  # overlays (either flavour alone is drawn too, with a WARN for the other).
  COMBF="$FORK_TEST/pO_fit_out${dsuf}/simfit/summary/comb_fitted_yields.root"
  HAVE_COMB=0; [ -f "$COMBF" ] && HAVE_COMB=1
  FLAVS=""
  for c in mu ele; do
    [ -f "$FORK_TEST/pO_fit_out${dsuf}/simfit_$c/summary/simfit_${c}_fitted_yields.root" ] && FLAVS="$FLAVS $c"
  done
  FLAVS="${FLAVS# }"
  HAVE_FLAV=0; [ -n "$FLAVS" ] && HAVE_FLAV=1

  if [ "$HAVE_COMB" -eq 0 ] && [ "$HAVE_FLAV" -eq 0 ]; then
    if [ "$strict" = "strict" ]; then
      echo "[ERROR] no fit outputs for disc=$disc:"
      echo "        simfit  : $COMBF"
      echo "        flavfit : pO_fit_out${dsuf}/simfit_{mu,ele}/summary/simfit_{mu,ele}_fitted_yields.root"
      echo "        Run the fit first: (fork) ./run_pO_fits.sh --disc $disc              (simfit, DEFAULT)"
      echo "        and/or:            (fork) ./run_pO_fits.sh both flavfit --disc $disc (mu-only + e-only)"
      echo "        (both at once: ./run_pO_fits.sh both all --disc $disc; then sync_lxplus.sh download if on lxplus)"
      exit 1
    fi
    echo "[SKIP]  disc=$disc: no fit outputs (neither simfit nor flavfit)"
    return 0
  fi

  echo "=================================================================="
  echo "[observables] disc=$disc  (out-tree: pO_fit_out${dsuf};  simfit: $([ "$HAVE_COMB" -eq 1 ] && echo yes || echo NO), flavfit: ${FLAVS:-NO})"
  echo "=================================================================="
  mkdir -p "$HERE/../skim/rootfile"

  # ==== PRIMARY: simfit (grand simultaneous fit) -> comb observables =========
  if [ "$HAVE_COMB" -eq 1 ]; then
    # FIDUCIAL (r x sigma_gen) charge asymmetry + F/B, never the raw counts
    fid_observables comb "$FORK_TEST/pO_fit_out${dsuf}/simfit" "$disc"
    (cd "$HERE/../plotting" && root -l -b -q -e "gROOT->LoadMacro(\"observables.C+\"); observables_comb(\"$disc\");") || exit 1
    require_file "$HERE/../plotting/plots/comb/charge_asym/$disc/chargeAsym.png" "observables_comb produced no plot"
    # comb fiducial sigma = r x sigma_gen-fid (per flavour, covariance-propagated)
    # + the per-flavour extraction diagnostic (gen / reco / r x gen / r x reco)
    (cd "$HERE/../plotting" && root -l -b -q -e "gROOT->LoadMacro(\"xsec_fiducial.C+\"); xsec_fiducial_comb(\"$disc\"); xsec_fiducial_diag(\"$disc\");") || exit 1
    require_file "$HERE/../plotting/plots/comb/xsec/$disc/W_fiducial.png" "xsec_fiducial_comb produced no plot"
    require_file "$HERE/../plotting/plots/comb/xsec/$disc/W_xsec_diag_mu.png" "xsec_fiducial_diag produced no plot"
    # (sigma_W, sigma_Z) confidence ellipse from the 25x25 POI covariance -- the
    # plane for comparing nPDF sets. Needs h_cov_poi (fork extractor, 2026-09-15)
    # and h_gen_sig_Z (skim/gen_xsec.C); the macro says so and bails if absent,
    # which is why this is not gated by require_file on a pre-09-15 tree.
    # BOTH binnings (2026-09-21): the fb variant sums a DIFFERENT lab window, so
    # it is a genuinely different number (sigma_W 93.22 vs 108.23 lab), not a
    # cosmetic twin. Until now only the default "lab" was run here, so
    # xsec_contour_WZ_fb.{png,csv} silently kept whatever the LAST manual run
    # left -- found stale by a full fit cycle on 2026-09-21 (it still showed the
    # previous fit's 92.01 nb). Both are cheap; run both.
    for bn in lab fb; do
      (cd "$HERE/../plotting" && root -l -b -q -e "gROOT->LoadMacro(\"xsec_contour.C+\"); xsec_contour_WZ(\"$disc\",\"$bn\");") || exit 1
      require_file "$HERE/../plotting/plots/comb/xsec/$disc/xsec_contour_WZ_$bn.png" \
                   "xsec_contour_WZ($disc,$bn) produced no plot"
    done
    # rapidity-INCLUSIVE postfit stacks (per flavour, Wp/Wm/W; exact because the
    # simfit parameters are pure normalizations -- see postfit_incl.C header)
    (cd "$HERE/../plotting" && root -l -b -q "postfit_incl.C+(\"$disc\")") || exit 1
    require_file "$HERE/../plotting/plots/comb/postfit_incl/$disc/postfit_mu_W.png" "postfit_incl produced no plot"
  else
    echo "[note] disc=$disc: no simfit output -> comb (primary) chain skipped"
  fi

  # ==== PER-FLAVOUR: the mu-only / e-only simultaneous fits (flavfit) =======
  if [ "$HAVE_FLAV" -eq 1 ]; then
    # ---- FIDUCIAL yields r x sigma_gen (+ r covariance) -> charge-asym / FB-ratio
    # graphs, per flavour fit (the grand fit's were made in the comb block and
    # serve as the overlays' reference). The count-based r x S yields would not
    # even be comparable here: their A x eps does not cancel in R_FB (F and B
    # are different |eta_lab| regions), so mu and e differ by up to 60% on
    # IDENTICAL physics -- see analysis/fiducial_yields.C.
    FIDTAGS=""
    for c in $FLAVS; do
      FIDTAGS="$FIDTAGS simfit_$c"
      fid_observables "simfit_$c" "$FORK_TEST/pO_fit_out${dsuf}/simfit_$c" "$disc"
    done
    [ "$HAVE_COMB" -eq 1 ] && FIDTAGS="$FIDTAGS comb"
    # ---- the mu-vs-e overlays: charge asymmetry + R_FB, cross sections --------
    (cd "$HERE/../plotting" && root -l -b -q -e "gROOT->LoadMacro(\"observables.C+\"); observables_flav(\"$disc\");") || exit 1
    require_file "$HERE/../plotting/plots/flavfit/charge_asym/$disc/chargeAsym.png" "observables_flav produced no plot"
    (cd "$HERE/../plotting" && root -l -b -q -e "gROOT->LoadMacro(\"xsec_fiducial.C+\"); xsec_fiducial_flav(\"$disc\");") || exit 1
    require_file "$HERE/../plotting/plots/flavfit/xsec/$disc/W_fiducial.png" "xsec_fiducial_flav produced no plot"
    # ---- the (sigma_W, sigma_Z) plane, BOTH binnings (see the comb block):
    # the overlay + each fit's own plot (Gaussian vs profiled, stat-only)
    for bn in lab fb; do
      (cd "$HERE/../plotting" && root -l -b -q -e "gROOT->LoadMacro(\"xsec_contour.C+\"); xsec_contour_WZ_flav(\"$disc\",\"$bn\");") || exit 1
      require_file "$HERE/../plotting/plots/flavfit/xsec/$disc/xsec_contour_WZ_$bn.png" \
                   "xsec_contour_WZ_flav($disc,$bn) produced no plot"
      for c in $FLAVS; do
        (cd "$HERE/../plotting" && root -l -b -q -e "gROOT->LoadMacro(\"xsec_contour.C+\"); xsec_contour_WZ_fit(\"$disc\",\"$bn\",\"simfit_$c\");") || exit 1
        require_file "$HERE/../plotting/plots/flavfit/xsec/$disc/xsec_contour_WZ_${bn}_$c.png" \
                     "xsec_contour_WZ_fit($disc,$bn,simfit_$c) produced no plot"
      done
    done
    # ---- each fit's rapidity-inclusive postfit stacks (its own flavour only) ---
    for c in $FLAVS; do
      (cd "$HERE/../plotting" && root -l -b -q -e "gROOT->LoadMacro(\"postfit_incl.C+\"); postfit_incl_fit(\"$disc\",\"simfit_$c\");") || exit 1
      require_file "$HERE/../plotting/plots/flavfit/postfit_incl/$disc/postfit_${c}_W.png" "postfit_incl_fit(simfit_$c) produced no plot"
    done
  else
    echo "[note] disc=$disc: no flavfit outputs (simfit_{mu,ele}/) -> mu-vs-e overlays skipped"
  fi

  echo "[observables] disc=$disc DONE. Outputs:"
  if [ "$HAVE_COMB" -eq 1 ]; then
    echo "    ../skim/rootfile/{fidyields,charge_asym_fid,FBratio_fid}_comb_${disc}.root"
    echo "    ../plotting/plots/comb/{charge_asym,FBratio}/$disc/     (PRIMARY: simfit mu+e)"
    echo "    ../plotting/plots/comb/xsec/$disc/                      (PRIMARY: simfit fiducial sigma vs MC)"
    echo "    ../plotting/plots/comb/postfit_incl/$disc/              (PRIMARY: y-inclusive postfit stacks)"
  fi
  if [ "$HAVE_FLAV" -eq 1 ]; then
    echo "    ../skim/rootfile/{fidyields,charge_asym_fid,FBratio_fid}_{$(echo $FIDTAGS | tr ' ' ',')}_${disc}.root"
    echo "    ../plotting/plots/flavfit/{charge_asym,FBratio}/$disc/  (mu-only vs e-only overlays + chi2 CSV)"
    echo "    ../plotting/plots/flavfit/xsec/$disc/                   (sigma + d(sigma)/d(eta) + contour overlays, per-fit contours)"
    echo "    ../plotting/plots/flavfit/postfit_incl/$disc/           (each flavour fit's y-inclusive postfit stacks)"
  fi
}

DISC_ARG="${1:-met}"
case "$DISC_ARG" in
  all)
    for d in met leppt_mt40; do run_one "$d" skip; done
    ;;
  met|leppt_mt40)
    run_one "$DISC_ARG" strict
    ;;
  leppt)
    echo "[ERROR] disc 'leppt' (plain W selection, no m_T cut) was dropped as a"
    echo "        discriminant 2026-08-16 and retired 2026-09-21: mtandmet.C no"
    echo "        longer writes combine_input_W_leppt.root or plots[/Elec]/leppt/."
    echo "        Use leppt_mt40 (primary) or met (backup)."
    exit 1
    ;;
  *)
    echo "usage: $0 [met|leppt_mt40|all]"
    exit 1
    ;;
esac
echo "[observables] all requested variants processed."
