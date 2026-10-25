#!/bin/bash
# =============================================================================
# run_observables.sh -- Module 5 driver: fitted yields -> final observable
# plots for ONE W-discriminant variant (or both), carrying the disc tag
# through every filename and output folder so variants never overwrite each
# other (2026-08-03).
#
#   ./run_observables.sh [leppt|leppt_mt40|met|all]        (default: leppt)
#   ./run_observables.sh leppt_mt40_abcdy                   (opt-in, 2026-10-19)
#   ./run_observables.sh met_mt40 | met_mt40_abcdy          (opt-in, 2026-10-20: PF MET
#                                                            on the leppt_mt40 events)
#
# leppt is THE NOMINAL since 2026-10-21 (user decision): the lepton-pT fit with
# NO m_T cut, its SR split at m_T = 30, the scheme-C multijet (fork
# QCD_MODE=mtsplit, out-tree pO_fit_out_leppt/, input
# plotting/combine_input_nomt.C), gen fiducial without m_T. leppt_mt40 (the
# nominal until then) is kept for later, met is the backup. ('leppt' was a
# retired plain no-m_T stack variant 2026-09-21..2026-10-21.)
#
# leppt_mt40_abcdy = the leppt_mt40 fits made with the fork's opt-in
# QCD_MODE=abcdy (rapidity-resolved in-fit ABCD), read from its own out-tree
# pO_fit_out_leppt_mt40_abcdy/; every output carries the full tag
# (plots/{comb,flavfit}/.../leppt_mt40_abcdy/, *_leppt_mt40_abcdy.root,
# logs/run_observables_leppt_mt40_abcdy.log), so the nominal outputs are never
# touched. Not part of 'all'. Comparison with the nominal fit:
# plotting/run_compare_fits.sh (plotting/compare_fits.C).
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
#                                     fit, d(sigma)/d(eta) per charge, e/mu ratios;
#                                     2026-10-22: also from the fb fits, *_fb, and
#                                     each fit's lab-vs-fb overlay xsec_fiducial_labfb)
#     xsec_summary(disc)              xsec/<disc>/xsec_summary_vertical  (2026-10-06:
#                                     W+/W-/W/Z panels, fits as rows + the mu+e band,
#                                     POWHEG x EPPS21 + theory_fid_predictions.txt rows)
#     xsec_contour_WZ_flav(disc, bn)  xsec/<disc>/xsec_contour_WZ_<bn>  (lab AND fb)
#   and each fit's own plots: xsec_contour_WZ_fit -> xsec/<disc>/..._<bn>_<flav>,
#   postfit_incl_fit -> postfit_incl/<disc>/postfit_<flav>_{Wp,Wm,W}.
#   r_per_bin(disc, ref) -> r/<disc>/ (2026-10-20): the POSTFIT r per rapidity
#   bin of mu-only | mu+e | e-only, r_e/r_mu, the e multijet each group would
#   need for r_e = r_mu against its fitted uncertainty, and for a fit-model
#   variant tag (leppt_mt40_abcdy) the base tag's fits overlaid (ref_tag).
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
# Log: analysis/logs/run_observables_<arg>.log = the whole console output
# (2026-09-29; the record of every quoted number, see the bottom of the file).
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
    # the fork's QCD_MODE=abcdy out-tree (2026-10-19, opt-in; see the header)
    leppt_mt40_abcdy) echo "_leppt_mt40_abcdy" ;;
    # the PF-MET counterpart of leppt_mt40 (2026-10-20): the same events, the
    # same in-fit ABCD and m_T > 40 gen fiducial -> pO_fit_out_met_mt40[_abcdy]/
    met_mt40)         echo "_met_mt40" ;;
    met_mt40_abcdy)   echo "_met_mt40_abcdy" ;;
    # the nominal since 2026-10-21: the no-m_T lepton-pT fit (QCD_MODE=mtsplit)
    leppt)      echo "_leppt" ;;
    *)          echo "[ERROR] unknown disc '$1' (leppt|met|leppt_mt40|leppt_mt40_abcdy|met_mt40|met_mt40_abcdy)" >&2; exit 1 ;;
  esac
}

# ref_tag <disc>: the tag whose fits differ from <disc>'s only in the fit model
# (same data, same inputs) -- r_per_bin.C overlays them. Empty for a plain
# discriminant (= pODisc::Base in plotting/disc_variants.h).
ref_tag() {
  case "$1" in
    leppt_mt40_abcdy) echo "leppt_mt40" ;;
    met_mt40_abcdy)   echo "met_mt40" ;;
    *)                echo "" ;;
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
  # $d also picks the gen fiducial set (pODisc::GenFid: m_T > 40 for leppt_mt40)
  (cd "$HERE" && root -l -b -q "fiducial_yields.C+(\"$SUMD/${T}_W_yields.csv\",\"$SUMD/${T}_fitted_yields.root\",\"$FY\",\"$d\")") || exit 1
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
      # a fit-model variant tag (<disc>_abcdy) is the fork's QCD_MODE=abcdy on <disc>
      local fdisc="${disc%_abcdy}" fenv=""
      [ "$fdisc" != "$disc" ] && fenv="QCD_MODE=abcdy "
      echo "        Run the fit first: (fork) ${fenv}./run_pO_fits.sh --disc $fdisc              (simfit, DEFAULT)"
      echo "        and/or:            (fork) ${fenv}./run_pO_fits.sh both flavfit --disc $fdisc (mu-only + e-only)"
      echo "        (both at once: ${fenv}./run_pO_fits.sh both all --disc $fdisc; then sync_lxplus.sh download if on lxplus)"
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
    # the same from the fb fit (2026-10-22: the fb-fit control) + the lab-vs-fb
    # d(sigma)/d(eta) overlay with r_eff of each (lab over the fb window too)
    (cd "$HERE/../plotting" && root -l -b -q -e "gROOT->LoadMacro(\"xsec_fiducial.C+\"); xsec_fiducial_comb(\"$disc\",\"fb\"); xsec_fiducial_labfb(\"$disc\",\"comb\");") || exit 1
    require_file "$HERE/../plotting/plots/comb/xsec/$disc/W_fiducial_fb.png" "xsec_fiducial_comb($disc,fb) produced no plot"
    require_file "$HERE/../plotting/plots/comb/xsec/$disc/W_dsigma_deta_labfb_W.png" "xsec_fiducial_labfb($disc,comb) produced no plot"
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
    # the fb fits (2026-10-22) + each flavour fit's lab-vs-fb overlay
    (cd "$HERE/../plotting" && root -l -b -q -e "gROOT->LoadMacro(\"xsec_fiducial.C+\"); xsec_fiducial_flav(\"$disc\",\"fb\");") || exit 1
    require_file "$HERE/../plotting/plots/flavfit/xsec/$disc/W_fiducial_fb.png" "xsec_fiducial_flav($disc,fb) produced no plot"
    for c in $FLAVS; do
      (cd "$HERE/../plotting" && root -l -b -q -e "gROOT->LoadMacro(\"xsec_fiducial.C+\"); xsec_fiducial_labfb(\"$disc\",\"simfit_$c\");") || exit 1
      require_file "$HERE/../plotting/plots/flavfit/xsec/$disc/W_dsigma_deta_labfb_${c}_W.png" "xsec_fiducial_labfb($disc,simfit_$c) produced no plot"
    done
    # ---- the vertical summary (ttbar-in-PbPb Fig. 9 style): W+/W-/W/Z panels,
    # fits as rows + the mu+e band, theory rows (POWHEG x EPPS21 + any extra
    # predictions in plotting/theory_fid_predictions.txt); reads xsec_flavfit.csv
    (cd "$HERE/../plotting" && root -l -b -q -e "gROOT->LoadMacro(\"xsec_summary.C+\"); xsec_summary(\"$disc\");") || exit 1
    require_file "$HERE/../plotting/plots/flavfit/xsec/$disc/xsec_summary_vertical.png" "xsec_summary produced no plot"
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
    # ---- the POSTFIT r per rapidity bin, mu vs e (2026-10-20): mu-only | mu+e |
    # e-only, r_e/r_mu, the multijet each group would need for r_e = r_mu, and --
    # for a fit-model variant -- its base tag's fits overlaid (the other model).
    # Needs both flavour fits.
    if [ "$FLAVS" = "mu ele" ]; then
      REF="$(ref_tag "$disc")"
      (cd "$HERE/../plotting" && root -l -b -q "r_per_bin.C+(\"$disc\",\"$REF\")") || exit 1
      for bn in lab fb; do
        require_file "$HERE/../plotting/plots/flavfit/r/$disc/r_Wm_$bn.png" "r_per_bin produced no plot"
        require_file "$HERE/../plotting/plots/flavfit/r/$disc/qcd_needed_Wm_$bn.png" "r_per_bin produced no qcd_needed plot"
      done
      if [ -n "$REF" ] && [ -f "$FORK_TEST/pO_fit_out$(disc_suffix "$REF")/simfit_ele/summary/simfit_ele_W_yields.csv" ]; then
        require_file "$HERE/../plotting/plots/flavfit/r/$disc/r_model_Wm_lab.png" "r_per_bin produced no model plot vs $REF"
      fi
    else
      echo "[note] disc=$disc: r_per_bin needs both flavour fits (have: $FLAVS) -> skipped"
    fi
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
    echo "    ../plotting/plots/flavfit/r/$disc/                      (postfit r per bin, mu vs e, multijet needed for r_e = r_mu)"
  fi
}

DISC_ARG="${1:-leppt}"
case "$DISC_ARG" in
  all|leppt|met|leppt_mt40|leppt_mt40_abcdy|met_mt40|met_mt40_abcdy) ;;
  *)
    echo "usage: $0 [leppt|met|leppt_mt40|all|leppt_mt40_abcdy|met_mt40|met_mt40_abcdy]"
    exit 1
    ;;
esac

main() {
  case "$DISC_ARG" in
    all)            for d in leppt met leppt_mt40; do run_one "$d" skip; done ;;
    leppt|met|leppt_mt40|leppt_mt40_abcdy|met_mt40|met_mt40_abcdy) run_one "$DISC_ARG" strict ;;
  esac
  echo "[observables] all requested variants processed."
}

# THE record (2026-09-29): the macros print the numbers the plots are quoted
# with -- sigma per fit with stat/syst, the e/mu ratios, the mu-vs-e chi2, the
# r-covariance source lines -- to stdout only, so the whole run is teed to a
# log named from the argument (one per variant; the next run of the same
# variant overwrites it). run_one's `exit 1` ends the pipeline's subshell, so
# the script's exit status is main's (PIPESTATUS), not tee's.
mkdir -p "$HERE/logs"
LOG="$HERE/logs/run_observables_${DISC_ARG}.log"
main 2>&1 | tee "$LOG"
st=${PIPESTATUS[0]}
echo "[observables] log: $LOG"
exit "$st"
