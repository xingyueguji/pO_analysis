#ifndef PO_DISC_VARIANTS_H
#define PO_DISC_VARIANTS_H

#include "TString.h"
#include "TSystem.h"
#include <iostream>

// =============================================================================
// disc_variants.h -- single source for the W-discriminant variant tags used
// DOWNSTREAM of the fit (2026-08-03).
//
// The fork fits two W discriminants (run_pO_fits.sh --disc met|leppt_mt40)
// into separate out-trees pO_fit_out[_leppt_mt40]/. The result files carry no
// internal marker of which discriminant produced them, so the SAME tag must be
// carried through the whole observables chain or variants silently overwrite /
// mislabel each other. Every consumer (observables.C, xsec_fiducial.C, ...)
// derives its default input paths and its output folder from the tag via this
// header (and WHICH fit -- grand / mu-only / e-only -- via fit_variants.h):
//   fitted yields : <fork>/test/pO_fit_out<suffix>/<fit dir>/summary/...
//   graphs        : ../skim/rootfile/{charge_asym,FBratio}_fid_<fit>_<disc>.root
//                   (FIDUCIAL, r x sigma_gen -- see GraphFile below)
//   plots         : plots/{comb,flavfit}/{charge_asym,FBratio,xsec,postfit_incl}/<disc>/
// One-command runner for the whole chain: analysis/run_observables.sh
// (see README Module 5).
//
// FIT-MODEL VARIANT TAGS (2026-10-19): "leppt_mt40_abcdy" = the leppt_mt40
// discriminant FITTED WITH THE RAPIDITY-RESOLVED IN-FIT ABCD (the fork's
// opt-in QCD_MODE=abcdy, whose out-tree is pO_fit_out_leppt_mt40_abcdy/).
// The full tag goes into every downstream file and folder name, so its
// outputs sit next to the nominal ones and never overwrite them; its
// discriminant -- Base(tag) -- picks what is the same for both: the Combine
// inputs the fit read (InputSuffix) and the gen fiducial set (GenFid). Only
// named explicitly: run_observables.sh 'all' runs met + leppt_mt40.
//
// met_mt40 (2026-10-20, user decision): PF MET as the discriminant on EXACTLY
// the leppt_mt40 events (pT > 25 && m_T > 40) -- the apples-to-apples MET
// counterpart, with the same in-fit ABCD regions and groups and the same
// m_T > 40 gen fiducial (so sigma, A_ch and R_FB compare directly). Fork:
// run_pO_fits.sh --disc met_mt40 -> pO_fit_out_met_mt40[_abcdy]/; tags
// met_mt40 and met_mt40_abcdy here.
// =============================================================================
namespace pODisc
{

// tag -> the W discriminant it was fitted on ("leppt_mt40_abcdy" ->
// "leppt_mt40", "met_mt40_abcdy" -> "met_mt40"; a plain discriminant tag maps
// to itself)
inline TString Base(const char *disc)
{
    const TString d(disc ? disc : "");
    if (d == "leppt_mt40_abcdy") return "leppt_mt40";
    if (d == "met_mt40_abcdy")   return "met_mt40";
    return d;
}

// the fitted variable is PF MET (met, met_mt40[_abcdy]) rather than lepton pT
inline bool IsMetDisc(const char *disc)
{
    const TString d = Base(disc);
    return d == "met" || d == "met_mt40";
}

// disc -> (out-tree suffix, short human label for the plot info box).
// Returns false (with an error) on an unknown tag so a typo cannot silently
// land outputs in a wrong folder.
inline bool Spec(const char *disc, TString &suffix, TString &label)
{
    const TString d(disc ? disc : "");
    if (d == "met")        { suffix = "";            label = "PF MET";               return true; }
    if (d == "leppt_mt40") { suffix = "_leppt_mt40"; label = "lep p_{T} (m_{T}>40)"; return true; }
    if (d == "leppt_mt40_abcdy")
    {
        suffix = "_leppt_mt40_abcdy"; label = "lep p_{T} (m_{T}>40), y-group ABCD"; return true;
    }
    if (d == "met_mt40")   { suffix = "_met_mt40";   label = "PF MET (m_{T}>40)";    return true; }
    if (d == "met_mt40_abcdy")
    {
        suffix = "_met_mt40_abcdy"; label = "PF MET (m_{T}>40), y-group ABCD"; return true;
    }
    // THE NOMINAL since 2026-10-21 (user decision): the lepton-pT fit with NO
    // m_T cut, its SR split at m_T = 30 and the scheme-C multijet (fork
    // QCD_MODE=mtsplit; input plotting/combine_input_nomt.C). Retired as a plain
    // no-m_T stack variant 2026-09-21 -- restored with its own model.
    if (d == "leppt")      { suffix = "_leppt";      label = "lep p_{T} (no m_{T} cut)"; return true; }
    std::cerr << "[ERROR] unknown discriminant tag '" << d
              << "' (use leppt|met|leppt_mt40|leppt_mt40_abcdy|met_mt40|met_mt40_abcdy)\n";
    return false;
}

// tag -> the suffix of the Combine inputs (plots[/Elec]/combine_input_W<suffix>.root)
// the fit read: the out-tree suffix of its discriminant (a fit-model variant
// reads the inputs of its discriminant; the abcdy objects live in that file)
inline TString InputSuffix(const char *disc)
{
    TString suffix, label;
    if (!Spec(Base(disc).Data(), suffix, label)) return "";
    return suffix;
}

// disc -> the GEN FIDUCIAL set of skim/gen_xsec.C it is measured in
// (2026-10-13). Each discriminant's sigma = r x sigma_gen uses the gen twin of
// its own reco selection: met has no m_T cut, leppt_mt40 cuts m_T > 40 (from
// the W's lepton and neutrino). `tag` is inserted after the charge in every
// gen histogram name -- h_gen_sig_Wp<tag>[_FB][_epps21] -- and in the sidecar
// name skim/output/gen_xsec_fid<tag>.txt; `label` is for plot text. The Z
// (h_gen_sig_Z) is the same for both. Keep in step with gen_xsec.C kFidTag.
inline bool GenFid(const char *disc, TString &tag, TString &label)
{
    const TString d = Base(disc); // a fit-model variant shares its discriminant's set
    // (plain p_{T}: the callers prefix "bare lepton", and TLatex stacks
    //  p_{T}^{l} badly)
    // met and (since 2026-10-21) the nominal no-m_T leppt: no m_T cut
    if (d == "met" || d == "leppt") { tag = ""; label = "p_{T} > 25 GeV, |#eta_{lab}| < 2.4";                   return true; }
    if (d == "leppt_mt40" || d == "met_mt40") // the same selection, hence the same fiducial set
                           { tag = "_mt40"; label = "p_{T} > 25 GeV, |#eta_{lab}| < 2.4, m_{T} > 40 GeV";   return true; }
    std::cerr << "[ERROR] GenFid: no gen fiducial set for discriminant '" << d
              << "' (use leppt|met|leppt_mt40|met_mt40)\n";
    return false;
}

// Default path of a charge_asym.C / FBratio.C output (stem = "charge_asym" or
// "FBratio") for one fit (comb | simfit_mu | simfit_ele), relative to
// plotting/: <stem>_fid_<fit>_<disc>.root -- the graphs built from the FIDUCIAL
// yields r x sigma_gen (analysis/fiducial_yields.C). Since 2026-09-22 that is
// the ONLY input any observable reads (user: never raw counts -- there is no
// dedicated efficiency or acceptance correction, r x sigma_gen applies it).
// The count-based <stem>_fit_<chan>_<disc>.root files of earlier versions,
// and their untagged pre-2026-08-03 met fallback, are deliberately NOT read.
inline TString GraphFile(const char *stem, const char *fit, const char *disc)
{
    return TString::Format("../skim/rootfile/%s_fid_%s_%s.root", stem, fit, disc);
}

} // namespace pODisc

#endif
