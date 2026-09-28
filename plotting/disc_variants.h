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
// =============================================================================
namespace pODisc
{

// disc -> (out-tree suffix, short human label for the plot info box).
// Returns false (with an error) on an unknown tag so a typo cannot silently
// land outputs in a wrong folder.
inline bool Spec(const char *disc, TString &suffix, TString &label)
{
    const TString d(disc ? disc : "");
    if (d == "met")        { suffix = "";            label = "PF MET";               return true; }
    if (d == "leppt_mt40") { suffix = "_leppt_mt40"; label = "lep p_{T} (m_{T}>40)"; return true; }
    // The plain no-m_T-cut "leppt" variant was dropped as a discriminant on
    // 2026-08-16 and fully retired on 2026-09-21 (mtandmet.C no longer writes
    // its stacks or combine_input_W_leppt.root). Named explicitly so a stale
    // command line gets an accurate message instead of "unknown tag".
    if (d == "leppt")
    {
        std::cerr << "[ERROR] discriminant 'leppt' (plain W selection) was retired "
                     "2026-09-21 -- use leppt_mt40 (primary) or met (backup)\n";
        return false;
    }
    std::cerr << "[ERROR] unknown discriminant tag '" << d
              << "' (use met|leppt_mt40)\n";
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
