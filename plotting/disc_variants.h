#ifndef PO_DISC_VARIANTS_H
#define PO_DISC_VARIANTS_H

#include "TString.h"
#include "TSystem.h"
#include <iostream>

// =============================================================================
// disc_variants.h -- single source for the W-discriminant variant tags used
// DOWNSTREAM of the fit (2026-08-03).
//
// The fork fits three W discriminants (run_pO_fits.sh --disc met|leppt|
// leppt_mt40) into separate out-trees pO_fit_out[_leppt[_mt40]]/. The result
// files carry no internal marker of which discriminant produced them, so the
// SAME tag must be carried through the whole observables chain or variants
// silently overwrite / mislabel each other. Every consumer (observables.C,
// xsec_fiducial.C) derives its default input paths and its output folder from
// the tag via this header:
//   fitted yields : <fork>/test/pO_fit_out<suffix>/<chan>/summary/...
//   graphs        : ../skim/rootfile/{charge_asym,FBratio}_fit_<chan>_<disc>.root
//   plots         : plots[/Elec]/{charge_asym,FBratio}/<disc>/,
//                   plots/merged/<disc>/, plots/xsec/<disc>/
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
    if (d == "leppt")      { suffix = "_leppt";      label = "lep p_{T}";            return true; }
    if (d == "leppt_mt40") { suffix = "_leppt_mt40"; label = "lep p_{T} (m_{T}>40)"; return true; }
    std::cerr << "[ERROR] unknown discriminant tag '" << d
              << "' (use met|leppt|leppt_mt40)\n";
    return false;
}

// Default path of a charge_asym.C / FBratio.C output (stem = "charge_asym" or
// "FBratio"), relative to plotting/. Tagged name first
// (<stem>_fit_<chan>_<disc>.root); for met ONLY, fall back to the legacy
// untagged <stem>_fit_<chan>.root when the tagged file does not exist (every
// pre-2026-08-03 untagged file came from the MET chain; a variant must never
// fall back, or an untagged file of unknown provenance re-enters the chain).
inline TString GraphFile(const char *stem, const char *chan, const char *disc)
{
    TString tagged = TString::Format("../skim/rootfile/%s_fit_%s_%s.root", stem, chan, disc);
    if (TString(disc) == "met" && gSystem->AccessPathName(tagged)) // true = missing
    {
        TString legacy = TString::Format("../skim/rootfile/%s_fit_%s.root", stem, chan);
        if (!gSystem->AccessPathName(legacy))
        {
            std::cout << "[INFO] " << tagged << " not found -> using legacy untagged "
                      << legacy << " (met chain only)\n";
            return legacy;
        }
    }
    return tagged;
}

} // namespace pODisc

#endif
