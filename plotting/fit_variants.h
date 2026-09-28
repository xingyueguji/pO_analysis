#ifndef PO_FIT_VARIANTS_H
#define PO_FIT_VARIANTS_H

#include "Rtypes.h" // EColor: kBlack, kBlue, kRed
#include "TString.h"
#include "TSystem.h"
#include <iostream>
#include <vector>

// =============================================================================
// fit_variants.h -- single source for WHICH simultaneous fit a downstream macro
// reads (2026-09-22), the companion of disc_variants.h (which picks the W
// discriminant). The fork's run_pO_fits.sh writes one work dir per fit into the
// discriminant's out-tree pO_fit_out<suffix>/:
//
//   tag         work dir     summary/ files                  the fit
//   comb        simfit/      comb_{W_yields.csv,             the GRAND fit (mode simfit):
//                            summary.csv,fitted_yields.root}  mu + e in one likelihood, r shared
//   simfit_mu   simfit_mu/   simfit_mu_{...}                 the MUON-ONLY fit (mode flavfit)
//   simfit_ele  simfit_ele/  simfit_ele_{...}                the ELECTRON-ONLY fit (mode flavfit)
//
// All three share the layout, the 25 POIs and the treatment (stat/syst split,
// 25x25 POI covariance, contour scan), so a consumer takes the tag and derives
// every path from here -- together with the marker/colour the fit wears in the
// mu-vs-e overlays (plots/flavfit/...). Those avoid every glyph + colour pair
// the CHARGE convention uses elsewhere (W black circle, W+ azure square, W- red
// diamond -- xsec_fiducial.C kMk*) and keep the legacy mu/e glyphs (mu circle,
// e square); sizes are the ones measured to rasterize centred on the error bar
// (circle any size, square 1.5, diamond 2.0 -- see the kMk* note).
//
// The fork location honours $FORK_TEST, which analysis/run_observables.sh
// exports (absolute), so a macro and its driver can never read different trees.
// =============================================================================
namespace pOFit
{

struct Spec
{
    TString tag;                   // comb | simfit_mu | simfit_ele (= the summary-file prefix)
    TString dir;                   // work dir inside pO_fit_out<suffix>/
    TString label;                 // legend / header text
    TString flav;                  // "mu" | "ele" for a per-flavour fit, "" for comb
    std::vector<TString> flavours; // the lepton flavours IN the fit
    int color = kBlack;
    int marker = 20;
    double markerSize = 1.3;
};

inline bool Get(const char *fit, Spec &s)
{
    const TString f(fit ? fit : "");
    s = Spec();
    s.tag = f;
    if (f == "comb")
    {
        s.dir = "simfit"; s.label = "#mu + e combined"; s.flavours = {"mu", "ele"};
        s.color = kBlack; s.marker = 33; s.markerSize = 2.0;
        return true;
    }
    if (f == "simfit_mu")
    {
        s.dir = "simfit_mu"; s.label = "#mu only"; s.flav = "mu"; s.flavours = {"mu"};
        s.color = kBlue + 1; s.marker = 20; s.markerSize = 1.3;
        return true;
    }
    if (f == "simfit_ele")
    {
        s.dir = "simfit_ele"; s.label = "e only"; s.flav = "ele"; s.flavours = {"ele"};
        s.color = kRed + 1; s.marker = 21; s.markerSize = 1.5;
        return true;
    }
    std::cerr << "[ERROR] unknown fit tag '" << f << "' (comb | simfit_mu | simfit_ele)\n";
    return false;
}

// <fork>/test -- $FORK_TEST when set, else the sibling checkout seen from plotting/
inline TString ForkTest()
{
    const char *e = gSystem->Getenv("FORK_TEST");
    return (e && *e) ? TString(e) : TString("../../HiggsAnalysis-CombinedLimit/test");
}

// <fork>/test/pO_fit_out<dsuf>/<fit work dir>
inline TString WorkDir(const TString &dsuf, const Spec &s)
{
    return ForkTest() + "/pO_fit_out" + dsuf + "/" + s.dir;
}

// the fit's summary files, e.g. SummaryFile(dsuf, s, "W_yields.csv")
inline TString SummaryFile(const TString &dsuf, const Spec &s, const char *what)
{
    return WorkDir(dsuf, s) + "/summary/" + s.tag + "_" + what;
}

// true = the fit's fitted-yields file exists for this discriminant
inline bool Available(const TString &dsuf, const Spec &s)
{
    return !gSystem->AccessPathName(SummaryFile(dsuf, s, "fitted_yields.root")); // false = exists
}

} // namespace pOFit

#endif
