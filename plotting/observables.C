#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <functional>

#include "plotting_helper.C"               // PlotStyle, SaveNiceGraph[_ErrorBand][_TwoData]
#include "../analysis/analysis_helpers.h"  // pOAnalysis::YieldInRange (Sumw2-aware yields)
#include "disc_variants.h"                 // pODisc::Spec/GraphFile (W-discriminant tags)

// =============================================================================
// observables.C -- final charge-asymmetry + forward/backward plots from the
// FITTED signal yields (Combine).
//
//   observables_comb(disc)           PRIMARY: the simfit grand fit's mu+e-combined
//                                    observables (plots/comb/..., 2026-08-04)
//   observables(isElec, disc)        per-channel plots from the LEGACY per-bin
//                                    fits (data + all-4-model theory bands)
//   observables_overlay(disc)        MERGED muon+electron overlay (legacy fits)
//
// `disc` = met|leppt|leppt_mt40 selects WHICH fit's yields are plotted
// (2026-08-03): it drives every default input path (the fork out-tree
// pO_fit_out<suffix>/, the tagged charge_asym/FBratio_fit_<chan>_<disc>.root)
// AND the output folders (plots[/Elec]/{charge_asym,FBratio}/<disc>/,
// plots/merged/<disc>/), so the three variants coexist without overwriting.
// The discriminant is also stamped into the plot info box. One-command runner
// for the whole chain: analysis/run_observables.sh (README Module 5).
//
// All theory bands now show ALL FOUR nPDF sets (EPPS21, nCTEQ15HQ, nNNPDF3.0,
// TUJU21nlo) and the obsolete "Projection with Electrons" pseudo-band is gone --
// the real electron measurement is overlaid instead (observables_overlay).
// =============================================================================

// -----------------------------------------------------------------------------
// Build the 4 count-weighted SUM-channel theory graphs:
//   R_FB^sum(|y|) = (N+ R+ + N- R-) / (N+ + N-)
// 'count(charge, signed-iy)' returns the (fitted) yield of one charge in one
// signed-rapidity bin; weights N+/N- are pooled over the forward bin + backward
// mirror (the same pairing FBratio.C uses). g_RFB_sum supplies the |y| binning.
// gWp/gWm are the per-charge theory graphs {EPPS21,nCTEQ15HQ,nNNPDF30,TUJU21}.
// NOTE: abundance-weighted mean of the two ratios, NOT the exact
// (F+ + F-)/(B+ + B-) identity (they coincide only when R+ = R-).
// -----------------------------------------------------------------------------
// Fetch a TGraphErrors by its discriminant-neutral name, falling back to the
// legacy "_mt"/"_met"-suffixed alias. The suffix never meant m_T for these
// objects: in the production path the yields come from the PF-MET-shape fit
// (see analysis/charge_asym.C and analysis/FBratio.C, which now write both).
static TGraphErrors *GetGraph(TFile *f, const char *name, const char *legacy)
{
    if (!f) return nullptr;
    auto *g = (TGraphErrors *)f->Get(name);
    if (!g) g = (TGraphErrors *)f->Get(legacy);
    return g;
}

static std::vector<TGraphErrors *> buildSumTheorySet(
    std::function<double(const char *, int)> count,
    TGraphErrors *g_RFB_sum,
    TGraphErrors *gWp[4], TGraphErrors *gWm[4])
{
    std::vector<TGraphErrors *> out(4, nullptr);
    if (!g_RFB_sum) { std::cerr << "[WARN] buildSumTheorySet: no g_RFB_sum -> skipped\n"; return out; }

    const int Nabs = g_RFB_sum->GetN(); // |y| bins (6 for NY=12)
    const int NY = 2 * Nabs;            // signed-rapidity bins (12)
    std::vector<double> wP(Nabs, 0.5), wM(Nabs, 0.5), xc(Nabs, 0.0), exc(Nabs, 0.0);

    for (int iabs = 0; iabs < Nabs; ++iabs)
    {
        const int iyB = iabs;          // backward bin (negative y_CM)
        const int iyF = NY - 1 - iabs; // forward mirror (positive y_CM)
        const double Np = count("Wp", iyF) + count("Wp", iyB); // N+ = F+ + B+
        const double Nm = count("Wm", iyF) + count("Wm", iyB); // N- = F- + B-
        const double S = Np + Nm;
        if (S > 0.0) { wP[iabs] = Np / S; wM[iabs] = Nm / S; }
        else std::cerr << "[WARN] zero W+/W- yield in |y| bin " << iabs << " -> equal weights\n";
        xc[iabs] = g_RFB_sum->GetPointX(iabs);
        exc[iabs] = g_RFB_sum->GetErrorX(iabs);
        std::cout << "[INFO] sum-theory weight |y|bin=" << iabs << " (center " << xc[iabs]
                  << ")  N+=" << Np << "  N-=" << Nm
                  << "  w+=" << wP[iabs] << "  w-=" << wM[iabs] << "\n";
    }

    auto findBin = [&](double x) -> int {
        for (int i = 0; i < Nabs; ++i)
            if (x >= xc[i] - exc[i] && x <= xc[i] + exc[i]) return i;
        int best = 0; double bd = 1e30;
        for (int i = 0; i < Nabs; ++i) { const double d = std::fabs(x - xc[i]); if (d < bd) { bd = d; best = i; } }
        return best;
    };

    auto build = [&](TGraphErrors *gP, TGraphErrors *gM, const char *name) -> TGraphErrors * {
        if (!gP || !gM) { std::cerr << "[WARN] " << name << ": missing W+/W- theory -> skipped\n"; return nullptr; }
        const int N = gP->GetN();
        const bool aligned = (gM->GetN() == N);
        if (!aligned) std::cerr << "[WARN] " << name << ": W+/W- theory point mismatch -> Eval-interpolating W-\n";
        std::vector<double> vx, vy, vex, vey;
        for (int i = 0; i < N; ++i) {
            const double x = gP->GetPointX(i);
            const double Rp = gP->GetPointY(i), eRp = gP->GetErrorY(i);
            const double Rm = aligned ? gM->GetPointY(i) : gM->Eval(x);
            const double eRm = aligned ? gM->GetErrorY(i) : 0.0;
            const int b = findBin(x);
            const double Rsum = wP[b] * Rp + wM[b] * Rm;
            const double eRsum = std::sqrt(wP[b] * eRp * wP[b] * eRp + wM[b] * eRm * wM[b] * eRm);
            vx.push_back(x); vy.push_back(Rsum); vex.push_back(gP->GetErrorX(i)); vey.push_back(eRsum);
        }
        TGraphErrors *g = new TGraphErrors((int)vx.size(), vx.data(), vy.data(), vex.data(), vey.data());
        g->SetName(name); g->SetTitle(name);
        return g;
    };

    const char *snames[4] = {"g_RFB_sum_EPPS21", "g_RFB_sum_nCTEQ15HQ", "g_RFB_sum_nNNPDF30", "g_RFB_sum_TUJU21"};
    for (int i = 0; i < 4; ++i) out[i] = build(gWp[i], gWm[i], snames[i]);
    return out;
}

// Read the 4 per-charge nPDF theory graphs (EPPS21,nCTEQ15HQ,nNNPDF30,TUJU21)
// for one boson charge ("WPlus" or "WMinus") into out[4]. Boson-level => the
// same graphs serve muon and electron.
static void readTheoryCharge(TFile *fT, const char *boson, TGraphErrors *out[4])
{
    const char *m[4] = {"EPPS21", "nCTEQ15HQ", "nNNPDF30", "TUJU21"};
    for (int i = 0; i < 4; ++i)
        out[i] = fT ? (TGraphErrors *)fT->Get(
                     Form("pQCDLightIon_mcfm_nuclearmodfactor_%sBoson_RpO_y_dep_%s_FB", boson, m[i]))
                    : nullptr;
}

// shared tuners (y-range + reference line)
static GraphTuner makeTuneCharge()
{
    return [](TCanvas *c, TGraphErrors *g) {
        g->SetMarkerStyle(20); g->SetMarkerSize(1.2); g->SetLineWidth(2);
        if (auto *h = g->GetHistogram()) { h->SetMinimum(-0.6); h->SetMaximum(0.6); }
        double xmin = g->GetXaxis()->GetXmin(), xmax = g->GetXaxis()->GetXmax();
        TLine *l0 = new TLine(xmin, 0.0, xmax, 0.0); l0->SetLineStyle(2); l0->Draw("same");
        c->Modified(); c->Update();
    };
}
static GraphTuner makeTuneRFB()
{
    return [](TCanvas *c, TGraphErrors *g) {
        g->SetMarkerStyle(20); g->SetMarkerSize(1.2); g->SetLineWidth(2);
        if (auto *h = g->GetHistogram()) { h->SetMinimum(0.3); h->SetMaximum(2.1); }
        double xmin = g->GetXaxis()->GetXmin(), xmax = g->GetXaxis()->GetXmax();
        TLine *l1 = new TLine(xmin, 1.0, xmax, 1.0); l1->SetLineStyle(2); l1->Draw("same");
        c->Modified(); c->Update();
    };
}

// Wrap a tuner so it also stamps the discriminant tag as a short 4th header
// line (headerX, headerY - 3*headerDy) -- the plot then self-identifies which
// fit variant produced it. Kept OFF the sub2 line: the long leppt_mt40 label
// would overflow the left-anchored DrawHeader there. The tuner runs after
// DrawHeader and before SaveAs in every SaveNiceGraph* variant, and the spot
// is clear of the bottom-left legends.
static GraphTuner withDiscTag(GraphTuner base, const TString &tagText, const PlotStyle &ps)
{
    const TString tag = tagText;
    const double x = ps.headerX, y = ps.headerY - 3.0 * ps.headerDy;
    const double size = ps.boxTextSize;
    const int font = ps.font;
    return [base, tag, x, y, size, font](TCanvas *c, TGraphErrors *g) {
        if (base) base(c, g);
        TLatex lat; lat.SetNDC(true); lat.SetTextFont(font);
        lat.SetTextAlign(13); lat.SetTextSize(size);
        lat.DrawLatex(x, y, tag.Data());
        c->Modified(); c->Update();
    };
}

// -----------------------------------------------------------------------------
// Shared implementation behind observables() / observables_comb(): everything
// channel-specific arrives as arguments.
//   chan      "mu" | "ele" | "comb"   (comb = the simfit grand fit, 2026-08-04)
//   lepSym    "#mu" | "e" | "l"       (lepton symbol used in the plot titles)
//   outBase   "./plots" | "./plots/Elec" | "./plots/comb"
//   sYieldDef default fitted-yields file (per-flavour out-tree, or the simfit
//             comb_fitted_yields.root) -- supplies the sum-theory weights
//   sub2      3rd header line ("post-fit signal yield" / the simfit label)
// -----------------------------------------------------------------------------
static void observables_run(const char *chan, const char *lepSym,
                            const char *outBase, const TString &sYieldDef,
                            const char *disc, const TString &discLabel,
                            const char *fittedYieldsFile,
                            const char *chargeFile, const char *fbFile,
                            const char *theoryFile, const char *sub2)
{
    gStyle->SetEndErrorSize(4);

    const TString sCharge = chargeFile ? TString(chargeFile)
                            : pODisc::GraphFile("charge_asym", chan, disc);
    const TString sFB = fbFile ? TString(fbFile)
                        : pODisc::GraphFile("FBratio", chan, disc);
    const TString sYield = fittedYieldsFile ? TString(fittedYieldsFile) : sYieldDef;

    TFile *fCharge = TFile::Open(sCharge, "READ");
    if (!fCharge || fCharge->IsZombie())
    { std::cerr << "[ERROR] Cannot open charge-asym file: " << sCharge << "\n        (run charge_asym.C on the fitted-yields file first.)\n"; return; }

    TFile *fFB = TFile::Open(sFB, "READ");
    if (!fFB || fFB->IsZombie())
    { std::cerr << "[ERROR] Cannot open FB-ratio file: " << sFB << "\n        (run FBratio.C on the fitted-yields file first.)\n"; return; }

    TFile *fFB_theory = TFile::Open(theoryFile, "READ");
    if (!fFB_theory || fFB_theory->IsZombie())
    { std::cerr << "[WARN] Theory file not found (" << theoryFile << ") -> plotting data only.\n"; fFB_theory = nullptr; }

    // Per-discriminant output folders so met / leppt / leppt_mt40 coexist.
    const std::string outB(outBase);
    const std::string outFBDir = outB + "/FBratio/" + disc;
    const std::string outChargeDir = outB + "/charge_asym/" + disc;
    gSystem->mkdir(outFBDir.c_str(), kTRUE);
    gSystem->mkdir(outChargeDir.c_str(), kTRUE);

    PlotStyle ps; ps.showStats = false; ps.logy = false;

    // Data graphs. Primary names are discriminant-neutral (g_chargeAsym,
    // g_RFB_*); the "_mt" aliases are legacy only -- they were NEVER an m_T
    // quantity here (nor a mu-tau tag): the yields come from the PF-MET fit.
    // Prefer the discriminant-neutral name; "_mt" is the legacy alias (the
    // yields are from the PF-MET-shape fit, so "_mt" never meant m_T here).
    auto *g_charge = GetGraph(fCharge, "g_chargeAsym", "g_chargeAsym_mt");
    auto *g_RFB_sum = GetGraph(fFB, "g_RFB_sum", "g_RFB_mt_sum");
    auto *g_RFB_Wp = GetGraph(fFB, "g_RFB_Wp", "g_RFB_mt_Wp");
    auto *g_RFB_Wm = GetGraph(fFB, "g_RFB_Wm", "g_RFB_mt_Wm");

    // theory graphs (boson-level -> same for e and mu)
    TGraphErrors *thWp[4], *thWm[4];
    readTheoryCharge(fFB_theory, "WPlus", thWp);
    readTheoryCharge(fFB_theory, "WMinus", thWm);

    if (!g_charge) std::cerr << "[ERROR] Missing g_chargeAsym(_mt) in " << sCharge << "\n";
    if (!g_RFB_sum) std::cerr << "[ERROR] Missing g_RFB_sum(_mt) in " << sFB << "\n";
    if (!g_RFB_Wp) std::cerr << "[ERROR] Missing g_RFB_Wp(_mt) in " << sFB << "\n";
    if (!g_RFB_Wm) std::cerr << "[ERROR] Missing g_RFB_Wm(_mt) in " << sFB << "\n";

    GraphTuner tuneCharge = makeTuneCharge();
    GraphTuner tuneRFB = makeTuneRFB();

    // Discriminant stamped as a 4th header line so a saved plot self-identifies.
    const TString fitTag = TString::Format("%s fit", discLabel.Data());
    GraphTuner tuneChargeTag = withDiscTag(tuneCharge, fitTag, ps);
    GraphTuner tuneRFBTag = withDiscTag(tuneRFB, fitTag, ps);

    // ---- charge asymmetry (data only) ----
    if (g_charge)
        SaveNiceGraph(g_charge, outChargeDir + "/chargeAsym",
                      Form("#eta^{%s}_{CM}", lepSym), "A_{ch}", "",
                      Form("W #rightarrow %s #nu", lepSym), sub2,
                      {}, ps, tuneChargeTag);

    // ---- sum-channel theory (count-weighted by the fitted yields) ----
    std::vector<TGraphErrors *> sumTheory(4, nullptr);
    {
        TFile *fW = TFile::Open(sYield, "READ");
        if (!fW || fW->IsZombie())
            std::cerr << "[WARN] Cannot open fitted-yields file for sum-theory weights: " << sYield << " -> sum theory skipped\n";
        else if (!fFB_theory)
            std::cerr << "[INFO] No theory file -> sum-theory curves skipped (data still plotted).\n";
        else
        {
            using pOAnalysis::YieldInRange;
            auto count = [&](const char *chg, int iy) -> double {
                // primary fitted-yield name; h_mt_* is the deprecated alias
                TH1D *h = (TH1D *)fW->Get(Form("h_yield_%s_y%d_FB", chg, iy));
                if (!h) h = (TH1D *)fW->Get(Form("h_mt_%s_y%d_FB", chg, iy));
                return YieldInRange(h, 30.0, 200.0, true).value;
            };
            sumTheory = buildSumTheorySet(count, g_RFB_sum, thWp, thWm);
        }
        if (fW) { fW->Close(); delete fW; }
    }

    // ---- R_FB plots (sum, W+, W-) with all-4-model bands ----
    auto plotRFB = [&](TGraphErrors *g, const std::string &tag, const std::string &subtitle,
                       TGraphErrors *t1, TGraphErrors *t2, TGraphErrors *t3, TGraphErrors *t4) {
        if (!g) return;
        SaveNiceGraph_ErrorBand(g, outFBDir + "/" + tag, Form("#eta^{%s}_{CM}", lepSym), "R_{FB}",
                                "", subtitle, sub2, {}, ps, tuneRFBTag, t1, t2, t3, t4);
    };

    plotRFB(g_RFB_sum, "RFB_sum", Form("W #rightarrow %s #nu", lepSym),
            sumTheory[0], sumTheory[1], sumTheory[2], sumTheory[3]);
    plotRFB(g_RFB_Wp, "RFB_Wp", Form("W^{+} #rightarrow %s^{+} #nu", lepSym),
            thWp[0], thWp[1], thWp[2], thWp[3]);
    plotRFB(g_RFB_Wm, "RFB_Wm", Form("W^{-} #rightarrow %s^{-} #bar{#nu}", lepSym),
            thWm[0], thWm[1], thWm[2], thWm[3]);

    fCharge->Close(); fFB->Close(); delete fCharge; delete fFB;
    if (fFB_theory) { fFB_theory->Close(); delete fFB_theory; }
    std::cout << "[OK] Saved " << chan << " observables (disc=" << disc << ") to: "
              << outChargeDir << " and " << outFBDir << "\n";
}

// =============================================================================
// observables -- per-flavour plots (legacy per-bin fit out-trees).
// =============================================================================
void observables(bool isElec = false,
                 const char *disc = "met",
                 const char *fittedYieldsFile = nullptr,
                 const char *chargeFile = nullptr,
                 const char *fbFile = nullptr,
                 const char *theoryFile = "./RpO_rootfile/RpO_FB_graphs.root")
{
    TString dsuf, discLabel;
    if (!pODisc::Spec(disc, dsuf, discLabel)) return;

    const char *chan = isElec ? "ele" : "mu";
    const TString sYieldDef = TString::Format(
        "../../HiggsAnalysis-CombinedLimit/test/pO_fit_out%s/%s/summary/%s_fitted_yields.root",
        dsuf.Data(), chan, chan);
    observables_run(chan, isElec ? "e" : "#mu",
                    isElec ? "./plots/Elec" : "./plots", sYieldDef,
                    disc, discLabel, fittedYieldsFile, chargeFile, fbFile,
                    theoryFile, "post-fit signal yield");
}

// =============================================================================
// observables_comb -- the PRIMARY observables of the simfit GRAND SIMULTANEOUS
// FIT (2026-08-04): mu/e-shared per-(charge, y-bin) signal strengths r_<C>_y<i>
// + global r_Z, one likelihood. Yields (and their covariance, used upstream in
// charge_asym.C / FBratio.C) come from the fork's
// pO_fit_out<suffix>/simfit/summary/comb_fitted_yields.root.
// Outputs: ./plots/comb/{charge_asym,FBratio}/<disc>/.
// NB with the shared r's the per-flavour observables of observables() are 100%
// correlated with these -- the comb plots are the result; per-flavour ones are
// only the legacy per-bin-fit cross-check.
// =============================================================================
void observables_comb(const char *disc = "met",
                      const char *fittedYieldsFile = nullptr,
                      const char *chargeFile = nullptr,
                      const char *fbFile = nullptr,
                      const char *theoryFile = "./RpO_rootfile/RpO_FB_graphs.root")
{
    TString dsuf, discLabel;
    if (!pODisc::Spec(disc, dsuf, discLabel)) return;

    const TString sYieldDef = TString::Format(
        "../../HiggsAnalysis-CombinedLimit/test/pO_fit_out%s/simfit/summary/comb_fitted_yields.root",
        dsuf.Data());
    observables_run("comb", "l", "./plots/comb", sYieldDef,
                    disc, discLabel, fittedYieldsFile, chargeFile, fbFile,
                    theoryFile, "#mu + e combined (simfit)");
}

// =============================================================================
// Merged muon + electron overlay (keeps the individual per-channel plots from
// observables(false)/(true); this just adds the combined view).  Run AFTER you
// have produced both channels' charge_asym_fit_<chan>_<disc>.root +
// FBratio_fit_<chan>_<disc>.root for the SAME disc.
//   disc                met|leppt|leppt_mt40 (drives defaults + output folder)
//   muYields/eleYields  fitted-yields files (combined abundance weights for sum)
//   *Charge/*FB         per-channel charge_asym/FBratio outputs (defaults below)
//   theoryFile          RpO_FB_graphs.root (missing -> data-only overlays)
// Output: ./plots/merged/<disc>/{chargeAsym,RFB_sum,RFB_Wp,RFB_Wm}_overlay.{png,pdf}
// =============================================================================
void observables_overlay(const char *disc = "met",
                         const char *muYields = nullptr,
                         const char *eleYields = nullptr,
                         const char *muCharge = nullptr,
                         const char *muFB = nullptr,
                         const char *eleCharge = nullptr,
                         const char *eleFB = nullptr,
                         const char *theoryFile = "./RpO_rootfile/RpO_FB_graphs.root")
{
    gStyle->SetEndErrorSize(4);

    TString dsuf, discLabel;
    if (!pODisc::Spec(disc, dsuf, discLabel)) return;

    const TString sMuCharge  = muCharge  ? TString(muCharge)  : pODisc::GraphFile("charge_asym", "mu", disc);
    const TString sMuFB      = muFB      ? TString(muFB)      : pODisc::GraphFile("FBratio", "mu", disc);
    const TString sElCharge  = eleCharge ? TString(eleCharge) : pODisc::GraphFile("charge_asym", "ele", disc);
    const TString sElFB      = eleFB     ? TString(eleFB)     : pODisc::GraphFile("FBratio", "ele", disc);
    const TString sMuYield   = muYields  ? TString(muYields)  : TString::Format("../../HiggsAnalysis-CombinedLimit/test/pO_fit_out%s/mu/summary/mu_fitted_yields.root", dsuf.Data());
    const TString sElYield   = eleYields ? TString(eleYields) : TString::Format("../../HiggsAnalysis-CombinedLimit/test/pO_fit_out%s/ele/summary/ele_fitted_yields.root", dsuf.Data());

    TFile *fCmu = TFile::Open(sMuCharge, "READ");
    TFile *fCel = TFile::Open(sElCharge, "READ");
    TFile *fFmu = TFile::Open(sMuFB, "READ");
    TFile *fFel = TFile::Open(sElFB, "READ");
    auto bad = [](TFile *f) { return !f || f->IsZombie(); };
    if (bad(fCmu) || bad(fCel) || bad(fFmu) || bad(fFel))
    {
        std::cerr << "[ERROR] overlay needs BOTH channels' charge_asym + FBratio outputs\n"
                  << "        for disc=" << disc << ". Run charge_asym.C / FBratio.C on both\n"
                  << "        channels' fitted yields first -- easiest via\n"
                  << "        analysis/run_observables.sh " << disc << " (README Module 5).\n";
        return;
    }

    TFile *fT = TFile::Open(theoryFile, "READ");
    if (!fT || fT->IsZombie()) { std::cerr << "[WARN] no theory file -> data-only overlays\n"; fT = nullptr; }

    const std::string outDir = std::string("./plots/merged/") + disc;
    gSystem->mkdir(outDir.c_str(), kTRUE);

    PlotStyle ps; ps.showStats = false; ps.logy = false;
    GraphTuner tuneCharge = makeTuneCharge();
    GraphTuner tuneRFB = makeTuneRFB();

    // data graphs, both channels
    auto *gC_mu = GetGraph(fCmu, "g_chargeAsym", "g_chargeAsym_mt");
    auto *gC_el = GetGraph(fCel, "g_chargeAsym", "g_chargeAsym_mt");
    auto *gS_mu = GetGraph(fFmu, "g_RFB_sum", "g_RFB_mt_sum");
    auto *gS_el = GetGraph(fFel, "g_RFB_sum", "g_RFB_mt_sum");
    auto *gP_mu = GetGraph(fFmu, "g_RFB_Wp", "g_RFB_mt_Wp");
    auto *gP_el = GetGraph(fFel, "g_RFB_Wp", "g_RFB_mt_Wp");
    auto *gM_mu = GetGraph(fFmu, "g_RFB_Wm", "g_RFB_mt_Wm");
    auto *gM_el = GetGraph(fFel, "g_RFB_Wm", "g_RFB_mt_Wm");

    // theory (shared, boson-level)
    TGraphErrors *thWp[4], *thWm[4];
    readTheoryCharge(fT, "WPlus", thWp);
    readTheoryCharge(fT, "WMinus", thWm);

    // sum-channel theory weighted by the COMBINED (mu+ele) fitted yields
    std::vector<TGraphErrors *> sumTheory(4, nullptr);
    {
        TFile *fYmu = TFile::Open(sMuYield, "READ");
        TFile *fYel = TFile::Open(sElYield, "READ");
        TGraphErrors *gS_ref = gS_mu ? gS_mu : gS_el;
        if (fT && gS_ref && !bad(fYmu) && !bad(fYel))
        {
            using pOAnalysis::YieldInRange;
            // primary fitted-yield name; h_mt_* is the deprecated alias
            auto getY = [](TFile *fy, const char *chg, int iy) -> TH1D * {
                TH1D *h = (TH1D *)fy->Get(Form("h_yield_%s_y%d_FB", chg, iy));
                if (!h) h = (TH1D *)fy->Get(Form("h_mt_%s_y%d_FB", chg, iy));
                return h;
            };
            auto count = [&](const char *chg, int iy) -> double {
                return YieldInRange(getY(fYmu, chg, iy), 30.0, 200.0, true).value
                     + YieldInRange(getY(fYel, chg, iy), 30.0, 200.0, true).value;
            };
            sumTheory = buildSumTheorySet(count, gS_ref, thWp, thWm);
        }
        else
            std::cerr << "[WARN] combined fitted-yields unavailable -> sum-overlay theory skipped\n";
        if (fYmu) { fYmu->Close(); delete fYmu; }
        if (fYel) { fYel->Close(); delete fYel; }
    }

    const char *muLab = "Muon", *elLab = "Electron";
    // Discriminant stamped as a 4th header line so a saved plot self-identifies.
    const TString fitTag = TString::Format("%s fit", discLabel.Data());
    GraphTuner tuneChargeTag = withDiscTag(tuneCharge, fitTag, ps);
    GraphTuner tuneRFBTag = withDiscTag(tuneRFB, fitTag, ps);

    // charge asymmetry overlay (no theory bands)
    SaveNiceGraph_ErrorBand_TwoData(gC_mu, muLab, gC_el, elLab,
        outDir + "/chargeAsym_overlay", "#eta_{CM}", "A_{ch}", "",
        "W charge asymmetry", "#mu + e (post-fit)", {}, ps, tuneChargeTag);

    // R_FB overlays (data both channels + all-4-model bands)
    SaveNiceGraph_ErrorBand_TwoData(gS_mu, muLab, gS_el, elLab,
        outDir + "/RFB_sum_overlay", "#eta_{CM}", "R_{FB}", "",
        "W #rightarrow l #nu", "#mu + e (post-fit)", {}, ps, tuneRFBTag,
        sumTheory[0], sumTheory[1], sumTheory[2], sumTheory[3]);

    SaveNiceGraph_ErrorBand_TwoData(gP_mu, muLab, gP_el, elLab,
        outDir + "/RFB_Wp_overlay", "#eta_{CM}", "R_{FB}", "",
        "W^{+} #rightarrow l^{+} #nu", "#mu + e (post-fit)", {}, ps, tuneRFBTag,
        thWp[0], thWp[1], thWp[2], thWp[3]);

    SaveNiceGraph_ErrorBand_TwoData(gM_mu, muLab, gM_el, elLab,
        outDir + "/RFB_Wm_overlay", "#eta_{CM}", "R_{FB}", "",
        "W^{-} #rightarrow l^{-} #bar{#nu}", "#mu + e (post-fit)", {}, ps, tuneRFBTag,
        thWm[0], thWm[1], thWm[2], thWm[3]);

    fCmu->Close(); fCel->Close(); fFmu->Close(); fFel->Close();
    delete fCmu; delete fCel; delete fFmu; delete fFel;
    if (fT) { fT->Close(); delete fT; }
    std::cout << "[OK] Saved merged mu+ele overlays (disc=" << disc << ") to: " << outDir << "\n";
}
