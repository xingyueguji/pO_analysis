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

// =============================================================================
// observables.C -- final charge-asymmetry + forward/backward plots from the
// FITTED signal yields (Combine).
//
//   observables(isElec)        per-channel plots (data + all-4-model theory bands)
//   observables_overlay()      MERGED muon+electron overlay on the same axes
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

// =============================================================================
// Per-channel plots (data + all-4-model theory bands).
//   isElec           false=muon, true=electron
//   fittedYieldsFile <chan>_fitted_yields.root (sum-theory abundance weights)
//   chargeFile/fbFile  charge_asym.C / FBratio.C outputs (default per-channel)
//   theoryFile       RpO_FB_graphs.root (missing -> data-only)
// =============================================================================
void observables(bool isElec = false,
                 const char *fittedYieldsFile = nullptr,
                 const char *chargeFile = nullptr,
                 const char *fbFile = nullptr,
                 const char *theoryFile = "./RpO_rootfile/RpO_FB_graphs.root")
{
    gStyle->SetEndErrorSize(4);

    const char *chan   = isElec ? "ele" : "mu";
    const char *lepSym = isElec ? "e" : "#mu";

    const TString sCharge = chargeFile ? TString(chargeFile)
                            : TString::Format("../skim/rootfile/charge_asym_fit_%s.root", chan);
    const TString sFB = fbFile ? TString(fbFile)
                        : TString::Format("../skim/rootfile/FBratio_fit_%s.root", chan);
    const TString sYield = fittedYieldsFile ? TString(fittedYieldsFile)
                           : TString::Format("../../HiggsAnalysis-CombinedLimit/test/pO_fit_out/%s/summary/%s_fitted_yields.root", chan, chan);

    TFile *fCharge = TFile::Open(sCharge, "READ");
    if (!fCharge || fCharge->IsZombie())
    { std::cerr << "[ERROR] Cannot open charge-asym file: " << sCharge << "\n        (run charge_asym.C on the fitted-yields file first.)\n"; return; }

    TFile *fFB = TFile::Open(sFB, "READ");
    if (!fFB || fFB->IsZombie())
    { std::cerr << "[ERROR] Cannot open FB-ratio file: " << sFB << "\n        (run FBratio.C on the fitted-yields file first.)\n"; return; }

    TFile *fFB_theory = TFile::Open(theoryFile, "READ");
    if (!fFB_theory || fFB_theory->IsZombie())
    { std::cerr << "[WARN] Theory file not found (" << theoryFile << ") -> plotting data only.\n"; fFB_theory = nullptr; }

    const std::string outBase = isElec ? "./plots/Elec" : "./plots";
    const std::string outFBDir = outBase + "/FBratio";
    const std::string outChargeDir = outBase + "/charge_asym";
    gSystem->mkdir(outFBDir.c_str(), kTRUE);
    gSystem->mkdir(outChargeDir.c_str(), kTRUE);

    PlotStyle ps; ps.showStats = false; ps.logy = false;

    // data graphs (names channel-independent: "_mt" = useMT flag, not flavour)
    auto *g_charge = (TGraphErrors *)fCharge->Get("g_chargeAsym_mt");
    auto *g_RFB_sum = (TGraphErrors *)fFB->Get("g_RFB_mt_sum");
    auto *g_RFB_Wp = (TGraphErrors *)fFB->Get("g_RFB_mt_Wp");
    auto *g_RFB_Wm = (TGraphErrors *)fFB->Get("g_RFB_mt_Wm");

    // theory graphs (boson-level -> same for e and mu)
    TGraphErrors *thWp[4], *thWm[4];
    readTheoryCharge(fFB_theory, "WPlus", thWp);
    readTheoryCharge(fFB_theory, "WMinus", thWm);

    if (!g_charge) std::cerr << "[ERROR] Missing g_chargeAsym_mt in " << sCharge << "\n";
    if (!g_RFB_sum) std::cerr << "[ERROR] Missing g_RFB_mt_sum in " << sFB << "\n";
    if (!g_RFB_Wp) std::cerr << "[ERROR] Missing g_RFB_mt_Wp in " << sFB << "\n";
    if (!g_RFB_Wm) std::cerr << "[ERROR] Missing g_RFB_mt_Wm in " << sFB << "\n";

    GraphTuner tuneCharge = makeTuneCharge();
    GraphTuner tuneRFB = makeTuneRFB();

    // ---- charge asymmetry (data only) ----
    if (g_charge)
        SaveNiceGraph(g_charge, outChargeDir + "/chargeAsym_mt",
                      Form("#eta^{%s}_{CM}", lepSym), "A_{ch}", "",
                      Form("W #rightarrow %s #nu", lepSym), "post-fit signal yield",
                      {}, ps, tuneCharge);

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
                TH1D *h = (TH1D *)fW->Get(Form("h_mt_%s_y%d_FB", chg, iy));
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
                                "", subtitle, "post-fit signal yield", {}, ps, tuneRFB, t1, t2, t3, t4);
    };

    plotRFB(g_RFB_sum, "RFB_mt_sum", Form("W #rightarrow %s #nu", lepSym),
            sumTheory[0], sumTheory[1], sumTheory[2], sumTheory[3]);
    plotRFB(g_RFB_Wp, "RFB_mt_Wp", Form("W^{+} #rightarrow %s^{+} #nu", lepSym),
            thWp[0], thWp[1], thWp[2], thWp[3]);
    plotRFB(g_RFB_Wm, "RFB_mt_Wm", Form("W^{-} #rightarrow %s^{-} #bar{#nu}", lepSym),
            thWm[0], thWm[1], thWm[2], thWm[3]);

    fCharge->Close(); fFB->Close(); delete fCharge; delete fFB;
    if (fFB_theory) { fFB_theory->Close(); delete fFB_theory; }
    std::cout << "[OK] Saved " << chan << " observables to: " << outBase << "\n";
}

// =============================================================================
// Merged muon + electron overlay (keeps the individual per-channel plots from
// observables(false)/(true); this just adds the combined view).  Run AFTER you
// have produced both channels' charge_asym_fit_<chan>.root + FBratio_fit_<chan>.root.
//   muYields/eleYields  fitted-yields files (combined abundance weights for sum)
//   *Charge/*FB         per-channel charge_asym/FBratio outputs (defaults below)
//   theoryFile          RpO_FB_graphs.root (missing -> data-only overlays)
// Output: ./plots/merged/{chargeAsym,RFB_sum,RFB_Wp,RFB_Wm}_overlay.{png,pdf}
// =============================================================================
void observables_overlay(const char *muYields = nullptr,
                         const char *eleYields = nullptr,
                         const char *muCharge = nullptr,
                         const char *muFB = nullptr,
                         const char *eleCharge = nullptr,
                         const char *eleFB = nullptr,
                         const char *theoryFile = "./RpO_rootfile/RpO_FB_graphs.root")
{
    gStyle->SetEndErrorSize(4);

    const TString sMuCharge  = muCharge  ? TString(muCharge)  : "../skim/rootfile/charge_asym_fit_mu.root";
    const TString sMuFB      = muFB      ? TString(muFB)      : "../skim/rootfile/FBratio_fit_mu.root";
    const TString sElCharge  = eleCharge ? TString(eleCharge) : "../skim/rootfile/charge_asym_fit_ele.root";
    const TString sElFB      = eleFB     ? TString(eleFB)     : "../skim/rootfile/FBratio_fit_ele.root";
    const TString sMuYield   = muYields  ? TString(muYields)  : "../../HiggsAnalysis-CombinedLimit/test/pO_fit_out/mu/summary/mu_fitted_yields.root";
    const TString sElYield   = eleYields ? TString(eleYields) : "../../HiggsAnalysis-CombinedLimit/test/pO_fit_out/ele/summary/ele_fitted_yields.root";

    TFile *fCmu = TFile::Open(sMuCharge, "READ");
    TFile *fCel = TFile::Open(sElCharge, "READ");
    TFile *fFmu = TFile::Open(sMuFB, "READ");
    TFile *fFel = TFile::Open(sElFB, "READ");
    auto bad = [](TFile *f) { return !f || f->IsZombie(); };
    if (bad(fCmu) || bad(fCel) || bad(fFmu) || bad(fFel))
    {
        std::cerr << "[ERROR] overlay needs BOTH channels' charge_asym + FBratio outputs.\n"
                  << "        Run charge_asym.C / FBratio.C on mu_fitted_yields.root AND\n"
                  << "        ele_fitted_yields.root first (see README Step 6).\n";
        return;
    }

    TFile *fT = TFile::Open(theoryFile, "READ");
    if (!fT || fT->IsZombie()) { std::cerr << "[WARN] no theory file -> data-only overlays\n"; fT = nullptr; }

    const std::string outDir = "./plots/merged";
    gSystem->mkdir(outDir.c_str(), kTRUE);

    PlotStyle ps; ps.showStats = false; ps.logy = false;
    GraphTuner tuneCharge = makeTuneCharge();
    GraphTuner tuneRFB = makeTuneRFB();

    // data graphs, both channels
    auto *gC_mu = (TGraphErrors *)fCmu->Get("g_chargeAsym_mt");
    auto *gC_el = (TGraphErrors *)fCel->Get("g_chargeAsym_mt");
    auto *gS_mu = (TGraphErrors *)fFmu->Get("g_RFB_mt_sum");
    auto *gS_el = (TGraphErrors *)fFel->Get("g_RFB_mt_sum");
    auto *gP_mu = (TGraphErrors *)fFmu->Get("g_RFB_mt_Wp");
    auto *gP_el = (TGraphErrors *)fFel->Get("g_RFB_mt_Wp");
    auto *gM_mu = (TGraphErrors *)fFmu->Get("g_RFB_mt_Wm");
    auto *gM_el = (TGraphErrors *)fFel->Get("g_RFB_mt_Wm");

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
            auto count = [&](const char *chg, int iy) -> double {
                double s = 0.0;
                TH1D *h = (TH1D *)fYmu->Get(Form("h_mt_%s_y%d_FB", chg, iy));
                s += YieldInRange(h, 30.0, 200.0, true).value;
                h = (TH1D *)fYel->Get(Form("h_mt_%s_y%d_FB", chg, iy));
                s += YieldInRange(h, 30.0, 200.0, true).value;
                return s;
            };
            sumTheory = buildSumTheorySet(count, gS_ref, thWp, thWm);
        }
        else
            std::cerr << "[WARN] combined fitted-yields unavailable -> sum-overlay theory skipped\n";
        if (fYmu) { fYmu->Close(); delete fYmu; }
        if (fYel) { fYel->Close(); delete fYel; }
    }

    const char *muLab = "Muon", *elLab = "Electron";

    // charge asymmetry overlay (no theory bands)
    SaveNiceGraph_ErrorBand_TwoData(gC_mu, muLab, gC_el, elLab,
        outDir + "/chargeAsym_overlay", "#eta_{CM}", "A_{ch}", "",
        "W charge asymmetry", "post-fit signal yield, #mu + e", {}, ps, tuneCharge);

    // R_FB overlays (data both channels + all-4-model bands)
    SaveNiceGraph_ErrorBand_TwoData(gS_mu, muLab, gS_el, elLab,
        outDir + "/RFB_sum_overlay", "#eta_{CM}", "R_{FB}", "",
        "W #rightarrow l #nu", "post-fit signal yield, #mu + e", {}, ps, tuneRFB,
        sumTheory[0], sumTheory[1], sumTheory[2], sumTheory[3]);

    SaveNiceGraph_ErrorBand_TwoData(gP_mu, muLab, gP_el, elLab,
        outDir + "/RFB_Wp_overlay", "#eta_{CM}", "R_{FB}", "",
        "W^{+} #rightarrow l^{+} #nu", "post-fit signal yield, #mu + e", {}, ps, tuneRFB,
        thWp[0], thWp[1], thWp[2], thWp[3]);

    SaveNiceGraph_ErrorBand_TwoData(gM_mu, muLab, gM_el, elLab,
        outDir + "/RFB_Wm_overlay", "#eta_{CM}", "R_{FB}", "",
        "W^{-} #rightarrow l^{-} #bar{#nu}", "post-fit signal yield, #mu + e", {}, ps, tuneRFB,
        thWm[0], thWm[1], thWm[2], thWm[3]);

    fCmu->Close(); fCel->Close(); fFmu->Close(); fFel->Close();
    delete fCmu; delete fCel; delete fFmu; delete fFel;
    if (fT) { fT->Close(); delete fT; }
    std::cout << "[OK] Saved merged mu+ele overlays to: " << outDir << "\n";
}
