#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <functional>

#include "plotting_helper.C"               // PlotStyle, SaveNiceGraph[_ErrorBand|_Overlay]
#include "../analysis/analysis_helpers.h"  // pOAnalysis::YieldInRange (Sumw2-aware yields)
#include "disc_variants.h"                 // pODisc::Spec/GraphFile (W-discriminant tags)
#include "fit_variants.h"                  // pOFit::Spec (which fit: comb / mu-only / e-only)
#include "TMath.h"

// =============================================================================
// observables.C -- final charge-asymmetry + forward/backward plots from the
// FITTED signal yields (Combine).
//
//   observables_comb(disc)           PRIMARY: the simfit grand fit's mu+e-combined
//                                    observables (plots/comb/..., 2026-08-04)
//   observables_flav(disc)           the MU-ONLY and E-ONLY simultaneous fits
//                                    (fork mode flavfit, 2026-09-22) OVERLAID,
//                                    each with stat bars + syst boxes
//                                    (plots/flavfit/...), + a mu-vs-e chi2
//   observables(disc)                both of the above
//
// `disc` = met|leppt_mt40 selects WHICH fit's yields are plotted
// (2026-08-03): it drives every default input path (the fork out-tree
// pO_fit_out<suffix>/, the tagged charge_asym/FBratio_fid_<fit>_<disc>.root --
// FIDUCIAL, r x sigma_gen, never raw counts: user directive 2026-09-22)
// AND the output folders (plots/{comb,flavfit}/{charge_asym,FBratio}/<disc>/),
// so the variants coexist without overwriting. The discriminant is also
// stamped into the plot info box. One-command runner for the whole chain:
// analysis/run_observables.sh (README Module 5).
//
// All theory bands show ALL FOUR nPDF sets (EPPS21, nCTEQ15HQ, nNNPDF3.0,
// TUJU21nlo). The LEGACY per-flavour plots (observables(isElec, disc) and the
// mu+e observables_overlay, from the per-bin legacy fits) were removed with
// those fits on 2026-09-22 -- observables_flav is the mu-vs-e comparison now.
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
// Implementation behind observables_comb(): everything fit-specific arrives as
// arguments (it served the legacy per-flavour plots too until 2026-09-22).
//   chan      the fit tag in the graph file names ("comb")
//   lepSym    "l"                     (lepton symbol used in the plot titles)
//   outBase   "./plots/comb"
//   sYieldDef default yields file = the fit's FIDUCIAL yields r x sigma_gen,
//             ../skim/rootfile/fidyields_comb_<disc>.root (analysis/
//             fiducial_yields.C) -- supplies the sum-theory weights
//   sub2      3rd header line (the simfit label)
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
    { std::cerr << "[ERROR] Cannot open charge-asym file: " << sCharge << "\n        (analysis/run_observables.sh makes it: fiducial_yields.C -> charge_asym.C.)\n"; return; }

    TFile *fFB = TFile::Open(sFB, "READ");
    if (!fFB || fFB->IsZombie())
    { std::cerr << "[ERROR] Cannot open FB-ratio file: " << sFB << "\n        (analysis/run_observables.sh makes it: fiducial_yields.C -> FBratio.C.)\n"; return; }

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

    // Statistical twins (2026-09-15): same points, statistical error only,
    // written by charge_asym.C / FBratio.C whenever their input yields file
    // carries the stat covariance h_cov_yield[_FB]_stat (fiducial_yields.C
    // writes it from the fit's h_cov_poi[_FB]_stat). When present, the
    // plotted error bars become the STATISTICAL ones and the systematic --
    // sqrt(total^2 - stat^2) -- is drawn as a TBox per point. Absent (raw
    // skim files, pre-2026-09-14 extractions) -> a single total-error bar,
    // exactly as before.
    auto *g_charge_stat = (TGraphErrors *)fCharge->Get("g_chargeAsym_stat");
    auto *g_RFB_sum_stat = (TGraphErrors *)fFB->Get("g_RFB_sum_stat");
    auto *g_RFB_Wp_stat = (TGraphErrors *)fFB->Get("g_RFB_Wp_stat");
    auto *g_RFB_Wm_stat = (TGraphErrors *)fFB->Get("g_RFB_Wm_stat");
    if (!g_charge_stat && !g_RFB_sum_stat)
        std::cerr << "[WARN] no *_stat graphs in " << sCharge << " / " << sFB
                  << "\n        -> single total-error bars (no systematic boxes)."
                  << "\n        Re-run the fork extraction (run_pO_fits.sh --extract-only) so the"
                  << " fit's <tag>_fitted_yields.root carries the stat covariance, then run_observables.sh.\n";

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
                      {}, ps, tuneChargeTag,
                      nullptr, nullptr, nullptr, nullptr, g_charge_stat);

    // ---- sum-channel theory (W+/W- weighted by the fit's FIDUCIAL yields r x sigma_gen) ----
    std::vector<TGraphErrors *> sumTheory(4, nullptr);
    {
        TFile *fW = TFile::Open(sYield, "READ");
        if (!fW || fW->IsZombie())
            std::cerr << "[WARN] Cannot open fiducial-yields file for sum-theory weights: " << sYield << " -> sum theory skipped\n";
        else if (!fFB_theory)
            std::cerr << "[INFO] No theory file -> sum-theory curves skipped (data still plotted).\n";
        else
        {
            using pOAnalysis::YieldInRange;
            auto count = [&](const char *chg, int iy) -> double {
                // h_yield_* (the name fiducial_yields.C writes); h_mt_* is the deprecated alias
                TH1D *h = (TH1D *)fW->Get(Form("h_yield_%s_y%d_FB", chg, iy));
                if (!h) h = (TH1D *)fW->Get(Form("h_mt_%s_y%d_FB", chg, iy));
                return YieldInRange(h, 30.0, 200.0, true).value;
            };
            sumTheory = buildSumTheorySet(count, g_RFB_sum, thWp, thWm);
        }
        if (fW) { fW->Close(); delete fW; }
    }

    // ---- R_FB plots (sum, W+, W-) with all-4-model bands ----
    auto plotRFB = [&](TGraphErrors *g, TGraphErrors *gStat,
                       const std::string &tag, const std::string &subtitle,
                       TGraphErrors *t1, TGraphErrors *t2, TGraphErrors *t3, TGraphErrors *t4) {
        if (!g) return;
        SaveNiceGraph_ErrorBand(g, outFBDir + "/" + tag, Form("#eta^{%s}_{CM}", lepSym), "R_{FB}",
                                "", subtitle, sub2, {}, ps, tuneRFBTag, t1, t2, t3, t4, gStat);
    };

    plotRFB(g_RFB_sum, g_RFB_sum_stat, "RFB_sum", Form("W #rightarrow %s #nu", lepSym),
            sumTheory[0], sumTheory[1], sumTheory[2], sumTheory[3]);
    plotRFB(g_RFB_Wp, g_RFB_Wp_stat, "RFB_Wp", Form("W^{+} #rightarrow %s^{+} #nu", lepSym),
            thWp[0], thWp[1], thWp[2], thWp[3]);
    plotRFB(g_RFB_Wm, g_RFB_Wm_stat, "RFB_Wm", Form("W^{-} #rightarrow %s^{-} #bar{#nu}", lepSym),
            thWm[0], thWm[1], thWm[2], thWm[3]);

    fCharge->Close(); fFB->Close(); delete fCharge; delete fFB;
    if (fFB_theory) { fFB_theory->Close(); delete fFB_theory; }
    std::cout << "[OK] Saved " << chan << " observables (disc=" << disc << ") to: "
              << outChargeDir << " and " << outFBDir << "\n";
}

// =============================================================================
// observables_comb -- the PRIMARY observables of the simfit GRAND SIMULTANEOUS
// FIT (2026-08-04): mu/e-shared per-(charge, y-bin) signal strengths r_<C>_y<i>
// + global r_Z, one likelihood. The r's and their covariance come from the
// fork's pO_fit_out<suffix>/simfit/summary/comb_{W_yields.csv,fitted_yields.root}
// and are turned into r x sigma_gen by analysis/fiducial_yields.C (below).
// Outputs: ./plots/comb/{charge_asym,FBratio}/<disc>/.
// NB with the shared r's the mu and e observables are 100% correlated inside
// this fit -- the independent per-flavour results are observables_flav's.
// FIDUCIAL since 2026-09-22 (user: every observable from r x sigma_gen, never
// raw counts): the graphs are charge_asym.C / FBratio.C on the fit's
// fiducial yields ../skim/rootfile/fidyields_comb_<disc>.root
// (analysis/fiducial_yields.C), which also supply the sum-theory weights.
// =============================================================================
void observables_comb(const char *disc = "met",
                      const char *fittedYieldsFile = nullptr,
                      const char *chargeFile = nullptr,
                      const char *fbFile = nullptr,
                      const char *theoryFile = "./RpO_rootfile/RpO_FB_graphs.root")
{
    TString dsuf, discLabel;
    if (!pODisc::Spec(disc, dsuf, discLabel)) return;

    const TString sYieldDef = TString::Format("../skim/rootfile/fidyields_comb_%s.root", disc);
    observables_run("comb", "l", "./plots/comb", sYieldDef,
                    disc, discLabel, fittedYieldsFile, chargeFile, fbFile,
                    theoryFile, "#mu + e simfit, fiducial");
}

// =============================================================================
// observables_flav -- the MU-ONLY and E-ONLY simultaneous fits (fork
// run_pO_fits.sh mode flavfit, 2026-09-22) OVERLAID: charge asymmetry and R_FB
// (sum, W+, W-). Every fit is drawn exactly as observables_comb draws the grand
// one -- point + statistical bar + the systematic, sqrt(total^2 - stat^2), as a
// box -- from the charge_asym.C / FBratio.C graphs and their *_stat twins, side
// by side in each bin (mu left, e right; the grand fit in the middle with
// withComb, off by default: in these per-bin plots it mostly crowds the bins,
// and it has its own plots in plots/comb/).
//
// FIDUCIAL (acceptance-corrected) inputs, NOT the count-based ones: the graphs
// are charge_asym.C / FBratio.C run on analysis/fiducial_yields.C's
// r_i x sigma_gen-fid,i yields (with the full r covariance, total + stat).
// The count-based R_FB of the two flavours is NOT comparable -- F and B are
// different |eta_lab| regions, so each flavour's A x eps (the electron ECAL
// crack above all) enters the ratio and does not cancel: on IDENTICAL r's the
// mu and e count-based R_FB differ by up to 60% (see fiducial_yields.C).
// In one fiducial volume the only difference left is the data.
//
// COMPARING THEM: the mu and e STATISTICAL errors are independent (disjoint
// event samples); their systematics are not -- lumi and the theory shapes act
// coherently on both (lumi cancels in A and R_FB anyway), muSF and the QCD
// nuisances on one flavour only. The console therefore prints, per observable,
//   chi2 = Sum_i (x_mu,i - x_e,i)^2 / (stat_mu,i^2 + stat_e,i^2)
// and the same with the TOTAL errors (bins treated as independent): the first
// is the tension against statistics alone, the second a conservative bound
// (it counts the common systematics as if they were independent). Same numbers
// in plots/flavfit/{charge_asym,FBratio}/<disc>/mu_vs_e_chi2.csv.
//
// Inputs:  ../skim/rootfile/{charge_asym,FBratio}_fid_simfit_{mu,ele}_<disc>.root
//          [+ _fid_comb_<disc>] (pODisc::GraphFile), made by
//          analysis/run_observables.sh (fiducial_yields.C -> charge_asym.C /
//          FBratio.C), and for the sum-theory weights
//          ../skim/rootfile/fidyields_<tag>_<disc>.root
// Outputs: ./plots/flavfit/charge_asym/<disc>/chargeAsym.{png,pdf}
//          ./plots/flavfit/FBratio/<disc>/RFB_{sum,Wp,Wm}.{png,pdf}
// =============================================================================
void observables_flav(const char *disc = "leppt_mt40",
                      bool withComb = false,
                      const char *theoryFile = "./RpO_rootfile/RpO_FB_graphs.root")
{
    gStyle->SetEndErrorSize(4);

    TString dsuf, discLabel;
    if (!pODisc::Spec(disc, dsuf, discLabel)) return;

    // the fits to overlay, in drawing order left -> right
    std::vector<std::string> tags = {"simfit_mu"};
    if (withComb) tags.push_back("comb");
    tags.push_back("simfit_ele");

    struct Fit { pOFit::Spec spec; TFile *fC = nullptr, *fF = nullptr; };
    std::vector<Fit> fits;
    for (const std::string &t : tags)
    {
        Fit f;
        if (!pOFit::Get(t.c_str(), f.spec)) continue;
        const TString sC = pODisc::GraphFile("charge_asym", t.c_str(), disc);
        const TString sF = pODisc::GraphFile("FBratio", t.c_str(), disc);
        if (gSystem->AccessPathName(sC) || gSystem->AccessPathName(sF)) // true = missing
        {
            std::cerr << "[WARN] observables_flav: no " << sC << " / " << sF << " -> '"
                      << f.spec.label << "' not drawn (fit not run, or its fiducial yields not"
                      << " made: analysis/run_observables.sh " << disc << ")\n";
            continue;
        }
        f.fC = TFile::Open(sC, "READ");
        f.fF = TFile::Open(sF, "READ");
        if (!f.fC || f.fC->IsZombie() || !f.fF || f.fF->IsZombie())
        {
            std::cerr << "[WARN] observables_flav: cannot open " << sC << " / " << sF << "\n";
            continue;
        }
        fits.push_back(f);
    }
    int nFlav = 0;
    for (const Fit &f : fits)
        if (f.spec.flav != "") ++nFlav;
    if (nFlav == 0)
    {
        std::cerr << "[ERROR] observables_flav: neither the mu-only nor the e-only fit is available for disc="
                  << disc << " -- fork: ./run_pO_fits.sh both flavfit --disc " << disc
                  << ", then analysis/run_observables.sh " << disc << ".\n";
        for (Fit &f : fits) { f.fC->Close(); f.fF->Close(); }
        return;
    }

    TFile *fT = TFile::Open(theoryFile, "READ");
    if (!fT || fT->IsZombie()) { std::cerr << "[WARN] no theory file -> data-only overlays\n"; fT = nullptr; }
    TGraphErrors *thWp[4], *thWm[4];
    readTheoryCharge(fT, "WPlus", thWp);
    readTheoryCharge(fT, "WMinus", thWm);

    // sum-channel theory weighted by the W+/W- abundances of the per-flavour
    // fits' FIDUCIAL yields (only the relative weights matter; the grand fit's
    // would double count them)
    std::vector<TGraphErrors *> sumTheory(4, nullptr);
    {
        TGraphErrors *gSumRef = nullptr;
        std::vector<TFile *> fY;
        for (const Fit &f : fits)
        {
            if (f.spec.flav == "") continue;
            if (!gSumRef) gSumRef = GetGraph(f.fF, "g_RFB_sum", "g_RFB_mt_sum");
            TFile *y = TFile::Open(TString::Format("../skim/rootfile/fidyields_%s_%s.root", f.spec.tag.Data(), disc), "READ");
            if (y && !y->IsZombie()) fY.push_back(y);
            else std::cerr << "[WARN] no fiducial yields of '" << f.spec.label << "' -> left out of the sum-theory weights\n";
        }
        if (fT && gSumRef && !fY.empty())
        {
            using pOAnalysis::YieldInRange;
            auto count = [&](const char *chg, int iy) -> double {
                double s = 0;
                for (TFile *y : fY)
                {
                    // h_yield_* (the name fiducial_yields.C writes); h_mt_* is the deprecated alias
                    TH1D *h = (TH1D *)y->Get(Form("h_yield_%s_y%d_FB", chg, iy));
                    if (!h) h = (TH1D *)y->Get(Form("h_mt_%s_y%d_FB", chg, iy));
                    s += YieldInRange(h, 30.0, 200.0, true).value;
                }
                return s;
            };
            sumTheory = buildSumTheorySet(count, gSumRef, thWp, thWm);
        }
        for (TFile *y : fY) { y->Close(); delete y; }
    }

    // ---- the overlaid series of one observable ------------------------------
    // shifts in units of the bin half-width: 2 fits at -+0.35, 3 at -0.6/0/+0.6,
    // with box half-widths (systBoxWidthFrac) that keep neighbours apart
    PlotStyle ps; ps.showStats = false; ps.logy = false;
    ps.systBoxWidthFrac = (fits.size() > 2) ? 0.17 : 0.20;
    ps.systBoxFillAlpha = 0.30;
    const double shift2[2] = {-0.35, 0.35}, shift3[3] = {-0.60, 0.0, 0.60};
    bool warnedStat = false;
    auto series = [&](bool isCharge, const char *name, const char *legacy) {
        std::vector<OverlaySeries> out;
        for (size_t k = 0; k < fits.size(); ++k)
        {
            TFile *f = isCharge ? fits[k].fC : fits[k].fF;
            OverlaySeries s;
            s.gTot = GetGraph(f, name, legacy);
            if (!s.gTot) { std::cerr << "[WARN] no " << name << " in " << f->GetName() << "\n"; continue; }
            s.gStat = (TGraphErrors *)f->Get(Form("%s_stat", name));
            if (!s.gStat && !warnedStat)
            {
                std::cerr << "[WARN] no " << name << "_stat in " << f->GetName()
                          << " -> that fit gets one total bar (no box); re-run the fork extraction"
                          << " so its fitted yields carry h_cov_yield[_FB]_stat\n";
                warnedStat = true;
            }
            s.label = std::string(fits[k].spec.label.Data()) + " (simfit)";
            s.color = fits[k].spec.color;
            s.marker = fits[k].spec.marker;
            s.markerSize = fits[k].spec.markerSize;
            s.xShift = (fits.size() == 1) ? 0.0 : (fits.size() == 2) ? shift2[k] : shift3[k];
            out.push_back(s);
        }
        return out;
    };

    const std::string outC = std::string("./plots/flavfit/charge_asym/") + disc;
    const std::string outF = std::string("./plots/flavfit/FBratio/") + disc;
    gSystem->mkdir(outC.c_str(), kTRUE);
    gSystem->mkdir(outF.c_str(), kTRUE);

    const TString fitTag = TString::Format("%s fit", discLabel.Data());
    GraphTuner tuneCharge = withDiscTag(makeTuneCharge(), fitTag, ps);
    GraphTuner tuneRFB = withDiscTag(makeTuneRFB(), fitTag, ps);
    // "fiducial": from r x sigma_gen, NOT the count-based yields (see the header)
    const char *sub2 = withComb ? "#mu / comb. / e, fiducial" : "#mu vs e fits, fiducial";

    SaveNiceGraph_Overlay(series(true, "g_chargeAsym", "g_chargeAsym_mt"), outC + "/chargeAsym",
                          "#eta^{l}_{CM}", "A_{ch}", "", "W #rightarrow l #nu", sub2, {}, ps, tuneCharge);
    SaveNiceGraph_Overlay(series(false, "g_RFB_sum", "g_RFB_mt_sum"), outF + "/RFB_sum",
                          "#eta^{l}_{CM}", "R_{FB}", "", "W #rightarrow l #nu", sub2, {}, ps, tuneRFB,
                          sumTheory[0], sumTheory[1], sumTheory[2], sumTheory[3]);
    SaveNiceGraph_Overlay(series(false, "g_RFB_Wp", "g_RFB_mt_Wp"), outF + "/RFB_Wp",
                          "#eta^{l}_{CM}", "R_{FB}", "", "W^{+} #rightarrow l^{+} #nu", sub2, {}, ps, tuneRFB,
                          thWp[0], thWp[1], thWp[2], thWp[3]);
    SaveNiceGraph_Overlay(series(false, "g_RFB_Wm", "g_RFB_mt_Wm"), outF + "/RFB_Wm",
                          "#eta^{l}_{CM}", "R_{FB}", "", "W^{-} #rightarrow l^{-} #bar{#nu}", sub2, {}, ps, tuneRFB,
                          thWm[0], thWm[1], thWm[2], thWm[3]);

    // ---- mu vs e compatibility (see the header) -----------------------------
    const Fit *fm = nullptr, *fe = nullptr;
    for (const Fit &f : fits)
    {
        if (f.spec.flav == "mu") fm = &f;
        if (f.spec.flav == "ele") fe = &f;
    }
    if (fm && fe)
    {
        auto compat = [&](bool isCharge, const char *name, const char *legacy, const char *what,
                          std::ofstream &out) {
            TFile *am = isCharge ? fm->fC : fm->fF, *ae = isCharge ? fe->fC : fe->fF;
            TGraphErrors *tm = GetGraph(am, name, legacy), *te = GetGraph(ae, name, legacy);
            TGraphErrors *sm = (TGraphErrors *)am->Get(Form("%s_stat", name));
            TGraphErrors *se = (TGraphErrors *)ae->Get(Form("%s_stat", name));
            if (!tm || !te || tm->GetN() != te->GetN()) return;
            const bool haveS = sm && se && sm->GetN() == tm->GetN() && se->GetN() == te->GetN();
            double c2s = 0, c2t = 0;
            int n = 0;
            for (int i = 0; i < tm->GetN(); ++i)
            {
                const double d = tm->GetPointY(i) - te->GetPointY(i);
                const double vt = std::pow(tm->GetErrorY(i), 2) + std::pow(te->GetErrorY(i), 2);
                if (vt <= 0) continue;
                c2t += d * d / vt;
                if (haveS)
                {
                    const double vs = std::pow(sm->GetErrorY(i), 2) + std::pow(se->GetErrorY(i), 2);
                    if (vs > 0) c2s += d * d / vs;
                }
                ++n;
            }
            if (n == 0) return;
            if (haveS)
                printf("[mu-vs-e] %-12s chi2/ndf = %5.1f/%d (stat only, p = %.3f)   %5.1f/%d (total, p = %.3f)\n",
                       what, c2s, n, TMath::Prob(c2s, n), c2t, n, TMath::Prob(c2t, n));
            else
                printf("[mu-vs-e] %-12s chi2/ndf = %5.1f/%d (total, p = %.3f)   [no stat twins]\n",
                       what, c2t, n, TMath::Prob(c2t, n));
            out << what << "," << n << "," << (haveS ? c2s : -1.0) << ","
                << (haveS ? TMath::Prob(c2s, n) : -1.0) << "," << c2t << "," << TMath::Prob(c2t, n) << "\n";
        };
        std::ofstream cC((outC + "/mu_vs_e_chi2.csv").c_str()), cF((outF + "/mu_vs_e_chi2.csv").c_str());
        cC << "observable,ndf,chi2_stat,p_stat,chi2_total,p_total\n";
        cF << "observable,ndf,chi2_stat,p_stat,chi2_total,p_total\n";
        compat(true, "g_chargeAsym", "g_chargeAsym_mt", "A_ch", cC);
        compat(false, "g_RFB_sum", "g_RFB_mt_sum", "R_FB (sum)", cF);
        compat(false, "g_RFB_Wp", "g_RFB_mt_Wp", "R_FB (W+)", cF);
        compat(false, "g_RFB_Wm", "g_RFB_mt_Wm", "R_FB (W-)", cF);
    }
    else
        std::cout << "[mu-vs-e] only one flavour available -> no compatibility numbers\n";

    for (Fit &f : fits) { f.fC->Close(); f.fF->Close(); delete f.fC; delete f.fF; }
    if (fT) { fT->Close(); delete fT; }
    std::cout << "[OK] Saved mu-vs-e overlays (disc=" << disc << ") to: " << outC << " and " << outF << "\n";
}

// =============================================================================
// observables(disc) -- both views of one discriminant: the grand fit
// (observables_comb) and the mu-vs-e overlay (observables_flav).
//   root -l -b -q 'observables.C+("leppt_mt40")'
// =============================================================================
void observables(const char *disc = "leppt_mt40")
{
    observables_comb(disc);
    observables_flav(disc);
}
