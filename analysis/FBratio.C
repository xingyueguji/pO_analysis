#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TString.h"
#include <iostream>
#include <vector>
#include <cmath>

#include "analysis_helpers.h" // YieldInRange, RatioErr, kPORapidityShift

// ---------- main ----------
// Production input (2026-09-22): the FIDUCIAL yields r x sigma_gen of a fit,
// ../skim/rootfile/fidyields_<fit>_<disc>.root (analysis/fiducial_yields.C,
// same histogram names) -- every observable is built from r x sigma_gen,
// never from the raw fitted yields (user directive: there is no dedicated
// efficiency/acceptance correction -- and here it matters most: F and B are
// different |eta_lab| regions, so a raw-count ratio carries the detector
// acceptance, up to 60% apart between mu and e). Any file with h_yield_*
// (+ h_cov_yield_FB*) still works.
void FBratio(
    const char *inFile = "../skim/rootfile/WToMuNu_pO_PFMet_hist.root",
    const char *outFile = "../skim/rootfile/FBratio.root", // update same
    bool useMT = true,                // raw skim files only (ignored if h_yield_* exists)
    bool integrateFull = true,                             // integrate full range or [xMin,xMax]
    double xMin = 30.0,
    double xMax = 200.0,
    int NY = 12,
    bool combineCharges = true,          // if true -> use (Wp+Wm) yields
    bool alsoWriteChargeSeparated = true // if combineCharges==true, optionally also write W+ and W-
)
{
    using pOAnalysis::RatioErr;
    using pOAnalysis::Yield;
    using pOAnalysis::YieldInRange;

    if (NY % 2 != 0)
    {
        std::cerr << "[ERROR] NY must be even to pair +/-y bins cleanly. NY=" << NY << "\n";
        return;
    }

    TFile *f = TFile::Open(inFile, "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open input file: " << inFile << "\n";
        return;
    }

    // yEdges: rapidity bin edges in the *lab frame*, used only to label the
    // x-axis of the output TGraph. The histogram integration itself uses the
    // bin contents: raw skim histos, or the fork's h_yield_* fitted yields.
    //
    // These 13 edges are chosen symmetric around deltaY = 0.3466 (the pO
    // rapidity shift), so that after the lab->CM shift below they become
    // symmetric around 0 in the CM frame. F/B pairing (iyB <-> NY-1-iyB)
    // then matches CM-frame |y| bins cleanly.
    std::vector<double> yEdges; // empty -> use index axis
    yEdges = {
        -1.7068,
        -1.3646,
        -1.0223,
        -0.6801,
        -0.3379,
        0.0044,
        0.3466,
        0.6888,
        1.0311,
        1.3733,
        1.7155,
        2.0578,
        2.4000};

    // pO rapidity shift (lab -> CM frame). Defined once in analysis_helpers.h.
    const double deltaY = pOAnalysis::kPORapidityShift;
    std::vector<double> yEdgesCM;
    yEdgesCM.reserve(yEdges.size());

    for (double y : yEdges)
    {
        yEdgesCM.push_back(y - deltaY);
    }

    const int Nabs = NY / 2;
    std::vector<double> x(Nabs), ex(Nabs);

    for (int iabs = 0; iabs < Nabs; ++iabs)
    {
        // backward bin index on negative side:
        int iyB = iabs;          // assumes bins ordered from negative to positive
        int iyF = NY - 1 - iabs; // paired positive bin

        if (!yEdgesCM.empty() && (int)yEdgesCM.size() == NY + 1)
        {
            // Use |y| bin center from the positive side bin edges
            double y1 = yEdgesCM[iyF];
            double y2 = yEdgesCM[iyF + 1];
            x[iabs] = 0.5 * (std::fabs(y1) + std::fabs(y2));
            ex[iabs] = 0.5 * (std::fabs(y2 - y1));
        }
        else
        {
            x[iabs] = iabs + 0.5; // |y| bin index axis
            ex[iabs] = 0.5;
        }
    }

    // Fitted-yield covariance from the simfit grand fit (2026-08-04): the
    // 24x24 matrix h_cov_yield_FB, fixed order [Wp_y0..11, Wm_y0..11], written
    // by the fork's extract_pO_simfit.C into comb_fitted_yields.root. Absent in
    // raw-skim and legacy per-flavour-fit files -> all covariances 0, which
    // reproduces the old independent-yield errors exactly.
    TH2D *hcov = (TH2D *)f->Get("h_cov_yield_FB");
    if (hcov && hcov->GetNbinsX() != 2 * NY)
    {
        std::cerr << "[WARN] h_cov_yield_FB has " << hcov->GetNbinsX()
                  << " rows, expected " << 2 * NY << " -> covariance ignored\n";
        hcov = nullptr;
    }
    if (hcov)
        std::cout << "[INFO] simfit covariance found -> R_FB errors include all "
                     "cross terms (within F, within B, and cov(F, B))\n";

    // STATISTICAL component (2026-09-15): the same matrix with the constrained
    // nuisances conditioned out (h_cov_yield_FB_stat, fork
    // extract_pO_simfit.C::ComputeStatCov). Present it and every graph below is
    // built TWICE -- once with the total matrix (the primary g_RFB_*) and once
    // with the stat one (g_RFB_*_stat) -- so observables.C can draw statistical
    // bars with the systematic, sqrt(tot^2 - stat^2), as a box per point.
    // Absent -> only the total graphs are written, exactly as before.
    TH2D *hcovS = (TH2D *)f->Get("h_cov_yield_FB_stat");
    if (hcovS && hcovS->GetNbinsX() != 2 * NY)
    {
        std::cerr << "[WARN] h_cov_yield_FB_stat has " << hcovS->GetNbinsX()
                  << " rows, expected " << 2 * NY << " -> stat component ignored\n";
        hcovS = nullptr;
    }
    if (hcovS)
        std::cout << "[INFO] stat-only covariance found -> also writing g_RFB_*_stat\n";

    // covariance-matrix index of one yield: Wp_yi -> iy, Wm_yi -> NY + iy
    auto covIdx = [&](int iy, bool wantWp) { return (wantWp ? 0 : NY) + iy; };
    auto covEl = [&](const TH2D *M, int a, int b) -> double {
        return M ? M->GetBinContent(a + 1, b + 1) : 0.0;
    };

    auto get_yield = [&](int iy, bool wantWp) -> Yield
    {
        // Prefer the fork's discriminant-neutral fitted-yield name; fall back
        // to the raw-skim h_mt_/h_met_ pair (see charge_asym.C for the why).
        TString name = Form("h_yield_%s_y%d_FB", wantWp ? "Wp" : "Wm", iy);
        if (!f->Get(name))
            name = useMT
                       ? Form("h_mt_%s_y%d_FB", wantWp ? "Wp" : "Wm", iy)
                       : Form("h_met_%s_y%d_FB", wantWp ? "Wp" : "Wm", iy);

        TH1D *h = (TH1D *)f->Get(name);
        if (!h)
        {
            std::cerr << "[WARN] Missing " << name << "\n";
            return Yield();
        }
        return YieldInRange(h, xMin, xMax, integrateFull);
    };

    // Make graphs (combined and/or separated)
    // M = the yield covariance the errors are taken from: hcov for the primary
    // (TOTAL) graphs, hcovS for the statistical twins. Null -> the legacy
    // independent-yield errors straight out of the histograms.
    auto build_graph = [&](const char *gname, const char *gtitle,
                           bool useWp, bool useWm, bool sumCharges,
                           const TH2D *M) -> TGraphErrors *
    {
        std::vector<double> yv(Nabs, 0.0), ey(Nabs, 0.0);

        for (int iabs = 0; iabs < Nabs; ++iabs)
        {
            int iyB = iabs;
            int iyF = NY - 1 - iabs;

            Yield FB, BB; // Forward yield, Backward yield (value + error)
            std::vector<int> idxF, idxB; // covariance-matrix indices of F / B parts

            if (sumCharges)
            {
                // F = W+ + W-, B = W+ + W-
                FB = get_yield(iyF, true) + get_yield(iyF, false);
                BB = get_yield(iyB, true) + get_yield(iyB, false);
                idxF = {covIdx(iyF, true), covIdx(iyF, false)};
                idxB = {covIdx(iyB, true), covIdx(iyB, false)};
            }
            else
            {
                // choose W+ or W-
                if (useWp)
                {
                    FB = get_yield(iyF, true);
                    BB = get_yield(iyB, true);
                    idxF = {covIdx(iyF, true)};
                    idxB = {covIdx(iyB, true)};
                }
                if (useWm)
                {
                    FB = get_yield(iyF, false);
                    BB = get_yield(iyB, false);
                    idxF = {covIdx(iyF, false)};
                    idxB = {covIdx(iyB, false)};
                }
            }

            // With the simfit covariance: replace the independent-quadrature
            // errors with the full Var(F) = sum_ab C_ab over the F parts (this
            // adds the 2*cov(Wp,Wm) term the Yield operator+ cannot know), and
            // build cov(F, B) for the ratio. Without hcov everything is 0/kept.
            double covFB = 0.0;
            if (M)
            {
                auto setVar = [&](Yield &Y, const std::vector<int> &idx) {
                    double var = 0.0;
                    for (int a : idx)
                        for (int b : idx)
                            var += covEl(M, a, b);
                    if (var > 0.0)
                        Y.error = std::sqrt(var);
                };
                setVar(FB, idxF);
                setVar(BB, idxB);
                for (int a : idxF)
                    for (int b : idxB)
                        covFB += covEl(M, a, b);
            }

            if (FB.value > 0.0 && BB.value > 0.0)
            {
                yv[iabs] = FB.value / BB.value;
                ey[iabs] = RatioErr(FB, BB, covFB);
            }
            else
            {
                yv[iabs] = 0.0;
                ey[iabs] = 0.0;
            }

            if (M == hcov) // the stat pass repeats the same values: report once
                std::cout << "[INFO] |y|bin=" << iabs
                          << "  F=" << FB.value << "  B=" << BB.value
                          << "  R_FB=" << yv[iabs] << " +/- " << ey[iabs] << "\n";
        }

        TGraphErrors *g = new TGraphErrors(Nabs, x.data(), yv.data(), ex.data(), ey.data());
        g->SetName(gname);
        g->SetTitle(gtitle);
        return g;
    };

    // Build requested graphs
    std::vector<TGraphErrors *> graphs;     // primary, TOTAL error
    std::vector<TGraphErrors *> graphsStat; // statistical twins, "<name>_stat"

    // one call -> the total graph plus (when the stat matrix exists) its twin
    auto addGraph = [&](const TString &gname, const TString &gtitle,
                        bool useWp, bool useWm, bool sumCharges) {
        graphs.push_back(build_graph(gname.Data(), gtitle.Data(), useWp, useWm, sumCharges, hcov));
        if (hcovS)
            graphsStat.push_back(build_graph((gname + "_stat").Data(), gtitle.Data(),
                                             useWp, useWm, sumCharges, hcovS));
    };

    if (combineCharges)
    {
        TString gname = "g_RFB_sum";   // legacy alias: g_RFB_mt_sum / g_RFB_met_sum (written too)
        TString gtitle = useMT
                             ? "R_{FB} (sum charges) from m_{T} yields; |y| bin; R_{FB}"
                             : "R_{FB} (sum charges) from MET yields; |y| bin; R_{FB}";
        cout << "[INFO] " << " Now producing Sum " << endl;
        addGraph(gname, gtitle, false, false, true);

        if (alsoWriteChargeSeparated)
        {
            TString gnameP = "g_RFB_Wp";   // legacy alias: g_RFB_mt_Wp / g_RFB_met_Wp (written too)
            TString gtitleP = useMT
                                  ? "R_{FB} (W^{+}) from m_{T} yields; |y| bin; R_{FB}"
                                  : "R_{FB} (W^{+}) from MET yields; |y| bin; R_{FB}";
            cout << "[INFO] " << " Now producing + " << endl;
            addGraph(gnameP, gtitleP, true, false, false);

            TString gnameM = "g_RFB_Wm";   // legacy alias: g_RFB_mt_Wm / g_RFB_met_Wm (written too)
            TString gtitleM = useMT
                                  ? "R_{FB} (W^{-}) from m_{T} yields; |y| bin; R_{FB}"
                                  : "R_{FB} (W^{-}) from MET yields; |y| bin; R_{FB}";
            cout << "[INFO] " << " Now producing - " << endl;
            addGraph(gnameM, gtitleM, false, true, false);
        }
    }
    else
    {
        // Only charge-separated
        TString gnameP = "g_RFB_Wp";   // legacy alias: g_RFB_mt_Wp / g_RFB_met_Wp (written too)
        TString gtitleP = useMT
                              ? "R_{FB} (W^{+}) from m_{T} yields; |y| bin; R_{FB}"
                              : "R_{FB} (W^{+}) from MET yields; |y| bin; R_{FB}";
        addGraph(gnameP, gtitleP, true, false, false);

        TString gnameM = "g_RFB_Wm";   // legacy alias: g_RFB_mt_Wm / g_RFB_met_Wm (written too)
        TString gtitleM = useMT
                              ? "R_{FB} (W^{-}) from m_{T} yields; |y| bin; R_{FB}"
                              : "R_{FB} (W^{-}) from MET yields; |y| bin; R_{FB}";
        addGraph(gnameM, gtitleM, false, true, false);
    }

    // Write out
    TFile *fout = TFile::Open(outFile, "UPDATE");
    if (!fout || fout->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open output file for UPDATE: " << outFile << "\n";
        for (auto *g : graphs)
            delete g;
        f->Close();
        delete f;
        return;
    }

    fout->cd();
    for (auto *g : graphs)
    {
        g->Write("", TObject::kOverwrite);
        // Deprecated alias under the old discriminant-flavoured name, so any
        // existing reader of g_RFB_{mt,met}_* keeps working. The primary name
        // (g_RFB_*) is discriminant-neutral because in the production path the
        // yields come from the PF-MET-shape fit, NOT from m_T.
        TString legacy = TString(g->GetName());
        legacy.ReplaceAll("g_RFB_", useMT ? "g_RFB_mt_" : "g_RFB_met_");
        g->Write(legacy, TObject::kOverwrite);
    }
    // statistical twins: primary name only (they are new, nothing reads an alias)
    for (auto *g : graphsStat)
        g->Write("", TObject::kOverwrite);

    // stat / syst breakdown, so the console record carries the split too
    for (size_t k = 0; k < graphsStat.size() && k < graphs.size(); ++k)
    {
        const TGraphErrors *gt = graphs[k], *gs = graphsStat[k];
        for (int i = 0; i < gt->GetN() && i < gs->GetN(); ++i)
        {
            const double et = gt->GetErrorY(i), es = gs->GetErrorY(i);
            printf("[INFO] %-12s |y|bin=%d  R_FB = %.5f +/- %.5f (stat) +/- %.5f (syst)  [total %.5f]\n",
                   gt->GetName(), i, gt->GetPointY(i), es,
                   std::sqrt(std::max(0.0, et * et - es * es)), et);
        }
    }

    fout->Close();
    delete fout;
    f->Close();
    delete f;

    std::cout << "[OK] Wrote " << graphs.size() << " R_FB graph(s)"
              << (graphsStat.empty() ? "" : Form(" + %d stat twin(s)", (int)graphsStat.size()))
              << " into " << outFile << "\n";
    for (auto *g : graphs)
        delete g;
    for (auto *g : graphsStat)
        delete g;
}