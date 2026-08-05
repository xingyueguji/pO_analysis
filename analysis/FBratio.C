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
    // covariance-matrix index of one yield: Wp_yi -> iy, Wm_yi -> NY + iy
    auto covIdx = [&](int iy, bool wantWp) { return (wantWp ? 0 : NY) + iy; };
    auto covEl = [&](int a, int b) -> double {
        return hcov ? hcov->GetBinContent(a + 1, b + 1) : 0.0;
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
    auto build_graph = [&](const char *gname, const char *gtitle,
                           bool useWp, bool useWm, bool sumCharges) -> TGraphErrors *
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
            if (hcov)
            {
                auto setVar = [&](Yield &Y, const std::vector<int> &idx) {
                    double var = 0.0;
                    for (int a : idx)
                        for (int b : idx)
                            var += covEl(a, b);
                    if (var > 0.0)
                        Y.error = std::sqrt(var);
                };
                setVar(FB, idxF);
                setVar(BB, idxB);
                for (int a : idxF)
                    for (int b : idxB)
                        covFB += covEl(a, b);
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
    std::vector<TGraphErrors *> graphs;

    if (combineCharges)
    {
        TString gname = "g_RFB_sum";   // legacy alias: g_RFB_mt_sum / g_RFB_met_sum (written too)
        TString gtitle = useMT
                             ? "R_{FB} (sum charges) from m_{T} yields; |y| bin; R_{FB}"
                             : "R_{FB} (sum charges) from MET yields; |y| bin; R_{FB}";
        cout << "[INFO] " << " Now producing Sum " << endl;
        graphs.push_back(build_graph(gname, gtitle, false, false, true));

        if (alsoWriteChargeSeparated)
        {
            TString gnameP = "g_RFB_Wp";   // legacy alias: g_RFB_mt_Wp / g_RFB_met_Wp (written too)
            TString gtitleP = useMT
                                  ? "R_{FB} (W^{+}) from m_{T} yields; |y| bin; R_{FB}"
                                  : "R_{FB} (W^{+}) from MET yields; |y| bin; R_{FB}";
            cout << "[INFO] " << " Now producing + " << endl;
            graphs.push_back(build_graph(gnameP, gtitleP, true, false, false));

            TString gnameM = "g_RFB_Wm";   // legacy alias: g_RFB_mt_Wm / g_RFB_met_Wm (written too)
            TString gtitleM = useMT
                                  ? "R_{FB} (W^{-}) from m_{T} yields; |y| bin; R_{FB}"
                                  : "R_{FB} (W^{-}) from MET yields; |y| bin; R_{FB}";
            cout << "[INFO] " << " Now producing - " << endl;
            graphs.push_back(build_graph(gnameM, gtitleM, false, true, false));
        }
    }
    else
    {
        // Only charge-separated
        TString gnameP = "g_RFB_Wp";   // legacy alias: g_RFB_mt_Wp / g_RFB_met_Wp (written too)
        TString gtitleP = useMT
                              ? "R_{FB} (W^{+}) from m_{T} yields; |y| bin; R_{FB}"
                              : "R_{FB} (W^{+}) from MET yields; |y| bin; R_{FB}";
        graphs.push_back(build_graph(gnameP, gtitleP, true, false, false));

        TString gnameM = "g_RFB_Wm";   // legacy alias: g_RFB_mt_Wm / g_RFB_met_Wm (written too)
        TString gtitleM = useMT
                              ? "R_{FB} (W^{-}) from m_{T} yields; |y| bin; R_{FB}"
                              : "R_{FB} (W^{-}) from MET yields; |y| bin; R_{FB}";
        graphs.push_back(build_graph(gnameM, gtitleM, false, true, false));
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

    fout->Close();
    delete fout;
    f->Close();
    delete f;

    std::cout << "[OK] Wrote " << graphs.size() << " R_FB graph(s) into " << outFile << "\n";
    for (auto *g : graphs)
        delete g;
}