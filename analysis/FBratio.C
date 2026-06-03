#include "TFile.h"
#include "TH1D.h"
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
    bool useMT = true,                                     // true: h_mt_*, false: h_met_*
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
    // bin contents already produced by the skim macros (h_mt_W{p,m}_y0..y11).
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

    auto get_yield = [&](int iy, bool wantWp) -> Yield
    {
        TString name = useMT
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

            if (sumCharges)
            {
                // F = W+ + W-, B = W+ + W-
                FB = get_yield(iyF, true) + get_yield(iyF, false);
                BB = get_yield(iyB, true) + get_yield(iyB, false);
            }
            else
            {
                // choose W+ or W-
                if (useWp)
                {
                    FB = get_yield(iyF, true);
                    BB = get_yield(iyB, true);
                }
                if (useWm)
                {
                    FB = get_yield(iyF, false);
                    BB = get_yield(iyB, false);
                }
            }

            if (FB.value > 0.0 && BB.value > 0.0)
            {
                yv[iabs] = FB.value / BB.value;
                ey[iabs] = RatioErr(FB, BB);
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
        TString gname = useMT ? "g_RFB_mt_sum" : "g_RFB_met_sum";
        TString gtitle = useMT
                             ? "R_{FB} (sum charges) from m_{T} yields; |y| bin; R_{FB}"
                             : "R_{FB} (sum charges) from MET yields; |y| bin; R_{FB}";
        cout << "[INFO] " << " Now producing Sum " << endl;
        graphs.push_back(build_graph(gname, gtitle, false, false, true));

        if (alsoWriteChargeSeparated)
        {
            TString gnameP = useMT ? "g_RFB_mt_Wp" : "g_RFB_met_Wp";
            TString gtitleP = useMT
                                  ? "R_{FB} (W^{+}) from m_{T} yields; |y| bin; R_{FB}"
                                  : "R_{FB} (W^{+}) from MET yields; |y| bin; R_{FB}";
            cout << "[INFO] " << " Now producing + " << endl;
            graphs.push_back(build_graph(gnameP, gtitleP, true, false, false));

            TString gnameM = useMT ? "g_RFB_mt_Wm" : "g_RFB_met_Wm";
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
        TString gnameP = useMT ? "g_RFB_mt_Wp" : "g_RFB_met_Wp";
        TString gtitleP = useMT
                              ? "R_{FB} (W^{+}) from m_{T} yields; |y| bin; R_{FB}"
                              : "R_{FB} (W^{+}) from MET yields; |y| bin; R_{FB}";
        graphs.push_back(build_graph(gnameP, gtitleP, true, false, false));

        TString gnameM = useMT ? "g_RFB_mt_Wm" : "g_RFB_met_Wm";
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
        g->Write("", TObject::kOverwrite);

    fout->Close();
    delete fout;
    f->Close();
    delete f;

    std::cout << "[OK] Wrote " << graphs.size() << " R_FB graph(s) into " << outFile << "\n";
    for (auto *g : graphs)
        delete g;
}