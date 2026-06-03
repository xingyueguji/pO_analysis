#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TString.h"
#include <iostream>
#include <vector>
#include <cmath>

#include "analysis_helpers.h" // YieldInRange, AsymErr, kPORapidityShift

// ---------- main ----------
void charge_asym(
    const char *inFile = "../skim/rootfile/WToMuNu_pO_PFMet_hist.root",
    const char *outFile = "../skim/rootfile/charge_asym.root", // can overwrite/update same
    bool useMT = true,                                         // true -> use h_mt_*, false -> use h_met_*
    bool integrateFull = true,                                 // true -> integrate full x-range (excluding under/overflow)
    double xMin = 30.0,                                        // if integrateFull=false, integrate [xMin, xMax]
    double xMax = 200.0,
    int NY = 12)
{
    using pOAnalysis::AsymErr;
    using pOAnalysis::YieldInRange;

    TFile *f = TFile::Open(inFile, "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open input file: " << inFile << "\n";
        return;
    }

    // yEdges: rapidity bin edges in the *lab frame*, used only to label the
    // x-axis of the output TGraph. The histogram integration uses bin
    // contents already produced by the skim macros (h_mt_W{p,m}_y0..y11).
    //
    // These edges are symmetric around 0 in the *lab* frame (unlike FBratio.C
    // which uses edges symmetric around deltaY for F/B pairing). Charge
    // asymmetry doesn't need a symmetric CM-frame binning, so this is fine,
    // but be aware the two analysis macros label their y-axis differently.
    std::vector<double> yEdges; // empty => use index axis
    yEdges = {
        -2.4, -2.0, -1.6, -1.2, -0.8, -0.4,
        0.0, 0.4, 0.8, 1.2, 1.6, 2.0,
        2.4};

    // pO rapidity shift (lab -> CM frame). Defined once in analysis_helpers.h.
    const double deltaY = pOAnalysis::kPORapidityShift;
    std::vector<double> yEdgesCM;
    yEdgesCM.reserve(yEdges.size());

    for (double y : yEdges)
    {
        yEdgesCM.push_back(y - deltaY);
    }

    std::vector<double> x(NY), ex(NY), y(NY), ey(NY);

    for (int iy = 0; iy < NY; ++iy)
    {
        TString hWpName = useMT ? Form("h_mt_Wp_y%d", iy) : Form("h_met_Wp_y%d", iy);
        TString hWmName = useMT ? Form("h_mt_Wm_y%d", iy) : Form("h_met_Wm_y%d", iy);

        TH1D *hWp = (TH1D *)f->Get(hWpName);
        TH1D *hWm = (TH1D *)f->Get(hWmName);

        if (!hWp || !hWm)
        {
            std::cerr << "[WARN] Missing hist(s) for iy=" << iy
                      << " : " << hWpName << " or " << hWmName << "\n";
            x[iy] = iy + 0.5;
            ex[iy] = 0.5;
            y[iy] = 0.0;
            ey[iy] = 0.0;
            continue;
        }

        const auto Np = YieldInRange(hWp, xMin, xMax, integrateFull);
        const auto Nm = YieldInRange(hWm, xMin, xMax, integrateFull);

        const double S = Np.value + Nm.value;
        double A = 0.0;
        double sA = 0.0;

        if (S > 0.0)
        {
            A = (Np.value - Nm.value) / S;
            sA = AsymErr(Np, Nm);
        }

        // x-axis: y-bin center
        if (!yEdgesCM.empty() && (int)yEdgesCM.size() == NY + 1)
        {
            x[iy] = 0.5 * (yEdgesCM[iy] + yEdgesCM[iy + 1]);
            ex[iy] = 0.5 * (yEdgesCM[iy + 1] - yEdgesCM[iy]);
        }
        else
        {
            x[iy] = iy + 0.5;
            ex[iy] = 0.5;
        }

        y[iy] = A;
        ey[iy] = sA;

        std::cout << "[INFO] iy=" << iy
                  << "  Np=" << Np.value << " Nm=" << Nm.value
                  << "  A=" << A << " +/- " << sA << "\n";
    }

    // Build graph
    TString gname = useMT ? "g_chargeAsym_mt" : "g_chargeAsym_met";
    TString gtitle = useMT
                         ? "Charge asymmetry from m_{T} yields; y-bin; (N^{+}-N^{-})/(N^{+}+N^{-})"
                         : "Charge asymmetry from MET yields; y-bin; (N^{+}-N^{-})/(N^{+}+N^{-})";

    TGraphErrors *g = new TGraphErrors(NY, x.data(), y.data(), ex.data(), ey.data());
    g->SetName(gname);
    g->SetTitle(gtitle);

    // Write out
    TFile *fout = TFile::Open(outFile, "UPDATE");
    if (!fout || fout->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open output file for UPDATE: " << outFile << "\n";
        delete g;
        f->Close();
        delete f;
        return;
    }

    fout->cd();
    g->Write("", TObject::kOverwrite);

    fout->Close();
    delete fout;
    f->Close();
    delete f;

    std::cout << "[OK] Wrote " << gname << " into " << outFile << "\n";
}