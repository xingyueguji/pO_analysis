#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TString.h"
#include <iostream>
#include <vector>
#include <cmath>

// ---------- helpers ----------
static double YieldInRange(const TH1D *h, double xmin, double xmax)
{
    if (!h)
        return 0.0;

    // If xmin/xmax are NaN -> integrate full histogram (excluding under/overflow by default).
    if (std::isnan(xmin) || std::isnan(xmax))
        return h->Integral(1, h->GetNbinsX());

    int b1 = h->GetXaxis()->FindBin(xmin);
    int b2 = h->GetXaxis()->FindBin(xmax);

    // protect edges
    if (b1 < 1)
        b1 = 1;
    if (b2 > h->GetNbinsX())
        b2 = h->GetNbinsX();

    return h->Integral(b1, b2);
}

static double AsymErr(double Np, double Nm)
{
    const double S = Np + Nm;
    if (S <= 0.0)
        return 0.0;
    // sigma_A^2 = 4 Np Nm / S^3
    return std::sqrt(4.0 * Np * Nm / (S * S * S));
}

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
    TFile *f = TFile::Open(inFile, "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open input file: " << inFile << "\n";
        return;
    }

    // If you have real y bin edges, put them here.
    // Otherwise we use bin centers = iy + 0.5 with x-error = 0.5.
    std::vector<double> yEdges; // empty => use index axis
                                // Example if you later know them:
    yEdges = {
        -2.4, -2.0, -1.6, -1.2, -0.8, -0.4,
        0.0, 0.4, 0.8, 1.2, 1.6, 2.0,
        2.4};

    const double deltaY = 0.3466; // pO rapidity shift
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

        double Np = 0.0, Nm = 0.0;
        if (integrateFull)
        {
            Np = YieldInRange(hWp, NAN, NAN);
            Nm = YieldInRange(hWm, NAN, NAN);
        }
        else
        {
            Np = YieldInRange(hWp, xMin, xMax);
            Nm = YieldInRange(hWm, xMin, xMax);
        }

        const double S = Np + Nm;
        double A = 0.0;
        double sA = 0.0;

        if (S > 0.0)
        {
            A = (Np - Nm) / S;
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
                  << "  Np=" << Np << " Nm=" << Nm
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