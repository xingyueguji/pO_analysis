#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
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
    bool useMT = true,                    // raw skim files only: h_mt_* vs h_met_*
                                          // (ignored when h_yield_* is present)
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
    // contents: either the raw skim histos (h_mt_/h_met_W{p,m}_y0..y11) or,
    // in the production path, the fork's fitted-yield containers h_yield_*.
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

    // Fitted-yield covariance from the simfit grand fit (2026-08-04): the
    // 24x24 matrix h_cov_yield, fixed order [Wp_y0..11, Wm_y0..11], written by
    // the fork's extract_pO_simfit.C into comb_fitted_yields.root. Absent in
    // raw-skim and legacy per-flavour-fit files -> cov = 0 reproduces the old
    // independent-yield errors exactly.
    TH2D *hcov = (TH2D *)f->Get("h_cov_yield");
    if (hcov && hcov->GetNbinsX() != 2 * NY)
    {
        std::cerr << "[WARN] h_cov_yield has " << hcov->GetNbinsX()
                  << " rows, expected " << 2 * NY << " -> covariance ignored\n";
        hcov = nullptr;
    }
    if (hcov)
        std::cout << "[INFO] simfit covariance found -> A errors include the "
                     "cov(N+, N-) cross term\n";

    // STATISTICAL component (2026-09-15): the same matrix with the constrained
    // nuisances conditioned out (h_cov_yield_stat, written by the fork's
    // extract_pO_simfit.C::ComputeStatCov -- the Gaussian-exact equivalent of
    // refitting with every nuisance frozen at its post-fit value). Present it
    // and this macro writes a SECOND graph g_chargeAsym_stat carrying the
    // statistical error only; observables.C then draws the bars from that one
    // and the systematic -- the quadratic difference -- as a box per point.
    // Absent (raw skim, legacy per-flavour fits, any pre-2026-09-14 extraction)
    // -> only the total-error graph is written, exactly as before.
    //
    // NB the errors stored in h_yield_* are the TOTAL ones (= sqrt of the
    // h_cov_yield diagonal), so the statistical A must take BOTH its diagonal
    // terms and its cross term from the stat matrix -- not the histograms.
    TH2D *hcovS = (TH2D *)f->Get("h_cov_yield_stat");
    if (hcovS && hcovS->GetNbinsX() != 2 * NY)
    {
        std::cerr << "[WARN] h_cov_yield_stat has " << hcovS->GetNbinsX()
                  << " rows, expected " << 2 * NY << " -> stat component ignored\n";
        hcovS = nullptr;
    }
    if (hcovS)
        std::cout << "[INFO] stat-only covariance found -> also writing "
                     "g_chargeAsym_stat (stat. error; syst = sqrt(tot^2 - stat^2))\n";

    std::vector<double> x(NY), ex(NY), y(NY), ey(NY), eyStat(NY, 0.0);

    for (int iy = 0; iy < NY; ++iy)
    {
        // The Combine fork's <chan>_fitted_yields.root stores 1-bin FITTED
        // yields under the discriminant-neutral name h_yield_* (the fit uses the
        // PF MET shape, so calling them h_mt_* was misleading -- that legacy
        // alias is still written and still accepted below). A RAW skim file
        // instead holds the real distributions, where h_mt_*/h_met_* do mean
        // m_T / MET and `useMT` genuinely chooses between them.
        TString hWpName = Form("h_yield_Wp_y%d", iy);
        TString hWmName = Form("h_yield_Wm_y%d", iy);
        if (!f->Get(hWpName) || !f->Get(hWmName))
        {
            hWpName = useMT ? Form("h_mt_Wp_y%d", iy) : Form("h_met_Wp_y%d", iy);
            hWmName = useMT ? Form("h_mt_Wm_y%d", iy) : Form("h_met_Wm_y%d", iy);
        }

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
        double sAstat = 0.0;

        if (S > 0.0)
        {
            A = (Np.value - Nm.value) / S;
            // cov(Wp_yi, Wm_yi): matrix order is [Wp_y0..11, Wm_y0..11]
            const double cov = hcov ? hcov->GetBinContent(iy + 1, NY + iy + 1) : 0.0;
            sA = AsymErr(Np, Nm, cov);

            // statistical A: same values, all three second moments from the
            // conditioned matrix. Note the 3% luminosity nuisance is fully
            // correlated across every channel, so it rescales N+ and N-
            // coherently and cancels in A -- expect sAstat ~ sA here (it is
            // the cross section, not the asymmetry, that lumi dominates).
            if (hcovS)
            {
                const pOAnalysis::Yield NpS(Np.value, std::sqrt(std::max(0.0, hcovS->GetBinContent(iy + 1, iy + 1))));
                const pOAnalysis::Yield NmS(Nm.value, std::sqrt(std::max(0.0, hcovS->GetBinContent(NY + iy + 1, NY + iy + 1))));
                sAstat = AsymErr(NpS, NmS, hcovS->GetBinContent(iy + 1, NY + iy + 1));
            }
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
        eyStat[iy] = sAstat;

        std::cout << "[INFO] iy=" << iy
                  << "  Np=" << Np.value << " Nm=" << Nm.value
                  << "  A=" << A << " +/- " << sA;
        if (hcovS)
            std::cout << " (total)  = +/- " << sAstat << " (stat) +/- "
                      << std::sqrt(std::max(0.0, sA * sA - sAstat * sAstat)) << " (syst)";
        std::cout << "\n";
    }

    // Build graph
    // Honest primary name + the legacy "_mt"/"_met" alias (observables.C
    // prefers the primary and falls back to the alias for old files).
    TString gname      = "g_chargeAsym";
    TString gnameLegacy = useMT ? "g_chargeAsym_mt" : "g_chargeAsym_met";
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
    g->Write(gnameLegacy, TObject::kOverwrite); // deprecated alias

    // the statistical-error twin (same points, stat-only bars); observables.C
    // picks it up by name and turns the difference into the systematic boxes
    if (hcovS)
    {
        TGraphErrors *gs = new TGraphErrors(NY, x.data(), y.data(), ex.data(), eyStat.data());
        gs->SetName(gname + "_stat");
        gs->SetTitle(gtitle);
        gs->Write("", TObject::kOverwrite);
        delete gs;
    }

    fout->Close();
    delete fout;
    f->Close();
    delete f;

    std::cout << "[OK] Wrote " << gname << " into " << outFile << "\n";
}