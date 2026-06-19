#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TLine.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "plotting_helper.C"               // where your PlotStyle + SaveNicePlotGraph live
#include "../analysis/analysis_helpers.h"  // pOAnalysis::YieldInRange (Sumw2-aware yields)

void observables()
{

    gStyle->SetEndErrorSize(4);

    const std::string inFileCharge = "../skim/rootfile/charge_asym.root";
    TFile *fCharge = TFile::Open(inFileCharge.c_str(), "READ");
    if (!fCharge || fCharge->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open file: " << inFileCharge << "\n";
        return;
    }

    const std::string inFileFB = "../skim/rootfile/FBratio.root";
    TFile *fFB = TFile::Open(inFileFB.c_str(), "READ");
    if (!fFB || fFB->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open file: " << inFileFB << "\n";
        return;
    }

    const std::string inFileFB_theory = "./RpO_rootfile/RpO_FB_graphs.root";
    TFile *fFB_theory = TFile::Open(inFileFB_theory.c_str(), "READ");
    if (!fFB_theory || fFB_theory->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open file: " << inFileFB_theory << "\n";
        return;
    }

    // output directories
    const std::string outBase = "./plots";
    const std::string outFBDir = outBase + "/FBratio";
    const std::string outChargeDir = outBase + "/charge_asym";
    gSystem->mkdir(outFBDir.c_str(), kTRUE);
    gSystem->mkdir(outChargeDir.c_str(), kTRUE);

    PlotStyle ps;
    ps.showStats = false;
    ps.logy = false;

    // -------------------------
    // Read graphs by name
    // -------------------------
    auto *g_charge = (TGraphErrors *)fCharge->Get("g_chargeAsym_mt");
    auto *g_RFB_sum = (TGraphErrors *)fFB->Get("g_RFB_mt_sum");
    auto *g_RFB_Wp = (TGraphErrors *)fFB->Get("g_RFB_mt_Wp");
    auto *g_RFB_Wm = (TGraphErrors *)fFB->Get("g_RFB_mt_Wm");

    // -------------------------
    // Read theory graphs by name
    // -------------------------

    auto *g_RFB_Wp_EPPS21 = (TGraphErrors *)fFB_theory->Get("pQCDLightIon_mcfm_nuclearmodfactor_WPlusBoson_RpO_y_dep_EPPS21_FB");
    auto *g_RFB_Wp_nCTEQ15HQ = (TGraphErrors *)fFB_theory->Get("pQCDLightIon_mcfm_nuclearmodfactor_WPlusBoson_RpO_y_dep_nCTEQ15HQ_FB");
    auto *g_RFB_Wp_nNNPDF30 = (TGraphErrors *)fFB_theory->Get("pQCDLightIon_mcfm_nuclearmodfactor_WPlusBoson_RpO_y_dep_nNNPDF30_FB");
    auto *g_RFB_Wp_TUJU21 = (TGraphErrors *)fFB_theory->Get("pQCDLightIon_mcfm_nuclearmodfactor_WPlusBoson_RpO_y_dep_TUJU21_FB");

    auto *g_RFB_Wm_EPPS21 = (TGraphErrors *)fFB_theory->Get("pQCDLightIon_mcfm_nuclearmodfactor_WMinusBoson_RpO_y_dep_EPPS21_FB");
    auto *g_RFB_Wm_nCTEQ15HQ = (TGraphErrors *)fFB_theory->Get("pQCDLightIon_mcfm_nuclearmodfactor_WMinusBoson_RpO_y_dep_nCTEQ15HQ_FB");
    auto *g_RFB_Wm_nNNPDF30 = (TGraphErrors *)fFB_theory->Get("pQCDLightIon_mcfm_nuclearmodfactor_WMinusBoson_RpO_y_dep_nNNPDF30_FB");
    auto *g_RFB_Wm_TUJU21 = (TGraphErrors *)fFB_theory->Get("pQCDLightIon_mcfm_nuclearmodfactor_WMinusBoson_RpO_y_dep_TUJU21_FB");

    if (!g_charge)
        std::cerr
            << "[ERROR] Missing g_chargeAsym_mt\n";
    if (!g_RFB_sum)
        std::cerr << "[ERROR] Missing g_RFB_mt_sum\n";
    if (!g_RFB_Wp)
        std::cerr << "[ERROR] Missing g_RFB_mt_Wp\n";
    if (!g_RFB_Wm)
        std::cerr << "[ERROR] Missing g_RFB_mt_Wm\n";

    // -------------------------
    // Tuners
    // -------------------------
    GraphTuner tuneCharge = [&](TCanvas *c, TGraphErrors *g)
    {
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.2);
        g->SetLineWidth(2);

        // set y-range after axes exist
        if (auto *h = g->GetHistogram())
        {
            h->SetMinimum(-0.6);
            h->SetMaximum(0.6);
        }

        // Optional: line at 0
        double xmin = g->GetXaxis()->GetXmin();
        double xmax = g->GetXaxis()->GetXmax();
        TLine *l0 = new TLine(xmin, 0.0, xmax, 0.0);
        l0->SetLineStyle(2);
        l0->Draw("same");

        c->Modified();
        c->Update();
    };

    GraphTuner tuneRFB = [&](TCanvas *c, TGraphErrors *g)
    {
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.2);
        g->SetLineWidth(2);

        if (auto *h = g->GetHistogram())
        {
            h->SetMinimum(0.3);
            h->SetMaximum(2.1);
        }

        // line at 1
        double xmin = g->GetXaxis()->GetXmin();
        double xmax = g->GetXaxis()->GetXmax();
        TLine *l1 = new TLine(xmin, 1.0, xmax, 1.0);
        l1->SetLineStyle(2);
        l1->Draw("same");

        c->Modified();
        c->Update();
    };

    // -------------------------
    // Charge asym plot
    // -------------------------
    if (g_charge)
    {
        SaveNiceGraph(
            g_charge,
            outChargeDir + "/chargeAsym_mt",
            "#eta^{#mu}_{CM}", // xTitle (change if you want)
            "A_{ch}",          // yTitle
            "",                // mainTitle
            "W #rightarrow #mu #nu",
            "m_{T}-based selection",
            {}, // box lines (or fill later)
            ps,
            tuneCharge);
    }

    // -------------------------
    // R_FB plots (sum, W+, W-)
    // -------------------------
    auto plotRFB = [&](TGraphErrors *g, const std::string &tag, const std::string &subtitle, TGraphErrors *g_theory_1 = nullptr,
                       TGraphErrors *g_theory_2 = nullptr, TGraphErrors *g_theory_3 = nullptr, TGraphErrors *g_theory_4 = nullptr)
    {
        if (!g)
            return;
        SaveNiceGraph_ErrorBand(
            g,
            outFBDir + "/" + tag,
            "#eta^{#mu}_{CM}",
            "R_{FB}",
            "",
            subtitle,
            "m_{T}-based selection",
            {},
            ps,
            tuneRFB, g_theory_1, g_theory_2, g_theory_3, g_theory_4);
    };

    // -------------------------------------------------------------------
    // Weighted-average theory curve for the SUM channel.
    //
    // Theory provides R_FB only for W+ and W- separately. We build the sum
    // theory curve as the DATA-COUNT-weighted average of the two per-charge
    // theory ratios, per |y| bin:
    //
    //     R_FB^sum(|y|) = (N+ R+ + N- R-) / (N+ + N-)
    //
    // where N+, N- are the TOTAL W+ / W- event counts DATA has in that |y| bin
    // -- i.e. forward bin + backward mirror, N+ = F+ + B+, read from the same
    // FB histograms / integration FBratio.C uses. Each charge's theory ratio is
    // weighted by how abundant that charge is in the data, so the sum sits
    // closer to whichever charge dominates. Each theory point's |y| is mapped
    // to the data |y| bin that contains it (nearest center if outside).
    //
    // NOTE: this is the abundance-weighted mean of the two ratios; it is NOT
    // the exact (F+ + F-)/(B+ + B-) identity (the two coincide only when
    // R+ = R-), but it is the intuitive data-driven combination requested.
    // -------------------------------------------------------------------
    std::vector<TGraphErrors *> sumTheory(4, nullptr);
    {
        using pOAnalysis::YieldInRange;

        const char *inFileW = "../skim/rootfile/WToMuNu_pO_PFMet_hist.root";
        TFile *fW = TFile::Open(inFileW, "READ");
        if (!fW || fW->IsZombie())
        {
            std::cerr << "[WARN] Cannot open W skim file for sum-theory weights: "
                      << inFileW << " -> sum theory curves skipped\n";
        }
        else if (!g_RFB_sum)
        {
            std::cerr << "[WARN] Missing g_RFB_mt_sum -> cannot map theory points to |y| bins\n";
        }
        else
        {
            const int Nabs = g_RFB_sum->GetN(); // |y| bins (6 for NY=12)
            const int NY = 2 * Nabs;            // signed-rapidity bins (12)

            std::vector<double> wP(Nabs, 0.5), wM(Nabs, 0.5); // per-bin weights
            std::vector<double> xc(Nabs, 0.0), exc(Nabs, 0.0); // bin center / half-width

            // Total DATA yield of one charge in one signed-rapidity bin.
            auto count = [&](const char *chg, int iy) -> double
            {
                TH1D *h = (TH1D *)fW->Get(Form("h_mt_%s_y%d_FB", chg, iy));
                return YieldInRange(h, 30.0, 200.0, true).value;
            };

            for (int iabs = 0; iabs < Nabs; ++iabs)
            {
                // A |y| bin pools its backward bin and its forward mirror
                // (same pairing as FBratio.C), so the total count of a charge is
                // N = (forward yield) + (backward yield).
                const int iyB = iabs;          // backward bin (negative y_CM)
                const int iyF = NY - 1 - iabs; // forward mirror (positive y_CM)

                const double Np = count("Wp", iyF) + count("Wp", iyB); // N+ = F+ + B+
                const double Nm = count("Wm", iyF) + count("Wm", iyB); // N- = F- + B-
                const double S = Np + Nm;
                if (S > 0.0)
                {
                    wP[iabs] = Np / S;
                    wM[iabs] = Nm / S;
                }
                else
                {
                    std::cerr << "[WARN] Zero W+/W- yield in |y| bin " << iabs
                              << " -> using equal W+/W- weights\n";
                }

                xc[iabs] = g_RFB_sum->GetPointX(iabs);
                exc[iabs] = g_RFB_sum->GetErrorX(iabs);

                std::cout << "[INFO] sum-theory weight |y|bin=" << iabs
                          << " (center " << xc[iabs] << ")"
                          << "  N+=" << Np << "  N-=" << Nm
                          << "  w+=" << wP[iabs] << "  w-=" << wM[iabs] << "\n";
            }

            // Map a theory point |y| to a data |y| bin: containing bin if any,
            // otherwise the nearest bin center.
            auto findBin = [&](double x) -> int
            {
                for (int i = 0; i < Nabs; ++i)
                    if (x >= xc[i] - exc[i] && x <= xc[i] + exc[i])
                        return i;
                int best = 0;
                double bd = 1e30;
                for (int i = 0; i < Nabs; ++i)
                {
                    const double d = std::fabs(x - xc[i]);
                    if (d < bd) { bd = d; best = i; }
                }
                return best;
            };

            // Build the count-weighted theory sum from a W+ / W- theory pair.
            auto buildSumTheory = [&](TGraphErrors *gWp, TGraphErrors *gWm,
                                      const char *name) -> TGraphErrors *
            {
                if (!gWp || !gWm)
                {
                    std::cerr << "[WARN] " << name << ": missing W+ or W- theory -> skipped\n";
                    return nullptr;
                }

                const int N = gWp->GetN();
                const bool aligned = (gWm->GetN() == N);
                if (!aligned)
                    std::cerr << "[WARN] " << name << ": W+ (" << N << ") and W- ("
                              << gWm->GetN() << ") theory have different point counts;"
                              << " interpolating W- with TGraph::Eval\n";

                std::vector<double> vx, vy, vex, vey;
                vx.reserve(N); vy.reserve(N); vex.reserve(N); vey.reserve(N);

                for (int i = 0; i < N; ++i)
                {
                    const double x = gWp->GetPointX(i);
                    const double Rp = gWp->GetPointY(i);
                    const double eRp = gWp->GetErrorY(i);

                    double Rm = aligned ? gWm->GetPointY(i) : gWm->Eval(x);
                    double eRm = aligned ? gWm->GetErrorY(i) : 0.0; // Eval gives no error

                    const int b = findBin(x);
                    const double Rsum = wP[b] * Rp + wM[b] * Rm;
                    const double eRsum = std::sqrt(wP[b] * eRp * wP[b] * eRp +
                                                   wM[b] * eRm * wM[b] * eRm);

                    vx.push_back(x);
                    vy.push_back(Rsum);
                    vex.push_back(gWp->GetErrorX(i));
                    vey.push_back(eRsum);
                }

                TGraphErrors *g = new TGraphErrors((int)vx.size(), vx.data(), vy.data(),
                                                   vex.data(), vey.data());
                g->SetName(name);
                g->SetTitle(name);
                return g;
            };

            sumTheory[0] = buildSumTheory(g_RFB_Wp_EPPS21,    g_RFB_Wm_EPPS21,    "g_RFB_sum_EPPS21");
            sumTheory[1] = buildSumTheory(g_RFB_Wp_nCTEQ15HQ, g_RFB_Wm_nCTEQ15HQ, "g_RFB_sum_nCTEQ15HQ");
            sumTheory[2] = buildSumTheory(g_RFB_Wp_nNNPDF30,  g_RFB_Wm_nNNPDF30,  "g_RFB_sum_nNNPDF30");
            sumTheory[3] = buildSumTheory(g_RFB_Wp_TUJU21,    g_RFB_Wm_TUJU21,    "g_RFB_sum_TUJU21");
        }

        if (fW)
        {
            fW->Close();
            delete fW;
        }
    }

    plotRFB(g_RFB_sum, "RFB_mt_sum", "W #rightarrow #mu #nu (sum)",
            sumTheory[0], sumTheory[1], sumTheory[2], sumTheory[3]);
    plotRFB(g_RFB_Wp, "RFB_mt_Wp", "W^{+} #rightarrow #mu^{+} #nu", g_RFB_Wp_EPPS21,
            g_RFB_Wp_nCTEQ15HQ, g_RFB_Wp_nNNPDF30, g_RFB_Wp_TUJU21);
    plotRFB(g_RFB_Wm, "RFB_mt_Wm", "W^{-} #rightarrow #mu^{-} #bar{#nu}", g_RFB_Wm_EPPS21,
            g_RFB_Wm_nCTEQ15HQ, g_RFB_Wm_nNNPDF30, g_RFB_Wm_TUJU21);

    fCharge->Close();
    fFB->Close();
    delete fCharge;
    delete fFB;

    std::cout << "[OK] Saved plots to: " << outBase << "\n";
}