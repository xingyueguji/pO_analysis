#include "TFile.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TLine.h"
#include <iostream>
#include <string>
#include <vector>

#include "plotting_helper.C" // where your PlotStyle + SaveNicePlotGraph live

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

    plotRFB(g_RFB_sum, "RFB_mt_sum", "W #rightarrow #mu #nu (sum)");
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