#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TSystem.h"
#include "plotting_helper.C"

#include <string>
#include <vector>
#include <iostream>

void mtandmet_overlay()
{
    std::string muFile = "../skim/rootfile/WToMuNu_pO_PFMet_hist.root";
    std::string eleFile = "../skim/rootfile/WToElecNu_pO_PFMet_hist.root";

    TFile *fMu = TFile::Open(muFile.c_str(), "READ");
    TFile *fEle = TFile::Open(eleFile.c_str(), "READ");

    if (!fMu || fMu->IsZombie() || !fEle || fEle->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open input files\n";
        return;
    }

    const int NY = 12;
    double yEdges[NY + 1] = {
        -2.4, -2.0, -1.6, -1.2, -0.8, -0.4,
        0.0, 0.4, 0.8, 1.2, 1.6, 2.0,
        2.4};
    double yEdges_FB[NY + 1] = {
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

    // If you already have your yEdges, you can format them here.
    // Otherwise just label by index (y0..y7) for now.
    // Example placeholders:
    std::vector<std::string> yLabel(NY);
    std::vector<std::string> yLabel_FB(NY);
    for (int iy = 0; iy < NY; ++iy)
    {
        yLabel[iy] = RangeLabel("#eta_{#mu}", iy, yEdges, NY);       // or "#eta_{CM}"
        yLabel_FB[iy] = RangeLabel("#eta_{#mu}", iy, yEdges_FB, NY); // or "#eta_{CM}"
    }

    std::string outBase = "./plots_both";
    std::string outMetDir = outBase + "/met";
    std::string outMtDir = outBase + "/mt";

    gSystem->mkdir(outMetDir.c_str(), kTRUE);
    gSystem->mkdir(outMtDir.c_str(), kTRUE);

    PlotStyle ps;
    ps.drawOpt = "hist";
    ps.showStats = false;
    ps.logy = false;
    ps.lm = 0.17;
    ps.yTitleOffset = 1.75;
    ps.boxY1 = 0.52;

    PlotStyle ps1;
    ps1.drawOpt = "hist SAME";
    ps1.showStats = false;
    ps1.logy = false;

    TH1D *h_met_mu_inclusive = nullptr;
    TH1D *h_met_ele_inclusive = nullptr;

    TH1D *h_mt_mu_inclusive = nullptr;
    TH1D *h_mt_ele_inclusive = nullptr;

    PlotTuner commonTuner = [&](TCanvas *c, TH1 *h)
    {
        (void)c;
        if (!h)
            return;

        // Example: ensure y-min is 0 for linear plots
        if (!ps.logy)
            h->SetMinimum(0.0);

        // Example: leave some headroom
        h->SetMaximum(1.25 * h->GetMaximum());
    };

    for (int iy = 0; iy < NY; iy++)
    {

        //--------------------------------
        // MET
        //--------------------------------

        TH1D *h_met_Wp_mu = (TH1D *)fMu->Get(Form("h_met_Wp_y%d", iy));
        TH1D *h_met_Wp_ele = (TH1D *)fEle->Get(Form("h_met_Wp_y%d", iy));

        TH1D *h_met_Wm_mu = (TH1D *)fMu->Get(Form("h_met_Wm_y%d", iy));
        TH1D *h_met_Wm_ele = (TH1D *)fEle->Get(Form("h_met_Wm_y%d", iy));

        //--------------------------------
        // mT
        //--------------------------------

        TH1D *h_mt_Wp_mu = (TH1D *)fMu->Get(Form("h_mt_Wp_y%d", iy));
        TH1D *h_mt_Wp_ele = (TH1D *)fEle->Get(Form("h_mt_Wp_y%d", iy));

        TH1D *h_mt_Wm_mu = (TH1D *)fMu->Get(Form("h_mt_Wm_y%d", iy));
        TH1D *h_mt_Wm_ele = (TH1D *)fEle->Get(Form("h_mt_Wm_y%d", iy));

        if (!h_met_Wp_mu || !h_met_Wp_ele || !h_mt_Wp_mu || !h_mt_Wp_ele)
        {
            std::cout << "[WARN] missing histogram in bin " << iy << "\n";
            continue;
        }

        if (!h_met_mu_inclusive)
        {
            h_met_mu_inclusive = (TH1D *)h_met_Wp_mu->Clone("h_met_mu_inclusive");
            h_met_mu_inclusive->Reset();
        }

        if (!h_met_ele_inclusive)
        {
            h_met_ele_inclusive = (TH1D *)h_met_Wp_ele->Clone("h_met_ele_inclusive");
            h_met_ele_inclusive->Reset();
        }

        if (!h_mt_mu_inclusive)
        {
            h_mt_mu_inclusive = (TH1D *)h_mt_Wp_mu->Clone("h_mt_mu_inclusive");
            h_mt_mu_inclusive->Reset();
        }

        if (!h_mt_ele_inclusive)
        {
            h_mt_ele_inclusive = (TH1D *)h_mt_Wp_ele->Clone("h_mt_ele_inclusive");
            h_mt_ele_inclusive->Reset();
        }

        h_met_mu_inclusive->Add(h_met_Wp_mu);
        h_met_mu_inclusive->Add(h_met_Wm_mu);

        h_met_ele_inclusive->Add(h_met_Wp_ele);
        h_met_ele_inclusive->Add(h_met_Wm_ele);

        h_mt_mu_inclusive->Add(h_mt_Wp_mu);
        h_mt_mu_inclusive->Add(h_mt_Wm_mu);

        h_mt_ele_inclusive->Add(h_mt_Wp_ele);
        h_mt_ele_inclusive->Add(h_mt_Wm_ele);

        //--------------------------------
        // Style
        //--------------------------------

        h_met_Wp_mu->SetLineColor(kBlue + 1);
        h_met_Wp_ele->SetLineColor(kRed + 1);

        h_met_Wm_mu->SetLineColor(kBlue + 1);
        h_met_Wm_ele->SetLineColor(kRed + 1);

        h_mt_Wp_mu->SetLineColor(kBlue + 1);
        h_mt_Wp_ele->SetLineColor(kRed + 1);

        h_mt_Wm_mu->SetLineColor(kBlue + 1);
        h_mt_Wm_ele->SetLineColor(kRed + 1);

        //--------------------------------
        // MET W+
        //--------------------------------
        {
            std::vector<std::string> box = {Form("Passing Events: %.0f, %.0f", h_met_Wp_mu->Integral(1, h_met_Wp_mu->GetNbinsX()),
                                                 h_met_Wp_ele->Integral(1, h_met_Wp_ele->GetNbinsX()))};

            h_met_Wp_mu->Scale(1.0 / h_met_Wp_mu->Integral("width"));
            h_met_Wp_ele->Scale(1.0 / h_met_Wp_ele->Integral("width"));

            SaveNicePlot1D_twoplots(h_met_Wp_mu, h_met_Wp_ele, Form("%s/met_Wp_y%d", outMetDir.c_str(), iy),
                                    "PF MET (GeV)", "Events / 2.0 GeV", "", "W^{+} #rightarrow l^{+} #nu", yLabel[iy], box, ps, ps1, commonTuner);
        }

        //--------------------------------
        // MET W-
        //--------------------------------
        {
            std::vector<std::string> box = {Form("Passing Events: %.0f, %.0f", h_met_Wm_mu->Integral(1, h_met_Wm_mu->GetNbinsX()),
                                                 h_met_Wm_ele->Integral(1, h_met_Wm_ele->GetNbinsX()))};

            h_met_Wm_mu->Scale(1.0 / h_met_Wm_mu->Integral("width"));
            h_met_Wm_ele->Scale(1.0 / h_met_Wm_ele->Integral("width"));

            SaveNicePlot1D_twoplots(h_met_Wm_mu, h_met_Wm_ele, Form("%s/met_Wm_y%d", outMetDir.c_str(), iy),
                                    "PF MET (GeV)", "Events / 2.0 GeV", "", "W^{-} #rightarrow l^{+}#bar{#nu}", yLabel[iy], box, ps, ps1, commonTuner);
        }

        //--------------------------------
        // mT W+
        //--------------------------------

        {
            std::vector<std::string> box = {Form("Passing Events: %.0f, %.0f", h_mt_Wp_mu->Integral(1, h_mt_Wp_mu->GetNbinsX()),
                                                 h_mt_Wp_ele->Integral(1, h_mt_Wp_ele->GetNbinsX()))};

            h_mt_Wp_mu->Scale(1.0 / h_mt_Wp_mu->Integral("width"));
            h_mt_Wp_ele->Scale(1.0 / h_mt_Wp_ele->Integral("width"));

            SaveNicePlot1D_twoplots(h_mt_Wp_mu, h_mt_Wp_ele, Form("%s/mt_Wp_y%d", outMtDir.c_str(), iy),
                                    "m_{T} (GeV)", "Events / 2.5 GeV", "", "W^{+} #rightarrow l^{+}#nu", yLabel[iy], box, ps, ps1, commonTuner);
        }

        //--------------------------------
        // mT W-
        //--------------------------------
        {
            std::vector<std::string> box = {Form("Passing Events: %.0f, %.0f", h_mt_Wm_mu->Integral(1, h_mt_Wm_mu->GetNbinsX()),
                                                 h_mt_Wm_ele->Integral(1, h_mt_Wm_ele->GetNbinsX()))};

            h_mt_Wm_mu->Scale(1.0 / h_mt_Wm_mu->Integral("width"));
            h_mt_Wm_ele->Scale(1.0 / h_mt_Wm_ele->Integral("width"));

            SaveNicePlot1D_twoplots(h_mt_Wm_mu, h_mt_Wm_ele, Form("%s/mt_Wm_y%d", outMtDir.c_str(), iy),
                                    "m_{T} (GeV)", "Events / 2.5 GeV", "", "W^{-} #rightarrow l^{+}#bar{#nu}", yLabel[iy], box, ps, ps1, commonTuner);
        }
    }

    h_met_mu_inclusive->SetLineColor(kBlue + 1);
    h_met_ele_inclusive->SetLineColor(kRed + 1);

    h_mt_mu_inclusive->SetLineColor(kBlue + 1);
    h_mt_ele_inclusive->SetLineColor(kRed + 1);

    {
        std::vector<std::string> box = {Form("Passing Events: %.0f, %.0f", h_met_mu_inclusive->Integral(1, h_met_mu_inclusive->GetNbinsX()),
                                             h_met_ele_inclusive->Integral(1, h_met_ele_inclusive->GetNbinsX()))};

        h_met_mu_inclusive->Scale(1.0 / h_met_mu_inclusive->Integral("width"));
        h_met_ele_inclusive->Scale(1.0 / h_met_ele_inclusive->Integral("width"));

        SaveNicePlot1D_twoplots(h_met_mu_inclusive, h_met_ele_inclusive, Form("%s/h_met_inclusive", outMetDir.c_str()),
                                "PF MET (GeV)", "Events / 2.0 GeV", "", "W #rightarrow l #nu", "inclusive", box, ps, ps1, commonTuner);
    }

    {
        std::vector<std::string> box = {Form("Passing Events: %.0f, %.0f", h_mt_mu_inclusive->Integral(1, h_mt_mu_inclusive->GetNbinsX()),
                                             h_mt_ele_inclusive->Integral(1, h_mt_ele_inclusive->GetNbinsX()))};

        h_mt_mu_inclusive->Scale(1.0 / h_mt_mu_inclusive->Integral("width"));
        h_mt_ele_inclusive->Scale(1.0 / h_mt_ele_inclusive->Integral("width"));

        SaveNicePlot1D_twoplots(h_mt_mu_inclusive, h_mt_ele_inclusive, Form("%s/h_mt_inclusive", outMtDir.c_str()),
                                "m_{T} (GeV)", "Events / 2.5 GeV", "", "W #rightarrow l #nu", "inclusive", box, ps, ps1, commonTuner);
    }

    fMu->Close();
    fEle->Close();

    std::cout << "[INFO] Done. Plots saved in " << outBase << "\n";
}