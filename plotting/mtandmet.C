#include "TCanvas.h"
#include "TH1.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TSystem.h"
#include "plotting_helper.C"

#include <string>
#include <vector>
#include <functional>

void mtandmet(bool isElec = 0)
{
    std::string inFile;
    const char *Channeltype;
    const char *Channeltypewplus;
    const char *Channeltypewminus;
    if (isElec)
    {
        inFile = "../skim/rootfile/WToElecNu_pO_PFMet_hist.root";
        Channeltype = "W #rightarrow e #nu";
        Channeltypewplus = "W^{+} #rightarrow e^{+} #nu";
        Channeltypewminus = "W^{-} #rightarrow e^{-}#bar{#nu}";
    }
    else
    {
        inFile = "../skim/rootfile/WToMuNu_pO_PFMet_hist.root";
        Channeltype = "W #rightarrow #mu #nu";
        Channeltypewplus = "W^{+} #rightarrow #mu^{+} #nu";
        Channeltypewminus = "W^{-} #rightarrow #mu^{-}#bar{#nu}";
    }

    TFile *f = TFile::Open(inFile.c_str(), "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open file: " << inFile << "\n";
        return;
    }

    // output directories
    std::string outBase;
    if (isElec)
    {
        outBase = "./plots/Elec";
    }
    else
    {
        outBase = "./plots";
    }
    const std::string outMetDir = outBase + "/met";
    const std::string outMtDir = outBase + "/mt";
    gSystem->mkdir(outMetDir.c_str(), kTRUE);
    gSystem->mkdir(outMtDir.c_str(), kTRUE);

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

    // Common style
    PlotStyle ps;
    ps.drawOpt = "hist";
    ps.showStats = false;
    ps.logy = false; // set true if you want log-y for MET tails

    // (Optional) one tuner you can reuse for all plots
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

    TH1D *h_met_inclusive = nullptr;
    TH1D *h_mt_inclusive = nullptr;

    // Loop over rapidity bins and charge
    for (int iy = 0; iy < NY; ++iy)
    {
        // ---- MET ----
        TH1D *h_met_Wp = (TH1D *)f->Get(Form("h_met_Wp_y%d", iy));
        TH1D *h_met_Wm = (TH1D *)f->Get(Form("h_met_Wm_y%d", iy));

        TH1D *h_met_Wp_FB = (TH1D *)f->Get(Form("h_met_Wp_y%d_FB", iy));
        TH1D *h_met_Wm_FB = (TH1D *)f->Get(Form("h_met_Wm_y%d_FB", iy));

        // ---- mT ----
        TH1D *h_mt_Wp = (TH1D *)f->Get(Form("h_mt_Wp_y%d", iy));
        TH1D *h_mt_Wm = (TH1D *)f->Get(Form("h_mt_Wm_y%d", iy));

        TH1D *h_mt_Wp_FB = (TH1D *)f->Get(Form("h_mt_Wp_y%d_FB", iy));
        TH1D *h_mt_Wm_FB = (TH1D *)f->Get(Form("h_mt_Wm_y%d_FB", iy));

        if (!h_met_Wp || !h_met_Wm || !h_mt_Wp || !h_mt_Wm)
        {
            std::cerr << "[WARN] Missing histogram(s) at y bin " << iy << "\n";
            continue;
        }

        if (!h_mt_inclusive)
        {
            h_mt_inclusive = (TH1D *)h_mt_Wp->Clone("h_mt_inclusive");
            h_mt_inclusive->Reset();
            h_mt_inclusive->SetDirectory(nullptr);
        }
        if (!h_met_inclusive)
        {
            h_met_inclusive = (TH1D *)h_met_Wp->Clone("h_met_inclusive");
            h_met_inclusive->Reset();
            h_met_inclusive->SetDirectory(nullptr); // avoid ROOT ownership issues
        }

        h_met_inclusive->Add(h_met_Wp);
        h_met_inclusive->Add(h_met_Wm);
        h_mt_inclusive->Add(h_mt_Wp);
        h_mt_inclusive->Add(h_mt_Wm);

        if (h_met_Wp->Integral(1,h_met_Wp->GetNbinsX()) != h_mt_Wp->Integral(1,h_mt_Wp->GetNbinsX()) ){
            cout << "ERROR check plot " << iy << endl;
            cout << h_met_Wp->Integral(1,h_met_Wp->GetNbinsX()) << " And " << h_mt_Wp->Integral(1,h_mt_Wp->GetNbinsX());
        }

        // If you want to close the file later safely, detach histograms:
        // (Not strictly required if you keep file open until end)
        // h_met_Wp->SetDirectory(0); h_met_Wm->SetDirectory(0);
        // h_mt_Wp->SetDirectory(0);  h_mt_Wm->SetDirectory(0);

        // --------- MET plots ----------
        {
            std::vector<std::string> box = {Form("Passing Events: %.0f", h_met_Wp->Integral(1, h_met_Wp->GetNbinsX()))};

            SaveNicePlot1D(
                h_met_Wp,
                outMetDir + Form("/met_Wp_y%d", iy), // outPathNoExt
                "PF MET (GeV)",                      // xTitle
                "Events / 2.0 GeV",                  // yTitle
                "",                                  // mainTitle (whatever you like)
                Channeltypewplus,                    // subTitle1
                yLabel[iy],                          // subTitle2
                box,                                 // info box lines
                ps,
                commonTuner);
        }

        {
            std::vector<std::string> box = {Form("Passing Events: %.0f", h_met_Wm->Integral(1, h_met_Wm->GetNbinsX()))};

            SaveNicePlot1D(
                h_met_Wm,
                outMetDir + Form("/met_Wm_y%d", iy),
                "PF MET (GeV)",
                "Events / 2.0 GeV",
                "",
                Channeltypewminus,
                yLabel[iy],
                box,
                ps,
                commonTuner);
        }

        // --------- Met with FB ratio --------

        {
            std::vector<std::string> box = {Form("Passing Events: %.0f", h_met_Wp_FB->Integral(1, h_met_Wp_FB->GetNbinsX()))};

            SaveNicePlot1D(
                h_met_Wp_FB,
                outMetDir + Form("/met_Wp_y%d_FB", iy), // outPathNoExt
                "PF MET (GeV)",                         // xTitle
                "Events / 2.0 GeV",                     // yTitle
                "",                                     // mainTitle (whatever you like)
                Channeltypewplus,                       // subTitle1
                yLabel_FB[iy],                          // subTitle2
                box,                                    // info box lines
                ps,
                commonTuner);
        }

        {
            std::vector<std::string> box = {Form("Passing Events: %.0f", h_met_Wm_FB->Integral(1, h_met_Wm_FB->GetNbinsX()))};

            SaveNicePlot1D(
                h_met_Wm_FB,
                outMetDir + Form("/met_Wm_y%d_FB", iy),
                "PF MET (GeV)",
                "Events / 2.0 GeV",
                "",
                Channeltypewminus,
                yLabel_FB[iy],
                box,
                ps,
                commonTuner);
        }

        // --------- mT plots ----------
        {
            std::vector<std::string> box = {Form("Passing Events: %.0f", h_mt_Wp->Integral(1, h_mt_Wp->GetNbinsX()))};

            SaveNicePlot1D(
                h_mt_Wp,
                outMtDir + Form("/mt_Wp_y%d", iy),
                "m_{T} (GeV)",
                "Events / 2.5 GeV",
                "",
                Channeltypewplus,
                yLabel[iy],
                box,
                ps,
                commonTuner);
        }

        {
            std::vector<std::string> box = {Form("Passing Events: %.0f", h_mt_Wm->Integral(1, h_mt_Wm->GetNbinsX()))};

            SaveNicePlot1D(
                h_mt_Wm,
                outMtDir + Form("/mt_Wm_y%d", iy),
                "m_{T} (GeV)",
                "Events / 2.5 GeV",
                "",
                Channeltypewminus,
                yLabel[iy],
                box,
                ps,
                commonTuner);
        }

        // --------- mT plots For FB ----------
        {
            std::vector<std::string> box = {Form("Passing Events: %.0f", h_mt_Wp_FB->Integral(1, h_mt_Wp_FB->GetNbinsX()))};

            SaveNicePlot1D(
                h_mt_Wp_FB,
                outMtDir + Form("/mt_Wp_y%d_FB", iy),
                "m_{T} (GeV)",
                "Events / 2.5 GeV",
                "",
                Channeltypewplus,
                yLabel_FB[iy],
                box,
                ps,
                commonTuner);
        }

        {
            std::vector<std::string> box = {Form("Passing Events: %.0f", h_mt_Wm_FB->Integral(1, h_mt_Wm_FB->GetNbinsX()))};

            SaveNicePlot1D(
                h_mt_Wm_FB,
                outMtDir + Form("/mt_Wm_y%d_FB", iy),
                "m_{T} (GeV)",
                "Events / 2.5 GeV",
                "",
                Channeltypewminus,
                yLabel_FB[iy],
                box,
                ps,
                commonTuner);
        }
    }

    {

        std::vector<std::string> box = {Form("Passing Events: %.0f", h_mt_inclusive->Integral(1, h_mt_inclusive->GetNbinsX()))};

        SaveNicePlot1D(
            h_mt_inclusive,
            outMtDir + Form("/h_mt_inclusive"),
            "m_{T} (GeV)",
            "Events / 2.5 GeV",
            "",
            Channeltype,
            "inclusive",
            box,
            ps,
            commonTuner);
    }

    {

        std::vector<std::string> box = {Form("Passing Events: %.0f", h_met_inclusive->Integral(1, h_met_inclusive->GetNbinsX()))};

        SaveNicePlot1D(
            h_met_inclusive,
            outMetDir + Form("/h_met_inclusive"),
            "PF MET (GeV)",
            "Events / 2.0 GeV",
            "",
            Channeltype,
            "inclusive",
            box,
            ps,
            commonTuner);
    }

    f->Close();
    delete f;

    std::cout << "[INFO] Done. Saved plots under: " << outBase << "\n";
}