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

void mtandmet()
{
    const std::string inFile = "../skim/rootfile/WToMuNu_pO_PFMet_hist.root";
    TFile *f = TFile::Open(inFile.c_str(), "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open file: " << inFile << "\n";
        return;
    }

    // output directories
    const std::string outBase = "./plots";
    const std::string outMetDir = outBase + "/met";
    const std::string outMtDir = outBase + "/mt";
    gSystem->mkdir(outMetDir.c_str(), kTRUE);
    gSystem->mkdir(outMtDir.c_str(), kTRUE);

    const int NY = 16;
    double yEdges[NY + 1] = {
        -2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3,
        0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4};

    double yEdges_FB[NY + 1] = {
        -1.7068, -1.4501, -1.1934, -0.9368,
        -0.6801, -0.4234, -0.1667, 0.0899,
        0.3466, 0.6033, 0.8599, 1.1166,
        1.3733, 1.6300, 1.8866, 2.1433,
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
                "PF MET [GeV]",                      // xTitle
                "Events",                            // yTitle
                "",                                  // mainTitle (whatever you like)
                "W^{+} #rightarrow #mu^{+}#nu",      // subTitle1
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
                "PF MET [GeV]",
                "Events",
                "",
                "W^{-} #rightarrow #mu^{-}#bar{#nu}",
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
                "PF MET [GeV]",                         // xTitle
                "Events",                               // yTitle
                "",                                     // mainTitle (whatever you like)
                "W^{+} #rightarrow #mu^{+}#nu",         // subTitle1
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
                "PF MET [GeV]",
                "Events",
                "",
                "W^{-} #rightarrow #mu^{-}#bar{#nu}",
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
                "m_{T} [GeV]",
                "Events",
                "",
                "W^{+} #rightarrow #mu^{+}#nu",
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
                "m_{T} [GeV]",
                "Events",
                "",
                "W^{-} #rightarrow #mu^{-}#bar{#nu}",
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
                "m_{T} [GeV]",
                "Events",
                "",
                "W^{+} #rightarrow #mu^{+}#nu",
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
                "m_{T} [GeV]",
                "Events",
                "",
                "W^{-} #rightarrow #mu^{-}#bar{#nu}",
                yLabel_FB[iy],
                box,
                ps,
                commonTuner);
        }
    }

    f->Close();
    delete f;

    std::cout << "[INFO] Done. Saved plots under: " << outBase << "\n";
}