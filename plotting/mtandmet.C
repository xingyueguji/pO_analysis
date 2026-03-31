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

void mtandmet(bool isElec = 1)
{
    std::string inFile;
    std::string inFile_MC_signal_Wp;
    std::string inFile_MC_signal_Wm;
    std::string inFile_MC_background_Z;

    const char *Channeltype;
    const char *Channeltypewplus;
    const char *Channeltypewminus;

    if (isElec)
    {
        inFile = "../skim/rootfile/WToElecNu_pO_PFMet_hist.root";
        inFile_MC_signal_Wp = "../skim/rootfile/WToElecNu_pO_PFMet_Wp_hist.root";
        inFile_MC_signal_Wm = "../skim/rootfile/WToElecNu_pO_PFMet_Wm_hist.root";
        inFile_MC_background_Z = "../skim/rootfile/WToElecNu_pO_PFMet_DY_hist.root";
        Channeltype = "W #rightarrow e #nu";
        Channeltypewplus = "W^{+} #rightarrow e^{+} #nu";
        Channeltypewminus = "W^{-} #rightarrow e^{-}#bar{#nu}";
    }
    else
    {
        inFile = "../skim/rootfile/WToMuNu_pO_PFMet_hist.root";
        inFile_MC_signal_Wp = "../skim/rootfile/WToMuNu_pO_PFMet_Wp_hist.root";
        inFile_MC_signal_Wm = "../skim/rootfile/WToMuNu_pO_PFMet_Wm_hist.root";
        inFile_MC_background_Z = "../skim/rootfile/WToMuNu_pO_PFMet_DY_hist.root";
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

    TFile *f_Wp = TFile::Open(inFile_MC_signal_Wp.c_str(), "READ");
    if (!f_Wp || f_Wp->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open file: " << inFile_MC_signal_Wp << "\n";
        return;
    }

    TFile *f_Wm = TFile::Open(inFile_MC_signal_Wm.c_str(), "READ");
    if (!f_Wm || f_Wm->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open file: " << inFile_MC_signal_Wm << "\n";
        return;
    }

    TFile *f_DY = TFile::Open(inFile_MC_background_Z.c_str(), "READ");
    if (!f_DY || f_DY->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open file: " << inFile_MC_background_Z << "\n";
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
    ps.boxY1 = 0.62;
    ps.boxY2 = 0.82;

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

    TH1D *h_met_inclusive_MC_signal = nullptr;
    TH1D *h_mt_inclusive_MC_signal = nullptr;

    TH1D *h_met_inclusive_MC_Z = nullptr;
    TH1D *h_mt_inclusive_MC_Z = nullptr;

    // Loop over rapidity bins and charge
    for (int iy = 0; iy < NY; ++iy)
    {
        // ---- MET ----
        TH1D *h_met_Wp = (TH1D *)f->Get(Form("h_met_Wp_y%d", iy));
        TH1D *h_met_Wm = (TH1D *)f->Get(Form("h_met_Wm_y%d", iy));

        TH1D *h_met_Wp_FB = (TH1D *)f->Get(Form("h_met_Wp_y%d_FB", iy));
        TH1D *h_met_Wm_FB = (TH1D *)f->Get(Form("h_met_Wm_y%d_FB", iy));

        // ---- MET MC Signal ----

        TH1D *h_met_Wp_MC_Wp = (TH1D *)f_Wp->Get(Form("h_met_Wp_y%d", iy));
        TH1D *h_met_Wm_MC_Wp = (TH1D *)f_Wp->Get(Form("h_met_Wm_y%d", iy));

        TH1D *h_met_Wp_FB_MC_Wp = (TH1D *)f_Wp->Get(Form("h_met_Wp_y%d_FB", iy));
        TH1D *h_met_Wm_FB_MC_Wp = (TH1D *)f_Wp->Get(Form("h_met_Wm_y%d_FB", iy));

        TH1D *h_met_Wp_MC_Wm = (TH1D *)f_Wm->Get(Form("h_met_Wp_y%d", iy));
        TH1D *h_met_Wm_MC_Wm = (TH1D *)f_Wm->Get(Form("h_met_Wm_y%d", iy));

        TH1D *h_met_Wp_FB_MC_Wm = (TH1D *)f_Wm->Get(Form("h_met_Wp_y%d_FB", iy));
        TH1D *h_met_Wm_FB_MC_Wm = (TH1D *)f_Wm->Get(Form("h_met_Wm_y%d_FB", iy));

        // ---- MET MC Z BK ----

        TH1D *h_met_Wp_MC_Z = (TH1D *)f_DY->Get(Form("h_met_Wp_y%d", iy));
        TH1D *h_met_Wm_MC_Z = (TH1D *)f_DY->Get(Form("h_met_Wm_y%d", iy));

        TH1D *h_met_Wp_FB_MC_Z = (TH1D *)f_DY->Get(Form("h_met_Wp_y%d_FB", iy));
        TH1D *h_met_Wm_FB_MC_Z = (TH1D *)f_DY->Get(Form("h_met_Wm_y%d_FB", iy));

        // ---- mT ----
        TH1D *h_mt_Wp = (TH1D *)f->Get(Form("h_mt_Wp_y%d", iy));
        TH1D *h_mt_Wm = (TH1D *)f->Get(Form("h_mt_Wm_y%d", iy));

        TH1D *h_mt_Wp_FB = (TH1D *)f->Get(Form("h_mt_Wp_y%d_FB", iy));
        TH1D *h_mt_Wm_FB = (TH1D *)f->Get(Form("h_mt_Wm_y%d_FB", iy));

        // ---- mT MC Signal ----

        TH1D *h_mt_Wp_MC_Wp = (TH1D *)f_Wp->Get(Form("h_mt_Wp_y%d", iy));
        TH1D *h_mt_Wm_MC_Wp = (TH1D *)f_Wp->Get(Form("h_mt_Wm_y%d", iy));

        TH1D *h_mt_Wp_FB_MC_Wp = (TH1D *)f_Wp->Get(Form("h_mt_Wp_y%d_FB", iy));
        TH1D *h_mt_Wm_FB_MC_Wp = (TH1D *)f_Wp->Get(Form("h_mt_Wm_y%d_FB", iy));

        TH1D *h_mt_Wp_MC_Wm = (TH1D *)f_Wm->Get(Form("h_mt_Wp_y%d", iy));
        TH1D *h_mt_Wm_MC_Wm = (TH1D *)f_Wm->Get(Form("h_mt_Wm_y%d", iy));

        TH1D *h_mt_Wp_FB_MC_Wm = (TH1D *)f_Wm->Get(Form("h_mt_Wp_y%d_FB", iy));
        TH1D *h_mt_Wm_FB_MC_Wm = (TH1D *)f_Wm->Get(Form("h_mt_Wm_y%d_FB", iy));

        // ---- mT MC Z BK ----

        TH1D *h_mt_Wp_MC_Z = (TH1D *)f_DY->Get(Form("h_mt_Wp_y%d", iy));
        TH1D *h_mt_Wm_MC_Z = (TH1D *)f_DY->Get(Form("h_mt_Wm_y%d", iy));

        TH1D *h_mt_Wp_FB_MC_Z = (TH1D *)f_DY->Get(Form("h_mt_Wp_y%d_FB", iy));
        TH1D *h_mt_Wm_FB_MC_Z = (TH1D *)f_DY->Get(Form("h_mt_Wm_y%d_FB", iy));

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

            h_mt_inclusive_MC_signal = (TH1D *)h_mt_Wp->Clone("h_mt_inclusive_MC_signal");
            h_mt_inclusive_MC_signal->Reset();
            h_mt_inclusive_MC_signal->SetDirectory(nullptr);

            h_mt_inclusive_MC_Z = (TH1D *)h_mt_Wp->Clone("h_mt_inclusive_MC_Z");
            h_mt_inclusive_MC_Z->Reset();
            h_mt_inclusive_MC_Z->SetDirectory(nullptr);
        }
        if (!h_met_inclusive)
        {
            h_met_inclusive = (TH1D *)h_met_Wp->Clone("h_met_inclusive");
            h_met_inclusive->Reset();
            h_met_inclusive->SetDirectory(nullptr); // avoid ROOT ownership issues

            h_met_inclusive_MC_signal = (TH1D *)h_met_Wp->Clone("h_met_inclusive_MC_signal");
            h_met_inclusive_MC_signal->Reset();
            h_met_inclusive_MC_signal->SetDirectory(nullptr); // avoid ROOT ownership issues

            h_met_inclusive_MC_Z = (TH1D *)h_met_Wp->Clone("h_met_inclusive_MC_Z");
            h_met_inclusive_MC_Z->Reset();
            h_met_inclusive_MC_Z->SetDirectory(nullptr); // avoid ROOT ownership issues
        }

        h_met_inclusive->Add(h_met_Wp);
        h_met_inclusive->Add(h_met_Wm);

        h_met_inclusive_MC_signal->Add(h_met_Wp_MC_Wp);
        h_met_inclusive_MC_signal->Add(h_met_Wm_MC_Wp);
        h_met_inclusive_MC_signal->Add(h_met_Wp_MC_Wm);
        h_met_inclusive_MC_signal->Add(h_met_Wm_MC_Wm);

        h_met_inclusive_MC_Z->Add(h_met_Wp_MC_Z);
        h_met_inclusive_MC_Z->Add(h_met_Wm_MC_Z);

        h_mt_inclusive->Add(h_mt_Wp);
        h_mt_inclusive->Add(h_mt_Wm);

        h_mt_inclusive_MC_signal->Add(h_mt_Wp_MC_Wp);
        h_mt_inclusive_MC_signal->Add(h_mt_Wm_MC_Wp);
        h_mt_inclusive_MC_signal->Add(h_mt_Wp_MC_Wm);
        h_mt_inclusive_MC_signal->Add(h_mt_Wm_MC_Wm);

        h_mt_inclusive_MC_Z->Add(h_mt_Wp_MC_Z);
        h_mt_inclusive_MC_Z->Add(h_mt_Wm_MC_Z);

        if (h_met_Wp->Integral(1, h_met_Wp->GetNbinsX()) != h_mt_Wp->Integral(1, h_mt_Wp->GetNbinsX()))
        {
            cout << "ERROR check plot " << iy << endl;
            cout << h_met_Wp->Integral(1, h_met_Wp->GetNbinsX()) << " And " << h_mt_Wp->Integral(1, h_mt_Wp->GetNbinsX());
        }

        // If you want to close the file later safely, detach histograms:
        // (Not strictly required if you keep file open until end)
        // h_met_Wp->SetDirectory(0); h_met_Wm->SetDirectory(0);
        // h_mt_Wp->SetDirectory(0);  h_mt_Wm->SetDirectory(0);

        // --------- MET plots ----------
        {
            std::vector<std::string> box = {Form("Passing Events: %.0f", h_met_Wp->Integral(1, h_met_Wp->GetNbinsX()))};

            std::vector<TH1 *> bkgs = {
                h_met_Wp_MC_Wp,
                h_met_Wp_MC_Wm,
                h_met_Wp_MC_Z};

            std::vector<std::string> names = {
                "W+",
                "W-",
                "DY"};

            SaveNicePlot1D_WithBkg(
                h_met_Wp,
                bkgs,
                names,
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

            std::vector<TH1 *> bkgs = {
                h_met_Wm_MC_Wp,
                h_met_Wm_MC_Wm,
                h_met_Wm_MC_Z};

            std::vector<std::string> names = {
                "W+",
                "W-",
                "DY"};

            SaveNicePlot1D_WithBkg(
                h_met_Wm,
                bkgs,
                names,
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

            std::vector<TH1 *> bkgs = {
                h_met_Wp_FB_MC_Wp,
                h_met_Wp_FB_MC_Wm,
                h_met_Wp_FB_MC_Z};

            std::vector<std::string> names = {
                "W+",
                "W-",
                "DY"};

            SaveNicePlot1D_WithBkg(
                h_met_Wp_FB,
                bkgs,
                names,
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

            std::vector<TH1 *> bkgs = {
                h_met_Wm_FB_MC_Wp,
                h_met_Wm_FB_MC_Wm,
                h_met_Wm_FB_MC_Z};

            std::vector<std::string> names = {
                "W+",
                "W-",
                "DY"};

            SaveNicePlot1D_WithBkg(
                h_met_Wm_FB,
                bkgs,
                names,
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

            std::vector<TH1 *> bkgs = {
                h_mt_Wp_MC_Wp,
                h_mt_Wp_MC_Wm,
                h_mt_Wp_MC_Z};

            std::vector<std::string> names = {
                "W+",
                "W-",
                "DY"};

            SaveNicePlot1D_WithBkg(
                h_mt_Wp,
                bkgs,
                names,
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

            std::vector<TH1 *> bkgs = {
                h_mt_Wm_MC_Wp,
                h_mt_Wm_MC_Wm,
                h_mt_Wm_MC_Z};

            std::vector<std::string> names = {
                "W+",
                "W-",
                "DY"};

            SaveNicePlot1D_WithBkg(
                h_mt_Wm,
                bkgs,
                names,
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

            std::vector<TH1 *> bkgs = {
                h_mt_Wp_FB_MC_Wp,
                h_mt_Wp_FB_MC_Wm,
                h_mt_Wp_FB_MC_Z};

            std::vector<std::string> names = {
                "W+",
                "W-",
                "DY"};

            SaveNicePlot1D_WithBkg(
                h_mt_Wp_FB,
                bkgs,
                names,
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

            std::vector<TH1 *> bkgs = {
                h_mt_Wm_FB_MC_Wp,
                h_mt_Wm_FB_MC_Wm,
                h_mt_Wm_FB_MC_Z};

            std::vector<std::string> names = {
                "W+",
                "W-",
                "DY"};

            SaveNicePlot1D_WithBkg(
                h_mt_Wm_FB,
                bkgs,
                names,
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

        std::vector<TH1 *> bkgs = {
            h_mt_inclusive_MC_signal,
            h_mt_inclusive_MC_Z};

        std::vector<std::string> names = {
            "W+/W-",
            "DY"};

        SaveNicePlot1D_WithBkg(
            h_mt_inclusive,
            bkgs,
            names,
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

        std::vector<TH1 *> bkgs = {
            h_met_inclusive_MC_signal,
            h_met_inclusive_MC_Z};

        std::vector<std::string> names = {
            "W+/W-",
            "DY"};

        SaveNicePlot1D_WithBkg(
            h_met_inclusive,
            bkgs,
            names,
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