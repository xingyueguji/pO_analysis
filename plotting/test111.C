#include "TCanvas.h"
#include "TH1.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TFile.h"
#include "plotting_helper.C"

#include <string>
#include <vector>
#include <iostream>

void plot_qcd_sideband(bool isElec = 0)
{
    if (isElec)
    {
        std::cout << "[INFO] Electron channel not implemented yet.\n";
        return;
    }

    std::string inFile = "../skim/rootfile/WToMuNu_pO_PFMet_hist.root";

    TFile *f = TFile::Open(inFile.c_str(), "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open file: " << inFile << "\n";
        return;
    }

    const int NISO = 5;
    double isoEdges[NISO + 1] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

    const std::string outBase = "./plots/qcd_sideband";
    gSystem->mkdir(outBase.c_str(), kTRUE);

    PlotStyle ps;
    ps.drawOpt = "E";
    ps.showStats = false;
    ps.logy = false;
    ps.boxY1 = 0.62;
    ps.boxY2 = 0.82;

    PlotTuner commonTuner = [&](TCanvas *c, TH1 *h)
    {
        (void)c;
        if (!h)
            return;

        h->SetMinimum(0.0);
        h->SetMaximum(1.30 * h->GetMaximum());
    };

    TH1D *h_met_iso_muPlus_inclusive = nullptr;
    TH1D *h_met_iso_muMinus_inclusive = nullptr;

    for (int ib = 0; ib < NISO; ++ib)
    {
        TH1D *h_muPlus  = (TH1D *)f->Get(Form("h_met_iso_muPlus_bin%d", ib));
        TH1D *h_muMinus = (TH1D *)f->Get(Form("h_met_iso_muMinus_bin%d", ib));

        if (!h_muPlus || !h_muMinus)
        {
            std::cerr << "[WARN] Missing sideband histograms in iso bin " << ib << "\n";
            continue;
        }

        if (!h_met_iso_muPlus_inclusive)
        {
            h_met_iso_muPlus_inclusive = (TH1D *)h_muPlus->Clone("h_met_iso_muPlus_inclusive");
            h_met_iso_muPlus_inclusive->Reset();
            h_met_iso_muPlus_inclusive->SetDirectory(nullptr);

            h_met_iso_muMinus_inclusive = (TH1D *)h_muMinus->Clone("h_met_iso_muMinus_inclusive");
            h_met_iso_muMinus_inclusive->Reset();
            h_met_iso_muMinus_inclusive->SetDirectory(nullptr);
        }

        h_met_iso_muPlus_inclusive->Add(h_muPlus);
        h_met_iso_muMinus_inclusive->Add(h_muMinus);

        std::string isoLabel = Form("Muon isolation #in [%.1f, %.1f)", isoEdges[ib], isoEdges[ib + 1]);

        // mu+
        {
            std::vector<std::string> box = {
                Form("Passing Events: %.0f", h_muPlus->Integral(1, h_muPlus->GetNbinsX()))
            };

            SaveNicePlot1D(
                h_muPlus,
                outBase + Form("/met_iso_muPlus_bin%d", ib),
                "PF MET (GeV)",
                "Events / 2.0 GeV",
                "",
                "QCD sideband, #mu^{+}",
                isoLabel,
                box,
                ps,
                commonTuner);
        }

        // mu-
        {
            std::vector<std::string> box = {
                Form("Passing Events: %.0f", h_muMinus->Integral(1, h_muMinus->GetNbinsX()))
            };

            SaveNicePlot1D(
                h_muMinus,
                outBase + Form("/met_iso_muMinus_bin%d", ib),
                "PF MET (GeV)",
                "Events / 2.0 GeV",
                "",
                "QCD sideband, #mu^{-}",
                isoLabel,
                box,
                ps,
                commonTuner);
        }
    }

    // inclusive mu+
    if (h_met_iso_muPlus_inclusive)
    {
        std::vector<std::string> box = {
            Form("Passing Events: %.0f", h_met_iso_muPlus_inclusive->Integral(1, h_met_iso_muPlus_inclusive->GetNbinsX()))
        };

        SaveNicePlot1D(
            h_met_iso_muPlus_inclusive,
            outBase + "/met_iso_muPlus_inclusive",
            "PF MET (GeV)",
            "Events / 2.0 GeV",
            "",
            "QCD sideband, #mu^{+}",
            "inclusive in isolation sideband",
            box,
            ps,
            commonTuner);
    }

    // inclusive mu-
    if (h_met_iso_muMinus_inclusive)
    {
        std::vector<std::string> box = {
            Form("Passing Events: %.0f", h_met_iso_muMinus_inclusive->Integral(1, h_met_iso_muMinus_inclusive->GetNbinsX()))
        };

        SaveNicePlot1D(
            h_met_iso_muMinus_inclusive,
            outBase + "/met_iso_muMinus_inclusive",
            "PF MET (GeV)",
            "Events / 2.0 GeV",
            "",
            "QCD sideband, #mu^{-}",
            "inclusive in isolation sideband",
            box,
            ps,
            commonTuner);
    }

    f->Close();
    delete f;

    std::cout << "[INFO] Done. Saved plots under: " << outBase << "\n";
}