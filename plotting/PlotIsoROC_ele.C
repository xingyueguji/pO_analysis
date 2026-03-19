#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TParameter.h"
#include "TAxis.h"
#include "TSystem.h"

#include <iostream>
#include <vector>
#include <string>

#include "plotting_helper.C"

// --------------------------------------------
// Convert TGraph -> TGraphErrors
// --------------------------------------------
static TGraphErrors *GraphToGraphErrors(const TGraph *g)
{
    if (!g)
        return nullptr;

    int n = g->GetN();
    TGraphErrors *ge = new TGraphErrors(n);

    const double *x = g->GetX();
    const double *y = g->GetY();

    for (int i = 0; i < n; i++)
    {
        ge->SetPoint(i, x[i], y[i]);
        ge->SetPointError(i, 0, 0);
    }

    return ge;
}

void PlotIsoROC_ele()
{
    //--------------------------------------------
    // Open ROOT file
    //--------------------------------------------
    TFile *h_graph = TFile::Open("../skim/rootfile/IsoStudyOutputs_electron.root", "READ");

    if (!h_graph || h_graph->IsZombie())
    {
        std::cerr << "Cannot open ROOT file\n";
        return;
    }

    gSystem->mkdir("plotsROC_ele", kTRUE);

    //--------------------------------------------
    // Electron ID names
    //--------------------------------------------
    std::vector<std::string> idNames =
        {
            "MVAIdWP80",
            "MVAIdWP85",
            "MVAIdWP90",
            "MVAIdWP95",
            "CutIdWP70",
            "CutIdWP80",
            "CutIdWP90",
            "CutIdWP95"};

    //--------------------------------------------
    // Isolation labels
    //--------------------------------------------
    std::vector<std::string> isoNames =
        {
            "IsoWP80",
            "IsoWP85",
            "IsoWP90",
            "IsoWP95"};

    const int nID = idNames.size();

    //--------------------------------------------
    // Containers
    //--------------------------------------------
    TGraph *g_eff[3][8] = {{nullptr}};
    TGraphErrors *ge_eff[3][8] = {{nullptr}};

    TGraph *roc_g[2][8] = {{nullptr}};
    TGraphErrors *roc_ge[2][8] = {{nullptr}};

    TParameter<double> *roc_integral[2][8] = {{nullptr}};

    TString type[3] = {"OS", "SS", "MET"};

    TGraph *g_passOS[8] = {nullptr};
    TGraph *g_totOS[8] = {nullptr};

    TGraph *g_passSS[8] = {nullptr};
    TGraph *g_totSS[8] = {nullptr};

    TGraph *g_passMET[8] = {nullptr};
    TGraph *g_totMET[8] = {nullptr};

    //--------------------------------------------
    // Load graphs
    //--------------------------------------------
    for (int id = 0; id < nID; id++)
    {
        TString tag = idNames[id];

        for (int i = 0; i < 3; i++)
        {
            g_eff[i][id] =
                (TGraph *)h_graph->Get(Form("eff%s_%s",
                                            type[i].Data(), tag.Data()));

            ge_eff[i][id] = GraphToGraphErrors(g_eff[i][id]);
        }

        g_passOS[id] = (TGraph *)h_graph->Get(Form("passOS_%s", tag.Data()));
        g_totOS[id] = (TGraph *)h_graph->Get(Form("totOS_%s", tag.Data()));

        g_passSS[id] = (TGraph *)h_graph->Get(Form("passSS_%s", tag.Data()));
        g_totSS[id] = (TGraph *)h_graph->Get(Form("totSS_%s", tag.Data()));

        g_passMET[id] = (TGraph *)h_graph->Get(Form("passMET_%s", tag.Data()));
        g_totMET[id] = (TGraph *)h_graph->Get(Form("totMET_%s", tag.Data()));

        roc_g[0][id] =
            (TGraph *)h_graph->Get(Form("rocSS_%s", tag.Data()));

        roc_g[1][id] =
            (TGraph *)h_graph->Get(Form("rocMET_%s", tag.Data()));

        roc_ge[0][id] = GraphToGraphErrors(roc_g[0][id]);
        roc_ge[1][id] = GraphToGraphErrors(roc_g[1][id]);

        roc_integral[0][id] =
            (TParameter<double> *)h_graph->Get(Form("aucSS_%s", tag.Data()));

        roc_integral[1][id] =
            (TParameter<double> *)h_graph->Get(Form("aucMET_%s", tag.Data()));
    }

    //--------------------------------------------
    // Plot styles
    //--------------------------------------------
    PlotStyle psEff;
    psEff.drawOptGraph = "AP";
    psEff.boxY2 = 0.7;
    psEff.boxY1 = 0.55;

    PlotStyle psRoc;
    psRoc.drawOptGraph = "AP";

    //--------------------------------------------
    // Efficiency plots
    //--------------------------------------------
    for (int id = 0; id < nID; id++)
    {

        int currentID = id;

        auto gPassOS = g_passOS[currentID];
        auto gTotOS = g_totOS[currentID];

        auto gPassSS = g_passSS[currentID];
        auto gTotSS = g_totSS[currentID];

        auto gPassMET = g_passMET[currentID];
        auto gTotMET = g_totMET[currentID];

        //--------------------------------------------
        // Efficiency tuner (replace x labels)
        //--------------------------------------------
        GraphTuner Eff_tuner = [&, currentID,
                                gPassOS, gTotOS,
                                gPassSS, gTotSS,
                                gPassMET, gTotMET](TCanvas *c, TGraphErrors *g)
        {
            if (!g)
                return;

            //--------------------------------------------
            // axis formatting (keep yours)
            //--------------------------------------------
            TAxis *ax = g->GetXaxis();
            ax->SetLimits(0.5, 4.5);
            ax->SetNdivisions(4);

            for (int i = 0; i < 4; i++)
                ax->ChangeLabel(i + 1, -1, -1, -1, -1, -1, isoNames[i].c_str());

            g->GetYaxis()->SetRangeUser(0, 1.05);

            g->SetMarkerStyle(20);
            g->SetMarkerSize(1.1);

            //--------------------------------------------
            // draw text
            //--------------------------------------------
            TLatex latex;
            latex.SetTextSize(0.028);

            for (int i = 0; i < g->GetN(); i++)
            {
                double x = g->GetX()[i];

                // --- OS ---
                if (gPassOS && gTotOS)
                {
                    double y = ge_eff[0][currentID]->GetY()[i];
                    double pass = gPassOS->GetY()[i];
                    double tot = gTotOS->GetY()[i];

                    latex.SetTextColor(kBlack);
                    latex.DrawLatex(x + 0.1, y + 0.02,
                                    Form("OS %.0f/%.0f", pass, tot));
                }

                // --- SS ---
                if (gPassSS && gTotSS)
                {
                    double y = ge_eff[1][currentID]->GetY()[i];
                    double pass = gPassSS->GetY()[i];
                    double tot = gTotSS->GetY()[i];

                    latex.SetTextColor(kRed + 1);
                    latex.DrawLatex(x + 0.1, y,
                                    Form("SS %.0f/%.0f", pass, tot));
                }

                // --- MET ---
                if (gPassMET && gTotMET)
                {
                    double y = ge_eff[2][currentID]->GetY()[i];
                    double pass = gPassMET->GetY()[i];
                    double tot = gTotMET->GetY()[i];

                    latex.SetTextColor(kBlue + 1);
                    latex.DrawLatex(x + 0.1, y - 0.02,
                                    Form("MET %.0f/%.0f", pass, tot));
                }
            }

            c->Modified();
            c->Update();
        };

        TGraphErrors *gBkg = nullptr;

        if (g_passSS[id] && g_passMET[id])
        {
            int n = g_passSS[id]->GetN();

            gBkg = new TGraphErrors(n);

            for (int i = 0; i < n; i++)
            {
                double x = g_passSS[id]->GetX()[i];

                double pass = g_passSS[id]->GetY()[i] + g_passMET[id]->GetY()[i];
                double tot = g_totSS[id]->GetY()[i] + g_totMET[id]->GetY()[i];

                double eff = (tot > 0 ? pass / tot : 0);

                gBkg->SetPoint(i, x, eff);
                gBkg->SetPointError(i, 0, 0);
            }
        }

        TGraphErrors *gOS = ge_eff[0][id];
        TGraphErrors *gSS = ge_eff[1][id];
        TGraphErrors *gMET = ge_eff[2][id];

        if (!gOS)
            continue;

        std::vector<std::string> boxLines =
            {
                std::string("Electron ID = ") + idNames[id],
                "Isolation scan"};

        SaveNiceGraph(
            gOS,
            Form("plotsROC_ele/eff_%s", idNames[id].c_str()),
            "Isolation WP",
            "Efficiency",
            "",
            "",
            "",
            boxLines,
            psEff,
            Eff_tuner,
            gSS,
            gMET,
            gBkg);
    }

    //--------------------------------------------
    // ROC plots
    //--------------------------------------------

    //--------------------------------------------
    // ROC tuner
    //--------------------------------------------
    GraphTuner Roc_tuner = [](TCanvas *c, TGraphErrors *g)
    {
        if (!g)
            return;

        g->GetXaxis()->SetRangeUser(0, 1);
        g->GetYaxis()->SetRangeUser(0, 1);

        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.0);

        c->Modified();
        c->Update();
    };

    for (int id = 0; id < nID; id++)
    {
        for (int t = 0; t < 2; t++)
        {
            TGraphErrors *ge = roc_ge[t][id];
            if (!ge)
                continue;

            std::string aucStr = "AUC = (missing)";

            if (roc_integral[t][id])
                aucStr = Form("AUC = %.4f",
                              roc_integral[t][id]->GetVal());

            std::vector<std::string> boxLines =
                {
                    std::string("Electron ID = ") + idNames[id],
                    aucStr};

            SaveNiceGraph(
                ge,
                Form("plotsROC_ele/roc_%s_%s",
                     (t == 0 ? "SS" : "MET"),
                     idNames[id].c_str()),
                "Signal efficiency",
                "Background rejection",
                "",
                "",
                "",
                boxLines,
                psRoc,
                Roc_tuner);
        }
    }

    h_graph->Close();
}