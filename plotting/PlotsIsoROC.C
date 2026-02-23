#include "TCanvas.h"
#include "TClass.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TKey.h"
#include "TLegend.h"
#include "TLine.h"
#include "TParameter.h"
#include "TStyle.h"
#include "TSystem.h"

#include <iostream>
#include <string>
#include <vector>

#include "plotting_helper.C" // <-- put your header (the code you pasted) into this

static TGraphErrors *GraphToGraphErrors(const TGraph *g)
{
    if (!g)
        return nullptr;

    const int n = g->GetN();

    TGraphErrors *ge = new TGraphErrors(n);
    ge->SetName(Form("%s", g->GetName()));
    ge->SetTitle(g->GetTitle());

    const double *x = g->GetX();
    const double *y = g->GetY();

    for (int i = 0; i < n; ++i)
    {
        ge->SetPoint(i, x[i], y[i]);
        ge->SetPointError(i, 0.0, 0.0);
    }

    return ge;
}

void PlotsIsoROC()
{
    TFile *h_graph = new TFile("../skim/rootfile/IsoStudyOutputs.root", "READ");
    if (!h_graph || h_graph->IsZombie())
    {
        std::cerr << "ERROR: cannot open IsoStudyOutputs.root\n";
        return;
    }

    TString type[3] = {"OS", "SS", "MET"};     // i = 0,1,2
    int cases[3] = {0, 1, 2};                  // j = 0,1,2
    TString Rvalue[3] = {"0.2", "0.3", "0.4"}; // k = 0,1,2

    TGraph *g_eff[3][3][3] = {{{nullptr}}};
    TGraphErrors *ge_eff[3][3][3] = {{{nullptr}}};

    TGraph *roc_g[2][3][3] = {{{nullptr}}}; // i-1 = 0,1 for SS/MET if i!=0
    TParameter<double> *roc_integral[2][3][3] = {{{nullptr}}};

    // --- Load graphs
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                g_eff[i][j][k] = (TGraph *)h_graph->Get(
                    Form("eff%s_case%i_R%s", type[i].Data(), cases[j], Rvalue[k].Data()));

                ge_eff[i][j][k] = GraphToGraphErrors(g_eff[i][j][k]);

                if (i != 0)
                {
                    roc_g[i - 1][j][k] = (TGraph *)h_graph->Get(
                        Form("roc%s_case%i_R%s", type[i].Data(), cases[j], Rvalue[k].Data()));

                    roc_integral[i - 1][j][k] = (TParameter<double> *)h_graph->Get(
                        Form("auc%s_case%i_R%s", type[i].Data(), cases[j], Rvalue[k].Data()));
                }
            }
        }
    }

    // --- Define styles (you can tune these using your plotting_helper.C)
    PlotStyle psEff;
    psEff.drawOptGraph = "AP"; // example: axis+line+points
    psEff.boxX1 = 0.5;
    psEff.boxY1 = 0.33;
    psEff.boxX2 = 0.9;
    psEff.boxY2 = 0.53;

    PlotStyle psRoc;
    psRoc.drawOptGraph = "AP";
    psRoc.boxX1 = 0.5;
    psRoc.boxY1 = 0.33;
    psRoc.boxX2 = 0.9;
    psRoc.boxY2 = 0.53;

    GraphTuner Roc_tuner = [](TCanvas *c, TGraphErrors *h)
    {
        // Axis ranges
        h->GetXaxis()->SetRangeUser(0.0, 1.0);
        h->GetYaxis()->SetRangeUser(0.0, 1.0);

        // Solid black markers
        h->SetMarkerColor(kBlack);
        h->SetLineColor(kBlack);

        // Filled marker (20 = solid circle)
        h->SetMarkerStyle(20);

        // Larger points
        h->SetMarkerSize(1.0);

        c->Modified();
        c->Update();
    };

    // --- Plot efficiency overlays: for each (j,k) overlay i=0,1,2
    for (int j = 0; j < 3; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            TGraphErrors *g0 = ge_eff[0][j][k];
            TGraphErrors *g1 = ge_eff[1][j][k];
            TGraphErrors *g2 = ge_eff[2][j][k];

            if (!g0)
                continue;

            std::vector<std::string> boxLines = {
                std::string("case = ") + std::to_string(cases[j]),
                std::string("R = ") + Rvalue[k].Data(),
                "overlay: OS / SS / MET"};

            std::string out = Form("plotsROC/eff_case%i_R%s_overlay",
                                   cases[j], Rvalue[k].Data());

            SaveNiceGraph(g0,
                          out,
                          "x", "Efficiency",
                          "",
                          "",
                          "",
                          boxLines,
                          psEff,
                          Roc_tuner,
                          g1,
                          g2);
        }
    }

    // --- Plot ROC separately: i=1,2 only (stored in roc_g[0]=SS, roc_g[1]=MET)
    for (int i = 1; i < 3; i++)
    {
        int ii = i - 1; // 0 or 1
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                TGraph *g = roc_g[ii][j][k];
                if (!g)
                    continue;

                // Convert to TGraphErrors for SaveNiceGraph
                TGraphErrors *ge = GraphToGraphErrors(g);

                std::string aucStr = "AUC = (missing)";
                if (roc_integral[ii][j][k])
                    aucStr = Form("AUC = %.4f", roc_integral[ii][j][k]->GetVal());

                std::vector<std::string> boxLines = {
                    std::string("type = ") + type[i].Data(),
                    std::string("case = ") + std::to_string(cases[j]),
                    std::string("R = ") + Rvalue[k].Data(),
                    aucStr};

                std::string out = Form("plotsROC/roc_%s_case%i_R%s",
                                       type[i].Data(), cases[j], Rvalue[k].Data());

                SaveNiceGraph(ge,
                              out,
                              "Signal efficiency", "Background rejection",
                              "",
                              "",
                              "",
                              boxLines,
                              psRoc,
                              Roc_tuner);

                delete ge; // since GraphToGraphErrors allocates
            }
        }
    }

    h_graph->Close();
    delete h_graph;
}