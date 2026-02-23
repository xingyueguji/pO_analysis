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

    TString type[3] = {"OS", "SS", "MET"};
    int case[3] = {0,1,2};
    TString Rvalue[3] = {"0.2","0.3","0.4"};

    TGraph* g_eff[3][3][3];
    TGraphErrors* ge_eff[3][3][3];

    TGraph* roc[2][3][3];
    TParameters* roc_integral[2][3][3];

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++){

        for (int j = 0; j < 3; j++){


            
        }
            
        }
    }
}