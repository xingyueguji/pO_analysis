// correction/plot_iso_summary.C
//
// Per-ID isolation summary plots in one consistent cosmetic:
//   LEFT  pad: continuous efficiency vs relIso cut (signal OS, QCD MET<5, SS)
//   RIGHT pad: the corresponding ROC (signal eff vs QCD rejection) + AUC,
//              with the operating-point marker at the chosen cut.
//
// Electron: one canvas per ID (8), reading the continuous-relIso scan added to
//           isolation_ele.C  -> plotsROC_ele/summary_<ID>.png
// Muon:     one canvas (single tight ID), reading the pure tight-ID scan from
//           isolation_mu_tight.C -> plotsROC/summary_muon.png
//           (muon needs no ID study; shown the same way for consistency).
//
// Both scans use the Delta-beta PU-corrected relIso (MuonPOG convention),
// (ch + max(0, neu + pho - 0.5*PU))/pT -- same as the skim's RelIsoPF
// (switched 2026-07-06; before that both were the uncorrected sum).
//
// Inputs (run isolation_mu_tight.C+ / isolation_ele.C+ first):
//   rootfile/IsoStudyOutputs_electron.root    (contEff*_<ID>, contRocMET_<ID>, contAucMET_<ID>)
//   rootfile/IsoStudyOutputs_muon_tight.root  (contEff*_dbeta, contRocMET_dbeta, contAucMET_dbeta)
//
// Run from correction/:  root -l -q 'plot_iso_summary.C+'

#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMarker.h"
#include "TAxis.h"
#include "TParameter.h"
#include "TStyle.h"
#include "TSystem.h"

#include <string>
#include <vector>

namespace
{
// linear interpolation of an (increasing-x) TGraph at x
double interp(TGraph *g, double x)
{
    if (!g) return -1;
    const double *X = g->GetX();
    const double *Y = g->GetY();
    const int n = g->GetN();
    if (x <= X[0]) return Y[0];
    if (x >= X[n - 1]) return Y[n - 1];
    for (int i = 1; i < n; ++i)
        if (X[i] >= x) { double t = (x - X[i - 1]) / (X[i] - X[i - 1]); return Y[i - 1] + t * (Y[i] - Y[i - 1]); }
    return Y[n - 1];
}

// One ID: eff-vs-cut (left) + ROC (right) on a 2-pad canvas, shared cosmetic.
void drawOne(const std::string &outNoExt, const std::string &title,
             TGraph *gOS, TGraph *gMET, TGraph *gSS, TGraph *gROC,
             double auc, double cut, double isoXmax)
{
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas(("c_" + outNoExt).c_str(), "", 1500, 640);
    c->Divide(2, 1);

    const int cSig = kGreen + 2, cQCD = kRed + 1, cSS = kGray + 2, cROC = kAzure + 2, cCut = kBlue + 1;

    // ---------------- LEFT: efficiency vs relIso cut ----------------
    c->cd(1);
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.08);
    gPad->SetGridx(); gPad->SetGridy();
    gOS->SetLineColor(cSig); gOS->SetLineWidth(3); gOS->SetMarkerColor(cSig);
    gMET->SetLineColor(cQCD); gMET->SetLineWidth(3);
    gOS->SetTitle(";relIso (#Delta#beta-corr.) cut;efficiency");
    gOS->GetXaxis()->SetLimits(0.0, isoXmax);
    gOS->GetYaxis()->SetRangeUser(0.0, 1.08);
    gOS->GetYaxis()->SetTitleOffset(1.4);
    gOS->Draw("AL");
    gMET->Draw("L SAME");
    if (gSS) { gSS->SetLineColor(cSS); gSS->SetLineStyle(2); gSS->SetLineWidth(2); gSS->Draw("L SAME"); }

    TLine *l = new TLine(cut, 0.0, cut, 1.08);
    l->SetLineColor(cCut); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw();

    const double eS = interp(gOS, cut), eB = interp(gMET, cut);
    TLatex tt; tt.SetNDC();
    tt.SetTextSize(0.040); tt.DrawLatex(0.13, 0.94, title.c_str());
    tt.SetTextSize(0.034); tt.SetTextColor(cCut);
    tt.DrawLatex(0.46, 0.86, Form("cut: relIso < %.3f", cut));
    tt.SetTextColor(kBlack);
    tt.DrawLatex(0.46, 0.81, Form("#varepsilon_{sig}=%.3f   #varepsilon_{QCD}=%.3f", eS, eB));

    TLegend *lg = new TLegend(0.44, 0.17, 0.90, 0.36);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.033);
    lg->AddEntry(gOS, "signal  #varepsilon (Z#rightarrow#font[12]{ll} OS)", "l");
    lg->AddEntry(gMET, "QCD  #varepsilon (single-#font[12]{l}, MET<5)", "l");
    if (gSS) lg->AddEntry(gSS, "SS bkg  #varepsilon", "l");
    lg->Draw();

    // ---------------- RIGHT: ROC + AUC ----------------
    c->cd(2);
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.08);
    gPad->SetGridx(); gPad->SetGridy();
    gROC->SetLineColor(cROC); gROC->SetLineWidth(3);
    gROC->SetTitle(";signal efficiency;QCD rejection (1 - #varepsilon_{MET})");
    gROC->GetXaxis()->SetLimits(0.5, 1.0);
    gROC->GetYaxis()->SetRangeUser(0.0, 1.0);
    gROC->GetYaxis()->SetTitleOffset(1.4);
    gROC->Draw("AL");

    TLine *diag = new TLine(0.5, 0.5, 1.0, 0.0); // random-classifier reference
    diag->SetLineColor(kGray + 1); diag->SetLineStyle(3); diag->Draw();

    TMarker *m = new TMarker(eS, 1.0 - eB, 29);
    m->SetMarkerColor(kBlack); m->SetMarkerSize(2.8); m->Draw();

    tt.SetTextSize(0.040); tt.DrawLatex(0.13, 0.94, title.c_str());
    TLegend *lg2 = new TLegend(0.16, 0.16, 0.66, 0.33);
    lg2->SetBorderSize(0); lg2->SetFillStyle(0); lg2->SetTextSize(0.033);
    lg2->AddEntry(gROC, Form("ROC vs QCD  (AUC = %.3f)", auc), "l");
    lg2->AddEntry(m, Form("operating pt (relIso < %.3f)", cut), "p");
    lg2->Draw();

    c->SaveAs((outNoExt + ".png").c_str());
    c->SaveAs((outNoExt + ".pdf").c_str());
    delete c;
}
} // namespace

void plot_iso_summary()
{
    gSystem->mkdir("plotsROC_ele", kTRUE);
    gSystem->mkdir("plotsROC", kTRUE);

    // ---------------- ELECTRON: per ID ----------------
    TFile *fe = TFile::Open("rootfile/IsoStudyOutputs_electron.root", "READ");
    if (fe && !fe->IsZombie())
    {
        const std::vector<std::string> ids = {
            "MVAIdWP80", "MVAIdWP85", "MVAIdWP90", "MVAIdWP95",
            "CutIdWP70", "CutIdWP80", "CutIdWP90", "CutIdWP95"};
        for (const auto &id : ids)
        {
            TGraph *gOS = (TGraph *)fe->Get(Form("contEffOS_%s", id.c_str()));
            TGraph *gMET = (TGraph *)fe->Get(Form("contEffMET_%s", id.c_str()));
            TGraph *gSS = (TGraph *)fe->Get(Form("contEffSS_%s", id.c_str()));
            TGraph *gROC = (TGraph *)fe->Get(Form("contRocMET_%s", id.c_str()));
            TParameter<double> *pa = (TParameter<double> *)fe->Get(Form("contAucMET_%s", id.c_str()));
            if (!gOS || !gMET || !gROC) { printf("[WARN] missing graphs for %s\n", id.c_str()); continue; }
            const std::string note = (id == "MVAIdWP95") ? "  [current]" : (id == "CutIdWP95") ? "  [old]" : "";
            drawOne("plotsROC_ele/summary_" + id, "W#rightarrowe#nu   ID = " + id + note,
                    gOS, gMET, gSS, gROC, pa ? pa->GetVal() : -1.0, 0.095, 0.50);
        }
        fe->Close();
    }
    else printf("[WARN] no rootfile/IsoStudyOutputs_electron.root (run isolation_ele.C+ first)\n");

    // ---------------- MUON: single tight ID, Delta-beta relIso (= skim RelIsoPF) ----------------
    TFile *fm = TFile::Open("rootfile/IsoStudyOutputs_muon_tight.root", "READ");
    if (fm && !fm->IsZombie())
    {
        TGraph *gOS = (TGraph *)fm->Get("contEffOS_dbeta");
        TGraph *gMET = (TGraph *)fm->Get("contEffMET_dbeta");
        TGraph *gSS = (TGraph *)fm->Get("contEffSS_dbeta");
        TGraph *gROC = (TGraph *)fm->Get("contRocMET_dbeta");
        TParameter<double> *pa = (TParameter<double> *)fm->Get("contAucMET_dbeta");
        if (gOS && gMET && gROC)
            drawOne("plotsROC/summary_muon", "W#rightarrow#mu#nu   ID = muIDTight",
                    gOS, gMET, gSS, gROC, pa ? pa->GetVal() : -1.0, 0.15, 0.50);
        else printf("[WARN] missing muon dbeta graphs\n");
        fm->Close();
    }
    else printf("[WARN] no rootfile/IsoStudyOutputs_muon_tight.root (run isolation_mu_tight.C+ first)\n");

    printf("[DONE] iso summary -> plotsROC_ele/summary_<ID>.png , plotsROC/summary_muon.png\n");
}
