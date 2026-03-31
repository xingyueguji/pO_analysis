#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "plotting_helper.C"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

struct QCDFitConfig
{
    double fitMin = 0.0;
    double fitMax = 80.0;

    double p0 = 1.0;
    double p1 = 20.0;
    double p2 = 10.0;
    double p3 = 0.0;

    double p0min = 0.0, p0max = 1e9;
    double p1min = -20.0, p1max = 30.0;
    double p2min = -20.0, p2max = 30.0;
    double p3min = -1.0, p3max = 1.0;
};

// ---------------------------------------------------------
// QCD shape
// p[0] = normalization
// p[1] = sigma0
// p[2] = sigma1
// p[3] = sigma2
// ---------------------------------------------------------
double QCDRayleighLike(double *x, double *p)
{
    const double xx = x[0];
    if (xx <= 0)
        return 0.0;

    const double xt = xx / 50.0 - 1.0;
    const double sigma = p[1] + p[2] * xt + p[3] * (2.0 * xt * xt - 1.0);

    if (sigma <= 0)
        return 0.0;

    return p[0] * xx * std::exp(-(xx * xx) / (2.0 * sigma * sigma));
}

void FitAndDrawOneHist(TH1D *h,
                       const std::string &outPathNoExt,
                       const std::string &sub1,
                       const std::string &sub2,
                       double &sigma0, double &sigma0err,
                       double &sigma1, double &sigma1err,
                       double &sigma2, double &sigma2err,
                       QCDFitConfig cfg,
                       double fitMin = 0.0,
                       double fitMax = 100.0)
{
    if (!h)
        return;

    PlotStyle ps;
    ps.drawOpt = "E";
    ps.showStats = false;
    ps.logy = true;
    ps.boxX1 = 0.56;
    ps.boxY1 = 0.52;
    ps.boxX2 = 0.92;
    ps.boxY2 = 0.74;
    ps.headerX = 0.52;
    ps.headerY = 0.92;

    PlotTuner tuner = [&](TCanvas *c, TH1 *hh)
    {
        (void)c;
        if (!hh)
            return;
        hh->SetMinimum(0.01);
        hh->SetMaximum(10.5 * hh->GetMaximum());
        hh->SetMarkerStyle(20);
        hh->SetMarkerSize(1.2);
        hh->SetLineColor(kBlack);
    };

    // -------------------------
    // Define fit function
    // -------------------------
    TF1 *fQCD = new TF1(Form("fQCD_%s", h->GetName()), QCDRayleighLike, fitMin, fitMax, 4);
    fQCD->SetParName(0, "Norm");
    fQCD->SetParName(1, "sigma0");
    fQCD->SetParName(2, "sigma1");
    fQCD->SetParName(3, "sigma2");

    // Initial values
    double peakx = h->GetBinCenter(h->GetMaximumBin());
    double maxy = h->GetMaximum();
    fQCD->SetParameters(maxy / std::max(peakx, 1.0), cfg.p1, cfg.p2, cfg.p3);

    // Loose parameter limits for stability
    fQCD->SetParLimits(0, 0.0, 1e9);
    fQCD->SetParLimits(1, cfg.p1min, cfg.p1max);
    fQCD->SetParLimits(2, cfg.p2min, cfg.p2max);
    fQCD->SetParLimits(3, cfg.p3min, cfg.p3max);

    TFitResultPtr fitRes = h->Fit(fQCD, "RSQ"); // R=range, S=result, Q=quiet

    sigma0 = fQCD->GetParameter(1);
    sigma0err = fQCD->GetParError(1);
    sigma1 = fQCD->GetParameter(2);
    sigma1err = fQCD->GetParError(2);
    sigma2 = fQCD->GetParameter(3);
    sigma2err = fQCD->GetParError(3);

    double chi2 = fQCD->GetChisquare();
    int ndf = fQCD->GetNDF();
    double chi2ndf = (ndf > 0 ? chi2 / ndf : 0.0);

    // -------------------------
    // Draw plot with fit
    // -------------------------
    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas(Form("c_%s", h->GetName()), "", ps.w, ps.h);
    ApplyCanvasStyle(c, ps);
    c->cd();

    ApplyHistStyle(h, ps, "PF MET (GeV)", "Events / 2.0 GeV");
    h->Draw(ps.drawOpt.c_str());

    fQCD->SetLineColor(kRed + 1);
    fQCD->SetLineWidth(3);
    fQCD->Draw("SAME");

    std::vector<std::string> boxLines;
    boxLines.push_back(Form("Yield = %.0f", h->Integral(1, h->GetNbinsX())));
    boxLines.push_back(Form("#sigma_{0} = %.3f #pm %.3f", sigma0, sigma0err));
    boxLines.push_back(Form("#sigma_{1} = %.3f #pm %.3f", sigma1, sigma1err));
    boxLines.push_back(Form("#sigma_{2} = %.3f #pm %.3f", sigma2, sigma2err));
    boxLines.push_back(Form("#chi^{2}/ndf = %.2f", chi2ndf));
    if ((int)fitRes != 0)
        boxLines.push_back(Form("fit status = %d", (int)fitRes));

    DrawHeader(ps, "", sub1, sub2);
    DrawInfoBox(ps, boxLines);

    if (tuner)
        tuner(c, h);

    CMS_lumi(c, 13, 10);
    c->Modified();
    c->Update();

    c->SaveAs((outPathNoExt + ".png").c_str());
    c->SaveAs((outPathNoExt + ".pdf").c_str());

    delete c;
    delete fQCD;
}

void MakeExtrapPlot(const std::string &name,
                    const std::string &outPathNoExt,
                    const std::string &sub1,
                    const std::vector<double> &x,
                    const std::vector<double> &y,
                    const std::vector<double> &ey,
                    const std::string &yTitle,
                    double xEval,
                    double &yEval,
                    double &yEvalErr)
{
    const int n = x.size();
    TGraphErrors *g = new TGraphErrors(n);

    for (int i = 0; i < n; ++i)
    {
        g->SetPoint(i, x[i], y[i]);
        g->SetPointError(i, 0.0, ey[i]);
    }

    g->SetName(name.c_str());
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.2);
    g->SetLineWidth(2);

    TF1 *fLin = new TF1(Form("fLin_%s", name.c_str()), "pol1", 0.01, 1.1);
    TFitResultPtr r = g->Fit(fLin, "RSQ");

    yEval = fLin->Eval(xEval);

    // error propagation using fit covariance for pol1
    // y = p0 + p1*x
    // Var(y) = Var(p0) + x^2 Var(p1) + 2x Cov(p0,p1)
    double var = 0.0;
    if ((int)r == 0)
    {
        const TMatrixDSym &cov = r->GetCovarianceMatrix();
        const double var0 = cov(0, 0);
        const double var1 = cov(1, 1);
        const double cov01 = cov(0, 1);
        var = var0 + xEval * xEval * var1 + 2.0 * xEval * cov01;
        if (var < 0)
            var = 0;
    }
    yEvalErr = std::sqrt(var);

    TGraphErrors *gExtr = new TGraphErrors(1);
    gExtr->SetPoint(0, xEval, yEval);
    gExtr->SetPointError(0, 0.0, yEvalErr);

    gExtr->SetMarkerStyle(29); // star / diamond-like
    gExtr->SetMarkerSize(2.0);
    gExtr->SetMarkerColor(kBlue + 2);
    gExtr->SetLineColor(kBlue + 2);
    gExtr->SetLineWidth(2);

    PlotStyle ps;
    ps.drawOptGraph = "AP";
    ps.showStats = false;
    ps.logy = false;
    ps.boxX1 = 0.18;
    ps.boxY1 = 0.52;
    ps.boxX2 = 0.52;
    ps.boxY2 = 0.74;

    GraphTuner tuner = [&](TCanvas *c, TGraphErrors *gg)
    {
        (void)c;
        if (!gg)
            return;
        gg->GetXaxis()->SetLimits(-0.1, 1.2);
        gg->SetMinimum(*std::min_element(y.begin(), y.end()) - 0.55 * std::fabs(*std::min_element(y.begin(), y.end())));
        gg->SetMaximum(*std::max_element(y.begin(), y.end()) + 0.75 * std::fabs(*std::max_element(y.begin(), y.end())));
        gPad->SetFrameLineWidth(3);
    };

    std::vector<std::string> boxLines;
    boxLines.push_back("Linear fit");
    boxLines.push_back(Form("value @ iso=%.2f", xEval));
    boxLines.push_back(Form("= %.4f #pm %.4f", yEval, yEvalErr));

    /*SaveNiceGraph(g,
                  outPathNoExt,
                  "Muon isolation",
                  yTitle,
                  "",
                  sub1,
                  "",
                  boxLines,
                  ps,
                  tuner);*/

    // draw line / point on top in a second pass
    TCanvas *c = new TCanvas(Form("c_extra_%s", name.c_str()), "", ps.w, ps.h);
    ApplyCanvasStyle(c, ps);
    c->cd();
    ApplyGraphStyle(g, ps, "Muon isolation", yTitle);
    g->Draw("AP");
    fLin->SetLineColor(kRed + 1);
    fLin->SetLineWidth(3);
    fLin->Draw("SAME");
    gExtr->Draw("PE SAME");

    TLine *l = new TLine(xEval, g->GetYaxis()->GetXmin(), xEval, g->GetYaxis()->GetXmax());
    l->SetLineStyle(2);
    l->SetLineColor(kBlue + 1);

    DrawHeader(ps, "", sub1, "");
    DrawInfoBox(ps, boxLines);

    tuner(c, g);

    const double ymin = gPad->GetUymin();
    const double ymax = gPad->GetUymax();
    TLine *lv = new TLine(xEval, ymin, xEval, ymax);
    lv->SetLineStyle(2);
    lv->SetLineColor(kBlue + 1);
    lv->SetLineWidth(2);
    lv->Draw("SAME");

    CMS_lumi(c, 13, 10);
    c->Modified();
    c->Update();

    c->SaveAs((outPathNoExt + "_fit.png").c_str());
    c->SaveAs((outPathNoExt + "_fit.pdf").c_str());

    delete lv;
    delete l;
    delete c;
    delete fLin;
    delete g;
}

void qcd_sideband_fit_and_extrapolate(bool isElec = 0)
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

    const int NISO = 3;
    double isoEdges[NISO + 1] = {0.1, 0.4, 0.7, 1.0};
    double isoCenter[NISO] = {0.25, 0.55, 0.85};
    const double isoSignal = 0.03;

    const std::string outBase = "./plots/qcd_sideband_fit";
    gSystem->mkdir(outBase.c_str(), kTRUE);

    std::vector<double> xIso;
    std::vector<double> sig0_p, sig1_p, sig2_p, esig0_p, esig1_p, esig2_p;
    std::vector<double> sig0_m, sig1_m, sig2_m, esig0_m, esig1_m, esig2_m;

    // Here I create intial value for different fits

    static const int NCHG = 2; // 0=muPlus, 1=muMinus

    QCDFitConfig fitCfg[NCHG][NISO];

    // muPlus, bin 1
    fitCfg[0][0].fitMax = 70.0;
    fitCfg[0][0].p1 = 18.0;
    fitCfg[0][0].p2 = 8.0;
    fitCfg[0][0].p3 = 0.0;

    fitCfg[0][0].p1min = 10.0;
    fitCfg[0][0].p1max = 30.0;

    fitCfg[0][0].p2min = 0.0;
    fitCfg[0][0].p2max = 20.0;

    fitCfg[0][0].p3min = -0.2;
    fitCfg[0][0].p3max = 0.2;

    // muMinus, bin 1
    fitCfg[1][0].fitMax = 70.0;
    fitCfg[1][0].p1 = 20.0;
    fitCfg[1][0].p2 = 10.0;
    fitCfg[1][0].p3 = 0.0;

    fitCfg[1][0].p1min = 12.0;
    fitCfg[1][0].p1max = 28.0;

    fitCfg[1][0].p2min = 0.0;
    fitCfg[1][0].p2max = 18.0;

    fitCfg[1][0].p3min = -0.2;
    fitCfg[1][0].p3max = 0.5;

    // muPlus, bin 2
    fitCfg[0][1].fitMax = 70.0;
    fitCfg[0][1].p1 = 18.0;
    fitCfg[0][1].p2 = 8.0;
    fitCfg[0][1].p3 = 0.0;

    fitCfg[0][1].p1min = 10.0;
    fitCfg[0][1].p1max = 30.0;

    fitCfg[0][1].p2min = 0.0;
    fitCfg[0][1].p2max = 20.0;

    fitCfg[0][1].p3min = -0.2;
    fitCfg[0][1].p3max = 0.2;

    // muMinus, bin 2
    fitCfg[1][1].fitMax = 70.0;
    fitCfg[1][1].p1 = 20.0;
    fitCfg[1][1].p2 = 10.0;
    fitCfg[1][1].p3 = 0.0;

    fitCfg[1][1].p1min = 12.0;
    fitCfg[1][1].p1max = 28.0;

    fitCfg[1][1].p2min = 0.0;
    fitCfg[1][1].p2max = 18.0;

    fitCfg[1][1].p3min = -0.2;
    fitCfg[1][1].p3max = 0.2;

    // muPlus, bin 3
    fitCfg[0][2].fitMax = 70.0;
    fitCfg[0][2].p1 = 18.0;
    fitCfg[0][2].p2 = 8.0;
    fitCfg[0][2].p3 = 0.0;

    fitCfg[0][2].p1min = 10.0;
    fitCfg[0][2].p1max = 30.0;

    fitCfg[0][2].p2min = 0.0;
    fitCfg[0][2].p2max = 20.0;

    fitCfg[0][2].p3min = -0.2;
    fitCfg[0][2].p3max = 0.2;

    // muMinus, bin 3
    fitCfg[1][2].fitMax = 70.0;
    fitCfg[1][2].p1 = 20.0;
    fitCfg[1][2].p2 = 10.0;
    fitCfg[1][2].p3 = 0.0;

    fitCfg[1][2].p1min = 12.0;
    fitCfg[1][2].p1max = 28.0;

    fitCfg[1][2].p2min = 0.0;
    fitCfg[1][2].p2max = 18.0;

    fitCfg[1][2].p3min = -0.2;
    fitCfg[1][2].p3max = 0.2;

    // ------------------------------------------------
    // Fit each sideband histogram and save overlay plot
    // ------------------------------------------------
    for (int ib = 0; ib < NISO; ++ib)
    {
        TH1D *h_muPlus = (TH1D *)f->Get(Form("h_met_iso_muPlus_bin%d", ib));
        TH1D *h_muMinus = (TH1D *)f->Get(Form("h_met_iso_muMinus_bin%d", ib));

        if (!h_muPlus || !h_muMinus)
        {
            std::cerr << "[WARN] Missing histogram in iso bin " << ib << "\n";
            continue;
        }

        double s0, e0, s1, e1, s2, e2;

        FitAndDrawOneHist(h_muPlus,
                          outBase + Form("/fit_met_iso_muPlus_bin%d", ib),
                          "QCD sideband, #mu^{+}",
                          Form("Muon isolation #in [%.1f, %.1f)", isoEdges[ib], isoEdges[ib + 1]),
                          s0, e0, s1, e1, s2, e2, fitCfg[0][ib]);

        xIso.push_back(isoCenter[ib]);
        sig0_p.push_back(s0);
        esig0_p.push_back(e0);
        sig1_p.push_back(s1);
        esig1_p.push_back(e1);
        sig2_p.push_back(s2);
        esig2_p.push_back(e2);

        FitAndDrawOneHist(h_muMinus,
                          outBase + Form("/fit_met_iso_muMinus_bin%d", ib),
                          "QCD sideband, #mu^{-}",
                          Form("Muon isolation #in [%.1f, %.1f)", isoEdges[ib], isoEdges[ib + 1]),
                          s0, e0, s1, e1, s2, e2, fitCfg[1][ib]);

        sig0_m.push_back(s0);
        esig0_m.push_back(e0);
        sig1_m.push_back(s1);
        esig1_m.push_back(e1);
        sig2_m.push_back(s2);
        esig2_m.push_back(e2);
    }

    // ------------------------------------------------
    // Extrapolation plots
    // ------------------------------------------------
    double yEval = 0.0, yEvalErr = 0.0;

    // mu+
    MakeExtrapPlot("g_sigma0_muPlus",
                   outBase + "/sigma0_vs_iso_muPlus",
                   "QCD sideband, #mu^{+}",
                   xIso, sig0_p, esig0_p,
                   "#sigma_{0}",
                   isoSignal, yEval, yEvalErr);
    std::cout << "[mu+] sigma0(0.03) = " << yEval << " +/- " << yEvalErr << "\n";

    MakeExtrapPlot("g_sigma1_muPlus",
                   outBase + "/sigma1_vs_iso_muPlus",
                   "QCD sideband, #mu^{+}",
                   xIso, sig1_p, esig1_p,
                   "#sigma_{1}",
                   isoSignal, yEval, yEvalErr);
    std::cout << "[mu+] sigma1(0.03) = " << yEval << " +/- " << yEvalErr << "\n";

    MakeExtrapPlot("g_sigma2_muPlus",
                   outBase + "/sigma2_vs_iso_muPlus",
                   "QCD sideband, #mu^{+}",
                   xIso, sig2_p, esig2_p,
                   "#sigma_{2}",
                   isoSignal, yEval, yEvalErr);
    std::cout << "[mu+] sigma2(0.03) = " << yEval << " +/- " << yEvalErr << "\n";

    // mu-
    MakeExtrapPlot("g_sigma0_muMinus",
                   outBase + "/sigma0_vs_iso_muMinus",
                   "QCD sideband, #mu^{-}",
                   xIso, sig0_m, esig0_m,
                   "#sigma_{0}",
                   isoSignal, yEval, yEvalErr);
    std::cout << "[mu-] sigma0(0.03) = " << yEval << " +/- " << yEvalErr << "\n";

    MakeExtrapPlot("g_sigma1_muMinus",
                   outBase + "/sigma1_vs_iso_muMinus",
                   "QCD sideband, #mu^{-}",
                   xIso, sig1_m, esig1_m,
                   "#sigma_{1}",
                   isoSignal, yEval, yEvalErr);
    std::cout << "[mu-] sigma1(0.03) = " << yEval << " +/- " << yEvalErr << "\n";

    MakeExtrapPlot("g_sigma2_muMinus",
                   outBase + "/sigma2_vs_iso_muMinus",
                   "QCD sideband, #mu^{-}",
                   xIso, sig2_m, esig2_m,
                   "#sigma_{2}",
                   isoSignal, yEval, yEvalErr);
    std::cout << "[mu-] sigma2(0.03) = " << yEval << " +/- " << yEvalErr << "\n";

    f->Close();
    delete f;

    std::cout << "[INFO] Done. Saved plots under: " << outBase << "\n";
}