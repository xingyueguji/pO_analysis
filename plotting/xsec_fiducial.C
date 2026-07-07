#include "TCanvas.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TLatex.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

#include "TGraph.h"
#include "TFile.h"
#include "plotting_helper.C"           // PlotStyle, ApplyCanvasStyle/HistStyle, DrawHeader, CMS_lumi
#include "../skim/mc_norm.h"           // pONorm::kLumi_invnb (single-source data lumi)

// =============================================================================
// xsec_fiducial.C -- fiducial W cross sections (W+, W-, W inclusive) from the
// fitted signal yields in the Combine summary CSVs, muon and electron overlaid.
//
//   sigma_fid = N_fit / L_int        (N_fit = fitted signal yield from the CSV)
//
// NOTE: this is the fiducial cross section BEFORE the lepton efficiency
// correction (eff/reco/ID/iso/trigger not yet applied), i.e. sigma_fid x eps.
// That is exactly why muon and electron disagree here -- the gap is the channel
// efficiency difference; once eps is applied they should converge (universality).
// Stat (fit) uncertainty only; lumi uncertainty not included.
//
//   root -l -q 'xsec_fiducial.C+'
// =============================================================================

// read fitted_yield (col 5) + fitted_yield_err (col 6) for a region row.
static bool readYield(const char *csv, const char *region, double &y, double &e)
{
    std::ifstream in(csv);
    if (!in) { std::cerr << "[ERROR] cannot open CSV: " << csv << "\n"; return false; }
    std::string line;
    std::getline(in, line); // header
    while (std::getline(in, line))
    {
        std::stringstream ss(line);
        std::string f;
        std::vector<std::string> c;
        while (std::getline(ss, f, ',')) c.push_back(f);
        if (c.size() >= 7 && c[0] == region && c[1] == "r")
        {
            y = std::atof(c[5].c_str());
            e = std::atof(c[6].c_str());
            return true;
        }
    }
    std::cerr << "[WARN] region '" << region << "' not found in " << csv << "\n";
    return false;
}

void xsec_fiducial(
    const char *muCsv  = "../../HiggsAnalysis-CombinedLimit/test/pO_fit_out/mu/summary/mu_summary.csv",
    const char *eleCsv = "../../HiggsAnalysis-CombinedLimit/test/pO_fit_out/ele/summary/ele_summary.csv",
    double Lint = pONorm::kLumi_invnb /* 46.5 nb^-1 */)
{
    gStyle->SetEndErrorSize(4);

    const char *regions[3] = {"Wp_incl", "Wm_incl", "W_incl"};
    const char *xlab[3]    = {"W^{+}", "W^{-}", "W"};

    double yMu[3] = {0, 0, 0}, eMu[3] = {0, 0, 0};
    double yEl[3] = {0, 0, 0}, eEl[3] = {0, 0, 0};
    bool haveMu = false, haveEl = false;
    for (int i = 0; i < 3; ++i)
    {
        double y, e;
        if (readYield(muCsv, regions[i], y, e))  { yMu[i] = y / Lint; eMu[i] = e / Lint; haveMu = true; }
        if (readYield(eleCsv, regions[i], y, e)) { yEl[i] = y / Lint; eEl[i] = e / Lint; haveEl = true; }
        printf("[xsec] %-8s  mu = %6.2f +/- %5.2f nb   e = %6.2f +/- %5.2f nb\n",
               regions[i], yMu[i], eMu[i], yEl[i], eEl[i]);
    }
    if (!haveMu && !haveEl) { std::cerr << "[ERROR] no yields read -- check the CSV paths.\n"; return; }

    // points, offset slightly in x so mu/e error bars don't overlap
    double xMu[3] = {1 - 0.08, 2 - 0.08, 3 - 0.08};
    double xEl[3] = {1 + 0.08, 2 + 0.08, 3 + 0.08};
    double ex[3]  = {0, 0, 0};
    TGraphErrors *gMu = new TGraphErrors(3, xMu, yMu, ex, eMu);
    TGraphErrors *gEl = new TGraphErrors(3, xEl, yEl, ex, eEl);

    PlotStyle ps;
    ps.logy = false;
    TCanvas *c = new TCanvas("c_xsec_fid", "", ps.w, ps.h);
    ApplyCanvasStyle(c, ps);
    c->cd();

    double ymax = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        if (haveMu) ymax = std::max(ymax, yMu[i] + eMu[i]);
        if (haveEl) ymax = std::max(ymax, yEl[i] + eEl[i]);
    }

    TH1F *fr = new TH1F("fr_xsec", "", 3, 0.5, 3.5); // category frame
    for (int i = 0; i < 3; ++i) fr->GetXaxis()->SetBinLabel(i + 1, xlab[i]);
    fr->SetStats(0);
    ApplyHistStyle(fr, ps, "", "#sigma^{fid}_{W} (nb)");
    fr->GetXaxis()->SetLabelSize(0.065);
    fr->GetXaxis()->SetNdivisions(3);
    fr->SetMinimum(0.0);
    fr->SetMaximum(1.45 * ymax);
    fr->Draw();

    if (haveMu)
    {
        gMu->SetMarkerStyle(20); gMu->SetMarkerSize(1.5);
        gMu->SetMarkerColor(kBlack); gMu->SetLineColor(kBlack); gMu->SetLineWidth(2);
        gMu->Draw("P SAME");
    }
    if (haveEl)
    {
        gEl->SetMarkerStyle(21); gEl->SetMarkerSize(1.5);
        gEl->SetMarkerColor(kRed + 1); gEl->SetLineColor(kRed + 1); gEl->SetLineWidth(2);
        gEl->Draw("P SAME");
    }

    DrawHeader(ps, "", "W #rightarrow l #nu", "fiducial cross section");

    // caveat note (this is sigma_fid x eps, no efficiency correction yet)
    TLatex tx;
    tx.SetNDC();
    tx.SetTextFont(42);
    tx.SetTextSize(0.030);
    tx.DrawLatex(0.18, 0.255, "stat. unc. only"); // L is already in the CMS_lumi label
    tx.DrawLatex(0.18, 0.21, "no eff and acc. correction");

    TLegend *leg = new TLegend(0.20, 0.66, 0.44, 0.82); // upper-left (header is upper-right)
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.038);
    if (haveMu) leg->AddEntry(gMu, "Muon", "lep");
    if (haveEl) leg->AddEntry(gEl, "Electron", "lep");
    leg->Draw();

    CMS_lumi(c, 13, 10);
    c->RedrawAxis();
    c->Modified();
    c->Update();

    const std::string outDir = "./plots/xsec";
    gSystem->mkdir(outDir.c_str(), kTRUE);
    c->SaveAs((outDir + "/W_fiducial.png").c_str());
    c->SaveAs((outDir + "/W_fiducial.pdf").c_str());
    std::cout << "[OK] Saved " << outDir << "/W_fiducial.{png,pdf}  (L = " << Lint << " nb^-1)\n";
}

// -----------------------------------------------------------------------------
// per-bin yields for one charge+binning from <chan>_W_yields.csv:
//   region,charge,binning,ybin,r,rErr,signal_prefit,fitted_yield(7),fitted_yield_err(8),...
static bool readBinYields(const char *csv, const char *charge, const char *binning,
                          double y[12], double e[12])
{
    for (int i = 0; i < 12; ++i) { y[i] = 0; e[i] = 0; }
    std::ifstream in(csv);
    if (!in) { std::cerr << "[ERROR] cannot open CSV: " << csv << "\n"; return false; }
    std::string line;
    std::getline(in, line); // header
    int found = 0;
    while (std::getline(in, line))
    {
        std::stringstream ss(line);
        std::string f;
        std::vector<std::string> c;
        while (std::getline(ss, f, ',')) c.push_back(f);
        if (c.size() < 9) continue;
        if (c[1] == charge && c[2] == binning)
        {
            int iy = std::atoi(c[3].c_str());
            if (iy >= 0 && iy < 12) { y[iy] = std::atof(c[7].c_str()); e[iy] = std::atof(c[8].c_str()); ++found; }
        }
    }
    return found > 0;
}

// =============================================================================
// xsec_fiducial_diff -- differential fiducial cross section d(sigma_fid)/d(eta_CM)
// for W+, W-, and W inclusive, vs the CM-frame lepton pseudorapidity, ONE channel.
// Per bin: d(sigma)/d(eta) = N_fit / (L * d-eta). The eta bin centers/widths are
// taken from g_chargeAsym_mt (charge_asym_fit_<chan>.root) so the axis matches the
// charge-asymmetry plot exactly (lab bins shifted to the CM frame). Same caveats
// as xsec_fiducial: sigma_fid x eps (no efficiency correction), stat unc only.
//
//   root -l -q 'xsec_fiducial.C+' -e 'xsec_fiducial_diff(false)'   // muon; (true)=electron
// =============================================================================
void xsec_fiducial_diff(bool isElec = false,
                        const char *csv = nullptr,
                        const char *chargeFile = nullptr,
                        double Lint = pONorm::kLumi_invnb)
{
    gStyle->SetEndErrorSize(3);
    const char *chan = isElec ? "ele" : "mu";
    const char *lepSym = isElec ? "e" : "#mu";

    const TString sCsv = csv ? TString(csv)
        : TString::Format("../../HiggsAnalysis-CombinedLimit/test/pO_fit_out/%s/summary/%s_W_yields.csv", chan, chan);
    const TString sCharge = chargeFile ? TString(chargeFile)
        : TString::Format("../skim/rootfile/charge_asym_fit_%s.root", chan);

    // eta bin centers + half-widths from the charge-asym graph (CM frame)
    TFile *fC = TFile::Open(sCharge, "READ");
    TGraphErrors *gA = (fC && !fC->IsZombie()) ? (TGraphErrors *)fC->Get("g_chargeAsym_mt") : nullptr;
    if (!gA) { std::cerr << "[ERROR] need g_chargeAsym_mt in " << sCharge
                         << " for the CM eta bin positions (run charge_asym.C on the fitted yields first).\n"; return; }
    const int NY = gA->GetN();
    std::vector<double> xc(NY), exc(NY);
    for (int i = 0; i < NY; ++i) { xc[i] = gA->GetPointX(i); exc[i] = gA->GetErrorX(i); }

    double yWp[12], eWp[12], yWm[12], eWm[12];
    if (!readBinYields(sCsv, "Wp", "lab", yWp, eWp) || !readBinYields(sCsv, "Wm", "lab", yWm, eWm))
    { std::cerr << "[ERROR] could not read per-bin yields from " << sCsv << "\n"; return; }

    std::vector<double> sP(NY), esP(NY), sM(NY), esM(NY), sI(NY), esI(NY);
    double ymax = 0.0;
    for (int i = 0; i < NY && i < 12; ++i)
    {
        double dEta = 2.0 * exc[i];
        if (dEta <= 0) dEta = 0.4;
        sP[i] = yWp[i] / (Lint * dEta);  esP[i] = eWp[i] / (Lint * dEta);
        sM[i] = yWm[i] / (Lint * dEta);  esM[i] = eWm[i] / (Lint * dEta);
        const double yi = yWp[i] + yWm[i];
        const double ei = std::sqrt(eWp[i] * eWp[i] + eWm[i] * eWm[i]);
        sI[i] = yi / (Lint * dEta);      esI[i] = ei / (Lint * dEta);
        ymax = std::max(ymax, std::max(sI[i] + esI[i], std::max(sP[i] + esP[i], sM[i] + esM[i])));
    }

    TGraphErrors *gI = new TGraphErrors(NY, xc.data(), sI.data(), exc.data(), esI.data());
    TGraphErrors *gP = new TGraphErrors(NY, xc.data(), sP.data(), exc.data(), esP.data());
    TGraphErrors *gM = new TGraphErrors(NY, xc.data(), sM.data(), exc.data(), esM.data());

    PlotStyle ps; ps.logy = false;
    TCanvas *c = new TCanvas(Form("c_dsig_%s", chan), "", ps.w, ps.h);
    ApplyCanvasStyle(c, ps);
    c->cd();

    const double xlo = xc[0] - exc[0] - 0.1, xhi = xc[NY - 1] + exc[NY - 1] + 0.1;
    TH1F *fr = new TH1F(Form("fr_dsig_%s", chan), "", 100, xlo, xhi);
    fr->SetStats(0);
    ApplyHistStyle(fr, ps, Form("#eta^{%s}_{CM}", lepSym), "d#sigma^{fid}_{W}/d#eta (nb)");
    fr->SetMinimum(0.0);
    fr->SetMaximum(1.55 * ymax);
    fr->Draw();

    gI->SetMarkerStyle(20); gI->SetMarkerSize(1.3); gI->SetMarkerColor(kBlack);    gI->SetLineColor(kBlack);    gI->SetLineWidth(2);
    gP->SetMarkerStyle(21); gP->SetMarkerSize(1.3); gP->SetMarkerColor(kAzure + 2); gP->SetLineColor(kAzure + 2); gP->SetLineWidth(2);
    gM->SetMarkerStyle(22); gM->SetMarkerSize(1.4); gM->SetMarkerColor(kRed + 1);   gM->SetLineColor(kRed + 1);   gM->SetLineWidth(2);
    gI->Draw("P SAME"); gP->Draw("P SAME"); gM->Draw("P SAME");

    DrawHeader(ps, "", Form("W #rightarrow %s #nu", lepSym), "fiducial cross section");

    TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.030);
    tx.DrawLatex(0.18, 0.235, "stat. unc. only"); // L is already in the CMS_lumi label
    tx.DrawLatex(0.18, 0.190, "no eff and acc. correction");

    TLegend *leg = new TLegend(0.18, 0.61, 0.40, 0.79); // upper-left, below CMS (header is upper-right)
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.036);
    leg->AddEntry(gI, "W (incl.)", "lep");
    leg->AddEntry(gP, "W^{+}", "lep");
    leg->AddEntry(gM, "W^{-}", "lep");
    leg->Draw();

    CMS_lumi(c, 13, 10);
    c->RedrawAxis(); c->Modified(); c->Update();

    const std::string outDir = "./plots/xsec";
    gSystem->mkdir(outDir.c_str(), kTRUE);
    c->SaveAs(Form("%s/W_dsigma_deta_%s.png", outDir.c_str(), chan));
    c->SaveAs(Form("%s/W_dsigma_deta_%s.pdf", outDir.c_str(), chan));
    std::cout << "[OK] Saved " << outDir << "/W_dsigma_deta_" << chan << ".{png,pdf}  (L = " << Lint << " nb^-1)\n";
}
