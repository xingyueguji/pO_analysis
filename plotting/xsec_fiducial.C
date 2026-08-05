#include "TCanvas.h"
#include "TH1F.h"
#include "TH2D.h"
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
#include "disc_variants.h"             // pODisc::Spec/GraphFile (W-discriminant tags)

// =============================================================================
// xsec_fiducial.C -- fiducial W cross sections (W+, W-, W inclusive) from the
// fitted signal yields in the Combine summary CSVs, muon and electron overlaid.
//
//   sigma_fid = N_fit / L_int        (N_fit = fitted signal yield from the CSV)
//
// `disc` = met|leppt|leppt_mt40 selects WHICH fit's yields are used
// (2026-08-03): it drives the default CSV / charge-asym paths (out-tree
// pO_fit_out<suffix>/) and the output folder plots/xsec/<disc>/, and is
// stamped into the plot header. Runner: analysis/run_observables.sh.
//
// NOTE: this is the fiducial cross section BEFORE the lepton efficiency
// correction (eff/reco/ID/iso/trigger not yet applied), i.e. sigma_fid x eps.
// That is exactly why muon and electron disagree here -- the gap is the channel
// efficiency difference; once eps is applied they should converge (universality).
// Stat (fit) uncertainty only; lumi uncertainty not included.
//
//   root -l -q 'xsec_fiducial.C+'              // met; +("leppt") etc. for variants
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
    const char *disc   = "met",    // met|leppt|leppt_mt40 (drives defaults + out folder)
    const char *muCsv  = nullptr,  // default: pO_fit_out<suffix>/mu/summary/mu_summary.csv
    const char *eleCsv = nullptr,  // default: pO_fit_out<suffix>/ele/summary/ele_summary.csv
    double Lint = pONorm::kLumi_invnb /* 46.5 nb^-1 */)
{
    gStyle->SetEndErrorSize(4);

    TString dsuf, discLabel;
    if (!pODisc::Spec(disc, dsuf, discLabel)) return;
    const TString sMuCsv = muCsv ? TString(muCsv)
        : TString::Format("../../HiggsAnalysis-CombinedLimit/test/pO_fit_out%s/mu/summary/mu_summary.csv", dsuf.Data());
    const TString sEleCsv = eleCsv ? TString(eleCsv)
        : TString::Format("../../HiggsAnalysis-CombinedLimit/test/pO_fit_out%s/ele/summary/ele_summary.csv", dsuf.Data());

    const char *regions[3] = {"Wp_incl", "Wm_incl", "W_incl"};
    const char *xlab[3]    = {"W^{+}", "W^{-}", "W"};

    double yMu[3] = {0, 0, 0}, eMu[3] = {0, 0, 0};
    double yEl[3] = {0, 0, 0}, eEl[3] = {0, 0, 0};
    bool haveMu = false, haveEl = false;
    for (int i = 0; i < 3; ++i)
    {
        double y, e;
        if (readYield(sMuCsv.Data(), regions[i], y, e))  { yMu[i] = y / Lint; eMu[i] = e / Lint; haveMu = true; }
        if (readYield(sEleCsv.Data(), regions[i], y, e)) { yEl[i] = y / Lint; eEl[i] = e / Lint; haveEl = true; }
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
    // Discriminant tag as a short 4th header line (long labels overflow sub2)
    TLatex ltag; ltag.SetNDC(); ltag.SetTextFont(ps.font); ltag.SetTextAlign(13);
    ltag.SetTextSize(ps.boxTextSize);
    ltag.DrawLatex(ps.headerX, ps.headerY - 3.0 * ps.headerDy,
                   Form("%s fit", discLabel.Data()));

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

    // Per-discriminant output folder so met / leppt / leppt_mt40 coexist.
    const std::string outDir = std::string("./plots/xsec/") + disc;
    gSystem->mkdir(outDir.c_str(), kTRUE);
    c->SaveAs((outDir + "/W_fiducial.png").c_str());
    c->SaveAs((outDir + "/W_fiducial.pdf").c_str());
    std::cout << "[OK] Saved " << outDir << "/W_fiducial.{png,pdf}  (disc=" << disc
              << ", L = " << Lint << " nb^-1)\n";
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
// taken from g_chargeAsym (charge_asym_fit_<chan>_<disc>.root; legacy alias
// g_chargeAsym_mt accepted) so the axis matches the charge-asymmetry plot
// exactly (lab bins shifted to the CM frame). Same caveats as xsec_fiducial:
// sigma_fid x eps (no efficiency correction), stat unc only.
//
//   root -l -q 'xsec_fiducial.C+' -e 'xsec_fiducial_diff(false)'   // muon; (true)=electron
//   (2nd arg = disc, e.g. xsec_fiducial_diff(false, "leppt"))
// =============================================================================
void xsec_fiducial_diff(bool isElec = false,
                        const char *disc = "met",
                        const char *csv = nullptr,
                        const char *chargeFile = nullptr,
                        double Lint = pONorm::kLumi_invnb)
{
    gStyle->SetEndErrorSize(3);
    const char *chan = isElec ? "ele" : "mu";
    const char *lepSym = isElec ? "e" : "#mu";

    TString dsuf, discLabel;
    if (!pODisc::Spec(disc, dsuf, discLabel)) return;
    const TString sCsv = csv ? TString(csv)
        : TString::Format("../../HiggsAnalysis-CombinedLimit/test/pO_fit_out%s/%s/summary/%s_W_yields.csv", dsuf.Data(), chan, chan);
    const TString sCharge = chargeFile ? TString(chargeFile)
        : pODisc::GraphFile("charge_asym", chan, disc);

    // eta bin centers + half-widths from the charge-asym graph (CM frame);
    // primary name first, legacy "_mt" alias as fallback (see naming audit).
    TFile *fC = TFile::Open(sCharge, "READ");
    TGraphErrors *gA = (fC && !fC->IsZombie()) ? (TGraphErrors *)fC->Get("g_chargeAsym") : nullptr;
    if (!gA && fC && !fC->IsZombie()) gA = (TGraphErrors *)fC->Get("g_chargeAsym_mt");
    if (!gA) { std::cerr << "[ERROR] need g_chargeAsym(_mt) in " << sCharge
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
    // Discriminant tag as a short 4th header line (long labels overflow sub2)
    TLatex ltag; ltag.SetNDC(); ltag.SetTextFont(ps.font); ltag.SetTextAlign(13);
    ltag.SetTextSize(ps.boxTextSize);
    ltag.DrawLatex(ps.headerX, ps.headerY - 3.0 * ps.headerDy,
                   Form("%s fit", discLabel.Data()));

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

    // Per-discriminant output folder so met / leppt / leppt_mt40 coexist.
    const std::string outDir = std::string("./plots/xsec/") + disc;
    gSystem->mkdir(outDir.c_str(), kTRUE);
    c->SaveAs(Form("%s/W_dsigma_deta_%s.png", outDir.c_str(), chan));
    c->SaveAs(Form("%s/W_dsigma_deta_%s.pdf", outDir.c_str(), chan));
    std::cout << "[OK] Saved " << outDir << "/W_dsigma_deta_" << chan << ".{png,pdf}  (disc="
              << disc << ", L = " << Lint << " nb^-1)\n";
}

// -----------------------------------------------------------------------------
// per-bin PREFIT MC signal (mu+e summed template integrals, col 6
// signal_prefit of comb_W_yields.csv) -- the fit's r = 1 expectation; the
// postfit/prefit ratio per bin IS the fitted r_<C>_y<i>.
static bool readBinPrefit(const char *csv, const char *charge, const char *binning,
                          double s[12])
{
    for (int i = 0; i < 12; ++i) s[i] = 0;
    std::ifstream in(csv);
    if (!in) return false;
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
            if (iy >= 0 && iy < 12) { s[iy] = std::atof(c[6].c_str()); ++found; }
        }
    }
    return found > 0;
}

// -----------------------------------------------------------------------------
// comb_summary.csv reader: fit,param,value,error,...  (extract_pO_simfit.C
// writes the covariance-propagated inclusive sums as e.g.
// "simfit_lab,W_sum_yield,<yield>,<err>,,,").
static bool readCombSum(const char *csv, const char *fit, const char *param,
                        double &y, double &e)
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
        if (c.size() >= 4 && c[0] == fit && c[1] == param)
        {
            y = std::atof(c[2].c_str());
            e = std::atof(c[3].c_str());
            return true;
        }
    }
    std::cerr << "[WARN] row '" << fit << "," << param << "' not found in " << csv << "\n";
    return false;
}

// =============================================================================
// xsec_fiducial_comb -- the simfit (grand simultaneous fit) analog (2026-08-05):
// mu+e-SUMMED fiducial cross sections from the comb summaries.
//   inclusive   : sigma = Y/L with Y = the COVARIANCE-propagated Wp/Wm/W sums
//                 written by extract_pO_simfit.C (simfit_lab rows).
//   differential: d(sigma)/d(eta_CM) from comb_W_yields.csv (lab bins); the
//                 per-bin W-inclusive error includes cov(Wp_yi, Wm_yi) from
//                 h_cov_yield in comb_fitted_yields.root (quadrature fallback).
// NB with the mu/e-SHARED r's this is sigma_fid(mu) + sigma_fid(e), still
// x eps (no eff/acc correction) and stat-only. The per-flavour mu-vs-e
// comparison plot remains a LEGACY-fit product by construction.
// Outputs: plots/comb/xsec/<disc>/{W_fiducial, W_dsigma_deta_comb}.{png,pdf}
//
//   root -l -q -e 'gROOT->LoadMacro("xsec_fiducial.C+"); xsec_fiducial_comb("met");'
// =============================================================================
void xsec_fiducial_comb(const char *disc = "met",
                        const char *summaryCsv = nullptr, // default: simfit/summary/comb_summary.csv
                        const char *binsCsv    = nullptr, // default: simfit/summary/comb_W_yields.csv
                        const char *chargeFile = nullptr, // default: charge_asym_fit_comb_<disc>.root
                        const char *yieldsRoot = nullptr, // default: simfit/summary/comb_fitted_yields.root
                        const char *genFile    = "../skim/rootfile/gen_xsec.root", // skim/gen_xsec.C output (disc-independent); missing -> overlay skipped
                        double Lint = pONorm::kLumi_invnb)
{
    gStyle->SetEndErrorSize(4);

    TString dsuf, discLabel;
    if (!pODisc::Spec(disc, dsuf, discLabel)) return;
    const TString base = TString::Format(
        "../../HiggsAnalysis-CombinedLimit/test/pO_fit_out%s/simfit/summary", dsuf.Data());
    const TString sSum    = summaryCsv ? TString(summaryCsv) : base + "/comb_summary.csv";
    const TString sBins   = binsCsv    ? TString(binsCsv)    : base + "/comb_W_yields.csv";
    const TString sRoot   = yieldsRoot ? TString(yieldsRoot) : base + "/comb_fitted_yields.root";
    const TString sCharge = chargeFile ? TString(chargeFile) : pODisc::GraphFile("charge_asym", "comb", disc);

    const std::string outDir = std::string("./plots/comb/xsec/") + disc;
    gSystem->mkdir(outDir.c_str(), kTRUE);

    // ---- 1) inclusive W+/W-/W (covariance-propagated sums from the lab fit) --
    const char *params[3] = {"Wp_sum_yield", "Wm_sum_yield", "W_sum_yield"};
    const char *xlab[3]   = {"W^{+}", "W^{-}", "W"};
    double yv[3] = {0, 0, 0}, ev[3] = {0, 0, 0};
    bool haveIncl = true;
    for (int i = 0; i < 3; ++i)
    {
        double y, e;
        if (readCombSum(sSum.Data(), "simfit_lab", params[i], y, e)) { yv[i] = y / Lint; ev[i] = e / Lint; }
        else haveIncl = false;
        printf("[xsec-comb] %-13s  mu+e = %6.2f +/- %5.2f nb\n", params[i], yv[i], ev[i]);
    }
    // GENERATOR level (skim/gen_xsec.C): per-flavour sigma_i per lab bin, drawn
    // x2 = mu+e summed so every series on these plots shares the same
    // convention. gen/reco gap = combined acceptance x efficiency.
    TFile *fG = TFile::Open(genFile, "READ");
    TH1D *hGp = (fG && !fG->IsZombie()) ? (TH1D *)fG->Get("h_gen_sig_Wp") : nullptr;
    TH1D *hGm = (fG && !fG->IsZombie()) ? (TH1D *)fG->Get("h_gen_sig_Wm") : nullptr;
    const bool haveGen = (hGp && hGm);
    double genTot[3] = {0, 0, 0};
    if (haveGen)
    {
        genTot[0] = 2.0 * hGp->Integral();
        genTot[1] = 2.0 * hGm->Integral();
        genTot[2] = genTot[0] + genTot[1];
        printf("[xsec-comb] gen level    mu+e = %6.2f / %6.2f / %6.2f nb (W+/W-/W, |eta_lab|<2.4)\n",
               genTot[0], genTot[1], genTot[2]);
    }
    else
        std::cerr << "[WARN] no usable " << genFile
                  << " -> generator-level overlay skipped (run skim/gen_xsec.C once)\n";

    // prefit MC expectation (r = 1): the mu+e signal-template integrals summed
    // over the lab bins -- the same S the fitted yield Y = r*S scales.
    double sPre[3] = {0, 0, 0};
    bool havePre = false;
    {
        double sp[12], sm[12];
        if (readBinPrefit(sBins.Data(), "Wp", "lab", sp) && readBinPrefit(sBins.Data(), "Wm", "lab", sm))
        {
            havePre = true;
            for (int i = 0; i < 12; ++i) { sPre[0] += sp[i]; sPre[1] += sm[i]; }
            sPre[2] = sPre[0] + sPre[1];
            for (int i = 0; i < 3; ++i) sPre[i] /= Lint;
            printf("[xsec-comb] prefit MC   mu+e = %6.2f / %6.2f / %6.2f nb (W+/W-/W)\n",
                   sPre[0], sPre[1], sPre[2]);
        }
    }

    if (haveIncl)
    {
        double x[3] = {1, 2, 3}, ex0[3] = {0, 0, 0};
        double xp[3] = {1.14, 2.14, 3.14}, ey0[3] = {0, 0, 0};
        double xg[3] = {0.86, 1.86, 2.86};
        TGraphErrors *g = new TGraphErrors(3, x, yv, ex0, ev);
        TGraphErrors *gPre = havePre ? new TGraphErrors(3, xp, sPre, ex0, ey0) : nullptr;
        TGraphErrors *gGen = haveGen ? new TGraphErrors(3, xg, genTot, ex0, ey0) : nullptr;

        PlotStyle ps; ps.logy = false;
        TCanvas *c = new TCanvas("c_xsec_fid_comb", "", ps.w, ps.h);
        ApplyCanvasStyle(c, ps);
        c->cd();
        double ymax = 0.0;
        for (int i = 0; i < 3; ++i)
        {
            ymax = std::max(ymax, yv[i] + ev[i]);
            if (havePre) ymax = std::max(ymax, sPre[i]);
            if (haveGen) ymax = std::max(ymax, genTot[i]);
        }
        TH1F *fr = new TH1F("fr_xsec_comb", "", 3, 0.5, 3.5);
        for (int i = 0; i < 3; ++i) fr->GetXaxis()->SetBinLabel(i + 1, xlab[i]);
        fr->SetStats(0);
        ApplyHistStyle(fr, ps, "", "#sigma^{fid}_{W} (nb)");
        fr->GetXaxis()->SetLabelSize(0.065);
        fr->GetXaxis()->SetNdivisions(3);
        fr->SetMinimum(0.0);
        fr->SetMaximum(1.45 * ymax);
        fr->Draw();

        if (gGen)
        {
            gGen->SetMarkerStyle(27); gGen->SetMarkerSize(1.9); // open diamond = gen
            gGen->SetMarkerColor(kGreen + 2); gGen->SetLineColor(kGreen + 2); gGen->SetLineWidth(2);
            gGen->Draw("P SAME");
        }
        if (gPre)
        {
            gPre->SetMarkerStyle(24); gPre->SetMarkerSize(1.6); // open circle = reco MC
            gPre->SetMarkerColor(kGray + 2); gPre->SetLineColor(kGray + 2); gPre->SetLineWidth(2);
            gPre->Draw("P SAME");
        }
        g->SetMarkerStyle(20); g->SetMarkerSize(1.5);
        g->SetMarkerColor(kBlack); g->SetLineColor(kBlack); g->SetLineWidth(2);
        g->Draw("P SAME");

        DrawHeader(ps, "", "W #rightarrow l #nu", "fiducial cross section");
        TLatex ltag; ltag.SetNDC(); ltag.SetTextFont(ps.font); ltag.SetTextAlign(13);
        ltag.SetTextSize(ps.boxTextSize);
        ltag.DrawLatex(ps.headerX, ps.headerY - 3.0 * ps.headerDy, Form("%s fit", discLabel.Data()));

        TLegend *leg = new TLegend(0.20, 0.60, 0.56, 0.80);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.032);
        leg->AddEntry(g, "post-fit (#mu + e)", "lep");
        if (gPre) leg->AddEntry(gPre, "MC reco expectation (prefit)", "p");
        if (gGen) leg->AddEntry(gGen, "MC generator level (|#eta^{l}_{lab}| < 2.4)", "p");
        leg->Draw();

        TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.030);
        tx.DrawLatex(0.18, 0.30, "#mu + e combined (simfit)");
        tx.DrawLatex(0.18, 0.255, "#sigma = [N(#mu) + N(e)] / L, stat. unc. only");
        tx.DrawLatex(0.18, 0.21, "no eff and acc. correction");

        CMS_lumi(c, 13, 10);
        c->RedrawAxis(); c->Modified(); c->Update();
        c->SaveAs((outDir + "/W_fiducial.png").c_str());
        c->SaveAs((outDir + "/W_fiducial.pdf").c_str());
        std::cout << "[OK] Saved " << outDir << "/W_fiducial.{png,pdf}  (disc=" << disc
                  << ", L = " << Lint << " nb^-1)\n";
    }

    // ---- 2) differential d(sigma)/d(eta_CM), mu+e summed ---------------------
    TFile *fC = TFile::Open(sCharge, "READ");
    TGraphErrors *gA = (fC && !fC->IsZombie()) ? (TGraphErrors *)fC->Get("g_chargeAsym") : nullptr;
    if (!gA && fC && !fC->IsZombie()) gA = (TGraphErrors *)fC->Get("g_chargeAsym_mt");
    if (!gA) { std::cerr << "[ERROR] need g_chargeAsym(_mt) in " << sCharge
                         << " for the CM eta bin positions (run charge_asym.C on the comb yields first).\n"; return; }
    const int NY = gA->GetN();
    std::vector<double> xc(NY), exc(NY);
    for (int i = 0; i < NY; ++i) { xc[i] = gA->GetPointX(i); exc[i] = gA->GetErrorX(i); }

    double yWp[12], eWp[12], yWm[12], eWm[12];
    if (!readBinYields(sBins, "Wp", "lab", yWp, eWp) || !readBinYields(sBins, "Wm", "lab", yWm, eWm))
    { std::cerr << "[ERROR] could not read per-bin yields from " << sBins << "\n"; return; }

    // cov(Wp_yi, Wm_yi) for the per-bin W-inclusive error: the simfit covariance
    // (order [Wp_y0..11, Wm_y0..11]); missing/mismatched -> quadrature fallback.
    TFile *fY = TFile::Open(sRoot, "READ");
    TH2D *hcov = (fY && !fY->IsZombie()) ? (TH2D *)fY->Get("h_cov_yield") : nullptr;
    if (hcov && hcov->GetNbinsX() != 2 * NY) hcov = nullptr;
    if (!hcov) std::cerr << "[WARN] no usable h_cov_yield in " << sRoot
                         << " -> W-inclusive errors use quadrature (no cross term)\n";

    std::vector<double> sP(NY), esP(NY), sM(NY), esM(NY), sI(NY), esI(NY);
    double ymax = 0.0;
    for (int i = 0; i < NY && i < 12; ++i)
    {
        double dEta = 2.0 * exc[i];
        if (dEta <= 0) dEta = 0.4;
        sP[i] = yWp[i] / (Lint * dEta);  esP[i] = eWp[i] / (Lint * dEta);
        sM[i] = yWm[i] / (Lint * dEta);  esM[i] = eWm[i] / (Lint * dEta);
        const double cov = hcov ? hcov->GetBinContent(i + 1, NY + i + 1) : 0.0;
        const double yi  = yWp[i] + yWm[i];
        double vi = eWp[i] * eWp[i] + eWm[i] * eWm[i] + 2.0 * cov;
        if (vi < 0) vi = 0;
        sI[i] = yi / (Lint * dEta);      esI[i] = std::sqrt(vi) / (Lint * dEta);
        ymax = std::max(ymax, std::max(sI[i] + esI[i], std::max(sP[i] + esP[i], sM[i] + esM[i])));
    }

    // prefit MC expectation per bin (r = 1): S/(L*deta), drawn as OPEN markers
    // under the post-fit points -- postfit/prefit per bin IS the fitted r.
    double pWp[12], pWm[12];
    const bool havePreD = readBinPrefit(sBins.Data(), "Wp", "lab", pWp) &&
                          readBinPrefit(sBins.Data(), "Wm", "lab", pWm);
    std::vector<double> pP(NY, 0), pM(NY, 0), pI(NY, 0), zeroe(NY, 0.0);
    if (havePreD)
        for (int i = 0; i < NY && i < 12; ++i)
        {
            double dEta = 2.0 * exc[i];
            if (dEta <= 0) dEta = 0.4;
            pP[i] = pWp[i] / (Lint * dEta);
            pM[i] = pWm[i] / (Lint * dEta);
            pI[i] = (pWp[i] + pWm[i]) / (Lint * dEta);
            ymax = std::max(ymax, std::max(pI[i], std::max(pP[i], pM[i])));
        }
    else
        std::cerr << "[WARN] no prefit signal columns in " << sBins << " -> MC overlay skipped\n";

    // generator level per bin (2 x per-flavour sigma_i / dEta = mu+e summed),
    // drawn as dashed lines; gen bins = the same lab edge set as the yields.
    const bool haveGenD = haveGen && hGp->GetNbinsX() == NY && hGm->GetNbinsX() == NY;
    std::vector<double> gP2(NY, 0), gM2(NY, 0), gI2(NY, 0);
    if (haveGenD)
        for (int i = 0; i < NY; ++i)
        {
            double dEta = 2.0 * exc[i];
            if (dEta <= 0) dEta = 0.4;
            gP2[i] = 2.0 * hGp->GetBinContent(i + 1) / dEta;
            gM2[i] = 2.0 * hGm->GetBinContent(i + 1) / dEta;
            gI2[i] = gP2[i] + gM2[i];
            ymax = std::max(ymax, gI2[i]);
        }

    TGraphErrors *gI = new TGraphErrors(NY, xc.data(), sI.data(), exc.data(), esI.data());
    TGraphErrors *gP = new TGraphErrors(NY, xc.data(), sP.data(), exc.data(), esP.data());
    TGraphErrors *gM = new TGraphErrors(NY, xc.data(), sM.data(), exc.data(), esM.data());
    TGraphErrors *gIp = havePreD ? new TGraphErrors(NY, xc.data(), pI.data(), zeroe.data(), zeroe.data()) : nullptr;
    TGraphErrors *gPp = havePreD ? new TGraphErrors(NY, xc.data(), pP.data(), zeroe.data(), zeroe.data()) : nullptr;
    TGraphErrors *gMp = havePreD ? new TGraphErrors(NY, xc.data(), pM.data(), zeroe.data(), zeroe.data()) : nullptr;
    TGraphErrors *gIg = haveGenD ? new TGraphErrors(NY, xc.data(), gI2.data(), zeroe.data(), zeroe.data()) : nullptr;
    TGraphErrors *gPg = haveGenD ? new TGraphErrors(NY, xc.data(), gP2.data(), zeroe.data(), zeroe.data()) : nullptr;
    TGraphErrors *gMg = haveGenD ? new TGraphErrors(NY, xc.data(), gM2.data(), zeroe.data(), zeroe.data()) : nullptr;

    PlotStyle ps; ps.logy = false;
    TCanvas *c = new TCanvas("c_dsig_comb", "", ps.w, ps.h);
    ApplyCanvasStyle(c, ps);
    c->cd();

    const double xlo = xc[0] - exc[0] - 0.1, xhi = xc[NY - 1] + exc[NY - 1] + 0.1;
    TH1F *fr = new TH1F("fr_dsig_comb", "", 100, xlo, xhi);
    fr->SetStats(0);
    ApplyHistStyle(fr, ps, "#eta^{l}_{CM}", "d#sigma^{fid}_{W}/d#eta (nb)");
    fr->SetMinimum(0.0);
    fr->SetMaximum(1.55 * ymax);
    fr->Draw();

    // gen-level dashed lines first (bottom layer)...
    if (haveGenD)
    {
        gIg->SetLineStyle(2); gIg->SetLineWidth(3); gIg->SetLineColor(kBlack);
        gPg->SetLineStyle(2); gPg->SetLineWidth(3); gPg->SetLineColor(kAzure + 2);
        gMg->SetLineStyle(2); gMg->SetLineWidth(3); gMg->SetLineColor(kRed + 1);
        gIg->Draw("L SAME"); gPg->Draw("L SAME"); gMg->Draw("L SAME");
    }
    // ...then open MC markers (underneath the fit), slightly larger so the ring
    // stays visible when the fitted point sits on top of the expectation (r ~ 1)
    if (havePreD)
    {
        gIp->SetMarkerStyle(24); gIp->SetMarkerSize(1.7); gIp->SetMarkerColor(kBlack);     gIp->SetLineColor(kBlack);
        gPp->SetMarkerStyle(25); gPp->SetMarkerSize(1.7); gPp->SetMarkerColor(kAzure + 2); gPp->SetLineColor(kAzure + 2);
        gMp->SetMarkerStyle(26); gMp->SetMarkerSize(1.8); gMp->SetMarkerColor(kRed + 1);   gMp->SetLineColor(kRed + 1);
        gIp->Draw("P SAME"); gPp->Draw("P SAME"); gMp->Draw("P SAME");
    }
    gI->SetMarkerStyle(20); gI->SetMarkerSize(1.3); gI->SetMarkerColor(kBlack);     gI->SetLineColor(kBlack);     gI->SetLineWidth(2);
    gP->SetMarkerStyle(21); gP->SetMarkerSize(1.3); gP->SetMarkerColor(kAzure + 2); gP->SetLineColor(kAzure + 2); gP->SetLineWidth(2);
    gM->SetMarkerStyle(22); gM->SetMarkerSize(1.4); gM->SetMarkerColor(kRed + 1);   gM->SetLineColor(kRed + 1);   gM->SetLineWidth(2);
    gI->Draw("P SAME"); gP->Draw("P SAME"); gM->Draw("P SAME");

    DrawHeader(ps, "", "W #rightarrow l #nu", "fiducial cross section");
    TLatex ltag2; ltag2.SetNDC(); ltag2.SetTextFont(ps.font); ltag2.SetTextAlign(13);
    ltag2.SetTextSize(ps.boxTextSize);
    ltag2.DrawLatex(ps.headerX, ps.headerY - 3.0 * ps.headerDy, Form("%s fit", discLabel.Data()));

    TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.030);
    tx.DrawLatex(0.18, 0.28, "#mu + e combined (simfit)");
    tx.DrawLatex(0.18, 0.235, "stat. unc. only");
    tx.DrawLatex(0.18, 0.190, "no eff and acc. correction");

    TLegend *leg = new TLegend(0.18, 0.52, 0.48, 0.79);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.031);
    leg->AddEntry(gI, "W (incl.) post-fit", "lep");
    leg->AddEntry(gP, "W^{+} post-fit", "lep");
    leg->AddEntry(gM, "W^{-} post-fit", "lep");
    if (havePreD) leg->AddEntry(gIp, "MC reco expectation (open)", "p");
    if (haveGenD) leg->AddEntry(gIg, "MC generator level (dashed)", "l");
    leg->Draw();

    CMS_lumi(c, 13, 10);
    c->RedrawAxis(); c->Modified(); c->Update();
    c->SaveAs((outDir + "/W_dsigma_deta_comb.png").c_str());
    c->SaveAs((outDir + "/W_dsigma_deta_comb.pdf").c_str());
    std::cout << "[OK] Saved " << outDir << "/W_dsigma_deta_comb.{png,pdf}  (disc="
              << disc << ", L = " << Lint << " nb^-1)\n";

    if (fC) { fC->Close(); delete fC; }
    if (fY) { fY->Close(); delete fY; }
    if (fG) { fG->Close(); delete fG; }
}
