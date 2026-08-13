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

// -----------------------------------------------------------------------------
// per-bin fitted signal strength r +/- rErr (cols r/rErr of comb_W_yields.csv)
// -- the POIs of the grand simultaneous fit; sigma_meas,i = r_i x sigma_gen-fid,i.
static bool readBinR(const char *csv, const char *charge, const char *binning,
                     double r[12], double e[12])
{
    for (int i = 0; i < 12; ++i) { r[i] = 0; e[i] = 0; }
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
            if (iy >= 0 && iy < 12) { r[iy] = std::atof(c[4].c_str()); e[iy] = std::atof(c[5].c_str()); ++found; }
        }
    }
    return found > 0;
}

// -----------------------------------------------------------------------------
// Ingredients shared by every simfit-based cross section, stacked in the
// h_cov_yield order [Wp_y0..11, Wm_y0..11]:
//   G    = gen FIDUCIAL sigma per bin, PER FLAVOUR (nb)  <- skim/gen_xsec.C
//   S    = prefit signal template integral, mu+e summed  <- comb_W_yields.csv
//   R,ER = fitted signal strength r and its error        <- comb_W_yields.csv
//   covr = r covariance = cov_Y_ij/(S_i S_j)             <- h_cov_yield
struct CombIngredients
{
    double G[24], S[24], R[24], ER[24], covr[24][24];
};

static bool loadCombIngredients(const char *binsCsv, const char *yieldsRoot,
                                const char *genFile, CombIngredients &in)
{
    TFile *fG = TFile::Open(genFile, "READ");
    TH1D *hGp = (fG && !fG->IsZombie()) ? (TH1D *)fG->Get("h_gen_sig_Wp") : nullptr;
    TH1D *hGm = (fG && !fG->IsZombie()) ? (TH1D *)fG->Get("h_gen_sig_Wm") : nullptr;
    if (!hGp || !hGm || hGp->GetNbinsX() != 12 || hGm->GetNbinsX() != 12)
    {
        std::cerr << "[ERROR] no usable 12-bin h_gen_sig_{Wp,Wm} in " << genFile
                  << " -- the r x sigma_gen extraction needs it (run skim/gen_xsec.C once).\n";
        if (fG) { fG->Close(); delete fG; }
        return false;
    }

    double rP[12], erP[12], rM[12], erM[12], sp[12], sm[12];
    if (!readBinR(binsCsv, "Wp", "lab", rP, erP) || !readBinR(binsCsv, "Wm", "lab", rM, erM))
    {
        std::cerr << "[ERROR] could not read per-bin r from " << binsCsv << "\n";
        fG->Close(); delete fG; return false;
    }
    const bool haveS = readBinPrefit(binsCsv, "Wp", "lab", sp) &&
                       readBinPrefit(binsCsv, "Wm", "lab", sm);
    for (int i = 0; i < 12; ++i)
    {
        in.G[i] = hGp->GetBinContent(i + 1);  in.G[12 + i] = hGm->GetBinContent(i + 1);
        in.S[i] = haveS ? sp[i] : 0;          in.S[12 + i] = haveS ? sm[i] : 0;
        in.R[i] = rP[i]; in.ER[i] = erP[i];   in.R[12 + i] = rM[i]; in.ER[12 + i] = erM[i];
    }
    fG->Close(); delete fG;

    // r covariance from the yield covariance; fallback = diagonal rErr^2
    // (sums then lose the cross terms).
    TFile *fY = TFile::Open(yieldsRoot, "READ");
    TH2D *hcov = (fY && !fY->IsZombie()) ? (TH2D *)fY->Get("h_cov_yield") : nullptr;
    if (hcov && hcov->GetNbinsX() != 24) hcov = nullptr;
    if (!hcov || !haveS)
    {
        hcov = nullptr;
        std::cerr << "[WARN] no usable h_cov_yield in " << yieldsRoot
                  << (haveS ? "" : " (and/or no signal_prefit columns)")
                  << " -> r-covariance falls back to DIAGONAL rErr^2\n";
    }
    for (int i = 0; i < 24; ++i)
        for (int j = 0; j < 24; ++j)
            in.covr[i][j] = (hcov && in.S[i] > 0 && in.S[j] > 0)
                                ? hcov->GetBinContent(i + 1, j + 1) / (in.S[i] * in.S[j])
                                : (i == j ? in.ER[i] * in.ER[i] : 0.0);
    if (fY) { fY->Close(); delete fY; }
    return true;
}

// Sum_i w_i r_i over bins [lo, hi) with Var = Sum_ij w_i w_j cov_r,ij -- the
// one place the r-covariance is propagated (w = sigma_gen for the measurement,
// w = per-flavour prefit S for the count-based version).
static void sumWithCov(const CombIngredients &in, const double *w, int lo, int hi,
                       double &val, double &err)
{
    val = 0;
    double var = 0;
    for (int i = lo; i < hi; ++i)
    {
        val += w[i] * in.R[i];
        for (int j = lo; j < hi; ++j) var += w[i] * w[j] * in.covr[i][j];
    }
    err = std::sqrt(std::max(0.0, var));
}

// =============================================================================
// xsec_fiducial_comb -- the simfit (grand simultaneous fit) analog (2026-08-05;
// r x sigma_gen EXTRACTION since 2026-08-12):
//
//     sigma_meas,i = r_i x sigma_gen-fid,i          [PER LEPTON FLAVOUR, nb]
//
// r_i = the mu/e-shared per-(charge, y-bin) POI (r/rErr of comb_W_yields.csv),
// sigma_gen-fid,i = the gen FIDUCIAL sigma from skim/gen_xsec.C (h_gen_sig_*:
// bare post-FSR lepton pT > 25, |eta_lab| < 2.4). Algebraically identical to
// N_fit/(L*(A*eps)_MC), but kA_O/kSigma_* cancel between r and sigma_gen and
// future template SFs propagate through r automatically. Inclusive sums
// Sum_i r_i*sigma_i carry the full r-covariance, recovered from h_cov_yield
// (yield space) via cov(r_i,r_j) = cov_Y_ij/(S_i*S_j), S = signal_prefit
// (diagonal-rErr fallback + WARN). Gen MC-stat error (~0.3%/bin) is neglected
// against the fit error (~6%/bin).
//
// PER-FLAVOUR convention since 2026-08-12: the result is sigma(W -> l nu) for
// ONE lepton flavour (lepton universality; the r's are mu/e-shared) -- no more
// "x2 = mu+e summed" display.
//
// The COUNT record is untouched: comb_W_yields.csv / comb_summary.csv keep the
// fitted yields (printed here as a console cross-check), and charge_asym /
// FBratio stay count-based by design. This macro additionally writes the sigma
// table plots/comb/xsec/<disc>/xsec_comb.csv (per-bin + inclusive).
// Outputs: plots/comb/xsec/<disc>/{W_fiducial, W_dsigma_deta_comb}.{png,pdf}
//          + xsec_comb.csv
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

    // ---- 0) ingredients for sigma = r x sigma_gen-fid ------------------------
    CombIngredients in;
    if (!loadCombIngredients(sBins.Data(), sRoot.Data(), genFile, in)) return;
    const double *G = in.G, *R = in.R, *ER = in.ER;
    auto &covr = in.covr;

    // ---- 1) inclusive W+/W-/W: sigma = Sum_i r_i sigma_gen,i (cov-propagated)
    const char *xlab[3] = {"W^{+}", "W^{-}", "W"};
    double yv[3] = {0, 0, 0}, ev[3] = {0, 0, 0}, genTot[3] = {0, 0, 0};
    const int klo[3] = {0, 12, 0}, khi[3] = {12, 24, 24};
    for (int k = 0; k < 3; ++k)
    {
        sumWithCov(in, G, klo[k], khi[k], yv[k], ev[k]);
        for (int i = klo[k]; i < khi[k]; ++i) genTot[k] += G[i];
        printf("[xsec-comb] %-6s measured = %6.2f +/- %4.2f nb   (gen fid %6.2f nb, eff. r = %.3f)  [per flavour]\n",
               xlab[k], yv[k], ev[k], genTot[k], genTot[k] > 0 ? yv[k] / genTot[k] : 0.0);
    }

    // console cross-check ONLY: the COUNT record (stays in the fork CSVs;
    // charge_asym / FBratio remain count-based).
    {
        const char *params[3] = {"Wp_sum_yield", "Wm_sum_yield", "W_sum_yield"};
        for (int i = 0; i < 3; ++i)
        {
            double y, e;
            if (readCombSum(sSum.Data(), "simfit_lab", params[i], y, e))
                printf("[xsec-comb] counts %-13s N/L = %6.2f +/- %4.2f nb (mu+e, uncorrected)\n",
                       params[i], y / Lint, e / Lint);
        }
    }

    {
        double x[3] = {1, 2, 3}, ex0[3] = {0, 0, 0};
        double xg[3] = {0.86, 1.86, 2.86}, ey0[3] = {0, 0, 0};
        TGraphErrors *g = new TGraphErrors(3, x, yv, ex0, ev);
        TGraphErrors *gGen = new TGraphErrors(3, xg, genTot, ex0, ey0);

        PlotStyle ps; ps.logy = false;
        TCanvas *c = new TCanvas("c_xsec_fid_comb", "", ps.w, ps.h);
        ApplyCanvasStyle(c, ps);
        c->cd();
        double ymax = 0.0;
        for (int i = 0; i < 3; ++i)
            ymax = std::max(ymax, std::max(yv[i] + ev[i], genTot[i]));
        TH1F *fr = new TH1F("fr_xsec_comb", "", 3, 0.5, 3.5);
        for (int i = 0; i < 3; ++i) fr->GetXaxis()->SetBinLabel(i + 1, xlab[i]);
        fr->SetStats(0);
        ApplyHistStyle(fr, ps, "", "#sigma^{fid}_{W#rightarrowl#nu} (nb)");
        fr->GetXaxis()->SetLabelSize(0.065);
        fr->GetXaxis()->SetNdivisions(3);
        fr->SetMinimum(0.0);
        fr->SetMaximum(1.45 * ymax);
        fr->Draw();

        gGen->SetMarkerStyle(27); gGen->SetMarkerSize(1.9); // open diamond = gen fid (r = 1)
        gGen->SetMarkerColor(kGreen + 2); gGen->SetLineColor(kGreen + 2); gGen->SetLineWidth(2);
        gGen->Draw("P SAME");
        g->SetMarkerStyle(20); g->SetMarkerSize(1.5);
        g->SetMarkerColor(kBlack); g->SetLineColor(kBlack); g->SetLineWidth(2);
        g->Draw("P SAME");

        DrawHeader(ps, "", "W #rightarrow l #nu", "fiducial cross section");
        TLatex ltag; ltag.SetNDC(); ltag.SetTextFont(ps.font); ltag.SetTextAlign(13);
        ltag.SetTextSize(ps.boxTextSize);
        ltag.DrawLatex(ps.headerX, ps.headerY - 3.0 * ps.headerDy, Form("%s fit", discLabel.Data()));

        TLegend *leg = new TLegend(0.20, 0.56, 0.56, 0.72);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.032);
        leg->AddEntry(g, "measured, #sigma = r #times #sigma^{gen}_{fid}", "lep");
        leg->AddEntry(gGen, "#sigma^{gen}_{fid} (POWHEG, r = 1)", "p");
        leg->Draw();

        TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.030);
        tx.DrawLatex(0.18, 0.30,  "#mu + e combined (simfit), per lepton flavour");
        tx.DrawLatex(0.18, 0.255, "fiducial: bare lepton p_{T} > 25 GeV, |#eta_{lab}| < 2.4");
        tx.DrawLatex(0.18, 0.21,  "#sigma_{i} = r_{i} #times #sigma^{gen}_{fid,i}, summed with r-covariance");
        tx.DrawLatex(0.18, 0.165, "stat. (fit) unc. only");

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

    if (NY != 12)
    { std::cerr << "[ERROR] charge-asym graph has " << NY << " points (expected 12) -- binning mismatch.\n"; return; }

    // measured d(sigma)/d(eta): r_i * sigma_gen,i / dEta (per flavour);
    // W-inclusive per bin carries cov(r_Wp_i, r_Wm_i) from the r-covariance.
    std::vector<double> sP(NY), esP(NY), sM(NY), esM(NY), sI(NY), esI(NY);
    std::vector<double> gP2(NY, 0), gM2(NY, 0), gI2(NY, 0), zeroe(NY, 0.0);
    double ymax = 0.0;
    for (int i = 0; i < NY; ++i)
    {
        double dEta = 2.0 * exc[i];
        if (dEta <= 0) dEta = 0.4;
        const double gp = G[i], gm = G[12 + i];
        sP[i] = R[i] * gp / dEta;       esP[i] = std::sqrt(covr[i][i]) * gp / dEta;
        sM[i] = R[12 + i] * gm / dEta;  esM[i] = std::sqrt(covr[12 + i][12 + i]) * gm / dEta;
        double vi = gp * gp * covr[i][i] + gm * gm * covr[12 + i][12 + i]
                  + 2.0 * gp * gm * covr[i][12 + i];
        if (vi < 0) vi = 0;
        sI[i]  = (R[i] * gp + R[12 + i] * gm) / dEta;
        esI[i] = std::sqrt(vi) / dEta;
        // gen fiducial (r = 1 expectation), dashed lines
        gP2[i] = gp / dEta;
        gM2[i] = gm / dEta;
        gI2[i] = gP2[i] + gM2[i];
        ymax = std::max(ymax, std::max(sI[i] + esI[i], gI2[i]));
    }

    TGraphErrors *gI = new TGraphErrors(NY, xc.data(), sI.data(), exc.data(), esI.data());
    TGraphErrors *gP = new TGraphErrors(NY, xc.data(), sP.data(), exc.data(), esP.data());
    TGraphErrors *gM = new TGraphErrors(NY, xc.data(), sM.data(), exc.data(), esM.data());
    TGraphErrors *gIg = new TGraphErrors(NY, xc.data(), gI2.data(), zeroe.data(), zeroe.data());
    TGraphErrors *gPg = new TGraphErrors(NY, xc.data(), gP2.data(), zeroe.data(), zeroe.data());
    TGraphErrors *gMg = new TGraphErrors(NY, xc.data(), gM2.data(), zeroe.data(), zeroe.data());

    PlotStyle ps; ps.logy = false;
    TCanvas *c = new TCanvas("c_dsig_comb", "", ps.w, ps.h);
    ApplyCanvasStyle(c, ps);
    c->cd();

    const double xlo = xc[0] - exc[0] - 0.1, xhi = xc[NY - 1] + exc[NY - 1] + 0.1;
    TH1F *fr = new TH1F("fr_dsig_comb", "", 100, xlo, xhi);
    fr->SetStats(0);
    ApplyHistStyle(fr, ps, "#eta^{l}_{CM}", "d#sigma^{fid}_{W#rightarrowl#nu}/d#eta (nb)");
    fr->SetMinimum(0.0);
    fr->SetMaximum(1.55 * ymax);
    fr->Draw();

    // gen-fiducial dashed lines first (bottom layer) = the r = 1 expectation;
    // measured point / dashed line per bin IS the fitted r_i.
    gIg->SetLineStyle(2); gIg->SetLineWidth(3); gIg->SetLineColor(kBlack);
    gPg->SetLineStyle(2); gPg->SetLineWidth(3); gPg->SetLineColor(kAzure + 2);
    gMg->SetLineStyle(2); gMg->SetLineWidth(3); gMg->SetLineColor(kRed + 1);
    gIg->Draw("L SAME"); gPg->Draw("L SAME"); gMg->Draw("L SAME");

    gI->SetMarkerStyle(20); gI->SetMarkerSize(1.3); gI->SetMarkerColor(kBlack);     gI->SetLineColor(kBlack);     gI->SetLineWidth(2);
    gP->SetMarkerStyle(21); gP->SetMarkerSize(1.3); gP->SetMarkerColor(kAzure + 2); gP->SetLineColor(kAzure + 2); gP->SetLineWidth(2);
    gM->SetMarkerStyle(22); gM->SetMarkerSize(1.4); gM->SetMarkerColor(kRed + 1);   gM->SetLineColor(kRed + 1);   gM->SetLineWidth(2);
    gI->Draw("P SAME"); gP->Draw("P SAME"); gM->Draw("P SAME");

    DrawHeader(ps, "", "W #rightarrow l #nu", "fiducial cross section");
    TLatex ltag2; ltag2.SetNDC(); ltag2.SetTextFont(ps.font); ltag2.SetTextAlign(13);
    ltag2.SetTextSize(ps.boxTextSize);
    ltag2.DrawLatex(ps.headerX, ps.headerY - 3.0 * ps.headerDy, Form("%s fit", discLabel.Data()));

    TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.030);
    tx.DrawLatex(0.52, 0.26,  "#mu + e combined (simfit), per flavour");
    tx.DrawLatex(0.52, 0.215, "#sigma_{i} = r_{i} #times #sigma^{gen}_{fid,i} (p_{T}^{l} > 25 GeV)");
    tx.DrawLatex(0.52, 0.17,  "stat. (fit) unc. only");

    TLegend *leg = new TLegend(0.18, 0.55, 0.48, 0.79);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.031);
    leg->AddEntry(gI, "W (incl.) measured", "lep");
    leg->AddEntry(gP, "W^{+} measured", "lep");
    leg->AddEntry(gM, "W^{-} measured", "lep");
    leg->AddEntry(gIg, "#sigma^{gen}_{fid} (POWHEG, r = 1, dashed)", "l");
    leg->Draw();

    CMS_lumi(c, 13, 10);
    c->RedrawAxis(); c->Modified(); c->Update();
    c->SaveAs((outDir + "/W_dsigma_deta_comb.png").c_str());
    c->SaveAs((outDir + "/W_dsigma_deta_comb.pdf").c_str());
    std::cout << "[OK] Saved " << outDir << "/W_dsigma_deta_comb.{png,pdf}  (disc="
              << disc << ", L = " << Lint << " nb^-1)\n";

    // ---- 3) sigma table: the r x sigma_gen record (per-bin + inclusive) ------
    // (the COUNT record stays in comb_W_yields.csv / comb_summary.csv, fork.)
    {
        const std::string csvPath = outDir + "/xsec_comb.csv";
        std::ofstream out(csvPath.c_str());
        out << "charge,ybin,etaCM_center,etaCM_halfwidth,r,rErr,sigma_gen_fid_nb,sigma_meas_nb,sigma_meas_err_nb\n";
        const char *cn2[2] = {"Wp", "Wm"};
        for (int q = 0; q < 2; ++q)
            for (int i = 0; i < NY; ++i)
            {
                const int k = 12 * q + i;
                out << Form("%s,%d,%.4f,%.4f,%.5f,%.5f,%.5f,%.5f,%.5f\n",
                            cn2[q], i, xc[i], exc[i], R[k], ER[k],
                            G[k], R[k] * G[k], ER[k] * G[k]);
            }
        const char *cn3[3] = {"Wp", "Wm", "W"};
        for (int k = 0; k < 3; ++k)  // r column = effective (gen-weighted mean) r
            out << Form("%s,incl,,,%.5f,%.5f,%.5f,%.5f,%.5f\n", cn3[k],
                        genTot[k] > 0 ? yv[k] / genTot[k] : 0.0,
                        genTot[k] > 0 ? ev[k] / genTot[k] : 0.0,
                        genTot[k], yv[k], ev[k]);
        std::cout << "[OK] Wrote " << csvPath << "\n";
    }

    if (fC) { fC->Close(); delete fC; }
}

// -----------------------------------------------------------------------------
// PER-FLAVOUR prefit signal per lab bin from THIS repo's structured Combine
// inputs (plotting/mtandmet.C): region dirs W{p,m}_lab_y{i}, absolute `signal`
// template. mu + e reproduce the CSV signal_prefit column exactly (checked).
// Filled into the stacked [Wp_y0..11, Wm_y0..11] order.
static bool readFlavourPrefit(const char *file, double s[24])
{
    TFile *f = TFile::Open(file, "READ");
    if (!f || f->IsZombie())
    { std::cerr << "[WARN] cannot open " << file << "\n"; if (f) delete f; return false; }
    bool ok = true;
    const char *cn[2] = {"Wp", "Wm"};
    for (int q = 0; q < 2 && ok; ++q)
        for (int i = 0; i < 12; ++i)
        {
            TH1 *h = (TH1 *)f->Get(Form("%s_lab_y%d/signal", cn[q], i));
            if (!h) { ok = false; break; }
            s[12 * q + i] = h->Integral();
        }
    if (!ok) std::cerr << "[WARN] missing W{p,m}_lab_y*/signal in " << file << "\n";
    f->Close(); delete f;
    return ok;
}

// =============================================================================
// xsec_fiducial_diag -- PER-FLAVOUR DIAGNOSTIC of the cross-section extraction
// (2026-08-12). One plot per lepton flavour, each with four series per charge:
//
//   sigma^gen        = sigma_gen-fid                  (POWHEG, r = 1)
//   sigma^reco       = S_f / L                        (MC reco expectation, r = 1)
//   r x sigma^gen    = Sum_i r_i sigma_gen,i          (THE measurement -- identical
//                                                      in both plots: r is shared)
//   r x sigma^reco   = Sum_i r_i S_f,i / L = N_fit,f/L (count-based, uncorrected)
//
// The vertical gap within each pair is the SAME factor (A x eps)_MC,f, so the
// plot reads directly as "how big is the MC correction, and does it differ
// between mu and e" -- the acceptance x efficiency the r x sigma_gen extraction
// applies implicitly. Errors on the two r-based series use the full
// r-covariance; the r = 1 series are MC and drawn without errors.
// Sum over flavours of r x sigma^reco reproduces the fork's *_sum_yield/L
// (printed as a cross-check).
//
// The same four series are drawn ETA-BINNED as d(sigma)/d(eta_CM), one canvas
// per (flavour, charge), plus a summary of (A x eps)(eta) for both flavours --
// where an ECAL-crack acceptance hole would show up as an electron-only dip.
//
// Purely diagnostic: the physics result stays in xsec_fiducial_comb.
// Outputs: plots/comb/xsec/<disc>/
//            W_xsec_diag_{mu,ele}.{png,pdf}                    (inclusive)
//            W_dsigma_diag_{mu,ele}_{Wp,Wm,W}.{png,pdf}        (eta-binned)
//            AxEps_vs_eta.{png,pdf}                            (summary)
//            xsec_diag.csv, xsec_diag_bins.csv
//
//   root -l -q -e 'gROOT->LoadMacro("xsec_fiducial.C+"); xsec_fiducial_diag("met");'
// =============================================================================
void xsec_fiducial_diag(const char *disc = "met",
                        const char *binsCsv    = nullptr, // default: simfit/summary/comb_W_yields.csv
                        const char *yieldsRoot = nullptr, // default: simfit/summary/comb_fitted_yields.root
                        const char *chargeFile = nullptr, // default: charge_asym_fit_comb_<disc>.root (eta bin positions)
                        const char *genFile    = "../skim/rootfile/gen_xsec.root",
                        double Lint = pONorm::kLumi_invnb)
{
    gStyle->SetEndErrorSize(4);

    TString dsuf, discLabel;
    if (!pODisc::Spec(disc, dsuf, discLabel)) return;
    const TString base = TString::Format(
        "../../HiggsAnalysis-CombinedLimit/test/pO_fit_out%s/simfit/summary", dsuf.Data());
    const TString sBins = binsCsv    ? TString(binsCsv)    : base + "/comb_W_yields.csv";
    const TString sRoot = yieldsRoot ? TString(yieldsRoot) : base + "/comb_fitted_yields.root";
    const TString sCharge = chargeFile ? TString(chargeFile) : pODisc::GraphFile("charge_asym", "comb", disc);

    const std::string outDir = std::string("./plots/comb/xsec/") + disc;
    gSystem->mkdir(outDir.c_str(), kTRUE);

    CombIngredients in;
    if (!loadCombIngredients(sBins.Data(), sRoot.Data(), genFile, in)) return;

    // eta_CM bin centers/widths for the eta-binned part (same source as
    // xsec_fiducial_comb: the charge-asymmetry graph); missing -> inclusive only.
    TFile *fC = TFile::Open(sCharge, "READ");
    TGraphErrors *gA = (fC && !fC->IsZombie()) ? (TGraphErrors *)fC->Get("g_chargeAsym") : nullptr;
    if (!gA && fC && !fC->IsZombie()) gA = (TGraphErrors *)fC->Get("g_chargeAsym_mt");
    const int NY = (gA && gA->GetN() == 12) ? 12 : 0;
    std::vector<double> xc(12, 0), exc(12, 0);
    if (NY == 12)
        for (int i = 0; i < 12; ++i) { xc[i] = gA->GetPointX(i); exc[i] = gA->GetErrorX(i); }
    else
        std::cerr << "[WARN] no 12-point g_chargeAsym(_mt) in " << sCharge
                  << " -> eta-binned diagnostic skipped (run charge_asym.C on the comb yields first)\n";

    const char *flav[2]   = {"mu", "ele"};
    const char *lepSym[2] = {"#mu", "e"};
    // this repo's Combine inputs: muon in plots/, electron in plots/Elec/
    const TString inFile[2] = {TString::Format("./plots/combine_input_W%s.root", dsuf.Data()),
                               TString::Format("./plots/Elec/combine_input_W%s.root", dsuf.Data())};
    const char *xlab[3]  = {"W^{+}", "W^{-}", "W"};
    const int klo[3] = {0, 12, 0}, khi[3] = {12, 24, 24};

    std::ofstream csv((outDir + "/xsec_diag.csv").c_str());
    csv << "flavour,charge,sigma_gen_nb,sigma_reco_nb,r_x_gen_nb,r_x_gen_err_nb,"
           "r_x_reco_nb,r_x_reco_err_nb,AxEps_reco_over_gen\n";
    std::ofstream csvB((outDir + "/xsec_diag_bins.csv").c_str());
    csvB << "flavour,charge,ybin,etaCM_center,etaCM_halfwidth,dsigma_gen,dsigma_reco,"
            "dsigma_r_gen,dsigma_r_gen_err,dsigma_r_reco,dsigma_r_reco_err,AxEps\n";
    double sumCount[3] = {0, 0, 0}; // mu + e of r x sigma^reco, for the cross-check
    double axeps[2][3][12];         // (A x eps)(eta) per flavour/charge, for the summary plot
    for (int a = 0; a < 2; ++a) for (int b = 0; b < 3; ++b) for (int cc = 0; cc < 12; ++cc) axeps[a][b][cc] = 0;

    // per-bin value with the covariance-correct error, charge k (0=Wp, 1=Wm, 2=W)
    auto binVal = [&](const double *wt, int k, int i, double &v, double &e)
    {
        double tmp[24];
        for (int t = 0; t < 24; ++t) tmp[t] = 0.0;
        if (k != 1) tmp[i]      = wt[i];
        if (k != 0) tmp[12 + i] = wt[12 + i];
        sumWithCov(in, tmp, 0, 24, v, e);
    };

    for (int fl = 0; fl < 2; ++fl)
    {
        double Sf[24], w[24];
        if (!readFlavourPrefit(inFile[fl].Data(), Sf))
        {
            std::cerr << "[ERROR] " << flav[fl] << " diagnostic skipped -- run plotting/mtandmet.C"
                      << (dsuf.Length() ? " for this discriminant" : "") << " first.\n";
            continue;
        }
        for (int i = 0; i < 24; ++i) w[i] = Sf[i] / Lint; // reco sigma per bin (nb)

        double gen[3] = {0, 0, 0}, reco[3] = {0, 0, 0};
        double rg[3], erg[3], rr[3], err[3];
        for (int k = 0; k < 3; ++k)
        {
            for (int i = klo[k]; i < khi[k]; ++i) { gen[k] += in.G[i]; reco[k] += w[i]; }
            sumWithCov(in, in.G, klo[k], khi[k], rg[k], erg[k]);
            sumWithCov(in, w,    klo[k], khi[k], rr[k], err[k]);
            sumCount[k] += rr[k];
            csv << Form("%s,%s,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n",
                        flav[fl], (k == 0 ? "Wp" : k == 1 ? "Wm" : "W"),
                        gen[k], reco[k], rg[k], erg[k], rr[k], err[k],
                        gen[k] > 0 ? reco[k] / gen[k] : 0.0);
            printf("[xsec-diag/%-3s] %-6s gen %6.2f | reco %6.2f (A*eps = %.3f) | "
                   "r*gen %6.2f +/- %4.2f | r*reco %6.2f +/- %4.2f nb\n",
                   flav[fl], xlab[k], gen[k], reco[k], gen[k] > 0 ? reco[k] / gen[k] : 0.0,
                   rg[k], erg[k], rr[k], err[k]);
        }

        // ---- draw ------------------------------------------------------------
        double xGen[3], xReco[3], xMeas[3], xCnt[3], ex0[3] = {0, 0, 0}, ey0[3] = {0, 0, 0};
        for (int k = 0; k < 3; ++k)
        { xGen[k] = k + 1 - 0.27; xMeas[k] = k + 1 - 0.09; xReco[k] = k + 1 + 0.09; xCnt[k] = k + 1 + 0.27; }
        TGraphErrors *gGen  = new TGraphErrors(3, xGen,  gen,  ex0, ey0);
        TGraphErrors *gMeas = new TGraphErrors(3, xMeas, rg,   ex0, erg);
        TGraphErrors *gReco = new TGraphErrors(3, xReco, reco, ex0, ey0);
        TGraphErrors *gCnt  = new TGraphErrors(3, xCnt,  rr,   ex0, err);

        PlotStyle ps; ps.logy = false;
        TCanvas *c = new TCanvas(Form("c_xsec_diag_%s", flav[fl]), "", ps.w, ps.h);
        ApplyCanvasStyle(c, ps);
        c->cd();
        double ymax = 0.0;
        for (int k = 0; k < 3; ++k)
            ymax = std::max(ymax, std::max(std::max(gen[k], reco[k]),
                                           std::max(rg[k] + erg[k], rr[k] + err[k])));
        TH1F *fr = new TH1F(Form("fr_diag_%s", flav[fl]), "", 3, 0.5, 3.5);
        for (int k = 0; k < 3; ++k) fr->GetXaxis()->SetBinLabel(k + 1, xlab[k]);
        fr->SetStats(0);
        ApplyHistStyle(fr, ps, "", Form("#sigma^{fid}_{W#rightarrow%s#nu} (nb)", lepSym[fl]));
        fr->GetXaxis()->SetLabelSize(0.065);
        fr->GetXaxis()->SetNdivisions(3);
        fr->SetMinimum(0.0);
        fr->SetMaximum(1.55 * ymax);
        fr->Draw();

        gGen->SetMarkerStyle(27); gGen->SetMarkerSize(2.0);  // open diamond  = gen, r = 1
        gGen->SetMarkerColor(kGreen + 2); gGen->SetLineColor(kGreen + 2); gGen->SetLineWidth(2);
        gMeas->SetMarkerStyle(20); gMeas->SetMarkerSize(1.4); // filled circle = r x gen (measured)
        gMeas->SetMarkerColor(kBlack); gMeas->SetLineColor(kBlack); gMeas->SetLineWidth(2);
        gReco->SetMarkerStyle(25); gReco->SetMarkerSize(1.6); // open square   = reco MC, r = 1
        gReco->SetMarkerColor(kAzure + 2); gReco->SetLineColor(kAzure + 2); gReco->SetLineWidth(2);
        gCnt->SetMarkerStyle(21); gCnt->SetMarkerSize(1.4);   // filled square = r x reco (counts)
        gCnt->SetMarkerColor(kAzure + 2); gCnt->SetLineColor(kAzure + 2); gCnt->SetLineWidth(2);
        gGen->Draw("P SAME"); gReco->Draw("P SAME"); gMeas->Draw("P SAME"); gCnt->Draw("P SAME");

        // (A x eps) label next to each gen/reco pair -- the size of the correction
        TLatex te; te.SetTextFont(42); te.SetTextSize(0.026); te.SetTextAlign(21);
        te.SetTextColor(kAzure + 2);
        for (int k = 0; k < 3; ++k)
            te.DrawLatex(k + 1 + 0.18, reco[k] + 0.045 * ymax,
                         Form("A#times#varepsilon = %.2f", gen[k] > 0 ? reco[k] / gen[k] : 0.0));

        DrawHeader(ps, "", Form("W #rightarrow %s #nu", lepSym[fl]), "extraction diagnostic");
        TLatex ltag; ltag.SetNDC(); ltag.SetTextFont(ps.font); ltag.SetTextAlign(13);
        ltag.SetTextSize(ps.boxTextSize);
        ltag.DrawLatex(ps.headerX, ps.headerY - 3.0 * ps.headerDy, Form("%s fit", discLabel.Data()));

        TLegend *leg = new TLegend(0.18, 0.50, 0.60, 0.72);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.030);
        leg->AddEntry(gGen,  "#sigma^{gen}_{fid} (POWHEG, r = 1)", "p");
        leg->AddEntry(gMeas, "r #times #sigma^{gen}_{fid}  (measured)", "lep");
        leg->AddEntry(gReco, "#sigma^{reco}_{MC} = S/L (r = 1)", "p");
        leg->AddEntry(gCnt,  "r #times #sigma^{reco} = N_{fit}/L  (counts)", "lep");
        leg->Draw();

        TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.029);
        tx.DrawLatex(0.18, 0.255, Form("simfit (r shared #mu/e), %s channel only", lepSym[fl]));
        tx.DrawLatex(0.18, 0.215, "fiducial: bare lepton p_{T} > 25 GeV, |#eta_{lab}| < 2.4");
        tx.DrawLatex(0.18, 0.175, "stat. (fit) unc. only");

        CMS_lumi(c, 13, 10);
        c->RedrawAxis(); c->Modified(); c->Update();
        c->SaveAs(Form("%s/W_xsec_diag_%s.png", outDir.c_str(), flav[fl]));
        c->SaveAs(Form("%s/W_xsec_diag_%s.pdf", outDir.c_str(), flav[fl]));
        std::cout << "[OK] Saved " << outDir << "/W_xsec_diag_" << flav[fl] << ".{png,pdf}\n";

        // ---- eta-binned: the same four series as d(sigma)/d(eta_CM) ----------
        if (NY != 12) continue;
        const char *cn[3] = {"Wp", "Wm", "W"};
        for (int k = 0; k < 3; ++k)
        {
            std::vector<double> dg(NY), dr(NY), dm(NY), edm(NY), dc(NY), edc(NY), zero(NY, 0.0);
            double ymaxD = 0;
            for (int i = 0; i < NY; ++i)
            {
                double dEta = 2.0 * exc[i];
                if (dEta <= 0) dEta = 0.4;
                double tg = 0, tr = 0;
                if (k != 1) { tg += in.G[i];      tr += w[i]; }
                if (k != 0) { tg += in.G[12 + i]; tr += w[12 + i]; }
                double v, e;
                binVal(in.G, k, i, v, e); dm[i] = v / dEta; edm[i] = e / dEta;
                binVal(w,    k, i, v, e); dc[i] = v / dEta; edc[i] = e / dEta;
                dg[i] = tg / dEta;  dr[i] = tr / dEta;
                axeps[fl][k][i] = tg > 0 ? tr / tg : 0.0;
                ymaxD = std::max(ymaxD, std::max(dg[i], dm[i] + edm[i]));
                csvB << Form("%s,%s,%d,%.4f,%.4f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n",
                             flav[fl], cn[k], i, xc[i], exc[i], dg[i], dr[i],
                             dm[i], edm[i], dc[i], edc[i], axeps[fl][k][i]);
            }

            TGraph *ggen = new TGraph(NY, xc.data(), dg.data());
            TGraph *grec = new TGraph(NY, xc.data(), dr.data());
            TGraphErrors *gmea = new TGraphErrors(NY, xc.data(), dm.data(), exc.data(), edm.data());
            TGraphErrors *gcnt = new TGraphErrors(NY, xc.data(), dc.data(), exc.data(), edc.data());

            PlotStyle ps2; ps2.logy = false;
            TCanvas *cd = new TCanvas(Form("c_diag_deta_%s_%s", flav[fl], cn[k]), "", ps2.w, ps2.h);
            ApplyCanvasStyle(cd, ps2);
            cd->cd();
            const double xlo = xc[0] - exc[0] - 0.1, xhi = xc[NY - 1] + exc[NY - 1] + 0.1;
            TH1F *frD = new TH1F(Form("fr_diag_deta_%s_%s", flav[fl], cn[k]), "", 100, xlo, xhi);
            frD->SetStats(0);
            ApplyHistStyle(frD, ps2, "#eta^{l}_{CM}",
                           Form("d#sigma^{fid}_{W#rightarrow%s#nu}/d#eta (nb)", lepSym[fl]));
            frD->SetMinimum(0.0);
            frD->SetMaximum(1.60 * ymaxD);
            frD->Draw();

            ggen->SetLineStyle(2); ggen->SetLineWidth(3); ggen->SetLineColor(kGreen + 2);
            grec->SetLineStyle(3); grec->SetLineWidth(3); grec->SetLineColor(kAzure + 2);
            ggen->Draw("L SAME"); grec->Draw("L SAME");
            gmea->SetMarkerStyle(20); gmea->SetMarkerSize(1.3); gmea->SetMarkerColor(kBlack);
            gmea->SetLineColor(kBlack); gmea->SetLineWidth(2);
            gcnt->SetMarkerStyle(21); gcnt->SetMarkerSize(1.3); gcnt->SetMarkerColor(kAzure + 2);
            gcnt->SetLineColor(kAzure + 2); gcnt->SetLineWidth(2);
            gmea->Draw("P SAME"); gcnt->Draw("P SAME");

            const char *procD = (k == 0) ? "W^{+} #rightarrow %s #nu"
                              : (k == 1) ? "W^{-} #rightarrow %s #nu"
                                         : "W #rightarrow %s #nu";
            DrawHeader(ps2, "", Form(procD, lepSym[fl]), "extraction diagnostic");
            TLatex ltD; ltD.SetNDC(); ltD.SetTextFont(ps2.font); ltD.SetTextAlign(13);
            ltD.SetTextSize(ps2.boxTextSize);
            ltD.DrawLatex(ps2.headerX, ps2.headerY - 3.0 * ps2.headerDy, Form("%s fit", discLabel.Data()));

            TLegend *legD = new TLegend(0.18, 0.55, 0.60, 0.78);
            legD->SetBorderSize(0); legD->SetFillStyle(0); legD->SetTextFont(42); legD->SetTextSize(0.030);
            legD->AddEntry(ggen, "#sigma^{gen}_{fid} (POWHEG, r = 1)", "l");
            legD->AddEntry(gmea, "r #times #sigma^{gen}_{fid}  (measured)", "lep");
            legD->AddEntry(grec, "#sigma^{reco}_{MC} = S/L (r = 1)", "l");
            legD->AddEntry(gcnt, "r #times #sigma^{reco} = N_{fit}/L  (counts)", "lep");
            legD->Draw();

            TLatex txD; txD.SetNDC(); txD.SetTextFont(42); txD.SetTextSize(0.029);
            txD.DrawLatex(0.18, 0.255, Form("simfit (r shared #mu/e), %s channel only", lepSym[fl]));
            txD.DrawLatex(0.18, 0.215, "gen #rightarrow reco gap = (A#times#varepsilon)_{MC}");
            txD.DrawLatex(0.18, 0.175, "stat. (fit) unc. only");

            CMS_lumi(cd, 13, 10);
            cd->RedrawAxis(); cd->Modified(); cd->Update();
            cd->SaveAs(Form("%s/W_dsigma_diag_%s_%s.png", outDir.c_str(), flav[fl], cn[k]));
            cd->SaveAs(Form("%s/W_dsigma_diag_%s_%s.pdf", outDir.c_str(), flav[fl], cn[k]));
        }
        std::cout << "[OK] Saved " << outDir << "/W_dsigma_diag_" << flav[fl]
                  << "_{Wp,Wm,W}.{png,pdf}\n";
    }

    // ---- (A x eps)(eta) summary: mu vs e, W+ and W- --------------------------
    // The MC correction the r x sigma_gen extraction applies, bin by bin. An
    // electron-only dip is the ECAL-crack acceptance hole (|eta_SC| 1.44-1.57).
    if (NY == 12 && axeps[0][0][0] > 0 && axeps[1][0][0] > 0)
    {
        PlotStyle ps3; ps3.logy = false;
        TCanvas *ca = new TCanvas("c_axeps", "", ps3.w, ps3.h);
        ApplyCanvasStyle(ca, ps3);
        ca->cd();
        double vlo = 1e9, vhi = 0;
        for (int fl = 0; fl < 2; ++fl)
            for (int k = 0; k < 2; ++k)
                for (int i = 0; i < NY; ++i)
                { vlo = std::min(vlo, axeps[fl][k][i]); vhi = std::max(vhi, axeps[fl][k][i]); }
        const double xlo = xc[0] - exc[0] - 0.1, xhi = xc[NY - 1] + exc[NY - 1] + 0.1;
        TH1F *frA = new TH1F("fr_axeps", "", 100, xlo, xhi);
        frA->SetStats(0);
        ApplyHistStyle(frA, ps3, "#eta^{l}_{CM}", "(A #times #varepsilon)_{MC} = #sigma^{reco}/#sigma^{gen}_{fid}");
        frA->SetMinimum(std::max(0.0, vlo - 0.45 * (vhi - vlo) - 0.05));
        frA->SetMaximum(vhi + 0.55 * (vhi - vlo) + 0.05);
        frA->Draw();

        TGraph *ga[2][2];
        const int col[2] = {kAzure + 2, kRed + 1};   // charge
        const int mst[2][2] = {{20, 22}, {24, 26}};  // [flavour][charge]: mu filled, e open
        for (int fl = 0; fl < 2; ++fl)
            for (int k = 0; k < 2; ++k)
            {
                ga[fl][k] = new TGraph(NY, xc.data(), axeps[fl][k]);
                ga[fl][k]->SetMarkerStyle(mst[fl][k]); ga[fl][k]->SetMarkerSize(1.4);
                ga[fl][k]->SetMarkerColor(col[k]);     ga[fl][k]->SetLineColor(col[k]);
                ga[fl][k]->SetLineWidth(2);            ga[fl][k]->SetLineStyle(fl == 0 ? 1 : 2);
                ga[fl][k]->Draw("PL SAME");
            }

        DrawHeader(ps3, "", "W #rightarrow l #nu", "acceptance #times efficiency");
        TLatex ltA; ltA.SetNDC(); ltA.SetTextFont(ps3.font); ltA.SetTextAlign(13);
        ltA.SetTextSize(ps3.boxTextSize);
        ltA.DrawLatex(ps3.headerX, ps3.headerY - 3.0 * ps3.headerDy, Form("%s fit", discLabel.Data()));

        TLegend *legA = new TLegend(0.18, 0.60, 0.52, 0.80);
        legA->SetBorderSize(0); legA->SetFillStyle(0); legA->SetTextFont(42); legA->SetTextSize(0.030);
        legA->SetNColumns(2);
        legA->AddEntry(ga[0][0], "#mu, W^{+}", "lp");
        legA->AddEntry(ga[1][0], "e, W^{+}", "lp");
        legA->AddEntry(ga[0][1], "#mu, W^{-}", "lp");
        legA->AddEntry(ga[1][1], "e, W^{-}", "lp");
        legA->Draw();

        TLatex txA; txA.SetNDC(); txA.SetTextFont(42); txA.SetTextSize(0.029);
        txA.DrawLatex(0.18, 0.255, "MC correction implicit in #sigma = r #times #sigma^{gen}_{fid}");
        txA.DrawLatex(0.18, 0.215, "fiducial: bare lepton p_{T} > 25 GeV, |#eta_{lab}| < 2.4");

        CMS_lumi(ca, 13, 10);
        ca->RedrawAxis(); ca->Modified(); ca->Update();
        ca->SaveAs((outDir + "/AxEps_vs_eta.png").c_str());
        ca->SaveAs((outDir + "/AxEps_vs_eta.pdf").c_str());
        std::cout << "[OK] Saved " << outDir << "/AxEps_vs_eta.{png,pdf}\n";
    }

    printf("[xsec-diag] cross-check: mu+e of r*sigma^reco = %.2f / %.2f / %.2f nb"
           " (should equal the fork's {Wp,Wm,W}_sum_yield / L)\n",
           sumCount[0], sumCount[1], sumCount[2]);
    std::cout << "[OK] Wrote " << outDir << "/xsec_diag.csv + xsec_diag_bins.csv\n";
    if (fC) { fC->Close(); delete fC; }
}
