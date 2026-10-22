// correction/open_issues_diag.C
//
// Diagnostics of the two open issues of the analysis note (AN-26-030,
// Section "Outstanding issues", label sec:openissues):
//
//   (A) the charge-dependent electron/muon difference of the W yields
//       (the electron W- yield low in the outermost rapidity bins);
//   (B) the Z->ee vs Z->mumu normalization behind the per-flavour r_Z.
//
// Reads ONLY stored outputs -- no ntuple loop, no fit:
//   ../plotting/plots[/Elec]/combine_input_W_leppt_mt40.root  (the fitted inputs)
//   ../plotting/plots[/Elec]/combine_input_W.root             (PF-MET inputs)
//   ../plotting/plots[/Elec]/combine_input_Z.root
//   ../skim/rootfile/ZTo{MuMu,EE}_pO2025_{Data,DY_MC}_hist.root
//
// Every number is PREFIT (absolute MC, r = 1) and every error is the data
// Poisson error only (sqrt(N)); the MC statistical error is neglected.
//
// Definitions per region (flavour f, charge c, lab rapidity bin i):
//   D  = data,  W = signal + wtau (what the fit's r scales),
//   B  = z + ztau + ttbar (every other MC),  Q = the multijet template.
//   r_Q = (D - B - Q) / W     (= prefit_r.C's r_pre[_abcd])
//   r_0 = (D - B) / W         (no multijet at all: an UPPER bound on r)
//   f_Q = (D - B - r_mu W) / Q   (electron rows only): the factor on the
//         electron multijet template that would give that electron bin the
//         signal strength r_Q of the muon bin of the same charge and rapidity.
//
// Part B normalizes the DY skim histograms to the absolute scale of the
// Combine input (k = Z_incl/signal / hMass integral, the same k_s), so the
// lepton-leg and mass distributions are compared with data in absolute terms.
// Non-DY MC (< 0.5% under the peaks) is not included in Part B's leg
// distributions.
//
// Run via ./run_open_issues_diag.sh -> logs/open_issues_diag.log (the record).

#include "../plotting/plotting_helper.C"

#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TSystem.h"
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

const int kNY = 12;
const char *kCharge[2] = {"Wp", "Wm"};
const char *kPlotDir = "plots/open_issues";

// per-bin electron/muon ratios of one Part-A variant, kept for the figures
struct EMuCurve {
  double rq[2][12] = {{0}}, erq[2][12] = {{0}}; // e/mu of r_Q
  double r0[2][12] = {{0}}, er0[2][12] = {{0}}; // e/mu of r_0
};

double Integ(TFile *f, const TString &dir, const char *proc, int b1, int b2)
{
  TH1 *h = (TH1 *)f->Get(dir + "/" + proc);
  if (!h) return 0.;
  if (b2 < 0) b2 = h->GetNbinsX();
  return h->Integral(b1, b2);
}

struct Reg { double D = 0, W = 0, B = 0, Q = 0; };

Reg ReadReg(TFile *f, const TString &dir, const char *qcdName, int b1, int b2)
{
  Reg r;
  r.D = Integ(f, dir, "data_obs", b1, b2);
  r.W = Integ(f, dir, "signal", b1, b2) + Integ(f, dir, "wtau", b1, b2);
  r.B = Integ(f, dir, "z", b1, b2) + Integ(f, dir, "ztau", b1, b2) + Integ(f, dir, "ttbar", b1, b2);
  r.Q = Integ(f, dir, qcdName, b1, b2);
  return r;
}

double RQ(const Reg &r) { return (r.D - r.B - r.Q) / r.W; }
double R0(const Reg &r) { return (r.D - r.B) / r.W; }
double RErr(const Reg &r) { return std::sqrt(r.D) / r.W; }

void Ratio(double a, double ea, double b, double eb, double &q, double &eq)
{
  q = a / b;
  eq = std::fabs(q) * std::sqrt(std::pow(ea / a, 2) + std::pow(eb / b, 2));
}

// Part A: one discriminant, bins [b1, b2] of every region.
void PartA(const char *title, const char *fileStem, const char *qcdName, int b1, int b2, EMuCurve &out)
{
  TFile *fm = TFile::Open(Form("../plotting/plots/%s.root", fileStem));
  TFile *fe = TFile::Open(Form("../plotting/plots/Elec/%s.root", fileStem));
  if (!fm || fm->IsZombie() || !fe || fe->IsZombie()) {
    printf("[ERROR] cannot open %s inputs\n", fileStem);
    return;
  }
  TH1 *href = (TH1 *)fm->Get("Wp_lab_y0/data_obs");
  printf("\n=== (A) %s : %s.root, multijet template '%s', bins %d..%d = [%g, %g) ===\n",
         title, fileStem, qcdName, b1, b2 < 0 ? href->GetNbinsX() : b2,
         href->GetXaxis()->GetBinLowEdge(b1),
         href->GetXaxis()->GetBinUpEdge(b2 < 0 ? href->GetNbinsX() : b2));
  printf("    r_Q = (D-B-Q)/W, r_0 = (D-B)/W (no multijet), errors = data Poisson only\n");

  for (int c = 0; c < 2; ++c) {
    printf("\n  %s   bin |   mu: D     Q/D    r_Q            r_0    |   e: D     Q/D    r_Q            r_0    "
           "| e/mu r_Q        e/mu r_0       | f_Q (e)\n", kCharge[c]);
    Reg sm[kNY], se[kNY];
    Reg edgeM, edgeE;
    for (int i = 0; i < kNY; ++i) {
      TString dir = Form("%s_lab_y%d", kCharge[c], i);
      Reg m = ReadReg(fm, dir, qcdName, b1, b2);
      Reg e = ReadReg(fe, dir, qcdName, b1, b2);
      sm[i] = m; se[i] = e;
      if (i == 0 || i == kNY - 1) {
        edgeM.D += m.D; edgeM.W += m.W; edgeM.B += m.B; edgeM.Q += m.Q;
        edgeE.D += e.D; edgeE.W += e.W; edgeE.B += e.B; edgeE.Q += e.Q;
      }
      double qq, eqq, q0, eq0;
      Ratio(RQ(e), RErr(e), RQ(m), RErr(m), qq, eqq);
      Ratio(R0(e), RErr(e), R0(m), RErr(m), q0, eq0);
      out.rq[c][i] = qq; out.erq[c][i] = eqq; out.r0[c][i] = q0; out.er0[c][i] = eq0;
      double fQ = (e.D - e.B - RQ(m) * e.W) / e.Q;
      double efQ = std::sqrt(e.D + std::pow(RErr(m) * e.W, 2)) / e.Q;
      printf("        y%-2d |  %5.0f  %5.3f  %5.3f+-%5.3f  %5.3f  |  %5.0f  %5.3f  %5.3f+-%5.3f  %5.3f  "
             "| %5.3f+-%5.3f  %5.3f+-%5.3f | %5.2f+-%4.2f\n",
             i, m.D, m.Q / m.D, RQ(m), RErr(m), R0(m), e.D, e.Q / e.D, RQ(e), RErr(e), R0(e),
             qq, eqq, q0, eq0, fQ, efQ);
    }
    double qq, eqq, q0, eq0;
    Ratio(RQ(edgeE), RErr(edgeE), RQ(edgeM), RErr(edgeM), qq, eqq);
    Ratio(R0(edgeE), RErr(edgeE), R0(edgeM), RErr(edgeM), q0, eq0);
    double fQ = (edgeE.D - edgeE.B - RQ(edgeM) * edgeE.W) / edgeE.Q;
    double efQ = std::sqrt(edgeE.D + std::pow(RErr(edgeM) * edgeE.W, 2)) / edgeE.Q;
    printf("    y0+y11 |  %5.0f  %5.3f  %5.3f+-%5.3f  %5.3f  |  %5.0f  %5.3f  %5.3f+-%5.3f  %5.3f  "
           "| %5.3f+-%5.3f  %5.3f+-%5.3f | %5.2f+-%4.2f\n",
           edgeM.D, edgeM.Q / edgeM.D, RQ(edgeM), RErr(edgeM), R0(edgeM), edgeE.D, edgeE.Q / edgeE.D,
           RQ(edgeE), RErr(edgeE), R0(edgeE), qq, eqq, q0, eq0, fQ, efQ);
    // inclusive over the 12 lab bins
    Reg im, ie;
    for (int i = 0; i < kNY; ++i) {
      im.D += sm[i].D; im.W += sm[i].W; im.B += sm[i].B; im.Q += sm[i].Q;
      ie.D += se[i].D; ie.W += se[i].W; ie.B += se[i].B; ie.Q += se[i].Q;
    }
    Ratio(RQ(ie), RErr(ie), RQ(im), RErr(im), qq, eqq);
    Ratio(R0(ie), RErr(ie), R0(im), RErr(im), q0, eq0);
    printf("      incl |  %5.0f  %5.3f  %5.3f+-%5.3f  %5.3f  |  %5.0f  %5.3f  %5.3f+-%5.3f  %5.3f  "
           "| %5.3f+-%5.3f  %5.3f+-%5.3f |\n",
           im.D, im.Q / im.D, RQ(im), RErr(im), R0(im), ie.D, ie.Q / ie.D, RQ(ie), RErr(ie), R0(ie),
           qq, eqq, q0, eq0);
  }
  fm->Close(); fe->Close();
}

// Part A3: the inclusive fitted cross sections of the three fits
// (../plotting/plots/flavfit/xsec/leppt_mt40/xsec_flavfit.csv, rows <fit>,W[pm],incl):
// the electron/muon double ratio (W+ e/mu)/(W- e/mu) and the inclusive charge
// asymmetry per fit, with STATISTICAL errors. The W+/W- statistical covariance
// within one fit is neglected (0.3% correlation in the combined fit's stat matrix,
// xsec_contour_WZ_lab.csv), and the e and mu fits are independent samples.
void PartA3()
{
  const char *csv = "../plotting/plots/flavfit/xsec/leppt_mt40/xsec_flavfit.csv";
  FILE *fp = std::fopen(csv, "r");
  if (!fp) { printf("[ERROR] cannot open %s\n", csv); return; }
  const char *fits[3] = {"simfit_mu", "simfit_ele", "comb"};
  double s[3][2] = {{0}}, es[3][2] = {{0}};
  char line[1024];
  while (std::fgets(line, sizeof(line), fp)) {
    TString l(line);
    // Tokenize collapses the empty etaCM fields of the incl rows, so they read
    // fit,charge,incl,r,rErr,sigma_gen,sigma_meas,err,stat,syst (10 tokens)
    TObjArray *tok = l.Tokenize(",");
    if (tok->GetEntries() == 10 && TString(((TObjString *)tok->At(2))->GetString()) == "incl") {
      TString fit = ((TObjString *)tok->At(0))->GetString();
      TString chg = ((TObjString *)tok->At(1))->GetString();
      for (int f = 0; f < 3; ++f)
        for (int c = 0; c < 2; ++c)
          if (fit == fits[f] && chg == kCharge[c]) {
            s[f][c] = TString(((TObjString *)tok->At(6))->GetString()).Atof();  // sigma_meas_nb
            es[f][c] = TString(((TObjString *)tok->At(8))->GetString()).Atof(); // sigma_meas_stat_nb
          }
    }
    delete tok;
  }
  std::fclose(fp);
  printf("\n=== (A3) inclusive fitted cross sections (%s), stat errors ===\n", csv);
  for (int f = 0; f < 3; ++f) {
    double p = s[f][0], m = s[f][1], S = p + m;
    double A = (p - m) / S;
    double eA = 2. / (S * S) * std::sqrt(m * m * es[f][0] * es[f][0] + p * p * es[f][1] * es[f][1]);
    printf("    %-10s sigma(W+) %6.2f +- %4.2f  sigma(W-) %6.2f +- %4.2f  A_ch(incl) = %.4f +- %.4f\n",
           fits[f], p, es[f][0], m, es[f][1], A, eA);
  }
  double rp, erp, rm, erm, dr, edr;
  Ratio(s[1][0], es[1][0], s[0][0], es[0][0], rp, erp);
  Ratio(s[1][1], es[1][1], s[0][1], es[0][1], rm, erm);
  Ratio(rp, erp, rm, erm, dr, edr);
  printf("    e/mu: W+ %.4f +- %.4f   W- %.4f +- %.4f   double ratio (W+ e/mu)/(W- e/mu) = %.3f +- %.3f  (%.1f sigma from 1)\n",
         rp, erp, rm, erm, dr, edr, (dr - 1) / edr);
  auto acc = [&](int f, double &A, double &eA) {
    double p = s[f][0], m = s[f][1], S = p + m;
    A = (p - m) / S;
    eA = 2. / (S * S) * std::sqrt(m * m * es[f][0] * es[f][0] + p * p * es[f][1] * es[f][1]);
  };
  double Am, eAm, Ae, eAe;
  acc(0, Am, eAm); acc(1, Ae, eAe);
  printf("    A_ch(incl) e - mu = %+.4f +- %.4f  (%.1f sigma)\n", Ae - Am, std::hypot(eAm, eAe),
         (Ae - Am) / std::hypot(eAm, eAe));
}

// Part B helpers
struct ZFiles { TFile *data = nullptr, *dy = nullptr; double k = 0; };

ZFiles OpenZ(bool isElec)
{
  ZFiles z;
  const char *stem = isElec ? "ZToEE_pO2025" : "ZToMuMu_pO2025";
  z.data = TFile::Open(Form("../skim/rootfile/%s_Data_hist.root", stem));
  z.dy = TFile::Open(Form("../skim/rootfile/%s_DY_MC_hist.root", stem));
  TFile *ci = TFile::Open(Form("../plotting/plots/%scombine_input_Z.root", isElec ? "Elec/" : ""));
  TH1 *sig = ci ? (TH1 *)ci->Get("Z_incl/signal") : nullptr;
  TH1 *hm = z.dy ? (TH1 *)z.dy->Get("hMass") : nullptr;
  if (sig && hm) z.k = sig->Integral() / hm->Integral();
  return z;
}

// prints data/MC per range and returns it as a graph (x = range centre, ex = half width)
TGraphErrors *PrintRatioBins(const char *flav, TH1 *hd, TH1 *hm, double k, const std::vector<double> &edges, bool absval)
{
  TGraphErrors *g = new TGraphErrors();
  for (size_t j = 0; j + 1 < edges.size(); ++j) {
    double d = 0, m = 0;
    for (int b = 1; b <= hd->GetNbinsX(); ++b) {
      double x = hd->GetXaxis()->GetBinCenter(b);
      if (absval) x = std::fabs(x);
      if (x >= edges[j] && x < edges[j + 1]) { d += hd->GetBinContent(b); m += hm->GetBinContent(b) * k; }
    }
    if (m <= 0) continue;
    printf("    %-3s [%6.2f, %6.2f)  data %5.0f  MC %7.1f  data/MC %5.3f +- %5.3f\n",
           flav, edges[j], edges[j + 1], d, m, d / m, std::sqrt(d) / m);
    const int n = g->GetN();
    g->SetPoint(n, 0.5 * (edges[j] + edges[j + 1]), d / m);
    g->SetPointError(n, 0.5 * (edges[j + 1] - edges[j]), std::sqrt(d) / m);
  }
  return g;
}

// ---- figures ---------------------------------------------------------------
void DrawRatioOverlay(std::vector<OverlaySeries> ser, const char *name, const char *xTitle, const char *yTitle,
                      const char *head1, const char *head2, const std::vector<std::string> &box,
                      double ymin, double ymax)
{
  PlotStyle ps;
  // header (two lines) and box right of the CMS label; the main-title slot is
  // empty, so the two sub-lines land at headerY - dy and headerY - 2 dy
  ps.headerX = 0.47;
  ps.headerY = 0.94;
  ps.subSize = 0.032;
  ps.boxTextSize = 0.026;
  ps.boxX1 = 0.47; ps.boxX2 = 0.93; ps.boxY1 = 0.71; ps.boxY2 = 0.77;
  ps.legX1 = 0.55; ps.legY1 = 0.15; ps.legX2 = 0.93; ps.legY2 = 0.27;
  GraphTuner tuner = [ymin, ymax](TCanvas *c, TGraphErrors *g) {
    g->SetMinimum(ymin);
    g->SetMaximum(ymax);
    c->cd();
    TLine *one = new TLine(g->GetXaxis()->GetXmin(), 1.0, g->GetXaxis()->GetXmax(), 1.0);
    one->SetLineStyle(2);
    one->SetLineColor(kGray + 2);
    one->Draw();
  };
  SaveNiceGraph_Overlay(ser, Form("%s/%s", kPlotDir, name), xTitle, yTitle, "", head1, head2, box, ps, tuner);
}

// e/mu vs eta_lab for W+ and W- (statistical errors of the data only)
void PlotEMu(const EMuCurve &cv, bool useR0, const char *name, const char *head1, const char *head2)
{
  std::vector<OverlaySeries> ser(2);
  for (int c = 0; c < 2; ++c) {
    TGraphErrors *g = new TGraphErrors(kNY);
    for (int i = 0; i < kNY; ++i) {
      g->SetPoint(i, -2.2 + 0.4 * i, useR0 ? cv.r0[c][i] : cv.rq[c][i]);
      g->SetPointError(i, 0.2, useR0 ? cv.er0[c][i] : cv.erq[c][i]);
    }
    ser[c].gTot = g;
    ser[c].label = c == 0 ? "W^{+}" : "W^{-}";
    ser[c].color = c == 0 ? kAzure + 2 : kRed + 1; // the repo's W+/W- convention (xsec_fiducial.C)
    ser[c].marker = c == 0 ? 21 : 33;
    ser[c].markerSize = c == 0 ? 1.5 : 2.0;
    ser[c].xShift = c == 0 ? -0.35 : 0.35;
  }
  DrawRatioOverlay(ser, name, "#eta_{lab} (proton-going > 0)", useR0 ? "e / #mu  of  (D - B) / W" : "e / #mu  of  (D - B - Q) / W",
                   head1, head2, {"prefit, absolute MC, data statistical errors only"}, 0.2, 2.2);
}

void PlotZ(TGraphErrors *gmu, TGraphErrors *ge, const char *name, const char *xTitle)
{
  std::vector<OverlaySeries> ser(2);
  ser[0].gTot = gmu; ser[0].label = "Z#rightarrow#mu#mu"; ser[0].color = kBlue + 1; ser[0].marker = 20;
  ser[0].markerSize = 1.3; ser[0].xShift = -0.25;
  ser[1].gTot = ge; ser[1].label = "Z#rightarrowee"; ser[1].color = kRed + 1; ser[1].marker = 21;
  ser[1].markerSize = 1.5; ser[1].xShift = 0.25;
  DrawRatioOverlay(ser, name, xTitle, "data / MC (absolute)", "Z#rightarrow#mu#mu and Z#rightarrowee, prefit",
                   "60 < m_{ll} < 120 GeV, both legs", {"data statistical errors only"}, 0.4, 2.2);
}

void PartB()
{
  printf("\n=== (B) Z->ee vs Z->mumu (prefit; MC = DY skim histograms x k, k from combine_input_Z.root) ===\n");
  ZFiles zz[2] = {OpenZ(false), OpenZ(true)};
  const char *fl[2] = {"mu", "e"};
  double ratio[2] = {0, 0}, eratio[2] = {0, 0};

  printf("\n  (B1) Z peak, 60 < m < 120 GeV, all MC of the Combine input\n");
  for (int f = 0; f < 2; ++f) {
    TFile *ci = TFile::Open(Form("../plotting/plots/%scombine_input_Z.root", f ? "Elec/" : ""));
    double D = Integ(ci, "Z_incl", "data_obs", 1, -1);
    double S = Integ(ci, "Z_incl", "signal", 1, -1);
    double O = Integ(ci, "Z_incl", "ztau", 1, -1) + Integ(ci, "Z_incl", "w", 1, -1) +
               Integ(ci, "Z_incl", "wtau", 1, -1) + Integ(ci, "Z_incl", "ttbar", 1, -1);
    ratio[f] = D / (S + O); eratio[f] = std::sqrt(D) / (S + O);
    printf("    %-3s data %4.0f  DY %7.2f  other MC %5.2f  data/MC %5.3f +- %5.3f   (k = %.4g)\n",
           fl[f], D, S, O, ratio[f], eratio[f], zz[f].k);
  }
  double q, eq;
  Ratio(ratio[1], eratio[1], ratio[0], eratio[0], q, eq);
  printf("    e/mu of data/MC = %5.3f +- %5.3f\n", q, eq);

  printf("\n  (B2) lepton legs (both legs, h_lepPt), data/MC per leg-pT range\n");
  std::vector<double> ptE = {10, 20, 26, 35, 45, 60, 100.0001};
  TGraphErrors *gPt[2];
  for (int f = 0; f < 2; ++f)
    gPt[f] = PrintRatioBins(fl[f], (TH1 *)zz[f].data->Get("h_lepPt"), (TH1 *)zz[f].dy->Get("h_lepPt"), zz[f].k, ptE, false);
  PlotZ(gPt[0], gPt[1], "z_datamc_vs_legpt", "lepton p_{T} (GeV)");

  printf("\n  (B3) lepton legs, data/MC per |eta| range (h_lepEta; the ECAL crack 1.4442-1.566 is vetoed for e)\n");
  std::vector<double> etaE = {0, 0.8, 1.4, 1.6, 2.0, 2.41};
  for (int f = 0; f < 2; ++f)
    PrintRatioBins(fl[f], (TH1 *)zz[f].data->Get("h_lepEta"), (TH1 *)zz[f].dy->Get("h_lepEta"), zz[f].k, etaE, true);

  printf("\n  (B4) dilepton rapidity, data/MC per |y_Z| range (h_Zy)\n");
  std::vector<double> yE = {0, 0.8, 1.6, 2.41};
  for (int f = 0; f < 2; ++f)
    PrintRatioBins(fl[f], (TH1 *)zz[f].data->Get("h_Zy"), (TH1 *)zz[f].dy->Get("h_Zy"), zz[f].k, yE, true);

  printf("\n  (B5) dilepton mass: data/MC per mass range, and the mean in [80, 100) GeV (hMass, 0.5 GeV bins)\n");
  std::vector<double> mE = {60, 80, 89, 92, 100, 120.0001};
  TGraphErrors *gM[2];
  for (int f = 0; f < 2; ++f)
    gM[f] = PrintRatioBins(fl[f], (TH1 *)zz[f].data->Get("hMass"), (TH1 *)zz[f].dy->Get("hMass"), zz[f].k, mE, false);
  PlotZ(gM[0], gM[1], "z_datamc_vs_mass", "m_{ll} (GeV)");
  for (int f = 0; f < 2; ++f) {
    TH1 *hd = (TH1 *)zz[f].data->Get("hMass");
    TH1 *hm = (TH1 *)zz[f].dy->Get("hMass");
    TH1 *cd = (TH1 *)hd->Clone(Form("cd_%d", f)); cd->SetDirectory(nullptr);
    TH1 *cm = (TH1 *)hm->Clone(Form("cm_%d", f)); cm->SetDirectory(nullptr);
    cd->GetXaxis()->SetRangeUser(80.0001, 99.9999);
    cm->GetXaxis()->SetRangeUser(80.0001, 99.9999);
    double nd = cd->Integral();
    printf("    %-3s mean[80,100): data %.3f +- %.3f GeV (N = %.0f, RMS %.2f)  MC %.3f GeV (RMS %.2f)  "
           "data - MC = %+.3f +- %.3f GeV\n",
           fl[f], cd->GetMean(), cd->GetRMS() / std::sqrt(nd), nd, cd->GetRMS(), cm->GetMean(), cm->GetRMS(),
           cd->GetMean() - cm->GetMean(), cd->GetRMS() / std::sqrt(nd));
  }
}

} // namespace

void open_issues_diag()
{
  printf("[open_issues_diag] prefit diagnostics of the W charge split and the Z->ee normalization\n");
  // (A1) the fitted discriminant: lepton pT at m_T > 40, the in-fit-ABCD multijet template
  gSystem->mkdir(kPlotDir, kTRUE);
  EMuCurve cvPt, cvMet;
  PartA("lepton pT, m_T > 40 GeV (fit window)", "combine_input_W_leppt_mt40", "qcd_abcd", 1, -1, cvPt);
  // (A2) PF MET >= 30 GeV: the multijet is suppressed independently of m_T
  //      (bins of 2 GeV: bin 16 = [30, 32)); Q = the ABCD MET template of the MET fit
  PartA("PF MET >= 30 GeV (PF-MET inputs)", "combine_input_W", "qcd", 16, -1, cvMet);
  PlotEMu(cvPt, false, "emu_ptmt40_abcd", "W#rightarrowl#nu, lepton p_{T}, m_{T} > 40 GeV",
          "Q = in-fit ABCD multijet template");
  PlotEMu(cvMet, true, "emu_met30_noqcd", "W#rightarrowl#nu, p_{T}^{miss} #geq 30 GeV",
          "no multijet subtracted (Q = 0)");
  PartA3();
  PartB();
  printf("[OUT] %s/{emu_ptmt40_abcd,emu_met30_noqcd,z_datamc_vs_legpt,z_datamc_vs_mass}.{png,pdf}\n", kPlotDir);
  printf("\n[open_issues_diag] done\n");
}
