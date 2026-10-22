// correction/qcd_ysplit.C
//
// Test of the rapidity split of the W signal-region multijet (QCD) template --
// step 1 of the charge-split follow-up of the analysis note (AN-26-030,
// Section "Outstanding issues", subsection sec:openissues:charge).
//
// The SR multijet template of the leppt_mt40 fit (qcd_abcd) is built
// y-INCLUSIVELY per (flavour, charge) by qcd_abcd.C (total A0 = B0*C40/D0, the
// electron one fake-factor weighted, total A0ff) and split over the 12 lab /
// 12 FB rapidity bins by plotting/mtandmet.C (qcdWeights): the per-bin
// iso-pass low-MET excess, data - EWK(r = 1) at MET < 30, normalized to 1.
// That ASSUMES the SR multijet sits in rapidity like the low-MET excess.
//
// This macro measures the ABCD regions PER RAPIDITY BIN from the y-binned
// twins of the (iso, m_T) plane that the W skim stores since 2026-10-19,
//   h_iso_mt_{mu,ele}{Plus,Minus}_y<i>[_FB]
// (B = iso-pass m_T < 30, C40 = anti-iso m_T > 40, D = anti-iso m_T < 30, the
// regions and windows of qcd_abcd.C) and compares three estimates of the SR
// multijet yield per rapidity bin and charge:
//   tmpl_y   the template the fit uses: Combine input <C>_<b>_y<i>/qcd_abcd
//   ABCD_y = B_y * C40_y / D_y   (the transfer factor flat WITHIN each y bin only)
//   C40s_y = A0 * C40_y / Sum C40 (the anti-iso SR shape: T flat in y)
// and the per-y transfer factors T_y = B_y/D_y (iso, at m_T < 30) and
// R_y = C40_y/D_y (m_T, in the anti-iso sideband). If the low-MET excess
// tracks B_y, then ABCD_y/tmpl_y is essentially R_y/R: the template assumes
// R flat in rapidity.
//
// Also printed per bin: Q_req = data - r_W*(signal + wtau) - (z + ztau + ttbar)
// in the fitted region itself (Combine input), i.e. the multijet the data
// would need for a common W scale r_W (default 1.2 = the fitted level).
// Then: the anti-iso window variation of the edge result, and the e/mu ratio
// of the prefit r_Q (open_issues_diag.C's convention) with the template and
// with the template total split in y like the per-y ABCD.
//
// FB binning: the FB window [-1.7068, 2.4] misses y < -1.7068, so the fb
// templates must sum to the window's share of the inclusive multijet. Until
// 2026-10-19 mtandmet.C normalized the FB weights to 1 over the FB bins (the
// full inclusive total in the window: x1.10 mu / x1.24-1.34 e too much),
// found by this macro; the "[FB window]" line compares the fb template total
// with the window's share (low-MET excess and per-y ABCD) and flags a mismatch.
//
// Conventions:
//  * QCD-only counts = data - r*W-EWK - other EWK, absolute MC (mc_norm.h k_s);
//    W-EWK = Wp, Wm, Wptau, Wmtau; other = DY, DYtau, TT (qcd_abcd.C's split).
//    Two subtractions: r = 1 (the prefit convention of the template) and r = r_W.
//  * errors: data Poisson + EWK MC stat. The ratio ABCD_y/tmpl_y carries the
//    error of ABCD_y only (the template is treated as exact; its own
//    statistical error comes from the low-MET iso-pass events, which largely
//    overlap B_y, so B_y's error stands in for it once).
//  * "abs" ratio = ABCD_y * (Sum tmpl / A0) / tmpl_y: the per-y estimate scaled
//    by the factor the template total carries beyond the flat A0 (1 for the
//    muon, A0ff/A0 = the fake-factor weighting for the electron);
//    "shape" ratio = (ABCD_y / Sum ABCD) / (tmpl_y / Sum tmpl).
//  * chi2 vs 1 over the 12 bins uses Gaussian errors -- only indicative where
//    the regions hold few events (the muon edge bins).
//  * x axis: the lab bins at y = -eta_lab (proton-going > 0, the skim's
//    convention; the figures label it eta_lab like open_issues_diag.C), the fb
//    bins at y_CM = y - 0.3466.
//
// Inputs : ../skim/rootfile/WTo{MuNu,ElecNu}_pO_PFMet[_<sample>]_hist.root
//          ../plotting/plots[/Elec]/combine_input_W_leppt_mt40.root
//          ./rootfile/qcd_abcd_{mu,ele}.root  (abcd_counts_*: cross-check only)
// Outputs: plots/qcd_ysplit/  (figures + qcd_ysplit.csv)
// Run via ./run_qcd_ysplit.sh -> logs/qcd_ysplit.log (the record).

#include "../plotting/plotting_helper.C"
#include "../skim/skim_common.h" // kYEdges/kYEdgesFB, the electron iso axis
#include "../skim/mc_norm.h"     // pONorm::MCScale

#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TSystem.h"
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

const int kNYs = pOSkim::kNY;
const char *kYsPlotDir = "plots/qcd_ysplit";
const double kYsMtLow = 30.0;  // B/D: m_T < 30   (qcd_abcd.C cfg.yCut)
const double kYsMtSR  = 40.0;  // C40/SR: m_T > 40 (qcd_abcd.C kSRMtCut)

struct YsCnt { double n = 0, v = 0; }; // value, variance
double YsErr(const YsCnt &c) { return std::sqrt(std::max(c.v, 0.0)); }

struct YsSample { std::string label, suffix; bool isW; TFile *f; double k; };

struct YsInputs {
  bool isElec = false;
  std::string lep, fileBase;
  double isoPassLo = 0, isoCut = 0, isoFailLo = 0, isoFailHi = 0;
  TFile *data = nullptr, *ci = nullptr, *abcd = nullptr;
  std::vector<YsSample> mc;
};

// Integral + variance over [xlo,xhi) x [ylo,yhi): qcd_abcd.C count2D's
// convention (no x under/overflow; a yhi beyond the axis includes the overflow).
YsCnt YsCount2D(TH2 *h, double xlo, double xhi, double ylo, double yhi)
{
  YsCnt c;
  if (!h) return c;
  const TAxis *ax = h->GetXaxis(), *ay = h->GetYaxis();
  const double eps = 1e-6;
  const int bxlo = (xlo <= ax->GetXmin()) ? 1 : ax->FindBin(xlo + eps);
  const int bxhi = (xhi >= ax->GetXmax()) ? ax->GetNbins() : ax->FindBin(xhi - eps);
  const int bylo = (ylo <= ay->GetXmin()) ? 1 : ay->FindBin(ylo + eps);
  const int byhi = (yhi >= ay->GetXmax()) ? ay->GetNbins() + 1 : ay->FindBin(yhi - eps);
  double e = 0;
  c.n = h->IntegralAndError(bxlo, bxhi, bylo, byhi, e);
  c.v = e * e;
  return c;
}

// one ABCD region of one plane: data, W-related EWK, other EWK (k_s applied)
struct YsReg {
  YsCnt d, w, z;
  bool ok = true;
  YsCnt qcd(double r) const { return {d.n - r * w.n - z.n, d.v + r * r * w.v + z.v}; }
};

YsReg YsReadRegion(const YsInputs &in, const std::string &hname, double xlo, double xhi, double ylo, double yhi)
{
  YsReg r;
  TH2 *hd = (TH2 *)in.data->Get(hname.c_str());
  if (!hd) { printf("[ERROR] %s missing in the data file -- re-skim the W channels\n", hname.c_str()); r.ok = false; return r; }
  r.d = YsCount2D(hd, xlo, xhi, ylo, yhi);
  for (const YsSample &s : in.mc) {
    TH2 *h = (TH2 *)s.f->Get(hname.c_str());
    if (!h) { printf("[ERROR] %s missing in %s\n", hname.c_str(), s.suffix.c_str()); r.ok = false; continue; }
    YsCnt c = YsCount2D(h, xlo, xhi, ylo, yhi);
    YsCnt &t = s.isW ? r.w : r.z;
    t.n += s.k * c.n;
    t.v += s.k * s.k * c.v;
  }
  return r;
}

// the three regions of one plane
struct YsABCD { YsReg B, C, D; bool ok = true; };

YsABCD YsReadPlane(const YsInputs &in, const std::string &hname)
{
  YsABCD p;
  p.B = YsReadRegion(in, hname, in.isoPassLo, in.isoCut, 0.0, kYsMtLow);
  p.C = YsReadRegion(in, hname, in.isoFailLo, in.isoFailHi, kYsMtSR, 1e9);
  p.D = YsReadRegion(in, hname, in.isoFailLo, in.isoFailHi, 0.0, kYsMtLow);
  p.ok = p.B.ok && p.C.ok && p.D.ok;
  return p;
}

// A = B*C/D with relative errors in quadrature (as qcd_abcd.C); ok=false when
// any QCD-only count is <= 0
struct YsPred { double a = 0, e = 0, T = 0, eT = 0, R = 0, eR = 0; bool ok = false; };

YsPred YsPredict(const YsABCD &p, double r)
{
  YsPred o;
  const YsCnt B = p.B.qcd(r), C = p.C.qcd(r), D = p.D.qcd(r);
  if (B.n <= 0 || C.n <= 0 || D.n <= 0) return o;
  const double rb = YsErr(B) / B.n, rc = YsErr(C) / C.n, rd = YsErr(D) / D.n;
  o.a = B.n * C.n / D.n;
  o.e = o.a * std::sqrt(rb * rb + rc * rc + rd * rd);
  o.T = B.n / D.n;  o.eT = o.T * std::sqrt(rb * rb + rd * rd);
  o.R = C.n / D.n;  o.eR = o.R * std::sqrt(rc * rc + rd * rd);
  o.ok = true;
  return o;
}

double YsInteg(TFile *f, const std::string &path)
{
  TH1 *h = f ? (TH1 *)f->Get(path.c_str()) : nullptr;
  return h ? h->Integral() : 0.0;
}

double YsAbcdLabel(TFile *f, const std::string &hname, const char *label)
{
  TH1 *h = f ? (TH1 *)f->Get(hname.c_str()) : nullptr;
  if (!h) return -1;
  for (int i = 1; i <= h->GetNbinsX(); ++i)
    if (std::string(h->GetXaxis()->GetBinLabel(i)) == label) return h->GetBinContent(i);
  return -1;
}

bool YsOpen(YsInputs &in, bool isElec)
{
  in.isElec = isElec;
  in.lep = isElec ? "ele" : "mu";
  in.fileBase = isElec ? "WToElecNu_pO_PFMet" : "WToMuNu_pO_PFMet";
  if (isElec) { // qcd_abcd.C electron config (iso axis of skim_common.h)
    in.isoPassLo = pOSkim::kEleIsoAxisLo; in.isoCut = 0.0;
    in.isoFailLo = pOSkim::kEleAntiIsoLo; in.isoFailHi = pOSkim::kEleAntiIsoHi;
  } else {      // qcd_abcd.C muon defaults
    in.isoPassLo = 0.0; in.isoCut = 0.15; in.isoFailLo = 0.30; in.isoFailHi = 1.00;
  }
  const std::string dir = "../skim/rootfile/";
  in.data = TFile::Open((dir + in.fileBase + "_hist.root").c_str());
  if (!in.data || in.data->IsZombie()) { printf("[ERROR] cannot open the %s data skim\n", in.lep.c_str()); return false; }
  const std::string dyLab = isElec ? "DYee" : "DYmu";
  struct Spec { std::string label, suffix; bool isW; };
  const std::vector<Spec> spec = {{"Wp_" + in.lep, "Wp", true}, {"Wm_" + in.lep, "Wm", true},
                                  {"Wp_tau", "Wptau", true},    {"Wm_tau", "Wmtau", true},
                                  {dyLab, "DY", false},         {"DYtau", "DYtau", false},
                                  {"TT", "TT", false}};
  for (const Spec &s : spec) {
    TFile *f = TFile::Open((dir + in.fileBase + "_" + s.suffix + "_hist.root").c_str());
    if (!f || f->IsZombie()) { printf("[ERROR] cannot open %s MC %s\n", in.lep.c_str(), s.suffix.c_str()); return false; }
    in.mc.push_back({s.label, s.suffix, s.isW, f, pONorm::MCScale(s.label)});
  }
  in.ci = TFile::Open(Form("../plotting/plots/%scombine_input_W_leppt_mt40.root", isElec ? "Elec/" : ""));
  if (!in.ci || in.ci->IsZombie()) { printf("[ERROR] cannot open the %s leppt_mt40 Combine input\n", in.lep.c_str()); return false; }
  in.abcd = TFile::Open(Form("./rootfile/qcd_abcd_%s.root", in.lep.c_str()));
  return true;
}

// per-y results of one (flavour, charge, binning)
struct YsCurve {
  double x[12] = {0}, ex[12] = {0};
  YsPred p1[12], pR[12];          // per-y ABCD at r = 1 and r = r_W
  double tmpl[12] = {0};          // the fit's template per y
  double c40s[12] = {0};          // A0 * C40_y / Sum C40 (r = 1)
  double wLow1[12] = {0}, wLowR[12] = {0}; // low-MET weights (r = 1 replica, r = r_W)
  double qreq[12] = {0}, eqreq[12] = {0};  // required multijet at r_W (Combine input)
  double A0 = 0, A0R = 0, T0 = 0, R0 = 0, sumT = 0, sumA1 = 0, sumAR = 0, kFF = 1;
  bool ok = false;
};

// the low-MET excess weights of mtandmet.C qcdWeights (data - r*W - other at
// MET < 30, per bin, clamped >= 0, normalized to 1); returns the unnormalized sum
double YsLowMetWeights(const YsInputs &in, const char *chgW, const char *suf, double r, double *w)
{
  double sum = 0;
  for (int i = 0; i < kNYs; ++i) {
    const std::string hn = Form("h_met_%s_y%d%s", chgW, i, suf);
    auto low = [&](TFile *f) -> double {
      TH1 *h = (TH1 *)f->Get(hn.c_str());
      if (!h) return 0.0;
      return h->Integral(1, h->GetXaxis()->FindBin(30.0 - 1e-6));
    };
    double e = low(in.data);
    for (const YsSample &s : in.mc) e -= (s.isW ? r : 1.0) * s.k * low(s.f);
    if (e < 0) e = 0;
    w[i] = e;
    sum += e;
  }
  for (int i = 0; i < kNYs; ++i) w[i] = sum > 0 ? w[i] / sum : 1.0 / kNYs;
  return sum;
}

FILE *gYsCsv = nullptr;

YsCurve YsAnalyze(const YsInputs &in, int ic, bool fb, double rW)
{
  YsCurve cv;
  const char *chgW = ic == 0 ? "Wp" : "Wm";
  const char *chgL = ic == 0 ? "Plus" : "Minus";
  const char *bn = fb ? "fb" : "lab";
  const char *suf = fb ? "_FB" : "";
  const double *edges = fb ? pOSkim::kYEdgesFB : pOSkim::kYEdges;
  const double shift = fb ? 0.3466 : 0.0; // fb bins drawn at y_CM

  // inclusive plane (the qcd_abcd.C normalization) and the per-y planes
  const YsABCD incl = YsReadPlane(in, Form("h_iso_mt_%s%s", in.lep.c_str(), chgL));
  YsABCD py[12];
  bool ok = incl.ok;
  for (int i = 0; i < kNYs; ++i) {
    py[i] = YsReadPlane(in, Form("h_iso_mt_%s%s_y%d%s", in.lep.c_str(), chgL, i, suf));
    ok = ok && py[i].ok;
  }
  if (!ok) return cv;
  const YsPred P1 = YsPredict(incl, 1.0), PR = YsPredict(incl, rW);
  cv.A0 = P1.a; cv.A0R = PR.a; cv.T0 = P1.T; cv.R0 = P1.R;

  printf("\n--- %s %s, %s binning ------------------------------------------------\n", in.lep.c_str(),
         chgW, bn);
  // closure: the y bins against the inclusive plane (data, every region)
  YsCnt sB, sC, sD;
  for (int i = 0; i < kNYs; ++i) { sB.n += py[i].B.d.n; sC.n += py[i].C.d.n; sD.n += py[i].D.d.n; }
  printf("  [closure] data Sum_y vs inclusive plane: B %.0f / %.0f  C40 %.0f / %.0f  D %.0f / %.0f  %s\n", sB.n,
         incl.B.d.n, sC.n, incl.C.d.n, sD.n, incl.D.d.n,
         fb ? "(fb window [-1.7068, 2.4] covers part of the lab range)"
            : ((sB.n == incl.B.d.n && sC.n == incl.C.d.n && sD.n == incl.D.d.n) ? "EXACT" : "MISMATCH"));
  // cross-check with the exported qcd_abcd.C counts (r = 1)
  const std::string ac = Form("abcd_counts_%s%s", in.lep.c_str(), chgL);
  const double xB = YsAbcdLabel(in.abcd, ac, "qcdB0"), xC = YsAbcdLabel(in.abcd, ac, "qcdC0"),
               xD = YsAbcdLabel(in.abcd, ac, "qcdD0");
  if (!fb && xB > 0)
    printf("  [check] QCD-only B/C40/D (r=1) here %.3f / %.3f / %.3f vs qcd_abcd.C %.3f / %.3f / %.3f\n",
           incl.B.qcd(1).n, incl.C.qcd(1).n, incl.D.qcd(1).n, xB, xC, xD);

  // the fit's template and the per-y quantities
  double sumC40 = 0;
  for (int i = 0; i < kNYs; ++i) {
    const std::string dir = Form("%s_%s_y%d/", chgW, bn, i);
    cv.tmpl[i] = YsInteg(in.ci, dir + "qcd_abcd");
    cv.sumT += cv.tmpl[i];
    cv.p1[i] = YsPredict(py[i], 1.0);
    cv.pR[i] = YsPredict(py[i], rW);
    if (cv.p1[i].ok) cv.sumA1 += cv.p1[i].a;
    if (cv.pR[i].ok) cv.sumAR += cv.pR[i].a;
    // C40 summed over the LAB bins (the whole acceptance), also for fb, so the
    // anti-iso split gives the FB window its share of A0, not all of it
    sumC40 += std::max(YsReadPlane(in, Form("h_iso_mt_%s%s_y%d", in.lep.c_str(), chgL, i)).C.qcd(1).n, 0.0);
    const double D = YsInteg(in.ci, dir + "data_obs");
    const double W = YsInteg(in.ci, dir + "signal") + YsInteg(in.ci, dir + "wtau");
    const double Z = YsInteg(in.ci, dir + "z") + YsInteg(in.ci, dir + "ztau") + YsInteg(in.ci, dir + "ttbar");
    cv.qreq[i] = D - rW * W - Z;
    cv.eqreq[i] = std::sqrt(D);
    cv.x[i] = 0.5 * (edges[i] + edges[i + 1]) - shift;
    cv.ex[i] = 0.5 * (edges[i + 1] - edges[i]);
  }
  for (int i = 0; i < kNYs; ++i) cv.c40s[i] = sumC40 > 0 ? cv.A0 * std::max(py[i].C.qcd(1).n, 0.0) / sumC40 : 0;
  const double lowSum = YsLowMetWeights(in, chgW, suf, 1.0, cv.wLow1);
  YsLowMetWeights(in, chgW, suf, rW, cv.wLowR);
  if (fb) {
    // The FB window [-1.7068, 2.4] does not cover y < -1.7068, but mtandmet.C
    // normalizes the FB weights to 1 over the FB bins, so the fb template sums
    // to the FULL y-inclusive total. Fraction of the multijet that is actually
    // inside the window, two ways: the low-MET excess, and the per-y ABCD.
    double wl[12], sumLab = 0;
    const double lowSumLab = YsLowMetWeights(in, chgW, "", 1.0, wl);
    for (int i = 0; i < kNYs; ++i) {
      const YsPred q = YsPredict(YsReadPlane(in, Form("h_iso_mt_%s%s_y%d", in.lep.c_str(), chgL, i)), 1.0);
      if (q.ok) sumLab += q.a;
    }
    double sumFB = 0;
    for (int i = 0; i < kNYs; ++i) {
      const YsPred q = YsPredict(py[i], 1.0);
      if (q.ok) sumFB += q.a;
    }
    double tLab = 0, dLab = 0, dFB = 0;
    for (int i = 0; i < kNYs; ++i) {
      tLab += YsInteg(in.ci, Form("%s_lab_y%d/qcd_abcd", chgW, i));
      dLab += YsInteg(in.ci, Form("%s_lab_y%d/data_obs", chgW, i));
      dFB  += YsInteg(in.ci, Form("%s_fb_y%d/data_obs", chgW, i));
    }
    double tFB = 0;
    for (int i = 0; i < kNYs; ++i) tFB += YsInteg(in.ci, Form("%s_fb_y%d/qcd_abcd", chgW, i));
    const double tShare = tLab > 0 ? tFB / tLab : 0;
    const double lShare = lowSumLab > 0 ? lowSum / lowSumLab : 0, aShare = sumLab > 0 ? sumFB / sumLab : 0;
    printf("  [FB window] Sum qcd_abcd fb/lab = %.4f; inside the window: SR data %.3f, low-MET excess %.3f, "
           "per-y ABCD %.3f of the lab total -> fb template / window share = x%.3f (low-MET) / x%.3f (per-y ABCD)%s\n",
           tShare, dLab > 0 ? dFB / dLab : 0, lShare, aShare, lShare > 0 ? tShare / lShare : 0,
           aShare > 0 ? tShare / aShare : 0,
           (lShare > 0 && std::fabs(tShare / lShare - 1) < 1e-3) ? "  (consistent: the 2026-10-19 normalization)"
                                                                  : "  <-- the fb template does not carry the window's share");
  }
  // the factor the template total carries beyond the flat A0 (the electron
  // fake-factor weighting): always from the LAB template, which holds the
  // whole inclusive total (the fb one holds only the FB window's share)
  double sumTLab = 0;
  for (int i = 0; i < kNYs; ++i) sumTLab += YsInteg(in.ci, Form("%s_lab_y%d/qcd_abcd", chgW, i));
  cv.kFF = cv.A0 > 0 ? sumTLab / cv.A0 : 1.0;
  double maxDev = 0;
  for (int i = 0; i < kNYs; ++i)
    maxDev = std::max(maxDev, std::fabs(cv.wLow1[i] - (cv.sumT > 0 ? cv.tmpl[i] / cv.sumT : 0)));
  printf("  [template] Sum_y qcd_abcd = %.2f  (flat A0 = B0*C40/D0 = %.2f at r=1, %.2f at r=%.2f; "
         "lab template/A0 = %.4f%s)\n", cv.sumT, cv.A0, cv.A0R, rW, cv.kFF,
         in.isElec ? " = the fake-factor weighting" : "");
  printf("  [template] low-MET weights replica vs template split: max |dw| = %.1e%s\n", maxDev,
         maxDev < 1e-6 ? " (identical)" : "  <-- template NOT the low-MET split?");
  printf("  [inclusive] T = B/D = %.4f +- %.4f   R = C40/D = %.4f +- %.4f   (r = 1)\n", P1.T, P1.eT, P1.R, P1.eR);

  printf("  bin    y   | QCD-only (r=1):  B      C40      D   |   T_y/T        R_y/R      | ABCD_y (r=1)    tmpl_y "
         " C40s_y | abs ratio      shape ratio  shape(r=%.1f)| w_low  w_low(r=%.1f) | Q_req(r=%.1f)\n", rW, rW, rW);
  double chiS = 0, chiR = 0; int nS = 0, nR = 0;
  for (int i = 0; i < kNYs; ++i) {
    const YsCnt B = py[i].B.qcd(1), C = py[i].C.qcd(1), D = py[i].D.qcd(1);
    const YsPred &p = cv.p1[i], &q = cv.pR[i];
    double rAbs = 0, eAbs = 0, rSh = 0, eSh = 0, rShR = 0, eShR = 0;
    if (p.ok && cv.tmpl[i] > 0) {
      rAbs = p.a * cv.kFF / cv.tmpl[i];  eAbs = rAbs * p.e / p.a;
      rSh = (p.a / cv.sumA1) / (cv.tmpl[i] / cv.sumT);  eSh = rSh * p.e / p.a;
      chiS += std::pow((rSh - 1) / eSh, 2); ++nS;
    }
    if (q.ok && cv.tmpl[i] > 0) {
      rShR = (q.a / cv.sumAR) / (cv.tmpl[i] / cv.sumT);  eShR = rShR * q.e / q.a;
    }
    if (p.ok) { chiR += std::pow((p.R - P1.R) / p.eR, 2); ++nR; }
    printf("  y%-2d %6.3f | %8.1f %8.1f %8.1f | %5.3f+-%5.3f %5.3f+-%5.3f | %6.1f+-%5.1f %7.1f %7.1f | "
           "%5.3f+-%5.3f %5.3f+-%5.3f %5.3f+-%5.3f | %5.3f %5.3f | %7.1f+-%4.1f\n",
           i, cv.x[i], B.n, C.n, D.n, p.ok ? p.T / P1.T : 0, p.ok ? p.eT / P1.T : 0, p.ok ? p.R / P1.R : 0,
           p.ok ? p.eR / P1.R : 0, p.a, p.e, cv.tmpl[i], cv.c40s[i], rAbs, eAbs, rSh, eSh, rShR, eShR, cv.wLow1[i],
           cv.wLowR[i], cv.qreq[i], cv.eqreq[i]);
    if (gYsCsv)
      fprintf(gYsCsv, "%s,%s,%s,%d,%.4f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.5f,%.5f,%.5f,%.5f,%.3f,%.3f,%.3f,%.3f,%.3f,"
                      "%.3f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.3f,%.3f\n",
              in.lep.c_str(), chgW, bn, i, cv.x[i], B.n, YsErr(B), C.n, YsErr(C), D.n, YsErr(D), p.T, p.eT, p.R,
              p.eR, p.a, p.e, q.a, q.e, cv.tmpl[i], cv.c40s[i], rAbs, eAbs, rSh, eSh, rShR, eShR, cv.wLow1[i],
              cv.wLowR[i], cv.qreq[i], cv.eqreq[i]);
  }
  printf("  chi2 vs 1 of the shape ratio (r=1): %.1f / %d    chi2 of R_y vs inclusive R: %.1f / %d  "
         "(Gaussian errors, indicative)\n", chiS, nS > 0 ? nS - 1 : 0, chiR, nR > 0 ? nR - 1 : 0);
  {
    double ea2 = 0;
    for (int i = 0; i < kNYs; ++i) ea2 += cv.p1[i].e * cv.p1[i].e;
    printf("  [sum] Sum_y ABCD_y = %.1f +- %.1f = %.3f x the inclusive A0 (r=1); %.1f = %.3f x A0 at r=%.1f\n",
           cv.sumA1, std::sqrt(ea2), cv.A0 > 0 ? cv.sumA1 / cv.A0 : 0, cv.sumAR, cv.A0R > 0 ? cv.sumAR / cv.A0R : 0, rW);
  }

  // the outermost bins y0 + y11 (where the electron W- deficit sits, lab)
  {
    double a = 0, ea2 = 0, aR = 0, eaR2 = 0, t = 0, qr = 0, eqr2 = 0, c40s = 0;
    YsReg Bs, Cs, Ds; // regions summed over the two bins (ABCD of the sum)
    for (int i : {0, kNYs - 1}) {
      a += cv.p1[i].a; ea2 += cv.p1[i].e * cv.p1[i].e;
      aR += cv.pR[i].a; eaR2 += cv.pR[i].e * cv.pR[i].e;
      t += cv.tmpl[i]; qr += cv.qreq[i]; eqr2 += cv.eqreq[i] * cv.eqreq[i]; c40s += cv.c40s[i];
      auto add = [](YsReg &s, const YsReg &r) {
        s.d.n += r.d.n; s.d.v += r.d.v; s.w.n += r.w.n; s.w.v += r.w.v; s.z.n += r.z.n; s.z.v += r.z.v;
      };
      add(Bs, py[i].B); add(Cs, py[i].C); add(Ds, py[i].D);
    }
    YsABCD es; es.B = Bs; es.C = Cs; es.D = Ds;
    const YsPred pe = YsPredict(es, 1.0);
    const double rAbs = a * cv.kFF / t, eAbs = std::sqrt(ea2) * cv.kFF / t;
    const double rSh = (a / cv.sumA1) / (t / cv.sumT), eSh = rSh * std::sqrt(ea2) / a;
    const double rShR = (aR / cv.sumAR) / (t / cv.sumT), eShR = rShR * std::sqrt(eaR2) / aR;
    printf("  y0+y11: ABCD (sum of bins) %.1f +- %.1f, ABCD of the summed regions %.1f +- %.1f, template %.1f, "
           "C40s %.1f, Q_req(r=%.1f) %.1f +- %.1f\n", a, std::sqrt(ea2), pe.a, pe.e, t, c40s, rW, qr, std::sqrt(eqr2));
    printf("  y0+y11: abs ratio ABCD/template = %.3f +- %.3f   shape ratio = %.3f +- %.3f (r=1), %.3f +- %.3f (r=%.1f)"
           "   T_edge/T = %.3f +- %.3f   R_edge/R = %.3f +- %.3f   Q_req/template = %.3f +- %.3f\n",
           rAbs, eAbs, rSh, eSh, rShR, eShR, rW, pe.T / P1.T, pe.eT / P1.T, pe.R / P1.R, pe.eR / P1.R, qr / t,
           std::sqrt(eqr2) / t);
    if (gYsCsv)
      fprintf(gYsCsv, "%s,%s,%s,edge,,,,,,,,%.5f,%.5f,%.5f,%.5f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.5f,%.5f,%.5f,%.5f,"
                      "%.5f,%.5f,,,%.3f,%.3f\n",
              in.lep.c_str(), chgW, bn, pe.T, pe.eT, pe.R, pe.eR, a, std::sqrt(ea2), aR, std::sqrt(eaR2), t, c40s, rAbs,
              eAbs, rSh, eSh, rShR, eShR, qr, std::sqrt(eqr2));
  }
  cv.ok = true;
  return cv;
}

// ---- figures ---------------------------------------------------------------
PlotStyle YsStyle()
{
  PlotStyle ps;
  ps.headerX = 0.47;
  ps.headerY = 0.94;
  ps.subSize = 0.032;
  ps.boxTextSize = 0.026;
  ps.boxX1 = 0.47; ps.boxX2 = 0.93; ps.boxY1 = 0.71; ps.boxY2 = 0.77;
  return ps;
}

TGraphErrors *YsGraph(const YsCurve &cv, const double *y, const double *ey)
{
  TGraphErrors *g = new TGraphErrors();
  for (int i = 0; i < kNYs; ++i) {
    const int n = g->GetN();
    g->SetPoint(n, cv.x[i], y[i]);
    g->SetPointError(n, cv.ex[i], ey ? ey[i] : 0.0);
  }
  return g;
}

const char *YsXTitle(bool fb)
{
  return fb ? "#eta_{CM} (proton-going > 0)" : "#eta_{lab} (proton-going > 0)";
}

// SR multijet yield per rapidity bin: template vs per-y ABCD vs anti-iso split
void YsPlotYields(const YsInputs &in, const YsCurve &cv, int ic, bool fb, double rW)
{
  const char *chgW = ic == 0 ? "Wp" : "Wm";
  double a1[12], e1[12], aR[12], eR[12];
  double ymax = 0;
  for (int i = 0; i < kNYs; ++i) {
    a1[i] = cv.p1[i].a * cv.kFF; e1[i] = cv.p1[i].e * cv.kFF;
    aR[i] = cv.pR[i].a * cv.kFF; eR[i] = cv.pR[i].e * cv.kFF;
    ymax = std::max({ymax, a1[i] + e1[i], aR[i] + eR[i], cv.tmpl[i], cv.c40s[i] * cv.kFF});
  }
  double c40k[12];
  for (int i = 0; i < kNYs; ++i) c40k[i] = cv.c40s[i] * cv.kFF;
  std::vector<OverlaySeries> ser(4);
  ser[0].gTot = YsGraph(cv, cv.tmpl, nullptr);
  ser[0].label = "fit template (low-MET split)";
  ser[0].color = kBlack; ser[0].marker = 20; ser[0].markerSize = 1.3; ser[0].xShift = -0.6;
  ser[1].gTot = YsGraph(cv, a1, e1);
  ser[1].label = "B_{y} C40_{y} / D_{y}  (r = 1)";
  ser[1].color = kRed + 1; ser[1].marker = 21; ser[1].markerSize = 1.5; ser[1].xShift = -0.2;
  ser[2].gTot = YsGraph(cv, aR, eR);
  ser[2].label = Form("B_{y} C40_{y} / D_{y}  (r_{W} = %.1f)", rW);
  ser[2].color = kOrange + 7; ser[2].marker = 25; ser[2].markerSize = 1.5; ser[2].xShift = 0.2;
  ser[3].gTot = YsGraph(cv, c40k, nullptr);
  ser[3].label = "A_{0} #times C40_{y} / #Sigma C40  (T flat in y)";
  ser[3].color = kGreen + 2; ser[3].marker = 33; ser[3].markerSize = 2.0; ser[3].xShift = 0.6;
  PlotStyle ps = YsStyle();
  // data in the lower ~45% of the frame, the legend above them, under the box
  ps.legX1 = 0.17; ps.legY1 = 0.51; ps.legX2 = 0.62; ps.legY2 = 0.69;
  GraphTuner tuner = [ymax](TCanvas *, TGraphErrors *g) { g->SetMinimum(0); g->SetMaximum(2.3 * ymax); };
  std::vector<std::string> box = {in.isElec ? "per-y ABCD scaled by the template's FF factor" : "prefit, absolute MC",
                                  "errors: data + MC stat. of the per-y ABCD"};
  SaveNiceGraph_Overlay(ser, Form("%s/qcd_ysplit_%s_%s_%s", kYsPlotDir, in.lep.c_str(), chgW, fb ? "fb" : "lab"),
                        YsXTitle(fb), "SR multijet events (p_{T} > 25, m_{T} > 40)", "",
                        Form("%s channel, %s, %s bins", in.isElec ? "e" : "#mu", ic == 0 ? "W^{+}" : "W^{-}",
                             fb ? "FB" : "lab"),
                        "multijet: fit template vs per-y ABCD", box, ps, tuner);
}

// (per-y ABCD / template), shape-normalized, W+ and W-
void YsPlotRatio(const YsInputs &in, const YsCurve cv[2], bool fb)
{
  std::vector<OverlaySeries> ser(2);
  for (int c = 0; c < 2; ++c) {
    double r[12], e[12];
    for (int i = 0; i < kNYs; ++i) {
      const YsPred &p = cv[c].p1[i];
      r[i] = (p.ok && cv[c].tmpl[i] > 0) ? (p.a / cv[c].sumA1) / (cv[c].tmpl[i] / cv[c].sumT) : 0;
      e[i] = p.ok ? r[i] * p.e / p.a : 0;
    }
    ser[c].gTot = YsGraph(cv[c], r, e);
    ser[c].label = c == 0 ? "W^{+}" : "W^{-}";
    ser[c].color = c == 0 ? kAzure + 2 : kRed + 1;
    ser[c].marker = c == 0 ? 21 : 33;
    ser[c].markerSize = c == 0 ? 1.5 : 2.0;
    ser[c].xShift = c == 0 ? -0.35 : 0.35;
  }
  PlotStyle ps = YsStyle();
  ps.legX1 = 0.55; ps.legY1 = 0.15; ps.legX2 = 0.93; ps.legY2 = 0.27;
  GraphTuner tuner = [](TCanvas *c, TGraphErrors *g) {
    g->SetMinimum(0.0);
    g->SetMaximum(2.6);
    c->cd();
    TLine *one = new TLine(g->GetXaxis()->GetXmin(), 1.0, g->GetXaxis()->GetXmax(), 1.0);
    one->SetLineStyle(2); one->SetLineColor(kGray + 2); one->Draw();
  };
  SaveNiceGraph_Overlay(ser, Form("%s/qcd_ysplit_ratio_%s_%s", kYsPlotDir, in.lep.c_str(), fb ? "fb" : "lab"),
                        YsXTitle(fb), "per-y ABCD / fit template  (shape)", "",
                        Form("%s channel, %s bins, r = 1", in.isElec ? "e" : "#mu", fb ? "FB" : "lab"),
                        "both normalized to their sum over y",
                        {"errors: per-y ABCD stat. only"}, ps, tuner);
}

// T_y/T and R_y/R per charge
void YsPlotTR(const YsInputs &in, const YsCurve &cv, int ic, bool fb)
{
  double t[12], et[12], r[12], er[12];
  for (int i = 0; i < kNYs; ++i) {
    const YsPred &p = cv.p1[i];
    t[i] = p.ok ? p.T / cv.T0 : 0;  et[i] = p.ok ? p.eT / cv.T0 : 0;
    r[i] = p.ok ? p.R / cv.R0 : 0;  er[i] = p.ok ? p.eR / cv.R0 : 0;
  }
  std::vector<OverlaySeries> ser(2);
  ser[0].gTot = YsGraph(cv, t, et);
  ser[0].label = Form("T_{y}/T,  T = B/D = %.3f", cv.T0);
  ser[0].color = kViolet + 1; ser[0].marker = 20; ser[0].markerSize = 1.3; ser[0].xShift = -0.35;
  ser[1].gTot = YsGraph(cv, r, er);
  ser[1].label = Form("R_{y}/R,  R = C40/D = %.3f", cv.R0);
  ser[1].color = kTeal + 3; ser[1].marker = 21; ser[1].markerSize = 1.5; ser[1].xShift = 0.35;
  PlotStyle ps = YsStyle();
  ps.legX1 = 0.25; ps.legY1 = 0.59; ps.legX2 = 0.68; ps.legY2 = 0.69;
  GraphTuner tuner = [](TCanvas *c, TGraphErrors *g) {
    g->SetMinimum(0.0);
    g->SetMaximum(3.6);
    c->cd();
    TLine *one = new TLine(g->GetXaxis()->GetXmin(), 1.0, g->GetXaxis()->GetXmax(), 1.0);
    one->SetLineStyle(2); one->SetLineColor(kGray + 2); one->Draw();
  };
  SaveNiceGraph_Overlay(ser, Form("%s/qcd_ysplit_TR_%s_%s_%s", kYsPlotDir, in.lep.c_str(), ic == 0 ? "Wp" : "Wm",
                                  fb ? "fb" : "lab"),
                        YsXTitle(fb), "per-y transfer factor / inclusive", "",
                        Form("%s channel, %s, %s bins, r = 1", in.isElec ? "e" : "#mu",
                             ic == 0 ? "W^{+}" : "W^{-}", fb ? "FB" : "lab"),
                        "T: iso (m_{T}<30), R: m_{T} (anti-iso)",
                        {"errors: data + MC stat."}, ps, tuner);
}

// ---- robustness: the anti-iso window --------------------------------------
// The edge (y0+y11, lab) shape ratio and R_edge/R for the three sideband
// windows of qcd_abcd.C's window scan ([lo, 1.0) nominal, [lo, mid), [mid, 1.0)
// with mid = 0.6), r = 1.
void YsWindowScan(const YsInputs &inNom)
{
  printf("\n=== anti-iso window variation: edge bins y0+y11 (lab), r = 1 ===\n");
  printf("  %-4s %-3s | window       | shape ratio ABCD/template | R_edge/R      | T_edge/T\n", "", "");
  const double wins[3][2] = {{inNom.isoFailLo, inNom.isoFailHi}, {inNom.isoFailLo, 0.6}, {0.6, inNom.isoFailHi}};
  for (int c = 0; c < 2; ++c)
    for (int k = 0; k < 3; ++k) {
      YsInputs in = inNom;
      in.isoFailLo = wins[k][0];
      in.isoFailHi = wins[k][1];
      const char *chgL = c == 0 ? "Plus" : "Minus";
      const YsABCD incl = YsReadPlane(in, Form("h_iso_mt_%s%s", in.lep.c_str(), chgL));
      const YsPred P = YsPredict(incl, 1.0);
      double sumA = 0, sumT = 0, aE = 0, eE2 = 0, tE = 0;
      YsABCD edge;
      bool ok = incl.ok && P.ok;
      for (int i = 0; i < kNYs; ++i) {
        const YsABCD p = YsReadPlane(in, Form("h_iso_mt_%s%s_y%d", in.lep.c_str(), chgL, i));
        const YsPred q = YsPredict(p, 1.0);
        const double t = YsInteg(in.ci, Form("%s_lab_y%d/qcd_abcd", c == 0 ? "Wp" : "Wm", i));
        ok = ok && p.ok && q.ok;
        sumA += q.a; sumT += t;
        if (i == 0 || i == kNYs - 1) {
          aE += q.a; eE2 += q.e * q.e; tE += t;
          auto add = [](YsReg &s, const YsReg &r) {
            s.d.n += r.d.n; s.d.v += r.d.v; s.w.n += r.w.n; s.w.v += r.w.v; s.z.n += r.z.n; s.z.v += r.z.v;
          };
          add(edge.B, p.B); add(edge.C, p.C); add(edge.D, p.D);
        }
      }
      if (!ok) { printf("  %-4s %-3s | [%.2f, %.2f) | (a QCD-only count <= 0)\n", in.lep.c_str(), c == 0 ? "Wp" : "Wm", wins[k][0], wins[k][1]); continue; }
      const YsPred pe = YsPredict(edge, 1.0);
      const double sh = (aE / sumA) / (tE / sumT), esh = sh * std::sqrt(eE2) / aE;
      printf("  %-4s %-3s | [%.2f, %.2f) | %5.3f +- %5.3f             | %5.3f +- %5.3f | %5.3f +- %5.3f\n",
             in.lep.c_str(), c == 0 ? "Wp" : "Wm", wins[k][0], wins[k][1], sh, esh, pe.R / P.R, pe.eR / P.R,
             pe.T / P.T, pe.eT / P.T);
    }
}

// ---- consequence for the electron/muon comparison ---------------------------
// open_issues_diag.C's prefit r_Q = (D - B - Q)/W per lab bin (D data, W =
// signal + wtau, B = z + ztau + ttbar, Q = multijet; data Poisson errors only)
// with Q = the fit template, and with Q = the template total redistributed in
// y like the per-y ABCD (Q'_y = Sum tmpl * ABCD_y / Sum ABCD, r = 1): the
// alternative rapidity split at the SAME inclusive multijet total, applied to
// both flavours. Also the "abs" variant Q''_y = ABCD_y * (Sum tmpl / A0).
struct YsRq { double D = 0, W = 0, B = 0, Qt = 0, Qs = 0, Qa = 0; };

YsRq YsReadRq(const YsInputs &in, const YsCurve &cv, int ic, int i)
{
  YsRq r;
  const std::string dir = Form("%s_lab_y%d/", ic == 0 ? "Wp" : "Wm", i);
  r.D = YsInteg(in.ci, dir + "data_obs");
  r.W = YsInteg(in.ci, dir + "signal") + YsInteg(in.ci, dir + "wtau");
  r.B = YsInteg(in.ci, dir + "z") + YsInteg(in.ci, dir + "ztau") + YsInteg(in.ci, dir + "ttbar");
  r.Qt = cv.tmpl[i];
  r.Qs = cv.p1[i].ok ? cv.sumT * cv.p1[i].a / cv.sumA1 : cv.tmpl[i];
  r.Qa = cv.p1[i].ok ? cv.kFF * cv.p1[i].a : cv.tmpl[i];
  return r;
}

void YsAdd(YsRq &s, const YsRq &r)
{
  s.D += r.D; s.W += r.W; s.B += r.B; s.Qt += r.Qt; s.Qs += r.Qs; s.Qa += r.Qa;
}

// e/mu of r for one multijet choice q (0 template, 1 shape split, 2 abs)
void YsEMuRatio(const YsRq &m, const YsRq &e, int q, double &rat, double &err, double &rm, double &re)
{
  auto Q = [q](const YsRq &r) { return q == 0 ? r.Qt : (q == 1 ? r.Qs : r.Qa); };
  rm = (m.D - m.B - Q(m)) / m.W;
  re = (e.D - e.B - Q(e)) / e.W;
  const double erm = std::sqrt(m.D) / m.W, ere = std::sqrt(e.D) / e.W;
  rat = re / rm;
  err = std::fabs(rat) * std::sqrt(std::pow(erm / rm, 2) + std::pow(ere / re, 2));
}

void YsEMu(const YsInputs &inM, const YsInputs &inE, const YsCurve cvM[2], const YsCurve cvE[2])
{
  printf("\n=== e/mu of the prefit r_Q = (D-B-Q)/W per lab bin (open_issues_diag.C convention, data Poisson errors) ===\n");
  printf("    Q = fit template | Q' = template total split in y like the per-y ABCD (r=1) | Q'' = per-y ABCD x FF "
         "factor\n");
  double emu[2][3][12] = {{{0}}}, eemu[2][3][12] = {{{0}}};
  for (int c = 0; c < 2; ++c) {
    printf("\n  %s  bin |  e: D    Q     Q'    Q''  |  r_e(Q)  r_e(Q')  r_mu(Q)  r_mu(Q') | e/mu (Q)       e/mu (Q')"
           "      e/mu (Q'')\n", c == 0 ? "Wp" : "Wm");
    YsRq em, ee, im, ie, bm, be; // edges y0+y11, inclusive, backward y0..y3 (eta_lab < -0.8)
    for (int i = 0; i < kNYs; ++i) {
      const YsRq m = YsReadRq(inM, cvM[c], c, i), e = YsReadRq(inE, cvE[c], c, i);
      YsAdd(im, m); YsAdd(ie, e);
      if (i == 0 || i == kNYs - 1) { YsAdd(em, m); YsAdd(ee, e); }
      if (i <= 3) { YsAdd(bm, m); YsAdd(be, e); }
      double rat[3], err[3], rm[3], re[3];
      for (int q = 0; q < 3; ++q) {
        YsEMuRatio(m, e, q, rat[q], err[q], rm[q], re[q]);
        emu[c][q][i] = rat[q]; eemu[c][q][i] = err[q];
      }
      printf("        y%-2d | %5.0f %5.1f %5.1f %5.1f | %6.3f  %6.3f   %6.3f   %6.3f  | %5.3f+-%5.3f  %5.3f+-%5.3f  "
             "%5.3f+-%5.3f\n", i, e.D, e.Qt, e.Qs, e.Qa, re[0], re[1], rm[0], rm[1], rat[0], err[0], rat[1], err[1],
             rat[2], err[2]);
    }
    for (int k = 0; k < 3; ++k) {
      const YsRq &m = k == 0 ? em : (k == 1 ? bm : im), &e = k == 0 ? ee : (k == 1 ? be : ie);
      double rat[3], err[3], rm[3], re[3];
      for (int q = 0; q < 3; ++q) YsEMuRatio(m, e, q, rat[q], err[q], rm[q], re[q]);
      printf("    %s | %5.0f %5.1f %5.1f %5.1f | %6.3f  %6.3f   %6.3f   %6.3f  | %5.3f+-%5.3f  %5.3f+-%5.3f  "
             "%5.3f+-%5.3f\n", k == 0 ? "y0+y11" : (k == 1 ? " y0..y3" : "  incl"), e.D, e.Qt, e.Qs, e.Qa, re[0], re[1],
             rm[0], rm[1], rat[0], err[0], rat[1], err[1], rat[2], err[2]);
    }
  }
  // figure: e/mu vs eta_lab per charge, template vs per-y ABCD split
  for (int c = 0; c < 2; ++c) {
    std::vector<OverlaySeries> ser(2);
    for (int q = 0; q < 2; ++q) {
      ser[q].gTot = YsGraph(cvE[c], emu[c][q], eemu[c][q]);
      ser[q].label = q == 0 ? "multijet = fit template (low-MET split)" : "multijet split like the per-y ABCD";
      ser[q].color = q == 0 ? kGray + 2 : (c == 0 ? kAzure + 2 : kRed + 1);
      ser[q].marker = q == 0 ? 24 : (c == 0 ? 21 : 33);
      ser[q].markerSize = q == 0 ? 1.3 : (c == 0 ? 1.5 : 2.0);
      ser[q].xShift = q == 0 ? -0.35 : 0.35;
    }
    PlotStyle ps = YsStyle();
    ps.legX1 = 0.17; ps.legY1 = 0.59; ps.legX2 = 0.70; ps.legY2 = 0.69;
    GraphTuner tuner = [](TCanvas *cc, TGraphErrors *g) {
      g->SetMinimum(0.2);
      g->SetMaximum(2.6);
      cc->cd();
      TLine *one = new TLine(g->GetXaxis()->GetXmin(), 1.0, g->GetXaxis()->GetXmax(), 1.0);
      one->SetLineStyle(2); one->SetLineColor(kGray + 2); one->Draw();
    };
    SaveNiceGraph_Overlay(ser, Form("%s/qcd_ysplit_emu_%s_lab", kYsPlotDir, c == 0 ? "Wp" : "Wm"), YsXTitle(false),
                          "e / #mu  of  (D - B - Q) / W", "",
                          Form("%s, lepton p_{T} at m_{T} > 40, prefit", c == 0 ? "W^{+}" : "W^{-}"),
                          "same inclusive multijet total in both",
                          {"errors: data statistical only"}, ps, tuner);
  }
}

} // namespace

void qcd_ysplit(double rW = 1.2)
{
  gROOT->SetBatch(kTRUE);
  gSystem->mkdir(kYsPlotDir, kTRUE);
  gYsCsv = std::fopen(Form("%s/qcd_ysplit.csv", kYsPlotDir), "w");
  if (gYsCsv)
    std::fprintf(gYsCsv, "flavour,charge,binning,ybin,y,B,B_err,C40,C40_err,D,D_err,T_y,T_y_err,R_y,R_y_err,"
                         "abcd_r1,abcd_r1_err,abcd_rW,abcd_rW_err,template,c40split,abs_ratio,abs_ratio_err,"
                         "shape_ratio,shape_ratio_err,shape_ratio_rW,shape_ratio_rW_err,w_low_r1,w_low_rW,"
                         "Q_req_rW,Q_req_rW_err\n");
  printf("=== qcd_ysplit: SR multijet per rapidity bin, fit template vs per-y ABCD (r_W = %.2f) ===\n", rW);
  printf("    regions (qcd_abcd.C): B iso-pass m_T<30, C40 anti-iso m_T>40, D anti-iso m_T<30\n");
  printf("    template = Combine input <C>_<bins>_y<i>/qcd_abcd (leppt_mt40); counts PREFIT, absolute MC\n");
  YsInputs in[2];
  YsCurve labCv[2][2]; // [flavour][charge], lab binning
  bool ok[2] = {false, false};
  for (int f = 0; f < 2; ++f) {
    if (!YsOpen(in[f], f == 1)) continue;
    printf("\n==================== %s : iso-pass [%g, %g), anti-iso [%g, %g) ====================\n",
           in[f].lep.c_str(), in[f].isoPassLo, in[f].isoCut, in[f].isoFailLo, in[f].isoFailHi);
    ok[f] = true;
    for (int b = 0; b < 2; ++b) {
      YsCurve cv[2];
      for (int c = 0; c < 2; ++c) {
        cv[c] = YsAnalyze(in[f], c, b == 1, rW);
        if (!cv[c].ok) { ok[f] = false; continue; }
        YsPlotYields(in[f], cv[c], c, b == 1, rW);
        YsPlotTR(in[f], cv[c], c, b == 1);
        if (b == 0) labCv[f][c] = cv[c];
      }
      if (cv[0].ok && cv[1].ok) YsPlotRatio(in[f], cv, b == 1);
    }
  }
  for (int f = 0; f < 2; ++f)
    if (ok[f]) YsWindowScan(in[f]);
  if (ok[0] && ok[1]) YsEMu(in[0], in[1], labCv[0], labCv[1]);
  if (gYsCsv) std::fclose(gYsCsv);
  printf("\n[qcd_ysplit] figures + qcd_ysplit.csv in %s/\n", kYsPlotDir);
}
