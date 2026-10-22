// correction/qcd_ysplit_pergroup.C
//
// BEFORE vs AFTER of the rapidity-resolved in-fit ABCD (QCD_MODE=abcdy,
// 2026-10-19): the SR multijet template of the leppt_mt40 fit, per rapidity
// GROUP (correction/abcd_ygroups.h: fwd = y10,y11 | bwd = y0,y1 | cen =
// y2..y9, the same bin indices lab and fb), per flavour x charge, lab and fb.
//
//   OLD = qcd_abcd   the y-inclusive in-fit ABCD (total A0 = B0*C40/D0, the
//                    electron one fake-factor weighted, total A0ff) split over
//                    the rapidity bins by the low-MET excess (mtandmet.C
//                    qcdWeights); one transfer factor per flavour x charge
//   NEW = qcd_abcdy  each group at ITS OWN ABCD total A_g = B_g*C40_g/D_g
//                    (x A0ff/A0 for the electron: the fake-factor weighting
//                    stays inclusive), split inside the group by the same
//                    low-MET proxy; in the fit each group floats with its own
//                    control regions (sB_g*sC_g/sD_g)
//
// Printed / drawn: per group and per bin the old and new yields with the group
// ABCD statistical error (data + MC stat of B_g, C40_g, D_g; fully correlated
// within a group), the group transfer factors T_g = B_g/D_g (iso, m_T < 30)
// and R_g = C40_g/D_g (m_T, anti-iso) against the inclusive T and R, and the
// e/mu ratio of the prefit r_Q = (D - B - Q)/W with Q = old vs new
// (open_issues_diag.C's convention: D data, W = signal + wtau, B = z + ztau +
// ttbar, data Poisson errors only) per bin and for the edge (y0+y11), the
// backward (y0..y3), the three groups and the inclusive sums.
// Self-checks (printed): each group's qcd_abcdy total = A_g x (template/A0);
// the group CR dirs of the Combine input hold exactly the exported counts;
// the lab groups' data add up to the inclusive control regions.
// All numbers PREFIT (r = 1 EWK subtraction), from stored outputs only.
//
// Inputs : ../plotting/plots[/Elec]/combine_input_W_leppt_mt40.root
//          ./rootfile/qcd_abcd_{mu,ele}.root   (abcd_counts_*, incl. the groups)
// Outputs: plots/qcd_ysplit_pergroup/  (figures + qcd_ysplit_pergroup.csv)
// Run via ./run_qcd_ysplit_pergroup.sh -> logs/qcd_ysplit_pergroup.log (the record).

#include "../plotting/plotting_helper.C"
#include "../skim/skim_common.h" // kYEdges / kYEdgesFB
#include "abcd_ygroups.h"

#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLine.h"
#include "TROOT.h"
#include "TSystem.h"
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

static_assert(pOAbcdY::kNYBins == pOSkim::kNY, "abcd_ygroups.h assumes the 12 rapidity bins of skim_common.h");

namespace {

const int kNPg = pOSkim::kNY;
const int kNG = pOAbcdY::kNGroups;
const char *kPgDir = "plots/qcd_ysplit_pergroup";
const double kPgShift = 0.3466; // fb bins drawn at y_CM
FILE *gPgCsv = nullptr;

double PgInteg(TFile *f, const std::string &path)
{
  TH1 *h = f ? (TH1 *)f->Get(path.c_str()) : nullptr;
  return h ? h->Integral() : 0.0;
}

bool PgHas(TFile *f, const std::string &path) { return f && f->Get(path.c_str()); }

// labeled-bin reader (loop + compare: TAxis::FindBin(label) would append)
double PgLabel(TFile *f, const std::string &hname, const char *label, double *err = nullptr, bool *ok = nullptr)
{
  TH1 *h = f ? (TH1 *)f->Get(hname.c_str()) : nullptr;
  if (h)
    for (int i = 1; i <= h->GetNbinsX(); ++i)
      if (std::string(h->GetXaxis()->GetBinLabel(i)) == label)
      {
        if (err) *err = h->GetBinError(i);
        return h->GetBinContent(i);
      }
  if (ok) *ok = false;
  return 0.0;
}

// in-fit ABCD counts of one abcd_counts_* object
struct PgCounts {
  double B = 0, eB = 0, C = 0, eC = 0, D = 0, eD = 0, dB = 0, dC = 0, dD = 0, wB = 0;
  bool ok = true;
  double A() const { return D > 0 ? B * C / D : 0.0; }
  double rA() const { return (B > 0 && C > 0 && D > 0) ? std::sqrt(eB * eB / (B * B) + eC * eC / (C * C) + eD * eD / (D * D)) : 0.0; }
  double T() const { return D > 0 ? B / D : 0.0; }
  double rT() const { return (B > 0 && D > 0) ? std::sqrt(eB * eB / (B * B) + eD * eD / (D * D)) : 0.0; }
  double R() const { return D > 0 ? C / D : 0.0; }
  double rR() const { return (C > 0 && D > 0) ? std::sqrt(eC * eC / (C * C) + eD * eD / (D * D)) : 0.0; }
};

PgCounts PgReadCounts(TFile *f, const std::string &nm)
{
  PgCounts c;
  c.B = PgLabel(f, nm, "qcdB0", &c.eB, &c.ok);
  c.C = PgLabel(f, nm, "qcdC0", &c.eC, &c.ok);
  c.D = PgLabel(f, nm, "qcdD0", &c.eD, &c.ok);
  c.dB = PgLabel(f, nm, "dataB", nullptr, &c.ok);
  c.dC = PgLabel(f, nm, "dataC40", nullptr, &c.ok);
  c.dD = PgLabel(f, nm, "dataD", nullptr, &c.ok);
  c.wB = PgLabel(f, nm, "wB", nullptr, &c.ok);
  if (!c.ok) printf("[ERROR] %s missing or incomplete -- re-run ./run_qcd_abcd.sh\n", nm.c_str());
  return c;
}

struct PgInputs {
  bool isElec = false;
  std::string lep;
  TFile *ci = nullptr, *abcd = nullptr;
};

bool PgOpen(PgInputs &in, bool isElec)
{
  in.isElec = isElec;
  in.lep = isElec ? "ele" : "mu";
  in.ci = TFile::Open(Form("../plotting/plots/%scombine_input_W_leppt_mt40.root", isElec ? "Elec/" : ""));
  in.abcd = TFile::Open(Form("./rootfile/qcd_abcd_%s.root", in.lep.c_str()));
  if (!in.ci || in.ci->IsZombie() || !in.abcd || in.abcd->IsZombie())
  {
    printf("[ERROR] cannot open the %s leppt_mt40 Combine input or qcd_abcd_%s.root\n", in.lep.c_str(), in.lep.c_str());
    return false;
  }
  return true;
}

// one flavour x charge x binning
struct PgSet {
  double x[12] = {0}, ex[12] = {0};
  double oldY[12] = {0}, newY[12] = {0}, eNewY[12] = {0};
  double gOld[3] = {0}, gNew[3] = {0}, eGNew[3] = {0}, gx[3] = {0}, gex[3] = {0};
  PgCounts g[3], incl;
  double kFF = 1.0;
  bool ok = false;
};

std::string PgBins(int g)
{
  const std::vector<int> b = pOAbcdY::BinsOf(g);
  return Form("y%d-y%d", b.front(), b.back());
}

PgSet PgAnalyze(const PgInputs &in, int ic, bool fb)
{
  PgSet s;
  const char *chgW = ic == 0 ? "Wp" : "Wm";
  const char *chgL = ic == 0 ? "Plus" : "Minus";
  const char *bn = fb ? "fb" : "lab";
  const double *edges = fb ? pOSkim::kYEdgesFB : pOSkim::kYEdges;
  const double shift = fb ? kPgShift : 0.0;
  printf("\n--- %s %s, %s binning ------------------------------------------------\n", in.lep.c_str(), chgW, bn);

  s.incl = PgReadCounts(in.abcd, Form("abcd_counts_%s%s", in.lep.c_str(), chgL));
  bool ok = s.incl.ok;
  for (int g = 0; g < kNG; ++g)
  {
    s.g[g] = PgReadCounts(in.abcd, Form("abcd_counts_%s%s_%s_%s", in.lep.c_str(), chgL, bn, pOAbcdY::kGroupName[g]));
    ok = ok && s.g[g].ok;
  }
  for (int i = 0; i < kNPg; ++i)
  {
    const std::string dir = Form("%s_%s_y%d/", chgW, bn, i);
    if (!PgHas(in.ci, dir + "qcd_abcdy"))
    {
      printf("[ERROR] %s missing in the %s Combine input -- run ../plotting/run_combine_inputs.sh\n",
             (dir + "qcd_abcdy").c_str(), in.lep.c_str());
      ok = false;
    }
    s.oldY[i] = PgInteg(in.ci, dir + "qcd_abcd");
    s.newY[i] = PgInteg(in.ci, dir + "qcd_abcdy");
    s.x[i] = 0.5 * (edges[i] + edges[i + 1]) - shift;
    s.ex[i] = 0.5 * (edges[i + 1] - edges[i]);
  }
  if (!ok) return s;
  // the factor the template carries beyond the flat A0 (the electron fake
  // factor): from the LAB template, which holds the whole inclusive total
  double sumLab = 0;
  for (int i = 0; i < kNPg; ++i) sumLab += PgInteg(in.ci, Form("%s_lab_y%d/qcd_abcd", chgW, i));
  s.kFF = s.incl.A() > 0 ? sumLab / s.incl.A() : 1.0;

  // per group sums + the self-checks
  double maxDevTot = 0, maxDevCR = 0, dB = 0, dC = 0, dD = 0;
  for (int g = 0; g < kNG; ++g)
  {
    const std::vector<int> gb = pOAbcdY::BinsOf(g);
    for (int i : gb)
    {
      s.gOld[g] += s.oldY[i];
      s.gNew[g] += s.newY[i];
      s.eNewY[i] = s.newY[i] * s.g[g].rA();
    }
    s.eGNew[g] = s.gNew[g] * s.g[g].rA();
    s.gx[g] = 0.5 * (edges[gb.front()] + edges[gb.back() + 1]) - shift;
    s.gex[g] = 0.5 * (edges[gb.back() + 1] - edges[gb.front()]);
    const double want = s.g[g].A() * s.kFF;
    maxDevTot = std::max(maxDevTot, want > 0 ? std::fabs(s.gNew[g] / want - 1.0) : 1.0);
    const std::string pre = Form("%s_%s_%s_", chgW, bn, pOAbcdY::kGroupName[g]);
    double wsum = 0;
    for (int i : gb) wsum += PgInteg(in.ci, pre + Form("CRB/w_%s_y%d", bn, i));
    const double crv[6] = {PgInteg(in.ci, pre + "CRB/qcd"), PgInteg(in.ci, pre + "CRC/qcd"), PgInteg(in.ci, pre + "CRD/qcd"),
                           PgInteg(in.ci, pre + "CRB/data_obs"), PgInteg(in.ci, pre + "CRC/data_obs"), wsum};
    const double cnt[6] = {s.g[g].B, s.g[g].C, s.g[g].D, s.g[g].dB, s.g[g].dC, s.g[g].wB};
    for (int k = 0; k < 6; ++k)
      maxDevCR = std::max(maxDevCR, cnt[k] != 0 ? std::fabs(crv[k] / cnt[k] - 1.0) : std::fabs(crv[k]));
    dB += s.g[g].dB; dC += s.g[g].dC; dD += s.g[g].dD;
  }
  printf("  [check] Sum_{y in g} qcd_abcdy = A_g x (template/A0 = %.4f): max |dev| %.1e %s\n", s.kFF, maxDevTot,
         maxDevTot < 1e-6 ? "OK" : "<-- MISMATCH (stale input?)");
  printf("  [check] group CR dirs (CRB/CRC/CRD qcd, CRB/CRC data, Sum CRB w) vs the exported counts: max |dev| %.1e %s\n",
         maxDevCR, maxDevCR < 1e-6 ? "OK" : "<-- MISMATCH (stale input?)");
  if (!fb)
    printf("  [check] lab groups' data vs the inclusive CRs: B %.0f / %.0f  C40 %.0f / %.0f  D %.0f / %.0f  %s\n", dB,
           s.incl.dB, dC, s.incl.dC, dD, s.incl.dD,
           (dB == s.incl.dB && dC == s.incl.dC && dD == s.incl.dD) ? "EXACT" : "MISMATCH");
  else
    printf("  [check] fb groups' data share of the inclusive CRs (the FB window): B %.3f  C40 %.3f  D %.3f\n",
           s.incl.dB > 0 ? dB / s.incl.dB : 0, s.incl.dC > 0 ? dC / s.incl.dC : 0, s.incl.dD > 0 ? dD / s.incl.dD : 0);
  printf("  inclusive: A0 = B0*C40/D0 = %.1f +- %.1f, T = %.4f +- %.4f, R = %.4f +- %.4f%s\n", s.incl.A(),
         s.incl.A() * s.incl.rA(), s.incl.T(), s.incl.T() * s.incl.rT(), s.incl.R(), s.incl.R() * s.incl.rR(),
         in.isElec ? Form("; template/A0 = %.4f (the fake-factor weighting)", s.kFF) : "");

  printf("  group bins     | OLD qcd_abcd | NEW qcd_abcdy +- ABCD stat | NEW/OLD        | T_g = B/D        T_g/T  | "
         "R_g = C40/D      R_g/R\n");
  double so = 0, sn = 0, sne2 = 0;
  for (int g = 0; g < kNG; ++g)
  {
    const PgCounts &c = s.g[g];
    const double r = s.gOld[g] > 0 ? s.gNew[g] / s.gOld[g] : 0, er = r * c.rA();
    printf("  %-4s %-9s | %12.1f | %10.1f +- %6.1f        | %5.3f +- %5.3f | %.4f+-%.4f  %5.3f  | %.4f+-%.4f  %5.3f\n",
           pOAbcdY::kGroupName[g], PgBins(g).c_str(), s.gOld[g], s.gNew[g], s.eGNew[g], r, er, c.T(), c.T() * c.rT(),
           s.incl.T() > 0 ? c.T() / s.incl.T() : 0, c.R(), c.R() * c.rR(), s.incl.R() > 0 ? c.R() / s.incl.R() : 0);
    so += s.gOld[g]; sn += s.gNew[g]; sne2 += s.eGNew[g] * s.eGNew[g];
    if (gPgCsv)
      fprintf(gPgCsv, "group,%s,%s,%s,%s,%.4f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,,,,,,,,\n",
              in.lep.c_str(), chgW, bn, pOAbcdY::kGroupName[g], s.gx[g], s.gOld[g], s.gNew[g], s.eGNew[g], r, er,
              c.T(), c.T() * c.rT(), s.incl.T() > 0 ? c.T() / s.incl.T() : 0, c.R(), c.R() * c.rR(),
              s.incl.R() > 0 ? c.R() / s.incl.R() : 0);
  }
  printf("  sum            | %12.1f | %10.1f +- %6.1f        | %5.3f  (Sum_g A_g / A0 = %.3f)\n", so, sn, std::sqrt(sne2),
         so > 0 ? sn / so : 0, s.incl.A() > 0 ? sn / s.kFF / s.incl.A() : 0);
  printf("  per bin:");
  for (int i = 0; i < kNPg; ++i) printf(" y%d %.1f->%.1f", i, s.oldY[i], s.newY[i]);
  printf("\n");
  for (int i = 0; i < kNPg; ++i)
    if (gPgCsv)
    {
      const double r = s.oldY[i] > 0 ? s.newY[i] / s.oldY[i] : 0;
      fprintf(gPgCsv, "bin,%s,%s,%s,y%d,%.4f,%.4f,%.4f,%.4f,%.5f,%.5f,,,,,,,,,,,,,,\n", in.lep.c_str(), chgW, bn, i, s.x[i],
              s.oldY[i], s.newY[i], s.eNewY[i], r, s.oldY[i] > 0 ? s.eNewY[i] / s.oldY[i] : 0);
    }
  s.ok = true;
  return s;
}

// ---- figures ---------------------------------------------------------------
PlotStyle PgStyle()
{
  PlotStyle ps;
  ps.headerX = 0.47;
  ps.headerY = 0.94;
  ps.subSize = 0.032;
  ps.boxTextSize = 0.026;
  ps.boxX1 = 0.47; ps.boxX2 = 0.93; ps.boxY1 = 0.71; ps.boxY2 = 0.77;
  return ps;
}

TGraphErrors *PgGraph(int n, const double *x, const double *ex, const double *y, const double *ey)
{
  TGraphErrors *g = new TGraphErrors();
  for (int i = 0; i < n; ++i)
  {
    g->SetPoint(i, x[i], y[i]);
    g->SetPointError(i, ex[i], ey ? ey[i] : 0.0);
  }
  return g;
}

const char *PgXTitle(bool fb) { return fb ? "#eta_{CM} (proton-going > 0)" : "#eta_{lab} (proton-going > 0)"; }

GraphTuner PgRatioTuner(double lo, double hi)
{
  return [lo, hi](TCanvas *c, TGraphErrors *g) {
    g->SetMinimum(lo);
    g->SetMaximum(hi);
    c->cd();
    TLine *one = new TLine(g->GetXaxis()->GetXmin(), 1.0, g->GetXaxis()->GetXmax(), 1.0);
    one->SetLineStyle(2); one->SetLineColor(kGray + 2); one->Draw();
  };
}

// per-bin yields, old vs new
void PgPlotYields(const PgInputs &in, const PgSet &s, int ic, bool fb)
{
  double ymax = 0;
  for (int i = 0; i < kNPg; ++i) ymax = std::max({ymax, s.oldY[i], s.newY[i] + s.eNewY[i]});
  std::vector<OverlaySeries> ser(2);
  ser[0].gTot = PgGraph(kNPg, s.x, s.ex, s.oldY, nullptr);
  ser[0].label = "before: qcd_abcd (one ABCD per charge)";
  ser[0].color = kBlack; ser[0].marker = 20; ser[0].markerSize = 1.3; ser[0].xShift = -0.35;
  ser[1].gTot = PgGraph(kNPg, s.x, s.ex, s.newY, s.eNewY);
  ser[1].label = "after: qcd_abcdy (one ABCD per group)";
  ser[1].color = kRed + 1; ser[1].marker = 21; ser[1].markerSize = 1.5; ser[1].xShift = 0.35;
  PlotStyle ps = PgStyle();
  ps.legX1 = 0.17; ps.legY1 = 0.55; ps.legX2 = 0.62; ps.legY2 = 0.67;
  GraphTuner tuner = [ymax](TCanvas *, TGraphErrors *g) { g->SetMinimum(0); g->SetMaximum(2.0 * ymax); };
  SaveNiceGraph_Overlay(ser, Form("%s/pergroup_yields_%s_%s_%s", kPgDir, in.lep.c_str(), ic == 0 ? "Wp" : "Wm",
                                  fb ? "fb" : "lab"),
                        PgXTitle(fb), "SR multijet events (p_{T} > 25, m_{T} > 40)", "",
                        Form("%s channel, %s, %s bins, prefit", in.isElec ? "e" : "#mu", ic == 0 ? "W^{+}" : "W^{-}",
                             fb ? "FB" : "lab"),
                        "groups: y0-1 | y2-9 | y10-11",
                        {"after: group ABCD stat. error (correlated)"}, ps, tuner);
}

// per-group NEW/OLD, W+ and W-
void PgPlotGroupRatio(const PgInputs &in, const PgSet s[2], bool fb)
{
  std::vector<OverlaySeries> ser(2);
  for (int c = 0; c < 2; ++c)
  {
    double r[3], e[3];
    for (int g = 0; g < kNG; ++g)
    {
      r[g] = s[c].gOld[g] > 0 ? s[c].gNew[g] / s[c].gOld[g] : 0;
      e[g] = r[g] * s[c].g[g].rA();
    }
    ser[c].gTot = PgGraph(kNG, s[c].gx, s[c].gex, r, e);
    ser[c].label = c == 0 ? "W^{+}" : "W^{-}";
    ser[c].color = c == 0 ? kAzure + 2 : kRed + 1;
    ser[c].marker = c == 0 ? 21 : 33;
    ser[c].markerSize = c == 0 ? 1.5 : 2.0;
    ser[c].xShift = c == 0 ? -0.15 : 0.15;
  }
  PlotStyle ps = PgStyle();
  ps.legX1 = 0.55; ps.legY1 = 0.15; ps.legX2 = 0.93; ps.legY2 = 0.27;
  SaveNiceGraph_Overlay(ser, Form("%s/pergroup_ratio_%s_%s", kPgDir, in.lep.c_str(), fb ? "fb" : "lab"), PgXTitle(fb),
                        "group multijet: after / before", "",
                        Form("%s channel, %s bins, prefit", in.isElec ? "e" : "#mu", fb ? "FB" : "lab"),
                        "qcd_abcdy / qcd_abcd per group",
                        {"errors: the group ABCD stat. only"}, ps, PgRatioTuner(0.0, 2.2));
}

// T_g/T and R_g/R per group
void PgPlotTR(const PgInputs &in, const PgSet &s, int ic, bool fb)
{
  double t[3], et[3], r[3], er[3];
  for (int g = 0; g < kNG; ++g)
  {
    t[g] = s.g[g].T() / s.incl.T();  et[g] = t[g] * s.g[g].rT();
    r[g] = s.g[g].R() / s.incl.R();  er[g] = r[g] * s.g[g].rR();
  }
  std::vector<OverlaySeries> ser(2);
  ser[0].gTot = PgGraph(kNG, s.gx, s.gex, t, et);
  ser[0].label = Form("T_{g}/T,  T = B/D = %.3f", s.incl.T());
  ser[0].color = kViolet + 1; ser[0].marker = 20; ser[0].markerSize = 1.3; ser[0].xShift = -0.15;
  ser[1].gTot = PgGraph(kNG, s.gx, s.gex, r, er);
  ser[1].label = Form("R_{g}/R,  R = C40/D = %.3f", s.incl.R());
  ser[1].color = kTeal + 3; ser[1].marker = 21; ser[1].markerSize = 1.5; ser[1].xShift = 0.15;
  PlotStyle ps = PgStyle();
  ps.legX1 = 0.25; ps.legY1 = 0.59; ps.legX2 = 0.68; ps.legY2 = 0.69;
  SaveNiceGraph_Overlay(ser, Form("%s/pergroup_TR_%s_%s_%s", kPgDir, in.lep.c_str(), ic == 0 ? "Wp" : "Wm",
                                  fb ? "fb" : "lab"),
                        PgXTitle(fb), "group transfer factor / inclusive", "",
                        Form("%s channel, %s, %s bins, r = 1", in.isElec ? "e" : "#mu", ic == 0 ? "W^{+}" : "W^{-}",
                             fb ? "FB" : "lab"),
                        "T: iso (m_{T}<30), R: m_{T} (anti-iso)", {"errors: data + MC stat."}, ps,
                        PgRatioTuner(0.0, 3.0));
}

// ---- e/mu of the prefit r_Q (open_issues_diag.C convention) ---------------
struct PgRq { double D = 0, W = 0, B = 0, Qo = 0, Qn = 0; };

PgRq PgReadRq(const PgInputs &in, int ic, bool fb, int i)
{
  PgRq r;
  const std::string dir = Form("%s_%s_y%d/", ic == 0 ? "Wp" : "Wm", fb ? "fb" : "lab", i);
  r.D = PgInteg(in.ci, dir + "data_obs");
  r.W = PgInteg(in.ci, dir + "signal") + PgInteg(in.ci, dir + "wtau");
  r.B = PgInteg(in.ci, dir + "z") + PgInteg(in.ci, dir + "ztau") + PgInteg(in.ci, dir + "ttbar");
  r.Qo = PgInteg(in.ci, dir + "qcd_abcd");
  r.Qn = PgInteg(in.ci, dir + "qcd_abcdy");
  return r;
}

void PgAdd(PgRq &s, const PgRq &r) { s.D += r.D; s.W += r.W; s.B += r.B; s.Qo += r.Qo; s.Qn += r.Qn; }

// e/mu of r for Q = old (q = 0) or new (q = 1); data Poisson errors only
void PgEMu(const PgRq &m, const PgRq &e, int q, double &rat, double &err, double &rm, double &re)
{
  const double Qm = q ? m.Qn : m.Qo, Qe = q ? e.Qn : e.Qo;
  rm = (m.D - m.B - Qm) / m.W;
  re = (e.D - e.B - Qe) / e.W;
  const double erm = std::sqrt(m.D) / m.W, ere = std::sqrt(e.D) / e.W;
  rat = re / rm;
  err = std::fabs(rat) * std::sqrt(std::pow(erm / rm, 2) + std::pow(ere / re, 2));
}

void PgEMuAll(const PgInputs &inM, const PgInputs &inE, bool fb)
{
  const char *bn = fb ? "fb" : "lab";
  const double *edges = fb ? pOSkim::kYEdgesFB : pOSkim::kYEdges;
  const double shift = fb ? kPgShift : 0.0;
  printf("\n=== e/mu of the prefit r_Q = (D-B-Q)/W, %s bins (open_issues_diag.C convention, data Poisson errors) ===\n", bn);
  printf("    Q = OLD qcd_abcd (one ABCD per charge) | NEW qcd_abcdy (one ABCD per rapidity group)\n");
  double x[12], ex[12], emu[2][2][12], eemu[2][2][12];
  for (int i = 0; i < kNPg; ++i) { x[i] = 0.5 * (edges[i] + edges[i + 1]) - shift; ex[i] = 0.5 * (edges[i + 1] - edges[i]); }
  for (int c = 0; c < 2; ++c)
  {
    const char *cg = c == 0 ? "Wp" : "Wm";
    printf("\n  %s  bin |  e: D    Qold   Qnew  |  r_e(old) r_e(new)  r_mu(old) r_mu(new) | e/mu (old)      e/mu (new)\n", cg);
    PgRq em, ee, bm, be, im, ie, gm[3], ge[3];
    auto row = [&](const char *lab, const PgRq &m, const PgRq &e, double xx, double *ro, double *eo, double *rn, double *en) {
      double rat[2], err[2], rm[2], re[2];
      for (int q = 0; q < 2; ++q) PgEMu(m, e, q, rat[q], err[q], rm[q], re[q]);
      printf("    %7s | %5.0f %6.1f %6.1f | %7.3f  %7.3f   %7.3f   %7.3f  | %5.3f +- %5.3f  %5.3f +- %5.3f\n", lab, e.D, e.Qo,
             e.Qn, re[0], re[1], rm[0], rm[1], rat[0], err[0], rat[1], err[1]);
      if (ro) { *ro = rat[0]; *eo = err[0]; *rn = rat[1]; *en = err[1]; }
      if (gPgCsv)
        fprintf(gPgCsv, "emu,e/mu,%s,%s,%s,%.4f,,,,,,,,,,,,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", cg, bn, lab, xx,
                re[0], re[1], rm[0], rm[1], rat[0], err[0], rat[1], err[1]);
    };
    for (int i = 0; i < kNPg; ++i)
    {
      const PgRq m = PgReadRq(inM, c, fb, i), e = PgReadRq(inE, c, fb, i);
      PgAdd(im, m); PgAdd(ie, e);
      if (i == 0 || i == kNPg - 1) { PgAdd(em, m); PgAdd(ee, e); }
      if (i <= 3) { PgAdd(bm, m); PgAdd(be, e); }
      PgAdd(gm[pOAbcdY::GroupOf(i)], m); PgAdd(ge[pOAbcdY::GroupOf(i)], e);
      row(Form("y%d", i), m, e, x[i], &emu[c][0][i], &eemu[c][0][i], &emu[c][1][i], &eemu[c][1][i]);
    }
    row("y0+y11", em, ee, 0.0, nullptr, nullptr, nullptr, nullptr);
    row("y0..y3", bm, be, 0.0, nullptr, nullptr, nullptr, nullptr);
    for (int g = 0; g < kNG; ++g) row(pOAbcdY::kGroupName[g], gm[g], ge[g], 0.0, nullptr, nullptr, nullptr, nullptr);
    row("incl", im, ie, 0.0, nullptr, nullptr, nullptr, nullptr);
  }
  for (int c = 0; c < 2; ++c)
  {
    std::vector<OverlaySeries> ser(2);
    for (int q = 0; q < 2; ++q)
    {
      ser[q].gTot = PgGraph(kNPg, x, ex, emu[c][q], eemu[c][q]);
      ser[q].label = q == 0 ? "multijet before: qcd_abcd" : "multijet after: qcd_abcdy";
      ser[q].color = q == 0 ? kGray + 2 : (c == 0 ? kAzure + 2 : kRed + 1);
      ser[q].marker = q == 0 ? 24 : (c == 0 ? 21 : 33);
      ser[q].markerSize = q == 0 ? 1.3 : (c == 0 ? 1.5 : 2.0);
      ser[q].xShift = q == 0 ? -0.35 : 0.35;
    }
    PlotStyle ps = PgStyle();
    ps.legX1 = 0.17; ps.legY1 = 0.59; ps.legX2 = 0.70; ps.legY2 = 0.69;
    SaveNiceGraph_Overlay(ser, Form("%s/pergroup_emu_%s_%s", kPgDir, c == 0 ? "Wp" : "Wm", bn), PgXTitle(fb),
                          "e / #mu  of  (D - B - Q) / W", "",
                          Form("%s, m_{T} > 40, prefit, %s bins", c == 0 ? "W^{+}" : "W^{-}", fb ? "FB" : "lab"),
                          "ABCD per group vs per charge", {"errors: data statistical only"}, ps,
                          PgRatioTuner(0.2, 2.6));
  }
}

} // namespace

void qcd_ysplit_pergroup()
{
  gROOT->SetBatch(kTRUE);
  gSystem->mkdir(kPgDir, kTRUE);
  gPgCsv = std::fopen(Form("%s/qcd_ysplit_pergroup.csv", kPgDir), "w");
  if (gPgCsv)
    std::fprintf(gPgCsv, "kind,flavour,charge,binning,item,x,old,new,new_err,new_over_old,new_over_old_err,T_g,T_g_err,"
                         "T_g_over_T,R_g,R_g_err,R_g_over_R,r_e_old,r_e_new,r_mu_old,r_mu_new,emu_old,emu_old_err,"
                         "emu_new,emu_new_err\n");
  printf("=== qcd_ysplit_pergroup: SR multijet BEFORE (qcd_abcd) vs AFTER (qcd_abcdy, QCD_MODE=abcdy) ===\n");
  printf("    groups (abcd_ygroups.h, same bin indices lab and fb):");
  for (int g = 0; g < kNG; ++g) printf(" %s = %s", pOAbcdY::kGroupName[g], PgBins(g).c_str());
  printf("\n    regions (qcd_abcd.C): B iso-pass m_T<30, C40 anti-iso m_T>40, D anti-iso m_T<30; PREFIT, r = 1\n");
  PgInputs in[2];
  bool ok[2] = {false, false};
  for (int f = 0; f < 2; ++f)
  {
    if (!PgOpen(in[f], f == 1)) continue;
    printf("\n==================== %s ====================\n", in[f].lep.c_str());
    ok[f] = true;
    for (int b = 0; b < 2; ++b)
    {
      PgSet s[2];
      for (int c = 0; c < 2; ++c)
      {
        s[c] = PgAnalyze(in[f], c, b == 1);
        if (!s[c].ok) { ok[f] = false; continue; }
        PgPlotYields(in[f], s[c], c, b == 1);
        PgPlotTR(in[f], s[c], c, b == 1);
      }
      if (s[0].ok && s[1].ok) PgPlotGroupRatio(in[f], s, b == 1);
    }
  }
  if (ok[0] && ok[1])
    for (int b = 0; b < 2; ++b) PgEMuAll(in[0], in[1], b == 1);
  else
    printf("[ERROR] the e/mu comparison needs both flavours\n");
  if (gPgCsv) std::fclose(gPgCsv);
  printf("\n[qcd_ysplit_pergroup] figures + qcd_ysplit_pergroup.csv in %s/\n", kPgDir);
}
