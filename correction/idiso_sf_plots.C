// correction/idiso_sf_plots.C
//
// The SF plots of the ELECTRON ID+ISO SCALE-FACTOR CROSS-CHECK (2026-09-26).
// Reads the Combine results downloaded from lxplus (`sync_lxplus.sh
// download-idiso`) -- the fork's test/pO_idiso_sf_out/<tag>/summary/:
// idiso_sf_<tag>.csv (the fitted eps / SF per bin) and the extraction log,
// whose [expected] lines carry the Asimov precision and whose "QCD in-fit ABCD"
// lines carry the SR multipliers -- plus the EGM wp90iso table and the skim's
// EGM-cell maps of passing MC electrons (the EGM prediction for our bins).
// Definitions: idiso_sf_common.h.
//
// Writes to plots/idiso_sf_ele/:
//   sf_summary             every measurement on one axis: the Z tag-and-probe and
//                          the 9 W variants (inclusive / |eta| / pT combined x 3
//                          sideband windows) with the fitted error bars, the
//                          Asimov-expected precision as boxes, and the EGM
//                          predictions for the same electrons as bands
//   sf_egm_abseta          the EGM overlay the study was built for. Top: our
//   sf_egm_pt              coarse-bin W SF (3 windows) with the expected precision,
//                          the EGM prediction per coarse bin, the Z tag-and-probe
//                          band. Bottom (zoom): the EGM fine-binned SF per pT /
//                          |eta| slice (+-eta averaged) with the same two references
//   abcd_multiplier_<scheme>  the in-fit ABCD SR multiplier kappa^theta sB sC / sD
//                          per channel and window: how far the fit moved each QCD
//   fail_qcdshape_<window> the inclusive fail SR (charges summed): data - EWK at
//                          r = 1 vs the prefit in-fit-ABCD QCD template -- the
//                          shape the fail W count has to be read off
// + plots/idiso_sf_ele/idiso_sf_summary.csv (every number drawn) and the record
// on stdout: [FIT] [SF] [EGM] [PULL] [OUT].
// Markers: sideband window 0.3-1.0 / 0.3-0.6 / 0.6-1.0 = circle / square /
// diamond (the repo's centred sizes); FILLED = measured, OPEN = a fitted fail r of
// the bin (either charge) sits on a LIMIT of its range [0, 10] -- at 0 (no fail W,
// eps_data = 1, SF = 1/eps_MC ~ 1.24) or at 10 (the upper edge of the POI range,
// SF ~ 0.41): the fitted error is then meaningless and no error bar is drawn --
// the grey box (the expected precision) is the honest one. (Until 2026-09-27 only
// the 0 side was flagged, via eps_data > 0.99; the 0.6-1.0 fits sit at 10.)
//
// Run from correction/ after the download:
//   ./run_idiso_sf.sh plots                  (keeps logs/idiso_sf_ele_plots.log)
// env FORK_TEST = the fork's test dir (default /Users/zhenghuang/HiggsAnalysis-CombinedLimit/test)

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TBox.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TSystem.h"
#include "TStyle.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "../skim/skim_common.h"
#include "../skim/mc_norm.h"
#include "../plotting/plotting_helper.C"
#include "idiso_sf_common.h"

using namespace pOIdIso;

namespace
{

const char *const kWinSfx[3]   = {"", "_sblo", "_sbhi"};
const char *const kWinText[3]  = {"0.3-1.0", "0.3-0.6", "0.6-1.0"};
const char *const kWinShort[3] = {"full", "sblo", "sbhi"};
const int    kWinMarker[3]     = {20, 21, 33};   // filled circle / square / diamond (measured)
const int    kWinMarkerOpen[3] = {24, 25, 27};   // their open twins (the fail W at 0)
const double kWinSize[3]       = {1.3, 1.5, 2.0}; // the repo's centred sizes (xsec_fiducial.C kMk*)
const int    kWinColor[3]      = {kBlack, kBlue + 1, kRed + 1};
const int    kSliceColor[4]    = {kAzure + 2, kOrange + 7, kMagenta + 1, kTeal + 3};
const int    kZColor           = kGreen + 2;
const int    kEgmPredColor     = kOrange + 1;

std::string FitDir()
{
  const char *e = gSystem->Getenv("FORK_TEST");
  return std::string(e ? e : "/Users/zhenghuang/HiggsAnalysis-CombinedLimit/test") + "/pO_idiso_sf_out";
}

// ------------------------------------------------------------------ the fit results
// the range of every r POI (the maps of make_pO_idiso_sf_cards.sh: [1,0,10])
const double kRMin = 0.0, kRMax = 10.0;

struct Row
{
  std::string bin, charge;
  double lo = 0, hi = 0, epsMC = 0, epsData = 0, epsDataErr = 0, sf = 0, sfErr = 0, nPass = 0, nFail = 0;
  double rPass = -1, rFail = -1; // single (pass, fail) pair rows only; -1 on the combined rows
  int status = -1, covQual = -1;
  int bnd = 0; // a fail r of the bin at a limit: 1 = at 0, 2 = at 10, 3 = both (OR over the bin's rows)
  bool Boundary() const { return bnd != 0; }
};
const char *BndText(int bnd)
{
  return bnd == 1 ? "fail r at 0" : bnd == 2 ? "fail r at its limit 10" : "fail r at 0 and at 10";
}
// the same on the summary plot, where the numbers column is narrow
const char *BndShort(int bnd) { return bnd == 1 ? "r_{fail} = 0" : bnd == 2 ? "r_{fail} = 10" : "r_{fail} = 0, 10"; }
const char *BndCode(int bnd) { return bnd == 0 ? "no" : bnd == 1 ? "r0" : bnd == 2 ? "r10" : "r0+r10"; }
struct Fit
{
  std::string tag;
  bool ok = false;
  std::vector<Row> rows;
  std::map<std::string, double> expErr;                  // bin -> expected SF error (the Asimov fit)
  std::map<std::string, std::pair<double, double>> mult; // SR channel -> in-fit ABCD multiplier +- error
  std::vector<std::string> multOrder;
  const Row *Get(const std::string &bin, const std::string &chg = "all") const
  {
    for (const auto &r : rows)
      if (r.bin == bin && r.charge == chg) return &r;
    return nullptr;
  }
  double Exp(const std::string &bin) const
  {
    auto it = expErr.find(bin);
    return it == expErr.end() ? -1.0 : it->second;
  }
};

std::vector<std::string> SplitCsv(const std::string &s)
{
  std::vector<std::string> v;
  std::stringstream ss(s);
  std::string t;
  while (std::getline(ss, t, ',')) v.push_back(t);
  return v;
}
std::vector<std::string> Tokens(const std::string &s)
{
  std::istringstream ss(s);
  std::vector<std::string> v;
  std::string t;
  while (ss >> t) v.push_back(t);
  return v;
}

Fit ReadFit(const std::string &tag)
{
  Fit f;
  f.tag = tag;
  const std::string dir = FitDir() + "/" + tag + "/summary/";
  std::ifstream csv(dir + "idiso_sf_" + tag + ".csv");
  if (!csv)
  {
    std::cout << "[WARN] no " << dir << "idiso_sf_" << tag << ".csv -- not downloaded?\n";
    return f;
  }
  std::string line;
  std::getline(csv, line); // header
  while (std::getline(csv, line))
  {
    const auto c = SplitCsv(line);
    if (c.size() < 24) continue;
    Row r;
    r.bin = c[4]; r.lo = std::stod(c[5]); r.hi = std::stod(c[6]); r.charge = c[7];
    r.epsMC = std::stod(c[10]); r.rPass = std::stod(c[11]); r.rFail = std::stod(c[13]);
    r.nPass = std::stod(c[16]); r.nFail = std::stod(c[17]);
    r.epsData = std::stod(c[18]); r.epsDataErr = std::stod(c[19]); r.sf = std::stod(c[20]); r.sfErr = std::stod(c[21]);
    r.status = std::stoi(c[22]); r.covQual = std::stoi(c[23]);
    f.rows.push_back(r);
  }
  // a fail r on a limit of its range makes the fitted error of every row built on it
  // meaningless: flag the bin (both charges and the charge-combined row) and, for the
  // scheme-combined row "all", every bin. Pass r's are checked too, and only warned.
  std::map<std::string, int> binFlag;
  int allFlag = 0;
  for (const auto &r : f.rows)
  {
    if (r.rFail < 0) continue; // the combined rows carry no r
    int fl = 0;
    if (r.rFail < kRMin + 1e-3) fl |= 1;
    if (r.rFail > kRMax - 1e-2) fl |= 2;
    binFlag[r.bin] |= fl;
    allFlag |= fl;
    if (r.rPass < kRMin + 1e-3 || r.rPass > kRMax - 1e-2)
      std::cout << Form("[WARN] %s bin %s %s: r_pass = %.3g sits on a limit of [%g, %g]\n", tag.c_str(), r.bin.c_str(),
                        r.charge.c_str(), r.rPass, kRMin, kRMax);
  }
  for (auto &r : f.rows) r.bnd = (r.bin == "all") ? allFlag : binFlag[r.bin];
  std::ifstream log(dir + "extract_idiso_sf_" + tag + ".log");
  while (std::getline(log, line))
  {
    const auto t = Tokens(line);
    if (t.size() >= 6 && t[0] == "[expected]" && t[2] == "bin") f.expErr[t[3]] = std::stod(t.back());
    if (t.size() >= 8 && t[0] == "[idiso-sf]" && t[2] == "QCD" && t[4] == "ABCD")
    {
      f.mult[t[5]] = {std::stod(t[t.size() - 3]), std::stod(t.back())};
      f.multOrder.push_back(t[5]);
    }
  }
  f.ok = !f.rows.empty();
  if (f.ok)
    std::cout << Form("[FIT] %-24s status %d, covQual %d, %zu rows, %zu expected errors, %zu ABCD multipliers\n",
                      tag.c_str(), f.rows[0].status, f.rows[0].covQual, f.rows.size(), f.expErr.size(), f.mult.size());
  return f;
}

// ------------------------------------------------------------------ the EGM predictions
TH2D *EgmMap(const std::vector<int> &samples, const std::vector<std::string> &names)
{
  TH2D *out = nullptr;
  for (int is : samples)
  {
    TFile *f = TFile::Open(SkimFile(kSampleName[is]).c_str(), "READ");
    if (!f || f->IsZombie()) { std::cerr << "[FATAL] cannot open " << SkimFile(kSampleName[is]) << "\n"; return nullptr; }
    const double k = pONorm::MCScale(kSampleNormLabel[is]);
    for (const auto &n : names)
    {
      TH2D *h = (TH2D *)f->Get(n.c_str());
      if (!h) { std::cerr << "[FATAL] " << kSampleName[is] << ": missing " << n << "\n"; return nullptr; }
      TH2D *c = (TH2D *)h->Clone(Form("%s_%s_egm", n.c_str(), kSampleName[is]));
      c->SetDirectory(nullptr);
      c->Scale(k);
      if (!out) out = c;
      else { out->Add(c); delete c; }
    }
    f->Close();
  }
  return out;
}

// the EGM SF of a folded (|eta|, pT) cell: the +eta and -eta cells averaged
bool EgmFold(double aLo, double aHi, double pLo, double pHi, double &sf, double &err)
{
  const double eta = 0.5 * (aLo + aHi), pt = 0.5 * (pLo + pHi);
  const EgmCell *a = FindEgm(eta, pt), *b = FindEgm(-eta, pt);
  if (!a || !b) return false;
  sf = 0.5 * (a->sf + b->sf);
  err = 0.5 * (a->err + b->err);
  return true;
}

TBox *Box(double x1, double y1, double x2, double y2, int col, double alpha, int style = 1001)
{
  TBox *b = new TBox(x1, y1, x2, y2);
  b->SetFillColorAlpha(col, alpha);
  b->SetFillStyle(style);
  b->SetLineColor(col);
  b->SetLineWidth(0);
  return b;
}
TLine *Line(double x1, double y1, double x2, double y2, int col, int style, int width)
{
  TLine *l = new TLine(x1, y1, x2, y2);
  l->SetLineColor(col);
  l->SetLineStyle(style);
  l->SetLineWidth(width);
  return l;
}
TGraphAsymmErrors *Point(double x, double y, double exl, double exh, double ey, int marker, int col, double size)
{
  TGraphAsymmErrors *g = new TGraphAsymmErrors(1);
  g->SetPoint(0, x, y);
  g->SetPointError(0, exl, exh, ey, ey);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(size);
  g->SetMarkerColor(col);
  g->SetLineColor(col);
  g->SetLineWidth(2);
  return g;
}
TLegend *Legend(double x1, double y1, double x2, double y2, int ncol, double size, bool opaque = false)
{
  TLegend *l = new TLegend(x1, y1, x2, y2);
  l->SetNColumns(ncol);
  l->SetBorderSize(0);
  l->SetFillStyle(opaque ? 1001 : 0);
  l->SetFillColor(0);
  l->SetTextFont(42);
  l->SetTextSize(size);
  return l;
}
// legend-only markers: every window gets its entry even when all its points are open
TGraphAsymmErrors *LegMarker(int marker, int col, double size)
{
  TGraphAsymmErrors *g = new TGraphAsymmErrors(1);
  g->SetMarkerStyle(marker); g->SetMarkerSize(size);
  g->SetMarkerColor(col); g->SetLineColor(col); g->SetLineWidth(2);
  return g;
}
void Save(TCanvas *c, const std::string &path)
{
  c->SaveAs((path + ".png").c_str());
  c->SaveAs((path + ".pdf").c_str());
  std::cout << "[OUT] " << path << ".{png,pdf}\n";
}

struct Pred { double sf = 0, err = 0; };

// ------------------------------------------------------------------ 1. the summary
void DrawSummary(const std::vector<Fit> &w, const Fit &z, const Pred &pz, const Pred &pw)
{
  struct Item { std::string label; double sf, fitErr, expErr; int bnd; bool isZ; int win; };
  std::vector<Item> items;
  if (const Row *r = z.Get("ztnp"))
    items.push_back({"Z tag-and-probe, inclusive", r->sf, r->sfErr, z.Exp("ztnp"), r->bnd, true, -1});
  const char *schemeLabel[kNScheme] = {"W inclusive", "W |#eta_{SC}| bins combined", "W p_{T} bins combined"};
  for (int s = 0; s < kNScheme; ++s)
    for (int iw = 0; iw < 3; ++iw)
    {
      const Fit &f = w[3 * s + iw];
      const Row *r = f.Get(s == 0 ? "0" : "all"); // "all" carries the OR over every bin of the scheme
      if (!r) continue;
      items.push_back({Form("%s, sideband relIso %s", schemeLabel[s], kWinText[iw]), r->sf, r->sfErr,
                       f.Exp(s == 0 ? "0" : "all"), r->bnd, false, iw});
    }
  const int n = (int)items.size();
  PlotStyle ps;
  ps.w = 1150; ps.h = 900; ps.lm = 0.35; ps.rm = 0.03; ps.tm = 0.07; ps.bm = 0.25;
  TCanvas *c = new TCanvas("c_sf_summary", "", ps.w, ps.h);
  c->SetLeftMargin(ps.lm); c->SetRightMargin(ps.rm); c->SetTopMargin(ps.tm); c->SetBottomMargin(ps.bm);
  c->SetTicks(1, 0);
  const double xlo = 0.2, xhi = 2.35;
  TH2F *fr = new TH2F("fr_sf_summary", "", 10, xlo, xhi, n, 0, n);
  fr->SetDirectory(nullptr);
  fr->SetStats(0);
  ApplyHistStyle(fr, ps, "data / MC scale factor", "");
  fr->GetXaxis()->SetTitleOffset(1.0);
  fr->GetYaxis()->SetLabelSize(0.027);
  fr->GetYaxis()->SetTickLength(0);
  for (int i = 0; i < n; ++i) fr->GetYaxis()->SetBinLabel(n - i, items[i].label.c_str());
  fr->Draw("AXIS");
  // the EGM predictions for the same electrons: bands over the Z row and the W rows
  Box(pz.sf - pz.err, n - 1, pz.sf + pz.err, n, kEgmPredColor, 0.45)->Draw();
  Line(pz.sf, n - 1, pz.sf, n, kEgmPredColor, 1, 2)->Draw();
  Box(pw.sf - pw.err, 0, pw.sf + pw.err, n - 1, kEgmPredColor, 0.30)->Draw();
  Line(pw.sf, 0, pw.sf, n - 1, kEgmPredColor, 1, 2)->Draw();
  Line(1.0, 0, 1.0, n, kGray + 2, 2, 1)->Draw();
  Line(xlo, n - 1, xhi, n - 1, kGray + 1, 3, 1)->Draw(); // Z | W separator
  TGraphAsymmErrors *gLeg[5] = {};                        // Z, W windows 0-2, boundary
  TBox *bExp = nullptr;
  for (int i = 0; i < n; ++i)
  {
    const Item &it = items[i];
    const double y = n - i - 0.5;
    if (it.expErr > 0)
    {
      TBox *b = Box(it.sf - it.expErr, y - 0.3, it.sf + it.expErr, y + 0.3, it.isZ ? kZColor : kGray + 1, 0.35);
      b->Draw();
      if (!it.isZ) bExp = b;
    }
    const int col = it.isZ ? kZColor : kWinColor[it.win];
    const bool open = it.bnd != 0;
    const int mk = it.isZ ? 21 : (open ? kWinMarkerOpen[it.win] : kWinMarker[it.win]);
    const double sz = it.isZ ? 1.5 : kWinSize[it.win];
    TGraphAsymmErrors *g = Point(it.sf, y, open ? 0.0 : it.fitErr, open ? 0.0 : it.fitErr, 0.0, mk, col, sz);
    g->Draw("P SAME");
    const int slot = it.isZ ? 0 : (open ? 4 : 1 + it.win);
    if (!gLeg[slot]) gLeg[slot] = g;
    TLatex *t = new TLatex(1.50, y, open ? Form("%.3f  [%s; exp. #pm %.3f]", it.sf, BndShort(it.bnd), it.expErr)
                                         : Form("%.3f #pm %.3f  [exp. #pm %.3f]", it.sf, it.fitErr, it.expErr));
    t->SetTextFont(42); t->SetTextSize(0.021); t->SetTextAlign(12); t->SetTextColor(col);
    t->Draw();
  }
  TLegend *leg = Legend(0.02, 0.01, 0.99, 0.13, 3, 0.022);
  if (gLeg[0]) leg->AddEntry(gLeg[0], "Z tag-and-probe fit", "lp");
  for (int iw = 0; iw < 3; ++iw)
    leg->AddEntry(LegMarker(kWinMarker[iw], kWinColor[iw], kWinSize[iw]), Form("W fit, sideband relIso %s", kWinText[iw]), "lp");
  leg->AddEntry(LegMarker(24, kGray + 2, 1.3), "open: a bin's r_{fail} at 0 or 10 (its limits)", "p");
  if (bExp) leg->AddEntry(bExp, "expected precision (Asimov fit)", "f");
  leg->AddEntry(Box(0, 0, 1, 1, kEgmPredColor, 0.40), "EGM wp90iso for the same electrons", "f");
  leg->Draw();
  (void)gLeg;
  TLatex hdr;
  hdr.SetNDC(); hdr.SetTextFont(42); hdr.SetTextSize(0.034);
  hdr.DrawLatex(0.02, 0.945, "#bf{CMS} #it{Work in Progress}");
  hdr.SetTextSize(0.027); hdr.SetTextAlign(31);
  hdr.DrawLatex(1 - ps.rm, 0.945, "Electron WP90 ID + WP90 iso, pO 46.5 nb^{-1}");
  c->RedrawAxis();
  Save(c, Form("%s/sf_summary", kPlotDir));
}

// ------------------------------------------------------------------ 2. the EGM overlay
void DrawEgmOverlay(int s, const std::vector<Fit> &w, const Row *zr, const std::map<std::string, Pred> &pred)
{
  const bool eta = (s == 1);
  const double xlo = eta ? 0.0 : 25.0, xhi = eta ? 2.5 : 120.0;
  const double ylo = 0.2, yhi = 2.9, zlo = 0.88, zhi = 1.04;
  const Scheme &sc = kSchemes[s];
  PlotStyle ps;
  ps.w = 900; ps.h = 1000; ps.lm = 0.13; ps.rm = 0.04;
  TCanvas *c = new TCanvas(Form("c_egm_%d", s), "", ps.w, ps.h);
  const double split = 0.30;
  TPad *pT = new TPad(Form("pT_egm_%d", s), "", 0, split, 1, 1);
  TPad *pB = new TPad(Form("pB_egm_%d", s), "", 0, 0, 1, split);
  pT->SetLeftMargin(ps.lm); pT->SetRightMargin(ps.rm); pT->SetTopMargin(0.07); pT->SetBottomMargin(0.02);
  pB->SetLeftMargin(ps.lm); pB->SetRightMargin(ps.rm); pB->SetTopMargin(0.04); pB->SetBottomMargin(0.33);
  pT->SetTicks(1, 1); pB->SetTicks(1, 1);
  pT->Draw(); pB->Draw();

  // the EGM fine-binned SF (+-eta averaged): slices in pT (|eta| plot) or |eta| (pT plot)
  const double etaLo[4] = {0.0, 0.8, 1.566, 2.0}, etaHi[4] = {0.8, 1.444, 2.0, 2.5};
  const double ptLo[5] = {25.0, 35.0, 50.0, 70.0, 100.0}, ptHi[5] = {35.0, 50.0, 70.0, 100.0, 120.0};
  auto slices = [&](bool draw, TLegend *leg) {
    for (int sl = 0; sl < 4; ++sl)
    {
      const int nCell = eta ? 4 : 5;
      bool any = false;
      for (int ce = 0; ce < nCell; ++ce)
      {
        double sf = 0, err = 0;
        const double aLo = eta ? etaLo[ce] : etaLo[sl], aHi = eta ? etaHi[ce] : etaHi[sl];
        const double pLo = eta ? ptLo[sl] : ptLo[ce], pHi = eta ? ptHi[sl] : ptHi[ce];
        if (!EgmFold(aLo, aHi, pLo, pHi, sf, err)) continue;
        any = true;
        if (!draw) continue;
        const double x1 = eta ? aLo : pLo, x2 = eta ? aHi : pHi;
        Box(x1, sf - err, x2, sf + err, kSliceColor[sl], 0.25)->Draw();
        Line(x1, sf, x2, sf, kSliceColor[sl], 1, 3)->Draw();
      }
      if (any && leg)
      {
        TBox *e = Box(0, 0, 1, 1, kSliceColor[sl], 0.5);
        e->SetLineColor(kSliceColor[sl]); e->SetLineWidth(3);
        leg->AddEntry(e, eta ? Form("EGM wp90iso, p_{T} %.0f-%.0f GeV", ptLo[sl], ptHi[sl])
                             : Form("EGM wp90iso, |#eta_{SC}| %.4g-%.4g", etaLo[sl], etaHi[sl]),
                      "lf");
      }
    }
  };
  // the references drawn in both pads: the EGM prediction per coarse bin, the Z band, 1
  auto references = [&](double y1, double y2) {
    if (eta) Box(pOSkim::kEcalGapLo, y1, pOSkim::kEcalGapHi, y2, kGray + 1, 0.6, 3004)->Draw();
    if (zr)
    {
      Box(xlo, zr->sf - zr->sfErr, xhi, zr->sf + zr->sfErr, kZColor, 0.22)->Draw();
      Line(xlo, zr->sf, xhi, zr->sf, kZColor, 1, 2)->Draw();
    }
    Line(xlo, 1.0, xhi, 1.0, kGray + 2, 2, 1)->Draw();
    for (int k = 0; k < sc.n; ++k)
    {
      auto it = pred.find(Form("%s_%d", sc.name, k));
      if (it == pred.end()) continue;
      const double lo = sc.lo[k], hi = (sc.hi[k] >= kPtInf) ? xhi : sc.hi[k];
      Line(lo, it->second.sf, hi, it->second.sf, kBlack, 2, 2)->Draw();
    }
  };

  // ---- top: our W SF ----
  pT->cd();
  TH1F *fr = new TH1F(Form("fr_egm_%d", s), "", 100, xlo, xhi);
  fr->SetDirectory(nullptr);
  fr->SetStats(0);
  ApplyHistStyle(fr, ps, "", "data / MC scale factor");
  fr->GetXaxis()->SetLabelSize(0);
  fr->GetYaxis()->SetTitleSize(0.050); fr->GetYaxis()->SetLabelSize(0.042); fr->GetYaxis()->SetTitleOffset(1.10);
  fr->SetMinimum(ylo); fr->SetMaximum(yhi);
  fr->Draw("AXIS");
  references(ylo, 1.55);
  TLegend *leg = Legend(0.15, 0.47, 0.95, 0.77, 2, 0.026, true);
  if (zr) leg->AddEntry(Box(0, 0, 1, 1, kZColor, 0.35), Form("Z tag-and-probe: %.3f #pm %.3f", zr->sf, zr->sfErr), "f");
  leg->AddEntry(Line(0, 0, 1, 1, kBlack, 2, 2), "EGM prediction for our W electrons", "l");
  TBox *bExp = nullptr;
  TGraphAsymmErrors *gw[3] = {}, *gOpen = nullptr;
  for (int k = 0; k < sc.n; ++k)
  {
    const double lo = sc.lo[k], hi = (sc.hi[k] >= kPtInf) ? xhi : sc.hi[k];
    const double half = 0.5 * (hi - lo), cen = 0.5 * (hi + lo);
    for (int iw = 0; iw < 3; ++iw)
    {
      const Row *r = w[3 * s + iw].Get(std::to_string(k));
      if (!r) continue;
      const double dx = (iw - 1) * 0.30 * half;
      if (iw == 0)
      {
        const double e = w[3 * s + iw].Exp(std::to_string(k));
        if (e > 0)
        {
          bExp = Box(cen + dx - 0.12 * half, std::max(ylo, r->sf - e), cen + dx + 0.12 * half, std::min(yhi, r->sf + e),
                     kGray + 1, 0.45);
          bExp->Draw();
        }
      }
      const bool b = r->Boundary();
      TGraphAsymmErrors *g = Point(cen + dx, r->sf, 0.0, 0.0, b ? 0.0 : r->sfErr,
                                   b ? kWinMarkerOpen[iw] : kWinMarker[iw], kWinColor[iw], kWinSize[iw]);
      g->Draw("P SAME");
      if (b) { if (!gOpen) gOpen = g; }
      else if (!gw[iw]) gw[iw] = g;
    }
    if (k > 0) Line(lo, ylo, lo, 1.55, kGray, 3, 1)->Draw(); // the coarse-bin boundaries
  }
  for (int iw = 0; iw < 3; ++iw)
    leg->AddEntry(LegMarker(kWinMarker[iw], kWinColor[iw], kWinSize[iw]), Form("W fit, sideband relIso %s", kWinText[iw]), "lp");
  if (gOpen) leg->AddEntry(LegMarker(24, kGray + 2, 1.3), "open: r_{fail} at 0 or 10 (no fit error)", "p");
  if (bExp) leg->AddEntry(bExp, "W expected precision (Asimov, 0.3-1.0)", "f");
  (void)gw;
  slices(false, leg);
  leg->Draw();
  CMS_lumi(pT, 13, 10);
  TLatex t;
  t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.034); t.SetTextAlign(31);
  t.DrawLatex(0.94, 0.86, "Electron WP90 ID + WP90 iso");
  t.SetTextSize(0.026);
  t.DrawLatex(0.94, 0.825, "W: m_{T} > 40, in-fit ABCD;  Z: tag-and-probe");
  pT->RedrawAxis();

  // ---- bottom (zoom): the EGM fine-binned SF + the same references ----
  pB->cd();
  TH1F *fz = new TH1F(Form("fz_egm_%d", s), "", 100, xlo, xhi);
  fz->SetDirectory(nullptr);
  fz->SetStats(0);
  ApplyHistStyle(fz, ps, eta ? "|#eta_{SC}|" : "p_{T}^{e} (GeV)", "SF (zoom)");
  const double sfac = (1.0 - split) / split;
  fz->GetXaxis()->SetTitleSize(0.045 * sfac); fz->GetXaxis()->SetLabelSize(0.040 * sfac);
  fz->GetXaxis()->SetTitleOffset(1.0);
  fz->GetYaxis()->SetTitleSize(0.050 * sfac * 0.85); fz->GetYaxis()->SetLabelSize(0.042 * sfac * 0.85);
  fz->GetYaxis()->SetTitleOffset(1.10 / sfac / 0.85); fz->GetYaxis()->SetNdivisions(504);
  fz->SetMinimum(zlo); fz->SetMaximum(zhi);
  fz->Draw("AXIS");
  references(zlo, zhi);
  slices(true, nullptr);
  pB->RedrawAxis();

  Save(c, Form("%s/sf_egm_%s", kPlotDir, sc.name));
}

// ------------------------------------------------------------------ 3. the ABCD multipliers
void DrawMultipliers(int s, const std::vector<Fit> &w)
{
  const Fit &f0 = w[3 * s];
  const int n = (int)f0.multOrder.size();
  if (n == 0) return;
  PlotStyle ps;
  ps.w = 1000; ps.h = 750; ps.bm = 0.19;
  TCanvas *c = new TCanvas(Form("c_mult_%d", s), "", ps.w, ps.h);
  c->SetLeftMargin(ps.lm); c->SetRightMargin(ps.rm); c->SetTopMargin(ps.tm); c->SetBottomMargin(ps.bm);
  c->SetTicks(0, 1);
  TH1F *fr = new TH1F(Form("fr_mult_%d", s), "", n, 0, n);
  fr->SetDirectory(nullptr);
  fr->SetStats(0);
  ApplyHistStyle(fr, ps, "", "SR QCD multiplier #kappa^{#theta} s_{B} s_{C} / s_{D}");
  for (int i = 0; i < n; ++i)
  {
    TString t(f0.multOrder[i].c_str()); // Wp_pass_k0 -> W^{+} pass k0
    t.ReplaceAll("Wp_", "W^{+} ").ReplaceAll("Wm_", "W^{-} ").ReplaceAll("_k", " k");
    fr->GetXaxis()->SetBinLabel(i + 1, t.Data());
  }
  fr->GetXaxis()->LabelsOption("v");
  fr->GetXaxis()->SetLabelSize(n > 8 ? 0.032 : 0.042);
  fr->SetMinimum(0.0); fr->SetMaximum(1.75);
  fr->Draw("AXIS");
  Line(0, 1.0, n, 1.0, kGray + 2, 2, 1)->Draw();
  TLegend *leg = Legend(0.55, 0.70, 0.93, 0.86, 1, 0.032);
  for (int iw = 0; iw < 3; ++iw)
  {
    const Fit &f = w[3 * s + iw];
    TGraphAsymmErrors *g = new TGraphAsymmErrors();
    int ip = 0;
    for (int i = 0; i < n; ++i)
    {
      auto it = f.mult.find(f0.multOrder[i]);
      if (it == f.mult.end()) continue;
      g->SetPoint(ip, i + 0.5 + (iw - 1) * 0.24, it->second.first);
      g->SetPointError(ip, 0, 0, it->second.second, it->second.second);
      ++ip;
    }
    g->SetMarkerStyle(kWinMarker[iw]); g->SetMarkerSize(kWinSize[iw]);
    g->SetMarkerColor(kWinColor[iw]); g->SetLineColor(kWinColor[iw]); g->SetLineWidth(2);
    g->Draw("P SAME");
    leg->AddEntry(g, Form("sideband relIso %s", kWinText[iw]), "lp");
  }
  leg->Draw();
  CMS_lumi(c, 13, 10);
  TLatex t;
  t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.036); t.SetTextAlign(11);
  t.DrawLatex(0.18, 0.76, Form("In-fit ABCD, %s", s == 0 ? "inclusive" : (s == 1 ? "|#eta_{SC}| bins" : "p_{T} bins")));
  t.SetTextSize(0.029);
  t.DrawLatex(0.18, 0.72, "postfit / prefit QCD in each W signal region");
  c->RedrawAxis();
  Save(c, Form("%s/abcd_multiplier_%s", kPlotDir, kSchemes[s].name));
}

// ------------------------------------------------------------------ 4. the fail-channel QCD shape
void DrawFailShape(int iw)
{
  const std::string fn = Form("%s/combine_input_idiso_incl_%s%s.root", kRootDir, kDiscs[kDiscFit].name, kWinSfx[iw]);
  TFile *f = TFile::Open(fn.c_str(), "READ");
  if (!f || f->IsZombie()) { std::cout << "[WARN] no " << fn << "\n"; return; }
  TH1D *res = nullptr, *qcd = nullptr;
  for (const char *ch : {"Wp_fail_k0", "Wm_fail_k0"})
  {
    TH1D *d = (TH1D *)f->Get(Form("%s/data_obs", ch));
    TH1D *q = (TH1D *)f->Get(Form("%s/qcd", ch));
    if (!d || !q) { std::cout << "[WARN] " << fn << ": missing " << ch << "\n"; return; }
    TH1D *r = (TH1D *)d->Clone(Form("res_%s_%d", ch, iw));
    r->SetDirectory(nullptr);
    for (const char *p : {"signal", "z", "ztau", "wtau"})
    {
      TH1D *h = (TH1D *)f->Get(Form("%s/%s", ch, p));
      if (h) r->Add(h, -1.0);
    }
    TH1D *qc = (TH1D *)q->Clone(Form("qcd_%s_%d", ch, iw));
    qc->SetDirectory(nullptr);
    if (!res) { res = r; qcd = qc; }
    else { res->Add(r); qcd->Add(qc); delete r; delete qc; }
  }
  PlotStyle ps;
  ps.logy = true;
  ps.yTitleOffset = 1.55;
  ps.headerX = 0.55; ps.headerY = 0.66; ps.titleSize = 0.036;
  ps.xRangeLo = 24.0; ps.xRangeHi = 200.0;
  ps.ratioLo = 0.4; ps.ratioHi = 3.4;
  SaveDataMCRatio(res, qcd, Form("%s/fail_qcdshape_%s", kPlotDir, kWinShort[iw]), kDiscs[kDiscFit].axisTitle,
                  "Events / 2.0 GeV", "fail, inclusive, W^{+} + W^{-}", Form("sideband relIso %s", kWinText[iw]),
                  "m_{T} > 40, prefit (r = 1)", ps, false, "Data #minus EWK MC", "QCD (in-fit ABCD, prefit)",
                  "(Data #minus EWK) / QCD");
  std::cout << "[OUT] " << kPlotDir << "/fail_qcdshape_" << kWinShort[iw] << ".{png,pdf}\n";
  f->Close();
}

} // namespace

int idiso_sf_plots()
{
  gErrorIgnoreLevel = kWarning;
  gStyle->SetOptStat(0);
  gSystem->mkdir(kPlotDir, kTRUE);
  std::cout << "[INPUT] fit results: " << FitDir() << "\n";

  // ---------------- the fits ----------------
  std::vector<Fit> w; // index 3 * scheme + window
  for (int s = 0; s < kNScheme; ++s)
    for (int iw = 0; iw < 3; ++iw) w.push_back(ReadFit(Form("%s_%s%s", kSchemes[s].name, kDiscs[kDiscFit].name, kWinSfx[iw])));
  const Fit z = ReadFit("ztnp");
  if (!z.ok) { std::cerr << "[FATAL] no Z tag-and-probe result -- run sync_lxplus.sh download-idiso\n"; return 2; }
  for (const auto &f : w)
    if (!f.ok) { std::cerr << "[FATAL] a W fit is missing (" << f.tag << ") -- run sync_lxplus.sh download-idiso\n"; return 2; }

  // ---------------- the EGM predictions ----------------
  if (!LoadEgm()) { std::cerr << "[FATAL] cannot read " << kEgmCsv << "\n"; return 2; }
  std::map<std::string, Pred> pred;
  TH2D *mz = EgmMap({SampleIndex("DY")}, {EgmTnpMapName(true)});
  TH2D *mw = EgmMap({SampleIndex("Wp"), SampleIndex("Wm")}, {EgmMapName(true, 0), EgmMapName(true, 1)});
  if (!mz || !mw) return 2;
  EgmPrediction(mz, -1, 0, pred["ztnp"].sf, pred["ztnp"].err);
  std::cout << Form("[EGM] Z TnP probes                     SF_EGM = %.4f +- %.4f\n", pred["ztnp"].sf, pred["ztnp"].err);
  for (int s = 0; s < kNScheme; ++s)
    for (int k = 0; k < kSchemes[s].n; ++k)
    {
      Pred &p = pred[Form("%s_%d", kSchemes[s].name, k)];
      EgmPrediction(mw, s, k, p.sf, p.err);
      std::cout << Form("[EGM] W %-6s k%d (m_T > 40)              SF_EGM = %.4f +- %.4f\n", kSchemes[s].name, k, p.sf, p.err);
    }

  // ---------------- the record + the CSV ----------------
  const std::string csv = std::string(kPlotDir) + "/idiso_sf_summary.csv";
  std::ofstream o(csv);
  // boundary = no | r0 | r10 | r0+r10: a fail r of the bin at 0 / at its limit 10 (the SF is then
  // set by the POI range, not by the data, and neither its fit error nor the pull means anything)
  o << "kind,tag,bin,lo,hi,sf,sf_fit_err,sf_expected_err,n_fail,boundary,egm_pred,egm_err,pull_expected\n";
  auto record = [&](const char *kind, const Fit &f, const Row &r, const Pred &p) {
    const double e = f.Exp(r.bin);
    const double sig = std::sqrt((e > 0 ? e * e : r.sfErr * r.sfErr) + p.err * p.err);
    const double pull = sig > 0 ? (r.sf - p.sf) / sig : 0.0;
    std::cout << Form("[SF] %-24s bin %-4s SF %.3f +- %.3f (expected +- %.3f)  N_fail %6.0f%s\n", f.tag.c_str(),
                      r.bin.c_str(), r.sf, r.sfErr, e, r.nFail, r.Boundary() ? Form("   [%s]", BndText(r.bnd)) : "");
    std::cout << Form("[PULL] %-24s bin %-4s (SF - SF_EGM) / sqrt(expected^2 + EGM^2) = %+.2f%s\n", f.tag.c_str(),
                      r.bin.c_str(), pull, r.Boundary() ? "   (not meaningful: a fail r at a limit)" : "");
    o << kind << "," << f.tag << "," << r.bin << "," << r.lo << "," << r.hi << "," << r.sf << "," << r.sfErr << "," << e
      << "," << r.nFail << "," << BndCode(r.bnd) << "," << p.sf << "," << p.err << "," << pull << "\n";
  };
  if (const Row *r = z.Get("ztnp")) record("Z", z, *r, pred["ztnp"]);
  for (int s = 0; s < kNScheme; ++s)
    for (int iw = 0; iw < 3; ++iw)
    {
      const Fit &f = w[3 * s + iw];
      for (int k = 0; k < kSchemes[s].n; ++k)
        if (const Row *r = f.Get(std::to_string(k))) record("W", f, *r, pred[Form("%s_%d", kSchemes[s].name, k)]);
    }
  o.close();
  std::cout << "[OUT] " << csv << "\n";

  // ---------------- the plots ----------------
  DrawSummary(w, z, pred["ztnp"], pred["incl_0"]);
  DrawEgmOverlay(1, w, z.Get("ztnp"), pred);
  DrawEgmOverlay(2, w, z.Get("ztnp"), pred);
  for (int s = 0; s < kNScheme; ++s) DrawMultipliers(s, w);
  for (int iw = 0; iw < 3; ++iw) DrawFailShape(iw);
  return 0;
}
