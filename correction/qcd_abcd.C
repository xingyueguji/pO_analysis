// correction/qcd_abcd.C
//
// Data-driven QCD (multijet / "low-MET") background for the W -> l nu selection
// via the ABCD method, MUON channel first.
//
// Replaces the old shape-extrapolation approach
// (correction/qcd_sideband_fit_and_extrapolate.C), which fitted a Rayleigh-like
// MET shape in three anti-iso bins and linearly extrapolated the shape
// parameters to the signal isolation -- it did not behave well at the available
// pO statistics.
//
// ----------------------------------------------------------------------------
// METHOD (per charge, inclusive in rapidity for this first look)
// ----------------------------------------------------------------------------
// Two ~uncorrelated handles: lepton isolation (relIso) x a recoil-side variable y.
// NAMING: the SAME 2x2 machinery is run over more than one y (see `planes` in the
// driver): y = PF MET (the plane whose T normalizes everything) and y = m_T. A
// third relIso x lepton-pT plane supplies SHAPE ONLY and is never counted with the
// 2x2. So read "y high/low" below as MET high/low on the met pass and m_T high/low
// on the mt pass; the config field that splits it is ABCDConfig::yCut.
//
//            y high (>= yCut)       y low (< yCut)
//   iso    +---------------------+---------------------+
//   pass   |          A          |          B          |   relIso < isoCut
//  (<isoCut)|   (signal region)   |   (QCD-enriched)    |
//          +---------------------+---------------------+
//   iso    |          C          |          D          |   isoFailLo <= relIso
//   fail   |    (QCD-enriched)   |   (QCD-dominated)   |             < isoFailHi
//          +---------------------+---------------------+
//
// If relIso _|_ y for QCD, the 2x2 table factorizes:  A = B * C / D.
// Each region count is QCD-ONLY: N_QCD = N_data - N_EWK(MC), where the EWK
// (W/DY/tau) prediction is the ABSOLUTELY normalized MC (k_s = A_O*sigma*L/N_gen
// from skim/mc_norm.h). The anti-iso (fail) window is opened with a GAP above
// isoCut to keep real-W leakage out of the sideband.
//
// Deliverables:
//   - The transfer factor  T = N_QCD(B) / N_QCD(D)  (iso pass/fail at low y,
//     where QCD dominates and the EWK subtraction is small).
//   - The QCD yield in the high-y signal region  A_pred = B * C / D.
//   - The QCD TEMPLATE in the iso-pass region, one per plane (y = MET and y =
//     m_T) = the anti-iso shape (data - EWK, all y) scaled by that plane's T. By
//     construction it reproduces B at low y and predicts A at high y.
//   - A closure/overlay: data(iso-pass) vs [EWK + ABCD-QCD] with a ratio pad.
//
// The 2D inputs come from the skim (h_iso_met_mu{Plus,Minus},
// h_iso_mt_mu{Plus,Minus}); regions are chosen here by PROJECTION, so the
// isoCut / isoFail window / metCut can be retuned WITHOUT re-skimming.
//
// LEPTON-pT plane (2026-07-29): a third template qcd_pt_<lep><chg> for the
// lepton-pT stacks in plotting/mtandmet.C, built from h_iso_pt_<lep><chg> as
// (anti-iso pT shape) x (the MET-plane T). relIso has pT in its denominator,
// so relIso x pT is correlated for QCD -- the pT plane is NEVER counted with
// the 2x2 factorization; it only contributes the shape, and an anti-iso
// slice-stability plot (antiiso_shape_pt_*) checks that shape for bias.
// 2026-07-30: the same is also produced WITH the lepton-pT-discriminant m_T >
// 40 cut -> qcd_pt_mt40_<lep><chg> (from h_iso_pt_mt40_*), i.e. the QCD
// template of the "pT discriminant, pT>25 && m_T>40" selection. Same MET-plane
// T (QCD isolation efficiency assumed independent of the recoil variable), so
// the template total IS the ABCD prediction for QCD surviving m_T > 40.
//
// Inputs : ../skim/rootfile/WToMuNu_pO_PFMet_hist.root        (data)
//          ../skim/rootfile/WToMuNu_pO_PFMet_{Wp,Wm,DY,DYtau,Wptau,Wmtau}_hist.root
// Output : ./plots/qcd_abcd_mu/   + ./rootfile/qcd_abcd_mu.root (templates)
//
// Run from correction/:  root -l -q 'qcd_abcd.C+'           (muon)
//                        root -l -q 'qcd_abcd.C+(true)'     (electron)

#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../plotting/plotting_helper.C"
#include "../skim/mc_norm.h" // pONorm::MCScale -> per-sample k_s = A*sigma*L/N_gen

namespace
{
// ---------------- ABCD configuration (all tunable; projections, no re-skim) --
struct ABCDConfig
{
  double isoCut    = 0.15; // signal: relIso < isoCut          (regions A, B)
  double isoFailLo = 0.30; // anti-iso sideband, low edge      (regions C, D)
  double isoFailHi = 1.00; // anti-iso sideband, high edge
  // y-axis split: high (>=) vs low (<). Named yCut, NOT metCut: the same value
  // splits whichever y the plane carries (MET on the met pass, m_T on the mt pass).
  double yCut      = 30.0;
};

struct Count
{
  double n = 0.0, e = 0.0;
};

// QCD-only counts (data - EWK) in the four ABCD regions, + the raw data and EWK
// totals in the signal region A for context.
struct ABCDResult
{
  Count A, B, C, D;      // QCD-only
  double T = 0, Terr = 0;        // transfer factor B/D
  double Apred = 0, Aprederr = 0; // predicted QCD in signal region = B*C/D
  Count dataA, ewkA;     // raw data and EWK in A (for the contamination fraction)
  double qcdTotPass = 0, qcdTotPassErr = 0; // total iso-pass QCD = template integral
};

// Integral + error of a TH2 over [xlo,xhi) x [ylo,yhi). A sentinel hi >= axis
// max includes the overflow on that axis (so the "y high" region reaches the tail
// of whichever y the plane carries: MET, m_T or lepton pT).
Count count2D(TH2 *h, double xlo, double xhi, double ylo, double yhi)
{
  Count c;
  if (!h) return c;
  const TAxis *ax = h->GetXaxis();
  const TAxis *ay = h->GetYaxis();
  const double eps = 1e-6;
  const int bxlo = (xlo <= ax->GetXmin()) ? 1 : ax->FindBin(xlo + eps);
  const int bxhi = (xhi >= ax->GetXmax()) ? ax->GetNbins() : ax->FindBin(xhi - eps);
  const int bylo = (ylo <= ay->GetXmin()) ? 1 : ay->FindBin(ylo + eps);
  const int byhi = (yhi >= ay->GetXmax()) ? ay->GetNbins() + 1 // include overflow
                                          : ay->FindBin(yhi - eps);
  double err = 0.0;
  c.n = h->IntegralAndError(bxlo, bxhi, bylo, byhi, err);
  c.e = err;
  return c;
}

// ProjectionY (the y distribution of the plane: MET, m_T or lepton pT) over the
// x-window (relIso) [xlo,xhi).
TH1D *projY(TH2 *h, double xlo, double xhi, const char *name)
{
  if (!h) return nullptr;
  const TAxis *ax = h->GetXaxis();
  const double eps = 1e-6;
  const int bxlo = (xlo <= ax->GetXmin()) ? 1 : ax->FindBin(xlo + eps);
  const int bxhi = (xhi >= ax->GetXmax()) ? ax->GetNbins() : ax->FindBin(xhi - eps);
  TH1D *p = h->ProjectionY(name, bxlo, bxhi);
  p->SetDirectory(nullptr);
  return p;
}

double relErr(const Count &c) { return (c.n != 0.0) ? c.e / c.n : 0.0; }

// EWK (W/DY/tau) absolute-normalized sum of one named 2D histogram across all MC
// files. Returns a fresh owned TH2D (caller deletes). Each file contributes
// k_s * its histogram, with k_s from mc_norm.h.
struct MCFile { std::string label; TFile *f; };

TH2D *sumEWK2D(const std::vector<MCFile> &mc, const char *hname, const char *outname)
{
  TH2D *sum = nullptr;
  for (const MCFile &m : mc)
  {
    TH2D *h = (TH2D *)m.f->Get(hname);
    if (!h)
    {
      std::cerr << "[WARN] sumEWK2D: missing " << hname << " in " << m.label << "\n";
      continue;
    }
    const double k = pONorm::MCScale(m.label);
    if (!sum)
    {
      sum = (TH2D *)h->Clone(outname);
      sum->SetDirectory(nullptr);
      sum->Reset("ICESM"); // zero contents+errors, keep binning & Sumw2
    }
    sum->Add(h, k); // adds k*h with error propagation
  }
  return sum;
}

// Draw a QCD-only 2D (colz) with the ABCD region boundaries overlaid.
// drawYSplit=false suppresses the horizontal yCut line (the lepton-pT plane has
// no y-axis region split -- its normalization comes from the MET plane).
void draw2D(TH2 *h, const ABCDConfig &cfg, const std::string &outNoExt,
            const std::string &title, const char *ytit, bool drawYSplit = true)
{
  if (!h) return;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0); // title drawn via TLatex below; avoid the doubled pavetext
  TCanvas *c = new TCanvas(Form("c2d_%s", h->GetName()), "", 850, 750);
  c->SetRightMargin(0.15);
  c->SetLeftMargin(0.13);
  h->SetTitle("");
  h->GetXaxis()->SetTitle("relIso");
  h->GetYaxis()->SetTitle(ytit);
  h->GetXaxis()->SetRangeUser(0.0, cfg.isoFailHi);
  // Clamp the z (color) range: a single over-subtracted W-peak bin (relIso~0,
  // y~40 -- the W peak of whichever y this plane carries) otherwise dominates the
  // scale and washes out the actual QCD structure. Floor at -0.3*max so that
  // outlier just saturates the low color.
  const double zmax = h->GetMaximum();
  if (zmax > 0) { h->SetMaximum(zmax); h->SetMinimum(-0.3 * zmax); }
  h->Draw("COLZ");

  auto vline = [&](double x) {
    TLine *l = new TLine(x, h->GetYaxis()->GetXmin(), x, h->GetYaxis()->GetXmax());
    l->SetLineColor(kRed + 1); l->SetLineWidth(2); l->SetLineStyle(2); l->Draw();
  };
  auto hline = [&](double y, double x1, double x2) {
    TLine *l = new TLine(x1, y, x2, y);
    l->SetLineColor(kRed + 1); l->SetLineWidth(2); l->SetLineStyle(2); l->Draw();
  };
  vline(cfg.isoCut);
  vline(cfg.isoFailLo);
  if (drawYSplit) hline(cfg.yCut, 0.0, cfg.isoFailHi);

  TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.03);
  t.DrawLatex(0.14, 0.92, title.c_str());
  c->SaveAs((outNoExt + ".png").c_str());
  c->SaveAs((outNoExt + ".pdf").c_str());
  delete c;
}

// ----------------------------------------------------------------------------
// One ABCD plane (MET or MT), one charge. Builds the QCD-only 2D, fills the
// ABCD result, builds the iso-pass QCD template (returned, owned), and draws
// the 2D + the closure overlay.
// ----------------------------------------------------------------------------
TH1D *runPlaneCharge(TFile *fData, const std::vector<MCFile> &mc,
                     const std::string &plane, const std::string &chg,
                     const ABCDConfig &cfg, const std::string &outDir,
                     const std::string &chgLatex, const std::string &ytit,
                     const std::string &lep, ABCDResult &res)
{
  const std::string hname = Form("h_iso_%s_%s%s", plane.c_str(), lep.c_str(), chg.c_str());
  TH2D *data2D = (TH2D *)fData->Get(hname.c_str());
  if (!data2D)
  {
    std::cerr << "[ERR] missing " << hname << " in data\n";
    return nullptr;
  }
  data2D = (TH2D *)data2D->Clone(Form("%s_data", hname.c_str()));
  data2D->SetDirectory(nullptr);

  TH2D *ewk2D = sumEWK2D(mc, hname.c_str(), Form("%s_ewk", hname.c_str()));

  // QCD-only 2D = data - EWK
  TH2D *qcd2D = (TH2D *)data2D->Clone(Form("%s_qcd", hname.c_str()));
  qcd2D->SetDirectory(nullptr);
  if (ewk2D) qcd2D->Add(ewk2D, -1.0);

  // --- region counts (QCD-only) ---
  res.A = count2D(qcd2D, 0.0, cfg.isoCut,            cfg.yCut, 1e9);
  res.B = count2D(qcd2D, 0.0, cfg.isoCut,            0.0,        cfg.yCut);
  res.C = count2D(qcd2D, cfg.isoFailLo, cfg.isoFailHi, cfg.yCut, 1e9);
  res.D = count2D(qcd2D, cfg.isoFailLo, cfg.isoFailHi, 0.0,        cfg.yCut);

  // context: raw data and EWK in the signal region A
  res.dataA = count2D(data2D, 0.0, cfg.isoCut, cfg.yCut, 1e9);
  res.ewkA  = ewk2D ? count2D(ewk2D, 0.0, cfg.isoCut, cfg.yCut, 1e9) : Count{};

  // transfer factor and prediction
  res.T = (res.D.n != 0.0) ? res.B.n / res.D.n : 0.0;
  res.Terr = (res.B.n != 0.0 && res.D.n != 0.0)
                 ? std::fabs(res.T) * std::sqrt(relErr(res.B) * relErr(res.B) +
                                                relErr(res.D) * relErr(res.D))
                 : 0.0;
  res.Apred = (res.D.n != 0.0) ? res.B.n * res.C.n / res.D.n : 0.0;
  res.Aprederr = (res.B.n && res.C.n && res.D.n)
                     ? std::fabs(res.Apred) *
                           std::sqrt(relErr(res.B) * relErr(res.B) +
                                     relErr(res.C) * relErr(res.C) +
                                     relErr(res.D) * relErr(res.D))
                     : 0.0;

  // --- QCD template in the iso-pass region = (anti-iso shape) x T ---
  TH1D *tmpl = projY(qcd2D, cfg.isoFailLo, cfg.isoFailHi,
                     Form("qcd_%s_%s%s", plane.c_str(), lep.c_str(), chg.c_str()));
  for (int i = 1; i <= tmpl->GetNbinsX(); ++i)
  {
    const double v = tmpl->GetBinContent(i);
    const double ev = tmpl->GetBinError(i);
    // Clamp tiny negative high-y tail bins (data-EWK fluctuations) to 0 so the
    // QCD template is a physical, non-negative background (also keeps THStack /
    // Combine happy). The clamped yield is negligible (~0.03% of the integral).
    tmpl->SetBinContent(i, std::max(0.0, v * res.T));
    // fold the shape stat error and the T uncertainty
    tmpl->SetBinError(i, std::sqrt((res.T * ev) * (res.T * ev) +
                                   (v * res.Terr) * (v * res.Terr)));
  }
  double tot_e = 0.0;
  res.qcdTotPass = tmpl->IntegralAndError(1, tmpl->GetNbinsX(), tot_e);
  res.qcdTotPassErr = tot_e;

  // --- diagnostics: 2D colz with boundaries ---
  draw2D(qcd2D, cfg, outDir + Form("/qcd2D_%s_%s%s", plane.c_str(), lep.c_str(), chg.c_str()),
         Form("QCD-only (data-EWK), %s", chgLatex.c_str()), ytit.c_str());

  // --- closure overlay in the iso-pass region: data vs EWK(+QCD) ---
  TH1D *passData = projY(data2D, 0.0, cfg.isoCut,
                         Form("pass_data_%s_%s%s", plane.c_str(), lep.c_str(), chg.c_str()));
  TH1D *passEWK = ewk2D ? projY(ewk2D, 0.0, cfg.isoCut,
                                Form("pass_ewk_%s_%s%s", plane.c_str(), lep.c_str(), chg.c_str()))
                        : nullptr;

  // Layout: SaveDataMCRatio puts CMS/WiP top-left and its legend top-right.
  // Unlike the Z-peak kinematic plots, these log-y W distributions fill the
  // top-LEFT too (the low-MET QCD excess), so the header goes BELOW the legend on
  // the right -- the one genuinely empty patch in the high-MET/m_T tail. Keep the
  // header short (2 lines); "Data vs ..." is already conveyed by the legend.
  PlotStyle ps;
  ps.logy = true;
  ps.yTitleOffset = 1.55;
  ps.headerX = 0.64;
  ps.headerY = 0.66;
  ps.titleSize = 0.04;

  if (passEWK)
  {
    // model WITHOUT QCD (shows the low-y excess the QCD fills)
    SaveDataMCRatio(passData, passEWK,
                    outDir + Form("/closure_noQCD_%s_%s%s", plane.c_str(), lep.c_str(), chg.c_str()),
                    ytit, "Events", Form("iso-pass, %s", chgLatex.c_str()),
                    "no QCD", "", ps, false, "Data", "EWK MC");

    // model WITH ABCD QCD added
    TH1D *model = (TH1D *)passEWK->Clone(Form("model_%s_%s%s", plane.c_str(), lep.c_str(), chg.c_str()));
    model->SetDirectory(nullptr);
    model->Add(tmpl);
    SaveDataMCRatio(passData, model,
                    outDir + Form("/closure_withQCD_%s_%s%s", plane.c_str(), lep.c_str(), chg.c_str()),
                    ytit, "Events", Form("iso-pass, %s", chgLatex.c_str()),
                    "EWK + ABCD QCD", "", ps, false, "Data", "EWK+QCD");
    delete model;
  }

  delete passData;
  if (passEWK) delete passEWK;
  delete data2D;
  if (ewk2D) delete ewk2D;
  delete qcd2D;
  return tmpl; // caller writes & deletes
}

// ----------------------------------------------------------------------------
// The lepton-pT plane, one charge (2026-07-29). relIso carries pT in its
// denominator, so relIso x pT is CORRELATED for QCD and must NOT be counted
// with the 2x2 factorization. This plane only supplies the anti-iso pT SHAPE;
// the normalization is the MET-plane transfer factor metRes.T. Since the
// anti-iso pT projection integrates the same events as the anti-iso MET
// projection, the pT template total equals the MET template total (same
// iso-pass QCD yield) by construction.
// Also draws the shape-STABILITY check the MET template never needed: the
// QCD-only anti-iso pT shape in a low vs high anti-iso slice (split at the
// window midpoint), shape-normalized -- if the two disagree, the anti-iso pT
// shape is biased by the iso-pT correlation and the window needs retuning.
// ----------------------------------------------------------------------------
// `tag` selects the phase space: "" = the nominal iso-pass selection,
// "mt40" = the same WITH the lepton-pT-discriminant m_T > 40 cut (reads
// h_iso_pt_mt40_*, writes qcd_pt_mt40_*). The transfer factor is the SAME
// MET-plane T in both cases -- the ABCD factorization assumes QCD isolation
// efficiency is independent of the recoil-side variable, so the m_T>40
// template total is the ABCD prediction for QCD in that region.
TH1D *runPtCharge(TFile *fData, const std::vector<MCFile> &mc,
                  const std::string &chg, const ABCDConfig &cfg,
                  const std::string &outDir, const std::string &chgLatex,
                  const std::string &lep, const ABCDResult &metRes,
                  const std::string &tag = "")
{
  const std::string sfx   = tag.empty() ? "" : ("_" + tag);
  const std::string hname = Form("h_iso_pt%s_%s%s", sfx.c_str(), lep.c_str(), chg.c_str());
  TH2D *data2D = (TH2D *)fData->Get(hname.c_str());
  if (!data2D)
  {
    std::cerr << "[WARN] missing " << hname << " in data (skim output predates the"
              << " lepton-pT ABCD plane; re-skim to enable the pT QCD template)\n";
    return nullptr;
  }
  data2D = (TH2D *)data2D->Clone(Form("%s_data", hname.c_str()));
  data2D->SetDirectory(nullptr);

  TH2D *ewk2D = sumEWK2D(mc, hname.c_str(), Form("%s_ewk", hname.c_str()));

  TH2D *qcd2D = (TH2D *)data2D->Clone(Form("%s_qcd", hname.c_str()));
  qcd2D->SetDirectory(nullptr);
  if (ewk2D) qcd2D->Add(ewk2D, -1.0);

  const std::string ytit = Form("p_{T}^{%s} [GeV]", lep == "ele" ? "e" : "#mu");

  // --- QCD template in the iso-pass region = (anti-iso pT shape) x T_MET ---
  TH1D *tmpl = projY(qcd2D, cfg.isoFailLo, cfg.isoFailHi,
                     Form("qcd_pt%s_%s%s", sfx.c_str(), lep.c_str(), chg.c_str()));
  for (int i = 1; i <= tmpl->GetNbinsX(); ++i)
  {
    const double v = tmpl->GetBinContent(i);
    const double ev = tmpl->GetBinError(i);
    tmpl->SetBinContent(i, std::max(0.0, v * metRes.T));
    tmpl->SetBinError(i, std::sqrt((metRes.T * ev) * (metRes.T * ev) +
                                   (v * metRes.Terr) * (v * metRes.Terr)));
  }

  // --- diagnostics: 2D colz (no y-split line -- normalization is MET-plane) ---
  draw2D(qcd2D, cfg, outDir + Form("/qcd2D_pt%s_%s%s", sfx.c_str(), lep.c_str(), chg.c_str()),
         Form("QCD-only (data-EWK), %s", chgLatex.c_str()), ytit.c_str(),
         /*drawYSplit=*/false);

  // --- closure overlay in the iso-pass region: data vs EWK(+QCD) ---
  TH1D *passData = projY(data2D, 0.0, cfg.isoCut,
                         Form("pass_data_pt%s_%s%s", sfx.c_str(), lep.c_str(), chg.c_str()));
  TH1D *passEWK = ewk2D ? projY(ewk2D, 0.0, cfg.isoCut,
                                Form("pass_ewk_pt%s_%s%s", sfx.c_str(), lep.c_str(), chg.c_str()))
                        : nullptr;

  PlotStyle ps;
  ps.logy = true;
  ps.yTitleOffset = 1.55;
  ps.headerX = 0.64;
  ps.headerY = 0.66;
  ps.titleSize = 0.04;

  if (passEWK)
  {
    SaveDataMCRatio(passData, passEWK,
                    outDir + Form("/closure_noQCD_pt%s_%s%s", sfx.c_str(), lep.c_str(), chg.c_str()),
                    ytit, "Events", Form("iso-pass, %s", chgLatex.c_str()),
                    "no QCD", "", ps, false, "Data", "EWK MC");

    TH1D *model = (TH1D *)passEWK->Clone(Form("model_pt%s_%s%s", sfx.c_str(), lep.c_str(), chg.c_str()));
    model->SetDirectory(nullptr);
    model->Add(tmpl);
    SaveDataMCRatio(passData, model,
                    outDir + Form("/closure_withQCD_pt%s_%s%s", sfx.c_str(), lep.c_str(), chg.c_str()),
                    ytit, "Events", Form("iso-pass, %s", chgLatex.c_str()),
                    "EWK + ABCD QCD", "", ps, false, "Data", "EWK+QCD");
    delete model;
  }

  // --- anti-iso slice stability: low vs high half of the sideband, shape-only ---
  const double mid = 0.5 * (cfg.isoFailLo + cfg.isoFailHi);
  TH1D *sLo = projY(qcd2D, cfg.isoFailLo, mid,
                    Form("antiiso_lo_pt%s_%s%s", sfx.c_str(), lep.c_str(), chg.c_str()));
  TH1D *sHi = projY(qcd2D, mid, cfg.isoFailHi,
                    Form("antiiso_hi_pt%s_%s%s", sfx.c_str(), lep.c_str(), chg.c_str()));
  if (sLo && sHi)
  {
    PlotStyle ps2;
    ps2.logy = false; // shape comparison reads better linear
    ps2.yTitleOffset = 1.55;
    ps2.headerX = 0.60;
    ps2.headerY = 0.66;
    ps2.titleSize = 0.04;
    // NB both inputs are DATA-derived (EWK-subtracted sideband slices) -- the
    // "data" point style / "MC" fill style are just the helper's cosmetics, so
    // relabel the ratio pad accordingly (low slice / high slice).
    SaveDataMCRatio(sLo, sHi,
                    outDir + Form("/antiiso_shape_pt%s_%s%s", sfx.c_str(), lep.c_str(), chg.c_str()),
                    ytit, "Events (a.u.)",
                    Form("anti-iso pT shape, %s", chgLatex.c_str()),
                    "slice stability check", "shape-normalized", ps2,
                    /*normToData=*/true,
                    Form("QCD, relIso [%.2f,%.2f)", cfg.isoFailLo, mid),
                    Form("QCD, relIso [%.2f,%.2f)", mid, cfg.isoFailHi),
                    "low / high slice");
  }
  delete sLo;
  delete sHi;

  delete passData;
  if (passEWK) delete passEWK;
  delete data2D;
  if (ewk2D) delete ewk2D;
  delete qcd2D;
  return tmpl; // caller writes & deletes
}

void printRow(const char *tag, const Count &c)
{
  std::cout << "    " << tag << " = " << Form("%10.1f +/- %7.1f", c.n, c.e) << "\n";
}

// Human-readable name of a plane's y variable, so the region labels below say
// what was ACTUALLY split (they used to hardcode "MET" even on the m_T pass).
std::string YName(const std::string &plane)
{
  if (plane == "met") return "MET";
  if (plane == "mt")  return "m_T";
  if (plane == "pt")  return "pT";
  return plane;
}

void printResult(const std::string &plane, const std::string &lep,
                 const std::string &chg, const ABCDResult &r)
{
  const std::string y = YName(plane);
  std::cout << "\n=== ABCD QCD estimate  [iso x " << plane << "]  " << lep << " " << chg << " ===\n";
  std::cout << "  QCD-only region counts (data - EWK), y = " << y << ":\n";
  printRow(Form("A  iso-pass, %s-high (SIGNAL)", y.c_str()), r.A);
  printRow(Form("B  iso-pass, %s-low         ", y.c_str()), r.B);
  printRow(Form("C  iso-fail, %s-high        ", y.c_str()), r.C);
  printRow(Form("D  iso-fail, %s-low         ", y.c_str()), r.D);
  std::cout << "  Transfer factor T = B/D            = "
            << Form("%.4f +/- %.4f", r.T, r.Terr) << "\n";
  std::cout << "  Predicted QCD in signal A = B*C/D  = "
            << Form("%.1f +/- %.1f", r.Apred, r.Aprederr) << "\n";
  std::cout << "  Total iso-pass QCD (template norm) = "
            << Form("%.1f +/- %.1f", r.qcdTotPass, r.qcdTotPassErr) << "\n";
  std::cout << "  --- context (signal region A) ---\n";
  std::cout << "    data(A) = " << Form("%.1f", r.dataA.n)
            << " ,  EWK(A) = " << Form("%.1f", r.ewkA.n)
            << " ,  QCD(A)/data(A) = "
            << Form("%.1f%%", r.dataA.n > 0 ? 100.0 * r.Apred / r.dataA.n : 0.0) << "\n";
}
} // namespace

void qcd_abcd(bool isElec = false)
{
  // ---- channel-specific config -------------------------------------------
  // Muon and electron share the machinery; only the file base, the MCScale
  // labels, the lepton tag in the histogram names, and the isolation window
  // differ. The electron signal iso cut is tighter (0.095 vs 0.15), so the
  // anti-iso sideband opens higher and the gap shifts accordingly. All of these
  // are projections -> retune freely with no re-skim.
  const std::string lep      = isElec ? "ele" : "mu";
  const std::string fileBase = isElec ? "WToElecNu_pO_PFMet" : "WToMuNu_pO_PFMet";
  const std::string lepP     = isElec ? "e^{+}" : "#mu^{+}";
  const std::string lepM     = isElec ? "e^{-}" : "#mu^{-}";

  ABCDConfig cfg; // muon defaults (isoCut 0.15, isoFail [0.30,1.0), yCut 30)
  if (isElec)
  {
    cfg.isoCut    = 0.095; // electron signal isolation cut (skim_Wel isoMax)
    cfg.isoFailLo = 0.20;  // anti-iso sideband low edge (gap 0.095-0.20)
    cfg.isoFailHi = 1.00;
    // cfg.yCut stays 30 GeV
  }

  const std::string dir = "../skim/rootfile/";
  const std::string dataFile = dir + fileBase + "_hist.root";
  TFile *fData = TFile::Open(dataFile.c_str(), "READ");
  if (!fData || fData->IsZombie())
  {
    std::cerr << "[ERR] cannot open data: " << dataFile << "\n";
    return;
  }

  // EWK MC: MCScale label -> per-sample skim suffix. The W/DY labels carry the
  // lepton flavour (Wp_mu/Wp_ele, DYmu/DYee); the tau labels are shared (the tau
  // N_gen is per physical file, channel-independent).
  const std::string dyLab = isElec ? "DYee" : "DYmu";
  struct Spec { std::string label, suffix; };
  const std::vector<Spec> spec = {
      {"Wp_" + lep, "Wp"},
      {"Wm_" + lep, "Wm"},
      {dyLab,       "DY"},
      {"DYtau",     "DYtau"},
      {"Wp_tau",    "Wptau"},
      {"Wm_tau",    "Wmtau"},
  };
  std::vector<MCFile> mc;
  for (const Spec &s : spec)
  {
    const std::string fn = dir + fileBase + "_" + s.suffix + "_hist.root";
    TFile *f = TFile::Open(fn.c_str(), "READ");
    if (!f || f->IsZombie())
    {
      std::cerr << "[ERR] cannot open MC: " << fn << "\n";
      return;
    }
    mc.push_back({s.label, f});
  }

  const std::string outDir  = "./plots/qcd_abcd_" + lep;
  const std::string outFile = "./rootfile/qcd_abcd_" + lep + ".root";
  gSystem->mkdir(outDir.c_str(), kTRUE);
  gSystem->mkdir("./rootfile", kTRUE);

  std::cout << "\n[CONFIG] channel=" << lep << "  isoCut=" << cfg.isoCut
            << "  isoFail=[" << cfg.isoFailLo << "," << cfg.isoFailHi << ")"
            << "  yCut=" << cfg.yCut << " GeV"
            << "  (the y split applies to EACH counted plane: MET and m_T)\n";

  TFile *fout = TFile::Open(outFile.c_str(), "RECREATE");

  struct Plane { std::string key; std::string ytit; };
  const std::vector<Plane> planes = {{"met", "PF MET (GeV)"}, {"mt", "m_{T} (GeV)"}};
  const std::vector<std::pair<std::string, std::string>> charges = {
      {"Plus", lepP}, {"Minus", lepM}};

  // The MET-plane results are kept per charge: the lepton-pT plane below
  // borrows its transfer factor (relIso x pT is correlated for QCD, so the pT
  // plane must not be counted with the 2x2 factorization itself).
  ABCDResult metResPlus, metResMinus;

  for (const Plane &pl : planes)
  {
    for (const auto &cq : charges)
    {
      ABCDResult res;
      TH1D *tmpl = runPlaneCharge(fData, mc, pl.key, cq.first, cfg, outDir,
                                  cq.second, pl.ytit, lep, res);
      if (tmpl)
      {
        printResult(pl.key, lep, cq.first, res);
        if (pl.key == "met")
          ((cq.first == "Plus") ? metResPlus : metResMinus) = res;
        fout->cd();
        tmpl->Write(tmpl->GetName(), TObject::kOverwrite);
        delete tmpl;
      }
    }
  }

  // --- lepton-pT plane: anti-iso SHAPE x the MET-plane T (template only) ---
  for (const auto &cq : charges)
  {
    const ABCDResult &mres = (cq.first == "Plus") ? metResPlus : metResMinus;
    if (mres.T == 0.0)
    {
      std::cerr << "[WARN] no MET-plane transfer factor for " << cq.first
                << "; skipping the lepton-pT QCD template\n";
      continue;
    }
    // "" = nominal iso-pass; "mt40" = + the lepton-pT-discriminant m_T > 40 cut
    for (const std::string &tag : {std::string(""), std::string("mt40")})
    {
      TH1D *tmpl = runPtCharge(fData, mc, cq.first, cfg, outDir, cq.second, lep, mres, tag);
      if (!tmpl) continue;
      double te = 0.0;
      const double tt = tmpl->IntegralAndError(1, tmpl->GetNbinsX(), te);
      std::cout << "\n=== ABCD QCD lepton-pT template  " << lep << " " << cq.first
                << (tag.empty() ? "  [no m_T cut]" : "  [m_T > 40]") << " ===\n"
                << "  (anti-iso pT shape x MET-plane T = "
                << Form("%.4f +/- %.4f", mres.T, mres.Terr) << ")\n"
                << "  Total iso-pass QCD (template norm) = "
                << Form("%.1f +/- %.1f", tt, te);
      if (tag.empty())
        std::cout << "  [MET-plane value: " << Form("%.1f", mres.qcdTotPass) << "]\n";
      else
        std::cout << "  -> m_{T}>40 keeps "
                  << Form("%.1f%%", mres.qcdTotPass > 0 ? 100.0 * tt / mres.qcdTotPass : 0.0)
                  << " of the QCD\n";
      fout->cd();
      tmpl->Write(tmpl->GetName(), TObject::kOverwrite);
      delete tmpl;
    }
  }

  fout->Close();
  delete fout;
  fData->Close();
  delete fData;
  for (MCFile &m : mc) { m.f->Close(); delete m.f; }

  std::cout << "\n[DONE] ABCD QCD (" << lep << ") -> plots: " << outDir
            << "/   templates: " << outFile << "\n";
}
