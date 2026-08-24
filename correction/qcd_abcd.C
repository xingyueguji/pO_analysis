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
// IN-FIT ABCD (2026-08-23, QCD_MODE=abcd in the Combine fork): the leppt_mt40
// QCD NORMALIZATION moves into the likelihood -- counting CR channels B
// (iso-pass, m_T < yCut), C40 (anti-iso, m_T > 40 = the SR's own cut) and D
// (anti-iso, m_T < yCut) on the m_T plane, SR total = (sB*sC/sD) * B0*C40/D0
// via a functional rateParam, so the EWK subtraction floats with the fitted
// r's instead of being frozen at r = 1 (the r-scan below shows that prefit
// assumption is what made T_MET move ~20% for muons). The SHAPE is unchanged:
// still the anti-iso m_T>40 pT distribution (qcd_pt_mt40_*). This macro
// supplies the machine-readable prefit inputs (abcd_counts_<lep><chg>,
// labeled-bin TH1D in the output rootfile), the in-fit window-variation
// systematic, the REDUCED lnN kappa (window + tilt only: the stat rows are
// profiled in-fit and the plane-transport row does not apply -- the m_T plane
// IS the counting plane), and the per-pT-bin fake-factor diagnostic ff_pt_*
// (the flat-T shape-assumption check; see runFFCheck). The MET-plane template
// machinery and the original kappa budget are UNCHANGED (QCD_MODE=lnN|free
// still consume them).
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
#include <sstream>
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

// The lepton-pT discriminant's m_T cut (= mtCutForPtDisc in the skim). The
// in-fit ABCD C region and the fake-factor diagnostic re-cut the y axis here;
// a bin edge on both the m_T (2.5 GeV) and MET (2 GeV) axes. The [yCut, 40)
// band stays an unused buffer (Jacobian protection for B/D).
const double kSRMtCut = 40.0;

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
  // Raw data and EWK in EVERY region -> the contamination (purity) table. The
  // EWK part is split W-related / DY-related, which the r-scan diagnostic needs
  // (only the W part scales with the signal strength).
  Count dataA, ewkA, dataB, ewkB, dataC, ewkC, dataD, ewkD;
  Count wB, zB, wD, zD;  // W- and DY-related EWK in B and D (for the r-scan)
  // In-fit ABCD (QCD_MODE=abcd): the C region re-cut at the SR's own m_T > 40
  // (kSRMtCut) so A40 = B*C40/D predicts the QCD of the leppt_mt40 selection
  // directly. Only meaningful on the mt pass (exported/printed there alone).
  Count C40, dataC40, wC40, zC40;   // QCD-only / raw data / W-rel / DY-rel EWK
  Count zB_dy, zB_dytau;            // DY vs DYtau split of zB (CR-B processes)
  double Apred40 = 0, Apred40err = 0; // B*C40/D = the abcd-mode template total
  double winSpread40 = 0; // max |dev| of A40 over the anti-iso window variations, %
  double qcdTotPass = 0, qcdTotPassErr = 0; // total iso-pass QCD = template integral
  double winSpread = 0;  // max |deviation| of Apred over the anti-iso window variations, %
  double subSlice  = 0;  // max |deviation| of T over the low-y sub-slices, %
  double clampedFrac = 0;// negative tail removed by the template clamp, %
  std::string diag;      // buffered diagnostics, emitted by printResult()
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

// Human-readable name of a plane's y variable, so the region labels below say
// what was ACTUALLY split (they used to hardcode "MET" even on the m_T pass).
std::string YName(const std::string &plane)
{
  if (plane == "met") return "MET";
  if (plane == "mt")  return "m_T";
  if (plane == "pt")  return "pT";
  return plane;
}

// ----------------------------------------------------------------------------
// Turn an (EWK-subtracted) anti-iso sideband shape into the iso-pass QCD
// template: clamp negative tail bins, multiply by the transfer factor T, and
// keep the per-bin errors as SHAPE STATISTICS ONLY.
//
// delta_T is a SINGLE, fully correlated normalisation uncertainty and must not
// be spread over the bins: TH1::IntegralAndError would then add it in quadrature
// across bins, i.e. use sum(v_i^2) where the correct total needs (sum v_i)^2,
// diluting it by the effective sqrt(N_bins) (a factor 1.8 on the mu+ m_T>40
// template: 4.3% quoted vs 7.8% correct).  It would also be picked up as N
// independent nuisances if autoMCStats is ever switched on in the datacards,
// double-counting the qcd_rate_<flav>_<chg> lnN that already carries it.
// So it is folded into the TOTAL only, and returned here.
struct TemplateNorm { double total = 0.0, err = 0.0; double clampedFrac = 0.0; };  // clampedFrac = negative tail removed, % of the template

TemplateNorm ScaleToIsoPass(TH1D *h, double T, double Terr)
{
  TemplateNorm out;
  if (!h) return out;
  const int n = h->GetNbinsX();
  // Clamp tiny negative high-y tail bins (data-EWK fluctuations) to 0 FIRST, so
  // the quoted total is the integral of the histogram actually written out. The
  // clamped yield is negligible (<=0.02%, printed per template). Bin errors are
  // kept (a clamped bin
  // still has a statistical uncertainty).
  double clamped = 0.0;
  for (int i = 1; i <= n; ++i)
    if (h->GetBinContent(i) < 0.0) { clamped += -h->GetBinContent(i); h->SetBinContent(i, 0.0); }

  // Integrate the SHAPE before scaling: V and its statistical error.
  double shapeErr = 0.0;
  const double V = h->IntegralAndError(1, n, shapeErr);

  for (int i = 1; i <= n; ++i)
  {
    h->SetBinContent(i, T * h->GetBinContent(i));
    h->SetBinError  (i, T * h->GetBinError(i)); // shape stat only
  }

  out.total = T * V;
  out.err   = std::sqrt((T * shapeErr) * (T * shapeErr) + (V * Terr) * (V * Terr));
  out.clampedFrac = (V > 0.0) ? 100.0 * clamped / V : 0.0;
  return out;
}

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
TH1D *runPlaneCharge(TFile *fData, const std::vector<MCFile> &mcW,
                     const std::vector<MCFile> &mcZ,
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

  TH2D *w2D   = sumEWK2D(mcW, hname.c_str(), Form("%s_w", hname.c_str()));
  TH2D *z2D   = sumEWK2D(mcZ, hname.c_str(), Form("%s_z", hname.c_str()));
  TH2D *ewk2D = nullptr;
  if (w2D) { ewk2D = (TH2D *)w2D->Clone(Form("%s_ewk", hname.c_str())); ewk2D->SetDirectory(nullptr); }
  if (ewk2D && z2D) ewk2D->Add(z2D);

  // QCD-only 2D = data - EWK
  TH2D *qcd2D = (TH2D *)data2D->Clone(Form("%s_qcd", hname.c_str()));
  qcd2D->SetDirectory(nullptr);
  if (ewk2D) qcd2D->Add(ewk2D, -1.0);

  // --- region counts (QCD-only) ---
  res.A = count2D(qcd2D, 0.0, cfg.isoCut,            cfg.yCut, 1e9);
  res.B = count2D(qcd2D, 0.0, cfg.isoCut,            0.0,        cfg.yCut);
  res.C = count2D(qcd2D, cfg.isoFailLo, cfg.isoFailHi, cfg.yCut, 1e9);
  res.D = count2D(qcd2D, cfg.isoFailLo, cfg.isoFailHi, 0.0,        cfg.yCut);

  // context: raw data and EWK in EVERY region (purity table + r-scan inputs)
  auto reg = [&](TH2 *h, int which) -> Count {
    if (!h) return Count{};
    const double xlo = (which < 2) ? 0.0 : cfg.isoFailLo;
    const double xhi = (which < 2) ? cfg.isoCut : cfg.isoFailHi;
    const bool   hi  = (which % 2 == 0);           // 0=A,1=B,2=C,3=D
    return count2D(h, xlo, xhi, hi ? cfg.yCut : 0.0, hi ? 1e9 : cfg.yCut);
  };
  res.dataA = reg(data2D, 0); res.ewkA = reg(ewk2D, 0);
  res.dataB = reg(data2D, 1); res.ewkB = reg(ewk2D, 1);
  res.dataC = reg(data2D, 2); res.ewkC = reg(ewk2D, 2);
  res.dataD = reg(data2D, 3); res.ewkD = reg(ewk2D, 3);
  res.wB = reg(w2D, 1); res.zB = reg(z2D, 1);
  res.wD = reg(w2D, 3); res.zD = reg(z2D, 3);

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

  // --- in-fit ABCD (QCD_MODE=abcd): C re-cut at the SR's own m_T > 40 -------
  // B and D stay at y < yCut (the [yCut, 40) band is a buffer), so
  // A40 = B*C40/D is the prefit QCD of the leppt_mt40 selection = the value
  // the CR scales sB*sC/sD multiply in the fit. Computed generically on every
  // plane; only the mt pass is exported (on the met pass "y > 40" is just a
  // different MET cut). The B-region DY is additionally split DY vs DYtau:
  // the in-fit CR-B carries them as separate processes (both ride r_Z; the
  // bookkeeping mirrors the SR process list).
  res.C40     = count2D(qcd2D,  cfg.isoFailLo, cfg.isoFailHi, kSRMtCut, 1e9);
  res.dataC40 = count2D(data2D, cfg.isoFailLo, cfg.isoFailHi, kSRMtCut, 1e9);
  res.wC40    = count2D(w2D,    cfg.isoFailLo, cfg.isoFailHi, kSRMtCut, 1e9);
  res.zC40    = count2D(z2D,    cfg.isoFailLo, cfg.isoFailHi, kSRMtCut, 1e9);
  res.Apred40 = (res.D.n != 0.0) ? res.B.n * res.C40.n / res.D.n : 0.0;
  res.Apred40err = (res.B.n && res.C40.n && res.D.n)
                       ? std::fabs(res.Apred40) *
                             std::sqrt(relErr(res.B) * relErr(res.B) +
                                       relErr(res.C40) * relErr(res.C40) +
                                       relErr(res.D) * relErr(res.D))
                       : 0.0;
  for (const MCFile &m : mcZ)
  {
    TH2D *h = (TH2D *)m.f->Get(hname.c_str());
    if (!h) continue;
    Count c = count2D(h, 0.0, cfg.isoCut, 0.0, cfg.yCut); // region B, raw MC
    const double k = pONorm::MCScale(m.label);
    c.n *= k;
    c.e *= k;
    Count &dst = (m.label == "DYtau") ? res.zB_dytau : res.zB_dy;
    dst.n += c.n;
    dst.e = std::sqrt(dst.e * dst.e + c.e * c.e);
  }

  // --- QCD template in the iso-pass region = (anti-iso shape) x T ---
  TH1D *tmpl = projY(qcd2D, cfg.isoFailLo, cfg.isoFailHi,
                     Form("qcd_%s_%s%s", plane.c_str(), lep.c_str(), chg.c_str()));
  const TemplateNorm tn = ScaleToIsoPass(tmpl, res.T, res.Terr);
  res.qcdTotPass    = tn.total;
  res.qcdTotPassErr = tn.err;
  res.clampedFrac   = tn.clampedFrac;

  // The two blocks below are BUFFERED, not printed here: runPlaneCharge runs
  // before printResult(), so writing straight to stdout would file them under
  // the PREVIOUS charge's header. printResult() emits res.diag at the end.
  std::ostringstream dg;
  // ------------------------------------------------------------------------
  // FACTORISATION TEST 1: is T stable WITHIN the low-y region?
  // Under exact relIso _|_ y factorisation, T measured in sub-slices of the
  // low-y band must be constant. The spread is a direct, data-driven measure of
  // the correlation the method assumes away.
  // ------------------------------------------------------------------------
  {
    dg << "  --- T = B/D in sub-slices of the " << YName(plane) << "-low region ---\n";
    const double h = cfg.yCut / 3.0, half = cfg.yCut / 2.0;
    const double edges[5][2] = {{0, h}, {h, 2 * h}, {2 * h, cfg.yCut},
                                {0, half}, {half, cfg.yCut}};
    double maxdev = 0.0;
    for (int i = 0; i < 5; ++i)
    {
      const Count b = count2D(qcd2D, 0.0, cfg.isoCut, edges[i][0], edges[i][1]);
      const Count d = count2D(qcd2D, cfg.isoFailLo, cfg.isoFailHi, edges[i][0], edges[i][1]);
      const double t = (d.n != 0.0) ? b.n / d.n : 0.0;
      const double dev = (res.T != 0.0) ? 100.0 * (t / res.T - 1.0) : 0.0;
      if (i < 3 && std::fabs(dev) > maxdev) maxdev = std::fabs(dev);
      dg << Form("    %s in [%5.1f,%5.1f) : B=%8.1f  D=%8.1f  T=%.4f  (%+.1f%% vs incl)\n",
                        YName(plane).c_str(), edges[i][0], edges[i][1], b.n, d.n, t, dev);
    }
    res.subSlice = maxdev;
    dg << Form("    -> max |deviation| over the thirds = %.1f%%\n", maxdev);
  }

  // ------------------------------------------------------------------------
  // SYSTEMATIC 1: how much does the answer move with the anti-iso window?
  // Both edges are varied. NB the TEMPLATE TOTAL is far more stable than
  // Apred, because the total is anchored by the directly measured N_B and only
  // its high-y tail is extrapolated (total = N_B + Apred by construction).
  // ------------------------------------------------------------------------
  {
    dg << "  --- anti-iso window variation (nominal ["
              << Form("%.3f,%.2f", cfg.isoFailLo, cfg.isoFailHi) << ")) ---\n";
    const double lo[4] = {cfg.isoFailLo - 0.05, cfg.isoFailLo, cfg.isoFailLo + 0.05,
                          cfg.isoFailLo + 0.10};
    const double hi[3] = {cfg.isoFailHi, 0.60, 0.50};
    double maxdev = 0.0, maxdev40 = 0.0;
    for (int il = 0; il < 4; ++il)
      for (int ih = 0; ih < 3; ++ih)
      {
        if (il != 1 && ih != 0) continue;          // vary one edge at a time
        if (hi[ih] <= lo[il] + 0.05) continue;     // keep a sane window
        const Count c = count2D(qcd2D, lo[il], hi[ih], cfg.yCut, 1e9);
        const Count d = count2D(qcd2D, lo[il], hi[ih], 0.0, cfg.yCut);
        const Count all = count2D(qcd2D, lo[il], hi[ih], 0.0, 1e9);
        if (d.n == 0.0) continue;
        const double t = res.B.n / d.n, ap = t * c.n, tot = t * all.n;
        const double dev = (res.Apred != 0.0) ? 100.0 * (ap / res.Apred - 1.0) : 0.0;
        const double devT = (res.qcdTotPass != 0.0) ? 100.0 * (tot / res.qcdTotPass - 1.0) : 0.0;
        if (std::fabs(dev) > maxdev) maxdev = std::fabs(dev);
        // the in-fit observable A40 = B*C40/D (the abcd-mode normalization):
        // its window dependence is the systematic that survives in the
        // reduced kappa, so it is scanned alongside Apred.
        const Count c40 = count2D(qcd2D, lo[il], hi[ih], kSRMtCut, 1e9);
        const double ap40 = t * c40.n;
        const double dev40 = (res.Apred40 != 0.0) ? 100.0 * (ap40 / res.Apred40 - 1.0) : 0.0;
        if (std::fabs(dev40) > maxdev40) maxdev40 = std::fabs(dev40);
        dg << Form("    relIso [%.3f,%.2f): T=%.4f  Apred=%7.1f (%+6.1f%%)"
                          "  template total=%7.1f (%+5.1f%%)",
                          lo[il], hi[ih], t, ap, dev, tot, devT);
        if (plane == "mt")
          dg << Form("  A40=%7.1f (%+6.1f%%)", ap40, dev40);
        dg << "\n";
      }
    res.winSpread = maxdev;
    res.winSpread40 = maxdev40;
    dg << Form("    -> max |deviation| of Apred = %.1f%%\n", maxdev);
    if (plane == "mt")
      dg << Form("    -> max |deviation| of the in-fit A40 = B*C(m_T>%.0f)/D = %.1f%%\n",
                 kSRMtCut, maxdev40);
  }

  res.diag = dg.str();

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
  if (w2D) delete w2D;
  if (z2D) delete z2D;
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
                  TemplateNorm &norm, const std::string &tag = "",
                  double *tiltOut = nullptr)  // [tilt %, <pT> lower half, <pT> upper half]
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
  norm = ScaleToIsoPass(tmpl, metRes.T, metRes.Terr);

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
    // SYSTEMATIC 2: quantify the shape tilt the plot shows. relIso carries pT in
    // its denominator, so a deeper anti-iso slice preferentially selects softer
    // leptons; if the two halves disagree, the borrowed pT shape is biased.
    const double mLo = sLo->GetMean(), mHi = sHi->GetMean();
    const double tilt = (mLo != 0.0) ? 100.0 * (mHi / mLo - 1.0) : 0.0;
    if (tiltOut) { tiltOut[0] = tilt; tiltOut[1] = mLo; tiltOut[2] = mHi; }
    // (printed by the caller, under its own header)

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


// ----------------------------------------------------------------------------
// Everything the AN quotes that is NOT a single-region number: the plane-to-
// plane transport, its sensitivity to the assumed signal strength, the
// composition of the fitted selections, and the assembled systematic budget.
// Collected per charge so the report can be printed once at the end.
// ----------------------------------------------------------------------------
struct ChargeInfo
{
  ABCDResult met, mt;
  double ptTot = 0, ptTotErr = 0;      // lepton-pT template, no m_T cut
  double ptTot40 = 0, ptTot40Err = 0;  // lepton-pT template, m_T > 40
  double tiltNom = 0, tilt40 = 0;      // anti-iso pT-shape tilt, %
  double dataAll = 0, data40 = 0;      // iso-pass data yields of the two selections
  double ffShift = 0;                  // FF-weighted vs flat-T total, % (runFFCheck)
};

// Iso-pass data yield of a lepton-pT plane (the denominator of the QCD fraction).
double isoPassData(TFile *fData, const std::string &hname, double isoCut)
{
  TH2D *h = (TH2D *)fData->Get(hname.c_str());
  if (!h) return 0.0;
  return count2D(h, 0.0, isoCut, 0.0, 1e9).n;
}

// ----------------------------------------------------------------------------
// Per-pT-bin fake-factor DIAGNOSTIC (2026-08-23). The in-fit ABCD (like the
// old flat-T template) assumes the QCD isolation transfer factor is
// pT-independent: the SR shape is the anti-iso pT distribution times ONE
// number. The group's alternative ("fake factor method") would instead measure
// F(pT) = QCD(iso-pass)/QCD(anti-iso) per pT bin in the dijet-dominated
// m_T < yCut region and apply it bin-by-bin to the anti-iso m_T > 40 sideband.
// Here that is run as a CHECK of the flat-T assumption, deliberately kept out
// of the fit (the per-bin prompt subtraction is still prefit-r and the B-region
// stats are thin; if F(pT) shows a genuine trend, the documented upgrade is
// Combine's RooParametricHist bin-by-bin ABCD): compare
//   P(pT)    = F(pT) x N_antiiso,mT>40(pT)     (fake-factor prediction)
//   Flat(pT) = T_mT  x N_antiiso,mT>40(pT)     (what the fit uses)
// and quote the per-bin F/T deviations and the total shift.
// Inputs: the (pT x m_T) scan planes h_pt_mt[_antiiso]_* (pT 1 GeV, m_T 2.5
// GeV, filled for pT > kPtScanFloor = 20; anti-iso window FROZEN AT SKIM TIME
// to the nominal one). No (pT x MET) plane exists, so the MET-sideband
// counterpart of F(pT) is only available pT-integrated (= T_MET). Everything
// below pT 25 (outside the W selection) is dropped.
void runFFCheck(TFile *fData, const std::vector<MCFile> &mc, const std::string &chg,
                const ABCDConfig &cfg, const std::string &outDir,
                const std::string &chgLatex, const std::string &lep,
                ChargeInfo &ci, TFile *fout) // non-const: stores ffShift
{
  const std::string hpass = Form("h_pt_mt_%s%s", lep.c_str(), chg.c_str());
  const std::string hanti = Form("h_pt_mt_antiiso_%s%s", lep.c_str(), chg.c_str());
  TH2D *dPass = (TH2D *)fData->Get(hpass.c_str());
  TH2D *dAnti = (TH2D *)fData->Get(hanti.c_str());
  if (!dPass || !dAnti)
  {
    std::cerr << "[WARN] " << hpass << " / " << hanti << " not in the skim output"
              << " (pre-2026-07-30 skim?) -> fake-factor diagnostic skipped\n";
    return;
  }
  dPass = (TH2D *)dPass->Clone(Form("%s_ffdata", hpass.c_str()));
  dPass->SetDirectory(nullptr);
  dAnti = (TH2D *)dAnti->Clone(Form("%s_ffdata", hanti.c_str()));
  dAnti->SetDirectory(nullptr);
  TH2D *ePass = sumEWK2D(mc, hpass.c_str(), Form("%s_ffewk", hpass.c_str()));
  TH2D *eAnti = sumEWK2D(mc, hanti.c_str(), Form("%s_ffewk", hanti.c_str()));
  if (ePass) dPass->Add(ePass, -1.0); // QCD-only (prefit r=1 subtraction)
  if (eAnti) dAnti->Add(eAnti, -1.0);

  const double kPtMin = 25.0; // W-selection threshold (scan planes reach 20)
  const int kRebin = 5;       // 1 GeV -> 5 GeV pT bins (B stats are thin)
  const double eps = 1e-6;

  // ProjectionX (pT) over an m_T window; hi >= axis max includes the overflow.
  auto projPt = [&](TH2D *h, double ylo, double yhi, const char *nm) -> TH1D * {
    const TAxis *ay = h->GetYaxis();
    const int bylo = (ylo <= ay->GetXmin()) ? 1 : ay->FindBin(ylo + eps);
    const int byhi = (yhi >= ay->GetXmax()) ? ay->GetNbins() + 1 : ay->FindBin(yhi - eps);
    TH1D *p = h->ProjectionX(nm, bylo, byhi);
    p->SetDirectory(nullptr);
    p->Rebin(kRebin);
    // drop the [20, 25) scan-floor slice: outside the fitted selection
    for (int i = 1; i <= p->GetNbinsX(); ++i)
      if (p->GetXaxis()->GetBinUpEdge(i) <= kPtMin + eps)
      {
        p->SetBinContent(i, 0.0);
        p->SetBinError(i, 0.0);
      }
    return p;
  };

  const std::string tagc = lep + chg;
  TH1D *numLow = projPt(dPass, 0.0, cfg.yCut, Form("ff_num_%s", tagc.c_str()));
  TH1D *denLow = projPt(dAnti, 0.0, cfg.yCut, Form("ff_den_%s", tagc.c_str()));
  TH1D *anti40 = projPt(dAnti, kSRMtCut, 1e9, Form("ff_anti40_%s", tagc.c_str()));

  TH1D *ff = (TH1D *)numLow->Clone(Form("ff_pt_%s", tagc.c_str()));
  ff->SetDirectory(nullptr);
  ff->SetTitle(Form("F(p_{T}) = QCD(iso-pass)/QCD(anti-iso) at m_{T}<%.0f, %s;"
                    "p_{T}^{%s} [GeV];F(p_{T})",
                    cfg.yCut, chgLatex.c_str(), lep == "ele" ? "e" : "#mu"));
  ff->Divide(denLow); // independent samples -> default (uncorrelated) errors

  // Bins where F is not measurable (EWK-subtracted denominator below kMinDen
  // events) fall back to the flat T in the prediction: a near-zero denominator
  // otherwise turns one empty tail bin into a wild multiplier that dominates
  // the total. Those ff bins are zeroed (clean plot), the table marks them,
  // and the fallback fraction of the anti-iso m_T>40 yield is reported.
  const double kMinDen = 3.0;
  const int nb = ff->GetNbinsX();
  TH1D *pred = (TH1D *)anti40->Clone(Form("ff_pred_%s", tagc.c_str()));
  pred->SetDirectory(nullptr);
  TH1D *flat = (TH1D *)anti40->Clone(Form("ff_flat_%s", tagc.c_str()));
  flat->SetDirectory(nullptr);
  flat->Scale(ci.mt.T);
  double fbYield = 0.0;
  for (int i = 1; i <= nb + 1; ++i) // include the pT overflow
  {
    const bool ok = denLow->GetBinContent(i) >= kMinDen;
    const double f  = ok ? ff->GetBinContent(i) : ci.mt.T;
    const double fe = ok ? ff->GetBinError(i)   : ci.mt.Terr;
    const double a  = anti40->GetBinContent(i), ae = anti40->GetBinError(i);
    pred->SetBinContent(i, f * a);
    pred->SetBinError(i, std::sqrt(f * f * ae * ae + a * a * fe * fe));
    if (!ok)
    {
      fbYield += a;
      ff->SetBinContent(i, 0.0);
      ff->SetBinError(i, 0.0);
    }
  }

  // ---- printed table + totals (the log is the AN's source) ----
  std::cout << "\n=== per-pT-bin fake factor (diagnostic)  " << lep << " " << chg
            << "  [m_T plane] ===\n"
            << Form("  F(pT) from m_T < %.0f (QCD-only, prefit EWK subtraction),"
                    " applied to anti-iso m_T > %.0f;  flat reference T_mT = %.4f\n",
                    cfg.yCut, kSRMtCut, ci.mt.T);
  std::cout << "     pT bin       N(iso,mTlow)  D(anti,mTlow)   F(pT)              F/T_mT-1\n";
  for (int i = 1; i <= nb; ++i)
  {
    if (ff->GetXaxis()->GetBinUpEdge(i) <= kPtMin + eps) continue;
    if (denLow->GetBinContent(i) <= 0.0 && numLow->GetBinContent(i) <= 0.0) continue;
    if (denLow->GetBinContent(i) >= kMinDen)
      std::cout << Form("   [%3.0f,%3.0f)    %10.1f    %10.1f     %.4f +/- %.4f   %+7.1f%%\n",
                        ff->GetXaxis()->GetBinLowEdge(i), ff->GetXaxis()->GetBinUpEdge(i),
                        numLow->GetBinContent(i), denLow->GetBinContent(i),
                        ff->GetBinContent(i), ff->GetBinError(i),
                        ci.mt.T != 0.0 ? 100.0 * (ff->GetBinContent(i) / ci.mt.T - 1.0) : 0.0);
    else
      std::cout << Form("   [%3.0f,%3.0f)    %10.1f    %10.1f     -- (D < %.0f: flat-T fallback)\n",
                        ff->GetXaxis()->GetBinLowEdge(i), ff->GetXaxis()->GetBinUpEdge(i),
                        numLow->GetBinContent(i), denLow->GetBinContent(i), kMinDen);
  }
  double eP = 0, eF = 0, eA = 0;
  const double sP = pred->IntegralAndError(1, nb + 1, eP); // incl pT>100 overflow
  const double sF = flat->IntegralAndError(1, nb + 1, eF);
  const double sA = anti40->IntegralAndError(1, nb + 1, eA);
  ci.ffShift = (sF != 0.0) ? 100.0 * (sP / sF - 1.0) : 0.0;
  std::cout << Form("  totals (pT >= %.0f, incl overflow):  FF-weighted = %.1f +/- %.1f"
                    "   flat-T = %.1f +/- %.1f   ->  FF/flat - 1 = %+.1f%%\n",
                    kPtMin, sP, eP, sF, eF, ci.ffShift);
  std::cout << Form("  (%.1f%% of the anti-iso m_T>%.0f yield used the flat-T fallback:"
                    " F unmeasurable there)\n",
                    sA != 0.0 ? 100.0 * fbYield / sA : 0.0, kSRMtCut);
  std::cout << Form("  consistency: anti-iso(m_T>%.0f) count here = %.1f +/- %.1f  vs the"
                    " iso x m_T plane C40 = %.1f  (same events, two binnings)\n",
                    kSRMtCut, sA, eA, ci.mt.C40.n);
  std::cout << Form("  NB no (pT x MET) skim plane exists -> the MET-sideband F(pT) is only"
                    " available pT-integrated: T_MET = %.4f (vs T_mT = %.4f)\n",
                    ci.met.T, ci.mt.T);

  // ---- F(pT) with the flat-T line ----
  {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TCanvas *c = new TCanvas(Form("cff_%s", tagc.c_str()), "", 800, 650);
    c->SetLeftMargin(0.13);
    ff->GetXaxis()->SetRangeUser(kPtMin, 100.0);
    ff->SetMinimum(0.0);
    ff->SetMaximum(std::max(3.0 * ci.mt.T, 1.3 * ff->GetMaximum()));
    ff->SetLineColor(kBlue + 1);
    ff->SetMarkerColor(kBlue + 1);
    ff->SetMarkerStyle(20);
    ff->Draw("E1");
    TLine *l = new TLine(kPtMin, ci.mt.T, 100.0, ci.mt.T);
    l->SetLineColor(kRed + 1);
    l->SetLineWidth(2);
    l->SetLineStyle(2);
    l->Draw();
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(0.035);
    t.DrawLatex(0.16, 0.86, Form("F(p_{T}) at m_{T} < %.0f GeV, %s", cfg.yCut, chgLatex.c_str()));
    t.DrawLatex(0.16, 0.80, Form("#color[633]{dashed: flat T_{m_{T}} = %.4f}", ci.mt.T));
    c->SaveAs((outDir + Form("/ff_pt_%s.png", tagc.c_str())).c_str());
    c->SaveAs((outDir + Form("/ff_pt_%s.pdf", tagc.c_str())).c_str());
    delete l;
    delete c;
  }

  // ---- FF-weighted vs flat-T prediction, with ratio pad ----
  {
    PlotStyle ps;
    ps.logy = true;
    ps.yTitleOffset = 1.55;
    ps.headerX = 0.60;
    ps.headerY = 0.66;
    ps.titleSize = 0.04;
    // Both inputs are data-derived (the helper's data/MC styles are cosmetics).
    SaveDataMCRatio(pred, flat,
                    outDir + Form("/ff_closure_pt_%s", tagc.c_str()),
                    Form("p_{T}^{%s} [GeV]", lep == "ele" ? "e" : "#mu"),
                    "QCD events",
                    Form("QCD prediction (anti-iso, m_{T}>%.0f), %s", kSRMtCut, chgLatex.c_str()),
                    "fake-factor vs flat-T", "", ps, /*normToData=*/false,
                    "F(p_{T})-weighted", Form("flat T_{m_{T}} = %.4f", ci.mt.T),
                    "FF / flat");
  }

  fout->cd();
  ff->Write(ff->GetName(), TObject::kOverwrite);

  delete numLow;
  delete denLow;
  delete anti40;
  delete ff;
  delete pred;
  delete flat;
  delete dPass;
  delete dAnti;
  if (ePass) delete ePass;
  if (eAnti) delete eAnti;
}

void printChannelReport(const std::string &lep, const ABCDConfig &cfg,
                        const ChargeInfo &P, const ChargeInfo &M)
{
  const char *cn[2] = {"+", "-"};
  const ChargeInfo *ci[2] = {&P, &M};

  std::cout << "\n"
            << "================================================================\n"
            << "  CHANNEL REPORT  (" << lep << ")\n"
            << "================================================================\n";

  // ---- 1. transfer factors on the two planes -----------------------------
  // T is meant to be a property of the multijet LEPTON alone, so measuring it
  // against two different recoil variables tests the factorisation directly.
  std::cout << "\n-- transfer factor, both planes --\n";
  std::cout << "   chg   T(MET)             T(m_T)             T(m_T)/T(MET)-1\n";
  for (int i = 0; i < 2; ++i)
    std::cout << Form("   %-4s  %.4f +/- %.4f   %.4f +/- %.4f   %+6.1f%%\n", cn[i],
                      ci[i]->met.T, ci[i]->met.Terr, ci[i]->mt.T, ci[i]->mt.Terr,
                      ci[i]->met.T != 0 ? 100.0 * (ci[i]->mt.T / ci[i]->met.T - 1.0) : 0.0);

  // ---- 2. r-scan: how much of that spread is just the prefit EWK scale? ---
  // The EWK subtraction is done at the PREFIT normalisation (r = 1). Region B
  // on the MET plane is EWK-rich, so if the true signal strength is not 1 the
  // subtraction is wrong there and T_MET moves a lot; the m_T-plane B is much
  // purer and barely moves. Where the two planes cross IS the r the sideband
  // prefers -- compare it with the fitted r before assigning a "transport"
  // systematic, or the signal strength gets counted twice.
  std::cout << "\n-- sensitivity of T to the assumed W normalisation r (EWK subtracted at r) --\n";
  for (int i = 0; i < 2; ++i)
  {
    const ABCDResult &me = ci[i]->met, &mt = ci[i]->mt;
    std::cout << Form("   %s%s :  W fraction of region B:  MET %.1f%% , m_T %.1f%%\n", lep.c_str(),
                      cn[i], me.dataB.n > 0 ? 100.0 * me.wB.n / me.dataB.n : 0.0,
                      mt.dataB.n > 0 ? 100.0 * mt.wB.n / mt.dataB.n : 0.0);
    std::cout << "        r      T(MET)   T(m_T)   diff\n";
    double cross = 0.0, prev = 0.0;
    for (int k = 0; k <= 12; ++k)
    {
      const double r = 1.0 + 0.025 * k;
      const double tm = (me.dataD.n - r * me.wD.n - me.zD.n) != 0.0
                            ? (me.dataB.n - r * me.wB.n - me.zB.n) /
                                  (me.dataD.n - r * me.wD.n - me.zD.n) : 0.0;
      const double tt = (mt.dataD.n - r * mt.wD.n - mt.zD.n) != 0.0
                            ? (mt.dataB.n - r * mt.wB.n - mt.zB.n) /
                                  (mt.dataD.n - r * mt.wD.n - mt.zD.n) : 0.0;
      const double d = (tm != 0.0) ? 100.0 * (tt / tm - 1.0) : 0.0;
      if (k > 0 && prev < 0.0 && d >= 0.0) cross = r - 0.025 * d / (d - prev);
      prev = d;
      if (k % 2 == 0)
        std::cout << Form("      %.3f   %.4f   %.4f   %+6.1f%%\n", r, tm, tt, d);
    }
    if (cross > 0.0)
      std::cout << Form("      -> the two planes agree at r = %.2f  (compare with the FITTED r:"
                        " if they match, the spread is the prefit normalisation, NOT"
                        " iso-recoil correlation)\n", cross);
    else
      std::cout << "      -> no crossing in r = [1.0, 1.3]: the spread is r-stable and more"
                   " likely a genuine correlation\n";
  }

  // ---- 3. composition of the two fitted selections -----------------------
  // Per charge AND charge-summed: the fit is per charge, but the AN quotes the
  // combined fraction when it describes how multijet-dominated a channel is.
  std::cout << "\n-- multijet fraction of the fitted selections --\n";
  std::cout << "   selection                          chg   data      QCD        frac\n";
  for (int i = 0; i < 2; ++i)
    std::cout << Form("   W sel, no MET/m_T cut (MET fit)    %-4s  %7.0f  %8.1f   %5.1f%%\n", cn[i],
                      ci[i]->dataAll, ci[i]->ptTot,
                      ci[i]->dataAll > 0 ? 100.0 * ci[i]->ptTot / ci[i]->dataAll : 0.0);
  for (int i = 0; i < 2; ++i)
    std::cout << Form("   W sel + m_T>40   (lepton-pT fit)   %-4s  %7.0f  %8.1f   %5.1f%%\n", cn[i],
                      ci[i]->data40, ci[i]->ptTot40,
                      ci[i]->data40 > 0 ? 100.0 * ci[i]->ptTot40 / ci[i]->data40 : 0.0);
  for (int i = 0; i < 2; ++i)
    std::cout << Form("   high-MET signal region A          %-4s  %7.0f  %8.1f   %5.1f%%\n", cn[i],
                      ci[i]->met.dataA.n, ci[i]->met.Apred,
                      ci[i]->met.dataA.n > 0 ? 100.0 * ci[i]->met.Apred / ci[i]->met.dataA.n : 0.0);
  auto sumrow = [](const char *tag, double d, double q) {
    std::cout << Form("   %-33s both  %7.0f  %8.1f   %5.1f%%\n", tag, d, q,
                      d > 0 ? 100.0 * q / d : 0.0);
  };
  std::cout << "   -- charge-summed --\n";
  sumrow("W sel, no MET/m_T cut (MET fit)", P.dataAll + M.dataAll, P.ptTot + M.ptTot);
  sumrow("W sel + m_T>40   (lepton-pT fit)", P.data40 + M.data40, P.ptTot40 + M.ptTot40);
  sumrow("high-MET signal region A", P.met.dataA.n + M.met.dataA.n, P.met.Apred + M.met.Apred);

  // ---- 4. the assembled systematic budget --------------------------------
  // Quadrature sum of the components measured above -> the log-normal width for
  // the qcd_rate_<flav>_<chg> nuisance in the datacards. Quoted per charge and
  // then symmetrised (the cards use one kappa per flavour).
  std::cout << "\n-- systematic budget on the multijet normalisation (-> lnN kappa) --\n";
  std::cout << "   component                        " << lep << "+       " << lep << "-\n";
  const double statT[2]  = {P.met.T != 0 ? 100 * P.met.Terr / P.met.T : 0,
                            M.met.T != 0 ? 100 * M.met.Terr / M.met.T : 0};
  const double statSh[2] = {P.ptTot40 != 0 ? 100 * std::sqrt(std::max(0.0,
                              P.ptTot40Err * P.ptTot40Err -
                              std::pow(P.ptTot40 * statT[0] / 100.0, 2))) / P.ptTot40 : 0,
                            M.ptTot40 != 0 ? 100 * std::sqrt(std::max(0.0,
                              M.ptTot40Err * M.ptTot40Err -
                              std::pow(M.ptTot40 * statT[1] / 100.0, 2))) / M.ptTot40 : 0};
  const double transp[2] = {P.met.T != 0 ? 100 * std::fabs(P.mt.T / P.met.T - 1.0) : 0,
                            M.met.T != 0 ? 100 * std::fabs(M.mt.T / M.met.T - 1.0) : 0};
  const double window[2] = {P.met.winSpread, M.met.winSpread};
  const double tilt[2]   = {std::fabs(P.tilt40), std::fabs(M.tilt40)};
  auto row = [&](const char *n, const double v[2]) {
    std::cout << Form("   %-30s %6.1f%%  %6.1f%%\n", n, v[0], v[1]);
  };
  row("stat on T", statT);
  row("sideband stat (m_T>40 tmpl)", statSh);
  row("plane transport (m_T vs MET)", transp);
  row("anti-iso window", window);
  row("anti-iso pT-shape tilt", tilt);
  double tot[2];
  for (int i = 0; i < 2; ++i)
    tot[i] = std::sqrt(statT[i] * statT[i] + statSh[i] * statSh[i] + transp[i] * transp[i] +
                       window[i] * window[i] + tilt[i] * tilt[i]);
  row("TOTAL (quadrature)", tot);
  const double kap = 1.0 + 0.01 * std::max(tot[0], tot[1]);
  std::cout << Form("   => lnN kappa (max of the two charges)            = %.2f\n", kap);
  // The same budget without the transport row, for when the r-scan shows that
  // row is the prefit normalisation rather than a factorisation violation.
  double totNoT[2];
  for (int i = 0; i < 2; ++i)
    totNoT[i] = std::sqrt(statT[i] * statT[i] + statSh[i] * statSh[i] +
                          window[i] * window[i] + tilt[i] * tilt[i]);
  std::cout << Form("   => lnN kappa EXCLUDING the transport row         = %.2f\n",
                    1.0 + 0.01 * std::max(totNoT[0], totNoT[1]));
  std::cout << "   NB the transport row is only a factorisation systematic if the r-scan above\n"
               "      does NOT cross near the fitted r; if it does, that row is the prefit\n"
               "      normalisation and including it double-counts the signal strength.\n";

  // ---- 5. the in-fit ABCD (QCD_MODE=abcd) numbers ------------------------
  // The 2026-08-23 migration moves the leppt_mt40 QCD normalisation into the
  // likelihood: counting CR channels B (iso-pass, m_T<yCut), C40 (anti-iso,
  // m_T>40) and D (anti-iso, m_T<yCut) with free QCD scales, SR total =
  // (sB*sC/sD)*B0*C40/D0 via a functional rateParam, EWK subtraction riding
  // the fitted r's. Printed here: the prefit inputs the CR templates carry
  // (also exported as abcd_counts_* in the rootfile) and the REDUCED lnN that
  // stays on the SR qcd process -- the two stat rows are profiled in-fit and
  // the plane-transport row does not exist (the m_T plane IS the counting
  // plane; T_MET becomes a pure cross-check).
  std::cout << Form("\n-- in-fit ABCD (m_T plane, C re-cut at m_T>%.0f; QCD_MODE=abcd) --\n",
                    kSRMtCut);
  std::cout << "   chg   B0        C40       D0        A40=B0*C40/D0      [T_MET-norm total]\n";
  for (int i = 0; i < 2; ++i)
    std::cout << Form("   %-4s  %8.1f  %8.1f  %8.1f  %7.1f +/- %5.1f   [%7.1f]\n", cn[i],
                      ci[i]->mt.B.n, ci[i]->mt.C40.n, ci[i]->mt.D.n,
                      ci[i]->mt.Apred40, ci[i]->mt.Apred40err, ci[i]->ptTot40);
  std::cout << "   (A40 = the abcd-mode SR template total; the bracketed T_MET x sideband\n"
               "    value is what QCD_MODE=lnN|free keep fitting -- the old qcd_pt_mt40 total)\n";
  std::cout << "\n-- in-fit ABCD reduced systematic budget (lnN on the SR qcd process only) --\n";
  std::cout << "   component                        " << lep << "+       " << lep << "-\n";
  const double window40[2] = {P.mt.winSpread40, M.mt.winSpread40};
  row(Form("anti-iso window (B*C%.0f/D)", kSRMtCut), window40);
  row("anti-iso pT-shape tilt", tilt);
  double totIF[2];
  for (int i = 0; i < 2; ++i)
    totIF[i] = std::sqrt(window40[i] * window40[i] + tilt[i] * tilt[i]);
  row("TOTAL (quadrature)", totIF);
  std::cout << Form("   => in-fit ABCD reduced lnN kappa (max of the two charges)  = %.2f\n",
                    1.0 + 0.01 * std::max(totIF[0], totIF[1]));
  // The fake-factor diagnostic measures the same iso-pT correlation MORE
  // DIRECTLY than the <pT> tilt: its FF-weighted vs flat-T total shift is the
  // normalization impact of letting F vary with pT. Quote the budget with that
  // row replacing the tilt -- use the LARGER kappa in the cards when they
  // differ (2026-08-23: mu ~identical, ele FF shift ~13% > tilt 7-9%).
  const double ffs[2] = {std::fabs(P.ffShift), std::fabs(M.ffShift)};
  row("alt: FF total shift (per-pT F)", ffs);
  double totFF[2];
  for (int i = 0; i < 2; ++i)
    totFF[i] = std::sqrt(window40[i] * window40[i] + ffs[i] * ffs[i]);
  std::cout << Form("   => reduced kappa with the FF shift replacing the tilt row  = %.2f\n",
                    1.0 + 0.01 * std::max(totFF[0], totFF[1]));
  // If the CR-B W-related content is FROZEN at absolute MC instead of floated
  // with the r POIs (QCD_WCR=frozen in the fork), the prefit-r subtraction
  // residual ~ (wB/dataB) * |r-1| comes back; quote it at |r-1| = 0.2 (the
  // r-scan crossing) so that card variant has its own kappa.
  const double resid[2] = {P.mt.dataB.n > 0 ? 100.0 * 0.2 * P.mt.wB.n / P.mt.dataB.n : 0.0,
                           M.mt.dataB.n > 0 ? 100.0 * 0.2 * M.mt.wB.n / M.mt.dataB.n : 0.0};
  double totIFz[2];
  for (int i = 0; i < 2; ++i)
    totIFz[i] = std::sqrt(totIF[i] * totIF[i] + resid[i] * resid[i]);
  std::cout << Form("   [QCD_WCR=frozen only] + prefit-W residual (wB/dataB x 0.2) = %.1f%% / %.1f%%"
                    "  => kappa = %.2f\n",
                    resid[0], resid[1], 1.0 + 0.01 * std::max(totIFz[0], totIFz[1]));
  std::cout << "================================================================\n";
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
            << Form("%.1f +/- %.1f", r.qcdTotPass, r.qcdTotPassErr)
            << Form("   (negative-tail clamp removed %.2f%%)\n", r.clampedFrac);
  if (plane == "mt")
  {
    // The in-fit ABCD numbers (QCD_MODE=abcd): C re-cut at the SR's own
    // m_T > 40; B/D untouched at y < yCut. A40 is the prefit value the CR
    // scales sB*sC/sD multiply in the fit = the abcd-mode SR template total.
    std::cout << Form("  --- in-fit ABCD (C re-cut at m_T > %.0f = the leppt_mt40 SR) ---\n",
                      kSRMtCut);
    printRow(Form("C40 iso-fail, m_T > %.0f     ", kSRMtCut), r.C40);
    std::cout << "    A40 = B*C40/D (abcd-mode norm)   = "
              << Form("%.1f +/- %.1f", r.Apred40, r.Apred40err) << "\n";
    std::cout << Form("    CR-B EWK split:  W-rel = %.1f   DY = %.1f   DYtau = %.1f"
                      "   (raw data in B = %.0f)\n",
                      r.wB.n, r.zB_dy.n, r.zB_dytau.n, r.dataB.n);
  }
  std::cout << "  --- context (signal region A) ---\n";
  std::cout << "    data(A) = " << Form("%.1f", r.dataA.n)
            << " ,  EWK(A) = " << Form("%.1f", r.ewkA.n)
            << " ,  QCD(A)/data(A) = "
            << Form("%.1f%%", r.dataA.n > 0 ? 100.0 * r.Apred / r.dataA.n : 0.0) << "\n";
  // Purity of each region: how much of the ABCD input is actually electroweak.
  // D (and C) must be QCD-dominated for the method to work; B is the sensitive
  // one -- where it is EWK-rich, T is a difference of two comparable numbers and
  // inherits the absolute MC normalisation (see the r-scan diagnostic).
  std::cout << "  --- region composition (data | EWK | QCD=data-EWK | EWK/data) ---\n";
  auto prow = [](const char *tag, const Count &d, const Count &e) {
    std::cout << Form("    %-24s %9.1f | %8.1f | %9.1f | %5.1f%%\n", tag, d.n, e.n,
                      d.n - e.n, d.n > 0 ? 100.0 * e.n / d.n : 0.0);
  };
  prow(Form("A iso-pass %s-high", y.c_str()), r.dataA, r.ewkA);
  prow(Form("B iso-pass %s-low",  y.c_str()), r.dataB, r.ewkB);
  prow(Form("C iso-fail %s-high", y.c_str()), r.dataC, r.ewkC);
  prow(Form("D iso-fail %s-low",  y.c_str()), r.dataD, r.ewkD);
  {
    const double model = r.ewkA.n + r.Apred;
    std::cout << Form("    closure in A: data/(EWK+QCD_pred) = %.3f   -> data sit %+.1f%%"
                      " (%+.0f events) vs EWK+QCD\n",
                      model > 0 ? r.dataA.n / model : 0.0,
                      model > 0 ? 100.0 * (r.dataA.n / model - 1.0) : 0.0,
                      r.dataA.n - model);
    std::cout << Form("       NB (data-EWK) in A = %.1f exceeds the ABCD prediction %.1f by"
                      " %.0f events: that is the SIGNAL STRENGTH the fit measures\n"
                      "       (region A is EWK-dominated, so data-EWK there is not a QCD"
                      " measurement), not an ABCD failure.\n",
                      r.A.n, r.Apred, r.A.n - r.Apred);
  }
  std::cout << r.diag;   // buffered factorisation test + window variation
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
  // W-related and DY-related are kept apart: only the W part scales with the
  // signal strength, which the r-sensitivity diagnostic in the channel report
  // needs. `mc` (their union) is what the pT plane and the closure plots use.
  std::vector<MCFile> mc, mcW, mcZ;
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
    const bool isW = (s.suffix.rfind("Wp", 0) == 0 || s.suffix.rfind("Wm", 0) == 0);
    (isW ? mcW : mcZ).push_back({s.label, f});
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
  ChargeInfo infoPlus, infoMinus;   // everything the channel report needs

  for (const Plane &pl : planes)
  {
    for (const auto &cq : charges)
    {
      ABCDResult res;
      TH1D *tmpl = runPlaneCharge(fData, mcW, mcZ, pl.key, cq.first, cfg, outDir,
                                  cq.second, pl.ytit, lep, res);
      if (tmpl)
      {
        printResult(pl.key, lep, cq.first, res);
        ChargeInfo &ci = (cq.first == "Plus") ? infoPlus : infoMinus;
        if (pl.key == "met") { ((cq.first == "Plus") ? metResPlus : metResMinus) = res; ci.met = res; }
        else                 { ci.mt = res; }
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
      TemplateNorm ptn;
      double tilt3[3] = {0.0, 0.0, 0.0};   // tilt %, <pT> lower half, <pT> upper half
      TH1D *tmpl = runPtCharge(fData, mc, cq.first, cfg, outDir, cq.second, lep, mres, ptn,
                               tag, tilt3);
      const double tilt = tilt3[0];
      if (!tmpl) continue;
      {
        ChargeInfo &ci = (cq.first == "Plus") ? infoPlus : infoMinus;
        const std::string sfx = tag.empty() ? "" : ("_" + tag);
        const double dat = isoPassData(fData,
            Form("h_iso_pt%s_%s%s", sfx.c_str(), lep.c_str(), cq.first.c_str()), cfg.isoCut);
        if (tag.empty()) { ci.ptTot = ptn.total;   ci.ptTotErr = ptn.err;   ci.tiltNom = tilt; ci.dataAll = dat; }
        else             { ci.ptTot40 = ptn.total; ci.ptTot40Err = ptn.err; ci.tilt40 = tilt;  ci.data40 = dat; }
      }
      const double tt = ptn.total, te = ptn.err;
      std::cout << "\n=== ABCD QCD lepton-pT template  " << lep << " " << cq.first
                << (tag.empty() ? "  [no m_T cut]" : "  [m_T > 40]") << " ===\n"
                << "  (anti-iso pT shape x MET-plane T = "
                << Form("%.4f +/- %.4f", mres.T, mres.Terr) << ")\n"
                << "  Total iso-pass QCD (template norm) = "
                << Form("%.1f +/- %.1f", tt, te)
                << Form("\n  anti-iso pT-shape tilt: <pT> = %.2f (relIso [%.2f,%.2f)) vs"
                        " %.2f (relIso [%.2f,%.2f)) GeV  ->  %+.1f%%",
                        tilt3[1], cfg.isoFailLo, 0.5 * (cfg.isoFailLo + cfg.isoFailHi),
                        tilt3[2], 0.5 * (cfg.isoFailLo + cfg.isoFailHi), cfg.isoFailHi, tilt3[0]);
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

  // --- per-pT-bin fake-factor diagnostic (flat-T shape-assumption check) ---
  runFFCheck(fData, mc, "Plus", cfg, outDir, lepP, lep, infoPlus, fout);
  runFFCheck(fData, mc, "Minus", cfg, outDir, lepM, lep, infoMinus, fout);

  printChannelReport(lep, cfg, infoPlus, infoMinus);

  // ---- machine-readable in-fit ABCD counts (m_T plane, C at m_T>40) --------
  // Consumed by plotting/mtandmet.C to build the abcd-normalized SR template
  // (qcd_abcd) and the CR-channel 1-bin templates for QCD_MODE=abcd in the
  // fork. Labeled bins: content = count, error = its stat uncertainty. qcdX0
  // are the r=1 EWK-subtracted QCD counts (the CR scales' init anchor); the
  // raw data and the W/DY components are carried separately so the CR channels
  // can float the subtraction in-fit (W-part per-y via the SR fractions, DY
  // via r_Z). T_met/T_mt ride along for provenance / staleness checks.
  for (int ic = 0; ic < 2; ++ic)
  {
    const ChargeInfo &ci = (ic == 0) ? infoPlus : infoMinus;
    const char *chg = (ic == 0) ? "Plus" : "Minus";
    const ABCDResult &r = ci.mt;
    const double bCheck = r.dataB.n - r.wB.n - r.zB_dy.n - r.zB_dytau.n;
    if (std::fabs(bCheck - r.B.n) > 0.1)
      std::cerr << Form("[WARN] abcd_counts_%s%s: B0 self-check off (%.2f vs %.2f)"
                        " -- exported pieces do not reassemble the QCD count\n",
                        lep.c_str(), chg, bCheck, r.B.n);
    TH1D *hc = new TH1D(Form("abcd_counts_%s%s", lep.c_str(), chg),
                        Form("in-fit ABCD counts (%s%s): m_T plane, B/D at m_T<%.0f,"
                             " C40 at m_T>%.0f, anti-iso [%.3f,%.2f)",
                             lep.c_str(), chg, cfg.yCut, kSRMtCut,
                             cfg.isoFailLo, cfg.isoFailHi),
                        15, 0.0, 15.0);
    hc->SetDirectory(nullptr);
    int b = 0;
    auto setb = [&](const char *lab, double v, double e) {
      ++b;
      hc->GetXaxis()->SetBinLabel(b, lab);
      hc->SetBinContent(b, v);
      hc->SetBinError(b, e);
    };
    setb("dataB", r.dataB.n, r.dataB.e);
    setb("qcdB0", r.B.n, r.B.e);
    setb("wB", r.wB.n, r.wB.e);
    setb("zB", r.zB_dy.n, r.zB_dy.e);
    setb("ztauB", r.zB_dytau.n, r.zB_dytau.e);
    setb("dataC40", r.dataC40.n, r.dataC40.e);
    setb("qcdC0", r.C40.n, r.C40.e);
    setb("wC40", r.wC40.n, r.wC40.e);
    setb("zC40", r.zC40.n, r.zC40.e);
    setb("dataD", r.dataD.n, r.dataD.e);
    setb("qcdD0", r.D.n, r.D.e);
    setb("wD", r.wD.n, r.wD.e);
    setb("zD", r.zD.n, r.zD.e);
    setb("T_met", ci.met.T, ci.met.Terr);
    setb("T_mt", r.T, r.Terr);
    fout->cd();
    hc->Write(hc->GetName(), TObject::kOverwrite);
    delete hc;
  }

  fout->Close();
  delete fout;
  fData->Close();
  delete fData;
  for (MCFile &m : mc) { m.f->Close(); delete m.f; }

  std::cout << "\n[DONE] ABCD QCD (" << lep << ") -> plots: " << outDir
            << "/   templates: " << outFile << "\n";
}
