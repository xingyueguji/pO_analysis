// ============================================================================
// qcd_abcd_nomt.h -- THE NO-m_T LEPTON-pT FIT's multijet (scheme C) and its
// diagnostics. Included by correction/qcd_abcd.C (not standalone: it uses that
// file's MCFile, ABCDConfig, kSRMtCut and the fake-factor likelihood ffLnL /
// ffMaximize / ffInterval); moved out of qcd_abcd.C on 2026-10-21b.
//
// Scheme C (2026-10-21, user decision; design study correction/qcd_nomt_study.C,
// records docs/iterations/2026-10-21_noreco-nomt-study.md and _leppt-nomt-schemeC.md).
// The signal region is iso-pass with NO m_T cut, fitted as two pT channels per
// rapidity bin split at m_T = 30 (skim h_leppt_mt{lo,hi}_*), sharing r: the
// m_T < 30 channel IS region B, so each bin's transfer factor is measured inside
// the likelihood from its own B pT spectrum. The multijet of bin y, channel c:
//   Q_{y,c}(pT) = T_y x f(pT; s) x Chat_{g,c}(pT) x N_{y,c} / N_{g,c}
// with g = the rapidity group of y (abcd_ygroups.h), Chat_{g,c} = the group's
// QCD-only (data - EWK at r = 1) anti-iso pT spectrum in channel c (per bin
// clamped >= 0) = the SHAPE (per group: per-bin shapes are template noise),
// N_{y,c} = that spectrum's total in bin y alone, f = exp(s (pT - 40)/10).
// PER-BIN LEVEL since 2026-10-21b (user decision, after the check that T varies
// within a group: e+ y10 0.62 vs y11 1.14, y6 0.12 vs y9 0.55): T_y is the
// bin's own fake factor at m_T < 30 at the shared slope (profile likelihood,
// every bin its own level; electron: one slope per charge shared by all bins;
// muon: flat, s = 0). In the fit the level floats freely per bin
// (qcd_sT_<F>_<C>_y<i>, prefit 1 = this template) and the bin's 1-bin
// anti-iso anchors D_y (m_T < 30) / H_y (m_T >= 30) carry the template
// normalization statistics (qcd_sD/qcd_sH). The electron slope is a shape
// nuisance at s +- kNoMtFFPriorK x the measurement's 68% interval (a wider
// prior: the B data that measured it are IN the fit), normalized so each
// group's m_T < 30 total is unchanged. The residual lnN on the m_T >= 30
// multijet only (the part extrapolated across m_T) is TWO terms with
// different correlations:
//   window   = the anti-iso window variation ([0.3,0.6) / [0.6,1.0)) of the
//              predicted hi/lo ratio with B fixed -- one method choice, so ONE
//              kappa per charge, COHERENT over all bins (fork qcd_rate_<F>_<C>);
//   transfer = F(m_T 30-X) / F(m_T < 30), X = 40 (muon) / 50 (electron),
//              measured POOLED over all 12 bins with the W MC at r_W = 1.2,
//              kappa = 1 + sqrt(dev^2 + stat^2), the same kappa for every
//              group but one UNCORRELATED theta per group (fork
//              qcd_xfer_<F>_<C>_<g>) -- since 2026-10-22 (user decision; it
//              was the group's own 30-40 value, which for the muon was pure
//              statistics, +-0.2-0.55); the per-group values are printed as
//              diagnostics, with their chi2 against the pooled one.
// The group shape's TAIL (pT >= 50) is the pooled all-bins anti-iso shape,
// scaled to the group's own integral above 50 (2026-10-22, user decision):
// every group's tail agrees with the pooled one (p 0.21-1.00), and the group
// tails were template noise -- single anti-iso events x f(pT) -- spiky in the
// forward groups and under-filled (negative bins clipped) in the backward ones.
// Outputs (only ADDED to the rootfile):
//   nomt_qcd_<lep><chg>_<lab|fb>_y<i>_<lo|hi>   (+ _ffUp / _ffDown, electron)
//   nomt_counts_<lep><chg>_<lab|fb>_y<i>         (labeled: the bin's anchors, T_y, slope)
//   nomt_counts_<lep><chg>_<lab|fb>_<group>      (labeled: the group's sums, common T_g -- diagnostics)
//   nomt_kappa_<lep><chg>                        (labeled: kappaWin + per binning/group kappaXfer)
// Plots: <plotDir> = plots/qcd_abcd_<lep>/nomt_fit/ (see NoMtDiag below).
// Inputs: the skim cubes h_iso_pt_mt_<lep><chg>_y<i>[_FB] (2026-10-21).
// ============================================================================
#ifndef PO_QCD_ABCD_NOMT_H
#define PO_QCD_ABCD_NOMT_H

const double kNoMtSplitMt  = pOSkim::kNoMtSplit; // 30 = cfg.yCut
const double kNoMtFFPriorK = 2.0;                // the slope templates at s +- 2 x its 68% interval
const double kNoMtWinMid   = 0.6;                // the anti-iso sub-window boundary (window variation)
// The W scale of the EWK subtraction when the transfer is measured above m_T 30,
// where the iso-pass muon data are 25-34% W: at r_W = 1 every W excess of the data
// (the muon m_T >= 30 channel prefers r ~ 1.20 / 1.22) would be read as multijet
// and inflate the transfer -- the signal strength the fit itself determines,
// counted again in a prior (the same trap as the 2026-08-19 ABCD transport row).
// The electron is insensitive (W 4-18% at 30-40). r_W = 1 is printed too.
const double kNoMtXferR    = 1.2;
// The upper edge of the transfer window (lower edge = the split, 30): the muon
// stops at 40 (the iso-pass W share is 50-57% at 40-50, and the ratio there
// moves with r_W); the electron goes to 50 (W share still 14-37% at 40-50, twice
// the multijet events, nearer the m_T >= 30 multijet's <m_T> ~ 48 GeV).
const double kNoMtXferMtHiMu  = 40.0;
const double kNoMtXferMtHiEle = 50.0;
// Above this lepton pT the group shape takes the pooled (all-bins) tail.
const double kNoMtTailPt = 50.0;

// ---- the diagnostic plots (2026-10-21b, user request: "diagnostic plots like
// we have for the abcd, to verify the potential issues"). One family per
// question the model answers or assumes:
//   nomt_regions_<lep><chg>              the (iso, m_T) plane with the four regions
//   nomt_T_y_<lep><chg>_<bin>            T_y per bin vs the one-level-per-group fit
//   nomt_lofit_<lep><chg>_<bin>          m_T < 30 channel data/prefit, group level vs per-bin level,
//                                        + the m_T >= 30 channel's needed multijet factor (W at r = 1.2)
//   nomt_ff_<lep><chg>_<bin>             F(pT) per group vs the shared slope (+ per-group slopes)
//   nomt_xfer_<lep><chg>_<bin>           F(m_T slice) / F(m_T < 30) per group + the pooled window value (the transfer kappa)
//   nomt_tail_<lep><chg>_<bin>_<g>_<lo|hi> the group's anti-iso pT shape vs the pooled one, and the tail taken above 50
//   nomt_shape_<lep><chg>_lab_<g>_<lo|hi> per-bin anti-iso pT shape / the group shape
//   nomt_shape_p_<lep><chg>_<bin>        p-value of each bin's shape vs its group shape
//   nomt_mt_<lep><chg>_<bin>_<g>         iso-pass m_T spectrum per group: data vs EWK (r = 1) + multijet
//   nomt_sr_<lep><chg>_<bin>_<g>_<lo|hi> per-group pT stacks of the two SR channels (W at r = 1)
namespace NoMtDiag
{
struct Series
{
  TGraphAsymmErrors *g;
  std::string label;
  int color, marker;
  double size;
  std::string opt = "P";
};

void Draw(const std::string &out, std::vector<Series> ser, double xlo, double xhi, double ylo, double yhi, const char *xt,
          const char *yt, const std::vector<std::string> &text, const std::function<void()> &under = nullptr,
          double legX1 = 0.62, double legY1 = 0.70, double legX2 = 0.93, double legY2 = 0.90)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas c(("c_" + out).c_str(), "", 850, 650);
  c.SetLeftMargin(0.13);
  c.SetRightMargin(0.04);
  c.SetTopMargin(0.06);
  TH1F *fr = c.DrawFrame(xlo, ylo, xhi, yhi);
  fr->GetXaxis()->SetTitle(xt);
  fr->GetYaxis()->SetTitle(yt);
  fr->GetYaxis()->SetTitleOffset(1.25);
  if (under) under();
  TLegend leg(legX1, legY1, legX2, legY2);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.030);
  for (Series &s : ser)
  {
    s.g->SetLineColor(s.color);
    s.g->SetMarkerColor(s.color);
    s.g->SetMarkerStyle(s.marker);
    s.g->SetMarkerSize(s.size);
    s.g->SetLineWidth(2);
    s.g->Draw((s.opt + " SAME").c_str());
    if (!s.label.empty()) leg.AddEntry(s.g, s.label.c_str(), s.opt == "L" ? "l" : "lp");
  }
  leg.Draw();
  TLatex t;
  t.SetNDC();
  t.SetTextFont(42);
  t.SetTextSize(0.029);
  double y = 0.89;
  for (const std::string &l : text) { t.DrawLatex(0.16, y, l.c_str()); y -= 0.045; }
  c.SaveAs((out + ".png").c_str());
  c.SaveAs((out + ".pdf").c_str());
}
} // namespace NoMtDiag

void runNoMt(TFile *fData, const std::vector<MCFile> &mcW, const std::vector<MCFile> &mcZ, const ABCDConfig &cfg,
             const std::string &lep, TFile *fout, const std::string &plotDir)
{
  std::cout << "\n=========================================================================\n"
            << " THE NO-m_T LEPTON-pT FIT (scheme C): SR split at m_T = " << kNoMtSplitMt
            << ", multijet = T_y x f(pT) x group anti-iso shape x per-bin anti-iso count (" << lep << ")\n"
            << "=========================================================================\n";
  gSystem->mkdir(plotDir.c_str(), kTRUE);
  const int errLevel0 = gErrorIgnoreLevel; // keep the TCanvas::Print "Info" lines out of the log
  gErrorIgnoreLevel = std::max(errLevel0, (int)kWarning);
  const std::vector<double> FE = {24, 28, 32, 36, 40, 46, 52, 60, 70, 80, 100, 1e9}; // F(pT) bins
  const double eps = 1e-6;
  const char *bins[2] = {"lab", "fb"};
  const char *chgLab[2] = {"+", "-"};
  const char *lepLatex = (lep == "ele") ? "e" : "#mu";
  const int gColor[pOAbcdY::kNGroups] = {kRed + 1, kAzure + 2, kGreen + 2};
  const int gMarker[pOAbcdY::kNGroups] = {21, 33, 20};
  const double gSize[pOAbcdY::kNGroups] = {1.5, 2.0, 1.3};
  const int kBinColor[8] = {kBlack, kRed + 1, kAzure + 2, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kViolet + 1};
  const int kBinMarker[8] = {20, 21, 33, 34, 24, 25, 27, 28};
  double kWinMax = 1.0, kXferMax = 1.0;

  // a pT spectrum (50 x 2 GeV + overflow) over iso [ilo, ihi) x m_T [mlo, mhi) of a cube;
  // every window edge must be a cube bin edge (checked)
  auto edgeBin = [](const TAxis *a, double v, bool low) {
    for (int b = 1; b <= a->GetNbins(); ++b)
      if (std::fabs((low ? a->GetBinLowEdge(b) : a->GetBinUpEdge(b)) - v) < 1e-9) return b;
    return -1;
  };
  auto slice = [&](TH3D *h, double ilo, double ihi, double mlo, double mhi, bool &ok) -> TH1D * {
    const int x1 = edgeBin(h->GetXaxis(), ilo, true), x2 = edgeBin(h->GetXaxis(), ihi, false);
    const int z1 = edgeBin(h->GetZaxis(), mlo, true);
    const int z2 = (mhi >= h->GetZaxis()->GetXmax()) ? h->GetZaxis()->GetNbins() + 1 : edgeBin(h->GetZaxis(), mhi, false);
    if (x1 < 0 || x2 < 0 || z1 < 0 || z2 < 0)
    {
      std::cerr << "[ERROR] runNoMt: window iso [" << ilo << "," << ihi << ") x m_T [" << mlo << "," << mhi
                << ") not on the bin edges of " << h->GetName() << "\n";
      ok = false;
      return nullptr;
    }
    static int uid = 0;
    TH1D *p = h->ProjectionY(Form("nomt_sl_%d", ++uid), x1, x2, z1, z2, "e");
    p->SetDirectory(nullptr);
    return p;
  };
  auto sumOver = [&](const TH1D *h, double lo, double hi) {
    const int b1 = h->GetXaxis()->FindBin(lo + eps);
    const int b2 = (hi >= 1e8) ? h->GetNbinsX() + 1 : h->GetXaxis()->FindBin(hi - eps);
    return h->Integral(b1, b2);
  };
  auto xOf = [](const TH1D *h, int b) { return (std::fabs(h->GetBinLowEdge(b) - 24.0) < 1e-6) ? 25.5 : h->GetBinCenter(b); };
  auto yCenter = [](int ib, int iy) {
    const double *e = ib ? pOSkim::kYEdgesFB : pOSkim::kYEdges;
    return 0.5 * (e[iy] + e[iy + 1]);
  };
  auto yHalf = [](int ib, int iy) {
    const double *e = ib ? pOSkim::kYEdgesFB : pOSkim::kYEdges;
    return 0.5 * (e[iy + 1] - e[iy]);
  };

  // ---- the region map (inclusive plane, data) ----
  for (int ic = 0; ic < 2; ++ic)
  {
    const char *chg = ic ? "Minus" : "Plus";
    TH2D *h2 = (TH2D *)fData->Get(Form("h_iso_mt_%s%s", lep.c_str(), chg));
    if (!h2) continue;
    gStyle->SetOptStat(0);
    TCanvas c(Form("c_nomt_regions_%s%s", lep.c_str(), chg), "", 900, 700);
    c.SetLogz();
    c.SetRightMargin(0.14);
    c.SetLeftMargin(0.12);
    TH2D *hd = (TH2D *)h2->Clone(Form("nomt_regions_%s%s", lep.c_str(), chg));
    hd->SetDirectory(nullptr);
    hd->RebinX(4);
    hd->RebinY(2);
    hd->GetXaxis()->SetRangeUser(cfg.isoPassLo, cfg.isoFailHi);
    hd->GetYaxis()->SetRangeUser(0, 120);
    hd->SetTitle(""); // axis by axis: the electron iso title itself contains a ';'
    hd->GetXaxis()->SetTitle(cfg.isoTitle.c_str());
    hd->GetYaxis()->SetTitle("m_{T} (GeV)");
    hd->GetZaxis()->SetTitle("events");
    hd->Draw("COLZ");
    TBox gap(cfg.isoCut, 0, cfg.isoFailLo, 120);
    gap.SetFillColorAlpha(kGray + 1, 0.45);
    gap.SetFillStyle(3004);
    gap.Draw();
    TLine lm(cfg.isoPassLo, kNoMtSplitMt, cfg.isoFailHi, kNoMtSplitMt);
    lm.SetLineColor(kRed + 1); lm.SetLineWidth(3); lm.SetLineStyle(2); lm.Draw();
    for (double xv : {cfg.isoCut, cfg.isoFailLo, kNoMtWinMid})
    {
      TLine *l = new TLine(xv, 0, xv, 120);
      l->SetLineColor(xv == kNoMtWinMid ? kGray + 2 : kBlack);
      l->SetLineWidth(xv == kNoMtWinMid ? 1 : 3);
      l->SetLineStyle(xv == kNoMtWinMid ? 3 : 1);
      l->Draw();
    }
    TLatex t;
    t.SetTextFont(62);
    t.SetTextSize(0.040);
    const double xp = 0.5 * (cfg.isoPassLo + cfg.isoCut), xa = 0.5 * (cfg.isoFailLo + cfg.isoFailHi);
    t.SetTextAlign(22);
    t.DrawLatex(xp, 15, "lo (= B)");
    t.DrawLatex(xp, 75, "hi (= A)");
    t.DrawLatex(xa, 15, "D");
    t.DrawLatex(xa, 75, "H");
    t.SetTextFont(42);
    t.SetTextSize(0.028);
    t.DrawLatex(0.5 * (cfg.isoCut + cfg.isoFailLo), 105, "gap");
    t.SetNDC();
    t.SetTextAlign(11);
    t.DrawLatex(0.13, 0.95, Form("%s^{%s} data: SR = iso-pass, all m_{T}, split at %.0f GeV; T measured lo/D, applied to H",
                                 lepLatex, chgLab[ic], kNoMtSplitMt));
    c.SaveAs((plotDir + Form("/nomt_regions_%s%s.png", lep.c_str(), chg)).c_str());
    c.SaveAs((plotDir + Form("/nomt_regions_%s%s.pdf", lep.c_str(), chg)).c_str());
    delete hd;
  }

  for (int ic = 0; ic < 2; ++ic)
  {
    const char *chg = ic ? "Minus" : "Plus";
    double kWin = 1.0, winDev[3] = {0, 0, 0};
    double xferDev[2][pOAbcdY::kNGroups] = {}, xferErr[2][pOAbcdY::kNGroups] = {}, kXfer[2][pOAbcdY::kNGroups] = {};
    double xferDevP[2] = {}, xferErrP[2] = {}; // the pooled transfer measurement (what sets kXfer)
    for (int ib = 0; ib < 2; ++ib)
    {
      const char *suf = ib ? "_FB" : "";
      const std::string tagCB = Form("%s%s_%s", lep.c_str(), chg, bins[ib]);
      // ---- per-bin slices: data, W-related MC x k_s, other MC x k_s ----
      // pT spectra: the SR channels (P), the anti-iso anchors (A), the two
      // window halves (A1, A2), m_T 30-40 (Pmid, Amid) and 40-50 (Pm4, Am4) =
      // the transfer window, and the m_T slices 0-10 / 10-20 / 20-30 (Pm*, Am*)
      enum { kPlo = 0, kPhi, kAlo, kAhi, kA1lo, kA1hi, kA2lo, kA2hi, kPmid, kAmid,
             kPm0, kPm1, kPm2, kAm0, kAm1, kAm2, kPm4, kAm4, kNSl };
      TH1D *S[pOSkim::kNY][kNSl][3] = {}; // [y][slice][0 data, 1 W, 2 other]
      // m_T spectra (all pT): [y][0 iso-pass, 1 anti-iso QCD weighted by T_y f(pT) (filled later)][sample]
      TH1D *M[pOSkim::kNY][3] = {};      // iso-pass m_T: data, W, other
      TH2D *Aqp[pOSkim::kNY] = {};       // anti-iso QCD-only (pT x m_T), for the multijet m_T template
      bool ok = true;
      for (int iy = 0; iy < pOSkim::kNY && ok; ++iy)
      {
        const std::string nm = Form("h_iso_pt_mt_%s%s_y%d%s", lep.c_str(), chg, iy, suf);
        TH3D *cube[3] = {(TH3D *)fData->Get(nm.c_str()), nullptr, nullptr};
        if (!cube[0])
        {
          std::cerr << "[WARN] runNoMt: " << nm << " not in the data skim (skim before 2026-10-21?) -> no no-m_T templates\n";
          return;
        }
        for (int k = 1; k <= 2; ++k)
          for (const MCFile &m : (k == 1 ? mcW : mcZ))
          {
            TH3D *h = (TH3D *)m.f->Get(nm.c_str());
            if (!h) { std::cerr << "[ERROR] runNoMt: " << nm << " missing in " << m.label << "\n"; return; }
            if (!cube[k]) { cube[k] = (TH3D *)h->Clone(Form("%s_sum%d", nm.c_str(), k)); cube[k]->SetDirectory(nullptr); cube[k]->Reset("ICESM"); }
            cube[k]->Add(h, pONorm::MCScale(m.label));
          }
        const double big = 1e9;
        const double win[kNSl][4] = {{cfg.isoPassLo, cfg.isoCut, 0, kNoMtSplitMt},        {cfg.isoPassLo, cfg.isoCut, kNoMtSplitMt, big},
                                     {cfg.isoFailLo, cfg.isoFailHi, 0, kNoMtSplitMt},     {cfg.isoFailLo, cfg.isoFailHi, kNoMtSplitMt, big},
                                     {cfg.isoFailLo, kNoMtWinMid, 0, kNoMtSplitMt},       {cfg.isoFailLo, kNoMtWinMid, kNoMtSplitMt, big},
                                     {kNoMtWinMid, cfg.isoFailHi, 0, kNoMtSplitMt},       {kNoMtWinMid, cfg.isoFailHi, kNoMtSplitMt, big},
                                     {cfg.isoPassLo, cfg.isoCut, kNoMtSplitMt, kSRMtCut}, {cfg.isoFailLo, cfg.isoFailHi, kNoMtSplitMt, kSRMtCut},
                                     {cfg.isoPassLo, cfg.isoCut, 0, 10},                  {cfg.isoPassLo, cfg.isoCut, 10, 20},
                                     {cfg.isoPassLo, cfg.isoCut, 20, 30},                 {cfg.isoFailLo, cfg.isoFailHi, 0, 10},
                                     {cfg.isoFailLo, cfg.isoFailHi, 10, 20},              {cfg.isoFailLo, cfg.isoFailHi, 20, 30},
                                     {cfg.isoPassLo, cfg.isoCut, kSRMtCut, 50},           {cfg.isoFailLo, cfg.isoFailHi, kSRMtCut, 50}};
        for (int s = 0; s < kNSl; ++s)
          for (int k = 0; k < 3; ++k) S[iy][s][k] = slice(cube[k], win[s][0], win[s][1], win[s][2], win[s][3], ok);
        // m_T spectra over all pT: iso-pass per sample, and the anti-iso QCD-only (pT x m_T) plane
        if (ok)
        {
          const TAxis *ax = cube[0]->GetXaxis(), *ay = cube[0]->GetYaxis(), *az = cube[0]->GetZaxis();
          const int ip1 = edgeBin(ax, cfg.isoPassLo, true), ip2 = edgeBin(ax, cfg.isoCut, false);
          const int ia1 = edgeBin(ax, cfg.isoFailLo, true), ia2 = edgeBin(ax, cfg.isoFailHi, false);
          const int nz = az->GetNbins(), ny = ay->GetNbins();
          for (int k = 0; k < 3; ++k)
          {
            M[iy][k] = new TH1D(Form("nomt_m_%s_%d_%d", tagCB.c_str(), iy, k), "", nz, az->GetXbins()->GetArray());
            M[iy][k]->SetDirectory(nullptr);
            M[iy][k]->Sumw2();
          }
          Aqp[iy] = new TH2D(Form("nomt_aq_%s_%d", tagCB.c_str(), iy), "", ny, ay->GetXmin(), ay->GetXmax(), nz, az->GetXbins()->GetArray());
          Aqp[iy]->SetDirectory(nullptr);
          for (int jz = 1; jz <= nz + 1; ++jz)
          {
            const int tz = std::min(jz, nz); // overflow into the last m_T bin
            for (int jy = 1; jy <= ny + 1; ++jy)
            {
              const int ty = std::min(jy, ny); // pT overflow into the last bin
              double a = 0;
              for (int k = 0; k < 3; ++k)
              {
                double vp = 0, ep2 = 0, va = 0;
                for (int jx = ip1; jx <= ip2; ++jx) { vp += cube[k]->GetBinContent(jx, jy, jz); ep2 += std::pow(cube[k]->GetBinError(jx, jy, jz), 2); }
                for (int jx = ia1; jx <= ia2; ++jx) va += cube[k]->GetBinContent(jx, jy, jz);
                M[iy][k]->SetBinContent(tz, M[iy][k]->GetBinContent(tz) + vp);
                M[iy][k]->SetBinError(tz, std::sqrt(std::pow(M[iy][k]->GetBinError(tz), 2) + ep2));
                a += (k == 0 ? va : -va);
              }
              Aqp[iy]->SetBinContent(ty, tz, Aqp[iy]->GetBinContent(ty, tz) + a);
            }
          }
        }
        delete cube[1]; // the MC sums (cube[0] is the file's)
        delete cube[2];
      }
      if (!ok) return;
      // QCD-only spectrum (data - W - other), per bin clamped >= 0 when asked
      auto qcdOf = [&](int iy, int sl, bool clamp) {
        TH1D *q = (TH1D *)S[iy][sl][0]->Clone(Form("nomt_q_%s%s%d_%d_%d", chg, suf, iy, sl, clamp ? 1 : 0));
        q->SetDirectory(nullptr);
        q->Add(S[iy][sl][1], -1.0);
        q->Add(S[iy][sl][2], -1.0);
        if (clamp)
          for (int b = 0; b <= q->GetNbinsX() + 1; ++b)
            if (q->GetBinContent(b) < 0.0) q->SetBinContent(b, 0.0);
        return q;
      };

      // ---- the fake-factor inputs at m_T < 30, per bin in the F(pT) bins ----
      struct BinF { std::vector<double> x, np, nf, bp, bf; double F0 = 0, F0lo = 0, F0hi = 0; };
      BinF Y[pOSkim::kNY];
      for (int iy = 0; iy < pOSkim::kNY; ++iy)
        for (size_t k = 0; k + 1 < FE.size(); ++k)
        {
          double sx = 0, sw = 0;
          if (FE[k + 1] < 1e8)
            for (int b = S[iy][kAlo][0]->FindBin(FE[k] + eps); b <= S[iy][kAlo][0]->FindBin(FE[k + 1] - eps); ++b)
            { sx += xOf(S[iy][kAlo][0], b) * S[iy][kAlo][0]->GetBinContent(b); sw += S[iy][kAlo][0]->GetBinContent(b); }
          Y[iy].x.push_back(FE[k + 1] >= 1e8 ? 105.0 : (sw > 0 ? sx / sw : 0.5 * (FE[k] + FE[k + 1])));
          Y[iy].np.push_back(sumOver(S[iy][kPlo][0], FE[k], FE[k + 1]));
          Y[iy].nf.push_back(sumOver(S[iy][kAlo][0], FE[k], FE[k + 1]));
          Y[iy].bp.push_back(std::max(sumOver(S[iy][kPlo][1], FE[k], FE[k + 1]) + sumOver(S[iy][kPlo][2], FE[k], FE[k + 1]), 1e-6));
          Y[iy].bf.push_back(std::max(sumOver(S[iy][kAlo][1], FE[k], FE[k + 1]) + sumOver(S[iy][kAlo][2], FE[k], FE[k + 1]), 1e-6));
        }
      auto lnY = [&](int iy, double F0, double s) {
        double t = 0;
        for (size_t k = 0; k < Y[iy].x.size(); ++k)
          t += ffLnL(F0 * std::exp(s * (Y[iy].x[k] - 40.0) / 10.0), Y[iy].np[k], Y[iy].nf[k], Y[iy].bp[k], Y[iy].bf[k]);
        return t;
      };
      auto maxY = [&](int iy, double s, double &best) { return ffMaximize([&](double f) { return lnY(iy, f, s); }, 0.0, 5.0, 60, best); };
      // the profile over per-bin levels of a set of bins at slope s
      auto profBins = [&](const std::vector<int> &ys, double s, double *F0out) {
        double tot = 0;
        for (int iy : ys)
        {
          double best = 0;
          const double F0 = maxY(iy, s, best);
          if (F0out) F0out[iy] = F0;
          tot += best;
        }
        return tot;
      };
      // a common level over a set of bins at slope s
      auto common = [&](const std::vector<int> &ys, double s, double &best) {
        return ffMaximize([&](double f) { double t = 0; for (int iy : ys) t += lnY(iy, f, s); return t; }, 0.0, 5.0, 60, best);
      };
      // the 68% profile interval of a 1D function's maximum (2 DlnL = 1)
      auto interval = [&](const std::function<double(double)> &L, double xhat, double step, double lo0, double &lo, double &hi) {
        const double L0 = L(xhat);
        for (int side = -1; side <= 1; side += 2)
        {
          double in = xhat, out = xhat;
          for (int k = 1; k <= 4000; ++k) { out = xhat + side * step * k; if (out <= lo0 || 2.0 * (L0 - L(out)) > 1.0) break; in = out; }
          for (int it = 0; it < 40 && out > lo0; ++it) { const double m = 0.5 * (in + out); (2.0 * (L0 - L(m)) > 1.0 ? out : in) = m; }
          (side < 0 ? lo : hi) = std::max(lo0, 0.5 * (in + out));
        }
      };
      std::vector<int> allY;
      for (int iy = 0; iy < pOSkim::kNY; ++iy) allY.push_back(iy);

      // ---- the shared slope, every bin its own level: USED for the electron;
      //      for the muon (flat model, s = 0) fitted as a diagnostic only ----
      double sFit = 0.0, sFitLo = 0.0, sFitHi = 0.0, F0hat[pOSkim::kNY] = {};
      {
        double Lmax = 0;
        sFit = ffMaximize([&](double s) { return profBins(allY, s, nullptr); }, -1.5, 1.5, 40, Lmax);
        interval([&](double s) { return profBins(allY, s, nullptr); }, sFit, 0.01, -9.0, sFitLo, sFitHi);
      }
      const double sHat = cfg.useFF ? sFit : 0.0, sLo = cfg.useFF ? sFitLo : 0.0, sHi = cfg.useFF ? sFitHi : 0.0;
      // flat vs the shared slope (the muon model's own test)
      const double lrtFlat = std::max(0.0, 2.0 * (profBins(allY, sFit, nullptr) - profBins(allY, 0.0, nullptr)));
      profBins(allY, sHat, F0hat);
      for (int iy = 0; iy < pOSkim::kNY; ++iy)
      {
        Y[iy].F0 = F0hat[iy];
        interval([&](double f) { return lnY(iy, f, sHat); }, F0hat[iy], 0.002, 0.0, Y[iy].F0lo, Y[iy].F0hi);
      }
      const double sUp = sHat + kNoMtFFPriorK * (sHi - sHat), sDn = sHat - kNoMtFFPriorK * (sHat - sLo);
      std::cout << Form("\n[NOMT] %s %-3s: F_y = T_y exp(s (pT-40)/10) at m_T < %.0f, every bin its own level T_y; %s\n",
                        (lep + chg).c_str(), bins[ib], kNoMtSplitMt,
                        cfg.useFF ? Form("shared slope s = %+.3f [%+.3f,%+.3f] (68%%); templates at s %+.3f / %+.3f (x%.0f)",
                                         sHat, sLo, sHi, sUp, sDn, kNoMtFFPriorK)
                                  : "flat (muon: s = 0)");

      // ---- the diagnostics of the level and slope assumptions ----
      // (a) one level per group? LRT of a common T_g vs per-bin T_y, per group, at the shared slope
      double F0g[pOAbcdY::kNGroups] = {}, F0glo[pOAbcdY::kNGroups] = {}, F0ghi[pOAbcdY::kNGroups] = {}, lrtT[pOAbcdY::kNGroups] = {};
      for (int g = 0; g < pOAbcdY::kNGroups; ++g)
      {
        const std::vector<int> ys = pOAbcdY::BinsOf(g);
        double Lc = 0;
        F0g[g] = common(ys, sHat, Lc);
        interval([&](double f) { double t = 0; for (int iy : ys) t += lnY(iy, f, sHat); return t; }, F0g[g], 0.002, 0.0, F0glo[g], F0ghi[g]);
        lrtT[g] = std::max(0.0, 2.0 * (profBins(ys, sHat, nullptr) - Lc));
        std::cout << Form("[NOMT-T] %s %-3s %-3s: common T_g %.4f [%.4f,%.4f] vs per-bin T_y:", (lep + chg).c_str(), bins[ib],
                          pOAbcdY::kGroupName[g], F0g[g], F0glo[g], F0ghi[g]);
        for (int iy : ys) std::cout << Form(" y%d %.3f[%.3f,%.3f]", iy, Y[iy].F0, Y[iy].F0lo, Y[iy].F0hi);
        std::cout << Form("  | LRT one-level-per-group: %.1f / %d dof, p = %.2g\n", lrtT[g], (int)ys.size() - 1,
                          ys.size() > 1 ? TMath::Prob(lrtT[g], (int)ys.size() - 1) : 1.0);
      }
      // (b) one slope for all groups? free slope per group (levels per bin), LRT vs the shared slope
      double sG[pOAbcdY::kNGroups] = {}, sGlo[pOAbcdY::kNGroups] = {}, sGhi[pOAbcdY::kNGroups] = {}, LgSum = 0;
      for (int g = 0; g < pOAbcdY::kNGroups; ++g)
      {
        const std::vector<int> ys = pOAbcdY::BinsOf(g);
        double Lg = 0;
        sG[g] = ffMaximize([&](double s) { return profBins(ys, s, nullptr); }, -1.5, 1.5, 40, Lg);
        interval([&](double s) { return profBins(ys, s, nullptr); }, sG[g], 0.01, -9.0, sGlo[g], sGhi[g]);
        LgSum += Lg;
      }
      const double lrtS = std::max(0.0, 2.0 * (LgSum - profBins(allY, sFit, nullptr)));
      std::cout << Form("[NOMT-SLOPE] %s %-3s: slope per group fwd %+.3f [%+.3f,%+.3f], bwd %+.3f [%+.3f,%+.3f], cen %+.3f [%+.3f,%+.3f]"
                        " vs shared %+.3f [%+.3f,%+.3f] | LRT one slope: %.1f / 2 dof, p = %.2g | LRT flat vs shared: %.1f / 1, p = %.2g%s\n",
                        (lep + chg).c_str(), bins[ib], sG[0], sGlo[0], sGhi[0], sG[1], sGlo[1], sGhi[1], sG[2], sGlo[2], sGhi[2],
                        sFit, sFitLo, sFitHi, lrtS, TMath::Prob(lrtS, 2), lrtFlat, TMath::Prob(lrtFlat, 1),
                        cfg.useFF ? "" : "  (the template uses flat, s = 0)");

      // ---- templates per bin and channel ----
      // group shapes Chat_{g,c} (clamped) and their f-weighted totals; per-bin counts N_{y,c}
      TH1D *Cg[pOAbcdY::kNGroups][2] = {};
      double Ng[pOAbcdY::kNGroups][2] = {}, Ny[pOSkim::kNY][2] = {};
      for (int iy = 0; iy < pOSkim::kNY; ++iy)
        for (int c = 0; c < 2; ++c)
        {
          TH1D *q = qcdOf(iy, c ? kAhi : kAlo, true);
          Ny[iy][c] = q->Integral(1, q->GetNbinsX());
          const int g = pOAbcdY::GroupOf(iy);
          if (!Cg[g][c]) { Cg[g][c] = (TH1D *)q->Clone(Form("nomt_Cg_%s%s%d%d", chg, suf, g, c)); Cg[g][c]->SetDirectory(nullptr); }
          else Cg[g][c]->Add(q);
          Ng[g][c] += Ny[iy][c];
          delete q;
        }
      // ---- the group shape's tail (pT >= kNoMtTailPt) = the pooled all-bins shape, scaled to
      //      the group's own integral above kNoMtTailPt (2026-10-22, user decision): the group's
      //      normalization and everything below the edge are unchanged ----
      TH1D *Cpool[2] = {}, *CgOwn[pOAbcdY::kNGroups][2] = {};
      for (int c = 0; c < 2; ++c)
      {
        for (int g = 0; g < pOAbcdY::kNGroups; ++g)
        {
          if (!Cpool[c]) { Cpool[c] = (TH1D *)Cg[g][c]->Clone(Form("nomt_Cpool_%s%s%d", chg, suf, c)); Cpool[c]->SetDirectory(nullptr); }
          else Cpool[c]->Add(Cg[g][c]);
        }
        const int bT = Cpool[c]->GetXaxis()->FindBin(kNoMtTailPt + eps), nb = Cpool[c]->GetNbinsX();
        const double ip = Cpool[c]->Integral(bT, nb);
        for (int g = 0; g < pOAbcdY::kNGroups; ++g)
        {
          CgOwn[g][c] = (TH1D *)Cg[g][c]->Clone(Form("nomt_CgOwn_%s%s%d%d", chg, suf, g, c));
          CgOwn[g][c]->SetDirectory(nullptr);
          const double ig = Cg[g][c]->Integral(bT, nb), itot = Cg[g][c]->Integral(1, nb);
          // the group's tail vs the pooled tail shape, in 10 GeV chunks (Pearson chi2)
          double chi = 0;
          int nd = 0, nEmpty = 0;
          for (int b = bT; b <= nb; b += 5)
          {
            const double vg = Cg[g][c]->Integral(b, std::min(b + 4, nb)), vp = Cpool[c]->Integral(b, std::min(b + 4, nb));
            const double ex = ip > 0 ? vp * ig / ip : 0.0;
            if (ex > 0) { chi += (vg - ex) * (vg - ex) / ex; ++nd; }
          }
          for (int b = bT; b <= nb; ++b) if (Cg[g][c]->GetBinContent(b) <= 0) ++nEmpty;
          if (ig > 0 && ip > 0)
            for (int b = bT; b <= nb; ++b)
            {
              Cg[g][c]->SetBinContent(b, Cpool[c]->GetBinContent(b) * ig / ip);
              Cg[g][c]->SetBinError(b, Cpool[c]->GetBinError(b) * ig / ip);
            }
          std::cout << Form("[NOMT-TAIL] %s %-3s %-3s m_T %s 30: anti-iso QCD-only above %.0f GeV %.1f (%.1f%% of the group), %d of %d"
                            " 2-GeV bins empty; group tail vs pooled shape chi2 %.1f / %d, p = %.2f -> pooled tail taken%s\n",
                            (lep + chg).c_str(), bins[ib], pOAbcdY::kGroupName[g], c ? ">=" : "<", kNoMtTailPt, ig,
                            itot > 0 ? 100 * ig / itot : 0.0, nEmpty, nb - bT + 1, chi, std::max(nd - 1, 1),
                            TMath::Prob(chi, std::max(nd - 1, 1)), (ig > 0 && ip > 0) ? "" : " (NOT: empty tail, kept as is)");
          // the plot: the group's own shape vs the pooled one, both per unit group integral (lab only)
          if (ib == 0)
          {
            TGraphAsymmErrors *go = new TGraphAsymmErrors(), *gp = new TGraphAsymmErrors();
            for (int b = 1; b <= nb; ++b)
            {
              const double x = xOf(CgOwn[g][c], b);
              if (itot <= 0) break;
              if (CgOwn[g][c]->GetBinContent(b) > 0)
              {
                const int n = go->GetN();
                go->SetPoint(n, x, CgOwn[g][c]->GetBinContent(b) / itot);
                go->SetPointError(n, 0, 0, CgOwn[g][c]->GetBinError(b) / itot, CgOwn[g][c]->GetBinError(b) / itot);
              }
              if (b >= bT && ip > 0) gp->SetPoint(gp->GetN(), x, Cpool[c]->GetBinContent(b) * ig / ip / itot);
            }
            auto under = [&]() {
              gPad->SetLogy();
              TLine *l = new TLine(kNoMtTailPt, 1e-4, kNoMtTailPt, 0.5);
              l->SetLineColor(kRed + 1); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw();
            };
            NoMtDiag::Draw(plotDir + Form("/nomt_tail_%s_%s_%s", tagCB.c_str(), pOAbcdY::kGroupName[g], c ? "hi" : "lo"),
                           {{go, Form("group %s: its own anti-iso shape", pOAbcdY::kGroupName[g]), kBlack, 20, 1.1},
                            {gp, "all-bins shape (used above 50)", kRed + 1, 1, 1.0, "L"}},
                           22, 102, 1e-4, 0.5, "lepton p_{T} (GeV)", "fraction of the group's anti-iso multijet / 2 GeV",
                           {Form("%s^{%s}, lab, group %s, anti-iso m_{T} %s 30", lepLatex, chgLab[ic], pOAbcdY::kGroupName[g], c ? "#geq" : "<"),
                            Form("above %.0f GeV: %.1f%% of the group; tail vs pooled shape p = %.2f", kNoMtTailPt,
                                 itot > 0 ? 100 * ig / itot : 0.0, TMath::Prob(chi, std::max(nd - 1, 1)))},
                           under, 0.50, 0.70, 0.93, 0.80);
          }
        }
      }
      // f-weighted total of a shape at slope s
      auto fTot = [&](const TH1D *C, double s) {
        double t = 0;
        for (int b = 1; b <= C->GetNbinsX(); ++b) t += C->GetBinContent(b) * std::exp(s * (xOf(C, b) - 40.0) / 10.0);
        return t;
      };
      TH1D *Tpl[pOSkim::kNY][2] = {}; // the nominal templates (kept for the stacks)
      double srQ[pOSkim::kNY][2] = {};
      for (int iy = 0; iy < pOSkim::kNY; ++iy)
      {
        const int g = pOAbcdY::GroupOf(iy);
        for (int c = 0; c < 2; ++c)
        {
          const std::string tn = Form("nomt_qcd_%s_y%d_%s", tagCB.c_str(), iy, c ? "hi" : "lo");
          for (int v = 0; v < (cfg.useFF ? 3 : 1); ++v) // 0 nominal, 1 slope Up, 2 slope Down
          {
            const double s = (v == 0) ? sHat : (v == 1 ? sUp : sDn);
            // the variations keep the group's m_T < 30 total (the level is the free T_y)
            const double norm = (v == 0) ? 1.0 : fTot(Cg[g][0], sHat) / std::max(1e-12, fTot(Cg[g][0], s));
            const double share = (Ng[g][c] > 0) ? Ny[iy][c] / Ng[g][c] : 0.0;
            TH1D *t = (TH1D *)Cg[g][c]->Clone((tn + (v == 1 ? "_ffUp" : (v == 2 ? "_ffDown" : ""))).c_str());
            t->SetDirectory(nullptr);
            t->SetTitle(Form("no-m_T multijet, %s y%d, m_T %s %.0f%s;p_{T} [GeV];Events / 2 GeV", tagCB.c_str(), iy,
                             c ? ">=" : "<", kNoMtSplitMt, v == 1 ? ", slope Up" : (v == 2 ? ", slope Down" : "")));
            for (int b = 0; b <= t->GetNbinsX() + 1; ++b)
            {
              const double fac = (b >= 1 && b <= t->GetNbinsX())
                                     ? Y[iy].F0 * norm * std::exp(s * (xOf(t, b) - 40.0) / 10.0) * share
                                     : 0.0; // no under/overflow: the fit window is pT 24-100
              t->SetBinContent(b, Cg[g][c]->GetBinContent(b) * fac);
              t->SetBinError(b, Cg[g][c]->GetBinError(b) * fac);
            }
            fout->cd();
            t->Write(t->GetName(), TObject::kOverwrite);
            if (v == 0) { Tpl[iy][c] = t; srQ[iy][c] = t->Integral(1, t->GetNbinsX()); }
            else delete t;
          }
        }
      }
      // ---- per-bin counts: the anchors D_y / H_y, the SR composition, T_y ----
      double dA[pOSkim::kNY][2] = {}, eA[pOSkim::kNY][2] = {}, dP[pOSkim::kNY][2] = {}, wP[pOSkim::kNY][2] = {}, oP[pOSkim::kNY][2] = {};
      auto writeCounts = [&](const std::string &cn, const char *what, const double *d, const double *q, const double *e,
                             const double *sq, double F0, double F0lo, double F0hi) {
        TH1D *hc = new TH1D(cn.c_str(), Form("no-m_T scheme C counts (%s, %s): anti-iso [%.2f,%.2f), m_T split %.0f, pT 0-100",
                                             tagCB.c_str(), what, cfg.isoFailLo, cfg.isoFailHi, kNoMtSplitMt),
                            15, 0.0, 15.0);
        hc->SetDirectory(nullptr);
        int b = 0;
        auto setb = [&](const char *lab, double v) { ++b; hc->GetXaxis()->SetBinLabel(b, lab); hc->SetBinContent(b, v); };
        setb("dataLo", d[0]); setb("qcdLo", q[0]); setb("ewkLo", e[0]);
        setb("dataHi", d[1]); setb("qcdHi", q[1]); setb("ewkHi", e[1]);
        setb("srQcdLo", sq[0]); setb("srQcdHi", sq[1]);
        setb("F0", F0); setb("F0lo", F0lo); setb("F0hi", F0hi);
        setb("slope", sHat); setb("slopeLo", sLo); setb("slopeHi", sHi); setb("slopeK", cfg.useFF ? kNoMtFFPriorK : 0.0);
        fout->cd();
        hc->Write(hc->GetName(), TObject::kOverwrite);
        delete hc;
      };
      double prefitTot[2] = {0, 0}, dataTot[2] = {0, 0};
      for (int iy = 0; iy < pOSkim::kNY; ++iy)
      {
        for (int c = 0; c < 2; ++c)
        {
          const int sa = c ? kAhi : kAlo, sp = c ? kPhi : kPlo;
          dA[iy][c] = S[iy][sa][0]->Integral(1, 50);
          eA[iy][c] = S[iy][sa][1]->Integral(1, 50) + S[iy][sa][2]->Integral(1, 50);
          dP[iy][c] = S[iy][sp][0]->Integral(1, 50);
          wP[iy][c] = S[iy][sp][1]->Integral(1, 50);
          oP[iy][c] = S[iy][sp][2]->Integral(1, 50);
          prefitTot[c] += srQ[iy][c];
          dataTot[c] += dP[iy][c];
        }
        writeCounts(Form("nomt_counts_%s_y%d", tagCB.c_str(), iy), Form("bin y%d", iy), dA[iy], Ny[iy], eA[iy], srQ[iy],
                    Y[iy].F0, Y[iy].F0lo, Y[iy].F0hi);
        std::cout << Form("[NOMT-BIN] %s %-3s y%-2d (%s): T_y %.4f [%.4f,%.4f] | D %.0f data (QCD %.1f) H %.0f (QCD %.1f) | lo: data %.0f, W %.1f,"
                          " other %.1f, multijet %.1f (data/prefit %.3f) | hi: data %.0f, W %.1f, other %.1f, multijet %.1f (%.0f%%)\n",
                          (lep + chg).c_str(), bins[ib], iy, pOAbcdY::kGroupName[pOAbcdY::GroupOf(iy)], Y[iy].F0, Y[iy].F0lo,
                          Y[iy].F0hi, dA[iy][0], Ny[iy][0], dA[iy][1], Ny[iy][1], dP[iy][0], wP[iy][0], oP[iy][0], srQ[iy][0],
                          (wP[iy][0] + oP[iy][0] + srQ[iy][0]) > 0 ? dP[iy][0] / (wP[iy][0] + oP[iy][0] + srQ[iy][0]) : 0.0,
                          dP[iy][1], wP[iy][1], oP[iy][1], srQ[iy][1], dP[iy][1] > 0 ? 100 * srQ[iy][1] / dP[iy][1] : 0.0);
      }
      for (int g = 0; g < pOAbcdY::kNGroups; ++g) // group sums (diagnostics; the fit uses the per-bin ones)
      {
        double d[2] = {0, 0}, q[2] = {0, 0}, e[2] = {0, 0}, sq[2] = {0, 0};
        for (int iy : pOAbcdY::BinsOf(g))
          for (int c = 0; c < 2; ++c) { d[c] += dA[iy][c]; q[c] += Ny[iy][c]; e[c] += eA[iy][c]; sq[c] += srQ[iy][c]; }
        writeCounts(Form("nomt_counts_%s_%s", tagCB.c_str(), pOAbcdY::kGroupName[g]), Form("group %s", pOAbcdY::kGroupName[g]),
                    d, q, e, sq, F0g[g], F0glo[g], F0ghi[g]);
      }
      std::cout << Form("[NOMT] %s %-3s all: SR multijet prefit m_T<30 %.1f of %.0f data (%.0f%%), m_T>=30 %.1f of %.0f (%.0f%%)\n",
                        (lep + chg).c_str(), bins[ib], prefitTot[0], dataTot[0], dataTot[0] > 0 ? 100 * prefitTot[0] / dataTot[0] : 0.0,
                        prefitTot[1], dataTot[1], dataTot[1] > 0 ? 100 * prefitTot[1] / dataTot[1] : 0.0);

      // ---- the transfer m_T < 30 -> >= 30: F (pT >= 24) in m_T slices / F(m_T < 30) ----
      // Measured POOLED over all 12 bins in the window m_T 30 - xferHi with the W
      // MC at r_W = kNoMtXferR (2026-10-22, user decision); kappa = 1 + sqrt(dev^2 +
      // stat^2) is given to every group, each with its own theta in the fit. The
      // groups' own values are diagnostics only (printed, and their chi2 vs the
      // pooled one): per group the muon measurement is pure statistics.
      const double xferHi = (lep == "ele") ? kNoMtXferMtHiEle : kNoMtXferMtHiMu;
      // slice index: 0-4 = m_T 0-10 / 10-20 / 20-30 / 30-40 / 40-50, 5 = the m_T < 30 reference
      const int slP[6] = {kPm0, kPm1, kPm2, kPmid, kPm4, kPlo}, slA[6] = {kAm0, kAm1, kAm2, kAmid, kAm4, kAlo};
      const std::vector<int> winK = (xferHi > 45.0) ? std::vector<int>{3, 4} : std::vector<int>{3};
      // F of a set of bins in a set of m_T slices with the W MC scaled by rW
      auto Fslices = [&](const std::vector<int> &ys, const std::vector<int> &ks, double rW, double *wFrac) {
        double np = 0, nf = 0, bp = 0, bf = 0, wp = 0;
        for (int iy : ys)
          for (int k : ks)
          {
            np += sumOver(S[iy][slP[k]][0], 24, 1e9); nf += sumOver(S[iy][slA[k]][0], 24, 1e9);
            wp += sumOver(S[iy][slP[k]][1], 24, 1e9);
            bp += rW * sumOver(S[iy][slP[k]][1], 24, 1e9) + sumOver(S[iy][slP[k]][2], 24, 1e9);
            bf += rW * sumOver(S[iy][slA[k]][1], 24, 1e9) + sumOver(S[iy][slA[k]][2], 24, 1e9);
          }
        if (wFrac) *wFrac = np > 0 ? wp / np : 0.0;
        return ffInterval(np, nf, std::max(bp, 1e-6), std::max(bf, 1e-6));
      };
      // sets 0-2 = the groups, 3 = all 12 bins (the pooled measurement)
      std::vector<std::vector<int>> xset;
      for (int g = 0; g < pOAbcdY::kNGroups; ++g) xset.push_back(pOAbcdY::BinsOf(g));
      xset.push_back(allY);
      const int kAll = pOAbcdY::kNGroups;
      FFInterval Isl[pOAbcdY::kNGroups + 1][6];             // per slice, at r_W = kNoMtXferR
      double xr[pOAbcdY::kNGroups + 1] = {}, xe[pOAbcdY::kNGroups + 1] = {}, xr1[pOAbcdY::kNGroups + 1] = {},
             xw[pOAbcdY::kNGroups + 1] = {};               // window ratio, its stat, at r_W = 1, W share
      for (size_t s = 0; s < xset.size(); ++s)
      {
        for (int k = 0; k < 6; ++k) Isl[s][k] = Fslices(xset[s], {k}, kNoMtXferR, nullptr);
        const FFInterval win = Fslices(xset[s], winK, kNoMtXferR, &xw[s]);
        const FFInterval win1 = Fslices(xset[s], winK, 1.0, nullptr), ref1 = Fslices(xset[s], {5}, 1.0, nullptr);
        const double ref = Isl[s][5].mle;
        xr[s] = ref > 0 ? win.mle / ref : 1.0;
        xe[s] = ref > 0 ? 0.5 * (win.hi - win.lo) / ref : 0.0;
        xr1[s] = ref1.mle > 0 ? win1.mle / ref1.mle : 1.0;
      }
      const double kPooled = 1.0 + std::sqrt((xr[kAll] - 1.0) * (xr[kAll] - 1.0) + xe[kAll] * xe[kAll]);
      double xchi2 = 0;
      for (int g = 0; g < pOAbcdY::kNGroups; ++g)
      {
        xferDev[ib][g] = xr[g] - 1.0; // the group's own value (diagnostic)
        xferErr[ib][g] = xe[g];
        kXfer[ib][g] = kPooled;       // what the fit uses: the pooled kappa, one theta per group
        xchi2 += std::pow((xr[g] - xr[kAll]) / std::max(xe[g], 1e-3), 2);
      }
      xferDevP[ib] = xr[kAll] - 1.0;
      xferErrP[ib] = xe[kAll];
      kXferMax = std::max(kXferMax, kPooled);
      const char *sll[5] = {"0-10", "10-20", "20-30", "30-40", "40-50"};
      for (size_t s = 0; s < xset.size(); ++s)
      {
        std::cout << Form("[NOMT-XFER] %s %-3s %-3s: F(m_T slice)/F(m_T<30) at r_W = %.1f:", (lep + chg).c_str(), bins[ib],
                          s < (size_t)kAll ? pOAbcdY::kGroupName[s] : "ALL", kNoMtXferR);
        for (int k = 0; k < 5; ++k)
          std::cout << Form(" [%s] %.3f", sll[k], Isl[s][5].mle > 0 ? Isl[s][k].mle / Isl[s][5].mle : 0.0);
        std::cout << Form(" | window 30-%.0f: %.3f +- %.3f (r_W = 1: %.3f; W MC = %.0f%% of the iso-pass data)%s\n", xferHi, xr[s], xe[s],
                          xr1[s], 100 * xw[s], s < (size_t)kAll ? "  [diagnostic]" : "");
      }
      std::cout << Form("[NOMT-XFER] %s %-3s POOLED: transfer kappa = 1 + sqrt(%.3f^2 + %.3f^2) = %.3f, given to every group (own theta);"
                        " groups vs pooled chi2 %.2f / 2, p = %.3f\n",
                        (lep + chg).c_str(), bins[ib], xr[kAll] - 1.0, xe[kAll], kPooled, xchi2, TMath::Prob(xchi2, 2));

      // ---- the window term (lab binning; coherent over the bins of the charge) ----
      if (ib == 0)
      {
        // with B fixed (= the per-bin m_T < 30 prediction), the predicted hi multijet
        // P(w) = Sum_g B_g x [f-weighted QCD_hi / QCD_lo of the group](w)
        auto Pw = [&](int sLo2, int sHi2) {
          double P = 0;
          for (int g = 0; g < pOAbcdY::kNGroups; ++g)
          {
            TH1D *cl = nullptr, *ch = nullptr;
            double Bg = 0;
            for (int iy : pOAbcdY::BinsOf(g))
            {
              Bg += srQ[iy][0];
              TH1D *a = qcdOf(iy, sLo2, true), *h = qcdOf(iy, sHi2, true);
              if (!cl) { cl = a; ch = h; } else { cl->Add(a); ch->Add(h); delete a; delete h; }
            }
            const double rl = fTot(cl, sHat);
            if (rl > 0) P += Bg * fTot(ch, sHat) / rl;
            delete cl; delete ch;
          }
          return P;
        };
        const double P0 = Pw(kAlo, kAhi), P1 = Pw(kA1lo, kA1hi), P2 = Pw(kA2lo, kA2hi);
        winDev[1] = P1 / P0 - 1.0;
        winDev[2] = P2 / P0 - 1.0;
        winDev[0] = std::max(std::fabs(winDev[1]), std::fabs(winDev[2]));
        kWin = 1.0 + winDev[0];
        kWinMax = std::max(kWinMax, kWin);
        std::cout << Form("[NOMT-WIN] %s: anti-iso window [%.2f,%.2f) / [%.2f,%.2f) moves the hi prediction (B fixed) by %+.1f%% / %+.1f%%"
                          " -> window kappa %.3f (coherent over the bins)\n",
                          (lep + chg).c_str(), cfg.isoFailLo, kNoMtWinMid, kNoMtWinMid, cfg.isoFailHi, 100 * winDev[1], 100 * winDev[2], kWin);
      }

      // =================== diagnostic plots ===================
      const std::string pre = plotDir + "/";
      const std::string head = Form("%s^{%s}, %s bins", lepLatex, chgLab[ic], bins[ib]);
      // (1) T_y per bin vs one level per group
      {
        TGraphAsymmErrors *g = new TGraphAsymmErrors();
        double ymax = 0;
        for (int iy = 0; iy < pOSkim::kNY; ++iy)
        {
          g->SetPoint(iy, yCenter(ib, iy), Y[iy].F0);
          g->SetPointError(iy, yHalf(ib, iy), yHalf(ib, iy), Y[iy].F0 - Y[iy].F0lo, Y[iy].F0hi - Y[iy].F0);
          ymax = std::max(ymax, Y[iy].F0hi);
        }
        std::vector<std::string> txt = {head + Form(": T_{y} at m_{T} < %.0f, slope %s", kNoMtSplitMt, cfg.useFF ? Form("%+.2f (shared)", sHat) : "0 (flat)"),
                                        "band: one level per group (the group-level model)"};
        for (int k = 0; k < pOAbcdY::kNGroups; ++k)
          txt.push_back(Form("%s: LRT one level %.1f / %d, p = %.2g", pOAbcdY::kGroupName[k], lrtT[k], (int)pOAbcdY::BinsOf(k).size() - 1,
                             TMath::Prob(lrtT[k], (int)pOAbcdY::BinsOf(k).size() - 1)));
        auto under = [&]() {
          for (int k = 0; k < pOAbcdY::kNGroups; ++k)
          {
            const std::vector<int> ys = pOAbcdY::BinsOf(k);
            const double x1 = yCenter(ib, ys.front()) - yHalf(ib, ys.front()), x2 = yCenter(ib, ys.back()) + yHalf(ib, ys.back());
            TBox *bx = new TBox(x1, F0glo[k], x2, F0ghi[k]);
            bx->SetFillColorAlpha(gColor[k], 0.25);
            bx->Draw();
            TLine *l = new TLine(x1, F0g[k], x2, F0g[k]);
            l->SetLineColor(gColor[k]); l->SetLineWidth(2); l->SetLineStyle(2); l->Draw();
          }
        };
        NoMtDiag::Draw(pre + "nomt_T_y_" + tagCB, {{g, "T_{y} per bin (68%)", kBlack, 20, 1.3}}, -2.6, 2.6, 0.0, 1.6 * ymax,
                       "y_{lab} = #minus#eta_{lab}", "T_{y} = multijet iso-pass / anti-iso at p_{T} = 40", txt, under, 0.62, 0.80, 0.93, 0.88);
      }
      // (2) the m_T < 30 channel: data/prefit with the group level vs the per-bin level;
      //     and the m_T >= 30 channel's needed multijet factor at W x 1.2
      {
        TGraphAsymmErrors *gOld = new TGraphAsymmErrors(), *gNew = new TGraphAsymmErrors(), *gHi = new TGraphAsymmErrors();
        for (int iy = 0; iy < pOSkim::kNY; ++iy)
        {
          const int g = pOAbcdY::GroupOf(iy);
          const double qN = srQ[iy][0], qO = (Y[iy].F0 > 0) ? qN * F0g[g] / Y[iy].F0 : 0.0;
          const double num = dP[iy][0] - wP[iy][0] - oP[iy][0], er = std::sqrt(std::max(dP[iy][0], 1.0));
          const double x = yCenter(ib, iy);
          if (qO > 0) { gOld->SetPoint(gOld->GetN(), x - 0.05, num / qO); gOld->SetPointError(gOld->GetN() - 1, 0, 0, er / qO, er / qO); }
          if (qN > 0) { gNew->SetPoint(gNew->GetN(), x, num / qN); gNew->SetPointError(gNew->GetN() - 1, 0, 0, er / qN, er / qN); }
          const double qH = srQ[iy][1], eh = std::sqrt(std::max(dP[iy][1], 1.0));
          if (qH > 0)
          {
            gHi->SetPoint(gHi->GetN(), x + 0.05, (dP[iy][1] - 1.2 * wP[iy][1] - oP[iy][1]) / qH);
            gHi->SetPointError(gHi->GetN() - 1, 0, 0, eh / qH, eh / qH);
          }
        }
        auto under = [&]() {
          TLine *l = new TLine(-2.6, 1, 2.6, 1);
          l->SetLineStyle(2); l->SetLineColor(kGray + 2); l->Draw();
        };
        NoMtDiag::Draw(pre + "nomt_lofit_" + tagCB,
                       {{gOld, "m_{T}<30: group level (data#minusEWK)/multijet", kGray + 2, 25, 1.4},
                        {gNew, "m_{T}<30: per-bin level", kBlack, 20, 1.3},
                        {gHi, "m_{T}#geq30: needed multijet factor (W #times 1.2)", kRed + 1, 33, 2.0}},
                       -2.6, 2.6, 0.0, 3.0, "y_{lab} = #minus#eta_{lab}", "(data #minus EWK) / prefit multijet",
                       {head + ": the per-bin level closes m_{T}<30 by construction",
                        "red: what the m_{T}#geq30 data need if W = 1.2 #times MC (the fit decides)"},
                       under, 0.40, 0.64, 0.93, 0.80);
      }
      // (3) F(pT) per group vs the shared slope (+ the per-group free slopes)
      {
        std::vector<NoMtDiag::Series> ser;
        for (int g = 0; g < pOAbcdY::kNGroups; ++g)
        {
          TGraphAsymmErrors *gr = new TGraphAsymmErrors();
          for (size_t k = 0; k + 1 < FE.size(); ++k)
          {
            double np = 0, nf = 0, bp = 0, bf = 0, xs = 0;
            for (int iy : pOAbcdY::BinsOf(g))
            {
              np += Y[iy].np[k]; nf += Y[iy].nf[k]; bp += Y[iy].bp[k]; bf += Y[iy].bf[k]; xs += Y[iy].x[k];
            }
            const FFInterval I = ffInterval(np, nf, bp, bf);
            if (F0g[g] <= 0 || nf < 3) continue;
            const double x = xs / pOAbcdY::BinsOf(g).size() + (g - 1) * 0.8;
            const int n = gr->GetN();
            gr->SetPoint(n, x, I.mle / F0g[g]);
            gr->SetPointError(n, 0, 0, (I.mle - I.lo) / F0g[g], (I.hi - I.mle) / F0g[g]);
          }
          ser.push_back({gr, Form("%s (free slope %+.2f [%+.2f,%+.2f])", pOAbcdY::kGroupName[g], sG[g], sGlo[g], sGhi[g]), gColor[g],
                         gMarker[g], gSize[g]});
        }
        auto curve = [&](double s, int col, int sty) {
          TGraphAsymmErrors *c = new TGraphAsymmErrors();
          for (int i = 0; i <= 80; ++i) { const double x = 24 + i; c->SetPoint(i, x, std::exp(s * (x - 40.0) / 10.0)); }
          c->SetLineColor(col); c->SetLineStyle(sty); c->SetLineWidth(sty == 1 ? 3 : 2);
          return c;
        };
        auto under = [&]() {
          if (!cfg.useFF) curve(sFit, kGray + 1, 1)->Draw("L SAME");
          curve(sHat, kBlack, 1)->Draw("L SAME");
          for (int g = 0; g < pOAbcdY::kNGroups; ++g) curve(sG[g], gColor[g], 2)->Draw("L SAME");
        };
        NoMtDiag::Draw(pre + "nomt_ff_" + tagCB, ser, 22, 102, 0.0, 3.5, "lepton p_{T} (GeV)", "F(p_{T}) / T_{g} at m_{T} < 30",
                       {head + (cfg.useFF ? Form(": black = the shared slope (used) %+.3f [%+.3f,%+.3f]", sFit, sFitLo, sFitHi)
                                          : Form(": black = flat (used); grey = best shared slope %+.3f [%+.3f,%+.3f]", sFit, sFitLo, sFitHi)),
                        Form("dashed: each group's own slope; LRT one slope for all: %.1f / 2, p = %.2g; flat vs shared p = %.2g", lrtS,
                             TMath::Prob(lrtS, 2), TMath::Prob(lrtFlat, 1))},
                       under, 0.16, 0.56, 0.70, 0.76);
      }
      // (4) the transfer across m_T: per group (diagnostic) and the pooled window value (sets the kappa)
      {
        std::vector<NoMtDiag::Series> ser;
        const double xm[5] = {5, 15, 25, 35, 45};
        for (int s = 0; s <= kAll; ++s)
        {
          TGraphAsymmErrors *gr = new TGraphAsymmErrors();
          const double ref = Isl[s][5].mle;
          if (ref <= 0) continue;
          for (int k = 0; k < 5; ++k)
          {
            gr->SetPoint(k, xm[k] + (s - 1.5) * 1.0, Isl[s][k].mle / ref);
            gr->SetPointError(k, 0, 0, (Isl[s][k].mle - Isl[s][k].lo) / ref, (Isl[s][k].hi - Isl[s][k].mle) / ref);
          }
          if (s < kAll)
            ser.push_back({gr, Form("%s (own window value %.2f #pm %.2f)", pOAbcdY::kGroupName[s], xr[s], xe[s]), gColor[s], gMarker[s], gSize[s]});
          else
            ser.push_back({gr, "all 12 bins", kBlack, 29, 2.0});
        }
        auto under = [&]() {
          TBox *bx = new TBox(kNoMtSplitMt, 0.0, 50.0, 3.0);
          bx->SetFillColorAlpha(kOrange - 9, 0.30); bx->Draw();
          // the pooled window measurement: a band over the window, kappa in the legend text
          TBox *pw = new TBox(kNoMtSplitMt, xr[kAll] - xe[kAll], xferHi, xr[kAll] + xe[kAll]);
          pw->SetFillColorAlpha(kGray + 2, 0.35); pw->Draw();
          TLine *pl = new TLine(kNoMtSplitMt, xr[kAll], xferHi, xr[kAll]); pl->SetLineWidth(3); pl->Draw();
          TLine *l = new TLine(0, 1, 50, 1); l->SetLineStyle(2); l->SetLineColor(kGray + 2); l->Draw();
          TLine *lm = new TLine(kNoMtSplitMt, 0.0, kNoMtSplitMt, 3.0); lm->SetLineColor(kRed + 1); lm->SetLineStyle(2); lm->SetLineWidth(2); lm->Draw();
        };
        NoMtDiag::Draw(pre + "nomt_xfer_" + tagCB, ser, 0, 50, 0.0, 3.0, "m_{T} slice (GeV)", "F(m_{T} slice) / F(m_{T} < 30), p_{T} #geq 24",
                       {head + Form(": is T the same below and above m_{T} = 30? (W MC subtracted at r_{W} = %.1f)", kNoMtXferR),
                        Form("grey band: all bins, m_{T} 30-%.0f = %.2f #pm %.2f #rightarrow #kappa_{transfer} = %.2f for every group",
                             xferHi, xr[kAll], xe[kAll], kPooled),
                        Form("groups vs the all-bins value: p = %.2f", TMath::Prob(xchi2, 2))},
                       under, 0.16, 0.56, 0.66, 0.80);
      }
      // (5) per-bin anti-iso shape vs the group shape (both channels), + the p-value summary
      {
        TGraphAsymmErrors *gpL = new TGraphAsymmErrors(), *gpH = new TGraphAsymmErrors();
        for (int c = 0; c < 2; ++c)
          for (int g = 0; g < pOAbcdY::kNGroups; ++g)
          {
            std::vector<NoMtDiag::Series> ser;
            const std::vector<int> ys = pOAbcdY::BinsOf(g);
            TH1D *qg = qcdOf(ys.front(), c ? kAhi : kAlo, false);
            for (size_t j = 1; j < ys.size(); ++j) { TH1D *q = qcdOf(ys[j], c ? kAhi : kAlo, false); qg->Add(q); delete q; }
            const double Qg = sumOver(qg, 24, 1e9);
            int iser = 0;
            for (int iy : ys)
            {
              TH1D *q = qcdOf(iy, c ? kAhi : kAlo, false), *dd = S[iy][c ? kAhi : kAlo][0];
              const double Qy = sumOver(q, 24, 1e9);
              TGraphAsymmErrors *gr = new TGraphAsymmErrors();
              double chi2 = 0;
              int ndf = 0;
              for (size_t k = 0; k + 1 < FE.size(); ++k)
              {
                const double vy = sumOver(q, FE[k], FE[k + 1]), vg = sumOver(qg, FE[k], FE[k + 1]);
                const double ny = sumOver(dd, FE[k], FE[k + 1]);
                if (Qy <= 0 || Qg <= 0 || vg <= 0 || ny < 1) continue;
                const double r = (vy / Qy) / (vg / Qg), er = std::sqrt(ny) / Qy / (vg / Qg);
                chi2 += (r - 1.0) * (r - 1.0) / (er * er);
                ++ndf;
                const int n = gr->GetN();
                gr->SetPoint(n, 0.5 * (FE[k] + std::min(FE[k + 1], 110.0)) + (iser - 0.5 * (ys.size() - 1)) * 0.5, r);
                gr->SetPointError(n, 0, 0, er, er);
              }
              const double p = ndf > 1 ? TMath::Prob(chi2, ndf - 1) : 1.0;
              (c ? gpH : gpL)->SetPoint((c ? gpH : gpL)->GetN(), yCenter(ib, iy) + (c ? 0.05 : -0.05), std::max(p, 1e-6));
              if (ib == 0)
                ser.push_back({gr, Form("y%d: #chi^{2} %.1f / %d, p = %.2g", iy, chi2, std::max(ndf - 1, 1), p),
                               kBinColor[iser % 8], kBinMarker[iser % 8], 1.2});
              else
                delete gr;
              ++iser;
              delete q;
            }
            delete qg;
            if (ib == 0)
            {
              auto under = [&]() { TLine *l = new TLine(22, 1, 112, 1); l->SetLineStyle(2); l->SetLineColor(kGray + 2); l->Draw(); };
              NoMtDiag::Draw(pre + Form("nomt_shape_%s_%s_%s", tagCB.c_str(), pOAbcdY::kGroupName[g], c ? "hi" : "lo"), ser, 22, 112, 0.0,
                             3.0, "lepton p_{T} (GeV)", "per-bin anti-iso shape / group shape",
                             {head + Form(", group %s, anti-iso m_{T} %s 30: the template uses the GROUP shape", pOAbcdY::kGroupName[g],
                                          c ? "#geq" : "<")},
                             under, 0.58, 0.48, 0.93, 0.84);
            }
          }
        auto under = [&]() {
          TLine *l = new TLine(-2.6, 0.05, 2.6, 0.05); l->SetLineStyle(2); l->SetLineColor(kRed + 1); l->Draw();
          gPad->SetLogy();
        };
        NoMtDiag::Draw(pre + "nomt_shape_p_" + tagCB,
                       {{gpL, "anti-iso m_{T} < 30", kBlack, 20, 1.3}, {gpH, "anti-iso m_{T} #geq 30", kAzure + 2, 21, 1.5}},
                       -2.6, 2.6, 1e-6, 5.0, "y_{lab} = #minus#eta_{lab}", "p-value: bin shape = group shape",
                       {head + ": below the red line = the bin's anti-iso shape differs from its group's"}, under, 0.62, 0.20, 0.93, 0.32);
      }
      // (6) per-group iso-pass m_T spectrum: data vs EWK (r = 1) + multijet (T_y f(pT) x anti-iso per bin)
      for (int g = 0; g < pOAbcdY::kNGroups; ++g)
      {
        TH1D *hd = nullptr, *hw = nullptr, *ho = nullptr, *hq = nullptr;
        for (int iy : pOAbcdY::BinsOf(g))
        {
          if (!hd)
          {
            hd = (TH1D *)M[iy][0]->Clone(Form("nomt_md_%s_%d", tagCB.c_str(), g)); hd->SetDirectory(nullptr); hd->Reset();
            hw = (TH1D *)hd->Clone(Form("nomt_mw_%s_%d", tagCB.c_str(), g)); hw->SetDirectory(nullptr);
            ho = (TH1D *)hd->Clone(Form("nomt_mo_%s_%d", tagCB.c_str(), g)); ho->SetDirectory(nullptr);
            hq = (TH1D *)hd->Clone(Form("nomt_mq_%s_%d", tagCB.c_str(), g)); hq->SetDirectory(nullptr);
          }
          hd->Add(M[iy][0]); hw->Add(M[iy][1]); ho->Add(M[iy][2]);
          for (int jz = 1; jz <= Aqp[iy]->GetNbinsY(); ++jz)
          {
            double v = 0;
            for (int jy = 1; jy <= Aqp[iy]->GetNbinsX(); ++jy)
            {
              const double x = Aqp[iy]->GetXaxis()->GetBinCenter(jy);
              const double xf = (x < 25.0) ? 25.5 : (x > 100 ? 105.0 : x);
              v += Aqp[iy]->GetBinContent(jy, jz) * Y[iy].F0 * std::exp(sHat * (xf - 40.0) / 10.0);
            }
            hq->SetBinContent(jz, hq->GetBinContent(jz) + v);
          }
        }
        for (int b = 1; b <= hq->GetNbinsX(); ++b) if (hq->GetBinContent(b) < 0) hq->SetBinContent(b, 0.0);
        // per GeV (the m_T bins are 10-80 GeV wide)
        for (TH1D *h : {hd, hw, ho, hq}) h->Scale(1.0, "width");
        const double dLo = hd->Integral(1, hd->FindBin(kNoMtSplitMt - eps), "width"), dHi = hd->Integral(hd->FindBin(kNoMtSplitMt + eps), hd->GetNbinsX(), "width");
        const double wLo = hw->Integral(1, hw->FindBin(kNoMtSplitMt - eps), "width"), wAll = hw->Integral("width");
        PlotStyle ps;
        ps.drawOpt = "hist";
        ps.logy = true;
        ps.normBkgToData = false;
        ps.pullPad = true;
        ps.boxX1 = 0.50; ps.boxY1 = 0.58; ps.boxY2 = 0.80;          // info box under the header
        ps.legX1 = 0.50; ps.legY1 = 0.36; ps.legX2 = 0.92; ps.legY2 = 0.56; // legend under the box
        ps.xRangeLo = 0.0;
        ps.xRangeHi = 200.0;
        PlotTuner tuner = [&](TCanvas *cc, TH1 *h) {
          if (!h) return;
          h->SetMinimum(0.05);
          h->SetMaximum(200.0 * h->GetMaximum());
          (void)cc;
        };
        SaveNicePlot1D_WithBkg(hd, {hw, ho, hq}, {"W (MC, r = 1)", "DY, #tau, t#bar{t} (MC)", "multijet (T_{y} f(p_{T}) #times anti-iso)"},
                               pre + Form("nomt_mt_%s_%s", tagCB.c_str(), pOAbcdY::kGroupName[g]), "m_{T} (GeV)", "Events / GeV", "",
                               head + Form(", group %s", pOAbcdY::kGroupName[g]), "iso-pass, split at m_{T} = 30",
                               {Form("data: m_{T}<30 %.0f, #geq30 %.0f", dLo, dHi), Form("W MC below 30: %.1f%%", wAll > 0 ? 100 * wLo / wAll : 0.0)},
                               ps, tuner);
        delete hd; delete hw; delete ho; delete hq;
      }
      // (7) per-group pT stacks of the two SR channels (W at r = 1; multijet = the templates)
      for (int g = 0; g < pOAbcdY::kNGroups; ++g)
        for (int c = 0; c < 2; ++c)
        {
          TH1D *hd = nullptr, *hw = nullptr, *ho = nullptr, *hq = nullptr;
          for (int iy : pOAbcdY::BinsOf(g))
          {
            const int sp = c ? kPhi : kPlo;
            if (!hd)
            {
              hd = (TH1D *)S[iy][sp][0]->Clone(Form("nomt_sd_%s_%d_%d", tagCB.c_str(), g, c)); hd->SetDirectory(nullptr); hd->Reset();
              hw = (TH1D *)hd->Clone(Form("nomt_sw_%s_%d_%d", tagCB.c_str(), g, c)); hw->SetDirectory(nullptr);
              ho = (TH1D *)hd->Clone(Form("nomt_so_%s_%d_%d", tagCB.c_str(), g, c)); ho->SetDirectory(nullptr);
              hq = (TH1D *)hd->Clone(Form("nomt_sq_%s_%d_%d", tagCB.c_str(), g, c)); hq->SetDirectory(nullptr);
            }
            hd->Add(S[iy][sp][0]); hw->Add(S[iy][sp][1]); ho->Add(S[iy][sp][2]); hq->Add(Tpl[iy][c]);
          }
          for (TH1D *h : {hd, hw, ho, hq}) { h->SetBinContent(h->GetNbinsX() + 1, 0.0); h->SetBinContent(0, 0.0); }
          const double D = hd->Integral(1, 50), W = hw->Integral(1, 50), O = ho->Integral(1, 50), Q = hq->Integral(1, 50);
          PlotStyle ps;
          ps.drawOpt = "hist";
          ps.logy = true;
          ps.normBkgToData = false;
          ps.pullPad = true;
          ps.boxX1 = 0.50; ps.boxY1 = 0.58; ps.boxY2 = 0.80;
          ps.legX1 = 0.50; ps.legY1 = 0.36; ps.legX2 = 0.92; ps.legY2 = 0.56;
          ps.xRangeLo = 24.0;
          ps.xRangeHi = 100.0;
          PlotTuner tuner = [&](TCanvas *cc, TH1 *h) {
            if (!h) return;
            h->SetMinimum(0.3);
            h->SetMaximum(50.0 * h->GetMaximum());
            (void)cc;
          };
          SaveNicePlot1D_WithBkg(hd, {hw, ho, hq}, {"W (MC, r = 1)", "DY, #tau, t#bar{t} (MC)", "multijet (per-bin T_{y})"},
                                 pre + Form("nomt_sr_%s_%s_%s", tagCB.c_str(), pOAbcdY::kGroupName[g], c ? "hi" : "lo"),
                                 "lepton p_{T} (GeV)", "Events / 2 GeV", "", head + Form(", group %s", pOAbcdY::kGroupName[g]),
                                 Form("SR m_{T} %s %.0f (prefit)", c ? "#geq" : "<", kNoMtSplitMt),
                                 {Form("data %.0f", D), Form("W %.1f, other %.1f", W, O), Form("multijet %.1f (%.0f%%)", Q, D > 0 ? 100 * Q / D : 0.0),
                                  Form("data / prefit %.3f", (W + O + Q) > 0 ? D / (W + O + Q) : 0.0)},
                                 ps, tuner);
          delete hd; delete hw; delete ho; delete hq;
        }

      for (int g = 0; g < pOAbcdY::kNGroups; ++g) for (int c = 0; c < 2; ++c) { delete Cg[g][c]; delete CgOwn[g][c]; }
      for (int c = 0; c < 2; ++c) delete Cpool[c];
      for (int iy = 0; iy < pOSkim::kNY; ++iy)
      {
        for (int s = 0; s < kNSl; ++s) for (int k = 0; k < 3; ++k) delete S[iy][s][k];
        for (int k = 0; k < 3; ++k) delete M[iy][k];
        delete Aqp[iy];
        for (int c = 0; c < 2; ++c) delete Tpl[iy][c];
      }
    }
    // ---- the residual kappas of this charge (labeled; read by plotting/combine_input_nomt.C) ----
    // kappaXfer_<bin>_<group> = the POOLED kappa (the same for the three groups; the fork
    // gives each group its own theta); xferDev/xferErr_<bin>_<group> = the group's own
    // window value (diagnostic); xferDevPooled/xferErrPooled_<bin> = what sets the kappa
    const int nk = 5 + 2 * (pOAbcdY::kNGroups * 3 + 2);
    TH1D *hk = new TH1D(Form("nomt_kappa_%s%s", lep.c_str(), chg), "no-m_T scheme C: the residual lnN of the m_T>=30 multijet", nk, 0, nk);
    hk->SetDirectory(nullptr);
    int b = 0;
    auto setk = [&](const std::string &lab, double v) { ++b; hk->GetXaxis()->SetBinLabel(b, lab.c_str()); hk->SetBinContent(b, v); };
    setk("kappaWin", kWin); setk("windowDev", winDev[0]); setk("winDevLo", winDev[1]); setk("winDevHi", winDev[2]);
    setk("xferMtHi", (lep == "ele") ? kNoMtXferMtHiEle : kNoMtXferMtHiMu);
    for (int ib = 0; ib < 2; ++ib)
    {
      for (int g = 0; g < pOAbcdY::kNGroups; ++g)
      {
        setk(Form("kappaXfer_%s_%s", bins[ib], pOAbcdY::kGroupName[g]), kXfer[ib][g]);
        setk(Form("xferDev_%s_%s", bins[ib], pOAbcdY::kGroupName[g]), xferDev[ib][g]);
        setk(Form("xferErr_%s_%s", bins[ib], pOAbcdY::kGroupName[g]), xferErr[ib][g]);
      }
      setk(Form("xferDevPooled_%s", bins[ib]), xferDevP[ib]);
      setk(Form("xferErrPooled_%s", bins[ib]), xferErrP[ib]);
    }
    fout->cd();
    hk->Write(hk->GetName(), TObject::kOverwrite);
    delete hk;
  }
  std::cout << Form("[NOMT] => residual lnN of the m_T>=30 multijet: window kappa up to %.3f (coherent per charge), transfer kappa up to %.3f"
                    " (measured on all bins, one theta per group); the values go to the _mtsplit.txt sidecar via plotting/combine_input_nomt.C\n",
                    kWinMax, kXferMax);
  gErrorIgnoreLevel = errLevel0;
}

#endif // PO_QCD_ABCD_NOMT_H
