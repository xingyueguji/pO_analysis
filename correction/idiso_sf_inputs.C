// correction/idiso_sf_inputs.C
//
// Combine inputs of the ELECTRON ID+ISO SF fit (2026-09-24) + the composition
// record. Reads the separate skim (rootfile/idiso_sf_ele/skim_<sample>.root,
// correction/idiso_sf_skim.C); definition and names in idiso_sf_common.h.
//
// One W input file per (scheme, sideband window) -- the fitted variable is the
// lepton pT at m_T > 40 (leppt_mt40), the QCD the nominal IN-FIT ABCD in the
// (m_T x relIso) plane (idiso_sf_common.h; user 2026-09-25):
//   rootfile/idiso_sf_ele/combine_input_idiso_<scheme>_leppt_mt40[_sblo|_sbhi].root
//     W{p,m}_{pass,fail}_k<k>/      the SR (m_T > 40, pT): data_obs signal z ztau wtau qcd
//     W{p,m}_{pass,fail}_k<k>_CRB/  1-bin, the category at m_T < 30: data_obs qcd signal wtau z ztau
//     W{p,m}_{pass,fail}_k<k>_CRC/  1-bin, the sideband at m_T > 40:  data_obs qcd ewk
//     W{p,m}_{pass,fail}_k<k>_CRD/  1-bin, the sideband at m_T < 30:  data_obs qcd ewk
//     Z_PP/                         data_obs signal w wtau ztau bkg (both TnP legs
//                                   pass: the DY anchor of the W fit, r_Z)
//   + <same>_meta.txt               scheme, disc, sbwin, nbins, bin edges, x title
//                                   (read by the fork's make_pO_idiso_sf_cards.sh)
// and ONE Z tag-and-probe input (the Z TnP is its OWN fit -- see below):
//   rootfile/idiso_sf_ele/combine_input_idiso_ztnp.root  Z_PP/, Z_PF/ (+ _meta.txt,
//                                   scheme ztnp)
// WHY a separate Z fit: sharing r_Z between the Z_PP peak and the DY of the W
// channels lets the W side pull it -- and with it the Z TnP efficiency: in one
// likelihood the fitted SF_Z moved 0.96 -> 1.01 across the W variants although
// their Z channels are identical (2026-09-25 stand-in).
// scheme = incl | abseta | pt ; window = full (relIso 0.3-1.0, no suffix) |
// sblo (0.3-0.6) | sbhi (0.6-1.0). The fit runs in the Combine fork:
// test/run_pO_idiso_sf.sh.
//
// Templates (all ABSOLUTE, k_s = A.sigma.L/N_gen from skim/mc_norm.h, no area
// normalization -- the nominal convention):
//   signal = W -> e nu (Wp + Wm samples, in the reco-charge region), z = DY -> ee,
//   ztau = DY -> tautau, wtau = W -> tau nu, ewk = all of them, data_obs = data.
//   CR qcd = the prefit QCD count, data - EWK at r = 1 (the anchor of the fit's
//            free scale sB / sC / sD -- the in-fit ABCD floats all three)
//   SR qcd = the sideband pT shape at m_T > 40 (data - EWK, negative bins
//            clamped to 0) scaled to A0 = B0 C0 / D0 -- the prefit ABCD; the fit
//            multiplies it by sB sC / sD (a formula rateParam), so nothing of
//            the prefit normalization is assumed. pass-like channels take the
//            sideband with eleMVAIdWP90 (sbid), fail-like ones any ID (sball).
//   Z_PP, Z_PF: signal = DY -> ee, w = W -> e nu, wtau, ztau (MC, absolute),
//            bkg = FLAT over 60-120 with the same-sign data count of the
//            category (>= 1 event) as its prefit size; the fit floats it.
//
// The Z tag-and-probe record ([ZTNP] lines): the per-category counts, the
// same-sign-subtracted cut-and-count efficiency in data and DY MC (the
// cross-check of the fitted one), and [EGMPRED] = the EGM wp90iso SF averaged
// over the passing MC probes / electrons, i.e. what the measured SF should be
// if the pp SF describes our selection (Z probes inclusive; W per coarse bin,
// m_T > 40).
//
// The ABCD record ([ABCD] lines, per W channel): the prefit B0, C0, D0 and A0,
// and in the SR the prefit closure (data - EWK at r = 1) / A0 and the W / A0 --
// the go/no-go of the fail-category W extraction (in the fail channels W is a
// few % of the QCD, so an ABCD non-closure of that size decides the W count).
//
// Run from correction/ (after the skim):
//   ./run_idiso_sf.sh inputs                   (keeps logs/idiso_sf_ele_inputs.log)
//   root -l -b -q 'idiso_sf_inputs.C+()'       (bare)

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"
#include "TCanvas.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
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

// sideband windows of the in-fit ABCD (regions C, D and the SR QCD shape)
const int         kNWin = 3;
const char *const kWinName[kNWin]   = {"full", "sblo", "sbhi"};
const char *const kWinSuffix[kNWin] = {"", "_sblo", "_sbhi"};
const char *const kWinLabel[kNWin]  = {"relIso 0.3-1.0", "relIso 0.3-0.6", "relIso 0.6-1.0"};

// sample indices (idiso_sf_common.h order: Data Wp Wm DY DYtau Wptau Wmtau)
const std::vector<int> kSigS  = {1, 2};
const std::vector<int> kZS    = {3};
const std::vector<int> kZtauS = {4};
const std::vector<int> kWtauS = {5, 6};
const std::vector<int> kEwkS  = {1, 2, 3, 4, 5, 6};

TFile *gF[kNSample] = {};
double gK[kNSample] = {};
int    gClone = 0;

// Clone (never mutate a file-owned histogram) and apply k_s to MC.
TH1D *Get(int is, const std::string &name)
{
  TH1D *src = gF[is] ? (TH1D *)gF[is]->Get(name.c_str()) : nullptr;
  if (!src)
  {
    std::cerr << "[FATAL] " << kSampleName[is] << ": missing " << name << "\n";
    return nullptr;
  }
  TH1D *h = (TH1D *)src->Clone(Form("%s_%s_c%d", name.c_str(), kSampleName[is], gClone++));
  h->SetDirectory(nullptr);
  if (is > 0) h->Scale(gK[is]);
  return h;
}
TH1D *Sum(const std::vector<int> &samples, const std::vector<std::string> &names)
{
  TH1D *out = nullptr;
  for (int is : samples)
    for (const auto &n : names)
    {
      TH1D *h = Get(is, n);
      if (!h) return nullptr;
      if (!out) out = h;
      else { out->Add(h); delete h; }
    }
  return out;
}

double IntegralBelow(const TH1D *h, double xMax)
{
  if (!h) return 0.0;
  const int b = h->GetXaxis()->FindBin(xMax - 1e-9);
  return h->Integral(1, b);
}
double IntegralAbove(const TH1D *h, double xMin)
{
  if (!h) return 0.0;
  const int b = h->GetXaxis()->FindBin(xMin + 1e-9);
  return h->Integral(b, h->GetNbinsX() + 1);
}

// a 1-bin counting template of a control region
TH1D *Count(const std::string &name, double v)
{
  TH1D *h = new TH1D(Form("%s_c%d", name.c_str(), gClone++), "", 1, 0.0, 1.0);
  h->SetDirectory(nullptr);
  h->SetBinContent(1, v);
  h->SetBinError(1, std::sqrt(std::fabs(v)));
  return h;
}

// The ABCD regions of one W channel: B = the category at m_T < 30, C = the
// sideband at m_T > 40, D = the sideband at m_T < 30
enum Reg { kRB = 0, kRC, kRD, kNReg };
const char *const kRegName[kNReg] = {"CRB", "CRC", "CRD"};

struct WChan
{
  std::string dir;
  // the SR: m_T > 40, the lepton pT (the fitted variable)
  TH1D *data = nullptr, *sig = nullptr, *z = nullptr, *ztau = nullptr, *wtau = nullptr, *qcd = nullptr;
  // the three control regions (Sum_w, MC at k_s); qcd0 = data - EWK at r = 1
  double cData[kNReg] = {}, cSig[kNReg] = {}, cWtau[kNReg] = {}, cZ[kNReg] = {}, cZtau[kNReg] = {},
         cEwk[kNReg] = {}, cQcd0[kNReg] = {};
  double A0 = 0.0;
  bool ok() const { return data && sig && z && ztau && wtau && qcd; }
  void Delete() { for (TH1D *h : {data, sig, z, ztau, wtau, qcd}) delete h; }
};

// The W channel (category c = kPass|kFail, charge q, scheme s, bin k) with its
// in-fit-ABCD control regions for sideband window iw (idiso_sf_common.h).
WChan MakeW(int c, int q, int s, int k, int iw)
{
  WChan ch;
  ch.dir = WDir(c, q, k);
  const std::string n = HName(kDiscFit, c, q, s, k);
  ch.data = Get(0, n);
  ch.sig  = Sum(kSigS, {n});
  ch.z    = Sum(kZS, {n});
  ch.ztau = Sum(kZtauS, {n});
  ch.wtau = Sum(kWtauS, {n});

  // the sideband categories of this window
  const int lo = (c == kPass) ? kSbIdLo : kSbAllLo;
  const int hi = (c == kPass) ? kSbIdHi : kSbAllHi;
  std::vector<int> sbc;
  if (iw == 0 || iw == 1) sbc.push_back(lo);
  if (iw == 0 || iw == 2) sbc.push_back(hi);
  std::vector<std::string> sbMt, sbPt;
  for (int cc : sbc) { sbMt.push_back(HName(kDiscMt, cc, q, s, k)); sbPt.push_back(HName(kDiscFit, cc, q, s, k)); }

  // the counts, from the m_T histograms (30 and 40 are bin edges of the 2.5 GeV axis)
  const std::string nMt = HName(kDiscMt, c, q, s, k);
  TH1D *bD = Get(0, nMt), *bS = Sum(kSigS, {nMt}), *bWt = Sum(kWtauS, {nMt}), *bZ = Sum(kZS, {nMt}),
       *bZt = Sum(kZtauS, {nMt});
  TH1D *sD = Sum({0}, sbMt), *sE = Sum(kEwkS, sbMt);
  TH1D *pD = Sum({0}, sbPt), *pE = Sum(kEwkS, sbPt);
  if (!ch.data || !ch.sig || !ch.z || !ch.ztau || !ch.wtau || !bD || !bS || !bWt || !bZ || !bZt || !sD || !sE ||
      !pD || !pE)
    return ch;
  ch.cData[kRB] = IntegralBelow(bD, kCRMtMax);
  ch.cSig[kRB]  = IntegralBelow(bS, kCRMtMax);
  ch.cWtau[kRB] = IntegralBelow(bWt, kCRMtMax);
  ch.cZ[kRB]    = IntegralBelow(bZ, kCRMtMax);
  ch.cZtau[kRB] = IntegralBelow(bZt, kCRMtMax);
  ch.cEwk[kRB]  = ch.cSig[kRB] + ch.cWtau[kRB] + ch.cZ[kRB] + ch.cZtau[kRB];
  ch.cData[kRC] = IntegralAbove(sD, kSRMtMin);
  ch.cEwk[kRC]  = IntegralAbove(sE, kSRMtMin);
  ch.cData[kRD] = IntegralBelow(sD, kCRMtMax);
  ch.cEwk[kRD]  = IntegralBelow(sE, kCRMtMax);
  for (int r = 0; r < kNReg; ++r)
  {
    ch.cQcd0[r] = ch.cData[r] - ch.cEwk[r];
    if (ch.cQcd0[r] <= 0)
    {
      std::cout << Form("[WARN] %s %s: prefit QCD count %.2f <= 0 -> 1e-3 (the fit's free scale takes it from there)\n",
                        ch.dir.c_str(), kRegName[r], ch.cQcd0[r]);
      ch.cQcd0[r] = 1e-3;
    }
  }
  ch.A0 = ch.cQcd0[kRB] * ch.cQcd0[kRC] / ch.cQcd0[kRD];

  // the SR QCD template: the sideband pT shape at m_T > 40, scaled to A0 over the fitted bins
  TH1D *t = (TH1D *)pD->Clone(Form("qcd_%s_c%d", ch.dir.c_str(), gClone++));
  t->SetDirectory(nullptr);
  t->Add(pE, -1.0);
  for (int b = 0; b <= t->GetNbinsX() + 1; ++b)
    if (t->GetBinContent(b) < 0) t->SetBinContent(b, 0.0);
  const double have = t->Integral(1, t->GetNbinsX());
  if (have > 0) t->Scale(ch.A0 / have);
  else std::cout << "[WARN] " << ch.dir << ": empty sideband pT shape -> the SR QCD template is empty\n";
  ch.qcd = t;
  for (TH1D *h : {bD, bS, bWt, bZ, bZt, sD, sE, pD, pE}) delete h;
  return ch;
}

// A Z tag-and-probe category (zc = kZPP | kZPF): OS templates + the SS data
struct ZChan
{
  std::string dir;
  TH1D *data = nullptr, *zsig = nullptr, *w = nullptr, *wtau = nullptr, *ztau = nullptr, *bkg = nullptr,
       *ssData = nullptr, *ssDY = nullptr;
  bool ok() const { return data && zsig && w && wtau && ztau && bkg && ssData && ssDY; }
  void Delete() { for (TH1D *h : {data, zsig, w, wtau, ztau, bkg, ssData, ssDY}) delete h; }
};
ZChan MakeZ(int zc)
{
  ZChan ch;
  ch.dir    = ZDir(zc);
  const std::string os = ZHistName(zc, false), ss = ZHistName(zc, true);
  ch.data   = Get(0, os);
  ch.zsig   = Sum(kZS, {os});
  ch.w      = Sum(kSigS, {os});
  ch.wtau   = Sum(kWtauS, {os});
  ch.ztau   = Sum(kZtauS, {os});
  ch.ssData = Get(0, ss);
  ch.ssDY   = Sum(kZS, {ss});
  if (!ch.data || !ch.ssData) return ch;
  // flat background, prefit size = the SS data count (>= 1): the fit floats it
  ch.bkg = new TH1D(Form("bkg_%s_c%d", ch.dir.c_str(), gClone++), "", kZNb, kZLo, kZHi);
  ch.bkg->SetDirectory(nullptr);
  const double nb = std::max(ch.ssData->Integral(), 1.0);
  for (int b = 1; b <= kZNb; ++b) ch.bkg->SetBinContent(b, nb / kZNb);
  return ch;
}

// (the EGM wp90iso table -- LoadEgm, FindEgm, EgmPrediction -- lives in idiso_sf_common.h)

void WriteHist(TDirectory *d, TH1D *h, const char *name, const std::string &file)
{
  if (!h) return;
  d->cd();
  TH1D *c = (TH1D *)h->Clone(name);
  c->SetDirectory(d);
  // text2workspace cannot build a pdf from an ALL-ZERO shape (the nominal
  // writer's floor): an empty MC template gets 1e-6 in its central bin
  if (std::string(name) != "data_obs" && c->Integral() <= 0.0)
  {
    std::cout << "[WARN] " << file << " " << d->GetName() << "/" << name
              << " is empty -> flooring central bin to 1e-6 for Combine\n";
    c->SetBinContent(c->GetNbinsX() / 2, 1e-6);
  }
  c->Write(name, TObject::kOverwrite);
}

// One Z TnP category into its dir of `fo`; `plot` draws the prefit stack once.
bool WriteZ(TFile &fo, int zc, const std::string &fout, bool plot)
{
  ZChan z = MakeZ(zc);
  if (!z.ok()) return false;
  TDirectory *zd = fo.mkdir(z.dir.c_str());
  WriteHist(zd, z.data, "data_obs", fout);
  WriteHist(zd, z.zsig, "signal", fout);
  WriteHist(zd, z.w, "w", fout);
  WriteHist(zd, z.wtau, "wtau", fout);
  WriteHist(zd, z.ztau, "ztau", fout);
  WriteHist(zd, z.bkg, "bkg", fout);
  if (plot)
  {
    PlotStyle ps;
    ps.drawOpt = "hist";
    ps.normBkgToData = false;
    ps.pullPad = true;
    ps.boxX1 = 0.17; ps.boxX2 = 0.54; ps.boxY1 = 0.46; ps.boxY2 = 0.76;
    ps.headerX = 0.58;  // the default 0.6 runs the TnP titles off the frame
    ps.titleSize = 0.042;
    std::vector<std::string> box = {Form("Data: %.0f (SS %.0f)", z.data->Integral(), z.ssData->Integral()),
                                    Form("DY MC: %.1f", z.zsig->Integral())};
    SaveNicePlot1D_WithBkg(z.data, {z.zsig, z.bkg, z.w, z.wtau, z.ztau},
                           {"Z #rightarrow ee", "flat bkg (SS count)", "W^{+}/W^{-}", "W #tau", "DY #tau"},
                           Form("%s/prefit/%s", kPlotDir, z.dir.c_str()), "m_{ee} (GeV)", "Events / 1.0 GeV",
                           zc == kZPP ? "Z TnP: both pass" : "Z TnP: probe fails", "WP90 ID && WP90 iso",
                           "prefit, r = 1", box, ps);
  }
  z.Delete();
  return true;
}

std::string BinLabel(int s, int k)
{
  const Scheme &sc = kSchemes[s];
  if (sc.var == kVarIncl) return "inclusive (p_{T} > 25, |#eta_{SC}| < 2.4)";
  if (sc.var == kVarAbsEta) return Form("%.4g #leq |#eta_{SC}| < %.4g", sc.lo[k], sc.hi[k]);
  if (sc.hi[k] >= kPtInf) return Form("p_{T} #geq %.0f GeV", sc.lo[k]);
  return Form("%.0f #leq p_{T} < %.0f GeV", sc.lo[k], sc.hi[k]);
}

void PrefitPlot(const WChan &ch, const std::string &out, const std::string &main, const std::string &sub)
{
  PlotStyle ps;
  ps.drawOpt = "hist";
  ps.logy = true;
  ps.boxY1 = 0.595; ps.boxY2 = 0.795;
  ps.normBkgToData = false; // absolute prefit (r = 1, all ABCD scales 1)
  ps.pullPad = true;
  ps.xRangeLo = 24.0; ps.xRangeHi = 200.0; // from the 2 GeV edge enclosing the 25 GeV cut, as the nominal
  ps.headerX = 0.50;  // the default 0.6 runs the bin label off the frame
  ps.titleSize = 0.042;
  PlotTuner tuner = [&](TCanvas *c, TH1 *h)
  {
    (void)c;
    if (!h) return;
    h->SetMinimum(1.0);
    h->SetMaximum(10.25 * h->GetMaximum());
  };
  double tot = ch.sig->Integral() + ch.z->Integral() + ch.ztau->Integral() + ch.wtau->Integral() + ch.qcd->Integral();
  std::vector<std::string> box = {Form("Data: %.0f", ch.data->Integral()),
                                  Form("W signal MC: %.1f", ch.sig->Integral()),
                                  Form("Stack: %.1f", tot)};
  SaveNicePlot1D_WithBkg(ch.data, {ch.sig, ch.z, ch.ztau, ch.wtau, ch.qcd},
                         {"W signal", "DY", "DY #tau", "W #tau", "QCD (in-fit ABCD, prefit)"}, out,
                         kDiscs[kDiscFit].axisTitle, "Events / 2.0 GeV", main, sub, "prefit, m_{T} > 40", box, ps, tuner);
}

} // namespace

int idiso_sf_inputs()
{
  gErrorIgnoreLevel = kWarning; // the ~70 "png file ... has been created" lines would bury the record
  // ---------------- inputs + k_s ----------------
  for (int is = 0; is < kNSample; ++is)
  {
    const std::string fn = SkimFile(kSampleName[is]);
    gF[is] = TFile::Open(fn.c_str(), "READ");
    if (!gF[is] || gF[is]->IsZombie())
    {
      std::cerr << "[FATAL] cannot open " << fn << " -- run ./run_idiso_sf.sh skim first\n";
      return 2;
    }
    gK[is] = (is == 0) ? 1.0 : pONorm::MCScale(kSampleNormLabel[is]);
    std::cout << Form("[INPUT] %-6s k_s = %.6g  <- %s\n", kSampleName[is], gK[is], fn.c_str());
  }
  gSystem->mkdir(kRootDir, kTRUE);
  gSystem->mkdir(Form("%s/prefit", kPlotDir), kTRUE);

  // ---------------- eps_MC: templates vs gen-matched, eps(relIso < 0.3) ----------------
  std::cout << "\n==== eps_MC = W(pass) / W(pass + fail) in W -> e nu MC (k_s-weighted Wp + Wm) ====\n";
  std::cout << "  SR = the fitted signal templates (m_T > 40: what the fit's eps_MC is); all m_T = every W MC\n"
               "  event of the category; gen-matched = the leading electron is a prompt W electron (all m_T)\n";
  for (int s = 0; s < kNScheme; ++s)
    for (int k = 0; k < kSchemes[s].n; ++k)
      for (int q = 0; q <= kNChg; ++q) // q == kNChg: both charges
      {
        double sp = 0, sf = 0, gp = 0, gf = 0, ga = 0, rp = 0, rf = 0;
        for (int qq = 0; qq < kNChg; ++qq)
        {
          if (q < kNChg && qq != q) continue;
          TH1D *a = Sum(kSigS, {CntName(kPass, qq, s)}), *b = Sum(kSigS, {CntName(kFail, qq, s)});
          TH1D *ag = Sum(kSigS, {CntGmName(kPass, qq, s)}), *bg = Sum(kSigS, {CntGmName(kFail, qq, s)});
          TH1D *an = Sum(kSigS, {CntGmAnyName(qq, s)});
          TH1D *ra = Sum(kSigS, {HName(kDiscFit, kPass, qq, s, k)}), *rb = Sum(kSigS, {HName(kDiscFit, kFail, qq, s, k)});
          if (!a || !b || !ag || !bg || !an || !ra || !rb) return 2;
          sp += a->GetBinContent(k + 1);  sf += b->GetBinContent(k + 1);
          gp += ag->GetBinContent(k + 1); gf += bg->GetBinContent(k + 1);
          ga += an->GetBinContent(k + 1);
          rp += ra->Integral(1, ra->GetNbinsX()); rf += rb->Integral(1, rb->GetNbinsX());
          delete a; delete b; delete ag; delete bg; delete an; delete ra; delete rb;
        }
        std::cout << Form("[EPS] %-6s k%d %-3s  SR %.4f   all m_T %.4f   gen-matched %.4f   eps(relIso<0.3 | gen-matched) %.4f   [%s]\n",
                          kSchemes[s].name, k, q < kNChg ? kChgName[q] : "all", (rp + rf) > 0 ? rp / (rp + rf) : 0.0,
                          (sp + sf) > 0 ? sp / (sp + sf) : 0.0, (gp + gf) > 0 ? gp / (gp + gf) : 0.0,
                          ga > 0 ? (gp + gf) / ga : 0.0, BinLabel(s, k).c_str());
      }

  // ---------------- the Z tag-and-probe record ----------------
  std::cout << "\n==== Z tag-and-probe (inclusive): legs pT > 25, |eta_SC| < 2.4, no crack, relIso < 0.3;"
               " tag = passing + EG10-matched; the OS pair closest to the Z mass ====\n";
  double cnt[kNZCat][2] = {}, cntDY[kNZCat][2] = {};
  for (int zc = 0; zc < kNZCat; ++zc)
  {
    ZChan z = MakeZ(zc);
    if (!z.ok()) return 2;
    cnt[zc][0] = z.data->Integral();   cnt[zc][1] = z.ssData->Integral();
    cntDY[zc][0] = z.zsig->Integral(); cntDY[zc][1] = z.ssDY->Integral();
    std::cout << Form("[ZTNP] %-5s OS data %4.0f   SS data %3.0f | DY MC %8.2f (SS %6.2f)   W MC %.3f   W tau %.3f   DY tau %.3f\n",
                      z.dir.c_str(), cnt[zc][0], cnt[zc][1], cntDY[zc][0], cntDY[zc][1], z.w->Integral(),
                      z.wtau->Integral(), z.ztau->Integral());
    z.Delete();
  }
  // the fit's estimator eps = 2 N_PP / (2 N_PP + N_PF) as a cut-and-count, SS-subtracted
  // in data AND in DY MC (the SS sample also holds charge-flipped DY; subtracting it
  // on both sides cancels that part). Errors: Poisson per EVENT category, as in the fit.
  auto EpsZ = [](double pp, double pf, double vpp, double vpf, double &err)
  {
    const double P = 2.0 * pp, F = pf, t = P + F;
    err = 0.0;
    if (t <= 0) return 0.0;
    err = std::sqrt(4.0 * F * F * vpp + P * P * vpf) / (t * t);
    return P / t;
  };
  double eDe = 0, eMe = 0, eM0e = 0;
  const double eD  = EpsZ(cnt[kZPP][0] - cnt[kZPP][1], cnt[kZPF][0] - cnt[kZPF][1],
                          cnt[kZPP][0] + cnt[kZPP][1], cnt[kZPF][0] + cnt[kZPF][1], eDe);
  const double eM  = EpsZ(cntDY[kZPP][0] - cntDY[kZPP][1], cntDY[kZPF][0] - cntDY[kZPF][1], 0, 0, eMe);
  const double eM0 = EpsZ(cntDY[kZPP][0], cntDY[kZPF][0], 0, 0, eM0e);
  std::cout << Form("[ZTNP] cut-and-count (2 PP / (2 PP + PF), SS-subtracted): eps_data %.4f +- %.4f, eps_DY-MC %.4f"
                    " (%.4f without SS subtraction) -> SF %.4f +- %.4f\n",
                    eD, eDe, eM, eM0, eM > 0 ? eD / eM : 0.0, eM > 0 ? eDe / eM : 0.0);
  {
    TH1D *pd = Get(0, kTnpProbeHist), *pm = Sum(kZS, {kTnpProbeHist});
    if (!pd || !pm) return 2;
    const double dP = pd->GetBinContent(1) - pd->GetBinContent(3), dF = pd->GetBinContent(2) - pd->GetBinContent(4);
    const double mP = pm->GetBinContent(1) - pm->GetBinContent(3), mF = pm->GetBinContent(2) - pm->GetBinContent(4);
    const double ed = (dP + dF) > 0 ? dP / (dP + dF) : 0.0, em = (mP + mF) > 0 ? mP / (mP + mF) : 0.0;
    std::cout << Form("[ZTNP] exact per-pair probes: data OS pass %.0f fail %.0f, SS pass %.0f fail %.0f -> eps %.4f;"
                      " DY MC -> eps %.4f; SF %.4f  (vs the 2-probes-per-PP estimator above)\n",
                      pd->GetBinContent(1), pd->GetBinContent(2), pd->GetBinContent(3), pd->GetBinContent(4), ed, em,
                      em > 0 ? ed / em : 0.0);
    delete pd; delete pm;
  }

  // ---------------- the EGM predictions: <SF_EGM> over the passing MC probes / electrons ----------------
  if (!LoadEgm())
  {
    std::cerr << "[FATAL] cannot read " << kEgmCsv << " -- run python3 ../skim/sf/extract_electron_sf.py\n";
    return 2;
  }
  auto Get2 = [&](int is, const std::string &name) -> TH2D *
  {
    TH2D *src = gF[is] ? (TH2D *)gF[is]->Get(name.c_str()) : nullptr;
    if (!src) { std::cerr << "[FATAL] " << kSampleName[is] << ": missing " << name << "\n"; return nullptr; }
    TH2D *h = (TH2D *)src->Clone(Form("%s_%s_c%d", name.c_str(), kSampleName[is], gClone++));
    h->SetDirectory(nullptr);
    if (is > 0) h->Scale(gK[is]);
    return h;
  };
  {
    double sfp = 0, sfe = 0;
    TH2D *mz = Get2(kZS[0], EgmTnpMapName(true));
    if (!mz) return 2;
    EgmPrediction(mz, -1, 0, sfp, sfe);
    std::cout << Form("[EGMPRED] Z TnP (passing DY-MC probes)             SF_EGM = %.4f +- %.4f\n", sfp, sfe);
    delete mz;
    TH2D *mw = nullptr;
    for (int is : kSigS)
      for (int q = 0; q < kNChg; ++q)
      {
        TH2D *h = Get2(is, EgmMapName(true, q));
        if (!h) return 2;
        if (!mw) mw = h;
        else { mw->Add(h); delete h; }
      }
    for (int s = 0; s < kNScheme; ++s)
      for (int k = 0; k < kSchemes[s].n; ++k)
      {
        EgmPrediction(mw, s, k, sfp, sfe);
        std::cout << Form("[EGMPRED] W %-6s k%d (passing W-MC e, m_T > 40)    SF_EGM = %.4f +- %.4f   [%s]\n",
                          kSchemes[s].name, k, sfp, sfe, BinLabel(s, k).c_str());
      }
    delete mw;
  }

  // ---------------- the W inputs: SR (pT at m_T > 40) + the in-fit ABCD control regions ----------------
  int nFiles = 0;
  const char *discName = kDiscs[kDiscFit].name;
  for (int s = 0; s < kNScheme; ++s)
  {
    std::cout << "\n==== scheme " << kSchemes[s].name << ", discriminant " << discName
              << " (in-fit ABCD: B/D at m_T < " << kCRMtMax << ", SR/C at m_T > " << kSRMtMin << ") ====\n";
    for (int iw = 0; iw < kNWin; ++iw)
    {
      const std::string base = Form("%s/combine_input_idiso_%s_%s%s", kRootDir, kSchemes[s].name, discName, kWinSuffix[iw]);
      const std::string fout = base + ".root";
      TFile fo(fout.c_str(), "RECREATE");
      if (fo.IsZombie()) { std::cerr << "[FATAL] cannot write " << fout << "\n"; return 2; }
      const bool report = (iw == 0); // the prefit plots once, for the nominal window
      for (int q = 0; q < kNChg; ++q)
        for (int k = 0; k < kSchemes[s].n; ++k)
          for (int c : {(int)kPass, (int)kFail})
          {
            WChan ch = MakeW(c, q, s, k, iw);
            if (!ch.ok()) return 2;
            TDirectory *dir = fo.mkdir(ch.dir.c_str());
            WriteHist(dir, ch.data, "data_obs", fout);
            WriteHist(dir, ch.sig, "signal", fout);
            WriteHist(dir, ch.z, "z", fout);
            WriteHist(dir, ch.ztau, "ztau", fout);
            WriteHist(dir, ch.wtau, "wtau", fout);
            WriteHist(dir, ch.qcd, "qcd", fout);
            // the three 1-bin control regions
            for (int r = 0; r < kNReg; ++r)
            {
              TDirectory *cd = fo.mkdir(Form("%s_%s", ch.dir.c_str(), kRegName[r]));
              std::vector<std::pair<const char *, double>> procs = {{"data_obs", ch.cData[r]}, {"qcd", ch.cQcd0[r]}};
              if (r == kRB)
                for (auto p : std::vector<std::pair<const char *, double>>{
                         {"signal", ch.cSig[r]}, {"wtau", ch.cWtau[r]}, {"z", ch.cZ[r]}, {"ztau", ch.cZtau[r]}})
                  procs.push_back(p);
              else
                procs.push_back({"ewk", ch.cEwk[r]});
              for (auto &p : procs)
              {
                TH1D *h = Count(p.first, p.second);
                WriteHist(cd, h, p.first, fout);
                delete h;
              }
            }
            // the record: the prefit ABCD and, in the SR, its closure at r = 1 and the W / QCD
            const double srD = ch.data->Integral(1, ch.data->GetNbinsX());
            const double srW = ch.sig->Integral(1, ch.sig->GetNbinsX()) + ch.wtau->Integral(1, ch.wtau->GetNbinsX());
            const double srDY = ch.z->Integral(1, ch.z->GetNbinsX()) + ch.ztau->Integral(1, ch.ztau->GetNbinsX());
            std::cout << Form("[ABCD] %-6s %-4s %-12s B %7.0f (QCD %8.1f) C %7.0f (QCD %8.1f) D %7.0f (QCD %8.1f)"
                              " -> A0 %8.1f | SR data %7.0f W %7.1f DY %6.1f: (data - EWK)/A0 = %.3f, W/A0 = %.3f\n",
                              kSchemes[s].name, kWinName[iw], ch.dir.c_str(), ch.cData[kRB], ch.cQcd0[kRB], ch.cData[kRC],
                              ch.cQcd0[kRC], ch.cData[kRD], ch.cQcd0[kRD], ch.A0, srD, srW, srDY,
                              ch.A0 > 0 ? (srD - srW - srDY) / ch.A0 : 0.0, ch.A0 > 0 ? srW / ch.A0 : 0.0);
            if (report)
              PrefitPlot(ch, Form("%s/prefit/%s_%s_%s", kPlotDir, kSchemes[s].name, discName, ch.dir.c_str()),
                         Form("W #rightarrow e#nu, %s %s", c == kPass ? "pass" : "fail", kChgName[q]), BinLabel(s, k));
            ch.Delete();
          }

      // the DY anchor of the W fit: Z_PP (the same in every input file)
      if (!WriteZ(fo, kZPP, fout, false)) return 2;
      fo.Close();

      // the sidecar the fork's card generator reads
      std::ofstream meta(base + "_meta.txt");
      meta << "# written by correction/idiso_sf_inputs.C -- the ID+iso SF fit input next to it\n";
      meta << "scheme " << kSchemes[s].name << "\n";
      meta << "disc " << discName << "\n";
      meta << "sbwin " << kWinName[iw] << "\n";
      meta << "nbins " << kSchemes[s].n << "\n";
      for (int k = 0; k < kSchemes[s].n; ++k)
        meta << "bin " << k << " " << kSchemes[s].lo[k] << " " << (kSchemes[s].hi[k] >= kPtInf ? -1 : kSchemes[s].hi[k]) << "\n";
      meta << "xtitle " << kDiscs[kDiscFit].axisTitle << "\n";
      meta << "ytitle Events / 2.0 GeV\n";
      meta.close();
      std::cout << "[OUT] " << fout << " (+ _meta.txt)\n";
      ++nFiles;
    }
  }
  // ---------------- the Z tag-and-probe input: its own fit ----------------
  {
    const std::string base = Form("%s/combine_input_idiso_ztnp", kRootDir);
    const std::string fout = base + ".root";
    TFile fo(fout.c_str(), "RECREATE");
    if (fo.IsZombie()) { std::cerr << "[FATAL] cannot write " << fout << "\n"; return 2; }
    for (int zc = 0; zc < kNZCat; ++zc)
      if (!WriteZ(fo, zc, fout, true)) return 2;
    fo.Close();
    std::ofstream meta(base + "_meta.txt");
    meta << "# written by correction/idiso_sf_inputs.C -- the Z tag-and-probe input next to it\n";
    meta << "scheme ztnp\n" << "disc mee\n" << "sbwin none\n" << "nbins 1\n" << "bin 0 60 120\n";
    meta << "xtitle m_{ee} (GeV)\n" << "ytitle Events / 1.0 GeV\n";
    meta.close();
    std::cout << "[OUT] " << fout << " (+ _meta.txt)\n";
    ++nFiles;
  }
  std::cout << "\n[RESULT] wrote " << nFiles << " Combine inputs to " << kRootDir << "/\n";
  for (int is = 0; is < kNSample; ++is) gF[is]->Close();
  return 0;
}
