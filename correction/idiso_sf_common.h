// correction/idiso_sf_common.h
//
// Single source for the ELECTRON ID+ISO SCALE-FACTOR CROSS-CHECK (2026-09-24).
// correction/idiso_sf_skim.C (the separate skim), idiso_sf_inputs.C (Combine
// inputs + the composition record) and idiso_sf_plots.C (efficiencies, SF, the
// EGM overlay) read the efficiency definition, the binning and every histogram
// name from here. The fit runs in the Combine fork (test/run_pO_idiso_sf.sh),
// SEPARATE from the nominal fit stream. Nothing in skim/skim.C or the nominal
// inputs is touched by this study.
//
// WHY. No electron SF is applied in the analysis. The EGM 2025Prompt SFs exist
// (skim/sf/EGM/electron.json.gz), but our ntuple IDs are heavy-ion forest
// classifiers, not EGM VID. The candidate switch is to eleMVAIdWP90 &&
// eleMVAIsoWP90 -- the closest analogue of EGM's `wp90iso` (the Run-3 MVA ID
// trained WITH isolation inputs, 90% WP: ONE SF for ID and isolation together,
// measured by EGM tag-and-probe relative to a reconstructed electron) -- and
// to apply the pp SF only if our own data show that it describes the forest
// selection. The two are not the same object (our W-MC efficiency is ~0.80 vs
// EGM effMC 0.85-0.93), so only the data/MC RATIO is compared.
//
// WHAT IS MEASURED -- two efficiencies of the same electron definition, run
// together by the fork's run_pO_idiso_sf.sh:
//  (1) a W-sample efficiency, built like the trigger SF of
//      correction/trig_eff_mb.C (user decision 2026-09-24), per coarse bin, in
//      an electron-only W + Z_PP fit (the Z_PP peak pins the DY);
//  (2) since 2026-09-25 (user: "an inclusive TnP style using Z, at the same
//      time"), an INCLUSIVE Z tag-and-probe efficiency -- see "Z tag-and-probe"
//      below -- in its OWN small fit (Z_PP + Z_PF). Sharing the DY scale with
//      the W channels let them pull the Z efficiency: in one likelihood SF_Z
//      moved 0.96 -> 1.01 across the W variants, whose Z channels are identical.
// The fail category of (1) is ~95-99% QCD, which makes its W count depend on the
// QCD model at the level of the W itself; the Z probes have a ~4% background.
// (1), the W sample:
//   total = event filter, |vz| < 15, HLT_OxyL1SingleEG10_v1 fired, DY veto, and
//           the LEADING electron (highest pT with pT > 25, |eta_SC| < 2.4,
//           outside the ECAL crack) matched to the EG10 trigger object and
//           with relIso < 0.3
//   pass  = total && that electron passes eleMVAIdWP90 && eleMVAIsoWP90
//   fail  = total && !pass
//   eps   = N_W(pass) / N_W(pass + fail); in data the W counts are the POST-FIT
//           W yields of an electron-only W+Z fit (the electron QCD is far too
//           large to count), in MC they are counted; SF = eps_data / eps_MC.
// The relIso < 0.3 cut on the total exists so that relIso 0.3-1.0 is OUTSIDE
// the measurement and can supply the data-driven QCD MET templates (without
// it the template events would be the fitted fail events themselves). It costs
// 0.8% of W electrons (W+ MC: 0.992 pass it) and 99.93% of passing electrons
// are below 0.3 anyway, so our SF differs from EGM's reco-electron definition
// only by SF(relIso < 0.3) ~ 1 (sub-percent). The skim measures that fraction
// on gen-matched W electrons (the "any" counters) so the log can quote it.
// The DY veto uses TOTAL-level legs (independent of the tested ID) over
// 60 < m < 120, and the Z-channel legs are a subset of them, so the W and Z
// channels of the fit never share an event.
//
// Other differences from EGM's definition, to be stated with the result: W
// kinematics instead of Z, and a trigger-matched electron.
//
// THE QCD of the W channels is the NOMINAL IN-FIT ABCD (QCD_MODE=abcd of the
// nominal simfit; correction/qcd_abcd.C) -- user 2026-09-25: "always the in-fit
// ABCD, so we don't assume the prefit normalization". Plane = (m_T x relIso),
// the FITTED variable is the lepton pT (never an ABCD axis), exactly like the
// nominal leppt_mt40 fit. Per W channel (charge, category, coarse bin):
//   SR  = the channel, m_T > 40, fitted in pT
//   CRB = the same category at m_T < 30         (1-bin counting channel)
//   CRC = the relIso sideband at m_T > 40       (1-bin; its pT shape = the SR
//                                                QCD template)
//   CRD = the relIso sideband at m_T < 30       (1-bin)
//   (30-40 = the nominal buffer). Free scales sB, sC, sD float the three QCD
//   counts; the SR QCD (total A0 = B0 C0 / D0 of the prefit EWK-subtracted
//   counts) is scaled by the formula sB sC / sD; the W in CRB rides the
//   channel's own r (the in-fit EWK subtraction), CRC/CRD EWK is frozen; the
//   residual lnN (kappa 1.15 = the nominal electron value) sits on the SR QCD
//   only. The sideband is 0.3-1.0 (not the nominal 0.2-1.0) because the total
//   already holds everything below 0.3; pass-like channels use the sideband
//   with eleMVAIdWP90 (sbid), fail-like ones any ID (sball).
// NB this makes eps the efficiency of W electrons with m_T > 40 (the nominal
// W signal region), in data and MC alike.
//
// (2), THE Z TAG-AND-PROBE (inclusive, charge-inclusive):
//   leg   = the total's electron without the trigger match: pT > 25,
//           |eta_SC| < 2.4, outside the crack, relIso < 0.3
//   tag   = a leg that passes AND is matched to an EG10 object (DR < 0.4)
//   pair  = the OS pair of legs with 60 < m < 120 and at least one tag,
//           closest to the Z mass; one per event
//   PP    = both legs pass -> two passing probes (each leg probes the other)
//   PF    = exactly one passes -- it is the tag -> one failing probe
//   eps   = 2 N_PP / (2 N_PP + N_PF), with N the POST-FIT DY yields of the
//           Z_PP and Z_PF channels of the Z TnP fit (DY template, own scales
//           r_ZPP / r_ZPF, + a flat background with free normalization); MC
//           counts the DY MC; SF = eps_data / eps_MC.
// A PP event with only ONE matched leg really holds a single probe; counting it
// twice biases eps by ~0.2% in data and MC alike (it cancels in the SF). The
// skim also records the exact per-pair probe counts (h_tnp_probes) and the
// same-sign twins, for the SS-subtracted cut-and-count cross-check in the log.
// The Z_PP events replace the old "both legs pass, lead 15 / sub 10" Z channel
// as the r_Z anchor of the DY in the W channels. Every Z leg is a DY-veto leg,
// so the W and Z channels still never share an event.

#ifndef PO_IDISO_SF_COMMON_H
#define PO_IDISO_SF_COMMON_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TH2D.h"
#include "TString.h"

#include "../skim/skim_common.h"

namespace pOIdIso
{

// ---------------------------------------------------------------- definition
const double kPtMin        = 25.0;   // leading electron (the W analysis cut)
const double kEtaSCMax     = 2.4;    // |eta_SC|; the ECAL crack is vetoed on eta_SC too
const double kIsoTotalMax  = 0.3;    // relIso cut of the TOTAL (see header)
const double kIsoSbLo      = 0.3;    // QCD-template sideband: relIso in [0.3, 1.0),
const double kIsoSbMid     = 0.6;    //   stored split at 0.6 for a window variation
const double kIsoSbHi      = 1.0;
const double kTrigMatchDR  = 0.4;    // skim.C's trigMatchDR
const double kVzMax        = 15.0;
const double kDyLegPtMin   = 10.0;   // DY-veto legs (total-level: relIso < 0.3, acceptance)
const double kMllLo        = 60.0;   // DY veto window == Z-channel window
const double kMllHi        = 120.0;
const double kZMass        = 91.1876; // the TnP pair closest to it

const char *const kPassIdBranch  = "eleMVAIdWP90";
const char *const kPassIsoBranch = "eleMVAIsoWP90";
const char *const kTrigPath      = "HLT_OxyL1SingleEG10_v1";
const char *const kTrigObjTree   = "hltobject/HLT_OxyL1SingleEG10_v";

// ---------------------------------------------------------------- categories
enum Cat { kPass = 0, kFail, kSbIdLo, kSbIdHi, kSbAllLo, kSbAllHi, kNCat };
const char *const kCatName[kNCat] = {"pass", "fail", "sbidlo", "sbidhi", "sballlo", "sballhi"};
// sbid*  = sideband && eleMVAIdWP90 (the pass-like QCD template)
// sball* = sideband, any ID          (the fail-like QCD template)

const int         kNChg = 2;
const char *const kChgName[kNChg] = {"Wp", "Wm"};
inline int ChgIndex(int q) { return q > 0 ? 0 : 1; }

// ---------------------------------------------------------------- binning
// Coarse bins = EGM's own edges, merged (user request: our coarse result is
// overlaid on the EGM fine-binned SF). |eta_SC| bins are EGM's eta edges folded
// in +-eta (the crack gap between 1.4442 and 1.566 belongs to no bin); pT bins
// merge EGM's 20-35 (from our 25 floor), 35-50, and everything above 50.
enum SchemeVar { kVarIncl = 0, kVarAbsEta, kVarPt };
struct Scheme
{
  const char *name;
  SchemeVar   var;
  int         n;
  double      lo[4], hi[4];
  const char *axisTitle;
};
const int    kNScheme  = 3;
const int    kMaxBins  = 4;
const double kPtInf    = 1.0e9;
const Scheme kSchemes[kNScheme] = {
    {"incl",   kVarIncl,   1, {0, 0, 0, 0},            {0, 0, 0, 0},                 "inclusive"},
    {"abseta", kVarAbsEta, 4, {0.0, 0.8, 1.566, 2.0},  {0.8, 1.4442, 2.0, 2.4},      "|#eta_{SC}|"},
    {"pt",     kVarPt,     3, {25.0, 35.0, 50.0, 0},   {35.0, 50.0, kPtInf, 0},      "p_{T} (GeV)"},
};

inline int SchemeIndex(const std::string &name)
{
  for (int s = 0; s < kNScheme; ++s)
    if (name == kSchemes[s].name) return s;
  return -1;
}

// bin of the leading electron in scheme s (-1 = in no bin, e.g. the crack gap)
inline int FindBin(int s, double scEta, double pt)
{
  const Scheme &sc = kSchemes[s];
  if (sc.var == kVarIncl) return 0;
  const double x = (sc.var == kVarAbsEta) ? std::fabs(scEta) : pt;
  for (int k = 0; k < sc.n; ++k)
    if (x >= sc.lo[k] && x < sc.hi[k]) return k;
  if (sc.var == kVarAbsEta && std::fabs(x - kEtaSCMax) < 1e-9) return sc.n - 1; // |eta| == 2.4 exactly
  return -1;
}

// ---------------------------------------------------------------- discriminants
struct Disc
{
  const char *name;
  int         nb;
  double      lo, hi;
  const char *axisTitle;
};
const int  kNDisc = 3;
const Disc kDiscs[kNDisc] = {
    {"met",        60, 0.0, 120.0, "PF MET (GeV)"},     // = the nominal h_met_* binning (diagnostic)
    {"mt",         80, 0.0, 200.0, "m_{T} (GeV)"},      // = the nominal h_mt_* binning: the ABCD m_T axis
    {"leppt_mt40", 100, 0.0, 200.0, "p_{T}^{e} (GeV)"}, // lepton pT at m_T > 40: THE FITTED variable
};
const int    kDiscFit  = 2;     // the fitted discriminant (leppt_mt40)
const int    kDiscMt   = 1;     // the ABCD axis
const double kCRMtMax  = 30.0;  // in-fit ABCD: B, D at m_T < 30   (the nominal yCut)
const double kSRMtMin  = 40.0;  //              SR, C at m_T > 40  (the nominal kSRMtCut)
inline int DiscIndex(const std::string &name)
{
  for (int d = 0; d < kNDisc; ++d)
    if (name == kDiscs[d].name) return d;
  return -1;
}

// Z tag-and-probe channels: m_ee of the chosen pair, one entry per event
const int    kZNb = 60;
const double kZLo = 60.0, kZHi = 120.0;
enum ZCat { kZPP = 0, kZPF, kNZCat };
const char *const kZCatName[kNZCat] = {"PP", "PF"};
// exact per-pair probe counts (Sum_w): bin 1 OS pass, 2 OS fail, 3 SS pass, 4 SS fail
const char *const kTnpProbeHist = "h_tnp_probes";

// EGM-cell map of the leading electron (signed eta_SC x pT) -- the input of the
// "EGM prediction for our coarse bin". Every coarse bin is a union of these
// cells, and the EGM SF is constant inside each, so nothing finer is needed.
// NB EGM's crack edge is 1.444, ours 1.4442: [0.8, 1.4442) maps to EGM's
// [0.8, 1.444) cell (the 0.0002-wide slice is not worth a separate cell).
const int    kNEgmEta = 10;
const double kEgmEta[kNEgmEta + 1] = {-2.4, -2.0, -1.566, -1.4442, -0.8, 0.0, 0.8, 1.4442, 1.566, 2.0, 2.4};
const int    kNEgmPt  = 7;
const double kEgmPt[kNEgmPt + 1]   = {25.0, 35.0, 50.0, 70.0, 100.0, 200.0, 300.0, 500.0};

// ---------------------------------------------------------------- the EGM table
// The wp90iso SF per (signed eta_SC, pT) cell, from the CSV that
// skim/sf/extract_electron_sf.py writes; err = (sfup - sfdown)/2 =
// sqrt(stat^2 + syst^2). Used by idiso_sf_inputs.C ([EGMPRED]) and
// idiso_sf_plots.C (the EGM overlay).
struct EgmCell { double etaLo, etaHi, ptLo, ptHi, sf, err; };
const char *const kEgmCsv = "../skim/sf/electron_sf_2025Prompt_wp90iso.csv"; // relative to correction/
inline std::vector<EgmCell> &EgmCells() { static std::vector<EgmCell> v; return v; }
inline bool LoadEgm()
{
  std::vector<EgmCell> &cells = EgmCells();
  if (!cells.empty()) return true;
  std::ifstream in(kEgmCsv);
  std::string line;
  while (std::getline(in, line))
  {
    if (line.empty() || line[0] == '#' || line.rfind("wp,", 0) == 0) continue;
    std::stringstream ss(line);
    std::string f[12];
    int n = 0;
    while (n < 12 && std::getline(ss, f[n], ',')) ++n;
    if (n < 8) continue;
    EgmCell c;
    c.etaLo = std::stod(f[1]); c.etaHi = std::stod(f[2]); c.ptLo = std::stod(f[3]); c.ptHi = std::stod(f[4]);
    c.sf = std::stod(f[5]);
    c.err = 0.5 * (std::stod(f[6]) - std::stod(f[7]));
    cells.push_back(c);
  }
  return !cells.empty();
}
inline const EgmCell *FindEgm(double eta, double pt)
{
  for (const auto &c : EgmCells())
    if (eta >= c.etaLo && eta < c.etaHi && pt >= c.ptLo && pt < c.ptHi) return &c;
  return nullptr;
}
// <SF_EGM> over the entries of a (signed eta_SC x pT) map of PASSING electrons,
// restricted to scheme s bin k (s < 0: everything). The EGM errors are averaged
// linearly (fully correlated -- a conservative band). Returns the weight used.
inline double EgmPrediction(const TH2D *m, int s, int k, double &sf, double &err)
{
  sf = err = 0.0;
  double wsum = 0.0;
  for (int ix = 1; ix <= m->GetNbinsX(); ++ix)
    for (int iy = 1; iy <= m->GetNbinsY() + 1; ++iy) // + the pT overflow (> 500)
    {
      const double wv = m->GetBinContent(ix, iy);
      if (wv <= 0) continue;
      const double eta = m->GetXaxis()->GetBinCenter(ix);
      const double pt  = iy <= m->GetNbinsY() ? m->GetYaxis()->GetBinCenter(iy) : 600.0;
      if (s >= 0 && FindBin(s, eta, pt) != k) continue;
      const EgmCell *c = FindEgm(eta, pt);
      if (!c) { std::cout << Form("[WARN] no EGM cell at eta %.3f, pT %.0f\n", eta, pt); continue; }
      sf += wv * c->sf; err += wv * c->err; wsum += wv;
    }
  if (wsum > 0) { sf /= wsum; err /= wsum; }
  return wsum;
}

// ---------------------------------------------------------------- names
inline std::string HName(int d, int c, int q, int s, int k)
{
  return Form("h_%s_%s_%s_%s_k%d", kDiscs[d].name, kCatName[c], kChgName[q], kSchemes[s].name, k);
}
// Sum_w per coarse bin (TH1D with n bins); gen-matched twin for MC
inline std::string CntName(int c, int q, int s)   { return Form("h_cnt_%s_%s_%s",   kCatName[c], kChgName[q], kSchemes[s].name); }
inline std::string CntGmName(int c, int q, int s) { return Form("h_cntgm_%s_%s_%s", kCatName[c], kChgName[q], kSchemes[s].name); }
// gen-matched leading electrons with NO relIso requirement (-> eps(relIso < 0.3))
inline std::string CntGmAnyName(int q, int s)     { return Form("h_cntgm_any_%s_%s", kChgName[q], kSchemes[s].name); }
// EGM-cell maps of the total and of pass
inline std::string EgmMapName(bool pass, int q)   { return Form("h2_egm_%s_%s", pass ? "pass" : "total", kChgName[q]); }
// Z TnP: the m_ee of a category, OS (fitted) or SS (the cut-and-count background)
inline std::string ZHistName(int zc, bool ss)     { return Form("h_mZ%s_%s", ss ? "ss" : "", kZCatName[zc]); }
// EGM-cell maps of the OS probes: every probe, and the passing ones
inline std::string EgmTnpMapName(bool pass)       { return Form("h2_egm_tnp_%s", pass ? "pass" : "total"); }

// Combine-input directory of one W channel, and of a Z TnP category
inline std::string WDir(int c, int q, int k) { return Form("%s_%s_k%d", kChgName[q], kCatName[c], k); }
inline std::string ZDir(int zc)              { return Form("Z_%s", kZCatName[zc]); }

// ---------------------------------------------------------------- samples
const int         kNSample = 7;
const char *const kSampleName[kNSample] = {"Data", "Wp", "Wm", "DY", "DYtau", "Wptau", "Wmtau"};
const SampleType  kSampleType[kNSample] = {kData, kWp, kWm, kDY, kDYtau, kWptau, kWmtau};
// pONorm::MCScale labels (skim/mc_norm.h); "" for data
const char *const kSampleNormLabel[kNSample] = {"", "Wp_ele", "Wm_ele", "DYee", "DYtau", "Wp_tau", "Wm_tau"};
inline int SampleIndex(const std::string &name)
{
  for (int i = 0; i < kNSample; ++i)
    if (name == kSampleName[i]) return i;
  return -1;
}

const char *const kRootDir = "rootfile/idiso_sf_ele";   // relative to correction/
const char *const kPlotDir = "plots/idiso_sf_ele";
inline std::string SkimFile(const std::string &sample) { return std::string(kRootDir) + "/skim_" + sample + ".root"; }

// ---------------------------------------------------------------- gen match
// Is the selected reco electron a prompt W electron? Reads the EventTree's OWN
// gen block (mcPID/mcStatus/mcPt/mcEta/mcPhi/mcMomPID/mcGMomPID) -- NOT
// HiGenParticleAna/hi, whose motherIdx is -999 for every W lepton. Charge-blind,
// FSR chains (l <- l <- W) accepted via mcGMomPID, DR < 0.5, |dpT|/pT < 0.5.
// Copied from correction/trig_eff_mb.C (MatchGenLeptonFromW, anonymous
// namespace there); used only for the gen-matched eps_MC cross-check here --
// the fit templates are NOT gen-matched, exactly like the nominal ones.
const double kGenMatchDR  = 0.5;
const double kGenMatchDPt = 0.5;
inline int MatchGenLeptonFromW(double pt, double eta, double phi, int flavPdg,
                               const std::vector<int> *gPid, const std::vector<int> *gStatus,
                               const std::vector<float> *gPt, const std::vector<float> *gEta,
                               const std::vector<float> *gPhi, const std::vector<int> *gMom,
                               const std::vector<int> *gGMom)
{
  if (!gPid || !gStatus || !gPt || !gEta || !gPhi || !gMom) return -1;
  const size_t ng = std::min({gPid->size(), gStatus->size(), gPt->size(), gEta->size(),
                              gPhi->size(), gMom->size()});
  int    best   = -1;
  double drBest = 1e9;
  for (size_t i = 0; i < ng; ++i)
  {
    if (std::abs(gPid->at(i)) != flavPdg) continue;
    if (gStatus->at(i) != 1) continue;
    if (gPt->at(i) <= 0) continue;
    const bool momW = std::abs(gMom->at(i)) == 24;
    const bool fsrW = (std::abs(gMom->at(i)) == flavPdg) && gGMom && i < gGMom->size() &&
                      std::abs(gGMom->at(i)) == 24;
    if (!momW && !fsrW) continue;
    const double d = pOSkim::DeltaR(gEta->at(i), gPhi->at(i), eta, phi);
    if (d >= kGenMatchDR) continue;
    if (std::fabs(gPt->at(i) - pt) / gPt->at(i) >= kGenMatchDPt) continue;
    if (d < drBest) { drBest = d; best = (int)i; }
  }
  return best;
}

} // namespace pOIdIso

#endif
