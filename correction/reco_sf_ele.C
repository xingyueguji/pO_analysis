// correction/reco_sf_ele.C
//
// ELECTRON RECO EFFICIENCY IN pO -- a Z->ee tag-and-probe of the FIRST link of
// the EGM factorization (2026-09-29). THE pO CROSS-CHECK of the applied
// electron Reco SF: since 2026-10-01 (user decision, as on 2026-09-29) the
// analysis applies the pp EGM Reco SF (Reco20to75 / RecoAbove75; <RECO> =
// 0.929 per W electron, 0.865 per Z->ee pair), which this macro keeps as its
// comparison (RecoEGM). On 2026-09-30 only, its inclusive result WAS the
// applied factor: skim/electron_sf.h reads recosf_noid_incl0 from
// rootfile/reco_sf_ele/reco_sf_ele_result.root and applies that one number to
// every electron under kEleRecoSource = kRecoTnP (one constant switches back).
// The chain, each SF relative to the step before it:
//
//   supercluster --[Reco]--> reconstructed electron --[wp90iso]--> WP90 ID+iso
//                --[our MB trigger SF]--> fired and matched
//
// so the Reco probe carries NO ID: it only has to be a supercluster. The
// ID+iso link was checked on 2026-09-25 (correction/idiso_sf_*: Z TnP SF
// 0.959 +- 0.025 vs EGM 0.956). This macro checks the Reco link.
//
// WHY. With the pp Reco SF, Z->ee data/MC is 1.31 against Z->mumu 1.07
// (sigma_Z e/mu = 1.27 +- 0.11 in the 09-29 per-flavour fits). The pp 2025
// Reco SF (data 0.89-0.93 vs MC 0.97-0.98, SF ~0.93) is the one electron factor
// never checked in pO. A scratch feasibility pass (2026-09-29) found ~550 probes
// in 80-100 GeV and a first estimate SF ~0.94-0.98.
//
// DEFINITIONS (single source: namespace pORecoSF below)
//   event = pO event filter, |vz| < 15, HLT_OxyL1SingleEG10_v1 fired
//   tag   = the tag of the ID+iso tag-and-probe (correction/idiso_sf_common.h):
//           pT > 25, |eta_SC| < 2.4 outside the ECAL crack, relIso < 0.3,
//           eleMVAIdWP90 && eleMVAIsoWP90, matched to an EG10 object (DR < 0.4)
//   probe = a PHOTON object, i.e. its supercluster: E_T^SC > 20 GeV (the lower
//           edge of EGM Reco20to75), |eta_SC| < 2.4 outside the crack,
//           DR_SC > 0.3 from the tag's SC, one probe per SC (photons within
//           DR_SC < 0.01 of an earlier one are dropped). The forest's photon
//           collection starts at E_T^SC ~ 5 GeV, and every electron SC also
//           appears as a photon (100.00% of DY-MC tags, 99.99% in data), so the
//           photons are a complete SC list above 20 GeV.
//           NO ID ON THE PROBE (nominal): H/E is part of the electron seeding,
//           so an H/E-like cut would remove precisely the SCs that fail the
//           reconstruction and bias eps high -- e.g. if the pO underlying event
//           raises H/E more in data than in MC. A loose ECAL-only probe ID
//           (H/E < 0.1, sigma_ietaieta < 0.012 EB / 0.035 EE) is filled as a
//           CROSS-CHECK variant only (it cuts the fake-photon background).
//   pass  = any forest electron on the probe's SC (DR_SC < 0.05), no ID. That
//           is the collection the analysis selects from, so the SF applies to
//           exactly what the analysis uses.
//   pair  = every (tag, probe) with 60 <= m < 120; m from the tag electron and
//           the probe PHOTON, for pass and fail alike (same m definition).
//   eps   = the pass fraction of the Z signal. DATA: a simultaneous binned ML
//           fit of the pass and fail m spectra: signal = the DY-MC pass / fail
//           templates convolved with ONE Gaussian (mass scale + resolution,
//           shared by both categories -- the Z->ee data peak sits ~0.45 GeV
//           above MC), background = one exponential per category, eps =
//           N_sig(pass) / N_sig. MC: counted in the same m range from the same
//           DY-MC histograms (the fit templates), so eps_data and eps_MC share
//           one definition, FSR-photon probes included in both.
//   SF    = eps_data / eps_MC, compared with the EGM Reco SF averaged over the
//           same probes (the passing DY-MC probes, at the matched electron's
//           eta_SC and pT) -- the "[EGMPRED]" of the ID+iso study.
// BINS: inclusive; barrel / endcap (|eta_SC| < 1.4442 / > 1.566); E_T^SC 20-45 /
// > 45 GeV (EGM's pT edges). Per-bin fits fix the Gaussian to the inclusive fit.
// SYSTEMATIC (nominal variant, every bin) = the largest |eps shift| among:
// background Chebychev(1) instead of the exponential, fit range 70-110 instead
// of 60-120, the Gaussian floated per bin (bins) / fixed at MC (inclusive) --
// counting only variations whose fit CONVERGED (status 0): in the small bins a
// narrower range lets the background swallow every fail and drives eps to the
// eps = 1 boundary, which is a failed fit, not a shift. A nominal fit that does
// not converge is flagged "[WARN] ... do not quote this bin".
// CROSS-CHECK "[CC]": cut-and-count in 80-100 GeV, raw pass over raw pass +
// (raw fail - the fitted fail background inside the window), MC in the same
// window.
//
// RUN: ./run_reco_sf.sh [skim|fit|all]  (logs: logs/reco_sf_ele_*.log = THE record)
//   root -l -b -q 'reco_sf_ele.C+("skim","Data")'   -> rootfile/reco_sf_ele/skim_Data.root
//   root -l -b -q 'reco_sf_ele.C+("skim","DY")'     -> rootfile/reco_sf_ele/skim_DY.root
//   root -l -b -q 'reco_sf_ele.C+("fit")'           -> plots/reco_sf_ele/fit_<variant>_<bin>_{pass,fail}
//                                                      + reco_sf_vs_{eta,et} + reco_sf_ele.csv
//                                                      + rootfile/reco_sf_ele/reco_sf_ele_result.root
//                                                        (labeled TH1D recosf_<variant>_<scheme><bin>;
//                                                        read by the electron MC skim under
//                                                        kEleRecoSource = kRecoTnP only -- then re-skim
//                                                        the electron MC after re-running the fit)
// PLOTS use the repo's standard functions (plotting_helper.C): the fits as
// SaveNicePlot1D_WithBkg (data over the stacked fitted components, pull pad),
// the SF vs |eta_SC| / E_T^SC as SaveNiceGraph_Overlay (stat bars + syst boxes).
// Run from correction/ (the EGM tables are found at ../skim/sf/).

#include "../skim/skim_common.h"
#include "../skim/electron_sf.h"
#include "../plotting/plotting_helper.C"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector2.h"

#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooCategory.h"
#include "RooChebychev.h"
#include "RooDataHist.h"
#include "RooExponential.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooMsgService.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace pORecoSF
{
// ---------------------------------------------------------------- definition
const double kTagPtMin    = 25.0;
const double kEtaSCMax    = 2.4;
const double kTrigMatchDR = 0.4;  // skim.C's trigMatchDR
const double kVzMax       = 15.0;
const double kProbeEtMin  = 20.0; // E_T^SC, the lower edge of EGM Reco20to75
const double kProbeTagDR  = 0.3;  // probe SC not within 0.3 of the tag's SC
const double kDupDR       = 0.01; // photons sharing one SC
const double kPassDR      = 0.05; // a forest electron on the probe's SC
const double kGenDR       = 0.1;  // MC: a prompt Z electron under the probe photon
const double kMLo = 60.0, kMHi = 120.0;
const int    kNM  = 60;           // 1 GeV bins
const double kSysLo = 70.0, kSysHi = 110.0;
const char *const kTrigPath    = "HLT_OxyL1SingleEG10_v1";
const char *const kTrigObjTree = "hltobject/HLT_OxyL1SingleEG10_v";
const char *const kOutDir      = "rootfile/reco_sf_ele";
const char *const kPlotDir     = "plots/reco_sf_ele";

// probe variants
enum Var { kNoId = 0, kLooseId, kNVar };
const char *const kVarName[kNVar]  = {"noid", "looseid"};
const char *const kVarLabel[kNVar] = {"probe: any SC (nominal)", "probe: loose ECAL ID (cross-check)"};
inline bool LooseProbeId(double hoe, double sieie, double absEtaSC)
{
  return hoe < 0.10 && sieie < (absEtaSC < 1.4442 ? 0.012 : 0.035);
}

// binning schemes
struct Scheme
{
  const char *name;
  int         n;
  double      lo[2], hi[2];
  bool        byEta;
  const char *label[2];
};
const int    kNScheme = 3;
const Scheme kSchemes[kNScheme] = {
    {"incl", 1, {0.0, 0.0}, {0.0, 0.0}, false, {"inclusive", ""}},
    {"eta", 2, {0.0, 1.566}, {1.4442, 2.4}, true, {"barrel", "endcap"}},
    {"et", 2, {20.0, 45.0}, {45.0, 1.0e9}, false, {"E_{T}^{SC} 20-45 GeV", "E_{T}^{SC} > 45 GeV"}},
};
inline int RecoBinOf(int s, double absEtaSC, double et)
{
  const Scheme &sc = kSchemes[s];
  if (sc.n == 1) return 0;
  const double x = sc.byEta ? absEtaSC : et;
  for (int k = 0; k < sc.n; ++k)
    if (x >= sc.lo[k] && x < sc.hi[k]) return k;
  return -1;
}
inline TString RecoHName(int v, int s, int k, bool pass, const char *sfx = "")
{
  return Form("h_m_%s_%s%d_%s%s", kVarName[v], kSchemes[s].name, k, pass ? "pass" : "fail", sfx);
}
inline double RecoDRsc(double e1, double p1, double e2, double p2)
{
  const double dp = TVector2::Phi_mpi_pi(p1 - p2);
  return std::sqrt((e1 - e2) * (e1 - e2) + dp * dp);
}
} // namespace pORecoSF

// ============================================================================
// SKIM: one ntuple (Data or DY->ee MC) -> m(tag, probe) histograms per
// variant x scheme x bin x {pass, fail}
// ============================================================================
int RecoSfSkim(const std::string &sampleName, Long64_t nmax)
{
  using namespace pORecoSF;
  const bool isMC = (sampleName != "Data");
  std::string fname;
  if (!isMC)
    fname = pOSkim::kDefaultDataFile;
  else if (sampleName == "DY")
    fname = pOSkim::ResolveMCSample(kDY, "ele").fname; // SampleType is global (skim_common.h)
  else
  {
    std::cerr << "[FATAL] sample must be Data or DY, got " << sampleName << "\n";
    return 2;
  }
  std::cout << "[INPUT] " << sampleName << " <- " << fname << "\n";
  std::cout << "[CONFIG] tag: WP90 ID+iso, pT > " << kTagPtMin << ", |eta_SC| < " << kEtaSCMax
            << " (no crack), relIso < " << pOSkim::kEleIsoRelMax << ", EG10-matched DR < " << kTrigMatchDR
            << " | probe: photon SC, E_T^SC > " << kProbeEtMin << ", |eta_SC| < " << kEtaSCMax
            << " (no crack), DR_SC(tag) > " << kProbeTagDR << ", NO ID (nominal) | pass: a forest electron with DR_SC < "
            << kPassDR << " | m in [" << kMLo << ", " << kMHi << ")\n";

  TFile *f = TFile::Open(fname.c_str());
  if (!f || f->IsZombie()) { std::cerr << "[FATAL] cannot open " << fname << "\n"; return 2; }
  TTree *tEle = (TTree *)f->Get("ggHiNtuplizer/EventTree");
  TTree *tHi  = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  TTree *tHLT = (TTree *)f->Get("hltanalysis/HltTree");
  TTree *tObj = (TTree *)f->Get(kTrigObjTree);
  TTree *tEvt = (TTree *)f->Get("skimanalysis/HltTree");
  if (!tEle || !tHi || !tHLT || !tObj || !tEvt) { std::cerr << "[FATAL] missing a required tree in " << fname << "\n"; return 2; }
  for (TTree *t : {tHi, tHLT, tObj, tEvt})
    if (t->GetEntries() != tEle->GetEntries())
    {
      std::cerr << "[FATAL] tree " << t->GetName() << " (" << t->GetEntries() << ") and EventTree ("
                << tEle->GetEntries() << ") entry counts differ\n";
      return 2;
    }

  // -------- EventTree: every branch read is enabled AND addressed --------
  Int_t nEle = 0, nPho = 0;
  std::vector<float> *elePt = nullptr, *eleEta = nullptr, *elePhi = nullptr, *eleSCEta = nullptr, *eleSCPhi = nullptr;
  std::vector<int>   *eleId = nullptr, *eleIso = nullptr;
  std::vector<float> *chIso = nullptr, *neuIso = nullptr, *phIso = nullptr, *puIso = nullptr;
  std::vector<float> *phoEt = nullptr, *phoEta = nullptr, *phoPhi = nullptr, *phoSCE = nullptr, *phoSCEta = nullptr,
                     *phoSCPhi = nullptr, *phoHoE = nullptr, *phoSieie = nullptr;
  std::vector<int>   *mcPID = nullptr, *mcStatus = nullptr, *mcMomPID = nullptr;
  std::vector<float> *mcPt = nullptr, *mcEta = nullptr, *mcPhi = nullptr;
  std::vector<const char *> need = {"nEle", "elePt", "eleEta", "elePhi", "eleSCEta", "eleSCPhi",
                                    pOSkim::kEleIDBranch, pOSkim::kEleIsoBranch, "elePFChIso", "elePFNeuIso",
                                    "elePFPhoIso", "elePFPUIso", "nPho", "phoEt", "phoEta", "phoPhi", "phoSCE",
                                    "phoSCEta", "phoSCPhi", "phoHoverE", "phoSigmaIEtaIEta"};
  if (isMC)
    for (const char *b : {"mcPID", "mcStatus", "mcMomPID", "mcPt", "mcEta", "mcPhi"}) need.push_back(b);
  for (const char *b : need)
    if (!pOSkim::HasBranch(tEle, b)) { std::cerr << "[FATAL] missing EventTree branch " << b << "\n"; return 2; }
  tEle->SetBranchStatus("*", 0);
  auto bind = [&](const char *b, auto *addr) { tEle->SetBranchStatus(b, 1); tEle->SetBranchAddress(b, addr); };
  bind("nEle", &nEle);             bind("elePt", &elePt);          bind("eleEta", &eleEta);
  bind("elePhi", &elePhi);         bind("eleSCEta", &eleSCEta);    bind("eleSCPhi", &eleSCPhi);
  bind(pOSkim::kEleIDBranch, &eleId);  bind(pOSkim::kEleIsoBranch, &eleIso);
  bind("elePFChIso", &chIso);      bind("elePFNeuIso", &neuIso);   bind("elePFPhoIso", &phIso);
  bind("elePFPUIso", &puIso);
  bind("nPho", &nPho);             bind("phoEt", &phoEt);          bind("phoEta", &phoEta);
  bind("phoPhi", &phoPhi);         bind("phoSCE", &phoSCE);        bind("phoSCEta", &phoSCEta);
  bind("phoSCPhi", &phoSCPhi);     bind("phoHoverE", &phoHoE);     bind("phoSigmaIEtaIEta", &phoSieie);
  if (isMC)
  {
    bind("mcPID", &mcPID);   bind("mcStatus", &mcStatus); bind("mcMomPID", &mcMomPID);
    bind("mcPt", &mcPt);     bind("mcEta", &mcEta);       bind("mcPhi", &mcPhi);
  }

  // -------- trigger bit (exact name: FindBranchContaining may return a _Prescale twin) --------
  const std::string hltName = pOSkim::HasBranch(tHLT, kTrigPath) ? std::string(kTrigPath)
                                                                 : pOSkim::FindBranchContaining(tHLT, kTrigPath);
  if (hltName.empty()) { std::cerr << "[FATAL] hltanalysis/HltTree lacks " << kTrigPath << "\n"; return 2; }
  Int_t hltBit = 0;
  tHLT->SetBranchStatus("*", 0);
  tHLT->SetBranchStatus(hltName.c_str(), 1); tHLT->SetBranchAddress(hltName.c_str(), &hltBit);

  // -------- trigger objects --------
  std::vector<double> *toPt = nullptr, *toEta = nullptr, *toPhi = nullptr;
  for (const char *b : {"pt", "eta", "phi"})
    if (!pOSkim::HasBranch(tObj, b)) { std::cerr << "[FATAL] no " << b << " in " << kTrigObjTree << "\n"; return 2; }
  tObj->SetBranchStatus("*", 0);
  tObj->SetBranchStatus("pt", 1);  tObj->SetBranchAddress("pt", &toPt);
  tObj->SetBranchStatus("eta", 1); tObj->SetBranchAddress("eta", &toEta);
  tObj->SetBranchStatus("phi", 1); tObj->SetBranchAddress("phi", &toPhi);

  // -------- event filters, vz, gen weight --------
  Int_t ppv = 1, pcc = 1;
  const bool has_ppv = pOSkim::HasBranch(tEvt, "pprimaryVertexFilter");
  const bool has_pcc = pOSkim::HasBranch(tEvt, "pclusterCompatibilityFilter");
  tEvt->SetBranchStatus("*", 0);
  if (has_ppv) { tEvt->SetBranchStatus("pprimaryVertexFilter", 1);        tEvt->SetBranchAddress("pprimaryVertexFilter", &ppv); }
  if (has_pcc) { tEvt->SetBranchStatus("pclusterCompatibilityFilter", 1); tEvt->SetBranchAddress("pclusterCompatibilityFilter", &pcc); }
  Float_t vz = 0.f, genWeight = 1.f;
  if (!pOSkim::HasBranch(tHi, "vz")) { std::cerr << "[FATAL] HiTree lacks vz\n"; return 2; }
  tHi->SetBranchStatus("*", 0);
  tHi->SetBranchStatus("vz", 1); tHi->SetBranchAddress("vz", &vz);
  const bool has_w = isMC && pOSkim::HasBranch(tHi, "weight");
  if (has_w) { tHi->SetBranchStatus("weight", 1); tHi->SetBranchAddress("weight", &genWeight); }
  else if (isMC) std::cout << "[WARN] no HiTree 'weight' -- MC filled unweighted.\n";

  // -------- the pp EGM Reco SF, for the prediction on the same probes (MC only) --------
  // Load(.., false): everything but this macro's OWN result (reco_sf_ele_result.root,
  // which the skim applies since 2026-09-30) -- it must not need what it produces.
  pOSF::ElectronSF egm;
  if (isMC && !egm.Load(std::cout, false)) { std::cerr << "[FATAL] cannot load the EGM electron SF tables\n"; return 2; }

  // -------- output --------
  gSystem->mkdir(kOutDir, kTRUE);
  const TString outName = Form("%s/skim_%s.root", kOutDir, sampleName.c_str());
  TFile *fo = TFile::Open(outName, "RECREATE");
  std::map<TString, TH1D *> h;
  auto book = [&](const TString &n)
  {
    TH1D *x = new TH1D(n, ";m(tag, probe) [GeV];pairs / GeV", kNM, kMLo, kMHi);
    x->Sumw2();
    h[n] = x;
  };
  for (int v = 0; v < kNVar; ++v)
    for (int s = 0; s < kNScheme; ++s)
      for (int k = 0; k < kSchemes[s].n; ++k)
        for (int p = 0; p < 2; ++p) book(RecoHName(v, s, k, p));
  if (isMC)
    for (int p = 0; p < 2; ++p) book(RecoHName(kNoId, 0, 0, p, "_gm"));
  // EGM Reco SF summed over the passing probes, one bin per bin of the scheme
  TH1D *hEgmSum[kNVar][kNScheme], *hEgmW[kNVar][kNScheme];
  for (int v = 0; v < kNVar; ++v)
    for (int s = 0; s < kNScheme; ++s)
    {
      const int n = kSchemes[s].n;
      hEgmSum[v][s] = new TH1D(Form("h_egmsum_%s_%s", kVarName[v], kSchemes[s].name), "sum w*SF_EGM (pass)", n, 0, n);
      hEgmW[v][s]   = new TH1D(Form("h_egmw_%s_%s", kVarName[v], kSchemes[s].name), "sum w (pass)", n, 0, n);
    }
  const char *cntLab[] = {"events", "filter + vz", "trigger", "events with a tag", "tags",
                          "pairs 60-120", "pairs 80-100", "MC: probes on a Z electron"};
  TH1D *hCnt = new TH1D("h_count", "counters (weighted)", 8, 0, 8);
  TH1D *hCntRaw = new TH1D("h_count_raw", "counters (raw)", 8, 0, 8);
  for (int i = 0; i < 8; ++i)
  {
    hCnt->GetXaxis()->SetBinLabel(i + 1, cntLab[i]);
    hCntRaw->GetXaxis()->SetBinLabel(i + 1, cntLab[i]);
  }
  auto count = [&](int i, double w) { hCnt->Fill(i + 0.5, w); hCntRaw->Fill(i + 0.5); };

  // -------- loop --------
  const Long64_t nEntries = tEle->GetEntries();
  const Long64_t nRun     = (nmax > 0 && nmax < nEntries) ? nmax : nEntries;
  std::cout << "Processing entries: " << nRun << " of " << nEntries << "\n";
  bool warnedFilters = false, warnedTrig = false;
  Long64_t nMultiPair = 0;
  for (Long64_t ie = 0; ie < nRun; ++ie)
  {
    if (ie % 1000000 == 0) std::cout << "  event " << ie << "/" << nRun << "\n";
    tEle->GetEntry(ie);
    tHi->GetEntry(ie);
    tEvt->GetEntry(ie);
    tHLT->GetEntry(ie);
    const double w = has_w ? (double)genWeight : 1.0;
    count(0, w);
    if (!elePt || !eleEta || !elePhi || !eleSCEta || !eleSCPhi || !eleId || !eleIso || !phoEt || !phoEta ||
        !phoPhi || !phoSCE || !phoSCEta || !phoSCPhi || !phoHoE || !phoSieie)
      continue;
    if (!pOSkim::PassEventSelection_pO(warnedFilters, has_ppv, ppv, has_pcc, pcc)) continue;
    if (std::fabs(vz) > kVzMax) continue;
    count(1, w);
    if (!pOSkim::TriggerFired(hltBit)) continue;
    count(2, w);

    const int ne = (int)elePt->size();
    tObj->GetEntry(ie);
    std::vector<int> tags;
    for (int i = 0; i < ne; ++i)
    {
      if (elePt->at(i) < kTagPtMin) continue;
      const double se = eleSCEta->at(i);
      if (std::fabs(se) >= kEtaSCMax || pOSkim::InEcalGap(se)) continue;
      if (!eleId->at(i) || !eleIso->at(i)) continue;
      if (pOSkim::RelIsoPF(i, elePt, chIso, neuIso, phIso, puIso) >= pOSkim::kEleIsoRelMax) continue;
      if (!pOSkim::PassLeadingLeptonTrigMatch(kTrigMatchDR, i, eleEta, elePhi, true, toPt, true, toEta, true, toPhi,
                                              warnedTrig, "electron"))
        continue;
      tags.push_back(i);
    }
    if (tags.empty()) continue;
    count(3, w);

    const int np = (int)phoEt->size();
    int pairsInEvent = 0;
    for (int it : tags)
    {
      count(4, w);
      TLorentzVector vt;
      vt.SetPtEtaPhiM(elePt->at(it), eleEta->at(it), elePhi->at(it), pOSkim::ELE_MASS);
      std::vector<int> used;
      for (int k = 0; k < np; ++k)
      {
        const double se = phoSCEta->at(k), sp = phoSCPhi->at(k), ase = std::fabs(se);
        const double et = phoSCE->at(k) / std::cosh(se);
        if (et < kProbeEtMin || ase >= kEtaSCMax || pOSkim::InEcalGap(se)) continue;
        if (RecoDRsc(se, sp, eleSCEta->at(it), eleSCPhi->at(it)) < kProbeTagDR) continue;
        bool dup = false;
        for (int u : used)
          if (RecoDRsc(se, sp, phoSCEta->at(u), phoSCPhi->at(u)) < kDupDR) { dup = true; break; }
        if (dup) continue;
        used.push_back(k);
        TLorentzVector vp;
        vp.SetPtEtaPhiM(phoEt->at(k), phoEta->at(k), phoPhi->at(k), 0.0);
        const double m = (vt + vp).M();
        if (m < kMLo || m >= kMHi) continue;
        count(5, w);
        if (m >= 80.0 && m < 100.0) count(6, w);
        ++pairsInEvent;
        int jPass = -1;
        for (int j = 0; j < ne && jPass < 0; ++j)
          if (j != it && RecoDRsc(eleSCEta->at(j), eleSCPhi->at(j), se, sp) < kPassDR) jPass = j;
        const bool   pass  = (jPass >= 0);
        const bool   loose = LooseProbeId(phoHoE->at(k), phoSieie->at(k), ase);
        const double sfEgm = (isMC && pass) ? egm.RecoEGM(eleSCEta->at(jPass), elePt->at(jPass), 0) : 0.0;
        for (int v = 0; v < kNVar; ++v)
        {
          if (v == kLooseId && !loose) continue;
          for (int s = 0; s < kNScheme; ++s)
          {
            const int b = RecoBinOf(s, ase, et);
            if (b < 0) continue;
            h[RecoHName(v, s, b, pass)]->Fill(m, w);
            if (isMC && pass)
            {
              hEgmSum[v][s]->Fill(b + 0.5, w * sfEgm);
              hEgmW[v][s]->Fill(b + 0.5, w);
            }
          }
        }
        if (isMC)
        { // is the probe photon on a prompt Z electron (mother 23, or an FSR e -> e chain)?
          bool gm = false;
          for (size_t g = 0; g < mcPID->size() && !gm; ++g)
          {
            if (std::abs(mcPID->at(g)) != 11 || mcStatus->at(g) != 1 || mcPt->at(g) < 10.0) continue;
            const int mom = std::abs(mcMomPID->at(g));
            if (mom != 23 && mom != 11) continue;
            if (RecoDRsc(mcEta->at(g), mcPhi->at(g), phoEta->at(k), phoPhi->at(k)) < kGenDR) gm = true;
          }
          if (gm)
          {
            h[RecoHName(kNoId, 0, 0, pass, "_gm")]->Fill(m, w);
            count(7, w);
          }
        }
      }
    }
    if (pairsInEvent > 1) ++nMultiPair;
  }

  // -------- report --------
  std::cout << "[COUNT] " << sampleName << " (raw | weighted):\n";
  for (int i = 0; i < 8; ++i)
    std::cout << Form("  %-28s %12.0f | %14.1f\n", cntLab[i], hCntRaw->GetBinContent(i + 1), hCnt->GetBinContent(i + 1));
  std::cout << "  events with more than one pair in 60-120: " << nMultiPair << "\n";
  for (int v = 0; v < kNVar; ++v)
  {
    const double p = h[RecoHName(v, 0, 0, true)]->Integral(), fl = h[RecoHName(v, 0, 0, false)]->Integral();
    TH1D *hp = h[RecoHName(v, 0, 0, true)], *hf = h[RecoHName(v, 0, 0, false)];
    const int b1 = hp->FindBin(80.5), b2 = hp->FindBin(99.5);
    const double pw = hp->Integral(b1, b2), fw = hf->Integral(b1, b2);
    std::cout << Form("[RAW] %-7s inclusive: 60-120 pass %.1f fail %.1f -> %.4f | 80-100 pass %.1f fail %.1f -> %.4f\n",
                      kVarName[v], p, fl, p / std::max(1e-9, p + fl), pw, fw, pw / std::max(1e-9, pw + fw));
  }
  if (isMC)
  {
    const double p = h[RecoHName(kNoId, 0, 0, true, "_gm")]->Integral(), fl = h[RecoHName(kNoId, 0, 0, false, "_gm")]->Integral();
    std::cout << Form("[RAW] noid, probe on a Z electron (MC truth): pass %.1f fail %.1f -> %.4f\n", p, fl,
                      p / std::max(1e-9, p + fl));
    egm.PrintCounters(std::cout, "reco_sf_ele (EGM prediction lookups)");
  }
  fo->Write();
  fo->Close();
  std::cout << "[OUT] " << outName << "\n";
  return 0;
}

// ============================================================================
// FIT: eps_data from the pass/fail m spectra, eps_MC counted, SF vs EGM
// ============================================================================
namespace pORecoSF
{
struct FitOut
{
  bool   ok = false;
  int    status = -1, covQual = -1;
  double eps = 0, lo = 0, hi = 0, nSig = 0, nBP = 0, nBF = 0, mu = 0, sg = 0;
  double bFwin = 0; // fitted fail background inside 80-100 GeV (the cut-and-count cross-check)
};

// the [lo, hi) part of a 1-GeV-binned histogram (lo, hi on bin edges)
TH1D *RecoSubHist(const TH1 *h, double lo, double hi, const char *name)
{
  const int b1 = h->GetXaxis()->FindBin(lo + 1e-6), b2 = h->GetXaxis()->FindBin(hi - 1e-6);
  TH1D *o = new TH1D(name, "", b2 - b1 + 1, h->GetXaxis()->GetBinLowEdge(b1), h->GetXaxis()->GetBinUpEdge(b2));
  o->SetDirectory(nullptr);
  o->Sumw2();
  for (int b = b1; b <= b2; ++b)
  {
    o->SetBinContent(b - b1 + 1, h->GetBinContent(b));
    o->SetBinError(b - b1 + 1, h->GetBinError(b));
  }
  return o;
}

// bkg: 0 = exponential, 1 = Chebychev(1). plot = "" -> no plot.
FitOut RecoFitEps(const TH1 *dP0, const TH1 *dF0, const TH1 *tP0, const TH1 *tF0, double lo, double hi,
                  bool floatConv, double mu0, double sg0, int bkg, const TString &plot, const TString &title)
{
  FitOut out;
  std::unique_ptr<TH1D> dP(RecoSubHist(dP0, lo, hi, "fit_dP")), dF(RecoSubHist(dF0, lo, hi, "fit_dF"));
  std::unique_ptr<TH1D> tP(RecoSubHist(tP0, lo, hi, "fit_tP")), tF(RecoSubHist(tF0, lo, hi, "fit_tF"));
  RooRealVar m("m", "m(tag, probe) [GeV]", lo, hi);
  m.setBins(dP->GetNbinsX());
  m.setBins(2000, "cache");
  RooCategory cat("cat", "cat");
  cat.defineType("pass");
  cat.defineType("fail");
  std::map<std::string, TH1 *> hmap = {{"pass", dP.get()}, {"fail", dF.get()}};
  RooDataHist data("data", "data", RooArgList(m), cat, hmap);
  RooDataHist rtP("rtP", "", RooArgList(m), tP.get()), rtF("rtF", "", RooArgList(m), tF.get());
  RooHistPdf  hP("hP", "", RooArgSet(m), rtP, 1), hF("hF", "", RooArgSet(m), rtF, 1);
  RooRealVar  mu("mu", "#Deltam", mu0, -5.0, 5.0), sg("sg", "#sigma", sg0, 0.05, 6.0);
  if (!floatConv) { mu.setConstant(true); sg.setConstant(true); }
  RooGaussian   gs("gs", "", m, mu, sg);
  RooFFTConvPdf sP("sP", "", m, hP, gs), sF("sF", "", m, hF, gs);
  RooRealVar     lP("lP", "", -0.02, -0.5, 0.2), lF("lF", "", -0.02, -0.5, 0.2);
  RooExponential eP("eP", "", m, lP), eF("eF", "", m, lF);
  RooRealVar     aP("aP", "", 0.0, -1.0, 1.0), aF("aF", "", 0.0, -1.0, 1.0);
  RooChebychev   cP("cP", "", m, RooArgList(aP)), cF("cF", "", m, RooArgList(aF));
  RooAbsPdf     *bP = (bkg == 0) ? (RooAbsPdf *)&eP : (RooAbsPdf *)&cP;
  RooAbsPdf     *bF = (bkg == 0) ? (RooAbsPdf *)&eF : (RooAbsPdf *)&cF;
  const double   nP = dP->Integral(), nF = dF->Integral();
  RooRealVar     nSig("nSig", "", 0.97 * (nP + nF), 0.0, 3.0 * (nP + nF) + 10.0);
  RooRealVar     eps("eps", "#epsilon", nP / std::max(1.0, nP + nF), 0.0, 1.0);
  RooFormulaVar  nSP("nSP", "@0*@1", RooArgList(nSig, eps)), nSF("nSF", "@0*(1-@1)", RooArgList(nSig, eps));
  RooRealVar     nBP("nBP", "", 0.02 * nP + 1.0, 0.0, nP + 10.0), nBF("nBF", "", 0.3 * nF + 1.0, 0.0, nF + 10.0);
  RooAddPdf      pP("pP", "", RooArgList(sP, *bP), RooArgList(nSP, nBP));
  RooAddPdf      pF("pF", "", RooArgList(sF, *bF), RooArgList(nSF, nBF));
  RooSimultaneous sim("sim", "", cat);
  sim.addPdf(pP, "pass");
  sim.addPdf(pF, "fail");
  std::unique_ptr<RooFitResult> r{sim.fitTo(data, RooFit::Extended(true), RooFit::Save(true),
                                            RooFit::Minos(RooArgSet(eps)), RooFit::Strategy(2),
                                            RooFit::PrintLevel(-1), RooFit::Verbose(false))};
  if (!r) return out;
  out.status  = r->status();
  out.covQual = r->covQual();
  const RooRealVar *e = (const RooRealVar *)r->floatParsFinal().find("eps");
  out.eps = e->getVal();
  out.lo  = std::fabs(e->getErrorLo());
  out.hi  = std::fabs(e->getErrorHi());
  if (out.lo == 0.0) out.lo = e->getError();
  if (out.hi == 0.0) out.hi = std::min(e->getError(), 1.0 - out.eps); // at the eps = 1 boundary
  out.nSig = nSig.getVal();
  out.nBP  = nBP.getVal();
  out.nBF  = nBF.getVal();
  out.mu   = mu.getVal();
  out.sg   = sg.getVal();
  out.ok   = (out.status == 0);
  if (lo <= 80.0 && hi >= 100.0)
  { // the fitted fail background inside 80-100 GeV, for the cut-and-count cross-check
    m.setRange("win8010", 80.0, 100.0);
    std::unique_ptr<RooAbsReal> fr(bF->createIntegral(m, RooFit::NormSet(m), RooFit::Range("win8010")));
    out.bFwin = out.nBF * fr->getVal();
  }

  if (plot.Length())
  { // the repo's standard postfit plot (plotting_helper.C SaveNicePlot1D_WithBkg:
    // data over the stacked fitted components, pull pad), one per category.
    // Component histograms = the fitted pdfs at the bin centres x bin width,
    // normalized to the fitted yields (the stack sums to the fit exactly).
    RooArgSet ns(m);
    const TString tag = TString(gSystem->BaseName(plot.Data()));
    auto comp = [&](RooAbsPdf &pdf, double yield, const TH1 *ref, const TString &name)
    {
      TH1D *hc = (TH1D *)ref->Clone(name);
      hc->SetDirectory(nullptr);
      hc->Reset();
      double sum = 0.0;
      for (int b = 1; b <= hc->GetNbinsX(); ++b)
      {
        m.setVal(hc->GetBinCenter(b));
        const double v = pdf.getVal(&ns) * hc->GetBinWidth(b);
        hc->SetBinContent(b, v);
        sum += v;
      }
      if (sum > 0.0) hc->Scale(yield / sum);
      for (int b = 1; b <= hc->GetNbinsX(); ++b) hc->SetBinError(b, 0.0);
      return hc;
    };
    const char *cn[2] = {"pass", "fail"};
    for (int i = 0; i < 2; ++i)
    {
      const bool   isP  = (i == 0);
      const double ySig = isP ? out.eps * out.nSig : (1.0 - out.eps) * out.nSig;
      const double yBkg = isP ? out.nBP : out.nBF;
      // kept alive: the helper's canvas references them (it is never deleted)
      TH1D *hd = (TH1D *)(isP ? dP.get() : dF.get())->Clone(Form("%s_%s_data", tag.Data(), cn[i]));
      hd->SetDirectory(nullptr);
      TH1D *hs = comp(isP ? (RooAbsPdf &)sP : (RooAbsPdf &)sF, ySig, hd, Form("%s_%s_sig", tag.Data(), cn[i]));
      TH1D *hb = comp(isP ? *bP : *bF, yBkg, hd, Form("%s_%s_bkg", tag.Data(), cn[i]));
      PlotStyle ps;
      ps.drawOpt       = "hist";
      ps.showStats     = false;
      ps.normBkgToData = false; // absolute fitted yields
      ps.pullPad       = true;
      ps.headerX       = 0.45;
      ps.boxTextSize   = 0.028;
      ps.boxX1 = 0.56; ps.boxX2 = 0.93; ps.boxY1 = 0.56; ps.boxY2 = 0.72;
      ps.legX1 = 0.56; ps.legY1 = 0.40; ps.legX2 = 0.93; ps.legY2 = 0.55;
      PlotTuner tuner = [](TCanvas *, TH1 *h)
      { // headroom from the tallest error bar: the fail spectra are 0-3 events
        // per bin, whose Poisson bars reach far above the bin contents
        if (!h) return;
        double top = 0.0;
        for (int b = 1; b <= h->GetNbinsX(); ++b) top = std::max(top, h->GetBinContent(b) + h->GetBinError(b));
        h->SetMinimum(0.0);
        h->SetMaximum(1.75 * top);
      };
      const std::vector<std::string> box = {
          Form("signal %.1f, background %.1f", ySig, yBkg),
          Form("#varepsilon_{data} = %.3f^{+%.3f}_{-%.3f}", out.eps, out.hi, out.lo),
          Form("fit status %d, covQual %d", out.status, out.covQual)};
      SaveNicePlot1D_WithBkg(hd, {hs, hb}, {"Z#rightarrowee (DY template)", "background"},
                             (plot + "_" + cn[i]).Data(), "m(tag, probe) [GeV]", "pairs / GeV", "",
                             "Z#rightarrowee tag-and-probe, electron reco",
                             Form("%s probes, %s", cn[i], title.Data()), box, ps, tuner);
    }
  }
  return out;
}
} // namespace pORecoSF

int RecoSfFit()
{
  using namespace pORecoSF;
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  TFile *fd = TFile::Open(Form("%s/skim_Data.root", kOutDir));
  TFile *fm = TFile::Open(Form("%s/skim_DY.root", kOutDir));
  if (!fd || fd->IsZombie() || !fm || fm->IsZombie())
  {
    std::cerr << "[FATAL] run the skim first: " << kOutDir << "/skim_{Data,DY}.root\n";
    return 2;
  }
  gSystem->mkdir(kPlotDir, kTRUE);
  std::ofstream csv(Form("%s/reco_sf_ele.csv", kPlotDir));
  csv << "variant,scheme,bin,label,n_pass_data,n_fail_data,eps_data,eps_data_errlo,eps_data_errhi,eps_data_syst,"
         "eps_mc,eps_mc_err,sf,sf_errlo,sf_errhi,sf_syst,egm_sf,nbkg_pass,nbkg_fail,fit_status,cov_qual,mu,sigma\n";

  struct Row { double sf, lo, hi, syst, egm, cc; bool ok; };
  std::vector<Row> rows[kNVar];
  // the machine-readable result (2026-09-30): one TH1D per (variant, scheme,
  // bin) with LABELED bins -> rootfile/reco_sf_ele/reco_sf_ele_result.root.
  // under kEleRecoSource = kRecoTnP skim/electron_sf.h applies recosf_noid_incl0
  // (nominal probe, inclusive) as THE electron Reco SF (2026-09-30 only; the pp
  // EGM tables since 2026-10-01); the CSV above stays the human-readable twin.
  const char *kResLab[] = {"sf", "sf_stat_lo", "sf_stat_hi", "sf_syst", "eps_data", "eps_data_lo", "eps_data_hi",
                           "eps_data_syst", "eps_mc", "eps_mc_err", "egm_sf", "sf_cc", "n_pass_data", "n_fail_data",
                           "fit_status", "cov_qual"};
  const int kNResLab = sizeof(kResLab) / sizeof(kResLab[0]);
  std::vector<TH1D *> results;

  for (int v = 0; v < kNVar; ++v)
  {
    double muIncl = 0.0, sgIncl = 0.5;
    std::cout << "\n==================== " << kVarLabel[v] << " ====================\n";
    for (int s = 0; s < kNScheme; ++s)
      for (int k = 0; k < kSchemes[s].n; ++k)
      {
        TH1 *dP = (TH1 *)fd->Get(RecoHName(v, s, k, true)), *dF = (TH1 *)fd->Get(RecoHName(v, s, k, false));
        TH1 *tP = (TH1 *)fm->Get(RecoHName(v, s, k, true)), *tF = (TH1 *)fm->Get(RecoHName(v, s, k, false));
        if (!dP || !dF || !tP || !tF) { std::cerr << "[ERROR] missing histograms for " << RecoHName(v, s, k, true) << "\n"; continue; }
        const TString lab = (kSchemes[s].n == 1) ? TString("inclusive") : TString(kSchemes[s].label[k]);

        // eps_MC: counted from the templates over the nominal range
        double eP = 0.0, eF = 0.0;
        const double mP = tP->IntegralAndError(1, tP->GetNbinsX(), eP), mF = tF->IntegralAndError(1, tF->GetNbinsX(), eF);
        const double epsMC = mP / (mP + mF);
        const double epsMCe = std::sqrt(std::pow((1 - epsMC) * eP, 2) + std::pow(epsMC * eF, 2)) / (mP + mF);

        // nominal fit: the Gaussian floats in the inclusive fit, is fixed to it per bin
        const bool isIncl = (s == 0);
        const TString plot = Form("%s/fit_%s_%s%d", kPlotDir, kVarName[v], kSchemes[s].name, k);
        const TString title = (v == 0) ? lab : TString(Form("%s, loose-ID probe", lab.Data()));
        const FitOut nom = RecoFitEps(dP, dF, tP, tF, kMLo, kMHi, isIncl, isIncl ? 0.0 : muIncl, isIncl ? 0.5 : sgIncl,
                                      0, plot, title);
        if (isIncl) { muIncl = nom.mu; sgIncl = nom.sg; }

        // systematics on the nominal variant
        double syst = 0.0;
        if (v == kNoId)
        {
          const FitOut sB = RecoFitEps(dP, dF, tP, tF, kMLo, kMHi, isIncl, isIncl ? 0.0 : muIncl, isIncl ? 0.5 : sgIncl, 1, "", "");
          const FitOut sR = RecoFitEps(dP, dF, tP, tF, kSysLo, kSysHi, isIncl, isIncl ? 0.0 : muIncl, isIncl ? 0.5 : sgIncl, 0, "", "");
          // the Gaussian: floated per bin (bins) / fixed at no smearing (inclusive: mu 0, sigma 0.05)
          const FitOut sG = isIncl ? RecoFitEps(dP, dF, tP, tF, kMLo, kMHi, false, 0.0, 0.05, 0, "", "")
                                   : RecoFitEps(dP, dF, tP, tF, kMLo, kMHi, true, muIncl, sgIncl, 0, "", "");
          // a variation counts only if its fit converged: a failed fit (e.g. eps
          // driven to the eps = 1 boundary when a narrower range lets the
          // background swallow every fail) is not a systematic shift
          const FitOut *vars[3] = {&sB, &sR, &sG};
          const char   *vnm[3]  = {"bkg Chebychev", "range 70-110", isIncl ? "Gaussian unsmeared" : "Gaussian floated"};
          std::string ignored;
          for (int iv = 0; iv < 3; ++iv)
          {
            if (vars[iv]->status == 0) syst = std::max(syst, std::fabs(vars[iv]->eps - nom.eps));
            else ignored += std::string(ignored.empty() ? "" : ", ") + vnm[iv];
          }
          std::cout << Form("  [SYST] %-24s %s %+.4f (st %d) | %s %+.4f (st %d) | %s %+.4f (st %d) -> %.4f%s\n",
                            lab.Data(), vnm[0], sB.eps - nom.eps, sB.status, vnm[1], sR.eps - nom.eps, sR.status, vnm[2],
                            sG.eps - nom.eps, sG.status, syst,
                            ignored.empty() ? "" : Form("   (not converged, ignored: %s)", ignored.c_str()));
        }
        if (nom.status != 0)
          std::cout << Form("  [WARN] %s %s: the NOMINAL fit did not converge (status %d) -- do not quote this bin\n",
                            kVarName[v], lab.Data(), nom.status);

        // cut-and-count cross-check, 80-100 GeV: raw pass, raw fail minus the fitted
        // fail background inside the window; MC counted in the same window
        double sfCC = 0.0;
        {
          const int w1 = dP->FindBin(80.5), w2 = dP->FindBin(99.5);
          const double pw = dP->Integral(w1, w2), fw = dF->Integral(w1, w2);
          const double mpw = tP->Integral(w1, w2), mfw = tF->Integral(w1, w2);
          const double fsig = std::max(0.0, fw - nom.bFwin);
          const double eCC = pw / (pw + fsig), eMCw = mpw / (mpw + mfw);
          sfCC = eCC / eMCw;
          std::cout << Form("  [CC]   %-7s %-22s 80-100: pass %.0f, fail %.0f - fitted bkg %.1f = %.1f -> eps %.4f | MC %.4f | SF %.4f\n",
                            kVarName[v], lab.Data(), pw, fw, nom.bFwin, fsig, eCC, eMCw, eCC / eMCw);
        }

        // EGM prediction: <SF_EGM> over the passing DY-MC probes of this bin
        TH1 *hs = (TH1 *)fm->Get(Form("h_egmsum_%s_%s", kVarName[v], kSchemes[s].name));
        TH1 *hw = (TH1 *)fm->Get(Form("h_egmw_%s_%s", kVarName[v], kSchemes[s].name));
        const double egmSF = (hs && hw && hw->GetBinContent(k + 1) > 0) ? hs->GetBinContent(k + 1) / hw->GetBinContent(k + 1) : 0.0;

        const double sf = nom.eps / epsMC;
        const double sfLo = sf * std::sqrt(std::pow(nom.lo / nom.eps, 2) + std::pow(epsMCe / epsMC, 2));
        const double sfHi = sf * std::sqrt(std::pow(nom.hi / nom.eps, 2) + std::pow(epsMCe / epsMC, 2));
        const double sfSyst = syst / epsMC;
        rows[v].push_back({sf, sfLo, sfHi, sfSyst, egmSF, sfCC, nom.ok});
        std::cout << Form("[RESULT] %-7s %-22s data pass %4.0f fail %3.0f | eps_data %.4f +%.4f -%.4f (syst %.4f) | eps_MC %.4f +- %.4f"
                          " | SF %.4f +%.4f -%.4f (syst %.4f) | EGM %.4f | bkg pass %.1f fail %.1f | status %d cov %d | dm %.2f sigma %.2f\n",
                          kVarName[v], lab.Data(), dP->Integral(), dF->Integral(), nom.eps, nom.hi, nom.lo, syst, epsMC,
                          epsMCe, sf, sfHi, sfLo, sfSyst, egmSF, nom.nBP, nom.nBF, nom.status, nom.covQual, nom.mu, nom.sg);
        csv << Form("%s,%s,%d,%s,%.0f,%.0f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.2f,%.2f,%d,%d,%.3f,%.3f\n",
                    kVarName[v], kSchemes[s].name, k, lab.Data(), dP->Integral(), dF->Integral(), nom.eps, nom.lo, nom.hi,
                    syst, epsMC, epsMCe, sf, sfLo, sfHi, sfSyst, egmSF, nom.nBP, nom.nBF, nom.status, nom.covQual, nom.mu,
                    nom.sg);
        {
          TH1D *hr = new TH1D(Form("recosf_%s_%s%d", kVarName[v], kSchemes[s].name, k),
                              Form("electron Reco SF, pO Z#rightarrowee tag-and-probe: %s, %s", kVarLabel[v], lab.Data()),
                              kNResLab, 0, kNResLab);
          hr->SetDirectory(nullptr);
          const double val[] = {sf, sfLo, sfHi, sfSyst, nom.eps, nom.lo, nom.hi, syst, epsMC, epsMCe, egmSF, sfCC,
                                dP->Integral(), dF->Integral(), (double)nom.status, (double)nom.covQual};
          for (int i = 0; i < kNResLab; ++i)
          {
            hr->GetXaxis()->SetBinLabel(i + 1, kResLab[i]);
            hr->SetBinContent(i + 1, val[i]);
          }
          results.push_back(hr);
        }
      }
  }
  // MC truth: the efficiency of probes that sit on a prompt Z electron
  {
    TH1 *gp = (TH1 *)fm->Get(RecoHName(kNoId, 0, 0, true, "_gm")), *gf = (TH1 *)fm->Get(RecoHName(kNoId, 0, 0, false, "_gm"));
    if (gp && gf)
      std::cout << Form("\n[MC-TRUTH] probes on a prompt Z electron: eps_MC = %.4f (vs %.4f for all probes, 60-120 GeV)\n",
                        gp->Integral() / (gp->Integral() + gf->Integral()),
                        ((TH1 *)fm->Get(RecoHName(kNoId, 0, 0, true)))->Integral() /
                            (((TH1 *)fm->Get(RecoHName(kNoId, 0, 0, true)))->Integral() + ((TH1 *)fm->Get(RecoHName(kNoId, 0, 0, false)))->Integral()));
  }
  csv.close();
  {
    const TString resName = Form("%s/reco_sf_ele_result.root", kOutDir);
    TDirectory::TContext ctx;
    TFile *fr = TFile::Open(resName, "RECREATE");
    if (!fr || fr->IsZombie()) { std::cerr << "[FATAL] cannot write " << resName << "\n"; return 2; }
    for (TH1D *hr : results) hr->Write();
    fr->Close();
    delete fr;
    std::cout << "[OUT] " << resName << " (" << results.size()
              << " labeled result histograms; skim/electron_sf.h applies recosf_noid_incl0)\n";
  }

  // ---- summary: the repo's standard overlay (plotting_helper.C
  //      SaveNiceGraph_Overlay), SF vs |eta_SC| and vs E_T^SC, one series each:
  //      the fit (stat bars + syst boxes; bins whose nominal fit did not
  //      converge are left out and named in the box), the 80-100 GeV
  //      cut-and-count, and the pp EGM Reco SF for the same probes. The
  //      inclusive numbers go to the info box. The fit points are drawn with
  //      their ASYMMETRIC MINOS intervals (OverlaySeries gTotAsym/gStatAsym,
  //      2026-09-30 -- until then the helper took symmetric errors only and the
  //      plot showed their mean, while the box printed the asymmetric
  //      inclusive interval); the syst box sits on top of each side.
  //      rows[v] follows the fit loop: [0] inclusive, [1,2] |eta| bins, [3,4] E_T bins.
  struct AxisDef { const char *name, *xTitle; int first; double x[2], ex[2]; const char *binNote[2]; };
  const AxisDef axes[2] = {
      {"eta", "|#eta_{SC}|", 1, {0.7221, 1.983}, {0.7221, 0.417}, {"barrel", "endcap"}},
      {"et", "E_{T}^{SC} (GeV)", 3, {32.5, 72.5}, {12.5, 27.5}, {"E_{T} 20-45 GeV", "E_{T} > 45 GeV"}}};
  const Row &inc = rows[kNoId][0];
  for (const AxisDef &ax : axes)
  {
    TGraphErrors *gFitT = new TGraphErrors(), *gFitS = new TGraphErrors();
    TGraphAsymmErrors *gFitTA = new TGraphAsymmErrors(), *gFitSA = new TGraphAsymmErrors();
    TGraphErrors *gCC = new TGraphErrors(2), *gEgm = new TGraphErrors(2);
    gFitT->SetName(Form("g_recosf_fit_%s", ax.name));
    gFitS->SetName(Form("g_recosf_fitstat_%s", ax.name));
    gFitTA->SetName(Form("g_recosf_fitasym_%s", ax.name));
    gFitSA->SetName(Form("g_recosf_fitstatasym_%s", ax.name));
    gCC->SetName(Form("g_recosf_cc_%s", ax.name));
    gEgm->SetName(Form("g_recosf_egm_%s", ax.name));
    std::string skipped;
    for (int k = 0; k < 2; ++k)
    {
      const Row   &r    = rows[kNoId][ax.first + k];
      const double stat = 0.5 * (r.lo + r.hi);
      if (r.ok)
      { // symmetric twins = the frame (points, x errors); the asymmetric ones are drawn
        const int ip = gFitT->GetN();
        gFitT->SetPoint(ip, ax.x[k], r.sf);
        gFitT->SetPointError(ip, ax.ex[k], std::hypot(stat, r.syst));
        gFitS->SetPoint(ip, ax.x[k], r.sf);
        gFitS->SetPointError(ip, ax.ex[k], stat);
        gFitTA->SetPoint(ip, ax.x[k], r.sf);
        gFitTA->SetPointError(ip, ax.ex[k], ax.ex[k], std::hypot(r.lo, r.syst), std::hypot(r.hi, r.syst));
        gFitSA->SetPoint(ip, ax.x[k], r.sf);
        gFitSA->SetPointError(ip, ax.ex[k], ax.ex[k], r.lo, r.hi);
      }
      else
        skipped += std::string(skipped.empty() ? "" : ", ") + ax.binNote[k];
      gCC->SetPoint(k, ax.x[k], r.cc);
      gCC->SetPointError(k, ax.ex[k], 0.0);
      gEgm->SetPoint(k, ax.x[k], r.egm);
      gEgm->SetPointError(k, ax.ex[k], 0.0);
    }
    std::vector<OverlaySeries> ser(3);
    ser[0].gTot = gFitT; ser[0].gStat = gFitS; ser[0].label = "pO data / MC, fit";
    ser[0].gTotAsym = gFitTA; ser[0].gStatAsym = gFitSA; // drawn: the MINOS intervals
    ser[0].color = kBlue + 1;  ser[0].marker = 20; ser[0].markerSize = 1.3; ser[0].xShift = -0.35;
    ser[1].gTot = gCC; ser[1].label = "pO data / MC, cut-and-count 80-100 GeV";
    ser[1].color = kGray + 2;  ser[1].marker = 24; ser[1].markerSize = 1.3; ser[1].xShift = 0.0;
    ser[2].gTot = gEgm; ser[2].label = "pp EGM Reco SF, same probes";
    ser[2].color = kGreen + 2; ser[2].marker = 33; ser[2].markerSize = 2.0; ser[2].xShift = 0.35;
    std::vector<std::string> box = {
        Form("inclusive: %.3f^{+%.3f}_{-%.3f} (stat) #pm %.3f (syst)", inc.sf, inc.hi, inc.lo, inc.syst),
        Form("cut-and-count %.3f, pp EGM %.3f", inc.cc, inc.egm)};
    if (!skipped.empty()) box.push_back(Form("fit not converged: %s", skipped.c_str()));
    // header + box right of the CMS label, box under the header, legend lower
    // left; the y range puts SF = 1 at ~45% of the frame so nothing meets a point
    PlotStyle ps;
    ps.headerX = 0.50;
    ps.headerY = 0.90;
    ps.subSize = 0.034;
    ps.boxTextSize = 0.028;
    ps.boxX1 = 0.40; ps.boxX2 = 0.93; ps.boxY1 = 0.56; ps.boxY2 = 0.72;
    ps.legX1 = 0.18; ps.legY1 = 0.15; ps.legX2 = 0.70; ps.legY2 = 0.33;
    GraphTuner tuner = [](TCanvas *c, TGraphErrors *g)
    {
      g->SetMinimum(0.75);
      g->SetMaximum(1.30);
      c->cd();
      TLine *one = new TLine(g->GetXaxis()->GetXmin(), 1.0, g->GetXaxis()->GetXmax(), 1.0);
      one->SetLineStyle(2);
      one->SetLineColor(kGray + 2);
      one->Draw();
    };
    SaveNiceGraph_Overlay(ser, Form("%s/reco_sf_vs_%s", kPlotDir, ax.name), ax.xTitle,
                          "Reco SF (#varepsilon_{data} / #varepsilon_{MC})", "",
                          "Z#rightarrowee tag-and-probe", "electron reco efficiency", box, ps, tuner);
  }
  std::cout << "[OUT] " << kPlotDir << "/{fit_*_{pass,fail},reco_sf_vs_{eta,et}}.{png,pdf}, reco_sf_ele.csv\n";
  return 0;
}

int reco_sf_ele(const char *mode = "all", const char *sample = "", Long64_t nmax = -1)
{
  const std::string md = mode, sm = sample;
  if (md == "skim")
  {
    if (!sm.empty()) return RecoSfSkim(sm, nmax);
    for (const char *s : {"Data", "DY"})
      if (const int rc = RecoSfSkim(s, nmax)) return rc;
    return 0;
  }
  if (md == "fit") return RecoSfFit();
  if (md == "all")
  {
    for (const char *s : {"Data", "DY"})
      if (const int rc = RecoSfSkim(s, nmax)) return rc;
    return RecoSfFit();
  }
  std::cerr << "usage: reco_sf_ele.C+(\"skim\"[, \"Data\"|\"DY\"] | \"fit\" | \"all\")\n";
  return 1;
}
