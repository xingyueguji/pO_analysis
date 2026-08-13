// skim/skim.C
//
// Unified entry point for the four pO skim channels.
// Dispatch through ChannelType:
//   kZmm -> Z -> mu+mu
//   kZee -> Z -> e+e
//   kWmu -> W -> mu+nu
//   kWel -> W -> e+nu
//
// Each per-channel function reproduces, byte-for-byte in selection logic
// and histogram output, the legacy macros in skim/legacy/:
//   skim_Wmu : legacy/DrawWToMuNu_PFMet.C
//   skim_Wel : legacy/DrawWToElecNu_PFMet.C
//   skim_Zmm : legacy/DrawDimuonPeak.C
//   skim_Zee : legacy/DrawDielectronPeak.C
//
// Driver invocation (from run_all.sh):
//   root -l -q -b 'skim.C+(<channel-enum>, "<root-file>", <sample-enum>)'
//
// Output filenames and histogram names are unchanged.

#include "skim_common.h"

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TVector2.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TBranch.h"
#include "TPaveText.h"
#include "TInterpreter.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace pOSkim;

// ============================================================
// Per-channel selection helpers (kept here, not in skim_common.h,
// because they are too entangled with each channel's branch set
// to factor cleanly into a shared header).
// ============================================================

// -------- W -> mu nu : leading-muon definition & DY veto -------------

namespace {

// Step 1 helper: exists a PF muon with pt > 25
static bool ExistsPFMuonPt25_W(double ptCut,
                               int nMu,
                               std::vector<float> *muPt,
                               bool has_muIsPF, std::vector<int> *muIsPF)
{
  if (!muPt) return false;
  for (int i = 0; i < nMu; ++i)
  {
    const double pt = muPt->at(i);
    if (pt <= ptCut) continue;
    if (has_muIsPF && muIsPF && muIsPF->at(i) == 0) continue;
    return true;
  }
  return false;
}

// Step 5 helper: exists Tight ID muon
static bool ExistsTightMuon_W(int nMu,
                              bool has_muIDTight, std::vector<int> *muIDTight)
{
  if (!has_muIDTight || !muIDTight) return false;
  for (int i = 0; i < nMu; ++i)
    if (muIDTight->at(i) != 0) return true;
  return false;
}

static int FindLeadingMuon_TightPF(int nMu,
                                   std::vector<float> *muPt,
                                   std::vector<float> *muEta,
                                   std::vector<float> *muPhi,
                                   bool has_muIDTight, std::vector<int> *muIDTight,
                                   bool has_muIsPF,    std::vector<int> *muIsPF)
{
  if (!muPt || !muEta || !muPhi) return -1;

  int    best   = -1;
  double bestPt = -1;

  for (int i = 0; i < nMu; ++i)
  {
    const double pt = muPt->at(i);
    if (has_muIDTight && muIDTight && muIDTight->at(i) == 0) continue;
    if (has_muIsPF    && muIsPF    && muIsPF->at(i)    == 0) continue;
    if (pt > bestPt) { bestPt = pt; best = i; }
  }
  return best;
}

// Step 4 (muon W): OS dimuon DY veto. Veto if exists OS pair where both
// muons are Tight+PF, pT > 15, relIso < isoMax, AND mll in (80, 110).
static bool PassDYVeto_Wmu(double dyMuPtMin, double isoMax,
                           double dyMassMin, double dyMassMax,
                           int nMu,
                           std::vector<float> *muPt,
                           std::vector<float> *muEta,
                           std::vector<float> *muPhi,
                           std::vector<int>   *muCharge,
                           bool has_muIDTight, std::vector<int> *muIDTight,
                           bool has_muIsPF,    std::vector<int> *muIsPF,
                           std::vector<float> *muPFChIso,
                           std::vector<float> *muPFNeuIso,
                           std::vector<float> *muPFPhoIso,
                           std::vector<float> *muPFPUIso)
{
  if (!muPt || !muEta || !muPhi || !muCharge) return true;
  if (!has_muIDTight || !muIDTight)            return true;

  std::vector<int> cand;
  cand.reserve(nMu);

  for (int i = 0; i < nMu; ++i)
  {
    if (muPt->at(i) <= dyMuPtMin) continue;
    if (muIDTight->at(i) == 0)     continue;
    if (has_muIsPF && muIsPF && muIsPF->at(i) == 0) continue;
    if (RelIsoPF(i, muPt, muPFChIso, muPFNeuIso, muPFPhoIso, muPFPUIso) >= isoMax) continue;
    cand.push_back(i);
  }

  for (size_t a = 0; a < cand.size(); ++a)
    for (size_t b = a + 1; b < cand.size(); ++b)
    {
      const int i1 = cand[a], i2 = cand[b];
      if (muCharge->at(i1) * muCharge->at(i2) >= 0) continue;
      TLorentzVector p1, p2;
      p1.SetPtEtaPhiM(muPt->at(i1), muEta->at(i1), muPhi->at(i1), MU_MASS);
      p2.SetPtEtaPhiM(muPt->at(i2), muEta->at(i2), muPhi->at(i2), MU_MASS);
      const double mll = (p1 + p2).M();
      if (mll > dyMassMin && mll < dyMassMax) return false;
    }
  return true;
}

// -------- W -> e nu --------------------------------------------------

// Step 1: exists electron with pt > 25 (no PF flag for electrons)
static bool ExistsPFElectronPt25_W(double ptCut,
                                   int nEle,
                                   std::vector<float> *elePt)
{
  if (!elePt) return false;
  for (int i = 0; i < nEle; ++i)
    if (elePt->at(i) > ptCut) return true;
  return false;
}

// Step 5: exists Tight ID electron. NOTE: legacy macro had no nullptr
// guard on eleIDTight — preserved exactly.
static bool ExistsTightElectron_W(int nEle, std::vector<int> *eleIDTight)
{
  for (int i = 0; i < nEle; ++i)
    if (eleIDTight->at(i) != 0) return true;
  return false;
}

static int FindLeadingElectron_TightPF(int nEle,
                                       std::vector<float> *elePt,
                                       std::vector<float> *eleEta,
                                       std::vector<float> *elePhi,
                                       std::vector<int>   *eleIDTight)
{
  if (!elePt || !eleEta || !elePhi) return -1;

  int    best   = -1;
  double bestPt = -1;

  for (int i = 0; i < nEle; ++i)
  {
    if (eleIDTight && eleIDTight->at(i) == 0) continue;
    const double pt = elePt->at(i);
    if (pt > bestPt) { bestPt = pt; best = i; }
  }
  return best;
}

// Step 4 (electron W): OS dielectron DY veto. Veto if exists OS pair where
// both electrons are ID'd (eleMVAIdWP95) with continuous PF relIso < isoMax,
// pT > dyElePtMin, AND mll in (dyMassMin, dyMassMax) — same surface as the
// muon veto. (Iso gate switched from the integer eleMVAIsoWP95 WP to the
// continuous RelIsoPF on 2026-07-02: the isolation study characterized the
// continuous cut, and the MVA-Iso WP was shown NOT to track it.)
static bool PassDYVeto_Wel(double dyElePtMin, double isoMax,
                           double dyMassMin, double dyMassMax,
                           int nEle,
                           std::vector<float> *elePt,
                           std::vector<float> *eleEta,
                           std::vector<float> *elePhi,
                           std::vector<int>   *eleCharge,
                           std::vector<int>   *eleIDTight,
                           std::vector<float> *elePFChIso,
                           std::vector<float> *elePFNeuIso,
                           std::vector<float> *elePFPhoIso,
                           std::vector<float> *elePFPUIso)
{
  if (!elePt || !eleEta || !elePhi || !eleCharge) return true;

  std::vector<int> cand;
  cand.reserve(nEle);

  for (int i = 0; i < nEle; ++i)
  {
    if (elePt->at(i) <= dyElePtMin) continue;
    if (eleIDTight->at(i) == 0)      continue;
    if (RelIsoPF(i, elePt, elePFChIso, elePFNeuIso, elePFPhoIso, elePFPUIso) >= isoMax) continue;
    cand.push_back(i);
  }

  for (size_t a = 0; a < cand.size(); ++a)
    for (size_t b = a + 1; b < cand.size(); ++b)
    {
      const int i1 = cand[a], i2 = cand[b];
      if (eleCharge->at(i1) * eleCharge->at(i2) >= 0) continue;
      TLorentzVector p1, p2;
      p1.SetPtEtaPhiM(elePt->at(i1), eleEta->at(i1), elePhi->at(i1), ELE_MASS);
      p2.SetPtEtaPhiM(elePt->at(i2), eleEta->at(i2), elePhi->at(i2), ELE_MASS);
      const double mll = (p1 + p2).M();
      if (mll > dyMassMin && mll < dyMassMax) return false;
    }
  return true;
}

} // anonymous namespace


// ============================================================
// W -> mu+nu
// Legacy: skim/legacy/DrawWToMuNu_PFMet.C
// ============================================================

int skim_Wmu(const char *fname, SampleType sample)
{
  gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");

  const bool isMC = IsMC(sample);

  gStyle->SetOptStat(0);

  // -------- Config (kept identical to legacy WConfig) --------
  const double muPt25      = 25.0;
  const double dyMuPtMin   = 15.0;
  const double isoMax      = 0.15;
  // m_T cut of the lepton-pT-DISCRIMINANT selection (pT > 25 && m_T > 40).
  // Applies ONLY to the `_mt40` histogram set below -- the nominal MET/m_T
  // histograms and the cutflow are untouched (no MET/m_T cut in the W skim).
  // NAMING CONTRACT: the literal "40" in the persisted `_mt40` histogram names
  // asserts this value. If you change it, RENAME those histograms too (and the
  // readers: correction/qcd_abcd.C runPtCharge tag, plotting/mtandmet.C
  // varStem[kVarMt40], correction/dataMC_kinematics.C's _mt40 channels) --
  // otherwise the stored names silently claim a cut that was not applied.
  const double mtCutForPtDisc     = 40.0;
  const double dyMassMin   = 80.0;
  const double dyMassMax   = 110.0;
  const double muEtaMax    = 2.4;
  // 0.1 -> 0.4 (2026-08-12): the saved trigger objects are L1-granularity
  // candidates whose phi is the muon's azimuth AT THE MUON STATION, while
  // {mu,ele}Phi is the vertex direction. The two differ by the solenoid
  // bending, Dphi = c + q*k/pT with k = 2.55 rad*GeV and a charge-blind
  // c = +0.006 rad (half of the 2pi/576 muGMT phi cell). At pT 25-30 in the
  // barrel that is ~0.10 rad, i.e. ON the old window -- and because c adds to
  // the mu+ bending and cancels part of the mu- one, the old cut kept only
  // 48% of barrel mu+ vs 91% of mu- there (a ~10% CHARGE-DEPENDENT loss riding
  // straight on the charge-asymmetry observable). 0.4 is comfortably beyond
  // the largest |Dphi| (~0.10) and restores 100%/100% for both charges.
  const double trigMatchDR = 0.4;
  const bool   applyVz     = true;
  const double vzMax       = 15.0;

  std::string outPrefix = "WToMuNu_pO_PFMet";

  // Sample -> filename / outPrefix override
  if (isMC)
  {
    const auto info = ResolveMCSample(sample, "mu");
    if (!info.fname.empty()) fname = Form("%s", info.fname.c_str());
    outPrefix += info.outSuffix;
  }

  std::cout << "[INPUT] " << fname << std::endl;
  TFile *f = TFile::Open(fname);
  if (!f || f->IsZombie()) { std::cerr << "ERROR opening file\n"; return 1; }

  TTree *tMu     = (TTree *)f->Get("ggHiNtuplizer/EventTree");
  TTree *tHi     = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  TTree *tPF     = (TTree *)f->Get("particleFlowAnalyser/pftree");
  TTree *tHLT    = (TTree *)f->Get("hltanalysis/HltTree");
  TTree *tHLTobj = (TTree *)f->Get("hltobject/HLT_OxyL1SingleMuOpen_v");
  TTree *tEvent  = (TTree *)f->Get("skimanalysis/HltTree");
  if (!tEvent)
  {
    std::cerr << "[FATAL] skim_Wmu: missing TTree skimanalysis/HltTree in " << fname << "\n";
    f->Close();
    return 2;
  }
  TTree *tGen    = nullptr;
  if (isMC) tGen = (TTree *)f->Get("HiGenParticleAna/hi");
  if (isMC && !tGen)
  {
    std::cerr << "[FATAL] skim_Wmu: isMC=true but no gen tree found (HiGenParticleAna/hi) in " << fname << "\n";
    f->Close();
    return 2;
  }

  if (!tMu || !tHi || !tPF || !tHLT || !tHLTobj)
  {
    std::cerr << "ERROR: missing one of trees: ggHiNtuplizer/EventTree, hiEvtAnalyzer/HiTree, particleFlowAnalyser/pftree\n";
    return 1;
  }

  // -------- Muon branches --------
  Int_t nMu = 0;
  std::vector<float> *muPt = nullptr, *muEta = nullptr, *muPhi = nullptr;
  std::vector<int>   *muCharge = nullptr;
  std::vector<int>   *muIDTight = nullptr, *muIsPF = nullptr;
  std::vector<float> *muPFChIso = nullptr, *muPFNeuIso = nullptr, *muPFPhoIso = nullptr;
  std::vector<float> *muPFPUIso = nullptr;

  tMu->SetBranchStatus("*", 0);
  tMu->SetBranchStatus("nMu", 1);
  tMu->SetBranchStatus("muPt", 1);
  tMu->SetBranchStatus("muEta", 1);
  tMu->SetBranchStatus("muPhi", 1);
  tMu->SetBranchStatus("muCharge", 1);
  if (HasBranch(tMu, "muIDTight"))   tMu->SetBranchStatus("muIDTight", 1);
  if (HasBranch(tMu, "muIsPF"))      tMu->SetBranchStatus("muIsPF", 1);
  if (HasBranch(tMu, "muPFChIso"))   tMu->SetBranchStatus("muPFChIso", 1);
  if (HasBranch(tMu, "muPFNeuIso"))  tMu->SetBranchStatus("muPFNeuIso", 1);
  if (HasBranch(tMu, "muPFPhoIso"))  tMu->SetBranchStatus("muPFPhoIso", 1);
  if (HasBranch(tMu, "muPFPUIso"))   tMu->SetBranchStatus("muPFPUIso", 1);

  tMu->SetBranchAddress("nMu",      &nMu);
  tMu->SetBranchAddress("muPt",     &muPt);
  tMu->SetBranchAddress("muEta",    &muEta);
  tMu->SetBranchAddress("muPhi",    &muPhi);
  tMu->SetBranchAddress("muCharge", &muCharge);

  const bool has_muIDTight  = HasBranch(tMu, "muIDTight");
  const bool has_muIsPF     = HasBranch(tMu, "muIsPF");
  if (has_muIDTight)  tMu->SetBranchAddress("muIDTight", &muIDTight);
  if (has_muIsPF)     tMu->SetBranchAddress("muIsPF",    &muIsPF);

  const bool has_muPFChIso  = HasBranch(tMu, "muPFChIso");
  const bool has_muPFNeuIso = HasBranch(tMu, "muPFNeuIso");
  const bool has_muPFPhoIso = HasBranch(tMu, "muPFPhoIso");
  const bool has_muPFPUIso  = HasBranch(tMu, "muPFPUIso");
  if (has_muPFChIso)  tMu->SetBranchAddress("muPFChIso",  &muPFChIso);
  if (has_muPFNeuIso) tMu->SetBranchAddress("muPFNeuIso", &muPFNeuIso);
  if (has_muPFPhoIso) tMu->SetBranchAddress("muPFPhoIso", &muPFPhoIso);
  if (has_muPFPUIso)  tMu->SetBranchAddress("muPFPUIso",  &muPFPUIso);
  if (!has_muPFPUIso)
    std::cout << "[WARN] skim_Wmu: no muPFPUIso branch; relIso falls back to uncorrected (no Delta-beta).\n";

  // -------- HLT bit (step 3) --------
  const std::string hltName = FindBranchContaining(tHLT, "HLT_OxyL1SingleMuOpen_v1");
  Int_t HLT_OxyL1SingleMuOpen_v1 = 0;
  bool  has_hlt = false;
  tHLT->SetBranchStatus("*", 0);
  if (!hltName.empty())
  {
    has_hlt = true;
    tHLT->SetBranchStatus(hltName.c_str(), 1);
    tHLT->SetBranchAddress(hltName.c_str(), &HLT_OxyL1SingleMuOpen_v1);
    std::cout << "[INFO] Using trigger bit branch: " << hltName << "\n";
  }
  else
  {
    std::cout << "[WARN] Could not find any branch containing 'HLT_OxyL1SingleMuOpen_v1' in hltanalysis/HltTree.\n"
              << "       Step 3 (trigger fired) will be treated as PASS.\n";
  }

  // -------- HLT objects (step 8) --------
  std::vector<double> *trgObjPt = nullptr, *trgObjEta = nullptr;
  std::vector<double> *trgObjPhi = nullptr, *trgObjId = nullptr;

  const bool has_trgObjPt  = HasBranch(tHLTobj, "pt");
  const bool has_trgObjEta = HasBranch(tHLTobj, "eta");
  const bool has_trgObjPhi = HasBranch(tHLTobj, "phi");
  const bool has_trgObjId  = HasBranch(tHLTobj, "TriggerObjID");

  tHLTobj->SetBranchStatus("*", false);
  if (has_trgObjPt && has_trgObjEta && has_trgObjPhi && has_trgObjId)
  {
    tHLTobj->SetBranchStatus("pt", 1);
    tHLTobj->SetBranchStatus("eta", 1);
    tHLTobj->SetBranchStatus("phi", 1);
    tHLTobj->SetBranchStatus("TriggerObjID", 1);

    tHLTobj->SetBranchAddress("pt",           &trgObjPt);
    tHLTobj->SetBranchAddress("eta",          &trgObjEta);
    tHLTobj->SetBranchAddress("phi",          &trgObjPhi);
    tHLTobj->SetBranchAddress("TriggerObjID", &trgObjId);
  }
  else
  {
    std::cout << "[WARN] Could not find HLT related Object." << std::endl;
  }

  // -------- HiTree + event filters (step 2) --------
  Float_t vz    = 999.f;
  Int_t   hiBin = -999;
  Int_t   pprimaryVertexFilter        = 1;
  Int_t   pclusterCompatibilityFilter = 1;

  const bool has_pprimaryVertexFilter        = HasBranch(tEvent, "pprimaryVertexFilter");
  const bool has_pclusterCompatibilityFilter = HasBranch(tEvent, "pclusterCompatibilityFilter");

  tEvent->SetBranchStatus("*", 0);
  if (has_pprimaryVertexFilter && has_pclusterCompatibilityFilter)
  {
    tEvent->SetBranchStatus("pprimaryVertexFilter", 1);
    tEvent->SetBranchStatus("pclusterCompatibilityFilter", 1);
    tEvent->SetBranchAddress("pprimaryVertexFilter",        &pprimaryVertexFilter);
    tEvent->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);
  }

  tHi->SetBranchStatus("*", 0);
  if (HasBranch(tHi, "vz"))    { tHi->SetBranchStatus("vz", 1);    tHi->SetBranchAddress("vz",    &vz); }
  if (HasBranch(tHi, "hiBin")) { tHi->SetBranchStatus("hiBin", 1); tHi->SetBranchAddress("hiBin", &hiBin); }

  // -------- Generator event weight (MC only) --------
  // hiEvtAnalyzer/HiTree::weight is the per-event generator weight. Applied to
  // every Fill below for MC; data is filled with weight 1.
  Float_t    genWeight     = 1.f;
  const bool has_genWeight = isMC && HasBranch(tHi, "weight");
  if (has_genWeight) { tHi->SetBranchStatus("weight", 1); tHi->SetBranchAddress("weight", &genWeight); }
  else if (isMC)     std::cout << "[WARN] skim_Wmu: MC sample but no 'weight' branch on HiTree; filling unweighted.\n";

  // -------- PF tree (for MET) --------
  Int_t nPF = 0;
  std::vector<int>   *pfId  = nullptr;
  std::vector<float> *pfPt  = nullptr;
  std::vector<float> *pfPhi = nullptr;
  tPF->SetBranchStatus("*", 0);
  tPF->SetBranchStatus("nPF", 1);
  tPF->SetBranchStatus("pfId", 1);
  tPF->SetBranchStatus("pfPt", 1);
  tPF->SetBranchStatus("pfPhi", 1);
  tPF->SetBranchAddress("nPF",  &nPF);
  tPF->SetBranchAddress("pfId", &pfId);
  tPF->SetBranchAddress("pfPt", &pfPt);
  tPF->SetBranchAddress("pfPhi", &pfPhi);

  // -------- Gen branches (MC only) --------
  std::vector<float> *genPt = nullptr, *genEta = nullptr, *genPhi = nullptr;
  std::vector<int>   *genChg = nullptr, *genPdg = nullptr;
  std::vector<int>   *genNMothers = nullptr;
  std::vector<std::vector<int>> *genMotherIdx = nullptr;

  if (isMC)
  {
    tGen->SetBranchStatus("*", 0);
    tGen->SetBranchStatus("pt", 1);
    tGen->SetBranchStatus("eta", 1);
    tGen->SetBranchStatus("phi", 1);
    tGen->SetBranchStatus("chg", 1);
    tGen->SetBranchStatus("pdg", 1);
    tGen->SetBranchStatus("nMothers", 1);
    tGen->SetBranchStatus("motherIdx", 1);

    tGen->SetBranchAddress("pt",        &genPt);
    tGen->SetBranchAddress("eta",       &genEta);
    tGen->SetBranchAddress("phi",       &genPhi);
    tGen->SetBranchAddress("chg",       &genChg);
    tGen->SetBranchAddress("pdg",       &genPdg);
    tGen->SetBranchAddress("nMothers",  &genNMothers);
    tGen->SetBranchAddress("motherIdx", &genMotherIdx);
  }

  // -------- Output histograms --------
  TH1D *h_met_Wp[kNY], *h_met_Wm[kNY], *h_mt_Wp[kNY], *h_mt_Wm[kNY];
  TH1D *h_met_Wp_FB[kNY], *h_met_Wm_FB[kNY], *h_mt_Wp_FB[kNY], *h_mt_Wm_FB[kNY];

  for (int b = 0; b < kNY; ++b)
  {
    const double y1 = kYEdges[b],     y2 = kYEdges[b + 1];
    const double y1FB = kYEdgesFB[b], y2FB = kYEdgesFB[b + 1];

    h_met_Wp[b] = new TH1D(Form("h_met_Wp_y%d", b),
        Form("W+ PF MET;MET [GeV];Events (%.2f<y<%.2f)", y1, y2), 60, 0, 120);
    h_met_Wm[b] = new TH1D(Form("h_met_Wm_y%d", b),
        Form("W- PF MET;MET [GeV];Events (%.2f<y<%.2f)", y1, y2), 60, 0, 120);
    h_mt_Wp[b]  = new TH1D(Form("h_mt_Wp_y%d", b),
        Form("W+ m_{T};m_{T} [GeV];Events (%.2f<y<%.2f)", y1, y2), 80, 0, 200);
    h_mt_Wm[b]  = new TH1D(Form("h_mt_Wm_y%d", b),
        Form("W- m_{T};m_{T} [GeV];Events (%.2f<y<%.2f)", y1, y2), 80, 0, 200);

    h_met_Wp_FB[b] = new TH1D(Form("h_met_Wp_y%d_FB", b),
        Form("W+ PF MET;MET [GeV];Events (%.2f<y<%.2f)", y1FB, y2FB), 60, 0, 120);
    h_met_Wm_FB[b] = new TH1D(Form("h_met_Wm_y%d_FB", b),
        Form("W- PF MET;MET [GeV];Events (%.2f<y<%.2f)", y1FB, y2FB), 60, 0, 120);
    h_mt_Wp_FB[b]  = new TH1D(Form("h_mt_Wp_y%d_FB", b),
        Form("W+ m_{T};m_{T} [GeV];Events (%.2f<y<%.2f)", y1FB, y2FB), 80, 0, 200);
    h_mt_Wm_FB[b]  = new TH1D(Form("h_mt_Wm_y%d_FB", b),
        Form("W- m_{T};m_{T} [GeV];Events (%.2f<y<%.2f)", y1FB, y2FB), 80, 0, 200);

    h_met_Wp[b]->Sumw2();    h_met_Wm[b]->Sumw2();
    h_mt_Wp[b]->Sumw2();     h_mt_Wm[b]->Sumw2();
    h_met_Wp_FB[b]->Sumw2(); h_met_Wm_FB[b]->Sumw2();
    h_mt_Wp_FB[b]->Sumw2();  h_mt_Wm_FB[b]->Sumw2();
  }

  // QCD sideband (anti-iso) histograms
  static const int NISO = 3;
  double isoEdges[NISO + 1] = {0.1, 0.4, 0.7, 1.0};

  TH1D *h_met_iso_muPlus[NISO], *h_met_iso_muMinus[NISO];
  for (int b = 0; b < NISO; ++b)
  {
    h_met_iso_muPlus[b]  = new TH1D(Form("h_met_iso_muPlus_bin%d", b),
        Form("MET, #mu^{+}, %.1f < iso < %.1f;MET [GeV];Events", isoEdges[b], isoEdges[b + 1]),
        60, 0, 120);
    h_met_iso_muMinus[b] = new TH1D(Form("h_met_iso_muMinus_bin%d", b),
        Form("MET, #mu^{-}, %.1f < iso < %.1f;MET [GeV];Events", isoEdges[b], isoEdges[b + 1]),
        60, 0, 120);
    h_met_iso_muPlus[b]->Sumw2();
    h_met_iso_muMinus[b]->Sumw2();
  }

  auto FindIsoBin = [&](double iso) -> int
  {
    for (int b = 0; b < NISO; ++b)
      if (iso >= isoEdges[b] && iso < isoEdges[b + 1]) return b;
    return -1;
  };

  // -------- ABCD inputs: 2D (relIso, MET) and (relIso, MT) per charge --------
  // Filled for EVERY event passing the full W selection EXCEPT the isolation
  // cut, so the relIso axis spans the signal (iso pass) and the QCD-enriched
  // anti-iso sideband (iso fail). The ABCD regions are defined DOWNSTREAM by
  // projecting these in correction/qcd_abcd.C, so the iso/MET boundaries can be
  // retuned without re-skimming (same philosophy as the recoil 2D histos).
  //   x = relIso : see the kNIsoAB note below (0.005-wide bins; every boundary
  //               used downstream must be a bin EDGE)
  //   y = MET    : 0..120 GeV / 2.0  (matches h_met_* binning)
  //   y = m_T    : 0..200 GeV / 2.5  (matches h_mt_*  binning)
  // MC is filled with the per-event gen weight w (data: w=1); the absolute
  // per-sample k_s is applied later in the ABCD macro (EWK subtraction).
  // 200 bins over [0,1] => 0.005 wide, so EVERY relIso boundary used
  // downstream lands on a bin EDGE: the electron cut 0.095 (= 19*0.005),
  // the muon 0.15, and the anti-iso window edges 0.20/0.30/0.60/0.65.
  // With the old 0.01 bins, 0.095 fell mid-bin and qcd_abcd.C's
  // FindBin(isoCut-eps) projection silently integrated relIso < 0.10 --
  // 170 electron events that FAIL the analysis cut leaked into the
  // ABCD iso-pass region, inflating T and every electron QCD template by
  // ~4% (found 2026-07-30; the muon was exact all along).
  const int    kNIsoAB  = 200;
  const double kIsoLoAB = 0.0, kIsoHiAB = 1.0;
  TH2D *h_iso_met_muPlus  = new TH2D("h_iso_met_muPlus",
      "#mu^{+} ABCD;relIso;PF MET [GeV]", kNIsoAB, kIsoLoAB, kIsoHiAB, 60, 0, 120);
  TH2D *h_iso_met_muMinus = new TH2D("h_iso_met_muMinus",
      "#mu^{-} ABCD;relIso;PF MET [GeV]", kNIsoAB, kIsoLoAB, kIsoHiAB, 60, 0, 120);
  TH2D *h_iso_mt_muPlus   = new TH2D("h_iso_mt_muPlus",
      "#mu^{+} ABCD;relIso;m_{T} [GeV]",  kNIsoAB, kIsoLoAB, kIsoHiAB, 80, 0, 200);
  TH2D *h_iso_mt_muMinus  = new TH2D("h_iso_mt_muMinus",
      "#mu^{-} ABCD;relIso;m_{T} [GeV]",  kNIsoAB, kIsoLoAB, kIsoHiAB, 80, 0, 200);
  //   y = lepton pT : 0..100 GeV / 2.0 (matches h_lepPt / h_leppt_* binning).
  // Third ABCD plane (2026-07-29): provides the anti-iso lepton-pT SHAPE for
  // the QCD template in the lepton-pT stacks (normalization still comes from
  // the relIso x MET counting -- relIso has pT in its denominator, so this
  // plane is NOT used for the 2x2 factorization itself).
  TH2D *h_iso_pt_muPlus   = new TH2D("h_iso_pt_muPlus",
      "#mu^{+} ABCD;relIso;p_{T}^{#mu} [GeV]", kNIsoAB, kIsoLoAB, kIsoHiAB, 50, 0, 100);
  TH2D *h_iso_pt_muMinus  = new TH2D("h_iso_pt_muMinus",
      "#mu^{-} ABCD;relIso;p_{T}^{#mu} [GeV]", kNIsoAB, kIsoLoAB, kIsoHiAB, 50, 0, 100);
  // Same plane WITH the discriminant selection's m_T cut: supplies the
  // anti-iso lepton-pT shape for the QCD template of the pT+m_T>40 stacks.
  TH2D *h_iso_pt_mt40_muPlus  = new TH2D("h_iso_pt_mt40_muPlus",
      Form("#mu^{+} ABCD, m_{T}>%.0f;relIso;p_{T}^{#mu} [GeV]", mtCutForPtDisc), kNIsoAB, kIsoLoAB, kIsoHiAB, 50, 0, 100);
  TH2D *h_iso_pt_mt40_muMinus = new TH2D("h_iso_pt_mt40_muMinus",
      Form("#mu^{-} ABCD, m_{T}>%.0f;relIso;p_{T}^{#mu} [GeV]", mtCutForPtDisc), kNIsoAB, kIsoLoAB, kIsoHiAB, 50, 0, 100);
  h_iso_pt_mt40_muPlus->Sumw2(); h_iso_pt_mt40_muMinus->Sumw2();
  h_iso_met_muPlus->Sumw2();  h_iso_met_muMinus->Sumw2();
  h_iso_mt_muPlus->Sumw2();   h_iso_mt_muMinus->Sumw2();
  h_iso_pt_muPlus->Sumw2();   h_iso_pt_muMinus->Sumw2();

  // -------- (lepton pT, m_T) joint 2Ds for the pT x mT cut scan --------
  // For correction/ptmt_scan.C (2026-07-30): scanning a (pT, mT) cut PAIR
  // needs the JOINT distribution, which the 1D h_leppt_*/h_mt_* cannot give.
  // Isolation stays at the nominal cut (this is NOT an iso scan):
  //   h_pt_mt_*         -- iso-pass events, no MET condition. W MC -> the
  //                        signal model S(pTcut, mTcut); data -> full sample.
  //   h_pt_mt_antiiso_* -- events in the ANTI-ISO sideband (relIso in
  //                        [0.30, 1.0), the qcd_abcd.C window), with their
  //                        ACTUAL MET/mT: the QCD background model of the
  //                        scan. Isolation is ~independent of (pT, mT), so
  //                        unlike a MET-conditioned proxy this does NOT
  //                        sculpt the mT axis (a MET<5 proxy caps mT at
  //                        ~2*sqrt(pT*MET) ~ 30 GeV by construction --
  //                        that earlier convention was removed 2026-07-30).
  // pT 1 GeV bins (integer-GeV cuts land on edges), mT 2.5 GeV (matches h_mt_*).
  // 2026-08-04: these four planes are filled down to leading-lepton
  // pT > kPtScanFloor (20 GeV, skim_common.h) instead of the nominal 25, so
  // ptmt_scan.C can probe LOWERING the pT cut. They are the ONLY histograms
  // with the relaxed floor -- everything else keeps pT > 25 (passPtNominal).
  const double antiIsoLo = 0.30, antiIsoHi = 1.00; // = qcd_abcd.C isoFail window (mu)
  TH2D *h_pt_mt_muPlus  = new TH2D("h_pt_mt_muPlus",
      "#mu^{+}, iso-pass;p_{T}^{#mu} [GeV];m_{T} [GeV]", 100, 0, 100, 80, 0, 200);
  TH2D *h_pt_mt_muMinus = new TH2D("h_pt_mt_muMinus",
      "#mu^{-}, iso-pass;p_{T}^{#mu} [GeV];m_{T} [GeV]", 100, 0, 100, 80, 0, 200);
  TH2D *h_pt_mt_antiiso_muPlus  = new TH2D("h_pt_mt_antiiso_muPlus",
      "#mu^{+}, anti-iso [0.30,1.0);p_{T}^{#mu} [GeV];m_{T} [GeV]", 100, 0, 100, 80, 0, 200);
  TH2D *h_pt_mt_antiiso_muMinus = new TH2D("h_pt_mt_antiiso_muMinus",
      "#mu^{-}, anti-iso [0.30,1.0);p_{T}^{#mu} [GeV];m_{T} [GeV]", 100, 0, 100, 80, 0, 200);
  h_pt_mt_muPlus->Sumw2();         h_pt_mt_muMinus->Sumw2();
  h_pt_mt_antiiso_muPlus->Sumw2(); h_pt_mt_antiiso_muMinus->Sumw2();

  // -------- Leading-lepton pT: per-charge, per-y (lab + FB) --------
  // pT twins of the h_met_*/h_mt_* rapidity-binned histos, for the lepton-pT
  // data/MC stacks in plotting/mtandmet.C (plots only -- NOT written into the
  // Combine input) and a possible future pT-discriminant fit. 2 GeV bins,
  // matching h_lepPt.
  //
  // The `_mt40` twins carry the SAME distributions with the m_T > kMtCutPt cut
  // applied (2026-07-30): the lepton-pT-DISCRIMINANT selection is pT > 25 AND
  // m_T > 40, the m_T cut replacing the QCD-suppression role that the MET
  // shape plays in the MET-discriminant fit. The cut cannot be applied
  // downstream (the pT histos have no m_T axis), hence the separate set.
  TH1D *h_leppt_Wp[kNY],    *h_leppt_Wm[kNY];
  TH1D *h_leppt_Wp_FB[kNY], *h_leppt_Wm_FB[kNY];
  TH1D *h_leppt_mt40_Wp[kNY],    *h_leppt_mt40_Wm[kNY];
  TH1D *h_leppt_mt40_Wp_FB[kNY], *h_leppt_mt40_Wm_FB[kNY];
  for (int b = 0; b < kNY; ++b)
  {
    const double y1 = kYEdges[b],     y2 = kYEdges[b + 1];
    const double y1FB = kYEdgesFB[b], y2FB = kYEdgesFB[b + 1];

    h_leppt_Wp[b] = new TH1D(Form("h_leppt_Wp_y%d", b),
        Form("W+ lepton p_{T};p_{T}^{#mu} [GeV];Events (%.2f<y<%.2f)", y1, y2), 50, 0, 100);
    h_leppt_Wm[b] = new TH1D(Form("h_leppt_Wm_y%d", b),
        Form("W- lepton p_{T};p_{T}^{#mu} [GeV];Events (%.2f<y<%.2f)", y1, y2), 50, 0, 100);
    h_leppt_Wp_FB[b] = new TH1D(Form("h_leppt_Wp_y%d_FB", b),
        Form("W+ lepton p_{T};p_{T}^{#mu} [GeV];Events (%.2f<y<%.2f)", y1FB, y2FB), 50, 0, 100);
    h_leppt_Wm_FB[b] = new TH1D(Form("h_leppt_Wm_y%d_FB", b),
        Form("W- lepton p_{T};p_{T}^{#mu} [GeV];Events (%.2f<y<%.2f)", y1FB, y2FB), 50, 0, 100);

    h_leppt_mt40_Wp[b] = new TH1D(Form("h_leppt_mt40_Wp_y%d", b),
        Form("W+ lepton p_{T}, m_{T}>%.0f;p_{T}^{#mu} [GeV];Events (%.2f<y<%.2f)", mtCutForPtDisc, y1, y2), 50, 0, 100);
    h_leppt_mt40_Wm[b] = new TH1D(Form("h_leppt_mt40_Wm_y%d", b),
        Form("W- lepton p_{T}, m_{T}>%.0f;p_{T}^{#mu} [GeV];Events (%.2f<y<%.2f)", mtCutForPtDisc, y1, y2), 50, 0, 100);
    h_leppt_mt40_Wp_FB[b] = new TH1D(Form("h_leppt_mt40_Wp_y%d_FB", b),
        Form("W+ lepton p_{T}, m_{T}>%.0f;p_{T}^{#mu} [GeV];Events (%.2f<y<%.2f)", mtCutForPtDisc, y1FB, y2FB), 50, 0, 100);
    h_leppt_mt40_Wm_FB[b] = new TH1D(Form("h_leppt_mt40_Wm_y%d_FB", b),
        Form("W- lepton p_{T}, m_{T}>%.0f;p_{T}^{#mu} [GeV];Events (%.2f<y<%.2f)", mtCutForPtDisc, y1FB, y2FB), 50, 0, 100);

    h_leppt_Wp[b]->Sumw2();         h_leppt_Wm[b]->Sumw2();
    h_leppt_Wp_FB[b]->Sumw2();      h_leppt_Wm_FB[b]->Sumw2();
    h_leppt_mt40_Wp[b]->Sumw2();    h_leppt_mt40_Wm[b]->Sumw2();
    h_leppt_mt40_Wp_FB[b]->Sumw2(); h_leppt_mt40_Wm_FB[b]->Sumw2();
  }

  // -------- Leading-lepton kinematics (full W selection) --------
  // Data/MC shape check for the W channel (correction/dataMC_kinematics.C),
  // mirroring the Z-channel h_lep* histos (same names, so the macro reads both
  // W and Z files uniformly). Filled once per event with the LEADING muon after
  // the full selection (isolation + trigger match included), charges combined.
  // The `_mt40` twins additionally require m_T > kMtCutPt (the lepton-pT
  // discriminant selection) -- with QCD suppressed there, the data/MC shape
  // check becomes quantitative rather than QCD-dominated.
  TH1D *h_lepPt  = new TH1D("h_lepPt",  Form("%s; p_{T}^{#mu} [GeV]; Events", outPrefix.c_str()), 50,  0,  100);
  TH1D *h_lepEta = new TH1D("h_lepEta", Form("%s; #eta_{#mu}; Events",        outPrefix.c_str()), 50, -2.5, 2.5);
  TH1D *h_lepPhi = new TH1D("h_lepPhi", Form("%s; #phi_{#mu}; Events",        outPrefix.c_str()), 32, -TMath::Pi(), TMath::Pi());
  TH1D *h_lepPt_mt40  = new TH1D("h_lepPt_mt40",  Form("%s, m_{T}>%.0f; p_{T}^{#mu} [GeV]; Events", outPrefix.c_str(), mtCutForPtDisc), 50,  0,  100);
  TH1D *h_lepEta_mt40 = new TH1D("h_lepEta_mt40", Form("%s, m_{T}>%.0f; #eta_{#mu}; Events", outPrefix.c_str(), mtCutForPtDisc), 50, -2.5, 2.5);
  TH1D *h_lepPhi_mt40 = new TH1D("h_lepPhi_mt40", Form("%s, m_{T}>%.0f; #phi_{#mu}; Events", outPrefix.c_str(), mtCutForPtDisc), 32, -TMath::Pi(), TMath::Pi());
  h_lepPt->Sumw2(); h_lepEta->Sumw2(); h_lepPhi->Sumw2();
  h_lepPt_mt40->Sumw2(); h_lepEta_mt40->Sumw2(); h_lepPhi_mt40->Sumw2();

  // -------- Cutflow loop --------
  unsigned long long N[9] = {0};
  const Long64_t nEntries = tMu->GetEntries();
  std::cout << "Entries: " << nEntries << "\n";

  bool warnedEventFiltersOnce = false;
  bool warnedNoTrigMatchInfo  = false;

  for (Long64_t ie = 0; ie < nEntries; ++ie)
  {
    if (ie % 200000 == 0) std::cout << "Event " << ie << "/" << nEntries << "\n";

    tMu->GetEntry(ie);
    tHi->GetEntry(ie);
    tPF->GetEntry(ie);
    tHLT->GetEntry(ie);
    tHLTobj->GetEntry(ie);
    tEvent->GetEntry(ie);
    if (isMC) tGen->GetEntry(ie);

    const double w = has_genWeight ? (double)genWeight : 1.0;

    N[0]++;

    if (!muPt || !muEta || !muPhi || !muCharge) continue;

    // Precut relaxed to the (pT, mT)-scan floor (pT > 20) so the scan planes
    // h_pt_mt[_antiiso]_* reach below the nominal cut. The cutflow keeps its
    // nominal meaning: N[1..5] count only events the pT > 25 precut passes.
    const bool hasPF25 = ExistsPFMuonPt25_W(muPt25, nMu, muPt, has_muIsPF, muIsPF);
    if (!ExistsPFMuonPt25_W(kPtScanFloor, nMu, muPt, has_muIsPF, muIsPF)) continue;
    if (hasPF25) N[1]++;

    if (!PassEventSelection_pO(warnedEventFiltersOnce,
                               has_pprimaryVertexFilter,        pprimaryVertexFilter,
                               has_pclusterCompatibilityFilter, pclusterCompatibilityFilter))
      continue;
    if (applyVz && TMath::Abs(vz) > vzMax) continue;
    if (hasPF25) N[2]++;

    if (has_hlt && !TriggerFired(HLT_OxyL1SingleMuOpen_v1)) continue;
    if (hasPF25) N[3]++;

    if (!PassDYVeto_Wmu(dyMuPtMin, isoMax, dyMassMin, dyMassMax,
                        nMu, muPt, muEta, muPhi, muCharge,
                        has_muIDTight, muIDTight,
                        has_muIsPF, muIsPF,
                        muPFChIso, muPFNeuIso, muPFPhoIso, muPFPUIso))
      continue;
    if (hasPF25) N[4]++;

    if (!ExistsTightMuon_W(nMu, has_muIDTight, muIDTight)) continue;
    if (hasPF25) N[5]++;

    const int iLead = FindLeadingMuon_TightPF(nMu, muPt, muEta, muPhi,
                                              has_muIDTight, muIDTight,
                                              has_muIsPF,    muIsPF);
    if (iLead < 0) continue;

    if (isMC)
    {
      bool passWmuMatch = PassGenRecoMatchingWithAncestor(
          iLead,
          muPt, muEta, muPhi, muCharge,
          genPt, genEta, genPhi, genChg,
          genPdg, genMotherIdx,
          13, // matched gen particle must be muon
          24  // ancestor must be W
      );
      if (!passWmuMatch)
      {
        // cout << "We see a gen-reco matching failed case" << endl;
        // continue;
      }
    }

    if (muPt->at(iLead) <= kPtScanFloor) continue; // scan floor; nominal pT cut via passPtNominal
    if (std::abs(muEta->at(iLead)) > muEtaMax) continue;
    const bool passPtNominal = (muPt->at(iLead) > muPt25);
    if (passPtNominal) N[6]++;

    const double isoLead = RelIsoPF(iLead, muPt, muPFChIso, muPFNeuIso, muPFPhoIso, muPFPUIso);
    const int    isoBin  = FindIsoBin(isoLead);
    const bool   isQCDSideband   = (isoBin >= 0);
    const bool   passIsoNominal  = (isoLead < isoMax);

    if (passIsoNominal && passPtNominal) N[7]++;

    const bool passMatch = PassLeadingLeptonTrigMatch(
        trigMatchDR, iLead, muEta, muPhi,
        has_trgObjPt, trgObjPt, has_trgObjEta, trgObjEta,
        has_trgObjPhi, trgObjPhi,
        warnedNoTrigMatchInfo, "muon");
    if (!passMatch) continue;
    if (passIsoNominal && passPtNominal) N[8]++;

    // -------- Fill final distributions (after step 8) --------
    TVector2 metv = ComputePFMET(pfId, pfPt, pfPhi);
    const double met = metv.Mod();
    const double mt  = TransverseMass(muPt->at(iLead), muPhi->at(iLead), metv);

    if (isQCDSideband && passPtNominal)
    {
      if      (muCharge->at(iLead) > 0) h_met_iso_muPlus [isoBin]->Fill(met, w);
      else if (muCharge->at(iLead) < 0) h_met_iso_muMinus[isoBin]->Fill(met, w);
    }

    // ABCD: fill the full-iso-range (relIso, MET) and (relIso, MT) planes.
    // Unconditional on the isolation cut -- iso pass AND fail both populated.
    // The h_iso_* ABCD planes keep the NOMINAL pT > 25 selection; only the
    // (pT, mT) scan planes h_pt_mt[_antiiso]_* go down to kPtScanFloor.
    const bool passMtCut = (mt > mtCutForPtDisc); // lepton-pT-discriminant m_T cut
    if (muCharge->at(iLead) > 0)
    {
      if (passPtNominal)
      {
        h_iso_met_muPlus->Fill(isoLead, met, w);
        h_iso_mt_muPlus ->Fill(isoLead, mt,  w);
        h_iso_pt_muPlus ->Fill(isoLead, muPt->at(iLead), w);
        if (passMtCut) h_iso_pt_mt40_muPlus->Fill(isoLead, muPt->at(iLead), w);
      }
      if (isoLead >= antiIsoLo && isoLead < antiIsoHi)
        h_pt_mt_antiiso_muPlus->Fill(muPt->at(iLead), mt, w);
    }
    else if (muCharge->at(iLead) < 0)
    {
      if (passPtNominal)
      {
        h_iso_met_muMinus->Fill(isoLead, met, w);
        h_iso_mt_muMinus ->Fill(isoLead, mt,  w);
        h_iso_pt_muMinus ->Fill(isoLead, muPt->at(iLead), w);
        if (passMtCut) h_iso_pt_mt40_muMinus->Fill(isoLead, muPt->at(iLead), w);
      }
      if (isoLead >= antiIsoLo && isoLead < antiIsoHi)
        h_pt_mt_antiiso_muMinus->Fill(muPt->at(iLead), mt, w);
    }

    if (!passIsoNominal) continue;

    // pT x mT scan signal plane (iso-pass; filled down to the pT > 20 scan floor)
    if      (muCharge->at(iLead) > 0) h_pt_mt_muPlus ->Fill(muPt->at(iLead), mt, w);
    else if (muCharge->at(iLead) < 0) h_pt_mt_muMinus->Fill(muPt->at(iLead), mt, w);

    if (!passPtNominal) continue; // everything below keeps the nominal pT > 25 cut

    h_lepPt ->Fill(muPt->at(iLead),  w);
    h_lepEta->Fill(muEta->at(iLead), w);
    h_lepPhi->Fill(muPhi->at(iLead), w);
    if (passMtCut)
    {
      h_lepPt_mt40 ->Fill(muPt->at(iLead),  w);
      h_lepEta_mt40->Fill(muEta->at(iLead), w);
      h_lepPhi_mt40->Fill(muPhi->at(iLead), w);
    }

    const int q = muCharge->at(iLead);
    const bool isWp = (q > 0), isWm = (q < 0);

    // [NOTICE] Sign flip: p-going (-Z) defined as forward
    const double y       = -muEta->at(iLead);
    const int    ybin    = FindYBin(y, kYEdges,   kNY);
    const int    ybin_FB = FindYBin(y, kYEdgesFB, kNY);

    if (ybin >= 0)
    {
      if      (isWp) { h_met_Wp[ybin]->Fill(met, w); h_mt_Wp[ybin]->Fill(mt, w); h_leppt_Wp[ybin]->Fill(muPt->at(iLead), w); }
      else if (isWm) { h_met_Wm[ybin]->Fill(met, w); h_mt_Wm[ybin]->Fill(mt, w); h_leppt_Wm[ybin]->Fill(muPt->at(iLead), w); }
      if (passMtCut)
      {
        if      (isWp) h_leppt_mt40_Wp[ybin]->Fill(muPt->at(iLead), w);
        else if (isWm) h_leppt_mt40_Wm[ybin]->Fill(muPt->at(iLead), w);
      }
    }
    if (ybin_FB >= 0)
    {
      if      (isWp) { h_met_Wp_FB[ybin_FB]->Fill(met, w); h_mt_Wp_FB[ybin_FB]->Fill(mt, w); h_leppt_Wp_FB[ybin_FB]->Fill(muPt->at(iLead), w); }
      else if (isWm) { h_met_Wm_FB[ybin_FB]->Fill(met, w); h_mt_Wm_FB[ybin_FB]->Fill(mt, w); h_leppt_Wm_FB[ybin_FB]->Fill(muPt->at(iLead), w); }
      if (passMtCut)
      {
        if      (isWp) h_leppt_mt40_Wp_FB[ybin_FB]->Fill(muPt->at(iLead), w);
        else if (isWm) h_leppt_mt40_Wm_FB[ybin_FB]->Fill(muPt->at(iLead), w);
      }
    }
  }

  // -------- Cutflow + outputs --------
  PrintCutflow(std::cout, N);
  WriteCutflowTxt(outPrefix, isMC, N);

  gSystem->mkdir("rootfile", kTRUE); // ensure ./rootfile exists (fresh checkout)
  TFile *fout = new TFile(("./rootfile/" + outPrefix + "_hist.root").c_str(), "RECREATE");
  for (int i = 0; i < kNY; ++i)
  {
    h_met_Wp[i]->Write("", 2);
    h_met_Wm[i]->Write("", 2);
    h_mt_Wp [i]->Write("", 2);
    h_mt_Wm [i]->Write("", 2);
    h_met_Wp_FB[i]->Write("", 2);
    h_met_Wm_FB[i]->Write("", 2);
    h_mt_Wp_FB [i]->Write("", 2);
    h_mt_Wm_FB [i]->Write("", 2);
    h_leppt_Wp[i]->Write("", 2);
    h_leppt_Wm[i]->Write("", 2);
    h_leppt_Wp_FB[i]->Write("", 2);
    h_leppt_Wm_FB[i]->Write("", 2);
    h_leppt_mt40_Wp[i]->Write("", 2);
    h_leppt_mt40_Wm[i]->Write("", 2);
    h_leppt_mt40_Wp_FB[i]->Write("", 2);
    h_leppt_mt40_Wm_FB[i]->Write("", 2);
  }
  for (int b = 0; b < NISO; ++b)
  {
    h_met_iso_muPlus [b]->Write("", TObject::kOverwrite);
    h_met_iso_muMinus[b]->Write("", TObject::kOverwrite);
  }
  // ABCD 2D planes (relIso vs MET / m_T / lepton pT), per charge
  h_iso_met_muPlus ->Write("", TObject::kOverwrite);
  h_iso_met_muMinus->Write("", TObject::kOverwrite);
  h_iso_mt_muPlus  ->Write("", TObject::kOverwrite);
  h_iso_mt_muMinus ->Write("", TObject::kOverwrite);
  h_iso_pt_muPlus  ->Write("", TObject::kOverwrite);
  h_iso_pt_muMinus ->Write("", TObject::kOverwrite);
  h_iso_pt_mt40_muPlus ->Write("", TObject::kOverwrite);
  h_iso_pt_mt40_muMinus->Write("", TObject::kOverwrite);
  // pT x mT scan inputs (correction/ptmt_scan.C)
  h_pt_mt_muPlus         ->Write("", TObject::kOverwrite);
  h_pt_mt_muMinus        ->Write("", TObject::kOverwrite);
  h_pt_mt_antiiso_muPlus ->Write("", TObject::kOverwrite);
  h_pt_mt_antiiso_muMinus->Write("", TObject::kOverwrite);
  // Leading-lepton kinematics (full selection; _mt40 = + the m_T>40 cut)
  h_lepPt ->Write("", TObject::kOverwrite);
  h_lepEta->Write("", TObject::kOverwrite);
  h_lepPhi->Write("", TObject::kOverwrite);
  h_lepPt_mt40 ->Write("", TObject::kOverwrite);
  h_lepEta_mt40->Write("", TObject::kOverwrite);
  h_lepPhi_mt40->Write("", TObject::kOverwrite);
  fout->Close();

  std::cout << "[INFO] Wrote outputs with prefix: " << outPrefix << "\n";
  return 0;
}


// ============================================================
// W -> e+nu
// Legacy: skim/legacy/DrawWToElecNu_PFMet.C
// ============================================================

int skim_Wel(const char *fname, SampleType sample)
{
  const bool isMC = IsMC(sample);

  // -------- Config (kept identical to legacy WConfig for the electron channel) --------
  const double elePt25     = 25.0;
  const double dyElePtMin  = 10.0;
  const double isoMax      = 0.095;
  // m_T cut of the lepton-pT-DISCRIMINANT selection (pT > 25 && m_T > 40);
  // applies ONLY to the `_mt40` histogram set (see skim_Wmu for the rationale).
  // NAMING CONTRACT: the literal "40" in the persisted `_mt40` histogram names
  // asserts this value. If you change it, RENAME those histograms too (and the
  // readers: correction/qcd_abcd.C runPtCharge tag, plotting/mtandmet.C
  // varStem[kVarMt40], correction/dataMC_kinematics.C's _mt40 channels) --
  // otherwise the stored names silently claim a cut that was not applied.
  const double mtCutForPtDisc     = 40.0;
  const double dyMassMin   = 80.0;
  const double dyMassMax   = 110.0;
  const double eleEtaMax   = 2.4;
  // 0.1 -> 0.4 (2026-08-12): the saved trigger objects are L1-granularity
  // candidates whose phi is the muon's azimuth AT THE MUON STATION, while
  // {mu,ele}Phi is the vertex direction. The two differ by the solenoid
  // bending, Dphi = c + q*k/pT with k = 2.55 rad*GeV and a charge-blind
  // c = +0.006 rad (half of the 2pi/576 muGMT phi cell). At pT 25-30 in the
  // barrel that is ~0.10 rad, i.e. ON the old window -- and because c adds to
  // the mu+ bending and cancels part of the mu- one, the old cut kept only
  // 48% of barrel mu+ vs 91% of mu- there (a ~10% CHARGE-DEPENDENT loss riding
  // straight on the charge-asymmetry observable). 0.4 is comfortably beyond
  // the largest |Dphi| (~0.10) and restores 100%/100% for both charges.
  const double trigMatchDR = 0.4;
  const bool   applyVz     = true;
  const double vzMax       = 15.0;

  std::string outPrefix = "WToElecNu_pO_PFMet";

  gStyle->SetOptStat(0);

  if (isMC)
  {
    const auto info = ResolveMCSample(sample, "ele");
    if (!info.fname.empty()) fname = Form("%s", info.fname.c_str());
    outPrefix += info.outSuffix;
  }

  std::cout << "[INPUT] " << fname << std::endl;
  TFile *f = TFile::Open(fname);
  if (!f || f->IsZombie()) { std::cerr << "ERROR opening file\n"; return 1; }

  TTree *tEle    = (TTree *)f->Get("ggHiNtuplizer/EventTree");
  TTree *tHi     = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  TTree *tPF     = (TTree *)f->Get("particleFlowAnalyser/pftree");
  TTree *tHLT    = (TTree *)f->Get("hltanalysis/HltTree");
  TTree *tHLTobj = (TTree *)f->Get("hltobject/HLT_OxyL1SingleEG10_v");
  TTree *tEvent  = (TTree *)f->Get("skimanalysis/HltTree");
  if (!tEvent)
  {
    std::cerr << "[FATAL] skim_Wel: missing TTree skimanalysis/HltTree in " << fname << "\n";
    f->Close();
    return 2;
  }
  TTree *tGen    = nullptr;
  if (isMC) tGen = (TTree *)f->Get("HiGenParticleAna/hi");

  if (!tEle || !tHi || !tPF || !tHLT || !tHLTobj)
  {
    std::cerr << "ERROR: missing one of trees: ggHiNtuplizer/EventTree, hiEvtAnalyzer/HiTree, particleFlowAnalyser/pftree\n";
    return 1;
  }

  // -------- Electron kinematics --------
  Int_t nEle = 0;
  std::vector<float> *elePt = nullptr, *eleEta = nullptr, *elePhi = nullptr;
  std::vector<int>   *eleCharge = nullptr;
  std::vector<float> *elePFChIso = nullptr, *elePFNeuIso = nullptr, *elePFPhoIso = nullptr;
  std::vector<float> *elePFPUIso = nullptr;

  tEle->SetBranchStatus("*", 0);
  tEle->SetBranchStatus("nEle", 1);
  tEle->SetBranchStatus("elePt", 1);
  tEle->SetBranchStatus("eleEta", 1);
  tEle->SetBranchStatus("elePhi", 1);
  tEle->SetBranchStatus("eleCharge", 1);
  tEle->SetBranchAddress("nEle",      &nEle);
  tEle->SetBranchAddress("elePt",     &elePt);
  tEle->SetBranchAddress("eleEta",    &eleEta);
  tEle->SetBranchAddress("elePhi",    &elePhi);
  tEle->SetBranchAddress("eleCharge", &eleCharge);

  // -------- Electron IDs (cut-based + MVA) --------
  std::vector<int> *eleMVAIdWP80 = nullptr, *eleMVAIdWP85 = nullptr;
  std::vector<int> *eleMVAIdWP90 = nullptr, *eleMVAIdWP95 = nullptr;
  std::vector<int> *eleCutIdWP70 = nullptr, *eleCutIdWP80 = nullptr;
  std::vector<int> *eleCutIdWP90 = nullptr, *eleCutIdWP95 = nullptr;

  tEle->SetBranchStatus("eleMVAIdWP80", 1); tEle->SetBranchStatus("eleMVAIdWP85", 1);
  tEle->SetBranchStatus("eleMVAIdWP90", 1); tEle->SetBranchStatus("eleMVAIdWP95", 1);
  tEle->SetBranchStatus("eleCutIdWP70", 1); tEle->SetBranchStatus("eleCutIdWP80", 1);
  tEle->SetBranchStatus("eleCutIdWP90", 1); tEle->SetBranchStatus("eleCutIdWP95", 1);
  tEle->SetBranchAddress("eleMVAIdWP80", &eleMVAIdWP80);
  tEle->SetBranchAddress("eleMVAIdWP85", &eleMVAIdWP85);
  tEle->SetBranchAddress("eleMVAIdWP90", &eleMVAIdWP90);
  tEle->SetBranchAddress("eleMVAIdWP95", &eleMVAIdWP95);
  tEle->SetBranchAddress("eleCutIdWP70", &eleCutIdWP70);
  tEle->SetBranchAddress("eleCutIdWP80", &eleCutIdWP80);
  tEle->SetBranchAddress("eleCutIdWP90", &eleCutIdWP90);
  tEle->SetBranchAddress("eleCutIdWP95", &eleCutIdWP95);

  // -------- Electron iso (PF + WP) --------
  if (HasBranch(tEle, "elePFChIso"))  tEle->SetBranchStatus("elePFChIso", 1);
  if (HasBranch(tEle, "elePFNeuIso")) tEle->SetBranchStatus("elePFNeuIso", 1);
  if (HasBranch(tEle, "elePFPhoIso")) tEle->SetBranchStatus("elePFPhoIso", 1);
  if (HasBranch(tEle, "elePFPUIso"))  tEle->SetBranchStatus("elePFPUIso", 1);

  const bool has_elePFChIso  = HasBranch(tEle, "elePFChIso");
  const bool has_elePFNeuIso = HasBranch(tEle, "elePFNeuIso");
  const bool has_elePFPhoIso = HasBranch(tEle, "elePFPhoIso");
  const bool has_elePFPUIso  = HasBranch(tEle, "elePFPUIso");
  if (has_elePFChIso)  tEle->SetBranchAddress("elePFChIso",  &elePFChIso);
  if (has_elePFNeuIso) tEle->SetBranchAddress("elePFNeuIso", &elePFNeuIso);
  if (has_elePFPhoIso) tEle->SetBranchAddress("elePFPhoIso", &elePFPhoIso);
  if (has_elePFPUIso)  tEle->SetBranchAddress("elePFPUIso",  &elePFPUIso);
  if (!has_elePFPUIso)
    std::cout << "[WARN] skim_Wel: no elePFPUIso branch; relIso falls back to uncorrected (no Delta-beta).\n";

  std::vector<int> *eleMVAIsoWP80 = nullptr, *eleMVAIsoWP85 = nullptr;
  std::vector<int> *eleMVAIsoWP90 = nullptr, *eleMVAIsoWP95 = nullptr;
  tEle->SetBranchStatus("eleMVAIsoWP80", 1); tEle->SetBranchStatus("eleMVAIsoWP85", 1);
  tEle->SetBranchStatus("eleMVAIsoWP90", 1); tEle->SetBranchStatus("eleMVAIsoWP95", 1);
  tEle->SetBranchAddress("eleMVAIsoWP80", &eleMVAIsoWP80);
  tEle->SetBranchAddress("eleMVAIsoWP85", &eleMVAIsoWP85);
  tEle->SetBranchAddress("eleMVAIsoWP90", &eleMVAIsoWP90);
  tEle->SetBranchAddress("eleMVAIsoWP95", &eleMVAIsoWP95);

  // -------- HLT bit (step 3) --------
  const std::string hltName = FindBranchContaining(tHLT, "HLT_OxyL1SingleEG10_v1");
  Int_t HLT_OxyL1SingleEG10_v1 = 0;
  bool  has_hlt = false;
  tHLT->SetBranchStatus("*", 0);
  if (!hltName.empty())
  {
    has_hlt = true;
    tHLT->SetBranchStatus(hltName.c_str(), 1);
    tHLT->SetBranchAddress(hltName.c_str(), &HLT_OxyL1SingleEG10_v1);
    std::cout << "[INFO] Using trigger bit branch: " << hltName << "\n";
  }
  else
  {
    std::cout << "[WARN] Could not find any branch containing 'HLT_OxyL1SingleEG10_v' in hltanalysis/HltTree.\n"
              << "       Step 3 (trigger fired) will be treated as PASS.\n";
  }

  // -------- HLT objects (step 8) --------
  std::vector<double> *trgObjPt = nullptr, *trgObjEta = nullptr;
  std::vector<double> *trgObjPhi = nullptr, *trgObjId = nullptr;

  const bool has_trgObjPt  = HasBranch(tHLTobj, "pt");
  const bool has_trgObjEta = HasBranch(tHLTobj, "eta");
  const bool has_trgObjPhi = HasBranch(tHLTobj, "phi");
  const bool has_trgObjId  = HasBranch(tHLTobj, "TriggerObjID");

  tHLTobj->SetBranchStatus("*", false);
  if (has_trgObjPt && has_trgObjEta && has_trgObjPhi && has_trgObjId)
  {
    tHLTobj->SetBranchStatus("pt", 1);
    tHLTobj->SetBranchStatus("eta", 1);
    tHLTobj->SetBranchStatus("phi", 1);
    tHLTobj->SetBranchStatus("TriggerObjID", 1);
    tHLTobj->SetBranchAddress("pt",           &trgObjPt);
    tHLTobj->SetBranchAddress("eta",          &trgObjEta);
    tHLTobj->SetBranchAddress("phi",          &trgObjPhi);
    tHLTobj->SetBranchAddress("TriggerObjID", &trgObjId);
  }
  else
  {
    std::cout << "[WARN] Could not find HLT related Object." << std::endl;
  }

  // -------- HiTree + event filters (step 2) --------
  Float_t vz    = 999.f;
  Int_t   hiBin = -999;
  Int_t   pprimaryVertexFilter        = 1;
  Int_t   pclusterCompatibilityFilter = 1;

  const bool has_pprimaryVertexFilter        = HasBranch(tEvent, "pprimaryVertexFilter");
  const bool has_pclusterCompatibilityFilter = HasBranch(tEvent, "pclusterCompatibilityFilter");

  tEvent->SetBranchStatus("*", 0);
  if (has_pprimaryVertexFilter && has_pclusterCompatibilityFilter)
  {
    tEvent->SetBranchStatus("pprimaryVertexFilter", 1);
    tEvent->SetBranchStatus("pclusterCompatibilityFilter", 1);
    tEvent->SetBranchAddress("pprimaryVertexFilter",        &pprimaryVertexFilter);
    tEvent->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);
  }

  tHi->SetBranchStatus("*", 0);
  if (HasBranch(tHi, "vz"))    { tHi->SetBranchStatus("vz", 1);    tHi->SetBranchAddress("vz",    &vz); }
  if (HasBranch(tHi, "hiBin")) { tHi->SetBranchStatus("hiBin", 1); tHi->SetBranchAddress("hiBin", &hiBin); }

  // -------- Generator event weight (MC only) --------
  Float_t    genWeight     = 1.f;
  const bool has_genWeight = isMC && HasBranch(tHi, "weight");
  if (has_genWeight) { tHi->SetBranchStatus("weight", 1); tHi->SetBranchAddress("weight", &genWeight); }
  else if (isMC)     std::cout << "[WARN] skim_Wel: MC sample but no 'weight' branch on HiTree; filling unweighted.\n";

  // -------- PF tree (MET) --------
  Int_t nPF = 0;
  std::vector<int>   *pfId  = nullptr;
  std::vector<float> *pfPt  = nullptr;
  std::vector<float> *pfPhi = nullptr;
  tPF->SetBranchStatus("*", 0);
  tPF->SetBranchStatus("nPF", 1);
  tPF->SetBranchStatus("pfId", 1);
  tPF->SetBranchStatus("pfPt", 1);
  tPF->SetBranchStatus("pfPhi", 1);
  tPF->SetBranchAddress("nPF",  &nPF);
  tPF->SetBranchAddress("pfId", &pfId);
  tPF->SetBranchAddress("pfPt", &pfPt);
  tPF->SetBranchAddress("pfPhi", &pfPhi);

  // -------- Gen (MC) --------
  std::vector<float> *genPt = nullptr, *genEta = nullptr, *genPhi = nullptr;
  std::vector<int>   *genChg = nullptr, *genPdg = nullptr;
  std::vector<std::vector<int>> *genMotherIdx = nullptr;

  if (isMC && !tGen)
  {
    std::cerr << "[FATAL] isMC=true but no gen tree found\n";
    return 1;
  }
  if (isMC)
  {
    gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");
    tGen->SetBranchStatus("*", 0);
    tGen->SetBranchStatus("pt", 1);
    tGen->SetBranchStatus("eta", 1);
    tGen->SetBranchStatus("phi", 1);
    tGen->SetBranchStatus("chg", 1);
    tGen->SetBranchStatus("pdg", 1);
    tGen->SetBranchStatus("motherIdx", 1);
    tGen->SetBranchAddress("pt",        &genPt);
    tGen->SetBranchAddress("eta",       &genEta);
    tGen->SetBranchAddress("phi",       &genPhi);
    tGen->SetBranchAddress("chg",       &genChg);
    tGen->SetBranchAddress("pdg",       &genPdg);
    tGen->SetBranchAddress("motherIdx", &genMotherIdx);
  }

  // -------- Output histograms --------
  TH1D *h_met_Wp[kNY], *h_met_Wm[kNY], *h_mt_Wp[kNY], *h_mt_Wm[kNY];
  TH1D *h_met_Wp_FB[kNY], *h_met_Wm_FB[kNY], *h_mt_Wp_FB[kNY], *h_mt_Wm_FB[kNY];

  for (int b = 0; b < kNY; ++b)
  {
    const double y1 = kYEdges[b],     y2 = kYEdges[b + 1];
    const double y1FB = kYEdgesFB[b], y2FB = kYEdgesFB[b + 1];

    h_met_Wp[b] = new TH1D(Form("h_met_Wp_y%d", b),
        Form("W+ PF MET;MET [GeV];Events (%.2f<y<%.2f)", y1, y2), 60, 0, 120);
    h_met_Wm[b] = new TH1D(Form("h_met_Wm_y%d", b),
        Form("W- PF MET;MET [GeV];Events (%.2f<y<%.2f)", y1, y2), 60, 0, 120);
    h_mt_Wp[b]  = new TH1D(Form("h_mt_Wp_y%d", b),
        Form("W+ m_{T};m_{T} [GeV];Events (%.2f<y<%.2f)", y1, y2), 80, 0, 200);
    h_mt_Wm[b]  = new TH1D(Form("h_mt_Wm_y%d", b),
        Form("W- m_{T};m_{T} [GeV];Events (%.2f<y<%.2f)", y1, y2), 80, 0, 200);

    h_met_Wp_FB[b] = new TH1D(Form("h_met_Wp_y%d_FB", b),
        Form("W+ PF MET;MET [GeV];Events (%.2f<y<%.2f)", y1FB, y2FB), 60, 0, 120);
    h_met_Wm_FB[b] = new TH1D(Form("h_met_Wm_y%d_FB", b),
        Form("W- PF MET;MET [GeV];Events (%.2f<y<%.2f)", y1FB, y2FB), 60, 0, 120);
    h_mt_Wp_FB[b]  = new TH1D(Form("h_mt_Wp_y%d_FB", b),
        Form("W+ m_{T};m_{T} [GeV];Events (%.2f<y<%.2f)", y1FB, y2FB), 80, 0, 200);
    h_mt_Wm_FB[b]  = new TH1D(Form("h_mt_Wm_y%d_FB", b),
        Form("W- m_{T};m_{T} [GeV];Events (%.2f<y<%.2f)", y1FB, y2FB), 80, 0, 200);

    h_met_Wp[b]->Sumw2();    h_met_Wm[b]->Sumw2();
    h_mt_Wp[b]->Sumw2();     h_mt_Wm[b]->Sumw2();
    h_met_Wp_FB[b]->Sumw2(); h_met_Wm_FB[b]->Sumw2();
    h_mt_Wp_FB[b]->Sumw2();  h_mt_Wm_FB[b]->Sumw2();
  }

  static const int NISO = 3;
  double isoEdges[NISO + 1] = {0.1, 0.4, 0.7, 1.0}; // tune for e if needed

  TH1D *h_met_iso_elePlus[NISO], *h_met_iso_eleMinus[NISO];
  for (int b = 0; b < NISO; ++b)
  {
    h_met_iso_elePlus[b]  = new TH1D(Form("h_met_iso_elePlus_bin%d", b),
        Form("MET, e^{+}, %.1f < iso < %.1f;MET [GeV];Events", isoEdges[b], isoEdges[b + 1]),
        60, 0, 120);
    h_met_iso_eleMinus[b] = new TH1D(Form("h_met_iso_eleMinus_bin%d", b),
        Form("MET, e^{-}, %.1f < iso < %.1f;MET [GeV];Events", isoEdges[b], isoEdges[b + 1]),
        60, 0, 120);
    h_met_iso_elePlus[b]->Sumw2();
    h_met_iso_eleMinus[b]->Sumw2();
  }

  auto FindIsoBin = [&](double iso) -> int
  {
    for (int b = 0; b < NISO; ++b)
      if (iso >= isoEdges[b] && iso < isoEdges[b + 1]) return b;
    return -1;
  };

  // -------- ABCD inputs: 2D (relIso, MET) and (relIso, MT) per charge --------
  // Electron twin of the skim_Wmu ABCD histos (consumed by correction/qcd_abcd.C
  // with isElec=true). Filled for every event passing the full W selection
  // EXCEPT the isolation cut, so the relIso axis spans iso-pass (signal) and the
  // anti-iso sideband. Regions are chosen downstream by projecting (no re-skim).
  //   x = relIso : see the kNIsoAB note below (0.005-wide bins so the electron
  //               signal cut 0.095 is exactly a bin edge)
  //   y = MET    : 0..120 GeV / 2.0  (matches h_met_* binning)
  //   y = m_T    : 0..200 GeV / 2.5  (matches h_mt_*  binning)
  // 200 bins over [0,1] => 0.005 wide, so EVERY relIso boundary used
  // downstream lands on a bin EDGE: the electron cut 0.095 (= 19*0.005),
  // the muon 0.15, and the anti-iso window edges 0.20/0.30/0.60/0.65.
  // With the old 0.01 bins, 0.095 fell mid-bin and qcd_abcd.C's
  // FindBin(isoCut-eps) projection silently integrated relIso < 0.10 --
  // 170 electron events that FAIL the analysis cut leaked into the
  // ABCD iso-pass region, inflating T and every electron QCD template by
  // ~4% (found 2026-07-30; the muon was exact all along).
  const int    kNIsoAB  = 200;
  const double kIsoLoAB = 0.0, kIsoHiAB = 1.0;
  TH2D *h_iso_met_elePlus  = new TH2D("h_iso_met_elePlus",
      "e^{+} ABCD;relIso;PF MET [GeV]", kNIsoAB, kIsoLoAB, kIsoHiAB, 60, 0, 120);
  TH2D *h_iso_met_eleMinus = new TH2D("h_iso_met_eleMinus",
      "e^{-} ABCD;relIso;PF MET [GeV]", kNIsoAB, kIsoLoAB, kIsoHiAB, 60, 0, 120);
  TH2D *h_iso_mt_elePlus   = new TH2D("h_iso_mt_elePlus",
      "e^{+} ABCD;relIso;m_{T} [GeV]",  kNIsoAB, kIsoLoAB, kIsoHiAB, 80, 0, 200);
  TH2D *h_iso_mt_eleMinus  = new TH2D("h_iso_mt_eleMinus",
      "e^{-} ABCD;relIso;m_{T} [GeV]",  kNIsoAB, kIsoLoAB, kIsoHiAB, 80, 0, 200);
  //   y = lepton pT : 0..100 GeV / 2.0 (matches h_lepPt / h_leppt_* binning).
  // Third ABCD plane (2026-07-29): anti-iso lepton-pT SHAPE for the QCD
  // template in the lepton-pT stacks (normalization from the relIso x MET
  // counting -- this plane is NOT used for the 2x2 factorization).
  TH2D *h_iso_pt_elePlus   = new TH2D("h_iso_pt_elePlus",
      "e^{+} ABCD;relIso;p_{T}^{e} [GeV]", kNIsoAB, kIsoLoAB, kIsoHiAB, 50, 0, 100);
  TH2D *h_iso_pt_eleMinus  = new TH2D("h_iso_pt_eleMinus",
      "e^{-} ABCD;relIso;p_{T}^{e} [GeV]", kNIsoAB, kIsoLoAB, kIsoHiAB, 50, 0, 100);
  // Same plane WITH the discriminant selection's m_T cut: supplies the
  // anti-iso lepton-pT shape for the QCD template of the pT+m_T>40 stacks.
  TH2D *h_iso_pt_mt40_elePlus  = new TH2D("h_iso_pt_mt40_elePlus",
      Form("e^{+} ABCD, m_{T}>%.0f;relIso;p_{T}^{e} [GeV]", mtCutForPtDisc), kNIsoAB, kIsoLoAB, kIsoHiAB, 50, 0, 100);
  TH2D *h_iso_pt_mt40_eleMinus = new TH2D("h_iso_pt_mt40_eleMinus",
      Form("e^{-} ABCD, m_{T}>%.0f;relIso;p_{T}^{e} [GeV]", mtCutForPtDisc), kNIsoAB, kIsoLoAB, kIsoHiAB, 50, 0, 100);
  h_iso_pt_mt40_elePlus->Sumw2(); h_iso_pt_mt40_eleMinus->Sumw2();
  h_iso_met_elePlus->Sumw2();  h_iso_met_eleMinus->Sumw2();
  h_iso_mt_elePlus->Sumw2();   h_iso_mt_eleMinus->Sumw2();
  h_iso_pt_elePlus->Sumw2();   h_iso_pt_eleMinus->Sumw2();

  // -------- (lepton pT, m_T) joint 2Ds for the pT x mT cut scan --------
  // Electron twins of the skim_Wmu pT x mT scan inputs (see the comment there;
  // consumed by correction/ptmt_scan.C+(true)). Iso stays at the nominal cut;
  // the QCD model is the ANTI-ISO sideband (relIso in [0.20, 1.0), the
  // qcd_abcd.C electron window) with its actual MET/mT -- no MET conditioning.
  const double antiIsoLo = 0.20, antiIsoHi = 1.00; // = qcd_abcd.C isoFail window (ele)
  TH2D *h_pt_mt_elePlus  = new TH2D("h_pt_mt_elePlus",
      "e^{+}, iso-pass;p_{T}^{e} [GeV];m_{T} [GeV]", 100, 0, 100, 80, 0, 200);
  TH2D *h_pt_mt_eleMinus = new TH2D("h_pt_mt_eleMinus",
      "e^{-}, iso-pass;p_{T}^{e} [GeV];m_{T} [GeV]", 100, 0, 100, 80, 0, 200);
  TH2D *h_pt_mt_antiiso_elePlus  = new TH2D("h_pt_mt_antiiso_elePlus",
      "e^{+}, anti-iso [0.20,1.0);p_{T}^{e} [GeV];m_{T} [GeV]", 100, 0, 100, 80, 0, 200);
  TH2D *h_pt_mt_antiiso_eleMinus = new TH2D("h_pt_mt_antiiso_eleMinus",
      "e^{-}, anti-iso [0.20,1.0);p_{T}^{e} [GeV];m_{T} [GeV]", 100, 0, 100, 80, 0, 200);
  h_pt_mt_elePlus->Sumw2();         h_pt_mt_eleMinus->Sumw2();
  h_pt_mt_antiiso_elePlus->Sumw2(); h_pt_mt_antiiso_eleMinus->Sumw2();

  // -------- Leading-lepton pT: per-charge, per-y (lab + FB) --------
  // Electron twins of the skim_Wmu h_leppt_* histos: pT versions of the
  // rapidity-binned h_met_*/h_mt_*, for the lepton-pT data/MC stacks in
  // plotting/mtandmet.C (plots only -- NOT in the Combine input). 2 GeV bins.
  TH1D *h_leppt_Wp[kNY],    *h_leppt_Wm[kNY];
  TH1D *h_leppt_Wp_FB[kNY], *h_leppt_Wm_FB[kNY];
  TH1D *h_leppt_mt40_Wp[kNY],    *h_leppt_mt40_Wm[kNY];
  TH1D *h_leppt_mt40_Wp_FB[kNY], *h_leppt_mt40_Wm_FB[kNY];
  for (int b = 0; b < kNY; ++b)
  {
    const double y1 = kYEdges[b],     y2 = kYEdges[b + 1];
    const double y1FB = kYEdgesFB[b], y2FB = kYEdgesFB[b + 1];

    h_leppt_Wp[b] = new TH1D(Form("h_leppt_Wp_y%d", b),
        Form("W+ lepton p_{T};p_{T}^{e} [GeV];Events (%.2f<y<%.2f)", y1, y2), 50, 0, 100);
    h_leppt_Wm[b] = new TH1D(Form("h_leppt_Wm_y%d", b),
        Form("W- lepton p_{T};p_{T}^{e} [GeV];Events (%.2f<y<%.2f)", y1, y2), 50, 0, 100);
    h_leppt_Wp_FB[b] = new TH1D(Form("h_leppt_Wp_y%d_FB", b),
        Form("W+ lepton p_{T};p_{T}^{e} [GeV];Events (%.2f<y<%.2f)", y1FB, y2FB), 50, 0, 100);
    h_leppt_Wm_FB[b] = new TH1D(Form("h_leppt_Wm_y%d_FB", b),
        Form("W- lepton p_{T};p_{T}^{e} [GeV];Events (%.2f<y<%.2f)", y1FB, y2FB), 50, 0, 100);

    h_leppt_mt40_Wp[b] = new TH1D(Form("h_leppt_mt40_Wp_y%d", b),
        Form("W+ lepton p_{T}, m_{T}>%.0f;p_{T}^{e} [GeV];Events (%.2f<y<%.2f)", mtCutForPtDisc, y1, y2), 50, 0, 100);
    h_leppt_mt40_Wm[b] = new TH1D(Form("h_leppt_mt40_Wm_y%d", b),
        Form("W- lepton p_{T}, m_{T}>%.0f;p_{T}^{e} [GeV];Events (%.2f<y<%.2f)", mtCutForPtDisc, y1, y2), 50, 0, 100);
    h_leppt_mt40_Wp_FB[b] = new TH1D(Form("h_leppt_mt40_Wp_y%d_FB", b),
        Form("W+ lepton p_{T}, m_{T}>%.0f;p_{T}^{e} [GeV];Events (%.2f<y<%.2f)", mtCutForPtDisc, y1FB, y2FB), 50, 0, 100);
    h_leppt_mt40_Wm_FB[b] = new TH1D(Form("h_leppt_mt40_Wm_y%d_FB", b),
        Form("W- lepton p_{T}, m_{T}>%.0f;p_{T}^{e} [GeV];Events (%.2f<y<%.2f)", mtCutForPtDisc, y1FB, y2FB), 50, 0, 100);

    h_leppt_Wp[b]->Sumw2();         h_leppt_Wm[b]->Sumw2();
    h_leppt_Wp_FB[b]->Sumw2();      h_leppt_Wm_FB[b]->Sumw2();
    h_leppt_mt40_Wp[b]->Sumw2();    h_leppt_mt40_Wm[b]->Sumw2();
    h_leppt_mt40_Wp_FB[b]->Sumw2(); h_leppt_mt40_Wm_FB[b]->Sumw2();
  }

  // -------- Leading-lepton kinematics (full W selection) --------
  // Electron twin of the skim_Wmu h_lep* histos (correction/dataMC_kinematics.C):
  // LEADING electron after the full selection (isolation + trigger match
  // included), charges combined. Same names as the Z-channel histos on purpose.
  TH1D *h_lepPt  = new TH1D("h_lepPt",  Form("%s; p_{T}^{e} [GeV]; Events", outPrefix.c_str()), 50,  0,  100);
  TH1D *h_lepEta = new TH1D("h_lepEta", Form("%s; #eta_{e}; Events",        outPrefix.c_str()), 50, -2.5, 2.5);
  TH1D *h_lepPhi = new TH1D("h_lepPhi", Form("%s; #phi_{e}; Events",        outPrefix.c_str()), 32, -TMath::Pi(), TMath::Pi());
  TH1D *h_lepPt_mt40  = new TH1D("h_lepPt_mt40",  Form("%s, m_{T}>%.0f; p_{T}^{e} [GeV]; Events", outPrefix.c_str(), mtCutForPtDisc), 50,  0,  100);
  TH1D *h_lepEta_mt40 = new TH1D("h_lepEta_mt40", Form("%s, m_{T}>%.0f; #eta_{e}; Events", outPrefix.c_str(), mtCutForPtDisc), 50, -2.5, 2.5);
  TH1D *h_lepPhi_mt40 = new TH1D("h_lepPhi_mt40", Form("%s, m_{T}>%.0f; #phi_{e}; Events", outPrefix.c_str(), mtCutForPtDisc), 32, -TMath::Pi(), TMath::Pi());
  h_lepPt->Sumw2(); h_lepEta->Sumw2(); h_lepPhi->Sumw2();
  h_lepPt_mt40->Sumw2(); h_lepEta_mt40->Sumw2(); h_lepPhi_mt40->Sumw2();

  // -------- Cutflow loop --------
  unsigned long long N[9] = {0};
  const Long64_t nEntries = tEle->GetEntries();
  std::cout << "Entries: " << nEntries << "\n";

  bool warnedEventFiltersOnce = false;
  bool warnedNoTrigMatchInfo  = false;

  for (Long64_t ie = 0; ie < nEntries; ++ie)
  {
    if (ie % 200000 == 0) std::cout << "Event " << ie << "/" << nEntries << "\n";

    tEle->GetEntry(ie);
    tHi->GetEntry(ie);
    tPF->GetEntry(ie);
    tHLT->GetEntry(ie);
    tHLTobj->GetEntry(ie);
    tEvent->GetEntry(ie);
    if (isMC) tGen->GetEntry(ie);

    const double w = has_genWeight ? (double)genWeight : 1.0;

    N[0]++;

    if (!elePt || !eleEta || !elePhi || !eleCharge) continue;

    // Precut relaxed to the (pT, mT)-scan floor (pT > 20) so the scan planes
    // h_pt_mt[_antiiso]_* reach below the nominal cut. The cutflow keeps its
    // nominal meaning: N[1..5] count only events the pT > 25 precut passes.
    const bool hasPF25 = ExistsPFElectronPt25_W(elePt25, nEle, elePt);
    if (!ExistsPFElectronPt25_W(kPtScanFloor, nEle, elePt)) continue;
    if (hasPF25) N[1]++;

    if (!PassEventSelection_pO(warnedEventFiltersOnce,
                               has_pprimaryVertexFilter,        pprimaryVertexFilter,
                               has_pclusterCompatibilityFilter, pclusterCompatibilityFilter))
      continue;
    if (applyVz && TMath::Abs(vz) > vzMax) continue;
    if (hasPF25) N[2]++;

    if (has_hlt && !TriggerFired(HLT_OxyL1SingleEG10_v1)) continue;
    if (hasPF25) N[3]++;

    // Electron ID: eleMVAIdWP95 (switched from eleCutIdWP95 on 2026-06-24 after the
    // isolation/ID study in correction/isolation_ele.C — the MVA WP is the best QCD
    // rejector while CutWP95 was the worst; continuous AUC_MET 0.911 vs 0.870). Used
    // consistently for the DY veto, the tight-ID gate, and the leading-electron pick.
    // ALL isolation gates are the continuous PF relIso < isoMax (0.095, ~optimal):
    // the leading electron AND (since 2026-07-02) the DY-veto legs, which
    // previously used the integer eleMVAIsoWP95 WP.
    if (!PassDYVeto_Wel(dyElePtMin, isoMax, dyMassMin, dyMassMax,
                        nEle, elePt, eleEta, elePhi, eleCharge,
                        eleMVAIdWP95, elePFChIso, elePFNeuIso, elePFPhoIso, elePFPUIso))
      continue;
    if (hasPF25) N[4]++;

    if (!ExistsTightElectron_W(nEle, eleMVAIdWP95)) continue;
    if (hasPF25) N[5]++;

    const int iLead = FindLeadingElectron_TightPF(nEle, elePt, eleEta, elePhi, eleMVAIdWP95);
    if (iLead < 0) continue;

    if (isMC)
    {
      bool passEleW = PassGenRecoMatchingWithAncestor(
          iLead,
          elePt, eleEta, elePhi, eleCharge,
          genPt, genEta, genPhi, genChg,
          genPdg, genMotherIdx,
          11, 24);
      if (!passEleW)
      {
        // cout << "We see a gen-reco matching failed case" << endl;
        // continue;
      }
    }

    if (elePt->at(iLead) <= kPtScanFloor) continue; // scan floor; nominal pT cut via passPtNominal
    if (std::abs(eleEta->at(iLead)) > eleEtaMax) continue;
    const bool passPtNominal = (elePt->at(iLead) > elePt25);
    if (passPtNominal) N[6]++;

    // (7) Leading electron Iso (continuous PF rel-iso; legacy WP path kept commented):
    // if (eleMVAIsoWP95->at(iLead) != 1) continue;
    const double isoLead = RelIsoPF(iLead, elePt, elePFChIso, elePFNeuIso, elePFPhoIso, elePFPUIso);
    const int    isoBin  = FindIsoBin(isoLead);
    const bool   isQCDSideband   = (isoBin >= 0);
    const bool   passIsoNominal  = (isoLead < isoMax);
    // Isolation study DONE (2026-06-24, correction/isolation_ele.C continuous scan):
    // continuous relIso < 0.095 is the Youden-J(QCD) optimum for the MVA IDs, so the
    // cut value stands. (2026-07-06: relIso now carries the Delta-beta PU correction,
    // consistent with the MuonPOG definition used in the muon channel.)

    if (passIsoNominal && passPtNominal) N[7]++;

    const bool passMatch = PassLeadingLeptonTrigMatch(
        trigMatchDR, iLead, eleEta, elePhi,
        has_trgObjPt, trgObjPt, has_trgObjEta, trgObjEta,
        has_trgObjPhi, trgObjPhi,
        warnedNoTrigMatchInfo, "electron");
    if (!passMatch) continue;
    if (passIsoNominal && passPtNominal) N[8]++;

    // -------- Fill final distributions (after step 8) --------
    TVector2 metv = ComputePFMET(pfId, pfPt, pfPhi);
    const double met = metv.Mod();
    const double mt  = TransverseMass(elePt->at(iLead), elePhi->at(iLead), metv);

    if (isQCDSideband && passPtNominal)
    {
      if      (eleCharge->at(iLead) > 0) h_met_iso_elePlus [isoBin]->Fill(met, w);
      else if (eleCharge->at(iLead) < 0) h_met_iso_eleMinus[isoBin]->Fill(met, w);
    }

    // ABCD: fill the full-iso-range (relIso, MET) and (relIso, MT) planes.
    // Unconditional on the isolation cut -- iso pass AND fail both populated.
    // The h_iso_* ABCD planes keep the NOMINAL pT > 25 selection; only the
    // (pT, mT) scan planes h_pt_mt[_antiiso]_* go down to kPtScanFloor.
    const bool passMtCut = (mt > mtCutForPtDisc); // lepton-pT-discriminant m_T cut
    if (eleCharge->at(iLead) > 0)
    {
      if (passPtNominal)
      {
        h_iso_met_elePlus->Fill(isoLead, met, w);
        h_iso_mt_elePlus ->Fill(isoLead, mt,  w);
        h_iso_pt_elePlus ->Fill(isoLead, elePt->at(iLead), w);
        if (passMtCut) h_iso_pt_mt40_elePlus->Fill(isoLead, elePt->at(iLead), w);
      }
      if (isoLead >= antiIsoLo && isoLead < antiIsoHi)
        h_pt_mt_antiiso_elePlus->Fill(elePt->at(iLead), mt, w);
    }
    else if (eleCharge->at(iLead) < 0)
    {
      if (passPtNominal)
      {
        h_iso_met_eleMinus->Fill(isoLead, met, w);
        h_iso_mt_eleMinus ->Fill(isoLead, mt,  w);
        h_iso_pt_eleMinus ->Fill(isoLead, elePt->at(iLead), w);
        if (passMtCut) h_iso_pt_mt40_eleMinus->Fill(isoLead, elePt->at(iLead), w);
      }
      if (isoLead >= antiIsoLo && isoLead < antiIsoHi)
        h_pt_mt_antiiso_eleMinus->Fill(elePt->at(iLead), mt, w);
    }

    if (!passIsoNominal) continue;

    // pT x mT scan signal plane (iso-pass; filled down to the pT > 20 scan floor)
    if      (eleCharge->at(iLead) > 0) h_pt_mt_elePlus ->Fill(elePt->at(iLead), mt, w);
    else if (eleCharge->at(iLead) < 0) h_pt_mt_eleMinus->Fill(elePt->at(iLead), mt, w);

    if (!passPtNominal) continue; // everything below keeps the nominal pT > 25 cut

    h_lepPt ->Fill(elePt->at(iLead),  w);
    h_lepEta->Fill(eleEta->at(iLead), w);
    h_lepPhi->Fill(elePhi->at(iLead), w);
    if (passMtCut)
    {
      h_lepPt_mt40 ->Fill(elePt->at(iLead),  w);
      h_lepEta_mt40->Fill(eleEta->at(iLead), w);
      h_lepPhi_mt40->Fill(elePhi->at(iLead), w);
    }

    const int q = eleCharge->at(iLead);
    const bool isWp = (q > 0), isWm = (q < 0);

    // [NOTICE] Sign flip: p-going (-Z) defined as forward
    const double y       = -eleEta->at(iLead);
    const int    ybin    = FindYBin(y, kYEdges,   kNY);
    const int    ybin_FB = FindYBin(y, kYEdgesFB, kNY);

    if (ybin >= 0)
    {
      if      (isWp) { h_met_Wp[ybin]->Fill(met, w); h_mt_Wp[ybin]->Fill(mt, w); h_leppt_Wp[ybin]->Fill(elePt->at(iLead), w); }
      else if (isWm) { h_met_Wm[ybin]->Fill(met, w); h_mt_Wm[ybin]->Fill(mt, w); h_leppt_Wm[ybin]->Fill(elePt->at(iLead), w); }
      if (passMtCut)
      {
        if      (isWp) h_leppt_mt40_Wp[ybin]->Fill(elePt->at(iLead), w);
        else if (isWm) h_leppt_mt40_Wm[ybin]->Fill(elePt->at(iLead), w);
      }
    }
    if (ybin_FB >= 0)
    {
      if      (isWp) { h_met_Wp_FB[ybin_FB]->Fill(met, w); h_mt_Wp_FB[ybin_FB]->Fill(mt, w); h_leppt_Wp_FB[ybin_FB]->Fill(elePt->at(iLead), w); }
      else if (isWm) { h_met_Wm_FB[ybin_FB]->Fill(met, w); h_mt_Wm_FB[ybin_FB]->Fill(mt, w); h_leppt_Wm_FB[ybin_FB]->Fill(elePt->at(iLead), w); }
      if (passMtCut)
      {
        if      (isWp) h_leppt_mt40_Wp_FB[ybin_FB]->Fill(elePt->at(iLead), w);
        else if (isWm) h_leppt_mt40_Wm_FB[ybin_FB]->Fill(elePt->at(iLead), w);
      }
    }
  }

  // -------- Cutflow + outputs --------
  PrintCutflow(std::cout, N);
  WriteCutflowTxt(outPrefix, isMC, N);

  gSystem->mkdir("rootfile", kTRUE); // ensure ./rootfile exists (fresh checkout)
  TFile *fout = new TFile(("./rootfile/" + outPrefix + "_hist.root").c_str(), "RECREATE");
  for (int i = 0; i < kNY; ++i)
  {
    h_met_Wp[i]->Write("", 2);
    h_met_Wm[i]->Write("", 2);
    h_mt_Wp [i]->Write("", 2);
    h_mt_Wm [i]->Write("", 2);
    h_met_Wp_FB[i]->Write("", 2);
    h_met_Wm_FB[i]->Write("", 2);
    h_mt_Wp_FB [i]->Write("", 2);
    h_mt_Wm_FB [i]->Write("", 2);
    h_leppt_Wp[i]->Write("", 2);
    h_leppt_Wm[i]->Write("", 2);
    h_leppt_Wp_FB[i]->Write("", 2);
    h_leppt_Wm_FB[i]->Write("", 2);
    h_leppt_mt40_Wp[i]->Write("", 2);
    h_leppt_mt40_Wm[i]->Write("", 2);
    h_leppt_mt40_Wp_FB[i]->Write("", 2);
    h_leppt_mt40_Wm_FB[i]->Write("", 2);
  }
  for (int b = 0; b < NISO; ++b)
  {
    h_met_iso_elePlus [b]->Write("", TObject::kOverwrite);
    h_met_iso_eleMinus[b]->Write("", TObject::kOverwrite);
  }
  // ABCD 2D planes (relIso vs MET / m_T / lepton pT), per charge
  h_iso_met_elePlus ->Write("", TObject::kOverwrite);
  h_iso_met_eleMinus->Write("", TObject::kOverwrite);
  h_iso_mt_elePlus  ->Write("", TObject::kOverwrite);
  h_iso_mt_eleMinus ->Write("", TObject::kOverwrite);
  h_iso_pt_elePlus  ->Write("", TObject::kOverwrite);
  h_iso_pt_eleMinus ->Write("", TObject::kOverwrite);
  h_iso_pt_mt40_elePlus ->Write("", TObject::kOverwrite);
  h_iso_pt_mt40_eleMinus->Write("", TObject::kOverwrite);
  // pT x mT scan inputs (correction/ptmt_scan.C)
  h_pt_mt_elePlus         ->Write("", TObject::kOverwrite);
  h_pt_mt_eleMinus        ->Write("", TObject::kOverwrite);
  h_pt_mt_antiiso_elePlus ->Write("", TObject::kOverwrite);
  h_pt_mt_antiiso_eleMinus->Write("", TObject::kOverwrite);
  // Leading-lepton kinematics (full selection; _mt40 = + the m_T>40 cut)
  h_lepPt ->Write("", TObject::kOverwrite);
  h_lepEta->Write("", TObject::kOverwrite);
  h_lepPhi->Write("", TObject::kOverwrite);
  h_lepPt_mt40 ->Write("", TObject::kOverwrite);
  h_lepEta_mt40->Write("", TObject::kOverwrite);
  h_lepPhi_mt40->Write("", TObject::kOverwrite);
  fout->Close();

  std::cout << "[INFO] Wrote outputs with prefix: " << outPrefix << "\n";
  return 0;
}


// ============================================================
// Z -> mu+mu
// Legacy: skim/legacy/DrawDimuonPeak.C
// ============================================================

int skim_Zmm(const char *fname, SampleType sample)
{
  // -------- Z dimuon config (kept identical to legacy DileptonConfig) --------
  double massMin = 60.0,  massMax = 120.0;
  int    nBins   = 120;
  double ptMin1  = 15.0; // leading
  double ptMin2  = 10.0; // subleading
  const double etaMax = 2.4;
  const bool   requireOS      = true;
  const bool   requireGood    = false;
  const bool   requireGlobal  = true;
  const bool   requirePF      = false;  // NOTE: original Z dimuon code sets this false
  const bool   requireTightID = true;
  const bool   applyVz        = true;
  const double vzMax          = 15.0;
  const bool   applyHiBin     = false;
  const int    hiBinMin       = 0;
  const int    hiBinMax       = 200;
  const double isoMax         = 0.2;

  std::string outPrefix = "ZToMuMu_pO2025";

  const bool isMC = IsMC(sample);

  if (isMC)
  {
    const auto info = ResolveMCSample(sample, "mu");
    if (!info.fname.empty()) fname = Form("%s", info.fname.c_str());
    outPrefix += info.outSuffix;
  }
  const std::string mcTag = isMC ? "MC" : "Data";

  std::cout << "[INPUT] " << fname << std::endl;
  TFile *f = TFile::Open(fname);
  if (!f || f->IsZombie()) { std::cerr << "ERROR: cannot open " << fname << "\n"; return 1; }

  TTree *tMu = (TTree *)f->Get("ggHiNtuplizer/EventTree");
  if (!tMu) { std::cerr << "ERROR: missing ggHiNtuplizer/EventTree\n"; return 1; }

  TTree *tHi = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  const bool haveHiTree = (tHi != nullptr);

  // -------- PF candidates (for the hadronic-recoil / MET correction) --------
  // Soft-gated: if pftree is absent the Z-mass plots still run, but the recoil
  // (u_par/u_perp) histograms stay empty -- warn loudly so it isn't silent.
  TTree *tPF = (TTree *)f->Get("particleFlowAnalyser/pftree");
  const bool haveRecoil = (tPF != nullptr);
  std::vector<int>   *pfId  = nullptr;
  std::vector<float> *pfPt  = nullptr;
  std::vector<float> *pfPhi = nullptr;
  if (haveRecoil)
  {
    tPF->SetBranchStatus("*", 0); // read only what ComputePFMET needs (speed)
    tPF->SetBranchStatus("pfPt", 1);
    tPF->SetBranchStatus("pfPhi", 1);
    tPF->SetBranchAddress("pfPt", &pfPt);
    tPF->SetBranchAddress("pfPhi", &pfPhi);
  }
  else
    std::cerr << "[WARN] skim_Zmm: missing particleFlowAnalyser/pftree in " << fname
              << " -- recoil (u_par/u_perp) histograms will be EMPTY.\n";

  // -------- Muon branches --------
  Int_t nMu = 0;
  std::vector<float> *muPt = nullptr, *muEta = nullptr, *muPhi = nullptr;
  std::vector<int>   *muCharge = nullptr;
  std::vector<int>   *muIsGood = nullptr, *muIsGlobal = nullptr;
  std::vector<int>   *muIsPF   = nullptr, *muIDTight  = nullptr;

  tMu->SetBranchStatus("*", 0); // speed: read only the muon branches we use
  tMu->SetBranchStatus("nMu", 1);
  tMu->SetBranchStatus("muPt", 1);
  tMu->SetBranchStatus("muEta", 1);
  tMu->SetBranchStatus("muPhi", 1);
  tMu->SetBranchStatus("muCharge", 1);
  tMu->SetBranchAddress("nMu",      &nMu);
  tMu->SetBranchAddress("muPt",     &muPt);
  tMu->SetBranchAddress("muEta",    &muEta);
  tMu->SetBranchAddress("muPhi",    &muPhi);
  tMu->SetBranchAddress("muCharge", &muCharge);

  const bool has_muIsGood   = HasBranch(tMu, "muIsGood");
  const bool has_muIsGlobal = HasBranch(tMu, "muIsGlobal");
  const bool has_muIsPF     = HasBranch(tMu, "muIsPF");
  const bool has_muIDTight  = HasBranch(tMu, "muIDTight");

  if (has_muIsGood)   { tMu->SetBranchStatus("muIsGood",   1); tMu->SetBranchAddress("muIsGood",   &muIsGood);   }
  if (has_muIsGlobal) { tMu->SetBranchStatus("muIsGlobal", 1); tMu->SetBranchAddress("muIsGlobal", &muIsGlobal); }
  if (has_muIsPF)     { tMu->SetBranchStatus("muIsPF",     1); tMu->SetBranchAddress("muIsPF",     &muIsPF);     }
  if (has_muIDTight)  { tMu->SetBranchStatus("muIDTight",  1); tMu->SetBranchAddress("muIDTight",  &muIDTight);  }

  // -------- Event branches --------
  Float_t vz    = 999.f;
  Int_t   hiBin = -999;

  const bool has_vz    = (haveHiTree && HasBranch(tHi, "vz"));
  const bool has_hiBin = (haveHiTree && HasBranch(tHi, "hiBin"));

  if (haveHiTree) tHi->SetBranchStatus("*", 0); // speed (HiTree is optional here)
  if (has_vz)    { tHi->SetBranchStatus("vz", 1);    tHi->SetBranchAddress("vz",    &vz); }
  if (has_hiBin) { tHi->SetBranchStatus("hiBin", 1); tHi->SetBranchAddress("hiBin", &hiBin); }

  // -------- Generator event weight (MC only) --------
  Float_t    genWeight     = 1.f;
  const bool has_genWeight = isMC && haveHiTree && HasBranch(tHi, "weight");
  if (has_genWeight) { tHi->SetBranchStatus("weight", 1); tHi->SetBranchAddress("weight", &genWeight); }
  else if (isMC)     std::cout << "[WARN] skim_Zmm: MC sample but no 'weight' branch on HiTree; filling unweighted.\n";

  // -------- PF iso --------
  std::vector<float> *muPFChIso = nullptr, *muPFNeuIso = nullptr, *muPFPhoIso = nullptr;
  std::vector<float> *muPFPUIso = nullptr;

  const bool has_muPFChIso  = HasBranch(tMu, "muPFChIso");
  const bool has_muPFNeuIso = HasBranch(tMu, "muPFNeuIso");
  const bool has_muPFPhoIso = HasBranch(tMu, "muPFPhoIso");
  const bool has_muPFPUIso  = HasBranch(tMu, "muPFPUIso");

  if (has_muPFChIso)  { tMu->SetBranchStatus("muPFChIso",  1); tMu->SetBranchAddress("muPFChIso",  &muPFChIso);  }
  if (has_muPFNeuIso) { tMu->SetBranchStatus("muPFNeuIso", 1); tMu->SetBranchAddress("muPFNeuIso", &muPFNeuIso); }
  if (has_muPFPhoIso) { tMu->SetBranchStatus("muPFPhoIso", 1); tMu->SetBranchAddress("muPFPhoIso", &muPFPhoIso); }
  if (has_muPFPUIso)  { tMu->SetBranchStatus("muPFPUIso",  1); tMu->SetBranchAddress("muPFPUIso",  &muPFPUIso);  }
  if (!has_muPFPUIso)
    std::cout << "[WARN] skim_Zmm: no muPFPUIso branch; relIso falls back to uncorrected (no Delta-beta).\n";

  // -------- HLT objects --------
  TTree *tHLTobj = (TTree *)f->Get("hltobject/HLT_OxyL1SingleMuOpen_v");
  if (!tHLTobj)
  {
    std::cerr << "[FATAL] skim_Zmm: missing TTree hltobject/HLT_OxyL1SingleMuOpen_v in " << fname << "\n";
    f->Close();
    return 2;
  }

  std::vector<double> *trgObjPt = nullptr, *trgObjEta = nullptr;
  std::vector<double> *trgObjPhi = nullptr, *trgObjId = nullptr;

  const bool has_trgObjPt  = HasBranch(tHLTobj, "pt");
  const bool has_trgObjEta = HasBranch(tHLTobj, "eta");
  const bool has_trgObjPhi = HasBranch(tHLTobj, "phi");
  const bool has_trgObjId  = HasBranch(tHLTobj, "TriggerObjID");

  tHLTobj->SetBranchStatus("*", false);
  if (has_trgObjPt && has_trgObjEta && has_trgObjPhi && has_trgObjId)
  {
    tHLTobj->SetBranchStatus("pt", 1);
    tHLTobj->SetBranchStatus("eta", 1);
    tHLTobj->SetBranchStatus("phi", 1);
    tHLTobj->SetBranchStatus("TriggerObjID", 1);
    tHLTobj->SetBranchAddress("pt",           &trgObjPt);
    tHLTobj->SetBranchAddress("eta",          &trgObjEta);
    tHLTobj->SetBranchAddress("phi",          &trgObjPhi);
    tHLTobj->SetBranchAddress("TriggerObjID", &trgObjId);
  }
  else
  {
    std::cout << "[WARN] Could not find HLT related Object." << std::endl;
  }

  // -------- HLT bit --------
  TTree *tHLT = (TTree *)f->Get("hltanalysis/HltTree");
  if (!tHLT)
  {
    std::cerr << "[FATAL] skim_Zmm: missing TTree hltanalysis/HltTree in " << fname << "\n";
    f->Close();
    return 2;
  }

  const std::string hltName = FindBranchContaining(tHLT, "HLT_OxyL1SingleMuOpen_v1");
  Int_t HLT_OxyL1SingleMuOpen_v1 = 0;
  bool  has_hlt = false;
  tHLT->SetBranchStatus("*", 0);
  if (!hltName.empty())
  {
    has_hlt = true;
    tHLT->SetBranchStatus(hltName.c_str(), 1);
    tHLT->SetBranchAddress(hltName.c_str(), &HLT_OxyL1SingleMuOpen_v1);
    std::cout << "[INFO] Using trigger bit branch: " << hltName << "\n";
  }
  else
  {
    std::cout << "[WARN] Could not find any branch containing 'HLT_OxyL1SingleMuOpen_v1' in hltanalysis/HltTree.\n"
              << "       Step 3 (trigger fired) will be treated as PASS.\n";
  }

  // -------- Event filters --------
  TTree *tEvent = (TTree *)f->Get("skimanalysis/HltTree");
  if (!tEvent)
  {
    std::cerr << "[FATAL] skim_Zmm: missing TTree skimanalysis/HltTree in " << fname << "\n";
    f->Close();
    return 2;
  }

  Int_t pprimaryVertexFilter        = 1;
  Int_t pclusterCompatibilityFilter = 1;

  const bool has_pprimaryVertexFilter        = HasBranch(tEvent, "pprimaryVertexFilter");
  const bool has_pclusterCompatibilityFilter = HasBranch(tEvent, "pclusterCompatibilityFilter");

  tEvent->SetBranchStatus("*", 0);
  if (has_pprimaryVertexFilter && has_pclusterCompatibilityFilter)
  {
    tEvent->SetBranchStatus("pprimaryVertexFilter", 1);
    tEvent->SetBranchStatus("pclusterCompatibilityFilter", 1);
    tEvent->SetBranchAddress("pprimaryVertexFilter",        &pprimaryVertexFilter);
    tEvent->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);
  }

  // -------- Histograms --------
  gStyle->SetOptStat(0);
  TH1D *hMass = new TH1D("hMass",
      Form("%s; m_{#mu#mu} [GeV]; Events", outPrefix.c_str()),
      nBins, massMin, massMax);
  TH1D *hMass_extended = new TH1D("hMass_extended",
      Form("%s; m_{#mu#mu} [GeV]; Events", outPrefix.c_str()),
      nBins, 15, 120);
  TH1D *hMass_vipul = new TH1D("hMass_vipul",
      Form("%s; m_{#mu#mu} [GeV]; Events", outPrefix.c_str()),
      nBins, massMin, massMax);
  hMass->Sumw2(); hMass_extended->Sumw2(); hMass_vipul->Sumw2();

  // -------- Kinematics of the dimuon (Z) system and its muons --------
  // Filled for iso-selected OS pairs inside the Z peak [60,120] GeV (below).
  // Used for the Data/MC kinematic (boson-pT) check in correction/.
  TH1D *h_Zpt   = new TH1D("h_Zpt",  Form("%s; p_{T}^{#mu#mu} [GeV]; Events", outPrefix.c_str()), 50,  0,  100);
  TH1D *h_Zeta  = new TH1D("h_Zeta", Form("%s; #eta_{#mu#mu}; Events",        outPrefix.c_str()), 50, -5,    5);
  TH1D *h_Zy    = new TH1D("h_Zy",   Form("%s; y_{#mu#mu}; Events",           outPrefix.c_str()), 50, -5,    5);
  TH1D *h_Zphi  = new TH1D("h_Zphi", Form("%s; #phi_{#mu#mu}; Events",        outPrefix.c_str()), 32, -TMath::Pi(), TMath::Pi());
  TH1D *h_lepPt  = new TH1D("h_lepPt",  Form("%s; p_{T}^{#mu} [GeV]; Muons", outPrefix.c_str()), 50,  0,  100);
  TH1D *h_lepEta = new TH1D("h_lepEta", Form("%s; #eta_{#mu}; Muons",        outPrefix.c_str()), 50, -2.5, 2.5);
  TH1D *h_lepPhi = new TH1D("h_lepPhi", Form("%s; #phi_{#mu}; Muons",        outPrefix.c_str()), 32, -TMath::Pi(), TMath::Pi());
  h_Zpt->Sumw2();   h_Zeta->Sumw2();   h_Zphi->Sumw2();   h_Zy->Sumw2();
  h_lepPt->Sumw2(); h_lepEta->Sumw2(); h_lepPhi->Sumw2();

  // -------- Hadronic recoil (for the MET recoil correction) --------
  // u_par / u_perp = recoil components parallel / perpendicular to q_T (the
  // dimuon pT), from u = -MET - q_T. Stored inclusively (1D) and vs q_T (2D)
  // so the q_T binning for the recoil fits can be chosen later by projecting
  // slices -- no need to re-run the skim. u-axis uses the AN's 2 GeV/c binning.
  // Filled for the same selection as the kinematics (iso OS pairs, Z peak).
  TH1D *h_uPar  = new TH1D("h_uPar",  Form("%s; u_{#parallel} [GeV]; Events", outPrefix.c_str()), 200, -200, 200);
  TH1D *h_uPerp = new TH1D("h_uPerp", Form("%s; u_{#perp} [GeV]; Events",     outPrefix.c_str()), 200, -200, 200);
  TH2D *h_uPar_qT  = new TH2D("h_uPar_qT",  Form("%s; q_{T} [GeV]; u_{#parallel} [GeV]", outPrefix.c_str()), 70, 0, 140, 200, -200, 200);
  TH2D *h_uPerp_qT = new TH2D("h_uPerp_qT", Form("%s; q_{T} [GeV]; u_{#perp} [GeV]",     outPrefix.c_str()), 70, 0, 140, 200, -200, 200);
  h_uPar->Sumw2(); h_uPerp->Sumw2(); h_uPar_qT->Sumw2(); h_uPerp_qT->Sumw2();

  // -------- Loop --------
  const Long64_t nEntries = tMu->GetEntries();
  std::cout << "EventTree entries = " << nEntries << "\n";
  std::cout << "Have HiTree? " << (haveHiTree ? "YES" : "NO") << "\n";
  if (haveHiTree)
  {
    std::cout << "HiTree has vz? "    << (has_vz    ? "YES" : "NO") << "\n";
    std::cout << "HiTree has hiBin? " << (has_hiBin ? "YES" : "NO") << "\n";
  }

  Long64_t nPassEvent = 0;
  Long64_t nPassPair  = 0;
  bool warnedEventFiltersOnce = false;

  for (Long64_t ie = 0; ie < nEntries; ++ie)
  {
    if (ie % 200000 == 0) std::cout << "Event " << ie << "/" << nEntries << "\n";

    tMu->GetEntry(ie);
    if (haveHiTree) tHi->GetEntry(ie);
    if (haveRecoil) tPF->GetEntry(ie);
    tHLT->GetEntry(ie);
    tHLTobj->GetEntry(ie);
    tEvent->GetEntry(ie);

    const double w = has_genWeight ? (double)genWeight : 1.0;

    if (applyVz && has_vz && TMath::Abs(vz) > vzMax) continue;

    if (!PassEventSelection_pO(warnedEventFiltersOnce,
                               has_pprimaryVertexFilter,        pprimaryVertexFilter,
                               has_pclusterCompatibilityFilter, pclusterCompatibilityFilter))
      continue;

    if (applyHiBin && has_hiBin && (hiBin < hiBinMin || hiBin > hiBinMax)) continue;

    if (nMu < 2) continue;

    if (has_hlt && !TriggerFired(HLT_OxyL1SingleMuOpen_v1)) continue;
    nPassEvent++;

    for (int i = 0; i < nMu; ++i)
    {
      bool passIsolead = true;

      if (!muPt || !muEta || !muPhi || !muCharge) continue;
      if (muPt->at(i) < ptMin2) continue;
      if (TMath::Abs(muEta->at(i)) > etaMax) continue;

      if (requireGood    && has_muIsGood   && muIsGood->at(i)   == 0) continue;
      if (requireGlobal  && has_muIsGlobal && muIsGlobal->at(i) == 0) continue;
      if (requirePF      && has_muIsPF     && muIsPF->at(i)     == 0) continue;
      if (requireTightID && has_muIDTight  && muIDTight->at(i)  == 0) continue;

      const double isoLead = RelIsoPF(i, muPt, muPFChIso, muPFNeuIso, muPFPhoIso, muPFPUIso);
      if (isoLead >= isoMax) passIsolead = false;

      for (int j = i + 1; j < nMu; ++j)
      {
        bool passIsosec = true;

        if (muPt->at(j) < ptMin2) continue;
        if (TMath::Abs(muEta->at(j)) > etaMax) continue;

        if (requireGood    && has_muIsGood   && muIsGood->at(j)   == 0) continue;
        if (requireGlobal  && has_muIsGlobal && muIsGlobal->at(j) == 0) continue;
        if (requirePF      && has_muIsPF     && muIsPF->at(j)     == 0) continue;
        if (requireTightID && has_muIDTight  && muIDTight->at(j)  == 0) continue;

        if (requireOS && muCharge->at(i) * muCharge->at(j) >= 0) continue;

        const double pt1  = muPt->at(i);
        const double pt2  = muPt->at(j);
        const double lead = (pt1 > pt2 ? pt1 : pt2);
        const double sub  = (pt1 > pt2 ? pt2 : pt1);
        if (lead < ptMin1) continue;
        if (sub  < ptMin2) continue;

        TLorentzVector v1, v2;
        v1.SetPtEtaPhiM(muPt->at(i), muEta->at(i), muPhi->at(i), MU_MASS);
        v2.SetPtEtaPhiM(muPt->at(j), muEta->at(j), muPhi->at(j), MU_MASS);
        const double m = (v1 + v2).M();

        const double isosub = RelIsoPF(j, muPt, muPFChIso, muPFNeuIso, muPFPhoIso, muPFPUIso);
        if (isosub >= isoMax) passIsosec = false;

        if (passIsolead && passIsosec)
        {
          hMass_extended->Fill(m, w);
          if (!(m < 60 || m > 120))
          {
            hMass->Fill(m, w);
            const TLorentzVector ll = v1 + v2;
            h_Zpt ->Fill(ll.Pt(),  w);
            h_Zeta->Fill(ll.Eta(), w);
            h_Zy  ->Fill(ll.Rapidity(), w);
            h_Zphi->Fill(ll.Phi(), w);
            h_lepPt ->Fill(v1.Pt(),  w); h_lepPt ->Fill(v2.Pt(),  w);
            h_lepEta->Fill(v1.Eta(), w); h_lepEta->Fill(v2.Eta(), w);
            h_lepPhi->Fill(v1.Phi(), w); h_lepPhi->Fill(v2.Phi(), w);

            // Hadronic recoil: u = -MET - q_T, with q_T = the dimuon (ll).
            // u_par is along q_T (should peak near -q_T), u_perp is transverse.
            if (haveRecoil)
            {
              const TVector2 metv = ComputePFMET(pfId, pfPt, pfPhi);
              const double qT = ll.Pt();
              if (qT > 0.0)
              {
                const double ux = -metv.X() - ll.Px();
                const double uy = -metv.Y() - ll.Py();
                const double cphi = ll.Px() / qT, sphi = ll.Py() / qT;
                const double uPar  =  ux * cphi + uy * sphi;
                const double uPerp = -ux * sphi + uy * cphi;
                h_uPar    ->Fill(uPar,  w);
                h_uPerp   ->Fill(uPerp, w);
                h_uPar_qT ->Fill(qT, uPar,  w);
                h_uPerp_qT->Fill(qT, uPerp, w);
              }
            }
          }
        }

        if (m < 60 || m > 120) continue;
        if (lead < 20 || sub < 20) continue;

        // Tighter pT cut, no iso requirement.
        hMass_vipul->Fill(m, w);
        nPassPair++;
      }
    }
  }

  std::cout << "Passed event preselection: " << nPassEvent << "\n";
  std::cout << "Found selected OS pair: "    << nPassPair  << "\n";

  gSystem->mkdir("rootfile", kTRUE); // ensure ./rootfile exists (fresh checkout)
  TFile *fout = new TFile(("./rootfile/" + outPrefix + "_" + mcTag + "_hist.root").c_str(), "RECREATE");
  hMass->Write("", 2);
  hMass_extended->Write("", 2);
  hMass_vipul->Write("", 2);
  h_Zpt->Write("", 2);   h_Zeta->Write("", 2);   h_Zphi->Write("", 2);   h_Zy->Write("", 2);
  h_lepPt->Write("", 2); h_lepEta->Write("", 2); h_lepPhi->Write("", 2);
  h_uPar->Write("", 2); h_uPerp->Write("", 2);
  h_uPar_qT->Write("", 2); h_uPerp_qT->Write("", 2);
  fout->Close();

  std::cout << "Wrote: " << outPrefix << "_mass.png/.pdf and _hist.root\n";
  return 0;
}


// ============================================================
// Z -> e+e
// Legacy: skim/legacy/DrawDielectronPeak.C
// ============================================================

int skim_Zee(const char *fname, SampleType sample)
{
  // -------- Z dielectron config (kept identical to legacy DileptonConfig) --------
  double massMin = 60.0,  massMax = 120.0;
  int    nBins   = 120;
  double ptMin1  = 15.0;
  double ptMin2  = 10.0;
  const double etaMax         = 2.4;   // FIXME: ECAL is 2.5 with 1.4442-1.566 gap
  const bool   requireOS      = true;
  const bool   requireTightID = true;  // active gate: eleMVAIdWP95 (switched from eleCutIdWP95 2026-07-02, matching skim_Wel + the iso/ID study)
  const bool   applyVz        = true;
  const double vzMax          = 15.0;
  const bool   applyHiBin     = false;
  const int    hiBinMin       = 0;
  const int    hiBinMax       = 200;
  const double isoMax         = 0.095; // FIXME: tune for electrons (PU correction missing)

  std::string outPrefix = "ZToEE_pO2025";

  const bool isMC = IsMC(sample);

  if (isMC)
  {
    const auto info = ResolveMCSample(sample, "ele");
    if (!info.fname.empty()) fname = Form("%s", info.fname.c_str());
    outPrefix += info.outSuffix;
  }
  const std::string mcTag = isMC ? "MC" : "Data";

  std::cout << "[INPUT] " << fname << std::endl;
  TFile *f = TFile::Open(fname);
  if (!f || f->IsZombie()) { std::cerr << "ERROR: cannot open " << fname << "\n"; return 1; }

  TTree *tEle = (TTree *)f->Get("ggHiNtuplizer/EventTree");
  if (!tEle) { std::cerr << "ERROR: missing ggHiNtuplizer/EventTree\n"; return 1; }

  TTree *tHi = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  const bool haveHiTree = (tHi != nullptr);

  // -------- PF candidates (for the hadronic-recoil / MET correction) --------
  // Soft-gated: if pftree is absent the Z-mass plots still run, but the recoil
  // (u_par/u_perp) histograms stay empty -- warn loudly so it isn't silent.
  TTree *tPF = (TTree *)f->Get("particleFlowAnalyser/pftree");
  const bool haveRecoil = (tPF != nullptr);
  std::vector<int>   *pfId  = nullptr;
  std::vector<float> *pfPt  = nullptr;
  std::vector<float> *pfPhi = nullptr;
  if (haveRecoil)
  {
    tPF->SetBranchStatus("*", 0); // read only what ComputePFMET needs (speed)
    tPF->SetBranchStatus("pfPt", 1);
    tPF->SetBranchStatus("pfPhi", 1);
    tPF->SetBranchAddress("pfPt", &pfPt);
    tPF->SetBranchAddress("pfPhi", &pfPhi);
  }
  else
    std::cerr << "[WARN] skim_Zee: missing particleFlowAnalyser/pftree in " << fname
              << " -- recoil (u_par/u_perp) histograms will be EMPTY.\n";

  // -------- Electron kinematics --------
  Int_t nEle = 0;
  std::vector<float> *elePt = nullptr, *eleEta = nullptr, *elePhi = nullptr;
  std::vector<int>   *eleCharge = nullptr;
  tEle->SetBranchStatus("*", 0); // speed: read only the electron branches we use
  tEle->SetBranchStatus("nEle", 1);
  tEle->SetBranchStatus("elePt", 1);
  tEle->SetBranchStatus("eleEta", 1);
  tEle->SetBranchStatus("elePhi", 1);
  tEle->SetBranchStatus("eleCharge", 1);
  tEle->SetBranchAddress("nEle",      &nEle);
  tEle->SetBranchAddress("elePt",     &elePt);
  tEle->SetBranchAddress("eleEta",    &eleEta);
  tEle->SetBranchAddress("elePhi",    &elePhi);
  tEle->SetBranchAddress("eleCharge", &eleCharge);

  // -------- Electron IDs (kept around; active gate uses eleMVAIdWP95) --------
  std::vector<int> *eleMVAIdWP80 = nullptr, *eleMVAIdWP85 = nullptr;
  std::vector<int> *eleMVAIdWP90 = nullptr, *eleMVAIdWP95 = nullptr;
  std::vector<int> *eleCutIdWP70 = nullptr, *eleCutIdWP80 = nullptr;
  std::vector<int> *eleCutIdWP90 = nullptr, *eleCutIdWP95 = nullptr;

  auto hookID = [&](const char *bname, std::vector<int> **addr)
  {
    if (HasBranch(tEle, bname))
    {
      tEle->SetBranchStatus(bname, 1);
      tEle->SetBranchAddress(bname, addr);
    }
  };
  hookID("eleMVAIdWP80", &eleMVAIdWP80);
  hookID("eleMVAIdWP85", &eleMVAIdWP85);
  hookID("eleMVAIdWP90", &eleMVAIdWP90);
  hookID("eleMVAIdWP95", &eleMVAIdWP95);
  hookID("eleCutIdWP70", &eleCutIdWP70);
  hookID("eleCutIdWP80", &eleCutIdWP80);
  hookID("eleCutIdWP90", &eleCutIdWP90);
  hookID("eleCutIdWP95", &eleCutIdWP95);

  const bool has_eleID = HasBranch(tEle, "eleMVAIdWP95");

  // -------- Iso WPs (kept around; active gate uses RelIsoPF) --------
  std::vector<int> *eleMVAIsoWP80 = nullptr, *eleMVAIsoWP85 = nullptr;
  std::vector<int> *eleMVAIsoWP90 = nullptr, *eleMVAIsoWP95 = nullptr;
  hookID("eleMVAIsoWP80", &eleMVAIsoWP80);
  hookID("eleMVAIsoWP85", &eleMVAIsoWP85);
  hookID("eleMVAIsoWP90", &eleMVAIsoWP90);
  hookID("eleMVAIsoWP95", &eleMVAIsoWP95);

  std::vector<float> *elePFChIso = nullptr, *elePFNeuIso = nullptr, *elePFPhoIso = nullptr;
  std::vector<float> *elePFPUIso = nullptr;
  const bool has_elePFChIso  = HasBranch(tEle, "elePFChIso");
  const bool has_elePFNeuIso = HasBranch(tEle, "elePFNeuIso");
  const bool has_elePFPhoIso = HasBranch(tEle, "elePFPhoIso");
  const bool has_elePFPUIso  = HasBranch(tEle, "elePFPUIso");
  if (has_elePFChIso)  { tEle->SetBranchStatus("elePFChIso", 1);  tEle->SetBranchAddress("elePFChIso",  &elePFChIso);  }
  if (has_elePFNeuIso) { tEle->SetBranchStatus("elePFNeuIso", 1); tEle->SetBranchAddress("elePFNeuIso", &elePFNeuIso); }
  if (has_elePFPhoIso) { tEle->SetBranchStatus("elePFPhoIso", 1); tEle->SetBranchAddress("elePFPhoIso", &elePFPhoIso); }
  if (has_elePFPUIso)  { tEle->SetBranchStatus("elePFPUIso", 1);  tEle->SetBranchAddress("elePFPUIso",  &elePFPUIso);  }
  if (!has_elePFPUIso)
    std::cout << "[WARN] skim_Zee: no elePFPUIso branch; relIso falls back to uncorrected (no Delta-beta).\n";

  // -------- Event branches --------
  Float_t vz    = 999.f;
  Int_t   hiBin = -999;
  const bool has_vz    = (haveHiTree && HasBranch(tHi, "vz"));
  const bool has_hiBin = (haveHiTree && HasBranch(tHi, "hiBin"));
  if (haveHiTree) tHi->SetBranchStatus("*", 0); // speed (HiTree is optional here)
  if (has_vz)    { tHi->SetBranchStatus("vz", 1);    tHi->SetBranchAddress("vz",    &vz); }
  if (has_hiBin) { tHi->SetBranchStatus("hiBin", 1); tHi->SetBranchAddress("hiBin", &hiBin); }

  // -------- Generator event weight (MC only) --------
  Float_t    genWeight     = 1.f;
  const bool has_genWeight = isMC && haveHiTree && HasBranch(tHi, "weight");
  if (has_genWeight) { tHi->SetBranchStatus("weight", 1); tHi->SetBranchAddress("weight", &genWeight); }
  else if (isMC)     std::cout << "[WARN] skim_Zee: MC sample but no 'weight' branch on HiTree; filling unweighted.\n";

  // -------- HLT objects --------
  TTree *tHLTobj = (TTree *)f->Get("hltobject/HLT_OxyL1SingleEG10_v");
  if (!tHLTobj)
  {
    std::cerr << "[FATAL] skim_Zee: missing TTree hltobject/HLT_OxyL1SingleEG10_v in " << fname << "\n";
    f->Close();
    return 2;
  }

  std::vector<double> *trgObjPt = nullptr, *trgObjEta = nullptr;
  std::vector<double> *trgObjPhi = nullptr, *trgObjId = nullptr;

  const bool has_trgObjPt  = HasBranch(tHLTobj, "pt");
  const bool has_trgObjEta = HasBranch(tHLTobj, "eta");
  const bool has_trgObjPhi = HasBranch(tHLTobj, "phi");
  const bool has_trgObjId  = HasBranch(tHLTobj, "TriggerObjID");

  tHLTobj->SetBranchStatus("*", false);
  if (has_trgObjPt && has_trgObjEta && has_trgObjPhi && has_trgObjId)
  {
    tHLTobj->SetBranchStatus("pt", 1);
    tHLTobj->SetBranchStatus("eta", 1);
    tHLTobj->SetBranchStatus("phi", 1);
    tHLTobj->SetBranchStatus("TriggerObjID", 1);
    tHLTobj->SetBranchAddress("pt",           &trgObjPt);
    tHLTobj->SetBranchAddress("eta",          &trgObjEta);
    tHLTobj->SetBranchAddress("phi",          &trgObjPhi);
    tHLTobj->SetBranchAddress("TriggerObjID", &trgObjId);
  }
  else
  {
    std::cout << "[WARN] Could not find HLT related Object." << std::endl;
  }

  // -------- HLT bit --------
  TTree *tHLT = (TTree *)f->Get("hltanalysis/HltTree");
  if (!tHLT)
  {
    std::cerr << "[FATAL] skim_Zee: missing TTree hltanalysis/HltTree in " << fname << "\n";
    f->Close();
    return 2;
  }

  const std::string hltName = FindBranchContaining(tHLT, "HLT_OxyL1SingleEG10_v1");
  Int_t HLT_OxyL1SingleEG10_v1 = 0;
  bool  has_hlt = false;
  tHLT->SetBranchStatus("*", 0);
  if (!hltName.empty())
  {
    has_hlt = true;
    tHLT->SetBranchStatus(hltName.c_str(), 1);
    tHLT->SetBranchAddress(hltName.c_str(), &HLT_OxyL1SingleEG10_v1);
    std::cout << "[INFO] Using trigger bit branch: " << hltName << "\n";
  }
  else
  {
    std::cout << "[WARN] Could not find any branch containing 'HLT_OxyL1SingleEG10_v1' in hltanalysis/HltTree.\n"
              << "       Trigger gate will be treated as PASS.\n";
  }

  // -------- Event filters --------
  TTree *tEvent = (TTree *)f->Get("skimanalysis/HltTree");
  if (!tEvent)
  {
    std::cerr << "[FATAL] skim_Zee: missing TTree skimanalysis/HltTree in " << fname << "\n";
    f->Close();
    return 2;
  }

  Int_t pprimaryVertexFilter        = 1;
  Int_t pclusterCompatibilityFilter = 1;

  const bool has_pprimaryVertexFilter        = HasBranch(tEvent, "pprimaryVertexFilter");
  const bool has_pclusterCompatibilityFilter = HasBranch(tEvent, "pclusterCompatibilityFilter");

  tEvent->SetBranchStatus("*", 0);
  if (has_pprimaryVertexFilter && has_pclusterCompatibilityFilter)
  {
    tEvent->SetBranchStatus("pprimaryVertexFilter", 1);
    tEvent->SetBranchStatus("pclusterCompatibilityFilter", 1);
    tEvent->SetBranchAddress("pprimaryVertexFilter",        &pprimaryVertexFilter);
    tEvent->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);
  }

  // -------- Histograms --------
  gStyle->SetOptStat(0);
  TH1D *hMass = new TH1D("hMass",
      Form("%s; m_{ee} [GeV]; Events", outPrefix.c_str()),
      nBins, massMin, massMax);
  TH1D *hMass_extended = new TH1D("hMass_extended",
      Form("%s; m_{ee} [GeV]; Events", outPrefix.c_str()),
      nBins, 15, 120);
  TH1D *hMass_vipul = new TH1D("hMass_vipul",
      Form("%s; m_{ee} [GeV]; Events", outPrefix.c_str()),
      nBins, massMin, massMax);
  hMass->Sumw2(); hMass_extended->Sumw2(); hMass_vipul->Sumw2();

  // -------- Kinematics of the dielectron (Z) system and its electrons --------
  // Filled for iso-selected OS pairs inside the Z peak [60,120] GeV (below).
  // Used for the Data/MC kinematic (boson-pT) check in correction/.
  TH1D *h_Zpt   = new TH1D("h_Zpt",  Form("%s; p_{T}^{ee} [GeV]; Events", outPrefix.c_str()), 50,  0,  100);
  TH1D *h_Zeta  = new TH1D("h_Zeta", Form("%s; #eta_{ee}; Events",        outPrefix.c_str()), 50, -5,    5);
  TH1D *h_Zy    = new TH1D("h_Zy",   Form("%s; y_{ee}; Events",           outPrefix.c_str()), 50, -5,    5);
  TH1D *h_Zphi  = new TH1D("h_Zphi", Form("%s; #phi_{ee}; Events",        outPrefix.c_str()), 32, -TMath::Pi(), TMath::Pi());
  TH1D *h_lepPt  = new TH1D("h_lepPt",  Form("%s; p_{T}^{e} [GeV]; Electrons", outPrefix.c_str()), 50,  0,  100);
  TH1D *h_lepEta = new TH1D("h_lepEta", Form("%s; #eta_{e}; Electrons",        outPrefix.c_str()), 50, -2.5, 2.5);
  TH1D *h_lepPhi = new TH1D("h_lepPhi", Form("%s; #phi_{e}; Electrons",        outPrefix.c_str()), 32, -TMath::Pi(), TMath::Pi());
  h_Zpt->Sumw2();   h_Zeta->Sumw2();   h_Zphi->Sumw2();   h_Zy->Sumw2();
  h_lepPt->Sumw2(); h_lepEta->Sumw2(); h_lepPhi->Sumw2();

  // -------- Hadronic recoil (for the MET recoil correction) --------
  // u_par / u_perp = recoil components parallel / perpendicular to q_T (the
  // dielectron pT), from u = -MET - q_T. Stored inclusively (1D) and vs q_T (2D)
  // so the q_T binning for the recoil fits can be chosen later by projecting
  // slices -- no need to re-run the skim. u-axis uses the AN's 2 GeV/c binning.
  // Filled for the same selection as the kinematics (iso OS pairs, Z peak).
  // Mirrors skim_Zmm (the muon recoil); see correction/recoil_raw_ele.C.
  TH1D *h_uPar  = new TH1D("h_uPar",  Form("%s; u_{#parallel} [GeV]; Events", outPrefix.c_str()), 200, -200, 200);
  TH1D *h_uPerp = new TH1D("h_uPerp", Form("%s; u_{#perp} [GeV]; Events",     outPrefix.c_str()), 200, -200, 200);
  TH2D *h_uPar_qT  = new TH2D("h_uPar_qT",  Form("%s; q_{T} [GeV]; u_{#parallel} [GeV]", outPrefix.c_str()), 70, 0, 140, 200, -200, 200);
  TH2D *h_uPerp_qT = new TH2D("h_uPerp_qT", Form("%s; q_{T} [GeV]; u_{#perp} [GeV]",     outPrefix.c_str()), 70, 0, 140, 200, -200, 200);
  h_uPar->Sumw2(); h_uPerp->Sumw2(); h_uPar_qT->Sumw2(); h_uPerp_qT->Sumw2();

  // -------- Loop --------
  const Long64_t nEntries = tEle->GetEntries();
  std::cout << "EventTree entries = " << nEntries << "\n";
  std::cout << "Have HiTree? " << (haveHiTree ? "YES" : "NO") << "\n";
  if (haveHiTree)
  {
    std::cout << "HiTree has vz? "    << (has_vz    ? "YES" : "NO") << "\n";
    std::cout << "HiTree has hiBin? " << (has_hiBin ? "YES" : "NO") << "\n";
  }

  Long64_t nPassEvent = 0;
  Long64_t nPassPair  = 0;
  bool warnedEventFiltersOnce = false;

  for (Long64_t ie = 0; ie < nEntries; ++ie)
  {
    if (ie % 200000 == 0) std::cout << "Event " << ie << "/" << nEntries << "\n";

    tEle->GetEntry(ie);
    if (haveHiTree) tHi->GetEntry(ie);
    if (haveRecoil) tPF->GetEntry(ie);
    tHLT->GetEntry(ie);
    tHLTobj->GetEntry(ie);
    tEvent->GetEntry(ie);

    const double w = has_genWeight ? (double)genWeight : 1.0;

    if (applyVz && has_vz && TMath::Abs(vz) > vzMax) continue;

    if (!PassEventSelection_pO(warnedEventFiltersOnce,
                               has_pprimaryVertexFilter,        pprimaryVertexFilter,
                               has_pclusterCompatibilityFilter, pclusterCompatibilityFilter))
      continue;

    if (applyHiBin && has_hiBin && (hiBin < hiBinMin || hiBin > hiBinMax)) continue;

    if (nEle < 2) continue;

    if (has_hlt && !TriggerFired(HLT_OxyL1SingleEG10_v1)) continue;
    nPassEvent++;

    for (int i = 0; i < nEle; ++i)
    {
      bool passIsolead = true;

      if (!elePt || !eleEta || !elePhi || !eleCharge) continue;
      if (elePt->at(i) < ptMin2) continue;
      if (TMath::Abs(eleEta->at(i)) > etaMax) continue;

      if (requireTightID && has_eleID && eleMVAIdWP95->at(i) == 0) continue;

      /* WP-based iso, kept for cross-check
      if (eleMVAIsoWP95->at(i) != 1) passIsolead = false;
      */
      const double isoLead = RelIsoPF(i, elePt, elePFChIso, elePFNeuIso, elePFPhoIso, elePFPUIso);
      if (isoLead >= isoMax) passIsolead = false;
      // FIXME: redo electron isolation study; not using a working point.

      for (int j = i + 1; j < nEle; ++j)
      {
        bool passIsosec = true;

        if (elePt->at(j) < ptMin2) continue;
        if (TMath::Abs(eleEta->at(j)) > etaMax) continue;
        if (requireTightID && has_eleID && eleMVAIdWP95->at(j) == 0) continue;

        if (requireOS && eleCharge->at(i) * eleCharge->at(j) >= 0) continue;

        const double pt1  = elePt->at(i);
        const double pt2  = elePt->at(j);
        const double lead = (pt1 > pt2 ? pt1 : pt2);
        const double sub  = (pt1 > pt2 ? pt2 : pt1);
        if (lead < ptMin1) continue;
        if (sub  < ptMin2) continue;

        TLorentzVector v1, v2;
        v1.SetPtEtaPhiM(elePt->at(i), eleEta->at(i), elePhi->at(i), ELE_MASS);
        v2.SetPtEtaPhiM(elePt->at(j), eleEta->at(j), elePhi->at(j), ELE_MASS);
        const double m = (v1 + v2).M();

        /* WP-based iso, kept for cross-check
        if (eleMVAIsoWP95->at(j) != 1) passIsosec = false;
        */
        const double isosub = RelIsoPF(j, elePt, elePFChIso, elePFNeuIso, elePFPhoIso, elePFPUIso);
        if (isosub >= isoMax) passIsosec = false;

        if (passIsolead && passIsosec)
        {
          hMass_extended->Fill(m, w);
          if (!(m < 60 || m > 120))
          {
            hMass->Fill(m, w);
            const TLorentzVector ll = v1 + v2;
            h_Zpt ->Fill(ll.Pt(),  w);
            h_Zeta->Fill(ll.Eta(), w);
            h_Zy  ->Fill(ll.Rapidity(), w);
            h_Zphi->Fill(ll.Phi(), w);
            h_lepPt ->Fill(v1.Pt(),  w); h_lepPt ->Fill(v2.Pt(),  w);
            h_lepEta->Fill(v1.Eta(), w); h_lepEta->Fill(v2.Eta(), w);
            h_lepPhi->Fill(v1.Phi(), w); h_lepPhi->Fill(v2.Phi(), w);

            // Hadronic recoil: u = -MET - q_T, with q_T = the dielectron (ll).
            // u_par is along q_T (should peak near -q_T), u_perp is transverse.
            if (haveRecoil)
            {
              const TVector2 metv = ComputePFMET(pfId, pfPt, pfPhi);
              const double qT = ll.Pt();
              if (qT > 0.0)
              {
                const double ux = -metv.X() - ll.Px();
                const double uy = -metv.Y() - ll.Py();
                const double cphi = ll.Px() / qT, sphi = ll.Py() / qT;
                const double uPar  =  ux * cphi + uy * sphi;
                const double uPerp = -ux * sphi + uy * cphi;
                h_uPar    ->Fill(uPar,  w);
                h_uPerp   ->Fill(uPerp, w);
                h_uPar_qT ->Fill(qT, uPar,  w);
                h_uPerp_qT->Fill(qT, uPerp, w);
              }
            }
          }
        }

        if (m < 60 || m > 120) continue;
        if (lead < 20 || sub < 20) continue;

        // Tighter pT cut, no iso requirement.
        hMass_vipul->Fill(m, w);
        nPassPair++;
      }
    }
  }

  std::cout << "Passed event preselection: " << nPassEvent << "\n";
  std::cout << "Found selected OS pair: "    << nPassPair  << "\n";

  gSystem->mkdir("rootfile", kTRUE); // ensure ./rootfile exists (fresh checkout)
  TFile *fout = new TFile(("./rootfile/" + outPrefix + "_" + mcTag + "_hist.root").c_str(), "RECREATE");
  hMass->Write("", 2);
  hMass_extended->Write("", 2);
  hMass_vipul->Write("", 2);
  h_Zpt->Write("", 2);   h_Zeta->Write("", 2);   h_Zphi->Write("", 2);   h_Zy->Write("", 2);
  h_lepPt->Write("", 2); h_lepEta->Write("", 2); h_lepPhi->Write("", 2);
  h_uPar->Write("", 2); h_uPerp->Write("", 2);
  h_uPar_qT->Write("", 2); h_uPerp_qT->Write("", 2);
  fout->Close();

  std::cout << "[INFO] Wrote outputs with prefix: " << outPrefix << "\n";
  return 0;
}


// ============================================================
// Dispatcher
// ============================================================

int skim(ChannelType channel,
         const char *dataFile = pOSkim::kDefaultDataFile,
         SampleType  sample   = kData)
{
  // Defensive: turn on Sumw2() for every TH1 created in this process.
  // Every active histogram below ALSO calls Sumw2() explicitly, so this is
  // belt-and-braces — guards against future additions forgetting. Idempotent.
  // Critical the moment any Fill() in this file starts passing an event
  // weight (lepton SFs, recoil correction, etc.): without Sumw2, errors on
  // weighted bins are silently wrong; with it, ROOT tracks sum-of-w^2 correctly.
  TH1::SetDefaultSumw2(kTRUE);

  switch (channel)
  {
    case kWmu: return skim_Wmu(dataFile, sample);
    case kWel: return skim_Wel(dataFile, sample);
    case kZmm: return skim_Zmm(dataFile, sample);
    case kZee: return skim_Zee(dataFile, sample);
  }
  std::cerr << "[ERR] skim(): unknown channel " << (int)channel << "\n";
  return 2;
}
