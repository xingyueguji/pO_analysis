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
                           std::vector<float> *muPFPhoIso)
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
    if (RelIsoPF(i, muPt, muPFChIso, muPFNeuIso, muPFPhoIso) >= isoMax) continue;
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

// Step 4 (electron W): OS dielectron DY veto. Uses integer WP gates
// (eleIDTight + eleIso) — different surface from the muon version.
static bool PassDYVeto_Wel(double dyElePtMin,
                           double dyMassMin, double dyMassMax,
                           int nEle,
                           std::vector<float> *elePt,
                           std::vector<float> *eleEta,
                           std::vector<float> *elePhi,
                           std::vector<int>   *eleCharge,
                           std::vector<int>   *eleIDTight,
                           std::vector<int>   *eleIso)
{
  if (!elePt || !eleEta || !elePhi || !eleCharge) return true;

  std::vector<int> cand;
  cand.reserve(nEle);

  for (int i = 0; i < nEle; ++i)
  {
    if (elePt->at(i) <= dyElePtMin) continue;
    if (eleIDTight->at(i) == 0)      continue;
    if (eleIso->at(i)    == 0)       continue;
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
  const double dyMassMin   = 80.0;
  const double dyMassMax   = 110.0;
  const double muEtaMax    = 2.4;
  const double trigMatchDR = 0.1;
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
  if (has_muPFChIso)  tMu->SetBranchAddress("muPFChIso",  &muPFChIso);
  if (has_muPFNeuIso) tMu->SetBranchAddress("muPFNeuIso", &muPFNeuIso);
  if (has_muPFPhoIso) tMu->SetBranchAddress("muPFPhoIso", &muPFPhoIso);

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

    if (!ExistsPFMuonPt25_W(muPt25, nMu, muPt, has_muIsPF, muIsPF)) continue;
    N[1]++;

    if (!PassEventSelection_pO(warnedEventFiltersOnce,
                               has_pprimaryVertexFilter,        pprimaryVertexFilter,
                               has_pclusterCompatibilityFilter, pclusterCompatibilityFilter))
      continue;
    if (applyVz && TMath::Abs(vz) > vzMax) continue;
    N[2]++;

    if (has_hlt && !TriggerFired(HLT_OxyL1SingleMuOpen_v1)) continue;
    N[3]++;

    if (!PassDYVeto_Wmu(dyMuPtMin, isoMax, dyMassMin, dyMassMax,
                        nMu, muPt, muEta, muPhi, muCharge,
                        has_muIDTight, muIDTight,
                        has_muIsPF, muIsPF,
                        muPFChIso, muPFNeuIso, muPFPhoIso))
      continue;
    N[4]++;

    if (!ExistsTightMuon_W(nMu, has_muIDTight, muIDTight)) continue;
    N[5]++;

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

    if (muPt->at(iLead) <= muPt25) continue;
    if (std::abs(muEta->at(iLead)) > muEtaMax) continue;
    N[6]++;

    const double isoLead = RelIsoPF(iLead, muPt, muPFChIso, muPFNeuIso, muPFPhoIso);
    const int    isoBin  = FindIsoBin(isoLead);
    const bool   isQCDSideband   = (isoBin >= 0);
    const bool   passIsoNominal  = (isoLead < isoMax);

    if (passIsoNominal) N[7]++;

    const bool passMatch = PassLeadingLeptonTrigMatch(
        trigMatchDR, iLead, muEta, muPhi,
        has_trgObjPt, trgObjPt, has_trgObjEta, trgObjEta,
        has_trgObjPhi, trgObjPhi,
        warnedNoTrigMatchInfo, "muon");
    if (!passMatch) continue;
    if (passIsoNominal) N[8]++;

    // -------- Fill final distributions (after step 8) --------
    TVector2 metv = ComputePFMET(pfId, pfPt, pfPhi);
    const double met = metv.Mod();
    const double mt  = TransverseMass(muPt->at(iLead), muPhi->at(iLead), metv);

    if (isQCDSideband)
    {
      if      (muCharge->at(iLead) > 0) h_met_iso_muPlus [isoBin]->Fill(met, w);
      else if (muCharge->at(iLead) < 0) h_met_iso_muMinus[isoBin]->Fill(met, w);
    }

    if (!passIsoNominal) continue;

    const int q = muCharge->at(iLead);
    const bool isWp = (q > 0), isWm = (q < 0);

    // [NOTICE] Sign flip: p-going (-Z) defined as forward
    const double y       = -muEta->at(iLead);
    const int    ybin    = FindYBin(y, kYEdges,   kNY);
    const int    ybin_FB = FindYBin(y, kYEdgesFB, kNY);

    if (ybin >= 0)
    {
      if      (isWp) { h_met_Wp[ybin]->Fill(met, w); h_mt_Wp[ybin]->Fill(mt, w); }
      else if (isWm) { h_met_Wm[ybin]->Fill(met, w); h_mt_Wm[ybin]->Fill(mt, w); }
    }
    if (ybin_FB >= 0)
    {
      if      (isWp) { h_met_Wp_FB[ybin_FB]->Fill(met, w); h_mt_Wp_FB[ybin_FB]->Fill(mt, w); }
      else if (isWm) { h_met_Wm_FB[ybin_FB]->Fill(met, w); h_mt_Wm_FB[ybin_FB]->Fill(mt, w); }
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
  }
  for (int b = 0; b < NISO; ++b)
  {
    h_met_iso_muPlus [b]->Write("", TObject::kOverwrite);
    h_met_iso_muMinus[b]->Write("", TObject::kOverwrite);
  }
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
  const double dyMassMin   = 80.0;
  const double dyMassMax   = 110.0;
  const double eleEtaMax   = 2.4;
  const double trigMatchDR = 0.1;
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

  const bool has_elePFChIso  = HasBranch(tEle, "elePFChIso");
  const bool has_elePFNeuIso = HasBranch(tEle, "elePFNeuIso");
  const bool has_elePFPhoIso = HasBranch(tEle, "elePFPhoIso");
  if (has_elePFChIso)  tEle->SetBranchAddress("elePFChIso",  &elePFChIso);
  if (has_elePFNeuIso) tEle->SetBranchAddress("elePFNeuIso", &elePFNeuIso);
  if (has_elePFPhoIso) tEle->SetBranchAddress("elePFPhoIso", &elePFPhoIso);

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

    if (!ExistsPFElectronPt25_W(elePt25, nEle, elePt)) continue;
    N[1]++;

    if (!PassEventSelection_pO(warnedEventFiltersOnce,
                               has_pprimaryVertexFilter,        pprimaryVertexFilter,
                               has_pclusterCompatibilityFilter, pclusterCompatibilityFilter))
      continue;
    if (applyVz && TMath::Abs(vz) > vzMax) continue;
    N[2]++;

    if (has_hlt && !TriggerFired(HLT_OxyL1SingleEG10_v1)) continue;
    N[3]++;

    if (!PassDYVeto_Wel(dyElePtMin, dyMassMin, dyMassMax,
                        nEle, elePt, eleEta, elePhi, eleCharge,
                        eleCutIdWP95, eleMVAIsoWP95))
      continue;
    N[4]++;

    if (!ExistsTightElectron_W(nEle, eleCutIdWP95)) continue;
    N[5]++;

    const int iLead = FindLeadingElectron_TightPF(nEle, elePt, eleEta, elePhi, eleCutIdWP95);
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

    if (elePt->at(iLead) <= elePt25) continue;
    if (std::abs(eleEta->at(iLead)) > eleEtaMax) continue;
    N[6]++;

    // (7) Leading electron Iso (continuous PF rel-iso; legacy WP path kept commented):
    // if (eleMVAIsoWP95->at(iLead) != 1) continue;
    const double isoLead = RelIsoPF(iLead, elePt, elePFChIso, elePFNeuIso, elePFPhoIso);
    const int    isoBin  = FindIsoBin(isoLead);
    const bool   isQCDSideband   = (isoBin >= 0);
    const bool   passIsoNominal  = (isoLead < isoMax);
    // FIXME: Need to redo electron isolation study, not using working point.

    if (passIsoNominal) N[7]++;

    const bool passMatch = PassLeadingLeptonTrigMatch(
        trigMatchDR, iLead, eleEta, elePhi,
        has_trgObjPt, trgObjPt, has_trgObjEta, trgObjEta,
        has_trgObjPhi, trgObjPhi,
        warnedNoTrigMatchInfo, "electron");
    if (!passMatch) continue;
    if (passIsoNominal) N[8]++;

    // -------- Fill final distributions (after step 8) --------
    TVector2 metv = ComputePFMET(pfId, pfPt, pfPhi);
    const double met = metv.Mod();
    const double mt  = TransverseMass(elePt->at(iLead), elePhi->at(iLead), metv);

    if (isQCDSideband)
    {
      if      (eleCharge->at(iLead) > 0) h_met_iso_elePlus [isoBin]->Fill(met, w);
      else if (eleCharge->at(iLead) < 0) h_met_iso_eleMinus[isoBin]->Fill(met, w);
    }

    if (!passIsoNominal) continue;

    const int q = eleCharge->at(iLead);
    const bool isWp = (q > 0), isWm = (q < 0);

    // [NOTICE] Sign flip: p-going (-Z) defined as forward
    const double y       = -eleEta->at(iLead);
    const int    ybin    = FindYBin(y, kYEdges,   kNY);
    const int    ybin_FB = FindYBin(y, kYEdgesFB, kNY);

    if (ybin >= 0)
    {
      if      (isWp) { h_met_Wp[ybin]->Fill(met, w); h_mt_Wp[ybin]->Fill(mt, w); }
      else if (isWm) { h_met_Wm[ybin]->Fill(met, w); h_mt_Wm[ybin]->Fill(mt, w); }
    }
    if (ybin_FB >= 0)
    {
      if      (isWp) { h_met_Wp_FB[ybin_FB]->Fill(met, w); h_mt_Wp_FB[ybin_FB]->Fill(mt, w); }
      else if (isWm) { h_met_Wm_FB[ybin_FB]->Fill(met, w); h_mt_Wm_FB[ybin_FB]->Fill(mt, w); }
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
  }
  for (int b = 0; b < NISO; ++b)
  {
    h_met_iso_elePlus [b]->Write("", TObject::kOverwrite);
    h_met_iso_eleMinus[b]->Write("", TObject::kOverwrite);
  }
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

  // -------- Muon branches --------
  Int_t nMu = 0;
  std::vector<float> *muPt = nullptr, *muEta = nullptr, *muPhi = nullptr;
  std::vector<int>   *muCharge = nullptr;
  std::vector<int>   *muIsGood = nullptr, *muIsGlobal = nullptr;
  std::vector<int>   *muIsPF   = nullptr, *muIDTight  = nullptr;

  tMu->SetBranchAddress("nMu",      &nMu);
  tMu->SetBranchAddress("muPt",     &muPt);
  tMu->SetBranchAddress("muEta",    &muEta);
  tMu->SetBranchAddress("muPhi",    &muPhi);
  tMu->SetBranchAddress("muCharge", &muCharge);

  const bool has_muIsGood   = HasBranch(tMu, "muIsGood");
  const bool has_muIsGlobal = HasBranch(tMu, "muIsGlobal");
  const bool has_muIsPF     = HasBranch(tMu, "muIsPF");
  const bool has_muIDTight  = HasBranch(tMu, "muIDTight");

  if (has_muIsGood)   tMu->SetBranchAddress("muIsGood",   &muIsGood);
  if (has_muIsGlobal) tMu->SetBranchAddress("muIsGlobal", &muIsGlobal);
  if (has_muIsPF)     tMu->SetBranchAddress("muIsPF",     &muIsPF);
  if (has_muIDTight)  tMu->SetBranchAddress("muIDTight",  &muIDTight);

  // -------- Event branches --------
  Float_t vz    = 999.f;
  Int_t   hiBin = -999;

  const bool has_vz    = (haveHiTree && HasBranch(tHi, "vz"));
  const bool has_hiBin = (haveHiTree && HasBranch(tHi, "hiBin"));

  if (has_vz)    tHi->SetBranchAddress("vz",    &vz);
  if (has_hiBin) tHi->SetBranchAddress("hiBin", &hiBin);

  // -------- Generator event weight (MC only) --------
  Float_t    genWeight     = 1.f;
  const bool has_genWeight = isMC && haveHiTree && HasBranch(tHi, "weight");
  if (has_genWeight) tHi->SetBranchAddress("weight", &genWeight);
  else if (isMC)     std::cout << "[WARN] skim_Zmm: MC sample but no 'weight' branch on HiTree; filling unweighted.\n";

  // -------- PF iso --------
  std::vector<float> *muPFChIso = nullptr, *muPFNeuIso = nullptr, *muPFPhoIso = nullptr;

  const bool has_muPFChIso  = HasBranch(tMu, "muPFChIso");
  const bool has_muPFNeuIso = HasBranch(tMu, "muPFNeuIso");
  const bool has_muPFPhoIso = HasBranch(tMu, "muPFPhoIso");

  if (has_muPFChIso)  { tMu->SetBranchStatus("muPFChIso",  1); tMu->SetBranchAddress("muPFChIso",  &muPFChIso);  }
  if (has_muPFNeuIso) { tMu->SetBranchStatus("muPFNeuIso", 1); tMu->SetBranchAddress("muPFNeuIso", &muPFNeuIso); }
  if (has_muPFPhoIso) { tMu->SetBranchStatus("muPFPhoIso", 1); tMu->SetBranchAddress("muPFPhoIso", &muPFPhoIso); }

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
  TH1D *h_Zphi  = new TH1D("h_Zphi", Form("%s; #phi_{#mu#mu}; Events",        outPrefix.c_str()), 32, -TMath::Pi(), TMath::Pi());
  TH1D *h_lepPt  = new TH1D("h_lepPt",  Form("%s; p_{T}^{#mu} [GeV]; Muons", outPrefix.c_str()), 50,  0,  100);
  TH1D *h_lepEta = new TH1D("h_lepEta", Form("%s; #eta_{#mu}; Muons",        outPrefix.c_str()), 50, -2.5, 2.5);
  TH1D *h_lepPhi = new TH1D("h_lepPhi", Form("%s; #phi_{#mu}; Muons",        outPrefix.c_str()), 32, -TMath::Pi(), TMath::Pi());
  h_Zpt->Sumw2();   h_Zeta->Sumw2();   h_Zphi->Sumw2();
  h_lepPt->Sumw2(); h_lepEta->Sumw2(); h_lepPhi->Sumw2();

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

      const double isoLead = RelIsoPF(i, muPt, muPFChIso, muPFNeuIso, muPFPhoIso);
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

        const double isosub = RelIsoPF(j, muPt, muPFChIso, muPFNeuIso, muPFPhoIso);
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
            h_Zphi->Fill(ll.Phi(), w);
            h_lepPt ->Fill(v1.Pt(),  w); h_lepPt ->Fill(v2.Pt(),  w);
            h_lepEta->Fill(v1.Eta(), w); h_lepEta->Fill(v2.Eta(), w);
            h_lepPhi->Fill(v1.Phi(), w); h_lepPhi->Fill(v2.Phi(), w);
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
  h_Zpt->Write("", 2);   h_Zeta->Write("", 2);   h_Zphi->Write("", 2);
  h_lepPt->Write("", 2); h_lepEta->Write("", 2); h_lepPhi->Write("", 2);
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
  const bool   requireTightID = true;  // active gate uses eleCutIdWP95 (FIXME: WP95 is loose)
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

  // -------- Electron kinematics --------
  Int_t nEle = 0;
  std::vector<float> *elePt = nullptr, *eleEta = nullptr, *elePhi = nullptr;
  std::vector<int>   *eleCharge = nullptr;
  tEle->SetBranchAddress("nEle",      &nEle);
  tEle->SetBranchAddress("elePt",     &elePt);
  tEle->SetBranchAddress("eleEta",    &eleEta);
  tEle->SetBranchAddress("elePhi",    &elePhi);
  tEle->SetBranchAddress("eleCharge", &eleCharge);

  // -------- Electron IDs (kept around; active gate uses eleCutIdWP95) --------
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

  const bool has_eleID = HasBranch(tEle, "eleCutIdWP95");

  // -------- Iso WPs (kept around; active gate uses RelIsoPF) --------
  std::vector<int> *eleMVAIsoWP80 = nullptr, *eleMVAIsoWP85 = nullptr;
  std::vector<int> *eleMVAIsoWP90 = nullptr, *eleMVAIsoWP95 = nullptr;
  hookID("eleMVAIsoWP80", &eleMVAIsoWP80);
  hookID("eleMVAIsoWP85", &eleMVAIsoWP85);
  hookID("eleMVAIsoWP90", &eleMVAIsoWP90);
  hookID("eleMVAIsoWP95", &eleMVAIsoWP95);

  std::vector<float> *elePFChIso = nullptr, *elePFNeuIso = nullptr, *elePFPhoIso = nullptr;
  const bool has_elePFChIso  = HasBranch(tEle, "elePFChIso");
  const bool has_elePFNeuIso = HasBranch(tEle, "elePFNeuIso");
  const bool has_elePFPhoIso = HasBranch(tEle, "elePFPhoIso");
  if (has_elePFChIso)  { tEle->SetBranchStatus("elePFChIso", 1);  tEle->SetBranchAddress("elePFChIso",  &elePFChIso);  }
  if (has_elePFNeuIso) { tEle->SetBranchStatus("elePFNeuIso", 1); tEle->SetBranchAddress("elePFNeuIso", &elePFNeuIso); }
  if (has_elePFPhoIso) { tEle->SetBranchStatus("elePFPhoIso", 1); tEle->SetBranchAddress("elePFPhoIso", &elePFPhoIso); }

  // -------- Event branches --------
  Float_t vz    = 999.f;
  Int_t   hiBin = -999;
  const bool has_vz    = (haveHiTree && HasBranch(tHi, "vz"));
  const bool has_hiBin = (haveHiTree && HasBranch(tHi, "hiBin"));
  if (has_vz)    tHi->SetBranchAddress("vz",    &vz);
  if (has_hiBin) tHi->SetBranchAddress("hiBin", &hiBin);

  // -------- Generator event weight (MC only) --------
  Float_t    genWeight     = 1.f;
  const bool has_genWeight = isMC && haveHiTree && HasBranch(tHi, "weight");
  if (has_genWeight) tHi->SetBranchAddress("weight", &genWeight);
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
  TH1D *h_Zphi  = new TH1D("h_Zphi", Form("%s; #phi_{ee}; Events",        outPrefix.c_str()), 32, -TMath::Pi(), TMath::Pi());
  TH1D *h_lepPt  = new TH1D("h_lepPt",  Form("%s; p_{T}^{e} [GeV]; Electrons", outPrefix.c_str()), 50,  0,  100);
  TH1D *h_lepEta = new TH1D("h_lepEta", Form("%s; #eta_{e}; Electrons",        outPrefix.c_str()), 50, -2.5, 2.5);
  TH1D *h_lepPhi = new TH1D("h_lepPhi", Form("%s; #phi_{e}; Electrons",        outPrefix.c_str()), 32, -TMath::Pi(), TMath::Pi());
  h_Zpt->Sumw2();   h_Zeta->Sumw2();   h_Zphi->Sumw2();
  h_lepPt->Sumw2(); h_lepEta->Sumw2(); h_lepPhi->Sumw2();

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

      if (requireTightID && has_eleID && eleCutIdWP95->at(i) == 0) continue;

      /* WP-based iso, kept for cross-check
      if (eleMVAIsoWP95->at(i) != 1) passIsolead = false;
      */
      const double isoLead = RelIsoPF(i, elePt, elePFChIso, elePFNeuIso, elePFPhoIso);
      if (isoLead >= isoMax) passIsolead = false;
      // FIXME: redo electron isolation study; not using a working point.

      for (int j = i + 1; j < nEle; ++j)
      {
        bool passIsosec = true;

        if (elePt->at(j) < ptMin2) continue;
        if (TMath::Abs(eleEta->at(j)) > etaMax) continue;
        if (requireTightID && has_eleID && eleCutIdWP95->at(j) == 0) continue;

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
        const double isosub = RelIsoPF(j, elePt, elePFChIso, elePFNeuIso, elePFPhoIso);
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
            h_Zphi->Fill(ll.Phi(), w);
            h_lepPt ->Fill(v1.Pt(),  w); h_lepPt ->Fill(v2.Pt(),  w);
            h_lepEta->Fill(v1.Eta(), w); h_lepEta->Fill(v2.Eta(), w);
            h_lepPhi->Fill(v1.Phi(), w); h_lepPhi->Fill(v2.Phi(), w);
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
  h_Zpt->Write("", 2);   h_Zeta->Write("", 2);   h_Zphi->Write("", 2);
  h_lepPt->Write("", 2); h_lepEta->Write("", 2); h_lepPhi->Write("", 2);
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
