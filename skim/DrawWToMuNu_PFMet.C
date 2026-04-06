// DrawWToMuNu_PFMet_CutflowTableStrict.C
//
// Goal: follow the cutflow logic EXACTLY as in the table:
//
//  0) All events in dataset
//  1) >=1 PF muon with pT > 25
//  2) pO collision event selection
//  3) Trigger PAL3Mu12 fired
//  4) Drell–Yan veto
//  5) >=1 Tight ID muon
//  6) Leading muon pT > 25
//  7) Leading muon Iso < 0.15
//  8) Leading muon matched to trigger
//
// Notes / assumptions (kept as “auto-detect” wherever possible):
// - “PF muon” is taken as muIsPF == 1 when the branch exists.
// - “Tight ID muon” is taken as muIDTight == 1 when the branch exists.
// - Event selection (step 2): we AND together common pO/pA filters IF they exist;
//   if a given filter branch is missing, we don’t apply it (and we print a warning once).
// - Trigger fired (step 3): we auto-find any branch name containing "HLT_PAL3Mu12"
//   in ggHiNtuplizer/EventTree and use it.
// - Trigger matching (step 8): we try, in order:
//    (A) trigger object collections (trgObj* style) if present
//    (B) per-muon trigger match bits (muTrig/muTrg style) if present
//   If neither exists, we can’t enforce step 8 and we will WARN and treat as pass.
// - DY veto (step 4): veto event if exists an OS muon pair with BOTH muons:
//      pT>15, PF, Tight, relIso<0.15, and (mll > 30 GeV)   <-- your requested m>30
//
// Run:
//   root -l -q 'DrawWToMuNu_PFMet_CutflowTableStrict.C("root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2025.root")'
//
// Output:
//   - Prints cutflow counts + efficiencies
//   - Saves MET and mT after the FINAL selection (step 8)

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

#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>

enum SampleType
{
  kData,
  kDY,
  kWp,
  kWm
};

// ---------------------------------------------
// Small utilities
// ---------------------------------------------
static bool HasBranch(TTree *t, const char *name)
{
  return (t && t->GetListOfBranches() && t->GetListOfBranches()->FindObject(name));
}

static std::string FindBranchContaining(TTree *t, const std::string &needle)
{
  if (!t || !t->GetListOfBranches())
    return "";
  TObjArray *arr = t->GetListOfBranches();
  const int n = arr->GetEntries();
  for (int i = 0; i < n; ++i)
  {
    TBranch *br = (TBranch *)arr->At(i);
    if (!br)
      continue;
    std::string bn = br->GetName();
    if (bn.find(needle) != std::string::npos)
      return bn;
  }
  return "";
}

static double DeltaR(double eta1, double phi1, double eta2, double phi2)
{
  const double dEta = eta1 - eta2;
  const double dPhi = TVector2::Phi_mpi_pi(phi1 - phi2);
  return std::sqrt(dEta * dEta + dPhi * dPhi);
}

// ---------------------------------------------
// Config
// ---------------------------------------------
struct WConfig
{
  // Table parameters
  double muPt25 = 25.0;
  double dyMuPtMin = 15.0;
  double isoMax = 0.15;
  double dyMassMin = 80.0; // your requested M>30 inside DY veto
  double dyMassMax = 110.0;
  double muEtaMax = 2.4;

  // For trigger-object matching
  double trigMatchDR = 0.1;
  bool applyVz = true;
  double vzMax = 15.0;

  // Histograms
  int nBinsMt = 100;
  double mtMin = 0.0;
  double mtMax = 200.0;

  std::string outPrefix = "WToMuNu_pO_PFMet";
};

// ---------------------------------------------
// PF MET from pftree
// ---------------------------------------------
static TVector2 ComputePFMET(std::vector<int> *pfId,
                             std::vector<float> *pfPt,
                             std::vector<float> *pfPhi)
{
  (void)pfId;
  double sumPx = 0.0, sumPy = 0.0;
  const int nPF = (pfPt ? (int)pfPt->size() : 0);

  for (int i = 0; i < nPF; ++i)
  {
    const double pt = pfPt->at(i);
    const double phi = pfPhi->at(i);
    sumPx += pt * std::cos(phi);
    sumPy += pt * std::sin(phi);
  }
  return TVector2(-sumPx, -sumPy);
}

static double RelIsoPF(int i,
                       std::vector<float> *muPt,
                       std::vector<float> *muPFChIso,
                       std::vector<float> *muPFNeuIso,
                       std::vector<float> *muPFPhoIso)
{
  if (!muPt || !muPFChIso || !muPFNeuIso || !muPFPhoIso)
    return 999.0;
  const double pt = muPt->at(i);
  if (pt <= 0)
    return 999.0;
  const double isoAbs = muPFChIso->at(i) + muPFNeuIso->at(i) + muPFPhoIso->at(i);
  return isoAbs / pt;
}

// ---------------------------------------------
// Step 2: event selection (auto-detect common filter flags if present)
// We AND all existing filters in the list.
// If a filter branch is missing, we do not apply it, but we warn once.
// ---------------------------------------------
static bool PassEventSelection_pO(TTree *tHi,
                                  bool &warnedOnce,
                                  // common skim/filter flags (Int_t or Bool_t stored as Int_t often)
                                  bool has_pprimaryVertexFilter, int pprimaryVertexFilter,
                                  bool has_pclusterCompatibilityFilter, int pclusterCompatibilityFilter)
{
  (void)tHi;

  auto requireIfExists = [&](bool has, int val, const char *name) -> bool
  {
    if (!has)
    {
      if (!warnedOnce)
        std::cout << "[WARN] Missing event filter branch: " << name << " (will not apply)\n";
      return true;
    }
    return (val != 0);
  };

  bool ok = true;
  ok = ok && requireIfExists(has_pprimaryVertexFilter, pprimaryVertexFilter, "pprimaryVertexFilter");
  ok = ok && requireIfExists(has_pclusterCompatibilityFilter, pclusterCompatibilityFilter, "pclusterCompatibilityFilter");

  warnedOnce = true; // after first call, suppress repeated warnings
  return ok;
}

// ---------------------------------------------
// Step 1 helper: exists PF muon with pt>25
// ---------------------------------------------
static bool ExistsPFMuonPt25(const WConfig &cfg,
                             int nMu,
                             std::vector<float> *muPt,
                             bool has_muIsPF, std::vector<int> *muIsPF)
{
  if (!muPt)
    return false;

  for (int i = 0; i < nMu; ++i)
  {
    const double pt = muPt->at(i);
    if (pt <= cfg.muPt25)
      continue;

    // PF requirement only if branch exists; if missing we can't enforce PF, so treat as pass
    if (has_muIsPF && muIsPF && muIsPF->at(i) == 0)
      continue;

    return true;
  }
  return false;
}

// ---------------------------------------------
// Step 5 helper: exists Tight ID muon (no pt requirement in table)
// ---------------------------------------------
static bool ExistsTightMuon(int nMu,
                            bool has_muIDTight, std::vector<int> *muIDTight)
{
  if (!has_muIDTight || !muIDTight)
    return false; // if table requires tight, and branch absent, we should fail rather than silently pass

  for (int i = 0; i < nMu; ++i)
  {
    if (muIDTight->at(i) != 0)
      return true;
  }
  return false;
}

// ---------------------------------------------
// Leading muon definition for steps 6-8:
// We'll define "leading muon" among muons that are Tight (if exists) and PF (if exists).
// This is the closest to typical W analysis and makes steps 6-8 consistent.
// ---------------------------------------------
static int FindLeadingMuon_TightPF(int nMu,
                                   std::vector<float> *muPt,
                                   std::vector<float> *muEta,
                                   std::vector<float> *muPhi,
                                   bool has_muIDTight, std::vector<int> *muIDTight,
                                   bool has_muIsPF, std::vector<int> *muIsPF)
{
  if (!muPt || !muEta || !muPhi)
    return -1;

  int best = -1;
  double bestPt = -1;

  for (int i = 0; i < nMu; ++i)
  {
    const double pt = muPt->at(i);

    // Tight requirement if branch exists; if missing, cannot define “tight leading muon”
    if (has_muIDTight && muIDTight && muIDTight->at(i) == 0)
      continue;

    // PF requirement if branch exists
    if (has_muIsPF && muIsPF && muIsPF->at(i) == 0)
      continue;

    if (pt > bestPt)
    {
      bestPt = pt;
      best = i;
    }
  }
  return best;
}

// ---------------------------------------------
// Step 4: DY veto (with Mll > 30 requirement)
// Veto event if exists OS pair where both muons satisfy:
//   pT > 15, PF, Tight, relIso < 0.15, AND mll > 30
// ---------------------------------------------
static bool PassDYVeto(const WConfig &cfg,
                       int nMu,
                       std::vector<float> *muPt,
                       std::vector<float> *muEta,
                       std::vector<float> *muPhi,
                       std::vector<int> *muCharge,
                       bool has_muIDTight, std::vector<int> *muIDTight,
                       bool has_muIsPF, std::vector<int> *muIsPF,
                       std::vector<float> *muPFChIso,
                       std::vector<float> *muPFNeuIso,
                       std::vector<float> *muPFPhoIso)
{
  if (!muPt || !muEta || !muPhi || !muCharge)
    return true; // can't apply -> do not veto

  // If tight branch missing, we can't apply the DY veto as defined (table says "DY veto" though),
  // but to avoid silently over-vetoing/under-vetoing, we choose: if missing, do not veto.
  if (!has_muIDTight || !muIDTight)
    return true;

  std::vector<int> cand;
  cand.reserve(nMu);

  for (int i = 0; i < nMu; ++i)
  {
    const double pt = muPt->at(i);
    if (pt <= cfg.dyMuPtMin)
      continue;

    if (muIDTight->at(i) == 0)
      continue;

    if (has_muIsPF && muIsPF && muIsPF->at(i) == 0)
      continue;

    const double isoRel = RelIsoPF(i, muPt, muPFChIso, muPFNeuIso, muPFPhoIso);
    if (isoRel >= cfg.isoMax)
      continue;

    cand.push_back(i);
  }

  const double muMass = 0.105658;

  for (size_t a = 0; a < cand.size(); ++a)
  {
    for (size_t b = a + 1; b < cand.size(); ++b)
    {
      const int i1 = cand[a];
      const int i2 = cand[b];

      if (muCharge->at(i1) * muCharge->at(i2) >= 0)
        continue; // require OS

      TLorentzVector mu1, mu2;
      mu1.SetPtEtaPhiM(muPt->at(i1), muEta->at(i1), muPhi->at(i1), muMass);
      mu2.SetPtEtaPhiM(muPt->at(i2), muEta->at(i2), muPhi->at(i2), muMass);

      const double mll = (mu1 + mu2).M();

      // DY-like OS pair with mll > 30 => veto event (fail DY veto)
      if (mll > cfg.dyMassMin && mll < cfg.dyMassMax)
        return false;
    }
  }

  return true;
}

// ---------------------------------------------
// Step 3: Trigger fired (PAL3Mu12) from HLT bit branch (auto-found)
// ---------------------------------------------
static bool TriggerFired(int hltBit)
{
  return (hltBit != 0);
}

// ---------------------------------------------
// Step 8: Leading muon matched to trigger
// Try trigger objects, else try per-muon match bits.
// If neither exists, warn and accept (cannot enforce).
// ---------------------------------------------
static bool PassLeadingMuonTrigMatch(const WConfig &cfg,
                                     int iLead,
                                     std::vector<float> *muEta,
                                     std::vector<float> *muPhi,

                                     // Option A: trigger object branches (if exist)
                                     bool has_trgObjPt, std::vector<double> *trgObjPt,
                                     bool has_trgObjEta, std::vector<double> *trgObjEta,
                                     bool has_trgObjPhi, std::vector<double> *trgObjPhi,
                                     bool has_trgObjId, std::vector<double> *trgObjId,
                                     bool &warnedNoTrigMatchInfo)
{
  if (iLead < 0 || !muEta || !muPhi)
    return false;

  const double eta = muEta->at(iLead);
  const double phi = muPhi->at(iLead);

  // A) Trigger objects
  if (has_trgObjPt && has_trgObjEta && has_trgObjPhi &&
      trgObjPt && trgObjEta && trgObjPhi)
  {
    const int nObj = (int)trgObjPt->size();
    for (int j = 0; j < nObj; ++j)
    {
      const double dr = DeltaR(eta, phi, trgObjEta->at(j), trgObjPhi->at(j));
      if (dr < cfg.trigMatchDR)
        return true;
    }
    return false;
  }

  // C) Unknown: cannot enforce
  if (!warnedNoTrigMatchInfo)
  {
    std::cout << "[WARN] No trigger-object collection and no per-muon trigger-match info found. "
              << "Cannot apply step 8 (leading muon matched to trigger). Treating as PASS.\n";
    warnedNoTrigMatchInfo = true;
  }
  return true;
}

// DecayID check

static bool HasAncestor(
    int idx,
    int targetPdg,
    std::vector<int> *genPdg,
    std::vector<std::vector<int>> *genMotherIdx)
{
  if (!genPdg || !genMotherIdx)
    return false;
  if (idx < 0 || (size_t)idx >= genPdg->size() || (size_t)idx >= genMotherIdx->size())
    return false;

  const auto &moms = genMotherIdx->at(idx);
  for (int midx : moms)
  {
    if (midx < 0 || (size_t)midx >= genPdg->size() || (size_t)midx >= genMotherIdx->size())
      continue;

    if (std::abs(genPdg->at(midx)) == std::abs(targetPdg))
      return true;

    if (HasAncestor(midx, targetPdg, genPdg, genMotherIdx))
      return true;
  }

  return false;
}

// ---------------------------------------------
// Extra: Gen-reco matching
// ---------------------------------------------

static bool PassGenRecoMatchingWithAncestor(
    int iLead,
    std::vector<float> *muPt, std::vector<float> *muEta, std::vector<float> *muPhi, std::vector<int> *muCharge,
    std::vector<float> *genPt, std::vector<float> *genEta, std::vector<float> *genPhi, std::vector<int> *genChg,
    std::vector<int> *genPdg,
    std::vector<std::vector<int>> *genMotherIdx,
    int requiredGenPdg,      // e.g. 13 for muon
    int requiredAncestorPdg) // e.g. 24 for W, 23 for Z
{
  if (!muPt || !muEta || !muPhi || !muCharge)
    return false;
  if (!genPt || !genEta || !genPhi || !genChg || !genPdg || !genMotherIdx)
    return false;
  if (iLead < 0 || (size_t)iLead >= muPt->size())
    return false;

  const double muonPt = muPt->at(iLead);
  const double muonEta = muEta->at(iLead);
  const double muonPhi = muPhi->at(iLead);
  const int muonChg = muCharge->at(iLead);

  const size_t ngen = std::min({genPt->size(),
                                genEta->size(),
                                genPhi->size(),
                                genChg->size(),
                                genPdg->size(),
                                genMotherIdx->size()});

  int bestIdx = -1;
  double bestDR = 1e9;

  for (size_t i = 0; i < ngen; ++i)
  {
    if (genPt->at(i) <= 0)
      continue;
    if (std::abs(genPdg->at(i)) != std::abs(requiredGenPdg))
      continue;
    if (genChg->at(i) != muonChg)
      continue;

    const double dEta = genEta->at(i) - muonEta;
    const double dPhi = TVector2::Phi_mpi_pi(genPhi->at(i) - muonPhi);
    const double dR = std::sqrt(dEta * dEta + dPhi * dPhi);
    const double dPtRel = std::fabs(genPt->at(i) - muonPt) / genPt->at(i);

    if (dR >= 0.5)
      continue;
    if (dPtRel >= 0.5)
      continue;

    if (dR < bestDR)
    {
      bestDR = dR;
      bestIdx = (int)i;
    }
  }

  if (bestIdx < 0)
    return false;

  /*std::cout << "[MatchDebug] matched reco muon iLead=" << iLead
            << " to gen idx=" << bestIdx
            << " genPdg=" << genPdg->at(bestIdx)
            << " genPt=" << genPt->at(bestIdx)
            << " dR=" << bestDR
            << "\n";

  const auto &moms = genMotherIdx->at(bestIdx);

  if (moms.empty())
  {
    std::cout << "[MatchDebug] matched gen idx=" << bestIdx
              << " has no stored mother\n";
    return true;
  }

  std::cout << "[MatchDebug] matched gen idx=" << bestIdx << " mothers:";
  for (int midx : moms)
  {
    std::cout << " " << midx;
    if (midx >= 0 && (size_t)midx < genPdg->size())
      std::cout << "(pdg=" << genPdg->at(midx) << ")";
    else
      std::cout << "(pdg=OUT_OF_LOCAL_RANGE)";
  }
  std::cout << "\n";*/

  return true;
}

// ---------------------------------------------
// Main
// ---------------------------------------------
void DrawWToMuNu_PFMet(const char *fname =
                           "root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_DATA_pass_2.root",
                       SampleType sample = kWm)
{

  gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");

  bool isMC = false;

  gStyle->SetOptStat(0);

  if (sample == kData)
  {
    isMC = false;
  }
  else
  {
    isMC = true;
  }

  WConfig cfg;

  // "root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2025.root"
  // root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_DATA_pass_2.root

  if (isMC)
  {
    if (sample == kDY)
    {
      fname = Form("root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_MC_DY_mu_Z.root");
      cfg.outPrefix += "_DY";
    }
    if (sample == kWp)
    {
      fname = Form("root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_MC_Wp_mu.root");
      cfg.outPrefix += "_Wp";
    }
    if (sample == kWm)
    {
      fname = Form("root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_MC_Wm_mu.root");
      cfg.outPrefix += "_Wm";
    }
  }

  TFile *f = TFile::Open(fname);
  if (!f || f->IsZombie())
  {
    std::cerr << "ERROR opening file\n";
    return;
  }

  TTree *tMu = (TTree *)f->Get("ggHiNtuplizer/EventTree");
  TTree *tHi = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  TTree *tPF = (TTree *)f->Get("particleFlowAnalyser/pftree");
  TTree *tHLT = (TTree *)f->Get("hltanalysis/HltTree");
  TTree *tHLTobj = (TTree *)f->Get("hltobject/HLT_OxyL1SingleMuOpen_v");
  TTree *tEvent = (TTree *)f->Get("skimanalysis/HltTree");
  TTree *tGen;

  if (isMC)
  {
    tGen = (TTree *)f->Get("HiGenParticleAna/hi");
  }

  if (!tMu || !tHi || !tPF || !tHLT || !tHLTobj)
  {
    std::cerr << "ERROR: missing one of trees: ggHiNtuplizer/EventTree, hiEvtAnalyzer/HiTree, particleFlowAnalyser/pftree\n";
    return;
  }

  // -------------------------
  // Muon branches
  // -------------------------
  Int_t nMu = 0;
  std::vector<float> *muPt = nullptr;
  std::vector<float> *muEta = nullptr;
  std::vector<float> *muPhi = nullptr;
  std::vector<int> *muCharge = nullptr;
  std::vector<int> *muIDTight = nullptr;
  std::vector<int> *muIsPF = nullptr;

  std::vector<float> *muPFChIso = nullptr;
  std::vector<float> *muPFNeuIso = nullptr;
  std::vector<float> *muPFPhoIso = nullptr;

  // Optional per-muon trigger match
  std::vector<int> *muTrig = nullptr;

  tMu->SetBranchStatus("*", 0);
  tMu->SetBranchStatus("nMu", 1);
  tMu->SetBranchStatus("muPt", 1);
  tMu->SetBranchStatus("muEta", 1);
  tMu->SetBranchStatus("muPhi", 1);
  tMu->SetBranchStatus("muCharge", 1);

  // ID/PF
  if (HasBranch(tMu, "muIDTight"))
    tMu->SetBranchStatus("muIDTight", 1);
  if (HasBranch(tMu, "muIsPF"))
    tMu->SetBranchStatus("muIsPF", 1);

  // Iso
  if (HasBranch(tMu, "muPFChIso"))
    tMu->SetBranchStatus("muPFChIso", 1);
  if (HasBranch(tMu, "muPFNeuIso"))
    tMu->SetBranchStatus("muPFNeuIso", 1);
  if (HasBranch(tMu, "muPFPhoIso"))
    tMu->SetBranchStatus("muPFPhoIso", 1);

  tMu->SetBranchAddress("nMu", &nMu);
  tMu->SetBranchAddress("muPt", &muPt);
  tMu->SetBranchAddress("muEta", &muEta);
  tMu->SetBranchAddress("muPhi", &muPhi);
  tMu->SetBranchAddress("muCharge", &muCharge);

  const bool has_muIDTight = HasBranch(tMu, "muIDTight");
  const bool has_muIsPF = HasBranch(tMu, "muIsPF");

  if (has_muIDTight)
    tMu->SetBranchAddress("muIDTight", &muIDTight);
  if (has_muIsPF)
    tMu->SetBranchAddress("muIsPF", &muIsPF);

  const bool has_muPFChIso = HasBranch(tMu, "muPFChIso");
  const bool has_muPFNeuIso = HasBranch(tMu, "muPFNeuIso");
  const bool has_muPFPhoIso = HasBranch(tMu, "muPFPhoIso");

  if (has_muPFChIso)
    tMu->SetBranchAddress("muPFChIso", &muPFChIso);
  if (has_muPFNeuIso)
    tMu->SetBranchAddress("muPFNeuIso", &muPFNeuIso);
  if (has_muPFPhoIso)
    tMu->SetBranchAddress("muPFPhoIso", &muPFPhoIso);

  // -------------------------
  // Trigger fired bit (step 3):
  // -------------------------
  const std::string hltName = FindBranchContaining(tHLT, "HLT_OxyL1SingleMuOpen_v1");
  Int_t HLT_OxyL1SingleMuOpen_v1 = 0;
  bool has_hlt = false;
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

  // -------------------------
  // Trigger objects (step 8): auto-detect common names
  // -------------------------
  std::vector<double> *trgObjPt = nullptr;
  std::vector<double> *trgObjEta = nullptr;
  std::vector<double> *trgObjPhi = nullptr;
  std::vector<double> *trgObjId = nullptr;

  const bool has_trgObjPt = HasBranch(tHLTobj, "pt");
  const bool has_trgObjEta = HasBranch(tHLTobj, "eta");
  const bool has_trgObjPhi = HasBranch(tHLTobj, "phi");
  const bool has_trgObjId = HasBranch(tHLTobj, "TriggerObjID");

  tHLTobj->SetBranchStatus("*", false);

  if (has_trgObjPt && has_trgObjEta && has_trgObjPhi && has_trgObjId)
  {
    tHLTobj->SetBranchStatus("pt", 1);
    tHLTobj->SetBranchStatus("eta", 1);
    tHLTobj->SetBranchStatus("phi", 1);
    tHLTobj->SetBranchStatus("TriggerObjID", 1);

    tHLTobj->SetBranchAddress("pt", &trgObjPt);
    tHLTobj->SetBranchAddress("eta", &trgObjEta);
    tHLTobj->SetBranchAddress("phi", &trgObjPhi);
    tHLTobj->SetBranchAddress("TriggerObjID", &trgObjId);
  }
  else
  {
    cout << "[WARN] Could not find HLT related Object." << endl;
  }

  // -------------------------
  // HiTree: event filters (step 2)
  // -------------------------
  Float_t vz = 999.f; // not used in table steps directly; keep available if you later want
  Int_t hiBin = -999; // not used

  // Common pO filters (names can vary by production)
  Int_t pprimaryVertexFilter = 1;
  Int_t pclusterCompatibilityFilter = 1;

  const bool has_pprimaryVertexFilter = HasBranch(tEvent, "pprimaryVertexFilter");
  const bool has_pclusterCompatibilityFilter = HasBranch(tEvent, "pclusterCompatibilityFilter");

  tEvent->SetBranchStatus("*", 0);

  if (has_pprimaryVertexFilter && has_pclusterCompatibilityFilter)
  {
    tEvent->SetBranchStatus("pprimaryVertexFilter", 1);
    tEvent->SetBranchStatus("pclusterCompatibilityFilter", 1);

    tEvent->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter);
    tEvent->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);
  }

  tHi->SetBranchStatus("*", 0);

  if (HasBranch(tHi, "vz"))
  {
    tHi->SetBranchStatus("vz", 1);
    tHi->SetBranchAddress("vz", &vz);
  }
  if (HasBranch(tHi, "hiBin"))
  {
    tHi->SetBranchStatus("hiBin", 1);
    tHi->SetBranchAddress("hiBin", &hiBin);
  }

  // -------------------------
  // PF tree
  // -------------------------
  Int_t nPF = 0;
  std::vector<int> *pfId = nullptr;
  std::vector<float> *pfPt = nullptr;
  std::vector<float> *pfPhi = nullptr;

  tPF->SetBranchStatus("*", 0);
  tPF->SetBranchStatus("nPF", 1);
  tPF->SetBranchStatus("pfId", 1);
  tPF->SetBranchStatus("pfPt", 1);
  tPF->SetBranchStatus("pfPhi", 1);

  tPF->SetBranchAddress("nPF", &nPF);
  tPF->SetBranchAddress("pfId", &pfId);
  tPF->SetBranchAddress("pfPt", &pfPt);
  tPF->SetBranchAddress("pfPhi", &pfPhi);

  // -------------------------
  // Gen Branches (Just for MC)
  // -------------------------

  std::vector<float> *genPt = nullptr;
  std::vector<float> *genEta = nullptr;
  std::vector<float> *genPhi = nullptr;
  std::vector<int> *genChg = nullptr;

  std::vector<int> *genPdg = nullptr;
  std::vector<int> *genNMothers = nullptr;
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

    tGen->SetBranchAddress("pt", &genPt);
    tGen->SetBranchAddress("eta", &genEta);
    tGen->SetBranchAddress("phi", &genPhi);
    tGen->SetBranchAddress("chg", &genChg);

    tGen->SetBranchAddress("pdg", &genPdg);
    tGen->SetBranchAddress("nMothers", &genNMothers);
    tGen->SetBranchAddress("motherIdx", &genMotherIdx);
  }

  // -------------------------
  // Output hist (after final step 8)
  // -------------------------
  static const int NY = 12;
  double yEdges[NY + 1] = {
      -2.4, -2.0, -1.6, -1.2, -0.8, -0.4,
      0.0, 0.4, 0.8, 1.2, 1.6, 2.0,
      2.4};
  double yEdges_FB[NY + 1] = {
      -1.7068,
      -1.3646,
      -1.0223,
      -0.6801,
      -0.3379,
      0.0044,
      0.3466,
      0.6888,
      1.0311,
      1.3733,
      1.7155,
      2.0578,
      2.4000};

  TH1D *h_met_Wp[NY];
  TH1D *h_met_Wm[NY];
  TH1D *h_mt_Wp[NY];
  TH1D *h_mt_Wm[NY];

  TH1D *h_met_Wp_FB[NY];
  TH1D *h_met_Wm_FB[NY];
  TH1D *h_mt_Wp_FB[NY];
  TH1D *h_mt_Wm_FB[NY];

  for (int b = 0; b < NY; ++b)
  {
    const double y1 = yEdges[b];
    const double y2 = yEdges[b + 1];
    const double y1_FB = yEdges_FB[b];
    const double y2_FB = yEdges_FB[b + 1];

    h_met_Wp[b] = new TH1D(Form("h_met_Wp_y%d", b),
                           Form("W+ PF MET;MET [GeV];Events (%.2f<y<%.2f)", y1, y2),
                           60, 0, 120);
    h_met_Wm[b] = new TH1D(Form("h_met_Wm_y%d", b),
                           Form("W- PF MET;MET [GeV];Events (%.2f<y<%.2f)", y1, y2),
                           60, 0, 120);

    h_mt_Wp[b] = new TH1D(Form("h_mt_Wp_y%d", b),
                          Form("W+ m_{T};m_{T} [GeV];Events (%.2f<y<%.2f)", y1, y2),
                          80, 0, 200);
    h_mt_Wm[b] = new TH1D(Form("h_mt_Wm_y%d", b),
                          Form("W- m_{T};m_{T} [GeV];Events (%.2f<y<%.2f)", y1, y2),
                          80, 0, 200);

    h_met_Wp_FB[b] = new TH1D(Form("h_met_Wp_y%d_FB", b),
                              Form("W+ PF MET;MET [GeV];Events (%.2f<y<%.2f)", y1_FB, y2_FB),
                              60, 0, 120);
    h_met_Wm_FB[b] = new TH1D(Form("h_met_Wm_y%d_FB", b),
                              Form("W- PF MET;MET [GeV];Events (%.2f<y<%.2f)", y1_FB, y2_FB),
                              60, 0, 120);

    h_mt_Wp_FB[b] = new TH1D(Form("h_mt_Wp_y%d_FB", b),
                             Form("W+ m_{T};m_{T} [GeV];Events (%.2f<y<%.2f)", y1_FB, y2_FB),
                             80, 0, 200);
    h_mt_Wm_FB[b] = new TH1D(Form("h_mt_Wm_y%d_FB", b),
                             Form("W- m_{T};m_{T} [GeV];Events (%.2f<y<%.2f)", y1_FB, y2_FB),
                             80, 0, 200);

    // optional but recommended if you later sum/scale:
    h_met_Wp[b]->Sumw2();
    h_met_Wm[b]->Sumw2();
    h_mt_Wp[b]->Sumw2();
    h_mt_Wm[b]->Sumw2();

    h_met_Wp_FB[b]->Sumw2();
    h_met_Wm_FB[b]->Sumw2();
    h_mt_Wp_FB[b]->Sumw2();
    h_mt_Wm_FB[b]->Sumw2();
  }

  // QCD BK histogram

  static const int NISO = 3;
  double isoEdges[NISO + 1] = {0.1, 0.4, 0.7, 1.0};

  TH1D *h_met_iso_muPlus[NISO];
  TH1D *h_met_iso_muMinus[NISO];

  for (int b = 0; b < NISO; ++b)
  {
    h_met_iso_muPlus[b] = new TH1D(
        Form("h_met_iso_muPlus_bin%d", b),
        Form("MET, #mu^{+}, %.1f < iso < %.1f;MET [GeV];Events", isoEdges[b], isoEdges[b + 1]),
        60, 0, 120);

    h_met_iso_muMinus[b] = new TH1D(
        Form("h_met_iso_muMinus_bin%d", b),
        Form("MET, #mu^{-}, %.1f < iso < %.1f;MET [GeV];Events", isoEdges[b], isoEdges[b + 1]),
        60, 0, 120);

    h_met_iso_muPlus[b]->Sumw2();
    h_met_iso_muMinus[b]->Sumw2();
  }

  // -------------------------
  // Cutflow counters
  // -------------------------
  // Ni in the table = events passing current and all previous cuts
  unsigned long long N[9] = {0};

  const Long64_t nEntries = tMu->GetEntries();
  std::cout << "Entries: " << nEntries << "\n";

  bool warnedEventFiltersOnce = false;
  bool warnedNoTrigMatchInfo = false;

  // Iso bin finder:

  auto FindIsoBin = [&](double iso) -> int
  {
    for (int b = 0; b < NISO; ++b)
    {
      if (iso >= isoEdges[b] && iso < isoEdges[b + 1])
        return b;
    }
    return -1;
  };

  for (Long64_t ie = 0; ie < nEntries; ++ie)
  {
    if (ie % 200000 == 0)
      std::cout << "Event " << ie << "/" << nEntries << "\n";

    tMu->GetEntry(ie);
    tHi->GetEntry(ie);
    tPF->GetEntry(ie);
    tHLT->GetEntry(ie);
    tHLTobj->GetEntry(ie);
    tEvent->GetEntry(ie);
    if (isMC)
      tGen->GetEntry(ie);

    // (0) all events
    N[0]++;

    // (1) >=1 PF muon with pT > 25
    if (!ExistsPFMuonPt25(cfg, nMu, muPt, has_muIsPF, muIsPF))
      continue;
    N[1]++;

    // (2) pO collision event selection
    if (!PassEventSelection_pO(tHi,
                               warnedEventFiltersOnce,
                               has_pprimaryVertexFilter, pprimaryVertexFilter,
                               has_pclusterCompatibilityFilter, pclusterCompatibilityFilter))
      continue;
    if (cfg.applyVz)
    {
      if (TMath::Abs(vz) > cfg.vzMax)
        continue;
    }
    N[2]++;

    // (3) Trigger PAL3Mu12 fired
    if (has_hlt)
    {
      if (!TriggerFired(HLT_OxyL1SingleMuOpen_v1))
        continue;
    }
    // if has_hlt==false, treat as pass (warned already)
    N[3]++;

    // (4) Drell–Yan veto
    if (!PassDYVeto(cfg,
                    nMu, muPt, muEta, muPhi, muCharge,
                    has_muIDTight, muIDTight,
                    has_muIsPF, muIsPF,
                    muPFChIso, muPFNeuIso, muPFPhoIso))
      continue;
    N[4]++;

    // (5) >=1 Tight ID muon
    if (!ExistsTightMuon(nMu, has_muIDTight, muIDTight))
      continue;
    N[5]++;

    // Define leading muon for steps 6-8 (leading among Tight+PF)
    const int iLead = FindLeadingMuon_TightPF(nMu, muPt, muEta, muPhi, has_muIDTight, muIDTight, has_muIsPF, muIsPF);
    if (iLead < 0)
      continue;

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

    // (6) Leading muon pT > 25 and |eta| < 2.4
    if (!muPt || muPt->at(iLead) <= cfg.muPt25)
      continue;
    if (!muEta || abs(muEta->at(iLead)) > cfg.muEtaMax)
      continue;
    N[6]++;

    // (7) Leading muon Iso < 0.15
    const double isoLead = RelIsoPF(iLead, muPt, muPFChIso, muPFNeuIso, muPFPhoIso);

    const int isoBin = FindIsoBin(isoLead);
    const bool isQCDSideband = (isoBin >= 0);

    // nominal signal-region bool
    const bool passIsoNominal = (isoLead < cfg.isoMax);

    if (passIsoNominal)
      N[7]++;

    // (8) Leading muon matched to trigger
    const bool passMatch = PassLeadingMuonTrigMatch(cfg, iLead, muEta, muPhi,
                                                    has_trgObjPt, trgObjPt,
                                                    has_trgObjEta, trgObjEta,
                                                    has_trgObjPhi, trgObjPhi,
                                                    has_trgObjId, trgObjId,
                                                    warnedNoTrigMatchInfo);

    if (!passMatch)
      continue;

    if (passIsoNominal)
      N[8]++;

    // Fill final distributions (after step 8)
    TVector2 metv = ComputePFMET(pfId, pfPt, pfPhi);
    const double met = metv.Mod();
    const double metPhi = metv.Phi();
    const double mu_phi = muPhi->at(iLead);
    const double dphi = TVector2::Phi_mpi_pi(mu_phi - metPhi);
    const double mt = std::sqrt(2.0 * muPt->at(iLead) * met * (1.0 - std::cos(dphi)));

    if (isQCDSideband)
    {
      if (muCharge->at(iLead) > 0)
        h_met_iso_muPlus[isoBin]->Fill(met);
      else if (muCharge->at(iLead) < 0)
        h_met_iso_muMinus[isoBin]->Fill(met);
    }

    if (!passIsoNominal)
      continue;

    auto FindBin = [&](double y) -> int
    {
      for (int b = 0; b < NY; b++)
      {
        if (y >= yEdges[b] && y < yEdges[b + 1])
          return b;
      }
      if (y == yEdges[NY])
        return NY - 1; // include right edge
      return -1;
    };

    auto FindBin_FB = [&](double y) -> int
    {
      for (int b = 0; b < NY; b++)
      {
        if (y >= yEdges_FB[b] && y < yEdges_FB[b + 1])
          return b;
      }
      if (y == yEdges_FB[NY])
        return NY - 1; // include right edge
      return -1;
    };

    int q = muCharge->at(iLead); // +1 or -1
    bool isWp = (q > 0);
    bool isWm = (q < 0);

    // [NOTICE] I filp the sign here for eta, define p going (-Z) as forward

    const double y = (-1) * muEta->at(iLead);
    const int ybin = FindBin(y);
    const int ybin_FB = FindBin_FB(y);

    if (ybin >= 0)
    {
      if (isWp)
      {
        h_met_Wp[ybin]->Fill(met);
        h_mt_Wp[ybin]->Fill(mt);
      }
      else if (isWm)
      {
        h_met_Wm[ybin]->Fill(met);
        h_mt_Wm[ybin]->Fill(mt);
      }
    }

    if (ybin_FB >= 0)
    {
      if (isWp)
      {
        h_met_Wp_FB[ybin_FB]->Fill(met);
        h_mt_Wp_FB[ybin_FB]->Fill(mt);
      }
      else if (isWm)
      {
        h_met_Wm_FB[ybin_FB]->Fill(met);
        h_mt_Wm_FB[ybin_FB]->Fill(mt);
      }
    }
  }

  // -------------------------
  // Print cutflow like the table (Ni and ratios)
  // -------------------------
  std::cout << "\n========== Cutflow (AN table logic) ==========\n";
  std::cout << " i   Ni              (Ni-Ni-1)/N1        Ni/Ni-1\n";
  std::cout << "--------------------------------------------------\n";

  const double N1 = (N[1] > 0 ? (double)N[1] : 1.0);

  for (int i = 0; i <= 8; ++i)
  {
    double fracStep = 0.0;
    double ratio = 0.0;

    if (i >= 1 && N[i - 1] > 0)
      ratio = (double)N[i] / (double)N[i - 1];

    if (i >= 2)
      fracStep = ((double)N[i] - (double)N[i - 1]) / N1;

    // format similar to table:
    // - for i=0/1: table shows Ni/Ni-1 only for i=1; step frac shown from i>=2
    std::cout << " " << i << "   " << N[i];

    if (i >= 2)
      std::cout << "        " << Form("%+7.2f%%", 100.0 * fracStep);
    else
      std::cout << "        " << "   (n/a) ";

    if (i >= 1)
      std::cout << "        " << Form("%7.2f%%", 100.0 * ratio);
    else
      std::cout << "        " << "  (n/a) ";

    std::cout << "\n";
  }

  std::cout << "==================================================\n\n";

  // ---------------------------------------------
  // Save cutflow table to text file (overwrite)
  // ---------------------------------------------

  std::string outDir = "./output/";
  std::string mcTag = isMC ? "MC" : "Data";

  std::string outFile =
      outDir + cfg.outPrefix + "_" + mcTag + ".txt";

  std::ofstream ftxtout(outFile.c_str(), std::ios::out | std::ios::trunc);
  if (!ftxtout.is_open())
  {
    std::cerr << "ERROR: cannot open " << outFile << std::endl;
    return;
  }
  else
  {
    ftxtout << "\n========== Cutflow (AN table logic) ==========\n";
    ftxtout << " i   Ni              (Ni-Ni-1)/N1        Ni/Ni-1\n";
    ftxtout << "--------------------------------------------------\n";

    const double N1 = (N[1] > 0 ? (double)N[1] : 1.0);

    for (int i = 0; i <= 8; ++i)
    {
      double fracStep = 0.0;
      double ratio = 0.0;

      if (i >= 1 && N[i - 1] > 0)
        ratio = (double)N[i] / (double)N[i - 1];

      if (i >= 2)
        fracStep = ((double)N[i] - (double)N[i - 1]) / N1;

      ftxtout << " " << i << "   " << N[i];

      if (i >= 2)
        ftxtout << "        " << Form("%+7.2f%%", 100.0 * fracStep);
      else
        ftxtout << "        " << "   (n/a) ";

      if (i >= 1)
        ftxtout << "        " << Form("%7.2f%%", 100.0 * ratio);
      else
        ftxtout << "        " << "  (n/a) ";

      ftxtout << "\n";
    }

    ftxtout << "==================================================\n\n";
    ftxtout.close();
  }

  // Save hists
  TFile *fout = new TFile(("./rootfile/" + cfg.outPrefix + "_hist.root").c_str(), "RECREATE");

  for (int i = 0; i < NY; i++)
  {
    h_met_Wp[i]->Write("", 2);
    h_met_Wm[i]->Write("", 2);
    h_mt_Wp[i]->Write("", 2);
    h_mt_Wm[i]->Write("", 2);

    h_met_Wp_FB[i]->Write("", 2);
    h_met_Wm_FB[i]->Write("", 2);
    h_mt_Wp_FB[i]->Write("", 2);
    h_mt_Wm_FB[i]->Write("", 2);
  }

  for (int b = 0; b < NISO; ++b)
  {
    h_met_iso_muPlus[b]->Write("", TObject::kOverwrite);
    h_met_iso_muMinus[b]->Write("", TObject::kOverwrite);
  }
  fout->Close();

  std::cout << "[INFO] Wrote outputs with prefix: " << cfg.outPrefix << "\n";
}