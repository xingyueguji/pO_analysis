// DrawWToEleNu_PFMet_CutflowTableStrict.C
//
// Goal: follow the cutflow logic EXACTLY as in the table:
//
//  0) All events in dataset
//  1) >=1 PF electron with pT > 25
//  2) pO collision event selection
//  3) Trigger PAL3Ele12 fired
//  4) Drell–Yan veto
//  5) >=1 Tight ID electron
//  6) Leading electron pT > 25
//  7) Leading electron Iso < 0.15
//  8) Leading electron matched to trigger
//
// Notes / assumptions (kept as “auto-detect” wherever possible):
// - “PF electron” is taken as eleIsPF == 1 when the branch exists.
// - “Tight ID electron” is taken as eleIDTight == 1 when the branch exists.
// - Event selection (step 2): we AND together common pO/pA filters IF they exist;
//   if a given filter branch is missing, we don’t apply it (and we print a warning once).
// - Trigger fired (step 3): we auto-find any branch name containing "HLT_PAL3Ele12"
//   in ggHiNtuplizer/EventTree and use it.
// - Trigger matching (step 8): we try, in order:
//    (A) trigger object collections (trgObj* style) if present
//    (B) per-electron trigger match bits (eleTrig/eleTrg style) if present
//   If neither exists, we can’t enforce step 8 and we will WARN and treat as pass.
// - DY veto (step 4): veto event if exists an OS electron pair with BOTH electrons:
//      pT>15, PF, Tight, relIso<0.15, and (mll > 30 GeV)   <-- your requested m>30
//
// Run:
//   root -l -q 'DrawWToEleNu_PFMet_CutflowTableStrict.C("root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2025.root")'
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
  kData = 0,
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
  double elePt25 = 25.0;
  double dyElePtMin = 10.0;
  double isoMax = 0.095;
  double dyMassMin = 80.0; // your requested M>30 inside DY veto
  double dyMassMax = 110.0;
  double eleEtaMax = 2.4;

  // For trigger-object matching
  double trigMatchDR = 0.1;
  bool applyVz = true;
  double vzMax = 15.0;

  // Histograms
  int nBinsMt = 100;
  double mtMin = 0.0;
  double mtMax = 200.0;

  std::string outPrefix = "WToElecNu_pO_PFMet";
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
                       std::vector<float> *elePt,
                       std::vector<float> *elePFChIso,
                       std::vector<float> *elePFNeuIso,
                       std::vector<float> *elePFPhoIso)
{
  if (!elePt || !elePFChIso || !elePFNeuIso || !elePFPhoIso)
    return 999.0;
  const double pt = elePt->at(i);
  if (pt <= 0)
    return 999.0;
  const double isoAbs = elePFChIso->at(i) + elePFNeuIso->at(i) + elePFPhoIso->at(i);
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
// Step 1 helper: exists PF electron with pt>25 // to be changed
// ---------------------------------------------
static bool ExistsPFElectronPt25(const WConfig &cfg,
                                 int nEle,
                                 std::vector<float> *elePt)
{
  if (!elePt)
    return false;

  for (int i = 0; i < nEle; ++i)
  {
    const double pt = elePt->at(i);
    if (pt <= cfg.elePt25)
      continue;

    return true;
  }
  return false;
}

// ---------------------------------------------
// Step 5 helper: exists Tight ID electron (no pt requirement in table) // to be changed
// ---------------------------------------------
static bool ExistsTightElectron(int nEle,
                                std::vector<int> *eleIDTight)
{

  for (int i = 0; i < nEle; ++i)
  {
    if (eleIDTight->at(i) != 0)
      return true;
  }
  return false;
}

// ---------------------------------------------
// Leading electron definition for steps 6-8:
// We'll define "leading electron" among electrons that are Tight (if exists) and PF (if exists).
// This is the closest to typical W analysis and makes steps 6-8 consistent.
// ---------------------------------------------
static int FindLeadingElectron_TightPF(int nEle,
                                       std::vector<float> *elePt,
                                       std::vector<float> *eleEta,
                                       std::vector<float> *elePhi,
                                       std::vector<int> *eleIDTight)
{
  if (!elePt || !eleEta || !elePhi)
    return -1;

  int best = -1;
  double bestPt = -1;

  for (int i = 0; i < nEle; ++i)
  {
    const double pt = elePt->at(i);

    // Tight requirement if branch exists; if missing, cannot define “tight leading electron”
    if (eleIDTight && eleIDTight->at(i) == 0)
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
// Veto event if exists OS pair where both electrons satisfy:
//   pT > 15, PF, Tight, relIso < 0.15, AND mll > 30
// ---------------------------------------------
static bool PassDYVeto(const WConfig &cfg,
                       int nEle,
                       std::vector<float> *elePt,
                       std::vector<float> *eleEta,
                       std::vector<float> *elePhi,
                       std::vector<int> *eleCharge,
                       std::vector<int> *eleIDTight,
                       std::vector<int> *eleIso)
{
  if (!elePt || !eleEta || !elePhi || !eleCharge)
    return true; // can't apply -> do not veto

  std::vector<int> cand;
  cand.reserve(nEle);

  for (int i = 0; i < nEle; ++i)
  {
    const double pt = elePt->at(i);
    if (pt <= cfg.dyElePtMin)
      continue;

    if (eleIDTight->at(i) == 0)
      continue;

    if (eleIso->at(i) == 0)
      continue;

    cand.push_back(i);
  }

  const double eleMass = 0.000511;

  for (size_t a = 0; a < cand.size(); ++a)
  {
    for (size_t b = a + 1; b < cand.size(); ++b)
    {
      const int i1 = cand[a];
      const int i2 = cand[b];

      if (eleCharge->at(i1) * eleCharge->at(i2) >= 0)
        continue; // require OS

      TLorentzVector ele1, ele2;
      ele1.SetPtEtaPhiM(elePt->at(i1), eleEta->at(i1), elePhi->at(i1), eleMass);
      ele2.SetPtEtaPhiM(elePt->at(i2), eleEta->at(i2), elePhi->at(i2), eleMass);

      const double mll = (ele1 + ele2).M();

      // DY-like OS pair with mll > 30 => veto event (fail DY veto)
      if (mll > cfg.dyMassMin && mll < cfg.dyMassMax)
        return false;
    }
  }

  return true;
}

// ---------------------------------------------
// Step 3: Trigger fired (PAL3Ele12) from HLT bit branch (auto-found)
// ---------------------------------------------
static bool TriggerFired(int hltBit)
{
  return (hltBit != 0);
}

// ---------------------------------------------
// Step 8: Leading electron matched to trigger
// Try trigger objects, else try per-electron match bits.
// If neither exists, warn and accept (cannot enforce).
// ---------------------------------------------
static bool PassLeadingElectronTrigMatch(const WConfig &cfg,
                                         int iLead,
                                         std::vector<float> *eleEta,
                                         std::vector<float> *elePhi,

                                         // Option A: trigger object branches (if exist)
                                         bool has_trgObjPt, std::vector<double> *trgObjPt,
                                         bool has_trgObjEta, std::vector<double> *trgObjEta,
                                         bool has_trgObjPhi, std::vector<double> *trgObjPhi,
                                         bool has_trgObjId, std::vector<double> *trgObjId,
                                         bool &warnedNoTrigMatchInfo)
{
  if (iLead < 0 || !eleEta || !elePhi)
    return false;

  const double eta = eleEta->at(iLead);
  const double phi = elePhi->at(iLead);

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
    std::cout << "[WARN] No trigger-object collection and no per-electron trigger-match info found. "
              << "Cannot apply step 8 (leading electron matched to trigger). Treating as PASS.\n";
    warnedNoTrigMatchInfo = true;
  }
  return true;
}

// ---------------------------------------------
// Extra: Gen-reco matching
// ---------------------------------------------

static bool PassGenRecoMatching(
    int iLead,
    std::vector<float> *elePt, std::vector<float> *eleEta, std::vector<float> *elePhi, std::vector<int> *eleCharge,
    std::vector<float> *genPt, std::vector<float> *genEta, std::vector<float> *genPhi, std::vector<int> *genChg)
{
  if (!elePt || !eleEta || !elePhi || !eleCharge || !genPt || !genEta || !genPhi || !genChg)
    return false;
  if (iLead < 0 || (size_t)iLead >= elePt->size())
    return false;

  const double electronPt = elePt->at(iLead);
  const double electronEta = eleEta->at(iLead);
  const double electronPhi = elePhi->at(iLead);
  const int electronChg = eleCharge->at(iLead);

  const size_t ngen = std::min({genPt->size(), genEta->size(), genPhi->size(), genChg->size()});
  for (size_t i = 0; i < ngen; ++i)
  {
    // avoid division by 0 in dPtRel
    if (genPt->at(i) <= 0)
      continue;

    const double dEta = genEta->at(i) - electronEta;
    const double dPhi = TVector2::Phi_mpi_pi(genPhi->at(i) - electronPhi);
    const double dR = std::sqrt(dEta * dEta + dPhi * dPhi);

    const double dPtRel = std::fabs(genPt->at(i) - electronPt) / genPt->at(i);

    if (genChg->at(i) != electronChg)
      continue; // same charge
    if (dR >= 0.5)
      continue; // ΔR < 0.5
    if (dPtRel >= 0.5)
      continue; // δpT < 0.5

    return true; // found a match satisfying all cuts
  }

  return false;
}

// ---------------------------------------------
// Main
// ---------------------------------------------
void DrawWToElecNu_PFMet(const char *fname =
                             "root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2025.root",
                         SampleType sample = kData)
{
  bool isMC = false;

  if (sample == kData)
  {
    isMC = false;
  }
  else
  {
    isMC = true;
  }

  gStyle->SetOptStat(0);

  WConfig cfg;

  if (isMC)
  {
    if (sample == kDY)
    {
      fname = Form("root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_MC_DY_ele_Z.root");
      cfg.outPrefix += "_DY";
    }
    if (sample == kWp)
    {
      fname = Form("root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_MC_Wp_ele.root");
      cfg.outPrefix += "_Wp";
    }
    if (sample == kWm)
    {
      fname = Form("root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_MC_Wm_ele.root");
      cfg.outPrefix += "_Wm";
    }
  }

  TFile *f = TFile::Open(fname);
  if (!f || f->IsZombie())
  {
    std::cerr << "ERROR opening file\n";
    return;
  }

  TTree *tEle = (TTree *)f->Get("ggHiNtuplizer/EventTree");
  TTree *tHi = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  TTree *tPF = (TTree *)f->Get("particleFlowAnalyser/pftree");
  TTree *tHLT = (TTree *)f->Get("hltanalysis/HltTree");
  TTree *tHLTobj = (TTree *)f->Get("hltobject/HLT_OxyL1SingleEG10_v");
  TTree *tEvent = (TTree *)f->Get("skimanalysis/HltTree");
  TTree *tGen;
  if (isMC)
  {
    tGen = (TTree *)f->Get("HiGenParticleAna/hi");
  }

  if (!tEle || !tHi || !tPF || !tHLT || !tHLTobj)
  {
    std::cerr << "ERROR: missing one of trees: ggHiNtuplizer/EventTree, hiEvtAnalyzer/HiTree, particleFlowAnalyser/pftree\n";
    return;
  }

  // -------------------------
  // Electron branches
  // -------------------------
  Int_t nEle = 0;
  std::vector<float> *elePt = nullptr;
  std::vector<float> *eleEta = nullptr;
  std::vector<float> *elePhi = nullptr;
  std::vector<int> *eleCharge = nullptr;

  std::vector<float> *elePFChIso = nullptr;
  std::vector<float> *elePFNeuIso = nullptr;
  std::vector<float> *elePFPhoIso = nullptr;
  std::vector<float> *elePFPUIso = nullptr;

  // Optional per-electron trigger match
  std::vector<int> *eleTrig = nullptr;

  tEle->SetBranchStatus("*", 0);
  tEle->SetBranchStatus("nEle", 1);
  tEle->SetBranchStatus("elePt", 1);
  tEle->SetBranchStatus("eleEta", 1);
  tEle->SetBranchStatus("elePhi", 1);
  tEle->SetBranchStatus("eleCharge", 1);

  tEle->SetBranchAddress("nEle", &nEle);
  tEle->SetBranchAddress("elePt", &elePt);
  tEle->SetBranchAddress("eleEta", &eleEta);
  tEle->SetBranchAddress("elePhi", &elePhi);
  tEle->SetBranchAddress("eleCharge", &eleCharge);

  // Here's the electron CUT BASED ID and MVA BASED ID

  std::vector<int> *eleMVAIdWP80 = nullptr;
  std::vector<int> *eleMVAIdWP85 = nullptr;
  std::vector<int> *eleMVAIdWP90 = nullptr;
  std::vector<int> *eleMVAIdWP95 = nullptr;

  std::vector<int> *eleCutIdWP70 = nullptr;
  std::vector<int> *eleCutIdWP80 = nullptr;
  std::vector<int> *eleCutIdWP90 = nullptr;
  std::vector<int> *eleCutIdWP95 = nullptr;

  tEle->SetBranchStatus("eleMVAIdWP80", 1);
  tEle->SetBranchStatus("eleMVAIdWP85", 1);
  tEle->SetBranchStatus("eleMVAIdWP90", 1);
  tEle->SetBranchStatus("eleMVAIdWP95", 1);
  tEle->SetBranchStatus("eleCutIdWP70", 1);
  tEle->SetBranchStatus("eleCutIdWP80", 1);
  tEle->SetBranchStatus("eleCutIdWP90", 1);
  tEle->SetBranchStatus("eleCutIdWP95", 1);

  tEle->SetBranchAddress("eleMVAIdWP80", &eleMVAIdWP80);
  tEle->SetBranchAddress("eleMVAIdWP85", &eleMVAIdWP85);
  tEle->SetBranchAddress("eleMVAIdWP90", &eleMVAIdWP90);
  tEle->SetBranchAddress("eleMVAIdWP95", &eleMVAIdWP95);
  tEle->SetBranchAddress("eleCutIdWP70", &eleCutIdWP70);
  tEle->SetBranchAddress("eleCutIdWP80", &eleCutIdWP80);
  tEle->SetBranchAddress("eleCutIdWP90", &eleCutIdWP90);
  tEle->SetBranchAddress("eleCutIdWP95", &eleCutIdWP95);

  // Iso
  if (HasBranch(tEle, "elePFChIso"))
    tEle->SetBranchStatus("elePFChIso", 1);
  if (HasBranch(tEle, "elePFNeuIso"))
    tEle->SetBranchStatus("elePFNeuIso", 1);
  if (HasBranch(tEle, "elePFPhoIso"))
    tEle->SetBranchStatus("elePFPhoIso", 1);

  const bool has_elePFChIso = HasBranch(tEle, "elePFChIso");
  const bool has_elePFNeuIso = HasBranch(tEle, "elePFNeuIso");
  const bool has_elePFPhoIso = HasBranch(tEle, "elePFPhoIso");

  if (has_elePFChIso)
    tEle->SetBranchAddress("elePFChIso", &elePFChIso);
  if (has_elePFNeuIso)
    tEle->SetBranchAddress("elePFNeuIso", &elePFNeuIso);
  if (has_elePFPhoIso)
    tEle->SetBranchAddress("elePFPhoIso", &elePFPhoIso);

  // Another region for isolation
  std::vector<int> *eleMVAIsoWP80 = nullptr;
  std::vector<int> *eleMVAIsoWP85 = nullptr;
  std::vector<int> *eleMVAIsoWP90 = nullptr;
  std::vector<int> *eleMVAIsoWP95 = nullptr;

  tEle->SetBranchStatus("eleMVAIsoWP80", 1);
  tEle->SetBranchStatus("eleMVAIsoWP85", 1);
  tEle->SetBranchStatus("eleMVAIsoWP90", 1);
  tEle->SetBranchStatus("eleMVAIsoWP95", 1);

  tEle->SetBranchAddress("eleMVAIsoWP80", &eleMVAIsoWP80);
  tEle->SetBranchAddress("eleMVAIsoWP85", &eleMVAIsoWP85);
  tEle->SetBranchAddress("eleMVAIsoWP90", &eleMVAIsoWP90);
  tEle->SetBranchAddress("eleMVAIsoWP95", &eleMVAIsoWP95);

  // -------------------------
  // Trigger fired bit (step 3):
  // -------------------------
  const std::string hltName = FindBranchContaining(tHLT, "HLT_OxyL1SingleEG10_v1");
  Int_t HLT_OxyL1SingleEG10_v1 = 0;
  bool has_hlt = false;
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

  if (isMC && !tGen)
  {
    std::cerr << "[FATAL] isMC=true but no gen tree found\n";
    return;
  }

  if (isMC)
  {
    tGen->SetBranchStatus("*", 0);
    tGen->SetBranchStatus("pt", 1);
    tGen->SetBranchStatus("eta", 1);
    tGen->SetBranchStatus("phi", 1);
    tGen->SetBranchStatus("chg", 1);

    tGen->SetBranchAddress("pt", &genPt);
    tGen->SetBranchAddress("eta", &genEta);
    tGen->SetBranchAddress("phi", &genPhi);
    tGen->SetBranchAddress("chg", &genChg);
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

  // -------------------------
  // Cutflow counters
  // -------------------------
  // Ni in the table = events passing current and all previous cuts
  unsigned long long N[9] = {0};

  const Long64_t nEntries = tEle->GetEntries();
  std::cout << "Entries: " << nEntries << "\n";

  bool warnedEventFiltersOnce = false;
  bool warnedNoTrigMatchInfo = false;

  for (Long64_t ie = 0; ie < nEntries; ++ie)
  {
    if (ie % 200000 == 0)
      std::cout << "Event " << ie << "/" << nEntries << "\n";

    tEle->GetEntry(ie);
    tHi->GetEntry(ie);
    tPF->GetEntry(ie);
    tHLT->GetEntry(ie);
    tHLTobj->GetEntry(ie);
    tEvent->GetEntry(ie);

    // (0) all events
    N[0]++;

    // (1) >=1 PF eletron with pT > 25
    if (!ExistsPFElectronPt25(cfg, nEle, elePt))
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

    // (3) Trigger PAL3Ele12 fired
    if (has_hlt)
    {
      if (!TriggerFired(HLT_OxyL1SingleEG10_v1))
        continue;
    }
    // if has_hlt==false, treat as pass (warned already)
    N[3]++;

    // (4) Drell–Yan veto
    if (!PassDYVeto(cfg,
                    nEle, elePt, eleEta, elePhi, eleCharge,
                    eleCutIdWP95, eleMVAIsoWP95))
      continue;
    N[4]++;

    // (5) >=1 Tight ID electron
    if (!ExistsTightElectron(nEle, eleCutIdWP95))
      continue;
    N[5]++;

    // Define leading electron for steps 6-8 (leading among Tight+PF)
    const int iLead = FindLeadingElectron_TightPF(nEle, elePt, eleEta, elePhi, eleCutIdWP95);
    if (iLead < 0)
      continue;

    if (isMC)
    {
      if (!PassGenRecoMatching(iLead, elePt, eleEta, elePhi, eleCharge, genPt, genEta, genPhi, genChg))
      {
        // cout << "We see a gen-reco matching failed case" << endl;
        // continue;
      }
    }

    // (6) Leading electron pT > 25 and |eta| < 2.4
    if (!elePt || elePt->at(iLead) <= cfg.elePt25)
      continue;
    if (!eleEta || abs(eleEta->at(iLead)) > cfg.eleEtaMax)
      continue;
    N[6]++;

    // (7) Leading electron Iso < 0.15
    if (eleMVAIsoWP95->at(iLead) != 1)
      continue;

    // const double isoLead = RelIsoPF(iLead, elePt, elePFChIso, elePFNeuIso, elePFPhoIso);
    // if (isoLead >= cfg.isoMax)
    // continue;
    N[7]++;

    // (8) Leading electron matched to trigger
    const bool passMatch = PassLeadingElectronTrigMatch(cfg, iLead, eleEta, elePhi,
                                                        has_trgObjPt, trgObjPt,
                                                        has_trgObjEta, trgObjEta,
                                                        has_trgObjPhi, trgObjPhi,
                                                        has_trgObjId, trgObjId,
                                                        warnedNoTrigMatchInfo);
    if (!passMatch)
      continue;
    N[8]++;

    // Fill final distributions (after step 8)
    TVector2 metv = ComputePFMET(pfId, pfPt, pfPhi);
    const double met = metv.Mod();
    const double metPhi = metv.Phi();
    const double ele_phi = elePhi->at(iLead);
    const double dphi = TVector2::Phi_mpi_pi(ele_phi - metPhi);
    const double mt = std::sqrt(2.0 * elePt->at(iLead) * met * (1.0 - std::cos(dphi)));

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

    int q = eleCharge->at(iLead); // +1 or -1
    bool isWp = (q > 0);
    bool isWm = (q < 0);

    // [NOTICE] I filp the sign here for eta, define p going (-Z) as forward

    const double y = (-1) * eleEta->at(iLead);
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
  fout->Close();

  std::cout << "[INFO] Wrote outputs with prefix: " << cfg.outPrefix << "\n";
}