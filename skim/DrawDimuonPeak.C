// DrawDimuonPeak.C
// Example:
//   root -l -q 'DrawDimuonPeak.C("root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2025.root")'

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TMath.h"
#include "TPaveText.h"
#include <vector>
#include <iostream>
#include <string>

// --------------------
// Small helpers
// --------------------
static double MU_MASS() { return 0.1056583745; }

static double ELE_MASS() { return 0.000511; }

static bool HasBranch(TTree *t, const char *name)
{
  return (t && t->GetListOfBranches() && t->GetListOfBranches()->FindObject(name));
}

static bool TriggerFired(int hltBit)
{
  return (hltBit != 0);
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

// --------------------
// User config struct
// --------------------
struct DileptonConfig
{
  // resonance window (plot range)
  double massMin = 60.0;
  double massMax = 120.0;
  int nBins = 120;

  // selection
  double ptMin1 = 20.0; // leading
  double ptMin2 = 10.0; // subleading
  double etaMax = 2.4;

  bool requireOS = true;

  bool requireGood = false;   // muIsGood
  bool requireGlobal = true;  // muIsGlobal
  bool requirePF = true;      // muIsPF
  bool requireTightID = true; // muIDTight

  // optional event cuts (only applied if branch exists and enabled)
  bool applyVz = true;
  double vzMax = 15.0;

  bool applyHiBin = false; // OFF by default for pO
  int hiBinMin = 0;
  int hiBinMax = 200;

  double isoMax = 0.2;

  // output
  std::string outPrefix = "dimuon";
};

// --------------------
// Main
// --------------------
void DrawDimuonPeak(const char *fname = "root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2025.root")
{
  // ==============
  // Configure here
  // ==============
  DileptonConfig cfg;

  bool warnedEventFiltersOnce = false;

  // Z defaults:
  cfg.massMin = 60;
  cfg.massMax = 120;
  cfg.nBins = 120;
  cfg.ptMin1 = 15;
  cfg.ptMin2 = 10;
  cfg.requireGood = false;
  cfg.requireGlobal = true;
  cfg.requirePF = false;
  cfg.requireTightID = true;
  cfg.applyVz = true;
  cfg.vzMax = 15;
  cfg.applyHiBin = false; // pO: leave OFF unless you confirm it's meaningful
  cfg.outPrefix = "ZToMuMu_pO2025";

  // If you later want J/psi quickly:
  // cfg.massMin=2.6; cfg.massMax=3.6; cfg.nBins=100;
  // cfg.ptMin1=4; cfg.ptMin2=3; cfg.outPrefix="JPsiToMuMu";

  // If you later want Upsilon quickly:
  // cfg.massMin=8.0; cfg.massMax=12.0; cfg.nBins=80;
  // cfg.ptMin1=6; cfg.ptMin2=4; cfg.outPrefix="UpsilonToMuMu";

  // --------------------
  // Open file / get trees
  // --------------------
  TFile *f = TFile::Open(fname);
  if (!f || f->IsZombie())
  {
    std::cerr << "ERROR: cannot open " << fname << "\n";
    return;
  }

  TTree *tMu = (TTree *)f->Get("ggHiNtuplizer/EventTree");
  if (!tMu)
  {
    std::cerr << "ERROR: missing ggHiNtuplizer/EventTree\n";
    return;
  }

  // HiTree is optional
  TTree *tHi = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  bool haveHiTree = (tHi != nullptr);

  // --------------------
  // Set up branches (muons)
  // --------------------
  Int_t nMu = 0;
  std::vector<float> *muPt = nullptr;
  std::vector<float> *muEta = nullptr;
  std::vector<float> *muPhi = nullptr;
  std::vector<int> *muCharge = nullptr;

  std::vector<int> *muIsGood = nullptr;
  std::vector<int> *muIsGlobal = nullptr;
  std::vector<int> *muIsPF = nullptr;
  std::vector<int> *muIDTight = nullptr;

  tMu->SetBranchAddress("nMu", &nMu);
  tMu->SetBranchAddress("muPt", &muPt);
  tMu->SetBranchAddress("muEta", &muEta);
  tMu->SetBranchAddress("muPhi", &muPhi);
  tMu->SetBranchAddress("muCharge", &muCharge);

  const bool has_muIsGood = HasBranch(tMu, "muIsGood");
  const bool has_muIsGlobal = HasBranch(tMu, "muIsGlobal");
  const bool has_muIsPF = HasBranch(tMu, "muIsPF");
  const bool has_muIDTight = HasBranch(tMu, "muIDTight");

  if (has_muIsGood)
    tMu->SetBranchAddress("muIsGood", &muIsGood);
  if (has_muIsGlobal)
    tMu->SetBranchAddress("muIsGlobal", &muIsGlobal);
  if (has_muIsPF)
    tMu->SetBranchAddress("muIsPF", &muIsPF);
  if (has_muIDTight)
    tMu->SetBranchAddress("muIDTight", &muIDTight);

  // --------------------
  // Optional event branches
  // --------------------
  Float_t vz = 999.f;
  Int_t hiBin = -999;

  const bool has_vz = (haveHiTree && HasBranch(tHi, "vz"));
  const bool has_hiBin = (haveHiTree && HasBranch(tHi, "hiBin"));

  if (has_vz)
    tHi->SetBranchAddress("vz", &vz);
  if (has_hiBin)
    tHi->SetBranchAddress("hiBin", &hiBin);

  // PF info

  std::vector<float> *muPFChIso = nullptr;
  std::vector<float> *muPFNeuIso = nullptr;
  std::vector<float> *muPFPhoIso = nullptr;

  tMu->SetBranchStatus("muPFChIso", 1);
  tMu->SetBranchStatus("muPFNeuIso", 1);
  tMu->SetBranchStatus("muPFPhoIso", 1);

  const bool has_muPFChIso = HasBranch(tMu, "muPFChIso");
  const bool has_muPFNeuIso = HasBranch(tMu, "muPFNeuIso");
  const bool has_muPFPhoIso = HasBranch(tMu, "muPFPhoIso");

  if (has_muPFChIso)
    tMu->SetBranchAddress("muPFChIso", &muPFChIso);
  if (has_muPFNeuIso)
    tMu->SetBranchAddress("muPFNeuIso", &muPFNeuIso);
  if (has_muPFPhoIso)
    tMu->SetBranchAddress("muPFPhoIso", &muPFPhoIso);

  // Here's the HLT info

  TTree *tHLTobj = (TTree *)f->Get("hltobject/HLT_OxyL1SingleMuOpen_v");

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

  TTree *tHLT = (TTree *)f->Get("hltanalysis/HltTree");

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

  TTree *tEvent = (TTree *)f->Get("skimanalysis/HltTree");

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

  // --------------------
  // Histogram
  // --------------------
  gStyle->SetOptStat(0);
  TH1D *hMass = new TH1D("hMass",
                         Form("%s; m_{#mu#mu} [GeV]; Events", cfg.outPrefix.c_str()),
                         cfg.nBins, cfg.massMin, cfg.massMax);

  TH1D *hMass_extended = new TH1D("hMass_extended",
                                  Form("%s; m_{#mu#mu} [GeV]; Events", cfg.outPrefix.c_str()),
                                  cfg.nBins, 15, 120);

  TH1D *hMass_vipul = new TH1D("hMass_vipul", Form("%s; m_{#mu#mu} [GeV]; Events", cfg.outPrefix.c_str()),
                               cfg.nBins, cfg.massMin, cfg.massMax);
  hMass->Sumw2();
  hMass_extended->Sumw2();
  hMass_vipul->Sumw2();

  // --------------------
  // Loop
  // --------------------
  const Long64_t nEntries = tMu->GetEntries();
  std::cout << "EventTree entries = " << nEntries << "\n";
  std::cout << "Have HiTree? " << (haveHiTree ? "YES" : "NO") << "\n";
  if (haveHiTree)
  {
    std::cout << "HiTree has vz? " << (has_vz ? "YES" : "NO") << "\n";
    std::cout << "HiTree has hiBin? " << (has_hiBin ? "YES" : "NO") << "\n";
  }

  Long64_t nPassEvent = 0;
  Long64_t nPassPair = 0;

  for (Long64_t ie = 0; ie < nEntries; ++ie)
  {

    if (ie % 200000 == 0)
      std::cout << "Event " << ie << "/" << nEntries << "\n";

    tMu->GetEntry(ie);
    tHi->GetEntry(ie);
    tHLT->GetEntry(ie);
    tHLTobj->GetEntry(ie);
    tEvent->GetEntry(ie);

    if (cfg.applyVz && has_vz)
    {
      if (TMath::Abs(vz) > cfg.vzMax)
        continue;
    }

    if (!PassEventSelection_pO(tHi,
                               warnedEventFiltersOnce,
                               has_pprimaryVertexFilter, pprimaryVertexFilter,
                               has_pclusterCompatibilityFilter, pclusterCompatibilityFilter))
      continue;

    if (cfg.applyHiBin && has_hiBin)
    {
      if (hiBin < cfg.hiBinMin || hiBin > cfg.hiBinMax)
        continue;
    }

    if (nMu < 2)
      continue;

    if (has_hlt)
    {
      if (!TriggerFired(HLT_OxyL1SingleMuOpen_v1))
        continue;
    }
    nPassEvent++; // skip all single muon events

    for (int i = 0; i < nMu; ++i)
    { // loop through all muon candiates

      bool passIsolead = true;

      if (!muPt || !muEta || !muPhi || !muCharge)
        continue; // if empty branch just skip

      if (muPt->at(i) < cfg.ptMin2)
        continue; // loose for looping
      if (TMath::Abs(muEta->at(i)) > cfg.etaMax)
        continue; // eta cuts

      if (cfg.requireGood && has_muIsGood && muIsGood->at(i) == 0)
        continue;
      if (cfg.requireGlobal && has_muIsGlobal && muIsGlobal->at(i) == 0)
        continue;
      if (cfg.requirePF && has_muIsPF && muIsPF->at(i) == 0)
        continue;
      if (cfg.requireTightID && has_muIDTight && muIDTight->at(i) == 0)
        continue;

      const double isoLead = RelIsoPF(i, muPt, muPFChIso, muPFNeuIso, muPFPhoIso);
      if (isoLead >= cfg.isoMax)
        passIsolead = false;

      for (int j = i + 1; j < nMu; ++j)
      {
        bool passIsosec = true;

        if (muPt->at(j) < cfg.ptMin2)
          continue;
        if (TMath::Abs(muEta->at(j)) > cfg.etaMax)
          continue;

        if (cfg.requireGood && has_muIsGood && muIsGood->at(j) == 0)
          continue;
        if (cfg.requireGlobal && has_muIsGlobal && muIsGlobal->at(j) == 0)
          continue;
        if (cfg.requirePF && has_muIsPF && muIsPF->at(j) == 0)
          continue;
        if (cfg.requireTightID && has_muIDTight && muIDTight->at(j) == 0)
          continue;

        if (cfg.requireOS)
        {
          if (muCharge->at(i) * muCharge->at(j) >= 0)
            continue;
        }

        // apply leading/subleading pT thresholds
        const double pt1 = muPt->at(i);
        const double pt2 = muPt->at(j);
        const double lead = (pt1 > pt2 ? pt1 : pt2);
        const double sub = (pt1 > pt2 ? pt2 : pt1);
        if (lead < cfg.ptMin1)
          continue;
        if (sub < cfg.ptMin2)
          continue;

        TLorentzVector v1, v2;
        v1.SetPtEtaPhiM(muPt->at(i), muEta->at(i), muPhi->at(i), MU_MASS());
        v2.SetPtEtaPhiM(muPt->at(j), muEta->at(j), muPhi->at(j), MU_MASS());
        const double m = (v1 + v2).M();

        const double isosub = RelIsoPF(j, muPt, muPFChIso, muPFNeuIso, muPFPhoIso);
        if (isosub >= cfg.isoMax)
          passIsosec = false;

        if (passIsolead && passIsosec)
        {
          hMass_extended->Fill(m);
          if (!(m < 60 || m > 120))
          {
            hMass->Fill(m);
          }
        }

        if (m < 60 || m > 120)
          continue;

        if (lead < 20 || sub < 20)
          continue;

        // Here's with tighter pT cut, without iso cut.

        hMass_vipul->Fill(m);
        nPassPair++;
      }
    }
  }

  std::cout << "Passed event preselection: " << nPassEvent << "\n";
  std::cout << "Found selected OS pair: " << nPassPair << "\n";

  // --------------------
  // Draw
  // --------------------
  /*TCanvas *c = new TCanvas("c", "c", 900, 700);
  hMass->SetLineWidth(2);
  hMass->Draw("hist");

  TPaveText *note = new TPaveText(0.55, 0.60, 0.88, 0.88, "NDC");
  note->SetFillStyle(0);
  note->SetBorderSize(0);
  note->AddText(cfg.outPrefix.c_str());
  note->AddText(Form("pT lead > %.1f, sub > %.1f GeV", cfg.ptMin1, cfg.ptMin2));
  note->AddText(Form("|#eta| < %.1f", cfg.etaMax));
  note->AddText(Form("Total # of Z : %lli", nPassPair));
  if (cfg.applyVz && has_vz)
    note->AddText(Form("|vz| < %.1f cm", cfg.vzMax));
  if (cfg.applyHiBin && has_hiBin)
    note->AddText(Form("hiBin in [%d,%d]", cfg.hiBinMin, cfg.hiBinMax));
  note->Draw();

  c->SaveAs((cfg.outPrefix + "_mass.png").c_str());
  c->SaveAs((cfg.outPrefix + "_mass.pdf").c_str());*/

  // optional: save histogram to a ROOT file
  TFile *fout = new TFile((cfg.outPrefix + "_hist.root").c_str(), "RECREATE");
  hMass->Write("", 2);
  hMass_extended->Write("", 2);
  hMass_vipul->Write("", 2);
  fout->Close();

  std::cout << "Wrote: " << cfg.outPrefix << "_mass.png/.pdf and _hist.root\n";
}