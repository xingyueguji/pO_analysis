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

static bool HasBranch(TTree* t, const char* name) {
  return (t && t->GetListOfBranches() && t->GetListOfBranches()->FindObject(name));
}

// --------------------
// User config struct
// --------------------
struct DileptonConfig {
  // resonance window (plot range)
  double massMin = 60.0;
  double massMax = 120.0;
  int    nBins   = 120;

  // selection
  double ptMin1  = 20.0;   // leading
  double ptMin2  = 10.0;   // subleading
  double etaMax  = 2.4;

  bool requireOS      = true;

  bool requireGood    = true;   // muIsGood
  bool requireGlobal  = true;   // muIsGlobal
  bool requirePF      = false;  // muIsPF
  bool requireTightID = false;  // muIDTight

  // optional event cuts (only applied if branch exists and enabled)
  bool   applyVz = true;
  double vzMax   = 15.0;

  bool   applyHiBin = false; // OFF by default for pO
  int    hiBinMin   = 0;
  int    hiBinMax   = 200;

  // output
  std::string outPrefix = "dimuon";
};

// --------------------
// Main
// --------------------
void DrawDimuonPeak(const char* fname = "root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2025.root")
{
  // ==============
  // Configure here
  // ==============
  DileptonConfig cfg;

  // Z defaults:
  cfg.massMin = 60; cfg.massMax = 120; cfg.nBins = 120;
  cfg.ptMin1 = 15; cfg.ptMin2 = 15;
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
  TFile* f = TFile::Open(fname);
  if (!f || f->IsZombie()) {
    std::cerr << "ERROR: cannot open " << fname << "\n";
    return;
  }

  TTree* tMu = (TTree*)f->Get("ggHiNtuplizer/EventTree");
  if (!tMu) {
    std::cerr << "ERROR: missing ggHiNtuplizer/EventTree\n";
    return;
  }

  // HiTree is optional
  TTree* tHi = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
  bool haveHiTree = (tHi != nullptr);

  // --------------------
  // Set up branches (muons)
  // --------------------
  Int_t nMu = 0;
  std::vector<float>* muPt  = nullptr;
  std::vector<float>* muEta = nullptr;
  std::vector<float>* muPhi = nullptr;
  std::vector<int>*   muCharge = nullptr;

  std::vector<int>* muIsGood   = nullptr;
  std::vector<int>* muIsGlobal = nullptr;
  std::vector<int>* muIsPF     = nullptr;
  std::vector<int>* muIDTight  = nullptr;

  tMu->SetBranchAddress("nMu", &nMu);
  tMu->SetBranchAddress("muPt", &muPt);
  tMu->SetBranchAddress("muEta", &muEta);
  tMu->SetBranchAddress("muPhi", &muPhi);
  tMu->SetBranchAddress("muCharge", &muCharge);

  const bool has_muIsGood   = HasBranch(tMu, "muIsGood");
  const bool has_muIsGlobal = HasBranch(tMu, "muIsGlobal");
  const bool has_muIsPF     = HasBranch(tMu, "muIsPF");
  const bool has_muIDTight  = HasBranch(tMu, "muIDTight");

  if (has_muIsGood)   tMu->SetBranchAddress("muIsGood", &muIsGood);
  if (has_muIsGlobal) tMu->SetBranchAddress("muIsGlobal", &muIsGlobal);
  if (has_muIsPF)     tMu->SetBranchAddress("muIsPF", &muIsPF);
  if (has_muIDTight)  tMu->SetBranchAddress("muIDTight", &muIDTight);

  // --------------------
  // Optional event branches
  // --------------------
  Float_t vz = 999.f;
  Int_t   hiBin = -999;

  const bool has_vz   = (haveHiTree && HasBranch(tHi, "vz"));
  const bool has_hiBin= (haveHiTree && HasBranch(tHi, "hiBin"));

  if (has_vz)    tHi->SetBranchAddress("vz", &vz);
  if (has_hiBin) tHi->SetBranchAddress("hiBin", &hiBin);

  // --------------------
  // Histogram
  // --------------------
  gStyle->SetOptStat(0);
  TH1D* hMass = new TH1D("hMass",
                         Form("%s; m_{#mu#mu} [GeV]; Events", cfg.outPrefix.c_str()),
                         cfg.nBins, cfg.massMin, cfg.massMax);
  hMass->Sumw2();

  // --------------------
  // Loop
  // --------------------
  const Long64_t nEntries = tMu->GetEntries();
  std::cout << "EventTree entries = " << nEntries << "\n";
  std::cout << "Have HiTree? " << (haveHiTree ? "YES" : "NO") << "\n";
  if (haveHiTree) {
    std::cout << "HiTree has vz? " << (has_vz ? "YES" : "NO") << "\n";
    std::cout << "HiTree has hiBin? " << (has_hiBin ? "YES" : "NO") << "\n";
  }

  Long64_t nPassEvent = 0;
  Long64_t nPassPair  = 0;

  for (Long64_t ie = 0; ie < nEntries; ++ie) {
    tMu->GetEntry(ie);
    if (haveHiTree) tHi->GetEntry(ie);

    if (cfg.applyVz && has_vz) {
      if (TMath::Abs(vz) > cfg.vzMax) continue;
    }

    if (cfg.applyHiBin && has_hiBin) {
      if (hiBin < cfg.hiBinMin || hiBin > cfg.hiBinMax) continue;
    }

    if (nMu < 2) continue;
    nPassEvent++; // skip all single muon events

    int bestI = -1, bestJ = -1;
    double bestScore = 1e30; // choose pair closest to center of window
    const double target = 0.5*(cfg.massMin + cfg.massMax);

    for (int i = 0; i < nMu; ++i) { // loop through all muon candiates
      if (!muPt || !muEta || !muPhi || !muCharge) continue; // if empty branch just skip

      if (muPt->at(i) < cfg.ptMin2) continue;              // loose for looping
      if (TMath::Abs(muEta->at(i)) > cfg.etaMax) continue; // eta cuts

      if (cfg.requireGood && has_muIsGood && muIsGood->at(i) == 0) continue; 
      if (cfg.requireGlobal && has_muIsGlobal && muIsGlobal->at(i) == 0) continue;
      if (cfg.requirePF && has_muIsPF && muIsPF->at(i) == 0) continue;
      if (cfg.requireTightID && has_muIDTight && muIDTight->at(i) == 0) continue;

      for (int j = i+1; j < nMu; ++j) {
        if (muPt->at(j) < cfg.ptMin2) continue;
        if (TMath::Abs(muEta->at(j)) > cfg.etaMax) continue;

        if (cfg.requireGood && has_muIsGood && muIsGood->at(j) == 0) continue;
        if (cfg.requireGlobal && has_muIsGlobal && muIsGlobal->at(j) == 0) continue;
        if (cfg.requirePF && has_muIsPF && muIsPF->at(j) == 0) continue;
        if (cfg.requireTightID && has_muIDTight && muIDTight->at(j) == 0) continue;

        if (cfg.requireOS) {
          if (muCharge->at(i) * muCharge->at(j) >= 0) continue;
        }

        // apply leading/subleading pT thresholds
        const double pt1 = muPt->at(i);
        const double pt2 = muPt->at(j);
        const double lead = (pt1 > pt2 ? pt1 : pt2);
        const double sub  = (pt1 > pt2 ? pt2 : pt1);
        if (lead < cfg.ptMin1) continue;
        if (sub  < cfg.ptMin2) continue;

        TLorentzVector v1, v2;
        v1.SetPtEtaPhiM(muPt->at(i), muEta->at(i), muPhi->at(i), MU_MASS());
        v2.SetPtEtaPhiM(muPt->at(j), muEta->at(j), muPhi->at(j), MU_MASS());
        const double m = (v1+v2).M();

        if (m < 60 || m > 120) continue;

        hMass->Fill(m);
        nPassPair++;

        // pick "best" pair closest to target (works for Z; still ok for others)
        /*const double score = TMath::Abs(m - target);
        if (score < bestScore) {
          bestScore = score;
          bestI = i; bestJ = j;
        }*/
      }
    }

    //if (bestI < 0 || bestJ < 0) continue;
    //nPassPair++;

    //TLorentzVector a, b;
    //a.SetPtEtaPhiM(muPt->at(bestI), muEta->at(bestI), muPhi->at(bestI), MU_MASS());
    //b.SetPtEtaPhiM(muPt->at(bestJ), muEta->at(bestJ), muPhi->at(bestJ), MU_MASS());
    //hMass->Fill((a+b).M());
  }

  std::cout << "Passed event preselection: " << nPassEvent << "\n";
  std::cout << "Found selected OS pair: " << nPassPair << "\n";

  // --------------------
  // Draw
  // --------------------
  TCanvas* c = new TCanvas("c","c",900,700);
  hMass->SetLineWidth(2);
  hMass->Draw("hist");

  TPaveText* note = new TPaveText(0.55,0.60,0.88,0.88,"NDC");
  note->SetFillStyle(0);
  note->SetBorderSize(0);
  note->AddText(cfg.outPrefix.c_str());
  note->AddText(Form("pT lead > %.1f, sub > %.1f GeV", cfg.ptMin1, cfg.ptMin2));
  note->AddText(Form("|#eta| < %.1f", cfg.etaMax));
  note->AddText(Form("Total # of Z : %lli", nPassPair));
  if (cfg.applyVz && has_vz) note->AddText(Form("|vz| < %.1f cm", cfg.vzMax));
  if (cfg.applyHiBin && has_hiBin) note->AddText(Form("hiBin in [%d,%d]", cfg.hiBinMin, cfg.hiBinMax));
  note->Draw();

  c->SaveAs((cfg.outPrefix + "_mass.png").c_str());
  c->SaveAs((cfg.outPrefix + "_mass.pdf").c_str());

  // optional: save histogram to a ROOT file
  TFile* fout = new TFile((cfg.outPrefix + "_hist.root").c_str(), "RECREATE");
  hMass->Write("",2);
  fout->Close();

  std::cout << "Wrote: " << cfg.outPrefix << "_mass.png/.pdf and _hist.root\n";
}