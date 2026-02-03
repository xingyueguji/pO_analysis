// DrawWToMuNu_PFMet_3Versions.C
// Run:
//   root -l -q 'DrawWToMuNu_PFMet_3Versions.C("root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2025.root")'

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TMath.h"
#include "TVector2.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TStyle.h"
#include <vector>
#include <iostream>
#include <string>
#include <cmath>

static bool HasBranch(TTree *t, const char *name)
{
  return (t && t->GetListOfBranches() && t->GetListOfBranches()->FindObject(name));
}

struct WConfig
{
  // muon selection
  double muPtMin = 25.0;
  double muEtaMax = 2.4;
  bool requireGlobal = true;
  bool requireGood = false;
  bool requireTight = true; // muIDTight
  bool requirePF = true;

  // event selection
  bool applyVz = true;
  double vzMax = 15.0;

  // do NOT use hiBin by default in pO
  bool applyHiBin = false;
  int hiBinMin = 0;
  int hiBinMax = 200;

  // simple W-shape cuts (optional)
  bool applyMetMin = false;
  double metMin = 20.0;

  // histogram
  int nBins = 100;
  double mtMin = 0.0;
  double mtMax = 200.0;

  std::string outPrefix = "WToMuNu_pO_PFMet";
};

// Compute PF MET from particleFlowAnalyser/pftree
static TVector2 ComputePFMET(std::vector<int> *pfId,
                             std::vector<float> *pfPt,
                             std::vector<float> *pfPhi)
{
  (void)pfId; // not used in this simple MET
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

static bool PassMuonID(int i,
                       const WConfig &cfg,
                       bool has_muIsGlobal, std::vector<int> *muIsGlobal,
                       bool has_muIsGood, std::vector<int> *muIsGood,
                       bool has_muIDTight, std::vector<int> *muIDTight,
                       bool has_muIsPF, std::vector<int> *muIsPF)
{
  if (cfg.requireGlobal && has_muIsGlobal && muIsGlobal && muIsGlobal->at(i) == 0)
    return false;
  if (cfg.requireGood && has_muIsGood && muIsGood && muIsGood->at(i) == 0)
    return false;
  if (cfg.requireTight && has_muIDTight && muIDTight && muIDTight->at(i) == 0)
    return false;
  if (cfg.requirePF && has_muIsPF && muIsPF && muIsPF->at(i) == 0)
    return false;
  return true;

  // This is actually redundant, same as just asking for TightID muon
}

// Find best muon (highest pT) passing kinematics+ID, and optionally isolation.
static int FindBestMuon(const WConfig &cfg,
                        bool applyIso, double isoMax,
                        Int_t nMu,
                        std::vector<float> *muPt,
                        std::vector<float> *muEta,
                        std::vector<float> *muPhi,
                        bool has_muIsGlobal, std::vector<int> *muIsGlobal,
                        bool has_muIsGood, std::vector<int> *muIsGood,
                        bool has_muIDTight, std::vector<int> *muIDTight,
                        bool has_muIsPF, std::vector<int> *muIsPF,
                        std::vector<float> *muPFChIso,
                        std::vector<float> *muPFNeuIso,
                        std::vector<float> *muPFPhoIso,
                        double &bestPtOut)
{
  int best = -1;
  double bestPt = -1.0;

  for (int i = 0; i < nMu; ++i)
  {
    if (!muPt || !muEta || !muPhi)
      continue;

    const double pt = muPt->at(i);
    const double eta = muEta->at(i);

    if (pt < cfg.muPtMin)
      continue;
    if (std::abs(eta) > cfg.muEtaMax)
      continue;

    if (!PassMuonID(i, cfg,
                    has_muIsGlobal, muIsGlobal,
                    has_muIsGood, muIsGood,
                    has_muIDTight, muIDTight,
                    has_muIsPF, muIsPF))
      continue;

    if (applyIso)
    {
      const double isoRel = RelIsoPF(i, muPt, muPFChIso, muPFNeuIso, muPFPhoIso);
      if (isoRel > isoMax)
        continue;
    }

    if (pt > bestPt)
    {
      bestPt = pt;
      best = i;
    }
  }

  bestPtOut = bestPt;
  return best;
}

// DY veto: reject event if there exists an OS pair where BOTH muons satisfy:
// pt>15, |eta|<cfg.muEtaMax, Tight+PF, relIso<isoMax (PF iso)
static bool PassDYVeto(const WConfig &cfg,
                       double isoMax,
                       Int_t nMu,
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
  if (!muPt || !muEta || !muCharge)
    return true; // can't apply -> don't veto

  std::vector<int> cand;
  cand.reserve(nMu);

  for (int i = 0; i < nMu; ++i)
  {
    const double pt = muPt->at(i);
    const double eta = muEta->at(i);

    if (pt < 15.0)
      continue;
    if (std::abs(eta) > cfg.muEtaMax)
      continue;

    if (cfg.requireTight && has_muIDTight && muIDTight && muIDTight->at(i) == 0)
      continue;
    if (cfg.requirePF && has_muIsPF && muIsPF && muIsPF->at(i) == 0)
      continue;

    const double isoRel = RelIsoPF(i, muPt, muPFChIso, muPFNeuIso, muPFPhoIso);
    if (isoRel > isoMax)
      continue;

    cand.push_back(i);
  }

  for (size_t a = 0; a < cand.size(); ++a)
  {
    for (size_t b = a + 1; b < cand.size(); ++b)
    {
      const int i1 = cand[a];
      const int i2 = cand[b];
      if (muCharge->at(i1) * muCharge->at(i2) < 0)
        return false; // found OS pair => FAIL veto

      TLorentzVector mu1, mu2;
      mu1.SetPtEtaPhiM(muPt->at(i1), muEta->at(i1), muPhi->at(i1), 0.105658);
      mu2.SetPtEtaPhiM(muPt->at(i2), muEta->at(i2), muPhi->at(i2), 0.105658);

      const double mll = (mu1 + mu2).M();

      // DY-like pair → veto event
      if (mll > 30.0)
        return false;
    }
  }

  return true; // no OS DY-like pair found => pass veto
}



void DrawWToMuNu_PFMet(const char *fname = "root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2025.root")
{
  gStyle->SetOptStat(0);

  WConfig cfg; // edit knobs here if you want
  const double isoMax = 0.15;

  TFile *f = TFile::Open(fname);
  if (!f || f->IsZombie())
  {
    std::cerr << "ERROR opening file\n";
    return;
  }

  TTree *tMu = (TTree *)f->Get("ggHiNtuplizer/EventTree");
  TTree *tHi = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  TTree *tPF = (TTree *)f->Get("particleFlowAnalyser/pftree");

  if (!tMu || !tHi || !tPF)
  {
    std::cerr << "ERROR: missing one of trees: ggHiNtuplizer/EventTree, hiEvtAnalyzer/HiTree, particleFlowAnalyser/pftree\n";
    return;
  }

  // -------- Muon branches
  Int_t nMu = 0;
  std::vector<float> *muPt = nullptr;
  std::vector<float> *muEta = nullptr;
  std::vector<float> *muPhi = nullptr;
  std::vector<int> *muCharge = nullptr;
  std::vector<int> *muIsGlobal = nullptr;
  std::vector<int> *muIsGood = nullptr;
  std::vector<int> *muIDTight = nullptr;
  std::vector<int> *muIsPF = nullptr;
  std::vector<float> *muPFChIso = nullptr;
  std::vector<float> *muPFNeuIso = nullptr;
  std::vector<float> *muPFPhoIso = nullptr;
  std::vector<float> *muPFPUIso = nullptr;

  tMu->SetBranchStatus("*", 0);
  tMu->SetBranchStatus("nMu", 1);
  tMu->SetBranchStatus("muPt", 1);
  tMu->SetBranchStatus("muEta", 1);
  tMu->SetBranchStatus("muPhi", 1);
  tMu->SetBranchStatus("muCharge", 1);
  tMu->SetBranchStatus("muIsGlobal", 1);
  tMu->SetBranchStatus("muIsGood", 1);
  tMu->SetBranchStatus("muIDTight", 1);
  tMu->SetBranchStatus("muIsPF", 1);
  tMu->SetBranchStatus("muPFChIso", 1);
  tMu->SetBranchStatus("muPFNeuIso", 1);
  tMu->SetBranchStatus("muPFPhoIso", 1);
  tMu->SetBranchStatus("muPFPUIso", 1);

  tMu->SetBranchAddress("nMu", &nMu);
  tMu->SetBranchAddress("muPt", &muPt);
  tMu->SetBranchAddress("muEta", &muEta);
  tMu->SetBranchAddress("muPhi", &muPhi);
  tMu->SetBranchAddress("muCharge", &muCharge);
  tMu->SetBranchAddress("muPFChIso", &muPFChIso);
  tMu->SetBranchAddress("muPFNeuIso", &muPFNeuIso);
  tMu->SetBranchAddress("muPFPhoIso", &muPFPhoIso);
  tMu->SetBranchAddress("muPFPUIso", &muPFPUIso);

  const bool has_muIsGlobal = HasBranch(tMu, "muIsGlobal");
  const bool has_muIsGood = HasBranch(tMu, "muIsGood");
  const bool has_muIDTight = HasBranch(tMu, "muIDTight");
  const bool has_muIsPF = HasBranch(tMu, "muIsPF");

  if (has_muIsGlobal)
    tMu->SetBranchAddress("muIsGlobal", &muIsGlobal);
  if (has_muIsGood)
    tMu->SetBranchAddress("muIsGood", &muIsGood);
  if (has_muIDTight)
    tMu->SetBranchAddress("muIDTight", &muIDTight);
  if (has_muIsPF)
    tMu->SetBranchAddress("muIsPF", &muIsPF);

  // -------- HiTree branches
  Float_t vz = 999.f;
  Int_t hiBin = -999;

  tHi->SetBranchStatus("*", 0);
  tHi->SetBranchStatus("vz", 1);
  tHi->SetBranchStatus("hiBin", 1);

  tHi->SetBranchAddress("vz", &vz);
  tHi->SetBranchAddress("hiBin", &hiBin);

  // -------- PF branches
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

  // -------- Output histograms (3 versions)
  TH1D *hMt_noIso_noDY = new TH1D("hMt_noIso_noDY",
                                  "; m_{T}(#mu, MET) [GeV]; Events",
                                  cfg.nBins, cfg.mtMin, cfg.mtMax);
  TH1D *hMt_isoOnly = new TH1D("hMt_isoOnly",
                               "; m_{T}(#mu, MET) [GeV]; Events",
                               cfg.nBins, cfg.mtMin, cfg.mtMax);
  TH1D *hMt_iso_DY = new TH1D("hMt_iso_DY",
                              "; m_{T}(#mu, MET) [GeV]; Events",
                              cfg.nBins, cfg.mtMin, cfg.mtMax);

  TH1D *hMet_noIso_noDY = new TH1D("hMet_noIso_noDY",
                                   "; PF MET [GeV]; Events",
                                   100, 0, 200);
  TH1D *hMet_isoOnly = new TH1D("hMet_isoOnly",
                                "; PF MET [GeV]; Events",
                                100, 0, 200);
  TH1D *hMet_iso_DY = new TH1D("hMet_iso_DY",
                               "; PF MET [GeV]; Events",
                               100, 0, 200);

  hMt_noIso_noDY->Sumw2();
  hMt_isoOnly->Sumw2();
  hMt_iso_DY->Sumw2();
  hMet_noIso_noDY->Sumw2();
  hMet_isoOnly->Sumw2();
  hMet_iso_DY->Sumw2();

  // Give different line styles (no explicit colors per your preference)
  hMt_noIso_noDY->SetLineWidth(2);
  hMt_isoOnly->SetLineWidth(2);
  hMt_iso_DY->SetLineWidth(2);
  hMt_noIso_noDY->SetLineStyle(1);
  hMt_isoOnly->SetLineStyle(2);
  hMt_iso_DY->SetLineStyle(3);
  hMt_noIso_noDY->SetLineColor(kAzure + 2);
  hMt_isoOnly->SetLineColor(kRed + 1);
  hMt_iso_DY->SetLineColor(kGreen + 2);

  hMet_noIso_noDY->SetLineWidth(2);
  hMet_isoOnly->SetLineWidth(2);
  hMet_iso_DY->SetLineWidth(2);
  hMet_noIso_noDY->SetLineStyle(1);
  hMet_isoOnly->SetLineStyle(2);
  hMet_iso_DY->SetLineStyle(3);
  hMet_noIso_noDY->SetLineColor(kAzure + 2);
  hMet_isoOnly->SetLineColor(kRed + 1);
  hMet_iso_DY->SetLineColor(kGreen + 2);

  // -------- Loop
  const Long64_t nEntries = tMu->GetEntries();
  std::cout << "Entries: " << nEntries << "\n";
  std::cout << "NOTE: using PF-MET computed from pftree (no dedicated MET branch found)\n";

  Long64_t nPassEvent = 0;
  Long64_t nPassMuon_noIso_noDY = 0, nPassMuon_isoOnly = 0, nPassMuon_iso_DY = 0;

  for (Long64_t ie = 0; ie < nEntries; ++ie)
  {
    if (ie % 100000 == 0)
      std::cout << "This is event " << ie << std::endl;

    tMu->GetEntry(ie);
    tHi->GetEntry(ie);
    tPF->GetEntry(ie);

    if (cfg.applyVz && std::abs(vz) > cfg.vzMax)
      continue;
    if (cfg.applyHiBin && (hiBin < cfg.hiBinMin || hiBin > cfg.hiBinMax))
      continue;

    nPassEvent++;

    // Compute PF MET once per event
    TVector2 metv = ComputePFMET(pfId, pfPt, pfPhi);
    const double met = metv.Mod();
    const double metPhi = metv.Phi();

    if (cfg.applyMetMin && met < cfg.metMin)
      continue;

    // ===============================
    // Version 1: No Iso, No DY veto
    // ===============================
    {
      double bestPt = -1.0;
      int best = FindBestMuon(cfg,
                              /*applyIso=*/false, isoMax,
                              nMu, muPt, muEta, muPhi,
                              has_muIsGlobal, muIsGlobal,
                              has_muIsGood, muIsGood,
                              has_muIDTight, muIDTight,
                              has_muIsPF, muIsPF,
                              muPFChIso, muPFNeuIso, muPFPhoIso,
                              bestPt);

      if (best >= 0)
      {
        nPassMuon_noIso_noDY++;
        const double mu_phi = muPhi->at(best);
        const double dphi = TVector2::Phi_mpi_pi(mu_phi - metPhi);
        const double mt = std::sqrt(2.0 * bestPt * met * (1.0 - std::cos(dphi)));

        hMet_noIso_noDY->Fill(met);
        hMt_noIso_noDY->Fill(mt);
      }
    }

    // ===============================
    // Version 2: Iso only (no DY veto)
    // ===============================
    {
      double bestPt = -1.0;
      int best = FindBestMuon(cfg,
                              /*applyIso=*/true, isoMax,
                              nMu, muPt, muEta, muPhi,
                              has_muIsGlobal, muIsGlobal,
                              has_muIsGood, muIsGood,
                              has_muIDTight, muIDTight,
                              has_muIsPF, muIsPF,
                              muPFChIso, muPFNeuIso, muPFPhoIso,
                              bestPt);

      if (best >= 0)
      {
        nPassMuon_isoOnly++;
        const double mu_phi = muPhi->at(best);
        const double dphi = TVector2::Phi_mpi_pi(mu_phi - metPhi);
        const double mt = std::sqrt(2.0 * bestPt * met * (1.0 - std::cos(dphi)));

        hMet_isoOnly->Fill(met);
        hMt_isoOnly->Fill(mt);
      }
    }

    // ===============================
    // Version 3: Iso + DY veto
    // ===============================
    {
      const bool passDY = PassDYVeto(cfg, isoMax,
                                     nMu, muPt, muEta, muPhi, muCharge,
                                     has_muIDTight, muIDTight,
                                     has_muIsPF, muIsPF,
                                     muPFChIso, muPFNeuIso, muPFPhoIso);

      if (passDY)
      {
        double bestPt = -1.0;
        int best = FindBestMuon(cfg,
                                /*applyIso=*/true, isoMax,
                                nMu, muPt, muEta, muPhi,
                                has_muIsGlobal, muIsGlobal,
                                has_muIsGood, muIsGood,
                                has_muIDTight, muIDTight,
                                has_muIsPF, muIsPF,
                                muPFChIso, muPFNeuIso, muPFPhoIso,
                                bestPt);

        if (best >= 0)
        {
          nPassMuon_iso_DY++;
          const double mu_phi = muPhi->at(best);
          const double dphi = TVector2::Phi_mpi_pi(mu_phi - metPhi);
          const double mt = std::sqrt(2.0 * bestPt * met * (1.0 - std::cos(dphi)));

          hMet_iso_DY->Fill(met);
          hMt_iso_DY->Fill(mt);
        }
      }
    }
  }

  std::cout << "Pass event cuts: " << nPassEvent << "\n";
  std::cout << "Have >=1 W muon candidate (NoIso, NoDY): " << nPassMuon_noIso_noDY << "\n";
  std::cout << "Have >=1 W muon candidate (Iso only):   " << nPassMuon_isoOnly << "\n";
  std::cout << "Have >=1 W muon candidate (Iso+DYVeto): " << nPassMuon_iso_DY << "\n";

  // -------- Drawing properties

  // Global style defaults (apply to all newly-created pads/canvases)
  gStyle->SetPadTopMargin(0.04);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.08);
  gStyle->SetPadBottomMargin(0.13);

  // ticks on both sides (x and y)
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  // Axis fonts (42 = Helvetica / ROOT default, not Times despite the comment)
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");

  // Axis label sizes
  gStyle->SetLabelSize(0.045, "XYZ");

  gStyle->SetTitleOffset(1.1, "X");
  gStyle->SetTitleOffset(1.6, "Y");

  gStyle->SetTitleSize(0.06, "XYZ");

  // -------- Draw mT overlay
  TCanvas *c1 = new TCanvas("c1", "mT (3 versions)", 800, 800);
  c1->SetLogy(false);

  // Make y-range nice
  double ymax = std::max({hMt_noIso_noDY->GetMaximum(), hMt_isoOnly->GetMaximum(), hMt_iso_DY->GetMaximum()});
  hMt_noIso_noDY->SetMaximum(1.25 * ymax);

  hMt_noIso_noDY->Draw("hist");
  hMt_isoOnly->Draw("hist same");
  hMt_iso_DY->Draw("hist same");

  TLegend *leg1 = new TLegend(0.55, 0.68, 0.88, 0.88);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->AddEntry(hMt_noIso_noDY, "No Iso, No DY veto", "l");
  leg1->AddEntry(hMt_isoOnly, "Iso < 0.15 only", "l");
  leg1->AddEntry(hMt_iso_DY, "Iso < 0.15 + DY veto", "l");
  leg1->Draw();

  TPaveText *p = new TPaveText(0.55, 0.43, 0.88, 0.68, "NDC");
  p->SetFillStyle(0);
  p->SetBorderSize(0);
  p->AddText(cfg.outPrefix.c_str());
  p->AddText(Form("mu p_{T}>%.0f, |#eta|<%.1f", cfg.muPtMin, cfg.muEtaMax));
  if (cfg.applyVz)
    p->AddText(Form("|v_{z}|<%.0f cm", cfg.vzMax));
  if (cfg.applyMetMin)
    p->AddText(Form("MET>%.0f GeV", cfg.metMin));
  p->AddText("MET: computed from PF candidates");
  p->Draw();

  c1->SaveAs(("./mT/" + cfg.outPrefix + "_mT_overlay.png").c_str());
  c1->SaveAs(("./mT/" + cfg.outPrefix + "_mT_overlay.pdf").c_str());

  // -------- Draw MET overlay
  TCanvas *c2 = new TCanvas("c2", "MET (3 versions)", 800, 800);
  double ymax2 = std::max({hMet_noIso_noDY->GetMaximum(), hMet_isoOnly->GetMaximum(), hMet_iso_DY->GetMaximum()});
  hMet_noIso_noDY->SetMaximum(1.25 * ymax2);

  hMet_noIso_noDY->Draw("hist");
  hMet_isoOnly->Draw("hist same");
  hMet_iso_DY->Draw("hist same");

  TLegend *leg2 = new TLegend(0.55, 0.68, 0.88, 0.88);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->AddEntry(hMet_noIso_noDY, "No Iso, No DY veto", "l");
  leg2->AddEntry(hMet_isoOnly, "Iso < 0.15 only", "l");
  leg2->AddEntry(hMet_iso_DY, "Iso < 0.15 + DY veto", "l");
  leg2->Draw();

  c2->SaveAs(("./met/" + cfg.outPrefix + "_met_overlay.png").c_str());
  c2->SaveAs(("./met/" + cfg.outPrefix + "_met_overlay.pdf").c_str());

  // -------- Save histograms
  TFile *fout = new TFile(("./rootfile/" + cfg.outPrefix + "_hist.root").c_str(), "RECREATE");
  hMt_noIso_noDY->Write("", 2);
  hMt_isoOnly->Write("", 2);
  hMt_iso_DY->Write("", 2);
  hMet_noIso_noDY->Write("", 2);
  hMet_isoOnly->Write("", 2);
  hMet_iso_DY->Write("", 2);
  fout->Close();

  std::cout << "Wrote outputs with prefix: " << cfg.outPrefix << "\n";
}