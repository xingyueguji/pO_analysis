// correction/njet_WZ.C
//
// Jet-multiplicity (njet) distributions for W-tagged and Z-tagged events,
// data vs signal MC.
//
// Jet definition (the point of this macro):
//   ak4PFJetAnalyzer/t jets with  jtpt > 30 GeV  &&  |jteta| < 2.1,
//   cleaned by DeltaR(jet, selected lepton) > 0.4 (AK4 cone) against the
//   selected lepton(s) -- essential for electrons, which are always also
//   reconstructed as a PF jet. jtpt is the calibrated jet pT (rawpt is the
//   uncorrected one). NO jet ID is applied beyond the kinematic cuts (the
//   tree carries the PF composition fractions jtPfCHF/... if one is wanted
//   later). The jet tree is entry-aligned with ggHiNtuplizer/EventTree
//   (verified run/evt-identical across the July-29 files; the macro enforces
//   equal entry counts and reads both trees in lockstep).
//
// Event selections replicate skim/skim.C exactly (the skim outputs carry no
// jet information, so this macro re-runs the selection on the ntuples):
//   W channels (Wmu/Wel): the full 8-step W selection at the NOMINAL cuts --
//     PF lepton pT>25, pO filters + |vz|<15, trigger, DY veto, Tight ID
//     (mu: muIDTight+PF, e: eleMVAIdWP95), leading-lepton pT>25 & |eta|<2.4,
//     Delta-beta relIso < 0.15/0.095, trigger match DR<0.4. The gen-reco
//     match is a no-op in the skim and is skipped here.
//     NB there is NO MET/m_T cut in the W selection, so the DATA W sample
//     contains sizeable QCD (~19% mu / ~50% e of the events, 2026-07-29
//     numbers) which peaks at HIGH njet -- read the nominal W data/MC njet
//     with that in mind. The `_mt40` twin (pT>25 && m_T>40, the lepton-pT-
//     discriminant selection, purity 84.5% mu / 69% e) is filled alongside
//     as the cleaner W-njet distribution.
//   Z channels (Zmm/Zee): iso-selected OS pair in the Z peak [60,120] GeV
//     (legs: pT>10, |eta|<2.4, ID'd (mu: Global+Tight, e: eleMVAIdWP95),
//     relIso < 0.2/0.095; leading leg pT>15; trigger bit). The event is
//     counted ONCE, on the first accepted pair (skim iteration order);
//     jets are cleaned against both legs.
//
// MC: signal only (W: Wp+Wm combined with the relative k_s from
// skim/mc_norm.h; Z: DY), gen-weighted like the skim; comparison is
// shape-normalized (SaveDataMCRatio), so only the W+/W- mix matters.
// Note the MC is UNEMBEDDED (no pO underlying event), but at jet pT > 30
// the UE contribution is small -- the dominant data/MC njet difference in
// the W channels is the QCD contamination.
//
// Outputs (run from correction/):
//   ./rootfile/njet_<chan>.root                     raw + combined histograms
//   ./plots/njet_<chan>/{njet,njet_mt40,jetpt,jeteta}.png/.pdf
//
// Run:
//   root -l -q 'njet_WZ.C+'                 # all four channels (full loops)
//   root -l -q 'njet_WZ.C+("Wmu")'          # or Wel / Zmm / Zee
//   root -l -q 'njet_WZ.C+("all", true)'    # plots-only: redraw from the
//                                           # stored rootfile/njet_<chan>.root
//                                           # (cosmetics iterations, no loops)

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../skim/skim_common.h"
#include "../skim/mc_norm.h"
#include "../plotting/plotting_helper.C"

using namespace pOSkim;

namespace
{

// -------- jet selection --------
const double kJetPtMin  = 30.0; // GeV (calibrated jtpt)
const double kJetEtaMax = 2.1;
const double kJetLepDR  = 0.4;  // AK4 cone: drop jets overlapping a selected lepton
const int    kNJetMax   = 10;   // njet axis [0, 10]; fills are capped at 10

// m_T cut of the `_mt40` W variant = the skim's lepton-pT-discriminant
// selection (naming contract: the literal "40" in the histogram names).
const double kMtCutW = 40.0;

// -------- jet-tree reader (ak4PFJetAnalyzer/t) --------
struct JetReader
{
  static const int kMaxJets = 200; // max nref in the July-29 files is 48
  Int_t   nref = 0;
  Float_t jtpt [kMaxJets];
  Float_t jteta[kMaxJets];
  Float_t jtphi[kMaxJets];
  bool warnedCap = false;

  bool Attach(TTree *t)
  {
    if (!t) return false;
    if (!(HasBranch(t, "nref") && HasBranch(t, "jtpt") &&
          HasBranch(t, "jteta") && HasBranch(t, "jtphi")))
      return false;
    t->SetBranchStatus("*", 0);
    t->SetBranchStatus("nref", 1);  t->SetBranchAddress("nref", &nref);
    t->SetBranchStatus("jtpt", 1);  t->SetBranchAddress("jtpt", jtpt);
    t->SetBranchStatus("jteta", 1); t->SetBranchAddress("jteta", jteta);
    t->SetBranchStatus("jtphi", 1); t->SetBranchAddress("jtphi", jtphi);
    return true;
  }

  // Count jets passing (pT, |eta|) and DR > kJetLepDR to EVERY lepton in
  // lepEtaPhi; optionally fill the per-jet pt/eta histos with weight w.
  int CountSelected(const std::vector<std::pair<double, double>> &lepEtaPhi,
                    double w, TH1D *hpt, TH1D *heta)
  {
    int n = nref;
    if (n > kMaxJets)
    {
      if (!warnedCap)
      {
        std::cout << "[WARN] nref=" << n << " > kMaxJets=" << kMaxJets
                  << "; capping (should never happen)\n";
        warnedCap = true;
      }
      n = kMaxJets;
    }
    int njet = 0;
    for (int i = 0; i < n; ++i)
    {
      if (jtpt[i] <= kJetPtMin) continue;
      if (std::abs(jteta[i]) >= kJetEtaMax) continue;
      bool overlap = false;
      for (const auto &lep : lepEtaPhi)
        if (DeltaR(jteta[i], jtphi[i], lep.first, lep.second) < kJetLepDR)
        {
          overlap = true;
          break;
        }
      if (overlap) continue;
      ++njet;
      if (hpt)  hpt ->Fill(jtpt[i],  w);
      if (heta) heta->Fill(jteta[i], w);
    }
    return njet;
  }
};

// -------- per-sample outputs --------
struct NJetOut
{
  TH1D *njet = nullptr, *jetpt = nullptr, *jeteta = nullptr;
  TH1D *njet_mt40 = nullptr; // W channels only
  long long nSelRaw = 0;     // selected events, raw
  double    nSelW   = 0.0;   // selected events, gen-weighted
};

NJetOut MakeOut(const std::string &tag, bool isW)
{
  NJetOut o;
  const std::string njTitle =
      Form(";N_{jets} (p_{T}>%.0f GeV, |#eta|<%.1f);Events", kJetPtMin, kJetEtaMax);
  o.njet   = new TH1D(Form("h_njet_%s", tag.c_str()),   njTitle.c_str(),
                      kNJetMax + 1, -0.5, kNJetMax + 0.5);
  o.jetpt  = new TH1D(Form("h_jetpt_%s", tag.c_str()),  ";p_{T}^{jet} [GeV];Jets",
                      25, 30, 130);
  o.jeteta = new TH1D(Form("h_jeteta_%s", tag.c_str()), ";#eta_{jet};Jets",
                      21, -kJetEtaMax, kJetEtaMax);
  o.njet->Sumw2(); o.jetpt->Sumw2(); o.jeteta->Sumw2();
  o.njet->SetDirectory(nullptr); o.jetpt->SetDirectory(nullptr);
  o.jeteta->SetDirectory(nullptr);
  if (isW)
  {
    o.njet_mt40 = new TH1D(Form("h_njet_mt40_%s", tag.c_str()), njTitle.c_str(),
                           kNJetMax + 1, -0.5, kNJetMax + 0.5);
    o.njet_mt40->Sumw2();
    o.njet_mt40->SetDirectory(nullptr);
  }
  return o;
}

// ============================================================
// W channels: full skim W selection (skim_Wmu/skim_Wel, nominal cuts),
// then count cleaned jets. Returns 0 on success.
// ============================================================
int RunW(bool isMu, SampleType sample, JetReader &jr, NJetOut &out)
{
  // -------- config = the skim's --------
  const double lepPt25     = 25.0;
  const double dyPtMin     = isMu ? 15.0 : 10.0;   // DY-veto leg pT (intentionally asymmetric)
  const double isoMax      = isMu ? 0.15 : 0.095;
  const double dyMassMin   = 80.0, dyMassMax = 110.0;
  const double etaMax      = 2.4;
  const double trigMatchDR = 0.4;  // 2026-08-12: tracks skim.C:283/963 (see there)
  const double vzMax       = 15.0;
  const double lepMass     = isMu ? MU_MASS : ELE_MASS;
  const char  *tagLep      = isMu ? "muon" : "electron";

  const bool isMC = IsMC(sample);
  std::string fname = kDefaultDataFile;
  if (isMC)
  {
    const auto info = ResolveMCSample(sample, isMu ? "mu" : "ele");
    if (info.fname.empty())
    {
      std::cerr << "[ERR] RunW: cannot resolve MC sample\n";
      return 1;
    }
    fname = info.fname;
  }

  std::cout << "[INPUT] " << fname << std::endl;
  TFile *f = TFile::Open(fname.c_str());
  if (!f || f->IsZombie()) { std::cerr << "[ERR] cannot open " << fname << "\n"; return 1; }

  TTree *tLep    = (TTree *)f->Get("ggHiNtuplizer/EventTree");
  TTree *tHi     = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  TTree *tPF     = (TTree *)f->Get("particleFlowAnalyser/pftree");
  TTree *tHLT    = (TTree *)f->Get("hltanalysis/HltTree");
  TTree *tHLTobj = (TTree *)f->Get(isMu ? "hltobject/HLT_OxyL1SingleMuOpen_v"
                                        : "hltobject/HLT_OxyL1SingleEG10_v");
  TTree *tEvent  = (TTree *)f->Get("skimanalysis/HltTree");
  TTree *tJet    = (TTree *)f->Get("ak4PFJetAnalyzer/t");
  if (!tLep || !tHi || !tPF || !tHLT || !tHLTobj || !tEvent)
  {
    std::cerr << "[FATAL] RunW: missing a required tree in " << fname << "\n";
    f->Close();
    return 2;
  }
  if (!tJet)
  {
    std::cerr << "[FATAL] RunW: missing ak4PFJetAnalyzer/t in " << fname << "\n";
    f->Close();
    return 2;
  }
  if (tJet->GetEntries() != tLep->GetEntries())
  {
    std::cerr << "[FATAL] RunW: jet tree (" << tJet->GetEntries()
              << ") and EventTree (" << tLep->GetEntries()
              << ") entry counts differ -- not entry-aligned\n";
    f->Close();
    return 2;
  }
  if (!jr.Attach(tJet))
  {
    std::cerr << "[FATAL] RunW: jet tree lacks nref/jtpt/jteta/jtphi\n";
    f->Close();
    return 2;
  }

  // -------- lepton branches (names differ by flavour, logic identical) --------
  const char *bN   = isMu ? "nMu"        : "nEle";
  const char *bPt  = isMu ? "muPt"       : "elePt";
  const char *bEta = isMu ? "muEta"      : "eleEta";
  const char *bPhi = isMu ? "muPhi"      : "elePhi";
  const char *bChg = isMu ? "muCharge"   : "eleCharge";
  const char *bID  = isMu ? "muIDTight"  : "eleMVAIdWP95";
  const char *bCh  = isMu ? "muPFChIso"  : "elePFChIso";
  const char *bNeu = isMu ? "muPFNeuIso" : "elePFNeuIso";
  const char *bPho = isMu ? "muPFPhoIso" : "elePFPhoIso";
  const char *bPU  = isMu ? "muPFPUIso"  : "elePFPUIso";

  Int_t nLep = 0;
  std::vector<float> *lepPt = nullptr, *lepEta = nullptr, *lepPhi = nullptr;
  std::vector<int>   *lepCharge = nullptr, *lepID = nullptr, *muIsPF = nullptr;
  std::vector<float> *chIso = nullptr, *neuIso = nullptr, *phoIso = nullptr, *puIso = nullptr;

  tLep->SetBranchStatus("*", 0);
  for (const char *bn : {bN, bPt, bEta, bPhi, bChg, bID, bCh, bNeu, bPho})
    if (!HasBranch(tLep, bn))
    {
      std::cerr << "[FATAL] RunW: missing EventTree branch " << bn << "\n";
      f->Close();
      return 2;
    }
  tLep->SetBranchStatus(bN, 1);   tLep->SetBranchAddress(bN,   &nLep);
  tLep->SetBranchStatus(bPt, 1);  tLep->SetBranchAddress(bPt,  &lepPt);
  tLep->SetBranchStatus(bEta, 1); tLep->SetBranchAddress(bEta, &lepEta);
  tLep->SetBranchStatus(bPhi, 1); tLep->SetBranchAddress(bPhi, &lepPhi);
  tLep->SetBranchStatus(bChg, 1); tLep->SetBranchAddress(bChg, &lepCharge);
  tLep->SetBranchStatus(bID, 1);  tLep->SetBranchAddress(bID,  &lepID);
  tLep->SetBranchStatus(bCh, 1);  tLep->SetBranchAddress(bCh,  &chIso);
  tLep->SetBranchStatus(bNeu, 1); tLep->SetBranchAddress(bNeu, &neuIso);
  tLep->SetBranchStatus(bPho, 1); tLep->SetBranchAddress(bPho, &phoIso);

  const bool has_puIso = HasBranch(tLep, bPU);
  if (has_puIso) { tLep->SetBranchStatus(bPU, 1); tLep->SetBranchAddress(bPU, &puIso); }
  else std::cout << "[WARN] RunW: no " << bPU << " branch; relIso uncorrected (no Delta-beta).\n";

  const bool has_muIsPF = isMu && HasBranch(tLep, "muIsPF");
  if (has_muIsPF) { tLep->SetBranchStatus("muIsPF", 1); tLep->SetBranchAddress("muIsPF", &muIsPF); }

  // -------- HLT bit --------
  const std::string hltNeedle = isMu ? "HLT_OxyL1SingleMuOpen_v1" : "HLT_OxyL1SingleEG10_v1";
  const std::string hltName   = FindBranchContaining(tHLT, hltNeedle);
  Int_t hltBit  = 0;
  bool  has_hlt = false;
  tHLT->SetBranchStatus("*", 0);
  if (!hltName.empty())
  {
    has_hlt = true;
    tHLT->SetBranchStatus(hltName.c_str(), 1);
    tHLT->SetBranchAddress(hltName.c_str(), &hltBit);
  }
  else
    std::cout << "[WARN] RunW: no branch containing " << hltNeedle
              << "; trigger treated as PASS.\n";

  // -------- HLT objects (trigger match) --------
  std::vector<double> *toPt = nullptr, *toEta = nullptr, *toPhi = nullptr;
  const bool has_toPt  = HasBranch(tHLTobj, "pt");
  const bool has_toEta = HasBranch(tHLTobj, "eta");
  const bool has_toPhi = HasBranch(tHLTobj, "phi");
  tHLTobj->SetBranchStatus("*", 0);
  if (has_toPt && has_toEta && has_toPhi)
  {
    tHLTobj->SetBranchStatus("pt", 1);  tHLTobj->SetBranchAddress("pt",  &toPt);
    tHLTobj->SetBranchStatus("eta", 1); tHLTobj->SetBranchAddress("eta", &toEta);
    tHLTobj->SetBranchStatus("phi", 1); tHLTobj->SetBranchAddress("phi", &toPhi);
  }
  else
    std::cout << "[WARN] RunW: could not find HLT objects.\n";

  // -------- event filters + vz + gen weight --------
  Int_t ppv = 1, pcc = 1;
  const bool has_ppv = HasBranch(tEvent, "pprimaryVertexFilter");
  const bool has_pcc = HasBranch(tEvent, "pclusterCompatibilityFilter");
  tEvent->SetBranchStatus("*", 0);
  if (has_ppv && has_pcc)
  {
    tEvent->SetBranchStatus("pprimaryVertexFilter", 1);
    tEvent->SetBranchStatus("pclusterCompatibilityFilter", 1);
    tEvent->SetBranchAddress("pprimaryVertexFilter", &ppv);
    tEvent->SetBranchAddress("pclusterCompatibilityFilter", &pcc);
  }

  Float_t vz = 999.f;
  tHi->SetBranchStatus("*", 0);
  const bool has_vz = HasBranch(tHi, "vz");
  if (has_vz) { tHi->SetBranchStatus("vz", 1); tHi->SetBranchAddress("vz", &vz); }

  Float_t    genWeight     = 1.f;
  const bool has_genWeight = isMC && HasBranch(tHi, "weight");
  if (has_genWeight) { tHi->SetBranchStatus("weight", 1); tHi->SetBranchAddress("weight", &genWeight); }
  else if (isMC) std::cout << "[WARN] RunW: MC but no HiTree 'weight'; filling unweighted.\n";

  // -------- PF tree (MET, for the _mt40 variant only) --------
  std::vector<float> *pfPt = nullptr, *pfPhi = nullptr;
  tPF->SetBranchStatus("*", 0);
  tPF->SetBranchStatus("pfPt", 1);  tPF->SetBranchAddress("pfPt",  &pfPt);
  tPF->SetBranchStatus("pfPhi", 1); tPF->SetBranchAddress("pfPhi", &pfPhi);

  // -------- event loop --------
  const Long64_t nEntries = tLep->GetEntries();
  std::cout << "Entries: " << nEntries << "\n";
  bool warnedFilters = false, warnedTrig = false;

  for (Long64_t ie = 0; ie < nEntries; ++ie)
  {
    if (ie % 500000 == 0) std::cout << "  event " << ie << "/" << nEntries << "\n";

    tLep->GetEntry(ie);
    tHi->GetEntry(ie);
    tHLT->GetEntry(ie);
    tHLTobj->GetEntry(ie);
    tEvent->GetEntry(ie);
    // tPF / tJet (the two big collections) are read only AFTER the event
    // passes the full selection -- same entry number, so lockstep is kept.

    const double w = has_genWeight ? (double)genWeight : 1.0;

    if (!lepPt || !lepEta || !lepPhi || !lepCharge || !lepID) continue;

    // (1) exists PF lepton pT > 25
    bool hasPF25 = false;
    for (int i = 0; i < nLep; ++i)
    {
      if (lepPt->at(i) <= lepPt25) continue;
      if (has_muIsPF && muIsPF && muIsPF->at(i) == 0) continue;
      hasPF25 = true;
      break;
    }
    if (!hasPF25) continue;

    // (2) pO event filters + vz
    if (!PassEventSelection_pO(warnedFilters, has_ppv, ppv, has_pcc, pcc)) continue;
    if (has_vz && TMath::Abs(vz) > vzMax) continue;

    // (3) trigger fired
    if (has_hlt && !TriggerFired(hltBit)) continue;

    // (4) DY veto: OS pair, both legs ID'd + isolated, mll in (80, 110)
    {
      std::vector<int> cand;
      cand.reserve(nLep);
      for (int i = 0; i < nLep; ++i)
      {
        if (lepPt->at(i) <= dyPtMin) continue;
        if (lepID->at(i) == 0) continue;
        if (has_muIsPF && muIsPF && muIsPF->at(i) == 0) continue;
        if (RelIsoPF(i, lepPt, chIso, neuIso, phoIso, puIso) >= isoMax) continue;
        cand.push_back(i);
      }
      bool veto = false;
      for (size_t a = 0; a < cand.size() && !veto; ++a)
        for (size_t b = a + 1; b < cand.size() && !veto; ++b)
        {
          const int i1 = cand[a], i2 = cand[b];
          if (lepCharge->at(i1) * lepCharge->at(i2) >= 0) continue;
          TLorentzVector p1, p2;
          p1.SetPtEtaPhiM(lepPt->at(i1), lepEta->at(i1), lepPhi->at(i1), lepMass);
          p2.SetPtEtaPhiM(lepPt->at(i2), lepEta->at(i2), lepPhi->at(i2), lepMass);
          const double mll = (p1 + p2).M();
          if (mll > dyMassMin && mll < dyMassMax) veto = true;
        }
      if (veto) continue;
    }

    // (5)+(6) leading ID'd (mu: +PF) lepton, pT > 25, |eta| < 2.4
    int iLead = -1;
    double bestPt = -1.0;
    for (int i = 0; i < nLep; ++i)
    {
      if (lepID->at(i) == 0) continue;
      if (has_muIsPF && muIsPF && muIsPF->at(i) == 0) continue;
      if (lepPt->at(i) > bestPt) { bestPt = lepPt->at(i); iLead = i; }
    }
    if (iLead < 0) continue;
    if (lepPt->at(iLead) <= lepPt25) continue;
    if (std::abs(lepEta->at(iLead)) > etaMax) continue;

    // (7) leading-lepton isolation
    if (RelIsoPF(iLead, lepPt, chIso, neuIso, phoIso, puIso) >= isoMax) continue;

    // (8) trigger match
    if (!PassLeadingLeptonTrigMatch(trigMatchDR, iLead, lepEta, lepPhi,
                                    has_toPt, toPt, has_toEta, toEta,
                                    has_toPhi, toPhi, warnedTrig, tagLep))
      continue;

    // -------- W-tagged: count jets --------
    tJet->GetEntry(ie);
    const std::vector<std::pair<double, double>> lepClean = {
        {(double)lepEta->at(iLead), (double)lepPhi->at(iLead)}};
    int njet = jr.CountSelected(lepClean, w, out.jetpt, out.jeteta);
    if (njet > kNJetMax) njet = kNJetMax; // keep everything in axis range
    out.njet->Fill(njet, w);
    out.nSelRaw++;
    out.nSelW += w;

    if (out.njet_mt40)
    {
      tPF->GetEntry(ie);
      const TVector2 metv = ComputePFMET(nullptr, pfPt, pfPhi);
      const double mt = TransverseMass(lepPt->at(iLead), lepPhi->at(iLead), metv);
      if (mt > kMtCutW) out.njet_mt40->Fill(njet, w);
    }
  }

  f->Close();
  delete f;
  std::cout << "[INFO] RunW selected " << out.nSelRaw << " events (weighted "
            << out.nSelW << ")\n";
  return 0;
}

// ============================================================
// Z channels: iso-selected OS pair in the Z peak [60,120] GeV
// (skim_Zmm/skim_Zee pair logic; event counted once, on the first
// accepted pair), then count jets cleaned against both legs.
// ============================================================
int RunZ(bool isMu, SampleType sample, JetReader &jr, NJetOut &out)
{
  // -------- config = the skim's --------
  const double ptMin1  = 15.0; // leading leg
  const double ptMin2  = 10.0; // both legs
  const double etaMax  = 2.4;
  const double isoMax  = isMu ? 0.2 : 0.095;
  const double vzMax   = 15.0;
  const double massMin = 60.0, massMax = 120.0;
  const double lepMass = isMu ? MU_MASS : ELE_MASS;

  const bool isMC = IsMC(sample);
  std::string fname = kDefaultDataFile;
  if (isMC)
  {
    const auto info = ResolveMCSample(sample, isMu ? "mu" : "ele");
    if (info.fname.empty())
    {
      std::cerr << "[ERR] RunZ: cannot resolve MC sample\n";
      return 1;
    }
    fname = info.fname;
  }

  std::cout << "[INPUT] " << fname << std::endl;
  TFile *f = TFile::Open(fname.c_str());
  if (!f || f->IsZombie()) { std::cerr << "[ERR] cannot open " << fname << "\n"; return 1; }

  TTree *tLep   = (TTree *)f->Get("ggHiNtuplizer/EventTree");
  TTree *tHLT   = (TTree *)f->Get("hltanalysis/HltTree");
  TTree *tEvent = (TTree *)f->Get("skimanalysis/HltTree");
  TTree *tJet   = (TTree *)f->Get("ak4PFJetAnalyzer/t");
  TTree *tHi    = (TTree *)f->Get("hiEvtAnalyzer/HiTree"); // optional (like the skim)
  const bool haveHiTree = (tHi != nullptr);
  if (!tLep || !tHLT || !tEvent)
  {
    std::cerr << "[FATAL] RunZ: missing a required tree in " << fname << "\n";
    f->Close();
    return 2;
  }
  if (!tJet)
  {
    std::cerr << "[FATAL] RunZ: missing ak4PFJetAnalyzer/t in " << fname << "\n";
    f->Close();
    return 2;
  }
  if (tJet->GetEntries() != tLep->GetEntries())
  {
    std::cerr << "[FATAL] RunZ: jet tree (" << tJet->GetEntries()
              << ") and EventTree (" << tLep->GetEntries()
              << ") entry counts differ -- not entry-aligned\n";
    f->Close();
    return 2;
  }
  if (!jr.Attach(tJet))
  {
    std::cerr << "[FATAL] RunZ: jet tree lacks nref/jtpt/jteta/jtphi\n";
    f->Close();
    return 2;
  }

  // -------- lepton branches --------
  const char *bN   = isMu ? "nMu"        : "nEle";
  const char *bPt  = isMu ? "muPt"       : "elePt";
  const char *bEta = isMu ? "muEta"      : "eleEta";
  const char *bPhi = isMu ? "muPhi"      : "elePhi";
  const char *bChg = isMu ? "muCharge"   : "eleCharge";
  const char *bID  = isMu ? "muIDTight"  : "eleMVAIdWP95";
  const char *bCh  = isMu ? "muPFChIso"  : "elePFChIso";
  const char *bNeu = isMu ? "muPFNeuIso" : "elePFNeuIso";
  const char *bPho = isMu ? "muPFPhoIso" : "elePFPhoIso";
  const char *bPU  = isMu ? "muPFPUIso"  : "elePFPUIso";

  Int_t nLep = 0;
  std::vector<float> *lepPt = nullptr, *lepEta = nullptr, *lepPhi = nullptr;
  std::vector<int>   *lepCharge = nullptr, *lepID = nullptr, *muIsGlobal = nullptr;
  std::vector<float> *chIso = nullptr, *neuIso = nullptr, *phoIso = nullptr, *puIso = nullptr;

  tLep->SetBranchStatus("*", 0);
  for (const char *bn : {bN, bPt, bEta, bPhi, bChg, bID, bCh, bNeu, bPho})
    if (!HasBranch(tLep, bn))
    {
      std::cerr << "[FATAL] RunZ: missing EventTree branch " << bn << "\n";
      f->Close();
      return 2;
    }
  tLep->SetBranchStatus(bN, 1);   tLep->SetBranchAddress(bN,   &nLep);
  tLep->SetBranchStatus(bPt, 1);  tLep->SetBranchAddress(bPt,  &lepPt);
  tLep->SetBranchStatus(bEta, 1); tLep->SetBranchAddress(bEta, &lepEta);
  tLep->SetBranchStatus(bPhi, 1); tLep->SetBranchAddress(bPhi, &lepPhi);
  tLep->SetBranchStatus(bChg, 1); tLep->SetBranchAddress(bChg, &lepCharge);
  tLep->SetBranchStatus(bID, 1);  tLep->SetBranchAddress(bID,  &lepID);
  tLep->SetBranchStatus(bCh, 1);  tLep->SetBranchAddress(bCh,  &chIso);
  tLep->SetBranchStatus(bNeu, 1); tLep->SetBranchAddress(bNeu, &neuIso);
  tLep->SetBranchStatus(bPho, 1); tLep->SetBranchAddress(bPho, &phoIso);

  const bool has_puIso = HasBranch(tLep, bPU);
  if (has_puIso) { tLep->SetBranchStatus(bPU, 1); tLep->SetBranchAddress(bPU, &puIso); }
  else std::cout << "[WARN] RunZ: no " << bPU << " branch; relIso uncorrected (no Delta-beta).\n";

  const bool has_muIsGlobal = isMu && HasBranch(tLep, "muIsGlobal");
  if (has_muIsGlobal) { tLep->SetBranchStatus("muIsGlobal", 1); tLep->SetBranchAddress("muIsGlobal", &muIsGlobal); }
  else if (isMu) std::cout << "[WARN] RunZ: no muIsGlobal branch; Global requirement not applied.\n";

  // -------- HLT bit --------
  const std::string hltNeedle = isMu ? "HLT_OxyL1SingleMuOpen_v1" : "HLT_OxyL1SingleEG10_v1";
  const std::string hltName   = FindBranchContaining(tHLT, hltNeedle);
  Int_t hltBit  = 0;
  bool  has_hlt = false;
  tHLT->SetBranchStatus("*", 0);
  if (!hltName.empty())
  {
    has_hlt = true;
    tHLT->SetBranchStatus(hltName.c_str(), 1);
    tHLT->SetBranchAddress(hltName.c_str(), &hltBit);
  }
  else
    std::cout << "[WARN] RunZ: no branch containing " << hltNeedle
              << "; trigger treated as PASS.\n";

  // -------- event filters + vz + gen weight --------
  Int_t ppv = 1, pcc = 1;
  const bool has_ppv = HasBranch(tEvent, "pprimaryVertexFilter");
  const bool has_pcc = HasBranch(tEvent, "pclusterCompatibilityFilter");
  tEvent->SetBranchStatus("*", 0);
  if (has_ppv && has_pcc)
  {
    tEvent->SetBranchStatus("pprimaryVertexFilter", 1);
    tEvent->SetBranchStatus("pclusterCompatibilityFilter", 1);
    tEvent->SetBranchAddress("pprimaryVertexFilter", &ppv);
    tEvent->SetBranchAddress("pclusterCompatibilityFilter", &pcc);
  }

  Float_t vz = 999.f;
  const bool has_vz = (haveHiTree && HasBranch(tHi, "vz"));
  if (haveHiTree) tHi->SetBranchStatus("*", 0);
  if (has_vz) { tHi->SetBranchStatus("vz", 1); tHi->SetBranchAddress("vz", &vz); }

  Float_t    genWeight     = 1.f;
  const bool has_genWeight = isMC && haveHiTree && HasBranch(tHi, "weight");
  if (has_genWeight) { tHi->SetBranchStatus("weight", 1); tHi->SetBranchAddress("weight", &genWeight); }
  else if (isMC) std::cout << "[WARN] RunZ: MC but no HiTree 'weight'; filling unweighted.\n";

  // -------- event loop --------
  const Long64_t nEntries = tLep->GetEntries();
  std::cout << "Entries: " << nEntries << "\n";
  bool warnedFilters = false;

  for (Long64_t ie = 0; ie < nEntries; ++ie)
  {
    if (ie % 500000 == 0) std::cout << "  event " << ie << "/" << nEntries << "\n";

    tLep->GetEntry(ie);
    if (haveHiTree) tHi->GetEntry(ie);
    tHLT->GetEntry(ie);
    tEvent->GetEntry(ie);
    // tJet read only after an accepted pair is found (lockstep entry number)

    const double w = has_genWeight ? (double)genWeight : 1.0;

    if (has_vz && TMath::Abs(vz) > vzMax) continue;
    if (!PassEventSelection_pO(warnedFilters, has_ppv, ppv, has_pcc, pcc)) continue;
    if (!lepPt || !lepEta || !lepPhi || !lepCharge || !lepID) continue;
    if (nLep < 2) continue;
    if (has_hlt && !TriggerFired(hltBit)) continue;

    // one ID'd + isolated OS pair in [60,120]; first accepted pair in the
    // skim's iteration order tags the event (>1 accepted pair is negligible)
    auto passLeg = [&](int i) -> bool
    {
      if (lepPt->at(i) < ptMin2) return false;
      if (TMath::Abs(lepEta->at(i)) > etaMax) return false;
      if (has_muIsGlobal && muIsGlobal->at(i) == 0) return false;
      if (lepID->at(i) == 0) return false;
      return true;
    };

    int iSel = -1, jSel = -1;
    for (int i = 0; i < nLep && iSel < 0; ++i)
    {
      if (!passLeg(i)) continue;
      if (RelIsoPF(i, lepPt, chIso, neuIso, phoIso, puIso) >= isoMax) continue;
      for (int j = i + 1; j < nLep; ++j)
      {
        if (!passLeg(j)) continue;
        if (lepCharge->at(i) * lepCharge->at(j) >= 0) continue;

        const double pt1  = lepPt->at(i), pt2 = lepPt->at(j);
        const double lead = (pt1 > pt2 ? pt1 : pt2);
        const double sub  = (pt1 > pt2 ? pt2 : pt1);
        if (lead < ptMin1 || sub < ptMin2) continue;

        TLorentzVector v1, v2;
        v1.SetPtEtaPhiM(lepPt->at(i), lepEta->at(i), lepPhi->at(i), lepMass);
        v2.SetPtEtaPhiM(lepPt->at(j), lepEta->at(j), lepPhi->at(j), lepMass);
        const double m = (v1 + v2).M();
        if (m < massMin || m > massMax) continue;

        if (RelIsoPF(j, lepPt, chIso, neuIso, phoIso, puIso) >= isoMax) continue;

        iSel = i;
        jSel = j;
        break;
      }
    }
    if (iSel < 0) continue;

    // -------- Z-tagged: count jets (cleaned against both legs) --------
    tJet->GetEntry(ie);
    const std::vector<std::pair<double, double>> lepClean = {
        {(double)lepEta->at(iSel), (double)lepPhi->at(iSel)},
        {(double)lepEta->at(jSel), (double)lepPhi->at(jSel)}};
    int njet = jr.CountSelected(lepClean, w, out.jetpt, out.jeteta);
    if (njet > kNJetMax) njet = kNJetMax; // keep everything in axis range
    out.njet->Fill(njet, w);
    out.nSelRaw++;
    out.nSelW += w;
  }

  f->Close();
  delete f;
  std::cout << "[INFO] RunZ selected " << out.nSelRaw << " events (weighted "
            << out.nSelW << ")\n";
  return 0;
}

// Combine two per-sample histos with relative scales (clone -- never mutate
// the raw inputs, they are written to the rootfile as-is).
TH1D *CombineMC(TH1D *hp, double kp, TH1D *hm, double km, const char *name)
{
  if (!hp) return nullptr;
  TH1D *h = (TH1D *)hp->Clone(name);
  h->SetDirectory(nullptr);
  h->Scale(kp);
  if (hm) h->Add(hm, km);
  return h;
}

// -------- per-channel plotting set (data + combined MC) --------
struct ChannelHists
{
  TH1D *dNjet = nullptr, *dNjetMt40 = nullptr, *dJetpt = nullptr, *dJeteta = nullptr;
  TH1D *mNjet = nullptr, *mNjetMt40 = nullptr, *mJetpt = nullptr, *mJeteta = nullptr;
};

// Read a histogram from an open file, cloned + detached (repo rule: never
// hold a file-owned histogram past the file's Close).
TH1D *GetCloned(TFile *f, const std::string &name, bool required)
{
  TH1D *h = (TH1D *)f->Get(name.c_str());
  if (!h)
  {
    if (required)
      std::cerr << "[ERR] missing histogram " << name << " in " << f->GetName() << "\n";
    return nullptr;
  }
  TH1D *c = (TH1D *)h->Clone((name + "_plot").c_str());
  c->SetDirectory(nullptr);
  return c;
}

// -------- plots + printed summary for one channel --------
void MakeChannelPlots(const std::string &ch, bool isW, bool isMu, const ChannelHists &H)
{
  const std::string outDir = "./plots/njet_" + ch;
  gSystem->mkdir(outDir.c_str(), kTRUE);

  const std::string lepSym = isMu ? "#mu" : "e";
  const std::string bosSym = isMu ? "#mu#mu" : "ee";

  // Header block at the LOWER RIGHT of the main pad (empty for these
  // steeply-falling distributions), clear of the data/MC content and of the
  // SaveDataMCRatio legend (top right); generous y-headroom keeps the frame
  // top open (user request 2026-08-05).
  PlotStyle ps;
  ps.yTitleOffset = 1.55;
  ps.headerX = 0.47;
  ps.headerY = 0.32;
  ps.titleSize = 0.04;
  ps.yHeadroom = 1.8;     // linear pads (jeteta)

  PlotStyle psLog = ps;
  psLog.logy = true;
  psLog.yHeadroom = 20.0; // log pads: ~1.3 decades of air above the peak

  const std::string header  = ch + (isW ? ": Data vs W MC" : ": Data vs DY MC");
  const std::string selNom  = isW ? "full W#rightarrow" + lepSym + "#nu selection"
                                  : "Z#rightarrow" + bosSym + " peak [60,120] GeV";
  const std::string jetDef  = Form("jets: p_{T}>%.0f, |#eta|<%.1f, #DeltaR(j,%s)>%.1f",
                                   kJetPtMin, kJetEtaMax, lepSym.c_str(), kJetLepDR);
  const std::string mcLabel = isW ? "W signal MC (W^{+}+W^{-})" : "DY signal MC";

  if (H.dNjet && H.mNjet)
  {
    SaveDataMCRatio(H.dNjet, H.mNjet, outDir + "/njet",
                    H.dNjet->GetXaxis()->GetTitle(), "Events (a.u.)",
                    header, selNom, jetDef, psLog,
                    /*normToData=*/true, "Data", mcLabel);
    std::cout << "[OK] wrote " << outDir << "/njet.png/.pdf\n";
  }
  if (isW && H.dNjetMt40 && H.mNjetMt40)
  {
    SaveDataMCRatio(H.dNjetMt40, H.mNjetMt40, outDir + "/njet_mt40",
                    H.dNjetMt40->GetXaxis()->GetTitle(), "Events (a.u.)",
                    header, selNom + ", m_{T}>40 GeV", jetDef, psLog,
                    /*normToData=*/true, "Data", mcLabel);
    std::cout << "[OK] wrote " << outDir << "/njet_mt40.png/.pdf\n";
  }
  if (H.dJetpt && H.mJetpt)
    SaveDataMCRatio(H.dJetpt, H.mJetpt, outDir + "/jetpt",
                    "p_{T}^{jet} [GeV]", "Jets (a.u.)",
                    header, selNom, jetDef, psLog,
                    /*normToData=*/true, "Data", mcLabel);
  if (H.dJeteta && H.mJeteta)
    SaveDataMCRatio(H.dJeteta, H.mJeteta, outDir + "/jeteta",
                    "#eta_{jet}", "Jets (a.u.)",
                    header, selNom, jetDef, ps,
                    /*normToData=*/true, "Data", mcLabel);
  std::cout << "[OK] wrote " << outDir << "/jetpt,jeteta .png/.pdf\n";

  // summary from the histograms (data is unweighted, so GetEntries = events;
  // the MC number is the k-mixed combined weighted total)
  auto frac1p = [](TH1D *h) -> double
  {
    const double tot = h->Integral();
    return (tot > 0) ? (tot - h->GetBinContent(1)) / tot : 0.0;
  };
  std::cout << "[SUMMARY] " << ch << " (jet: pT>" << kJetPtMin << ", |eta|<"
            << kJetEtaMax << ", DR>" << kJetLepDR << ")\n";
  std::cout << Form("  data    : %8.0f events   <njet> = %.4f   frac(>=1 jet) = %.4f\n",
                    H.dNjet->GetEntries(), H.dNjet->GetMean(), frac1p(H.dNjet));
  std::cout << Form("  sig MC  : %8.1f (wgt, k-mixed)   <njet> = %.4f   frac(>=1 jet) = %.4f\n",
                    H.mNjet->Integral(), H.mNjet->GetMean(), frac1p(H.mNjet));
  if (isW && H.dNjetMt40)
    std::cout << Form("  data mT>40: %6.0f events   <njet> = %.4f   frac(>=1 jet) = %.4f\n",
                      H.dNjetMt40->Integral(), H.dNjetMt40->GetMean(),
                      frac1p(H.dNjetMt40));
}

} // anonymous namespace

// ============================================================
// Driver
// ============================================================
void njet_WZ(const char *channel = "all", bool plotsOnly = false)
{
  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetOptStat(0);

  const std::string chArg = channel;
  std::vector<std::string> chans;
  if (chArg == "all")
    chans = {"Wmu", "Wel", "Zmm", "Zee"};
  else if (chArg == "Wmu" || chArg == "Wel" || chArg == "Zmm" || chArg == "Zee")
    chans = {chArg};
  else
  {
    std::cerr << "[ERR] channel must be \"Wmu\", \"Wel\", \"Zmm\", \"Zee\" or \"all\" (got \""
              << chArg << "\")\n";
    return;
  }

  gSystem->mkdir("rootfile", kTRUE);

  for (const auto &ch : chans)
  {
    const bool isW  = (ch[0] == 'W');
    const bool isMu = (ch == "Wmu" || ch == "Zmm");
    std::cout << "\n========== " << ch << " ==========\n";

    ChannelHists H;
    const std::string outRoot = "./rootfile/njet_" + ch + ".root";

    if (plotsOnly)
    {
      // redraw from the stored histograms -- no ntuple loops
      TFile *fin = TFile::Open(outRoot.c_str(), "READ");
      if (!fin || fin->IsZombie())
      {
        std::cerr << "[ERR] " << ch << ": cannot open " << outRoot
                  << " (run without plotsOnly first); skipping\n";
        continue;
      }
      H.dNjet   = GetCloned(fin, "h_njet_"   + ch + "_data", true);
      H.dJetpt  = GetCloned(fin, "h_jetpt_"  + ch + "_data", true);
      H.dJeteta = GetCloned(fin, "h_jeteta_" + ch + "_data", true);
      H.mNjet   = GetCloned(fin, "h_njet_"   + ch + "_mc", true);
      H.mJetpt  = GetCloned(fin, "h_jetpt_"  + ch + "_mc", true);
      H.mJeteta = GetCloned(fin, "h_jeteta_" + ch + "_mc", true);
      if (isW)
      {
        H.dNjetMt40 = GetCloned(fin, "h_njet_mt40_" + ch + "_data", true);
        H.mNjetMt40 = GetCloned(fin, "h_njet_mt40_" + ch + "_mc", true);
      }
      fin->Close();
      if (!H.dNjet || !H.mNjet)
      {
        std::cerr << "[ERR] " << ch << ": incomplete rootfile; skipping\n";
        continue;
      }
      MakeChannelPlots(ch, isW, isMu, H);
      continue;
    }

    JetReader jr; // one reader reused across the channel's files

    // ---- data ----
    NJetOut oData = MakeOut(ch + "_data", isW);
    if ((isW ? RunW(isMu, kData, jr, oData) : RunZ(isMu, kData, jr, oData)) != 0)
    {
      std::cerr << "[ERR] " << ch << ": data loop failed; skipping channel\n";
      continue;
    }

    // ---- signal MC (+ combined) ----
    NJetOut oWp, oWm, oDY;
    TH1D *hMC_njet = nullptr, *hMC_njet_mt40 = nullptr;
    TH1D *hMC_jetpt = nullptr, *hMC_jeteta = nullptr;
    if (isW)
    {
      oWp = MakeOut(ch + "_Wp", isW);
      oWm = MakeOut(ch + "_Wm", isW);
      if (RunW(isMu, kWp, jr, oWp) != 0 || RunW(isMu, kWm, jr, oWm) != 0)
      {
        std::cerr << "[ERR] " << ch << ": W MC loop failed; skipping channel\n";
        continue;
      }
      // relative W+/W- mix (shape norm downstream cancels the absolute scale)
      const std::string flav = isMu ? "mu" : "ele";
      const double kWp = pONorm::MCScale("Wp_" + flav);
      const double kWm = pONorm::MCScale("Wm_" + flav);
      std::cout << "[INFO] " << ch << " W+/W- mix: k_Wp=" << kWp << ", k_Wm=" << kWm << "\n";
      hMC_njet      = CombineMC(oWp.njet, kWp, oWm.njet, kWm, Form("h_njet_%s_mc", ch.c_str()));
      hMC_njet_mt40 = CombineMC(oWp.njet_mt40, kWp, oWm.njet_mt40, kWm, Form("h_njet_mt40_%s_mc", ch.c_str()));
      hMC_jetpt     = CombineMC(oWp.jetpt, kWp, oWm.jetpt, kWm, Form("h_jetpt_%s_mc", ch.c_str()));
      hMC_jeteta    = CombineMC(oWp.jeteta, kWp, oWm.jeteta, kWm, Form("h_jeteta_%s_mc", ch.c_str()));
    }
    else
    {
      oDY = MakeOut(ch + "_DY", isW);
      if (RunZ(isMu, kDY, jr, oDY) != 0)
      {
        std::cerr << "[ERR] " << ch << ": DY MC loop failed; skipping channel\n";
        continue;
      }
      // absolute k_s, so the *_mc histos carry the same (absolute) meaning
      // in the W and Z channels (shape norm downstream cancels it anyway)
      const double kDYs = pONorm::MCScale(isMu ? "DYmu" : "DYee");
      std::cout << "[INFO] " << ch << " DY scale: k_DY=" << kDYs << "\n";
      hMC_njet   = CombineMC(oDY.njet, kDYs, nullptr, 0.0, Form("h_njet_%s_mc", ch.c_str()));
      hMC_jetpt  = CombineMC(oDY.jetpt, kDYs, nullptr, 0.0, Form("h_jetpt_%s_mc", ch.c_str()));
      hMC_jeteta = CombineMC(oDY.jeteta, kDYs, nullptr, 0.0, Form("h_jeteta_%s_mc", ch.c_str()));
    }

    // ---- rootfile ----
    TFile *fout = TFile::Open(outRoot.c_str(), "RECREATE");
    auto writeOut = [&](const NJetOut &o)
    {
      if (o.njet)      o.njet->Write("", TObject::kOverwrite);
      if (o.jetpt)     o.jetpt->Write("", TObject::kOverwrite);
      if (o.jeteta)    o.jeteta->Write("", TObject::kOverwrite);
      if (o.njet_mt40) o.njet_mt40->Write("", TObject::kOverwrite);
    };
    writeOut(oData);
    if (isW) { writeOut(oWp); writeOut(oWm); }
    else     { writeOut(oDY); }
    if (hMC_njet)      hMC_njet->Write("", TObject::kOverwrite);
    if (hMC_njet_mt40) hMC_njet_mt40->Write("", TObject::kOverwrite);
    if (hMC_jetpt)     hMC_jetpt->Write("", TObject::kOverwrite);
    if (hMC_jeteta)    hMC_jeteta->Write("", TObject::kOverwrite);
    fout->Close();
    std::cout << "[OK] wrote " << outRoot << "\n";

    // ---- plots + summary ----
    H.dNjet  = oData.njet;   H.dNjetMt40 = oData.njet_mt40;
    H.dJetpt = oData.jetpt;  H.dJeteta   = oData.jeteta;
    H.mNjet  = hMC_njet;     H.mNjetMt40 = hMC_njet_mt40;
    H.mJetpt = hMC_jetpt;    H.mJeteta   = hMC_jeteta;
    MakeChannelPlots(ch, isW, isMu, H);
  }

  std::cout << "\n[DONE] njet_WZ\n";
}
