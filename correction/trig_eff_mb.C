// correction/trig_eff_mb.C
//
// SINGLE-LEPTON TRIGGER EFFICIENCY (turn-on) of the W analysis path, measured
// on the MINIMUM-BIAS-triggered sample, DATA vs W SIGNAL MC on the same plot.
//
// Definition (per the group's request, 2026-09-09):
//   den = event passes the analysis W selection  &&  MB path fired
//   num = den  &&  analysis single-lepton path fired
//              &&  leading lepton matched to a trigger object (DR < 0.4)
//   eps(pT) = num / den,   SF(pT) = eps_data / eps_MC
//
// "Analysis W selection" = the skim's 8-step W cutflow (skim/skim.C) WITHOUT
// steps 3 (trigger fired) and 8 (trigger match) -- those two ARE the quantity
// being measured. Steps 1/2/4/5/6/7 are replicated exactly as in
// charge_flip.C / njet_WZ.C (both verified event-identical with the skim).
// The lepton-pT floor of steps 1 and 6 is LOWERED to kPtFloor (10 GeV) so the
// actual turn-on is visible; the analysis selection is the pT > 25 sub-range
// (25 is a bin edge, so the two are the same histogram -- the leading tight
// lepton is the same object whenever it has pT > 25).
//
// Why the MB denominator works:
//   HLT_MinimumBiasHF_OR_BptxAND_v1 is a global HF-activity path. It knows
//   nothing about the lepton, so requiring it selects an UNBIASED subsample of
//   the offline-selected events -- including the ones the lepton path missed,
//   which are exactly the inefficiency. There is no object to match for the
//   MB path ("MB trigger matched" == the MB bit fired). In DATA the path was
//   prescaled run by run (fraction of lepton-triggered events with the MB bit:
//   393952 100%, 393953 91%, 393974/393975 0%, 393976 67%, 394004 51%,
//   394005/6 68%, 394007 78%; 49% overall, identical for muon- and electron-
//   triggered events within each run -> a pure prescale, uncorrelated with the
//   lepton; measured 2026-09-09). In MC it fires on 99.93% of W events (no
//   prescale). The prescale only costs statistics; the two MB-FRACTION control
//   plots below (den/all and num/trg vs pT) verify it is flat in lepton pT.
//
// Two selection variants, both filled in one loop:
//   nom   -- the plain W selection (NO MET / m_T cut; in data the low-pT part
//            is QCD-dominated: ~19% mu / ~50% e of the pT > 25 sample, more
//            below -- so the data turn-on is that of the MIXTURE)
//   mt40  -- + m_T > 40 GeV (PF MET from particleFlowAnalyser/pftree, read
//            only for selected events) = the leppt_mt40 discriminant selection,
//            QCD ~5% mu / ~28% e -> the cleaner data/MC comparison.
//
// What this is NOT: a tag-and-probe measurement. T&P (Z->ll) measures the
// per-LEPTON efficiency on a background-subtracted pure-lepton probe sample
// with Z kinematics; this measures the per-EVENT trigger efficiency of the W
// selection itself (W kinematics, W charge mix, the real fake content in
// data), with no tag bias, at the cost of the MB prescale and the data
// background composition. See the discussion in the README / AN.
//
// Outputs (run from correction/):
//   rootfile/trig_eff_mb_<mu|ele>.root  raw count histos per sample (+ MC
//                                       combined), TEfficiency objects, SF graphs
//   plots/trig_eff_mb_<mu|ele>/         turnon_pt_<sel>[_zoom25], eff_y_<sel>,
//                                       eff_y_<sel>_charge, eff_pt_<sel>_charge,
//                                       mbfrac_pt_<sel>, bit_vs_match_pt_<sel>,
//                                       trig_eff_<sel>.csv (per-bin table)
//   stdout                              the tables -- run through
//                                       ./run_trig_eff_mb.sh to keep the log
//
// Run:
//   ./run_trig_eff_mb.sh [mu|ele|both]           (keeps logs/trig_eff_mb_<chan>.log)
//   root -l -b -q 'trig_eff_mb.C+("mu")'          (bare; ~5 min data + ~1 min per MC file)
//   root -l -b -q 'trig_eff_mb.C+("both", true)'  (tables + plots only, from the rootfile)

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "../skim/skim_common.h"
#include "../skim/mc_norm.h"
#include "../plotting/plotting_helper.C"

using namespace pOSkim;

namespace
{

// -------- selection constants (tracking skim.C) --------
const double kPtFloor     = 10.0;  // loop floor for steps 1 and 6 (analysis = 25)
const double kPtNominal   = 25.0;  // the analysis cut; MUST be a bin edge of kPtEdges
const double kEtaMax      = 2.4;
const double kVzMax       = 15.0;
const double kTrigMatchDR = 0.4;   // skim.C:284/1008
const double kMtCut       = 40.0;  // the leppt_mt40 discriminant selection
const double kDyMassMin   = 80.0, kDyMassMax = 110.0;

// -------- trigger paths --------
const char *kMBPath = "HLT_MinimumBiasHF_OR_BptxAND_v1";

// -------- binning --------
const int    kNPt = 18;
const double kPtEdges[kNPt + 1] = {10, 12, 14, 16, 18, 20, 22, 25, 27.5, 30, 32.5, 35, 37.5, 40, 45, 50, 60, 80, 120};
const double kPtFoldMax = 119.9; // pT above the last edge is folded into the last bin

const char *kSel[2]      = {"nom", "mt40"};
const char *kSelLabel[2] = {"W selection (no m_{T} cut)", "W selection, m_{T} > 40 GeV"};

const double kCL = 0.6827;

// ============================================================
// Per-sample histogram set
// ============================================================
struct TrigOut
{
  std::string tag;
  std::map<std::string, TH1D *> h1;
  std::map<std::string, TH2D *> h2;
  double nSel[2] = {0, 0}, nMB[2] = {0, 0}, nTrg[2] = {0, 0}, nNum[2] = {0, 0}; // pT > 25 counts
  std::map<int, std::vector<double>> byRun;                                    // run -> {den, num} (nom, pT>25)

  TH1D *H1(const std::string &n) const
  {
    auto it = h1.find(n);
    if (it == h1.end()) { std::cerr << "[FATAL] TrigOut: no histo " << n << "\n"; std::abort(); }
    return it->second;
  }
  TH2D *H2(const std::string &n) const
  {
    auto it = h2.find(n);
    if (it == h2.end()) { std::cerr << "[FATAL] TrigOut: no histo2 " << n << "\n"; std::abort(); }
    return it->second;
  }

  void Book(const std::string &t)
  {
    tag = t;
    auto b1 = [&](const std::string &n, int nb, const double *e, const char *xt)
    {
      TH1D *h = new TH1D(Form("h_%s_%s", n.c_str(), tag.c_str()), Form(";%s;events", xt), nb, e);
      h->Sumw2();
      h->SetDirectory(nullptr);
      h1[n] = h;
    };
    auto b2 = [&](const std::string &n)
    {
      TH2D *h = new TH2D(Form("h_%s_%s", n.c_str(), tag.c_str()), ";p_{T} [GeV];y = -#eta_{lab}",
                         kNPt, kPtEdges, kNY, kYEdges);
      h->Sumw2();
      h->SetDirectory(nullptr);
      h2[n] = h;
    };
    for (int s = 0; s < 2; ++s)
    {
      const std::string S = kSel[s];
      // stages: all (analysis cuts only), den (+MB), bit (+MB+HLT bit), num (+MB+HLT+match),
      //         trg (analysis cuts + HLT + match, NO MB) -- for the MB-fraction control
      for (const char *st : {"all", "den", "bit", "num", "trg"})
      {
        b1(std::string(st) + "_pt_" + S, kNPt, kPtEdges, "leading-lepton p_{T} [GeV]");
        b1(std::string(st) + "_pt_" + S + "_plus", kNPt, kPtEdges, "leading-lepton p_{T} [GeV]");
        b1(std::string(st) + "_pt_" + S + "_minus", kNPt, kPtEdges, "leading-lepton p_{T} [GeV]");
        b1(std::string(st) + "_y_" + S, kNY, kYEdges, "y = -#eta_{lab}");          // pT > 25 only
        b1(std::string(st) + "_y_" + S + "_plus", kNY, kYEdges, "y = -#eta_{lab}");
        b1(std::string(st) + "_y_" + S + "_minus", kNY, kYEdges, "y = -#eta_{lab}");
      }
      // gen-weighted den/num (MC cross-check of the raw-count efficiency)
      b1("den_pt_" + S + "_w", kNPt, kPtEdges, "leading-lepton p_{T} [GeV]");
      b1("num_pt_" + S + "_w", kNPt, kPtEdges, "leading-lepton p_{T} [GeV]");
      b2("den_2d_" + S);
      b2("num_2d_" + S);
    }
  }

  // Fill one selected event. stage flags: mb, bit (HLT path fired), match.
  void Fill(int s, double pt, double y, int q, bool mb, bool bit, bool match, double w)
  {
    const std::string S = kSel[s];
    const double ptF = std::min(pt, kPtFoldMax);
    const bool   nomPt = (pt > kPtNominal);
    const char  *qs = (q > 0) ? "_plus" : "_minus";
    auto fill = [&](const char *st)
    {
      H1(std::string(st) + "_pt_" + S)->Fill(ptF);
      H1(std::string(st) + "_pt_" + S + qs)->Fill(ptF);
      if (nomPt)
      {
        H1(std::string(st) + "_y_" + S)->Fill(y);
        H1(std::string(st) + "_y_" + S + qs)->Fill(y);
      }
    };
    fill("all");
    if (bit && match) fill("trg");
    if (mb)
    {
      fill("den");
      H1("den_pt_" + S + "_w")->Fill(ptF, w);
      H2("den_2d_" + S)->Fill(ptF, y);
      if (bit) fill("bit");
      if (bit && match)
      {
        fill("num");
        H1("num_pt_" + S + "_w")->Fill(ptF, w);
        H2("num_2d_" + S)->Fill(ptF, y);
      }
    }
    if (nomPt)
    {
      nSel[s] += 1;
      if (bit && match) nTrg[s] += 1;
      if (mb) { nMB[s] += 1; if (bit && match) nNum[s] += 1; }
    }
  }

  void Write(TFile *f) const
  {
    f->cd();
    for (auto &p : h1) p.second->Write("", TObject::kOverwrite);
    for (auto &p : h2) p.second->Write("", TObject::kOverwrite);
  }
};

// ============================================================
// The event loop: W selection minus the trigger steps, on one file
// ============================================================
int RunTrig(bool isMu, SampleType sample, TrigOut &out)
{
  const bool   isMC      = IsMC(sample);
  const double dyPtMin   = isMu ? 15.0 : 10.0;   // DY-veto leg pT (intentionally asymmetric)
  const double isoMax    = isMu ? 0.15 : 0.095;
  const double lepMass   = isMu ? MU_MASS : ELE_MASS;
  const char  *tagLep    = isMu ? "muon" : "electron";

  std::string fname;
  if (isMC)
  {
    const auto info = ResolveMCSample(sample, isMu ? "mu" : "ele");
    if (info.fname.empty()) { std::cerr << "[ERR] RunTrig: cannot resolve MC sample\n"; return 1; }
    fname = info.fname;
  }
  else
    fname = kDefaultDataFile;

  std::cout << "[INPUT] " << out.tag << " <- " << fname << std::endl;
  TFile *f = TFile::Open(fname.c_str());
  if (!f || f->IsZombie()) { std::cerr << "[ERR] cannot open " << fname << "\n"; return 1; }

  TTree *tLep    = (TTree *)f->Get("ggHiNtuplizer/EventTree");
  TTree *tHi     = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  TTree *tPF     = (TTree *)f->Get("particleFlowAnalyser/pftree");
  TTree *tHLT    = (TTree *)f->Get("hltanalysis/HltTree");
  TTree *tHLTobj = (TTree *)f->Get(isMu ? "hltobject/HLT_OxyL1SingleMuOpen_v"
                                        : "hltobject/HLT_OxyL1SingleEG10_v");
  TTree *tEvent  = (TTree *)f->Get("skimanalysis/HltTree");
  if (!tLep || !tHi || !tPF || !tHLT || !tHLTobj || !tEvent)
  {
    std::cerr << "[FATAL] RunTrig: missing a required tree in " << fname << "\n";
    f->Close();
    return 2;
  }
  for (TTree *t : {tHi, tPF, tHLT, tHLTobj, tEvent})
    if (t->GetEntries() != tLep->GetEntries())
    {
      std::cerr << "[FATAL] RunTrig: tree " << t->GetName() << " (" << t->GetEntries()
                << ") and EventTree (" << tLep->GetEntries() << ") entry counts differ\n";
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
  std::vector<int>   *lepCharge = nullptr, *lepID = nullptr, *muIsPF = nullptr;
  std::vector<float> *chIso = nullptr, *neuIso = nullptr, *phoIso = nullptr, *puIso = nullptr;

  tLep->SetBranchStatus("*", 0);
  for (const char *bn : {bN, bPt, bEta, bPhi, bChg, bID, bCh, bNeu, bPho})
    if (!HasBranch(tLep, bn))
    {
      std::cerr << "[FATAL] RunTrig: missing EventTree branch " << bn << "\n";
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
  else std::cout << "[WARN] RunTrig: no " << bPU << " branch; relIso uncorrected (no Delta-beta).\n";

  const bool has_muIsPF = isMu && HasBranch(tLep, "muIsPF");
  if (has_muIsPF) { tLep->SetBranchStatus("muIsPF", 1); tLep->SetBranchAddress("muIsPF", &muIsPF); }

  // -------- HLT bits: the analysis path AND the MB path (both mandatory here) --------
  const std::string hltNeedle = isMu ? "HLT_OxyL1SingleMuOpen_v1" : "HLT_OxyL1SingleEG10_v1";
  const std::string hltName   = FindBranchContaining(tHLT, hltNeedle);
  const std::string mbName    = FindBranchContaining(tHLT, kMBPath);
  Int_t hltBit = 0, mbBit = 0;
  tHLT->SetBranchStatus("*", 0);
  if (hltName.empty() || mbName.empty())
  {
    std::cerr << "[FATAL] RunTrig: hltanalysis/HltTree lacks " << (hltName.empty() ? hltNeedle : std::string(kMBPath))
              << " -- the MB-denominator efficiency cannot be defined on this file.\n";
    f->Close();
    return 2;
  }
  // FindBranchContaining may return the _Prescale* twin; insist on the exact bit name.
  auto exact = [&](const std::string &needle) -> std::string
  {
    return HasBranch(tHLT, needle.c_str()) ? needle : FindBranchContaining(tHLT, needle);
  };
  const std::string hltExact = exact(hltNeedle), mbExact = exact(kMBPath);
  tHLT->SetBranchStatus(hltExact.c_str(), 1); tHLT->SetBranchAddress(hltExact.c_str(), &hltBit);
  tHLT->SetBranchStatus(mbExact.c_str(), 1);  tHLT->SetBranchAddress(mbExact.c_str(),  &mbBit);
  Int_t run = 0;
  const bool has_run = HasBranch(tHLT, "Run");
  if (has_run) { tHLT->SetBranchStatus("Run", 1); tHLT->SetBranchAddress("Run", &run); }
  std::cout << "[CONFIG] analysis path: " << hltExact << "   MB path: " << mbExact
            << "   match DR < " << kTrigMatchDR << "   pT floor " << kPtFloor << " (analysis " << kPtNominal << ")\n";

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
  {
    std::cerr << "[FATAL] RunTrig: no pt/eta/phi in " << tHLTobj->GetName() << " -- no trigger objects to match.\n";
    f->Close();
    return 2;
  }

  // -------- event filters + vz + gen weight --------
  Int_t ppv = 1, pcc = 1;
  const bool has_ppv = HasBranch(tEvent, "pprimaryVertexFilter");
  const bool has_pcc = HasBranch(tEvent, "pclusterCompatibilityFilter");
  tEvent->SetBranchStatus("*", 0);
  if (has_ppv) { tEvent->SetBranchStatus("pprimaryVertexFilter", 1);        tEvent->SetBranchAddress("pprimaryVertexFilter", &ppv); }
  if (has_pcc) { tEvent->SetBranchStatus("pclusterCompatibilityFilter", 1); tEvent->SetBranchAddress("pclusterCompatibilityFilter", &pcc); }

  Float_t vz = 999.f;
  tHi->SetBranchStatus("*", 0);
  const bool has_vz = HasBranch(tHi, "vz");
  if (has_vz) { tHi->SetBranchStatus("vz", 1); tHi->SetBranchAddress("vz", &vz); }

  Float_t    genWeight     = 1.f;
  const bool has_genWeight = isMC && HasBranch(tHi, "weight");
  if (has_genWeight) { tHi->SetBranchStatus("weight", 1); tHi->SetBranchAddress("weight", &genWeight); }
  else if (isMC) std::cout << "[WARN] RunTrig: no HiTree 'weight'; weighted sums unweighted.\n";

  // -------- PF tree (MET for the mt40 variant; read only for selected events) --------
  std::vector<float> *pfPt = nullptr, *pfPhi = nullptr;
  tPF->SetBranchStatus("*", 0);
  if (!HasBranch(tPF, "pfPt") || !HasBranch(tPF, "pfPhi"))
  {
    std::cerr << "[FATAL] RunTrig: pftree lacks pfPt/pfPhi\n";
    f->Close();
    return 2;
  }
  tPF->SetBranchStatus("pfPt", 1);  tPF->SetBranchAddress("pfPt",  &pfPt);
  tPF->SetBranchStatus("pfPhi", 1); tPF->SetBranchAddress("pfPhi", &pfPhi);

  // -------- event loop --------
  const Long64_t nEntries = tLep->GetEntries();
  std::cout << "Entries: " << nEntries << "\n";
  bool warnedFilters = false, warnedTrig = false;
  Long64_t nPass = 0;

  for (Long64_t ie = 0; ie < nEntries; ++ie)
  {
    if (ie % 500000 == 0) std::cout << "  event " << ie << "/" << nEntries << "\n";

    tLep->GetEntry(ie);
    tHi->GetEntry(ie);
    tHLT->GetEntry(ie);
    tHLTobj->GetEntry(ie);
    tEvent->GetEntry(ie);

    const double w = has_genWeight ? (double)genWeight : 1.0;

    if (!lepPt || !lepEta || !lepPhi || !lepCharge || !lepID) continue;

    // (1) exists PF lepton pT > floor (analysis: 25; the leading tight lepton
    //     with pT > 25 below implies the analysis version of this step)
    bool hasPF = false;
    for (int i = 0; i < nLep; ++i)
    {
      if (lepPt->at(i) <= kPtFloor) continue;
      if (has_muIsPF && muIsPF && muIsPF->at(i) == 0) continue;
      hasPF = true;
      break;
    }
    if (!hasPF) continue;

    // (2) pO event filters + vz
    if (!PassEventSelection_pO(warnedFilters, has_ppv, ppv, has_pcc, pcc)) continue;
    if (has_vz && TMath::Abs(vz) > kVzMax) continue;

    // (3) trigger fired  -- NOT applied: this is the numerator condition

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
          if (mll > kDyMassMin && mll < kDyMassMax) veto = true;
        }
      if (veto) continue;
    }

    // (5)+(6) leading ID'd (mu: +PF) lepton, pT > floor, |eta| < 2.4
    int iLead = -1;
    double bestPt = -1.0;
    for (int i = 0; i < nLep; ++i)
    {
      if (lepID->at(i) == 0) continue;
      if (has_muIsPF && muIsPF && muIsPF->at(i) == 0) continue;
      if (lepPt->at(i) > bestPt) { bestPt = lepPt->at(i); iLead = i; }
    }
    if (iLead < 0) continue;
    if (lepPt->at(iLead) <= kPtFloor) continue;
    if (std::abs(lepEta->at(iLead)) > kEtaMax) continue;

    // (7) leading-lepton isolation
    if (RelIsoPF(iLead, lepPt, chIso, neuIso, phoIso, puIso) >= isoMax) continue;

    // (8) trigger match -- NOT applied: numerator condition, evaluated below
    ++nPass;

    const bool mb    = TriggerFired(mbBit);
    const bool bit   = TriggerFired(hltBit);
    const bool match = bit && PassLeadingLeptonTrigMatch(kTrigMatchDR, iLead, lepEta, lepPhi,
                                                         has_toPt, toPt, has_toEta, toEta,
                                                         has_toPhi, toPhi, warnedTrig, tagLep);

    const double pt  = lepPt->at(iLead);
    const double y   = -lepEta->at(iLead);          // skim convention: p-going (-Z) = forward
    const int    q   = lepCharge->at(iLead);

    out.Fill(0, pt, y, q, mb, bit, match, w);

    // mt40 twin: PF MET only now (selected events are a tiny fraction)
    tPF->GetEntry(ie);
    const TVector2 metv = ComputePFMET(nullptr, pfPt, pfPhi);
    const double   mt   = TransverseMass(pt, lepPhi->at(iLead), metv);
    if (mt > kMtCut) out.Fill(1, pt, y, q, mb, bit, match, w);

    if (!isMC && has_run && pt > kPtNominal)
    {
      auto &v = out.byRun[run];
      if (v.empty()) v.assign(4, 0.0);
      v[0] += 1;                    // selected
      if (mb) v[1] += 1;            // den
      if (mb && match) v[2] += 1;   // num
      if (match) v[3] += 1;         // triggered (no MB)
    }
  }

  f->Close();
  delete f;
  std::cout << Form("[INFO] %s: %lld events pass the trigger-free W selection (pT > %.0f); pT > %.0f: sel %.0f, MB %.0f, num %.0f  |  mt40: sel %.0f, MB %.0f, num %.0f\n",
                    out.tag.c_str(), nPass, kPtFloor, kPtNominal,
                    out.nSel[0], out.nMB[0], out.nNum[0], out.nSel[1], out.nMB[1], out.nNum[1]);
  return 0;
}

// ============================================================
// Efficiency / ratio helpers
// ============================================================
TEfficiency *MakeEff(TH1 *pass, TH1 *tot, const std::string &name, const char *title = "")
{
  if (!pass || !tot) return nullptr;
  if (!TEfficiency::CheckConsistency(*pass, *tot))
  {
    std::cerr << "[ERR] TEfficiency inconsistency for " << name << "\n";
    return nullptr;
  }
  TEfficiency *e = new TEfficiency(*pass, *tot);
  e->SetName(name.c_str());
  e->SetTitle(title);
  e->SetStatisticOption(TEfficiency::kFCP);
  e->SetConfidenceLevel(kCL);
  e->SetDirectory(nullptr);
  return e;
}

// a / b per bin with the CP intervals propagated (a_up with b_dn and vice versa).
TGraphAsymmErrors *RatioGraph(TEfficiency *a, TEfficiency *b, const std::string &name)
{
  if (!a || !b) return nullptr;
  TGraphAsymmErrors *g = new TGraphAsymmErrors();
  g->SetName(name.c_str());
  const TH1 *ha = a->GetTotalHistogram();
  int n = 0;
  for (int i = 1; i <= ha->GetNbinsX(); ++i)
  {
    const double ea = a->GetEfficiency(i), eb = b->GetEfficiency(i);
    if (a->GetTotalHistogram()->GetBinContent(i) <= 0 || b->GetTotalHistogram()->GetBinContent(i) <= 0) continue;
    if (eb <= 0 || ea <= 0) continue;
    const double r  = ea / eb;
    const double ru = r * std::sqrt(std::pow(a->GetEfficiencyErrorUp(i) / ea, 2) + std::pow(b->GetEfficiencyErrorLow(i) / eb, 2));
    const double rd = r * std::sqrt(std::pow(a->GetEfficiencyErrorLow(i) / ea, 2) + std::pow(b->GetEfficiencyErrorUp(i) / eb, 2));
    const double x  = ha->GetBinCenter(i);
    const double xl = x - ha->GetBinLowEdge(i), xh = ha->GetBinLowEdge(i + 1) - x;
    g->SetPoint(n, x, r);
    g->SetPointError(n, xl, xh, rd, ru);
    ++n;
  }
  return g;
}

// "no vertical marker line" sentinel for DrawEffRatio -- must lie outside every
// x range used (the rapidity frames span [-2.4, 2.4], so -1 would NOT do).
const double kNoLine = -999.0;

struct Series
{
  TEfficiency *eff;
  std::string  label;
  int color, marker;
};
struct RatioSpec
{
  int iNum, iDen; // indices into the Series vector
  std::string label;
  int color, marker;
};

// Two-pad canvas: efficiencies on top, the requested ratios (SFs) below.
void DrawEffRatio(const std::vector<Series> &series, const std::vector<RatioSpec> &ratios,
                  const std::string &outPath, const std::string &xTitle, const std::string &ratioTitle,
                  const std::string &hdr, const std::string &sub1, const std::string &sub2,
                  const std::vector<std::string> &box,
                  double xlo, double xhi, double ylo, double yhi, double rlo, double rhi,
                  double xLine = kNoLine, const std::string &yTitle = "Trigger efficiency",
                  double legY1 = 0.37)
{
  // Layout: CMS_lumi(pad, 13, 10) paints "CMS / Work in Progress" top-LEFT
  // inside the frame (NDC y ~ 0.83-0.92), so the header block goes directly
  // below it across the full width; the legend and the info box sit lower
  // RIGHT, where the plateau points (eps ~ 1 -> NDC y ~ 0.6 for yhi = 1.6, or
  // ~0.55 on the zoomed 0.6-1.25 frames) never reach. Four-series plots pass a
  // lower legY1 and no box.
  gStyle->SetOptStat(0);
  PlotStyle ps;
  ps.headerX = 0.17; ps.headerY = 0.80; ps.headerDy = 0.045;
  ps.titleSize = 0.038; ps.subSize = 0.029;
  ps.boxX1 = 0.55; ps.boxX2 = 0.93; ps.boxY1 = 0.14; ps.boxY2 = 0.34; ps.boxTextSize = 0.027;

  static int nFrame = 0;
  TCanvas *c = new TCanvas(Form("c_eff_%d", nFrame), "", ps.w, ps.h);
  const bool haveRatio = !ratios.empty();
  const double split = haveRatio ? 0.30 : 0.0;
  TPad *pTop = new TPad("pTop", "", 0.0, split, 1.0, 1.0);
  TPad *pBot = haveRatio ? new TPad("pBot", "", 0.0, 0.0, 1.0, split) : nullptr;
  pTop->SetTopMargin(ps.tm); pTop->SetBottomMargin(haveRatio ? 0.02 : ps.bm);
  pTop->SetLeftMargin(ps.lm); pTop->SetRightMargin(ps.rm);
  pTop->SetTicks(1, 1);
  pTop->Draw();
  if (pBot)
  {
    pBot->SetTopMargin(0.04); pBot->SetBottomMargin(0.35);
    pBot->SetLeftMargin(ps.lm); pBot->SetRightMargin(ps.rm);
    pBot->SetTicks(1, 1);
    pBot->Draw();
  }

  // ---- top ----
  pTop->cd();
  TH1F *frame = new TH1F(Form("frame_eff_%d", nFrame++), "", 100, xlo, xhi);
  frame->SetDirectory(nullptr);
  ApplyHistStyle(frame, ps, haveRatio ? "" : xTitle, yTitle);
  if (haveRatio) { frame->GetXaxis()->SetLabelSize(0.0); frame->GetXaxis()->SetTitleSize(0.0); }
  frame->SetMinimum(ylo); frame->SetMaximum(yhi);
  frame->Draw("AXIS");

  TLegend *leg = new TLegend(0.55, legY1, 0.93, legY1 + 0.045 * series.size());
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(ps.font); leg->SetTextSize(0.030);
  std::vector<TGraphAsymmErrors *> gs;
  for (const auto &s : series)
  {
    TGraphAsymmErrors *g = s.eff ? s.eff->CreateGraph() : nullptr;
    gs.push_back(g);
    if (!g) continue;
    g->SetMarkerStyle(s.marker); g->SetMarkerSize(1.2);
    g->SetMarkerColor(s.color);  g->SetLineColor(s.color); g->SetLineWidth(2);
    g->Draw("P SAME");
    leg->AddEntry(g, s.label.c_str(), "lp");
  }
  leg->Draw();
  TLine *l1 = new TLine(xlo, 1.0, xhi, 1.0);
  l1->SetLineStyle(2); l1->SetLineColor(kGray + 2); l1->Draw();
  if (xLine > xlo && xLine < xhi)
  {
    TLine *lx = new TLine(xLine, ylo, xLine, yhi);
    lx->SetLineStyle(3); lx->SetLineColor(kBlue + 1); lx->SetLineWidth(2); lx->Draw();
  }
  DrawHeader(ps, hdr, sub1, sub2);
  DrawInfoBox(ps, box);
  CMS_lumi(pTop, 13, 10);
  pTop->RedrawAxis();

  // ---- bottom ----
  if (pBot)
  {
    pBot->cd();
    const double sf = 0.7 / 0.3;
    TH1F *fr = new TH1F(Form("frame_rat_%d", nFrame++), "", 100, xlo, xhi);
    fr->SetDirectory(nullptr);
    ApplyHistStyle(fr, ps, xTitle, ratioTitle);
    fr->GetXaxis()->SetTitleSize(ps.xTitleSize * sf); fr->GetYaxis()->SetTitleSize(ps.yTitleSize * sf);
    fr->GetXaxis()->SetLabelSize(ps.xLabelSize * sf); fr->GetYaxis()->SetLabelSize(ps.yLabelSize * sf);
    fr->GetXaxis()->SetTitleOffset(1.0); fr->GetYaxis()->SetTitleOffset(ps.yTitleOffset / sf);
    fr->GetYaxis()->SetNdivisions(505);
    fr->SetMinimum(rlo); fr->SetMaximum(rhi);
    fr->Draw("AXIS");
    TLine *lr = new TLine(xlo, 1.0, xhi, 1.0);
    lr->SetLineStyle(2); lr->SetLineColor(kRed + 1); lr->Draw();
    if (xLine > xlo && xLine < xhi)
    {
      TLine *lx = new TLine(xLine, rlo, xLine, rhi);
      lx->SetLineStyle(3); lx->SetLineColor(kBlue + 1); lx->SetLineWidth(2); lx->Draw();
    }
    TLegend *lg2 = ratios.size() > 1 ? new TLegend(0.55, 0.72, 0.93, 0.72 + 0.10 * ratios.size()) : nullptr;
    if (lg2) { lg2->SetBorderSize(0); lg2->SetFillStyle(0); lg2->SetTextFont(ps.font); lg2->SetTextSize(0.030 * sf); }
    for (const auto &r : ratios)
    {
      if (r.iNum >= (int)series.size() || r.iDen >= (int)series.size()) continue;
      TGraphAsymmErrors *g = RatioGraph(series[r.iNum].eff, series[r.iDen].eff, "r_tmp");
      if (!g) continue;
      g->SetMarkerStyle(r.marker); g->SetMarkerSize(1.2);
      g->SetMarkerColor(r.color);  g->SetLineColor(r.color); g->SetLineWidth(2);
      g->Draw("P SAME");
      if (lg2) lg2->AddEntry(g, r.label.c_str(), "lp");
    }
    if (lg2) lg2->Draw();
    pBot->RedrawAxis();
  }

  c->SaveAs((outPath + ".png").c_str());
  c->SaveAs((outPath + ".pdf").c_str());
  delete c;
}

// Sum of raw-count histograms of several samples (MC combined = Wp + Wm raw).
TH1D *SumH1(const std::vector<const TrigOut *> &outs, const std::string &key, const std::string &name)
{
  TH1D *h = nullptr;
  for (const TrigOut *o : outs)
  {
    TH1D *src = o->H1(key);
    if (!h) { h = (TH1D *)src->Clone(name.c_str()); h->SetDirectory(nullptr); }
    else h->Add(src);
  }
  return h;
}
TH2D *SumH2(const std::vector<const TrigOut *> &outs, const std::string &key, const std::string &name)
{
  TH2D *h = nullptr;
  for (const TrigOut *o : outs)
  {
    TH2D *src = o->H2(key);
    if (!h) { h = (TH2D *)src->Clone(name.c_str()); h->SetDirectory(nullptr); }
    else h->Add(src);
  }
  return h;
}

std::string FmtEff(TEfficiency *e, int bin)
{
  if (!e || e->GetTotalHistogram()->GetBinContent(bin) <= 0) return "      --          ";
  return Form("%.4f -%.4f +%.4f", e->GetEfficiency(bin), e->GetEfficiencyErrorLow(bin), e->GetEfficiencyErrorUp(bin));
}

// Inclusive efficiency over pT > 25 from the pT histograms (bins >= the 25 edge).
struct Incl { double den = 0, num = 0, eff = 0, lo = 0, hi = 0; };
Incl InclusiveAbove(TH1 *num, TH1 *den, double ptMin)
{
  Incl r;
  const int b0 = den->GetXaxis()->FindBin(ptMin + 1e-6);
  for (int i = b0; i <= den->GetNbinsX(); ++i) { r.den += den->GetBinContent(i); r.num += num->GetBinContent(i); }
  if (r.den > 0)
  {
    r.eff = r.num / r.den;
    r.lo  = r.eff - TEfficiency::ClopperPearson((int)r.den, (int)r.num, kCL, false);
    r.hi  = TEfficiency::ClopperPearson((int)r.den, (int)r.num, kCL, true) - r.eff;
  }
  return r;
}

// ============================================================
// Report: tables, TEfficiency objects, SF graphs, plots
// ============================================================
void Report(bool isMu, const TrigOut &oData, const TrigOut &oWp, const TrigOut &oWm, TFile *fout)
{
  const std::string ch      = isMu ? "mu" : "ele";
  const std::string lepSym  = isMu ? "#mu" : "e";
  const std::string procLbl = isMu ? "W #rightarrow #mu#nu" : "W #rightarrow e#nu";
  const std::string pathLbl = isMu ? "HLT_OxyL1SingleMuOpen" : "HLT_OxyL1SingleEG10";
  const std::string outDir  = "plots/trig_eff_mb_" + ch;
  gSystem->mkdir(outDir.c_str(), kTRUE);

  const std::vector<const TrigOut *> mcs = {&oWp, &oWm};

  std::cout << "\n==================== TRIGGER EFFICIENCY REPORT (" << ch << ") ====================\n";
  std::cout << "den = W selection (steps 1,2,4,5,6,7) && " << kMBPath << "\n";
  std::cout << "num = den && " << pathLbl << "_v1 && leading-lepton match DR < " << kTrigMatchDR << "\n";
  std::cout << "MC = Wp + Wm signal, RAW counts (gen weight ~ constant; weighted eff printed as a check)\n";

  for (int s = 0; s < 2; ++s)
  {
    const std::string S = kSel[s];
    std::cout << "\n---------- selection: " << S << "  (" << kSelLabel[s] << ") ----------\n";

    // ---- combined MC histograms ----
    std::map<std::string, TH1D *> M; // key -> summed MC histo
    for (const char *st : {"all", "den", "bit", "num", "trg"})
      for (const char *suf : {"", "_plus", "_minus"})
      {
        const std::string kpt = std::string(st) + "_pt_" + S + suf;
        const std::string ky  = std::string(st) + "_y_" + S + suf;
        M[kpt] = SumH1(mcs, kpt, "h_" + kpt + "_mc");
        M[ky]  = SumH1(mcs, ky, "h_" + ky + "_mc");
      }
    M["den_pt_" + S + "_w"] = SumH1(mcs, "den_pt_" + S + "_w", "h_den_pt_" + S + "_w_mc");
    M["num_pt_" + S + "_w"] = SumH1(mcs, "num_pt_" + S + "_w", "h_num_pt_" + S + "_w_mc");
    TH2D *mDen2 = SumH2(mcs, "den_2d_" + S, "h_den_2d_" + S + "_mc");
    TH2D *mNum2 = SumH2(mcs, "num_2d_" + S, "h_num_2d_" + S + "_mc");

    // ---- efficiencies ----
    auto effD = [&](const std::string &num, const std::string &den, const std::string &n)
    { return MakeEff(oData.H1(num), oData.H1(den), n); };
    auto effM = [&](const std::string &num, const std::string &den, const std::string &n)
    { return MakeEff(M[num], M[den], n); };

    TEfficiency *eD_pt  = effD("num_pt_" + S, "den_pt_" + S, "eff_pt_" + S + "_data");
    TEfficiency *eM_pt  = effM("num_pt_" + S, "den_pt_" + S, "eff_pt_" + S + "_mc");
    TEfficiency *eWp_pt = MakeEff(oWp.H1("num_pt_" + S), oWp.H1("den_pt_" + S), "eff_pt_" + S + "_Wp");
    TEfficiency *eWm_pt = MakeEff(oWm.H1("num_pt_" + S), oWm.H1("den_pt_" + S), "eff_pt_" + S + "_Wm");
    TEfficiency *eD_ptP = effD("num_pt_" + S + "_plus",  "den_pt_" + S + "_plus",  "eff_pt_" + S + "_plus_data");
    TEfficiency *eD_ptM = effD("num_pt_" + S + "_minus", "den_pt_" + S + "_minus", "eff_pt_" + S + "_minus_data");
    TEfficiency *eM_ptP = effM("num_pt_" + S + "_plus",  "den_pt_" + S + "_plus",  "eff_pt_" + S + "_plus_mc");
    TEfficiency *eM_ptM = effM("num_pt_" + S + "_minus", "den_pt_" + S + "_minus", "eff_pt_" + S + "_minus_mc");
    TEfficiency *eD_y   = effD("num_y_" + S, "den_y_" + S, "eff_y_" + S + "_data");
    TEfficiency *eM_y   = effM("num_y_" + S, "den_y_" + S, "eff_y_" + S + "_mc");
    TEfficiency *eD_yP  = effD("num_y_" + S + "_plus",  "den_y_" + S + "_plus",  "eff_y_" + S + "_plus_data");
    TEfficiency *eD_yM  = effD("num_y_" + S + "_minus", "den_y_" + S + "_minus", "eff_y_" + S + "_minus_data");
    TEfficiency *eM_yP  = effM("num_y_" + S + "_plus",  "den_y_" + S + "_plus",  "eff_y_" + S + "_plus_mc");
    TEfficiency *eM_yM  = effM("num_y_" + S + "_minus", "den_y_" + S + "_minus", "eff_y_" + S + "_minus_mc");
    // HLT bit only (no object match) -- separates the path from the matching
    TEfficiency *eD_bit = effD("bit_pt_" + S, "den_pt_" + S, "effbit_pt_" + S + "_data");
    TEfficiency *eM_bit = effM("bit_pt_" + S, "den_pt_" + S, "effbit_pt_" + S + "_mc");
    // MB-fraction controls: den/all (no trigger requirement) and num/trg (triggered)
    TEfficiency *fD_all = effD("den_pt_" + S, "all_pt_" + S, "mbfrac_all_pt_" + S + "_data");
    TEfficiency *fD_trg = effD("num_pt_" + S, "trg_pt_" + S, "mbfrac_trg_pt_" + S + "_data");
    TEfficiency *fM_all = effM("den_pt_" + S, "all_pt_" + S, "mbfrac_all_pt_" + S + "_mc");
    TEfficiency *fM_trg = effM("num_pt_" + S, "trg_pt_" + S, "mbfrac_trg_pt_" + S + "_mc");
    TEfficiency *e2D_D  = MakeEff(oData.H2("num_2d_" + S), oData.H2("den_2d_" + S), "eff_2d_" + S + "_data");
    TEfficiency *e2D_M  = MakeEff(mNum2, mDen2, "eff_2d_" + S + "_mc");

    TGraphAsymmErrors *sf_pt = RatioGraph(eD_pt, eM_pt, "sf_pt_" + S);
    TGraphAsymmErrors *sf_y  = RatioGraph(eD_y, eM_y, "sf_y_" + S);
    TGraphAsymmErrors *sf_yP = RatioGraph(eD_yP, eM_yP, "sf_y_" + S + "_plus");
    TGraphAsymmErrors *sf_yM = RatioGraph(eD_yM, eM_yM, "sf_y_" + S + "_minus");

    // ---- per-bin table (pT) ----
    std::cout << Form("%-12s | %9s %9s %-24s | %9s %9s %-24s | %-22s | %s\n",
                      "pT bin", "data den", "data num", "eff_data", "MC den", "MC num", "eff_MC", "SF = data/MC", "MB frac (data, all|trg)");
    std::ofstream csv(outDir + "/trig_eff_" + S + ".csv");
    csv << "ptlo,pthi,data_den,data_num,eff_data,eff_data_lo,eff_data_hi,mc_den,mc_num,eff_mc,eff_mc_lo,eff_mc_hi,sf,sf_lo,sf_hi\n";
    for (int i = 1; i <= kNPt; ++i)
    {
      const double lo = kPtEdges[i - 1], hi = kPtEdges[i];
      std::string sfStr = "        --          ";
      double sfv = 0, sfl = 0, sfh = 0;
      if (sf_pt)
        for (int k = 0; k < sf_pt->GetN(); ++k)
          if (std::abs(sf_pt->GetX()[k] - 0.5 * (lo + hi)) < 1e-6)
          {
            sfv = sf_pt->GetY()[k]; sfl = sf_pt->GetEYlow()[k]; sfh = sf_pt->GetEYhigh()[k];
            sfStr = Form("%.4f -%.4f +%.4f", sfv, sfl, sfh);
          }
      const double fa = fD_all ? fD_all->GetEfficiency(i) : 0, ft = fD_trg ? fD_trg->GetEfficiency(i) : 0;
      std::cout << Form("%5.1f-%-6.1f | %9.0f %9.0f %-24s | %9.0f %9.0f %-24s | %-22s | %.3f | %.3f\n",
                        lo, hi,
                        oData.H1("den_pt_" + S)->GetBinContent(i), oData.H1("num_pt_" + S)->GetBinContent(i), FmtEff(eD_pt, i).c_str(),
                        M["den_pt_" + S]->GetBinContent(i), M["num_pt_" + S]->GetBinContent(i), FmtEff(eM_pt, i).c_str(),
                        sfStr.c_str(), fa, ft);
      csv << Form("%g,%g,%.0f,%.0f,%.6f,%.6f,%.6f,%.0f,%.0f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n", lo, hi,
                  oData.H1("den_pt_" + S)->GetBinContent(i), oData.H1("num_pt_" + S)->GetBinContent(i),
                  eD_pt ? eD_pt->GetEfficiency(i) : 0, eD_pt ? eD_pt->GetEfficiencyErrorLow(i) : 0, eD_pt ? eD_pt->GetEfficiencyErrorUp(i) : 0,
                  M["den_pt_" + S]->GetBinContent(i), M["num_pt_" + S]->GetBinContent(i),
                  eM_pt ? eM_pt->GetEfficiency(i) : 0, eM_pt ? eM_pt->GetEfficiencyErrorLow(i) : 0, eM_pt ? eM_pt->GetEfficiencyErrorUp(i) : 0,
                  sfv, sfl, sfh);
    }
    csv.close();

    // ---- per-bin table (y, pT > 25) ----
    std::cout << Form("\n%-14s | %9s %9s %-24s | %9s %9s %-24s | %-22s\n",
                      "y bin (pT>25)", "data den", "data num", "eff_data", "MC den", "MC num", "eff_MC", "SF = data/MC");
    for (int i = 1; i <= kNY; ++i)
    {
      std::string sfStr = "        --          ";
      if (sf_y)
        for (int k = 0; k < sf_y->GetN(); ++k)
          if (std::abs(sf_y->GetX()[k] - 0.5 * (kYEdges[i - 1] + kYEdges[i])) < 1e-6)
            sfStr = Form("%.4f -%.4f +%.4f", sf_y->GetY()[k], sf_y->GetEYlow()[k], sf_y->GetEYhigh()[k]);
      std::cout << Form("%5.2f..%-6.2f | %9.0f %9.0f %-24s | %9.0f %9.0f %-24s | %-22s\n",
                        kYEdges[i - 1], kYEdges[i],
                        oData.H1("den_y_" + S)->GetBinContent(i), oData.H1("num_y_" + S)->GetBinContent(i), FmtEff(eD_y, i).c_str(),
                        M["den_y_" + S]->GetBinContent(i), M["num_y_" + S]->GetBinContent(i), FmtEff(eM_y, i).c_str(), sfStr.c_str());
    }

    // ---- inclusive pT > 25 ----
    const Incl iD  = InclusiveAbove(oData.H1("num_pt_" + S), oData.H1("den_pt_" + S), kPtNominal);
    const Incl iM  = InclusiveAbove(M["num_pt_" + S], M["den_pt_" + S], kPtNominal);
    const Incl iDP = InclusiveAbove(oData.H1("num_pt_" + S + "_plus"),  oData.H1("den_pt_" + S + "_plus"),  kPtNominal);
    const Incl iDM = InclusiveAbove(oData.H1("num_pt_" + S + "_minus"), oData.H1("den_pt_" + S + "_minus"), kPtNominal);
    const Incl iMP = InclusiveAbove(M["num_pt_" + S + "_plus"],  M["den_pt_" + S + "_plus"],  kPtNominal);
    const Incl iMM = InclusiveAbove(M["num_pt_" + S + "_minus"], M["den_pt_" + S + "_minus"], kPtNominal);
    const Incl iWp = InclusiveAbove(oWp.H1("num_pt_" + S), oWp.H1("den_pt_" + S), kPtNominal);
    const Incl iWm = InclusiveAbove(oWm.H1("num_pt_" + S), oWm.H1("den_pt_" + S), kPtNominal);
    const Incl iBitD = InclusiveAbove(oData.H1("bit_pt_" + S), oData.H1("den_pt_" + S), kPtNominal);
    const Incl iBitM = InclusiveAbove(M["bit_pt_" + S], M["den_pt_" + S], kPtNominal);
    // gen-weighted MC (normal-approx binomial error)
    double wDen = 0, wNum = 0;
    {
      TH1D *hd = M["den_pt_" + S + "_w"], *hn = M["num_pt_" + S + "_w"];
      const int b0 = hd->GetXaxis()->FindBin(kPtNominal + 1e-6);
      for (int i = b0; i <= hd->GetNbinsX(); ++i) { wDen += hd->GetBinContent(i); wNum += hn->GetBinContent(i); }
    }
    // k_s-weighted Wp/Wm combination (physical charge mix)
    const double kWpS = pONorm::MCScale(isMu ? "Wp_mu" : "Wp_ele");
    const double kWmS = pONorm::MCScale(isMu ? "Wm_mu" : "Wm_ele");
    const double ksDen = kWpS * iWp.den + kWmS * iWm.den, ksNum = kWpS * iWp.num + kWmS * iWm.num;
    const double sfIncl = (iM.eff > 0) ? iD.eff / iM.eff : 0;
    const double sfErr  = (iM.eff > 0 && iD.eff > 0) ? sfIncl * std::sqrt(std::pow(0.5 * (iD.lo + iD.hi) / iD.eff, 2) + std::pow(0.5 * (iM.lo + iM.hi) / iM.eff, 2)) : 0;
    // MB-bit fractions from the stored histograms (so plotsOnly reproduces them):
    // of ALL selected events (den/all) and of the TRIGGERED ones (num/trg).
    const Incl fSelD = InclusiveAbove(oData.H1("den_pt_" + S), oData.H1("all_pt_" + S), kPtNominal);
    const Incl fTrgD = InclusiveAbove(oData.H1("num_pt_" + S), oData.H1("trg_pt_" + S), kPtNominal);
    const Incl fSelM = InclusiveAbove(M["den_pt_" + S], M["all_pt_" + S], kPtNominal);
    const double mbFracSel = fSelD.eff, mbFracTrg = fTrgD.eff;

    std::cout << "\n[INCL] " << S << " pT > " << kPtNominal << ":\n";
    std::cout << Form("  data : den %.0f  num %.0f  eff = %.4f -%.4f +%.4f   (bit only %.4f;  W+ %.4f -%.4f +%.4f [%.0f]  W- %.4f -%.4f +%.4f [%.0f])\n",
                      iD.den, iD.num, iD.eff, iD.lo, iD.hi, iBitD.eff, iDP.eff, iDP.lo, iDP.hi, iDP.den, iDM.eff, iDM.lo, iDM.hi, iDM.den);
    std::cout << Form("  MC   : den %.0f  num %.0f  eff = %.4f -%.4f +%.4f   (bit only %.4f;  W+ %.4f -%.4f +%.4f [%.0f]  W- %.4f -%.4f +%.4f [%.0f])\n",
                      iM.den, iM.num, iM.eff, iM.lo, iM.hi, iBitM.eff, iMP.eff, iMP.lo, iMP.hi, iMP.den, iMM.eff, iMM.lo, iMM.hi, iMM.den);
    std::cout << Form("  MC   : per sample Wp %.4f (den %.0f)  Wm %.4f (den %.0f);  gen-weighted %.4f;  k_s-weighted Wp+Wm %.4f\n",
                      iWp.eff, iWp.den, iWm.eff, iWm.den, wDen > 0 ? wNum / wDen : 0, ksDen > 0 ? ksNum / ksDen : 0);
    std::cout << Form("  data MB fraction of the selected events: all %.4f (%.0f/%.0f)   triggered %.4f (%.0f/%.0f)   [prescale check: should agree];  MC %.4f\n",
                      mbFracSel, fSelD.num, fSelD.den, mbFracTrg, fTrgD.num, fTrgD.den, fSelM.eff);
    std::cout << Form("[RESULT] %s %s pT>%.0f: eff_data = %.4f -%.4f +%.4f  eff_MC = %.4f -%.4f +%.4f  SF = %.4f +- %.4f\n",
                      ch.c_str(), S.c_str(), kPtNominal, iD.eff, iD.lo, iD.hi, iM.eff, iM.lo, iM.hi, sfIncl, sfErr);

    // ---- write objects ----
    fout->cd();
    for (auto &p : M) p.second->Write("", TObject::kOverwrite);
    mDen2->Write("", TObject::kOverwrite); mNum2->Write("", TObject::kOverwrite);
    for (TEfficiency *e : {eD_pt, eM_pt, eWp_pt, eWm_pt, eD_ptP, eD_ptM, eM_ptP, eM_ptM, eD_y, eM_y, eD_yP, eD_yM, eM_yP, eM_yM,
                           eD_bit, eM_bit, fD_all, fD_trg, fM_all, fM_trg, e2D_D, e2D_M})
      if (e) e->Write("", TObject::kOverwrite);
    for (TGraphAsymmErrors *g : {sf_pt, sf_y, sf_yP, sf_yM})
      if (g) g->Write("", TObject::kOverwrite);

    // ---- plots ----
    const std::string hdr  = procLbl + ", " + lepSym + " trigger turn-on";
    const std::string sub1 = std::string(kSelLabel[s]) + ", |#eta| < 2.4";
    const std::string sub2 = "den: MB-triggered;  num: + " + pathLbl + ", #DeltaR < 0.4 match";
    const std::vector<std::string> box = {
        Form("p_{T} > %.0f: #varepsilon_{data} = %.3f_{-%.3f}^{+%.3f}", kPtNominal, iD.eff, iD.lo, iD.hi),
        Form("#varepsilon_{MC} = %.3f_{-%.3f}^{+%.3f}", iM.eff, iM.lo, iM.hi),
        Form("SF = %.3f #pm %.3f", sfIncl, sfErr),
        Form("N_{den}: data %.0f, MC %.0f", iD.den, iM.den)};

    // full turn-on (from the loop floor) + zoom on the analysis range
    DrawEffRatio({{eD_pt, "Data (MB-triggered)", kBlack, 20}, {eM_pt, "W signal MC", kRed + 1, 24}},
                 {{0, 1, "Data / MC", kBlack, 20}},
                 outDir + "/turnon_pt_" + S, "leading-lepton p_{T} [GeV]", "Data / MC",
                 hdr, sub1, sub2, box, kPtFloor, 120.0, 0.0, 1.6, 0.7, 1.3, kPtNominal);
    DrawEffRatio({{eD_pt, "Data (MB-triggered)", kBlack, 20}, {eM_pt, "W signal MC", kRed + 1, 24}},
                 {{0, 1, "Data / MC", kBlack, 20}},
                 outDir + "/turnon_pt_" + S + "_zoom25", "leading-lepton p_{T} [GeV]", "Data / MC",
                 hdr, sub1, sub2, box, kPtNominal, 120.0, 0.6, 1.25, 0.85, 1.15);
    // per charge
    DrawEffRatio({{eD_ptP, "Data " + lepSym + "^{+}", kBlack, 20}, {eD_ptM, "Data " + lepSym + "^{-}", kBlue + 1, 21},
                  {eM_ptP, "MC " + lepSym + "^{+}", kRed + 1, 24}, {eM_ptM, "MC " + lepSym + "^{-}", kOrange + 7, 25}},
                 {{0, 2, "Data / MC, " + lepSym + "^{+}", kBlack, 20}, {1, 3, "Data / MC, " + lepSym + "^{-}", kBlue + 1, 21}},
                 outDir + "/eff_pt_" + S + "_charge", "leading-lepton p_{T} [GeV]", "Data / MC",
                 hdr, sub1, sub2, {}, kPtNominal, 120.0, 0.6, 1.25, 0.85, 1.15, kNoLine, "Trigger efficiency", 0.16);
    // vs y (pT > 25)
    DrawEffRatio({{eD_y, "Data (MB-triggered)", kBlack, 20}, {eM_y, "W signal MC", kRed + 1, 24}},
                 {{0, 1, "Data / MC", kBlack, 20}},
                 outDir + "/eff_y_" + S, "y = -#eta_{lab}", "Data / MC",
                 hdr, sub1 + ", p_{T} > 25 GeV", sub2, box, kYEdges[0], kYEdges[kNY], 0.6, 1.25, 0.85, 1.15);
    DrawEffRatio({{eD_yP, "Data " + lepSym + "^{+}", kBlack, 20}, {eD_yM, "Data " + lepSym + "^{-}", kBlue + 1, 21},
                  {eM_yP, "MC " + lepSym + "^{+}", kRed + 1, 24}, {eM_yM, "MC " + lepSym + "^{-}", kOrange + 7, 25}},
                 {{0, 2, "Data / MC, " + lepSym + "^{+}", kBlack, 20}, {1, 3, "Data / MC, " + lepSym + "^{-}", kBlue + 1, 21}},
                 outDir + "/eff_y_" + S + "_charge", "y = -#eta_{lab}", "Data / MC",
                 hdr, sub1 + ", p_{T} > 25 GeV", sub2, {}, kYEdges[0], kYEdges[kNY], 0.6, 1.25, 0.85, 1.15, kNoLine, "Trigger efficiency", 0.16);
    // HLT bit vs bit+match
    DrawEffRatio({{eD_bit, "Data: path fired", kGray + 2, 20}, {eD_pt, "Data: path fired + matched", kBlack, 21},
                  {eM_bit, "MC: path fired", kRed - 7, 24}, {eM_pt, "MC: path fired + matched", kRed + 1, 25}},
                 {{1, 0, "Data: matched / fired", kBlack, 20}, {3, 2, "MC: matched / fired", kRed + 1, 24}},
                 outDir + "/bit_vs_match_pt_" + S, "leading-lepton p_{T} [GeV]", "matched / fired",
                 hdr, sub1, "den: MB-triggered", {}, kPtFloor, 120.0, 0.0, 1.6, 0.8, 1.1, kPtNominal, "Trigger efficiency", 0.16);
    // MB-fraction control (prescale independence): den/all and num/trg, data only
    // (MC fires the MB path on 99.9% of W events -- quoted in the box).
    DrawEffRatio({{fD_all, "all selected events", kBlack, 20}, {fD_trg, "triggered events only", kBlue + 1, 21}},
                 {},
                 outDir + "/mbfrac_pt_" + S, "leading-lepton p_{T} [GeV]", "",
                 procLbl + ", MB-tag control (data)", sub1, "MB-bit fraction of the selected events (a pure prescale: flat, identical)",
                 {Form("p_{T} > 25: all %.4f, triggered %.4f", mbFracSel, mbFracTrg),
                  Form("MC: MB fires on %.2f%% of W events", 100.0 * fSelM.eff)},
                 kPtFloor, 120.0, 0.0, 1.6, 0.8, 1.2, kPtNominal, "MB-bit fraction");
  }

  // ---- data per-run table (nom, pT > 25) ----
  std::cout << "\n[RUNS] data, nom selection, pT > 25:  run | selected | MB (den) | num | eff | MB frac | triggered-no-MB\n";
  for (auto &p : oData.byRun)
  {
    const auto &v = p.second;
    std::cout << Form("   %d | %8.0f | %8.0f | %8.0f | %s | %.3f | %.0f\n", p.first, v[0], v[1], v[2],
                      v[1] > 0 ? Form("%.4f", v[2] / v[1]) : "  --  ", v[0] > 0 ? v[1] / v[0] : 0.0, v[3]);
  }
  std::cout << "[OK] plots in " << outDir << "/\n";
}

} // anonymous namespace

// ============================================================
// Entry point
// ============================================================
void trig_eff_mb(const char *channel = "both", bool plotsOnly = false)
{
  gSystem->mkdir("rootfile", kTRUE);
  gSystem->mkdir("plots", kTRUE);
  gSystem->mkdir("logs", kTRUE);

  const std::string chArg = channel;
  std::vector<bool> flavours;
  if (chArg == "both") flavours = {true, false};
  else if (chArg == "mu") flavours = {true};
  else if (chArg == "ele") flavours = {false};
  else { std::cerr << "[ERR] channel must be mu | ele | both\n"; return; }

  for (bool isMu : flavours)
  {
    const std::string ch = isMu ? "mu" : "ele";
    const std::string rootName = "rootfile/trig_eff_mb_" + ch + ".root";

    TrigOut oData, oWp, oWm;
    oData.Book("data_" + ch);
    oWp.Book("Wp_" + ch);
    oWm.Book("Wm_" + ch);

    if (plotsOnly)
    {
      TFile *fin = TFile::Open(rootName.c_str());
      if (!fin || fin->IsZombie()) { std::cerr << "[ERR] plotsOnly: cannot open " << rootName << "\n"; continue; }
      for (TrigOut *o : {&oData, &oWp, &oWm})
      {
        for (auto &p : o->h1)
        {
          TH1D *h = (TH1D *)fin->Get(p.second->GetName());
          if (!h) { std::cerr << "[ERR] plotsOnly: missing " << p.second->GetName() << "\n"; continue; }
          p.second->Reset(); p.second->Add(h);
        }
        for (auto &p : o->h2)
        {
          TH2D *h = (TH2D *)fin->Get(p.second->GetName());
          if (!h) { std::cerr << "[ERR] plotsOnly: missing " << p.second->GetName() << "\n"; continue; }
          p.second->Reset(); p.second->Add(h);
        }
      }
      fin->Close();
      // NB the pT>25 counters and the per-run table are not stored; they are only in the full-run log.
      TFile *fout = TFile::Open(rootName.c_str(), "UPDATE");
      Report(isMu, oData, oWp, oWm, fout);
      fout->Close();
      continue;
    }

    if (RunTrig(isMu, kData, oData) != 0) { std::cerr << "[ERR] data run failed (" << ch << ")\n"; continue; }
    if (RunTrig(isMu, kWp, oWp) != 0)     { std::cerr << "[ERR] Wp run failed (" << ch << ")\n"; continue; }
    if (RunTrig(isMu, kWm, oWm) != 0)     { std::cerr << "[ERR] Wm run failed (" << ch << ")\n"; continue; }

    TFile *fout = TFile::Open(rootName.c_str(), "RECREATE");
    oData.Write(fout); oWp.Write(fout); oWm.Write(fout);
    Report(isMu, oData, oWp, oWm, fout);
    fout->Close();
    std::cout << "[OK] wrote " << rootName << "\n";
  }
}
