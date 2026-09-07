// correction/charge_flip.C
//
// Lepton CHARGE-MISIDENTIFICATION (charge-flip) rate of the selected W lepton,
// measured on the W SIGNAL MC (Wp + Wm, per flavor), for muons and electrons.
// A flip moves a W+ event into the reco W- region (and vice versa), i.e. it
// dilutes the charge asymmetry directly; this macro measures how often that
// happens for the exact object the analysis uses.
//
// Method
//   1. Replicate the skim's full 8-step W selection (skim/skim.C, nominal
//      cuts; the same replication as correction/njet_WZ.C, which is verified
//      event-identical with the skim) on the W MC -> the leading selected
//      lepton = the object whose reco charge decides the W+/W- region.
//   2. Match it to the generator lepton CHARGE-BLIND: gen |pdg| = 13 (11),
//      DeltaR(gen, reco) < 0.5 and |pT(gen) - pT(reco)| / pT(gen) < 0.5 --
//      the AN's gen-reco criteria (= skim_common.h::
//      PassGenRecoMatchingWithAncestor) WITHOUT its same-charge requirement.
//      Same-charge matching cannot see a flip: the flipped lepton just fails
//      to match. Among candidates, W-ancestor ones would be preferred, ties
//      by smallest DeltaR -- but the filtered gen collection (pT > 5,
//      |eta| < 2.5) stores NO |pdg| = 24 entry at all (measured 2026-09-01:
//      0 in 100k events), so in practice the DR-nearest candidate is always
//      used; the log reports the count for transparency (every stored W
//      lepton has nMothers = 1 with motherIdx = -999, i.e. "mother not in
//      the stored list"). NB the ntuplizer's own mu_/ele_genMatchedIndex is
//      an index into the EventTree's OWN gen block (nMC/mcPID/mcStatus/mcPt/
//      mcMomPID/..., which DOES store the W as a status-62 entry with direct
//      parentage), NOT into HiGenParticleAna/hi -- read against `hi` it
//      lands on the neutrino. Against mc* it is valid (99.96% of pT>25 reco
//      muons -> a status-1 muon with mcMomPID = 24) but CHARGE-AWARE: 0
//      charge disagreements in 45k matched electrons, so a flipped lepton
//      comes back as idx = -1 -- the same blindness as the same-charge
//      helper, which is why it is not used here either.
//   3. Classify the selected lepton:
//        CORRECT    matched, reco charge == gen charge
//        FLIP       matched, reco charge != gen charge
//        UNMATCHED  no candidate (fakes / non-prompt / hard-brem electrons
//                   outside the pT window) -- reported, NOT in the denominator
//      f = N_flip / (N_correct + N_flip), Clopper-Pearson 68.3% intervals
//      (TEfficiency) on RAW counts: the gen weight is ~constant (= sigma, with
//      ~0.9% negative entries), so raw and weighted rates agree to ~2%; the
//      weighted rate is printed next to it.
//   4. The "so what": the asymmetry bias. On the same matched leptons the MC
//      gives the charge asymmetry once with the RECO charge and once with the
//      GEN charge (Wp/Wm combined with the absolute k_s from skim/mc_norm.h),
//      per analysis y bin: DeltaA = A_reco - A_gen is the flip-induced bias
//      the measurement would carry uncorrected (= -2 f A to first order for a
//      charge-symmetric f).
//
// Binning: y = -eta_lab in the analysis bins (pOSkim::kYEdges; the skim's
// p-going = forward convention, so bin i here IS analysis bin i), and lepton
// pT in {25,30,35,40,45,50,60,80,120} (pT > 120 folded into the last bin).
//
// Electron extras: the ntuple carries eleCharge (the analysis charge) AND
// eleTrkCharge; the rate is also given (a) using eleTrkCharge as the reco
// charge, (b) on the CHARGE-CONSISTENT subset eleCharge == eleTrkCharge
// (the standard CMS handle against flips -- shows what a consistency
// requirement would buy), plus the inconsistent fraction itself.
//
// Caveats: the MC is UNEMBEDDED (no pO underlying event), so flips from track
// confusion in the UE are not represented (small at pO multiplicities), and
// no lepton SFs / momentum corrections are applied (irrelevant for a rate).
// The data-driven cross-check is the same-sign / opposite-sign Z->ee ratio
// (Z->mumu has ~1 SS pair in the whole dataset).
//
// Outputs (run from correction/):
//   rootfile/charge_flip_<mu|ele>.root   raw count histos per sample + combined,
//                                        TEfficiency objects, labeled counters
//   plots/charge_flip_<mu|ele>/          fliprate_vs_y, fliprate_vs_pt,
//                                        fliprate_2D_pt_y, unmatched_vs_y,
//                                        match_dR, match_dpt (+ electron extras
//                                        fliprate_vs_y_chargedef, incons_vs_y)
//   stdout                               the tables -- run through
//                                        ./run_charge_flip.sh to keep the log
//
// Run:
//   ./run_charge_flip.sh [mu|ele|both]            (keeps logs/charge_flip_<chan>.log)
//   root -l -b -q 'charge_flip.C+("mu")'          (bare; ~2 min per MC file)
//   root -l -b -q 'charge_flip.C+("both", true)'  (tables + plots only, from
//                                                  the stored rootfile)

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TInterpreter.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"

#include <cmath>
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

// -------- gen-reco matching window (the AN criteria, charge-blind) --------
const double kMatchDR  = 0.5;
const double kMatchDPt = 0.5;

// -------- binning --------
const int    kNPt = 8;
const double kPtEdges[kNPt + 1] = {25, 30, 35, 40, 45, 50, 60, 80, 120};
const double kPtFoldMax = 119.9; // pT above the last edge is folded into the last bin

// log10(DR) storage for the match-quality shapes: [-4, -0.3] (DR 1e-4 .. 0.5)
const int    kNLogDR  = 74;
const double kLogDRLo = -4.0;
const double kLogDRHi = -0.3;
const double kDRFloor = 1e-4; // DR below this is stored at the floor

// Clopper-Pearson level for every quoted interval.
const double kCL = 0.6827;

// ============================================================
// Per-sample bookkeeping
// ============================================================
struct Counters
{
  double nSel = 0, nMatched = 0, nFlip = 0, nUnm = 0, nNoWAnc = 0, nGenChgNe = 0;
  double nIncons = 0, nInconsFlip = 0; // electron: eleCharge != eleTrkCharge (matched)
  double wSel = 0, wMatched = 0, wFlip = 0; // gen-weighted
};

const char *kCounterLabels[] = {"nSel", "nMatched", "nFlip", "nUnm", "nNoWAnc",
                                "nGenChgNe", "nIncons", "nInconsFlip",
                                "wSel", "wMatched", "wFlip"};
const int kNCounters = 11;

struct FlipOut
{
  std::string tag;   // "mu_Wp", "ele_Wm", "mu_all", ...
  bool isMu = true;
  Counters c;
  std::map<std::string, TH1 *> H; // short name -> histogram (owned, detached)

  TH1D *H1(const std::string &k) const
  {
    auto it = H.find(k);
    if (it == H.end()) { std::cerr << "[BUG] FlipOut::H1 missing " << k << "\n"; return nullptr; }
    return (TH1D *)it->second;
  }
  TH2D *H2(const std::string &k) const
  {
    auto it = H.find(k);
    if (it == H.end()) { std::cerr << "[BUG] FlipOut::H2 missing " << k << "\n"; return nullptr; }
    return (TH2D *)it->second;
  }

  // Histogram key list (electron extras only when !isMu).
  std::vector<std::string> Keys1D() const
  {
    std::vector<std::string> k = {
        "sel_y", "tot_y", "flip_y", "unm_y",
        "sel_pt", "tot_pt", "flip_pt", "unm_pt",
        "dR_ok", "dR_flip", "dpt_ok", "dpt_flip",
        "w_recoP_y", "w_recoM_y", "w_genP_y", "w_genM_y"};
    if (!isMu)
      for (const char *e : {"tot_y_trk", "flip_y_trk", "tot_y_cons", "flip_y_cons", "incons_y",
                            "tot_pt_trk", "flip_pt_trk", "tot_pt_cons", "flip_pt_cons", "incons_pt"})
        k.push_back(e);
    return k;
  }
  std::vector<std::string> Keys2D() const { return {"tot2", "flip2"}; }

  void Book(const std::string &t, bool mu)
  {
    tag  = t;
    isMu = mu;
    const char *yT  = "y^{l}_{lab} = -#eta^{l}_{lab}";
    const char *ptT = "p_{T}^{l} [GeV]";
    auto mk1y = [&](const std::string &k, const char *title)
    {
      TH1D *h = new TH1D(Form("h_%s_%s", k.c_str(), tag.c_str()),
                         Form("%s;%s;Leptons", title, yT), kNY, kYEdges);
      h->SetDirectory(nullptr);
      H[k] = h;
    };
    auto mk1pt = [&](const std::string &k, const char *title)
    {
      TH1D *h = new TH1D(Form("h_%s_%s", k.c_str(), tag.c_str()),
                         Form("%s;%s;Leptons", title, ptT), kNPt, kPtEdges);
      h->SetDirectory(nullptr);
      H[k] = h;
    };
    for (const auto &k : Keys1D())
    {
      if (k.find("_pt") != std::string::npos)      mk1pt(k, k.c_str());
      else if (k.find("_y") != std::string::npos)  mk1y(k, k.c_str());
    }
    // match quality (raw counts). DR spans decades (mu ~3e-4, e brem tails),
    // so it is stored as log10(DR) on a 0.05 grid: -2 (DR = 0.01) and -1
    // (DR = 0.1) are bin edges for the fractions quoted in the log.
    H["dR_ok"]    = new TH1D(Form("h_dR_ok_%s", tag.c_str()),    ";log_{10} #DeltaR(gen, reco);Leptons", kNLogDR, kLogDRLo, kLogDRHi);
    H["dR_flip"]  = new TH1D(Form("h_dR_flip_%s", tag.c_str()),  ";log_{10} #DeltaR(gen, reco);Leptons", kNLogDR, kLogDRLo, kLogDRHi);
    H["dpt_ok"]   = new TH1D(Form("h_dpt_ok_%s", tag.c_str()),   ";(p_{T}^{reco} - p_{T}^{gen}) / p_{T}^{gen};Leptons", 100, -kMatchDPt, kMatchDPt);
    H["dpt_flip"] = new TH1D(Form("h_dpt_flip_%s", tag.c_str()), ";(p_{T}^{reco} - p_{T}^{gen}) / p_{T}^{gen};Leptons", 100, -kMatchDPt, kMatchDPt);
    for (const char *k : {"dR_ok", "dR_flip", "dpt_ok", "dpt_flip"}) H[k]->SetDirectory(nullptr);
    // 2D (pT x y) raw counts
    H["tot2"]  = new TH2D(Form("h2_tot_pt_y_%s", tag.c_str()),  Form(";%s;%s", ptT, yT), kNPt, kPtEdges, kNY, kYEdges);
    H["flip2"] = new TH2D(Form("h2_flip_pt_y_%s", tag.c_str()), Form(";%s;%s", ptT, yT), kNPt, kPtEdges, kNY, kYEdges);
    H["tot2"]->SetDirectory(nullptr);
    H["flip2"]->SetDirectory(nullptr);
  }

  // Clone-and-add combination (P + M), raw and weighted alike.
  void CombineFrom(const FlipOut &a, const FlipOut &b, const std::string &t)
  {
    tag  = t;
    isMu = a.isMu;
    for (const auto &kv : a.H)
    {
      TH1 *h = (TH1 *)kv.second->Clone(Form("%s_%s", kv.second->GetName(), "all"));
      h->SetDirectory(nullptr);
      // rename to the canonical <what>_<tag> name
      std::string nm = kv.second->GetName();
      const size_t p = nm.rfind("_" + a.tag);
      if (p != std::string::npos) nm = nm.substr(0, p) + "_" + tag;
      h->SetName(nm.c_str());
      auto it = b.H.find(kv.first);
      if (it != b.H.end()) h->Add(it->second);
      H[kv.first] = h;
    }
    c.nSel = a.c.nSel + b.c.nSel;             c.nMatched = a.c.nMatched + b.c.nMatched;
    c.nFlip = a.c.nFlip + b.c.nFlip;          c.nUnm = a.c.nUnm + b.c.nUnm;
    c.nNoWAnc = a.c.nNoWAnc + b.c.nNoWAnc;    c.nGenChgNe = a.c.nGenChgNe + b.c.nGenChgNe;
    c.nIncons = a.c.nIncons + b.c.nIncons;    c.nInconsFlip = a.c.nInconsFlip + b.c.nInconsFlip;
    c.wSel = a.c.wSel + b.c.wSel;             c.wMatched = a.c.wMatched + b.c.wMatched;
    c.wFlip = a.c.wFlip + b.c.wFlip;
  }

  TH1D *CountersHist() const
  {
    TH1D *h = new TH1D(Form("h_counters_%s", tag.c_str()), "labeled counters", kNCounters, 0, kNCounters);
    h->SetDirectory(nullptr);
    const double v[kNCounters] = {c.nSel, c.nMatched, c.nFlip, c.nUnm, c.nNoWAnc, c.nGenChgNe,
                                  c.nIncons, c.nInconsFlip, c.wSel, c.wMatched, c.wFlip};
    for (int i = 0; i < kNCounters; ++i)
    {
      h->GetXaxis()->SetBinLabel(i + 1, kCounterLabels[i]);
      h->SetBinContent(i + 1, v[i]);
    }
    return h;
  }

  void Write(TFile *f) const
  {
    f->cd();
    for (const auto &kv : H) kv.second->Write("", TObject::kOverwrite);
    TH1D *hc = CountersHist();
    hc->Write("", TObject::kOverwrite);
    delete hc;
  }

  // Read back everything written by Write() (plots-only mode). Labels are
  // resolved by looping GetBinLabel -- TAxis::FindBin(label) would APPEND on a
  // miss (repo rule).
  bool Read(TFile *f, const std::string &t, bool mu)
  {
    tag  = t;
    isMu = mu;
    for (const auto &k : Keys1D())
    {
      TH1 *h = (TH1 *)f->Get(Form("h_%s_%s", k.c_str(), tag.c_str()));
      if (!h) { std::cerr << "[ERR] missing h_" << k << "_" << tag << " in " << f->GetName() << "\n"; return false; }
      TH1 *cl = (TH1 *)h->Clone(Form("%s_rd", h->GetName()));
      cl->SetDirectory(nullptr);
      cl->SetName(h->GetName());
      H[k] = cl;
    }
    const char *n2[2] = {"h2_tot_pt_y", "h2_flip_pt_y"};
    const char *k2[2] = {"tot2", "flip2"};
    for (int i = 0; i < 2; ++i)
    {
      TH1 *h = (TH1 *)f->Get(Form("%s_%s", n2[i], tag.c_str()));
      if (!h) { std::cerr << "[ERR] missing " << n2[i] << "_" << tag << "\n"; return false; }
      TH1 *cl = (TH1 *)h->Clone(Form("%s_rd", h->GetName()));
      cl->SetDirectory(nullptr);
      cl->SetName(h->GetName());
      H[k2[i]] = cl;
    }
    TH1D *hc = (TH1D *)f->Get(Form("h_counters_%s", tag.c_str()));
    if (!hc) { std::cerr << "[ERR] missing h_counters_" << tag << "\n"; return false; }
    double *v[kNCounters] = {&c.nSel, &c.nMatched, &c.nFlip, &c.nUnm, &c.nNoWAnc, &c.nGenChgNe,
                             &c.nIncons, &c.nInconsFlip, &c.wSel, &c.wMatched, &c.wFlip};
    for (int i = 0; i < kNCounters; ++i)
    {
      bool found = false;
      for (int b = 1; b <= hc->GetNbinsX(); ++b)
        if (std::string(hc->GetXaxis()->GetBinLabel(b)) == kCounterLabels[i])
        { *v[i] = hc->GetBinContent(b); found = true; break; }
      if (!found) { std::cerr << "[ERR] counter " << kCounterLabels[i] << " missing\n"; return false; }
    }
    return true;
  }
};

// ============================================================
// Charge-blind gen match: best gen lepton of |pdg| = flavPdg inside the
// (DR, dpT/pT) window; W-ancestor candidates preferred, then smallest DR.
// Returns the gen index or -1.
// ============================================================
int MatchGenChargeBlind(double pt, double eta, double phi, int flavPdg,
                        std::vector<float> *gPt, std::vector<float> *gEta,
                        std::vector<float> *gPhi, std::vector<int> *gPdg,
                        std::vector<std::vector<int>> *gMom,
                        double &dR, double &dPtRel, bool &hasWAnc)
{
  dR = 1e9; dPtRel = 0; hasWAnc = false;
  if (!gPt || !gEta || !gPhi || !gPdg || !gMom) return -1;
  const size_t ng = std::min({gPt->size(), gEta->size(), gPhi->size(), gPdg->size(), gMom->size()});
  int best = -1, bestAnc = -1;
  double drBest = 1e9, drAnc = 1e9;
  for (size_t i = 0; i < ng; ++i)
  {
    if (std::abs(gPdg->at(i)) != flavPdg) continue;
    if (gPt->at(i) <= 0) continue;
    const double d = DeltaR(gEta->at(i), gPhi->at(i), eta, phi);
    if (d >= kMatchDR) continue;
    if (std::fabs(gPt->at(i) - pt) / gPt->at(i) >= kMatchDPt) continue;
    if (d < drBest) { drBest = d; best = (int)i; }
    if (HasAncestor((int)i, 24, gPdg, gMom))
      if (d < drAnc) { drAnc = d; bestAnc = (int)i; }
  }
  const int use = (bestAnc >= 0) ? bestAnc : best;
  if (use < 0) return -1;
  hasWAnc = (bestAnc >= 0);
  dR      = (bestAnc >= 0) ? drAnc : drBest;
  dPtRel  = (pt - gPt->at(use)) / gPt->at(use);
  return use;
}

// ============================================================
// W selection (the skim's 8 steps, nominal cuts -- copied from
// correction/njet_WZ.C::RunW, which is verified event-identical with skim.C)
// + charge-blind gen match + classification. Returns 0 on success.
// ============================================================
int RunFlip(bool isMu, SampleType sample, FlipOut &out)
{
  const double lepPt25     = 25.0;
  const double dyPtMin     = isMu ? 15.0 : 10.0;   // DY-veto leg pT (intentionally asymmetric)
  const double isoMax      = isMu ? 0.15 : 0.095;
  const double dyMassMin   = 80.0, dyMassMax = 110.0;
  const double etaMax      = 2.4;
  const double trigMatchDR = 0.4;  // 2026-08-12: tracks skim.C:283/963
  const double vzMax       = 15.0;
  const double lepMass     = isMu ? MU_MASS : ELE_MASS;
  const char  *tagLep      = isMu ? "muon" : "electron";
  const int    flavPdg     = isMu ? 13 : 11;
  const int    sampleChg   = (sample == kWp) ? +1 : -1;

  const auto info = ResolveMCSample(sample, isMu ? "mu" : "ele");
  if (info.fname.empty()) { std::cerr << "[ERR] RunFlip: cannot resolve MC sample\n"; return 1; }
  const std::string fname = info.fname;

  std::cout << "[INPUT] " << fname << std::endl;
  TFile *f = TFile::Open(fname.c_str());
  if (!f || f->IsZombie()) { std::cerr << "[ERR] cannot open " << fname << "\n"; return 1; }

  TTree *tLep    = (TTree *)f->Get("ggHiNtuplizer/EventTree");
  TTree *tHi     = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  TTree *tHLT    = (TTree *)f->Get("hltanalysis/HltTree");
  TTree *tHLTobj = (TTree *)f->Get(isMu ? "hltobject/HLT_OxyL1SingleMuOpen_v"
                                        : "hltobject/HLT_OxyL1SingleEG10_v");
  TTree *tEvent  = (TTree *)f->Get("skimanalysis/HltTree");
  TTree *tGen    = (TTree *)f->Get("HiGenParticleAna/hi");
  if (!tLep || !tHi || !tHLT || !tHLTobj || !tEvent || !tGen)
  {
    std::cerr << "[FATAL] RunFlip: missing a required tree in " << fname << "\n";
    f->Close();
    return 2;
  }
  if (tGen->GetEntries() != tLep->GetEntries())
  {
    std::cerr << "[FATAL] RunFlip: gen tree (" << tGen->GetEntries() << ") and EventTree ("
              << tLep->GetEntries() << ") entry counts differ\n";
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
  std::vector<int>   *lepCharge = nullptr, *lepID = nullptr, *muIsPF = nullptr, *trkCharge = nullptr;
  std::vector<float> *chIso = nullptr, *neuIso = nullptr, *phoIso = nullptr, *puIso = nullptr;

  tLep->SetBranchStatus("*", 0);
  for (const char *bn : {bN, bPt, bEta, bPhi, bChg, bID, bCh, bNeu, bPho})
    if (!HasBranch(tLep, bn))
    {
      std::cerr << "[FATAL] RunFlip: missing EventTree branch " << bn << "\n";
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
  else std::cout << "[WARN] RunFlip: no " << bPU << " branch; relIso uncorrected (no Delta-beta).\n";

  const bool has_muIsPF = isMu && HasBranch(tLep, "muIsPF");
  if (has_muIsPF) { tLep->SetBranchStatus("muIsPF", 1); tLep->SetBranchAddress("muIsPF", &muIsPF); }

  const bool has_trkChg = !isMu && HasBranch(tLep, "eleTrkCharge");
  if (has_trkChg) { tLep->SetBranchStatus("eleTrkCharge", 1); tLep->SetBranchAddress("eleTrkCharge", &trkCharge); }
  else if (!isMu) std::cout << "[WARN] RunFlip: no eleTrkCharge branch; electron charge-consistency extras skipped.\n";

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
    std::cout << "[WARN] RunFlip: no branch containing " << hltNeedle << "; trigger treated as PASS.\n";

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
    std::cout << "[WARN] RunFlip: could not find HLT objects.\n";

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
  const bool has_genWeight = HasBranch(tHi, "weight");
  if (has_genWeight) { tHi->SetBranchStatus("weight", 1); tHi->SetBranchAddress("weight", &genWeight); }
  else std::cout << "[WARN] RunFlip: no HiTree 'weight'; weighted sums unweighted.\n";

  // -------- gen collection --------
  std::vector<float> *gPt = nullptr, *gEta = nullptr, *gPhi = nullptr;
  std::vector<int>   *gChg = nullptr, *gPdg = nullptr;
  std::vector<std::vector<int>> *gMom = nullptr;
  tGen->SetBranchStatus("*", 0);
  for (const char *bn : {"pt", "eta", "phi", "chg", "pdg", "motherIdx"})
    if (!HasBranch(tGen, bn))
    {
      std::cerr << "[FATAL] RunFlip: missing gen branch " << bn << "\n";
      f->Close();
      return 2;
    }
  tGen->SetBranchStatus("pt", 1);        tGen->SetBranchAddress("pt",        &gPt);
  tGen->SetBranchStatus("eta", 1);       tGen->SetBranchAddress("eta",       &gEta);
  tGen->SetBranchStatus("phi", 1);       tGen->SetBranchAddress("phi",       &gPhi);
  tGen->SetBranchStatus("chg", 1);       tGen->SetBranchAddress("chg",       &gChg);
  tGen->SetBranchStatus("pdg", 1);       tGen->SetBranchAddress("pdg",       &gPdg);
  tGen->SetBranchStatus("motherIdx", 1); tGen->SetBranchAddress("motherIdx", &gMom);

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
    // tGen is read only for events passing the full selection (lockstep entry)

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

    // -------- selected: charge-blind gen match + classification --------
    tGen->GetEntry(ie);
    const double pt  = lepPt->at(iLead);
    const double eta = lepEta->at(iLead);
    const double phi = lepPhi->at(iLead);
    const int    qReco = lepCharge->at(iLead);
    const double y   = -eta;                       // skim convention: p-going (-Z) = forward
    const double ptF = std::min(pt, kPtFoldMax);   // overflow folded into the last pT bin

    out.c.nSel += 1; out.c.wSel += w;
    out.H1("sel_y")->Fill(y);
    out.H1("sel_pt")->Fill(ptF);

    double dR = 0, dPtRel = 0;
    bool hasWAnc = false;
    const int ig = MatchGenChargeBlind(pt, eta, phi, flavPdg, gPt, gEta, gPhi, gPdg, gMom,
                                       dR, dPtRel, hasWAnc);
    if (ig < 0)
    {
      out.c.nUnm += 1;
      out.H1("unm_y")->Fill(y);
      out.H1("unm_pt")->Fill(ptF);
      continue;
    }

    const int qGen = gChg->at(ig);
    out.c.nMatched += 1; out.c.wMatched += w;
    if (!hasWAnc) out.c.nNoWAnc += 1;
    if (qGen != sampleChg) out.c.nGenChgNe += 1;

    out.H1("tot_y")->Fill(y);
    out.H1("tot_pt")->Fill(ptF);
    out.H2("tot2")->Fill(ptF, y);
    out.H1(qReco > 0 ? "w_recoP_y" : "w_recoM_y")->Fill(y, w);
    out.H1(qGen  > 0 ? "w_genP_y"  : "w_genM_y")->Fill(y, w);

    const bool flip = (qReco != qGen);
    const double logDR = std::log10(std::max(dR, kDRFloor));
    if (flip)
    {
      out.c.nFlip += 1; out.c.wFlip += w;
      out.H1("flip_y")->Fill(y);
      out.H1("flip_pt")->Fill(ptF);
      out.H2("flip2")->Fill(ptF, y);
      out.H1("dR_flip")->Fill(logDR);
      out.H1("dpt_flip")->Fill(dPtRel);
    }
    else
    {
      out.H1("dR_ok")->Fill(logDR);
      out.H1("dpt_ok")->Fill(dPtRel);
    }

    if (has_trkChg && trkCharge)
    {
      const int qTrk = trkCharge->at(iLead);
      out.H1("tot_y_trk")->Fill(y);
      out.H1("tot_pt_trk")->Fill(ptF);
      if (qTrk != qGen) { out.H1("flip_y_trk")->Fill(y); out.H1("flip_pt_trk")->Fill(ptF); }
      if (qTrk == qReco)
      {
        out.H1("tot_y_cons")->Fill(y);
        out.H1("tot_pt_cons")->Fill(ptF);
        if (flip) { out.H1("flip_y_cons")->Fill(y); out.H1("flip_pt_cons")->Fill(ptF); }
      }
      else
      {
        out.c.nIncons += 1;
        if (flip) out.c.nInconsFlip += 1;
        out.H1("incons_y")->Fill(y);
        out.H1("incons_pt")->Fill(ptF);
      }
    }
  }

  f->Close();
  delete f;
  std::cout << Form("[INFO] %s: selected %.0f, matched %.0f, flip %.0f, unmatched %.0f\n",
                    out.tag.c_str(), out.c.nSel, out.c.nMatched, out.c.nFlip, out.c.nUnm);
  return 0;
}

// ============================================================
// Reporting helpers
// ============================================================
struct CP
{
  double f, lo, hi; // rate and its CP-68 interval half-widths (down, up)
};
CP Rate(double pass, double tot)
{
  CP r{0, 0, 0};
  if (tot <= 0) return r;
  r.f  = pass / tot;
  const int t = (int)std::lround(tot), p = (int)std::lround(pass);
  r.lo = r.f - TEfficiency::ClopperPearson(t, p, kCL, false);
  r.hi = TEfficiency::ClopperPearson(t, p, kCL, true) - r.f;
  return r;
}
std::string FmtRate(const CP &r) { return Form("%.3e -%.1e +%.1e", r.f, r.lo, r.hi); }

TEfficiency *MakeEff(TH1 *pass, TH1 *tot, const std::string &name, const char *title)
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

// Draw several TEfficiency's as points with CP error bars on one frame.
void DrawRates(const std::vector<TEfficiency *> &effs,
               const std::vector<std::string> &labels,
               const std::vector<int> &colors,
               const std::vector<int> &markers,
               const std::string &outPath,
               const std::string &xTitle, const std::string &yTitle,
               const std::string &hdr, const std::string &sub1, const std::string &sub2,
               const std::vector<std::string> &box,
               double xlo, double xhi, bool logy)
{
  // Layout: CMS_lumi(…, 10) paints "CMS / Work in Progress" at the top LEFT
  // inside the frame, so the header block sits top RIGHT, the legend under it,
  // and the info box on the left below the CMS label. Generous y headroom
  // keeps the points in the lower ~60% of the frame.
  PlotStyle ps;
  ps.logy = logy;
  ps.headerX = 0.47; ps.headerY = 0.89; ps.headerDy = 0.05;
  ps.titleSize = 0.040; ps.subSize = 0.031;
  ps.boxX1 = 0.17; ps.boxX2 = 0.47; ps.boxY1 = 0.64; ps.boxY2 = 0.76; ps.boxTextSize = 0.029;

  std::vector<TGraphAsymmErrors *> gs;
  double ymax = 0, ymin = 1e30;
  for (auto *e : effs)
  {
    TGraphAsymmErrors *g = e ? e->CreateGraph() : nullptr;
    gs.push_back(g);
    if (!g) continue;
    for (int i = 0; i < g->GetN(); ++i)
    {
      const double yv = g->GetY()[i];
      ymax = std::max(ymax, yv + g->GetEYhigh()[i]);
      if (yv > 0) ymin = std::min(ymin, yv - 0.9 * g->GetEYlow()[i] > 0 ? yv - 0.9 * g->GetEYlow()[i] : yv);
    }
  }
  if (ymax <= 0) ymax = 1e-3;
  if (ymin >= ymax) ymin = ymax / 100.0;

  static int nFrame = 0;
  TCanvas *c = new TCanvas(Form("c_rate_%d", nFrame), "", ps.w, ps.h);
  ApplyCanvasStyle(c, ps);
  TH1F *frame = new TH1F(Form("frame_rate_%d", nFrame++), "", 100, xlo, xhi);
  frame->SetDirectory(nullptr);
  ApplyHistStyle(frame, ps, xTitle, yTitle);
  frame->SetStats(0);
  if (logy) { frame->SetMinimum(ymin / 3.0); frame->SetMaximum(ymax * 200.0); }
  else      { frame->SetMinimum(0.0);        frame->SetMaximum(ymax * 2.4); }
  frame->Draw("AXIS");

  TLegend *leg = new TLegend(0.50, 0.60, 0.93, 0.74);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(ps.font);
  leg->SetTextSize(0.030);
  for (size_t i = 0; i < gs.size(); ++i)
  {
    if (!gs[i]) continue;
    gs[i]->SetMarkerStyle(markers[i]);
    gs[i]->SetMarkerSize(1.3);
    gs[i]->SetMarkerColor(colors[i]);
    gs[i]->SetLineColor(colors[i]);
    gs[i]->SetLineWidth(2);
    gs[i]->Draw("P SAME");
    leg->AddEntry(gs[i], labels[i].c_str(), "lp");
  }
  leg->Draw();
  DrawHeader(ps, hdr, sub1, sub2);
  DrawInfoBox(ps, box);
  CMS_lumi(c, 13, 10);
  c->RedrawAxis();
  c->SaveAs((outPath + ".png").c_str());
  c->SaveAs((outPath + ".pdf").c_str());
  delete c;
}

// Two unit-normalized shapes (correct vs flipped) on one log frame.
void DrawTwoShapes(TH1 *a, TH1 *b, const std::string &la, const std::string &lb,
                   const std::string &outPath, const std::string &xTitle,
                   const std::string &hdr, const std::string &sub1, const std::string &sub2)
{
  if (!a || !b) return;
  PlotStyle ps;
  ps.logy = true;
  ps.headerX = 0.47; ps.headerY = 0.89; ps.headerDy = 0.05; // right of the CMS label
  ps.titleSize = 0.040; ps.subSize = 0.031;

  TH1 *ha = (TH1 *)a->Clone(Form("%s_shape", a->GetName())); ha->SetDirectory(nullptr);
  TH1 *hb = (TH1 *)b->Clone(Form("%s_shape", b->GetName())); hb->SetDirectory(nullptr);
  if (ha->Integral() > 0) ha->Scale(1.0 / ha->Integral());
  if (hb->Integral() > 0) hb->Scale(1.0 / hb->Integral());

  static int n = 0;
  TCanvas *c = new TCanvas(Form("c_shape_%d", n++), "", ps.w, ps.h);
  ApplyCanvasStyle(c, ps);
  ApplyHistStyle(ha, ps, xTitle, "Fraction of leptons");
  ha->SetStats(0);
  const double mx = std::max(ha->GetMaximum(), hb->GetMaximum());
  ha->SetMaximum(mx * 300.0);
  ha->SetMinimum(1e-5);
  ha->SetLineColor(kBlue + 1); ha->SetLineWidth(2); ha->SetMarkerColor(kBlue + 1); ha->SetMarkerStyle(20);
  hb->SetLineColor(kRed + 1);  hb->SetLineWidth(2); hb->SetMarkerColor(kRed + 1);  hb->SetMarkerStyle(24);
  ha->Draw("HIST E");
  hb->Draw("HIST E SAME");
  TLegend *leg = new TLegend(0.50, 0.62, 0.93, 0.74);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(ps.font); leg->SetTextSize(0.030);
  leg->AddEntry(ha, Form("%s, N = %.0f", la.c_str(), a->GetEntries()), "l");
  leg->AddEntry(hb, Form("%s, N = %.0f", lb.c_str(), b->GetEntries()), "l");
  leg->Draw();
  DrawHeader(ps, hdr, sub1, sub2);
  CMS_lumi(c, 13, 10);
  c->RedrawAxis();
  c->SaveAs((outPath + ".png").c_str());
  c->SaveAs((outPath + ".pdf").c_str());
  delete c;
}

// 2D flip rate map (pT x y), text + color.
void Draw2DRate(TH2D *flip, TH2D *tot, const std::string &outPath,
                const std::string &hdr, const std::string &sub1)
{
  if (!flip || !tot) return;
  TH2D *r = (TH2D *)flip->Clone(Form("%s_rate", flip->GetName()));
  r->SetDirectory(nullptr);
  r->Divide(flip, tot, 1, 1, "B");
  // Color only (the per-bin numbers are in the log tables; text labels do not
  // fit the 5 GeV pT bins), log z since the rate spans 1e-4 .. 2e-2.
  PlotStyle ps;
  ps.lm = 0.15; ps.rm = 0.20; ps.yTitleOffset = 1.55;
  ps.headerX = 0.15; ps.headerY = 0.97; ps.headerDy = 0.04;
  ps.titleSize = 0.038; ps.subSize = 0.030;
  static int n = 0;
  TCanvas *c = new TCanvas(Form("c_2d_%d", n++), "", ps.w, ps.h);
  ApplyCanvasStyle(c, ps);
  c->SetTopMargin(0.11);
  c->SetLogz();
  ApplyHistStyle(r, ps, r->GetXaxis()->GetTitle(), r->GetYaxis()->GetTitle());
  r->SetStats(0);
  r->GetZaxis()->SetTitle("charge-flip rate");
  r->GetZaxis()->SetTitleFont(ps.font);
  r->GetZaxis()->SetLabelFont(ps.font);
  r->GetZaxis()->SetTitleOffset(1.6);
  r->Draw("COLZ");
  DrawHeader(ps, hdr, sub1, "");
  c->SaveAs((outPath + ".png").c_str());
  c->SaveAs((outPath + ".pdf").c_str());
  delete c;
}

// ============================================================
// Tables (stdout = the record) + plots + rootfile, per flavor
// ============================================================
void Report(bool isMu, FlipOut &oP, FlipOut &oM, TFile *fout)
{
  const std::string flav = isMu ? "mu" : "ele";
  const std::string lepSym = isMu ? "#mu" : "e";
  const std::string outDir = "./plots/charge_flip_" + flav;
  gSystem->mkdir(outDir.c_str(), kTRUE);

  FlipOut oA;
  oA.CombineFrom(oP, oM, flav + "_all");

  FlipOut *S[3] = {&oP, &oM, &oA};
  const char *sName[3]  = {"W+ (Wp MC)", "W- (Wm MC)", "W+ + W-"};
  const char *sShort[3] = {"Wp", "Wm", "all"};

  std::cout << "\n==================== CHARGE-FLIP REPORT: " << flav << " ====================\n";
  std::cout << Form("[CONFIG] %s: gen match |pdg|=%d, DR<%.2f, |dpT|/pT_gen<%.2f, CHARGE-BLIND; "
                    "selection = skim W 8-step nominal; CP %.1f%% intervals on raw counts\n",
                    flav.c_str(), isMu ? 13 : 11, kMatchDR, kMatchDPt, 100 * kCL);

  // ---- inclusive ----
  for (int s = 0; s < 3; ++s)
  {
    const Counters &c = S[s]->c;
    const CP f  = Rate(c.nFlip, c.nMatched);
    const CP fu = Rate(c.nUnm, c.nSel);
    std::cout << Form("[RESULT] %-3s %-11s sel=%7.0f matched=%7.0f flip=%5.0f  f = %s   (weighted f = %.3e)   unmatched = %.0f (%.2f%%)\n",
                      flav.c_str(), sName[s], c.nSel, c.nMatched, c.nFlip, FmtRate(f).c_str(),
                      c.wMatched > 0 ? c.wFlip / c.wMatched : 0.0, c.nUnm, 100 * fu.f);
  }
  std::cout << Form("[INFO] %s: matches with a W in the gen mother chain: %.0f of %.0f (expected 0: the filtered gen "
                    "collection stores no |pdg|=24 entry, so the DR-nearest candidate is used); "
                    "matched gen charge != sample charge: %.0f\n",
                    flav.c_str(), oA.c.nMatched - oA.c.nNoWAnc, oA.c.nMatched, oA.c.nGenChgNe);

  // ---- vs y ----
  std::cout << Form("\n[TABLE] %s flip rate vs y (analysis bins, y = -eta_lab)\n", flav.c_str());
  std::cout << Form("%4s %6s %6s | %7s %5s %-28s | %7s %5s %-28s | %7s %5s %-28s\n",
                    "bin", "ylo", "yhi", "N+", "flip", "f+ (CP68)", "N-", "flip", "f- (CP68)",
                    "Nall", "flip", "f_all (CP68)");
  for (int b = 1; b <= kNY; ++b)
  {
    std::cout << Form("%4d %6.2f %6.2f", b - 1, kYEdges[b - 1], kYEdges[b]);
    for (int s = 0; s < 3; ++s)
    {
      const double t = S[s]->H1("tot_y")->GetBinContent(b), p = S[s]->H1("flip_y")->GetBinContent(b);
      std::cout << Form(" | %7.0f %5.0f %-28s", t, p, FmtRate(Rate(p, t)).c_str());
    }
    std::cout << "\n";
  }

  // ---- vs pT ----
  std::cout << Form("\n[TABLE] %s flip rate vs lepton pT (last bin includes pT > %.0f)\n", flav.c_str(), kPtEdges[kNPt]);
  std::cout << Form("%4s %6s %6s | %7s %5s %-28s | %7s %5s %-28s | %7s %5s %-28s\n",
                    "bin", "ptlo", "pthi", "N+", "flip", "f+ (CP68)", "N-", "flip", "f- (CP68)",
                    "Nall", "flip", "f_all (CP68)");
  for (int b = 1; b <= kNPt; ++b)
  {
    std::cout << Form("%4d %6.0f %6.0f", b - 1, kPtEdges[b - 1], kPtEdges[b]);
    for (int s = 0; s < 3; ++s)
    {
      const double t = S[s]->H1("tot_pt")->GetBinContent(b), p = S[s]->H1("flip_pt")->GetBinContent(b);
      std::cout << Form(" | %7.0f %5.0f %-28s", t, p, FmtRate(Rate(p, t)).c_str());
    }
    std::cout << "\n";
  }

  // ---- match quality ----
  {
    TH1D *dOk = oA.H1("dR_ok"), *dFl = oA.H1("dR_flip");
    auto fracBelow = [](TH1D *h, double logDR) -> double
    {
      const double tot = h->Integral(0, h->GetNbinsX() + 1);
      return tot > 0 ? h->Integral(0, h->FindFixBin(logDR - 1e-6)) / tot : 0.0;
    };
    std::cout << Form("\n[INFO] %s match quality (all matched): <log10 DR> correct %.2f / flip %.2f; "
                      "fraction with DR < 0.01: correct %.4f / flip %.4f; DR < 0.1: correct %.4f / flip %.4f; "
                      "<dpT/pT> correct %+.4f / flip %+.4f\n",
                      flav.c_str(), dOk->GetMean(), dFl->GetMean(),
                      fracBelow(dOk, -2.0), fracBelow(dFl, -2.0), fracBelow(dOk, -1.0), fracBelow(dFl, -1.0),
                      oA.H1("dpt_ok")->GetMean(), oA.H1("dpt_flip")->GetMean());
  }

  // ---- asymmetry bias: A_reco - A_gen on the matched leptons, k_s-combined ----
  const double kP = pONorm::MCScale("Wp_" + flav);
  const double kM = pONorm::MCScale("Wm_" + flav);
  std::cout << Form("\n[TABLE] %s charge-asymmetry bias from flips (matched leptons, Wp/Wm combined with k_s = %.4g / %.4g)\n",
                    flav.c_str(), kP, kM);
  std::cout << Form("%4s %6s %6s | %9s %9s %9s | %9s %9s | %9s\n", "bin", "ylo", "yhi",
                    "A_gen", "A_reco", "dA", "f_all", "-2fA_gen", "dA/A_gen");
  double sRP = 0, sRM = 0, sGP = 0, sGM = 0;
  for (int b = 1; b <= kNY; ++b)
  {
    const double rp = kP * oP.H1("w_recoP_y")->GetBinContent(b) + kM * oM.H1("w_recoP_y")->GetBinContent(b);
    const double rm = kP * oP.H1("w_recoM_y")->GetBinContent(b) + kM * oM.H1("w_recoM_y")->GetBinContent(b);
    const double gp = kP * oP.H1("w_genP_y")->GetBinContent(b)  + kM * oM.H1("w_genP_y")->GetBinContent(b);
    const double gm = kP * oP.H1("w_genM_y")->GetBinContent(b)  + kM * oM.H1("w_genM_y")->GetBinContent(b);
    sRP += rp; sRM += rm; sGP += gp; sGM += gm;
    const double aR = (rp + rm > 0) ? (rp - rm) / (rp + rm) : 0.0;
    const double aG = (gp + gm > 0) ? (gp - gm) / (gp + gm) : 0.0;
    const double t = oA.H1("tot_y")->GetBinContent(b), p = oA.H1("flip_y")->GetBinContent(b);
    const double fA = (t > 0) ? p / t : 0.0;
    std::cout << Form("%4d %6.2f %6.2f | %9.4f %9.4f %+9.5f | %9.2e %+9.5f | %+9.4f\n",
                      b - 1, kYEdges[b - 1], kYEdges[b], aG, aR, aR - aG, fA, -2.0 * fA * aG,
                      aG != 0 ? (aR - aG) / aG : 0.0);
  }
  {
    const double aR = (sRP + sRM > 0) ? (sRP - sRM) / (sRP + sRM) : 0.0;
    const double aG = (sGP + sGM > 0) ? (sGP - sGM) / (sGP + sGM) : 0.0;
    const double fA = oA.c.nMatched > 0 ? oA.c.nFlip / oA.c.nMatched : 0.0;
    std::cout << Form("%4s %6s %6s | %9.4f %9.4f %+9.5f | %9.2e %+9.5f | %+9.4f\n", "incl", "", "",
                      aG, aR, aR - aG, fA, -2.0 * fA * aG, aG != 0 ? (aR - aG) / aG : 0.0);
    std::cout << Form("[RESULT] %s inclusive asymmetry bias: A_gen = %.4f, A_reco = %.4f, dA = %+.5f (relative %+.3f%%)\n",
                      flav.c_str(), aG, aR, aR - aG, aG != 0 ? 100 * (aR - aG) / aG : 0.0);
  }

  // ---- electron extras ----
  if (!isMu)
  {
    const Counters &c = oA.c;
    const double totTrk = oA.H1("tot_y_trk")->Integral(), flTrk = oA.H1("flip_y_trk")->Integral();
    const double totCon = oA.H1("tot_y_cons")->Integral(), flCon = oA.H1("flip_y_cons")->Integral();
    std::cout << "\n[RESULT] ele charge definitions (matched electrons):\n";
    std::cout << Form("         eleCharge (analysis)            f = %s\n", FmtRate(Rate(c.nFlip, c.nMatched)).c_str());
    std::cout << Form("         eleTrkCharge                     f = %s\n", FmtRate(Rate(flTrk, totTrk)).c_str());
    std::cout << Form("         eleCharge == eleTrkCharge subset f = %s   (subset = %.0f of %.0f, %.2f%%)\n",
                      FmtRate(Rate(flCon, totCon)).c_str(), totCon, c.nMatched,
                      c.nMatched > 0 ? 100 * totCon / c.nMatched : 0.0);
    std::cout << Form("         inconsistent (eleCharge != eleTrkCharge): %.0f (%.2f%% of matched), of which flipped %.0f  ->  f(incons) = %s\n",
                      c.nIncons, c.nMatched > 0 ? 100 * c.nIncons / c.nMatched : 0.0, c.nInconsFlip,
                      FmtRate(Rate(c.nInconsFlip, c.nIncons)).c_str());
    std::cout << Form("         a consistency requirement would drop %.2f%% of matched electrons and remove %.1f%% of the flips\n",
                      c.nMatched > 0 ? 100 * c.nIncons / c.nMatched : 0.0,
                      c.nFlip > 0 ? 100 * c.nInconsFlip / c.nFlip : 0.0);
  }

  // ---- rootfile ----
  if (fout)
  {
    oP.Write(fout); oM.Write(fout); oA.Write(fout);
  }

  // ---- TEfficiency objects + plots ----
  std::vector<TEfficiency *> eY, ePt, eUnm;
  const int cols[3] = {kRed + 1, kBlue + 1, kBlack};
  const int mks[3]  = {20, 21, 24};
  for (int s = 0; s < 3; ++s)
  {
    eY.push_back(MakeEff(S[s]->H1("flip_y"), S[s]->H1("tot_y"), Form("eff_flip_y_%s_%s", flav.c_str(), sShort[s]), "flip rate vs y"));
    ePt.push_back(MakeEff(S[s]->H1("flip_pt"), S[s]->H1("tot_pt"), Form("eff_flip_pt_%s_%s", flav.c_str(), sShort[s]), "flip rate vs pT"));
    eUnm.push_back(MakeEff(S[s]->H1("unm_y"), S[s]->H1("sel_y"), Form("eff_unm_y_%s_%s", flav.c_str(), sShort[s]), "unmatched fraction vs y"));
  }
  if (fout)
  {
    fout->cd();
    for (auto *e : eY)   if (e) e->Write("", TObject::kOverwrite);
    for (auto *e : ePt)  if (e) e->Write("", TObject::kOverwrite);
    for (auto *e : eUnm) if (e) e->Write("", TObject::kOverwrite);
  }

  const std::string hdr  = "W#rightarrow" + lepSym + "#nu signal MC";
  const std::string sub1 = "full W selection, charge-blind match";
  const std::string sub2 = Form("#DeltaR<%.1f, |#Deltap_{T}|/p_{T}^{gen}<%.1f", kMatchDR, kMatchDPt);
  const CP fIncl = Rate(oA.c.nFlip, oA.c.nMatched);
  const std::vector<std::string> box = {
      Form("f_{incl} = %.2e_{-%.1e}^{+%.1e}", fIncl.f, fIncl.lo, fIncl.hi),
      Form("N_{matched} = %.0f, N_{flip} = %.0f", oA.c.nMatched, oA.c.nFlip)};
  const std::vector<std::string> labs = {"W^{+}", "W^{-}", "W^{+}+W^{-}"};
  const std::vector<int> vcols(cols, cols + 3), vmks(mks, mks + 3);

  DrawRates(eY, labs, vcols, vmks, outDir + "/fliprate_vs_y",
            "y^{l}_{lab} = -#eta^{l}_{lab}", "Charge-flip rate", hdr, sub1, sub2, box,
            kYEdges[0], kYEdges[kNY], /*logy=*/!isMu);
  DrawRates(ePt, labs, vcols, vmks, outDir + "/fliprate_vs_pt",
            "p_{T}^{l} [GeV]", "Charge-flip rate", hdr, sub1, sub2, box,
            kPtEdges[0], kPtEdges[kNPt], /*logy=*/!isMu);
  DrawRates(eUnm, labs, vcols, vmks, outDir + "/unmatched_vs_y",
            "y^{l}_{lab} = -#eta^{l}_{lab}", "Unmatched fraction", hdr, sub1, sub2, {},
            kYEdges[0], kYEdges[kNY], /*logy=*/false);
  Draw2DRate(oA.H2("flip2"), oA.H2("tot2"), outDir + "/fliprate_2D_pt_y", hdr, sub1 + ", W^{+}+W^{-}");
  DrawTwoShapes(oA.H1("dR_ok"), oA.H1("dR_flip"), "correct charge", "flipped charge",
                outDir + "/match_dR", "log_{10} #DeltaR(gen, reco)", hdr, sub1, sub2);
  DrawTwoShapes(oA.H1("dpt_ok"), oA.H1("dpt_flip"), "correct charge", "flipped charge",
                outDir + "/match_dpt", "(p_{T}^{reco} - p_{T}^{gen}) / p_{T}^{gen}", hdr, sub1, sub2);

  if (!isMu)
  {
    std::vector<TEfficiency *> eDef = {
        MakeEff(oA.H1("flip_y"),      oA.H1("tot_y"),      Form("eff_flip_y_%s_all_eleCharge", flav.c_str()), "eleCharge"),
        MakeEff(oA.H1("flip_y_trk"),  oA.H1("tot_y_trk"),  Form("eff_flip_y_%s_all_trkCharge", flav.c_str()), "eleTrkCharge"),
        MakeEff(oA.H1("flip_y_cons"), oA.H1("tot_y_cons"), Form("eff_flip_y_%s_all_consistent", flav.c_str()), "consistent subset")};
    std::vector<TEfficiency *> eDefPt = {
        MakeEff(oA.H1("flip_pt"),      oA.H1("tot_pt"),      Form("eff_flip_pt_%s_all_eleCharge", flav.c_str()), "eleCharge"),
        MakeEff(oA.H1("flip_pt_trk"),  oA.H1("tot_pt_trk"),  Form("eff_flip_pt_%s_all_trkCharge", flav.c_str()), "eleTrkCharge"),
        MakeEff(oA.H1("flip_pt_cons"), oA.H1("tot_pt_cons"), Form("eff_flip_pt_%s_all_consistent", flav.c_str()), "consistent subset")};
    std::vector<TEfficiency *> eInc = {
        MakeEff(oA.H1("incons_y"), oA.H1("tot_y"), Form("eff_incons_y_%s_all", flav.c_str()), "inconsistent fraction")};
    if (fout)
    {
      fout->cd();
      for (auto *e : eDef)   if (e) e->Write("", TObject::kOverwrite);
      for (auto *e : eDefPt) if (e) e->Write("", TObject::kOverwrite);
      for (auto *e : eInc)   if (e) e->Write("", TObject::kOverwrite);
    }
    const std::vector<std::string> dl = {"eleCharge (analysis)", "eleTrkCharge", "eleCharge = eleTrkCharge"};
    const std::vector<int> dc = {kBlack, kGreen + 2, kMagenta + 1};
    const std::vector<int> dm = {24, 22, 23};
    const std::string sub2c = sub2 + ", W^{+}+W^{-}";
    DrawRates(eDef, dl, dc, dm, outDir + "/fliprate_vs_y_chargedef",
              "y^{e}_{lab} = -#eta^{e}_{lab}", "Charge-flip rate", hdr, sub1, sub2c, {},
              kYEdges[0], kYEdges[kNY], /*logy=*/true);
    DrawRates(eDefPt, dl, dc, dm, outDir + "/fliprate_vs_pt_chargedef",
              "p_{T}^{e} [GeV]", "Charge-flip rate", hdr, sub1, sub2c, {},
              kPtEdges[0], kPtEdges[kNPt], /*logy=*/true);
    DrawRates(eInc, {"eleCharge #neq eleTrkCharge"}, {kBlack}, {20}, outDir + "/incons_vs_y",
              "y^{e}_{lab} = -#eta^{e}_{lab}", "Fraction of matched electrons", hdr, sub1, sub2c, {},
              kYEdges[0], kYEdges[kNY], /*logy=*/false);
  }
  std::cout << "[OK] plots in " << outDir << "/\n";
}

} // anonymous namespace

// ============================================================
// Driver
// ============================================================
void charge_flip(const char *channel = "both", bool plotsOnly = false)
{
  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetOptStat(0);
  // motherIdx is a vector<vector<int>> branch -- same dictionary the skim needs
  gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");

  const std::string ch = channel;
  std::vector<bool> flavs;
  if (ch == "both" || ch == "all") flavs = {true, false};
  else if (ch == "mu"  || ch == "Wmu") flavs = {true};
  else if (ch == "ele" || ch == "Wel") flavs = {false};
  else
  {
    std::cerr << "[ERR] channel must be \"mu\", \"ele\" or \"both\" (got \"" << ch << "\")\n";
    return;
  }

  gSystem->mkdir("rootfile", kTRUE);
  gSystem->mkdir("plots", kTRUE);

  for (bool isMu : flavs)
  {
    const std::string flav = isMu ? "mu" : "ele";
    const std::string outRoot = "./rootfile/charge_flip_" + flav + ".root";
    std::cout << "\n========== charge_flip " << flav << " ==========\n";

    FlipOut oP, oM;
    if (plotsOnly)
    {
      TFile *fin = TFile::Open(outRoot.c_str(), "READ");
      if (!fin || fin->IsZombie())
      {
        std::cerr << "[ERR] cannot open " << outRoot << " (run without plotsOnly first); skipping\n";
        continue;
      }
      const bool ok = oP.Read(fin, flav + "_Wp", isMu) && oM.Read(fin, flav + "_Wm", isMu);
      fin->Close();
      if (!ok) { std::cerr << "[ERR] incomplete " << outRoot << "; skipping\n"; continue; }
      Report(isMu, oP, oM, nullptr);
      continue;
    }

    oP.Book(flav + "_Wp", isMu);
    oM.Book(flav + "_Wm", isMu);
    if (RunFlip(isMu, kWp, oP) != 0 || RunFlip(isMu, kWm, oM) != 0)
    {
      std::cerr << "[ERR] " << flav << ": MC loop failed; skipping\n";
      continue;
    }

    TFile *fout = TFile::Open(outRoot.c_str(), "RECREATE");
    Report(isMu, oP, oM, fout);
    fout->Close();
    std::cout << "[OK] wrote " << outRoot << "\n";
  }
  std::cout << "\n[DONE] charge_flip\n";
}
