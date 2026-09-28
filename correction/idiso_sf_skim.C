// correction/idiso_sf_skim.C
//
// THE SEPARATE SKIM of the electron ID+ISO scale-factor cross-check
// (2026-09-24). The efficiency definition, binning and histogram names live
// in correction/idiso_sf_common.h -- read its header first. skim/skim.C is NOT
// touched: this macro replicates what it needs (as trig_eff_mb.C and
// charge_flip.C do) because the study changes the electron definition itself.
//
// One job per sample (Data, Wp, Wm, DY, DYtau, Wptau, Wmtau), all from the
// single-sourced inputs (pOSkim::kDefaultDataFile / ResolveMCSample(s, "ele")).
// Per event:
//   (a) pO event filter (PassEventSelection_pO) and |vz| < 15
//   (b) HLT_OxyL1SingleEG10_v1 fired
//   Z TAG-AND-PROBE (idiso_sf_common.h, since 2026-09-25): legs = pT > 25,
//       |eta_SC| < 2.4, no crack, relIso < 0.3; tag = a passing leg matched to
//       an EG10 object; the OS pair with >= 1 tag and 60 <= m <= 120 closest to
//       the Z mass -> h_mZ_PP (both pass) or h_mZ_PF (one passes), once per
//       event; the same-sign twin -> h_mZss_*; the exact per-pair probe counts
//       -> h_tnp_probes; the OS probes' EGM-cell maps -> h2_egm_tnp_*
//   (c) a leading electron exists: highest pT with pT > 25, |eta_SC| < 2.4,
//       outside the crack -- NO ID, NO isolation
//   (d) DY veto: no OS pair of total-level legs (pT >= 10, |eta_SC| < 2.4, no
//       crack, relIso < 0.3) with 60 <= m <= 120 -- a superset of the Z legs,
//       so no event reaches both the W and the Z channels
//   (e) the leading electron is matched to an EG10 trigger object (DR < 0.4)
//   [MC: gen-matched "any" counter -- the denominator of eps(relIso < 0.3)]
//   (f) relIso < 1.0: relIso < 0.3 -> TOTAL (pass | fail),
//                     0.3-1.0      -> sideband (QCD MET templates)
// Filled per category x charge x scheme bin: PF MET, m_T (the in-fit ABCD axis)
// and the lepton pT at m_T > 40 (the fitted variable); plus Sum_w counters
// (gen-matched twins in MC) and the EGM-cell maps. No SFs, no LHE twins: this
// study only needs the nominal MC and the generator weight.
//
// Run from correction/:
//   ./run_idiso_sf.sh skim [samples]                  (keeps the logs)
//   root -l -b -q 'idiso_sf_skim.C+("Data")'          (bare)
//   root -l -b -q 'idiso_sf_skim.C+("Wp", 200000)'    (first N events, a quick test)
// Output: rootfile/idiso_sf_ele/skim_<sample>.root

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../skim/skim_common.h"
#include "idiso_sf_common.h"

using namespace pOSkim;
using namespace pOIdIso;

namespace
{

// cutflow steps (event counts, raw and weighted)
enum Step { kStAll = 0, kStFilt, kStTrig, kStLead, kStDyVeto, kStMatch, kStIso1, kStTotal, kStPass, kStFail, kStSb, kNStep };
const char *const kStepName[kNStep] = {
    "all events",
    "(a) event filter + |vz| < 15",
    "(b) HLT_OxyL1SingleEG10_v1 fired",
    "(c) leading e: pT > 25, |eta_SC| < 2.4, no crack",
    "(d) DY veto (total-level legs, 60-120)",
    "(e) leading e trigger-matched (DR < 0.4)",
    "(f) leading e relIso < 1.0",
    "    TOTAL: relIso < 0.3",
    "      pass: eleMVAIdWP90 && eleMVAIsoWP90",
    "      fail",
    "    sideband: relIso 0.3-1.0"};

bool InAcceptance(double pt, double scEta, double ptMin)
{
  return pt > ptMin && std::fabs(scEta) < kEtaSCMax && !InEcalGap(scEta);
}

} // namespace

int idiso_sf_skim(const char *sampleName = "Data", Long64_t nmax = -1)
{
  const int is = SampleIndex(sampleName);
  if (is < 0)
  {
    std::cerr << "[FATAL] idiso_sf_skim: unknown sample '" << sampleName
              << "' (Data Wp Wm DY DYtau Wptau Wmtau)\n";
    return 2;
  }
  const SampleType sample = kSampleType[is];
  const bool isMC = IsMC(sample);

  std::string fname = kDefaultDataFile;
  if (isMC)
  {
    const auto info = ResolveMCSample(sample, "ele");
    if (info.fname.empty()) { std::cerr << "[FATAL] cannot resolve MC sample " << sampleName << "\n"; return 2; }
    fname = info.fname;
  }
  std::cout << "[INPUT] " << sampleName << " <- " << fname << std::endl;
  std::cout << "[CONFIG] total: leading e pT > " << kPtMin << ", |eta_SC| < " << kEtaSCMax
            << " (no crack), trigger-matched DR < " << kTrigMatchDR << ", relIso < " << kIsoTotalMax
            << "; pass = " << kPassIdBranch << " && " << kPassIsoBranch
            << "; sideband relIso [" << kIsoSbLo << ", " << kIsoSbMid << ") + [" << kIsoSbMid << ", " << kIsoSbHi << ")"
            << "; DY veto / Z window " << kMllLo << "-" << kMllHi << "\n";

  TFile *f = TFile::Open(fname.c_str());
  if (!f || f->IsZombie()) { std::cerr << "[FATAL] cannot open " << fname << "\n"; return 2; }

  TTree *tEle    = (TTree *)f->Get("ggHiNtuplizer/EventTree");
  TTree *tHi     = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  TTree *tPF     = (TTree *)f->Get("particleFlowAnalyser/pftree");
  TTree *tHLT    = (TTree *)f->Get("hltanalysis/HltTree");
  TTree *tHLTobj = (TTree *)f->Get(kTrigObjTree);
  TTree *tEvent  = (TTree *)f->Get("skimanalysis/HltTree");
  if (!tEle || !tHi || !tPF || !tHLT || !tHLTobj || !tEvent)
  {
    std::cerr << "[FATAL] missing a required tree in " << fname << "\n";
    return 2;
  }
  for (TTree *t : {tHi, tPF, tHLT, tHLTobj, tEvent})
    if (t->GetEntries() != tEle->GetEntries())
    {
      std::cerr << "[FATAL] tree " << t->GetName() << " (" << t->GetEntries() << ") and EventTree ("
                << tEle->GetEntries() << ") entry counts differ\n";
      return 2;
    }

  // -------- EventTree: every branch read is enabled AND addressed (the pair
  //          is load-bearing after SetBranchStatus("*", 0)) --------
  Int_t nEle = 0;
  std::vector<float> *elePt = nullptr, *eleEta = nullptr, *elePhi = nullptr, *eleSCEta = nullptr;
  std::vector<int>   *eleCharge = nullptr, *eleIdWP90 = nullptr, *eleIsoWP90 = nullptr;
  std::vector<float> *chIso = nullptr, *neuIso = nullptr, *phoIso = nullptr, *puIso = nullptr;
  tEle->SetBranchStatus("*", 0);
  for (const char *bn : {"nEle", "elePt", "eleEta", "elePhi", "eleSCEta", "eleCharge", kPassIdBranch,
                         kPassIsoBranch, "elePFChIso", "elePFNeuIso", "elePFPhoIso", "elePFPUIso"})
    if (!HasBranch(tEle, bn))
    {
      std::cerr << "[FATAL] missing EventTree branch " << bn << "\n";
      return 2;
    }
  tEle->SetBranchStatus("nEle", 1);          tEle->SetBranchAddress("nEle", &nEle);
  tEle->SetBranchStatus("elePt", 1);         tEle->SetBranchAddress("elePt", &elePt);
  tEle->SetBranchStatus("eleEta", 1);        tEle->SetBranchAddress("eleEta", &eleEta);
  tEle->SetBranchStatus("elePhi", 1);        tEle->SetBranchAddress("elePhi", &elePhi);
  tEle->SetBranchStatus("eleSCEta", 1);      tEle->SetBranchAddress("eleSCEta", &eleSCEta);
  tEle->SetBranchStatus("eleCharge", 1);     tEle->SetBranchAddress("eleCharge", &eleCharge);
  tEle->SetBranchStatus(kPassIdBranch, 1);   tEle->SetBranchAddress(kPassIdBranch, &eleIdWP90);
  tEle->SetBranchStatus(kPassIsoBranch, 1);  tEle->SetBranchAddress(kPassIsoBranch, &eleIsoWP90);
  tEle->SetBranchStatus("elePFChIso", 1);    tEle->SetBranchAddress("elePFChIso", &chIso);
  tEle->SetBranchStatus("elePFNeuIso", 1);   tEle->SetBranchAddress("elePFNeuIso", &neuIso);
  tEle->SetBranchStatus("elePFPhoIso", 1);   tEle->SetBranchAddress("elePFPhoIso", &phoIso);
  tEle->SetBranchStatus("elePFPUIso", 1);    tEle->SetBranchAddress("elePFPUIso", &puIso);

  // gen block for the gen-matched eps_MC cross-check (MC; mandatory there)
  std::vector<int>   *mcPID = nullptr, *mcStatus = nullptr, *mcMomPID = nullptr, *mcGMomPID = nullptr;
  std::vector<float> *mcPt = nullptr, *mcEta = nullptr, *mcPhi = nullptr;
  if (isMC)
  {
    for (const char *bn : {"mcPID", "mcStatus", "mcPt", "mcEta", "mcPhi", "mcMomPID"})
      if (!HasBranch(tEle, bn))
      {
        std::cerr << "[FATAL] MC file lacks EventTree branch " << bn << " (gen-matched eps_MC)\n";
        return 2;
      }
    tEle->SetBranchStatus("mcPID", 1);    tEle->SetBranchAddress("mcPID", &mcPID);
    tEle->SetBranchStatus("mcStatus", 1); tEle->SetBranchAddress("mcStatus", &mcStatus);
    tEle->SetBranchStatus("mcPt", 1);     tEle->SetBranchAddress("mcPt", &mcPt);
    tEle->SetBranchStatus("mcEta", 1);    tEle->SetBranchAddress("mcEta", &mcEta);
    tEle->SetBranchStatus("mcPhi", 1);    tEle->SetBranchAddress("mcPhi", &mcPhi);
    tEle->SetBranchStatus("mcMomPID", 1); tEle->SetBranchAddress("mcMomPID", &mcMomPID);
    if (HasBranch(tEle, "mcGMomPID")) { tEle->SetBranchStatus("mcGMomPID", 1); tEle->SetBranchAddress("mcGMomPID", &mcGMomPID); }
    else std::cout << "[WARN] no mcGMomPID; FSR chains (e <- e <- W) will not be gen-matched.\n";
  }

  // -------- trigger bit (exact name: FindBranchContaining may return a _Prescale twin) --------
  std::string hltName = HasBranch(tHLT, kTrigPath) ? std::string(kTrigPath) : FindBranchContaining(tHLT, kTrigPath);
  if (hltName.empty()) { std::cerr << "[FATAL] hltanalysis/HltTree lacks " << kTrigPath << "\n"; return 2; }
  Int_t hltBit = 0;
  tHLT->SetBranchStatus("*", 0);
  tHLT->SetBranchStatus(hltName.c_str(), 1); tHLT->SetBranchAddress(hltName.c_str(), &hltBit);

  // -------- trigger objects --------
  std::vector<double> *toPt = nullptr, *toEta = nullptr, *toPhi = nullptr;
  if (!HasBranch(tHLTobj, "pt") || !HasBranch(tHLTobj, "eta") || !HasBranch(tHLTobj, "phi"))
  {
    std::cerr << "[FATAL] no pt/eta/phi in " << kTrigObjTree << "\n";
    return 2;
  }
  tHLTobj->SetBranchStatus("*", 0);
  tHLTobj->SetBranchStatus("pt", 1);  tHLTobj->SetBranchAddress("pt", &toPt);
  tHLTobj->SetBranchStatus("eta", 1); tHLTobj->SetBranchAddress("eta", &toEta);
  tHLTobj->SetBranchStatus("phi", 1); tHLTobj->SetBranchAddress("phi", &toPhi);

  // -------- event filters, vz, gen weight --------
  Int_t ppv = 1, pcc = 1;
  const bool has_ppv = HasBranch(tEvent, "pprimaryVertexFilter");
  const bool has_pcc = HasBranch(tEvent, "pclusterCompatibilityFilter");
  tEvent->SetBranchStatus("*", 0);
  if (has_ppv) { tEvent->SetBranchStatus("pprimaryVertexFilter", 1);        tEvent->SetBranchAddress("pprimaryVertexFilter", &ppv); }
  if (has_pcc) { tEvent->SetBranchStatus("pclusterCompatibilityFilter", 1); tEvent->SetBranchAddress("pclusterCompatibilityFilter", &pcc); }

  Float_t vz = 0.f, genWeight = 1.f;
  tHi->SetBranchStatus("*", 0);
  if (!HasBranch(tHi, "vz")) { std::cerr << "[FATAL] HiTree lacks vz\n"; return 2; }
  tHi->SetBranchStatus("vz", 1); tHi->SetBranchAddress("vz", &vz);
  const bool has_genWeight = isMC && HasBranch(tHi, "weight");
  if (has_genWeight) { tHi->SetBranchStatus("weight", 1); tHi->SetBranchAddress("weight", &genWeight); }
  else if (isMC) std::cout << "[WARN] no HiTree 'weight' -- MC filled unweighted.\n";

  // -------- PF tree: MET, read only for selected events --------
  std::vector<float> *pfPt = nullptr, *pfPhi = nullptr;
  if (!HasBranch(tPF, "pfPt") || !HasBranch(tPF, "pfPhi")) { std::cerr << "[FATAL] pftree lacks pfPt/pfPhi\n"; return 2; }
  tPF->SetBranchStatus("*", 0);
  tPF->SetBranchStatus("pfPt", 1);  tPF->SetBranchAddress("pfPt", &pfPt);
  tPF->SetBranchStatus("pfPhi", 1); tPF->SetBranchAddress("pfPhi", &pfPhi);

  // -------- histograms --------
  TH1::SetDefaultSumw2(kTRUE);
  TH1D *hV[kNDisc][kNCat][kNChg][kNScheme][kMaxBins] = {};
  for (int d = 0; d < kNDisc; ++d)
    for (int c = 0; c < kNCat; ++c)
      for (int q = 0; q < kNChg; ++q)
        for (int s = 0; s < kNScheme; ++s)
          for (int k = 0; k < kSchemes[s].n; ++k)
          {
            const std::string n = HName(d, c, q, s, k);
            hV[d][c][q][s][k] = new TH1D(n.c_str(), Form(";%s;Events", kDiscs[d].axisTitle),
                                         kDiscs[d].nb, kDiscs[d].lo, kDiscs[d].hi);
            hV[d][c][q][s][k]->SetDirectory(nullptr);
          }
  TH1D *hCnt[kNCat][kNChg][kNScheme] = {}, *hCntGm[kNCat][kNChg][kNScheme] = {}, *hCntGmAny[kNChg][kNScheme] = {};
  for (int q = 0; q < kNChg; ++q)
    for (int s = 0; s < kNScheme; ++s)
    {
      const int nb = kSchemes[s].n;
      for (int c = 0; c < kNCat; ++c)
      {
        hCnt[c][q][s]   = new TH1D(CntName(c, q, s).c_str(),   ";coarse bin;#Sigma w", nb, -0.5, nb - 0.5);
        hCntGm[c][q][s] = new TH1D(CntGmName(c, q, s).c_str(), ";coarse bin;#Sigma w (gen-matched W e)", nb, -0.5, nb - 0.5);
        hCnt[c][q][s]->SetDirectory(nullptr);
        hCntGm[c][q][s]->SetDirectory(nullptr);
      }
      hCntGmAny[q][s] = new TH1D(CntGmAnyName(q, s).c_str(), ";coarse bin;#Sigma w (gen-matched, any relIso)", nb, -0.5, nb - 0.5);
      hCntGmAny[q][s]->SetDirectory(nullptr);
    }
  TH2D *hEgm[2][kNChg] = {};
  for (int p = 0; p < 2; ++p)
    for (int q = 0; q < kNChg; ++q)
    {
      hEgm[p][q] = new TH2D(EgmMapName(p == 1, q).c_str(), ";#eta_{SC};p_{T} (GeV)",
                            kNEgmEta, kEgmEta, kNEgmPt, kEgmPt);
      hEgm[p][q]->SetDirectory(nullptr);
    }
  TH1D *hZ[kNZCat][2] = {}; // [category][0 = OS, 1 = SS]
  for (int zc = 0; zc < kNZCat; ++zc)
    for (int ss = 0; ss < 2; ++ss)
    {
      hZ[zc][ss] = new TH1D(ZHistName(zc, ss == 1).c_str(), ";m_{ee} (GeV);Events", kZNb, kZLo, kZHi);
      hZ[zc][ss]->SetDirectory(nullptr);
    }
  TH1D *hTnpProbes = new TH1D(kTnpProbeHist, ";probe;#Sigma w", 4, -0.5, 3.5);
  hTnpProbes->SetDirectory(nullptr);
  const char *const probeLab[4] = {"OS pass", "OS fail", "SS pass", "SS fail"};
  for (int i = 0; i < 4; ++i) hTnpProbes->GetXaxis()->SetBinLabel(i + 1, probeLab[i]);
  TH2D *hEgmTnp[2] = {};
  for (int p = 0; p < 2; ++p)
  {
    hEgmTnp[p] = new TH2D(EgmTnpMapName(p == 1).c_str(), ";#eta_{SC};p_{T} (GeV)", kNEgmEta, kEgmEta, kNEgmPt, kEgmPt);
    hEgmTnp[p]->SetDirectory(nullptr);
  }
  TH1D *hFlow = new TH1D("h_cutflow", ";;events", kNStep, -0.5, kNStep - 0.5);
  TH1D *hFlowW = new TH1D("h_cutflow_w", ";;#Sigma w", kNStep, -0.5, kNStep - 0.5);
  hFlow->SetDirectory(nullptr);
  hFlowW->SetDirectory(nullptr);
  for (int i = 0; i < kNStep; ++i)
  {
    hFlow->GetXaxis()->SetBinLabel(i + 1, kStepName[i]);
    hFlowW->GetXaxis()->SetBinLabel(i + 1, kStepName[i]);
  }

  // -------- event loop --------
  const Long64_t nEntries = tEle->GetEntries();
  const Long64_t nRun = (nmax > 0 && nmax < nEntries) ? nmax : nEntries;
  std::cout << "Processing entries: " << nRun << " of " << nEntries << "\n";
  bool warnedFilters = false, warnedTrig = false;
  Long64_t nGmTried = 0, nGmFail = 0, nZ[kNZCat][2] = {};
  auto step = [&](int i, double w) { hFlow->Fill(i); hFlowW->Fill(i, w); };

  for (Long64_t ie = 0; ie < nRun; ++ie)
  {
    if (ie % 500000 == 0) std::cout << "  event " << ie << "/" << nRun << "\n";
    tEle->GetEntry(ie);
    tHi->GetEntry(ie);
    tEvent->GetEntry(ie);
    tHLT->GetEntry(ie);

    const double w = has_genWeight ? (double)genWeight : 1.0;
    step(kStAll, w);
    if (!elePt || !eleEta || !elePhi || !eleSCEta || !eleCharge || !eleIdWP90 || !eleIsoWP90) continue;

    // (a) event filter + vz
    if (!PassEventSelection_pO(warnedFilters, has_ppv, ppv, has_pcc, pcc)) continue;
    if (TMath::Abs(vz) > kVzMax) continue;
    step(kStFilt, w);
    // (b) trigger bit
    if (!TriggerFired(hltBit)) continue;
    step(kStTrig, w);

    // relIso once per electron (the DY veto, the Z legs and the leading electron all use it)
    const int ne = (int)elePt->size();
    std::vector<double> iso(ne, 999.0);
    for (int i = 0; i < ne; ++i) iso[i] = RelIsoPF(i, elePt, chIso, neuIso, phoIso, puIso);
    auto pairMass = [&](int i, int j)
    {
      TLorentzVector a, b;
      a.SetPtEtaPhiM(elePt->at(i), eleEta->at(i), elePhi->at(i), ELE_MASS);
      b.SetPtEtaPhiM(elePt->at(j), eleEta->at(j), elePhi->at(j), ELE_MASS);
      return (a + b).M();
    };
    // total-level leg: acceptance (pT >= 10) + relIso < 0.3, no ID -- the DY-veto leg
    auto vetoLeg = [&](int i)
    {
      return elePt->at(i) >= kDyLegPtMin && std::fabs(eleSCEta->at(i)) < kEtaSCMax &&
             !InEcalGap(eleSCEta->at(i)) && iso[i] < kIsoTotalMax;
    };
    auto passes = [&](int i) { return eleIdWP90->at(i) != 0 && eleIsoWP90->at(i) != 0; };
    // the trigger objects, read once: the Z tag and step (e) both match against them
    tHLTobj->GetEntry(ie);
    auto matched = [&](int i)
    {
      return PassLeadingLeptonTrigMatch(kTrigMatchDR, i, eleEta, elePhi, true, toPt, true, toEta, true, toPhi,
                                        warnedTrig, "electron");
    };

    // Z TAG-AND-PROBE -- before the W steps (every TnP leg is a DY-veto leg, so
    // step (d) removes these events from the W channels). Per charge
    // combination (OS = the fit, SS = the cut-and-count background): the pair
    // of legs with 60 <= m <= 120 and at least one tag, closest to the Z mass.
    std::vector<int> leg;
    for (int i = 0; i < ne; ++i)
      if (InAcceptance(elePt->at(i), eleSCEta->at(i), kPtMin) && iso[i] < kIsoTotalMax) leg.push_back(i);
    if (leg.size() >= 2)
    {
      std::vector<char> pas(leg.size(), 0), tag(leg.size(), 0);
      for (size_t a = 0; a < leg.size(); ++a)
      {
        pas[a] = passes(leg[a]);
        tag[a] = pas[a] && matched(leg[a]);
      }
      for (int ss = 0; ss < 2; ++ss)
      {
        int ba = -1, bb = -1;
        double bm = 0.0, bd = 1e9;
        for (size_t a = 0; a < leg.size(); ++a)
          for (size_t b = a + 1; b < leg.size(); ++b)
          {
            const bool os = eleCharge->at(leg[a]) * eleCharge->at(leg[b]) < 0;
            if (os != (ss == 0)) continue;
            if (!tag[a] && !tag[b]) continue;
            const double m = pairMass(leg[a], leg[b]);
            if (m < kMllLo || m > kMllHi) continue;
            if (std::fabs(m - kZMass) < bd) { bd = std::fabs(m - kZMass); ba = (int)a; bb = (int)b; bm = m; }
          }
        if (ba < 0) continue;
        const int zc = (pas[ba] && pas[bb]) ? kZPP : kZPF; // PF: the passing leg is the tag
        hZ[zc][ss]->Fill(bm, w);
        ++nZ[zc][ss];
        // exact per-pair probes: a leg is a probe when its partner is a tag
        for (int side = 0; side < 2; ++side)
        {
          const int self = side ? bb : ba, other = side ? ba : bb;
          if (!tag[other]) continue;
          hTnpProbes->Fill(2 * ss + (pas[self] ? 0 : 1), w);
          if (ss == 0)
          {
            hEgmTnp[0]->Fill(eleSCEta->at(leg[self]), elePt->at(leg[self]), w);
            if (pas[self]) hEgmTnp[1]->Fill(eleSCEta->at(leg[self]), elePt->at(leg[self]), w);
          }
        }
      }
    }

    // (c) leading electron in acceptance -- no ID, no isolation
    int iLead = -1;
    for (int i = 0; i < ne; ++i)
    {
      if (!InAcceptance(elePt->at(i), eleSCEta->at(i), kPtMin)) continue;
      if (iLead < 0 || elePt->at(i) > elePt->at(iLead)) iLead = i;
    }
    if (iLead < 0) continue;
    step(kStLead, w);

    // (d) DY veto on total-level legs
    bool veto = false;
    for (int i = 0; i < ne && !veto; ++i)
    {
      if (!vetoLeg(i)) continue;
      for (int j = i + 1; j < ne && !veto; ++j)
      {
        if (!vetoLeg(j) || eleCharge->at(i) * eleCharge->at(j) >= 0) continue;
        const double m = pairMass(i, j);
        if (m >= kMllLo && m <= kMllHi) veto = true;
      }
    }
    if (veto) continue;
    step(kStDyVeto, w);

    // (e) trigger match of the leading electron
    if (!matched(iLead)) continue;
    step(kStMatch, w);

    const double pt    = elePt->at(iLead);
    const double scEta = eleSCEta->at(iLead);
    const double isoL  = iso[iLead];
    const int    q     = ChgIndex(eleCharge->at(iLead));
    int bins[kNScheme];
    for (int s = 0; s < kNScheme; ++s) bins[s] = FindBin(s, scEta, pt);

    // MC: is the leading electron a prompt W electron? (eps_MC cross-check only)
    bool gm = false;
    if (isMC)
    {
      gm = MatchGenLeptonFromW(pt, eleEta->at(iLead), elePhi->at(iLead), 11, mcPID, mcStatus, mcPt, mcEta,
                               mcPhi, mcMomPID, mcGMomPID) >= 0;
      if (gm)
        for (int s = 0; s < kNScheme; ++s)
          if (bins[s] >= 0) hCntGmAny[q][s]->Fill(bins[s], w);
    }

    // (f) relIso < 1.0 -> total (< 0.3) or sideband (0.3-1.0)
    if (isoL >= kIsoSbHi) continue;
    step(kStIso1, w);

    tPF->GetEntry(ie);
    const TVector2 metv = ComputePFMET(nullptr, pfPt, pfPhi);
    const double   mt   = TransverseMass(pt, elePhi->at(iLead), metv);
    const double   val[kNDisc] = {metv.Mod(), mt, pt};
    const bool     use[kNDisc] = {true, true, mt > kSRMtMin}; // leppt_mt40: the SR / CRC pT at m_T > 40

    std::vector<int> cats;
    if (isoL < kIsoTotalMax)
    {
      const bool pass = passes(iLead);
      cats.push_back(pass ? kPass : kFail);
      step(kStTotal, w);
      step(pass ? kStPass : kStFail, w);
      if (isMC) { ++nGmTried; if (!gm) ++nGmFail; }
      if (mt > kSRMtMin) // the EGM prediction of the FITTED efficiency: the in-fit-ABCD SR, m_T > 40
      {
        hEgm[0][q]->Fill(scEta, pt, w);
        if (pass) hEgm[1][q]->Fill(scEta, pt, w);
      }
    }
    else
    {
      const bool lo = isoL < kIsoSbMid;
      cats.push_back(lo ? kSbAllLo : kSbAllHi);
      if (eleIdWP90->at(iLead) != 0) cats.push_back(lo ? kSbIdLo : kSbIdHi);
      step(kStSb, w);
    }
    for (int c : cats)
      for (int s = 0; s < kNScheme; ++s)
      {
        const int k = bins[s];
        if (k < 0) continue;
        for (int d = 0; d < kNDisc; ++d)
          if (use[d]) hV[d][c][q][s][k]->Fill(val[d], w);
        hCnt[c][q][s]->Fill(k, w);
        if (gm) hCntGm[c][q][s]->Fill(k, w);
      }
  }

  // -------- record --------
  std::cout << "\n[CUTFLOW] " << sampleName << "   (events | sum of weights)\n";
  for (int i = 0; i < kNStep; ++i)
    std::cout << Form("  %-50s %10.0f   %14.1f\n", kStepName[i], hFlow->GetBinContent(i + 1), hFlowW->GetBinContent(i + 1));
  std::cout << Form("[RESULT] %s: total %.0f (pass %.0f, fail %.0f), sideband %.0f\n", sampleName,
                    hFlow->GetBinContent(kStTotal + 1), hFlow->GetBinContent(kStPass + 1), hFlow->GetBinContent(kStFail + 1),
                    hFlow->GetBinContent(kStSb + 1));
  {
    const double pP = hTnpProbes->GetBinContent(1), pF = hTnpProbes->GetBinContent(2);
    const double sP = hTnpProbes->GetBinContent(3), sF = hTnpProbes->GetBinContent(4);
    std::cout << Form("[RESULT] %s: Z TnP events OS PP %lld, PF %lld | SS PP %lld, PF %lld;"
                      " probes (Sum w) OS pass %.1f, fail %.1f | SS pass %.1f, fail %.1f;"
                      " OS probe eps %.4f (no SS subtraction)\n",
                      sampleName, nZ[kZPP][0], nZ[kZPF][0], nZ[kZPP][1], nZ[kZPF][1], pP, pF, sP, sF,
                      (pP + pF) > 0 ? pP / (pP + pF) : 0.0);
  }
  if (sample == kWp || sample == kWm) // the gen match (prompt W electron) means something only here
  {
    const double tot = hFlowW->GetBinContent(kStTotal + 1), pas = hFlowW->GetBinContent(kStPass + 1);
    std::cout << Form("[RESULT] %s: weighted eps(pass | total) = %.4f   (all events of the sample, not gen-matched)\n",
                      sampleName, tot > 0 ? pas / tot : 0.0);
    double gAny = 0, gTot = 0, gPass = 0;
    for (int q = 0; q < kNChg; ++q)
    {
      gAny  += hCntGmAny[q][0]->GetBinContent(1);
      gTot  += hCntGm[kPass][q][0]->GetBinContent(1) + hCntGm[kFail][q][0]->GetBinContent(1);
      gPass += hCntGm[kPass][q][0]->GetBinContent(1);
    }
    std::cout << Form("[GENMATCH] %s: %lld total-category events, %lld not a prompt W electron (%.3f%%);"
                      " gen-matched eps(pass | total) = %.4f, eps(relIso < 0.3) = %.4f\n",
                      sampleName, nGmTried, nGmFail, nGmTried > 0 ? 100.0 * nGmFail / nGmTried : 0.0,
                      gTot > 0 ? gPass / gTot : 0.0, gAny > 0 ? gTot / gAny : 0.0);
  }

  gSystem->mkdir(kRootDir, kTRUE);
  const std::string out = SkimFile(sampleName);
  TFile fo(out.c_str(), "RECREATE");
  if (fo.IsZombie()) { std::cerr << "[FATAL] cannot write " << out << "\n"; return 2; }
  for (int d = 0; d < kNDisc; ++d)
    for (int c = 0; c < kNCat; ++c)
      for (int q = 0; q < kNChg; ++q)
        for (int s = 0; s < kNScheme; ++s)
          for (int k = 0; k < kSchemes[s].n; ++k) hV[d][c][q][s][k]->Write();
  for (int q = 0; q < kNChg; ++q)
    for (int s = 0; s < kNScheme; ++s)
    {
      for (int c = 0; c < kNCat; ++c) { hCnt[c][q][s]->Write(); if (isMC) hCntGm[c][q][s]->Write(); }
      if (isMC) hCntGmAny[q][s]->Write();
    }
  for (int p = 0; p < 2; ++p)
    for (int q = 0; q < kNChg; ++q) hEgm[p][q]->Write();
  for (int zc = 0; zc < kNZCat; ++zc)
    for (int ss = 0; ss < 2; ++ss) hZ[zc][ss]->Write();
  hTnpProbes->Write();
  for (int p = 0; p < 2; ++p) hEgmTnp[p]->Write();
  hFlow->Write();
  hFlowW->Write();
  fo.Close();
  f->Close();
  std::cout << "[INFO] wrote " << out << "\n";
  return 0;
}
