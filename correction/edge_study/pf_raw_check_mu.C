// SCRATCH (2026-10-20): the muon twin of pf_raw_check.C -- does the PF muon carry the
// same momentum as muPt in data and MC? For every analysis muon (muIDTight && muIsPF,
// pT > 25, |eta| < 2.4, Delta-beta relIso < 0.15; event filter, |vz| < 15,
// HLT_OxyL1SingleMuOpen fired), per region of |eta|, data vs W -> mu nu MC, the MEDIAN of
//   PF / muPt   (PF = sum pT of the PF candidates within DR < 0.1)
//   PFmu / muPt (the PF muon candidate, id 3, closest within DR < 0.05)
// Run from correction/edge_study: root -l -b -q pf_raw_check_mu.C+ | tee pf_raw_check_mu.log
#include "/Users/zhenghuang/pO_analysis/skim/skim_common.h"

#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

using namespace pOSkim;

namespace PRCM
{
const int kNReg = 3, kNV = 2;
const char *const kReg[kNReg] = {"|eta|>2.0", "1.2-2.0", "|eta|<1.2"};
const char *const kV[kNV] = {"PF/muPt", "PFmu/muPt"};
std::vector<double> val[2][kNReg][kNV];
long long nNoPFmu[2][kNReg];

inline int RegionOf(double eta)
{
  const double a = std::fabs(eta);
  return a > 2.0 ? 0 : (a > 1.2 ? 1 : 2);
}

void Loop(const std::string &fname, int g, Long64_t nMax)
{
  TFile *f = TFile::Open(fname.c_str());
  TTree *t = (TTree *)f->Get("ggHiNtuplizer/EventTree"), *h = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  TTree *pf = (TTree *)f->Get("particleFlowAnalyser/pftree"), *hlt = (TTree *)f->Get("hltanalysis/HltTree");
  TTree *ev = (TTree *)f->Get("skimanalysis/HltTree");
  std::vector<float> *pt = nullptr, *eta = nullptr, *phi = nullptr, *ch = nullptr, *neu = nullptr, *pho = nullptr, *pu = nullptr;
  std::vector<int> *id = nullptr, *isPF = nullptr;
  t->SetBranchStatus("*", 0);
  auto bind = [&](const char *b, auto *a) { t->SetBranchStatus(b, 1); t->SetBranchAddress(b, a); };
  bind("muPt", &pt); bind("muEta", &eta); bind("muPhi", &phi); bind("muIDTight", &id); bind("muIsPF", &isPF);
  bind("muPFChIso", &ch); bind("muPFNeuIso", &neu); bind("muPFPhoIso", &pho); bind("muPFPUIso", &pu);
  std::vector<int> *pfId = nullptr;
  std::vector<float> *pfPt = nullptr, *pfEta = nullptr, *pfPhi = nullptr;
  pf->SetBranchStatus("*", 0);
  for (const char *b : {"pfId", "pfPt", "pfEta", "pfPhi"}) pf->SetBranchStatus(b, 1);
  pf->SetBranchAddress("pfId", &pfId); pf->SetBranchAddress("pfPt", &pfPt);
  pf->SetBranchAddress("pfEta", &pfEta); pf->SetBranchAddress("pfPhi", &pfPhi);
  Int_t bit = 0, ppv = 1, pcc = 1;
  hlt->SetBranchStatus("*", 0);
  hlt->SetBranchStatus("HLT_OxyL1SingleMuOpen_v1", 1); hlt->SetBranchAddress("HLT_OxyL1SingleMuOpen_v1", &bit);
  ev->SetBranchStatus("*", 0);
  ev->SetBranchStatus("pprimaryVertexFilter", 1); ev->SetBranchAddress("pprimaryVertexFilter", &ppv);
  ev->SetBranchStatus("pclusterCompatibilityFilter", 1); ev->SetBranchAddress("pclusterCompatibilityFilter", &pcc);
  Float_t vz = 0;
  h->SetBranchStatus("*", 0);
  h->SetBranchStatus("vz", 1); h->SetBranchAddress("vz", &vz);
  bool warned = false;
  long long nMu = 0;
  const Long64_t n = (nMax > 0) ? std::min(nMax, t->GetEntries()) : t->GetEntries();
  for (Long64_t ie = 0; ie < n; ++ie)
  {
    t->GetEntry(ie); h->GetEntry(ie); ev->GetEntry(ie); hlt->GetEntry(ie);
    if (!pt || !id || !isPF) continue;
    if (!PassEventSelection_pO(warned, true, ppv, true, pcc) || std::fabs(vz) > 15.0 || !TriggerFired(bit)) continue;
    std::vector<int> sel;
    for (size_t i = 0; i < pt->size(); ++i)
      if (pt->at(i) > 25.0 && std::fabs(eta->at(i)) < 2.4 && id->at(i) != 0 && isPF->at(i) != 0 &&
          RelIsoPF((int)i, pt, ch, neu, pho, pu) < 0.15)
        sel.push_back((int)i);
    if (sel.empty()) continue;
    pf->GetEntry(ie);
    for (int i : sel)
    {
      double sum = 0, pfm = -1, best = 0.05;
      for (size_t k = 0; k < pfPt->size(); ++k)
      {
        const double dr = DeltaR(pfEta->at(k), pfPhi->at(k), eta->at(i), phi->at(i));
        if (dr < 0.1) sum += pfPt->at(k);
        if (pfId->at(k) == 3 && dr < best) { best = dr; pfm = pfPt->at(k); }
      }
      const int r = RegionOf(eta->at(i));
      val[g][r][0].push_back(sum / pt->at(i));
      if (pfm > 0) val[g][r][1].push_back(pfm / pt->at(i));
      else ++nNoPFmu[g][r];
      ++nMu;
    }
  }
  std::cout << "[LOOP] " << fname << ": " << nMu << " analysis muons\n";
  f->Close();
}

double Median(std::vector<double> v)
{
  if (v.empty()) return 0;
  std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
  return v[v.size() / 2];
}
// the fraction of values within +-0.1% of 1 (an exact copy of muPt)
double FracExact(const std::vector<double> &v)
{
  if (v.empty()) return 0;
  long long k = 0;
  for (double x : v) if (std::fabs(x - 1.0) < 1e-3) ++k;
  return double(k) / v.size();
}
} // namespace PRCM

void pf_raw_check_mu()
{
  using namespace PRCM;
  Loop(kDefaultDataFile, 0, -1);
  Loop("/Users/zhenghuang/pO_2026_Aug_20/Aug_20_MC_Wp_mu.root", 1, 300000);
  std::cout << "\nMEDIANS (data | MC), analysis muons; in brackets the fraction within 0.1% of 1\n";
  std::cout << Form("  %-10s", "");
  for (int r = 0; r < kNReg; ++r) std::cout << Form(" | %-34s", kReg[r]);
  std::cout << "\n";
  for (int v = 0; v < kNV; ++v)
  {
    std::cout << Form("  %-10s", kV[v]);
    for (int r = 0; r < kNReg; ++r)
      std::cout << Form(" | %.4f [%.2f] | %.4f [%.2f] (N %zu)", Median(val[0][r][v]), FracExact(val[0][r][v]), Median(val[1][r][v]),
                        FracExact(val[1][r][v]), val[0][r][v].size());
    std::cout << "\n";
  }
  for (int r = 0; r < kNReg; ++r)
    std::cout << Form("  %s: muons without a PF muon within DR < 0.05: data %lld, MC %lld\n", kReg[r], nNoPFmu[0][r], nNoPFmu[1][r]);
}
