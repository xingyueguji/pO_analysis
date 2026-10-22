// SCRATCH (2026-10-20): which electron energy does the PF electron carry?
// For every analysis electron (WP90 ID outside the crack, WP90 iso && relIso < 0.3,
// pT > 25, |eta| < 2.4; event filter, |vz| < 15, EG10 fired), per region of its track
// eta, data vs W -> e nu MC, the MEDIAN of
//   PF / elePt, PF / eleRawPt, elePt / eleRawPt, PF / (eleSCEn / cosh eta_SC), PF / (eleEn / cosh eta)
// with PF = the sum of pT of the PF candidates within DR < 0.1 of the electron, and
// the median of the PF electron candidate (id 2) closest to the electron (DR < 0.05)
// over eleRawPt and over elePt. If the PF electron carries the uncorrected energy,
// PF / eleRawPt ~ 1 in data and MC alike.
// Run from correction/edge_study: root -l -b -q pf_raw_check.C+ | tee pf_raw_check.log
#include "/Users/zhenghuang/pO_analysis/skim/skim_common.h"

#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

using namespace pOSkim;

namespace PRC
{
const int kNReg = 3, kNV = 7;
const char *const kReg[kNReg] = {"|eta|>2.0", "1.566-2.0", "barrel"};
const char *const kV[kNV] = {"PF/elePt", "PF/eleRawPt", "elePt/eleRawPt", "PF/E_SC_T", "PF/eleEn_T", "PFe/eleRawPt", "PFe/elePt"};
std::vector<double> val[2][kNReg][kNV];

inline int RegionOf(double eta)
{
  const double a = std::fabs(eta);
  return a > 2.0 ? 0 : (a > kEcalGapHi ? 1 : 2);
}

void Loop(const std::string &fname, int g, Long64_t nMax)
{
  TFile *f = TFile::Open(fname.c_str());
  TTree *t = (TTree *)f->Get("ggHiNtuplizer/EventTree"), *h = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  TTree *pf = (TTree *)f->Get("particleFlowAnalyser/pftree"), *hlt = (TTree *)f->Get("hltanalysis/HltTree");
  TTree *ev = (TTree *)f->Get("skimanalysis/HltTree");
  std::vector<float> *pt = nullptr, *raw = nullptr, *en = nullptr, *scEn = nullptr, *eta = nullptr, *phi = nullptr, *scEta = nullptr,
                     *ch = nullptr, *neu = nullptr, *pho = nullptr, *pu = nullptr;
  std::vector<int> *id = nullptr, *iso = nullptr;
  t->SetBranchStatus("*", 0);
  auto bind = [&](const char *b, auto *a) { t->SetBranchStatus(b, 1); t->SetBranchAddress(b, a); };
  bind("elePt", &pt); bind("eleRawPt", &raw); bind("eleEn", &en); bind("eleSCEn", &scEn); bind("eleEta", &eta);
  bind("elePhi", &phi); bind("eleSCEta", &scEta); bind(kEleIDBranch, &id); bind(kEleIsoBranch, &iso);
  bind("elePFChIso", &ch); bind("elePFNeuIso", &neu); bind("elePFPhoIso", &pho); bind("elePFPUIso", &pu);
  std::vector<int> *pfId = nullptr;
  std::vector<float> *pfPt = nullptr, *pfEta = nullptr, *pfPhi = nullptr;
  pf->SetBranchStatus("*", 0);
  for (const char *b : {"pfId", "pfPt", "pfEta", "pfPhi"}) pf->SetBranchStatus(b, 1);
  pf->SetBranchAddress("pfId", &pfId); pf->SetBranchAddress("pfPt", &pfPt);
  pf->SetBranchAddress("pfEta", &pfEta); pf->SetBranchAddress("pfPhi", &pfPhi);
  Int_t bit = 0, ppv = 1, pcc = 1;
  hlt->SetBranchStatus("*", 0);
  hlt->SetBranchStatus("HLT_OxyL1SingleEG10_v1", 1); hlt->SetBranchAddress("HLT_OxyL1SingleEG10_v1", &bit);
  ev->SetBranchStatus("*", 0);
  ev->SetBranchStatus("pprimaryVertexFilter", 1); ev->SetBranchAddress("pprimaryVertexFilter", &ppv);
  ev->SetBranchStatus("pclusterCompatibilityFilter", 1); ev->SetBranchAddress("pclusterCompatibilityFilter", &pcc);
  Float_t vz = 0;
  h->SetBranchStatus("*", 0);
  h->SetBranchStatus("vz", 1); h->SetBranchAddress("vz", &vz);
  bool warned = false;
  long long nEle = 0;
  const Long64_t n = (nMax > 0) ? std::min(nMax, t->GetEntries()) : t->GetEntries();
  for (Long64_t ie = 0; ie < n; ++ie)
  {
    t->GetEntry(ie); h->GetEntry(ie); ev->GetEntry(ie); hlt->GetEntry(ie);
    if (!pt || !raw || !id) continue;
    if (!PassEventSelection_pO(warned, true, ppv, true, pcc) || std::fabs(vz) > 15.0 || !TriggerFired(bit)) continue;
    std::vector<int> sel;
    for (size_t i = 0; i < pt->size(); ++i)
      if (pt->at(i) > 25.0 && std::fabs(eta->at(i)) < 2.4 && id->at(i) != 0 && !InEcalGap(scEta->at(i)) && iso->at(i) != 0 &&
          RelIsoPF((int)i, pt, ch, neu, pho, pu) < kEleIsoRelMax)
        sel.push_back((int)i);
    if (sel.empty()) continue;
    pf->GetEntry(ie);
    for (int i : sel)
    {
      double sum = 0, pfe = -1, best = 0.05;
      for (size_t k = 0; k < pfPt->size(); ++k)
      {
        const double dr = DeltaR(pfEta->at(k), pfPhi->at(k), eta->at(i), phi->at(i));
        if (dr < 0.1) sum += pfPt->at(k);
        if (pfId->at(k) == 2 && dr < best) { best = dr; pfe = pfPt->at(k); }
      }
      const int r = RegionOf(eta->at(i));
      const double x[kNV] = {sum / pt->at(i), sum / raw->at(i), pt->at(i) / raw->at(i),
                             sum / (scEn->at(i) / std::cosh(scEta->at(i))), sum / (en->at(i) / std::cosh(eta->at(i))),
                             pfe > 0 ? pfe / raw->at(i) : -1, pfe > 0 ? pfe / pt->at(i) : -1};
      for (int v = 0; v < kNV; ++v)
        if (x[v] > 0) val[g][r][v].push_back(x[v]);
      ++nEle;
    }
  }
  std::cout << "[LOOP] " << fname << ": " << nEle << " analysis electrons\n";
  f->Close();
}

double Median(std::vector<double> v)
{
  if (v.empty()) return 0;
  std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
  return v[v.size() / 2];
}
} // namespace PRC

void pf_raw_check()
{
  using namespace PRC;
  Loop(kDefaultDataFile, 0, -1);
  Loop("/Users/zhenghuang/pO_2026_Aug_20/Aug_20_MC_Wp_ele.root", 1, 300000);
  std::cout << "\nMEDIANS (data | MC), analysis electrons. PF = sum pT of PF candidates in DR < 0.1; PFe = the matched PF electron (id 2, DR < 0.05);\n"
               "E_SC_T = eleSCEn / cosh(eta_SC); eleEn_T = eleEn / cosh(eta)\n";
  std::cout << Form("  %-16s", "");
  for (int r = 0; r < kNReg; ++r) std::cout << Form(" | %-21s", kReg[r]);
  std::cout << "\n";
  for (int v = 0; v < kNV; ++v)
  {
    std::cout << Form("  %-16s", kV[v]);
    for (int r = 0; r < kNReg; ++r) std::cout << Form(" | %.3f | %.3f (N %zu)", Median(val[0][r][v]), Median(val[1][r][v]), val[0][r][v].size());
    std::cout << "\n";
  }
}
