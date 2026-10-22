// SCRATCH STUDY (2026-10-20): does the PF electron's missing data-only energy
// calibration explain the electron W deficit at high MET / m_T?
// ep_charge.C found PF pT (DR < 0.1) / elePt = 0.916-0.92 (|eta| > 1.566) and 0.98
// (barrel) in DATA but 1.01 in MC: the PF electron (hence the PF MET) carries the
// uncorrected energy, elePt the forest's data-only correction (elePt/eleRawPt =
// 1.094 endcap / 1.021 barrel). Here the MET is recomputed with the PF candidates
// within DR < 0.1 of the leading electron replaced by the electron itself:
//   MET_corr = MET_PF + (sum_PF(DR<0.1) pT - elePt) * e_T(electron)
// (data and MC alike) and the W selection (step 8) is followed by m_T > 40 and
// MET >= 30 with the PF and with the corrected MET, per charge and region:
// edges |eta| > 2.0 (bwd + fwd), mid 1.566-2.0, cen |eta| < 1.4442.
// MC = W+ and W- -> e nu + DY + tau + ttbar (k_s, gen weight, no SFs).
// Run: root -l -b -q 'met_corr.C+'
#include "/Users/zhenghuang/pO_analysis/skim/skim_common.h"
#include "/Users/zhenghuang/pO_analysis/skim/mc_norm.h"

#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TVector2.h"

#include <iostream>
#include <string>
#include <vector>

using namespace pOSkim;

namespace MC2
{
const int kNReg = 3, kNQ = 2, kNStep = 5; // step 8, +mT>40 (PF), +mT>40 (corr), +MET>=30 (PF), +MET>=30 (corr)
const char *const kReg[kNReg] = {"edges |eta|>2.0", "mid 1.566-2.0", "cen |eta|<1.4442"};
const char *const kStep[kNStep] = {"step 8 (W sel.)", "+ m_T>40, PF MET", "+ m_T>40, corrected MET", "+ MET>=30, PF MET",
                                   "+ MET>=30, corrected MET"};
double D[kNReg][kNQ][kNStep], M[kNReg][kNQ][kNStep], MW[kNReg][kNQ][kNStep], M2[kNReg][kNQ][kNStep];
double shiftSum[2][kNReg][kNQ], shiftW[2][kNReg][kNQ]; // <|MET_corr| - |MET_PF|> per (data/mc)
const char *const kNgen = "/Users/zhenghuang/pO_analysis/skim/rootfile/ngen.root";

inline int RegionOf(double eta)
{
  const double a = std::fabs(eta);
  if (a > 2.0) return 0;
  if (a > kEcalGapHi) return 1;
  return 2;
}

int Loop(const std::string &fname, bool isMC, double ks, bool isW)
{
  TFile *f = TFile::Open(fname.c_str());
  TTree *t = (TTree *)f->Get("ggHiNtuplizer/EventTree"), *h = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  TTree *pf = (TTree *)f->Get("particleFlowAnalyser/pftree"), *hlt = (TTree *)f->Get("hltanalysis/HltTree");
  TTree *ev = (TTree *)f->Get("skimanalysis/HltTree"), *obj = (TTree *)f->Get("hltobject/HLT_OxyL1SingleEG10_v");
  std::vector<float> *pt = nullptr, *eta = nullptr, *phi = nullptr, *ch = nullptr, *neu = nullptr, *pho = nullptr, *pu = nullptr,
                     *scEta = nullptr;
  std::vector<int> *q = nullptr, *id = nullptr, *iso = nullptr;
  t->SetBranchStatus("*", 0);
  auto bind = [&](const char *b, auto *a) { t->SetBranchStatus(b, 1); t->SetBranchAddress(b, a); };
  bind("elePt", &pt); bind("eleEta", &eta); bind("elePhi", &phi); bind("eleCharge", &q); bind(kEleIDBranch, &id);
  bind(kEleIsoBranch, &iso); bind("elePFChIso", &ch); bind("elePFNeuIso", &neu); bind("elePFPhoIso", &pho);
  bind("elePFPUIso", &pu); bind("eleSCEta", &scEta);
  std::vector<float> *pfPt = nullptr, *pfEta = nullptr, *pfPhi = nullptr;
  pf->SetBranchStatus("*", 0);
  for (const char *b : {"pfPt", "pfEta", "pfPhi"}) pf->SetBranchStatus(b, 1);
  pf->SetBranchAddress("pfPt", &pfPt); pf->SetBranchAddress("pfEta", &pfEta); pf->SetBranchAddress("pfPhi", &pfPhi);
  Int_t bit = 0;
  hlt->SetBranchStatus("*", 0);
  hlt->SetBranchStatus("HLT_OxyL1SingleEG10_v1", 1); hlt->SetBranchAddress("HLT_OxyL1SingleEG10_v1", &bit);
  std::vector<double> *toPt = nullptr, *toEta = nullptr, *toPhi = nullptr;
  obj->SetBranchStatus("*", 0);
  for (const char *b : {"pt", "eta", "phi"}) obj->SetBranchStatus(b, 1);
  obj->SetBranchAddress("pt", &toPt); obj->SetBranchAddress("eta", &toEta); obj->SetBranchAddress("phi", &toPhi);
  Int_t ppv = 1, pcc = 1;
  ev->SetBranchStatus("*", 0);
  ev->SetBranchStatus("pprimaryVertexFilter", 1); ev->SetBranchAddress("pprimaryVertexFilter", &ppv);
  ev->SetBranchStatus("pclusterCompatibilityFilter", 1); ev->SetBranchAddress("pclusterCompatibilityFilter", &pcc);
  Float_t vz = 0, w0 = 1;
  h->SetBranchStatus("*", 0);
  h->SetBranchStatus("vz", 1); h->SetBranchAddress("vz", &vz);
  if (isMC) { h->SetBranchStatus("weight", 1); h->SetBranchAddress("weight", &w0); }
  bool warnedF = false, warnedT = false;
  long long nSel = 0;
  const Long64_t n = t->GetEntries();
  for (Long64_t ie = 0; ie < n; ++ie)
  {
    t->GetEntry(ie); h->GetEntry(ie); ev->GetEntry(ie); hlt->GetEntry(ie);
    if (!pt || !id || !iso || !scEta) continue;
    if (!PassEventSelection_pO(warnedF, true, ppv, true, pcc) || std::fabs(vz) > 15.0) continue;
    if (!TriggerFired(bit)) continue;
    const int ne = (int)pt->size();
    bool has25 = false;
    for (int i = 0; i < ne && !has25; ++i) has25 = pt->at(i) > 25.0;
    if (!has25) continue;
    auto idOk = [&](int i) { return id->at(i) != 0 && !InEcalGap(scEta->at(i)); };
    auto isoOk = [&](int i) { return iso->at(i) != 0 && RelIsoPF(i, pt, ch, neu, pho, pu) < kEleIsoRelMax; };
    std::vector<int> c;
    for (int i = 0; i < ne; ++i)
      if (pt->at(i) > 10.0 && idOk(i) && isoOk(i)) c.push_back(i);
    bool veto = false;
    for (size_t a = 0; a < c.size() && !veto; ++a)
      for (size_t b = a + 1; b < c.size() && !veto; ++b)
      {
        if (q->at(c[a]) * q->at(c[b]) >= 0) continue;
        TLorentzVector pa, pb;
        pa.SetPtEtaPhiM(pt->at(c[a]), eta->at(c[a]), phi->at(c[a]), ELE_MASS);
        pb.SetPtEtaPhiM(pt->at(c[b]), eta->at(c[b]), phi->at(c[b]), ELE_MASS);
        const double m = (pa + pb).M();
        veto = m > 80.0 && m < 110.0;
      }
    if (veto) continue;
    int iL = -1;
    for (int i = 0; i < ne; ++i)
      if (idOk(i) && (iL < 0 || pt->at(i) > pt->at(iL))) iL = i;
    if (iL < 0 || pt->at(iL) <= 25.0 || std::fabs(eta->at(iL)) > 2.4 || !isoOk(iL)) continue;
    obj->GetEntry(ie);
    if (!PassLeadingLeptonTrigMatch(0.4, iL, eta, phi, true, toPt, true, toEta, true, toPhi, warnedT, "electron")) continue;
    pf->GetEntry(ie);
    const TVector2 met = ComputePFMET(nullptr, pfPt, pfPhi);
    double sum = 0;
    for (size_t k = 0; k < pfPt->size(); ++k)
      if (DeltaR(pfEta->at(k), pfPhi->at(k), eta->at(iL), phi->at(iL)) < 0.1) sum += pfPt->at(k);
    const TVector2 e(std::cos(phi->at(iL)), std::sin(phi->at(iL)));
    const TVector2 metC = met + (sum - pt->at(iL)) * e;
    const double mt = TransverseMass(pt->at(iL), phi->at(iL), met), mtC = TransverseMass(pt->at(iL), phi->at(iL), metC);
    const int r = RegionOf(eta->at(iL)), qi = q->at(iL) > 0 ? 1 : 0;
    const double w = (isMC ? (double)w0 : 1.0) * ks;
    const bool pass[kNStep] = {true, mt > 40.0, mtC > 40.0, met.Mod() >= 30.0, metC.Mod() >= 30.0};
    for (int s = 0; s < kNStep; ++s)
      if (pass[s])
      {
        if (!isMC) D[r][qi][s] += w;
        else { M[r][qi][s] += w; M2[r][qi][s] += w * w; if (isW) MW[r][qi][s] += w; }
      }
    shiftSum[isMC][r][qi] += w * (metC.Mod() - met.Mod());
    shiftW[isMC][r][qi] += w;
    ++nSel;
  }
  std::cout << Form("[LOOP] %s: %lld W candidates (raw)\n", fname.c_str(), nSel);
  f->Close();
  return 0;
}
} // namespace MC2

void met_corr()
{
  using namespace MC2;
  const std::string b = "/Users/zhenghuang/pO_2026_Aug_20/Aug_20_MC_";
  Loop(kDefaultDataFile, false, 1.0, false);
  const std::vector<std::pair<std::string, const char *>> mc = {{"Wp_ele", "Wp_ele"}, {"Wm_ele", "Wm_ele"}, {"DY_ele_Z", "DYee"},
                                                               {"DY_tau_Z", "DYtau"}, {"Wp_tau", "Wp_tau"}, {"Wm_tau", "Wm_tau"},
                                                               {"TT", "TT"}};
  for (auto &s : mc) Loop(b + s.first + ".root", true, pONorm::MCScale(s.second, kNgen), s.first.rfind("W", 0) == 0 && s.first.find("tau") == std::string::npos);
  for (int r = 0; r < kNReg; ++r)
  {
    std::cout << Form("\n=== %s ===  <|MET_corr| - |MET_PF|>: data e- %+.2f e+ %+.2f GeV | MC e- %+.2f e+ %+.2f GeV\n", kReg[r],
                      shiftSum[0][r][0] / shiftW[0][r][0], shiftSum[0][r][1] / shiftW[0][r][1], shiftSum[1][r][0] / shiftW[1][r][0],
                      shiftSum[1][r][1] / shiftW[1][r][1]);
    std::cout << "  step                       |   D-     D+   |   M-       M+   | r- = D-/M-     r+ = D+/M+     | r-/r+\n";
    for (int s = 0; s < kNStep; ++s)
    {
      double rr[2], re[2];
      for (int qi = 0; qi < 2; ++qi)
      {
        rr[qi] = D[r][qi][s] / M[r][qi][s];
        re[qi] = rr[qi] * std::sqrt(1.0 / std::max(1.0, D[r][qi][s]) + M2[r][qi][s] / (M[r][qi][s] * M[r][qi][s]));
      }
      const double dr = rr[0] / rr[1], de = dr * std::hypot(re[0] / rr[0], re[1] / rr[1]);
      std::cout << Form("  %-26s | %5.0f %5.0f | %7.1f %7.1f | %.3f +- %.3f  %.3f +- %.3f  | %.3f +- %.3f\n", kStep[s], D[r][0][s],
                        D[r][1][s], M[r][0][s], M[r][1][s], rr[0], re[0], rr[1], re[1], dr, de);
    }
    // the symmetric-QCD prediction of the e- count from the e+ side (r_W = 1.2 x <SF> 0.887 on the MC)
    const double rW = 1.2 * 0.887;
    for (int s : {1, 2, 3, 4})
    {
      const double Q = D[r][1][s] - rW * M[r][1][s];
      const double pred = rW * M[r][0][s] + Q;
      std::cout << Form("  [sym-QCD] %-26s: e- predicted %.0f (QCD per charge %.0f) vs observed %.0f -> ratio %.3f +- %.3f\n",
                        kStep[s], pred, Q, D[r][0][s], D[r][0][s] / pred,
                        D[r][0][s] / pred * std::sqrt(1.0 / std::max(1.0, D[r][0][s]) + D[r][1][s] / (pred * pred)));
    }
  }
}
