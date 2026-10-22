// SCRATCH STUDY (2026-10-20): is there a CHARGE-DEPENDENT electron momentum or MET
// bias at the lab edges in data? For every analysis electron (WP90 ID outside the
// crack, WP90 iso && relIso < 0.3, |eta| < 2.4), per charge and region of its track eta
//   bwd  eta > 2.0 (= lab y0) | fwd eta < -2.0 (= y11) | mid 1.566 < |eta| <= 2.0 | cen |eta| < 1.4442
// in three selections -- Z leg (Z->ee: two such electrons pT > 20, OS, 80 < m < 100,
// EG10 fired), W (the skim's W selection, step 8, leading electron pT > 25), W with
// PF MET >= 30 -- data vs MC (DY->ee for Z, W+ and W- -> e nu for W; gen weight x k_s):
//   EoP   = eleEoverP (CMSSW eSuperClusterOverP)
//   rTrk  = eleTrkPt / elePt
//   rPF   = (sum pT of the PF candidates within DR < 0.1 of the electron) / elePt
//   fPFe  = fraction whose leading PF candidate in DR < 0.1 is a PF electron (id 2)
//   uPar  = the recoil projected on the electron: u = -MET - (sum of the selected
//           electrons' pT vectors), uPar = u . e_T(electron)   [Z: per leg; W: the lead]
// The test is the CHARGE DIFFERENCE <x>(e-) - <x>(e+) in data minus the same in MC:
// energy-scale offsets common to both charges cancel. Truncated means (windows
// below) so the brem / fake tails do not dominate.
// Run: root -l -b -q 'ep_charge.C+("loop")' then 'ep_charge.C+("report")'
#include "/Users/zhenghuang/pO_analysis/skim/skim_common.h"
#include "/Users/zhenghuang/pO_analysis/skim/mc_norm.h"

#include "TFile.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TVector2.h"

#include <iostream>
#include <string>
#include <vector>

using namespace pOSkim;

namespace EP
{
const int kNG = 2, kNSel = 3, kNReg = 4, kNQ = 2, kNVar = 5;
const char *const kG[kNG]     = {"data", "mc"};
const char *const kSel[kNSel] = {"Zleg", "W", "W_met30"};
const char *const kReg[kNReg] = {"bwd(y0)", "fwd(y11)", "mid", "cen"};
const char *const kVar[kNVar] = {"EoP", "rTrk", "rPF", "fPFe", "uPar"};
const double kLo[kNVar] = {0.5, 0.5, 0.5, -0.5, -80.0}, kHi[kNVar] = {2.0, 1.5, 1.5, 1.5, 80.0};
const char *const kNgen = "/Users/zhenghuang/pO_analysis/skim/rootfile/ngen.root";
const char *const kOut = "ep_charge.root";

inline int RegionOf(double eta)
{
  const double a = std::fabs(eta);
  if (eta > 2.0) return 0;
  if (eta < -2.0) return 1;
  if (a > kEcalGapHi) return 2;
  return 3;
}
// [g][sel][reg][q][var] -> sum w, sum w x, sum w x^2, sum w^2 (within the window)
double S[kNG][kNSel][kNReg][kNQ][kNVar][4];
void Acc(int g, int sel, int reg, int q, int v, double x, double w)
{
  if (x < kLo[v] || x > kHi[v]) return;
  double *s = S[g][sel][reg][q][v];
  s[0] += w; s[1] += w * x; s[2] += w * x * x; s[3] += w * w;
}

int Loop(const std::string &fname, int g, double ks, bool doZ, bool doW)
{
  TFile *f = TFile::Open(fname.c_str());
  TTree *t = (TTree *)f->Get("ggHiNtuplizer/EventTree"), *h = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  TTree *pf = (TTree *)f->Get("particleFlowAnalyser/pftree"), *hlt = (TTree *)f->Get("hltanalysis/HltTree");
  TTree *ev = (TTree *)f->Get("skimanalysis/HltTree"), *obj = (TTree *)f->Get("hltobject/HLT_OxyL1SingleEG10_v");
  if (!t || !h || !pf || !hlt || !ev || !obj) { std::cerr << "[FATAL] trees missing in " << fname << "\n"; return 2; }
  std::vector<float> *pt = nullptr, *eta = nullptr, *phi = nullptr, *ch = nullptr, *neu = nullptr, *pho = nullptr, *pu = nullptr,
                     *scEta = nullptr, *eop = nullptr, *trkPt = nullptr;
  std::vector<int> *q = nullptr, *id = nullptr, *iso = nullptr;
  t->SetBranchStatus("*", 0);
  auto bind = [&](const char *b, auto *a) { t->SetBranchStatus(b, 1); t->SetBranchAddress(b, a); };
  bind("elePt", &pt); bind("eleEta", &eta); bind("elePhi", &phi); bind("eleCharge", &q); bind(kEleIDBranch, &id);
  bind(kEleIsoBranch, &iso); bind("elePFChIso", &ch); bind("elePFNeuIso", &neu); bind("elePFPhoIso", &pho);
  bind("elePFPUIso", &pu); bind("eleSCEta", &scEta); bind("eleEoverP", &eop); bind("eleTrkPt", &trkPt);
  std::vector<int> *pfId = nullptr;
  std::vector<float> *pfPt = nullptr, *pfEta = nullptr, *pfPhi = nullptr;
  pf->SetBranchStatus("*", 0);
  for (const char *b : {"pfId", "pfPt", "pfEta", "pfPhi"}) pf->SetBranchStatus(b, 1);
  pf->SetBranchAddress("pfId", &pfId); pf->SetBranchAddress("pfPt", &pfPt);
  pf->SetBranchAddress("pfEta", &pfEta); pf->SetBranchAddress("pfPhi", &pfPhi);
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
  if (g == 1) { h->SetBranchStatus("weight", 1); h->SetBranchAddress("weight", &w0); }
  bool warnedF = false, warnedT = false;
  long long nZ = 0, nW = 0;
  const Long64_t n = t->GetEntries();
  std::cout << "[LOOP] " << fname << " (" << n << ") k_s " << ks << " Z " << doZ << " W " << doW << "\n";
  for (Long64_t ie = 0; ie < n; ++ie)
  {
    t->GetEntry(ie); h->GetEntry(ie); ev->GetEntry(ie); hlt->GetEntry(ie);
    if (!pt || !id || !iso || !scEta || !eop || !trkPt) continue;
    if (!PassEventSelection_pO(warnedF, true, ppv, true, pcc) || std::fabs(vz) > 15.0) continue;
    if (!TriggerFired(bit)) continue;
    const double w = (g == 1 ? (double)w0 : 1.0) * ks;
    const int ne = (int)pt->size();
    auto idOk = [&](int i) { return id->at(i) != 0 && !InEcalGap(scEta->at(i)); };
    auto isoOk = [&](int i) { return iso->at(i) != 0 && RelIsoPF(i, pt, ch, neu, pho, pu) < kEleIsoRelMax; };
    bool pfRead = false;
    TVector2 met;
    auto readPF = [&]() { if (!pfRead) { pf->GetEntry(ie); met = ComputePFMET(nullptr, pfPt, pfPhi); pfRead = true; } };
    // per-electron quantities; uRec = the recoil vector of the event
    auto fillEle = [&](int sel, int i, const TVector2 &uRec)
    {
      readPF();
      const int r = RegionOf(eta->at(i)), qi = q->at(i) > 0 ? 1 : 0;
      Acc(g, sel, r, qi, 0, eop->at(i), w);
      Acc(g, sel, r, qi, 1, trkPt->at(i) / pt->at(i), w);
      double sum = 0, lead = 0;
      int leadId = -1;
      for (size_t k = 0; k < pfPt->size(); ++k)
      {
        if (DeltaR(pfEta->at(k), pfPhi->at(k), eta->at(i), phi->at(i)) >= 0.1) continue;
        sum += pfPt->at(k);
        if (pfPt->at(k) > lead) { lead = pfPt->at(k); leadId = pfId->at(k); }
      }
      Acc(g, sel, r, qi, 2, sum / pt->at(i), w);
      Acc(g, sel, r, qi, 3, leadId == 2 ? 1.0 : 0.0, w);
      const TVector2 e(std::cos(phi->at(i)), std::sin(phi->at(i)));
      Acc(g, sel, r, qi, 4, uRec * e, w);
    };
    // ---- Z->ee legs
    if (doZ)
    {
      int ba = -1, bb = -1;
      double bd = 1e9;
      for (int a = 0; a < ne; ++a)
      {
        if (pt->at(a) <= 20.0 || std::fabs(eta->at(a)) >= 2.4 || !idOk(a) || !isoOk(a)) continue;
        for (int b = a + 1; b < ne; ++b)
        {
          if (pt->at(b) <= 20.0 || std::fabs(eta->at(b)) >= 2.4 || !idOk(b) || !isoOk(b)) continue;
          if (q->at(a) * q->at(b) >= 0) continue;
          TLorentzVector pa, pb;
          pa.SetPtEtaPhiM(pt->at(a), eta->at(a), phi->at(a), ELE_MASS);
          pb.SetPtEtaPhiM(pt->at(b), eta->at(b), phi->at(b), ELE_MASS);
          const double m = (pa + pb).M();
          if (m > 80 && m < 100 && std::fabs(m - 91.19) < bd) { bd = std::fabs(m - 91.19); ba = a; bb = b; }
        }
      }
      if (ba >= 0)
      {
        readPF();
        TVector2 qT(pt->at(ba) * std::cos(phi->at(ba)) + pt->at(bb) * std::cos(phi->at(bb)),
                    pt->at(ba) * std::sin(phi->at(ba)) + pt->at(bb) * std::sin(phi->at(bb)));
        const TVector2 u = -1.0 * met - qT;
        fillEle(0, ba, u);
        fillEle(0, bb, u);
        ++nZ;
      }
    }
    // ---- W: the skim's W selection (step 8) on the leading ID'd electron
    if (doW)
    {
      bool has25 = false;
      for (int i = 0; i < ne && !has25; ++i) has25 = pt->at(i) > 25.0;
      if (!has25) continue;
      std::vector<int> c; // DY veto legs
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
      readPF();
      const TVector2 pe(pt->at(iL) * std::cos(phi->at(iL)), pt->at(iL) * std::sin(phi->at(iL)));
      const TVector2 u = -1.0 * met - pe;
      fillEle(1, iL, u);
      if (met.Mod() >= 30.0) fillEle(2, iL, u);
      ++nW;
    }
  }
  std::cout << Form("[LOOP]   Z pairs %lld, W candidates %lld (raw)\n", nZ, nW);
  f->Close();
  return 0;
}

int DoLoop()
{
  for (auto &a : S) for (auto &b : a) for (auto &c : b) for (auto &d : c) for (auto &e : d) for (double &x : e) x = 0;
  const std::string base = "/Users/zhenghuang/pO_2026_Aug_20/";
  if (Loop(kDefaultDataFile, 0, 1.0, true, true)) return 2;
  if (Loop(base + "Aug_20_MC_DY_ele_Z.root", 1, pONorm::MCScale("DYee", kNgen), true, false)) return 2;
  if (Loop(base + "Aug_20_MC_Wp_ele.root", 1, pONorm::MCScale("Wp_ele", kNgen), false, true)) return 2;
  if (Loop(base + "Aug_20_MC_Wm_ele.root", 1, pONorm::MCScale("Wm_ele", kNgen), false, true)) return 2;
  TFile *fo = TFile::Open(kOut, "RECREATE");
  const int nTot = kNG * kNSel * kNReg * kNQ * kNVar * 4;
  TH1D *hs = new TH1D("h_sums", "flattened sums [g][sel][reg][q][var][w,wx,wx2,w2]", nTot, 0, nTot);
  int k = 0;
  for (auto &a : S) for (auto &b : a) for (auto &c : b) for (auto &d : c) for (auto &e : d) for (double x : e) hs->SetBinContent(++k, x);
  fo->Write();
  fo->Close();
  std::cout << "[OUT] " << kOut << "\n";
  return 0;
}

struct M { double m = 0, e = 0, n = 0; };
M Mean(int g, int sel, int reg, int q, int v)
{
  const double *s = S[g][sel][reg][q][v];
  M r;
  if (s[0] <= 0) return r;
  r.m = s[1] / s[0];
  const double var = std::max(0.0, s[2] / s[0] - r.m * r.m);
  r.n = s[0] * s[0] / s[3]; // effective entries
  r.e = std::sqrt(var / r.n);
  return r;
}

int DoReport()
{
  TFile *fi = TFile::Open(kOut);
  if (!fi) return 2;
  TH1 *hs = (TH1 *)fi->Get("h_sums");
  int k = 0;
  for (auto &a : S) for (auto &b : a) for (auto &c : b) for (auto &d : c) for (auto &e : d) for (double &x : e) x = hs->GetBinContent(++k);
  const char *vlab[kNVar] = {"E/p (eleEoverP)", "track pT / electron pT", "PF pT (DR<0.1) / electron pT",
                             "lead PF cand. in DR<0.1 is a PF electron", "recoil along the electron [GeV]"};
  for (int sel = 0; sel < kNSel; ++sel)
  {
    std::cout << "\n==================== " << kSel[sel] << " ====================\n";
    for (int v = 0; v < kNVar; ++v)
    {
      std::cout << Form("--- %s (truncated to [%g, %g]) ---\n", vlab[v], kLo[v], kHi[v]);
      std::cout << "  region    |  data e-           data e+           Neff-/+   |  MC e-     MC e+    | (e- - e+)data - (e- - e+)MC\n";
      for (int r = 0; r < kNReg; ++r)
      {
        const M dm = Mean(0, sel, r, 0, v), dp = Mean(0, sel, r, 1, v), mm = Mean(1, sel, r, 0, v), mp = Mean(1, sel, r, 1, v);
        const double dd = (dm.m - dp.m) - (mm.m - mp.m);
        const double de = std::sqrt(dm.e * dm.e + dp.e * dp.e + mm.e * mm.e + mp.e * mp.e);
        std::cout << Form("  %-9s | %8.4f+-%.4f  %8.4f+-%.4f  %5.0f/%5.0f | %8.4f  %8.4f | %+8.4f +- %.4f (%+.1f sigma)\n", kReg[r],
                          dm.m, dm.e, dp.m, dp.e, dm.n, dp.n, mm.m, mp.m, dd, de, de > 0 ? dd / de : 0.0);
      }
    }
  }
  return 0;
}
} // namespace EP

int ep_charge(const char *mode = "loop")
{
  const std::string m = mode;
  if (m == "loop") return EP::DoLoop();
  if (m == "report") return EP::DoReport();
  return 1;
}
