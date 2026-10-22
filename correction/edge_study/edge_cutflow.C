// SCRATCH STUDY (2026-10-20): where does the W- -> e nu edge deficit appear?
// A per-charge cutflow at the lab edges, data vs MC, electrons with muons as the
// control, plus e-/e+ (mu-/mu+) vs phi and vs run in data.
//
// The selection replicates skim/skim.C's W selection (via correction/njet_WZ.C,
// event-identical with the skim): steps (1) a lepton pT > 25, (2) pO filter +
// |vz| < 15, (3) the flavour's trigger fired, (4) DY veto, then the leading
// lepton. Two chains on top of the event steps:
//   A1 leading ID'd lepton (e: MVA WP90 outside the crack; mu: Tight + PF), pT > 25, |eta| < 2.4
//   A2 + isolation (e: MVA iso WP90 && relIso < 0.3; mu: relIso < 0.15)
//   A3 + trigger match DR < 0.4                       = the skim's W selection (N[8])
//   A4 + m_T > 40                                     = the leppt_mt40 selection
//   A5 A3 + PF MET >= 30
//   B0 leading lepton WITHOUT ID (crack veto / PF), pT > 25, |eta| < 2.4, m_T > 40
//   B1 B0 + that lepton passes the ID                  (the ID step on a W-enriched sample)
// Regions of the selected lepton (track eta, y = -eta): "bwd" eta in (2.0, 2.4]
// (= lab bin y0), "fwd" eta in [-2.4, -2.0) (= y11), "cen" |eta| <= 2.0.
// MC = absolute (k_s = A sigma L / N_gen, gen weight), NO efficiency SFs: every
// quantity quoted is a charge ratio, in which the (charge-blind) SFs cancel.
// QCD is not subtracted: it is (close to) charge-symmetric, so the charge
// DIFFERENCE of data is QCD-free (see the analysis part).
//
// Run (from anywhere): root -l -b -q 'edge_cutflow.C+("loop")'  -> edge_cutflow.root
//                      root -l -b -q 'edge_cutflow.C+("report")' -> console + plots/ + csv
#include "/Users/zhenghuang/pO_analysis/skim/skim_common.h"
#include "/Users/zhenghuang/pO_analysis/skim/mc_norm.h"
#include "/Users/zhenghuang/pO_analysis/plotting/plotting_helper.C"

#include "TFile.h"
#include "TH1D.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector2.h"

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace pOSkim;

namespace EC
{
const int kNFl = 2, kNGrp = 3, kNStep = 7, kNQ = 2, kNReg = 3, kNPhi = 12;
const char *const kFl[kNFl]    = {"mu", "e"};
const char *const kGrp[kNGrp]  = {"data", "mcW", "mcOther"};
const char *const kStep[kNStep] = {"A1 lead ID'd lepton",  "A2 + iso",          "A3 + trig match (W sel.)", "A4 + m_T > 40",
                                   "A5 A3 + MET >= 30",    "B0 no-ID lead, m_T>40", "B1 B0 + passes ID"};
const char *const kQ[kNQ]      = {"minus", "plus"};
const char *const kReg[kNReg]  = {"bwd", "fwd", "cen"};
const char *const kOut = "edge_cutflow.root";
const char *const kNgen = "/Users/zhenghuang/pO_analysis/skim/rootfile/ngen.root";

inline int RegionOf(double eta) { return eta > 2.0 ? 0 : (eta < -2.0 ? 1 : 2); }
inline TString HName(int fl, int g, int r) { return Form("h_cf_%s_%s_%s", kFl[fl], kGrp[g], kReg[r]); }    // x = step*2 + q
inline TString PName(int fl, int g, int q, int r) { return Form("h_phi_%s_%s_%s_%s", kFl[fl], kGrp[g], kQ[q], kReg[r]); }

// one flavour's branches on the shared EventTree
struct Lep
{
  bool isMu = true;
  std::vector<float> *pt = nullptr, *eta = nullptr, *phi = nullptr, *ch = nullptr, *neu = nullptr, *pho = nullptr, *pu = nullptr,
                     *scEta = nullptr;
  std::vector<int> *q = nullptr, *id = nullptr, *isPF = nullptr, *mvaIso = nullptr;
  Int_t hlt = 0;
  std::vector<double> *toPt = nullptr, *toEta = nullptr, *toPhi = nullptr;
  TTree *tObj = nullptr;
  bool warnedTrig = false;

  bool Bind(TTree *t, TTree *tHLT, TFile *f)
  {
    auto bind = [&](const char *b, auto *addr) -> bool
    {
      if (!HasBranch(t, b)) { std::cerr << "[FATAL] missing " << b << "\n"; return false; }
      t->SetBranchStatus(b, 1);
      t->SetBranchAddress(b, addr);
      return true;
    };
    bool ok = true;
    if (isMu)
    {
      ok &= bind("muPt", &pt); ok &= bind("muEta", &eta); ok &= bind("muPhi", &phi); ok &= bind("muCharge", &q);
      ok &= bind("muIDTight", &id); ok &= bind("muPFChIso", &ch); ok &= bind("muPFNeuIso", &neu);
      ok &= bind("muPFPhoIso", &pho); ok &= bind("muPFPUIso", &pu); ok &= bind("muIsPF", &isPF);
    }
    else
    {
      ok &= bind("elePt", &pt); ok &= bind("eleEta", &eta); ok &= bind("elePhi", &phi); ok &= bind("eleCharge", &q);
      ok &= bind(kEleIDBranch, &id); ok &= bind("elePFChIso", &ch); ok &= bind("elePFNeuIso", &neu);
      ok &= bind("elePFPhoIso", &pho); ok &= bind("elePFPUIso", &pu); ok &= bind("eleSCEta", &scEta);
      ok &= bind(kEleIsoBranch, &mvaIso);
    }
    const std::string hn = isMu ? "HLT_OxyL1SingleMuOpen_v1" : "HLT_OxyL1SingleEG10_v1";
    if (!HasBranch(tHLT, hn.c_str())) { std::cerr << "[FATAL] no " << hn << "\n"; return false; }
    tHLT->SetBranchStatus(hn.c_str(), 1);
    tHLT->SetBranchAddress(hn.c_str(), &hlt);
    tObj = (TTree *)f->Get(isMu ? "hltobject/HLT_OxyL1SingleMuOpen_v" : "hltobject/HLT_OxyL1SingleEG10_v");
    if (!tObj) { std::cerr << "[FATAL] no trigger-object tree\n"; return false; }
    tObj->SetBranchStatus("*", 0);
    for (const char *b : {"pt", "eta", "phi"}) tObj->SetBranchStatus(b, 1);
    tObj->SetBranchAddress("pt", &toPt); tObj->SetBranchAddress("eta", &toEta); tObj->SetBranchAddress("phi", &toPhi);
    return ok;
  }
  bool Valid() const { return pt && eta && phi && q && id && (isMu ? (isPF != nullptr) : (scEta && mvaIso)); }
  bool Base(int i) const // the object requirements before ID: e outside the crack, mu PF
  {
    return isMu ? isPF->at(i) != 0 : !InEcalGap(scEta->at(i));
  }
  bool Iso(int i) const
  {
    if (RelIsoPF(i, pt, ch, neu, pho, pu) >= (isMu ? 0.15 : kEleIsoRelMax)) return false;
    return isMu || mvaIso->at(i) != 0;
  }
};

// one sample: fill the cutflow / phi / run histograms of the enabled flavours
int Loop(const std::string &fname, bool doMu, bool doEle, int grp, const char *labelMu, const char *labelEle, TFile *fo,
         std::map<unsigned, std::vector<double>> &runCnt)
{
  TFile *f = TFile::Open(fname.c_str());
  if (!f || f->IsZombie()) { std::cerr << "[FATAL] cannot open " << fname << "\n"; return 2; }
  TTree *tLep = (TTree *)f->Get("ggHiNtuplizer/EventTree"), *tHi = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  TTree *tPF = (TTree *)f->Get("particleFlowAnalyser/pftree"), *tHLT = (TTree *)f->Get("hltanalysis/HltTree");
  TTree *tEvt = (TTree *)f->Get("skimanalysis/HltTree");
  if (!tLep || !tHi || !tPF || !tHLT || !tEvt) { std::cerr << "[FATAL] missing tree in " << fname << "\n"; return 2; }
  tLep->SetBranchStatus("*", 0);
  tHLT->SetBranchStatus("*", 0);
  Lep L[2];
  L[0].isMu = true; L[1].isMu = false;
  const bool on[2] = {doMu, doEle};
  for (int fl = 0; fl < 2; ++fl)
    if (on[fl] && !L[fl].Bind(tLep, tHLT, f)) return 2;
  Int_t ppv = 1, pcc = 1;
  const bool has_ppv = HasBranch(tEvt, "pprimaryVertexFilter"), has_pcc = HasBranch(tEvt, "pclusterCompatibilityFilter");
  tEvt->SetBranchStatus("*", 0);
  if (has_ppv && has_pcc)
  {
    tEvt->SetBranchStatus("pprimaryVertexFilter", 1); tEvt->SetBranchAddress("pprimaryVertexFilter", &ppv);
    tEvt->SetBranchStatus("pclusterCompatibilityFilter", 1); tEvt->SetBranchAddress("pclusterCompatibilityFilter", &pcc);
  }
  Float_t vz = 0.f, w0 = 1.f;
  UInt_t run = 0;
  tHi->SetBranchStatus("*", 0);
  tHi->SetBranchStatus("vz", 1); tHi->SetBranchAddress("vz", &vz);
  tHi->SetBranchStatus("run", 1); tHi->SetBranchAddress("run", &run);
  const bool isMC = (grp != 0);
  if (isMC) { tHi->SetBranchStatus("weight", 1); tHi->SetBranchAddress("weight", &w0); }
  std::vector<float> *pfPt = nullptr, *pfPhi = nullptr;
  tPF->SetBranchStatus("*", 0);
  tPF->SetBranchStatus("pfPt", 1); tPF->SetBranchAddress("pfPt", &pfPt);
  tPF->SetBranchStatus("pfPhi", 1); tPF->SetBranchAddress("pfPhi", &pfPhi);

  double ks[2] = {1.0, 1.0};
  if (isMC)
  {
    if (doMu) ks[0] = pONorm::MCScale(labelMu, kNgen);
    if (doEle) ks[1] = pONorm::MCScale(labelEle, kNgen);
  }
  TH1D *hcf[2][kNReg], *hphi[2][kNQ][kNReg];
  for (int fl = 0; fl < 2; ++fl)
    for (int r = 0; r < kNReg; ++r)
    {
      hcf[fl][r] = (TH1D *)fo->Get(HName(fl, grp, r));
      for (int q = 0; q < kNQ; ++q) hphi[fl][q][r] = (TH1D *)fo->Get(PName(fl, grp, q, r));
    }
  long long nA3[2] = {0, 0}, nA4[2] = {0, 0};
  bool warnedFilters = false;
  const Long64_t n = tLep->GetEntries();
  std::cout << "[LOOP] " << fname << " (" << n << " entries) mu " << doMu << " e " << doEle << " k_s " << ks[0] << " / " << ks[1]
            << "\n";
  for (Long64_t ie = 0; ie < n; ++ie)
  {
    tLep->GetEntry(ie);
    tHi->GetEntry(ie);
    tEvt->GetEntry(ie);
    tHLT->GetEntry(ie);
    if (!PassEventSelection_pO(warnedFilters, has_ppv, ppv, has_pcc, pcc)) continue;
    if (std::fabs(vz) > 15.0) continue;
    bool pfRead = false;
    double metX = 0, metY = 0;
    auto met = [&]() -> TVector2
    {
      if (!pfRead) { tPF->GetEntry(ie); const TVector2 m = ComputePFMET(nullptr, pfPt, pfPhi); metX = m.X(); metY = m.Y(); pfRead = true; }
      return TVector2(metX, metY);
    };
    for (int fl = 0; fl < 2; ++fl)
    {
      if (!on[fl]) continue;
      Lep &l = L[fl];
      if (!l.Valid()) continue;
      const double w = (isMC ? (double)w0 : 1.0) * ks[fl];
      const int nl = (int)l.pt->size();
      // (1) a lepton pT > 25 (mu: PF)
      bool has25 = false;
      for (int i = 0; i < nl && !has25; ++i) has25 = l.pt->at(i) > 25.0 && (!l.isMu || l.isPF->at(i) != 0);
      if (!has25) continue;
      // (3) trigger
      if (!TriggerFired(l.hlt)) continue;
      // (4) DY veto
      {
        std::vector<int> c;
        for (int i = 0; i < nl; ++i)
          if (l.pt->at(i) > (l.isMu ? 15.0 : 10.0) && l.id->at(i) != 0 && l.Base(i) && l.Iso(i)) c.push_back(i);
        bool veto = false;
        for (size_t a = 0; a < c.size() && !veto; ++a)
          for (size_t b = a + 1; b < c.size() && !veto; ++b)
          {
            if (l.q->at(c[a]) * l.q->at(c[b]) >= 0) continue;
            TLorentzVector p1, p2;
            const double m = l.isMu ? MU_MASS : ELE_MASS;
            p1.SetPtEtaPhiM(l.pt->at(c[a]), l.eta->at(c[a]), l.phi->at(c[a]), m);
            p2.SetPtEtaPhiM(l.pt->at(c[b]), l.eta->at(c[b]), l.phi->at(c[b]), m);
            const double mll = (p1 + p2).M();
            veto = mll > 80.0 && mll < 110.0;
          }
        if (veto) continue;
      }
      auto fill = [&](int step, int i)
      {
        const int qi = l.q->at(i) > 0 ? 1 : 0, r = RegionOf(l.eta->at(i));
        hcf[fl][r]->Fill(2 * step + qi + 0.5, w);
      };
      // chain A: the leading ID'd lepton (the skim's)
      int iA = -1;
      for (int i = 0; i < nl; ++i)
        if (l.id->at(i) != 0 && l.Base(i) && (iA < 0 || l.pt->at(i) > l.pt->at(iA))) iA = i;
      if (iA >= 0 && l.pt->at(iA) > 25.0 && std::fabs(l.eta->at(iA)) <= 2.4)
      {
        fill(0, iA);
        if (l.Iso(iA))
        {
          fill(1, iA);
          l.tObj->GetEntry(ie);
          if (PassLeadingLeptonTrigMatch(0.4, iA, l.eta, l.phi, true, l.toPt, true, l.toEta, true, l.toPhi, l.warnedTrig,
                                         l.isMu ? "muon" : "electron"))
          {
            fill(2, iA);
            ++nA3[fl];
            const TVector2 m = met();
            const double mt = TransverseMass(l.pt->at(iA), l.phi->at(iA), m);
            if (mt > 40.0)
            {
              fill(3, iA);
              ++nA4[fl];
              const int qi = l.q->at(iA) > 0 ? 1 : 0, r = RegionOf(l.eta->at(iA));
              hphi[fl][qi][r]->Fill(l.phi->at(iA), w);
              if (!isMC)
              {
                std::vector<double> &v = runCnt[run];
                if (v.empty()) v.assign(kNFl * kNQ * kNReg, 0.0);
                v[(fl * kNQ + qi) * kNReg + r] += 1.0;
              }
            }
            if (m.Mod() >= 30.0) fill(4, iA);
          }
        }
      }
      // chain B: the leading lepton before the ID, m_T > 40
      int iB = -1;
      for (int i = 0; i < nl; ++i)
        if (l.Base(i) && (iB < 0 || l.pt->at(i) > l.pt->at(iB))) iB = i;
      if (iB >= 0 && l.pt->at(iB) > 25.0 && std::fabs(l.eta->at(iB)) <= 2.4)
      {
        const double mt = TransverseMass(l.pt->at(iB), l.phi->at(iB), met());
        if (mt > 40.0)
        {
          fill(5, iB);
          if (l.id->at(iB) != 0) fill(6, iB);
        }
      }
    }
  }
  for (int fl = 0; fl < 2; ++fl)
    if (on[fl]) std::cout << Form("[LOOP]   %s: A3 (W selection) %lld events, A4 (+ m_T > 40) %lld (raw)\n", kFl[fl], nA3[fl], nA4[fl]);
  f->Close();
  return 0;
}

int DoLoop()
{
  TFile *fo = TFile::Open(kOut, "RECREATE");
  for (int fl = 0; fl < kNFl; ++fl)
    for (int g = 0; g < kNGrp; ++g)
      for (int r = 0; r < kNReg; ++r)
      {
        TH1D *h = new TH1D(HName(fl, g, r), "cutflow: bin = 2*step + (q>0)", 2 * kNStep, 0, 2 * kNStep);
        h->Sumw2();
        for (int q = 0; q < kNQ; ++q)
        {
          TH1D *p = new TH1D(PName(fl, g, q, r), ";#phi;events", kNPhi, -TMath::Pi(), TMath::Pi());
          p->Sumw2();
        }
      }
  std::map<unsigned, std::vector<double>> runCnt;
  const std::string base = "/Users/zhenghuang/pO_2026_Aug_20/";
  struct S { std::string file; bool mu, ele; int grp; const char *lm, *le; };
  const std::vector<S> samples = {
      {kDefaultDataFile, true, true, 0, "", ""},
      {base + "Aug_20_MC_Wp_mu.root", true, false, 1, "Wp_mu", ""},
      {base + "Aug_20_MC_Wm_mu.root", true, false, 1, "Wm_mu", ""},
      {base + "Aug_20_MC_Wp_ele.root", false, true, 1, "", "Wp_ele"},
      {base + "Aug_20_MC_Wm_ele.root", false, true, 1, "", "Wm_ele"},
      {base + "Aug_20_MC_DY_mu_Z.root", true, false, 2, "DYmu", ""},
      {base + "Aug_20_MC_DY_ele_Z.root", false, true, 2, "", "DYee"},
      {base + "Aug_20_MC_DY_tau_Z.root", true, true, 2, "DYtau", "DYtau"},
      {base + "Aug_20_MC_Wp_tau.root", true, true, 2, "Wp_tau", "Wp_tau"},
      {base + "Aug_20_MC_Wm_tau.root", true, true, 2, "Wm_tau", "Wm_tau"},
      {base + "Aug_20_MC_TT.root", true, true, 2, "TT", "TT"},
  };
  for (const S &s : samples)
    if (Loop(s.file, s.mu, s.ele, s.grp, s.lm, s.le, fo, runCnt)) return 2;
  // runs: one TH1D per (flavour, charge, region), bin = run index; labels = run numbers
  fo->cd();
  const int nRun = (int)runCnt.size();
  for (int fl = 0; fl < kNFl; ++fl)
    for (int q = 0; q < kNQ; ++q)
      for (int r = 0; r < kNReg; ++r)
      {
        TH1D *h = new TH1D(Form("h_run_%s_%s_%s", kFl[fl], kQ[q], kReg[r]), "A4 data per run", nRun, 0, nRun);
        int k = 0;
        for (const auto &kv : runCnt)
        {
          h->GetXaxis()->SetBinLabel(k + 1, Form("%u", kv.first));
          h->SetBinContent(k + 1, kv.second[(fl * kNQ + q) * kNReg + r]);
          h->SetBinError(k + 1, std::sqrt(kv.second[(fl * kNQ + q) * kNReg + r]));
          ++k;
        }
      }
  fo->Write();
  fo->Close();
  std::cout << "[OUT] " << kOut << "\n";
  return 0;
}

// ---------------------------------------------------------------- report
struct V { double v = 0, e2 = 0; };
inline V Get(TFile *f, int fl, int g, int r, int step, int q)
{
  TH1 *h = (TH1 *)f->Get(HName(fl, g, r));
  V x;
  x.v = h->GetBinContent(2 * step + q + 1);
  x.e2 = std::pow(h->GetBinError(2 * step + q + 1), 2);
  return x;
}
inline V Add(V a, V b) { a.v += b.v; a.e2 += b.e2; return a; }

int DoReport()
{
  TFile *f = TFile::Open(kOut);
  if (!f || f->IsZombie()) { std::cerr << "[FATAL] run the loop first\n"; return 2; }
  gSystem->mkdir("plots", kTRUE);
  std::ofstream csv("edge_cutflow.csv");
  csv << "flavour,region,step,D_minus,D_plus,M_minus,M_plus,Wmc_minus,Wmc_plus,r_minus,r_minus_err,r_plus,r_plus_err,DR,DR_err,"
         "kDelta,kDelta_err\n";
  // region sets: bwd, fwd, both edges, central
  const char *rsName[4] = {"bwd (y0)", "fwd (y11)", "edges", "central"};
  const std::vector<std::vector<int>> rs = {{0}, {1}, {0, 1}, {2}};
  double DRtab[kNFl][4][kNStep], DRerr[kNFl][4][kNStep];
  for (int fl = 0; fl < kNFl; ++fl)
  {
    std::cout << "\n================ " << kFl[fl] << " ================\n";
    for (int s4 = 0; s4 < 4; ++s4)
    {
      std::cout << Form("--- %s ---   D = data, M = all MC (absolute, no SFs), W = W MC only; r = D/M; DR = r-/r+;"
                        " kDelta = [(D+ - O+) - (D- - O-)] / (W+ - W-) (QCD-free if QCD is charge-symmetric)\n", rsName[s4]);
      std::cout << "  step                       |   D-     D+   |    M-       M+   |  r-            r+            | DR = r-/r+      | kDelta\n";
      for (int st = 0; st < kNStep; ++st)
      {
        V D[2], M[2], W[2], O[2];
        for (int q = 0; q < 2; ++q)
          for (int r : rs[s4])
          {
            D[q] = Add(D[q], Get(f, fl, 0, r, st, q));
            W[q] = Add(W[q], Get(f, fl, 1, r, st, q));
            O[q] = Add(O[q], Get(f, fl, 2, r, st, q));
          }
        for (int q = 0; q < 2; ++q) M[q] = Add(W[q], O[q]);
        double rr[2], re[2];
        for (int q = 0; q < 2; ++q)
        {
          rr[q] = D[q].v / M[q].v;
          re[q] = rr[q] * std::sqrt(1.0 / std::max(D[q].v, 1.0) + M[q].e2 / (M[q].v * M[q].v));
        }
        const double dr = rr[0] / rr[1], dre = dr * std::hypot(re[0] / rr[0], re[1] / rr[1]);
        const double num = (D[1].v - O[1].v) - (D[0].v - O[0].v), den = W[1].v - W[0].v;
        const double kd = num / den, kde = std::sqrt(D[0].v + D[1].v) / std::fabs(den);
        DRtab[fl][s4][st] = dr;
        DRerr[fl][s4][st] = dre;
        std::cout << Form("  %-26s | %6.0f %6.0f | %8.1f %8.1f | %.3f+-%.3f  %.3f+-%.3f  | %.3f +- %.3f | %6.2f +- %.2f\n", kStep[st],
                          D[0].v, D[1].v, M[0].v, M[1].v, rr[0], re[0], rr[1], re[1], dr, dre, kd, kde);
        csv << Form("%s,%s,%d,%.0f,%.0f,%.3f,%.3f,%.3f,%.3f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.4f,%.4f\n", kFl[fl], rsName[s4], st, D[0].v,
                    D[1].v, M[0].v, M[1].v, W[0].v, W[1].v, rr[0], re[0], rr[1], re[1], dr, dre, kd, kde);
      }
    }
  }
  // the electron-over-muon double ratio of DR per step and region
  std::cout << "\n================ e/mu: DR_e / DR_mu = [(D/M)-/(D/M)+]_e / [(D/M)-/(D/M)+]_mu ================\n";
  std::cout << "  step                       |  bwd (y0)        fwd (y11)        edges            central\n";
  for (int st = 0; st < kNStep; ++st)
  {
    std::cout << Form("  %-26s |", kStep[st]);
    for (int s4 = 0; s4 < 4; ++s4)
    {
      const double v = DRtab[1][s4][st] / DRtab[0][s4][st];
      const double e = v * std::hypot(DRerr[1][s4][st] / DRtab[1][s4][st], DRerr[0][s4][st] / DRtab[0][s4][st]);
      std::cout << Form("  %.3f +- %.3f ", v, e);
    }
    std::cout << "\n";
  }
  // the step efficiencies per charge, data vs MC, where both are W-dominated (A3 -> A4: the m_T > 40 cut) and
  // for the ID (B0 -> B1) -- data contain QCD at B0, so the data "efficiency" mixes W and fakes
  std::cout << "\n================ step efficiencies per charge (data | MC), edges and central ================\n";
  for (int fl = 0; fl < kNFl; ++fl)
    for (int s4 : {2, 3})
      for (auto pr : std::vector<std::pair<int, int>>{{0, 1}, {1, 2}, {2, 3}, {5, 6}})
      {
        std::cout << Form("  %-2s %-8s %-26s -> %-26s:", kFl[fl], rsName[s4], kStep[pr.first], kStep[pr.second]);
        for (int q = 0; q < 2; ++q)
        {
          V d0, d1, m0, m1;
          for (int r : rs[s4])
          {
            d0 = Add(d0, Get(f, fl, 0, r, pr.first, q)); d1 = Add(d1, Get(f, fl, 0, r, pr.second, q));
            for (int g = 1; g < 3; ++g) { m0 = Add(m0, Get(f, fl, g, r, pr.first, q)); m1 = Add(m1, Get(f, fl, g, r, pr.second, q)); }
          }
          const double ed = d1.v / d0.v, ee = std::sqrt(ed * (1 - ed) / d0.v), em = m1.v / m0.v;
          std::cout << Form("  %s: data %.4f +- %.4f  MC %.4f", q ? "+" : "-", ed, ee, em);
        }
        std::cout << "\n";
      }

  // ---- phi: data/MC per charge at A4, both edges and central; the e-/e+ double ratio vs phi
  for (int fl = 0; fl < kNFl; ++fl)
    for (int s4 : {2, 3})
    {
      TGraphErrors *g[2] = {new TGraphErrors(), new TGraphErrors()}, *gdr = new TGraphErrors();
      double chi2 = 0;
      int ndf = 0;
      // normalize each charge's data/MC to its phi-integrated value: the plot shows the SHAPE in phi
      double tot[2];
      for (int q = 0; q < 2; ++q)
      {
        double d = 0, m = 0;
        for (int r : rs[s4])
        {
          d += ((TH1 *)f->Get(PName(fl, 0, q, r)))->Integral();
          for (int gg = 1; gg < 3; ++gg) m += ((TH1 *)f->Get(PName(fl, gg, q, r)))->Integral();
        }
        tot[q] = d / m;
      }
      for (int b = 1; b <= kNPhi; ++b)
      {
        double rr[2], re[2], x = 0, ex = 0;
        for (int q = 0; q < 2; ++q)
        {
          double d = 0, m = 0, m2 = 0;
          for (int r : rs[s4])
          {
            TH1 *hd = (TH1 *)f->Get(PName(fl, 0, q, r));
            d += hd->GetBinContent(b);
            x = hd->GetXaxis()->GetBinCenter(b);
            ex = 0.5 * hd->GetXaxis()->GetBinWidth(b);
            for (int gg = 1; gg < 3; ++gg)
            {
              TH1 *hm = (TH1 *)f->Get(PName(fl, gg, q, r));
              m += hm->GetBinContent(b);
              m2 += std::pow(hm->GetBinError(b), 2);
            }
          }
          rr[q] = (d / m) / tot[q];
          re[q] = rr[q] * std::sqrt(1.0 / std::max(d, 1.0) + m2 / (m * m));
          g[q]->SetPoint(b - 1, x, rr[q]);
          g[q]->SetPointError(b - 1, ex, re[q]);
        }
        const double dr = rr[0] / rr[1], dre = dr * std::hypot(re[0] / rr[0], re[1] / rr[1]);
        gdr->SetPoint(b - 1, x, dr);
        gdr->SetPointError(b - 1, ex, dre);
        chi2 += std::pow((dr - 1.0) / dre, 2);
        ++ndf;
      }
      std::cout << Form("[PHI] %s %s: (data/MC)-/(data/MC)+ vs phi, each normalized to its phi average: chi2 vs flat = %.1f / %d\n",
                        kFl[fl], rsName[s4], chi2, ndf - 1);
      auto ser = [](TGraphErrors *gg, const char *lab, int col, int mk, double sh)
      {
        OverlaySeries s;
        s.gTot = gg; s.label = lab; s.color = col; s.marker = mk; s.markerSize = 1.4; s.xShift = sh;
        return s;
      };
      auto tuner = [](double lo, double hi) {
        return GraphTuner([lo, hi](TCanvas *c, TGraphErrors *gg)
        {
          gg->SetMinimum(lo); gg->SetMaximum(hi);
          c->cd();
          TLine *l = new TLine(gg->GetXaxis()->GetXmin(), 1.0, gg->GetXaxis()->GetXmax(), 1.0);
          l->SetLineStyle(2); l->SetLineColor(kGray + 2); l->Draw();
        });
      };
      PlotStyle ps;
      ps.headerX = 0.55; ps.headerY = 0.90; ps.subSize = 0.032; ps.boxTextSize = 0.028;
      ps.boxX1 = 0.45; ps.boxX2 = 0.93; ps.boxY1 = 0.66; ps.boxY2 = 0.74;
      ps.legX1 = 0.50; ps.legY1 = 0.15; ps.legX2 = 0.93; ps.legY2 = 0.26;
      const std::string lepM = fl ? "e^{-}" : "#mu^{-}", lepP = fl ? "e^{+}" : "#mu^{+}";
      SaveNiceGraph_Overlay({ser(g[0], (lepM + " (data/MC)").c_str(), kRed + 1, 33, -0.25),
                             ser(g[1], (lepP + " (data/MC)").c_str(), kAzure + 2, 21, 0.25)},
                            Form("plots/phi_%s_%s", kFl[fl], s4 == 2 ? "edges" : "central"), "#phi of the lepton",
                            "data/MC, normalized to its #phi average", "", Form("W, m_{T} > 40, %s", rsName[s4]),
                            fl ? "electrons" : "muons", {Form("#chi^{2}(ratio of the two flat) = %.1f / %d", chi2, ndf - 1)}, ps,
                            tuner(0.0, 2.5));
    }

  // ---- runs: the negative-charge fraction per run at A4, edges, data
  for (int fl = 0; fl < kNFl; ++fl)
    for (int s4 : {2, 3})
    {
      std::vector<double> nm, np;
      TH1 *ref = (TH1 *)f->Get(Form("h_run_%s_minus_bwd", kFl[fl]));
      const int nb = ref->GetNbinsX();
      double Nm = 0, Np = 0;
      for (int b = 1; b <= nb; ++b)
      {
        double a = 0, c = 0;
        for (int r : rs[s4])
        {
          a += ((TH1 *)f->Get(Form("h_run_%s_minus_%s", kFl[fl], kReg[r])))->GetBinContent(b);
          c += ((TH1 *)f->Get(Form("h_run_%s_plus_%s", kFl[fl], kReg[r])))->GetBinContent(b);
        }
        nm.push_back(a); np.push_back(c); Nm += a; Np += c;
      }
      const double p = Nm / (Nm + Np);
      double chi2 = 0;
      int ndf = 0;
      std::cout << Form("[RUN] %s %s: overall %s-fraction %.3f (%0.f / %.0f); per run:", kFl[fl], rsName[s4], fl ? "e" : "mu", p, Nm,
                        Nm + Np);
      for (int b = 0; b < nb; ++b)
      {
        const double n = nm[b] + np[b];
        if (n < 1) continue;
        const double fr = nm[b] / n, e = std::sqrt(p * (1 - p) / n);
        chi2 += std::pow((fr - p) / e, 2);
        ++ndf;
        std::cout << Form(" %s:%.0f/%.0f", ref->GetXaxis()->GetBinLabel(b + 1), nm[b], n);
      }
      std::cout << Form("\n[RUN] %s %s: chi2 of the per-run fraction vs constant = %.1f / %d\n", kFl[fl], rsName[s4], chi2, ndf - 1);
    }
  csv.close();
  std::cout << "[OUT] edge_cutflow.csv, plots/phi_{mu,e}_{edges,central}\n";
  return 0;
}
} // namespace EC

int edge_cutflow(const char *mode = "loop")
{
  const std::string m = mode;
  if (m == "loop") return EC::DoLoop();
  if (m == "report") return EC::DoReport();
  std::cerr << "usage: edge_cutflow.C+(\"loop\"|\"report\")\n";
  return 1;
}
