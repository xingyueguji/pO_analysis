// SCRATCH STUDY (2026-10-20): e+ vs e- MET and m_T distributions at the W selection
// (skim step 8, NO MET / m_T cut), data vs MC, per region of the electron:
//   edges |eta| > 2.0 (lab y0 + y11) | mid 1.566-2.0 | cen |eta| < 1.4442
// for the PF MET and for the electron-corrected MET (the PF candidates within
// DR < 0.1 of the electron replaced by elePt, met_corr.C), and the m_T of each.
// Per charge: data vs a stack of W MC + other MC (DY, tau, ttbar) scaled by
// kRW = 1.2 (the muon prefit r) x 0.887 (the electron <SF>; no SFs in this loop)
// + a multijet template = that charge's anti-iso data (not iso-pass, relIso in
// [0.3, 1.0)) minus the anti-iso MC, normalized to iso-pass data - MC in the
// low bins (MET < 20 / m_T < 30) of the same variable.
// QCD-free: the charge difference data(e+) - data(e-) vs kRW x [MC(e+) - MC(e-)].
// Run: root -l -b -q dist_charge.C+   -> dist_charge.root, plots/dist/*, dist_charge.log (console)
#include "/Users/zhenghuang/pO_analysis/skim/skim_common.h"
#include "/Users/zhenghuang/pO_analysis/skim/mc_norm.h"
#include "/Users/zhenghuang/pO_analysis/plotting/plotting_helper.C"

#include "TFile.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector2.h"

#include <iostream>
#include <string>
#include <vector>

using namespace pOSkim;

namespace DC
{
const int kNReg = 3, kNQ = 2, kNVar = 4, kNGrp = 5;
const char *const kReg[kNReg]  = {"edges", "mid", "cen"};
const char *const kRegLab[kNReg] = {"|#eta| > 2.0", "1.566 < |#eta| < 2.0", "|#eta| < 1.4442"};
const char *const kQ[kNQ]      = {"em", "ep"};
const char *const kQLab[kNQ]   = {"e^{-}", "e^{+}"};
const char *const kVar[kNVar]  = {"met", "metcorr", "mt", "mtcorr"};
const char *const kVarLab[kNVar] = {"PF MET [GeV]", "corrected MET [GeV]", "m_{T} (PF MET) [GeV]", "m_{T} (corrected MET) [GeV]"};
const double kLowCut[kNVar] = {20.0, 20.0, 30.0, 30.0}; // the multijet normalization window
const char *const kGrp[kNGrp]  = {"data", "data_anti", "mcW", "mcO", "mc_anti"};
const double kRW = 1.2 * 0.887;
const char *const kNgen = "/Users/zhenghuang/pO_analysis/skim/rootfile/ngen.root";
TH1D *H[kNGrp][kNReg][kNQ][kNVar];
void Plot();

inline int RegionOf(double eta)
{
  const double a = std::fabs(eta);
  if (a > 2.0) return 0;
  if (a > kEcalGapHi) return 1;
  return 2;
}

void Book()
{
  for (int g = 0; g < kNGrp; ++g)
    for (int r = 0; r < kNReg; ++r)
      for (int q = 0; q < kNQ; ++q)
        for (int v = 0; v < kNVar; ++v)
        {
          const bool isMt = v >= 2;
          H[g][r][q][v] = new TH1D(Form("h_%s_%s_%s_%s", kGrp[g], kReg[r], kQ[q], kVar[v]), Form(";%s;events", kVarLab[v]),
                                   isMt ? 14 : 20, 0, isMt ? 140 : 100);
          H[g][r][q][v]->Sumw2();
        }
}

// g0 = the iso-pass group of the sample (0 data, 2 mcW, 3 mcO); its anti-iso group is 1 (data) or 4 (MC)
int Loop(const std::string &fname, int gIso, double ks)
{
  const bool isMC = gIso != 0;
  const int gAnti = isMC ? 4 : 1;
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
  long long nIso = 0, nAnti = 0;
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
    auto relIso = [&](int i) { return RelIsoPF(i, pt, ch, neu, pho, pu); };
    auto isoOk = [&](int i) { return iso->at(i) != 0 && relIso(i) < kEleIsoRelMax; };
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
    if (iL < 0 || pt->at(iL) <= 25.0 || std::fabs(eta->at(iL)) > 2.4) continue;
    const bool pass = isoOk(iL);
    const double ri = relIso(iL);
    const bool anti = !pass && ri >= 0.3 && ri < 1.0;
    if (!pass && !anti) continue;
    obj->GetEntry(ie);
    if (!PassLeadingLeptonTrigMatch(0.4, iL, eta, phi, true, toPt, true, toEta, true, toPhi, warnedT, "electron")) continue;
    pf->GetEntry(ie);
    const TVector2 met = ComputePFMET(nullptr, pfPt, pfPhi);
    double sum = 0;
    for (size_t k = 0; k < pfPt->size(); ++k)
      if (DeltaR(pfEta->at(k), pfPhi->at(k), eta->at(iL), phi->at(iL)) < 0.1) sum += pfPt->at(k);
    const TVector2 metC = met + (sum - pt->at(iL)) * TVector2(std::cos(phi->at(iL)), std::sin(phi->at(iL)));
    const double x[kNVar] = {met.Mod(), metC.Mod(), TransverseMass(pt->at(iL), phi->at(iL), met),
                             TransverseMass(pt->at(iL), phi->at(iL), metC)};
    const int r = RegionOf(eta->at(iL)), qi = q->at(iL) > 0 ? 1 : 0, g = pass ? gIso : gAnti;
    const double w = (isMC ? (double)w0 : 1.0) * ks;
    for (int v = 0; v < kNVar; ++v)
    {
      TH1D *hh = H[g][r][qi][v];
      hh->Fill(std::min(x[v], hh->GetXaxis()->GetXmax() - 1e-6), w); // overflow into the last bin
    }
    (pass ? nIso : nAnti)++;
  }
  std::cout << Form("[LOOP] %s: iso-pass %lld, anti-iso %lld (raw)\n", fname.c_str(), nIso, nAnti);
  f->Close();
  return 0;
}
} // namespace DC

void dist_charge()
{
  using namespace DC;
  Book();
  const std::string b = "/Users/zhenghuang/pO_2026_Aug_20/Aug_20_MC_";
  Loop(kDefaultDataFile, 0, 1.0);
  Loop(b + "Wp_ele.root", 2, pONorm::MCScale("Wp_ele", kNgen));
  Loop(b + "Wm_ele.root", 2, pONorm::MCScale("Wm_ele", kNgen));
  for (auto &s : std::vector<std::pair<std::string, const char *>>{{"DY_ele_Z", "DYee"}, {"DY_tau_Z", "DYtau"}, {"Wp_tau", "Wp_tau"},
                                                                    {"Wm_tau", "Wm_tau"}, {"TT", "TT"}})
    Loop(b + s.first + ".root", 3, pONorm::MCScale(s.second, kNgen));
  TFile *fo = TFile::Open("dist_charge.root", "RECREATE");
  for (auto &a : H) for (auto &bb : a) for (auto &c : bb) for (TH1D *x : c) x->Write();
  fo->Close();
  DC::Plot();
}

// redraw from dist_charge.root (no loops)
void dist_charge_replot()
{
  using namespace DC;
  TFile *fi = TFile::Open("dist_charge.root");
  for (int g = 0; g < kNGrp; ++g)
    for (int r = 0; r < kNReg; ++r)
      for (int q = 0; q < kNQ; ++q)
        for (int v = 0; v < kNVar; ++v)
        {
          H[g][r][q][v] = (TH1D *)fi->Get(Form("h_%s_%s_%s_%s", kGrp[g], kReg[r], kQ[q], kVar[v]));
          H[g][r][q][v]->SetDirectory(nullptr);
        }
  DC::Plot();
}

void DC::Plot()
{
  using namespace DC;
  gSystem->mkdir("plots/dist", kTRUE);
  for (int r = 0; r < kNReg; ++r)
    for (int v = 0; v < kNVar; ++v)
    {
      TH1D *mcTot[kNQ], *dat[kNQ];
      for (int qi = 0; qi < kNQ; ++qi)
      {
        TH1D *d = (TH1D *)H[0][r][qi][v]->Clone(Form("d_%d_%d_%d", r, qi, v));
        TH1D *w = (TH1D *)H[2][r][qi][v]->Clone(Form("w_%d_%d_%d", r, qi, v));
        TH1D *o = (TH1D *)H[3][r][qi][v]->Clone(Form("o_%d_%d_%d", r, qi, v));
        w->Scale(kRW);
        o->Scale(kRW);
        // the multijet template: anti-iso data - anti-iso MC, normalized in the low window
        TH1D *qcd = (TH1D *)H[1][r][qi][v]->Clone(Form("q_%d_%d_%d", r, qi, v));
        TH1D *ma = (TH1D *)H[4][r][qi][v]->Clone(Form("ma_%d_%d_%d", r, qi, v));
        ma->Scale(kRW);
        qcd->Add(ma, -1.0);
        for (int i = 1; i <= qcd->GetNbinsX(); ++i)
          if (qcd->GetBinContent(i) < 0) qcd->SetBinContent(i, 0.0);
        const int bLo = 1, bHi = qcd->GetXaxis()->FindBin(kLowCut[v] - 1e-6);
        const double need = d->Integral(bLo, bHi) - w->Integral(bLo, bHi) - o->Integral(bLo, bHi);
        const double have = qcd->Integral(bLo, bHi);
        qcd->Scale(have > 0 ? std::max(0.0, need) / have : 0.0);
        const double nD = d->Integral(), nW = w->Integral(), nO = o->Integral(), nQ = qcd->Integral();
        const int bHiCut = d->GetXaxis()->FindBin((v >= 2 ? 40.0 : 30.0) + 1e-6);
        const double dHi = d->Integral(bHiCut, d->GetNbinsX()), mHi = w->Integral(bHiCut, d->GetNbinsX()) +
                                                                           o->Integral(bHiCut, d->GetNbinsX()) +
                                                                           qcd->Integral(bHiCut, d->GetNbinsX());
        std::cout << Form("[DIST] %-5s %-7s %-3s: data %5.0f | W %7.1f other %6.1f QCD %7.1f -> total %7.1f | above %s: data %4.0f vs model %6.1f (%.3f)\n",
                          kReg[r], kVar[v], kQ[qi], nD, nW, nO, nQ, nW + nO + nQ, v >= 2 ? "m_T 40" : "MET 30", dHi, mHi, dHi / mHi);
        PlotStyle ps;
        ps.normBkgToData = false;
        ps.pullPad = true;
        ps.headerX = 0.55;
        ps.legX1 = 0.58; ps.legY1 = 0.60; ps.legX2 = 0.93; ps.legY2 = 0.78;
        ps.boxX1 = 0.50; ps.boxX2 = 0.93; ps.boxY1 = 0.48; ps.boxY2 = 0.58; ps.boxTextSize = 0.026;
        SaveNicePlot1D_WithBkg(d, {w, o, qcd}, {"W#rightarrowe#nu MC", "DY, #tau, t#bar{t} MC", "multijet (anti-iso)"},
                               Form("plots/dist/%s_%s_%s", kVar[v], kReg[r], kQ[qi]), kVarLab[v], "events", "",
                               Form("%s, %s", kQLab[qi], kRegLab[r]), "W selection, no MET / m_{T} cut",
                               {Form("MC #times %.3f; multijet norm. below %.0f", kRW, kLowCut[v]),
                                Form("above %s: data %.0f, model %.1f", v >= 2 ? "m_{T} 40" : "MET 30", dHi, mHi)},
                               ps);
        dat[qi] = d;
        mcTot[qi] = (TH1D *)w->Clone(Form("mt_%d_%d_%d", r, qi, v));
        mcTot[qi]->Add(o);
      }
      // the QCD-free charge difference
      TH1D *dd = (TH1D *)dat[1]->Clone(Form("dd_%d_%d", r, v));
      dd->Add(dat[0], -1.0);
      TH1D *md = (TH1D *)mcTot[1]->Clone(Form("md_%d_%d", r, v));
      md->Add(mcTot[0], -1.0);
      double chi2 = 0;
      int nb = 0;
      for (int i = 1; i <= dd->GetNbinsX(); ++i)
      {
        const double e2 = std::pow(dd->GetBinError(i), 2) + std::pow(md->GetBinError(i), 2);
        if (e2 <= 0) continue;
        chi2 += std::pow(dd->GetBinContent(i) - md->GetBinContent(i), 2) / e2;
        ++nb;
      }
      const int bHiCut = dd->GetXaxis()->FindBin((v >= 2 ? 40.0 : 30.0) + 1e-6);
      double eHi = 0, eLo = 0;
      const double dHi = dd->IntegralAndError(bHiCut, dd->GetNbinsX(), eHi), dLo = dd->IntegralAndError(1, bHiCut - 1, eLo);
      const double mHi = md->Integral(bHiCut, md->GetNbinsX()), mLo = md->Integral(1, bHiCut - 1);
      std::cout << Form("[DIFF] %-5s %-7s: data(e+)-data(e-) below %s %+6.0f +- %4.0f vs MC %+6.1f | above %+6.0f +- %4.0f vs MC %+6.1f | chi2 %.1f/%d\n",
                        kReg[r], kVar[v], v >= 2 ? "m_T 40" : "MET 30", dLo, eLo, mLo, dHi, eHi, mHi, chi2, nb);
      PlotStyle ps;
      ps.normBkgToData = false;
      ps.pullPad = false; // a pull against sqrt(MC) is meaningless for a difference
      ps.headerX = 0.55;
      ps.legX1 = 0.55; ps.legY1 = 0.66; ps.legX2 = 0.93; ps.legY2 = 0.76;
      ps.boxX1 = 0.50; ps.boxX2 = 0.93; ps.boxY1 = 0.50; ps.boxY2 = 0.64; ps.boxTextSize = 0.026;
      SaveNicePlot1D_WithBkg(dd, {md}, {Form("MC (e^{+} #minus e^{-}) #times %.3f", kRW)}, Form("plots/dist/diff_%s_%s", kVar[v], kReg[r]),
                             kVarLab[v], "events (e^{+} #minus e^{-})", "", Form("e^{+} #minus e^{-}, %s", kRegLab[r]),
                             "W selection, no MET / m_{T} cut",
                             {Form("below: data %+.0f #pm %.0f, MC %+.1f", dLo, eLo, mLo), Form("above: data %+.0f #pm %.0f, MC %+.1f", dHi, eHi, mHi),
                              Form("#chi^{2} = %.1f / %d", chi2, nb)},
                             ps);
    }
}
