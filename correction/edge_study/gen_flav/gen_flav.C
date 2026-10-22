// scratch (2026-10-20): are the W->e nu and W->mu nu MC samples the same physics at
// generator level? Per sample, the W decay lepton from the EventTree mc* block
// (status 1, |pdg| = flavour, mother W -- or grandmother for an FSR copy), bare,
// fiducial pT > 25, |eta| < 2.4; sigma per lab bin y = -eta (kYEdges) =
// kA_O * kSigma * Sum w(bin) / Sum w(all events). Then the reco signal template of
// the leppt_mt40 Combine input / (L * sigma_gen) = A.eps per flavour, charge, bin.
// Run from skim/: root -l -b -q <this file>+
#include "/Users/zhenghuang/pO_analysis/skim/skim_common.h"
#include "/Users/zhenghuang/pO_analysis/skim/mc_norm.h"

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include <cmath>
#include <iostream>
#include <vector>

struct GenOut
{
  double sig[12] = {}, tot = 0, nolep = 0;
};

GenOut LoopGen(const char *fname, int flav, int q, double sigmaNN)
{
  GenOut o;
  TFile *f = TFile::Open(fname);
  TTree *t = (TTree *)f->Get("ggHiNtuplizer/EventTree"), *h = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  std::vector<int> *pid = nullptr, *st = nullptr, *mom = nullptr, *gmom = nullptr;
  std::vector<float> *pt = nullptr, *eta = nullptr;
  Float_t w = 0;
  t->SetBranchStatus("*", 0);
  for (const char *b : {"mcPID", "mcStatus", "mcMomPID", "mcGMomPID", "mcPt", "mcEta"}) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("mcPID", &pid);  t->SetBranchAddress("mcStatus", &st);
  t->SetBranchAddress("mcMomPID", &mom); t->SetBranchAddress("mcGMomPID", &gmom);
  t->SetBranchAddress("mcPt", &pt);    t->SetBranchAddress("mcEta", &eta);
  h->SetBranchStatus("*", 0); h->SetBranchStatus("weight", 1); h->SetBranchAddress("weight", &w);
  const int want = (q > 0) ? -flav : flav; // l+ has pdg -11/-13
  double sumAll = 0, sumBin[12] = {};
  const Long64_t n = t->GetEntries();
  for (Long64_t i = 0; i < n; ++i)
  {
    t->GetEntry(i);
    h->GetEntry(i);
    sumAll += w;
    int best = -1;
    for (size_t k = 0; k < pid->size(); ++k)
    {
      if (pid->at(k) != want || st->at(k) != 1) continue;
      if (std::abs(mom->at(k)) != 24 && !(std::abs(mom->at(k)) == flav && std::abs(gmom->at(k)) == 24)) continue;
      if (best < 0 || pt->at(k) > pt->at(best)) best = (int)k;
    }
    if (best < 0) { o.nolep += 1; continue; }
    if (pt->at(best) <= 25.0 || std::fabs(eta->at(best)) >= 2.4) continue;
    const double y = -eta->at(best);
    for (int b = 0; b < 12; ++b)
      if (y >= pOSkim::kYEdges[b] && y < pOSkim::kYEdges[b + 1]) { sumBin[b] += w; break; }
  }
  for (int b = 0; b < 12; ++b) o.sig[b] = pONorm::kA_O * sigmaNN * sumBin[b] / sumAll;
  for (int b = 0; b < 12; ++b) o.tot += o.sig[b];
  f->Close();
  return o;
}

void gen_flav()
{
  const char *base = "/Users/zhenghuang/pO_2026_Aug_20/Aug_20_MC_";
  const char *cn[2] = {"Wp", "Wm"};
  const double sig[2] = {pONorm::kSigma_Wp, pONorm::kSigma_Wm};
  GenOut g[2][2]; // [charge][flav: 0 mu, 1 e]
  for (int c = 0; c < 2; ++c)
    for (int fl = 0; fl < 2; ++fl)
    {
      g[c][fl] = LoopGen(Form("%s%s_%s.root", base, cn[c], fl ? "ele" : "mu"), fl ? 11 : 13, c == 0 ? +1 : -1, sig[c]);
      std::cout << Form("[GEN] %s %-3s fiducial sigma %.3f nb, events without the W lepton %.0f\n", cn[c], fl ? "e" : "mu",
                        g[c][fl].tot, g[c][fl].nolep);
    }
  // reco signal templates (leppt_mt40 input, absolute, both W samples in the reco-charge region)
  TFile *fin[2] = {TFile::Open("/Users/zhenghuang/pO_analysis/plotting/plots/combine_input_W_leppt_mt40.root"),
                   TFile::Open("/Users/zhenghuang/pO_analysis/plotting/plots/Elec/combine_input_W_leppt_mt40.root")};
  const double L = pONorm::kLumi_invnb;
  for (int c = 0; c < 2; ++c)
  {
    std::cout << Form("\n[GEN] %s per lab bin: sigma_gen (nb) mu | e | e/mu (gen MC stat ~0.3%%)   ||   A.eps = S_reco / (L sigma_gen): mu | e | e/mu\n", cn[c]);
    for (int b = 0; b < 12; ++b)
    {
      const double gm = g[c][0].sig[b], ge = g[c][1].sig[b];
      double S[2];
      for (int fl = 0; fl < 2; ++fl)
      {
        TH1 *hs = (TH1 *)fin[fl]->Get(Form("%s_lab_y%d/signal", cn[c], b));
        S[fl] = hs ? hs->Integral(0, hs->GetNbinsX() + 1) : -1;
      }
      const double am = S[0] / (L * gm), ae = S[1] / (L * ge);
      std::cout << Form("  y%-2d [%5.2f,%5.2f) | %.4f %.4f %.4f || %.4f %.4f %.4f\n", b, pOSkim::kYEdges[b], pOSkim::kYEdges[b + 1], gm,
                        ge, ge / gm, am, ae, ae / am);
    }
  }
}
