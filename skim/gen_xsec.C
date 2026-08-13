// =============================================================================
// gen_xsec.C -- GENERATOR-LEVEL W cross sections per rapidity bin (2026-08-05;
// FIDUCIAL since 2026-08-12).
//
// Loops the four W MC files (Wp/Wm x mu/ele) over ALL generated events -- no
// reco, no selection -- and histograms the gen charged lepton's LAB eta in the
// analysis binning (pOSkim::kYEdges, + the FB edge set), in the SAME
// y = -eta_lab convention as the reco skim (p-going = forward; see the
// [NOTICE] at the fill). The pO-scaled
// PER-FLAVOUR cross section in bin i is the weighted FRACTION times the
// single-source cross section (unit-proof: the July-29 weights are sigma in
// pb, <w> ~ 6376, while mc_norm's kSigma_* are nb -- the fraction cancels the
// weight unit exactly like k_s = A*sigma*L/Sum(w) does in the reco chain):
//
//     sigma_i = kA_O * kSigma_{Wp,Wm} * (Sum_i w) / (Sum_all w)      [nb]
//
// (the gen-level analog of the reco normalization: reco-template integral
// / L / sigma_i = acceptance x efficiency of bin i).
// The mu and e samples are two statistically independent estimates of the SAME
// distribution (lepton universality), so both files are accumulated together.
//
// Gen lepton choice: highest-pT gen lepton of the sample's flavour and charge
// WITH a W (|pdg|=24) ancestor; if the ntuple's mother chain never reaches a W
// (counted + warned), fall back to the highest-pT flavour+charge match --
// for exclusive W->l nu samples that is the decay lepton either way.
//
// FIDUCIAL definition (2026-08-12): the primary histograms apply the gen-level
// twin of the reco selection cuts -- lepton pT > kFidPtMin (= the skim's
// nominal leading-lepton cut, 25 GeV) with |eta_lab| < 2.4 implicit in the
// binning window. That makes sigma_i a FIDUCIAL cross section, so the fitted
// signal strength converts directly: sigma_meas,i = r_i * sigma_i^gen,fid
// (= N_fit,i / L / (A*eps)_MC,i -- algebraically identical to the yield-based
// extraction with MC eff/acc, but kA_O and kSigma_* cancel between r and
// sigma_gen). The lepton is the BARE post-FSR gen lepton (ntuple gen
// collection); no m_T cut in the fiducial even for the leppt_mt40 variant --
// ONE fiducial definition, the m_T-cut efficiency stays inside eps.
//
// Output: rootfile/gen_xsec.root
//   h_gen_sig_{Wp,Wm}      (lab edges;  bin content = per-flavour FIDUCIAL
//                           sigma_i in nb, gen lepton pT > 25)
//   h_gen_sig_{Wp,Wm}_FB   (FB edges;   same content convention)
//   h_gen_tot_{Wp,Wm}[_FB] (reference:  NO pT cut -- the pre-2026-08-12
//                           definition; NB "tot" is NOT a true no-cut sigma:
//                           the ntuple gen collection is itself FILTERED at
//                           pT > 5 GeV and |eta| < 2.5 (measured 2026-08-12),
//                           so fid/tot is a pT>25 / pT>5 ratio, not acceptance)
//
// The gen filter is why 25-37% of events report "no gen lepton" below: those
// leptons are below 5 GeV or beyond |eta| = 2.5 -- outside the fiducial either
// way, so sigma_fid is unaffected (every pT > 25, |eta| < 2.4 lepton IS stored)
// and the Sum(w) denominator still runs over ALL events.
// Sumw2 carries the MC-stat error. Consumed by
// plotting/xsec_fiducial.C::xsec_fiducial_comb, which draws 2 x sigma_i
// (= mu+e summed) next to the reco-level expectation and the post-fit points.
// Discriminant-independent: produce ONCE, serves met/leppt/leppt_mt40 alike.
//
//   cd skim/ && root -l -b -q 'gen_xsec.C+'
// =============================================================================
#include "skim_common.h"
#include "mc_norm.h"

#include "TFile.h"
#include "TInterpreter.h"
#include "TTree.h"
#include "TH1D.h"
#include "TString.h"
#include "TSystem.h"
#include <algorithm>
#include <iostream>
#include <vector>

namespace {

// Gen-level fiducial lepton-pT cut = the skim's nominal leading-lepton cut.
constexpr double kFidPtMin = 25.0;

// Accumulate one W MC file into the (shared) weighted gen-eta histograms:
// hLab/hFB get the FIDUCIAL fill (pT > kFidPtMin), hLabTot/hFBTot every lepton.
bool AccumulateGen(const char *fname, int flavPdg, int wantChg,
                   TH1D *hLab, TH1D *hFB, TH1D *hLabTot, TH1D *hFBTot,
                   double &sumw, long long &nraw,
                   long long &nNoLep, long long &nNoWAnc)
{
  TFile *f = TFile::Open(fname, "READ");
  if (!f || f->IsZombie()) { std::cerr << "[ERROR] cannot open " << fname << "\n"; return false; }
  TTree *tGen = (TTree *)f->Get("HiGenParticleAna/hi");
  TTree *tHi  = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  if (!tGen || !tHi)
  { std::cerr << "[ERROR] missing HiGenParticleAna/hi or hiEvtAnalyzer/HiTree in " << fname << "\n"; f->Close(); return false; }
  const Long64_t n = tGen->GetEntries();
  if (tHi->GetEntries() != n)
  { std::cerr << "[ERROR] gen/HiTree entry mismatch in " << fname << "\n"; f->Close(); return false; }

  // fast-skim discipline: disable everything, enable only what is read
  std::vector<float> *pt = nullptr, *eta = nullptr;
  std::vector<int> *chg = nullptr, *pdg = nullptr;
  std::vector<std::vector<int>> *motherIdx = nullptr;
  tGen->SetBranchStatus("*", 0);
  tGen->SetBranchStatus("pt", 1);        tGen->SetBranchAddress("pt", &pt);
  tGen->SetBranchStatus("eta", 1);       tGen->SetBranchAddress("eta", &eta);
  tGen->SetBranchStatus("chg", 1);       tGen->SetBranchAddress("chg", &chg);
  tGen->SetBranchStatus("pdg", 1);       tGen->SetBranchAddress("pdg", &pdg);
  tGen->SetBranchStatus("motherIdx", 1); tGen->SetBranchAddress("motherIdx", &motherIdx);
  Float_t weight = 1.f;
  tHi->SetBranchStatus("*", 0);
  tHi->SetBranchStatus("weight", 1);     tHi->SetBranchAddress("weight", &weight);

  for (Long64_t ie = 0; ie < n; ++ie)
  {
    tGen->GetEntry(ie);
    tHi->GetEntry(ie);
    const double w = (double)weight;
    sumw += w;
    ++nraw;
    if (!pt || !eta || !chg || !pdg || !motherIdx) continue;
    const size_t ng = std::min({pt->size(), eta->size(), chg->size(),
                                pdg->size(), motherIdx->size()});
    int best = -1, bestAnc = -1;
    for (size_t i = 0; i < ng; ++i)
    {
      if (std::abs(pdg->at(i)) != flavPdg) continue;
      if (chg->at(i) != wantChg) continue;
      if (best < 0 || pt->at(i) > pt->at(best)) best = (int)i;
      if (pOSkim::HasAncestor((int)i, 24, pdg, motherIdx))
        if (bestAnc < 0 || pt->at(i) > pt->at(bestAnc)) bestAnc = (int)i;
    }
    const int use = (bestAnc >= 0) ? bestAnc : best;
    if (bestAnc < 0 && best >= 0) ++nNoWAnc; // fallback (mother chain w/o W)
    if (use < 0) { ++nNoLep; continue; }     // no gen lepton of this flavour+charge
    // [NOTICE] Sign flip: p-going (-Z) defined as forward -- MUST match the
    // reco convention in skim.C (`const double y = -{mu,ele}Eta->at(iLead)`),
    // or gen bin i pairs with the MIRRORED reco region y_i. Caught 2026-08-12
    // by the per-bin (A x eps) diagnostic: mirrored pairing made A x eps
    // charge-dependent (mu W- ran 1.27 -> 0.73 across eta) instead of the
    // charge-symmetric detector response it must be.
    const double y = -eta->at(use);
    hLabTot->Fill(y, w);
    hFBTot->Fill(y, w);
    if (pt->at(use) > kFidPtMin)
    {
      hLab->Fill(y, w);
      hFB->Fill(y, w);
    }
  }
  f->Close();
  delete f;
  return true;
}

} // namespace

void gen_xsec()
{
  using namespace pOSkim;
  TH1::SetDefaultSumw2(kTRUE);
  gSystem->mkdir("rootfile", kTRUE);
  // motherIdx is a vector<vector<int>> branch -- same dictionary the skim needs
  gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");

  const char *cname[2]     = {"Wp", "Wm"};
  const SampleType samp[2] = {kWp, kWm};
  const int wantChg[2]     = {+1, -1};
  const double sigNN[2]    = {pONorm::kSigma_Wp, pONorm::kSigma_Wm}; // nb, single source

  TFile *fout = TFile::Open("rootfile/gen_xsec.root", "RECREATE");
  printf("[gen_xsec] kA_O = %.1f, sigma_NN(Wp/Wm) = %.3f/%.3f nb (mc_norm.h);"
         " sigma_i = A * sigma * Sumw_i/Sumw_all\n",
         pONorm::kA_O, sigNN[0], sigNN[1]);

  for (int ic = 0; ic < 2; ++ic)
  {
    fout->cd();
    TH1D *hLab = new TH1D(Form("h_gen_sig_%s", cname[ic]),
                          Form("gen fiducial d#sigma bins (p_{T}>%.0f), %s (per flavour);#eta^{l}_{lab};#sigma_{i} (nb)", kFidPtMin, cname[ic]),
                          kNY, kYEdges);
    TH1D *hFB = new TH1D(Form("h_gen_sig_%s_FB", cname[ic]),
                         Form("gen fiducial d#sigma bins (FB edges, p_{T}>%.0f), %s (per flavour);#eta^{l}_{lab};#sigma_{i} (nb)", kFidPtMin, cname[ic]),
                         kNY, kYEdgesFB);
    TH1D *hLabTot = new TH1D(Form("h_gen_tot_%s", cname[ic]),
                             Form("gen d#sigma bins (no p_{T} cut), %s (per flavour);#eta^{l}_{lab};#sigma_{i} (nb)", cname[ic]),
                             kNY, kYEdges);
    TH1D *hFBTot = new TH1D(Form("h_gen_tot_%s_FB", cname[ic]),
                            Form("gen d#sigma bins (FB edges, no p_{T} cut), %s (per flavour);#eta^{l}_{lab};#sigma_{i} (nb)", cname[ic]),
                            kNY, kYEdgesFB);

    double sumw = 0; long long nraw = 0, nNoLep = 0, nNoWAnc = 0;
    const char *flavs[2] = {"mu", "ele"};
    for (int fl = 0; fl < 2; ++fl)
    {
      const int flavPdg = (fl == 0) ? 13 : 11;
      SampleFileInfo info = ResolveMCSample(samp[ic], flavs[fl]);
      printf("[gen_xsec] %s %-3s <- %s\n", cname[ic], flavs[fl], info.fname.c_str());
      AccumulateGen(info.fname.c_str(), flavPdg, wantChg[ic], hLab, hFB,
                    hLabTot, hFBTot, sumw, nraw, nNoLep, nNoWAnc);
    }
    if (nraw <= 0 || sumw <= 0) { std::cerr << "[ERROR] no events accumulated for " << cname[ic] << "\n"; continue; }

    const double sigAll = pONorm::kA_O * sigNN[ic]; // all eta, per flavour, pO-scaled
    const double scale  = sigAll / sumw;            // weighted-fraction normalization
    hLab->Scale(scale);
    hFB->Scale(scale);
    hLabTot->Scale(scale);
    hFBTot->Scale(scale);

    printf("[gen_xsec] %s: Nraw = %lld, sigma(all eta) = %.2f nb,"
           " sigma(|eta_lab|<2.4) = %.2f nb, sigma_fid(pT>%.0f, |eta_lab|<2.4) = %.2f nb"
           " (per flavour, pO-scaled; pT-acceptance = %.3f)\n",
           cname[ic], nraw, sigAll, hLabTot->Integral(), kFidPtMin, hLab->Integral(),
           hLabTot->Integral() > 0 ? hLab->Integral() / hLabTot->Integral() : 0.0);
    if (nNoLep > 0)
      printf("[gen_xsec][WARN] %s: %lld events without a gen %s of charge %+d (skipped)\n",
             cname[ic], nNoLep, "lepton", wantChg[ic]);
    if (nNoWAnc > 0)
      printf("[gen_xsec][WARN] %s: %lld events used the highest-pT fallback (no W ancestor in the mother chain)\n",
             cname[ic], nNoWAnc);

    fout->cd();
    hLab->Write("", TObject::kOverwrite);
    hFB->Write("", TObject::kOverwrite);
    hLabTot->Write("", TObject::kOverwrite);
    hFBTot->Write("", TObject::kOverwrite);
  }
  fout->Close();
  delete fout;
  std::cout << "[gen_xsec] wrote rootfile/gen_xsec.root\n";
}
