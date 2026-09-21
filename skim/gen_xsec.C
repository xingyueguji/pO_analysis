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
//   h_gen_sig_Z            (2026-09-15; ONE bin = the Z->ll FIDUCIAL sigma in
//                           nb, per flavour: gen twin of the skim's dilepton
//                           selection -- lead pT > 15, sub > 10, |eta| < 2.4,
//                           60 < m_ll < 120. The fit has ONE global DY scale
//                           r_Z, so sigma_Z = r_Z * h_gen_sig_Z is the Z-axis
//                           of the (sigma_W, sigma_Z) plane; nothing to bin.)
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
#include "lhe_index.h"   // pOLhe: the EPPS21 member weights (gen-level variations)

#include "TFile.h"
#include "TInterpreter.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TSystem.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

namespace {

// Gen-level fiducial lepton-pT cut = the skim's nominal leading-lepton cut.
constexpr double kFidPtMin = 25.0;

// ---- Z -> ll gen fiducial (2026-09-15) --------------------------------------
// The gen twin of the skim's Z dilepton selection (skim.C: skim_Zmm/skim_Zee
// `ptMin1/ptMin2/etaMax/massMin/massMax`) -- keep these in step with it, they
// are what makes sigma_Z = r_Z * sigma_gen-fid,Z the analog of the W's
// sigma_i = r_i * sigma_gen-fid,i. Both legs sit inside the ntuple gen filter
// (pT > 5, |eta| < 2.5), so the fiducial is fully populated.
constexpr double kZPtLead  = 15.0;
constexpr double kZPtSub   = 10.0;
constexpr double kZEtaMax  = 2.4;
constexpr double kZMassLo  = 60.0;
constexpr double kZMassHi  = 120.0;

// ---- EPPS21 member variations of the GEN cross sections (2026-09-15c) -------
// The same 107 members the skim stores as reco twins (skim/lhe_index.h:
// the LHAPDF set EPPS21nlo_CT18Anlo_O16, member for member), applied here at
// GEN level so the (sigma_W, sigma_Z) plane can show the nPDF coverage of the
// PREDICTION -- what plotting/xsec_contour.C scatters next to the r = 1 point.
//
// NORMALIZATION -- the one choice that matters. The nominal uses
//     sigma_i = kA_O * kSigma * Sumw_i / Sumw_all,
// a weighted FRACTION of a FIXED total. For a member that would be wrong: a
// PDF variation changes the total cross section too, and dividing by that
// member's own Sumw_all would divide exactly that change out, leaving only the
// shape -- i.e. every member would land at the same sigma and the scatter
// would collapse onto a line. So the member cross sections are normalized to
// the NOMINAL total:
//     sigma_i(m) = kA_O * kSigma * Sumw_i(m) / Sumw_all(0),
// which makes sigma_total(m)/sigma_total(0) = Sumw_all(m)/Sumw_all(0) exactly,
// and reproduces the nominal identically for m = 0 (checked and printed).
constexpr int kNEp = pOLhe::kNMembers[pOLhe::kEpps21];   // 107

// Per-member accumulators, in plain arrays: 7M events x 107 members is far too
// many TH2D::Fill calls, so sum here and fill the histograms once at the end.
struct MemAccum
{
  double lab[pOSkim::kNY][kNEp];
  double fb [pOSkim::kNY][kNEp];
  double fid[kNEp];                 // the Z's single fiducial bin
  double all[kNEp];                 // over ALL events (the member total sigma)
  double allNom = 0;                // Sumw over all events, member 0 == nominal
  long long nNoLhe = 0;             // events without a usable ttbar_w
  unsigned long long nNeutral = 0;  // events the kMaxMemberRatio guard flattened
  bool warned = false;
  MemAccum()
  {
    for (int i = 0; i < pOSkim::kNY; ++i)
      for (int m = 0; m < kNEp; ++m) { lab[i][m] = 0; fb[i][m] = 0; }
    for (int m = 0; m < kNEp; ++m) { fid[m] = 0; all[m] = 0; }
  }
};

// Accumulate one W MC file into the (shared) weighted gen-eta histograms:
// hLab/hFB get the FIDUCIAL fill (pT > kFidPtMin), hLabTot/hFBTot every lepton.
bool AccumulateGen(const char *fname, int flavPdg, int wantChg,
                   TH1D *hLab, TH1D *hFB, TH1D *hLabTot, TH1D *hFBTot,
                   double &sumw, long long &nraw,
                   long long &nNoLep, long long &nNoWAnc,
                   MemAccum *ma)
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
  std::vector<float> *ttbar_w = nullptr;
  tHi->SetBranchStatus("*", 0);
  tHi->SetBranchStatus("weight", 1);     tHi->SetBranchAddress("weight", &weight);
  const bool hasLhe = ma && pOSkim::HasBranch(tHi, "ttbar_w");
  if (hasLhe) { tHi->SetBranchStatus("ttbar_w", 1); tHi->SetBranchAddress("ttbar_w", &ttbar_w); }
  else if (ma)
    std::cout << "[gen_xsec][WARN] no ttbar_w in " << fname
              << " -- EPPS21 member cross sections not filled from this file\n";

  pOLhe::MemberWeights mw;
  for (Long64_t ie = 0; ie < n; ++ie)
  {
    tGen->GetEntry(ie);
    tHi->GetEntry(ie);
    const double w = (double)weight;
    sumw += w;
    ++nraw;
    // member weights once per event; `okMem` false -> this event contributes to
    // the nominal only (counted), never a silent zero
    bool okMem = false;
    if (hasLhe)
    {
      okMem = pOLhe::ComputeMemberWeights(w, ttbar_w, mw, ma->warned, "gen_xsec", &ma->nNeutral);
      if (!okMem) ++ma->nNoLhe;
      else { ma->allNom += w; for (int m = 0; m < kNEp; ++m) ma->all[m] += mw.w[pOLhe::kEpps21][m]; }
    }
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
      if (okMem)
      {
        // bin index from the histogram itself -- never a second copy of the edges
        const int bl = hLab->GetXaxis()->FindFixBin(y), bf = hFB->GetXaxis()->FindFixBin(y);
        for (int m = 0; m < kNEp; ++m)
        {
          const double wm = mw.w[pOLhe::kEpps21][m];
          if (bl >= 1 && bl <= pOSkim::kNY) ma->lab[bl - 1][m] += wm;
          if (bf >= 1 && bf <= pOSkim::kNY) ma->fb [bf - 1][m] += wm;
        }
      }
    }
  }
  f->Close();
  delete f;
  return true;
}

// Accumulate one DY MC file into the Z gen-fiducial weighted sums: sumwFid over
// events with an OS same-flavour gen pair inside the fiducial box above, sumw
// over ALL events (the cross-section denominator, exactly as for the W).
// Pair choice: the highest-pT positive and highest-pT negative gen lepton of
// the sample's flavour. No Z-ancestor preference is attempted -- the filtered
// gen collection `HiGenParticleAna/hi` stores no boson (every lepton has
// motherIdx = -999, so pOSkim::HasAncestor is a structural no-op there; see
// CLAUDE.md, correction/charge_flip.C), and for an exclusive DY sample the two
// hardest same-flavour leptons are the decay products anyway.
bool AccumulateGenZ(const char *fname, int flavPdg, double lepMass,
                    double &sumwFid, double &sumw, long long &nraw,
                    long long &nNoPair, MemAccum *ma)
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

  std::vector<float> *pt = nullptr, *eta = nullptr, *phi = nullptr;
  std::vector<int> *chg = nullptr, *pdg = nullptr;
  tGen->SetBranchStatus("*", 0);
  tGen->SetBranchStatus("pt", 1);   tGen->SetBranchAddress("pt", &pt);
  tGen->SetBranchStatus("eta", 1);  tGen->SetBranchAddress("eta", &eta);
  tGen->SetBranchStatus("phi", 1);  tGen->SetBranchAddress("phi", &phi);
  tGen->SetBranchStatus("chg", 1);  tGen->SetBranchAddress("chg", &chg);
  tGen->SetBranchStatus("pdg", 1);  tGen->SetBranchAddress("pdg", &pdg);
  Float_t weight = 1.f;
  std::vector<float> *ttbar_w = nullptr;
  tHi->SetBranchStatus("*", 0);
  tHi->SetBranchStatus("weight", 1); tHi->SetBranchAddress("weight", &weight);
  const bool hasLhe = ma && pOSkim::HasBranch(tHi, "ttbar_w");
  if (hasLhe) { tHi->SetBranchStatus("ttbar_w", 1); tHi->SetBranchAddress("ttbar_w", &ttbar_w); }
  else if (ma)
    std::cout << "[gen_xsec][WARN] no ttbar_w in " << fname
              << " -- EPPS21 member cross sections not filled from this file\n";

  pOLhe::MemberWeights mw;
  for (Long64_t ie = 0; ie < n; ++ie)
  {
    tGen->GetEntry(ie);
    tHi->GetEntry(ie);
    const double w = (double)weight;
    sumw += w;
    ++nraw;
    bool okMem = false;
    if (hasLhe)
    {
      okMem = pOLhe::ComputeMemberWeights(w, ttbar_w, mw, ma->warned, "gen_xsec Z", &ma->nNeutral);
      if (!okMem) ++ma->nNoLhe;
      else { ma->allNom += w; for (int m = 0; m < kNEp; ++m) ma->all[m] += mw.w[pOLhe::kEpps21][m]; }
    }
    if (!pt || !eta || !phi || !chg || !pdg) continue;
    const size_t ng = std::min({pt->size(), eta->size(), phi->size(),
                                chg->size(), pdg->size()});
    int iPos = -1, iNeg = -1;
    for (size_t i = 0; i < ng; ++i)
    {
      if (std::abs(pdg->at(i)) != flavPdg) continue;
      if (chg->at(i) > 0) { if (iPos < 0 || pt->at(i) > pt->at(iPos)) iPos = (int)i; }
      else if (chg->at(i) < 0) { if (iNeg < 0 || pt->at(i) > pt->at(iNeg)) iNeg = (int)i; }
    }
    if (iPos < 0 || iNeg < 0) { ++nNoPair; continue; }
    const double pt1 = pt->at(iPos), pt2 = pt->at(iNeg);
    const double lead = std::max(pt1, pt2), sub = std::min(pt1, pt2);
    if (lead < kZPtLead || sub < kZPtSub) continue;
    if (std::abs(eta->at(iPos)) > kZEtaMax || std::abs(eta->at(iNeg)) > kZEtaMax) continue;
    TLorentzVector v1, v2;
    v1.SetPtEtaPhiM(pt->at(iPos), eta->at(iPos), phi->at(iPos), lepMass);
    v2.SetPtEtaPhiM(pt->at(iNeg), eta->at(iNeg), phi->at(iNeg), lepMass);
    const double m = (v1 + v2).M();
    if (m < kZMassLo || m > kZMassHi) continue;
    sumwFid += w;
    if (okMem)
      for (int k = 0; k < kNEp; ++k) ma->fid[k] += mw.w[pOLhe::kEpps21][k];
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

  // kept for the plain-text sidecar written at the end (see below)
  TH1D *keepLab[2] = {nullptr, nullptr}, *keepFB[2] = {nullptr, nullptr};

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
    MemAccum ma;
    const char *flavs[2] = {"mu", "ele"};
    for (int fl = 0; fl < 2; ++fl)
    {
      const int flavPdg = (fl == 0) ? 13 : 11;
      SampleFileInfo info = ResolveMCSample(samp[ic], flavs[fl]);
      printf("[gen_xsec] %s %-3s <- %s\n", cname[ic], flavs[fl], info.fname.c_str());
      AccumulateGen(info.fname.c_str(), flavPdg, wantChg[ic], hLab, hFB,
                    hLabTot, hFBTot, sumw, nraw, nNoLep, nNoWAnc, &ma);
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
    keepLab[ic] = hLab; keepFB[ic] = hFB;

    // ---- EPPS21 member twins of the gen fiducial sigma (2026-09-15c) -------
    // x = the rapidity bin (axis copied from the nominal), y = member index.
    // Normalized to the NOMINAL total (see the header of MemAccum) so a
    // member's own total-sigma change survives -- that spread IS the coverage
    // plotting/xsec_contour.C scatters.
    if (ma.allNom > 0)
    {
      const double sc = sigAll / ma.allNom;
      TH2D *tw[2];
      TH1D *src[2] = {hLab, hFB};
      const char *tag[2] = {"", "_FB"};
      for (int v = 0; v < 2; ++v)
      {
        tw[v] = new TH2D(Form("h_gen_sig_%s%s_epps21", cname[ic], tag[v]),
                         Form("gen fiducial #sigma_{i} per EPPS21 member, %s%s"
                              " (%s);#eta^{l}_{lab};member;#sigma_{i} (nb)",
                              cname[ic], tag[v], pOLhe::kEpps21SetName),
                         kNY, src[v]->GetXaxis()->GetXbins()->GetArray(),
                         kNEp, -0.5, kNEp - 0.5);
        for (int i = 0; i < kNY; ++i)
          for (int m = 0; m < kNEp; ++m)
            tw[v]->SetBinContent(i + 1, m + 1,
                                 sc * (v == 0 ? ma.lab[i][m] : ma.fb[i][m]));
        fout->cd();
        tw[v]->Write("", TObject::kOverwrite);
      }
      // member 0 must reproduce the nominal histogram exactly
      double dmax = 0;
      for (int i = 0; i < kNY; ++i)
        dmax = std::max(dmax, std::fabs(tw[0]->GetBinContent(i + 1, 1) - hLab->GetBinContent(i + 1)));
      double lo = 1e30, hi = -1e30;
      for (int m = 0; m < kNEp; ++m) { lo = std::min(lo, ma.all[m]); hi = std::max(hi, ma.all[m]); }
      printf("[gen_xsec] %s: EPPS21 members -> sigma_total spread %.2f..%.2f nb"
             " (nominal %.2f); max|member0 - nominal| per bin = %.3e\n",
             cname[ic], sigAll * lo / ma.allNom, sigAll * hi / ma.allNom, sigAll, dmax);
      if (ma.nNoLhe > 0)
        printf("[gen_xsec][WARN] %s: %lld events without a usable ttbar_w (member sums miss them)\n",
               cname[ic], ma.nNoLhe);
      if (ma.nNeutral > 0)
        printf("[gen_xsec] %s: %llu events flattened to the nominal by the"
               " kMaxMemberRatio = %.0f guard\n", cname[ic], ma.nNeutral, pOLhe::kMaxMemberRatio);
    }
    else
      printf("[gen_xsec][WARN] %s: no LHE member weights -> h_gen_sig_%s*_epps21 NOT written\n",
             cname[ic], cname[ic]);
  }

  double sigZfid = -1.0; // for the plain-text sidecar below

  // ---- Z -> ll gen fiducial sigma (2026-09-15) ------------------------------
  // sigma_Z^fid = kA_O * kSigma_DY * Sumw_fid / Sumw_all, per lepton flavour --
  // the exact analog of the W above, so the fitted global DY scale converts as
  // sigma_Z = r_Z * sigma_gen-fid,Z. The mu and e files are pooled like the W's
  // (lepton universality; their individual acceptances are printed, and they
  // differ slightly because bare post-FSR electrons lose more energy).
  // One bin, inclusive: the fit has ONE r_Z, so there is nothing to bin in.
  {
    fout->cd();
    TH1D *hZ = new TH1D("h_gen_sig_Z",
                        Form("gen fiducial #sigma(Z#rightarrow ll), p_{T}>%.0f/%.0f,"
                             " |#eta|<%.1f, %.0f<m_{ll}<%.0f (per flavour);;#sigma (nb)",
                             kZPtLead, kZPtSub, kZEtaMax, kZMassLo, kZMassHi),
                        1, 0., 1.);
    const double sigAllZ = pONorm::kA_O * pONorm::kSigma_DY; // DY->ll, m>50, pO-scaled
    double sumwFidTot = 0, sumwTot = 0;
    long long nrawTot = 0, nNoPairTot = 0;
    MemAccum maZ;
    const char *flavs[2] = {"mu", "ele"};
    const int   flavPdg[2] = {13, 11};
    const double lepM[2] = {MU_MASS, ELE_MASS};
    for (int fl = 0; fl < 2; ++fl)
    {
      SampleFileInfo info = ResolveMCSample(kDY, flavs[fl]);
      printf("[gen_xsec] Z  %-3s <- %s\n", flavs[fl], info.fname.c_str());
      double sumwFid = 0, sumw = 0; long long nraw = 0, nNoPair = 0;
      if (!AccumulateGenZ(info.fname.c_str(), flavPdg[fl], lepM[fl],
                          sumwFid, sumw, nraw, nNoPair, &maZ))
        continue;
      printf("[gen_xsec] Z  %-3s: Nraw = %lld, acceptance = %.4f,"
             " sigma_fid = %.4f nb%s\n",
             flavs[fl], nraw, sumw > 0 ? sumwFid / sumw : 0.0,
             sumw > 0 ? sigAllZ * sumwFid / sumw : 0.0,
             nNoPair > 0 ? Form(" (%lld events without an OS gen pair)", nNoPair) : "");
      sumwFidTot += sumwFid; sumwTot += sumw; nrawTot += nraw; nNoPairTot += nNoPair;
    }
    if (nrawTot > 0 && sumwTot > 0)
    {
      const double acc = sumwFidTot / sumwTot;
      // MC-stat error on the weighted fraction (binomial on the pooled Sumw)
      const double accErr = std::sqrt(std::max(0.0, acc * (1.0 - acc) / nrawTot));
      hZ->SetBinContent(1, sigAllZ * acc);
      hZ->SetBinError(1, sigAllZ * accErr);
      sigZfid = sigAllZ * acc;
      printf("[gen_xsec] Z  : Nraw = %lld, sigma(DY->ll, m>50) = %.3f nb,"
             " acceptance = %.4f, sigma_fid = %.4f +/- %.4f nb (per flavour, pO-scaled)\n",
             nrawTot, sigAllZ, acc, sigAllZ * acc, sigAllZ * accErr);
      if (nNoPairTot > 0)
        printf("[gen_xsec][WARN] Z: %lld events without an OS same-flavour gen pair (skipped)\n",
               nNoPairTot);
      fout->cd();
      hZ->Write("", TObject::kOverwrite);

      // EPPS21 member twins of sigma_Z (1D over members; same normalization)
      if (maZ.allNom > 0)
      {
        TH1D *hZm = new TH1D("h_gen_sig_Z_epps21",
                             Form("gen fiducial #sigma(Z#rightarrow ll) per EPPS21 member"
                                  " (%s);member;#sigma (nb)", pOLhe::kEpps21SetName),
                             kNEp, -0.5, kNEp - 0.5);
        const double scZ = sigAllZ / maZ.allNom;
        for (int m = 0; m < kNEp; ++m) hZm->SetBinContent(m + 1, scZ * maZ.fid[m]);
        double lo = 1e30, hi = -1e30;
        for (int m = 0; m < kNEp; ++m) { lo = std::min(lo, maZ.all[m]); hi = std::max(hi, maZ.all[m]); }
        printf("[gen_xsec] Z  : EPPS21 members -> sigma_fid spread %.4f..%.4f nb"
               " (nominal %.4f); |member0 - nominal| = %.3e\n",
               hZm->GetMinimum() > 0 ? hZm->GetMinimum() : 0.0, hZm->GetMaximum(),
               sigZfid, std::fabs(hZm->GetBinContent(1) - sigZfid));
        printf("[gen_xsec] Z  : member total-sigma spread %.3f..%.3f nb (nominal %.3f)\n",
               sigAllZ * lo / maZ.allNom, sigAllZ * hi / maZ.allNom, sigAllZ);
        if (maZ.nNoLhe > 0)
          printf("[gen_xsec][WARN] Z: %lld events without a usable ttbar_w\n", maZ.nNoLhe);
        if (maZ.nNeutral > 0)
          printf("[gen_xsec] Z  : %llu events flattened by the kMaxMemberRatio guard\n", maZ.nNeutral);
        fout->cd();
        hZm->Write("", TObject::kOverwrite);
      }
      else
        printf("[gen_xsec][WARN] Z: no LHE member weights -> h_gen_sig_Z_epps21 NOT written\n");
    }
    else
      std::cerr << "[ERROR] no DY events accumulated -- h_gen_sig_Z not written\n";
  }

  // ---- plain-text sidecar output/gen_xsec_fid.txt (2026-09-15) -------------
  // Same numbers as the histograms, in a form the FORK's card generator can
  // read with awk: the reparametrized workspace (sigma_W promoted to a POI for
  // the profiled (sigma_W, sigma_Z) contour) has to bake the per-bin gen
  // fiducial cross sections in as constants, because sigma_W = Sum_i r_i
  // sigma_gen,i is a specific linear combination of the POIs. Uploaded next to
  // the Combine inputs by the fork's sync_lxplus.sh, like the _systs.txt
  // sidecars. NB the weights are the NOMINAL-PDF gen cross sections -- the same
  // approximation xsec_fiducial_comb already makes.
  {
    gSystem->mkdir("output", kTRUE);
    std::ofstream sc("output/gen_xsec_fid.txt");
    sc << "# gen FIDUCIAL cross sections, nb, PER LEPTON FLAVOUR -- skim/gen_xsec.C\n"
       << "# W: bare post-FSR lepton pT > " << kFidPtMin << ", |eta_lab| < 2.4 (bin edges = pOSkim::kYEdges[_FB])\n"
       << "# Z: lead pT > " << kZPtLead << ", sub > " << kZPtSub << ", |eta| < " << kZEtaMax
       << ", " << kZMassLo << " < m_ll < " << kZMassHi << "\n"
       << "# columns: binning charge ybin sigma_nb   (binning: lab | fb | incl)\n";
    for (int ic = 0; ic < 2; ++ic)
    {
      if (keepLab[ic])
        for (int i = 0; i < kNY; ++i)
          sc << "lab " << cname[ic] << " " << i << " "
             << Form("%.6f", keepLab[ic]->GetBinContent(i + 1)) << "\n";
      if (keepFB[ic])
        for (int i = 0; i < kNY; ++i)
          sc << "fb " << cname[ic] << " " << i << " "
             << Form("%.6f", keepFB[ic]->GetBinContent(i + 1)) << "\n";
    }
    if (sigZfid > 0) sc << "incl Z 0 " << Form("%.6f", sigZfid) << "\n";
    sc.close();
    std::cout << "[gen_xsec] wrote output/gen_xsec_fid.txt (the fork's card generator reads this)\n";
  }

  fout->Close();
  delete fout;
  std::cout << "[gen_xsec] wrote rootfile/gen_xsec.root\n";
}
