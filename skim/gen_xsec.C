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
// W GEN SOURCE (2026-10-13): the EventTree mc* block (ggHiNtuplizer), NOT
// HiGenParticleAna/hi. The mc* block stores the W itself (status 62) and its
// decay lepton and neutrino with direct parentage (mcMomPID = +-24) and NO
// pT/eta filter; hi is filtered at pT > 5, |eta| < 2.5 and stores no boson
// (every lepton has motherIdx = -999). For the pT/eta fiducial the two give the
// same lepton (verified 2026-09-01; the base set below must reproduce the
// pre-2026-10-13 numbers), but the m_T cut needs the NEUTRINO, which hi drops
// whenever it is forward: 22% (W+) / 9% (W-) of the events with a fiducial
// m_T > 40 lepton. Lepton = the status-1 (bare, post-FSR) lepton of the
// sample's flavour and charge whose mother is the W (or, for an FSR copy the
// ntuplizer keeps separately, ~0.1% of events, whose grandmother is);
// neutrino = the same with |pdg| = flavour + 1 (exactly one of each per event,
// measured). The Z below stays on hi: it needs no neutrino, and hi is complete
// inside its filter, which contains the whole Z fiducial.
//
// FIDUCIAL definition (2026-08-12; two SETS since 2026-10-13). Both sets apply
// the gen-level twin of the reco selection cuts -- bare lepton pT > kFidPtMin
// (= the skim's nominal leading-lepton cut, 25 GeV) with |eta_lab| < 2.4
// implicit in the binning window -- and differ only in the m_T cut:
//     ""      no m_T cut                                    -> the met fit
//     "_mt40" m_T(l, nu) > pOSkim::kPtDiscMtCut (40 GeV)    -> the leppt_mt40 fit
// with m_T from the W's own lepton and neutrino (not a "gen MET": that would
// add neutrinos from heavy-flavour decays and the UE, which no fixed-order
// calculation contains). Each discriminant is measured in the fiducial volume
// of its own reco selection (plotting/disc_variants.h::GenFid is the map), so
// the m_T > 40 acceptance (~0.96) is no longer extrapolated with the MC for
// leppt_mt40; what remains is the reco-gen migration across 40 GeV inside eps.
// sigma_i is a FIDUCIAL cross section, so the fitted signal strength converts
// directly: sigma_meas,i = r_i * sigma_i^gen,fid (= N_fit,i / L / (A*eps)_MC,i
// -- algebraically identical to the yield-based extraction with MC eff/acc,
// but kA_O and kSigma_* cancel between r and sigma_gen).
//
// Output: rootfile/gen_xsec.root, <set> = "" | "_mt40"
//   h_gen_sig_{Wp,Wm}<set>      (lab edges; bin content = per-flavour FIDUCIAL
//                                sigma_i in nb)
//   h_gen_sig_{Wp,Wm}<set>_FB   (FB edges;  same content convention)
//   h_gen_sig_{Wp,Wm}<set>[_FB]_epps21  (the 107 EPPS21 members, see MemAccum)
//   h_gen_tot_{Wp,Wm}[_FB]      (reference: NO pT and no m_T cut, the lepton in
//                                the |eta_lab| < 2.4 window. Before 2026-10-13
//                                this came from the FILTERED hi tree, i.e. it
//                                was a pT > 5 sigma; now it is the true one)
//   h_gen_sig_Z                 (2026-09-15; ONE bin = the Z->ll FIDUCIAL sigma
//                                in nb, per flavour: gen twin of the skim's
//                                dilepton selection -- lead pT > 15, sub > 10,
//                                |eta| < 2.4, 60 < m_ll < 120. The fit has ONE
//                                global DY scale r_Z, so sigma_Z = r_Z *
//                                h_gen_sig_Z is the Z-axis of the
//                                (sigma_W, sigma_Z) plane; nothing to bin.
//                                The same for both W sets.)
// plus the plain-text sidecars output/gen_xsec_fid<set>.txt (see the end).
//
// The Sum(w) denominator runs over ALL events. Sumw2 carries the MC-stat
// error. Consumed by analysis/fiducial_yields.C and plotting/xsec_{fiducial,
// contour,summary}.C, each picking the set of its discriminant.
//
//   cd skim/ && ./run_gen_xsec.sh        (log: logs/gen_xsec.log)
// =============================================================================
#include "skim_common.h"
#include "mc_norm.h"
#include "lhe_index.h"   // pOLhe: the EPPS21 member weights (gen-level variations)

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TSystem.h"
#include "TVector2.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

namespace {

// Gen-level fiducial lepton-pT cut = the skim's nominal leading-lepton cut.
constexpr double kFidPtMin = 25.0;

// The W fiducial SETS (2026-10-13; see the header). The tag is appended to the
// charge in every histogram name (h_gen_sig_Wp_mt40, h_gen_sig_Wp_mt40_FB,
// h_gen_sig_Wp_mt40_FB_epps21) and to the sidecar name; set s applies the
// m_T cut iff kFidMt[s]. plotting/disc_variants.h::GenFid maps a discriminant
// to its tag -- keep the two in step.
constexpr int kNFid = 2;
const char *const kFidTag[kNFid] = {"", "_mt40"};
const bool        kFidMt[kNFid]  = {false, true};

// the cuts of set s, for titles and printouts (from the constants, never typed)
TString FidCuts(int s, bool latex)
{
  TString t = TString::Format(latex ? "p_{T}>%.0f" : "pT > %.0f", kFidPtMin);
  if (kFidMt[s]) t += TString::Format(latex ? ", m_{T}>%.0f" : ", mT > %.0f", pOSkim::kPtDiscMtCut);
  return t;
}

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
  double lab[kNFid][pOSkim::kNY][kNEp]; // W, per fiducial set
  double fb [kNFid][pOSkim::kNY][kNEp];
  double fid[kNEp];                 // the Z's single fiducial bin
  double all[kNEp];                 // over ALL events (the member total sigma)
  double allNom = 0;                // Sumw over all events, member 0 == nominal
  long long nNoLhe = 0;             // events without a usable ttbar_w
  unsigned long long nNeutral = 0;  // events the kMaxMemberRatio guard flattened
  bool warned = false;
  MemAccum()
  {
    for (int s = 0; s < kNFid; ++s)
      for (int i = 0; i < pOSkim::kNY; ++i)
        for (int m = 0; m < kNEp; ++m) { lab[s][i][m] = 0; fb[s][i][m] = 0; }
    for (int m = 0; m < kNEp; ++m) { fid[m] = 0; all[m] = 0; }
  }
};

// The W daughter of |pdg| == absPdg in the EventTree mc* block (2026-10-13):
// status 1 (for a charged lepton = the BARE, post-FSR one), mother the W, or a
// same-flavour copy whose mother is the W (the ntuplizer keeps the pre-FSR
// lepton as a separate entry in ~0.1% of events; the status-1 copy then has
// mcMomPID = the lepton and mcGMomPID = the W). signedPdg != 0 also fixes the
// charge (a charged lepton with pdg > 0 is NEGATIVE). Highest pT if several;
// -1 if none.
int FindWDaughter(const std::vector<int> &pid, const std::vector<int> &st,
                  const std::vector<int> &mom, const std::vector<int> &gmom,
                  const std::vector<float> &pt, int absPdg, int signedPdg)
{
  const size_t n = std::min({pid.size(), st.size(), mom.size(), gmom.size(), pt.size()});
  int best = -1;
  for (size_t i = 0; i < n; ++i)
  {
    if (std::abs(pid[i]) != absPdg) continue;
    if (signedPdg != 0 && pid[i] != signedPdg) continue;
    if (st[i] != 1) continue;
    const int m = std::abs(mom[i]);
    if (m != 24 && !(m == absPdg && std::abs(gmom[i]) == 24)) continue;
    if (best < 0 || pt[i] > pt[best]) best = (int)i;
  }
  return best;
}

// Accumulate one W MC file into the (shared) weighted gen-eta histograms:
// hLab[s]/hFB[s] get the FIDUCIAL fill of set s, hLabTot/hFBTot every lepton
// (no pT, no m_T cut). Reads the EventTree mc* block -- see the header.
bool AccumulateGen(const char *fname, int flavPdg, int wantChg,
                   TH1D *const *hLab, TH1D *const *hFB, TH1D *hLabTot, TH1D *hFBTot,
                   double &sumw, long long &nraw,
                   long long &nNoLep, long long &nNoNu,
                   MemAccum *ma)
{
  TFile *f = TFile::Open(fname, "READ");
  if (!f || f->IsZombie()) { std::cerr << "[ERROR] cannot open " << fname << "\n"; return false; }
  TTree *tGen = (TTree *)f->Get("ggHiNtuplizer/EventTree");
  TTree *tHi  = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  if (!tGen || !tHi)
  { std::cerr << "[ERROR] missing ggHiNtuplizer/EventTree or hiEvtAnalyzer/HiTree in " << fname << "\n"; f->Close(); return false; }
  const Long64_t n = tGen->GetEntries();
  if (tHi->GetEntries() != n)
  { std::cerr << "[ERROR] EventTree/HiTree entry mismatch in " << fname << "\n"; f->Close(); return false; }

  // fast-skim discipline: disable everything, enable only what is read (a
  // missing mc* branch is FATAL -- never a silently empty gen collection)
  std::vector<int> *mcPID = nullptr, *mcStatus = nullptr, *mcMomPID = nullptr, *mcGMomPID = nullptr;
  std::vector<float> *mcPt = nullptr, *mcEta = nullptr, *mcPhi = nullptr;
  tGen->SetBranchStatus("*", 0);
  const char *bInt[4] = {"mcPID", "mcStatus", "mcMomPID", "mcGMomPID"};
  std::vector<int> **pInt[4] = {&mcPID, &mcStatus, &mcMomPID, &mcGMomPID};
  const char *bFlt[3] = {"mcPt", "mcEta", "mcPhi"};
  std::vector<float> **pFlt[3] = {&mcPt, &mcEta, &mcPhi};
  for (int b = 0; b < 7; ++b)
  {
    const char *bn = (b < 4) ? bInt[b] : bFlt[b - 4];
    if (!pOSkim::HasBranch(tGen, bn))
    { std::cerr << "[FATAL] gen_xsec: no branch " << bn << " in the EventTree of " << fname << "\n"; f->Close(); return false; }
    tGen->SetBranchStatus(bn, 1);
    if (b < 4) tGen->SetBranchAddress(bn, pInt[b]);
    else       tGen->SetBranchAddress(bn, pFlt[b - 4]);
  }
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
    if (!mcPID || !mcStatus || !mcMomPID || !mcGMomPID || !mcPt || !mcEta || !mcPhi) { ++nNoLep; continue; }
    // the W's charged lepton (pdg < 0 is the POSITIVE one) and its neutrino
    const int iL = FindWDaughter(*mcPID, *mcStatus, *mcMomPID, *mcGMomPID, *mcPt,
                                 flavPdg, -wantChg * flavPdg);
    if (iL < 0) { ++nNoLep; continue; } // no W lepton of this flavour+charge (expected: never)
    const int iN = FindWDaughter(*mcPID, *mcStatus, *mcMomPID, *mcGMomPID, *mcPt,
                                 flavPdg + 1, 0);
    if (iN < 0) ++nNoNu;                // then it fails the m_T sets only (expected: never)
    const double lepPt = mcPt->at(iL);
    // [NOTICE] Sign flip: p-going (-Z) defined as forward -- MUST match the
    // reco convention in skim.C (`const double y = -{mu,ele}Eta->at(iLead)`),
    // or gen bin i pairs with the MIRRORED reco region y_i. Caught 2026-08-12
    // by the per-bin (A x eps) diagnostic: mirrored pairing made A x eps
    // charge-dependent (mu W- ran 1.27 -> 0.73 across eta) instead of the
    // charge-symmetric detector response it must be.
    const double y = -mcEta->at(iL);
    hLabTot->Fill(y, w);
    hFBTot->Fill(y, w);
    if (lepPt <= kFidPtMin) continue;
    // m_T of the W's own lepton and neutrino (the same formula as the reco
    // m_T in skim.C, with the neutrino in place of the PF MET)
    TVector2 nu;
    if (iN >= 0) nu.SetMagPhi(mcPt->at(iN), mcPhi->at(iN));
    const double mt = (iN >= 0) ? pOSkim::TransverseMass(lepPt, mcPhi->at(iL), nu) : -1.0;
    // bin index from the histogram itself -- never a second copy of the edges
    const int bl = hLab[0]->GetXaxis()->FindFixBin(y), bf = hFB[0]->GetXaxis()->FindFixBin(y);
    for (int s = 0; s < kNFid; ++s)
    {
      if (kFidMt[s] && !(mt > pOSkim::kPtDiscMtCut)) continue;
      hLab[s]->Fill(y, w);
      hFB[s]->Fill(y, w);
      if (okMem)
        for (int m = 0; m < kNEp; ++m)
        {
          const double wm = mw.w[pOLhe::kEpps21][m];
          if (bl >= 1 && bl <= pOSkim::kNY) ma->lab[s][bl - 1][m] += wm;
          if (bf >= 1 && bf <= pOSkim::kNY) ma->fb [s][bf - 1][m] += wm;
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
  // (no vector<vector<int>> dictionary any more: the W reads the flat mc*
  //  block since 2026-10-13, and the Z never read hi's motherIdx)

  const char *cname[2]     = {"Wp", "Wm"};
  const SampleType samp[2] = {kWp, kWm};
  const int wantChg[2]     = {+1, -1};
  const double sigNN[2]    = {pONorm::kSigma_Wp, pONorm::kSigma_Wm}; // nb, single source

  TFile *fout = TFile::Open("rootfile/gen_xsec.root", "RECREATE");
  printf("[gen_xsec] kA_O = %.1f, sigma_NN(Wp/Wm) = %.3f/%.3f nb (mc_norm.h);"
         " sigma_i = A * sigma * Sumw_i/Sumw_all\n",
         pONorm::kA_O, sigNN[0], sigNN[1]);

  // kept for the plain-text sidecars written at the end (see below)
  TH1D *keepLab[2][kNFid] = {}, *keepFB[2][kNFid] = {};

  for (int ic = 0; ic < 2; ++ic)
  {
    fout->cd();
    TH1D *hLab[kNFid], *hFB[kNFid];
    for (int s = 0; s < kNFid; ++s)
    {
      const TString cuts = FidCuts(s, true);
      hLab[s] = new TH1D(Form("h_gen_sig_%s%s", cname[ic], kFidTag[s]),
                         Form("gen fiducial d#sigma bins (%s), %s (per flavour);#eta^{l}_{lab};#sigma_{i} (nb)", cuts.Data(), cname[ic]),
                         kNY, kYEdges);
      hFB[s] = new TH1D(Form("h_gen_sig_%s%s_FB", cname[ic], kFidTag[s]),
                        Form("gen fiducial d#sigma bins (FB edges, %s), %s (per flavour);#eta^{l}_{lab};#sigma_{i} (nb)", cuts.Data(), cname[ic]),
                        kNY, kYEdgesFB);
    }
    TH1D *hLabTot = new TH1D(Form("h_gen_tot_%s", cname[ic]),
                             Form("gen d#sigma bins (no p_{T}, no m_{T} cut), %s (per flavour);#eta^{l}_{lab};#sigma_{i} (nb)", cname[ic]),
                             kNY, kYEdges);
    TH1D *hFBTot = new TH1D(Form("h_gen_tot_%s_FB", cname[ic]),
                            Form("gen d#sigma bins (FB edges, no p_{T}, no m_{T} cut), %s (per flavour);#eta^{l}_{lab};#sigma_{i} (nb)", cname[ic]),
                            kNY, kYEdgesFB);

    double sumw = 0; long long nraw = 0, nNoLep = 0, nNoNu = 0;
    MemAccum ma;
    const char *flavs[2] = {"mu", "ele"};
    for (int fl = 0; fl < 2; ++fl)
    {
      const int flavPdg = (fl == 0) ? 13 : 11;
      SampleFileInfo info = ResolveMCSample(samp[ic], flavs[fl]);
      printf("[gen_xsec] %s %-3s <- %s\n", cname[ic], flavs[fl], info.fname.c_str());
      if (!AccumulateGen(info.fname.c_str(), flavPdg, wantChg[ic], hLab, hFB,
                         hLabTot, hFBTot, sumw, nraw, nNoLep, nNoNu, &ma))
      { std::cerr << "[ERROR] " << cname[ic] << " " << flavs[fl] << " not accumulated -- " << cname[ic] << " NOT written\n"; nraw = 0; break; }
    }
    if (nraw <= 0 || sumw <= 0) { std::cerr << "[ERROR] no events accumulated for " << cname[ic] << "\n"; continue; }

    const double sigAll = pONorm::kA_O * sigNN[ic]; // all eta, per flavour, pO-scaled
    const double scale  = sigAll / sumw;            // weighted-fraction normalization
    for (int s = 0; s < kNFid; ++s) { hLab[s]->Scale(scale); hFB[s]->Scale(scale); }
    hLabTot->Scale(scale);
    hFBTot->Scale(scale);

    printf("[gen_xsec] %s: Nraw = %lld, sigma(all eta) = %.2f nb,"
           " sigma(|eta_lab|<2.4) = %.2f nb (per flavour, pO-scaled)\n",
           cname[ic], nraw, sigAll, hLabTot->Integral());
    for (int s = 0; s < kNFid; ++s)
      printf("[gen_xsec] %s: set '%s' sigma_fid(%s, |eta_lab|<2.4) = %.2f nb"
             " (lab), %.2f nb (FB window); /sigma(|eta_lab|<2.4) = %.4f%s\n",
             cname[ic], kFidTag[s], FidCuts(s, false).Data(), hLab[s]->Integral(), hFB[s]->Integral(),
             hLabTot->Integral() > 0 ? hLab[s]->Integral() / hLabTot->Integral() : 0.0,
             s > 0 ? Form(", /set '%s' = %.4f", kFidTag[0], hLab[s]->Integral() / hLab[0]->Integral()) : "");
    if (nNoLep > 0)
      printf("[gen_xsec][WARN] %s: %lld events without a W lepton of charge %+d in the mc* block (skipped)\n",
             cname[ic], nNoLep, wantChg[ic]);
    if (nNoNu > 0)
      printf("[gen_xsec][WARN] %s: %lld events without the W neutrino in the mc* block (fail every m_T set)\n",
             cname[ic], nNoNu);

    fout->cd();
    for (int s = 0; s < kNFid; ++s)
    {
      hLab[s]->Write("", TObject::kOverwrite);
      hFB[s]->Write("", TObject::kOverwrite);
      keepLab[ic][s] = hLab[s]; keepFB[ic][s] = hFB[s];
    }
    hLabTot->Write("", TObject::kOverwrite);
    hFBTot->Write("", TObject::kOverwrite);

    // ---- EPPS21 member twins of the gen fiducial sigma (2026-09-15c) -------
    // x = the rapidity bin (axis copied from the nominal), y = member index.
    // Normalized to the NOMINAL total (see the header of MemAccum) so a
    // member's own total-sigma change survives -- that spread IS the coverage
    // plotting/xsec_contour.C scatters. One pair per fiducial set.
    if (ma.allNom > 0)
    {
      const double sc = sigAll / ma.allNom;
      const char *tag[2] = {"", "_FB"};
      double dmax = 0; // member 0 must reproduce the nominal histograms exactly
      for (int s = 0; s < kNFid; ++s)
        for (int v = 0; v < 2; ++v)
        {
          TH1D *src = (v == 0) ? hLab[s] : hFB[s];
          TH2D *tw = new TH2D(Form("h_gen_sig_%s%s%s_epps21", cname[ic], kFidTag[s], tag[v]),
                              Form("gen fiducial #sigma_{i} per EPPS21 member, %s%s (%s)"
                                   " (%s);#eta^{l}_{lab};member;#sigma_{i} (nb)",
                                   cname[ic], tag[v], FidCuts(s, true).Data(), pOLhe::kEpps21SetName),
                              kNY, src->GetXaxis()->GetXbins()->GetArray(),
                              kNEp, -0.5, kNEp - 0.5);
          for (int i = 0; i < kNY; ++i)
            for (int m = 0; m < kNEp; ++m)
              tw->SetBinContent(i + 1, m + 1,
                                sc * (v == 0 ? ma.lab[s][i][m] : ma.fb[s][i][m]));
          for (int i = 0; i < kNY; ++i)
            dmax = std::max(dmax, std::fabs(tw->GetBinContent(i + 1, 1) - src->GetBinContent(i + 1)));
          fout->cd();
          tw->Write("", TObject::kOverwrite);
        }
      double lo = 1e30, hi = -1e30;
      for (int m = 0; m < kNEp; ++m) { lo = std::min(lo, ma.all[m]); hi = std::max(hi, ma.all[m]); }
      printf("[gen_xsec] %s: EPPS21 members -> sigma_total spread %.2f..%.2f nb"
             " (nominal %.2f); max|member0 - nominal| per bin (all sets, lab+FB) = %.3e\n",
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

  // ---- plain-text sidecars output/gen_xsec_fid<set>.txt (2026-09-15) -------
  // Same numbers as the histograms, in a form the FORK's card generator can
  // read with awk: the reparametrized workspace (sigma_W promoted to a POI for
  // the profiled (sigma_W, sigma_Z) contour) has to bake the per-bin gen
  // fiducial cross sections in as constants, because sigma_W = Sum_i r_i
  // sigma_gen,i is a specific linear combination of the POIs. Uploaded next to
  // the Combine inputs by the fork's sync_lxplus.sh, like the _systs.txt
  // sidecars. NB the weights are the NOMINAL-PDF gen cross sections -- the same
  // approximation xsec_fiducial_comb already makes.
  // ONE FILE PER FIDUCIAL SET since 2026-10-13 (gen_xsec_fid.txt = no m_T cut,
  // gen_xsec_fid_mt40.txt = m_T > 40), identical layout; the fork's
  // run_pO_fits.sh picks the file of the fitted discriminant. The Z row is the
  // same in both.
  gSystem->mkdir("output", kTRUE);
  for (int s = 0; s < kNFid; ++s)
  {
    const TString path = TString::Format("output/gen_xsec_fid%s.txt", kFidTag[s]);
    std::ofstream sc(path.Data());
    sc << "# gen FIDUCIAL cross sections, nb, PER LEPTON FLAVOUR -- skim/gen_xsec.C\n"
       << "# fiducial set '" << kFidTag[s] << "'\n"
       << "# W: bare post-FSR lepton (EventTree mc* block) " << FidCuts(s, false).Data()
       << ", |eta_lab| < 2.4 (bin edges = pOSkim::kYEdges[_FB])"
       << (kFidMt[s] ? "; m_T from the W's lepton and neutrino" : "") << "\n"
       << "# Z: lead pT > " << kZPtLead << ", sub > " << kZPtSub << ", |eta| < " << kZEtaMax
       << ", " << kZMassLo << " < m_ll < " << kZMassHi << "\n"
       << "# columns: binning charge ybin sigma_nb   (binning: lab | fb | incl)\n";
    for (int ic = 0; ic < 2; ++ic)
    {
      if (keepLab[ic][s])
        for (int i = 0; i < kNY; ++i)
          sc << "lab " << cname[ic] << " " << i << " "
             << Form("%.6f", keepLab[ic][s]->GetBinContent(i + 1)) << "\n";
      if (keepFB[ic][s])
        for (int i = 0; i < kNY; ++i)
          sc << "fb " << cname[ic] << " " << i << " "
             << Form("%.6f", keepFB[ic][s]->GetBinContent(i + 1)) << "\n";
    }
    if (sigZfid > 0) sc << "incl Z 0 " << Form("%.6f", sigZfid) << "\n";
    sc.close();
    std::cout << "[gen_xsec] wrote " << path << " (the fork's card generator reads this)\n";
  }

  fout->Close();
  delete fout;
  std::cout << "[gen_xsec] wrote rootfile/gen_xsec.root\n";
}
