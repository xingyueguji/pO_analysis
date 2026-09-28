// correction/trig_eff_mb.C
//
// SINGLE-LEPTON TRIGGER EFFICIENCY (turn-on) of the W analysis path, measured
// on the MINIMUM-BIAS-triggered sample, DATA vs W SIGNAL MC on the same plot.
//
// Definition (per the group's request, 2026-09-09):
//   den = event passes the analysis W selection  &&  MB path fired
//              [MC only, since 2026-09-15: && the leading lepton gen-matches a
//               PROMPT W LEPTON -- see MatchGenLeptonFromW]
//   num = den  &&  analysis single-lepton path fired
//              &&  leading lepton matched to a trigger object (DR < 0.4)
//   eps(pT) = num / den,   SF(pT) = eps_data / eps_MC
//
// The prompt-W gen match is the ONLY difference between the data and MC legs
// (user decision 2026-09-15). It makes eps_MC "the efficiency for a real W
// lepton" by definition instead of by sample purity -- for muons it is a ~0.04%
// effect, since the Wp_mu/Wm_mu samples are W -> mu nu by construction and
// correction/charge_flip.C measured only 0.001% unmatched muons in this same
// selection, but it self-documents the definition and matters for electrons
// (0.7% unmatched there). NB it does NOT symmetrize the two legs: the DATA
// denominator is still a mixture, ~5% QCD fakes with m_T > 40 and ~19% without,
// which is the real asymmetry and the reason mt40 is the nominal variant.
//
// "Analysis W selection" = the skim's 8-step W cutflow (skim/skim.C) WITHOUT
// steps 3 (trigger fired) and 8 (trigger match) -- those two ARE the quantity
// being measured. Steps 1/2/4/5/6/7 are replicated exactly as in
// charge_flip.C / njet_WZ.C (both verified event-identical with the skim).
// The lepton-pT floor of steps 1 and 6 is LOWERED to kPtFloor (10 GeV) so the
// actual turn-on is visible; the analysis selection is the pT > 25 sub-range
// (25 is a bin edge, so the two are the same histogram -- the leading tight
// lepton is the same object whenever it has pT > 25).
//
// Why the MB denominator works:
//   HLT_MinimumBiasHF_OR_BptxAND_v1 is a global HF-activity path. It knows
//   nothing about the lepton, so requiring it selects an UNBIASED subsample of
//   the offline-selected events -- including the ones the lepton path missed,
//   which are exactly the inefficiency. There is no object to match for the
//   MB path ("MB trigger matched" == the MB bit fired). In DATA the path was
//   prescaled run by run (fraction of lepton-triggered events with the MB bit:
//   393952 100%, 393953 91%, 393974/393975 0%, 393976 67%, 394004 51%,
//   394005/6 68%, 394007 78%; 49% overall, identical for muon- and electron-
//   triggered events within each run -> a pure prescale, uncorrelated with the
//   lepton; measured 2026-09-09). In MC it fires on 99.93% of W events (no
//   prescale). The prescale only costs statistics; the two MB-FRACTION control
//   plots below (den/all and num/trg vs pT) verify it is flat in lepton pT.
//
// Two selection variants, both filled in one loop:
//   nom   -- the plain W selection (NO MET / m_T cut; in data the low-pT part
//            is QCD-dominated: ~19% mu / ~50% e of the pT > 25 sample, more
//            below -- so the data turn-on is that of the MIXTURE)
//   mt40  -- + m_T > 40 GeV (PF MET from particleFlowAnalyser/pftree, read
//            only for selected events) = the leppt_mt40 discriminant selection,
//            QCD ~5% mu / ~28% e -> the cleaner data/MC comparison.
//
// What this is NOT: a tag-and-probe measurement. T&P (Z->ll) measures the
// per-LEPTON efficiency on a background-subtracted pure-lepton probe sample
// with Z kinematics; this measures the per-EVENT trigger efficiency of the W
// selection itself (W kinematics, W charge mix, the real fake content in
// data), with no tag bias, at the cost of the MB prescale and the data
// background composition. See the discussion in the README / AN.
//
// Outputs (run from correction/):
//   rootfile/trig_eff_mb_<mu|ele>.root  raw count histos per sample (+ MC
//                                       combined), TEfficiency objects, SF graphs
//   plots/trig_eff_mb_<mu|ele>/         turnon_pt_<sel>[_zoom25], eff_y_<sel>,
//                                       eff_y_<sel>_charge, eff_pt_<sel>_charge,
//                                       eff_absy_<sel>[_charge], mbfrac_pt_<sel>,
//                                       bit_vs_match_pt_<sel>,
//                                       trig_eff_<sel>.csv (per-bin table)
//
// WHAT IS APPLIED (2026-09-21, user decision): ONE INCLUSIVE SF -- skim/muon_sf.h
// runs with kTrigBinning = kTrigInclusive, i.e. 0.9971 (+0.0020 -0.0024) for
// every muon. This macro nevertheless measures, prints and plots the full
// rapidity dependence: that record is what JUSTIFIES the choice, and it is the
// drop-in input if the decision is revisited (flip that one constant back to
// kTrigPerAbsY and the sf_absy_<sel> + h_{den,num,bit}_absy_<sel>_* histograms
// below are applied instead, unchanged).
//
// IF a binned SF is applied, the binning is |y| = |eta_lab|, 6 bins of 0.4,
// charge-inclusive -- established by the likelihood-ratio tests this macro
// prints under "BINNING DECISION" (2026-09-15, mt40 / nom):
//   flat -> pT (2 groups, </>35)  p = 0.71  / 0.80   => NO pT dependence
//   flat -> |y| 6 bins            p = 0.0006/ 0.0004 => the eta dependence is real
//   |y| 3 -> |y| 6 bins           p = 0.0043/ 0.0132 => 3 coarse bins are not enough
//   |y| 6 -> signed y 12          p = 0.39  / 0.14   => folding to |y| is justified
//   |y|3 -> |y|3 x pT2            p = 0.68  / 0.68   => NO pT x eta interaction, 1D
//   |y|6 -> |y|6 x charge         p = 0.13           => charge-inclusive
// The per-bin pT scan (11 bins) gives p = 0.087 / 0.0054, but that answers a
// different question -- "is any single pT bin anomalous?" -- and is driven
// entirely by the [50,60) bin (SF 0.964 / 0.962 on ~150 events, the same ~2
// sigma downward fluctuation in both selections, no trend). The decision test
// is the coarse one above.
// The SF runs 0.982 (|y| < 0.4) to 1.006 (|y| > 1.2), a 2.4% spread; the MC
// itself dips at |y| < 0.4 (eps_MC 0.985 vs 0.994) = the eta ~ 0 barrel wheel
// gap, and the data dips further (0.967), i.e. the L1 emulation under-models a
// real, localized detector feature. So the structure is real and is NOT being
// corrected for: leaving it out is a deliberate choice, whose cost is a pure
// rapidity SHAPE effect on the W templates of -1.5% (|y| < 0.4) to +1.0%
// (|y| > 2.0), |y|- and charge-symmetric, with the inclusive normalization
// (and hence sigma_incl and the charge asymmetry) unchanged by construction --
// only dsigma/deta and R_FB would move. Do not quote the flatness chi2 of the
// [SF] log block as the justification; see the note there.
//   stdout                              the tables -- run through
//                                       ./run_trig_eff_mb.sh to keep the log
//
// Run:
//   ./run_trig_eff_mb.sh [mu|ele|both]           (keeps logs/trig_eff_mb_<chan>.log)
//   root -l -b -q 'trig_eff_mb.C+("mu")'          (bare; ~5 min data + ~1 min per MC file)
//   root -l -b -q 'trig_eff_mb.C+("both", true)'  (tables + plots only, from the rootfile)

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "../skim/skim_common.h"
#include "../skim/mc_norm.h"
#include "../plotting/plotting_helper.C"

using namespace pOSkim;

namespace
{

// -------- selection constants (tracking skim.C) --------
const double kPtFloor     = 10.0;  // loop floor for steps 1 and 6 (analysis = 25)
const double kPtNominal   = 25.0;  // the analysis cut; MUST be a bin edge of kPtEdges
const double kEtaMax      = 2.4;
const double kVzMax       = 15.0;
const double kTrigMatchDR = 0.4;   // skim.C:284/1008
const double kMtCut       = 40.0;  // the leppt_mt40 discriminant selection
const double kDyMassMin   = 80.0, kDyMassMax = 110.0;

// -------- trigger paths --------
const char *kMBPath = "HLT_MinimumBiasHF_OR_BptxAND_v1";

// -------- binning --------
const int    kNPt = 18;
const double kPtEdges[kNPt + 1] = {10, 12, 14, 16, 18, 20, 22, 25, 27.5, 30, 32.5, 35, 37.5, 40, 45, 50, 60, 80, 120};
const double kPtFoldMax = 119.9; // pT above the last edge is folded into the last bin

const char *kSel[2]      = {"nom", "mt40"};
const char *kSelLabel[2] = {"W selection (no m_{T} cut)", "W selection, m_{T} > 40 GeV"};

// |y| = |eta_lab| binning -- the folded twin of kYEdges. This is the binning a
// rapidity-dependent trigger SF WOULD use (folding is justified by the
// likelihood-ratio test printed below -- |y| 6 -> signed y 12 is not
// significant -- and doubles the data statistics per bin). skim/muon_sf.h
// currently applies the INCLUSIVE SF instead (kTrigBinning = kTrigInclusive,
// user decision 2026-09-21), but still reads this table and prints it as the
// measured-but-not-applied record; switching back needs only that constant.
// The signed 12-bin table is kept as the left/right asymmetry cross-check.
const int    kNAbsY = 6;
const double kAbsYEdges[kNAbsY + 1] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};

// Gen matching of the MC leg (2026-09-15, user decision). Charge-blind window,
// identical to correction/charge_flip.C's kMatchDR / kMatchDPt.
const double kGenMatchDR  = 0.5;
const double kGenMatchDPt = 0.5;

const double kCL = 0.6827;

// ============================================================
// Gen match of the MC leg: is the selected reco lepton a prompt W lepton?
//
// Reads the EventTree's OWN gen block (nMC/mcPID/mcStatus/mcPt/mcEta/mcPhi/
// mcMomPID/mcGMomPID) -- NOT HiGenParticleAna/hi, whose motherIdx is -999 for
// every W lepton so that a W-ancestor test there is always false (measured in
// correction/charge_flip.C). Verified on July_29_MC_Wp_mu: one |pdg| = 24 entry
// per event, and 99.98% of gen muons with pT > 20, |eta| < 2.4 carry
// mcMomPID = +-24 directly, all with mcStatus = 1.
//
// CHARGE-BLIND (the sign of mcPID is ignored) so that a charge-misidentified
// lepton still counts as a prompt W lepton -- charge misID is not a trigger
// inefficiency, and folding it in here would double-count correction/charge_flip.C.
// FSR chains (l <- l <- W) are accepted via mcGMomPID.
// Returns the gen index, or -1 if no prompt-W lepton matches.
// ============================================================
int MatchGenLeptonFromW(double pt, double eta, double phi, int flavPdg,
                        const std::vector<int> *gPid, const std::vector<int> *gStatus,
                        const std::vector<float> *gPt, const std::vector<float> *gEta,
                        const std::vector<float> *gPhi, const std::vector<int> *gMom,
                        const std::vector<int> *gGMom)
{
  if (!gPid || !gStatus || !gPt || !gEta || !gPhi || !gMom) return -1;
  const size_t ng = std::min({gPid->size(), gStatus->size(), gPt->size(), gEta->size(),
                              gPhi->size(), gMom->size()});
  int    best   = -1;
  double drBest = 1e9;
  for (size_t i = 0; i < ng; ++i)
  {
    if (std::abs(gPid->at(i)) != flavPdg) continue;
    if (gStatus->at(i) != 1) continue;
    if (gPt->at(i) <= 0) continue;
    const bool momW  = std::abs(gMom->at(i)) == 24;
    const bool fsrW  = (std::abs(gMom->at(i)) == flavPdg) && gGMom && i < gGMom->size() &&
                       std::abs(gGMom->at(i)) == 24;
    if (!momW && !fsrW) continue;
    const double d = DeltaR(gEta->at(i), gPhi->at(i), eta, phi);
    if (d >= kGenMatchDR) continue;
    if (std::fabs(gPt->at(i) - pt) / gPt->at(i) >= kGenMatchDPt) continue;
    if (d < drBest) { drBest = d; best = (int)i; }
  }
  return best;
}

// ============================================================
// Is the SF flat across a binning?  LIKELIHOOD-RATIO test.
//
// Model A: one common SF, eps_data,i = SF * eps_MC,i.  Model B: a free SF per
// bin.  Data counts are binomial, k_i ~ B(n_i, SF * eps_MC,i); eps_MC is taken
// as exact (MC denominators are ~200x the data ones, so its error is ~0.03%).
// -2 dlnL is then chi2-distributed with (nbins - 1) dof.
//
// This REPLACES the Clopper-Pearson pull chi2 as the decision test. Toys under
// a true flat null with the real denominators (2026-09-15) show the CP pull
// chi2 gives <chi2> = 3.15 for 5 dof and rejects at only 0.3% when it should
// reject at 5% -- it is ~16x under-powered, because CP intervals over-cover
// badly as eps -> 1 (and 2 of the 6 |y| bins have eps_data = 1 exactly). The
// LRT is calibrated on the same toys: <-2dlnL> = 5.48 vs 5.00, 7.2% vs 5%.
// The old p = 0.41 "consistent with a flat SF" was an artifact of that test.
// ============================================================
struct FlatFit { double sf = 1; double lrt = 0; int ndf = 0; };

double BinomLnL(double k, double n, double p)
{
  if (p <= 0) return (k > 0) ? -1e30 : 0.0;
  if (p >= 1) return (k < n) ? -1e30 : 0.0;
  return k * std::log(p) + (n - k) * std::log(1.0 - p);
}

// One cell of the test: data den/num and MC den/num.
struct LrtCell { double nd = 0, kd = 0, nm = 0, km = 0; };

// Cells with an empty data or MC denominator are skipped.
FlatFit FlatnessLRT(const std::vector<LrtCell> &c)
{
  FlatFit r;
  double best = -1e30;
  for (int g = 0; g <= 20000; ++g) // scan the common SF
  {
    const double s = 0.80 + g * (0.40 / 20000.0);
    double L = 0; bool ok = true;
    for (const LrtCell &x : c)
    {
      if (x.nd <= 0 || x.nm <= 0) continue;
      const double p = s * (x.km / x.nm);
      if (p > 1) { ok = false; break; }
      L += BinomLnL(x.kd, x.nd, p);
    }
    if (ok && L > best) { best = L; r.sf = s; }
  }
  double lFree = 0;
  r.ndf = -1;
  for (const LrtCell &x : c)
  {
    if (x.nd <= 0 || x.nm <= 0) continue;
    lFree += BinomLnL(x.kd, x.nd, std::min(x.kd / x.nd, 1.0));
    ++r.ndf;
  }
  r.lrt = 2.0 * (lFree - best);
  if (r.ndf < 0) r.ndf = 0;
  return r;
}

FlatFit FlatnessLRT(const TH1 *dDen, const TH1 *dNum, const TH1 *mDen, const TH1 *mNum)
{
  if (!dDen || !dNum || !mDen || !mNum) return FlatFit();
  std::vector<LrtCell> c;
  for (int i = 1; i <= dDen->GetNbinsX(); ++i)
    c.push_back({dDen->GetBinContent(i), dNum->GetBinContent(i),
                 mDen->GetBinContent(i), mNum->GetBinContent(i)});
  return FlatnessLRT(c);
}

// Coarse (pT group) x (|y| group) cells from the 2D maps -- the interaction test.
// nPtGrp = 1 pools pT; nAbsGrp = 1 pools |y|. ptSplit is the pT bin index (1-based,
// inclusive) ending the first pT group; ptLo is the first bin of the analysis range.
std::vector<LrtCell> Cells2D(const TH2 *dD, const TH2 *dN, const TH2 *mD, const TH2 *mN,
                             int ptLo, int ptSplit, int nPtGrp, int nAbsGrp)
{
  const int nx = dD->GetNbinsX(), ny = dD->GetNbinsY();
  std::vector<LrtCell> c(nPtGrp * nAbsGrp);
  for (int ix = ptLo; ix <= nx; ++ix)
  {
    const int gp = (nPtGrp == 1) ? 0 : (ix <= ptSplit ? 0 : 1);
    for (int iy = 1; iy <= ny; ++iy)
    {
      const int ay = (iy <= ny / 2) ? (ny / 2 - iy) : (iy - ny / 2 - 1); // 0..5 = |y| index
      const int ga = (nAbsGrp == 1) ? 0 : std::min(ay * nAbsGrp / (ny / 2), nAbsGrp - 1);
      LrtCell &x = c[gp * nAbsGrp + ga];
      x.nd += dD->GetBinContent(ix, iy); x.kd += dN->GetBinContent(ix, iy);
      x.nm += mD->GetBinContent(ix, iy); x.km += mN->GetBinContent(ix, iy);
    }
  }
  return c;
}

void PrintFlatness(const char *label, const TH1 *dDen, const TH1 *dNum, const TH1 *mDen, const TH1 *mNum)
{
  const FlatFit r = FlatnessLRT(dDen, dNum, mDen, mNum);
  std::cout << Form("  %-30s common SF %.4f   -2dlnL = %6.2f / %2d dof   p = %.4f  %s\n",
                    label, r.sf, r.lrt, r.ndf, TMath::Prob(r.lrt, r.ndf),
                    TMath::Prob(r.lrt, r.ndf) < 0.05 ? "<-- NOT flat" : "");
}

// ============================================================
// Per-sample histogram set
// ============================================================
struct TrigOut
{
  std::string tag;
  std::map<std::string, TH1D *> h1;
  std::map<std::string, TH2D *> h2;
  double nSel[2] = {0, 0}, nMB[2] = {0, 0}, nTrg[2] = {0, 0}, nNum[2] = {0, 0}; // pT > 25 counts
  std::map<int, std::vector<double>> byRun;                                    // run -> {den, num} (nom, pT>25)

  TH1D *H1(const std::string &n) const
  {
    auto it = h1.find(n);
    if (it == h1.end()) { std::cerr << "[FATAL] TrigOut: no histo " << n << "\n"; std::abort(); }
    return it->second;
  }
  TH2D *H2(const std::string &n) const
  {
    auto it = h2.find(n);
    if (it == h2.end()) { std::cerr << "[FATAL] TrigOut: no histo2 " << n << "\n"; std::abort(); }
    return it->second;
  }

  void Book(const std::string &t)
  {
    tag = t;
    auto b1 = [&](const std::string &n, int nb, const double *e, const char *xt)
    {
      TH1D *h = new TH1D(Form("h_%s_%s", n.c_str(), tag.c_str()), Form(";%s;events", xt), nb, e);
      h->Sumw2();
      h->SetDirectory(nullptr);
      h1[n] = h;
    };
    auto b2 = [&](const std::string &n)
    {
      TH2D *h = new TH2D(Form("h_%s_%s", n.c_str(), tag.c_str()), ";p_{T} [GeV];y = -#eta_{lab}",
                         kNPt, kPtEdges, kNY, kYEdges);
      h->Sumw2();
      h->SetDirectory(nullptr);
      h2[n] = h;
    };
    for (int s = 0; s < 2; ++s)
    {
      const std::string S = kSel[s];
      // stages: all (analysis cuts only), den (+MB), bit (+MB+HLT bit), num (+MB+HLT+match),
      //         trg (analysis cuts + HLT + match, NO MB) -- for the MB-fraction control
      for (const char *st : {"all", "den", "bit", "num", "trg"})
      {
        b1(std::string(st) + "_pt_" + S, kNPt, kPtEdges, "leading-lepton p_{T} [GeV]");
        b1(std::string(st) + "_pt_" + S + "_plus", kNPt, kPtEdges, "leading-lepton p_{T} [GeV]");
        b1(std::string(st) + "_pt_" + S + "_minus", kNPt, kPtEdges, "leading-lepton p_{T} [GeV]");
        b1(std::string(st) + "_y_" + S, kNY, kYEdges, "y = -#eta_{lab}");          // pT > 25 only
        b1(std::string(st) + "_y_" + S + "_plus", kNY, kYEdges, "y = -#eta_{lab}");
        b1(std::string(st) + "_y_" + S + "_minus", kNY, kYEdges, "y = -#eta_{lab}");
        // |y| = |eta_lab|, folded -- the binning a rapidity-dependent trigger
        // SF would use; read by skim/muon_sf.h either way (see the header)
        b1(std::string(st) + "_absy_" + S, kNAbsY, kAbsYEdges, "|y| = |#eta_{lab}|");
        b1(std::string(st) + "_absy_" + S + "_plus", kNAbsY, kAbsYEdges, "|y| = |#eta_{lab}|");
        b1(std::string(st) + "_absy_" + S + "_minus", kNAbsY, kAbsYEdges, "|y| = |#eta_{lab}|");
      }
      // gen-weighted den/num (MC cross-check of the raw-count efficiency)
      b1("den_pt_" + S + "_w", kNPt, kPtEdges, "leading-lepton p_{T} [GeV]");
      b1("num_pt_" + S + "_w", kNPt, kPtEdges, "leading-lepton p_{T} [GeV]");
      b2("den_2d_" + S);
      b2("num_2d_" + S);
    }
  }

  // Fill one selected event. stage flags: mb, bit (HLT path fired), match.
  void Fill(int s, double pt, double y, int q, bool mb, bool bit, bool match, double w)
  {
    const std::string S = kSel[s];
    const double ptF = std::min(pt, kPtFoldMax);
    const bool   nomPt = (pt > kPtNominal);
    const char  *qs = (q > 0) ? "_plus" : "_minus";
    auto fill = [&](const char *st)
    {
      H1(std::string(st) + "_pt_" + S)->Fill(ptF);
      H1(std::string(st) + "_pt_" + S + qs)->Fill(ptF);
      if (nomPt)
      {
        H1(std::string(st) + "_y_" + S)->Fill(y);
        H1(std::string(st) + "_y_" + S + qs)->Fill(y);
        H1(std::string(st) + "_absy_" + S)->Fill(std::fabs(y));
        H1(std::string(st) + "_absy_" + S + qs)->Fill(std::fabs(y));
      }
    };
    fill("all");
    if (bit && match) fill("trg");
    if (mb)
    {
      fill("den");
      H1("den_pt_" + S + "_w")->Fill(ptF, w);
      H2("den_2d_" + S)->Fill(ptF, y);
      if (bit) fill("bit");
      if (bit && match)
      {
        fill("num");
        H1("num_pt_" + S + "_w")->Fill(ptF, w);
        H2("num_2d_" + S)->Fill(ptF, y);
      }
    }
    if (nomPt)
    {
      nSel[s] += 1;
      if (bit && match) nTrg[s] += 1;
      if (mb) { nMB[s] += 1; if (bit && match) nNum[s] += 1; }
    }
  }

  void Write(TFile *f) const
  {
    f->cd();
    for (auto &p : h1) p.second->Write("", TObject::kOverwrite);
    for (auto &p : h2) p.second->Write("", TObject::kOverwrite);
  }
};

// ============================================================
// The event loop: W selection minus the trigger steps, on one file
// ============================================================
int RunTrig(bool isMu, SampleType sample, TrigOut &out)
{
  const bool   isMC      = IsMC(sample);
  const double dyPtMin   = isMu ? 15.0 : 10.0;   // DY-veto leg pT (intentionally asymmetric)
  const double isoMax    = isMu ? 0.15 : 0.095;
  const double lepMass   = isMu ? MU_MASS : ELE_MASS;
  const char  *tagLep    = isMu ? "muon" : "electron";

  std::string fname;
  if (isMC)
  {
    const auto info = ResolveMCSample(sample, isMu ? "mu" : "ele");
    if (info.fname.empty()) { std::cerr << "[ERR] RunTrig: cannot resolve MC sample\n"; return 1; }
    fname = info.fname;
  }
  else
    fname = kDefaultDataFile;

  std::cout << "[INPUT] " << out.tag << " <- " << fname << std::endl;
  TFile *f = TFile::Open(fname.c_str());
  if (!f || f->IsZombie()) { std::cerr << "[ERR] cannot open " << fname << "\n"; return 1; }

  TTree *tLep    = (TTree *)f->Get("ggHiNtuplizer/EventTree");
  TTree *tHi     = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  TTree *tPF     = (TTree *)f->Get("particleFlowAnalyser/pftree");
  TTree *tHLT    = (TTree *)f->Get("hltanalysis/HltTree");
  TTree *tHLTobj = (TTree *)f->Get(isMu ? "hltobject/HLT_OxyL1SingleMuOpen_v"
                                        : "hltobject/HLT_OxyL1SingleEG10_v");
  TTree *tEvent  = (TTree *)f->Get("skimanalysis/HltTree");
  if (!tLep || !tHi || !tPF || !tHLT || !tHLTobj || !tEvent)
  {
    std::cerr << "[FATAL] RunTrig: missing a required tree in " << fname << "\n";
    f->Close();
    return 2;
  }
  for (TTree *t : {tHi, tPF, tHLT, tHLTobj, tEvent})
    if (t->GetEntries() != tLep->GetEntries())
    {
      std::cerr << "[FATAL] RunTrig: tree " << t->GetName() << " (" << t->GetEntries()
                << ") and EventTree (" << tLep->GetEntries() << ") entry counts differ\n";
      f->Close();
      return 2;
    }

  // -------- lepton branches --------
  const char *bN   = isMu ? "nMu"        : "nEle";
  const char *bPt  = isMu ? "muPt"       : "elePt";
  const char *bEta = isMu ? "muEta"      : "eleEta";
  const char *bPhi = isMu ? "muPhi"      : "elePhi";
  const char *bChg = isMu ? "muCharge"   : "eleCharge";
  const char *bID  = isMu ? "muIDTight"  : "eleMVAIdWP95";
  const char *bCh  = isMu ? "muPFChIso"  : "elePFChIso";
  const char *bNeu = isMu ? "muPFNeuIso" : "elePFNeuIso";
  const char *bPho = isMu ? "muPFPhoIso" : "elePFPhoIso";
  const char *bPU  = isMu ? "muPFPUIso"  : "elePFPUIso";

  Int_t nLep = 0;
  std::vector<float> *lepPt = nullptr, *lepEta = nullptr, *lepPhi = nullptr;
  std::vector<int>   *lepCharge = nullptr, *lepID = nullptr, *muIsPF = nullptr;
  std::vector<float> *chIso = nullptr, *neuIso = nullptr, *phoIso = nullptr, *puIso = nullptr;

  tLep->SetBranchStatus("*", 0);
  for (const char *bn : {bN, bPt, bEta, bPhi, bChg, bID, bCh, bNeu, bPho})
    if (!HasBranch(tLep, bn))
    {
      std::cerr << "[FATAL] RunTrig: missing EventTree branch " << bn << "\n";
      f->Close();
      return 2;
    }
  tLep->SetBranchStatus(bN, 1);   tLep->SetBranchAddress(bN,   &nLep);
  tLep->SetBranchStatus(bPt, 1);  tLep->SetBranchAddress(bPt,  &lepPt);
  tLep->SetBranchStatus(bEta, 1); tLep->SetBranchAddress(bEta, &lepEta);
  tLep->SetBranchStatus(bPhi, 1); tLep->SetBranchAddress(bPhi, &lepPhi);
  tLep->SetBranchStatus(bChg, 1); tLep->SetBranchAddress(bChg, &lepCharge);
  tLep->SetBranchStatus(bID, 1);  tLep->SetBranchAddress(bID,  &lepID);
  tLep->SetBranchStatus(bCh, 1);  tLep->SetBranchAddress(bCh,  &chIso);
  tLep->SetBranchStatus(bNeu, 1); tLep->SetBranchAddress(bNeu, &neuIso);
  tLep->SetBranchStatus(bPho, 1); tLep->SetBranchAddress(bPho, &phoIso);

  // ECAL-gap veto input (electrons only, 2026-09-14) -- tracks the skim's eleIDNoGap
  std::vector<float> *scEta = nullptr;
  if (!isMu)
  {
    if (!HasBranch(tLep, "eleSCEta")) { std::cerr << "[FATAL] RunTrig: missing EventTree branch eleSCEta (ECAL-gap veto)\n"; f->Close(); return 2; }
    tLep->SetBranchStatus("eleSCEta", 1); tLep->SetBranchAddress("eleSCEta", &scEta);
  }

  const bool has_puIso = HasBranch(tLep, bPU);
  if (has_puIso) { tLep->SetBranchStatus(bPU, 1); tLep->SetBranchAddress(bPU, &puIso); }
  else std::cout << "[WARN] RunTrig: no " << bPU << " branch; relIso uncorrected (no Delta-beta).\n";

  const bool has_muIsPF = isMu && HasBranch(tLep, "muIsPF");
  if (has_muIsPF) { tLep->SetBranchStatus("muIsPF", 1); tLep->SetBranchAddress("muIsPF", &muIsPF); }

  // -------- gen block for the MC leg's prompt-W requirement (2026-09-15) --------
  // MANDATORY in MC: a silent fall-back to "no gen match required" would change
  // the definition of eps_MC without any message.
  std::vector<int>   *mcPID = nullptr, *mcStatus = nullptr, *mcMomPID = nullptr, *mcGMomPID = nullptr;
  std::vector<float> *mcPt = nullptr, *mcEta = nullptr, *mcPhi = nullptr;
  if (isMC)
  {
    for (const char *bn : {"mcPID", "mcStatus", "mcPt", "mcEta", "mcPhi", "mcMomPID"})
      if (!HasBranch(tLep, bn))
      {
        std::cerr << "[FATAL] RunTrig: MC file lacks EventTree branch " << bn
                  << " -- the prompt-W gen match of the MC leg cannot be applied.\n";
        f->Close();
        return 2;
      }
    tLep->SetBranchStatus("mcPID", 1);    tLep->SetBranchAddress("mcPID",    &mcPID);
    tLep->SetBranchStatus("mcStatus", 1); tLep->SetBranchAddress("mcStatus", &mcStatus);
    tLep->SetBranchStatus("mcPt", 1);     tLep->SetBranchAddress("mcPt",     &mcPt);
    tLep->SetBranchStatus("mcEta", 1);    tLep->SetBranchAddress("mcEta",    &mcEta);
    tLep->SetBranchStatus("mcPhi", 1);    tLep->SetBranchAddress("mcPhi",    &mcPhi);
    tLep->SetBranchStatus("mcMomPID", 1); tLep->SetBranchAddress("mcMomPID", &mcMomPID);
    if (HasBranch(tLep, "mcGMomPID")) { tLep->SetBranchStatus("mcGMomPID", 1); tLep->SetBranchAddress("mcGMomPID", &mcGMomPID); }
    else std::cout << "[WARN] RunTrig: no mcGMomPID; FSR chains (l <- l <- W) will not be accepted.\n";
    std::cout << "[CONFIG] MC leg: leading lepton required to gen-match a prompt W lepton"
              << " (|pdg| = " << (isMu ? 13 : 11) << ", status 1, |mcMomPID| = 24 or FSR via mcGMomPID,"
              << " charge-blind, DR < " << kGenMatchDR << ", |dpT|/pT < " << kGenMatchDPt << ")\n";
  }

  // -------- HLT bits: the analysis path AND the MB path (both mandatory here) --------
  const std::string hltNeedle = isMu ? "HLT_OxyL1SingleMuOpen_v1" : "HLT_OxyL1SingleEG10_v1";
  const std::string hltName   = FindBranchContaining(tHLT, hltNeedle);
  const std::string mbName    = FindBranchContaining(tHLT, kMBPath);
  Int_t hltBit = 0, mbBit = 0;
  tHLT->SetBranchStatus("*", 0);
  if (hltName.empty() || mbName.empty())
  {
    std::cerr << "[FATAL] RunTrig: hltanalysis/HltTree lacks " << (hltName.empty() ? hltNeedle : std::string(kMBPath))
              << " -- the MB-denominator efficiency cannot be defined on this file.\n";
    f->Close();
    return 2;
  }
  // FindBranchContaining may return the _Prescale* twin; insist on the exact bit name.
  auto exact = [&](const std::string &needle) -> std::string
  {
    return HasBranch(tHLT, needle.c_str()) ? needle : FindBranchContaining(tHLT, needle);
  };
  const std::string hltExact = exact(hltNeedle), mbExact = exact(kMBPath);
  tHLT->SetBranchStatus(hltExact.c_str(), 1); tHLT->SetBranchAddress(hltExact.c_str(), &hltBit);
  tHLT->SetBranchStatus(mbExact.c_str(), 1);  tHLT->SetBranchAddress(mbExact.c_str(),  &mbBit);
  Int_t run = 0;
  const bool has_run = HasBranch(tHLT, "Run");
  if (has_run) { tHLT->SetBranchStatus("Run", 1); tHLT->SetBranchAddress("Run", &run); }
  std::cout << "[CONFIG] analysis path: " << hltExact << "   MB path: " << mbExact
            << "   match DR < " << kTrigMatchDR << "   pT floor " << kPtFloor << " (analysis " << kPtNominal << ")\n";

  // -------- HLT objects (trigger match) --------
  std::vector<double> *toPt = nullptr, *toEta = nullptr, *toPhi = nullptr;
  const bool has_toPt  = HasBranch(tHLTobj, "pt");
  const bool has_toEta = HasBranch(tHLTobj, "eta");
  const bool has_toPhi = HasBranch(tHLTobj, "phi");
  tHLTobj->SetBranchStatus("*", 0);
  if (has_toPt && has_toEta && has_toPhi)
  {
    tHLTobj->SetBranchStatus("pt", 1);  tHLTobj->SetBranchAddress("pt",  &toPt);
    tHLTobj->SetBranchStatus("eta", 1); tHLTobj->SetBranchAddress("eta", &toEta);
    tHLTobj->SetBranchStatus("phi", 1); tHLTobj->SetBranchAddress("phi", &toPhi);
  }
  else
  {
    std::cerr << "[FATAL] RunTrig: no pt/eta/phi in " << tHLTobj->GetName() << " -- no trigger objects to match.\n";
    f->Close();
    return 2;
  }

  // -------- event filters + vz + gen weight --------
  Int_t ppv = 1, pcc = 1;
  const bool has_ppv = HasBranch(tEvent, "pprimaryVertexFilter");
  const bool has_pcc = HasBranch(tEvent, "pclusterCompatibilityFilter");
  tEvent->SetBranchStatus("*", 0);
  if (has_ppv) { tEvent->SetBranchStatus("pprimaryVertexFilter", 1);        tEvent->SetBranchAddress("pprimaryVertexFilter", &ppv); }
  if (has_pcc) { tEvent->SetBranchStatus("pclusterCompatibilityFilter", 1); tEvent->SetBranchAddress("pclusterCompatibilityFilter", &pcc); }

  Float_t vz = 999.f;
  tHi->SetBranchStatus("*", 0);
  const bool has_vz = HasBranch(tHi, "vz");
  if (has_vz) { tHi->SetBranchStatus("vz", 1); tHi->SetBranchAddress("vz", &vz); }

  Float_t    genWeight     = 1.f;
  const bool has_genWeight = isMC && HasBranch(tHi, "weight");
  if (has_genWeight) { tHi->SetBranchStatus("weight", 1); tHi->SetBranchAddress("weight", &genWeight); }
  else if (isMC) std::cout << "[WARN] RunTrig: no HiTree 'weight'; weighted sums unweighted.\n";

  // -------- PF tree (MET for the mt40 variant; read only for selected events) --------
  std::vector<float> *pfPt = nullptr, *pfPhi = nullptr;
  tPF->SetBranchStatus("*", 0);
  if (!HasBranch(tPF, "pfPt") || !HasBranch(tPF, "pfPhi"))
  {
    std::cerr << "[FATAL] RunTrig: pftree lacks pfPt/pfPhi\n";
    f->Close();
    return 2;
  }
  tPF->SetBranchStatus("pfPt", 1);  tPF->SetBranchAddress("pfPt",  &pfPt);
  tPF->SetBranchStatus("pfPhi", 1); tPF->SetBranchAddress("pfPhi", &pfPhi);

  // -------- event loop --------
  const Long64_t nEntries = tLep->GetEntries();
  std::cout << "Entries: " << nEntries << "\n";
  bool warnedFilters = false, warnedTrig = false;
  Long64_t nPass = 0;
  Long64_t nGenTried = 0, nGenFail = 0, nGenFail25 = 0; // prompt-W gen match bookkeeping (MC)

  for (Long64_t ie = 0; ie < nEntries; ++ie)
  {
    if (ie % 500000 == 0) std::cout << "  event " << ie << "/" << nEntries << "\n";

    tLep->GetEntry(ie);
    tHi->GetEntry(ie);
    tHLT->GetEntry(ie);
    tHLTobj->GetEntry(ie);
    tEvent->GetEntry(ie);

    const double w = has_genWeight ? (double)genWeight : 1.0;

    if (!lepPt || !lepEta || !lepPhi || !lepCharge || !lepID) continue;

    // (1) exists PF lepton pT > floor (analysis: 25; the leading tight lepton
    //     with pT > 25 below implies the analysis version of this step)
    bool hasPF = false;
    for (int i = 0; i < nLep; ++i)
    {
      if (lepPt->at(i) <= kPtFloor) continue;
      if (has_muIsPF && muIsPF && muIsPF->at(i) == 0) continue;
      hasPF = true;
      break;
    }
    if (!hasPF) continue;

    // (2) pO event filters + vz
    if (!PassEventSelection_pO(warnedFilters, has_ppv, ppv, has_pcc, pcc)) continue;
    if (has_vz && TMath::Abs(vz) > kVzMax) continue;

    // (3) trigger fired  -- NOT applied: this is the numerator condition

    // (4) DY veto: OS pair, both legs ID'd + isolated, mll in (80, 110)
    {
      std::vector<int> cand;
      cand.reserve(nLep);
      for (int i = 0; i < nLep; ++i)
      {
        if (lepPt->at(i) <= dyPtMin) continue;
        if (lepID->at(i) == 0) continue;
        if (scEta && pOSkim::InEcalGap(scEta->at(i))) continue; // e: ECAL crack veto
        if (has_muIsPF && muIsPF && muIsPF->at(i) == 0) continue;
        if (RelIsoPF(i, lepPt, chIso, neuIso, phoIso, puIso) >= isoMax) continue;
        cand.push_back(i);
      }
      bool veto = false;
      for (size_t a = 0; a < cand.size() && !veto; ++a)
        for (size_t b = a + 1; b < cand.size() && !veto; ++b)
        {
          const int i1 = cand[a], i2 = cand[b];
          if (lepCharge->at(i1) * lepCharge->at(i2) >= 0) continue;
          TLorentzVector p1, p2;
          p1.SetPtEtaPhiM(lepPt->at(i1), lepEta->at(i1), lepPhi->at(i1), lepMass);
          p2.SetPtEtaPhiM(lepPt->at(i2), lepEta->at(i2), lepPhi->at(i2), lepMass);
          const double mll = (p1 + p2).M();
          if (mll > kDyMassMin && mll < kDyMassMax) veto = true;
        }
      if (veto) continue;
    }

    // (5)+(6) leading ID'd (mu: +PF) lepton, pT > floor, |eta| < 2.4
    int iLead = -1;
    double bestPt = -1.0;
    for (int i = 0; i < nLep; ++i)
    {
      if (lepID->at(i) == 0) continue;
      if (scEta && pOSkim::InEcalGap(scEta->at(i))) continue; // e: ECAL crack veto
      if (has_muIsPF && muIsPF && muIsPF->at(i) == 0) continue;
      if (lepPt->at(i) > bestPt) { bestPt = lepPt->at(i); iLead = i; }
    }
    if (iLead < 0) continue;
    if (lepPt->at(iLead) <= kPtFloor) continue;
    if (std::abs(lepEta->at(iLead)) > kEtaMax) continue;

    // (7) leading-lepton isolation
    if (RelIsoPF(iLead, lepPt, chIso, neuIso, phoIso, puIso) >= isoMax) continue;

    // (8) trigger match -- NOT applied: numerator condition, evaluated below

    // (MC only) prompt-W gen match of the selected lepton. Applied BEFORE any
    // den/num fill, so it moves eps_MC's denominator and numerator together and
    // makes the MC leg "the efficiency for a real W lepton" by construction
    // rather than by sample purity. Charge-blind; see MatchGenLeptonFromW.
    if (isMC)
    {
      ++nGenTried;
      if (MatchGenLeptonFromW(lepPt->at(iLead), lepEta->at(iLead), lepPhi->at(iLead),
                              isMu ? 13 : 11, mcPID, mcStatus, mcPt, mcEta, mcPhi,
                              mcMomPID, mcGMomPID) < 0)
      {
        ++nGenFail;
        if (lepPt->at(iLead) > kPtNominal) ++nGenFail25;
        continue;
      }
    }
    ++nPass;

    const bool mb    = TriggerFired(mbBit);
    const bool bit   = TriggerFired(hltBit);
    const bool match = bit && PassLeadingLeptonTrigMatch(kTrigMatchDR, iLead, lepEta, lepPhi,
                                                         has_toPt, toPt, has_toEta, toEta,
                                                         has_toPhi, toPhi, warnedTrig, tagLep);

    const double pt  = lepPt->at(iLead);
    const double y   = -lepEta->at(iLead);          // skim convention: p-going (-Z) = forward
    const int    q   = lepCharge->at(iLead);

    out.Fill(0, pt, y, q, mb, bit, match, w);

    // mt40 twin: PF MET only now (selected events are a tiny fraction)
    tPF->GetEntry(ie);
    const TVector2 metv = ComputePFMET(nullptr, pfPt, pfPhi);
    const double   mt   = TransverseMass(pt, lepPhi->at(iLead), metv);
    if (mt > kMtCut) out.Fill(1, pt, y, q, mb, bit, match, w);

    if (!isMC && has_run && pt > kPtNominal)
    {
      auto &v = out.byRun[run];
      if (v.empty()) v.assign(4, 0.0);
      v[0] += 1;                    // selected
      if (mb) v[1] += 1;            // den
      if (mb && match) v[2] += 1;   // num
      if (match) v[3] += 1;         // triggered (no MB)
    }
  }

  f->Close();
  delete f;
  std::cout << Form("[INFO] %s: %lld events pass the trigger-free W selection (pT > %.0f); pT > %.0f: sel %.0f, MB %.0f, num %.0f  |  mt40: sel %.0f, MB %.0f, num %.0f\n",
                    out.tag.c_str(), nPass, kPtFloor, kPtNominal,
                    out.nSel[0], out.nMB[0], out.nNum[0], out.nSel[1], out.nMB[1], out.nNum[1]);
  if (isMC)
    std::cout << Form("[GENMATCH] %s: %lld selected, %lld fail the prompt-W match (%.3f%%; pT > %.0f: %lld)\n",
                      out.tag.c_str(), nGenTried, nGenFail,
                      nGenTried > 0 ? 100.0 * nGenFail / nGenTried : 0.0, kPtNominal, nGenFail25);
  return 0;
}

// ============================================================
// Efficiency / ratio helpers
// ============================================================
TEfficiency *MakeEff(TH1 *pass, TH1 *tot, const std::string &name, const char *title = "")
{
  if (!pass || !tot) return nullptr;
  if (!TEfficiency::CheckConsistency(*pass, *tot))
  {
    std::cerr << "[ERR] TEfficiency inconsistency for " << name << "\n";
    return nullptr;
  }
  TEfficiency *e = new TEfficiency(*pass, *tot);
  e->SetName(name.c_str());
  e->SetTitle(title);
  e->SetStatisticOption(TEfficiency::kFCP);
  e->SetConfidenceLevel(kCL);
  e->SetDirectory(nullptr);
  return e;
}

// a / b per bin with the CP intervals propagated (a_up with b_dn and vice versa).
TGraphAsymmErrors *RatioGraph(TEfficiency *a, TEfficiency *b, const std::string &name)
{
  if (!a || !b) return nullptr;
  TGraphAsymmErrors *g = new TGraphAsymmErrors();
  g->SetName(name.c_str());
  const TH1 *ha = a->GetTotalHistogram();
  int n = 0;
  for (int i = 1; i <= ha->GetNbinsX(); ++i)
  {
    const double ea = a->GetEfficiency(i), eb = b->GetEfficiency(i);
    if (a->GetTotalHistogram()->GetBinContent(i) <= 0 || b->GetTotalHistogram()->GetBinContent(i) <= 0) continue;
    if (eb <= 0 || ea <= 0) continue;
    const double r  = ea / eb;
    const double ru = r * std::sqrt(std::pow(a->GetEfficiencyErrorUp(i) / ea, 2) + std::pow(b->GetEfficiencyErrorLow(i) / eb, 2));
    const double rd = r * std::sqrt(std::pow(a->GetEfficiencyErrorLow(i) / ea, 2) + std::pow(b->GetEfficiencyErrorUp(i) / eb, 2));
    const double x  = ha->GetBinCenter(i);
    const double xl = x - ha->GetBinLowEdge(i), xh = ha->GetBinLowEdge(i + 1) - x;
    g->SetPoint(n, x, r);
    g->SetPointError(n, xl, xh, rd, ru);
    ++n;
  }
  return g;
}

// "no vertical marker line" sentinel for DrawEffRatio -- must lie outside every
// x range used (the rapidity frames span [-2.4, 2.4], so -1 would NOT do).
const double kNoLine = -999.0;

struct Series
{
  TEfficiency *eff;
  std::string  label;
  int color, marker;
};
struct RatioSpec
{
  int iNum, iDen; // indices into the Series vector
  std::string label;
  int color, marker;
};

// Two-pad canvas: efficiencies on top, the requested ratios (SFs) below.
void DrawEffRatio(const std::vector<Series> &series, const std::vector<RatioSpec> &ratios,
                  const std::string &outPath, const std::string &xTitle, const std::string &ratioTitle,
                  const std::string &hdr, const std::string &sub1, const std::string &sub2,
                  const std::vector<std::string> &box,
                  double xlo, double xhi, double ylo, double yhi, double rlo, double rhi,
                  double xLine = kNoLine, const std::string &yTitle = "Trigger efficiency",
                  double legY1 = 0.37)
{
  // Layout: CMS_lumi(pad, 13, 10) paints "CMS / Work in Progress" top-LEFT
  // inside the frame (NDC y ~ 0.83-0.92), so the header block goes directly
  // below it across the full width; the legend and the info box sit lower
  // RIGHT, where the plateau points (eps ~ 1 -> NDC y ~ 0.6 for yhi = 1.6, or
  // ~0.55 on the zoomed 0.6-1.25 frames) never reach. Four-series plots pass a
  // lower legY1 and no box.
  gStyle->SetOptStat(0);
  PlotStyle ps;
  ps.headerX = 0.17; ps.headerY = 0.80; ps.headerDy = 0.045;
  ps.titleSize = 0.038; ps.subSize = 0.029;
  ps.boxX1 = 0.55; ps.boxX2 = 0.93; ps.boxY1 = 0.14; ps.boxY2 = 0.34; ps.boxTextSize = 0.027;

  static int nFrame = 0;
  TCanvas *c = new TCanvas(Form("c_eff_%d", nFrame), "", ps.w, ps.h);
  const bool haveRatio = !ratios.empty();
  const double split = haveRatio ? 0.30 : 0.0;
  TPad *pTop = new TPad("pTop", "", 0.0, split, 1.0, 1.0);
  TPad *pBot = haveRatio ? new TPad("pBot", "", 0.0, 0.0, 1.0, split) : nullptr;
  pTop->SetTopMargin(ps.tm); pTop->SetBottomMargin(haveRatio ? 0.02 : ps.bm);
  pTop->SetLeftMargin(ps.lm); pTop->SetRightMargin(ps.rm);
  pTop->SetTicks(1, 1);
  pTop->Draw();
  if (pBot)
  {
    pBot->SetTopMargin(0.04); pBot->SetBottomMargin(0.35);
    pBot->SetLeftMargin(ps.lm); pBot->SetRightMargin(ps.rm);
    pBot->SetTicks(1, 1);
    pBot->Draw();
  }

  // ---- top ----
  pTop->cd();
  TH1F *frame = new TH1F(Form("frame_eff_%d", nFrame++), "", 100, xlo, xhi);
  frame->SetDirectory(nullptr);
  ApplyHistStyle(frame, ps, haveRatio ? "" : xTitle, yTitle);
  if (haveRatio) { frame->GetXaxis()->SetLabelSize(0.0); frame->GetXaxis()->SetTitleSize(0.0); }
  frame->SetMinimum(ylo); frame->SetMaximum(yhi);
  frame->Draw("AXIS");

  TLegend *leg = new TLegend(0.55, legY1, 0.93, legY1 + 0.045 * series.size());
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(ps.font); leg->SetTextSize(0.030);
  std::vector<TGraphAsymmErrors *> gs;
  for (const auto &s : series)
  {
    TGraphAsymmErrors *g = s.eff ? s.eff->CreateGraph() : nullptr;
    gs.push_back(g);
    if (!g) continue;
    g->SetMarkerStyle(s.marker); g->SetMarkerSize(1.2);
    g->SetMarkerColor(s.color);  g->SetLineColor(s.color); g->SetLineWidth(2);
    g->Draw("P SAME");
    leg->AddEntry(g, s.label.c_str(), "lp");
  }
  leg->Draw();
  TLine *l1 = new TLine(xlo, 1.0, xhi, 1.0);
  l1->SetLineStyle(2); l1->SetLineColor(kGray + 2); l1->Draw();
  if (xLine > xlo && xLine < xhi)
  {
    TLine *lx = new TLine(xLine, ylo, xLine, yhi);
    lx->SetLineStyle(3); lx->SetLineColor(kBlue + 1); lx->SetLineWidth(2); lx->Draw();
  }
  DrawHeader(ps, hdr, sub1, sub2);
  DrawInfoBox(ps, box);
  CMS_lumi(pTop, 13, 10);
  pTop->RedrawAxis();

  // ---- bottom ----
  if (pBot)
  {
    pBot->cd();
    const double sf = 0.7 / 0.3;
    TH1F *fr = new TH1F(Form("frame_rat_%d", nFrame++), "", 100, xlo, xhi);
    fr->SetDirectory(nullptr);
    ApplyHistStyle(fr, ps, xTitle, ratioTitle);
    fr->GetXaxis()->SetTitleSize(ps.xTitleSize * sf); fr->GetYaxis()->SetTitleSize(ps.yTitleSize * sf);
    fr->GetXaxis()->SetLabelSize(ps.xLabelSize * sf); fr->GetYaxis()->SetLabelSize(ps.yLabelSize * sf);
    fr->GetXaxis()->SetTitleOffset(1.0); fr->GetYaxis()->SetTitleOffset(ps.yTitleOffset / sf);
    fr->GetYaxis()->SetNdivisions(505);
    fr->SetMinimum(rlo); fr->SetMaximum(rhi);
    fr->Draw("AXIS");
    TLine *lr = new TLine(xlo, 1.0, xhi, 1.0);
    lr->SetLineStyle(2); lr->SetLineColor(kRed + 1); lr->Draw();
    if (xLine > xlo && xLine < xhi)
    {
      TLine *lx = new TLine(xLine, rlo, xLine, rhi);
      lx->SetLineStyle(3); lx->SetLineColor(kBlue + 1); lx->SetLineWidth(2); lx->Draw();
    }
    TLegend *lg2 = ratios.size() > 1 ? new TLegend(0.55, 0.72, 0.93, 0.72 + 0.10 * ratios.size()) : nullptr;
    if (lg2) { lg2->SetBorderSize(0); lg2->SetFillStyle(0); lg2->SetTextFont(ps.font); lg2->SetTextSize(0.030 * sf); }
    for (const auto &r : ratios)
    {
      if (r.iNum >= (int)series.size() || r.iDen >= (int)series.size()) continue;
      TGraphAsymmErrors *g = RatioGraph(series[r.iNum].eff, series[r.iDen].eff, "r_tmp");
      if (!g) continue;
      g->SetMarkerStyle(r.marker); g->SetMarkerSize(1.2);
      g->SetMarkerColor(r.color);  g->SetLineColor(r.color); g->SetLineWidth(2);
      g->Draw("P SAME");
      if (lg2) lg2->AddEntry(g, r.label.c_str(), "lp");
    }
    if (lg2) lg2->Draw();
    pBot->RedrawAxis();
  }

  c->SaveAs((outPath + ".png").c_str());
  c->SaveAs((outPath + ".pdf").c_str());
  delete c;
}

// Sum of raw-count histograms of several samples (MC combined = Wp + Wm raw).
TH1D *SumH1(const std::vector<const TrigOut *> &outs, const std::string &key, const std::string &name)
{
  TH1D *h = nullptr;
  for (const TrigOut *o : outs)
  {
    TH1D *src = o->H1(key);
    if (!h) { h = (TH1D *)src->Clone(name.c_str()); h->SetDirectory(nullptr); }
    else h->Add(src);
  }
  return h;
}
TH2D *SumH2(const std::vector<const TrigOut *> &outs, const std::string &key, const std::string &name)
{
  TH2D *h = nullptr;
  for (const TrigOut *o : outs)
  {
    TH2D *src = o->H2(key);
    if (!h) { h = (TH2D *)src->Clone(name.c_str()); h->SetDirectory(nullptr); }
    else h->Add(src);
  }
  return h;
}

std::string FmtEff(TEfficiency *e, int bin)
{
  if (!e || e->GetTotalHistogram()->GetBinContent(bin) <= 0) return "      --          ";
  return Form("%.4f -%.4f +%.4f", e->GetEfficiency(bin), e->GetEfficiencyErrorLow(bin), e->GetEfficiencyErrorUp(bin));
}

// Inclusive efficiency over pT > 25 from the pT histograms (bins >= the 25 edge).
struct Incl { double den = 0, num = 0, eff = 0, lo = 0, hi = 0; };
Incl InclusiveAbove(TH1 *num, TH1 *den, double ptMin)
{
  Incl r;
  const int b0 = den->GetXaxis()->FindBin(ptMin + 1e-6);
  for (int i = b0; i <= den->GetNbinsX(); ++i) { r.den += den->GetBinContent(i); r.num += num->GetBinContent(i); }
  if (r.den > 0)
  {
    r.eff = r.num / r.den;
    r.lo  = r.eff - TEfficiency::ClopperPearson((int)r.den, (int)r.num, kCL, false);
    r.hi  = TEfficiency::ClopperPearson((int)r.den, (int)r.num, kCL, true) - r.eff;
  }
  return r;
}

// ============================================================
// Report: tables, TEfficiency objects, SF graphs, plots
// ============================================================
void Report(bool isMu, const TrigOut &oData, const TrigOut &oWp, const TrigOut &oWm, TFile *fout)
{
  const std::string ch      = isMu ? "mu" : "ele";
  const std::string lepSym  = isMu ? "#mu" : "e";
  const std::string procLbl = isMu ? "W #rightarrow #mu#nu" : "W #rightarrow e#nu";
  const std::string pathLbl = isMu ? "HLT_OxyL1SingleMuOpen" : "HLT_OxyL1SingleEG10";
  const std::string outDir  = "plots/trig_eff_mb_" + ch;
  gSystem->mkdir(outDir.c_str(), kTRUE);

  const std::vector<const TrigOut *> mcs = {&oWp, &oWm};

  std::cout << "\n==================== TRIGGER EFFICIENCY REPORT (" << ch << ") ====================\n";
  std::cout << "den = W selection (steps 1,2,4,5,6,7) && " << kMBPath << "\n";
  std::cout << "num = den && " << pathLbl << "_v1 && leading-lepton match DR < " << kTrigMatchDR << "\n";
  std::cout << "MC = Wp + Wm signal, RAW counts (gen weight ~ constant; weighted eff printed as a check)\n";

  for (int s = 0; s < 2; ++s)
  {
    const std::string S = kSel[s];
    std::cout << "\n---------- selection: " << S << "  (" << kSelLabel[s] << ") ----------\n";

    // ---- combined MC histograms ----
    std::map<std::string, TH1D *> M; // key -> summed MC histo
    for (const char *st : {"all", "den", "bit", "num", "trg"})
      for (const char *suf : {"", "_plus", "_minus"})
      {
        const std::string kpt = std::string(st) + "_pt_" + S + suf;
        const std::string ky  = std::string(st) + "_y_" + S + suf;
        const std::string ka  = std::string(st) + "_absy_" + S + suf;
        M[kpt] = SumH1(mcs, kpt, "h_" + kpt + "_mc");
        M[ky]  = SumH1(mcs, ky, "h_" + ky + "_mc");
        M[ka]  = SumH1(mcs, ka, "h_" + ka + "_mc");
      }
    M["den_pt_" + S + "_w"] = SumH1(mcs, "den_pt_" + S + "_w", "h_den_pt_" + S + "_w_mc");
    M["num_pt_" + S + "_w"] = SumH1(mcs, "num_pt_" + S + "_w", "h_num_pt_" + S + "_w_mc");
    TH2D *mDen2 = SumH2(mcs, "den_2d_" + S, "h_den_2d_" + S + "_mc");
    TH2D *mNum2 = SumH2(mcs, "num_2d_" + S, "h_num_2d_" + S + "_mc");

    // ---- efficiencies ----
    auto effD = [&](const std::string &num, const std::string &den, const std::string &n)
    { return MakeEff(oData.H1(num), oData.H1(den), n); };
    auto effM = [&](const std::string &num, const std::string &den, const std::string &n)
    { return MakeEff(M[num], M[den], n); };

    TEfficiency *eD_pt  = effD("num_pt_" + S, "den_pt_" + S, "eff_pt_" + S + "_data");
    TEfficiency *eM_pt  = effM("num_pt_" + S, "den_pt_" + S, "eff_pt_" + S + "_mc");
    TEfficiency *eWp_pt = MakeEff(oWp.H1("num_pt_" + S), oWp.H1("den_pt_" + S), "eff_pt_" + S + "_Wp");
    TEfficiency *eWm_pt = MakeEff(oWm.H1("num_pt_" + S), oWm.H1("den_pt_" + S), "eff_pt_" + S + "_Wm");
    TEfficiency *eD_ptP = effD("num_pt_" + S + "_plus",  "den_pt_" + S + "_plus",  "eff_pt_" + S + "_plus_data");
    TEfficiency *eD_ptM = effD("num_pt_" + S + "_minus", "den_pt_" + S + "_minus", "eff_pt_" + S + "_minus_data");
    TEfficiency *eM_ptP = effM("num_pt_" + S + "_plus",  "den_pt_" + S + "_plus",  "eff_pt_" + S + "_plus_mc");
    TEfficiency *eM_ptM = effM("num_pt_" + S + "_minus", "den_pt_" + S + "_minus", "eff_pt_" + S + "_minus_mc");
    TEfficiency *eD_y   = effD("num_y_" + S, "den_y_" + S, "eff_y_" + S + "_data");
    TEfficiency *eM_y   = effM("num_y_" + S, "den_y_" + S, "eff_y_" + S + "_mc");
    TEfficiency *eD_yP  = effD("num_y_" + S + "_plus",  "den_y_" + S + "_plus",  "eff_y_" + S + "_plus_data");
    TEfficiency *eD_yM  = effD("num_y_" + S + "_minus", "den_y_" + S + "_minus", "eff_y_" + S + "_minus_data");
    TEfficiency *eM_yP  = effM("num_y_" + S + "_plus",  "den_y_" + S + "_plus",  "eff_y_" + S + "_plus_mc");
    TEfficiency *eM_yM  = effM("num_y_" + S + "_minus", "den_y_" + S + "_minus", "eff_y_" + S + "_minus_mc");
    // HLT bit only (no object match) -- separates the path from the matching
    TEfficiency *eD_bit = effD("bit_pt_" + S, "den_pt_" + S, "effbit_pt_" + S + "_data");
    TEfficiency *eM_bit = effM("bit_pt_" + S, "den_pt_" + S, "effbit_pt_" + S + "_mc");
    // MB-fraction controls: den/all (no trigger requirement) and num/trg (triggered)
    TEfficiency *fD_all = effD("den_pt_" + S, "all_pt_" + S, "mbfrac_all_pt_" + S + "_data");
    TEfficiency *fD_trg = effD("num_pt_" + S, "trg_pt_" + S, "mbfrac_trg_pt_" + S + "_data");
    TEfficiency *fM_all = effM("den_pt_" + S, "all_pt_" + S, "mbfrac_all_pt_" + S + "_mc");
    TEfficiency *fM_trg = effM("num_pt_" + S, "trg_pt_" + S, "mbfrac_trg_pt_" + S + "_mc");
    TEfficiency *e2D_D  = MakeEff(oData.H2("num_2d_" + S), oData.H2("den_2d_" + S), "eff_2d_" + S + "_data");
    TEfficiency *e2D_M  = MakeEff(mNum2, mDen2, "eff_2d_" + S + "_mc");

    // |y|-folded twins -- the rapidity record skim/muon_sf.h reads and prints
    TEfficiency *eD_ay  = effD("num_absy_" + S, "den_absy_" + S, "eff_absy_" + S + "_data");
    TEfficiency *eM_ay  = effM("num_absy_" + S, "den_absy_" + S, "eff_absy_" + S + "_mc");
    TEfficiency *eD_ayP = effD("num_absy_" + S + "_plus",  "den_absy_" + S + "_plus",  "eff_absy_" + S + "_plus_data");
    TEfficiency *eD_ayM = effD("num_absy_" + S + "_minus", "den_absy_" + S + "_minus", "eff_absy_" + S + "_minus_data");
    TEfficiency *eM_ayP = effM("num_absy_" + S + "_plus",  "den_absy_" + S + "_plus",  "eff_absy_" + S + "_plus_mc");
    TEfficiency *eM_ayM = effM("num_absy_" + S + "_minus", "den_absy_" + S + "_minus", "eff_absy_" + S + "_minus_mc");

    TGraphAsymmErrors *sf_pt = RatioGraph(eD_pt, eM_pt, "sf_pt_" + S);
    TGraphAsymmErrors *sf_y  = RatioGraph(eD_y, eM_y, "sf_y_" + S);
    TGraphAsymmErrors *sf_yP = RatioGraph(eD_yP, eM_yP, "sf_y_" + S + "_plus");
    TGraphAsymmErrors *sf_yM = RatioGraph(eD_yM, eM_yM, "sf_y_" + S + "_minus");
    TGraphAsymmErrors *sf_ay  = RatioGraph(eD_ay, eM_ay, "sf_absy_" + S);
    TGraphAsymmErrors *sf_ayP = RatioGraph(eD_ayP, eM_ayP, "sf_absy_" + S + "_plus");
    TGraphAsymmErrors *sf_ayM = RatioGraph(eD_ayM, eM_ayM, "sf_absy_" + S + "_minus");

    // ---- per-bin table (pT) ----
    std::cout << Form("%-12s | %9s %9s %-24s | %9s %9s %-24s | %-22s | %s\n",
                      "pT bin", "data den", "data num", "eff_data", "MC den", "MC num", "eff_MC", "SF = data/MC", "MB frac (data, all|trg)");
    std::ofstream csv(outDir + "/trig_eff_" + S + ".csv");
    csv << "ptlo,pthi,data_den,data_num,eff_data,eff_data_lo,eff_data_hi,mc_den,mc_num,eff_mc,eff_mc_lo,eff_mc_hi,sf,sf_lo,sf_hi\n";
    for (int i = 1; i <= kNPt; ++i)
    {
      const double lo = kPtEdges[i - 1], hi = kPtEdges[i];
      std::string sfStr = "        --          ";
      double sfv = 0, sfl = 0, sfh = 0;
      if (sf_pt)
        for (int k = 0; k < sf_pt->GetN(); ++k)
          if (std::abs(sf_pt->GetX()[k] - 0.5 * (lo + hi)) < 1e-6)
          {
            sfv = sf_pt->GetY()[k]; sfl = sf_pt->GetEYlow()[k]; sfh = sf_pt->GetEYhigh()[k];
            sfStr = Form("%.4f -%.4f +%.4f", sfv, sfl, sfh);
          }
      const double fa = fD_all ? fD_all->GetEfficiency(i) : 0, ft = fD_trg ? fD_trg->GetEfficiency(i) : 0;
      std::cout << Form("%5.1f-%-6.1f | %9.0f %9.0f %-24s | %9.0f %9.0f %-24s | %-22s | %.3f | %.3f\n",
                        lo, hi,
                        oData.H1("den_pt_" + S)->GetBinContent(i), oData.H1("num_pt_" + S)->GetBinContent(i), FmtEff(eD_pt, i).c_str(),
                        M["den_pt_" + S]->GetBinContent(i), M["num_pt_" + S]->GetBinContent(i), FmtEff(eM_pt, i).c_str(),
                        sfStr.c_str(), fa, ft);
      csv << Form("%g,%g,%.0f,%.0f,%.6f,%.6f,%.6f,%.0f,%.0f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n", lo, hi,
                  oData.H1("den_pt_" + S)->GetBinContent(i), oData.H1("num_pt_" + S)->GetBinContent(i),
                  eD_pt ? eD_pt->GetEfficiency(i) : 0, eD_pt ? eD_pt->GetEfficiencyErrorLow(i) : 0, eD_pt ? eD_pt->GetEfficiencyErrorUp(i) : 0,
                  M["den_pt_" + S]->GetBinContent(i), M["num_pt_" + S]->GetBinContent(i),
                  eM_pt ? eM_pt->GetEfficiency(i) : 0, eM_pt ? eM_pt->GetEfficiencyErrorLow(i) : 0, eM_pt ? eM_pt->GetEfficiencyErrorUp(i) : 0,
                  sfv, sfl, sfh);
    }
    csv.close();

    // ---- per-bin table (y, pT > 25) ----
    std::cout << Form("\n%-14s | %9s %9s %-24s | %9s %9s %-24s | %-22s\n",
                      "y bin (pT>25)", "data den", "data num", "eff_data", "MC den", "MC num", "eff_MC", "SF = data/MC");
    for (int i = 1; i <= kNY; ++i)
    {
      std::string sfStr = "        --          ";
      if (sf_y)
        for (int k = 0; k < sf_y->GetN(); ++k)
          if (std::abs(sf_y->GetX()[k] - 0.5 * (kYEdges[i - 1] + kYEdges[i])) < 1e-6)
            sfStr = Form("%.4f -%.4f +%.4f", sf_y->GetY()[k], sf_y->GetEYlow()[k], sf_y->GetEYhigh()[k]);
      std::cout << Form("%5.2f..%-6.2f | %9.0f %9.0f %-24s | %9.0f %9.0f %-24s | %-22s\n",
                        kYEdges[i - 1], kYEdges[i],
                        oData.H1("den_y_" + S)->GetBinContent(i), oData.H1("num_y_" + S)->GetBinContent(i), FmtEff(eD_y, i).c_str(),
                        M["den_y_" + S]->GetBinContent(i), M["num_y_" + S]->GetBinContent(i), FmtEff(eM_y, i).c_str(), sfStr.c_str());
    }

    // ---- per-bin table (|y|, pT > 25) -- the rapidity record; skim/muon_sf.h
    //      applies the INCLUSIVE SF, this is what justifies that (see header) ----
    std::cout << Form("\n%-14s | %9s %9s %-24s | %9s %9s %-24s | %-22s\n",
                      "|y| bin (pT>25)", "data den", "data num", "eff_data", "MC den", "MC num", "eff_MC", "SF = data/MC (not appl.)");
    for (int i = 1; i <= kNAbsY; ++i)
    {
      std::string sfStr = "        --          ";
      if (sf_ay && i - 1 < sf_ay->GetN())
        sfStr = Form("%.4f -%.4f +%.4f", sf_ay->GetY()[i - 1], sf_ay->GetEYlow()[i - 1], sf_ay->GetEYhigh()[i - 1]);
      std::cout << Form("%5.2f..%-6.2f | %9.0f %9.0f %-24s | %9.0f %9.0f %-24s | %-22s\n",
                        kAbsYEdges[i - 1], kAbsYEdges[i],
                        oData.H1("den_absy_" + S)->GetBinContent(i), oData.H1("num_absy_" + S)->GetBinContent(i), FmtEff(eD_ay, i).c_str(),
                        M["den_absy_" + S]->GetBinContent(i), M["num_absy_" + S]->GetBinContent(i), FmtEff(eM_ay, i).c_str(), sfStr.c_str());
    }

    // ---- BINNING DECISION: likelihood-ratio flatness tests (see FlatnessLRT) ----
    std::cout << "\n  BINNING DECISION (likelihood-ratio tests vs one flat SF; pT bins restricted to the pT > 25 analysis range)\n";
    {
      // pT: restrict both legs to bins at/above the 25 GeV edge by zeroing the rest
      auto above25 = [&](const TH1 *h, const char *nm)
      {
        TH1D *c = (TH1D *)h->Clone(Form("%s_a25_%s", nm, S.c_str()));
        c->SetDirectory(nullptr);
        for (int i = 1; i <= c->GetNbinsX(); ++i)
          if (c->GetXaxis()->GetBinLowEdge(i) < kPtNominal - 1e-6) c->SetBinContent(i, 0);
        return c;
      };
      TH1D *dDp = above25(oData.H1("den_pt_" + S), "dd"), *dNp = above25(oData.H1("num_pt_" + S), "dn");
      TH1D *mDp = above25(M["den_pt_" + S], "md"),        *mNp = above25(M["num_pt_" + S], "mn");
      PrintFlatness("SF vs pT (11 bins, pT>25)", dDp, dNp, mDp, mNp);
      delete dDp; delete dNp; delete mDp; delete mNp;

      // Coarse pT (< 35 vs > 35) -- THE decision-relevant pT test. The 11-bin
      // scan above answers a different question ("is any single pT bin
      // anomalous?") and is driven by the [50,60) bin alone, a ~2 sigma
      // downward fluctuation on ~150 events present in both selections.
      const TH2 *dD2 = oData.H2("den_2d_" + S), *dN2 = oData.H2("num_2d_" + S);
      int ptLo = 1;
      while (ptLo <= dD2->GetNbinsX() && dD2->GetXaxis()->GetBinLowEdge(ptLo) < kPtNominal - 1e-6) ++ptLo;
      int ptSplit = ptLo;
      while (ptSplit <= dD2->GetNbinsX() && dD2->GetXaxis()->GetBinUpEdge(ptSplit) < 35.0 - 1e-6) ++ptSplit;
      const FlatFit fPt2 = FlatnessLRT(Cells2D(dD2, dN2, mDen2, mNum2, ptLo, ptSplit, 2, 1));
      std::cout << Form("  %-30s common SF %.4f   -2dlnL = %6.2f / %2d dof   p = %.4f  %s\n",
                        "SF vs pT (2 groups, </> 35)", fPt2.sf, fPt2.lrt, fPt2.ndf,
                        TMath::Prob(fPt2.lrt, fPt2.ndf),
                        TMath::Prob(fPt2.lrt, fPt2.ndf) < 0.05 ? "<-- NOT flat" : "(no pT dependence)");
      // Granularity: are 3 coarse |y| bins enough, or is the 6-bin structure real?
      const FlatFit fA3 = FlatnessLRT(Cells2D(dD2, dN2, mDen2, mNum2, ptLo, ptSplit, 1, 3));
      const FlatFit fA6 = FlatnessLRT(Cells2D(dD2, dN2, mDen2, mNum2, ptLo, ptSplit, 1, 6));
      {
        const double d = fA6.lrt - fA3.lrt;
        const int    nd = std::max(fA6.ndf - fA3.ndf, 1);
        std::cout << Form("  %-30s %35s -2dlnL = %6.2f / %2d dof   p = %.4f  %s\n",
                          "|y| 3 -> |y| 6 bins", "", d, nd, TMath::Prob(std::max(d, 0.0), nd),
                          TMath::Prob(std::max(d, 0.0), nd) < 0.05 ? "<-- 3 bins NOT enough" : "(3 bins would do)");
      }
      // Interaction: does |y| structure depend on pT?  |y|3 -> |y|3 x pT2
      const FlatFit fA3P2 = FlatnessLRT(Cells2D(dD2, dN2, mDen2, mNum2, ptLo, ptSplit, 2, 3));
      const double dI = fA3P2.lrt - fA3.lrt;
      const int    nI = std::max(fA3P2.ndf - fA3.ndf, 1);
      std::cout << Form("  %-30s %35s -2dlnL = %6.2f / %2d dof   p = %.4f  %s\n",
                        "|y|3 -> |y|3 x pT2 (2D)", "", dI, nI, TMath::Prob(std::max(dI, 0.0), nI),
                        TMath::Prob(std::max(dI, 0.0), nI) < 0.05 ? "<-- 2D needed" : "(1D in |y| suffices)");
    }
    PrintFlatness("SF vs |y| (6 bins)  MEASURED, NOT APPLIED", oData.H1("den_absy_" + S), oData.H1("num_absy_" + S),
                  M["den_absy_" + S], M["num_absy_" + S]);
    PrintFlatness("SF vs signed y (12 bins)", oData.H1("den_y_" + S), oData.H1("num_y_" + S),
                  M["den_y_" + S], M["num_y_" + S]);
    {
      // left/right asymmetry: is the signed-y structure more than the folded one?
      const FlatFit f6 = FlatnessLRT(oData.H1("den_absy_" + S), oData.H1("num_absy_" + S), M["den_absy_" + S], M["num_absy_" + S]);
      const FlatFit f12 = FlatnessLRT(oData.H1("den_y_" + S), oData.H1("num_y_" + S), M["den_y_" + S], M["num_y_" + S]);
      const double d = f12.lrt - f6.lrt;
      const int    nd = f12.ndf - f6.ndf;
      std::cout << Form("  %-30s %35s -2dlnL = %6.2f / %2d dof   p = %.4f  %s\n",
                        "|y| 6 -> signed y 12", "", d, nd, TMath::Prob(std::max(d, 0.0), std::max(nd, 1)),
                        TMath::Prob(std::max(d, 0.0), std::max(nd, 1)) < 0.05 ? "<-- folding NOT justified" : "(folding justified)");
      // charge: same test on the per-charge |y| tables
      const FlatFit fp = FlatnessLRT(oData.H1("den_absy_" + S + "_plus"), oData.H1("num_absy_" + S + "_plus"),
                                     M["den_absy_" + S + "_plus"], M["num_absy_" + S + "_plus"]);
      const FlatFit fm = FlatnessLRT(oData.H1("den_absy_" + S + "_minus"), oData.H1("num_absy_" + S + "_minus"),
                                     M["den_absy_" + S + "_minus"], M["num_absy_" + S + "_minus"]);
      const char *lp = isMu ? "mu" : "e ";
      std::cout << Form("  %-30s SF(%s+) %.4f  SF(%s-) %.4f  (the SF is charge-inclusive; see the charge split above)\n",
                        "per-charge |y| common SFs", lp, fp.sf, lp, fm.sf);
    }

    // ---- inclusive pT > 25 ----
    const Incl iD  = InclusiveAbove(oData.H1("num_pt_" + S), oData.H1("den_pt_" + S), kPtNominal);
    const Incl iM  = InclusiveAbove(M["num_pt_" + S], M["den_pt_" + S], kPtNominal);
    const Incl iDP = InclusiveAbove(oData.H1("num_pt_" + S + "_plus"),  oData.H1("den_pt_" + S + "_plus"),  kPtNominal);
    const Incl iDM = InclusiveAbove(oData.H1("num_pt_" + S + "_minus"), oData.H1("den_pt_" + S + "_minus"), kPtNominal);
    const Incl iMP = InclusiveAbove(M["num_pt_" + S + "_plus"],  M["den_pt_" + S + "_plus"],  kPtNominal);
    const Incl iMM = InclusiveAbove(M["num_pt_" + S + "_minus"], M["den_pt_" + S + "_minus"], kPtNominal);
    const Incl iWp = InclusiveAbove(oWp.H1("num_pt_" + S), oWp.H1("den_pt_" + S), kPtNominal);
    const Incl iWm = InclusiveAbove(oWm.H1("num_pt_" + S), oWm.H1("den_pt_" + S), kPtNominal);
    const Incl iBitD = InclusiveAbove(oData.H1("bit_pt_" + S), oData.H1("den_pt_" + S), kPtNominal);
    const Incl iBitM = InclusiveAbove(M["bit_pt_" + S], M["den_pt_" + S], kPtNominal);
    // gen-weighted MC (normal-approx binomial error)
    double wDen = 0, wNum = 0;
    {
      TH1D *hd = M["den_pt_" + S + "_w"], *hn = M["num_pt_" + S + "_w"];
      const int b0 = hd->GetXaxis()->FindBin(kPtNominal + 1e-6);
      for (int i = b0; i <= hd->GetNbinsX(); ++i) { wDen += hd->GetBinContent(i); wNum += hn->GetBinContent(i); }
    }
    // k_s-weighted Wp/Wm combination (physical charge mix)
    const double kWpS = pONorm::MCScale(isMu ? "Wp_mu" : "Wp_ele");
    const double kWmS = pONorm::MCScale(isMu ? "Wm_mu" : "Wm_ele");
    const double ksDen = kWpS * iWp.den + kWmS * iWm.den, ksNum = kWpS * iWp.num + kWmS * iWm.num;
    const double sfIncl = (iM.eff > 0) ? iD.eff / iM.eff : 0;
    const double sfErr  = (iM.eff > 0 && iD.eff > 0) ? sfIncl * std::sqrt(std::pow(0.5 * (iD.lo + iD.hi) / iD.eff, 2) + std::pow(0.5 * (iM.lo + iM.hi) / iM.eff, 2)) : 0;
    // MB-bit fractions from the stored histograms (so plotsOnly reproduces them):
    // of ALL selected events (den/all) and of the TRIGGERED ones (num/trg).
    const Incl fSelD = InclusiveAbove(oData.H1("den_pt_" + S), oData.H1("all_pt_" + S), kPtNominal);
    const Incl fTrgD = InclusiveAbove(oData.H1("num_pt_" + S), oData.H1("trg_pt_" + S), kPtNominal);
    const Incl fSelM = InclusiveAbove(M["den_pt_" + S], M["all_pt_" + S], kPtNominal);
    const double mbFracSel = fSelD.eff, mbFracTrg = fTrgD.eff;

    std::cout << "\n[INCL] " << S << " pT > " << kPtNominal << ":\n";
    std::cout << Form("  data : den %.0f  num %.0f  eff = %.4f -%.4f +%.4f   (bit only %.4f;  W+ %.4f -%.4f +%.4f [%.0f]  W- %.4f -%.4f +%.4f [%.0f])\n",
                      iD.den, iD.num, iD.eff, iD.lo, iD.hi, iBitD.eff, iDP.eff, iDP.lo, iDP.hi, iDP.den, iDM.eff, iDM.lo, iDM.hi, iDM.den);
    std::cout << Form("  MC   : den %.0f  num %.0f  eff = %.4f -%.4f +%.4f   (bit only %.4f;  W+ %.4f -%.4f +%.4f [%.0f]  W- %.4f -%.4f +%.4f [%.0f])\n",
                      iM.den, iM.num, iM.eff, iM.lo, iM.hi, iBitM.eff, iMP.eff, iMP.lo, iMP.hi, iMP.den, iMM.eff, iMM.lo, iMM.hi, iMM.den);
    std::cout << Form("  MC   : per sample Wp %.4f (den %.0f)  Wm %.4f (den %.0f);  gen-weighted %.4f;  k_s-weighted Wp+Wm %.4f\n",
                      iWp.eff, iWp.den, iWm.eff, iWm.den, wDen > 0 ? wNum / wDen : 0, ksDen > 0 ? ksNum / ksDen : 0);
    std::cout << Form("  data MB fraction of the selected events: all %.4f (%.0f/%.0f)   triggered %.4f (%.0f/%.0f)   [prescale check: should agree];  MC %.4f\n",
                      mbFracSel, fSelD.num, fSelD.den, mbFracTrg, fTrgD.num, fTrgD.den, fSelM.eff);
    std::cout << Form("[RESULT] %s %s pT>%.0f: eff_data = %.4f -%.4f +%.4f  eff_MC = %.4f -%.4f +%.4f  SF = %.4f +- %.4f\n",
                      ch.c_str(), S.c_str(), kPtNominal, iD.eff, iD.lo, iD.hi, iM.eff, iM.lo, iM.hi, sfIncl, sfErr);

    // ---- write objects ----
    fout->cd();
    for (auto &p : M) p.second->Write("", TObject::kOverwrite);
    mDen2->Write("", TObject::kOverwrite); mNum2->Write("", TObject::kOverwrite);
    for (TEfficiency *e : {eD_pt, eM_pt, eWp_pt, eWm_pt, eD_ptP, eD_ptM, eM_ptP, eM_ptM, eD_y, eM_y, eD_yP, eD_yM, eM_yP, eM_yM,
                           eD_ay, eM_ay, eD_ayP, eD_ayM, eM_ayP, eM_ayM,
                           eD_bit, eM_bit, fD_all, fD_trg, fM_all, fM_trg, e2D_D, e2D_M})
      if (e) e->Write("", TObject::kOverwrite);
    for (TGraphAsymmErrors *g : {sf_pt, sf_y, sf_yP, sf_yM, sf_ay, sf_ayP, sf_ayM})
      if (g) g->Write("", TObject::kOverwrite);

    // ---- plots ----
    const std::string hdr  = procLbl + ", " + lepSym + " trigger turn-on";
    const std::string sub1 = std::string(kSelLabel[s]) + ", |#eta| < 2.4";
    const std::string sub2 = "den: MB-triggered;  num: + " + pathLbl + ", #DeltaR < 0.4 match";
    const std::vector<std::string> box = {
        Form("p_{T} > %.0f: #varepsilon_{data} = %.3f_{-%.3f}^{+%.3f}", kPtNominal, iD.eff, iD.lo, iD.hi),
        Form("#varepsilon_{MC} = %.3f_{-%.3f}^{+%.3f}", iM.eff, iM.lo, iM.hi),
        Form("SF = %.3f #pm %.3f", sfIncl, sfErr),
        Form("N_{den}: data %.0f, MC %.0f", iD.den, iM.den)};

    // full turn-on (from the loop floor) + zoom on the analysis range
    DrawEffRatio({{eD_pt, "Data (MB-triggered)", kBlack, 20}, {eM_pt, "W signal MC", kRed + 1, 24}},
                 {{0, 1, "Data / MC", kBlack, 20}},
                 outDir + "/turnon_pt_" + S, "leading-lepton p_{T} [GeV]", "Data / MC",
                 hdr, sub1, sub2, box, kPtFloor, 120.0, 0.0, 1.6, 0.7, 1.3, kPtNominal);
    DrawEffRatio({{eD_pt, "Data (MB-triggered)", kBlack, 20}, {eM_pt, "W signal MC", kRed + 1, 24}},
                 {{0, 1, "Data / MC", kBlack, 20}},
                 outDir + "/turnon_pt_" + S + "_zoom25", "leading-lepton p_{T} [GeV]", "Data / MC",
                 hdr, sub1, sub2, box, kPtNominal, 120.0, 0.6, 1.25, 0.85, 1.15);
    // per charge
    DrawEffRatio({{eD_ptP, "Data " + lepSym + "^{+}", kBlack, 20}, {eD_ptM, "Data " + lepSym + "^{-}", kBlue + 1, 21},
                  {eM_ptP, "MC " + lepSym + "^{+}", kRed + 1, 24}, {eM_ptM, "MC " + lepSym + "^{-}", kOrange + 7, 25}},
                 {{0, 2, "Data / MC, " + lepSym + "^{+}", kBlack, 20}, {1, 3, "Data / MC, " + lepSym + "^{-}", kBlue + 1, 21}},
                 outDir + "/eff_pt_" + S + "_charge", "leading-lepton p_{T} [GeV]", "Data / MC",
                 hdr, sub1, sub2, {}, kPtNominal, 120.0, 0.6, 1.25, 0.85, 1.15, kNoLine, "Trigger efficiency", 0.16);
    // vs y (pT > 25)
    DrawEffRatio({{eD_y, "Data (MB-triggered)", kBlack, 20}, {eM_y, "W signal MC", kRed + 1, 24}},
                 {{0, 1, "Data / MC", kBlack, 20}},
                 outDir + "/eff_y_" + S, "y = -#eta_{lab}", "Data / MC",
                 hdr, sub1 + ", p_{T} > 25 GeV", sub2, box, kYEdges[0], kYEdges[kNY], 0.6, 1.25, 0.85, 1.15);
    DrawEffRatio({{eD_yP, "Data " + lepSym + "^{+}", kBlack, 20}, {eD_yM, "Data " + lepSym + "^{-}", kBlue + 1, 21},
                  {eM_yP, "MC " + lepSym + "^{+}", kRed + 1, 24}, {eM_yM, "MC " + lepSym + "^{-}", kOrange + 7, 25}},
                 {{0, 2, "Data / MC, " + lepSym + "^{+}", kBlack, 20}, {1, 3, "Data / MC, " + lepSym + "^{-}", kBlue + 1, 21}},
                 outDir + "/eff_y_" + S + "_charge", "y = -#eta_{lab}", "Data / MC",
                 hdr, sub1 + ", p_{T} > 25 GeV", sub2, {}, kYEdges[0], kYEdges[kNY], 0.6, 1.25, 0.85, 1.15, kNoLine, "Trigger efficiency", 0.16);
    // vs |y| (pT > 25) -- the ratio pad IS the rapidity-dependent SF. Measured
    // and kept as the justification for applying the inclusive one instead.
    DrawEffRatio({{eD_ay, "Data (MB-triggered)", kBlack, 20}, {eM_ay, "W signal MC", kRed + 1, 24}},
                 {{0, 1, "Data / MC = SF", kBlack, 20}},
                 outDir + "/eff_absy_" + S, "|y| = |#eta_{lab}|", "Data / MC = SF",
                 hdr, sub1 + ", p_{T} > 25 GeV", "rapidity-dependent SF (measured; an inclusive SF is applied)", box,
                 kAbsYEdges[0], kAbsYEdges[kNAbsY], 0.6, 1.25, 0.90, 1.10);
    DrawEffRatio({{eD_ayP, "Data " + lepSym + "^{+}", kBlack, 20}, {eD_ayM, "Data " + lepSym + "^{-}", kBlue + 1, 21},
                  {eM_ayP, "MC " + lepSym + "^{+}", kRed + 1, 24}, {eM_ayM, "MC " + lepSym + "^{-}", kOrange + 7, 25}},
                 {{0, 2, "Data / MC, " + lepSym + "^{+}", kBlack, 20}, {1, 3, "Data / MC, " + lepSym + "^{-}", kBlue + 1, 21}},
                 outDir + "/eff_absy_" + S + "_charge", "|y| = |#eta_{lab}|", "Data / MC",
                 hdr, sub1 + ", p_{T} > 25 GeV", sub2, {}, kAbsYEdges[0], kAbsYEdges[kNAbsY], 0.6, 1.25, 0.85, 1.15, kNoLine, "Trigger efficiency", 0.16);
    // HLT bit vs bit+match
    DrawEffRatio({{eD_bit, "Data: path fired", kGray + 2, 20}, {eD_pt, "Data: path fired + matched", kBlack, 21},
                  {eM_bit, "MC: path fired", kRed - 7, 24}, {eM_pt, "MC: path fired + matched", kRed + 1, 25}},
                 {{1, 0, "Data: matched / fired", kBlack, 20}, {3, 2, "MC: matched / fired", kRed + 1, 24}},
                 outDir + "/bit_vs_match_pt_" + S, "leading-lepton p_{T} [GeV]", "matched / fired",
                 hdr, sub1, "den: MB-triggered", {}, kPtFloor, 120.0, 0.0, 1.6, 0.8, 1.1, kPtNominal, "Trigger efficiency", 0.16);
    // MB-fraction control (prescale independence): den/all and num/trg, data only
    // (MC fires the MB path on 99.9% of W events -- quoted in the box).
    DrawEffRatio({{fD_all, "all selected events", kBlack, 20}, {fD_trg, "triggered events only", kBlue + 1, 21}},
                 {},
                 outDir + "/mbfrac_pt_" + S, "leading-lepton p_{T} [GeV]", "",
                 procLbl + ", MB-tag control (data)", sub1, "MB-bit fraction of the selected events (a pure prescale: flat, identical)",
                 {Form("p_{T} > 25: all %.4f, triggered %.4f", mbFracSel, mbFracTrg),
                  Form("MC: MB fires on %.2f%% of W events", 100.0 * fSelM.eff)},
                 kPtFloor, 120.0, 0.0, 1.6, 0.8, 1.2, kPtNominal, "MB-bit fraction");
  }

  // ---- data per-run table (nom, pT > 25) ----
  std::cout << "\n[RUNS] data, nom selection, pT > 25:  run | selected | MB (den) | num | eff | MB frac | triggered-no-MB\n";
  for (auto &p : oData.byRun)
  {
    const auto &v = p.second;
    std::cout << Form("   %d | %8.0f | %8.0f | %8.0f | %s | %.3f | %.0f\n", p.first, v[0], v[1], v[2],
                      v[1] > 0 ? Form("%.4f", v[2] / v[1]) : "  --  ", v[0] > 0 ? v[1] / v[0] : 0.0, v[3]);
  }
  std::cout << "[OK] plots in " << outDir << "/\n";
}

} // anonymous namespace

// ============================================================
// Entry point
// ============================================================
void trig_eff_mb(const char *channel = "both", bool plotsOnly = false)
{
  gSystem->mkdir("rootfile", kTRUE);
  gSystem->mkdir("plots", kTRUE);
  gSystem->mkdir("logs", kTRUE);

  const std::string chArg = channel;
  std::vector<bool> flavours;
  if (chArg == "both") flavours = {true, false};
  else if (chArg == "mu") flavours = {true};
  else if (chArg == "ele") flavours = {false};
  else { std::cerr << "[ERR] channel must be mu | ele | both\n"; return; }

  for (bool isMu : flavours)
  {
    const std::string ch = isMu ? "mu" : "ele";
    const std::string rootName = "rootfile/trig_eff_mb_" + ch + ".root";

    TrigOut oData, oWp, oWm;
    oData.Book("data_" + ch);
    oWp.Book("Wp_" + ch);
    oWm.Book("Wm_" + ch);

    if (plotsOnly)
    {
      TFile *fin = TFile::Open(rootName.c_str());
      if (!fin || fin->IsZombie()) { std::cerr << "[ERR] plotsOnly: cannot open " << rootName << "\n"; continue; }
      for (TrigOut *o : {&oData, &oWp, &oWm})
      {
        for (auto &p : o->h1)
        {
          TH1D *h = (TH1D *)fin->Get(p.second->GetName());
          if (!h) { std::cerr << "[ERR] plotsOnly: missing " << p.second->GetName() << "\n"; continue; }
          p.second->Reset(); p.second->Add(h);
        }
        for (auto &p : o->h2)
        {
          TH2D *h = (TH2D *)fin->Get(p.second->GetName());
          if (!h) { std::cerr << "[ERR] plotsOnly: missing " << p.second->GetName() << "\n"; continue; }
          p.second->Reset(); p.second->Add(h);
        }
      }
      fin->Close();
      // NB the pT>25 counters and the per-run table are not stored; they are only in the full-run log.
      TFile *fout = TFile::Open(rootName.c_str(), "UPDATE");
      Report(isMu, oData, oWp, oWm, fout);
      fout->Close();
      continue;
    }

    if (RunTrig(isMu, kData, oData) != 0) { std::cerr << "[ERR] data run failed (" << ch << ")\n"; continue; }
    if (RunTrig(isMu, kWp, oWp) != 0)     { std::cerr << "[ERR] Wp run failed (" << ch << ")\n"; continue; }
    if (RunTrig(isMu, kWm, oWm) != 0)     { std::cerr << "[ERR] Wm run failed (" << ch << ")\n"; continue; }

    TFile *fout = TFile::Open(rootName.c_str(), "RECREATE");
    oData.Write(fout); oWp.Write(fout); oWm.Write(fout);
    Report(isMu, oData, oWp, oWm, fout);
    fout->Close();
    std::cout << "[OK] wrote " << rootName << "\n";
  }
}
