// correction/gen_boost_check.C
//
// Is the W/Z hard process of the MC in the LAB frame? (2026-10-21)
//
// pO beams are asymmetric: p 6800 GeV on O 3400 GeV per nucleon (same magnetic
// rigidity, Z/A = 1/2), so the nucleon-nucleon CM frame moves along the proton
// by y = 1/2 ln(6800/3400) = 0.347 (= pOAnalysis::kPORapidityShift). A correct
// pO sample has every hard process boosted by that amount; one generated with
// symmetric 4810 + 4810 GeV beams and never boosted sits 0.347 too far toward
// the oxygen. The Aug-20 production (POWHEG + Angantyr) FAILS this check: the
// Born-level events are not boosted at all and the others only ~0.10 of the
// way (CLAUDE.md, TODO "MC rapidity boost"). Run this on every new production
// BEFORE switching to it.
//
// Modes (first argument; "all" = the three in order):
//
//   gen   GENERATOR LEVEL, MC only. hiEvtAnalyzer/HiTree::pdfX = (x1, x2), the
//         momentum fractions of the incoming partons relative to the beams
//         the generator used (x1 = the +z beam); the boson = the EventTree
//         mc* entry with |mcPID| = 23/24 and mcStatus = 62.
//         (a) Born-level events, r = sqrt(x1 x2 s)/m_V < 1.005 (exact 2->1
//             kinematics: the two partons ARE the boson): light-cone momentum
//             is conserved per side, so the beam energies the hard process saw
//             are E(+z) = (E+pz)_V/2x1 and E(-z) = (E-pz)_V/2x2, in the frame of
//             the final MC record. Lab frame = 3400 (O) / 6800 (p) GeV; the
//             symmetric NN CM frame = 4810 / 4810.
//         (b) Delta = y_V - 1/2 ln(x1/x2) = the rapidity of the NN CM frame as
//             seen in the MC record (1/2 ln(x1/x2) is the y of the partonic
//             system in the frame where both beams have equal energy per
//             nucleon, whatever frame the generator ran in). Expected
//             1/2 ln(E(+z)/E(-z)) = -0.347 with the oxygen on +z, exactly for
//             Born events, within +-ln r for events with an emission.
//         VERDICT = PASS when the Born-event Delta is within 0.02 of the
//         expectation. Also prints the x range per side (true momentum
//         fractions stay below 1), the r<1 fraction (s = (9.62 TeV)^2), the
//         oxygen side from the gen record (particles whose mother is the O16
//         nucleus) and the W+/W- lean (W+ toward the u-rich proton).
//   side  BEAM DIRECTION in data and MC. Data: hiZDCplus/hiZDCminus (the
//         oxygen's spectator neutrons hit its own ZDC). Data and MC: central
//         charged PF multiplicity per side (0.5 < |eta| < 2.4, more on the
//         oxygen side) and the HF PF sum pT per side (3 < |eta| < 5.2). The
//         analysis convention y = -eta ("p-going (-Z) = forward", skim.C)
//         needs the oxygen on +z in BOTH. VERDICT = PASS when the data ZDC and
//         the MC gen record put the oxygen on the same side.
//   scan  DO THE DATA SEE A MISSING BOOST? The prefit per-lab-bin W yields of
//         the current leppt_mt40 Combine inputs (data - z - ztau - ttbar -
//         qcd_abcd vs signal + wtau, one free normalization per flavour and
//         charge, data Poisson errors) against the MC W templates reweighted
//         per bin by the gen-level effect of an extra rapidity shift dy
//         (analysis convention, + = toward the proton). The reweighting comes
//         from the gen loop (needs mode gen in the same run). Meaningful only
//         when the Combine inputs were built from the MC under test, i.e.
//         AFTER the production switch: a correct sample gives dy ~ 0. On the
//         Aug-20 inputs (2026-10-20) mu gave +0.43 +- 0.10. Degenerate with
//         any other rapidity-dependent data/MC difference (nPDF shape,
//         efficiency vs eta), so it is supporting evidence, not a test.
//
// Inputs: MC files from pOSkim::ResolveMCSample (W+, W-, DY, muon flavour),
// data from pOSkim::kDefaultDataFile -- the production in use. To check a NEW
// production before switching, pass its files as arguments 2-4 (wrapper:
// MC_WP / MC_WM / MC_DY env vars); argument 5 overrides the data file.
//
// Outputs (run from correction/):
//   rootfile/gen_boost_check.root   per-sample histograms (Born beam energies,
//                                   Delta, Delta vs r, boson y, shift response)
//                                   + the scan graphs
//   plots/gen_boost/                born_beam_energy, delta, delta_vs_r,
//                                   data_shift_scan
//   stdout                          [RESULT]/[VERDICT] lines -- run through
//                                   ./run_gen_boost_check.sh to keep the log
//
// Run:
//   ./run_gen_boost_check.sh [gen|side|scan|all]          (default all)
//   MC_WP=<file> MC_WM=<file> MC_DY=<file> ./run_gen_boost_check.sh gen
//   root -l -b -q 'gen_boost_check.C+("gen")'             (bare)
//
// Branch gotcha: under SetBranchStatus("*", 0), pdfX needs "pdf*" AND the
// sub-branches "first" and "second" enabled -- "pdfX*" alone leaves x = 0.

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "../skim/skim_common.h"
#include "../analysis/analysis_helpers.h"
#include "../plotting/plotting_helper.C"

namespace
{
const double kGbcEp = 6800.0;                   // proton beam energy (GeV)
const double kGbcEN = 3400.0;                   // oxygen beam energy per nucleon
const double kGbcS = 4.0 * kGbcEp * kGbcEN;     // s_NN = (9.62 TeV)^2
const double kGbcShift = 0.5 * std::log(kGbcEp / kGbcEN); // 0.3466
const int kGbcOxygenPdg = 1000080160;
const double kGbcBornR = 1.005;                 // Born-level: r below this
const double kGbcSoftR = 1.05, kGbcSoftPt = 3.0;
const double kGbcPassTol = 0.02;                // |<Delta_Born> - expected| for PASS
const Long64_t kGbcSideEvents = 200000;         // events per file in mode side
// extra-shift grid of the scan (analysis convention, + = toward the proton)
const int kGbcNShift = 57;
const double kGbcShift0 = -0.4, kGbcShiftStep = 0.025;
const char *kGbcRoot = "rootfile/gen_boost_check.root";
const char *kGbcPlotDir = "plots/gen_boost";

struct GbcGen
{
  std::string tag, label;
  int pdg = 0;
  TH1D *hEp = nullptr, *hEm = nullptr;          // Born-event beam energies (+z, -z)
  TH1D *hDb = nullptr, *hDs = nullptr, *hDa = nullptr; // Delta: Born, soft emission, all
  TH1D *hR = nullptr, *hY = nullptr;            // r, boson y (analysis convention)
  TH2D *hDr = nullptr;                          // Delta vs r, pT_V < 3 GeV
  TH2D *hShift = nullptr;                       // W only: fiducial lepton y bin vs extra shift
  long nEv = 0, nBorn = 0, nNoBoson = 0, nBadX = 0, nRbelow = 0;
  long nOplus = 0, nOminus = 0;                 // particles with an O16 mother, eta>0 / eta<0
  double xMax1 = 0, xMax2 = 0;
  bool ok = false;
};

TH1D *GbcH1(const char *n, int nb, double lo, double hi)
{
  TH1D *h = new TH1D(n, "", nb, lo, hi);
  h->SetDirectory(nullptr);
  h->Sumw2();
  return h;
}

double GbcMedian(TH1D *h)
{
  if (!h || h->Integral() <= 0) return 0;
  double q = 0.5, m = 0;
  h->GetQuantiles(1, &m, &q);
  return m;
}

// One MC file: boson kinematics vs pdfX, Born-event beam energies, the oxygen
// side of the gen record and (W) the fiducial lepton response to a rapidity shift.
bool GbcLoopGen(GbcGen &g, const std::string &fname)
{
  TFile *f = TFile::Open(fname.c_str(), "READ");
  if (!f || f->IsZombie()) { printf("[ERROR] cannot open %s\n", fname.c_str()); return false; }
  TTree *t = (TTree *)f->Get("ggHiNtuplizer/EventTree");
  TTree *h = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  if (!t || !h) { printf("[ERROR] %s: missing EventTree or HiTree\n", fname.c_str()); return false; }
  if (t->GetEntries() != h->GetEntries())
  { printf("[ERROR] %s: EventTree (%lld) and HiTree (%lld) not entry-aligned\n", fname.c_str(), t->GetEntries(), h->GetEntries()); return false; }
  for (const char *b : {"mcPID", "mcStatus", "mcMomPID", "mcGMomPID", "mcMomEta", "mcPt", "mcEta", "mcPhi", "mcE", "mcMass"})
    if (!t->GetBranch(b)) { printf("[ERROR] %s: EventTree has no %s\n", fname.c_str(), b); return false; }
  if (!h->GetBranch("pdfX") || !h->GetBranch("weight"))
  { printf("[ERROR] %s: HiTree has no pdfX / weight\n", fname.c_str()); return false; }

  t->SetBranchStatus("*", 0);
  h->SetBranchStatus("*", 0);
  std::vector<int> *pid = nullptr, *st = nullptr, *mom = nullptr, *gmom = nullptr;
  std::vector<float> *pt = nullptr, *eta = nullptr, *phi = nullptr, *E = nullptr, *m = nullptr, *momEta = nullptr;
  t->SetBranchStatus("mcPID", 1);     t->SetBranchAddress("mcPID", &pid);
  t->SetBranchStatus("mcMomEta", 1);  t->SetBranchAddress("mcMomEta", &momEta);
  t->SetBranchStatus("mcStatus", 1);  t->SetBranchAddress("mcStatus", &st);
  t->SetBranchStatus("mcMomPID", 1);  t->SetBranchAddress("mcMomPID", &mom);
  t->SetBranchStatus("mcGMomPID", 1); t->SetBranchAddress("mcGMomPID", &gmom);
  t->SetBranchStatus("mcPt", 1);      t->SetBranchAddress("mcPt", &pt);
  t->SetBranchStatus("mcEta", 1);     t->SetBranchAddress("mcEta", &eta);
  t->SetBranchStatus("mcPhi", 1);     t->SetBranchAddress("mcPhi", &phi);
  t->SetBranchStatus("mcE", 1);       t->SetBranchAddress("mcE", &E);
  t->SetBranchStatus("mcMass", 1);    t->SetBranchAddress("mcMass", &m);
  std::pair<float, float> *px = nullptr;
  float w = 1;
  h->SetBranchStatus("pdf*", 1);   // pdfX is split: the sub-branches are "first"/"second"
  h->SetBranchStatus("first", 1);
  h->SetBranchStatus("second", 1);
  h->SetBranchAddress("pdfX", &px);
  h->SetBranchStatus("weight", 1);
  h->SetBranchAddress("weight", &w);

  const TString s = g.tag.c_str();
  g.hEp = GbcH1("born_Eplus_" + s, 1600, 0, 8000);   // 5 GeV: the medians come from these
  g.hEm = GbcH1("born_Eminus_" + s, 1600, 0, 8000);
  g.hDb = GbcH1("delta_born_" + s, 160, -0.6, 0.2);
  g.hDs = GbcH1("delta_soft_" + s, 160, -0.6, 0.2);
  g.hDa = GbcH1("delta_all_" + s, 200, -1.5, 1.0);
  g.hR = GbcH1("r_" + s, 220, 0.9, 3.1);
  g.hY = GbcH1("bosonY_" + s, 96, -6, 6);
  g.hDr = new TH2D("delta_vs_r_" + s, "", 60, 0.98, 1.28, 125, -0.65, 0.60);   // headroom for the header
  g.hDr->SetDirectory(nullptr);
  const bool isW = (g.pdg == 24);
  if (isW)
  {
    g.hShift = new TH2D("shift_response_" + s, "", kGbcNShift, kGbcShift0 - 0.5 * kGbcShiftStep,
                        kGbcShift0 + (kGbcNShift - 0.5) * kGbcShiftStep, pOSkim::kNY, pOSkim::kYEdges);
    g.hShift->SetDirectory(nullptr);
  }

  const Long64_t n = t->GetEntries();
  for (Long64_t i = 0; i < n; ++i)
  {
    t->GetEntry(i);
    h->GetEntry(i);
    const double ws = (w < 0) ? -1.0 : 1.0;   // POWHEG: |w| constant, sign = sign of B-bar
    ++g.nEv;
    int ib = -1, il = -1, inu = -1;
    for (size_t k = 0; k < pid->size(); ++k)
    {
      const int a = std::abs((*pid)[k]);
      if (a == g.pdg && (*st)[k] == 62) ib = (int)k;
      // the nucleus itself (pT = 0) carries its direction in the MOTHER's eta;
      // its daughters can go anywhere
      if (std::abs((*mom)[k]) == kGbcOxygenPdg) { if ((*momEta)[k] > 0) ++g.nOplus; else ++g.nOminus; }
      if (isW && (*st)[k] == 1)
      {
        const bool fromW = std::abs((*mom)[k]) == 24 || (std::abs((*gmom)[k]) == 24 && std::abs((*mom)[k]) == a);
        if (fromW && (a == 11 || a == 13) && (il < 0 || (*pt)[k] > (*pt)[il])) il = (int)k;
        if (fromW && (a == 12 || a == 14) && (inu < 0 || (*pt)[k] > (*pt)[inu])) inu = (int)k;
      }
    }
    if (ib < 0) { ++g.nNoBoson; continue; }
    const double x1 = px->first, x2 = px->second;
    if (!(x1 > 0 && x2 > 0)) { ++g.nBadX; continue; }
    g.xMax1 = std::max(g.xMax1, x1);
    g.xMax2 = std::max(g.xMax2, x2);

    const double PT = (*pt)[ib], pz = PT * std::sinh((*eta)[ib]), EE = (*E)[ib], M = (*m)[ib];
    const double yDet = 0.5 * std::log((EE + pz) / (EE - pz));
    const double D = yDet - 0.5 * std::log(x1 / x2);
    const double r = std::sqrt(x1 * x2 * kGbcS) / M;
    g.hY->Fill(-yDet, ws);                     // analysis convention y = -y_det
    g.hR->Fill(r, ws);
    g.hDa->Fill(D, ws);
    if (r < 0.995) ++g.nRbelow;
    if (PT < kGbcSoftPt) g.hDr->Fill(r, D, ws);
    if (r < kGbcBornR)
    {
      ++g.nBorn;
      g.hEp->Fill((EE + pz) / (2.0 * x1), ws);
      g.hEm->Fill((EE - pz) / (2.0 * x2), ws);
      g.hDb->Fill(D, ws);
    }
    else if (r < kGbcSoftR && PT < kGbcSoftPt)
      g.hDs->Fill(D, ws);

    // W: the gen fiducial lepton (pT > 25, m_T(l,nu) > 40, |y| < 2.4) per lab
    // bin for every extra shift dy (a longitudinal boost shifts the lepton eta
    // by dy and leaves pT, phi and m_T alone).
    if (isW && il >= 0 && inu >= 0)
    {
      const double lpt = (*pt)[il], npt = (*pt)[inu];
      const double mt = std::sqrt(2.0 * lpt * npt * (1.0 - std::cos((*phi)[il] - (*phi)[inu])));
      if (lpt > 25.0 && mt > 40.0)
        for (int d = 0; d < kGbcNShift; ++d)
        {
          const double dy = kGbcShift0 + d * kGbcShiftStep, y = -(*eta)[il] + dy;
          if (std::fabs(y) < 2.4) g.hShift->Fill(dy, y, ws);
        }
    }
  }
  f->Close();
  g.ok = true;
  return true;
}

// A plot frame with the repo's cosmetics.
TCanvas *GbcCanvas(const char *name, PlotStyle &ps)
{
  TCanvas *c = new TCanvas(name, "", ps.w, ps.h);
  ApplyCanvasStyle(c, ps);
  return c;
}

void GbcVLine(double x, double y0, double y1, int col, int sty)
{
  TLine *l = new TLine(x, y0, x, y1);
  l->SetLineColor(col);
  l->SetLineStyle(sty);
  l->SetLineWidth(2);
  l->Draw();
}

void GbcText(double x, double y, const char *s, double sz, int col = kBlack)
{
  TLatex *t = new TLatex(x, y, s);
  t->SetNDC();
  t->SetTextFont(42);
  t->SetTextSize(sz);
  t->SetTextColor(col);
  t->Draw();
}

void GbcSave(TCanvas *c, const char *stem)
{
  CMS_lumi(c, 13, 10);
  c->RedrawAxis();
  c->SaveAs(Form("%s/%s.png", kGbcPlotDir, stem));
  c->SaveAs(Form("%s/%s.pdf", kGbcPlotDir, stem));
  delete c;
}

// ---------------------------------------------------------------- mode gen
bool GbcRunGen(std::vector<GbcGen> &gens, int &oxygenSide)
{
  printf("\n================ MODE gen: where does the W/Z hard process sit? ================\n");
  printf("[CONFIG] beams p %.0f GeV, O %.0f GeV/nucleon; sqrt(s_NN) = %.0f GeV; boost 1/2 ln(Ep/EN) = %.4f"
         " (pOAnalysis::kPORapidityShift = %.4f)\n",
         kGbcEp, kGbcEN, std::sqrt(kGbcS), kGbcShift, pOAnalysis::kPORapidityShift);
  for (auto &g : gens)
  {
    printf("[INPUT] %-3s %s\n", g.tag.c_str(), g.label.c_str());
    if (!GbcLoopGen(g, g.label)) return false;
  }

  // the oxygen side from the gen record (all samples)
  long op = 0, om = 0;
  for (auto &g : gens) { op += g.nOplus; om += g.nOminus; }
  oxygenSide = (op > 10 * om) ? +1 : (om > 10 * op) ? -1 : 0;
  printf("[RESULT] MC gen record: particles whose O16 mother has eta>0: %ld, eta<0: %ld -> oxygen on %s\n",
         op, om, oxygenSide > 0 ? "+z" : oxygenSide < 0 ? "-z" : "UNDETERMINED");
  const double eExpPlus = oxygenSide >= 0 ? kGbcEN : kGbcEp, eExpMinus = oxygenSide >= 0 ? kGbcEp : kGbcEN;
  const double dExp = 0.5 * std::log(eExpPlus / eExpMinus);
  printf("[RESULT] expected in the lab frame: E(+z) = %.0f, E(-z) = %.0f GeV, Delta = %+.4f"
         "  (unboosted, symmetric NN CM frame: %.0f / %.0f GeV, Delta = 0)\n",
         eExpPlus, eExpMinus, dExp, std::sqrt(kGbcS) / 2, std::sqrt(kGbcS) / 2);

  bool pass = true;
  for (auto &g : gens)
  {
    const double fBorn = g.nEv ? 100.0 * g.nBorn / g.nEv : 0;
    printf("[RESULT] %-3s events %ld (no status-62 boson %ld, x<=0 %ld); x max +z %.3f / -z %.3f;"
           " r<0.995 %.2f%%\n", g.tag.c_str(), g.nEv, g.nNoBoson, g.nBadX, g.xMax1, g.xMax2,
           g.nEv ? 100.0 * g.nRbelow / g.nEv : 0);
    printf("[RESULT] %-3s Born events %ld (%.2f%%): E(+z) median %.0f, E(-z) median %.0f GeV;"
           " Delta(Born) %+.4f rms %.4f; Delta(soft emission) %+.3f; <Delta>(all) %+.3f; <y_V> %+.3f (analysis y)\n",
           g.tag.c_str(), g.nBorn, fBorn, GbcMedian(g.hEp), GbcMedian(g.hEm), g.hDb->GetMean(), g.hDb->GetRMS(),
           g.hDs->GetMean(), g.hDa->GetMean(), g.hY->GetMean());
    if (g.nBorn < 100) { printf("[WARN] %s: only %ld Born events -- the check is weak\n", g.tag.c_str(), g.nBorn); }
    if (g.hDb->GetRMS() > 0.05)
      printf("[WARN] %s: Born-event Delta rms %.3f > 0.05 -- x1 may not be the +z beam in this production\n",
             g.tag.c_str(), g.hDb->GetRMS());
    if (std::max(g.xMax1, g.xMax2) > 1.0)
      printf("[WARN] %s: x > 1 seen -- the stored x are not momentum fractions of the generator beams\n", g.tag.c_str());
    if (std::fabs(g.hDb->GetMean() - dExp) > kGbcPassTol) pass = false;
  }
  // W lean: W+ goes with the u-rich proton
  const GbcGen *wp = nullptr, *wm = nullptr;
  for (auto &g : gens) { if (g.tag == "Wp") wp = &g; if (g.tag == "Wm") wm = &g; }
  if (wp && wm)
    printf("[RESULT] W lean (analysis y = -y_det, proton side +): <y> W+ %+.3f, W- %+.3f -> W+ %s of W-%s\n",
           wp->hY->GetMean(), wm->hY->GetMean(), wp->hY->GetMean() > wm->hY->GetMean() ? "forward" : "BACKWARD",
           wp->hY->GetMean() > wm->hY->GetMean() ? " (proton forward, as expected)" : " -- proton NOT on the forward side?");
  printf("[VERDICT] gen boost: %s  (Born-event Delta vs %+.3f within %.2f; 0 = not boosted at all)\n",
         pass ? "PASS" : "FAIL", dExp, kGbcPassTol);

  // ---- plots
  gSystem->mkdir(kGbcPlotDir, kTRUE);
  TH1D *ep = (TH1D *)gens[0].hEp->Clone("born_Eplus_sum"), *em = (TH1D *)gens[0].hEm->Clone("born_Eminus_sum");
  TH1D *db = (TH1D *)gens[0].hDb->Clone("delta_born_sum"), *ds = (TH1D *)gens[0].hDs->Clone("delta_soft_sum");
  for (size_t i = 1; i < gens.size(); ++i)
  { ep->Add(gens[i].hEp); em->Add(gens[i].hEm); db->Add(gens[i].hDb); ds->Add(gens[i].hDs); }
  std::string samples;
  for (auto &g : gens) samples += (samples.empty() ? "" : " + ") + g.tag;

  // Layout: CMS_lumi(..., 10) paints "CMS / Work in Progress" top LEFT inside
  // the frame, so the header block sits top RIGHT with the legend under it.
  auto style = [](PlotStyle &ps) {
    ps.lm = 0.15; ps.yTitleOffset = 1.55; ps.xTitleSize = 0.040;
    ps.headerX = 0.50; ps.headerY = 0.89; ps.headerDy = 0.042; ps.titleSize = 0.032; ps.subSize = 0.027;
  };
  {
    PlotStyle ps;
    style(ps);
    ps.lm = 0.18; ps.yTitleOffset = 1.85; ps.headerX = 0.53;   // five-digit y labels push the CMS label right
    TCanvas *c = GbcCanvas("c_gbc_E", ps);
    const double medP = GbcMedian(ep), medM = GbcMedian(em);   // before the display rebin
    ep->Rebin(10); em->Rebin(10);
    const double mx = 1.8 * std::max(ep->GetMaximum(), em->GetMaximum());
    ApplyHistStyle(ep, ps, "beam energy per nucleon seen by the hard process [GeV]", "Born-level events");
    ep->GetXaxis()->SetRangeUser(2800, 7400);
    ep->GetXaxis()->SetNdivisions(505);
    ep->SetMinimum(0); ep->SetMaximum(mx);
    ep->SetLineColor(kAzure + 2); ep->SetFillColorAlpha(kAzure + 2, 0.35); ep->SetLineWidth(2);
    em->SetLineColor(kRed + 1); em->SetFillColorAlpha(kRed + 1, 0.25); em->SetLineWidth(2); em->SetLineStyle(2);
    ep->Draw("hist"); em->Draw("hist same");
    GbcVLine(eExpPlus, 0, mx * 0.55, kAzure + 2, 7);
    GbcVLine(eExpMinus, 0, mx * 0.55, kRed + 1, 7);
    GbcVLine(std::sqrt(kGbcS) / 2, 0, mx * 0.55, kGray + 2, 3);
    TLegend *L = new TLegend(0.53, 0.63, 0.94, 0.74);
    L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextFont(42); L->SetTextSize(0.027);
    L->AddEntry(ep,"+z beam: (E+p_{z})_{V} / 2x_{1}", "f");
    L->AddEntry(em, "#minusz beam: (E#minusp_{z})_{V} / 2x_{2}", "f");
    L->Draw();
    DrawHeader(ps, Form("MC %s, Born level", samples.c_str()),
               Form("lab: %.0f (+z) / %.0f (#minusz) GeV", eExpPlus, eExpMinus),
               Form("MC: %.0f / %.0f GeV (NN CM: %.0f)", medP, medM, std::sqrt(kGbcS) / 2));
    GbcText(0.55, 0.45, pass ? "PASS: boosted to the lab" : "FAIL: not boosted", 0.032, pass ? kGreen + 2 : kRed + 1);
    GbcSave(c, "born_beam_energy");
  }
  {
    PlotStyle ps;
    style(ps);
    TCanvas *c = GbcCanvas("c_gbc_D", ps);
    if (db->Integral() > 0) db->Scale(1.0 / db->Integral());
    if (ds->Integral() > 0) ds->Scale(1.0 / ds->Integral());
    const double mx = 1.8 * std::max(db->GetMaximum(), ds->GetMaximum());
    ApplyHistStyle(db, ps, "#Delta = y_{V} #minus #frac{1}{2} ln(x_{1}/x_{2})", "normalized");
    db->SetMinimum(0); db->SetMaximum(mx);
    db->SetLineColor(kBlack); db->SetLineWidth(2);
    ds->SetLineColor(kOrange + 7); ds->SetLineWidth(2); ds->SetFillColorAlpha(kOrange + 7, 0.3);
    db->Draw("hist"); ds->Draw("hist same"); db->Draw("hist same");
    GbcVLine(dExp, 0, mx * 0.55, kGreen + 2, 7);
    GbcVLine(0, 0, mx * 0.55, kGray + 2, 3);
    TLegend *L = new TLegend(0.50, 0.63, 0.93, 0.74);
    L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextFont(42); L->SetTextSize(0.027);
    L->AddEntry(db, Form("Born level: %+.3f", db->GetMean()), "l");
    L->AddEntry(ds, Form("soft emission: %+.3f", ds->GetMean()), "f");
    L->Draw();
    DrawHeader(ps, Form("MC %s", samples.c_str()),
               "y of the NN CM frame (detector y)",
               Form("lab %+.3f (green), unboosted 0 (grey)", dExp));
    GbcText(0.18, 0.60, Form("soft emission: r < %.2f,", kGbcSoftR), 0.025);
    GbcText(0.18, 0.565, Form("p_{T}^{V} < %.0f GeV", kGbcSoftPt), 0.025);
    GbcSave(c, "delta");
  }
  {
    PlotStyle ps;
    style(ps);
    ps.rm = 0.16; ps.headerX = 0.48; ps.subSize = 0.025;
    TCanvas *c = GbcCanvas("c_gbc_Dr", ps);
    c->SetLogz();
    TH2D *h2 = (TH2D *)gens[0].hDr->Clone("delta_vs_r_sum");
    for (size_t i = 1; i < gens.size(); ++i) h2->Add(gens[i].hDr);
    h2->SetMinimum(1);
    ApplyHistStyle(h2, ps, "r = #sqrt{x_{1}x_{2}s} / m_{V}", "#Delta");
    h2->Draw("colz");
    for (double d0 : {0.0, dExp})
    {
      TGraph *up = new TGraph(), *lo = new TGraph();
      for (double r = 1.0; r <= 1.2801; r += 0.005) { up->SetPoint(up->GetN(), r, d0 + std::log(r)); lo->SetPoint(lo->GetN(), r, d0 - std::log(r)); }
      for (TGraph *gr : {up, lo}) { gr->SetLineColor(d0 == 0 ? kGray + 2 : kGreen + 2); gr->SetLineStyle(d0 == 0 ? 3 : 7); gr->SetLineWidth(2); gr->Draw("L same"); }
    }
    DrawHeader(ps, Form("MC %s, p_{T}^{V} < %.0f GeV", samples.c_str(), kGbcSoftPt),
               Form("lab: #Delta = %+.3f #pm ln r (green)", dExp), "unboosted: #Delta = 0 #pm ln r (grey)");
    GbcSave(c, "delta_vs_r");
  }
  return true;
}

// ---------------------------------------------------------------- mode side
struct GbcSide { double hfP = 0, hfM = 0, nP = 0, nM = 0, zP = 0, zM = 0; long n = 0, zMore = 0, zLess = 0; bool zdc = false, ok = false; };

GbcSide GbcLoopSide(const std::string &fname, bool isData)
{
  GbcSide s;
  TFile *f = TFile::Open(fname.c_str(), "READ");
  if (!f || f->IsZombie()) { printf("[ERROR] cannot open %s\n", fname.c_str()); return s; }
  TTree *p = (TTree *)f->Get("particleFlowAnalyser/pftree");
  TTree *h = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  if (!p || !h) { printf("[ERROR] %s: missing pftree or HiTree\n", fname.c_str()); return s; }
  p->SetBranchStatus("*", 0);
  h->SetBranchStatus("*", 0);
  std::vector<int> *id = nullptr;
  std::vector<float> *eta = nullptr, *pt = nullptr;
  p->SetBranchStatus("pfId", 1);  p->SetBranchAddress("pfId", &id);
  p->SetBranchStatus("pfEta", 1); p->SetBranchAddress("pfEta", &eta);
  p->SetBranchStatus("pfPt", 1);  p->SetBranchAddress("pfPt", &pt);
  float zp = 0, zm = 0;
  s.zdc = isData && h->GetBranch("hiZDCplus") && h->GetBranch("hiZDCminus");
  if (s.zdc)
  {
    h->SetBranchStatus("hiZDCplus", 1);  h->SetBranchAddress("hiZDCplus", &zp);
    h->SetBranchStatus("hiZDCminus", 1); h->SetBranchAddress("hiZDCminus", &zm);
  }
  const Long64_t n = std::min(kGbcSideEvents, p->GetEntries());
  for (Long64_t i = 0; i < n; ++i)
  {
    p->GetEntry(i);
    if (s.zdc) h->GetEntry(i);
    for (size_t k = 0; k < eta->size(); ++k)
    {
      const double e = (*eta)[k];
      if (e > 3.0 && e < 5.2) s.hfP += (*pt)[k];
      if (e < -3.0 && e > -5.2) s.hfM += (*pt)[k];
      const int t = (*id)[k];       // PF: 1 h+-, 2 e, 3 mu
      if (t >= 1 && t <= 3 && (*pt)[k] > 0.3)
      {
        if (e > 0.5 && e < 2.4) s.nP += 1;
        if (e < -0.5 && e > -2.4) s.nM += 1;
      }
    }
    if (s.zdc) { s.zP += zp; s.zM += zm; if (zp > zm) ++s.zMore; else if (zm > zp) ++s.zLess; }
  }
  s.n = n;
  s.ok = true;
  f->Close();
  return s;
}

bool GbcRunSide(const std::string &dataFile, const std::vector<GbcGen> &gens, int mcOxygenSide)
{
  printf("\n================ MODE side: which way do the beams go, in data and MC? ================\n");
  printf("[INPUT] data %s (first %lld events)\n", dataFile.c_str(), kGbcSideEvents);
  GbcSide d = GbcLoopSide(dataFile, true);
  if (!d.ok) return false;
  int dataOxygenSide = 0;
  if (d.zdc)
  {
    dataOxygenSide = d.zP > 3 * d.zM ? +1 : d.zM > 3 * d.zP ? -1 : 0;
    printf("[RESULT] data ZDC energy +z/-z = %.2f (ZDC+ larger in %.1f%% of events) -> oxygen on %s\n",
           d.zM > 0 ? d.zP / d.zM : 0, 100.0 * d.zMore / std::max(1L, d.zMore + d.zLess),
           dataOxygenSide > 0 ? "+z" : dataOxygenSide < 0 ? "-z" : "UNDETERMINED");
  }
  else
    printf("[WARN] data HiTree has no hiZDCplus/hiZDCminus -- direction from the multiplicity only\n");
  printf("[RESULT] data: charged PF multiplicity 0.5<|eta|<2.4 +z/-z = %.3f; HF PF sum pT +z/-z = %.3f\n",
         d.nP / d.nM, d.hfP / d.hfM);
  if (!d.zdc) dataOxygenSide = d.nP > d.nM ? +1 : -1;

  for (auto &g : gens)
  {
    GbcSide m = GbcLoopSide(g.label, false);
    if (!m.ok) return false;
    printf("[RESULT] MC %-3s: charged PF multiplicity 0.5<|eta|<2.4 +z/-z = %.3f; HF PF sum pT +z/-z = %.3f\n",
           g.tag.c_str(), m.nP / m.nM, m.hfP / m.hfM);
  }
  if (mcOxygenSide == 0)
    printf("[WARN] MC oxygen side unknown (mode gen not run, or no O16 mothers in the gen record)\n");
  const bool agree = mcOxygenSide != 0 && mcOxygenSide == dataOxygenSide;
  printf("[VERDICT] beam direction: %s  (oxygen data %s, MC %s; the analysis needs the oxygen on +z in both, "
         "i.e. the proton forward after y = -eta)\n",
         agree && dataOxygenSide > 0 ? "PASS" : "FAIL",
         dataOxygenSide > 0 ? "+z" : dataOxygenSide < 0 ? "-z" : "?", mcOxygenSide > 0 ? "+z" : mcOxygenSide < 0 ? "-z" : "?");
  return true;
}

// ---------------------------------------------------------------- mode scan
TGraph *GbcScanOne(const char *inFile, const char *charge, TH2D *resp, bool &ok)
{
  ok = false;
  TFile f(inFile, "READ");
  if (f.IsZombie()) { printf("[WARN] cannot open %s\n", inFile); return nullptr; }
  double D[pOSkim::kNY], B[pOSkim::kNY], W[pOSkim::kNY];
  for (int i = 0; i < pOSkim::kNY; ++i)
  {
    const TString dir = Form("%s_lab_y%d/", charge, i);
    auto I = [&](const char *proc) { TH1 *h = (TH1 *)f.Get(dir + proc); return h ? h->Integral(1, h->GetNbinsX()) : -1.0; };
    D[i] = I("data_obs");
    W[i] = I("signal") + std::max(0.0, I("wtau"));
    double q = I("qcd_abcd");
    if (q < 0) q = I("qcd");
    B[i] = std::max(0.0, I("z")) + std::max(0.0, I("ztau")) + std::max(0.0, I("ttbar")) + std::max(0.0, q);
    if (D[i] <= 0 || W[i] <= 0) { printf("[WARN] %s %s: empty region %s\n", inFile, charge, dir.Data()); return nullptr; }
  }
  const int i0 = resp->GetXaxis()->FindBin(0.0);
  TGraph *g = new TGraph();
  for (int d = 1; d <= resp->GetNbinsX(); ++d)
  {
    double Wd[pOSkim::kNY], sxy = 0, sxx = 0;
    for (int i = 0; i < pOSkim::kNY; ++i)
    {
      Wd[i] = W[i] * resp->GetBinContent(d, i + 1) / resp->GetBinContent(i0, i + 1);
      sxy += (D[i] - B[i]) * Wd[i] / D[i];
      sxx += Wd[i] * Wd[i] / D[i];
    }
    const double s = sxy / sxx;
    double c2 = 0;
    for (int i = 0; i < pOSkim::kNY; ++i) c2 += std::pow(D[i] - B[i] - s * Wd[i], 2) / D[i];
    g->SetPoint(g->GetN(), resp->GetXaxis()->GetBinCenter(d), c2);
  }
  ok = true;
  return g;
}

// Shift to the minimum; best value and the outermost Delta chi2 = 1 crossings
// (err < 0 = a crossing lies outside the scanned range).
TGraph *GbcNorm(TGraph *g, double &best, double &err, double &chi2min)
{
  int im = 0;
  for (int k = 1; k < g->GetN(); ++k) if (g->GetY()[k] < g->GetY()[im]) im = k;
  chi2min = g->GetY()[im];
  best = g->GetX()[im];
  TGraph *o = new TGraph();
  for (int k = 0; k < g->GetN(); ++k) o->SetPoint(k, g->GetX()[k], g->GetY()[k] - chi2min);
  double lo = -99, hi = 99;
  for (int k = o->GetN() - 1; k > im; --k)
    if (o->GetY()[k - 1] <= 1 && o->GetY()[k] > 1)
    { hi = o->GetX()[k - 1] + (1 - o->GetY()[k - 1]) / (o->GetY()[k] - o->GetY()[k - 1]) * (o->GetX()[k] - o->GetX()[k - 1]); break; }
  for (int k = 0; k < im; ++k)
    if (o->GetY()[k] > 1 && o->GetY()[k + 1] <= 1)
    { lo = o->GetX()[k] + (o->GetY()[k] - 1) / (o->GetY()[k] - o->GetY()[k + 1]) * (o->GetX()[k + 1] - o->GetX()[k]); break; }
  err = (lo < -90 || hi > 90) ? -1 : 0.5 * (hi - lo);
  return o;
}

bool GbcRunScan(const std::vector<GbcGen> &gens, TFile *out)
{
  printf("\n================ MODE scan: do the data prefer the MC shifted in rapidity? ================\n");
  TH2D *rp = nullptr, *rm = nullptr;
  for (auto &g : gens) { if (g.tag == "Wp") rp = g.hShift; if (g.tag == "Wm") rm = g.hShift; }
  if (!rp || !rm) { printf("[ERROR] scan needs the W+ and W- gen loops (mode gen in the same run)\n"); return false; }
  const char *inMu = "../plotting/plots/combine_input_W_leppt_mt40.root";
  const char *inEl = "../plotting/plots/Elec/combine_input_W_leppt_mt40.root";
  printf("[INPUT] %s\n[INPUT] %s\n", inMu, inEl);
  printf("[CONFIG] per lab bin: (data - z - ztau - ttbar - qcd_abcd) vs s * (signal + wtau) * R(dy), s free per"
         " flavour and charge, data Poisson errors; R(dy) = gen fiducial response of the MC under test."
         " Meaningful only if the Combine inputs were built from that same MC.\n");

  struct Ch { const char *in, *charge, *label; TH2D *resp; int col, sty; };
  const std::vector<Ch> chans = {
      {inMu, "Wp", "W^{+}#rightarrow#mu#nu", rp, kAzure + 2, 1}, {inMu, "Wm", "W^{#minus}#rightarrow#mu#nu", rm, kRed + 1, 1},
      {inEl, "Wp", "W^{+}#rightarrowe#nu", rp, kAzure - 9, 7}, {inEl, "Wm", "W^{#minus}#rightarrowe#nu", rm, kPink + 6, 7}};
  std::vector<TGraph *> raw;
  for (auto &c : chans)
  {
    bool ok = false;
    TGraph *g = GbcScanOne(c.in, c.charge, c.resp, ok);
    if (!ok) return false;
    raw.push_back(g);
  }
  auto sum = [&](std::vector<int> idx) {
    TGraph *o = new TGraph();
    for (int k = 0; k < raw[0]->GetN(); ++k) { double s = 0; for (int i : idx) s += raw[i]->GetY()[k]; o->SetPoint(k, raw[0]->GetX()[k], s); }
    return o; };

  gSystem->mkdir(kGbcPlotDir, kTRUE);
  PlotStyle ps;
  ps.lm = 0.15; ps.xTitleSize = 0.040;
  ps.headerX = 0.50; ps.headerY = 0.89; ps.headerDy = 0.042; ps.titleSize = 0.032; ps.subSize = 0.027;
  const double yMax = 60;
  TCanvas *c = GbcCanvas("c_gbc_scan", ps);
  TH1D *fr = new TH1D("gbc_scan_frame", "", 10, kGbcShift0, kGbcShift0 + (kGbcNShift - 1) * kGbcShiftStep);
  fr->SetDirectory(nullptr);
  ApplyHistStyle(fr, ps, "extra MC rapidity shift #deltay (+ = proton side)", "#Delta#chi^{2}");
  fr->SetMinimum(0); fr->SetMaximum(yMax);
  fr->Draw("AXIS");
  TLegend *L = new TLegend(0.33, 0.50, 0.80, 0.76);
  L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextFont(42); L->SetTextSize(0.026);
  auto draw = [&](TGraph *g, const char *name, const char *label, int col, int sty, int wid) {
    double best, err, cmin;
    TGraph *o = GbcNorm(g, best, err, cmin);
    o->SetLineColor(col); o->SetLineStyle(sty); o->SetLineWidth(wid); o->Draw("L same");
    o->Write(Form("scan_%s", name));
    printf("[RESULT] scan %-8s best dy %+.3f %s  chi2min %.1f  Delta chi2 at dy=0: %.1f\n", name, best,
           err < 0 ? "(range edge)" : Form("+- %.3f", err), cmin, o->Eval(0.0));
    if (err < 0) L->AddEntry(o, Form("%s: %+.2f (range edge)", label, best), "l");
    else L->AddEntry(o, Form("%s: %+.2f #pm %.2f, #Delta#chi^{2}(0) = %.1f", label, best, err, o->Eval(0.0)), "l");
    return std::make_pair(best, err); };
  out->cd();
  const char *nm[4] = {"mu_Wp", "mu_Wm", "e_Wp", "e_Wm"};
  for (size_t i = 0; i < chans.size(); ++i) draw(raw[i], nm[i], chans[i].label, chans[i].col, chans[i].sty, 2);
  auto mu = draw(sum({0, 1}), "mu", "#mu, W^{+}+W^{#minus}", kBlack, 1, 4);
  draw(sum({0, 1, 2, 3}), "all", "all four", kGray + 1, 2, 3);
  L->Draw();
  GbcVLine(0, 0, 0.45 * yMax, kGray + 2, 3);
  GbcVLine(kGbcShift, 0, 0.45 * yMax, kGreen + 2, 7);
  DrawHeader(ps, "Data vs MC W yields per lab bin",
             "correct MC: #deltay #approx 0", Form("full pO boost: %+.3f (green)", kGbcShift));
  printf("[VERDICT] data shift scan (muons): dy = %+.2f %s -- %s\n", mu.first,
         mu.second < 0 ? "(range edge)" : Form("+- %.2f", mu.second),
         (mu.second > 0 && std::fabs(mu.first) < 2 * mu.second) ? "consistent with no missing shift"
                                                                 : "the data prefer a shifted MC (degenerate with other y-dependent mismodeling)");
  GbcSave(c, "data_shift_scan");
  return true;
}
} // namespace

void gen_boost_check(TString mode = "all", TString wpFile = "", TString wmFile = "", TString dyFile = "",
                     TString dataFile = "")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  if (mode != "gen" && mode != "side" && mode != "scan" && mode != "all")
  { printf("[ERROR] mode must be gen | side | scan | all (got %s)\n", mode.Data()); return; }
  // The gen loop always runs: side needs its oxygen side, scan its shift response.
  const bool doSide = (mode == "side" || mode == "all");
  const bool doScan = (mode == "scan" || mode == "all");

  // inputs: the production in use, unless overridden
  std::vector<GbcGen> gens(3);
  gens[0].tag = "Wp"; gens[0].pdg = 24;
  gens[1].tag = "Wm"; gens[1].pdg = 24;
  gens[2].tag = "Z";  gens[2].pdg = 23;
  gens[0].label = wpFile.Length() ? wpFile.Data() : pOSkim::ResolveMCSample(kWp, "mu").fname;
  gens[1].label = wmFile.Length() ? wmFile.Data() : pOSkim::ResolveMCSample(kWm, "mu").fname;
  gens[2].label = dyFile.Length() ? dyFile.Data() : pOSkim::ResolveMCSample(kDY, "mu").fname;
  const std::string data = dataFile.Length() ? dataFile.Data() : pOSkim::kDefaultDataFile;
  printf("[CONFIG] mode %s; MC %s%s\n", mode.Data(),
         (wpFile.Length() || wmFile.Length() || dyFile.Length()) ? "OVERRIDDEN by arguments" : "from pOSkim::ResolveMCSample",
         doSide ? Form("; data %s", dataFile.Length() ? "OVERRIDDEN" : "from pOSkim::kDefaultDataFile") : "");

  gSystem->mkdir("rootfile", kTRUE);
  gSystem->mkdir(kGbcPlotDir, kTRUE);
  int oxygenSide = 0;
  if (!GbcRunGen(gens, oxygenSide)) { printf("[FATAL] gen loop failed\n"); return; }

  TFile out(kGbcRoot, "RECREATE");
  for (auto &g : gens)
    for (TH1 *h : std::vector<TH1 *>{g.hEp, g.hEm, g.hDb, g.hDs, g.hDa, g.hR, g.hY, g.hDr, g.hShift})
      if (h) h->Write();
  if (doSide && !GbcRunSide(data, gens, oxygenSide)) printf("[ERROR] mode side failed\n");
  if (doScan && !GbcRunScan(gens, &out)) printf("[ERROR] mode scan failed\n");
  out.Close();
  printf("\n[DONE] %s, %s/\n", kGbcRoot, kGbcPlotDir);
}
