// correction/ptmt_scan.C
//
// 2D scan of the (lepton pT, m_T) cut pair for the W selection -- which
// combination gives the best S-to-noise -- in the SAME format/conventions as
// the isolation cut studies (isolation_mu_tight.C / isolation_ele.C):
//
//   SIGNAL     = W signal MC (Wp + Wm), absolute k_s = A*sigma*L/N_gen from
//                skim/mc_norm.h, iso-pass, full W selection.
//   BACKGROUND = the ANTI-ISO sideband: DATA passing the full W selection
//                except isolation, with relIso in [0.30,1.0) mu / [0.20,1.0)
//                e (the qcd_abcd.C windows), keeping the ACTUAL MET and m_T,
//                minus the tiny EWK leakage (absolute MC, ~0.5-0.7%). The
//                sideband provides the QCD SHAPE in (pT, m_T); its
//                normalization is set to the total iso-pass QCD = data - EWK
//                over the full plane (the ~19% mu / ~50% e contamination
//                measured 2026-07-29). On top of the QCD, the non-W EWK
//                backgrounds (DY, DYtau, W/DY->tau) above the cuts are added
//                to the denominator.
//
// WHY anti-iso (2026-07-30, replaces the earlier MET<5 "low-MET proxy"):
// the background model for a scan must be selected with a variable that is
// independent of the SCANNED variables. m_T is built from MET, so a
// MET-conditioned proxy sculpts the m_T axis mechanically (MET<5 caps the
// proxy at m_T <~ 2*sqrt(pT*MET) ~ 30 GeV -> fake total rejection above ~35).
// Isolation is ~independent of (pT, m_T), so the sideband QCD keeps its true
// m_T reach. Residual caveat: relIso carries pT in its denominator, so the
// sideband (pT, m_T) shape has a mild iso-pT correlation bias -- bounded by
// the qcd_abcd.C anti-iso slice-stability checks (antiiso_shape_pt_*: flat
// for mu, ~15% tilt for e) -- a small, MEASURED bias instead of the
// structural m_T cap of the MET<5 convention.
//
// Isolation is NOT scanned here -- it stays at the nominal cut baked into the
// skim histograms (mu 0.15 / e 0.095).
//
// Inputs : ../skim/rootfile/<base>[_<sample>]_hist.root with the joint 2Ds
//          h_pt_mt_{mu,ele}{Plus,Minus}           (iso-pass, all events)
//          h_pt_mt_antiiso_{mu,ele}{Plus,Minus}   (anti-iso sideband)
// Output : ./plots/ptmt_scan_{mu,ele}/ + ./rootfile/ptmt_scan_{mu,ele}.root
//
// Run from correction/:  root -l -q 'ptmt_scan.C+'        (muon)
//                        root -l -q 'ptmt_scan.C+(true)'  (electron)

#include "TFile.h"
#include "TH2.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../skim/mc_norm.h" // pONorm::MCScale -> per-sample k_s

namespace
{
// ---- scan grid: cut values land on histogram bin edges (pT 1 GeV, mT 2.5) --
const double kPtLo = 25.0, kPtHi = 45.0, kPtStep = 1.0;
const double kMtLo = 0.0,  kMtHi = 80.0, kMtStep = 2.5;
const int kNPt = (int)((kPtHi - kPtLo) / kPtStep) + 1; // 21
const int kNMt = (int)((kMtHi - kMtLo) / kMtStep) + 1; // 33

// Integral of h over pT >= ptc AND mT >= mtc, overflow included on both axes.
double CumAbove(TH2 *h, double ptc, double mtc)
{
  if (!h) return 0.0;
  const int bx = h->GetXaxis()->FindBin(ptc + 1e-6);
  const int by = h->GetYaxis()->FindBin(mtc + 1e-6);
  return h->Integral(bx, h->GetNbinsX() + 1, by, h->GetNbinsY() + 1);
}

// Clone (detached) the Plus+Minus pair of `stem` from f, scaled by k.
TH2D *SumCharges(TFile *f, const std::string &stem, double k, const char *out)
{
  TH2D *hp = (TH2D *)f->Get((stem + "Plus").c_str());
  TH2D *hm = (TH2D *)f->Get((stem + "Minus").c_str());
  if (!hp || !hm)
  {
    std::cerr << "[ERR] missing " << stem << "Plus/Minus in " << f->GetName()
              << " (re-skim with the updated skim.C first)\n";
    return nullptr;
  }
  TH2D *s = (TH2D *)hp->Clone(out); // clone before mutating (repo rule)
  s->SetDirectory(nullptr);
  s->Add(hm);
  s->Scale(k);
  return s;
}

void ClampNonNegative(TH2 *h)
{
  if (!h) return;
  for (int i = 0; i <= h->GetNbinsX() + 1; ++i)
    for (int j = 0; j <= h->GetNbinsY() + 1; ++j)
      if (h->GetBinContent(i, j) < 0) h->SetBinContent(i, j, 0.0);
}

// MC sample spec and ROC sweep definition. File-scope (not function-local):
// local structs have no linkage and cling refuses std::vector<local-type>.
struct Spec { std::string label, suffix; bool isSig; };
struct RocDef { const char *name; bool sweepMt; double fixed; int color; };

// AUC of a (effB, effS) curve via trapezoids, points sorted in effB.
double AucOf(std::vector<std::pair<double, double>> pts)
{
  std::sort(pts.begin(), pts.end());
  double a = 0.0;
  for (size_t i = 1; i < pts.size(); ++i)
    a += 0.5 * (pts[i].first - pts[i - 1].first) * (pts[i].second + pts[i - 1].second);
  // extend to (0,0) and (1,1) endpoints if the sweep does not reach them
  if (!pts.empty())
  {
    a += 0.5 * pts.front().first * pts.front().second;                      // from (0,0)
    a += 0.5 * (1.0 - pts.back().first) * (1.0 + pts.back().second);       // to (1,1)
  }
  return a;
}
} // namespace

void ptmt_scan(bool isElec = false)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const std::string lep      = isElec ? "ele" : "mu";
  const std::string lepSym   = isElec ? "e" : "#mu";
  const std::string fileBase = isElec ? "WToElecNu_pO_PFMet" : "WToMuNu_pO_PFMet";
  const std::string dir      = "../skim/rootfile/";

  // tentative working point to benchmark against
  const double kPtTent = 25.0, kMtTent = 40.0;

  std::cout << "\n[NOTE] background = anti-iso sideband QCD (relIso "
            << (isElec ? "[0.20,1.0)" : "[0.30,1.0)")
            << ", actual MET/m_T -- no m_T sculpting). Residual caveat: mild"
            << " iso-pT shape correlation, bounded by the qcd_abcd.C"
            << " antiiso_shape_pt_* slice checks.\n";

  // ---- open files ----------------------------------------------------------
  TFile *fData = TFile::Open((dir + fileBase + "_hist.root").c_str(), "READ");
  if (!fData || fData->IsZombie()) { std::cerr << "[ERR] no data file\n"; return; }

  const std::string dyLab = isElec ? "DYee" : "DYmu";
  const std::vector<Spec> spec = {
      {"Wp_" + lep, "Wp",    true},
      {"Wm_" + lep, "Wm",    true},
      {dyLab,       "DY",    false},
      {"DYtau",     "DYtau", false},
      {"Wp_tau",    "Wptau", false},
      {"Wm_tau",    "Wmtau", false},
  };
  std::vector<TFile *> fMC;
  for (const Spec &s : spec)
  {
    TFile *f = TFile::Open((dir + fileBase + "_" + s.suffix + "_hist.root").c_str(), "READ");
    if (!f || f->IsZombie()) { std::cerr << "[ERR] no MC file " << s.suffix << "\n"; return; }
    fMC.push_back(f);
  }

  const std::string stem        = "h_pt_mt_" + lep;
  const std::string stemAntiIso = "h_pt_mt_antiiso_" + lep;

  // ---- build the ingredients ----------------------------------------------
  // S: absolute W signal MC. EWKbkg: absolute non-W EWK. EWKall: for the QCD
  // total. hSide*: the ANTI-ISO SIDEBAND (data - EWK leakage) = the QCD shape.
  // (Until 2026-07-30 this was a MET<5 "proxy"; see the header for why that
  //  convention was removed. Any remaining `Side` naming means sideband.)
  TH2D *hS = nullptr, *hEwkBkg = nullptr, *hEwkAll = nullptr, *hSideEwk = nullptr;
  for (size_t i = 0; i < spec.size(); ++i)
  {
    const double k = pONorm::MCScale(spec[i].label);
    TH2D *h  = SumCharges(fMC[i], stem,        k, Form("mc_%zu", i));
    TH2D *hl = SumCharges(fMC[i], stemAntiIso, k, Form("mcanti_%zu", i));
    if (!h || !hl) return;
    if (!hEwkAll)   { hEwkAll = (TH2D *)h->Clone("ewk_all");     hEwkAll->SetDirectory(nullptr); }
    else            hEwkAll->Add(h);
    if (!hSideEwk) { hSideEwk = (TH2D *)hl->Clone("sideband_ewk"); hSideEwk->SetDirectory(nullptr); }
    else            hSideEwk->Add(hl);
    if (spec[i].isSig) { if (!hS) { hS = h; h = nullptr; } else { hS->Add(h); } }
    else               { if (!hEwkBkg) { hEwkBkg = h; h = nullptr; } else { hEwkBkg->Add(h); } }
    delete h;
    delete hl;
  }

  TH2D *hData      = SumCharges(fData, stem,        1.0, "data_all");
  TH2D *hSideData = SumCharges(fData, stemAntiIso, 1.0, "sideband_data");
  if (!hData || !hSideData || !hS) return;

  // QCD shape = anti-iso data - anti-iso EWK leakage, clamped >= 0.
  // (The `Proxy` in these names is historical -- it was the MET<5 low-MET proxy
  //  until 2026-07-30. These objects now hold the ANTI-ISO SIDEBAND.)
  TH2D *hSide = (TH2D *)hSideData->Clone("sideband_qcd");
  hSide->SetDirectory(nullptr);
  hSide->Add(hSideEwk, -1.0);
  ClampNonNegative(hSide);

  // absolute iso-pass QCD normalization = data - EWK over the FULL plane
  const double nData   = CumAbove(hData, 0, 0);
  const double nEwkAll = CumAbove(hEwkAll, 0, 0);
  const double qcdTot  = std::max(0.0, nData - nEwkAll);
  const double sideTot    = CumAbove(hSide, 0, 0);   // = sideband total
  const double sideRawTot = CumAbove(hSideData, 0, 0);

  std::cout << "[INFO] " << lep << ": data=" << nData << "  EWK(abs MC)=" << nEwkAll
            << "  -> iso-pass QCD total=" << qcdTot << "\n"
            << "[INFO] anti-iso sideband: " << sideRawTot << " data events, "
            << sideTot << " after EWK subtraction -> sideband shape scaled to "
            << qcdTot << "\n";

  // ---- the scan ------------------------------------------------------------
  // maps binned so each cut value is a bin center
  TH2D *mSton = new TH2D("m_ston", ";p_{T} cut [GeV];m_{T} cut [GeV];S/#sqrt{S+B}",
                         kNPt, kPtLo - 0.5 * kPtStep, kPtHi + 0.5 * kPtStep,
                         kNMt, kMtLo - 0.5 * kMtStep, kMtHi + 0.5 * kMtStep);
  TH2D *mSoB  = (TH2D *)mSton->Clone("m_sob");  mSoB->SetDirectory(nullptr);  mSoB->SetTitle(";p_{T} cut [GeV];m_{T} cut [GeV];S/B");
  TH2D *mEffS = (TH2D *)mSton->Clone("m_effS"); mEffS->SetDirectory(nullptr); mEffS->SetTitle(";p_{T} cut [GeV];m_{T} cut [GeV];#varepsilon_{S}");
  mSton->SetDirectory(nullptr);

  const double S0 = CumAbove(hS, kPtLo, 0.0); // baseline: current selection pT>25, no mT cut

  double bestSton = -1, bestPt = kPtLo, bestMt = 0;
  auto evalPoint = [&](double ptc, double mtc, double &S, double &B, double &ston, double &sob)
  {
    S = CumAbove(hS, ptc, mtc);
    const double bq = qcdTot * (sideTot > 0 ? CumAbove(hSide, ptc, mtc) / sideTot : 0.0);
    const double be = hEwkBkg ? CumAbove(hEwkBkg, ptc, mtc) : 0.0;
    B = bq + be;
    ston = (S + B > 0) ? S / std::sqrt(S + B) : 0.0;
    sob  = (B > 0) ? S / B : 0.0;
  };

  for (int ip = 0; ip < kNPt; ++ip)
    for (int im = 0; im < kNMt; ++im)
    {
      const double ptc = kPtLo + ip * kPtStep;
      const double mtc = kMtLo + im * kMtStep;
      double S, B, ston, sob;
      evalPoint(ptc, mtc, S, B, ston, sob);
      mSton->SetBinContent(ip + 1, im + 1, ston);
      mSoB->SetBinContent(ip + 1, im + 1, sob);
      mEffS->SetBinContent(ip + 1, im + 1, S0 > 0 ? S / S0 : 0.0);
      if (ston > bestSton) { bestSton = ston; bestPt = ptc; bestMt = mtc; }
    }

  // ---- ROC curves + AUC (normalization-free view) --------------------------
  // sweep the mT cut at fixed pT cuts, and the pT cut at mT=0; efficiencies
  // relative to the (25, 0) baseline. Background eff uses QCD only (the anti-iso
  // sideband), excluding the non-W EWK that the FOM's B does include.
  // (named qcdBase, not "B0": B0 is a termios.h baud-rate macro on macOS)
  const double qcdBase = qcdTot * (sideTot > 0 ? CumAbove(hSide, kPtLo, 0.0) / sideTot : 0.0);
  const std::vector<RocDef> rocs = {
      {"m_{T} sweep, p_{T}>25", true, 25.0, 600 + 1}, // kBlue+1
      {"m_{T} sweep, p_{T}>28", true, 28.0, 632 + 1}, // kRed+1
      {"m_{T} sweep, p_{T}>30", true, 30.0, 416 + 2}, // kGreen+2
      {"p_{T} sweep, no m_{T} cut", false, 0.0, 880 + 1}, // kViolet+1
  };
  std::vector<TGraph *> gRocs;
  std::vector<double> aucs;
  for (const RocDef &rd : rocs)
  {
    TGraph *g = new TGraph();
    std::vector<std::pair<double, double>> pts;
    const int n = rd.sweepMt ? kNMt : kNPt;
    for (int i = 0; i < n; ++i)
    {
      const double ptc = rd.sweepMt ? rd.fixed : (kPtLo + i * kPtStep);
      const double mtc = rd.sweepMt ? (kMtLo + i * kMtStep) : rd.fixed;
      const double eS = S0 > 0 ? CumAbove(hS, ptc, mtc) / S0 : 0.0;
      const double eB = qcdBase > 0 ? qcdTot * (sideTot > 0 ? CumAbove(hSide, ptc, mtc) / sideTot : 0.0) / qcdBase : 0.0;
      g->SetPoint(g->GetN(), eB, eS);
      pts.push_back({eB, eS});
    }
    g->SetLineColor(rd.color);
    g->SetMarkerColor(rd.color);
    g->SetLineWidth(2);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.5);
    gRocs.push_back(g);
    aucs.push_back(AucOf(pts));
  }

  // ---- report --------------------------------------------------------------
  auto printPoint = [&](const char *tag, double ptc, double mtc)
  {
    double S, B, ston, sob;
    evalPoint(ptc, mtc, S, B, ston, sob);
    std::cout << Form("  %-28s pT>%4.1f mT>%5.1f : S=%7.1f  B=%7.1f  S/B=%6.2f"
                      "  S/sqrt(S+B)=%6.2f  effS=%5.1f%%\n",
                      tag, ptc, mtc, S, B, sob, ston, S0 > 0 ? 100.0 * S / S0 : 0.0);
  };
  std::cout << "\n=== pT x mT scan  [" << lep << "]  (B = anti-iso-sideband QCD + non-W EWK) ===\n";
  printPoint("baseline (current)", kPtLo, 0.0);
  printPoint("tentative", kPtTent, kMtTent);
  printPoint("OPTIMUM (max S/sqrt(S+B))", bestPt, bestMt);
  for (size_t i = 0; i < rocs.size(); ++i)
    std::cout << Form("  AUC  %-24s = %.4f\n", rocs[i].name, aucs[i]);

  // ---- draw ----------------------------------------------------------------
  const std::string outDir = "./plots/ptmt_scan_" + lep;
  gSystem->mkdir(outDir.c_str(), kTRUE);
  gSystem->mkdir("./rootfile", kTRUE);

  auto drawMap = [&](TH2D *m, const char *ztit, const std::string &outNoExt)
  {
    TCanvas *c = new TCanvas(Form("c_%s", m->GetName()), "", 850, 750);
    c->SetRightMargin(0.16);
    c->SetLeftMargin(0.12);
    m->GetZaxis()->SetTitle(ztit);
    m->Draw("COLZ");
    TMarker *opt = new TMarker(bestPt, bestMt, 29);  // star
    opt->SetMarkerColor(1); opt->SetMarkerSize(2.6); opt->Draw();
    TMarker *tnt = new TMarker(kPtTent, kMtTent, 34); // cross
    tnt->SetMarkerColor(2); tnt->SetMarkerSize(2.0); tnt->Draw();
    TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.032);
    t.DrawLatex(0.12, 0.93, Form("%s channel   #star optimum (%.0f, %.1f)   red cross: tentative (%.0f, %.0f)",
                                 lepSym.c_str(), bestPt, bestMt, kPtTent, kMtTent));
    c->SaveAs((outNoExt + ".png").c_str());
    c->SaveAs((outNoExt + ".pdf").c_str());
    delete c;
  };
  drawMap(mSton, "S/#sqrt{S+B}",  outDir + "/fom_ston");
  drawMap(mSoB,  "S/B",           outDir + "/fom_sob");
  drawMap(mEffS, "#varepsilon_{S} (rel. to p_{T}>25, no m_{T} cut)", outDir + "/eff_signal");

  {
    TCanvas *c = new TCanvas("c_roc", "", 800, 750);
    c->SetLeftMargin(0.12);
    TH2D *frame = new TH2D("frame", ";#varepsilon_{B} (anti-iso QCD);#varepsilon_{S}", 10, 0, 1.02, 10, 0, 1.05);
    frame->Draw();
    TLegend *leg = new TLegend(0.45, 0.16, 0.88, 0.42);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.028);
    for (size_t i = 0; i < gRocs.size(); ++i)
    {
      gRocs[i]->Draw("PL same");
      leg->AddEntry(gRocs[i], Form("%s  (AUC %.3f)", rocs[i].name, aucs[i]), "pl");
    }
    leg->Draw();
    TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.032);
    t.DrawLatex(0.14, 0.93, Form("%s channel: cut-sweep ROC (anti-iso sideband QCD)", lepSym.c_str()));
    c->SaveAs((outDir + "/roc.png").c_str());
    c->SaveAs((outDir + "/roc.pdf").c_str());
    delete c;
  }

  {
    // 1D profile: S/sqrt(S+B) vs mT cut at pT>25 and pT>30
    TCanvas *c = new TCanvas("c_prof", "", 800, 650);
    c->SetLeftMargin(0.12);
    TGraph *g25 = new TGraph(), *g30 = new TGraph();
    for (int im = 0; im < kNMt; ++im)
    {
      const double mtc = kMtLo + im * kMtStep;
      double S, B, ston, sob;
      evalPoint(25.0, mtc, S, B, ston, sob); g25->SetPoint(im, mtc, ston);
      evalPoint(30.0, mtc, S, B, ston, sob); g30->SetPoint(im, mtc, ston);
    }
    g25->SetLineColor(600 + 1); g25->SetLineWidth(3);
    g30->SetLineColor(632 + 1); g30->SetLineWidth(3); g30->SetLineStyle(2);
    g25->SetTitle(";m_{T} cut [GeV];S/#sqrt{S+B}");
    g25->Draw("AL");
    g30->Draw("L same");
    TLegend *leg = new TLegend(0.15, 0.18, 0.5, 0.34);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
    leg->AddEntry(g25, "p_{T} > 25 GeV", "l");
    leg->AddEntry(g30, "p_{T} > 30 GeV", "l");
    leg->Draw();
    TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.032);
    t.DrawLatex(0.14, 0.93, Form("%s channel: S/#sqrt{S+B} vs m_{T} cut", lepSym.c_str()));
    c->SaveAs((outDir + "/profile_mt.png").c_str());
    c->SaveAs((outDir + "/profile_mt.pdf").c_str());
    delete c;
  }

  // ---- persist -------------------------------------------------------------
  TFile *fout = TFile::Open(("./rootfile/ptmt_scan_" + lep + ".root").c_str(), "RECREATE");
  mSton->Write(); mSoB->Write(); mEffS->Write();
  hS->Write("h_signal"); hSide->Write("h_qcd_antiiso");
  if (hEwkBkg) hEwkBkg->Write("h_ewk_bkg");
  for (size_t i = 0; i < gRocs.size(); ++i) gRocs[i]->Write(Form("roc_%zu", i));
  fout->Close();

  std::cout << "[DONE] " << lep << ": plots -> " << outDir
            << "/   maps -> rootfile/ptmt_scan_" << lep << ".root\n";

  fData->Close();
  for (TFile *f : fMC) f->Close();
}
