// SANDBOX ONLY (2026-10-20): the prefit signal strength per lab rapidity bin,
// muon vs electron, with the electron MC weighted by the pp EGM SFs (the
// repo's nominal) and by the pO charge-dependent SFs of the Z->ee
// tag-and-probe split by probe charge (this sandbox's re-skim). Reads two
// prefit_r.csv files written by plotting/prefit_r.C (r_pre = (data - DY -
// QCD)/W, QCD as drawn; data Poisson errors only) and draws, per charge,
//   rpre_<C>   : mu | e (EGM SF) | e (pO charge SF)
//   emu_<C>    : e/mu with each electron SF choice
// + a console table and compare_chargesf.csv.
// Run from plotting/: root -l -b -q 'compare_chargesf.C+("<egm csv>", "<pO csv>")'
#include "plotting_helper.C"

#include "TSystem.h"

#include <fstream>
#include <map>
#include <sstream>

namespace
{
struct Row { double x = 0, r = 0, e = 0, rAbcd = 0; };
// key: flavour|charge|ybin (lab only)
std::map<std::string, Row> ReadPrefit(const char *path)
{
  std::map<std::string, Row> out;
  std::ifstream in(path);
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line))
  {
    std::vector<std::string> f;
    std::stringstream ss(line);
    std::string t;
    while (std::getline(ss, t, ',')) f.push_back(t);
    if (f.size() < 14 || f[1] != "lab") continue;
    Row r;
    r.x = std::stod(f[4]); r.r = std::stod(f[10]); r.e = std::stod(f[11]); r.rAbcd = std::stod(f[12]);
    out[f[0] + "|" + f[2] + "|" + f[3]] = r;
  }
  return out;
}
} // namespace

int compare_chargesf(const char *egmCsv, const char *poCsv, const char *outDir = "plots/chargesf")
{
  const auto egm = ReadPrefit(egmCsv), po = ReadPrefit(poCsv);
  if (egm.empty() || po.empty()) { std::cerr << "[FATAL] cannot read " << egmCsv << " / " << poCsv << "\n"; return 2; }
  gSystem->mkdir(outDir, kTRUE);
  std::ofstream csv(Form("%s/compare_chargesf.csv", outDir));
  csv << "charge,ybin,x,r_mu,r_mu_err,r_e_egm,r_e_egm_err,r_e_po,r_e_po_err,emu_egm,emu_egm_err,emu_po,emu_po_err,"
         "r_mu_abcd,r_e_egm_abcd,r_e_po_abcd\n";
  const char *chg[3] = {"Wp", "Wm", "W"}, *chgLab[3] = {"W^{+}#rightarrowl^{+}#nu", "W^{-}#rightarrowl^{-}#bar{#nu}", "W#rightarrowl#nu"};
  for (int c = 0; c < 3; ++c)
  {
    // the muon rows of the two files must agree (the muon input is common)
    TGraphErrors *gMu = new TGraphErrors(), *gEgm = new TGraphErrors(), *gPo = new TGraphErrors();
    TGraphErrors *rEgm = new TGraphErrors(), *rPo = new TGraphErrors();
    double chiEgm = 0, chiPo = 0;
    std::cout << Form("\n[CHARGESF] %s, lab bins: prefit r (QCD as drawn)  |  e/mu  (data stat errors only)\n", chg[c])
              << "   bin      x |   mu              | e (EGM SF)        | e (pO charge SF)  | e/mu EGM         | e/mu pO\n";
    for (int i = 0; i <= 12; ++i)
    {
      const std::string yb = (i < 12) ? std::to_string(i) : "incl";
      const std::string kMu = std::string("mu|") + chg[c] + "|" + yb, kE = std::string("ele|") + chg[c] + "|" + yb;
      if (!egm.count(kMu) || !egm.count(kE) || !po.count(kMu) || !po.count(kE))
      { std::cerr << "[FATAL] missing row " << kE << "\n"; return 2; }
      const Row &m = egm.at(kMu), &m2 = po.at(kMu), &a = egm.at(kE), &b = po.at(kE);
      if (std::fabs(m.r - m2.r) > 1e-9) std::cout << "  [WARN] the muon r differs between the two files in bin " << yb << "\n";
      const double qa = a.r / m.r, qaE = qa * std::hypot(a.e / a.r, m.e / m.r);
      const double qb = b.r / m.r, qbE = qb * std::hypot(b.e / b.r, m.e / m.r);
      std::cout << Form("  %5s %6.2f | %.3f +- %.3f   | %.3f +- %.3f   | %.3f +- %.3f   | %.3f +- %.3f  | %.3f +- %.3f\n",
                        yb.c_str(), m.x, m.r, m.e, a.r, a.e, b.r, b.e, qa, qaE, qb, qbE);
      csv << Form("%s,%s,%.4f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", chg[c], yb.c_str(), m.x, m.r,
                  m.e, a.r, a.e, b.r, b.e, qa, qaE, qb, qbE, m.rAbcd, a.rAbcd, b.rAbcd);
      if (i == 12) continue;
      // the e/mu chi2 against 1 per binning: (e/mu - 1)^2 / err^2 over the 12 bins (stat)
      chiEgm += std::pow((qa - 1.0) / qaE, 2);
      chiPo += std::pow((qb - 1.0) / qbE, 2);
      const int k = gMu->GetN();
      gMu->SetPoint(k, m.x, m.r);   gMu->SetPointError(k, 0.2, m.e);
      gEgm->SetPoint(k, a.x, a.r);  gEgm->SetPointError(k, 0.2, a.e);
      gPo->SetPoint(k, b.x, b.r);   gPo->SetPointError(k, 0.2, b.e);
      rEgm->SetPoint(k, m.x, qa);   rEgm->SetPointError(k, 0.2, qaE);
      rPo->SetPoint(k, m.x, qb);    rPo->SetPointError(k, 0.2, qbE);
    }
    std::cout << Form("  chi2 of e/mu vs 1 over the 12 lab bins (stat): EGM SF %.1f/12, pO charge SF %.1f/12\n", chiEgm, chiPo);
    const Row &ia = egm.at(std::string("ele|") + chg[c] + "|incl"), &ib = po.at(std::string("ele|") + chg[c] + "|incl"),
              &im = egm.at(std::string("mu|") + chg[c] + "|incl");

    auto series = [](TGraphErrors *g, const char *lab, int col, int mk, double sz, double sh)
    {
      OverlaySeries s;
      s.gTot = g; s.label = lab; s.color = col; s.marker = mk; s.markerSize = sz; s.xShift = sh;
      return s;
    };
    auto tuner = [](double lo, double hi, double ref)
    {
      return GraphTuner([lo, hi, ref](TCanvas *cv, TGraphErrors *g)
      {
        g->SetMinimum(lo);
        g->SetMaximum(hi);
        cv->cd();
        TLine *l = new TLine(g->GetXaxis()->GetXmin(), ref, g->GetXaxis()->GetXmax(), ref);
        l->SetLineStyle(2);
        l->SetLineColor(kGray + 2);
        l->Draw();
      });
    };
    PlotStyle ps;
    ps.headerX = 0.55; ps.headerY = 0.90; ps.subSize = 0.032; ps.boxTextSize = 0.028;
    ps.boxX1 = 0.40; ps.boxX2 = 0.93; ps.boxY1 = 0.61; ps.boxY2 = 0.71;
    ps.legX1 = 0.50; ps.legY1 = 0.14; ps.legX2 = 0.93; ps.legY2 = 0.28;
    {
      std::vector<OverlaySeries> s = {series(gMu, "#mu", kAzure + 2, 20, 1.3, -0.5),
                                      series(gEgm, "e, pp EGM SF", kRed + 1, 25, 1.5, 0.0),
                                      series(gPo, "e, pO charge-dependent SF", kRed + 1, 21, 1.5, 0.5)};
      SaveNiceGraph_Overlay(s, Form("%s/rpre_%s", outDir, chg[c]), "y_{lab} = #minus#eta_{lab}",
                            "prefit r = (data #minus DY #minus QCD) / W", "", chgLab[c], "lepton p_{T}, m_{T} > 40 GeV",
                            {Form("inclusive: #mu %.3f, e %.3f (EGM) / %.3f (pO)", im.r, ia.r, ib.r),
                             "data stat. errors only"},
                            ps, tuner(0.2, 2.8, 1.0));
    }
    {
      std::vector<OverlaySeries> s = {series(rEgm, "e/#mu, pp EGM SF", kRed + 1, 25, 1.5, -0.25),
                                      series(rPo, "e/#mu, pO charge-dependent SF", kRed + 1, 21, 1.5, 0.25)};
      PlotStyle p2 = ps;
      p2.legX1 = 0.40; p2.legX2 = 0.85; p2.legY2 = 0.23;
      SaveNiceGraph_Overlay(s, Form("%s/emu_%s", outDir, chg[c]), "y_{lab} = #minus#eta_{lab}", "prefit r_{e} / r_{#mu}", "",
                            chgLab[c], "lepton p_{T}, m_{T} > 40 GeV",
                            {Form("inclusive e/#mu: %.3f (EGM), %.3f (pO)", ia.r / im.r, ib.r / im.r),
                             Form("#chi^{2}(e/#mu = 1), stat: %.1f / %.1f (12 bins)", chiEgm, chiPo)},
                            p2, tuner(0.2, 2.2, 1.0));
    }
  }
  csv.close();
  std::cout << "[OUT] " << outDir << "/{rpre,emu}_{Wp,Wm,W}.{png,pdf}, compare_chargesf.csv\n";
  return 0;
}
