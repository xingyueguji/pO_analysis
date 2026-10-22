// SCRATCH (2026-10-20): per charge, a 2-parameter fit of the iso-pass m_T spectrum
//   N_i = r * W_i + O_i + s * Q_i
// W = W MC (signal), O = other MC (DY, tau, ttbar), both absolute (k_s, skim SFs);
// Q = the SAME charge's anti-iso m_T shape from data minus EWK MC (unit area), s = its
// yield. Inputs = the skim's y-binned planes h_iso_mt_<lep><Plus|Minus>_y<i> (no new
// loop). Regions: bwd edge y0, fwd edge y11, both edges, central y2..y9. Fit on m_T
// 0-120 in 10 GeV bins, Poisson likelihood. Also r from the m_T > 40 bins with s
// fixed to the full fit, and the per-bin residuals -- a MIGRATION in m_T shows up as
// a bad fit with an excess at low m_T and a deficit at high m_T.
// Run: root -l -b -q fit_mt_charge.C+
#include "/Users/zhenghuang/pO_analysis/skim/skim_common.h"
#include "/Users/zhenghuang/pO_analysis/skim/mc_norm.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TFile.h"
#include "TH2D.h"

#include <iostream>
#include <memory>
#include <vector>

namespace
{
const char *kSkim = "/Users/zhenghuang/pO_analysis/skim/rootfile/";
const char *kNgen = "/Users/zhenghuang/pO_analysis/skim/rootfile/ngen.root";
const int kNb = 12; // m_T 0-120 in 10 GeV
struct Spec { std::vector<double> d, w, o, q; };

// project the iso-pass / anti-iso windows of the summed y bins onto m_T (10 GeV)
void Proj(TFile *f, const char *lep, const char *chg, const std::vector<int> &ys, double scale, bool isMu,
          std::vector<double> &pass, std::vector<double> &anti)
{
  pass.assign(kNb, 0.0); anti.assign(kNb, 0.0);
  const double pLo = isMu ? 0.0 : -0.3, pHi = isMu ? 0.15 : 0.0, aLo = 0.3, aHi = 1.0;
  for (int y : ys)
  {
    TH2 *h = (TH2 *)f->Get(Form("h_iso_mt_%s%s_y%d", lep, chg, y));
    if (!h) { std::cerr << "[FATAL] missing " << Form("h_iso_mt_%s%s_y%d", lep, chg, y) << " in " << f->GetName() << "\n"; exit(2); }
    for (int ix = 1; ix <= h->GetNbinsX(); ++ix)
    {
      const double xc = h->GetXaxis()->GetBinCenter(ix);
      const bool inP = xc > pLo && xc < pHi, inA = xc > aLo && xc < aHi;
      if (!inP && !inA) continue;
      for (int iy = 1; iy <= h->GetNbinsY(); ++iy)
      {
        const int b = (int)(h->GetYaxis()->GetBinCenter(iy) / 10.0);
        if (b < 0 || b >= kNb) continue;
        (inP ? pass : anti)[b] += scale * h->GetBinContent(ix, iy);
      }
    }
  }
}

double NLL(const Spec &s, double r, double q, int b0)
{
  double v = 0;
  for (int i = b0; i < kNb; ++i)
  {
    const double mu = std::max(1e-9, r * s.w[i] + s.o[i] + q * s.q[i]);
    v += mu - s.d[i] * std::log(mu);
  }
  return v;
}
} // namespace

void fit_mt_charge()
{
  struct Fl { const char *tag, *lep; bool isMu; const char *data, *wp, *wm; std::vector<std::pair<const char *, const char *>> oth; };
  const std::vector<Fl> fls = {
      {"mu", "mu", true, "WToMuNu_pO_PFMet_hist.root", "Wp_mu", "Wm_mu",
       {{"DY", "DYmu"}, {"DYtau", "DYtau"}, {"Wptau", "Wp_tau"}, {"Wmtau", "Wm_tau"}, {"TT", "TT"}}},
      {"e", "ele", false, "WToElecNu_pO_PFMet_hist.root", "Wp_ele", "Wm_ele",
       {{"DY", "DYee"}, {"DYtau", "DYtau"}, {"Wptau", "Wp_tau"}, {"Wmtau", "Wm_tau"}, {"TT", "TT"}}}};
  const std::vector<std::pair<const char *, std::vector<int>>> regs = {
      {"bwd y0", {0}}, {"fwd y11", {11}}, {"edges", {0, 11}}, {"central", {2, 3, 4, 5, 6, 7, 8, 9}}};
  for (const Fl &fl : fls)
  {
    const std::string base = std::string(kSkim) + (fl.isMu ? "WToMuNu_pO_PFMet_" : "WToElecNu_pO_PFMet_");
    TFile *fd = TFile::Open((std::string(kSkim) + fl.data).c_str());
    TFile *fwp = TFile::Open((base + "Wp_hist.root").c_str()), *fwm = TFile::Open((base + "Wm_hist.root").c_str());
    std::vector<std::pair<TFile *, double>> fo;
    for (auto &o : fl.oth) fo.push_back({TFile::Open((base + o.first + "_hist.root").c_str()), pONorm::MCScale(o.second, kNgen)});
    const double kp = pONorm::MCScale(fl.wp, kNgen), km = pONorm::MCScale(fl.wm, kNgen);
    std::cout << "\n================ " << fl.tag << " ================\n"
              << "  region   chg | data(iso,mT<120) | fit: r          s(QCD)      chi2/ndf | r(mT>40, s fixed)  | QCD/data(mT>40) | residuals (data-fit)/sqrt(fit) per 10 GeV, mT 0..120\n";
    for (auto &rg : regs)
      for (const char *chg : {"Minus", "Plus"})
      {
        Spec s;
        std::vector<double> p, a, pw, aw, po, ao;
        Proj(fd, fl.lep, chg, rg.second, 1.0, fl.isMu, s.d, a);
        s.w.assign(kNb, 0); s.o.assign(kNb, 0); s.q.assign(kNb, 0);
        std::vector<double> aE(kNb, 0.0); // anti-iso EWK
        for (auto fw : {std::make_pair(fwp, kp), std::make_pair(fwm, km)})
        {
          Proj(fw.first, fl.lep, chg, rg.second, fw.second, fl.isMu, pw, aw);
          for (int i = 0; i < kNb; ++i) { s.w[i] += pw[i]; aE[i] += aw[i]; }
        }
        for (auto &o : fo)
        {
          Proj(o.first, fl.lep, chg, rg.second, o.second, fl.isMu, po, ao);
          for (int i = 0; i < kNb; ++i) { s.o[i] += po[i]; aE[i] += ao[i]; }
        }
        double qt = 0;
        for (int i = 0; i < kNb; ++i) { s.q[i] = std::max(0.0, a[i] - aE[i]); qt += s.q[i]; }
        for (int i = 0; i < kNb; ++i) s.q[i] /= qt;
        // fit r, s
        auto f2 = [&](const double *x) { return NLL(s, x[0], x[1], 0); };
        std::unique_ptr<ROOT::Math::Minimizer> mn(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
        ROOT::Math::Functor fn(f2, 2);
        mn->SetFunction(fn);
        mn->SetErrorDef(0.5);
        mn->SetLimitedVariable(0, "r", 1.2, 0.05, 0.0, 5.0);
        mn->SetLimitedVariable(1, "s", 100, 10, 0.0, 1e5);
        mn->Minimize();
        mn->Hesse();
        const double r = mn->X()[0], q = mn->X()[1], re = mn->Errors()[0], qe = mn->Errors()[1];
        double chi2 = 0, nd = 0, dt = 0;
        std::string res;
        for (int i = 0; i < kNb; ++i)
        {
          const double mu = r * s.w[i] + s.o[i] + q * s.q[i];
          chi2 += (mu > 0 && s.d[i] > 0) ? 2 * (mu - s.d[i] + s.d[i] * std::log(s.d[i] / mu)) : 2 * mu;
          res += Form(" %+4.1f", (s.d[i] - mu) / std::sqrt(std::max(mu, 1e-9)));
          dt += s.d[i];
        }
        // r from m_T > 40 with s fixed
        auto f1 = [&](const double *x) { return NLL(s, x[0], q, 4); };
        std::unique_ptr<ROOT::Math::Minimizer> m1(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
        ROOT::Math::Functor fn1(f1, 1);
        m1->SetFunction(fn1);
        m1->SetErrorDef(0.5);
        m1->SetLimitedVariable(0, "r", 1.2, 0.05, 0.0, 5.0);
        m1->Minimize();
        m1->Hesse();
        double dHi = 0, qHi = 0;
        for (int i = 4; i < kNb; ++i) { dHi += s.d[i]; qHi += q * s.q[i]; }
        std::cout << Form("  %-8s %-5s | %6.0f           | %.3f +- %.3f  %6.1f +- %5.1f  %5.1f/%d | %.3f +- %.3f     | %.3f          |%s\n",
                          rg.first, chg, dt, r, re, q, qe, chi2, kNb - 2, m1->X()[0], m1->Errors()[0], qHi / dHi, res.c_str());
      }
  }
}
