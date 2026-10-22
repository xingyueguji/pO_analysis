// plotting/combine_input_nomt.C
//
// The Combine input of the NO-m_T lepton-pT fit = the nominal since 2026-10-21
// (user decision: drop the m_T > 40 cut because of the electron PF-MET issue;
// leppt_mt40 stays for later). Multijet model = scheme C of
// correction/qcd_nomt_study.C (user confirmation 2026-10-21; record
// docs/iterations/2026-10-21_noreco-nomt-study.md):
//
//   * the signal region is iso-pass, pT 24-100, NO m_T cut, fitted per
//     (charge, rapidity bin) as TWO pT channels split at m_T = 30
//     (skim h_leppt_mt{lo,hi}_*) that share r -- every event once; the m_T < 30
//     channel is the ABCD region B, so the multijet level is measured inside
//     the fit;
//   * multijet per bin and channel = T_y x f(pT) x the group's anti-iso pT
//     shape x the bin's anti-iso count (correction/qcd_abcd_nomt.h runNoMt,
//     nomt_qcd_*), the level T_y floating PER RAPIDITY BIN in the fit
//     (since 2026-10-21b; it was one level per group before);
//   * the anti-iso anchors per bin, D_y (m_T < 30) and H_y (m_T >= 30), as
//     1-bin control regions.
//
// Writes plots[/Elec]/combine_input_W_leppt.root:
//   <C>_<B>_y<i>_lo, <C>_<B>_y<i>_hi   the SR channels: data_obs signal z ztau wtau
//                                      ttbar qcd + <proc>_<syst>Up/Down (signal, z,
//                                      ztau, wtau, ttbar; electron also
//                                      qcd_qcdFF_ele_<C>Up/Down)
//   <C>_<B>_y<i>                       lo + hi (display, the extractor's prefit S,
//                                      xsec_fiducial_diag) -- no systematics
//   Wp_incl, Wm_incl, W_incl           the lab sums (display)
//   <C>_<B>_y<i>_CRD, <C>_<B>_y<i>_CRH the anchors: data_obs, qcd (the QCD-only
//                                      count), ewk (frozen EWK MC), 1 bin each
// + the sidecars <input minus .root>_systs.txt (the shape systematics, as
//   plotting/mtandmet.C writes them) and _mtsplit.txt (the rapidity groups of
//   the shape and the transfer term, and the measured residual kappas of the
//   m_T >= 30 multijet: '#! kappaWin <C> <k>' = the anti-iso window variation,
//   coherent over the bins of a charge, and '#! kappaXfer <B> <group> <C> <k>' =
//   the m_T < 30 -> >= 30 transfer, measured on all 12 bins (the same k for the
//   three groups since 2026-10-22), one uncorrelated nuisance per group --
//   the fork's QCD_MODE=mtsplit reads both sidecars and refuses without them).
// Systematics: the LHE families, the lepton SF(s) of the flavour (electron:
// eleSF + eleReco since 2026-10-21) and the MC MET resolution of the flavour
// (metResMu / metResEle: on these templates it moves W across m_T = 30, the
// only place MET enters this fit), plus the electron fake-factor slope.
//
// Run through plotting/run_combine_inputs.sh (after correction/run_qcd_abcd.sh).

#include "TFile.h"
#include "TH1.h"
#include "TStyle.h"
#include "TSystem.h"
#include "plotting_helper.C"
#include "../skim/mc_norm.h"
#include "../skim/lhe_index.h"
#include "../skim/muon_sf.h"
#include "../skim/electron_sf.h"
#include "../skim/met_res.h"
#include "../correction/abcd_ygroups.h"

#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace nomtin
{
const int kNY = 12;
const char *kBin[2] = {"lab", "fb"};
const char *kPart[2] = {"lo", "hi"};
const char *kPartHist[2] = {"mtlo", "mthi"};

struct Region
{
  std::string dir;
  TH1D *data = nullptr, *sig = nullptr, *z = nullptr, *ztau = nullptr, *wtau = nullptr, *ttbar = nullptr, *qcd = nullptr;
  std::map<std::string, TH1D *> syst; // "<proc>_<syst><Up|Down>"
};

TH1D *Detach(TH1 *h, const std::string &nm)
{
  if (!h) return nullptr;
  TH1D *c = (TH1D *)h->Clone(nm.c_str());
  c->SetDirectory(nullptr);
  return c;
}

// k-scaled copy of one MC histogram (nullptr when absent)
TH1D *GetScaled(TFile *f, double k, const std::string &name, const std::string &out)
{
  TH1 *h = (TH1 *)f->Get(name.c_str());
  if (!h) return nullptr;
  TH1D *c = Detach(h, out);
  c->Scale(k);
  return c;
}

TH1D *Sum(TH1D *a, TH1D *b, const std::string &nm)
{
  if (!a && !b) return nullptr;
  TH1D *r = Detach(a ? a : b, nm);
  if (a && b) r->Add(b);
  return r;
}

void AddTo(TH1D *&acc, const TH1D *h, const std::string &nm)
{
  if (!h) return;
  if (!acc) acc = Detach((TH1 *)h, nm);
  else acc->Add(h);
}
} // namespace nomtin

void combine_input_nomt(bool isElec = false)
{
  using namespace nomtin;
  const std::string base = isElec ? "WToElecNu_pO_PFMet" : "WToMuNu_pO_PFMet";
  const std::string lep = isElec ? "ele" : "mu";
  const std::string outBase = isElec ? "./plots/Elec" : "./plots";
  const std::string inName = "combine_input_W_leppt";
  const std::string outFile = outBase + "/" + inName + ".root";
  gSystem->mkdir((outBase + "/leppt").c_str(), kTRUE);

  // ---- inputs ----
  const std::string dir = "../skim/rootfile/";
  TFile *fD = TFile::Open((dir + base + "_hist.root").c_str(), "READ");
  struct MC { std::string suffix, label; TFile *f = nullptr; double k = 1; };
  std::vector<MC> mc = {{"Wp", "Wp_" + lep}, {"Wm", "Wm_" + lep}, {"DY", isElec ? "DYee" : "DYmu"},
                        {"DYtau", "DYtau"}, {"Wptau", "Wp_tau"}, {"Wmtau", "Wm_tau"}, {"TT", "TT"}};
  if (!fD || fD->IsZombie()) { std::cerr << "[ERROR] cannot open the data skim\n"; return; }
  for (MC &m : mc)
  {
    m.f = TFile::Open((dir + base + "_" + m.suffix + "_hist.root").c_str(), "READ");
    if (!m.f || m.f->IsZombie()) { std::cerr << "[ERROR] cannot open the " << m.suffix << " skim\n"; return; }
    m.k = pONorm::MCScale(m.label);
  }
  TFile *fQ = TFile::Open(("../correction/rootfile/qcd_abcd_" + lep + ".root").c_str(), "READ");
  if (!fQ || fQ->IsZombie()) { std::cerr << "[ERROR] cannot open qcd_abcd_" << lep << ".root (run correction/run_qcd_abcd.sh)\n"; return; }

  // the shape systematics of this flavour's MC templates
  std::vector<std::string> systs(pOLhe::kLheSystNames, pOLhe::kLheSystNames + pOLhe::kNLheSysts);
  if (isElec) systs.insert(systs.end(), pOSF::kEleSFSystNames, pOSF::kEleSFSystNames + pOSF::kNEleSFSysts);
  else        systs.insert(systs.end(), pOSF::kMuonSFSystNames, pOSF::kMuonSFSystNames + pOSF::kNMuonSFSysts);
  systs.push_back(pOMetRes::SystName(isElec));
  const std::string ffName[2] = {"qcdFF_ele_Wp", "qcdFF_ele_Wm"};

  // ---- the SR channels ----
  std::vector<Region> regs;        // the 96 SR channels
  std::vector<Region> sums;        // the 48 lo + hi sums
  bool missing = false;
  std::map<std::string, int> systComplete; // syst -> number of SR channels carrying it fully
  int ffComplete[2] = {0, 0};
  for (int ic = 0; ic < 2; ++ic)
  {
    const std::string C = ic ? "Wm" : "Wp";
    const std::string chg = ic ? "Minus" : "Plus";
    for (int ib = 0; ib < 2; ++ib)
      for (int iy = 0; iy < kNY; ++iy)
      {
        Region sum;
        sum.dir = Form("%s_%s_y%d", C.c_str(), kBin[ib], iy);
        for (int p = 0; p < 2; ++p)
        {
          Region r;
          r.dir = sum.dir + "_" + kPart[p];
          const std::string hn = Form("h_leppt_%s_%s_y%d%s", kPartHist[p], C.c_str(), iy, ib ? "_FB" : "");
          r.data = Detach((TH1 *)fD->Get(hn.c_str()), r.dir + "_data");
          TH1D *s[7];
          for (int k = 0; k < 7; ++k) s[k] = GetScaled(mc[k].f, mc[k].k, hn, r.dir + "_mc" + mc[k].suffix);
          r.sig = Sum(s[0], s[1], r.dir + "_signal");
          r.z = s[2]; r.ztau = s[3];
          r.wtau = Sum(s[4], s[5], r.dir + "_wtau");
          r.ttbar = s[6];
          const std::string qn = Form("nomt_qcd_%s%s_%s_y%d_%s", lep.c_str(), chg.c_str(), kBin[ib], iy, kPart[p]);
          r.qcd = Detach((TH1 *)fQ->Get(qn.c_str()), r.dir + "_qcd");
          if (!r.data || !r.sig || !r.z || !r.ztau || !r.wtau || !r.ttbar || !r.qcd)
          {
            std::cerr << "[ERROR] " << r.dir << ": missing " << hn << " (skim before 2026-10-21?) or " << qn
                      << " (run correction/run_qcd_abcd.sh)\n";
            missing = true;
            continue;
          }
          // the MC shape systematics: k-scaled and combined exactly like the nominals
          for (const std::string &sy : systs)
          {
            std::map<std::string, TH1D *> got;
            bool full = true;
            for (const char *d : {"Up", "Down"})
            {
              const std::string vn = hn + "_" + sy + d;
              TH1D *v[7];
              for (int k = 0; k < 7; ++k) { v[k] = GetScaled(mc[k].f, mc[k].k, vn, r.dir + "_" + sy + d + mc[k].suffix); full &= (v[k] != nullptr); }
              if (!full) break;
              got["signal_" + sy + d] = Sum(v[0], v[1], r.dir + "_signal_" + sy + d);
              got["z_" + sy + d] = v[2];
              got["ztau_" + sy + d] = v[3];
              got["wtau_" + sy + d] = Sum(v[4], v[5], r.dir + "_wtau_" + sy + d);
              got["ttbar_" + sy + d] = v[6];
            }
            if (!full) continue;
            r.syst.insert(got.begin(), got.end());
            ++systComplete[sy];
          }
          if (isElec) // the fake-factor slope twins of the multijet template
          {
            TH1D *up = Detach((TH1 *)fQ->Get((qn + "_ffUp").c_str()), r.dir + "_qffUp");
            TH1D *dn = Detach((TH1 *)fQ->Get((qn + "_ffDown").c_str()), r.dir + "_qffDown");
            if (up && dn)
            {
              r.syst["qcd_" + ffName[ic] + "Up"] = up;
              r.syst["qcd_" + ffName[ic] + "Down"] = dn;
              ++ffComplete[ic];
            }
          }
          AddTo(sum.data, r.data, sum.dir + "_data");
          AddTo(sum.sig, r.sig, sum.dir + "_signal");
          AddTo(sum.z, r.z, sum.dir + "_z");
          AddTo(sum.ztau, r.ztau, sum.dir + "_ztau");
          AddTo(sum.wtau, r.wtau, sum.dir + "_wtau");
          AddTo(sum.ttbar, r.ttbar, sum.dir + "_ttbar");
          AddTo(sum.qcd, r.qcd, sum.dir + "_qcd");
          regs.push_back(r);
        }
        // consistency: lo + hi data = the no-m_T template of the same bin
        TH1 *hAll = (TH1 *)fD->Get(Form("h_leppt_%s_y%d%s", C.c_str(), iy, ib ? "_FB" : ""));
        if (hAll && sum.data && std::fabs(hAll->Integral(1, 50) - sum.data->Integral(1, 50)) > 1e-6)
          std::cerr << "[WARN] " << sum.dir << ": data lo + hi = " << sum.data->Integral(1, 50) << " != h_leppt "
                    << hAll->Integral(1, 50) << "\n";
        sums.push_back(sum);
      }
  }
  if (missing) { std::cerr << "[ERROR] incomplete inputs -> " << outFile << " NOT written\n"; return; }

  // ---- the anchors per bin + the residual kappas ----
  struct CR { std::string dir; double data = 0, qcd = 0, ewk = 0; };
  std::vector<CR> crs;
  auto countBin = [](TH1 *h, const char *lab) -> double {
    for (int b = 1; b <= h->GetNbinsX(); ++b)
      if (std::string(h->GetXaxis()->GetBinLabel(b)) == lab) return h->GetBinContent(b);
    std::cerr << "[ERROR] label " << lab << " not in " << h->GetName() << "\n";
    return NAN;
  };
  double kWin[2] = {1.0, 1.0}, kXfer[2][2][pOAbcdY::kNGroups] = {};
  for (int ic = 0; ic < 2; ++ic)
  {
    const std::string C = ic ? "Wm" : "Wp", chg = ic ? "Minus" : "Plus";
    TH1 *hk = (TH1 *)fQ->Get(Form("nomt_kappa_%s%s", lep.c_str(), chg.c_str()));
    if (!hk) { std::cerr << "[ERROR] nomt_kappa_" << lep << chg << " missing (run correction/run_qcd_abcd.sh)\n"; missing = true; continue; }
    kWin[ic] = countBin(hk, "kappaWin");
    for (int ib = 0; ib < 2; ++ib)
      for (int g = 0; g < pOAbcdY::kNGroups; ++g)
        kXfer[ic][ib][g] = countBin(hk, Form("kappaXfer_%s_%s", kBin[ib], pOAbcdY::kGroupName[g]));
    if (std::isnan(kWin[ic])) missing = true;
    for (int ib = 0; ib < 2; ++ib)
      for (int iy = 0; iy < kNY; ++iy)
      {
        TH1 *hc = (TH1 *)fQ->Get(Form("nomt_counts_%s%s_%s_y%d", lep.c_str(), chg.c_str(), kBin[ib], iy));
        if (!hc) { std::cerr << "[ERROR] nomt_counts missing for " << C << " " << kBin[ib] << " y" << iy << " (qcd_abcd.C before the per-bin level?)\n"; missing = true; continue; }
        for (int p = 0; p < 2; ++p)
        {
          CR cr;
          cr.dir = Form("%s_%s_y%d_CR%s", C.c_str(), kBin[ib], iy, p ? "H" : "D");
          cr.data = countBin(hc, p ? "dataHi" : "dataLo");
          cr.qcd = countBin(hc, p ? "qcdHi" : "qcdLo");
          cr.ewk = countBin(hc, p ? "ewkHi" : "ewkLo");
          crs.push_back(cr);
          // the SR template of the bin must add up to what qcd_abcd_nomt.h wrote
          double sr = 0;
          for (const Region &r : regs)
            if (r.dir == Form("%s_%s_y%d_%s", C.c_str(), kBin[ib], iy, kPart[p])) sr += r.qcd->Integral(1, 50);
          const double srRef = countBin(hc, p ? "srQcdHi" : "srQcdLo");
          if (std::fabs(sr - srRef) > 1e-6 * std::max(1.0, srRef))
            std::cerr << "[WARN] " << cr.dir << ": SR multijet template sums to " << sr << ", qcd_abcd_nomt.h says " << srRef << "\n";
        }
      }
  }
  if (missing) { std::cerr << "[ERROR] incomplete multijet inputs -> " << outFile << " NOT written\n"; return; }

  // ---- write ----
  TFile *fo = TFile::Open(outFile.c_str(), "RECREATE");
  auto writeRegion = [&](const Region &r, bool withSyst) {
    TDirectory *d = fo->mkdir(r.dir.c_str());
    d->cd();
    auto wn = [&](TH1D *h, const std::string &nm) {
      if (!h) return;
      TH1D *c = (TH1D *)h->Clone(nm.c_str());
      c->SetDirectory(d);
      if (nm != "data_obs" && c->Integral() <= 0.0)
      {
        std::cerr << "[WARN] " << outFile << " " << r.dir << "/" << nm << " is empty -> flooring central bin to 1e-6 for Combine\n";
        c->SetBinContent(c->GetNbinsX() / 2, 1e-6);
      }
      c->Write(nm.c_str(), TObject::kOverwrite);
    };
    wn(r.data, "data_obs"); wn(r.sig, "signal"); wn(r.z, "z"); wn(r.ztau, "ztau");
    wn(r.wtau, "wtau"); wn(r.ttbar, "ttbar"); wn(r.qcd, "qcd");
    if (withSyst) for (const auto &kv : r.syst) wn(kv.second, kv.first);
  };
  for (const Region &r : regs) writeRegion(r, true);
  for (const Region &r : sums) writeRegion(r, false);
  Region inc[3];
  inc[0].dir = "Wp_incl"; inc[1].dir = "Wm_incl"; inc[2].dir = "W_incl";
  for (const Region &s : sums)
  {
    if (s.dir.find("_lab_") == std::string::npos) continue;
    for (int k = 0; k < 3; ++k)
    {
      if (k == 0 && s.dir.compare(0, 2, "Wp") != 0) continue;
      if (k == 1 && s.dir.compare(0, 2, "Wm") != 0) continue;
      AddTo(inc[k].data, s.data, inc[k].dir + "_data"); AddTo(inc[k].sig, s.sig, inc[k].dir + "_signal");
      AddTo(inc[k].z, s.z, inc[k].dir + "_z"); AddTo(inc[k].ztau, s.ztau, inc[k].dir + "_ztau");
      AddTo(inc[k].wtau, s.wtau, inc[k].dir + "_wtau"); AddTo(inc[k].ttbar, s.ttbar, inc[k].dir + "_ttbar");
      AddTo(inc[k].qcd, s.qcd, inc[k].dir + "_qcd");
    }
  }
  for (int k = 0; k < 3; ++k) writeRegion(inc[k], false);
  for (const CR &cr : crs)
  {
    TDirectory *d = fo->mkdir(cr.dir.c_str());
    d->cd();
    auto one = [&](const char *nm, double v, double e) {
      TH1D *h = new TH1D(nm, nm, 1, 0, 1);
      h->SetDirectory(d);
      h->SetBinContent(1, (std::string(nm) != "data_obs" && v <= 0) ? 1e-6 : v);
      h->SetBinError(1, e);
      h->Write(nm, TObject::kOverwrite);
    };
    one("data_obs", cr.data, std::sqrt(std::max(cr.data, 0.0)));
    one("qcd", cr.qcd, std::sqrt(std::max(cr.data, 0.0)));
    one("ewk", cr.ewk, 0.0);
  }
  fo->Close();
  delete fo;
  std::cout << "[INFO] Saved structured Combine input: " << outFile << "  (" << regs.size() << " SR channels (lo/hi), "
            << sums.size() << " lo+hi sums, Wp_incl/Wm_incl/W_incl, " << crs.size() << " anchor CR dirs)\n";

  // ---- sidecars ----
  {
    const std::string side = outBase + "/" + inName + "_systs.txt";
    std::ofstream sf(side.c_str());
    sf << "# Shape systematics in " << outFile << " (plotting/combine_input_nomt.C)\n"
       << "# <systematic> <processes carrying <process>_<systematic>Up/Down in every SR channel>\n"
       << "# a process token <hist>@<charge> restricts the row to that charge's channels\n";
    int n = 0;
    for (const std::string &sy : systs)
    {
      if (systComplete[sy] == (int)regs.size()) { sf << sy << " signal z ztau wtau ttbar\n"; ++n; }
      else if (systComplete[sy] > 0)
        std::cerr << "[WARN] shape syst " << sy << " complete in only " << systComplete[sy] << "/" << regs.size()
                  << " SR channels -> NOT listed\n";
      else std::cerr << "[WARN] shape syst " << sy << " absent (run skim/run_lhe_updown.sh after the skim?) -> NOT listed\n";
    }
    if (isElec)
      for (int ic = 0; ic < 2; ++ic)
      {
        if (ffComplete[ic] == (int)regs.size() / 2) { sf << ffName[ic] << " qcd@" << (ic ? "Wm" : "Wp") << "\n"; ++n; }
        else std::cerr << "[WARN] " << ffName[ic] << " complete in only " << ffComplete[ic] << " channels -> NOT listed\n";
      }
    std::cout << "[INFO] shape systematics in " << outFile << ": " << n << " listed in " << side << "\n";
  }
  {
    const std::string gside = outBase + "/" + inName + "_mtsplit.txt";
    std::ofstream gs(gside.c_str());
    gs << "# The no-m_T lepton-pT fit, scheme C (QCD_MODE=mtsplit) in " << outFile << " (plotting/combine_input_nomt.C)\n"
       << "# SR channels <C>_<binning>_y<i>_lo (m_T < 30) and _hi (m_T >= 30); per rapidity bin the anchors\n"
       << "# <C>_<binning>_y<i>_CRD (anti-iso, m_T < 30) and _CRH (anti-iso, m_T >= 30). Lines: <binning> <group> <bins>\n"
       << "# = the groups of the multijet pT shape and of the transfer lnN. Residual lnN of the m_T >= 30 multijet\n"
       << "# (correction/qcd_abcd_nomt.h runNoMt): kappaWin <C> = anti-iso window, coherent over the bins of the charge;\n"
       << "# kappaXfer <binning> <group> <C> = the m_T < 30 -> >= 30 transfer, measured on all bins (same k per group),\n"
       << "# one uncorrelated nuisance per group\n";
    for (int ic = 0; ic < 2; ++ic) gs << "#! kappaWin " << (ic ? "Wm" : "Wp") << " " << Form("%.3f", kWin[ic]) << "\n";
    for (int ib = 0; ib < 2; ++ib)
      for (int g = 0; g < pOAbcdY::kNGroups; ++g)
        for (int ic = 0; ic < 2; ++ic)
          gs << "#! kappaXfer " << kBin[ib] << " " << pOAbcdY::kGroupName[g] << " " << (ic ? "Wm" : "Wp") << " "
             << Form("%.3f", kXfer[ic][ib][g]) << "\n";
    for (int ib = 0; ib < 2; ++ib)
      for (int g = 0; g < pOAbcdY::kNGroups; ++g)
      {
        gs << kBin[ib] << " " << pOAbcdY::kGroupName[g];
        for (int iy : pOAbcdY::BinsOf(g)) gs << " " << iy;
        gs << "\n";
      }
    std::cout << "[INFO] group sidecar " << gside << " (kappaWin W+ " << Form("%.3f", kWin[0]) << ", W- " << Form("%.3f", kWin[1])
              << "; kappaXfer lab fwd/bwd/cen W+ " << Form("%.3f/%.3f/%.3f", kXfer[0][0][0], kXfer[0][0][1], kXfer[0][0][2])
              << ", W- " << Form("%.3f/%.3f/%.3f", kXfer[1][0][0], kXfer[1][0][1], kXfer[1][0][2]) << ")\n";
  }

  // ---- composition summary + prefit stacks (display) ----
  for (int ic = 0; ic < 2; ++ic)
  {
    const std::string C = ic ? "Wm" : "Wp";
    for (int p = 0; p < 3; ++p) // lo, hi, all
    {
      TH1D *acc[7] = {};
      for (const Region &r : (p < 2 ? regs : sums))
      {
        if (r.dir.compare(0, 2, C) != 0 || r.dir.find("_lab_") == std::string::npos) continue;
        if (p < 2 && r.dir.substr(r.dir.size() - 2) != kPart[p]) continue;
        TH1D *in[7] = {r.data, r.sig, r.z, r.ztau, r.wtau, r.ttbar, r.qcd};
        for (int k = 0; k < 7; ++k) AddTo(acc[k], in[k], Form("nomt_acc_%s_%d_%d", C.c_str(), p, k));
      }
      const char *pn[3] = {"m_{T} < 30", "m_{T} #geq 30", "all m_{T}"};
      const char *ps[3] = {"lo", "hi", "all"};
      double tot = 0;
      for (int k = 1; k < 7; ++k) tot += acc[k]->Integral(1, 50);
      std::cout << Form("[NOMT-IN] %s %-3s: data %.0f | W %.1f, DY %.1f, DYtau %.1f, Wtau %.1f, ttbar %.1f, multijet %.1f (%.0f%%) | data/prefit %.3f\n",
                        C.c_str(), ps[p], acc[0]->Integral(1, 50), acc[1]->Integral(1, 50), acc[2]->Integral(1, 50),
                        acc[3]->Integral(1, 50), acc[4]->Integral(1, 50), acc[5]->Integral(1, 50), acc[6]->Integral(1, 50),
                        100 * acc[6]->Integral(1, 50) / acc[0]->Integral(1, 50), acc[0]->Integral(1, 50) / tot);
      PlotStyle st;
      st.drawOpt = "hist";
      st.logy = true;
      st.normBkgToData = false;
      st.pullPad = true;
      st.xRangeLo = 24.0; st.xRangeHi = 100.0;
      SaveNicePlot1D_WithBkg(acc[0], {acc[1], acc[2], acc[3], acc[4], acc[6], acc[5]},
                             {"W signal", "DY", "DY #tau", "W #tau", "QCD (no-m_{T} ABCD)", "t#bar{t}"},
                             outBase + "/leppt/h_leppt_" + ps[p] + "_" + C + "_inclusive",
                             isElec ? "p_{T}^{e} (GeV)" : "p_{T}^{#mu} (GeV)", "Events / 2 GeV",
                             (C == "Wp" ? "W^{+}" : "W^{-}") + std::string(isElec ? " #rightarrow e#nu" : " #rightarrow #mu#nu"),
                             std::string("no-m_{T} pT fit, ") + pn[p], "prefit, lab bins summed",
                             {Form("Data: %.0f", acc[0]->Integral(1, 50)), Form("Total prefit: %.0f", tot),
                              Form("multijet: %.0f%%", 100 * acc[6]->Integral(1, 50) / acc[0]->Integral(1, 50))},
                             st, [](TCanvas *, TH1 *h) { if (h) { h->SetMinimum(1.0); h->SetMaximum(10.25 * h->GetMaximum()); } });
    }
  }
  std::cout << "[INFO] Done. " << outFile << "\n";
}
