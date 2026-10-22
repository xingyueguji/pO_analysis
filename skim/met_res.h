// skim/met_res.h -- the MET-RESOLUTION variation of the MC (2026-10-20).
//
// The Z recoil says the MC MET is too narrow: the u_perp width (the recoil
// resolution, no q_T response in it) is k = sigma_data/sigma_MC = 1.10 +- 0.04
// in Z->mumu and 1.26 +- 0.06 in Z->ee (correction/recoil_res.C, log
// correction/logs/recoil_res.log). This header turns that into a shape
// nuisance for the MET-discriminant fit (user decision 2026-10-20: one
// nuisance per flavour, the two k's are independent measurements):
//   +1 sigma (Up)  = every MC event's MET vector smeared by an independent
//                    Gaussian of width sigma_add = sigma_MC sqrt(k^2 - 1) in
//                    EACH transverse component, i.e. the MC MET resolution
//                    brought to the data one; m_T is recomputed with the
//                    smeared MET, so events migrate across the m_T > 40 cut
//                    of the _mt40 selections.
//   -1 sigma (Down) = the mirror, 2 x nominal - Up per bin (floored at 0) --
//                    a resolution cannot be un-smeared.
// The smearing is deterministic: the two Gaussians of an event are drawn from
// a TRandom3 seeded with the event's entry number, so a re-skim reproduces the
// templates bit for bit and an event keeps its draw whatever else changes.
// Twins <h>_metResMuUp/Down (muon files) / <h>_metResEleUp/Down (electron
// files) are booked for h_met_*, h_met_mt40_* and h_leppt_mt40_* (lab + FB),
// MC only, filled with the full event weight (gen x SF). NOT area-normalized:
// a migration across the m_T cut is part of the effect. plotting/mtandmet.C
// attaches them to the met_mt40 Combine input only (user decision 2026-10-20:
// the MET fit only; the leppt_mt40 twins are stored so the lepton-pT fit can
// switch them on later).
//
// Input  : correction/rootfile/recoil_res.root (TParameter sigmaAdd_<flav>,
//          k_<flav>, kErr_<flav>, sigmaMC_<flav>). Missing => FATAL in the
//          W MC skim (never a silent zero smearing). Run
//          correction/run_recoil_res.sh first.
// Env    : PO_MET_RES=off skips the twins (regression checks); the job says so.

#ifndef PO_MET_RES_H
#define PO_MET_RES_H

#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TParameter.h"
#include "TRandom3.h"
#include "TString.h"
#include "TSystem.h"
#include "TVector2.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace pOMetRes
{
inline const char *kInputCandidates[] = {
    "../correction/rootfile/recoil_res.root",
    "correction/rootfile/recoil_res.root",
    nullptr};
// the per-event seed = kSeedBase + entry (TRandom3 treats seed 0 as "use the clock")
inline const UInt_t kSeedBase = 20261020u;

inline const char *SystName(bool isElec) { return isElec ? "metResEle" : "metResMu"; }

inline bool Enabled()
{
  const char *e = gSystem->Getenv("PO_MET_RES");
  return !(e && (TString(e) == "off" || TString(e) == "0" || TString(e) == "OFF"));
}

struct MetRes
{
  std::string flav, name;
  double k = 0.0, kErr = 0.0, sigmaMC = 0.0, sigmaAdd = 0.0;
  std::string file;
  mutable TRandom3 rng;

  bool Load(bool isElec, std::ostream &log)
  {
    flav = isElec ? "ele" : "mu";
    name = SystName(isElec);
    for (int i = 0; kInputCandidates[i]; ++i)
      if (!gSystem->AccessPathName(kInputCandidates[i])) { file = kInputCandidates[i]; break; }
    if (file.empty())
    {
      log << "[MetRes] ERR recoil_res.root not found (run correction/run_recoil_res.sh)\n";
      return false;
    }
    TFile *f = TFile::Open(file.c_str(), "READ");
    if (!f || f->IsZombie()) { log << "[MetRes] ERR cannot open " << file << "\n"; return false; }
    auto get = [&](const std::string &nm, double &v) -> bool {
      auto *p = dynamic_cast<TParameter<double> *>(f->Get(nm.c_str()));
      if (!p) return false;
      v = p->GetVal();
      return true;
    };
    const bool ok = get("sigmaAdd_" + flav, sigmaAdd) && get("k_" + flav, k) && get("kErr_" + flav, kErr) &&
                    get("sigmaMC_" + flav, sigmaMC);
    f->Close();
    delete f;
    if (!ok) { log << "[MetRes] ERR " << file << " lacks the " << flav << " parameters\n"; return false; }
    log << Form("[MetRes] %s: k = %.4f +- %.4f (Z u_perp width data/MC), sigma_MC = %.3f GeV -> +1 sigma = MET"
                " smeared by sigma_add = %.3f GeV per component (seed = %u + entry); twins <h>_%sUp/Down; input %s\n",
                flav.c_str(), k, kErr, sigmaMC, sigmaAdd, kSeedBase, name.c_str(), file.c_str());
    return true;
  }

  // the smeared MET vector of the event at `entry` (deterministic)
  TVector2 Smear(const TVector2 &met, Long64_t entry) const
  {
    if (sigmaAdd <= 0.0) return met;
    rng.SetSeed(kSeedBase + (UInt_t)entry);
    const double dx = rng.Gaus(0.0, sigmaAdd), dy = rng.Gaus(0.0, sigmaAdd);
    return TVector2(met.X() + dx, met.Y() + dy);
  }
};

struct Twin
{
  const TH1 *nom = nullptr;
  TH1D *up = nullptr;
};

// Book the Up twin of `nom` (same axes, empty); the Down is built at the end.
inline Twin Book(const TH1 *nom, const std::string &systName, std::vector<Twin> &registry)
{
  Twin t;
  t.nom = nom;
  const std::string title = nom->GetTitle();
  const size_t semi = title.find(';');
  const std::string note = " [" + systName + " +1 sigma: MC MET smeared to the data resolution]";
  t.up = (TH1D *)nom->Clone((std::string(nom->GetName()) + "_" + systName + "Up").c_str());
  t.up->Reset("ICES");
  t.up->SetTitle(semi == std::string::npos ? (title + note).c_str()
                                           : (title.substr(0, semi) + note + title.substr(semi)).c_str());
  if (t.up->GetSumw2N() == 0) t.up->Sumw2();
  registry.push_back(t);
  return t;
}

inline void Fill(const Twin &t, double x, double w)
{
  if (t.up) t.up->Fill(x, w);
}

// Down = 2 x nominal - Up per bin (under/overflow included), floored at 0;
// bin errors of both = the nominal's (Combine reads the contents). Writes
// Up and Down into the current directory; returns the number written.
inline int FinalizeAndWrite(std::vector<Twin> &registry, const std::string &systName)
{
  int n = 0;
  for (Twin &t : registry)
  {
    if (!t.nom || !t.up) continue;
    TH1D *dn = (TH1D *)t.up->Clone((std::string(t.nom->GetName()) + "_" + systName + "Down").c_str());
    dn->SetTitle(TString(t.up->GetTitle()).ReplaceAll("+1 sigma: MC MET smeared to the data resolution",
                                                     "-1 sigma: the mirror 2 x nominal - Up").Data());
    const int nb = t.nom->GetNbinsX();
    for (int b = 0; b <= nb + 1; ++b)
    {
      const double v0 = t.nom->GetBinContent(b), e0 = t.nom->GetBinError(b);
      dn->SetBinContent(b, std::max(0.0, 2.0 * v0 - t.up->GetBinContent(b)));
      dn->SetBinError(b, e0);
      t.up->SetBinError(b, e0);
    }
    t.up->Write("", TObject::kOverwrite);
    dn->Write("", TObject::kOverwrite);
    delete dn;
    n += 2;
  }
  return n;
}

// Weighted bookkeeping of the m_T > 40 migration (the record of each job)
struct Stats
{
  double wNom40 = 0.0, wIn = 0.0, wOut = 0.0, wAll = 0.0, dMet = 0.0;
  void Add(double w, bool pass40, bool pass40s, double met, double metS)
  {
    wAll += w;
    dMet += w * (metS - met);
    if (pass40) wNom40 += w;
    if (!pass40 && pass40s) wIn += w;
    if (pass40 && !pass40s) wOut += w;
  }
  void Print(std::ostream &log, const std::string &who) const
  {
    log << Form("[MetRes] %s: smeared %.6g (weighted) iso-pass events; <MET_smeared - MET> = %+.3f GeV;"
                " m_T>40: nominal %.6g, in %.6g (%.2f%%), out %.6g (%.2f%%) -> net %+.2f%%\n",
                who.c_str(), wAll, wAll > 0 ? dMet / wAll : 0.0, wNom40, wIn, wNom40 > 0 ? 100.0 * wIn / wNom40 : 0.0,
                wOut, wNom40 > 0 ? 100.0 * wOut / wNom40 : 0.0, wNom40 > 0 ? 100.0 * (wIn - wOut) / wNom40 : 0.0);
  }
};
} // namespace pOMetRes

#endif // PO_MET_RES_H
