// skim_common.h
//
// Shared utilities for the four pO skim channels (Wmu, Wel, Zmm, Zee).
// Everything is inline / constexpr — no .cc file needed.
// Wrapped in namespace pOSkim to avoid polluting the ROOT CLING global scope.
//
// Physics conventions, cut values, histogram binning, etc. are kept identical
// to the original per-channel macros (see skim/legacy/).
//
// This file is included from skim/skim.C.

#ifndef PO_SKIM_COMMON_H
#define PO_SKIM_COMMON_H

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TString.h"
#include "TVector2.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// ============================================================
// Enums shared with run_all.sh
// ============================================================

// NOTE: kept at file scope (not inside the namespace) so that
// `root -l -q 'skim.C+(kZmm, "...", kData)'` from the shell driver
// resolves the bare names without needing a namespace prefix.
enum SampleType
{
  kData,
  kDY,
  kWp,
  kWm,
  kDYtau,
  kWptau,
  kWmtau,
  kTT     // inclusive ttbar background (2026-10-16); appended so the values above are unchanged
};

enum ChannelType
{
  kZmm,
  kZee,
  kWmu,
  kWel
};

namespace pOSkim
{

// ============================================================
// Constants
// ============================================================

inline constexpr double MU_MASS  = 0.1056583745;
inline constexpr double ELE_MASS = 0.000511;

// pO -> CM-frame rapidity boost (used downstream; defined here for reference)
inline constexpr double Y_SHIFT_pO_CM = 0.3466;

// Default data file. SINGLE SOURCE OF TRUTH for the data-file path:
// run_all.sh greps this value rather than hardcoding its own (override a
// single run by exporting DATA_FILE). Change the data file here, and the MC
// files in ResolveMCSample below -- both live in this header.
//
// Everything lives under ~/pO_2026_Aug_20/ since 2026-09-21: the MC is the
// Aug-20 EMBEDDED production (POWHEG + Angantyr underlying event), and the
// data file sits beside it as a HARD LINK to the July-29 copy (same inode,
// no second 45 GB on disk; `ls -li` shows link count 2). The Aug-20
// production ships MC only, so the DATA IS STILL THE JULY-29
// RECONSTRUCTION -- the filename keeps saying so on purpose.
// Absolute paths so ROOT's TFile::Open and run_all.sh resolve them without
// ~ expansion. To read from EOS/lxplus again, point this and `inputBase` in
// ResolveMCSample at
//   root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2026_Aug_20/
// (the data still lives under .../pO_2026_July_29/ on EOS -- there is no
// hard link there, so the data path needs the July-29 prefix remotely).
inline const char *kDefaultDataFile =
    "/Users/zhenghuang/pO_2026_Aug_20/July_29_DATA.root";

// ============================================================
// Branch / tree introspection helpers
// ============================================================

inline bool HasBranch(TTree *t, const char *name)
{
  return (t && t->GetListOfBranches() && t->GetListOfBranches()->FindObject(name));
}

inline std::string FindBranchContaining(TTree *t, const std::string &needle)
{
  if (!t || !t->GetListOfBranches())
    return "";
  TObjArray *arr = t->GetListOfBranches();
  const int n = arr->GetEntries();
  for (int i = 0; i < n; ++i)
  {
    TBranch *br = (TBranch *)arr->At(i);
    if (!br)
      continue;
    std::string bn = br->GetName();
    if (bn.find(needle) != std::string::npos)
      return bn;
  }
  return "";
}

// ============================================================
// Geometry / kinematics helpers
// ============================================================

inline double DeltaR(double eta1, double phi1, double eta2, double phi2)
{
  const double dEta = eta1 - eta2;
  const double dPhi = TVector2::Phi_mpi_pi(phi1 - phi2);
  return std::sqrt(dEta * dEta + dPhi * dPhi);
}

// PF MET from particleFlowAnalyser/pftree. pfId is taken but unused (kept
// for backwards-compatible signature with the legacy macros).
inline TVector2 ComputePFMET(std::vector<int> *pfId,
                             std::vector<float> *pfPt,
                             std::vector<float> *pfPhi)
{
  (void)pfId;
  double sumPx = 0.0, sumPy = 0.0;
  const int nPF = (pfPt ? (int)pfPt->size() : 0);

  for (int i = 0; i < nPF; ++i)
  {
    const double pt  = pfPt->at(i);
    const double phi = pfPhi->at(i);
    sumPx += pt * std::cos(phi);
    sumPy += pt * std::sin(phi);
  }
  return TVector2(-sumPx, -sumPy);
}

// The CHARGED-only missing transverse momentum (2026-10-21, user question:
// would the multijet ABCD behave better with a MET built from charged PF
// candidates only?): minus the vector sum of the packed PF candidates of type
// 1 (charged hadron), 2 (electron), 3 (muon) -- i.e. tracks, the lepton's own
// candidate included; no eta cut (charged candidates end at |eta| ~ 3), no
// vertex association (none stored; pO pileup is small). Diagnostic only: it
// fills the h_iso_pt_mttrk_* cubes of the W skims, nothing else.
inline TVector2 ComputeChargedPFMET(const std::vector<int> *pfId,
                                    const std::vector<float> *pfPt,
                                    const std::vector<float> *pfPhi)
{
  double sumPx = 0.0, sumPy = 0.0;
  const int nPF = (pfPt && pfId) ? (int)std::min(pfPt->size(), pfId->size()) : 0;
  for (int i = 0; i < nPF; ++i)
  {
    const int id = pfId->at(i);
    if (id < 1 || id > 3) continue;
    const double pt = pfPt->at(i), phi = pfPhi->at(i);
    sumPx += pt * std::cos(phi);
    sumPy += pt * std::sin(phi);
  }
  return TVector2(-sumPx, -sumPy);
}

// The e/gamma-CALIBRATED PF MET (2026-10-21, user decision) -- the MET every
// skim channel and the MET-using correction macros fill since that date.
// In DATA the forest corrects the electron and photon energies (elePt/eleRawPt
// = 1.021 barrel / 1.094 endcap, phoEt/phoRawEt = 1.025 / 1.101; in MC both
// ~1.00), but the PF candidates keep the RAW energies (the PF electron ==
// eleRawPt, median ratio 1.000 in data and MC), so the plain PF MET of data is
// biased along every electron and photon (for W -> e nu |MET| and m_T come out
// low, most in the endcap). Correction: the PF candidate of type 2 (electron)
// closest to a reco electron within DR < kMetCalibDR gets pT x elePt/eleRawPt,
// and the type-4 (photon) candidate closest to a reco photon gets
// pT x phoEt/phoRawEt; each candidate at most once (electrons first), ratios
// outside [kMetCalibRatioLo, kMetCalibRatioHi] are not applied and counted.
// Muons need nothing (the PF muon == muPt, correction/edge_study/pf_raw_check_mu).
// NB the matched PF photon carries only 45-80% of phoRawEt (PF splits the
// cluster), so the photon part corrects that candidate's share: the energy
// SCALE, not the whole cluster. Applied to data AND MC. PO_MET_EGCALIB=off ->
// the plain PF MET (regression checks; the job log says so).
inline constexpr double kMetCalibDR       = 0.05;
inline constexpr double kMetCalibRatioLo  = 0.5;
inline constexpr double kMetCalibRatioHi  = 2.0;
inline bool MetEGCalibEnabled()
{
  const char *e = std::getenv("PO_MET_EGCALIB");
  return !(e && std::string(e) == "off");
}

// Counters of the correction over the events a job passes to it (printed as
// the job's [MET-EG] block)
struct MetCalibStats
{
  unsigned long long nEvt = 0, nEle = 0, nEleMatched = 0, nPho = 0, nPhoMatched = 0, nRejected = 0;
  double sumRatioE = 0, sumRatioG = 0;  // over the applied corrections
  double sumDMet = 0, sumDMetVec = 0;   // |MET_calib| - |MET_PF| and |MET_calib - MET_PF|
  void Print(std::ostream &os, const char *who) const
  {
    os << "[MET-EG] " << who << ": e/gamma-calibrated PF MET "
       << (MetEGCalibEnabled() ? "ON" : "OFF (PO_MET_EGCALIB=off -> plain PF MET)") << " over " << nEvt << " events\n";
    if (!nEvt) return;
    os << Form("[MET-EG]   electrons: %llu, with a matched PF electron (DR < %.2f) %.1f%%, <elePt/eleRawPt> %.4f\n", nEle,
               kMetCalibDR, nEle ? 100.0 * nEleMatched / nEle : 0.0, nEleMatched ? sumRatioE / nEleMatched : 0.0);
    os << Form("[MET-EG]   photons:   %llu, with a matched PF photon   (DR < %.2f) %.1f%%, <phoEt/phoRawEt> %.4f\n", nPho,
               kMetCalibDR, nPho ? 100.0 * nPhoMatched / nPho : 0.0, nPhoMatched ? sumRatioG / nPhoMatched : 0.0);
    os << Form("[MET-EG]   ratios outside [%.1f, %.1f] not applied: %llu;  <|MET_calib| - |MET_PF|> = %+.3f GeV, <|dMET vector|> = %.3f GeV\n",
               kMetCalibRatioLo, kMetCalibRatioHi, nRejected, sumDMet / nEvt, sumDMetVec / nEvt);
  }
};

// The vectors are the caller's own bound branches (a branch can hold ONE
// address: never bind elePt etc. a second time for this -- see
// BindEGammaMetBranches). The PF vectors must be entry-aligned with the
// EventTree (same entry read).
inline TVector2 ComputePFMETEGCalib(const std::vector<int> *pfId, const std::vector<float> *pfPt,
                                    const std::vector<float> *pfEta, const std::vector<float> *pfPhi,
                                    const std::vector<float> *elePt, const std::vector<float> *eleRawPt,
                                    const std::vector<float> *eleEta, const std::vector<float> *elePhi,
                                    const std::vector<float> *phoEt, const std::vector<float> *phoRawEt,
                                    const std::vector<float> *phoEta, const std::vector<float> *phoPhi,
                                    MetCalibStats *st = nullptr)
{
  const int nPF = (pfPt && pfPhi) ? (int)std::min(pfPt->size(), pfPhi->size()) : 0;
  std::vector<double> scale(nPF, 1.0);
  std::vector<char> claimed(nPF, 0);
  const bool on = MetEGCalibEnabled() && pfId && pfEta && (int)pfId->size() >= nPF && (int)pfEta->size() >= nPF;
  // the closest unclaimed PF candidate of type `id` within kMetCalibDR of (eta, phi), or -1
  auto closest = [&](int id, double eta, double phi) {
    int best = -1;
    double bdr = kMetCalibDR;
    for (int j = 0; j < nPF; ++j)
    {
      if (pfId->at(j) != id || claimed[j]) continue;
      const double d = DeltaR(eta, phi, pfEta->at(j), pfPhi->at(j));
      if (d < bdr) { bdr = d; best = j; }
    }
    return best;
  };
  auto correct = [&](int id, const std::vector<float> *pt, const std::vector<float> *raw, const std::vector<float> *eta,
                     const std::vector<float> *phi, unsigned long long &nObj, unsigned long long &nMatched, double &sumRatio) {
    if (!pt || !raw || !eta || !phi) return;
    const size_t n = std::min(std::min(pt->size(), raw->size()), std::min(eta->size(), phi->size()));
    for (size_t k = 0; k < n; ++k)
    {
      ++nObj;
      const int j = closest(id, eta->at(k), phi->at(k));
      if (j < 0) continue;
      const double r = (raw->at(k) > 0) ? pt->at(k) / raw->at(k) : -1.0;
      if (!(r >= kMetCalibRatioLo && r <= kMetCalibRatioHi)) { if (st) ++st->nRejected; continue; }
      ++nMatched;
      sumRatio += r;
      scale[j] = r;
      claimed[j] = 1;
    }
  };
  MetCalibStats dummy;
  MetCalibStats &s = st ? *st : dummy;
  if (on)
  {
    correct(2, elePt, eleRawPt, eleEta, elePhi, s.nEle, s.nEleMatched, s.sumRatioE);
    correct(4, phoEt, phoRawEt, phoEta, phoPhi, s.nPho, s.nPhoMatched, s.sumRatioG);
  }
  double px = 0.0, py = 0.0, px0 = 0.0, py0 = 0.0;
  for (int i = 0; i < nPF; ++i)
  {
    const double c = std::cos(pfPhi->at(i)), sn = std::sin(pfPhi->at(i)), pt = pfPt->at(i);
    px0 -= pt * c; py0 -= pt * sn;
    px -= scale[i] * pt * c; py -= scale[i] * pt * sn;
  }
  ++s.nEvt;
  s.sumDMet += std::hypot(px, py) - std::hypot(px0, py0);
  s.sumDMetVec += std::hypot(px - px0, py - py0);
  return TVector2(px, py);
}

// Status AND address together (the fast-skim pairing rule); false (nothing
// bound) if the branch is absent.
template <class T>
inline bool BindBranch(TTree *t, const char *name, T *addr)
{
  if (!HasBranch(t, name)) return false;
  t->SetBranchStatus(name, 1);
  t->SetBranchAddress(name, addr);
  return true;
}

// Binds what ComputePFMETEGCalib needs and the caller has not bound already:
// pfId + pfEta on the PF tree, eleRawPt and the four photon branches on the
// EventTree, and elePt/eleEta/elePhi ONLY when the caller passes their
// pointers (the muon channels; the electron channels bind those themselves and
// pass nullptr -- binding them again would steal the caller's binding).
// Returns the name of the first missing branch ("" = all bound); the caller
// treats a missing one as FATAL (a silent plain-MET fallback would change the
// physics without a message).
inline std::string BindEGammaMetBranches(TTree *tEv, TTree *tPF, std::vector<int> **pfId, std::vector<float> **pfEta,
                                         std::vector<float> **eleRawPt, std::vector<float> **phoEt,
                                         std::vector<float> **phoRawEt, std::vector<float> **phoEta,
                                         std::vector<float> **phoPhi, std::vector<float> **elePt = nullptr,
                                         std::vector<float> **eleEta = nullptr, std::vector<float> **elePhi = nullptr)
{
  if (!BindBranch(tPF, "pfId", pfId)) return "pfId";
  if (!BindBranch(tPF, "pfEta", pfEta)) return "pfEta";
  if (!BindBranch(tEv, "eleRawPt", eleRawPt)) return "eleRawPt";
  if (!BindBranch(tEv, "phoEt", phoEt)) return "phoEt";
  if (!BindBranch(tEv, "phoRawEt", phoRawEt)) return "phoRawEt";
  if (!BindBranch(tEv, "phoEta", phoEta)) return "phoEta";
  if (!BindBranch(tEv, "phoPhi", phoPhi)) return "phoPhi";
  if (elePt && !BindBranch(tEv, "elePt", elePt)) return "elePt";
  if (eleEta && !BindBranch(tEv, "eleEta", eleEta)) return "eleEta";
  if (elePhi && !BindBranch(tEv, "elePhi", elePhi)) return "elePhi";
  return "";
}

// Relative PF isolation with the Delta-beta pileup correction (MuonPOG
// convention), (chHad + max(0, neuHad + pho - 0.5*PU)) / pT.
// Generic over the lepton flavour: pass the lepton pT vector + four iso vectors.
// Verified 2026-07-06: this formula applied to the ele component branches
// reproduces the ntuplizer's elePFRelIsoWithDBeta exactly.
// If the PU vector is null (branch absent), degrades to the uncorrected
// (ch + neu + pho)/pT sum.
inline double RelIsoPF(int i,
                       std::vector<float> *lepPt,
                       std::vector<float> *lepPFChIso,
                       std::vector<float> *lepPFNeuIso,
                       std::vector<float> *lepPFPhoIso,
                       std::vector<float> *lepPFPUIso)
{
  if (!lepPt || !lepPFChIso || !lepPFNeuIso || !lepPFPhoIso)
    return 999.0;
  const double pt = lepPt->at(i);
  if (pt <= 0)
    return 999.0;
  const double pu      = (lepPFPUIso ? (double)lepPFPUIso->at(i) : 0.0);
  const double neutral = std::max(0.0, (double)lepPFNeuIso->at(i)
                                     + (double)lepPFPhoIso->at(i) - 0.5 * pu);
  return (lepPFChIso->at(i) + neutral) / pt;
}

// MT for a (lepton, MET) pair, given lepton phi & MET TVector2.
inline double TransverseMass(double lepPt, double lepPhi, const TVector2 &met)
{
  const double dphi = TVector2::Phi_mpi_pi(lepPhi - met.Phi());
  return std::sqrt(2.0 * lepPt * met.Mod() * (1.0 - std::cos(dphi)));
}

// ============================================================
// Event-level selection (pO filters)
// ============================================================

// Step 2 in the W cutflow; reused by Z macros as a preselection.
// AND of the existing pO filter flags. Missing branches are treated as
// "do not apply" with a single warning.
inline bool PassEventSelection_pO(bool &warnedOnce,
                                  bool has_pprimaryVertexFilter,        int pprimaryVertexFilter,
                                  bool has_pclusterCompatibilityFilter, int pclusterCompatibilityFilter)
{
  auto requireIfExists = [&](bool has, int val, const char *name) -> bool
  {
    if (!has)
    {
      if (!warnedOnce)
        std::cout << "[WARN] Missing event filter branch: " << name << " (will not apply)\n";
      return true;
    }
    return (val != 0);
  };

  bool ok = true;
  ok = ok && requireIfExists(has_pprimaryVertexFilter,        pprimaryVertexFilter,        "pprimaryVertexFilter");
  // pclusterCompatibilityFilter DISABLED 2026-08-18: the filter is problematic
  // in pO, so it is not applied. Branch wiring in the callers is kept so this
  // is a one-line re-enable.
  // ok = ok && requireIfExists(has_pclusterCompatibilityFilter, pclusterCompatibilityFilter, "pclusterCompatibilityFilter");
  (void)has_pclusterCompatibilityFilter; (void)pclusterCompatibilityFilter;

  warnedOnce = true;
  return ok;
}

// ============================================================
// Trigger helpers
// ============================================================

inline bool TriggerFired(int hltBit)
{
  return (hltBit != 0);
}

// Step 8 in the W cutflows: leading-lepton trigger-object DR match.
// Identical structure in the muon and electron W macros (only arg names differed).
inline bool PassLeadingLeptonTrigMatch(double matchDR,
                                       int iLead,
                                       std::vector<float> *lepEta,
                                       std::vector<float> *lepPhi,
                                       bool has_trgObjPt,  std::vector<double> *trgObjPt,
                                       bool has_trgObjEta, std::vector<double> *trgObjEta,
                                       bool has_trgObjPhi, std::vector<double> *trgObjPhi,
                                       bool &warnedNoTrigMatchInfo,
                                       const char *flavourTag = "lepton")
{
  if (iLead < 0 || !lepEta || !lepPhi)
    return false;

  const double eta = lepEta->at(iLead);
  const double phi = lepPhi->at(iLead);

  if (has_trgObjPt && has_trgObjEta && has_trgObjPhi &&
      trgObjPt && trgObjEta && trgObjPhi)
  {
    const int nObj = (int)trgObjPt->size();
    for (int j = 0; j < nObj; ++j)
    {
      const double dr = DeltaR(eta, phi, trgObjEta->at(j), trgObjPhi->at(j));
      if (dr < matchDR)
        return true;
    }
    return false;
  }

  if (!warnedNoTrigMatchInfo)
  {
    std::cout << "[WARN] No trigger-object collection and no per-" << flavourTag
              << " trigger-match info found. Cannot apply step 8 "
              << "(leading " << flavourTag << " matched to trigger). Treating as PASS.\n";
    warnedNoTrigMatchInfo = true;
  }
  return true;
}

// ============================================================
// Gen-matching helpers (W macros only)
// ============================================================

// Recursively check whether the gen particle at `idx` has an ancestor with |pdg| == |targetPdg|.
inline bool HasAncestor(int idx,
                        int targetPdg,
                        std::vector<int> *genPdg,
                        std::vector<std::vector<int>> *genMotherIdx)
{
  if (!genPdg || !genMotherIdx)
    return false;
  if (idx < 0 || (size_t)idx >= genPdg->size() || (size_t)idx >= genMotherIdx->size())
    return false;

  const auto &moms = genMotherIdx->at(idx);
  for (int midx : moms)
  {
    if (midx < 0 || (size_t)midx >= genPdg->size() || (size_t)midx >= genMotherIdx->size())
      continue;

    if (std::abs(genPdg->at(midx)) == std::abs(targetPdg))
      return true;

    if (HasAncestor(midx, targetPdg, genPdg, genMotherIdx))
      return true;
  }

  return false;
}

// Reco lepton -> gen lepton matching with dR<0.5 and |dpT/genPt|<0.5.
// `requiredAncestorPdg` is accepted for legacy-signature compatibility; the
// original macros do not gate the return value on it (the matched-mother
// printout was commented out). Behavior preserved here: we look up the
// best gen match by dR among gen particles of the requested |pdg| and
// matching charge, and return true if any such match exists.
inline bool PassGenRecoMatchingWithAncestor(
    int iLead,
    std::vector<float> *lepPt, std::vector<float> *lepEta,
    std::vector<float> *lepPhi, std::vector<int>  *lepCharge,
    std::vector<float> *genPt, std::vector<float> *genEta,
    std::vector<float> *genPhi, std::vector<int> *genChg,
    std::vector<int> *genPdg,
    std::vector<std::vector<int>> *genMotherIdx,
    int requiredGenPdg,
    int requiredAncestorPdg)
{
  (void)requiredAncestorPdg;
  if (!lepPt || !lepEta || !lepPhi || !lepCharge)
    return false;
  if (!genPt || !genEta || !genPhi || !genChg || !genPdg || !genMotherIdx)
    return false;
  if (iLead < 0 || (size_t)iLead >= lepPt->size())
    return false;

  const double lepPtV  = lepPt->at(iLead);
  const double lepEtaV = lepEta->at(iLead);
  const double lepPhiV = lepPhi->at(iLead);
  const int    lepChg  = lepCharge->at(iLead);

  const size_t ngen = std::min({genPt->size(),
                                genEta->size(),
                                genPhi->size(),
                                genChg->size(),
                                genPdg->size(),
                                genMotherIdx->size()});

  int    bestIdx = -1;
  double bestDR  = 1e9;

  for (size_t i = 0; i < ngen; ++i)
  {
    if (genPt->at(i) <= 0)
      continue;
    if (std::abs(genPdg->at(i)) != std::abs(requiredGenPdg))
      continue;
    if (genChg->at(i) != lepChg)
      continue;

    const double dEta   = genEta->at(i) - lepEtaV;
    const double dPhi   = TVector2::Phi_mpi_pi(genPhi->at(i) - lepPhiV);
    const double dR     = std::sqrt(dEta * dEta + dPhi * dPhi);
    const double dPtRel = std::fabs(genPt->at(i) - lepPtV) / genPt->at(i);

    if (dR     >= 0.5) continue;
    if (dPtRel >= 0.5) continue;

    if (dR < bestDR)
    {
      bestDR  = dR;
      bestIdx = (int)i;
    }
  }

  if (bestIdx < 0)
    return false;
  return true;
}

// ============================================================
// Sample -> filename / outPrefix-suffix mapping
// ============================================================

// The W skims overwrite the input file path per MC sample. The DY MC file
// name is flavour-dependent: muon channels read Aug_20_MC_DY_mu_Z.root,
// electron channels read Aug_20_MC_DY_ele_Z.root.
//
// `flavour` should be "mu" or "ele".
//
// For non-MC (kData) this is a no-op (caller leaves `fname` untouched and
// `outSuffix` empty).
struct SampleFileInfo
{
  std::string fname;     // empty => leave caller's input file alone
  std::string outSuffix; // appended to outPrefix
};

inline SampleFileInfo ResolveMCSample(SampleType sample, const char *flavour)
{
  SampleFileInfo info;
  // Local copy of the Aug-20 EMBEDDED production (POWHEG + Angantyr) under
  // ~/pO_2026_Aug_20/ (absolute path; see the kDefaultDataFile note above for
  // how to switch back to EOS). MC ONLY -- the data stays on July-29.
  // Verified 2026-09-21 before the switch: all 7 trees present and
  // entry-aligned, hltobject carries the Oxy paths, ttbar_w is 217 long with
  // the July-29 layout, and <weight> = 6376.1 / 5463 / 1174.5 reproduces
  // sigma = 6.376 / 5.464 / 1.175 nb, so mc_norm.h is unchanged.
  const std::string inputBase =
      "/Users/zhenghuang/pO_2026_Aug_20/";
  // Production tag carried by every MC filename (mirrored by LabelFromFname
  // in count_ngen.C, which strips it to recover the canonical ngen labels).
  const std::string mcPrefix = "Aug_20_MC_";

  const std::string flav = flavour;

  // The filenames use the SAME flavour token for DY and W ("mu"/"ele"/"tau");
  // normalize the legacy "ee"/"e" spellings to "ele". (The May-26 production
  // used "ee" for DY -- that translation is gone with the new filenames.)
  const std::string tok = (flav == "ee" || flav == "e") ? "ele" : flav;

  switch (sample)
  {
    case kData:
      return info;
    case kDY:
      info.fname     = inputBase + mcPrefix + "DY_" + tok + "_Z.root";  // Aug_20_MC_DY_ele_Z.root
      info.outSuffix = "_DY";
      return info;
    case kWp:
      info.fname     = inputBase + mcPrefix + "Wp_" + tok + ".root";
      info.outSuffix = "_Wp";
      return info;
    case kWm:
      info.fname     = inputBase + mcPrefix + "Wm_" + tok + ".root";
      info.outSuffix = "_Wm";
      return info;
    case kDYtau:
      info.fname     = inputBase + mcPrefix + "DY_tau_Z.root";
      info.outSuffix = "_DYtau";
      return info;
    case kWptau:
      info.fname     = inputBase + mcPrefix + "Wp_tau.root";
      info.outSuffix = "_Wptau";
      return info;
    case kWmtau:
      info.fname     = inputBase + mcPrefix + "Wm_tau.root";
      info.outSuffix = "_Wmtau";
      return info;
    case kTT:
      // Inclusive ttbar (POWHEG + Angantyr, same Aug-20 production, added
      // 2026-10-16): ONE file for every lepton channel, like the tau samples.
      // Verified before wiring it in: the same 17 trees, entry-aligned
      // (1172000), ttbar_w 217 long with element 0 == weight and the W/DY
      // layout fingerprints exact (idx 0 == 44 == 110, 38 == 39, 9 == 107),
      // <weight> = 338.98 = sigma(ttbar) in pb (|w| = 340.72 constant, 0.26%
      // negative) -- see kSigma_TT in mc_norm.h.
      info.fname     = inputBase + mcPrefix + "TT.root";
      info.outSuffix = "_TT";
      return info;
  }
  // Every MC SampleType must have a case above. Falling through would return
  // an empty fname, the caller would keep the DATA input path with an empty
  // output suffix, and the W skims would overwrite the data skim file with
  // data filled as MC -- so stop instead (2026-10-16, when kTT was added).
  std::cerr << "[FATAL] ResolveMCSample: no input file for SampleType "
            << static_cast<int>(sample) << " -- add its case in skim_common.h\n";
  std::exit(2);
  return info;
}

inline bool IsMC(SampleType sample) { return sample != kData; }

// ============================================================
// (pT, mT) cut-scan floor (W macros, 2026-08-04)
// ============================================================
// The W skims fill the joint scan planes h_pt_mt_* / h_pt_mt_antiiso_* for
// leading-lepton pT > kPtScanFloor instead of the nominal 25 GeV, so
// correction/ptmt_scan.C can probe LOWERING the pT cut. Every other W
// histogram and the printed cutflow keep the nominal pT > 25 selection.
// Keep equal to ptmt_scan.C's kPtLo; must land on a h_pt_mt_* bin edge
// (1 GeV bins -- projection cuts must be bin edges).
inline constexpr double kPtScanFloor = 20.0;

// ============================================================
// m_T cut of the lepton-pT discriminant (2026-10-13)
// ============================================================
// m_T > kPtDiscMtCut defines the leppt_mt40 selection: at RECO in skim.C (the
// `_mt40` histograms, m_T from the PF MET) and at GEN in skim/gen_xsec.C (the
// `_mt40` fiducial set, m_T from the W's own lepton and neutrino). One number
// for both, so the fiducial volume is the gen twin of the reco cut. The `mt40`
// in the histogram names states this value -- rename them if it changes.
inline constexpr double kPtDiscMtCut = 40.0;

// The no-m_T lepton-pT fit (2026-10-21, user decision: the nominal; scheme C of
// correction/qcd_nomt_study.C): the signal region has NO m_T cut and is fitted
// as two pT channels split at m_T = kNoMtSplit (h_leppt_mtlo_* = m_T < 30,
// h_leppt_mthi_* = m_T >= 30, sharing r); the m_T < 30 channel is the ABCD
// region B. = qcd_abcd.C yCut (30), a bin edge of every m_T axis.
inline constexpr double kNoMtSplit = 30.0;

// ============================================================
// ECAL barrel-endcap transition ("crack") veto for electrons (2026-09-14)
// ============================================================
// 1.4442 < |eta_SC| < 1.566 (CMS EGamma convention), evaluated on the
// SUPERCLUSTER eta (eleSCEta = the calorimeter position seen from the
// origin). eleEta is the track direction at the vertex and differs by
// ~z_vtx/(R_ECAL cosh eta) ~ 0.02 at |eta| 1.5 -- comparable to the crack
// half-width, so the veto must use the SC position.
// RECO level only, data and MC alike, on EVERY electron that enters an
// ID'd-electron definition: the W leading electron, the DY-veto legs, both
// Z legs, and the macros replicating the selection (njet_WZ, charge_flip,
// trig_eff_mb, isolation_ele). The gen fiducial stays |eta| < 2.4, common
// with the muon channel: the crack is an acceptance hole inside the fiducial
// (0.122 of each 0.4-wide |eta| in [1.2,1.6] bin), filled by the MC's
// intra-bin rapidity shape. Why veto at all: the crack efficiency (A.eps
// ~0.62 there) cannot be measured with tag-and-probe, which excludes the
// crack as well -- removing the electrons turns an unmeasurable efficiency
// into a small, smooth acceptance extrapolation.
inline constexpr double kEcalGapLo = 1.4442;
inline constexpr double kEcalGapHi = 1.566;
inline bool InEcalGap(double scEta)
{
  const double a = scEta < 0 ? -scEta : scEta;
  return a > kEcalGapLo && a < kEcalGapHi;
}
// The effective electron ID used by the W skim helpers: the MVA WP
// (kEleIDBranch below) AND not in the crack. Built once per event; passed
// wherever the bare ID used to go (DY-veto legs, tight-ID gate,
// leading-electron pick) -- so the veto lives in one place and the helpers'
// signatures stay unchanged.
inline void BuildEleIDNoGap(int nEle,
                            const std::vector<int>   *eleID,
                            const std::vector<float> *eleSCEta,
                            std::vector<int>         &out)
{
  out.assign(nEle > 0 ? nEle : 0, 0);
  if (!eleID || !eleSCEta) return;
  for (int i = 0; i < nEle; ++i)
    out[i] = (eleID->at(i) != 0 && !InEcalGap(eleSCEta->at(i))) ? 1 : 0;
}

// ============================================================
// THE electron ID + isolation (2026-09-29): forest WP90 MVA ID + WP90 MVA iso
// ============================================================
// Switched from eleMVAIdWP95 + continuous relIso < 0.095 (user decision
// 2026-09-29) so that the EGM `wp90iso` scale factor -- ONE SF for ID and
// isolation together, relative to a reconstructed electron -- can be applied
// (skim/electron_sf.h). The forest's eleMVAIdWP90 / eleMVAIsoWP90 are the
// heavy-ion XGBoost classifiers of HIElectronInfoProducer, NOT EGM VID; the
// ID+iso SF cross-check (correction/idiso_sf_*, 2026-09-25/26) measured
// SF = 0.959 +- 0.025 for exactly this selection with an inclusive Z
// tag-and-probe in our pO data, vs the EGM prediction 0.956 +- 0.009 for the
// same probes -- the pp SF describes the forest selection at the 2.6% level.
//   ID      : eleMVAIdWP90, crack-vetoed (BuildEleIDNoGap)
//   iso-pass: eleMVAIsoWP90 AND Delta-beta relIso < kEleIsoRelMax (0.3)
// The relIso < 0.3 requirement is the cross-check's definition too: it keeps
// relIso 0.3-1.0 OUTSIDE the selection so that range can serve as the ABCD
// anti-iso sideband (kEleAntiIsoLo/Hi) with no event in both regions; it
// costs ~0.1% of passing prompt electrons (99.9% of WP90-iso-passing W
// electrons have relIso < 0.3), so the EGM SF applies up to SF(relIso<0.3) ~ 1.
// Used by EVERY electron definition: the W lepton, the DY-veto legs, both Z
// legs and the macros that replicate them (njet_WZ, charge_flip,
// trig_eff_mb). Muons are untouched (TightID + relIso < 0.15).
inline constexpr const char *kEleIDBranch  = "eleMVAIdWP90";
inline constexpr const char *kEleIsoBranch = "eleMVAIsoWP90";
inline constexpr double      kEleIsoRelMax = 0.3;
// ABCD anti-iso sideband of the electron channel: relIso in [0.3, 1.0) with
// the WP90 ID (= the cross-check's `sbid`). The "gap" between the regions is
// the set of MVA-iso-FAILING electrons with relIso < 0.3.
inline constexpr double      kEleAntiIsoLo = 0.3;
inline constexpr double      kEleAntiIsoHi = 1.0;

// The per-electron iso decision (vector aligned with the electron
// collection). relIso via RelIsoPF (Delta-beta), the SAME value the skim
// histograms carry.
inline void BuildEleIsoPass(int nEle,
                            std::vector<int>   *eleMVAIso,
                            std::vector<float> *elePt,
                            std::vector<float> *ch,
                            std::vector<float> *neu,
                            std::vector<float> *pho,
                            std::vector<float> *pu,
                            std::vector<int>   &out)
{
  out.assign(nEle > 0 ? nEle : 0, 0);
  if (!eleMVAIso || !elePt) return;
  for (int i = 0; i < nEle; ++i)
    out[i] = (eleMVAIso->at(i) != 0 && RelIsoPF(i, elePt, ch, neu, pho, pu) < kEleIsoRelMax) ? 1 : 0;
}

// The x value of an electron in the skim's ABCD planes h_iso_{met,mt,pt,
// pt_mt40}_ele* (axis [kEleIsoAxisLo, 1.0), 0.005 bins). The iso-pass is a
// boolean MVA decision, not a relIso window, so it cannot be recovered by
// projecting a relIso axis; the planes therefore store ISO-PASS electrons at
// relIso - 0.3 (i.e. in [-0.3, 0)) and every other electron at its relIso
// ([0, 1)). Downstream (correction/qcd_abcd.C) the regions are then plain
// projections again: iso-pass = [-0.3, 0), the gap = [0, 0.3), the anti-iso
// sideband = [0.3, 1.0). Every edge in play (-0.3, 0, 0.3, 0.6, 1.0) is a bin
// edge. An iso-passing electron with relIso >= 0.3 cannot exist (kEleIsoRelMax).
inline constexpr double kEleIsoAxisLo = -kEleIsoRelMax;
inline constexpr int    kEleIsoAxisNb = 260; // (1.0 - (-0.3)) / 0.005
inline double EleIsoAxis(bool isoPass, double relIso)
{
  return isoPass ? relIso - kEleIsoRelMax : relIso;
}

// ============================================================
// Rapidity binning (W macros)
// ============================================================

inline constexpr int kNY = 12;

// Lab-frame rapidity edges (12 bins, |y| < 2.4)
inline const double kYEdges[kNY + 1] = {
    -2.4, -2.0, -1.6, -1.2, -0.8, -0.4,
     0.0,  0.4,  0.8,  1.2,  1.6,  2.0,
     2.4};

// CM-frame ("FB") rapidity edges (shifted by +0.3466)
inline const double kYEdgesFB[kNY + 1] = {
    -1.7068,
    -1.3646,
    -1.0223,
    -0.6801,
    -0.3379,
     0.0044,
     0.3466,
     0.6888,
     1.0311,
     1.3733,
     1.7155,
     2.0578,
     2.4000};

// Inclusive-or-right-edge bin finder over an [N+1]-array of edges.
// Returns -1 if y is outside the full range.
inline int FindYBin(double y, const double *edges, int n)
{
  for (int b = 0; b < n; ++b)
    if (y >= edges[b] && y < edges[b + 1])
      return b;
  if (y == edges[n])
    return n - 1; // include right edge
  return -1;
}

// ============================================================
// Cutflow table (W macros)
// ============================================================

// Pretty-print the 0..8 cutflow exactly as the original macros did, to a
// stream (cout for stdout, ofstream for the txt file in ./output/).
template <typename Stream>
inline void PrintCutflow(Stream &os, const unsigned long long N[9])
{
  os << "\n========== Cutflow (AN table logic) ==========\n";
  os << " i   Ni              (Ni-Ni-1)/N1        Ni/Ni-1\n";
  os << "--------------------------------------------------\n";

  const double N1 = (N[1] > 0 ? (double)N[1] : 1.0);

  for (int i = 0; i <= 8; ++i)
  {
    double fracStep = 0.0;
    double ratio    = 0.0;

    if (i >= 1 && N[i - 1] > 0)
      ratio = (double)N[i] / (double)N[i - 1];

    if (i >= 2)
      fracStep = ((double)N[i] - (double)N[i - 1]) / N1;

    os << " " << i << "   " << N[i];

    if (i >= 2)
      os << "        " << Form("%+7.2f%%", 100.0 * fracStep);
    else
      os << "        " << "   (n/a) ";

    if (i >= 1)
      os << "        " << Form("%7.2f%%", 100.0 * ratio);
    else
      os << "        " << "  (n/a) ";

    os << "\n";
  }

  os << "==================================================\n\n";
}

// The W cutflow split by the CHARGE of the leading lepton (2026-10-20, user
// request, while chasing the W- -> e nu deficit at the lab edges). The charge
// exists only once the leading lepton is picked, so the split starts at step 5
// (attributed to the leading tight lepton) and follows N[6..8] exactly; two
// lines that are NOT skim cuts are added for information: step 8 + m_T > 40
// (the leppt_mt40 selection) and step 8 + PF MET >= 30. Raw counts like N[],
// for all |eta| and for the edges |eta| > 2.0 (lab bins y0 and y11).
struct ChargeCutflow
{
  static const int kNStep = 6; // steps 5, 6, 7, 8, 8+m_T>40, 8+MET>=30
  unsigned long long n[kNStep][2][2] = {}; // [step][q: 0 = -, 1 = +][region: 0 = all, 1 = |eta| > 2.0]
  void Fill(int step, int charge, double eta)
  {
    const int q = (charge > 0) ? 1 : 0;
    ++n[step][q][0];
    if (std::fabs(eta) > 2.0) ++n[step][q][1];
  }
};

template <typename Stream>
inline void PrintChargeCutflow(Stream &os, const ChargeCutflow &c, const char *lep)
{
  const char *name[ChargeCutflow::kNStep] = {"5 tight lepton (leading)", "6 lead pT>25, |eta|<2.4", "7 + isolation",
                                             "8 + trigger match",        "  8 + m_T > 40 [info]",   "  8 + MET >= 30 [info]"};
  os << "========== Cutflow by the charge of the leading " << lep << " (raw counts) ==========\n";
  os << "  step                      |          all |eta|                   |       edges |eta| > 2.0\n";
  os << "                            |    N-      N+    N-/N+  eff-   eff+ |   N-     N+   N-/N+  eff-   eff+\n";
  for (int s = 0; s < ChargeCutflow::kNStep; ++s)
  {
    os << Form("  %-25s |", name[s]);
    for (int r = 0; r < 2; ++r)
    {
      const double m = (double)c.n[s][0][r], p = (double)c.n[s][1][r];
      // step efficiency w.r.t. the previous line (the two [info] lines w.r.t. step 8)
      const int prev = (s == 0) ? -1 : (s >= 4 ? 3 : s - 1);
      const double em = (prev >= 0 && c.n[prev][0][r] > 0) ? m / c.n[prev][0][r] : 0.0;
      const double ep = (prev >= 0 && c.n[prev][1][r] > 0) ? p / c.n[prev][1][r] : 0.0;
      os << Form(r == 0 ? " %8.0f %8.0f %6.4f" : " %6.0f %6.0f %6.4f", m, p, p > 0 ? m / p : 0.0);
      if (prev >= 0) os << Form(" %6.4f %6.4f |", em, ep);
      else           os << "    --     -- |";
    }
    os << "\n";
  }
  os << "  (eff = this line / the line above; the [info] lines / step 8. Not skim cuts: m_T, MET.)\n\n";
}

// Write the cutflow to ./output/<outPrefix>_<Data|MC>.txt (overwrite), with the
// per-charge table appended when given.
inline void WriteCutflowTxt(const std::string &outPrefix, bool isMC, const unsigned long long N[9],
                            const ChargeCutflow *byCharge = nullptr, const char *lep = "lepton")
{
  const std::string outDir  = "./output/";
  const std::string mcTag   = isMC ? "MC" : "Data";
  const std::string outFile = outDir + outPrefix + "_" + mcTag + ".txt";

  std::ofstream ftxtout(outFile.c_str(), std::ios::out | std::ios::trunc);
  if (!ftxtout.is_open())
  {
    std::cerr << "ERROR: cannot open " << outFile << std::endl;
    return;
  }
  PrintCutflow(ftxtout, N);
  if (byCharge) PrintChargeCutflow(ftxtout, *byCharge, lep);
  ftxtout.close();
}

} // namespace pOSkim

#endif // PO_SKIM_COMMON_H
