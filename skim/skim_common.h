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
#include "TVector2.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <cmath>
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
  kWmtau
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
// Currently the LOCAL copy under ~/pO_2026_May_26/ (written as an absolute path
// so ROOT's TFile::Open and run_all.sh resolve it without needing ~ expansion).
// To read from EOS/lxplus again, point this and `inputBase` in ResolveMCSample
// back at the prefix
//   root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2026_May_26/
inline const char *kDefaultDataFile =
    "/Users/zhenghuang/pO_2026_May_26/Data_May_26.root";

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
  ok = ok && requireIfExists(has_pclusterCompatibilityFilter, pclusterCompatibilityFilter, "pclusterCompatibilityFilter");

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
// name is flavour-dependent: muon channels read pO_MC_DY_mu_Z.root,
// electron channels read pO_MC_DY_ele_Z.root.
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
  // Local copy under ~/pO_2026_May_26/ (absolute path; see the kDefaultDataFile
  // note above for how to switch back to EOS).
  const std::string inputBase =
      "/Users/zhenghuang/pO_2026_May_26/";

  const std::string flav = flavour;

  // DY uses ee/mu/tau; W uses ele/mu/tau. Translate the electron token.
  const std::string dyFlav = (flav == "ele" || flav == "e") ? "ee"  : flav;
  const std::string wFlav  = (flav == "ee"  || flav == "e") ? "ele" : flav;

  switch (sample)
  {
    case kData:
      return info;
    case kDY:
      info.fname     = inputBase + "MC_DY" + dyFlav + "_May26.root";  // MC_DYee_May26.root
      info.outSuffix = "_DY";
      return info;
    case kWp:
      info.fname     = inputBase + "MC_Wp_" + wFlav + "_May26.root";
      info.outSuffix = "_Wp";
      return info;
    case kWm:
      info.fname     = inputBase + "MC_Wm_" + wFlav + "_May26.root";
      info.outSuffix = "_Wm";
      return info;
    case kDYtau:
      info.fname     = inputBase + "MC_DYtau_May26.root";
      info.outSuffix = "_DYtau";
      return info;
    case kWptau:
      info.fname     = inputBase + "MC_Wp_tau_May26.root";
      info.outSuffix = "_Wptau";
      return info;
    case kWmtau:
      info.fname     = inputBase + "MC_Wm_tau_May26.root";
      info.outSuffix = "_Wmtau";
      return info;
  }
  return info;
}

inline bool IsMC(SampleType sample) { return sample != kData; }

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

// Write the cutflow to ./output/<outPrefix>_<Data|MC>.txt (overwrite).
inline void WriteCutflowTxt(const std::string &outPrefix, bool isMC, const unsigned long long N[9])
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
  ftxtout.close();
}

} // namespace pOSkim

#endif // PO_SKIM_COMMON_H
