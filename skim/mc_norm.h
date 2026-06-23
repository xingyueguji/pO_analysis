// skim/mc_norm.h
//
// Single source of truth for the MC -> data ABSOLUTE normalization.
//
//     k_s = sigma_s * L_int / N_gen,s        (per physical MC sample s)
//
// Applied DOWNSTREAM as a TH1::Scale(k_s) on each per-sample MC template, so the
// template carries the absolute predicted yield at the data luminosity:
//
//     N_pred,s(bin) = k_s * sum_{selected} w_gen  =  sigma_s * L_int * (A x eps)
//
// This produces FIXED, absolutely-normalized templates. In the downstream fit
// (the separate Combine fork, intentionally NOT touched here) the backgrounds
// stay frozen at this prediction and only the signal RATE (signal strength)
// floats. Because the backgrounds are frozen, getting this absolute scale RIGHT
// is REQUIRED for the fit to be unbiased -- it is not merely cosmetic.
//
// Ingredients:
//   - N_gen,s : skim/count_ngen.C -> skim/rootfile/ngen.root, per-sample scalar
//               TParameter<double> named "ngen_<label>" (label = file-derived,
//               e.g. "Wp_mu", "DYee", "DYtau", "DYmu_low").
//   - L_int   : data luminosity (46.5 nb^-1).   [the only data-derived number]
//   - sigma_s : per-process effective cross section. FILLED (2026-06-23) from the
//               MC weights — this production's per-event weight IS a cross-section
//               weight with <w> = sigma; see kSigma_* below.
//
// CONVENTION (must be self-consistent or the normalization is silently wrong):
//   * L_int and every sigma_s MUST be in the SAME frame -- per-pO OR
//     per-nucleon-nucleon (the A=16 factor lives in whichever you pick) --
//     matching however the 46.5 nb^-1 was measured.   TODO: confirm.
//   * Units below: L in nb^-1, sigma in nb  =>  k_s dimensionless.
//   * sigma_s is the EFFECTIVE cross section for the sample's final state: it
//     already includes the branching ratio into that channel AND any
//     generator-level filter efficiency.
//   * Lepton universality: sigma(W+ -> e nu) = sigma(W+ -> mu nu) =
//     sigma(W+ -> tau nu), and likewise DY -> ee/mumu/tautau -- so ONE sigma per
//     process serves all three lepton-flavour FILES (their N_gen differ and are
//     handled per file).   TODO: confirm the samples are exclusive-per-lepton.
//
// The kSigma_* are now filled from the MC weights (below). MCScale() only falls
// back to 1.0 + warns for a label whose sigma is left kUnset (e.g. low-mass DY,
// which is skipped). Re-tuning = editing the numbers here and re-running the
// template macros (NO re-skim needed).
//
// KEY FACT for this production: because <w> = sigma, the formula k_s = sigma*L/N_gen
// reduces algebraically to k_s = L / N_raw,s  (sigma cancels: sigma*L/Sum(w) =
// sigma*L/(N_raw*sigma) = L/N_raw). Both give the same answer; the sigma form is
// kept because it is robust to the weight unit (pb vs nb cancels in the
// efficiency ratio Sum_sel(w)/Sum_all(w)), so only sigma[nb] + L[nb^-1] matter.

#ifndef PO_MC_NORM_H
#define PO_MC_NORM_H

#include "TFile.h"
#include "TParameter.h"

#include <iostream>
#include <string>

namespace pONorm
{

// -------- data luminosity (the only data-derived number) --------
inline constexpr double kLumi_invnb = 46.5; // [nb^-1] -- the pO (hadronic) luminosity

// -------- proton-Oxygen A-scaling -------------------------------
// The kSigma_* below are per-NUCLEON-NUCLEON (<w> ~ 6 nb, pp-scale; the nuclear
// PDF modification EPPS21 is already folded into the POWHEG weights). For a hard
// process the proton can scatter off any of Oxygen's A nucleons, so the cross
// section A-scales: sigma_pO = A * sigma_NN, with A = 16 for O-16. Multiplying
// by A (with the pO luminosity above) gives the absolute predicted pO yield.
// If your luminosity is instead the per-nucleon-nucleon L_NN, set kA_O = 1.0.
inline constexpr double kA_O = 16.0; // Oxygen mass number (per-NN sigma -> per-pO)

// -------- per-process effective cross sections [nb] --------
// POWHEG generator cross sections, read straight from the MC weights: this
// production's per-event weight is a cross-section weight with <w> = sigma (pb),
// so sigma_s = <w>_s (skim/output/ngen.txt, averaged over the e/mu/tau files):
//   Wp: <w> ~ 6375.6 pb,  Wm: ~ 5463.6 pb,  DY(m>50): ~ 1174.7 pb.
// Consistent across the three lepton files (universality); the ~0.9% NLO
// negative-weight fraction is already folded into <w> and into N_gen = Sum w.
//
// FRAME: these are ~per-nucleon-nucleon magnitudes (W+ -> l nu ~ 6 nb, pp-scale).
// If the 46.5 nb^-1 is the per-pO luminosity, the absolute MC sits ~A=16 (or
// ~N_coll) below data -- that factor IS the frame, and for a stack normalized to
// data it cancels. Confirm against data before using these absolutely in the fit.
inline constexpr double kUnset       = -1.0;
inline constexpr double kSigma_Wp    = 6.376; // sigma(W+ -> l nu)      = <w>_Wp  (POWHEG, nb)
inline constexpr double kSigma_Wm    = 5.464; // sigma(W- -> l nu)      = <w>_Wm
inline constexpr double kSigma_DY    = 1.175; // sigma(DY -> l l, m>50) = <w>_DY
inline constexpr double kSigma_DYlow = kUnset; // low-mass DY -- SKIPPED (the _low files are left unwired)

// Map a count_ngen file-label ("Wp_mu","DYee","DYtau","DYmu_low",...) to its
// effective cross section. One sigma per process (universality); _low overrides.
inline double SigmaForLabel(const std::string &label)
{
  if (label.find("_low") != std::string::npos) return kSigma_DYlow; // any DY low-mass file
  if (label.rfind("Wp", 0) == 0)               return kSigma_Wp;    // Wp_mu / Wp_ele / Wp_tau
  if (label.rfind("Wm", 0) == 0)               return kSigma_Wm;    // Wm_mu / Wm_ele / Wm_tau
  if (label.rfind("DY", 0) == 0)               return kSigma_DY;    // DYmu / DYee / DYtau
  std::cout << "[WARN] pONorm::SigmaForLabel: unknown label '" << label
            << "' -> sigma unset\n";
  return kUnset;
}

// Read N_gen,s = ngen_<label> (TParameter<double>) from ngen.root.
// ngenFile default assumes a caller one directory below skim/ (plotting/,
// correction/, analysis/). Returns <=0 if the file/param is missing.
inline double NgenForLabel(const std::string &label,
                           const char *ngenFile = "../skim/rootfile/ngen.root")
{
  TFile *f = TFile::Open(ngenFile);
  if (!f || f->IsZombie())
  {
    std::cout << "[WARN] pONorm::NgenForLabel: cannot open " << ngenFile
              << " (run skim/run_ngen.sh first)\n";
    if (f) { f->Close(); delete f; }
    return -1.0;
  }
  TParameter<double> *p =
      (TParameter<double> *)f->Get((std::string("ngen_") + label).c_str());
  const double ngen = p ? p->GetVal() : -1.0;
  if (!p)
    std::cout << "[WARN] pONorm::NgenForLabel: no ngen_" << label
              << " in " << ngenFile << "\n";
  f->Close();
  delete f;
  return ngen;
}

// k_s = A_O * sigma_s * L_int / N_gen,s  for the given file-label (absolute pO
// event yield). Returns 1.0 (no-op) + warns if sigma is unset or N_gen
// unavailable, so it is safe to wire in before the cross sections are provided.
inline double MCScale(const std::string &label,
                      const char *ngenFile = "../skim/rootfile/ngen.root")
{
  const double sigma = SigmaForLabel(label);
  if (sigma <= 0.0)
  {
    std::cout << "[WARN] pONorm::MCScale: sigma unset for '" << label
              << "' -> k=1.0 (no-op; fill kSigma_* in skim/mc_norm.h)\n";
    return 1.0;
  }
  const double ngen = NgenForLabel(label, ngenFile);
  if (ngen <= 0.0)
  {
    std::cout << "[WARN] pONorm::MCScale: N_gen unavailable for '" << label
              << "' -> k=1.0 (no-op)\n";
    return 1.0;
  }
  return kA_O * sigma * kLumi_invnb / ngen; // absolute pO yield scale
}

} // namespace pONorm

#endif // PO_MC_NORM_H
