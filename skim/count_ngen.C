// skim/count_ngen.C
//
// Compute N_gen for every MC sample = the sum of the per-event generator
// weight (hiEvtAnalyzer/HiTree::weight) over ALL events in the file, with NO
// analysis selection applied. This is the generator-level denominator used for
// absolute cross-section normalization:
//
//     sigma = N_selected / (N_gen * acceptance * ...);   N_gen = Sum_i w_i
//
// where the sum runs over the entire generated sample. Per the project's
// gen-weight convention (see skim/skim.C), data carries no generator weight,
// so only MC files are counted here.
//
// N_gen is a property of a *physical MC file*, not of a channel/selection:
//   - the W->mu / Z->mu skims both read MC_Wp_mu / MC_Wm_mu / MC_DYmu,
//   - the W->e  / Z->e  skims both read MC_Wp_ele / MC_Wm_ele / MC_DYee,
//   - the tau backgrounds (MC_*tau) are shared by every channel.
// So we enumerate the samples through ResolveMCSample(sample, flavour) and
// dedupe by the resolved filename -> 9 distinct MC files, each counted once.
//
// Run from skim/:
//     root -l -q -b 'count_ngen.C+'
//   (or ./run_ngen.sh)
//
// Outputs (both written, dirs created if absent):
//   output/ngen.txt   -- human-readable table (label, file, N_raw, N_gen, ...)
//   rootfile/ngen.root -- h_ngen / h_nraw (one labeled bin per file) plus a
//                         TParameter<double> ngen_<label> / nraw_<label> each,
//                         so downstream code can Get() a value by name.
//
// The per-file mean/min/max weight and the negative-weight count are also
// reported -- this is the diagnostic the CLAUDE.md TODO asks for ("verify the
// weight distribution: all == 1 -> weighting is a no-op; spread or negative
// weights -> it matters").

//k_s = σ · L / N_gen = (N_gen / N_raw) · L / N_gen = L / N_raw

#include "skim_common.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TParameter.h"
#include "TSystem.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

using namespace pOSkim;

namespace
{

// One physical MC file to count. flavour ("mu"/"ele") only matters for the
// flavour-split DY/W samples; it is ignored by ResolveMCSample for tau.
struct McJob
{
  SampleType  sample;
  const char *flavour;
};

// Result of counting one file.
struct NgenResult
{
  std::string label;   // short, file-derived (e.g. "Wp_mu", "DYtau")
  std::string fname;   // full resolved path that was opened
  bool        ok;      // file opened and HiTree found
  bool        hasWeight;
  Long64_t    nraw;    // raw event count (HiTree entries)
  double      sumw;    // N_gen = Sum w_i   (= nraw when no weight branch)
  double      sumw2;   // Sum w_i^2  (-> stat error sqrt(sumw2) on N_gen)
  double      wmin;
  double      wmax;
  Long64_t    nNeg;    // # events with w < 0
};

// Derive the CANONICAL sample label from the resolved filename basename. The
// labels ("Wp_mu", "Wm_ele", "DYmu", "DYee", "DYtau", "Wp_tau", "Wm_tau") are
// load-bearing: mc_norm.h consumers request ngen_<label> by these exact
// strings (MCScale("DYee") etc.), so they must stay STABLE across productions
// even though the filenames change.
//   May-26:  MC_Wp_mu_May26.root      -> "Wp_mu",  MC_DYee_May26.root  -> "DYee"
//   July-29: July_29_MC_Wp_mu.root    -> "Wp_mu",  July_29_MC_DY_ele_Z.root -> "DYee"
std::string LabelFromFname(const std::string &fname)
{
  std::string base = fname;
  const size_t slash = base.find_last_of('/');
  if (slash != std::string::npos)
    base = base.substr(slash + 1);

  // strip known prefixes (July-29 first: it contains the May-era "MC_" too)
  for (const std::string &pre : {std::string("July_29_MC_"), std::string("MC_")})
    if (base.rfind(pre, 0) == 0)
    {
      base = base.substr(pre.size());
      break;
    }

  // strip known suffixes, longest first
  for (const std::string &suf : {std::string("_May26.root"), std::string("_Z.root"),
                                 std::string(".root")})
    if (base.size() >= suf.size() &&
        base.compare(base.size() - suf.size(), suf.size(), suf) == 0)
    {
      base = base.substr(0, base.size() - suf.size());
      break;
    }

  // canonicalize the July-29 DY spellings to the legacy labels
  if (base == "DY_mu")  return "DYmu";
  if (base == "DY_ele") return "DYee";
  if (base == "DY_tau") return "DYtau";
  return base;
}

// Count Sum(weight) over all HiTree entries in one already-resolved MC file.
NgenResult CountOne(const std::string &label, const std::string &fname)
{
  NgenResult r;
  r.label     = label;
  r.fname     = fname;
  r.ok        = false;
  r.hasWeight = false;
  r.nraw      = 0;
  r.sumw      = 0.0;
  r.sumw2     = 0.0;
  r.wmin      = 0.0;
  r.wmax      = 0.0;
  r.nNeg      = 0;

  std::cout << "[INPUT] " << fname << std::endl;
  TFile *f = TFile::Open(fname.c_str());
  if (!f || f->IsZombie())
  {
    std::cerr << "[ERR] count_ngen: cannot open " << fname << " (skipping)\n";
    if (f) { f->Close(); delete f; }
    return r;
  }

  TTree *tHi = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  if (!tHi)
  {
    std::cerr << "[ERR] count_ngen: no hiEvtAnalyzer/HiTree in " << fname << " (skipping)\n";
    f->Close();
    delete f;
    return r;
  }

  // Same wiring as skim.C: Float_t weight, soft-gated on the branch existing.
  // Fast read: disable everything, enable only 'weight'.
  Float_t    genWeight = 1.f;
  const bool hasWeight = HasBranch(tHi, "weight");
  tHi->SetBranchStatus("*", 0);
  if (hasWeight)
  {
    tHi->SetBranchStatus("weight", 1);
    tHi->SetBranchAddress("weight", &genWeight);
  }
  else
  {
    std::cout << "[WARN] count_ngen: '" << label
              << "' has no 'weight' branch on HiTree; counting unweighted (w=1).\n";
  }

  const Long64_t n = tHi->GetEntries();
  double sumw = 0.0, sumw2 = 0.0;
  double wmin = 1e300, wmax = -1e300;
  Long64_t nNeg = 0;

  for (Long64_t i = 0; i < n; ++i)
  {
    tHi->GetEntry(i);
    const double w = hasWeight ? (double)genWeight : 1.0;
    sumw  += w;
    sumw2 += w * w;
    if (w < wmin) wmin = w;
    if (w > wmax) wmax = w;
    if (w < 0.0)  ++nNeg;
  }

  r.ok        = true;
  r.hasWeight = hasWeight;
  r.nraw      = n;
  r.sumw      = sumw;
  r.sumw2     = sumw2;
  r.wmin      = (n > 0 ? wmin : 0.0);
  r.wmax      = (n > 0 ? wmax : 0.0);
  r.nNeg      = nNeg;

  std::cout << "       " << std::left << std::setw(10) << label
            << "  N_raw=" << n
            << "  N_gen=" << std::fixed << std::setprecision(2) << sumw
            << "  <w>=" << (n > 0 ? sumw / n : 0.0)
            << "  w in [" << r.wmin << ", " << r.wmax << "]"
            << "  Nneg=" << nNeg << "\n";

  f->Close();
  delete f;
  return r;
}

void WriteTxt(const std::vector<NgenResult> &res)
{
  const std::string outDir  = "./output/";
  const std::string outFile = outDir + "ngen.txt";
  std::ofstream os(outFile.c_str(), std::ios::out | std::ios::trunc);
  if (!os.is_open())
  {
    std::cerr << "[ERR] cannot open " << outFile << "\n";
    return;
  }

  os << "# N_gen per MC sample = Sum of hiEvtAnalyzer/HiTree::weight over ALL events (no selection).\n";
  os << "# N_gen is the cross-section normalization denominator. err(N_gen) = sqrt(Sum w^2).\n";
  os << "# Columns: label  N_raw  N_gen  err(N_gen)  <w>  w_min  w_max  N_neg  file\n";
  os << "#-------------------------------------------------------------------------------\n";
  os << std::left << std::setw(10) << "# label"
     << std::right
     << std::setw(12) << "N_raw"
     << std::setw(16) << "N_gen"
     << std::setw(14) << "err(N_gen)"
     << std::setw(10) << "<w>"
     << std::setw(10) << "w_min"
     << std::setw(10) << "w_max"
     << std::setw(10) << "N_neg"
     << "  file\n";

  for (const auto &r : res)
  {
    const double err  = std::sqrt(r.sumw2);
    const double mean = (r.nraw > 0 ? r.sumw / r.nraw : 0.0);
    os << std::left << std::setw(10) << r.label
       << std::right << std::fixed
       << std::setw(12) << r.nraw
       << std::setprecision(2) << std::setw(16) << r.sumw
       << std::setprecision(2) << std::setw(14) << err
       << std::setprecision(4) << std::setw(10) << mean
       << std::setprecision(4) << std::setw(10) << r.wmin
       << std::setprecision(4) << std::setw(10) << r.wmax
       << std::setw(10) << r.nNeg
       << "  " << r.fname
       << (r.ok ? "" : "   [FAILED TO COUNT]")
       << (r.hasWeight ? "" : "   [no weight branch -> unweighted]")
       << "\n";
  }
  os.close();
  std::cout << "[INFO] wrote " << outFile << "\n";
}

void WriteRoot(const std::vector<NgenResult> &res)
{
  gSystem->mkdir("rootfile", kTRUE);
  TFile *fout = new TFile("./rootfile/ngen.root", "RECREATE");

  const int nb = (int)res.size();
  TH1D *h_ngen = new TH1D("h_ngen", "N_{gen}=#Sigma w per MC sample; sample; N_{gen}", nb, 0, nb);
  TH1D *h_nraw = new TH1D("h_nraw", "Raw event count per MC sample; sample; N_{raw}", nb, 0, nb);
  h_ngen->Sumw2();
  h_nraw->Sumw2();

  for (int i = 0; i < nb; ++i)
  {
    const auto &r = res[i];
    const int   b = i + 1; // ROOT bins are 1-based

    h_ngen->GetXaxis()->SetBinLabel(b, r.label.c_str());
    h_ngen->SetBinContent(b, r.sumw);
    h_ngen->SetBinError(b, std::sqrt(r.sumw2)); // err(N_gen) = sqrt(Sum w^2)

    h_nraw->GetXaxis()->SetBinLabel(b, r.label.c_str());
    h_nraw->SetBinContent(b, (double)r.nraw);
    h_nraw->SetBinError(b, std::sqrt((double)r.nraw)); // raw Poisson on the count

    // Per-sample scalars for direct Get() by name downstream.
    TParameter<double> pNgen((std::string("ngen_") + r.label).c_str(), r.sumw);
    TParameter<double> pNraw((std::string("nraw_") + r.label).c_str(), (double)r.nraw);
    pNgen.Write();
    pNraw.Write();
  }

  h_ngen->Write();
  h_nraw->Write();
  fout->Close();
  delete fout;
  std::cout << "[INFO] wrote ./rootfile/ngen.root\n";
}

} // anonymous namespace

int count_ngen()
{
  // The 9 distinct physical MC files, enumerated through the single source of
  // truth (ResolveMCSample in skim_common.h). flavour only splits DY/W; tau is
  // flavour-independent. Dedupe by resolved filename so nothing is counted twice.
  const std::vector<McJob> jobs = {
      {kDY,    "mu"},  {kWp,    "mu"},  {kWm,    "mu"},   // muon-channel MC
      {kDY,    "ele"}, {kWp,    "ele"}, {kWm,    "ele"},  // electron-channel MC
      {kDYtau, "mu"},  {kWptau, "mu"},  {kWmtau, "mu"},   // tau backgrounds (shared)
  };

  std::cout << "============================================================\n";
  std::cout << "N_gen = Sum(HiTree::weight) over ALL events, per MC sample\n";
  std::cout << "============================================================\n";

  std::vector<NgenResult> res;
  std::vector<std::string> seen; // resolved filenames already counted

  for (const auto &j : jobs)
  {
    const auto info = ResolveMCSample(j.sample, j.flavour);
    if (info.fname.empty()) // kData (never in this list) or unmapped
      continue;

    // Dedupe: same physical file reached via two (sample, flavour) routes.
    bool dup = false;
    for (const auto &s : seen)
      if (s == info.fname) { dup = true; break; }
    if (dup)
      continue;
    seen.push_back(info.fname);

    res.push_back(CountOne(LabelFromFname(info.fname), info.fname));
  }

  gSystem->mkdir("output", kTRUE);
  WriteTxt(res);
  WriteRoot(res);

  std::cout << "------------------------------------------------------------\n";
  std::cout << "[DONE] counted " << res.size() << " MC files. "
            << "See output/ngen.txt and rootfile/ngen.root\n";
  std::cout << "       (TODO check: if every <w> == 1 and N_neg == 0, the gen\n";
  std::cout << "        weight is a no-op; otherwise it matters for the xsec.)\n";
  return 0;
}
