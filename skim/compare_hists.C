// skim/compare_hists.C -- bit-identity check of the histograms in two ROOT files.
//
// For every TH1-derived key of file A (recursing into TDirectories, e.g. the
// per-region dirs of a Combine input), finds the same path in file B and
// compares axes, every cell content (under/overflow included), the Sumw2
// arrays (when both carry them) and the entry count with exact equality.
// Keys only in B are reported as NEW (e.g. the LHE member twins after the
// 2026-09-07 re-skim, or the <process>_<syst>Up/Down Combine templates). Use
// it after any change that is supposed to leave the existing histograms
// untouched:
//
//   root -l -b -q 'compare_hists.C+("rootfile_pre_lhe/X_hist.root","rootfile/X_hist.root")'
//   ./compare_reskim.sh rootfile_pre_lhe        # all files of a backup dir
//
// Returns the number of DIFF + MISSING keys (0 = the old content is identical).
#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TKey.h"
#include "TList.h"
#include "TObject.h"
#include "TArrayD.h"
#include <cstdio>
#include <set>
#include <string>

namespace {

bool SameAxis(const TAxis *a, const TAxis *b)
{
  if (a->GetNbins() != b->GetNbins()) return false;
  for (int i = 0; i <= a->GetNbins(); ++i)
    if (a->GetBinLowEdge(i + 1) != b->GetBinLowEdge(i + 1)) return false;
  return true;
}

// "" when identical, otherwise the first reason
std::string Compare(const TH1 *a, const TH1 *b)
{
  if (std::string(a->ClassName()) != b->ClassName()) return "class";
  if (!SameAxis(a->GetXaxis(), b->GetXaxis())) return "x axis";
  if (a->GetDimension() >= 2 && !SameAxis(a->GetYaxis(), b->GetYaxis())) return "y axis";
  if (a->GetNcells() != b->GetNcells()) return "ncells";
  for (int i = 0; i < a->GetNcells(); ++i)
    if (a->GetBinContent(i) != b->GetBinContent(i)) return "content";
  const bool sa = a->GetSumw2N() > 0, sb = b->GetSumw2N() > 0;
  if (sa != sb) return "sumw2 presence";
  if (sa)
    for (int i = 0; i < a->GetNcells(); ++i)
      if (a->GetSumw2()->At(i) != b->GetSumw2()->At(i)) return "sumw2";
  if (a->GetEntries() != b->GetEntries()) return "entries";
  return "";
}

struct Counts { int ident = 0, diff = 0, missing = 0, fresh = 0; };

// Walk directory A (and its subdirectories); look every histogram up in B.
void CompareDir(TDirectory *A, TDirectory *B, const std::string &prefix, Counts &c, bool listNew)
{
  std::set<std::string> seen;
  TIter itA(A->GetListOfKeys());
  while (TKey *k = (TKey *)itA())
  {
    if (!seen.insert(k->GetName()).second) continue; // one cycle per name
    const std::string path = prefix + k->GetName();
    if (std::string(k->GetClassName()).rfind("TDirectory", 0) == 0)
    {
      TDirectory *dA = (TDirectory *)k->ReadObj();
      TDirectory *dB = B ? B->GetDirectory(k->GetName()) : nullptr;
      if (!dB) { ++c.missing; printf("  MISSING  %s/ (whole directory)\n", path.c_str()); continue; }
      CompareDir(dA, dB, path + "/", c, listNew);
      continue;
    }
    TObject *oa = k->ReadObj();
    if (!oa || !oa->InheritsFrom("TH1")) continue;
    TH1 *ha = (TH1 *)oa;
    TH1 *hb = B ? (TH1 *)B->Get(k->GetName()) : nullptr;
    if (!hb) { ++c.missing; printf("  MISSING  %s\n", path.c_str()); continue; }
    const std::string why = Compare(ha, hb);
    if (why.empty()) ++c.ident;
    else { ++c.diff; printf("  DIFF     %s  (%s)\n", path.c_str(), why.c_str()); }
  }
  if (!B) return;
  std::set<std::string> seenB;
  TIter itB(B->GetListOfKeys());
  while (TKey *k = (TKey *)itB())
  {
    if (!seenB.insert(k->GetName()).second) continue;
    if (seen.count(k->GetName())) continue;
    if (std::string(k->GetClassName()).rfind("TDirectory", 0) == 0)
    {
      // a whole new directory: count its histograms as NEW
      Counts sub;
      CompareDir((TDirectory *)k->ReadObj(), nullptr, prefix + k->GetName() + "/", sub, listNew);
      c.fresh += sub.ident + sub.diff + sub.missing; // all of them are "missing in A"
      continue;
    }
    ++c.fresh;
    if (listNew) printf("  NEW      %s%s (%s)\n", prefix.c_str(), k->GetName(), k->GetClassName());
  }
}

} // namespace

int compare_hists(const char *fileA, const char *fileB, bool listNew = false)
{
  TH1::AddDirectory(kFALSE);
  TFile *A = TFile::Open(fileA, "READ");
  TFile *B = TFile::Open(fileB, "READ");
  if (!A || A->IsZombie() || !B || B->IsZombie())
  {
    printf("[ERR] compare_hists: cannot open %s or %s\n", fileA, fileB);
    return 999;
  }
  Counts c;
  CompareDir(A, B, "", c, listNew);
  printf("[SUMMARY] %s vs %s: IDENTICAL %d  DIFF %d  MISSING %d  NEW %d\n",
         fileA, fileB, c.ident, c.diff, c.missing, c.fresh);
  A->Close(); B->Close();
  return c.diff + c.missing;
}
