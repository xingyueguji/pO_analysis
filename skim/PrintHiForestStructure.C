// File: PrintHiForestStructure.C
// Usage:
//   root -l -q 'PrintHiForestStructure.C("pO_2025.root")'
//   root -l -q 'PrintHiForestStructure.C("pO_2025.root", false, "structure.txt")'
//   root -l -q 'PrintHiForestStructure.C("pO_2025.root", true,  "structure_with_leaves.txt")'

#include <iostream>
#include <cstdio>
#include <string>

#include "TFile.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TObjArray.h"
#include "TClass.h"

static void indent(int n) {
  for (int i = 0; i < n; ++i) std::cout << "  ";
}

static std::string GetBranchTypeSummary(TBranch* br)
{
  if (!br) return "unknown";

  // 1) If this is an object branch (e.g. std::vector<float>), ROOT often fills GetClassName()
  const char* cls = br->GetClassName();
  if (cls && std::string(cls).size() > 0) {
    return std::string("brClass=") + cls;
  }

  // 2) Otherwise fall back to leaf type(s)
  TObjArray* leaves = br->GetListOfLeaves();
  if (!leaves || leaves->GetEntries() == 0) return "noLeaves";

  // If there is exactly one leaf and leaf name == branch name, print just that type.
  if (leaves->GetEntries() == 1) {
    TLeaf* lf = (TLeaf*)leaves->At(0);
    if (lf) {
      const char* t = lf->GetTypeName();
      if (t && std::string(t).size() > 0) return std::string("leafType=") + t;
    }
  }

  // Otherwise, print all leaves (split branches / structs)
  std::string out = "leafTypes=[";
  for (int j = 0; j < leaves->GetEntries(); ++j) {
    TLeaf* lf = (TLeaf*)leaves->At(j);
    if (!lf) continue;
    out += lf->GetName();
    out += ":";
    out += (lf->GetTypeName() ? lf->GetTypeName() : "unknown");
    if (j != leaves->GetEntries() - 1) out += ", ";
  }
  out += "]";
  return out;
}

static void PrintTreeBranches(TTree *t, int level, bool printLeaves, bool printBranchType)
{
  if (!t) return;

  TObjArray *branches = t->GetListOfBranches();
  if (!branches) return;

  for (int i = 0; i < branches->GetEntries(); ++i) {
    TBranch *br = (TBranch*)branches->At(i);
    if (!br) continue;

    indent(level);
    std::cout << "  * " << br->GetName();

    // NEW: branch type summary
    if (printBranchType) {
      std::cout << "  {" << GetBranchTypeSummary(br) << "}";
    }

    // Optional: leaf list (your original behavior)
    if (printLeaves) {
      TObjArray *leaves = br->GetListOfLeaves();
      if (leaves && leaves->GetEntries() > 0) {
        std::cout << "  [";
        for (int j = 0; j < leaves->GetEntries(); ++j) {
          TLeaf *lf = (TLeaf*)leaves->At(j);
          if (!lf) continue;
          std::cout << lf->GetName() << ":" << lf->GetTypeName();
          if (j != leaves->GetEntries() - 1) std::cout << ", ";
        }
        std::cout << "]";
      }
    }

    std::cout << "\n";
  }
}

static void PrintDirRecursive(TDirectory *dir, int level, bool printLeaves, bool printBranchType)
{
  if (!dir) return;

  TIter next(dir->GetListOfKeys());
  TKey *key = nullptr;

  while ((key = (TKey*)next())) {
    TObject *obj = key->ReadObj();
    if (!obj) continue;

    indent(level);
    std::cout << "- " << obj->GetName()
              << "  [" << obj->ClassName() << "]";

    if (obj->InheritsFrom(TTree::Class())) {
      TTree *t = (TTree*)obj;
      std::cout << "  entries=" << t->GetEntries() << "\n";
      PrintTreeBranches(t, level + 1, printLeaves, printBranchType);
    } else {
      std::cout << "\n";
    }

    if (obj->InheritsFrom(TDirectory::Class())) {
      PrintDirRecursive((TDirectory*)obj, level + 1, printLeaves, printBranchType);
    }
  }
}

void PrintHiForestStructure(const char *filename = "root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_mc_version_2_2025.root",
                            bool printLeaves = false,
                            const char *outFile = "./mc_version_2.txt",
                            bool printBranchType = true)
{
  // Redirect stdout to file if requested
  FILE *old_stdout = nullptr;
  if (outFile && std::string(outFile).size() > 0) {
    old_stdout = stdout;
    freopen(outFile, "w", stdout);
  }

  std::cout << "File: " << filename << "\n";

  TFile *f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    std::cerr << "ERROR: cannot open file: " << filename << "\n";
    if (outFile) {
      fflush(stdout);
      fclose(stdout);
      stdout = old_stdout;
    }
    return;
  }

  PrintDirRecursive(f, 0, printLeaves, printBranchType);

  f->Close();
  delete f;

  // Restore stdout
  if (outFile) {
    fflush(stdout);
    fclose(stdout);
    stdout = old_stdout;
  }
}