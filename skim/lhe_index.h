// skim/lhe_index.h -- the decoded layout of hiEvtAnalyzer/HiTree::ttbar_w
// (217 LHE reweighting weights per MC event) and the STORED-MEMBER families
// the skim fills as per-variation histogram twins.
//
// SINGLE SOURCE for every macro that touches ttbar_w:
//   skim/lhe_weights.C   inclusive decoding/diagnostics (the index tables,
//                        Hessian() and SetLabel() below were moved out of its
//                        anonymous namespace on 2026-09-07 -- byte-identical
//                        lhe_weights.txt after the move)
//   skim/skim.C          fills the <hist>_epps21 / _scale / _alphas TH2D twins
//                        of every fit-template histogram (MC only)
//   skim/lhe_updown.py   combines the _epps21 members per bin with LHAPDF's
//                        PDFSet.uncertainty() into <hist>_nPDFUp/Down
//
// The full provenance of the 217 entries (producing script:line per index)
// is in skim/output/lhe_weights.txt; the physics of the two blocks and the
// agreed combination recipe are in CLAUDE.md (Stage-1, lhe_weights.C bullet).
// In short: f_O = R(nPDFerrSet) x f_p^LHAPDF -- idx 45-102 vary f_p (CT18ANLO
// eigenvectors, both beams, R central), idx 111-158 vary R (EPPS21 nuclear
// error sets), idx 159-216 vary R along the CT18A-BASELINE directions only, so
// the coherent EPPS21 baseline member (paper Eq. 40) is the PRODUCT of the
// per-event ratios idx 44+k and 158+k' -- built per event in
// ComputeMemberWeights() below, never from histograms.
//
// STORED-MEMBER FAMILIES (the y axis of each TH2D twin = member index):
//   kEpps21  107 members = the LHAPDF set EPPS21nlo_CT18Anlo_O16, member for
//            member, so a bin's column feeds lhapdf.PDFSet.uncertainty() 1:1:
//              0      nominal                       = ttbar_w[0]
//              1-48   nuclear error sets 2..49      = ttbar_w[111..158]
//              49-106 baseline pairs, k = m-48:       (ttbar_w[44+k]/w0) x (ttbar_w[158+k']/w0)
//                     k' = k (kStraight) or the other member of k's pair
//                     (kCrossed = coherent if EPPS21.f even=S- / CTEQ odd="+")
//   kScale   9 members = ttbar_w[0..8], the (muR,muF) grid (inner loop muF)
//   kAlphaS  5 members = ttbar_w[0], ttbar_w[103..106] (CT18ANLO alpha_s
//            0.116 / 0.117 / 0.119 / 0.120)
// Per-event member weight = w * (ttbar_w[a]/ttbar_w[0]) [* ttbar_w[b]/ttbar_w[0]]
// with w the skim's nominal event weight -- so member 0 == the nominal
// histogram bit for bit, and any SF later folded into w propagates.
#ifndef PO_LHE_INDEX_H
#define PO_LHE_INDEX_H

#include "TAxis.h"
#include "TH1.h"
#include "TH2D.h"
#include "TString.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace pOLhe
{

// ============================================================================
// Part 1 -- the decoded ttbar_w layout (moved verbatim from lhe_weights.C)
// ============================================================================
inline constexpr int kNW = 217; // expected vector length (checked at run time)

// ---- the decoded layout (see header comment) -------------------------------
enum BlockKind { kCentral, kEnvelope, kList, kHessian };
struct Block
{
  int         lo, hi;
  BlockKind   kind;
  const char *name;   // one-line description (block summaries)
  const char *tag;    // short role, printed on every row of the per-index table
  const char *where;  // which script/line PRODUCES these weights (paths relative to kRefDir)
  const char *why;    // why the production includes them
};

inline const std::vector<Block> kBlocks = {
    {  0,   0, kCentral,  "nominal: id 1001 muR=muF=1 (== HiTree::weight)",
       "nominal (reference)",
       "bin/Powheg/make_rwl.py L43-54: the (1d0,1d0) corner of the scale loop, `lhapdf=CentralPDF` (14600)",
       "the nominal configuration re-evaluated by the reweighter; the reference every other weight is compared to (== HiTree::weight up to float rounding)"},
    {  1,   8, kEnvelope, "QCD scale grid ids 1002-1009, (muR,muF) in {1,2,0.5}^2 minus (1,1)",
       "QCD scale variation",
       "bin/Powheg/make_rwl.py L43-54: `for m_rensc in m_factor: for m_facsc in m_factor`, m_factor = [1, 2, 0.5]",
       "the QCD scale uncertainty: renormalization and factorization scales varied by 2 and 1/2; the envelope is the scale systematic (7-point convention drops the (2,1/2) and (1/2,2) corners = idx 5 and 7)"},
    {  9,  43, kList,     "35 proton-PDF centrals (NNPDF3.x/4.0 + alpha_s, CT18*, MSHT20, PDF4LHC21, HERAPDF, ABMP16)",
       "PDF-set comparison central",
       "bin/Powheg/make_rwl.py L490-524: the 35 n=1 entries of the `PDF_variation1 , hessian` group of the `\"EPPS21\" in Period` branch (L481)",
       "the CMS Run-3 standard list of alternative proton PDFs and alpha_s(M_Z) values, reduced to their central members: for comparing/choosing the proton PDF and quantifying the alpha_s dependence -- a comparison set, NOT an uncertainty band (mostly NNLO sets in an NLO calculation, hence the +2..5% offsets)"},
    { 44,  44, kCentral,  "CT18ANLO member 0 (== the nominal LHAPDF set, both beams)",
       "nominal proton PDF central",
       "bin/Powheg/make_rwl.py L526: `[4400, 14600, 'CT18ANLO', 59]` member 0 (comment L525: \"Reference proton PDF of EPPS21nlo nuclear PDF\")",
       "the central member of the nominal proton PDF (CT18ANLO = the baseline EPPS21 is fitted on); reproduces the nominal exactly"},
    { 45, 102, kHessian,  "CT18ANLO eigenvectors 1-58 on both beams, R central (29 pairs, 90% CL)",
       "proton-PDF uncertainty (CT18ANLO eig.)",
       "bin/Powheg/make_rwl.py L526: `[4400, 14600, 'CT18ANLO', 59]` members 1-58",
       "the proton-PDF uncertainty: the 29 Hessian eigenvector pairs of CT18ANLO (90% CL), applied to the LHAPDF set of both beams while the nuclear R-factor stays at its central set"},
    {103, 106, kList,     "CT18ANLO alpha_s 0.116 / 0.117 / 0.119 / 0.120",
       "proton alpha_s variation",
       "bin/Powheg/make_rwl.py L527-530: `[4500..4503, 14666/14667/14669/14670, 'CT18ANLO_as_0116..0120', 1]`",
       "the alpha_s(M_Z) uncertainty of the nominal proton PDF: CT18ANLO refits at 0.116/0.117/0.119/0.120 (nominal 0.118)"},
    {107, 109, kList,     "replica-group centrals: NNPDF31 mc, NNPDF40 pdfas, NNPDF40 pch",
       "PDF-set comparison central (replica)",
       "bin/Powheg/make_rwl.py L535-537: the 3 n=1 entries of the `PDF_variation2 , replica` group (written after the hessian group by `sorted(pdf_sets.items())`, L548)",
       "the same comparison list, Monte-Carlo-replica representations of the NNPDF fits (central members only); 107 duplicates idx 9 because the replica mean equals the mc_hessian central"},
    {110, 110, kCentral,  "EPPS21 R-factor set 1 = central (== nominal), id 9001",
       "nominal nPDF R central",
       "bin/Powheg/runcmsgrid_powheg.sh L254-262: `for iset in {1..107}`, `<weight id='$((9000+iset))'> nPDFerrSet=iset`, iset=1",
       "the central EPPS21 nuclear modification (EPPS21.f pset 1) on the oxygen beam; reproduces the nominal exactly"},
    {111, 158, kHessian,  "EPPS21 R-factor sets 2-49: nuclear eigen-directions (24 pairs, 90% CL)",
       "NUCLEAR-PDF uncertainty (EPPS21 R)",
       "bin/Powheg/runcmsgrid_powheg.sh L254-262 with iset = 2..49; R applied in patches/EPPS21/lhapdf6if_nPDF.patch L38-40 (evolvePDFhi -> EPPS21.f)",
       "the nuclear-PDF uncertainty: EPPS21 error sets S-+1..S-+24 = the 24 nuclear eigen-directions, plus and minus (90% CL) -- the systematic this analysis needs"},
    {159, 216, kHessian,  "EPPS21 R-factor sets 50-107: CT18A-baseline directions, R only (29 pairs, 90% CL)",
       "R baseline dependence (not a PDF var.)",
       "bin/Powheg/runcmsgrid_powheg.sh L254-262 with iset = 50..107",
       "how the nuclear modification R responds to the 29 CT18A-baseline eigen-directions; ONLY R changes here, the LHAPDF baseline stays central -- the baseline PDF variation itself is idx 45-102 (combine k with 44+k as the product of the two ratios for the coherent EPPS21 baseline variation)"},
};
inline const Block &BlockOf(int i)
{
  for (const Block &b : kBlocks)
    if (i >= b.lo && i <= b.hi) return b;
  return kBlocks.back();
}

// ---- the <initrwgt> map (make_rwl.py "EPPS21" branch + the runcmsgrid EPPS21 loop) --
// One entry per contiguous run of weights from the same source: first ttbar_w
// index, number of members, first LHE weight id, first LHAPDF id (-1 = the
// scale grid, -2 = EPPS21 R-factor error sets selected by nPDFerrSet, no LHAPDF id).
struct SetRun
{
  int         lo, n, wid0, lha0;
  const char *name;
  const char *src;  // the line of the reference script that writes this run (relative to kRefDir)
};
inline const std::vector<SetRun> kSets = {
    {  0, 9, 1001,     -1, "scale",                                  "bin/Powheg/make_rwl.py:52-54 scale loop"},
    {  9, 1, 2000, 325300, "NNPDF31_nnlo_as_0118_mc_hessian_pdfas",  "bin/Powheg/make_rwl.py:490"},
    { 10, 1, 2200, 306000, "NNPDF31_nnlo_hessian_pdfas",             "bin/Powheg/make_rwl.py:491"},
    { 11, 1, 2201, 322500, "NNPDF31_nnlo_as_0108",                   "bin/Powheg/make_rwl.py:492"},
    { 12, 1, 2202, 322700, "NNPDF31_nnlo_as_0110",                   "bin/Powheg/make_rwl.py:493"},
    { 13, 1, 2203, 322900, "NNPDF31_nnlo_as_0112",                   "bin/Powheg/make_rwl.py:494"},
    { 14, 1, 2204, 323100, "NNPDF31_nnlo_as_0114",                   "bin/Powheg/make_rwl.py:495"},
    { 15, 1, 2205, 323300, "NNPDF31_nnlo_as_0117",                   "bin/Powheg/make_rwl.py:496"},
    { 16, 1, 2206, 323500, "NNPDF31_nnlo_as_0119",                   "bin/Powheg/make_rwl.py:497"},
    { 17, 1, 2207, 323700, "NNPDF31_nnlo_as_0122",                   "bin/Powheg/make_rwl.py:498"},
    { 18, 1, 2208, 323900, "NNPDF31_nnlo_as_0124",                   "bin/Powheg/make_rwl.py:499"},
    { 19, 1, 2300, 305800, "NNPDF31_nlo_hessian_pdfas",              "bin/Powheg/make_rwl.py:500"},
    { 20, 1, 2500, 303200, "NNPDF30_nnlo_as_0118_hessian",           "bin/Powheg/make_rwl.py:501"},
    { 21, 1, 2501, 292200, "NNPDF30_nlo_nf_5_pdfas",                 "bin/Powheg/make_rwl.py:502"},
    { 22, 1, 2600, 331600, "NNPDF40_nnlo_hessian_pdfas",             "bin/Powheg/make_rwl.py:503"},
    { 23, 1, 2700, 332700, "NNPDF40_nnlo_as_01160",                  "bin/Powheg/make_rwl.py:504"},
    { 24, 1, 2701, 332900, "NNPDF40_nnlo_as_01170",                  "bin/Powheg/make_rwl.py:505"},
    { 25, 1, 2702, 333100, "NNPDF40_nnlo_as_01175",                  "bin/Powheg/make_rwl.py:506"},
    { 26, 1, 2703, 333300, "NNPDF40_nnlo_as_01185",                  "bin/Powheg/make_rwl.py:507"},
    { 27, 1, 2704, 333500, "NNPDF40_nnlo_as_01190",                  "bin/Powheg/make_rwl.py:508"},
    { 28, 1, 2705, 333700, "NNPDF40_nnlo_as_01200",                  "bin/Powheg/make_rwl.py:509"},
    { 29, 1, 2800, 332300, "NNPDF40_nlo_pch_as_01180",               "bin/Powheg/make_rwl.py:510"},
    { 30, 1, 4000,  14000, "CT18NNLO",                               "bin/Powheg/make_rwl.py:511"},
    { 31, 1, 4100,  14066, "CT18NNLO_as_0116",                       "bin/Powheg/make_rwl.py:512"},
    { 32, 1, 4101,  14067, "CT18NNLO_as_0117",                       "bin/Powheg/make_rwl.py:513"},
    { 33, 1, 4102,  14069, "CT18NNLO_as_0119",                       "bin/Powheg/make_rwl.py:514"},
    { 34, 1, 4103,  14070, "CT18NNLO_as_0120",                       "bin/Powheg/make_rwl.py:515"},
    { 35, 1, 4200,  14100, "CT18ZNNLO",                              "bin/Powheg/make_rwl.py:516"},
    { 36, 1, 4300,  14200, "CT18ANNLO",                              "bin/Powheg/make_rwl.py:517"},
    { 37, 1, 4301,  14300, "CT18XNNLO",                              "bin/Powheg/make_rwl.py:518"},
    { 38, 1, 5000,  27400, "MSHT20nnlo_as118",                       "bin/Powheg/make_rwl.py:519"},
    { 39, 1, 5100,  27500, "MSHT20nnlo_as_smallrange",               "bin/Powheg/make_rwl.py:520"},
    { 40, 1, 5101,  27550, "MSHT20nnlo_as_largerange",               "bin/Powheg/make_rwl.py:521"},
    { 41, 1, 6000,  93300, "PDF4LHC21_40_pdfas",                     "bin/Powheg/make_rwl.py:522"},
    { 42, 1, 7000,  61200, "HERAPDF20_NNLO_EIG",                     "bin/Powheg/make_rwl.py:523"},
    { 43, 1, 8000,  42780, "ABMP16als118_5_nnlo",                    "bin/Powheg/make_rwl.py:524"},
    { 44, 59, 4400, 14600, "CT18ANLO",                               "bin/Powheg/make_rwl.py:526"},
    {103, 1, 4500,  14666, "CT18ANLO_as_0116",                       "bin/Powheg/make_rwl.py:527"},
    {104, 1, 4501,  14667, "CT18ANLO_as_0117",                       "bin/Powheg/make_rwl.py:528"},
    {105, 1, 4502,  14669, "CT18ANLO_as_0119",                       "bin/Powheg/make_rwl.py:529"},
    {106, 1, 4503,  14670, "CT18ANLO_as_0120",                       "bin/Powheg/make_rwl.py:530"},
    {107, 1, 3000, 316200, "NNPDF31_nnlo_as_0118_mc",                "bin/Powheg/make_rwl.py:535"},
    {108, 1, 3200, 331300, "NNPDF40_nnlo_pdfas",                     "bin/Powheg/make_rwl.py:536"},
    {109, 1, 3400, 332100, "NNPDF40_nnlo_pch_as_01180",              "bin/Powheg/make_rwl.py:537"},
    {110, 107, 9001,     -2, "EPPS21 R-factor (O16)",                 "bin/Powheg/runcmsgrid_powheg.sh:254-262 nPDFerrSet loop"},
};
inline const SetRun &SetOf(int i)
{
  for (const SetRun &s : kSets)
    if (i >= s.lo && i < s.lo + s.n) return s;
  return kSets.back();
}
inline const char *const kScaleLabel[9] = {"muR=1 muF=1", "muR=1 muF=2", "muR=1 muF=0.5", "muR=2 muF=1", "muR=2 muF=2",
                              "muR=2 muF=0.5", "muR=0.5 muF=1", "muR=0.5 muF=2", "muR=0.5 muF=0.5"};

// "CT18ANLO m12  (id 4412, lhapdf 14612)" / "muR=2 muF=0.5  (id 1006)" /
// "EPPS21 R-factor (O16) nPDFerrSet=17  (id 9017)"
inline std::string SetLabel(int i)
{
  for (const SetRun &s : kSets)
    if (i >= s.lo && i < s.lo + s.n)
    {
      const int k = i - s.lo;
      if (s.lha0 == -1) return std::string(kScaleLabel[k]) + Form("  (id %d)", s.wid0 + k);
      if (s.lha0 == -2) return std::string(s.name) + Form(" nPDFerrSet=%d  (id %d)", k + 1, s.wid0 + k);
      return std::string(s.name) + Form(" m%d  (id %d, lhapdf %d)", k, s.wid0 + k, s.lha0 + k);
    }
  return "(beyond the 217-entry map)";
}

// ---- Hessian combination (moved verbatim from lhe_weights.C) -----------------
struct HessRes { double sym, up, dn; };

// Hessian sums over consecutive pairs (lo,lo+1),(lo+2,lo+3),... of the
// relative deviations d_i = x_i/x_0 - 1, where x is any functional (sigma,
// ratio, asymmetry) evaluated per member; passed in as a function of index.
template <class F>
HessRes Hessian(const Block &b, F dev)
{
  double hs = 0, up = 0, dn = 0;
  for (int i = b.lo; i + 1 <= b.hi; i += 2)
  {
    const double a = dev(i), c = dev(i + 1);
    hs += (a - c) * (a - c);
    const double u = std::max(std::max(a, c), 0.0);
    const double d = std::min(std::min(a, c), 0.0);
    up += u * u; dn += d * d;
  }
  return {0.5 * std::sqrt(hs), std::sqrt(up), std::sqrt(dn)};
}

// ============================================================================
// Part 2 -- stored-member families for the skim twins (2026-09-07)
// ============================================================================
enum Family { kEpps21 = 0, kScale = 1, kAlphaS = 2, kNFamilies = 3 };
inline constexpr int         kNMembers[kNFamilies]    = {107, 9, 5};
inline constexpr int         kNMaxMembers             = 107;
inline constexpr const char *kFamilySuffix[kNFamilies] = {"_epps21", "_scale", "_alphas"};
inline constexpr const char *kEpps21SetName           = "EPPS21nlo_CT18Anlo_O16";

// The Up/Down template families skim/lhe_updown.py writes into every MC skim
// file (<hist>_<name>Up / <hist>_<name>Down) and plotting/mtandmet.C +
// dileptonpeak.C carry into the Combine inputs as <process>_<name>Up/Down.
// MUST match FAMILIES in lhe_updown.py (nPDF = LHAPDF PDFSet.uncertainty() on
// the _epps21 members; qcdScale = per-bin max/min over the _scale members;
// alphaS = the 0.119 / 0.117 _alphas members).
inline constexpr const char *kLheSystNames[] = {"nPDF", "qcdScale", "alphaS"};
inline constexpr int         kNLheSysts      = 3;

// Sign pairing of the 29 baseline products (open verification item, CLAUDE.md):
// kCrossed pairs CT18ANLO member 2i-1 with R set 49+2i (and 2i with 48+2i).
enum BaselinePairing { kStraight = 0, kCrossed = 1 };
inline constexpr BaselinePairing kBaselinePairing = kCrossed;

// One stored member: weight = w * (ttbar_w[a]/w0) * (b < 0 ? 1 : ttbar_w[b]/w0)
struct MemberDef { int a; int b; };

inline MemberDef MemberIndices(Family f, int m)
{
  switch (f)
  {
    case kEpps21:
      if (m <= 0)  return {0, -1};
      if (m <= 48) return {110 + m, -1};                       // nuclear sets 2..49
      {
        const int k  = m - 48;                                  // 1..58 = CT18ANLO member
        const int kp = (kBaselinePairing == kCrossed) ? ((k % 2) ? k + 1 : k - 1) : k;
        return {44 + k, 158 + kp};                              // x R baseline set 49+kp
      }
    case kScale:  return {std::max(0, std::min(m, 8)), -1};    // idx 0..8
    case kAlphaS: return {m <= 0 ? 0 : 102 + std::min(m, 4), -1}; // idx 0, 103..106
    default:      return {0, -1};
  }
}

// Precomputed per family (called ~120 times per event otherwise).
inline const std::vector<MemberDef> &MemberTable(Family f)
{
  static std::vector<MemberDef> tab[kNFamilies];
  if (tab[f].empty())
    for (int m = 0; m < kNMembers[f]; ++m) tab[f].push_back(MemberIndices(f, m));
  return tab[f];
}

inline std::string MemberName(Family f, int m)
{
  const MemberDef d = MemberIndices(f, m);
  if (d.b < 0) return SetLabel(d.a);
  return SetLabel(d.a) + "  x  " + SetLabel(d.b);
}

// Stamped into every twin's title (no ';' -- ROOT would split it into axis titles).
inline std::string FamilyTitle(Family f)
{
  switch (f)
  {
    case kEpps21:
      return std::string(kEpps21SetName) + " members: 0 nominal, 1-48 = ttbar_w[111..158] (nuclear sets 2-49),"
             " 49-106 = (ttbar_w[44+k]/w0)x(ttbar_w[158+k']/w0), k=1..58, "
             + (kBaselinePairing == kCrossed ? "k'=k odd?k+1:k-1 (crossed pairing)" : "k'=k (straight pairing)");
    case kScale:
      return "QCD scale members = ttbar_w[0..8], (muR,muF): 0=(1,1) 1=(1,2) 2=(1,0.5) 3=(2,1) 4=(2,2)"
             " 5=(2,0.5) 6=(0.5,1) 7=(0.5,2) 8=(0.5,0.5)";
    case kAlphaS:
      return "alpha_s members: 0 nominal (0.118), 1-4 = ttbar_w[103..106] = CT18ANLO alpha_s 0.116, 0.117, 0.119, 0.120";
    default: return "";
  }
}

// The per-event member weights of all families (computed once per event).
struct MemberWeights { double w[kNFamilies][kNMaxMembers]; };

// false (and nothing usable in `out`) when the vector is absent / not 217 long
// (warned once per job) or when ttbar_w[0] == 0 (the nominal adds 0 anyway).
inline bool ComputeMemberWeights(double w, const std::vector<float> *ww, MemberWeights &out,
                                 bool &warnedOnce, const char *who)
{
  if (!ww || (int)ww->size() != kNW)
  {
    if (!warnedOnce)
    {
      std::cout << "[WARN] " << who << ": ttbar_w missing or has " << (ww ? ww->size() : 0)
                << " entries (expected " << kNW << ") -- LHE member twins not filled for such events\n";
      warnedOnce = true;
    }
    return false;
  }
  const double w0 = (*ww)[0];
  if (w0 == 0.0) return false;
  for (int f = 0; f < kNFamilies; ++f)
  {
    const std::vector<MemberDef> &tab = MemberTable((Family)f);
    for (int m = 0; m < kNMembers[f]; ++m)
    {
      double r = (*ww)[tab[m].a] / w0;
      if (tab[m].b >= 0) r *= (*ww)[tab[m].b] / w0;
      out.w[f][m] = w * r;
    }
  }
  return true;
}

// One twin: x axis copied from the nominal (fixed or variable bins), y = member index.
inline TH2D *BookTwin(const TH1 *nom, Family f)
{
  const int   n     = kNMembers[f];
  const TAxis *ax   = nom->GetXaxis();
  const std::string name  = std::string(nom->GetName()) + kFamilySuffix[f];
  const std::string title = std::string(nom->GetTitle()) + " [" + FamilyTitle(f) + "];"
                            + ax->GetTitle() + ";LHE member";
  TH2D *h = (ax->GetXbins()->GetSize() > 0)
              ? new TH2D(name.c_str(), title.c_str(), ax->GetNbins(), ax->GetXbins()->GetArray(), n, -0.5, n - 0.5)
              : new TH2D(name.c_str(), title.c_str(), ax->GetNbins(), ax->GetXmin(), ax->GetXmax(), n, -0.5, n - 0.5);
  if (h->GetSumw2N() == 0) h->Sumw2();
  return h;
}

// The three twins of one fit-template histogram.
struct Twins { TH2D *h[kNFamilies] = {nullptr, nullptr, nullptr}; };

inline Twins BookTwins(const TH1 *nom, std::vector<TH2D *> &registry)
{
  Twins t;
  for (int f = 0; f < kNFamilies; ++f)
  {
    t.h[f] = BookTwin(nom, (Family)f);
    registry.push_back(t.h[f]);
  }
  return t;
}

// Member m of family f gets weight mw.w[f][m] at x (bin centre m = bin m+1).
inline void FillTwins(const Twins &t, double x, const MemberWeights &mw)
{
  for (int f = 0; f < kNFamilies; ++f)
  {
    TH2D *h = t.h[f];
    if (!h) continue;
    for (int m = 0; m < kNMembers[f]; ++m) h->Fill(x, (double)m, mw.w[f][m]);
  }
}

} // namespace pOLhe

#endif // PO_LHE_INDEX_H
