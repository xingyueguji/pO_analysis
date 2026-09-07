// skim/lhe_weights.C
//
// Decode the per-event generator systematic weights stored in
//   hiEvtAnalyzer/HiTree::ttbar_w   (vector<float>, 217 entries/event)
//
// WHAT THE BRANCH IS (HiEvtAnalyzer.cc, pO_HiForest EventAnalysis):
//   ttbar_w[i] = (genInfo->weight() / LHE originalXWGTUP) * LHEEventProduct::weights()[i]
// i.e. the FULL LHE reweighting block (<rwgt>) of the POWHEG event, rescaled so
// that element 0 reproduces HiTree::weight (the nominal, = sigma in pb for this
// production, <w> = sigma). The name "ttbar_w" is historical -- it has nothing
// to do with ttbar; it is the generic "weights for systematics" container.
//
// The ntuple does NOT carry the LHE <initrwgt> header (weight ids / PDF set
// names). What IS certain from the numbers alone: three members are identical
// to the nominal (i = 0, 44, 110 -- a PDF central that reproduces the generation
// PDF), which splits the vector into three REGIONS A = 0-43, B = 44-109,
// C = 110-216; the default FIGURE shows exactly that.
//
// IDENTITIES -- CONFIRMED 2026-09-02 from the producer's generator scripts
// (CMS genproductions bin/Powheg, archive a1a26254; the producer confirmed the
// make_rwl.py `"EPPS21" in Period` branch). TWO mechanisms write the 217:
//  (1) make_rwl.py -> pwg-rwl.dat <initrwgt>, read by POWHEG at generation:
//      the 3x3 scale grid (ids 1001-1009, lhapdf=14600 = run_pwg_condor.py's
//      defaultPDF for an "EPPS21_*" ion) + the "hessian" PDF group + the
//      "replica" PDF group = 110 weights = idx 0-109. `lhapdf=X` swaps the
//      LHAPDF set of BOTH beams (the card must have lhans1 == lhans2 = 14600).
//  (2) runcmsgrid_powheg.sh, AFTER generation: if the card has `nPDFerrSet`,
//      it reruns pwhg_main 107 times with `rwl_add 1`, once per
//      nPDFerrSet = 1..107, appending weight ids 9001-9107 (weightgroup
//      'EPPS21_variation', combine=hessian) = idx 110-216.
// The nuclear PDF itself is a POWHEG source patch (patches/EPPS21/*): for the
// beam with ia >= 16 the LHAPDF proton PDF is multiplied by the EPPS21
// R-factors of error set nPDFerrSet (EPPS21.f + the EPPS21NLOR_16 grid) and
// isospin-averaged (Z/A mixing u<->d, ubar<->dbar). So f_O = R(set) x f_p^LHAPDF
// and the two weight kinds vary the two factors SEPARATELY: (1) varies f_p on
// both beams with R at set 1 (central); (2) varies R with f_p fixed at CT18ANLO
// central. Consequences:
//   * idx 44 (lhapdf=14600) and idx 110 (nPDFerrSet=1) are both exactly the
//     nominal configuration -> the two exact ==nominal markers;
//   * idx 159-216 (EPPS21 sets 50-107 = the CT18A-baseline eigen-directions)
//     move ONLY R's response to the baseline, NOT the baseline PDF itself; the
//     proton-PDF uncertainty proper is idx 45-102 (CT18ANLO eigenvectors on
//     both beams). The coherent EPPS21 baseline variation f_A,k = R_k x f_p,k
//     is, to first order, the PRODUCT of the idx 44+k and idx 158+k per-event
//     ratios (k = 1..58).
// Every testable feature of the map is reproduced by the data: the alpha_s
// scans are monotonic (idx 11-18, 23-28, 30-34, 103-106), 38 == 39 (MSHT20
// as118 vs as_smallrange member 0 are the same grid), 107 == 9 (NNPDF31 replica
// mean == mc_hessian central), 108 ~ 22 (NNPDF40 replica vs hessian central),
// the scale grid obeys the additive/cross-term algebra with idx 1,2,3,6 =
// single-scale variations, and the EPPS21 nuclear/baseline sub-blocks change
// character at 158/159 (EPPS21.f sets 2-49 = 24 nuclear pairs S-+1..S-+24,
// 50-107 = the 29 CT18A-baseline pairs). The table kSets below IS the map:
//
//   idx   0        id 1001  muR=1 muF=1, lhapdf=14600  (== HiTree::weight)
//   idx   1-  8    ids 1002-1009: (muR,muF) in {1,2,0.5}^2 minus (1,1), inner
//                  loop over muF:  1=(1,2) 2=(1,.5) 3=(2,1) 4=(2,2) 5=(2,.5)
//                  6=(.5,1) 7=(.5,2) 8=(.5,.5)
//   idx   9- 43    35 single proton-PDF centrals (hessian group, ids 2000-8000):
//                  NNPDF31 (+alpha_s 0.108..0.124), NNPDF30, NNPDF40 (+alpha_s
//                  0.116..0.120), CT18NNLO (+alpha_s), CT18Z/A/X, MSHT20 (x3),
//                  PDF4LHC21, HERAPDF20, ABMP16 -- swapped on both beams, R kept
//   idx  44        CT18ANLO member 0 (== the nominal LHAPDF set)
//   idx  45-102    CT18ANLO members 1-58 (29 eigenvector pairs, 90 % CL)
//   idx 103-106    CT18ANLO alpha_s 0.116 / 0.117 / 0.119 / 0.120
//   idx 107-109    replica-group centrals: NNPDF31_nnlo_as_0118_mc,
//                  NNPDF40_nnlo_pdfas, NNPDF40_nnlo_pch_as_01180
//   idx 110        id 9001  nPDFerrSet=1: central EPPS21 R (== nominal)
//   idx 111-158    ids 9002-9049  nPDFerrSet 2-49: nuclear eigen-directions (24 pairs)
//   idx 159-216    ids 9050-9107  nPDFerrSet 50-107: CT18A-baseline directions
//                  of R (29 pairs; R only, see above)
//   9 + 35 + 59 + 4 + 3 + 107 = 217.
// EPPS21 / CT18 error sets are delivered at 90 % CL -> divide by 1.645 for 68 %.
//
// CROSS-CHECKED 2026-09-02 against LHAPDF (gitlab hepcedar/lhapdf, main):
// PDFSet::uncertainty(values, cl = CL1SIGMA) for ErrorType "hessian" loops
// ieigen = 1..nmem/2 over members (2i-1, 2i) and computes exactly the Hessian()
// below -- errplus/errminus with the max(.,0) per pair, errsymm = 0.5*sqrt(sum
// (v_{2i-1}-v_{2i})^2) -- then rescales by sqrt(chi2quantile(68.27%)/
// chi2quantile(90%)) = 1/1.645 because the default cl is 1 sigma (a negative cl
// keeps the set's own CL). So a colleague's `pset.uncertainty(val)` on the 107
// EPPS21 values == this macro's asymmetric Hessian / 1.645, bin by bin.
// EPPS21 paper (arXiv:2112.12462, Sect. 4.2-4.3): 24 free parameters, tolerance
// Delta chi2 ~ 33 = 90 % CL; baseline sets = the nuclear fit REPEATED with each
// CT18ANLO error set S_i^+- as baseline (so their sign labels follow CT18A's);
// the full uncertainty = nuclear (+) baseline in quadrature, Eq. (39). Eq. (40)
// DEFINES the proton-error cross sections as
//     sigma(S_i^+-), i = 25..53  =  f^p_{i-24,+-} (x) sigma_hat (x) f^A_{i,+-}
// i.e. the PROTON-BEAM PDF is the CT18A error set i-24 AND the nuclear PDF is the
// refit made with that baseline (f^A = R_i x f^p_{i-24}). Our idx 159-216 supply
// only R_i (f^p fixed at central on both beams); the CT18A member on both beams is
// idx 44+m. Their per-event product is exactly the Eq. (40) configuration -- so
// the EPPS21 block alone is NOT "the CT18A variation", block B is required. EPPS21.f numbers
// even psets = S-, odd = S+, while the CTEQ convention numbers odd LHAPDF members
// = "+" direction: if both hold, the coherent product pairing is CROSSED,
// idx 44+(2i-1) <-> idx 158+2i and idx 44+2i <-> idx 158+(2i-1), giving a
// baseline term of 3.0 % (W+) instead of 3.7 % for the naive pairing. The CTEQ
// sign convention still has to be verified (CT18 paper / CT18ANLO .info).
//
// For each requested MC sample this macro loops ALL events (no selection, like
// count_ngen.C) and reports, per weight index i:
//   S_i/S_0 - 1   = the relative change of the INCLUSIVE cross section
//                   (S_i = Sum_events ttbar_w[i]),
//   rms(w_i/w_0)  = the standard deviation over events of the per-event ratio
//                   (how much the member re-shapes the sample, as opposed to
//                   re-scaling it). Plain unweighted std. dev. about the mean;
//                   tail-dominated for members whose ratio blows up on the
//                   near-zero-weight (NLO sign-flip) events -- see WriteTable.
// plus per-block summaries (envelope / Hessian sym+asym), and, when both W
// charges are present, the same for the W+/W- ratio and the charge asymmetry.
//
// NB for the analysis: sigma_meas = r x sigma_gen-fid and the count-based
// observables depend on the PDF only through the MC acceptance x efficiency
// and template shapes -- NOT through the inclusive sigma shift quoted here
// (kSigma in mc_norm.h is fixed; the total-sigma variation cancels in r x
// sigma_gen). The per-bin A x eps variation is the relevant systematic and is
// the natural next use of these weights (fill the skim histograms with
// w * ttbar_w[i]/ttbar_w[0]).
//
// Run from skim/ (or ./run_lhe_weights.sh, which tees the log):
//     root -l -b -q 'lhe_weights.C+'                      // Wp, Wm, DY (mu files)
//     root -l -b -q 'lhe_weights.C+("Wp,Wm,DY", "ele")'   // electron files
//     root -l -b -q 'lhe_weights.C+("Wp", "mu", 200000)'   // quick look
//     root -l -b -q 'lhe_weights.C+("draw")'              // redraw the figure from
//                                                          // rootfile/lhe_weights.root (no loop)
// Input paths come from ResolveMCSample (skim_common.h) -- never hardcoded.
//
// Outputs (dirs created if absent):
//   output/lhe_weights.txt                  -- per-index table (with the set/member/id of every
//                                              weight from the confirmed map) + block summaries
//   output/lhe_weights_structure.png        -- the two-pad figure, regions A/B/C only (slide-ready)
//   output/lhe_weights_structure_named.png  -- same + the confirmed set content per region and
//                                              dashed sub-block boundaries
//   rootfile/lhe_weights.root         -- h_wrel_<label> (S_i/S_0), h_wrms_<label>
//                                        (RMS of the per-event ratio), h_sumw_<label>

#include "skim_common.h"
#include "mc_norm.h" // pONorm::kSigma_Wp / kSigma_Wm -- the nominal W+/W- ratio

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TCollection.h"
#include "TKey.h"
#include "TPad.h"
#include "TBox.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace pOSkim;

namespace
{

const int kNW = 217; // expected vector length (checked at run time)

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
// The reference copy of the generator scripts every `where` below points into.
const char *kRefDir = "skim/reference/genproductions_a1a26254/";
const char *kRefZip = "/Users/zhenghuang/Downloads/genproductions_scripts-a1a2625485abc436ed10809327323d5d5daefc5e.zip";

const std::vector<Block> kBlocks = {
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
const Block &BlockOf(int i)
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
const std::vector<SetRun> kSets = {
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
const SetRun &SetOf(int i)
{
  for (const SetRun &s : kSets)
    if (i >= s.lo && i < s.lo + s.n) return s;
  return kSets.back();
}
const char *kScaleLabel[9] = {"muR=1 muF=1", "muR=1 muF=2", "muR=1 muF=0.5", "muR=2 muF=1", "muR=2 muF=2",
                              "muR=2 muF=0.5", "muR=0.5 muF=1", "muR=0.5 muF=2", "muR=0.5 muF=0.5"};

// "CT18ANLO m12  (id 4412, lhapdf 14612)" / "muR=2 muF=0.5  (id 1006)" /
// "EPPS21 R-factor (O16) nPDFerrSet=17  (id 9017)"
std::string SetLabel(int i)
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

// ---- per-sample accumulator -------------------------------------------------
struct WSample
{
  std::string           label;   // "Wp_mu", "DY_mu", ...
  std::string           fname;
  Long64_t              nEvents = 0;
  Long64_t              nUsed   = 0;
  Long64_t              nBadLen = 0;
  size_t                len     = 0;
  double                sumW0   = 0; // Sum HiTree::weight
  std::vector<double>   sumW;        // Sum ttbar_w[i]
  std::vector<double>   sumR, sumR2; // Sum (ttbar_w[i]/ttbar_w[0]), squared
  double                maxDev0 = 0; // max |ttbar_w[0]/weight - 1|
  bool                  ok      = false;

  double rel(int i) const { return sumW[i] / sumW[0]; }              // S_i/S_0
  double rms(int i) const                                              // RMS of per-event ratio
  {
    const double m = sumR[i] / nUsed;
    return std::sqrt(std::max(0.0, sumR2[i] / nUsed - m * m));
  }
};

SampleType SampleFromToken(const std::string &tok)
{
  if (tok == "Wp")    return kWp;
  if (tok == "Wm")    return kWm;
  if (tok == "DY")    return kDY;
  if (tok == "Wptau") return kWptau;
  if (tok == "Wmtau") return kWmtau;
  if (tok == "DYtau") return kDYtau;
  std::cerr << "[ERR] lhe_weights: unknown sample token '" << tok
            << "' (use Wp, Wm, DY, Wptau, Wmtau, DYtau)\n";
  return kData;
}

WSample LoopOne(const std::string &tok, const char *flavour, Long64_t nmax)
{
  WSample s;
  const SampleType st = SampleFromToken(tok);
  if (st == kData) return s;
  const SampleFileInfo info = ResolveMCSample(st, flavour);
  s.fname = info.fname;
  const bool tauFile = (st == kWptau || st == kWmtau || st == kDYtau);
  s.label = tok + (tauFile ? "" : std::string("_") + flavour);

  std::cout << "[INPUT] " << s.fname << std::endl;
  TFile *f = TFile::Open(s.fname.c_str());
  if (!f || f->IsZombie())
  {
    std::cerr << "[ERR] lhe_weights: cannot open " << s.fname << "\n";
    if (f) { f->Close(); delete f; }
    return s;
  }
  TTree *t = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  if (!t || !HasBranch(t, "weight") || !HasBranch(t, "ttbar_w"))
  {
    std::cerr << "[ERR] lhe_weights: HiTree with 'weight' + 'ttbar_w' not found in " << s.fname << "\n";
    f->Close(); delete f;
    return s;
  }

  Float_t             weight = 0.f;
  std::vector<float> *ww     = nullptr;
  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("weight", 1);
  t->SetBranchStatus("ttbar_w", 1);
  t->SetBranchAddress("weight", &weight);
  t->SetBranchAddress("ttbar_w", &ww);

  Long64_t n = t->GetEntries();
  if (nmax > 0 && nmax < n) n = nmax;
  s.nEvents = n;

  for (Long64_t ie = 0; ie < n; ++ie)
  {
    t->GetEntry(ie);
    const size_t L = ww ? ww->size() : 0;
    if (s.len == 0 && L > 0)
    {
      s.len = L;
      s.sumW.assign(L, 0.0); s.sumR.assign(L, 0.0); s.sumR2.assign(L, 0.0);
      if ((int)L != kNW)
        std::cout << "[WARN] " << s.label << ": ttbar_w has " << L << " entries, expected " << kNW
                  << " -- block table may not apply\n";
    }
    if (L != s.len || L == 0) { ++s.nBadLen; continue; }
    const double w0 = (*ww)[0];
    if (w0 == 0.0) { ++s.nBadLen; continue; }
    ++s.nUsed;
    s.sumW0 += weight;
    const double dev0 = std::fabs(w0 / weight - 1.0);
    if (dev0 > s.maxDev0) s.maxDev0 = dev0;
    for (size_t i = 0; i < L; ++i)
    {
      const double wi = (*ww)[i];
      const double r  = wi / w0;
      s.sumW[i]  += wi;
      s.sumR[i]  += r;
      s.sumR2[i] += r * r;
    }
    if (ie % 200000 == 0 && ie > 0)
      std::cout << "       " << s.label << ": " << ie << " / " << n << "\n";
  }
  s.ok = (s.nUsed > 0);
  std::cout << "       " << s.label << ": events " << n << ", used " << s.nUsed
            << ", length " << s.len << ", length/zero rejects " << s.nBadLen
            << ", Sum(weight) = " << std::setprecision(8) << s.sumW0
            << ", Sum(ttbar_w[0]) = " << s.sumW[0]
            << ", max|ttbar_w[0]/weight - 1| = " << std::scientific << std::setprecision(2)
            << s.maxDev0 << std::defaultfloat << "\n";
  f->Close(); delete f;
  return s;
}

// ---- summaries ---------------------------------------------------------------
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

template <class F>
void PrintBlockSummary(std::ostream &os, const char *what, F dev)
{
  os << "  " << what << "\n";
  for (const Block &b : kBlocks)
  {
    const int n = b.hi - b.lo + 1;
    double mean = 0, s2 = 0, mn = 1e9, mx = -1e9;
    for (int i = b.lo; i <= b.hi; ++i)
    {
      const double d = dev(i);
      mean += d; s2 += d * d; mn = std::min(mn, d); mx = std::max(mx, d);
    }
    mean /= n;
    const double sd = std::sqrt(std::max(0.0, s2 / n - mean * mean));
    os << "    idx " << std::setw(3) << b.lo << "-" << std::setw(3) << b.hi
       << "  n=" << std::setw(3) << n << "  " << std::left << std::setw(62) << b.name << std::right;
    if (b.kind == kCentral)
      os << "  dev = " << std::showpos << std::fixed << std::setprecision(3) << 100 * mean << std::noshowpos << " %\n";
    else
    {
      os << std::fixed << std::setprecision(2)
         << "  envelope " << std::showpos << 100 * mx << " / " << 100 * mn << std::noshowpos << " %"
         << "  mean " << std::showpos << 100 * mean << std::noshowpos << " %  std " << 100 * sd << " %";
      if (b.kind == kHessian)
      {
        const HessRes h = Hessian(b, dev);
        os << "  | Hessian sym " << 100 * h.sym << " %, asym +" << 100 * h.up << " / -" << 100 * h.dn
           << " %  (68% CL: sym " << 100 * h.sym / 1.645 << " %)";
      }
      os << "\n";
    }
  }
}

void WriteTable(const std::vector<WSample> &S, const std::string &outFile)
{
  std::ofstream os(outFile.c_str(), std::ios::out | std::ios::trunc);
  if (!os.is_open()) { std::cerr << "[ERR] cannot open " << outFile << "\n"; return; }

  os << "# hiEvtAnalyzer/HiTree::ttbar_w = the LHE reweighting block, rescaled so ttbar_w[0] == HiTree::weight\n"
     << "# (HiEvtAnalyzer.cc: ttbar_w[i] = weight/originalXWGTUP * LHEEventProduct::weights()[i]).\n"
     << "#\n"
     << "# REFERENCE SCRIPTS -- every file:line below points into the repo copy\n"
     << "#   " << kRefDir << "   (README.md there)\n"
     << "# copied 2026-09-02 from the CMS genproductions archive used for this production:\n"
     << "#   " << kRefZip << "\n"
     << "#   bin/Powheg/make_rwl.py                         writes pwg-rwl.dat = the <initrwgt> header POWHEG reads at generation -> idx 0-109\n"
     << "#                                                  (\"EPPS21\" in Period branch, L481: scale grid L43-54 + hessian group L487-530 + replica group L532-537;\n"
     << "#                                                  groups written in sorted-key order L548, one weight per member L559)\n"
     << "#   bin/Powheg/runcmsgrid_powheg.sh                appends the EPPS21 nPDF weights AFTER generation -> idx 110-216 (L252-266: one pwhg_main rerun\n"
     << "#                                                  with rwl_add 1 per nPDFerrSet = 1..107, weight id 9000+iset, weightgroup 'EPPS21_variation')\n"
     << "#   bin/Powheg/run_pwg_condor.py                   L324-326: for an EPPS21_* ion, defaultPDF = 14600 (CT18ANLO) and period = \"Run3_\"+ion;\n"
     << "#                                                  L648-667: the card must have lhans1 == lhans2 == 14600\n"
     << "#   bin/Powheg/Templates/runGetSource_template.sh  L55 calls make_rwl.py; L110-116 apply the EPPS21 patches + fetch the EPPS21NLOR_16 R grid\n"
     << "#   bin/Powheg/patches/EPPS21/*.patch, EPPS21.f    the nuclear PDF: LHAPDF proton PDF x EPPS21 R(nPDFerrSet) on the ia>=16 beam (lhapdf6if_nPDF.patch\n"
     << "#                                                  L38-40, L80-95), isospin-averaged (L63-67); pdfcalls_nPDF.patch passes ia1/ia2 per beam\n"
     << "#   MetaData/npdflist_O_5f_run3.dat                the Run3_O branch's LHAPDF-grid nPDF list -- NOT used by this production, kept for contrast\n"
     << "# Mechanism: f_O = R_EPPS21(nPDFerrSet) x f_p^LHAPDF. `lhapdf=X` weights (idx 9-109) swap f_p on BOTH beams with R at set 1;\n"
     << "# `nPDFerrSet` weights (idx 110-216) vary R with f_p fixed at CT18ANLO central. Hence idx 44 and idx 110 are both exactly the nominal.\n"
     << "# Data checks of the map: exact ==nominal at 0/44/110, monotonic alpha_s scans (11-18, 23-28, 30-34, 103-106), 38==39, 107==9,\n"
     << "# scale-grid algebra (idx 1,2,3,6 single-scale), nuclear/baseline character change at 158/159.\n"
     << "#\n"
     << "# PROVENANCE BY BLOCK  (idx range | produced by | why it is in the list)\n";
  for (const Block &b : kBlocks)
    os << "#   idx " << std::setw(3) << b.lo << "-" << std::setw(3) << b.hi << "  [" << b.tag << "]  " << b.name << "\n"
       << "#              where: " << b.where << "\n"
       << "#              why:   " << b.why << "\n";
  os << "#\n# Samples:\n";
  for (const WSample &s : S)
    os << "#   " << std::left << std::setw(8) << s.label << std::right
       << "  events " << s.nEvents << "  used " << s.nUsed << "  vector length " << s.len
       << "  Sum(weight) " << std::fixed << std::setprecision(1) << s.sumW0
       << "  max|w[0]/weight-1| " << std::scientific << std::setprecision(1) << s.maxDev0
       << std::defaultfloat << "  " << s.fname << "\n";
  os << "#\n# Per-index table. rel = 100*(S_i/S_0 - 1) [%] = relative change of the inclusive cross section;\n"
     << "# rms = standard deviation over events of the per-event ratio ttbar_w[i]/ttbar_w[0] (ROOT's 'RMS'\n"
     << "#       convention, i.e. about the mean). NB it is a plain unweighted std. dev., so it is dominated by\n"
     << "#       the few events whose nominal weight is close to zero (NLO sign flips): for the region-A members\n"
     << "#       the robust MAD-based spread is ~6x smaller than the std. dev. (0.025 vs 0.157 at idx 5), while\n"
     << "#       for region B/C members the two agree.\n"
     << "#  idx";
  for (const WSample &s : S) os << std::setw(11) << ("rel_" + s.label) << std::setw(10) << ("rms_" + s.label);
  os << "   " << std::left << std::setw(70) << "LHE weight (set, member, id, lhapdf)"
     << std::setw(62) << "written by (file:line in the reference copy)" << "role\n" << std::right;
  for (int i = 0; i < kNW; ++i)
  {
    os << std::setw(5) << i;
    for (const WSample &s : S)
    {
      if (!s.ok || i >= (int)s.len) { os << std::setw(11) << "-" << std::setw(10) << "-"; continue; }
      os << std::fixed << std::setprecision(3) << std::showpos << std::setw(11) << 100 * (s.rel(i) - 1.0)
         << std::noshowpos << std::setprecision(4) << std::setw(10) << s.rms(i);
    }
    os << "   " << std::left << std::setw(70) << SetLabel(i) << std::setw(62) << SetOf(i).src
       << BlockOf(i).tag << std::right;
    for (const Block &b : kBlocks)
      if (i == b.lo) os << "   <-- block start: " << b.name;
    os << "\n";
  }

  os << "\n# ---------------- BLOCK SUMMARY (per sample, inclusive cross section) ----------------\n";
  for (const WSample &s : S)
  {
    if (!s.ok) continue;
    PrintBlockSummary(os, (s.label + ":  d_i = S_i/S_0 - 1").c_str(),
                      [&](int i) { return s.rel(i) - 1.0; });
  }

  // W+/W- ratio and charge asymmetry, when both charges are present (same flavour).
  const WSample *wp = nullptr, *wm = nullptr;
  for (const WSample &s : S)
  {
    if (s.ok && s.label.rfind("Wp", 0) == 0) wp = &s;
    if (s.ok && s.label.rfind("Wm", 0) == 0) wm = &s;
  }
  if (wp && wm)
  {
    os << "\n# ---------------- W+/W- RATIO and CHARGE ASYMMETRY (weights only; the absolute sigma_W+/sigma_W- from mc_norm.h) ----------------\n";
    const double R0 = pONorm::kSigma_Wp / pONorm::kSigma_Wm; // nominal ratio; member i scales it by rel_p/rel_m
    PrintBlockSummary(os, "W+/W- ratio:  d_i = R_i/R_0 - 1",
                      [&](int i) { return wp->rel(i) / wm->rel(i) - 1.0; });
    PrintBlockSummary(os, "charge asymmetry A = (s+ - s-)/(s+ + s-):  d_i = A_i - A_0  (ABSOLUTE shift x100)",
                      [&](int i)
                      {
                        const double Ri = R0 * wp->rel(i) / wm->rel(i);
                        return (Ri - 1) / (Ri + 1) - (R0 - 1) / (R0 + 1);
                      });
    os << "    (A_0 = " << std::fixed << std::setprecision(4) << (R0 - 1) / (R0 + 1)
       << " from kSigma_Wp/kSigma_Wm = " << R0 << ")\n";
  }
  os.close();
  std::cout << "[INFO] wrote " << outFile << "\n";
}

void WriteRoot(const std::vector<WSample> &S)
{
  gSystem->mkdir("rootfile", kTRUE);
  TFile *fout = new TFile("./rootfile/lhe_weights.root", "RECREATE");
  for (const WSample &s : S)
  {
    if (!s.ok) continue;
    const int L = (int)s.len;
    TH1D *hrel  = new TH1D(("h_wrel_" + s.label).c_str(),
                           (s.label + ";ttbar_w index;S_{i}/S_{0}").c_str(), L, -0.5, L - 0.5);
    TH1D *hrms  = new TH1D(("h_wrms_" + s.label).c_str(),
                           (s.label + ";ttbar_w index;std. dev. over events of w_{i}/w_{0}").c_str(), L, -0.5, L - 0.5);
    TH1D *hsumw = new TH1D(("h_sumw_" + s.label).c_str(),
                           (s.label + ";ttbar_w index;#Sigma_{events} ttbar_w[i]").c_str(), L, -0.5, L - 0.5);
    for (int i = 0; i < L; ++i)
    {
      hrel->SetBinContent(i + 1, s.rel(i));
      hrms->SetBinContent(i + 1, s.rms(i));
      hsumw->SetBinContent(i + 1, s.sumW[i]);
    }
    hrel->Write(); hrms->Write(); hsumw->Write();
  }
  fout->Close(); delete fout;
  std::cout << "[INFO] wrote rootfile/lhe_weights.root\n";
}

// ---- drawing -----------------------------------------------------------------
// The FIGURE is deliberately interpretation-free: it shows the numbers and
// splits the index axis into REGIONS that start at every member identical to
// the nominal weight (S_i/S_0 == 1 with zero event-by-event spread, i.e. a PDF
// central that reproduces the generation PDF), labelled A, B, C, ... The block
// names in kBlocks (scale / EPPS21 / ...) are a working hypothesis pending the
// LHE header from the MC producer and appear ONLY in the text outputs.

// Light-weight input so the figure can be redrawn from rootfile/lhe_weights.root
// without re-looping the ntuples (mode "draw").
struct DrawInput
{
  std::string         label;
  std::vector<double> rel; // S_i/S_0
  std::vector<double> rms; // RMS over events of w_i/w_0
};

DrawInput ToDrawInput(const WSample &s)
{
  DrawInput d;
  d.label = s.label;
  for (size_t i = 0; i < s.len; ++i) { d.rel.push_back(s.rel(i)); d.rms.push_back(s.rms(i)); }
  return d;
}

std::vector<DrawInput> ReadDrawInputs(const char *rootfile)
{
  std::vector<DrawInput> out;
  TFile *f = TFile::Open(rootfile);
  if (!f || f->IsZombie())
  {
    std::cerr << "[ERR] lhe_weights: cannot open " << rootfile << " (run the event loop first)\n";
    if (f) { f->Close(); delete f; }
    return out;
  }
  TIter next(f->GetListOfKeys());
  while (TKey *key = (TKey *)next())
  {
    const std::string n = key->GetName();
    if (n.rfind("h_wrel_", 0) != 0) continue;
    const std::string label = n.substr(std::string("h_wrel_").size());
    TH1D *hrel = (TH1D *)f->Get(n.c_str());
    TH1D *hrms = (TH1D *)f->Get(("h_wrms_" + label).c_str());
    if (!hrel || !hrms) continue;
    DrawInput d;
    d.label = label;
    for (int i = 1; i <= hrel->GetNbinsX(); ++i)
    {
      d.rel.push_back(hrel->GetBinContent(i));
      d.rms.push_back(hrms->GetBinContent(i));
    }
    out.push_back(d);
    std::cout << "[INFO] draw: read " << label << " (" << d.rel.size() << " weights) from " << rootfile << "\n";
  }
  f->Close(); delete f;
  return out;
}

// Region starts = every index at which ALL samples have S_i/S_0 == 1 and no
// event-by-event spread (tolerance 1e-3 covers the float rounding of the LHE
// text). Index 0 always opens region A.
std::vector<int> RegionStarts(const std::vector<DrawInput> &D, int L)
{
  std::vector<int> starts;
  for (int i = 0; i < L; ++i)
  {
    bool nominalLike = !D.empty();
    for (const DrawInput &d : D)
      if (i >= (int)d.rel.size() || std::fabs(d.rel[i] - 1.0) > 1e-3 || d.rms[i] > 1e-3)
      {
        nominalLike = false;
        break;
      }
    if (nominalLike) starts.push_back(i);
  }
  if (starts.empty() || starts.front() != 0) starts.insert(starts.begin(), 0);
  return starts;
}

// named = false: regions A/B/C only (interpretation-free).
// named = true : additionally the confirmed set content per region and dashed
//                separators at the sub-block boundaries of kBlocks.
void Draw(const std::vector<DrawInput> &D, const std::string &outPng, bool named = false)
{
  int L = 0;
  for (const DrawInput &d : D) L = std::max(L, (int)d.rel.size());
  if (L == 0) { std::cerr << "[ERR] lhe_weights: nothing to draw\n"; return; }
  const std::vector<int> starts = RegionStarts(D, L);
  auto regionEnd = [&](size_t r) { return (r + 1 < starts.size() ? starts[r + 1] - 1 : L - 1); };
  const char *sfx = named ? "_named" : "";

  std::cout << "[INFO] regions (start = member identical to the nominal):";
  for (size_t r = 0; r < starts.size(); ++r)
    std::cout << "  " << (char)('A' + (int)r) << " = idx " << starts[r] << "-" << regionEnd(r);
  std::cout << "\n";

  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas(Form("c_lhe%s", sfx), "ttbar_w structure", 1500, 800);
  TPad *p1 = new TPad(Form("p1%s", sfx), "", 0, 0.36, 1, 1);
  TPad *p2 = new TPad(Form("p2%s", sfx), "", 0, 0, 1, 0.36);
  p1->SetBottomMargin(0.02); p1->SetTopMargin(0.09); p1->SetLeftMargin(0.07); p1->SetRightMargin(0.02);
  p2->SetTopMargin(0.03); p2->SetBottomMargin(0.28); p2->SetLeftMargin(0.07); p2->SetRightMargin(0.02);
  p1->Draw(); p2->Draw();

  const int    cols[3]  = {kRed + 1, kBlue + 1, kGreen + 2};
  const int    mks[3]   = {20, 21, 22};
  const int    shade[2] = {kGray, kAzure - 9};
  const double yLo = -8.5, yHi = (named ? 12.5 : 10.8); // named: room for a third label line
  const double yLab = yHi - 1.3;
  // sub-block boundaries (kBlocks starts that are not region starts), named mode only
  std::vector<int> subStarts;
  if (named)
    for (const Block &b : kBlocks)
      if (std::find(starts.begin(), starts.end(), b.lo) == starts.end()) subStarts.push_back(b.lo);
  auto drawSubSeps = [&](double lo, double hi)
  {
    for (int s0 : subStarts)
    {
      TLine *lb = new TLine(s0 - 0.5, lo, s0 - 0.5, hi); lb->SetLineColor(kGray + 2); lb->SetLineStyle(2); lb->Draw("same");
    }
  };

  auto shadeRegions = [&](double lo, double hi)
  {
    for (size_t r = 0; r < starts.size(); ++r)
    {
      TBox *bx = new TBox(starts[r] - 0.5, lo, regionEnd(r) + 0.5, hi);
      bx->SetFillColorAlpha(shade[r % 2], 0.25); bx->SetLineWidth(0);
      bx->Draw("same");
    }
  };

  // ---- top: relative change of the inclusive cross section
  p1->cd();
  TH1D *fr = new TH1D(Form("fr_lhe%s", sfx), ";;100 #times (S_{i}/S_{0} #minus 1)  [%]", L, -0.5, L - 0.5);
  fr->SetMinimum(yLo); fr->SetMaximum(yHi);
  fr->GetXaxis()->SetLabelSize(0); fr->GetYaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetLabelSize(0.045);
  fr->GetYaxis()->SetTitleOffset(0.6);
  fr->Draw("axis");
  shadeRegions(yLo, yHi);
  drawSubSeps(yLo, yHi);
  TLine *l0 = new TLine(-0.5, 0, L - 0.5, 0); l0->SetLineStyle(2); l0->Draw("same");
  for (int s0 : starts)
  {
    TLine *lc = new TLine(s0, yLo, s0, yHi); lc->SetLineColor(kBlack); lc->SetLineStyle(3); lc->Draw("same");
  }
  TLatex tx; tx.SetTextAlign(21);
  std::string startsTxt;
  // confirmed content per region (make_rwl.py map, see header) -- named mode only
  const char *content[3] = {"nominal | 8 scale | 35 PDF centrals",
                            "CT18ANLO m0 | 58 eig. | 4 #alpha_{s} | 3 replica centrals",
                            "EPPS21 R-factor sets: 1 central | 2-49 nuclear | 50-107 CT18A-baseline (R only)"};
  for (size_t r = 0; r < starts.size(); ++r)
  {
    const int a = starts[r], b = regionEnd(r);
    const double xc = 0.5 * (a + b);
    tx.SetTextSize(0.046);
    tx.DrawLatex(xc, yLab, Form("#bf{Region %c}", 'A' + (int)r));
    tx.SetTextSize(0.036);
    tx.DrawLatex(xc, yLab - 1.1, Form("idx %d #minus %d   (%d weights)", a, b, b - a + 1));
    if (named && r < 3)
    {
      tx.SetTextSize(0.030);
      tx.DrawLatex(xc, yLab - 2.3, content[r]);
    }
    startsTxt += (startsTxt.empty() ? "" : ", ") + std::to_string(a);
  }
  tx.SetTextAlign(11); tx.SetTextSize(0.034);
  tx.DrawLatex(112, -7.6, ("dotted: member identical to the nominal weight = region start (idx " + startsTxt + ")").c_str());
  if (named) tx.DrawLatex(112, -6.7, "dashed: sub-block boundaries of the <initrwgt> map (make_rwl.py)");

  // legend sits in the lower part of region B, where the points stay within +-1.5 %
  TLegend *leg = new TLegend(0.27, 0.07, 0.49, 0.27);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.04);
  for (size_t j = 0; j < D.size() && j < 3; ++j)
  {
    const DrawInput &d = D[j];
    TH1D *h = new TH1D(("g_rel_" + d.label + sfx).c_str(), "", L, -0.5, L - 0.5);
    for (int i = 0; i < (int)d.rel.size() && i < L; ++i) h->SetBinContent(i + 1, 100 * (d.rel[i] - 1.0));
    h->SetMarkerStyle(mks[j]); h->SetMarkerColor(cols[j]); h->SetLineColor(cols[j]); h->SetMarkerSize(0.9);
    h->Draw("P same");
    leg->AddEntry(h, (d.label + "  (#Sigma w over all events)").c_str(), "p");
  }
  leg->Draw();
  tx.SetTextAlign(11); tx.SetTextSize(named ? 0.044 : 0.05);
  tx.DrawLatex(-0.5, yHi + 0.3 * (yHi - yLo) / 19.3,
               Form("#bf{HiTree::ttbar_w}  =  LHE reweighting block, %d weights / event  (POWHEG pO 9.62 TeV MC)%s",
                    L, named ? "  |  map: genproductions scripts" : ""));

  // ---- bottom: event-by-event RMS of the ratio (log)
  p2->cd(); p2->SetLogy();
  TH1D *fr2 = new TH1D(Form("fr2_lhe%s", sfx), ";ttbar_w index;std. dev._{events}(w_{i}/w_{0})", L, -0.5, L - 0.5);
  fr2->SetMinimum(1e-3); fr2->SetMaximum(2.0);
  fr2->GetXaxis()->SetTitleSize(0.10); fr2->GetXaxis()->SetLabelSize(0.085); fr2->GetXaxis()->SetTitleOffset(1.1);
  fr2->GetYaxis()->SetTitleSize(0.075); fr2->GetYaxis()->SetLabelSize(0.075); fr2->GetYaxis()->SetTitleOffset(0.42);
  fr2->GetXaxis()->SetNdivisions(522);
  fr2->Draw("axis");
  shadeRegions(1e-3, 2.0);
  drawSubSeps(1e-3, 2.0);
  for (int s0 : starts)
  {
    TLine *lc = new TLine(s0, 1e-3, s0, 2.0); lc->SetLineColor(kBlack); lc->SetLineStyle(3); lc->Draw("same");
  }
  for (size_t j = 0; j < D.size() && j < 3; ++j)
  {
    const DrawInput &d = D[j];
    TH1D *h = new TH1D(("g_rms_" + d.label + sfx).c_str(), "", L, -0.5, L - 0.5);
    for (int i = 0; i < (int)d.rms.size() && i < L; ++i) h->SetBinContent(i + 1, std::max(d.rms[i], 1.1e-3));
    h->SetMarkerStyle(mks[j]); h->SetMarkerColor(cols[j]); h->SetLineColor(cols[j]); h->SetMarkerSize(0.8);
    h->Draw("P same");
  }

  c->SaveAs(outPng.c_str());
  std::cout << "[INFO] wrote " << outPng << "\n";
}

} // namespace

// ============================================================================
void lhe_weights(const char *samples = "Wp,Wm,DY", const char *flavour = "mu", Long64_t nmax = -1)
{
  gSystem->mkdir("output", kTRUE);
  gSystem->mkdir("rootfile", kTRUE);

  // Redraw-only mode: rebuild the figure from the stored rootfile (seconds).
  if (std::string(samples) == "draw")
  {
    const std::vector<DrawInput> D = ReadDrawInputs("./rootfile/lhe_weights.root");
    if (D.empty()) { std::cerr << "[ERR] lhe_weights: no h_wrel_* histograms to draw\n"; return; }
    Draw(D, "./output/lhe_weights_structure.png");             // regions A/B/C only
    Draw(D, "./output/lhe_weights_structure_named.png", true); // + the confirmed set content
    return;
  }

  std::vector<WSample> S;
  std::stringstream ss(samples);
  std::string tok;
  while (std::getline(ss, tok, ','))
  {
    if (tok.empty()) continue;
    S.push_back(LoopOne(tok, flavour, nmax));
  }
  bool any = false;
  for (const WSample &s : S) any |= s.ok;
  if (!any) { std::cerr << "[ERR] lhe_weights: nothing processed\n"; return; }

  const std::string tag = (S.size() == 3 && std::string(flavour) == "mu" && nmax <= 0)
                              ? "" : std::string("_") + flavour + (nmax > 0 ? "_partial" : "");
  WriteTable(S, "./output/lhe_weights" + tag + ".txt");
  WriteRoot(S);
  std::vector<DrawInput> D;
  for (const WSample &s : S)
    if (s.ok) D.push_back(ToDrawInput(s));
  Draw(D, "./output/lhe_weights_structure" + tag + ".png");
  Draw(D, "./output/lhe_weights_structure" + tag + "_named.png", true);

  // headline echo for the log
  std::cout << "\n=== ttbar_w headline (inclusive sigma, S_i/S_0 - 1) ===\n";
  for (const WSample &s : S)
  {
    if (!s.ok) continue;
    PrintBlockSummary(std::cout, s.label.c_str(), [&](int i) { return s.rel(i) - 1.0; });
  }
}
