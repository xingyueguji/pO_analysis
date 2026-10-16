#include "TCanvas.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TLatex.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

#include "TGraph.h"
#include "TFile.h"
#include "TLine.h"
#include "TMath.h"
#include "plotting_helper.C"           // PlotStyle, ApplyCanvasStyle/HistStyle, DrawHeader, CMS_lumi
#include "../skim/mc_norm.h"           // pONorm::kLumi_invnb (single-source data lumi)
#include "disc_variants.h"             // pODisc::Spec/GraphFile (W-discriminant tags)
#include "fit_variants.h"              // pOFit::Spec (which fit: comb / mu-only / e-only)

// -----------------------------------------------------------------------------
// Charge-series markers (2026-09-15) -- one place, so the W / W+ / W- glyphs
// stay consistent across every plot in this macro.
//
// TWO rules, both from a pixel measurement of ROOT's rasterizer (markers drawn
// on top of a width-2 line, glyph centre vs line centre, sizes 1.1-2.0):
//
//  (1) POINT-SYMMETRIC GLYPHS ONLY -- no triangles. A triangle's visual centre
//      sits ~1/3 of the height above its base, so the eye reads the data point
//      systematically low, and style 22 additionally rasterizes +0.5..+1.0 px
//      off in y (the open triangle 26 is worse: +1.0..+1.5 px). Circle (20/24)
//      and diamond (33/27) are symmetric under a half-turn, so their centre is
//      unambiguous whatever the size.
//
//  (2) SIZES THAT RASTERIZE CENTRED. ROOT draws each glyph as an integer pixel
//      block; when that block comes out an even number of pixels wide it cannot
//      straddle the error bar exactly, and the marker lands up to 1.5 px off.
//      Measured offsets (dx, dy in px): circle 20 = (+0.5, 0.0) at ANY size;
//      square 21 = (-1.0, -1.5) at 1.3 but (0.0, -0.5) at 1.5; diamond 33 =
//      (-0.5, 0.0) at 2.0 (width 11 px, matching the other two).
//      The square at 1.3 is what made the W+ point look shifted in the legend,
//      where a long horizontal line runs through the glyph and makes a 1.5 px
//      offset obvious. Re-measure before changing a size.
// -----------------------------------------------------------------------------
static const int    kMkW  = 20, kMkWp = 21, kMkWm = 33;  // filled: circle / square / diamond
static const double kMkWSz = 1.3, kMkWpSz = 1.5, kMkWmSz = 2.0;
static const int    kMkWOpen = 24, kMkWmOpen = 27;       // open circle / open diamond

// -----------------------------------------------------------------------------
// pp REFERENCE for the inclusive plot of xsec_fiducial_flav (2026-10-05, user
// request: a COARSE comparison -- "we don't need to be precise here").
//
// CMS SMP-20-004, "Measurement of the inclusive cross sections for W and Z
// boson production in pp collisions at sqrt(s) = 5.02 and 13 TeV",
// JHEP 04 (2025) 162, arXiv:2408.03744, HEPData ins2816048 -- ONE paper, ONE
// fiducial definition at both energies, which bracket our 9.62 TeV (there was
// NO pp run at 9.62 TeV in 2025: the LPC special-run page lists pO 9.62, OO
// 5.36 and NeNe 5.36 TeV only). Combined e+mu FIDUCIAL cross sections, PER
// LEPTON FLAVOUR, in nb, from the paper's Tables 1 (5.02 TeV) and 2 (13 TeV):
//   W: Born lepton pT > 25, |eta| < 2.4, mT > 40        (since 2026-10-13 our
//      leppt_mt40 gen fiducial has the same mT cut; ours is bare, not Born --
//      ~2% (mu) / ~4.5% (e) lower -- NOT corrected)
//   Z: both Born leptons pT > 25, |eta| < 2.4, 60 < m < 120 (ours: 15/10 GeV legs)
// Drawn x kA_O (= 16): sigma_pO = A x sigma_pp for NO nuclear modification AND
// NO isospin effect -- the latter is why W+ and W- separately are NOT expected
// to match even at the right energy (oxygen is 8p + 8n: pp W+/W- ~ 1.4, the
// isospin-averaged POWHEG 1.17); the W+ + W- sum and the Z cancel it to a few %.
// The energy difference is not corrected either (a log-log interpolation of the
// two points to 9.62 TeV gives W 7.02 nb, Z 0.575 nb per flavour).
// Errors drawn = stat (+) syst (+) lumi in quadrature.
// -----------------------------------------------------------------------------
struct PPRef
{
    const char *label;     // legend text
    double v[4];           // W+, W-, W, Z  [nb, per flavour, pp fiducial]
    double stat[4], syst[4], lumi[4];
    int marker; double markerSize; int color;
};
static const PPRef kPPRef[2] = {
    {"CMS pp 5.02 TeV",
     {2.475, 1.525, 4.000, 0.3198}, {0.002, 0.002, 0.003, 0.0009},
     {0.008, 0.005, 0.011, 0.0012}, {0.047, 0.029, 0.076, 0.0062},
     24, 1.3, kOrange + 7},                              // open circle
    {"CMS pp 13 TeV",
     {5.170, 3.932, 9.110, 0.754}, {0.000, 0.004, 0.010, 0.002},
     {0.020, 0.014, 0.030, 0.003}, {0.120, 0.089, 0.210, 0.017},
     25, 1.5, kMagenta + 2},                             // open square (1.5 = centred, see kMk* note)
};
static const char *kPPRefCite = "CMS, JHEP 04 (2025) 162";
static double ppRefErr(const PPRef &p, int k)
{
    return std::sqrt(p.stat[k] * p.stat[k] + p.syst[k] * p.syst[k] + p.lumi[k] * p.lumi[k]);
}

// =============================================================================
// xsec_fiducial.C -- fiducial W cross sections (W+, W-, W inclusive, and
// d(sigma)/d(eta_CM)) from the simultaneous fits:
//
//   sigma_meas,i = r_i x sigma_gen-fid,i        (r_i = the fitted POI of bin i)
//
//   xsec_fiducial_comb(disc)   the GRAND fit (mu + e, r shared)  -> plots/comb/xsec/<disc>/
//   xsec_fiducial_diag(disc)   its per-flavour extraction diagnostic (A x eps)
//   xsec_fiducial_flav(disc)   the MU-ONLY and E-ONLY fits OVERLAID (fork mode
//                              flavfit, 2026-09-22)       -> plots/flavfit/xsec/<disc>/
//   xsec_fiducial(disc)        all three
//
// The LEGACY N_fit/L cross sections of the per-flavour per-bin fits
// (xsec_fiducial(disc, muCsv, eleCsv) and xsec_fiducial_diff) were removed
// with those fits on 2026-09-22; xsec_fiducial_flav is the mu-vs-e view now.
// POLICY (user, 2026-09-22): every observable is built from r x sigma_gen,
// never from raw counts (no dedicated efficiency/acceptance correction
// exists); the count-based r x S = N_fit numbers survive only as console
// cross-checks and in xsec_fiducial_diag, whose point is to SHOW the A x eps
// that r x sigma_gen applies.
// Runner: analysis/run_observables.sh.
// =============================================================================

// -----------------------------------------------------------------------------
// per-bin PREFIT MC signal (mu+e summed template integrals, col 6
// signal_prefit of comb_W_yields.csv) -- the fit's r = 1 expectation; the
// postfit/prefit ratio per bin IS the fitted r_<C>_y<i>.
static bool readBinPrefit(const char *csv, const char *charge, const char *binning,
                          double s[12])
{
    for (int i = 0; i < 12; ++i) s[i] = 0;
    std::ifstream in(csv);
    if (!in) return false;
    std::string line;
    std::getline(in, line); // header
    int found = 0;
    while (std::getline(in, line))
    {
        std::stringstream ss(line);
        std::string f;
        std::vector<std::string> c;
        while (std::getline(ss, f, ',')) c.push_back(f);
        if (c.size() < 9) continue;
        if (c[1] == charge && c[2] == binning)
        {
            int iy = std::atoi(c[3].c_str());
            if (iy >= 0 && iy < 12) { s[iy] = std::atof(c[6].c_str()); ++found; }
        }
    }
    return found > 0;
}

// -----------------------------------------------------------------------------
// comb_summary.csv reader: fit,param,value,error,...  (extract_pO_simfit.C
// writes the covariance-propagated inclusive sums as e.g.
// "simfit_lab,W_sum_yield,<yield>,<err>,,,").
static bool readCombSum(const char *csv, const char *fit, const char *param,
                        double &y, double &e)
{
    std::ifstream in(csv);
    if (!in) { std::cerr << "[ERROR] cannot open CSV: " << csv << "\n"; return false; }
    std::string line;
    std::getline(in, line); // header
    while (std::getline(in, line))
    {
        std::stringstream ss(line);
        std::string f;
        std::vector<std::string> c;
        while (std::getline(ss, f, ',')) c.push_back(f);
        if (c.size() >= 4 && c[0] == fit && c[1] == param)
        {
            y = std::atof(c[2].c_str());
            e = std::atof(c[3].c_str());
            return true;
        }
    }
    std::cerr << "[WARN] row '" << fit << "," << param << "' not found in " << csv << "\n";
    return false;
}

// -----------------------------------------------------------------------------
// per-bin fitted signal strength r +/- rErr (cols r/rErr of comb_W_yields.csv)
// -- the POIs of the grand simultaneous fit; sigma_meas,i = r_i x sigma_gen-fid,i.
// es (optional) = the STATISTICAL r error (19th column rErr_stat, 2026-09-14/15:
// the fit covariance conditioned on the constrained nuisances = the frozen-
// nuisance error without a refit, fork extract_pO_simfit.C); -1 when the
// column is absent. rErr itself is the TOTAL (profiled) error since the
// 2026-08-17 nuisances.
static bool readBinR(const char *csv, const char *charge, const char *binning,
                     double r[12], double e[12], double es[12] = nullptr)
{
    for (int i = 0; i < 12; ++i) { r[i] = 0; e[i] = 0; if (es) es[i] = -1.0; }
    std::ifstream in(csv);
    if (!in) { std::cerr << "[ERROR] cannot open CSV: " << csv << "\n"; return false; }
    std::string line;
    std::getline(in, line); // header
    int found = 0;
    while (std::getline(in, line))
    {
        std::stringstream ss(line);
        std::string f;
        std::vector<std::string> c;
        while (std::getline(ss, f, ',')) c.push_back(f);
        if (c.size() < 9) continue;
        if (c[1] == charge && c[2] == binning)
        {
            int iy = std::atoi(c[3].c_str());
            if (iy >= 0 && iy < 12)
            {
                r[iy] = std::atof(c[4].c_str()); e[iy] = std::atof(c[5].c_str()); ++found;
                if (es) es[iy] = (c.size() >= 19) ? std::atof(c[18].c_str()) : -1.0;
            }
        }
    }
    return found > 0;
}

// -----------------------------------------------------------------------------
// Ingredients shared by every simfit-based cross section, stacked in the
// h_cov_yield order [Wp_y0..11, Wm_y0..11]:
//   G    = gen FIDUCIAL sigma per bin, PER FLAVOUR (nb)  <- skim/gen_xsec.C
//   S    = prefit signal template integral, mu+e summed  <- comb_W_yields.csv
//   R,ER = fitted signal strength r and its TOTAL error  <- comb_W_yields.csv
//   covr = r covariance = cov_Y_ij/(S_i S_j)             <- h_cov_yield
//   ERs, covrS = the STATISTICAL part of the same (rErr_stat column,
//          h_cov_yield_stat; 2026-09-14/15) -- haveStat false when absent, then
//          ERs/covrS just mirror the total ones and the plots draw one bar.
struct CombIngredients
{
    double G[24], S[24], R[24], ER[24], covr[24][24];
    double ERs[24], covrS[24][24];
    bool   haveStat = false;
};

// The fiducial line drawn on the plots, from the discriminant's gen set.
static TString fidText(const char *disc)
{
    TString tag, label;
    pODisc::GenFid(disc, tag, label);
    return "fiducial: bare lepton " + label;
}

// `disc` picks the gen fiducial set (pODisc::GenFid -- since 2026-10-13 the
// leppt_mt40 fit is measured in the m_T > 40 fiducial volume, met in none).
static bool loadCombIngredients(const char *binsCsv, const char *yieldsRoot,
                                const char *genFile, const char *disc, CombIngredients &in)
{
    TString gtag, glabel;
    if (!pODisc::GenFid(disc, gtag, glabel)) return false;
    TFile *fG = TFile::Open(genFile, "READ");
    TH1D *hGp = (fG && !fG->IsZombie()) ? (TH1D *)fG->Get("h_gen_sig_Wp" + gtag) : nullptr;
    TH1D *hGm = (fG && !fG->IsZombie()) ? (TH1D *)fG->Get("h_gen_sig_Wm" + gtag) : nullptr;
    if (!hGp || !hGm || hGp->GetNbinsX() != 12 || hGm->GetNbinsX() != 12)
    {
        std::cerr << "[ERROR] no usable 12-bin h_gen_sig_{Wp,Wm}" << gtag << " in " << genFile
                  << " -- the r x sigma_gen extraction needs it (run skim/run_gen_xsec.sh).\n";
        if (fG) { fG->Close(); delete fG; }
        return false;
    }

    double rP[12], erP[12], rM[12], erM[12], sp[12], sm[12], esP[12], esM[12];
    if (!readBinR(binsCsv, "Wp", "lab", rP, erP, esP) || !readBinR(binsCsv, "Wm", "lab", rM, erM, esM))
    {
        std::cerr << "[ERROR] could not read per-bin r from " << binsCsv << "\n";
        fG->Close(); delete fG; return false;
    }
    const bool haveS = readBinPrefit(binsCsv, "Wp", "lab", sp) &&
                       readBinPrefit(binsCsv, "Wm", "lab", sm);
    in.haveStat = true;
    for (int i = 0; i < 12; ++i)
    {
        in.G[i] = hGp->GetBinContent(i + 1);  in.G[12 + i] = hGm->GetBinContent(i + 1);
        in.S[i] = haveS ? sp[i] : 0;          in.S[12 + i] = haveS ? sm[i] : 0;
        in.R[i] = rP[i]; in.ER[i] = erP[i];   in.R[12 + i] = rM[i]; in.ER[12 + i] = erM[i];
        in.ERs[i] = esP[i];                   in.ERs[12 + i] = esM[i];
        if (esP[i] < 0 || esM[i] < 0) in.haveStat = false;
    }
    if (!in.haveStat)
    {
        std::cerr << "[WARN] no stat-only r errors (rErr_stat column) in " << binsCsv
                  << " -> plots draw ONE bar = the total (profiled) fit uncertainty\n";
        for (int i = 0; i < 24; ++i) in.ERs[i] = in.ER[i];
    }
    fG->Close(); delete fG;

    // r covariance from the yield covariance; fallback = diagonal rErr^2
    // (sums then lose the cross terms). Same for the stat-only matrix.
    TFile *fY = TFile::Open(yieldsRoot, "READ");
    TH2D *hcov  = (fY && !fY->IsZombie()) ? (TH2D *)fY->Get("h_cov_yield") : nullptr;
    TH2D *hcovS = (fY && !fY->IsZombie() && in.haveStat) ? (TH2D *)fY->Get("h_cov_yield_stat") : nullptr;
    if (hcov && hcov->GetNbinsX() != 24) hcov = nullptr;
    if (hcovS && hcovS->GetNbinsX() != 24) hcovS = nullptr;
    if (!hcov || !haveS)
    {
        hcov = nullptr; hcovS = nullptr;
        std::cerr << "[WARN] no usable h_cov_yield in " << yieldsRoot
                  << (haveS ? "" : " (and/or no signal_prefit columns)")
                  << " -> r-covariance falls back to DIAGONAL rErr^2\n";
    }
    else if (in.haveStat && !hcovS)
        std::cerr << "[WARN] no h_cov_yield_stat in " << yieldsRoot
                  << " -> stat-only r-covariance falls back to DIAGONAL rErr_stat^2\n";
    for (int i = 0; i < 24; ++i)
        for (int j = 0; j < 24; ++j)
        {
            in.covr[i][j] = (hcov && in.S[i] > 0 && in.S[j] > 0)
                                ? hcov->GetBinContent(i + 1, j + 1) / (in.S[i] * in.S[j])
                                : (i == j ? in.ER[i] * in.ER[i] : 0.0);
            in.covrS[i][j] = !in.haveStat ? in.covr[i][j]
                             : (hcovS && in.S[i] > 0 && in.S[j] > 0)
                                   ? hcovS->GetBinContent(i + 1, j + 1) / (in.S[i] * in.S[j])
                                   : (i == j ? in.ERs[i] * in.ERs[i] : 0.0);
        }
    if (fY) { fY->Close(); delete fY; }
    return true;
}

// Sum_i w_i r_i over bins [lo, hi) with Var = Sum_ij w_i w_j cov_r,ij -- the
// one place the r-covariance is propagated (w = sigma_gen for the measurement,
// w = per-flavour prefit S for the count-based version). stat = true uses the
// stat-only covariance (identical to the total when no companion fit exists).
static void sumWithCov(const CombIngredients &in, const double *w, int lo, int hi,
                       double &val, double &err, bool stat = false)
{
    val = 0;
    double var = 0;
    for (int i = lo; i < hi; ++i)
    {
        val += w[i] * in.R[i];
        for (int j = lo; j < hi; ++j) var += w[i] * w[j] * (stat ? in.covrS[i][j] : in.covr[i][j]);
    }
    err = std::sqrt(std::max(0.0, var));
}

// syst = sqrt(total^2 - stat^2) (0 when no stat-only companion exists)
static double systFrom(double tot, double stat, bool haveStat)
{
    return haveStat ? std::sqrt(std::max(0.0, tot * tot - stat * stat)) : 0.0;
}

// =============================================================================
// xsec_fiducial_comb -- the simfit (grand simultaneous fit) analog (2026-08-05;
// r x sigma_gen EXTRACTION since 2026-08-12):
//
//     sigma_meas,i = r_i x sigma_gen-fid,i          [PER LEPTON FLAVOUR, nb]
//
// r_i = the mu/e-shared per-(charge, y-bin) POI (r/rErr of comb_W_yields.csv),
// sigma_gen-fid,i = the gen FIDUCIAL sigma from skim/gen_xsec.C (h_gen_sig_*:
// bare post-FSR lepton pT > 25, |eta_lab| < 2.4). Algebraically identical to
// N_fit/(L*(A*eps)_MC), but kA_O/kSigma_* cancel between r and sigma_gen and
// future template SFs propagate through r automatically. Inclusive sums
// Sum_i r_i*sigma_i carry the full r-covariance, recovered from h_cov_yield
// (yield space) via cov(r_i,r_j) = cov_Y_ij/(S_i*S_j), S = signal_prefit
// (diagonal-rErr fallback + WARN). Gen MC-stat error (~0.3%/bin) is neglected
// against the fit error (~6%/bin).
//
// PER-FLAVOUR convention since 2026-08-12: the result is sigma(W -> l nu) for
// ONE lepton flavour (lepton universality; the r's are mu/e-shared) -- no more
// "x2 = mu+e summed" display.
//
// The COUNT record is untouched: comb_W_yields.csv / comb_summary.csv keep the
// fitted yields (printed here as a console cross-check), and charge_asym /
// FBratio stay count-based by design. This macro additionally writes the sigma
// table plots/comb/xsec/<disc>/xsec_comb.csv (per-bin + inclusive).
// Outputs: plots/comb/xsec/<disc>/{W_fiducial, W_dsigma_deta_comb}.{png,pdf}
//          + xsec_comb.csv
//
// UNCERTAINTIES (2026-09-14/15): the fitted rErr is the TOTAL profiled error
// (stat (+) every nuisance in the fit: QCD lnN, lumi 3%, the theory shapes,
// the lepton SFs). The fork's extractor also writes the STATISTICAL part
// (CSV column rErr_stat + h_cov_yield_stat: the post-fit covariance
// conditioned on the constrained nuisances -- what a refit with them frozen
// at their post-fit values gives, without the refit; since 2026-09-15b the
// --statonly companion fit itself); the plots then draw the point + the
// STATISTICAL bar + the systematic as a TBox (MakeSystBoxes), and
// xsec_comb.csv splits the error into sigma_meas_stat_nb / sigma_meas_syst_nb
// (syst = sqrt(total^2 - stat^2)); without the column a single (total) bar is
// drawn and labelled as such.
//
//   root -l -q -e 'gROOT->LoadMacro("xsec_fiducial.C+"); xsec_fiducial_comb("met");'
// =============================================================================
void xsec_fiducial_comb(const char *disc = "met",
                        const char *summaryCsv = nullptr, // default: simfit/summary/comb_summary.csv
                        const char *binsCsv    = nullptr, // default: simfit/summary/comb_W_yields.csv
                        const char *chargeFile = nullptr, // default: charge_asym_fid_comb_<disc>.root (eta bin positions)
                        const char *yieldsRoot = nullptr, // default: simfit/summary/comb_fitted_yields.root
                        const char *genFile    = "../skim/rootfile/gen_xsec.root", // skim/gen_xsec.C output; the set (m_T cut or not) follows disc
                        double Lint = pONorm::kLumi_invnb)
{
    gStyle->SetEndErrorSize(4);

    TString dsuf, discLabel;
    if (!pODisc::Spec(disc, dsuf, discLabel)) return;
    pOFit::Spec comb;
    pOFit::Get("comb", comb);
    const TString base = pOFit::WorkDir(dsuf, comb) + "/summary";   // honours $FORK_TEST
    const TString sSum    = summaryCsv ? TString(summaryCsv) : base + "/comb_summary.csv";
    const TString sBins   = binsCsv    ? TString(binsCsv)    : base + "/comb_W_yields.csv";
    const TString sRoot   = yieldsRoot ? TString(yieldsRoot) : base + "/comb_fitted_yields.root";
    const TString sCharge = chargeFile ? TString(chargeFile) : pODisc::GraphFile("charge_asym", "comb", disc);

    const std::string outDir = std::string("./plots/comb/xsec/") + disc;
    gSystem->mkdir(outDir.c_str(), kTRUE);

    // ---- 0) ingredients for sigma = r x sigma_gen-fid ------------------------
    CombIngredients in;
    if (!loadCombIngredients(sBins.Data(), sRoot.Data(), genFile, disc, in)) return;
    const double *G = in.G, *R = in.R, *ER = in.ER;
    auto &covr = in.covr;

    // ---- 1) inclusive W+/W-/W: sigma = Sum_i r_i sigma_gen,i (cov-propagated)
    // ev = TOTAL (profiled) uncertainty, evS = statistical (stat-only companion
    // fit), syst = sqrt(total^2 - stat^2)
    const char *xlab[3] = {"W^{+}", "W^{-}", "W"};
    double yv[3] = {0, 0, 0}, ev[3] = {0, 0, 0}, evS[3] = {0, 0, 0}, genTot[3] = {0, 0, 0};
    const int klo[3] = {0, 12, 0}, khi[3] = {12, 24, 24};
    for (int k = 0; k < 3; ++k)
    {
        sumWithCov(in, G, klo[k], khi[k], yv[k], ev[k]);
        double dummy; sumWithCov(in, G, klo[k], khi[k], dummy, evS[k], true);
        for (int i = klo[k]; i < khi[k]; ++i) genTot[k] += G[i];
        if (in.haveStat)
            printf("[xsec-comb] %-6s measured = %6.2f +/- %4.2f (stat) +/- %4.2f (syst) = +/- %4.2f (total) nb   (gen fid %6.2f nb, eff. r = %.3f)  [per flavour]\n",
                   xlab[k], yv[k], evS[k], systFrom(ev[k], evS[k], true), ev[k], genTot[k], genTot[k] > 0 ? yv[k] / genTot[k] : 0.0);
        else
            printf("[xsec-comb] %-6s measured = %6.2f +/- %4.2f nb (stat+syst profiled)   (gen fid %6.2f nb, eff. r = %.3f)  [per flavour]\n",
                   xlab[k], yv[k], ev[k], genTot[k], genTot[k] > 0 ? yv[k] / genTot[k] : 0.0);
    }

    // console cross-check ONLY: the COUNT record (stays in the fork CSVs;
    // charge_asym / FBratio remain count-based).
    {
        const char *params[3] = {"Wp_sum_yield", "Wm_sum_yield", "W_sum_yield"};
        for (int i = 0; i < 3; ++i)
        {
            double y, e;
            if (readCombSum(sSum.Data(), "simfit_lab", params[i], y, e))
                printf("[xsec-comb] counts %-13s N/L = %6.2f +/- %4.2f nb (mu+e, uncorrected)\n",
                       params[i], y / Lint, e / Lint);
        }
    }

    {
        double x[3] = {1, 2, 3}, ex0[3] = {0, 0, 0};
        double xg[3] = {0.86, 1.86, 2.86}, ey0[3] = {0, 0, 0};
        TGraphErrors *g = new TGraphErrors(3, x, yv, ex0, ev);      // TOTAL error
        TGraphErrors *gS = new TGraphErrors(3, x, yv, ex0, evS);    // STATISTICAL error
        TGraphErrors *gGen = new TGraphErrors(3, xg, genTot, ex0, ey0);

        PlotStyle ps; ps.logy = false;
        // these three points sit at discrete x with no x error -> the syst
        // boxes need an absolute half-width (in x units of the 3-bin frame)
        ps.systBoxHalfWidthAbs = 0.13;
        TCanvas *c = new TCanvas("c_xsec_fid_comb", "", ps.w, ps.h);
        ApplyCanvasStyle(c, ps);
        c->cd();
        double ymax = 0.0;
        for (int i = 0; i < 3; ++i)
            ymax = std::max(ymax, std::max(yv[i] + ev[i], genTot[i]));
        TH1F *fr = new TH1F("fr_xsec_comb", "", 3, 0.5, 3.5);
        for (int i = 0; i < 3; ++i) fr->GetXaxis()->SetBinLabel(i + 1, xlab[i]);
        fr->SetStats(0);
        ApplyHistStyle(fr, ps, "", "#sigma^{fid}_{W#rightarrowl#nu} (nb)");
        fr->GetXaxis()->SetLabelSize(0.065);
        fr->GetXaxis()->SetNdivisions(3);
        fr->SetMinimum(0.0);
        fr->SetMaximum(1.45 * ymax);
        fr->Draw();

        gGen->SetMarkerStyle(27); gGen->SetMarkerSize(1.9); // open diamond = gen fid (r = 1)
        gGen->SetMarkerColor(kGreen + 2); gGen->SetLineColor(kGreen + 2); gGen->SetLineWidth(2);
        gGen->Draw("P SAME");
        // point + STATISTICAL bar, with the systematic -- the quadratic
        // difference from the total (profiled) error -- as a TBox per point.
        // Without the stat component (in.haveStat false) the single bar is the
        // total and no box is drawn.
        std::vector<TBox *> sboxes = MakeSystBoxes(g, in.haveStat ? gS : nullptr, ps);
        for (TBox *b : sboxes) b->Draw();
        TGraphErrors *gDraw = in.haveStat ? gS : g;
        gDraw->SetMarkerStyle(20); gDraw->SetMarkerSize(1.5);
        gDraw->SetMarkerColor(kBlack); gDraw->SetLineColor(kBlack); gDraw->SetLineWidth(2);
        gDraw->Draw("P SAME");

        DrawHeader(ps, "", "W #rightarrow l #nu", "fiducial cross section");
        TLatex ltag; ltag.SetNDC(); ltag.SetTextFont(ps.font); ltag.SetTextAlign(13);
        ltag.SetTextSize(ps.boxTextSize);
        ltag.DrawLatex(ps.headerX, ps.headerY - 3.0 * ps.headerDy, Form("%s fit", discLabel.Data()));

        TLegend *leg = new TLegend(0.20, 0.54, 0.56, 0.72);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.032);
        leg->AddEntry(gDraw, in.haveStat ? "measured, #sigma = r #times #sigma^{gen}_{fid} (bars: stat.)"
                                         : "measured, #sigma = r #times #sigma^{gen}_{fid}", "lep");
        if (!sboxes.empty()) leg->AddEntry(sboxes[0], "syst. (fit nuisances)", "f");
        leg->AddEntry(gGen, "#sigma^{gen}_{fid} (POWHEG, r = 1)", "p");
        leg->Draw();

        TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.030);
        tx.DrawLatex(0.18, 0.30,  "#mu + e combined (simfit), per lepton flavour");
        tx.DrawLatex(0.18, 0.255, fidText(disc));
        tx.DrawLatex(0.18, 0.21,  "#sigma_{i} = r_{i} #times #sigma^{gen}_{fid,i}, summed with r-covariance");
        tx.DrawLatex(0.18, 0.165, in.haveStat ? "syst. = profiled fit nuisances (QCD, lumi 3%, theory, lepton SF)"
                                              : "one bar = stat. #oplus syst. (all fit nuisances profiled)");

        CMS_lumi(c, 13, 10);
        c->RedrawAxis(); c->Modified(); c->Update();
        c->SaveAs((outDir + "/W_fiducial.png").c_str());
        c->SaveAs((outDir + "/W_fiducial.pdf").c_str());
        std::cout << "[OK] Saved " << outDir << "/W_fiducial.{png,pdf}  (disc=" << disc
                  << ", L = " << Lint << " nb^-1)\n";
    }

    // ---- 2) differential d(sigma)/d(eta_CM), mu+e summed ---------------------
    TFile *fC = TFile::Open(sCharge, "READ");
    TGraphErrors *gA = (fC && !fC->IsZombie()) ? (TGraphErrors *)fC->Get("g_chargeAsym") : nullptr;
    if (!gA && fC && !fC->IsZombie()) gA = (TGraphErrors *)fC->Get("g_chargeAsym_mt");
    if (!gA) { std::cerr << "[ERROR] need g_chargeAsym(_mt) in " << sCharge
                         << " for the CM eta bin positions (run charge_asym.C on the comb yields first).\n"; return; }
    const int NY = gA->GetN();
    std::vector<double> xc(NY), exc(NY);
    for (int i = 0; i < NY; ++i) { xc[i] = gA->GetPointX(i); exc[i] = gA->GetErrorX(i); }

    if (NY != 12)
    { std::cerr << "[ERROR] charge-asym graph has " << NY << " points (expected 12) -- binning mismatch.\n"; return; }

    // measured d(sigma)/d(eta): r_i * sigma_gen,i / dEta (per flavour);
    // W-inclusive per bin carries cov(r_Wp_i, r_Wm_i) from the r-covariance.
    // (es* = total errors from covr; esS* = stat-only from covrS)
    std::vector<double> sP(NY), esP(NY), sM(NY), esM(NY), sI(NY), esI(NY);
    std::vector<double> esSP(NY), esSM(NY), esSI(NY);
    std::vector<double> gP2(NY, 0), gM2(NY, 0), gI2(NY, 0), zeroe(NY, 0.0);
    const auto &covrS = in.covrS;
    double ymax = 0.0;
    for (int i = 0; i < NY; ++i)
    {
        double dEta = 2.0 * exc[i];
        if (dEta <= 0) dEta = 0.4;
        const double gp = G[i], gm = G[12 + i];
        sP[i] = R[i] * gp / dEta;       esP[i] = std::sqrt(covr[i][i]) * gp / dEta;
        sM[i] = R[12 + i] * gm / dEta;  esM[i] = std::sqrt(covr[12 + i][12 + i]) * gm / dEta;
        esSP[i] = std::sqrt(covrS[i][i]) * gp / dEta;
        esSM[i] = std::sqrt(covrS[12 + i][12 + i]) * gm / dEta;
        double vi = gp * gp * covr[i][i] + gm * gm * covr[12 + i][12 + i]
                  + 2.0 * gp * gm * covr[i][12 + i];
        double viS = gp * gp * covrS[i][i] + gm * gm * covrS[12 + i][12 + i]
                   + 2.0 * gp * gm * covrS[i][12 + i];
        if (vi < 0) vi = 0;
        if (viS < 0) viS = 0;
        sI[i]  = (R[i] * gp + R[12 + i] * gm) / dEta;
        esI[i] = std::sqrt(vi) / dEta;
        esSI[i] = std::sqrt(viS) / dEta;
        // gen fiducial (r = 1 expectation), dashed lines
        gP2[i] = gp / dEta;
        gM2[i] = gm / dEta;
        gI2[i] = gP2[i] + gM2[i];
        ymax = std::max(ymax, std::max(sI[i] + esI[i], gI2[i]));
    }

    TGraphErrors *gI = new TGraphErrors(NY, xc.data(), sI.data(), exc.data(), esI.data());
    TGraphErrors *gP = new TGraphErrors(NY, xc.data(), sP.data(), exc.data(), esP.data());
    TGraphErrors *gM = new TGraphErrors(NY, xc.data(), sM.data(), exc.data(), esM.data());
    // statistical twins: same points, stat error, and the bin-width x error kept
    // so MakeSystBoxes can size the systematic boxes from it
    TGraphErrors *gIS = new TGraphErrors(NY, xc.data(), sI.data(), exc.data(), esSI.data());
    TGraphErrors *gPS = new TGraphErrors(NY, xc.data(), sP.data(), exc.data(), esSP.data());
    TGraphErrors *gMS = new TGraphErrors(NY, xc.data(), sM.data(), exc.data(), esSM.data());
    TGraphErrors *gIg = new TGraphErrors(NY, xc.data(), gI2.data(), zeroe.data(), zeroe.data());
    TGraphErrors *gPg = new TGraphErrors(NY, xc.data(), gP2.data(), zeroe.data(), zeroe.data());
    TGraphErrors *gMg = new TGraphErrors(NY, xc.data(), gM2.data(), zeroe.data(), zeroe.data());

    PlotStyle ps; ps.logy = false;
    TCanvas *c = new TCanvas("c_dsig_comb", "", ps.w, ps.h);
    ApplyCanvasStyle(c, ps);
    c->cd();

    const double xlo = xc[0] - exc[0] - 0.1, xhi = xc[NY - 1] + exc[NY - 1] + 0.1;
    TH1F *fr = new TH1F("fr_dsig_comb", "", 100, xlo, xhi);
    fr->SetStats(0);
    ApplyHistStyle(fr, ps, "#eta^{l}_{CM}", "d#sigma^{fid}_{W#rightarrowl#nu}/d#eta (nb)");
    fr->SetMinimum(0.0);
    // Annotations live in the upper band of the frame (NDC y >= kAnnoY): the
    // legend top-left and the three info lines on the right under the header.
    // The frame maximum is raised so that every point (+ error) and gen line
    // stays below that band -- with a fixed 1.55 x ymax the 2026-09-14 values
    // (W incl. 22-25 nb) put the 4th legend entry on the eta_CM = -2.1 ... -1.4
    // points, and info lines at the bottom right always collided with the W-
    // points and the falling gen line at forward eta.
    const double kAnnoY = 0.55;                                     // lower edge of the annotation band (NDC)
    const double legX1 = 0.18, legY1 = kAnnoY, legX2 = 0.48, legY2 = 0.79;
    const double fAnno = (kAnnoY - ps.bm) / (1.0 - ps.bm - ps.tm);  // that edge as a fraction of the y axis
    double yMaxFrame = 1.55 * ymax;
    if (fAnno > 0.05) yMaxFrame = std::max(yMaxFrame, 1.06 * ymax / fAnno);
    fr->SetMaximum(yMaxFrame);
    fr->Draw();

    // gen-fiducial dashed lines first (bottom layer) = the r = 1 expectation;
    // measured point / dashed line per bin IS the fitted r_i.
    gIg->SetLineStyle(2); gIg->SetLineWidth(3); gIg->SetLineColor(kBlack);
    gPg->SetLineStyle(2); gPg->SetLineWidth(3); gPg->SetLineColor(kAzure + 2);
    gMg->SetLineStyle(2); gMg->SetLineWidth(3); gMg->SetLineColor(kRed + 1);
    gIg->Draw("L SAME"); gPg->Draw("L SAME"); gMg->Draw("L SAME");

    // Point + STATISTICAL bar, systematic as a TBox per point in the series
    // colour. The three series share the x positions but sit at different y
    // (W ~ W+ + W-), so their boxes never overlap.
    //
    // MARKER CHOICE (measured 2026-09-15, see kMk* in the header comment):
    // point-symmetric glyphs only -- a TRIANGLE (22/26) has its visual centre
    // ~1/3 of the height above the base, so the eye reads the point low, and
    // ROOT rasterizes it 1.0-1.5 px off the error bar. Sizes are the ones that
    // rasterize CENTRED on the bar: square 21 at 1.3 lands 1.0 px left and
    // 1.5 px high (this is the offset that showed up in the legend), at 1.5 it
    // is exact in x and 0.5 px in y.
    gI->SetMarkerStyle(kMkW);  gI->SetMarkerSize(kMkWSz);  gI->SetMarkerColor(kBlack);     gI->SetLineColor(kBlack);     gI->SetLineWidth(2);
    gP->SetMarkerStyle(kMkWp); gP->SetMarkerSize(kMkWpSz); gP->SetMarkerColor(kAzure + 2); gP->SetLineColor(kAzure + 2); gP->SetLineWidth(2);
    gM->SetMarkerStyle(kMkWm); gM->SetMarkerSize(kMkWmSz); gM->SetMarkerColor(kRed + 1);   gM->SetLineColor(kRed + 1);   gM->SetLineWidth(2);
    gIS->SetMarkerStyle(kMkW);  gIS->SetMarkerSize(kMkWSz);  gIS->SetMarkerColor(kBlack);     gIS->SetLineColor(kBlack);     gIS->SetLineWidth(2);
    gPS->SetMarkerStyle(kMkWp); gPS->SetMarkerSize(kMkWpSz); gPS->SetMarkerColor(kAzure + 2); gPS->SetLineColor(kAzure + 2); gPS->SetLineWidth(2);
    gMS->SetMarkerStyle(kMkWm); gMS->SetMarkerSize(kMkWmSz); gMS->SetMarkerColor(kRed + 1);   gMS->SetLineColor(kRed + 1);   gMS->SetLineWidth(2);

    std::vector<TBox *> boxI, boxP, boxM;
    if (in.haveStat)
    {
        PlotStyle pb = ps;
        pb.systBoxWidthFrac = 0.60;
        pb.systBoxFillAlpha = 0.30;
        pb.systBoxFillColor = kBlack;     pb.systBoxLineColor = kBlack;     boxI = MakeSystBoxes(gI, gIS, pb);
        pb.systBoxFillColor = kAzure + 2; pb.systBoxLineColor = kAzure + 2; boxP = MakeSystBoxes(gP, gPS, pb);
        pb.systBoxFillColor = kRed + 1;   pb.systBoxLineColor = kRed + 1;   boxM = MakeSystBoxes(gM, gMS, pb);
        for (TBox *b : boxI) b->Draw();
        for (TBox *b : boxP) b->Draw();
        for (TBox *b : boxM) b->Draw();
        gIS->Draw("P SAME"); gPS->Draw("P SAME"); gMS->Draw("P SAME");
    }
    else
    {
        gI->Draw("P SAME"); gP->Draw("P SAME"); gM->Draw("P SAME");
    }

    DrawHeader(ps, "", "W #rightarrow l #nu", "fiducial cross section");
    TLatex ltag2; ltag2.SetNDC(); ltag2.SetTextFont(ps.font); ltag2.SetTextAlign(13);
    ltag2.SetTextSize(ps.boxTextSize);
    ltag2.DrawLatex(ps.headerX, ps.headerY - 3.0 * ps.headerDy, Form("%s fit", discLabel.Data()));

    // info lines: right column under the header, inside the annotation band
    // (x = 0.55 keeps them clear of the legend's longest entry on the same row)
    TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.029);
    tx.DrawLatex(0.55, kAnnoY + 0.095, "#mu + e combined (simfit), per flavour");
    tx.DrawLatex(0.55, kAnnoY + 0.050, "#sigma_{i} = r_{i} #times #sigma^{gen}_{fid,i} (p_{T}^{l} > 25 GeV)");
    tx.DrawLatex(0.55, kAnnoY + 0.005, in.haveStat ? "bars: stat.,  boxes: syst. (profiled)"
                                                    : "one bar: stat. #oplus syst. (profiled)");

    TLegend *leg = new TLegend(legX1, legY1, legX2, legY2);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.031);
    leg->AddEntry(in.haveStat ? gIS : gI, "W (incl.) measured", "lep");
    leg->AddEntry(in.haveStat ? gPS : gP, "W^{+} measured", "lep");
    leg->AddEntry(in.haveStat ? gMS : gM, "W^{-} measured", "lep");
    leg->AddEntry(gIg, "#sigma^{gen}_{fid} (POWHEG, r = 1)", "l");  // dashed sample shows the style
    leg->Draw();

    CMS_lumi(c, 13, 10);
    c->RedrawAxis(); c->Modified(); c->Update();
    c->SaveAs((outDir + "/W_dsigma_deta_comb.png").c_str());
    c->SaveAs((outDir + "/W_dsigma_deta_comb.pdf").c_str());
    std::cout << "[OK] Saved " << outDir << "/W_dsigma_deta_comb.{png,pdf}  (disc="
              << disc << ", L = " << Lint << " nb^-1)\n";

    // ---- 3) sigma table: the r x sigma_gen record (per-bin + inclusive) ------
    // (the COUNT record stays in comb_W_yields.csv / comb_summary.csv, fork.)
    {
        // sigma_meas_err_nb = TOTAL (profiled); the two trailing columns split
        // it into stat (the fit covariance conditioned on the constrained
        // nuisances) and syst = sqrt(total^2 - stat^2); both -1 when the
        // fitted-yields file carries no h_cov_yield_stat / rErr_stat.
        const std::string csvPath = outDir + "/xsec_comb.csv";
        std::ofstream out(csvPath.c_str());
        out << "charge,ybin,etaCM_center,etaCM_halfwidth,r,rErr,sigma_gen_fid_nb,sigma_meas_nb,sigma_meas_err_nb,"
               "sigma_meas_stat_nb,sigma_meas_syst_nb\n";
        const char *cn2[2] = {"Wp", "Wm"};
        for (int q = 0; q < 2; ++q)
            for (int i = 0; i < NY; ++i)
            {
                const int k = 12 * q + i;
                const double st = in.haveStat ? in.ERs[k] * G[k] : -1.0;
                out << Form("%s,%d,%.4f,%.4f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n",
                            cn2[q], i, xc[i], exc[i], R[k], ER[k],
                            G[k], R[k] * G[k], ER[k] * G[k],
                            st, in.haveStat ? systFrom(ER[k] * G[k], st, true) : -1.0);
            }
        const char *cn3[3] = {"Wp", "Wm", "W"};
        for (int k = 0; k < 3; ++k)  // r column = effective (gen-weighted mean) r
            out << Form("%s,incl,,,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", cn3[k],
                        genTot[k] > 0 ? yv[k] / genTot[k] : 0.0,
                        genTot[k] > 0 ? ev[k] / genTot[k] : 0.0,
                        genTot[k], yv[k], ev[k],
                        in.haveStat ? evS[k] : -1.0,
                        in.haveStat ? systFrom(ev[k], evS[k], true) : -1.0);
        std::cout << "[OK] Wrote " << csvPath << "\n";
    }

    if (fC) { fC->Close(); delete fC; }
}

// -----------------------------------------------------------------------------
// PER-FLAVOUR prefit signal per lab bin from THIS repo's structured Combine
// inputs (plotting/mtandmet.C): region dirs W{p,m}_lab_y{i}, absolute `signal`
// template. mu + e reproduce the CSV signal_prefit column exactly (checked).
// Filled into the stacked [Wp_y0..11, Wm_y0..11] order.
static bool readFlavourPrefit(const char *file, double s[24])
{
    TFile *f = TFile::Open(file, "READ");
    if (!f || f->IsZombie())
    { std::cerr << "[WARN] cannot open " << file << "\n"; if (f) delete f; return false; }
    bool ok = true;
    const char *cn[2] = {"Wp", "Wm"};
    for (int q = 0; q < 2 && ok; ++q)
        for (int i = 0; i < 12; ++i)
        {
            TH1 *h = (TH1 *)f->Get(Form("%s_lab_y%d/signal", cn[q], i));
            if (!h) { ok = false; break; }
            s[12 * q + i] = h->Integral();
        }
    if (!ok) std::cerr << "[WARN] missing W{p,m}_lab_y*/signal in " << file << "\n";
    f->Close(); delete f;
    return ok;
}

// =============================================================================
// xsec_fiducial_diag -- PER-FLAVOUR DIAGNOSTIC of the cross-section extraction
// (2026-08-12). One plot per lepton flavour, each with four series per charge:
//
//   sigma^gen        = sigma_gen-fid                  (POWHEG, r = 1)
//   sigma^reco       = S_f / L                        (MC reco expectation, r = 1)
//   r x sigma^gen    = Sum_i r_i sigma_gen,i          (THE measurement -- identical
//                                                      in both plots: r is shared)
//   r x sigma^reco   = Sum_i r_i S_f,i / L = N_fit,f/L (count-based, uncorrected)
//
// The vertical gap within each pair is the SAME factor (A x eps)_MC,f, so the
// plot reads directly as "how big is the MC correction, and does it differ
// between mu and e" -- the acceptance x efficiency the r x sigma_gen extraction
// applies implicitly. Errors on the two r-based series use the full
// r-covariance; the r = 1 series are MC and drawn without errors.
// Sum over flavours of r x sigma^reco reproduces the fork's *_sum_yield/L
// (printed as a cross-check).
//
// The same four series are drawn ETA-BINNED as d(sigma)/d(eta_CM), one canvas
// per (flavour, charge), plus a summary of (A x eps)(eta) for both flavours --
// where an ECAL-crack acceptance hole would show up as an electron-only dip.
//
// Purely diagnostic: the physics result stays in xsec_fiducial_comb.
// Outputs: plots/comb/xsec/<disc>/
//            W_xsec_diag_{mu,ele}.{png,pdf}                    (inclusive)
//            W_dsigma_diag_{mu,ele}_{Wp,Wm,W}.{png,pdf}        (eta-binned)
//            AxEps_vs_eta.{png,pdf}                            (summary)
//            xsec_diag.csv, xsec_diag_bins.csv
//
//   root -l -q -e 'gROOT->LoadMacro("xsec_fiducial.C+"); xsec_fiducial_diag("met");'
// =============================================================================
void xsec_fiducial_diag(const char *disc = "met",
                        const char *binsCsv    = nullptr, // default: simfit/summary/comb_W_yields.csv
                        const char *yieldsRoot = nullptr, // default: simfit/summary/comb_fitted_yields.root
                        const char *chargeFile = nullptr, // default: charge_asym_fid_comb_<disc>.root (eta bin positions)
                        const char *genFile    = "../skim/rootfile/gen_xsec.root",
                        double Lint = pONorm::kLumi_invnb)
{
    gStyle->SetEndErrorSize(4);

    TString dsuf, discLabel;
    if (!pODisc::Spec(disc, dsuf, discLabel)) return;
    pOFit::Spec comb;
    pOFit::Get("comb", comb);
    const TString base = pOFit::WorkDir(dsuf, comb) + "/summary";   // honours $FORK_TEST
    const TString sBins = binsCsv    ? TString(binsCsv)    : base + "/comb_W_yields.csv";
    const TString sRoot = yieldsRoot ? TString(yieldsRoot) : base + "/comb_fitted_yields.root";
    const TString sCharge = chargeFile ? TString(chargeFile) : pODisc::GraphFile("charge_asym", "comb", disc);

    const std::string outDir = std::string("./plots/comb/xsec/") + disc;
    gSystem->mkdir(outDir.c_str(), kTRUE);

    CombIngredients in;
    if (!loadCombIngredients(sBins.Data(), sRoot.Data(), genFile, disc, in)) return;

    // eta_CM bin centers/widths for the eta-binned part (same source as
    // xsec_fiducial_comb: the charge-asymmetry graph); missing -> inclusive only.
    TFile *fC = TFile::Open(sCharge, "READ");
    TGraphErrors *gA = (fC && !fC->IsZombie()) ? (TGraphErrors *)fC->Get("g_chargeAsym") : nullptr;
    if (!gA && fC && !fC->IsZombie()) gA = (TGraphErrors *)fC->Get("g_chargeAsym_mt");
    const int NY = (gA && gA->GetN() == 12) ? 12 : 0;
    std::vector<double> xc(12, 0), exc(12, 0);
    if (NY == 12)
        for (int i = 0; i < 12; ++i) { xc[i] = gA->GetPointX(i); exc[i] = gA->GetErrorX(i); }
    else
        std::cerr << "[WARN] no 12-point g_chargeAsym(_mt) in " << sCharge
                  << " -> eta-binned diagnostic skipped (run charge_asym.C on the comb yields first)\n";

    const char *flav[2]   = {"mu", "ele"};
    const char *lepSym[2] = {"#mu", "e"};
    // this repo's Combine inputs: muon in plots/, electron in plots/Elec/
    const TString inFile[2] = {TString::Format("./plots/combine_input_W%s.root", dsuf.Data()),
                               TString::Format("./plots/Elec/combine_input_W%s.root", dsuf.Data())};
    const char *xlab[3]  = {"W^{+}", "W^{-}", "W"};
    const int klo[3] = {0, 12, 0}, khi[3] = {12, 24, 24};

    std::ofstream csv((outDir + "/xsec_diag.csv").c_str());
    csv << "flavour,charge,sigma_gen_nb,sigma_reco_nb,r_x_gen_nb,r_x_gen_err_nb,"
           "r_x_reco_nb,r_x_reco_err_nb,AxEps_reco_over_gen\n";
    std::ofstream csvB((outDir + "/xsec_diag_bins.csv").c_str());
    csvB << "flavour,charge,ybin,etaCM_center,etaCM_halfwidth,dsigma_gen,dsigma_reco,"
            "dsigma_r_gen,dsigma_r_gen_err,dsigma_r_reco,dsigma_r_reco_err,AxEps\n";
    double sumCount[3] = {0, 0, 0}; // mu + e of r x sigma^reco, for the cross-check
    double axeps[2][3][12];         // (A x eps)(eta) per flavour/charge, for the summary plot
    for (int a = 0; a < 2; ++a) for (int b = 0; b < 3; ++b) for (int cc = 0; cc < 12; ++cc) axeps[a][b][cc] = 0;

    // per-bin value with the covariance-correct error, charge k (0=Wp, 1=Wm, 2=W)
    auto binVal = [&](const double *wt, int k, int i, double &v, double &e)
    {
        double tmp[24];
        for (int t = 0; t < 24; ++t) tmp[t] = 0.0;
        if (k != 1) tmp[i]      = wt[i];
        if (k != 0) tmp[12 + i] = wt[12 + i];
        sumWithCov(in, tmp, 0, 24, v, e);
    };

    for (int fl = 0; fl < 2; ++fl)
    {
        double Sf[24], w[24];
        if (!readFlavourPrefit(inFile[fl].Data(), Sf))
        {
            std::cerr << "[ERROR] " << flav[fl] << " diagnostic skipped -- run plotting/mtandmet.C"
                      << (dsuf.Length() ? " for this discriminant" : "") << " first.\n";
            continue;
        }
        for (int i = 0; i < 24; ++i) w[i] = Sf[i] / Lint; // reco sigma per bin (nb)

        double gen[3] = {0, 0, 0}, reco[3] = {0, 0, 0};
        double rg[3], erg[3], rr[3], err[3];
        for (int k = 0; k < 3; ++k)
        {
            for (int i = klo[k]; i < khi[k]; ++i) { gen[k] += in.G[i]; reco[k] += w[i]; }
            sumWithCov(in, in.G, klo[k], khi[k], rg[k], erg[k]);
            sumWithCov(in, w,    klo[k], khi[k], rr[k], err[k]);
            sumCount[k] += rr[k];
            csv << Form("%s,%s,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n",
                        flav[fl], (k == 0 ? "Wp" : k == 1 ? "Wm" : "W"),
                        gen[k], reco[k], rg[k], erg[k], rr[k], err[k],
                        gen[k] > 0 ? reco[k] / gen[k] : 0.0);
            printf("[xsec-diag/%-3s] %-6s gen %6.2f | reco %6.2f (A*eps = %.3f) | "
                   "r*gen %6.2f +/- %4.2f | r*reco %6.2f +/- %4.2f nb\n",
                   flav[fl], xlab[k], gen[k], reco[k], gen[k] > 0 ? reco[k] / gen[k] : 0.0,
                   rg[k], erg[k], rr[k], err[k]);
        }

        // ---- draw ------------------------------------------------------------
        double xGen[3], xReco[3], xMeas[3], xCnt[3], ex0[3] = {0, 0, 0}, ey0[3] = {0, 0, 0};
        for (int k = 0; k < 3; ++k)
        { xGen[k] = k + 1 - 0.27; xMeas[k] = k + 1 - 0.09; xReco[k] = k + 1 + 0.09; xCnt[k] = k + 1 + 0.27; }
        TGraphErrors *gGen  = new TGraphErrors(3, xGen,  gen,  ex0, ey0);
        TGraphErrors *gMeas = new TGraphErrors(3, xMeas, rg,   ex0, erg);
        TGraphErrors *gReco = new TGraphErrors(3, xReco, reco, ex0, ey0);
        TGraphErrors *gCnt  = new TGraphErrors(3, xCnt,  rr,   ex0, err);

        PlotStyle ps; ps.logy = false;
        TCanvas *c = new TCanvas(Form("c_xsec_diag_%s", flav[fl]), "", ps.w, ps.h);
        ApplyCanvasStyle(c, ps);
        c->cd();
        double ymax = 0.0;
        for (int k = 0; k < 3; ++k)
            ymax = std::max(ymax, std::max(std::max(gen[k], reco[k]),
                                           std::max(rg[k] + erg[k], rr[k] + err[k])));
        TH1F *fr = new TH1F(Form("fr_diag_%s", flav[fl]), "", 3, 0.5, 3.5);
        for (int k = 0; k < 3; ++k) fr->GetXaxis()->SetBinLabel(k + 1, xlab[k]);
        fr->SetStats(0);
        ApplyHistStyle(fr, ps, "", Form("#sigma^{fid}_{W#rightarrow%s#nu} (nb)", lepSym[fl]));
        fr->GetXaxis()->SetLabelSize(0.065);
        fr->GetXaxis()->SetNdivisions(3);
        fr->SetMinimum(0.0);
        fr->SetMaximum(1.55 * ymax);
        fr->Draw();

        gGen->SetMarkerStyle(27); gGen->SetMarkerSize(2.0);  // open diamond  = gen, r = 1
        gGen->SetMarkerColor(kGreen + 2); gGen->SetLineColor(kGreen + 2); gGen->SetLineWidth(2);
        gMeas->SetMarkerStyle(20); gMeas->SetMarkerSize(1.4); // filled circle = r x gen (measured)
        gMeas->SetMarkerColor(kBlack); gMeas->SetLineColor(kBlack); gMeas->SetLineWidth(2);
        gReco->SetMarkerStyle(25); gReco->SetMarkerSize(1.5); // open square   = reco MC, r = 1 (1.5 = the centred size)
        gReco->SetMarkerColor(kAzure + 2); gReco->SetLineColor(kAzure + 2); gReco->SetLineWidth(2);
        gCnt->SetMarkerStyle(21); gCnt->SetMarkerSize(1.5);   // filled square = r x reco (counts; 1.5 = the centred size, see kMk*)
        gCnt->SetMarkerColor(kAzure + 2); gCnt->SetLineColor(kAzure + 2); gCnt->SetLineWidth(2);
        gGen->Draw("P SAME"); gReco->Draw("P SAME"); gMeas->Draw("P SAME"); gCnt->Draw("P SAME");

        // (A x eps) label next to each gen/reco pair -- the size of the correction
        TLatex te; te.SetTextFont(42); te.SetTextSize(0.026); te.SetTextAlign(21);
        te.SetTextColor(kAzure + 2);
        for (int k = 0; k < 3; ++k)
            te.DrawLatex(k + 1 + 0.18, reco[k] + 0.045 * ymax,
                         Form("A#times#varepsilon = %.2f", gen[k] > 0 ? reco[k] / gen[k] : 0.0));

        DrawHeader(ps, "", Form("W #rightarrow %s #nu", lepSym[fl]), "extraction diagnostic");
        TLatex ltag; ltag.SetNDC(); ltag.SetTextFont(ps.font); ltag.SetTextAlign(13);
        ltag.SetTextSize(ps.boxTextSize);
        ltag.DrawLatex(ps.headerX, ps.headerY - 3.0 * ps.headerDy, Form("%s fit", discLabel.Data()));

        TLegend *leg = new TLegend(0.18, 0.50, 0.60, 0.72);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.030);
        leg->AddEntry(gGen,  "#sigma^{gen}_{fid} (POWHEG, r = 1)", "p");
        leg->AddEntry(gMeas, "r #times #sigma^{gen}_{fid}  (measured)", "lep");
        leg->AddEntry(gReco, "#sigma^{reco}_{MC} = S/L (r = 1)", "p");
        leg->AddEntry(gCnt,  "r #times #sigma^{reco} = N_{fit}/L  (counts)", "lep");
        leg->Draw();

        TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.029);
        tx.DrawLatex(0.18, 0.255, Form("simfit (r shared #mu/e), %s channel only", lepSym[fl]));
        tx.DrawLatex(0.18, 0.215, fidText(disc));
        tx.DrawLatex(0.18, 0.175, "stat. (fit) unc. only");

        CMS_lumi(c, 13, 10);
        c->RedrawAxis(); c->Modified(); c->Update();
        c->SaveAs(Form("%s/W_xsec_diag_%s.png", outDir.c_str(), flav[fl]));
        c->SaveAs(Form("%s/W_xsec_diag_%s.pdf", outDir.c_str(), flav[fl]));
        std::cout << "[OK] Saved " << outDir << "/W_xsec_diag_" << flav[fl] << ".{png,pdf}\n";

        // ---- eta-binned: the same four series as d(sigma)/d(eta_CM) ----------
        if (NY != 12) continue;
        const char *cn[3] = {"Wp", "Wm", "W"};
        for (int k = 0; k < 3; ++k)
        {
            std::vector<double> dg(NY), dr(NY), dm(NY), edm(NY), dc(NY), edc(NY), zero(NY, 0.0);
            double ymaxD = 0;
            for (int i = 0; i < NY; ++i)
            {
                double dEta = 2.0 * exc[i];
                if (dEta <= 0) dEta = 0.4;
                double tg = 0, tr = 0;
                if (k != 1) { tg += in.G[i];      tr += w[i]; }
                if (k != 0) { tg += in.G[12 + i]; tr += w[12 + i]; }
                double v, e;
                binVal(in.G, k, i, v, e); dm[i] = v / dEta; edm[i] = e / dEta;
                binVal(w,    k, i, v, e); dc[i] = v / dEta; edc[i] = e / dEta;
                dg[i] = tg / dEta;  dr[i] = tr / dEta;
                axeps[fl][k][i] = tg > 0 ? tr / tg : 0.0;
                ymaxD = std::max(ymaxD, std::max(dg[i], dm[i] + edm[i]));
                csvB << Form("%s,%s,%d,%.4f,%.4f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n",
                             flav[fl], cn[k], i, xc[i], exc[i], dg[i], dr[i],
                             dm[i], edm[i], dc[i], edc[i], axeps[fl][k][i]);
            }

            TGraph *ggen = new TGraph(NY, xc.data(), dg.data());
            TGraph *grec = new TGraph(NY, xc.data(), dr.data());
            TGraphErrors *gmea = new TGraphErrors(NY, xc.data(), dm.data(), exc.data(), edm.data());
            TGraphErrors *gcnt = new TGraphErrors(NY, xc.data(), dc.data(), exc.data(), edc.data());

            PlotStyle ps2; ps2.logy = false;
            TCanvas *cd = new TCanvas(Form("c_diag_deta_%s_%s", flav[fl], cn[k]), "", ps2.w, ps2.h);
            ApplyCanvasStyle(cd, ps2);
            cd->cd();
            const double xlo = xc[0] - exc[0] - 0.1, xhi = xc[NY - 1] + exc[NY - 1] + 0.1;
            TH1F *frD = new TH1F(Form("fr_diag_deta_%s_%s", flav[fl], cn[k]), "", 100, xlo, xhi);
            frD->SetStats(0);
            ApplyHistStyle(frD, ps2, "#eta^{l}_{CM}",
                           Form("d#sigma^{fid}_{W#rightarrow%s#nu}/d#eta (nb)", lepSym[fl]));
            frD->SetMinimum(0.0);
            frD->SetMaximum(1.60 * ymaxD);
            frD->Draw();

            ggen->SetLineStyle(2); ggen->SetLineWidth(3); ggen->SetLineColor(kGreen + 2);
            grec->SetLineStyle(3); grec->SetLineWidth(3); grec->SetLineColor(kAzure + 2);
            ggen->Draw("L SAME"); grec->Draw("L SAME");
            gmea->SetMarkerStyle(20); gmea->SetMarkerSize(1.3); gmea->SetMarkerColor(kBlack);
            gmea->SetLineColor(kBlack); gmea->SetLineWidth(2);
            gcnt->SetMarkerStyle(21); gcnt->SetMarkerSize(1.5); gcnt->SetMarkerColor(kAzure + 2); // 1.5 = the centred square size
            gcnt->SetLineColor(kAzure + 2); gcnt->SetLineWidth(2);
            gmea->Draw("P SAME"); gcnt->Draw("P SAME");

            const char *procD = (k == 0) ? "W^{+} #rightarrow %s #nu"
                              : (k == 1) ? "W^{-} #rightarrow %s #nu"
                                         : "W #rightarrow %s #nu";
            DrawHeader(ps2, "", Form(procD, lepSym[fl]), "extraction diagnostic");
            TLatex ltD; ltD.SetNDC(); ltD.SetTextFont(ps2.font); ltD.SetTextAlign(13);
            ltD.SetTextSize(ps2.boxTextSize);
            ltD.DrawLatex(ps2.headerX, ps2.headerY - 3.0 * ps2.headerDy, Form("%s fit", discLabel.Data()));

            TLegend *legD = new TLegend(0.18, 0.55, 0.60, 0.78);
            legD->SetBorderSize(0); legD->SetFillStyle(0); legD->SetTextFont(42); legD->SetTextSize(0.030);
            legD->AddEntry(ggen, "#sigma^{gen}_{fid} (POWHEG, r = 1)", "l");
            legD->AddEntry(gmea, "r #times #sigma^{gen}_{fid}  (measured)", "lep");
            legD->AddEntry(grec, "#sigma^{reco}_{MC} = S/L (r = 1)", "l");
            legD->AddEntry(gcnt, "r #times #sigma^{reco} = N_{fit}/L  (counts)", "lep");
            legD->Draw();

            TLatex txD; txD.SetNDC(); txD.SetTextFont(42); txD.SetTextSize(0.029);
            txD.DrawLatex(0.18, 0.255, Form("simfit (r shared #mu/e), %s channel only", lepSym[fl]));
            txD.DrawLatex(0.18, 0.215, "gen #rightarrow reco gap = (A#times#varepsilon)_{MC}");
            txD.DrawLatex(0.18, 0.175, "stat. (fit) unc. only");

            CMS_lumi(cd, 13, 10);
            cd->RedrawAxis(); cd->Modified(); cd->Update();
            cd->SaveAs(Form("%s/W_dsigma_diag_%s_%s.png", outDir.c_str(), flav[fl], cn[k]));
            cd->SaveAs(Form("%s/W_dsigma_diag_%s_%s.pdf", outDir.c_str(), flav[fl], cn[k]));
        }
        std::cout << "[OK] Saved " << outDir << "/W_dsigma_diag_" << flav[fl]
                  << "_{Wp,Wm,W}.{png,pdf}\n";
    }

    // ---- (A x eps)(eta) summary: mu vs e, W+ and W- --------------------------
    // The MC correction the r x sigma_gen extraction applies, bin by bin. An
    // electron-only dip is the ECAL-crack acceptance hole (|eta_SC| 1.44-1.57).
    if (NY == 12 && axeps[0][0][0] > 0 && axeps[1][0][0] > 0)
    {
        PlotStyle ps3; ps3.logy = false;
        TCanvas *ca = new TCanvas("c_axeps", "", ps3.w, ps3.h);
        ApplyCanvasStyle(ca, ps3);
        ca->cd();
        double vlo = 1e9, vhi = 0;
        for (int fl = 0; fl < 2; ++fl)
            for (int k = 0; k < 2; ++k)
                for (int i = 0; i < NY; ++i)
                { vlo = std::min(vlo, axeps[fl][k][i]); vhi = std::max(vhi, axeps[fl][k][i]); }
        const double xlo = xc[0] - exc[0] - 0.1, xhi = xc[NY - 1] + exc[NY - 1] + 0.1;
        TH1F *frA = new TH1F("fr_axeps", "", 100, xlo, xhi);
        frA->SetStats(0);
        ApplyHistStyle(frA, ps3, "#eta^{l}_{CM}", "(A #times #varepsilon)_{MC} = #sigma^{reco}/#sigma^{gen}_{fid}");
        frA->SetMinimum(std::max(0.0, vlo - 0.45 * (vhi - vlo) - 0.05));
        frA->SetMaximum(vhi + 0.55 * (vhi - vlo) + 0.05);
        frA->Draw();

        TGraph *ga[2][2];
        const int col[2] = {kAzure + 2, kRed + 1};   // charge
        // [flavour][charge]: mu filled, e open -- circle/diamond, no triangles
        // (see the kMk* block at the top); diamonds need the larger size to
        // match the circles' visual weight.
        const int mst[2][2] = {{kMkW, kMkWm}, {kMkWOpen, kMkWmOpen}};
        const double msz[2] = {1.4, 2.0};  // [charge]: circle, diamond
        for (int fl = 0; fl < 2; ++fl)
            for (int k = 0; k < 2; ++k)
            {
                ga[fl][k] = new TGraph(NY, xc.data(), axeps[fl][k]);
                ga[fl][k]->SetMarkerStyle(mst[fl][k]); ga[fl][k]->SetMarkerSize(msz[k]);
                ga[fl][k]->SetMarkerColor(col[k]);     ga[fl][k]->SetLineColor(col[k]);
                ga[fl][k]->SetLineWidth(2);            ga[fl][k]->SetLineStyle(fl == 0 ? 1 : 2);
                ga[fl][k]->Draw("PL SAME");
            }

        DrawHeader(ps3, "", "W #rightarrow l #nu", "acceptance #times efficiency");
        TLatex ltA; ltA.SetNDC(); ltA.SetTextFont(ps3.font); ltA.SetTextAlign(13);
        ltA.SetTextSize(ps3.boxTextSize);
        ltA.DrawLatex(ps3.headerX, ps3.headerY - 3.0 * ps3.headerDy, Form("%s fit", discLabel.Data()));

        TLegend *legA = new TLegend(0.18, 0.60, 0.52, 0.80);
        legA->SetBorderSize(0); legA->SetFillStyle(0); legA->SetTextFont(42); legA->SetTextSize(0.030);
        legA->SetNColumns(2);
        legA->AddEntry(ga[0][0], "#mu, W^{+}", "lp");
        legA->AddEntry(ga[1][0], "e, W^{+}", "lp");
        legA->AddEntry(ga[0][1], "#mu, W^{-}", "lp");
        legA->AddEntry(ga[1][1], "e, W^{-}", "lp");
        legA->Draw();

        TLatex txA; txA.SetNDC(); txA.SetTextFont(42); txA.SetTextSize(0.029);
        txA.DrawLatex(0.18, 0.255, "MC correction implicit in #sigma = r #times #sigma^{gen}_{fid}");
        txA.DrawLatex(0.18, 0.215, fidText(disc));

        CMS_lumi(ca, 13, 10);
        ca->RedrawAxis(); ca->Modified(); ca->Update();
        ca->SaveAs((outDir + "/AxEps_vs_eta.png").c_str());
        ca->SaveAs((outDir + "/AxEps_vs_eta.pdf").c_str());
        std::cout << "[OK] Saved " << outDir << "/AxEps_vs_eta.{png,pdf}\n";
    }

    printf("[xsec-diag] cross-check: mu+e of r*sigma^reco = %.2f / %.2f / %.2f nb"
           " (should equal the fork's {Wp,Wm,W}_sum_yield / L)\n",
           sumCount[0], sumCount[1], sumCount[2]);
    std::cout << "[OK] Wrote " << outDir << "/xsec_diag.csv + xsec_diag_bins.csv\n";
    if (fC) { fC->Close(); delete fC; }
}

// =============================================================================
// xsec_fiducial_flav -- the MU-ONLY and E-ONLY simultaneous fits (fork
// run_pO_fits.sh mode flavfit, 2026-09-22) OVERLAID, each extracted exactly as
// xsec_fiducial_comb extracts the grand fit:
//
//     sigma_meas,i = r_i x sigma_gen-fid,i ,  inclusive = Sum_i with the full
//     r-covariance (h_cov_yield), the stat part from h_cov_yield_stat /
//     rErr_stat, drawn as point + stat bar + syst box (MakeSystBoxes)
//
// with the SAME (mu+e pooled) gen fiducial cross sections for every fit, so all
// results are quoted in ONE fiducial volume (bare lepton pT > 25 GeV,
// |eta_lab| < 2.4, and m_T > 40 for leppt_mt40 -- pODisc::GenFid) and
// sigma_e / sigma_mu = r_e / r_mu exactly -- a pure
// data-vs-MC comparison; the ~1-2% bare-lepton FSR difference between the
// flavours sits in each flavour's MC acceptance, where it belongs.
//
//   W_fiducial.{png,pdf}              W+, W-, W AND Z (x10): mu only, e only
//                                     [+ the grand fit, withComb], the r = 1
//                                     gen-fiducial expectation as a dashed
//                                     segment, and (2026-10-05) the CMS pp
//                                     5.02 / 13 TeV fiducial points x A
//                                     (kPPRef, a COARSE comparison: no sqrt(s),
//                                     isospin or fiducial-definition correction).
//                                     Z = r_Z x sigma_gen-fid,Z from the fit's
//                                     summary.csv (total + stat rows). Wide
//                                     canvas, all text in a right-hand panel.
//   W_dsigma_deta_{Wp,Wm,W}.{png,pdf} d(sigma)/d(eta_CM) per charge, mu vs e
//                                     [+ the grand fit, withCombDiff]
//   xsec_flavfit.csv                  per fit: per-bin + inclusive rows (the
//                                     xsec_comb.csv columns + a fit column)
//   xsec_flavfit_ratio.csv            sigma_e / sigma_mu, per bin + inclusive
//
// COMPARING: the two fits' STATISTICAL errors are independent (disjoint event
// samples); the systematics are largely common -- the lumi 3% (most of the
// inclusive syst box) and the theory shapes move both flavours together and
// cancel in sigma_e / sigma_mu. The ratio is therefore quoted with the
// stat-only error, and the per-bin mu-vs-e chi2 is printed both against the
// stat errors and, as a conservative bound, against the total ones. The
// flavour-specific systematics (in-fit ABCD QCD, muSF) are in neither ratio
// error.
// Outputs: ./plots/flavfit/xsec/<disc>/
//   root -l -b -q -e 'gROOT->LoadMacro("xsec_fiducial.C+"); xsec_fiducial_flav("leppt_mt40");'
// =============================================================================
void xsec_fiducial_flav(const char *disc = "leppt_mt40",
                        bool withComb = true,       // the grand fit in the inclusive plot
                        bool withCombDiff = false,  // ... and in the d(sigma)/d(eta) ones
                        const char *genFile = "../skim/rootfile/gen_xsec.root")
{
    gStyle->SetEndErrorSize(4);

    TString dsuf, discLabel;
    if (!pODisc::Spec(disc, dsuf, discLabel)) return;
    const std::string outDir = std::string("./plots/flavfit/xsec/") + disc;

    // ---- the fits, left to right: mu, [comb], e -- each through the SAME
    // ingredient loader and covariance propagation as the grand fit
    // index 3 = Z (2026-10-05): sigma_Z = r_Z x sigma_gen-fid,Z, the fit's one
    // global DY scale (summary.csv rows simfit_lab,r_Z = total error and
    // simfit_lab_stat,r_Z = statistical), the same product xsec_contour.C uses
    struct FitX
    {
        pOFit::Spec spec;
        CombIngredients in;
        double yv[4] = {0, 0, 0, 0}, ev[4] = {0, 0, 0, 0}, evS[4] = {0, 0, 0, 0}; // W+, W-, W, Z
        bool haveZ = false;     // r_Z row found
        bool haveZStat = false; // ... and its stat-only row
    };
    const int klo[3] = {0, 12, 0}, khi[3] = {12, 24, 24};
    // the Z gen fiducial (one bin; skim/gen_xsec.C since 2026-09-15)
    double gZ = 0.0;
    {
        TFile *fG = TFile::Open(genFile, "READ");
        TH1D *hGZ = (fG && !fG->IsZombie()) ? (TH1D *)fG->Get("h_gen_sig_Z") : nullptr;
        if (hGZ) gZ = hGZ->GetBinContent(1);
        else std::cerr << "[WARN] xsec_fiducial_flav: no h_gen_sig_Z in " << genFile
                       << " -> the Z column shows the pp reference only (re-run skim/gen_xsec.C)\n";
        if (fG) { fG->Close(); delete fG; }
    }
    std::vector<FitX> fits;
    const char *tags[3] = {"simfit_mu", "comb", "simfit_ele"};
    for (int t = 0; t < 3; ++t)
    {
        FitX f;
        if (!pOFit::Get(tags[t], f.spec)) continue;
        const bool isComb = (f.spec.flav == "");
        if (isComb && !withComb && !withCombDiff) continue;
        if (!pOFit::Available(dsuf, f.spec))
        {
            std::cerr << "[WARN] xsec_fiducial_flav: no " << pOFit::SummaryFile(dsuf, f.spec, "fitted_yields.root")
                      << " -> '" << f.spec.label << "' not drawn"
                      << (isComb ? "\n" : TString::Format(" (fork: ./run_pO_fits.sh both flavfit --disc %s)\n", disc).Data());
            continue;
        }
        if (!loadCombIngredients(pOFit::SummaryFile(dsuf, f.spec, "W_yields.csv").Data(),
                                 pOFit::SummaryFile(dsuf, f.spec, "fitted_yields.root").Data(),
                                 genFile, disc, f.in))
            continue;
        for (int k = 0; k < 3; ++k)
        {
            sumWithCov(f.in, f.in.G, klo[k], khi[k], f.yv[k], f.ev[k]);
            double dummy;
            sumWithCov(f.in, f.in.G, klo[k], khi[k], dummy, f.evS[k], true);
        }
        if (gZ > 0)
        {
            const TString sumCsv = pOFit::SummaryFile(dsuf, f.spec, "summary.csv");
            double rz = 0, erz = 0, rzS = 0, erzS = 0;
            if (readCombSum(sumCsv.Data(), "simfit_lab", "r_Z", rz, erz))
            {
                f.haveZ = true;
                f.yv[3] = rz * gZ; f.ev[3] = erz * gZ;
                f.haveZStat = f.in.haveStat && readCombSum(sumCsv.Data(), "simfit_lab_stat", "r_Z", rzS, erzS);
                f.evS[3] = f.haveZStat ? erzS * gZ : f.ev[3]; // no stat row -> one (total) bar, no box
            }
        }
        fits.push_back(f);
    }
    const FitX *fm = nullptr, *fe = nullptr;
    for (const FitX &f : fits)
    {
        if (f.spec.flav == "mu") fm = &f;
        if (f.spec.flav == "ele") fe = &f;
    }
    if (!fm && !fe)
    {
        std::cerr << "[ERROR] xsec_fiducial_flav: neither the mu-only nor the e-only fit is available for disc="
                  << disc << " -- fork: ./run_pO_fits.sh both flavfit --disc " << disc << " (+ sync_lxplus.sh download).\n";
        return;
    }
    gSystem->mkdir(outDir.c_str(), kTRUE);

    // the gen fiducial cross sections are the same for every fit (pooled mu+e)
    const double *G = fits[0].in.G;
    double genTot[4] = {0, 0, 0, gZ};
    for (int k = 0; k < 3; ++k)
        for (int i = klo[k]; i < khi[k]; ++i) genTot[k] += G[i];

    // ---- console: every fit, then the e/mu ratios ----------------------------
    const char *cn3[4] = {"Wp", "Wm", "W", "Z"};
    const char *xlab[4] = {"W^{+}", "W^{-}", "W", "Z"};
    for (const FitX &f : fits)
        for (int k = 0; k < 4; ++k)
        {
            if (k == 3 && !f.haveZ) continue;
            const bool st = (k == 3) ? f.haveZStat : f.in.haveStat;
            if (st)
                printf("[xsec-flav] %-10s %-3s measured = %6.2f +/- %4.2f (stat) +/- %4.2f (syst) = +/- %4.2f (total) nb   (gen fid %6.2f nb, eff. r = %.3f)\n",
                       f.spec.tag.Data(), cn3[k], f.yv[k], f.evS[k], systFrom(f.ev[k], f.evS[k], true), f.ev[k],
                       genTot[k], genTot[k] > 0 ? f.yv[k] / genTot[k] : 0.0);
            else
                printf("[xsec-flav] %-10s %-3s measured = %6.2f +/- %4.2f nb (total; no stat split)   (eff. r = %.3f)\n",
                       f.spec.tag.Data(), cn3[k], f.yv[k], f.ev[k], genTot[k] > 0 ? f.yv[k] / genTot[k] : 0.0);
        }
    for (int p = 0; p < 2; ++p)
        printf("[xsec-flav] pp reference %-16s x %.0f: W+ %6.2f  W- %6.2f  W %6.2f  Z %6.2f nb  (%s; drawn as-is, no energy/isospin/fiducial correction)\n",
               kPPRef[p].label, pONorm::kA_O, kPPRef[p].v[0] * pONorm::kA_O, kPPRef[p].v[1] * pONorm::kA_O,
               kPPRef[p].v[2] * pONorm::kA_O, kPPRef[p].v[3] * pONorm::kA_O, kPPRefCite);
    // sigma_e / sigma_mu with the STAT error (independent samples; lumi and the
    // theory shapes cancel) -- the total-error version only when no stat split
    const bool haveRatio = fm && fe;
    const bool ratioStat = haveRatio && fm->in.haveStat && fe->in.haveStat;
    const bool haveRatioZ = haveRatio && fm->haveZ && fe->haveZ;
    double rat[4] = {0, 0, 0, 0}, ratE[4] = {0, 0, 0, 0};
    if (haveRatio)
        for (int k = 0; k < 4; ++k)
        {
            if (k == 3 && !haveRatioZ) continue;
            if (fm->yv[k] <= 0 || fe->yv[k] <= 0) continue;
            rat[k] = fe->yv[k] / fm->yv[k];
            const double em = ratioStat ? fm->evS[k] : fm->ev[k], ee = ratioStat ? fe->evS[k] : fe->ev[k];
            ratE[k] = rat[k] * std::sqrt(std::pow(em / fm->yv[k], 2) + std::pow(ee / fe->yv[k], 2));
            printf("[xsec-flav] sigma_e / sigma_mu  %-3s = %.4f +/- %.4f (%s)  -> %+.1f sigma from 1\n",
                   cn3[k], rat[k], ratE[k], ratioStat ? "stat" : "total", ratE[k] > 0 ? (rat[k] - 1.0) / ratE[k] : 0.0);
        }

    // ---- 1) inclusive W+ / W- / W / Z, with the CMS pp reference --------------
    // (2026-10-05, user request) FOUR categories. In each: the fits mu | [comb]
    // | e on the left, the two pp reference points (kPPRef x A) on the right,
    // the r = 1 gen segment under the fit cluster only. The Z column is drawn
    // x kZScale so all four sit on one linear axis. WIDE canvas with every
    // text in a panel right of the frame (the xsec_contour_WZ_flav layout), so
    // the frame holds nothing but the points and the CMS banner.
    std::vector<const FitX *> incl;
    for (const FitX &f : fits)
        if (f.spec.flav != "" || withComb) incl.push_back(&f);
    {
        const int nF = (int)incl.size();
        const double kZScale = 10.0;
        const double zsc[4] = {1.0, 1.0, 1.0, kZScale};
        const double offF1[1] = {-0.22}, offF2[2] = {-0.32, -0.14}, offF3[3] = {-0.37, -0.22, -0.07};
        const double offPP[2] = {+0.17, +0.33};
        const double genLo = -0.45, genHi = +0.02;
        auto offF = [&](int j) { return (nF == 1) ? offF1[0] : (nF == 2) ? offF2[j] : offF3[j]; };

        PlotStyle ps; ps.logy = false;
        ps.w = 1300; ps.h = 850;                                   // frame ~660 x 705 px, panel ~520 px
        ps.lm = 0.09; ps.rm = 0.40; ps.bm = 0.11; ps.tm = 0.06;
        ps.xTitleSize = ps.yTitleSize = 0.040;
        ps.xLabelSize = 0.050; ps.yLabelSize = 0.034;
        ps.yTitleOffset = 1.05;
        ps.systBoxHalfWidthAbs = 0.05;
        ps.systBoxFillAlpha = 0.30;
        TCanvas *c = new TCanvas("c_xsec_fid_flav", "", ps.w, ps.h);
        ApplyCanvasStyle(c, ps);
        c->cd();

        double ymax = 0.0;
        for (const FitX *f : incl)
            for (int k = 0; k < 4; ++k)
                if (k < 3 || f->haveZ) ymax = std::max(ymax, (f->yv[k] + f->ev[k]) * zsc[k]);
        for (int k = 0; k < 4; ++k) ymax = std::max(ymax, genTot[k] * zsc[k]);
        for (int p = 0; p < 2; ++p)
            for (int k = 0; k < 4; ++k)
                ymax = std::max(ymax, (kPPRef[p].v[k] + ppRefErr(kPPRef[p], k)) * pONorm::kA_O * zsc[k]);
        const double yFrameMax = 1.32 * ymax; // headroom for the in-frame CMS banner only
        TH1F *fr = new TH1F("fr_xsec_flav", "", 4, 0.5, 4.5);
        for (int k = 0; k < 4; ++k)
            fr->GetXaxis()->SetBinLabel(k + 1, k == 3 ? Form("Z (#times%.0f)", kZScale) : xlab[k]);
        fr->SetStats(0);
        ApplyHistStyle(fr, ps, "", "#sigma^{fid} (nb)");
        fr->GetXaxis()->SetNdivisions(4);
        fr->SetMinimum(0.0);
        fr->SetMaximum(yFrameMax);
        fr->Draw();

        // the r = 1 expectation (same for every fit): a dashed segment under
        // the fit cluster of each category
        TLine *lGen = nullptr;
        for (int k = 0; k < 4; ++k)
        {
            if (genTot[k] <= 0) continue;
            lGen = new TLine(k + 1 + genLo, genTot[k] * zsc[k], k + 1 + genHi, genTot[k] * zsc[k]);
            lGen->SetLineStyle(2); lGen->SetLineWidth(3); lGen->SetLineColor(kGreen + 2);
            lGen->Draw();
        }
        // every fit's boxes first, then all points on top
        std::vector<TGraphErrors *> gDraw;
        bool anyBox = false;
        for (int j = 0; j < nF; ++j)
        {
            const FitX *f = incl[j];
            const double o = offF(j);
            std::vector<double> x, y, e, es, ex0;
            for (int k = 0; k < 4; ++k)
            {
                if (k == 3 && !f->haveZ) continue;
                x.push_back(k + 1 + o);
                y.push_back(f->yv[k] * zsc[k]);
                e.push_back(f->ev[k] * zsc[k]);
                // an entry without a stat split gets es = e: a zero-height box,
                // which MakeSystBoxes skips, and the total bar is what is drawn
                const bool st = (k == 3) ? f->haveZStat : f->in.haveStat;
                es.push_back((st ? f->evS[k] : f->ev[k]) * zsc[k]);
                ex0.push_back(0.0);
            }
            const int n = (int)x.size();
            TGraphErrors *gT = new TGraphErrors(n, x.data(), y.data(), ex0.data(), e.data());  // TOTAL error
            TGraphErrors *gS = new TGraphErrors(n, x.data(), y.data(), ex0.data(), es.data()); // STATISTICAL error
            PlotStyle pb = ps;
            pb.systBoxFillColor = f->spec.color; pb.systBoxLineColor = f->spec.color;
            for (TBox *b : MakeSystBoxes(gT, gS, pb)) { b->Draw(); anyBox = true; }
            TGraphErrors *gd = gS; // == gT wherever no split exists (see above)
            gd->SetMarkerStyle(f->spec.marker); gd->SetMarkerSize(f->spec.markerSize);
            gd->SetMarkerColor(f->spec.color); gd->SetLineColor(f->spec.color); gd->SetLineWidth(2);
            gDraw.push_back(gd);
        }
        // the pp reference x A: open markers, one TOTAL bar (stat (+) syst (+) lumi)
        TGraphErrors *gPP[2] = {nullptr, nullptr};
        for (int p = 0; p < 2; ++p)
        {
            double x[4], y[4], e[4], ex0[4] = {0, 0, 0, 0};
            for (int k = 0; k < 4; ++k)
            {
                x[k] = k + 1 + offPP[p];
                y[k] = kPPRef[p].v[k] * pONorm::kA_O * zsc[k];
                e[k] = ppRefErr(kPPRef[p], k) * pONorm::kA_O * zsc[k];
            }
            gPP[p] = new TGraphErrors(4, x, y, ex0, e);
            gPP[p]->SetMarkerStyle(kPPRef[p].marker); gPP[p]->SetMarkerSize(kPPRef[p].markerSize);
            gPP[p]->SetMarkerColor(kPPRef[p].color); gPP[p]->SetLineColor(kPPRef[p].color);
            gPP[p]->SetLineWidth(2);
            gPP[p]->Draw("P SAME");
        }
        for (TGraphErrors *gd : gDraw) gd->Draw("P SAME");

        // ---- the text panel right of the frame -----------------------------------
        const double tx0 = 1.0 - ps.rm + 0.02; // left edge of the panel (NDC)
        TLatex hd; hd.SetNDC(true); hd.SetTextFont(ps.font); hd.SetTextAlign(13);
        hd.SetTextSize(0.034);
        hd.DrawLatex(tx0, 0.935, "W #rightarrow l #nu,  Z #rightarrow l l");
        hd.SetTextSize(0.029);
        hd.DrawLatex(tx0, 0.893, "fiducial cross sections per lepton flavour");
        hd.DrawLatex(tx0, 0.858, Form("#mu-only vs e-only simfits, %s fit", discLabel.Data()));

        const int nLeg = nF + 1 + (anyBox ? 1 : 0) + 2;
        const double legTop = 0.822, rowH = 0.035;
        TLegend *leg = new TLegend(tx0 - 0.01, legTop - rowH * nLeg, 0.995, legTop);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.027);
        for (int j = 0; j < nF; ++j)
            leg->AddEntry(gDraw[j], Form("%s (simfit, pO)", incl[j]->spec.label.Data()), "lep");
        leg->AddEntry(lGen, "#sigma^{gen}_{fid} (POWHEG + EPPS21, r = 1)", "l");
        if (anyBox)
        {
            TBox *key = new TBox(0, 0, 1, 1); // legend glyph only
            key->SetFillColorAlpha(kGray + 1, ps.systBoxFillAlpha); key->SetLineColor(kGray + 3);
            leg->AddEntry(key, "bars: stat.,  boxes: syst.", "f");
        }
        for (int p = 0; p < 2; ++p)
            leg->AddEntry(gPP[p], Form("%s, #times A (A = %.0f)", kPPRef[p].label, pONorm::kA_O), "lep");
        leg->Draw();

        TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextAlign(13);
        double yc = legTop - rowH * nLeg - 0.015;
        if (haveRatio)
        {
            const double dy = 0.033;
            tx.SetTextSize(0.027);
            tx.DrawLatex(tx0, yc, Form("#sigma_{e} / #sigma_{#mu} (%s unc.):", ratioStat ? "stat." : "total"));
            yc -= dy;
            for (int k = 0; k < 4; ++k)
            {
                if (k == 3 && !haveRatioZ) continue;
                tx.DrawLatex(tx0 + 0.02, yc, Form("%s:  %.3f #pm %.3f", xlab[k], rat[k], ratE[k]));
                yc -= dy;
            }
            // mu vs e chi2 over W+ and W- (W is their sum, not independent):
            // (sigma_mu - sigma_e)^2 / (delta_mu^2 + delta_e^2), stat and total
            double c2s = 0, c2t = 0;
            for (int k = 0; k < 2; ++k)
            {
                const double d2 = std::pow(fm->yv[k] - fe->yv[k], 2);
                const double vs = fm->evS[k] * fm->evS[k] + fe->evS[k] * fe->evS[k];
                const double vt = fm->ev[k] * fm->ev[k] + fe->ev[k] * fe->ev[k];
                if (vs > 0) c2s += d2 / vs;
                if (vt > 0) c2t += d2 / vt;
            }
            printf("[mu-vs-e] inclusive W+,W- chi2/ndf = %.2f/2 (stat only, p = %.3f)   %.2f/2 (total, p = %.3f)\n",
                   c2s, TMath::Prob(c2s, 2), c2t, TMath::Prob(c2t, 2));
            yc -= 0.010;
            tx.DrawLatex(tx0, yc, "#mu vs e (W^{+}, W^{-}):");
            yc -= dy;
            tx.DrawLatex(tx0 + 0.02, yc, "#chi^{2} = #Sigma(#sigma_{#mu}#minus#sigma_{e})^{2}/(#delta_{#mu}^{2}+#delta_{e}^{2})");
            yc -= dy;
            tx.DrawLatex(tx0 + 0.02, yc, ratioStat ? Form("#chi^{2}/ndf = %.2f/2 (#delta = stat.),  %.2f/2 (total)", c2s, c2t)
                                                 : Form("#chi^{2}/ndf = %.2f/2 (#delta = total)", c2t));
            yc -= dy;
        }
        // the pp reference, then the pO definitions
        const double dy2 = 0.030;
        tx.SetTextSize(0.024);
        yc -= 0.010;
        tx.SetTextColor(kGray + 3);
        tx.DrawLatex(tx0, yc, Form("pp reference: %s, pp fiducial #times A", kPPRefCite)); yc -= dy2;
        tx.DrawLatex(tx0 + 0.02, yc, "W: p_{T}^{l} > 25 GeV, |#eta| < 2.4, m_{T} > 40 GeV"); yc -= dy2;
        tx.DrawLatex(tx0 + 0.02, yc, "Z: p_{T}^{l} > 25 GeV, |#eta| < 2.4, m_{ll} 60#minus120 GeV"); yc -= dy2;
        tx.DrawLatex(tx0 + 0.02, yc, "no #sqrt{s}, isospin or nuclear correction (coarse)"); yc -= dy2;
        yc -= 0.010;
        tx.SetTextColor(kBlack);
        {
            TString gtag, glabel;
            pODisc::GenFid(disc, gtag, glabel);
            tx.DrawLatex(tx0, yc, "pO fiducial (bare leptons):"); yc -= dy2;
            tx.DrawLatex(tx0 + 0.02, yc, "W: lepton " + glabel); yc -= dy2;
        }
        tx.DrawLatex(tx0 + 0.02, yc, Form("Z: p_{T}^{l} > 15/10 GeV, m_{ll} 60#minus120 GeV, drawn #times%.0f", kZScale)); yc -= dy2;
        tx.DrawLatex(tx0, yc, "#sigma_{i} = r_{i} #times #sigma^{gen}_{fid,i}, summed with the r-covariance"); yc -= dy2;
        tx.DrawLatex(tx0, yc, "syst. = profiled fit nuisances (lumi 3% common)"); yc -= dy2;

        CMS_lumi(c, 13, 10);
        c->RedrawAxis(); c->Modified(); c->Update();
        c->SaveAs((outDir + "/W_fiducial.png").c_str());
        c->SaveAs((outDir + "/W_fiducial.pdf").c_str());
        std::cout << "[OK] Saved " << outDir << "/W_fiducial.{png,pdf}  (disc=" << disc << ")\n";
    }

    // ---- 2) d(sigma)/d(eta_CM) per charge --------------------------------------
    // eta_CM bin positions from the (fiducial) charge-asymmetry graph of any fit
    // (the same lab edges shifted to the CM frame in every one of them)
    std::vector<double> xc(12, 0.0), exc(12, 0.0);
    bool haveEta = false;
    for (const FitX &f : fits)
    {
        const TString sC = pODisc::GraphFile("charge_asym", f.spec.tag.Data(), disc);
        if (gSystem->AccessPathName(sC)) continue;
        TFile *fC = TFile::Open(sC, "READ");
        TGraphErrors *gA = (fC && !fC->IsZombie()) ? (TGraphErrors *)fC->Get("g_chargeAsym") : nullptr;
        if (gA && gA->GetN() == 12)
        {
            for (int i = 0; i < 12; ++i) { xc[i] = gA->GetPointX(i); exc[i] = gA->GetErrorX(i); }
            haveEta = true;
        }
        if (fC) { fC->Close(); delete fC; }
        if (haveEta) break;
    }
    // per-bin value of fit f, charge k (0 W+, 1 W-, 2 W = both, with their
    // covariance) through the same sumWithCov as the inclusive numbers
    auto binVal = [&](const FitX &f, int k, int i, double &v, double &e, double &es) {
        double w[24];
        for (int t = 0; t < 24; ++t) w[t] = 0.0;
        if (k != 1) w[i] = f.in.G[i];
        if (k != 0) w[12 + i] = f.in.G[12 + i];
        sumWithCov(f.in, w, 0, 24, v, e);
        double dummy;
        sumWithCov(f.in, w, 0, 24, dummy, es, true);
    };
    std::ofstream rcsv((outDir + "/xsec_flavfit_ratio.csv").c_str());
    rcsv << "charge,ybin,etaCM_center,ratio_e_over_mu,err_stat,err_total\n";
    if (!haveEta)
        std::cerr << "[WARN] xsec_fiducial_flav: no 12-point g_chargeAsym in any charge_asym_fid_<fit>_" << disc
                  << ".root -> d(sigma)/d(eta) plots skipped (analysis/run_observables.sh makes them first)\n";
    else
    {
        std::vector<const FitX *> diff;
        for (const FitX &f : fits)
            if (f.spec.flav != "" || withCombDiff) diff.push_back(&f);
        const int nF = (int)diff.size();
        const double sh2[2] = {-0.35, 0.35}, sh3[3] = {-0.60, 0.0, 0.60};
        for (int k = 0; k < 3; ++k)
        {
            std::vector<std::vector<double>> v(nF, std::vector<double>(12)), e(nF, std::vector<double>(12)),
                es(nF, std::vector<double>(12));
            std::vector<double> gen(12), zero(12, 0.0);
            double ymax = 0.0;
            for (int i = 0; i < 12; ++i)
            {
                double dEta = 2.0 * exc[i];
                if (dEta <= 0) dEta = 0.4;
                gen[i] = ((k != 1 ? G[i] : 0.0) + (k != 0 ? G[12 + i] : 0.0)) / dEta;
                ymax = std::max(ymax, gen[i]);
                for (int j = 0; j < nF; ++j)
                {
                    binVal(*diff[j], k, i, v[j][i], e[j][i], es[j][i]);
                    v[j][i] /= dEta; e[j][i] /= dEta; es[j][i] /= dEta;
                    ymax = std::max(ymax, v[j][i] + e[j][i]);
                }
            }

            PlotStyle ps; ps.logy = false;
            ps.systBoxWidthFrac = (nF >= 3) ? 0.17 : 0.20;
            ps.systBoxFillAlpha = 0.30;
            TCanvas *c = new TCanvas(Form("c_dsig_flav_%s", cn3[k]), "", ps.w, ps.h);
            ApplyCanvasStyle(c, ps);
            c->cd();
            const double xlo = xc[0] - exc[0] - 0.1, xhi = xc[11] + exc[11] + 0.1;
            TH1F *fr = new TH1F(Form("fr_dsig_flav_%s", cn3[k]), "", 100, xlo, xhi);
            fr->SetStats(0);
            ApplyHistStyle(fr, ps, "#eta^{l}_{CM}", "d#sigma^{fid}_{W#rightarrowl#nu}/d#eta (nb)");
            fr->SetMinimum(0.0);
            // annotations in the band above NDC kAnnoY (as xsec_fiducial_comb)
            const double kAnnoY = 0.55;
            const double fAnno = (kAnnoY - ps.bm) / (1.0 - ps.bm - ps.tm);
            double yMaxFrame = 1.55 * ymax;
            if (fAnno > 0.05) yMaxFrame = std::max(yMaxFrame, 1.06 * ymax / fAnno);
            fr->SetMaximum(yMaxFrame);
            fr->Draw();

            TGraph *gGen = new TGraph(12, xc.data(), gen.data());
            gGen->SetLineStyle(2); gGen->SetLineWidth(3); gGen->SetLineColor(kGreen + 2);
            gGen->Draw("L SAME");

            std::vector<TGraphErrors *> gDraw;
            bool anyBox = false;
            for (int j = 0; j < nF; ++j)
            {
                const FitX *f = diff[j];
                const double s = (nF == 1) ? 0.0 : (nF == 2) ? sh2[j] : sh3[j];
                std::vector<double> xs(12);
                for (int i = 0; i < 12; ++i) xs[i] = xc[i] + s * exc[i];
                // box sources keep the bin half-widths (MakeSystBoxes sizes the
                // boxes from them); the drawn graph has no horizontal bars
                TGraphErrors *gT = new TGraphErrors(12, xs.data(), v[j].data(), exc.data(), e[j].data());
                TGraphErrors *gS = new TGraphErrors(12, xs.data(), v[j].data(), exc.data(), es[j].data());
                PlotStyle pb = ps;
                pb.systBoxFillColor = f->spec.color; pb.systBoxLineColor = f->spec.color;
                if (f->in.haveStat)
                    for (TBox *b : MakeSystBoxes(gT, gS, pb)) { b->Draw(); anyBox = true; }
                TGraphErrors *gd = new TGraphErrors(12, xs.data(), v[j].data(), zero.data(),
                                                    f->in.haveStat ? es[j].data() : e[j].data());
                gd->SetMarkerStyle(f->spec.marker); gd->SetMarkerSize(f->spec.markerSize);
                gd->SetMarkerColor(f->spec.color); gd->SetLineColor(f->spec.color); gd->SetLineWidth(2);
                gDraw.push_back(gd);
            }
            for (TGraphErrors *gd : gDraw) gd->Draw("P SAME");

            const char *proc = (k == 0) ? "W^{+} #rightarrow l^{+} #nu" : (k == 1) ? "W^{-} #rightarrow l^{-} #bar{#nu}"
                                                                               : "W #rightarrow l #nu";
            DrawHeader(ps, "", proc, "fiducial cross section");
            TLatex ltag; ltag.SetNDC(); ltag.SetTextFont(ps.font); ltag.SetTextAlign(13);
            ltag.SetTextSize(ps.boxTextSize);
            ltag.DrawLatex(ps.headerX, ps.headerY - 3.0 * ps.headerDy, Form("%s fit", discLabel.Data()));

            TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.029);
            // right column under the header; the band just above kAnnoY holds
            // the full-width mu-vs-e chi2 lines
            tx.DrawLatex(0.55, kAnnoY + 0.120, "#sigma_{i} = r_{i} #times #sigma^{gen}_{fid,i} (p_{T}^{l} > 25 GeV)");
            tx.DrawLatex(0.55, kAnnoY + 0.077, "bars: stat.,  boxes: syst. (profiled)");

            const int nLeg = nF + 1;
            TLegend *leg = new TLegend(0.18, 0.79 - 0.045 * nLeg, 0.50, 0.79);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.031);
            for (int j = 0; j < nF; ++j)
                leg->AddEntry(gDraw[j], Form("%s (simfit)", diff[j]->spec.label.Data()), "lep");
            leg->AddEntry(gGen, "#sigma^{gen}_{fid} (POWHEG, r = 1)", "l");
            leg->Draw();
            (void)anyBox;

            // mu vs e, bin by bin: chi2 = sum (x_mu - x_e)^2 / (sigma_mu^2 + sigma_e^2)
            // against the stat and the total errors (the two fits are
            // independent samples; the total version double-counts the common
            // lumi 3%, so it is the conservative bound), printed on the plot
            // and the console, plus the per-bin e/mu ratio for the CSV
            if (fm && fe)
            {
                double vm, emT, emS, ve, eeT, eeS, c2s = 0, c2t = 0;
                int n = 0;
                for (int i = 0; i < 12; ++i)
                {
                    binVal(*fm, k, i, vm, emT, emS);
                    binVal(*fe, k, i, ve, eeT, eeS);
                    const double vt = emT * emT + eeT * eeT, vs = emS * emS + eeS * eeS;
                    if (vt <= 0 || vm <= 0 || ve <= 0) continue;
                    c2t += (vm - ve) * (vm - ve) / vt;
                    if (vs > 0) c2s += (vm - ve) * (vm - ve) / vs;
                    ++n;
                    const double rr = ve / vm;
                    rcsv << Form("%s,%d,%.4f,%.5f,%.5f,%.5f\n", cn3[k], i, xc[i], rr,
                                 rr * std::sqrt(std::pow(emS / vm, 2) + std::pow(eeS / ve, 2)),
                                 rr * std::sqrt(std::pow(emT / vm, 2) + std::pow(eeT / ve, 2)));
                }
                if (n > 0)
                {
                    printf("[mu-vs-e] d(sigma)/d(eta) %-3s chi2/ndf = %5.1f/%d (stat only, p = %.3f)   %5.1f/%d (total, p = %.3f)\n",
                           cn3[k], c2s, n, TMath::Prob(c2s, n), c2t, n, TMath::Prob(c2t, n));
                    // under the legend, left of the right-hand annotation band
                    TLatex tc; tc.SetNDC(); tc.SetTextFont(42); tc.SetTextSize(0.027);
                    tc.DrawLatex(0.19, kAnnoY + 0.045, "#mu vs e: #chi^{2} = #Sigma(#sigma_{#mu}#minus#sigma_{e})^{2}/(#delta_{#mu}^{2}+#delta_{e}^{2})");
                    tc.DrawLatex(0.19, kAnnoY + 0.005, (fm->in.haveStat && fe->in.haveStat)
                        ? Form("#chi^{2}/ndf = %.1f/%d (stat.),  %.1f/%d (total)", c2s, n, c2t, n)
                        : Form("#chi^{2}/ndf = %.1f/%d (total)", c2t, n));
                }
            }

            CMS_lumi(c, 13, 10);
            c->RedrawAxis(); c->Modified(); c->Update();
            c->SaveAs(Form("%s/W_dsigma_deta_%s.png", outDir.c_str(), cn3[k]));
            c->SaveAs(Form("%s/W_dsigma_deta_%s.pdf", outDir.c_str(), cn3[k]));
        }
        std::cout << "[OK] Saved " << outDir << "/W_dsigma_deta_{Wp,Wm,W}.{png,pdf}\n";
    }
    if (haveRatio)
        for (int k = 0; k < 4; ++k)
        {
            if (k == 3 && !haveRatioZ) continue;
            rcsv << Form("%s,incl,,%.5f,%.5f,%.5f\n", cn3[k], rat[k], ratioStat ? ratE[k] : -1.0,
                         rat[k] * std::sqrt(std::pow(fm->ev[k] / fm->yv[k], 2) + std::pow(fe->ev[k] / fe->yv[k], 2)));
        }
    rcsv.close();

    // ---- 3) the sigma table, every fit (the xsec_comb.csv columns + fit) ------
    {
        std::ofstream out((outDir + "/xsec_flavfit.csv").c_str());
        out << "fit,charge,ybin,etaCM_center,etaCM_halfwidth,r,rErr,sigma_gen_fid_nb,sigma_meas_nb,sigma_meas_err_nb,"
               "sigma_meas_stat_nb,sigma_meas_syst_nb\n";
        const char *cn2[2] = {"Wp", "Wm"};
        for (const FitX &f : fits)
        {
            for (int q = 0; q < 2; ++q)
                for (int i = 0; i < 12; ++i)
                {
                    const int k = 12 * q + i;
                    const double st = f.in.haveStat ? f.in.ERs[k] * G[k] : -1.0;
                    out << Form("%s,%s,%d,%.4f,%.4f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n",
                                f.spec.tag.Data(), cn2[q], i, haveEta ? xc[i] : -1.0, haveEta ? exc[i] : -1.0,
                                f.in.R[k], f.in.ER[k], G[k], f.in.R[k] * G[k], f.in.ER[k] * G[k],
                                st, f.in.haveStat ? systFrom(f.in.ER[k] * G[k], st, true) : -1.0);
                }
            for (int k = 0; k < 4; ++k) // r column = effective (gen-weighted mean) r; k = 3: Z = r_Z x sigma_gen,Z
            {
                if (k == 3 && !f.haveZ) continue;
                const bool st = (k == 3) ? f.haveZStat : f.in.haveStat;
                out << Form("%s,%s,incl,,,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", f.spec.tag.Data(), cn3[k],
                            genTot[k] > 0 ? f.yv[k] / genTot[k] : 0.0, genTot[k] > 0 ? f.ev[k] / genTot[k] : 0.0,
                            genTot[k], f.yv[k], f.ev[k], st ? f.evS[k] : -1.0,
                            st ? systFrom(f.ev[k], f.evS[k], true) : -1.0);
            }
        }
        // the pp reference as drawn (x A), one row per energy and process
        for (int p = 0; p < 2; ++p)
            for (int k = 0; k < 4; ++k)
                out << Form("pp_%s,%s,incl,,,,,,%.5f,%.5f,,\n", p == 0 ? "5.02TeV" : "13TeV", cn3[k],
                            kPPRef[p].v[k] * pONorm::kA_O, ppRefErr(kPPRef[p], k) * pONorm::kA_O);
    }
    std::cout << "[OK] Wrote " << outDir << "/xsec_flavfit.csv + xsec_flavfit_ratio.csv\n";
}

// =============================================================================
// xsec_fiducial(disc) -- every cross-section view of one discriminant: the
// grand fit (xsec_fiducial_comb + its per-flavour diagnostic
// xsec_fiducial_diag) and the mu-vs-e overlay (xsec_fiducial_flav).
//   root -l -b -q 'xsec_fiducial.C+("leppt_mt40")'
// =============================================================================
void xsec_fiducial(const char *disc = "leppt_mt40")
{
    xsec_fiducial_comb(disc);
    xsec_fiducial_diag(disc);
    xsec_fiducial_flav(disc);
}
