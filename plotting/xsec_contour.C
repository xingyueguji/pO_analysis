// =============================================================================
// xsec_contour.C -- the (sigma_W, sigma_Z) covariance ellipse from the grand
// simultaneous fit (2026-09-15).
//
// WHY THIS EXISTS
// ---------------
// Both inclusive cross sections are LINEAR functions of the simfit POIs:
//
//     sigma_W = Sum_i r_i * sigma_gen-fid,i        (i = 24 charge x rapidity bins)
//     sigma_Z = r_Z * sigma_gen-fid,Z              (ONE global DY scale)
//
// so their joint uncertainty is J V J^T with V the POI covariance and J the gen
// fiducial cross sections. The 24x24 h_cov_yield written for xsec_fiducial.C
// covers the W POIs only and has no r_Z row, so the fork's extractor also
// writes the full **25x25** h_cov_poi[_FB][_stat] in PARAMETER space, order
// [r_Wp_y0..11, r_Wm_y0..11, r_Z] -- that matrix is what this macro turns into
// the 2x2 (sigma_W, sigma_Z) covariance and its confidence ellipses.
//
// WHAT THE ELLIPSE IS (and is not)
// --------------------------------
// It is the GAUSSIAN (Hesse) confidence region: the exact contour if the
// profile likelihood is parabolic in these two directions. Here that is a good
// approximation -- the dominant systematic is the 3% lumi log-normal, whose
// exact profiled interval differs from the linearized one by ~0.1 nb on
// sigma_W (~2% of the error), and every r_i sits >10 sigma from its r >= 0
// boundary. The EXACT profiled contour needs sigma_W promoted to a real POI
// (reparametrized workspace) and `combine -M MultiDimFit --algo grid -P
// <sigmaW> -P r_Z`; this macro is both the fast answer and the cross-check
// that the reparametrized scan must reproduce.
//
// NOTE the correlation is dominated by LUMI (3%, fully correlated between the
// two, positive) fighting the DY-background anticorrelation in the W channels.
// The ellipse is therefore elongated along the diagonal, i.e. the sigma_W /
// sigma_Z RATIO is much better measured than either cross section -- that
// short axis is where nPDF sets separate, and is the point of the 2D plot.
//
// THEORY POINT: the open marker at (Sum_i sigma_gen,i, sigma_gen,Z) is r = 1,
// i.e. POWHEG with the generation nPDF (EPPS21nlo_CT18Anlo_O16 central, see
// skim/lhe_index.h). Its nPDF uncertainty ellipse needs the gen-level member
// twins in skim/gen_xsec.C (the deferred LHE TODO) and is NOT drawn yet; other
// nPDF sets need external absolute fiducial predictions.
//
// Outputs: plots/comb/xsec/<disc>/{xsec_contour_WZ_<binning>}.{png,pdf}
//          + xsec_contour_WZ_<binning>.csv (the point, the 2x2, and the 3x3
//            over (W+, W-, Z) for whoever wants the charge-split version)
//
//   root -l -q -e 'gROOT->LoadMacro("xsec_contour.C+"); xsec_contour_WZ("leppt_mt40");'
// =============================================================================
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLeaf.h"
#include "TList.h"
#include "TObjArray.h"
#include "TROOT.h"
#include "TTree.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TString.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "plotting_helper.C"  // PlotStyle, ApplyCanvasStyle, DrawHeader, DrawInfoBox
#include "disc_variants.h"    // pODisc::Spec
#include "../skim/mc_norm.h"  // pONorm::kLumi_invnb

namespace {

// ---- comb_W_yields.csv reader (self-contained; the file also carries r_Z) ----
// layout: region,charge,binning,ybin,r,rErr,signal_prefit,...,r_Z(13),r_ZErr(14),...
struct FitR
{
    double r[24] = {0};   // [Wp_y0..11, Wm_y0..11]
    double rZ = 0, rZErr = -1;
    bool ok = false;
};

FitR readFitR(const char *csv, const char *binning)
{
    FitR out;
    std::ifstream in(csv);
    if (!in) { std::cerr << "[ERROR] cannot open CSV: " << csv << "\n"; return out; }
    std::string line;
    std::getline(in, line); // header
    int found = 0;
    while (std::getline(in, line))
    {
        std::stringstream ss(line);
        std::string f;
        std::vector<std::string> c;
        while (std::getline(ss, f, ',')) c.push_back(f);
        if (c.size() < 15) continue;
        if (c[2] != binning) continue;
        const int iy = std::atoi(c[3].c_str());
        if (iy < 0 || iy >= 12) continue;
        const int off = (c[1] == "Wp") ? 0 : (c[1] == "Wm") ? 12 : -1;
        if (off < 0) continue;
        out.r[off + iy] = std::atof(c[4].c_str());
        out.rZ    = std::atof(c[13].c_str());
        out.rZErr = std::atof(c[14].c_str());
        ++found;
    }
    out.ok = (found == 24);
    if (!out.ok)
        std::cerr << "[ERROR] found " << found << "/24 '" << binning << "' rows in " << csv << "\n";
    return out;
}

// ---- the 25x25 POI covariance -----------------------------------------------
// Read by AXIS LABEL, never by position: a future model change that inserts a
// POI would otherwise silently shuffle the matrix into the wrong J.
bool readPoiCov(TFile *f, const char *name, double C[25][25])
{
    TH2D *h = (TH2D *)f->Get(name);
    if (!h) return false;
    if (h->GetNbinsX() != 25 || h->GetNbinsY() != 25)
    {
        std::cerr << "[ERROR] " << name << " is " << h->GetNbinsX() << "x" << h->GetNbinsY()
                  << ", expected 25x25\n";
        return false;
    }
    // expected label order
    std::vector<TString> want;
    for (int i = 0; i < 12; ++i) want.push_back(TString::Format("r_Wp_y%d", i));
    for (int i = 0; i < 12; ++i) want.push_back(TString::Format("r_Wm_y%d", i));
    want.push_back("r_Z");
    for (int a = 0; a < 25; ++a)
    {
        const TString got = h->GetXaxis()->GetBinLabel(a + 1);
        if (got != want[a])
        {
            std::cerr << "[ERROR] " << name << " bin " << a + 1 << " is labelled '" << got
                      << "', expected '" << want[a] << "' -- refusing to guess the order\n";
            return false;
        }
    }
    for (int a = 0; a < 25; ++a)
        for (int b = 0; b < 25; ++b) C[a][b] = h->GetBinContent(a + 1, b + 1);
    return true;
}

// ---- the EPPS21 member scatter of the PREDICTION (2026-09-15c) --------------
// skim/gen_xsec.C stores the gen fiducial cross sections for all 107 members of
// EPPS21nlo_CT18Anlo_O16 (h_gen_sig_{Wp,Wm}[_FB]_epps21, x = rapidity bin,
// y = member; h_gen_sig_Z_epps21). Member 0 IS the nominal (verified there to
// 0), so points 1..106 are the variations and the r = 1 marker is member 0.
// Scattering them shows the nPDF coverage of the prediction directly -- no
// Hessian combination, just where the members actually land.
// Returns the number of variation points, 0 when the twins are absent.
int epps21Scatter(const char *genFile, bool isFB, TGraph *&g,
                  double &wLo, double &wHi, double &zLo, double &zHi)
{
    g = nullptr;
    TFile *f = TFile::Open(genFile, "READ");
    if (!f || f->IsZombie()) { if (f) { f->Close(); delete f; } return 0; }
    TH2D *hp = (TH2D *)f->Get(isFB ? "h_gen_sig_Wp_FB_epps21" : "h_gen_sig_Wp_epps21");
    TH2D *hm = (TH2D *)f->Get(isFB ? "h_gen_sig_Wm_FB_epps21" : "h_gen_sig_Wm_epps21");
    TH1D *hz = (TH1D *)f->Get("h_gen_sig_Z_epps21");
    if (!hp || !hm || !hz)
    {
        std::cout << "[INFO] no EPPS21 member twins in " << genFile
                  << " -> no member scatter (re-run skim/gen_xsec.C, which gained them 2026-09-15c)\n";
        f->Close(); delete f; return 0;
    }
    const int nm = hz->GetNbinsX();
    if (hp->GetNbinsY() != nm || hm->GetNbinsY() != nm)
    {
        std::cerr << "[WARN] EPPS21 member counts disagree (W " << hp->GetNbinsY()
                  << "/" << hm->GetNbinsY() << " vs Z " << nm << ") -> scatter skipped\n";
        f->Close(); delete f; return 0;
    }
    g = new TGraph();
    g->SetName("g_epps21_members");
    wLo = zLo = 1e30; wHi = zHi = -1e30;
    int n = 0;
    for (int m = 2; m <= nm; ++m)   // ROOT bin 1 = member 0 = the nominal
    {
        double sw = 0;
        for (int i = 1; i <= hp->GetNbinsX(); ++i)
            sw += hp->GetBinContent(i, m) + hm->GetBinContent(i, m);
        const double sz = hz->GetBinContent(m);
        if (sw <= 0 || sz <= 0) continue;
        g->SetPoint(n++, sw, sz);
        wLo = std::min(wLo, sw); wHi = std::max(wHi, sw);
        zLo = std::min(zLo, sz); zHi = std::max(zHi, sz);
    }
    f->Close(); delete f;
    if (n == 0) { delete g; g = nullptr; }
    return n;
}

// ---- ellipse of the Delta(chi2) = s contour of a 2x2 covariance -------------
// {x : (x-mu)^T V^-1 (x-mu) = s}, drawn through the eigen-decomposition
// x = mu + sqrt(s) * (sqrt(l1) v1 cos t + sqrt(l2) v2 sin t).
TGraph *ellipse(double mx, double my, double vxx, double vyy, double vxy,
                double s, int npt = 400)
{
    const double tr = vxx + vyy, det = vxx * vyy - vxy * vxy;
    const double disc = std::sqrt(std::max(0.0, 0.25 * tr * tr - det));
    const double l1 = 0.5 * tr + disc, l2 = std::max(0.0, 0.5 * tr - disc);
    // eigenvector of l1
    double ex = 1, ey = 0;
    if (std::fabs(vxy) > 1e-12) { ex = l1 - vyy; ey = vxy; }
    else if (vyy > vxx)         { ex = 0; ey = 1; }
    const double norm = std::sqrt(ex * ex + ey * ey);
    ex /= norm; ey /= norm;
    const double a = std::sqrt(s * l1), b = std::sqrt(s * std::max(0.0, l2));
    TGraph *g = new TGraph(npt + 1);
    for (int i = 0; i <= npt; ++i)
    {
        const double t = 2.0 * TMath::Pi() * i / npt;
        const double u = a * std::cos(t), v = b * std::sin(t);
        g->SetPoint(i, mx + u * ex - v * ey, my + u * ey + v * ex);
    }
    return g;
}

// ---- the PROFILED contour from a `MultiDimFit --algo grid` scan -------------
// Reads the `limit` tree of higgsCombine_contour_<B>.MultiDimFit.mH120.root
// (run_pO_fits.sh --contour: POIs sigmaW and r_Z, everything else profiled),
// interpolates 2*deltaNLL over (sigma_W, sigma_Z = r_Z * gZ) and returns the
// TGraphs of the Delta(chi2) = 2.296 / 5.991 contours. This is the EXACT
// region; the covariance ellipses drawn alongside are its Gaussian
// approximation, and the two lying on top of each other is the statement that
// the Hesse error on the inclusive cross section is trustworthy.
// Locate the scan file without hardcoding Combine's mass tag. `combine -M
// MultiDimFit -n _contour_<B>` writes higgsCombine_contour_<B>.MultiDimFit.mH<M>.root
// with M = 120 by default -- but that default is not ours to rely on, so try
// the usual name and otherwise glob the directory for the one file with the
// right prefix. Returns "" when the contour pass has not been run.
TString resolveScanFile(const TString &dir, const TString &bn, const char *stem = "contour")
{
    const TString pre = TString::Format("higgsCombine_%s_%s.MultiDimFit.mH", stem, bn.Data());
    const TString usual = dir + "/" + pre + "120.root";
    if (!gSystem->AccessPathName(usual)) return usual;     // false = exists
    TSystemDirectory sd("scan", dir);
    TList *files = sd.GetListOfFiles();
    if (!files) return "";
    TString hit;
    int nhit = 0;
    TIter next(files);
    while (TObject *o = next())
    {
        const TString fn = o->GetName();
        if (!fn.BeginsWith(pre) || !fn.EndsWith(".root")) continue;
        ++nhit;
        if (fn > hit) hit = fn;                            // deterministic pick
    }
    delete files;
    if (nhit > 1)
        std::cerr << "[WARN] " << nhit << " scan files match " << pre << "*.root in " << dir
                  << " -> using " << hit << " (delete the stale ones)\n";
    return hit.IsNull() ? TString("") : dir + "/" + hit;
}

bool profiledContours(const char *scanFile, double gZ,
                      std::vector<TGraph *> &c68, std::vector<TGraph *> &c95,
                      double &bestW, double &bestZ)
{
    // check first: TFile::Open on a missing path prints a scary ROOT error, and
    // "no scan yet" is the NORMAL state until someone runs --contour
    if (gSystem->AccessPathName(scanFile)) return false;   // true = missing
    TFile *f = TFile::Open(scanFile, "READ");
    if (!f || f->IsZombie()) { if (f) { f->Close(); delete f; } return false; }
    TTree *t = (TTree *)f->Get("limit");
    if (!t) { std::cerr << "[WARN] no 'limit' tree in " << scanFile << "\n"; f->Close(); delete f; return false; }
    // Read through TLeaf, NOT SetBranchAddress: combine stores sigmaW, r_Z AND
    // deltaNLL as Float_t, and binding a double* to a float branch makes
    // SetBranchAddress return a mismatch code -- which is exactly how this
    // silently reported "no scan" on the first REAL file (2026-09-15c). The
    // synthetic test fixture could not catch it: it was written with the types
    // this reader assumed, not the ones combine actually writes.
    // TLeaf::GetValue() returns Double_t whatever the storage type is, so this
    // also survives combine ever switching to /D.
    TLeaf *lw = t->GetLeaf("sigmaW"), *lz = t->GetLeaf("r_Z"), *ld = t->GetLeaf("deltaNLL");
    if (!lw || !lz || !ld)
    {
        std::cerr << "[WARN] " << scanFile << " has no sigmaW/r_Z/deltaNLL leaves"
                  << " -- was it made with --contour?\n";
        f->Close(); delete f; return false;
    }
    // NB MultiDimFit writes quantileExpected = -1 on EVERY row (it is a fit, not
    // a limit), so it cannot mark the leading best-fit row -- and need not: that
    // row is simply the point at the minimum.

    TGraph2D *g2 = new TGraph2D();
    g2->SetName("g_scan2d");
    int n = 0, nNeg = 0, nBad = 0;
    bestW = bestZ = 0;
    double bestNll = 1e300;
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
        t->GetEntry(i);
        const double sw = lw->GetValue(), rz = lz->GetValue();
        double dnll = ld->GetValue();
        if (!std::isfinite(dnll)) { ++nBad; continue; }   // failed grid points
        // a slightly negative deltaNLL means that grid point found a marginally
        // better minimum than the reported best fit: keep it, clamped
        if (dnll < 0) { ++nNeg; dnll = 0; }
        g2->SetPoint(n++, sw, rz * gZ, 2.0 * dnll);
        if (dnll < bestNll) { bestNll = dnll; bestW = sw; bestZ = rz * gZ; }
    }
    f->Close(); delete f;
    if (nBad) std::cerr << "[WARN] " << nBad << " scan points with non-finite deltaNLL (dropped)\n";
    if (nNeg) std::cerr << "[WARN] " << nNeg << " scan points with deltaNLL < 0 (clamped to 0 --"
                        << " the grid beat the reported minimum; check the fit)\n";
    if (n < 10) { std::cerr << "[WARN] only " << n << " usable scan points in " << scanFile << "\n"; return false; }

    g2->SetNpx(250); g2->SetNpy(250);
    TH2D *h2 = (TH2D *)g2->GetHistogram();
    if (!h2) { std::cerr << "[WARN] could not interpolate the scan grid\n"; return false; }
    double levels[2] = {2.296, 5.991};                        // 68%, 95% for 2 dof
    h2->SetContour(2, levels);
    // CONT LIST publishes the contour TGraphs under gROOT's specials; it needs
    // a real (if throwaway) pad and an Update to actually run.
    TCanvas tmp("c_contour_tmp", "", 10, 10);
    h2->Draw("CONT Z LIST");
    tmp.Update();
    TObjArray *conts = (TObjArray *)gROOT->GetListOfSpecials()->FindObject("contours");
    if (!conts || conts->GetSize() < 2) { std::cerr << "[WARN] contour extraction returned nothing\n"; return false; }
    for (int lev = 0; lev < 2; ++lev)
    {
        TList *l = (TList *)conts->At(lev);
        if (!l) continue;
        for (int k = 0; k < l->GetSize(); ++k)
        {
            TGraph *g = (TGraph *)((TGraph *)l->At(k))->Clone(Form("c%d_%d", lev, k));
            (lev == 0 ? c68 : c95).push_back(g);
        }
    }
    std::cout << Form("  [profiled] %d scan points -> %d/%d contour piece(s) at 68%%/95%%;"
                      " scan minimum at (%.2f, %.3f)\n",
                      n, (int)c68.size(), (int)c95.size(), bestW, bestZ);
    return !c68.empty();
}

} // namespace

// =============================================================================
void xsec_contour_WZ(const char *disc    = "leppt_mt40",
                     const char *binning = "lab",   // lab | fb (fb = the shifted
                                                    // FB edge set, a DIFFERENT
                                                    // lab window -- gen sigmas
                                                    // follow automatically)
                     const char *binsCsv    = nullptr,
                     const char *yieldsRoot = nullptr,
                     const char *genFile    = "../skim/rootfile/gen_xsec.root",
                     const char *scanFile   = nullptr)  // MultiDimFit grid from
                                                        // run_pO_fits.sh --contour;
                                                        // absent -> ellipse only
{
    TString dsuf, discLabel;
    if (!pODisc::Spec(disc, dsuf, discLabel)) return;
    const TString bn = binning;
    if (bn != "lab" && bn != "fb")
    { std::cerr << "[ERROR] binning must be 'lab' or 'fb', got '" << bn << "'\n"; return; }
    const bool isFB = (bn == "fb");

    const TString base = TString::Format(
        "../../HiggsAnalysis-CombinedLimit/test/pO_fit_out%s/simfit/summary", dsuf.Data());
    const TString sBins = binsCsv    ? TString(binsCsv)    : base + "/comb_W_yields.csv";
    const TString sRoot = yieldsRoot ? TString(yieldsRoot) : base + "/comb_fitted_yields.root";
    const std::string outDir = std::string("./plots/comb/xsec/") + disc;
    gSystem->mkdir(outDir.c_str(), kTRUE);

    // ---- 1) gen fiducial cross sections --------------------------------------
    TFile *fG = TFile::Open(genFile, "READ");
    if (!fG || fG->IsZombie())
    { std::cerr << "[ERROR] cannot open " << genFile << " (run skim/gen_xsec.C)\n"; return; }
    TH1D *hGp = (TH1D *)fG->Get(isFB ? "h_gen_sig_Wp_FB" : "h_gen_sig_Wp");
    TH1D *hGm = (TH1D *)fG->Get(isFB ? "h_gen_sig_Wm_FB" : "h_gen_sig_Wm");
    TH1D *hGZ = (TH1D *)fG->Get("h_gen_sig_Z");
    if (!hGp || !hGm || hGp->GetNbinsX() != 12 || hGm->GetNbinsX() != 12)
    { std::cerr << "[ERROR] no usable 12-bin h_gen_sig_{Wp,Wm}" << (isFB ? "_FB" : "") << " in " << genFile << "\n"; return; }
    if (!hGZ)
    {
        std::cerr << "[ERROR] no h_gen_sig_Z in " << genFile
                  << " -- re-run skim/gen_xsec.C (it gained the Z fiducial on 2026-09-15)\n";
        return;
    }
    double G[25] = {0};
    for (int i = 0; i < 12; ++i) { G[i] = hGp->GetBinContent(i + 1); G[12 + i] = hGm->GetBinContent(i + 1); }
    const double gZ = hGZ->GetBinContent(1), gZErr = hGZ->GetBinError(1);
    G[24] = gZ;
    fG->Close(); delete fG;

    // ---- 2) fitted POIs + the 25x25 covariance -------------------------------
    const FitR fit = readFitR(sBins.Data(), bn.Data());
    if (!fit.ok) return;
    TFile *fY = TFile::Open(sRoot.Data(), "READ");
    if (!fY || fY->IsZombie())
    { std::cerr << "[ERROR] cannot open " << sRoot << "\n"; return; }
    const TString cname  = isFB ? "h_cov_poi_FB" : "h_cov_poi";
    double C[25][25] = {{0}}, Cs[25][25] = {{0}};
    if (!readPoiCov(fY, cname, C))
    {
        std::cerr << "[ERROR] no usable " << cname << " in " << sRoot
                  << " -- re-run the fork's extraction (`run_pO_fits.sh ... --extract-only`)\n";
        fY->Close(); delete fY; return;
    }
    const bool haveStat = readPoiCov(fY, cname + "_stat", Cs);
    if (!haveStat)
        std::cerr << "[WARN] no " << cname << "_stat -> only the TOTAL ellipse is drawn\n";
    fY->Close(); delete fY;

    // ---- 3) propagate: sigma_W = Sum_i G_i r_i, sigma_Z = G_Z r_Z ------------
    // J rows: W+ (0..11), W- (12..23), Z (24). The 3x3 is written to the CSV;
    // the plot uses the (W = W+ + W-, Z) 2x2.
    double J[3][25] = {{0}};
    for (int i = 0; i < 12; ++i) { J[0][i] = G[i]; J[1][12 + i] = G[12 + i]; }
    J[2][24] = G[24];
    double mu3[3] = {0, 0, 0};
    for (int k = 0; k < 3; ++k)
        for (int a = 0; a < 25; ++a) mu3[k] += J[k][a] * (a < 24 ? fit.r[a] : fit.rZ);
    auto propagate = [&](const double CC[25][25], double V3[3][3]) {
        for (int k = 0; k < 3; ++k)
            for (int l = 0; l < 3; ++l)
            {
                double s = 0;
                for (int a = 0; a < 25; ++a)
                    for (int b = 0; b < 25; ++b) s += J[k][a] * CC[a][b] * J[l][b];
                V3[k][l] = s;
            }
    };
    double V3[3][3], V3s[3][3];
    propagate(C, V3);
    if (haveStat) propagate(Cs, V3s);

    // (W, Z) 2x2: W = W+ + W-
    auto wz = [&](const double V[3][3], double &vWW, double &vZZ, double &vWZ) {
        vWW = V[0][0] + V[1][1] + 2 * V[0][1];
        vZZ = V[2][2];
        vWZ = V[0][2] + V[1][2];
    };
    double vWW, vZZ, vWZ, vWWs = 0, vZZs = 0, vWZs = 0;
    wz(V3, vWW, vZZ, vWZ);
    if (haveStat) wz(V3s, vWWs, vZZs, vWZs);

    const double sW = mu3[0] + mu3[1], sZ = mu3[2];
    const double eW = std::sqrt(std::max(0.0, vWW)), eZ = std::sqrt(std::max(0.0, vZZ));
    const double rho = (eW > 0 && eZ > 0) ? vWZ / (eW * eZ) : 0.0;
    const double eWs = std::sqrt(std::max(0.0, vWWs)), eZs = std::sqrt(std::max(0.0, vZZs));
    const double sysW = haveStat ? std::sqrt(std::max(0.0, vWW - vWWs)) : 0.0;
    const double sysZ = haveStat ? std::sqrt(std::max(0.0, vZZ - vZZs)) : 0.0;

    // theory point = r = 1 (POWHEG + the generation nPDF, EPPS21 central)
    double tW = 0;
    for (int i = 0; i < 24; ++i) tW += G[i];
    const double tZ = gZ;

    // the 106 EPPS21 member variations of that same prediction
    TGraph *gMem = nullptr;
    double mWlo = 0, mWhi = 0, mZlo = 0, mZhi = 0;
    const int nMem = epps21Scatter(genFile, isFB, gMem, mWlo, mWhi, mZlo, mZhi);
    if (nMem > 0)
        std::cout << Form("  [EPPS21] %d member variations: sigma_W %.2f..%.2f nb (%+.1f/%+.1f%%),"
                          " sigma_Z %.3f..%.3f nb (%+.1f/%+.1f%%) around the central point\n",
                          nMem, mWlo, mWhi, 100 * (mWlo / tW - 1), 100 * (mWhi / tW - 1),
                          mZlo, mZhi, 100 * (mZlo / tZ - 1), 100 * (mZhi / tZ - 1));

    std::cout << "\n===== xsec_contour_WZ (" << disc << ", " << bn << ") =====\n";
    std::cout << Form("  sigma_W = %.2f +/- %.2f nb", sW, eW);
    if (haveStat) std::cout << Form("  (%.2f stat, %.2f syst)", eWs, sysW);
    std::cout << Form("   [r_eff = %.4f]\n", tW > 0 ? sW / tW : 0.0);
    std::cout << Form("  sigma_Z = %.3f +/- %.3f nb", sZ, eZ);
    if (haveStat) std::cout << Form("  (%.3f stat, %.3f syst)", eZs, sysZ);
    std::cout << Form("   [r_Z = %.4f]\n", fit.rZ);
    std::cout << Form("  corr(sigma_W, sigma_Z) = %+.4f   cov = %+.4f nb^2\n", rho, vWZ);
    std::cout << Form("  theory (r=1, EPPS21 central): sigma_W = %.2f, sigma_Z = %.3f nb"
                      " (gen MC-stat on sigma_Z: %.4f nb, neglected)\n", tW, tZ, gZErr);
    // The ratio direction. The 3% lumi is fully correlated and CANCELS here
    // (that is what the +rho tilt of the ellipse is), so the ratio's systematic
    // collapses -- but sigma_Z is statistics-limited (~360 Z events), so the
    // ratio's TOTAL error is not smaller than sigma_W's. Report the split, not
    // just the total, or the cancellation is invisible.
    auto ratErr = [&](double vww, double vzz, double vwz) {
        return (sZ > 0 && sW > 0)
                   ? (sW / sZ) * std::sqrt(std::max(0.0, vww / (sW * sW) + vzz / (sZ * sZ)
                                                             - 2 * vwz / (sW * sZ)))
                   : 0.0;
    };
    const double ratio  = (sZ > 0) ? sW / sZ : 0;
    const double eRatio = ratErr(vWW, vZZ, vWZ);
    const double eRatioS = haveStat ? ratErr(vWWs, vZZs, vWZs) : 0.0;
    const double sysRatio = haveStat ? std::sqrt(std::max(0.0, eRatio * eRatio - eRatioS * eRatioS)) : 0.0;
    std::cout << Form("  sigma_W / sigma_Z = %.3f +/- %.3f (%.2f%%)", ratio, eRatio,
                      100 * eRatio / ratio);
    if (haveStat)
        std::cout << Form("  = %.3f stat (+) %.3f syst  -- syst %.2f%% vs %.2f%% on sigma_W"
                          " and %.2f%% on sigma_Z (the lumi cancellation)",
                          eRatioS, sysRatio, 100 * sysRatio / ratio, 100 * sysW / sW, 100 * sysZ / sZ);
    std::cout << "\n";

    // ---- 4) the plot ---------------------------------------------------------
    // Text is deliberately smaller than the repo default here (user request
    // 2026-09-15c): this panel carries the IN-FRAME CMS banner, a legend of up
    // to 7 entries, the numbers AND the EPPS21 member cloud, so the axis
    // titles/labels are shrunk to buy the vertical room.
    PlotStyle ps;
    ps.w = 800; ps.h = 800;
    ps.lm = 0.14; ps.rm = 0.05; ps.bm = 0.11; ps.tm = 0.06;
    ps.xTitleSize = ps.yTitleSize = 0.036;
    ps.xLabelSize = ps.yLabelSize = 0.030;
    ps.xTitleOffset = 1.20; ps.yTitleOffset = 1.60;

    TCanvas *c = new TCanvas("c_xsec_contour", "c_xsec_contour", ps.w, ps.h);
    ApplyCanvasStyle(c, ps);

    // axis range from the 95% total ellipse plus room for the theory point
    const double pad = 1.30;
    double xlo = std::min(sW - pad * 2.45 * eW, tW - 0.35 * eW);
    double xhi = std::max(sW + pad * 2.45 * eW, tW + 0.35 * eW);
    double ylo = std::min(sZ - pad * 2.45 * eZ, tZ - 0.35 * eZ);
    double yhi = std::max(sZ + pad * 2.45 * eZ, tZ + 0.35 * eZ);
    if (nMem > 0)   // the member cloud must fit too
    {
        xlo = std::min(xlo, mWlo); xhi = std::max(xhi, mWhi);
        ylo = std::min(ylo, mZlo); yhi = std::max(yhi, mZhi);
    }
    // 6% margin so a theory point that lands outside the ellipse is not clipped
    // by the frame (r_eff ~ 1.19 puts the r = 1 point ~5 sigma below in W)
    const double mx = 0.06 * (xhi - xlo), my = 0.06 * (yhi - ylo);
    xlo -= mx; xhi += mx; ylo -= my; yhi += my;
    TH1D *fr = new TH1D("fr_contour", "", 100, xlo, xhi);
    fr->SetMinimum(ylo); fr->SetMaximum(yhi);
    fr->SetStats(0);
    ApplyHistStyle(fr, ps, "#sigma^{fid}(W#rightarrow l#nu) [nb]",
                   "#sigma^{fid}(Z#rightarrow ll) [nb]");
    fr->Draw("AXIS");

    // 95% then 68% (2 dof: Delta chi2 = 5.991, 2.296)
    TGraph *e95 = ellipse(sW, sZ, vWW, vZZ, vWZ, 5.991);
    TGraph *e68 = ellipse(sW, sZ, vWW, vZZ, vWZ, 2.296);
    e95->SetFillColorAlpha(kAzure + 7, 0.25); e95->SetLineColor(kAzure + 3); e95->SetLineWidth(2);
    e68->SetFillColorAlpha(kAzure + 7, 0.55); e68->SetLineColor(kAzure + 3); e68->SetLineWidth(2);
    e95->Draw("F SAME"); e95->Draw("L SAME");
    e68->Draw("F SAME"); e68->Draw("L SAME");

    // Colour convention on this panel: BLACK = total, DARK RED = stat only;
    // DASHED = the Gaussian (covariance) region, SOLID = the profiled scan.
    TGraph *e68s = nullptr;
    if (haveStat)
    {
        e68s = ellipse(sW, sZ, vWWs, vZZs, vWZs, 2.296);
        e68s->SetLineColor(kRed + 2); e68s->SetLineWidth(2); e68s->SetLineStyle(7);
        e68s->Draw("L SAME");
    }
    // (the profiled stat-only twin is drawn with the other scan contours below,
    //  where it is in scope)

    TGraph *gBest = new TGraph(1); gBest->SetPoint(0, sW, sZ);
    gBest->SetMarkerStyle(20); gBest->SetMarkerSize(1.3); gBest->SetMarkerColor(kBlack);
    gBest->Draw("P SAME");

    // the EXACT profiled region, when a --contour scan exists
    std::vector<TGraph *> c68, c95;
    double pbW = 0, pbZ = 0;
    const TString scanDir = base + TString::Format("/../contour/contour_%s", bn.Data());
    const TString sScan = scanFile ? TString(scanFile) : resolveScanFile(scanDir, bn);
    bool haveScan = !sScan.IsNull() && profiledContours(sScan.Data(), gZ, c68, c95, pbW, pbZ);

    // ---- STALENESS GATE ------------------------------------------------------
    // The scan and the ellipse MUST come from the same fit. A reparametrization
    // cannot move the minimum, and the first row of a MultiDimFit grid is the
    // exact best fit (deltaNLL = 0), so the scan's minimum must reproduce
    // Sum_i r_i sigma_gen,i and r_Z sigma_gen,Z to float precision. If it does
    // not, the contour file belongs to a DIFFERENT fit (or to a different
    // gen_xsec.C run, since the sigma_gen weights are baked into the workspace)
    // -- and a contour drawn from one fit over an ellipse from another is worse
    // than no contour at all, so REFUSE to draw it and say exactly why.
    // Tolerance 1e-3 relative: the tree stores the POIs as Float_t (~1e-7) and
    // the CSV r's carry 6 digits, so anything real is orders of magnitude below.
    bool scanRejected = false;
    if (haveScan)
    {
        const double dW = (sW != 0) ? std::fabs(pbW / sW - 1.0) : 0.0;
        const double dZ = (sZ != 0) ? std::fabs(pbZ / sZ - 1.0) : 0.0;
        if (dW > 1e-3 || dZ > 1e-3)
        {
            std::cerr << Form("[ERROR] the MultiDimFit scan does NOT match this fit:"
                              " its minimum is (%.4f, %.5f) but the covariance best fit is"
                              " (%.4f, %.5f) -- relative %.2e / %.2e, tolerance 1e-3.\n",
                              pbW, pbZ, sW, sZ, dW, dZ)
                      << "        A reparametrization cannot move the minimum, so this scan is from a"
                         " DIFFERENT fit (or a different gen_xsec.C run: the sigma_gen weights are baked"
                         " into the contour workspace).\n"
                      << "        File: " << sScan << "\n"
                      << "        -> REFUSING to draw it. Re-run the fit with --contour and"
                         " `sync_lxplus.sh download`, or delete the stale file.\n";
            haveScan = false;
            scanRejected = true;
            c68.clear(); c95.clear();
        }
    }
    // the STAT-ONLY profiled twin (nuisances frozen at their post-fit values --
    // the same recipe as the --statonly companion fit), gated the same way
    std::vector<TGraph *> s68, s95;
    double sbW = 0, sbZ = 0;
    const TString sScanStat = resolveScanFile(scanDir, bn, "contourstat");
    bool haveScanStat = !sScanStat.IsNull() &&
                        profiledContours(sScanStat.Data(), gZ, s68, s95, sbW, sbZ);
    if (haveScanStat)
    {
        const double dW = (sW != 0) ? std::fabs(sbW / sW - 1.0) : 0.0;
        const double dZ = (sZ != 0) ? std::fabs(sbZ / sZ - 1.0) : 0.0;
        if (dW > 1e-3 || dZ > 1e-3)
        {
            std::cerr << Form("[ERROR] the STAT-ONLY scan does not match this fit: minimum"
                              " (%.4f, %.5f) vs (%.4f, %.5f), relative %.2e / %.2e"
                              " -> REFUSING to draw it.\n", sbW, sbZ, sW, sZ, dW, dZ);
            haveScanStat = false; s68.clear(); s95.clear();
        }
    }

    if (haveScan)
    {
        std::cout << "  [profiled] using " << sScan << "\n";
        for (size_t k = 0; k < c95.size(); ++k)
        { c95[k]->SetLineColor(kBlack); c95[k]->SetLineWidth(2); c95[k]->SetLineStyle(3); c95[k]->Draw("L SAME"); }
        for (size_t k = 0; k < c68.size(); ++k)
        { c68[k]->SetLineColor(kBlack); c68[k]->SetLineWidth(3); c68[k]->Draw("L SAME"); }
    }
    if (haveScanStat)
    {
        for (size_t k = 0; k < s68.size(); ++k)
        { s68[k]->SetLineColor(kRed + 2); s68[k]->SetLineWidth(3); s68[k]->Draw("L SAME"); }
        std::cout << Form("  [profiled] stat-only scan: minimum (%.4f, %.5f), d = %+.4f / %+.5f nb\n",
                          sbW, sbZ, sbW - sW, sbZ - sZ);
    }
    if (haveScan)
    {
        std::cout << Form("  [profiled] scan min vs covariance best fit:"
                          " d(sigma_W) = %+.4f nb, d(sigma_Z) = %+.5f nb"
                          " (must be ~0: a reparametrization cannot move the minimum)\n",
                          pbW - sW, pbZ - sZ);
        // HOW GAUSSIAN IS IT: every point of the profiled 68% contour should sit
        // at Mahalanobis distance^2 = 2.296 from the best fit if the likelihood
        // is exactly parabolic. Radius goes as sqrt(d), so sqrt(d/2.296) - 1 is
        // the fractional error the covariance ellipse makes in that direction --
        // i.e. the fractional error on the quoted Hesse uncertainty.
        // FLOOR: fed an exactly Gaussian scan on a 60x60 grid this reports
        // -0.4% .. -0.1% (measured 2026-09-15), which is the TGraph2D
        // interpolation + CONT LIST resolution, not physics. Do not read
        // anything below ~0.5% as non-Gaussianity; raise --points if you need to.
        const double detV = vWW * vZZ - vWZ * vWZ;
        if (detV > 0)
        {
            double dmin = 1e300, dmax = -1e300;
            long npts = 0;
            for (size_t k = 0; k < c68.size(); ++k)
                for (int p = 0; p < c68[k]->GetN(); ++p)
                {
                    double x, y; c68[k]->GetPoint(p, x, y);
                    const double dx = x - sW, dy = y - sZ;
                    const double d = (dx * dx * vZZ - 2 * dx * dy * vWZ + dy * dy * vWW) / detV;
                    dmin = std::min(dmin, d); dmax = std::max(dmax, d); ++npts;
                }
            if (npts > 0)
                std::cout << Form("  [profiled] Gaussianity: over %ld points of the 68%% contour the"
                                  " Mahalanobis d^2 runs %.3f..%.3f (exactly 2.296 if parabolic)"
                                  " -> the ellipse misstates the 1 sigma radius by %+.1f%% .. %+.1f%%\n",
                                  npts, dmin, dmax,
                                  100 * (std::sqrt(dmin / 2.296) - 1), 100 * (std::sqrt(dmax / 2.296) - 1));
        }
    }
    else if (!scanRejected)
        std::cout << "  [profiled] no MultiDimFit scan in " << scanDir
                  << " -> Gaussian ellipse only; the fork's `run_pO_fits.sh --disc " << disc
                  << "` writes one there (the --contour pass is on by default)\n";

    if (gMem)   // under the central marker, so the nominal stays readable
    {
        gMem->SetMarkerStyle(20);
        gMem->SetMarkerSize(0.6);
        gMem->SetMarkerColorAlpha(kGreen + 2, 0.7);
        gMem->Draw("P SAME");
    }
    TGraph *gTh = new TGraph(1); gTh->SetPoint(0, tW, tZ);
    gTh->SetMarkerStyle(33); gTh->SetMarkerSize(2.2); gTh->SetMarkerColor(kGreen + 3);
    gTh->Draw("P SAME");

    // Layout: the r = 1 prediction and its EPPS21 member cloud sit ~5 sigma
    // below the data in sigma_W, so the ellipses live in the upper right and
    // the left column is free. Stacked top to bottom: the in-frame CMS banner
    // (drawn last by CMS_lumi), one header line, the legend, the numbers --
    // stopping above the member cloud. Energy + luminosity come from CMS_lumi
    // (period 13 = "46.5 nb^{-1} (9.62 TeV pO)"), the repo's single source.
    ps.headerX = 0.17; ps.headerY = 0.845; ps.headerDy = 0.048;
    DrawHeader(ps, "", Form("simfit, %s, %s binning", discLabel.Data(), bn.Data()), "");

    const int nLeg = 4 + (e68s ? 1 : 0) + (haveScan ? 1 : 0) + (haveScanStat ? 1 : 0) + (gMem ? 1 : 0);
    TLegend *leg = new TLegend(0.16, 0.765 - 0.0335 * nLeg, 0.60, 0.765);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(ps.font); leg->SetTextSize(0.0245);
    leg->AddEntry(gBest, "Data (best fit)", "p");
    // label the ellipses as the GAUSSIAN approximation they are, so that a plot
    // WITHOUT the profiled contour cannot be mistaken for one that has it
    leg->AddEntry(e68, "68% CL (Gaussian)", "f");
    leg->AddEntry(e95, "95% CL (Gaussian)", "f");
    if (e68s) leg->AddEntry(e68s, "68% CL (Gaussian, stat. only)", "l");
    if (haveScan) leg->AddEntry(c68[0], "68%, 95% CL (profiled)", "l");
    if (haveScanStat) leg->AddEntry(s68[0], "68% CL (profiled, stat. only)", "l");
    leg->AddEntry(gTh, "POWHEG + EPPS21 (central)", "p");
    if (gMem) leg->AddEntry(gMem, Form("EPPS21 members (%d)", nMem), "p");
    leg->Draw();

    TLatex num;
    num.SetNDC(true); num.SetTextFont(ps.font); num.SetTextAlign(13); num.SetTextSize(0.026);
    std::vector<std::string> info;
    info.push_back(Form("#sigma_{W} = %.1f #pm %.1f nb", sW, eW));
    info.push_back(Form("#sigma_{Z} = %.2f #pm %.2f nb", sZ, eZ));
    info.push_back(Form("#rho = %+.3f", rho));
    info.push_back(Form("#sigma_{W}/#sigma_{Z} = %.2f #pm %.2f", ratio, eRatio));
    // hang the numbers off the legend's own bottom edge so they follow it when
    // an entry appears or disappears (stat ellipse, profiled contour, members)
    const double numTop = 0.765 - 0.0335 * nLeg - 0.028, numDy = 0.0345;
    for (size_t i = 0; i < info.size(); ++i)
        num.DrawLatex(0.18, numTop - numDy * i, info[i].c_str());
    // PROVENANCE STAMP: the plot must say, on its face, whether the profiled
    // contour is really there -- a missing legend entry is weak evidence, and a
    // rejected (stale) scan must not look like "never ran".
    num.SetTextSize(0.023);
    num.SetTextColor(haveScan ? kGreen + 3 : (scanRejected ? kRed + 1 : kGray + 2));
    num.DrawLatex(0.18, numTop - numDy * info.size(),
                  haveScan      ? "profiled contour: from the scan"
                  : scanRejected ? "profiled contour: REJECTED (other fit)"
                                 : "profiled contour: not available");
    num.SetTextColor(kBlack);

    // iPosX = 0 -> the CMS banner goes ABOVE the frame (out of frame). The
    // in-frame variant used elsewhere (10) would sit on the header/legend
    // stack: here the r = 1 point pins the ellipses to the upper right, so the
    // whole left column is needed for text.
    // iPosX = 10 -> the CMS banner sits INSIDE the frame, top left, as in every
    // other plot of the repo (out-of-frame was the 2026-09-15 stopgap). It fits
    // only because the text above was shrunk; the header/legend/numbers stack
    // starts below it.
    CMS_lumi(c, 13, 10);
    c->Update();

    const TString stem = TString::Format("%s/xsec_contour_WZ_%s", outDir.c_str(), bn.Data());
    c->SaveAs(stem + ".png");
    c->SaveAs(stem + ".pdf");

    // ---- 5) the machine-readable record --------------------------------------
    std::ofstream out((stem + ".csv").Data());
    out << "# (sigma_W, sigma_Z) from the simfit POIs via h_cov_poi; disc=" << disc
        << " binning=" << bn << "\n";
    out << "quantity,value_nb,err_total_nb,err_stat_nb,err_syst_nb\n";
    out << Form("sigma_Wp,%.5f,%.5f,%.5f,%.5f\n", mu3[0], std::sqrt(V3[0][0]),
                haveStat ? std::sqrt(V3s[0][0]) : 0.0,
                haveStat ? std::sqrt(std::max(0.0, V3[0][0] - V3s[0][0])) : 0.0);
    out << Form("sigma_Wm,%.5f,%.5f,%.5f,%.5f\n", mu3[1], std::sqrt(V3[1][1]),
                haveStat ? std::sqrt(V3s[1][1]) : 0.0,
                haveStat ? std::sqrt(std::max(0.0, V3[1][1] - V3s[1][1])) : 0.0);
    out << Form("sigma_W,%.5f,%.5f,%.5f,%.5f\n", sW, eW, eWs, sysW);
    out << Form("sigma_Z,%.5f,%.5f,%.5f,%.5f\n", sZ, eZ, eZs, sysZ);
    out << "\n# 3x3 covariance (nb^2), order [Wp, Wm, Z] -- total then stat\n";
    out << "matrix,row,Wp,Wm,Z\n";
    const char *rn[3] = {"Wp", "Wm", "Z"};
    for (int k = 0; k < 3; ++k)
        out << Form("total,%s,%.6f,%.6f,%.6f\n", rn[k], V3[k][0], V3[k][1], V3[k][2]);
    if (haveStat)
        for (int k = 0; k < 3; ++k)
            out << Form("stat,%s,%.6f,%.6f,%.6f\n", rn[k], V3s[k][0], V3s[k][1], V3s[k][2]);
    out << "\n# theory (r = 1: POWHEG + EPPS21nlo_CT18Anlo_O16 central)\n";
    out << Form("theory,sigma_W,%.5f\ntheory,sigma_Z,%.5f\n", tW, tZ);
    if (nMem > 0)
    {
        out << Form("theory,n_epps21_members,%d\n", nMem);
        out << Form("theory,sigma_W_member_min,%.5f\ntheory,sigma_W_member_max,%.5f\n", mWlo, mWhi);
        out << Form("theory,sigma_Z_member_min,%.5f\ntheory,sigma_Z_member_max,%.5f\n", mZlo, mZhi);
    }
    // provenance: which region the companion plot actually shows
    out << "\n# profiled contour (run_pO_fits.sh --contour)\n";
    out << "contour,status," << (haveScan ? "drawn" : (scanRejected ? "rejected_stale" : "absent")) << "\n";
    out << "contour,file," << (sScan.IsNull() ? TString("(none)") : sScan) << "\n";
    out << "contour,stat_status," << (haveScanStat ? "drawn" : (sScanStat.IsNull() ? "absent" : "rejected_stale")) << "\n";
    out << "contour,stat_file," << (sScanStat.IsNull() ? TString("(none)") : sScanStat) << "\n";
    if (haveScan || scanRejected)
        out << Form("contour,scan_min_sigma_W,%.5f\ncontour,scan_min_sigma_Z,%.5f\n", pbW, pbZ);
    out.close();
    std::cout << "[xsec_contour] wrote " << stem << ".{png,pdf,csv}\n";
}
