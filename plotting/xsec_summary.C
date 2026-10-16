// =============================================================================
// xsec_summary.C -- the fiducial cross sections W+, W-, W and Z, ONE PANEL EACH,
// stacked vertically, in the summary style of the CMS ttbar-in-PbPb paper
// (arXiv:2604.27091, Fig. 9): a horizontal sigma axis per panel, the mu+e
// combined measurement as a VERTICAL BAND (outer = total, inner = statistical,
// dashed line = central value) and every entry as a ROW, top to bottom:
//
//   mu only / e only / mu + e      the three simultaneous fits (fit_variants.h
//                                  styles): point + THICK stat bar + THIN total
//                                  bar -- the same numbers as W_fiducial, read
//                                  from xsec_flavfit.csv (xsec_fiducial_flav)
//   ---- dashed separator ----
//   POWHEG NLO x EPPS21            sigma_gen-fid of gen_xsec.root (member 0) with
//                                  its nPDF band: asymmetric Hessian over the
//                                  106 stored EPPS21 members (pOLhe::Hessian,
//                                  = LHAPDF PDFSet::uncertainty), 90% -> 68% CL
//   <extra theory rows>            optional, from theory_fid_predictions.txt
//                                  (other nPDF sets / NNLO -- to be added; the
//                                  format is documented in that file)
//   ---- dashed separator ----
//   CMS pp 5.02 / 13 TeV x A       (withPP) the pp fiducial measurements x A,
//                                  the coarse experimental "no nuclear effect"
//                                  reference of W_fiducial (no sqrt(s), isospin
//                                  or fiducial-definition correction)
//
// The measurement never moves when a theory row is added: sigma = r x
// sigma_gen-fid is a property of the data (the acceptance correction is taken
// from the EPPS21 MC once), so other nPDF sets enter ONLY as prediction rows.
//
// Inputs:  plots/flavfit/xsec/<disc>/xsec_flavfit.csv   (xsec_fiducial_flav)
//          ../skim/rootfile/gen_xsec.root               (skim/gen_xsec.C)
//          theory_fid_predictions.txt                   (optional extra rows)
// Outputs: plots/flavfit/xsec/<disc>/xsec_summary_vertical.{png,pdf,csv}
//   root -l -b -q -e 'gROOT->LoadMacro("xsec_summary.C+"); xsec_summary("leppt_mt40");'
// Run by analysis/run_observables.sh (flavfit block), after xsec_fiducial_flav.
// =============================================================================
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TBox.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TMath.h"
#include "TColor.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "plotting_helper.C"   // PlotStyle, CMS_lumi (via CMS_lumi.C)
#include "../skim/mc_norm.h"   // pONorm::kA_O
#include "../skim/lhe_index.h" // pOLhe::Hessian / Block (the EPPS21 member combination)
#include "disc_variants.h"     // pODisc::Spec
#include "fit_variants.h"      // pOFit::Spec (fit labels / markers / colours)

namespace
{
enum Kind { kFitRow = 0, kTheoryRow = 1, kPPRow = 2 };

struct Row
{
    TString label;
    Kind kind = kFitRow;
    double v[4] = {0, 0, 0, 0};
    double lo[4] = {0, 0, 0, 0}, hi[4] = {0, 0, 0, 0};   // the (outer / only) bar, asymmetric
    double slo[4] = {0, 0, 0, 0}, shi[4] = {0, 0, 0, 0}; // the stat bar (fits only)
    bool has[4] = {false, false, false, false};
    bool hasStat = false;
    int color = kBlack, marker = 20;
    double msize = 1.3;
    TString note; // CSV: where the numbers come from
};

const char *kProc[4]    = {"Wp", "Wm", "W", "Z"};
const char *kProcTex[4] = {"W^{+} #rightarrow l^{+}#nu", "W^{-} #rightarrow l^{-}#bar{#nu}",
                           "W^{#pm} #rightarrow l^{#pm}#nu", "Z #rightarrow l^{+}l^{-}"};
const char *kAxisTex[4] = {"#sigma^{fid}_{W^{+}} (nb)", "#sigma^{fid}_{W^{-}} (nb)",
                           "#sigma^{fid}_{W} (nb)", "#sigma^{fid}_{Z} (nb)"};

// the exact LHAPDF 90% -> 68.27% CL factor (= 0.607957; "/1.645" is a rounding)
const double kCL90to68 = std::sqrt(TMath::ChisquareQuantile(0.682689492137, 1) /
                                   TMath::ChisquareQuantile(0.90, 1));

std::vector<std::string> splitKeep(const std::string &line, char sep)
{
    std::vector<std::string> out;
    std::string cur;
    std::istringstream ss(line);
    while (std::getline(ss, cur, sep)) out.push_back(cur);
    if (!line.empty() && line.back() == sep) out.push_back("");
    return out;
}
std::string trim(const std::string &s)
{
    const size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
    return (a == std::string::npos) ? "" : s.substr(a, b - a + 1);
}
int procIndex(const std::string &c)
{
    for (int k = 0; k < 4; ++k)
        if (c == kProc[k]) return k;
    return -1;
}

// ---- the inclusive rows of xsec_flavfit.csv for one tag ----------------------
// columns: fit,charge,ybin,etaCM_center,etaCM_halfwidth,r,rErr,sigma_gen_fid_nb,
//          sigma_meas_nb,sigma_meas_err_nb,sigma_meas_stat_nb,sigma_meas_syst_nb
// (pp reference rows carry only sigma_meas_nb and sigma_meas_err_nb)
bool readInclRows(const char *csv, const char *tag, Row &r)
{
    std::ifstream in(csv);
    if (!in.is_open()) return false;
    std::string line;
    bool any = false;
    while (std::getline(in, line))
    {
        const std::vector<std::string> f = splitKeep(line, ',');
        if (f.size() < 10 || f[0] != tag || f[2] != "incl") continue;
        const int k = procIndex(f[1]);
        if (k < 0) continue;
        const double v = std::atof(f[8].c_str()), e = std::atof(f[9].c_str());
        r.v[k] = v; r.lo[k] = r.hi[k] = e; r.has[k] = true;
        if (f.size() >= 11 && !trim(f[10]).empty())
        {
            r.slo[k] = r.shi[k] = std::atof(f[10].c_str());
            r.hasStat = true;
        }
        any = true;
    }
    return any;
}

// ---- POWHEG NLO x EPPS21: nominal + the nPDF band from the stored members ----
// `disc` picks the gen fiducial set (pODisc::GenFid), so the prediction sits
// in the same volume as the measurement rows (m_T > 40 for leppt_mt40).
bool theoryEpps21(const char *genFile, const char *disc, Row &r, std::vector<std::string> &report)
{
    TString gtag, glabel;
    if (!pODisc::GenFid(disc, gtag, glabel)) return false;
    TFile *f = TFile::Open(genFile, "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "[ERROR] xsec_summary: cannot open " << genFile << " (run skim/run_gen_xsec.sh)\n";
        if (f) { f->Close(); delete f; }
        return false;
    }
    TH1D *hp = (TH1D *)f->Get("h_gen_sig_Wp" + gtag), *hm = (TH1D *)f->Get("h_gen_sig_Wm" + gtag), *hz = (TH1D *)f->Get("h_gen_sig_Z");
    TH2D *hp2 = (TH2D *)f->Get("h_gen_sig_Wp" + gtag + "_epps21"), *hm2 = (TH2D *)f->Get("h_gen_sig_Wm" + gtag + "_epps21");
    TH1D *hz2 = (TH1D *)f->Get("h_gen_sig_Z_epps21");
    if (!hp || !hm)
    {
        std::cerr << "[ERROR] xsec_summary: no h_gen_sig_Wp/Wm" << gtag << " in " << genFile << "\n";
        f->Close(); delete f; return false;
    }
    // nominal
    double nom[4] = {0, 0, 0, 0};
    for (int i = 1; i <= hp->GetNbinsX(); ++i) { nom[0] += hp->GetBinContent(i); nom[1] += hm->GetBinContent(i); }
    nom[2] = nom[0] + nom[1];
    nom[3] = hz ? hz->GetBinContent(1) : 0.0;
    const int nm = (hp2 && hm2) ? hp2->GetNbinsY() : 0; // 107 members, ROOT bin 1 = member 0
    for (int k = 0; k < 4; ++k)
    {
        if (nom[k] <= 0) continue;
        r.v[k] = nom[k]; r.has[k] = true;
        if (nm < 3 || (k == 3 && !hz2)) { r.lo[k] = r.hi[k] = 0.0; continue; }
        std::vector<double> S(nm, 0.0);
        for (int m = 0; m < nm; ++m)
        {
            if (k == 3) S[m] = hz2->GetBinContent(m + 1);
            else
                for (int i = 1; i <= hp2->GetNbinsX(); ++i)
                    S[m] += (k != 1 ? hp2->GetBinContent(i, m + 1) : 0.0) + (k != 0 ? hm2->GetBinContent(i, m + 1) : 0.0);
        }
        auto dev = [&](int m) { return S[0] > 0 ? S[m] / S[0] - 1.0 : 0.0; };
        pOLhe::Block bAll{1, nm - 1, pOLhe::kHessian, "", "", "", ""};
        pOLhe::Block bNuc{1, 48, pOLhe::kHessian, "", "", "", ""};
        pOLhe::Block bBas{49, nm - 1, pOLhe::kHessian, "", "", "", ""};
        const pOLhe::HessRes hA = pOLhe::Hessian(bAll, dev), hN = pOLhe::Hessian(bNuc, dev), hB = pOLhe::Hessian(bBas, dev);
        r.lo[k] = nom[k] * hA.dn * kCL90to68;
        r.hi[k] = nom[k] * hA.up * kCL90to68;
        const double closure = std::fabs(S[0] - nom[k]) / nom[k];
        report.push_back(Form("[theory] %-2s POWHEG NLO x EPPS21: %8.3f nb  nPDF 68%% CL +%.3f -%.3f (+%.2f%% -%.2f%%)  [nuclear +%.2f/-%.2f%%, baseline +%.2f/-%.2f%%; 90%% CL x %.4f]  member0-vs-nominal %.1e",
                              kProc[k], nom[k], r.hi[k], r.lo[k], 100 * hA.up * kCL90to68, 100 * hA.dn * kCL90to68,
                              100 * hN.up * kCL90to68, 100 * hN.dn * kCL90to68, 100 * hB.up * kCL90to68, 100 * hB.dn * kCL90to68,
                              kCL90to68, closure));
        if (closure > 1e-6)
            std::cerr << "[WARN] xsec_summary: EPPS21 member 0 differs from the nominal gen sigma for " << kProc[k]
                      << " by " << closure << " (expected 0 -- stale gen_xsec.root?)\n";
    }
    if (nm < 3)
        std::cerr << "[WARN] xsec_summary: no EPPS21 member twins in " << genFile
                  << " -> the theory row has no nPDF band (re-run skim/gen_xsec.C)\n";
    f->Close(); delete f;
    return r.has[0] && r.has[1];
}

// ---- extra theory rows: "label | Wp dn up | Wm dn up | W dn up | Z dn up" -----
int readExtraTheory(const char *file, std::vector<Row> &rows)
{
    std::ifstream in(file);
    if (!in.is_open()) return 0;
    const int pal[4] = {kGray + 3, kViolet + 1, kOrange + 7, kTeal + 3};
    const int mk[4] = {24, 25, 27, 28}; // open circle / square / diamond / cross (no triangles)
    const double ms[4] = {1.3, 1.5, 2.0, 1.5};
    std::string line;
    int n = 0;
    while (std::getline(in, line))
    {
        const std::string t = trim(line);
        if (t.empty() || t[0] == '#') continue;
        const std::vector<std::string> blk = splitKeep(t, '|');
        if (blk.size() < 5)
        {
            std::cerr << "[WARN] xsec_summary: skipping malformed theory line (need label | Wp | Wm | W | Z): " << t << "\n";
            continue;
        }
        Row r;
        r.kind = kTheoryRow;
        r.label = trim(blk[0]).c_str();
        r.note = file;
        for (int k = 0; k < 4; ++k)
        {
            std::istringstream ss(trim(blk[k + 1]));
            double v, dn, up;
            if (ss >> v >> dn >> up) { r.v[k] = v; r.lo[k] = std::fabs(dn); r.hi[k] = std::fabs(up); r.has[k] = true; }
        }
        r.color = pal[n % 4]; r.marker = mk[n % 4]; r.msize = ms[n % 4];
        rows.push_back(r);
        ++n;
    }
    return n;
}
} // namespace

void xsec_summary(const char *disc = "leppt_mt40",
                  bool withPP = true,
                  const char *extraTheoryFile = "theory_fid_predictions.txt",
                  const char *genFile = "../skim/rootfile/gen_xsec.root")
{
    TString dsuf, discLabel;
    if (!pODisc::Spec(disc, dsuf, discLabel)) return;
    const std::string outDir = std::string("./plots/flavfit/xsec/") + disc;
    const std::string csv = outDir + "/xsec_flavfit.csv";
    if (gSystem->AccessPathName(csv.c_str()))
    {
        std::cerr << "[ERROR] xsec_summary: " << csv << " not found -- run xsec_fiducial_flav(\"" << disc
                  << "\") first (analysis/run_observables.sh does)\n";
        return;
    }

    // ---- the rows, in drawing order ------------------------------------------
    std::vector<Row> rows;
    int nFit = 0, nTheory = 0, nPP = 0;
    const char *fitTags[3] = {"simfit_mu", "simfit_ele", "comb"};
    int iComb = -1;
    for (int t = 0; t < 3; ++t)
    {
        pOFit::Spec s;
        if (!pOFit::Get(fitTags[t], s)) continue;
        Row r;
        if (!readInclRows(csv.c_str(), fitTags[t], r))
        {
            std::cerr << "[WARN] xsec_summary: no '" << fitTags[t] << "' inclusive rows in " << csv << " -> row not drawn\n";
            continue;
        }
        r.kind = kFitRow;
        r.label = s.label; // "mu only" / "e only" / "mu + e combined" -- "simfit" is in the legend
        r.color = s.color; r.marker = s.marker; r.msize = s.markerSize;
        r.note = "xsec_flavfit.csv";
        if (TString(fitTags[t]) == "comb") iComb = (int)rows.size();
        rows.push_back(r);
        ++nFit;
    }
    if (nFit == 0)
    {
        std::cerr << "[ERROR] xsec_summary: no fit rows in " << csv << "\n";
        return;
    }
    std::vector<std::string> report;
    {
        Row r;
        if (theoryEpps21(genFile, disc, r, report))
        {
            r.kind = kTheoryRow;
            r.label = "POWHEG NLO #otimes EPPS21"; // the CT18ANLO baseline is named in the footnote
            r.color = kGreen + 2; r.marker = 20; r.msize = 1.3;
            r.note = "gen_xsec.root member 0, nPDF 68% CL (53 Hessian pairs)";
            rows.push_back(r);
            ++nTheory;
        }
    }
    nTheory += readExtraTheory(extraTheoryFile, rows);
    if (withPP)
    {
        const char *ppTag[2] = {"pp_5.02TeV", "pp_13TeV"};
        const char *ppLab[2] = {"pp 5.02 TeV #times A", "pp 13 TeV #times A"}; // CMS + the citation in the footnote
        const int ppCol[2] = {kOrange + 7, kMagenta + 2}, ppMk[2] = {24, 25}; // = kPPRef of xsec_fiducial.C
        const double ppMs[2] = {1.3, 1.5};
        for (int p = 0; p < 2; ++p)
        {
            Row r;
            if (!readInclRows(csv.c_str(), ppTag[p], r)) continue;
            r.kind = kPPRow; r.hasStat = false;
            r.label = ppLab[p];
            r.color = ppCol[p]; r.marker = ppMk[p]; r.msize = ppMs[p];
            r.note = "CMS JHEP 04 (2025) 162, pp fiducial x A (total unc.)";
            rows.push_back(r);
            ++nPP;
        }
    }
    const int nRows = (int)rows.size();

    // ---- console record --------------------------------------------------------
    for (const std::string &l : report) std::cout << l << "\n";
    for (int k = 0; k < 4; ++k)
        for (const Row &r : rows)
        {
            if (!r.has[k]) continue;
            if (r.hasStat)
                printf("[summary] %-2s %-40s %8.3f  +%.3f -%.3f (total)  +/- %.3f (stat) nb\n", kProc[k], r.label.Data(), r.v[k], r.hi[k], r.lo[k], r.shi[k]);
            else
                printf("[summary] %-2s %-40s %8.3f  +%.3f -%.3f nb\n", kProc[k], r.label.Data(), r.v[k], r.hi[k], r.lo[k]);
        }

    // ---- the canvas: header | legend | 4 panels | footnote ---------------------
    gSystem->mkdir(outDir.c_str(), kTRUE);
    gStyle->SetEndErrorSize(3);
    PlotStyle ps;
    const int cw = 800, ch = 1500;
    TCanvas *c = new TCanvas("c_xsec_summary", "", cw, ch);
    c->SetFillColor(0); c->SetBorderMode(0); c->SetFrameBorderMode(0);
    const double lmP = 0.035, rmP = 0.035;   // panel margins (the frame spans the width, as in the paper)
    const double yHdr0 = 0.935, yLeg0 = 0.880, yFoot = 0.075;
    const double panelH = (yLeg0 - yFoot) / 4.0;

    const int colTot = TColor::GetColor("#4F81BD"), colStat = TColor::GetColor("#F2A93B");
    const int colSep = kAzure + 1;

    // the band = the mu+e combined fit
    const Row *comb = (iComb >= 0) ? &rows[iComb] : nullptr;

    std::vector<TObject *> keep; // ROOT owns nothing drawn here; keep the pointers alive
    for (int k = 0; k < 4; ++k)
    {
        // rows present in this panel, in order
        std::vector<int> idx;
        for (int j = 0; j < nRows; ++j)
            if (rows[j].has[k]) idx.push_back(j);
        const int nR = (int)idx.size();
        if (nR == 0) continue;

        const double y1 = yFoot + (3 - k) * panelH, y2 = y1 + panelH;
        TPad *pad = new TPad(Form("pad_%s", kProc[k]), "", 0.0, y1, 1.0, y2);
        pad->SetMargin(lmP, rmP, 0.19, 0.03);
        pad->SetTicks(1, 0);
        pad->SetFillColor(0); pad->SetBorderMode(0); pad->SetFrameBorderMode(0);
        pad->SetFrameLineWidth(2);
        c->cd();
        pad->Draw();
        pad->cd();
        keep.push_back(pad);

        // x range: the data occupy the right 60% of the frame, the row labels the
        // left 40% (the paper's layout: labels inside the frame, over the low end
        // of the axis); a cross section axis never starts below zero
        double dlo = 1e300, dhi = -1e300;
        for (int j : idx)
        {
            dlo = std::min(dlo, rows[j].v[k] - rows[j].lo[k]);
            dhi = std::max(dhi, rows[j].v[k] + rows[j].hi[k]);
        }
        double span = dhi - dlo;
        if (span < 0.04 * std::fabs(dhi)) span = 0.04 * std::fabs(dhi);
        const double dataLo = dlo - 0.08 * span, dataHi = dhi + 0.08 * span;
        const double xhi = dataHi;
        const double xlo = std::max(0.0, dataHi - (dataHi - dataLo) / 0.60);
        const double yTop = nR + 1.6, yBot = 0.4; // room above the first row for the process label

        TH1F *fr = new TH1F(Form("fr_summary_%s", kProc[k]), "", 10, xlo, xhi);
        fr->SetStats(0);
        fr->SetMinimum(yBot); fr->SetMaximum(yTop);
        fr->GetXaxis()->SetTitle(kAxisTex[k]);
        fr->GetXaxis()->SetTitleFont(ps.font); fr->GetXaxis()->SetLabelFont(ps.font);
        fr->GetXaxis()->SetTitleSize(0.075); fr->GetXaxis()->SetLabelSize(0.065);
        fr->GetXaxis()->SetTitleOffset(1.00);
        fr->GetXaxis()->SetNdivisions(506);
        fr->GetYaxis()->SetNdivisions(0);
        fr->GetYaxis()->SetLabelSize(0); fr->GetYaxis()->SetTickLength(0);
        fr->Draw("AXIS");
        keep.push_back(fr);

        // the combined measurement as a band: total (outer), stat (inner), central value
        if (comb && comb->has[k])
        {
            const double yb1 = yBot, yb2 = nR + 0.6;
            TBox *bT = new TBox(comb->v[k] - comb->lo[k], yb1, comb->v[k] + comb->hi[k], yb2);
            bT->SetFillColor(colTot); bT->SetLineColor(colTot); bT->Draw();
            TBox *bS = new TBox(comb->v[k] - comb->slo[k], yb1, comb->v[k] + comb->shi[k], yb2);
            bS->SetFillColor(colStat); bS->SetLineColor(colStat);
            if (comb->hasStat) bS->Draw();
            TLine *lc = new TLine(comb->v[k], yb1, comb->v[k], yb2);
            lc->SetLineStyle(2); lc->SetLineWidth(2); lc->SetLineColor(kBlack); lc->Draw();
            keep.push_back(bT); keep.push_back(bS); keep.push_back(lc);
        }
        // separators between the row groups
        {
            int nF = 0, nT = 0;
            for (int j : idx) { if (rows[j].kind == kFitRow) ++nF; if (rows[j].kind == kTheoryRow) ++nT; }
            const int nP = nR - nF - nT;
            std::vector<double> ys;
            if (nF > 0 && nT + nP > 0) ys.push_back(nR - nF + 0.5);
            if (nT > 0 && nP > 0) ys.push_back(nR - nF - nT + 0.5);
            for (double ysep : ys)
            {
                TLine *ls = new TLine(xlo, ysep, xhi, ysep);
                ls->SetLineStyle(2); ls->SetLineWidth(1); ls->SetLineColor(colSep); ls->Draw();
                keep.push_back(ls);
            }
        }
        // the rows: thin total bar (caps), thick stat bar, marker, label
        TLatex lab; lab.SetTextFont(ps.font); lab.SetTextAlign(12); lab.SetTextSize(0.066);
        const double xLab = xlo + 0.015 * (xhi - xlo);
        for (int q = 0; q < nR; ++q)
        {
            const Row &r = rows[idx[q]];
            const double y = nR - q;
            double xv = r.v[k], yv = y, exl = r.lo[k], exh = r.hi[k], zero = 0.0;
            TGraphAsymmErrors *gT = new TGraphAsymmErrors(1, &xv, &yv, &exl, &exh, &zero, &zero);
            gT->SetLineColor(r.color); gT->SetMarkerColor(r.color);
            gT->SetMarkerStyle(r.marker); gT->SetMarkerSize(r.msize);
            if (r.hasStat)
            {
                gT->SetLineWidth(1);
                gT->Draw("P SAME"); // thin total bar with end caps, marker on top later
                double sl = r.slo[k], sh = r.shi[k];
                TGraphAsymmErrors *gS = new TGraphAsymmErrors(1, &xv, &yv, &sl, &sh, &zero, &zero);
                gS->SetLineColor(r.color); gS->SetMarkerColor(r.color);
                gS->SetMarkerStyle(r.marker); gS->SetMarkerSize(r.msize);
                gS->SetLineWidth(4);
                gS->Draw("PZ SAME"); // thick stat bar, no caps
                keep.push_back(gS);
            }
            else
            {
                gT->SetLineWidth(r.kind == kTheoryRow ? 4 : 2);
                gT->Draw(r.kind == kTheoryRow ? "PZ SAME" : "P SAME");
            }
            keep.push_back(gT);
            lab.SetTextColor(r.kind == kPPRow ? kGray + 3 : kBlack);
            lab.DrawLatex(xLab, y, r.label);
        }
        // the process, above the first row
        TLatex tp; tp.SetTextFont(62); tp.SetTextAlign(12); tp.SetTextSize(0.080);
        tp.DrawLatex(xLab, nR + 1.08, kProcTex[k]);
        pad->RedrawAxis();
    }

    // ---- legend strip -----------------------------------------------------------
    c->cd();
    TPad *pLeg = new TPad("pad_leg", "", 0.0, yLeg0, 1.0, yHdr0);
    pLeg->SetMargin(0, 0, 0, 0); pLeg->SetFillColor(0); pLeg->SetBorderMode(0);
    pLeg->Draw(); pLeg->cd();
    keep.push_back(pLeg);
    {
        TLegend *leg = new TLegend(lmP, 0.02, 1.0 - rmP, 0.98);
        leg->SetNColumns(2); leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->SetTextFont(ps.font); leg->SetTextSize(0.26);
        TBox *kT = new TBox(0, 0, 1, 1); kT->SetFillColor(colTot); kT->SetLineColor(colTot);
        TBox *kS = new TBox(0, 0, 1, 1); kS->SetFillColor(colStat); kS->SetLineColor(colStat);
        double x0 = 0, y0 = 0, e0 = 0;
        TGraphAsymmErrors *kB = new TGraphAsymmErrors(1, &x0, &y0, &e0, &e0, &e0, &e0);
        kB->SetLineWidth(4); kB->SetLineColor(kBlack); kB->SetMarkerStyle(20); kB->SetMarkerSize(1.2);
        TGraphAsymmErrors *kN = new TGraphAsymmErrors(1, &x0, &y0, &e0, &e0, &e0, &e0);
        kN->SetLineWidth(4); kN->SetLineColor(kGreen + 2); kN->SetMarkerStyle(20); kN->SetMarkerColor(kGreen + 2); kN->SetMarkerSize(1.2);
        leg->AddEntry(kT, "#mu + e simfit: total unc.", "f");
        leg->AddEntry(kB, "fits: stat. (thick), total (thin) bars", "lp");
        leg->AddEntry(kS, "#mu + e simfit: stat. unc.", "f");
        leg->AddEntry(kN, "theory: nPDF unc. (EPPS21, 68% CL)", "lp");
        leg->Draw();
        keep.push_back(leg); keep.push_back(kT); keep.push_back(kS); keep.push_back(kB); keep.push_back(kN);
    }

    // ---- header: CMS + lumi above everything (the paper's banner) -----------------
    c->cd();
    TPad *pHdr = new TPad("pad_hdr", "", 0.0, yHdr0, 1.0, 1.0);
    pHdr->SetMargin(lmP, rmP, 0.0, 0.55); pHdr->SetFillColor(0); pHdr->SetBorderMode(0);
    pHdr->Draw(); pHdr->cd();
    keep.push_back(pHdr);
    {
        const float savedRelPosX = relPosX;
        relPosX = 0.135; // keep "Work in Progress" clear of "CMS" in the out-of-frame layout
        CMS_lumi(pHdr, 13, 0);
        relPosX = savedRelPosX;
    }

    // ---- footnote: the fiducial definitions --------------------------------------
    c->cd();
    TPad *pFoot = new TPad("pad_foot", "", 0.0, 0.0, 1.0, yFoot);
    pFoot->SetMargin(lmP, rmP, 0, 0); pFoot->SetFillColor(0); pFoot->SetBorderMode(0);
    pFoot->Draw(); pFoot->cd();
    keep.push_back(pFoot);
    {
        TLatex ft; ft.SetNDC(); ft.SetTextFont(ps.font); ft.SetTextAlign(13); ft.SetTextSize(0.15); ft.SetTextColor(kGray + 3);
        const double dyF = 0.185;
        double yf = 0.96;
        TString gtag, glabel;
        pODisc::GenFid(disc, gtag, glabel);   // the W fiducial set of this discriminant
        ft.DrawLatex(lmP, yf, "pO fiducial (per lepton flavour):  W: bare lepton " + glabel); yf -= dyF;
        ft.DrawLatex(lmP, yf, "Z: p_{T}^{l} > 15 / 10 GeV, |#eta_{lab}| < 2.4, 60 < m_{ll} < 120 GeV"); yf -= dyF;
        ft.DrawLatex(lmP, yf, Form("#sigma = r #times #sigma^{gen}_{fid} (POWHEG NLO, EPPS21nlo_CT18Anlo_O16), %s fit", discLabel.Data())); yf -= dyF;
        ft.DrawLatex(lmP, yf, "syst. = profiled fit nuisances incl. lumi 3%;  theory band = EPPS21 nPDF, 68% CL (53 Hessian directions)"); yf -= dyF;
        if (nPP > 0)
            ft.DrawLatex(lmP, yf, Form("pp: CMS JHEP 04 (2025) 162, pp fiducial #times A = %.0f; no #sqrt{s}, isospin or fiducial-definition correction", pONorm::kA_O));
    }

    c->cd(); c->Modified(); c->Update();
    c->SaveAs((outDir + "/xsec_summary_vertical.png").c_str());
    c->SaveAs((outDir + "/xsec_summary_vertical.pdf").c_str());
    std::cout << "[OK] Saved " << outDir << "/xsec_summary_vertical.{png,pdf}  (disc=" << disc << ", rows: "
              << nFit << " fits + " << nTheory << " theory + " << nPP << " pp)\n";

    // ---- CSV: every number drawn (row labels in plain text, no TLatex markup) -------
    auto plainLabel = [](const TString &l) {
        TString p = l;
        p.ReplaceAll("#mu", "mu"); p.ReplaceAll("#otimes", "x"); p.ReplaceAll("#times", "x");
        p.ReplaceAll(",", ";");
        return p;
    };
    std::ofstream out((outDir + "/xsec_summary_vertical.csv").c_str());
    out << "process,row,kind,value_nb,err_lo_nb,err_hi_nb,stat_lo_nb,stat_hi_nb,source\n";
    for (int k = 0; k < 4; ++k)
        for (const Row &r : rows)
        {
            if (!r.has[k]) continue;
            out << kProc[k] << "," << plainLabel(r.label) << "," << (r.kind == kFitRow ? "fit" : r.kind == kTheoryRow ? "theory" : "pp")
                << "," << Form("%.5f,%.5f,%.5f,", r.v[k], r.lo[k], r.hi[k])
                << (r.hasStat ? Form("%.5f,%.5f", r.slo[k], r.shi[k]) : ",") << "," << r.note << "\n";
        }
    out.close();
    std::cout << "[OK] Wrote " << outDir << "/xsec_summary_vertical.csv\n";
}
