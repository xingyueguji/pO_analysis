// =============================================================================
// plotting/compare_fits.C -- TWO FITS OF THE SAME DATA SIDE BY SIDE (2026-10-19)
//
//   root -l -b -q 'compare_fits.C+("leppt_mt40","leppt_mt40_abcdy")'
//   (run it through plotting/run_compare_fits.sh -> logs/compare_fits_<ref>_vs_<alt>.log)
//
// First use: the rapidity-resolved in-fit ABCD (the fork's opt-in
// QCD_MODE=abcdy, downstream tag leppt_mt40_abcdy) against the nominal in-fit
// ABCD (QCD_MODE=abcd, tag leppt_mt40). Both arguments are downstream tags of
// plotting/disc_variants.h: each names an out-tree pO_fit_out<suffix>/ of the
// fork AND the observables analysis/run_observables.sh <tag> built from it --
// run that for both tags first. For every fit (grand mu+e, mu-only, e-only)
// and both binnings (lab, fb) the log prints, and compare_fits.csv records:
//   - fit status / covQual (summary CSV), the Asimov closure (extraction log)
//     and the multijet model (datacards/qcd_lnn_kappas.txt, qcdMode)
//   - the SR MULTIJET per rapidity group (correction/abcd_ygroups.h), prefit
//     and postfit, summed from the fit's OWN fitDiagnostics shapes
//     (shapes_prefit / shapes_fit_s, process qcd, the SR channels
//     <F>_<C>_<bn>_y<i>). That is the like-for-like quantity: the models'
//     multipliers act on different prefit templates (qcd_abcd vs qcd_abcdy),
//     so they are printed as well but compared through the yields. The
//     postfit error is the multiplier's relative error (the qcdFF_ele shape
//     morphing and the residual qcd_rate lnN are inside the multiplier only
//     as far as Combine's κ^θ·sB·sC/sD bookkeeping puts them there).
//   - every nuisance pull (<name>_theta), r per bin and r_Z (total, stat)
//   - sigma W+/W-/W/Z (the xsec_contour_WZ_<bn>[_<flav>].csv of each tag:
//     r x sigma_gen with the 25x25 POI covariance), the e/mu ratios with stat
//     errors and the (e/mu)_W+ / (e/mu)_W- double ratio
//   - A_ch and R_FB (sum, W+, W-) per bin from the FIDUCIAL graphs
//     ../skim/rootfile/{charge_asym,FBratio}_fid_<fit>_<tag>.root -- every
//     observable from r x sigma_gen, never from raw fitted counts.
// Differences are delta = alt - ref and delta / sigma_stat(ref). The data are
// the same, so that is the model shift in units of the statistical error, not
// a compatibility test.
//
// PROVENANCE CHECK: per fit and (flavour, charge) the fb prefit multijet total
// over the lab one is printed. = 1.000 means the fb templates carry the full
// y-inclusive multijet, i.e. the fit predates the 2026-10-19 FB-window fix
// (docs/iterations/2026-10-19_fbqcd.md; the FB window [-1.7068, 2.4] holds
// only part of it); an fb comparison against such a fit mixes that fix with
// the model change. The lab inputs were not touched by that fix.
//
// Outputs: plots/compare_fits/<ref>_vs_<alt>/
//   compare_fits.csv                       every number of the log (long format)
//   qcd_groups_<fit>_<bn>.{png,pdf}        multijet per group / ref prefit:
//                                          ref postfit, alt prefit, alt postfit
//   chargeAsym_<fit>, RFB_{sum,Wp,Wm}_<fit>   ref vs alt overlays (stat bars + syst boxes)
// Lines starting with [HEAD] are the headline numbers (the wrapper echoes them).
// =============================================================================
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"
#include "TSystem.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "plotting_helper.C"            // PlotStyle, SaveNiceGraph_Overlay
#include "disc_variants.h"              // pODisc::Spec / GraphFile
#include "fit_variants.h"               // pOFit::Spec / WorkDir (honours $FORK_TEST)
#include "../correction/abcd_ygroups.h" // pOAbcdY: the rapidity groups

namespace
{

struct Val
{
    double v = 0, e = 0, s = -1; // s < 0: no statistical error known
    bool ok = false;
};

const char *kBn[2] = {"lab", "fb"};
const char *kCh[2] = {"Wp", "Wm"};
const char *kSigQ[4] = {"sigma_Wp", "sigma_Wm", "sigma_W", "sigma_Z"};
const char *kRfb[3] = {"sum", "Wp", "Wm"};
// the groups in RAPIDITY order (bwd y0-1, cen y2-9, fwd y10-11); indices into pOAbcdY::kGroupName
const int kGOrd[3] = {1, 2, 0};

struct FitData
{
    bool ok = false;
    TString tag, fitTag, work, model = "?";
    pOFit::Spec spec;
    std::map<std::string, Val> par[2];
    int status[2] = {-1, -1}, covq[2] = {-1, -1};
    std::string asimov[2] = {"absent", "absent"};
    bool haveShapes[2] = {false, false};
    double qPre[2][2][2][3] = {}, qPost[2][2][2][3] = {}; // [bn][flav][charge][group]
    double qPreTot[2][2][2] = {};                         // [bn][flav][charge], all 12 bins
    Val sig[2][4];                                        // [bn][Wp, Wm, W, Z]
    TGraphErrors *gA = nullptr, *gAs = nullptr, *gR[3] = {}, *gRs[3] = {};
};

int FlavIndex(const TString &f) { return f == "mu" ? 0 : 1; }

std::vector<std::string> Split(const std::string &line)
{
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, ',')) out.push_back(tok);
    return out;
}

TString ModelLabel(const TString &m)
{
    if (m == "abcd") return "in-fit ABCD per charge";
    if (m == "abcdy") return "in-fit ABCD per y group";
    if (m == "mtsplit") return "m_T split, multijet level per y bin";
    if (m == "lnN") return "QCD template + lnN";
    if (m == "free") return "QCD free per region";
    return m;
}

// the multijet SR multiplier the extractor reports for (flavour, charge, group)
// (mtsplit: one multiplier per BIN since 2026-10-21b, no group multiplier ->
//  n/a here; the per-group multijet YIELDS below are the like-for-like quantity)
Val Multiplier(const FitData &f, int bn, const char *fl, const char *ch, int g)
{
    if (f.model == "mtsplit") return Val();
    const std::string base = std::string("qcd_abcd_mult_") + fl + "_" + ch;
    const bool perGroup = (f.model == "abcdy");
    const std::string key = perGroup ? base + "_" + pOAbcdY::kGroupName[g] : base;
    auto it = f.par[bn].find(key);
    return it == f.par[bn].end() ? Val() : it->second;
}

TGraphErrors *CloneGraph(TFile *f, const char *name, const TString &uniq)
{
    if (!f) return nullptr;
    TGraphErrors *g = dynamic_cast<TGraphErrors *>(f->Get(name));
    if (!g) return nullptr;
    g = (TGraphErrors *)g->Clone(Form("%s_%s", name, uniq.Data()));
    return g;
}

bool LoadFit(const char *tag, const char *fitTag, FitData &f)
{
    f = FitData();
    f.tag = tag;
    f.fitTag = fitTag;
    TString dsuf, dlab;
    if (!pODisc::Spec(tag, dsuf, dlab)) return false;
    if (!pOFit::Get(fitTag, f.spec)) return false;
    f.work = pOFit::WorkDir(dsuf, f.spec);
    const TString sumCsv = pOFit::SummaryFile(dsuf, f.spec, "summary.csv");
    std::ifstream in(sumCsv.Data());
    if (!in)
    {
        std::cerr << "[WARN] compare_fits: no " << sumCsv << " -> fit '" << fitTag << "' of " << tag << " skipped\n";
        return false;
    }
    std::string line;
    std::getline(in, line); // header
    while (std::getline(in, line))
    {
        const std::vector<std::string> t = Split(line);
        if (t.size() < 4) continue;
        for (int bn = 0; bn < 2; ++bn)
        {
            const std::string nom = std::string("simfit_") + kBn[bn];
            if (t[0] == nom)
            {
                if (t[1] == "fit_status") { f.status[bn] = std::atoi(t[2].c_str()); f.covq[bn] = std::atoi(t[3].c_str()); continue; }
                Val &v = f.par[bn][t[1]]; // NB keep .s: the _stat row may come first
                v.v = std::atof(t[2].c_str());
                v.e = std::atof(t[3].c_str());
                v.ok = true;
            }
            else if (t[0] == nom + "_stat")
                f.par[bn][t[1]].s = std::atof(t[3].c_str());
        }
    }
    // the Asimov closure verdict, from the extraction log
    const TString logName = (f.spec.tag == "comb") ? "extract_simfit.log" : TString("extract_") + f.spec.tag + ".log";
    std::ifstream lg((f.work + "/summary/" + logName).Data());
    while (lg && std::getline(lg, line))
        for (int bn = 0; bn < 2; ++bn)
        {
            const std::string key = std::string("[asimov] simfit_") + kBn[bn] + " closure";
            if (line.compare(0, key.size(), key) == 0)
            {
                const size_t p = line.find("=>");
                f.asimov[bn] = (p == std::string::npos) ? line : line.substr(p + 3);
            }
        }
    // the multijet model the cards were built with
    std::ifstream kp((f.work + "/datacards/qcd_lnn_kappas.txt").Data());
    while (kp && std::getline(kp, line))
        if (line.compare(0, 8, "qcdMode ") == 0) f.model = line.substr(8);
    // the SR multijet per group, prefit and postfit, from the fit's own shapes
    for (int bn = 0; bn < 2; ++bn)
    {
        const TString fd = f.work + "/fits/simfit_" + kBn[bn] + "/fitDiagnostics_simfit_" + kBn[bn] + ".root";
        TFile *fF = TFile::Open(fd, "READ");
        if (!fF || fF->IsZombie() || !fF->Get("shapes_fit_s") || !fF->Get("shapes_prefit"))
        {
            std::cerr << "[WARN] compare_fits: no shapes in " << fd << " -> no multijet-per-group numbers for "
                      << tag << " / " << fitTag << " / " << kBn[bn] << "\n";
            if (fF) { fF->Close(); delete fF; }
            continue;
        }
        // mtsplit (2026-10-21): a bin is the two channels <ch>_lo / <ch>_hi
        std::vector<TString> parts = {""};
        if (f.model == "mtsplit") parts = {"_lo", "_hi"};
        bool all = true;
        for (const TString &fl : f.spec.flavours)
            for (int ic = 0; ic < 2; ++ic)
                for (int iy = 0; iy < pOAbcdY::kNYBins; ++iy)
                    for (const TString &pt : parts)
                    {
                        const TString ch = TString::Format("%s_%s_%s_y%d", fl.Data(), kCh[ic], kBn[bn], iy) + pt;
                        TH1 *hPre = (TH1 *)fF->Get("shapes_prefit/" + ch + "/qcd");
                        TH1 *hPost = (TH1 *)fF->Get("shapes_fit_s/" + ch + "/qcd");
                        if (!hPre || !hPost) { all = false; continue; }
                        const int il = FlavIndex(fl), g = pOAbcdY::GroupOf(iy);
                        f.qPre[bn][il][ic][g] += hPre->Integral();
                        f.qPost[bn][il][ic][g] += hPost->Integral();
                        f.qPreTot[bn][il][ic] += hPre->Integral();
                    }
        if (!all)
            std::cerr << "[WARN] compare_fits: some SR qcd shapes missing in " << fd << "\n";
        f.haveShapes[bn] = all;
        fF->Close();
        delete fF;
    }
    // sigma W+/W-/W/Z per binning: the contour CSVs of run_observables.sh <tag>
    for (int bn = 0; bn < 2; ++bn)
    {
        const TString csv = (f.spec.tag == "comb")
                                ? TString::Format("./plots/comb/xsec/%s/xsec_contour_WZ_%s.csv", tag, kBn[bn])
                                : TString::Format("./plots/flavfit/xsec/%s/xsec_contour_WZ_%s_%s.csv", tag, kBn[bn], f.spec.flav.Data());
        std::ifstream sc(csv.Data());
        if (!sc) { std::cerr << "[WARN] compare_fits: no " << csv << " (run analysis/run_observables.sh " << tag << ")\n"; continue; }
        while (std::getline(sc, line))
        {
            const std::vector<std::string> t = Split(line);
            if (t.size() < 4) continue;
            for (int q = 0; q < 4; ++q)
                if (t[0] == kSigQ[q])
                {
                    f.sig[bn][q].v = std::atof(t[1].c_str());
                    f.sig[bn][q].e = std::atof(t[2].c_str());
                    f.sig[bn][q].s = std::atof(t[3].c_str());
                    f.sig[bn][q].ok = true;
                }
        }
    }
    // A_ch and R_FB: the FIDUCIAL graphs (r x sigma_gen)
    const TString uniq = TString(tag) + "_" + fitTag;
    TFile *fC = TFile::Open(pODisc::GraphFile("charge_asym", fitTag, tag), "READ");
    TFile *fB = TFile::Open(pODisc::GraphFile("FBratio", fitTag, tag), "READ");
    if (fC && !fC->IsZombie())
    {
        f.gA = CloneGraph(fC, "g_chargeAsym", uniq);
        f.gAs = CloneGraph(fC, "g_chargeAsym_stat", uniq);
    }
    if (fB && !fB->IsZombie())
        for (int k = 0; k < 3; ++k)
        {
            f.gR[k] = CloneGraph(fB, Form("g_RFB_%s", kRfb[k]), uniq);
            f.gRs[k] = CloneGraph(fB, Form("g_RFB_%s_stat", kRfb[k]), uniq);
        }
    if (!f.gA || !f.gR[0])
        std::cerr << "[WARN] compare_fits: no fiducial A_ch / R_FB graphs for " << tag << " / " << fitTag
                  << " (run analysis/run_observables.sh " << tag << ")\n";
    if (fC) { fC->Close(); delete fC; }
    if (fB) { fB->Close(); delete fB; }
    f.ok = true;
    return true;
}

// one comparison row: printed + written to the CSV
struct Out
{
    std::ofstream csv;
    void Row(const char *sec, const TString &fit, const char *bn, const std::string &q,
             const Val &a, const Val &b, bool print, const char *unit = "")
    {
        const double d = b.v - a.v;
        const double den = (a.s > 0) ? a.s : a.e;
        const double z = (den > 0) ? d / den : 0.0;
        csv << sec << "," << fit << "," << bn << "," << q << "," << a.v << "," << a.e << "," << a.s << ","
            << b.v << "," << b.e << "," << b.s << "," << d << "," << ((a.s > 0) ? z : 0.0) << "\n";
        if (!print) return;
        if (a.s > 0)
            printf("  %-26s ref %10.4f +- %.4f (stat %.4f)   alt %10.4f +- %.4f (stat %.4f)   d %+9.4f = %+5.2f stat%s\n",
                   q.c_str(), a.v, a.e, a.s, b.v, b.e, b.s, d, z, unit);
        else
            printf("  %-26s ref %10.4f +- %.4f   alt %10.4f +- %.4f   d %+9.4f%s\n",
                   q.c_str(), a.v, a.e, b.v, b.e, d, unit);
    }
};

// graph point -> Val (total and statistical error)
Val GraphVal(const TGraphErrors *g, const TGraphErrors *gs, int i)
{
    Val v;
    if (!g || i >= g->GetN()) return v;
    v.v = g->GetPointY(i);
    v.e = g->GetErrorY(i);
    v.s = (gs && i < gs->GetN()) ? gs->GetErrorY(i) : -1;
    v.ok = true;
    return v;
}

// one graph observable (A_ch or one R_FB) of one fit: every bin + a summary
void CompareGraph(Out &o, const char *sec, const TString &fit, TGraphErrors *ga, TGraphErrors *gas,
                  TGraphErrors *gb, TGraphErrors *gbs, const char *what, double &maxZ)
{
    maxZ = 0;
    if (!ga || !gb || ga->GetN() != gb->GetN())
    {
        printf("  %s: graphs missing or of different size -> skipped\n", what);
        return;
    }
    double sum2 = 0;
    int imax = -1;
    for (int i = 0; i < ga->GetN(); ++i)
    {
        const Val a = GraphVal(ga, gas, i), b = GraphVal(gb, gbs, i);
        o.Row(sec, fit, "", Form("%s_bin%d_eta%+.4f", what, i, ga->GetPointX(i)), a, b, true);
        if (a.s > 0)
        {
            const double z = (b.v - a.v) / a.s;
            sum2 += z * z;
            if (std::fabs(z) > std::fabs(maxZ)) { maxZ = z; imax = i; }
        }
    }
    printf("  -> %s: max |d|/stat = %.2f (bin %d, eta_CM %+.3f), sum (d/stat)^2 = %.2f over %d bins\n",
           what, std::fabs(maxZ), imax, imax >= 0 ? ga->GetPointX(imax) : 0.0, sum2, ga->GetN());
}

// ref/alt overlay of one graph observable
void OverlayPlot(TGraphErrors *ga, TGraphErrors *gas, TGraphErrors *gb, TGraphErrors *gbs,
                 const TString &labA, const TString &labB, const std::string &out,
                 const char *yTitle, const char *sub1, const TString &sub2, double ylo, double yhi, double yline,
                 const TString &note = "")
{
    if (!ga || !gb) return;
    PlotStyle ps;
    ps.showStats = false;
    ps.logy = false;
    ps.systBoxWidthFrac = 0.20;
    ps.systBoxFillAlpha = 0.30;
    std::vector<OverlaySeries> ser(2);
    ser[0].gTot = ga; ser[0].gStat = gas; ser[0].label = labA.Data();
    ser[0].color = kBlack; ser[0].marker = 20; ser[0].markerSize = 1.3; ser[0].xShift = -0.35;
    ser[1].gTot = gb; ser[1].gStat = gbs; ser[1].label = labB.Data();
    ser[1].color = kMagenta + 2; ser[1].marker = 21; ser[1].markerSize = 1.5; ser[1].xShift = 0.35;
    ps.headerX = 0.55; // room for the "<fit>, fiducial" line
    // optional 4th header line (e.g. the provenance warning), as observables.C's withDiscTag
    const double nx = ps.headerX, ny = ps.headerY - 3.0 * ps.headerDy;
    GraphTuner tune = [ylo, yhi, yline, note, nx, ny](TCanvas *c, TGraphErrors *g) {
        if (auto *h = g->GetHistogram()) { h->SetMinimum(ylo); h->SetMaximum(yhi); }
        const double x1 = g->GetXaxis()->GetXmin(), x2 = g->GetXaxis()->GetXmax();
        TLine *l = new TLine(x1, yline, x2, yline);
        l->SetLineStyle(2);
        l->Draw("same");
        if (note.Length())
        {
            TLatex lat;
            lat.SetNDC(true);
            lat.SetTextFont(42);
            lat.SetTextAlign(13);
            lat.SetTextSize(0.028);
            lat.SetTextColor(kRed + 1);
            lat.DrawLatex(nx, ny, note.Data());
        }
        c->Modified();
        c->Update();
    };
    SaveNiceGraph_Overlay(ser, out, "#eta^{l}_{CM}", yTitle, "", sub1, sub2.Data(), {}, ps, tune);
}

// the multijet-per-group plot of one fit and binning: every category divided
// by the REF prefit of that group -> ref postfit, alt prefit, alt postfit
void GroupPlot(const FitData &A, const FitData &B, int bn, const TString &labA, const TString &labB,
               const std::string &out, const TString &sub2)
{
    if (!A.haveShapes[bn] || !B.haveShapes[bn]) return;
    std::vector<std::string> catLab;
    std::vector<std::string> blockLab;
    std::vector<double> yA, eA, yBp, yB, eB;
    for (const TString &fl : A.spec.flavours)
        for (int ic = 0; ic < 2; ++ic)
        {
            const int il = FlavIndex(fl);
            blockLab.push_back(std::string(il == 0 ? "#mu" : "e") + (ic == 0 ? "^{+}" : "^{-}"));
            for (int k = 0; k < 3; ++k)
            {
                const int g = kGOrd[k];
                const double ref = A.qPre[bn][il][ic][g];
                const Val mA = Multiplier(A, bn, fl.Data(), kCh[ic], g), mB = Multiplier(B, bn, fl.Data(), kCh[ic], g);
                const double pa = A.qPost[bn][il][ic][g] / ref, pb = B.qPost[bn][il][ic][g] / ref;
                catLab.push_back(pOAbcdY::kGroupName[g]);
                yA.push_back(pa);
                eA.push_back((mA.ok && mA.v > 0) ? pa * mA.e / mA.v : 0.0);
                yBp.push_back(B.qPre[bn][il][ic][g] / ref);
                yB.push_back(pb);
                eB.push_back((mB.ok && mB.v > 0) ? pb * mB.e / mB.v : 0.0);
            }
        }
    const int n = (int)catLab.size();
    TGraphErrors *gA = new TGraphErrors(n), *gBp = new TGraphErrors(n), *gB = new TGraphErrors(n);
    gA->SetName(Form("g_qgrp_ref_%s_%s", A.fitTag.Data(), kBn[bn]));
    gBp->SetName(Form("g_qgrp_altpre_%s_%s", A.fitTag.Data(), kBn[bn]));
    gB->SetName(Form("g_qgrp_alt_%s_%s", A.fitTag.Data(), kBn[bn]));
    double lo = 1e9, hi = -1e9;
    for (int i = 0; i < n; ++i)
    {
        gA->SetPoint(i, i, yA[i]);   gA->SetPointError(i, 0.5, eA[i]);
        gBp->SetPoint(i, i, yBp[i]); gBp->SetPointError(i, 0.5, 0.0);
        gB->SetPoint(i, i, yB[i]);   gB->SetPointError(i, 0.5, eB[i]);
        lo = std::min({lo, yA[i] - eA[i], yBp[i], yB[i] - eB[i]});
        hi = std::max({hi, yA[i] + eA[i], yBp[i], yB[i] + eB[i]});
    }
    const double ylo = std::min(0.2, lo - 0.30), yhi = hi + 0.55;
    std::vector<OverlaySeries> ser(3);
    ser[0].gTot = gA;  ser[0].label = (labA + ", postfit").Data();
    ser[0].color = kBlack; ser[0].marker = 20; ser[0].markerSize = 1.3; ser[0].xShift = -0.6;
    ser[1].gTot = gBp; ser[1].label = (labB + ", prefit").Data();
    ser[1].color = kMagenta + 2; ser[1].marker = 25; ser[1].markerSize = 1.5; ser[1].xShift = 0.0;
    ser[2].gTot = gB;  ser[2].label = (labB + ", postfit").Data();
    ser[2].color = kMagenta + 2; ser[2].marker = 21; ser[2].markerSize = 1.5; ser[2].xShift = 0.6;
    PlotStyle ps;
    ps.showStats = false;
    ps.logy = false;
    ps.headerX = 0.55; // the header lines are longer than the usual "W -> l nu"
    GraphTuner tune = [=](TCanvas *c, TGraphErrors *g) {
        if (auto *h = g->GetHistogram()) { h->SetMinimum(ylo); h->SetMaximum(yhi); }
        g->GetXaxis()->SetLabelSize(0);
        g->GetXaxis()->SetTickLength(0);
        c->Modified();
        c->Update();
        const double ux1 = gPad->GetUxmin(), ux2 = gPad->GetUxmax();
        const double L = gPad->GetLeftMargin(), R = gPad->GetRightMargin(), Bm = gPad->GetBottomMargin();
        auto xndc = [&](double x) { return L + (x - ux1) / (ux2 - ux1) * (1.0 - L - R); };
        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextFont(42);
        lat.SetTextAlign(23);
        lat.SetTextSize(0.026);
        for (int i = 0; i < n; ++i) lat.DrawLatex(xndc(i), Bm - 0.012, catLab[i].c_str());
        lat.SetTextSize(0.036);
        for (size_t b = 0; b < blockLab.size(); ++b) lat.DrawLatex(xndc(3.0 * b + 1.0), Bm - 0.048, blockLab[b].c_str());
        for (size_t b = 1; b < blockLab.size(); ++b)
        {
            TLine *s = new TLine(3.0 * b - 0.5, ylo, 3.0 * b - 0.5, yhi);
            s->SetLineStyle(3);
            s->SetLineColor(kGray + 2);
            s->Draw("same");
        }
        TLine *l1 = new TLine(ux1, 1.0, ux2, 1.0);
        l1->SetLineStyle(2);
        l1->Draw("same");
        c->Modified();
        c->Update();
    };
    SaveNiceGraph_Overlay(ser, out, "", "SR multijet / ref. prefit", "", "SR multijet per y group",
                          sub2.Data(), {}, ps, tune);
}

} // namespace

void compare_fits(const char *refTag = "leppt_mt40", const char *altTag = "leppt_mt40_abcdy")
{
    gStyle->SetEndErrorSize(4);
    TString dsuf, refLab, altLab;
    if (!pODisc::Spec(refTag, dsuf, refLab) || !pODisc::Spec(altTag, dsuf, altLab)) return;
    if (pODisc::Base(refTag) != pODisc::Base(altTag))
    {
        // leppt_mt40 and met_mt40 (2026-10-20) are the SAME events (pT > 25 &&
        // m_T > 40), the same in-fit ABCD regions and the same gen fiducial:
        // one data set, two discriminants -- every row compares directly
        TString f1, l1, f2, l2;
        const bool sameFid = pODisc::GenFid(refTag, f1, l1) && pODisc::GenFid(altTag, f2, l2) && f1 == f2;
        const bool mt40Pair = pODisc::Base(refTag).EndsWith("_mt40") && pODisc::Base(altTag).EndsWith("_mt40");
        if (sameFid && mt40Pair)
            std::cout << "[NOTE] compare_fits: " << refTag << " and " << altTag << " fit DIFFERENT discriminants of the"
                      << " SAME events (pT > 25 && m_T > 40) in the same fiducial volume\n";
        else
            std::cerr << "[WARN] compare_fits: " << refTag << " and " << altTag
                      << " are fits of DIFFERENT discriminants -- the comparison is not of one data set\n";
    }

    const char *fits[3] = {"comb", "simfit_mu", "simfit_ele"};
    FitData F[2][3]; // [ref, alt][fit]
    for (int k = 0; k < 3; ++k)
    {
        LoadFit(refTag, fits[k], F[0][k]);
        LoadFit(altTag, fits[k], F[1][k]);
    }

    const std::string outDir = std::string("./plots/compare_fits/") + refTag + "_vs_" + altTag;
    gSystem->mkdir(outDir.c_str(), kTRUE);
    Out o;
    o.csv.open((outDir + "/compare_fits.csv").c_str());
    o.csv << "section,fit,binning,quantity,ref,ref_err,ref_stat,alt,alt_err,alt_stat,delta,delta_over_ref_stat\n";

    // ---- 1. what is compared ----------------------------------------------
    printf("==================================================================\n");
    printf("[compare_fits] ref = %s, alt = %s\n", refTag, altTag);
    printf("==================================================================\n");
    for (int k = 0; k < 3; ++k)
        for (int m = 0; m < 2; ++m)
        {
            const FitData &f = F[m][k];
            if (!f.ok) { printf("[fit] %-4s %-10s : ABSENT\n", m ? "alt" : "ref", fits[k]); continue; }
            printf("[fit] %-4s %-10s : %s  qcdMode %-6s (%s)\n", m ? "alt" : "ref", fits[k], f.work.Data(),
                   f.model.Data(), ModelLabel(f.model).Data());
            for (int bn = 0; bn < 2; ++bn)
                printf("        %-3s status %d covQual %d  Asimov closure: %s  r_Z %.4f +- %.4f\n", kBn[bn],
                       f.status[bn], f.covq[bn], f.asimov[bn].c_str(), f.par[bn].count("r_Z") ? f.par[bn].at("r_Z").v : 0.0,
                       f.par[bn].count("r_Z") ? f.par[bn].at("r_Z").e : 0.0);
        }
    const TString labA = ModelLabel(F[0][0].model), labB = ModelLabel(F[1][0].model);

    // provenance of the fb multijet templates (see the header)
    printf("\n[provenance] fb / lab prefit multijet total per (flavour, charge); 1.000 = fb templates carry the\n"
           "             full y-inclusive multijet (pre-2026-10-19 FB-window fix):\n");
    bool refPreFix = false;
    for (int m = 0; m < 2; ++m)
        for (int k = 0; k < 3; ++k)
        {
            const FitData &f = F[m][k];
            if (!f.ok || !f.haveShapes[0] || !f.haveShapes[1]) continue;
            printf("  %-4s %-10s", m ? "alt" : "ref", fits[k]);
            bool all1 = true;
            for (const TString &fl : f.spec.flavours)
                for (int ic = 0; ic < 2; ++ic)
                {
                    const int il = FlavIndex(fl);
                    const double r = f.qPreTot[1][il][ic] / f.qPreTot[0][il][ic];
                    printf("  %s_%s %.4f", fl.Data(), kCh[ic], r);
                    all1 &= std::fabs(r - 1.0) < 1e-3;
                    o.csv << "provenance," << fits[k] << ",fb," << (m ? "alt_" : "ref_") << fl << "_" << kCh[ic]
                          << "_fb_over_lab_prefit_qcd," << r << ",,,,,,,\n";
                }
            printf("%s\n", all1 ? "   <- PRE-FIX fb templates" : "");
            if (m == 0 && all1) refPreFix = true;
        }
    if (refPreFix)
        printf("[WARN] the REF fits used the pre-fix fb multijet templates: every fb comparison below mixes the\n"
               "       2026-10-19 FB-window fix with the model change (refit the ref model on the current inputs\n"
               "       to separate them). The lab comparison is clean (lab inputs unchanged by the fix).\n");

    // ---- 2. the multijet per rapidity group -------------------------------
    for (int k = 0; k < 3; ++k)
    {
        const FitData &A = F[0][k], &B = F[1][k];
        if (!A.ok || !B.ok) continue;
        for (int bn = 0; bn < 2; ++bn)
        {
            if (!A.haveShapes[bn] || !B.haveShapes[bn]) continue;
            printf("\n[qcd-groups] %s %s: SR multijet per group (events), ref = %s, alt = %s\n", fits[k], kBn[bn],
                   A.model.Data(), B.model.Data());
            printf("  %-12s %9s %9s %8s | %9s %9s %16s | %10s\n", "", "ref pre", "ref post", "ref mult", "alt pre",
                   "alt post", "alt mult", "post alt/ref");
            for (const TString &fl : A.spec.flavours)
                for (int ic = 0; ic < 2; ++ic)
                {
                    const int il = FlavIndex(fl);
                    double tA[2] = {0, 0}, tB[2] = {0, 0};
                    for (int kk = 0; kk < 3; ++kk)
                    {
                        const int g = kGOrd[kk];
                        const Val mA = Multiplier(A, bn, fl.Data(), kCh[ic], g), mB = Multiplier(B, bn, fl.Data(), kCh[ic], g);
                        const double aPre = A.qPre[bn][il][ic][g], aPost = A.qPost[bn][il][ic][g];
                        const double bPre = B.qPre[bn][il][ic][g], bPost = B.qPost[bn][il][ic][g];
                        tA[0] += aPre; tA[1] += aPost; tB[0] += bPre; tB[1] += bPost;
                        const std::string nm = std::string(fl.Data()) + "_" + kCh[ic] + "_" + pOAbcdY::kGroupName[g];
                        printf("  %-12s %9.1f %9.1f %8.3f | %9.1f %9.1f %8.3f +- %5.3f | %10.3f\n", nm.c_str(), aPre,
                               aPost, mA.v, bPre, bPost, mB.v, mB.e, aPost > 0 ? bPost / aPost : 0.0);
                        Val vA, vB;
                        vA.v = aPre; vB.v = bPre;
                        o.Row("qcd_group", fits[k], kBn[bn], nm + "_prefit", vA, vB, false);
                        vA.v = aPost; vA.e = (mA.ok && mA.v > 0) ? aPost * mA.e / mA.v : 0;
                        vB.v = bPost; vB.e = (mB.ok && mB.v > 0) ? bPost * mB.e / mB.v : 0;
                        o.Row("qcd_group", fits[k], kBn[bn], nm + "_postfit", vA, vB, false);
                        o.Row("qcd_group", fits[k], kBn[bn], nm + "_multiplier", mA, mB, false);
                    }
                    const std::string nm = std::string(fl.Data()) + "_" + kCh[ic] + "_all";
                    printf("  %-12s %9.1f %9.1f %8s | %9.1f %9.1f %16s | %10.3f\n", nm.c_str(), tA[0], tA[1], "", tB[0],
                           tB[1], "", tA[1] > 0 ? tB[1] / tA[1] : 0.0);
                    Val vA, vB;
                    vA.v = tA[0]; vB.v = tB[0];
                    o.Row("qcd_group", fits[k], kBn[bn], nm + "_prefit", vA, vB, false);
                    vA.v = tA[1]; vB.v = tB[1];
                    o.Row("qcd_group", fits[k], kBn[bn], nm + "_postfit", vA, vB, false);
                    if (k == 0 && bn == 0 && fl == "ele")
                        printf("[HEAD] comb lab e %s multijet postfit %.1f -> %.1f (x%.3f)\n", kCh[ic], tA[1], tB[1],
                               tA[1] > 0 ? tB[1] / tA[1] : 0.0);
                }
            GroupPlot(A, B, bn, labA, labB, outDir + "/qcd_groups_" + fits[k] + "_" + kBn[bn],
                      TString::Format("%s fit, %s", A.spec.label.Data(), kBn[bn]));
        }
    }

    // ---- 3. nuisance pulls ------------------------------------------------
    for (int k = 0; k < 3; ++k)
    {
        const FitData &A = F[0][k], &B = F[1][k];
        if (!A.ok || !B.ok) continue;
        for (int bn = 0; bn < 2; ++bn)
        {
            printf("\n[pulls] %s %s: theta (value +- post-fit constraint)\n", fits[k], kBn[bn]);
            std::vector<std::string> names;
            for (const auto &kv : A.par[bn]) if (kv.first.find("_theta") != std::string::npos) names.push_back(kv.first);
            for (const auto &kv : B.par[bn])
                if (kv.first.find("_theta") != std::string::npos && !A.par[bn].count(kv.first)) names.push_back(kv.first);
            for (const std::string &nm : names)
            {
                const Val a = A.par[bn].count(nm) ? A.par[bn].at(nm) : Val();
                const Val b = B.par[bn].count(nm) ? B.par[bn].at(nm) : Val();
                Val a2 = a, b2 = b;
                a2.s = -1; b2.s = -1; // a pull is already in prior-sigma units
                o.Row("pull", fits[k], kBn[bn], nm, a2, b2, true);
            }
        }
    }

    // ---- 4. r per bin and r_Z ---------------------------------------------
    for (int k = 0; k < 3; ++k)
    {
        const FitData &A = F[0][k], &B = F[1][k];
        if (!A.ok || !B.ok) continue;
        for (int bn = 0; bn < 2; ++bn)
        {
            printf("\n[r] %s %s: signal strengths (total, stat)\n", fits[k], kBn[bn]);
            double maxZ = 0, sumZ = 0;
            std::string where;
            for (int ic = 0; ic < 2; ++ic)
                for (int iy = 0; iy < 12; ++iy)
                {
                    const std::string nm = Form("r_%s_y%d", kCh[ic], iy);
                    const Val a = A.par[bn].count(nm) ? A.par[bn].at(nm) : Val();
                    const Val b = B.par[bn].count(nm) ? B.par[bn].at(nm) : Val();
                    o.Row("r", fits[k], kBn[bn], nm, a, b, true);
                    if (a.s > 0)
                    {
                        const double z = (b.v - a.v) / a.s;
                        sumZ += z;
                        if (std::fabs(z) > std::fabs(maxZ)) { maxZ = z; where = nm; }
                    }
                }
            const Val a = A.par[bn].count("r_Z") ? A.par[bn].at("r_Z") : Val();
            const Val b = B.par[bn].count("r_Z") ? B.par[bn].at("r_Z") : Val();
            o.Row("r", fits[k], kBn[bn], "r_Z", a, b, true);
            printf("  -> 24 r: mean d/stat %+.3f, max |d|/stat %.2f (%s)\n", sumZ / 24.0, std::fabs(maxZ), where.c_str());
        }
    }

    // ---- 5. cross sections, e/mu ratios -----------------------------------
    for (int bn = 0; bn < 2; ++bn)
    {
        printf("\n[sigma] %s binning: fiducial sigma = r x sigma_gen (nb)\n", kBn[bn]);
        for (int k = 0; k < 3; ++k)
        {
            const FitData &A = F[0][k], &B = F[1][k];
            if (!A.ok || !B.ok) continue;
            printf(" %s\n", fits[k]);
            for (int q = 0; q < 4; ++q)
            {
                o.Row("sigma", fits[k], kBn[bn], kSigQ[q], A.sig[bn][q], B.sig[bn][q], true);
                if (k == 0 && (q == 2 || q == 3))
                    printf("[HEAD] comb %s %s %.3f -> %.3f nb (d %+.3f = %+.2f stat)\n", kBn[bn], kSigQ[q],
                           A.sig[bn][q].v, B.sig[bn][q].v, B.sig[bn][q].v - A.sig[bn][q].v,
                           A.sig[bn][q].s > 0 ? (B.sig[bn][q].v - A.sig[bn][q].v) / A.sig[bn][q].s : 0.0);
            }
        }
        // e/mu from the per-flavour fits (independent samples -> stat errors in quadrature)
        const FitData *mu[2] = {&F[0][1], &F[1][1]}, *el[2] = {&F[0][2], &F[1][2]};
        if (!mu[0]->ok || !mu[1]->ok || !el[0]->ok || !el[1]->ok) continue;
        printf(" e/mu (per-flavour fits, stat errors)\n");
        Val R[2][4];
        for (int m = 0; m < 2; ++m)
            for (int q = 0; q < 4; ++q)
            {
                const Val &a = mu[m]->sig[bn][q], &b = el[m]->sig[bn][q];
                if (!a.ok || !b.ok || a.v <= 0) continue;
                R[m][q].v = b.v / a.v;
                R[m][q].s = R[m][q].v * std::sqrt(std::pow(a.s / a.v, 2) + std::pow(b.s / b.v, 2));
                R[m][q].e = R[m][q].s;
                R[m][q].ok = true;
            }
        for (int q = 0; q < 4; ++q)
        {
            o.Row("e_over_mu", "flavfit", kBn[bn], kSigQ[q], R[0][q], R[1][q], true);
            if (q == 1 || q == 3)
                printf("[HEAD] %s e/mu %s %.4f +- %.4f -> %.4f +- %.4f (stat)\n", kBn[bn], kSigQ[q], R[0][q].v,
                       R[0][q].s, R[1][q].v, R[1][q].s);
        }
        Val D[2];
        for (int m = 0; m < 2; ++m)
            if (R[m][0].ok && R[m][1].ok)
            {
                D[m].v = R[m][0].v / R[m][1].v;
                D[m].s = D[m].v * std::sqrt(std::pow(R[m][0].s / R[m][0].v, 2) + std::pow(R[m][1].s / R[m][1].v, 2));
                D[m].e = D[m].s;
                D[m].ok = true;
            }
        o.Row("e_over_mu", "flavfit", kBn[bn], "double_ratio_Wp_over_Wm", D[0], D[1], true);
        printf("[HEAD] %s (e/mu)_W+ / (e/mu)_W- %.4f +- %.4f -> %.4f +- %.4f (stat)\n", kBn[bn], D[0].v, D[0].s, D[1].v,
               D[1].s);
    }

    // ---- 6. A_ch and R_FB ------------------------------------------------
    for (int k = 0; k < 3; ++k)
    {
        const FitData &A = F[0][k], &B = F[1][k];
        if (!A.ok || !B.ok) continue;
        double z = 0;
        printf("\n[A_ch] %s (fiducial, lab bins)\n", fits[k]);
        CompareGraph(o, "A_ch", fits[k], A.gA, A.gAs, B.gA, B.gAs, "A_ch", z);
        printf("[HEAD] %s A_ch max |d|/stat %.2f\n", fits[k], std::fabs(z));
        for (int r = 0; r < 3; ++r)
        {
            printf("\n[R_FB] %s %s (fiducial, |eta_CM| bins)\n", fits[k], kRfb[r]);
            CompareGraph(o, "R_FB", fits[k], A.gR[r], A.gRs[r], B.gR[r], B.gRs[r], Form("R_FB_%s", kRfb[r]), z);
            printf("[HEAD] %s R_FB_%s max |d|/stat %.2f%s\n", fits[k], kRfb[r], std::fabs(z),
                   refPreFix ? "  (ref fb inputs pre-fix: mixes the FB-window fix)" : "");
        }
        // the overlays
        const char *lepOf = (A.spec.flav == "mu") ? "W #rightarrow #mu #nu" : (A.spec.flav == "ele") ? "W #rightarrow e #nu"
                                                                                                    : "W #rightarrow l #nu";
        const TString sub2 = TString::Format("%s, fiducial", A.spec.label.Data());
        OverlayPlot(A.gA, A.gAs, B.gA, B.gAs, labA, labB, outDir + "/chargeAsym_" + fits[k], "A_{ch}", lepOf, sub2,
                    -0.6, 0.6, 0.0);
        // R_FB comes from the fb fit: say so on the plot when the ref fb fit predates the FB-window fix
        const TString note = refPreFix ? "ref. fb fit: pre-FB-window-fix inputs" : "";
        for (int r = 0; r < 3; ++r)
            OverlayPlot(A.gR[r], A.gRs[r], B.gR[r], B.gRs[r], labA, labB,
                        outDir + "/RFB_" + kRfb[r] + "_" + fits[k], "R_{FB}",
                        r == 0 ? lepOf : (r == 1 ? "W^{+} #rightarrow l^{+} #nu" : "W^{-} #rightarrow l^{-} #bar{#nu}"),
                        sub2, 0.3, 2.1, 1.0, note);
    }

    o.csv.close();
    printf("\n[DONE] compare_fits: %s vs %s -> %s/{compare_fits.csv, qcd_groups_*, chargeAsym_*, RFB_*}\n", refTag,
           altTag, outDir.c_str());
}
