// =============================================================================
// plotting/r_per_bin.C -- the POSTFIT signal strength per rapidity bin, mu vs e
// (2026-10-20, user request: "a per rapidity bin r plot, based on the post-fit
// r, from which we can see why the ABCD cannot help").
//
//   root -l -b -q 'r_per_bin.C+("leppt_mt40_abcdy","leppt_mt40")'
//   (run by analysis/run_observables.sh for every tag; the second argument is
//    given for a fit-model variant tag = its base tag, see below)
//
// r_i = the fitted POI r_<C>_y<i> of each simultaneous fit (comb = mu + e with
// shared r, simfit_mu, simfit_ele), read from the fit's <tag>_W_yields.csv:
// value r, total error rErr, statistical error rErr_stat (the --statonly
// companion fit). Drawn as point + stat bar + syst box (syst = sqrt(tot^2 -
// stat^2)). These are the fitted strengths themselves -- not r x sigma_gen --
// i.e. data / absolute POWHEG MC per bin: for the per-flavour fits the mu/e
// difference of r IS the mu/e difference of the measured cross section.
//
// Plots, plots/flavfit/r/<disc>/:
//   r_{Wp,Wm}_{lab,fb}        mu-only | mu+e | e-only per bin
//   r_emu_{lab,fb}            r_e / r_mu per bin, W+ and W- (STATISTICAL errors:
//                             the two fits are independent, lumi cancels)
//   qcd_needed_{Wp,Wm}_{lab,fb}   THE "WHY" PLOT: per multijet group of the fit
//                             (abcdy: fwd y10-11 | bwd y0-1 | cen y2-9; any other
//                             model: one group = all 12 bins), the factor on the
//                             e-only POSTFIT multijet that would make r_e = r_mu,
//                               f_g = 1 - Sum_{i in g} (r_mu,i - r_e,i) S_e,i / Q_g
//                             (S_e,i = postfit W / r_e,i = the e W template the r
//                             scales; Q_g = postfit e multijet of the group; data
//                             and every other process held at their postfit
//                             values), against the grey band 1 +- the fitted
//                             group multiplier's relative error = how far the
//                             multijet control regions let that group move.
//                             Error on f_g: r_mu stat (+) the e DATA Poisson of
//                             the group's bins, nothing from the multijet itself
//                             (that is the band) -- bins treated as independent.
// With a reference tag (the same data fitted with another multijet model; for
// leppt_mt40_abcdy the reference is leppt_mt40 = one ABCD per charge):
//   r_model_{Wp,Wm}_{lab,fb}  mu and e, reference (open) vs this model (filled)
//   r_emu_model_{lab,fb}      r_e / r_mu, W+ and W-, reference vs this model
//   + the e postfit multijet per group, this model / reference, in the box
// r_per_bin.csv: every number of every plot. The console carries the per-bin
// tables (mu vs e pulls, the per-bin and per-group "needed" numbers) and goes
// to analysis/logs/run_observables_<tag>.log.
//
// PROVENANCE: a reference fit whose fb prefit multijet equals its lab one
// (ratio 1.000 per flavour and charge) predates the 2026-10-19 FB-window fix
// (see compare_fits.C); the fb model plots then say so and the console WARNs.
// =============================================================================
#include "TBox.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLine.h"
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

#include "plotting_helper.C"              // PlotStyle, SaveNiceGraph_Overlay, CMS_lumi
#include "disc_variants.h"                // pODisc::Spec / Base
#include "fit_variants.h"                 // pOFit::Spec / WorkDir / SummaryFile (honours $FORK_TEST)
#include "../correction/abcd_ygroups.h"   // pOAbcdY: the rapidity groups of QCD_MODE=abcdy
#include "../skim/skim_common.h"          // pOSkim::kNY, kYEdges, kYEdgesFB
#include "../analysis/analysis_helpers.h" // pOAnalysis::kPORapidityShift

namespace
{
const int kNYr = pOSkim::kNY;
const char *const kBnR[2] = {"lab", "fb"};
const char *const kChR[2] = {"Wp", "Wm"};
const char *const kChLab[2] = {"W^{+} #rightarrow l^{+} #nu", "W^{-} #rightarrow l^{-} #bar{#nu}"};
const char *const kFitTag[3] = {"simfit_mu", "comb", "simfit_ele"};

struct RVal
{
    double v = NAN, e = NAN, s = NAN; // value, total error, statistical error
    bool ok = false;
    double syst() const { return std::sqrt(std::max(0.0, e * e - s * s)); }
};

struct FitR
{
    bool ok = false;
    TString disc, model, work;
    pOFit::Spec spec;
    RVal r[2][2][kNYr];                     // [binning][charge][y]
    std::map<std::string, RVal> par[2];     // summary.csv rows of simfit_<bn>
    // postfit shapes of the fit's OWN flavour (per-flavour fits only)
    bool haveShapes[2] = {false, false};
    double W[2][2][kNYr] = {}, Q[2][2][kNYr] = {}, D[2][2][kNYr] = {};
    double qPreTot[2][2] = {};              // prefit multijet total per (binning, charge)
};

std::vector<std::string> SplitCsv(const std::string &line)
{
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, ',')) out.push_back(tok);
    return out;
}

TString ModelName(const TString &m)
{
    if (m == "abcd") return "in-fit ABCD per charge";
    if (m == "abcdy") return "in-fit ABCD per y group";
    if (m == "mtsplit") return "m_{T} split, multijet level per y bin";
    if (m == "lnN") return "QCD template + lnN";
    if (m == "free") return "QCD free per region";
    return m;
}

double XCenter(int bn, int y)
{
    if (bn == 0) return 0.5 * (pOSkim::kYEdges[y] + pOSkim::kYEdges[y + 1]);
    return 0.5 * (pOSkim::kYEdgesFB[y] + pOSkim::kYEdgesFB[y + 1]) - pOAnalysis::kPORapidityShift;
}
double XHalf(int bn, int y)
{
    const double *e = (bn == 0) ? pOSkim::kYEdges : pOSkim::kYEdgesFB;
    return 0.5 * (e[y + 1] - e[y]);
}

bool LoadFitR(const char *disc, const char *fitTag, FitR &f)
{
    f = FitR();
    f.disc = disc;
    TString dsuf, dlab;
    if (!pODisc::Spec(disc, dsuf, dlab)) return false;
    if (!pOFit::Get(fitTag, f.spec)) return false;
    f.work = pOFit::WorkDir(dsuf, f.spec);

    // ---- r per bin: <tag>_W_yields.csv, columns found by header name
    const TString yCsv = pOFit::SummaryFile(dsuf, f.spec, "W_yields.csv");
    std::ifstream in(yCsv.Data());
    if (!in)
    {
        std::cerr << "[WARN] r_per_bin: no " << yCsv << " -> fit '" << fitTag << "' of " << disc << " skipped\n";
        return false;
    }
    std::string line;
    std::getline(in, line);
    const std::vector<std::string> hdr = SplitCsv(line);
    auto col = [&](const char *name) -> int
    {
        for (size_t k = 0; k < hdr.size(); ++k)
            if (hdr[k] == name) return (int)k;
        return -1;
    };
    const int cCh = col("charge"), cBn = col("binning"), cY = col("ybin"), cR = col("r"), cE = col("rErr"),
              cS = col("rErr_stat"), cM = col("qcd_model");
    if (cCh < 0 || cBn < 0 || cY < 0 || cR < 0 || cE < 0)
    {
        std::cerr << "[ERROR] r_per_bin: unexpected header in " << yCsv << "\n";
        return false;
    }
    if (cS < 0)
        std::cerr << "[WARN] r_per_bin: no rErr_stat column in " << yCsv << " -> total error bars only\n";
    int nRead = 0;
    while (std::getline(in, line))
    {
        const std::vector<std::string> t = SplitCsv(line);
        if ((int)t.size() <= std::max({cCh, cBn, cY, cR, cE})) continue;
        const int ic = (t[cCh] == "Wp") ? 0 : (t[cCh] == "Wm" ? 1 : -1);
        const int bn = (t[cBn] == "lab") ? 0 : (t[cBn] == "fb" ? 1 : -1);
        const int iy = std::atoi(t[cY].c_str());
        if (ic < 0 || bn < 0 || iy < 0 || iy >= kNYr) continue;
        RVal &v = f.r[bn][ic][iy];
        v.v = std::atof(t[cR].c_str());
        v.e = std::atof(t[cE].c_str());
        v.s = (cS >= 0 && cS < (int)t.size()) ? std::atof(t[cS].c_str()) : v.e;
        v.ok = std::isfinite(v.v) && v.e > 0;
        if (cM >= 0 && cM < (int)t.size() && f.model.IsNull()) f.model = t[cM];
        ++nRead;
    }
    if (nRead < 2 * 2 * kNYr)
        std::cerr << "[WARN] r_per_bin: only " << nRead << " of 48 r rows in " << yCsv << "\n";

    // ---- the multijet multipliers: <tag>_summary.csv rows simfit_<bn>
    std::ifstream sc(pOFit::SummaryFile(dsuf, f.spec, "summary.csv").Data());
    while (sc && std::getline(sc, line))
    {
        const std::vector<std::string> t = SplitCsv(line);
        if (t.size() < 4) continue;
        for (int bn = 0; bn < 2; ++bn)
            if (t[0] == std::string("simfit_") + kBnR[bn])
            {
                RVal &v = f.par[bn][t[1]];
                v.v = std::atof(t[2].c_str());
                v.e = std::atof(t[3].c_str());
                v.ok = true;
            }
    }
    // the model as the cards were built (the CSV column says the same)
    std::ifstream kp((f.work + "/datacards/qcd_lnn_kappas.txt").Data());
    while (kp && std::getline(kp, line))
        if (line.compare(0, 8, "qcdMode ") == 0) f.model = line.substr(8);

    // ---- postfit shapes of a per-flavour fit (its own flavour's SR channels)
    if (!f.spec.flav.IsNull())
        for (int bn = 0; bn < 2; ++bn)
        {
            const TString fd = f.work + "/fits/simfit_" + kBnR[bn] + "/fitDiagnostics_simfit_" + kBnR[bn] + ".root";
            TFile *fF = TFile::Open(fd, "READ");
            if (!fF || fF->IsZombie() || !fF->Get("shapes_fit_s") || !fF->Get("shapes_prefit"))
            {
                std::cerr << "[WARN] r_per_bin: no shapes in " << fd << "\n";
                if (fF) { fF->Close(); delete fF; }
                continue;
            }
            // mtsplit (2026-10-21): each bin is two channels, <ch>_lo (m_T < 30)
            // and <ch>_hi, sharing r -- the bin's W / multijet / data = their sum
            std::vector<TString> parts = {""};
            if (f.model == "mtsplit") parts = {"_lo", "_hi"};
            bool all = true;
            for (int ic = 0; ic < 2; ++ic)
                for (int iy = 0; iy < kNYr; ++iy)
                    for (const TString &pt : parts)
                    {
                        const TString ch = TString::Format("%s_%s_%s_y%d", f.spec.flav.Data(), kChR[ic], kBnR[bn], iy) + pt;
                        TH1 *hS = (TH1 *)fF->Get("shapes_fit_s/" + ch + "/signal");
                        TH1 *hT = (TH1 *)fF->Get("shapes_fit_s/" + ch + "/wtau");
                        TH1 *hQ = (TH1 *)fF->Get("shapes_fit_s/" + ch + "/qcd");
                        TH1 *hQp = (TH1 *)fF->Get("shapes_prefit/" + ch + "/qcd");
                        TGraphAsymmErrors *gD = (TGraphAsymmErrors *)fF->Get("shapes_fit_s/" + ch + "/data");
                        if (!hS || !hQ || !hQp || !gD) { all = false; continue; }
                        f.W[bn][ic][iy] += hS->Integral() + (hT ? hT->Integral() : 0.0);
                        f.Q[bn][ic][iy] += hQ->Integral();
                        double d = 0;
                        for (int k = 0; k < gD->GetN(); ++k) d += gD->GetPointY(k);
                        f.D[bn][ic][iy] += d;
                        f.qPreTot[bn][ic] += hQp->Integral();
                    }
            if (!all) std::cerr << "[WARN] r_per_bin: some SR shapes missing in " << fd << "\n";
            f.haveShapes[bn] = all;
            fF->Close();
            delete fF;
        }
    f.ok = nRead > 0;
    return f.ok;
}

// the multijet "groups" of a fit model = the units with their own fitted
// multijet level: abcdy -> the three y groups; mtsplit (per-bin level since
// 2026-10-21b) -> every rapidity bin; else one
bool PerGroup(const TString &model) { return model == "abcdy"; }
bool PerBin(const TString &model) { return model == "mtsplit"; }
std::vector<std::vector<int>> GroupsOf(const TString &model)
{
    std::vector<std::vector<int>> g;
    if (PerGroup(model))
        for (int k = 0; k < pOAbcdY::kNGroups; ++k) g.push_back(pOAbcdY::BinsOf(k));
    else if (PerBin(model))
        for (int iy = 0; iy < kNYr; ++iy) g.push_back({iy});
    else
    {
        std::vector<int> all;
        for (int iy = 0; iy < kNYr; ++iy) all.push_back(iy);
        g.push_back(all);
    }
    return g;
}
TString GroupName(const TString &model, int k)
{
    if (PerBin(model)) return TString::Format("y%d", k);
    return PerGroup(model) ? TString(pOAbcdY::kGroupName[k]) : TString("all");
}
// the fitted SR multiplier of (e, charge, group) as the extractor reports it
// (mtsplit: the bin's m_T >= 30 multiplier, which the extractor also writes
// under this key -- the channel where the W, and so the needed factor, sits)
RVal Multiplier(const FitR &f, int bn, int ic, int k)
{
    std::string key = std::string("qcd_abcd_mult_") + f.spec.flav.Data() + "_" + kChR[ic];
    if (PerGroup(f.model)) key += std::string("_") + pOAbcdY::kGroupName[k];
    if (PerBin(f.model)) key += Form("_y%d", k);
    auto it = f.par[bn].find(key);
    return it == f.par[bn].end() ? RVal() : it->second;
}

TGraphErrors *MakeGraph(const char *name, int bn, const double *v, const double *e)
{
    TGraphErrors *g = new TGraphErrors(kNYr);
    g->SetName(name);
    for (int iy = 0; iy < kNYr; ++iy)
    {
        g->SetPoint(iy, XCenter(bn, iy), v[iy]);
        g->SetPointError(iy, XHalf(bn, iy), e[iy]);
    }
    return g;
}

// one series of r per bin of fit f (total + stat graphs)
OverlaySeries RSeries(const FitR &f, int bn, int ic, const char *uniq, const std::string &label, int marker,
                      double size, double shift)
{
    double v[kNYr], e[kNYr], s[kNYr];
    for (int iy = 0; iy < kNYr; ++iy)
    {
        const RVal &x = f.r[bn][ic][iy];
        v[iy] = x.v; e[iy] = x.e; s[iy] = x.s;
    }
    OverlaySeries o;
    o.gTot = MakeGraph(Form("g_r_%s_%s_%s_%s", uniq, f.spec.tag.Data(), kChR[ic], kBnR[bn]), bn, v, e);
    o.gStat = MakeGraph(Form("g_rs_%s_%s_%s_%s", uniq, f.spec.tag.Data(), kChR[ic], kBnR[bn]), bn, v, s);
    o.label = label;
    o.color = f.spec.color;
    o.marker = marker;
    o.markerSize = size;
    o.xShift = shift;
    return o;
}

// r_e / r_mu per bin, statistical errors (independent fits)
OverlaySeries RatioSeries(const FitR &fm, const FitR &fe, int bn, int ic, const char *uniq, const std::string &label,
                          int color, int marker, double size, double shift)
{
    double v[kNYr], e[kNYr];
    for (int iy = 0; iy < kNYr; ++iy)
    {
        const RVal &m = fm.r[bn][ic][iy], &x = fe.r[bn][ic][iy];
        v[iy] = x.v / m.v;
        e[iy] = v[iy] * std::hypot(x.s / x.v, m.s / m.v);
    }
    OverlaySeries o;
    o.gTot = MakeGraph(Form("g_emu_%s_%s_%s", uniq, kChR[ic], kBnR[bn]), bn, v, e);
    o.label = label;
    o.color = color;
    o.marker = marker;
    o.markerSize = size;
    o.xShift = shift;
    return o;
}

// sum_i (a_i - b_i)^2 / (s_a^2 + s_b^2) over the 12 bins, with the stat or total errors
double Chi2(const FitR &a, const FitR &b, int bn, int ic, bool stat)
{
    double c = 0;
    for (int iy = 0; iy < kNYr; ++iy)
    {
        const RVal &x = a.r[bn][ic][iy], &y = b.r[bn][ic][iy];
        const double sa = stat ? x.s : x.e, sb = stat ? y.s : y.e;
        c += std::pow(x.v - y.v, 2) / (sa * sa + sb * sb);
    }
    return c;
}

// layout: legend upper left under "CMS"; header (two lines) upper right from
// x = 0.53, the info box right under it (y 0.58-0.74); the y ranges below keep
// the data under ~0.57 NDC
PlotStyle RStyle(int nLegRows, double yLegTop = 0.80)
{
    PlotStyle ps;
    ps.showStats = false;
    ps.logy = false;
    ps.legX1 = 0.17; ps.legX2 = 0.52;
    ps.legY2 = yLegTop; ps.legY1 = yLegTop - 0.042 * nLegRows;
    ps.headerX = 0.53;
    ps.subSize = 0.036;
    ps.boxTextSize = 0.027;
    ps.boxX1 = 0.57; ps.boxX2 = 0.95; ps.boxY1 = 0.58; ps.boxY2 = 0.74;
    return ps;
}

GraphTuner RangeTuner(double lo, double hi, double ref)
{
    return [lo, hi, ref](TCanvas *cv, TGraphErrors *g)
    {
        if (auto *h = g->GetHistogram()) { h->SetMinimum(lo); h->SetMaximum(hi); }
        g->SetMinimum(lo);
        g->SetMaximum(hi);
        TLine *l = new TLine(g->GetXaxis()->GetXmin(), ref, g->GetXaxis()->GetXmax(), ref);
        l->SetLineStyle(2);
        l->SetLineColor(kGray + 2);
        l->Draw("same");
        cv->Modified();
        cv->Update();
    };
}
} // namespace

void r_per_bin(const char *disc = "leppt_mt40_abcdy", const char *refDisc = "")
{
    TString dsuf, dlab;
    if (!pODisc::Spec(disc, dsuf, dlab)) return;
    const bool withRef = refDisc && *refDisc && TString(refDisc) != disc;

    FitR F[3], R[3]; // [mu, comb, e]: this tag, the reference tag
    for (int k = 0; k < 3; ++k) LoadFitR(disc, kFitTag[k], F[k]);
    if (!F[0].ok || !F[2].ok)
    {
        std::cerr << "[ERROR] r_per_bin: need both per-flavour fits (simfit_mu, simfit_ele) of " << disc
                  << " -- run the fork's run_pO_fits.sh both flavfit first\n";
        return;
    }
    bool haveRef = false;
    TString refLab;
    if (withRef)
    {
        TString rsuf;
        if (pODisc::Spec(refDisc, rsuf, refLab))
        {
            LoadFitR(refDisc, "simfit_mu", R[0]);
            LoadFitR(refDisc, "simfit_ele", R[2]);
            haveRef = R[0].ok && R[2].ok;
        }
        if (!haveRef)
            std::cerr << "[WARN] r_per_bin: reference '" << refDisc << "' has no per-flavour fits -> model plots skipped\n";
    }
    const TString model = F[2].model;
    const TString outDir = TString("./plots/flavfit/r/") + disc;
    gSystem->mkdir(outDir, kTRUE);
    std::ofstream csv((outDir + "/r_per_bin.csv").Data());
    csv << "quantity,tag,fit,binning,charge,bin,x,value,err_total,err_stat\n";

    printf("\n=== r_per_bin: postfit r per rapidity bin, disc %s (%s), multijet model of the e-only fit: %s ===\n", disc,
           dlab.Data(), ModelName(model).Data());
    for (int k = 0; k < 3; ++k)
        printf("  %-10s %s\n", kFitTag[k], F[k].ok ? F[k].work.Data() : "ABSENT");

    // reference provenance (see the header)
    bool refPreFix = false;
    if (haveRef && R[2].haveShapes[0] && R[2].haveShapes[1])
    {
        bool all1 = true;
        for (int k : {0, 2})
            for (int ic = 0; ic < 2; ++ic)
                if (R[k].qPreTot[0][ic] > 0)
                    all1 &= std::fabs(R[k].qPreTot[1][ic] / R[k].qPreTot[0][ic] - 1.0) < 1e-3;
        refPreFix = all1;
        if (refPreFix)
            printf("[WARN] reference %s: fb prefit multijet == lab one -> its fb fits predate the 2026-10-19 FB-window\n"
                   "       fix; the fb model plots mix that fix with the model change (lab is clean)\n", refDisc);
    }

    for (int bn = 0; bn < 2; ++bn)
    {
        const char *xT = bn == 0 ? "#eta^{l}_{lab}" : "#eta^{l}_{CM}";
        // ==== 1. r per bin: mu-only | mu+e | e-only ===========================
        for (int ic = 0; ic < 2; ++ic)
        {
            printf("\n[r] %s %s  (r: value +- total (stat))\n", kChR[ic], kBnR[bn]);
            printf("   %-4s %7s | %-22s | %-22s | %-22s | %-15s | %s\n", "bin", "x", "mu only", "mu + e", "e only",
                   "e/mu (stat)", "(e - mu)/sigma stat, total");
            for (int iy = 0; iy < kNYr; ++iy)
            {
                std::string line = Form("   y%-3d %7.3f", iy, XCenter(bn, iy));
                for (int k = 0; k < 3; ++k)
                {
                    const RVal &v = F[k].r[bn][ic][iy];
                    line += F[k].ok ? Form(" | %6.3f +- %.3f (%.3f)", v.v, v.e, v.s) : " | -                     ";
                    if (F[k].ok)
                        csv << "r," << disc << "," << kFitTag[k] << "," << kBnR[bn] << "," << kChR[ic] << "," << iy
                            << "," << XCenter(bn, iy) << "," << v.v << "," << v.e << "," << v.s << "\n";
                }
                const RVal &m = F[0].r[bn][ic][iy], &e = F[2].r[bn][ic][iy];
                const double q = e.v / m.v, qe = q * std::hypot(e.s / e.v, m.s / m.v);
                line += Form(" | %.3f +- %.3f | %+5.2f  %+5.2f", q, qe, (e.v - m.v) / std::hypot(e.s, m.s),
                             (e.v - m.v) / std::hypot(e.e, m.e));
                csv << "r_e_over_mu," << disc << ",ratio," << kBnR[bn] << "," << kChR[ic] << "," << iy << ","
                    << XCenter(bn, iy) << "," << q << ",," << qe << "\n";
                printf("%s\n", line.c_str());
            }
            const double c2s = Chi2(F[0], F[2], bn, ic, true), c2t = Chi2(F[0], F[2], bn, ic, false);
            printf("   mu vs e: chi2/ndf = %.1f/12 (stat), %.1f/12 (total)\n", c2s, c2t);
            csv << "chi2_mu_vs_e_stat," << disc << ",," << kBnR[bn] << "," << kChR[ic] << ",,," << c2s << ",,\n";
            csv << "chi2_mu_vs_e_total," << disc << ",," << kBnR[bn] << "," << kChR[ic] << ",,," << c2t << ",,\n";

            std::vector<OverlaySeries> ser;
            const double sh[3] = {-0.6, 0.0, 0.6};
            for (int k = 0; k < 3; ++k)
                if (F[k].ok)
                    ser.push_back(RSeries(F[k], bn, ic, "fit", F[k].spec.label.Data(), F[k].spec.marker,
                                          F[k].spec.markerSize, sh[k]));
            PlotStyle ps = RStyle((int)ser.size() + 1);
            SaveNiceGraph_Overlay(ser, Form("%s/r_%s_%s", outDir.Data(), kChR[ic], kBnR[bn]), xT,
                                  "postfit r  (data / POWHEG)", "", kChLab[ic], "postfit signal strength",
                                  {dlab.Data(), Form("multijet: %s", ModelName(model).Data()),
                                   Form("#mu vs e #chi^{2}: %.1f/12 stat, %.1f/12 tot.", c2s, c2t)},
                                  ps, RangeTuner(0.0, 3.2, 1.0));
        }

        // ==== 2. r_e / r_mu, W+ and W- ======================================
        {
            std::vector<OverlaySeries> ser;
            ser.push_back(RatioSeries(F[0], F[2], bn, 0, "nom", "W^{+}", kAzure + 2, 21, 1.5, -0.35));
            ser.push_back(RatioSeries(F[0], F[2], bn, 1, "nom", "W^{-}", kRed + 1, 33, 2.0, 0.35));
            PlotStyle ps = RStyle(2);
            SaveNiceGraph_Overlay(ser, Form("%s/r_emu_%s", outDir.Data(), kBnR[bn]), xT, "r_{e} / r_{#mu}  (postfit)",
                                  "", "e-only fit / #mu-only fit", "postfit signal strength",
                                  {dlab.Data(), Form("multijet: %s", ModelName(model).Data()),
                                   "errors: statistical (independent fits)"},
                                  ps, RangeTuner(0.2, 2.6, 1.0));
        }

        // ==== 3. the multijet each group would need for r_e = r_mu ===========
        if (F[2].haveShapes[bn])
        {
            const std::vector<std::vector<int>> grp = GroupsOf(model);
            for (int ic = 0; ic < 2; ++ic)
            {
                printf("\n[needed] %s %s: factor on the e-only POSTFIT multijet that makes r_e = r_mu\n", kChR[ic],
                       kBnR[bn]);
                printf("   per bin: %-4s | %8s %8s %7s | %8s %8s | %s\n", "bin", "Q_post", "S_e", "data", "dQ", "f_bin",
                       "group");
                for (int iy = 0; iy < kNYr; ++iy)
                {
                    const RVal &m = F[0].r[bn][ic][iy], &e = F[2].r[bn][ic][iy];
                    const double S = F[2].W[bn][ic][iy] / e.v, Q = F[2].Q[bn][ic][iy];
                    const double dQ = (m.v - e.v) * S;
                    const double ef = std::sqrt(S * S * m.s * m.s + F[2].D[bn][ic][iy]) / Q;
                    int g = 0;
                    for (size_t k = 0; k < grp.size(); ++k)
                        for (int b : grp[k]) if (b == iy) g = (int)k;
                    printf("            y%-3d | %8.1f %8.1f %7.0f | %+8.1f %5.2f+-%.2f | %s\n", iy, Q, S,
                           F[2].D[bn][ic][iy], dQ, 1.0 - dQ / Q, ef, GroupName(model, g).Data());
                }
                double gx[3], gxe[3], gf[3], gfe[3], band[3], pull[3];
                const int nG = (int)grp.size();
                printf("   per group (the granularity of the fitted multijet); pull = (f - 1) / sqrt(err_f^2 + band^2):\n");
                for (int k = 0; k < nG; ++k)
                {
                    double Qg = 0, dQ = 0, var = 0;
                    for (int iy : grp[k])
                    {
                        const RVal &m = F[0].r[bn][ic][iy], &e = F[2].r[bn][ic][iy];
                        const double S = F[2].W[bn][ic][iy] / e.v;
                        Qg += F[2].Q[bn][ic][iy];
                        dQ += (m.v - e.v) * S;
                        var += S * S * m.s * m.s + F[2].D[bn][ic][iy];
                    }
                    const RVal mult = Multiplier(F[2], bn, ic, k);
                    band[k] = mult.ok ? mult.e / mult.v : NAN;
                    gf[k] = 1.0 - dQ / Qg;
                    gfe[k] = std::sqrt(var) / Qg;
                    const double *ed = bn == 0 ? pOSkim::kYEdges : pOSkim::kYEdgesFB;
                    const double sft = bn == 0 ? 0.0 : pOAnalysis::kPORapidityShift;
                    double lo = 1e9, hi = -1e9;
                    for (int iy : grp[k]) { lo = std::min(lo, ed[iy] - sft); hi = std::max(hi, ed[iy + 1] - sft); }
                    gx[k] = 0.5 * (lo + hi);
                    gxe[k] = 0.5 * (hi - lo);
                    pull[k] = (gf[k] - 1.0) / std::hypot(gfe[k], std::isfinite(band[k]) ? band[k] : 0.0);
                    printf("     %-4s Q_post %7.1f  needed dQ %+7.1f -> f = %.3f +- %.3f   fitted multiplier %.3f +- %.3f"
                           " (band +-%.1f%%)  pull %+.2f\n",
                           GroupName(model, k).Data(), Qg, dQ, gf[k], gfe[k], mult.v, mult.e, 100 * band[k], pull[k]);
                    csv << "qcd_needed_factor," << disc << ",simfit_ele," << kBnR[bn] << "," << kChR[ic] << ","
                        << GroupName(model, k) << "," << gx[k] << "," << gf[k] << "," << gfe[k] << ",\n";
                    csv << "qcd_fitted_relerr," << disc << ",simfit_ele," << kBnR[bn] << "," << kChR[ic] << ","
                        << GroupName(model, k) << "," << gx[k] << "," << band[k] << ",,\n";
                }
                TGraphErrors *gN = new TGraphErrors(nG);
                gN->SetName(Form("g_needed_%s_%s", kChR[ic], kBnR[bn]));
                double ymax = 2.4;
                std::string pulls;
                for (int k = 0; k < nG; ++k)
                {
                    gN->SetPoint(k, gx[k], gf[k]);
                    gN->SetPointError(k, gxe[k], gfe[k]);
                    ymax = std::max(ymax, 2.3 * (gf[k] + gfe[k]));
                    pulls += Form("%s%s %+.1f#sigma", k ? ",  " : "", GroupName(model, k).Data(), pull[k]);
                }
                OverlaySeries s;
                s.gTot = gN;
                s.label = "needed for r_{e} = r_{#mu}";
                s.color = ic == 0 ? kAzure + 2 : kRed + 1;
                s.marker = ic == 0 ? 21 : 33;
                s.markerSize = ic == 0 ? 1.5 : 2.0;
                const std::vector<double> bLo(gx, gx + nG), bW(gxe, gxe + nG), bB(band, band + nG);
                GraphTuner tune = [bLo, bW, bB, ymax](TCanvas *cv, TGraphErrors *g)
                {
                    g->SetMinimum(0.0);
                    g->SetMaximum(ymax);
                    if (auto *h = g->GetHistogram()) { h->SetMinimum(0.0); h->SetMaximum(ymax); }
                    for (size_t k = 0; k < bLo.size(); ++k)
                    {
                        if (!std::isfinite(bB[k])) continue;
                        TBox *b = new TBox(bLo[k] - bW[k], 1.0 - bB[k], bLo[k] + bW[k], 1.0 + bB[k]);
                        b->SetFillColorAlpha(kGray + 1, 0.45);
                        b->SetLineColor(kGray + 2);
                        b->Draw();
                    }
                    TLine *l = new TLine(g->GetXaxis()->GetXmin(), 1.0, g->GetXaxis()->GetXmax(), 1.0);
                    l->SetLineStyle(2);
                    l->SetLineColor(kGray + 3);
                    l->Draw("same");
                    cv->Modified();
                    cv->Update();
                };
                PlotStyle ps = RStyle(1);
                ps.boxY1 = 0.49; ps.boxY2 = 0.75; // five lines
                SaveNiceGraph_Overlay({s}, Form("%s/qcd_needed_%s_%s", outDir.Data(), kChR[ic], kBnR[bn]), xT,
                                      "factor on the e postfit multijet", "", kChLab[ic], "e-only fit, postfit",
                                      {Form("multijet: %s", ModelName(model).Data()),
                                       "grey: fitted multijet #pm 1#sigma per group",
                                       "points: factor that gives r_{e} = r_{#mu}",
                                       "(need #minus fitted) / error:", pulls},
                                      ps, tune);
            }
        }
        else
            printf("\n[needed] %s: no e-only fitDiagnostics shapes -> qcd_needed plots skipped\n", kBnR[bn]);

        // ==== 4. the reference model vs this one =============================
        if (!haveRef) continue;
        const bool note = (bn == 1 && refPreFix);
        const TString refName = ModelName(R[2].model), nomName = ModelName(model);
        // short legend tags: the part after "in-fit ABCD " when both are in-fit ABCD models
        auto shortName = [](const TString &m) -> TString
        {
            TString s = ModelName(m);
            return s.BeginsWith("in-fit ABCD ") ? TString(s(12, s.Length() - 12)) : s;
        };
        const TString refShort = shortName(R[2].model), nomShort = shortName(model);
        for (int ic = 0; ic < 2; ++ic)
        {
            // e postfit multijet per group, this model / reference (both from their own shapes)
            std::vector<std::string> box = {Form("open: %s", refName.Data()), Form("filled: %s", nomName.Data())};
            if (F[2].haveShapes[bn] && R[2].haveShapes[bn])
            {
                box.push_back("e postfit multijet, filled / open:");
                std::string l;
                for (int k = 0; k < pOAbcdY::kNGroups; ++k)
                {
                    double a = 0, b = 0;
                    for (int iy : pOAbcdY::BinsOf(k)) { a += F[2].Q[bn][ic][iy]; b += R[2].Q[bn][ic][iy]; }
                    l += Form("%s%s %.2f", k ? ",  " : "", pOAbcdY::kGroupName[k], a / b);
                    csv << "e_qcd_post_this_over_ref," << disc << ",simfit_ele," << kBnR[bn] << "," << kChR[ic] << ","
                        << pOAbcdY::kGroupName[k] << ",," << a / b << ",,\n";
                }
                box.push_back(l);
            }
            double maxShift = 0;
            for (int iy = 0; iy < kNYr; ++iy)
                maxShift = std::max(maxShift, std::fabs(F[2].r[bn][ic][iy].v - R[2].r[bn][ic][iy].v));
            box.push_back(Form("max |#Delta r_{e}| between models: %.3f", maxShift));
            if (note) box.push_back("ref. fb fit: pre-FB-window-fix inputs");
            printf("\n[model] %s %s: e-only r, %s -> %s: max |delta r_e| = %.3f\n", kChR[ic], kBnR[bn], refName.Data(),
                   nomName.Data(), maxShift);
            std::vector<OverlaySeries> ser;
            ser.push_back(RSeries(R[0], bn, ic, "ref", Form("#mu only, %s", refShort.Data()), 24, 1.3, -0.66));
            ser.push_back(RSeries(F[0], bn, ic, "nom", Form("#mu only, %s", nomShort.Data()), 20, 1.3, -0.22));
            ser.push_back(RSeries(R[2], bn, ic, "ref", Form("e only, %s", refShort.Data()), 25, 1.5, 0.22));
            ser.push_back(RSeries(F[2], bn, ic, "nom", Form("e only, %s", nomShort.Data()), 21, 1.5, 0.66));
            PlotStyle ps = RStyle(5);
            ps.boxY1 = 0.50; ps.boxY2 = 0.76; // five or six lines
            SaveNiceGraph_Overlay(ser, Form("%s/r_model_%s_%s", outDir.Data(), kChR[ic], kBnR[bn]), xT,
                                  "postfit r  (data / POWHEG)", "", kChLab[ic], "postfit, two multijet models", box, ps,
                                  RangeTuner(0.0, 3.6, 1.0));
        }
        {
            std::vector<OverlaySeries> ser;
            ser.push_back(RatioSeries(R[0], R[2], bn, 0, "ref", Form("W^{+}, %s", refShort.Data()), kAzure + 2, 25, 1.5, -0.66));
            ser.push_back(RatioSeries(F[0], F[2], bn, 0, "nom", Form("W^{+}, %s", nomShort.Data()), kAzure + 2, 21, 1.5, -0.22));
            ser.push_back(RatioSeries(R[0], R[2], bn, 1, "ref", Form("W^{-}, %s", refShort.Data()), kRed + 1, 27, 2.0, 0.22));
            ser.push_back(RatioSeries(F[0], F[2], bn, 1, "nom", Form("W^{-}, %s", nomShort.Data()), kRed + 1, 33, 2.0, 0.66));
            std::vector<std::string> box = {Form("open: %s", refName.Data()), Form("filled: %s", nomName.Data()),
                                            "errors: statistical (independent fits)"};
            if (note) box.push_back("ref. fb fit: pre-FB-window-fix inputs");
            PlotStyle ps = RStyle(4);
            SaveNiceGraph_Overlay(ser, Form("%s/r_emu_model_%s", outDir.Data(), kBnR[bn]), xT,
                                  "r_{e} / r_{#mu}  (postfit)", "", "e-only fit / #mu-only fit",
                                  "postfit, two multijet models", box, ps, RangeTuner(0.2, 2.6, 1.0));
        }
    }
    csv.close();
    printf("\n[r_per_bin] wrote %s/{r_{Wp,Wm},r_emu,qcd_needed_{Wp,Wm}%s}_{lab,fb}.{png,pdf} + r_per_bin.csv\n",
           outDir.Data(), haveRef ? ",r_model_{Wp,Wm},r_emu_model" : "");
}
