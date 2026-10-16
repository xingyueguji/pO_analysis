// =============================================================================
// plotting/prefit_r.C -- the PREFIT signal strength per rapidity bin, mu vs e
// (2026-10-13).
//
//   r_pre = (data - DY - QCD) / W     per region of the lepton-pT Combine inputs
//     W   = signal + wtau              (what the fit's r scales)
//     DY  = z + ztau
//     QCD = the DISPLAY template `qcd` that the leppt_mt40 stacks draw, so every
//           number here equals the r_pre line of that stack's info box
//           (mtandmet.C::prefitLines). The in-fit-ABCD prefit `qcd_abcd` (what
//           the abcd-mode fit starts from) is printed and written to the CSV
//           alongside: the two share the per-y split and differ only in the
//           per-charge total (inclusive r differs by <= 0.03).
//   r_pre > 1: the data sit above the absolute MC, so the fit will pull r up.
//   MC/data = (W + DY + QCD)/data is printed too; it is NOT 1/r_pre, since r
//   leaves the backgrounds alone.
//   Error: data Poisson only, sqrt(N)/W (the template stat is far smaller).
//
// Plots (mu and e overlaid, plotting_helper.C::SaveNiceGraph_Overlay with the
// fit_variants.h flavour styles, the inclusive r of each flavour in the
// legend), for W+, W- and W (= the two charge regions of a bin summed), in the
// lab and the fb binning:
//   plots/prefit_r/<disc>/prefit_r_{Wp,Wm,W}_{lab,fb}.{png,pdf}
//   plots/prefit_r/<disc>/prefit_r.csv   (every number, both QCD templates)
// The inclusive rows are the sum of the 12 lab regions; the macro checks them
// against the inputs' own <C>_incl directories.
//
// Bin positions come from the skim's single-source binning (pOSkim::kYEdges[FB]);
// fb bins are drawn at eta_CM = y_lab - kPORapidityShift (symmetric).
//
// Usage (from plotting/):  root -l -b -q 'prefit_r.C+("leppt_mt40")'
// Run by run_combine_inputs.sh after the inputs are written (log
// logs/prefit_r_<disc>.log). Needs both flavours' inputs for the overlay; one
// alone is plotted with a [WARN].
// =============================================================================
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLine.h"
#include "TString.h"
#include "TSystem.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "plotting_helper.C"              // PlotStyle, SaveNiceGraph_Overlay, CMS_lumi
#include "fit_variants.h"                 // pOFit::Spec (mu / e overlay styles)
#include "../skim/skim_common.h"          // pOSkim::kNY, kYEdges, kYEdgesFB (single-source binning)
#include "../analysis/analysis_helpers.h" // pOAnalysis::kPORapidityShift

namespace
{
// Integrated templates of one region (bins 1..N, as in the stacks' info boxes)
struct Counts
{
    double D = 0, W = 0, Y = 0, Q = 0, QA = 0;
    bool haveQA = true;
    void add(const Counts &o)
    {
        D += o.D; W += o.W; Y += o.Y; Q += o.Q; QA += o.QA;
        haveQA = haveQA && o.haveQA;
    }
    double r() const { return (D - Y - Q) / W; }
    double rErr() const { return std::sqrt(D) / W; }
    double rAbcd() const { return haveQA ? (D - Y - QA) / W : NAN; }
    double mcOverData() const { return (W + Y + Q) / D; }
};

bool ReadRegion(TFile *f, const std::string &dir, Counts &c)
{
    auto integ = [&](const char *name, bool required, bool &ok) -> double
    {
        TH1 *h = (TH1 *)f->Get((dir + "/" + name).c_str());
        if (!h)
        {
            if (required)
            {
                std::cerr << "[ERROR] prefit_r: no " << dir << "/" << name << " in " << f->GetName() << "\n";
                ok = false;
            }
            return 0.0;
        }
        return h->Integral(1, h->GetNbinsX());
    };
    bool ok = true, okQA = true;
    c.D = integ("data_obs", true, ok);
    c.W = integ("signal", true, ok) + integ("wtau", true, ok);
    c.Y = integ("z", true, ok) + integ("ztau", true, ok);
    c.Q = integ("qcd", true, ok);
    c.QA = integ("qcd_abcd", false, okQA);
    c.haveQA = (f->Get((dir + "/qcd_abcd").c_str()) != nullptr);
    return ok && c.W > 0.0;
}
} // namespace

void prefit_r(const char *disc = "leppt_mt40")
{
    // Only the lepton-pT inputs carry the stacks these numbers annotate. (A MET
    // version would need a decision: the per-y MET QCD is anchored at MET < 30
    // with r = 1, so r_pre over the full MET range is pulled to 1 by design.)
    if (TString(disc) != "leppt_mt40")
    {
        std::cerr << "[ERROR] prefit_r: only disc 'leppt_mt40' is supported (got '" << disc << "')\n";
        return;
    }
    const int NY = pOSkim::kNY;
    const char *flav[2] = {"mu", "ele"};
    const char *fdir[2] = {"./plots", "./plots/Elec"};
    const char *bnName[2] = {"lab", "fb"};
    const char *chg[3] = {"Wp", "Wm", "W"};

    // c[flavour][binning][charge][y], y = NY -> inclusive
    std::vector<Counts> c(2 * 2 * 3 * (NY + 1));
    auto C = [&](int f, int b, int q, int y) -> Counts & { return c[((f * 2 + b) * 3 + q) * (NY + 1) + y]; };
    bool haveFlav[2] = {false, false};

    for (int f = 0; f < 2; ++f)
    {
        const TString fn = TString::Format("%s/combine_input_W_%s.root", fdir[f], disc);
        TFile *fin = TFile::Open(fn, "READ");
        if (!fin || fin->IsZombie())
        {
            std::cerr << "[WARN] prefit_r: cannot open " << fn << " -> " << flav[f] << " not plotted\n";
            continue;
        }
        bool ok = true;
        for (int b = 0; b < 2; ++b)
            for (int q = 0; q < 2; ++q)
                for (int y = 0; y < NY; ++y)
                {
                    Counts &r = C(f, b, q, y);
                    ok &= ReadRegion(fin, Form("%s_%s_y%d", chg[q], bnName[b], y), r);
                    C(f, b, 2, y).add(r); // W = W+ + W- of the bin
                    C(f, b, q, NY).add(r);
                    C(f, b, 2, NY).add(r);
                }
        // consistency: the summed lab regions vs the inputs' own inclusive dirs
        for (int q = 0; q < 3 && ok; ++q)
        {
            Counts inc;
            if (!ReadRegion(fin, Form("%s_incl", chg[q]), inc)) continue;
            const Counts &s = C(f, 0, q, NY);
            const double dd = std::fabs(inc.D - s.D), dw = std::fabs(inc.W - s.W) / s.W;
            std::cout << "[prefit_r] " << flav[f] << " " << chg[q] << "_incl vs sum of lab regions: data "
                      << inc.D << " vs " << s.D << ", W rel. diff " << dw
                      << ((dd > 0.5 || dw > 1e-6) ? "   [WARN] mismatch" : "   ok") << "\n";
        }
        fin->Close();
        delete fin;
        if (!ok)
        {
            std::cerr << "[ERROR] prefit_r: incomplete regions in " << fn << " -> " << flav[f] << " not plotted\n";
            continue;
        }
        haveFlav[f] = true;
    }
    if (!haveFlav[0] && !haveFlav[1])
    {
        std::cerr << "[ERROR] prefit_r: no usable input -- run run_combine_inputs.sh first\n";
        return;
    }

    // ---- console + CSV ------------------------------------------------------
    const std::string outDir = std::string("./plots/prefit_r/") + disc;
    gSystem->mkdir(outDir.c_str(), kTRUE);
    std::ofstream csv(outDir + "/prefit_r.csv");
    csv << "flavour,binning,charge,ybin,x_center,data,W,DY,QCD_display,QCD_abcd,"
           "r_pre,r_pre_err,r_pre_abcd,MC_over_data\n";
    auto xCenter = [&](int b, int y) -> double
    {
        if (b == 0) return 0.5 * (pOSkim::kYEdges[y] + pOSkim::kYEdges[y + 1]);
        return 0.5 * (pOSkim::kYEdgesFB[y] + pOSkim::kYEdgesFB[y + 1]) - pOAnalysis::kPORapidityShift;
    };
    auto xHalf = [&](int b, int y) -> double
    {
        const double *e = (b == 0) ? pOSkim::kYEdges : pOSkim::kYEdgesFB;
        return 0.5 * (e[y + 1] - e[y]);
    };
    for (int b = 0; b < 2; ++b)
        for (int q = 0; q < 3; ++q)
        {
            std::cout << "\n[prefit_r] " << chg[q] << ", " << bnName[b]
                      << " binning  (r_pre = (data - DY - QCD)/W, QCD as drawn; [abcd] = with the in-fit prefit qcd_abcd)\n";
            std::cout << Form("   %-5s %7s | %-17s %-8s | %-17s %-8s | %s\n", "bin", "x",
                              "mu r_pre", "[abcd]", "e r_pre", "[abcd]", "e/mu");
            for (int y = 0; y <= NY; ++y)
            {
                const bool inc = (y == NY);
                if (inc && b == 1) continue; // the inclusive is the lab sum
                std::string line = inc ? "   incl          " : Form("   y%-4d %7.3f ", y, xCenter(b, y));
                double rv[2] = {NAN, NAN}, re[2] = {NAN, NAN};
                for (int f = 0; f < 2; ++f)
                {
                    if (!haveFlav[f]) { line += Form("| %-17s %-8s ", "-", "-"); continue; }
                    const Counts &k = C(f, b, q, y);
                    rv[f] = k.r(); re[f] = k.rErr();
                    line += Form("| %.3f +- %.3f   %-8s ", rv[f], re[f],
                                 std::isnan(k.rAbcd()) ? "n/a" : Form("%.3f", k.rAbcd()));
                    csv << flav[f] << "," << bnName[b] << "," << chg[q] << ","
                        << (inc ? std::string("incl") : std::to_string(y)) << ","
                        << (inc ? 0.0 : xCenter(b, y)) << "," << k.D << "," << k.W << "," << k.Y << ","
                        << k.Q << "," << (k.haveQA ? k.QA : NAN) << "," << rv[f] << "," << re[f] << ","
                        << k.rAbcd() << "," << k.mcOverData() << "\n";
                }
                if (haveFlav[0] && haveFlav[1])
                {
                    const double ratio = rv[1] / rv[0];
                    line += Form("| %.3f +- %.3f", ratio, ratio * std::hypot(re[0] / rv[0], re[1] / rv[1]));
                }
                std::cout << line << "\n";
            }
        }
    csv.close();

    // ---- plots ----------------------------------------------------------------
    pOFit::Spec spec[2];
    pOFit::Get("simfit_mu", spec[0]);
    pOFit::Get("simfit_ele", spec[1]);
    const char *flavLabel[2] = {"#mu", "e"};
    const char *sub1[3] = {"W^{+} #rightarrow l^{+} #nu", "W^{-} #rightarrow l^{-} #bar{#nu}", "W #rightarrow l #nu"};
    const int nFlav = (haveFlav[0] ? 1 : 0) + (haveFlav[1] ? 1 : 0);
    if (nFlav == 1)
        std::cerr << "[WARN] prefit_r: only one flavour available -> no mu/e overlay\n";

    PlotStyle ps;
    ps.showStats = false;
    ps.logy = false;
    // legend lower left, wide enough for the inclusive value in each entry
    ps.legX1 = 0.18; ps.legX2 = 0.62;
    ps.legY1 = 0.15; ps.legY2 = 0.15 + 0.042 * nFlav;
    GraphTuner tune = [](TCanvas *cv, TGraphErrors *g)
    {
        if (auto *h = g->GetHistogram()) { h->SetMinimum(0.0); h->SetMaximum(2.4); }
        const double xmin = g->GetXaxis()->GetXmin(), xmax = g->GetXaxis()->GetXmax();
        TLine *l1 = new TLine(xmin, 1.0, xmax, 1.0);
        l1->SetLineStyle(2);
        l1->Draw("same");
        cv->Modified();
        cv->Update();
    };

    for (int b = 0; b < 2; ++b)
        for (int q = 0; q < 3; ++q)
        {
            std::vector<OverlaySeries> series;
            for (int f = 0; f < 2; ++f)
            {
                if (!haveFlav[f]) continue;
                TGraphErrors *g = new TGraphErrors(NY);
                g->SetName(Form("g_prefit_r_%s_%s_%s", flav[f], chg[q], bnName[b]));
                for (int y = 0; y < NY; ++y)
                {
                    g->SetPoint(y, xCenter(b, y), C(f, b, q, y).r());
                    g->SetPointError(y, xHalf(b, y), C(f, b, q, y).rErr());
                }
                const Counts &inc = C(f, 0, q, NY); // inclusive = the lab sum
                OverlaySeries s;
                s.gTot = g;
                s.label = Form("%s,  inclusive r_{pre} = %.2f #pm %.2f", flavLabel[f], inc.r(), inc.rErr());
                s.color = spec[f].color;
                s.marker = spec[f].marker;
                s.markerSize = spec[f].markerSize;
                s.xShift = (nFlav == 2) ? (f == 0 ? -0.35 : 0.35) : 0.0;
                series.push_back(s);
            }
            SaveNiceGraph_Overlay(series, outDir + Form("/prefit_r_%s_%s", chg[q], bnName[b]),
                                  b == 0 ? "#eta^{l}_{lab}" : "#eta^{l}_{CM}",
                                  "r_{pre} = (data#minusbkg)/W", "", sub1[q],
                                  "prefit, m_{T} > 40 GeV", {}, ps, tune);
        }
    std::cout << "\n[prefit_r] wrote " << outDir << "/prefit_r_{Wp,Wm,W}_{lab,fb}.{png,pdf} + prefit_r.csv\n";
}
