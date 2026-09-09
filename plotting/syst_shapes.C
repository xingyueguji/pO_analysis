// =============================================================================
// syst_shapes.C -- diagnostics of the LHE shape systematics carried by the
// Combine inputs (2026-09-07): nPDF / qcdScale / alphaS, <process>_<syst>Up/Down
// next to every MC template (written by mtandmet.C / dileptonpeak.C from the
// skim twins, see CLAUDE.md "Structured inputs").
//
//   (1) PER-REGION SHAPE PLOTS  plots[/Elec]/syst_shapes/<disc>/perbin/<region>_<process>
//       Up/nominal and Down/nominal bin by bin for the three systematics on one
//       canvas (solid = Up, dashed = Down; nPDF red, qcdScale blue, alphaS
//       green), info box = integral shifts + the 68% range of the per-bin
//       ratios + the number of one-sided bins (Up and Down on the same side of
//       the nominal -- Combine only warns about those; alphaS can have some
//       because its two members are used as they are).
//   (2) SUMMARIES  plots[/Elec]/syst_shapes/<disc>/summary_<syst>_<binning>
//       the integral shifts of `signal` vs rapidity bin (W+ / W- x Up / Down),
//       lab and fb.
//   (3) INCLUSIVE CONSISTENCY (printed, [INCL] lines -- the log is the record):
//       (a) the TRUE inclusive uncertainty of the selected events: the member
//           twins (<h>_epps21/_scale/_alphas in the skim files) summed over the
//           24 lab templates, then combined -- asymmetric Hessian over the 53
//           EPPS21 pairs / 1.645 (= LHAPDF's PDFSet.uncertainty, pOLhe::Hessian),
//           the scale envelope, the alpha_s members -- per sample and, with the
//           k_s of mc_norm.h, per Combine process;
//       (b) the SUM OF THE PER-BIN Up/Down templates of the Combine inputs over
//           the 24 lab SR regions = what ONE collapsed nuisance implies in the
//           fit (a single theta moves every bin coherently). (b) >= (a) for the
//           Hessian and the envelope (triangle inequality: sum of per-bin
//           combinations vs combination of the sum); the ratio (b)/(a) is the
//           inflation of the inclusive uncertainty caused by collapsing the
//           eigen-directions into one nuisance;
//       (c) the ALL-EVENTS reference of skim/output/lhe_weights.txt (per-index
//           S_i/S_0-1 of the mu-flavour files): nuclear Hessian (idx 111-158),
//           scale envelope (idx 1-8), alpha_s (idx 105/104) -- compared with (a)
//           per sample. Agreement is expected only in pattern and magnitude:
//           (a) is the selected fiducial events, (c) all generated events, and
//           the tau / electron files have no reference row of their own (the
//           W+/W-/DY rows are used).
//   (4) LHAPDF-vs-pOLhe::Hessian CLOSURE ([CLOSURE] lines): the stored
//       <h>_nPDFUp/Down (written by lhe_updown.py with LHAPDF) recomputed bin by
//       bin from the _epps21 twin with pOLhe::Hessian -- max |difference| over
//       every template and file must be ~0 (two implementations, same data).
//   (5) MEMBER OVERLAYS  plots[/Elec]/syst_shapes/<disc>/members/<region>_signal
//       (2026-09-08): one canvas per W fit region (charge x lab/fb x y0..11) with
//       the nominal `signal` template (black) and ALL 106 EPPS21 member templates
//       (grey, transparent) on the same axes, rebuilt from the skim twins
//       <h>_epps21 with the k_s and W+/W- sample sum of mtandmet.C (member 0 is
//       checked against the Combine-input `signal` bin by bin); ratio pad =
//       member / nominal (grey) with the LHAPDF nPDF Up/Down of the region (red
//       solid / dashed) -- i.e. what PDFSet.uncertainty() made of the cloud.
//
//   root -l -b -q 'syst_shapes.C+("leppt_mt40")'      (or "met")
//   ./run_syst_shapes.sh [met|leppt_mt40|all]        -> logs/syst_shapes_<disc>.log
// =============================================================================
#include "plotting_helper.C"
#include "disc_variants.h"
#include "../skim/lhe_index.h" // pOLhe::kLheSystNames, kNMembers, Block, Hessian()
#include "../skim/mc_norm.h"   // pONorm::MCScale -> k_s per sample

#include "TBox.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TKey.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

const int kNSyst = pOLhe::kNLheSysts;
const int kColor[3] = {kRed + 1, kBlue + 1, kGreen + 2};

struct Shift { double up = 0, dn = 0; bool ok = false; }; // relative integral shifts

double Integ(const TH1 *h) { return h ? h->Integral(1, h->GetNbinsX()) : 0.0; }

Shift IntegralShift(const TH1 *nom, const TH1 *up, const TH1 *dn)
{
    Shift s;
    const double i0 = Integ(nom);
    if (!nom || !up || !dn || i0 <= 0) return s;
    s.ok = true;
    s.up = Integ(up) / i0 - 1.0;
    s.dn = Integ(dn) / i0 - 1.0;
    return s;
}

// ratio histogram (variation / nominal), 1 where the nominal is empty, no errors
TH1D *Ratio(const TH1 *var, const TH1 *nom, const char *name)
{
    TH1D *r = (TH1D *)nom->Clone(name);
    r->SetDirectory(nullptr);
    // the Combine-input templates carry the stack plots' fill attributes; a
    // filled "hist" would paint over the curves drawn before it
    r->SetFillStyle(0);
    r->SetFillColor(0);
    for (int b = 0; b <= nom->GetNbinsX() + 1; ++b)
    {
        const double n = nom->GetBinContent(b);
        r->SetBinContent(b, n > 0 ? var->GetBinContent(b) / n : 1.0);
        r->SetBinError(b, 0.0);
    }
    return r;
}

// median and central-68% range of the per-bin ratios (bins with nominal > 0 in [b1,b2])
void RatioSpread(const TH1 *nom, const TH1 *up, const TH1 *dn, int b1, int b2,
                 double &med, double &lo, double &hi)
{
    std::vector<double> v;
    for (int b = b1; b <= b2; ++b)
        if (nom->GetBinContent(b) > 0)
        {
            v.push_back(up->GetBinContent(b) / nom->GetBinContent(b));
            v.push_back(dn->GetBinContent(b) / nom->GetBinContent(b));
        }
    med = lo = hi = 1.0;
    if (v.empty()) return;
    std::sort(v.begin(), v.end());
    auto q = [&](double p) { return v[std::min((size_t)(p * v.size()), v.size() - 1)]; };
    med = q(0.5); lo = q(0.16); hi = q(0.84);
}

int OneSided(const TH1 *nom, const TH1 *up, const TH1 *dn, int b1, int b2)
{
    int n = 0;
    for (int b = b1; b <= b2; ++b)
    {
        const double x = nom->GetBinContent(b);
        if (x <= 0) continue;
        const double du = up->GetBinContent(b) - x, dd = dn->GetBinContent(b) - x;
        if (du * dd > 0) ++n; // both above or both below the nominal
    }
    return n;
}

// ---------------------------------------------------------------------------
// (1) one canvas per (region, process): the six ratio curves
// ---------------------------------------------------------------------------
void DrawRatioSet(const TH1 *nom, const std::vector<TH1 *> &ups, const std::vector<TH1 *> &dns,
                  const std::string &outNoExt, const char *xTitle, const std::string &head,
                  const std::string &sub1, double xlo, double xhi, bool logx = false)
{
    // layout: CMS_lumi owns the top-left corner -> header top-right, legend
    // right below it, info box lower-left; the y range is stretched downward
    // so the curves (near 1) sit above the box
    PlotStyle ps;
    ps.boxTextSize = 0.026;
    ps.headerX = 0.47; ps.headerY = 0.90; ps.headerDy = 0.05;
    ps.titleSize = 0.042; ps.subSize = 0.034;
    ps.boxX1 = 0.17; ps.boxX2 = 0.60; ps.boxY1 = 0.14; ps.boxY2 = 0.40;
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("c_syst", "", ps.w, ps.h);
    ApplyCanvasStyle(c, ps);
    if (logx) c->SetLogx();

    const int b1 = (xhi > xlo) ? nom->GetXaxis()->FindBin(xlo + 1e-6) : 1;
    const int b2 = (xhi > xlo) ? nom->GetXaxis()->FindBin(xhi - 1e-6) : nom->GetNbinsX();

    std::vector<TH1D *> rat;
    double ymin = 1.0, ymax = 1.0;
    for (int is = 0; is < kNSyst; ++is)
    {
        if (!ups[is] || !dns[is]) { rat.push_back(nullptr); rat.push_back(nullptr); continue; }
        TH1D *ru = Ratio(ups[is], nom, Form("r_up_%d", is));
        TH1D *rd = Ratio(dns[is], nom, Form("r_dn_%d", is));
        for (int b = b1; b <= b2; ++b)
        {
            ymin = std::min(ymin, std::min(ru->GetBinContent(b), rd->GetBinContent(b)));
            ymax = std::max(ymax, std::max(ru->GetBinContent(b), rd->GetBinContent(b)));
        }
        rat.push_back(ru); rat.push_back(rd);
    }
    double pad = std::max(0.02, 0.30 * (ymax - ymin));
    ymin = std::max(0.0, ymin - pad);
    ymax = ymax + 2.2 * pad;             // top: header + legend
    ymin -= 0.8 * (ymax - ymin);         // bottom: the info box
    if (ymin < 0) ymin = 0;

    TH1D *frame = (TH1D *)nom->Clone("frame_syst");
    frame->SetDirectory(nullptr);
    frame->Reset();
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle(xTitle);
    frame->GetYaxis()->SetTitle("variation / nominal");
    frame->GetXaxis()->SetTitleSize(ps.xTitleSize); frame->GetYaxis()->SetTitleSize(ps.yTitleSize);
    frame->GetXaxis()->SetLabelSize(ps.xLabelSize); frame->GetYaxis()->SetLabelSize(ps.yLabelSize);
    frame->GetXaxis()->SetTitleOffset(ps.xTitleOffset); frame->GetYaxis()->SetTitleOffset(ps.yTitleOffset);
    frame->SetMinimum(ymin); frame->SetMaximum(ymax);
    if (xhi > xlo) frame->GetXaxis()->SetRangeUser(xlo, xhi);
    frame->Draw("axis");

    TLine *one = new TLine(frame->GetXaxis()->GetBinLowEdge(b1), 1.0, frame->GetXaxis()->GetBinUpEdge(b2), 1.0);
    one->SetLineColor(kGray + 2); one->SetLineStyle(2); one->SetLineWidth(2);
    one->Draw();

    TLegend *leg = new TLegend(0.66, 0.15, 0.93, 0.41); // lower right, beside the info box
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(ps.font); leg->SetTextSize(0.028);
    std::vector<std::string> box;
    for (int is = 0; is < kNSyst; ++is)
    {
        TH1D *ru = rat[2 * is], *rd = rat[2 * is + 1];
        if (!ru) continue;
        ru->SetLineColor(kColor[is]); ru->SetLineWidth(2); ru->SetLineStyle(1);
        rd->SetLineColor(kColor[is]); rd->SetLineWidth(2); rd->SetLineStyle(2);
        ru->Draw("hist same"); rd->Draw("hist same");
        leg->AddEntry(ru, Form("%s Up", pOLhe::kLheSystNames[is]), "l");
        leg->AddEntry(rd, Form("%s Down", pOLhe::kLheSystNames[is]), "l");
        const Shift s = IntegralShift(nom, ups[is], dns[is]);
        double med, lo, hi;
        RatioSpread(nom, ups[is], dns[is], b1, b2, med, lo, hi);
        const int nos = OneSided(nom, ups[is], dns[is], b1, b2);
        box.push_back(Form("%s: integral %+.2f%% / %+.2f%%", pOLhe::kLheSystNames[is], 100 * s.up, 100 * s.dn));
        box.push_back(Form("   bins: med %.3f, 68%% [%.3f, %.3f]%s", med, lo, hi,
                           nos ? Form(", %d one-sided", nos) : ""));
    }
    leg->Draw();
    DrawHeader(ps, head, sub1, "LHE shape systematics");
    DrawInfoBox(ps, box);
    CMS_lumi(c, 13, 10);
    c->SaveAs((outNoExt + ".png").c_str());
    c->SaveAs((outNoExt + ".pdf").c_str());
    for (TH1D *r : rat) delete r;
    delete frame; delete one; delete leg; delete c;
}

// ---------------------------------------------------------------------------
// (2) integral shifts of `signal` vs rapidity bin
// ---------------------------------------------------------------------------
void DrawSummary(const std::vector<Shift> sh[2], // [charge][ybin]
                 const std::string &outNoExt, const char *syst, const char *binning,
                 const std::string &head)
{
    PlotStyle ps;
    ps.headerX = 0.47; ps.headerY = 0.90; ps.headerDy = 0.05;
    ps.titleSize = 0.042; ps.subSize = 0.034;
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("c_sum", "", ps.w, ps.h);
    ApplyCanvasStyle(c, ps);
    const int NY = 12;
    double ymin = 0, ymax = 0;
    TGraph *g[4];
    const int col[2] = {kRed + 1, kBlue + 1};
    const int mk[2][2] = {{20, 24}, {21, 25}}; // [charge][up/dn]: filled / open
    for (int ic = 0; ic < 2; ++ic)
        for (int ud = 0; ud < 2; ++ud)
        {
            TGraph *gg = new TGraph();
            for (int iy = 0; iy < NY; ++iy)
            {
                if (!sh[ic][iy].ok) continue;
                const double v = 100.0 * (ud == 0 ? sh[ic][iy].up : sh[ic][iy].dn);
                gg->SetPoint(gg->GetN(), iy + (ic == 0 ? -0.12 : 0.12), v);
                ymin = std::min(ymin, v); ymax = std::max(ymax, v);
            }
            gg->SetMarkerStyle(mk[ic][ud]); gg->SetMarkerColor(col[ic]); gg->SetLineColor(col[ic]);
            gg->SetMarkerSize(1.4);
            g[2 * ic + ud] = gg;
        }
    const double pad = std::max(1.0, 0.25 * (ymax - ymin));
    TH1D *frame = new TH1D("frame_sum", "", NY, -0.5, NY - 0.5);
    frame->SetDirectory(nullptr);
    for (int iy = 0; iy < NY; ++iy) frame->GetXaxis()->SetBinLabel(iy + 1, Form("y%d", iy));
    frame->GetXaxis()->SetTitle(Form("rapidity bin (%s)", binning));
    frame->GetYaxis()->SetTitle(Form("%s: integral shift of signal [%%]", syst));
    frame->GetXaxis()->SetLabelSize(0.045); frame->GetYaxis()->SetLabelSize(ps.yLabelSize);
    frame->GetXaxis()->SetTitleSize(ps.xTitleSize); frame->GetYaxis()->SetTitleSize(ps.yTitleSize);
    frame->GetYaxis()->SetTitleOffset(ps.yTitleOffset);
    frame->SetMinimum(ymin - pad - 0.6 * (ymax - ymin + 2 * pad)); // room for the legend
    frame->SetMaximum(ymax + pad + 0.45 * (ymax - ymin + 2 * pad)); // room for the header
    frame->Draw("axis");
    TLine *zero = new TLine(-0.5, 0, NY - 0.5, 0);
    zero->SetLineColor(kGray + 2); zero->SetLineStyle(2); zero->SetLineWidth(2); zero->Draw();
    TLegend *leg = new TLegend(0.17, 0.16, 0.60, 0.36);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(ps.font); leg->SetTextSize(0.03);
    const char *lab[4] = {"W^{+} Up", "W^{+} Down", "W^{-} Up", "W^{-} Down"};
    for (int i = 0; i < 4; ++i) { g[i]->Draw("P same"); leg->AddEntry(g[i], lab[i], "p"); }
    leg->Draw();
    DrawHeader(ps, head, Form("signal, %s binning", binning), Form("%s: (Up, Down) / nominal - 1", syst));
    CMS_lumi(c, 13, 10);
    c->SaveAs((outNoExt + ".png").c_str());
    c->SaveAs((outNoExt + ".pdf").c_str());
    for (int i = 0; i < 4; ++i) delete g[i];
    delete frame; delete zero; delete leg; delete c;
}

// ---------------------------------------------------------------------------
// (5) member overlays: the nominal + all 106 EPPS21 member templates of one
//     region on one canvas
// ---------------------------------------------------------------------------
struct MemberSet
{
    TH1D *nom = nullptr;          // member 0 = the nominal (== the Combine-input `signal`)
    std::vector<TH1D *> mem;      // members 1..106, k_s-scaled, W+ and W- samples summed
    double maxDev0 = -1.0;        // max |member 0 - Combine signal| / signal over the bins
    double intLo = 0, intHi = 0;  // min / max over the members of (integral / nominal - 1)
    void Clear() { delete nom; nom = nullptr; for (TH1D *h : mem) delete h; mem.clear(); }
};

// signal member m = k_Wp x (W+ sample twin row m) + k_Wm x (W- sample twin row m):
// the arithmetic of mtandmet.C (getScaled, then sum2) applied to every row of
// <stem>_epps21 (x = discriminant, y = member 0..106), so member 0 reproduces
// the Combine-input `signal` bit for bit (checked -> maxDev0)
bool BuildSignalMembers(TFile *fWp, TFile *fWm, double kWp, double kWm, const TString &stem,
                        const TH1 *combineSig, MemberSet &out)
{
    TH2D *a = (TH2D *)fWp->Get(stem + "_epps21");
    TH2D *b = (TH2D *)fWm->Get(stem + "_epps21");
    if (!a || !b || a->GetNbinsY() != 107 || b->GetNbinsY() != 107) return false;
    out.Clear();
    for (int m = 0; m < 107; ++m)
    {
        TH1D *h  = a->ProjectionX(Form("mem_%s_%d", stem.Data(), m), m + 1, m + 1);
        TH1D *hb = b->ProjectionX(Form("memb_%s_%d", stem.Data(), m), m + 1, m + 1);
        h->SetDirectory(nullptr);
        h->Scale(kWp);
        h->Add(hb, kWm);
        delete hb;
        h->SetFillStyle(0); h->SetFillColor(0);
        if (m == 0) out.nom = h; else out.mem.push_back(h);
    }
    const double i0 = Integ(out.nom);
    out.intLo = out.intHi = 0.0;
    for (TH1D *h : out.mem)
    {
        const double d = i0 > 0 ? Integ(h) / i0 - 1.0 : 0.0;
        out.intLo = std::min(out.intLo, d); out.intHi = std::max(out.intHi, d);
    }
    if (combineSig)
    {
        out.maxDev0 = 0.0;
        for (int i = 1; i <= out.nom->GetNbinsX(); ++i)
        {
            const double s = combineSig->GetBinContent(i);
            if (s > 0) out.maxDev0 = std::max(out.maxDev0, std::fabs(out.nom->GetBinContent(i) - s) / s);
        }
    }
    return true;
}

// top pad: nominal (black) + 106 members (grey, alpha); bottom pad: member /
// nominal for all of them + the LHAPDF Up/Down ratio (red solid / dashed)
void DrawMemberOverlay(const MemberSet &ms, const TH1 *up, const TH1 *dn, const std::string &outNoExt,
                       const char *xTitle, const std::string &head, const std::string &sub1,
                       double xlo, double xhi, bool logy = true)
{
    PlotStyle ps;
    ps.headerX = 0.42; ps.headerY = 0.90; ps.headerDy = 0.06;
    ps.titleSize = 0.048; ps.subSize = 0.038; ps.boxTextSize = 0.027;
    gStyle->SetOptStat(0);
    const double split = 0.32; // ratio-pad fraction of the canvas (pull-pad convention)
    TCanvas *c = new TCanvas("c_members", "", ps.w, (int)((double)ps.h * ps.pullCanvasScale + 0.5));
    TPad *pTop = new TPad("pTop_members", "", 0.0, split, 1.0, 1.0);
    TPad *pBot = new TPad("pBot_members", "", 0.0, 0.0, 1.0, split);
    pTop->SetTopMargin(ps.tm); pTop->SetBottomMargin(0.025); pTop->SetLeftMargin(ps.lm); pTop->SetRightMargin(ps.rm);
    pBot->SetTopMargin(0.04); pBot->SetBottomMargin(0.36); pBot->SetLeftMargin(ps.lm); pBot->SetRightMargin(ps.rm);
    for (TPad *p : {pTop, pBot}) { p->SetTicks(1, 1); p->SetFrameLineWidth(ps.FrameLineWidth); p->Draw(); }
    pTop->SetLogy(logy);

    const TH1D *nom = ms.nom;
    const int b1 = (xhi > xlo) ? nom->GetXaxis()->FindBin(xlo + 1e-6) : 1;
    const int b2 = (xhi > xlo) ? nom->GetXaxis()->FindBin(xhi - 1e-6) : nom->GetNbinsX();
    double peak = 0.0, floorPos = 1e300; // max and smallest positive content of the displayed bins
    for (int b = b1; b <= b2; ++b)
    {
        const double v = nom->GetBinContent(b);
        peak = std::max(peak, v);
        if (v > 0) floorPos = std::min(floorPos, v);
    }

    // ratios first (the legend on the top pad references the Up/Down curves)
    std::vector<TH1D *> rat;
    for (size_t m = 0; m < ms.mem.size(); ++m) rat.push_back(Ratio(ms.mem[m], nom, Form("r_mem_%zu", m)));
    TH1D *rUp = up ? Ratio(up, nom, "r_up_members") : nullptr;
    TH1D *rDn = dn ? Ratio(dn, nom, "r_dn_members") : nullptr;
    const double alpha = 0.25;

    // ---- top pad: the distributions ---------------------------------------
    pTop->cd();
    TH1D *frame = (TH1D *)nom->Clone("frame_members");
    frame->SetDirectory(nullptr); frame->Reset(); frame->SetTitle("");
    ApplyHistStyle(frame, ps, "", Form("Events / %g GeV", nom->GetXaxis()->GetBinWidth(1)));
    frame->GetXaxis()->SetLabelSize(0.0); frame->GetXaxis()->SetTitleSize(0.0);
    if (xhi > xlo) frame->GetXaxis()->SetRangeUser(xlo, xhi);
    if (logy)
    {
        // peak at 55% of the frame height (header + legend above it); the info
        // box goes lower-left, under the peak and the falling edge, and is kept
        // narrow so it stays clear of the descending tail
        const double ymin = std::max(0.5 * floorPos, 1e-3 * peak);
        frame->SetMinimum(ymin);
        frame->SetMaximum(ymin * std::pow(peak / ymin, 1.0 / 0.55));
        ps.boxX1 = 0.17; ps.boxX2 = 0.60; ps.boxY1 = 0.06; ps.boxY2 = 0.24;
    }
    else
    {
        frame->SetMinimum(0.0); frame->SetMaximum(1.8 * peak); // top 45%: header, legend, box
        ps.boxX1 = 0.42; ps.boxX2 = 0.93; ps.boxY1 = 0.33; ps.boxY2 = 0.49;
    }
    frame->Draw("axis");
    for (TH1D *h : ms.mem) { h->SetLineColorAlpha(kGray + 2, alpha); h->SetLineWidth(1); h->Draw("hist same"); }
    TH1D *nomC = (TH1D *)nom->Clone("nom_members_draw");
    nomC->SetDirectory(nullptr); nomC->SetLineColor(kBlack); nomC->SetLineWidth(3); nomC->SetFillStyle(0);
    nomC->Draw("hist same");
    TH1D *legMem = (TH1D *)nomC->Clone("leg_members"); // opaque proxy for the legend
    legMem->SetDirectory(nullptr); legMem->SetLineColor(kGray + 1); legMem->SetLineWidth(2);
    TLegend *leg = new TLegend(0.42, 0.53, 0.93, 0.75);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(ps.font); leg->SetTextSize(0.031);
    leg->AddEntry(nomC, "nominal (member 0)", "l");
    leg->AddEntry(legMem, "106 EPPS21 members", "l");
    if (rUp) { rUp->SetLineColor(kRed + 1); rUp->SetLineWidth(2); leg->AddEntry(rUp, "nPDF Up (LHAPDF, 68% CL)", "l"); }
    if (rDn) { rDn->SetLineColor(kRed + 1); rDn->SetLineWidth(2); rDn->SetLineStyle(2); leg->AddEntry(rDn, "nPDF Down", "l"); }
    leg->Draw();
    DrawHeader(ps, head, sub1, "EPPS21 member templates");
    std::vector<std::string> box; // integral shifts w.r.t. the nominal; short lines (see logy above)
    box.push_back("24 nuclear + 29 baseline pairs");
    box.push_back(Form("member integrals: %+.2f%%..%+.2f%%", 100 * ms.intLo, 100 * ms.intHi));
    const Shift s = IntegralShift(nom, up, dn);
    if (s.ok) box.push_back(Form("LHAPDF Up/Down: %+.2f%% / %+.2f%%", 100 * s.up, 100 * s.dn));
    DrawInfoBox(ps, box);
    CMS_lumi(pTop, 13, 10);
    pTop->RedrawAxis();

    // ---- zoom inset (log layout): the peak window on a linear y axis spanning
    // only a few % around the nominal, where the band IS visible -- at the scale
    // of the distribution a +-3% band is 1-2 px whatever the y axis does. The
    // window = the bins around the peak holding >= 85% of it (at least 5), its
    // y range = the nominal +- max(7%, 1.25 x the largest Up/Down deviation);
    // a dashed rectangle on the main plot marks it.
    TBox *zbox = nullptr; TH1D *zframe = nullptr, *upC = nullptr, *dnC = nullptr;
    if (logy)
    {
        int pk = b1;
        for (int b = b1; b <= b2; ++b) if (nom->GetBinContent(b) > nom->GetBinContent(pk)) pk = b;
        int z1 = pk, z2 = pk;
        while (z1 - 1 >= b1 && nom->GetBinContent(z1 - 1) >= 0.85 * peak) --z1;
        while (z2 + 1 <= b2 && nom->GetBinContent(z2 + 1) >= 0.85 * peak) ++z2;
        z1 = std::max(b1, std::min(z1, pk - 2)); z2 = std::min(b2, std::max(z2, pk + 2));
        double zmin = 1e300, zmax = 0.0, dev = 0.0;
        for (int b = z1; b <= z2; ++b)
        {
            zmin = std::min(zmin, nom->GetBinContent(b)); zmax = std::max(zmax, nom->GetBinContent(b));
            for (TH1D *r : {rUp, rDn}) if (r) dev = std::max(dev, std::fabs(r->GetBinContent(b) - 1.0));
            for (TH1D *r : rat) dev = std::max(dev, std::fabs(r->GetBinContent(b) - 1.0));
        }
        const double margin = std::max(0.07, 1.25 * dev);
        const double zxlo = nom->GetXaxis()->GetBinLowEdge(z1), zxhi = nom->GetXaxis()->GetBinUpEdge(z2);
        const double zlo = (1.0 - margin) * zmin, zhi = (1.0 + margin) * zmax;
        zbox = new TBox(zxlo, zlo, zxhi, zhi);
        zbox->SetFillStyle(0); zbox->SetLineColor(kBlue + 1); zbox->SetLineStyle(2); zbox->SetLineWidth(1);
        zbox->Draw("l");
        TPad *pIns = new TPad("pIns_members", "", 0.58, 0.30, 0.93, 0.52); // above the tail, below the legend
        pIns->SetLeftMargin(0.17); pIns->SetRightMargin(0.03); pIns->SetTopMargin(0.05); pIns->SetBottomMargin(0.22);
        pIns->SetTicks(1, 1); pIns->SetFrameLineWidth(1);
        pIns->Draw();
        pIns->cd();
        zframe = (TH1D *)frame->Clone("frame_members_zoom");
        zframe->SetDirectory(nullptr);
        zframe->GetXaxis()->SetRangeUser(zxlo, zxhi);
        zframe->SetMinimum(zlo); zframe->SetMaximum(zhi);
        zframe->GetXaxis()->SetLabelSize(0.10); zframe->GetYaxis()->SetLabelSize(0.10);
        zframe->GetXaxis()->SetTitleSize(0.0); zframe->GetYaxis()->SetTitleSize(0.0);
        zframe->GetXaxis()->SetNdivisions(505); zframe->GetYaxis()->SetNdivisions(404);
        zframe->Draw("axis");
        for (TH1D *h : ms.mem) h->Draw("hist same");
        if (up)
        {
            upC = (TH1D *)up->Clone("up_members_zoom"); upC->SetDirectory(nullptr);
            upC->SetFillStyle(0); upC->SetLineColor(kRed + 1); upC->SetLineWidth(2); upC->SetLineStyle(1);
            upC->Draw("hist same");
        }
        if (dn)
        {
            dnC = (TH1D *)dn->Clone("dn_members_zoom"); dnC->SetDirectory(nullptr);
            dnC->SetFillStyle(0); dnC->SetLineColor(kRed + 1); dnC->SetLineWidth(2); dnC->SetLineStyle(2);
            dnC->Draw("hist same");
        }
        nomC->Draw("hist same");
        TLatex zl; zl.SetNDC(); zl.SetTextFont(ps.font); zl.SetTextSize(0.11); zl.SetTextColor(kBlue + 1);
        zl.DrawLatex(0.21, 0.82, "zoom");
        pIns->RedrawAxis();
        pTop->cd();
    }

    // ---- bottom pad: member / nominal ----------------------------------------
    pBot->cd();
    const double sf = (1.0 - split) / split; // fonts are pad-relative
    TH1D *rframe = (TH1D *)frame->Clone("frame_members_ratio");
    rframe->SetDirectory(nullptr);
    rframe->GetXaxis()->SetTitle(xTitle); rframe->GetYaxis()->SetTitle("member / nominal");
    rframe->GetXaxis()->SetTitleSize(ps.xTitleSize * sf); rframe->GetYaxis()->SetTitleSize(ps.yTitleSize * sf);
    rframe->GetXaxis()->SetLabelSize(ps.xLabelSize * sf); rframe->GetYaxis()->SetLabelSize(ps.yLabelSize * sf);
    rframe->GetXaxis()->SetTitleOffset(1.0); rframe->GetYaxis()->SetTitleOffset(ps.yTitleOffset / sf);
    rframe->GetYaxis()->SetNdivisions(505);
    // robust y range: the members and Up/Down in bins holding >= 0.5% of the
    // peak (a few-event tail bin with sign-flipping weights would set the scale)
    double lo = 1.0, hi = 1.0;
    for (int b = b1; b <= b2; ++b)
    {
        if (nom->GetBinContent(b) < 0.005 * peak) continue;
        for (TH1D *r : rat) { lo = std::min(lo, r->GetBinContent(b)); hi = std::max(hi, r->GetBinContent(b)); }
        for (TH1D *r : {rUp, rDn})
            if (r) { lo = std::min(lo, r->GetBinContent(b)); hi = std::max(hi, r->GetBinContent(b)); }
    }
    const double half = std::min(0.5, std::max(0.03, 1.3 * std::max(hi - 1.0, 1.0 - lo)));
    rframe->SetMinimum(1.0 - half); rframe->SetMaximum(1.0 + half);
    rframe->Draw("axis");
    // "][" = no vertical line at the first/last drawn bin edge (108 of them would
    // otherwise pile up into a dark bar on the left frame edge)
    for (TH1D *r : rat) { r->SetLineColorAlpha(kGray + 2, alpha); r->SetLineWidth(1); r->Draw("hist ][ same"); }
    TLine *one = new TLine(frame->GetXaxis()->GetBinLowEdge(b1), 1.0, frame->GetXaxis()->GetBinUpEdge(b2), 1.0);
    one->SetLineColor(kBlack); one->SetLineStyle(2); one->SetLineWidth(1);
    one->Draw();
    if (rUp) rUp->Draw("hist ][ same");
    if (rDn) rDn->Draw("hist ][ same");
    pBot->RedrawAxis();

    c->SaveAs((outNoExt + ".png").c_str());
    c->SaveAs((outNoExt + ".pdf").c_str());
    for (TH1D *r : rat) delete r;
    delete rUp; delete rDn; delete one; delete leg; delete legMem; delete nomC; delete frame; delete rframe;
    delete zbox; delete zframe; delete upC; delete dnC;
    delete c;
}

// ---------------------------------------------------------------------------
// (3)/(4) helpers on the member twins
// ---------------------------------------------------------------------------
struct Incl { double nuc_up = 0, nuc_dn = 0, base_up = 0, base_dn = 0, tot_up = 0, tot_dn = 0,
                     sc_up = 0, sc_dn = 0, as_up = 0, as_dn = 0; bool ok = false; };

// asymmetric Hessian (LHAPDF formula) over the member pairs [lo,hi] of I[],
// rescaled from the set's 90% CL to 68.27% exactly as LHAPDF does:
// sqrt(chi2_quantile(0.6827, 1) / chi2_quantile(0.90, 1)) = 0.607957 = 1/1.64485
// (the "1.645" quoted in the notes is this number rounded).
const double kCL90to68 = std::sqrt(TMath::ChisquareQuantile(0.682689492137, 1) /
                                   TMath::ChisquareQuantile(0.90, 1));
void HessOn(const std::vector<double> &I, int lo, int hi, double &up, double &dn)
{
    pOLhe::Block b{lo, hi, pOLhe::kHessian, "", "", "", ""};
    auto dev = [&](int m) { return I[0] > 0 ? I[m] / I[0] - 1.0 : 0.0; };
    const pOLhe::HessRes h = pOLhe::Hessian(b, dev);
    up = h.up * kCL90to68; dn = h.dn * kCL90to68;
}

// combine member integrals (epps21 107, scale 9, alphas 5) into the relative shifts
Incl Combine(const std::vector<double> &E, const std::vector<double> &S, const std::vector<double> &A)
{
    Incl r;
    if (E.size() != 107 || S.size() != 9 || A.size() != 5 || E[0] <= 0) return r;
    r.ok = true;
    HessOn(E, 1, 48, r.nuc_up, r.nuc_dn);
    HessOn(E, 49, 106, r.base_up, r.base_dn);
    HessOn(E, 1, 106, r.tot_up, r.tot_dn);
    double mx = S[0], mn = S[0];
    for (int m = 0; m < 9; ++m) { mx = std::max(mx, S[m]); mn = std::min(mn, S[m]); }
    r.sc_up = mx / S[0] - 1.0; r.sc_dn = mn / S[0] - 1.0;
    r.as_up = A[3] / A[0] - 1.0; r.as_dn = A[2] / A[0] - 1.0;
    return r;
}

// add the member integrals of twin `name` (over x bins 1..nx) into acc[m]
bool AddMembers(TFile *f, const TString &name, std::vector<double> &acc, int nmem)
{
    TH2D *h = (TH2D *)f->Get(name);
    if (!h) return false;
    if (acc.empty()) acc.assign(nmem, 0.0);
    const int nx = h->GetNbinsX();
    for (int m = 0; m < nmem; ++m) acc[m] += h->Integral(1, nx, m + 1, m + 1);
    return true;
}

// (4): max |LHAPDF Up/Down - Hessian Up/Down| / nominal over the bins of one template
double ClosureOne(TFile *f, const TString &stem)
{
    TH1D *nom = (TH1D *)f->Get(stem);
    TH2D *tw  = (TH2D *)f->Get(stem + "_epps21");
    TH1D *up  = (TH1D *)f->Get(stem + "_nPDFUp");
    TH1D *dn  = (TH1D *)f->Get(stem + "_nPDFDown");
    if (!nom || !tw || !up || !dn) return -1.0;
    double worst = 0.0;
    std::vector<double> I(107);
    for (int ix = 1; ix <= nom->GetNbinsX(); ++ix)
    {
        for (int m = 0; m < 107; ++m) I[m] = tw->GetBinContent(ix, m + 1);
        if (I[0] <= 0) continue;
        double hu, hd;
        HessOn(I, 1, 106, hu, hd);
        const double eu = std::fabs(up->GetBinContent(ix) - std::max(I[0] * (1 + hu), 0.0)) / I[0];
        const double ed = std::fabs(dn->GetBinContent(ix) - std::max(I[0] * (1 - hd), 0.0)) / I[0];
        worst = std::max(worst, std::max(eu, ed));
    }
    return worst;
}

// (c): the all-events per-index S_i/S_0 - 1 (fraction) of lhe_weights.txt, per column
bool ReadReference(const char *path, std::vector<double> &wp, std::vector<double> &wm, std::vector<double> &dy)
{
    std::ifstream in(path);
    if (!in) return false;
    wp.assign(217, 0); wm.assign(217, 0); dy.assign(217, 0);
    std::string line;
    int nrow = 0;
    while (std::getline(in, line))
    {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        int idx; double a, ra, b, rb, c, rc;
        if (!(ss >> idx >> a >> ra >> b >> rb >> c >> rc)) continue;
        if (idx < 0 || idx > 216) continue;
        wp[idx] = a / 100.0; wm[idx] = b / 100.0; dy[idx] = c / 100.0;
        ++nrow;
    }
    return nrow == 217;
}

Incl ReferenceIncl(const std::vector<double> &rel)
{
    // rebuild member-ordered relative integrals: E[m] = 1 + rel[idx(m)] for the
    // nuclear members (1..48 <- idx 111..158); the baseline PRODUCTS have no
    // all-events counterpart in the txt -> left at 1 (base term = 0)
    std::vector<double> E(107, 1.0), S(9, 1.0), A(5, 1.0);
    for (int m = 1; m <= 48; ++m) E[m] = 1.0 + rel[110 + m];
    for (int m = 0; m <= 8; ++m) S[m] = 1.0 + rel[m];
    A[2] = 1.0 + rel[104]; A[3] = 1.0 + rel[105]; A[1] = 1.0 + rel[103]; A[4] = 1.0 + rel[106];
    return Combine(E, S, A);
}

std::string Fmt(const Incl &r)
{
    if (!r.ok) return "   (n/a)";
    return Form("nPDF %+.2f/%+.2f%% (nucl %+.2f/%+.2f, base %+.2f/%+.2f)  qcdScale %+.2f/%+.2f%%  alphaS %+.2f/%+.2f%%",
                100 * r.tot_up, -100 * r.tot_dn, 100 * r.nuc_up, -100 * r.nuc_dn, 100 * r.base_up, -100 * r.base_dn,
                100 * r.sc_up, 100 * r.sc_dn, 100 * r.as_up, 100 * r.as_dn);
}

} // namespace

void syst_shapes(const char *disc = "leppt_mt40")
{
    gROOT->SetBatch(kTRUE);
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);
    TString dsuf, discLabel;
    if (!pODisc::Spec(disc, dsuf, discLabel)) return;
    const bool isMET = (TString(disc) == "met");
    const char *xTitleW = isMET ? "PF MET (GeV)" : "Lepton p_{T} (GeV)";
    const double xlo = isMET ? 0.0 : 24.0, xhi = isMET ? -1.0 : 100.0; // pT floor convention
    const char *stem = isMET ? "h_met" : "h_leppt_mt40";               // skim twin stems

    // all-events reference (mu-flavour files; the only rows lhe_weights.txt has)
    std::vector<double> refWp, refWm, refDY;
    const bool haveRef = ReadReference("../skim/output/lhe_weights.txt", refWp, refWm, refDY);
    if (!haveRef) std::cerr << "[WARN] ../skim/output/lhe_weights.txt not readable -> no all-events reference\n";

    double worstClosure = 0.0, worstMember0 = -1.0;
    int nPlots = 0, nOneSided = 0, nTemplates = 0, nMemberPlots = 0;

    for (int fl = 0; fl < 2; ++fl)
    {
        const bool isElec = (fl == 1);
        const std::string plotsBase = isElec ? "./plots/Elec" : "./plots";
        const std::string outDir = plotsBase + "/syst_shapes/" + disc;
        gSystem->mkdir((outDir + "/perbin").c_str(), kTRUE);
        const char *flav = isElec ? "ele" : "mu";
        const std::string headW = isElec ? "W #rightarrow e #nu" : "W #rightarrow #mu #nu";
        const std::string headZ = isElec ? "Z #rightarrow e e" : "Z #rightarrow #mu #mu";
        const std::string wPrefix = isElec ? "WToElecNu_pO_PFMet" : "WToMuNu_pO_PFMet";
        const std::string zPrefix = isElec ? "ZToEE_pO2025" : "ZToMuMu_pO2025";

        std::cout << "\n===================== " << flav << ", " << disc << " =====================\n";

        // ---- (1)+(2) per-region plots from the Combine inputs ------------------
        struct Input { std::string file; std::string head; const char *xTitle; double lo, hi; bool isZ; };
        const Input inputs[2] = {
            {plotsBase + "/combine_input_W" + dsuf.Data() + ".root", headW, xTitleW, xlo, xhi, false},
            {plotsBase + "/combine_input_Z.root", headZ, "m_{ll} (GeV)", 0.0, -1.0, true}};
        std::map<std::string, std::vector<Shift>> sigShift; // key "lab|fb_Wp_<syst>" -> 12
        // (b) sum of the per-bin Up/Down over the 24 lab SR regions, per process
        std::map<std::string, double> sumNom, sumUp[3], sumDn[3];

        for (const Input &in : inputs)
        {
            TFile *f = TFile::Open(in.file.c_str(), "READ");
            if (!f || f->IsZombie()) { std::cerr << "[WARN] cannot open " << in.file << "\n"; continue; }
            TIter it(f->GetListOfKeys());
            while (TKey *k = (TKey *)it())
            {
                if (TString(k->GetClassName()).BeginsWith("TDirectory") == kFALSE) continue;
                const TString region = k->GetName();
                TDirectory *d = f->GetDirectory(region);
                if (!d) continue;
                // MC processes = those with an Up template
                std::vector<std::string> procs;
                TIter kt(d->GetListOfKeys());
                while (TKey *kk = (TKey *)kt())
                {
                    TString n = kk->GetName();
                    if (n.EndsWith("_nPDFUp"))
                    {
                        // NB TString::operator()(start,len) returns a TSubString whose
                        // Data() points into the FULL string -- copy into a TString first
                        const TString base = n(0, n.Length() - 7);
                        procs.push_back(std::string(base.Data()));
                    }
                }
                if (procs.empty()) continue; // CR dirs, data-only
                std::sort(procs.begin(), procs.end());
                const bool isSR = region.BeginsWith("Wp_") || region.BeginsWith("Wm_");
                const bool isLab = region.Contains("_lab_");
                const bool isFb  = region.Contains("_fb_");
                int iy = -1, ic = -1;
                if (isSR && (isLab || isFb))
                {
                    ic = region.BeginsWith("Wp_") ? 0 : 1;
                    iy = TString(region(region.Last('y') + 1, region.Length())).Atoi();
                }
                for (const std::string &p : procs)
                {
                    TH1 *nom = (TH1 *)d->Get(p.c_str());
                    if (!nom) continue;
                    std::vector<TH1 *> ups, dns;
                    for (int is = 0; is < kNSyst; ++is)
                    {
                        ups.push_back((TH1 *)d->Get(Form("%s_%sUp", p.c_str(), pOLhe::kLheSystNames[is])));
                        dns.push_back((TH1 *)d->Get(Form("%s_%sDown", p.c_str(), pOLhe::kLheSystNames[is])));
                    }
                    const int b1 = (in.hi > in.lo) ? nom->GetXaxis()->FindBin(in.lo + 1e-6) : 1;
                    const int b2 = (in.hi > in.lo) ? nom->GetXaxis()->FindBin(in.hi - 1e-6) : nom->GetNbinsX();
                    ++nTemplates;
                    std::string line = Form("[SHAPE] %-4s %-14s %-7s nominal %10.2f", flav, region.Data(), p.c_str(), Integ(nom));
                    for (int is = 0; is < kNSyst; ++is)
                    {
                        if (!ups[is] || !dns[is]) continue;
                        const Shift s = IntegralShift(nom, ups[is], dns[is]);
                        const int nos = OneSided(nom, ups[is], dns[is], b1, b2);
                        nOneSided += nos;
                        line += Form("  %s %+.2f/%+.2f%%%s", pOLhe::kLheSystNames[is], 100 * s.up, 100 * s.dn,
                                     nos ? Form(" (%d 1-sided)", nos) : "");
                        if (isSR && iy >= 0 && p == "signal")
                            sigShift[Form("%s_%s_%s", isLab ? "lab" : "fb", ic == 0 ? "Wp" : "Wm",
                                          pOLhe::kLheSystNames[is])].resize(12),
                            sigShift[Form("%s_%s_%s", isLab ? "lab" : "fb", ic == 0 ? "Wp" : "Wm",
                                          pOLhe::kLheSystNames[is])][iy] = s;
                        if (isSR && isLab)
                        {
                            sumUp[is][p] += Integ(ups[is]);
                            sumDn[is][p] += Integ(dns[is]);
                        }
                    }
                    if (isSR && isLab) sumNom[p] += Integ(nom);
                    std::cout << line << "\n";
                    DrawRatioSet(nom, ups, dns, outDir + "/perbin/" + region.Data() + "_" + p, in.xTitle, in.head,
                                 Form("%s, %s", region.Data(), p.c_str()), in.lo, in.hi);
                    ++nPlots;
                }
            }
            f->Close(); delete f;
        }

        // (2) summaries: signal vs y bin
        for (int is = 0; is < kNSyst; ++is)
            for (int ib = 0; ib < 2; ++ib)
            {
                const char *B = ib == 0 ? "lab" : "fb";
                std::vector<Shift> sh[2];
                sh[0] = sigShift[Form("%s_Wp_%s", B, pOLhe::kLheSystNames[is])];
                sh[1] = sigShift[Form("%s_Wm_%s", B, pOLhe::kLheSystNames[is])];
                if (sh[0].size() != 12 || sh[1].size() != 12) continue;
                DrawSummary(sh, outDir + "/summary_" + pOLhe::kLheSystNames[is] + "_" + B,
                            pOLhe::kLheSystNames[is], B, headW);
                ++nPlots;
            }

        // ---- (5) member overlays: nominal + all 106 EPPS21 members per region ----
        {
            gSystem->mkdir((outDir + "/members").c_str(), kTRUE);
            const std::string fWpName = "../skim/rootfile/" + wPrefix + "_Wp_hist.root";
            const std::string fWmName = "../skim/rootfile/" + wPrefix + "_Wm_hist.root";
            TFile *fWp = TFile::Open(fWpName.c_str(), "READ");
            TFile *fWm = TFile::Open(fWmName.c_str(), "READ");
            TFile *fC  = TFile::Open(inputs[0].file.c_str(), "READ");
            if (!fWp || fWp->IsZombie() || !fWm || fWm->IsZombie() || !fC || fC->IsZombie())
                std::cerr << "[WARN] member overlays skipped: cannot open " << fWpName << " / " << fWmName
                          << " / " << inputs[0].file << "\n";
            else
            {
                const double kWp = pONorm::MCScale(isElec ? "Wp_ele" : "Wp_mu");
                const double kWm = pONorm::MCScale(isElec ? "Wm_ele" : "Wm_mu");
                for (int ic = 0; ic < 2; ++ic)
                    for (int ib = 0; ib < 2; ++ib)
                        for (int iy = 0; iy < 12; ++iy)
                        {
                            const char *chg = ic == 0 ? "Wp" : "Wm";
                            const TString region = Form("%s_%s_y%d", chg, ib == 0 ? "lab" : "fb", iy);
                            const TString hstem  = Form("%s_%s_y%d%s", stem, chg, iy, ib == 0 ? "" : "_FB");
                            TDirectory *d = fC->GetDirectory(region);
                            TH1 *sig = d ? (TH1 *)d->Get("signal") : nullptr;
                            TH1 *up  = d ? (TH1 *)d->Get("signal_nPDFUp") : nullptr;
                            TH1 *dn  = d ? (TH1 *)d->Get("signal_nPDFDown") : nullptr;
                            MemberSet ms;
                            if (!BuildSignalMembers(fWp, fWm, kWp, kWm, hstem, sig, ms))
                            {
                                std::cerr << "[WARN] " << region << ": " << hstem << "_epps21 missing in a W skim file"
                                          << " -> no member overlay\n";
                                continue;
                            }
                            const Shift s = IntegralShift(ms.nom, up, dn);
                            std::cout << Form("[MEMBER] %-4s %-10s signal: member 0 vs Combine signal max|diff|/signal %.1e;"
                                              " member integrals %+.2f%% .. %+.2f%%; LHAPDF Up/Down %+.2f%% / %+.2f%%",
                                              flav, region.Data(), ms.maxDev0, 100 * ms.intLo, 100 * ms.intHi,
                                              100 * s.up, 100 * s.dn) << "\n";
                            worstMember0 = std::max(worstMember0, ms.maxDev0);
                            DrawMemberOverlay(ms, up, dn, outDir + "/members/" + region.Data() + "_signal", xTitleW, headW,
                                              Form("W^{%s}, %s binning, y%d, signal", ic == 0 ? "+" : "-",
                                                   ib == 0 ? "lab" : "fb", iy), xlo, xhi);
                            ++nPlots; ++nMemberPlots;
                            ms.Clear();
                        }
            }
            for (TFile *f : {fWp, fWm, fC}) if (f) { f->Close(); delete f; }
        }

        // ---- (3)(a) member-level inclusive from the skim twins ------------------
        // per sample: sum the 24 lab templates' twins; per process: k_s-weighted
        struct Samp { const char *tok; const char *label; const char *proc; };
        const Samp samps[6] = {{"Wp", isElec ? "Wp_ele" : "Wp_mu", "signal"},
                               {"Wm", isElec ? "Wm_ele" : "Wm_mu", "signal"},
                               {"DY", isElec ? "DYee" : "DYmu", "z"},
                               {"DYtau", "DYtau", "ztau"},
                               {"Wptau", "Wp_tau", "wtau"},
                               {"Wmtau", "Wm_tau", "wtau"}};
        std::map<std::string, std::vector<double>> procE, procS, procA; // k_s-weighted, per process
        std::cout << "\n[INCL] " << flav << " " << disc << ": inclusive (24 lab templates) -- (a) member-level"
                  << " combination of the skim twins vs (c) the all-events reference of lhe_weights.txt (mu rows)\n";
        for (const Samp &sp : samps)
        {
            const std::string fn = "../skim/rootfile/" + wPrefix + "_" + sp.tok + "_hist.root";
            TFile *f = TFile::Open(fn.c_str(), "READ");
            if (!f || f->IsZombie()) { std::cerr << "[WARN] cannot open " << fn << "\n"; continue; }
            std::vector<double> E, S, A;
            int nfound = 0;
            for (int ic = 0; ic < 2; ++ic)
                for (int iy = 0; iy < 12; ++iy)
                {
                    const TString nm = Form("%s_%s_y%d", stem, ic == 0 ? "Wp" : "Wm", iy);
                    const bool ok = AddMembers(f, nm + "_epps21", E, 107) && AddMembers(f, nm + "_scale", S, 9) &&
                                    AddMembers(f, nm + "_alphas", A, 5);
                    if (ok) ++nfound;
                    const double cl = ClosureOne(f, nm);
                    if (cl >= 0) worstClosure = std::max(worstClosure, cl);
                }
            if (nfound != 24) std::cerr << "[WARN] " << fn << ": only " << nfound << "/24 twins found\n";
            const Incl a = Combine(E, S, A);
            const std::vector<double> *ref = (TString(sp.tok).BeginsWith("Wp")) ? &refWp
                                             : (TString(sp.tok).BeginsWith("Wm")) ? &refWm : &refDY;
            const Incl c = haveRef ? ReferenceIncl(*ref) : Incl();
            std::cout << "[INCL]   sample " << Form("%-6s", sp.tok) << " (a) selected " << Fmt(a) << "\n";
            std::cout << "[INCL]   " << "       " << " (c) all-evts " << Fmt(c) << "   [ref row: "
                      << (ref == &refWp ? "Wp_mu" : ref == &refWm ? "Wm_mu" : "DY_mu") << "; no baseline products in the txt]\n";
            // per process, k_s-weighted
            const double ks = pONorm::MCScale(sp.label);
            std::vector<double> &pe = procE[sp.proc], &ps = procS[sp.proc], &pa = procA[sp.proc];
            if (pe.empty()) { pe.assign(107, 0); ps.assign(9, 0); pa.assign(5, 0); }
            for (size_t m = 0; m < E.size(); ++m) pe[m] += ks * E[m];
            for (size_t m = 0; m < S.size(); ++m) ps[m] += ks * S[m];
            for (size_t m = 0; m < A.size(); ++m) pa[m] += ks * A[m];
            f->Close(); delete f;
        }
        // Z peak: the single hMass template per sample (closure only, no reference)
        for (const Samp &sp : samps)
        {
            const std::string fn = "../skim/rootfile/" + zPrefix + "_" + sp.tok + "_MC_hist.root";
            TFile *f = TFile::Open(fn.c_str(), "READ");
            if (!f || f->IsZombie()) continue;
            const double cl = ClosureOne(f, "hMass");
            if (cl >= 0) worstClosure = std::max(worstClosure, cl);
            f->Close(); delete f;
        }

        // (a) vs (b) per process
        std::cout << "[INCL] " << flav << " " << disc << ": per Combine process (24 lab SR regions): (a) combination of the"
                  << " summed members = the true inclusive uncertainty; (b) sum of the per-bin Up/Down templates ="
                  << " what ONE collapsed nuisance implies; (b)/(a) = the collapse inflation\n";
        for (const char *p : {"signal", "z", "ztau", "wtau"})
        {
            if (!procE.count(p) || !sumNom.count(p) || sumNom[p] <= 0) continue;
            const Incl a = Combine(procE[p], procS[p], procA[p]);
            const double n0 = sumNom[p];
            std::cout << "[INCL]   " << Form("%-7s", p) << " (a) " << Fmt(a) << "\n";
            std::cout << "[INCL]   " << "        (b) "
                      << Form("nPDF %+.2f/%+.2f%%                              qcdScale %+.2f/%+.2f%%  alphaS %+.2f/%+.2f%%",
                              100 * (sumUp[0][p] / n0 - 1), 100 * (sumDn[0][p] / n0 - 1),
                              100 * (sumUp[1][p] / n0 - 1), 100 * (sumDn[1][p] / n0 - 1),
                              100 * (sumUp[2][p] / n0 - 1), 100 * (sumDn[2][p] / n0 - 1)) << "\n";
            if (a.ok && a.tot_up > 0 && a.tot_dn > 0)
                std::cout << "[INCL]   " << "        (b)/(a): nPDF up " << Form("%.2f", (sumUp[0][p] / n0 - 1) / a.tot_up)
                          << ", down " << Form("%.2f", (1 - sumDn[0][p] / n0) / a.tot_dn)
                          << "; qcdScale up " << Form("%.2f", (sumUp[1][p] / n0 - 1) / a.sc_up)
                          << ", down " << Form("%.2f", (sumDn[1][p] / n0 - 1) / a.sc_dn) << "\n";
        }
    }

    std::cout << "\n[CLOSURE] LHAPDF (stored <h>_nPDFUp/Down) vs pOLhe::Hessian recomputed from the _epps21 twins:"
              << " max |difference| / nominal over all W lab templates + hMass, both flavours = "
              << Form("%.2e", worstClosure) << (worstClosure < 1e-6 ? "  (identical)" : "  (CHECK)") << "\n";
    std::cout << "[CLOSURE] member 0 of the rebuilt signal members vs the Combine-input `signal`, max |difference| /"
              << " signal over all W regions, both flavours = " << Form("%.2e", worstMember0)
              << (worstMember0 < 0 ? "  (no overlays built)" : worstMember0 < 1e-9 ? "  (identical)" : "  (CHECK)") << "\n";
    std::cout << "[SUMMARY] " << disc << ": " << nTemplates << " MC templates checked, " << nPlots
              << " plots written (of which " << nMemberPlots << " member overlays), one-sided (Up and Down on the"
              << " same side) bins in total: " << nOneSided << " -> plots[/Elec]/syst_shapes/" << disc << "/\n";
    std::cout << "[OK] syst_shapes(" << disc << ")\n";
}
