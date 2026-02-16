#include "CMS_lumi.C"

// -----------------------------
// 1) Generic style/options blob
// -----------------------------
struct PlotStyle
{
    // canvas
    int w = 800, h = 800;
    double lm = 0.14, rm = 0.05, bm = 0.12, tm = 0.07;
    bool ticks = true;

    // fonts
    int font = 42;
    double titleSize = 0.050;
    double subSize = 0.040;
    double boxTextSize = 0.034;

    // header positions (NDC)
    double headerX = 0.6;
    double headerY = 0.93;  // main title top
    double headerDy = 0.06; // spacing between lines

    // box (NDC)
    double boxX1 = 0.55, boxY1 = 0.42, boxX2 = 0.93, boxY2 = 0.88;

    // axis sizes
    double xTitleSize = 0.045, yTitleSize = 0.045;
    double xLabelSize = 0.040, yLabelSize = 0.040;
    double xTitleOffset = 1.10, yTitleOffset = 1.35;
    double FrameLineWidth = 3;

    // draw options
    std::string drawOpt = "hist"; // e.g. "E", "hist", "E1", etc.
    std::string drawOptGraph = "AP";
    bool logy = false;

    // stat box
    bool showStats = false;
};

// -----------------------------
// 2) Basic helpers
// -----------------------------
static void ApplyCanvasStyle(TCanvas *c, const PlotStyle &ps)
{
    c->SetLeftMargin(ps.lm);
    c->SetRightMargin(ps.rm);
    c->SetBottomMargin(ps.bm);
    c->SetTopMargin(ps.tm);
    c->SetFrameLineWidth(ps.FrameLineWidth);
    if (ps.ticks)
        c->SetTicks(1, 1);
    c->SetLogy(ps.logy);
}

static void ApplyHistStyle(TH1 *h, const PlotStyle &ps,
                           const std::string &xtitle,
                           const std::string &ytitle)
{
    h->SetTitle("");

    h->GetXaxis()->SetTitle(xtitle.c_str());
    h->GetYaxis()->SetTitle(ytitle.c_str());

    h->GetXaxis()->SetTitleFont(ps.font);
    h->GetYaxis()->SetTitleFont(ps.font);
    h->GetXaxis()->SetLabelFont(ps.font);
    h->GetYaxis()->SetLabelFont(ps.font);

    h->GetXaxis()->SetTitleSize(ps.xTitleSize);
    h->GetYaxis()->SetTitleSize(ps.yTitleSize);
    h->GetXaxis()->SetLabelSize(ps.xLabelSize);
    h->GetYaxis()->SetLabelSize(ps.yLabelSize);

    h->GetXaxis()->SetTitleOffset(ps.xTitleOffset);
    h->GetYaxis()->SetTitleOffset(ps.yTitleOffset);
}

static void DrawHeader(const PlotStyle &ps,
                       const std::string &mainTitle,
                       const std::string &sub1,
                       const std::string &sub2)
{
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextFont(ps.font);
    lat.SetTextAlign(13); // left-top

    lat.SetTextSize(ps.titleSize);
    lat.DrawLatex(ps.headerX, ps.headerY, mainTitle.c_str());

    lat.SetTextSize(ps.subSize);
    if (!sub1.empty())
        lat.DrawLatex(ps.headerX, ps.headerY - ps.headerDy, sub1.c_str());
    if (!sub2.empty())
        lat.DrawLatex(ps.headerX, ps.headerY - 2.0 * ps.headerDy, sub2.c_str());
}

static void DrawInfoBox(const PlotStyle &ps,
                        const std::vector<std::string> &lines)
{
    if (lines.empty())
        return;

    TPaveText *p = new TPaveText(ps.boxX1, ps.boxY1, ps.boxX2, ps.boxY2, "NDC");
    p->SetFillStyle(0);
    p->SetBorderSize(0);
    p->SetTextFont(ps.font);
    p->SetTextSize(ps.boxTextSize);
    p->SetTextAlign(12);
    for (const auto &s : lines)
        p->AddText(s.c_str());
    p->Draw();
}

// ---------------------------------------------------------
// 3) The ONE generic plot function (works for any TH1*)
//    - h can be something you read from file
//    - optional "tuner" lets you do plot-specific tweaks
// ---------------------------------------------------------
using PlotTuner = std::function<void(TCanvas *, TH1 *)>;
using GraphTuner = std::function<void(TCanvas*, TGraphErrors*)>;


/* Example of how to use a lambda version of type defined Plot Tuner
PlotTuner mtTuner = [](TCanvas* c, TH1* h) {
    c->SetLogy();
    h->GetXaxis()->SetRangeUser(40, 200);
};
*/

// The following main function takes genearal input of
// Saving Path
// Axis Title
// ... etc check it yourself

static void SaveNicePlot1D(TH1 *h,
                           const std::string &outPathNoExt, // e.g. "./mT/foo"
                           const std::string &xTitle,
                           const std::string &yTitle,
                           const std::string &mainTitle,
                           const std::string &subTitle1,
                           const std::string &subTitle2,
                           const std::vector<std::string> &boxLines,
                           const PlotStyle &ps = PlotStyle(),
                           PlotTuner tuner = nullptr)
{
    if (!h)
        return;

    gStyle->SetOptStat(ps.showStats ? 1110 : 0);

    TCanvas *c = new TCanvas(Form("c_%s", h->GetName()), "", ps.w, ps.h);
    ApplyCanvasStyle(c, ps);

    c->cd();

    ApplyHistStyle(h, ps, xTitle, yTitle);

    h->Draw(ps.drawOpt.c_str());

    DrawHeader(ps, mainTitle, subTitle1, subTitle2);
    DrawInfoBox(ps, boxLines);

    if (tuner)
        tuner(c, h); // e.g. change ranges, add extra lines, etc.

    // ensure parent directory exists
    // outPathNoExt could include dirs, so you can mkdir manually outside too.

    CMS_lumi(c, 13, 10);
    c->Modified();
    c->Update();

    c->SaveAs((outPathNoExt + ".png").c_str());
    c->SaveAs((outPathNoExt + ".pdf").c_str());

    delete c;
}

static std::string RangeLabel(const char *var, int iy, const double *edges, int nY)
{
    const double a = edges[iy];
    const double b = edges[iy + 1];

    // Use [a, b) except last bin [a, b]
    if (iy == nY - 1)
        return Form("%s #in [%.1f, %.1f]", var, a, b);
    return Form("%s #in [%.1f, %.1f)", var, a, b);
}

static void ApplyGraphStyle(TGraphErrors* g, const PlotStyle& ps,
                            const std::string& xtitle,
                            const std::string& ytitle)
{
    g->SetTitle("");
    g->GetXaxis()->SetTitle(xtitle.c_str());
    g->GetYaxis()->SetTitle(ytitle.c_str());

    g->GetXaxis()->SetTitleFont(ps.font);
    g->GetYaxis()->SetTitleFont(ps.font);
    g->GetXaxis()->SetLabelFont(ps.font);
    g->GetYaxis()->SetLabelFont(ps.font);

    g->GetXaxis()->SetTitleSize(ps.xTitleSize);
    g->GetYaxis()->SetTitleSize(ps.yTitleSize);
    g->GetXaxis()->SetLabelSize(ps.xLabelSize);
    g->GetYaxis()->SetLabelSize(ps.yLabelSize);

    g->GetXaxis()->SetTitleOffset(ps.xTitleOffset);
    g->GetYaxis()->SetTitleOffset(ps.yTitleOffset);
}

static void SaveNiceGraph(TGraphErrors* g,
                          const std::string& outPathNoExt,
                          const std::string& xTitle,
                          const std::string& yTitle,
                          const std::string& mainTitle,
                          const std::string& subTitle1,
                          const std::string& subTitle2,
                          const std::vector<std::string>& boxLines,
                          const PlotStyle& ps = PlotStyle(),
                          GraphTuner tuner = nullptr)
{
    if (!g) return;

    gStyle->SetOptStat(ps.showStats ? 1110 : 0);

    TCanvas* c = new TCanvas(Form("c_%s", g->GetName()), "", ps.w, ps.h);
    ApplyCanvasStyle(c, ps);
    c->cd();

    ApplyGraphStyle(g, ps, xTitle, yTitle);

    g->Draw(ps.drawOptGraph.c_str());

    DrawHeader(ps, mainTitle, subTitle1, subTitle2);
    DrawInfoBox(ps, boxLines);

    if (tuner) tuner(c, g);

    CMS_lumi(c, 13, 10);
    c->Modified();
    c->Update();

    c->SaveAs((outPathNoExt + ".png").c_str());
    c->SaveAs((outPathNoExt + ".pdf").c_str());

    delete c;
}