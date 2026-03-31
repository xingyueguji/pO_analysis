#include "plotting_helper.C" // where your PlotStyle + SaveNicePlotGraph live

void plotZcurve()
{
    TFile *f1 = TFile::Open("../skim/ZToMuMu_pO2025_hist.root", "READ");
    auto *h_Z = (TH1D *)f1->Get("hMass");
    auto *h_extended = (TH1D *)f1->Get("hMass_extended");
    auto *h_Z_vipul = (TH1D *)f1->Get("hMass_vipul");

    PlotStyle ps;
    ps.logy = true;
    ps.boxX1 = 0.15;
    ps.boxY1 = 0.51;
    ps.boxTextSize = 0.025;
    ps.boxY2 = 0.78;
    ps.boxTextSize = 0.022;
    ps.drawOpt = "E";

    PlotTuner tunemass = [&](TCanvas *c, TH1 *h)
    {
        (void)c;
        if (!h)
            return;

        // Example: ensure y-min is 0 for linear plots
        if (!ps.logy)
            h->SetMinimum(0.0);

        // Example: leave some headroom
        h->SetMaximum(2.5 * h->GetMaximum());
        h->SetMinimum(0.8);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(1.2);
        h->SetLineWidth(2);
        h->SetLineColor(kBlack);
        h->SetMarkerColor(kBlack);
    };

    const std::string outBase = "./plots";
    const std::string outDir = outBase + "/ZPeak";
    gSystem->mkdir(outDir.c_str(), kTRUE);

    h_extended->Rebin(2);
    h_Z->Rebin(2);
    h_Z_vipul->Rebin(2);

    double nZ = h_Z->Integral(1, h_Z->GetNbinsX());

    std::vector<std::string> boxLines = {
        "Z #rightarrow #mu^{+}#mu^{-}",
        Form("N_{Z} = %.0f", nZ),
        "p_{T}^{lead} > 15 GeV p_{T}^{sub} > 10 GeV",
        "|#eta| < 2.4",
        "standard iso < 0.2",
        "TightID",
        "HLT_OxyL1SingleMuOpen_v"};
    double nZ_extended = h_extended->Integral(1, h_extended->GetNbinsX());

    std::vector<std::string> boxLines_extended = {
        "Z #rightarrow #mu^{+}#mu^{-}",
        Form("N_{Z} = %.0f", nZ_extended),
        "p_{T}^{lead} > 15 GeV p_{T}^{sub} > 10 GeV",
        "|#eta| < 2.4",
        "standard iso < 0.2",
        "TightID",
        "HLT_OxyL1SingleMuOpen_v"};

    double nZ_vipul = h_Z_vipul->Integral(1, h_Z_vipul->GetNbinsX());

    std::vector<std::string> boxLines_vipul = {
        "Z #rightarrow #mu^{+}#mu^{-}",
        Form("N_{Z} = %.0f", nZ_vipul),
        "p_{T}^{lead} > 20 GeV p_{T}^{sub} > 20 GeV",
        "|#eta| < 2.4",
        "TightID",
        "HLT_OxyL1SingleMuOpen_v"};

    SaveNicePlot1D(h_Z,
                   outDir + "/ZMass",
                   "m_{#mu#mu} (GeV)", // xTitle (change if you want)
                   "Events / 1.0 GeV", // yTitle
                   "",                 // mainTitle
                   "",
                   "",
                   boxLines, // box lines (or fill later)
                   ps,
                   tunemass);

    SaveNicePlot1D(h_extended,
                   outDir + "/ZMass_extended",
                   "m_{#mu#mu} (GeV)", // xTitle (change if you want)
                   "Events / 1.0 GeV", // yTitle
                   "",                 // mainTitle
                   "",
                   "",
                   boxLines_extended, // box lines (or fill later)
                   ps,
                   tunemass);

    SaveNicePlot1D(h_Z_vipul,
                   outDir + "/ZMass_vipul",
                   "m_{#mu#mu} (GeV)", // xTitle (change if you want)
                   "Events / 1.0 GeV", // yTitle
                   "",                 // mainTitle
                   "",
                   "",
                   boxLines_vipul, // box lines (or fill later)
                   ps,
                   tunemass);
}