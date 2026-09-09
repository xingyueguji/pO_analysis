// =============================================================================
// postfit_incl.C -- rapidity-INCLUSIVE postfit stacks of the simfit GRAND
// SIMULTANEOUS FIT (2026-08-05): data vs the fitted model summed over the 12
// LAB rapidity bins, per flavour, for W+, W-, and W (both charges).
//
// The grand fit has no inclusive channel (that would double-count the per-bin
// ones), but every simfit parameter is a PURE NORMALIZATION (no shape
// nuisances, no autoMCStats) -- including the 2026-08-17 lnN nuisances (QCD
// kappa^theta per flavour x charge, global lumi on all MC), which are rate-only
// too -- so each channel's postfit shape is EXACTLY its prefit template times
// its fitted scale:
//
//     postfit(bin i) = lumi * [ r_<C>_y<i> * (signal_i + wtau_i)
//                             + r_Z        * (z_i + ztau_i) ]
//                    + qcdScale_<F>_<C> * qcd_i
//
// where lumi = kLumi^theta_lumi and qcdScale = kappa^theta are read as
// MULTIPLIERS from the CSV (extract_pO_simfit.C writes them into the qcd_norm
// columns and the appended lumi column; older CSVs without the lumi column
// fall back to lumi = 1, and free-rateParam CSVs carry the per-channel
// qcd_norm as before -- both layouts work).
// QCD_MODE=abcd (2026-08-23): still pure normalizations -- qcdScale is then
// the evaluated formula (sB*sC/sD x kappa^theta) and applies to the
// A0-normalized `qcd_abcd` template instead of `qcd`; selected via the CSV's
// 18th column qcd_model ("abcd"); column absent -> plain qcd, unchanged.
// Summing those over i reproduces combine's postfit sum identically -- and
// needs only the structured inputs (combine_input_W*.root) plus the fitted
// parameters (simfit/summary/comb_W_yields.csv), i.e. no fitDiagnostics file.
//
// WITH SHAPE NUISANCES (2026-09-07: the LHE nPDF / qcdScale / alphaS rows,
// recorded in the fitted cards' sidecar datacards/qcd_lnn_kappas.txt as
// `lheSysts a,b,c`) the shortcut is no longer exact: the fitted thetas morph
// every MC template. Then the per-channel POSTFIT shapes are read from
// shapes_fit_s/<F>_<C>_lab_y<i>/<process> of fits/simfit_lab/
// fitDiagnostics_simfit_lab.root (pulled by `sync_lxplus.sh download`),
// remapped onto the input's physical axis (FitDiagnostics stores unit-width
// bins) and summed over the 12 y bins -- same procedure as the fork's
// draw_postfit_pO.C. If that file is missing the macro falls back to the
// prefit x scale reconstruction WITH A LOUD WARNING and labels the plot
// "approx" (the theta-shape effects are ignored there).
//
//   postfit_incl(disc)  -> plots/comb/postfit_incl/<disc>/postfit_{mu,ele}_{Wp,Wm,W}.{png,pdf}
//   postfit_incl(disc, csv, muW, eleW, lheSysts, fitDiag): explicit paths;
//     lheSysts "" = read the sidecar, "none" = force the shortcut path
//
// Same cosmetics as the per-bin postfit plots (SaveNicePlot1D_WithBkg, pull
// pad). Run by analysis/run_observables.sh (comb chain), or by hand:
//   root -l -b -q 'postfit_incl.C+("leppt_mt40")'
// =============================================================================
#include "plotting_helper.C"
#include "disc_variants.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

// fitted parameters of ONE binning variant (lab), parsed from comb_W_yields.csv:
// region,charge,binning,ybin,r,rErr,signal_prefit,fitted_yield,fitted_yield_err,
// qcd_norm_mu,qcd_norm_muErr,qcd_norm_ele,qcd_norm_eleErr,r_Z,r_ZErr[,lumi,lumiErr]
struct SimfitPars
{
    double r[2][12];      // [Wp,Wm][ybin]
    double qcdMu[2][12];  // factor on the muon qcd template (rateParam or kappa^theta)
    double qcdEle[2][12]; // factor on the electron qcd template
    double rZ = 1.0, rZe = 0.0;
    double lumi = 1.0;    // kLumi^theta multiplier on all MC (1 = column absent)
    bool qcdShared = false; // lnN mode: one QCD param per (flavour, charge)
    // 18th CSV column (2026-08-23): "abcd" -> the multiplier applies to the
    // qcd_abcd template (in-fit ABCD normalization), else the plain qcd.
    // Column absent (pre-2026-08-23 CSVs) -> "" -> qcd, unchanged behavior.
    std::string qcdModel;
    bool ok = false;
};

SimfitPars ReadPars(const TString &csv)
{
    SimfitPars p;
    std::ifstream in(csv.Data());
    if (!in) { std::cerr << "[ERROR] cannot open CSV: " << csv << "\n"; return p; }
    std::string line;
    std::getline(in, line); // header
    int found = 0;
    while (std::getline(in, line))
    {
        std::stringstream ss(line);
        std::string f;
        std::vector<std::string> c;
        while (std::getline(ss, f, ',')) c.push_back(f);
        if (c.size() < 15 || c[2] != "lab") continue;
        const int ic = (c[1] == "Wp") ? 0 : (c[1] == "Wm") ? 1 : -1;
        const int iy = std::atoi(c[3].c_str());
        if (ic < 0 || iy < 0 || iy >= 12) continue;
        p.r[ic][iy]      = std::atof(c[4].c_str());
        p.qcdMu[ic][iy]  = std::atof(c[9].c_str());
        p.qcdEle[ic][iy] = std::atof(c[11].c_str());
        p.rZ  = std::atof(c[13].c_str());
        p.rZe = std::atof(c[14].c_str());
        if (c.size() >= 17) p.lumi = std::atof(c[15].c_str()); // 2026-08-17 column
        if (c.size() >= 18) p.qcdModel = c[17];                // 2026-08-23 column
        ++found;
    }
    p.ok = (found == 24);
    if (!p.ok) std::cerr << "[ERROR] expected 24 lab rows in " << csv << ", found " << found << "\n";
    // lnN mode leaves one shared QCD factor per (flavour, charge): detect it as
    // all 12 y values identical (used only to count params for the chi2 ndf)
    p.qcdShared = p.ok;
    for (int ic = 0; ic < 2 && p.qcdShared; ++ic)
        for (int iy = 1; iy < 12 && p.qcdShared; ++iy)
            if (std::fabs(p.qcdMu[ic][iy] - p.qcdMu[ic][0]) > 1e-9 ||
                std::fabs(p.qcdEle[ic][iy] - p.qcdEle[ic][0]) > 1e-9)
                p.qcdShared = false;
    return p;
}

// Poisson (Baker-Cousins) chi2, as in the fork's draw_postfit_pO.C
double BCChi2(const TH1 *d, const TH1 *t, int &nUsed)
{
    nUsed = 0;
    double chi2 = 0.0;
    const int n = std::min(d->GetNbinsX(), t->GetNbinsX());
    for (int i = 1; i <= n; ++i)
    {
        double di = d->GetBinContent(i), ti = t->GetBinContent(i);
        if (ti <= 1e-9 && di <= 0.0) continue;
        if (ti < 1e-9) ti = 1e-9;
        double term = ti - di;
        if (di > 0.0) term += di * std::log(di / ti);
        chi2 += 2.0 * term;
        ++nUsed;
    }
    return chi2;
}

// clone-a-file-histogram (never mutate file-owned objects) or add into acc
void AccTemplate(TFile *f, const TString &path, TH1D *&acc, const char *nm)
{
    TH1 *h = (TH1 *)f->Get(path);
    if (!h) { std::cerr << "[WARN] missing " << path << "\n"; return; }
    if (!acc)
    {
        acc = (TH1D *)h->Clone(nm);
        acc->SetDirectory(nullptr);
        acc->Reset();
    }
    acc->Add(h);
}

// The `lheSysts` line of the fitted cards' sidecar ("" when absent/none).
std::string ReadLheSysts(const TString &sidecar)
{
    std::ifstream in(sidecar.Data());
    if (!in) return "";
    std::string key, val;
    while (in >> key >> val)
        if (key == "lheSysts") return (val == "none") ? "" : val;
    return "";
}

// Add a FitDiagnostics postfit shape (unit-width bins 1..N, contents = events)
// into acc on the physical axis of ref (the input's data_obs), as the fork's
// draw_postfit_pO.C::RemapToRef does; bin errors = the postfit uncertainties.
bool AccPostfitShape(TFile *fd, const TString &path, const TH1 *ref, TH1D *&acc, const char *nm)
{
    TH1 *src = (TH1 *)fd->Get(path);
    if (!src) return false;
    if (!acc)
    {
        acc = (TH1D *)ref->Clone(nm);
        acc->SetDirectory(nullptr);
        acc->Reset();
    }
    const int n = std::min(src->GetNbinsX(), acc->GetNbinsX());
    if (src->GetNbinsX() != acc->GetNbinsX())
        std::cerr << "[WARN] " << path << ": postfit nbins " << src->GetNbinsX()
                  << " != input " << acc->GetNbinsX() << " (copying " << n << ")\n";
    for (int b = 1; b <= n; ++b)
    {
        acc->SetBinContent(b, acc->GetBinContent(b) + src->GetBinContent(b));
        const double e0 = acc->GetBinError(b), e1 = src->GetBinError(b);
        acc->SetBinError(b, std::sqrt(e0 * e0 + e1 * e1));
    }
    return true;
}

} // namespace

void postfit_incl(const char *disc = "met",
                  const char *csv = nullptr,       // default: pO_fit_out<suffix>/simfit/summary/comb_W_yields.csv
                  const char *muWFile = nullptr,   // default: plots/combine_input_W<suffix>.root
                  const char *eleWFile = nullptr,  // default: plots/Elec/combine_input_W<suffix>.root
                  const char *lheSysts = "",       // "" = read the fitted cards' sidecar; "none" = shortcut path
                  const char *fitDiag = nullptr)   // default: pO_fit_out<suffix>/simfit/fits/simfit_lab/fitDiagnostics_simfit_lab.root
{
    TString dsuf, discLabel;
    if (!pODisc::Spec(disc, dsuf, discLabel)) return;

    const TString fork = TString::Format("../../HiggsAnalysis-CombinedLimit/test/pO_fit_out%s/simfit", dsuf.Data());
    const TString sCsv = csv ? TString(csv) : fork + "/summary/comb_W_yields.csv";
    const TString sMuW  = muWFile  ? TString(muWFile)  : TString::Format("./plots/combine_input_W%s.root", dsuf.Data());
    const TString sEleW = eleWFile ? TString(eleWFile) : TString::Format("./plots/Elec/combine_input_W%s.root", dsuf.Data());
    const TString sFD   = fitDiag ? TString(fitDiag) : fork + "/fits/simfit_lab/fitDiagnostics_simfit_lab.root";

    SimfitPars pars = ReadPars(sCsv);
    if (!pars.ok) return;

    // shape nuisances in the fit? -> the postfit shapes must come from shapes_fit_s
    std::string lhe = (TString(lheSysts) == "") ? ReadLheSysts(fork + "/datacards/qcd_lnn_kappas.txt")
                                                : ((TString(lheSysts) == "none") ? "" : std::string(lheSysts));
    int nLhe = 0;
    for (size_t i = 0; i < lhe.size(); ++i) if (lhe[i] == ',') ++nLhe;
    if (!lhe.empty()) ++nLhe;
    TFile *fFD = nullptr;
    bool useShapes = false;
    if (nLhe > 0)
    {
        fFD = TFile::Open(sFD, "READ");
        if (fFD && !fFD->IsZombie() && fFD->Get("shapes_fit_s"))
        {
            useShapes = true;
            std::cout << "[postfit-incl] shape nuisances in the fit (" << lhe
                      << ") -> postfit shapes from " << sFD << "\n";
        }
        else
        {
            if (fFD) { fFD->Close(); delete fFD; fFD = nullptr; }
            std::cerr << "[WARN] postfit_incl: the fit has shape nuisances (" << lhe << ") but " << sFD
                      << " is missing or has no shapes_fit_s -> falling back to prefit x fitted scale,"
                         " which IGNORES the fitted shape nuisances (approximate; run"
                         " sync_lxplus.sh download to fetch the fitDiagnostics files)\n";
        }
    }

    const bool isMET = (TString(disc) == "met");
    const char *xTitle = isMET ? "PF MET (GeV)" : "Lepton p_{T} (GeV)";
    const char *yTitle = "Events / 2.0 GeV";

    const std::string outDir = std::string("./plots/comb/postfit_incl/") + disc;
    gSystem->mkdir(outDir.c_str(), kTRUE);

    const char *flavs[2] = {"mu", "ele"};
    for (int fl = 0; fl < 2; ++fl)
    {
        const TString fpath = (fl == 0) ? sMuW : sEleW;
        TFile *fW = TFile::Open(fpath, "READ");
        if (!fW || fW->IsZombie()) { std::cerr << "[ERROR] cannot open " << fpath << "\n"; continue; }
        const char *lepLab = (fl == 0) ? "W #rightarrow #mu #nu" : "W #rightarrow e #nu";

        // charge sets: Wp only, Wm only, both
        const char *tag[3] = {"Wp", "Wm", "W"};
        for (int is = 0; is < 3; ++is)
        {
            TH1D *hData = nullptr, *hSig = nullptr, *hZ = nullptr, *hZt = nullptr, *hWt = nullptr, *hQ = nullptr;
            for (int ic = 0; ic < 2; ++ic)
            {
                if (is == 0 && ic == 1) continue; // Wp only
                if (is == 1 && ic == 0) continue; // Wm only
                const char *C = (ic == 0) ? "Wp" : "Wm";
                for (int iy = 0; iy < 12; ++iy)
                {
                    const TString R = TString::Format("%s_lab_y%d", C, iy);
                    // data (unscaled)
                    AccTemplate(fW, R + "/data_obs", hData, Form("pfincl_data_%s_%s", flavs[fl], tag[is]));
                    if (useShapes)
                    {
                        // the fit's own postfit templates (thetas applied), per
                        // grand-fit channel <F>_<C>_lab_y<i>, on the input's axis
                        const TH1 *ref = (TH1 *)fW->Get(R + "/data_obs");
                        const TString ch = TString::Format("shapes_fit_s/%s_%s", flavs[fl], R.Data());
                        bool okc = ref != nullptr;
                        okc &= AccPostfitShape(fFD, ch + "/signal", ref, hSig, Form("pfincl_sig_%s_%s", flavs[fl], tag[is]));
                        okc &= AccPostfitShape(fFD, ch + "/wtau",   ref, hWt,  Form("pfincl_wt_%s_%s",  flavs[fl], tag[is]));
                        okc &= AccPostfitShape(fFD, ch + "/z",      ref, hZ,   Form("pfincl_z_%s_%s",   flavs[fl], tag[is]));
                        okc &= AccPostfitShape(fFD, ch + "/ztau",   ref, hZt,  Form("pfincl_zt_%s_%s",  flavs[fl], tag[is]));
                        okc &= AccPostfitShape(fFD, ch + "/qcd",    ref, hQ,   Form("pfincl_q_%s_%s",   flavs[fl], tag[is]));
                        if (!okc) std::cerr << "[WARN] postfit_incl: incomplete " << ch << " in " << sFD << "\n";
                        continue;
                    }
                    // postfit-scaled templates: clone region histos, scale, add
                    auto addScaled = [&](const char *proc, TH1D *&acc, double scale, const char *nm) {
                        TH1 *h = (TH1 *)fW->Get(R + "/" + proc);
                        if (!h) { std::cerr << "[WARN] missing " << R << "/" << proc << "\n"; return; }
                        TH1D *tmp = (TH1D *)h->Clone(Form("tmp_%s_%s", nm, R.Data()));
                        tmp->SetDirectory(nullptr);
                        tmp->Scale(scale);
                        if (!acc) { acc = (TH1D *)tmp->Clone(nm); acc->SetDirectory(nullptr); }
                        else acc->Add(tmp);
                        delete tmp;
                    };
                    // lumi nuisance multiplies every MC template, not the
                    // data-driven qcd (whose factor already IS its full scale)
                    const double rB   = pars.r[ic][iy] * pars.lumi;
                    const double rZB  = pars.rZ * pars.lumi;
                    const double qcdB = (fl == 0) ? pars.qcdMu[ic][iy] : pars.qcdEle[ic][iy];
                    // abcd mode: the multiplier (sB*sC/sD x kappa^theta) applies
                    // to the A0-normalized qcd_abcd template, not the plain qcd
                    const char *qcdProc = (pars.qcdModel == "abcd") ? "qcd_abcd" : "qcd";
                    addScaled("signal", hSig, rB,   Form("pfincl_sig_%s_%s", flavs[fl], tag[is]));
                    addScaled("wtau",   hWt,  rB,   Form("pfincl_wt_%s_%s",  flavs[fl], tag[is]));
                    addScaled("z",      hZ,   rZB,  Form("pfincl_z_%s_%s",   flavs[fl], tag[is]));
                    addScaled("ztau",   hZt,  rZB,  Form("pfincl_zt_%s_%s",  flavs[fl], tag[is]));
                    addScaled(qcdProc,  hQ,   qcdB, Form("pfincl_q_%s_%s",   flavs[fl], tag[is]));
                }
            }
            if (!hData || !hSig) { std::cerr << "[ERROR] no histograms accumulated for " << tag[is] << "\n"; continue; }

            // stack order matches the per-bin postfit plots (draw_postfit_pO.C)
            std::vector<TH1 *> bkgs = {hSig, hZ, hZt, hWt, hQ};
            std::vector<std::string> names = {"W signal", "DY", "DY #tau", "W #tau", "QCD (ABCD)"};

            TH1D *hTot = (TH1D *)hSig->Clone(Form("pfincl_tot_%s_%s", flavs[fl], tag[is]));
            hTot->SetDirectory(nullptr);
            if (hZ) hTot->Add(hZ);
            if (hZt) hTot->Add(hZt);
            if (hWt) hTot->Add(hWt);
            if (hQ) hTot->Add(hQ);

            // info box: same content style as the per-bin plots; the per-bin r's
            // cannot be quoted here, so point at the y-binned plots instead.
            int nUsed = 0;
            const double chi2 = BCChi2(hData, hTot, nUsed);
            // approx params shaping this sum: 12(24) r's + r_Z + the QCD params
            // (lnN mode: 1 shared per flavour x charge; free mode: 12(24)) + lumi.
            // abcd mode: the multiplier is y-identical so qcdShared fires ->
            // 1 QCD param per (flavour, charge) -- the right count here, since
            // the 3 CR scales are constrained by CR data outside this sum.
            // + the LHE shape nuisances (one theta each, shared by everything)
            const int ndfPars = (pars.qcdShared ? ((is == 2) ? 28 : 15)
                                                : ((is == 2) ? 49 : 25)) + nLhe;
            int ndf = nUsed - ndfPars;
            if (ndf < 1) ndf = (nUsed > 0 ? nUsed : 1);
            std::vector<std::string> box = {
                Form("Data: %.0f", hData->Integral()),
                Form("Postfit total: %.0f", hTot->Integral()),
                Form("#chi^{2}/ndf = %.2f, p = %.2f", chi2 / ndf, TMath::Prob(chi2, ndf)),
                Form("DY norm = %.3f #pm %.3f", pars.rZ, pars.rZe),
                useShapes ? "sum of 12 lab y bins (shapes_fit_s)"
                          : (nLhe > 0 ? "sum of 12 lab y bins (APPROX: shape #theta's ignored)"
                                      : "sum of 12 lab y bins (per-bin r's)")};

            PlotStyle ps;
            ps.drawOpt = "hist";
            ps.showStats = false;
            ps.logy = true;
            ps.normBkgToData = false; // ABSOLUTE postfit yields
            ps.pullPad = true;
            ps.headerX = 0.56;
            ps.boxTextSize = 0.028;
            ps.boxX1 = 0.56; ps.boxX2 = 0.93;
            ps.boxY1 = 0.46; ps.boxY2 = 0.76;
            ps.legX1 = 0.72; ps.legY1 = 0.275;
            ps.legX2 = 0.93; ps.legY2 = 0.455;
            if (!isMET) // pT axes start at the selection floor (edge 24
            {           // encloses the 25 GeV cut on the 2 GeV grid)
                ps.xRangeLo = 24.0;
                ps.xRangeHi = 100.0;
            }
            PlotTuner tuner = [&](TCanvas *c, TH1 *h) {
                (void)c;
                if (!h) return;
                h->SetMinimum(1.0);
                h->SetMaximum(10.25 * h->GetMaximum());
            };

            const std::string out = outDir + "/postfit_" + flavs[fl] + "_" + tag[is];
            SaveNicePlot1D_WithBkg(hData, bkgs, names, out, xTitle, yTitle,
                                   "", lepLab,
                                   Form("%s lab incl (simfit postfit)", tag[is]),
                                   box, ps, tuner);
            std::cout << "[postfit-incl] " << out << ".png\n";
        }
        fW->Close();
        delete fW;
    }
    if (fFD) { fFD->Close(); delete fFD; }
    std::cout << "[OK] inclusive simfit postfit stacks (disc=" << disc << ", "
              << (useShapes ? "shapes_fit_s" : (nLhe > 0 ? "APPROX prefit x scale" : "prefit x scale"))
              << ") -> " << outDir << "\n";
}
