#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TSystem.h"
#include "plotting_helper.C"
#include "../skim/mc_norm.h"   // pONorm::MCScale -> per-sample k_s = sigma*L/N_gen
#include "../skim/lhe_index.h" // pOLhe::kLheSystNames: the <hist>_<syst>Up/Down twins written by skim/lhe_updown.py
#include "../skim/muon_sf.h"   // pOSF::kMuonSFSystNames: the <hist>_mu{ID,Iso,Trig}Up/Down twins written by skim.C (muon MC only)
#include "../skim/electron_sf.h" // pOSF::kEleSFSystNames: the <hist>_eleSFUp/Down twins written by skim.C (electron MC only, 2026-09-29)
#include "../skim/met_res.h"     // pOMetRes::SystName: the <hist>_metRes{Mu,Ele}Up/Down twins written by skim.C (MC only, 2026-10-20)
#include "../correction/abcd_ygroups.h" // pOAbcdY: the rapidity groups of the abcdy in-fit ABCD (QCD_MODE=abcdy, 2026-10-19)

#include <algorithm>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <functional>
#include <initializer_list>

void mtandmet(bool isElec = 1)
{
    std::string inFile;
    std::string inFile_MC_signal_Wp;
    std::string inFile_MC_signal_Wm;
    std::string inFile_MC_background_Z;
    std::string inFile_MC_background_Ztau;
    std::string inFile_MC_background_Wptau;
    std::string inFile_MC_background_Wmtau;
    std::string inFile_MC_background_TT; // ttbar (2026-10-16): one skim file per channel, like the taus

    const char *Channeltype;
    const char *Channeltypewplus;
    const char *Channeltypewminus;

    if (isElec)
    {
        inFile = "../skim/rootfile/WToElecNu_pO_PFMet_hist.root";
        inFile_MC_signal_Wp = "../skim/rootfile/WToElecNu_pO_PFMet_Wp_hist.root";
        inFile_MC_signal_Wm = "../skim/rootfile/WToElecNu_pO_PFMet_Wm_hist.root";
        inFile_MC_background_Z = "../skim/rootfile/WToElecNu_pO_PFMet_DY_hist.root";
        inFile_MC_background_Ztau  = "../skim/rootfile/WToElecNu_pO_PFMet_DYtau_hist.root";
        inFile_MC_background_Wptau = "../skim/rootfile/WToElecNu_pO_PFMet_Wptau_hist.root";
        inFile_MC_background_Wmtau = "../skim/rootfile/WToElecNu_pO_PFMet_Wmtau_hist.root";
        inFile_MC_background_TT    = "../skim/rootfile/WToElecNu_pO_PFMet_TT_hist.root";
        Channeltype = "W #rightarrow e #nu";
        Channeltypewplus = "W^{+} #rightarrow e^{+} #nu";
        Channeltypewminus = "W^{-} #rightarrow e^{-}#bar{#nu}";
    }
    else
    {
        inFile = "../skim/rootfile/WToMuNu_pO_PFMet_hist.root";
        inFile_MC_signal_Wp = "../skim/rootfile/WToMuNu_pO_PFMet_Wp_hist.root";
        inFile_MC_signal_Wm = "../skim/rootfile/WToMuNu_pO_PFMet_Wm_hist.root";
        inFile_MC_background_Z = "../skim/rootfile/WToMuNu_pO_PFMet_DY_hist.root";
        inFile_MC_background_Ztau = "../skim/rootfile/WToMuNu_pO_PFMet_DYtau_hist.root";
        inFile_MC_background_Wptau = "../skim/rootfile/WToMuNu_pO_PFMet_Wptau_hist.root";
        inFile_MC_background_Wmtau = "../skim/rootfile/WToMuNu_pO_PFMet_Wmtau_hist.root";
        inFile_MC_background_TT    = "../skim/rootfile/WToMuNu_pO_PFMet_TT_hist.root";
        Channeltype = "W #rightarrow #mu #nu";
        Channeltypewplus = "W^{+} #rightarrow #mu^{+} #nu";
        Channeltypewminus = "W^{-} #rightarrow #mu^{-}#bar{#nu}";
    }

    TFile *f = TFile::Open(inFile.c_str(), "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open file: " << inFile << "\n";
        return;
    }

    TFile *f_Wp = TFile::Open(inFile_MC_signal_Wp.c_str(), "READ");
    if (!f_Wp || f_Wp->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open file: " << inFile_MC_signal_Wp << "\n";
        return;
    }

    TFile *f_Wm = TFile::Open(inFile_MC_signal_Wm.c_str(), "READ");
    if (!f_Wm || f_Wm->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open file: " << inFile_MC_signal_Wm << "\n";
        return;
    }

    TFile *f_DY = TFile::Open(inFile_MC_background_Z.c_str(), "READ");
    if (!f_DY || f_DY->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open file: " << inFile_MC_background_Z << "\n";
        return;
    }

    TFile *f_DYtau = TFile::Open(inFile_MC_background_Ztau.c_str(), "READ");
    if (!f_DYtau || f_DYtau->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open file: " << inFile_MC_background_Ztau << "\n";
        return;
    }

    TFile *f_Wptau = TFile::Open(inFile_MC_background_Wptau.c_str(), "READ");
    if (!f_Wptau || f_Wptau->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open file: " << inFile_MC_background_Wptau << "\n";
        return;
    }

    TFile *f_Wmtau = TFile::Open(inFile_MC_background_Wmtau.c_str(), "READ");
    if (!f_Wmtau || f_Wmtau->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open file: " << inFile_MC_background_Wmtau << "\n";
        return;
    }

    TFile *f_TT = TFile::Open(inFile_MC_background_TT.c_str(), "READ");
    if (!f_TT || f_TT->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open file: " << inFile_MC_background_TT
                  << " (skim the TT sample: skim/run_all.sh <Wmu|Wel> TT)\n";
        return;
    }

    // output directories
    std::string outBase;
    if (isElec)
    {
        outBase = "./plots/Elec";
    }
    else
    {
        outBase = "./plots";
    }
    // Three display families: met (the backup discriminant), mt (the m_T shape
    // check), leppt_mt40 (the primary discriminant). The fourth, leppt (plain
    // W selection), was retired 2026-09-21 -- see the variant table below.
    const std::string outMetDir = outBase + "/met";
    const std::string outMtDir = outBase + "/mt";
    const std::string outLepPtMt40Dir = outBase + "/leppt_mt40";
    const std::string outMetMt40Dir = outBase + "/met_mt40"; // the MET counterpart (2026-10-20)
    gSystem->mkdir(outMetDir.c_str(), kTRUE);
    gSystem->mkdir(outMtDir.c_str(), kTRUE);
    gSystem->mkdir(outLepPtMt40Dir.c_str(), kTRUE);
    gSystem->mkdir(outMetMt40Dir.c_str(), kTRUE);

    // Lepton-pT stacks (2026-07-29) are DISPLAY ONLY: they are never written
    // into combine_input_W.root -- the fit discriminant stays PF MET.
    const char *ptTitle = isElec ? "p_{T}^{e} (GeV)" : "p_{T}^{#mu} (GeV)";

    const int NY = 12;
    double yEdges[NY + 1] = {
        -2.4, -2.0, -1.6, -1.2, -0.8, -0.4,
        0.0, 0.4, 0.8, 1.2, 1.6, 2.0,
        2.4};
    double yEdges_FB[NY + 1] = {
        -1.7068,
        -1.3646,
        -1.0223,
        -0.6801,
        -0.3379,
        0.0044,
        0.3466,
        0.6888,
        1.0311,
        1.3733,
        1.7155,
        2.0578,
        2.4000};

    // If you already have your yEdges, you can format them here.
    // Otherwise just label by index (y0..y7) for now.
    // Example placeholders:
    std::vector<std::string> yLabel(NY);
    std::vector<std::string> yLabel_FB(NY);
    const char *etaSym = isElec ? "#eta_{e}" : "#eta_{#mu}"; // lepton-correct axis label
    for (int iy = 0; iy < NY; ++iy)
    {
        yLabel[iy] = RangeLabel(etaSym, iy, yEdges, NY);       // or "#eta_{CM}"
        yLabel_FB[iy] = RangeLabel(etaSym, iy, yEdges_FB, NY); // or "#eta_{CM}"
    }

    // Common style
    PlotStyle ps;
    ps.drawOpt = "hist";
    ps.showStats = false;
    ps.logy = true; // set true if you want log-y for MET tails
    ps.boxY1 = 0.62;
    ps.boxY2 = 0.82;
    ps.normBkgToData = false; // ABSOLUTE pO scaling (k_s incl. A=16) -- no area norm
    ps.pullPad = true;        // (data-MC)/sigma bars in a sub-pad under each plot

    // The MET and lepton-pT stacks carry THREE count lines (data, W signal MC,
    // total stack -- see totalLine) against the m_T stacks' one. DrawInfoBox
    // centres the box in [boxY1, boxY2], so a third line would push the first
    // one up into the rapidity header: shift the band (and the auto legend
    // under it) down by half a line (DrawInfoBox's line pitch is 0.05), which
    // keeps the first line exactly where it sat on the two-line box.
    PlotStyle psCounts = ps;
    psCounts.boxY1 -= 0.025;
    psCounts.boxY2 -= 0.025;
    // The lepton-pT stacks add two more lines (MC/data and r_pre, see
    // prefitLines): shift the band down by one more line so the first line
    // still sits where it does on the two-line box (top at 0.77 either way),
    // and put the legend under the taller box in TWO columns: 0.1 high (the
    // same row height as the one-column 0.2), so its bottom edge stays where
    // the auto legend's was and it covers no more of the stack than before.
    PlotStyle psLepPt = psCounts;
    psLepPt.boxY1 -= 0.05;
    psLepPt.boxY2 -= 0.05;
    {
        const double boxBottom = 0.5 * (psLepPt.boxY1 + psLepPt.boxY2) - 0.5 * 0.05 * 5;
        psLepPt.legX1 = psLepPt.boxX1;
        psLepPt.legX2 = psLepPt.boxX2;
        psLepPt.legY2 = boxBottom - 0.025;
        psLepPt.legY1 = psLepPt.legY2 - 0.1;
        psLepPt.legNColumns = 2;
    }

    // (Optional) one tuner you can reuse for all plots
    PlotTuner commonTuner = [&](TCanvas *c, TH1 *h)
    {
        (void)c;
        if (!h)
            return;

        // Example: ensure y-min is 0 for linear plots
        if (ps.logy)
            h->SetMinimum(1.0);

        // Example: leave some headroom
        h->SetMaximum(10.25 * h->GetMaximum());
    };

    TH1D *h_met_inclusive = nullptr;
    TH1D *h_mt_inclusive = nullptr;

    TH1D *h_met_inclusive_MC_signal = nullptr;
    TH1D *h_mt_inclusive_MC_signal = nullptr;

    TH1D *h_met_inclusive_MC_Z = nullptr;
    TH1D *h_mt_inclusive_MC_Z = nullptr;

    TH1D *h_met_inclusive_MC_Ztau = nullptr;
    TH1D *h_mt_inclusive_MC_Ztau = nullptr;

    TH1D *h_met_inclusive_MC_Wtau = nullptr;
    TH1D *h_mt_inclusive_MC_Wtau = nullptr;

    TH1D *h_met_inclusive_MC_TT = nullptr;
    TH1D *h_mt_inclusive_MC_TT = nullptr;

    // --- Absolute per-sample MC normalization k_s = A*sigma*L/N_gen -----------
    // Scale each MC sample by its own k_s so the samples sit in the correct
    // PHYSICAL proportion AND at the absolute pO yield. sigma is already inside
    // the gen weights (<w>=sigma, per-NN); MCScale multiplies by the Oxygen
    // A-scaling (A=16) and L to give the predicted pO event count. The tau
    // samples share the W/DY cross sections but carry their own N_gen (per-file
    // label). The stack is drawn ABSOLUTE (ps.normBkgToData=false above) -- no
    // area normalization -- so any data/MC gap (e.g. the low-MET QCD excess,
    // missing reco corrections) is physical, not hidden.
    const double k_Wp    = pONorm::MCScale(isElec ? "Wp_ele" : "Wp_mu");
    const double k_Wm    = pONorm::MCScale(isElec ? "Wm_ele" : "Wm_mu");
    const double k_DY    = pONorm::MCScale(isElec ? "DYee"   : "DYmu");
    const double k_DYtau = pONorm::MCScale("DYtau");
    const double k_Wptau = pONorm::MCScale("Wp_tau");
    const double k_Wmtau = pONorm::MCScale("Wm_tau");
    const double k_TT    = pONorm::MCScale("TT");    // ttbar background (2026-10-16)

    auto scaleAll = [](double k, std::initializer_list<TH1D *> hs)
    {
        for (TH1D *h : hs)
            if (h) h->Scale(k);
    };

    // ttbar (2026-10-16) goes LAST in every stack, after the QCD: the stack
    // colours are assigned by position (plotting_helper.C), so every entry that
    // existed before keeps its colour and ttbar takes the 8th one.
    auto addTT = [](std::vector<TH1 *> &bkgs, std::vector<std::string> &names, TH1 *h)
    {
        if (h) { bkgs.push_back(h); names.push_back("t#bar{t}"); }
    };

    // Sum the absolute integrals of the W signal-MC histos (skips nulls). Used to
    // annotate the MET plots with the real-W yield (W+ + W-) under "Passing Events".
    auto sigInt = [](std::initializer_list<TH1D *> hs) -> double
    {
        double s = 0.0;
        for (TH1D *h : hs)
            if (h) s += h->Integral(1, h->GetNbinsX());
        return s;
    };

    // Third info-box line of the MET / lepton-pT stacks: the TOTAL of the stack
    // as drawn -- every MC sample plus the ABCD QCD when the stack has one, i.e.
    // the same total the pull pad compares the data to. Next to "Passing
    // Events" it separates a normalization offset (data != total) from a shape
    // mismatch (totals agree, pulls still structured). Takes the SAME vectors
    // handed to SaveNicePlot1D_WithBkg, so the label cannot disagree with the
    // legend; integral over bins 1..N like the other two lines (no
    // under/overflow), so the three numbers compare directly.
    auto totalLine = [](const std::vector<TH1 *> &stack,
                        const std::vector<std::string> &names) -> std::string
    {
        double s = 0.0;
        for (TH1 *h : stack)
            if (h) s += h->Integral(1, h->GetNbinsX());
        const bool withQcd = std::any_of(names.begin(), names.end(),
            [](const std::string &n) { return n.compare(0, 3, "QCD") == 0; });
        return Form("%s: %.0f", withQcd ? "Total MC+QCD" : "Total MC", s);
    };

    // Two more info-box lines of the lepton-pT stacks (2026-10-13): where the
    // prefit MC sits against the data, from the stack AS DRAWN (bins 1..N, like
    // the count lines above):
    //   MC/data = (W + DY + QCD)/data      < 1: the stack is BELOW the data
    //   r_pre   = (data - DY - QCD)/W      W = W signal + W->tau nu, i.e. what
    //                                      the fit's r scales; > 1: the data
    //                                      want more W (the fit's r starts here)
    // The two are not inverses of each other: r leaves the backgrounds alone.
    // Error = data Poisson only, sqrt(N)/W. NB the QCD is the DISPLAY template
    // the stack draws (qcd_pt_mt40), not the fit's in-fit-ABCD prefit
    // qcd_abcd -- plotting/prefit_r.C prints both (inclusive r differs by
    // <= 0.03).
    auto prefitLines = [](TH1 *hData, const std::vector<TH1 *> &wParts,
                          const std::vector<TH1 *> &stack) -> std::vector<std::string>
    {
        auto integ = [](TH1 *h) { return h ? h->Integral(1, h->GetNbinsX()) : 0.0; };
        const double D = integ(hData);
        double W = 0.0, tot = 0.0;
        for (TH1 *h : wParts) W += integ(h);
        for (TH1 *h : stack) tot += integ(h);
        if (D <= 0.0 || W <= 0.0) return {};
        return {Form("MC/data: %.3f", tot / D),
                Form("r_{pre} = (D#minusB)/W: %.2f#pm%.2f", (D - (tot - W)) / W, std::sqrt(D) / W)};
    };

    // --- ABCD QCD (low-MET) background, MET stacks ----------------------------
    // Inclusive-in-rapidity per-charge QCD MET template from
    // correction/qcd_abcd.C (ABCD: data-EWK anti-iso shape x transfer factor),
    // muon OR electron (qcd_abcd_{mu,ele}.root). Added to the MET stacks below.
    // The ABCD normalization is only measured INCLUSIVELY (low pO statistics), so
    // for the per-rapidity-bin plots the one inclusive template is split across y
    // by a data-driven proxy: the per-bin low-MET (MET < qcdMetCut) excess
    // (data - EWK), i.e. where the QCD sits in y, keeping the inclusive template
    // SHAPE. The lab per-y QCD therefore SUMS back to the inclusive template by
    // construction; the FB per-y QCD sums to the FB window's share of it (see
    // qcdWeights). This is an approximation for the per-bin plots; the inclusive
    // MET plot is the rigorous one.
    const double qcdMetCut = 30.0; // matches metCut in correction/qcd_abcd.C
    const std::string qlep = isElec ? "ele" : "mu";
    const std::string qcdFile = "../correction/rootfile/qcd_abcd_" + qlep + ".root";

    // Selection variants of the m_T > 40 fits (shared by the stacks, the QCD
    // templates and the Combine input writer below):
    //   kVarMt40    = pT>25 && m_T>40, LEPTON pT discriminant (h_leppt_mt40_*,
    //                 qcd_pt_mt40_*) -> combine_input_W_leppt_mt40.root, PRIMARY
    //   kVarMetMt40 = the SAME events, PF MET discriminant (h_met_mt40_*,
    //                 qcd_met_mt40_*; 2026-10-20, user decision: the apples-to-
    //                 apples MET counterpart of the lepton-pT fit) ->
    //                 combine_input_W_met_mt40.root. Same in-fit ABCD counts,
    //                 CRs and rapidity groups; only the SR shapes change, plus
    //                 the MET-resolution nuisance (metResMu/metResEle).
    // The plain no-m_T-cut lepton-pT variant (kVarNom, "h_leppt_*" ->
    // plots/leppt/ + combine_input_W_leppt.root) was RETIRED 2026-09-21: the
    // discriminant itself was dropped on 2026-08-16, its out-tree
    // pO_fit_out_leppt/ never had a simfit, and the skim stores no LHE/SF
    // systematic twins for h_leppt_* -- so its Combine input carried no shape
    // systematics and could not be fitted by the current card generator.
    enum { kVarMt40 = 0, kVarMetMt40 = 1, kNVar = 2 };
    const char *varStem[kNVar] = {"h_leppt_mt40", "h_met_mt40"};
    const char *varTag [kNVar] = {"_mt40", "_metmt40"};
    const char *varDisc[kNVar] = {"leppt_mt40", "met_mt40"}; // file / plot-folder stem
    const char *varXTitle[kNVar] = {ptTitle, "PF MET (GeV)"};
    // display x range (pT starts at the 24 GeV bin edge enclosing the 25 GeV
    // cut; MET keeps its full 0-120 axis): active when hi > lo
    const double varXLo[kNVar] = {24.0, 0.0};
    const double varXHi[kNVar] = {100.0, 0.0};
    TH1D *qcdPlusBase = nullptr, *qcdMinusBase = nullptr, *qcd_met_incl = nullptr;
    // No-m_T-cut lepton-pT QCD: since the plain-leppt variant was retired
    // (2026-09-21) these feed only the m_T-cut-efficiency line printed below.
    TH1D *qcdPtPlusBase = nullptr, *qcdPtMinusBase = nullptr;
    TH1D *qcdPtMt40PlusBase = nullptr, *qcdPtMt40MinusBase = nullptr, *qcd_pt_mt40_incl = nullptr;
    // In-fit ABCD (QCD_MODE=abcd, 2026-08-23): the exported m_T-plane counts and
    // the A0-renormalized SR template bases built from them (same anti-iso
    // m_T>40 pT SHAPE, total = A0 = qcdB0*qcdC0/qcdD0 instead of T_MET x V).
    TH1D *abcdCountsPlus = nullptr, *abcdCountsMinus = nullptr;
    TH1D *qcdAbcdPlusBase = nullptr, *qcdAbcdMinusBase = nullptr;
    // ELECTRONS since 2026-09-29 (user decision): the in-fit SR template is the
    // FAKE-FACTOR-WEIGHTED one of correction/qcd_abcd.C (buildFFTemplate:
    // T0 x f(pT) x C40(pT), f = the measured exponential pT dependence of the
    // fake factor), plus its slope +-1 sigma twins = the qcdFF_ele_<C> shape
    // nuisance (one per charge). Read as-is -- the total is the FF prediction,
    // NOT renormalized to the flat A0.
    TH1D *ffPlus = nullptr, *ffMinus = nullptr;              // qcd_pt_mt40_ff_ele{Plus,Minus}
    TH1D *ffPlusUp = nullptr, *ffPlusDn = nullptr, *ffMinusUp = nullptr, *ffMinusDn = nullptr;
    // the slope twins actually used (set only when the FF base was accepted)
    TH1D *qcdAbcdPlusUpBase = nullptr, *qcdAbcdPlusDnBase = nullptr;
    TH1D *qcdAbcdMinusUpBase = nullptr, *qcdAbcdMinusDnBase = nullptr;
    // The met_mt40 twins of everything above (2026-10-20): the display template
    // qcd_met_mt40_* (anti-iso m_T>40 MET shape x T_MET), its in-fit base
    // (muon: renormalized to A0; electron: the FF-weighted qcd_met_mt40_ff_*,
    // total checked against the exported A0ffMet) and the slope twins.
    TH1D *qcdMetMt40PlusBase = nullptr, *qcdMetMt40MinusBase = nullptr, *qcd_met_mt40_incl = nullptr;
    TH1D *qcdAbcdMetPlusBase = nullptr, *qcdAbcdMetMinusBase = nullptr;
    TH1D *ffMetPlus = nullptr, *ffMetMinus = nullptr;
    TH1D *ffMetPlusUp = nullptr, *ffMetPlusDn = nullptr, *ffMetMinusUp = nullptr, *ffMetMinusDn = nullptr;
    TH1D *qcdAbcdMetPlusUpBase = nullptr, *qcdAbcdMetPlusDnBase = nullptr;
    TH1D *qcdAbcdMetMinusUpBase = nullptr, *qcdAbcdMetMinusDnBase = nullptr;
    const std::string ffSystName[2] = {"qcdFF_" + std::string(isElec ? "ele" : "mu") + "_Wp",
                                       "qcdFF_" + std::string(isElec ? "ele" : "mu") + "_Wm"};
    // RAPIDITY-RESOLVED in-fit ABCD (QCD_MODE=abcdy, 2026-10-19, OPT-IN in the
    // fork): the per-group counts abcd_counts_<lep><chg>_<lab|fb>_<g> of
    // correction/qcd_abcd.C (runGroupABCD), [charge][binning][group], groups
    // from correction/abcd_ygroups.h. Read here, used for the qcd_abcdy SR
    // template and the group CR dirs below; nothing else depends on them.
    TH1D *abcdyCounts[2][2][pOAbcdY::kNGroups] = {};
    const char *abcdyBin[2] = {"lab", "fb"};
    {
        TFile *fqcd = TFile::Open(qcdFile.c_str(), "READ");
        if (fqcd && !fqcd->IsZombie())
        {
            TH1D *qp = (TH1D *)fqcd->Get(Form("qcd_met_%sPlus", qlep.c_str()));
            TH1D *qm = (TH1D *)fqcd->Get(Form("qcd_met_%sMinus", qlep.c_str()));
            if (qp) { qcdPlusBase  = (TH1D *)qp->Clone("qcd_met_Plus_base");  qcdPlusBase->SetDirectory(nullptr); }
            if (qm) { qcdMinusBase = (TH1D *)qm->Clone("qcd_met_Minus_base"); qcdMinusBase->SetDirectory(nullptr); }
            // lepton-pT QCD template (same iso-pass yield, anti-iso pT shape)
            TH1D *pp = (TH1D *)fqcd->Get(Form("qcd_pt_%sPlus", qlep.c_str()));
            TH1D *pm = (TH1D *)fqcd->Get(Form("qcd_pt_%sMinus", qlep.c_str()));
            // ... and the m_T>40 twin (the pT-discriminant selection's QCD)
            TH1D *pp40 = (TH1D *)fqcd->Get(Form("qcd_pt_mt40_%sPlus", qlep.c_str()));
            TH1D *pm40 = (TH1D *)fqcd->Get(Form("qcd_pt_mt40_%sMinus", qlep.c_str()));
            if (pp40) { qcdPtMt40PlusBase  = (TH1D *)pp40->Clone("qcd_pt_mt40_Plus_base");  qcdPtMt40PlusBase->SetDirectory(nullptr); }
            if (pm40) { qcdPtMt40MinusBase = (TH1D *)pm40->Clone("qcd_pt_mt40_Minus_base"); qcdPtMt40MinusBase->SetDirectory(nullptr); }
            // ... and its PF-MET twin (the met_mt40 fit's QCD, 2026-10-20)
            TH1D *mp40 = (TH1D *)fqcd->Get(Form("qcd_met_mt40_%sPlus", qlep.c_str()));
            TH1D *mm40 = (TH1D *)fqcd->Get(Form("qcd_met_mt40_%sMinus", qlep.c_str()));
            if (mp40) { qcdMetMt40PlusBase  = (TH1D *)mp40->Clone("qcd_met_mt40_Plus_base");  qcdMetMt40PlusBase->SetDirectory(nullptr); }
            if (mm40) { qcdMetMt40MinusBase = (TH1D *)mm40->Clone("qcd_met_mt40_Minus_base"); qcdMetMt40MinusBase->SetDirectory(nullptr); }
            if (pp) { qcdPtPlusBase  = (TH1D *)pp->Clone("qcd_pt_Plus_base");  qcdPtPlusBase->SetDirectory(nullptr); }
            if (pm) { qcdPtMinusBase = (TH1D *)pm->Clone("qcd_pt_Minus_base"); qcdPtMinusBase->SetDirectory(nullptr); }
            // in-fit ABCD counts (m_T plane; written by qcd_abcd.C since 2026-08-23)
            TH1D *cp = (TH1D *)fqcd->Get(Form("abcd_counts_%sPlus", qlep.c_str()));
            TH1D *cm = (TH1D *)fqcd->Get(Form("abcd_counts_%sMinus", qlep.c_str()));
            if (cp) { abcdCountsPlus  = (TH1D *)cp->Clone("abcd_counts_Plus_loc");  abcdCountsPlus->SetDirectory(nullptr); }
            if (cm) { abcdCountsMinus = (TH1D *)cm->Clone("abcd_counts_Minus_loc"); abcdCountsMinus->SetDirectory(nullptr); }
            // the per-group counts (abcdy); absent in a pre-2026-10-19 qcd_abcd run
            for (int ic = 0; ic < 2; ++ic)
                for (int ib = 0; ib < 2; ++ib)
                    for (int g = 0; g < pOAbcdY::kNGroups; ++g)
                    {
                        const std::string nm = Form("abcd_counts_%s%s_%s_%s", qlep.c_str(), ic == 0 ? "Plus" : "Minus",
                                                    abcdyBin[ib], pOAbcdY::kGroupName[g]);
                        TH1D *h = (TH1D *)fqcd->Get(nm.c_str());
                        if (!h) continue;
                        abcdyCounts[ic][ib][g] = (TH1D *)h->Clone((nm + "_loc").c_str());
                        abcdyCounts[ic][ib][g]->SetDirectory(nullptr);
                    }
            if (isElec) // the fake-factor-weighted SR template + its slope twins (2026-09-29)
            {
                auto getFF = [&](const char *nm) -> TH1D * {
                    TH1D *h = (TH1D *)fqcd->Get(nm);
                    if (!h) return nullptr;
                    h = (TH1D *)h->Clone(Form("%s_loc", nm));
                    h->SetDirectory(nullptr);
                    return h;
                };
                ffPlus    = getFF("qcd_pt_mt40_ff_elePlus");
                ffMinus   = getFF("qcd_pt_mt40_ff_eleMinus");
                ffPlusUp  = getFF("qcd_pt_mt40_ff_elePlus_ffUp");
                ffPlusDn  = getFF("qcd_pt_mt40_ff_elePlus_ffDown");
                ffMinusUp = getFF("qcd_pt_mt40_ff_eleMinus_ffUp");
                ffMinusDn = getFF("qcd_pt_mt40_ff_eleMinus_ffDown");
                // the same fake-factor model on PF MET (met_mt40, 2026-10-20)
                ffMetPlus    = getFF("qcd_met_mt40_ff_elePlus");
                ffMetMinus   = getFF("qcd_met_mt40_ff_eleMinus");
                ffMetPlusUp  = getFF("qcd_met_mt40_ff_elePlus_ffUp");
                ffMetPlusDn  = getFF("qcd_met_mt40_ff_elePlus_ffDown");
                ffMetMinusUp = getFF("qcd_met_mt40_ff_eleMinus_ffUp");
                ffMetMinusDn = getFF("qcd_met_mt40_ff_eleMinus_ffDown");
            }
            fqcd->Close();
            delete fqcd;
        }
        if (qcdPlusBase && qcdMinusBase)
        {
            qcd_met_incl = (TH1D *)qcdPlusBase->Clone("qcd_met_incl"); // both charges
            qcd_met_incl->SetDirectory(nullptr);
            qcd_met_incl->Add(qcdMinusBase);
        }
        else
            std::cerr << "[WARN] ABCD QCD template not found (" << qcdFile
                      << "); run qcd_abcd.C" << (isElec ? "+(true)" : "+")
                      << " first. MET stacks will omit QCD.\n";
        if (!qcdPtPlusBase || !qcdPtMinusBase)
            std::cerr << "[WARN] no-m_T-cut lepton-pT QCD template (qcd_pt_*) not in " << qcdFile
                      << "; re-run qcd_abcd.C on a skim with the h_iso_pt planes."
                      << " Only the m_T-cut-efficiency line below is affected.\n";
        if (qcdPtMt40PlusBase && qcdPtMt40MinusBase)
        {
            qcd_pt_mt40_incl = (TH1D *)qcdPtMt40PlusBase->Clone("qcd_pt_mt40_incl");
            qcd_pt_mt40_incl->SetDirectory(nullptr);
            qcd_pt_mt40_incl->Add(qcdPtMt40MinusBase);
            std::cout << "[QCD] m_T>40 lepton-pT QCD: W+ = "
                      << Form("%.1f", qcdPtMt40PlusBase->Integral()) << ", W- = "
                      << Form("%.1f", qcdPtMt40MinusBase->Integral())
                      << " events (vs no-m_T-cut "
                      << Form("%.1f", qcdPtPlusBase ? qcdPtPlusBase->Integral() : 0.0) << " + "
                      << Form("%.1f", qcdPtMinusBase ? qcdPtMinusBase->Integral() : 0.0) << ")\n";
        }
        else
            std::cerr << "[WARN] m_T>40 lepton-pT QCD template (qcd_pt_mt40_*) not in "
                      << qcdFile << "; re-run qcd_abcd.C on a skim with the"
                      << " h_iso_pt_mt40 planes. m_T>40 stacks will omit QCD.\n";
        if (qcdMetMt40PlusBase && qcdMetMt40MinusBase)
        {
            qcd_met_mt40_incl = (TH1D *)qcdMetMt40PlusBase->Clone("qcd_met_mt40_incl");
            qcd_met_mt40_incl->SetDirectory(nullptr);
            qcd_met_mt40_incl->Add(qcdMetMt40MinusBase);
            std::cout << "[QCD] m_T>40 PF-MET QCD (met_mt40): W+ = "
                      << Form("%.1f", qcdMetMt40PlusBase->Integral()) << ", W- = "
                      << Form("%.1f", qcdMetMt40MinusBase->Integral()) << " events\n";
        }
        else
            std::cerr << "[WARN] m_T>40 PF-MET QCD template (qcd_met_mt40_*) not in " << qcdFile
                      << "; re-run correction/run_qcd_abcd.sh on a skim with the h_iso_met_mt40 planes."
                      << " The met_mt40 stacks omit QCD and its Combine input gets no qcd_abcd.\n";
    }

    // Labeled-bin reader for abcd_counts_* (loop the labels and string-compare;
    // TAxis::FindBin(label) would APPEND a new label to the mutable axis on a
    // miss instead of failing).
    auto countBin = [](TH1D *h, const char *lab, double *err = nullptr) -> double
    {
        if (!h) return 0.0;
        for (int i = 1; i <= h->GetNbinsX(); ++i)
            if (strcmp(h->GetXaxis()->GetBinLabel(i), lab) == 0)
            {
                if (err) *err = h->GetBinError(i);
                return h->GetBinContent(i);
            }
        std::cerr << "[WARN] abcd_counts: no bin labeled '" << lab << "'\n";
        return 0.0;
    };

    // In-fit ABCD SR template bases: clone of qcd_pt_mt40_* rescaled so the
    // total is EXACTLY A0 = qcdB0*qcdC0/qcdD0 (then the card's functional
    // rateParam (sB*sC/sD) needs no baked constants and Asimov closes at 1).
    // Staleness guard: the template's own sideband count (total / T_MET) must
    // match the exported C40 to ~1% (pT>100 overflow + negative-bin clamp make
    // it not exact) -- a >2% mismatch means qcd_abcd_<lep>.root mixes objects
    // from different qcd_abcd.C runs.
    {
        auto mkAbcdBase = [&](TH1D *base, TH1D *counts, const char *nm,
                              const char *cg) -> TH1D *
        {
            if (!base || !counts) return nullptr;
            const double b0 = countBin(counts, "qcdB0"), c0 = countBin(counts, "qcdC0"),
                         d0 = countBin(counts, "qcdD0"), tmet = countBin(counts, "T_met");
            const double v = base->Integral();
            if (b0 <= 0.0 || c0 <= 0.0 || d0 <= 0.0 || v <= 0.0) return nullptr;
            const double A0 = b0 * c0 / d0;
            if (tmet > 0.0 && std::fabs(v / tmet / c0 - 1.0) > 0.02)
                std::cerr << Form("[WARN] qcd_abcd staleness (%s): template sideband %.1f vs"
                                  " exported C40 %.1f -- re-run correction/run_qcd_abcd.sh\n",
                                  cg, v / tmet, c0);
            TH1D *h = (TH1D *)base->Clone(nm);
            h->SetDirectory(nullptr);
            h->Scale(A0 / v);
            std::cout << Form("[QCD-ABCD] %s: B0=%.1f C40=%.1f D0=%.1f -> A0=%.1f"
                              "  (rescale x%.3f from the T_MET-normalized total %.1f)\n",
                              cg, b0, c0, d0, A0, A0 / v, v);
            return h;
        };
        if (!isElec)
        {
            qcdAbcdPlusBase  = mkAbcdBase(qcdPtMt40PlusBase,  abcdCountsPlus,  "qcd_abcd_Plus_base",  "W+");
            qcdAbcdMinusBase = mkAbcdBase(qcdPtMt40MinusBase, abcdCountsMinus, "qcd_abcd_Minus_base", "W-");
            // met_mt40: the SAME A0 (same events, same counts), MET shape
            qcdAbcdMetPlusBase  = mkAbcdBase(qcdMetMt40PlusBase,  abcdCountsPlus,  "qcd_abcd_met_Plus_base",  "W+ (met_mt40)");
            qcdAbcdMetMinusBase = mkAbcdBase(qcdMetMt40MinusBase, abcdCountsMinus, "qcd_abcd_met_Minus_base", "W- (met_mt40)");
        }
        else
        {
            // Electrons: the FF-weighted template IS the base (total = the FF
            // prediction A0ff, printed next to the flat A0 it replaces). It is
            // taken only if it is complete (nominal + both slope twins) and
            // matches the A0ff exported with the SAME qcd_abcd.C run -- a stale
            // or missing FF template must never silently fall back to flat T.
            // `lab` = the abcd_counts label of its total: A0ff (lepton pT) or,
            // for the met_mt40 MET template (2026-10-20), A0ffMet.
            auto mkFFBase = [&](TH1D *ff, TH1D *up, TH1D *dn, TH1D *counts, const char *cg,
                                const char *lab = "A0ff") -> bool {
                if (!ff || !up || !dn || !counts) return false;
                const double a0ff = countBin(counts, lab);
                const double b0 = countBin(counts, "qcdB0"), c0 = countBin(counts, "qcdC0"), d0 = countBin(counts, "qcdD0");
                const double v = ff->Integral();
                if (a0ff <= 0.0 || std::fabs(v / a0ff - 1.0) > 1e-6)
                {
                    std::cerr << Form("[ERROR] qcd_abcd FF template (%s): total %.3f vs exported %s %.3f"
                                      " -- stale qcd_abcd_ele.root, re-run correction/run_qcd_abcd.sh ele\n",
                                      cg, v, lab, a0ff);
                    return false;
                }
                std::cout << Form("[QCD-ABCD] %s: FAKE-FACTOR template, total %s = %.1f (flat A0 = B0*C40/D0 ="
                                  " %.1f -> FF/flat %.3f); slope twins %.1f / %.1f (%+.1f%% / %+.1f%%)\n",
                                  cg, lab, v, d0 > 0 ? b0 * c0 / d0 : 0.0, (d0 > 0 && b0 * c0 > 0) ? v / (b0 * c0 / d0) : 0.0,
                                  up->Integral(), dn->Integral(), 100.0 * (up->Integral() / v - 1.0),
                                  100.0 * (dn->Integral() / v - 1.0));
                return true;
            };
            if (mkFFBase(ffPlus, ffPlusUp, ffPlusDn, abcdCountsPlus, "W+") &&
                mkFFBase(ffMinus, ffMinusUp, ffMinusDn, abcdCountsMinus, "W-"))
            {
                qcdAbcdPlusBase  = (TH1D *)ffPlus->Clone("qcd_abcd_Plus_base");  qcdAbcdPlusBase->SetDirectory(nullptr);
                qcdAbcdMinusBase = (TH1D *)ffMinus->Clone("qcd_abcd_Minus_base"); qcdAbcdMinusBase->SetDirectory(nullptr);
                qcdAbcdPlusUpBase = ffPlusUp;   qcdAbcdPlusDnBase = ffPlusDn;
                qcdAbcdMinusUpBase = ffMinusUp; qcdAbcdMinusDnBase = ffMinusDn;
            }
            else
                std::cerr << "[ERROR] electron fake-factor SR template (qcd_pt_mt40_ff_ele*) missing, incomplete"
                          << " or stale in " << qcdFile << " -- qcd_abcd + CR templates will NOT be written.\n";
            if (mkFFBase(ffMetPlus, ffMetPlusUp, ffMetPlusDn, abcdCountsPlus, "W+ (met_mt40)", "A0ffMet") &&
                mkFFBase(ffMetMinus, ffMetMinusUp, ffMetMinusDn, abcdCountsMinus, "W- (met_mt40)", "A0ffMet"))
            {
                qcdAbcdMetPlusBase  = (TH1D *)ffMetPlus->Clone("qcd_abcd_met_Plus_base");  qcdAbcdMetPlusBase->SetDirectory(nullptr);
                qcdAbcdMetMinusBase = (TH1D *)ffMetMinus->Clone("qcd_abcd_met_Minus_base"); qcdAbcdMetMinusBase->SetDirectory(nullptr);
                qcdAbcdMetPlusUpBase = ffMetPlusUp;   qcdAbcdMetPlusDnBase = ffMetPlusDn;
                qcdAbcdMetMinusUpBase = ffMetMinusUp; qcdAbcdMetMinusDnBase = ffMetMinusDn;
            }
            else
                std::cerr << "[ERROR] electron fake-factor MET template (qcd_met_mt40_ff_ele*) missing, incomplete"
                          << " or stale in " << qcdFile << " -- the met_mt40 qcd_abcd + CR templates will NOT be written.\n";
        }
        if (!qcdAbcdPlusBase || !qcdAbcdMinusBase)
            std::cerr << "[WARN] abcd_counts_* missing/incomplete in " << qcdFile
                      << " -- re-run correction/run_qcd_abcd.sh. qcd_abcd + CR templates"
                      << " will NOT be written (QCD_MODE=abcd unavailable for this input).\n";
        if (!qcdAbcdMetPlusBase || !qcdAbcdMetMinusBase)
            std::cerr << "[WARN] no met_mt40 in-fit QCD base -- combine_input_W_met_mt40.root gets no"
                      << " qcd_abcd / CR templates (QCD_MODE=abcd unavailable for --disc met_mt40).\n";
    }
    // per-variant views of the in-fit bases [var][0 = W+, 1 = W-] (null = absent)
    TH1D *abcdBaseV[kNVar][2]   = {{qcdAbcdPlusBase, qcdAbcdMinusBase}, {qcdAbcdMetPlusBase, qcdAbcdMetMinusBase}};
    TH1D *abcdUpBaseV[kNVar][2] = {{qcdAbcdPlusUpBase, qcdAbcdMinusUpBase}, {qcdAbcdMetPlusUpBase, qcdAbcdMetMinusUpBase}};
    TH1D *abcdDnBaseV[kNVar][2] = {{qcdAbcdPlusDnBase, qcdAbcdMinusDnBase}, {qcdAbcdMetPlusDnBase, qcdAbcdMetMinusDnBase}};
    TH1D *dispBaseV[kNVar][2]   = {{qcdPtMt40PlusBase, qcdPtMt40MinusBase}, {qcdMetMt40PlusBase, qcdMetMt40MinusBase}};
    TH1D *dispInclV[kNVar]      = {qcd_pt_mt40_incl, qcd_met_mt40_incl};

    // Per-bin low-MET excess (data - EWK at MET < qcdMetCut, clamped >= 0) of
    // the lab (suf "") or FB (suf "_FB") rapidity bins.
    auto lowMetExcess = [&](const char *chg, const char *suf) -> std::vector<double>
    {
        std::vector<double> e(NY, 0.0);
        for (int iy = 0; iy < NY; ++iy)
        {
            auto lowInt = [&](TFile *ff, double k) -> double
            {
                TH1D *hh = (TH1D *)ff->Get(Form("h_met_%s_y%d%s", chg, iy, suf));
                if (!hh) return 0.0;
                const int bc = hh->GetXaxis()->FindBin(qcdMetCut - 1e-6);
                return k * hh->Integral(1, bc); // read-only: no clone needed
            };
            double d = lowInt(f, 1.0); // data (k=1)
            double ewk = lowInt(f_Wp, k_Wp) + lowInt(f_Wm, k_Wm) + lowInt(f_DY, k_DY) +
                         lowInt(f_DYtau, k_DYtau) + lowInt(f_Wptau, k_Wptau) + lowInt(f_Wmtau, k_Wmtau) +
                         lowInt(f_TT, k_TT);
            e[iy] = d - ewk;
            if (e[iy] < 0) e[iy] = 0; // EWK over-subtraction in a sparse bin -> no QCD there
        }
        return e;
    };
    // Per-y QCD weights = the bin's share of the low-MET excess of the WHOLE
    // acceptance, i.e. always divided by the LAB sum (the lab bins tile
    // |y| < 2.4, where the inclusive template lives). Lab weights therefore sum
    // to 1. The FB bins tile only [-1.7068, 2.4], so the FB weights sum to the
    // window's share f = Sum_FB / Sum_lab and the FB templates to f x the
    // inclusive total. (Until 2026-10-19 the FB weights were normalized to 1
    // over the FB bins, which pushed the multijet of y < -1.7068 into the
    // window: x1.10-1.14 (mu) / x1.24-1.34 (e) too much multijet in the fb
    // fit -- found by correction/qcd_ysplit.C.)
    auto qcdWeights = [&](const char *chg, const char *suf) -> std::vector<double>
    {
        std::vector<double> w = lowMetExcess(chg, suf);
        const std::vector<double> eLab = lowMetExcess(chg, "");
        double sum = 0.0;
        for (int iy = 0; iy < NY; ++iy) sum += eLab[iy];
        for (int iy = 0; iy < NY; ++iy) w[iy] = (sum > 0) ? w[iy] / sum : 1.0 / NY;
        return w;
    };

    const std::vector<double> wQcdWp   = qcdPlusBase  ? qcdWeights("Wp", "")    : std::vector<double>(NY, 0.0);
    const std::vector<double> wQcdWm   = qcdMinusBase ? qcdWeights("Wm", "")    : std::vector<double>(NY, 0.0);
    const std::vector<double> wQcdWpFB = qcdPlusBase  ? qcdWeights("Wp", "_FB") : std::vector<double>(NY, 0.0);
    const std::vector<double> wQcdWmFB = qcdMinusBase ? qcdWeights("Wm", "_FB") : std::vector<double>(NY, 0.0);
    if (qcdPlusBase && qcdMinusBase)
    {
        double sp = 0, sm = 0, spFB = 0, smFB = 0;
        for (int iy = 0; iy < NY; ++iy)
        {
            sp += wQcdWp[iy]; sm += wQcdWm[iy]; spFB += wQcdWpFB[iy]; smFB += wQcdWmFB[iy];
        }
        std::cout << Form("[QCD] per-y split weights sum: lab W+ %.4f W- %.4f; FB window share f: W+ %.4f W- %.4f"
                          " (FB templates = f x the inclusive total)\n", sp, sm, spFB, smFB);
    }

    // --- abcdy: the per-y factors of the qcd_abcdy template (2026-10-19) -----
    // QCD_MODE=abcdy in the fork (OPT-IN) runs the in-fit ABCD per rapidity
    // GROUP (correction/abcd_ygroups.h). Its SR template, per y bin of group g:
    //   qcd_abcdy_y = qcdAbcdBase x (A_g / A0) x w_in,y
    //   A0     = qcdB0*qcdC0/qcdD0 of the INCLUSIVE counts. The base is today's
    //            qcd_abcd base (total A0; the electron one is fake-factor
    //            weighted, total A0ff), so a group carries A_g (muon) or
    //            A_g x A0ff/A0 (electron: the FF weighting stays inclusive).
    //   A_g    = the group's own B_g*C40_g/D_g (abcd_counts_*_<bin>_<g>, r = 1).
    //   w_in,y = the bin's low-MET excess / the excess of its GROUP, from the
    //            same binning: the qcdWeights proxy, now only INSIDE a group.
    // In the fit each group gets its own sB_g*sC_g/sD_g. Built only when every
    // group count is present and positive and the lab groups' DATA close
    // exactly on the inclusive counts (= the same qcd_abcd.C run); otherwise
    // abcdyOk = false: no qcd_abcdy, no group CR dirs, no group sidecar (any
    // stale sidecar is deleted, so the fork cannot build abcdy cards).
    bool abcdyOk = (qcdAbcdPlusBase && qcdAbcdMinusBase);
    double abcdyFac[2][2][NY] = {};                   // [charge][lab|fb][y]: the factor on qcdAbcdBase
    double abcdyGrpFac[2][2][pOAbcdY::kNGroups] = {}; // A_g / A0
    if (abcdyOk)
    {
        for (int ic = 0; ic < 2 && abcdyOk; ++ic)
        {
            const char *cg = (ic == 0) ? "Wp" : "Wm";
            TH1D *inc = (ic == 0) ? abcdCountsPlus : abcdCountsMinus;
            TH1D *base = (ic == 0) ? qcdAbcdPlusBase : qcdAbcdMinusBase;
            const double b0 = countBin(inc, "qcdB0"), c0 = countBin(inc, "qcdC0"), d0 = countBin(inc, "qcdD0");
            const double A0 = (d0 > 0.0) ? b0 * c0 / d0 : 0.0;
            if (A0 <= 0.0) { abcdyOk = false; break; }
            for (int ib = 0; ib < 2; ++ib)
            {
                const std::vector<double> ex = lowMetExcess(cg, ib == 0 ? "" : "_FB");
                const std::vector<double> &wOld = (ic == 0) ? (ib == 0 ? wQcdWp : wQcdWpFB)
                                                            : (ib == 0 ? wQcdWm : wQcdWmFB);
                double sB = 0, sC = 0, sD = 0;
                for (int g = 0; g < pOAbcdY::kNGroups; ++g)
                {
                    TH1D *gc = abcdyCounts[ic][ib][g];
                    if (!gc)
                    {
                        std::cerr << "[WARN] abcdy: abcd_counts_" << qlep << (ic == 0 ? "Plus_" : "Minus_") << abcdyBin[ib]
                                  << "_" << pOAbcdY::kGroupName[g] << " not in " << qcdFile
                                  << " -- re-run correction/run_qcd_abcd.sh; no abcdy objects written\n";
                        abcdyOk = false;
                        continue;
                    }
                    double eB = 0, eC = 0, eD = 0;
                    const double bg = countBin(gc, "qcdB0", &eB), cgr = countBin(gc, "qcdC0", &eC),
                                 dg = countBin(gc, "qcdD0", &eD);
                    sB += countBin(gc, "dataB"); sC += countBin(gc, "dataC40"); sD += countBin(gc, "dataD");
                    if (bg <= 0.0 || cgr <= 0.0 || dg <= 0.0)
                    {
                        std::cerr << Form("[WARN] abcdy %s %s %s: QCD-only count <= 0 (B %.1f C40 %.1f D %.1f) -- no "
                                          "abcdy objects written\n", cg, abcdyBin[ib], pOAbcdY::kGroupName[g], bg, cgr, dg);
                        abcdyOk = false;
                        continue;
                    }
                    const double Ag = bg * cgr / dg;
                    const double eAg = Ag * std::sqrt(eB * eB / (bg * bg) + eC * eC / (cgr * cgr) + eD * eD / (dg * dg));
                    abcdyGrpFac[ic][ib][g] = Ag / A0;
                    const std::vector<int> gb = pOAbcdY::BinsOf(g);
                    double es = 0.0, wo = 0.0;
                    for (int iy : gb) { es += ex[iy]; wo += wOld[iy]; }
                    for (int iy : gb)
                        abcdyFac[ic][ib][iy] = (Ag / A0) * ((es > 0.0) ? ex[iy] / es : 1.0 / gb.size());
                    std::cout << Form("[QCD-ABCDY] %s %-3s %s (y%d-y%d): B0 %.1f C40 %.1f D0 %.1f -> A_g = %.1f +- %.1f"
                                      " = %.3f x A0; qcd_abcdy %.1f +- %.1f vs qcd_abcd %.1f in these bins (x%.3f)\n",
                                      cg, abcdyBin[ib], pOAbcdY::kGroupName[g], gb.front(), gb.back(), bg, cgr, dg,
                                      Ag, eAg, Ag / A0, base->Integral() * Ag / A0, base->Integral() * eAg / A0,
                                      base->Integral() * wo, wo > 0.0 ? (Ag / A0) / wo : 0.0);
                }
                // staleness: the lab groups tile the inclusive plane -> their
                // DATA counts (integers) must add up exactly
                if (ib == 0 && abcdyOk &&
                    (sB != countBin(inc, "dataB") || sC != countBin(inc, "dataC40") || sD != countBin(inc, "dataD")))
                {
                    std::cerr << Form("[ERROR] abcdy %s: the lab groups' data (B %.0f C40 %.0f D %.0f) do not add up to the "
                                      "inclusive counts (%.0f %.0f %.0f) -- stale %s, re-run correction/run_qcd_abcd.sh\n",
                                      cg, sB, sC, sD, countBin(inc, "dataB"), countBin(inc, "dataC40"),
                                      countBin(inc, "dataD"), qcdFile.c_str());
                    abcdyOk = false;
                }
            }
        }
        std::cout << "[QCD-ABCDY] rapidity-resolved in-fit ABCD inputs "
                  << (abcdyOk ? "COMPLETE -> qcd_abcdy + group CR dirs + the group sidecar will be written"
                              : "INCOMPLETE -> no abcdy objects (QCD_MODE=abcdy unavailable for this input)")
                  << "\n";
    }

    // abcdy per selection variant: the group factors above are COUNT-level and
    // shared; a variant needs its own in-fit base on top
    const bool abcdyOkV[kNVar] = {abcdyOk && abcdBaseV[kVarMt40][0] && abcdBaseV[kVarMt40][1],
                                  abcdyOk && abcdBaseV[kVarMetMt40][0] && abcdBaseV[kVarMetMt40][1]};

    if (qcd_met_incl)
    {
        std::cout << "[QCD] inclusive ABCD QCD in MET stack: W+ = "
                  << Form("%.1f", qcdPlusBase->Integral()) << ", W- = "
                  << Form("%.1f", qcdMinusBase->Integral()) << " events\n";
        std::cout << "[QCD] per-y split (std binning) W+ yields:";
        for (int iy = 0; iy < NY; ++iy)
            std::cout << " " << Form("%.0f", wQcdWp[iy] * qcdPlusBase->Integral());
        std::cout << "\n";
    }

    // Per-y QCD histo helper (inclusive template shape x per-bin weight).
    auto qcdPerY = [&](TH1D *base, double w, const char *name) -> TH1D *
    {
        if (!base) return nullptr;
        TH1D *h = (TH1D *)base->Clone(name);
        h->SetDirectory(nullptr);
        h->Scale(w);
        return h;
    };

    // (the selection-variant table kVarMt40 / kVarMetMt40 is declared above,
    // before the QCD templates are read)

    // --- Structured per-region Combine input -----------------------------------
    // For every fit region we emit the 7 absolute templates the downstream
    // Combine datacards consume -- signal, z, ztau, wtau, ttbar, qcd, data_obs
    // (ttbar since 2026-10-16) -- each in its own TDirectory. Regions:
    //   Wp_lab_y{0..11} / Wm_lab_y{0..11}   standard lab-frame rapidity bins
    //   Wp_fb_y{0..11}  / Wm_fb_y{0..11}    FB-symmetric bins (forward/backward)
    //   Wp_incl / Wm_incl                   per-charge, summed over the lab bins
    //   W_incl                              both charges (inclusive W fit / Z-control link)
    // "signal" = both W MC samples in that reco-charge region (Wp- and Wm-sample
    // events reconstructed with the region's charge); "wtau" = Wptau+Wmtau. Every
    // template is already at its absolute pO yield (k_s applied below) -- NO area
    // normalization, so the fit floats real normalizations, not shapes-to-data.
    //
    // TWO input files, one per W discriminant (three until 2026-09-21, when the
    // plain-leppt one was retired -- see the variant table above), same region
    // names and template roles in both so the fork's datacard generator applies
    // unchanged:
    //   combine_input_W.root            PF MET               (the backup fit)
    //   combine_input_W_leppt_mt40.root lepton pT, pT>25 && m_T>40 (PRIMARY)
    // NB for the lepton-pT variant the QCD-anchoring low-MET region is absent
    // from the discriminant, which is why its QCD is constrained by the in-fit
    // ABCD control regions (QCD_MODE=abcd) instead of by the shape.
    struct RegionTemplates {
        std::string dir;
        TH1D *data, *sig, *z, *ztau, *wtau, *qcd;
        // 7th template (leppt_mt40 only): the in-fit-ABCD-normalized QCD
        // (same shape, total = A0). NSDMI required -- accumulate()/W_incl
        // construct RegionTemplates without makeRegion.
        TH1D *qcdAbcd = nullptr;
        // 8th SR template (leppt_mt40 only, 2026-10-19): the rapidity-resolved
        // in-fit ABCD's qcd_abcdy = the qcd_abcd shape with each rapidity
        // group carrying its own ABCD total (QCD_MODE=abcdy; see abcdyFac)
        TH1D *qcdAbcdy = nullptr;
        // the ttbar background (2026-10-16), absolute MC like z/ztau/wtau,
        // in EVERY region of both inputs (NSDMI for the same reason)
        TH1D *ttbar = nullptr;
        // Shape systematics (2026-09-07 LHE, 2026-09-14 muon SFs): "<process>_<syst>Up|Down"
        // -> the absolute varied template, process in {signal, z, ztau, wtau, ttbar},
        // syst in `systNames` below (nPDF, qcdScale, alphaS from skim/lhe_updown.py;
        // muID, muIso, muTrig from skim.C, muon channel only), built by
        // attachSysts() from the <hist>_<syst>Up/Down twins in every MC skim
        // file. data_obs and the data-driven qcd/qcd_abcd have none.
        std::map<std::string, TH1D *> syst;
    };
    std::vector<RegionTemplates> regions;          // PF MET (the nominal fit input)
    std::vector<RegionTemplates> regionsPt[kNVar]; // lepton-pT: [kVarMt40]

    // `tag` disambiguates the detached-clone names across the three collections
    // (same region dir appears once per discriminant).
    auto cloneDetached = [](TH1D *h, const char *nm) -> TH1D * {
        if (!h) return nullptr;
        TH1D *c = (TH1D *)h->Clone(nm);
        c->SetDirectory(nullptr);
        return c;
    };
    auto sum2 = [&](TH1D *a, TH1D *b, const char *nm) -> TH1D * {
        TH1D *r = cloneDetached(a, nm);
        if (r) { if (b) r->Add(b); }
        else if (b) r = cloneDetached(b, nm);
        return r;
    };
    auto makeRegion = [&](const std::string &dir, const std::string &tag,
                          TH1D *data, TH1D *sigA, TH1D *sigB,
                          TH1D *z, TH1D *ztau, TH1D *wtA, TH1D *wtB, TH1D *tt,
                          TH1D *qcd, TH1D *qcdAbcd = nullptr) -> RegionTemplates {
        RegionTemplates r;
        r.dir  = dir;
        r.data = cloneDetached(data, (dir + tag + "_data_obs").c_str());
        r.sig  = sum2(sigA, sigB,    (dir + tag + "_signal").c_str());
        r.z    = cloneDetached(z,    (dir + tag + "_z").c_str());
        r.ztau = cloneDetached(ztau, (dir + tag + "_ztau").c_str());
        r.wtau = sum2(wtA, wtB,      (dir + tag + "_wtau").c_str());
        r.ttbar = cloneDetached(tt,  (dir + tag + "_ttbar").c_str());
        r.qcd  = cloneDetached(qcd,  (dir + tag + "_qcd").c_str());
        r.qcdAbcd = cloneDetached(qcdAbcd, (dir + tag + "_qcd_abcd").c_str());
        return r;
    };

    // --- Lepton-pT stacks -------------------------------------------------------
    // Absolute data-vs-stack comparisons in the leading-lepton pT (the pT twins
    // of the MET stacks: h_leppt_mt40_* from the skim + the qcd_pt_mt40_* ABCD
    // template). The PLOTS are display only; the same per-region templates are
    // collected here and written to combine_input_W_leppt_mt40.root, which IS
    // the primary fit input (combine_input_W.root itself stays PF-MET-only).
    //
    // ONE variant since 2026-09-21 (see the variant table above):
    //   kVarMt40 -- the lepton-pT-DISCRIMINANT selection pT>25 && m_T>40,
    //               h_leppt_mt40_* + qcd_pt_mt40_* (2026-07-30). The m_T cut is
    //               what suppresses QCD when the fit no longer uses the MET
    //               shape to do it.
    TH1D *h_leppt_inclusive[kNVar]           = {nullptr};
    TH1D *h_leppt_inclusive_MC_signal[kNVar] = {nullptr};
    TH1D *h_leppt_inclusive_MC_Z[kNVar]      = {nullptr};
    TH1D *h_leppt_inclusive_MC_Ztau[kNVar]   = {nullptr};
    TH1D *h_leppt_inclusive_MC_Wtau[kNVar]   = {nullptr};
    TH1D *h_leppt_inclusive_MC_TT[kNVar]     = {nullptr};

    // Get + absolute-scale one MC histo (fresh Get per name/file, scaled in
    // place and used once -- same convention as scaleAll above).
    auto getScaled = [](TFile *ff, double k, const char *name) -> TH1D *
    {
        TH1D *h = (TH1D *)ff->Get(name);
        if (h) h->Scale(k);
        return h;
    };

    // The shape systematics this channel's MC skim files carry: the LHE families
    // (both flavours; skim/lhe_updown.py) + the lepton-SF family of THIS
    // flavour (muon: muSF, skim/muon_sf.h; electron since 2026-09-29: eleSF,
    // skim/electron_sf.h -- both written by skim.C). Single list for
    // attachSysts() and the sidecar below. (The electron qcdFF shape nuisance
    // sits on the data-driven qcd_abcd template and is handled separately.)
    std::vector<std::string> systNames(pOLhe::kLheSystNames, pOLhe::kLheSystNames + pOLhe::kNLheSysts);
    if (!isElec)
        systNames.insert(systNames.end(), pOSF::kMuonSFSystNames, pOSF::kMuonSFSystNames + pOSF::kNMuonSFSysts);
    else
        systNames.insert(systNames.end(), pOSF::kEleSFSystNames, pOSF::kEleSFSystNames + pOSF::kNEleSFSysts);
    // Per selection variant: the met_mt40 input additionally carries the MET-
    // resolution nuisance of this flavour (metResMu / metResEle, skim/met_res.h;
    // user decision 2026-10-20: the MET fit only -- the leppt_mt40 twins exist
    // in the skim but are NOT attached, so that input is unchanged).
    std::vector<std::string> systNamesV[kNVar] = {systNames, systNames};
    systNamesV[kVarMetMt40].push_back(pOMetRes::SystName(isElec));

    // Shape systematics for one region (2026-09-07). For every syst in
    // `systNames` and direction Up/Down, read the per-sample
    // "<baseName>_<syst><dir>" twins (skim/lhe_updown.py for the LHE families,
    // skim.C for the muon SFs) from the MC skim files, scale each by its own
    // k_s and combine them EXACTLY like the nominal templates: signal = W+
    // sample + W- sample reconstructed with this region's charge, wtau = both
    // tau samples, z / ztau / ttbar single samples. A syst is attached only when all
    // seven inputs exist in BOTH directions (otherwise [WARN] once per syst and
    // skipped); the sidecar written next to each Combine input lists what
    // every region carries, and the fork's card generator must not reference
    // more. Each varied name is Get()'d exactly once per file (region names
    // are unique), so the in-place Scale is safe.
    std::set<std::string> lheWarned;
    auto attachSysts = [&](RegionTemplates &r, const std::string &tag, const char *baseName,
                           const std::vector<std::string> *names = nullptr)
    {
        for (const std::string &systStr : (names ? *names : systNames))
        {
            const char *syst = systStr.c_str();
            std::map<std::string, TH1D *> got;
            bool complete = true;
            for (const char *dir : {"Up", "Down"})
            {
                const std::string nm = Form("%s_%s%s", baseName, syst, dir);
                TH1D *hWp  = getScaled(f_Wp,    k_Wp,    nm.c_str());
                TH1D *hWm  = getScaled(f_Wm,    k_Wm,    nm.c_str());
                TH1D *hZ   = getScaled(f_DY,    k_DY,    nm.c_str());
                TH1D *hZt  = getScaled(f_DYtau, k_DYtau, nm.c_str());
                TH1D *hWpt = getScaled(f_Wptau, k_Wptau, nm.c_str());
                TH1D *hWmt = getScaled(f_Wmtau, k_Wmtau, nm.c_str());
                TH1D *hTT  = getScaled(f_TT,    k_TT,    nm.c_str());
                if (!hWp || !hWm || !hZ || !hZt || !hWpt || !hWmt || !hTT) { complete = false; break; }
                const std::string key = std::string("_") + syst + dir; // e.g. "_nPDFUp"
                const std::string pre = r.dir + tag;
                got["signal" + key] = sum2(hWp, hWm,   (pre + "_signal" + key).c_str());
                got["z" + key]      = cloneDetached(hZ,  (pre + "_z" + key).c_str());
                got["ztau" + key]   = cloneDetached(hZt, (pre + "_ztau" + key).c_str());
                got["wtau" + key]   = sum2(hWpt, hWmt, (pre + "_wtau" + key).c_str());
                got["ttbar" + key]  = cloneDetached(hTT, (pre + "_ttbar" + key).c_str());
            }
            if (!complete)
            {
                if (lheWarned.insert(syst).second)
                    std::cerr << "[WARN] shape syst " << syst << ": " << baseName << "_" << syst
                              << "Up/Down missing in at least one MC skim file -> not attached"
                              << " (LHE families: run skim/run_lhe_updown.sh after the skim; muon SFs: re-skim)\n";
                continue;
            }
            r.syst.insert(got.begin(), got.end());
        }
    };

    // One lepton-pT stack region. accumulate=true (lab bins) also feeds the
    // inclusive accumulators, mirroring how h_met_inclusive* are built.
    // qcdAbcdUp/Dn + qcdSyst (electrons, 2026-09-29): the per-y slope twins of
    // the FF-weighted qcd_abcd, attached as <qcd_abcd>_<qcdSyst>Up/Down.
    auto lepPtStack = [&](int var, int iy, const char *chg, const char *suf,
                          TH1D *qcdH, const char *sub1, const std::string &sub2,
                          const std::string &outPath, bool accumulate,
                          TH1D *qcdAbcdH = nullptr, TH1D *qcdAbcdUp = nullptr,
                          TH1D *qcdAbcdDn = nullptr, const std::string &qcdSyst = "",
                          TH1D *qcdAbcdyH = nullptr, TH1D *qcdAbcdyUp = nullptr,
                          TH1D *qcdAbcdyDn = nullptr) -> bool
    {
        const std::string nm = Form("%s_%s_y%d%s", varStem[var], chg, iy, suf);
        TH1D *hD = (TH1D *)f->Get(nm.c_str());
        if (!hD) return false;
        TH1D *hWp    = getScaled(f_Wp,    k_Wp,    nm.c_str());
        TH1D *hWm    = getScaled(f_Wm,    k_Wm,    nm.c_str());
        TH1D *hZ     = getScaled(f_DY,    k_DY,    nm.c_str());
        TH1D *hZtau  = getScaled(f_DYtau, k_DYtau, nm.c_str());
        TH1D *hWptau = getScaled(f_Wptau, k_Wptau, nm.c_str());
        TH1D *hWmtau = getScaled(f_Wmtau, k_Wmtau, nm.c_str());
        TH1D *hTT    = getScaled(f_TT,    k_TT,    nm.c_str());

        std::vector<TH1 *> bkgs;
        std::vector<std::string> names;
        auto pushIf = [&](TH1 *h, const char *n) { if (h) { bkgs.push_back(h); names.push_back(n); } };
        pushIf(hWp, "W+");
        pushIf(hWm, "W-");
        pushIf(hZ, "DY");
        pushIf(hZtau, "DY Tau");
        pushIf(hWptau, "Wp Tau");
        pushIf(hWmtau, "Wm Tau");
        pushIf(qcdH, "QCD");
        addTT(bkgs, names, hTT);

        std::vector<std::string> box = {
            Form("Passing Events: %.0f", hD->Integral(1, hD->GetNbinsX())),
            Form("W signal MC: %.0f", sigInt({hWp, hWm})),
            totalLine(bkgs, names)};
        for (const std::string &l : prefitLines(hD, {hWp, hWm, hWptau, hWmtau}, bkgs))
            box.push_back(l);

        // pT axes start at the selection floor (bin edge 24 encloses the
        // 25 GeV cut on the 2 GeV grid) -- no empty [0,25) band; the MET
        // variant keeps its full axis (varXLo/Hi). The helper restores the
        // full range after saving, so the Combine-input writing below sees
        // untouched histograms. Both variants have 2 GeV bins.
        PlotStyle psPt = psLepPt;
        if (varXHi[var] > varXLo[var])
        {
            psPt.xRangeLo = varXLo[var];
            psPt.xRangeHi = varXHi[var];
        }
        SaveNicePlot1D_WithBkg(hD, bkgs, names, outPath,
                               varXTitle[var], "Events / 2.0 GeV", "",
                               sub1, sub2, box, psPt, commonTuner);

        // Collect this region's templates for the variant's Combine input
        // (combine_input_W_<leppt_mt40|met_mt40>.root below -- combine_input_W.root
        // itself stays the no-m_T-cut PF-MET one).
        regionsPt[var].push_back(makeRegion(
            Form("%s_%s_y%d", chg, (suf[0] ? "fb" : "lab"), iy),
            Form("_lp%s", varTag[var]),
            hD, hWp, hWm, hZ, hZtau, hWptau, hWmtau, hTT, qcdH, qcdAbcdH));
        // LHE / lepton-SF (/ met_mt40: MET-resolution) shape systematics: the
        // skim stores twins for h_met_*, h_met_mt40_* and h_leppt_mt40_* only --
        // which is also why the plain leppt variant was retired (2026-09-21).
        attachSysts(regionsPt[var].back(), Form("_lp%s", varTag[var]), nm.c_str(), &systNamesV[var]);
        // The electron fake-factor slope nuisance (one per charge): the slope
        // +-1 sigma twins of the FF-weighted qcd_abcd, split in y exactly like
        // the nominal. Stored under Combine's <process>_$SYSTEMATIC naming for
        // the qcd_abcd histogram; the sidecar lists it for that charge only.
        if (!qcdSyst.empty() && qcdAbcdH && qcdAbcdUp && qcdAbcdDn)
        {
            RegionTemplates &rr = regionsPt[var].back();
            const std::string pre = rr.dir + Form("_lp%s", varTag[var]);
            rr.syst["qcd_abcd_" + qcdSyst + "Up"]   = cloneDetached(qcdAbcdUp, (pre + "_qcd_abcd_" + qcdSyst + "Up").c_str());
            rr.syst["qcd_abcd_" + qcdSyst + "Down"] = cloneDetached(qcdAbcdDn, (pre + "_qcd_abcd_" + qcdSyst + "Down").c_str());
        }
        // the abcdy SR template (2026-10-19) + its electron slope twins (same
        // per-y factor as the nominal, so the slope moves every group alike)
        if (qcdAbcdyH)
        {
            RegionTemplates &rr = regionsPt[var].back();
            const std::string pre = rr.dir + Form("_lp%s", varTag[var]);
            rr.qcdAbcdy = cloneDetached(qcdAbcdyH, (pre + "_qcd_abcdy").c_str());
            if (!qcdSyst.empty() && qcdAbcdyUp && qcdAbcdyDn)
            {
                rr.syst["qcd_abcdy_" + qcdSyst + "Up"]   = cloneDetached(qcdAbcdyUp, (pre + "_qcd_abcdy_" + qcdSyst + "Up").c_str());
                rr.syst["qcd_abcdy_" + qcdSyst + "Down"] = cloneDetached(qcdAbcdyDn, (pre + "_qcd_abcdy_" + qcdSyst + "Down").c_str());
            }
        }

        if (accumulate)
        {
            auto acc = [](TH1D *&dst, const char *dnm, std::initializer_list<TH1D *> hs)
            {
                for (TH1D *h : hs)
                {
                    if (!h) continue;
                    if (!dst)
                    {
                        dst = (TH1D *)h->Clone(dnm);
                        dst->Reset();
                        dst->SetDirectory(nullptr);
                    }
                    dst->Add(h);
                }
            };
            const char *t = varTag[var];
            acc(h_leppt_inclusive[var],           Form("h_leppt_inclusive%s", t),           {hD});
            acc(h_leppt_inclusive_MC_signal[var], Form("h_leppt_inclusive%s_MC_signal", t), {hWp, hWm});
            acc(h_leppt_inclusive_MC_Z[var],      Form("h_leppt_inclusive%s_MC_Z", t),      {hZ});
            acc(h_leppt_inclusive_MC_Ztau[var],   Form("h_leppt_inclusive%s_MC_Ztau", t),   {hZtau});
            acc(h_leppt_inclusive_MC_Wtau[var],   Form("h_leppt_inclusive%s_MC_Wtau", t),   {hWptau, hWmtau});
            acc(h_leppt_inclusive_MC_TT[var],     Form("h_leppt_inclusive%s_MC_TT", t),     {hTT});
        }
        return true;
    };

    // Loop over rapidity bins and charge
    for (int iy = 0; iy < NY; ++iy)
    {
        // ---- MET ----
        TH1D *h_met_Wp = (TH1D *)f->Get(Form("h_met_Wp_y%d", iy));
        TH1D *h_met_Wm = (TH1D *)f->Get(Form("h_met_Wm_y%d", iy));

        TH1D *h_met_Wp_FB = (TH1D *)f->Get(Form("h_met_Wp_y%d_FB", iy));
        TH1D *h_met_Wm_FB = (TH1D *)f->Get(Form("h_met_Wm_y%d_FB", iy));

        // ---- MET MC Signal ----

        TH1D *h_met_Wp_MC_Wp = (TH1D *)f_Wp->Get(Form("h_met_Wp_y%d", iy));
        TH1D *h_met_Wm_MC_Wp = (TH1D *)f_Wp->Get(Form("h_met_Wm_y%d", iy));

        TH1D *h_met_Wp_FB_MC_Wp = (TH1D *)f_Wp->Get(Form("h_met_Wp_y%d_FB", iy));
        TH1D *h_met_Wm_FB_MC_Wp = (TH1D *)f_Wp->Get(Form("h_met_Wm_y%d_FB", iy));

        TH1D *h_met_Wp_MC_Wm = (TH1D *)f_Wm->Get(Form("h_met_Wp_y%d", iy));
        TH1D *h_met_Wm_MC_Wm = (TH1D *)f_Wm->Get(Form("h_met_Wm_y%d", iy));

        TH1D *h_met_Wp_FB_MC_Wm = (TH1D *)f_Wm->Get(Form("h_met_Wp_y%d_FB", iy));
        TH1D *h_met_Wm_FB_MC_Wm = (TH1D *)f_Wm->Get(Form("h_met_Wm_y%d_FB", iy));

        // ---- MET MC Z BK ----

        TH1D *h_met_Wp_MC_Z = (TH1D *)f_DY->Get(Form("h_met_Wp_y%d", iy));
        TH1D *h_met_Wm_MC_Z = (TH1D *)f_DY->Get(Form("h_met_Wm_y%d", iy));

        TH1D *h_met_Wp_FB_MC_Z = (TH1D *)f_DY->Get(Form("h_met_Wp_y%d_FB", iy));
        TH1D *h_met_Wm_FB_MC_Z = (TH1D *)f_DY->Get(Form("h_met_Wm_y%d_FB", iy));

        // ---- MET MC Ztau BK ----

        TH1D *h_met_Wp_MC_Ztau = (TH1D *)f_DYtau->Get(Form("h_met_Wp_y%d", iy));
        TH1D *h_met_Wm_MC_Ztau = (TH1D *)f_DYtau->Get(Form("h_met_Wm_y%d", iy));

        TH1D *h_met_Wp_FB_MC_Ztau = (TH1D *)f_DYtau->Get(Form("h_met_Wp_y%d_FB", iy));
        TH1D *h_met_Wm_FB_MC_Ztau = (TH1D *)f_DYtau->Get(Form("h_met_Wm_y%d_FB", iy));

        // ---- MET MC Wtau BK ----

        TH1D *h_met_Wp_MC_Wptau = (TH1D *)f_Wptau->Get(Form("h_met_Wp_y%d", iy));
        TH1D *h_met_Wm_MC_Wptau = (TH1D *)f_Wptau->Get(Form("h_met_Wm_y%d", iy));

        TH1D *h_met_Wp_FB_MC_Wptau = (TH1D *)f_Wptau->Get(Form("h_met_Wp_y%d_FB", iy));
        TH1D *h_met_Wm_FB_MC_Wptau = (TH1D *)f_Wptau->Get(Form("h_met_Wm_y%d_FB", iy));

        TH1D *h_met_Wp_MC_Wmtau = (TH1D *)f_Wmtau->Get(Form("h_met_Wp_y%d", iy));
        TH1D *h_met_Wm_MC_Wmtau = (TH1D *)f_Wmtau->Get(Form("h_met_Wm_y%d", iy));

        TH1D *h_met_Wp_FB_MC_Wmtau = (TH1D *)f_Wmtau->Get(Form("h_met_Wp_y%d_FB", iy));
        TH1D *h_met_Wm_FB_MC_Wmtau = (TH1D *)f_Wmtau->Get(Form("h_met_Wm_y%d_FB", iy));

        // ---- MET MC ttbar BK ----

        TH1D *h_met_Wp_MC_TT = (TH1D *)f_TT->Get(Form("h_met_Wp_y%d", iy));
        TH1D *h_met_Wm_MC_TT = (TH1D *)f_TT->Get(Form("h_met_Wm_y%d", iy));

        TH1D *h_met_Wp_FB_MC_TT = (TH1D *)f_TT->Get(Form("h_met_Wp_y%d_FB", iy));
        TH1D *h_met_Wm_FB_MC_TT = (TH1D *)f_TT->Get(Form("h_met_Wm_y%d_FB", iy));

        // ---- mT ----
        TH1D *h_mt_Wp = (TH1D *)f->Get(Form("h_mt_Wp_y%d", iy));
        TH1D *h_mt_Wm = (TH1D *)f->Get(Form("h_mt_Wm_y%d", iy));

        TH1D *h_mt_Wp_FB = (TH1D *)f->Get(Form("h_mt_Wp_y%d_FB", iy));
        TH1D *h_mt_Wm_FB = (TH1D *)f->Get(Form("h_mt_Wm_y%d_FB", iy));

        // ---- mT MC Signal ----

        TH1D *h_mt_Wp_MC_Wp = (TH1D *)f_Wp->Get(Form("h_mt_Wp_y%d", iy));
        TH1D *h_mt_Wm_MC_Wp = (TH1D *)f_Wp->Get(Form("h_mt_Wm_y%d", iy));

        TH1D *h_mt_Wp_FB_MC_Wp = (TH1D *)f_Wp->Get(Form("h_mt_Wp_y%d_FB", iy));
        TH1D *h_mt_Wm_FB_MC_Wp = (TH1D *)f_Wp->Get(Form("h_mt_Wm_y%d_FB", iy));

        TH1D *h_mt_Wp_MC_Wm = (TH1D *)f_Wm->Get(Form("h_mt_Wp_y%d", iy));
        TH1D *h_mt_Wm_MC_Wm = (TH1D *)f_Wm->Get(Form("h_mt_Wm_y%d", iy));

        TH1D *h_mt_Wp_FB_MC_Wm = (TH1D *)f_Wm->Get(Form("h_mt_Wp_y%d_FB", iy));
        TH1D *h_mt_Wm_FB_MC_Wm = (TH1D *)f_Wm->Get(Form("h_mt_Wm_y%d_FB", iy));

        // ---- mT MC Z BK ----

        TH1D *h_mt_Wp_MC_Z = (TH1D *)f_DY->Get(Form("h_mt_Wp_y%d", iy));
        TH1D *h_mt_Wm_MC_Z = (TH1D *)f_DY->Get(Form("h_mt_Wm_y%d", iy));

        TH1D *h_mt_Wp_FB_MC_Z = (TH1D *)f_DY->Get(Form("h_mt_Wp_y%d_FB", iy));
        TH1D *h_mt_Wm_FB_MC_Z = (TH1D *)f_DY->Get(Form("h_mt_Wm_y%d_FB", iy));

        // ---- mT MC Ztau BK ----

        TH1D *h_mt_Wp_MC_Ztau = (TH1D *)f_DYtau->Get(Form("h_mt_Wp_y%d", iy));
        TH1D *h_mt_Wm_MC_Ztau = (TH1D *)f_DYtau->Get(Form("h_mt_Wm_y%d", iy));

        TH1D *h_mt_Wp_FB_MC_Ztau = (TH1D *)f_DYtau->Get(Form("h_mt_Wp_y%d_FB", iy));
        TH1D *h_mt_Wm_FB_MC_Ztau = (TH1D *)f_DYtau->Get(Form("h_mt_Wm_y%d_FB", iy));

        // ---- mT MC Wtau BK ----

        TH1D *h_mt_Wp_MC_Wptau = (TH1D *)f_Wptau->Get(Form("h_mt_Wp_y%d", iy));
        TH1D *h_mt_Wm_MC_Wptau = (TH1D *)f_Wptau->Get(Form("h_mt_Wm_y%d", iy));

        TH1D *h_mt_Wp_FB_MC_Wptau = (TH1D *)f_Wptau->Get(Form("h_mt_Wp_y%d_FB", iy));
        TH1D *h_mt_Wm_FB_MC_Wptau = (TH1D *)f_Wptau->Get(Form("h_mt_Wm_y%d_FB", iy));

        TH1D *h_mt_Wp_MC_Wmtau = (TH1D *)f_Wmtau->Get(Form("h_mt_Wp_y%d", iy));
        TH1D *h_mt_Wm_MC_Wmtau = (TH1D *)f_Wmtau->Get(Form("h_mt_Wm_y%d", iy));

        TH1D *h_mt_Wp_FB_MC_Wmtau = (TH1D *)f_Wmtau->Get(Form("h_mt_Wp_y%d_FB", iy));
        TH1D *h_mt_Wm_FB_MC_Wmtau = (TH1D *)f_Wmtau->Get(Form("h_mt_Wm_y%d_FB", iy));

        // ---- mT MC ttbar BK ----

        TH1D *h_mt_Wp_MC_TT = (TH1D *)f_TT->Get(Form("h_mt_Wp_y%d", iy));
        TH1D *h_mt_Wm_MC_TT = (TH1D *)f_TT->Get(Form("h_mt_Wm_y%d", iy));

        TH1D *h_mt_Wp_FB_MC_TT = (TH1D *)f_TT->Get(Form("h_mt_Wp_y%d_FB", iy));
        TH1D *h_mt_Wm_FB_MC_TT = (TH1D *)f_TT->Get(Form("h_mt_Wm_y%d_FB", iy));

        if (!h_met_Wp || !h_met_Wm || !h_mt_Wp || !h_mt_Wm)
        {
            std::cerr << "[WARN] Missing histogram(s) at y bin " << iy << "\n";
            continue;
        }

        // Apply per-sample absolute normalization (sets the relative composition).
        // In-place Scale is safe: each MC histo is Get()'d fresh for this iy and
        // used once (the name carries iy), so there is no repeated-Get hazard.
        // Data histos (h_met_Wp/Wm, h_mt_Wp/Wm and their _FB) are NOT scaled.
        scaleAll(k_Wp,    {h_met_Wp_MC_Wp, h_met_Wm_MC_Wp, h_met_Wp_FB_MC_Wp, h_met_Wm_FB_MC_Wp,
                           h_mt_Wp_MC_Wp,  h_mt_Wm_MC_Wp,  h_mt_Wp_FB_MC_Wp,  h_mt_Wm_FB_MC_Wp});
        scaleAll(k_Wm,    {h_met_Wp_MC_Wm, h_met_Wm_MC_Wm, h_met_Wp_FB_MC_Wm, h_met_Wm_FB_MC_Wm,
                           h_mt_Wp_MC_Wm,  h_mt_Wm_MC_Wm,  h_mt_Wp_FB_MC_Wm,  h_mt_Wm_FB_MC_Wm});
        scaleAll(k_DY,    {h_met_Wp_MC_Z, h_met_Wm_MC_Z, h_met_Wp_FB_MC_Z, h_met_Wm_FB_MC_Z,
                           h_mt_Wp_MC_Z,  h_mt_Wm_MC_Z,  h_mt_Wp_FB_MC_Z,  h_mt_Wm_FB_MC_Z});
        scaleAll(k_DYtau, {h_met_Wp_MC_Ztau, h_met_Wm_MC_Ztau, h_met_Wp_FB_MC_Ztau, h_met_Wm_FB_MC_Ztau,
                           h_mt_Wp_MC_Ztau,  h_mt_Wm_MC_Ztau,  h_mt_Wp_FB_MC_Ztau,  h_mt_Wm_FB_MC_Ztau});
        scaleAll(k_Wptau, {h_met_Wp_MC_Wptau, h_met_Wm_MC_Wptau, h_met_Wp_FB_MC_Wptau, h_met_Wm_FB_MC_Wptau,
                           h_mt_Wp_MC_Wptau,  h_mt_Wm_MC_Wptau,  h_mt_Wp_FB_MC_Wptau,  h_mt_Wm_FB_MC_Wptau});
        scaleAll(k_Wmtau, {h_met_Wp_MC_Wmtau, h_met_Wm_MC_Wmtau, h_met_Wp_FB_MC_Wmtau, h_met_Wm_FB_MC_Wmtau,
                           h_mt_Wp_MC_Wmtau,  h_mt_Wm_MC_Wmtau,  h_mt_Wp_FB_MC_Wmtau,  h_mt_Wm_FB_MC_Wmtau});
        scaleAll(k_TT,    {h_met_Wp_MC_TT, h_met_Wm_MC_TT, h_met_Wp_FB_MC_TT, h_met_Wm_FB_MC_TT,
                           h_mt_Wp_MC_TT,  h_mt_Wm_MC_TT,  h_mt_Wp_FB_MC_TT,  h_mt_Wm_FB_MC_TT});

        // Per-y ABCD QCD (inclusive MET template shape x per-bin weight), per
        // charge and per binning. Null when the template is absent (e/missing).
        TH1D *qcd_met_Wp    = qcdPerY(qcdPlusBase,  wQcdWp[iy],   Form("qcd_met_Wp_y%d", iy));
        TH1D *qcd_met_Wm    = qcdPerY(qcdMinusBase, wQcdWm[iy],   Form("qcd_met_Wm_y%d", iy));
        TH1D *qcd_met_Wp_FB = qcdPerY(qcdPlusBase,  wQcdWpFB[iy], Form("qcd_met_Wp_y%d_FB", iy));
        TH1D *qcd_met_Wm_FB = qcdPerY(qcdMinusBase, wQcdWmFB[iy], Form("qcd_met_Wm_y%d_FB", iy));

        // Collect the absolute per-region Combine templates (MET discriminant).
        // signal carries BOTH W samples reconstructed with this charge; the EWK
        // backgrounds and the per-y ABCD QCD are projected the same way.
        // Each region also gets its LHE shape systematics (attachSysts: the
        // <hist>_<syst>Up/Down twins of the same per-sample histograms).
        regions.push_back(makeRegion(Form("Wp_lab_y%d", iy), "",
            h_met_Wp, h_met_Wp_MC_Wp, h_met_Wp_MC_Wm,
            h_met_Wp_MC_Z, h_met_Wp_MC_Ztau, h_met_Wp_MC_Wptau, h_met_Wp_MC_Wmtau, h_met_Wp_MC_TT,
            qcd_met_Wp));
        attachSysts(regions.back(), "", Form("h_met_Wp_y%d", iy));
        regions.push_back(makeRegion(Form("Wm_lab_y%d", iy), "",
            h_met_Wm, h_met_Wm_MC_Wp, h_met_Wm_MC_Wm,
            h_met_Wm_MC_Z, h_met_Wm_MC_Ztau, h_met_Wm_MC_Wptau, h_met_Wm_MC_Wmtau, h_met_Wm_MC_TT,
            qcd_met_Wm));
        attachSysts(regions.back(), "", Form("h_met_Wm_y%d", iy));
        regions.push_back(makeRegion(Form("Wp_fb_y%d", iy), "",
            h_met_Wp_FB, h_met_Wp_FB_MC_Wp, h_met_Wp_FB_MC_Wm,
            h_met_Wp_FB_MC_Z, h_met_Wp_FB_MC_Ztau, h_met_Wp_FB_MC_Wptau, h_met_Wp_FB_MC_Wmtau, h_met_Wp_FB_MC_TT,
            qcd_met_Wp_FB));
        attachSysts(regions.back(), "", Form("h_met_Wp_y%d_FB", iy));
        regions.push_back(makeRegion(Form("Wm_fb_y%d", iy), "",
            h_met_Wm_FB, h_met_Wm_FB_MC_Wp, h_met_Wm_FB_MC_Wm,
            h_met_Wm_FB_MC_Z, h_met_Wm_FB_MC_Ztau, h_met_Wm_FB_MC_Wptau, h_met_Wm_FB_MC_Wmtau, h_met_Wm_FB_MC_TT,
            qcd_met_Wm_FB));
        attachSysts(regions.back(), "", Form("h_met_Wm_y%d_FB", iy));

        if (!h_mt_inclusive)
        {
            h_mt_inclusive = (TH1D *)h_mt_Wp->Clone("h_mt_inclusive");
            h_mt_inclusive->Reset();
            h_mt_inclusive->SetDirectory(nullptr);

            h_mt_inclusive_MC_signal = (TH1D *)h_mt_Wp->Clone("h_mt_inclusive_MC_signal");
            h_mt_inclusive_MC_signal->Reset();
            h_mt_inclusive_MC_signal->SetDirectory(nullptr);

            h_mt_inclusive_MC_Z = (TH1D *)h_mt_Wp->Clone("h_mt_inclusive_MC_Z");
            h_mt_inclusive_MC_Z->Reset();
            h_mt_inclusive_MC_Z->SetDirectory(nullptr);

            h_mt_inclusive_MC_Ztau = (TH1D *)h_mt_Wp->Clone("h_mt_inclusive_MC_Ztau");
            h_mt_inclusive_MC_Ztau->Reset();
            h_mt_inclusive_MC_Ztau->SetDirectory(nullptr);

            h_mt_inclusive_MC_Wtau = (TH1D *)h_mt_Wp->Clone("h_mt_inclusive_MC_Wtau");
            h_mt_inclusive_MC_Wtau->Reset();
            h_mt_inclusive_MC_Wtau->SetDirectory(nullptr);

            h_mt_inclusive_MC_TT = (TH1D *)h_mt_Wp->Clone("h_mt_inclusive_MC_TT");
            h_mt_inclusive_MC_TT->Reset();
            h_mt_inclusive_MC_TT->SetDirectory(nullptr);
        }
        if (!h_met_inclusive)
        {
            h_met_inclusive = (TH1D *)h_met_Wp->Clone("h_met_inclusive");
            h_met_inclusive->Reset();
            h_met_inclusive->SetDirectory(nullptr); // avoid ROOT ownership issues

            h_met_inclusive_MC_signal = (TH1D *)h_met_Wp->Clone("h_met_inclusive_MC_signal");
            h_met_inclusive_MC_signal->Reset();
            h_met_inclusive_MC_signal->SetDirectory(nullptr); // avoid ROOT ownership issues

            h_met_inclusive_MC_Z = (TH1D *)h_met_Wp->Clone("h_met_inclusive_MC_Z");
            h_met_inclusive_MC_Z->Reset();
            h_met_inclusive_MC_Z->SetDirectory(nullptr); // avoid ROOT ownership issues

            h_met_inclusive_MC_Ztau = (TH1D *)h_met_Wp->Clone("h_met_inclusive_MC_Ztau");
            h_met_inclusive_MC_Ztau->Reset();
            h_met_inclusive_MC_Ztau->SetDirectory(nullptr);

            h_met_inclusive_MC_Wtau = (TH1D *)h_met_Wp->Clone("h_met_inclusive_MC_Wtau");
            h_met_inclusive_MC_Wtau->Reset();
            h_met_inclusive_MC_Wtau->SetDirectory(nullptr);

            h_met_inclusive_MC_TT = (TH1D *)h_met_Wp->Clone("h_met_inclusive_MC_TT");
            h_met_inclusive_MC_TT->Reset();
            h_met_inclusive_MC_TT->SetDirectory(nullptr);
        }

        h_met_inclusive->Add(h_met_Wp);
        h_met_inclusive->Add(h_met_Wm);

        h_met_inclusive_MC_signal->Add(h_met_Wp_MC_Wp);
        h_met_inclusive_MC_signal->Add(h_met_Wm_MC_Wp);
        h_met_inclusive_MC_signal->Add(h_met_Wp_MC_Wm);
        h_met_inclusive_MC_signal->Add(h_met_Wm_MC_Wm);

        h_met_inclusive_MC_Z->Add(h_met_Wp_MC_Z);
        h_met_inclusive_MC_Z->Add(h_met_Wm_MC_Z);

        h_met_inclusive_MC_Ztau->Add(h_met_Wp_MC_Ztau);
        h_met_inclusive_MC_Ztau->Add(h_met_Wm_MC_Ztau);

        h_met_inclusive_MC_Wtau->Add(h_met_Wp_MC_Wptau);
        h_met_inclusive_MC_Wtau->Add(h_met_Wm_MC_Wptau);
        h_met_inclusive_MC_Wtau->Add(h_met_Wp_MC_Wmtau);
        h_met_inclusive_MC_Wtau->Add(h_met_Wm_MC_Wmtau);

        h_met_inclusive_MC_TT->Add(h_met_Wp_MC_TT);
        h_met_inclusive_MC_TT->Add(h_met_Wm_MC_TT);

        h_mt_inclusive->Add(h_mt_Wp);
        h_mt_inclusive->Add(h_mt_Wm);

        h_mt_inclusive_MC_signal->Add(h_mt_Wp_MC_Wp);
        h_mt_inclusive_MC_signal->Add(h_mt_Wm_MC_Wp);
        h_mt_inclusive_MC_signal->Add(h_mt_Wp_MC_Wm);
        h_mt_inclusive_MC_signal->Add(h_mt_Wm_MC_Wm);

        h_mt_inclusive_MC_Z->Add(h_mt_Wp_MC_Z);
        h_mt_inclusive_MC_Z->Add(h_mt_Wm_MC_Z);

        h_mt_inclusive_MC_Ztau->Add(h_mt_Wp_MC_Ztau);
        h_mt_inclusive_MC_Ztau->Add(h_mt_Wm_MC_Ztau);

        h_mt_inclusive_MC_Wtau->Add(h_mt_Wp_MC_Wptau);
        h_mt_inclusive_MC_Wtau->Add(h_mt_Wm_MC_Wptau);
        h_mt_inclusive_MC_Wtau->Add(h_mt_Wp_MC_Wmtau);
        h_mt_inclusive_MC_Wtau->Add(h_mt_Wm_MC_Wmtau);

        h_mt_inclusive_MC_TT->Add(h_mt_Wp_MC_TT);
        h_mt_inclusive_MC_TT->Add(h_mt_Wm_MC_TT);

        if (h_met_Wp->Integral(1, h_met_Wp->GetNbinsX()) != h_mt_Wp->Integral(1, h_mt_Wp->GetNbinsX()))
        {
            cout << "ERROR check plot " << iy << endl;
            cout << h_met_Wp->Integral(1, h_met_Wp->GetNbinsX()) << " And " << h_mt_Wp->Integral(1, h_mt_Wp->GetNbinsX());
        }

        // If you want to close the file later safely, detach histograms:
        // (Not strictly required if you keep file open until end)
        // h_met_Wp->SetDirectory(0); h_met_Wm->SetDirectory(0);
        // h_mt_Wp->SetDirectory(0);  h_mt_Wm->SetDirectory(0);

        // --------- MET plots ----------
        {
            std::vector<std::string> box = {
                Form("Passing Events: %.0f", h_met_Wp->Integral(1, h_met_Wp->GetNbinsX())),
                Form("W signal MC: %.0f", sigInt({h_met_Wp_MC_Wp, h_met_Wp_MC_Wm}))};

            std::vector<TH1 *> bkgs = {
                h_met_Wp_MC_Wp,
                h_met_Wp_MC_Wm,
                h_met_Wp_MC_Z,
                h_met_Wp_MC_Ztau,
                h_met_Wp_MC_Wptau,
                h_met_Wp_MC_Wmtau};

            std::vector<std::string> names = {
                "W+",
                "W-",
                "DY",
                "DY Tau",
                "Wp Tau",
                "Wm Tau"};

            if (qcd_met_Wp) { bkgs.push_back(qcd_met_Wp); names.push_back("QCD"); }
            addTT(bkgs, names, h_met_Wp_MC_TT);
            box.push_back(totalLine(bkgs, names));

            SaveNicePlot1D_WithBkg(
                h_met_Wp,
                bkgs,
                names,
                outMetDir + Form("/met_Wp_y%d", iy), // outPathNoExt
                "PF MET (GeV)",                      // xTitle
                "Events / 2.0 GeV",                  // yTitle
                "",                                  // mainTitle (whatever you like)
                Channeltypewplus,                    // subTitle1
                yLabel[iy],                          // subTitle2
                box,                                 // info box lines
                psCounts,
                commonTuner);
        }

        {
            std::vector<std::string> box = {
                Form("Passing Events: %.0f", h_met_Wm->Integral(1, h_met_Wm->GetNbinsX())),
                Form("W signal MC: %.0f", sigInt({h_met_Wm_MC_Wp, h_met_Wm_MC_Wm}))};

            std::vector<TH1 *> bkgs = {
                h_met_Wm_MC_Wp,
                h_met_Wm_MC_Wm,
                h_met_Wm_MC_Z,
                h_met_Wm_MC_Ztau,
                h_met_Wm_MC_Wptau,
                h_met_Wm_MC_Wmtau};

            std::vector<std::string> names = {
                "W+",
                "W-",
                "DY",
                "DY Tau",
                "Wp Tau",
                "Wm Tau"};

            if (qcd_met_Wm) { bkgs.push_back(qcd_met_Wm); names.push_back("QCD"); }
            addTT(bkgs, names, h_met_Wm_MC_TT);
            box.push_back(totalLine(bkgs, names));

            SaveNicePlot1D_WithBkg(
                h_met_Wm,
                bkgs,
                names,
                outMetDir + Form("/met_Wm_y%d", iy),
                "PF MET (GeV)",
                "Events / 2.0 GeV",
                "",
                Channeltypewminus,
                yLabel[iy],
                box,
                psCounts,
                commonTuner);
        }

        // --------- Met with FB ratio --------

        {
            std::vector<std::string> box = {
                Form("Passing Events: %.0f", h_met_Wp_FB->Integral(1, h_met_Wp_FB->GetNbinsX())),
                Form("W signal MC: %.0f", sigInt({h_met_Wp_FB_MC_Wp, h_met_Wp_FB_MC_Wm}))};

            std::vector<TH1 *> bkgs = {
                h_met_Wp_FB_MC_Wp,
                h_met_Wp_FB_MC_Wm,
                h_met_Wp_FB_MC_Z,
                h_met_Wp_FB_MC_Ztau,
                h_met_Wp_FB_MC_Wptau,
                h_met_Wp_FB_MC_Wmtau};

            std::vector<std::string> names = {
                "W+",
                "W-",
                "DY",
                "DY Tau",
                "Wp Tau",
                "Wm Tau"};

            if (qcd_met_Wp_FB) { bkgs.push_back(qcd_met_Wp_FB); names.push_back("QCD"); }
            addTT(bkgs, names, h_met_Wp_FB_MC_TT);
            box.push_back(totalLine(bkgs, names));

            SaveNicePlot1D_WithBkg(
                h_met_Wp_FB,
                bkgs,
                names,
                outMetDir + Form("/met_Wp_y%d_FB", iy), // outPathNoExt
                "PF MET (GeV)",                         // xTitle
                "Events / 2.0 GeV",                     // yTitle
                "",                                     // mainTitle (whatever you like)
                Channeltypewplus,                       // subTitle1
                yLabel_FB[iy],                          // subTitle2
                box,                                    // info box lines
                psCounts,
                commonTuner);
        }

        {
            std::vector<std::string> box = {
                Form("Passing Events: %.0f", h_met_Wm_FB->Integral(1, h_met_Wm_FB->GetNbinsX())),
                Form("W signal MC: %.0f", sigInt({h_met_Wm_FB_MC_Wp, h_met_Wm_FB_MC_Wm}))};

            std::vector<TH1 *> bkgs = {
                h_met_Wm_FB_MC_Wp,
                h_met_Wm_FB_MC_Wm,
                h_met_Wm_FB_MC_Z,
                h_met_Wm_FB_MC_Ztau,
                h_met_Wm_FB_MC_Wptau,
                h_met_Wm_FB_MC_Wmtau};

            std::vector<std::string> names = {
                "W+",
                "W-",
                "DY",
                "DY Tau",
                "Wp Tau",
                "Wm Tau"};

            if (qcd_met_Wm_FB) { bkgs.push_back(qcd_met_Wm_FB); names.push_back("QCD"); }
            addTT(bkgs, names, h_met_Wm_FB_MC_TT);
            box.push_back(totalLine(bkgs, names));

            SaveNicePlot1D_WithBkg(
                h_met_Wm_FB,
                bkgs,
                names,
                outMetDir + Form("/met_Wm_y%d_FB", iy),
                "PF MET (GeV)",
                "Events / 2.0 GeV",
                "",
                Channeltypewminus,
                yLabel_FB[iy],
                box,
                psCounts,
                commonTuner);
        }

        // --------- mT plots ----------
        {
            std::vector<std::string> box = {Form("Passing Events: %.0f", h_mt_Wp->Integral(1, h_mt_Wp->GetNbinsX()))};

            std::vector<TH1 *> bkgs = {
                h_mt_Wp_MC_Wp,
                h_mt_Wp_MC_Wm,
                h_mt_Wp_MC_Z,
                h_mt_Wp_MC_Ztau,
                h_mt_Wp_MC_Wptau,
                h_mt_Wp_MC_Wmtau};

            std::vector<std::string> names = {
                "W+",
                "W-",
                "DY",
                "DY Tau",
                "Wp Tau",
                "Wm Tau"};
            addTT(bkgs, names, h_mt_Wp_MC_TT);

            SaveNicePlot1D_WithBkg(
                h_mt_Wp,
                bkgs,
                names,
                outMtDir + Form("/mt_Wp_y%d", iy),
                "m_{T} (GeV)",
                "Events / 2.5 GeV",
                "",
                Channeltypewplus,
                yLabel[iy],
                box,
                ps,
                commonTuner);
        }

        {
            std::vector<std::string> box = {Form("Passing Events: %.0f", h_mt_Wm->Integral(1, h_mt_Wm->GetNbinsX()))};

            std::vector<TH1 *> bkgs = {
                h_mt_Wm_MC_Wp,
                h_mt_Wm_MC_Wm,
                h_mt_Wm_MC_Z,
                h_mt_Wm_MC_Ztau,
                h_mt_Wm_MC_Wptau,
                h_mt_Wm_MC_Wmtau};

            std::vector<std::string> names = {
                "W+",
                "W-",
                "DY",
                "DY Tau",
                "Wp Tau",
                "Wm Tau"};
            addTT(bkgs, names, h_mt_Wm_MC_TT);

            SaveNicePlot1D_WithBkg(
                h_mt_Wm,
                bkgs,
                names,
                outMtDir + Form("/mt_Wm_y%d", iy),
                "m_{T} (GeV)",
                "Events / 2.5 GeV",
                "",
                Channeltypewminus,
                yLabel[iy],
                box,
                ps,
                commonTuner);
        }

        // --------- mT plots For FB ----------
        {
            std::vector<std::string> box = {Form("Passing Events: %.0f", h_mt_Wp_FB->Integral(1, h_mt_Wp_FB->GetNbinsX()))};

            std::vector<TH1 *> bkgs = {
                h_mt_Wp_FB_MC_Wp,
                h_mt_Wp_FB_MC_Wm,
                h_mt_Wp_FB_MC_Z,
                h_mt_Wp_FB_MC_Ztau,
                h_mt_Wp_FB_MC_Wptau,
                h_mt_Wp_FB_MC_Wmtau};

            std::vector<std::string> names = {
                "W+",
                "W-",
                "DY",
                "DY Tau",
                "Wp Tau",
                "Wm Tau"};
            addTT(bkgs, names, h_mt_Wp_FB_MC_TT);

            SaveNicePlot1D_WithBkg(
                h_mt_Wp_FB,
                bkgs,
                names,
                outMtDir + Form("/mt_Wp_y%d_FB", iy),
                "m_{T} (GeV)",
                "Events / 2.5 GeV",
                "",
                Channeltypewplus,
                yLabel_FB[iy],
                box,
                ps,
                commonTuner);
        }

        {
            std::vector<std::string> box = {Form("Passing Events: %.0f", h_mt_Wm_FB->Integral(1, h_mt_Wm_FB->GetNbinsX()))};

            std::vector<TH1 *> bkgs = {
                h_mt_Wm_FB_MC_Wp,
                h_mt_Wm_FB_MC_Wm,
                h_mt_Wm_FB_MC_Z,
                h_mt_Wm_FB_MC_Ztau,
                h_mt_Wm_FB_MC_Wptau,
                h_mt_Wm_FB_MC_Wmtau};

            std::vector<std::string> names = {
                "W+",
                "W-",
                "DY",
                "DY Tau",
                "Wp Tau",
                "Wm Tau"};
            addTT(bkgs, names, h_mt_Wm_FB_MC_TT);

            SaveNicePlot1D_WithBkg(
                h_mt_Wm_FB,
                bkgs,
                names,
                outMtDir + Form("/mt_Wm_y%d_FB", iy),
                "m_{T} (GeV)",
                "Events / 2.5 GeV",
                "",
                Channeltypewminus,
                yLabel_FB[iy],
                box,
                ps,
                commonTuner);
        }

        // --------- m_T > 40 plots (+ the leppt_mt40 / met_mt40 Combine inputs) ---
        // One pass per selection variant (kVarMt40 = lepton pT, kVarMetMt40 =
        // PF MET of the same events): the display QCD, the in-fit bases and the
        // abcdy factors are the variant's own; the per-y weights are shared.
        for (int var = 0; var < kNVar; ++var)
        {
            const std::string vs = (var == kVarMt40) ? "" : "_met"; // unique clone names per variant
            TH1D *dP = dispBaseV[var][0], *dM = dispBaseV[var][1];
            TH1D *aP = abcdBaseV[var][0], *aM = abcdBaseV[var][1];
            TH1D *aPU = abcdUpBaseV[var][0], *aMU = abcdUpBaseV[var][1];
            TH1D *aPD = abcdDnBaseV[var][0], *aMD = abcdDnBaseV[var][1];
            // Per-y QCD = inclusive template shape x the SAME per-bin low-MET
            // excess weights as the MET stacks (the QCD y-distribution does not
            // depend on which variable is plotted). The m_T cut does not change
            // where QCD sits in y, so the m_T>40 split reuses those weights.
            TH1D *qcd_pt40_Wp    = qcdPerY(dP, wQcdWp[iy],   Form("qcd%s_pt_mt40_Wp_y%d", vs.c_str(), iy));
            TH1D *qcd_pt40_Wm    = qcdPerY(dM, wQcdWm[iy],   Form("qcd%s_pt_mt40_Wm_y%d", vs.c_str(), iy));
            TH1D *qcd_pt40_Wp_FB = qcdPerY(dP, wQcdWpFB[iy], Form("qcd%s_pt_mt40_Wp_y%d_FB", vs.c_str(), iy));
            TH1D *qcd_pt40_Wm_FB = qcdPerY(dM, wQcdWmFB[iy], Form("qcd%s_pt_mt40_Wm_y%d_FB", vs.c_str(), iy));

            // in-fit ABCD (abcd-mode) per-y templates: same shape, same per-y
            // weights (Sum_y = A0 per charge+binning since Sum w = 1), total
            // renormalized to A0. Written as the 7th SR template `qcd_abcd`;
            // deliberately NOT drawn in the display stacks.
            TH1D *qcdA_Wp    = qcdPerY(aP, wQcdWp[iy],   Form("qcd%s_abcd_Wp_y%d", vs.c_str(), iy));
            TH1D *qcdA_Wm    = qcdPerY(aM, wQcdWm[iy],   Form("qcd%s_abcd_Wm_y%d", vs.c_str(), iy));
            TH1D *qcdA_Wp_FB = qcdPerY(aP, wQcdWpFB[iy], Form("qcd%s_abcd_Wp_y%d_FB", vs.c_str(), iy));
            TH1D *qcdA_Wm_FB = qcdPerY(aM, wQcdWmFB[iy], Form("qcd%s_abcd_Wm_y%d_FB", vs.c_str(), iy));
            // electrons: the qcdFF slope twins, same per-y weights (null otherwise)
            TH1D *qA_Wp_U    = qcdPerY(aPU, wQcdWp[iy],   Form("qcd%s_abcd_Wp_y%d_ffU", vs.c_str(), iy));
            TH1D *qA_Wp_D    = qcdPerY(aPD, wQcdWp[iy],   Form("qcd%s_abcd_Wp_y%d_ffD", vs.c_str(), iy));
            TH1D *qA_Wm_U    = qcdPerY(aMU, wQcdWm[iy],   Form("qcd%s_abcd_Wm_y%d_ffU", vs.c_str(), iy));
            TH1D *qA_Wm_D    = qcdPerY(aMD, wQcdWm[iy],   Form("qcd%s_abcd_Wm_y%d_ffD", vs.c_str(), iy));
            TH1D *qA_Wp_FB_U = qcdPerY(aPU, wQcdWpFB[iy], Form("qcd%s_abcd_Wp_y%d_FB_ffU", vs.c_str(), iy));
            TH1D *qA_Wp_FB_D = qcdPerY(aPD, wQcdWpFB[iy], Form("qcd%s_abcd_Wp_y%d_FB_ffD", vs.c_str(), iy));
            TH1D *qA_Wm_FB_U = qcdPerY(aMU, wQcdWmFB[iy], Form("qcd%s_abcd_Wm_y%d_FB_ffU", vs.c_str(), iy));
            TH1D *qA_Wm_FB_D = qcdPerY(aMD, wQcdWmFB[iy], Form("qcd%s_abcd_Wm_y%d_FB_ffD", vs.c_str(), iy));
            // abcdy (2026-10-19): the same bases x the per-group factor x the
            // within-group low-MET weight (abcdyFac); null when !abcdyOkV[var].
            // [0..3] = Wp lab, Wm lab, Wp fb, Wm fb; nominal / slope Up / Down
            TH1D *qY[4] = {nullptr, nullptr, nullptr, nullptr};
            TH1D *qYU[4] = {nullptr, nullptr, nullptr, nullptr}, *qYD[4] = {nullptr, nullptr, nullptr, nullptr};
            if (abcdyOkV[var])
                for (int k = 0; k < 4; ++k)
                {
                    const int ic = k % 2, ib = k / 2;
                    const double fy = abcdyFac[ic][ib][iy];
                    const char *sfx = ib ? "_FB" : "";
                    const char *cgk = ic ? "Wm" : "Wp";
                    qY[k]  = qcdPerY(ic ? aM : aP, fy, Form("qcd%s_abcdy_%s_y%d%s", vs.c_str(), cgk, iy, sfx));
                    qYU[k] = qcdPerY(ic ? aMU : aPU, fy, Form("qcd%s_abcdy_%s_y%d%s_ffU", vs.c_str(), cgk, iy, sfx));
                    qYD[k] = qcdPerY(ic ? aMD : aPD, fy, Form("qcd%s_abcdy_%s_y%d%s_ffD", vs.c_str(), cgk, iy, sfx));
                }

            const std::string mtLab = " , m_{T} > 40 GeV";
            const std::string od = outBase + "/" + varDisc[var] + "/" + varDisc[var];
            bool ok40 = true;
            ok40 &= lepPtStack(var, iy, "Wp", "",    qcd_pt40_Wp,    Channeltypewplus,  yLabel[iy] + mtLab,
                               od + Form("_Wp_y%d", iy),    /*accumulate=*/true,  qcdA_Wp,
                               qA_Wp_U, qA_Wp_D, ffSystName[0], qY[0], qYU[0], qYD[0]);
            ok40 &= lepPtStack(var, iy, "Wm", "",    qcd_pt40_Wm,    Channeltypewminus, yLabel[iy] + mtLab,
                               od + Form("_Wm_y%d", iy),    /*accumulate=*/true,  qcdA_Wm,
                               qA_Wm_U, qA_Wm_D, ffSystName[1], qY[1], qYU[1], qYD[1]);
            ok40 &= lepPtStack(var, iy, "Wp", "_FB", qcd_pt40_Wp_FB, Channeltypewplus,  yLabel_FB[iy] + mtLab,
                               od + Form("_Wp_y%d_FB", iy), /*accumulate=*/false, qcdA_Wp_FB,
                               qA_Wp_FB_U, qA_Wp_FB_D, ffSystName[0], qY[2], qYU[2], qYD[2]);
            ok40 &= lepPtStack(var, iy, "Wm", "_FB", qcd_pt40_Wm_FB, Channeltypewminus, yLabel_FB[iy] + mtLab,
                               od + Form("_Wm_y%d_FB", iy), /*accumulate=*/false, qcdA_Wm_FB,
                               qA_Wm_FB_U, qA_Wm_FB_D, ffSystName[1], qY[3], qYU[3], qYD[3]);
            if (!ok40 && iy == 0)
                std::cerr << "[WARN] " << varStem[var] << "_* not found (skim predates them); "
                          << varDisc[var] << " stacks skipped.\n";
        }
    }

    {

        std::vector<std::string> box = {Form("Passing Events: %.0f", h_mt_inclusive->Integral(1, h_mt_inclusive->GetNbinsX()))};

        // Inclusive m_T: full data/MC stack (signal + DY + tau backgrounds) for
        // BOTH channels. (The muon channel previously used a signal-only,
        // peak-normalized comparison; reverted to all components for consistency
        // with the MET inclusive plot and the per-bin plots.)
        std::vector<TH1 *> bkgs = {
            h_mt_inclusive_MC_signal,
            h_mt_inclusive_MC_Z,
            h_mt_inclusive_MC_Ztau,
            h_mt_inclusive_MC_Wtau};

        std::vector<std::string> names = {
            "W+/W-",
            "DY",
            "DY tau",
            "W+/W- tau"};
        addTT(bkgs, names, h_mt_inclusive_MC_TT);

        SaveNicePlot1D_WithBkg(
            h_mt_inclusive,
            bkgs,
            names,
            outMtDir + Form("/h_mt_inclusive"),
            "m_{T} (GeV)",
            "Events / 2.5 GeV",
            "",
            Channeltype,
            "inclusive",
            box,
            ps,
            commonTuner);
    }

    {

        std::vector<std::string> box = {
            Form("Passing Events: %.0f", h_met_inclusive->Integral(1, h_met_inclusive->GetNbinsX())),
            Form("W signal MC: %.0f", sigInt({h_met_inclusive_MC_signal}))};

        std::vector<TH1 *> bkgs = {
            h_met_inclusive_MC_signal,
            h_met_inclusive_MC_Z,
            h_met_inclusive_MC_Ztau,
            h_met_inclusive_MC_Wtau};

        std::vector<std::string> names = {
            "W+/W-",
            "DY",
            "DY tau",
            "W+/W- tau"};

        if (qcd_met_incl) { bkgs.push_back(qcd_met_incl); names.push_back("QCD (ABCD)"); }
        addTT(bkgs, names, h_met_inclusive_MC_TT);
        box.push_back(totalLine(bkgs, names));

        SaveNicePlot1D_WithBkg(
            h_met_inclusive,
            bkgs,
            names,
            outMetDir + Form("/h_met_inclusive"),
            "PF MET (GeV)",
            "Events / 2.0 GeV",
            "",
            Channeltype,
            "inclusive",
            box,
            psCounts,
            commonTuner);
    }

    for (int var = 0; var < kNVar; ++var)
    {
        // Inclusive lepton-pT stack (display only): the absolute-normalization
        // closure test in the fit's own variable, with the discriminant's
        // m_T > 40 cut applied -- QCD is strongly suppressed there and the
        // Jacobian peak should dominate.
        if (!h_leppt_inclusive[var]) continue;

        std::vector<std::string> box = {
            Form("Passing Events: %.0f", h_leppt_inclusive[var]->Integral(1, h_leppt_inclusive[var]->GetNbinsX())),
            Form("W signal MC: %.0f", sigInt({h_leppt_inclusive_MC_signal[var]}))};

        std::vector<TH1 *> bkgs;
        std::vector<std::string> names;
        auto pushIf = [&](TH1 *h, const char *n) { if (h) { bkgs.push_back(h); names.push_back(n); } };
        pushIf(h_leppt_inclusive_MC_signal[var], "W+/W-");
        pushIf(h_leppt_inclusive_MC_Z[var], "DY");
        pushIf(h_leppt_inclusive_MC_Ztau[var], "DY tau");
        pushIf(h_leppt_inclusive_MC_Wtau[var], "W+/W- tau");
        pushIf(dispInclV[var], "QCD (ABCD)");
        addTT(bkgs, names, h_leppt_inclusive_MC_TT[var]);
        box.push_back(totalLine(bkgs, names));
        for (const std::string &l : prefitLines(h_leppt_inclusive[var],
                                                {h_leppt_inclusive_MC_signal[var], h_leppt_inclusive_MC_Wtau[var]}, bkgs))
            box.push_back(l);

        // pT axis starts at the selection floor (see lepPtStack); MET keeps its own
        PlotStyle psPt = psLepPt;
        if (varXHi[var] > varXLo[var])
        {
            psPt.xRangeLo = varXLo[var];
            psPt.xRangeHi = varXHi[var];
        }
        SaveNicePlot1D_WithBkg(
            h_leppt_inclusive[var],
            bkgs,
            names,
            outBase + "/" + varDisc[var] + "/" + varStem[var] + "_inclusive",
            varXTitle[var],
            "Events / 2.0 GeV",
            "",
            Channeltype,
            "inclusive, m_{T} > 40 GeV",
            box,
            psPt,
            commonTuner);
    }

    // --- in-fit ABCD control-region channels (QCD_MODE=abcd, leppt_mt40 only) --
    // Three counting CRs per charge on the (relIso x m_T) plane, from the
    // exported abcd_counts_*: CRB (iso-pass, m_T<30), CRC (anti-iso, m_T>40 --
    // the template's own source region), CRD (anti-iso, m_T<30). All 1-bin
    // templates: data_obs = raw data count; qcd = the prefit EWK-subtracted
    // count (the free CR scale's x1 anchor); z/ztau in CRB ride r_Z via the
    // fork's existing map regex; the CRB W-related content is written BOTH ways
    // -- per-y w_lab_y*/w_fb_y* (split by the SR per-y signal fractions and
    // mapped to the r POIs = the floating subtraction) AND a single frozen
    // wfix -- the card generator picks one (QCD_WCR=float|frozen). CRC/CRD
    // carry one frozen `ewk` (<= 2.5% of those regions). ttbar (2026-10-16):
    // its own frozen `ttbar` in CRB (from abcd_counts "ttB" -- like the SR
    // ttbar it must not ride r_Z, so it is NOT inside zB), and inside `ewk`
    // in CRC/CRD (zC40/zD hold every non-W sample, ttbar included).
    struct CRTemplates {
        std::string dir;
        TH1D *data = nullptr, *qcd = nullptr, *z = nullptr, *ztau = nullptr,
             *wfix = nullptr, *ewk = nullptr, *ttbar = nullptr;
        std::vector<TH1D *> wLab, wFb;
        // abcdy group CRB (2026-10-19): the per-y W parts of THIS group only,
        // written under their own names (w_<bin>_y<i>), not by vector position
        std::vector<std::pair<std::string, TH1D *>> wNamed;
    };
    // One set per selection variant (2026-10-20): the counts are the SAME for
    // leppt_mt40 and met_mt40 (same events, same m_T-plane regions); only the
    // CRB W split, which follows that variant's own SR signal fractions, differs.
    std::vector<CRTemplates> crRegionsV[kNVar];
    std::vector<CRTemplates> crRegionsYV[kNVar]; // abcdy: {Wp,Wm}_{lab,fb}_{fwd,bwd,cen}_CR{B,C,D}
    for (int var = 0; var < kNVar; ++var) // (the body is the single-variant block, unchanged)
    if (abcdBaseV[var][0] && abcdBaseV[var][1] && !regionsPt[var].empty())
    {
        std::vector<CRTemplates> &crRegions  = crRegionsV[var];
        std::vector<CRTemplates> &crRegionsY = crRegionsYV[var];
        const std::string crFile = std::string("combine_input_W_") + varDisc[var] + ".root";
        auto mk1 = [](const std::string &nm, double v, double e) -> TH1D * {
            TH1D *h = new TH1D(nm.c_str(), nm.c_str(), 1, 0.0, 1.0);
            h->SetDirectory(nullptr);
            h->SetBinContent(1, v);
            h->SetBinError(1, e);
            return h;
        };
        // Per-y (signal+wtau) fractions of the SR, per charge x binning, from
        // the already-collected mt40 regions (the same integrals the fit sees).
        // Approximates the y composition of the CRB W content by the SR's --
        // second-order on a ~12% component, see the qcd_abcd.C channel report.
        auto srFracs = [&](const std::string &chg, const std::string &bin) -> std::vector<double> {
            std::vector<double> fr(NY, 0.0);
            double sum = 0.0;
            for (int iy = 0; iy < NY; ++iy)
            {
                const std::string want = Form("%s_%s_y%d", chg.c_str(), bin.c_str(), iy);
                for (const auto &r : regionsPt[var])
                {
                    if (r.dir != want) continue;
                    double v = 0.0;
                    if (r.sig)  v += r.sig->Integral();
                    if (r.wtau) v += r.wtau->Integral();
                    fr[iy] = std::max(0.0, v);
                    sum += fr[iy];
                    break;
                }
            }
            for (int iy = 0; iy < NY; ++iy)
                fr[iy] = (sum > 0.0) ? fr[iy] / sum : 1.0 / NY;
            return fr;
        };
        for (int ic = 0; ic < 2; ++ic)
        {
            const std::string cg = (ic == 0) ? "Wp" : "Wm";
            TH1D *cnt = (ic == 0) ? abcdCountsPlus : abcdCountsMinus;
            auto cb = [&](const char *lab, double *err) { return countBin(cnt, lab, err); };
            double e = 0.0, v = 0.0;

            CRTemplates crb;
            crb.dir = cg + "_CRB";
            v = cb("dataB", &e);
            crb.data = mk1(crb.dir + "_data_obs", v, e);
            v = cb("qcdB0", &e);
            if (v < 0.0) { std::cerr << "[WARN] " << crb.dir << ": negative qcdB0 clamped to 0\n"; v = 0.0; }
            crb.qcd = mk1(crb.dir + "_qcd", v, e);
            v = cb("zB", &e);
            crb.z = mk1(crb.dir + "_z", v, e);
            v = cb("ztauB", &e);
            crb.ztau = mk1(crb.dir + "_ztau", v, e);
            v = cb("ttB", &e);
            crb.ttbar = mk1(crb.dir + "_ttbar", v, e);
            double wBe = 0.0;
            const double wBv = cb("wB", &wBe);
            crb.wfix = mk1(crb.dir + "_wfix", wBv, wBe);
            const std::vector<double> fLab = srFracs(cg, "lab");
            const std::vector<double> fFb  = srFracs(cg, "fb");
            for (int iy = 0; iy < NY; ++iy)
            {
                crb.wLab.push_back(mk1(Form("%s_w_lab_y%d", crb.dir.c_str(), iy), wBv * fLab[iy], wBe * fLab[iy]));
                crb.wFb.push_back(mk1(Form("%s_w_fb_y%d", crb.dir.c_str(), iy), wBv * fFb[iy], wBe * fFb[iy]));
            }
            crRegions.push_back(crb);

            CRTemplates crc;
            crc.dir = cg + "_CRC";
            v = cb("dataC40", &e);
            crc.data = mk1(crc.dir + "_data_obs", v, e);
            v = cb("qcdC0", &e);
            if (v < 0.0) { std::cerr << "[WARN] " << crc.dir << ": negative qcdC0 clamped to 0\n"; v = 0.0; }
            crc.qcd = mk1(crc.dir + "_qcd", v, e);
            {
                double e1 = 0.0, e2 = 0.0;
                const double v1 = cb("wC40", &e1), v2 = cb("zC40", &e2);
                crc.ewk = mk1(crc.dir + "_ewk", v1 + v2, std::sqrt(e1 * e1 + e2 * e2));
            }
            crRegions.push_back(crc);

            CRTemplates crd;
            crd.dir = cg + "_CRD";
            v = cb("dataD", &e);
            crd.data = mk1(crd.dir + "_data_obs", v, e);
            v = cb("qcdD0", &e);
            if (v < 0.0) { std::cerr << "[WARN] " << crd.dir << ": negative qcdD0 clamped to 0\n"; v = 0.0; }
            crd.qcd = mk1(crd.dir + "_qcd", v, e);
            {
                double e1 = 0.0, e2 = 0.0;
                const double v1 = cb("wD", &e1), v2 = cb("zD", &e2);
                crd.ewk = mk1(crd.dir + "_ewk", v1 + v2, std::sqrt(e1 * e1 + e2 * e2));
            }
            crRegions.push_back(crd);
        }
        std::cout << "[QCD-ABCD] built " << crRegions.size()
                  << " CR channel dirs (CRB/CRC/CRD x Wp/Wm) for " << crFile << "\n";

        // --- abcdy (2026-10-19): the SAME three CRs per rapidity GROUP and
        // binning, {Wp,Wm}_{lab,fb}_{fwd,bwd,cen}_CR{B,C,D}, from the group
        // counts abcd_counts_*_<bin>_<g>. Contents mirror the inclusive ones:
        // CRB data_obs / qcd / z / ztau / ttbar / wfix + the group's per-y W
        // parts w_<bin>_y<i> (wB_g split by the SR signal+wtau fractions WITHIN
        // the group, same binning -> mapped to r_<C>_y<i> by the fork);
        // CRC / CRD data_obs / qcd / ewk. The FB groups use the FB-binned
        // planes, i.e. the FB window's own control regions.
        if (abcdyOkV[var])
        {
            for (int ic = 0; ic < 2; ++ic)
                for (int ib = 0; ib < 2; ++ib)
                {
                    const std::string cg = (ic == 0) ? "Wp" : "Wm";
                    const std::string bn = abcdyBin[ib];
                    const std::vector<double> fr = srFracs(cg, bn); // over all 12 bins
                    for (int g = 0; g < pOAbcdY::kNGroups; ++g)
                    {
                        TH1D *cnt = abcdyCounts[ic][ib][g];
                        auto cb = [&](const char *lab, double *err) { return countBin(cnt, lab, err); };
                        const std::string pre = cg + "_" + bn + "_" + pOAbcdY::kGroupName[g];
                        double e = 0.0, v = 0.0;

                        CRTemplates crb;
                        crb.dir = pre + "_CRB";
                        v = cb("dataB", &e);
                        crb.data = mk1(crb.dir + "_data_obs", v, e);
                        v = cb("qcdB0", &e);
                        crb.qcd = mk1(crb.dir + "_qcd", v, e);   // > 0 (abcdyOk)
                        v = cb("zB", &e);
                        crb.z = mk1(crb.dir + "_z", v, e);
                        v = cb("ztauB", &e);
                        crb.ztau = mk1(crb.dir + "_ztau", v, e);
                        v = cb("ttB", &e);
                        crb.ttbar = mk1(crb.dir + "_ttbar", v, e);
                        double wBe = 0.0;
                        const double wBv = cb("wB", &wBe);
                        crb.wfix = mk1(crb.dir + "_wfix", wBv, wBe);
                        const std::vector<int> gb = pOAbcdY::BinsOf(g);
                        double fsum = 0.0;
                        for (int iy : gb) fsum += fr[iy];
                        for (int iy : gb)
                        {
                            const double fi = (fsum > 0.0) ? fr[iy] / fsum : 1.0 / gb.size();
                            const std::string wn = Form("w_%s_y%d", bn.c_str(), iy);
                            crb.wNamed.push_back({wn, mk1(crb.dir + "_" + wn, wBv * fi, wBe * fi)});
                        }
                        crRegionsY.push_back(crb);

                        for (int rg = 0; rg < 2; ++rg) // CRC (anti-iso, m_T>40), CRD (anti-iso, m_T<30)
                        {
                            CRTemplates cr;
                            cr.dir = pre + (rg == 0 ? "_CRC" : "_CRD");
                            v = cb(rg == 0 ? "dataC40" : "dataD", &e);
                            cr.data = mk1(cr.dir + "_data_obs", v, e);
                            v = cb(rg == 0 ? "qcdC0" : "qcdD0", &e);
                            cr.qcd = mk1(cr.dir + "_qcd", v, e);
                            double e1 = 0.0, e2 = 0.0;
                            const double v1 = cb(rg == 0 ? "wC40" : "wD", &e1), v2 = cb(rg == 0 ? "zC40" : "zD", &e2);
                            cr.ewk = mk1(cr.dir + "_ewk", v1 + v2, std::sqrt(e1 * e1 + e2 * e2));
                            crRegionsY.push_back(cr);
                        }
                    }
                }
            std::cout << "[QCD-ABCDY] built " << crRegionsY.size()
                      << " group CR dirs ({Wp,Wm} x {lab,fb} x {fwd,bwd,cen} x CRB/CRC/CRD) for " << crFile << "\n";
        }
    }

    {
        // --- Structured per-region Combine input files -------------------------
        // One TDirectory per fit region; each holds the 7 absolute templates
        // (8 for leppt_mt40: + qcd_abcd) and, for leppt_mt40, the 6 CR dirs.
        // Written once per discriminant: PF MET (the nominal fit input) and the
        // two lepton-pT variants (SEPARATE files -- same region names with a
        // different discriminant must never share a file).
        // crsY (2026-10-19): the abcdy group CR dirs, written after the 6
        // inclusive ones (leppt_mt40 only; nullptr = none)
        // snames (2026-10-20): the systematics the sidecar may list for this
        // input (default: the channel's systNames; met_mt40 adds metRes)
        auto writeCombineInput = [&](const std::string &combineOut,
                                     std::vector<RegionTemplates> &regs,
                                     const std::string &tag,
                                     const std::vector<CRTemplates> *crs,
                                     const std::vector<CRTemplates> *crsY = nullptr,
                                     const std::vector<std::string> *snames = nullptr) {
            const std::vector<std::string> &sideNames = snames ? *snames : systNames;
            TFile *fcomb = TFile::Open(combineOut.c_str(), "RECREATE");
            if (!fcomb || fcomb->IsZombie())
            {
                std::cerr << "[ERROR] Cannot create output file: " << combineOut << "\n";
                return;
            }

            // Per-charge inclusive = sum over the lab-frame bins; W_incl = both charges.
            auto accumulate = [&](const std::string &prefix, const std::string &dir) -> RegionTemplates {
                RegionTemplates s;
                s.dir = dir;
                s.data = s.sig = s.z = s.ztau = s.wtau = s.qcd = nullptr;
                for (auto &r : regs) {
                    if (r.dir.compare(0, prefix.size(), prefix) != 0) continue;
                    auto add = [&](TH1D *&dst, TH1D *src, const char *suf) {
                        if (!src) return;
                        if (!dst) { dst = (TH1D *)src->Clone((dir + tag + suf).c_str()); dst->SetDirectory(nullptr); }
                        else dst->Add(src);
                    };
                    add(s.data, r.data, "_data_obs");
                    add(s.sig,  r.sig,  "_signal");
                    add(s.z,    r.z,    "_z");
                    add(s.ztau, r.ztau, "_ztau");
                    add(s.wtau, r.wtau, "_wtau");
                    add(s.ttbar, r.ttbar, "_ttbar");
                    add(s.qcd,  r.qcd,  "_qcd");
                    add(s.qcdAbcd, r.qcdAbcd, "_qcd_abcd");
                    add(s.qcdAbcdy, r.qcdAbcdy, "_qcd_abcdy");
                    for (auto &kv : r.syst) add(s.syst[kv.first], kv.second, ("_" + kv.first).c_str());
                }
                return s;
            };
            RegionTemplates Wp_incl = accumulate("Wp_lab_", "Wp_incl");
            RegionTemplates Wm_incl = accumulate("Wm_lab_", "Wm_incl");
            RegionTemplates W_incl;
            W_incl.dir  = "W_incl";
            W_incl.data = sum2(Wp_incl.data, Wm_incl.data, ("W_incl" + tag + "_data_obs").c_str());
            W_incl.sig  = sum2(Wp_incl.sig,  Wm_incl.sig,  ("W_incl" + tag + "_signal").c_str());
            W_incl.z    = sum2(Wp_incl.z,    Wm_incl.z,    ("W_incl" + tag + "_z").c_str());
            W_incl.ztau = sum2(Wp_incl.ztau, Wm_incl.ztau, ("W_incl" + tag + "_ztau").c_str());
            W_incl.wtau = sum2(Wp_incl.wtau, Wm_incl.wtau, ("W_incl" + tag + "_wtau").c_str());
            W_incl.ttbar = sum2(Wp_incl.ttbar, Wm_incl.ttbar, ("W_incl" + tag + "_ttbar").c_str());
            W_incl.qcd  = sum2(Wp_incl.qcd,  Wm_incl.qcd,  ("W_incl" + tag + "_qcd").c_str());
            W_incl.qcdAbcd = sum2(Wp_incl.qcdAbcd, Wm_incl.qcdAbcd, ("W_incl" + tag + "_qcd_abcd").c_str());
            W_incl.qcdAbcdy = sum2(Wp_incl.qcdAbcdy, Wm_incl.qcdAbcdy, ("W_incl" + tag + "_qcd_abcdy").c_str());
            for (auto &kv : Wp_incl.syst)
            {
                auto it = Wm_incl.syst.find(kv.first);
                W_incl.syst[kv.first] = sum2(kv.second, it != Wm_incl.syst.end() ? it->second : nullptr,
                                             ("W_incl" + tag + "_" + kv.first).c_str());
            }

            auto writeRegionTemplates = [&](const RegionTemplates &r) {
                TDirectory *d = fcomb->mkdir(r.dir.c_str());
                if (!d) { std::cerr << "[WARN] mkdir failed: " << r.dir << "\n"; return; }
                auto wn = [&](TH1D *h, const char *nm) {
                    if (!h) { std::cerr << "[WARN] " << r.dir << ": missing " << nm << "\n"; return; }
                    d->cd();
                    TH1D *hc = (TH1D *)h->Clone(nm);
                    hc->SetDirectory(d);
                    // text2workspace cannot build a pdf from an ALL-ZERO shape; floor
                    // empty MC templates (possible in sparse tail bins) to a
                    // negligible epsilon so the region's card stays fittable.
                    if (strcmp(nm, "data_obs") != 0 && hc->Integral() <= 0.0) {
                        std::cerr << "[WARN] " << combineOut << " " << r.dir << "/" << nm
                                  << " is empty -> flooring central bin to 1e-6 for Combine\n";
                        hc->SetBinContent(hc->GetNbinsX() / 2, 1e-6);
                    }
                    hc->Write(nm, TObject::kOverwrite);
                };
                wn(r.data, "data_obs");
                wn(r.sig,  "signal");
                wn(r.z,    "z");
                wn(r.ztau, "ztau");
                wn(r.wtau, "wtau");
                wn(r.ttbar, "ttbar");
                if (r.qcd) wn(r.qcd, "qcd");
                if (r.qcdAbcd) wn(r.qcdAbcd, "qcd_abcd");
                if (r.qcdAbcdy) wn(r.qcdAbcdy, "qcd_abcdy");
                // LHE shape systematics: <process>_<syst>Up / Down (Combine's
                // `shapes ... $PROCESS_$SYSTEMATIC` convention)
                for (auto &kv : r.syst) wn(kv.second, kv.first.c_str());
            };

            if ((int)regs.size() != 4 * NY)
                std::cerr << "[WARN] " << combineOut << ": expected " << 4 * NY
                          << " per-(charge,y) regions, got " << regs.size()
                          << " -- input skim histograms incomplete? File may be unfittable.\n";
            for (auto &r : regs) writeRegionTemplates(r);
            writeRegionTemplates(Wp_incl);
            writeRegionTemplates(Wm_incl);
            writeRegionTemplates(W_incl);

            // --- in-fit ABCD CR channel dirs (leppt_mt40 only) ---
            // (+ since 2026-10-19 the abcdy group CR dirs, written by the SAME
            // lambda after the 6 inclusive ones; nCR counts the inclusive
            // ones as before, nCRY the group ones)
            int nCR = 0, nCRY = 0;
            auto writeCR = [&](const CRTemplates &cr) -> bool
                {
                    TDirectory *d = fcomb->mkdir(cr.dir.c_str());
                    if (!d) { std::cerr << "[WARN] mkdir failed: " << cr.dir << "\n"; return false; }
                    auto wn1 = [&](TH1D *h, const char *nm) {
                        if (!h) return;
                        d->cd();
                        TH1D *hc = (TH1D *)h->Clone(nm);
                        hc->SetDirectory(d);
                        // same all-zero-shape guard as wn() above
                        if (strcmp(nm, "data_obs") != 0 && hc->Integral() <= 0.0)
                        {
                            std::cerr << "[WARN] " << combineOut << " " << cr.dir << "/" << nm
                                      << " is empty -> flooring to 1e-6 for Combine\n";
                            hc->SetBinContent(1, 1e-6);
                        }
                        hc->Write(nm, TObject::kOverwrite);
                    };
                    wn1(cr.data, "data_obs");
                    wn1(cr.qcd, "qcd");
                    wn1(cr.z, "z");
                    wn1(cr.ztau, "ztau");
                    wn1(cr.ttbar, "ttbar");
                    wn1(cr.wfix, "wfix");
                    wn1(cr.ewk, "ewk");
                    for (int iy = 0; iy < (int)cr.wLab.size(); ++iy) wn1(cr.wLab[iy], Form("w_lab_y%d", iy));
                    for (int iy = 0; iy < (int)cr.wFb.size(); ++iy)  wn1(cr.wFb[iy],  Form("w_fb_y%d", iy));
                    for (const auto &p : cr.wNamed) wn1(p.second, p.first.c_str()); // abcdy group CRB only
                    return true;
                };
            if (crs)
                for (const auto &cr : *crs)
                    if (writeCR(cr)) ++nCR;
            if (crsY)
                for (const auto &cr : *crsY)
                    if (writeCR(cr)) ++nCRY;

            fcomb->Close();
            delete fcomb;
            std::string crMsg = nCR ? Form(" + %d CR dirs", nCR) : "";
            if (nCRY) crMsg += Form(" + %d abcdy group CR dirs", nCRY);
            std::cout << "[INFO] Saved structured Combine input: " << combineOut
                      << "  (" << regs.size()
                      << " per-(charge,y) regions + Wp_incl/Wm_incl/W_incl"
                      << crMsg << ")\n";

            // --- shape-systematics sidecar (read by the fork's card generator) ---
            // <combineOut minus .root>_systs.txt: one line per systematic that is
            // COMPLETE (signal/z/ztau/wtau/ttbar x Up/Down) in EVERY per-(charge,y)
            // region, followed by the processes carrying it. An incomplete syst
            // is warned about and left out, so a card can never reference a
            // shape that is missing somewhere. Flavour-specific systematics (the
            // muon SFs) appear only in that flavour's sidecar -- the generator
            // puts them on that flavour's channels alone.
            {
                const std::string side = combineOut.substr(0, combineOut.size() - 5) + "_systs.txt";
                std::ofstream sf(side.c_str());
                sf << "# Shape systematics in " << combineOut << " (plotting/mtandmet.C)\n"
                   << "# <systematic> <processes carrying <process>_<systematic>Up/Down in every SR region>\n"
                   << "# a process token <hist>@<charge> restricts the row to that charge's channels\n";
                // directive for the fork's card generator, only when a separate muTrig
                // nuisance is shipped (today the three sources go in as ONE combined
                // `muSF`, so this is dormant): its correlation model follows the
                // trigger-SF binning chosen in skim/muon_sf.h (inclusive -> coherent)
                if (std::find(sideNames.begin(), sideNames.end(), std::string("muTrig")) != sideNames.end())
                    sf << "#! muTrig corr " << pOSF::kMuTrigCorr << "\n";
                int nListed = 0;
                for (const std::string &systStr : sideNames)
                {
                    const char *syst = systStr.c_str();
                    int nFull = 0;
                    for (auto &r : regs)
                    {
                        int have = 0;
                        for (const char *p : {"signal", "z", "ztau", "wtau", "ttbar"})
                            for (const char *d : {"Up", "Down"})
                                have += (int)r.syst.count(std::string(p) + "_" + syst + d);
                        if (have == 10) ++nFull;
                    }
                    if (nFull > 0 && nFull == (int)regs.size()) { sf << syst << " signal z ztau wtau ttbar\n"; ++nListed; }
                    else if (nFull > 0)
                        std::cerr << "[WARN] " << combineOut << ": shape syst " << syst << " complete in only "
                                  << nFull << "/" << regs.size() << " regions -> NOT listed in " << side << "\n";
                }
                // The electron fake-factor slope nuisances (2026-09-29): on the
                // data-driven qcd_abcd template only, ONE PER CHARGE -- the token
                // qcd_abcd@<charge> puts the fork's `shape` row on the qcd column
                // of that charge's channels alone (the two slopes are measured
                // on disjoint samples). Listed only when every region of the
                // charge carries both directions.
                // Since 2026-10-19 the same slope also rides the abcdy SR
                // template: when its twins are complete the row names BOTH
                // histograms, 'qcdFF_ele_Wp qcd_abcd@Wp qcd_abcdy@Wp' -- the
                // card generator puts the entry on whichever of the two is the
                // card's qcd column, so an abcd card is unchanged by it.
                for (int ic = 0; ic < 2; ++ic)
                {
                    const std::string cg  = (ic == 0) ? "Wp" : "Wm";
                    const std::string key = "qcd_abcd_" + ffSystName[ic];
                    const std::string keyY = "qcd_abcdy_" + ffSystName[ic];
                    int nCh = 0, nFull = 0, nFullY = 0;
                    for (auto &r : regs)
                    {
                        if (r.dir.compare(0, cg.size() + 1, cg + "_") != 0) continue;
                        ++nCh;
                        if (r.syst.count(key + "Up") && r.syst.count(key + "Down")) ++nFull;
                        if (r.syst.count(keyY + "Up") && r.syst.count(keyY + "Down")) ++nFullY;
                    }
                    if (nFullY > 0 && nFullY != nCh)
                        std::cerr << "[WARN] " << combineOut << ": " << ffSystName[ic] << " on qcd_abcdy complete in only "
                                  << nFullY << "/" << nCh << " " << cg << " regions -> qcd_abcdy NOT listed for it in "
                                  << side << "\n";
                    if (nFull == 0) continue;
                    if (nFull == nCh)
                    {
                        sf << ffSystName[ic] << " qcd_abcd@" << cg;
                        if (nFullY == nCh) sf << " qcd_abcdy@" << cg;
                        sf << "\n";
                        ++nListed;
                    }
                    else
                        std::cerr << "[WARN] " << combineOut << ": " << ffSystName[ic] << " complete in only "
                                  << nFull << "/" << nCh << " " << cg << " regions -> NOT listed in " << side << "\n";
                }
                std::cout << "[INFO] shape systematics in " << combineOut << ": " << nListed
                          << " listed in " << side << "\n";
            }
        };

        writeCombineInput(outBase + "/combine_input_W.root", regions, "", nullptr);
        // (combine_input_W_leppt.root -- the no-m_T lepton-pT fit, the nominal
        // since 2026-10-21 -- is written by plotting/combine_input_nomt.C; the
        // guard that deleted the retired plain-leppt input here is gone)
        // the m_T > 40 inputs, one per selection variant (met_mt40 since 2026-10-20)
        for (int var = 0; var < kNVar; ++var)
        {
        const std::string inName = std::string("combine_input_W_") + varDisc[var];
        std::vector<CRTemplates> &crRegions  = crRegionsV[var];
        std::vector<CRTemplates> &crRegionsY = crRegionsYV[var];
        if (!regionsPt[var].empty())
            writeCombineInput(outBase + "/" + inName + ".root",
                              regionsPt[var], std::string("_lp") + varTag[var],
                              crRegions.empty() ? nullptr : &crRegions,
                              crRegionsY.empty() ? nullptr : &crRegionsY,
                              &systNamesV[var]);
        else
        {
            gSystem->Unlink((outBase + "/" + inName + ".root").c_str());
            std::cerr << "[WARN] no " << varDisc[var] << " regions collected -> " << inName << ".root"
                      << " not written (any stale copy deleted so the fork cannot fit old templates)\n";
        }

        // --- the abcdy GROUP SIDECAR (2026-10-19) -------------------------------
        // combine_input_W_<leppt_mt40|met_mt40>_abcdy.txt, next to the input: one line per
        // (binning, group) with its bin indices -- the fork's card generator
        // reads it for QCD_MODE=abcdy (group CR channels, the per-group formula
        // rateParams and the CRB w_y maps) and refuses abcdy without it. Written
        // ONLY when everything abcdy is in the file: qcd_abcdy in all 48 SR
        // regions, the 36 group CR dirs, and (electron) the slope twins of
        // qcd_abcdy in every region; otherwise any stale copy is deleted.
        {
            const std::string ySide = outBase + "/" + inName + "_abcdy.txt";
            int nSR = 0, nFF = 0;
            for (const auto &r : regionsPt[var])
            {
                if (r.qcdAbcdy) ++nSR;
                const std::string cg = r.dir.substr(0, 2);
                const std::string key = "qcd_abcdy_" + ffSystName[cg == "Wp" ? 0 : 1];
                if (r.syst.count(key + "Up") && r.syst.count(key + "Down")) ++nFF;
            }
            const int nGrpCR = 2 * 2 * pOAbcdY::kNGroups * 3;
            const bool abcdyOk = abcdyOkV[var];
            const bool complete = abcdyOk && nSR == 4 * NY && (int)regionsPt[var].size() == 4 * NY &&
                                  (int)crRegionsY.size() == nGrpCR && (!isElec || nFF == 4 * NY);
            if (complete)
            {
                std::ofstream ys(ySide.c_str());
                ys << "# Rapidity groups of the abcdy in-fit ABCD (QCD_MODE=abcdy) in " << inName << ".root\n"
                   << "# (plotting/mtandmet.C; groups from correction/abcd_ygroups.h, counts from correction/qcd_abcd.C)\n"
                   << "# <binning> <group> <bin indices>: SR template <C>_<binning>_y<i>/qcd_abcdy, the group's control\n"
                   << "# regions <C>_<binning>_<group>_CR{B,C,D} (CRB W parts w_<binning>_y<i> of the group's bins)\n";
                for (int ib = 0; ib < 2; ++ib)
                    for (int g = 0; g < pOAbcdY::kNGroups; ++g)
                    {
                        ys << abcdyBin[ib] << " " << pOAbcdY::kGroupName[g];
                        for (int iy : pOAbcdY::BinsOf(g)) ys << " " << iy;
                        ys << "\n";
                    }
                std::cout << "[QCD-ABCDY] wrote the group sidecar " << ySide << " (" << nSR << " SR qcd_abcdy, "
                          << crRegionsY.size() << " group CR dirs" << (isElec ? Form(", %d slope-twin pairs", nFF) : "")
                          << ")\n";
            }
            else
            {
                gSystem->Unlink(ySide.c_str());
                std::cerr << "[WARN] abcdy inputs incomplete (abcdyOk " << abcdyOk << ", SR qcd_abcdy " << nSR << "/"
                          << 4 * NY << ", group CR dirs " << crRegionsY.size() << "/" << nGrpCR
                          << (isElec ? Form(", slope twins %d/%d", nFF, 4 * NY) : "") << ") -> no group sidecar "
                          << ySide << " (any stale copy deleted; QCD_MODE=abcdy cannot be built from this input)\n";
            }
        }
        } // per variant
    }

    f->Close();
    delete f;

    std::cout << "[INFO] Done. Saved plots under: " << outBase << "\n";
}