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
    gSystem->mkdir(outMetDir.c_str(), kTRUE);
    gSystem->mkdir(outMtDir.c_str(), kTRUE);
    gSystem->mkdir(outLepPtMt40Dir.c_str(), kTRUE);

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

    auto scaleAll = [](double k, std::initializer_list<TH1D *> hs)
    {
        for (TH1D *h : hs)
            if (h) h->Scale(k);
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

    // --- ABCD QCD (low-MET) background, MET stacks ----------------------------
    // Inclusive-in-rapidity per-charge QCD MET template from
    // correction/qcd_abcd.C (ABCD: data-EWK anti-iso shape x transfer factor),
    // muon OR electron (qcd_abcd_{mu,ele}.root). Added to the MET stacks below.
    // The ABCD normalization is only measured INCLUSIVELY (low pO statistics), so
    // for the per-rapidity-bin plots the one inclusive template is split across y
    // by a data-driven proxy: the per-bin low-MET (MET < qcdMetCut) excess
    // (data - EWK), i.e. where the QCD sits in y, keeping the inclusive template
    // SHAPE. The per-y QCD therefore SUMS back to the inclusive template by
    // construction. This is an approximation for the per-bin plots; the inclusive
    // MET plot is the rigorous one.
    const double qcdMetCut = 30.0; // matches metCut in correction/qcd_abcd.C
    const std::string qlep = isElec ? "ele" : "mu";
    const std::string qcdFile = "../correction/rootfile/qcd_abcd_" + qlep + ".root";
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
            if (pp) { qcdPtPlusBase  = (TH1D *)pp->Clone("qcd_pt_Plus_base");  qcdPtPlusBase->SetDirectory(nullptr); }
            if (pm) { qcdPtMinusBase = (TH1D *)pm->Clone("qcd_pt_Minus_base"); qcdPtMinusBase->SetDirectory(nullptr); }
            // in-fit ABCD counts (m_T plane; written by qcd_abcd.C since 2026-08-23)
            TH1D *cp = (TH1D *)fqcd->Get(Form("abcd_counts_%sPlus", qlep.c_str()));
            TH1D *cm = (TH1D *)fqcd->Get(Form("abcd_counts_%sMinus", qlep.c_str()));
            if (cp) { abcdCountsPlus  = (TH1D *)cp->Clone("abcd_counts_Plus_loc");  abcdCountsPlus->SetDirectory(nullptr); }
            if (cm) { abcdCountsMinus = (TH1D *)cm->Clone("abcd_counts_Minus_loc"); abcdCountsMinus->SetDirectory(nullptr); }
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
        qcdAbcdPlusBase  = mkAbcdBase(qcdPtMt40PlusBase,  abcdCountsPlus,  "qcd_abcd_Plus_base",  "W+");
        qcdAbcdMinusBase = mkAbcdBase(qcdPtMt40MinusBase, abcdCountsMinus, "qcd_abcd_Minus_base", "W-");
        if (!qcdAbcdPlusBase || !qcdAbcdMinusBase)
            std::cerr << "[WARN] abcd_counts_* missing/incomplete in " << qcdFile
                      << " -- re-run correction/run_qcd_abcd.sh. qcd_abcd + CR templates"
                      << " will NOT be written (QCD_MODE=abcd unavailable for this input).\n";
    }

    // Per-y QCD weights (normalized to 1, clamped >=0) from the low-MET excess.
    auto qcdWeights = [&](const char *chg, const char *suf) -> std::vector<double>
    {
        std::vector<double> w(NY, 0.0);
        double sum = 0.0;
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
                         lowInt(f_DYtau, k_DYtau) + lowInt(f_Wptau, k_Wptau) + lowInt(f_Wmtau, k_Wmtau);
            double e = d - ewk;
            if (e < 0) e = 0; // EWK over-subtraction in a sparse bin -> no QCD there
            w[iy] = e;
            sum += e;
        }
        for (int iy = 0; iy < NY; ++iy) w[iy] = (sum > 0) ? w[iy] / sum : 1.0 / NY;
        return w;
    };

    const std::vector<double> wQcdWp   = qcdPlusBase  ? qcdWeights("Wp", "")    : std::vector<double>(NY, 0.0);
    const std::vector<double> wQcdWm   = qcdMinusBase ? qcdWeights("Wm", "")    : std::vector<double>(NY, 0.0);
    const std::vector<double> wQcdWpFB = qcdPlusBase  ? qcdWeights("Wp", "_FB") : std::vector<double>(NY, 0.0);
    const std::vector<double> wQcdWmFB = qcdMinusBase ? qcdWeights("Wm", "_FB") : std::vector<double>(NY, 0.0);

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

    // Lepton-pT selection variants (shared by the stacks below and the Combine
    // input writer). ONE variant since 2026-09-21: kVarMt40 = pT>25 && m_T>40.
    // The plain no-m_T-cut variant (kVarNom, "h_leppt_*" -> plots/leppt/ +
    // combine_input_W_leppt.root) was RETIRED that day: the discriminant itself
    // was dropped on 2026-08-16 (only met and leppt_mt40 are fitted), its
    // out-tree pO_fit_out_leppt/ never had a simfit, and the skim stores no LHE
    // /SF systematic twins for h_leppt_* -- so its Combine input carried no
    // shape systematics and could not be fitted by the current card generator.
    // The table is kept so re-adding a variant is a two-entry edit.
    enum { kVarMt40 = 0, kNVar = 1 };
    const char *varStem[kNVar] = {"h_leppt_mt40"};
    const char *varTag [kNVar] = {"_mt40"};

    // --- Structured per-region Combine input -----------------------------------
    // For every fit region we emit the 6 absolute templates the downstream
    // Combine datacards consume -- signal, z, ztau, wtau, qcd, data_obs -- each
    // in its own TDirectory. Regions:
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
        // Shape systematics (2026-09-07 LHE, 2026-09-14 muon SFs): "<process>_<syst>Up|Down"
        // -> the absolute varied template, process in {signal, z, ztau, wtau},
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
                          TH1D *z, TH1D *ztau, TH1D *wtA, TH1D *wtB,
                          TH1D *qcd, TH1D *qcdAbcd = nullptr) -> RegionTemplates {
        RegionTemplates r;
        r.dir  = dir;
        r.data = cloneDetached(data, (dir + tag + "_data_obs").c_str());
        r.sig  = sum2(sigA, sigB,    (dir + tag + "_signal").c_str());
        r.z    = cloneDetached(z,    (dir + tag + "_z").c_str());
        r.ztau = cloneDetached(ztau, (dir + tag + "_ztau").c_str());
        r.wtau = sum2(wtA, wtB,      (dir + tag + "_wtau").c_str());
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

    // Get + absolute-scale one MC histo (fresh Get per name/file, scaled in
    // place and used once -- same convention as scaleAll above).
    auto getScaled = [](TFile *ff, double k, const char *name) -> TH1D *
    {
        TH1D *h = (TH1D *)ff->Get(name);
        if (h) h->Scale(k);
        return h;
    };

    // The shape systematics this channel's MC skim files carry: the LHE families
    // (both flavours; skim/lhe_updown.py) + the lepton-SF families of THIS
    // flavour (muon: skim/muon_sf.h, written by skim.C; electron: none yet).
    // Single list for attachSysts() and the sidecar below.
    std::vector<std::string> systNames(pOLhe::kLheSystNames, pOLhe::kLheSystNames + pOLhe::kNLheSysts);
    if (!isElec)
        systNames.insert(systNames.end(), pOSF::kMuonSFSystNames, pOSF::kMuonSFSystNames + pOSF::kNMuonSFSysts);

    // Shape systematics for one region (2026-09-07). For every syst in
    // `systNames` and direction Up/Down, read the per-sample
    // "<baseName>_<syst><dir>" twins (skim/lhe_updown.py for the LHE families,
    // skim.C for the muon SFs) from the MC skim files, scale each by its own
    // k_s and combine them EXACTLY like the nominal templates: signal = W+
    // sample + W- sample reconstructed with this region's charge, wtau = both
    // tau samples, z / ztau single samples. A syst is attached only when all
    // six inputs exist in BOTH directions (otherwise [WARN] once per syst and
    // skipped); the sidecar written next to each Combine input lists what
    // every region carries, and the fork's card generator must not reference
    // more. Each varied name is Get()'d exactly once per file (region names
    // are unique), so the in-place Scale is safe.
    std::set<std::string> lheWarned;
    auto attachSysts = [&](RegionTemplates &r, const std::string &tag, const char *baseName)
    {
        for (const std::string &systStr : systNames)
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
                if (!hWp || !hWm || !hZ || !hZt || !hWpt || !hWmt) { complete = false; break; }
                const std::string key = std::string("_") + syst + dir; // e.g. "_nPDFUp"
                const std::string pre = r.dir + tag;
                got["signal" + key] = sum2(hWp, hWm,   (pre + "_signal" + key).c_str());
                got["z" + key]      = cloneDetached(hZ,  (pre + "_z" + key).c_str());
                got["ztau" + key]   = cloneDetached(hZt, (pre + "_ztau" + key).c_str());
                got["wtau" + key]   = sum2(hWpt, hWmt, (pre + "_wtau" + key).c_str());
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
    auto lepPtStack = [&](int var, int iy, const char *chg, const char *suf,
                          TH1D *qcdH, const char *sub1, const std::string &sub2,
                          const std::string &outPath, bool accumulate,
                          TH1D *qcdAbcdH = nullptr) -> bool
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

        std::vector<std::string> box = {
            Form("Passing Events: %.0f", hD->Integral(1, hD->GetNbinsX())),
            Form("W signal MC: %.0f", sigInt({hWp, hWm})),
            totalLine(bkgs, names)};

        // pT axes start at the selection floor (bin edge 24 encloses the
        // 25 GeV cut on the 2 GeV grid) -- no empty [0,25) band. The helper
        // restores the full range after saving, so the Combine-input writing
        // below sees untouched histograms.
        PlotStyle psPt = psCounts;
        psPt.xRangeLo = 24.0;
        psPt.xRangeHi = 100.0;
        SaveNicePlot1D_WithBkg(hD, bkgs, names, outPath,
                               ptTitle, "Events / 2.0 GeV", "",
                               sub1, sub2, box, psPt, commonTuner);

        // Collect this region's templates for the lepton-pT Combine input
        // (written to the SEPARATE combine_input_W_leppt_mt40.root below --
        // combine_input_W.root itself stays PF-MET-only).
        regionsPt[var].push_back(makeRegion(
            Form("%s_%s_y%d", chg, (suf[0] ? "fb" : "lab"), iy),
            Form("_lp%s", varTag[var]),
            hD, hWp, hWm, hZ, hZtau, hWptau, hWmtau, qcdH, qcdAbcdH));
        // LHE / lepton-SF shape systematics: the skim stores twins for h_met_*
        // and h_leppt_mt40_* only -- which is also why the plain leppt variant
        // was retired (2026-09-21): it could carry none.
        attachSysts(regionsPt[var].back(), Form("_lp%s", varTag[var]), nm.c_str());

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
            h_met_Wp_MC_Z, h_met_Wp_MC_Ztau, h_met_Wp_MC_Wptau, h_met_Wp_MC_Wmtau,
            qcd_met_Wp));
        attachSysts(regions.back(), "", Form("h_met_Wp_y%d", iy));
        regions.push_back(makeRegion(Form("Wm_lab_y%d", iy), "",
            h_met_Wm, h_met_Wm_MC_Wp, h_met_Wm_MC_Wm,
            h_met_Wm_MC_Z, h_met_Wm_MC_Ztau, h_met_Wm_MC_Wptau, h_met_Wm_MC_Wmtau,
            qcd_met_Wm));
        attachSysts(regions.back(), "", Form("h_met_Wm_y%d", iy));
        regions.push_back(makeRegion(Form("Wp_fb_y%d", iy), "",
            h_met_Wp_FB, h_met_Wp_FB_MC_Wp, h_met_Wp_FB_MC_Wm,
            h_met_Wp_FB_MC_Z, h_met_Wp_FB_MC_Ztau, h_met_Wp_FB_MC_Wptau, h_met_Wp_FB_MC_Wmtau,
            qcd_met_Wp_FB));
        attachSysts(regions.back(), "", Form("h_met_Wp_y%d_FB", iy));
        regions.push_back(makeRegion(Form("Wm_fb_y%d", iy), "",
            h_met_Wm_FB, h_met_Wm_FB_MC_Wp, h_met_Wm_FB_MC_Wm,
            h_met_Wm_FB_MC_Z, h_met_Wm_FB_MC_Ztau, h_met_Wm_FB_MC_Wptau, h_met_Wm_FB_MC_Wmtau,
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

        // --------- lepton-pT plots (+ the leppt_mt40 Combine input) -------------
        {
            // Per-y QCD = inclusive pT template shape x the SAME per-bin low-MET
            // excess weights as the MET stacks (the QCD y-distribution does not
            // depend on which variable is plotted). The m_T cut does not change
            // where QCD sits in y, so the m_T>40 split reuses those weights.
            TH1D *qcd_pt40_Wp    = qcdPerY(qcdPtMt40PlusBase,  wQcdWp[iy],   Form("qcd_pt_mt40_Wp_y%d", iy));
            TH1D *qcd_pt40_Wm    = qcdPerY(qcdPtMt40MinusBase, wQcdWm[iy],   Form("qcd_pt_mt40_Wm_y%d", iy));
            TH1D *qcd_pt40_Wp_FB = qcdPerY(qcdPtMt40PlusBase,  wQcdWpFB[iy], Form("qcd_pt_mt40_Wp_y%d_FB", iy));
            TH1D *qcd_pt40_Wm_FB = qcdPerY(qcdPtMt40MinusBase, wQcdWmFB[iy], Form("qcd_pt_mt40_Wm_y%d_FB", iy));

            // in-fit ABCD (abcd-mode) per-y templates: same shape, same per-y
            // weights (Sum_y = A0 per charge+binning since Sum w = 1), total
            // renormalized to A0. Written as the 7th SR template `qcd_abcd`;
            // deliberately NOT drawn in the display stacks.
            TH1D *qcdA_Wp    = qcdPerY(qcdAbcdPlusBase,  wQcdWp[iy],   Form("qcd_abcd_Wp_y%d", iy));
            TH1D *qcdA_Wm    = qcdPerY(qcdAbcdMinusBase, wQcdWm[iy],   Form("qcd_abcd_Wm_y%d", iy));
            TH1D *qcdA_Wp_FB = qcdPerY(qcdAbcdPlusBase,  wQcdWpFB[iy], Form("qcd_abcd_Wp_y%d_FB", iy));
            TH1D *qcdA_Wm_FB = qcdPerY(qcdAbcdMinusBase, wQcdWmFB[iy], Form("qcd_abcd_Wm_y%d_FB", iy));

            const std::string mtLab = " , m_{T} > 40 GeV";
            bool ok40 = true;
            ok40 &= lepPtStack(kVarMt40, iy, "Wp", "",    qcd_pt40_Wp,    Channeltypewplus,  yLabel[iy] + mtLab,
                               outLepPtMt40Dir + Form("/leppt_mt40_Wp_y%d", iy),    /*accumulate=*/true,  qcdA_Wp);
            ok40 &= lepPtStack(kVarMt40, iy, "Wm", "",    qcd_pt40_Wm,    Channeltypewminus, yLabel[iy] + mtLab,
                               outLepPtMt40Dir + Form("/leppt_mt40_Wm_y%d", iy),    /*accumulate=*/true,  qcdA_Wm);
            ok40 &= lepPtStack(kVarMt40, iy, "Wp", "_FB", qcd_pt40_Wp_FB, Channeltypewplus,  yLabel_FB[iy] + mtLab,
                               outLepPtMt40Dir + Form("/leppt_mt40_Wp_y%d_FB", iy), /*accumulate=*/false, qcdA_Wp_FB);
            ok40 &= lepPtStack(kVarMt40, iy, "Wm", "_FB", qcd_pt40_Wm_FB, Channeltypewminus, yLabel_FB[iy] + mtLab,
                               outLepPtMt40Dir + Form("/leppt_mt40_Wm_y%d_FB", iy), /*accumulate=*/false, qcdA_Wm_FB);
            if (!ok40 && iy == 0)
                std::cerr << "[WARN] h_leppt_mt40_* not found (skim predates the m_T>40"
                          << " lepton-pT histos); m_T>40 stacks skipped.\n";
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
        pushIf(qcd_pt_mt40_incl, "QCD (ABCD)");
        box.push_back(totalLine(bkgs, names));

        // pT axis starts at the selection floor (see lepPtStack)
        PlotStyle psPt = psCounts;
        psPt.xRangeLo = 24.0;
        psPt.xRangeHi = 100.0;
        SaveNicePlot1D_WithBkg(
            h_leppt_inclusive[var],
            bkgs,
            names,
            outLepPtMt40Dir + "/h_leppt_mt40_inclusive",
            ptTitle,
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
    // carry one frozen `ewk` (<= 2.5% of those regions).
    struct CRTemplates {
        std::string dir;
        TH1D *data = nullptr, *qcd = nullptr, *z = nullptr, *ztau = nullptr,
             *wfix = nullptr, *ewk = nullptr;
        std::vector<TH1D *> wLab, wFb;
    };
    std::vector<CRTemplates> crRegions;
    if (qcdAbcdPlusBase && qcdAbcdMinusBase && !regionsPt[kVarMt40].empty())
    {
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
                for (const auto &r : regionsPt[kVarMt40])
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
                  << " CR channel dirs (CRB/CRC/CRD x Wp/Wm) for combine_input_W_leppt_mt40.root\n";
    }

    {
        // --- Structured per-region Combine input files -------------------------
        // One TDirectory per fit region; each holds the 6 absolute templates
        // (7 for leppt_mt40: + qcd_abcd) and, for leppt_mt40, the 6 CR dirs.
        // Written once per discriminant: PF MET (the nominal fit input) and the
        // two lepton-pT variants (SEPARATE files -- same region names with a
        // different discriminant must never share a file).
        auto writeCombineInput = [&](const std::string &combineOut,
                                     std::vector<RegionTemplates> &regs,
                                     const std::string &tag,
                                     const std::vector<CRTemplates> *crs) {
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
                    add(s.qcd,  r.qcd,  "_qcd");
                    add(s.qcdAbcd, r.qcdAbcd, "_qcd_abcd");
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
            W_incl.qcd  = sum2(Wp_incl.qcd,  Wm_incl.qcd,  ("W_incl" + tag + "_qcd").c_str());
            W_incl.qcdAbcd = sum2(Wp_incl.qcdAbcd, Wm_incl.qcdAbcd, ("W_incl" + tag + "_qcd_abcd").c_str());
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
                if (r.qcd) wn(r.qcd, "qcd");
                if (r.qcdAbcd) wn(r.qcdAbcd, "qcd_abcd");
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
            int nCR = 0;
            if (crs)
                for (const auto &cr : *crs)
                {
                    TDirectory *d = fcomb->mkdir(cr.dir.c_str());
                    if (!d) { std::cerr << "[WARN] mkdir failed: " << cr.dir << "\n"; continue; }
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
                    wn1(cr.wfix, "wfix");
                    wn1(cr.ewk, "ewk");
                    for (int iy = 0; iy < (int)cr.wLab.size(); ++iy) wn1(cr.wLab[iy], Form("w_lab_y%d", iy));
                    for (int iy = 0; iy < (int)cr.wFb.size(); ++iy)  wn1(cr.wFb[iy],  Form("w_fb_y%d", iy));
                    ++nCR;
                }

            fcomb->Close();
            delete fcomb;
            std::cout << "[INFO] Saved structured Combine input: " << combineOut
                      << "  (" << regs.size()
                      << " per-(charge,y) regions + Wp_incl/Wm_incl/W_incl"
                      << (nCR ? Form(" + %d CR dirs", nCR) : "") << ")\n";

            // --- shape-systematics sidecar (read by the fork's card generator) ---
            // <combineOut minus .root>_systs.txt: one line per systematic that is
            // COMPLETE (signal/z/ztau/wtau x Up/Down) in EVERY per-(charge,y)
            // region, followed by the processes carrying it. An incomplete syst
            // is warned about and left out, so a card can never reference a
            // shape that is missing somewhere. Flavour-specific systematics (the
            // muon SFs) appear only in that flavour's sidecar -- the generator
            // puts them on that flavour's channels alone.
            {
                const std::string side = combineOut.substr(0, combineOut.size() - 5) + "_systs.txt";
                std::ofstream sf(side.c_str());
                sf << "# Shape systematics in " << combineOut << " (plotting/mtandmet.C)\n"
                   << "# <systematic> <processes carrying <process>_<systematic>Up/Down in every SR region>\n";
                // directive for the fork's card generator, only when a separate muTrig
                // nuisance is shipped (today the three sources go in as ONE combined
                // `muSF`, so this is dormant): its correlation model follows the
                // trigger-SF binning chosen in skim/muon_sf.h (inclusive -> coherent)
                if (std::find(systNames.begin(), systNames.end(), std::string("muTrig")) != systNames.end())
                    sf << "#! muTrig corr " << pOSF::kMuTrigCorr << "\n";
                int nListed = 0;
                for (const std::string &systStr : systNames)
                {
                    const char *syst = systStr.c_str();
                    int nFull = 0;
                    for (auto &r : regs)
                    {
                        int have = 0;
                        for (const char *p : {"signal", "z", "ztau", "wtau"})
                            for (const char *d : {"Up", "Down"})
                                have += (int)r.syst.count(std::string(p) + "_" + syst + d);
                        if (have == 8) ++nFull;
                    }
                    if (nFull > 0 && nFull == (int)regs.size()) { sf << syst << " signal z ztau wtau\n"; ++nListed; }
                    else if (nFull > 0)
                        std::cerr << "[WARN] " << combineOut << ": shape syst " << syst << " complete in only "
                                  << nFull << "/" << regs.size() << " regions -> NOT listed in " << side << "\n";
                }
                std::cout << "[INFO] shape systematics in " << combineOut << ": " << nListed
                          << " listed in " << side << "\n";
            }
        };

        writeCombineInput(outBase + "/combine_input_W.root", regions, "", nullptr);
        // The plain-leppt input was retired 2026-09-21 (see the variant table).
        // Delete any pre-retirement copy so the fork can never fit templates
        // that nothing regenerates -- the same guard the two live inputs use.
        for (const char *stale : {"/combine_input_W_leppt.root",
                                  "/combine_input_W_leppt_systs.txt"})
            gSystem->Unlink((outBase + stale).c_str());
        if (!regionsPt[kVarMt40].empty())
            writeCombineInput(outBase + "/combine_input_W_leppt_mt40.root",
                              regionsPt[kVarMt40], "_lp_mt40",
                              crRegions.empty() ? nullptr : &crRegions);
        else
        {
            gSystem->Unlink((outBase + "/combine_input_W_leppt_mt40.root").c_str());
            std::cerr << "[WARN] no m_T>40 lepton-pT regions collected -> combine_input_W_leppt_mt40.root"
                      << " not written (any stale copy deleted so the fork cannot fit old templates)\n";
        }
    }

    f->Close();
    delete f;

    std::cout << "[INFO] Done. Saved plots under: " << outBase << "\n";
}