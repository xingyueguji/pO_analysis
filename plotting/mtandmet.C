#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TSystem.h"
#include "plotting_helper.C"
#include "../skim/mc_norm.h"   // pONorm::MCScale -> per-sample k_s = sigma*L/N_gen

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
    const std::string outMetDir = outBase + "/met";
    const std::string outMtDir = outBase + "/mt";
    gSystem->mkdir(outMetDir.c_str(), kTRUE);
    gSystem->mkdir(outMtDir.c_str(), kTRUE);

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
    {
        TFile *fqcd = TFile::Open(qcdFile.c_str(), "READ");
        if (fqcd && !fqcd->IsZombie())
        {
            TH1D *qp = (TH1D *)fqcd->Get(Form("qcd_met_%sPlus", qlep.c_str()));
            TH1D *qm = (TH1D *)fqcd->Get(Form("qcd_met_%sMinus", qlep.c_str()));
            if (qp) { qcdPlusBase  = (TH1D *)qp->Clone("qcd_met_Plus_base");  qcdPlusBase->SetDirectory(nullptr); }
            if (qm) { qcdMinusBase = (TH1D *)qm->Clone("qcd_met_Minus_base"); qcdMinusBase->SetDirectory(nullptr); }
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
                ps,
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
                ps,
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
                ps,
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
                ps,
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
            ps,
            commonTuner);
    }

    {
        std::string combineOut = outBase + "/combine_input_inclusive.root";
        TFile *fout = TFile::Open(combineOut.c_str(), "RECREATE");
        if (!fout || fout->IsZombie())
        {
            std::cerr << "[ERROR] Cannot create output file: " << combineOut << "\n";
            return;
        }

        auto write_clone = [&](TH1D *h, const char *name)
        {
            if (!h)
            {
                std::cerr << "[WARN] Missing histogram: " << name << "\n";
                return;
            }
            fout->cd();
            TH1D *hc = (TH1D *)h->Clone(name);
            hc->SetDirectory(fout);
            hc->Write(name, TObject::kOverwrite);
        };

        // Data
        write_clone(h_met_inclusive, "data_obs");

        // Main templates
        write_clone(h_met_inclusive_MC_signal, "signal");
        write_clone(h_met_inclusive_MC_Z, "z");
        write_clone(h_met_inclusive_MC_Ztau, "ztau");
        write_clone(h_met_inclusive_MC_Wtau, "wtau");
        // Data-driven ABCD QCD (inclusive MET template; null/skipped for electron)
        if (qcd_met_incl) write_clone(qcd_met_incl, "qcd");

        fout->Close();
        delete fout;

        std::cout << "[INFO] Saved Combine hist file: " << combineOut << "\n";
    }

    f->Close();
    delete f;

    std::cout << "[INFO] Done. Saved plots under: " << outBase << "\n";
}