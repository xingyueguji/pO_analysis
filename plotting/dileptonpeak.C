#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "plotting_helper.C"
#include "../skim/mc_norm.h"   // pONorm::MCScale -> per-sample k_s = sigma*L/N_gen

#include <string>
#include <vector>
#include <functional>

void dileptonpeak(bool isElec = 0)
{
    const std::string baseDir = "../skim/rootfile/";
    std::string prefix;
    const char *Channeltype;

    if (isElec) {
        prefix      = "ZToEE_pO2025";
        Channeltype = "Z #rightarrow e e";
    } else {
        prefix      = "ZToMuMu_pO2025";
        Channeltype = "Z #rightarrow #mu #mu";
    }

    auto openFile = [](const std::string &p) -> TFile* {
        TFile *f = TFile::Open(p.c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "[ERROR] Cannot open file: " << p << "\n";
            return nullptr;
        }
        return f;
    };

    const int rebinFactor = 2;   // 120 -> 60 bins

    auto getRebinned = [&](TFile *file, const char *hname) -> TH1D* {
        if (!file) return nullptr;
        TH1D *src = (TH1D*)file->Get(hname);
        if (!src) { std::cerr << "[WARN] missing " << hname
                              << " in " << file->GetName() << "\n"; return nullptr; }
        static int counter = 0;
        TH1D *h = (TH1D*)src->Clone(Form("%s_copy%d", hname, counter++));
        h->SetDirectory(nullptr);                 // we own this copy now
        if (rebinFactor > 1) h->Rebin(rebinFactor);
        return h;
    };

    TFile *f       = openFile(baseDir + prefix + "_Data_hist.root");
    TFile *f_DY    = openFile(baseDir + prefix + "_DY_MC_hist.root");
    TFile *f_Wp    = openFile(baseDir + prefix + "_Wp_MC_hist.root");
    TFile *f_Wm    = openFile(baseDir + prefix + "_Wm_MC_hist.root");
    TFile *f_DYtau = openFile(baseDir + prefix + "_DYtau_MC_hist.root");
    TFile *f_Wptau = openFile(baseDir + prefix + "_Wptau_MC_hist.root");
    TFile *f_Wmtau = openFile(baseDir + prefix + "_Wmtau_MC_hist.root");
    if (!f || !f_DY || !f_Wp || !f_Wm || !f_DYtau || !f_Wptau || !f_Wmtau) return;

    const std::string outBase = isElec ? "./plots/Elec" : "./plots";
    const std::string outDir  = outBase + "/dilepton";
    gSystem->mkdir(outDir.c_str(), kTRUE);

    PlotStyle ps;
    ps.drawOpt   = "hist";
    ps.showStats = false;
    ps.logy      = true;
    ps.boxY1     = 0.62;
    ps.boxY2     = 0.82;
    ps.normBkgToData = false; // ABSOLUTE pO scaling (k_s incl. A=16) -- no area norm

    PlotTuner commonTuner = [&](TCanvas *c, TH1 *h) {
        (void)c;
        if (!h) return;
        if (ps.logy) h->SetMinimum(1.0);
        h->SetMaximum(10.25 * h->GetMaximum());
    };

    const char *xTitle = isElec ? "m_{ee} (GeV)" : "m_{#mu#mu} (GeV)";

    struct Spec { const char *hname; const char *yTitle; const char *suffix; };
    std::vector<Spec> specs = {
        {"hMass",          "Events / 1.0 GeV",   "mass"},
        {"hMass_extended", "Events / 1.75 GeV", "mass_extended"},
        {"hMass_vipul",    "Events / 1.0 GeV",   "mass_vipul"},
    };

    // --- Absolute per-sample MC normalization k_s = A*sigma*L/N_gen -----------
    // Scale each MC sample by its own k_s so the components sit in the correct
    // physical proportion AND at the absolute pO yield (drawn absolute, no area
    // norm -- see ps.normBkgToData=false above). sigma is already in the gen
    // weights (<w>=sigma, per-NN); MCScale multiplies by the Oxygen A-scaling
    // (A=16) and L. The tau samples share the W/DY cross sections but carry their
    // own N_gen (per-file label).
    const double k_Wp    = pONorm::MCScale(isElec ? "Wp_ele" : "Wp_mu");
    const double k_Wm    = pONorm::MCScale(isElec ? "Wm_ele" : "Wm_mu");
    const double k_DY    = pONorm::MCScale(isElec ? "DYee"   : "DYmu");
    const double k_DYtau = pONorm::MCScale("DYtau");
    const double k_Wptau = pONorm::MCScale("Wp_tau");
    const double k_Wmtau = pONorm::MCScale("Wm_tau");
    auto kScale = [](TH1D *h, double k) { if (h) h->Scale(k); };

        for (const auto &s : specs) {
        TH1D *h_data  = getRebinned(f,s.hname);
        TH1D *h_DYh   = getRebinned(f_DY,s.hname);
        TH1D *h_Wph   = getRebinned(f_Wp,s.hname);
        TH1D *h_Wmh   = getRebinned(f_Wm,s.hname);
        TH1D *h_DYth  = getRebinned(f_DYtau,s.hname);
        TH1D *h_Wpth  = getRebinned(f_Wptau,s.hname);
        TH1D *h_Wmth  = getRebinned(f_Wmtau,s.hname);

        if (!h_data) { std::cerr << "[WARN] Missing data hist " << s.hname << "\n"; continue; }

        // Apply per-sample absolute normalization (relative composition).
        // getRebinned already returned owned clones, so scaling in place is safe.
        kScale(h_Wph, k_Wp);     kScale(h_Wmh, k_Wm);
        kScale(h_DYh, k_DY);     kScale(h_DYth, k_DYtau);
        kScale(h_Wpth, k_Wptau); kScale(h_Wmth, k_Wmtau);

        // --- combine W+ and W- into a single "W+/W-" template ---
        TH1D *h_W = nullptr;
        if (h_Wph) {
            h_W = (TH1D*)h_Wph->Clone(Form("%s_WpWm_combined", s.hname));
            h_W->SetDirectory(nullptr);
            if (h_Wmh) h_W->Add(h_Wmh);
        } else if (h_Wmh) {
            h_W = (TH1D*)h_Wmh->Clone(Form("%s_WpWm_combined", s.hname));
            h_W->SetDirectory(nullptr);
        }

        // --- combine Wp tau and Wm tau into a single "W+/W- tau" template ---
        TH1D *h_Wtau = nullptr;
        if (h_Wpth) {
            h_Wtau = (TH1D*)h_Wpth->Clone(Form("%s_WpWmTau_combined", s.hname));
            h_Wtau->SetDirectory(nullptr);
            if (h_Wmth) h_Wtau->Add(h_Wmth);
        } else if (h_Wmth) {
            h_Wtau = (TH1D*)h_Wmth->Clone(Form("%s_WpWmTau_combined", s.hname));
            h_Wtau->SetDirectory(nullptr);
        }

        // (No area normalization: the MC is already at its absolute pO yield from
        // k_s above, and ps.normBkgToData=false keeps SaveNicePlot1D_WithBkg from
        // rescaling the stack to data.)

        std::vector<std::string> box = {
            Form("Passing Events: %.0f", h_data->Integral(1, h_data->GetNbinsX()))
        };
        std::vector<TH1*> bkgs         = { h_W,      h_DYh, h_DYth,   h_Wtau };
        std::vector<std::string> names = { "W+/W-",  "DY",  "DY tau", "W+/W- tau" };

        SaveNicePlot1D_WithBkg(
            h_data, bkgs, names,
            outDir + "/" + s.suffix,
            xTitle, s.yTitle,
            "",
            Channeltype,
            "inclusive",
            box, ps, commonTuner);
    }

        // --- Combine input ROOT file for Z fit ---
    {
        const std::string combineOut = outBase + "/combine_input_dilepton.root";
        TFile *fout = TFile::Open(combineOut.c_str(), "RECREATE");
        if (!fout || fout->IsZombie()) {
            std::cerr << "[ERROR] Cannot create output file: " << combineOut << "\n";
        } else {
            auto write_clone = [&](TH1D *h, const char *name) {
                if (!h) { std::cerr << "[WARN] Missing histogram: " << name << "\n"; return; }
                fout->cd();
                TH1D *hc = (TH1D*)h->Clone(name);
                hc->SetDirectory(fout);
                hc->Write(name, TObject::kOverwrite);
            };

            // Pick the histogram you want to fit (60-120 GeV signal window).
            // Swap to "hMass_vipul" or "hMass_extended" if that's the preferred fit input.
            const char *hname = "hMass";

            TH1D *h_data  = getRebinned(f,hname);
            TH1D *h_DYh   = getRebinned(f_DY,hname);
            TH1D *h_Wph   = getRebinned(f_Wp,hname);
            TH1D *h_Wmh   = getRebinned(f_Wm,hname);
            TH1D *h_DYth  = getRebinned(f_DYtau,hname);
            TH1D *h_Wpth  = getRebinned(f_Wptau,hname);
            TH1D *h_Wmth  = getRebinned(f_Wmtau,hname);

            // Same per-sample normalization as the plots above (physical
            // composition; harmless for the shape-only datacard, correct later).
            kScale(h_Wph, k_Wp);     kScale(h_Wmh, k_Wm);
            kScale(h_DYh, k_DY);     kScale(h_DYth, k_DYtau);
            kScale(h_Wpth, k_Wptau); kScale(h_Wmth, k_Wmtau);

            // Combine W+/W- into one template
            TH1D *h_W = nullptr;
            if (h_Wph) {
                h_W = (TH1D*)h_Wph->Clone("W_combined");
                h_W->SetDirectory(nullptr);
                if (h_Wmh) h_W->Add(h_Wmh);
            } else if (h_Wmh) {
                h_W = (TH1D*)h_Wmh->Clone("W_combined");
                h_W->SetDirectory(nullptr);
            }

            // Combine Wp tau + Wm tau into one template
            TH1D *h_Wtau = nullptr;
            if (h_Wpth) {
                h_Wtau = (TH1D*)h_Wpth->Clone("Wtau_combined");
                h_Wtau->SetDirectory(nullptr);
                if (h_Wmth) h_Wtau->Add(h_Wmth);
            } else if (h_Wmth) {
                h_Wtau = (TH1D*)h_Wmth->Clone("Wtau_combined");
                h_Wtau->SetDirectory(nullptr);
            }

            // Templates
            write_clone(h_data, "data_obs");
            write_clone(h_DYh,  "signal"); // DY -> ll is the Z signal
            write_clone(h_W,    "w");      // W+/W- -> l nu
            write_clone(h_Wtau, "wtau");   // W+/W- -> tau nu
            write_clone(h_DYth, "ztau");   // DY -> tau tau

            fout->Close();
            delete fout;
            std::cout << "[INFO] Saved Combine hist file: " << combineOut << "\n";
        }
    }

    f->Close();         delete f;
    f_DY->Close();      delete f_DY;
    f_Wp->Close();      delete f_Wp;
    f_Wm->Close();      delete f_Wm;
    f_DYtau->Close();   delete f_DYtau;
    f_Wptau->Close();   delete f_Wptau;
    f_Wmtau->Close();   delete f_Wmtau;

    std::cout << "[INFO] Done. Saved plots under: " << outBase << "\n";
}