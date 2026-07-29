// correction/dataMC_kinematics.C
//
// Data vs signal-MC kinematic check, with a ratio pad.
//
//   Z channels (Zmm/Zee): boson pT / eta / y / phi, lepton pT / eta / phi, and
//     the dilepton mass. Signal MC = DY.
//   W channels (Wmu/Wel): leading-lepton pT / eta / phi after the FULL W
//     selection (isolation + trigger match included). Signal MC = Wp + Wm,
//     combined with their relative absolute scales k_s from skim/mc_norm.h
//     (only the W+/W- RATIO matters here since the sum is shape-normalized).
//
// This is the "Data/MC pT check" used to validate (and later derive) the boson
// pT reweighting; the W lepton-pT comparison also probes whether the lepton pT
// spectrum is modeled well enough to serve as a fit discriminant (Jacobian
// peak at ~m_W/2, insensitive to the unembedded-MC MET mismodeling).
//
// The comparison is shape-only: the signal MC is scaled to the data integral,
// so the ratio pad probes shape differences. Backgrounds are negligible inside
// the Z peak; in the W channels they are NOT subtracted, and because there is
// NO MET cut in the W selection the QCD fraction here is much larger than in
// the high-MET fit region: data - absolute EWK MC (2026-07-29, Delta-beta skim)
// gives QCD ~19% (mu) / ~50% (e) of the plotted sample, plus DY+tau ~7% / ~5%.
// (The familiar "QCD ~3% mu / ~10% e" figures are HIGH-MET-region numbers.)
// QCD peaks at low lepton pT, so expect a data excess at pT ~ 25-32 GeV; read
// the W ratio as a shape check, not a precision statement. A QCD-subtracted
// version needs the anti-iso lepton-pT shape (relIso x lepPt 2D, not stored
// yet -- see the lepton-pT-discriminant plan).
//
// Inputs  : Z: ../skim/rootfile/<base>_Data_hist.root   (data)
//              ../skim/rootfile/<base>_DY_MC_hist.root  (DY signal MC)
//              where <base> = ZToMuMu_pO2025 (Zmm) or ZToEE_pO2025 (Zee).
//           W: ../skim/rootfile/<base>_hist.root        (data)
//              ../skim/rootfile/<base>_{Wp,Wm}_hist.root (W signal MC)
//              where <base> = WToMuNu_pO_PFMet (Wmu) or WToElecNu_pO_PFMet (Wel).
// Outputs : ./plots/dataMC_<channel>/<hist>.png/.pdf   (i.e. correction/plots/)
//
// Histograms (produced by skim.C): Z: hMass, h_Zpt, h_Zeta, h_Zy, h_Zphi,
// h_lepPt, h_lepEta, h_lepPhi (iso-selected OS pairs in the Z peak [60,120]);
// W: h_lepPt, h_lepEta, h_lepPhi (leading lepton, full W selection).
//
// Run from this directory (correction/):
//   root -l -q 'dataMC_kinematics.C+("Zmm")'
//   root -l -q 'dataMC_kinematics.C+("Zee")'
//   root -l -q 'dataMC_kinematics.C+("Wmu")'
//   root -l -q 'dataMC_kinematics.C+("Wel")'

#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"

#include <iostream>
#include <string>
#include <vector>

#include "../plotting/plotting_helper.C"
#include "../skim/mc_norm.h" // pONorm::MCScale -> per-sample k_s (W+/W- mix)

void dataMC_kinematics(const char *channel = "Zmm")
{
    const std::string ch = channel;
    const bool isW = (ch == "Wmu" || ch == "Wel");

    std::string base, lepSym, bosSym;
    if (ch == "Zmm")
    {
        base = "ZToMuMu_pO2025";
        lepSym = "#mu";
        bosSym = "#mu#mu";
    }
    else if (ch == "Zee")
    {
        base = "ZToEE_pO2025";
        lepSym = "e";
        bosSym = "ee";
    }
    else if (ch == "Wmu")
    {
        base = "WToMuNu_pO_PFMet";
        lepSym = "#mu";
    }
    else if (ch == "Wel")
    {
        base = "WToElecNu_pO_PFMet";
        lepSym = "e";
    }
    else
    {
        std::cerr << "[ERR] channel must be \"Zmm\", \"Zee\", \"Wmu\" or \"Wel\" (got \"" << ch << "\")\n";
        return;
    }

    // Data file: the W skim writes the data output with NO sample suffix.
    const std::string dataFile = isW
        ? "../skim/rootfile/" + base + "_hist.root"
        : "../skim/rootfile/" + base + "_Data_hist.root";

    TFile *fData = TFile::Open(dataFile.c_str(), "READ");
    if (!fData || fData->IsZombie())
    {
        std::cerr << "[ERR] cannot open data file: " << dataFile << "\n";
        return;
    }

    // Signal MC: one DY file for Z; Wp + Wm for W (combined below).
    TFile *fMC = nullptr, *fWp = nullptr, *fWm = nullptr;
    if (!isW)
    {
        const std::string mcFile = "../skim/rootfile/" + base + "_DY_MC_hist.root";
        fMC = TFile::Open(mcFile.c_str(), "READ");
        if (!fMC || fMC->IsZombie())
        {
            std::cerr << "[ERR] cannot open MC file: " << mcFile << "\n";
            return;
        }
    }
    else
    {
        const std::string wpFile = "../skim/rootfile/" + base + "_Wp_hist.root";
        const std::string wmFile = "../skim/rootfile/" + base + "_Wm_hist.root";
        fWp = TFile::Open(wpFile.c_str(), "READ");
        fWm = TFile::Open(wmFile.c_str(), "READ");
        if (!fWp || fWp->IsZombie()) { std::cerr << "[ERR] cannot open MC file: " << wpFile << "\n"; return; }
        if (!fWm || fWm->IsZombie()) { std::cerr << "[ERR] cannot open MC file: " << wmFile << "\n"; return; }
    }

    // W+/W- relative normalization (k_s = A*sigma*L/N_gen from mc_norm.h; the
    // shape-normalization downstream cancels the absolute scale, only the
    // W+/W- mix survives). Falls back to 1.0 (+[WARN]) if ngen.root is absent.
    double kWp = 1.0, kWm = 1.0;
    if (isW)
    {
        const std::string flav = (ch == "Wmu") ? "mu" : "ele";
        kWp = pONorm::MCScale("Wp_" + flav);
        kWm = pONorm::MCScale("Wm_" + flav);
        std::cout << "[INFO] " << ch << " W+/W- mix: k_Wp=" << kWp
                  << ", k_Wm=" << kWm << " (ratio " << (kWm > 0 ? kWp / kWm : -1)
                  << ")\n";
    }

    const std::string outDir = "./plots/dataMC_" + ch;
    gSystem->mkdir(outDir.c_str(), kTRUE);

    struct Obs
    {
        const char *hist;
        std::string xTitle;
    };
    std::vector<Obs> obs;
    if (!isW)
    {
        obs = {
            {"h_Zpt", "p_{T}^{" + bosSym + "} [GeV]"},
            {"h_Zeta", "#eta_{" + bosSym + "}"},
            {"h_Zy", "y_{" + bosSym + "}"},
            {"h_Zphi", "#phi_{" + bosSym + "}"},
            {"h_lepPt", "p_{T}^{" + lepSym + "} [GeV]"},
            {"h_lepEta", "#eta_{" + lepSym + "}"},
            {"h_lepPhi", "#phi_{" + lepSym + "}"},
            {"hMass", "m_{" + bosSym + "} [GeV]"},
        };
    }
    else
    {
        obs = {
            {"h_lepPt", "p_{T}^{" + lepSym + "} [GeV]"},
            {"h_lepEta", "#eta_{" + lepSym + "}"},
            {"h_lepPhi", "#phi_{" + lepSym + "}"},
        };
    }

    // Labels / subtitles (Z cosmetics kept identical; W follows the same style)
    const std::string header = ch + (isW ? ": Data vs W MC" : ": Data vs DY MC");
    const std::string sub1 = isW
        ? "full W#rightarrow" + lepSym + "#nu selection"
        : "Z peak [60,120] GeV";
    const std::string sub2 = isW ? "shape-norm., bkg not subtracted"
                                 : "shape-normalized";
    const std::string mcLabel = isW ? "W signal MC (W^{+}+W^{-})" : "DY signal MC";

    PlotStyle ps;            // defaults from plotting_helper.C
    ps.yTitleOffset = 1.55;  // a touch more room for "Events (a.u.)"
    ps.boxY2 = 0.65;
    ps.headerY = 0.75;
    ps.headerX = 0.20;
    ps.titleSize = 0.04;

    int nOK = 0;
    for (const auto &o : obs)
    {
        TH1 *hD = (TH1 *)fData->Get(o.hist);

        TH1 *hM = nullptr;
        if (!isW)
        {
            hM = (TH1 *)fMC->Get(o.hist);
        }
        else
        {
            TH1 *hWp = (TH1 *)fWp->Get(o.hist);
            TH1 *hWm = (TH1 *)fWm->Get(o.hist);
            if (hWp && hWm)
            {
                // Clone before mutating (repo rule: never Scale/Add a
                // file-owned histogram in place).
                TH1 *sum = (TH1 *)hWp->Clone(Form("%s_Wsum", o.hist));
                sum->SetDirectory(nullptr);
                sum->Scale(kWp);
                sum->Add(hWm, kWm);
                hM = sum;
            }
        }

        if (!hD || !hM)
        {
            std::cerr << "[WARN] missing '" << o.hist << "' in data or MC file; skipping\n";
            if (isW)
                std::cerr << "       (W lepton histos need a re-skim with the updated skim.C)\n";
            continue;
        }

        const std::string out = outDir + "/" + o.hist;
        SaveDataMCRatio(hD, hM, out,
                        o.xTitle, "Events (a.u.)",
                        header,
                        sub1,
                        sub2,
                        ps,
                        /*normToData=*/true,
                        "Data",
                        mcLabel);
        std::cout << "[OK] wrote " << out << ".png/.pdf\n";
        ++nOK;
    }

    fData->Close();
    if (fMC) fMC->Close();
    if (fWp) fWp->Close();
    if (fWm) fWm->Close();
    std::cout << "[DONE] " << ch << ": wrote " << nOK << " plot(s) under " << outDir << "/\n";
}
