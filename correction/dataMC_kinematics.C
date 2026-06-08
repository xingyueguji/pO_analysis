// correction/dataMC_kinematics.C
//
// Data vs signal-MC kinematic check for the Z channels (boson pT / eta / phi,
// lepton pT / eta / phi, and the dilepton mass), with a ratio pad.
//
// This is the "Data/MC pT check" used to validate (and later derive) the boson
// pT reweighting. The comparison is shape-only: the signal MC (DY) is scaled to
// the data integral, so the ratio pad probes shape differences (backgrounds are
// negligible inside the Z peak, so signal-vs-data is the right comparison here).
//
// Inputs  : ../skim/rootfile/<base>_Data_hist.root   (data)
//           ../skim/rootfile/<base>_DY_MC_hist.root  (DY signal MC)
//           where <base> = ZToMuMu_pO2025 (Zmm) or ZToEE_pO2025 (Zee).
// Outputs : ./plots/dataMC_<channel>/<hist>.png/.pdf   (i.e. correction/plots/)
//
// Histograms (produced by skim.C, filled for iso-selected OS pairs in the Z
// peak [60,120] GeV): hMass, h_Zpt, h_Zeta, h_Zphi, h_lepPt, h_lepEta, h_lepPhi.
//
// Run from this directory (correction/):
//   root -l -q 'dataMC_kinematics.C+("Zmm")'
//   root -l -q 'dataMC_kinematics.C+("Zee")'

#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"

#include <iostream>
#include <string>
#include <vector>

#include "../plotting/plotting_helper.C"

void dataMC_kinematics(const char *channel = "Zmm")
{
    const std::string ch = channel;

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
    else
    {
        std::cerr << "[ERR] channel must be \"Zmm\" or \"Zee\" (got \"" << ch << "\")\n";
        return;
    }

    const std::string dataFile = "../skim/rootfile/" + base + "_Data_hist.root";
    const std::string mcFile = "../skim/rootfile/" + base + "_DY_MC_hist.root";

    TFile *fData = TFile::Open(dataFile.c_str(), "READ");
    if (!fData || fData->IsZombie())
    {
        std::cerr << "[ERR] cannot open data file: " << dataFile << "\n";
        return;
    }
    TFile *fMC = TFile::Open(mcFile.c_str(), "READ");
    if (!fMC || fMC->IsZombie())
    {
        std::cerr << "[ERR] cannot open MC file: " << mcFile << "\n";
        return;
    }

    const std::string outDir = "./plots/dataMC_" + ch;
    gSystem->mkdir(outDir.c_str(), kTRUE);

    struct Obs
    {
        const char *hist;
        std::string xTitle;
    };
    const std::vector<Obs> obs = {
        {"h_Zpt", "p_{T}^{" + bosSym + "} [GeV]"},
        {"h_Zeta", "#eta_{" + bosSym + "}"},
        {"h_Zphi", "#phi_{" + bosSym + "}"},
        {"h_lepPt", "p_{T}^{" + lepSym + "} [GeV]"},
        {"h_lepEta", "#eta_{" + lepSym + "}"},
        {"h_lepPhi", "#phi_{" + lepSym + "}"},
        {"hMass", "m_{" + bosSym + "} [GeV]"},
    };

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
        TH1 *hM = (TH1 *)fMC->Get(o.hist);
        if (!hD || !hM)
        {
            std::cerr << "[WARN] missing '" << o.hist << "' in data or MC file; skipping\n";
            continue;
        }

        const std::string out = outDir + "/" + o.hist;
        SaveDataMCRatio(hD, hM, out,
                        o.xTitle, "Events (a.u.)",
                        ch + ": Data vs DY MC",
                        "Z peak [60,120] GeV",
                        "shape-normalized",
                        ps,
                        /*normToData=*/true,
                        "Data",
                        "DY signal MC");
        std::cout << "[OK] wrote " << out << ".png/.pdf\n";
        ++nOK;
    }

    fData->Close();
    fMC->Close();
    std::cout << "[DONE] " << ch << ": wrote " << nOK << " plot(s) under " << outDir << "/\n";
}
