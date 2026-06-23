// correction/recoil_raw_ele.C
//
// Raw look at the Z->ee hadronic-recoil components (u_par, u_perp) used for
// the MET recoil correction: data vs DY signal MC, inclusive and in q_T slices.
// Electron-channel twin of correction/recoil_raw.C (Z->mumu).
// NO FIT yet -- the point is to judge the available statistics and pick a
// sensible q_T binning before the double-Gaussian fits go in (Sec 6.1 of the AN).
//
// Inputs : ../skim/rootfile/ZToEE_pO2025_Data_hist.root   (data)
//          ../skim/rootfile/ZToEE_pO2025_DY_MC_hist.root  (DY signal MC)
//          histos: h_uPar, h_uPerp (1D inclusive); h_uPar_qT, h_uPerp_qT (2D vs q_T)
//          (these are produced by skim_Zee, mirroring skim_Zmm -- re-run the
//           Zee skim after the recoil histos were added or they'll be absent.)
// Output : ./plots/recoil_Zee/  (i.e. correction/plots/recoil_Zee/)
//          + a printed table of entries per q_T slice.
//
// The q_T slicing is done by projecting the 2D histograms, so you can change
// kQtEdges below and re-run WITHOUT re-running the skim.
//
// Run from correction/:  root -l -q 'recoil_raw_ele.C+'

#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"

#include <iostream>
#include <string>
#include <vector>

#include "../plotting/plotting_helper.C"

namespace
{
// q_T slice edges [GeV] for projecting the 2D recoil histograms. Keep them on
// the 2 GeV/c grid of the skim histograms. Coarse to start (low pO statistics,
// and the ee channel has fewer events than mumu); refine once you see the entry
// table. Edit freely -- these are projections.
const std::vector<double> kQtEdges = {0, 10, 20, 30, 50, 140};
} // namespace

// ProjectionY (the u distribution) of a 2D recoil histo for q_T in [qlo, qhi).
static TH1D *projU(TH2 *h2, double qlo, double qhi, const char *name)
{
    if (!h2)
        return nullptr;
    const int bxlo = h2->GetXaxis()->FindBin(qlo + 1e-6);
    const int bxhi = h2->GetXaxis()->FindBin(qhi - 1e-6);
    TH1D *p = h2->ProjectionY(name, bxlo, bxhi);
    p->SetDirectory(nullptr);
    return p;
}

void recoil_raw_ele()
{
    const std::string dir = "../skim/rootfile/";
    TFile *fData = TFile::Open((dir + "ZToEE_pO2025_Data_hist.root").c_str(), "READ");
    if (!fData || fData->IsZombie())
    {
        std::cerr << "[ERR] cannot open " << dir << "ZToEE_pO2025_Data_hist.root\n";
        return;
    }
    TFile *fMC = TFile::Open((dir + "ZToEE_pO2025_DY_MC_hist.root").c_str(), "READ");
    if (!fMC || fMC->IsZombie())
    {
        std::cerr << "[ERR] cannot open " << dir << "ZToEE_pO2025_DY_MC_hist.root\n";
        return;
    }

    const std::string outDir = "./plots/recoil_Zee";
    gSystem->mkdir(outDir.c_str(), kTRUE);

    PlotStyle ps;            // defaults from plotting_helper.C
    ps.yTitleOffset = 1.55;  // a touch more room for "Events (a.u.)"
    ps.boxY2 = 0.65;
    ps.headerY = 0.75;
    ps.headerX = 0.20;
    ps.titleSize = 0.04;

    // ---- inclusive u_par / u_perp (data vs MC, shape-normalized, with ratio) ----
    struct Incl { const char *name; std::string xt; };
    for (const Incl &o : std::vector<Incl>{{"h_uPar", "u_{#parallel} [GeV]"},
                                           {"h_uPerp", "u_{#perp} [GeV]"}})
    {
        TH1 *hD = (TH1 *)fData->Get(o.name);
        TH1 *hM = (TH1 *)fMC->Get(o.name);
        if (!hD || !hM)
        {
            std::cerr << "[WARN] missing " << o.name << " in data or MC; skipping\n";
            continue;
        }
        SaveDataMCRatio(hD, hM, outDir + "/" + o.name,
                        o.xt, "Events (a.u.)",
                        "Z#rightarrowee recoil (inclusive)",
                        "Data vs DY MC", "shape-normalized",
                        ps, true, "Data", "DY MC");
    }

    // ---- u_par / u_perp in q_T slices + entry table ----
    struct Comp { const char *h2; const char *tag; std::string xt; };
    std::cout << "\n=== recoil entries per q_T slice (use this to judge binning) ===\n";
    for (const Comp &c : std::vector<Comp>{{"h_uPar_qT", "uPar", "u_{#parallel} [GeV]"},
                                           {"h_uPerp_qT", "uPerp", "u_{#perp} [GeV]"}})
    {
        TH2 *d2 = (TH2 *)fData->Get(c.h2);
        TH2 *m2 = (TH2 *)fMC->Get(c.h2);
        if (!d2 || !m2)
        {
            std::cerr << "[WARN] missing " << c.h2 << " in data or MC; skipping\n";
            continue;
        }
        std::cout << "-- " << c.tag << " --\n";
        for (size_t i = 0; i + 1 < kQtEdges.size(); ++i)
        {
            const double qlo = kQtEdges[i], qhi = kQtEdges[i + 1];
            TH1D *pd = projU(d2, qlo, qhi, Form("d_%s_q%.0f_%.0f", c.tag, qlo, qhi));
            TH1D *pm = projU(m2, qlo, qhi, Form("m_%s_q%.0f_%.0f", c.tag, qlo, qhi));
            if (!pd || !pm)
                continue;

            std::cout << "   q_T [" << qlo << ", " << qhi << ")  data = "
                      << pd->Integral() << "   MC = " << pm->Integral() << "\n";

            if (pd->Integral() > 0 || pm->Integral() > 0)
                SaveDataMCRatio(pd, pm,
                                Form("%s/%s_qT%.0f_%.0f", outDir.c_str(), c.tag, qlo, qhi),
                                c.xt, "Events (a.u.)",
                                Form("Z#rightarrowee: q_{T} #in [%.0f, %.0f) GeV", qlo, qhi),
                                "Data vs DY MC", "shape-normalized",
                                ps, true, "Data", "DY MC");
            delete pd;
            delete pm;
        }
    }

    fData->Close();
    fMC->Close();
    std::cout << "[DONE] recoil raw plots -> " << outDir << "/\n";
}
