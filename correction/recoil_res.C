// correction/recoil_res.C
//
// The MET-RESOLUTION factor of the MC, measured on the Z recoil (2026-10-20).
// It feeds the metRes shape nuisance of the MET-discriminant fit (met_mt40):
// skim/met_res.h reads the result and fills <h>_metRes{Mu,Ele}Up/Down twins
// with the MC MET smeared to the data resolution.
//
// Measurement, per flavour (Z->mumu for the muon channel, Z->ee for the
// electron one -- the two are measured on different events and enter as two
// independent nuisances):
//   u_perp = the hadronic recoil component perpendicular to the dilepton q_T
//            (skim h_uPerp: u = -MET - q_T, iso-selected OS pairs in [60,120]).
//            It carries no q_T response, so its width is the recoil resolution.
//   sigma  = RMS of h_uPerp over its full range [-200, 200] GeV (2 GeV bins).
//   k      = sigma_data / sigma_MC (DY MC = the Z peak's MC; the non-DY
//            backgrounds under the peak are < 0.5% and left out).
//   sigma_add = sigma_MC x sqrt(k^2 - 1): the Gaussian smearing that, added
//            to EACH transverse MET component of the MC, brings the MC u_perp
//            width to the data one (0 if k <= 1).
// The u_par mean (the recoil RESPONSE) is printed as a record only: it is not
// varied (data - MC = +1.9 +- 1.1 GeV mu, +0.9 +- 1.7 GeV e on 2026-10-20).
//
// Inputs : ../skim/rootfile/ZTo{MuMu,EE}_pO2025_{Data,DY_MC}_hist.root
// Output : rootfile/recoil_res.root  -- TParameter<double> per flavour f in
//          {mu, ele}: k_<f>, kErr_<f>, sigmaMC_<f>, sigmaData_<f>, sigmaAdd_<f>
//          (GeV) -- read by skim/met_res.h (missing => the W MC skim is FATAL)
//          plots/recoil_res/uperp_{mu,ele}{,_smeared}  data vs MC (vs the MC
//          convolved with the Gaussian sigma_add: the closure of the smearing)
//
// Run via ./run_recoil_res.sh -> logs/recoil_res.log (the record). Run it
// BEFORE the W MC skim (and after any Z skim).

#include "../plotting/plotting_helper.C"

#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TParameter.h"
#include "TSystem.h"

#include <cmath>
#include <iostream>
#include <string>

namespace
{
struct RecoilRes
{
    double sigD = 0, sigDErr = 0, sigM = 0, sigMErr = 0, k = 0, kErr = 0, add = 0;
    double parD = 0, parDErr = 0, parM = 0, nData = 0;
};

// MC u_perp convolved with a Gaussian of width sAdd (bin-centre sum on the
// histogram's own grid): what the smeared MC u_perp looks like
TH1D *convolve(const TH1 *h, double sAdd, const char *name)
{
    TH1D *c = (TH1D *)h->Clone(name);
    c->SetDirectory(nullptr);
    c->Reset();
    const int nb = h->GetNbinsX();
    for (int i = 1; i <= nb; ++i)
    {
        const double v = h->GetBinContent(i);
        if (v == 0.0) continue;
        const double x = h->GetBinCenter(i);
        if (sAdd <= 0.0) { c->AddBinContent(i, v); continue; }
        for (int j = 1; j <= nb; ++j)
        {
            const double lo = h->GetXaxis()->GetBinLowEdge(j), hi = h->GetXaxis()->GetBinUpEdge(j);
            const double p = 0.5 * (TMath::Erf((hi - x) / (std::sqrt(2.0) * sAdd)) -
                                    TMath::Erf((lo - x) / (std::sqrt(2.0) * sAdd)));
            c->AddBinContent(j, v * p);
        }
    }
    for (int j = 1; j <= nb; ++j) c->SetBinError(j, 0.0);
    c->ResetStats(); // AddBinContent does not update the stored moments (GetRMS)
    return c;
}

bool measure(const std::string &flav, const std::string &prefix, const std::string &outDir,
             RecoilRes &r)
{
    const std::string dir = "../skim/rootfile/";
    TFile *fd = TFile::Open((dir + prefix + "_Data_hist.root").c_str(), "READ");
    TFile *fm = TFile::Open((dir + prefix + "_DY_MC_hist.root").c_str(), "READ");
    if (!fd || fd->IsZombie() || !fm || fm->IsZombie())
    {
        std::cerr << "[ERR] recoil_res " << flav << ": cannot open " << dir << prefix << "_{Data,DY_MC}_hist.root\n";
        return false;
    }
    TH1 *dPerp = (TH1 *)fd->Get("h_uPerp"), *mPerp = (TH1 *)fm->Get("h_uPerp");
    TH1 *dPar = (TH1 *)fd->Get("h_uPar"), *mPar = (TH1 *)fm->Get("h_uPar");
    if (!dPerp || !mPerp || !dPar || !mPar)
    {
        std::cerr << "[ERR] recoil_res " << flav << ": h_uPerp / h_uPar missing (Z skim without the PF tree?)\n";
        return false;
    }
    r.nData = dPerp->Integral(1, dPerp->GetNbinsX());
    if (r.nData <= 0.0 || mPerp->Integral() <= 0.0)
    {
        std::cerr << "[ERR] recoil_res " << flav << ": empty h_uPerp\n";
        return false;
    }
    const double ovD = dPerp->GetBinContent(0) + dPerp->GetBinContent(dPerp->GetNbinsX() + 1);
    const double ovM = mPerp->GetBinContent(0) + mPerp->GetBinContent(mPerp->GetNbinsX() + 1);
    r.sigD = dPerp->GetRMS();  r.sigDErr = dPerp->GetRMSError();
    r.sigM = mPerp->GetRMS();  r.sigMErr = mPerp->GetRMSError();
    r.k = r.sigD / r.sigM;
    r.kErr = r.k * std::sqrt(std::pow(r.sigDErr / r.sigD, 2) + std::pow(r.sigMErr / r.sigM, 2));
    r.add = (r.k > 1.0) ? r.sigM * std::sqrt(r.k * r.k - 1.0) : 0.0;
    r.parD = dPar->GetMean(); r.parDErr = dPar->GetMeanError(); r.parM = mPar->GetMean();

    std::cout << "\n=== MET resolution from the Z recoil: " << flav << " (" << prefix << ") ===\n"
              << Form("  u_perp RMS  data %.3f +- %.3f GeV (%.0f pairs, %.0f outside +-200)   DY MC %.3f +- %.3f GeV"
                      " (%.2f%% outside)\n",
                      r.sigD, r.sigDErr, r.nData, ovD, r.sigM, r.sigMErr,
                      mPerp->Integral() > 0 ? 100.0 * ovM / mPerp->Integral(0, mPerp->GetNbinsX() + 1) : 0.0)
              << Form("  k = sigma_data / sigma_MC = %.4f +- %.4f   (%.1f sigma from 1)\n", r.k, r.kErr,
                      r.kErr > 0 ? (r.k - 1.0) / r.kErr : 0.0)
              << Form("  sigma_add = sigma_MC sqrt(k^2-1) = %.3f GeV per MET component%s\n", r.add,
                      r.k > 1.0 ? "" : "  (k <= 1: no smearing possible, the nuisance is flat)")
              << Form("  u_par mean (response, NOT varied): data %.2f +- %.2f GeV, MC %.2f GeV -> data - MC = %+.2f +- %.2f\n",
                      r.parD, r.parDErr, r.parM, r.parD - r.parM, r.parDErr);

    // closure: data vs MC, and data vs MC (x) Gaussian(sigma_add)
    TH1D *mSm = convolve(mPerp, r.add, Form("uperp_smeared_%s", flav.c_str()));
    const double rmsSm = mSm->GetRMS();
    std::cout << Form("  closure: RMS of the MC convolved with Gauss(sigma_add) = %.3f GeV (data %.3f)\n", rmsSm, r.sigD);

    PlotStyle ps;
    ps.logy = false;
    ps.yTitleOffset = 1.55;
    ps.headerX = 0.62;
    ps.headerY = 0.66;
    ps.titleSize = 0.04;
    ps.xRangeLo = -60.0;
    ps.xRangeHi = 60.0;
    const std::string chan = (flav == "mu") ? "Z #rightarrow #mu#mu" : "Z #rightarrow ee";
    SaveDataMCRatio(dPerp, mPerp, outDir + "/uperp_" + flav, "u_{#perp} (GeV)", "Events (shape-normalized)",
                    chan, Form("RMS: data %.2f, MC %.2f GeV", r.sigD, r.sigM), Form("k = %.3f #pm %.3f", r.k, r.kErr),
                    ps, true, "Data", "DY MC");
    SaveDataMCRatio(dPerp, mSm, outDir + "/uperp_" + flav + "_smeared", "u_{#perp} (GeV)", "Events (shape-normalized)",
                    chan, Form("MC #otimes Gauss(%.2f GeV)", r.add),
                    Form("RMS: data %.2f, MC %.2f GeV", r.sigD, rmsSm), ps, true, "Data", "DY MC smeared");
    delete mSm;
    fd->Close(); fm->Close();
    delete fd; delete fm;
    return true;
}
} // namespace

void recoil_res()
{
    const std::string outDir = "./plots/recoil_res";
    gSystem->mkdir(outDir.c_str(), kTRUE);
    gSystem->mkdir("./rootfile", kTRUE);

    RecoilRes mu, ele;
    const bool okMu = measure("mu", "ZToMuMu_pO2025", outDir, mu);
    const bool okEle = measure("ele", "ZToEE_pO2025", outDir, ele);
    if (!okMu || !okEle)
    {
        std::cerr << "[ERR] recoil_res: incomplete -> rootfile/recoil_res.root NOT written\n";
        return;
    }

    TFile *fo = TFile::Open("./rootfile/recoil_res.root", "RECREATE");
    for (const auto &p : {std::make_pair(std::string("mu"), &mu), std::make_pair(std::string("ele"), &ele)})
    {
        const RecoilRes &r = *p.second;
        const std::string f = p.first;
        TParameter<double>(("k_" + f).c_str(), r.k).Write();
        TParameter<double>(("kErr_" + f).c_str(), r.kErr).Write();
        TParameter<double>(("sigmaMC_" + f).c_str(), r.sigM).Write();
        TParameter<double>(("sigmaData_" + f).c_str(), r.sigD).Write();
        TParameter<double>(("sigmaAdd_" + f).c_str(), r.add).Write();
    }
    fo->Close();
    delete fo;

    std::cout << "\n[RESULT] metRes nuisance inputs (+1 sigma = MC MET smeared by sigma_add per component):\n"
              << Form("  mu : k = %.4f +- %.4f  sigma_add = %.3f GeV\n", mu.k, mu.kErr, mu.add)
              << Form("  ele: k = %.4f +- %.4f  sigma_add = %.3f GeV\n", ele.k, ele.kErr, ele.add)
              << "[DONE] rootfile/recoil_res.root, plots in " << outDir << "/\n";
}
