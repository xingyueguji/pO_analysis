#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TSystem.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// =============================================================================
// fiducial_yields.C -- the ACCEPTANCE-CORRECTED yields of one simultaneous fit
// (2026-09-22):
//
//     Y_i = r_i x sigma_gen-fid,i     per (charge, rapidity bin), both binnings
//
// with the full r covariance (total and stat) from the fit's h_cov_poi[_FB][_stat],
// written under EXACTLY the names analysis/charge_asym.C and analysis/FBratio.C
// read -- h_yield_W{p,m}_y{i}[_FB] + h_cov_yield[_FB][_stat] -- so those two
// macros, unchanged, return the FIDUCIAL charge asymmetry and forward/backward
// ratios (Y in nb; both observables are ratios, so the unit drops out).
//
// WHY. The fitted yields r x S that the fork's extractor writes are RECO-level
// counts: S = L sigma_gen (A x eps)_MC. In the charge asymmetry the A x eps of
// W+ and W- leptons in the SAME eta bin nearly cancels, but R_FB divides two
// DIFFERENT |eta_lab| regions (eta_lab = eta_CM + 0.3466), so there the
// detector acceptance does not cancel at all: the electron ECAL crack
// (|eta_SC| 1.44-1.57) sits in the FORWARD bin at |eta_CM| ~ 1.2 and in the
// BACKWARD one at ~ 1.9. Measured on the grand fit's own r's (i.e. IDENTICAL
// physics for both flavours): the count-based R_FB of the mu and e yields
// differ by up to 60% (0.67 vs 1.11, 1.38 vs 0.84), and A_ch by up to 0.015.
// sigma_gen is the SAME pooled mu+e gen fiducial for every fit (skim/gen_xsec.C:
// bare lepton pT > 25 GeV, |eta_lab| < 2.4), so these yields sit in ONE
// fiducial volume: comparable between the mu-only, e-only and grand fits, and
// the level at which the theory curves are defined.
//
// POLICY (user, 2026-09-22): EVERY observable -- sigma, A_ch, R_FB, for the
// grand fit and the per-flavour ones -- is built from r x sigma_gen, NEVER
// from the raw (count-based) fitted yields: there is no dedicated efficiency
// or acceptance correction, and r x sigma_gen is what applies it (from MC).
// The count-based yields stay only as the fit's record and in diagnostics.
//
// Inputs: the fit's <tag>_W_yields.csv (r per bin, both binnings) and
//         <tag>_fitted_yields.root (h_cov_poi[_FB][_stat], the 25x25 POI
//         covariance of extract_pO_simfit.C -- read by AXIS LABEL; older
//         extractions without it fall back to h_cov_yield[_FB][_stat]/(S S),
//         then to the diagonal, and say so), and skim/gen_xsec.root
//         (h_gen_sig_W{p,m}[_FB]).
// Output: one ROOT file per fit, e.g. ../skim/rootfile/fidyields_simfit_mu_<disc>.root
// Driven by analysis/run_observables.sh; by hand:
//   root -l -b -q 'fiducial_yields.C+("<tag>_W_yields.csv","<tag>_fitted_yields.root","out.root")'
// =============================================================================

namespace
{

// one binning of the CSV, per (charge, y bin) in the [Wp_y0..11, Wm_y0..11]
// order: r (col 4), rErr (col 5), the prefit S (col 6) and, when the 19th
// column exists, rErr_stat (col 18; -1 otherwise); true when all 24 found
struct CsvR
{
    double r[24], e[24], S[24], es[24];
};
bool ReadR(const char *csv, const char *binning, CsvR &out)
{
    for (int k = 0; k < 24; ++k) { out.r[k] = out.e[k] = out.S[k] = 0.0; out.es[k] = -1.0; }
    std::ifstream in(csv);
    if (!in) { std::cerr << "[ERROR] cannot open " << csv << "\n"; return false; }
    std::string line;
    std::getline(in, line); // header
    int found = 0;
    while (std::getline(in, line))
    {
        std::stringstream ss(line);
        std::string f;
        std::vector<std::string> c;
        while (std::getline(ss, f, ',')) c.push_back(f);
        if (c.size() < 7 || c[2] != binning) continue;
        const int off = (c[1] == "Wp") ? 0 : (c[1] == "Wm") ? 12 : -1;
        const int iy = std::atoi(c[3].c_str());
        if (off < 0 || iy < 0 || iy >= 12) continue;
        out.r[off + iy] = std::atof(c[4].c_str());
        out.e[off + iy] = std::atof(c[5].c_str());
        out.S[off + iy] = std::atof(c[6].c_str());
        if (c.size() >= 19) out.es[off + iy] = std::atof(c[18].c_str());
        ++found;
    }
    if (found != 24) std::cerr << "[ERROR] " << found << "/24 '" << binning << "' rows in " << csv << "\n";
    return found == 24;
}

// the 24x24 W block of a 25x25 h_cov_poi*, by axis label (never by position)
bool ReadCovR(TFile *f, const char *name, double C[24][24])
{
    TH2D *h = (TH2D *)f->Get(name);
    if (!h || h->GetNbinsX() != 25 || h->GetNbinsY() != 25) return false;
    for (int a = 0; a < 24; ++a)
    {
        const TString want = TString::Format("r_%s_y%d", a < 12 ? "Wp" : "Wm", a % 12);
        if (TString(h->GetXaxis()->GetBinLabel(a + 1)) != want)
        {
            std::cerr << "[ERROR] " << name << " bin " << a + 1 << " is '" << h->GetXaxis()->GetBinLabel(a + 1)
                      << "', expected '" << want << "' -- refusing to guess the order\n";
            return false;
        }
    }
    for (int a = 0; a < 24; ++a)
        for (int b = 0; b < 24; ++b) C[a][b] = h->GetBinContent(a + 1, b + 1);
    return true;
}

// FALLBACK for extractions older than h_cov_poi (2026-09-15, e.g. the Aug-6
// met tree): the r covariance recovered from the 24x24 yield covariance,
// cov(r_a, r_b) = h_cov_yield[a][b] / (S_a S_b) -- exact, h_cov_yield being
// S_a S_b cov(r_a, r_b) by construction (the extractor, and
// xsec_fiducial.C::loadCombIngredients, use the same identity)
bool ReadCovFromYield(TFile *f, const char *name, const double S[24], double C[24][24])
{
    TH2D *h = (TH2D *)f->Get(name);
    if (!h || h->GetNbinsX() != 24 || h->GetNbinsY() != 24) return false;
    for (int a = 0; a < 24; ++a)
    {
        const TString want = TString::Format("%s_y%d", a < 12 ? "Wp" : "Wm", a % 12);
        const TString got = h->GetXaxis()->GetBinLabel(a + 1);
        if (got != "" && got != want)
        {
            std::cerr << "[ERROR] " << name << " bin " << a + 1 << " is '" << got << "', expected '" << want << "'\n";
            return false;
        }
        if (S[a] <= 0) { std::cerr << "[ERROR] no prefit S for bin " << a << " -- cannot convert " << name << "\n"; return false; }
    }
    for (int a = 0; a < 24; ++a)
        for (int b = 0; b < 24; ++b) C[a][b] = h->GetBinContent(a + 1, b + 1) / (S[a] * S[b]);
    return true;
}

} // namespace

void fiducial_yields(const char *yieldsCsv,        // <tag>_W_yields.csv of the fit
                     const char *fittedYieldsRoot, // <tag>_fitted_yields.root of the fit
                     const char *outFile,
                     const char *genFile = "../skim/rootfile/gen_xsec.root")
{
    TFile *fG = TFile::Open(genFile, "READ");
    TFile *fY = TFile::Open(fittedYieldsRoot, "READ");
    if (!fG || fG->IsZombie() || !fY || fY->IsZombie())
    {
        std::cerr << "[ERROR] fiducial_yields: cannot open " << genFile << " / " << fittedYieldsRoot << "\n";
        return;
    }
    gSystem->mkdir(gSystem->DirName(outFile), kTRUE);
    TFile *fo = TFile::Open(outFile, "RECREATE");
    if (!fo || fo->IsZombie()) { std::cerr << "[ERROR] cannot create " << outFile << "\n"; return; }

    int nWritten = 0;
    const char *bins[2] = {"lab", "fb"};
    for (int ib = 0; ib < 2; ++ib)
    {
        const bool fb = (ib == 1);
        const TString sfx = fb ? "_FB" : "";
        // gen fiducial sigma of this binning, [Wp_y0..11, Wm_y0..11]
        TH1D *hp = (TH1D *)fG->Get(TString("h_gen_sig_Wp") + sfx);
        TH1D *hm = (TH1D *)fG->Get(TString("h_gen_sig_Wm") + sfx);
        if (!hp || !hm || hp->GetNbinsX() != 12 || hm->GetNbinsX() != 12)
        {
            std::cerr << "[ERROR] no 12-bin h_gen_sig_W{p,m}" << sfx << " in " << genFile << "\n";
            continue;
        }
        double G[24], C[24][24], Cs[24][24];
        for (int i = 0; i < 12; ++i) { G[i] = hp->GetBinContent(i + 1); G[12 + i] = hm->GetBinContent(i + 1); }
        CsvR in;
        if (!ReadR(yieldsCsv, bins[ib], in)) continue;
        const double *r = in.r;
        // the r covariance, best source first: the 25x25 POI covariance
        // (2026-09-15+) -> the 24x24 yield covariance / (S S) (older
        // extractions) -> the diagonal (r errors from the CSV), each with its
        // stat twin when the fit has one; the source is printed
        const TString poi = TString("h_cov_poi") + sfx, yld = TString("h_cov_yield") + sfx;
        TString src;
        if (ReadCovR(fY, poi, C)) src = poi;
        else if (ReadCovFromYield(fY, yld, in.S, C)) src = yld + " / (S S)";
        else
        {
            for (int a = 0; a < 24; ++a)
                for (int b = 0; b < 24; ++b) C[a][b] = (a == b) ? in.e[a] * in.e[a] : 0.0;
            src = "the DIAGONAL rErr^2 (no covariance in the file -- bin correlations lost)";
        }
        bool haveStat = ReadCovR(fY, poi + "_stat", Cs) || ReadCovFromYield(fY, yld + "_stat", in.S, Cs);
        if (!haveStat && in.es[0] >= 0)
        {
            for (int a = 0; a < 24; ++a)
                for (int b = 0; b < 24; ++b) Cs[a][b] = (a == b) ? in.es[a] * in.es[a] : 0.0;
            haveStat = true;
            std::cerr << "[WARN] " << bins[ib] << ": stat covariance from the DIAGONAL rErr_stat only\n";
        }
        std::cout << "[fiducial_yields] " << bins[ib] << ": r covariance from " << src
                  << (haveStat ? "" : "; NO stat component in this extraction (single total bars downstream)") << "\n";

        fo->cd();
        for (int k = 0; k < 24; ++k)
        {
            const TString nm = TString::Format("h_yield_%s_y%d%s", k < 12 ? "Wp" : "Wm", k % 12, sfx.Data());
            TH1D *h = new TH1D(nm, nm + " (fiducial: r x sigma_gen, nb)", 1, 0.0, 1.0);
            h->Sumw2();
            h->SetBinContent(1, r[k] * G[k]);
            h->SetBinError(1, std::sqrt(std::max(0.0, C[k][k])) * G[k]);
            h->Write();
            ++nWritten;
        }
        for (int pass = 0; pass < 2; ++pass)
        {
            if (pass == 1 && !haveStat) continue;
            const TString nm = TString("h_cov_yield") + sfx + (pass ? "_stat" : "");
            TH2D *hc = new TH2D(nm, nm + " (fiducial yields, nb^2);index;index", 24, 0, 24, 24, 0, 24);
            for (int a = 0; a < 24; ++a)
            {
                const TString lab = TString::Format("%s_y%d", a < 12 ? "Wp" : "Wm", a % 12);
                hc->GetXaxis()->SetBinLabel(a + 1, lab);
                hc->GetYaxis()->SetBinLabel(a + 1, lab);
                for (int b = 0; b < 24; ++b)
                    hc->SetBinContent(a + 1, b + 1, G[a] * G[b] * (pass ? Cs[a][b] : C[a][b]));
            }
            hc->Write();
        }
    }
    fo->Close();
    fG->Close();
    fY->Close();
    if (nWritten != 48)
    {
        std::cerr << "[ERROR] fiducial_yields: wrote " << nWritten << "/48 yields -> " << outFile << " is incomplete\n";
        gSystem->Unlink(outFile); // never leave a half file for the next macro to read
        return;
    }
    std::cout << "[OK] fiducial yields (r x sigma_gen, lab + fb, total + stat covariance) -> " << outFile << "\n";
}
