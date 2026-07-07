// correction/isolation_mu_tight.C
//
// PURE muon isolation study under the tight muon ID -- the muon twin of the
// electron continuous-relIso scan in isolation_ele.C, in the same presentation.
//
// This is a fresh, minimal macro: the original isolation.C (kept untouched)
// also scanned PF-candidate cones (0.2/0.3/0.4), several iso types and a pT
// scan; here there is exactly ONE ID (muIDTight) and ONE isolation variable,
// the branch-based PF relIso with the Delta-beta pileup correction (MuonPOG
// convention), identical to the skim's RelIsoPF (skim/skim_common.h):
//
//     relIso = (muPFChIso + max(0, muPFNeuIso + muPFPhoIso - 0.5*muPFPUIso)) / muPt
//
// The uncorrected (ch+neu+pho)/pT scan is stored alongside purely as a
// reference, to quantify how little the Delta-beta term changes at pO pileup.
//
// Samples (all from the SAME data file, as in isolation_ele.C):
//   signal   = OS dimuon Z-window pairs (80 < m < 100), both legs tight ID
//   bkg SS   = same-sign pairs in the same window (~empty for muons)
//   bkg QCD  = single tight-ID muon events with PF MET < 5 GeV
//
// Output: rootfile/IsoStudyOutputs_muon_tight.root
//   contEffOS_<tag>, contEffSS_<tag>, contEffMET_<tag>   (eff vs relIso cut)
//   contPassOS_<tag>, contTotOS_<tag>, contPassSS_<tag>, contPassMET_<tag>
//   contRocSS_<tag>, contRocMET_<tag>  + TParameter contAucSS_/contAucMET_<tag>
//   with <tag> = dbeta (primary) and nodbeta (reference)
// -- the same object names/shapes plot_iso_summary.C consumes for the electron.
//
// Run from correction/:
//   root -l -q 'isolation_mu_tight.C+("/Users/zhenghuang/pO_2026_May_26/Data_May_26.root")'

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TGraph.h>
#include <TParameter.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <algorithm>

namespace
{

struct MuPairInfo
{
    int i;
    int j;
    double mass;
};

TVector2 ComputePFMET_mu(std::vector<int> *pfId,
                         std::vector<float> *pfPt,
                         std::vector<float> *pfPhi)
{
    (void)pfId;
    double sumPx = 0.0, sumPy = 0.0;
    const int nPF = (pfPt ? (int)pfPt->size() : 0);
    for (int i = 0; i < nPF; ++i)
    {
        sumPx += pfPt->at(i) * std::cos(pfPhi->at(i));
        sumPy += pfPt->at(i) * std::sin(pfPhi->at(i));
    }
    return TVector2(-sumPx, -sumPy);
}

double rocAUC_mu(const std::vector<double> &x, const std::vector<double> &y)
{
    std::vector<std::pair<double, double>> xy;
    xy.reserve(x.size());
    for (size_t i = 0; i < x.size(); ++i)
        xy.push_back({x[i], y[i]});
    std::sort(xy.begin(), xy.end(),
              [](const std::pair<double, double> &a, const std::pair<double, double> &b)
              { return a.first < b.first; });
    double auc = 0.0;
    for (size_t i = 1; i < xy.size(); ++i)
        auc += 0.5 * (xy[i - 1].second + xy[i].second) * (xy[i].first - xy[i - 1].first);
    return auc;
}

} // namespace

void isolation_mu_tight(const char *fname = "/Users/zhenghuang/pO_2026_May_26/Data_May_26.root")
{
    TFile *f = TFile::Open(fname, "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "Cannot open: " << fname << "\n";
        return;
    }

    gSystem->mkdir("rootfile", kTRUE);

    TTree *tMu = (TTree *)f->Get("ggHiNtuplizer/EventTree");
    TTree *tPF = (TTree *)f->Get("particleFlowAnalyser/pftree");
    if (!tMu || !tPF)
    {
        std::cerr << "Missing ggHiNtuplizer/EventTree or particleFlowAnalyser/pftree.\n";
        return;
    }

    // ---------- muon branches (fast skim: disable all, enable used) ----------
    Int_t nMu = 0;
    std::vector<float> *muPt = nullptr, *muEta = nullptr, *muPhi = nullptr;
    std::vector<int>   *muCharge = nullptr, *muIDTight = nullptr;
    std::vector<float> *muPFChIso = nullptr, *muPFNeuIso = nullptr, *muPFPhoIso = nullptr;
    std::vector<float> *muPFPUIso = nullptr;

    tMu->SetBranchStatus("*", 0);
    tMu->SetBranchStatus("nMu", 1);      tMu->SetBranchAddress("nMu", &nMu);
    tMu->SetBranchStatus("muPt", 1);     tMu->SetBranchAddress("muPt", &muPt);
    tMu->SetBranchStatus("muEta", 1);    tMu->SetBranchAddress("muEta", &muEta);
    tMu->SetBranchStatus("muPhi", 1);    tMu->SetBranchAddress("muPhi", &muPhi);
    tMu->SetBranchStatus("muCharge", 1); tMu->SetBranchAddress("muCharge", &muCharge);
    if (tMu->GetBranch("muIDTight"))
    {
        tMu->SetBranchStatus("muIDTight", 1);
        tMu->SetBranchAddress("muIDTight", &muIDTight);
    }
    else
    {
        std::cerr << "[FATAL] no muIDTight branch -- cannot run the tight-ID study.\n";
        return;
    }
    const bool hasCh  = (tMu->GetBranch("muPFChIso")  != nullptr);
    const bool hasNeu = (tMu->GetBranch("muPFNeuIso") != nullptr);
    const bool hasPho = (tMu->GetBranch("muPFPhoIso") != nullptr);
    const bool hasPU  = (tMu->GetBranch("muPFPUIso")  != nullptr);
    if (hasCh)  { tMu->SetBranchStatus("muPFChIso", 1);  tMu->SetBranchAddress("muPFChIso",  &muPFChIso);  }
    if (hasNeu) { tMu->SetBranchStatus("muPFNeuIso", 1); tMu->SetBranchAddress("muPFNeuIso", &muPFNeuIso); }
    if (hasPho) { tMu->SetBranchStatus("muPFPhoIso", 1); tMu->SetBranchAddress("muPFPhoIso", &muPFPhoIso); }
    if (hasPU)  { tMu->SetBranchStatus("muPFPUIso", 1);  tMu->SetBranchAddress("muPFPUIso",  &muPFPUIso);  }
    if (!hasCh || !hasNeu || !hasPho)
    {
        std::cerr << "[FATAL] muPFChIso/NeuIso/PhoIso branch missing -- no relIso to scan.\n";
        return;
    }
    if (!hasPU)
        std::cerr << "[WARN] no muPFPUIso branch; dbeta scan degenerates to the uncorrected one.\n";

    Int_t nPF = 0;
    std::vector<int> *pfId = nullptr;
    std::vector<float> *pfPt = nullptr, *pfEta = nullptr, *pfPhi = nullptr;
    tPF->SetBranchStatus("*", 0);
    tPF->SetBranchStatus("nPF", 1);  tPF->SetBranchAddress("nPF", &nPF);
    tPF->SetBranchStatus("pfId", 1); tPF->SetBranchAddress("pfId", &pfId);
    tPF->SetBranchStatus("pfPt", 1); tPF->SetBranchAddress("pfPt", &pfPt);
    tPF->SetBranchStatus("pfEta", 1); tPF->SetBranchAddress("pfEta", &pfEta);
    tPF->SetBranchStatus("pfPhi", 1); tPF->SetBranchAddress("pfPhi", &pfPhi);

    // ---------- selection (same surface as the old study / skim gates) ----------
    const double ptMin    = 15.0;
    const double etaMax   = 2.4;
    const double zMassLo  = 80.0, zMassHi = 100.0; // same window as isolation_ele.C
    const double metBkgMax = 5.0;                  // QCD-like: single muon + no MET
    const double cutSkim  = 0.15;                  // skim operating point

    auto passTight = [&](int i) -> bool
    {
        if (muPt->at(i) < ptMin) return false;
        if (std::abs(muEta->at(i)) > etaMax) return false;
        return (muIDTight->at(i) != 0);
    };

    // The two scanned definitions. idef: 0 = dbeta (primary), 1 = nodbeta (reference).
    auto relIso = [&](int i, int idef) -> double
    {
        const double pt = muPt->at(i);
        if (pt <= 0) return 999.0;
        if (idef == 0)
        {
            const double pu      = (muPFPUIso ? (double)muPFPUIso->at(i) : 0.0);
            const double neutral = std::max(0.0, (double)muPFNeuIso->at(i)
                                               + (double)muPFPhoIso->at(i) - 0.5 * pu);
            return (muPFChIso->at(i) + neutral) / pt;
        }
        return (muPFChIso->at(i) + muPFNeuIso->at(i) + muPFPhoIso->at(i)) / pt;
    };

    // ---------- continuous scan containers ----------
    const int nDef = 2;
    const char *defTags[nDef] = {"dbeta", "nodbeta"};
    const int nCuts = 200;
    std::vector<double> cuts(nCuts);
    for (int i = 0; i < nCuts; ++i)
        cuts[i] = (double)i / (double)(nCuts - 1); // 0..1

    std::vector<std::vector<double>> passOS(nDef, std::vector<double>(nCuts, 0.0));
    std::vector<std::vector<double>> totOS(nDef, std::vector<double>(nCuts, 0.0));
    std::vector<std::vector<double>> passSS(nDef, std::vector<double>(nCuts, 0.0));
    std::vector<std::vector<double>> totSS(nDef, std::vector<double>(nCuts, 0.0));
    std::vector<std::vector<double>> passMET(nDef, std::vector<double>(nCuts, 0.0));
    std::vector<std::vector<double>> totMET(nDef, std::vector<double>(nCuts, 0.0));

    long nOSpairs = 0, nSSpairs = 0, nQCDmu = 0;

    // ---------- event loop ----------
    const Long64_t nEv = std::min(tMu->GetEntries(), tPF->GetEntries());
    std::cout << "Processing entries: " << nEv << "\n";

    for (Long64_t e = 0; e < nEv; ++e)
    {
        if (e % 100000 == 0)
            std::cout << "Entries = " << e << std::endl;

        tMu->GetEntry(e);
        tPF->GetEntry(e);

        if (!muPt || !muEta || !muPhi || !muCharge || !muIDTight)
            continue;

        std::vector<int> idx;
        for (int i = 0; i < nMu; ++i)
            if (passTight(i))
                idx.push_back(i);
        if (idx.empty())
            continue;

        // ---- QCD-like background: exactly one tight muon + PF MET < 5 ----
        if (idx.size() == 1)
        {
            const double met = ComputePFMET_mu(pfId, pfPt, pfPhi).Mod();
            if (met < metBkgMax)
            {
                ++nQCDmu;
                const int i = idx[0];
                for (int idef = 0; idef < nDef; ++idef)
                {
                    const double riso = relIso(i, idef);
                    for (int ic = 0; ic < nCuts; ++ic)
                    {
                        totMET[idef][ic] += 1.0;
                        if (riso < cuts[ic])
                            passMET[idef][ic] += 1.0;
                    }
                }
            }
        }

        // ---- OS (signal) / SS (background) Z-window pairs ----
        std::vector<MuPairInfo> osPairs, ssPairs;
        for (size_t a = 0; a < idx.size(); ++a)
            for (size_t b = a + 1; b < idx.size(); ++b)
            {
                const int i = idx[a], j = idx[b];
                TLorentzVector v1, v2;
                v1.SetPtEtaPhiM(muPt->at(i), muEta->at(i), muPhi->at(i), 0.1056583745);
                v2.SetPtEtaPhiM(muPt->at(j), muEta->at(j), muPhi->at(j), 0.1056583745);
                const double m = (v1 + v2).M();
                if (m < zMassLo || m > zMassHi)
                    continue;
                if (muCharge->at(i) * muCharge->at(j) < 0)
                    osPairs.push_back({i, j, m});
                else
                    ssPairs.push_back({i, j, m});
            }
        nOSpairs += (long)osPairs.size();
        nSSpairs += (long)ssPairs.size();

        for (int idef = 0; idef < nDef; ++idef)
        {
            for (size_t p = 0; p < osPairs.size(); ++p)
            {
                const double ri = relIso(osPairs[p].i, idef);
                const double rj = relIso(osPairs[p].j, idef);
                for (int ic = 0; ic < nCuts; ++ic)
                {
                    totOS[idef][ic] += 2.0;
                    if (ri < cuts[ic]) passOS[idef][ic] += 1.0;
                    if (rj < cuts[ic]) passOS[idef][ic] += 1.0;
                }
            }
            for (size_t p = 0; p < ssPairs.size(); ++p)
            {
                const double ri = relIso(ssPairs[p].i, idef);
                const double rj = relIso(ssPairs[p].j, idef);
                for (int ic = 0; ic < nCuts; ++ic)
                {
                    totSS[idef][ic] += 2.0;
                    if (ri < cuts[ic]) passSS[idef][ic] += 1.0;
                    if (rj < cuts[ic]) passSS[idef][ic] += 1.0;
                }
            }
        }
    } // event loop

    std::cout << "Selected: " << nOSpairs << " OS Z-window pairs, "
              << nSSpairs << " SS pairs, "
              << nQCDmu << " single-mu MET<5 (QCD-like) muons\n";

    // ---------- outputs ----------
    TFile *fout = TFile::Open("./rootfile/IsoStudyOutputs_muon_tight.root", "RECREATE");
    if (!fout || fout->IsZombie())
    {
        std::cerr << "Cannot create rootfile/IsoStudyOutputs_muon_tight.root\n";
        return;
    }

    auto safeEff = [](double pass, double tot) -> double
    { return (tot > 0.0 ? pass / tot : 0.0); };

    std::cout << "\n============== CONTINUOUS relIso scan (muon, muIDTight) ==============\n";
    std::cout << "Signal = OS Z-window pairs; bkg(SS) = same-sign pairs; bkg(MET) = single-mu, MET<5.\n";
    std::cout << "Youden J = sigEff - bkgEff (maximize); current skim cut = " << cutSkim << ".\n";
    std::cout << Form("%-9s %8s %8s | %7s @cut %6s | %7s @cut %6s | %9s %9s\n",
                      "def", "aucSS", "aucMET", "J_SS", "", "J_MET", "", "effOS.15", "effMET.15");

    for (int idef = 0; idef < nDef; ++idef)
    {
        const char *tag = defTags[idef];
        std::vector<double> xcut(nCuts), effOS(nCuts), effSS(nCuts), effMET(nCuts);
        std::vector<double> pOS(nCuts), tOS(nCuts), pSS(nCuts), pMET(nCuts);
        for (int ic = 0; ic < nCuts; ++ic)
        {
            xcut[ic]   = cuts[ic];
            effOS[ic]  = safeEff(passOS[idef][ic],  totOS[idef][ic]);
            effSS[ic]  = safeEff(passSS[idef][ic],  totSS[idef][ic]);
            effMET[ic] = safeEff(passMET[idef][ic], totMET[idef][ic]);
            pOS[ic] = passOS[idef][ic];  tOS[ic]  = totOS[idef][ic];
            pSS[ic] = passSS[idef][ic];  pMET[ic] = passMET[idef][ic];
        }

        auto wg = [&](const char *pre, std::vector<double> &y)
        {
            TGraph *g = new TGraph(nCuts, xcut.data(), y.data());
            g->SetName(Form("%s_%s", pre, tag));
            g->Write();
        };
        wg("contEffOS", effOS);   wg("contEffSS", effSS);   wg("contEffMET", effMET);
        wg("contPassOS", pOS);    wg("contTotOS", tOS);     wg("contPassSS", pSS);
        wg("contPassMET", pMET);

        std::vector<double> xRocSS(nCuts), yRocSS(nCuts), xRocMET(nCuts), yRocMET(nCuts);
        for (int ic = 0; ic < nCuts; ++ic)
        {
            xRocSS[ic]  = effOS[ic]; yRocSS[ic]  = 1.0 - effSS[ic];
            xRocMET[ic] = effOS[ic]; yRocMET[ic] = 1.0 - effMET[ic];
        }
        TGraph *gRS = new TGraph(nCuts, xRocSS.data(), yRocSS.data());
        gRS->SetName(Form("contRocSS_%s", tag)); gRS->Write();
        TGraph *gRM = new TGraph(nCuts, xRocMET.data(), yRocMET.data());
        gRM->SetName(Form("contRocMET_%s", tag)); gRM->Write();

        const double aucSS  = rocAUC_mu(xRocSS, yRocSS);
        const double aucMET = rocAUC_mu(xRocMET, yRocMET);
        TParameter<double>(Form("contAucSS_%s", tag),  aucSS).Write();
        TParameter<double>(Form("contAucMET_%s", tag), aucMET).Write();

        int ibJSS = 0, ibJMET = 0;
        double JSS = -1.0, JMET = -1.0;
        for (int ic = 0; ic < nCuts; ++ic)
        {
            const double jss  = effOS[ic] - effSS[ic];
            const double jmet = effOS[ic] - effMET[ic];
            if (jss  > JSS)  { JSS  = jss;  ibJSS  = ic; }
            if (jmet > JMET) { JMET = jmet; ibJMET = ic; }
        }
        const int icCut = (int)(cutSkim * (nCuts - 1) + 0.5);
        std::cout << Form("%-9s %8.4f %8.4f | %7.4f %6.3f | %7.4f %6.3f | %9.4f %9.4f\n",
                          tag, aucSS, aucMET,
                          JSS, cuts[ibJSS], JMET, cuts[ibJMET],
                          effOS[icCut], effMET[icCut]);
    }
    std::cout << "=======================================================================\n";

    fout->Write();
    fout->Close();
    f->Close();

    std::cout << "Saved graphs to rootfile/IsoStudyOutputs_muon_tight.root\n";
    std::cout << "Now draw the summary: root -l -q 'plot_iso_summary.C+'\n";
}
