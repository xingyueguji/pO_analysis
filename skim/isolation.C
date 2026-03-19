// IsoStudy_pO.C
// ROOT macro to reproduce isolation efficiency + ROC/AUC
// using muonAnalyzer muonTree + particleFlowAnalyser pftree.
//
// Run: root -l -q 'IsoStudy_pO.C("HiForestMiniAOD.root")'

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <iostream>
#include <vector>
#include <algorithm>

struct PairInfo
{
    int i;
    int j;
    double mass;
};

static TVector2 ComputePFMET(std::vector<int> *pfId,
                             std::vector<float> *pfPt,
                             std::vector<float> *pfPhi)
{
    (void)pfId;
    double sumPx = 0.0, sumPy = 0.0;
    const int nPF = (pfPt ? (int)pfPt->size() : 0);

    for (int i = 0; i < nPF; ++i)
    {
        const double pt = pfPt->at(i);
        const double phi = pfPhi->at(i);
        sumPx += pt * std::cos(phi);
        sumPy += pt * std::sin(phi);
    }
    return TVector2(-sumPx, -sumPy);
}

static double dPhi(double a, double b)
{
    double d = a - b;
    while (d > TMath::Pi())
        d -= 2 * TMath::Pi();
    while (d < -TMath::Pi())
        d += 2 * TMath::Pi();
    return d;
}
static double dR(double eta1, double phi1, double eta2, double phi2)
{
    double de = eta1 - eta2;
    double dp = dPhi(phi1, phi2);
    return std::sqrt(de * de + dp * dp);
}

static TTree *getTree(TFile *f, const std::vector<TString> &paths)
{
    for (auto &p : paths)
    {
        TTree *t = (TTree *)f->Get(p);
        if (t)
        {
            std::cout << "Found tree: " << p << "\n";
            return t;
        }
    }
    return nullptr;
}

// A “tight-ish” muon selection you can tune.
// Uses branches available in your list.
static bool passTightMuon(int i, const std::vector<int> *muIDTight,
                          const std::vector<float> *muPt,
                          const std::vector<float> *muEta)
{
    if (!muIDTight || !muPt || !muEta)
        return false;
    if (muPt->at(i) < 15.0)
        return false;
    if (std::abs(muEta->at(i)) > 2.4)
        return false;
    return (muIDTight->at(i) != 0);
}

static bool passSoftMuon(int i, const std::vector<int> *muIDSoft,
                         const std::vector<float> *muPt,
                         const std::vector<float> *muEta)
{
    if (!muIDSoft || !muPt || !muEta)
        return false;
    if (muPt->at(i) < 15.0)
        return false;
    if (std::abs(muEta->at(i)) > 2.4)
        return false;
    return (muIDSoft->at(i) != 0);
}

// PF-ID meaning in HiForest particleFlowAnalyser is typically:
// 1=charged hadron, 2=e, 3=mu, 4=gamma, 5=neutral hadron, (6/7 HF stuff)
// If your mapping differs, adjust these predicates.
static bool isCH(int id) { return id == 1; }
static bool isNH(int id) { return id == 5; }
static bool isPH(int id) { return id == 4; }
static bool isMU(int id) { return id == 3; }
static bool isEL(int id) { return id == 2; }

/*
Particle::pdgId_	PFCandidate::particleId_	PFCandidate::ParticleType	Particle
0	0	X	unknown, or dummy
+211, -211	1	h	charged hadron
+11, -11	2	e	electron
+13, -13	3	mu	muon
22	4	gamma	photon
130	5	h0	neutral hadron
130	6	h_HF	hadronic energy in an HF tower
22	7	egamma_HF	electromagnetic energy in an HF tower
*/

// Compute isolation from PF candidates with variable cone R.
// caseMode:
// 0: ch + nh + pho
// 1: sum everything (except mu footprint)
// 2: deltaBeta PU correction (needs puAbs term provided externally)
static double computeIsoPF(double mu_eta, double mu_phi, double mu_pt,
                           const std::vector<int> *pfId,
                           const std::vector<float> *pfPt,
                           const std::vector<float> *pfEta,
                           const std::vector<float> *pfPhi,
                           double Rcone,
                           double vetoDR,
                           int caseMode,
                           double puAbs /*absolute PU pT in cone*/)
{
    if (!pfId || !pfPt || !pfEta || !pfPhi)
        return 999.0;
    if (mu_pt <= 0)
        return 999.0;

    double sumCH = 0, sumNH = 0, sumPH = 0, sumAll = 0;

    const int n = (int)pfId->size();
    for (int j = 0; j < n; j++)
    {
        const double dr = dR(mu_eta, mu_phi, pfEta->at(j), pfPhi->at(j));
        if (dr >= Rcone)
            continue;
        if (dr <= vetoDR)
            continue; // footprint veto (remove mu itself etc.)

        const int id = pfId->at(j);
        const double pt = pfPt->at(j);

        // always avoid counting PF muons in isolation unless you want them
        if (isMU(id))
            continue;
        if (isEL(id))
            continue;
        if (isCH(id))
            sumCH += pt;
        if (isNH(id))
            sumNH += pt;
        if (isPH(id))
            sumPH += pt;
        sumAll += pt;
    }

    if (caseMode == 0)
    {
        return (sumCH + sumNH + sumPH) / mu_pt;
    }
    else if (caseMode == 1)
    {
        return sumAll / mu_pt;
    }
    else
    {
        // Δβ: ch + max(0, nh+ph - 0.5*PU)
        double corr = (sumNH + sumPH) - 0.5 * puAbs;
        if (corr < 0)
            corr = 0;
        return (sumCH + corr) / mu_pt;
    }
}

// Main driver
void isolation(const char *fname = "root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2025.root", bool isSoft = true)
{
    TFile *f = TFile::Open(fname, "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "Cannot open: " << fname << "\n";
        return;
    }

    // Try common HiForest paths; adjust if your file uses different names
    TTree *tMu = getTree(f, {"ggHiNtuplizer/EventTree"});
    TTree *tPF = getTree(f, {"particleFlowAnalyser/pftree"});

    if (!tMu || !tPF)
    {
        std::cerr << "Missing muon or PF tree. Open the file in ROOT and check tree names.\n";
        return;
    }

    // ---------- Muon branches ----------
    Int_t nMu = 0;
    std::vector<float> *muPt = nullptr, *muEta = nullptr, *muPhi = nullptr, *muChi2NDF = nullptr, *muD0 = nullptr, *muDz = nullptr;
    std::vector<int> *muCharge = nullptr, *muIsGlobal = nullptr, *muIsPF = nullptr, *muMuonHits = nullptr, *muStations = nullptr, *muTrkLayers = nullptr, *muPixelHits = nullptr, *muIDTight = nullptr;
    std::vector<int> *muIDSoft = nullptr;
    std::vector<float> *muPFChIso = nullptr, *muPFPhoIso = nullptr, *muPFNeuIso = nullptr;

    // Optional PU iso (absolute) if present
    std::vector<float> *muPFPUIso = nullptr;
    bool hasMuPFPUIso = (tMu->GetBranch("muPFPUIso") != nullptr);

    tMu->SetBranchStatus("*", 0);

    tMu->SetBranchStatus("nMu", 1);
    tMu->SetBranchAddress("nMu", &nMu);
    tMu->SetBranchStatus("muPt", 1);
    tMu->SetBranchAddress("muPt", &muPt);
    tMu->SetBranchStatus("muEta", 1);
    tMu->SetBranchAddress("muEta", &muEta);
    tMu->SetBranchStatus("muPhi", 1);
    tMu->SetBranchAddress("muPhi", &muPhi);
    tMu->SetBranchStatus("muCharge", 1);
    tMu->SetBranchAddress("muCharge", &muCharge);

    tMu->SetBranchStatus("muIDTight", 1);
    tMu->SetBranchStatus("muIDSoft", 1);

    if (tMu->GetBranch("muIDTight"))
        tMu->SetBranchAddress("muIDTight", &muIDTight);

    if (tMu->GetBranch("muIDSoft"))
        tMu->SetBranchAddress("muIDSoft", &muIDSoft);

    tMu->SetBranchStatus("muPFPUIso", 1);
    if (hasMuPFPUIso)
        tMu->SetBranchAddress("muPFPUIso", &muPFPUIso);

    // Feb 22: added PF from gg branch to check dR, added soft variant.
    tMu->SetBranchStatus("muPFChIso", 1);
    tMu->SetBranchStatus("muPFPhoIso", 1);
    tMu->SetBranchStatus("muPFNeuIso", 1);

    tMu->SetBranchAddress("muPFChIso", &muPFChIso);
    tMu->SetBranchAddress("muPFPhoIso", &muPFPhoIso);
    tMu->SetBranchAddress("muPFNeuIso", &muPFNeuIso);

    // ---------- PF branches ----------
    Int_t nPF = 0;
    std::vector<int> *pfId = nullptr;
    std::vector<float> *pfPt = nullptr, *pfEta = nullptr, *pfPhi = nullptr;

    tPF->SetBranchStatus("*", 0);

    tPF->SetBranchStatus("nPF", 1);
    tPF->SetBranchAddress("nPF", &nPF);
    tPF->SetBranchStatus("pfId", 1);
    tPF->SetBranchAddress("pfId", &pfId);
    tPF->SetBranchStatus("pfPt", 1);
    tPF->SetBranchAddress("pfPt", &pfPt);
    tPF->SetBranchStatus("pfEta", 1);
    tPF->SetBranchAddress("pfEta", &pfEta);
    tPF->SetBranchStatus("pfPhi", 1);
    tPF->SetBranchAddress("pfPhi", &pfPhi);

    // ---------- Study setup ----------
    const std::vector<double> cones = {0.2, 0.3, 0.4};
    const int nCuts = 200;
    std::vector<double> cuts(nCuts);
    for (int i = 0; i < nCuts; i++)
        cuts[i] = (double)i / (double)(nCuts - 1); // 0..1

    std::vector<double> cuts_pT;
    for (double pt = 0.0; pt <= 100.0; pt += 5.0)
    {
        cuts_pT.push_back(pt);
    }
    int nCuts_pT = cuts_pT.size();

    // case names
    const int nCases = 3;
    const char *caseName[nCases] = {"CH+NH+PHO", "ALL", "DeltaBetaPU"};
    const double vetoDR = 0.01;

    // storage: effOS[case][cone][cut], effSS[case][cone][cut]
    std::vector<std::vector<std::vector<double>>> passOS(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts, 0)));
    std::vector<std::vector<std::vector<double>>> totOS(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts, 0)));
    std::vector<std::vector<std::vector<double>>> passSS(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts, 0)));
    std::vector<std::vector<std::vector<double>>> totSS(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts, 0)));
    std::vector<std::vector<std::vector<double>>> passMET(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts, 0)));
    std::vector<std::vector<std::vector<double>>> totMET(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts, 0)));

    std::vector<std::vector<std::vector<double>>> passOS_pT(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts_pT, 0)));
    std::vector<std::vector<std::vector<double>>> totOS_pT(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts_pT, 0)));
    std::vector<std::vector<std::vector<double>>> passSS_pT(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts_pT, 0)));
    std::vector<std::vector<std::vector<double>>> totSS_pT(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts_pT, 0)));
    std::vector<std::vector<std::vector<double>>> passMET_pT(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts_pT, 0)));
    std::vector<std::vector<std::vector<double>>> totMET_pT(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts_pT, 0)));

    std::vector<std::vector<std::vector<double>>> passOS_ggbranch(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts, 0)));
    std::vector<std::vector<std::vector<double>>> totOS_ggbranch(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts, 0)));
    std::vector<std::vector<std::vector<double>>> passSS_ggbranch(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts, 0)));
    std::vector<std::vector<std::vector<double>>> totSS_ggbranch(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts, 0)));
    std::vector<std::vector<std::vector<double>>> passMET_ggbranch(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts, 0)));
    std::vector<std::vector<std::vector<double>>> totMET_ggbranch(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts, 0)));

    // I do not need these for now

    /*std::vector<std::vector<std::vector<double>>> passOS_ggbranch_pT(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts_pT, 0)));
    std::vector<std::vector<std::vector<double>>> totOS_ggbranch_pT(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts_pT, 0)));
    std::vector<std::vector<std::vector<double>>> passSS_ggbranch_pT(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts_pT, 0)));
    std::vector<std::vector<std::vector<double>>> totSS_ggbranch_pT(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts_pT, 0)));
    std::vector<std::vector<std::vector<double>>> passMET_ggbranch_pT(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts_pT, 0)));
    std::vector<std::vector<std::vector<double>>> totMET_ggbranch_pT(nCases, std::vector<std::vector<double>>(cones.size(), std::vector<double>(nCuts_pT, 0)));*/

    // ---------- Event loop ----------
    Long64_t nEv = std::min(tMu->GetEntries(), tPF->GetEntries());
    std::cout << "Processing entries: " << nEv << "\n";

    auto relIsoFromBranches_case0 = [&](int imu) -> double
    {
        if (!muPt || !muPFChIso || !muPFNeuIso || !muPFPhoIso)
            return 999.0;
        const double pt = muPt->at(imu);
        if (pt <= 0)
            return 999.0;

        const double ch = muPFChIso->at(imu);
        const double neu = muPFNeuIso->at(imu);
        const double pho = muPFPhoIso->at(imu);

        return (ch + neu + pho) / pt;
    };

    for (Long64_t e = 0; e < nEv; e++)
    {
        if (e % 100000 == 0)
            cout << "Entries = " << e << endl;

        tMu->GetEntry(e);
        tPF->GetEntry(e);

        if (!muPt || !muEta || !muPhi || !muCharge)
            continue;

        TVector2 metv = ComputePFMET(pfId, pfPt, pfPhi);
        double met = metv.Mod();

        bool isBkgMET = (met < 5);

        // collect tight muons / or soft muon
        std::vector<int> idx;
        for (int i = 0; i < nMu; i++)
        {
            if (isSoft)
            {
                if (passSoftMuon(i, muIDSoft, muPt, muEta))
                    idx.push_back(i);
            }
            else
            {
                if (passTightMuon(i, muIDTight, muPt, muEta))
                    idx.push_back(i);
            }
        }
        if (idx.size() == 0)
            continue;

        // ---------------- MET background (define it exactly how you want) ----------------
        // Option A (what you said last): "single muon events with MET < 5 GeV"
        const bool useSingleMuonMETBkg = true;

        if (isBkgMET)
        {
            if (!useSingleMuonMETBkg || (idx.size() == 1))
            {
                // fill MET background using ALL tight muons in idx (for idx.size()==1 that's just one)
                for (int icase = 0; icase < nCases; ++icase)
                {
                    for (size_t ir = 0; ir < cones.size(); ++ir)
                    {
                        const double Rcone = cones[ir];

                        for (int imu : idx)
                        {
                            double puAbs = 0.0;
                            if (icase == 2 && hasMuPFPUIso)
                                puAbs = muPFPUIso->at(imu);

                            const double iso = computeIsoPF(
                                muEta->at(imu), muPhi->at(imu), muPt->at(imu),
                                pfId, pfPt, pfEta, pfPhi,
                                Rcone, vetoDR,
                                icase,
                                puAbs);

                            const double iso_br = relIsoFromBranches_case0(imu);

                            for (int icut = 0; icut < nCuts; ++icut)
                            {
                                const double cut = cuts[icut];
                                totMET[icase][ir][icut] += 1.0;
                                totMET_ggbranch[icase][ir][icut] += 1.0;
                                if (iso < cut)
                                    passMET[icase][ir][icut] += 1.0;
                                if (iso_br < cut)
                                    passMET_ggbranch[icase][ir][icut] += 1.0;
                            }

                            for (int icut_pT = 0; icut_pT < nCuts_pT; ++icut_pT)
                            {
                                // For now this has nothing to do with isolation, I am just testing efficiency of pT cuts, not sure if this is correct thing to do.
                                // Index for icase and ir is completely useless since I am not using those
                                const double cut_pT = cuts_pT[icut_pT];
                                totMET_pT[icase][ir][icut_pT] += 1.0;
                                if (muPt->at(imu) > cut_pT)
                                    passMET_pT[icase][ir][icut_pT] += 1.0;
                            }
                        }
                    }
                }
            }
        }

        // ---------------- OS / SS Z-like pairs ----------------
        if (idx.size() >= 2)
        {
            // collect Z-window pairs
            std::vector<PairInfo> osPairs;
            std::vector<PairInfo> ssPairs;

            for (size_t a = 0; a < idx.size(); ++a)
            {
                for (size_t b = a + 1; b < idx.size(); ++b)
                {
                    int i = idx[a];
                    int j = idx[b];

                    int qprod = muCharge->at(i) * muCharge->at(j);

                    TLorentzVector v1, v2;
                    v1.SetPtEtaPhiM(muPt->at(i), muEta->at(i), muPhi->at(i), 0.105658);
                    v2.SetPtEtaPhiM(muPt->at(j), muEta->at(j), muPhi->at(j), 0.105658);

                    double m = (v1 + v2).M();
                    if (m < 80.0 || m > 100.0)
                        continue;

                    if (qprod < 0)
                        osPairs.push_back({i, j, m});
                    else if (qprod > 0)
                        ssPairs.push_back({i, j, m});
                }
            }

            for (int icase = 0; icase < nCases; ++icase)
            {
                for (size_t ir = 0; ir < cones.size(); ++ir)
                {
                    const double Rcone = cones[ir];

                    // OS
                    for (const auto &p : osPairs)
                    {
                        double puAbs_i = 0.0, puAbs_j = 0.0;
                        if (icase == 2 && hasMuPFPUIso)
                        {
                            puAbs_i = muPFPUIso->at(p.i);
                            puAbs_j = muPFPUIso->at(p.j);
                        }

                        const double iso_i = computeIsoPF(muEta->at(p.i), muPhi->at(p.i), muPt->at(p.i),
                                                          pfId, pfPt, pfEta, pfPhi,
                                                          Rcone, vetoDR, icase, puAbs_i);

                        const double iso_j = computeIsoPF(muEta->at(p.j), muPhi->at(p.j), muPt->at(p.j),
                                                          pfId, pfPt, pfEta, pfPhi,
                                                          Rcone, vetoDR, icase, puAbs_j);

                        const double iso_i_br = relIsoFromBranches_case0(p.i);
                        const double iso_j_br = relIsoFromBranches_case0(p.j);

                        // ---- muon i ----
                        for (int icut = 0; icut < nCuts; ++icut)
                        {
                            const double cut = cuts[icut];
                            totOS[icase][ir][icut] += 1.0;
                            totOS_ggbranch[icase][ir][icut] += 1.0;
                            if (iso_i < cut)
                                passOS[icase][ir][icut] += 1.0;
                            if (iso_i_br < cut)
                                passOS_ggbranch[icase][ir][icut] += 1.0;
                        }

                        // ---- muon j ----
                        for (int icut = 0; icut < nCuts; ++icut)
                        {
                            const double cut = cuts[icut];
                            totOS[icase][ir][icut] += 1.0;
                            totOS_ggbranch[icase][ir][icut] += 1.0;
                            if (iso_j < cut)
                                passOS[icase][ir][icut] += 1.0;
                            if (iso_j_br < cut)
                                passOS_ggbranch[icase][ir][icut] += 1.0;
                        }

                        // ---- muon i ---- // pT
                        for (int icut_pT = 0; icut_pT < nCuts_pT; ++icut_pT)
                        {
                            const double cut_pT = cuts_pT[icut_pT];
                            totOS_pT[icase][ir][icut_pT] += 1.0;
                            if (muPt->at(p.i) > cut_pT)
                                passOS_pT[icase][ir][icut_pT] += 1.0;
                        }

                        // ---- muon j ---- // pT
                        for (int icut_pT = 0; icut_pT < nCuts_pT; ++icut_pT)
                        {
                            const double cut_pT = cuts_pT[icut_pT];
                            totOS_pT[icase][ir][icut_pT] += 1.0;
                            if (muPt->at(p.j) > cut_pT)
                                passOS_pT[icase][ir][icut_pT] += 1.0;
                        }
                    }

                    // SS
                    for (const auto &p : ssPairs)
                    {
                        double puAbs_i = 0.0, puAbs_j = 0.0;
                        if (icase == 2 && hasMuPFPUIso)
                        {
                            puAbs_i = muPFPUIso->at(p.i);
                            puAbs_j = muPFPUIso->at(p.j);
                        }

                        const double iso_i = computeIsoPF(muEta->at(p.i), muPhi->at(p.i), muPt->at(p.i),
                                                          pfId, pfPt, pfEta, pfPhi,
                                                          Rcone, vetoDR, icase, puAbs_i);

                        const double iso_j = computeIsoPF(muEta->at(p.j), muPhi->at(p.j), muPt->at(p.j),
                                                          pfId, pfPt, pfEta, pfPhi,
                                                          Rcone, vetoDR, icase, puAbs_j);

                        const double iso_i_br = relIsoFromBranches_case0(p.i);
                        const double iso_j_br = relIsoFromBranches_case0(p.j);

                        // muon i contribution
                        for (int icut = 0; icut < nCuts; ++icut)
                        {
                            const double cut = cuts[icut];
                            totSS[icase][ir][icut] += 1.0;
                            totSS_ggbranch[icase][ir][icut] += 1.0;
                            if (iso_i < cut)
                                passSS[icase][ir][icut] += 1.0;
                            if (iso_i_br < cut)
                                passSS_ggbranch[icase][ir][icut] += 1.0;
                        }

                        // muon j contribution
                        for (int icut = 0; icut < nCuts; ++icut)
                        {
                            const double cut = cuts[icut];
                            totSS[icase][ir][icut] += 1.0;
                            totSS_ggbranch[icase][ir][icut] += 1.0;
                            if (iso_j < cut)
                                passSS[icase][ir][icut] += 1.0;
                            if (iso_j_br < cut)
                                passSS_ggbranch[icase][ir][icut] += 1.0;
                        }

                        // muon i contribution // pT
                        for (int icut_pT = 0; icut_pT < nCuts_pT; ++icut_pT)
                        {
                            const double cut_pT = cuts_pT[icut_pT];
                            totSS_pT[icase][ir][icut_pT] += 1.0;
                            if (muPt->at(p.i) > cut_pT)
                                passSS_pT[icase][ir][icut_pT] += 1.0;
                        }

                        // muon j contribution // pT
                        for (int icut_pT = 0; icut_pT < nCuts_pT; ++icut_pT)
                        {
                            const double cut_pT = cuts_pT[icut_pT];
                            totSS_pT[icase][ir][icut_pT] += 1.0;
                            if (muPt->at(p.j) > cut_pT)
                                passSS_pT[icase][ir][icut_pT] += 1.0;
                        }
                    }
                }
            }
        }
    } // End of Event Loop

    // =========================
    // Save plots to ROOT file
    // =========================
    TFile *fout;
    if (isSoft)
    {
        fout = TFile::Open("./rootfile/IsoStudyOutputs_soft.root", "RECREATE");
    }
    else
    {
        fout = TFile::Open("./rootfile/IsoStudyOutputs.root", "RECREATE");
    }
    if (!fout || fout->IsZombie())
    {
        std::cerr << "Cannot create output ROOT file IsoStudyOutputs.root\n";
        return;
    }

    auto safeEff = [](double pass, double tot) -> double
    {
        return (tot > 0.0 ? pass / tot : 0.0);
    };

    auto rocAUC = [](const std::vector<double> &x, const std::vector<double> &y) -> double
    {
        // trapezoid integral in x (assumes x is monotonic-ish; we will sort by x)
        std::vector<std::pair<double, double>> xy;
        xy.reserve(x.size());
        for (size_t i = 0; i < x.size(); ++i)
            xy.push_back({x[i], y[i]});
        std::sort(xy.begin(), xy.end(),
                  [](auto &a, auto &b)
                  { return a.first < b.first; });

        double auc = 0.0;
        for (size_t i = 1; i < xy.size(); ++i)
        {
            const double x0 = xy[i - 1].first, x1 = xy[i].first;
            const double y0 = xy[i - 1].second, y1 = xy[i].second;
            auc += 0.5 * (y0 + y1) * (x1 - x0);
        }
        return auc; // in [0,1] typically
    };

    for (int icase = 0; icase < nCases; ++icase)
    {
        for (size_t ir = 0; ir < cones.size(); ++ir)
        {
            const double Rcone = cones[ir];

            // Prepare arrays for graphs
            std::vector<double> xCut(nCuts), yEffOS(nCuts), yEffSS(nCuts), yEffMET(nCuts);

            for (int icut = 0; icut < nCuts; ++icut)
            {
                xCut[icut] = cuts[icut];
                yEffOS[icut] = safeEff(passOS[icase][ir][icut], totOS[icase][ir][icut]);
                yEffSS[icut] = safeEff(passSS[icase][ir][icut], totSS[icase][ir][icut]);
                yEffMET[icut] = safeEff(passMET[icase][ir][icut], totMET[icase][ir][icut]);
            }

            // Name tag: case + cone
            TString tag;
            tag.Form("case%d_R%.1f", icase, Rcone); // e.g. case0_R0.3

            // Efficiency vs cut graphs
            TGraph *gEffOS = new TGraph(nCuts, xCut.data(), yEffOS.data());
            TGraph *gEffSS = new TGraph(nCuts, xCut.data(), yEffSS.data());
            TGraph *gEffMET = new TGraph(nCuts, xCut.data(), yEffMET.data());

            gEffOS->SetName(TString("effOS_") + tag);
            gEffSS->SetName(TString("effSS_") + tag);
            gEffMET->SetName(TString("effMET_") + tag);

            gEffOS->Write();
            gEffSS->Write();
            gEffMET->Write();

            // ---- ROC from SS background: x=bkg eff (SS), y=sig eff (OS) ----
            std::vector<double> xROC_SS(nCuts), yROC_SS(nCuts);
            for (int icut = 0; icut < nCuts; ++icut)
            {
                xROC_SS[icut] = yEffOS[icut];
                yROC_SS[icut] = 1.0 - yEffSS[icut];
            }
            TGraph *gROC_SS = new TGraph(nCuts, xROC_SS.data(), yROC_SS.data());
            gROC_SS->SetName(TString("rocSS_") + tag);

            // ---- ROC from MET background: x=bkg eff (MET), y=sig eff (OS) ----
            std::vector<double> xROC_MET(nCuts), yROC_MET(nCuts);
            for (int icut = 0; icut < nCuts; ++icut)
            {
                xROC_MET[icut] = yEffOS[icut];
                yROC_MET[icut] = 1.0 - yEffMET[icut];
            }
            TGraph *gROC_MET = new TGraph(nCuts, xROC_MET.data(), yROC_MET.data());
            gROC_MET->SetName(TString("rocMET_") + tag);

            gROC_SS->Write();
            gROC_MET->Write();

            // AUC numbers
            const double aucSS = rocAUC(xROC_SS, yROC_SS);
            const double aucMET = rocAUC(xROC_MET, yROC_MET);

            TParameter<double> pAucSS(TString("aucSS_") + tag, aucSS);
            TParameter<double> pAucMET(TString("aucMET_") + tag, aucMET);

            pAucSS.Write();
            pAucMET.Write();
        }
    }

    fout->Write();
    fout->Close();

    TFile *fout_pT;
    if (isSoft)
    {
        fout_pT = TFile::Open("./rootfile/PtcutStudyOutputs_soft.root", "RECREATE");
    }
    else
    {
        fout_pT = TFile::Open("./rootfile/PtcutStudyOutputs.root", "RECREATE");
    }
    if (!fout_pT || fout_pT->IsZombie())
    {
        std::cerr << "Cannot create output ROOT file PtcutStudyOutputs.root\n";
        return;
    }

    for (int icase = 0; icase < nCases; ++icase)
    {
        for (size_t ir = 0; ir < cones.size(); ++ir)
        {
            const double Rcone = cones[ir];

            // Prepare arrays for graphs
            std::vector<double> xCut(nCuts_pT), yEffOS(nCuts_pT), yEffSS(nCuts_pT), yEffMET(nCuts_pT);

            for (int icut = 0; icut < nCuts_pT; ++icut)
            {
                xCut[icut] = cuts_pT[icut];
                yEffOS[icut] = safeEff(passOS_pT[icase][ir][icut], totOS_pT[icase][ir][icut]);
                yEffSS[icut] = safeEff(passSS_pT[icase][ir][icut], totSS_pT[icase][ir][icut]);
                yEffMET[icut] = safeEff(passMET_pT[icase][ir][icut], totMET_pT[icase][ir][icut]);
            }

            // Name tag: case + cone
            TString tag;
            tag.Form("case%d_R%.1f", icase, Rcone); // e.g. case0_R0.3

            // Efficiency vs cut graphs
            TGraph *gEffOS = new TGraph(nCuts_pT, xCut.data(), yEffOS.data());
            TGraph *gEffSS = new TGraph(nCuts_pT, xCut.data(), yEffSS.data());
            TGraph *gEffMET = new TGraph(nCuts_pT, xCut.data(), yEffMET.data());

            gEffOS->SetName(TString("effOSpT_") + tag);
            gEffSS->SetName(TString("effSSpT_") + tag);
            gEffMET->SetName(TString("effMETpT_") + tag);

            gEffOS->Write();
            gEffSS->Write();
            gEffMET->Write();

            // ---- ROC from SS background: x=bkg eff (SS), y=sig eff (OS) ----
            std::vector<double> xROC_SS(nCuts_pT), yROC_SS(nCuts_pT);
            for (int icut = 0; icut < nCuts_pT; ++icut)
            {
                xROC_SS[icut] = yEffOS[icut];
                yROC_SS[icut] = 1.0 - yEffSS[icut];
            }
            TGraph *gROC_SS = new TGraph(nCuts_pT, xROC_SS.data(), yROC_SS.data());
            gROC_SS->SetName(TString("rocSSpT_") + tag);

            // ---- ROC from MET background: x=bkg eff (MET), y=sig eff (OS) ----
            std::vector<double> xROC_MET(nCuts_pT), yROC_MET(nCuts_pT);
            for (int icut = 0; icut < nCuts_pT; ++icut)
            {
                xROC_MET[icut] = yEffOS[icut];
                yROC_MET[icut] = 1.0 - yEffMET[icut];
            }
            TGraph *gROC_MET = new TGraph(nCuts_pT, xROC_MET.data(), yROC_MET.data());
            gROC_MET->SetName(TString("rocMETpT_") + tag);

            gROC_SS->Write();
            gROC_MET->Write();

            // AUC numbers
            const double aucSS = rocAUC(xROC_SS, yROC_SS);
            const double aucMET = rocAUC(xROC_MET, yROC_MET);

            TParameter<double> pAucSS(TString("aucSSpT_") + tag, aucSS);
            TParameter<double> pAucMET(TString("aucMETpT_") + tag, aucMET);

            pAucSS.Write();
            pAucMET.Write();
        }
    }

    fout_pT->Write();
    fout_pT->Close();

    TFile *fout_ggbranch;

    if (isSoft)
    {
        fout_ggbranch = TFile::Open("./rootfile/ggbranchStudyOutputs_soft.root", "RECREATE");
    }
    else
    {
        fout_ggbranch = TFile::Open("./rootfile/ggbranchStudyOutputs.root", "RECREATE");
    }

    if (!fout_ggbranch || fout_ggbranch->IsZombie())
    {
        std::cerr << "Cannot create output ROOT file IsoStudyOutputs.root\n";
        return;
    }

    for (int icase = 0; icase < nCases; ++icase)
    {
        for (size_t ir = 0; ir < cones.size(); ++ir)
        {
            const double Rcone = cones[ir];

            // Prepare arrays for graphs
            std::vector<double> xCut(nCuts), yEffOS(nCuts), yEffSS(nCuts), yEffMET(nCuts);

            for (int icut = 0; icut < nCuts; ++icut)
            {
                xCut[icut] = cuts[icut];
                yEffOS[icut] = safeEff(passOS_ggbranch[icase][ir][icut], totOS_ggbranch[icase][ir][icut]);
                yEffSS[icut] = safeEff(passSS_ggbranch[icase][ir][icut], totSS_ggbranch[icase][ir][icut]);
                yEffMET[icut] = safeEff(passMET_ggbranch[icase][ir][icut], totMET_ggbranch[icase][ir][icut]);
            }

            // Name tag: case + cone
            TString tag;
            tag.Form("case%d_R%.1f", icase, Rcone); // e.g. case0_R0.3

            // Efficiency vs cut graphs
            TGraph *gEffOS = new TGraph(nCuts, xCut.data(), yEffOS.data());
            TGraph *gEffSS = new TGraph(nCuts, xCut.data(), yEffSS.data());
            TGraph *gEffMET = new TGraph(nCuts, xCut.data(), yEffMET.data());

            gEffOS->SetName(TString("effOSgg_") + tag);
            gEffSS->SetName(TString("effSSgg_") + tag);
            gEffMET->SetName(TString("effMETgg_") + tag);

            gEffOS->Write();
            gEffSS->Write();
            gEffMET->Write();

            // ---- ROC from SS background: x=bkg eff (SS), y=sig eff (OS) ----
            std::vector<double> xROC_SS(nCuts), yROC_SS(nCuts);
            for (int icut = 0; icut < nCuts; ++icut)
            {
                xROC_SS[icut] = yEffOS[icut];
                yROC_SS[icut] = 1.0 - yEffSS[icut];
            }
            TGraph *gROC_SS = new TGraph(nCuts, xROC_SS.data(), yROC_SS.data());
            gROC_SS->SetName(TString("rocSSgg_") + tag);

            // ---- ROC from MET background: x=bkg eff (MET), y=sig eff (OS) ----
            std::vector<double> xROC_MET(nCuts), yROC_MET(nCuts);
            for (int icut = 0; icut < nCuts; ++icut)
            {
                xROC_MET[icut] = yEffOS[icut];
                yROC_MET[icut] = 1.0 - yEffMET[icut];
            }
            TGraph *gROC_MET = new TGraph(nCuts, xROC_MET.data(), yROC_MET.data());
            gROC_MET->SetName(TString("rocMETgg_") + tag);

            gROC_SS->Write();
            gROC_MET->Write();

            // AUC numbers
            const double aucSS = rocAUC(xROC_SS, yROC_SS);
            const double aucMET = rocAUC(xROC_MET, yROC_MET);

            TParameter<double> pAucSS(TString("aucSSgg_") + tag, aucSS);
            TParameter<double> pAucMET(TString("aucMETgg_") + tag, aucMET);

            pAucSS.Write();
            pAucMET.Write();
        }
    }

    fout_ggbranch->Write();
    fout_ggbranch->Close();

    std::cout << "Saved graphs to IsoStudyOutputs.root\n";
}