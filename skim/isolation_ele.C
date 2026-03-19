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
                          const std::vector<float> *elePt,
                          const std::vector<float> *eleEta)
{
    if (!muIDTight || !elePt || !eleEta)
        return false;
    if (elePt->at(i) < 15.0)
        return false;
    if (std::abs(eleEta->at(i)) > 2.4)
        return false;
    return (muIDTight->at(i) != 0);
}

static bool passSoftMuon(int i, const std::vector<int> *muIDSoft,
                         const std::vector<float> *elePt,
                         const std::vector<float> *eleEta)
{
    if (!muIDSoft || !elePt || !eleEta)
        return false;
    if (elePt->at(i) < 15.0)
        return false;
    if (std::abs(eleEta->at(i)) > 2.4)
        return false;
    return (muIDSoft->at(i) != 0);
}

static bool passEleID(int i, const std::vector<int> *eleID,
                      const std::vector<float> *elePt,
                      const std::vector<float> *eleEta)
{
    if (!eleID || !elePt || !eleEta)
        return false;
    if (elePt->at(i) < 15.0)
        return false;
    if (std::abs(eleEta->at(i)) > 2.4)
        return false;
    return (eleID->at(i) != 0);
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
void isolation_ele(const char *fname = "root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2025.root")
{
    TFile *f = TFile::Open(fname, "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "Cannot open: " << fname << "\n";
        return;
    }

    // Try common HiForest paths; adjust if your file uses different names
    TTree *tEle = getTree(f, {"ggHiNtuplizer/EventTree"});
    TTree *tPF = getTree(f, {"particleFlowAnalyser/pftree"});

    if (!tEle || !tPF)
    {
        std::cerr << "Missing muon or PF tree. Open the file in ROOT and check tree names.\n";
        return;
    }

    // ---------- electron branches ----------
    Int_t nEle = 0;
    std::vector<float> *elePt = nullptr, *eleEta = nullptr, *elePhi = nullptr, *eleDz = nullptr;
    std::vector<int> *eleCharge = nullptr, *eleMVAIdWP80 = nullptr, *eleMVAIdWP85 = nullptr, *eleMVAIdWP90 = nullptr, *eleMVAIdWP95 = nullptr;
    std::vector<int> *eleCutIdWP70 = nullptr, *eleCutIdWP80 = nullptr, *eleCutIdWP90 = nullptr, *eleCutIdWP95 = nullptr;
    std::vector<int> *eleMVAIsoWP80 = nullptr, *eleMVAIsoWP85 = nullptr, *eleMVAIsoWP90 = nullptr, *eleMVAIsoWP95 = nullptr;

    tEle->SetBranchStatus("*", 0);

    tEle->SetBranchStatus("nEle", 1);
    tEle->SetBranchAddress("nEle", &nEle);
    tEle->SetBranchStatus("elePt", 1);
    tEle->SetBranchAddress("elePt", &elePt);
    tEle->SetBranchStatus("eleEta", 1);
    tEle->SetBranchAddress("eleEta", &eleEta);
    tEle->SetBranchStatus("elePhi", 1);
    tEle->SetBranchAddress("elePhi", &elePhi);
    tEle->SetBranchStatus("eleCharge", 1);
    tEle->SetBranchAddress("eleCharge", &eleCharge);

    tEle->SetBranchStatus("eleMVAIdWP80", 1);
    tEle->SetBranchStatus("eleMVAIdWP85", 1);
    tEle->SetBranchStatus("eleMVAIdWP90", 1);
    tEle->SetBranchStatus("eleMVAIdWP95", 1);
    tEle->SetBranchAddress("eleMVAIdWP80", &eleMVAIdWP80);
    tEle->SetBranchAddress("eleMVAIdWP85", &eleMVAIdWP85);
    tEle->SetBranchAddress("eleMVAIdWP90", &eleMVAIdWP90);
    tEle->SetBranchAddress("eleMVAIdWP95", &eleMVAIdWP95);

    tEle->SetBranchStatus("eleCutIdWP70", 1);
    tEle->SetBranchStatus("eleCutIdWP80", 1);
    tEle->SetBranchStatus("eleCutIdWP90", 1);
    tEle->SetBranchStatus("eleCutIdWP95", 1);

    tEle->SetBranchAddress("eleCutIdWP70", &eleCutIdWP70);
    tEle->SetBranchAddress("eleCutIdWP80", &eleCutIdWP80);
    tEle->SetBranchAddress("eleCutIdWP90", &eleCutIdWP90);
    tEle->SetBranchAddress("eleCutIdWP95", &eleCutIdWP95);

    tEle->SetBranchStatus("eleMVAIsoWP80", 1);
    tEle->SetBranchStatus("eleMVAIsoWP85", 1);
    tEle->SetBranchStatus("eleMVAIsoWP90", 1);
    tEle->SetBranchStatus("eleMVAIsoWP95", 1);

    tEle->SetBranchAddress("eleMVAIsoWP80", &eleMVAIsoWP80);
    tEle->SetBranchAddress("eleMVAIsoWP85", &eleMVAIsoWP85);
    tEle->SetBranchAddress("eleMVAIsoWP90", &eleMVAIsoWP90);
    tEle->SetBranchAddress("eleMVAIsoWP95", &eleMVAIsoWP95);

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

    // ---------- Study setup --------- //

    // Different ID (4 + 4 ID), Different MVA (4)

    const int nID = 8;
    const int nIso = 4;

    std::vector<std::vector<double>> table(
        nID,
        std::vector<double>(nIso, 0.0));

    // Total 8 * 4 = 32 cases

    // Ordering:
    // ID: MVA ID first, then CUT ID

    std::vector<std::vector<double>> passOS(nID, std::vector<double>(nIso, 0.0));
    std::vector<std::vector<double>> totOS(nID, std::vector<double>(nIso, 0.0));
    std::vector<std::vector<double>> passSS(nID, std::vector<double>(nIso, 0.0));
    std::vector<std::vector<double>> totSS(nID, std::vector<double>(nIso, 0.0));
    std::vector<std::vector<double>> passMET(nID, std::vector<double>(nIso, 0.0));
    std::vector<std::vector<double>> totMET(nID, std::vector<double>(nIso, 0.0));

    // ---------- Event loop ----------
    Long64_t nEv = std::min(tEle->GetEntries(), tPF->GetEntries());
    std::cout << "Processing entries: " << nEv << "\n";

    for (Long64_t e = 0; e < nEv; e++)
    {
        if (e % 100000 == 0)
            cout << "Entries = " << e << endl;

        tEle->GetEntry(e);
        tPF->GetEntry(e);

        if (!elePt || !eleEta || !elePhi || !eleCharge)
            continue;

        TVector2 metv = ComputePFMET(pfId, pfPt, pfPhi);
        double met = metv.Mod();

        bool isBkgMET = (met < 5);

        // collect tight muons / or soft muon
        std::vector<int> idx_eleMVAIdWP80;
        std::vector<int> idx_eleMVAIdWP85;
        std::vector<int> idx_eleMVAIdWP90;
        std::vector<int> idx_eleMVAIdWP95;

        std::vector<int> idx_eleCutIdWP70;
        std::vector<int> idx_eleCutIdWP80;
        std::vector<int> idx_eleCutIdWP90;
        std::vector<int> idx_eleCutIdWP95;

        for (int i = 0; i < nEle; i++)
        {
            if (passEleID(i, eleMVAIdWP80, elePt, eleEta))
                idx_eleMVAIdWP80.push_back(i);
            if (passEleID(i, eleMVAIdWP85, elePt, eleEta))
                idx_eleMVAIdWP85.push_back(i);
            if (passEleID(i, eleMVAIdWP90, elePt, eleEta))
                idx_eleMVAIdWP90.push_back(i);
            if (passEleID(i, eleMVAIdWP95, elePt, eleEta))
                idx_eleMVAIdWP95.push_back(i);
            if (passEleID(i, eleCutIdWP70, elePt, eleEta))
                idx_eleCutIdWP70.push_back(i);
            if (passEleID(i, eleCutIdWP80, elePt, eleEta))
                idx_eleCutIdWP80.push_back(i);
            if (passEleID(i, eleCutIdWP90, elePt, eleEta))
                idx_eleCutIdWP90.push_back(i);
            if (passEleID(i, eleCutIdWP95, elePt, eleEta))
                idx_eleCutIdWP95.push_back(i);
        }

        // collect all ID selections into one container
        std::vector<std::vector<int> *> idLists = {
            &idx_eleMVAIdWP80,
            &idx_eleMVAIdWP85,
            &idx_eleMVAIdWP90,
            &idx_eleMVAIdWP95,
            &idx_eleCutIdWP70,
            &idx_eleCutIdWP80,
            &idx_eleCutIdWP90,
            &idx_eleCutIdWP95};

        // Option A: "single ID cut events with MET < 5 GeV"
        if (isBkgMET)
        {
            for (int id = 0; id < nID; id++)
            {
                auto &idxList = *idLists[id];

                // require exactly one electron passing this ID
                if (idxList.size() != 1)
                    continue;

                int i = idxList[0];

                for (int iso = 0; iso < nIso; iso++)
                {
                    bool passIso = false;

                    if (iso == 0)
                        passIso = (eleMVAIsoWP80->at(i));
                    if (iso == 1)
                        passIso = (eleMVAIsoWP85->at(i));
                    if (iso == 2)
                        passIso = (eleMVAIsoWP90->at(i));
                    if (iso == 3)
                        passIso = (eleMVAIsoWP95->at(i));

                    // total events
                    totMET[id][iso]++;

                    if (passIso)
                        passMET[id][iso]++;
                }
            }
        }

        // ---------------- OS / SS Z-like pairs ----------------

        for (int id = 0; id < nID; ++id)
        {
            auto &idx = *idLists[id];

            std::vector<PairInfo> osPairs;
            std::vector<PairInfo> ssPairs;

            // build Z-window pairs
            for (size_t a = 0; a < idx.size(); ++a)
            {
                for (size_t b = a + 1; b < idx.size(); ++b)
                {
                    int i = idx[a];
                    int j = idx[b];

                    int qprod = eleCharge->at(i) * eleCharge->at(j);

                    TLorentzVector v1, v2;
                    v1.SetPtEtaPhiM(elePt->at(i), eleEta->at(i), elePhi->at(i), 0.000511);
                    v2.SetPtEtaPhiM(elePt->at(j), eleEta->at(j), elePhi->at(j), 0.000511);

                    double m = (v1 + v2).M();

                    if (m < 80.0 || m > 100.0)
                        continue;

                    if (qprod < 0)
                        osPairs.push_back({i, j, m});
                    else
                        ssPairs.push_back({i, j, m});
                }
            }

            // evaluate isolation WPs
            for (int iso = 0; iso < nIso; ++iso)
            {
                // ---- OS pairs ----
                for (const auto &p : osPairs)
                {
                    bool passIso_i = false;
                    bool passIso_j = false;

                    if (iso == 0)
                    {
                        passIso_i = eleMVAIsoWP80->at(p.i);
                        passIso_j = eleMVAIsoWP80->at(p.j);
                    }
                    if (iso == 1)
                    {
                        passIso_i = eleMVAIsoWP85->at(p.i);
                        passIso_j = eleMVAIsoWP85->at(p.j);
                    }
                    if (iso == 2)
                    {
                        passIso_i = eleMVAIsoWP90->at(p.i);
                        passIso_j = eleMVAIsoWP90->at(p.j);
                    }
                    if (iso == 3)
                    {
                        passIso_i = eleMVAIsoWP95->at(p.i);
                        passIso_j = eleMVAIsoWP95->at(p.j);
                    }

                    totOS[id][iso] += 2.0;

                    if (passIso_i)
                        passOS[id][iso] += 1.0;
                    if (passIso_j)
                        passOS[id][iso] += 1.0;
                }

                // ---- SS pairs ----
                for (const auto &p : ssPairs)
                {
                    bool passIso_i = false;
                    bool passIso_j = false;

                    if (iso == 0)
                    {
                        passIso_i = eleMVAIsoWP80->at(p.i);
                        passIso_j = eleMVAIsoWP80->at(p.j);
                    }
                    if (iso == 1)
                    {
                        passIso_i = eleMVAIsoWP85->at(p.i);
                        passIso_j = eleMVAIsoWP85->at(p.j);
                    }
                    if (iso == 2)
                    {
                        passIso_i = eleMVAIsoWP90->at(p.i);
                        passIso_j = eleMVAIsoWP90->at(p.j);
                    }
                    if (iso == 3)
                    {
                        passIso_i = eleMVAIsoWP95->at(p.i);
                        passIso_j = eleMVAIsoWP95->at(p.j);
                    }

                    totSS[id][iso] += 2.0;

                    if (passIso_i)
                        passSS[id][iso] += 1.0;
                    if (passIso_j)
                        passSS[id][iso] += 1.0;
                }
            }
        }

    } // End of Event Loop

    // =========================
    // Save plots to ROOT file
    // =========================
    TFile *fout;

    fout = TFile::Open("./rootfile/IsoStudyOutputs_electron.root", "RECREATE");

    if (!fout || fout->IsZombie())
    {
        std::cerr << "Cannot create output ROOT file IsoStudyOutputs_electron.root\n";
        return;
    }

    auto safeEff = [](double pass, double tot) -> double
    {
        return (tot > 0.0 ? pass / tot : 0.0);
    };

    auto rocAUC = [](const std::vector<double> &x, const std::vector<double> &y) -> double
    {
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
            const double x0 = xy[i - 1].first;
            const double x1 = xy[i].first;
            const double y0 = xy[i - 1].second;
            const double y1 = xy[i].second;

            auc += 0.5 * (y0 + y1) * (x1 - x0);
        }

        return auc;
    };

    std::vector<std::string> idNames = {
        "MVAIdWP80",
        "MVAIdWP85",
        "MVAIdWP90",
        "MVAIdWP95",
        "CutIdWP70",
        "CutIdWP80",
        "CutIdWP90",
        "CutIdWP95"};

    std::vector<std::string> isoNames = {
        "IsoWP80",
        "IsoWP85",
        "IsoWP90",
        "IsoWP95"};

    for (int id = 0; id < nID; ++id)
    {
        std::vector<double> xIso(nIso);
        std::vector<double> yEffOS(nIso);
        std::vector<double> yEffSS(nIso);
        std::vector<double> yEffMET(nIso);

        std::vector<double> yPassOS(nIso), yTotOS(nIso);
        std::vector<double> yPassSS(nIso), yTotSS(nIso);
        std::vector<double> yPassMET(nIso), yTotMET(nIso);

        std::vector<double> yEffBkg(nIso);
        std::vector<double> xROC_comb(nIso), yROC_comb(nIso);

        TString tag;

        for (int iso = 0; iso < nIso; ++iso)
        {
            xIso[iso] = iso + 1;

            yEffOS[iso] = safeEff(passOS[id][iso], totOS[id][iso]);
            yEffSS[iso] = safeEff(passSS[id][iso], totSS[id][iso]);
            yEffMET[iso] = safeEff(passMET[id][iso], totMET[id][iso]);

            double passBkg = passSS[id][iso] + passMET[id][iso];
            double totBkg = totSS[id][iso] + totMET[id][iso];

            yEffBkg[iso] = safeEff(passBkg, totBkg);

            xROC_comb[iso] = yEffOS[iso];
            yROC_comb[iso] = 1.0 - yEffBkg[iso];

            yPassOS[iso] = passOS[id][iso];
            yTotOS[iso] = totOS[id][iso];

            yPassSS[iso] = passSS[id][iso];
            yTotSS[iso] = totSS[id][iso];

            yPassMET[iso] = passMET[id][iso];
            yTotMET[iso] = totMET[id][iso];

            tag.Form("%s", idNames[id].c_str());
        }

        // Efficiency graphs
        TGraph *gEffOS = new TGraph(nIso, xIso.data(), yEffOS.data());
        TGraph *gEffSS = new TGraph(nIso, xIso.data(), yEffSS.data());
        TGraph *gEffMET = new TGraph(nIso, xIso.data(), yEffMET.data());

        gEffOS->SetName(TString("effOS_") + tag);
        gEffSS->SetName(TString("effSS_") + tag);
        gEffMET->SetName(TString("effMET_") + tag);

        gEffOS->Write();
        gEffSS->Write();
        gEffMET->Write();

        TGraph *gPassOS = new TGraph(nIso, xIso.data(), yPassOS.data());
        TGraph *gTotOS = new TGraph(nIso, xIso.data(), yTotOS.data());

        gPassOS->SetName(TString("passOS_") + tag);
        gTotOS->SetName(TString("totOS_") + tag);

        gPassOS->Write();
        gTotOS->Write();

        TGraph *gPassSS = new TGraph(nIso, xIso.data(), yPassSS.data());
        TGraph *gTotSS = new TGraph(nIso, xIso.data(), yTotSS.data());

        gPassSS->SetName(TString("passSS_") + tag);
        gTotSS->SetName(TString("totSS_") + tag);

        gPassSS->Write();
        gTotSS->Write();

        TGraph *gPassMET = new TGraph(nIso, xIso.data(), yPassMET.data());
        TGraph *gTotMET = new TGraph(nIso, xIso.data(), yTotMET.data());

        gPassMET->SetName(TString("passMET_") + tag);
        gTotMET->SetName(TString("totMET_") + tag);

        gPassMET->Write();
        gTotMET->Write();

        // ---- ROC from SS background ----
        std::vector<double> xROC_SS(nIso), yROC_SS(nIso);

        for (int iso = 0; iso < nIso; ++iso)
        {
            xROC_SS[iso] = yEffOS[iso];
            yROC_SS[iso] = 1.0 - yEffSS[iso];
        }

        TGraph *gROC_SS = new TGraph(nIso, xROC_SS.data(), yROC_SS.data());
        gROC_SS->SetName(TString("rocSS_") + tag);
        gROC_SS->Write();

        // ---- ROC from MET background ----
        std::vector<double> xROC_MET(nIso), yROC_MET(nIso);

        for (int iso = 0; iso < nIso; ++iso)
        {
            xROC_MET[iso] = yEffOS[iso];
            yROC_MET[iso] = 1.0 - yEffMET[iso];
        }

        TGraph *gROC_MET = new TGraph(nIso, xROC_MET.data(), yROC_MET.data());
        gROC_MET->SetName(TString("rocMET_") + tag);
        gROC_MET->Write();

        TGraph *gROC_comb = new TGraph(nIso, xROC_comb.data(), yROC_comb.data());
        gROC_comb->SetName(TString("rocBkg_") + tag);
        gROC_comb->Write();

        // ---- AUC ----
        const double aucSS = rocAUC(xROC_SS, yROC_SS);
        const double aucMET = rocAUC(xROC_MET, yROC_MET);

        TParameter<double> pAucSS(TString("aucSS_") + tag, aucSS);
        TParameter<double> pAucMET(TString("aucMET_") + tag, aucMET);

        pAucSS.Write();
        pAucMET.Write();

        double aucBkg = rocAUC(xROC_comb, yROC_comb);
        TParameter<double> pAucBkg(TString("aucBkg_") + tag, aucBkg);
        pAucBkg.Write();
    }

    fout->Write();
    fout->Close();

    std::cout << "Saved graphs to IsoStudyOutputs.root\n";
}