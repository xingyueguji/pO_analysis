// DrawDielectronPeak.C
// Example:
//   root -l -q 'DrawDielectronPeak.C("root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2025.root", kData)'
//
// Sample selection mirrors DrawWToElecNu_PFMet.C:
//   kData / kDY / kWp / kWm / kDYtau / kWptau / kWmtau

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TMath.h"
#include "TPaveText.h"
#include <vector>
#include <iostream>
#include <string>

enum SampleType
{
    kData,
    kDY,
    kWp,
    kWm,
    kDYtau,
    kWptau,
    kWmtau
};

// --------------------
// Small helpers
// --------------------
static double MU_MASS() { return 0.1056583745; }
static double ELE_MASS() { return 0.000511; }

static bool HasBranch(TTree *t, const char *name)
{
    return (t && t->GetListOfBranches() && t->GetListOfBranches()->FindObject(name));
}

static bool TriggerFired(int hltBit)
{
    return (hltBit != 0);
}

static std::string FindBranchContaining(TTree *t, const std::string &needle)
{
    if (!t || !t->GetListOfBranches())
        return "";
    TObjArray *arr = t->GetListOfBranches();
    const int n = arr->GetEntries();
    for (int i = 0; i < n; ++i)
    {
        TBranch *br = (TBranch *)arr->At(i);
        if (!br)
            continue;
        std::string bn = br->GetName();
        if (bn.find(needle) != std::string::npos)
            return bn;
    }
    return "";
}

static double RelIsoPF(int i,
                       std::vector<float> *elePt,
                       std::vector<float> *elePFChIso,
                       std::vector<float> *elePFNeuIso,
                       std::vector<float> *elePFPhoIso)
{
    if (!elePt || !elePFChIso || !elePFNeuIso || !elePFPhoIso)
        return 999.0;
    const double pt = elePt->at(i);
    if (pt <= 0)
        return 999.0;
    return (elePFChIso->at(i) + elePFNeuIso->at(i) + elePFPhoIso->at(i)) / pt;
}

static bool PassEventSelection_pO(TTree *tHi,
                                  bool &warnedOnce,
                                  bool has_pprimaryVertexFilter, int pprimaryVertexFilter,
                                  bool has_pclusterCompatibilityFilter, int pclusterCompatibilityFilter)
{
    (void)tHi;

    auto requireIfExists = [&](bool has, int val, const char *name) -> bool
    {
        if (!has)
        {
            if (!warnedOnce)
                std::cout << "[WARN] Missing event filter branch: " << name << " (will not apply)\n";
            return true;
        }
        return (val != 0);
    };

    bool ok = true;
    ok = ok && requireIfExists(has_pprimaryVertexFilter, pprimaryVertexFilter, "pprimaryVertexFilter");
    ok = ok && requireIfExists(has_pclusterCompatibilityFilter, pclusterCompatibilityFilter, "pclusterCompatibilityFilter");

    warnedOnce = true;
    return ok;
}

// --------------------
// User config struct
// --------------------
struct DileptonConfig
{
    double massMin = 60.0, massMax = 120.0;
    int    nBins   = 120;
    double ptMin1 = 20.0, ptMin2 = 10.0;
    double etaMax = 2.4;          // FIXME: ECAL is 2.5 with 1.4442–1.566 gap

    bool requireOS      = true;
    bool requireTightID = true;   // active gate uses eleCutIdWP95 — FIXME WP95 is loose

    bool   applyVz    = true;
    double vzMax      = 15.0;
    bool   applyHiBin = false;
    int    hiBinMin   = 0;
    int    hiBinMax   = 200;

    double isoMax = 0.095;          // FIXME: tune for electrons (PU correction missing)

    std::string outPrefix = "ZToEE_pO2025";
};

// --------------------
// Main
// --------------------
void DrawDielectronPeak(const char *fname = "root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2025.root",
                        SampleType sample = kData)
{
    DileptonConfig cfg;
    bool warnedEventFiltersOnce = false;

    // Z defaults
    cfg.massMin = 60;
    cfg.massMax = 120;
    cfg.nBins   = 120;
    cfg.ptMin1  = 15;
    cfg.ptMin2  = 10;
    cfg.requireTightID = true;
    cfg.applyVz   = true;
    cfg.vzMax     = 15;
    cfg.applyHiBin = false;
    cfg.outPrefix  = "ZToEE_pO2025";

    // -------------------------
    // Sample / file switching (mirrors DrawWToElecNu_PFMet.C)
    // -------------------------
    bool isMC = (sample != kData);
    if (isMC)
    {
        if (sample == kDY)    { fname = Form("root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_MC_DY_ele_Z.root"); cfg.outPrefix += "_DY"; }
        if (sample == kWp)    { fname = Form("root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_MC_Wp_ele.root");   cfg.outPrefix += "_Wp"; }
        if (sample == kWm)    { fname = Form("root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_MC_Wm_ele.root");   cfg.outPrefix += "_Wm"; }
        if (sample == kDYtau) { fname = Form("root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_MC_DY_tau_Z.root"); cfg.outPrefix += "_DYtau"; }
        if (sample == kWptau) { fname = Form("root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_MC_Wp_tau.root");   cfg.outPrefix += "_Wptau"; }
        if (sample == kWmtau) { fname = Form("root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_MC_Wm_tau.root");   cfg.outPrefix += "_Wmtau"; }
    }
    std::string mcTag = isMC ? "MC" : "Data";

    // --------------------
    // Open file / get trees
    // --------------------
    TFile *f = TFile::Open(fname);
    if (!f || f->IsZombie())
    {
        std::cerr << "ERROR: cannot open " << fname << "\n";
        return;
    }

    TTree *tEle = (TTree *)f->Get("ggHiNtuplizer/EventTree");
    if (!tEle)
    {
        std::cerr << "ERROR: missing ggHiNtuplizer/EventTree\n";
        return;
    }

    TTree *tHi = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
    bool haveHiTree = (tHi != nullptr);

    // --------------------
    // Electron kinematics
    // --------------------
    Int_t nEle = 0;
    std::vector<float> *elePt = nullptr, *eleEta = nullptr, *elePhi = nullptr;
    std::vector<int>   *eleCharge = nullptr;

    tEle->SetBranchAddress("nEle",      &nEle);
    tEle->SetBranchAddress("elePt",     &elePt);
    tEle->SetBranchAddress("eleEta",    &eleEta);
    tEle->SetBranchAddress("elePhi",    &elePhi);
    tEle->SetBranchAddress("eleCharge", &eleCharge);

    // --------------------
    // ID WPs (kept around; active gate uses eleCutIdWP95 below)
    // --------------------
    std::vector<int> *eleMVAIdWP80 = nullptr, *eleMVAIdWP85 = nullptr;
    std::vector<int> *eleMVAIdWP90 = nullptr, *eleMVAIdWP95 = nullptr;
    std::vector<int> *eleCutIdWP70 = nullptr, *eleCutIdWP80 = nullptr;
    std::vector<int> *eleCutIdWP90 = nullptr, *eleCutIdWP95 = nullptr;

    auto hookID = [&](const char *bname, std::vector<int> **addr)
    {
        if (HasBranch(tEle, bname))
        {
            tEle->SetBranchStatus(bname, 1);
            tEle->SetBranchAddress(bname, addr);
        }
    };
    hookID("eleMVAIdWP80", &eleMVAIdWP80);
    hookID("eleMVAIdWP85", &eleMVAIdWP85);
    hookID("eleMVAIdWP90", &eleMVAIdWP90);
    hookID("eleMVAIdWP95", &eleMVAIdWP95);
    hookID("eleCutIdWP70", &eleCutIdWP70);
    hookID("eleCutIdWP80", &eleCutIdWP80);
    hookID("eleCutIdWP90", &eleCutIdWP90);
    hookID("eleCutIdWP95", &eleCutIdWP95);

    const bool has_eleID = HasBranch(tEle, "eleCutIdWP95");

    // --------------------
    // Iso WPs (kept around; active gate uses continuous RelIsoPF below)
    // --------------------
    std::vector<int> *eleMVAIsoWP80 = nullptr, *eleMVAIsoWP85 = nullptr;
    std::vector<int> *eleMVAIsoWP90 = nullptr, *eleMVAIsoWP95 = nullptr;

    hookID("eleMVAIsoWP80", &eleMVAIsoWP80);
    hookID("eleMVAIsoWP85", &eleMVAIsoWP85);
    hookID("eleMVAIsoWP90", &eleMVAIsoWP90);
    hookID("eleMVAIsoWP95", &eleMVAIsoWP95);

    // Continuous PF iso (active path)
    std::vector<float> *elePFChIso = nullptr, *elePFNeuIso = nullptr, *elePFPhoIso = nullptr;
    const bool has_elePFChIso  = HasBranch(tEle, "elePFChIso");
    const bool has_elePFNeuIso = HasBranch(tEle, "elePFNeuIso");
    const bool has_elePFPhoIso = HasBranch(tEle, "elePFPhoIso");
    if (has_elePFChIso)  { tEle->SetBranchStatus("elePFChIso", 1);  tEle->SetBranchAddress("elePFChIso",  &elePFChIso);  }
    if (has_elePFNeuIso) { tEle->SetBranchStatus("elePFNeuIso", 1); tEle->SetBranchAddress("elePFNeuIso", &elePFNeuIso); }
    if (has_elePFPhoIso) { tEle->SetBranchStatus("elePFPhoIso", 1); tEle->SetBranchAddress("elePFPhoIso", &elePFPhoIso); }

    // --------------------
    // Optional event branches
    // --------------------
    Float_t vz    = 999.f;
    Int_t   hiBin = -999;

    const bool has_vz    = (haveHiTree && HasBranch(tHi, "vz"));
    const bool has_hiBin = (haveHiTree && HasBranch(tHi, "hiBin"));

    if (has_vz)    tHi->SetBranchAddress("vz",    &vz);
    if (has_hiBin) tHi->SetBranchAddress("hiBin", &hiBin);

    // --------------------
    // HLT objects
    // --------------------
    TTree *tHLTobj = (TTree *)f->Get("hltobject/HLT_OxyL1SingleEG10_v");

    std::vector<double> *trgObjPt  = nullptr;
    std::vector<double> *trgObjEta = nullptr;
    std::vector<double> *trgObjPhi = nullptr;
    std::vector<double> *trgObjId  = nullptr;

    const bool has_trgObjPt  = HasBranch(tHLTobj, "pt");
    const bool has_trgObjEta = HasBranch(tHLTobj, "eta");
    const bool has_trgObjPhi = HasBranch(tHLTobj, "phi");
    const bool has_trgObjId  = HasBranch(tHLTobj, "TriggerObjID");

    tHLTobj->SetBranchStatus("*", false);
    if (has_trgObjPt && has_trgObjEta && has_trgObjPhi && has_trgObjId)
    {
        tHLTobj->SetBranchStatus("pt", 1);
        tHLTobj->SetBranchStatus("eta", 1);
        tHLTobj->SetBranchStatus("phi", 1);
        tHLTobj->SetBranchStatus("TriggerObjID", 1);

        tHLTobj->SetBranchAddress("pt", &trgObjPt);
        tHLTobj->SetBranchAddress("eta", &trgObjEta);
        tHLTobj->SetBranchAddress("phi", &trgObjPhi);
        tHLTobj->SetBranchAddress("TriggerObjID", &trgObjId);
    }
    else
    {
        cout << "[WARN] Could not find HLT related Object." << endl;
    }

    // --------------------
    // HLT bit
    // --------------------
    TTree *tHLT = (TTree *)f->Get("hltanalysis/HltTree");

    const std::string hltName = FindBranchContaining(tHLT, "HLT_OxyL1SingleEG10_v1");
    Int_t HLT_OxyL1SingleEG10_v1 = 0;
    bool has_hlt = false;
    tHLT->SetBranchStatus("*", 0);

    if (!hltName.empty())
    {
        has_hlt = true;
        tHLT->SetBranchStatus(hltName.c_str(), 1);
        tHLT->SetBranchAddress(hltName.c_str(), &HLT_OxyL1SingleEG10_v1);
        std::cout << "[INFO] Using trigger bit branch: " << hltName << "\n";
    }
    else
    {
        std::cout << "[WARN] Could not find any branch containing 'HLT_OxyL1SingleEG10_v1' in hltanalysis/HltTree.\n"
                  << "       Trigger gate will be treated as PASS.\n";
    }

    // --------------------
    // Event filters
    // --------------------
    TTree *tEvent = (TTree *)f->Get("skimanalysis/HltTree");

    Int_t pprimaryVertexFilter = 1;
    Int_t pclusterCompatibilityFilter = 1;

    const bool has_pprimaryVertexFilter        = HasBranch(tEvent, "pprimaryVertexFilter");
    const bool has_pclusterCompatibilityFilter = HasBranch(tEvent, "pclusterCompatibilityFilter");

    tEvent->SetBranchStatus("*", 0);

    if (has_pprimaryVertexFilter && has_pclusterCompatibilityFilter)
    {
        tEvent->SetBranchStatus("pprimaryVertexFilter", 1);
        tEvent->SetBranchStatus("pclusterCompatibilityFilter", 1);

        tEvent->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter);
        tEvent->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);
    }

    // --------------------
    // Histograms
    // --------------------
    gStyle->SetOptStat(0);
    TH1D *hMass = new TH1D("hMass",
                           Form("%s; m_{ee} [GeV]; Events", cfg.outPrefix.c_str()),
                           cfg.nBins, cfg.massMin, cfg.massMax);

    TH1D *hMass_extended = new TH1D("hMass_extended",
                                    Form("%s; m_{ee} [GeV]; Events", cfg.outPrefix.c_str()),
                                    cfg.nBins, 15, 120);

    TH1D *hMass_vipul = new TH1D("hMass_vipul",
                                 Form("%s; m_{ee} [GeV]; Events", cfg.outPrefix.c_str()),
                                 cfg.nBins, cfg.massMin, cfg.massMax);
    hMass->Sumw2();
    hMass_extended->Sumw2();
    hMass_vipul->Sumw2();

    // --------------------
    // Loop
    // --------------------
    const Long64_t nEntries = tEle->GetEntries();
    std::cout << "EventTree entries = " << nEntries << "\n";
    std::cout << "Have HiTree? " << (haveHiTree ? "YES" : "NO") << "\n";
    if (haveHiTree)
    {
        std::cout << "HiTree has vz? " << (has_vz ? "YES" : "NO") << "\n";
        std::cout << "HiTree has hiBin? " << (has_hiBin ? "YES" : "NO") << "\n";
    }

    Long64_t nPassEvent = 0;
    Long64_t nPassPair  = 0;

    for (Long64_t ie = 0; ie < nEntries; ++ie)
    {
        if (ie % 200000 == 0)
            std::cout << "Event " << ie << "/" << nEntries << "\n";

        tEle->GetEntry(ie);
        if (haveHiTree) tHi->GetEntry(ie);
        tHLT->GetEntry(ie);
        tHLTobj->GetEntry(ie);
        tEvent->GetEntry(ie);

        if (cfg.applyVz && has_vz)
        {
            if (TMath::Abs(vz) > cfg.vzMax)
                continue;
        }

        if (!PassEventSelection_pO(tHi,
                                   warnedEventFiltersOnce,
                                   has_pprimaryVertexFilter, pprimaryVertexFilter,
                                   has_pclusterCompatibilityFilter, pclusterCompatibilityFilter))
            continue;

        if (cfg.applyHiBin && has_hiBin)
        {
            if (hiBin < cfg.hiBinMin || hiBin > cfg.hiBinMax)
                continue;
        }

        if (nEle < 2)
            continue;

        if (has_hlt)
        {
            if (!TriggerFired(HLT_OxyL1SingleEG10_v1))
                continue;
        }
        nPassEvent++;

        for (int i = 0; i < nEle; ++i)
        {
            bool passIsolead = true;

            if (!elePt || !eleEta || !elePhi || !eleCharge)
                continue;

            if (elePt->at(i) < cfg.ptMin2)
                continue;
            if (TMath::Abs(eleEta->at(i)) > cfg.etaMax)
                continue;

            if (cfg.requireTightID && has_eleID && eleCutIdWP95->at(i) == 0)
                continue;

            /* WP-based iso, kept for cross-check
            if (eleMVAIsoWP95->at(i) != 1) passIsolead = false;
            */
            const double isoLead = RelIsoPF(i, elePt, elePFChIso, elePFNeuIso, elePFPhoIso);
            if (isoLead >= cfg.isoMax)
                passIsolead = false;
            // FIXME: redo electron isolation study; not using a working point.

            for (int j = i + 1; j < nEle; ++j)
            {
                bool passIsosec = true;

                if (elePt->at(j) < cfg.ptMin2)
                    continue;
                if (TMath::Abs(eleEta->at(j)) > cfg.etaMax)
                    continue;

                if (cfg.requireTightID && has_eleID && eleCutIdWP95->at(j) == 0)
                    continue;

                if (cfg.requireOS)
                {
                    if (eleCharge->at(i) * eleCharge->at(j) >= 0)
                        continue;
                }

                const double pt1  = elePt->at(i);
                const double pt2  = elePt->at(j);
                const double lead = (pt1 > pt2 ? pt1 : pt2);
                const double sub  = (pt1 > pt2 ? pt2 : pt1);
                if (lead < cfg.ptMin1) continue;
                if (sub  < cfg.ptMin2) continue;

                TLorentzVector v1, v2;
                v1.SetPtEtaPhiM(elePt->at(i), eleEta->at(i), elePhi->at(i), ELE_MASS());
                v2.SetPtEtaPhiM(elePt->at(j), eleEta->at(j), elePhi->at(j), ELE_MASS());
                const double m = (v1 + v2).M();

                /* WP-based iso, kept for cross-check
                if (eleMVAIsoWP95->at(j) != 1) passIsosec = false;
                */
                const double isosub = RelIsoPF(j, elePt, elePFChIso, elePFNeuIso, elePFPhoIso);
                if (isosub >= cfg.isoMax)
                    passIsosec = false;

                if (passIsolead && passIsosec)
                {
                    hMass_extended->Fill(m);
                    if (!(m < 60 || m > 120))
                    {
                        hMass->Fill(m);
                    }
                }

                if (m < 60 || m > 120) continue;
                if (lead < 20 || sub < 20) continue;

                // Tighter pT cut, no iso requirement.
                hMass_vipul->Fill(m);
                nPassPair++;
            }
        }
    }

    std::cout << "Passed event preselection: " << nPassEvent << "\n";
    std::cout << "Found selected OS pair: "    << nPassPair  << "\n";

    // --------------------
    // Save
    // --------------------
    TFile *fout = new TFile(("./rootfile/" + cfg.outPrefix + "_" + mcTag + "_hist.root").c_str(), "RECREATE");
    hMass->Write("", 2);
    hMass_extended->Write("", 2);
    hMass_vipul->Write("", 2);
    fout->Close();

    std::cout << "[INFO] Wrote outputs with prefix: " << cfg.outPrefix << "\n";
}