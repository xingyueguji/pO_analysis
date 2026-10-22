// correction/tnp_charge_ele.C
//
// THE ELECTRON TAG-AND-PROBE EFFICIENCIES SPLIT BY PROBE CHARGE (2026-10-20,
// user request: "split the Z->ee tag-and-probe (reco and ID+iso) by probe
// charge in |eta_SC| > 1.6 vs < 1.6 and report whether e- loses efficiency at
// high |eta| in data").
//
// WHY. The W- -> e nu yield sits ~25% below the muon one at the lab edges
// (|eta| 2.0-2.4 on BOTH sides, y = -eta), the W+ one does not:
// correction/logs/open_issues_diag.log, PF MET >= 30 with NO multijet
// subtracted, e/mu of (data - other MC)/W = 0.755 +- 0.086 (W-) vs 0.969 +-
// 0.088 (W+). Every efficiency the analysis applies or checks is
// charge-inclusive, so a loss of e- alone would go unnoticed.
//
// WHAT. The two electron tag-and-probe links of the analysis, with the probe
// split into six cells: charge (-, +) x {|eta_SC| < 1.6, eta_SC <= -1.6,
// eta_SC >= +1.6}. In the analysis convention y = -eta (proton-going =
// forward), eta_SC <= -1.6 is the FORWARD edge (lab y >= 1.6, bins y10-y11) and
// eta_SC >= +1.6 the BACKWARD one (y <= -1.6, bins y0-y1); "|eta_SC| > 1.6" =
// both (2026-10-20: the forward/backward split was added, user request).
//   reco  = correction/reco_sf_ele.C, its nominal probe: a photon object = a
//           supercluster, E_T^SC > 20, NO ID; pass = any forest electron on the
//           SC (DR_SC < 0.05). A failing probe has no electron, so the probe
//           charge is MINUS THE TAG'S (the pair is a Z), for pass and fail alike.
//           The probe four-vector in m(tag, probe) is reco_sf_ele.C's
//           RecoProbeP4 = EGM's (SC energy, direction from the origin, electron
//           mass) since 2026-10-20.
//   recopho = the same probes with the PREVIOUS four-vector, the photon
//           object's (phoEt, phoEta, phoPhi, massless: the calibrated energy,
//           as the electron's) -- kept side by side with reco since 2026-10-20
//           (user request), because the SC energy is not calibrated in pO data
//           (see CLAUDE.md, reco_sf_ele.C). Only m differs, so only which pairs
//           fall in 60-120 / 80-100 and the fitted shapes change. Its cells add
//           up to the pre-SC reco_sf_ele.C skim ("[CHECK]" against
//           backups/pre_scmass/correction/rootfile/reco_sf_ele/).
//   idiso = the Z tag-and-probe of correction/idiso_sf_skim.C: leg = an
//           electron with pT > 25, |eta_SC| < 2.4 outside the crack, relIso <
//           0.3; pass = eleMVAIdWP90 && eleMVAIsoWP90 (= the analysis
//           electron); tag = a passing leg matched to an EG10 object; one pair
//           per event (>= 1 tag, 60 <= m <= 120, closest to the Z mass), a leg
//           is a probe when its partner is a tag. The probe's OWN charge. The
//           same-sign pairs are filled too: the cut-and-count background, and a
//           data-vs-MC charge-flip cross-check.
// Tags, probes, pass definitions and the pair choice are copied from the two
// studies (the tag is the same in both), so the four cells add up to their
// inclusive numbers; the fit stage checks it ("[CHECK]") against
// electron_sf_check/reco_sf/rootfile/skim_<S>.root (h_m_noid_incl0_*, bin by
// bin) and electron_sf_check/idiso_sf/rootfile/skim_<S>.root (h_tnp_probes),
// i.e. the parent studies' kOutDir / kRootDir. No scale factor is
// applied to the MC, only the generator weight (as in both studies).
// The trigger link is not repeated: correction/trig_eff_mb.C measures it on the
// W sample itself, data efficiency 0.97-0.98 in both |y| > 1.6 bins and W+
// 0.9745 / W- 0.9740 inclusive (m_T > 40), which leaves no room for a ~25%
// loss of one charge.
//
// EXTRACTION per bin (inclusive, per charge, per charge x |eta| half, per charge
// x side at |eta| > 1.6):
//   fit = reco_sf_ele.C's RecoFitEps: a simultaneous binned ML fit of the pass
//         and fail m spectra, 60-120 GeV, the DY-MC templates of the same bin
//         convolved with one Gaussian (floated in the inclusive fit, fixed to
//         it in the others) + one exponential background per category; MINOS
//         on eps.
//   cc  = cut-and-count: reco 80-100 GeV with the fitted fail background in the
//         window subtracted (reco_sf_ele.C's [CC]); idiso OS minus SS over
//         60-120, data and MC alike.
//   raw = no background subtraction at all (OS 60-120): every failing probe is
//         counted as a lost Z electron, so this is (up to the ~2% background of
//         the pass category) the LOWEST efficiency the data allow -- the robust
//         number where a fit does not converge (few fails: eps runs into 1).
//   eps_MC is counted from the same DY-MC histograms (fit: OS, 60-120; cc: the
//   cc recipe; raw: OS 60-120), SF = eps_data / eps_MC.
//   ERRORS (2026-10-20, user request): the cc and raw data efficiencies are
//   68% PROFILE-LIKELIHOOD intervals capped at eps = 1 (EffInterval: the
//   signal profiled, the reco fail background known, the idiso SS counts a
//   measured background), so asymmetric -- one-sided [eps_lo, 1] when no fail
//   is left -- replacing the linear (Wald) propagation and its 1.841
//   zero-count patch; central values unchanged. eps_MC stays linear (~1e-4).
//   Asymmetric errors are combined per side in the SF, the e-/e+ ratios and
//   the reco x idiso products. CAVEAT, measured with toys (2026-10-20): with
//   ~35 probes the likelihood interval covers 52% at a true eps = 0.97 (36% of
//   the toys have no fail and get [0.986, 1]), 90% at 0.99, 67-68% at 0.95 /
//   0.90 and with 200 probes -- the reco side cells sit in the under-covering
//   regime.
// THE ANSWER is SF(e-) / SF(e+) per |eta| half and per side ("[RATIO]"): the e-
// and e+ probes are disjoint, so their errors are independent (the MINOS errors
// are symmetrized for the ratio). Its value at |eta_SC| > 1.6 is to be compared
// with the per-electron e-/e+ factor the W data would need.
//
// OUTPUTS live in correction/electron_sf_check/charge_split/{rootfile,plots,logs}/
// (since 2026-10-20; until then rootfile/tnp_charge_ele/, plots/tnp_charge_ele/, logs/).
// RUN: ./run_tnp_charge_ele.sh [skim|fit|all]  (<outputs>/logs/tnp_charge_ele_*.log = THE record)
//   root -l -b -q 'tnp_charge_ele.C+("skim","Data")' -> <outputs>/rootfile/skim_Data.root
//   root -l -b -q 'tnp_charge_ele.C+("skim","DY")'   -> <outputs>/rootfile/skim_DY.root
//   root -l -b -q 'tnp_charge_ele.C+("fit")'         -> <outputs>/plots/fit_<link>_<bin>_{pass,fail},
//                                                       sf_charge_<link> (vs |eta_SC|), sf_charge_signed_<link>
//                                                       (vs signed eta_SC) (.png/.pdf) + tnp_charge_ele.csv
// The two reco plots (sf_charge[_signed]_reco) show THREE cases (user request
// 2026-10-20): the cut-and-count SF with the SC-energy probe (filled) and with
// the photon-object probe (open), per charge, and the pp EGM 2025 Reco SF
// (Reco20to75, charge-inclusive) in ITS OWN (eta_SC, pT) cells -- a box per
// eta bin per pT bin, not the average over our probes (the "[EGMPRED]" of
// reco_sf_ele.C). The idiso plots keep fit + cut-and-count.
// Run from correction/. This file #includes reco_sf_ele.C (the definitions and
// the fit are single-sourced there): never load its library and
// reco_sf_ele_C.so in the same ROOT session.

#include "reco_sf_ele.C"
#include "idiso_sf_common.h"

#include "TBox.h"
#include "TGraphAsymmErrors.h"

#include <tuple>

namespace pOTnpQ
{
const double kEtaSplit = 1.6;
const char *const kOutDir  = "electron_sf_check/charge_split/rootfile";
const char *const kPlotDir = "electron_sf_check/charge_split/plots";

enum Link { kReco = 0, kRecoPho, kIdIso, kNLink };
const char *const kLinkName[kNLink]   = {"reco", "recopho", "idiso"};
const char *const kLinkLabel[kNLink]  = {"electron reco (SC #rightarrow electron), probe p4 = SC energy (EGM)",
                                         "electron reco (SC #rightarrow electron), probe p4 = photon object (previous)",
                                         "WP90 ID+iso (electron #rightarrow analysis e)"};
const char *const kLinkHeader[kNLink] = {"Z#rightarrowee tag-and-probe, electron reco (SC energy)",
                                         "Z#rightarrowee tag-and-probe, electron reco (photon p4)",
                                         "Z#rightarrowee tag-and-probe, WP90 ID+iso"};
const char *const kLinkShort[kNLink]  = {"electron reco", "electron reco (photon p4)", "WP90 ID+iso"};
inline bool IsReco(int link) { return link == kReco || link == kRecoPho; }
// the pre-2026-10-20 reco_sf_ele.C skim (photon-object probe p4): the [CHECK] reference of recopho
const char *const kPhoRefDir = "../backups/pre_scmass/correction/rootfile/reco_sf_ele";

// the six cells: probe charge x {|eta_SC| < 1.6, forward eta_SC <= -1.6,
// backward eta_SC >= +1.6} (y = -eta: forward = proton-going = lab y > 0)
const int         kNCell = 6;
const char *const kCellName[kNCell]  = {"em_lo", "em_fwd", "em_bwd", "ep_lo", "ep_fwd", "ep_bwd"};
const char *const kCellLabel[kNCell] = {"e^{-}, |#eta_{SC}| < 1.6", "e^{-}, fwd #eta_{SC} < -1.6", "e^{-}, bwd #eta_{SC} > 1.6",
                                        "e^{+}, |#eta_{SC}| < 1.6", "e^{+}, fwd #eta_{SC} < -1.6", "e^{+}, bwd #eta_{SC} > 1.6"};
inline int CellOf(int q, double etaSC)
{
  const int base = (q < 0) ? 0 : 3;
  if (std::fabs(etaSC) < kEtaSplit) return base;
  return base + (etaSC < 0.0 ? 1 : 2);
}

// the reported bins = unions of cells (the first seven = the original |eta|-half bins)
struct Bin
{
  const char *name, *label;
  int         n, cells[kNCell];
};
const int kNBin = 11;
const Bin kBins[kNBin] = {
    {"incl", "inclusive", 6, {0, 1, 2, 3, 4, 5}},
    {"em", "e^{-}", 3, {0, 1, 2}},
    {"ep", "e^{+}", 3, {3, 4, 5}},
    {"em_lo", "e^{-}, |#eta_{SC}| < 1.6", 1, {0}},
    {"em_hi", "e^{-}, |#eta_{SC}| > 1.6", 2, {1, 2}},
    {"ep_lo", "e^{+}, |#eta_{SC}| < 1.6", 1, {3}},
    {"ep_hi", "e^{+}, |#eta_{SC}| > 1.6", 2, {4, 5}},
    {"em_fwd", "e^{-}, fwd #eta_{SC} < -1.6", 1, {1}},
    {"em_bwd", "e^{-}, bwd #eta_{SC} > 1.6", 1, {2}},
    {"ep_fwd", "e^{+}, fwd #eta_{SC} < -1.6", 1, {4}},
    {"ep_bwd", "e^{+}, bwd #eta_{SC} > 1.6", 1, {5}},
};
const int kBinEm = 1, kBinEp = 2, kBinEmLo = 3, kBinEmHi = 4, kBinEpLo = 5, kBinEpHi = 6;
const int kBinEmFwd = 7, kBinEmBwd = 8, kBinEpFwd = 9, kBinEpBwd = 10;

inline TString HName(int link, int cell, bool pass, bool ss = false)
{
  return Form("h_m_%s_%s_%s%s", kLinkName[link], kCellName[cell], pass ? "pass" : "fail", ss ? "_ss" : "");
}
// the passing probes' electron pT per cell (reco: the SC-energy link, 60-120
// GeV; idiso: OS): which pT cells of the EGM tables our probes populate
inline TString PtName(int link, int cell) { return Form("h_pt_%s_%s_pass", kLinkName[link], kCellName[cell]); }
} // namespace pOTnpQ

// ============================================================================
// SKIM: one ntuple (Data or DY->ee MC) -> m(tag, probe) per link x cell x
// {pass, fail} (+ same-sign twins for idiso)
// ============================================================================
int TnpQSkim(const std::string &sampleName, Long64_t nmax)
{
  using namespace pOTnpQ;
  namespace R = pORecoSF;
  namespace I = pOIdIso;
  const bool isMC = (sampleName != "Data");
  std::string fname;
  if (!isMC)
    fname = pOSkim::kDefaultDataFile;
  else if (sampleName == "DY")
    fname = pOSkim::ResolveMCSample(kDY, "ele").fname; // SampleType is global (skim_common.h)
  else
  {
    std::cerr << "[FATAL] sample must be Data or DY, got " << sampleName << "\n";
    return 2;
  }
  std::cout << "[INPUT] " << sampleName << " <- " << fname << "\n";
  std::cout << "[CONFIG] cells: probe charge x |eta_SC| < / >= " << kEtaSplit << "\n"
            << "  reco : = reco_sf_ele.C nominal (tag WP90 ID+iso, pT >= " << R::kTagPtMin << ", |eta_SC| < " << R::kEtaSCMax
            << " no crack, relIso < " << pOSkim::kEleIsoRelMax << ", EG10 DR < " << R::kTrigMatchDR << " | probe photon SC E_T > "
            << R::kProbeEtMin << ", NO ID | pass: electron DR_SC < " << R::kPassDR << "), probe charge = -tag charge\n"
            << "  idiso: = idiso_sf_skim.C Z TnP (legs pT > " << I::kPtMin << ", relIso < " << I::kIsoTotalMax
            << ", pass " << I::kPassIdBranch << " && " << I::kPassIsoBranch << ", one OS (+ one SS) pair per event),"
            << " probe charge = its own\n";

  TFile *f = TFile::Open(fname.c_str());
  if (!f || f->IsZombie()) { std::cerr << "[FATAL] cannot open " << fname << "\n"; return 2; }
  TTree *tEle = (TTree *)f->Get("ggHiNtuplizer/EventTree");
  TTree *tHi  = (TTree *)f->Get("hiEvtAnalyzer/HiTree");
  TTree *tHLT = (TTree *)f->Get("hltanalysis/HltTree");
  TTree *tObj = (TTree *)f->Get(R::kTrigObjTree);
  TTree *tEvt = (TTree *)f->Get("skimanalysis/HltTree");
  if (!tEle || !tHi || !tHLT || !tObj || !tEvt) { std::cerr << "[FATAL] missing a required tree in " << fname << "\n"; return 2; }
  for (TTree *t : {tHi, tHLT, tObj, tEvt})
    if (t->GetEntries() != tEle->GetEntries())
    {
      std::cerr << "[FATAL] tree " << t->GetName() << " (" << t->GetEntries() << ") and EventTree ("
                << tEle->GetEntries() << ") entry counts differ\n";
      return 2;
    }

  // -------- EventTree: every branch read is enabled AND addressed --------
  std::vector<float> *elePt = nullptr, *eleEta = nullptr, *elePhi = nullptr, *eleSCEta = nullptr, *eleSCPhi = nullptr;
  std::vector<int>   *eleCharge = nullptr, *eleId = nullptr, *eleIso = nullptr;
  std::vector<float> *chIso = nullptr, *neuIso = nullptr, *phIso = nullptr, *puIso = nullptr;
  std::vector<float> *phoEt = nullptr, *phoEta = nullptr, *phoPhi = nullptr, *phoSCE = nullptr, *phoSCEta = nullptr,
                     *phoSCPhi = nullptr;
  const std::vector<const char *> need = {"elePt", "eleEta", "elePhi", "eleSCEta", "eleSCPhi", "eleCharge",
                                          pOSkim::kEleIDBranch, pOSkim::kEleIsoBranch, "elePFChIso", "elePFNeuIso",
                                          "elePFPhoIso", "elePFPUIso", "phoEt", "phoEta", "phoPhi", "phoSCE",
                                          "phoSCEta", "phoSCPhi"};
  for (const char *b : need)
    if (!pOSkim::HasBranch(tEle, b)) { std::cerr << "[FATAL] missing EventTree branch " << b << "\n"; return 2; }
  tEle->SetBranchStatus("*", 0);
  auto bind = [&](const char *b, auto *addr) { tEle->SetBranchStatus(b, 1); tEle->SetBranchAddress(b, addr); };
  bind("elePt", &elePt);           bind("eleEta", &eleEta);        bind("elePhi", &elePhi);
  bind("eleSCEta", &eleSCEta);     bind("eleSCPhi", &eleSCPhi);    bind("eleCharge", &eleCharge);
  bind(pOSkim::kEleIDBranch, &eleId);  bind(pOSkim::kEleIsoBranch, &eleIso);
  bind("elePFChIso", &chIso);      bind("elePFNeuIso", &neuIso);   bind("elePFPhoIso", &phIso);
  bind("elePFPUIso", &puIso);
  bind("phoEt", &phoEt);           bind("phoEta", &phoEta);        bind("phoPhi", &phoPhi);
  bind("phoSCE", &phoSCE);         bind("phoSCEta", &phoSCEta);    bind("phoSCPhi", &phoSCPhi);

  // -------- trigger bit (exact name: FindBranchContaining may return a _Prescale twin) --------
  const std::string hltName = pOSkim::HasBranch(tHLT, R::kTrigPath) ? std::string(R::kTrigPath)
                                                                    : pOSkim::FindBranchContaining(tHLT, R::kTrigPath);
  if (hltName.empty()) { std::cerr << "[FATAL] hltanalysis/HltTree lacks " << R::kTrigPath << "\n"; return 2; }
  Int_t hltBit = 0;
  tHLT->SetBranchStatus("*", 0);
  tHLT->SetBranchStatus(hltName.c_str(), 1); tHLT->SetBranchAddress(hltName.c_str(), &hltBit);

  // -------- trigger objects --------
  std::vector<double> *toPt = nullptr, *toEta = nullptr, *toPhi = nullptr;
  for (const char *b : {"pt", "eta", "phi"})
    if (!pOSkim::HasBranch(tObj, b)) { std::cerr << "[FATAL] no " << b << " in " << R::kTrigObjTree << "\n"; return 2; }
  tObj->SetBranchStatus("*", 0);
  tObj->SetBranchStatus("pt", 1);  tObj->SetBranchAddress("pt", &toPt);
  tObj->SetBranchStatus("eta", 1); tObj->SetBranchAddress("eta", &toEta);
  tObj->SetBranchStatus("phi", 1); tObj->SetBranchAddress("phi", &toPhi);

  // -------- event filters, vz, gen weight --------
  Int_t ppv = 1, pcc = 1;
  const bool has_ppv = pOSkim::HasBranch(tEvt, "pprimaryVertexFilter");
  const bool has_pcc = pOSkim::HasBranch(tEvt, "pclusterCompatibilityFilter");
  tEvt->SetBranchStatus("*", 0);
  if (has_ppv) { tEvt->SetBranchStatus("pprimaryVertexFilter", 1);        tEvt->SetBranchAddress("pprimaryVertexFilter", &ppv); }
  if (has_pcc) { tEvt->SetBranchStatus("pclusterCompatibilityFilter", 1); tEvt->SetBranchAddress("pclusterCompatibilityFilter", &pcc); }
  Float_t vz = 0.f, genWeight = 1.f;
  if (!pOSkim::HasBranch(tHi, "vz")) { std::cerr << "[FATAL] HiTree lacks vz\n"; return 2; }
  tHi->SetBranchStatus("*", 0);
  tHi->SetBranchStatus("vz", 1); tHi->SetBranchAddress("vz", &vz);
  const bool has_w = isMC && pOSkim::HasBranch(tHi, "weight");
  if (has_w) { tHi->SetBranchStatus("weight", 1); tHi->SetBranchAddress("weight", &genWeight); }
  else if (isMC) std::cout << "[WARN] no HiTree 'weight' -- MC filled unweighted.\n";

  // -------- output --------
  gSystem->mkdir(kOutDir, kTRUE);
  const TString outName = Form("%s/skim_%s.root", kOutDir, sampleName.c_str());
  TFile *fo = TFile::Open(outName, "RECREATE");
  if (!fo || fo->IsZombie()) { std::cerr << "[FATAL] cannot write " << outName << "\n"; return 2; }
  std::map<TString, TH1D *> h;
  for (int l = 0; l < kNLink; ++l)
    for (int c = 0; c < kNCell; ++c)
      for (int p = 0; p < 2; ++p)
        for (int ss = 0; ss < (l == kIdIso ? 2 : 1); ++ss)
        {
          const TString n = HName(l, c, p == 1, ss == 1);
          h[n] = new TH1D(n, ";m(tag, probe) [GeV];pairs / GeV", R::kNM, R::kMLo, R::kMHi);
          h[n]->Sumw2();
        }
  for (int l : {(int)kReco, (int)kIdIso})
    for (int c = 0; c < kNCell; ++c)
    {
      const TString n = PtName(l, c);
      h[n] = new TH1D(n, ";passing probe electron p_{T} [GeV];probes / GeV", 500, 0.0, 500.0);
      h[n]->Sumw2();
    }

  // -------- loop --------
  const Long64_t nEntries = tEle->GetEntries();
  const Long64_t nRun     = (nmax > 0 && nmax < nEntries) ? nmax : nEntries;
  std::cout << "Processing entries: " << nRun << " of " << nEntries << "\n";
  bool warnedFilters = false, warnedTrig = false;
  for (Long64_t ie = 0; ie < nRun; ++ie)
  {
    if (ie % 1000000 == 0) std::cout << "  event " << ie << "/" << nRun << "\n";
    tEle->GetEntry(ie);
    tHi->GetEntry(ie);
    tEvt->GetEntry(ie);
    tHLT->GetEntry(ie);
    const double w = has_w ? (double)genWeight : 1.0;
    if (!elePt || !eleEta || !elePhi || !eleSCEta || !eleSCPhi || !eleCharge || !eleId || !eleIso || !phoEt ||
        !phoEta || !phoPhi || !phoSCE || !phoSCEta || !phoSCPhi)
      continue;
    if (!pOSkim::PassEventSelection_pO(warnedFilters, has_ppv, ppv, has_pcc, pcc)) continue;
    if (std::fabs(vz) > R::kVzMax) continue;
    if (!pOSkim::TriggerFired(hltBit)) continue;

    const int ne = (int)elePt->size();
    tObj->GetEntry(ie);
    std::vector<double> iso(ne);
    for (int i = 0; i < ne; ++i) iso[i] = pOSkim::RelIsoPF(i, elePt, chIso, neuIso, phIso, puIso);
    auto passes  = [&](int i) { return eleId->at(i) != 0 && eleIso->at(i) != 0; };
    auto matched = [&](int i)
    {
      return pOSkim::PassLeadingLeptonTrigMatch(R::kTrigMatchDR, i, eleEta, elePhi, true, toPt, true, toEta, true, toPhi,
                                                warnedTrig, "electron");
    };

    // ---- reco links: reco_sf_ele.C's skim (nominal probe), cell from -(tag
    //      charge); one probe list, two four-vectors (reco = SC energy,
    //      recopho = the photon object, the pre-2026-10-20 definition)
    std::vector<int> tags;
    for (int i = 0; i < ne; ++i)
    {
      if (elePt->at(i) < R::kTagPtMin) continue;
      const double se = eleSCEta->at(i);
      if (std::fabs(se) >= R::kEtaSCMax || pOSkim::InEcalGap(se)) continue;
      if (!passes(i)) continue;
      if (iso[i] >= pOSkim::kEleIsoRelMax) continue;
      if (!matched(i)) continue;
      tags.push_back(i);
    }
    const int np = (int)phoEt->size();
    for (int it : tags)
    {
      const int qProbe = -eleCharge->at(it);
      TLorentzVector vt;
      vt.SetPtEtaPhiM(elePt->at(it), eleEta->at(it), elePhi->at(it), pOSkim::ELE_MASS);
      std::vector<int> used;
      for (int k = 0; k < np; ++k)
      {
        const double se = phoSCEta->at(k), sp = phoSCPhi->at(k), ase = std::fabs(se);
        const double et = phoSCE->at(k) / std::cosh(se);
        if (et < R::kProbeEtMin || ase >= R::kEtaSCMax || pOSkim::InEcalGap(se)) continue;
        if (R::RecoDRsc(se, sp, eleSCEta->at(it), eleSCPhi->at(it)) < R::kProbeTagDR) continue;
        bool dup = false;
        for (int u : used)
          if (R::RecoDRsc(se, sp, phoSCEta->at(u), phoSCPhi->at(u)) < R::kDupDR) { dup = true; break; }
        if (dup) continue;
        used.push_back(k);
        const double m = (vt + R::RecoProbeP4(phoSCE->at(k), se, sp)).M();
        TLorentzVector vp;
        vp.SetPtEtaPhiM(phoEt->at(k), phoEta->at(k), phoPhi->at(k), 0.0);
        const double mPho = (vt + vp).M();
        const bool   inSC = (m >= R::kMLo && m < R::kMHi), inPho = (mPho >= R::kMLo && mPho < R::kMHi);
        if (!inSC && !inPho) continue;
        int jPass = -1;
        for (int j = 0; j < ne && jPass < 0; ++j)
          if (j != it && R::RecoDRsc(eleSCEta->at(j), eleSCPhi->at(j), se, sp) < R::kPassDR) jPass = j;
        const bool pass = (jPass >= 0);
        const int  cell = CellOf(qProbe, se);
        if (inSC)
        {
          h[HName(kReco, cell, pass)]->Fill(m, w);
          if (pass) h[PtName(kReco, cell)]->Fill(elePt->at(jPass), w);
        }
        if (inPho) h[HName(kRecoPho, cell, pass)]->Fill(mPho, w);
      }
    }

    // ---- idiso link: idiso_sf_skim.C's Z tag-and-probe, per probe, cell from the probe's charge
    std::vector<int> leg;
    for (int i = 0; i < ne; ++i)
      if (elePt->at(i) > I::kPtMin && std::fabs(eleSCEta->at(i)) < I::kEtaSCMax && !pOSkim::InEcalGap(eleSCEta->at(i)) &&
          iso[i] < I::kIsoTotalMax)
        leg.push_back(i);
    if (leg.size() >= 2)
    {
      std::vector<char> pas(leg.size(), 0), tag(leg.size(), 0);
      for (size_t a = 0; a < leg.size(); ++a)
      {
        pas[a] = passes(leg[a]);
        tag[a] = pas[a] && matched(leg[a]);
      }
      for (int ss = 0; ss < 2; ++ss)
      {
        int ba = -1, bb = -1;
        double bm = 0.0, bd = 1e9;
        for (size_t a = 0; a < leg.size(); ++a)
          for (size_t b = a + 1; b < leg.size(); ++b)
          {
            const bool os = eleCharge->at(leg[a]) * eleCharge->at(leg[b]) < 0;
            if (os != (ss == 0)) continue;
            if (!tag[a] && !tag[b]) continue;
            TLorentzVector va, vb;
            va.SetPtEtaPhiM(elePt->at(leg[a]), eleEta->at(leg[a]), elePhi->at(leg[a]), pOSkim::ELE_MASS);
            vb.SetPtEtaPhiM(elePt->at(leg[b]), eleEta->at(leg[b]), elePhi->at(leg[b]), pOSkim::ELE_MASS);
            const double m = (va + vb).M();
            if (m < I::kMllLo || m > I::kMllHi) continue;
            if (std::fabs(m - I::kZMass) < bd) { bd = std::fabs(m - I::kZMass); ba = (int)a; bb = (int)b; bm = m; }
          }
        if (ba < 0) continue;
        for (int side = 0; side < 2; ++side)
        {
          const int self = side ? bb : ba, other = side ? ba : bb;
          if (!tag[other]) continue;
          const int ip = leg[self], cell = CellOf(eleCharge->at(ip), eleSCEta->at(ip));
          h[HName(kIdIso, cell, pas[self] != 0, ss == 1)]->Fill(bm, w);
          if (ss == 0 && pas[self]) h[PtName(kIdIso, cell)]->Fill(elePt->at(ip), w);
        }
      }
    }
  }

  // -------- report: raw counts per cell --------
  std::cout << "[RAW] " << sampleName << ": pass / fail per cell (60-120 GeV, weighted; idiso SS in brackets)\n";
  for (int l = 0; l < kNLink; ++l)
    for (int c = 0; c < kNCell; ++c)
    {
      const double p = h[HName(l, c, true)]->Integral(0, R::kNM + 1), fl = h[HName(l, c, false)]->Integral(0, R::kNM + 1);
      TString line = Form("[RAW]   %-5s %-6s pass %12.1f fail %10.1f -> %.4f", kLinkName[l], kCellName[c], p, fl,
                          p / std::max(1e-9, p + fl));
      if (l == kIdIso)
        line += Form("   [SS pass %.1f fail %.1f]", h[HName(l, c, true, true)]->Integral(0, R::kNM + 1),
                     h[HName(l, c, false, true)]->Integral(0, R::kNM + 1));
      std::cout << line << "\n";
    }
  fo->Write();
  fo->Close();
  std::cout << "[OUT] " << outName << "\n";
  return 0;
}

// ============================================================================
// FIT: eps per bin (fit + cut-and-count), SF, and the e-/e+ ratios
// ============================================================================
namespace pOTnpQ
{
// an asymmetric value: lo / hi = the distances below / above v (>= 0)
struct Asym
{
  double v = 0, lo = 0, hi = 0;
};

struct Res
{
  pORecoSF::FitOut fit;
  double nP = 0, nF = 0, nPss = 0, nFss = 0;   // data, 60-120 (OS; SS for idiso)
  double epsMC = 0, epsMCe = 0;                // DY MC counted, OS 60-120 (the fit's eps_MC)
  double epsMCcc = 0, epsMCcce = 0;            // DY MC, the cut-and-count recipe
  Asym   epsCC, epsRaw;                        // data, profile-likelihood intervals (raw: OS 60-120, no subtraction)
  double sf = 0, sfLo = 0, sfHi = 0;           // the fit (MINOS)
  Asym   sfCC, sfRaw;
};

// ---- a data efficiency as a PROFILE-LIKELIHOOD interval (2026-10-20, user
//      request: eps <= 1, so its uncertainty is asymmetric; until then the
//      cut-and-count and raw errors were linear (Wald) propagation,
//      symmetric, ~2x too small on the low side with 0-12 fails per cell and
//      unphysical above eps = 1). Counts are the DATA's (unweighted):
//        nP ~ Pois(eps S + bP),  nF ~ Pois((1 - eps) S + bF),  0 <= eps <= 1, S >= 0,
//      the backgrounds either KNOWN (reco cut-and-count: bP = 0, bF = the fitted
//      fail background in 80-100; raw: none) or MEASURED by the same-sign pairs
//      (idiso: mP ~ Pois(bP), mF ~ Pois(bF), the OS-minus-SS assumption of equal
//      SS and OS background). S and the measured backgrounds are profiled (the
//      log-likelihood is concave in them: exact coordinate ascent); the 68%
//      interval is {eps : -2 ln[L(eps) / L(eps_hat)] <= 1}, one-sided [eps_lo, 1]
//      when eps_hat = 1. With no background it is the binomial likelihood
//      interval of nP / (nP + nF). eps_hat = the old point estimate (the
//      background-subtracted ratio, capped at 1), checked by the caller.
struct EffModel
{
  double nP = 0, nF = 0;   // signal-region counts
  double bP = 0, bF = 0;   // known backgrounds (measured = false)
  bool   measured = false; // backgrounds from the auxiliary counts below
  double mP = 0, mF = 0;   // same-sign counts (measured = true)
};

inline double PoisLn(double n, double mu) // n ln mu - mu, 0 ln 0 = 0
{
  if (n <= 0.0) return -mu;
  return (mu > 0.0) ? n * std::log(mu) - mu : -1e300;
}

// ln L at fixed eps, maximized over S (and bP, bF when measured)
double ProfLnL(const EffModel &m, double eps)
{
  const double a = eps, c = 1.0 - eps;
  double bP = m.measured ? std::max(m.mP, 0.5) : m.bP, bF = m.measured ? std::max(m.mF, 0.5) : m.bF;
  double S = std::max(m.nP + m.nF - bP - bF, 1.0), prev = -1e301;
  auto lnL = [&]()
  {
    double v = PoisLn(m.nP, a * S + bP) + PoisLn(m.nF, c * S + bF);
    if (m.measured) v += PoisLn(m.mP, bP) + PoisLn(m.mF, bF);
    return v;
  };
  for (int it = 0; it < 20000; ++it)
  {
    // S: dlnL/dS = nP a/(aS + bP) + nF c/(cS + bF) - 1 = 0, decreasing in S
    auto g = [&](double s)
    {
      double v = -1.0;
      if (a > 0.0 && m.nP > 0.0) v += m.nP * a / std::max(a * s + bP, 1e-300);
      if (c > 0.0 && m.nF > 0.0) v += m.nF * c / std::max(c * s + bF, 1e-300);
      return v;
    };
    if (g(0.0) <= 0.0)
      S = 0.0;
    else
    {
      double lo = 0.0, hi = m.nP + m.nF + 10.0; // g(hi) < 0
      for (int k = 0; k < 200 && hi - lo > 1e-13 * hi; ++k) (g(0.5 * (lo + hi)) > 0.0 ? lo : hi) = 0.5 * (lo + hi);
      S = 0.5 * (lo + hi);
    }
    // measured backgrounds: d/db [n ln(s + b) - (s + b) + mm ln b - b] = 0
    //   -> 2 b^2 + (2 s - n - mm) b - mm s = 0
    if (m.measured)
    {
      auto bq = [](double n, double mm, double s)
      {
        const double B = 2.0 * s - n - mm;
        return std::max(0.0, (-B + std::sqrt(B * B + 8.0 * mm * s)) / 4.0);
      };
      bP = bq(m.nP, m.mP, a * S);
      bF = bq(m.nF, m.mF, c * S);
    }
    const double v = lnL();
    if (!m.measured || std::fabs(v - prev) < 1e-13) return v;
    prev = v;
  }
  return lnL();
}

// the MLE (golden section: the profile of a concave likelihood is unimodal)
// and the 68% interval (bisection on -2 Delta lnL = 1)
Asym EffInterval(const EffModel &m)
{
  const double gr = 0.5 * (std::sqrt(5.0) - 1.0);
  double x1 = 0.0, x2 = 1.0;
  for (int k = 0; k < 200 && x2 - x1 > 1e-12; ++k)
  {
    const double l = x2 - gr * (x2 - x1), r = x1 + gr * (x2 - x1);
    if (ProfLnL(m, l) < ProfLnL(m, r)) x1 = l; // the maximum lies in [l, x2]
    else                               x2 = r;
  }
  double e = 0.5 * (x1 + x2);
  if (ProfLnL(m, 1.0) >= ProfLnL(m, e)) e = 1.0; // the boundary MLE (no fail beyond the background)
  if (ProfLnL(m, 0.0) >= ProfLnL(m, e)) e = 0.0;
  const double lmax = ProfLnL(m, e);
  auto q = [&](double x) { return 2.0 * (lmax - ProfLnL(m, x)); };
  Asym r;
  r.v = e;
  if (e > 0.0)
  {
    if (q(0.0) <= 1.0) r.lo = e;
    else
    {
      double lo = 0.0, hi = e;
      for (int k = 0; k < 100; ++k) (q(0.5 * (lo + hi)) > 1.0 ? lo : hi) = 0.5 * (lo + hi);
      r.lo = e - 0.5 * (lo + hi);
    }
  }
  if (e < 1.0)
  {
    if (q(1.0) <= 1.0) r.hi = 1.0 - e;
    else
    {
      double lo = e, hi = 1.0;
      for (int k = 0; k < 100; ++k) (q(0.5 * (lo + hi)) > 1.0 ? hi : lo) = 0.5 * (lo + hi);
      r.hi = 0.5 * (lo + hi) - e;
    }
  }
  return r;
}

// x / y with y Gaussian: relative errors added in quadrature per side
inline Asym DivGauss(const Asym &x, double y, double ey)
{
  Asym r;
  r.v  = x.v / y;
  r.lo = r.v * std::hypot(x.lo / x.v, ey / y);
  r.hi = r.v * std::hypot(x.hi / x.v, ey / y);
  return r;
}
// x / y, both asymmetric: the ratio moves up with x's upper and y's lower side
inline Asym DivAsym(const Asym &x, const Asym &y)
{
  Asym r;
  r.v  = x.v / y.v;
  r.lo = r.v * std::hypot(x.lo / x.v, y.hi / y.v);
  r.hi = r.v * std::hypot(x.hi / x.v, y.lo / y.v);
  return r;
}
// x * y, both asymmetric (independent): per side in quadrature
inline Asym MulAsym(const Asym &x, const Asym &y)
{
  Asym r;
  r.v  = x.v * y.v;
  r.lo = r.v * std::hypot(x.lo / x.v, y.lo / y.v);
  r.hi = r.v * std::hypot(x.hi / x.v, y.hi / y.v);
  return r;
}
// the distance from 1 in units of the error on 1's side
inline double SigmaFromOne(const Asym &x)
{
  const double e = (x.v >= 1.0) ? x.lo : x.hi;
  return (e > 0.0) ? (x.v - 1.0) / e : 0.0;
}
// "v^{+hi}_{-lo}" for the plots' info boxes (TLatex)
inline std::string AsymTxt(const Asym &x)
{
  return Form("%.3f^{+%.3f}_{#minus%.3f}", x.v, x.hi, x.lo);
}

// the sum of a bin's cells (a clone: never mutate file-owned histograms)
TH1D *SumCells(TFile *f, int link, int bin, bool pass, bool ss, const char *tag)
{
  TH1D *o = nullptr;
  for (int i = 0; i < kBins[bin].n; ++i)
  {
    TH1 *x = (TH1 *)f->Get(HName(link, kBins[bin].cells[i], pass, ss));
    if (!x)
    {
      std::cerr << "[ERROR] missing " << HName(link, kBins[bin].cells[i], pass, ss) << " in " << f->GetName() << "\n";
      delete o;
      return nullptr;
    }
    if (!o)
    {
      o = (TH1D *)x->Clone(Form("%s_%s_%s_%s%s", tag, kLinkName[link], kBins[bin].name, pass ? "pass" : "fail", ss ? "_ss" : ""));
      o->SetDirectory(nullptr);
    }
    else
      o->Add(x);
  }
  return o;
}

// integral over [lo, hi) of a 1-GeV-binned histogram (lo, hi on bin edges)
inline double IntRange(const TH1 *h, double lo, double hi, double &err)
{
  const int b1 = h->GetXaxis()->FindBin(lo + 1e-6), b2 = h->GetXaxis()->FindBin(hi - 1e-6);
  return h->IntegralAndError(b1, b2, err);
}

// eps = P / (P + F) with independent errors on P and F
inline void EpsProp(double P, double eP, double F, double eF, double &eps, double &err)
{
  const double T = P + F;
  eps = (T > 0) ? P / T : 0.0;
  err = (T > 0) ? std::sqrt(std::pow(F * eP, 2) + std::pow(P * eF, 2)) / (T * T) : 0.0;
}

// ---- the pp EGM 2025 SFs in their OWN cells (user request 2026-10-20: the
//      original binning, not the average over our probes)
const double kEgmMinFrac = 0.05; // a pT cell is drawn if it holds >= 5% of the passing DY-MC probes
const double kEgmEtaMax  = 2.5;  // the tables' outer edges are +-inf: drawn at the tracker edge
struct EgmCell
{
  const pOSF::Table2D *t;
  int                  ipt;   // pT cell of t
  std::string          label; // legend
  int                  color;
};

// one working point from skim/sf/ (pOSF::LoadEgmCsv checks the grid is complete)
bool LoadEgm(const char *wp, pOSF::Table2D &t)
{
  std::string err, prov;
  for (int i = 0; pOSF::kEleCsvDirs[i]; ++i)
  {
    const std::string path = std::string(pOSF::kEleCsvDirs[i]) + "electron_sf_" + pOSF::kEgmYear + "_" + wp + ".csv";
    if (gSystem->AccessPathName(path.c_str(), kReadPermission)) continue;
    if (!pOSF::LoadEgmCsv(path, wp, t, err, prov)) { std::cerr << "[ERROR] " << err << "\n"; return false; }
    std::cout << "[EGM] " << wp << " <- " << path << " (" << prov << ")\n";
    return true;
  }
  std::cerr << "[ERROR] electron_sf_" << pOSF::kEgmYear << "_" << wp << ".csv not found -- run skim/sf/extract_electron_sf.py\n";
  return false;
}

// the pT cells of the tables (each only within [ptMin, ptMax)) that hold >=
// kEgmMinFrac of the passing probes of a link (DY MC, all six cells summed);
// every candidate's data and MC fraction is printed
std::vector<EgmCell> PickEgmCells(TFile *fd, TFile *fm, int link,
                                  const std::vector<std::tuple<const pOSF::Table2D *, double, double>> &tables)
{
  std::unique_ptr<TH1D> pd, pm;
  for (int c = 0; c < kNCell; ++c)
    for (int s = 0; s < 2; ++s)
    {
      TH1 *x = (TH1 *)(s ? fm : fd)->Get(PtName(link, c));
      if (!x) { std::cerr << "[ERROR] missing " << PtName(link, c) << " -- re-run the skim\n"; return {}; }
      std::unique_ptr<TH1D> &dst = s ? pm : pd;
      if (!dst) { dst.reset((TH1D *)x->Clone(Form("egmpt_%s_%d", kLinkName[link], s))); dst->SetDirectory(nullptr); }
      else dst->Add(x);
    }
  const int col[4] = {kGreen + 2, kOrange + 1, kViolet + 1, kGray + 2};
  std::vector<EgmCell> out;
  auto frac = [](const TH1D *h, double lo, double hi)
  {
    const int b1 = h->GetXaxis()->FindBin(lo + 1e-6), b2 = std::isinf(hi) ? h->GetNbinsX() + 1 : h->GetXaxis()->FindBin(hi - 1e-6);
    return h->Integral(b1, b2) / std::max(1e-12, h->Integral(0, h->GetNbinsX() + 1));
  };
  for (const auto &tb : tables)
  {
    const pOSF::Table2D *t = std::get<0>(tb);
    for (int ip = 0; ip < t->NPt(); ++ip)
    {
      const double lo = t->ptEdges[ip], hi = t->ptEdges[ip + 1];
      if (lo < std::get<1>(tb) || hi > std::get<2>(tb)) continue;
      const double fM = frac(pm.get(), lo, hi), fD = frac(pd.get(), lo, hi);
      const bool   draw = fM >= kEgmMinFrac && out.size() < 4;
      std::cout << Form("[EGM] %-7s %-12s pT [%5.0f, %5.0f): passing probes DY MC %5.1f%%, data %5.1f%% -> %s\n", kLinkName[link],
                        t->name.c_str(), lo, hi, 100.0 * fM, 100.0 * fD, draw ? "drawn" : "not drawn");
      if (!draw) continue;
      const std::string lab = std::isinf(hi) ? Form("p_{T} > %.0f GeV", lo) : Form("p_{T} %.0f#minus%.0f GeV", lo, hi);
      out.push_back({t, ip, lab, col[out.size()]});
    }
  }
  return out;
}

// one translucent box (sfdown..sfup = the EGM +-1 sigma) and the central line
// per eta cell per chosen pT cell; signed eta_SC, or folded onto |eta_SC| with
// both signs drawn (Reco tables are eta-symmetric, wp90iso is not, quite). The
// crack cells are skipped (no probes there; wp90iso holds EGM's placeholder 1)
void DrawEgmCells(const std::vector<EgmCell> &cells, bool folded)
{
  for (const EgmCell &c : cells)
  {
    const pOSF::Table2D &t = *c.t;
    for (int ie = 0; ie < t.NEta(); ++ie)
    {
      double lo = std::max(-kEgmEtaMax, t.etaEdges[ie]), hi = std::min(kEgmEtaMax, t.etaEdges[ie + 1]);
      const double amid = std::fabs(0.5 * (lo + hi));
      if (amid > pOSkim::kEcalGapLo && amid < pOSkim::kEcalGapHi) continue;
      if (folded && hi <= 0.0) { const double l = -hi; hi = -lo; lo = l; }
      const int k = ie * t.NPt() + c.ipt;
      TBox *b = new TBox(lo, t.dn[k], hi, t.up[k]);
      b->SetFillColorAlpha(c.color, 0.30);
      b->SetFillStyle(1001);
      b->Draw();
      TLine *ln = new TLine(lo, t.nom[k], hi, t.nom[k]);
      ln->SetLineColor(c.color);
      ln->SetLineWidth(2);
      ln->Draw();
    }
  }
}

// their legend: a header line + one entry per pT cell
void DrawEgmLegend(const std::vector<EgmCell> &cells, const char *header, double x1, double y1, double x2, double y2)
{
  if (cells.empty()) return;
  TLegend *leg = new TLegend(x1, y1, x2, y2);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.030);
  leg->SetHeader(header);
  for (const EgmCell &c : cells)
  {
    TBox *key = new TBox(0, 0, 1, 1); // legend glyph only
    key->SetFillColorAlpha(c.color, 0.30);
    key->SetFillStyle(1001);
    key->SetLineColor(c.color);
    key->SetLineWidth(2);
    leg->AddEntry(key, c.label.c_str(), "lf");
  }
  leg->Draw();
}
} // namespace pOTnpQ

int TnpQFit()
{
  using namespace pOTnpQ;
  namespace R = pORecoSF;
  namespace I = pOIdIso;
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  TFile *fd = TFile::Open(Form("%s/skim_Data.root", kOutDir));
  TFile *fm = TFile::Open(Form("%s/skim_DY.root", kOutDir));
  if (!fd || fd->IsZombie() || !fm || fm->IsZombie())
  {
    std::cerr << "[FATAL] run the skim first: " << kOutDir << "/skim_{Data,DY}.root\n";
    return 2;
  }
  gSystem->mkdir(kPlotDir, kTRUE);

  // ---- [CHECK] the four cells add up to the parent studies' inclusive numbers
  int nCheckFail = 0;
  for (const char *S : {"Data", "DY"})
  {
    TFile *mine = (std::string(S) == "Data") ? fd : fm;
    // reco vs the current reco_sf_ele.C skim (SC energy); recopho vs its
    // pre-2026-10-20 skim (the photon-object four-vector)
    for (int l : {(int)kReco, (int)kRecoPho})
    {
      const char   *refDir = (l == kReco) ? R::kOutDir : kPhoRefDir;
      const TString refFile = Form("%s/skim_%s.root", refDir, S);
      TFile *fr = gSystem->AccessPathName(refFile) ? nullptr : TFile::Open(refFile);
      if (!fr || fr->IsZombie())
      {
        std::cout << "[CHECK] " << kLinkName[l] << " " << S << ": " << refFile << " absent, not checked\n";
        continue;
      }
      for (int p = 1; p >= 0; --p)
      {
        TH1 *ref = (TH1 *)fr->Get(R::RecoHName(R::kNoId, 0, 0, p == 1));
        std::unique_ptr<TH1D> sum(SumCells(mine, l, 0, p == 1, false, "chk"));
        if (!ref || !sum) { std::cout << "[CHECK] " << kLinkName[l] << " " << S << ": reference histogram missing\n"; ++nCheckFail; continue; }
        double maxRel = 0.0;
        for (int b = 0; b <= ref->GetNbinsX() + 1; ++b)
          maxRel = std::max(maxRel, std::fabs(sum->GetBinContent(b) - ref->GetBinContent(b)) /
                                        std::max(1e-9, std::fabs(ref->GetBinContent(b))));
        const bool ok = maxRel < 1e-9;
        if (!ok) ++nCheckFail;
        std::cout << Form("[CHECK] %-7s %-4s %s: sum of cells %12.1f vs %s h_m_noid_incl0 %12.1f, max rel. bin diff %.1e -> %s\n",
                          kLinkName[l], S, p ? "pass" : "fail", sum->Integral(0, sum->GetNbinsX() + 1),
                          l == kReco ? "reco_sf_ele" : "pre-SC reco_sf_ele", ref->Integral(0, ref->GetNbinsX() + 1), maxRel,
                          ok ? "IDENTICAL" : "DIFFERENT");
      }
      fr->Close();
    }
    TFile *fi = TFile::Open(Form("%s/skim_%s.root", I::kRootDir, S));
    TH1 *ref = (fi && !fi->IsZombie()) ? (TH1 *)fi->Get(I::kTnpProbeHist) : nullptr;
    if (ref)
    {
      const char *lab[4] = {"OS pass", "OS fail", "SS pass", "SS fail"};
      for (int i = 0; i < 4; ++i)
      {
        std::unique_ptr<TH1D> sum(SumCells(mine, kIdIso, 0, i % 2 == 0, i >= 2, "chk"));
        const double a = sum ? sum->Integral(0, sum->GetNbinsX() + 1) : -1.0, b = ref->GetBinContent(i + 1);
        const bool ok = std::fabs(a - b) <= 1e-9 * std::max(1.0, std::fabs(b));
        if (!ok) ++nCheckFail;
        std::cout << Form("[CHECK] idiso   %-4s %s: sum of cells %12.1f vs idiso_sf_skim h_tnp_probes %12.1f -> %s\n", S, lab[i],
                          a, b, ok ? "IDENTICAL" : "DIFFERENT");
      }
    }
    else
      std::cout << "[CHECK] idiso " << S << ": " << I::kRootDir << "/skim_" << S << ".root (" << I::kTnpProbeHist
                << ") absent, not checked\n";
    if (fi) fi->Close();
  }
  if (nCheckFail) std::cout << "[WARN] " << nCheckFail << " consistency check(s) failed -- the split does not reproduce its parent study\n";

  std::ofstream csv(Form("%s/tnp_charge_ele.csv", kPlotDir));
  csv << "link,bin,label,n_pass_data,n_fail_data,n_pass_data_ss,n_fail_data_ss,eps_data,eps_data_errlo,eps_data_errhi,"
         "fit_status,cov_qual,eps_data_cc,eps_data_cc_errlo,eps_data_cc_errhi,eps_mc,eps_mc_err,eps_mc_cc,eps_mc_cc_err,"
         "sf,sf_errlo,sf_errhi,sf_cc,sf_cc_errlo,sf_cc_errhi,eps_data_raw,eps_data_raw_errlo,eps_data_raw_errhi,"
         "sf_raw,sf_raw_errlo,sf_raw_errhi\n";
  // the high-|eta| e-/e+ ratios of each link: [link][side: |eta| > 1.6, fwd, bwd][fit, cut-and-count, raw]
  Asym rSide[kNLink][3][3];
  // every link's results + its cut-and-count e-/e+ ratios (the reco plots
  // show both four-vectors, so they are drawn after the loop)
  const int kNPair = 5;
  Res  resAll[kNLink][kNBin];
  Asym rCCAll[kNLink][kNPair];

  // the pp EGM tables, drawn in their own cells (Reco: Reco20to75 below 75 GeV,
  // RecoAbove75 above, as skim/electron_sf.h applies them)
  pOSF::Table2D egmReco20, egmReco75, egmIdIso;
  if (!LoadEgm(pOSF::kEleReco20WP, egmReco20) || !LoadEgm(pOSF::kEleReco75WP, egmReco75) ||
      !LoadEgm(pOSF::kEleIdIsoWP, egmIdIso))
    return 2;
  const double kInf = std::numeric_limits<double>::infinity();
  const std::vector<EgmCell> egmRecoCells =
      PickEgmCells(fd, fm, kReco, {{&egmReco20, 0.0, pOSF::kEleRecoSplitPt}, {&egmReco75, pOSF::kEleRecoSplitPt, kInf}});
  const std::vector<EgmCell> egmIdIsoCells = PickEgmCells(fd, fm, kIdIso, {{&egmIdIso, 0.0, kInf}});
  if (egmRecoCells.empty() || egmIdIsoCells.empty()) return 2;

  // the SF plots' tuner: the y range, the SF = 1 line, and the EGM cells (under
  // the data points, which SaveNiceGraph_Overlay draws after the tuner) with
  // their legend at egmLeg = {x1, y1, x2, y2} (NDC)
  auto sfTuner = [](double yMin, double yMax, const std::vector<EgmCell> &egm, bool folded, std::string header,
                    std::array<double, 4> egmLeg) -> GraphTuner
  {
    return [=](TCanvas *c, TGraphErrors *g)
    {
      g->SetMinimum(yMin);
      g->SetMaximum(yMax);
      c->cd();
      TLine *one = new TLine(g->GetXaxis()->GetXmin(), 1.0, g->GetXaxis()->GetXmax(), 1.0);
      one->SetLineStyle(2);
      one->SetLineColor(kGray + 2);
      one->Draw();
      DrawEgmCells(egm, folded);
      DrawEgmLegend(egm, header.c_str(), egmLeg[0], egmLeg[1], egmLeg[2], egmLeg[3]);
    };
  };

  for (int l = 0; l < kNLink; ++l)
  {
    std::cout << "\n==================== " << kLinkName[l] << ": " << kLinkLabel[l] << " ====================\n";
    Res *res = resAll[l];
    double mu = 0.0, sg = 0.5;
    for (int b = 0; b < kNBin; ++b)
    {
      Res &r = res[b];
      std::unique_ptr<TH1D> dP(SumCells(fd, l, b, true, false, "d")), dF(SumCells(fd, l, b, false, false, "d"));
      std::unique_ptr<TH1D> tP(SumCells(fm, l, b, true, false, "t")), tF(SumCells(fm, l, b, false, false, "t"));
      if (!dP || !dF || !tP || !tF) return 2;
      std::unique_ptr<TH1D> dPs, dFs, tPs, tFs;
      if (l == kIdIso)
      {
        dPs.reset(SumCells(fd, l, b, true, true, "d"));  dFs.reset(SumCells(fd, l, b, false, true, "d"));
        tPs.reset(SumCells(fm, l, b, true, true, "t"));  tFs.reset(SumCells(fm, l, b, false, true, "t"));
        if (!dPs || !dFs || !tPs || !tFs) return 2;
      }
      double e1 = 0, e2 = 0;
      r.nP = dP->Integral();
      r.nF = dF->Integral();
      if (l == kIdIso) { r.nPss = dPs->Integral(); r.nFss = dFs->Integral(); }

      // eps_MC as in reco_sf_ele.C: counted from the templates over 60-120
      // (10^8 weighted DY pairs: its error stays linear, it is ~1e-4)
      const double mP = tP->IntegralAndError(1, tP->GetNbinsX(), e1), mF = tF->IntegralAndError(1, tF->GetNbinsX(), e2);
      EpsProp(mP, e1, mF, e2, r.epsMC, r.epsMCe);
      // the profile-likelihood interval, and a check that its MLE is the
      // background-subtracted ratio (capped at 1) quoted until now
      auto interval = [&](const EffModel &em, double sP, double sF, const char *what)
      {
        const Asym   a   = EffInterval(em);
        const double old = sP / (sP + std::max(0.0, sF));
        if (std::fabs(a.v - old) > 1e-6)
          std::cout << Form("  [WARN] %s %s %s: likelihood MLE %.6f differs from the subtracted ratio %.6f\n", kLinkName[l],
                            kBins[b].name, what, a.v, old);
        return a;
      };
      // no background subtraction: the lowest efficiency the data allow
      {
        EffModel em;
        em.nP = dP->Integral(1, dP->GetNbinsX());
        em.nF = dF->Integral(1, dF->GetNbinsX());
        r.epsRaw = interval(em, em.nP, em.nF, "raw");
      }

      // the fit: the Gaussian floats in the inclusive fit, fixed to it in every other bin
      const bool isIncl = (b == 0);
      r.fit = R::RecoFitEps(dP.get(), dF.get(), tP.get(), tF.get(), R::kMLo, R::kMHi, isIncl, isIncl ? 0.0 : mu,
                            isIncl ? 0.5 : sg, 0, Form("%s/fit_%s_%s", kPlotDir, kLinkName[l], kBins[b].name),
                            kBins[b].label, kLinkHeader[l]);
      if (isIncl) { mu = r.fit.mu; sg = r.fit.sg; }
      if (r.fit.status != 0)
        std::cout << Form("  [WARN] %s %s: the fit did not converge (status %d) -- use the cut-and-count\n", kLinkName[l],
                          kBins[b].name, r.fit.status);

      // the cut-and-count (a cell with no fail left in the window now simply
      // gets a one-sided interval [eps_lo, 1]: the 1.841 zero-count patch of
      // the linear propagation is gone)
      if (IsReco(l))
      { // 80-100 GeV, the fitted fail background in the window as a known background (reco_sf_ele.C's [CC])
        double eP = 0, eF = 0, mePw = 0, meFw = 0;
        EffModel em;
        em.nP = IntRange(dP.get(), 80.0, 100.0, eP);
        em.nF = IntRange(dF.get(), 80.0, 100.0, eF);
        em.bF = r.fit.bFwin;
        r.epsCC = interval(em, em.nP, em.nF - em.bF, "cc");
        const double mpw = IntRange(tP.get(), 80.0, 100.0, mePw), mfw = IntRange(tF.get(), 80.0, 100.0, meFw);
        EpsProp(mpw, mePw, mfw, meFw, r.epsMCcc, r.epsMCcce);
      }
      else
      { // OS minus SS over 60-120: the SS counts measure the OS background (data); MC subtracted
        auto osMinusSs = [](const TH1 *os, const TH1 *ss, double &err)
        {
          double eo = 0, es = 0;
          const double v = os->IntegralAndError(1, os->GetNbinsX(), eo) - ss->IntegralAndError(1, ss->GetNbinsX(), es);
          err = std::hypot(eo, es);
          return v;
        };
        EffModel em;
        em.nP = dP->Integral(1, dP->GetNbinsX());
        em.nF = dF->Integral(1, dF->GetNbinsX());
        em.measured = true;
        em.mP = dPs->Integral(1, dPs->GetNbinsX());
        em.mF = dFs->Integral(1, dFs->GetNbinsX());
        r.epsCC = interval(em, em.nP - em.mP, em.nF - em.mF, "cc");
        double meP = 0, meF = 0;
        const double MP = osMinusSs(tP.get(), tPs.get(), meP), MF = osMinusSs(tF.get(), tFs.get(), meF);
        EpsProp(MP, meP, MF, meF, r.epsMCcc, r.epsMCcce);
      }

      r.sf    = r.fit.eps / r.epsMC;
      r.sfLo  = r.sf * std::hypot(r.fit.lo / r.fit.eps, r.epsMCe / r.epsMC);
      r.sfHi  = r.sf * std::hypot(r.fit.hi / r.fit.eps, r.epsMCe / r.epsMC);
      r.sfCC  = DivGauss(r.epsCC, r.epsMCcc, r.epsMCcce);
      r.sfRaw = DivGauss(r.epsRaw, r.epsMC, r.epsMCe);

      TString ssPart = (l == kIdIso) ? TString(Form(" (SS %3.0f/%2.0f)", r.nPss, r.nFss)) : TString("");
      std::cout << Form("[RESULT] %-5s %-6s %-26s data pass %4.0f fail %3.0f%s | fit eps %.4f +%.4f -%.4f (status %d)"
                        " | cc %.4f +%.4f -%.4f | raw %.4f +%.4f -%.4f | MC %.4f (cc %.4f) | SF %.4f +%.4f -%.4f"
                        " | SF_cc %.4f +%.4f -%.4f | SF_raw %.4f +%.4f -%.4f\n",
                        kLinkName[l], kBins[b].name, kBins[b].label, r.nP, r.nF, ssPart.Data(), r.fit.eps, r.fit.hi, r.fit.lo,
                        r.fit.status, r.epsCC.v, r.epsCC.hi, r.epsCC.lo, r.epsRaw.v, r.epsRaw.hi, r.epsRaw.lo, r.epsMC,
                        r.epsMCcc, r.sf, r.sfHi, r.sfLo, r.sfCC.v, r.sfCC.hi, r.sfCC.lo, r.sfRaw.v, r.sfRaw.hi, r.sfRaw.lo);
      csv << Form("%s,%s,\"%s\",%.1f,%.1f,%.1f,%.1f,%.5f,%.5f,%.5f,%d,%d,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,"
                  "%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n",
                  kLinkName[l], kBins[b].name, kBins[b].label, r.nP, r.nF, r.nPss, r.nFss, r.fit.eps, r.fit.lo, r.fit.hi,
                  r.fit.status, r.fit.covQual, r.epsCC.v, r.epsCC.lo, r.epsCC.hi, r.epsMC, r.epsMCe, r.epsMCcc, r.epsMCcce,
                  r.sf, r.sfLo, r.sfHi, r.sfCC.v, r.sfCC.lo, r.sfCC.hi, r.epsRaw.v, r.epsRaw.lo, r.epsRaw.hi, r.sfRaw.v,
                  r.sfRaw.lo, r.sfRaw.hi);
    }

    // ---- the e-/e+ ratios. The probes are disjoint, so the errors are
    //      independent; asymmetric errors are combined per side (the ratio's
    //      lower side from e-'s lower and e+'s upper side, and vice versa) --
    //      the fit's MINOS errors, the cut-and-count and raw likelihood
    //      intervals alike. "sigma" = the distance from 1 in units of the error
    //      on 1's side.
    std::cout << "  ---- e- / e+ (disjoint probes: independent errors, asymmetric errors combined per side) ----\n";
    const int   pairs[kNPair][2] = {{kBinEm, kBinEp}, {kBinEmLo, kBinEpLo}, {kBinEmHi, kBinEpHi},
                                    {kBinEmFwd, kBinEpFwd}, {kBinEmBwd, kBinEpBwd}};
    const char *pname[kNPair] = {"all |eta_SC|", "|eta_SC| < 1.6", "|eta_SC| > 1.6", "fwd eta_SC<-1.6", "bwd eta_SC>+1.6"};
    Asym rFit[kNPair], rCC[kNPair], rRaw[kNPair];
    bool fitOk[kNPair];
    for (int i = 0; i < kNPair; ++i)
    {
      const Res &m = res[pairs[i][0]], &p = res[pairs[i][1]];
      rFit[i]  = DivAsym({m.sf, m.sfLo, m.sfHi}, {p.sf, p.sfLo, p.sfHi});
      rCC[i]   = DivAsym(m.sfCC, p.sfCC);
      rRaw[i]  = DivAsym(m.sfRaw, p.sfRaw);
      fitOk[i] = (m.fit.status == 0 && p.fit.status == 0);
      rCCAll[l][i] = rCC[i];
      const Asym rD = DivAsym({m.fit.eps, m.fit.lo, m.fit.hi}, {p.fit.eps, p.fit.lo, p.fit.hi});
      const double rM = m.epsMC / p.epsMC, eM = rM * std::hypot(m.epsMCe / m.epsMC, p.epsMCe / p.epsMC);
      std::cout << Form("[RATIO] %-5s %-15s SF(e-)/SF(e+) = fit %.4f +%.4f -%.4f (%+.1f sigma)%s | cc %.4f +%.4f -%.4f (%+.1f sigma)"
                        " | raw %.4f +%.4f -%.4f (%+.1f sigma) | eps_data(e-)/eps_data(e+) %.4f +%.4f -%.4f"
                        " | eps_MC(e-)/eps_MC(e+) %.4f +- %.4f\n",
                        kLinkName[l], pname[i], rFit[i].v, rFit[i].hi, rFit[i].lo, SigmaFromOne(rFit[i]),
                        fitOk[i] ? "" : " [NOT CONVERGED]", rCC[i].v, rCC[i].hi, rCC[i].lo, SigmaFromOne(rCC[i]), rRaw[i].v,
                        rRaw[i].hi, rRaw[i].lo, SigmaFromOne(rRaw[i]), rD.v, rD.hi, rD.lo, rM, eM);
    }
    {
      const Asym dr = DivAsym(rFit[2], rFit[1]);
      std::cout << Form("[RATIO] %-5s high / low |eta_SC|: [SF(e-)/SF(e+)]_{>1.6} / [SF(e-)/SF(e+)]_{<1.6} = %.4f +%.4f -%.4f (fit)\n",
                        kLinkName[l], dr.v, dr.hi, dr.lo);
    }
    for (int s = 0; s < 3; ++s)
    {
      rSide[l][s][0] = rFit[2 + s];
      rSide[l][s][1] = rCC[2 + s];
      rSide[l][s][2] = rRaw[2 + s];
    }

    // ---- the same-sign fraction of the passing idiso probes, data vs MC (charge flips + background)
    if (l == kIdIso)
      for (int c = 0; c < kNCell; ++c)
      {
        TH1 *dO = (TH1 *)fd->Get(HName(l, c, true)), *dS = (TH1 *)fd->Get(HName(l, c, true, true));
        TH1 *mO = (TH1 *)fm->Get(HName(l, c, true)), *mS = (TH1 *)fm->Get(HName(l, c, true, true));
        const double a = dS->Integral(), b = dO->Integral(), x = mS->Integral(), y = mO->Integral();
        std::cout << Form("[SSFRAC] idiso %-27s passing probes SS/(OS+SS): data %3.0f/%4.0f = %.4f | DY MC %.4f\n",
                          kCellLabel[c], a, a + b, a / std::max(1e-9, a + b), x / std::max(1e-9, x + y));
      }

    // ---- the ID+iso plots: SF vs |eta_SC| and vs signed eta_SC, e- vs e+ (the
    //      repo's standard overlay): the fit (filled; a bin whose fit did not
    //      converge is left out and named in the box) and the cut-and-count
    //      (open), statistical errors, over the EGM wp90iso cells. The box
    //      quotes the fit ratio where both fits converged, else the
    //      cut-and-count. The reco plots (both four-vectors) follow the loop.
    if (l != kIdIso) continue;
    auto rLine = [&](int i, const char *what)
    {
      return fitOk[i] ? std::string(Form("%s: %s (fit)", what, AsymTxt(rFit[i]).c_str()))
                      : std::string(Form("%s: %s (c&c)", what, AsymTxt(rCC[i]).c_str()));
    };
    // ps = the layout (header, box, legend); yMin/yMax = the SF range;
    // shortLab = "e^{-} fit" instead of "e^{-} probes, fit" (narrow legends);
    // folded = x is |eta_SC| (the EGM cells of both signs drawn on it);
    // egmLeg = the NDC corners of the EGM legend
    auto drawSF = [&](int npt, const double *x, const double *ex, const int *idxM, const int *idxP,
                      const char *const *ptLab, const char *stem, const char *xTitle, std::vector<std::string> box,
                      PlotStyle ps, double yMin, double yMax, bool shortLab, bool folded, std::array<double, 4> egmLeg)
    {
      const int   *idx[2] = {idxM, idxP};
      const char  *qn[2] = {"em", "ep"}, *qLab[2] = {"e^{-}", "e^{+}"};
      const int    qCol[2] = {kRed + 1, kAzure + 2}, qMk[2] = {33, 21}, qMkOpen[2] = {27, 25};
      const double qSz[2] = {2.0, 1.5}, qShift[2][2] = {{-0.6, -0.2}, {0.2, 0.6}};
      std::vector<OverlaySeries> ser(4);
      std::string skipped;
      int nSkipped = 0;
      for (int q = 0; q < 2; ++q)
      {
        TGraphErrors      *g  = new TGraphErrors(), *gc = new TGraphErrors(npt);
        TGraphAsymmErrors *gA = new TGraphAsymmErrors(), *gcA = new TGraphAsymmErrors(npt);
        g->SetName(Form("g_%s_%s_%s", stem, kLinkName[l], qn[q]));
        gA->SetName(Form("g_%s_asym_%s_%s", stem, kLinkName[l], qn[q]));
        gc->SetName(Form("g_%s_cc_%s_%s", stem, kLinkName[l], qn[q]));
        gcA->SetName(Form("g_%s_cc_asym_%s_%s", stem, kLinkName[l], qn[q]));
        for (int k = 0; k < npt; ++k)
        {
          const Res &r = res[idx[q][k]];
          if (r.fit.status == 0)
          {
            const int ip = g->GetN();
            g->SetPoint(ip, x[k], r.sf);
            g->SetPointError(ip, ex[k], 0.5 * (r.sfLo + r.sfHi));
            gA->SetPoint(ip, x[k], r.sf);
            gA->SetPointError(ip, ex[k], ex[k], r.sfLo, r.sfHi);
          }
          else
          {
            skipped += std::string(skipped.empty() ? "" : ", ") + Form("%s %s", qLab[q], ptLab[k]);
            ++nSkipped;
          }
          gc->SetPoint(k, x[k], r.sfCC.v);
          gc->SetPointError(k, ex[k], 0.5 * (r.sfCC.lo + r.sfCC.hi));
          gcA->SetPoint(k, x[k], r.sfCC.v);
          gcA->SetPointError(k, ex[k], ex[k], r.sfCC.lo, r.sfCC.hi);
        }
        OverlaySeries &sf = ser[2 * q], &sc = ser[2 * q + 1];
        sf.gTot = g;  sf.gTotAsym = gA;
        sf.label = Form(shortLab ? "%s fit" : "%s probes, fit", qLab[q]);
        sf.color = qCol[q]; sf.marker = qMk[q]; sf.markerSize = qSz[q]; sf.xShift = qShift[q][0];
        sc.gTot = gc; sc.gTotAsym = gcA;
        sc.label = Form(shortLab ? "%s cut-and-count" : "%s probes, cut-and-count", qLab[q]);
        sc.color = qCol[q]; sc.marker = qMkOpen[q]; sc.markerSize = qSz[q]; sc.xShift = qShift[q][1];
      }
      if (nSkipped > 2)
        box.push_back(Form("fit not converged in %d of %d points (c&c drawn)", nSkipped, 2 * npt));
      else if (!skipped.empty())
        box.push_back("fit not converged: " + skipped);
      SaveNiceGraph_Overlay(ser, Form("%s/%s_%s", kPlotDir, stem, kLinkName[l]), xTitle,
                            "SF (#varepsilon_{data} / #varepsilon_{MC})", "", "Z#rightarrowee tag-and-probe",
                            Form("%s, by probe charge", kLinkShort[l]), box, ps,
                            sfTuner(yMin, yMax, egmIdIsoCells, folded, "EGM pp 2025 wp90iso (e^{+} + e^{-})", egmLeg));
    };
    PlotStyle psBase;
    psBase.headerX = 0.50;
    psBase.headerY = 0.90;
    psBase.subSize = 0.034;
    psBase.boxTextSize = 0.026;
    {
      const double x[2] = {0.8, 2.0}, ex[2] = {0.8, 0.4};
      const int    iM[2] = {kBinEmLo, kBinEmHi}, iP[2] = {kBinEpLo, kBinEpHi};
      const char  *lab[2] = {"< 1.6", "> 1.6"};
      PlotStyle ps = psBase;
      ps.boxX1 = 0.40; ps.boxX2 = 0.93; ps.boxY1 = 0.56; ps.boxY2 = 0.74;
      ps.legX1 = 0.17; ps.legY1 = 0.15; ps.legX2 = 0.48; ps.legY2 = 0.32;
      drawSF(2, x, ex, iM, iP, lab, "sf_charge", "|#eta_{SC}|",
             {"SF(e^{-}) / SF(e^{+})", rLine(1, "|#eta_{SC}| < 1.6"), rLine(2, "|#eta_{SC}| > 1.6")}, ps, 0.40, 1.70,
             true, true, {0.52, 0.15, 0.93, 0.32});
    }
    { // signed eta_SC: the ID+iso points reach down to ~0.7, so the range opens
      // downward and the legend goes to the upper left, under the CMS label
      const double x[3] = {-2.0, 0.0, 2.0}, ex[3] = {0.4, 1.6, 0.4};
      const int    iM[3] = {kBinEmFwd, kBinEmLo, kBinEmBwd}, iP[3] = {kBinEpFwd, kBinEpLo, kBinEpBwd};
      const char  *lab[3] = {"fwd", "|#eta| < 1.6", "bwd"};
      PlotStyle ps = psBase;
      ps.boxX1 = 0.45; ps.boxX2 = 0.93; ps.boxY1 = 0.53; ps.boxY2 = 0.71; // centred band, grows with the lines
      ps.legX1 = 0.17; ps.legY1 = 0.57; ps.legX2 = 0.44; ps.legY2 = 0.76;
      drawSF(3, x, ex, iM, iP, lab, "sf_charge_signed", "#eta_{SC}   (forward = proton-going = #eta_{SC} < 0)",
             {"SF(e^{-}) / SF(e^{+})", rLine(3, "fwd #eta_{SC} < #minus1.6"), rLine(1, "|#eta_{SC}| < 1.6"),
              rLine(4, "bwd #eta_{SC} > 1.6")},
             ps, 0.10, 2.20, true, false, {0.17, 0.14, 0.55, 0.31});
    }
  }

  // ---- the reco plots: THREE cases (user request 2026-10-20) -- the
  //      cut-and-count SF with the SC-energy probe four-vector (filled) and with
  //      the photon-object one (open), per probe charge, over the EGM Reco cells
  //      (charge-inclusive, in their own eta x pT binning). The cut-and-count is
  //      the one extraction that exists in every cell (the fits rarely converge
  //      with 1-11 fails per cell; they stay in the log and the CSV).
  {
    auto drawReco = [&](int npt, const double *x, const double *ex, const int *idxM, const int *idxP,
                        const int *pairIdx, const char *const *pairLab, const char *stem, const char *xTitle,
                        PlotStyle ps, double yMin, double yMax, bool folded, std::array<double, 4> egmLeg)
    {
      const int   *idx[2] = {idxM, idxP};
      const char  *qn[2] = {"em", "ep"}, *qLab[2] = {"e^{-}", "e^{+}"};
      const int    qCol[2] = {kRed + 1, kAzure + 2}, qMk[2] = {33, 21}, qMkOpen[2] = {27, 25};
      const double qSz[2] = {2.0, 1.5}, qShift[2][2] = {{-0.6, -0.2}, {0.2, 0.6}};
      const int    lk[2] = {kReco, kRecoPho};
      const char  *lkLab[2] = {"SC energy", "photon p4"};
      std::vector<OverlaySeries> ser;
      for (int q = 0; q < 2; ++q)
        for (int v = 0; v < 2; ++v)
        {
          TGraphErrors      *g  = new TGraphErrors(npt);
          TGraphAsymmErrors *gA = new TGraphAsymmErrors(npt);
          g->SetName(Form("g_%s_reco_%s_%s", stem, kLinkName[lk[v]], qn[q]));
          gA->SetName(Form("g_%s_reco_asym_%s_%s", stem, kLinkName[lk[v]], qn[q]));
          for (int k = 0; k < npt; ++k)
          {
            const Asym &sfc = resAll[lk[v]][idx[q][k]].sfCC;
            g->SetPoint(k, x[k], sfc.v);
            g->SetPointError(k, ex[k], 0.5 * (sfc.lo + sfc.hi));
            gA->SetPoint(k, x[k], sfc.v);
            gA->SetPointError(k, ex[k], ex[k], sfc.lo, sfc.hi);
          }
          OverlaySeries s;
          s.gTot = g;
          s.gTotAsym = gA;
          s.label = Form("%s, %s", qLab[q], lkLab[v]);
          s.color = qCol[q];
          s.marker = v == 0 ? qMk[q] : qMkOpen[q];
          s.markerSize = qSz[q];
          s.xShift = qShift[q][v];
          ser.push_back(s);
        }
      std::vector<std::string> box = {"SF(e^{-}) / SF(e^{+}):  SC energy | photon p4"};
      for (int i = 0; i < npt; ++i)
      {
        const int p = pairIdx[i];
        box.push_back(Form("%s:  %s | %s", pairLab[i], AsymTxt(rCCAll[kReco][p]).c_str(),
                           AsymTxt(rCCAll[kRecoPho][p]).c_str()));
      }
      SaveNiceGraph_Overlay(ser, Form("%s/%s_reco", kPlotDir, stem), xTitle, "SF (#varepsilon_{data} / #varepsilon_{MC})", "",
                            "Z#rightarrowee tag-and-probe", "electron reco, cut-and-count", box, ps,
                            sfTuner(yMin, yMax, egmRecoCells, folded, "EGM pp 2025 Reco (e^{+} + e^{-})", egmLeg));
    };
    PlotStyle psBase;
    psBase.headerX = 0.50;
    psBase.headerY = 0.90;
    psBase.subSize = 0.030;
    psBase.boxTextSize = 0.026;
    {
      const double x[2] = {0.8, 2.0}, ex[2] = {0.8, 0.4};
      const int    iM[2] = {kBinEmLo, kBinEmHi}, iP[2] = {kBinEpLo, kBinEpHi}, pr[2] = {1, 2};
      const char  *pl[2] = {"|#eta_{SC}| < 1.6", "|#eta_{SC}| > 1.6"};
      PlotStyle ps = psBase;
      ps.boxX1 = 0.40; ps.boxX2 = 0.93; ps.boxY1 = 0.56; ps.boxY2 = 0.74;
      ps.legX1 = 0.18; ps.legY1 = 0.15; ps.legX2 = 0.52; ps.legY2 = 0.33;
      drawReco(2, x, ex, iM, iP, pr, pl, "sf_charge", "|#eta_{SC}|", ps, 0.70, 1.45, true, {0.54, 0.15, 0.93, 0.31});
    }
    {
      const double x[3] = {-2.0, 0.0, 2.0}, ex[3] = {0.4, 1.6, 0.4};
      const int    iM[3] = {kBinEmFwd, kBinEmLo, kBinEmBwd}, iP[3] = {kBinEpFwd, kBinEpLo, kBinEpBwd}, pr[3] = {3, 1, 4};
      const char  *pl[3] = {"fwd #eta_{SC} < #minus1.6", "|#eta_{SC}| < 1.6", "bwd #eta_{SC} > 1.6"};
      PlotStyle ps = psBase;
      ps.boxX1 = 0.36; ps.boxX2 = 0.93; ps.boxY1 = 0.52; ps.boxY2 = 0.72;
      ps.legX1 = 0.18; ps.legY1 = 0.15; ps.legX2 = 0.52; ps.legY2 = 0.33;
      drawReco(3, x, ex, iM, iP, pr, pl, "sf_charge_signed", "#eta_{SC}   (forward = proton-going = #eta_{SC} < 0)", ps,
               0.70, 1.45, false, {0.54, 0.15, 0.93, 0.31});
    }
  }

  {
    std::cout << "\n";
    const char *how[3]  = {"fit", "cc ", "raw"};
    const char *side[3] = {"|eta_SC| > 1.6", "fwd eta_SC<-1.6", "bwd eta_SC>+1.6"};
    for (int lr : {(int)kReco, (int)kRecoPho})
      for (int s = 0; s < 3; ++s)
        for (int i = 0; i < 3; ++i)
        {
          const Asym v = MulAsym(rSide[lr][s][i], rSide[kIdIso][s][i]);
          std::cout << Form("[RATIO] %s x idiso, %s, %s: SF(e-)/SF(e+) = %.4f +%.4f -%.4f (%+.1f sigma from 1;"
                            " the two links treated as independent)\n",
                            kLinkName[lr], side[s], how[i], v.v, v.hi, v.lo, SigmaFromOne(v));
        }
  }
  csv.close();
  std::cout << "[OUT] " << kPlotDir
            << "/{fit_<link>_<bin>_{pass,fail},sf_charge_{reco,idiso},sf_charge_signed_{reco,idiso}}.{png,pdf}, tnp_charge_ele.csv\n";
  return 0;
}

int tnp_charge_ele(const char *mode = "all", const char *sample = "", Long64_t nmax = -1)
{
  const std::string md = mode, sm = sample;
  if (md == "skim")
  {
    if (!sm.empty()) return TnpQSkim(sm, nmax);
    for (const char *s : {"Data", "DY"})
      if (const int rc = TnpQSkim(s, nmax)) return rc;
    return 0;
  }
  if (md == "fit") return TnpQFit();
  if (md == "all")
  {
    for (const char *s : {"Data", "DY"})
      if (const int rc = TnpQSkim(s, nmax)) return rc;
    return TnpQFit();
  }
  std::cerr << "usage: tnp_charge_ele.C+(\"skim\"[, \"Data\"|\"DY\"] | \"fit\" | \"all\")\n";
  return 1;
}
