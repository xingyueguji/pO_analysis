// skim/electron_sf.h -- electron efficiency SCALE FACTORS (data/MC) for the MC
// event weight of skim_Wel / skim_Zee (2026-09-29, the electron step of the
// SF-application phase; user decision: EGM Reco + ID+Iso SFs + a uniform
// trigger SF from our own MB study; 2026-09-30: the Reco factor replaced by our
// own pO tag-and-probe, one inclusive number; 2026-10-01: back to the EGM Reco
// tables, user decision -- the pO tag-and-probe stays as the cross-check).
//
// SINGLE SOURCE for: which electron SFs are applied, where they come from, how
// they are looked up, how the per-event factor and its +-1 sigma variations are
// formed, and the name of the resulting shape systematic. The flavour-blind
// machinery -- the (eta x pT) table with its lookup counters, EventSF, the
// per-source / combined twins and the <SF> record -- lives in muon_sf.h and is
// shared; skim.C only calls W() / Z() and fills the twins.
//
// WHAT IS APPLIED (MC only; data is always weight 1):
//   RECO  : the efficiency for a supercluster to be reconstructed as an
//           electron, ID-independent. kEleRecoSource selects it.
//           *** SINCE 2026-10-21 NOT APPLIED (kRecoUnit, user decision): factor
//           1 for every electron, with a flat +-10% per electron as its OWN
//           nuisance `eleReco` (one inclusive number, Z->ee legs coherent), so
//           the fit measures the electron Reco efficiency -- the pO
//           tag-and-probe gives ~1 (0.988 photon-p4 / 0.998 SC-p4 probe) where
//           the pp tables give ~0.93. ***
//           2026-10-01..2026-10-20 (user decision: the EGamma official SF), as on
//           2026-09-29: the pp EGM 2025Prompt tables Reco20to75 (pT < 75) /
//           RecoAbove75, per electron in (signed SC eta, pT) -- <RECO> = 0.929
//           per W electron, 0.865 per Z->ee pair; the table is eta-symmetric,
//           0.91-0.95 at 20-75 GeV, errors 0.4-1.1% per cell. There is no
//           RecoBelow20 table: the 10-20 GeV Z legs take the 20-45 GeV cell
//           (counted, cntReco.ptClamped).
//           THE pO CROSS-CHECK (kRecoTnP; the applied factor on 2026-09-30
//           only): ONE INCLUSIVE NUMBER for every electron from our own pO
//           Z->ee tag-and-probe, correction/reco_sf_ele.C (probe = a photon
//           supercluster with no ID, pass = any forest electron on it; nominal
//           probe variant, inclusive fit), read from
//           correction/electron_sf_check/reco_sf/rootfile/reco_sf_ele_result.root:
//           SF = 0.9882 +0.0181 -0.0226 (stat, MINOS) +- 0.0016 (syst), data
//           547 pass / 53 fail pairs, vs 0.929 for the pp tables on the same
//           probes (~2.6 sigma apart); flat within its errors (barrel 0.993,
//           endcap 0.970, E_T^SC 20-45 GeV 0.984). Its fit (2026-09-30) moved
//           sigma_Z e/mu 1.27 -> 1.14 and sigma_W e/mu 1.01 -> 0.96
//           (docs/iterations/). The skim reads that file only under kRecoTnP.
//   IDISO : EGM wp90iso -- ONE SF for the MVA ID trained with isolation inputs
//           at 90% efficiency, relative to a reconstructed electron. Applied to
//           electrons passing THE electron selection, eleMVAIdWP90 &&
//           eleMVAIsoWP90 && relIso < 0.3 (skim_common.h kEleIDBranch /
//           kEleIsoBranch / kEleIsoRelMax). Those are the forest's heavy-ion
//           classifiers, NOT EGM VID; the pp SF is applied because the ID+iso SF
//           cross-check (correction/idiso_sf_*, 2026-09-25/26) measured SF =
//           0.959 +- 0.025 for exactly this selection in pO data (inclusive Z
//           tag-and-probe) against the EGM prediction 0.956 +- 0.009 for the
//           same probes (pull +0.11).
//           IDISO (and the EGM Reco record) are read from
//           skim/sf/electron_sf_2025Prompt_<WP>.csv, flattened from
//           skim/sf/EGM/electron.json.gz by
//           skim/sf/extract_electron_sf.py (source md5 in each CSV header, echoed
//           in every job log). Binned in SIGNED supercluster eta x pT;
//           sfup/sfdown = sf +- sqrt(stat^2 + syst^2) (the script verifies it
//           per cell). The crack cells hold EGM placeholders (sf = 1) and are
//           never looked up: the crack is vetoed (skim_common.h InEcalGap).
//   TRIG  : the analysis path HLT_OxyL1SingleEG10_v1 fired AND the leading
//           electron matched within dR < 0.4 (= skim step 8), measured by
//           correction/trig_eff_mb.C with the minimum-bias denominator, data vs
//           the W signal MC, pT > 25, on the m_T > 40 selection variant (the
//           purer W sample: the plain W selection is ~50% QCD fakes in data,
//           whose lower efficiency drags the plain-selection SF down -- the
//           muon made the same choice). Applied as ONE INCLUSIVE number (user:
//           "a uniform electron trigger SF"), Clopper-Pearson 68% on the data
//           count. The per-|y| table of the same study is read and printed in
//           the [SF] block as the record, marked NOT APPLIED (the muon pattern).
//           Run correction/run_trig_eff_mb.sh ele (and, under kRecoTnP only,
//           correction/run_reco_sf.sh) BEFORE the electron MC skim.
//
// *** ETA CONVENTION -- read before touching a call site. ***
// W() / Z() take the RAW SUPERCLUSTER eta (eleSCEta). The analysis bins in
// y = -eta_lab, but that flip is a BINNING convention and must never reach the
// EGM lookup: the tables are in signed SC eta and wp90iso is genuinely
// ASYMMETRIC (20-35 GeV: 0.9216 at -2.0 < eta < -1.566 vs 0.9362 at +1.566 <
// eta < +2.0; 0.9584 vs 0.9677 beyond |eta| = 2), so passing y would mirror the
// correction with no error message -- the gen_xsec.C 2026-08 bug class, and the
// check the muon memory asked for "if a future SF (electrons!) is added".
// Verified by unit test 2026-09-29 (W(35, eta) returns the table at +eta, never
// at -eta, at eta = +-1.0, +-1.8, +-2.2). The parameters are named scEtaLab so
// the convention is visible at the call site; skim.C passes eleSCEta->at(i).
// (The EGM Reco table is eta-symmetric and the kRecoTnP factor does not depend
// on eta at all; wp90iso is the one that would be mirrored.)
//
// HOW IT ENTERS (skim.C):
//   W -> e nu : one factor per event for the LEADING electron after step 8,
//               SF = RECO(eta_SC, pT) x [passIso ? IDISO(eta_SC, pT) : 1] x TRIG,
//               folded into the weight of EVERY fill. The anti-iso sideband (and
//               the MVA-iso-failing gap) gets RECO x TRIG only: its electrons
//               pass the ID but not the iso, for which the combined wp90iso SF is
//               not defined -- and its EWK content is a few %.
//   Z -> e e  : both legs, SF = prod_legs RECO x [iso ? IDISO : 1], times the
//               per-EVENT trigger factor [1 - (1 - eps_data)^2] / [1 - (1 -
//               eps_MC)^2] with the inclusive per-lepton PATH-FIRED
//               efficiencies of the same study (skim_Zee requires the bit, no
//               matching). NB those are measured at pT > 25 and are applied to
//               the 10-25 GeV legs too, where EG10 is still turning on (MC eps
//               0.67 / 0.85 / 0.93 at 10-16 GeV): the event factor is 1 - O(1e-3)
//               either way, carried by the harder leg.
//   NOT covered : the DY-veto second electron (no SF, as for muons).
//
// UNCERTAINTIES: three independent sources (Reco, ID+iso, trigger), each
// varied +-1 sigma with the other two nominal (EventSF::var), stored per
// fit-template histogram as the per-SOURCE twins <h>_eleReco/eleIDIso/
// eleTrigUp|Down (diagnostics) and combined bin by bin in quadrature into THE
// ONE nuisance carried into the fit, <h>_eleSFUp/Down -- the muon muSF recipe
// (muon_sf.h FinalizeSFTwins). The EGM Reco and ID+iso errors are applied
// COHERENTLY (every electron at +1 sigma of its own cell at once): they are
// syst-dominated at 25-100 GeV (Reco ~+-0.7%, ID+iso ~+-1.0% on the W
// templates). Under kRecoTnP the Reco factor is one measurement with one
// (asymmetric) error, coherent by definition: +1.8% / -2.3% per electron,
// twice that per Z->ee pair -- statistics-dominated (53 failing probes).
// Z legs move together likewise.
// Under kRecoUnit (since 2026-10-21) the Reco source is +-10% per electron
// (W: x1.10 / x0.90; Z->ee: x1.21 / x0.81) and is NOT added into eleSF: its
// per-source twins <h>_eleRecoUp/Down are shipped as a SECOND nuisance,
// `eleReco` (kEleSFSystNames), so its pull is the fitted Reco efficiency.
//
// Because the SF is part of the nominal weight, the LHE member weights are
// rescaled by it (pOLhe::ScaleMemberWeights), as for the muon.
//
// Env switch: PO_ELE_SF=off fills MC without SFs (regression checks); the job
// log always says which. Missing inputs are FATAL in the skim -- never a
// silent SF = 1.
#ifndef PO_ELECTRON_SF_H
#define PO_ELECTRON_SF_H

#include "muon_sf.h" // shared: Table2D, LookupCounters, FindEdgeBin, EventSF, SFTwins, SFStats

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace pOSF
{

// ============================================================================
// Configuration -- edit here, nowhere else
// ============================================================================
// CSV directories tried in order (cwd = skim/ under run_all.sh, the repo root
// or a sibling directory otherwise); file = electron_sf_<kEgmYear>_<WP>.csv
inline const char *kEleCsvDirs[] = {"sf/", "skim/sf/", "../skim/sf/", nullptr};
inline const char *kEgmYear      = "2025Prompt";
inline const char *kEleIdIsoWP   = "wp90iso";
inline const char *kEleReco20WP  = "Reco20to75";
inline const char *kEleReco75WP  = "RecoAbove75";
inline constexpr double kEleRecoSplitPt = 75.0; // Reco20to75 below, RecoAbove75 at/above
// WHICH RECO factor is applied (switching = this constant + an electron MC
// re-skim + the downstream re-run):
//   kRecoEgm  the pp EGM 2025Prompt tables = the EGamma official SF (APPLIED:
//             2026-09-29, and again since 2026-10-01, user decision) -- see
//             the header
//   kRecoTnP  the pO tag-and-probe, ONE inclusive number (applied 2026-09-30
//             only; the pO cross-check: 0.988 vs 0.929 for the pp tables)
//   kRecoOff  factor 1 for every electron (the eleReco twins then equal the
//             nominal) -- the old kApplyEleRecoSF = false
//   kRecoUnit factor 1 for every electron (NO Reco SF) with a flat
//             +-kEleRecoUnitUnc per electron, shipped as its OWN nuisance
//             `eleReco` and left out of the combined eleSF (APPLIED since
//             2026-10-21, user decision: the pO data support a Reco SF ~ 1
//             (tag-and-probe 0.988-0.998) vs ~0.93 for the pp tables, so the
//             fit is to tell; one inclusive nuisance, no eta/pT dependence).
//             Z->ee: both legs move together, +-21% per pair.
enum EleRecoSource { kRecoOff = 0, kRecoEgm = 1, kRecoTnP = 2, kRecoUnit = 3 };
inline constexpr EleRecoSource kEleRecoSource = kRecoUnit;
inline constexpr double kEleRecoUnitUnc = 0.10; // kRecoUnit: +-1 sigma per electron (user: "~10%")
// the tag-and-probe result (correction/reco_sf_ele.C "fit"; run
// correction/run_reco_sf.sh) and which of its rows is applied: the nominal
// probe variant (noid = any supercluster), inclusive scheme, bin 0
// (paths since 2026-10-20: correction/electron_sf_check/<study>/rootfile/; until
// then correction/rootfile/reco_sf_ele/ and correction/rootfile/)
inline const char *kEleRecoTnPCandidates[] = {
    "../correction/electron_sf_check/reco_sf/rootfile/reco_sf_ele_result.root",
    "correction/electron_sf_check/reco_sf/rootfile/reco_sf_ele_result.root",
    nullptr};
inline const char *kEleRecoTnPHist = "recosf_noid_incl0";
inline constexpr double kEleRecoTnPEtMin = 20.0; // the probes' E_T^SC threshold (legs below it are counted)
inline const char *kEleTrigCandidates[] = {
    "../correction/electron_sf_check/trigger/rootfile/trig_eff_mb_ele.root",
    "correction/electron_sf_check/trigger/rootfile/trig_eff_mb_ele.root",
    nullptr};
inline const char *kEleTrigSel = "mt40"; // trig_eff_mb.C selection variant: nom | mt40

// The three sources map onto muon_sf.h's source slots 0/1/2 (EventSF::part/var,
// SFTwins::h) -- the machinery only needs three independent factors.
enum EleSource { kEReco = 0, kEIdIso = 1, kETrig = 2 };
static_assert(kNSources == 3, "electron_sf.h assumes three SF sources per event");
inline constexpr const char *kEleSFSourceNames[kNSources] = {"eleReco", "eleIDIso", "eleTrig"};
inline constexpr const char *kEleSFPartLabels[kNSources]  = {"<RECO>", "<IDISO>", "<TRIG>"};
inline constexpr const char *kEleSFCombinedName = "eleSF";
// Under kRecoUnit the Reco source is its own nuisance and stays OUT of the
// combined twins (muon_sf.h BookSFTwins combMask).
inline constexpr bool kEleRecoApart = (kEleRecoSource == kRecoUnit);
inline constexpr bool kEleSFCombMask[kNSources] = {!kEleRecoApart, true, true};
inline constexpr const char *kEleSFCombinedNote = kEleRecoApart
    ? "ID+iso, trigger shifts added in quadrature per bin (Reco = its own nuisance eleReco)"
    : "Reco, ID+iso, trigger shifts added in quadrature per bin";
// The nuisance(s) carried into the Combine inputs (plotting/mtandmet.C,
// dileptonpeak.C: <process>_<name>Up/Down) and the fork's `shape` rows
// (electron channels only): the combined electron-SF nuisance, plus -- under
// kRecoUnit -- the Reco one (= the per-source twins <h>_eleRecoUp/Down).
inline constexpr const char *kEleSFSystNames[2] = {"eleSF", "eleReco"};
inline constexpr int         kNEleSFSysts       = kEleRecoApart ? 2 : 1;

inline bool ElectronSFEnabled()
{
  const char *e = gSystem->Getenv("PO_ELE_SF");
  return !(e && (TString(e) == "off" || TString(e) == "0" || TString(e) == "OFF"));
}

// ============================================================================
// One EGM working point from its flattened CSV -> a Table2D (muon_sf.h)
// ============================================================================
// Rows: wp,eta_lo,eta_hi,pt_lo,pt_hi,sf,sfup,sfdown,... ('#' comments, one
// header row). The (eta, pT) grid must be COMPLETE -- every cell present
// exactly once, edges adjacent -- or the load fails: a hole would otherwise be
// filled silently by the clamping lookup. 'inf'/'-inf' edges are std::strtod's.
inline bool LoadEgmCsv(const std::string &path, const char *wp, Table2D &t,
                       std::string &err, std::string &provenance)
{
  std::ifstream in(path.c_str());
  if (!in) { err = "cannot open " + path; return false; }
  struct Cell { double elo, ehi, plo, phi, sf, up, dn; };
  std::vector<Cell> cells;
  std::string line;
  provenance.clear();
  while (std::getline(in, line))
  {
    if (line.empty()) continue;
    if (line[0] == '#')
    {
      if (line.rfind("# source:", 0) == 0) provenance = line.substr(2);
      continue;
    }
    if (line.rfind("wp,", 0) == 0) continue; // the column header
    std::stringstream ss(line);
    std::string f[12];
    int n = 0;
    while (n < 12 && std::getline(ss, f[n], ',')) ++n;
    if (n < 8) { err = path + ": short row '" + line + "'"; return false; }
    if (f[0] != wp) { err = path + ": row for working point '" + f[0] + "', expected '" + wp + "'"; return false; }
    Cell c;
    c.elo = std::strtod(f[1].c_str(), nullptr); c.ehi = std::strtod(f[2].c_str(), nullptr);
    c.plo = std::strtod(f[3].c_str(), nullptr); c.phi = std::strtod(f[4].c_str(), nullptr);
    c.sf  = std::strtod(f[5].c_str(), nullptr); c.up  = std::strtod(f[6].c_str(), nullptr);
    c.dn  = std::strtod(f[7].c_str(), nullptr);
    if (!(c.ehi > c.elo) || !(c.phi > c.plo)) { err = path + ": empty cell in row '" + line + "'"; return false; }
    cells.push_back(c);
  }
  if (cells.empty()) { err = path + ": no cells"; return false; }
  std::vector<double> ee, pe;
  for (const Cell &c : cells) { ee.push_back(c.elo); ee.push_back(c.ehi); pe.push_back(c.plo); pe.push_back(c.phi); }
  std::sort(ee.begin(), ee.end()); ee.erase(std::unique(ee.begin(), ee.end()), ee.end());
  std::sort(pe.begin(), pe.end()); pe.erase(std::unique(pe.begin(), pe.end()), pe.end());
  t = Table2D();
  t.name = wp;
  t.etaEdges = ee;
  t.ptEdges  = pe;
  const int ne = t.NEta(), np = t.NPt();
  if ((int)cells.size() != ne * np)
  { err = path + Form(": %d cells for a %d x %d grid", (int)cells.size(), ne, np); return false; }
  t.nom.assign(ne * np, NAN); t.up.assign(ne * np, NAN); t.dn.assign(ne * np, NAN);
  t.stat.assign(ne * np, NAN); t.syst.assign(ne * np, NAN);
  for (const Cell &c : cells)
  {
    const int ie = (int)(std::lower_bound(ee.begin(), ee.end(), c.elo) - ee.begin());
    const int ip = (int)(std::lower_bound(pe.begin(), pe.end(), c.plo) - pe.begin());
    if (ie + 1 >= (int)ee.size() || ee[ie + 1] != c.ehi || ip + 1 >= (int)pe.size() || pe[ip + 1] != c.phi)
    { err = path + Form(": cell eta [%g,%g) x pT [%g,%g) does not span adjacent edges", c.elo, c.ehi, c.plo, c.phi); return false; }
    const int k = ie * np + ip;
    if (!std::isnan(t.nom[k])) { err = path + Form(": duplicate cell eta [%g,%g) x pT [%g,%g)", c.elo, c.ehi, c.plo, c.phi); return false; }
    t.nom[k] = c.sf; t.up[k] = c.up; t.dn[k] = c.dn;
  }
  return true;
}

// ============================================================================
// The electron SF: Reco x ID+iso x trigger
// ============================================================================
class ElectronSF
{
public:
  Table2D reco20, reco75, idiso;
  std::string pathReco20, pathReco75, pathIdIso, provEgm;
  // the pO tag-and-probe Reco SF (kRecoTnP): the applied number, its +-1 sigma
  // and the record read with it (correction/reco_sf_ele.C, labeled bins)
  bool   recoTnPLoaded = false;
  double recoTnP = 1, recoTnPUp = 1, recoTnPDn = 1;           // applied: sf, sf + err_up, sf - err_dn
  double recoTnPStatLo = 0, recoTnPStatHi = 0, recoTnPSyst = 0;
  double recoTnPEpsData = 0, recoTnPEpsMC = 0, recoTnPEgm = 0, recoTnPCC = 0;
  double recoTnPNPass = 0, recoTnPNFail = 0;
  int    recoTnPStatus = -1, recoTnPCovQual = -1;
  std::string recoTnPSource;
  // the INCLUSIVE trigger numbers that are applied (pT > 25, all y)
  double sfI = 1, sfIUp = 1, sfIDn = 1;                  // W: fired && matched, data/MC
  double effDI = 1, effDIUp = 1, effDIDn = 1, effMI = 1; // per-lepton path-fired eff (Z formula)
  double nDenData = 0, nNumData = 0;                     // the data counts behind sfI
  // the per-|y| record (NOT applied): edges, SF and its CP errors
  std::vector<double> yEdges, sfY, sfYUp, sfYDn;
  std::string trigSource;
  mutable LookupCounters cntReco, cntIdIso, cntTrig;
  bool loaded = false;

  // Returns false (with [SF] ERR lines on `log`) when any input is missing or
  // malformed -- the caller treats that as FATAL, never as "SF = 1".
  // withRecoTnP = false loads everything BUT the tag-and-probe result: only for
  // correction/reco_sf_ele.C, which PRODUCES that result (and uses RecoEGM()
  // for its pp prediction); Reco() then refuses to run under kRecoTnP.
  bool Load(std::ostream &log, bool withRecoTnP = true)
  {
    loaded = false;
    std::string err, prov;
    auto load = [&](const char *wp, Table2D &t, std::string &path) -> bool
    {
      const std::string file = std::string("electron_sf_") + kEgmYear + "_" + wp + ".csv";
      path.clear();
      for (int i = 0; kEleCsvDirs[i]; ++i)
        if (!gSystem->AccessPathName((std::string(kEleCsvDirs[i]) + file).c_str(), kReadPermission))
        { path = std::string(kEleCsvDirs[i]) + file; break; }
      if (path.empty())
      { log << "[SF] ERR " << file << " not found (tried sf/, skim/sf/, ../skim/sf/); run skim/sf/extract_electron_sf.py\n"; return false; }
      if (!LoadEgmCsv(path, wp, t, err, prov)) { log << "[SF] ERR " << err << "\n"; return false; }
      if (provEgm.empty()) provEgm = prov;
      else if (prov != provEgm)
      { log << "[SF] ERR " << path << " was extracted from a different EGM file than the other tables:\n[SF]     "
            << prov << "\n[SF]  vs " << provEgm << "\n"; return false; }
      return true;
    };
    if (!load(kEleReco20WP, reco20, pathReco20)) return false;
    if (!load(kEleReco75WP, reco75, pathReco75)) return false;
    if (!load(kEleIdIsoWP, idiso, pathIdIso))    return false;
    // ---- the trigger study ----
    const char *tp = nullptr;
    for (int i = 0; kEleTrigCandidates[i]; ++i)
      if (!gSystem->AccessPathName(kEleTrigCandidates[i], kReadPermission)) { tp = kEleTrigCandidates[i]; break; }
    if (!tp) { log << "[SF] ERR trig_eff_mb_ele.root not found (tried ../correction/electron_sf_check/trigger/rootfile, correction/electron_sf_check/trigger/rootfile); run correction/run_trig_eff_mb.sh ele\n"; return false; }
    if (!LoadTrig(tp, log)) return false;
    // ---- the pO Reco tag-and-probe (the applied RECO factor under kRecoTnP) ----
    recoTnPLoaded = false;
    if (withRecoTnP && kEleRecoSource == kRecoTnP)
    {
      const char *rp = nullptr;
      for (int i = 0; kEleRecoTnPCandidates[i]; ++i)
        if (!gSystem->AccessPathName(kEleRecoTnPCandidates[i], kReadPermission)) { rp = kEleRecoTnPCandidates[i]; break; }
      if (!rp) { log << "[SF] ERR reco_sf_ele_result.root not found (tried ../correction/electron_sf_check/reco_sf/rootfile, correction/electron_sf_check/reco_sf/rootfile); run correction/run_reco_sf.sh\n"; return false; }
      if (!LoadRecoTnP(rp, log)) return false;
    }
    loaded = true;
    Print(log);
    return true;
  }

  // RECO, the APPLIED factor (kEleRecoSource). var: 0 nominal, +1 / -1 = +-1 sigma.
  //   kRecoTnP: one inclusive number, whatever eta and pT (legs below the probes'
  //             E_T^SC threshold are counted in cntReco.ptClamped)
  //   kRecoEgm: the pp tables (RecoEGM)
  double Reco(double scEta, double pt, int var) const
  {
    if (kEleRecoSource == kRecoOff) return 1.0;
    if (kEleRecoSource == kRecoUnit)
    {
      ++cntReco.calls;
      return var > 0 ? 1.0 + kEleRecoUnitUnc : (var < 0 ? 1.0 - kEleRecoUnitUnc : 1.0);
    }
    if (kEleRecoSource == kRecoEgm) return RecoEGM(scEta, pt, var);
    if (!recoTnPLoaded)
      throw std::runtime_error("pOSF::ElectronSF::Reco: kRecoTnP but the tag-and-probe result was not loaded "
                               "(Load(log, withRecoTnP = false) is for correction/reco_sf_ele.C only)");
    ++cntReco.calls;
    if (pt < kEleRecoTnPEtMin) ++cntReco.ptClamped;
    return var > 0 ? recoTnPUp : (var < 0 ? recoTnPDn : recoTnP);
  }

  // RECO from the pp EGM tables: Reco20to75 below 75 GeV, RecoAbove75 at/above
  // (pT < 20 clamps to the 20-45 GeV cell of Reco20to75 and is counted in
  // cntReco.ptClamped). The applied factor under kRecoEgm; otherwise the pp
  // PREDICTION of correction/reco_sf_ele.C, not a weight.
  double RecoEGM(double scEta, double pt, int var) const
  {
    return (pt < kEleRecoSplitPt ? reco20 : reco75).Eval(scEta, pt, var, cntReco);
  }

  // ---- W -> e nu: the leading electron (WP90 ID, matched); IDISO only if it passes the iso
  //
  // *** scEtaLab IS THE RAW SUPERCLUSTER eta (eleSCEta), NOT the analysis
  // rapidity y = -eta. *** See the header: the EGM tables are asymmetric in eta.
  EventSF W(double pt, double scEtaLab, bool passIso) const
  {
    const double eta = scEtaLab; // detector SC eta -> the EGM tables (NOT y = -eta)
    const double r  = Reco(eta, pt, 0), ru = Reco(eta, pt, +1), rd = Reco(eta, pt, -1);
    const double i  = passIso ? idiso.Eval(eta, pt, 0, cntIdIso) : 1.0;
    const double iu = passIso ? idiso.Eval(eta, pt, +1, cntIdIso) : 1.0;
    const double id = passIso ? idiso.Eval(eta, pt, -1, cntIdIso) : 1.0;
    ++cntTrig.calls;
    const double t = sfI, tu = sfIUp, td = sfIDn; // ONE inclusive number
    EventSF e;
    e.part[kEReco] = r; e.part[kEIdIso] = i; e.part[kETrig] = t;
    e.nom = r * i * t;
    e.var[kEReco][0]  = ru * i * t;  e.var[kEReco][1]  = rd * i * t;
    e.var[kEIdIso][0] = r * iu * t;  e.var[kEIdIso][1] = r * id * t;
    e.var[kETrig][0]  = r * i * tu;  e.var[kETrig][1]  = r * i * td;
    return e;
  }

  // ---- Z -> e e: Reco on both legs, IDISO on the legs passing the iso, the
  // standard OR trigger formula with the inclusive per-lepton path-fired
  // efficiencies (both legs share them: one inclusive number).
  // scEta1/scEta2 are RAW SUPERCLUSTER eta, not y -- see W().
  EventSF Z(double pt1, double scEtaLab1, bool iso1, double pt2, double scEtaLab2, bool iso2) const
  {
    const double eta1 = scEtaLab1, eta2 = scEtaLab2;
    const double a = Reco(eta1, pt1, 0), au = Reco(eta1, pt1, +1), ad = Reco(eta1, pt1, -1);
    const double b = Reco(eta2, pt2, 0), bu = Reco(eta2, pt2, +1), bd = Reco(eta2, pt2, -1);
    const double c  = iso1 ? idiso.Eval(eta1, pt1, 0, cntIdIso) : 1.0;
    const double cu = iso1 ? idiso.Eval(eta1, pt1, +1, cntIdIso) : 1.0;
    const double cd = iso1 ? idiso.Eval(eta1, pt1, -1, cntIdIso) : 1.0;
    const double d  = iso2 ? idiso.Eval(eta2, pt2, 0, cntIdIso) : 1.0;
    const double du = iso2 ? idiso.Eval(eta2, pt2, +1, cntIdIso) : 1.0;
    const double dd = iso2 ? idiso.Eval(eta2, pt2, -1, cntIdIso) : 1.0;
    ++cntTrig.calls;
    const double eM  = 1.0 - (1.0 - effMI) * (1.0 - effMI);
    const double eD  = 1.0 - (1.0 - effDI) * (1.0 - effDI);
    const double eDu = 1.0 - (1.0 - effDIUp) * (1.0 - effDIUp);
    const double eDd = 1.0 - (1.0 - effDIDn) * (1.0 - effDIDn);
    const double fT = eM > 0 ? eD / eM : 1.0, fTu = eM > 0 ? eDu / eM : 1.0, fTd = eM > 0 ? eDd / eM : 1.0;
    const double fR = a * b, fI = c * d;
    EventSF e;
    e.part[kEReco] = fR; e.part[kEIdIso] = fI; e.part[kETrig] = fT;
    e.nom = fR * fI * fT;
    e.var[kEReco][0]  = au * bu * fI * fT; e.var[kEReco][1]  = ad * bd * fI * fT;
    e.var[kEIdIso][0] = fR * cu * du * fT; e.var[kEIdIso][1] = fR * cd * dd * fT;
    e.var[kETrig][0]  = fR * fI * fTu;     e.var[kETrig][1]  = fR * fI * fTd;
    return e;
  }

  void PrintCounters(std::ostream &log, const char *who) const
  {
    log << "[SF] " << who << " lookups: RECO " << cntReco.calls;
    if (kEleRecoSource == kRecoTnP && recoTnPLoaded)
      log << " (pO TnP, inclusive; " << cntReco.ptClamped << " below the probes' E_T^SC > " << kEleRecoTnPEtMin
          << " GeV -- extrapolated)";
    else if (kEleRecoSource == kRecoUnit)
      log << Form(" (no Reco SF: factor 1 +- %.0f%%, inclusive)", 100.0 * kEleRecoUnitUnc);
    else
      log << " (eta clamped " << cntReco.etaClamped << ", pt clamped " << cntReco.ptClamped
          << " -- pT < 20 legs on the 20-45 GeV Reco cell)";
    log << ", IDISO " << cntIdIso.calls << " (eta " << cntIdIso.etaClamped << ", pt " << cntIdIso.ptClamped
        << "), TRIG " << cntTrig.calls << " (inclusive)\n";
  }

private:
  // The tag-and-probe result: one TH1D per (probe variant, scheme, bin) with
  // LABELED bins (correction/reco_sf_ele.C RecoSfFit). Labels are looked up by
  // looping GetBinLabel -- TAxis::FindBin(label) would APPEND a bin on a miss.
  // Refuses a nominal fit that did not converge.
  bool LoadRecoTnP(const char *path, std::ostream &log)
  {
    TDirectory::TContext ctx; // restore gDirectory afterwards
    TFile *f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) { log << "[SF] ERR cannot open " << path << "\n"; return false; }
    TH1 *h = (TH1 *)f->Get(kEleRecoTnPHist);
    if (!h) { log << "[SF] ERR " << path << " lacks " << kEleRecoTnPHist << " (re-run correction/run_reco_sf.sh fit)\n"; f->Close(); delete f; return false; }
    std::string missing;
    auto get = [&](const char *lab) -> double
    {
      for (int b = 1; b <= h->GetNbinsX(); ++b)
        if (std::string(h->GetXaxis()->GetBinLabel(b)) == lab) return h->GetBinContent(b);
      missing += std::string(missing.empty() ? "" : ", ") + lab;
      return NAN;
    };
    const double sf = get("sf"), lo = get("sf_stat_lo"), hi = get("sf_stat_hi"), sy = get("sf_syst");
    recoTnPEpsData = get("eps_data"); recoTnPEpsMC = get("eps_mc"); recoTnPEgm = get("egm_sf"); recoTnPCC = get("sf_cc");
    recoTnPNPass = get("n_pass_data"); recoTnPNFail = get("n_fail_data");
    const double st = get("fit_status"), cq = get("cov_qual");
    f->Close(); delete f;
    if (!missing.empty()) { log << "[SF] ERR " << path << ":" << kEleRecoTnPHist << " lacks the label(s) " << missing << "\n"; return false; }
    recoTnPStatus = (int)std::lround(st); recoTnPCovQual = (int)std::lround(cq);
    if (recoTnPStatus != 0)
    { log << "[SF] ERR " << path << ":" << kEleRecoTnPHist << ": the nominal fit did not converge (status " << recoTnPStatus << ") -- not applying it\n"; return false; }
    if (!(sf > 0.0) || !(lo > 0.0) || !(hi > 0.0) || !(sy >= 0.0))
    { log << "[SF] ERR " << path << ":" << kEleRecoTnPHist << ": unphysical SF/errors (" << sf << ", -" << lo << ", +" << hi << ", syst " << sy << ")\n"; return false; }
    if (recoTnPCovQual != 3) log << "[SF] WARN " << path << ":" << kEleRecoTnPHist << ": covQual " << recoTnPCovQual << " (not 3)\n";
    recoTnP = sf;
    recoTnPStatLo = lo; recoTnPStatHi = hi; recoTnPSyst = sy;
    recoTnPUp = sf + std::hypot(hi, sy);
    recoTnPDn = sf - std::hypot(lo, sy);
    recoTnPSource = std::string(path) + " : " + kEleRecoTnPHist + " (nominal probe = any supercluster, inclusive fit)";
    recoTnPLoaded = true;
    return true;
  }

  bool LoadTrig(const char *path, std::ostream &log)
  {
    TDirectory::TContext ctx; // restore gDirectory afterwards (the skim books histograms next)
    TFile *f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) { log << "[SF] ERR cannot open " << path << "\n"; return false; }
    const std::string sel = kEleTrigSel;
    auto get = [&](const std::string &n) { return (TH1 *)f->Get(n.c_str()); };
    TH1 *dDen = get("h_den_absy_" + sel + "_data_ele"), *dNum = get("h_num_absy_" + sel + "_data_ele"),
        *dBit = get("h_bit_absy_" + sel + "_data_ele");
    TH1 *pDen = get("h_den_absy_" + sel + "_Wp_ele"), *pNum = get("h_num_absy_" + sel + "_Wp_ele"),
        *pBit = get("h_bit_absy_" + sel + "_Wp_ele");
    TH1 *mDen = get("h_den_absy_" + sel + "_Wm_ele"), *mNum = get("h_num_absy_" + sel + "_Wm_ele"),
        *mBit = get("h_bit_absy_" + sel + "_Wm_ele");
    TGraphAsymmErrors *g = (TGraphAsymmErrors *)f->Get(("sf_absy_" + sel).c_str());
    if (!dDen || !dNum || !dBit || !pDen || !pNum || !pBit || !mDen || !mNum || !mBit || !g)
    {
      log << "[SF] ERR " << path << " lacks sf_absy_" << sel << " and/or the h_{den,num,bit}_absy_" << sel
          << "_{data,Wp,Wm}_ele histograms (re-run correction/run_trig_eff_mb.sh ele)\n";
      f->Close(); delete f; return false;
    }
    const int n = dDen->GetNbinsX();
    double sDd = 0, sDn = 0, sDb = 0, sMd = 0, sMn = 0, sMb = 0;
    for (int i = 1; i <= n; ++i)
    {
      sDd += dDen->GetBinContent(i); sDn += dNum->GetBinContent(i); sDb += dBit->GetBinContent(i);
      sMd += pDen->GetBinContent(i) + mDen->GetBinContent(i);
      sMn += pNum->GetBinContent(i) + mNum->GetBinContent(i);
      sMb += pBit->GetBinContent(i) + mBit->GetBinContent(i);
    }
    if (sDd <= 0 || sMd <= 0 || sMn <= 0) { log << "[SF] ERR empty trigger denominators in " << path << "\n"; f->Close(); delete f; return false; }
    const double eD = sDn / sDd, eM = sMn / sMd;
    nDenData = sDd; nNumData = sDn;
    sfI   = eD / eM;
    sfIUp = TEfficiency::ClopperPearson((int)sDd, (int)sDn, 0.683, true)  / eM;
    sfIDn = TEfficiency::ClopperPearson((int)sDd, (int)sDn, 0.683, false) / eM;
    effDI   = sDb / sDd;
    effDIUp = TEfficiency::ClopperPearson((int)sDd, (int)sDb, 0.683, true);
    effDIDn = TEfficiency::ClopperPearson((int)sDd, (int)sDb, 0.683, false);
    effMI   = sMb / sMd;
    yEdges.clear(); sfY.clear(); sfYUp.clear(); sfYDn.clear();
    for (int i = 0; i < g->GetN(); ++i)
    {
      yEdges.push_back(g->GetX()[i] - g->GetEXlow()[i]);
      if (i == g->GetN() - 1) yEdges.push_back(g->GetX()[i] + g->GetEXhigh()[i]);
      sfY.push_back(g->GetY()[i]);
      sfYUp.push_back(g->GetEYhigh()[i]);
      sfYDn.push_back(g->GetEYlow()[i]);
    }
    trigSource = std::string(path) + " : h_{den,num,bit}_absy_" + sel + "_{data,Wp,Wm}_ele summed over all |y| bins (inclusive); sf_absy_" + sel + " = the per-|y| record";
    f->Close(); delete f;
    return true;
  }

  void Print(std::ostream &log) const
  {
    auto rng = [](const Table2D &t, double ptMin, double ptMax) -> std::string
    {
      double lo = 9, hi = -9;
      for (int ie = 0; ie < t.NEta(); ++ie)
        for (int ip = 0; ip < t.NPt(); ++ip)
        {
          const double mid = 0.5 * (std::max(t.etaEdges[ie], -3.0) + std::min(t.etaEdges[ie + 1], 3.0));
          if (std::fabs(mid) > 1.4442 && std::fabs(mid) < 1.566) continue; // crack placeholders
          if (t.ptEdges[ip + 1] <= ptMin || t.ptEdges[ip] >= ptMax) continue;
          const double v = t.nom[ie * t.NPt() + ip];
          lo = std::min(lo, v); hi = std::max(hi, v);
        }
      return Form("%.4f..%.4f", lo, hi);
    };
    log << "[SF] electron scale factors ON (MC event weight x RECO x ID+ISO x TRIG; PO_ELE_SF=off disables)\n"
        << "[SF]   EGM   : " << provEgm << "\n";
    if (kEleRecoSource == kRecoOff)
      log << "[SF]   RECO  : *** NOT APPLIED (kEleRecoSource = kRecoOff): factor 1 for every electron ***\n";
    else if (kEleRecoSource == kRecoUnit)
      log << "[SF]   RECO  : *** NO Reco SF (kEleRecoSource = kRecoUnit): factor 1 for every electron, "
          << Form("+-%.0f%% per electron (Z->ee: both legs coherently)", 100.0 * kEleRecoUnitUnc)
          << " = its OWN nuisance '" << kEleSFSourceNames[kEReco] << "', NOT in the combined '"
          << kEleSFCombinedName << "' ***\n";
    else if (kEleRecoSource == kRecoTnP && recoTnPLoaded)
      log << "[SF]   RECO  : APPLIED = the pO Z->ee TAG-AND-PROBE, ONE INCLUSIVE NUMBER for every electron (kEleRecoSource = kRecoTnP):\n"
          << "[SF]           " << recoTnPSource << "\n"
          << "[SF]           SF = " << Form("%.5f +%.5f -%.5f (stat) +- %.5f (syst)", recoTnP, recoTnPStatHi, recoTnPStatLo, recoTnPSyst)
          << Form("  ->  applied +-1 sigma: %.5f / %.5f", recoTnPUp, recoTnPDn) << "\n"
          << "[SF]           eps_data " << Form("%.4f", recoTnPEpsData) << " vs eps_MC " << Form("%.4f", recoTnPEpsMC)
          << Form(" [%.0f pass / %.0f fail pairs in data; fit status %d, covQual %d]", recoTnPNPass, recoTnPNFail, recoTnPStatus, recoTnPCovQual)
          << "; cut-and-count 80-100 GeV " << Form("%.4f", recoTnPCC) << "; pp EGM on the same probes " << Form("%.4f", recoTnPEgm)
          << " (NOT applied)\n";
    else if (kEleRecoSource == kRecoTnP)
      log << "[SF]   RECO  : kRecoTnP, tag-and-probe result NOT LOADED (Load(log, false): correction/reco_sf_ele.C) -- Reco() unusable, RecoEGM() only\n";
    log << "[SF]   RECO  : " << (kEleRecoSource == kRecoEgm ? "APPLIED = the pp EGM tables: " : "pp EGM tables (the record, NOT applied): ")
        << pathReco20 << " (pT < " << kEleRecoSplitPt << "; nominal range 20-75 GeV outside the crack: "
        << rng(reco20, 20.0, kEleRecoSplitPt) << ")\n"
        << "[SF]           " << pathReco75 << " (pT >= " << kEleRecoSplitPt << "; 75-500 GeV: " << rng(reco75, kEleRecoSplitPt, 500.0) << ")\n"
        << "[SF]   IDISO : " << pathIdIso << " (eleMVAIdWP90 && eleMVAIsoWP90 && relIso < 0.3; nominal range at 25-100 GeV: "
        << rng(idiso, 25.0, 100.0) << ", 10-25 GeV (Z legs): " << rng(idiso, 10.0, 25.0) << ")\n"
        << "[SF]   TRIG  : " << trigSource << "\n"
        << "[SF]           APPLIED INCLUSIVE (one number for every electron): SF(W: fired&&matched) = "
        << Form("%.4f (+%.4f -%.4f)", sfI, sfIUp - sfI, sfI - sfIDn)
        << Form("  [%.0f/%.0f in data]", nNumData, nDenData)
        << ";  per-lepton eps(fired) data " << Form("%.4f (+%.4f -%.4f)", effDI, effDIUp - effDI, effDI - effDIDn)
        << ", MC " << Form("%.4f", effMI) << "\n";
    if (!sfY.empty())
    {
      log << "[SF]           per-|y| record (MEASURED, NOT APPLIED):";
      for (size_t i = 0; i < sfY.size(); ++i)
        log << Form("%s[%.1f,%.1f) %.4f +%.4f -%.4f", i % 3 == 0 ? "\n[SF]             " : "   ",
                    yEdges[i], yEdges[i + 1], sfY[i], sfYUp[i], sfYDn[i]);
      log << "\n";
    }
    log << "[SF]   per-source twins <hist>_" << kEleSFSourceNames[0] << "/" << kEleSFSourceNames[1] << "/"
        << kEleSFSourceNames[2] << "Up|Down (one source at +-1 sigma, others nominal) + the COMBINED <hist>_"
        << kEleSFCombinedName << "Up|Down (per-bin quadrature of "
        << (kEleRecoApart ? "ID+iso and trigger" : "the three") << ")";
    if (kEleRecoApart)
      log << "; the fit carries TWO electron nuisances: " << kEleSFCombinedName << " and "
          << kEleSFSourceNames[kEReco] << " (= the per-source twins)\n";
    else
      log << " = the one nuisance carried into the fit\n";
  }
};

} // namespace pOSF

#endif // PO_ELECTRON_SF_H
