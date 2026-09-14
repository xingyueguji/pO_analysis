// skim/muon_sf.h -- muon efficiency SCALE FACTORS (data/MC) for the MC event
// weight of skim_Wmu / skim_Zmm (2026-09-14, the SF-application phase).
//
// SINGLE SOURCE for: which SFs are applied, where they come from, how they are
// looked up, how the per-event factor and its +-1 sigma variations are formed,
// and the names of the resulting shape systematics. Nothing about the SFs is
// hardcoded anywhere else -- skim.C only calls W() / Z() and fills the twins.
//
// WHAT IS APPLIED (MC only; data is always weight 1):
//   ID   : NUM_TightID_DEN_TrackerMuons  -- the skim's muIDTight (CutBasedIdTight)
//   ISO  : NUM_TightPFIso_DEN_TightID    -- PF relIso(dBeta, R=0.4) < 0.15 given
//          TightID = the isolation cut of the W skim AND, since 2026-09-14, of the
//          Z skim (skim_Zmm harmonized from 0.20; skim_common.h RelIsoPF)
//          Both from the Muon POG pp 2025 correctionlib file (schema v2), reduced
//          to the used corrections by sf/extract_muon_sf.py -> sf/muon_sf_2025_*.json
//          (provenance paragraph inside the file's description). Binned in SIGNED
//          eta (24 x 0.2) and pT [10,15,20,25,30,40,50,60,120,inf); 'nominal' and
//          the XPOG 'systup'/'systdown' = nominal +- sqrt(stat^2 + syst^2).
//          The pp SFs are a TENTATIVE stand-in for pO tag-and-probe (user decision
//          2026-09-14): pO has no pile-up, so the ISO SF in particular corrects a
//          different mismodelling than it was derived for.
//   TRIG : the analysis path HLT_OxyL1SingleMuOpen_v1 (fired AND the leading muon
//          matched within dR<0.4 = skim step 8), measured by correction/trig_eff_mb.C
//          from the minimum-bias-triggered data vs the W signal MC, pT > 25. The
//          `mt40` selection variant is the nominal (the purer W sample: ~5% fakes
//          vs ~19% in the plain selection; 0.9971 vs 0.9962 inclusive). Applied
//          INCLUSIVELY in rapidity (kTrigBinning = kTrigInclusive, user decision
//          2026-09-14): the study's 12 per-y values rest on 170-245 data events
//          each (+-0.8-1.9% per bin vs +-0.2% inclusive) and are consistent with
//          one flat SF -- chi2/ndf = 11.5/11 (p = 0.41) for mt40, the |y| < 0.4 dip
//          is a 1.7 sigma effect -- so a per-y SF would only inject its own
//          statistical noise into the per-bin r's. The per-y table and the chi2
//          are still loaded and printed by every job as the standing check;
//          kTrigPerY switches back (the fork then decorrelates muTrig per y bin --
//          the correlation model travels in the sidecar directive `#! muTrig corr`,
//          see kMuTrigCorr). Errors = Clopper-Pearson 68% on the data count (pure
//          statistics; the MC error is negligible) -> slightly ASYMMETRIC. Read
//          from the study's rootfile (h_num/h_den_y_mt40 summed for the inclusive,
//          sf_y_mt40 for the per-y check) -- run correction/run_trig_eff_mb.sh mu
//          first.
//
// HOW IT ENTERS (skim.C):
//   W -> mu nu : one factor per event for the LEADING muon after step 8,
//                SF = ID(pt,eta) x [passIso ? ISO(pt,eta) : 1] x TRIG(y).
//                Folded into the weight of EVERY fill (all regions of the ABCD
//                planes have the leading muon TightID'd + matched; the anti-iso
//                sideband has no measured SF -> 1, its EWK content is a few %).
//   Z -> mu mu : both legs, SF = prod_legs ID x [iso ? ISO : 1], times the
//                per-EVENT trigger factor [1 - prod(1 - eps_data,i)] /
//                [1 - prod(1 - eps_MC,i)] with the per-lepton PATH-FIRED
//                efficiencies of the same study (skim_Zmm requires the bit, no
//                matching) -- numerically 1 +- 1e-4 with eps ~ 0.99, kept for
//                consistency. skim_Zmm's iso cut was 0.20 (the POG Medium WP, for
//                which the pp file has no TightID-denominator table) and was
//                HARMONIZED to the W's 0.15 on 2026-09-14 (user decision), so the
//                one Tight|TightID ISO SF serves both channels. (Checked then:
//                Tight|TightID vs Tight|MediumID agree to <= 0.2% per bin, so the
//                denominator would not have mattered either way.)
//   NOT covered : the reco/tracking efficiency (no SF supplied, ~1), the DY-veto
//                second lepton in the W selection, and the 1e-4 residual
//                correlation of the Z trigger factor with the W bins.
//
// UNCERTAINTIES: three independent sources (ID, ISO, trigger), each varied
// +-1 sigma with the other two nominal (EventSF::var), stored per fit-template
// histogram as the per-SOURCE twins <h>_muIDUp/Down, <h>_muIsoUp/Down,
// <h>_muTrigUp/Down (filled with w_gen x var, NOT area-normalized: the
// normalization IS the effect) -- diagnostics and the ingredients of
// **the ONE nuisance that reaches the fit (user decision 2026-09-14): <h>_muSFUp/Down,
// built in FinalizeSFTwins() bin by bin as nominal +- sqrt(sum_s (shift_s)^2)** --
// the quadrature sum of the three independent per-bin shifts = the exact 1 sigma
// of the product of independent factors (adding them linearly would assume full
// correlation). One parameter then moves every bin, process, charge and channel
// of the muon inputs coherently (the three factors are near-flat normalization
// factors, so a single coherent nuisance loses nothing; ID/ISO use the XPOG
// systup/systdown totals, muTrig the inclusive CP error). kMuonSFSystNames =
// {"muSF"} is what mtandmet.C / dileptonpeak.C carry into the Combine inputs
// and the fork turns into one `shape` row (group lepsf). kMuTrigCorr (the
// sidecar directive `#! muTrig corr`) only matters if the three sources are ever
// shipped separately again (kMuonSFSystNames = kMuonSFSourceNames): the fork
// would then split a per-y muTrig into 12 nuisances via `nuisance edit rename`.
//
// Because the SF is part of the nominal weight, the LHE member weights (and
// hence nPDF/qcdScale/alphaS Up/Down) are computed on the SF-weighted nominal
// (pOLhe::ScaleMemberWeights) -- every variation is one-at-a-time around the
// same nominal, which is how Combine composes them.
//
// Env switch: PO_MUON_SF=off fills MC without SFs (regression checks); the job
// log always says which.
#ifndef PO_MUON_SF_H
#define PO_MUON_SF_H

#include "TDirectory.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace pOSF
{

// ============================================================================
// Configuration -- edit here, nowhere else
// ============================================================================
// Candidate paths (tried in order; relative to the cwd -- skim/ when run through
// run_all.sh, the repo root or a sibling directory otherwise).
inline const char *kJsonCandidates[] = {
    "sf/muon_sf_2025_TightID_PFIso_schemaV2.json",
    "skim/sf/muon_sf_2025_TightID_PFIso_schemaV2.json",
    "../skim/sf/muon_sf_2025_TightID_PFIso_schemaV2.json",
    nullptr};
inline const char *kIDCorrection  = "NUM_TightID_DEN_TrackerMuons";
inline const char *kIsoCorrection = "NUM_TightPFIso_DEN_TightID"; // relIso < 0.15 | TightID = the W AND the Z cut
inline const char *kTrigCandidates[] = {
    "../correction/rootfile/trig_eff_mb_mu.root",
    "correction/rootfile/trig_eff_mb_mu.root",
    nullptr};
inline const char *kTrigSel = "mt40"; // trig_eff_mb.C selection variant: nom | mt40
// Rapidity binning of the trigger SF (header comment): inclusive by decision.
enum TrigBinning { kTrigInclusive = 0, kTrigPerY = 1 };
inline constexpr TrigBinning kTrigBinning = kTrigInclusive;
// The correlation model the fork must use for the muTrig nuisance, written into
// the Combine-input sidecars as `#! muTrig corr <coherent|perbin>` by
// plotting/mtandmet.C and dileptonpeak.C: one inclusive factor = one coherent
// nuisance; a per-y factor = per-y-bin nuisances (pure-statistics errors).
inline constexpr const char *kMuTrigCorr = (kTrigBinning == kTrigInclusive) ? "coherent" : "perbin";

enum Source { kID = 0, kIso = 1, kTrig = 2, kNSources = 3 };
// Per-SOURCE twins in the MC skim files: <hist>_<name>Up/Down (diagnostics +
// the ingredients of the combined one).
inline constexpr const char *kMuonSFSourceNames[kNSources] = {"muID", "muIso", "muTrig"};
// The nuisance(s) carried into the Combine inputs (plotting/mtandmet.C,
// dileptonpeak.C: <process>_<name>Up/Down) and the fork's `shape` rows (muon
// channels only): ONE combined muon-SF nuisance (see the header comment).
inline constexpr const char *kMuonSFCombinedName = "muSF";
inline constexpr const char *kMuonSFSystNames[1] = {"muSF"};
inline constexpr int         kNMuonSFSysts       = 1;

inline bool MuonSFEnabled()
{
  const char *e = gSystem->Getenv("PO_MUON_SF");
  return !(e && (TString(e) == "off" || TString(e) == "0" || TString(e) == "OFF"));
}

// ============================================================================
// Minimal JSON reader (the correctionlib schema-v2 subset we need; no external
// dependency -- ROOT 6.32 ships no generic JSON parser and correctionlib is not
// installed for the ROOT python).
// ============================================================================
struct JVal
{
  enum Type { kNull, kBool, kNumber, kString, kArray, kObject };
  Type                                    type = kNull;
  bool                                    b    = false;
  double                                  num  = 0.0;
  std::string                             str;
  std::vector<JVal>                       arr;
  std::vector<std::pair<std::string, JVal>> obj;

  const JVal *Get(const char *key) const
  {
    if (type != kObject) return nullptr;
    for (const auto &kv : obj)
      if (kv.first == key) return &kv.second;
    return nullptr;
  }
  size_t Size() const { return type == kArray ? arr.size() : (type == kObject ? obj.size() : 0); }
  // A number, or the strings "inf"/"-inf" (correctionlib bin edges)
  double Num() const
  {
    if (type == kNumber) return num;
    if (type == kString)
    {
      if (str == "inf" || str == "+inf" || str == "Infinity") return std::numeric_limits<double>::infinity();
      if (str == "-inf" || str == "-Infinity") return -std::numeric_limits<double>::infinity();
    }
    return std::numeric_limits<double>::quiet_NaN();
  }
  bool IsStr(const char *s) const { return type == kString && str == s; }
};

class JParser
{
public:
  explicit JParser(const std::string &text) : s_(text) {}
  bool Parse(JVal &root, std::string &err)
  {
    p_ = 0;
    if (!Value(root)) { err = err_ + Form(" (at byte %zu)", p_); return false; }
    Ws();
    if (p_ != s_.size()) { err = Form("trailing characters at byte %zu", p_); return false; }
    return true;
  }

private:
  const std::string &s_;
  size_t             p_ = 0;
  std::string        err_;

  void Ws() { while (p_ < s_.size() && (s_[p_] == ' ' || s_[p_] == '\t' || s_[p_] == '\n' || s_[p_] == '\r')) ++p_; }
  bool Fail(const char *m) { if (err_.empty()) err_ = m; return false; }
  bool Value(JVal &v)
  {
    Ws();
    if (p_ >= s_.size()) return Fail("unexpected end of input");
    const char c = s_[p_];
    if (c == '{') return Object(v);
    if (c == '[') return Array(v);
    if (c == '"') { v.type = JVal::kString; return String(v.str); }
    if (c == 't' && s_.compare(p_, 4, "true") == 0)  { v.type = JVal::kBool; v.b = true;  p_ += 4; return true; }
    if (c == 'f' && s_.compare(p_, 5, "false") == 0) { v.type = JVal::kBool; v.b = false; p_ += 5; return true; }
    if (c == 'n' && s_.compare(p_, 4, "null") == 0)  { v.type = JVal::kNull; p_ += 4; return true; }
    if (c == '-' || (c >= '0' && c <= '9')) { v.type = JVal::kNumber; return Number(v.num); }
    return Fail("unexpected character");
  }
  bool Object(JVal &v)
  {
    v.type = JVal::kObject;
    ++p_; // '{'
    Ws();
    if (p_ < s_.size() && s_[p_] == '}') { ++p_; return true; }
    for (;;)
    {
      Ws();
      if (p_ >= s_.size() || s_[p_] != '"') return Fail("expected object key");
      std::string key;
      if (!String(key)) return false;
      Ws();
      if (p_ >= s_.size() || s_[p_] != ':') return Fail("expected ':'");
      ++p_;
      v.obj.emplace_back(std::move(key), JVal());
      if (!Value(v.obj.back().second)) return false;
      Ws();
      if (p_ >= s_.size()) return Fail("unterminated object");
      if (s_[p_] == ',') { ++p_; continue; }
      if (s_[p_] == '}') { ++p_; return true; }
      return Fail("expected ',' or '}'");
    }
  }
  bool Array(JVal &v)
  {
    v.type = JVal::kArray;
    ++p_; // '['
    Ws();
    if (p_ < s_.size() && s_[p_] == ']') { ++p_; return true; }
    for (;;)
    {
      v.arr.emplace_back();
      if (!Value(v.arr.back())) return false;
      Ws();
      if (p_ >= s_.size()) return Fail("unterminated array");
      if (s_[p_] == ',') { ++p_; continue; }
      if (s_[p_] == ']') { ++p_; return true; }
      return Fail("expected ',' or ']'");
    }
  }
  bool String(std::string &out)
  {
    ++p_; // opening quote
    out.clear();
    while (p_ < s_.size())
    {
      const char c = s_[p_++];
      if (c == '"') return true;
      if (c != '\\') { out += c; continue; }
      if (p_ >= s_.size()) return Fail("bad escape");
      const char e = s_[p_++];
      switch (e)
      {
        case '"': out += '"'; break;
        case '\\': out += '\\'; break;
        case '/': out += '/'; break;
        case 'b': out += '\b'; break;
        case 'f': out += '\f'; break;
        case 'n': out += '\n'; break;
        case 'r': out += '\r'; break;
        case 't': out += '\t'; break;
        case 'u': // keep the escape verbatim (only appears in free text here)
          out += "\\u";
          break;
        default: return Fail("unknown escape");
      }
    }
    return Fail("unterminated string");
  }
  bool Number(double &out)
  {
    const char *start = s_.c_str() + p_;
    char       *end   = nullptr;
    out = std::strtod(start, &end);
    if (end == start) return Fail("bad number");
    p_ += (size_t)(end - start);
    return true;
  }
};

// ============================================================================
// One correction = a (signed eta) x (pT) table of nominal / systup / systdown
// ============================================================================
struct LookupCounters
{
  unsigned long long calls = 0, etaClamped = 0, ptClamped = 0;
};

// lower-inclusive bin [e_i, e_{i+1}); below range -> 0, at/above the last
// edge -> last bin (counted as clamped -- correctionlib would raise here)
inline int FindEdgeBin(const std::vector<double> &edges, double x, unsigned long long &clamped)
{
  const int n = (int)edges.size() - 1;
  if (x < edges[0]) { ++clamped; return 0; }
  for (int i = 0; i < n; ++i)
    if (x >= edges[i] && x < edges[i + 1]) return i;
  ++clamped;
  return n - 1;
}

struct Table2D
{
  std::string         name;
  std::vector<double> etaEdges, ptEdges;
  std::vector<double> nom, up, dn, stat, syst; // [ie * npt + ip]

  int  NPt() const { return (int)ptEdges.size() - 1; }
  int  NEta() const { return (int)etaEdges.size() - 1; }
  bool Empty() const { return nom.empty(); }

  // Parse one member of the file's "corrections" list (binning eta -> binning pt
  // -> category scale_factors). Any deviation from that layout is an error.
  bool Load(const JVal &corrections, const char *cname, std::string &err)
  {
    const JVal *corr = nullptr;
    for (const JVal &c : corrections.arr)
      if (const JVal *n = c.Get("name")) if (n->IsStr(cname)) { corr = &c; break; }
    if (!corr) { err = std::string("correction '") + cname + "' not in the file"; return false; }
    name = cname;
    // inputs must be (eta, pt, scale_factors) in this order
    const JVal *inputs = corr->Get("inputs");
    const char *want[3] = {"eta", "pt", "scale_factors"};
    bool inOk = inputs && inputs->Size() == 3;
    for (int i = 0; inOk && i < 3; ++i)
    { const JVal *nm = inputs->arr[i].Get("name"); inOk = nm && nm->IsStr(want[i]); }
    if (!inOk) { err = name + ": inputs are not (eta, pt, scale_factors)"; return false; }
    // node type / input / edges of a binning node (nullptr-safe)
    auto isBinningOn = [](const JVal &node, const char *input) -> bool
    {
      const JVal *nt = node.Get("nodetype"), *in = node.Get("input"), *ed = node.Get("edges"), *co = node.Get("content");
      return nt && nt->IsStr("binning") && in && in->IsStr(input) && ed && ed->type == JVal::kArray && co && co->type == JVal::kArray;
    };
    const JVal *d = corr->Get("data");
    if (!d || !isBinningOn(*d, "eta")) { err = name + ": data is not a binning on eta"; return false; }
    for (const JVal &e : d->Get("edges")->arr) etaEdges.push_back(e.Num());
    const JVal *econt = d->Get("content");
    if (!econt || (int)econt->Size() != NEta()) { err = name + ": eta content size mismatch"; return false; }
    for (int ie = 0; ie < NEta(); ++ie)
    {
      const JVal &pb = econt->arr[ie];
      if (!isBinningOn(pb, "pt")) { err = name + Form(": eta bin %d is not a binning on pt", ie); return false; }
      std::vector<double> pe;
      for (const JVal &e : pb.Get("edges")->arr) pe.push_back(e.Num());
      if (ie == 0) ptEdges = pe;
      else if (pe != ptEdges) { err = name + ": pt edges differ between eta bins"; return false; }
      const JVal *pcont = pb.Get("content");
      if (!pcont || (int)pcont->Size() != NPt()) { err = name + ": pt content size mismatch"; return false; }
      for (int ip = 0; ip < NPt(); ++ip)
      {
        const JVal &cat = pcont->arr[ip];
        const JVal *cnt = cat.Get("nodetype"), *cin = cat.Get("input"), *cco = cat.Get("content");
        if (!cnt || !cnt->IsStr("category") || !cin || !cin->IsStr("scale_factors") || !cco || cco->type != JVal::kArray)
        { err = name + Form(": (eta %d, pt %d) is not a category on scale_factors", ie, ip); return false; }
        double vn = NAN, vu = NAN, vd = NAN, vs = NAN, vy = NAN;
        for (const JVal &kv : cco->arr)
        {
          const JVal *k = kv.Get("key"), *v = kv.Get("value");
          if (!k || !v) continue;
          if (k->IsStr("nominal"))  vn = v->Num();
          if (k->IsStr("systup"))   vu = v->Num();
          if (k->IsStr("systdown")) vd = v->Num();
          if (k->IsStr("stat"))     vs = v->Num();
          if (k->IsStr("syst"))     vy = v->Num();
        }
        if (std::isnan(vn) || std::isnan(vu) || std::isnan(vd))
        { err = name + Form(": (eta %d, pt %d) lacks nominal/systup/systdown", ie, ip); return false; }
        nom.push_back(vn); up.push_back(vu); dn.push_back(vd); stat.push_back(vs); syst.push_back(vy);
      }
    }
    return true;
  }

  // var: 0 = nominal, +1 = systup, -1 = systdown
  double Eval(double eta, double pt, int var, LookupCounters &c) const
  {
    ++c.calls;
    const int ie = FindEdgeBin(etaEdges, eta, c.etaClamped);
    const int ip = FindEdgeBin(ptEdges, pt, c.ptClamped);
    const int k  = ie * NPt() + ip;
    return var > 0 ? up[k] : (var < 0 ? dn[k] : nom[k]);
  }

  // min..max of the nominal over the pT bins at/above ptMin (for the log)
  std::pair<double, double> Range(double ptMin) const
  {
    double lo = 9, hi = -9;
    for (int ie = 0; ie < NEta(); ++ie)
      for (int ip = 0; ip < NPt(); ++ip)
        if (ptEdges[ip + 1] > ptMin)
        { const double v = nom[ie * NPt() + ip]; lo = std::min(lo, v); hi = std::max(hi, v); }
    return {lo, hi};
  }
};

// ============================================================================
// Trigger SF per rapidity bin (W: fired && matched) + per-lepton path-fired
// efficiencies (Z), from correction/rootfile/trig_eff_mb_mu.root
// ============================================================================
struct TrigTable
{
  std::vector<double> yEdges;            // kNY + 1
  std::vector<double> sf, sfUp, sfDn;    // per-y W factor (fired && matched), asymmetric errors applied
  std::vector<double> effD, effDUp, effDDn, effM; // per-y per-lepton path-fired eff (Z formula)
  double sfI = 1, sfIUp = 1, sfIDn = 1;           // the INCLUSIVE W factor (pT > 25, all y)
  double effDI = 1, effDIUp = 1, effDIDn = 1, effMI = 1; // inclusive per-lepton path-fired eff
  double chi2 = 0; int ndf = 0;                   // per-y SF vs the inclusive one (the standing check)
  std::string         source;

  int NBins() const { return (int)sf.size(); }
  int Bin(double y, unsigned long long &clamped) const { return FindEdgeBin(yEdges, y, clamped); }
  // what the skim applies, by kTrigBinning
  double SF(int b) const    { return kTrigBinning == kTrigInclusive ? sfI    : sf[b]; }
  double SFUp(int b) const  { return kTrigBinning == kTrigInclusive ? sfIUp  : sfUp[b]; }
  double SFDn(int b) const  { return kTrigBinning == kTrigInclusive ? sfIDn  : sfDn[b]; }
  double EffD(int b) const  { return kTrigBinning == kTrigInclusive ? effDI  : effD[b]; }
  double EffDUp(int b) const{ return kTrigBinning == kTrigInclusive ? effDIUp: effDUp[b]; }
  double EffDDn(int b) const{ return kTrigBinning == kTrigInclusive ? effDIDn: effDDn[b]; }
  double EffM(int b) const  { return kTrigBinning == kTrigInclusive ? effMI  : effM[b]; }
};

// ============================================================================
// The per-event factor and its one-source-at-a-time variations
// ============================================================================
struct EventSF
{
  double nom = 1.0;                 // the factor folded into the event weight
  double part[kNSources]  = {1, 1, 1}; // nominal factor of each source (diagnostics)
  double var[kNSources][2] = {{1, 1}, {1, 1}, {1, 1}}; // TOTAL factor with source s Up / Down
  static EventSF Unit() { return EventSF(); }
};

class MuonSF
{
public:
  Table2D        id, iso; // the ID and ISO tables (the same 0.15 iso cut in W and Z)
  TrigTable      trig;
  std::string    jsonPath, jsonProvenance;
  mutable LookupCounters cntID, cntIso, cntTrig;
  bool           loaded = false;

  // Returns false (with [SF] ERR lines on `log`) when any input is missing or
  // malformed -- the caller should treat that as FATAL, never as "SF = 1".
  bool Load(std::ostream &log)
  {
    loaded = false;
    // ---- the POG JSON ----
    const char *jp = FindFile(kJsonCandidates);
    if (!jp) { log << "[SF] ERR muon SF JSON not found (tried sf/, skim/sf/, ../skim/sf/); run skim/sf/extract_muon_sf.py\n"; return false; }
    std::ifstream in(jp, std::ios::binary);
    std::stringstream ss; ss << in.rdbuf();
    const std::string text = ss.str();
    JVal root; std::string err;
    if (!JParser(text).Parse(root, err)) { log << "[SF] ERR cannot parse " << jp << ": " << err << "\n"; return false; }
    const JVal *sv = root.Get("schema_version");
    if (!sv || sv->Num() != 2) { log << "[SF] ERR " << jp << ": schema_version != 2\n"; return false; }
    const JVal *corr = root.Get("corrections");
    if (!corr || corr->type != JVal::kArray) { log << "[SF] ERR " << jp << ": no corrections list\n"; return false; }
    if (!id.Load(*corr, kIDCorrection, err))   { log << "[SF] ERR " << err << "\n"; return false; }
    if (!iso.Load(*corr, kIsoCorrection, err)) { log << "[SF] ERR " << err << "\n"; return false; }
    jsonPath = jp;
    if (const JVal *desc = root.Get("description"))
    {
      const size_t k = desc->str.find("[pO_analysis provenance]");
      jsonProvenance = (k == std::string::npos) ? "(no provenance paragraph -- not the reduced copy?)" : desc->str.substr(k);
    }
    // ---- the trigger study ----
    const char *tp = FindFile(kTrigCandidates);
    if (!tp) { log << "[SF] ERR trig_eff_mb_mu.root not found (tried ../correction/rootfile, correction/rootfile); run correction/run_trig_eff_mb.sh mu\n"; return false; }
    if (!LoadTrig(tp, log)) return false;
    loaded = true;
    Print(log);
    return true;
  }

  // ---- W -> mu nu: the leading muon (TightID, matched); ISO only if it passes the cut
  EventSF W(double pt, double eta, bool passIso) const
  {
    const double fID  = id.Eval(eta, pt, 0, cntID),   fIDu = id.Eval(eta, pt, +1, cntID),   fIDd = id.Eval(eta, pt, -1, cntID);
    const double fIso = passIso ? iso.Eval(eta, pt, 0, cntIso) : 1.0;
    const double fIsou = passIso ? iso.Eval(eta, pt, +1, cntIso) : 1.0;
    const double fIsod = passIso ? iso.Eval(eta, pt, -1, cntIso) : 1.0;
    const int    b    = trig.Bin(-eta, cntTrig.etaClamped); ++cntTrig.calls; // y = -eta_lab (skim convention)
    const double fT = trig.SF(b), fTu = trig.SFUp(b), fTd = trig.SFDn(b);    // inclusive or per-y (kTrigBinning)
    EventSF e;
    e.part[kID] = fID; e.part[kIso] = fIso; e.part[kTrig] = fT;
    e.nom = fID * fIso * fT;
    e.var[kID][0]   = fIDu * fIso * fT;  e.var[kID][1]   = fIDd * fIso * fT;
    e.var[kIso][0]  = fID * fIsou * fT;  e.var[kIso][1]  = fID * fIsod * fT;
    e.var[kTrig][0] = fID * fIso * fTu;  e.var[kTrig][1] = fID * fIso * fTd;
    return e;
  }

  // ---- Z -> mu mu: two TightID legs; the same Tight ISO table as the W per leg,
  // only if that leg passes the (0.15) cut; trigger = per-event path-fired
  // probability from the per-lepton efficiencies (inclusive or per-y, kTrigBinning)
  EventSF Z(double pt1, double eta1, bool iso1, double pt2, double eta2, bool iso2) const
  {
    const double a = id.Eval(eta1, pt1, 0, cntID),  au = id.Eval(eta1, pt1, +1, cntID),  ad = id.Eval(eta1, pt1, -1, cntID);
    const double b = id.Eval(eta2, pt2, 0, cntID),  bu = id.Eval(eta2, pt2, +1, cntID),  bd = id.Eval(eta2, pt2, -1, cntID);
    const double c = iso1 ? iso.Eval(eta1, pt1, 0, cntIso) : 1, cu = iso1 ? iso.Eval(eta1, pt1, +1, cntIso) : 1, cd = iso1 ? iso.Eval(eta1, pt1, -1, cntIso) : 1;
    const double d = iso2 ? iso.Eval(eta2, pt2, 0, cntIso) : 1, du = iso2 ? iso.Eval(eta2, pt2, +1, cntIso) : 1, dd = iso2 ? iso.Eval(eta2, pt2, -1, cntIso) : 1;
    const int    b1 = trig.Bin(-eta1, cntTrig.etaClamped), b2 = trig.Bin(-eta2, cntTrig.etaClamped); ++cntTrig.calls;
    const double eM = 1.0 - (1.0 - trig.EffM(b1)) * (1.0 - trig.EffM(b2));
    const double eD = 1.0 - (1.0 - trig.EffD(b1)) * (1.0 - trig.EffD(b2));
    const double eDu = 1.0 - (1.0 - trig.EffDUp(b1)) * (1.0 - trig.EffDUp(b2));
    const double eDd = 1.0 - (1.0 - trig.EffDDn(b1)) * (1.0 - trig.EffDDn(b2));
    const double fT = eM > 0 ? eD / eM : 1.0, fTu = eM > 0 ? eDu / eM : 1.0, fTd = eM > 0 ? eDd / eM : 1.0;
    const double fID = a * b, fIso = c * d;
    EventSF e;
    e.part[kID] = fID; e.part[kIso] = fIso; e.part[kTrig] = fT;
    e.nom = fID * fIso * fT;
    e.var[kID][0]   = au * bu * fIso * fT; e.var[kID][1]   = ad * bd * fIso * fT;
    e.var[kIso][0]  = fID * cu * du * fT;  e.var[kIso][1]  = fID * cd * dd * fT;
    e.var[kTrig][0] = fID * fIso * fTu;    e.var[kTrig][1] = fID * fIso * fTd;
    return e;
  }

  void PrintCounters(std::ostream &log, const char *who) const
  {
    log << "[SF] " << who << " lookups: ID " << cntID.calls << " (eta clamped " << cntID.etaClamped
        << ", pt clamped " << cntID.ptClamped << "), ISO " << cntIso.calls << " (eta " << cntIso.etaClamped
        << ", pt " << cntIso.ptClamped << "), TRIG " << cntTrig.calls << " (y clamped " << cntTrig.etaClamped << ")\n";
  }

private:
  static const char *FindFile(const char *const *cands)
  {
    for (int i = 0; cands[i]; ++i)
      if (!gSystem->AccessPathName(cands[i], kReadPermission)) return cands[i];
    return nullptr;
  }

  bool LoadTrig(const char *path, std::ostream &log)
  {
    TDirectory::TContext ctx; // restore gDirectory afterwards (the skim books histograms next)
    TFile *f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) { log << "[SF] ERR cannot open " << path << "\n"; return false; }
    const std::string sel = kTrigSel;
    TGraphAsymmErrors *g = (TGraphAsymmErrors *)f->Get(("sf_y_" + sel).c_str());
    TH1 *dDen = (TH1 *)f->Get(("h_den_y_" + sel + "_data_mu").c_str());
    TH1 *dBit = (TH1 *)f->Get(("h_bit_y_" + sel + "_data_mu").c_str());
    TH1 *pDen = (TH1 *)f->Get(("h_den_y_" + sel + "_Wp_mu").c_str());
    TH1 *pBit = (TH1 *)f->Get(("h_bit_y_" + sel + "_Wp_mu").c_str());
    TH1 *mDen = (TH1 *)f->Get(("h_den_y_" + sel + "_Wm_mu").c_str());
    TH1 *mBit = (TH1 *)f->Get(("h_bit_y_" + sel + "_Wm_mu").c_str());
    if (!g || !dDen || !dBit || !pDen || !pBit || !mDen || !mBit)
    {
      log << "[SF] ERR " << path << " lacks sf_y_" << sel << " and/or the h_{den,bit}_y_" << sel
          << "_{data,Wp,Wm}_mu histograms (re-run correction/run_trig_eff_mb.sh mu)\n";
      f->Close(); delete f; return false;
    }
    const int n = g->GetN();
    if (n != dDen->GetNbinsX()) { log << "[SF] ERR sf_y graph and den histogram disagree on the y binning\n"; f->Close(); delete f; return false; }
    trig = TrigTable();
    trig.source = std::string(path) + " : h_num/h_den_y_" + sel + " summed (inclusive W), h_bit/h_den_y_" + sel
                  + " summed (Z per-lepton); sf_y_" + sel + " = the per-y check";
    TH1 *dNum = (TH1 *)f->Get(("h_num_y_" + sel + "_data_mu").c_str());
    TH1 *pNum = (TH1 *)f->Get(("h_num_y_" + sel + "_Wp_mu").c_str());
    TH1 *mNum = (TH1 *)f->Get(("h_num_y_" + sel + "_Wm_mu").c_str());
    if (!dNum || !pNum || !mNum)
    { log << "[SF] ERR " << path << " lacks the h_num_y_" << sel << "_{data,Wp,Wm}_mu histograms\n"; f->Close(); delete f; return false; }
    double sDd = 0, sDn = 0, sMd = 0, sMn = 0, sDb = 0, sMb = 0; // summed den / num(matched) / bit(fired), data and MC
    for (int i = 0; i < n; ++i)
    {
      sDd += dDen->GetBinContent(i + 1); sDn += dNum->GetBinContent(i + 1); sDb += dBit->GetBinContent(i + 1);
      sMd += pDen->GetBinContent(i + 1) + mDen->GetBinContent(i + 1);
      sMn += pNum->GetBinContent(i + 1) + mNum->GetBinContent(i + 1);
      sMb += pBit->GetBinContent(i + 1) + mBit->GetBinContent(i + 1);
    }
    if (sDd <= 0 || sMd <= 0 || sMn <= 0) { log << "[SF] ERR empty trigger denominators in " << path << "\n"; f->Close(); delete f; return false; }
    {
      const double eD = sDn / sDd, eM = sMn / sMd;
      trig.sfI   = eD / eM;
      trig.sfIUp = TEfficiency::ClopperPearson((int)sDd, (int)sDn, 0.683, true)  / eM;
      trig.sfIDn = TEfficiency::ClopperPearson((int)sDd, (int)sDn, 0.683, false) / eM;
      trig.effDI   = sDb / sDd;
      trig.effDIUp = TEfficiency::ClopperPearson((int)sDd, (int)sDb, 0.683, true);
      trig.effDIDn = TEfficiency::ClopperPearson((int)sDd, (int)sDb, 0.683, false);
      trig.effMI   = sMb / sMd;
    }
    for (int i = 0; i < n; ++i)
    {
      const double x = g->GetX()[i];
      trig.yEdges.push_back(x - g->GetEXlow()[i]);
      if (i == n - 1) trig.yEdges.push_back(x + g->GetEXhigh()[i]);
      trig.sf.push_back(g->GetY()[i]);
      trig.sfUp.push_back(g->GetY()[i] + g->GetEYhigh()[i]);
      trig.sfDn.push_back(g->GetY()[i] - g->GetEYlow()[i]);
      // per-lepton path-fired efficiencies (Clopper-Pearson 68% on the data count)
      const double nd = dDen->GetBinContent(i + 1), kd = dBit->GetBinContent(i + 1);
      const double nm = pDen->GetBinContent(i + 1) + mDen->GetBinContent(i + 1);
      const double km = pBit->GetBinContent(i + 1) + mBit->GetBinContent(i + 1);
      const double ed = nd > 0 ? kd / nd : 1.0, em = nm > 0 ? km / nm : 1.0;
      trig.effD.push_back(ed);
      trig.effM.push_back(em);
      trig.effDUp.push_back(nd > 0 ? TEfficiency::ClopperPearson((int)nd, (int)kd, 0.683, true)  : 1.0);
      trig.effDDn.push_back(nd > 0 ? TEfficiency::ClopperPearson((int)nd, (int)kd, 0.683, false) : 1.0);
    }
    // bin-edge consistency between the graph (x +- ex) and the histogram axis
    for (int i = 0; i <= n; ++i)
    {
      const double hEdge = (i < n) ? dDen->GetXaxis()->GetBinLowEdge(i + 1) : dDen->GetXaxis()->GetBinUpEdge(n);
      if (std::fabs(trig.yEdges[i] - hEdge) > 1e-6)
      { log << "[SF] ERR trigger y edges inconsistent between graph and histogram at edge " << i << "\n"; f->Close(); delete f; return false; }
    }
    // the standing check: are the per-y SFs consistent with one flat (inclusive) SF?
    // (pull per bin with the asymmetric error on the side facing the inclusive value)
    trig.chi2 = 0.0; trig.ndf = n - 1;
    for (int i = 0; i < n; ++i)
    {
      const double s = (trig.sfI < trig.sf[i]) ? (trig.sf[i] - trig.sfDn[i]) : (trig.sfUp[i] - trig.sf[i]);
      if (s > 0) trig.chi2 += std::pow((trig.sf[i] - trig.sfI) / s, 2);
    }
    f->Close(); delete f;
    return true;
  }

  void Print(std::ostream &log) const
  {
    const auto rID = id.Range(25.0), rIso = iso.Range(25.0), rIsoZ = iso.Range(10.0);
    log << "[SF] muon scale factors ON (MC event weight x ID x ISO x TRIG; PO_MUON_SF=off disables)\n"
        << "[SF]   JSON  : " << jsonPath << "\n"
        << "[SF]           " << jsonProvenance << "\n"
        << "[SF]   ID    : " << id.name << "  (" << id.NEta() << " signed-eta x " << id.NPt()
        << " pT bins; nominal range at pT>25: " << Form("%.4f..%.4f", rID.first, rID.second) << ")\n"
        << "[SF]   ISO   : " << iso.name << "  (relIso < 0.15 in W and Z; nominal range at pT>25: "
        << Form("%.4f..%.4f", rIso.first, rIso.second) << ", at pT>10 (Z legs): " << Form("%.4f..%.4f", rIsoZ.first, rIsoZ.second) << ")\n"
        << "[SF]   TRIG  : " << trig.source << "\n"
        << "[SF]           APPLIED " << (kTrigBinning == kTrigInclusive ? "INCLUSIVE" : "PER y BIN")
        << " (kTrigBinning); muTrig nuisance correlation for the fork: " << kMuTrigCorr << "\n"
        << "[SF]           inclusive SF(W: fired&&matched) = " << Form("%.4f (+%.4f -%.4f)", trig.sfI, trig.sfIUp - trig.sfI, trig.sfI - trig.sfIDn)
        << ";  per-lepton eps(fired) data " << Form("%.4f (+%.4f -%.4f)", trig.effDI, trig.effDIUp - trig.effDI, trig.effDI - trig.effDIDn)
        << ", MC " << Form("%.4f", trig.effMI) << "\n"
        << "[SF]           per-y check (" << (kTrigBinning == kTrigInclusive ? "NOT applied" : "applied") << "):  y bin        SF(W)                    eps_data(fired)          eps_MC   pull\n";
    for (int i = 0; i < trig.NBins(); ++i)
    {
      const double s = (trig.sfI < trig.sf[i]) ? (trig.sf[i] - trig.sfDn[i]) : (trig.sfUp[i] - trig.sf[i]);
      log << "[SF]                                       " << Form("[%5.2f,%5.2f]  %.4f (+%.4f -%.4f)   %.4f (+%.4f -%.4f)   %.4f   %+.2f",
                                       trig.yEdges[i], trig.yEdges[i + 1], trig.sf[i], trig.sfUp[i] - trig.sf[i],
                                       trig.sf[i] - trig.sfDn[i], trig.effD[i], trig.effDUp[i] - trig.effD[i],
                                       trig.effD[i] - trig.effDDn[i], trig.effM[i], s > 0 ? (trig.sf[i] - trig.sfI) / s : 0.0) << "\n";
    }
    log << "[SF]           per-y SFs vs the inclusive one: chi2/ndf = " << Form("%.1f / %d (p = %.2f)", trig.chi2, trig.ndf, TMath::Prob(trig.chi2, trig.ndf))
        << " -> " << (TMath::Prob(trig.chi2, trig.ndf) > 0.05 ? "consistent with one flat SF" : "TENSION with a flat SF (revisit kTrigBinning)") << "\n"
        << "[SF]   per-source twins <hist>_" << kMuonSFSourceNames[0] << "/" << kMuonSFSourceNames[1] << "/"
        << kMuonSFSourceNames[2] << "Up|Down (one source at +-1 sigma, others nominal) + the COMBINED <hist>_"
        << kMuonSFCombinedName << "Up|Down (per-bin quadrature of the three) = the one nuisance carried into the fit\n";
  }
};

// ============================================================================
// Twins of a fit-template histogram: the per-source <h>_<src>Up/Down (filled
// per event) and the combined <h>_muSFUp/Down (built once at the end)
// ============================================================================
struct SFTwins
{
  const TH1 *nom = nullptr;                                           // the nominal (filled with w_gen x SF)
  TH1D *h[kNSources][2] = {{nullptr, nullptr}, {nullptr, nullptr}, {nullptr, nullptr}}; // per source, Up / Down
  TH1D *comb[2] = {nullptr, nullptr};                                 // combined Up / Down
};

inline TH1D *CloneEmptyTwin(const TH1 *nom, const std::string &name, const std::string &note)
{
  const std::string title = nom->GetTitle();
  const size_t      semi  = title.find(';');
  TH1D *h = (TH1D *)nom->Clone(name.c_str());
  h->Reset("ICES");
  h->SetTitle(semi == std::string::npos ? (title + note).c_str() : (title.substr(0, semi) + note + title.substr(semi)).c_str());
  if (h->GetSumw2N() == 0) h->Sumw2();
  return h;
}

// Books 3 x 2 per-source twins + the 2 combined ones; the set is appended to
// `registry` for FinalizeAndWriteSFTwins() at the end of the job.
inline SFTwins BookSFTwins(const TH1 *nom, std::vector<SFTwins> &registry)
{
  SFTwins t;
  t.nom = nom;
  for (int s = 0; s < kNSources; ++s)
    for (int d = 0; d < 2; ++d)
      t.h[s][d] = CloneEmptyTwin(nom, std::string(nom->GetName()) + "_" + kMuonSFSourceNames[s] + (d == 0 ? "Up" : "Down"),
                                 std::string(" [") + kMuonSFSourceNames[s] + (d == 0 ? " +1 sigma" : " -1 sigma") + ", others nominal]");
  for (int d = 0; d < 2; ++d)
    t.comb[d] = CloneEmptyTwin(nom, std::string(nom->GetName()) + "_" + kMuonSFCombinedName + (d == 0 ? "Up" : "Down"),
                               std::string(" [") + kMuonSFCombinedName + (d == 0 ? " +1 sigma" : " -1 sigma")
                                   + ": ID, ISO, trigger shifts added in quadrature per bin]");
  registry.push_back(t);
  return t;
}

// wBase = the event weight WITHOUT the SF (gen weight [x LHE member 0]); the twin
// receives wBase x (total factor with source s varied)
inline void FillSFTwins(const SFTwins &t, double x, double wBase, const EventSF &e)
{
  for (int s = 0; s < kNSources; ++s)
    for (int d = 0; d < 2; ++d)
      if (t.h[s][d]) t.h[s][d]->Fill(x, wBase * e.var[s][d]);
}

// The combined twins: per bin (under/overflow included)
//   Up   = nominal + sqrt( sum_s (Up_s   - nominal)^2 )
//   Down = nominal - sqrt( sum_s (nominal - Down_s)^2 )   (floored at 0 for non-negative bins)
// = the 1 sigma of the product of three independent factors, one coherent
// parameter. Bin errors = the nominal's (Combine reads the contents only).
inline void FinalizeSFTwins(SFTwins &t)
{
  if (!t.nom || !t.comb[0] || !t.comb[1]) return;
  const int nb = t.nom->GetNbinsX();
  for (int b = 0; b <= nb + 1; ++b)
  {
    const double n = t.nom->GetBinContent(b);
    double up2 = 0.0, dn2 = 0.0;
    for (int s = 0; s < kNSources; ++s)
    {
      const double du = t.h[s][0]->GetBinContent(b) - n;
      const double dd = n - t.h[s][1]->GetBinContent(b);
      up2 += du * du;
      dn2 += dd * dd;
    }
    const double dn = n - std::sqrt(dn2);
    t.comb[0]->SetBinContent(b, n + std::sqrt(up2));
    t.comb[1]->SetBinContent(b, (n >= 0.0 && dn < 0.0) ? 0.0 : dn);
    t.comb[0]->SetBinError(b, t.nom->GetBinError(b));
    t.comb[1]->SetBinError(b, t.nom->GetBinError(b));
  }
  t.comb[0]->SetEntries(t.nom->GetEntries());
  t.comb[1]->SetEntries(t.nom->GetEntries());
}

// Finalize every set and write all 8 histograms per set into the current
// directory (the output file). Returns the number of histograms written.
inline int FinalizeAndWriteSFTwins(std::vector<SFTwins> &sets)
{
  int n = 0;
  for (SFTwins &t : sets)
  {
    FinalizeSFTwins(t);
    for (int s = 0; s < kNSources; ++s)
      for (int d = 0; d < 2; ++d)
        if (t.h[s][d]) { t.h[s][d]->Write("", TObject::kOverwrite); ++n; }
    for (int d = 0; d < 2; ++d)
      if (t.comb[d]) { t.comb[d]->Write("", TObject::kOverwrite); ++n; }
  }
  return n;
}

// ============================================================================
// Running <SF> bookkeeping for the job log
// ============================================================================
struct SFStats
{
  double n = 0, sw = 0, swsf = 0, part[kNSources] = {0, 0, 0}, varUp[kNSources] = {0, 0, 0}, varDn[kNSources] = {0, 0, 0};
  void Add(double wBase, const EventSF &e)
  {
    n += 1; sw += wBase; swsf += wBase * e.nom;
    for (int s = 0; s < kNSources; ++s)
    { part[s] += wBase * e.part[s]; varUp[s] += wBase * e.var[s][0]; varDn[s] += wBase * e.var[s][1]; }
  }
  void Print(std::ostream &log, const char *what) const
  {
    if (sw <= 0) { log << "[SF] " << what << ": no weighted events\n"; return; }
    log << "[SF] " << what << ": N = " << (long long)n << ", sum w = " << Form("%.1f", sw)
        << ", sum w*SF = " << Form("%.1f", swsf) << "  ->  <SF> = " << Form("%.5f", swsf / sw)
        << "  (<ID> " << Form("%.5f", part[kID] / sw) << ", <ISO> " << Form("%.5f", part[kIso] / sw)
        << ", <TRIG> " << Form("%.5f", part[kTrig] / sw) << ")\n";
    double up2 = 0.0, dn2 = 0.0;
    for (int s = 0; s < kNSources; ++s)
    {
      log << "[SF]     " << Form("%-7s", kMuonSFSourceNames[s]) << " Up/Down on sum w*SF: "
          << Form("%+.3f%% / %+.3f%%", 100.0 * (varUp[s] / swsf - 1.0), 100.0 * (varDn[s] / swsf - 1.0)) << "\n";
      up2 += std::pow(varUp[s] / swsf - 1.0, 2);
      dn2 += std::pow(varDn[s] / swsf - 1.0, 2);
    }
    log << "[SF]     " << Form("%-7s", kMuonSFCombinedName) << " Up/Down (quadrature of the three): "
        << Form("%+.3f%% / %+.3f%%", 100.0 * std::sqrt(up2), -100.0 * std::sqrt(dn2)) << "  <- the one nuisance in the fit\n";
  }
};

} // namespace pOSF

#endif // PO_MUON_SF_H
