// correction/qcd_nomt_study.C
//
// STUDY (2026-10-21, user decisions + questions): the multijet background of
// the lepton-pT fit WITHOUT the m_T > 40 cut, which becomes the nominal (the
// cut is dropped because of the known electron PF-MET issue; leppt_mt40 stays
// for future use). The question is how the in-fit ABCD should look when the
// signal region spans all m_T, the user's proposal being
//   SR multijet(pT, y) = F_g(pT) x C_y(pT)
// with C_y(pT) = the anti-iso pT spectrum of EACH rapidity bin at all m_T
// ("per rapidity bin shape") and F_g(pT) = the pT-dependent transfer factor
// QCD(iso-pass)/QCD(anti-iso), measured at m_T < 30 per rapidity group
// (fwd y10-11 / bwd y0-1 / cen y2-9, abcd_ygroups.h).
//
// What it measures, per flavour and charge (everything from the y-binned
// (iso, pT, m_T) cubes of the W skims, written since 2026-10-21):
//  [COMP]  the composition of the no-m_T SR per rapidity bin (data, W MC,
//          other MC, the F x C multijet estimate) next to the m_T > 40 SR;
//  [FPT]   F(pT) per group and inclusive at m_T < 30: profile-likelihood
//          intervals per bin, the exponential slope s of F = F0 exp(s(pT-40)/10)
//          with its 68% interval, the flatness LRTs, and whether the level F0
//          and the slope s depend on the group (nested LRTs);
//          + the anti-iso window variation ([0.3,0.6) / [0.6,1.0) of the
//          relIso axis) of F0 and s;
//  [TY]    the pT-integrated T_y = QCD(iso-pass)/QCD(anti-iso) at m_T < 30
//          per rapidity bin;
//  [SHAPE] is the anti-iso pT shape of a rapidity bin compatible with the rest
//          of the acceptance (chi2 homogeneity test, bin vs complement) and
//          with its group;
//  [MTFAC] factorization: F (pT-integrated) in m_T slices 0-10/10-20/20-30/
//          30-40 relative to m_T < 30, at the prefit W scale and at r_W = 1.2,
//          and the fraction of the SR multijet sitting at m_T >= 30 (the part
//          that relies on carrying F across m_T);
//  [FISHER] expected precision (Asimov, prefit) of three in-fit schemes:
//            A  SR = iso-pass, all m_T; per-y anti-iso count CRs (anchors);
//               the level T_g of each group FREE -- fixed by the SR pT shape
//               alone: no m_T anywhere in the likelihood
//            C  the same SR split at m_T = 30 into two channels (shared r):
//               the low-m_T part plays region B and pins T_g, but the W MC
//               m_T spectrum then decides how W splits between the two
//            R  reference = today's leppt_mt40 in-fit ABCD per group (SR at
//               m_T > 40, 1-bin B / C40 / D per group) -- a different
//               fiducial volume, quoted for the precision scale only
//          (slope s Gaussian-constrained to its m_T < 30 measurement);
//  [FIT]   the data stand-in fits of the same schemes (Poisson likelihood,
//          Minuit2): r per rapidity bin, the fitted T_g against the m_T < 30
//          measurement F0_g (THE closure test of scheme A), the postfit
//          Baker-Cousins chi2; scheme A also with the group / inclusive
//          anti-iso shape instead of the per-bin one, and with a flat F;
//  [EMU]   the e/mu ratios of the fitted r per rapidity bin.
// Every m_T-dependent item is done twice (user question 2026-10-21): with m_T
// from the PF MET (cube h_iso_pt_mt_*) and from the CHARGED-only MET (cube
// h_iso_pt_mttrk_*, skim_common.h ComputeChargedPFMET).
//
// The stand-in fits are a cross-check of the design, NOT the measurement: no
// shape systematics, no lumi / SF nuisances, the Z peak not included (DY
// fixed at its MC), r per bin per flavour (no mu/e sharing). The electron MC
// carries NO Reco SF since 2026-10-21 (eleReco = +-10% in the real fit), so
// the electron r here are ~7% below the pre-2026-10-21 ones by construction.
//
// Run through the wrapper:  ./run_qcd_nomt_study.sh  (log = the record)
// Output: plots/qcd_nomt_study/  (figures + qcd_nomt_study.csv)

#include "qcd_abcd.C" // MCFile, ffLnL / ffInterval / ffFlatness / ffSlopeInterval / ffMaximize,
                      // plotting_helper.C, mc_norm.h, skim_common.h, abcd_ygroups.h

#include "TH3D.h"
#include "TMatrixDSym.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

#include <fstream>
#include <map>
#include <memory>

namespace nomt
{
const int kNY = pOSkim::kNY;
const int kNG = pOAbcdY::kNGroups;
const double kMtB   = 30.0;  // B / D: m_T < 30 (= qcd_abcd.C yCut)
const double kMtSR  = 40.0;  // the leppt_mt40 cut (reference scheme R)
const double kBig   = 1e9;   // "to the overflow"
const double kRAlt  = 1.2;   // the alternative W scale of the subtractions (the fitted r ~ 1.2)
const int    kFitLo = 13, kFitHi = 50; // the fit window in the 2-GeV pT histograms: [24, 100), overflow dropped
// F(pT) bins (2-GeV-grid edges; [24,28) holds 25-28 GeV, the last one includes the > 100 overflow)
const std::vector<double> kFE = {24, 28, 32, 36, 40, 46, 52, 60, 70, 80, 100, 1e9};
const std::vector<double> kMtSl = {0, 10, 20, 30, 40}; // the m_T slices of the factorization test

std::ofstream gCsv;
int gUid = 0;
std::string Uid(const std::string &s) { return s + Form("_%d", ++gUid); }

struct LepCfg
{
  std::string lep, fileBase, sym;
  bool isElec = false;
  double pLo = 0, pHi = 0, aLo = 0, aMid = 0.6, aHi = 1.0; // iso-pass / anti-iso windows on the cube's iso axis
};

LepCfg MakeCfg(bool isElec)
{
  LepCfg c;
  c.isElec = isElec;
  c.lep = isElec ? "ele" : "mu";
  c.sym = isElec ? "e" : "#mu";
  c.fileBase = isElec ? "WToElecNu_pO_PFMet" : "WToMuNu_pO_PFMet";
  if (isElec) { c.pLo = pOSkim::kEleIsoAxisLo; c.pHi = 0.0; c.aLo = pOSkim::kEleAntiIsoLo; c.aHi = pOSkim::kEleAntiIsoHi; }
  else        { c.pLo = 0.0; c.pHi = 0.15; c.aLo = 0.30; c.aHi = 1.0; } // = qcd_abcd.C ABCDConfig (muon)
  return c;
}

struct Files
{
  TFile *data = nullptr;
  std::vector<MCFile> w, o; // W-related (W, W->tau nu) and every other EWK sample (DY, DY->tautau, ttbar)
};

bool OpenFiles(const LepCfg &c, Files &F)
{
  const std::string dir = "../skim/rootfile/";
  F.data = TFile::Open((dir + c.fileBase + "_hist.root").c_str(), "READ");
  if (!F.data || F.data->IsZombie()) { std::cerr << "[ERROR] cannot open the " << c.lep << " data skim\n"; return false; }
  const std::string dyLab = c.isElec ? "DYee" : "DYmu";
  const std::vector<std::pair<std::string, std::string>> spec = {
      {"Wp_" + c.lep, "Wp"}, {"Wm_" + c.lep, "Wm"}, {"Wp_tau", "Wptau"}, {"Wm_tau", "Wmtau"},
      {dyLab, "DY"}, {"DYtau", "DYtau"}, {"TT", "TT"}};
  for (const auto &s : spec)
  {
    TFile *f = TFile::Open((dir + c.fileBase + "_" + s.second + "_hist.root").c_str(), "READ");
    if (!f || f->IsZombie()) { std::cerr << "[ERROR] cannot open the " << s.second << " skim\n"; return false; }
    const bool isW = (s.second.rfind("Wp", 0) == 0 || s.second.rfind("Wm", 0) == 0);
    (isW ? F.w : F.o).push_back({s.first, f});
  }
  return true;
}

// the three cubes of one rapidity bin: data, W-related MC x k_s, other MC x k_s
struct Cube { TH3D *d = nullptr, *w = nullptr, *o = nullptr; };

TH3D *SumCube(const std::vector<MCFile> &v, const std::string &name, bool &ok)
{
  TH3D *sum = nullptr;
  for (const MCFile &m : v)
  {
    TH3D *h = (TH3D *)m.f->Get(name.c_str());
    if (!h) { std::cerr << "[ERROR] " << name << " missing in " << m.label << " (skim before 2026-10-21?)\n"; ok = false; return nullptr; }
    if (!sum) { sum = (TH3D *)h->Clone(Uid(name + "_sum").c_str()); sum->SetDirectory(nullptr); sum->Reset("ICESM"); }
    sum->Add(h, pONorm::MCScale(m.label));
  }
  return sum;
}

// a pT spectrum of each component, and the QCD-only one at a W scale rW
struct Comp
{
  TH1D *d = nullptr, *w = nullptr, *o = nullptr;
  TH1D *Qcd(double rW, const char *nm) const
  {
    TH1D *q = (TH1D *)d->Clone(Uid(nm).c_str());
    q->SetDirectory(nullptr);
    q->Add(w, -rW);
    q->Add(o, -1.0);
    return q;
  }
};

// cube bin whose low (or high) edge is exactly v -- every window edge must be a bin edge
int EdgeBin(const TAxis *a, double v, bool lowEdge)
{
  for (int b = 1; b <= a->GetNbins(); ++b)
    if (std::fabs((lowEdge ? a->GetBinLowEdge(b) : a->GetBinUpEdge(b)) - v) < 1e-9) return b;
  return -1;
}

TH1D *Slice1(const TH3D *h, double ilo, double ihi, double mlo, double mhi, bool &ok)
{
  const TAxis *ax = h->GetXaxis(), *az = h->GetZaxis();
  const int bx1 = EdgeBin(ax, ilo, true), bx2 = EdgeBin(ax, ihi, false);
  const int bz1 = EdgeBin(az, mlo, true);
  const int bz2 = (mhi >= az->GetXmax()) ? az->GetNbins() + 1 : EdgeBin(az, mhi, false);
  if (bx1 < 0 || bx2 < 0 || bz1 < 0 || bz2 < 0)
  {
    std::cerr << "[ERROR] window iso [" << ilo << "," << ihi << ") x m_T [" << mlo << "," << mhi
              << ") is not on the cube's bin edges (" << h->GetName() << ")\n";
    ok = false;
    return nullptr;
  }
  TH1D *p = ((TH3D *)h)->ProjectionY(Uid("sl").c_str(), bx1, bx2, bz1, bz2, "e");
  p->SetDirectory(nullptr);
  return p;
}

Comp Slice(const Cube &c, double ilo, double ihi, double mlo, double mhi, bool &ok)
{
  Comp r;
  r.d = Slice1(c.d, ilo, ihi, mlo, mhi, ok);
  r.w = Slice1(c.w, ilo, ihi, mlo, mhi, ok);
  r.o = Slice1(c.o, ilo, ihi, mlo, mhi, ok);
  return r;
}

// sum (and variance) of a pT histogram over [lo, hi); hi >= 1e8 includes the overflow
double Cnt(const TH1D *h, double lo, double hi, double *var = nullptr)
{
  const double eps = 1e-6;
  const int b1 = h->GetXaxis()->FindBin(lo + eps);
  const int b2 = (hi >= 1e8) ? h->GetNbinsX() + 1 : h->GetXaxis()->FindBin(hi - eps);
  double e = 0.0;
  const double s = h->IntegralAndError(b1, b2, e);
  if (var) *var = e * e;
  return s;
}

// the regions of one rapidity bin
struct YReg
{
  Comp P_lo, P_hi, P_40, P_all, A_lo, A_hi, A_40, A_all, A1_lo, A2_lo;
  Comp P_m[4], A_m[4]; // the m_T slices kMtSl
};

YReg MakeYReg(const Cube &c, const LepCfg &L, bool &ok)
{
  YReg r;
  r.P_lo  = Slice(c, L.pLo, L.pHi, 0, kMtB, ok);
  r.P_hi  = Slice(c, L.pLo, L.pHi, kMtB, kBig, ok);
  r.P_40  = Slice(c, L.pLo, L.pHi, kMtSR, kBig, ok);
  r.P_all = Slice(c, L.pLo, L.pHi, 0, kBig, ok);
  r.A_lo  = Slice(c, L.aLo, L.aHi, 0, kMtB, ok);
  r.A_hi  = Slice(c, L.aLo, L.aHi, kMtB, kBig, ok);
  r.A_40  = Slice(c, L.aLo, L.aHi, kMtSR, kBig, ok);
  r.A_all = Slice(c, L.aLo, L.aHi, 0, kBig, ok);
  r.A1_lo = Slice(c, L.aLo, L.aMid, 0, kMtB, ok);
  r.A2_lo = Slice(c, L.aMid, L.aHi, 0, kMtB, ok);
  for (int m = 0; m < 4; ++m)
  {
    r.P_m[m] = Slice(c, L.pLo, L.pHi, kMtSl[m], kMtSl[m + 1], ok);
    r.A_m[m] = Slice(c, L.aLo, L.aHi, kMtSl[m], kMtSl[m + 1], ok);
  }
  return r;
}

// add the components of several rapidity bins (a group or the acceptance)
Comp SumComp(const std::vector<const Comp *> &v)
{
  Comp s;
  for (const Comp *c : v)
  {
    if (!s.d)
    {
      s.d = (TH1D *)c->d->Clone(Uid("sd").c_str()); s.d->SetDirectory(nullptr);
      s.w = (TH1D *)c->w->Clone(Uid("sw").c_str()); s.w->SetDirectory(nullptr);
      s.o = (TH1D *)c->o->Clone(Uid("so").c_str()); s.o->SetDirectory(nullptr);
    }
    else { s.d->Add(c->d); s.w->Add(c->w); s.o->Add(c->o); }
  }
  return s;
}

std::vector<int> GroupBins(int g) { return g < 0 ? std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11} : pOAbcdY::BinsOf(g); }
std::string GroupName(int g) { return g < 0 ? "incl" : pOAbcdY::kGroupName[g]; }

// ============================================================================
// F(pT) at low m_T: the qcd_abcd.C fake-factor likelihood (per bin n_p ~
// Pois(F mu + b_p), n_f ~ Pois(mu + b_f), mu profiled, the EWK b known)
// ============================================================================
struct FMeas
{
  std::vector<double> x, np, nf, bp, bf;
  FFFlatness fl;
  double sLo = 0, sHi = 0;
  bool sOk = false;
  std::vector<FFInterval> bin;
  FFInterval Tint; // pT-integrated
};

FMeas MeasureF(const Comp &pass, const Comp &anti, double rW)
{
  FMeas m;
  for (size_t k = 0; k + 1 < kFE.size(); ++k)
  {
    const double lo = kFE[k], hi = kFE[k + 1];
    const double np = Cnt(pass.d, lo, hi), nf = Cnt(anti.d, lo, hi);
    const double bp = std::max(rW * Cnt(pass.w, lo, hi) + Cnt(pass.o, lo, hi), 1e-6);
    const double bf = std::max(rW * Cnt(anti.w, lo, hi) + Cnt(anti.o, lo, hi), 1e-6);
    // x = the anti-iso data mean pT in the bin (the first 2-GeV bin holds 25-26 GeV; > 100 at 105)
    double sx = 0, sw = 0;
    if (hi < 1e8)
      for (int b = anti.d->FindBin(lo + 1e-6); b <= anti.d->FindBin(hi - 1e-6); ++b)
      {
        const double xc = (std::fabs(anti.d->GetBinLowEdge(b) - 24.0) < 1e-6) ? 25.5 : anti.d->GetBinCenter(b);
        sx += xc * anti.d->GetBinContent(b);
        sw += anti.d->GetBinContent(b);
      }
    const double x = (hi >= 1e8) ? 105.0 : (sw > 0 ? sx / sw : 0.5 * (lo + hi));
    m.x.push_back(x); m.np.push_back(np); m.nf.push_back(nf); m.bp.push_back(bp); m.bf.push_back(bf);
    m.bin.push_back(ffInterval(np, nf, bp, bf));
  }
  m.fl = ffFlatness(m.x, m.np, m.nf, m.bp, m.bf, true);
  m.sOk = ffSlopeInterval(m.x, m.np, m.nf, m.bp, m.bf, m.fl.s, m.sLo, m.sHi);
  double Np = 0, Nf = 0, Bp = 0, Bf = 0;
  for (size_t k = 0; k < m.x.size(); ++k) { Np += m.np[k]; Nf += m.nf[k]; Bp += m.bp[k]; Bf += m.bf[k]; }
  m.Tint = ffInterval(Np, Nf, Bp, Bf);
  return m;
}

double LnGroupSum(const FMeas &m, double F0, double s)
{
  double t = 0;
  for (size_t k = 0; k < m.x.size(); ++k) t += ffLnL(F0 * std::exp(s * (m.x[k] - 40.0) / 10.0), m.np[k], m.nf[k], m.bp[k], m.bf[k]);
  return t;
}

// the nested group tests: common (F0, s) | per-group F0, common s | per-group (F0, s)
struct GroupTest
{
  double lnCommon = 0, lnF0g = 0, lnFull = 0, sShared = 0, F0shared[kNG] = {0, 0, 0};
  double F0lo[kNG] = {0, 0, 0}, F0hi[kNG] = {0, 0, 0}; // 68% profile interval of F0_g at sShared
};

GroupTest TestGroups(const FMeas g[kNG])
{
  GroupTest t;
  for (int i = 0; i < kNG; ++i) t.lnFull += g[i].fl.lnSlope;
  auto inner = [&](double s, double *F0out) {
    double tot = 0;
    for (int i = 0; i < kNG; ++i)
    {
      double best = 0;
      const double F0 = ffMaximize([&](double f) { return LnGroupSum(g[i], f, s); }, 0.0, 3.0, 60, best);
      if (F0out) F0out[i] = F0;
      tot += best;
    }
    return tot;
  };
  t.sShared = ffMaximize([&](double s) { return inner(s, nullptr); }, -1.5, 1.5, 40, t.lnF0g);
  inner(t.sShared, t.F0shared);
  auto common = [&](double s) {
    double best = 0;
    ffMaximize([&](double f) { double v = 0; for (int i = 0; i < kNG; ++i) v += LnGroupSum(g[i], f, s); return v; },
               0.0, 3.0, 60, best);
    return best;
  };
  ffMaximize(common, -1.5, 1.5, 40, t.lnCommon);
  t.lnF0g = std::max(t.lnF0g, t.lnCommon);
  t.lnFull = std::max(t.lnFull, t.lnF0g);
  for (int i = 0; i < kNG; ++i) // F0_g interval at the shared slope (s fixed), 2 DlnL = 1
  {
    const double Lmax = LnGroupSum(g[i], t.F0shared[i], t.sShared);
    auto q = [&](double f) { return 2.0 * (Lmax - LnGroupSum(g[i], f, t.sShared)) - 1.0; };
    for (int side = -1; side <= 1; side += 2)
    {
      double in = t.F0shared[i], out = in;
      for (int k = 1; k <= 3000; ++k) { out = t.F0shared[i] + side * 0.001 * k; if (out <= 0 || q(out) > 0) break; in = out; }
      for (int it = 0; it < 40 && out > 0; ++it) { const double mm = 0.5 * (in + out); (q(mm) > 0 ? out : in) = mm; }
      (side < 0 ? t.F0lo[i] : t.F0hi[i]) = std::max(0.0, 0.5 * (in + out));
    }
  }
  return t;
}

// ============================================================================
// The stand-in likelihood (schemes A / C / R): Poisson SR pT bins + 1-bin
// count CRs, one Gaussian-constrained slope s, everything else free.
// ============================================================================
struct SRChan
{
  int iy = 0;
  int part = 0; // scheme C: 1 = the m_T < 30 channel, 2 = m_T >= 30; 0 otherwise
  std::vector<int> mPar;     // the QCD multiplier = prod p[mPar_i]^mExp_i x mConst
  std::vector<double> mExp;
  double mConst = 1.0;
  std::vector<double> n, W, O, Q, x; // per pT bin of the fit window
};
struct CRChan
{
  double n = 0, base = 0, ewk = 0; // nu = p[par] * base + ewk + sum_i p[rPar_i] * rW_i
  int par = -1;
  std::vector<int> rPar;
  std::vector<double> rW;
};
struct Model
{
  std::string name;
  std::vector<std::string> pname;
  std::vector<double> p0, plo, phi;
  std::vector<SRChan> sr;
  std::vector<CRChan> cr;
  int sPar = -1;
  double sMean = 0, sSig = 1;
  bool fixS = false;
  int sPar2 = -1; // scheme C_2s: the slope of the m_T >= 30 channels (same prior as sPar)
  double Slope(const SRChan &c, const double *p) const { return (c.part == 2 && sPar2 >= 0) ? p[sPar2] : p[sPar]; }
  int AddPar(const std::string &nm, double v, double lo, double hi)
  {
    pname.push_back(nm); p0.push_back(v); plo.push_back(lo); phi.push_back(hi);
    return (int)pname.size() - 1;
  }
  double Mult(const SRChan &c, const double *p) const
  {
    double m = c.mConst;
    for (size_t i = 0; i < c.mPar.size(); ++i) m *= std::pow(p[c.mPar[i]], c.mExp[i]);
    return m;
  }
  // every expected count, SR bins first then CRs (the order of Obs())
  void Nu(const double *p, std::vector<double> &nu) const
  {
    nu.clear();
    for (const SRChan &c : sr)
    {
      const double m = Mult(c, p), r = p[c.iy];
      for (size_t k = 0; k < c.n.size(); ++k)
        nu.push_back(r * c.W[k] + m * std::exp(Slope(c, p) * (c.x[k] - 40.0) / 10.0) * c.Q[k] + c.O[k]);
    }
    for (const CRChan &c : cr)
    {
      double v = p[c.par] * c.base + c.ewk;
      for (size_t i = 0; i < c.rPar.size(); ++i) v += p[c.rPar[i]] * c.rW[i];
      nu.push_back(v);
    }
  }
  void Obs(std::vector<double> &n) const
  {
    n.clear();
    for (const SRChan &c : sr) n.insert(n.end(), c.n.begin(), c.n.end());
    for (const CRChan &c : cr) n.push_back(c.n);
  }
  double NLL(const double *p, const std::vector<double> &n) const
  {
    std::vector<double> nu;
    Nu(p, nu);
    double t = 0;
    for (size_t i = 0; i < nu.size(); ++i)
    {
      const double v = std::max(nu[i], 1e-9);
      t += v - (n[i] > 0 ? n[i] * std::log(v) : 0.0);
    }
    if (!fixS) t += 0.5 * std::pow((p[sPar] - sMean) / sSig, 2);
    if (sPar2 >= 0) t += 0.5 * std::pow((p[sPar2] - sMean) / sSig, 2);
    return t;
  }
};

// Asimov (n = nu(p0)) covariance from the Fisher information
bool FisherCov(const Model &M, TMatrixDSym &cov)
{
  const int np = (int)M.p0.size();
  std::vector<double> nu0, nuP, nuM, p = M.p0;
  M.Nu(p.data(), nu0);
  const int nb = (int)nu0.size();
  std::vector<std::vector<double>> J(np, std::vector<double>(nb, 0.0));
  for (int a = 0; a < np; ++a)
  {
    if (a == M.sPar && M.fixS) continue;
    const double h = 1e-4 * std::max(std::fabs(M.p0[a]), 0.1);
    p[a] = M.p0[a] + h; M.Nu(p.data(), nuP);
    p[a] = M.p0[a] - h; M.Nu(p.data(), nuM);
    p[a] = M.p0[a];
    for (int i = 0; i < nb; ++i) J[a][i] = (nuP[i] - nuM[i]) / (2 * h);
  }
  TMatrixDSym I(np);
  for (int a = 0; a < np; ++a)
    for (int b = a; b < np; ++b)
    {
      double s = 0;
      for (int i = 0; i < nb; ++i) s += J[a][i] * J[b][i] / std::max(nu0[i], 1e-9);
      I(a, b) = I(b, a) = s;
    }
  I(M.sPar, M.sPar) += M.fixS ? 1e12 : 1.0 / (M.sSig * M.sSig);
  if (M.sPar2 >= 0) I(M.sPar2, M.sPar2) += 1.0 / (M.sSig * M.sSig);
  double det = 0;
  cov.ResizeTo(np, np);
  cov = I;
  cov.Invert(&det);
  return det != 0;
}

struct FitRes
{
  bool ok = false;
  int status = -1;
  std::vector<double> p, e;
  TMatrixDSym cov;
  double chi2BC = 0;
  int nbins = 0, npar = 0;
};

// the SR multijet of the channels with sel(chan) = true, at parameters p
double QcdYield(const Model &M, const double *p, const std::function<bool(const SRChan &)> &sel)
{
  double t = 0;
  for (const SRChan &c : M.sr)
  {
    if (!sel(c)) continue;
    const double m = M.Mult(c, p);
    for (size_t k = 0; k < c.n.size(); ++k) t += m * std::exp(M.Slope(c, p) * (c.x[k] - 40.0) / 10.0) * c.Q[k];
  }
  return t;
}

// its postfit value and error (linear propagation with the fit covariance)
void QcdYieldErr(const Model &M, const FitRes &R, const std::function<bool(const SRChan &)> &sel, double &val, double &err)
{
  const int np = (int)R.p.size();
  std::vector<double> p = R.p, g(np, 0.0);
  val = QcdYield(M, p.data(), sel);
  for (int a = 0; a < np; ++a)
  {
    const double h = 1e-4 * std::max(std::fabs(R.p[a]), 0.1);
    p[a] = R.p[a] + h; const double up = QcdYield(M, p.data(), sel);
    p[a] = R.p[a] - h; const double dn = QcdYield(M, p.data(), sel);
    p[a] = R.p[a];
    g[a] = (up - dn) / (2 * h);
  }
  double v = 0;
  for (int a = 0; a < np; ++a) for (int b = 0; b < np; ++b) v += g[a] * g[b] * R.cov(a, b);
  err = std::sqrt(std::max(v, 0.0));
}

FitRes FitData(const Model &M)
{
  FitRes R;
  std::vector<double> n;
  M.Obs(n);
  const int np = (int)M.p0.size();
  std::unique_ptr<ROOT::Math::Minimizer> mn(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
  mn->SetMaxFunctionCalls(200000);
  mn->SetTolerance(0.01);
  mn->SetPrintLevel(0);
  mn->SetStrategy(1);
  ROOT::Math::Functor fcn([&](const double *p) { return M.NLL(p, n); }, np);
  mn->SetFunction(fcn);
  for (int a = 0; a < np; ++a)
  {
    if (a == M.sPar && M.fixS) mn->SetFixedVariable(a, M.pname[a], M.p0[a]);
    else mn->SetLimitedVariable(a, M.pname[a], M.p0[a], 0.01 * std::max(std::fabs(M.p0[a]), 0.1), M.plo[a], M.phi[a]);
  }
  R.ok = mn->Minimize();
  mn->Hesse();
  R.status = mn->Status();
  R.p.assign(mn->X(), mn->X() + np);
  R.e.assign(mn->Errors(), mn->Errors() + np);
  R.cov.ResizeTo(np, np);
  for (int a = 0; a < np; ++a) for (int b = 0; b < np; ++b) R.cov(a, b) = mn->CovMatrix(a, b);
  // Baker-Cousins chi2 of the SR bins
  std::vector<double> nu;
  M.Nu(R.p.data(), nu);
  size_t nsr = 0;
  for (const SRChan &c : M.sr) nsr += c.n.size();
  for (size_t i = 0; i < nsr; ++i)
  {
    const double v = std::max(nu[i], 1e-9);
    R.chi2BC += 2.0 * (v - n[i] + (n[i] > 0 ? n[i] * std::log(n[i] / v) : 0.0));
  }
  R.nbins = (int)nsr;
  R.npar = np - (M.fixS ? 1 : 0);
  return R;
}

// per pT bin of the fit window
void FillBins(const TH1D *h, std::vector<double> &v, double scale = 1.0, bool clamp = false)
{
  v.clear();
  for (int b = kFitLo; b <= kFitHi; ++b) v.push_back(clamp ? std::max(0.0, scale * h->GetBinContent(b)) : scale * h->GetBinContent(b));
}
std::vector<double> FitX(const TH1D *h)
{
  std::vector<double> x;
  for (int b = kFitLo; b <= kFitHi; ++b) x.push_back(b == kFitLo ? 25.5 : h->GetBinCenter(b));
  return x;
}

// one charge of one flavour, one m_T family: everything the schemes need
struct ChargeData
{
  std::vector<YReg> y;
  FMeas fG[kNG], fIncl;
  GroupTest gt;
};

enum Shape { kPerY = 0, kGroupShape = 1, kInclShape = 2 };

// QCD-only anti-iso template of bin iy (fit window, clamped >= 0) in one of the shape choices
void QShape(const ChargeData &cd, int iy, Shape sh, const char *reg, std::vector<double> &Q, double &tot)
{
  auto pick = [&](int j) -> const Comp & {
    const YReg &r = cd.y[j];
    return std::string(reg) == "all" ? r.A_all : (std::string(reg) == "lo" ? r.A_lo : (std::string(reg) == "hi" ? r.A_hi : r.A_40));
  };
  auto qOf = [&](int j) { TH1D *q = pick(j).Qcd(1.0, "qs"); return q; };
  TH1D *own = qOf(iy);
  std::vector<double> qy;
  FillBins(own, qy, 1.0, true);
  tot = 0;
  for (double v : qy) tot += v;
  if (sh == kPerY) { Q = qy; delete own; return; }
  const std::vector<int> bins = (sh == kGroupShape) ? GroupBins(pOAbcdY::GroupOf(iy)) : GroupBins(-1);
  std::vector<double> s(qy.size(), 0.0);
  double st = 0;
  for (int j : bins)
  {
    TH1D *q = qOf(j);
    std::vector<double> v;
    FillBins(q, v, 1.0, true);
    for (size_t k = 0; k < v.size(); ++k) { s[k] += v[k]; st += v[k]; }
    delete q;
  }
  Q.assign(qy.size(), 0.0);
  for (size_t k = 0; k < s.size(); ++k) Q[k] = st > 0 ? s[k] * tot / st : 0.0;
  delete own;
}

double Tot(const std::vector<double> &v) { double t = 0; for (double x : v) t += x; return t; }

// scheme A: SR all m_T, per-y anchor CRs (anti-iso all m_T), T_g free
Model BuildA(const ChargeData &cd, Shape sh, bool flat)
{
  Model M;
  M.name = std::string("A") + (sh == kPerY ? "" : (sh == kGroupShape ? "_groupshape" : "_inclshape")) + (flat ? "_flatF" : "");
  for (int iy = 0; iy < kNY; ++iy) M.AddPar(Form("r_y%d", iy), 1.0, 0.0, 5.0);
  int Tp[kNG];
  // prefit level: F0_g at the shared slope; flat F: the group's pT-integrated T
  for (int g = 0; g < kNG; ++g) Tp[g] = M.AddPar(Form("T_%s", pOAbcdY::kGroupName[g]), flat ? cd.fG[g].Tint.mle : cd.gt.F0shared[g], 0.0, 10.0);
  M.sPar = M.AddPar("s", flat ? 0.0 : cd.gt.sShared, -2.0, 2.0);
  M.sMean = flat ? 0.0 : cd.gt.sShared;
  M.sSig = std::max(0.02, 0.5 * (cd.fIncl.sHi - cd.fIncl.sLo));
  M.fixS = flat;
  for (int iy = 0; iy < kNY; ++iy)
  {
    const YReg &r = cd.y[iy];
    SRChan c;
    c.iy = iy;
    FillBins(r.P_all.d, c.n); FillBins(r.P_all.w, c.W); FillBins(r.P_all.o, c.O);
    double qt = 0;
    QShape(cd, iy, sh, "all", c.Q, qt);
    c.x = FitX(r.P_all.d);
    const int a = M.AddPar(Form("sC_y%d", iy), 1.0, 0.0, 5.0);
    c.mPar = {Tp[pOAbcdY::GroupOf(iy)], a};
    c.mExp = {1.0, 1.0};
    M.sr.push_back(c);
    CRChan k;
    std::vector<double> dd, ww, oo;
    FillBins(r.A_all.d, dd); FillBins(r.A_all.w, ww); FillBins(r.A_all.o, oo);
    k.n = Tot(dd);
    k.ewk = Tot(ww) + Tot(oo);
    k.base = qt;
    k.par = a;
    M.cr.push_back(k);
  }
  return M;
}

// scheme C: the SR split at m_T = 30 (lo/hi channels, shared r), anchors per y for both;
// variants: the group anti-iso shape, and an independent slope for the m_T >= 30 channels
Model BuildC(const ChargeData &cd, Shape sh = kPerY, bool twoSlopes = false)
{
  Model M;
  M.name = std::string("C") + (sh == kPerY ? "" : (sh == kGroupShape ? "_groupshape" : "_inclshape")) + (twoSlopes ? "_2s" : "");
  for (int iy = 0; iy < kNY; ++iy) M.AddPar(Form("r_y%d", iy), 1.0, 0.0, 5.0);
  int Tp[kNG];
  for (int g = 0; g < kNG; ++g) Tp[g] = M.AddPar(Form("T_%s", pOAbcdY::kGroupName[g]), cd.gt.F0shared[g], 0.0, 10.0);
  M.sPar = M.AddPar("s", cd.gt.sShared, -2.0, 2.0);
  M.sMean = cd.gt.sShared;
  M.sSig = std::max(0.02, 0.5 * (cd.fIncl.sHi - cd.fIncl.sLo));
  if (twoSlopes) M.sPar2 = M.AddPar("s_hi", cd.gt.sShared, -2.0, 2.0);
  for (int iy = 0; iy < kNY; ++iy)
  {
    const YReg &r = cd.y[iy];
    for (int part = 0; part < 2; ++part)
    {
      const Comp &P = part == 0 ? r.P_lo : r.P_hi;
      const Comp &A = part == 0 ? r.A_lo : r.A_hi;
      SRChan c;
      c.iy = iy;
      c.part = part + 1;
      FillBins(P.d, c.n); FillBins(P.w, c.W); FillBins(P.o, c.O);
      double qt = 0;
      QShape(cd, iy, sh, part == 0 ? "lo" : "hi", c.Q, qt);
      c.x = FitX(P.d);
      const int a = M.AddPar(Form("s%s_y%d", part == 0 ? "D" : "H", iy), 1.0, 0.0, 5.0);
      c.mPar = {Tp[pOAbcdY::GroupOf(iy)], a};
      c.mExp = {1.0, 1.0};
      M.sr.push_back(c);
      CRChan k;
      std::vector<double> dd, ww, oo;
      FillBins(A.d, dd); FillBins(A.w, ww); FillBins(A.o, oo);
      k.n = Tot(dd); k.ewk = Tot(ww) + Tot(oo); k.base = qt; k.par = a;
      M.cr.push_back(k);
    }
  }
  return M;
}

// scheme R: today's leppt_mt40 in-fit ABCD per group (SR m_T > 40; 1-bin B / C40 / D per group)
Model BuildR(const ChargeData &cd)
{
  Model M;
  M.name = "R";
  for (int iy = 0; iy < kNY; ++iy) M.AddPar(Form("r_y%d", iy), 1.0, 0.0, 5.0);
  M.sPar = M.AddPar("s", cd.gt.sShared, -2.0, 2.0);
  M.sMean = cd.gt.sShared;
  M.sSig = std::max(0.02, 0.5 * (cd.fIncl.sHi - cd.fIncl.sLo));
  for (int g = 0; g < kNG; ++g)
  {
    const int pB = M.AddPar(Form("sB_%s", pOAbcdY::kGroupName[g]), 1.0, 0.0, 5.0);
    const int pC = M.AddPar(Form("sC_%s", pOAbcdY::kGroupName[g]), 1.0, 0.0, 5.0);
    const int pD = M.AddPar(Form("sD_%s", pOAbcdY::kGroupName[g]), 1.0, 0.0, 5.0);
    // group B / C40 / D counts (pT >= 25 incl. overflow, the qcd_abcd.C convention)
    double dB = 0, oB = 0, dC = 0, eC = 0, dD = 0, eD = 0, qD = 0, fD = 0, qC = 0;
    std::vector<double> wB;
    for (int iy : GroupBins(g))
    {
      const YReg &r = cd.y[iy];
      dB += Cnt(r.P_lo.d, 24, kBig); oB += Cnt(r.P_lo.o, 24, kBig); wB.push_back(Cnt(r.P_lo.w, 24, kBig));
      dC += Cnt(r.A_40.d, 24, kBig); eC += Cnt(r.A_40.w, 24, kBig) + Cnt(r.A_40.o, 24, kBig);
      dD += Cnt(r.A_lo.d, 24, kBig); eD += Cnt(r.A_lo.w, 24, kBig) + Cnt(r.A_lo.o, 24, kBig);
      TH1D *q = r.A_lo.Qcd(1.0, "qD");
      for (int b = q->FindBin(24.5); b <= q->GetNbinsX() + 1; ++b)
      {
        const double xc = (b == q->FindBin(24.5)) ? 25.5 : (b > q->GetNbinsX() ? 105.0 : q->GetBinCenter(b));
        const double v = std::max(0.0, q->GetBinContent(b));
        qD += v;
        fD += v * std::exp(cd.gt.sShared * (xc - 40.0) / 10.0);
      }
      delete q;
    }
    double qB = dB - oB; // QCD-only B / C40 / D at r = 1 (not "B0": a termios.h macro)
    for (double v : wB) qB -= v;
    const double qC0 = dC - eC, qD0 = dD - eD;
    const double fbar = qD > 0 ? fD / qD : 1.0; // <f>_D: the FF template is normalized on D (buildFFTemplate)
    for (int iy : GroupBins(g))
    {
      const YReg &r = cd.y[iy];
      SRChan c;
      c.iy = iy;
      FillBins(r.P_40.d, c.n); FillBins(r.P_40.w, c.W); FillBins(r.P_40.o, c.O);
      TH1D *q = r.A_40.Qcd(1.0, "q40");
      FillBins(q, c.Q, 1.0, true);
      delete q;
      qC += Tot(c.Q);
      c.x = FitX(r.P_40.d);
      c.mPar = {pB, pC, pD};
      c.mExp = {1.0, 1.0, -1.0};
      c.mConst = (qD0 > 0 ? qB / qD0 : 0.0) / fbar;
      M.sr.push_back(c);
    }
    CRChan kB; kB.n = dB; kB.base = qB; kB.ewk = oB; kB.par = pB;
    std::vector<int> ys = GroupBins(g);
    for (size_t i = 0; i < ys.size(); ++i) { kB.rPar.push_back(ys[i]); kB.rW.push_back(wB[i]); }
    CRChan kC; kC.n = dC; kC.base = qC0; kC.ewk = eC; kC.par = pC;
    CRChan kD; kD.n = dD; kD.base = qD0; kD.ewk = eD; kD.par = pD;
    M.cr.push_back(kB); M.cr.push_back(kC); M.cr.push_back(kD);
    (void)qC;
  }
  return M;
}

int ParIndex(const Model &M, const std::string &nm)
{
  for (size_t i = 0; i < M.pname.size(); ++i) if (M.pname[i] == nm) return (int)i;
  return -1;
}

// ============================================================================
// Driver pieces
// ============================================================================
struct Run
{
  // fitted r_y and errors per scheme name, for [EMU]
  std::map<std::string, std::vector<double>> r, e;
  // per-charge numbers for the e+ / e- comparison [CHG]: key "<fam>_<what>_<group>_<chg>"
  std::map<std::string, double> v;
};

void PrintComp(const LepCfg &L, const std::string &chg, const std::string &fam, const ChargeData &cd)
{
  std::cout << "\n[COMP] " << L.lep << chg << " (" << fam << "): no-m_T SR (iso-pass, pT 24-100) vs the m_T>40 SR, prefit\n"
            << "[COMP]   bin |   data      W     other   QCD(FxC) QCD/data |  data40  QCD40/data40\n";
  double S[6] = {0, 0, 0, 0, 0, 0};
  for (int iy = 0; iy < kNY; ++iy)
  {
    const YReg &r = cd.y[iy];
    const int g = pOAbcdY::GroupOf(iy);
    auto qcdPred = [&](const Comp &A) {
      TH1D *q = A.Qcd(1.0, "qp");
      double t = 0;
      for (int b = kFitLo; b <= kFitHi; ++b)
        t += std::max(0.0, q->GetBinContent(b)) * cd.gt.F0shared[g] *
             std::exp(cd.gt.sShared * (((b == kFitLo) ? 25.5 : q->GetBinCenter(b)) - 40.0) / 10.0);
      delete q;
      return t;
    };
    const double d = Cnt(r.P_all.d, 24, 100), w = Cnt(r.P_all.w, 24, 100), o = Cnt(r.P_all.o, 24, 100), q = qcdPred(r.A_all);
    const double d4 = Cnt(r.P_40.d, 24, 100), q4 = qcdPred(r.A_40);
    std::cout << Form("[COMP]   y%-2d | %6.0f %7.1f %7.1f %8.1f   %5.3f  | %6.0f   %5.3f\n", iy, d, w, o, q, d > 0 ? q / d : 0., d4, d4 > 0 ? q4 / d4 : 0.);
    S[0] += d; S[1] += w; S[2] += o; S[3] += q; S[4] += d4; S[5] += q4;
    gCsv << "comp," << L.lep << chg << "," << fam << ",y" << iy << "," << d << "," << w << "," << o << "," << q << "," << d4 << "," << q4 << "\n";
  }
  std::cout << Form("[COMP]   all | %6.0f %7.1f %7.1f %8.1f   %5.3f  | %6.0f   %5.3f   (data/(W+other+QCD) = %.3f)\n",
                    S[0], S[1], S[2], S[3], S[3] / S[0], S[4], S[5] / S[4], S[0] / (S[1] + S[2] + S[3]));
  // region B (iso-pass, m_T < 30) per group, and how much of the SR W the MC puts below m_T = 30
  // (what an m_T-split SR would take from the MC m_T spectrum)
  for (int g = -1; g < kNG; ++g)
  {
    double dB = 0, wB = 0, oB = 0, wAll = 0;
    for (int iy : GroupBins(g))
    {
      dB += Cnt(cd.y[iy].P_lo.d, 24, kBig); wB += Cnt(cd.y[iy].P_lo.w, 24, kBig); oB += Cnt(cd.y[iy].P_lo.o, 24, kBig);
      wAll += Cnt(cd.y[iy].P_all.w, 24, kBig);
    }
    std::cout << Form("[BCOMP] %s%s %-4s B (iso-pass, m_T<30): data %5.0f, W MC %6.1f (%4.1f%%), other %5.1f, QCD %6.1f;  W MC at m_T<30 = %4.1f%% of the SR W\n",
                      L.lep.c_str(), chg.c_str(), GroupName(g).c_str(), dB, wB, 100 * wB / dB, oB, dB - wB - oB, 100 * wB / wAll);
  }
}

void PrintF(const LepCfg &L, const std::string &chg, const std::string &fam, const ChargeData &cd)
{
  std::cout << "\n[FPT] " << L.lep << chg << " (" << fam << "): F(pT) = QCD(iso-pass)/QCD(anti-iso) at m_T < " << kMtB
            << ", 68% profile-likelihood intervals; F = F0 exp(s (pT-40)/10)\n";
  std::cout << "[FPT]   pT bin     ";
  for (int g = -1; g < kNG; ++g) std::cout << Form("| %-22s", GroupName(g).c_str());
  std::cout << "\n";
  for (size_t k = 0; k + 1 < kFE.size(); ++k)
  {
    const std::string hiEdge = kFE[k + 1] > 1e8 ? std::string("inf") : std::string(Form("%3.0f", kFE[k + 1]));
    std::cout << Form("[FPT]   [%3.0f,%s) ", kFE[k], hiEdge.c_str());
    for (int g = -1; g < kNG; ++g)
    {
      const FMeas &m = g < 0 ? cd.fIncl : cd.fG[g];
      std::cout << Form("| %5.3f [%5.3f,%5.3f] %4.0f ", m.bin[k].mle, m.bin[k].lo, m.bin[k].hi, m.np[k]);
    }
    std::cout << "\n";
  }
  for (int g = -1; g < kNG; ++g)
  {
    const FMeas &m = g < 0 ? cd.fIncl : cd.fG[g];
    const int n = (int)m.x.size();
    std::cout << Form("[FPT]   %-4s: T(pT-integrated) %.4f [%.4f,%.4f]; slope s = %+.3f [%+.3f,%+.3f] (x%.2f per 10 GeV, %.1f sigma from flat),"
                      " F0 = %.4f; LRT flat: slope p = %s, free p = %s\n",
                      GroupName(g).c_str(), m.Tint.mle, m.Tint.lo, m.Tint.hi, m.fl.s, m.sLo, m.sHi, std::exp(m.fl.s),
                      std::sqrt(std::max(0.0, m.fl.qSlope())), m.fl.F0s, ffFmtP(TMath::Prob(m.fl.qSlope(), 1)).c_str(),
                      ffFmtP(TMath::Prob(m.fl.qFree(), n - 1)).c_str());
    gCsv << "fpt," << L.lep << chg << "," << fam << "," << GroupName(g) << "," << m.Tint.mle << "," << m.fl.s << "," << m.sLo << ","
         << m.sHi << "," << m.fl.F0s << "," << TMath::Prob(m.fl.qSlope(), 1) << "\n";
  }
  const GroupTest &t = cd.gt;
  const double qN = 2 * (t.lnF0g - t.lnCommon), qS = 2 * (t.lnFull - t.lnF0g);
  std::cout << Form("[FPT]   group dependence: level F0 per group (common s) vs common: q = %.2f, p = %s;  slope per group vs common: q = %.2f, p = %s\n",
                    qN, ffFmtP(TMath::Prob(qN, 2)).c_str(), qS, ffFmtP(TMath::Prob(qS, 2)).c_str());
  std::cout << Form("[FPT]   shared slope s = %+.3f; F0_g (68%%):", t.sShared);
  for (int g = 0; g < kNG; ++g) std::cout << Form("  %s %.4f [%.4f,%.4f]", pOAbcdY::kGroupName[g], t.F0shared[g], t.F0lo[g], t.F0hi[g]);
  std::cout << Form("   (F0_g / F0_cen: fwd %.2f, bwd %.2f)\n", t.F0shared[0] / t.F0shared[2], t.F0shared[1] / t.F0shared[2]);
}

void RunCharge(const LepCfg &L, const Files &F, const std::string &chg, Run &run, const std::string &outDir)
{
  const std::string chgLatex = (chg == "Plus") ? (L.isElec ? "e^{+}" : "#mu^{+}") : (L.isElec ? "e^{-}" : "#mu^{-}");
  for (const std::string fam : {"mt", "mttrk"})
  {
    const std::string famName = fam == "mt" ? "m_T from PF MET" : "m_T from charged-only MET";
    ChargeData cd;
    bool ok = true;
    for (int iy = 0; iy < kNY; ++iy)
    {
      const std::string nm = Form("h_iso_pt_%s_%s%s_y%d", fam.c_str(), L.lep.c_str(), chg.c_str(), iy);
      Cube c;
      TH3D *hd = (TH3D *)F.data->Get(nm.c_str());
      if (!hd) { std::cerr << "[ERROR] " << nm << " missing in the data skim (skim before 2026-10-21?)\n"; return; }
      c.d = (TH3D *)hd->Clone(Uid(nm + "_d").c_str());
      c.d->SetDirectory(nullptr);
      c.w = SumCube(F.w, nm, ok);
      c.o = SumCube(F.o, nm, ok);
      if (!ok) return;
      cd.y.push_back(MakeYReg(c, L, ok));
      if (!ok) return;
    }
    // F(pT) at low m_T per group and inclusive
    for (int g = -1; g < kNG; ++g)
    {
      std::vector<const Comp *> pv, av;
      for (int iy : GroupBins(g)) { pv.push_back(&cd.y[iy].P_lo); av.push_back(&cd.y[iy].A_lo); }
      const FMeas m = MeasureF(SumComp(pv), SumComp(av), 1.0);
      if (g < 0) cd.fIncl = m; else cd.fG[g] = m;
    }
    cd.gt = TestGroups(cd.fG);
    for (int g = -1; g < kNG; ++g)
    {
      const FMeas &m = g < 0 ? cd.fIncl : cd.fG[g];
      const std::string k = "_" + GroupName(g) + "_" + chg;
      run.v[fam + "_T" + k] = m.Tint.mle; run.v[fam + "_Tlo" + k] = m.Tint.lo; run.v[fam + "_Thi" + k] = m.Tint.hi;
      run.v[fam + "_s" + k] = m.fl.s;     run.v[fam + "_slo" + k] = m.sLo;     run.v[fam + "_shi" + k] = m.sHi;
    }
    std::cout << "\n=============== " << L.lep << chg << " -- " << famName << " ===============\n";
    PrintComp(L, chg, fam, cd);
    PrintF(L, chg, fam, cd);

    // anti-iso window variation of F (inclusive, m_T < 30)
    {
      std::vector<const Comp *> pv, a1, a2;
      for (int iy = 0; iy < kNY; ++iy) { pv.push_back(&cd.y[iy].P_lo); a1.push_back(&cd.y[iy].A1_lo); a2.push_back(&cd.y[iy].A2_lo); }
      const FMeas m1 = MeasureF(SumComp(pv), SumComp(a1), 1.0), m2 = MeasureF(SumComp(pv), SumComp(a2), 1.0);
      std::cout << Form("[FPT]   anti-iso window: [%.2f,%.2f) T %.4f, s %+.3f [%+.3f,%+.3f] | [%.2f,%.2f) T %.4f, s %+.3f [%+.3f,%+.3f] | nominal s %+.3f\n",
                        L.aLo, L.aMid, m1.Tint.mle, m1.fl.s, m1.sLo, m1.sHi, L.aMid, L.aHi, m2.Tint.mle, m2.fl.s, m2.sLo, m2.sHi, cd.fIncl.fl.s);
    }
    // prefit W-scale dependence of the slope
    {
      std::vector<const Comp *> pv, av;
      for (int iy = 0; iy < kNY; ++iy) { pv.push_back(&cd.y[iy].P_lo); av.push_back(&cd.y[iy].A_lo); }
      const FMeas mr = MeasureF(SumComp(pv), SumComp(av), kRAlt);
      std::cout << Form("[FPT]   at r_W = %.1f: T %.4f (%+.1f%%), s %+.3f (%+.2f sigma)\n", kRAlt, mr.Tint.mle,
                        100 * (mr.Tint.mle / cd.fIncl.Tint.mle - 1), mr.fl.s,
                        (mr.fl.s - cd.fIncl.fl.s) / std::max(1e-6, 0.5 * (cd.fIncl.sHi - cd.fIncl.sLo)));
    }

    // [TY] per-bin pT-integrated T_y and [SHAPE]
    if (fam == "mt")
    {
      std::cout << "\n[TY] " << L.lep << chg << ": pT-integrated T_y = QCD(iso)/QCD(anti) at m_T < 30 per rapidity bin (68%)\n[TY]  ";
      for (int iy = 0; iy < kNY; ++iy)
      {
        const FMeas m = MeasureF(cd.y[iy].P_lo, cd.y[iy].A_lo, 1.0);
        std::cout << Form(" y%d %.3f[%.3f,%.3f]", iy, m.Tint.mle, m.Tint.lo, m.Tint.hi) << (iy == 5 ? "\n[TY]  " : "");
        gCsv << "ty," << L.lep << chg << ",y" << iy << "," << m.Tint.mle << "," << m.Tint.lo << "," << m.Tint.hi << "\n";
      }
      std::cout << "\n";
      // anti-iso all-m_T QCD-only pT shapes: bin vs the rest of the acceptance and vs the rest of its group
      std::cout << "[SHAPE] " << L.lep << chg << ": anti-iso (all m_T) pT shape of each bin vs the rest (chi2 homogeneity, WW, kFE bins):\n[SHAPE]  ";
      auto coarse = [&](const TH1D *h) {
        TH1D *c = new TH1D(Uid("cs").c_str(), "", (int)kFE.size() - 1, 0, (double)kFE.size() - 1);
        c->Sumw2();
        for (size_t k = 0; k + 1 < kFE.size(); ++k) { double v2 = 0; const double v = Cnt(h, kFE[k], kFE[k + 1], &v2); c->SetBinContent(k + 1, std::max(0.0, v)); c->SetBinError(k + 1, std::sqrt(v2)); }
        return c;
      };
      std::vector<TH1D *> qy(kNY);
      for (int iy = 0; iy < kNY; ++iy) qy[iy] = coarse(cd.y[iy].A_all.Qcd(1.0, "qa"));
      int nSmall = 0;
      for (int iy = 0; iy < kNY; ++iy)
      {
        TH1D *rest = (TH1D *)qy[iy]->Clone(Uid("rest").c_str()); rest->Reset();
        TH1D *restG = (TH1D *)qy[iy]->Clone(Uid("restg").c_str()); restG->Reset();
        for (int j = 0; j < kNY; ++j)
        {
          if (j == iy) continue;
          rest->Add(qy[j]);
          if (pOAbcdY::GroupOf(j) == pOAbcdY::GroupOf(iy)) restG->Add(qy[j]);
        }
        const double p = qy[iy]->Chi2Test(rest, "WW");
        const double pg = (restG->Integral() > 0) ? qy[iy]->Chi2Test(restG, "WW") : -1;
        if (p < 0.05) ++nSmall;
        const std::string pgs = pg < 0 ? std::string("-") : std::string(Form("%.3f", pg));
        std::cout << Form(" y%d p=%.3f/%s", iy, p, pgs.c_str()) << (iy == 5 ? "\n[SHAPE]  " : "");
        gCsv << "shape," << L.lep << chg << ",y" << iy << "," << p << "," << pg << "\n";
      }
      std::cout << "\n[SHAPE]   (p vs the rest of the acceptance / vs the rest of its group; " << nSmall << " of 12 below 0.05)\n";
    }

    // [MTFAC] F vs m_T per group (pT-integrated, pT >= 25) and the SR multijet fraction at m_T >= 30
    std::cout << "[MTFAC] " << L.lep << chg << " (" << fam << "): F(m_T slice)/F(m_T<30), pT-integrated, at r_W = 1 | " << kRAlt << "\n";
    for (int g = -1; g < kNG; ++g)
    {
      std::cout << Form("[MTFAC]   %-4s:", GroupName(g).c_str());
      std::vector<const Comp *> pl, al;
      for (int iy : GroupBins(g)) { pl.push_back(&cd.y[iy].P_lo); al.push_back(&cd.y[iy].A_lo); }
      const Comp Pl = SumComp(pl), Al = SumComp(al);
      for (double rW : {1.0, kRAlt})
      {
        const FMeas ref = MeasureF(Pl, Al, rW);
        for (int m = 0; m < 4; ++m)
        {
          std::vector<const Comp *> pv, av;
          for (int iy : GroupBins(g)) { pv.push_back(&cd.y[iy].P_m[m]); av.push_back(&cd.y[iy].A_m[m]); }
          const Comp P = SumComp(pv), A = SumComp(av);
          const double np = Cnt(P.d, 24, kBig), nf = Cnt(A.d, 24, kBig);
          const double bp = std::max(1e-6, rW * Cnt(P.w, 24, kBig) + Cnt(P.o, 24, kBig));
          const double bf = std::max(1e-6, rW * Cnt(A.w, 24, kBig) + Cnt(A.o, 24, kBig));
          const FFInterval I = ffInterval(np, nf, bp, bf);
          const double R0 = ref.Tint.mle;
          if (m == 3 && rW == 1.0)
          {
            const std::string k = "_" + GroupName(g) + "_" + chg;
            run.v[fam + "_f34" + k] = I.mle / R0; run.v[fam + "_f34lo" + k] = I.lo / R0; run.v[fam + "_f34hi" + k] = I.hi / R0;
          }
          std::cout << Form(" [%2.0f,%2.0f) %.2f[%.2f,%.2f]", kMtSl[m], kMtSl[m + 1], I.mle / R0, I.lo / R0, I.hi / R0);
          gCsv << "mtfac," << L.lep << chg << "," << fam << "," << GroupName(g) << "," << rW << "," << kMtSl[m] << "," << I.mle / R0 << ","
               << I.lo / R0 << "," << I.hi / R0 << "," << (np > 0 ? bp / np : 0) << "\n";
        }
        std::cout << (rW == 1.0 ? "  |" : "");
      }
      std::cout << "\n";
    }
    {
      double lo = 0, hi = 0, all = 0;
      for (int iy = 0; iy < kNY; ++iy)
      {
        const int g = pOAbcdY::GroupOf(iy);
        TH1D *qa = cd.y[iy].A_all.Qcd(1.0, "qa"), *qh = cd.y[iy].A_hi.Qcd(1.0, "qh");
        for (int b = kFitLo; b <= kFitHi; ++b)
        {
          const double f = cd.gt.F0shared[g] * std::exp(cd.gt.sShared * ((b == kFitLo ? 25.5 : qa->GetBinCenter(b)) - 40.0) / 10.0);
          all += f * std::max(0.0, qa->GetBinContent(b));
          hi += f * std::max(0.0, qh->GetBinContent(b));
        }
        delete qa; delete qh;
      }
      lo = all - hi;
      std::cout << Form("[MTFAC]   SR multijet (F x C, fit window): %.1f, of which m_T < 30 %.1f (%.0f%%), m_T >= 30 %.1f (%.0f%%)\n",
                        all, lo, 100 * lo / all, hi, 100 * hi / all);
      run.v[fam + "_fhi_incl_" + chg] = hi / all;
    }

    // [FISHER] + [FIT]
    std::vector<Model> models = {BuildA(cd, kPerY, false), BuildA(cd, kGroupShape, false), BuildA(cd, kInclShape, false),
                                 BuildA(cd, kPerY, true), BuildC(cd), BuildC(cd, kGroupShape, false), BuildC(cd, kPerY, true),
                                 BuildR(cd)};
    std::cout << "[FISHER] " << L.lep << chg << " (" << fam << "): expected (Asimov, prefit) errors -- sigma(r_y) per bin, then T_g / scales\n";
    std::map<std::string, std::vector<double>> expR;
    for (const Model &M : models)
    {
      if (M.name != "A" && M.name != "C" && M.name != "R") continue;
      if (fam == "mttrk" && M.name == "A") continue; // scheme A has no m_T in it: identical for both families
      TMatrixDSym cov;
      if (!FisherCov(M, cov)) { std::cout << "[FISHER]   " << M.name << ": singular Fisher matrix\n"; continue; }
      std::cout << Form("[FISHER]   %-2s r:", M.name.c_str());
      double sw = 0, se2 = 0;
      for (int iy = 0; iy < kNY; ++iy)
      {
        const double e = std::sqrt(cov(iy, iy));
        expR[M.name].push_back(e);
        std::cout << Form(" %.3f", e);
      }
      // the W-weighted average r (~ sigma_W / sigma_gen): its expected error with the full covariance
      std::vector<double> wts(kNY, 0.0);
      for (const SRChan &c : M.sr) wts[c.iy] += Tot(c.W);
      for (double v : wts) sw += v;
      for (int a = 0; a < kNY; ++a) for (int b = 0; b < kNY; ++b) se2 += wts[a] * wts[b] * cov(a, b) / (sw * sw);
      std::cout << Form("  | <r>_W %.4f |", std::sqrt(se2));
      for (size_t a = kNY; a < M.pname.size(); ++a)
        if (M.pname[a].rfind("T_", 0) == 0 || M.pname[a].rfind("sB", 0) == 0 || M.pname[a] == "s")
          std::cout << Form(" %s %.1f%%", M.pname[a].c_str(), 100 * std::sqrt(cov(a, a)) / std::max(1e-9, std::fabs(M.p0[a])));
      if (M.name == "A") // the r-T correlation of scheme A, per group (mean over its bins)
        for (int g = 0; g < kNG; ++g)
        {
          const int t = ParIndex(M, Form("T_%s", pOAbcdY::kGroupName[g]));
          double rho = 0; int n = 0;
          for (int iy : GroupBins(g)) { rho += cov(iy, t) / std::sqrt(cov(iy, iy) * cov(t, t)); ++n; }
          std::cout << Form(" rho(r,T_%s) %+.2f", pOAbcdY::kGroupName[g], rho / n);
        }
      std::cout << "\n";
      gCsv << "fisher," << L.lep << chg << "," << fam << "," << M.name << ",avgr," << std::sqrt(se2) << "\n";
    }
    std::cout << "[FIT] " << L.lep << chg << " (" << fam << "): data stand-in fits -- r_y, then the fitted T_g against the m_T<30 measurement F0_g\n";
    for (const Model &M : models)
    {
      if (fam == "mttrk" && M.name.rfind("A", 0) == 0) continue; // scheme A has no m_T in it
      const FitRes R = FitData(M);
      std::cout << Form("[FIT]   %-13s status %d  chi2_BC/ndf %.1f/%d  r:", M.name.c_str(), R.status, R.chi2BC, R.nbins - R.npar);
      for (int iy = 0; iy < kNY; ++iy) std::cout << Form(" %.3f", R.p[iy]);
      std::cout << "\n[FIT]   " << std::string(13, ' ') << "                         +-";
      for (int iy = 0; iy < kNY; ++iy) std::cout << Form(" %.3f", R.e[iy]);
      std::cout << "\n";
      if (M.name.rfind("A", 0) == 0 || M.name.rfind("C", 0) == 0)
      {
        // THE closure test: the SR multijet the fit wants (its level free, scheme A; pinned by the
        // m_T < 30 channel, scheme C) against the prefit F x C = the m_T < 30 B/D measurement,
        // per group; the prefit error = the relative error of that measurement (F0_g at the shared
        // slope; the pT-integrated T_g for the flat variant)
        for (int part = 0; part <= (M.name.rfind("C", 0) == 0 ? 1 : 0); ++part)
        {
          std::cout << "[FIT]   " << std::string(13, ' ') << (part == 0 ? " SR multijet postfit/prefit:" : " ... its m_T<30 part:    ");
          for (int g = 0; g < kNG; ++g)
          {
            auto sel = [&](const SRChan &c) { return pOAbcdY::GroupOf(c.iy) == g && (part == 0 || c.part == 1); };
            double post = 0, postE = 0;
            QcdYieldErr(M, R, sel, post, postE);
            const double pre = QcdYield(M, M.p0.data(), sel);
            const double relRef = M.fixS ? 0.5 * (cd.fG[g].Tint.hi - cd.fG[g].Tint.lo) / cd.fG[g].Tint.mle
                                         : 0.5 * (cd.gt.F0hi[g] - cd.gt.F0lo[g]) / cd.gt.F0shared[g];
            const double preE = relRef * pre;
            std::cout << Form("  %s %.0f+-%.0f vs %.0f+-%.0f (x%.2f, %+.1f sigma)", pOAbcdY::kGroupName[g], post, postE, pre, preE,
                              post / pre, (post - pre) / std::sqrt(postE * postE + preE * preE));
            gCsv << "closure," << L.lep << chg << "," << fam << "," << M.name << "," << part << "," << pOAbcdY::kGroupName[g] << ","
                 << post << "," << postE << "," << pre << "," << preE << "\n";
          }
          if (part == 0) std::cout << Form("  s %+.3f (prior %+.3f+-%.3f)", R.p[M.sPar], M.sMean, M.fixS ? 0.0 : M.sSig);
          if (part == 0 && M.sPar2 >= 0) std::cout << Form("  s(m_T>=30) %+.3f", R.p[M.sPar2]);
          std::cout << "\n";
        }
      }
      else
      {
        std::cout << "[FIT]   " << std::string(13, ' ') << " SR(m_T>40) multijet postfit/prefit:";
        for (int g = 0; g < kNG; ++g)
        {
          auto sel = [&](const SRChan &c) { return pOAbcdY::GroupOf(c.iy) == g; };
          double post = 0, postE = 0;
          QcdYieldErr(M, R, sel, post, postE);
          const double pre = QcdYield(M, M.p0.data(), sel);
          std::cout << Form("  %s %.0f+-%.0f vs %.0f (x%.2f)", pOAbcdY::kGroupName[g], post, postE, pre, post / pre);
        }
        std::cout << "\n[FIT]   " << std::string(13, ' ') << " group multipliers sB*sC/sD:";
        for (int g = 0; g < kNG; ++g)
        {
          const int b = ParIndex(M, Form("sB_%s", pOAbcdY::kGroupName[g])), c = ParIndex(M, Form("sC_%s", pOAbcdY::kGroupName[g])),
                    d = ParIndex(M, Form("sD_%s", pOAbcdY::kGroupName[g]));
          std::cout << Form("  %s %.3f", pOAbcdY::kGroupName[g], R.p[b] * R.p[c] / R.p[d]);
        }
        std::cout << Form("  s %+.3f\n", R.p[M.sPar]);
      }
      const std::string key = M.name + (fam == "mttrk" ? "_trk" : "");
      run.r[key + "_" + chg] = std::vector<double>(R.p.begin(), R.p.begin() + kNY);
      run.e[key + "_" + chg] = std::vector<double>(R.e.begin(), R.e.begin() + kNY);
      gCsv << "fit," << L.lep << chg << "," << fam << "," << M.name << "," << R.status << "," << R.chi2BC << "," << R.nbins - R.npar;
      for (int iy = 0; iy < kNY; ++iy) gCsv << "," << R.p[iy] << "," << R.e[iy];
      gCsv << "\n";

      // the postfit rapidity-inclusive stack of scheme A (per-y shape)
      if (M.name == "A")
      {
        TH1D *hD = (TH1D *)cd.y[0].P_all.d->Clone(Uid("pfD").c_str()); hD->Reset(); hD->SetDirectory(nullptr);
        TH1D *hW = (TH1D *)hD->Clone(Uid("pfW").c_str()), *hQ = (TH1D *)hD->Clone(Uid("pfQ").c_str()), *hO = (TH1D *)hD->Clone(Uid("pfO").c_str());
        for (const SRChan &c : M.sr)
        {
          const double mlt = M.Mult(c, R.p.data());
          for (size_t k = 0; k < c.n.size(); ++k)
          {
            const int b = kFitLo + (int)k;
            hD->SetBinContent(b, hD->GetBinContent(b) + c.n[k]);
            hW->SetBinContent(b, hW->GetBinContent(b) + R.p[c.iy] * c.W[k]);
            hQ->SetBinContent(b, hQ->GetBinContent(b) + mlt * std::exp(R.p[M.sPar] * (c.x[k] - 40.0) / 10.0) * c.Q[k]);
            hO->SetBinContent(b, hO->GetBinContent(b) + c.O[k]);
          }
        }
        for (int b = 1; b <= hD->GetNbinsX(); ++b) { hD->SetBinError(b, std::sqrt(hD->GetBinContent(b))); hW->SetBinError(b, 0); hQ->SetBinError(b, 0); hO->SetBinError(b, 0); }
        PlotStyle ps;
        ps.normBkgToData = false;
        ps.pullPad = true;
        ps.xRangeLo = 24; ps.xRangeHi = 100;
        ps.legX1 = 0.55; ps.legY1 = 0.50; ps.legX2 = 0.92; ps.legY2 = 0.66;
        SaveNicePlot1D_WithBkg(hD, {hW, hO, hQ}, {"W (r per bin)", "DY, #tau, t#bar{t}", "multijet F_{g}(p_{T}) #times C_{y}(p_{T})"},
                               outDir + "/postfitA_" + L.lep + chg, Form("p_{T}^{%s} [GeV]", L.sym.c_str()), "Events / 2 GeV",
                               "no-m_{T} lepton-p_{T} stand-in fit (scheme A)", chgLatex + ", all rapidity bins",
                               "stat. only, no shape systematics",
                               {Form("data %.0f", hD->Integral()), Form("W %.1f, multijet %.1f", hW->Integral(), hQ->Integral()),
                                Form("BC #chi^{2}/ndf = %.0f/%d", R.chi2BC, R.nbins - R.npar)},
                               ps);
      }
    }

    // the F(pT) plot per group (PF-MET m_T only)
    if (fam == "mt")
    {
      std::vector<OverlaySeries> ser;
      const int col[kNG + 1] = {kBlack, kAzure + 2, kRed + 1, kGreen + 2};
      const int mk[kNG + 1] = {20, 21, 33, 34};
      const double ms[kNG + 1] = {1.3, 1.5, 2.0, 1.5};
      const double sh[kNG + 1] = {0.0, -0.45, 0.45, 0.0};
      for (int g = -1; g < kNG; ++g)
      {
        if (g == 2) continue; // cen ~ incl; keep the plot readable: incl, fwd, bwd
        const FMeas &m = g < 0 ? cd.fIncl : cd.fG[g];
        const int n = (int)m.x.size() - 1; // drop the > 100 overflow bin from the plot
        TGraphErrors *gt = new TGraphErrors(n);
        TGraphAsymmErrors *ga = new TGraphAsymmErrors(n);
        for (int k = 0; k < n; ++k)
        {
          const double xc = 0.5 * (kFE[k] + kFE[k + 1]), ex = 0.5 * (kFE[k + 1] - kFE[k]);
          gt->SetPoint(k, xc, m.bin[k].mle); gt->SetPointError(k, ex, 0);
          ga->SetPoint(k, xc, m.bin[k].mle); ga->SetPointError(k, 0, 0, m.bin[k].mle - m.bin[k].lo, m.bin[k].hi - m.bin[k].mle);
        }
        OverlaySeries s;
        s.gTot = gt; s.gTotAsym = ga;
        s.label = Form("%s: s = %+.2f^{%+.2f}_{%+.2f}", GroupName(g).c_str(), m.fl.s, m.sHi - m.fl.s, m.sLo - m.fl.s);
        s.color = col[g + 1]; s.marker = mk[g + 1]; s.markerSize = ms[g + 1]; s.xShift = sh[g + 1];
        ser.push_back(s);
      }
      PlotStyle ps;
      ps.headerX = 0.47; ps.headerY = 0.89; ps.titleSize = 0.036; ps.subSize = 0.030; ps.headerDy = 0.042;
      ps.boxX1 = 0.47; ps.boxX2 = 0.93; ps.boxY1 = 0.66; ps.boxY2 = 0.76; ps.boxTextSize = 0.027;
      ps.legX1 = 0.47; ps.legX2 = 0.92; ps.legY1 = 0.49; ps.legY2 = 0.64;
      const double yMax = L.isElec ? 2.6 : 0.9;
      SaveNiceGraph_Overlay(ser, outDir + "/ff_groups_" + L.lep + chg, Form("p_{T}^{%s} [GeV]", L.sym.c_str()),
                            "F = QCD(iso-pass) / QCD(anti-iso)", "F(p_{T}) at m_{T} < 30 GeV", chgLatex + ", PF MET",
                            "68% profile-likelihood intervals",
                            {Form("level F0_{g}/F0_{cen}: fwd %.2f, bwd %.2f", cd.gt.F0shared[0] / cd.gt.F0shared[2], cd.gt.F0shared[1] / cd.gt.F0shared[2]),
                             Form("group LRT: level p = %s, slope p = %s", ffFmtP(TMath::Prob(2 * (cd.gt.lnF0g - cd.gt.lnCommon), 2)).c_str(),
                                  ffFmtP(TMath::Prob(2 * (cd.gt.lnFull - cd.gt.lnF0g), 2)).c_str())},
                            ps, [yMax](TCanvas *, TGraphErrors *g) { g->SetMinimum(0.0); g->SetMaximum(yMax); });
    }
  }
}

// [CHG] the + / - comparison of the ABCD ingredients, PF MET vs charged-only MET (user question
// 2026-10-21): pT-integrated T at m_T < 30, the F(pT) slope, F at m_T 30-40 over m_T < 30, and the
// share of the SR multijet at m_T >= 30
void PrintChg(const std::string &lep, const Run &run)
{
  std::cout << "\n=============== [CHG] " << lep << ": + vs - (pT >= 25; 68% intervals) -- PF MET | charged-only MET ===============\n";
  auto get = [&](const std::string &k) { auto it = run.v.find(k); return it == run.v.end() ? NAN : it->second; };
  for (const std::string fam : {"mt", "mttrk"})
  {
    std::cout << "[CHG] " << (fam == "mt" ? "PF MET" : "charged-only MET") << "\n";
    for (int g = -1; g < kNG; ++g)
    {
      const std::string G = GroupName(g);
      auto kk = [&](const std::string &w, const std::string &c) { return fam + "_" + w + "_" + G + "_" + c; };
      const double Tp = get(kk("T", "Plus")), Tm = get(kk("T", "Minus"));
      const double eTp = 0.5 * (get(kk("Thi", "Plus")) - get(kk("Tlo", "Plus"))), eTm = 0.5 * (get(kk("Thi", "Minus")) - get(kk("Tlo", "Minus")));
      const double rat = Tp / Tm, eRat = rat * std::hypot(eTp / Tp, eTm / Tm);
      const double sp = get(kk("s", "Plus")), sm = get(kk("s", "Minus"));
      const double esp = 0.5 * (get(kk("shi", "Plus")) - get(kk("slo", "Plus"))), esm = 0.5 * (get(kk("shi", "Minus")) - get(kk("slo", "Minus")));
      const double fp = get(kk("f34", "Plus")), fm = get(kk("f34", "Minus"));
      const double efp = 0.5 * (get(kk("f34hi", "Plus")) - get(kk("f34lo", "Plus"))), efm = 0.5 * (get(kk("f34hi", "Minus")) - get(kk("f34lo", "Minus")));
      std::cout << Form("[CHG]   %-4s T+ %.3f T- %.3f  T+/T- %.2f+-%.2f (%+.1f sigma) | slope s+ %+.2f+-%.2f s- %+.2f+-%.2f (diff %+.1f sigma)"
                        " | F(30-40)/F(<30): + %.2f+-%.2f  - %.2f+-%.2f",
                        G.c_str(), Tp, Tm, rat, eRat, (rat - 1) / eRat, sp, esp, sm, esm, (sp - sm) / std::hypot(esp, esm), fp, efp, fm, efm);
      if (g < 0)
        std::cout << Form(" | SR multijet at m_T>=30: + %.0f%%  - %.0f%%", 100 * get(fam + "_fhi_incl_Plus"), 100 * get(fam + "_fhi_incl_Minus"));
      std::cout << "\n";
      gCsv << "chg," << lep << "," << fam << "," << G << "," << rat << "," << eRat << "," << sp << "," << esp << "," << sm << "," << esm
           << "," << fp << "," << efp << "," << fm << "," << efm << "\n";
    }
  }
}

void PrintEmu(const Run &mu, const Run &el, const std::string &outDir)
{
  std::cout << "\n=============== [EMU] e/mu of the stand-in r per rapidity bin (independent fits, stat. only) ===============\n";
  for (const std::string chg : {"Plus", "Minus"})
    for (const std::string sch : {"A", "A_groupshape", "A_inclshape", "A_flatF", "C", "C_groupshape", "C_2s", "C_trk", "R", "R_trk"})
    {
      const std::string key = sch + "_" + chg;
      if (!mu.r.count(key) || !el.r.count(key)) continue;
      const std::vector<double> &rm = mu.r.at(key), &em = mu.e.at(key), &re = el.r.at(key), &ee = el.e.at(key);
      std::cout << Form("[EMU] W%s %-13s:", chg == "Plus" ? "+" : "-", sch.c_str());
      double chi2 = 0;
      for (int iy = 0; iy < kNY; ++iy)
      {
        const double q = re[iy] / rm[iy], dq = q * std::sqrt(std::pow(ee[iy] / re[iy], 2) + std::pow(em[iy] / rm[iy], 2));
        std::cout << Form(" %.2f", q);
        chi2 += std::pow(re[iy] - rm[iy], 2) / (ee[iy] * ee[iy] + em[iy] * em[iy]);
        gCsv << "emu," << chg << "," << sch << ",y" << iy << "," << q << "," << dq << "\n";
      }
      std::cout << Form("   (edges y0 %.2f+-%.2f, y11 %.2f+-%.2f; chi2(r_mu = r_e) %.1f/12)\n", re[0] / rm[0],
                        re[0] / rm[0] * std::hypot(ee[0] / re[0], em[0] / rm[0]), re[11] / rm[11],
                        re[11] / rm[11] * std::hypot(ee[11] / re[11], em[11] / rm[11]), chi2);
    }
  // r per bin, mu vs e: the candidate scheme C (m_T-split SR, group anti-iso shape) and the reference R (lab y)
  for (const std::string chg : {"Plus", "Minus"})
  {
    std::vector<OverlaySeries> ser;
    struct S { const Run *r; std::string sch, lab; int col, mk; double ms, sh; };
    const std::vector<S> spec = {{&mu, "C_groupshape", "#mu, no m_{T} cut (C, group shape)", kAzure + 2, 20, 1.3, -0.5},
                                 {&el, "C_groupshape", "e, no m_{T} cut (C, group shape)", kRed + 1, 21, 1.5, -0.17},
                                 {&mu, "R", "#mu, m_{T}>40 (R)", kAzure - 9, 24, 1.3, 0.17},
                                 {&el, "R", "e, m_{T}>40 (R)", kRed - 7, 25, 1.5, 0.5}};
    for (const S &s : spec)
    {
      const std::string key = s.sch + "_" + chg;
      if (!s.r->r.count(key)) continue;
      TGraphErrors *g = new TGraphErrors(kNY);
      for (int iy = 0; iy < kNY; ++iy)
      {
        g->SetPoint(iy, 0.5 * (pOSkim::kYEdges[iy] + pOSkim::kYEdges[iy + 1]), s.r->r.at(key)[iy]);
        g->SetPointError(iy, 0.5 * (pOSkim::kYEdges[iy + 1] - pOSkim::kYEdges[iy]), s.r->e.at(key)[iy]);
      }
      OverlaySeries os;
      os.gTot = g; os.label = s.lab; os.color = s.col; os.marker = s.mk; os.markerSize = s.ms; os.xShift = s.sh;
      ser.push_back(os);
    }
    PlotStyle ps;
    ps.headerX = 0.47; ps.headerY = 0.89; ps.titleSize = 0.036; ps.subSize = 0.030; ps.headerDy = 0.042;
    ps.legX1 = 0.40; ps.legX2 = 0.92; ps.legY1 = 0.15; ps.legY2 = 0.33;
    SaveNiceGraph_Overlay(ser, outDir + "/r_nomt_W" + (chg == "Plus" ? "p" : "m"), "y = #minus#eta_{lab}", "r = data / MC",
                          Form("%s, stand-in fits (stat. only)", chg == "Plus" ? "W^{+}" : "W^{-}"),
                          "one fit per flavour, no #mu/e sharing", "e MC without Reco SF", {}, ps,
                          [](TCanvas *, TGraphErrors *g) { g->SetMinimum(0.4); g->SetMaximum(2.4); });
  }
}
} // namespace nomt

void qcd_nomt_study()
{
  using namespace nomt;
  gErrorIgnoreLevel = kWarning; // TH1::Chi2Test prints an Info per low-count bin
  const std::string outDir = "./plots/qcd_nomt_study";
  gSystem->mkdir(outDir.c_str(), kTRUE);
  gCsv.open((outDir + "/qcd_nomt_study.csv").c_str());
  gCsv << "# qcd_nomt_study.C -- rows: comp|fpt|ty|shape|mtfac|fisher|fit|emu (see the macro header)\n";
  Run runs[2];
  for (int e = 0; e < 2; ++e)
  {
    const LepCfg L = MakeCfg(e == 1);
    Files F;
    if (!OpenFiles(L, F)) { std::cerr << "[ERROR] inputs missing for " << L.lep << "\n"; return; }
    std::cout << "\n######## " << L.lep << ": iso-pass [" << L.pLo << "," << L.pHi << "), anti-iso [" << L.aLo << "," << L.aHi
              << ") (sub-windows at " << L.aMid << "), B/D at m_T < " << kMtB << ", fit window pT [24,100) ########\n";
    for (const std::string chg : {"Plus", "Minus"}) RunCharge(L, F, chg, runs[e], outDir);
  }
  PrintChg("mu", runs[0]);
  PrintChg("ele", runs[1]);
  PrintEmu(runs[0], runs[1], outDir);
  gCsv.close();
  std::cout << "\n[DONE] plots + csv in " << outDir << "\n";
}
