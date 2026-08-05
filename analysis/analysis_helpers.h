// analysis_helpers.h
//
// Shared utilities for analysis/ stage macros (charge_asym.C, FBratio.C, ...).
// Pure-header, no link-time dependencies. Include with:
//
//     #include "analysis_helpers.h"
//
// All helpers live in the `pOAnalysis` namespace to avoid polluting the
// global / ROOT CLING scope.

#ifndef PO_ANALYSIS_ANALYSIS_HELPERS_H
#define PO_ANALYSIS_ANALYSIS_HELPERS_H

#include "TH1.h"
#include <cmath>

namespace pOAnalysis
{

/// An integrated yield together with its uncertainty.
///
/// `error` is the statistical uncertainty taken from the histogram's stored
/// sum-of-weights-squared (via TH1::IntegralAndError). For an *unweighted*
/// histogram this is exactly sqrt(N) (Poisson); for a weighted histogram it is
/// sqrt(sum w^2), which is the only correct thing — recomputing sqrt(N) from
/// the bin contents would be wrong as soon as event weights enter the skim.
struct Yield
{
    double value = 0.0;
    double error = 0.0;

    Yield() = default;
    Yield(double v, double e) : value(v), error(e) {}

    /// Combine two *independent* yields: values add, errors add in quadrature.
    /// Used e.g. for F = (W+ yield) + (W- yield) in the F/B ratio.
    Yield operator+(const Yield &o) const
    {
        return Yield(value + o.value,
                     std::sqrt(error * error + o.error * o.error));
    }
};

/// Integrate a TH1 between [xmin, xmax] in x-axis units, returning the yield
/// and its uncertainty (from the histogram's Sumw2 array).
///
/// If `fullRange` is true, integrate bins [1, Nbins] (excluding under/overflow)
/// and ignore xmin/xmax. Otherwise integrate the bins containing xmin and xmax
/// inclusive, clamped to [1, Nbins].
///
/// Returns a zero Yield if `h` is null.
inline Yield YieldInRange(const TH1 *h, double xmin, double xmax, bool fullRange = false)
{
    if (!h)
        return Yield();

    double err = 0.0;
    if (fullRange || std::isnan(xmin) || std::isnan(xmax))
    {
        const double val = h->IntegralAndError(1, h->GetNbinsX(), err);
        return Yield(val, err);
    }

    int b1 = h->GetXaxis()->FindBin(xmin);
    int b2 = h->GetXaxis()->FindBin(xmax);
    if (b1 < 1)
        b1 = 1;
    if (b2 > h->GetNbinsX())
        b2 = h->GetNbinsX();
    const double val = h->IntegralAndError(b1, b2, err);
    return Yield(val, err);
}

/// Uncertainty on the charge asymmetry A = (Np - Nm) / (Np + Nm),
/// by linear error propagation from Np, Nm, their errors, and (optionally)
/// their covariance:
///
///     dA/dNp =  2*Nm / S^2 ,  dA/dNm = -2*Np / S^2 ,  S = Np + Nm
///     sigma_A^2 = 4 / S^4 * ( Nm^2 sigma_Np^2 + Np^2 sigma_Nm^2
///                             - 2 Np Nm cov(Np, Nm) )
///
/// `cov` defaults to 0 (independent yields -- the pre-2026-08 behavior, exact
/// for raw-skim and legacy per-bin-fit inputs). The simfit grand fit supplies
/// a real cov(Np, Nm) (the shared r_Z and QCD correlate the per-charge fitted
/// yields) via the h_cov_yield matrix in comb_fitted_yields.root.
/// For unweighted independent yields this reduces to the old Poisson form
/// 4*Np*Nm/S^3. Returns 0 if Np + Nm <= 0 (or if rounding drives s2 < 0).
inline double AsymErr(const Yield &Np, const Yield &Nm, double cov = 0.0)
{
    const double S = Np.value + Nm.value;
    if (S <= 0.0)
        return 0.0;
    const double s2 = 4.0 / (S * S * S * S) *
                      (Nm.value * Nm.value * Np.error * Np.error +
                       Np.value * Np.value * Nm.error * Nm.error -
                       2.0 * Np.value * Nm.value * cov);
    return s2 > 0.0 ? std::sqrt(s2) : 0.0;
}

/// Uncertainty on the forward/backward ratio R = F / B,
/// by linear error propagation from F, B, their errors, and (optionally) their
/// covariance:
///
///     sigma_R^2 = R^2 * ( (sigma_F/F)^2 + (sigma_B/B)^2 - 2 cov(F,B)/(F*B) )
///
/// `cov` defaults to 0 (independent F, B -- the pre-2026-08 behavior); the
/// simfit grand fit supplies cov(F, B) via h_cov_yield_FB (see AsymErr).
/// For unweighted independent yields this reduces to the old Poisson form
/// R*sqrt(1/F + 1/B). Returns 0 if F or B is non-positive (or if rounding
/// drives s2 < 0).
inline double RatioErr(const Yield &F, const Yield &B, double cov = 0.0)
{
    if (F.value <= 0.0 || B.value <= 0.0)
        return 0.0;
    const double R = F.value / B.value;
    const double s2 = (F.error * F.error) / (F.value * F.value) +
                      (B.error * B.error) / (B.value * B.value) -
                      2.0 * cov / (F.value * B.value);
    return s2 > 0.0 ? R * std::sqrt(s2) : 0.0;
}

/// pO rapidity shift (lab -> nucleon-nucleon CM frame).
///
/// The pO beam is asymmetric (different per-nucleon beam energies because
/// the oxygen ion has Z/A != 1), so the lab and CM frames differ by a fixed
/// longitudinal boost:
///
///     y_CM = y_lab - kPORapidityShift
///
/// Hardcoded — update here (single place) if the beam configuration changes.
constexpr double kPORapidityShift = 0.3466;

} // namespace pOAnalysis

#endif // PO_ANALYSIS_ANALYSIS_HELPERS_H
