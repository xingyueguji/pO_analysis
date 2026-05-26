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

/// Integrate a TH1 between [xmin, xmax] in x-axis units.
///
/// If `fullRange` is true, integrate bins [1, Nbins] (excluding under/overflow)
/// and ignore xmin/xmax. Otherwise integrate the bins containing xmin and xmax
/// inclusive, clamped to [1, Nbins].
///
/// Returns 0 if `h` is null.
inline double YieldInRange(const TH1 *h, double xmin, double xmax, bool fullRange = false)
{
    if (!h)
        return 0.0;

    if (fullRange || std::isnan(xmin) || std::isnan(xmax))
        return h->Integral(1, h->GetNbinsX());

    int b1 = h->GetXaxis()->FindBin(xmin);
    int b2 = h->GetXaxis()->FindBin(xmax);
    if (b1 < 1)
        b1 = 1;
    if (b2 > h->GetNbinsX())
        b2 = h->GetNbinsX();
    return h->Integral(b1, b2);
}

/// Uncertainty on the charge asymmetry A = (Np - Nm) / (Np + Nm),
/// assuming Poisson-independent Np, Nm.
///
///     sigma_A^2 = 4 * Np * Nm / (Np + Nm)^3
///
/// Returns 0 if Np + Nm <= 0.
inline double AsymErr(double Np, double Nm)
{
    const double S = Np + Nm;
    if (S <= 0.0)
        return 0.0;
    return std::sqrt(4.0 * Np * Nm / (S * S * S));
}

/// Uncertainty on the forward/backward ratio R = F / B,
/// assuming Poisson-independent F, B:
///
///     sigma_R = R * sqrt(1/F + 1/B)
///
/// Returns 0 if F or B is non-positive.
inline double RatioErr(double F, double B)
{
    if (F <= 0.0 || B <= 0.0)
        return 0.0;
    const double R = F / B;
    return R * std::sqrt(1.0 / F + 1.0 / B);
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
