// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>

#include <string>
#include <utility>
#include <vector>

namespace ODIA::Calibration
{

  /// Fit the map from library retention time onto this run's, from anchors.
  ///
  /// Adopted from okohlbacher/OpenDIAlyzer, tag odia-v0.3.0 (0ed5fb2), where it
  /// lived as private statics of the analyzer tool. Lifted out unchanged in
  /// substance so ODIA can use it without the rest of that tool.
  ///
  /// Robust by construction, because the anchors are dirty by construction:
  /// they are first-pass peak groups found in windows narrower than the RT
  /// error, so residuals of thousands of seconds are guaranteed in the set.
  /// Binned medians, optionally Cleveland's bisquare-reweighted LOESS, then
  /// PAVA isotonic regression. The last is not cosmetic -- elution order is
  /// physics, and a non-monotone RT map is always wrong however well it fits.
  ///
  /// "Optionally" is load-bearing and was not said here before: LOESS is gated
  /// on `loess_span > 0`, every production call site passed 0, so the middle
  /// step has never run. `-rt_loess_span` exists to measure what it is worth.
  ///
  /// @param anchors (library RT, observed RT) pairs.
  /// @param p95_resid if non-null, receives the 95th percentile residual.
  /// @param loess_span 0 uses BINNED MEDIANS ALONE -- it does not "let the fit
  ///        choose", which is what this comment claimed while every production
  ///        call site passed 0. LOESS runs only for span > 0 and >= 50 anchors.
  /// @param interpolation how the knots are joined: "akima" (default,
  ///        nonlinear and outlier-resistant), "cspline" (nonlinear, rings
  ///        around outliers) or "linear" (what this used to do -- a
  ///        retention-time map has no reason to be piecewise-linear, and the
  ///        derivative jumps at every knot). OpenSWATH aligns iRT with LOWESS
  ///        and DIA-NN fits a nonlinear monotone regression; neither joins
  ///        knots with line segments.
  OpenMS::TransformationDescription fit(std::vector<std::pair<double, double>> anchors,
                                        double* p95_resid = nullptr,
                                        double loess_span = 0.0,
                                        const std::string& interpolation = "akima");

  /// A robust straight line through dirty anchors, by Theil-Sen plus one
  /// inlier-restricted least-squares refit.
  ///
  /// LINEAR, and deliberately so. doc/20 measured both options from a blind
  /// CiRT search and a monotone PCHIP LOST 5-7 percentage points at +/-60 s on
  /// both instruments (Astral 65.5% -> 58.9%, IH1 81.2% -> 76.7%): ~30 anchors
  /// cannot constrain a curve, so the fit chases anchor noise. The
  /// nonlinearity is real and worth +5.1 pp (Astral) to +12.6 pp (IH1) -- but
  /// only from THOUSANDS of identifications, which is fine-tuning's job, not
  /// the seed's. Fitting the curve here is the measured mistake.
  ///
  /// RANSAC rather than Theil-Sen, decided by measurement rather than by the
  /// slash in doc/19 §3's "RANSAC / Theil-Sen". A blind search over an
  /// uncalibrated run places apexes on interference, and doc/20 records how
  /// many survive: 24 inliers on Astral, 47 on IH1, from ~149 precursors. That
  /// is 68-84% wrong, well past Theil-Sen's 29.3% breakdown -- and past it
  /// Theil-Sen degrades QUIETLY (3.7% slope error, 221 s p95 on the synthetic
  /// 50% case) rather than failing, which is the worst behaviour for a number
  /// the seed then gates on.
  struct Line
  {
    double slope = 0.0;
    double intercept = 0.0;
    /// Residuals of the INLIERS in the units of the anchors' second element.
    double p95_residual = 0.0;
    double median_abs_residual = 0.0;
    std::size_t inliers = 0;
    /// Consensus size as a fraction of the anchors offered. A tight residual
    /// over a small consensus is what coincidence looks like.
    double inlier_fraction = 0.0;
    bool ok = false;
  };

  /// @param anchors (library RT, observed RT) pairs.
  /// @param max_pairs cap on the Theil-Sen pair count; above it the pairs are
  ///        strided deterministically rather than sampled at random, so the
  ///        answer does not move between runs.
  /// @param inlier_tolerance consensus half-width in the anchors' y units;
  ///        0 uses 2% of the observed y span.
  Line fitRobustLine(const std::vector<std::pair<double, double>>& anchors,
                     double inlier_tolerance = 0.0,
                     std::size_t max_pairs = 4000000);

  /// The identity map, for when there are no anchors to fit.
  OpenMS::TransformationDescription identity();

  /// Invert a monotone retention-time map: given RT, return the iRT that maps
  /// to it.
  ///
  /// NOT `TransformationDescription::invert()`. That swaps the control points
  /// and refits, so for an akima or cspline model the "inverse" is a NEW
  /// interpolation through the swapped knots -- close to the true inverse where
  /// the knots are dense and free to diverge from it where they are not, with
  /// no guarantee that `invert(apply(x)) == x` anywhere. Fitting a curve to
  /// approximate a function we already have exactly is the wrong move.
  ///
  /// This does the direct thing instead: the forward map is monotone by
  /// construction (PAVA), so bisect on it. `apply` is a few interpolation
  /// lookups, ~60 iterations pins the answer to machine precision over any
  /// plausible gradient, and the result is the inverse of the map actually in
  /// use rather than of a refit approximation to it.
  ///
  /// @param lo,hi the iRT bracket to search. Outside it the answer is clamped,
  ///        which is correct: a retention time beyond the map's range has no
  ///        iRT preimage and inventing one by extrapolation is how a refinement
  ///        puts windows off the gradient.
  double invertAt(const OpenMS::TransformationDescription& map, double rt,
                  double lo, double hi);

} // namespace ODIA::Calibration
