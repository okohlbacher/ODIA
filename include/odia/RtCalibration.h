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
