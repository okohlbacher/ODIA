// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>

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
  OpenMS::TransformationDescription fit(std::vector<std::pair<double, double>> anchors,
                                        double* p95_resid = nullptr,
                                        double loess_span = 0.0);

  /// The identity map, for when there are no anchors to fit.
  OpenMS::TransformationDescription identity();

} // namespace ODIA::Calibration
