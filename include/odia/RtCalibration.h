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
  /// Binned medians, then Cleveland's bisquare-reweighted LOESS, then PAVA
  /// isotonic regression. The last is not cosmetic -- elution order is
  /// physics, and a non-monotone RT map is always wrong however well it fits.
  ///
  /// @param anchors (library RT, observed RT) pairs.
  /// @param p95_resid if non-null, receives the 95th percentile residual.
  /// @param loess_span 0 lets the fit choose; the floor of ~15 neighbours
  ///        means the effective span grows as the anchor count shrinks.
  OpenMS::TransformationDescription fit(std::vector<std::pair<double, double>> anchors,
                                        double* p95_resid = nullptr,
                                        double loess_span = 0.0);

  /// The identity map, for when there are no anchors to fit.
  OpenMS::TransformationDescription identity();

} // namespace ODIA::Calibration
