// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// The CiRT seed's line fit, against the regime it was chosen for.
//
// doc/19 §3 budgets for ~50% WRONG apexes: a blind search over an uncalibrated
// run puts some standards on interference, and those land anywhere in the
// gradient. Least squares cannot survive that; Theil-Sen is the estimator for
// exactly this case, and this test is what says so with a number.
//
// It also pins the other half of the contract, which is the half that failed in
// production: the fit must REFUSE noise. The shipped seed accepted whatever the
// prefilter's contiguity statistic handed it, and that statistic passed 81% of
// decoys -- peptides that are not in the sample. A fit that reports a tight
// residual on noise is worse than no fit, because pass 1 then trusts it.

#include <odia/RtCalibration.h>

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

namespace
{
  // Deterministic: every number in this suite is compared against a previous
  // run's, so a seeded generator is not a nicety.
  std::mt19937 rng(20260818u);

  constexpr double TRUE_SLOPE = 1800.0;   // rt_norm [0,1] -> a 1800 s gradient
  constexpr double TRUE_INTERCEPT = 120.0;

  int check(bool ok, const char* what)
  {
    std::printf("  %-58s %s\n", what, ok ? "yes" : "NO");
    return ok ? 0 : 1;
  }
}

int main()
{
  int bad = 0;

  // ---- 1. half the anchors are wrong, and the line must still come back ----
  {
    std::uniform_real_distribution<double> x01(0.0, 1.0);
    std::normal_distribution<double> jitter(0.0, 8.0);      // real peak jitter
    std::uniform_real_distribution<double> anywhere(0.0, 1920.0);
    std::vector<std::pair<double, double>> anchors;
    for (int i = 0; i < 120; ++i)
    {
      const double x = x01(rng);
      // Every second anchor is an apex placed on interference: uncorrelated
      // with the library value, uniform over the gradient.
      const double y = (i % 2 == 0)
                         ? TRUE_SLOPE * x + TRUE_INTERCEPT + jitter(rng)
                         : anywhere(rng);
      anchors.push_back({x, y});
    }
    const auto line = ODIA::Calibration::fitRobustLine(anchors);
    std::printf("50%% wrong apexes: slope %.1f (true %.1f), intercept %.1f "
                "(true %.1f), p95 %.1f s, %zu inliers (%.0f%%)\n",
                line.slope, TRUE_SLOPE, line.intercept, TRUE_INTERCEPT,
                line.p95_residual, line.inliers, 100.0 * line.inlier_fraction);
    bad += check(line.ok, "fit succeeded");
    // 2% of the gradient. A least-squares line through this set lands hundreds
    // of seconds out, so the bar separates the two estimators rather than
    // merely being satisfiable.
    bad += check(std::fabs(line.slope - TRUE_SLOPE) < 0.02 * TRUE_SLOPE,
                 "slope within 2% of truth");
    bad += check(std::fabs(line.intercept - TRUE_INTERCEPT) < 40.0,
                 "intercept within 40 s of truth");
    bad += check(line.p95_residual < 40.0, "p95 residual under 40 s");
  }

  // ---- 2. pure noise must NOT produce a confident line --------------------
  // This is the decoy control's job in the seed, reproduced here without the
  // instrument: anchors with no relation between library RT and observed RT.
  {
    std::uniform_real_distribution<double> x01(0.0, 1.0);
    std::uniform_real_distribution<double> anywhere(0.0, 1920.0);
    std::vector<std::pair<double, double>> anchors;
    for (int i = 0; i < 120; ++i) { anchors.push_back({x01(rng), anywhere(rng)}); }
    const auto line = ODIA::Calibration::fitRobustLine(anchors);
    std::printf("pure noise: ok=%d, %zu inliers (%.0f%%), p95 %.1f s\n",
                int(line.ok), line.inliers, 100.0 * line.inlier_fraction,
                line.p95_residual);
    // The RESIDUAL cannot be the guard here, and assuming it could was the
    // first version of this test. RANSAC maximises consensus over every
    // candidate line, so noise still yields the luckiest handful of points --
    // measured, ~10 of 120 with a 32 s p95, which would sail through a
    // 10%-of-run gate. What must fail is the fit itself, on consensus size
    // against chance.
    bad += check(!line.ok, "noise refused outright, not fitted tightly");
  }

  // ---- 3. too few anchors is a refusal, not a guess -----------------------
  {
    std::vector<std::pair<double, double>> few{{0.1, 300.0}, {0.5, 1000.0},
                                               {0.9, 1700.0}};
    bad += check(!ODIA::Calibration::fitRobustLine(few).ok,
                 "3 anchors refused rather than fitted");
  }

  // ---- 4. a degenerate abscissa cannot define a slope ---------------------
  {
    std::vector<std::pair<double, double>> flat;
    for (int i = 0; i < 40; ++i) { flat.push_back({0.5, 100.0 + i}); }
    bad += check(!ODIA::Calibration::fitRobustLine(flat).ok,
                 "zero library-RT span refused");
  }

  if (bad) { std::fprintf(stderr, "%d robust-line check(s) failed\n", bad); return 1; }
  std::printf("robust line fit OK\n");
  return 0;
}
