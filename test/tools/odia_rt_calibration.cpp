// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Does the adopted RT calibration recover a known warp from dirty anchors?
//
// The properties that matter are not "it fits". They are:
//   * it recovers a NON-LINEAR warp, since that is what a predicted library
//     against a real gradient actually looks like;
//   * it survives gross outliers, because the anchors are first-pass peak
//     groups found in windows narrower than the RT error, so residuals of
//     thousands of seconds are guaranteed rather than possible;
//   * the map is MONOTONE, because elution order is physics and a map that
//     reorders the run is wrong however well it fits.
//
// A test that only checked the residual would pass on a fit that silently
// reordered the gradient.
#include <odia/RtCalibration.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

int main()
{
  // The truth: a smooth, monotone, distinctly non-linear compression of the
  // library scale onto a 300..3300 s gradient.
  const auto truth = [](double x) { return 300.0 + 2600.0 * (x + 0.35 * x * x) / 1.35; };

  std::mt19937 rng(20260804);
  std::normal_distribution<double> jitter(0.0, 12.0);
  std::uniform_real_distribution<double> anywhere(300.0, 3300.0);

  std::vector<std::pair<double, double>> anchors;
  int planted_outliers = 0;
  for (int i = 0; i < 400; ++i)
  {
    const double x = static_cast<double>(i) / 399.0;
    if (i % 12 == 0)
    {
      anchors.emplace_back(x, anywhere(rng));   // a gross outlier, as pass 1 produces
      ++planted_outliers;
    }
    else
    {
      anchors.emplace_back(x, truth(x) + jitter(rng));
    }
  }
  std::printf("%zu anchors, %d of them gross outliers (%.0f%%)\n",
              anchors.size(), planted_outliers,
              100.0 * planted_outliers / static_cast<double>(anchors.size()));

  double p95 = 0.0;
  const auto trafo = ODIA::Calibration::fit(anchors, &p95, 0.0);
  std::printf("p95 anchor residual: %.2f s\n", p95);

  double worst = 0.0, sum = 0.0;
  double previous = -1e300;
  int non_monotone = 0;
  const int probes = 200;
  for (int k = 0; k <= probes; ++k)
  {
    const double x = static_cast<double>(k) / probes;
    const double got = trafo.apply(x);
    const double want = truth(x);
    const double err = std::abs(got - want);
    worst = std::max(worst, err);
    sum += err;
    if (got < previous) { ++non_monotone; }
    previous = got;
  }
  const double mean = sum / (probes + 1);
  std::printf("recovered the warp: mean |error| %.2f s, worst %.2f s\n", mean, worst);
  std::printf("monotone: %s\n", non_monotone ? "NO" : "yes");

  if (non_monotone)
  {
    std::fprintf(stderr, "the fitted map is not monotone at %d of %d probes; it "
                         "reorders the gradient\n", non_monotone, probes + 1);
    return 1;
  }
  // Generous on purpose: the point is that a third of the range is not lost to
  // the outliers, not that the fit is exact. The gradient spans 3000 s, so
  // 60 s is 2% of it -- and a fit that ignored the outliers entirely would
  // land near the 12 s jitter, while one that took them seriously would be
  // hundreds of seconds out.
  if (!(mean < 60.0))
  {
    std::fprintf(stderr, "mean error %.2f s: the outliers were not resisted\n", mean);
    return 1;
  }
  std::printf("adopted RT calibration OK\n");
  return 0;
}
