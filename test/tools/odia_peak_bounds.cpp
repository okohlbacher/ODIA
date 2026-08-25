// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

/// Does a peak group cover the peak, or the whole extraction window?
///
/// It covered the window. Measured on 400k production peak groups: the median
/// candidate spanned 129 of 130 cycles and 77.4% covered more than 80% of the
/// window, against a peak whose fitted FWHM on this data is 3.5 s -- about 2.5
/// cycles. The default co-elution picker descended while above 10% of the apex
/// with no rebound guard and no span bound, so on a noisy baseline the walk
/// never terminated. The amplitude picker had all three guards, inline and only
/// there.
///
/// That is invisible to every other test in this suite, because nothing else
/// asserts the GEOMETRY of a candidate -- only that the plumbing around it
/// holds. Three cases, each a failure mode a reviewer named:
///
///   1. a peak on a noisy baseline does not run to the window edge
///   2. a weak peak beside a stronger neighbour does not swallow the neighbour
///   3. a very sharp peak is widened to the minimum, because MS1_COELUTION
///      needs 5 cycles and the mass and mobility blocks need more than one --
///      below that width five sub-scores go NaN at once
#include <odia/PeakGroupScorer.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace
{
  int failures = 0;

  void check(bool ok, const char* what, long got, long lo, long hi)
  {
    std::printf("%-52s %4ld  expected %ld..%ld  %s\n", what, got, lo, hi,
                ok ? "ok" : "FAIL");
    if (!ok) { ++failures; }
  }

  /// A Gaussian of the measured width (sigma 1.07 cycles) on a flat baseline.
  std::vector<double> gaussian(std::size_t n, double centre, double sigma,
                               double height, double baseline)
  {
    std::vector<double> v(n, baseline);
    for (std::size_t i = 0; i < n; ++i)
    {
      const double d = (double(i) - centre) / sigma;
      v[i] += height * std::exp(-0.5 * d * d);
    }
    return v;
  }
}

int main()
{
  using Opt = ODIA::PeakGroupScorer::Options;
  const Opt d;
  std::printf("defaults: min_cycles %zu  max_half %zu  boundary_fraction %.2f\n\n",
              d.peak_min_cycles, d.peak_max_half_cycles, d.boundary_fraction);

  // 1. A peak on a baseline that never falls below a tenth of the apex.
  //    This is the production case: baseline 12% of apex, so the old rule's
  //    floor is never reached and it walks to both edges.
  {
    auto t = gaussian(130, 65.0, 1.07, 100.0, 12.0);
    const auto b = ODIA::PeakGroupScorer::peakBoundsForTest(
      t, 65, 65, d.boundary_fraction, d.peak_min_cycles, d.peak_max_half_cycles);
    const long w = long(b.second - b.first + 1);
    check(w <= 41, "peak on a 12%-of-apex baseline: width", w, 1, 41);
  }

  // 2. A weak peak with a stronger neighbour 12 cycles away. The boundary must
  //    not cross the valley; if it does, every sub-score integrates both.
  {
    auto t = gaussian(130, 60.0, 1.07, 40.0, 2.0);
    const auto big = gaussian(130, 72.0, 1.07, 200.0, 0.0);
    for (std::size_t i = 0; i < t.size(); ++i) { t[i] += big[i]; }
    const auto b = ODIA::PeakGroupScorer::peakBoundsForTest(
      t, 60, 60, d.boundary_fraction, d.peak_min_cycles, d.peak_max_half_cycles);
    check(long(b.second) <= 66, "weak peak beside a 5x neighbour: right edge",
          long(b.second), 60, 66);
  }

  // 3. A sharp isolated peak must still be widened to the minimum, or
  //    MS1_COELUTION, MASS_ACCURACY, MASS_SPREAD, IM_DELTA and IM_SPREAD all
  //    return NaN together.
  {
    auto t = gaussian(130, 65.0, 0.6, 100.0, 0.01);
    const auto b = ODIA::PeakGroupScorer::peakBoundsForTest(
      t, 65, 65, d.boundary_fraction, d.peak_min_cycles, d.peak_max_half_cycles);
    const long w = long(b.second - b.first + 1);
    check(w >= long(d.peak_min_cycles), "sharp peak widened to the minimum", w,
          long(d.peak_min_cycles), 130);
  }

  std::printf("\n%s\n", failures ? "FAILURES" : "all peak-boundary cases pass");
  return failures ? 1 : 0;
}
