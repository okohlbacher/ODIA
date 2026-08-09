// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

/// Does the fragment window sized from identifications measure the right thing,
/// and refuse to answer when it cannot?
///
/// Two failures this suite exists to prevent, both of which have already
/// happened once on real data:
///
///  1. THE WRONG UNIT. `PeakGroup::mass_ppm` is a median over (fragment x
///     cycle) cells, so its scatter across groups is sigma/sqrt(N_eff) and is
///     several times narrower than the fragment scatter a window must admit.
///     Sizing a window from it once turned a 1.6 ppm quantity into an 18 ppm
///     claim, and the same confusion in the other direction would produce an
///     absurdly narrow window that looked like a beautifully calibrated
///     instrument. The estimator must track the PER-FRAGMENT sigma and be
///     insensitive to how many groups there are.
///
///  2. CENSORING. A deviation can only be observed inside the window it was
///     extracted through, so the observed distribution is truncated and every
///     estimate from it is biased narrow. Narrowing on that estimate and
///     re-measuring spirals inward. A measurement taken through a window it is
///     pressed against must be refused, not scaled.

#include <odia/MassWidth.h>

#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

namespace
{
  int failures = 0;
  void check(bool ok, const std::string& what)
  {
    std::printf("  [%s] %s\n", ok ? "ok" : "FAIL", what.c_str());
    if (!ok) { ++failures; }
  }

  /// `n` accepted precursors whose per-fragment deviations are N(centre, sigma),
  /// each group summarising `frags` fragments the way the scorer does.
  ODIA::PeakGroupScorer::Result plant(std::size_t n, double centre, double sigma,
                                      std::size_t frags, unsigned seed,
                                      double decoy_centre = 0.0)
  {
    std::mt19937 rng(seed);
    std::normal_distribution<double> draw(centre, sigma);
    ODIA::PeakGroupScorer::Result r;
    r.fdr_valid = true;
    for (std::size_t i = 0; i < n; ++i)
    {
      std::vector<double> per_fragment;
      for (std::size_t k = 0; k < frags; ++k) { per_fragment.push_back(draw(rng)); }
      std::sort(per_fragment.begin(), per_fragment.end());
      const double median = per_fragment[per_fragment.size() / 2];
      std::vector<double> abs_dev;
      for (const double v : per_fragment) { abs_dev.push_back(std::abs(v - median)); }
      std::sort(abs_dev.begin(), abs_dev.end());

      ODIA::PeakGroupScorer::PeakGroup g;
      g.precursor = static_cast<std::uint32_t>(i);
      g.qvalue = 0.001;
      g.dscore = 1.0;
      g.decoy = false;
      g.mass_ppm = static_cast<float>(median);
      g.mass_ppm_n = static_cast<std::uint16_t>(frags * 8);
      g.mass_ppm_spread = static_cast<float>(1.4826 * abs_dev[abs_dev.size() / 2]);
      g.mass_ppm_frags = static_cast<std::uint8_t>(frags);
      r.groups.push_back(g);

      // A decoy at the same q must never contribute; it is planted off-centre
      // so that it would visibly move the answer if it did.
      ODIA::PeakGroupScorer::PeakGroup d = g;
      d.precursor = static_cast<std::uint32_t>(n + i);
      d.decoy = true;
      d.mass_ppm = static_cast<float>(decoy_centre + 40.0);
      d.mass_ppm_spread = 30.0f;
      r.groups.push_back(d);
    }
    return r;
  }
}

int main()
{
  std::printf("mass width from identifications\n");

  // Sigma is the PER-FRAGMENT scatter, and does not shrink as groups are added.
  {
    const auto small = ODIA::MassWidth::measure(plant(300, -8.0, 4.0, 6, 1), 0.01,
                                                50.0, 3.0, 100);
    const auto large = ODIA::MassWidth::measure(plant(3000, -8.0, 4.0, 6, 2), 0.01,
                                                50.0, 3.0, 100);
    check(small.valid && large.valid, "both estimates valid");
    check(std::abs(small.sigma_ppm - 4.0) < 1.2, "sigma recovers the planted 4 ppm");
    check(std::abs(large.sigma_ppm - small.sigma_ppm) < 1.0,
          "sigma does not shrink with 10x the groups (it is not a standard error)");
    check(std::abs(small.centre_ppm + 8.0) < 0.8, "centre recovers the planted -8 ppm");
    check(std::abs(small.width_ppm - 3.0 * small.sigma_ppm) < 1e-9,
          "width is sigmas x sigma");
    check(!small.censored, "4 ppm scatter measured through 50 ppm is not censored");
  }

  // The censoring guard: the same scatter through a window it fills.
  {
    const auto e = ODIA::MassWidth::measure(plant(500, 0.0, 12.0, 6, 3), 0.01,
                                            15.0, 3.0, 100);
    check(e.valid, "estimate is produced");
    check(e.censored, "12 ppm scatter measured through a 15 ppm window is censored");
  }
  {
    const auto e = ODIA::MassWidth::measure(plant(500, 0.0, 2.0, 6, 4), 0.01,
                                            50.0, 3.0, 100);
    check(!e.censored, "2 ppm scatter through 50 ppm is not censored");
  }

  // Refusals.
  {
    const auto few = ODIA::MassWidth::measure(plant(50, 0.0, 4.0, 6, 5), 0.01,
                                              50.0, 3.0, 200);
    check(!few.valid, "50 accepted groups against a floor of 200 is refused");

    auto invalid = plant(500, 0.0, 4.0, 6, 6);
    invalid.fdr_valid = false;
    check(!ODIA::MassWidth::measure(invalid, 0.01, 50.0, 3.0, 100).valid,
          "a pass whose FDR is not valid supplies nothing");

    // Two fragments is below `min_fragments`, so no group qualifies.
    auto thin = plant(500, 0.0, 4.0, 6, 7);
    for (auto& g : thin.groups) { g.mass_ppm_frags = 2; }
    check(!ODIA::MassWidth::measure(thin, 0.01, 50.0, 3.0, 100).valid,
          "groups with fewer than 3 fragments do not contribute");
  }

  // Decoys and rejected targets are excluded, so a wildly off-centre decoy
  // population cannot move the centre.
  {
    const auto e = ODIA::MassWidth::measure(plant(500, -8.0, 4.0, 6, 8), 0.01,
                                            50.0, 3.0, 100);
    check(std::abs(e.centre_ppm + 8.0) < 0.8,
          "decoys at +40 ppm do not move the centre");
    check(e.groups == 500, "exactly the accepted targets contribute");
  }

  // One row per precursor: extra candidates of the same precursor must not
  // weight it more heavily.
  {
    auto r = plant(400, -8.0, 4.0, 6, 9);
    const std::size_t original = r.groups.size();
    for (std::size_t i = 0; i < original; ++i)
    {
      if (r.groups[i].decoy) { continue; }
      auto extra = r.groups[i];
      extra.dscore = 0.1;              // a worse candidate of the same precursor
      extra.mass_ppm = 25.0f;
      extra.mass_ppm_spread = 25.0f;
      r.groups.push_back(extra);
    }
    const auto e = ODIA::MassWidth::measure(r, 0.01, 50.0, 3.0, 100);
    check(e.groups == 400, "each precursor contributes exactly one row");
    check(std::abs(e.centre_ppm + 8.0) < 0.8,
          "the worse candidate of each precursor is not the one used");
  }

  std::printf("%s\n", failures ? "FAILED" : "all good");
  return failures ? 1 : 0;
}
