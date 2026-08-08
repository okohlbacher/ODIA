// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

/// Does the RT-blocked mass recalibration recover a drift that is actually there,
/// and refuse one that is not?
///
/// The existing MassCalibration fits ONE global correction and cannot express
/// drift. A mass spectrometer drifts over an hour of gradient -- source
/// contamination, temperature, space charge -- so a single number is fitted to
/// a moving target and is wrong at both ends by construction.
///
/// The failure this suite must prevent is the one that has already happened
/// twice on real data: a statistic computed over a large search space that looks
/// like a measurement and is a property of the search. An unanchored m/z
/// residual read -8.44 ppm on S08 while the same probe 300 s away from the
/// peptide read -4.98 -- so most of it was the +/-50 ppm window, not the
/// instrument. Pinning the contaminated figure cost 90 identifications.

#include <odia/MassRecalibration.h>

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

  /// Residuals with a planted drift: ppm = a(rt) + b * log(mz/500).
  std::vector<ODIA::MassResidual> plant(std::size_t n, double a0, double a1,
                                        double slope, double noise_sd, unsigned seed)
  {
    std::mt19937 rng(seed);
    std::normal_distribution<double> noise(0.0, noise_sd);
    std::uniform_real_distribution<double> rt(0.0, 2400.0);
    std::uniform_real_distribution<double> mz(200.0, 1800.0);
    std::vector<ODIA::MassResidual> out;
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
    {
      ODIA::MassResidual r;
      r.rt = static_cast<float>(rt(rng));
      r.mz = static_cast<float>(mz(rng));
      const double drift = a0 + a1 * (r.rt / 2400.0);         // linear in time
      r.ppm = static_cast<float>(drift + slope * std::log(r.mz / 500.0) + noise(rng));
      r.intensity = 1000.0f;
      r.decoy = false;
      out.push_back(r);
    }
    return out;
  }
} // namespace

int main()
{
  ODIA::MassRecalibration::Options opt;

  // ---- 1. a planted drift is recovered, at both ends -----------------------
  {
    const auto res = plant(20000, -8.0, +12.0, 0.0, 1.5, 0xC0FFEEu);   // -8 -> +4 ppm
    const auto m = ODIA::MassRecalibration::fit(res, opt);
    std::printf("1. planted drift -8 -> +4 ppm over the gradient\n%s", m.describe().c_str());
    check(m.fitted(), "a drift over 20,000 anchors is fitted");
    const double early = m.ppmAt(500.0, 100.0), late = m.ppmAt(500.0, 2300.0);
    check(std::abs(early - (-7.5)) < 1.5, "the EARLY block recovers about -8 ppm");
    check(std::abs(late - (+3.5)) < 1.5, "the LATE block recovers about +4 ppm");
    check(late - early > 8.0,
          "and the model expresses the DRIFT -- a single global constant could not");
  }

  // ---- 2. the m/z slope is recovered independently of the drift ------------
  {
    const auto res = plant(20000, 0.0, 0.0, 6.0, 1.5, 0xBEEFu);
    const auto m = ODIA::MassRecalibration::fit(res, opt);
    std::printf("\n2. planted slope +6 ppm per e-fold, no drift\n");
    const double lo = m.ppmAt(200.0, 1200.0), hi = m.ppmAt(1800.0, 1200.0);
    const double got = (hi - lo) / std::log(1800.0 / 200.0);
    std::printf("   recovered slope %.2f ppm per e-fold\n", got);
    check(std::abs(got - 6.0) < 1.5, "the m/z slope comes back at about +6");
    check(std::abs(m.ppmAt(500.0, 1200.0)) < 1.5,
          "and the intercept stays near zero -- the two are not confounded");
  }

  // ---- 3. a FLAT run is not given a correction it does not need ------------
  {
    const auto res = plant(20000, 0.0, 0.0, 0.0, 1.5, 0x1234u);
    const auto m = ODIA::MassRecalibration::fit(res, opt);
    std::printf("\n3. no drift, no slope, pure noise\n");
    check(m.maxAbsCorrection() < 1.5,
          "a flat run yields a near-zero correction everywhere, not a fitted wobble");
  }

  // ---- 4. controls are ignored --------------------------------------------
  {
    auto res = plant(10000, -6.0, 0.0, 0.0, 1.0, 0x99u);
    auto bad = plant(10000, +30.0, 0.0, 0.0, 1.0, 0xAAu);
    for (auto& r : bad) { r.decoy = true; }        // an m/z-shifted control
    res.insert(res.end(), bad.begin(), bad.end());
    const auto m = ODIA::MassRecalibration::fit(res, opt);
    std::printf("\n4. half the input is a shifted control at +30 ppm\n");
    check(std::abs(m.ppmAt(500.0, 1200.0) - (-6.0)) < 1.5,
          "controls are excluded -- a null cannot inform a correction");
  }

  // ---- 5. too little to fit is refused ------------------------------------
  {
    const auto res = plant(20, -8.0, 0.0, 0.0, 1.0, 0x55u);
    const auto m = ODIA::MassRecalibration::fit(res, opt);
    std::printf("\n5. twenty anchors\n");
    check(!m.fitted() && m.ppmAt(500.0, 1200.0) == 0.0,
          "twenty anchors is refused outright rather than fitted badly");
  }

  // ---- 6. the correction is bounded ---------------------------------------
  {
    const auto res = plant(20000, -500.0, 0.0, 0.0, 1.0, 0x77u);
    const auto m = ODIA::MassRecalibration::fit(res, opt);
    std::printf("\n6. an absurd -500 ppm input\n");
    check(m.maxAbsCorrection() <= opt.max_correction_ppm + 1e-6,
          "the correction is clamped -- a fit that wants 500 ppm is not a calibration");
  }

  std::printf("\n%s (%d failed check%s)\n", failures ? "FAILED" : "PASSED", failures,
              failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
