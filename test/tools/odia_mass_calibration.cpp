// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Does the per-run fragment mass calibration recover a KNOWN error, and does it
// refuse to invent one?
//
// With no arguments this runs entirely on synthesised residuals, so it is part
// of the ctest suite and depends on no data at all. Four properties are checked,
// and each of them is a way the estimator could be broken while still looking
// like it works:
//
//   1. A constant offset on a uniform background is recovered. This is the S08
//      case: a narrow true peak at about -10 ppm sitting on interference that is
//      flat across the +/-50 ppm search window.
//   2. A purely uniform sample is REJECTED. A uniform sample still has a
//      perfectly well-defined mode and width; they are just meaningless, and
//      applying them is worse than doing nothing.
//   3. An m/z-DEPENDENT offset is recovered as a trend, not as its average. A
//      mean-only estimator passes test 1 and fails here, which is the point --
//      a TOF's calibration error is characteristically a function of m/z.
//   4. A residual that is flat in m/z yields the CONSTANT model, not a spurious
//      slope. This is the guard against over-fitting: a line fitted to a flat
//      residual contorts at the ends of the m/z range, where a tryptic library
//      has fewest fragments.
//
//   4c. The two scale estimators are compared head to head on a planted sigma,
//      which is what makes "the reference's estimator does not work at this
//      purity" a measurement rather than an assertion.
//
//   4b. A trend that is REAL and statistically overwhelming, but which barely
//      moves the residual, still yields the constant. This is S08's own case,
//      and it is the one where a criterion based on significance -- or on how
//      many ppm the trend swings -- would ship a curve that buys nothing.
//
// Given arguments it instead runs the real thing on a real run, which is how the
// number in the commit message was measured. That path is deliberately NOT a
// ctest: it needs 14 GB of mzPeak and a minute of decode.
//
//   odia_mass_calibration <library.tsv> <run.mzpeak> [precursors] [cycles]

#include <odia/DIANNLibraryFile.h>
#include <odia/MassCalibration.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

namespace
{
  /// A sample shaped like the real thing: `n_true` matches drawn from a narrow
  /// Gaussian about a (possibly m/z-dependent) offset, plus `n_noise` matches
  /// drawn uniformly across the whole search window.
  ///
  /// The background is uniform on purpose. That is what an unrelated centroid
  /// inside a wide window actually looks like, and it is what makes a mean
  /// useless and a mode necessary.
  std::vector<ODIA::MassResidual> synthesise(std::mt19937& rng, std::size_t n_true,
                                             std::size_t n_noise, double offset_at_ref,
                                             double log_slope, double sigma,
                                             double search_ppm, bool with_control = true)
  {
    std::vector<ODIA::MassResidual> out;
    std::uniform_real_distribution<double> mz_of(300.0, 1500.0);
    std::uniform_real_distribution<double> rt_of(60.0, 1800.0);
    std::uniform_real_distribution<double> flat(-search_ppm, search_ppm);
    std::normal_distribution<double> jitter(0.0, sigma);
    std::uniform_real_distribution<double> bright(1e3, 1e5);

    for (std::size_t i = 0; i < n_true; ++i)
    {
      ODIA::MassResidual r;
      r.mz = static_cast<float>(mz_of(rng));
      r.rt = static_cast<float>(rt_of(rng));
      r.intensity = static_cast<float>(bright(rng));
      const double truth = offset_at_ref + log_slope * std::log(r.mz / 700.0);
      r.ppm = static_cast<float>(truth + jitter(rng));
      if (std::abs(r.ppm) > search_ppm) { continue; }   // the search would not have seen it
      out.push_back(r);
    }
    for (std::size_t i = 0; i < n_noise; ++i)
    {
      ODIA::MassResidual r;
      r.mz = static_cast<float>(mz_of(rng));
      r.rt = static_cast<float>(rt_of(rng));
      r.intensity = static_cast<float>(bright(rng));
      r.ppm = static_cast<float>(flat(rng));
      out.push_back(r);
    }
    if (with_control)
    {
      // The control cells hold no real ions, so their residuals are uniform by
      // construction. Supplying them exercises the null-referenced half of the
      // gate rather than leaving it dead in the test.
      for (std::size_t i = 0; i < n_noise; ++i)
      {
        ODIA::MassResidual r;
        r.mz = static_cast<float>(mz_of(rng));
        r.rt = static_cast<float>(rt_of(rng));
        r.intensity = static_cast<float>(bright(rng));
        r.ppm = static_cast<float>(flat(rng));
        r.decoy = true;
        out.push_back(r);
      }
    }
    return out;
  }

  int failures = 0;

  void check(bool ok, const std::string& what)
  {
    std::printf("  %s %s\n", ok ? "ok  " : "FAIL", what.c_str());
    if (!ok) { ++failures; }
  }

  int selftest()
  {
    ODIA::MassCalibration::Options opt;
    opt.search_ppm = 50.0;
    opt.min_residuals = 200;

    // ---- 1. a constant offset on a uniform background ---------------------
    // 5,000 true matches at -10.5 ppm against 10,000 uniform ones: 33% purity.
    // That is above the gate's operating point, which is worth knowing exactly.
    // For a narrow peak on a uniform background the peakedness statistic works
    // out at 1 + 5N/M, so the ported threshold of 3 is the assertion that at
    // least ~29% of the sample is real. Case 1b sits below it on purpose.
    {
      std::mt19937 rng(20260805);
      const auto r = synthesise(rng, 5000, 10000, -10.5, 0.0, 2.0, opt.search_ppm);
      const auto m = ODIA::MassCalibration::fit(r, opt);
      std::printf("1. constant offset, 33%% purity\n%s",
                  ODIA::MassCalibration::report(m, nullptr).c_str());
      check(m.fitted, "the gate passes on a peaked sample");
      check(std::abs(m.intercept_ppm - (-10.5)) < 0.5,
            "offset recovered within 0.5 ppm (got " + std::to_string(m.intercept_ppm) + ")");
      check(m.form == "constant", "the model chosen is constant");
      check(std::abs(m.sigma_after - 2.0) < 0.5,
            "sigma recovered within 0.5 ppm of the planted 2.0 (got " +
              std::to_string(m.sigma_after) + ")");
      check(m.window_ppm > 0.0 && m.window_ppm < 10.0,
            "the window is sized from sigma, not from the search width (got " +
              std::to_string(m.window_ppm) + ")");
      // A mean over the same sample is the thing the mode replaces; showing it
      // is what makes "use a robust estimator" a measurement rather than a
      // preference.
      double sum = 0.0; std::size_t n = 0;
      for (const auto& x : r) { if (!x.decoy) { sum += x.ppm; ++n; } }
      std::printf("     for comparison, the MEAN of the same sample is %.2f ppm\n",
                  n ? sum / static_cast<double>(n) : 0.0);
    }

    // ---- 1b. the same peak, too dilute to trust ---------------------------
    // 3,000 true against 15,000 uniform: 17% purity, peakedness ~2. The offset
    // is still there and the mode would still find it, which is the point --
    // the gate refuses on the strength of the evidence, not on whether an
    // answer can be computed.
    {
      std::mt19937 rng(20260809);
      const auto r = synthesise(rng, 3000, 15000, -10.5, 0.0, 2.0, opt.search_ppm);
      const auto m = ODIA::MassCalibration::fit(r, opt);
      std::printf("1b. the same offset at 17%% purity\n%s",
                  ODIA::MassCalibration::report(m, nullptr).c_str());
      check(!m.fitted, "the gate refuses a peak too dilute to be evidence");
      check(m.ppmAt(700.0) == 0.0, "and applies nothing");
    }

    // ---- 2. a purely uniform sample must be refused ------------------------
    {
      std::mt19937 rng(20260806);
      const auto r = synthesise(rng, 0, 16000, 0.0, 0.0, 2.0, opt.search_ppm);
      const auto m = ODIA::MassCalibration::fit(r, opt);
      std::printf("2. no signal at all\n%s", ODIA::MassCalibration::report(m, nullptr).c_str());
      check(!m.fitted, "the gate REJECTS a flat residual distribution");
      check(m.window_ppm <= 0.0, "no window is inferred from it");
      check(m.ppmAt(500.0) == 0.0 && m.ppmAt(1200.0) == 0.0,
            "the model is a no-op at every m/z when it did not fit");
    }

    // ---- 3. an m/z-dependent offset must come back as a trend --------------
    // -6 ppm at 700 Th with -8 ppm per e-fold in m/z: -1.5 ppm at 400 Th and
    // -11.5 ppm at 1400. A constant fitted to this would be ~-6 ppm and wrong by
    // about 5 ppm at both ends, which on a +/-10 ppm window is the difference
    // between catching a fragment and not -- and, decisively, it leaves a
    // visibly wider residual, which is what the model choice actually tests.
    {
      std::mt19937 rng(20260807);
      const auto r = synthesise(rng, 8000, 12000, -6.0, -8.0, 2.0, opt.search_ppm);
      ODIA::MassCalibration::Diagnostics d;
      const auto m = ODIA::MassCalibration::fit(r, opt, &d);
      std::printf("3. offset shaped in m/z\n%s", ODIA::MassCalibration::report(m, &d).c_str());
      check(m.fitted, "the gate passes");
      check(m.form == "log_mz", "the model chosen is the shaped one");
      check(std::abs(m.log_slope_ppm - (-8.0)) < 1.5,
            "slope recovered within 1.5 ppm per e-fold (got " +
              std::to_string(m.log_slope_ppm) + ")");
      const double truth400 = -6.0 - 8.0 * std::log(400.0 / 700.0);
      const double truth1400 = -6.0 - 8.0 * std::log(1400.0 / 700.0);
      check(std::abs(m.ppmAt(400.0) - truth400) < 1.0 &&
            std::abs(m.ppmAt(1400.0) - truth1400) < 1.0,
            "the correction is right at BOTH ends of the m/z range, not just on average");
      check(m.sigma_after < m.sigma_constant,
            "and it leaves a tighter residual than the constant (" +
              std::to_string(m.sigma_constant) + " -> " + std::to_string(m.sigma_after) + " ppm)");
    }

    // ---- 4. a flat residual must not acquire a slope -----------------------
    // Same size and purity as case 3, so the only difference is the truth.
    {
      std::mt19937 rng(20260808);
      const auto r = synthesise(rng, 8000, 12000, -10.0, 0.0, 2.0, opt.search_ppm);
      ODIA::MassCalibration::Diagnostics d;
      const auto m = ODIA::MassCalibration::fit(r, opt, &d);
      std::printf("4. flat in m/z, same sample size as 3\n%s",
                  ODIA::MassCalibration::report(m, &d).c_str());
      check(m.fitted, "the gate passes");
      check(m.form == "constant", "a flat residual yields the CONSTANT model");
      check(m.log_slope_ppm == 0.0, "and no shape is applied at all");
      check(std::abs(m.ppmAt(350.0) - (-10.0)) < 0.6 && std::abs(m.ppmAt(1450.0) - (-10.0)) < 0.6,
            "so the correction stays put at the sparse ends of the range");
    }

    // ---- 4b. a REAL trend that does not pay for itself ---------------------
    // This is S08 as measured: a genuine log-m/z term of -2.79 ppm per e-fold
    // sitting under a per-fragment scatter of ~7.8 ppm. The slope is hugely
    // significant -- with this many residuals it is tens of standard errors from
    // zero -- and modelling it still barely moves the residual. A criterion built
    // on significance, or on how many ppm the trend swings, would ship the curve
    // here. The criterion built on what the residual actually does, does not.
    {
      std::mt19937 rng(20260810);
      const auto r = synthesise(rng, 24000, 24000, -9.8, -2.79, 7.8, opt.search_ppm);
      ODIA::MassCalibration::Diagnostics d;
      const auto m = ODIA::MassCalibration::fit(r, opt, &d);
      std::printf("4b. a real but immaterial m/z trend (S08's own numbers)\n%s",
                  ODIA::MassCalibration::report(m, &d).c_str());
      check(m.fitted, "the gate passes");
      check(m.slope_t >= opt.min_slope_t,
            "the slope IS statistically significant (t=" + std::to_string(m.slope_t) + ")");
      check(m.form == "constant",
            "and the constant is chosen anyway, because the shape does not pay");
    }

    // ---- 4c. why the reference's scale estimator was replaced -------------
    // Not an opinion: the two estimators are run on the same sample with the
    // same known answer. `localScaleAboutMode` shrinks its neighbourhood by
    // 3 sigma per pass and stops when it stops shrinking, so when the first MAD
    // is already dominated by background it cannot contract at all and returns
    // the MIXTURE's scale. That is a window sized from the interference, which
    // is the exact failure the reference itself warns about from the other
    // direction. Pinned here so that if someone reverts to it, a test says why.
    {
      std::mt19937 rng(20260811);
      const auto r = synthesise(rng, 5000, 10000, 0.0, 0.0, 2.0, 50.0, false);
      std::vector<double> ppm;
      for (const auto& x : r) { ppm.push_back(x.ppm); }
      std::sort(ppm.begin(), ppm.end());
      std::vector<double> dev;
      for (double x : ppm) { dev.push_back(std::abs(x)); }
      std::sort(dev.begin(), dev.end());
      const double ported = ODIA::MassCalibration::localScaleAboutMode(ppm, 0.0, 50.0);
      const double used = ODIA::MassCalibration::backgroundCorrectedScale(dev, 50.0);
      std::printf("4c. the two scale estimators on a planted sigma of 2.0 ppm, 33%% purity\n"
                  "    localScaleAboutMode (reference): %.2f ppm\n"
                  "    backgroundCorrectedScale (used): %.2f ppm\n", ported, used);
      check(std::abs(used - 2.0) < 0.5, "the background-subtracting estimator recovers 2.0 ppm");
      check(ported > 2.0 * used,
            "while the ported one returns the mixture's scale, more than twice as wide");
    }

    // ---- 5. nothing at all ------------------------------------------------
    {
      const auto m = ODIA::MassCalibration::fit({}, opt);
      check(!m.fitted && m.window_ppm <= 0.0, "an empty sample is a no-op, not a crash");
    }

    std::printf("\n%s (%d failed check%s)\n", failures ? "FAILED" : "PASSED", failures,
                failures == 1 ? "" : "s");
    return failures ? 1 : 0;
  }
} // namespace

int main(int argc, char** argv)
{
  if (argc < 3) { return selftest(); }

  try
  {
    ODIA::Library library;
    ODIA::DIANNLibraryFile::load(argv[1], library);
    std::printf("library: %zu precursors, %zu transitions\n",
                library.precursorCount(), library.transitionCount());

    auto run = ODIA::openRun(argv[2]);
    std::printf("%s\n", run->describe().c_str());

    ODIA::MassCalibration::Options opt;
    if (argc > 3) { opt.max_precursors = std::strtoul(argv[3], nullptr, 10); }
    if (argc > 4) { opt.cycles = std::strtoul(argv[4], nullptr, 10); }
    if (const char* s = std::getenv("ODIA_MZCAL_SEARCH_PPM")) { opt.search_ppm = std::atof(s); }
    if (const char* s = std::getenv("ODIA_MZCAL_MIN_FRAGMENTS"))
    { opt.min_fragments_matched = std::strtoul(s, nullptr, 10); }
    if (const char* s = std::getenv("ODIA_MZCAL_IM_WINDOW")) { opt.im_window = std::atof(s); }

    ODIA::MassCalibration::Diagnostics d;
    const auto m = ODIA::MassCalibration::calibrate(library, *run, opt, &d);
    std::printf("%s", ODIA::MassCalibration::report(m, &d).c_str());

    if (const char* path = std::getenv("ODIA_MZCAL_DUMP"))
    {
      // The raw residuals, so the shape can be checked by something other than
      // the code that fitted it.
      const auto r = ODIA::MassCalibration::collect(library, *run, opt, nullptr);
      std::FILE* f = std::fopen(path, "wb");
      if (f)
      {
        std::fprintf(f, "mz\trt\tppm\tintensity\tdecoy\n");
        for (const auto& x : r)
        {
          std::fprintf(f, "%.5f\t%.2f\t%.4f\t%.6g\t%d\n", x.mz, x.rt, x.ppm, x.intensity,
                       x.decoy ? 1 : 0);
        }
        std::fclose(f);
        std::printf("wrote %s\n", path);
      }
    }
    return 0;
  }
  catch (const std::exception& e)
  {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
}
