// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// The held-out evaluator for mass models, and the invariants the anchor
// harvest has to hold.
//
// `MassCalibration::fit` already reports a systematic residual, but against its
// own sample -- that measures how well a model describes the data it was fitted
// from. `systematicResidualPpm` scores a model on residuals it never saw, which
// is a different question and the one doc/15 step 1 turns on. These checks
// exist because the difference between the two is invisible on a good fit and
// decisive on an overfitted one.

#include <odia/MassCalibration.h>

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <random>
#include <vector>

namespace
{
  int failures = 0;

  void check(bool ok, const char* what)
  {
    std::printf("%s  %s\n", ok ? "ok  " : "FAIL", what);
    if (!ok) { ++failures; }
  }

  void checkNear(double got, double want, double tol, const char* what)
  {
    const bool ok = std::abs(got - want) <= tol;
    std::printf("%s  %s (got %.3f, want %.3f +/- %.3f)\n",
                ok ? "ok  " : "FAIL", what, got, want, tol);
    if (!ok) { ++failures; }
  }

  /// Residuals whose systematic error is a known function of m/z, plus noise.
  /// `n` is deliberately large: the evaluator bins 8-deep at 40 minimum, so a
  /// small sample returns 0.0 and would make every check pass vacuously.
  std::vector<ODIA::MassResidual> synth(std::mt19937& rng, std::size_t n,
                                        double intercept, double log_slope,
                                        double noise_ppm, bool decoy = false)
  {
    std::uniform_real_distribution<double> mz(250.0, 1500.0);
    std::normal_distribution<double> eps(0.0, noise_ppm);
    std::vector<ODIA::MassResidual> out;
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
    {
      ODIA::MassResidual r;
      r.mz = static_cast<float>(mz(rng));
      r.ppm = static_cast<float>(intercept +
                                 log_slope * std::log(r.mz / 700.0) + eps(rng));
      r.rt = 1000.0f;
      r.intensity = 1000.0f;
      r.im = 0.9f;
      r.decoy = decoy;
      out.push_back(r);
    }
    return out;
  }
}

int main()
{
  std::mt19937 rng(20260810);

  // --- the evaluator measures what it claims to measure -------------------
  {
    // A pure -9 ppm offset with a +3 ppm/e-fold slope, the S08 shape.
    const auto rs = synth(rng, 6000, -9.0, 3.0, 4.0);

    const double none = ODIA::MassCalibration::systematicResidualPpm(
      rs, [](double) { return 0.0; });
    const double constant = ODIA::MassCalibration::systematicResidualPpm(
      rs, [](double) { return -9.0; });
    const double truth = ODIA::MassCalibration::systematicResidualPpm(
      rs, [](double mz) { return -9.0 + 3.0 * std::log(mz / 700.0); });

    check(none > constant, "uncorrected leaves more systematic error than a constant");
    check(constant > truth, "a constant leaves more than the true shape");
    checkNear(truth, 0.0, 0.35, "the true correction leaves ~nothing");
  }

  // --- it is not fooled by a correction that shifts everything equally ----
  {
    // The modes agree with EACH OTHER but sit off zero. A statistic taken about
    // a refitted constant would call this perfect; this one must not, because a
    // uniform offset is exactly the error a window mis-centres on.
    const auto rs = synth(rng, 6000, -9.0, 0.0, 4.0);
    const double biased = ODIA::MassCalibration::systematicResidualPpm(
      rs, [](double) { return 0.0; });
    check(biased > 5.0, "a uniform -9 ppm offset is reported as systematic error");
  }

  // --- a model fitted on one sample, scored on another --------------------
  {
    const auto train = synth(rng, 6000, -9.0, 3.0, 4.0);
    const auto test  = synth(rng, 6000, -9.0, 3.0, 4.0);

    const ODIA::MassCalibration::Options opt;
    const auto m = ODIA::MassCalibration::fit(train, opt, nullptr);
    check(m.fitted, "a clean shaped sample is fitted");

    const double held = ODIA::MassCalibration::systematicResidualPpm(
      test, [&](double mz) { return m.ppmAt(mz); });
    const double held_none = ODIA::MassCalibration::systematicResidualPpm(
      test, [](double) { return 0.0; });
    check(held < held_none, "the fitted model beats no correction on HELD-OUT data");
    checkNear(held, 0.0, 0.60, "and leaves little on data it never saw");
  }

  // --- a model fitted to the WRONG shape does not transfer ----------------
  {
    // Train where the slope runs one way, test where it runs the other. An
    // in-sample statistic would look fine on the training set; the held-out one
    // has to notice. This is the overfitting case the whole helper exists for.
    const auto train = synth(rng, 6000, -9.0, +6.0, 3.0);
    const auto test  = synth(rng, 6000, -9.0, -6.0, 3.0);

    const ODIA::MassCalibration::Options opt;
    const auto m = ODIA::MassCalibration::fit(train, opt, nullptr);
    if (m.fitted)
    {
      const double held = ODIA::MassCalibration::systematicResidualPpm(
        test, [&](double mz) { return m.ppmAt(mz); });
      const double held_const = ODIA::MassCalibration::systematicResidualPpm(
        test, [&](double) { return m.intercept_ppm; });
      check(held > held_const,
            "a shape fitted to the wrong trend is WORSE than its own constant held out");
    }
    else { check(true, "wrong-shape model refused outright, which is also correct"); }
  }

  // --- control cells are excluded -----------------------------------------
  {
    // Decoy residuals are uniform across the window by construction. If they
    // reached the bins they would drag every mode towards the centre and make a
    // bad correction look good.
    auto rs = synth(rng, 3000, -9.0, 3.0, 4.0);
    std::uniform_real_distribution<double> flat(-50.0, 50.0);
    for (std::size_t i = 0; i < 9000; ++i)
    {
      ODIA::MassResidual d;
      d.mz = static_cast<float>(250.0 + (i % 1250));
      d.ppm = static_cast<float>(flat(rng));
      d.decoy = true;
      rs.push_back(d);
    }
    const double with_decoys = ODIA::MassCalibration::systematicResidualPpm(
      rs, [](double mz) { return -9.0 + 3.0 * std::log(mz / 700.0); });
    checkNear(with_decoys, 0.0, 0.50,
              "3x as many control cells do not move the statistic");
  }

  // --- failure is NaN, never zero -----------------------------------------
  {
    // Zero is the BEST possible score for this statistic, so returning it on
    // insufficient data reports failure as perfection -- and every
    // "smaller is better" comparison then silently prefers the failed model.
    // This is not hypothetical: three checks in this file passed vacuously
    // against an early version that returned 0.0 from a bailout path.
    const auto few = synth(rng, 20, -9.0, 3.0, 4.0);
    const double v = ODIA::MassCalibration::systematicResidualPpm(
      few, [](double) { return 0.0; });
    check(std::isnan(v), "below the binning floor it returns NaN, not 0");
  }

  // --- every candidate is scored on the SAME fragments ---------------------
  {
    // A model that shifts residuals must not be able to change WHICH residuals
    // it is judged on. The old version re-cut a 4-sigma core after applying
    // each correction, so a model could improve its score by pushing awkward
    // fragments out of its own core.
    //
    // Construct exactly that: a clean population plus a shoulder at +14 ppm.
    // A correction of -14 would move the shoulder onto zero and the bulk to
    // -14; if membership were recomputed per model it could keep the shoulder
    // and discard the bulk, and score well. With a fixed population it cannot.
    auto rs = synth(rng, 4000, 0.0, 0.0, 1.0);
    for (const auto& r : synth(rng, 800, 14.0, 0.0, 1.0)) { rs.push_back(r); }

    const double honest = ODIA::MassCalibration::systematicResidualPpm(
      rs, [](double) { return 0.0; });
    const double gaming = ODIA::MassCalibration::systematicResidualPpm(
      rs, [](double) { return 14.0; });
    check(std::isfinite(honest) && std::isfinite(gaming),
          "both candidates are scorable on the shouldered population");
    check(gaming > honest,
          "shifting onto the minority shoulder scores WORSE, not better");
  }

  // --- `im` survives a round trip through the struct -----------------------
  {
    ODIA::MassResidual r;
    check(std::isnan(r.im), "MassResidual::im defaults to NaN, not 0");
    r.im = 1.25f;
    check(r.im == 1.25f, "and carries a value when set");
  }

  // --- anchors survive the reordering finish() does --------------------
  {
    // finish() stable-sorts groups from extraction order into library order.
    // MassAnchor::group indexes that vector, so the sort must be applied
    // through a permutation and the anchors remapped. This models the
    // permutation directly: it is the invariant, independent of how the sort
    // is spelled. Two bugs of this exact family have already shipped here --
    // a staging buffer whose lifetime did not match its guard, and this.
    struct G { std::uint32_t precursor; float apex_rt; bool decoy; };
    std::vector<G> groups = {                       // extraction (RT) order
      {7, 10.0f, true }, {2, 20.0f, false}, {7, 30.0f, false}, {2, 40.0f, true }};
    std::vector<std::uint32_t> anchor_group = {0, 1, 2, 3};

    std::vector<std::uint32_t> order(groups.size());
    std::iota(order.begin(), order.end(), 0u);
    std::stable_sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t b) {
      if (groups[a].precursor != groups[b].precursor)
      { return groups[a].precursor < groups[b].precursor; }
      return groups[a].apex_rt < groups[b].apex_rt; });

    std::vector<std::uint32_t> moved_to(order.size());
    for (std::size_t n = 0; n < order.size(); ++n) { moved_to[order[n]] = std::uint32_t(n); }

    std::vector<G> before = groups;
    std::vector<G> reordered;
    for (const std::uint32_t o : order) { reordered.push_back(groups[o]); }
    for (std::uint32_t& g : anchor_group) { g = moved_to[g]; }

    check(reordered[0].precursor == 2 && reordered[3].precursor == 7,
          "the sort really does move groups (otherwise this proves nothing)");
    bool intact = true;
    for (std::size_t i = 0; i < anchor_group.size(); ++i)
    {
      const G& was = before[i];
      const G& now = reordered[anchor_group[i]];
      if (was.precursor != now.precursor || was.apex_rt != now.apex_rt ||
          was.decoy != now.decoy)
      { intact = false; }
    }
    check(intact, "every anchor still points at the group it was harvested from");

    // And the failure it is guarding against: without the remap, anchor 0 --
    // harvested from a decoy -- would read group 0 of the sorted vector, which
    // is a target.
    check(before[0].decoy && !reordered[0].decoy,
          "without the remap, a decoy's anchor would be read as a target");
  }

  std::printf("%s\n", failures == 0 ? "ALL PASSED" : "FAILURES");
  return failures == 0 ? 0 : 1;
}
