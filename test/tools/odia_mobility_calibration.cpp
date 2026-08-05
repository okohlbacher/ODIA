// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Does the per-run ion-mobility calibration recover a KNOWN 1/K0 error, does it
// refuse to invent one, and does it leave a run that has no mobility axis
// alone?
//
// With no arguments this runs on a SCRIPTED run -- every spectrum, peak,
// mobility and isolation window written out here -- so it is part of the ctest
// suite and needs no data file. It exercises collect() and fit() together
// rather than fit() alone, because the half that can be silently wrong is the
// half that decides what an observation of 1/K0 even is.
//
// Six properties, each a way this could be broken while still looking like it
// works:
//
//   1. A per-charge CONSTANT offset is recovered, separately for each charge,
//      and does not acquire an m/z shape it was not given.
//   2. An m/z-DEPENDENT offset comes back as a shape, right at both ends of the
//      range rather than only on average.
//   3. Fragments that agree with each other but sit at a RANDOM mobility are
//      refused. That is the flat case: the deltas still have a perfectly good
//      mode and width, they are just meaningless.
//   3b. A control that is as peaked as the data is refused too, which is the
//      other half of the gate and the half a target-only test leaves dead.
//   4. A run whose peaks carry NO ion mobility is a named no-op, not a failed
//      fit -- 12_80 is SCIEX SWATH and has no 1/K0 at all. Same for a library
//      that carries none.
//   5. The correction really is OUT OF SAMPLE. Two folds are planted with
//      DIFFERENT offsets, and each fold's anchors must come back corrected by
//      the other fold's number. An implementation that fitted on everything and
//      applied it to everything gives both of them the average, and fails.
//
// Given arguments it instead runs the real thing on a real run, which is how
// the numbers in the commit messages were measured. That path is deliberately
// NOT a ctest: it needs 14 GB of mzPeak and a minute of decode.
//
//   odia_mobility_calibration <library.tsv> <run.mzpeak> [precursors] [cycles]

#include <odia/DIANNLibraryFile.h>
#include <odia/LibraryGenerator.h>
#include <odia/MobilityCalibration.h>
#include <odia/SpectrumSource.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <limits>
#include <random>
#include <string>
#include <vector>

namespace
{
  const float NA = std::numeric_limits<float>::quiet_NaN();

  /// A run written out spectrum by spectrum, in acquisition order.
  class ScriptedRun : public ODIA::SpectrumSource
  {
  public:
    std::size_t addWindow(double mz_low, double mz_high, double im_low, double im_high)
    {
      ODIA::IsolationWindow w;
      w.mz_low = mz_low; w.mz_high = mz_high;
      w.im_low = im_low; w.im_high = im_high;
      windows_.push_back(w);
      return windows_.size() - 1;
    }

    std::size_t addSpectrum(std::size_t window, double rt)
    {
      ODIA::SpectrumInfo s;
      s.index = info_.size();
      s.retention_time = rt;
      s.window = windows_[window];
      info_.push_back(s);
      peaks_.emplace_back();
      return info_.size() - 1;
    }

    void addPeak(std::size_t spectrum, double mz, float intensity, float im)
    {
      auto& pk = peaks_[spectrum];
      pk.mz.push_back(mz);
      pk.intensity.push_back(intensity);
      // Either the run separates by mobility or it does not, so a spectrum
      // carries a mobility for every peak or for none.
      if (!std::isnan(im)) { pk.ion_mobility.push_back(im); }
    }

    /// Peaks are handed out in m/z order, which is what a reader guarantees and
    /// what the query index below assumes.
    void sortPeaks()
    {
      for (auto& pk : peaks_)
      {
        std::vector<std::size_t> order(pk.mz.size());
        for (std::size_t i = 0; i < order.size(); ++i) { order[i] = i; }
        std::sort(order.begin(), order.end(),
                  [&](std::size_t a, std::size_t b) { return pk.mz[a] < pk.mz[b]; });
        const auto permute = [&](auto& v) {
          if (v.empty()) { return; }
          std::decay_t<decltype(v)> tmp(v.size());
          for (std::size_t i = 0; i < order.size(); ++i) { tmp[i] = v[order[i]]; }
          v.swap(tmp);
        };
        permute(pk.mz); permute(pk.intensity); permute(pk.ion_mobility);
      }
    }

    const std::vector<ODIA::SpectrumInfo>& spectra() const override { return info_; }
    const std::vector<ODIA::IsolationWindow>& windows() const override { return windows_; }
    std::string describe() const override { return "scripted run"; }

    void peaks(std::size_t begin, std::size_t end,
               std::vector<ODIA::SpectrumPeaks>& out) override
    {
      out.assign(peaks_.begin() + begin, peaks_.begin() + end);
    }

  private:
    std::vector<ODIA::SpectrumInfo> info_;
    std::vector<ODIA::IsolationWindow> windows_;
    std::vector<ODIA::SpectrumPeaks> peaks_;
  };

  /// A library written out precursor by precursor.
  class ScriptedLibrary
  {
  public:
    std::size_t addPrecursor(double mz, std::uint8_t charge, float im)
    {
      auto& p = lib_.precursors();
      p.mz.push_back(ODIA::toFixed(mz));
      p.irt.push_back(0.0f);
      p.im.push_back(im);
      p.ccs.push_back(NA);
      p.charge.push_back(charge);
      p.decoy.push_back(0);
      p.modified_sequence.push_back(0);
      p.protein_group.push_back(0);
      p.transition_begin.push_back(
        static_cast<std::uint32_t>(lib_.transitions().product_mz.size()));
      p.transition_count.push_back(0);
      return p.mz.size() - 1;
    }

    void addTransition(double product_mz)
    {
      auto& t = lib_.transitions();
      t.product_mz.push_back(ODIA::toFixed(product_mz));
      t.library_intensity.push_back(1.0f);
      t.type.push_back(ODIA::FragmentType::Y);
      t.ordinal.push_back(4);
      t.charge.push_back(1);
      t.loss.push_back(ODIA::LossType::None);
      ++lib_.precursors().transition_count.back();
    }

    ODIA::Library& library() { return lib_; }

  private:
    ODIA::Library lib_;
  };

  int failures = 0;

  void check(bool ok, const std::string& what)
  {
    std::printf("  %s %s\n", ok ? "ok  " : "FAIL", what.c_str());
    if (!ok) { ++failures; }
  }

  // ------------------------------------------------------------------------
  // The synthetic run.
  //
  // Four isolation windows over 400-800 Th, each with a mobility band 0.60
  // wide -- roughly a diaPASEF frame, and comfortably wider than the +/-0.06
  // the calibration insists on having clear around a precursor's library 1/K0.
  // Each precursor gets eight fragments at m/z nobody else uses, emitted over
  // five mobility scans around its TRUE 1/K0 the way a TIMS ion actually
  // arrives, plus a flat background of unrelated peaks.
  // ------------------------------------------------------------------------
  struct Plan
  {
    std::size_t precursors = 1200;
    std::size_t cycles = 30;
    double present = 0.75;             ///< fraction actually in the run
    bool mobility = true;              ///< false = a run with no 1/K0 at all
    bool library_im = true;            ///< false = a library with no 1/K0
    bool random_mobility = false;      ///< a 1/K0 unrelated to the library's, but STABLE
    /// With `random_mobility`, the fraction that nevertheless sits where the
    /// library says. 0 is pure nonsense; a small value is a real peak too
    /// dilute to be evidence, which is what the peakedness gate is for.
    double true_fraction = 0.0;
    bool wandering = false;            ///< a 1/K0 that is different in every cycle
    bool peaked_control = false;       ///< also emit the m/z-shifted control, tighter
    std::size_t background = 2000;     ///< unrelated peaks per spectrum
    /// The truth: what to add to the library 1/K0 for this precursor.
    std::function<double(std::size_t /*index*/, std::uint8_t /*charge*/, double /*mz*/,
                         std::uint8_t /*charge2*/)> offset;
  };

  struct World
  {
    ScriptedRun run;
    ScriptedLibrary lib;
    std::vector<double> mz, im_library, planted;
    std::vector<std::uint8_t> charge;
    std::vector<char> present;
  };

  void build(World& w, const Plan& plan,
             const std::function<double(std::size_t, std::uint8_t, double)>& offset_of)
  {
    std::mt19937 rng(20260805u);
    const int NW = 4;
    for (int i = 0; i < NW; ++i)
    {
      w.run.addWindow(400.0 + 100.0 * i, 500.0 + 100.0 * i, 0.70, 1.30);
    }

    std::uniform_real_distribution<double> u01(0.0, 1.0);
    for (std::size_t i = 0; i < plan.precursors; ++i)
    {
      const double mz = 402.0 + 396.0 * double(i % 400) / 399.0;
      const std::uint8_t z = static_cast<std::uint8_t>(2 + (i % 3 == 2 ? 1 : 0));
      const double im = 0.80 + 0.40 * u01(rng);
      w.mz.push_back(mz);
      w.charge.push_back(z);
      w.im_library.push_back(im);
      w.planted.push_back(offset_of(i, z, mz));
      w.present.push_back(u01(rng) < plan.present ? 1 : 0);
      w.lib.addPrecursor(mz, z, plan.library_im ? static_cast<float>(im) : NA);
      for (int k = 0; k < 8; ++k)
      {
        // Distinct per (precursor, fragment) by construction, so nothing
        // matches by accident and every cluster found is the one planted.
        w.lib.addTransition(250.0 + 1.10 * double(i) + 137.0 * k);
      }
    }

    std::normal_distribution<double> jitter(0.0, 0.002);
    std::uniform_real_distribution<double> band(0.75, 1.25);
    // Drawn ONCE per precursor: a wrong 1/K0 that is nevertheless the same
    // every time the precursor is seen, which is the case the peakedness gate
    // has to catch. `wandering` is the other one, and the block rule catches it.
    std::vector<double> nonsense(plan.precursors);
    std::vector<char> real(plan.precursors, 0);
    for (std::size_t i = 0; i < plan.precursors; ++i)
    {
      nonsense[i] = band(rng);
      real[i] = u01(rng) < plan.true_fraction ? 1 : 0;
    }
    std::uniform_real_distribution<double> bg_mz(240.0, 1500.0);
    std::uniform_real_distribution<double> bright(2e3, 2e5);
    const double scan[5] = {-0.006, -0.003, 0.0, 0.003, 0.006};
    const float shape[5] = {0.3f, 0.7f, 1.0f, 0.7f, 0.3f};

    for (std::size_t c = 0; c < plan.cycles; ++c)
    {
      for (int win = 0; win < NW; ++win)
      {
        const std::size_t si = w.run.addSpectrum(win, 60.0 + 12.0 * double(c));
        const double lo = 400.0 + 100.0 * win, hi = 500.0 + 100.0 * win;
        for (std::size_t i = 0; i < plan.precursors; ++i)
        {
          if (!w.present[i]) { continue; }
          if (w.mz[i] < lo || w.mz[i] >= hi) { continue; }
          const double truth = plan.wandering ? band(rng)
                             : (plan.random_mobility && !real[i]) ? nonsense[i] + jitter(rng)
                             : w.im_library[i] + w.planted[i] + jitter(rng);
          const float I = static_cast<float>(bright(rng));
          for (int k = 0; k < 8; ++k)
          {
            const double fmz = 250.0 + 1.10 * double(i) + 137.0 * k;
            for (int s = 0; s < 5; ++s)
            {
              const double m = truth + scan[s];
              if (m < 0.70 || m >= 1.30) { continue; }
              w.run.addPeak(si, fmz, I * shape[s], plan.mobility ? float(m) : NA);
              if (plan.peaked_control)
              {
                // The same coherence at the control's m/z, and TIGHTER, so the
                // null is unambiguously at least as peaked as the data.
                w.run.addPeak(si, fmz + 7.33,
                              I * shape[s], plan.mobility ? float(truth + 0.25 * scan[s]) : NA);
              }
            }
          }
        }
        for (std::size_t b = 0; b < plan.background; ++b)
        {
          w.run.addPeak(si, bg_mz(rng), static_cast<float>(bright(rng)) * 0.2f,
                        plan.mobility ? float(band(rng)) : NA);
        }
      }
    }
    w.run.sortPeaks();
  }

  ODIA::MobilityCalibration::Options baseOptions()
  {
    ODIA::MobilityCalibration::Options opt;
    opt.fragment_ppm = 15.0;
    opt.min_per_bin = 30;
    opt.min_intensity_quantile = 0.0;   // the synthetic brightness carries no information
    return opt;
  }

  int selftest()
  {
    using MC = ODIA::MobilityCalibration;

    // ---- 1. a per-charge CONSTANT offset -----------------------------------
    {
      Plan plan;
      World w;
      build(w, plan, [](std::size_t, std::uint8_t z, double) {
        return z == 2 ? 0.020 : -0.015;
      });
      MC::Diagnostics d;
      const auto opt = baseOptions();
      const auto m = MC::calibrate(w.lib.library(), w.run, opt, &d);
      std::printf("1. a per-charge constant offset (+0.020 at 2+, -0.015 at 3+)\n%s",
                  MC::report(m, &d).c_str());
      check(m.fitted, "the gate passes");
      check(m.form == "constant", "and the model is CONSTANT, not a shape it was not given");
      double got2 = 0.0, got3 = 0.0;
      for (const auto& c : m.by_charge)
      {
        if (c.charge == 2) { got2 = c.constant; }
        if (c.charge == 3) { got3 = c.constant; }
      }
      check(std::abs(got2 - 0.020) < 0.002,
            "charge 2 recovered within 0.002 (got " + std::to_string(got2) + ")");
      check(std::abs(got3 - (-0.015)) < 0.002,
            "charge 3 recovered within 0.002 (got " + std::to_string(got3) + ")");
      check(got2 - got3 > 0.030,
            "the two charges are corrected SEPARATELY, not pooled");
      // A precursor that never appeared in the run was never an anchor, and is
      // the population this whole class exists for.
      std::size_t absent = 0;
      bool any = false;
      for (std::size_t i = 0; i < w.present.size(); ++i)
      {
        if (w.present[i] || w.charge[i] != 2) { continue; }
        absent = i;
        any = true;
        break;
      }
      check(any && std::abs(m.offsetFor(static_cast<std::uint32_t>(absent), w.mz[absent], 2)
                            - 0.020) < 0.002,
            "and a precursor the run never showed gets the same correction");
      check(m.squared_error_removed > 0.8,
            "most of the mean squared 1/K0 error is removed out of fold (" +
              std::to_string(m.squared_error_removed) + ")");
    }

    // ---- 2. an m/z-DEPENDENT offset ----------------------------------------
    {
      Plan plan;
      World w;
      const auto truth = [](double mz) { return -0.012 + 0.030 * (mz - 400.0) / 400.0; };
      build(w, plan, [&](std::size_t, std::uint8_t, double mz) { return truth(mz); });
      MC::Diagnostics d;
      auto opt = baseOptions();
      const auto m = MC::calibrate(w.lib.library(), w.run, opt, &d);
      std::printf("2. an offset that is a function of m/z (-0.012 at 400 Th, +0.018 at 800)\n%s",
                  MC::report(m, &d).c_str());
      check(m.fitted, "the gate passes");
      check(m.form == "mz_shaped", "the model is SHAPED, which a constant-only fit would miss");
      // Read the applied correction at both ends, through a precursor that was
      // never an anchor so the answer comes from the all-anchor model.
      double err_low = 0.0, err_high = 0.0;
      std::size_t n_low = 0, n_high = 0;
      for (std::size_t i = 0; i < w.present.size(); ++i)
      {
        if (w.present[i]) { continue; }
        const double got = m.offsetFor(static_cast<std::uint32_t>(i), w.mz[i], w.charge[i]);
        if (w.mz[i] < 460.0) { err_low += std::abs(got - truth(w.mz[i])); ++n_low; }
        if (w.mz[i] > 740.0) { err_high += std::abs(got - truth(w.mz[i])); ++n_high; }
      }
      err_low = n_low ? err_low / double(n_low) : 1.0;
      err_high = n_high ? err_high / double(n_high) : 1.0;
      std::printf("     mean |error| of the applied correction: %.5f below 460 Th (n=%zu), "
                  "%.5f above 740 Th (n=%zu)\n", err_low, n_low, err_high, n_high);
      check(err_low < 0.004 && err_high < 0.004,
            "the correction is right at BOTH ends of the range, not just on average");
    }

    // ---- 3. coherent, reproducible fragments at a MEANINGLESS 1/K0 --------
    {
      Plan plan;
      plan.random_mobility = true;
      World w;
      build(w, plan, [](std::size_t, std::uint8_t, double) { return 0.0; });
      MC::Diagnostics d;
      const auto m = MC::calibrate(w.lib.library(), w.run, baseOptions(), &d);
      std::printf("3. fragments that agree with each other at a RANDOM mobility\n%s",
                  MC::report(m, &d).c_str());
      check(!m.fitted, "the gate REFUSES a flat 1/K0 residual");
      check(m.run_has_mobility, "and says so as a failed fit, not as a missing axis");
      check(m.offsetFor(0, 600.0, 2) == 0.0 && m.offsetFor(7, 450.0, 3) == 0.0,
            "nothing is applied at any m/z or charge");
    }

    // ---- 3b. a control as peaked as the data must be refused ---------------
    {
      Plan plan;
      plan.peaked_control = true;
      World w;
      build(w, plan, [](std::size_t, std::uint8_t, double) { return 0.020; });
      MC::Diagnostics d;
      const auto m = MC::calibrate(w.lib.library(), w.run, baseOptions(), &d);
      std::printf("3b. the m/z-shifted control given the same structure, only tighter\n%s",
                  MC::report(m, &d).c_str());
      check(!m.fitted, "the gate REFUSES when the null is as peaked as the data");
      check(m.decoy_peakedness >= m.peakedness,
            "and that is why (" + std::to_string(m.decoy_peakedness) + " vs " +
              std::to_string(m.peakedness) + ")");
    }

    // ---- 3d. a real peak, too dilute to be evidence ------------------------
    // One precursor in twelve is where the library says; the rest are somewhere
    // stable and unrelated. There IS a peak, the mode would find it, and a
    // scale can be computed -- so nothing degenerate refuses this. Only the
    // peakedness gate does, which is the point: it refuses on the STRENGTH of
    // the evidence, not on whether an answer can be computed.
    {
      Plan plan;
      plan.random_mobility = true;
      plan.true_fraction = 0.04;
      plan.precursors = 2000;
      plan.present = 1.0;
      World w;
      build(w, plan, [](std::size_t, std::uint8_t, double) { return 0.015; });
      MC::Diagnostics d;
      const auto m = MC::calibrate(w.lib.library(), w.run, baseOptions(), &d);
      std::printf("3d. a real +0.015 offset on one precursor in twenty-five\n%s",
                  MC::report(m, &d).c_str());
      check(!m.fitted, "the gate refuses a peak too dilute to be evidence");
      check(m.peakedness > 1.2 && m.peakedness < baseOptions().min_peakedness,
            "and it is the PEAKEDNESS that refuses it (" + std::to_string(m.peakedness) +
              "), not a degenerate scale");
      check(m.offsetFor(3, 500.0, 2) == 0.0, "so nothing is applied");
    }

    // ---- 3c. a 1/K0 that does not survive the next cycle -------------------
    // The block rule's own test. These fragments agree with each other inside
    // every spectrum -- they would have passed the co-occurrence filter on its
    // own -- but the agreement is somewhere else in the next cycle, which is
    // what a coincidence looks like and what an eluting precursor does not.
    {
      Plan plan;
      plan.wandering = true;
      World w;
      build(w, plan, [](std::size_t, std::uint8_t, double) { return 0.0; });
      MC::Diagnostics d;
      const auto m = MC::calibrate(w.lib.library(), w.run, baseOptions(), &d);
      std::printf("3c. clusters that move to a different 1/K0 every cycle\n%s",
                  MC::report(m, &d).c_str());
      check(!m.fitted, "refused");
      check(m.residuals <= 5,
            "and refused because almost nothing was REPRODUCED across a block -- 913 cells "
            "when the same fragments hold still, " + std::to_string(m.residuals) + " when "
            "they do not");
    }

    // ---- 4. no mobility axis at all ---------------------------------------
    {
      Plan plan;
      plan.mobility = false;
      plan.cycles = 30;
      World w;
      build(w, plan, [](std::size_t, std::uint8_t, double) { return 0.020; });
      MC::Diagnostics d;
      const auto m = MC::calibrate(w.lib.library(), w.run, baseOptions(), &d);
      std::printf("4. a run whose peaks carry no 1/K0 (the SCIEX SWATH case)\n%s",
                  MC::report(m, &d).c_str());
      check(!m.run_has_mobility, "the run is reported as having NO mobility axis");
      check(!m.fitted && m.form == "none", "nothing is fitted");
      check(m.offsetFor(0, 600.0, 2) == 0.0, "and nothing is applied");
      check(m.reason.find("no ion mobility") != std::string::npos,
            "and it SAYS so, rather than reporting a failed fit");
      check(d.spectra_decoded <= 24 && d.spectra_decoded < 120,
            "it also stops after the first block rather than probing the whole run (" +
              std::to_string(d.spectra_decoded) + " of 120 spectra)");
    }

    // ---- 4b. a library with no 1/K0 ---------------------------------------
    {
      Plan plan;
      plan.library_im = false;
      plan.cycles = 10;
      World w;
      build(w, plan, [](std::size_t, std::uint8_t, double) { return 0.0; });
      MC::Diagnostics d;
      const auto m = MC::calibrate(w.lib.library(), w.run, baseOptions(), &d);
      std::printf("4b. a library that carries no 1/K0\n%s", MC::report(m, &d).c_str());
      check(!m.library_has_mobility && !m.fitted, "reported as nothing to correct");
      check(m.offsetFor(0, 600.0, 2) == 0.0, "and nothing is applied");
    }

    // ---- 5. the correction is OUT OF SAMPLE -------------------------------
    // Two folds, planted with DIFFERENT offsets. Each fold's anchors must come
    // back corrected by the OTHER fold's number, because the model that touches
    // them is the one fitted without them. Fit-on-everything gives both the
    // average and fails both checks; fit-per-fold-but-apply-your-own gives each
    // its own number and fails them the other way.
    {
      Plan plan;
      plan.present = 1.0;
      World w;
      auto opt = baseOptions();
      opt.folds = 2;
      const auto planted = [&](std::size_t i, std::uint8_t, double) {
        return ODIA::MobilityCalibration::foldIndex(static_cast<std::uint32_t>(i), 2) == 0
                 ? 0.005 : 0.025;
      };
      build(w, plan, planted);
      MC::Diagnostics d;
      const auto m = MC::calibrate(w.lib.library(), w.run, opt, &d);
      std::printf("5. two folds planted at +0.005 and +0.025\n%s", MC::report(m, &d).c_str());
      check(m.fitted, "the gate passes");
      double e0 = 0.0, e1 = 0.0;
      std::size_t n0 = 0, n1 = 0;
      for (std::size_t i = 0; i < w.mz.size(); ++i)
      {
        const auto id = static_cast<std::uint32_t>(i);
        if (m.foldOf(id) >= m.folds) { continue; }         // never an anchor
        const double got = m.offsetFor(id, w.mz[i], w.charge[i]);
        if (MC::foldIndex(id, 2) == 0) { e0 += got; ++n0; } else { e1 += got; ++n1; }
      }
      e0 = n0 ? e0 / double(n0) : 0.0;
      e1 = n1 ? e1 / double(n1) : 0.0;
      std::printf("     fold 0 anchors (planted +0.005) are corrected by %+.5f\n"
                  "     fold 1 anchors (planted +0.025) are corrected by %+.5f\n", e0, e1);
      check(n0 > 100 && n1 > 100, "both folds have anchors");
      check(std::abs(e0 - 0.025) < 0.003,
            "fold 0 is corrected by fold 1's number -- the fit that saw it did not");
      check(std::abs(e1 - 0.005) < 0.003,
            "fold 1 is corrected by fold 0's number");
    }

    // ---- 6. nothing at all -------------------------------------------------
    {
      const auto m = MC::fit({}, baseOptions());
      check(!m.fitted && m.offsetFor(0, 700.0, 2) == 0.0,
            "an empty sample is a no-op, not a crash");
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
    // The null is the library's own decoys, so this path has to build them the
    // way the tool does before extracting -- otherwise the measurement is made
    // against a different control than the one that runs.
    if (library.decoyCount() == 0)
    {
      ODIA::LibraryGenerator::appendDecoys(library, ODIA::DecoyMethod::Mutate, nullptr);
    }
    std::printf("library: %zu precursors (%zu decoys), %zu transitions\n",
                library.precursorCount(), library.decoyCount(), library.transitionCount());

    auto run = ODIA::openRun(argv[2]);
    std::printf("%s\n", run->describe().c_str());

    ODIA::MobilityCalibration::Options opt;
    if (argc > 3) { opt.max_precursors = std::strtoul(argv[3], nullptr, 10); }
    if (argc > 4) { opt.cycles = std::strtoul(argv[4], nullptr, 10); }
    if (const char* s = std::getenv("ODIA_IMCAL_SEARCH")) { opt.search_im = std::atof(s); }
    if (const char* s = std::getenv("ODIA_IMCAL_GATE")) { opt.gate_im = std::atof(s); }
    if (const char* s = std::getenv("ODIA_IMCAL_CLUSTER")) { opt.cluster_im = std::atof(s); }
    if (const char* s = std::getenv("ODIA_IMCAL_FRAGMENTS"))
    { opt.min_fragments_matched = std::strtoul(s, nullptr, 10); }
    if (const char* s = std::getenv("ODIA_IMCAL_BLOCK"))
    { opt.cycle_block = std::strtoul(s, nullptr, 10); }
    if (const char* s = std::getenv("ODIA_IMCAL_MIN_CYCLES"))
    { opt.min_cycles_matched = std::strtoul(s, nullptr, 10); }
    if (const char* s = std::getenv("ODIA_IMCAL_PPM")) { opt.fragment_ppm = std::atof(s); }
    if (const char* s = std::getenv("ODIA_IMCAL_PPM_OFFSET"))
    { opt.fragment_ppm_offset = std::atof(s); }
    if (const char* s = std::getenv("ODIA_IMCAL_PPM_LOG_SLOPE"))
    { opt.fragment_ppm_log_slope = std::atof(s); }
    if (const char* s = std::getenv("ODIA_IMCAL_PPM_REF_MZ"))
    { opt.fragment_ppm_ref_mz = std::atof(s); }
    if (const char* s = std::getenv("ODIA_IMCAL_QUANTILE"))
    { opt.min_intensity_quantile = std::atof(s); }

    ODIA::MobilityCalibration::Diagnostics d;
    const auto m = ODIA::MobilityCalibration::calibrate(library, *run, opt, &d);
    std::printf("%s", ODIA::MobilityCalibration::report(m, &d).c_str());

    if (const char* path = std::getenv("ODIA_IMCAL_DUMP"))
    {
      const auto r = ODIA::MobilityCalibration::collect(library, *run, opt, nullptr);
      std::FILE* f = std::fopen(path, "wb");
      if (f)
      {
        std::fprintf(f, "precursor\tmz\tcharge\tim_library\tim_observed\tdelta\tintensity\t"
                        "rt\tfragments\tcycles\tdecoy\n");
        for (const auto& x : r)
        {
          std::fprintf(f, "%u\t%.5f\t%u\t%.6f\t%.6f\t%.6f\t%.6g\t%.2f\t%u\t%u\t%d\n",
                       x.precursor, x.mz, unsigned(x.charge), x.im_library, x.im_observed,
                       x.delta, x.intensity, x.rt, unsigned(x.fragments), unsigned(x.cycles),
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
