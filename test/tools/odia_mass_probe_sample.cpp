// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// -mass_probe_source scored: the fragment mass probe on a sample drawn from
// SCORED peak groups instead of from the library stride.
//
// Runs on a scripted run and library written out inside the test, so it needs
// no data. Four properties:
//
//   1. MassProbeSample::fromGroups picks the right rows: the best row per
//      target precursor by q (then d-score), only at q <= the bar; the same
//      NUMBER of best-scoring decoys; non-finite scores never; and a cap that
//      strides over the precursor-sorted list rather than keeping the top.
//   2. Identity: handing collect() the library stride's own list as an
//      explicit sample, with no RT and no 1/K0, returns the SAME residuals,
//      bit for bit. The explicit path is a different source, not a different
//      probe.
//   3. The defect and the fix, on a run where 2% of a library is present
//      and half of the absent precursors' cells fill with stable
//      interference at random per-fragment mass errors (what an absent precursor looks like
//      in a dense diaPASEF frame): the library stride does NOT recover the
//      planted -4 ppm, and the scored sample -- the present precursors at their
//      apex RT and observed 1/K0 -- does, through the unchanged gate.
//   4. The control: the same explicit path on ABSENT precursors (interference
//      only), at arbitrary "apexes", must not pass. The explicit path cannot
//      manufacture a mode.
//   5. A malformed explicit sample throws.
//   6. `use_sample` with an EMPTY sample yields no residuals and an unfitted
//      model -- never the library stride it used to fall through to -- and is
//      a no-op on a non-empty sample.
//   7. The cycle cap (`sample_max_cycles`): at or above the reachable cycles
//      it is bit-identical to no cap; binding, it visits at most that many,
//      drops the entries no drawn cycle reaches, still recovers the planted
//      offset, and the precursor cap strides what it kept.

#include <odia/DIANNLibraryFile.h>
#include <odia/LibraryGenerator.h>
#include <odia/MassCalibration.h>
#include <odia/MassProbeSample.h>
#include <odia/PeakGroupScorer.h>
#include <odia/SpectrumSource.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
  const float NA = std::numeric_limits<float>::quiet_NaN();

  int failures = 0;

  void check(bool ok, const std::string& what)
  {
    std::printf("  %s %s\n", ok ? "ok  " : "FAIL", what.c_str());
    if (!ok) { ++failures; }
  }

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
      pk.ion_mobility.push_back(im);
    }

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

  void addPrecursor(ODIA::Library& lib, double mz, float im, bool decoy)
  {
    auto& p = lib.precursors();
    p.mz.push_back(ODIA::toFixed(mz));
    p.irt.push_back(0.0f);
    p.im.push_back(im);
    p.ccs.push_back(NA);
    p.charge.push_back(2);
    p.decoy.push_back(decoy ? 1 : 0);
    p.modified_sequence.push_back(0);
    p.protein_group.push_back(0);
    p.transition_begin.push_back(static_cast<std::uint32_t>(lib.transitions().product_mz.size()));
    p.transition_count.push_back(0);
  }

  void addTransition(ODIA::Library& lib, double product_mz)
  {
    auto& t = lib.transitions();
    t.product_mz.push_back(ODIA::toFixed(product_mz));
    t.library_intensity.push_back(1.0f);
    t.type.push_back(ODIA::FragmentType::Y);
    t.ordinal.push_back(4);
    t.charge.push_back(1);
    t.loss.push_back(ODIA::LossType::None);
    ++lib.precursors().transition_count.back();
  }

  ODIA::PeakGroupScorer::PeakGroup group(std::uint32_t precursor, bool decoy, double q,
                                         double dscore, float rt, float im)
  {
    ODIA::PeakGroupScorer::PeakGroup g;
    g.precursor = precursor;
    g.decoy = decoy;
    g.qvalue = q;
    g.dscore = dscore;
    g.apex_rt = rt;
    g.observed_im = im;
    return g;
  }

  // ---- 1. the selection ----------------------------------------------------
  void selection()
  {
    std::printf("1. MassProbeSample::fromGroups\n");
    std::vector<ODIA::PeakGroupScorer::PeakGroup> g;
    g.push_back(group(5, false, 0.005, 2.0, 100.0f, 0.90f));   // 0
    g.push_back(group(5, false, 0.002, 1.0, 200.0f, 0.91f));   // 1: best of 5 (lower q)
    g.push_back(group(7, false, 0.020, 9.0, 300.0f, 0.92f));   // 2: above the bar
    g.push_back(group(3, false, 0.010, 1.5, 400.0f, 0.93f));   // 3: at the bar -> in
    g.push_back(group(3, false, 0.010, 1.7, 450.0f, 0.94f));   // 4: same q, higher d -> best
    g.push_back(group(100, true, 1.0, 3.0, 500.0f, 1.00f));    // 5
    g.push_back(group(101, true, 1.0, 5.0, 600.0f, 1.01f));    // 6
    g.push_back(group(101, true, 1.0, 4.0, 650.0f, 1.02f));    // 7: worse row of 101
    g.push_back(group(102, true, 1.0, 1.0, 700.0f, 1.03f));    // 8: third decoy, cut by matching
    g.push_back(group(103, true, 1.0, std::numeric_limits<double>::quiet_NaN(), 800.0f, NA));
    g.push_back(group(9, false, 0.001, std::numeric_limits<double>::quiet_NaN(), 900.0f, NA));

    const auto s = ODIA::MassProbeSample::fromGroups(g, 0.01, 0);
    check(s.eligible_targets == 2, "two targets pass q <= 0.01 (got " +
                                     std::to_string(s.eligible_targets) + ")");
    check(s.targets.precursor == std::vector<std::uint32_t>({3, 5}),
          "targets are precursors 3 and 5, sorted by library index");
    check(s.targets.size() == 2 && s.targets.rt[0] == 450.0f && s.targets.rt[1] == 200.0f,
          "each target carries its BEST row's apex (q first, then d-score)");
    check(s.targets.size() == 2 && s.targets.im[0] == 0.94f && s.targets.im[1] == 0.91f,
          "... and that row's observed 1/K0");
    check(s.available_decoys == 3, "three decoy precursors have a finite-scored row (got " +
                                     std::to_string(s.available_decoys) + ")");
    check(s.decoys.precursor == std::vector<std::uint32_t>({100, 101}),
          "the control is the 2 best-scoring decoys (as many as eligible targets)");
    check(s.decoys.size() == 2 && s.decoys.rt[1] == 600.0f,
          "a decoy carries its best-scoring row");

    // The cap strides over the precursor-sorted list.
    std::vector<ODIA::PeakGroupScorer::PeakGroup> h;
    for (std::uint32_t i = 0; i < 10; ++i)
    {
      h.push_back(group(100 - 3 * i, false, 0.001, double(i), float(i), 1.0f));
    }
    const auto c = ODIA::MassProbeSample::fromGroups(h, 0.01, 4);
    // sorted precursors: 73,76,...,100; stride k*10/4 -> 0,2,5,7 -> 73,79,88,94
    check(c.eligible_targets == 10 && c.targets.precursor ==
                                        std::vector<std::uint32_t>({73, 79, 88, 94}),
          "a cap keeps the library stride over precursor order, not the top scores");
    check(c.decoys.size() == 0, "no decoy rows, no control");
  }

  // ---- the world for 2-4 ------------------------------------------------------
  struct World
  {
    ScriptedRun run;
    ODIA::Library lib;
    std::vector<double> mz, im;
    std::vector<char> present, interferes;
    std::vector<double> apex_rt;
  };

  constexpr double PLANTED_PPM = -4.0;

  void build(World& w)
  {
    std::mt19937 rng(20261003u);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    const int NW = 4;
    const std::size_t N = 8000, CYCLES = 300;
    const double CYCLE_S = 2.0;
    for (int i = 0; i < NW; ++i) { w.run.addWindow(400.0 + 100.0 * i, 500.0 + 100.0 * i, 0.6, 1.4); }

    for (std::size_t i = 0; i < N; ++i)
    {
      const double mz = 401.0 + 398.0 * u01(rng);
      const double im = 0.75 + 0.5 * u01(rng);
      w.mz.push_back(mz);
      w.im.push_back(im);
      w.present.push_back(i % 50 == 7 ? 1 : 0);    // 160 of 8,000 = 2%
      w.interferes.push_back(!w.present.back() && u01(rng) < 0.5 ? 1 : 0);
      w.apex_rt.push_back(CYCLE_S * (10.0 + double(CYCLES - 20) * u01(rng)));
      addPrecursor(w.lib, mz, static_cast<float>(im), false);
      for (int k = 0; k < 8; ++k) { addTransition(w.lib, 250.0 + 0.2913 * double(i) + 141.37 * k); }
    }

    std::normal_distribution<double> jitter(0.0, 1.5);   // ppm, per peak
    std::uniform_real_distribution<double> wrong(-40.0, 40.0);
    std::uniform_real_distribution<double> bright(2e3, 2e5);
    std::uniform_real_distribution<double> bg_mz(240.0, 1500.0);
    std::uniform_real_distribution<double> bg_im(0.6, 1.4);
    // One stable WRONG mass error per (interfering absent precursor, fragment):
    // reproducible in every cycle, at the precursor's 1/K0, and not the
    // precursor -- each fragment's slot holds some unrelated ion, so the errors
    // are independent across fragments.
    //
    // NOT one error per precursor shared by all its fragments. That variant was
    // tried first and the UNCHANGED gate passed 160 such precursors at +28 ppm
    // (peakedness 3.35, control 2.00): 8 identical residuals per cell are 160
    // point masses, not 1,280 independent draws, and the gate counts them as
    // the latter. A property of the gate's evidence counting, not of where the
    // sample comes from; recorded here so the model choice is not mistaken
    // for a free parameter.
    std::vector<double> wrong_ppm(N * 8);
    for (auto& x : wrong_ppm) { x = wrong(rng); }

    const auto& tr = w.lib.transitions();
    const auto& p = w.lib.precursors();
    for (std::size_t c = 0; c < CYCLES; ++c)
    {
      for (int win = 0; win < NW; ++win)
      {
        const double rt = CYCLE_S * double(c) + 0.4 * win;
        const std::size_t si = w.run.addSpectrum(win, rt);
        const double lo = 400.0 + 100.0 * win, hi = 500.0 + 100.0 * win;
        for (std::size_t i = 0; i < N; ++i)
        {
          if (w.mz[i] < lo || w.mz[i] >= hi) { continue; }
          double scale = 0.0, ppm = 0.0;
          if (w.present[i])
          {
            const double d = (rt - w.apex_rt[i]) / 4.0;
            if (std::abs(d) > 3.0) { continue; }
            scale = std::exp(-0.5 * d * d);
            ppm = PLANTED_PPM;
          }
          else if (w.interferes[i])
          {
            scale = 1.0;
            ppm = 0.0;   // per fragment, below
          }
          else { continue; }
          for (std::size_t k = 0; k < p.transition_count[i]; ++k)
          {
            const double f = ODIA::fromFixed(tr.product_mz[p.transition_begin[i] + k]);
            const double e = w.present[i] ? ppm + jitter(rng) : wrong_ppm[i * 8 + k];
            w.run.addPeak(si, f * (1.0 + e * 1e-6), static_cast<float>(bright(rng) * scale),
                          static_cast<float>(w.im[i]));
          }
        }
        for (int b = 0; b < 300; ++b)
        {
          w.run.addPeak(si, bg_mz(rng), static_cast<float>(100.0 + 900.0 * u01(rng)),
                        static_cast<float>(bg_im(rng)));
        }
      }
    }
    w.run.sortPeaks();
  }

  ODIA::MassCalibration::Options baseOptions()
  {
    ODIA::MassCalibration::Options o;
    o.threads = 1;
    return o;
  }

  // ---- 2. identity -------------------------------------------------------------
  void identity(World& w)
  {
    std::printf("2. an explicit sample equal to the stride is the stride, bit for bit\n");
    auto o = baseOptions();
    const auto a = ODIA::MassCalibration::collect(w.lib, w.run, o, nullptr);

    // RunProbe::samplePrecursors' stride, rebuilt from its definition: targets,
    // valid m/z, >= min_fragments_matched transitions, eligible[k * n / want].
    std::vector<std::uint32_t> eligible;
    const auto& p = w.lib.precursors();
    for (std::size_t i = 0; i < w.lib.precursorCount(); ++i)
    {
      if (p.decoy[i]) { continue; }
      if (p.transition_count[i] < o.min_fragments_matched) { continue; }
      eligible.push_back(static_cast<std::uint32_t>(i));
    }
    const std::size_t want = std::min(o.max_precursors, eligible.size());
    for (std::size_t k = 0; k < want; ++k) { o.sample.push_back(eligible[k * eligible.size() / want]); }
    const auto b = ODIA::MassCalibration::collect(w.lib, w.run, o, nullptr);

    bool same = a.size() == b.size();
    for (std::size_t i = 0; same && i < a.size(); ++i)
    {
      same = std::memcmp(&a[i].mz, &b[i].mz, sizeof a[i].mz) == 0 &&
             std::memcmp(&a[i].ppm, &b[i].ppm, sizeof a[i].ppm) == 0 &&
             std::memcmp(&a[i].rt, &b[i].rt, sizeof a[i].rt) == 0 &&
             std::memcmp(&a[i].intensity, &b[i].intensity, sizeof a[i].intensity) == 0 &&
             a[i].decoy == b[i].decoy;
    }
    check(!a.empty() && same, "identical residual vectors (" + std::to_string(a.size()) +
                                " vs " + std::to_string(b.size()) + ")");
  }

  // ---- 3 and 4. the defect, the fix, the control -------------------------------
  void probe(World& w)
  {
    std::printf("3. library stride vs scored sample, 2%% present\n");
    {
      auto o = baseOptions();
      ODIA::MassCalibration::Diagnostics d;
      const auto m = ODIA::MassCalibration::calibrate(w.lib, w.run, o, &d);
      std::printf("   stride:\n%s", ODIA::MassCalibration::report(m, &d).c_str());
      check(!(m.fitted && std::abs(m.ppmAt(700.0) - PLANTED_PPM) < 1.0),
            "the library stride does NOT recover the planted offset");
    }
    {
      auto o = baseOptions();
      for (std::size_t i = 0; i < w.mz.size(); ++i)
      {
        if (!w.present[i]) { continue; }
        o.sample.push_back(static_cast<std::uint32_t>(i));
        // An apex as pass 1 would report it: off by up to a cycle.
        o.sample_rt.push_back(static_cast<float>(w.apex_rt[i] + (i % 3 == 0 ? 1.7 : -0.9)));
        o.sample_im.push_back(static_cast<float>(w.im[i] + 0.002));
      }
      o.sample_rt_window = 5.0;
      ODIA::MassCalibration::Diagnostics d;
      const auto m = ODIA::MassCalibration::calibrate(w.lib, w.run, o, &d);
      std::printf("   scored:\n%s", ODIA::MassCalibration::report(m, &d).c_str());
      check(d.sample_entries == 160, "all 160 sampled entries are eligible");
      check(d.sample_cycles > 0 && d.sample_cycles < 300,
            "only cycles near an apex are visited (" + std::to_string(d.sample_cycles) + ")");
      check(m.fitted, "the scored sample PASSES the unchanged gate");
      check(std::abs(m.ppmAt(700.0) - PLANTED_PPM) < 0.7,
            "and recovers the planted " + std::to_string(PLANTED_PPM) + " ppm (got " +
              std::to_string(m.ppmAt(700.0)) + ")");
      // Not "more peaked than its shifted control": the control here is a
      // handful of residuals, and peakednessRatio returns its 1e9 sentinel for
      // an empty edge band -- the comparison would read a sentinel.
      check(m.residuals >= 200 && m.peakedness >= 3.0,
            "from >= 200 residuals at peakedness >= 3 (got " + std::to_string(m.residuals) +
              ", " + std::to_string(m.peakedness) + ")");
    }
    std::printf("4. the explicit path on absent precursors must not pass\n");
    {
      auto o = baseOptions();
      std::size_t n = 0;
      for (std::size_t i = 0; i < w.mz.size() && n < 160; ++i)
      {
        if (w.present[i] || !w.interferes[i]) { continue; }
        o.sample.push_back(static_cast<std::uint32_t>(i));
        o.sample_rt.push_back(static_cast<float>(w.apex_rt[i]));
        o.sample_im.push_back(NA);   // falls back to the library column
        ++n;
      }
      o.sample_rt_window = 5.0;
      ODIA::MassCalibration::Diagnostics d;
      const auto m = ODIA::MassCalibration::calibrate(w.lib, w.run, o, &d);
      std::printf("   absent:\n%s", ODIA::MassCalibration::report(m, &d).c_str());
      check(!m.fitted, "the gate FAILS on stable interference at random per-fragment mass errors");
    }
    std::printf("5. a malformed explicit sample is refused, not misread\n");
    {
      auto o = baseOptions();
      o.sample = {1, 2, 3};
      o.sample_rt = {1.0f, 2.0f};
      bool threw = false;
      try { (void)ODIA::MassCalibration::collect(w.lib, w.run, o, nullptr); }
      catch (const std::invalid_argument&) { threw = true; }
      check(threw, "sample_rt not parallel to sample throws");
    }
  }

  // ---- 6. an empty explicit sample is nothing to probe -------------------------
  //
  // The defect: with the explicit path keyed on `!sample.empty()`, a scored
  // sample that came out EMPTY (pass 1 scored no target at the q bar) took the
  // library stride -- the noise probe -- and was logged as measured at pass 1's
  // groups. `use_sample` makes the empty sample mean what it says.
  void emptySample(World& w)
  {
    std::printf("6. use_sample with an EMPTY sample yields nothing, not the stride\n");
    auto o = baseOptions();
    const auto stride = ODIA::MassCalibration::collect(w.lib, w.run, o, nullptr);
    check(!stride.empty(), "control: the same options without use_sample take the stride (" +
                             std::to_string(stride.size()) + " residuals)");

    o.use_sample = true;   // sample, sample_rt, sample_im all left empty
    ODIA::MassCalibration::Diagnostics d;
    const auto r = ODIA::MassCalibration::collect(w.lib, w.run, o, &d);
    check(r.empty(), "0 residuals (got " + std::to_string(r.size()) + ")");
    check(d.spectra_decoded == 0 && d.sample_entries == 0,
          "no spectrum decoded and no entry sampled");
    const auto m = ODIA::MassCalibration::calibrate(w.lib, w.run, o, nullptr);
    check(!m.fitted && m.residuals == 0 && m.form == "none",
          "an unfitted model from 0 residuals (" + m.reason + ")");

    // ...and use_sample on a NON-empty sample changes nothing.
    auto a = baseOptions();
    for (std::uint32_t i = 7; i < 8000; i += 50) { a.sample.push_back(i); }
    auto b = a;
    b.use_sample = true;
    const auto ra = ODIA::MassCalibration::collect(w.lib, w.run, a, nullptr);
    const auto rb = ODIA::MassCalibration::collect(w.lib, w.run, b, nullptr);
    bool same = !ra.empty() && ra.size() == rb.size();
    for (std::size_t i = 0; same && i < ra.size(); ++i)
    {
      same = std::memcmp(&ra[i].ppm, &rb[i].ppm, sizeof ra[i].ppm) == 0 &&
             std::memcmp(&ra[i].mz, &rb[i].mz, sizeof ra[i].mz) == 0;
    }
    check(same, "use_sample on a non-empty sample is bit-identical to leaving it unset");
  }

  // ---- 7. the cycle cap ---------------------------------------------------------
  //
  // Uncapped, the probe visits every cycle near any apex, serially; on a long
  // run that is nearly the whole run. The cap draws the library probe's
  // stratified sample from the reachable cycles and keeps only the entries a
  // drawn cycle reaches; the precursor cap is applied after it.
  void cycleCap(World& w)
  {
    std::printf("7. -mass_probe_max_cycles: the cap bounds the visited cycles\n");
    const auto scored = [&]() {
      auto o = baseOptions();
      for (std::size_t i = 0; i < w.mz.size(); ++i)
      {
        if (!w.present[i]) { continue; }
        o.sample.push_back(static_cast<std::uint32_t>(i));
        o.sample_rt.push_back(static_cast<float>(w.apex_rt[i] + (i % 3 == 0 ? 1.7 : -0.9)));
        o.sample_im.push_back(static_cast<float>(w.im[i] + 0.002));
      }
      o.use_sample = true;
      o.sample_rt_window = 5.0;
      return o;
    };
    const auto same = [](const std::vector<ODIA::MassResidual>& a,
                         const std::vector<ODIA::MassResidual>& b) {
      if (a.empty() || a.size() != b.size()) { return false; }
      for (std::size_t i = 0; i < a.size(); ++i)
      {
        if (std::memcmp(&a[i].ppm, &b[i].ppm, sizeof a[i].ppm) != 0 ||
            std::memcmp(&a[i].rt, &b[i].rt, sizeof a[i].rt) != 0) { return false; }
      }
      return true;
    };

    auto o0 = scored();
    ODIA::MassCalibration::Diagnostics d0;
    const auto r0 = ODIA::MassCalibration::collect(w.lib, w.run, o0, &d0);
    check(d0.sample_cycles == d0.sample_cycles_available && d0.sample_cycles > 100,
          "uncapped: every reachable cycle is visited (" + std::to_string(d0.sample_cycles) +
            " of " + std::to_string(d0.sample_cycles_available) + ")");
    check(d0.sample_entries_in_reach == 160 && d0.sample_entries_probed == 160,
          "uncapped: all 160 entries probed");

    // A cap at or above the reachable count is no cap: bit-identical.
    auto o1 = scored();
    o1.sample_max_cycles = d0.sample_cycles_available;
    const auto r1 = ODIA::MassCalibration::collect(w.lib, w.run, o1, nullptr);
    check(same(r0, r1), "a cap >= the reachable cycles changes nothing, bit for bit");

    // A binding cap.
    auto o2 = scored();
    o2.sample_max_cycles = 40;
    ODIA::MassCalibration::Diagnostics d2;
    const auto m2 = ODIA::MassCalibration::calibrate(w.lib, w.run, o2, &d2);
    std::printf("   capped at 40 cycles:\n%s", ODIA::MassCalibration::report(m2, &d2).c_str());
    check(d2.sample_cycles_available == d0.sample_cycles_available,
          "the cap is drawn from the same reachable set");
    check(d2.sample_cycles > 0 && d2.sample_cycles <= 40,
          "at most 40 cycles visited (" + std::to_string(d2.sample_cycles) + ")");
    check(d2.spectra_decoded == 4 * d2.sample_cycles,
          "and only their spectra decoded (" + std::to_string(d2.spectra_decoded) + ")");
    check(d2.sample_entries_in_reach > 0 && d2.sample_entries_in_reach < 160 &&
            d2.sample_entries_probed == d2.sample_entries_in_reach,
          "entries no drawn cycle reaches are dropped (" +
            std::to_string(d2.sample_entries_in_reach) + " of 160 kept)");
    check(m2.fitted && std::abs(m2.ppmAt(700.0) - PLANTED_PPM) < 0.7,
          "and the capped probe still recovers the planted " + std::to_string(PLANTED_PPM) +
            " ppm (got " + std::to_string(m2.ppmAt(700.0)) + ")");

    // The precursor cap comes AFTER the cycle cap, and the visited cycles are
    // re-derived from what it leaves.
    auto o3 = scored();
    o3.sample_max_cycles = 40;
    o3.max_precursors = 10;
    ODIA::MassCalibration::Diagnostics d3;
    (void)ODIA::MassCalibration::collect(w.lib, w.run, o3, &d3);
    check(d3.sample_entries_in_reach == d2.sample_entries_in_reach &&
            d3.sample_entries_probed == 10,
          "max_precursors strides the entries the cycle cap kept (" +
            std::to_string(d3.sample_entries_probed) + " of " +
            std::to_string(d3.sample_entries_in_reach) + ")");
    check(d3.sample_cycles > 0 && d3.sample_cycles <= d2.sample_cycles,
          "and only cycles that still reach a probed entry are visited (" +
            std::to_string(d3.sample_cycles) + ")");
  }
} // namespace

int main()
{
  selection();
  World w;
  build(w);
  identity(w);
  probe(w);
  emptySample(w);
  cycleCap(w);
  std::printf("%s (%d failure%s)\n", failures ? "FAILED" : "PASSED", failures,
              failures == 1 ? "" : "s");
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
