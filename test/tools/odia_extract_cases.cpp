// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// ChromatogramExtractor against a scripted run: every spectrum, peak, mobility
// and isolation window written out in the test, so each case is exactly the
// hazard it names and nothing else. No data file, no reader, milliseconds.
//
// The suite had 44 tests and not one of them constructed a Chromatograms or
// called extract(). Four defects lived in that gap, and each case below is one
// of them:
//
//   invalid_mz    a transition whose product m/z cannot be represented is
//                 counted but never matched
//   aggregate     Sum against Max over several peaks in one tolerance
//   mobility      the per-precursor 1/K0 window, admitting and rejecting
//   im_gating     -no_ion_mobility must not also switch that window off
//   band_edge     a peak exactly on the boundary between two co-packed
//                 windows belongs to one of them, not to both
//   wide_csr      more than 2^32 chromatogram points, which the 32-bit CSR
//                 offsets refused outright (needs ~17.2 GiB)
//   plane_identity  the -debug_plane_identity counter fires where the residual
//                 denominators differ from the intensity plane (positive
//                 control) and not where they do not
//
// Usage: odia_extract_cases <case>

#include <odia/ChromatogramExtractor.h>
#include <odia/SpectrumSource.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

namespace
{
  const float NA = std::numeric_limits<float>::quiet_NaN();

  /// A run written out spectrum by spectrum.
  ///
  /// Spectra are added in acquisition order, which is the contract
  /// SpectrumSource states; a spectrum names the window it was acquired in, so
  /// several windows can share one cycle (and, as on IH1, one peak list).
  class ScriptedRun : public ODIA::SpectrumSource
  {
  public:
    std::size_t addWindow(double mz_low, double mz_high,
                          double im_low = -std::numeric_limits<double>::infinity(),
                          double im_high = std::numeric_limits<double>::infinity())
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

    void addPeak(std::size_t spectrum, double mz, float intensity, float im = NA)
    {
      auto& pk = peaks_[spectrum];
      pk.mz.push_back(mz);
      pk.intensity.push_back(intensity);
      // Either the run separates by mobility or it does not, so a spectrum
      // carries a mobility for every peak or for none.
      if (!std::isnan(im)) { pk.ion_mobility.push_back(im); }
    }

    /// Like addPeak, but the mobility is stored even when it is NaN: a
    /// spectrum that carries mobilities with one peak's missing, which
    /// addPeak cannot express.
    void addPeakIm(std::size_t spectrum, double mz, float intensity, float im)
    {
      auto& pk = peaks_[spectrum];
      pk.mz.push_back(mz);
      pk.intensity.push_back(intensity);
      pk.ion_mobility.push_back(im);
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
    /// @p mz in Th, @p im the library 1/K0 or NA for "the library has none".
    std::size_t addPrecursor(double mz, float im = NA)
    {
      auto& p = lib_.precursors();
      p.mz.push_back(ODIA::toFixed(mz));
      p.irt.push_back(0.0f);
      p.im.push_back(im);
      p.ccs.push_back(NA);
      p.charge.push_back(2);
      p.decoy.push_back(0);
      p.modified_sequence.push_back(0);
      p.protein_group.push_back(0);
      p.transition_begin.push_back(
        static_cast<std::uint32_t>(lib_.transitions().product_mz.size()));
      p.transition_count.push_back(0);
      return p.mz.size() - 1;
    }

    /// A product m/z of 0 (or anything toFixed refuses) is the unrepresentable
    /// case the library reports as invalidMzTransitionCount().
    std::size_t addTransition(double product_mz, float library_intensity = 1.0f)
    {
      auto& t = lib_.transitions();
      t.product_mz.push_back(ODIA::toFixed(product_mz));
      t.library_intensity.push_back(library_intensity);
      t.type.push_back(ODIA::FragmentType::Y);
      t.ordinal.push_back(4);
      t.charge.push_back(1);
      t.loss.push_back(ODIA::LossType::None);
      ++lib_.precursors().transition_count.back();
      return t.product_mz.size() - 1;
    }

    ODIA::Library& library() { return lib_; }

  private:
    ODIA::Library lib_;
  };

  int failures = 0;

  void check(bool ok, const std::string& what)
  {
    std::printf("%s %s\n", ok ? "  ok  " : "  FAIL", what.c_str());
    if (!ok) { ++failures; }
  }

  void checkNear(double got, double want, const std::string& what)
  {
    const bool ok = std::abs(got - want) <= 1e-3;
    std::printf("%s %s (want %.4g, got %.4g)\n", ok ? "  ok  " : "  FAIL",
                what.c_str(), want, got);
    if (!ok) { ++failures; }
  }

  /// The invariant every consumer of a Chromatograms depends on: each
  /// transition's points lie on the axis it names, and the CSR index covers the
  /// intensity array exactly once.
  void checkRepresentation(const ODIA::Chromatograms& x)
  {
    std::size_t running = 0;
    bool consistent = true, in_bounds = true;
    for (std::size_t j = 0; j < x.count.size(); ++j)
    {
      if (x.begin[j] != running) { consistent = false; }
      running += x.count[j];
      if (x.count[j] == 0) { continue; }
      const std::size_t a = x.axisOf(j);
      if (a >= x.axes.size() ||
          std::size_t(x.axisBegin(j)) + x.count[j] > x.axes[a].size())
      {
        in_bounds = false;
      }
    }
    check(consistent && running == x.intensity.size(), "CSR index covers the intensities");
    check(in_bounds, "every transition's points lie on the axis it names");
  }

  double sumTrace(const ODIA::Chromatograms& x, std::uint32_t transition)
  {
    double total = 0.0;
    for (std::uint32_t j = 0; j < x.count[transition]; ++j)
    {
      total += x.intensity[x.begin[transition] + j];
    }
    return total;
  }

  ODIA::ChromatogramExtractor::Options plainOptions()
  {
    ODIA::ChromatogramExtractor::Options opt;
    opt.threads = 1;
    opt.progress_every = 0;
    return opt;
  }

  // ----------------------------------------------------------------- cases

  /// A transition whose product m/z could not be represented.
  ///
  /// It cannot be matched -- there is no mass to match against -- so it must
  /// get no points. It used to be COUNTED in the first pass and skipped in the
  /// second, which left it with a full-length point run it never filled, the
  /// default axis_of = 0, and therefore a read off the end of window 0's axis
  /// whenever window 0 was shorter than its own. So the two windows here have
  /// deliberately different lengths.
  void caseInvalidMz()
  {
    ScriptedRun run;
    const auto w0 = run.addWindow(500.0, 510.0);
    const auto w1 = run.addWindow(520.0, 530.0);
    for (int c = 0; c < 5; ++c)
    {
      if (c < 2)                                   // window 0 stops early: 2 cycles
      {
        const auto s = run.addSpectrum(w0, 100.0 + 10.0 * c);
        run.addPeak(s, 400.0, 3.0f);
      }
      const auto s = run.addSpectrum(w1, 100.0 + 10.0 * c + 1.0);   // 5 cycles
      run.addPeak(s, 400.0, 6.0f);
    }

    ScriptedLibrary lib;
    lib.addPrecursor(505.0);
    const auto short_window_transition = lib.addTransition(400.0);
    lib.addPrecursor(525.0);
    const auto valid = lib.addTransition(400.0);
    const auto invalid = lib.addTransition(0.0);   // unrepresentable
    const auto after = lib.addTransition(400.0);

    check(lib.library().invalidMzTransitionCount() == 1,
          "the library reports the unrepresentable transition");

    const auto x = ODIA::ChromatogramExtractor::extract(lib.library(), run, plainOptions());
    checkRepresentation(x);

    check(x.count[short_window_transition] == 2, "window 0 has two cycles");
    check(x.count[valid] == 5, "window 1 has five cycles");
    check(x.count[invalid] == 0, "the unrepresentable transition gets no points");
    check(x.count[after] == 5, "its neighbours are unaffected");
    checkNear(sumTrace(x, valid), 5 * 6.0, "the valid transition is extracted");
    checkNear(sumTrace(x, after), 5 * 6.0, "the transition after it is extracted");
  }

  /// Sum against Max over several peaks inside one transition's tolerance.
  ///
  /// The default is Sum and has been since the aggregate became an option: a
  /// diaPASEF frame's peak array is the concatenation of 600-810 mobility
  /// scans, so peaks inside one m/z tolerance are the same ion at many mobility
  /// steps and Max returns the interference envelope instead of the ion. Both
  /// arms are asserted here because a default that silently reverted, or an
  /// option that stopped being read, looks like nothing at all downstream.
  void caseAggregate()
  {
    ScriptedRun run;
    const auto w = run.addWindow(500.0, 510.0);
    for (int c = 0; c < 3; ++c)
    {
      const auto s = run.addSpectrum(w, 100.0 + 10.0 * c);
      // +/- 10 ppm of 400 Th is +/- 0.004: the first two are inside one
      // tolerance, the third is outside it and must never be counted.
      run.addPeak(s, 400.0, 10.0f);
      run.addPeak(s, 400.001, 4.0f);
      run.addPeak(s, 400.010, 99.0f);
    }

    ScriptedLibrary lib;
    lib.addPrecursor(505.0);
    const auto tr = lib.addTransition(400.0);

    auto opt = plainOptions();
    opt.aggregate = ODIA::ChromatogramExtractor::Options::Aggregate::Sum;
    const auto summed = ODIA::ChromatogramExtractor::extract(lib.library(), run, opt);
    checkRepresentation(summed);
    check(summed.count[tr] == 3, "three cycles");
    checkNear(sumTrace(summed, tr), 3 * 14.0, "Sum integrates both peaks in the tolerance");

    opt.aggregate = ODIA::ChromatogramExtractor::Options::Aggregate::Max;
    const auto maxed = ODIA::ChromatogramExtractor::extract(lib.library(), run, opt);
    checkNear(sumTrace(maxed, tr), 3 * 10.0, "Max takes the larger of the two");

    const auto& defaulted = ODIA::ChromatogramExtractor::Options{}.aggregate;
    check(defaulted == ODIA::ChromatogramExtractor::Options::Aggregate::Sum,
          "Sum is the default");
  }

  /// The per-precursor ion-mobility window, admitting and rejecting.
  ///
  /// Three peaks at the same m/z, one per mobility position: the precursor's
  /// own, one elsewhere in the frame's band, and one outside the band
  /// altogether. Which of them survives says which of the two mobility tests
  /// did the work, and the case is run four ways so that neither test can be
  /// removed without a failure.
  void caseMobility()
  {
    ScriptedRun run;
    const auto w = run.addWindow(500.0, 510.0, 0.80, 1.20);
    for (int c = 0; c < 3; ++c)
    {
      const auto s = run.addSpectrum(w, 100.0 + 10.0 * c);
      // Descending 1/K0, which is the order a merged TIMS frame arrives in.
      run.addPeak(s, 400.0, 3.0f, 1.30f);    // outside the frame's band
      run.addPeak(s, 400.0, 7.0f, 1.10f);    // inside the band, wrong precursor
      run.addPeak(s, 400.0, 10.0f, 0.91f);   // this precursor, 0.01 off its library 1/K0
    }

    ScriptedLibrary lib;
    lib.addPrecursor(505.0, 0.90f);
    const auto known = lib.addTransition(400.0);
    lib.addPrecursor(505.0, NA);            // library carries no 1/K0
    const auto unknown = lib.addTransition(400.0);

    auto opt = plainOptions();
    const auto both = ODIA::ChromatogramExtractor::extract(lib.library(), run, opt);
    checkRepresentation(both);
    checkNear(sumTrace(both, known), 3 * 10.0,
              "both tests: only the precursor's own mobility survives");
    checkNear(sumTrace(both, unknown), 3 * 17.0,
              "no library 1/K0 disables only the per-precursor test, not the band");

    opt.precursor_im_window = 0.0;
    const auto band_only = ODIA::ChromatogramExtractor::extract(lib.library(), run, opt);
    checkNear(sumTrace(band_only, known), 3 * 17.0,
              "band only: the frame's band admits the whole same-window axis");

    // What the per-precursor window is worth on this case: 17 -> 10, i.e. the
    // same-window interference is 0.7x the signal and the band cannot see it.
    opt.precursor_im_window = 0.025;
    opt.use_ion_mobility = true;
    ScriptedRun no_im;
    const auto w2 = no_im.addWindow(500.0, 510.0, 0.80, 1.20);
    for (int c = 0; c < 3; ++c)
    {
      const auto s = no_im.addSpectrum(w2, 100.0 + 10.0 * c);
      no_im.addPeak(s, 400.0, 10.0f);       // a run that reports no mobility at all
    }
    const auto without = ODIA::ChromatogramExtractor::extract(lib.library(), no_im, opt);
    checkNear(sumTrace(without, known), 3 * 10.0,
              "a run without mobility is not filtered by mobility");
  }

  /// `use_ion_mobility` gates the frame's band and NOTHING else.
  ///
  /// The two tests are distinct -- the header says so, and they answer
  /// different questions: the band separates co-packed windows, the
  /// per-precursor window separates a peptide from its same-window
  /// neighbours. Switching the band off (-no_ion_mobility) also switched the
  /// per-precursor window off, because the peak's mobility was only READ
  /// inside the band's branch and stayed NaN otherwise -- and a NaN mobility
  /// skips the per-precursor test by the "absent information is not evidence"
  /// rule that is meant for a run that carries no mobility at all.
  ///
  /// That is not a small thing: -no_ion_mobility is the control arm used to
  /// measure what mobility filtering buys, and it was silently measuring both
  /// filters against neither.
  void caseImGating()
  {
    ScriptedRun run;
    const auto w = run.addWindow(500.0, 510.0, 0.80, 1.20);
    for (int c = 0; c < 3; ++c)
    {
      const auto s = run.addSpectrum(w, 100.0 + 10.0 * c);
      run.addPeak(s, 400.0, 3.0f, 1.30f);    // outside the frame's band
      run.addPeak(s, 400.0, 7.0f, 1.10f);    // inside the band, wrong precursor
      run.addPeak(s, 400.0, 10.0f, 0.91f);   // this precursor
    }

    ScriptedLibrary lib;
    lib.addPrecursor(505.0, 0.90f);
    const auto tr = lib.addTransition(400.0);

    auto opt = plainOptions();
    opt.use_ion_mobility = false;            // band off
    opt.precursor_im_window = 0.025;         // per-precursor window still on
    const auto x = ODIA::ChromatogramExtractor::extract(lib.library(), run, opt);
    checkRepresentation(x);
    // Only the peak at 0.91 survives: 1.10 and 1.30 are 0.2 and 0.4 from the
    // precursor's 0.90, so the per-precursor window rejects both whether or
    // not the frame's band would have. Before the fix this returned 20 per
    // cycle -- everything at that m/z, at every mobility in the frame.
    checkNear(sumTrace(x, tr), 3 * 10.0,
              "the band is off, the per-precursor window is not");

    // And with both off, everything at that m/z arrives: 20 per cycle. This is
    // what the case above must NOT return.
    opt.precursor_im_window = 0.0;
    const auto none = ODIA::ChromatogramExtractor::extract(lib.library(), run, opt);
    checkNear(sumTrace(none, tr), 3 * 20.0, "with both off, nothing is filtered");
  }

  /// A peak exactly on the boundary between two co-packed windows.
  ///
  /// IH1 packs two isolation windows into one frame and the reader hands both
  /// the SAME merged peak list, separated only by a mobility band -- and that
  /// band is derived, not stated: the split is the midpoint between the two
  /// windows' mobility positions, so the bands are adjacent and share their
  /// boundary exactly. A test inclusive at both ends therefore puts a peak
  /// sitting on the midpoint into BOTH windows, where under Sum it is
  /// integrated twice into two different precursors' traces.
  ///
  /// Measure zero on real data and impossible to see downstream, which is
  /// exactly why it needs a test rather than a measurement: the band is
  /// half-open, [im_low, im_high), so every peak of a frame lands in one
  /// window.
  void caseBandEdge()
  {
    ScriptedRun run;
    const auto w0 = run.addWindow(500.0, 510.0, 0.80, 1.00);
    const auto w1 = run.addWindow(510.001, 520.0, 1.00, 1.20);
    for (int c = 0; c < 3; ++c)
    {
      // One frame, two windows, one peak list -- as the reader serves IH1.
      for (const auto w : {w0, w1})
      {
        const auto s = run.addSpectrum(w, 100.0 + 10.0 * c);
        run.addPeak(s, 400.0, 3.0f, 1.10f);    // window 1's band
        run.addPeak(s, 400.0, 5.0f, 1.00f);    // exactly on the boundary
        run.addPeak(s, 400.0, 2.0f, 0.90f);    // window 0's band
      }
    }

    ScriptedLibrary lib;
    lib.addPrecursor(505.0);                   // window 0, no library 1/K0
    const auto lower = lib.addTransition(400.0);
    lib.addPrecursor(515.0);                   // window 1
    const auto upper = lib.addTransition(400.0);

    auto opt = plainOptions();
    opt.precursor_im_window = 0.0;             // the frame's band is the only test
    const auto x = ODIA::ChromatogramExtractor::extract(lib.library(), run, opt);
    checkRepresentation(x);

    checkNear(sumTrace(x, lower), 3 * 2.0, "the lower window stops below the boundary");
    checkNear(sumTrace(x, upper), 3 * 8.0, "the upper window starts at the boundary");
    checkNear(sumTrace(x, lower) + sumTrace(x, upper), 3 * 10.0,
              "every peak of the frame is counted exactly once");
  }

  /// More than 2^32 chromatogram points.
  ///
  /// The CSR offsets used to be uint32 over a flat point array, so extraction
  /// threw above 2^32 points -- 266,664 precursors at IH1's 12 transitions x
  /// 1,342 cycles, against the 4,255,113 precursors Phase 1's own human
  /// library holds. ODIA could not extract against the library it generated,
  /// and the failure was a hard throw rather than a wrong number.
  ///
  /// Reproduced at the smallest geometry that crosses the line: 1,048,577
  /// transitions of 4,096 cycles each, so the last transition's `begin` is
  /// exactly 2^32 -- the first offset a uint32 cannot hold, and the one it
  /// would silently wrap to 0. The peak below is placed on that transition, so
  /// this asserts that the points ABOVE the old ceiling are written and read
  /// back at the right place, not merely that nothing threw.
  ///
  /// Costs ~17.2 GiB, because the point array is preallocated as transitions x
  /// cycles. That is the whole remaining limit and the reason the header says
  /// running at proteome scale needs chunking, not a wider type. The test is
  /// registered only where the memory exists.
  void caseWideCsr()
  {
    constexpr std::uint32_t CYCLES = 4096;
    // ceil(2^32 / CYCLES) + 1, so the last transition starts at or above 2^32.
    constexpr std::size_t TRANSITIONS = (std::size_t(1) << 32) / CYCLES + 1;

    ScriptedRun run;
    const auto w = run.addWindow(500.0, 510.0);

    ScriptedLibrary lib;
    // One transition per precursor, which also puts the precursor count 3.9x
    // over the old 266,664 ceiling.
    for (std::size_t k = 0; k < TRANSITIONS; ++k)
    {
      lib.addPrecursor(505.0);
      lib.addTransition(200.0 + double(k) * 0.01);
    }
    const std::uint32_t last = static_cast<std::uint32_t>(TRANSITIONS - 1);
    const double last_mz = 200.0 + double(TRANSITIONS - 1) * 0.01;

    for (std::uint32_t c = 0; c < CYCLES; ++c)
    {
      const auto s = run.addSpectrum(w, 100.0 + double(c));
      run.addPeak(s, last_mz, 7.0f);
    }

    const auto x = ODIA::ChromatogramExtractor::extract(lib.library(), run, plainOptions());
    checkRepresentation(x);

    const std::uint64_t ceiling = std::uint64_t(std::numeric_limits<std::uint32_t>::max()) + 1;
    check(x.points() > ceiling, "more than 2^32 points were extracted");
    check(x.begin[last] >= ceiling,
          "the last transition starts past what a 32-bit offset can hold");
    check(x.count[last] == CYCLES, "and it has one point per cycle");
    checkNear(sumTrace(x, last), double(CYCLES) * 7.0,
              "every point above the ceiling holds the intensity written to it");
  }

  /// A sink that keeps a private copy of everything it is handed.
  ///
  /// The storage behind a `PrecursorChromatogram` is valid only for the
  /// duration of `accept`, so copying is the only way to check afterwards what
  /// was handed over -- and a sink that reads it later would be exactly the
  /// defect this checks against.
  class RecordingSink final : public ODIA::ChromatogramSink
  {
  public:
    struct Trace
    {
      std::uint32_t precursor = 0, transition_begin = 0;
      std::uint32_t axis = 0, axis_begin = 0, cycles = 0;
      std::vector<float> rt;
      std::vector<std::vector<float>> points;   ///< one per transition
    };

    void accept(const ODIA::PrecursorChromatogram& c) override
    {
      Trace t;
      t.precursor = c.precursor;
      t.transition_begin = c.transition_begin;
      t.axis = c.axis;
      t.axis_begin = c.axis_begin;
      t.cycles = c.cycles;
      t.rt.assign(c.rt, c.rt + c.cycles);
      for (std::uint32_t k = 0; k < c.transition_count; ++k)
      {
        const std::uint32_t n = c.pointCount(k);
        t.points.emplace_back(n ? c.trace(k) : nullptr, n ? c.trace(k) + n : nullptr);
      }
      traces.push_back(std::move(t));
    }

    std::vector<Trace> traces;
  };

  /// The sliding window: a precursor is allocated when the pass reaches the
  /// start of its retention-time window and freed when the pass leaves it.
  ///
  /// Three claims, each of which the old extractor could not have satisfied
  /// because it allocated the whole library before reading a peak:
  ///
  ///   * only the precursors whose windows OVERLAP are ever resident at once;
  ///   * what a streaming consumer is handed is exactly what the flat array
  ///     holds -- same points, same axis, same cycle range;
  ///   * capping the live set and chunking the library changes nothing about
  ///     the answer, only about how many passes it takes.
  ///
  /// The geometry is deliberately extreme -- twenty precursors tiled across the
  /// gradient with windows a fifteenth of it -- because at a 600 s window on a
  /// 1,859 s run about 65% of a library is live at once and a test built at
  /// that ratio would pass with the sliding window removed.
  void caseSlidingWindow()
  {
    // Long enough to span many match batches. The sliding window can only move
    // between them, so a run shorter than one batch has everything live at once
    // whatever the windows say -- which is a real property of the design and
    // the reason this is not a four-spectrum test like the others.
    constexpr std::uint32_t CYCLES = 2048;
    constexpr std::uint32_t PRECURSORS = 32;
    constexpr double STEP = 1.0;               // seconds per cycle

    ScriptedRun run;
    const auto w = run.addWindow(500.0, 510.0);

    ScriptedLibrary lib;
    std::vector<double> product;
    for (std::uint32_t i = 0; i < PRECURSORS; ++i)
    {
      const std::size_t at = lib.addPrecursor(505.0);
      lib.library().precursors().irt[at] =
        static_cast<float>(STEP * double(CYCLES) * double(i) / double(PRECURSORS));
      product.push_back(200.0 + double(i) * 0.5);
      product.push_back(600.0 + double(i) * 0.5);
      lib.addTransition(product[2 * i]);
      lib.addTransition(product[2 * i + 1]);
    }

    // Every transition is present in every spectrum, with an intensity that
    // identifies the cycle. A precursor's trace is then a slice of the cycle
    // index, so a chromatogram assembled at the wrong offset is visible rather
    // than merely smaller.
    for (std::uint32_t c = 0; c < CYCLES; ++c)
    {
      const auto s = run.addSpectrum(w, STEP * double(c));
      for (const double mz : product) { run.addPeak(s, mz, float(c) + 1.0f); }
    }

    auto opt = plainOptions();
    opt.irt_slope = 1.0;                       // the library's iRT is run seconds
    opt.irt_intercept = 0.0;
    opt.rt_window_seconds = 100.0;             // +/- 100 cycles, so windows overlap

    ODIA::ChromatogramExtractor::Stats flat_stats;
    const auto x = ODIA::ChromatogramExtractor::extract(lib.library(), run, opt,
                                                        &flat_stats);
    checkRepresentation(x);

    // Every precursor is extracted, and over a window far shorter than the run.
    check(flat_stats.precursors_extracted == PRECURSORS,
          "every precursor was extracted");
    check(x.count[0] > 0 && x.count[0] < CYCLES,
          "a precursor's window is a slice of the run, not the whole of it");

    // The claim the design rests on. Windows are 5-6 cycles wide and start 2
    // cycles apart, so a handful overlap -- never all twenty.
    check(flat_stats.peak_live_precursors < PRECURSORS,
          "not every precursor was resident at once");
    check(flat_stats.peak_live_points < std::uint64_t(flat_stats.points),
          "peak resident points are below what the whole library would need");
    std::printf("       peak live %zu of %u precursors, %llu of %zu points\n",
                flat_stats.peak_live_precursors, PRECURSORS,
                static_cast<unsigned long long>(flat_stats.peak_live_points),
                flat_stats.points);

    // The value, not just the shape: the first live cycle of precursor 5 is
    // where its window starts, and its intensity is that cycle's own number.
    {
      const std::uint32_t tr = lib.library().precursors().transition_begin[5];
      const std::uint32_t first = x.axisBegin(tr);
      checkNear(x.intensity[x.begin[tr]], double(first) + 1.0,
                "a precursor's first point is the cycle its window starts at");
      checkNear(x.retentionTime(tr, 0), STEP * double(first),
                "and it carries that cycle's retention time");
    }

    // What a streaming consumer sees is what the flat array holds.
    RecordingSink recorded;
    ODIA::ChromatogramExtractor::Stats stream_stats;
    ODIA::ChromatogramExtractor::extract(lib.library(), run, opt, recorded, &stream_stats);

    std::vector<int> seen(PRECURSORS, 0);
    bool same = true, ordered = true;
    float previous_rt = -1.0f;
    for (const auto& t : recorded.traces)
    {
      if (t.precursor < PRECURSORS) { ++seen[t.precursor]; }
      if (t.cycles == 0) { continue; }
      // Handed over in the order the pass finishes them, which is the order
      // their windows end.
      if (!t.rt.empty())
      {
        if (t.rt.back() < previous_rt) { ordered = false; }
        previous_rt = t.rt.back();
      }
      for (std::size_t k = 0; k < t.points.size(); ++k)
      {
        const std::uint32_t tr = t.transition_begin + static_cast<std::uint32_t>(k);
        if (x.count[tr] != t.points[k].size()) { same = false; continue; }
        if (x.axisOf(tr) != t.axis || x.axisBegin(tr) != t.axis_begin) { same = false; }
        for (std::size_t j = 0; j < t.points[k].size(); ++j)
        {
          if (x.intensity[x.begin[tr] + j] != t.points[k][j]) { same = false; }
        }
      }
    }
    check(same, "the streaming sink is handed exactly what the flat array holds");
    check(ordered, "precursors are handed over in the order the pass leaves them");
    check(std::count(seen.begin(), seen.end(), 1) == int(PRECURSORS),
          "every precursor is handed over exactly once");

    // Capping the live set splits the library into chunks. The answer must not
    // move; only the number of passes does.
    auto capped = opt;
    // Below what the windows overlap by on their own, or the cap is not a cap
    // and the test asserts nothing -- exactly the way a threshold set above the
    // data silently disables the branch it guards.
    capped.max_live_precursors = 2;
    ODIA::ChromatogramExtractor::Stats chunk_stats;
    const auto y = ODIA::ChromatogramExtractor::extract(lib.library(), run, capped,
                                                        &chunk_stats);
    checkRepresentation(y);
    check(chunk_stats.chunks > 1, "the cap split the library into chunks");
    // Not "<= the cap": the cap partitions by retention-time overlap, and a
    // batch's worth of early activation can still put one or two more live than
    // the partition says. The claim that holds is that it is strictly tighter
    // than the uncapped run, which is what the mechanism is for.
    check(chunk_stats.peak_live_precursors < flat_stats.peak_live_precursors,
          "and the live set is tighter than without it");
    check(chunk_stats.spectra_decoded >= flat_stats.spectra_decoded,
          "which is paid for in spectra decoded more than once");
    check(y.count == x.count && y.precursor_axis == x.precursor_axis &&
          y.precursor_axis_begin == x.precursor_axis_begin &&
          y.precursor_cycles == x.precursor_cycles && y.intensity == x.intensity,
          "chunking changes how many passes it takes, not the answer");
    std::printf("       %zu chunks, %zu spectra decoded against %zu in the run\n",
                chunk_stats.chunks, chunk_stats.spectra_decoded, chunk_stats.spectra_read);

    // The per-chunk resource ledger on the multi-chunk path: each chunk adds
    // its index bracket and its snapshot exactly once, and no bracket starts
    // before the previous one ended.
    auto labelled = capped;
    labelled.stage_label = "test";
    ODIA::ChromatogramExtractor::Stats ledger_stats;
    const auto z = ODIA::ChromatogramExtractor::extract(lib.library(), run, labelled,
                                                        &ledger_stats);
    check(z.intensity == x.intensity, "the stage ledger changes nothing extracted");
    check(ledger_stats.chunks > 1 && ledger_stats.index.brackets == ledger_stats.chunks,
          "one index bracket per chunk, summed over the chunks");
    check(ledger_stats.probe.brackets == ledger_stats.chunks,
          "one memory snapshot per chunk, charged to the probe stage");
    check(ledger_stats.bracket_overlaps == 0, "no bracket overlaps another");
  }

  /// -debug_plane_identity: the comparator has to be able to fire.
  ///
  /// A count of zero on the fixture is only evidence for aliasing the residual
  /// denominators to the intensity plane if the same counter returns nonzero
  /// where they DO differ. Under Sum they differ exactly where a contributing
  /// peak was non-positive (ppm and im) or had no finite 1/K0 (im only); under
  /// Max wherever two peaks share a cell (base keeps one, the dens sum both).
  /// Each cell below is one of those, and the expected counts are exact.
  void casePlaneIdentity()
  {
    ScriptedRun run;
    const auto w = run.addWindow(500.0, 510.0);
    for (int c = 0; c < 3; ++c)
    {
      const auto s = run.addSpectrum(w, 100.0 + 10.0 * c);
      // Ascending m/z, as a spectrum arrives.
      if (c == 1)                                                // C: a negative peak
      {
        run.addPeakIm(s, 300.0, -3.0f, 1.0f);
        run.addPeakIm(s, 300.0, 10.0f, 1.0f);
      }
      run.addPeakIm(s, 400.0, 10.0f, 1.0f);                      // A: clean
      if (c == 2) { run.addPeakIm(s, 400.001, 4.0f, 1.0f); }    // A: two peaks, one cell
      if (c == 0) { run.addPeakIm(s, 450.0, 5.0f, NA); }         // B: no 1/K0
      else { run.addPeakIm(s, 450.0, 6.0f, 1.0f); }              // B: clean
    }
    ScriptedLibrary lib;
    lib.addPrecursor(505.0, 1.0f);
    lib.addTransition(400.0);
    lib.addTransition(450.0);
    lib.addTransition(300.0);

    auto opt = plainOptions();
    opt.collect_mass_residuals = true;
    opt.collect_im_residuals = true;
    const auto off = ODIA::ChromatogramExtractor::extract(lib.library(), run, opt);

    opt.count_plane_identity = true;
    ODIA::ChromatogramExtractor::Stats sum;
    const auto on = ODIA::ChromatogramExtractor::extract(lib.library(), run, opt, &sum);
    check(on.intensity == off.intensity, "counting changes nothing extracted");
    check(sum.plane_cells == 9, "3 transitions x 3 cycles were compared");
    // Sum: C@1 (-3 in base, not in either den), B@0 (no 1/K0: im_den only).
    check(sum.ppm_den_differs == 1, "Sum: ppm_den differs at the negative peak only");
    check(sum.im_den_differs == 2, "Sum: im_den differs at the negative and the NaN-1/K0 peak");
    check(sum.im_den_zero == 1, "Sum: one cell had no finite 1/K0 at all");
    std::printf("       Sum: cells %llu, ppm %llu, im %llu (zero %llu)\n",
                static_cast<unsigned long long>(sum.plane_cells),
                static_cast<unsigned long long>(sum.ppm_den_differs),
                static_cast<unsigned long long>(sum.im_den_differs),
                static_cast<unsigned long long>(sum.im_den_zero));

    opt.aggregate = ODIA::ChromatogramExtractor::Options::Aggregate::Max;
    ODIA::ChromatogramExtractor::Stats mx;
    ODIA::ChromatogramExtractor::extract(lib.library(), run, opt, &mx);
    // Max: A@2 (base 10, dens 14). C@1 now agrees: the -3 never raises the
    // max, and the dens skip it. B@0 is still im only.
    check(mx.ppm_den_differs == 1, "Max: ppm_den differs where two peaks share a cell");
    check(mx.im_den_differs == 2, "Max: im_den differs there and at the NaN-1/K0 peak");

    // Negative control: the clean transition alone, under Sum, counts zero.
    ScriptedLibrary clean;
    clean.addPrecursor(505.0, 1.0f);
    clean.addTransition(450.0);
    opt.aggregate = ODIA::ChromatogramExtractor::Options::Aggregate::Sum;
    ScriptedRun run2;
    const auto w2 = run2.addWindow(500.0, 510.0);
    for (int c = 0; c < 3; ++c)
    {
      const auto s = run2.addSpectrum(w2, 100.0 + 10.0 * c);
      run2.addPeakIm(s, 450.0, 6.0f, 1.0f);
    }
    ODIA::ChromatogramExtractor::Stats zero;
    ODIA::ChromatogramExtractor::extract(clean.library(), run2, opt, &zero);
    check(zero.plane_cells == 3 && zero.ppm_den_differs == 0 && zero.im_den_differs == 0,
          "a clean run counts zero");
  }
}

int main(int argc, char** argv)
{
  const std::string which = argc > 1 ? argv[1] : "";
  std::printf("case: %s\n", which.c_str());

  if (which == "invalid_mz") { caseInvalidMz(); }
  else if (which == "aggregate") { caseAggregate(); }
  else if (which == "mobility") { caseMobility(); }
  else if (which == "im_gating") { caseImGating(); }
  else if (which == "band_edge") { caseBandEdge(); }
  else if (which == "wide_csr") { caseWideCsr(); }
  else if (which == "sliding") { caseSlidingWindow(); }
  else if (which == "plane_identity") { casePlaneIdentity(); }
  else
  {
    std::fprintf(stderr,
                 "usage: odia_extract_cases "
                 "<invalid_mz|aggregate|mobility|im_gating|band_edge|wide_csr|sliding|"
                 "plane_identity>\n");
    return 2;
  }

  std::printf(failures ? "FAILED (%d)\n" : "ok (%d)\n", failures);
  return failures ? 1 : 0;
}
