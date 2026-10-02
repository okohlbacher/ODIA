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
//   chunk_invariant  chunked extraction must hand the sink the unchunked
//                 stream: the same cells (a float chunk bound skipped edge
//                 frames), the same order (Gate C calibrates on arrival
//                 order), and no more live than the cap
//   cap_zero_rows an explicit cap survives a byte budget that binds nothing
//   chunk_matrix  every decode block x chunk count x thread count against an
//                 INDEPENDENT brute-force oracle, all five planes bitwise, the
//                 stream order, restricted RT bounds on both sides of float
//                 rounding, chunks cut between co-packed windows, empty
//                 spectra, missing acquisitions, no-window / outside-run /
//                 zero-transition / zero-valid precursors, both caps at once
//   no_assignments  nothing assigned: every precursor still handed over, empty
//   plane_identity  the -debug_plane_identity counter fires where the residual
//                 denominators differ from the intensity plane (positive
//                 control) and not where they do not
//   instr_matrix  the stage ledger and plane counter on the chunked path: the
//                 stream unchanged, every chunk/block/batch/part bracketed once
//
// Usage: odia_extract_cases <case>

#include <odia/ChromatogramExtractor.h>
#include <odia/SpectrumSource.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <limits>
#include <set>
#include <string>
#include <tuple>
#include <utility>
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

    /// What a spectrum holds, read directly -- for an oracle that must not go
    /// through the extractor's own decode path.
    const ODIA::SpectrumPeaks& peakList(std::size_t spectrum) const { return peaks_[spectrum]; }

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
      std::uint32_t transition_count = 0;
      std::vector<float> rt;
      std::vector<std::vector<float>> points;   ///< one per transition
      /// The residual planes, one per transition, empty when not collected.
      std::vector<std::vector<float>> ppm_num, ppm_den, im_num, im_den;
    };

    void accept(const ODIA::PrecursorChromatogram& c) override
    {
      Trace t;
      t.precursor = c.precursor;
      t.transition_begin = c.transition_begin;
      t.axis = c.axis;
      t.axis_begin = c.axis_begin;
      t.cycles = c.cycles;
      t.transition_count = c.transition_count;
      if (c.cycles) { t.rt.assign(c.rt, c.rt + c.cycles); }
      const auto plane = [&](const float* base, std::uint32_t k,
                             std::vector<std::vector<float>>& into) {
        const std::uint32_t n = c.pointCount(k);
        if (base == nullptr || n == 0) { into.emplace_back(); return; }
        into.emplace_back(base + c.offset[k], base + c.offset[k] + n);
      };
      for (std::uint32_t k = 0; k < c.transition_count; ++k)
      {
        const std::uint32_t n = c.pointCount(k);
        t.points.emplace_back(n ? c.trace(k) : nullptr, n ? c.trace(k) + n : nullptr);
        plane(c.ppm_num, k, t.ppm_num);
        plane(c.ppm_den, k, t.ppm_den);
        plane(c.im_num, k, t.im_num);
        plane(c.im_den, k, t.im_den);
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
    // The planner counts liveness on the match-batch grid the pass moves on,
    // so the cap is exact. (It used to count retention-time overlap, and a
    // batch's worth of early activation put more live than the cap allowed --
    // 734,070 under a cap of 400,000 on the PXD fixture.)
    check(chunk_stats.peak_live_precursors <= capped.max_live_precursors,
          "and no more precursors are live than the cap allows");
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

  /// Chunking is invisible to the sink: the same traces, in the same order.
  ///
  /// `sliding` asserts this too, and passed for as long as it was false,
  /// because its retention times were whole seconds -- exact in float. The
  /// extractor kept its cycle axis in float and looked each chunk's spectrum
  /// range up in the run's DOUBLE times with those float keys, so whenever a
  /// time rounded the wrong way the chunk's first or last frame was never
  /// read and its cells stayed zero. A real run's times are never whole: on
  /// the PXD fixture 2,992 of 6,020 rounded up and 3,028 down, and every chunk
  /// edge lost its frame with probability one half. Here the times are
  /// 1000.37 + 0.1 c, which is the same hazard.
  ///
  /// Two further things are asserted, because the sink is not order-free and
  /// the cap is a memory promise:
  ///
  ///   * ORDER. Gate C sets its threshold from the first decoys to ARRIVE, so
  ///     a chunk plan that hands precursors over in a different order hands
  ///     the scorer a different calibration sample (13 chunks moved tau from
  ///     6.61 to 6.12 and peak groups by +17% on the fixture). The stream must
  ///     be the unchunked stream, element for element.
  ///   * CAP. No more precursors live than the cap: the planner must count on
  ///     the batch grid the pass moves on, not on retention-time overlap.
  ///
  /// The geometry has the features that break each: two co-packed windows
  /// sharing every frame time (they lose a frame together); windows clipped at
  /// the run's ends and at a 100 s gap in the middle, so many precursors start
  /// on one cycle and end on different ones and start order is not hand-over
  /// order; that gap falling inside one match batch, which then holds the
  /// precursors on both sides of it live at once although their retention
  /// times never overlap -- as the PXD fixture's 60 s slices do; and two caps
  /// that give two different chunk plans.
  void caseChunkInvariant()
  {
    constexpr std::uint32_t CYCLES = 1500, GAP_AT = 750;
    constexpr std::uint32_t PRECURSORS = 64;
    constexpr double T0 = 1000.37, STEP = 0.1, HALF = 25.0, GAP = 100.0;

    ScriptedRun run;
    const auto w0 = run.addWindow(500.0, 510.0);
    const auto w1 = run.addWindow(510.001, 520.0);

    ScriptedLibrary lib;
    std::vector<double> product;
    const double span = STEP * double(CYCLES) + GAP + 2.0 * 15.0;
    for (std::uint32_t i = 0; i < PRECURSORS; ++i)
    {
      const std::size_t at = lib.addPrecursor(i % 2 ? 515.0 : 505.0);
      lib.library().precursors().irt[at] =
        static_cast<float>(T0 - 15.0 + span * double(i) / double(PRECURSORS - 1));
      for (int k = 0; k < 3; ++k)
      {
        product.push_back(200.0 + double(i) * 0.5 + 300.0 * k);
        lib.addTransition(product.back());
      }
    }

    // One frame per cycle carrying both windows at the SAME time, every
    // transition present with an intensity that names the cycle and window.
    for (std::uint32_t c = 0; c < CYCLES; ++c)
    {
      const double rt = T0 + STEP * double(c) + (c >= GAP_AT ? GAP : 0.0);
      for (const auto w : {w0, w1})
      {
        const auto s = run.addSpectrum(w, rt);
        for (const double mz : product)
        { run.addPeak(s, mz, float(c) + 1.0f + (w == w1 ? 0.5f : 0.0f)); }
      }
    }
    std::size_t rounds_up = 0, rounds_down = 0;
    for (const auto& s : run.spectra())
    {
      const double f = double(static_cast<float>(s.retention_time));
      rounds_up += f > s.retention_time;
      rounds_down += f < s.retention_time;
    }
    check(rounds_up > 0 && rounds_down > 0,
          "the run's times are not exact in float, in both directions");

    auto opt = plainOptions();
    opt.irt_slope = 1.0;
    opt.irt_intercept = 0.0;
    opt.rt_window_seconds = HALF;

    RecordingSink reference;
    ODIA::ChromatogramExtractor::Stats ref_stats;
    ODIA::ChromatogramExtractor::extract(lib.library(), run, opt, reference, &ref_stats);
    check(ref_stats.chunks == 1, "uncapped, the run is one chunk");

    const auto same_trace = [](const RecordingSink::Trace& a, const RecordingSink::Trace& b) {
      return a.precursor == b.precursor && a.transition_begin == b.transition_begin &&
             a.axis == b.axis && a.axis_begin == b.axis_begin && a.cycles == b.cycles &&
             a.rt == b.rt && a.points == b.points;
    };
    std::size_t plans = 0, last_chunks = 0;
    for (const std::size_t cap : {std::size_t(2), std::size_t(5)})
    {
      auto capped = opt;
      capped.max_live_precursors = cap;
      RecordingSink chunked;
      ODIA::ChromatogramExtractor::Stats st;
      ODIA::ChromatogramExtractor::extract(lib.library(), run, capped, chunked, &st);
      const std::string at = " (cap " + std::to_string(cap) + ", " +
                             std::to_string(st.chunks) + " chunks)";
      check(st.chunks > 1, "the cap splits the library" + at);
      if (st.chunks > 1 && st.chunks != last_chunks) { ++plans; }
      last_chunks = st.chunks;

      // Values, whatever the order: every precursor's trace, cell for cell.
      std::size_t differing_cells = 0, missing = 0;
      for (const auto& r : reference.traces)
      {
        const RecordingSink::Trace* found = nullptr;
        for (const auto& c : chunked.traces)
        { if (c.precursor == r.precursor) { found = &c; break; } }
        if (found == nullptr || found->points.size() != r.points.size() ||
            found->axis != r.axis || found->axis_begin != r.axis_begin)
        { ++missing; continue; }
        for (std::size_t k = 0; k < r.points.size(); ++k)
        {
          if (found->points[k].size() != r.points[k].size()) { ++missing; continue; }
          for (std::size_t j = 0; j < r.points[k].size(); ++j)
          { differing_cells += found->points[k][j] != r.points[k][j]; }
        }
      }
      std::printf("       cap %zu: %zu chunks, %zu cells differ, %zu traces unmatched, "
                  "peak live %zu\n", cap, st.chunks, differing_cells, missing,
                  st.peak_live_precursors);
      check(missing == 0 && differing_cells == 0,
            "every trace holds what the single pass put in it" + at);

      // Order: the stream itself.
      bool same_stream = chunked.traces.size() == reference.traces.size();
      for (std::size_t i = 0; same_stream && i < reference.traces.size(); ++i)
      { same_stream = same_trace(chunked.traces[i], reference.traces[i]); }
      check(same_stream, "the sink is handed the single pass's stream, in its order" + at);

      check(st.peak_live_precursors <= cap, "no more precursors live than the cap" + at);
    }
    check(plans == 2, "the two caps give two different chunk plans");

    // The decode block decides how many spectra are held, nothing else. It
    // used to set the match-batch boundaries too (a block of 200 matched
    // [128, 200) and [200, 328) where 256 matched [128, 256)), which moved the
    // hand-over order.
    std::size_t block_mismatches = 0;
    for (const std::size_t block : {1, 64, 127, 128, 129, 200, 256})
    {
      for (const std::size_t cap : {std::size_t(0), std::size_t(2), std::size_t(5)})
      {
        auto o = opt;
        o.decode_block = block;
        o.max_live_precursors = cap;
        RecordingSink got;
        ODIA::ChromatogramExtractor::Stats st;
        ODIA::ChromatogramExtractor::extract(lib.library(), run, o, got, &st);
        bool same = got.traces.size() == reference.traces.size();
        for (std::size_t i = 0; same && i < reference.traces.size(); ++i)
        { same = same_trace(got.traces[i], reference.traces[i]); }
        if (!same || (cap != 0 && st.peak_live_precursors > cap))
        {
          ++block_mismatches;
          std::printf("       decode_block %zu cap %zu (%zu chunks): stream differs\n",
                      block, cap, st.chunks);
        }
      }
    }
    check(block_mismatches == 0,
          "the stream is the same at every decode block, chunked or not");
  }

  /// An explicit cap survives a byte budget that constrains nothing.
  ///
  /// A precursor all of whose transitions lack a product m/z is extracted
  /// with zero rows, so it costs zero bytes; when every assigned precursor is
  /// like that the budget inverts to "0 precursors", which used to be read as
  /// unlimited -- and `min(explicit, 0)` erased the caller's cap with it. Two
  /// such precursors over the same cycles, `-max_live_precursors 1` and a
  /// budget: they must go through one at a time.
  void caseCapZeroRows()
  {
    ScriptedRun run;
    const auto w = run.addWindow(500.0, 510.0);
    for (int c = 0; c < 300; ++c) { run.addSpectrum(w, 100.0 + 0.7 * c); }

    ScriptedLibrary lib;
    for (int i = 0; i < 2; ++i)
    {
      lib.addPrecursor(505.0);
      lib.addTransition(0.0);                  // unrepresentable: no row
      lib.addTransition(0.0);
    }

    auto opt = plainOptions();
    opt.max_live_precursors = 1;
    opt.live_memory_budget_bytes = std::size_t(1) << 30;
    RecordingSink sink;
    ODIA::ChromatogramExtractor::Stats st;
    ODIA::ChromatogramExtractor::extract(lib.library(), run, opt, sink, &st);
    std::printf("       %zu chunks, peak live %zu, %zu traces; %s\n", st.chunks,
                st.peak_live_precursors, sink.traces.size(), st.live_budget_note.c_str());
    check(st.precursors_extracted == 2, "both zero-row precursors are extracted");
    check(st.peak_live_precursors <= 1,
          "the explicit cap holds when the byte budget binds nothing");
    check(st.chunks == 2, "so the two overlapping precursors take two chunks");
    check(sink.traces.size() == 2 && sink.traces[0].cycles == 300 &&
          sink.traces[1].cycles == 300 && sink.traces[0].precursor == 0 &&
          sink.traces[1].precursor == 1,
          "and both are handed over, whole, in the unchunked order");

    // Without an explicit cap the zero-cost budget leaves the run unbounded.
    opt.max_live_precursors = 0;
    ODIA::ChromatogramExtractor::Stats free_st;
    RecordingSink free_sink;
    ODIA::ChromatogramExtractor::extract(lib.library(), run, opt, free_sink, &free_st);
    check(free_st.chunks == 1 && free_st.peak_live_precursors == 2,
          "a budget alone that binds nothing leaves one chunk");
  }

  // ------------------------------------------------- the boundary matrix

  /// A scripted run that logs every decode call, so a test can see which
  /// spectrum ranges the chunks actually read.
  class LoggingRun final : public ScriptedRun
  {
  public:
    void peaks(std::size_t begin, std::size_t end,
               std::vector<ODIA::SpectrumPeaks>& out) override
    {
      calls.emplace_back(begin, end);
      ScriptedRun::peaks(begin, end, out);
    }
    std::vector<std::pair<std::size_t, std::size_t>> calls;
  };

  bool sameFloats(const std::vector<float>& a, const std::vector<float>& b)
  {
    return a.size() == b.size() &&
           (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(float)) == 0);
  }
  bool samePlanes(const std::vector<std::vector<float>>& a,
                  const std::vector<std::vector<float>>& b)
  {
    if (a.size() != b.size()) { return false; }
    for (std::size_t k = 0; k < a.size(); ++k) { if (!sameFloats(a[k], b[k])) { return false; } }
    return true;
  }
  /// Every field and every BIT of every plane.
  bool sameBits(const RecordingSink::Trace& a, const RecordingSink::Trace& b)
  {
    return a.precursor == b.precursor && a.transition_begin == b.transition_begin &&
           a.transition_count == b.transition_count && a.axis == b.axis &&
           a.axis_begin == b.axis_begin && a.cycles == b.cycles && sameFloats(a.rt, b.rt) &&
           samePlanes(a.points, b.points) && samePlanes(a.ppm_num, b.ppm_num) &&
           samePlanes(a.ppm_den, b.ppm_den) && samePlanes(a.im_num, b.im_num) &&
           samePlanes(a.im_den, b.im_den);
  }
  bool sameStream(const std::vector<RecordingSink::Trace>& a,
                  const std::vector<RecordingSink::Trace>& b)
  {
    if (a.size() != b.size()) { return false; }
    for (std::size_t i = 0; i < a.size(); ++i) { if (!sameBits(a[i], b[i])) { return false; } }
    return true;
  }

  /// What extraction must produce, worked out by brute force from the
  /// documented contract and NOT from the extractor's code paths: no m/z
  /// index, no sliding window, no chunks, no batches, no decode blocks. Each
  /// assigned precursor's cells are read straight off its own (window, cycle)
  /// spectra.
  ///
  /// The contract, as the header states it:
  ///   * a window's spectra, in acquisition order, are its cycles; their
  ///     FLOAT times are its axis;
  ///   * a precursor goes to the covering window with a non-empty axis whose
  ///     centre is nearest (strict, so a tie keeps the lower index);
  ///   * its cycle range is [first cycle at or after centre - half, first at or
  ///     after centre + half) on that float axis, clipped to the same lookup
  ///     of rt_low / rt_high when those are set;
  ///   * a peak counts toward a transition when |m - mz| <= mz * ppm, it lies
  ///     in the spectrum's mobility band [im_low, im_high), and it is within
  ///     precursor_im_window of the library 1/K0 when both are known; Sum
  ///     accumulates in peak order, the residual planes accumulate
  ///     intensity * deviation and intensity for positive intensities;
  ///   * the stream hands over the assigned precursors ordered by (the batch
  ///     of 128 spectra from the pass's first spectrum that holds their last
  ///     cycle, window, last cycle), then every unassigned precursor that has
  ///     transitions, empty, in library order.
  struct Oracle
  {
    std::vector<RecordingSink::Trace> assigned;   ///< by library index
    std::vector<std::uint32_t> tail;              ///< handed over empty, in order
    std::vector<std::tuple<std::size_t, std::uint32_t, std::uint32_t>> key;   ///< per assigned
    std::vector<std::uint32_t> valid;             ///< valid transitions per assigned
    std::size_t first_spectrum = 0;
    /// Nonzero cells on an edge cycle that the float axis admitted although
    /// its spectrum's double time lies outside [rt_low, rt_high) -- the cells
    /// a spectrum range searched in the double times never reads.
    std::size_t live_edge_cells = 0;
  };

  Oracle bruteForce(const ScriptedRun& run, const ODIA::Library& lib,
                    const ODIA::ChromatogramExtractor::Options& opt)
  {
    const auto& info = run.spectra();
    const auto& windows = run.windows();
    const auto& p = lib.precursors();
    const auto& t = lib.transitions();
    const std::size_t W = windows.size();

    std::vector<std::vector<std::size_t>> pos(W);
    std::vector<std::vector<float>> axis(W);
    for (std::size_t si = 0; si < info.size(); ++si)
    {
      for (std::size_t w = 0; w < W; ++w)
      {
        if (std::abs(info[si].window.mz_low - windows[w].mz_low) < 1e-6 &&
            std::abs(info[si].window.mz_high - windows[w].mz_high) < 1e-6)
        {
          pos[w].push_back(si);
          axis[w].push_back(static_cast<float>(info[si].retention_time));
          break;
        }
      }
    }
    const auto firstAtOrAfter = [&](std::size_t w, double when) {
      std::size_t i = 0;
      while (i < axis[w].size() && double(axis[w][i]) < when) { ++i; }
      return i;
    };
    const bool restricted = opt.rt_high > opt.rt_low;
    std::vector<std::size_t> glo(W, 0), ghi(W, 0);
    for (std::size_t w = 0; w < W; ++w)
    {
      ghi[w] = axis[w].size();
      if (restricted) { glo[w] = firstAtOrAfter(w, opt.rt_low); ghi[w] = firstAtOrAfter(w, opt.rt_high); }
    }

    Oracle o;
    // The pass starts at the first spectrum at or after rt_low in the run's
    // own times, or earlier if a window's float axis admitted an earlier one.
    if (restricted)
    {
      o.first_spectrum = info.size();
      for (std::size_t si = 0; si < info.size(); ++si)
      { if (info[si].retention_time >= opt.rt_low) { o.first_spectrum = si; break; } }
      for (std::size_t w = 0; w < W; ++w)
      { if (glo[w] < ghi[w]) { o.first_spectrum = std::min(o.first_spectrum, pos[w][glo[w]]); } }
    }

    const double ppm = opt.fragment_ppm * 1e-6;
    for (std::size_t i = 0; i < lib.precursorCount(); ++i)
    {
      const double mz = ODIA::fromFixed(p.mz[i]);
      std::size_t best = W;
      double best_offset = std::numeric_limits<double>::infinity();
      for (std::size_t w = 0; w < W; ++w)
      {
        if (!windows[w].contains(mz) || axis[w].empty()) { continue; }
        const double off = std::abs(mz - windows[w].centre());
        if (off < best_offset) { best_offset = off; best = w; }
      }
      std::size_t lo = 0, hi = 0;
      if (best < W)
      {
        lo = glo[best]; hi = ghi[best];
        if (opt.irt_slope != 0.0)
        {
          const double centre = opt.irt_slope * double(p.irt[i]) + opt.irt_intercept;
          if (std::isnan(centre)) { hi = lo; }
          else
          {
            lo = std::max(lo, firstAtOrAfter(best, centre - opt.rt_window_seconds));
            hi = std::min(hi, firstAtOrAfter(best, centre + opt.rt_window_seconds));
          }
        }
      }
      if (best == W || lo >= hi)
      {
        if (p.transition_count[i] != 0) { o.tail.push_back(static_cast<std::uint32_t>(i)); }
        continue;
      }

      RecordingSink::Trace tr;
      tr.precursor = static_cast<std::uint32_t>(i);
      tr.transition_begin = p.transition_begin[i];
      tr.transition_count = p.transition_count[i];
      tr.axis = static_cast<std::uint32_t>(best);
      tr.axis_begin = static_cast<std::uint32_t>(lo);
      tr.cycles = static_cast<std::uint32_t>(hi - lo);
      tr.rt.assign(axis[best].begin() + std::ptrdiff_t(lo), axis[best].begin() + std::ptrdiff_t(hi));
      std::uint32_t valid = 0;
      for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
      {
        const std::size_t j = p.transition_begin[i] + k;
        const bool ok = t.product_mz[j] != ODIA::MZ_INVALID;
        const std::size_t n = ok ? hi - lo : 0;
        const bool mass = opt.collect_mass_residuals && ok, im = opt.collect_im_residuals && ok;
        tr.points.emplace_back(n, 0.0f);
        tr.ppm_num.emplace_back(mass ? n : 0, 0.0f);
        tr.ppm_den.emplace_back(mass ? n : 0, 0.0f);
        tr.im_num.emplace_back(im ? n : 0, 0.0f);
        tr.im_den.emplace_back(im ? n : 0, 0.0f);
        if (!ok) { continue; }
        ++valid;
        const double target = ODIA::fromFixed(t.product_mz[j]);
        for (std::size_t c = lo; c < hi; ++c)
        {
          const std::size_t si = pos[best][c];
          const auto& pk = run.peakList(si);
          const bool has_im = pk.hasIonMobility();
          for (std::size_t q = 0; q < pk.mz.size(); ++q)
          {
            const double m = pk.mz[q];
            const double peak_im = has_im ? double(pk.ion_mobility[q])
                                          : std::numeric_limits<double>::quiet_NaN();
            if (opt.use_ion_mobility && has_im &&
                (peak_im < info[si].window.im_low || peak_im >= info[si].window.im_high))
            { continue; }
            if (std::abs(m - target) > target * ppm) { continue; }
            if (opt.precursor_im_window > 0.0 && !std::isnan(peak_im))
            {
              const float want = p.im[i];
              if (!std::isnan(want) && std::abs(peak_im - want) > opt.precursor_im_window) { continue; }
            }
            const float intensity = pk.intensity[q];
            const std::size_t at = c - lo;
            tr.points.back()[at] += intensity;
            if (mass && intensity > 0.0f)
            {
              tr.ppm_num.back()[at] += intensity * static_cast<float>((m - target) / target * 1e6);
              tr.ppm_den.back()[at] += intensity;
            }
            if (im && intensity > 0.0f && !std::isnan(peak_im))
            {
              tr.im_num.back()[at] += intensity * static_cast<float>(peak_im);
              tr.im_den.back()[at] += intensity;
            }
          }
        }
        // The cells a range searched in the run's DOUBLE times would lose: an
        // edge cycle the float axis admitted whose spectrum's own time lies
        // outside [rt_low, rt_high).
        if (restricted &&
            ((lo == glo[best] && info[pos[best][lo]].retention_time < opt.rt_low &&
              tr.points.back().front() != 0.0f) ||
             (hi == ghi[best] && info[pos[best][hi - 1]].retention_time >= opt.rt_high &&
              tr.points.back().back() != 0.0f)))
        { ++o.live_edge_cells; }
      }
      o.key.emplace_back((pos[best][hi - 1] - o.first_spectrum) / 128, std::uint32_t(best),
                         std::uint32_t(hi));
      o.valid.push_back(valid);
      o.assigned.push_back(std::move(tr));
    }
    return o;
  }

  /// How a stream departs from the oracle, or empty when it does not.
  std::string againstOracle(const Oracle& o, const std::vector<RecordingSink::Trace>& got)
  {
    std::size_t wrong = 0, unknown = 0, duplicate = 0, out_of_order = 0;
    std::vector<int> seen(o.assigned.size(), 0);
    std::size_t at = 0;
    std::tuple<std::size_t, std::uint32_t, std::uint32_t> previous{0, 0, 0};
    // The assigned precursors come first ...
    for (; at < got.size() && got[at].cycles != 0; ++at)
    {
      const auto& g = got[at];
      std::size_t idx = o.assigned.size();
      for (std::size_t e = 0; e < o.assigned.size(); ++e)
      { if (o.assigned[e].precursor == g.precursor) { idx = e; break; } }
      if (idx == o.assigned.size()) { ++unknown; continue; }
      if (seen[idx]++) { ++duplicate; }
      if (!sameBits(g, o.assigned[idx])) { ++wrong; }
      if (o.key[idx] < previous) { ++out_of_order; }
      previous = o.key[idx];
    }
    const std::size_t missing =
      std::size_t(std::count(seen.begin(), seen.end(), 0));
    // ... then the empty ones, in library order.
    std::vector<std::uint32_t> tail;
    bool tail_empty = true;
    for (; at < got.size(); ++at)
    {
      tail.push_back(got[at].precursor);
      tail_empty = tail_empty && got[at].cycles == 0;
    }
    if (!wrong && !unknown && !duplicate && !out_of_order && !missing &&
        tail == o.tail && tail_empty)
    { return {}; }
    return std::to_string(wrong) + " traces differ from the brute force, " +
           std::to_string(missing) + " missing, " + std::to_string(unknown) + " unexpected, " +
           std::to_string(duplicate) + " duplicated, " + std::to_string(out_of_order) +
           " out of order" + (tail == o.tail && tail_empty ? "" : ", wrong empty tail");
  }

  /// The run every matrix arm reads, built so each boundary the extractor has
  /// is crossed somewhere:
  ///
  ///   * times 1000.37 + 0.1 c, which round both ways in float;
  ///   * two co-packed windows sharing every frame time and one peak list,
  ///     split by mobility band, so a chunk can start or end BETWEEN the two
  ///     spectra of one frame; and a third window acquired between frames;
  ///   * a 30 s gap inside the run (a batch straddling it holds both sides
  ///     live), empty spectra, and missing acquisitions on the third window
  ///     (its cycle index then runs ahead of the frame index);
  ///   * per transition, an on-target peak and a jittered one that is
  ///     sometimes outside the 10 ppm tolerance, the mobility band or the
  ///     precursor's 1/K0 window, plus occasional zero-intensity peaks, so
  ///     all five planes carry sums whose float order matters;
  ///   * precursors spread over the run with windows clipped at its ends and
  ///     at the gap, interleaved with the awkward kinds: no covering window,
  ///     predicted outside the run, NaN iRT, no transitions at all (inside and
  ///     outside a window), and transitions none of which has a product m/z.
  struct MatrixRun
  {
    LoggingRun run;
    ScriptedLibrary lib;
    static constexpr double T0 = 1000.37, STEP = 0.1, GAP = 30.0, HALF = 6.0;
    static constexpr std::uint32_t CYCLES = 700, GAP_AT = 350;

    MatrixRun()
    {
      const auto w0 = run.addWindow(500.0, 510.0, 0.80, 1.00);
      const auto w1 = run.addWindow(510.001, 520.0, 1.00, 1.20);
      const auto w2 = run.addWindow(520.001, 530.0, 0.70, 1.30);
      const double window_mz[3] = {505.0, 515.0, 525.0};
      const float window_im[3] = {0.90f, 1.10f, 1.00f};

      struct Product { double mz; float im; std::uint32_t window; std::uint32_t id; };
      std::vector<Product> products;
      const double end = T0 + STEP * CYCLES + GAP;
      constexpr std::uint32_t REGULAR = 60;
      std::uint32_t regular = 0;
      for (std::uint32_t slot = 0; regular < REGULAR; ++slot)
      {
        // The awkward kinds, interleaved so the empty tail's ORDER is tested.
        switch (slot)
        {
          case 7:  lib.addPrecursor(600.0); lib.addTransition(300.0); lib.addTransition(400.0);
                   continue;                                  // no window covers it
          case 15: { const auto i = lib.addPrecursor(505.0); lib.library().precursors().irt[i] =
                     float(T0 - 5000.0); lib.addTransition(310.0); continue; }   // before the run
          case 22: { const auto i = lib.addPrecursor(515.0); lib.library().precursors().irt[i] =
                     std::numeric_limits<float>::quiet_NaN(); lib.addTransition(320.0); continue; }
          case 30: { const auto i = lib.addPrecursor(505.0); lib.library().precursors().irt[i] =
                     float(T0 + 40.0); continue; }              // in a window, no transitions
          case 37: lib.addPrecursor(700.0); continue;           // no window, no transitions
          case 44: { const auto i = lib.addPrecursor(515.0); lib.library().precursors().irt[i] =
                     float(T0 + 60.0); lib.addTransition(0.0); lib.addTransition(0.0);
                     continue; }                                // no valid transition
          default: break;
        }
        const std::uint32_t w = regular % 3;
        const float im = (regular % 10 == 3) ? NA : window_im[w];
        const auto i = lib.addPrecursor(window_mz[w], im);
        lib.library().precursors().irt[i] =
          static_cast<float>(T0 - 8.0 + (end - T0 + 16.0) * double(regular) / double(REGULAR - 1));
        for (std::uint32_t k = 0; k < 3; ++k)
        {
          const double mz = 200.0 + 0.37 * regular + 310.0 * k;
          lib.addTransition(mz);
          products.push_back({mz, window_im[w], w, regular * 3 + k});
        }
        ++regular;
      }

      std::uint64_t state = 0x9E3779B97F4A7C15ull;
      const auto uniform = [&]() {
        state = state * 6364136223846793005ull + 1442695040888963407ull;
        return double(state >> 11) / double(1ull << 53);
      };
      struct Peak { double mz; float intensity, im; };
      const auto frame = [&](std::uint32_t c, std::uint32_t window_a, std::uint32_t window_b) {
        std::vector<Peak> list;
        for (const auto& pr : products)
        {
          if (pr.window != window_a && pr.window != window_b) { continue; }
          const float base = 1.0f + float((c * 7 + pr.id * 3) % 50) + 0.25f;
          // On target, inside its band and its precursor's 1/K0 window.
          list.push_back({pr.mz * (1.0 + (uniform() * 16.0 - 8.0) * 1e-6), base,
                          pr.im + float(uniform() * 0.03 - 0.015)});
          // Jittered: sometimes outside the tolerance, the band, or the window.
          list.push_back({pr.mz * (1.0 + (uniform() * 24.0 - 12.0) * 1e-6),
                          0.5f + float(uniform()) * 3.0f, float(0.75 + uniform() * 0.5)});
          if ((c + pr.id) % 17 == 0) { list.push_back({pr.mz, 0.0f, pr.im}); }
        }
        return list;
      };
      const auto fill = [&](std::size_t s, const std::vector<Peak>& list) {
        for (const auto& pk : list) { run.addPeak(s, pk.mz, pk.intensity, pk.im); }
      };
      for (std::uint32_t c = 0; c < CYCLES; ++c)
      {
        const double rt = T0 + STEP * c + (c >= GAP_AT ? GAP : 0.0);
        // One frame, two windows, ONE peak list -- as the reader serves
        // diaPASEF -- empty on every eleventh frame.
        const auto shared = c % 11 == 4 ? std::vector<Peak>{} : frame(c, 0, 1);
        fill(run.addSpectrum(w0, rt), shared);
        fill(run.addSpectrum(w1, rt), shared);
        if (c % 37 == 5) { continue; }                          // a missing acquisition
        const auto s2 = run.addSpectrum(w2, rt + STEP / 2.0);
        if (c % 13 != 7) { fill(s2, frame(c, 2, 2)); }
      }
    }
  };

  ODIA::ChromatogramExtractor::Options matrixOptions()
  {
    auto opt = plainOptions();
    opt.irt_slope = 1.0;
    opt.irt_intercept = 0.0;
    opt.rt_window_seconds = MatrixRun::HALF;
    opt.collect_mass_residuals = true;      // the engine's default planes
    opt.collect_im_residuals = true;
    return opt;
  }

  /// Every decode block x chunk count x thread count, each against the
  /// independent oracle AND against the default arm, bitwise in all five
  /// planes and in order; then restricted RT bounds on both sides of float
  /// rounding, and the explicit and byte caps together.
  void caseChunkMatrix()
  {
    MatrixRun m;
    auto& run = m.run;
    const auto& lib = m.lib.library();
    const auto base = matrixOptions();
    const Oracle oracle = bruteForce(run, lib, base);
    std::printf("       %zu spectra, %zu precursors: %zu assigned, %zu handed over empty\n",
                run.spectra().size(), lib.precursorCount(), oracle.assigned.size(),
                oracle.tail.size());

    struct Arm
    {
      std::vector<RecordingSink::Trace> traces;
      ODIA::ChromatogramExtractor::Stats st;
      std::vector<std::pair<std::size_t, std::size_t>> calls;
      std::string error;
    };
    const auto extract = [&](const ODIA::ChromatogramExtractor::Options& o) {
      Arm a;
      RecordingSink sink;
      run.calls.clear();
      try { ODIA::ChromatogramExtractor::extract(lib, run, o, sink, &a.st); }
      catch (const std::exception& e) { a.error = e.what(); }
      a.traces = std::move(sink.traces);
      a.calls = run.calls;
      return a;
    };

    const Arm reference = extract(base);
    check(reference.error.empty() && reference.st.chunks == 1, "the default arm is one chunk");
    const std::string ref_vs_oracle = againstOracle(oracle, reference.traces);
    check(ref_vs_oracle.empty(), "the default arm equals the brute force" +
                                 (ref_vs_oracle.empty() ? "" : ": " + ref_vs_oracle));

    // A chunk boundary falls between the two spectra of one frame when the
    // spectrum before it has the same time.
    const auto& info = run.spectra();
    const auto between_frame = [&](std::size_t b) {
      return b > 0 && b < info.size() && info[b - 1].retention_time == info[b].retention_time;
    };

    std::size_t arms = 0, failed = 0, split_frames = 0, printed = 0;
    std::set<std::size_t> plans;
    for (const std::size_t block : {1, 64, 127, 128, 129, 200, 256})
    {
      for (const std::size_t cap : {0, 2, 4, 7})
      {
        for (const unsigned threads : {1u, 4u})
        {
          auto o = base;
          o.decode_block = block;
          o.max_live_precursors = cap;
          o.threads = threads;
          const Arm a = extract(o);
          ++arms;
          if (a.st.chunks > 1) { plans.insert(a.st.chunks); }
          std::string why = a.error;
          if (why.empty())
          {
            why = againstOracle(oracle, a.traces);
            if (why.empty() && !sameStream(a.traces, reference.traces))
            { why = "not the default arm's stream"; }
            if (why.empty() && cap != 0 && a.st.peak_live_precursors > cap)
            { why = "peak live " + std::to_string(a.st.peak_live_precursors) + " over the cap"; }
          }
          if (!why.empty())
          {
            ++failed;
            if (printed++ < 12)
            {
              std::printf("       decode_block %zu, cap %zu, %u threads (%zu chunks): %s\n",
                          block, cap, threads, a.st.chunks, why.c_str());
            }
          }
          // Where the chunks began and ended: a call that does not continue
          // the previous one starts a chunk.
          for (std::size_t c = 0; c < a.calls.size(); ++c)
          {
            const bool starts = c == 0 || a.calls[c].first != a.calls[c - 1].second;
            const bool ends = c + 1 == a.calls.size() || a.calls[c + 1].first != a.calls[c].second;
            if (a.st.chunks > 1 && ((starts && between_frame(a.calls[c].first)) ||
                                    (ends && between_frame(a.calls[c].second))))
            { ++split_frames; }
          }
        }
      }
    }
    std::printf("       %zu arms, %zu failed; chunk plans of", arms, failed);
    for (const auto n : plans) { std::printf(" %zu", n); }
    std::printf(" chunks; %zu chunk edges between the two windows of one frame\n", split_frames);
    check(failed == 0, "every decode block x cap x thread count equals the brute force and "
                       "the default arm, bit for bit and in order");
    check(plans.size() >= 3, "the caps give at least three different chunk plans");
    check(split_frames > 0, "some chunk begins or ends between the two windows of one frame");

    // Restricted RT bounds that fall between a spectrum's double time and its
    // float rounding, on both sides. In the first pair each bound admits, on
    // the float axis, a spectrum the double times exclude -- a range searched
    // in the double times never reads it and its cells stay zero. In the
    // second each bound excludes one the double times admit.
    std::size_t up_lo = 0, down_lo = 0, up_hi = 0, down_hi = 0;
    for (std::size_t si = 0; si < info.size(); ++si)
    {
      const double rt = info[si].retention_time, f = double(static_cast<float>(rt));
      const bool has_peaks = !run.peakList(si).mz.empty();
      const bool early = has_peaks && rt > MatrixRun::T0 + 10.0 && rt < MatrixRun::T0 + 20.0;
      const bool late = has_peaks && rt > MatrixRun::T0 + 80.0 && rt < MatrixRun::T0 + 90.0;
      if (early && f > rt && !up_lo) { up_lo = si; }
      if (early && f < rt && !down_lo) { down_lo = si; }
      if (late && f > rt && !up_hi) { up_hi = si; }
      if (late && f < rt && !down_hi) { down_hi = si; }
    }
    const auto between = [&](std::size_t si) {
      const double rt = info[si].retention_time;
      return 0.5 * (rt + double(static_cast<float>(rt)));
    };
    std::size_t restricted_failed = 0, edge_cells = 0;
    for (int pair = 0; pair < 2; ++pair)
    {
      auto o = base;
      o.rt_low = between(pair == 0 ? up_lo : down_lo);
      o.rt_high = between(pair == 0 ? down_hi : up_hi);
      const Oracle ro = bruteForce(run, lib, o);
      edge_cells += pair == 0 ? ro.live_edge_cells : 0;
      const Arm rref = extract(o);
      for (const std::size_t block : {256, 129})
      {
        for (const std::size_t cap : {0, 3})
        {
          auto oc = o;
          oc.decode_block = block;
          oc.max_live_precursors = cap;
          const Arm a = extract(oc);
          std::string why = a.error;
          if (why.empty()) { why = againstOracle(ro, a.traces); }
          if (why.empty() && !sameStream(a.traces, rref.traces)) { why = "not the 1-chunk stream"; }
          if (why.empty() && cap != 0 && a.st.peak_live_precursors > cap) { why = "over the cap"; }
          if (!why.empty())
          {
            ++restricted_failed;
            std::printf("       rt [%.9f, %.9f) %s, block %zu, cap %zu (%zu chunks): %s\n",
                        o.rt_low, o.rt_high, pair == 0 ? "float admits more" : "float admits less",
                        block, cap, a.st.chunks, why.c_str());
          }
        }
      }
    }
    std::printf("       restricted bounds: %zu precursors hold a nonzero cell on a float-"
                "admitted edge cycle\n", edge_cells);
    check(edge_cells > 0, "the float-admitted edge spectra carry signal, so losing one shows");
    check(restricted_failed == 0,
          "restricted RT bounds on both sides of float rounding equal the brute force, "
          "chunked or not");

    // Both caps at once. The byte budget is inverted with the mean block over
    // the assigned precursors, all five planes; the tighter cap must bind.
    double mean_cells = 0.0;
    for (std::size_t e = 0; e < oracle.assigned.size(); ++e)
    { mean_cells += double(oracle.valid[e]) * double(oracle.assigned[e].cycles); }
    mean_cells /= double(oracle.assigned.size());
    const double per = mean_cells * 4.0 * 5.0;
    const std::size_t budget = std::size_t(per * 3.5);
    const std::size_t derived = std::size_t(double(budget) / per);
    std::size_t both_failed = 0;
    for (const std::size_t explicit_cap : {std::size_t(2), std::size_t(5), std::size_t(0)})
    {
      auto o = base;
      o.max_live_precursors = explicit_cap;
      o.live_memory_budget_bytes = budget;
      o.decode_block = 200;
      const std::size_t want = explicit_cap == 0 ? derived : std::min(explicit_cap, derived);
      const Arm a = extract(o);
      const std::string note_tail = "cap " + std::to_string(want);
      const bool ok = a.error.empty() && againstOracle(oracle, a.traces).empty() &&
                      sameStream(a.traces, reference.traces) &&
                      a.st.peak_live_precursors <= want && a.st.chunks > 1 &&
                      a.st.live_budget_note.size() >= note_tail.size() &&
                      a.st.live_budget_note.compare(a.st.live_budget_note.size() - note_tail.size(),
                                                    note_tail.size(), note_tail) == 0;
      if (!ok)
      {
        ++both_failed;
        std::printf("       explicit %zu + budget (derived %zu): %zu chunks, peak live %zu, "
                    "note \"%s\"%s\n", explicit_cap, derived, a.st.chunks,
                    a.st.peak_live_precursors, a.st.live_budget_note.c_str(),
                    a.error.empty() ? "" : (" -- " + a.error).c_str());
      }
    }
    check(derived == 3, "the byte budget inverts to a cap of 3");
    check(both_failed == 0, "with both caps the tighter binds, and the stream is unchanged");
  }

  /// The instrumentation on the chunked path: the per-chunk stage ledger and
  /// the plane-identity counter must change nothing the sink sees, at any
  /// decode block x cap x thread count, and the ledger must account for the
  /// pass's actual structure -- one index bracket and one snapshot per chunk,
  /// one decode bracket per decoded block, one activate bracket per batch, a
  /// match bracket per part of a batch (a batch straddling two decode blocks is
  /// matched in two parts), one emit bracket per batch plus the chunk's flush,
  /// and no bracket overlapping another. The plane-identity counts are a
  /// property of the cells, so they must not move with the chunk plan either.
  void caseInstrumentMatrix()
  {
    MatrixRun m;
    auto& run = m.run;
    const auto& lib = m.lib.library();
    const auto base = matrixOptions();
    const Oracle oracle = bruteForce(run, lib, base);

    const auto extract = [&](const ODIA::ChromatogramExtractor::Options& o,
                             ODIA::ChromatogramExtractor::Stats& st, std::size_t& decodes,
                             std::string& error) {
      RecordingSink sink;
      run.calls.clear();
      try { ODIA::ChromatogramExtractor::extract(lib, run, o, sink, &st); }
      catch (const std::exception& e) { error = e.what(); }
      decodes = run.calls.size();
      return std::move(sink.traces);
    };

    ODIA::ChromatogramExtractor::Stats ref_st;
    std::size_t ref_decodes = 0;
    std::string ref_error;
    const auto reference = extract(base, ref_st, ref_decodes, ref_error);
    check(ref_error.empty() && againstOracle(oracle, reference).empty(),
          "the uninstrumented default arm equals the brute force");

    std::size_t arms = 0, failed = 0, printed = 0, multi = 0;
    bool counted = false;
    std::uint64_t cells = 0, ppm = 0, im = 0, im_zero = 0, nonfinite = 0, negative = 0;
    for (const std::size_t block : {1, 128, 129, 200, 256})
    {
      for (const std::size_t cap : {0, 2, 7})
      {
        for (const unsigned threads : {1u, 4u})
        {
          auto o = base;
          o.decode_block = block;
          o.max_live_precursors = cap;
          o.threads = threads;
          o.stage_label = "instr_matrix";
          o.count_plane_identity = true;
          ODIA::ChromatogramExtractor::Stats st;
          std::size_t decodes = 0;
          std::string why;
          const auto traces = extract(o, st, decodes, why);
          ++arms;
          if (st.chunks > 1) { ++multi; }
          const auto bad = [&](bool ok, const std::string& what) {
            if (why.empty() && !ok) { why = what; }
          };
          if (why.empty()) { why = againstOracle(oracle, traces); }
          bad(sameStream(traces, reference), "not the uninstrumented default stream");
          bad(st.index.brackets == st.chunks, "index brackets " +
              std::to_string(st.index.brackets) + " != chunks " + std::to_string(st.chunks));
          bad(st.probe.brackets == st.chunks, "probe brackets != chunks");
          bad(st.decode.brackets == decodes, "decode brackets " +
              std::to_string(st.decode.brackets) + " != decode calls " + std::to_string(decodes));
          bad(st.emit.brackets == st.activate.brackets + st.chunks,
              "emit brackets != batches + chunk flushes");
          bad(st.match.brackets >= st.activate.brackets &&
              st.match.brackets <= st.activate.brackets + st.decode.brackets,
              "match brackets outside [batches, batches + decodes]");
          bad(st.bracket_overlaps == 0, "a bracket overlaps another");
          bad(st.setup.brackets == 1 && st.empty.brackets == 1, "setup/empty not one bracket each");
          bad(st.plane_cells > 0, "the plane counter saw no cells");
          if (!counted)
          {
            counted = true;
            cells = st.plane_cells; ppm = st.ppm_den_differs; im = st.im_den_differs;
            im_zero = st.im_den_zero; nonfinite = st.base_nonfinite; negative = st.base_negative;
          }
          bad(st.plane_cells == cells && st.ppm_den_differs == ppm && st.im_den_differs == im &&
              st.im_den_zero == im_zero && st.base_nonfinite == nonfinite &&
              st.base_negative == negative, "plane-identity counts move with the chunk plan");
          if (!why.empty())
          {
            ++failed;
            if (printed++ < 12)
            {
              std::printf("       decode_block %zu, cap %zu, %u threads (%zu chunks): %s\n",
                          block, cap, threads, st.chunks, why.c_str());
            }
          }
        }
      }
    }
    std::printf("       %zu instrumented arms (%zu chunked), %zu failed; plane cells %llu, "
                "ppm_den differs %llu, im_den differs %llu\n", arms, multi, failed,
                static_cast<unsigned long long>(cells), static_cast<unsigned long long>(ppm),
                static_cast<unsigned long long>(im));
    check(multi > 0, "some instrumented arm is chunked");
    check(failed == 0, "the stage ledger and the plane counter change nothing extracted, "
                       "account for every chunk, block, batch and part, and do not move "
                       "with the chunk plan");
  }

  /// Nothing assigned: no precursor has a covering window. Every precursor
  /// with transitions is still handed over, empty, in library order -- chunked
  /// or not, budgeted or not -- and nothing throws.
  void caseNoAssignments()
  {
    ScriptedRun run;
    const auto w = run.addWindow(500.0, 510.0);
    for (int c = 0; c < 300; ++c)
    {
      const auto s = run.addSpectrum(w, 100.37 + 0.1 * c);
      if (c % 3) { run.addPeak(s, 400.0, 5.0f); }
    }
    ScriptedLibrary lib;
    lib.addPrecursor(600.0); lib.addTransition(400.0);
    lib.addPrecursor(700.0);                                   // no transitions: omitted
    lib.addPrecursor(450.0); lib.addTransition(400.0); lib.addTransition(0.0);

    std::size_t failed = 0;
    for (const std::size_t cap : {0, 1})
    {
      for (const std::size_t budget : {std::size_t(0), std::size_t(1) << 20})
      {
        auto opt = plainOptions();
        opt.max_live_precursors = cap;
        opt.live_memory_budget_bytes = budget;
        RecordingSink sink;
        ODIA::ChromatogramExtractor::Stats st;
        std::string error;
        try { ODIA::ChromatogramExtractor::extract(lib.library(), run, opt, sink, &st); }
        catch (const std::exception& e) { error = e.what(); }
        const bool ok = error.empty() && st.precursors_extracted == 0 &&
                        st.precursors_without_window == 3 && sink.traces.size() == 2 &&
                        sink.traces[0].precursor == 0 && sink.traces[1].precursor == 2 &&
                        sink.traces[0].cycles == 0 && sink.traces[1].cycles == 0;
        if (!ok)
        {
          ++failed;
          std::printf("       cap %zu budget %zu: %zu traces, %zu extracted %s\n", cap, budget,
                      sink.traces.size(), st.precursors_extracted, error.c_str());
        }
      }
    }
    check(failed == 0, "with nothing assigned every precursor with transitions is handed "
                       "over empty, in library order, under any cap");
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
  else if (which == "chunk_invariant") { caseChunkInvariant(); }
  else if (which == "cap_zero_rows") { caseCapZeroRows(); }
  else if (which == "chunk_matrix") { caseChunkMatrix(); }
  else if (which == "no_assignments") { caseNoAssignments(); }
  else if (which == "plane_identity") { casePlaneIdentity(); }
  else if (which == "instr_matrix") { caseInstrumentMatrix(); }
  else
  {
    std::fprintf(stderr,
                 "usage: odia_extract_cases "
                 "<invalid_mz|aggregate|mobility|im_gating|band_edge|wide_csr|sliding|"
                 "chunk_invariant|cap_zero_rows|chunk_matrix|no_assignments|"
                 "plane_identity|instr_matrix>\n");
    return 2;
  }

  std::printf(failures ? "FAILED (%d)\n" : "ok (%d)\n", failures);
  return failures ? 1 : 0;
}
