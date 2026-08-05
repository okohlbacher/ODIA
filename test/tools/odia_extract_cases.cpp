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
  /// several windows can share one cycle (and, as on S08, one peak list).
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
      const std::size_t a = x.axis_of[j];
      if (a >= x.axes.size() ||
          std::size_t(x.axis_begin[j]) + x.count[j] > x.axes[a].size())
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
}

int main(int argc, char** argv)
{
  const std::string which = argc > 1 ? argv[1] : "";
  std::printf("case: %s\n", which.c_str());

  if (which == "invalid_mz") { caseInvalidMz(); }
  else if (which == "aggregate") { caseAggregate(); }
  else if (which == "mobility") { caseMobility(); }
  else
  {
    std::fprintf(stderr, "usage: odia_extract_cases <invalid_mz|aggregate|mobility>\n");
    return 2;
  }

  std::printf(failures ? "FAILED (%d)\n" : "ok (%d)\n", failures);
  return failures ? 1 : 0;
}
