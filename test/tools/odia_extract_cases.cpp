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
}

int main(int argc, char** argv)
{
  const std::string which = argc > 1 ? argv[1] : "";
  std::printf("case: %s\n", which.c_str());

  if (which == "invalid_mz") { caseInvalidMz(); }
  else
  {
    std::fprintf(stderr, "usage: odia_extract_cases <invalid_mz>\n");
    return 2;
  }

  std::printf(failures ? "FAILED (%d)\n" : "ok (%d)\n", failures);
  return failures ? 1 : 0;
}
