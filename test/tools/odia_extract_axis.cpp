// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// The retention-time axis of a precursor that falls in MORE THAN ONE isolation
// window -- 1.1% of them on 12_80 and on IH1, whose schemes overlap adjacent
// windows by 1.0 Th.
//
// The representation (2be7d84) stores one (axis_of, axis_begin) pair per
// transition, which says a transition's points are one contiguous run on one
// window's axis. The extractor used to CONCATENATE the windows' cycles into
// that one run, so the pair could only name one of them: the last window won
// both fields, every point of the first window was read from the second's axis
// starting at the second's first live cycle, and `axis_begin + j` ran off the
// end of that axis -- a heap over-read through retentionTime(), confirmed under
// ASan.
//
// The extractor now extracts each precursor from ONE window, the one whose
// centre it is nearest, and counts the rest in `precursors_in_several_windows`.
// This asserts that choice, the values that follow from it, and the invariant
// the representation rests on.
//
// No real run is needed to show any of it: two overlapping windows, three
// precursors, four cycles each.

#include <odia/ChromatogramExtractor.h>
#include <odia/SpectrumSource.h>

#include <cmath>
#include <cstdio>
#include <vector>

namespace
{
  /// Two overlapping windows, interleaved one spectrum each per cycle.
  class TwoWindowRun : public ODIA::SpectrumSource
  {
  public:
    TwoWindowRun()
    {
      windows_.push_back({500.0, 510.0});   // centre 505
      windows_.push_back({509.0, 519.0});   // centre 514
      for (int c = 0; c < CYCLES; ++c)
      {
        for (int w = 0; w < 2; ++w)
        {
          ODIA::SpectrumInfo s;
          s.index = info_.size();
          // Distinct per (window, cycle) so a wrong axis lookup is visible.
          s.retention_time = 100.0 + 10.0 * c + 5.0 * w;
          s.window = windows_[w];
          info_.push_back(s);
        }
      }
    }

    static constexpr int CYCLES = 4;

    const std::vector<ODIA::SpectrumInfo>& spectra() const override { return info_; }
    const std::vector<ODIA::IsolationWindow>& windows() const override { return windows_; }
    std::string describe() const override { return "two overlapping windows"; }

    void peaks(std::size_t begin, std::size_t end,
               std::vector<ODIA::SpectrumPeaks>& out) override
    {
      out.assign(end - begin, {});
      for (std::size_t i = begin; i < end; ++i)
      {
        out[i - begin].mz.push_back(400.0);
        // Intensity names the SPECTRUM, so the trace says which window's
        // spectra were read and not merely how many points came back.
        out[i - begin].intensity.push_back(1.0f + static_cast<float>(i));
      }
    }

  private:
    std::vector<ODIA::SpectrumInfo> info_;
    std::vector<ODIA::IsolationWindow> windows_;
  };

  /// Three precursors: two inside BOTH windows and one inside only window 0.
  ///
  /// 509.2 is nearer window 0's centre and 509.8 nearer window 1's, so the two
  /// overlap cases resolve to different windows and a fix that simply always
  /// took the first (or the last) window would fail one of them. 505.0 is the
  /// ordinary single-window case, and it is what says the overlap count counts
  /// overlaps rather than precursors.
  ODIA::Library sharedPrecursors()
  {
    ODIA::Library lib;
    auto& p = lib.precursors();
    auto& t = lib.transitions();
    const double mz[3] = {509.2, 509.8, 505.0};
    for (int i = 0; i < 3; ++i)
    {
      p.mz.push_back(ODIA::toFixed(mz[i]));
      p.irt.push_back(0.0f);
      p.im.push_back(std::nanf(""));       // disables the mobility test
      p.ccs.push_back(std::nanf(""));
      p.charge.push_back(2);
      p.decoy.push_back(0);
      p.modified_sequence.push_back(0);
      p.protein_group.push_back(0);
      p.transition_begin.push_back(static_cast<std::uint32_t>(i));
      p.transition_count.push_back(1);
      t.product_mz.push_back(ODIA::toFixed(400.0));
      t.library_intensity.push_back(1.0f);
      t.type.push_back(ODIA::FragmentType::Y);
      t.ordinal.push_back(4);
      t.charge.push_back(1);
      t.loss.push_back(ODIA::LossType::None);
    }
    return lib;
  }

  int failures = 0;

  void check(bool ok, const char* what)
  {
    if (!ok) { std::printf("FAIL: %s\n", what); ++failures; }
  }
}

int main()
{
  TwoWindowRun run;
  ODIA::Library lib = sharedPrecursors();

  ODIA::ChromatogramExtractor::Options opt;
  opt.threads = 1;
  opt.progress_every = 0;
  opt.precursor_im_window = 0.0;

  const auto x = ODIA::ChromatogramExtractor::extract(lib, run, opt);

  // Two of the three sit in both windows, and exactly one window is extracted.
  check(x.precursors_in_several_windows == 2,
        "the two precursors covered by both windows must be counted");
  check(x.precursors_without_window == 0, "every precursor is covered by a window");

  // Which window each transition was extracted from.
  const std::uint32_t want_window[3] = {0, 1, 0};

  for (std::uint32_t tr = 0; tr < 3; ++tr)
  {
    const std::uint32_t n = x.count[tr];
    const std::size_t axis_size =
      x.axisOf(tr) < x.axes.size() ? x.axes[x.axisOf(tr)].size() : 0;
    std::printf("transition %u: %u points, axis_of=%u axis_begin=%u, axes[%u].size()=%zu\n",
                tr, n, x.axisOf(tr), x.axisBegin(tr), x.axisOf(tr), axis_size);

    // One window's cycles, not two windows' concatenated.
    if (n != TwoWindowRun::CYCLES)
    {
      std::printf("FAIL: transition %u: expected %d points from one window, got %u\n",
                  tr, TwoWindowRun::CYCLES, n);
      ++failures;
      continue;
    }

    // The representation's own invariant: retentionTime(tr, j) indexes
    // axes[axis_of][axis_begin + j] with no bound of its own, so the whole run
    // must fit on the named axis. If it does not, every read past the end is
    // out of bounds. (The extractor also checks this itself and throws; this
    // asserts the result rather than trusting that check.)
    if (std::size_t(x.axisBegin(tr)) + n > axis_size)
    {
      std::printf("FAIL: transition %u: axis_begin(%u) + count(%u) = %zu exceeds axis size %zu "
                  "-- retentionTime() reads out of bounds\n",
                  tr, x.axisBegin(tr), n, std::size_t(x.axisBegin(tr)) + n, axis_size);
      ++failures;
      continue;
    }

    const std::uint32_t w = want_window[tr];
    if (x.axisOf(tr) != w)
    {
      std::printf("FAIL: transition %u: expected window %u (nearest centre), got %u\n",
                  tr, w, x.axisOf(tr));
      ++failures;
    }

    for (int c = 0; c < TwoWindowRun::CYCLES; ++c)
    {
      // Window w's cycle c is spectrum 2c + w, at 100 + 10c + 5w seconds and
      // carrying intensity 1 + (2c + w).
      const double want_rt = 100.0 + 10.0 * c + 5.0 * double(w);
      const double got_rt = x.retentionTime(tr, static_cast<std::uint32_t>(c));
      if (std::abs(got_rt - want_rt) > 1e-3)
      {
        std::printf("FAIL: transition %u point %d: want rt %.1f, got %.1f\n",
                    tr, c, want_rt, got_rt);
        ++failures;
      }
      const double want_intensity = 1.0 + double(2 * c) + double(w);
      const double got_intensity = x.intensity[x.begin[tr] + c];
      if (std::abs(got_intensity - want_intensity) > 1e-3)
      {
        std::printf("FAIL: transition %u point %d: want intensity %.1f from window %u, got %.1f\n",
                    tr, c, want_intensity, w, got_intensity);
        ++failures;
      }
    }

    // Times must ascend. A concatenated run does not, which is the reason the
    // scorer cannot be handed one.
    for (std::uint32_t j = 1; j < n; ++j)
    {
      check(x.retentionTime(tr, j) > x.retentionTime(tr, j - 1),
            "retention times must ascend along a transition's points");
    }
  }

  std::printf(failures ? "FAILED (%d)\n" : "ok (%d)\n", failures);
  return failures ? 1 : 0;
}
