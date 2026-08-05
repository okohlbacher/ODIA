// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/SpectrumSource.h>

#include <limits>
#include <utility>

#include <mzpeak.h>

#include <algorithm>
#include <map>
#include <numeric>
#include <stdexcept>

namespace ODIA
{

  namespace
  {
    /// Quantise a window bound so that the same window acquired thousands of
    /// times collapses to one entry. Floating-point bounds recorded per
    /// spectrum differ in the last bits, and without this a 40-window scheme
    /// reports tens of thousands of distinct windows.
    long long quantise(double mz) { return static_cast<long long>(mz * 1000.0 + 0.5); }
  }

  /// mzPeak-backed run.
  ///
  /// Peak decode is currently ~284 ms/spectrum, ~2,450x slower than the file
  /// layout allows -- see doc/05-mzpeak-batched-reader-handoff.md. That is a
  /// property of the reader, not of ODIA or of the format, and it is being
  /// fixed elsewhere. This implementation is written against the batch entry
  /// point that the fix will make fast, so nothing here changes when it lands.
  ///
  /// Every extraction timing taken through this class before that happens
  /// measures the reader.
  class MzPeakSource : public SpectrumSource
  {
  public:
    explicit MzPeakSource(const std::string& filename)
      : filename_(filename), index_(MzPeak::open(filename.c_str())),
        spectra_(index_.spectra())
    {
      // Metadata only. mzPeak serves this faster than mzML and without
      // decoding a peak, which is what lets the whole extraction be planned
      // before anything expensive happens.
      std::map<long long, IsolationWindow> distinct;
      info_.reserve(spectra_.size());

      for (std::size_t i = 0; i < spectra_.size(); ++i)
      {
        const auto s = spectra_[i];
        if (s.ms_level() != 2) { continue; }

        SpectrumInfo info;
        info.index = i;
        info.ms_level = s.ms_level();
        // Seconds. The file stores minutes and the reader converts; passing
        // minutes on would be a 60x error that silently selects wrong scans.
        info.retention_time = s.retention_time().value_or(0.0);

        // diaPASEF window attribution.
        //
        // S08 packs TWO isolation windows into one frame, separated only in ion
        // mobility, and the reader hands both the same merged peak list. The
        // per-window mobility band would separate them -- but measured over all
        // 32,210 entries of that file, `ion_mobility_lower_limit` and
        // `_upper_limit` are NULL in every one, so there is no band to read.
        //
        // The assignment IS present, mis-attached: both SelectedIonInfo entries
        // land on precursor #0 and precursor #1 gets none. Each carries a
        // `selected_ion_mz` equal to one of the two window centres, and an
        // `ion_mobility_value` -- the window's mobility position, not a band.
        //
        // So: collect every selected ion of the spectrum, match each to the
        // window whose centre it names, and take the midpoint between adjacent
        // mobility positions as the split. That is a derived boundary, not one
        // the file states, which is why it is only applied when the ions
        // actually resolve to distinct windows.
        struct Ion { double mz; double im; };
        std::vector<Ion> ions;
        for (const auto& prec : s.precursors())
        {
          for (const auto& sel : prec.selected_ions)
          {
            if (sel.selected_ion_mz && sel.ion_mobility_value)
            {
              ions.push_back({*sel.selected_ion_mz, *sel.ion_mobility_value});
            }
          }
        }

        std::vector<SpectrumInfo> here;
        for (const auto& prec : s.precursors())
        {
          const auto& w = prec.isolation_window;
          if (!w.target_mz) { continue; }
          SpectrumInfo one = info;
          one.window.mz_low = *w.target_mz - (w.lower_offset ? *w.lower_offset : 0.0f);
          one.window.mz_high = *w.target_mz + (w.upper_offset ? *w.upper_offset : 0.0f);

          for (const auto& sel : prec.selected_ions)
          {
            if (sel.ion_mobility_lower_limit) { one.window.im_low = *sel.ion_mobility_lower_limit; }
            if (sel.ion_mobility_upper_limit) { one.window.im_high = *sel.ion_mobility_upper_limit; }
          }

          // "lower" and "upper" are the writer's SCAN order, not an ordering of
          // 1/K0: on a timsTOF 1/K0 decreases with scan number, so the file's
          // "lower limit" can be the numerically larger value. Taken verbatim
          // the band test rejects every peak, silently, because the precursor
          // is still placed and still gets a full-length all-zero trace.
          if (one.window.im_low > one.window.im_high)
          {
            std::swap(one.window.im_low, one.window.im_high);
          }
          here.push_back(one);
        }

        // No stated band: derive one from the mobility positions, but only when
        // they actually distinguish the windows. Attributing by nearest centre
        // and splitting at the midpoint is a guess about the instrument, and a
        // wrong guess here rejects real signal -- so it is taken only when each
        // window claims exactly one ion.
        if (here.size() > 1 && ions.size() == here.size())
        {
          std::vector<double> at(here.size(), std::numeric_limits<double>::quiet_NaN());
          std::vector<char> used(ions.size(), 0);
          for (std::size_t k = 0; k < here.size(); ++k)
          {
            const double centre = 0.5 * (here[k].window.mz_low + here[k].window.mz_high);
            std::size_t best = ions.size();
            double best_d = std::numeric_limits<double>::infinity();
            for (std::size_t j = 0; j < ions.size(); ++j)
            {
              if (used[j]) { continue; }
              const double d = std::abs(ions[j].mz - centre);
              if (d < best_d) { best_d = d; best = j; }
            }
            // Half a window's width is generous; beyond that the ion is not
            // naming this window and the whole attribution is abandoned.
            const double tol = 0.5 * (here[k].window.mz_high - here[k].window.mz_low) + 1.0;
            if (best < ions.size() && best_d <= tol) { used[best] = 1; at[k] = ions[best].im; }
          }

          bool all_known = true;
          for (const double v : at) { if (!std::isfinite(v)) { all_known = false; } }
          if (all_known)
          {
            std::vector<std::size_t> order(here.size());
            for (std::size_t k = 0; k < order.size(); ++k) { order[k] = k; }
            std::sort(order.begin(), order.end(),
                      [&](std::size_t a, std::size_t b) { return at[a] < at[b]; });
            // Adjacent bands SHARE their boundary: one window's upper limit is
            // the next one's lower limit, the same double. That is deliberate
            // -- there is no gap to leave between two halves of a midpoint --
            // and it is safe only because the band is half-open, [low, high).
            // See IsolationWindow: with a closed interval a peak landing
            // exactly on a split enters both windows and is counted twice.
            for (std::size_t r = 0; r < order.size(); ++r)
            {
              const std::size_t k = order[r];
              if (!std::isfinite(here[k].window.im_low) || !std::isfinite(here[k].window.im_high) ||
                  here[k].window.im_low >= here[k].window.im_high)
              {
                const double lo = r == 0 ? -std::numeric_limits<double>::infinity()
                                         : 0.5 * (at[order[r - 1]] + at[k]);
                const double hi = r + 1 == order.size() ? std::numeric_limits<double>::infinity()
                                                        : 0.5 * (at[k] + at[order[r + 1]]);
                here[k].window.im_low = lo;
                here[k].window.im_high = hi;
              }
            }
          }
        }

        // A band that is still degenerate means the writer recorded only a
        // midpoint. Read as a band that is exact float equality -- the most
        // aggressive filter possible -- so it means "no band known" instead.
        for (auto& one : here)
        {
          if (!(one.window.im_low < one.window.im_high))
          {
            one.window.im_low = -std::numeric_limits<double>::infinity();
            one.window.im_high = std::numeric_limits<double>::infinity();
          }
          distinct.emplace(quantise(one.window.mz_low), one.window);
          info_.push_back(one);
        }
        if (here.empty()) { info_.push_back(info); }
      }

      // Ascending in retention time, which the extractor relies on. Stated as
      // a sort rather than assumed: a run whose spectra are not written in
      // acquisition order would otherwise produce chromatograms that zig-zag.
      std::stable_sort(info_.begin(), info_.end(),
                       [](const SpectrumInfo& a, const SpectrumInfo& b) {
                         return a.retention_time < b.retention_time;
                       });

      windows_.reserve(distinct.size());
      for (auto& [key, w] : distinct) { windows_.push_back(w); }
    }

    const std::vector<SpectrumInfo>& spectra() const override { return info_; }
    const std::vector<IsolationWindow>& windows() const override { return windows_; }

    void peaks(std::size_t begin, std::size_t end,
               std::vector<SpectrumPeaks>& out) override
    {
      if (begin > end || end > info_.size())
      {
        throw std::out_of_range("peaks(): range outside the run");
      }
      out.assign(end - begin, SpectrumPeaks{});

      // Through the batch entry point, so this call gets faster for free when
      // the reader decodes each row group once instead of once per spectrum.
      std::vector<std::size_t> want;
      want.reserve(end - begin);
      for (std::size_t i = begin; i < end; ++i) { want.push_back(info_[i].index); }
      auto batch = spectra_.get_spectra_batch(want);

      for (std::size_t k = 0; k < batch.size(); ++k)
      {
        auto& dst = out[k];
        const auto& mz = batch[k].mz();
        const auto& intensity = batch[k].intensity();
        // A short intensity array against a long m/z array would silently
        // pair the wrong values, so the arrays are trusted only to their
        // common length and the shortfall is visible as missing peaks.
        const std::size_t n = std::min(mz.size(), intensity.size());
        dst.mz.assign(mz.begin(), mz.begin() + n);
        dst.intensity.assign(intensity.begin(), intensity.begin() + n);

        const auto& im = batch[k].ion_mobility_array();
        if (im.size() >= n) { dst.ion_mobility.assign(im.begin(), im.begin() + n); }
      }
    }

    std::string describe() const override
    {
      return "mzPeak " + filename_ + " (" + std::to_string(info_.size()) + " MS2 spectra, " +
             std::to_string(windows_.size()) + " isolation windows)";
    }

  private:
    std::string filename_;
    MzPeak::Index index_;
    MzPeak::Spectra spectra_;
    std::vector<SpectrumInfo> info_;
    std::vector<IsolationWindow> windows_;
  };

  void SpectrumSource::peaks(std::size_t index, SpectrumPeaks& out)
  {
    std::vector<SpectrumPeaks> one;
    peaks(index, index + 1, one);
    out = one.empty() ? SpectrumPeaks{} : std::move(one.front());
  }

  std::unique_ptr<SpectrumSource> openRun(const std::string& filename)
  {
    return std::make_unique<MzPeakSource>(filename);
  }

} // namespace ODIA
