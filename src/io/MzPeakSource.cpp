// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/SpectrumSource.h>

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

        for (const auto& prec : s.precursors())
        {
          const auto& w = prec.isolation_window;
          if (!w.target_mz) { continue; }
          info.window.mz_low = *w.target_mz - (w.lower_offset ? *w.lower_offset : 0.0f);
          info.window.mz_high = *w.target_mz + (w.upper_offset ? *w.upper_offset : 0.0f);

          // diaPASEF: one frame carries several windows over disjoint mobility
          // ranges, so the mobility limits are part of the window's identity,
          // not a property of the spectrum.
          for (const auto& ion : prec.selected_ions)
          {
            if (ion.ion_mobility_lower_limit)
            {
              info.window.im_low = *ion.ion_mobility_lower_limit;
            }
            if (ion.ion_mobility_upper_limit)
            {
              info.window.im_high = *ion.ion_mobility_upper_limit;
            }
          }
          distinct.emplace(quantise(info.window.mz_low), info.window);
          break;   // one precursor per spectrum in every DIA scheme seen here
        }
        info_.push_back(info);
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
