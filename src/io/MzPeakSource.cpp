// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/MobilityBands.h>
#include <odia/SpectrumSource.h>

#include <limits>
#include <utility>

#include <mzpeak.h>

#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

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
        std::vector<MobilityPosition> ions;
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

        // No stated band: derive one from the mobility positions. The rule and
        // its refusals live in MobilityBands, where they can be tested without
        // a file; the outcome is COUNTED here, because a derivation that
        // quietly does not happen leaves two co-packed windows sharing one
        // merged peak list with nothing to separate them.
        if (here.size() > 1)
        {
          std::vector<IsolationWindow> band(here.size());
          for (std::size_t k = 0; k < here.size(); ++k) { band[k] = here[k].window; }
          const auto result = deriveMobilityBands(ions, band);
          ++derivation_[std::size_t(result)];
          if (result == MobilityBandResult::Derived)
          {
            for (std::size_t k = 0; k < here.size(); ++k) { here[k].window = band[k]; }
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

      // Say so when the band could not be derived. A refusal is the safe
      // outcome -- the windows keep the full mobility range, which admits
      // everything -- but it is not a harmless one: on a run that co-packs
      // windows into one frame it means the co-packed windows are no longer
      // separated at all, and the only visible symptom would be interference
      // that looks like the instrument's.
      std::size_t refused = 0;
      for (std::size_t r = 0; r < derivation_.size(); ++r)
      {
        if (r != std::size_t(MobilityBandResult::Derived) &&
            r != std::size_t(MobilityBandResult::NotNeeded))
        {
          refused += derivation_[r];
        }
      }
      if (refused)
      {
        std::cerr << "warning: no ion-mobility band derived for " << refused << " of "
                  << (refused + derivation_[std::size_t(MobilityBandResult::Derived)])
                  << " co-packed spectra;"
                  << " those windows keep the full mobility range and are not separated ("
                  << toString(MobilityBandResult::TooFewIons) << ": "
                  << derivation_[std::size_t(MobilityBandResult::TooFewIons)] << ", "
                  << toString(MobilityBandResult::Unattributable) << ": "
                  << derivation_[std::size_t(MobilityBandResult::Unattributable)] << ", "
                  << toString(MobilityBandResult::Ambiguous) << ": "
                  << derivation_[std::size_t(MobilityBandResult::Ambiguous)] << ")\n";
      }
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

      // One request per DISTINCT physical spectrum, not one per entry.
      //
      // `info_` holds one entry per (spectrum, isolation window), because a
      // diaPASEF frame co-packs several windows into one physical spectrum and
      // each of them is a separate thing to extract. Those entries share a
      // `index`. mzPeak does no deduplication and no inter-request caching, so
      // asking for the same index twice costs twice: measured on S08,
      // 8.500 ms/request for 4,096 distinct indices against 8.248 ms/request
      // for 2,048 indices asked for twice each -- flat per request, regardless
      // of whether the frame was just decoded.
      //
      // S08 has 17,448 physical spectra and 32,210 entries, so the un-deduped
      // request list asked for every MS2 frame twice and spent ~297 s of its
      // ~594 s decode re-decoding what it already had.
      std::vector<std::size_t> want;
      want.reserve(end - begin);
      // Position in [begin, end) -> which entry of `want`, and therefore of the
      // batch, carries its peaks. `out` is indexed by position in the range,
      // NOT by file index, so this indirection is what keeps the two apart.
      std::vector<std::size_t> slot(end - begin, 0);
      {
        std::unordered_map<std::size_t, std::size_t> first_request;
        first_request.reserve((end - begin) * 2);
        for (std::size_t i = begin; i < end; ++i)
        {
          const auto [it, fresh] = first_request.emplace(info_[i].index, want.size());
          if (fresh) { want.push_back(info_[i].index); }
          slot[i - begin] = it->second;
        }
      }
      auto batch = spectra_.get_spectra_batch(want);

      for (std::size_t k = 0; k < out.size(); ++k)
      {
        const std::size_t s = slot[k];
        if (s >= batch.size()) { continue; }
        auto& dst = out[k];
        const auto& mz = batch[s].mz();
        const auto& intensity = batch[s].intensity();
        // A short intensity array against a long m/z array would silently
        // pair the wrong values, so the arrays are trusted only to their
        // common length and the shortfall is visible as missing peaks.
        const std::size_t n = std::min(mz.size(), intensity.size());
        // Copies, one per entry, rather than a shared buffer: the caller keeps
        // the block alive across the whole match and several threads read it,
        // and a later stage may retain it. Two entries of a co-packed frame
        // hold equal peak lists and are separated by their mobility bands, not
        // by their storage.
        dst.mz.assign(mz.begin(), mz.begin() + n);
        dst.intensity.assign(intensity.begin(), intensity.begin() + n);

        const auto& im = batch[s].ion_mobility_array();
        if (im.size() >= n) { dst.ion_mobility.assign(im.begin(), im.begin() + n); }
      }
    }

    std::string describe() const override
    {
      std::string out = "mzPeak " + filename_ + " (" + std::to_string(info_.size()) +
                        " MS2 spectra, " + std::to_string(windows_.size()) +
                        " isolation windows";
      const std::size_t derived = derivation_[std::size_t(MobilityBandResult::Derived)];
      std::size_t refused = 0;
      for (std::size_t r = 0; r < derivation_.size(); ++r)
      {
        if (r != std::size_t(MobilityBandResult::Derived) &&
            r != std::size_t(MobilityBandResult::NotNeeded))
        {
          refused += derivation_[r];
        }
      }
      if (derived || refused)
      {
        out += ", mobility band derived for " + std::to_string(derived) +
               " co-packed spectra and refused for " + std::to_string(refused);
      }
      return out + ")";
    }

  private:
    std::string filename_;
    MzPeak::Index index_;
    MzPeak::Spectra spectra_;
    std::vector<SpectrumInfo> info_;
    std::vector<IsolationWindow> windows_;
    /// One counter per MobilityBandResult, so a refusal is reported rather
    /// than inferred from chromatograms that came out worse.
    std::array<std::size_t, 5> derivation_{};
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
