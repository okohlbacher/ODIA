// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace ODIA
{

  /// One DIA isolation window, as the run actually acquired it.
  ///
  /// Held by value rather than by scan index because that is the question
  /// extraction asks: "which spectra could contain a precursor at this m/z".
  /// A source that can answer it without scanning should be free to.
  struct IsolationWindow
  {
    double mz_low = 0.0;
    double mz_high = 0.0;

    /// Ion-mobility bounds, for a run that separates by it. A run without
    /// mobility reports the full range rather than a sentinel, so callers
    /// never have to branch on "does this run have mobility".
    double im_low = -std::numeric_limits<double>::infinity();
    double im_high = std::numeric_limits<double>::infinity();

    bool contains(double mz) const { return mz >= mz_low && mz <= mz_high; }

    /// Windows overlap in most schemes, so a precursor belongs to several.
    bool overlaps(const IsolationWindow& other) const
    {
      return mz_low <= other.mz_high && other.mz_low <= mz_high;
    }

    double centre() const { return 0.5 * (mz_low + mz_high); }
    double width() const { return mz_high - mz_low; }
  };

  /// The peaks of one spectrum, as parallel arrays.
  ///
  /// Struct-of-arrays for the same reason the library is (D3): extraction
  /// binary-searches m/z and then reads one intensity, and an array of pairs
  /// would pull an intensity into cache for every m/z compared.
  struct SpectrumPeaks
  {
    std::vector<double> mz;         ///< ascending
    std::vector<float> intensity;

    /// Ion mobility per peak, empty when the run has none. Not per-peak
    /// optional: either the run separates by mobility or it does not.
    std::vector<float> ion_mobility;

    std::size_t size() const { return mz.size(); }
    bool empty() const { return mz.empty(); }
    bool hasIonMobility() const { return !ion_mobility.empty(); }
    void clear() { mz.clear(); intensity.clear(); ion_mobility.clear(); }
  };

  /// What a run says about one MS2 spectrum without decoding its peaks.
  ///
  /// Separated from the peaks because the costs differ by four orders of
  /// magnitude: mzPeak serves metadata faster than mzML and peaks ~1000x
  /// slower, so a caller that only needs to plan must be able to plan.
  struct SpectrumInfo
  {
    std::size_t index = 0;          ///< position in the run, for batch requests
    double retention_time = 0.0;    ///< seconds
    IsolationWindow window;
    std::uint8_t ms_level = 2;
  };

  /// A DIA run, reached through what extraction needs rather than through a
  /// particular file format.
  ///
  /// Deliberately not mzML-shaped. Peaks are requested for a *range* of spectra
  /// in one call, because that is a columnar read, and forcing it through a
  /// per-spectrum accessor is what makes the mzPeak reader 1000x slower than
  /// parsing mzML today (R2 in doc/03-mzpeak-streaming-requirements.md). A
  /// format that can serve a range cheaply must not be penalised for it here.
  class SpectrumSource
  {
  public:
    virtual ~SpectrumSource() = default;

    /// Every MS2 spectrum's metadata, ascending in retention time.
    ///
    /// Materialised whole: the largest run here is 307,590 spectra, which is a
    /// few MB of this, and every planning decision needs all of it.
    virtual const std::vector<SpectrumInfo>& spectra() const = 0;

    /// The distinct isolation windows the run cycles through.
    ///
    /// Derived from the spectra rather than from the stated method, because a
    /// run's actual window scheme and its method routinely differ.
    virtual const std::vector<IsolationWindow>& windows() const = 0;

    /// Peaks for spectra [begin, end), in order. @p out is resized to the
    /// range. Requesting a range rather than a spectrum is the point.
    virtual void peaks(std::size_t begin, std::size_t end,
                       std::vector<SpectrumPeaks>& out) = 0;

    /// One spectrum, in terms of the range form, so a source only has to get
    /// the range right.
    void peaks(std::size_t index, SpectrumPeaks& out);

    virtual std::string describe() const = 0;

    /// Total peaks in the run when the source knows it cheaply, else 0. For
    /// reporting only -- never for allocation.
    virtual std::size_t peakCount() const { return 0; }
  };

  /// Open a run, choosing the implementation by extension.
  ///
  /// mzML today, mzPeak when its reader can serve peaks at a usable rate. The
  /// data directory carries both for the same three runs, so one is the
  /// reference the other is checked against.
  std::unique_ptr<SpectrumSource> openRun(const std::string& filename);

} // namespace ODIA
