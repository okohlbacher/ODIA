// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/Library.h>
#include <odia/ResourceProbe.h>
#include <odia/SpectrumSource.h>

#include <cstdint>
#include <string>
#include <vector>

namespace ODIA
{

  /// The MS1 precursor trace for every precursor, on the run's MS1 time grid.
  ///
  /// Kept SEPARATE from the MS2 extraction rather than threaded through it. The
  /// streaming extractor is organised around isolation windows and MS2 cycles;
  /// MS1 has neither, and forcing it through that loop would complicate the one
  /// piece of this pipeline whose memory behaviour is understood and measured.
  /// The MS1 arm is small enough to stand alone: one value per precursor per MS1
  /// spectrum, against ~12 fragments per precursor per MS2 cycle.
  ///
  /// WHY THIS EXISTS AT ALL, measured 2026-08-08 on IH1 + `v6_50k` (the search
  /// benchmark, 670 true of 50,000, base 13.4 per 1000):
  ///
  ///     MS1/MS2 co-elution correlation    13.7x enrichment in the top bin,
  ///                                       1.4-1.5x in bulk
  ///
  /// against 0.5x-1.7x for every PRESENCE statistic tried (fragment depth,
  /// mobility-sliced depth, qualifying spectra, total matches, MS1 isotope
  /// depth, MS1 intensity). It survives stratification by MS1 intensity, so it
  /// is shape and not brightness.
  ///
  /// It is the first genuinely ORTHOGONAL evidence in the scorer: all fifteen
  /// existing sub-scores read MS2 fragment traces, so one co-eluting interferent
  /// corrupts them together. This asks whether the PRECURSOR rises and falls
  /// with its fragments, which is a different measurement with an independent
  /// failure mode.
  ///
  /// NOT stored: `ms1_iso` (isotope presence, saturated at 99.8% of precursors)
  /// and `ms1_max` (intensity, an abundance proxy -- and the benchmark is
  /// DIA-NN's confident set, which rewards predicting what DIA-NN finds rather
  /// than what is true). See doc/13.
  class Ms1Traces
  {
  public:
    /// Retention times of the MS1 grid, ascending, seconds.
    const std::vector<float>& times() const { return times_; }

    std::size_t bins() const { return times_.size(); }
    std::size_t precursors() const { return bins_ ? values_.size() / bins_ : 0; }
    bool empty() const { return times_.empty() || values_.empty(); }

    /// Monoisotopic intensity of @p precursor at MS1 bin @p b, or 0.
    float at(std::size_t precursor, std::size_t b) const
    {
      const std::size_t i = precursor * bins_ + b;
      return i < values_.size() ? values_[i] : 0.0f;
    }

    /// The MS1 bin nearest a retention time. Binary search; the grid is sorted.
    std::size_t binFor(double rt) const;

    std::size_t footprintBytes() const
    {
      return times_.capacity() * sizeof(float) + values_.capacity() * sizeof(float);
    }

    std::string describe() const;

    /// Where build() spends its time, stage by stage.
    ///
    /// It was 16.3% of a full IH1 run's wall (doc/83) and is serial, and
    /// the obvious fix -- parallelise over frames (doc/83 F05) -- only pays if
    /// the time is in the MATCH and not in the driver's decode, which the
    /// single "in X s" line could not say. `alloc` is the zero-fill of the dense
    /// matrix (tens of GiB on a full run, so its page faults are a term of their own),
    /// `index` the sorted target table, `decode` the source.ms1Peaks() calls,
    /// `match` the per-peak lookup and max, `median` the residual's
    /// nth_element.
    struct BuildCost
    {
      StageCost alloc, index, decode, match, median;
    };

    /// Build from a run: monoisotopic intensity per precursor per MS1 spectrum.
    ///
    /// @p im_window restricts a match to the precursor's mobility slice. On
    /// diaPASEF this matters -- an MS1 frame is a merged stack of TIMS scans, so
    /// without it the trace collects whatever shares the m/z anywhere in the
    /// mobility range.
    /// @param ppm_offset  the fitted mass deviation to centre the search on,
    ///        in ppm, exactly as the fragment axis is centred. 0 reproduces the
    ///        old uncalibrated behaviour.
    /// @param observed_ppm_median  optional out-param receiving the median
    ///        residual against the CALIBRATED target, so the applied offset can
    ///        be checked rather than assumed: a correct offset centres it on 0.
    /// @param isotope_offset_da  isotope index times the neutron mass delta,
    ///        divided by the precursor's own charge at build time: the target
    ///        m/z becomes `mz + isotope_offset_da / charge`, BEFORE the ppm
    ///        calibration scaling. 0 is the monoisotopic trace. Used by the
    ///        -out_ms1_iso cohort export (M+1, M+2 evidence channels).
    /// @param keep  optional per-precursor mask (size precursorCount). When
    ///        given, only masked-in precursors get a dense row -- the matrix is
    ///        `kept x bins` instead of `np x bins`, which is what makes a 3x
    ///        isotope build affordable (the full-library dense matrix is 77.6
    ///        GB at 4.99M precursors and is charged before extraction begins).
    ///        With @p keep, `at()` takes ROW indices, not library indices;
    ///        @p kept_indices receives the library index of each row.
    /// @param cost  optional out-param receiving the stage costs; read-only
    ///        instrumentation, the traces do not depend on it.
    static Ms1Traces build(const Library& library, SpectrumSource& source,
                           double fragment_ppm, double im_window,
                           double ppm_offset = 0.0,
                           double* observed_ppm_median = nullptr,
                           double isotope_offset_da = 0.0,
                           const std::vector<std::uint8_t>* keep = nullptr,
                           std::vector<std::uint32_t>* kept_indices = nullptr,
                           BuildCost* cost = nullptr);

  private:
    std::vector<float> times_;
    std::vector<float> values_;      ///< precursor-major, `bins_` per precursor
    std::size_t bins_ = 0;
  };

} // namespace ODIA
