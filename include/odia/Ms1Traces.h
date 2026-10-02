// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/Library.h>
#include <odia/SpectrumSource.h>

#include <atomic>
#include <cassert>
#include <cstdint>
#include <memory>
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
    /// How many precursors `at()` answers for: the dense matrix's rows, or,
    /// with a library->row map (@p library_rows of build), the LIBRARY size,
    /// so a caller's `i < precursors()` guard means what it always meant.
    std::size_t precursors() const
    {
      if (!row_of_.empty()) { return row_of_.size(); }
      return bins_ ? values_.size() / bins_ : 0;
    }
    /// Rows actually stored.
    std::size_t rows() const { return bins_ ? values_.size() / bins_ : 0; }
    bool empty() const { return times_.empty() || values_.empty(); }

    /// Monoisotopic intensity of @p precursor at MS1 bin @p b, or 0.
    ///
    /// With a library->row map the index is still a LIBRARY index and goes
    /// through the map. A precursor the map gave no row (no isolation window
    /// covers it, so it is never extracted and never scored) must never be
    /// asked for: debug builds assert it, and every build counts it
    /// (`droppedRowReads()`), returning the 0.0f an all-zero dense row would.
    /// The count is the release-mode guard: nonzero means the drop changed a
    /// read, and the run says so.
    float at(std::size_t precursor, std::size_t b) const
    {
      std::size_t row = precursor;
      if (!row_of_.empty())
      {
        if (precursor >= row_of_.size()) { return 0.0f; }
        row = row_of_[precursor];
        if (row == NO_ROW)
        {
          assert(!"Ms1Traces::at() read a row -ms1_drop_uncovered removed");
          dropped_reads_->fetch_add(1, std::memory_order_relaxed);
          return 0.0f;
        }
      }
      const std::size_t i = row * bins_ + b;
      return i < values_.size() ? values_[i] : 0.0f;
    }

    /// Precursors that have no row because @p library_rows excluded them.
    std::size_t droppedRows() const { return dropped_rows_; }
    /// Reads of such a row since build. Must stay 0; see at().
    std::size_t droppedRowReads() const
    { return dropped_reads_ ? dropped_reads_->load(std::memory_order_relaxed) : 0; }

    /// The MS1 bin nearest a retention time. Binary search; the grid is sorted.
    std::size_t binFor(double rt) const;

    std::size_t footprintBytes() const
    {
      return times_.capacity() * sizeof(float) + values_.capacity() * sizeof(float) +
             row_of_.capacity() * sizeof(std::uint32_t);
    }

    std::string describe() const;

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
    /// @param library_rows  optional per-precursor mask (size precursorCount),
    ///        the scorer's counterpart of @p keep: only masked-in precursors
    ///        get a row, but `at()` keeps taking LIBRARY indices, through an
    ///        explicit library->row map built here and consulted on every
    ///        read -- the scorer never learns that rows moved. A kept row holds
    ///        exactly what the dense matrix's row held: the matrix is a max per
    ///        (precursor, spectrum) over that precursor's own matches, so
    ///        removing other precursors' targets changes no kept cell. Only the
    ///        @p observed_ppm_median diagnostic sees fewer matches. Used by
    ///        -ms1_drop_uncovered with ChromatogramExtractor::windowCoverage.
    ///        Not combinable with @p keep (throws).
    static Ms1Traces build(const Library& library, SpectrumSource& source,
                           double fragment_ppm, double im_window,
                           double ppm_offset = 0.0,
                           double* observed_ppm_median = nullptr,
                           double isotope_offset_da = 0.0,
                           const std::vector<std::uint8_t>* keep = nullptr,
                           std::vector<std::uint32_t>* kept_indices = nullptr,
                           const std::vector<std::uint8_t>* library_rows = nullptr);

  private:
    static constexpr std::uint32_t NO_ROW = UINT32_MAX;

    std::vector<float> times_;
    std::vector<float> values_;      ///< row-major, `bins_` per row
    std::size_t bins_ = 0;
    /// library index -> row, NO_ROW where @p library_rows excluded it. Empty
    /// without a map: row == library index (or == mask row under @p keep).
    std::vector<std::uint32_t> row_of_;
    std::size_t dropped_rows_ = 0;
    /// Shared so the object stays copyable and movable; the count belongs to
    /// the matrix, not to one copy of the handle.
    std::shared_ptr<std::atomic<std::size_t>> dropped_reads_;
  };

} // namespace ODIA
