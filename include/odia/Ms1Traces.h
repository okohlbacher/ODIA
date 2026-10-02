// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/Library.h>
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
    /// @param threads  matching threads (decode stays on the caller's thread).
    ///        0 and 1 are both serial. The matrix, the retained residual list
    ///        and so its median are identical for every value -- see the loop
    ///        for why, and `odia_ms1_parallel` for the test that holds it to
    ///        that.
    /// @param block_budget_bytes  how much MORE than the serial path the
    ///        parallel match may hold at once: decoded peaks of a match block
    ///        beyond its first 64-frame decode, plus the per-unit residual
    ///        buffers. Further capped at MAX_BLOCK_BYTES, so the held block
    ///        never grows with @p threads. The caller passes what the memory
    ///        budget has left after the matrix; SIZE_MAX means "no budget",
    ///        i.e. MAX_BLOCK_BYTES alone. 0 restricts every match block to one
    ///        64-frame decode -- what the serial path holds -- and still
    ///        matches it on @p threads (at most four units). Output is
    ///        identical for every value: the block boundaries never reach a
    ///        cell or the residual list.
    static Ms1Traces build(const Library& library, SpectrumSource& source,
                           double fragment_ppm, double im_window,
                           double ppm_offset = 0.0,
                           double* observed_ppm_median = nullptr,
                           double isotope_offset_da = 0.0,
                           const std::vector<std::uint8_t>* keep = nullptr,
                           std::vector<std::uint32_t>* kept_indices = nullptr,
                           unsigned threads = 1,
                           std::size_t block_budget_bytes = SIZE_MAX);

    /// Ceiling on the decoded bytes a parallel match block may hold beyond one
    /// 64-frame decode, whatever the budget and -threads say. ~6 serial-size
    /// blocks at the ~11 MiB per frame measured on PXD047793 (fixture, 2026-10-02),
    /// i.e. ~20 sixteen-frame units per match block.
    static constexpr std::size_t MAX_BLOCK_BYTES = std::size_t(4) << 30;

    /// Where the build's wall went. Reported, never used for a decision.
    ///
    /// Exists because nothing split the 5,590 s the reference arm spent here
    /// into decode and match (doc/83 F05): the decode is serial by
    /// construction, so it is the floor a parallel match runs down to, and it
    /// was only ever inferred.
    struct BuildStats
    {
      std::size_t frames = 0;           ///< MS1 spectra matched
      std::size_t blocks = 0;           ///< driver decode calls, 64 frames each
      std::size_t match_blocks = 0;     ///< match blocks (one or more decodes each)
      std::size_t units = 0;            ///< contiguous frame runs matched
      std::size_t peaks = 0;            ///< decoded MS1 peaks, all frames
      std::size_t max_block_bytes = 0;  ///< largest match block's decoded peaks
      /// Most decoded bytes alive at once: a match block plus the decode that
      /// did not fit it and is carried to the next one.
      std::size_t max_held_bytes = 0;
      /// Peak of the per-unit residual buffers' worst case charged to a block.
      std::size_t max_resid_bytes = 0;
      std::size_t block_cap_bytes = 0;  ///< the cap in force (budget vs MAX_BLOCK_BYTES)
      std::size_t carried = 0;          ///< decodes that outgrew their block and moved to the next
      unsigned threads = 1;             ///< most threads any block matched on
      double decode_s = 0.0;            ///< wall in source.ms1Peaks(), driver only
      double match_s = 0.0;             ///< wall in matching, all blocks
      /// The residual prefix behind the median, as kept (at most 2M samples)
      /// and hashed in kept order before the median reorders it: the identity
      /// gate is on the list, since two lists can share a median.
      std::size_t resid_kept = 0;
      std::uint64_t resid_hash = 0;
    };
    const BuildStats& buildStats() const { return stats_; }
    std::string describeBuild() const;

    /// A 64-bit hash of the shape, the time grid and every cell's bit pattern.
    ///
    /// For identity gates on the MATRIX itself rather than through the scorer,
    /// which reads only rows inside a precursor's window and could hide a
    /// difference outside it. A full pass over the matrix (50 GiB on IH1), so
    /// the caller prints it only when asked (ODIA_MS1_CHECKSUM).
    std::uint64_t checksum() const;

  private:
    std::vector<float> times_;
    std::vector<float> values_;      ///< precursor-major, `bins_` per precursor
    std::size_t bins_ = 0;
    BuildStats stats_;
  };

} // namespace ODIA
