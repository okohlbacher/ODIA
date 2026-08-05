// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/Library.h>
#include <odia/SpectrumSource.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace ODIA
{

  /// Extracted-ion chromatograms for a library's transitions.
  ///
  /// Stored as one flat point array with a CSR index, for the same reason the
  /// library is: a vector-per-transition would be millions of small
  /// allocations, which is exactly the arena fragmentation D3 exists to avoid.
  /// Extracted-ion chromatograms for a library's transitions.
  ///
  /// Two space decisions, both taken from the reference implementation:
  ///
  /// **Retention time is an index, not a number.** A transition's points sit on
  /// its window's acquisition grid, so storing a float per point stored the
  /// same ~1,000 timestamps once per transition. The axis is held once per
  /// window and a transition records where on it its run starts. At 157 M
  /// points that is 630 MB of duplicated timestamps not stored. The deeper
  /// reason is correctness rather than size: recalibrating a run rewrites one
  /// array of ~1,000 doubles and every point follows, instead of rewriting
  /// every structure that cached a time and silently missing one.
  ///
  /// Intensity stays float32 for now; one-byte log quantisation is in the
  /// backlog. It is a further 4x and costs nothing real, because the scale is
  /// per transition rather than global -- 255 log steps over one transition's
  /// own two-decade range is ~2.7% per step.
  ///
  /// The axis is behind an accessor, so callers index by (transition, position)
  /// and never see the representation.
  struct Chromatograms
  {
    /// One retention-time axis per isolation window, seconds, ascending.
    std::vector<std::vector<float>> axes;

    /// transition -> which axis, and where on it point 0 sits.
    std::vector<std::uint32_t> axis_of;
    std::vector<std::uint32_t> axis_begin;

    /// transition index -> [begin, begin + count) into the intensity arrays.
    std::vector<std::uint32_t> begin;
    std::vector<std::uint32_t> count;

    std::vector<float> intensity;

    /// Precursors that no isolation window covered. They cannot be extracted
    /// and are reported rather than silently absent from the output.
    std::size_t precursors_without_window = 0;

    std::size_t points() const { return intensity.size(); }
    std::size_t footprintBytes() const;

    /// Retention time of point @p j of @p transition, in seconds.
    float retentionTime(std::uint32_t transition, std::uint32_t j) const
    {
      const std::uint32_t a = axis_of[transition];
      return axes[a][axis_begin[transition] + j];
    }

  };

  /// Reads a run once and extracts every requested transition from it.
  ///
  /// One forward pass, not one pass per precursor. The cost of decoding a
  /// spectrum is paid once and shared by every transition that lands in it,
  /// which is what makes the extraction scale with the run rather than with
  /// the library. With the current mzPeak reader a spectrum costs ~294 ms to
  /// decode, so this ordering is the difference between one pass and none.
  class ChromatogramExtractor
  {
  public:
    struct Options
    {
      /// Fragment mass tolerance. Relative, because that is how a mass
      /// spectrometer's accuracy behaves; an absolute window would be far too
      /// wide at 200 Th and far too narrow at 1800.
      /// Fragment mass tolerance, as a HALF-width in ppm: a peak matches when
      /// |peak - transition| <= transition * fragment_ppm * 1e-6.
      ///
      /// 15 rather than 10, and only because the window is not yet CENTRED.
      ///
      /// This instrument has a systematic fragment mass offset of about
      /// -10 ppm, measured on S08 as the centroid of the retention-time-
      /// specific excess over a local decoy-cell null: -11.2 ppm on precursors
      /// we recover, -12.6 ppm on those we miss, and present at off-peak times
      /// too, so it is a genuine instrument term and not an artefact of peak
      /// selection. Fraction of true fragments captured:
      ///
      ///     half-width   centred on 0   centred on -10 ppm
      ///        +/- 5        0.22             0.52
      ///        +/-10        0.51             0.81
      ///        +/-15        0.78             0.92
      ///
      /// So narrowing a window that is centred on the wrong place THROWS AWAY
      /// signal: +/-10 about zero keeps half the evidence, and the half it
      /// keeps is the tail rather than the peak. Until `fragment_ppm_offset` is
      /// fitted per run, the safe default is the wider window. With a fitted
      /// offset, +/-10 captures 0.81 at a quarter of the interference area of
      /// +/-20, which is the configuration to aim at.
      double fragment_ppm = 15.0;

      /// Systematic fragment mass offset, ppm, added to the theoretical m/z
      /// before matching. 0 means uncalibrated.
      ///
      /// Should be fitted from the run rather than set by hand -- that is what
      /// per-run mass calibration is for, and it is the single measured
      /// difference most likely to account for the recovery gap. -10 is the
      /// measured value for S08 and is NOT a default, because it is a property
      /// of that instrument and that acquisition.
      double fragment_ppm_offset = 0.0;

      /// Half-width of the ion-mobility window around the PRECURSOR's own
      /// library 1/K0, in 1/K0 units. 0 disables it.
      ///
      /// Distinct from the isolation window's band, which is what
      /// `use_ion_mobility` gates and which only separates co-packed windows.
      /// A frame's band is ~0.40 wide on S08; a precursor occupies ~0.05 of it.
      /// Filtering only by the band therefore admits the entire same-window
      /// mobility axis -- measured as ~8.5x more mobility than the reference
      /// accepts, and the reason a band-only fix bought 1.12x while the
      /// same-window interference it left behind is what dominates.
      ///
      /// 0.025 half-width is ~2.6 sigma on the measured library-vs-observed
      /// agreement (SD 0.019, and 0.0186 on precursors we currently miss, so
      /// not a selection effect).
      double precursor_im_window = 0.025;

      /// How several peaks inside one transition's box become one number.
      ///
      /// `Sum` integrates; `Max` takes the largest. Max was the original and is
      /// wrong for this data: a frame's peak array is the concatenation of
      /// 600-810 TIMS mobility scans, so peaks inside one m/z tolerance are the
      /// same ion across many mobility steps PLUS every co-isolated interferent
      /// at every other step. Max over that returns the interference envelope,
      /// and -- being an order statistic -- acts as a hard threshold rather
      /// than a graded penalty: a peptide below the envelope contributes
      /// nothing at all. That is the shape of the measured failure, where the
      /// bottom four abundance deciles sat exactly at the decoy null.
      ///
      /// Sum is only correct once the box is small. Summing over the old
      /// oversized box makes the trace worse, not better.
      enum class Aggregate { Sum, Max };
      Aggregate aggregate = Aggregate::Sum;

      /// Extract only the first N precursors of the library, 0 for all.
      std::size_t max_precursors = 0;

      /// Maps the library's iRT onto this run's retention time, in seconds:
      /// rt = irt_slope * irt + irt_intercept.
      ///
      /// This is what lets a precursor be extracted around where it should
      /// elute instead of across the whole run, and it is the difference
      /// between 1.7e10 points and something that fits in memory. It cannot be
      /// derived here -- it needs the run to have been searched once -- so the
      /// caller supplies it. A slope of 0 disables the restriction and every
      /// precursor is extracted over the full range, which is the old
      /// behaviour and is kept because it is the only option before a first
      /// pass exists.
      double irt_slope = 0.0;
      double irt_intercept = 0.0;

      /// Half-width of the extraction window, seconds. Sized from the
      /// retention-time residual at a stated quantile -- NOT from its
      /// standard deviation, because a window has to cover the tail it is
      /// meant to catch.
      double rt_window_seconds = 60.0;

      /// Worker threads for the matching. 0 uses the hardware concurrency.
      /// The decode stays serial: the reader's thread-safety is unverified,
      /// and decode is a shared cost per spectrum rather than per transition.
      unsigned threads = 0;

      /// Restrict to this retention-time range, seconds. Both zero means all.
      double rt_low = 0.0;
      double rt_high = 0.0;

      /// Match a peak's ion mobility against the window's limits when the run
      /// carries it. A diaPASEF frame holds several windows over disjoint
      /// mobility ranges, so ignoring this mixes them.
      bool use_ion_mobility = true;

      /// Report progress every this many spectra, 0 to stay quiet.
      std::size_t progress_every = 500;
    };

    struct Stats
    {
      std::size_t spectra_read = 0;
      std::size_t precursors = 0;
      std::size_t transitions = 0;
      std::size_t points = 0;
      std::size_t nonzero_points = 0;
      double decode_seconds = 0.0;
      double match_seconds = 0.0;
      double index_seconds = 0.0;

      /// Precursors whose predicted elution fell outside the run entirely.
      std::size_t outside_rt_range = 0;
      /// Mean transitions live at one cycle, which is what the inverted match
      /// actually costs per spectrum.
      double mean_live_transitions = 0.0;
    };

    static Chromatograms extract(const Library& library, SpectrumSource& source,
                                 const Options& options, Stats* stats = nullptr);
  };

} // namespace ODIA
