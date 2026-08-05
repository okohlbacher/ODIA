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
      double fragment_ppm = 20.0;

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
