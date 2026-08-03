// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/Library.h>
#include <odia/SpectrumSource.h>

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
  struct Chromatograms
  {
    /// Retention times, seconds, ascending within each transition's run.
    std::vector<float> retention_time;
    std::vector<float> intensity;

    /// transition index -> [begin, begin + count) into the arrays above.
    std::vector<std::uint32_t> begin;
    std::vector<std::uint32_t> count;

    /// Precursors that no isolation window covered. They cannot be extracted
    /// and are reported rather than silently absent from the output.
    std::size_t precursors_without_window = 0;

    std::size_t points() const { return retention_time.size(); }
    std::size_t footprintBytes() const;
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
      ///
      /// Present because a whole library over a whole run is
      /// (transitions x spectra-per-window) points -- billions -- and nothing
      /// downstream can use that yet. Retention-time restriction, which is the
      /// real answer, needs an iRT-to-RT calibration that does not exist until
      /// the run has been searched once.
      std::size_t max_precursors = 0;

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
    };

    static Chromatograms extract(const Library& library, SpectrumSource& source,
                                 const Options& options, Stats* stats = nullptr);
  };

} // namespace ODIA
