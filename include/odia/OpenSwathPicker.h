// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <vector>

namespace ODIA
{
  /// OpenSWATH's own chromatographic peak picker, as a third transition finder.
  ///
  /// ODIA has two pickers of its own: an AMPLITUDE detector (local maxima on the
  /// summed trace) and a CO-ELUTION detector ported from DIA-NN's
  /// `Searcher::peaks`, which only calls a position a peak when the fragments
  /// already correlate there. Replacing the first with the second fixed
  /// retention time, FDR, recovery and memory at once.
  ///
  /// This is the reference implementation of the third approach, and it belongs
  /// in the comparison for a specific reason: `doc/07-scoring-plan.md` step 2
  /// mandates cross-checking against OpenSWATH on the same input, and that step
  /// has never been done. Two sub-scores were later found worthless
  /// (`USABLE_FRAGMENTS` constant 12.000 for a week, `IM_DELTA` all-NaN) --
  /// precisely what an independent implementation is there to catch.
  ///
  /// OpenSWATH picks on the SUMMED chromatogram: smooth (Savitzky-Golay or
  /// Gaussian), find local maxima, then set boundaries by signal-to-noise
  /// rather than by a fixed fraction of apex height. Co-elution enters
  /// afterwards, as `xcorr_shape`/`xcorr_coelution` scores on an
  /// already-chosen peak. That ordering is the substantive difference from the
  /// co-elution detector, not the smoothing kernel.
  ///
  /// Uses OpenMS 3.6's `PeakPickerChromatogram` unmodified -- the project uses
  /// OpenMS as a library and does not patch it. (The class was `PeakPickerMRM`
  /// before OpenMS 3.x.)
  struct OpenSwathCandidate
  {
    std::size_t apex = 0, left = 0, right = 0;
    double apex_value = 0.0;
  };

  /// Pick @p total (the summed trace) at the retention times @p rt.
  ///
  /// @param sn_threshold      signal-to-noise below which a peak is discarded
  /// @param use_gauss         Gaussian smoothing instead of Savitzky-Golay
  /// @param peak_width        expected peak width in seconds, or <=0 to disable
  ///                          the width-based filtering
  /// @param max_candidates    cap, applied after ranking by apex intensity
  std::vector<OpenSwathCandidate> pickOpenSwath(const std::vector<double>& total,
                                                const std::vector<float>& rt,
                                                double sn_threshold, bool use_gauss,
                                                double peak_width,
                                                std::size_t max_candidates);

} // namespace ODIA
