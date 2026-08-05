// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// The part of "probe the run before extracting" that is the same whichever axis
// is being calibrated.
//
// Two calibrations now sample a run before the first pass -- the fragment mass
// error (MassCalibration) and the 1/K0 prediction error (MobilityCalibration).
// They measure different things through different filters, but they reach their
// data the same way: cut the run into acquisition cycles, draw a stratified
// random subset of those, take a spread of library precursors, and give each one
// the isolation window that transmits it best. That machinery has no opinion
// about ppm or 1/K0, and two copies of it would drift -- a second copy of "a
// cycle boundary is where a window repeats" that disagrees with the first is a
// bug nothing would catch.
//
// Internal to src/calibrate. Not installed: it is a shared implementation detail
// of two classes, not part of ODIA's interface.
#pragma once

#include <odia/Library.h>
#include <odia/SpectrumSource.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <utility>
#include <vector>

namespace ODIA::RunProbe
{
  inline constexpr std::uint32_t NO_WINDOW = std::numeric_limits<std::uint32_t>::max();

  /// Which isolation window each spectrum belongs to, by exact bounds match.
  inline std::vector<std::uint32_t> windowOfSpectrum(const SpectrumSource& source)
  {
    const auto& windows = source.windows();
    const auto& info = source.spectra();
    std::vector<std::uint32_t> out(info.size(), NO_WINDOW);
    for (std::size_t si = 0; si < info.size(); ++si)
    {
      for (std::size_t w = 0; w < windows.size(); ++w)
      {
        if (std::abs(info[si].window.mz_low - windows[w].mz_low) < 1e-6 &&
            std::abs(info[si].window.mz_high - windows[w].mz_high) < 1e-6)
        {
          out[si] = static_cast<std::uint32_t>(w);
          break;
        }
      }
    }
    return out;
  }

  /// The run's acquisition cycles as half-open spectrum ranges.
  ///
  /// A DIA cycle is a contiguous run of spectra covering each window once, so a
  /// cycle boundary is where a window repeats. Contiguity is what makes the
  /// sample cheap: one cycle is ONE range request, and every window is probed at
  /// the same retention time.
  inline std::vector<std::pair<std::size_t, std::size_t>>
  cycles(const SpectrumSource& source, const std::vector<std::uint32_t>& window_of)
  {
    std::vector<std::pair<std::size_t, std::size_t>> out;
    const std::size_t n_windows = source.windows().size();
    if (n_windows == 0) { return out; }
    std::vector<char> seen(n_windows, 0);
    std::size_t begin = 0;
    for (std::size_t si = 0; si < window_of.size(); ++si)
    {
      const std::uint32_t w = window_of[si];
      if (w == NO_WINDOW) { continue; }
      if (seen[w])
      {
        out.emplace_back(begin, si);
        std::fill(seen.begin(), seen.end(), 0);
        begin = si;
      }
      seen[w] = 1;
    }
    if (begin < window_of.size()) { out.emplace_back(begin, window_of.size()); }
    return out;
  }

  /// STRATIFIED RANDOM draw of @p want cycle indices out of @p available.
  ///
  /// The run is cut into `want` equal strata and one cycle is drawn from each.
  /// Guaranteed coverage of the whole gradient -- so a drift is visible -- without
  /// the aliasing a fixed stride would risk against a periodic acquisition, and a
  /// DIA run is periodic by construction. A PREFIX measures one stretch of the
  /// gradient and cannot see a drift at all.
  ///
  /// Seeded, not clock-derived: a calibration that returns a different number
  /// each run cannot be checked against a previous one, and "the estimate moved"
  /// would be indistinguishable from "the sample moved".
  inline std::vector<std::size_t> stratifiedCycles(std::size_t available, std::size_t want,
                                                   std::uint64_t seed)
  {
    std::vector<std::size_t> chosen;
    want = std::min(want, available);
    if (want == 0) { return chosen; }
    chosen.reserve(want);
    std::mt19937_64 rng(seed);
    for (std::size_t b = 0; b < want; ++b)
    {
      const std::size_t lo = b * available / want;
      const std::size_t hi = std::max(lo + 1, (b + 1) * available / want);
      std::uniform_int_distribution<std::size_t> pick(lo, hi - 1);
      chosen.push_back(pick(rng));
    }
    std::sort(chosen.begin(), chosen.end());
    return chosen;
  }

  /// STRATIFIED RANDOM draw of @p blocks CONTIGUOUS runs of @p length cycles.
  ///
  /// The same stratification as above, drawing a short BLOCK rather than a
  /// single cycle. What a block buys is time: consecutive cycles are ~1.4 s
  /// apart on a diaPASEF run and a chromatographic peak is tens of seconds
  /// wide, so a precursor that is really there is present in every cycle of a
  /// block that overlaps its elution, while a chance coincidence is not. That
  /// is the only pre-identification handle on "is this signal real" that does
  /// not need a retention-time map, and MobilityCalibration cannot work without
  /// one -- measured, see its header.
  ///
  /// Returns the START cycle of each block; blocks never overlap and are
  /// returned ascending.
  inline std::vector<std::size_t> stratifiedBlocks(std::size_t available, std::size_t blocks,
                                                   std::size_t length, std::uint64_t seed)
  {
    std::vector<std::size_t> out;
    if (available == 0 || blocks == 0 || length == 0) { return out; }
    length = std::min(length, available);
    blocks = std::min(blocks, available / length);
    if (blocks == 0) { return out; }
    std::mt19937_64 rng(seed);
    for (std::size_t b = 0; b < blocks; ++b)
    {
      const std::size_t lo = b * available / blocks;
      const std::size_t hi = std::max(lo + 1, (b + 1) * available / blocks);
      // Keep the block inside its own stratum where it fits, so two blocks can
      // never overlap and re-measure the same cycles.
      const std::size_t last = hi > lo + length ? hi - length : lo;
      std::uniform_int_distribution<std::size_t> pick(lo, last);
      std::size_t start = pick(rng);
      if (start + length > available) { start = available - length; }
      out.push_back(start);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
  }

  /// The most CENTRAL isolation window containing @p mz, or NO_WINDOW.
  ///
  /// An edge precursor is transmitted with reduced efficiency, so when several
  /// windows contain it the central one is the one whose spectra actually hold
  /// its fragments.
  inline std::uint32_t centralWindow(const std::vector<IsolationWindow>& windows, double mz)
  {
    std::uint32_t best = NO_WINDOW;
    double best_margin = -1.0;
    for (std::size_t w = 0; w < windows.size(); ++w)
    {
      if (!windows[w].contains(mz)) { continue; }
      const double margin = std::min(mz - windows[w].mz_low, windows[w].mz_high - mz);
      if (margin > best_margin) { best_margin = margin; best = static_cast<std::uint32_t>(w); }
    }
    return best;
  }

  /// Up to @p want target precursors, spread evenly through the library.
  ///
  /// A stride rather than the first N: the library may be sorted by m/z, and
  /// taking a prefix would then measure the bottom of the mass range only --
  /// which is precisely the axis a shape fit has to see.
  ///
  /// @param require_im keep only precursors whose library 1/K0 is finite. A
  ///        precondition for calibrating the mobility axis; irrelevant to the
  ///        mass axis, which is why it is a parameter and not a rule.
  /// @param decoys take the library's DECOY precursors instead of its targets,
  ///        which is how MobilityCalibration builds a null that is made of real
  ///        peptide-shaped fragment masses rather than of an m/z shift.
  inline std::vector<std::uint32_t> samplePrecursors(const Library& library, std::size_t want,
                                                     std::size_t min_transitions,
                                                     bool require_im, bool decoys = false)
  {
    const auto& p = library.precursors();
    std::vector<std::uint32_t> eligible;
    eligible.reserve(library.precursorCount());
    for (std::size_t i = 0; i < library.precursorCount(); ++i)
    {
      if (static_cast<bool>(p.decoy[i]) != decoys) { continue; }
      if (p.mz[i] == MZ_INVALID) { continue; }
      if (p.transition_count[i] < min_transitions) { continue; }
      if (require_im && !std::isfinite(p.im[i])) { continue; }
      eligible.push_back(static_cast<std::uint32_t>(i));
    }
    std::vector<std::uint32_t> out;
    if (eligible.empty()) { return out; }
    // 0 means "all of them". A probe whose cost is the spectra it decodes, not
    // the queries it evaluates, has no reason to subsample the library.
    want = want == 0 ? eligible.size() : std::min(want, eligible.size());
    out.reserve(want);
    for (std::size_t k = 0; k < want; ++k)
    {
      out.push_back(eligible[k * eligible.size() / want]);
    }
    return out;
  }

} // namespace ODIA::RunProbe
