// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/MobilityBands.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace ODIA
{

  const char* toString(MobilityBandResult r)
  {
    switch (r)
    {
      case MobilityBandResult::Derived:        return "derived";
      case MobilityBandResult::NotNeeded:      return "not needed";
      case MobilityBandResult::TooFewIons:     return "fewer mobility positions than windows";
      case MobilityBandResult::Unattributable: return "a position too far from any window centre";
      case MobilityBandResult::Ambiguous:      return "two equally good attributions";
    }
    return "unknown";
  }

  namespace
  {
    /// Half a window's width is generous; beyond that the position is not
    /// naming this window.
    double tolerance(const IsolationWindow& w) { return 0.5 * w.width() + 1.0; }

    /// Two distances are the same distance when nothing in the data separates
    /// them. 1e-6 Th is far below the precision any file states a window centre
    /// to, so this catches the genuine tie without inventing one.
    constexpr double TIE = 1e-6;
  }

  MobilityBandResult deriveMobilityBands(const std::vector<MobilityPosition>& ions,
                                         std::vector<IsolationWindow>& windows)
  {
    if (windows.size() < 2) { return MobilityBandResult::NotNeeded; }
    // Collected over every precursor of the spectrum, because the file
    // mis-attaches them all to precursor 0 -- so there can be MORE positions
    // than windows (a precursor carrying no isolation window still contributes
    // one). Demanding equality made one such precursor abandon the whole
    // derivation without a word.
    if (ions.size() < windows.size()) { return MobilityBandResult::TooFewIons; }

    std::vector<double> at(windows.size(), std::numeric_limits<double>::quiet_NaN());
    std::vector<char> claimed(ions.size(), 0);
    std::vector<char> matched(windows.size(), 0);

    for (std::size_t done = 0; done < windows.size(); ++done)
    {
      // The globally smallest remaining distance. That pair is each other's
      // nearest among what is left, which is what makes the matching mutual
      // and independent of the order the windows arrive in.
      double best = std::numeric_limits<double>::infinity();
      std::size_t best_w = windows.size(), best_i = ions.size();
      for (std::size_t w = 0; w < windows.size(); ++w)
      {
        if (matched[w]) { continue; }
        for (std::size_t i = 0; i < ions.size(); ++i)
        {
          if (claimed[i]) { continue; }
          const double d = std::abs(ions[i].mz - windows[w].centre());
          if (d < best) { best = d; best_w = w; best_i = i; }
        }
      }
      if (best_w == windows.size()) { return MobilityBandResult::Unattributable; }

      // Ambiguous only when a pair at the SAME distance CONFLICTS with the one
      // chosen -- the same window with another position, or the same position
      // for another window. Two equally near pairs that share neither index do
      // not compete: two windows each naming their own centre are both at
      // distance 0, and that is the ordinary case, not a tie.
      for (std::size_t w = 0; w < windows.size() && best_w < windows.size(); ++w)
      {
        if (matched[w]) { continue; }
        for (std::size_t i = 0; i < ions.size(); ++i)
        {
          if (claimed[i] || (w == best_w && i == best_i)) { continue; }
          if ((w != best_w) == (i != best_i)) { continue; }   // shares neither, or both
          if (std::abs(std::abs(ions[i].mz - windows[w].centre()) - best) <= TIE)
          {
            return MobilityBandResult::Ambiguous;
          }
        }
      }
      if (best > tolerance(windows[best_w])) { return MobilityBandResult::Unattributable; }
      matched[best_w] = 1;
      claimed[best_i] = 1;
      at[best_w] = ions[best_i].im;
    }

    // Windows in mobility order, so a band's neighbours are the windows either
    // side of it on the mobility axis and not on the m/z axis.
    std::vector<std::size_t> order(windows.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return at[a] < at[b]; });

    // Adjacent bands SHARE their boundary: one window's upper limit is the
    // next one's lower limit, the same double. There is no gap to leave
    // between two halves of a midpoint, and none is needed, because the band
    // is half-open -- see IsolationWindow.
    for (std::size_t r = 0; r < order.size(); ++r)
    {
      const std::size_t k = order[r];
      // A stated band wins. Only a missing or degenerate one is derived.
      if (std::isfinite(windows[k].im_low) && std::isfinite(windows[k].im_high) &&
          windows[k].im_low < windows[k].im_high)
      {
        continue;
      }
      windows[k].im_low = r == 0 ? -std::numeric_limits<double>::infinity()
                                 : 0.5 * (at[order[r - 1]] + at[k]);
      windows[k].im_high = r + 1 == order.size() ? std::numeric_limits<double>::infinity()
                                                 : 0.5 * (at[k] + at[order[r + 1]]);
    }
    return MobilityBandResult::Derived;
  }

} // namespace ODIA
