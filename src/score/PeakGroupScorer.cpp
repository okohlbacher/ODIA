// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/PeakGroupScorer.h>
#include <mutex>
#include <odia/scoring/PercolatorEngine.h>

#include <odia/scoring/gbt.h>
#include <odia/scoring/lda.h>
#include <odia/scoring/score.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numeric>
#include <sstream>

#include <atomic>
#include <chrono>
#include <stdexcept>

namespace ODIA
{
namespace
{
  // MS1_COELUTION diagnosis -- see the block in scoreCandidate. Removed once the
  // question is answered; until then these are the only evidence about WHY the
  // one feature measured to discriminate arrives constant.
  std::atomic<std::size_t> ms1_unavailable{0};   ///< options.ms1 null/empty, or index out of range
  std::atomic<std::size_t> ms1_no_signal{0};     ///< every MS1 point in the candidate window is 0
  std::atomic<std::size_t> ms1_too_short{0};     ///< fewer than 5 cycles
  std::atomic<std::size_t> ms1_flat{0};          ///< MS1 leg has zero variance -> Pearson undefined
  std::atomic<std::size_t> ms1_ok{0};            ///< a finite correlation was produced
  std::atomic<std::size_t> ms1_bins_spanned{0};  ///< sum of distinct MS1 bins per candidate
  std::atomic<std::size_t> ms1_spans{0};         ///< candidates contributing to the above
  std::atomic<std::size_t> ms1_null{0};          ///< options.ms1 was a null pointer
  std::atomic<std::size_t> ms1_empty{0};         ///< traces object present but empty
  std::atomic<std::size_t> ms1_index_oob{0};     ///< precursor index past the matrix
  std::atomic<std::size_t> ms1_degenerate_span{0};  ///< hi <= lo, candidate spans nothing

  /// Distinct MS1 bins a candidate's cycles map onto. If this is ~1 the MS1 grid
  /// is too coarse for the candidate and no correlation is definable.
  std::size_t distinctBins_(const PrecursorChromatogram& c, std::size_t lo, std::size_t hi,
                            const Ms1Traces& ms1)
  {
    std::size_t n = 0, prev = static_cast<std::size_t>(-1);
    for (std::size_t j = lo; j <= hi; ++j)
    {
      const std::size_t b = ms1.binFor(c.retentionTime(static_cast<std::uint32_t>(j)));
      if (b != prev) { ++n; prev = b; }
    }
    return n;
  }
}


  namespace
  {
    /// The summed trace over a precursor's transitions, on the cycle axis.
    ///
    /// Summing rather than taking any single transition: a peak group is a
    /// coelution of all of them, and picking on one transition would find that
    /// transition's noise as readily as the precursor's peak.
    std::vector<double> summedTrace(const PrecursorChromatogram& c,
                                    std::size_t points,
                                    const std::vector<double>* weights = nullptr)
    {
      std::vector<double> total(points, 0.0);
      for (std::uint32_t t = 0; t < c.transition_count; ++t)
      {
        const std::uint32_t n = c.pointCount(t);
        if (n == 0) { continue; }
        const float* at = c.trace(t);
        const double w = weights ? (*weights)[t] : 1.0;
        for (std::uint32_t i = 0; i < n && i < points; ++i)
        {
          total[i] += w * at[i];
        }
      }
      return total;
    }

    /// The sum of each transition standardised by its own median and MAD.
    ///
    /// MAD rather than standard deviation because a trace containing a real
    /// peak has that peak in its own noise estimate, and the SD would be
    /// inflated by exactly the feature being looked for. MAD is unmoved by a
    /// few large points, so a strong fragment does not suppress itself.
    /// Gate C: the co-elution evidence statistic.
    ///
    /// M_p = max over cycles of a smoothed sum, ACROSS transitions, of each
    /// transition's variance-stabilised z-score. Three properties the two
    /// previous gates lacked:
    ///
    ///  * summing across transitions AT THE SAME CYCLE makes co-elution a
    ///    requirement. Gate B counted excursions anywhere in any trace, so 12
    ///    traces x 200 cycles gave ~2400 independent chances and it admitted
    ///    81.5% of pure noise -- worse than the 28.2% of the gate it replaced.
    ///  * the maximum is over ~200 cycles, not ~2400 trace-cycle pairs, so the
    ///    multiple-testing burden is an order of magnitude smaller.
    ///  * sqrt() stabilises the variance of counting noise, so the MAD estimate
    ///    is not dominated by the few bright cycles.
    ///
    /// The threshold is NOT set here. It is calibrated from the decoy null,
    /// because assuming a Gaussian null is exactly what made Gate B's 3 sigma
    /// mean 81% instead of 0.1%.
    /// @param contributing if non-null, receives how many transitions actually
    ///        contributed a z-scored trace. The statistic is a SUM of robust
    ///        z-scores over those transitions, so its null scale is
    ///        sqrt(contributing / w) -- which is what lets a threshold be set
    ///        from a noise model instead of from a run-wide decoy quantile.
    double coelutionEvidence(const PrecursorChromatogram& c, std::size_t points,
                             std::size_t half, std::size_t* contributing = nullptr)
    {
      if (points == 0) { return 0.0; }
      std::vector<double> s(points, 0.0);
      std::vector<double> y, scratch;
      for (std::uint32_t t = 0; t < c.transition_count; ++t)
      {
        const std::uint32_t n = c.pointCount(t);
        if (n == 0) { continue; }
        const float* at = c.trace(t);
        y.assign(n, 0.0);
        for (std::uint32_t i = 0; i < n; ++i)
        { y[i] = std::sqrt(std::max(0.0, double(at[i]))); }

        scratch = y;
        std::sort(scratch.begin(), scratch.end());
        const double median = scratch[scratch.size() / 2];
        for (auto& v : scratch) { v = std::abs(v - median); }
        std::sort(scratch.begin(), scratch.end());
        const double mad = scratch[scratch.size() / 2];
        // A transition with no spread contributes NOTHING rather than being
        // given an invented scale. It cannot vote for or against co-elution.
        if (!(mad > 0.0)) { continue; }
        if (contributing != nullptr) { ++*contributing; }
        const double scale = 1.0 / (1.4826 * mad);
        for (std::uint32_t i = 0; i < n && i < points; ++i)
        { s[i] += (y[i] - median) * scale; }
      }

      // Smooth, so a single bright cycle in one transition cannot carry the
      // precursor: a real peak spans several cycles.
      double best = 0.0;
      const std::size_t w = 2 * half + 1;
      for (std::size_t i = 0; i < points; ++i)
      {
        double acc = 0.0;
        std::size_t used = 0;
        for (std::size_t j = (i > half ? i - half : 0);
             j <= i + half && j < points; ++j) { acc += s[j]; ++used; }
        if (used) { best = std::max(best, acc / double(w)); }
      }
      return best;
    }

    /// Summed noise-normalised trace, and how many transitions actually rose
    /// above `sigma` at any point.
    ///
    /// The count is the point. Summing median-subtracted traces gives a
    /// ZERO-MEAN quantity for a precursor with no real peak, so testing that
    /// sum against zero tests the sign of noise and admits about half of all
    /// absent precursors. Counting transitions with a positive excursion tests
    /// whether anything is actually there.
    std::vector<double> noiseNormalisedTrace(const PrecursorChromatogram& c,
                                             std::size_t points,
                                             double sigma = 0.0,
                                             std::size_t* excursions = nullptr)
    {
      std::vector<double> total(points, 0.0);
      std::vector<double> scratch;
      for (std::uint32_t t = 0; t < c.transition_count; ++t)
      {
        const std::uint32_t n = c.pointCount(t);
        if (n == 0) { continue; }

        const float* at = c.trace(t);
        scratch.assign(at, at + n);
        std::sort(scratch.begin(), scratch.end());
        const double median = scratch[scratch.size() / 2];
        for (auto& v : scratch) { v = std::abs(v - median); }
        std::sort(scratch.begin(), scratch.end());
        const double mad = scratch[scratch.size() / 2];

        // A trace with no variation contributes nothing rather than dividing by
        // a floor and injecting a scaled copy of its own rounding.
        //
        // But ZERO MAD WITH A RISE ABOVE THE MEDIAN IS SIGNAL, and the cleanest
        // kind: a peak on a flat baseline. Skipping it outright made the
        // excursion count zero for a noiseless trace, so the new gate rejected
        // a perfect peak -- caught by picker_prefers_coelution_over_amplitude,
        // whose fixture is deliberately noise-free. The old sum-based gate
        // tolerated this by accident; an explicit signal test must not.
        if (!(mad > 0.0))
        {
          // MAD zero means the trace has no measurable spread, which happens
          // two very different ways and they must not be treated alike:
          //
          //   a clean peak on a flat baseline -- the strongest signal there is
          //   a SPARSE trace that is mostly zeros -- the weakest
          //
          // Counting "any rise above the median" made the second qualify, and
          // in real chromatograms almost every trace is mostly zeros, so the
          // gate stopped rejecting anything at all: measured on Astral it fired
          // 0 times and let 11.5x more precursors through than the old one.
          // Unbiased, but not a test.
          //
          // With no spread to scale by, the only defensible bar is a RELATIVE
          // one: the peak must stand well above the baseline it sits on. A
          // clean Gaussian on a flat baseline clears it easily; a sparse trace
          // whose single nonzero point is a lone count does not.
          if (excursions && sigma > 0.0)
          {
            const float* mx = std::max_element(at, at + n);
            if (mx)
            {
              // Requiring the same sigma-multiple of the baseline keeps one
              // knob rather than inventing a second.
              const double floor_level = std::max(median, 1.0);
              if (double(*mx) >= floor_level * (1.0 + sigma)) { ++*excursions; }
            }
          }
          continue;
        }
        const double scale = 1.0 / (1.4826 * mad);   // MAD -> sigma for normal noise
        double peak = 0.0;
        for (std::uint32_t i = 0; i < n && i < points; ++i)
        {
          const double z = (at[i] - median) * scale;
          total[i] += z;
          if (z > peak) { peak = z; }
        }
        // Already in sigma units, so the threshold is read directly.
        if (excursions && sigma > 0.0 && peak >= sigma) { ++*excursions; }
      }
      return total;
    }

    /// Effective number of fragments a weight vector really spreads across,
    /// 1 / sum(w^2) on normalised weights. Reported rather than assumed: a
    /// dominant library fragment at 0.6 caps this below 2.7 however many
    /// transitions are nominally present, which is not multi-fragment picking.
    double effectiveFragments(const std::vector<double>& w)
    {
      double sum = 0.0, sq = 0.0;
      for (const double v : w) { sum += v; sq += v * v; }
      if (!(sum > 0.0) || !(sq > 0.0)) { return 0.0; }
      const double norm = sq / (sum * sum);
      return norm > 0.0 ? 1.0 / norm : 0.0;
    }

    /// A robust local background for one trace: the median of its points
    /// outside the candidate, falling back to the low quantile of the whole
    /// trace when the candidate spans everything.
    /// Median of the trace in the FLANKS beside a peak group.
    ///
    /// It used to take the median of everything outside [lo,hi], which is a
    /// whole-trace estimate wearing a local name. Measured over 30,586
    /// per-fragment comparisons: a flank estimate and a whole-trace estimate
    /// differ by more than 20% of the local level for **48.7%** of fragments
    /// (median 19.3%, p90 84%), with a median difference of zero -- so the
    /// whole-trace version was not biased, it was noisy, and that variance went
    /// straight into `corrected[]` and from there into LIBRARY_CORR,
    /// LIBRARY_DOTPROD, INTENSITY_SCORE, LIBRARY_RMSD and the RT_SPREAD
    /// weights. Chromatographic background drifts across a +-90 s window; the
    /// baseline under a peak is the baseline BESIDE it, not the run's average.
    ///
    /// The flanks widen until there are enough points, so a group near an edge
    /// degrades to a one-sided estimate and then to the old whole-trace one
    /// rather than failing.
    double localBackground(const std::vector<double>& trace,
                           std::size_t lo, std::size_t hi,
                           std::size_t flank = 20)
    {
      const std::size_t n = trace.size();
      if (n == 0) { return 0.0; }
      // `span *= 2` never terminates from 0, and `hi` past the end silently
      // yields an empty right flank. Neither is reachable from today's callers;
      // both are one line to make unreachable from any caller.
      if (flank == 0) { flank = 1; }
      if (lo > hi || hi >= n) { hi = n - 1; lo = std::min(lo, hi); }
      std::vector<double> near;
      near.reserve(2 * flank + 8);
      for (std::size_t span = flank; span <= n; span *= 2)
      {
        near.clear();
        const std::size_t l0 = lo > span ? lo - span : 0;
        for (std::size_t i = l0; i < lo && i < n; ++i) { near.push_back(trace[i]); }
        for (std::size_t i = hi + 1; i < n && i <= hi + span; ++i)
        { near.push_back(trace[i]); }
        if (near.size() >= 8) { break; }
      }
      if (near.size() < 3)
      {
        // CLEAR first. Without this the surviving flank points are counted
        // twice, so the fallback was neither the flank estimate nor the
        // whole-trace one it claimed to degrade to.
        near.clear();
        for (std::size_t i = 0; i < n; ++i)
        { if (i < lo || i > hi) { near.push_back(trace[i]); } }
      }
      if (near.empty()) { return 0.0; }
      std::sort(near.begin(), near.end());
      return near[near.size() / 2];
    }

    /// Defined below with the other correlation helpers; declared here because
    /// the selection-time library check belongs beside the boundary code it is
    /// used from, not 300 lines away next to its arithmetic.
    double pearson(const std::vector<double>& a, const std::vector<double>& b);

    /// Do the fragments at `k` have the RELATIVE INTENSITIES the library
    /// predicts? Pearson of background-subtracted per-fragment area against
    /// library intensity, over the same window the sub-scores use.
    ///
    /// Separate from the sub-score of the same name because it is computed at
    /// candidate-SELECTION time, before boundaries exist, for a position that
    /// may never become a candidate. -1 when it cannot be formed.
    double libraryCorrelationAt(const std::vector<std::vector<double>>& tr,
                                const float* lib_intensity, std::size_t k,
                                std::size_t half, std::size_t n)
    {
      if (!lib_intensity || tr.empty() || n == 0) { return -1.0; }
      const std::size_t lo = k > half ? k - half : 0;
      const std::size_t hi = std::min(k + half, n - 1);
      std::vector<double> area, lib;
      area.reserve(tr.size()); lib.reserve(tr.size());
      for (std::size_t f = 0; f < tr.size(); ++f)
      {
        const double w = lib_intensity[f];
        if (!(w > 0.0)) { continue; }
        std::vector<double> flank;
        const std::size_t fl = 20;
        const std::size_t l0 = lo > fl ? lo - fl : 0;
        for (std::size_t j = l0; j < lo; ++j) { flank.push_back(tr[f][j]); }
        for (std::size_t j = hi + 1; j < n && j <= hi + fl; ++j) { flank.push_back(tr[f][j]); }
        double base = 0.0;
        if (!flank.empty())
        {
          std::sort(flank.begin(), flank.end());
          base = flank[flank.size() / 2];
        }
        double a = 0.0;
        for (std::size_t j = lo; j <= hi; ++j) { a += std::max(0.0, tr[f][j] - base); }
        area.push_back(a); lib.push_back(w);
      }
      if (area.size() < 3) { return -1.0; }
      const double r = pearson(area, lib);
      return std::isfinite(r) ? r : -1.0;
    }

    /// Walk out from an apex to a peak's boundaries on a SMOOTHED trace.
    ///
    /// Three guards, and the co-elution picker had none of them. A bare
    /// "descend while above a fraction of the apex" never terminates on a noisy
    /// baseline: measured on the corpus, the median candidate spanned 129 of
    /// 130 cycles and 77.4% covered more than 80% of the extraction window,
    /// against a peak whose measured FWHM is 3.5 s -- about 2.5 cycles. The
    /// "peak group" was the whole window for most precursors, so every
    /// sub-score integrated over [left,right] was integrating the interference
    /// as well as the peak, and `localBackground` -- which estimates from the
    /// points OUTSIDE the group -- had almost nothing left to estimate from.
    ///
    ///  * floor: stop below `boundary_fraction` of the apex.
    ///  * rebound: stop if the trace climbs more than a quarter of the apex
    ///    back above its running minimum, which is a valley into a neighbour.
    ///  * span: never walk further than `max_half` cycles either way -- 20 by
    ///    default, so about 15% of a 130-cycle window per side, not the
    ///    "quarter of the trace" an earlier version of this comment claimed.
    ///    It is a fuse, not a peak-width constraint: it fires on 0.2% of
    ///    confident positives, so it is not what makes the rule work.
    ///
    /// This is the amplitude picker's rule, which had it all along; the two now
    /// share one implementation rather than one having a good rule and the
    /// DEFAULT picker a naive one.
    ///
    /// `apex_out`, when given, receives the SNAPPED apex. It is not optional
    /// bookkeeping: the snap moved the position the floor and the rebound limit
    /// are computed from, and a caller that keeps its original seed ends up
    /// with boundaries derived from one position and an apex reported at
    /// another. Downstream that desynchronises the at-apex fragment count, the
    /// retention time the feature reports, and the one-candidate-per-apex
    /// deduplication -- two seeds two cycles apart snap to the same maximum and
    /// both survive as separate candidates.
    std::pair<std::size_t, std::size_t> peakBounds(const std::vector<double>& sm,
                                                   std::size_t left_from,
                                                   std::size_t right_from,
                                                   double boundary_fraction,
                                                   std::size_t min_cycles,
                                                   std::size_t max_half,
                                                   double boundary_sigmas,
                                                   std::size_t* apex_out = nullptr)
    {
      const std::size_t n = sm.size();
      if (apex_out) { *apex_out = left_from; }
      if (n == 0) { return {left_from, right_from}; }
      left_from = std::min(left_from, n - 1);
      right_from = std::min(right_from, n - 1);
      // Snap to the local maximum of THIS trace near the seed.
      //
      // The co-elution picker chooses its apex from one reference fragment and
      // then walks the SUMMED trace, on which that position need not be a
      // maximum at all. Starting a descent part-way up a slope makes the whole
      // rule meaningless -- the apex value that sets the floor and the rebound
      // limit is not the peak's -- and it is how a weak correct peak acquires a
      // boundary containing its stronger neighbour.
      {
        const std::size_t lo = left_from > 2 ? left_from - 2 : 0;
        const std::size_t hi = std::min(right_from + 2, n - 1);
        std::size_t best = left_from;
        for (std::size_t k = lo; k <= hi; ++k)
        { if (sm[k] > sm[best]) { best = k; } }
        if (best < left_from) { left_from = best; }
        if (best > right_from) { right_from = best; }
        if (apex_out) { *apex_out = best; }
      }
      const double apex = std::max(sm[left_from], sm[right_from]);

      // The floor is NOISE-relative, not apex-relative.
      //
      // Measured on 3,427 confident positives: the local baseline is a median
      // of 36% of the apex, and for 87.6% of them it EXCEEDS a tenth of the
      // apex. An apex-relative 10% floor is therefore unreachable for most real
      // peptides -- which is why the walk ran to the window edge and the median
      // peak group covered 129 of 130 cycles. The fraction was never a
      // threshold; it was a formality.
      //
      // baseline + k*sigma is where the peak has returned to noise, which is
      // what a boundary is. Robust statistics from OUTSIDE a core region so the
      // peak cannot raise its own floor, and never above the apex -- that would
      // give a zero-width group for a real but shallow peak, and shallow is
      // exactly the regime this project is short in.
      // The FALLBACK floor. The normal path is noise-relative and replaces this
      // via the max() below; this term survives for the case that defeats the
      // robust statistics -- a zero-inflated trace where more than half the far
      // points are exactly 0, giving base = 0 and sigma = 0, so `base + k*sigma`
      // is 0 and no k can rescue it. Intensities are non-negative, so a floor of
      // 0 is never crossed and the walk would run to the span bound: the
      // original whole-window failure, returning for precisely the faintest
      // precursors. A reachable fraction of the apex is the right degradation.
      double floor_value = boundary_fraction * apex;
      {
        const std::size_t core = 6;
        std::vector<double> far;
        far.reserve(n);
        for (std::size_t k = 0; k < n; ++k)
        {
          const bool near = (k + core >= left_from) && (k <= right_from + core);
          if (!near) { far.push_back(sm[k]); }
        }
        if (far.size() >= 8)
        {
          std::sort(far.begin(), far.end());
          const double base = far[far.size() / 2];
          std::vector<double> dev(far.size());
          for (std::size_t k = 0; k < far.size(); ++k)
          { dev[k] = std::fabs(far[k] - base); }
          std::sort(dev.begin(), dev.end());
          const double sigma = 1.4826 * dev[dev.size() / 2];
          // CAPPED, not discarded. Letting the noise floor through only when it
          // sits below the apex made the rule discontinuous exactly where it
          // matters: a candidate whose apex is at or below baseline+k*sigma
          // reverted to the fractional floor, which is the floor measured to be
          // unreachable for 87.6% of real peptides. The weakest peaks -- p10 of
          // apex height is 0.9 sigma -- got the original failure mode back, and
          // two near-identical low-S/N peaks could receive radically different
          // boundaries depending on which side of the apex the floor fell.
          //
          // Capping just below the apex is continuous in the signal-to-noise
          // ratio: a sub-sigma candidate gets the narrowest boundary the rule
          // can express, then `min_cycles` widens it to a scoreable width.
          // Whether such a candidate should exist at all is a question for the
          // admission gate, not for a boundary rule.
          // nextafter, not a 1e-9 absolute-scale epsilon: intensities here span
          // several orders of magnitude, and `1e-9 * max(1, |apex|)` is either
          // far larger than one representable step for a bright peak or, for a
          // normalised trace, large enough to matter on its own.
          const double noise_floor =
            std::min(base + boundary_sigmas * sigma,
                     std::nextafter(apex, -std::numeric_limits<double>::infinity()));
          floor_value = std::max(floor_value, noise_floor);
        }
      }
      // The rebound limit is apex-absolute, and DIA-NN's is not: it breaks at a
      // valley only when the valley is DEEP (below apex/3) and the signal has
      // since doubled off it, which is valley-relative and therefore fires on
      // the shallow saddle between co-eluting isomers that an absolute
      // threshold walks straight through. That argument is sound and the rule
      // still lost: swept over four minimum widths and three bands of
      // wrong-apex distance, the valley-relative guard scored lower on
      // separation in all 24 cells (far band, min 5: 0.769 against 0.778;
      // close band: 0.569 against 0.587). Kept as measured, not as preferred.
      const double rebound_limit = 0.25 * apex;
      const std::size_t max_span = std::max<std::size_t>(2, max_half);

      std::size_t l = std::min(left_from, n - 1), guard = 0;
      double run_min = sm[l];
      while (l > 0 && guard++ < max_span)
      {
        const double v = sm[l - 1];
        if (v <= floor_value) { break; }
        if (v > run_min + rebound_limit) { break; }
        run_min = std::min(run_min, v);
        --l;
      }
      std::size_t r = std::min(right_from, n - 1);
      guard = 0;
      run_min = sm[r];
      while (r + 1 < n && guard++ < max_span)
      {
        const double v = sm[r + 1];
        if (v <= floor_value) { break; }
        if (v > run_min + rebound_limit) { break; }
        run_min = std::min(run_min, v);
        ++r;
      }
      // Widen symmetrically to the minimum width. Below it, MS1_COELUTION (5
      // cycles), the mass and mobility blocks (hi > lo) and RT_SPREAD all
      // either vanish or report agreement they did not measure.
      while (r - l + 1 < min_cycles && (l > 0 || r + 1 < n))
      {
        const bool can_l = l > 0, can_r = r + 1 < n;
        if (can_l && (!can_r || sm[l - 1] >= sm[r + 1])) { --l; }
        else if (can_r) { ++r; }
        else { break; }
      }
      return {l, r};
    }

    std::vector<double> smooth(const std::vector<double>& x, std::size_t half)
    {
      if (half == 0 || x.size() < 2 * half + 1) { return x; }
      std::vector<double> out(x.size(), 0.0);
      for (std::size_t i = 0; i < x.size(); ++i)
      {
        const std::size_t lo = i >= half ? i - half : 0;
        const std::size_t hi = std::min(x.size() - 1, i + half);
        double s = 0.0;
        for (std::size_t k = lo; k <= hi; ++k) { s += x[k]; }
        out[i] = s / static_cast<double>(hi - lo + 1);
      }
      return out;
    }

    /// Why the co-elution detector rejected a scan position.
    ///
    /// Six criteria reject silently, so "19,150 precursors yielded no candidate"
    /// said nothing about WHICH test was responsible. On a realistic library
    /// 14% of DIA-NN's confident precursors get no candidate at all, and those
    /// lost are heavier (median m/z 894 vs 661), higher-mobility and ~34% less
    /// abundant -- a pattern that points at a threshold rather than at absence
    /// of signal, but only counting can say which threshold.
    struct PickerRejects
    {
      std::size_t too_few_present[2] = {0, 0};   ///< [0] target, [1] decoy   ///< <2 fragments in {k-1,k,k+1}
      /// Precursors that never entered the correlation loop at all.
      ///
      /// `n < 2*S+4 || tc < 2` returns before anything is counted, so these were
      /// invisible: they landed in `precursors_without_candidate` with no way to
      /// tell them from a precursor the correlation test rejected. On Astral
      /// 5,187 precursors yield no candidate while OpenSWATH features 100% of
      /// the library, and the census could not say which stage lost them.
      std::size_t too_few_cycles = 0;
      /// Never got a usable chromatogram OUT OF THE EXTRACTOR at all. These are
      /// extraction losses, not picking losses, and they were pooled with
      /// picker rejects under `precursors_without_candidate`.
      std::size_t no_points[2] = {0, 0};   ///< [0] target, [1] decoy        ///< pointCount(0) < 3
      /// THREE different gates shared this counter and only the last is what
      /// the name says. Every log before 2026-08-18 therefore reported Gate C
      /// rejections as "empty trace" -- doc/34's "84.5% of targets had an
      /// all-zero trace" was in fact 84.5% rejected by Gate C, a completely
      /// different statement, and Gate C is ON by default (gate_alpha 0.05).
      /// Kept as the sum so older numbers stay comparable; read the three below.
      std::size_t empty_trace[2] = {0, 0};   ///< [0] target, [1] decoy
      /// Co-elution evidence below the (1-alpha) quantile of the decoy null.
      std::size_t gate_c[2] = {0, 0};
      /// Too few transitions showing a noise excursion (the -gate_alpha 0 path).
      std::size_t few_excursions[2] = {0, 0};
      /// The summed trace really is zero -- nothing extracted at all.
      std::size_t zero_trace[2] = {0, 0};
      std::size_t too_few_transitions[2] = {0, 0};   ///< [0] target, [1] decoy
      std::size_t masked_candidates = 0;   ///< v1.15: candidates scored with a transition mask
      std::size_t masked_fragments = 0;    ///< v1.15: transitions removed over those candidates
      /// Entered the loop, computed correlations, and found no qualifying
      /// position anywhere in the window.
      std::size_t no_hit_anywhere = 0;
      std::size_t below_corr[2] = {0, 0};   ///< [0] target, [1] decoy        ///< reference corr sum < min_corr_score
      std::size_t reference_zero[2] = {0, 0};   ///< [0] target, [1] decoy    ///< smoothed reference not positive
      std::size_t not_local_max[2] = {0, 0};   ///< [0] target, [1] decoy     ///< k is not the local maximum
      std::size_t below_apex_evidence[2] = {0, 0};   ///< [0] target, [1] decoy
      std::size_t outside_margin[2] = {0, 0};   ///< [0] target, [1] decoy
      /// Precursors that reached Session::add, per class. The code says these
      /// MUST be equal -- assignment is unconditional, the only precursor-level
      /// drop is covering==0 which depends solely on m/z, and every decoy shares
      /// its target's m/z (verified: 4,991,888 of 4,991,888). The scan counters
      /// say otherwise (1.52x). One of those is wrong and this is the counter
      /// that decides which.
      std::size_t reached[2] = {0, 0};    ///< beyond MaxCorrDiff of the best
      std::size_t too_few_at_apex = 0;   ///< candidate emitted, then dropped by the scorer
      std::size_t scans[2] = {0, 0};   ///< [0] target, [1] decoy             ///< positions examined
    };

    struct Candidate
    {
      std::size_t apex = 0, left = 0, right = 0;
      double apex_value = 0.0;
      /// Summed pairwise fragment correlation at this position, from the
      /// co-elution detector. 0 from the amplitude detector, which never
      /// computes it.
      double corr_sum = 0.0;
    };

    /// Local maxima of the smoothed trace, strongest first, with boundaries
    /// descending to a fraction of the apex.
    std::vector<Candidate> findCandidates(const std::vector<double>& smoothed,
                                          std::size_t max_candidates,
                                          double boundary_fraction,
                                          std::size_t min_cycles,
                                          std::size_t max_half,
                                          double boundary_sigmas)
    {
      std::vector<Candidate> found;
      const std::size_t n = smoothed.size();
      if (n < 3) { return found; }

      for (std::size_t i = 1; i + 1 < n; ++i)
      {
        if (smoothed[i] <= 0.0) { continue; }
        if (smoothed[i] < smoothed[i - 1] || smoothed[i] < smoothed[i + 1]) { continue; }

        // A plateau is a peak, not a thing to skip. The old `== previous`
        // guard rejected it outright, so a broad real peak -- exactly what a
        // moving average produces from a strong one -- yielded no candidate at
        // all while a narrow noise spike passed. Take the plateau's extent and
        // report its midpoint, and emit it once by skipping to its end.
        std::size_t plateau_end = i;
        while (plateau_end + 1 < n && smoothed[plateau_end + 1] == smoothed[i]) { ++plateau_end; }
        if (plateau_end + 1 < n && smoothed[plateau_end + 1] > smoothed[i]) { continue; }

        Candidate c;
        c.apex = (i + plateau_end) / 2;
        c.apex_value = smoothed[i];
        // The same scale-aware descent the co-elution picker now uses. It lived
        // here, inline, and only here -- which is how the DEFAULT picker came to
        // run without it.
        const auto b = peakBounds(smoothed, i, plateau_end, boundary_fraction,
                                  min_cycles, max_half, boundary_sigmas);
        i = plateau_end;
        c.left = b.first;
        c.right = b.second;
        found.push_back(c);
      }

      std::sort(found.begin(), found.end(),
                [](const Candidate& a, const Candidate& b) { return a.apex_value > b.apex_value; });
      if (found.size() > max_candidates) { found.resize(max_candidates); }
      return found;
    }

    double pearson(const std::vector<double>& a, const std::vector<double>& b)
    {
      const std::size_t n = std::min(a.size(), b.size());
      if (n < 2) { return 0.0; }
      double ma = 0, mb = 0;
      for (std::size_t i = 0; i < n; ++i) { ma += a[i]; mb += b[i]; }
      ma /= n; mb /= n;
      double num = 0, da = 0, db = 0;
      for (std::size_t i = 0; i < n; ++i)
      {
        const double x = a[i] - ma, y = b[i] - mb;
        num += x * y; da += x * x; db += y * y;
      }
      const double d = std::sqrt(da) * std::sqrt(db);
      return d > 0.0 ? num / d : 0.0;
    }

    /// Library correlation, refusing to report one it cannot estimate.
    ///
    /// `pearson` guards only n < 2, so at n = 2 it returns +-1 ALWAYS -- two
    /// points are collinear by construction. The IH1 library carries 366,084
    /// targets with 0-2 fragments, and MEASURED on the full run the score was
    /// monotonically INVERTED in the evidence behind it:
    ///
    ///     usable fragments      2       3       4      12
    ///     var_library_corr   1.0000  1.0000  1.0000  0.6667
    ///     frac exactly 1.0    79.7%   63.8%   51.8%    5.6%
    ///     median DScore      14.389   1.760   0.734  -0.273
    ///
    /// The first attempt at this shrank r^2 by its null expectation 1/(n-1).
    /// That is the textbook correction and it was WRONG HERE, because it
    /// assumes the feature's working range sits well above the null. It does
    /// not: on the baseline run the median library correlation at twelve
    /// fragments is 0.0137, far below the 1/sqrt(n-1) = 0.302 cutoff, so the
    /// shrinkage zeroed 66.7% of the twelve-fragment population -- the very
    /// population that carries the identifications -- while claiming not to.
    ///
    /// What is actually degenerate is the SAMPLE SIZE, not the correlation.
    /// `corrected` clamps a fragment below its own local background to exactly
    /// 0, so the informative points are the non-zero ones; a precursor with two
    /// of those cannot support a correlation however many transitions it
    /// nominally has. Below four, report no evidence rather than arithmetic.
    /// Above it, the correlation is passed through untouched.
    double libraryCorrelation(const std::vector<double>& observed,
                              const std::vector<double>& library_intensity)
    {
      const std::size_t n = std::min(observed.size(), library_intensity.size());
      std::size_t informative = 0;
      for (std::size_t i = 0; i < n; ++i)
      {
        if (observed[i] > 0.0) { ++informative; }
      }
      if (informative < 4) { return 0.0; }
      return pearson(observed, library_intensity);
    }

    /// Candidates detected by CO-ELUTION rather than by amplitude.
    ///
    /// This is DIA-NN's `Searcher::peaks` (diann.cpp:7423) as described in
    /// DIA-NN-workflow-handoff.md §5, and it differs from `findCandidates` in
    /// the one way that matters: the acceptance criterion is the pairwise
    /// correlation among the precursor's own fragments, evaluated AT detection
    /// time, not the height of a summed trace with correlation checked later.
    ///
    /// Why that ordering is the whole point. A standardised sum is dominated by
    /// whatever is loud: one bright interfering ion, or a dozen unrelated
    /// species that happen to share the window, out-sum a real but weak
    /// peptide. So an amplitude-detected candidate list is ordered by "how much
    /// signal is here" and the correct peak lands at rank 2, 3, 7. Measured on
    /// IH1 against DIA-NN's confident set, our amplitude picker put the true
    /// peak first for 24.4% of precursors, inside the top 3 for 47.2%, and
    /// inside the top 25 for 89.8% -- a slow decay, which is the signature of a
    /// detector that finds the peak but cannot tell it from its neighbours.
    ///
    /// Fragments of one peptide come from one eluting molecule, so they
    /// correlate. Interference does not.
    ///
    /// Three parameters are DIA-NN's, named as it names them:
    ///  * `min_corr_score` (0.5) -- the reference fragment's summed correlation
    ///    to the others must reach this, or the position is not a peak at all;
    ///  * `max_corr_diff` (2.0) -- candidates are kept by MARGIN from the best
    ///    correlation sum, so a precursor with one convincing peak yields one
    ///    candidate and an ambiguous one yields several. That is why a fixed
    ///    top-N hurts: at depth 25 it manufactures 24 competitors whether or
    ///    not any of them is plausible.
    ///  * `apex_evidence` (0.99) -- the apex must be essentially the local
    ///    maximum of the reference fragment's own smoothed trace, which is a
    ///    far stricter shape test than a maximum of a sum.
    std::vector<Candidate> findCandidatesByCorrelation(
      const PrecursorChromatogram& c, PickerRejects& rej, bool is_decoy,
      std::size_t half_window,
      double min_corr_score, double max_corr_diff, double apex_evidence,
      std::size_t smooth_half_width, double boundary_fraction,
      std::size_t max_candidates, std::size_t min_cycles, std::size_t max_half,
      std::size_t boundary_smooth_half, double boundary_sigmas,
      const float* lib_intensity, std::size_t score_half,
      std::size_t min_separation, double lib_weight)
    {
      std::vector<Candidate> found;
      const std::uint32_t tc = c.transition_count;
      const std::size_t n = c.cycles;
      const std::size_t S = std::max<std::size_t>(1, half_window);
      if (tc < 2) { ++rej.too_few_transitions[is_decoy]; return found; }
      if (n < 2 * S + 4) { ++rej.too_few_cycles; return found; }

      // Traces once, smoothed once. DIA-NN smooths the reference trace before
      // the local-maximum test (kernel 1/4-1/2-1/4); `smooth` here is the
      // moving average the rest of this file uses, which is the same idea.
      std::vector<std::vector<double>> tr(tc), sm(tc);
      for (std::uint32_t k = 0; k < tc; ++k)
      {
        const std::uint32_t m = c.pointCount(k);
        const float* pk = m ? c.trace(k) : nullptr;
        tr[k].assign(n, 0.0);
        for (std::uint32_t j = 0; j < m && j < n; ++j) { tr[k][j] = pk[j]; }
        sm[k] = smooth(tr[k], smooth_half_width);
      }

      struct Hit { std::size_t apex; double corr_sum; };
      std::vector<Hit> hits;
      std::vector<double> a, b;
      // The bound is what the WINDOW requires, stated directly. It used to be
      // `k = S + 1; k + S + 2 < n`, which is one position tighter at the start
      // and two at the end than the correlation window actually needs: the
      // window is [k-S, k+S+1) so it fits for every k from S to n-S-1. Three
      // scan positions per precursor were therefore never examined, and they
      // did not land in any reject counter either -- they were not rejected,
      // they were never considered, which is the kind of loss no amount of
      // staring at the rejection table would have found.
      //
      // Small but not nothing: 198 of 33,337 confident positives, 0.59%, have
      // their true apex on exactly those three cycles. Worth stating in a
      // project whose emission gap is the thing being chased.
      for (std::size_t k = S; k + S + 1 <= n; ++k)
      {
        // Cheap rejects first: something must be here, and it must persist
        // across neighbouring cycles rather than being a single spike.
        std::size_t present = 0;
        for (std::uint32_t f = 0; f < tc; ++f)
        {
          if (tr[f][k - 1] > 0.0 || tr[f][k] > 0.0 || tr[f][k + 1] > 0.0) { ++present; }
        }
        ++rej.scans[is_decoy];
        if (present < 2) { ++rej.too_few_present[is_decoy]; continue; }

        // Pairwise correlation over [k-S, k+S]; each fragment scores the sum of
        // its correlations to the others.
        const std::size_t lo = k - S, hi = k + S + 1;
        std::vector<double> score(tc, 0.0);
        for (std::uint32_t i = 0; i < tc; ++i)
        {
          a.assign(tr[i].begin() + lo, tr[i].begin() + hi);
          for (std::uint32_t j = i + 1; j < tc; ++j)
          {
            b.assign(tr[j].begin() + lo, tr[j].begin() + hi);
            const double r = pearson(a, b);
            if (std::isfinite(r)) { score[i] += r; score[j] += r; }
          }
        }

        std::vector<std::uint32_t> order(tc);
        for (std::uint32_t i = 0; i < tc; ++i) { order[i] = i; }
        std::stable_sort(order.begin(), order.end(),
                         [&](std::uint32_t x, std::uint32_t y) { return score[x] > score[y]; });

        // Walk the ranked fragments; the first that satisfies every test makes
        // this position a peak, and we stop -- at most one candidate per cycle.
        for (const std::uint32_t ref : order)
        {
          if (score[ref] < min_corr_score) { ++rej.below_corr[is_decoy]; break; }
          if (!(sm[ref][k] > 0.0)) { ++rej.reference_zero[is_decoy]; continue; }

          const std::size_t half = std::max<std::size_t>(S / 3, 1);
          bool is_max = true;
          for (std::size_t j = (k > half ? k - half : 0); j <= k + half && j < n; ++j)
          {
            if (sm[ref][j] > sm[ref][k]) { is_max = false; break; }
          }
          if (!is_max) { ++rej.not_local_max[is_decoy]; continue; }

          double best_near = 0.0;
          const std::size_t e = S > 1 ? S - 1 : 1;
          for (std::size_t j = (k > e ? k - e : 0); j <= k + e && j < n; ++j)
          {
            best_near = std::max(best_near, sm[ref][j]);
          }
          if (best_near > 0.0 && sm[ref][k] < apex_evidence * best_near)
          { ++rej.below_apex_evidence[is_decoy]; continue; }

          // `score[ref]`, not `score[argmax]`: this position's evidence is
          // credited to the fragment that actually WITNESSED an apex here, and
          // the argmax fragment is, by construction, one that just failed to.
          //
          // The choice is consequential and was never measured until now. The
          // stored value is not the position's best correlation for 67.7% of
          // hit positions, deflated by a median of 0.709 against a margin of
          // 2.0, and ranking by the argmax instead would change the top-ranked
          // position for 69.3% of precursors and the margin-surviving set for
          // 85.3%.
          //
          // It costs nothing. Target fraction among the top-N with decoys as
          // control is 98.6% / 92.3% / 86.2% at N = 1,000 / 2,000 by the stored
          // value against 98.2% / 91.6% / 84.4% by the argmax. Two orderings
          // that disagree about which position to name on two thirds of
          // precursors separate targets from decoys equally well -- which says
          // the positions inside a margin set are largely interchangeable, and
          // is why every change to this ordering has had a small effect.
          hits.push_back({k, score[ref]});
          break;
        }
      }
      if (hits.empty()) { ++rej.no_hit_anywhere; return found; }

      // Keep by MARGIN from the best, not by rank.
      double best = 0.0;
      for (const auto& h : hits) { best = std::max(best, h.corr_sum); }
      std::stable_sort(hits.begin(), hits.end(),
                       [](const Hit& x, const Hit& y) { return x.corr_sum > y.corr_sum; });

      // ORDER the survivors by library agreement as well as by co-elution.
      //
      // corr_sum is a statement about SHAPE: do these fragments rise and fall
      // together. It is the same kind of evidence at every position, and it is
      // the only kind this detector has ever had -- which is why every attempt
      // to fix candidate emission by moving its thresholds failed. Measured:
      // raising the cap 3 -> 20 buys 15.9 points of recall and loses 3.4 of
      // selection accuracy; relaxing apex_evidence loses on both axes; ranking
      // by the co-elution sub-score instead of corr_sum is worse than corr_sum.
      // Shape had been exhausted.
      //
      // The library says something corr_sum cannot: which fragments should be
      // BRIGHT. Two co-eluting species have equally good shape and different
      // relative intensities, so this separates exactly the case shape cannot.
      // Measured over 11,731 confident positives, ranking the margin survivors
      // by normalised corr_sum + library correlation rather than corr_sum:
      //
      //     statistic            recall@3   selection accuracy
      //     corr_sum               69.5%          59.9%
      //     library correlation    71.3%          63.6%
      //     both                   72.7%          63.3%
      //
      // The first change to this detector that improves BOTH -- on that ruler.
      //
      // OFF BY DEFAULT, because that ruler is agreement with DIA-NN, and a
      // reference-free measurement does not confirm it. Target fraction among
      // the top-N with decoys as the control, no external tool involved:
      //
      //     discriminant     corr_sum   corr+lib      (top-500)
      //     lib                 96.8%      91.0%
      //     coelution           93.8%      94.8%
      //
      // The pattern is the winner's curse, not a defect in the ordering: an arm
      // that selects the position maximising library correlation, and is then
      // ranked by something containing library correlation, lets DECOYS shop
      // for their best value too. The more the discriminant overlaps the
      // selection statistic the worse it looks; on a discriminant it does not
      // touch, it is mildly ahead.
      //
      // That matters here specifically because LIBRARY_CORR is one of the
      // nineteen features the classifier sees, so the shipped pipeline is the
      // coupled case rather than the independent one -- diluted one-in-nineteen,
      // but in the direction the table warns about. Recall rising on a
      // DIA-NN-agreement metric is, by this project's own rule, never
      // sufficient evidence: it has already risen 85.1% to 93.8% once while
      // identifications fell by 319.
      //
      // So this ships as a knob at 0. The measurement that can settle it is a
      // full run gated on entrapment FDP, and this exists so that run is a flag
      // rather than a patch.
      //
      // Applied AFTER the margin filter and BEFORE the cap, which is where it
      // was measured. The gates and the margin still run on corr_sum alone, so
      // nothing new is admitted -- only the order in which the cap keeps what
      // corr_sum already accepted.
      if (lib_weight > 0.0 && lib_intensity && hits.size() > 1)
      {
        std::vector<Hit> kept;
        kept.reserve(hits.size());
        for (const auto& h : hits)
        {
          if (h.corr_sum >= best - max_corr_diff) { kept.push_back(h); }
        }
        if (kept.size() > 1)
        {
          double lo_c = kept.front().corr_sum, hi_c = kept.front().corr_sum;
          for (const auto& h : kept)
          { lo_c = std::min(lo_c, h.corr_sum); hi_c = std::max(hi_c, h.corr_sum); }
          const double span = hi_c - lo_c;
          std::vector<std::pair<double, Hit>> keyed;
          keyed.reserve(kept.size());
          for (const auto& h : kept)
          {
            // -1 when the correlation cannot be formed, which is the worst a
            // correlation can be: a position that cannot be checked against the
            // library must not outrank one that was checked and agreed.
            const double lc = libraryCorrelationAt(tr, lib_intensity, h.apex,
                                                   score_half, n);
            const double norm = span > 0.0 ? (h.corr_sum - lo_c) / span : 0.0;
            keyed.emplace_back(norm + lib_weight * lc, h);
          }
          std::stable_sort(keyed.begin(), keyed.end(),
                           [](const auto& x, const auto& y) { return x.first > y.first; });
          hits.clear();
          for (auto& k : keyed) { hits.push_back(k.second); }
        }
      }

      // Boundaries from the summed trace, as before: the extent of a peak is
      // not what changed here, only which positions are peaks.
      std::vector<double> total(n, 0.0);
      for (std::uint32_t f = 0; f < tc; ++f)
      {
        for (std::size_t j = 0; j < n; ++j) { total[j] += tr[f][j]; }
      }
      // Boundaries are found on a smoothed copy, as the amplitude picker does:
      // on the raw sum a single noise point terminates or extends the walk.
      // apex_value and the traces themselves stay raw.
      const std::vector<double> total_smoothed = smooth(total, boundary_smooth_half);
      for (const auto& h : hits)
      {
        if (h.corr_sum < best - max_corr_diff) { ++rej.outside_margin[is_decoy]; break; }
        if (found.size() >= max_candidates) { break; }
        Candidate cd;
        cd.apex = h.apex;
        cd.apex_value = total[h.apex];
        cd.corr_sum = h.corr_sum;
        // Boundaries on the SMOOTHED sum and with the amplitude picker's
        // guards. The previous rule -- descend while above 10% of the apex,
        // with no rebound guard and no span bound -- gave a median candidate
        // width of 129 of 130 cycles on a peak whose FWHM is 2.5 cycles.
        // The snapped apex, not the seed. This picker chooses its apex from ONE
        // reference fragment, and on the summed trace that position need not be
        // a maximum -- which is exactly the case the snap exists for, and the
        // reason the desynchronisation matters here and not in the amplitude
        // picker, whose seed is a local maximum by construction.
        std::size_t snapped = h.apex;
        const auto b = peakBounds(total_smoothed, h.apex, h.apex, boundary_fraction,
                                  min_cycles, max_half, boundary_sigmas, &snapped);
        cd.left = b.first; cd.right = b.second;
        cd.apex = snapped;
        cd.apex_value = total[snapped];
        // One candidate per BASIN, not per scan position.
        //
        // At a separation of 1 this is the old rule -- one candidate per apex,
        // as DIA-NN does -- and that is the default, because widening it is
        // measured to trade one thing for another rather than to be free.
        //
        // What widening buys: a cap of three slots currently means three
        // positions, which can be three samples of ONE broad basin, so the true
        // peak is evicted by near-duplicates of a wrong one. That eviction is
        // not hypothetical: relaxing `apex_evidence` makes recall FALL, and
        // relaxing an admission gate cannot shrink a pre-cap superset, so the
        // extra positions must be consuming slots. Suppressing within a basin
        // raises recall@3 from 72.7% to 75.5%.
        //
        // What it costs: the freed slots go to genuinely DIFFERENT peaks, which
        // are genuinely able to outscore the true one. Selection accuracy, on a
        // single sub-score as selector, falls 63.3% to 57.7%.
        //
        // Off by default for that reason. The measurement that could justify it
        // needs the full nineteen-feature model and target-decoy competition,
        // not one feature -- the model would have to be 7.4 points better
        // conditionally than the proxy selector, which is a real question and
        // not one a proxy can answer. The option exists so that experiment can
        // be run without a code change.
        bool dup = false;
        const std::size_t sep = std::max<std::size_t>(1, min_separation);
        for (const auto& g : found)
        {
          const std::size_t d = g.apex > cd.apex ? g.apex - cd.apex : cd.apex - g.apex;
          if (d < sep) { dup = true; break; }
        }
        if (!dup) { found.push_back(cd); }
      }
      return found;
    }

    /// The best fragment's summed correlation to the others at position `k`,
    /// with NO gate applied. This is the same quantity
    /// `findCandidatesByCorrelation` thresholds on, computed for a position
    /// somebody else chose.
    double bestCorrSumAt(const std::vector<std::vector<double>>& tr, std::size_t k,
                         std::size_t S, std::size_t n)
    {
      const std::size_t tc = tr.size();
      if (tc < 2 || k < S || k + S + 1 > n) { return 0.0; }
      const std::size_t lo = k - S, hi = k + S + 1;
      std::vector<double> score(tc, 0.0), a, b;
      for (std::size_t i = 0; i < tc; ++i)
      {
        a.assign(tr[i].begin() + lo, tr[i].begin() + hi);
        for (std::size_t j = i + 1; j < tc; ++j)
        {
          b.assign(tr[j].begin() + lo, tr[j].begin() + hi);
          const double r = pearson(a, b);
          if (std::isfinite(r)) { score[i] += r; score[j] += r; }
        }
      }
      return *std::max_element(score.begin(), score.end());
    }

    /// Co-elution as a FEATURE rather than a GATE: the hybrid.
    ///
    /// The two pickers disagree about what a candidate is. DIA-NN's asks
    /// whether the fragments rise and fall together and refuses the position if
    /// they do not; OpenSWATH's asks only whether there is a peak in the summed
    /// trace and leaves co-elution to the scores. A gate cannot be undone by a
    /// classifier and a feature can, so the gate is only correct if it is never
    /// wrong -- and it is a fixed threshold on an unnormalised sum, which is
    /// not the kind of thing that is never wrong.
    ///
    /// So take BOTH candidate sets, and give every candidate the co-elution
    /// evidence as a number. An amplitude candidate the gate would have thrown
    /// away now arrives with a low `corr_sum`, and the semi-supervised
    /// classifier and the FDR decide what that is worth. This also repairs the
    /// comparison: run with either alternative picker alone, `var_corr_sum` and
    /// `var_candidate_margin` are constant columns and are dropped, so those
    /// arms were scoring with 15 features against the co-elution arm's 17.
    ///
    /// Ordered by `corr_sum`, so when the cap bites it is the positions with
    /// the least co-elution evidence that are dropped -- the union can only add
    /// candidates the gate refused, never displace ones it accepted.
    std::vector<Candidate> unionCandidates(const PrecursorChromatogram& c,
                                           std::vector<Candidate> found,
                                           const std::vector<Candidate>& amplitude,
                                           std::size_t half_window,
                                           std::size_t max_candidates)
    {
      const std::uint32_t tc = c.transition_count;
      const std::size_t n = c.cycles;
      const std::size_t S = std::max<std::size_t>(1, half_window);
      if (tc < 2 || n < 2 * S + 4) { return found; }

      std::vector<std::vector<double>> tr(tc);
      for (std::uint32_t k = 0; k < tc; ++k)
      {
        const std::uint32_t m = c.pointCount(k);
        const float* pk = m ? c.trace(k) : nullptr;
        tr[k].assign(n, 0.0);
        for (std::uint32_t j = 0; j < m && j < n; ++j) { tr[k][j] = pk[j]; }
      }

      // Two apexes within S/3 cycles are the same peak seen by two detectors --
      // the neighbourhood the co-elution picker's own local-maximum test uses.
      // Admitting both would let one peak occupy two of the capped slots and
      // would give the classifier a duplicate row to compete against itself.
      const std::size_t near = std::max<std::size_t>(S / 3, 1);
      std::vector<Candidate> added;
      for (const auto& a : amplitude)
      {
        bool dup = false;
        for (const auto& g : found)
        {
          const std::size_t d = g.apex > a.apex ? g.apex - a.apex : a.apex - g.apex;
          if (d <= near) { dup = true; break; }
        }
        if (dup) { continue; }
        Candidate cd = a;
        cd.corr_sum = bestCorrSumAt(tr, cd.apex, S, n);
        added.push_back(cd);
      }

      // THE ADDED ONES FILL REMAINING SLOTS; THEY NEVER DISPLACE.
      //
      // The first version sorted the whole union by corr_sum and truncated,
      // which was wrong and measurably so: an amplitude candidate's corr_sum is
      // computed ungated and can exceed that of a co-elution candidate the
      // detector accepted, so once the co-elution set reached the cap the union
      // both ADDED weak candidates and DROPPED accepted ones -- strictly the
      // worst of the two designs. On IH1 that took 1,306 identifications to 818.
      //
      // The fixture missed it because it produces two candidates against a cap
      // of five, so the truncation never ran.
      if (found.size() >= max_candidates) { return found; }
      std::stable_sort(added.begin(), added.end(),
                       [](const Candidate& x, const Candidate& y)
                       { return x.corr_sum > y.corr_sum; });
      const std::size_t room = max_candidates - found.size();
      if (added.size() > room) { added.resize(room); }
      found.insert(found.end(), added.begin(), added.end());
      return found;
    }

    double dotProduct(std::vector<double> a, std::vector<double> b)
    {
      const std::size_t n = std::min(a.size(), b.size());
      double na = 0, nb = 0, dot = 0;
      for (std::size_t i = 0; i < n; ++i)
      {
        na += a[i] * a[i]; nb += b[i] * b[i]; dot += a[i] * b[i];
      }
      const double d = std::sqrt(na) * std::sqrt(nb);
      return d > 0.0 ? dot / d : 0.0;
    }
  } // namespace

  std::pair<std::size_t, std::size_t> PeakGroupScorer::peakBoundsForTest(
    const std::vector<double>& sm, std::size_t l, std::size_t r, double bf,
    std::size_t min_cycles, std::size_t max_half, double sigmas)
  { return peakBounds(sm, l, r, bf, min_cycles, max_half, sigmas); }

  const std::vector<std::string>& PeakGroupScorer::subScoreNames()
  {
    static const std::vector<std::string> names{
      "var_xcorr_shape", "var_xcorr_coelution", "var_library_corr",
      "var_library_dotprod", "var_intensity_score", "var_log_sn",
      "var_usable_fragments", "var_library_rmsd", "var_yseries_score",
      "var_fragment_coverage",
      "var_corr_sum", "var_candidate_margin", "var_peak_width_ratio",
      "var_im_delta", "var_ms1_coelution",
      "var_mass_accuracy", "var_mass_spread", "var_im_spread",
      "var_rt_spread", "var_mass_survival", "var_null_control",
      "var_ref_corr_sum",
      "var_ref_corr_1", "var_ref_corr_2", "var_ref_corr_3", "var_ref_corr_4",
      "var_ref_corr_5", "var_ref_corr_6", "var_ref_corr_7", "var_ref_corr_8",
      "var_ref_corr_9", "var_ref_corr_10", "var_ref_corr_11", "var_ref_corr_12",
      "var_sig_share_1", "var_sig_share_2", "var_sig_share_3",
      "var_sig_share_4", "var_sig_share_5", "var_sig_share_6",
      "var_cand_rank", "var_cand_count"};
    return names;
  }

  const std::vector<std::string>& PeakGroupScorer::fragvecNames()
  {
    // Built rather than spelled out: 72 of the 78 are `<base>_<k>` for k = 1..12
    // and writing them by hand is 72 chances to transpose a digit in a column
    // order that a comparison against the reference builder would then report
    // as an arithmetic disagreement.
    static const std::vector<std::string> names = [] {
      std::vector<std::string> n;
      n.reserve(N_FRAGVEC);
      for (const char* base : {"R1_LOGAREA", "R1_SHARE", "R1_LOGRATIO",
                               "R1_ATAPEX", "R1_MEAS", "R1_ABSENT"})
      {
        for (int k = 1; k <= 12; ++k)
        { n.push_back(std::string(base) + "_" + std::to_string(k)); }
      }
      for (const char* scalar : {"R1_N_MEAS", "R1_N_ABSENT", "R1_N_PRESENT",
                                 "R1_LOGTOT", "R1_LIB_CORR", "R1_LIB_CORR_LOO"})
      { n.emplace_back(scalar); }
      return n;
    }();
    return names;
  }

  namespace
  {
    // thread_local for speed -- these increment once per scan position, so a
    // shared atomic would serialise the hottest loop in the scorer. But a
    // thread_local counter that is REPORTED from one thread reports one
    // thread's slice of the work, and nothing distributes chromatograms to
    // threads in a class-balanced way.
    //
    // That is not hypothetical: reading the unaggregated counters produced a
    // "1.52x more scan positions for decoys" that was taken as evidence the
    // decoy excess arises upstream of the picker. It arose from thread
    // scheduling. The peak-group excess itself is real -- it comes from the
    // scored result, not from here -- but its LOCATION was wrong.
    //
    // So each thread's instance registers itself once and the reporter sums
    // them. The hot path stays a plain increment.
    /// Gate C's threshold, calibrated from the run's OWN decoy null.
    ///
    /// It cannot be known up front -- the null does not exist until decoys have
    /// been seen -- and assuming a distribution is exactly what made Gate B's
    /// "3 sigma" mean 81% instead of 0.1%. So the first `calibration_n` decoys
    /// are admitted unconditionally while their statistics accumulate, then tau
    /// is the (1-alpha) quantile of those and applies from that point on.
    ///
    /// Admitting the first tranche costs nothing at this scale: 20,000 against
    /// a 9,983,789-precursor library is 0.2%, and they are scored normally.
    struct NullCalibrationBody
    {
      std::mutex mu;
      std::vector<double> decoy_stats;
      double tau = 0.0;
      bool ready = false;
      std::FILE* log = nullptr;
      std::mutex log_mu;
      std::size_t admitted_uncalibrated = 0;
      // v1.17: what the calibration sample was made of, for the 'null armed' line.
      std::size_t zeros = 0;                     ///< sample statistics exactly 0
      double rt_first = std::numeric_limits<double>::quiet_NaN();
      double rt_last = std::numeric_limits<double>::quiet_NaN();

      /// Returns true if the precursor should be admitted. @p rt_centre is the
      /// window's position in run seconds; with @p rt_min > 0 a decoy before it
      /// is admitted (as every precursor is while the null is built) but does
      /// not enter the calibration sample (v1.17). rt_min 0 = off.
      bool admit(double stat, bool is_decoy, std::size_t n_needed, double alpha,
                 double rt_centre, double rt_min)
      {
        std::lock_guard<std::mutex> g(mu);
        if (!ready)
        {
          if (is_decoy && !(rt_min > 0.0 && rt_centre < rt_min))
          {
            decoy_stats.push_back(stat);
            if (stat == 0.0) { ++zeros; }
            if (decoy_stats.size() == 1) { rt_first = rt_centre; }
            rt_last = rt_centre;
          }
          if (decoy_stats.size() >= n_needed)
          {
            std::sort(decoy_stats.begin(), decoy_stats.end());
            const std::size_t k = std::min(decoy_stats.size() - 1,
              std::size_t((1.0 - alpha) * double(decoy_stats.size())));
            tau = decoy_stats[k];
            ready = true;
            // v1.17: say so, once per Session, on the channel the picker census
            // uses (stderr -> the run log). Nothing was logged before, which is
            // how tau = 0 from pre-gradient windows went unnoticed.
            std::fprintf(stderr,
                         "gate C null armed: n=%zu decoys, tau=%.6g, zeros=%.4f (%zu), "
                         "first/last calibration RT %.1f-%.1f s, rt_min %.1f s\n",
                         decoy_stats.size(), tau,
                         double(zeros) / double(decoy_stats.size()), zeros,
                         rt_first, rt_last, rt_min);
          }
          ++admitted_uncalibrated;
          return true;                 // admit while the null is being built
        }
        return stat >= tau;
      }

      /// D6: a threshold decided before any precursor arrives. With ready set
      /// the warm-up branch above never runs; nothing is admitted unchecked.
      void freeze(double t)
      {
        std::lock_guard<std::mutex> g(mu);
        tau = t;
        ready = true;
      }

      /// Record one decision. The path is passed in rather than stored because
      /// this struct is a file-scope singleton declared before Options is in
      /// scope here; the file is opened on first use and closed at exit.
      void note(const std::string& path, std::uint32_t precursor, bool is_decoy,
                double stat, bool admitted, bool was_ready)
      {
        if (path.empty()) { return; }
        std::lock_guard<std::mutex> g(log_mu);
        if (log == nullptr)
        {
          // APPEND, not truncate. Each Session owns its own calibration and so
          // opens this file independently; with "w" the second pass wiped the
          // first pass's decisions and the log appeared to contain a single tau.
          // The header is written only into an empty file.
          log = std::fopen(path.c_str(), "a");
          if (log == nullptr) { return; }
          if (std::ftell(log) == 0)
          { std::fprintf(log, "precursor\tdecoy\tstatistic\ttau\ttau_ready\tadmitted\n"); }
        }
        std::fprintf(log, "%u\t%d\t%.6g\t%.6g\t%d\t%d\n",
                     precursor, is_decoy ? 1 : 0, stat, tau,
                     was_ready ? 1 : 0, admitted ? 1 : 0);
        std::fflush(log);
      }
    };

    std::mutex rejects_registry_mutex_;
    std::vector<PickerRejects*> rejects_registry_;

    struct RegisteredRejects : PickerRejects
    {
      RegisteredRejects()
      {
        std::lock_guard<std::mutex> g(rejects_registry_mutex_);
        rejects_registry_.push_back(this);
      }
    };
    thread_local RegisteredRejects rejects_;

    /// Every thread's counters, summed. The only correct way to read them.
    PickerRejects totalRejects()
    {
      PickerRejects t;
      std::lock_guard<std::mutex> g(rejects_registry_mutex_);
      for (const PickerRejects* r : rejects_registry_)
      {
        for (int c = 0; c < 2; ++c)
        {
          t.scans[c] += r->scans[c];
          t.too_few_present[c] += r->too_few_present[c];
          t.too_few_transitions[c] += r->too_few_transitions[c];
          t.below_corr[c] += r->below_corr[c];
          t.reference_zero[c] += r->reference_zero[c];
          t.not_local_max[c] += r->not_local_max[c];
          t.below_apex_evidence[c] += r->below_apex_evidence[c];
          t.outside_margin[c] += r->outside_margin[c];
          t.reached[c] += r->reached[c];
        }
        for (int c = 0; c < 2; ++c) { t.no_points[c] += r->no_points[c]; }
        for (int c = 0; c < 2; ++c) { t.empty_trace[c] += r->empty_trace[c]; }
        for (int c = 0; c < 2; ++c) { t.gate_c[c] += r->gate_c[c]; }
        for (int c = 0; c < 2; ++c) { t.few_excursions[c] += r->few_excursions[c]; }
        for (int c = 0; c < 2; ++c) { t.zero_trace[c] += r->zero_trace[c]; }
        t.too_few_at_apex += r->too_few_at_apex;
        t.masked_candidates += r->masked_candidates;
        t.masked_fragments += r->masked_fragments;
      }
      return t;
    }
  }

  struct PeakGroupScorer::Session::GateNull : NullCalibrationBody {};

  PeakGroupScorer::Session::Session(const Library& library, const Options& options)
    : gate_null_(std::make_shared<GateNull>()), library_(&library), options_(options)
  {
  }

  void PeakGroupScorer::Session::freezeGate(double tau) { gate_null_->freeze(tau); }

  void PeakGroupScorer::Session::add(const PrecursorChromatogram& chromatogram_in)
  {
    Result& result = result_;
    const Options& options = options_;
    const Library& library = *library_;
    const auto& p = library.precursors();
    const auto& t = library.transitions();
    const PrecursorChromatogram& chromatogram = chromatogram_in;   // rebound per candidate under a mask (v1.15)

    const std::size_t i = chromatogram.precursor;
    const std::uint32_t tb = chromatogram.transition_begin;
    const std::uint32_t tc = chromatogram.transition_count;
    // Every return below records why. `mark` is a no-op unless
    // -out_terminal_reasons asked for the table.
    // The EXTRACTOR's reasons win. It knows things this function cannot -- that
    // no isolation window covers the precursor, or that the prefilter dropped
    // it -- and if the scorer overwrites them the specific reason is replaced
    // by a vaguer one that is also true.
    //
    // That is not hypothetical. On full_v9, `few_points` came out at 1,119,490
    // and `no_window_coverage` at ZERO, and 100.0% of the `few_points` targets
    // turned out to lie outside every window's m/z range (median m/z 1468.3
    // against 525.8 for scored precursors; the windows stop at 1400.62 Th).
    // The whole bucket was the uncovered population wearing the wrong label,
    // which is exactly the kind of misattribution this table exists to prevent.
    const auto mark = [&](TerminalReason r) {
      if (!options.terminal_reason) { return; }
      const std::uint8_t prior = options.terminal_reason[i];
      if (prior == static_cast<std::uint8_t>(TerminalReason::NoWindowCoverage) ||
          prior == static_cast<std::uint8_t>(TerminalReason::PrefilterExcluded))
      { return; }
      options.terminal_reason[i] = static_cast<std::uint8_t>(r);
    };
    if (tc == 0) { mark(TerminalReason::NoTransitions); return; }

    const std::size_t points = chromatogram.pointCount(0);
    // The class is needed BEFORE these gates, not after: `reached` is
    // incremented past them, so an asymmetry here is invisible in every counter
    // downstream. One thread's numbers already put the entire target/decoy
    // imbalance in `reached` (137/254) while scans PER PRECURSOR were identical
    // (2014 both), so the divergence happens at exactly these two returns.
    const bool is_decoy = p.decoy[chromatogram.precursor] != 0;
    if (points < 3)
    { mark(TerminalReason::FewPoints); ++rejects_.no_points[is_decoy]; ++result.precursors_without_candidate; return; }

    // D8: each transition standardised against its own local noise before
    // summing, so no transition dominates by being loud and none is boosted
    // by what the library expects. See the option's comment for why library
    // weighting was rejected.
    std::size_t excursions = 0;
    const auto total = options.noise_normalised_picking
      ? noiseNormalisedTrace(chromatogram, points,
                             options.empty_trace_sigma, &excursions)
      : summedTrace(chromatogram, points);

    // Does anything rise above this precursor's own noise?
    //
    // The old test was `sum(median-subtracted trace) <= 0.0`. That sum is
    // ZERO-MEAN for a precursor with no real peak, so it was a coin flip on the
    // sign of noise and admitted about half of every absent precursor.
    // Measured on Astral: of precursors with usable points, 7.4% of targets and
    // 13.6% of decoys survived it -- the entire 1.85x decoy excess, from a gate
    // sitting exactly on the centre of the distribution where a 2.35 Th mean
    // m/z difference between the classes is enough to tip it.
    //
    // An inflated decoy null raises the 1% threshold above the few real
    // targets, which is what "the classifier trained and the threshold rejected
    // everything" looks like.
    //
    // -empty_trace_sigma 0 restores the old behaviour for comparison.
    // The oracle bypasses admission entirely -- that is the whole point of it.
    const bool oracled =
      options.oracle_rt != nullptr && std::isfinite(options.oracle_rt[i]);
    if (options.gate_alpha > 0.0 && !oracled)
    {
      // Gate C. Co-elution evidence against a threshold calibrated from the
      // run's own decoy null. Measured at real dimensions (12 traces x 200
      // cycles, Poisson noise) against its two predecessors on identical data:
      //
      //     gate           ABSENT      PRESENT
      //     A sum>0         90.1%       100.0%
      //     B 3sigma        99.0%       100.0%
      //     C co-elution     4.8%       100.0%
      std::size_t contributing = 0;
      const double m = coelutionEvidence(chromatogram, points,
                                         options.gate_smooth_half, &contributing);
      bool was_ready = true;
      bool admitted = true;
      if (options.gate_mode == "prominence")
      {
        // PER-PRECURSOR admission, from the statistic's own noise model rather
        // than from a run-wide decoy quantile.
        //
        // coelutionEvidence sums robust z-scores across `contributing`
        // transitions and averages over a window of w = 2*half+1 cycles, so
        // under white noise its null is N(0, contributing/w) and the scale is
        // sqrt(contributing/w). A threshold of k sigma therefore needs no null,
        // no calibration sample, and no dependence on library composition or
        // arrival order -- the three things measured wrong about the quantile
        // mode (alpha 0.05 delivering 15.3% decoy admission because tau is
        // fixed from the earliest-eluting 20,000 decoys).
        //
        // k is a LOOK-ELSEWHERE threshold, not a significance level: the
        // statistic is a maximum over ~M_eff independent positions, so for a
        // per-precursor false-proposal rate alpha_P,
        //     k = Phi^-1( (1-alpha_P)^(1/M_eff) )
        // which for alpha_P = 0.05 gives 3.1 at M_eff 50, 3.3 at 100, 3.5 at
        // 200. The default is the M_eff = 100 value. It is a starting point for
        // a sweep, not a derived constant -- prominence after smoothing does
        // not follow the point-height Gaussian model and real interference is
        // heavy-tailed and structured.
        //
        // This admits candidates; it does NOT assert the precursor is present.
        // Structured interference that is itself peak-shaped passes at any k,
        // and must be separated downstream by co-elution, library agreement,
        // mass accuracy, retention time and mobility. `max_candidates` remains
        // the resource bound.
        const double w = double(2 * options.gate_smooth_half + 1);
        const double sigma = contributing > 0
                               ? std::sqrt(double(contributing) / w) : 0.0;
        admitted = sigma > 0.0 && m >= options.gate_k * sigma;
      }
      else if (options.gate_alpha > 0.0)
      {
        was_ready = gate_null_->ready;
        // v1.17: the window's position in run seconds -- the midpoint of the
        // extracted cycle range, (lo+hi)/2 of the calibrated RT window after
        // clipping to the axis; the extractor's own centre is not carried here.
        // NaN (no axis) compares false and keeps the pre-v1.17 path.
        const double rt_centre = chromatogram.rt != nullptr
          ? 0.5 * (double(chromatogram.rt[0]) + double(chromatogram.rt[points - 1]))
          : std::numeric_limits<double>::quiet_NaN();
        admitted = gate_null_->admit(m, is_decoy,
                                     options.gate_calibration_n,
                                     options.gate_alpha,
                                     rt_centre, options.gate_calibration_rt_min);
      }
      gate_null_->note(options.gate_log_path, static_cast<std::uint32_t>(i),
                       is_decoy, m, admitted, was_ready);
      if (!admitted)
      { mark(TerminalReason::GateC); ++rejects_.empty_trace[is_decoy]; ++rejects_.gate_c[is_decoy]; ++result.precursors_without_candidate; return; }
    }
    else if (options.noise_normalised_picking && options.empty_trace_sigma > 0.0)
    {
      if (excursions < options.empty_trace_min_transitions)
      { mark(TerminalReason::FewExcursions); ++rejects_.empty_trace[is_decoy]; ++rejects_.few_excursions[is_decoy]; ++result.precursors_without_candidate; return; }
    }
    else
    {
      const double window_total = std::accumulate(total.begin(), total.end(), 0.0);
      if (window_total <= 0.0)
      { mark(TerminalReason::ZeroTrace); ++rejects_.empty_trace[is_decoy]; ++rejects_.zero_trace[is_decoy]; ++result.precursors_without_candidate; return; }
    }

    std::vector<MassAnchor> staged_anchors;
    const std::size_t first_group = result.groups.size();
    // Three pickers now. OpenSWATH's is the independent implementation
    // doc/07 step 2 requires; it picks on the summed trace by amplitude and
    // sets boundaries by signal-to-noise, so co-elution enters only later as a
    // score -- the opposite ordering to findCandidatesByCorrelation.
    std::vector<Candidate> openswath_candidates;
    if (options.openswath_picking)
    {
      std::vector<float> rt(chromatogram.cycles);
      for (std::uint32_t j = 0; j < chromatogram.cycles; ++j)
      { rt[j] = chromatogram.retentionTime(j); }
      for (const auto& c : pickOpenSwath(total, rt, options.openswath_sn,
                                         options.openswath_gauss,
                                         options.openswath_peak_width,
                                         options.max_candidates))
      {
        Candidate cd;
        cd.apex = c.apex; cd.left = c.left; cd.right = c.right;
        cd.apex_value = c.apex_value;
        // corr_sum stays 0: this picker never computes it, exactly as the
        // amplitude picker does not. CORR_SUM is then a constant column and the
        // constant-column guard drops it, which is the honest outcome -- a
        // feature this picker cannot supply must not be faked.
        openswath_candidates.push_back(cd);
      }
    }
    const auto amplitude_candidates = [&] {
      return findCandidates(smooth(total, options.boundary_smooth_half),
                            options.max_candidates, options.boundary_fraction,
                            options.peak_min_cycles, options.peak_max_half_cycles,
                            options.boundary_sigmas);
    };
    const auto coelution_candidates = [&] {
      // The picker's rejection counters are split by class so the stage at
      // which targets and decoys diverge can be SEEN. Astral yields 1.47x more
      // decoy peak groups than target ones from a balanced library and the
      // cause is unknown; totals cannot localise it, and the previous counters
      // could not even be summed -- they overcounted scan positions by 27.8%.
      ++rejects_.reached[is_decoy];
      return findCandidatesByCorrelation(chromatogram, rejects_,
                                         p.decoy[chromatogram.precursor] != 0,
                                         options.corr_half_window,
                                         options.min_corr_score, options.max_corr_diff,
                                         options.apex_evidence, options.smooth_half_width,
                                         options.boundary_fraction, options.max_candidates,
                                         options.peak_min_cycles,
                                         options.peak_max_half_cycles,
                                         options.boundary_smooth_half,
                                         options.boundary_sigmas,
                                         // null disables the library ordering, which is what
                                         // a library with no intensities should get: an
                                         // all-equal vector correlates with nothing and would
                                         // hand every position the same -1.
                                         t.library_intensity.empty()
                                           ? nullptr
                                           : t.library_intensity.data() + chromatogram.transition_begin,
                                         options.score_half_cycles,
                                         options.candidate_min_separation,
                                         options.select_library_weight);
    };
    const auto candidates = options.union_picking
      ? unionCandidates(chromatogram, coelution_candidates(),
                        options.openswath_picking ? openswath_candidates
                                                  : amplitude_candidates(),
                        options.corr_half_window, options.max_candidates)
      : options.openswath_picking
      ? openswath_candidates
      : options.coelution_picking
      ? coelution_candidates()
      : amplitude_candidates();
    // `candidates` is const above; the oracle needs to add to it.
    std::vector<Candidate> forced;
    const std::vector<Candidate>* use = &candidates;
    if (oracled)
    {
      const double want = options.oracle_rt[i];
      bool have = false;
      for (const auto& c : candidates)
      {
        if (c.apex < chromatogram.cycles &&
            std::abs(double(chromatogram.retentionTime(
                       static_cast<std::uint32_t>(c.apex))) - want)
              <= options.oracle_rt_tol)
        { have = true; break; }
      }
      if (!have)
      {
        // Nearest cycle to the oracle time, with the same boundary rule the
        // picker uses, so the synthesised group is shaped like a real one and
        // every sub-score downstream sees what it expects.
        std::size_t at = 0;
        double best = std::numeric_limits<double>::infinity();
        for (std::uint32_t j = 0; j < chromatogram.cycles; ++j)
        {
          const double d = std::abs(double(chromatogram.retentionTime(j)) - want);
          if (d < best) { best = d; at = j; }
        }
        if (best <= options.oracle_rt_tol && at < points)
        {
          // The RAW sum, not add()'s `total`: with noise-normalised picking
          // that one is in units of sigma, while findCandidatesByCorrelation
          // takes its boundaries and apex_value from raw intensity. A
          // synthesised group has to be measured on the same scale as a picked
          // one or every apex-valued sub-score reads it wrong.
          std::vector<double> raw(points, 0.0);
          for (std::uint32_t k = 0; k < tc; ++k)
          {
            const std::uint32_t m = chromatogram.pointCount(k);
            const float* pk = m ? chromatogram.trace(k) : nullptr;
            for (std::uint32_t j = 0; j < m && j < points; ++j) { raw[j] += pk[j]; }
          }
          Candidate cd;
          cd.apex = at;
          cd.apex_value = raw[at];
          cd.corr_sum = 0.0;
          const double floor_value = options.boundary_fraction * raw[at];
          std::size_t l = at, r = at;
          while (l > 0 && raw[l - 1] > floor_value) { --l; }
          while (r + 1 < points && raw[r + 1] > floor_value) { ++r; }
          cd.left = l; cd.right = r;
          forced = candidates;
          forced.push_back(cd);
          use = &forced;
        }
      }
    }
    if (use->empty())
    { mark(TerminalReason::NoCandidate); ++result.precursors_without_candidate; return; }

    // Library intensities, in the transition order the chromatograms use.
    std::vector<double> library_intensity_all(tc, 0.0);
    for (std::uint32_t k = 0; k < tc; ++k)
    {
      library_intensity_all[k] = t.library_intensity[tb + k];
    }
    // v1.15: transition index map. Identity unless a mask compacts this candidate's transitions.
    std::vector<std::uint32_t> kmap_identity(tc);
    for (std::uint32_t k = 0; k < tc; ++k) { kmap_identity[k] = k; }
    const std::vector<Options::MaskEntry>* mask_entries = nullptr;
    if (options.transition_mask != nullptr)
    {
      const auto it_m = options.transition_mask->find(static_cast<std::uint32_t>(i));
      if (it_m != options.transition_mask->end()) { mask_entries = &it_m->second; }
    }

    for (const auto& cand : *use)
    {
      // v1.15 INSTRUMENT: a candidate-scoped transition mask. When this candidate's apex RT lies
      // inside a listed interval, the listed transitions are ABSENT for its scoring: the block
      // below runs on a compacted VIEW of the chromatogram (same storage, fewer transitions) and a
      // compacted library-intensity vector; every library lookup goes through kmap. Off, or for an
      // unlisted candidate, the view IS the input and kmap is the identity: identical arithmetic.
      std::vector<std::uint64_t> mask_off;
      std::vector<std::uint32_t> mask_cnt, kmap_masked;
      std::vector<double> lib_masked;
      PrecursorChromatogram view = chromatogram_in;
      bool mask_active = false;
      if (mask_entries != nullptr && cand.apex < chromatogram_in.cycles)
      {
        const float apex_rt = chromatogram_in.retentionTime(static_cast<std::uint32_t>(cand.apex));
        std::vector<char> drop(tc, 0);
        for (const auto& e : *mask_entries)
        { if (e.k < tc && apex_rt >= e.rt_lo && apex_rt <= e.rt_hi) { drop[e.k] = 1; } }
        std::uint32_t kept = 0;
        for (std::uint32_t k = 0; k < tc; ++k) { if (!drop[k]) { ++kept; } }
        if (kept > 0 && kept < tc)
        {
          mask_active = true;
          for (std::uint32_t k = 0; k < tc; ++k)
          {
            if (drop[k]) { continue; }
            kmap_masked.push_back(k);
            mask_off.push_back(chromatogram_in.offset[k]);
            mask_cnt.push_back(chromatogram_in.count ? chromatogram_in.count[k] : 0u);
            lib_masked.push_back(library_intensity_all[k]);
          }
          view.transition_count = kept;
          view.offset = mask_off.data();
          view.count = chromatogram_in.count ? mask_cnt.data() : nullptr;
          ++rejects_.masked_candidates;
          rejects_.masked_fragments += tc - kept;
        }
      }
      const PrecursorChromatogram& chromatogram = mask_active ? view : chromatogram_in;
      const std::uint32_t tc = chromatogram.transition_count;
      const std::vector<double>& library_intensity = mask_active ? lib_masked : library_intensity_all;
      const std::vector<std::uint32_t>& kmap = mask_active ? kmap_masked : kmap_identity;
      // Per CANDIDATE, not per mass block. The mass block below is guarded by
      // `hi > lo`, so a single-cycle candidate skips it -- and used to inherit
      // whatever the previous candidate staged, committing it under this
      // candidate's group index and decoy flag. Measured on IH1 before the fix:
      // 26,868 anchors passed a target-group filter against 13,778 that were
      // actually flagged target, and the mismatched pairs put 876 control
      // residuals into a set that contained none, which made the fit refuse.
      staged_anchors.clear();
      // TWO intervals, because one interval was being asked to do two jobs that
      // want opposite things.
      //
      // The walked boundaries answer "how far does this peak extend", which is
      // the right question for QUANTIFICATION and for the retention-time range
      // the group reports. The sub-scores are asking something else: "does the
      // evidence at THIS position look like this peptide". For that, extent is
      // a liability -- every cycle the interval gains past the peak is a cycle
      // of neighbouring signal diluting the correlation, and the wider the
      // interval the less the score depends on where the candidate actually is.
      //
      // Measured on 11,728 paired candidates: separating a correct candidate
      // from a wrong one, a fixed apex-centred window beats the walked bounds
      // by 0.0047 of AUC, 95% CI [0.0011, 0.0082] bootstrapped over precursors,
      // 1,988 of 2,000 resamples favouring fixed. The sign reverses where the
      // walk found a genuinely broad peak (walked width >= 9: 0.838 walked
      // against 0.832 fixed), which is the population whose area the walked
      // bounds are still used for.
      //
      // Corrected from 0.0136, which was measured while `peak_min_cycles` was
      // briefly 5. It is a small effect, and the reason to keep the split is
      // not its size: it is that the two intervals answer different questions,
      // that DIA-NN separates them for the same reason, and that a regression
      // test can hold them apart. The change this refines -- a boundary floor
      // that was unreachable for 87.6% of real peptides -- is worth +0.253 on
      // the same measurement, and that is where the result actually lives.
      //
      // This is also what DIA-NN does, arrived at independently: its
      // discriminating correlations are computed over a fixed W = 2S+1 window
      // and its descent-derived borders are reported as RT_start/RT_stop for
      // quantification, not fed to the scores.
      const std::size_t quant_lo = cand.left, quant_hi = cand.right;
      const std::size_t lo = options.score_half_cycles
        ? (cand.apex > options.score_half_cycles ? cand.apex - options.score_half_cycles : 0)
        : cand.left;
      std::size_t hi = options.score_half_cycles
        ? std::min(cand.apex + options.score_half_cycles,
                   chromatogram.cycles ? std::size_t(chromatogram.cycles) - 1 : cand.right)
        : cand.right;
      // `width` is unsigned, so an inverted pair does not produce a small
      // window, it produces a 2^64 one. Reachable only if the apex is past the
      // last cycle, which nothing should produce -- and which is the reason to
      // spend a line on it rather than reason about whether anything does.
      if (hi < lo) { hi = lo; }
      const std::size_t width = hi - lo + 1;
      // retentionTime() indexes cycles; candidate bounds are already clamped to
      // them, but the clamp is kept explicit so a future boundary rule cannot
      // walk off the axis silently.
      const std::size_t n_cycles_minus1 =
        chromatogram.cycles ? chromatogram.cycles - 1 : 0;

      // Per-transition traces over the candidate's own boundaries, and the
      // observed intensity of each transition as its area there.
      std::vector<std::vector<double>> traces;
      std::vector<double> observed(tc, 0.0);
      traces.reserve(tc);
      for (std::uint32_t k = 0; k < tc; ++k)
      {
        const std::uint32_t n = chromatogram.pointCount(k);
        const float* points_k = n ? chromatogram.trace(k) : nullptr;
        std::vector<double> tr(width, 0.0);
        for (std::size_t j = 0; j < width; ++j)
        {
          const std::size_t at = lo + j;
          if (at < n) { tr[j] = points_k[at]; }
          observed[k] += tr[j];
        }
        traces.push_back(std::move(tr));
      }

      // D4: a per-transition local background, subtracted before anything
      // compares observed intensities to the library. Without it `observed`
      // is a raw area over a 60 s window and is dominated by baseline and
      // interference, which is why it correlated with nothing. Clamped at 0
      // rather than allowed negative: a weak real fragment sitting below its
      // own local median is absent evidence, not negative evidence.
      std::vector<double> corrected(tc, 0.0);
      // Retained per fragment: RT_SPREAD needs each fragment's own baseline to
      // place its centroid, and recomputing localBackground for that would be
      // the same scan twice. Named for the fragment axis -- a plain
      // `background` collides with the scalar one the sub-scores below use.
      std::vector<double> frag_background(tc, 0.0);
      std::size_t at_apex = 0;
      // -out_fragvec only. R1_ATAPEX_k is the PER-FRAGMENT form of the counter
      // below, and this is the only place it exists: `at_apex` collapses it to
      // a scalar in the same statement that computes it, so the 12 ATAPEX
      // columns cannot be recovered at the per-candidate block downstream.
      // Empty, and never touched, with the flag off.
      std::vector<double> frag_at_apex;
      if (options.fragvec) { frag_at_apex.assign(tc, 0.0); }
      for (std::uint32_t k = 0; k < tc; ++k)
      {
        const std::uint32_t n = chromatogram.pointCount(k);
        const float* points_k = n ? chromatogram.trace(k) : nullptr;
        std::vector<double> whole(n, 0.0);
        for (std::uint32_t j = 0; j < n; ++j) { whole[j] = points_k[j]; }
        const double bg = localBackground(whole, lo, hi);
        frag_background[k] = bg;
        corrected[k] = std::max(0.0, observed[k] - bg * static_cast<double>(width));
        if (cand.apex < n && points_k[cand.apex] > bg)
        {
          ++at_apex;
          if (options.fragvec) { frag_at_apex[k] = 1.0; }
        }
      }

      // RT_SPREAD: do this group's fragments agree about WHEN they elute?
      //
      // Each informative fragment gets a background-subtracted, intensity-
      // weighted retention-time centroid over the candidate's own boundaries;
      // the feature is the weighted scatter of those centroids about the
      // group's apex time. A peptide's fragments come from one ion packet and
      // agree; an interferent belongs to a different species and does not.
      //
      // The scatter is about the fragments' OWN weighted mean -- inter-fragment
      // dispersion, not displacement from `cand.apex`, which cancels out of the
      // arithmetic. Deliberate: the apex comes from the SUMMED trace, so one
      // loud interferent would move the reference and the measurement together,
      // and a feature about internal consistency must not be anchored to
      // something the interference can drag. "Fragments disagree with the
      // picked apex" is a different statistic and needs its own slot.
      //
      // Width floor: with one or two cycles every centroid is forced to the
      // same value and the scatter is 0 -- the BEST possible score -- so a
      // three-fragment noise spike clearing min_fragments_at_apex would look
      // like perfect co-elution.
      double rt_spread = std::numeric_limits<double>::quiet_NaN();
      if (width >= 3)
      {
        std::vector<double> offset, weight;
        offset.reserve(tc); weight.reserve(tc);
        for (std::uint32_t k = 0; k < tc; ++k)
        {
          if (!(corrected[k] > 0.0)) { continue; }
          double num = 0.0, den = 0.0;
          for (std::size_t j = 0; j < width; ++j)
          {
            const double v = traces[k][j] - frag_background[k];
            if (!(v > 0.0)) { continue; }
            const double rt = static_cast<double>(chromatogram.retentionTime(
              static_cast<std::uint32_t>(std::min<std::size_t>(lo + j, n_cycles_minus1))));
            num += rt * v; den += v;
          }
          // A fragment entirely at or below its own baseline across the peak
          // has no opinion about when it eluted. Excluded rather than given the
          // apex by default, which would fake agreement. (Unreachable while the
          // corrected[k] > 0 gate holds -- that implies sum(tr) > bg*width and
          // hence one positive v -- but kept so this loop is correct on its own
          // terms rather than on a gate twenty lines away.)
          if (!(den > 0.0)) { continue; }
          offset.push_back(num / den);
          // Area weighting estimates the dominant ion packet well and is poorly
          // matched to noticing ONE weak interfering fragment: for two clusters
          // the weighted variance scales as W1*W2/(W1+W2)^2, so a small
          // interferer's contribution vanishes with its area -- precisely the
          // case this feature exists to catch. Unweighted has the opposite
          // failure, giving a barely-detected noisy fragment an equal vote.
          // Sqrt is between. Swept, not assumed.
          const double a = corrected[k];
          weight.push_back(options.rt_spread_weight == RtSpreadWeight::None ? 1.0
                         : options.rt_spread_weight == RtSpreadWeight::Sqrt ? std::sqrt(a)
                         : a);
        }
        if (offset.size() >= options.min_rt_spread_fragments)
        {
          double w = 0.0, m = 0.0;
          for (std::size_t k = 0; k < offset.size(); ++k)
          { w += weight[k]; m += weight[k] * offset[k]; }
          if (w > 0.0)
          {
            m /= w;
            double var = 0.0;
            for (std::size_t k = 0; k < offset.size(); ++k)
            { const double d = offset[k] - m; var += weight[k] * d * d; }
            // Negated: lower scatter is better, and every other sub-score is
            // higher-is-better. Same convention as MASS_SPREAD and IM_SPREAD.
            rt_spread = -std::sqrt(var / w);
          }
        }
      }

      // D6/D8 gate: a peak group is a co-elution. One transition above its
      // own background is a spike, and emitting it as a candidate is what
      // let single-fragment interference into the score matrix.
      if (at_apex < options.min_fragments_at_apex) { ++rejects_.too_few_at_apex; continue; }

      // D1/D2/D3: self-pairs excluded, shape selected on the signed maximum,
      // lag capped to the trace, degenerate traces dropped. See score.h.
      Scoring::PairOptions popt;
      std::size_t usable = 0;
      const auto pairs = Scoring::allpairs_xcorr_ex(traces, options.max_delay, popt, &usable);
      double shape = 0.0, coelution = 0.0;
      for (const auto& pr : pairs)
      {
        shape += pr.value;
        coelution += std::abs(static_cast<double>(pr.delay));
      }
      if (!pairs.empty())
      {
        shape /= static_cast<double>(pairs.size());
        coelution /= static_cast<double>(pairs.size());
      }

      const double group_total = std::accumulate(observed.begin(), observed.end(), 0.0);

      // Background from outside the candidate, so a trace that is peak
      // everywhere does not report an impressive signal-to-noise.
      double background = 0.0;
      std::size_t background_n = 0;
      for (std::size_t j = 0; j < total.size(); ++j)
      {
        if (j < lo || j > hi) { background += total[j]; ++background_n; }
      }
      background = background_n ? background / static_cast<double>(background_n) : 0.0;

      PeakGroup g;
      g.precursor = static_cast<std::uint32_t>(i);
      g.decoy = p.decoy[i] != 0;
      g.apex_rt = chromatogram.retentionTime(static_cast<std::uint32_t>(cand.apex));
      // The WALKED bounds, not the scoring window: this pair is the group's
      // reported retention-time range, which downstream quantification
      // integrates over. Narrowing it to the scoring window would truncate the
      // area of every peak broader than 2*score_half_cycles+1 cycles.
      g.left_rt = chromatogram.retentionTime(static_cast<std::uint32_t>(quant_lo));
      g.right_rt = chromatogram.retentionTime(static_cast<std::uint32_t>(quant_hi));
      g.apex_intensity = static_cast<float>(cand.apex_value);

      g.sub_scores.assign(N_SUB_SCORES, 0.0);
      g.sub_scores[XCORR_SHAPE] = shape;
      // Negated so that, like every other column, larger is better. A
      // classifier would learn the sign either way, but a human reading a
      // weight vector should not have to remember which column is inverted.
      g.sub_scores[XCORR_COELUTION] = -coelution;
      g.sub_scores[LIBRARY_CORR] = libraryCorrelation(corrected, library_intensity);
      g.sub_scores[LIBRARY_DOTPROD] = dotProduct(corrected, library_intensity);

      // D6: the old group/window area ratio carried no library or
      // co-elution information, and a narrow decoy spike in an empty window
      // approached 1.0. This is the share of the group's background-corrected
      // area sitting in the fragments the library says are brightest -- which
      // a single-transition spike cannot satisfy however tall it is.
      {
        std::vector<std::size_t> order(tc);
        for (std::size_t k = 0; k < tc; ++k) { order[k] = k; }
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b_) {
          return library_intensity[a] > library_intensity[b_]; });
        const std::size_t top = std::max<std::size_t>(1, tc / 3);
        double top_area = 0.0, all_area = 0.0;
        for (std::size_t r = 0; r < tc; ++r)
        {
          all_area += corrected[order[r]];
          if (r < top) { top_area += corrected[order[r]]; }
        }
        g.sub_scores[INTENSITY_SCORE] = all_area > 0.0 ? top_area / all_area : 0.0;
      }

      // D5: the background floor is derived from the data, not from 1e-6.
      // With the constant, a decoy in an empty window got
      // log(apex / 1e-6) ~ 16-25 and outranked a real target on a real
      // baseline -- the floor did not guard the ratio, it inverted it.
      // The 0.01 fraction caps the ratio at 100:1, and MEASURED on IH1 every
      // real peak is above that: of 21,055 known-present precursors, log_sn
      // took 10 distinct values and its p10, p50 and p90 were all exactly
      // log(100) = 4.605. The floor that stopped decoys scoring high also
      // removed the entire dynamic range where present precursors live, making
      // this a constant for the class it is supposed to discriminate.
      // `std::min(10.0, ...)` below already bounds the runaway case, so the
      // fraction is now a knob rather than a hard 1%.
      const double floor_bg =
        std::max(background, options.log_sn_floor_frac * cand.apex_value);
      g.sub_scores[LOG_SN] = floor_bg > 0.0
        ? std::min(10.0, std::log(std::max(1e-12, cand.apex_value) / floor_bg))
        : 0.0;
      // NOTE: this is a CONSTANT 12.000 for right and wrong answers alike --
      // a trace is degenerate only if constant, and in a window where 45-60%
      // of points are non-zero none ever is. Replacing it with co-elution
      // depth (fragments whose own maximum coincides with the group apex) was
      // tried and measured: best-ranked-right moved 47.0->47.0, 43.2->43.0,
      // 40.7->40.8 at depths 3/10/25 -- nothing, because it is collinear with
      // XCORR_SHAPE -- while on-RT identifications fell 530->332. Reverted.
      // A replacement has to be orthogonal to the pair correlations, and
      // co-elution timing is not.
      g.sub_scores[USABLE_FRAGMENTS] = static_cast<double>(usable);

      // Both vectors normalised to unit sum first: RMSD on raw areas would
      // measure how intense the precursor is, not how well it matches.
      {
        double so = 0.0, sl = 0.0;
        for (std::uint32_t k = 0; k < tc; ++k) { so += corrected[k]; sl += library_intensity[k]; }
        double rmsd = 0.0;
        if (so > 0.0 && sl > 0.0)
        {
          for (std::uint32_t k = 0; k < tc; ++k)
          {
            const double d = corrected[k] / so - library_intensity[k] / sl;
            rmsd += d * d;
          }
          rmsd = std::sqrt(rmsd / static_cast<double>(tc));
        }
        // Negated so larger is better, as every other column is.
        g.sub_scores[LIBRARY_RMSD] = -rmsd;
      }

      {
        double y_area = 0.0, all_area = 0.0;
        for (std::uint32_t k = 0; k < tc; ++k)
        {
          all_area += corrected[k];
          if (t.type[tb + kmap[k]] == FragmentType::Y) { y_area += corrected[k]; }
        }
        g.sub_scores[YSERIES_SCORE] = all_area > 0.0 ? y_area / all_area : 0.0;
        g.sub_scores[FRAGMENT_COVERAGE] =
          tc > 0 ? static_cast<double>(at_apex) / static_cast<double>(tc) : 0.0;
      }

      // ---- features that were already computed and thrown away ----
      g.sub_scores[CORR_SUM] = cand.corr_sum;

      // P1 engine port (doc/80 retrospective): per-fragment evidence against
      // the interference-robust best-fragment reference, at THIS candidate.
      // Identical arithmetic for targets and decoys; zero-filled below the
      // width floor, like every other windowed sub-score.
      if (width >= 3 && tc >= 3)
      {
        std::vector<std::uint32_t> area_order(tc);
        std::iota(area_order.begin(), area_order.end(), 0u);
        std::sort(area_order.begin(), area_order.end(),
                  [&](std::uint32_t a, std::uint32_t b)
                  { return corrected[a] > corrected[b]; });
        const std::size_t top = std::min<std::size_t>(6, tc);
        const auto pearson0 = [&](const std::vector<double>& a,
                                  const std::vector<double>& b)
        {
          double ma = 0.0, mb = 0.0;
          for (std::size_t j = 0; j < width; ++j) { ma += a[j]; mb += b[j]; }
          ma /= static_cast<double>(width); mb /= static_cast<double>(width);
          double num = 0.0, va = 0.0, vb = 0.0;
          for (std::size_t j = 0; j < width; ++j)
          {
            const double x = a[j] - ma, y = b[j] - mb;
            num += x * y; va += x * x; vb += y * y;
          }
          return (va > 0.0 && vb > 0.0) ? num / std::sqrt(va * vb) : 0.0;
        };
        // The reference: the area-top-6 member best correlated with the other
        // five. A contaminated fragment cannot become the reference unless it
        // out-correlates the clean majority -- the robustness this feature
        // family exists for.
        std::size_t best_ref = 0; double best_sum = -1e18;
        for (std::size_t a = 0; a < top; ++a)
        {
          double s = 0.0;
          for (std::size_t b = 0; b < top; ++b)
          { if (a != b) { s += pearson0(traces[area_order[a]], traces[area_order[b]]); } }
          if (s > best_sum) { best_sum = s; best_ref = a; }
        }
        const auto& rref = traces[area_order[best_ref]];
        std::vector<double> refv(rref.begin(), rref.begin() + width);
        for (std::size_t j = 1; j + 1 < width; ++j)
        { refv[j] = 0.5 * rref[j] + 0.25 * (rref[j - 1] + rref[j + 1]); }
        std::vector<double> pc(tc, 0.0);
        for (std::uint32_t k = 0; k < tc; ++k) { pc[k] = pearson0(refv, traces[k]); }
        double refsum = 0.0;
        for (std::size_t a = 0; a < top; ++a) { refsum += pc[area_order[a]]; }
        g.sub_scores[REF_CORR_SUM] = refsum;
        std::sort(pc.begin(), pc.end(), [](double a, double b) { return a > b; });
        for (std::size_t k = 0; k < 12; ++k)
        {
          g.sub_scores[static_cast<std::size_t>(REF_CORR_1) + k] =
            k < pc.size() ? pc[k] : 0.0;
        }
        double tot_area = 0.0;
        for (std::uint32_t k = 0; k < tc; ++k)
        { tot_area += std::max(corrected[k], 0.0); }
        if (tot_area > 0.0)
        {
          std::vector<double> sh(tc, 0.0);
          for (std::uint32_t k = 0; k < tc; ++k)
          { sh[k] = std::max(corrected[k], 0.0) / tot_area; }
          std::sort(sh.begin(), sh.end(), [](double a, double b) { return a > b; });
          for (std::size_t k = 0; k < 6; ++k)
          {
            g.sub_scores[static_cast<std::size_t>(SIG_SHARE_1) + k] =
              k < sh.size() ? sh[k] : 0.0;
          }
        }
      }

      g.sub_scores[PEAK_WIDTH_RATIO] = static_cast<double>(width);

      // Library 1/K0 against what the run observed for this precursor. NaN
      // wherever either side lacks mobility -- an absent measurement is not a
      // zero deviation, and a placeholder here would be a feature the
      // classifier weights as evidence.
      // The run's OBSERVED 1/K0 for this candidate, from the extractor's
      // intensity-weighted mobility planes.
      //
      // Until now this read `options.observed_im`, which the tool never
      // assigned, so IM_DELTA was all-NaN on EVERY run and the constant-column
      // guard dropped it every time. Nothing in the pipeline measured observed
      // mobility per precursor: `MobilityAnchor` carries only {precursor, rt,
      // decoy}, and the mobility calibration's probe keeps its residuals to
      // itself. The planes are that missing measurement.
      double im_obs = std::numeric_limits<double>::quiet_NaN();
      if (chromatogram.im_num != nullptr && chromatogram.im_den != nullptr && hi > lo)
      {
        double num = 0.0, den = 0.0;
        for (std::uint32_t k = 0; k < tc; ++k)
        {
          const std::uint32_t n = chromatogram.pointCount(k);
          const std::ptrdiff_t off = chromatogram.trace(k) - chromatogram.points;
          const float* inum = chromatogram.im_num + off;
          const float* iden = chromatogram.im_den + off;
          for (std::size_t j = lo; j <= hi && j < n; ++j)
          {
            if (iden[j] > 0.0f) { num += inum[j]; den += iden[j]; }
          }
        }
        if (den > 0.0) { im_obs = num / den; }

        // Per-FRAGMENT observed mobility, and its scatter across the group.
        // Same construction as mass_ppm_spread: one number per fragment (the
        // intensity-weighted mean over the cycles it was seen in), then a MAD
        // across fragments. Averaging within a fragment is legitimate -- those
        // cells measure one ion -- whereas averaging across fragments is what
        // would hide the very disagreement being looked for.
        std::vector<double> per_fragment_im;
        per_fragment_im.reserve(tc);
        for (std::uint32_t k = 0; k < tc; ++k)
        {
          const std::uint32_t n = chromatogram.pointCount(k);
          const std::ptrdiff_t off = chromatogram.trace(k) - chromatogram.points;
          const float* inum = chromatogram.im_num + off;
          const float* iden = chromatogram.im_den + off;
          double fn = 0.0, fd = 0.0;
          for (std::size_t j = lo; j <= hi && j < n; ++j)
          { if (iden[j] > 0.0f) { fn += inum[j]; fd += iden[j]; } }
          if (fd > 0.0) { per_fragment_im.push_back(fn / fd); }
        }
        if (per_fragment_im.size() >= 3)
        {
          const std::size_t c = per_fragment_im.size() / 2;
          std::nth_element(per_fragment_im.begin(), per_fragment_im.begin() + c,
                           per_fragment_im.end());
          const double centre = per_fragment_im[c];
          std::vector<double> ad;
          ad.reserve(per_fragment_im.size());
          for (const double v : per_fragment_im) { ad.push_back(std::abs(v - centre)); }
          const std::size_t m = ad.size() / 2;
          std::nth_element(ad.begin(), ad.begin() + m, ad.end());
          g.im_spread = static_cast<float>(1.4826 * ad[m]);
          g.im_frags = static_cast<std::uint8_t>(
            std::min<std::size_t>(per_fragment_im.size(), 255));
        }
      }
      // The external vector still wins when a caller supplies one, so the
      // existing plumbing point is not silently ignored.
      if (options.observed_im != nullptr && i < options.observed_im->size() &&
          std::isfinite((*options.observed_im)[i]))
      { im_obs = static_cast<double>((*options.observed_im)[i]); }

      g.observed_im = static_cast<float>(im_obs);
      if (i < p.im.size() && std::isfinite(p.im[i]) && std::isfinite(im_obs))
      {
        g.sub_scores[IM_DELTA] = std::fabs(static_cast<double>(p.im[i]) - im_obs);
      }
      else { g.sub_scores[IM_DELTA] = std::numeric_limits<double>::quiet_NaN(); }

      // The peak group's own fragment mass deviation, when the extractor kept it.
      //
      // Recorded on the group rather than turned into a sub-score: a mass
      // deviation is a property of the INSTRUMENT, and feeding it to the
      // discriminant would let the classifier learn "this run's fragments sit
      // at -3 ppm" and reject correct identifications for being well
      // calibrated. It exists to fit a recalibration, and is read only for
      // groups the FDR has already accepted.
      // MASS_SURVIVAL: how much of each fragment's MATCHED intensity is still
      // there when the mass window is tightened.
      //
      // `ppm_den` is sum(intensity) over the peaks that matched this cell, so
      // it is the matched intensity itself and needs no separate weight. A cell
      // whose intensity-weighted mean deviation falls outside the tightened
      // window is treated as lost. That is an approximation -- a cell holding
      // both an on-mass and an off-mass peak has a mean that hides the split --
      // and it is the same approximation the two-plane representation already
      // forces on MASS_ACCURACY, so it introduces no new inaccuracy here.
      // `hi > lo` for the same reason the block below it requires it: at a
      // single cycle the per-fragment ratio degenerates to 0 or 1 and enters
      // the classifier alongside fractions taken over five. The feature
      // reverted yesterday lacked exactly this guard.
      if (options.mass_survival_ppm > 0.0 && hi > lo &&
          chromatogram.ppm_num != nullptr && chromatogram.ppm_den != nullptr)
      {
        double acc = 0.0;
        std::size_t n_frag = 0;
        for (std::uint32_t k = 0; k < tc; ++k)
        {
          const std::uint32_t n = chromatogram.pointCount(k);
          const std::ptrdiff_t off = chromatogram.trace(k) - chromatogram.points;
          const float* num = chromatogram.ppm_num + off;
          const float* den = chromatogram.ppm_den + off;
          double full = 0.0, tight = 0.0;
          for (std::size_t j = lo; j <= hi && j < n; ++j)
          {
            if (!(den[j] > 0.0f)) { continue; }
            const double d = double(den[j]);
            full += d;
            // Centred on the run's fragment deviation, not on zero. The
            // calibration offset is normally folded into the query m/z, which
            // centres the residual -- but when the mass-calibration gate FAILS
            // the offset is set to 0 with a 50 ppm window and the distribution
            // is left where the instrument put it. On IH1 that is about -9 ppm,
            // and an absolute window about zero would then score a perfectly
            // real fragment near 0 while a decoy matching noise uniformly
            // across the window scores higher: the column INVERTS. MASS_ACCURACY
            // re-centres for precisely this reason.
            const double dev = double(num[j]) / d - options.mass_survival_centre;
            if (std::fabs(dev) <= options.mass_survival_ppm)
            { tight += d; }
          }
          // A fragment that matched nothing abstains rather than scoring 0:
          // "no peak to lose" is not "the peak did not survive".
          if (full > 0.0) { acc += tight / full; ++n_frag; }
        }
        if (n_frag > 0)
        { g.sub_scores[MASS_SURVIVAL] = acc / static_cast<double>(n_frag); }
        else
        { g.sub_scores[MASS_SURVIVAL] = std::numeric_limits<double>::quiet_NaN(); }
      }
      else
      {
        g.sub_scores[MASS_SURVIVAL] = std::numeric_limits<double>::quiet_NaN();
      }

      // The control. splitmix64 on (precursor, apex): deterministic, so a
      // repeat run reproduces it exactly, and uniform in [0,1) with no relation
      // to whether this candidate is a target or a decoy.
      // Two kinds of null, and the second is the one that tests anything.
      //
      // `hash` is a column uniform in label. It is index 20 of 21, and split
      // ties break toward the LOWEST feature index (gbt.h), so it is maximally
      // disadvantaged for winning a tie -- and at depth 2 it won zero splits in
      // any fold. Bit-identical output then proves only that a column which is
      // never selected changes nothing, which is a fixed point rather than
      // stability.
      //
      // `dup` is the test that matters: an epsilon-jittered copy of CORR_SUM,
      // the highest-attribution feature in the model. It carries NO information
      // the model does not already have, and it WILL win splits, because it is
      // very nearly as good as the column it copies. If output moves under it,
      // then "insensitive to an added column" was never true -- it was
      // "insensitive to a column that loses every gain comparison".
      if (options.null_feature)
      {
        std::uint64_t z = (std::uint64_t(i) << 20) ^ std::uint64_t(cand.apex)
                        ^ (options.null_feature_seed * 0xD1B54A32D192ED03ULL);
        z += 0x9E3779B97F4A7C15ULL;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        z ^= z >> 31;
        const double u =
          static_cast<double>(z >> 11) / static_cast<double>(1ULL << 53);
        if (options.null_feature_dup)
        {
          // Relative jitter, so the copy tracks CORR_SUM across its whole range
          // rather than being swamped by an absolute epsilon at large values.
          const double base = g.sub_scores[CORR_SUM];
          g.sub_scores[NULL_CONTROL] = base * (1.0 + 1e-6 * (u - 0.5));
        }
        else
        {
          g.sub_scores[NULL_CONTROL] = u;
        }
      }
      else
      {
        g.sub_scores[NULL_CONTROL] = std::numeric_limits<double>::quiet_NaN();
      }

      if (chromatogram.ppm_num != nullptr && chromatogram.ppm_den != nullptr && hi > lo)
      {
        std::vector<double> dev, per_fragment, cell;
        dev.reserve((hi - lo + 1) * tc);
        per_fragment.reserve(tc);
        for (std::uint32_t k = 0; k < tc; ++k)
        {
          const std::uint32_t n = chromatogram.pointCount(k);
          const std::ptrdiff_t off = chromatogram.trace(k) - chromatogram.points;
          const float* num = chromatogram.ppm_num + off;
          const float* den = chromatogram.ppm_den + off;
          cell.clear();
          for (std::size_t j = lo; j <= hi && j < n; ++j)
          {
            // A zero denominator is the unambiguous "no peak matched here" --
            // which the single-plane version could not distinguish from a real
            // deviation of exactly 0.000 ppm.
            if (den[j] > 0.0f) { cell.push_back(num[j] / den[j]); }
          }
          dev.insert(dev.end(), cell.begin(), cell.end());
          // One number per FRAGMENT: the median over the cycles it was seen in.
          // Averaging within a fragment is legitimate -- those cells measure the
          // same ion at the same m/z -- whereas averaging ACROSS fragments is
          // what would hide the scatter the window has to accommodate.
          if (!cell.empty())
          {
            const std::size_t c = cell.size() / 2;
            std::nth_element(cell.begin(), cell.begin() + c, cell.end());
            per_fragment.push_back(cell[c]);

            // The same median, kept WITH the m/z it belongs to. `per_fragment`
            // above collapses to a scatter; a calibration has to be fitted
            // against the axis it varies on, and that axis is destroyed the
            // moment the fragment identity is dropped.
            //
            // Staged, not committed: this candidate can still be rejected below
            // (min_library_corr), and an anchor pointing at a group that was
            // never pushed would index the wrong group after finish().
            if (options.collect_mass_anchors)
            {
              const double fmz = fromFixed(t.product_mz[tb + kmap[k]]);
              if (fmz > 0.0)
              {
                // Undo the correction the extractor applied to this
                // transition's target, so the residual is against the
                // UNCORRECTED theoretical m/z. See `applied_ppm_offset`.
                double applied = options.applied_ppm_offset;
                if (options.applied_ppm_log_slope != 0.0 && options.applied_ppm_ref_mz > 0.0)
                { applied += options.applied_ppm_log_slope * std::log(fmz / options.applied_ppm_ref_mz); }
                if (options.applied_ppm_slope_per_1000 != 0.0)
                { applied += options.applied_ppm_slope_per_1000 * (fmz - options.applied_ppm_ref_mz) / 1000.0; }

                MassAnchor a;
                a.residual.mz = static_cast<float>(fmz);
                a.residual.rt = g.apex_rt;
                // First order: the exact inverse carries a cross term of
                // r*applied/1e6, which at 10 ppm each is 1e-4 ppm.
                a.residual.ppm = static_cast<float>(cell[c] + applied);
                // The OBSERVED apex intensity of this fragment, not the library
                // intensity: an intensity-dependent mass error is a property of
                // how many ions arrived, and the library's number is a
                // prediction about a different run.
                a.residual.intensity = n ? chromatogram.trace(k)[
                  std::min<std::size_t>(cand.apex, n - 1)] : 0.0f;
                // The precursor's mobility. A fragment has none of its own.
                a.residual.im = (i < p.im.size()) ? p.im[i]
                                : std::numeric_limits<float>::quiet_NaN();
                // NOT set here, and not at commit either: `decoy` is a property
                // of the GROUP, and storing it a second time on the residual is
                // what let the two disagree. Filled from the group at read time.
                staged_anchors.push_back(a);
              }
            }
          }
        }
        if (!dev.empty())
        {
          const std::size_t h = dev.size() / 2;
          std::nth_element(dev.begin(), dev.begin() + h, dev.end());
          g.mass_ppm = static_cast<float>(dev[h]);
          g.mass_ppm_n = static_cast<std::uint16_t>(std::min<std::size_t>(dev.size(), 65535));
        }
        if (per_fragment.size() >= 3)
        {
          const std::size_t c = per_fragment.size() / 2;
          std::nth_element(per_fragment.begin(), per_fragment.begin() + c, per_fragment.end());
          const double centre = per_fragment[c];
          std::vector<double> abs_dev;
          abs_dev.reserve(per_fragment.size());
          for (const double v : per_fragment) { abs_dev.push_back(std::abs(v - centre)); }
          const std::size_t m = abs_dev.size() / 2;
          std::nth_element(abs_dev.begin(), abs_dev.begin() + m, abs_dev.end());
          // 1.4826 makes the MAD an estimate of sigma for a Gaussian, so the
          // width the driver derives from it can be stated in sigmas.
          g.mass_ppm_spread = static_cast<float>(1.4826 * abs_dev[m]);
          g.mass_ppm_frags = static_cast<std::uint8_t>(
            std::min<std::size_t>(per_fragment.size(), 255));
        }
      }

      // MASS_SPREAD is final here: a within-group scatter needs no context.
      // MASS_ACCURACY holds the RAW deviation for now and is re-centred against
      // the run's median in finish(), once every group has been seen.
      // Lower is better, so negate -- same convention as MASS_SPREAD.
      g.sub_scores[IM_SPREAD] = std::isfinite(g.im_spread)
        ? -static_cast<double>(g.im_spread)
        : std::numeric_limits<double>::quiet_NaN();
      g.sub_scores[MASS_SPREAD] = std::isfinite(g.mass_ppm_spread)
        ? -static_cast<double>(g.mass_ppm_spread)
        : std::numeric_limits<double>::quiet_NaN();
      g.sub_scores[MASS_ACCURACY] = std::isfinite(g.mass_ppm)
        ? static_cast<double>(g.mass_ppm)
        : std::numeric_limits<double>::quiet_NaN();
      // Already negated where it was computed; NaN passes through so the
      // constant-column guard can drop it on a run where it never fires.
      g.sub_scores[RT_SPREAD] = rt_spread;

      // MS1_COELUTION: does the PRECURSOR rise and fall with its fragments?
      //
      // Over the candidate's own cycles, pair each cycle's summed fragment
      // intensity with the MS1 monoisotopic intensity at the nearest MS1 bin,
      // and correlate. The MS1 grid is ~1.8 s on IH1 against a ~0.4 s MS2
      // cycle, so several cycles map to one bin -- that is a real resolution
      // limit of the survey scan, not an approximation to be apologised for,
      // and a 20-30 s peak still spans ~15 bins.
      //
      // NaN, not zero, when there is no MS1 or no signal: zero is a legitimate
      // correlation (a precursor whose trace is flat where the fragments peak
      // is EVIDENCE AGAINST), and collapsing the two would feed the classifier
      // a placeholder dressed as a measurement.
      {
        double r = std::numeric_limits<double>::quiet_NaN();
        if (options.ms1 != nullptr && !options.ms1->empty() &&
            i < options.ms1->precursors() && hi > lo)
        {
          // `total` is already the per-cycle summed fragment intensity and
          // [lo, hi] the candidate's own cycle bounds -- reuse both rather than
          // recomputing a second, subtly different fragment sum.
          std::vector<double> f, m;
          f.reserve(hi - lo + 1);
          m.reserve(hi - lo + 1);
          for (std::size_t j = lo; j <= hi && j < total.size(); ++j)
          {
            const std::size_t b = options.ms1->binFor(chromatogram.retentionTime(
              static_cast<std::uint32_t>(j)));
            f.push_back(total[j]);
            m.push_back(static_cast<double>(options.ms1->at(i, b)));
          }
          // Refuse a correlation that would be computed from a handful of
          // points: over 3 or 4 cycles almost anything correlates, and the
          // classifier cannot tell a well-supported 0.9 from a lucky one.
          bool any = false;
          for (const double v : m) { if (v > 0.0) { any = true; break; } }
          if (any && f.size() >= 5) { const double c = pearson(f, m);
                                      if (std::isfinite(c)) { r = c; } }

          // WHICH guard kills it? On Astral this sub-score came out constant
          // across all 96,259 rows and was dropped as carrying no information --
          // and it is the one feature measured to discriminate (13.7x top bin),
          // so "it is NaN" is not a good enough answer. Atomics, not
          // thread_local: a thread_local counter read from one thread already
          // produced a retracted conclusion in this file once.
          ms1_no_signal.fetch_add(any ? 0 : 1, std::memory_order_relaxed);
          ms1_too_short.fetch_add(f.size() >= 5 ? 0 : 1, std::memory_order_relaxed);
          if (any && f.size() >= 5)
          {
            // Zero variance in the MS1 leg makes Pearson undefined however much
            // signal is present -- distinct from "no signal", and the case that
            // fires if several MS2 cycles share one MS1 bin.
            bool flat = true;
            for (std::size_t k = 1; k < m.size(); ++k)
            { if (m[k] != m[0]) { flat = false; break; } }
            ms1_flat.fetch_add(flat ? 1 : 0, std::memory_order_relaxed);
            ms1_ok.fetch_add(std::isfinite(pearson(f, m)) ? 1 : 0, std::memory_order_relaxed);
            ms1_bins_spanned.fetch_add(distinctBins_(chromatogram, lo, hi, *options.ms1),
                                       std::memory_order_relaxed);
            ms1_spans.fetch_add(1, std::memory_order_relaxed);
          }
        }
        else
        {
          // Split the guard. The first cut said "11,877 no MS1 available, 0
          // everything else" while the run had just built a 148 GiB MS1 matrix
          // over 9,983,789 precursors -- so the composite condition is useless
          // and each term has to be counted on its own.
          ms1_unavailable.fetch_add(1, std::memory_order_relaxed);
          if (options.ms1 == nullptr) { ms1_null.fetch_add(1, std::memory_order_relaxed); }
          else
          {
            if (options.ms1->empty()) { ms1_empty.fetch_add(1, std::memory_order_relaxed); }
            if (i >= options.ms1->precursors())
            { ms1_index_oob.fetch_add(1, std::memory_order_relaxed); }
            if (!(hi > lo)) { ms1_degenerate_span.fetch_add(1, std::memory_order_relaxed); }
          }
        }
        g.sub_scores[MS1_COELUTION] = r;
      }
      // ---- -out_fragvec: rung (i) of the sealed fragment-evidence contract ----
      //
      // 78 float32 per candidate: six 12-long vectors by LIBRARY-INTENSITY rank
      // (LOGAREA, SHARE, LOGRATIO, ATAPEX, MEAS, ABSENT) and six scalars
      // (N_MEAS, N_ABSENT, N_PRESENT, LOGTOT, LIB_CORR, LIB_CORR_LOO). Every
      // input is a local the scorer already computed and already discards; the
      // definition is analysis77/pick/wf_v33_fragvec_contract.md s.5.1 and the
      // arithmetic is wf_v33_fragvec_features.py's `row_features`, which is
      // what the measured result was produced with.
      //
      // Staged on the stack and appended at the push_back below, not here: the
      // `min_library_corr` gate a few lines down `continue`s AFTER this point,
      // and a row emitted for a candidate that never becomes a group would put
      // the export permanently out of step with -out.
      //
      // TWO PLACES THIS AND THE REFERENCE BUILDER CAN DIVERGE. Both are
      // MEASURED UNREACHABLE on the library the arm runs (dn_pred_cam.parquet,
      // 58,576,095 transitions), so they are recorded here to make a future
      // gate-A4 disagreement diagnosable in one step rather than fixed blind:
      //   (1) NaN. `corrected[k]` is built with std::max(0.0, x) above, which
      //       returns 0.0 for a NaN x, while the builder's np.maximum returns
      //       NaN. A single non-finite exported intensity would therefore give
      //       LOGAREA 0 here and NaN there (ABSENT agrees by accident: NaN > 0
      //       is false either way). That line is EXISTING scoring arithmetic --
      //       LIBRARY_CORR and the mass residuals read the same `corrected` --
      //       so it must not be changed for this export's convenience; the
      //       reference is what would have to move. Measured: 0 non-numeric and
      //       0 negative Intensity over 4,294,710 exported points.
      //   (2) A transition the extractor could not place (Product.Mz invalid)
      //       still occupies a rank HERE -- tc counts it, it enters tot, libsum
      //       and both Pearsons with corrected == 0 -- but contributes no row to
      //       -out_chrom, so a reference that rebuilds tc from the export would
      //       be shifted by one from that rank on. Measured: 0 such transitions
      //       (Product.Mz in [200.015427, 1799.998901], no NaN, none <= 0).
      float fv[N_FRAGVEC];
      if (options.fragvec)
      {
        // Its OWN permutation, and a STABLE one. The D6 block above also orders
        // by library intensity, but that order is `std::sort` (unstable) and is
        // destroyed with its braced block. Reusing it would assign ranks
        // differently from the reference builder's stable sort on every
        // Relative.Intensity tie, and the resulting permuted columns would read
        // as an arithmetic disagreement rather than as a different order.
        std::vector<std::uint32_t> lib_order(tc);
        std::iota(lib_order.begin(), lib_order.end(), 0u);
        std::stable_sort(lib_order.begin(), lib_order.end(),
                         [&](std::uint32_t a, std::uint32_t b_)
                         { return library_intensity[a] > library_intensity[b_]; });
        const std::size_t m12 = std::min<std::size_t>(tc, 12);
        double tot = 0.0, libsum = 0.0;
        for (std::uint32_t k = 0; k < tc; ++k)
        { tot += corrected[k]; libsum += library_intensity[k]; }
        const float nanf = std::numeric_limits<float>::quiet_NaN();
        double n_absent = 0.0;
        for (std::size_t r = 0; r < 12; ++r)
        {
          if (r >= m12)
          {
            // Rank not measured. NaN in the four value columns, ZERO in the two
            // indicators -- the one place the contract's two missing-data
            // conventions differ, and the reason MEAS/ABSENT are never NaN.
            fv[r] = nanf; fv[12 + r] = nanf; fv[24 + r] = nanf; fv[36 + r] = nanf;
            fv[48 + r] = 0.0f; fv[60 + r] = 0.0f;
            continue;
          }
          const std::uint32_t k = lib_order[r];
          const double share = tot > 0.0 ? corrected[k] / tot : 0.0;
          const double libshare = libsum > 0.0
            ? library_intensity[k] / libsum
            : 1.0 / static_cast<double>(tc);
          // A5: "measured absent" is a zero background-corrected 5-cycle AREA.
          // It says nothing about the raw trace, which is why ATAPEX above is
          // not forced to 0 here.
          const double absent = corrected[k] > 0.0 ? 0.0 : 1.0;
          n_absent += absent;
          fv[r] = static_cast<float>(std::log1p(corrected[k]));
          fv[12 + r] = static_cast<float>(share);
          fv[24 + r] = static_cast<float>(std::log((share + 1e-3) / (libshare + 1e-3)));
          fv[36 + r] = static_cast<float>(frag_at_apex[k]);
          fv[48 + r] = 1.0f;
          fv[60 + r] = static_cast<float>(absent);
        }
        fv[72] = static_cast<float>(m12);
        fv[73] = static_cast<float>(n_absent);
        fv[74] = static_cast<float>(static_cast<double>(m12) - n_absent);
        fv[75] = static_cast<float>(std::log1p(tot));
        // The same helper the LIBRARY_CORR sub-score uses: fewer than four
        // fragments with a positive corrected area report 0, not a Pearson over
        // three points.
        fv[76] = static_cast<float>(libraryCorrelation(corrected, library_intensity));
        // Leave-one-transition-out, FLOORED AT 0 by starting the max there --
        // the builder initialises at 0.0 and never lowers it, so a group whose
        // every LOO correlation is negative reports 0 rather than its maximum.
        // O(tc^2): tc Pearsons of tc-1 points, ~144 multiplies at tc = 12.
        double loo = 0.0;
        if (tc > 1)
        {
          std::vector<double> obs_loo(tc - 1), lib_loo(tc - 1);
          for (std::uint32_t drop = 0; drop < tc; ++drop)
          {
            std::size_t at = 0;
            for (std::uint32_t k = 0; k < tc; ++k)
            {
              if (k == drop) { continue; }
              obs_loo[at] = corrected[k];
              lib_loo[at] = library_intensity[k];
              ++at;
            }
            loo = std::max(loo, libraryCorrelation(obs_loo, lib_loo));
          }
        }
        fv[77] = static_cast<float>(loo);
      }

      // A candidate whose spectrum does not resemble the library is not this
      // peptide, wherever it eluted.
      //
      // Measured on IH1 against DIA-NN's confident set, splitting our own
      // q<=0.01 calls by whether they land within 30 s of the true apex:
      //
      //   group             n     median library_corr   frac > 0.5
      //   on-RT targets   1360             0.582           54.1%
      //   off-RT targets  5724            -0.036           12.5%
      //   decoys          7007            -0.032           12.0%
      //
      // Off-RT targets and decoys are the SAME population on this feature. That
      // is why target-decoy FDR cannot see them: every target is trained as a
      // potential positive, so the classifier learns whatever separates off-RT
      // targets from decoys -- signal presence -- rather than library
      // agreement. A target at the wrong retention time still sits on real
      // co-eluting ions from a real peptide, so it looks alive on log_sn and
      // intensity_score in a way a shuffled decoy never does.
      //
      // This gate is safe for the FDR by the criterion doc/08 states for any
      // selection upstream of it: it must be label-symmetric. It is, and that
      // is measured rather than assumed -- 12.5% of off-RT targets and 12.0% of
      // decoys survive a 0.5 cut, while 54.1% of on-RT targets do. Targets and
      // decoys traverse identical code here.
      //
      // Off by default: it changes which candidates exist, so it must be turned
      // on deliberately and its effect measured, not inherited.
      if (options.min_library_corr > -1.0 &&
          g.sub_scores[LIBRARY_CORR] < options.min_library_corr)
      {
        ++result.candidates_below_library_corr;
        continue;
      }
      // Commit this candidate's staged mass anchors, now that it is certain to
      // become a group and its index is known. `decoy` is taken from the group
      // rather than the fragment: a residual's status is the status of the
      // sequence it was matched against.
      if (options.collect_mass_anchors && !staged_anchors.empty())
      {
        const std::uint32_t gi = static_cast<std::uint32_t>(result.groups.size());
        for (MassAnchor& a : staged_anchors)
        {
          if (result.mass_anchors.size() >= options.max_mass_anchors)
          { ++result.mass_anchors_dropped; continue; }
          a.group = gi;
          result.mass_anchors.push_back(a);
        }
        staged_anchors.clear();
      }
      // In lockstep with `groups`, and appended at the SAME statement so no
      // path can push one without the other.
      if (options.fragvec)
      { result.fragvec.insert(result.fragvec.end(), fv, fv + N_FRAGVEC); }
      result.groups.push_back(std::move(g));
    }

    // Margin over this precursor's own runner-up, on the detector's correlation
    // sum. Computed here because it is the only place all of one precursor's
    // candidates are in hand; a per-candidate score cannot express "this
    // precursor had one obvious answer" versus "three equally plausible ones".
    mark(result.groups.size() > first_group ? TerminalReason::Scored
                                           : TerminalReason::AllCandidatesDropped);
    if (result.groups.size() > first_group)
    {
      double best = 0.0, second = 0.0;
      for (std::size_t g = first_group; g < result.groups.size(); ++g)
      {
        const double v = result.groups[g].sub_scores[CORR_SUM];
        if (v > best) { second = best; best = v; }
        else if (v > second) { second = v; }
      }
      const bool alone = (result.groups.size() - first_group) == 1;
      for (std::size_t g = first_group; g < result.groups.size(); ++g)
      {
        const double v = result.groups[g].sub_scores[CORR_SUM];
        // The sole candidate is compared against nothing, not against zero:
        // giving it its full corr_sum as a margin would make "only one peak
        // found" look like overwhelming evidence.
        result.groups[g].sub_scores[CANDIDATE_MARGIN] =
          alone ? 0.0 : (v >= best ? best - second : v - best);
        // P1: the rest of the competition context. Rank on the same CORR_SUM
        // this block already reads (1 = best), and the block size itself --
        // "the only peak found" and "best of seven near-ties" are different
        // evidence, and no per-candidate score can express either.
        std::size_t rank = 1;
        for (std::size_t h = first_group; h < result.groups.size(); ++h)
        { if (result.groups[h].sub_scores[CORR_SUM] > v) { ++rank; } }
        result.groups[g].sub_scores[CAND_RANK] = static_cast<double>(rank);
        result.groups[g].sub_scores[CAND_COUNT] =
          static_cast<double>(result.groups.size() - first_group);
      }
    }
  }

  /// Fit the discriminant over the retained groups and assign q-values.
  ///
  /// Split out of `finish()` so it can be run more than once over the SAME
  /// groups. That is what makes iteration cheap: the candidate picker is
  /// Refit the discriminant over the groups already scored.
  ///
  /// This used to recompute RT_DELTA first -- the ONLY sub-score that depended
  /// on the fitted retention-time map, which is what made "refit the map, then
  /// rescore" cheap. RT_DELTA is gone (see the note where it was declared), so
  /// NO sub-score depends on the map any more and re-running this after a new
  /// map produces a bit-identical result. `refitsChangeScores()` says so, and
  /// the retention-time refinement loop asks before spending a round.
  ///
  /// The map still matters upstream: it centres pass 2's extraction window.
  /// That has already run by the time anything calls this.
  void PeakGroupScorer::refit(const Library& library, Result& result,
                              const Options& options)
  {
    if (result.groups.empty()) { return; }
    // Counters are recomputed by the fit; reset so they do not accumulate
    // across iterations.
    result.identified_at_1pct = 0;
    fitAndAssign_(library, result, options);
  }

  /// Whether refitting after a new retention-time map can change any score.
  ///
  /// False since RT_DELTA was removed. Kept as a function rather than deleted
  /// at the call sites so that adding a map-dependent sub-score re-enables the
  /// loop by flipping one return, instead of by remembering that a loop was
  /// deleted somewhere.
  bool PeakGroupScorer::refitsChangeScores() { return false; }

  void PeakGroupScorer::fitAndAssign_(const Library& library, Result& result,
                                      const Options& options)
  {
    const std::size_t n_precursors = library.precursorCount();
    // PART 4 step 1: "Columns that are constant or all-missing are dropped."
    //
    // BEFORE the feature matrix is built, not after. This block only zeroes
    // `result.groups`, so building `features` first meant the classifier
    // received the original constants and NaNs while the log announced they had
    // been dropped -- and a later refit(), which rebuilds from the now-zeroed
    // groups, silently trained on different inputs from the first fit.
    //
    // Not implemented until now, and it bites today: USABLE_FRAGMENTS is a
    // constant 12.000 (a trace is degenerate only if constant, and in a window
    // where 45-60% of points are non-zero none ever is) and IM_DELTA is
    // all-NaN (Options::observed_im was added as a plumbing point and never
    // connected). A constant column makes the within-class covariance singular,
    // which the ridge then papers over silently; for the tree it is a wasted
    // split candidate. Either way the discriminant is asked to learn from a
    // column carrying no information.
    {
      std::vector<char> useful(N_SUB_SCORES, 0);
      std::vector<double> first(N_SUB_SCORES, 0.0);
      std::vector<char> seen(N_SUB_SCORES, 0);
      for (const auto& g : result.groups)
      {
        for (std::size_t j = 0; j < N_SUB_SCORES && j < g.sub_scores.size(); ++j)
        {
          const double v = g.sub_scores[j];
          if (!std::isfinite(v)) { continue; }
          if (!seen[j]) { seen[j] = 1; first[j] = v; }
          else if (v != first[j]) { useful[j] = 1; }
        }
      }
      std::vector<std::size_t> dropped;
      for (std::size_t j = 0; j < N_SUB_SCORES; ++j)
      {
        if (!useful[j]) { dropped.push_back(j); }
      }
      if (!dropped.empty())
      {
        std::ostringstream d;
        d << "dropping " << dropped.size() << " sub-score(s) carrying no information:";
        for (const std::size_t j : dropped) { d << ' ' << subScoreNames()[j]; }
        std::fprintf(stderr, "%s\n", d.str().c_str());
        // Zeroed rather than removed: the column indices are a contract with
        // subScoreNames(), nonpositive_features and seed_mask, and renumbering
        // them here would silently misalign all three.
        for (auto& g : result.groups)
        {
          for (const std::size_t j : dropped) { g.sub_scores[j] = 0.0; }
        }
      }
    }

    std::vector<std::vector<double>> features;
    std::vector<int> labels;
    std::vector<long long> group;
    features.reserve(result.groups.size());
    labels.reserve(result.groups.size());
    group.reserve(result.groups.size());
    for (const auto& g : result.groups)
    {
      features.push_back(g.sub_scores);
      labels.push_back(g.decoy ? 0 : 1);
      group.push_back(static_cast<long long>(g.precursor));
    }

    Scoring::LDAParams params;
    // These three were registered, parsed, and stored in Options -- and never
    // copied here, so `-use_pi0`, `-train_fdr` and `-train_fdr_initial` silently
    // returned the defaults. `acd000f` added all four of this family and touched
    // nothing in this file; only `-classifier_iterations` was later connected
    // (`c4b1416`), and it was found because an arm came back byte-identical.
    // The consequence is not merely a dead knob: triage arms on record claim to
    // have tested `use_pi0` and a relaxed `train_fdr`, and did not.
    //
    // Wiring rather than deprecating, because the intent is unambiguous -- all
    // three carry detailed help text describing behaviour they never had -- and
    // because the defaults are unchanged, so only a caller who asks for
    // something different sees any difference.
    if (options.train_fdr_initial > 0.0) { params.train_fdr_initial = options.train_fdr_initial; }
    if (options.train_fdr > 0.0) { params.train_fdr = options.train_fdr; }
    params.use_pi0 = options.use_pi0;
    // The classifier's stability knobs, exposed because it turned out to need
    // them. Adding a column of PURE NOISE moves identifications by up to 11.9%
    // at the operating point, and the churn is in the DISCRIMINANT rather than
    // the threshold: Spearman rho between a run and the same run with a noise
    // column is 0.79-0.82, and agreement in the top 500 is 0.38-0.42. A model
    // with 120 trees of depth 4 over ~3,400 positives has enough freedom to
    // find a different-but-equally-good solution whenever the feature set
    // changes, which makes every feature-level measurement on this pipeline
    // unreadable. These let that be tested rather than argued about.
    if (options.gbt_max_depth > 0) { params.gbt.max_depth = options.gbt_max_depth; }
    if (options.gbt_min_child_rows > 0)
    { params.gbt.min_child_rows = options.gbt_min_child_rows; }
    if (options.gbt_lambda > 0.0) { params.gbt.lambda = options.gbt_lambda; }
    if (options.gbt_max_delta_step > 0.0) { params.gbt.max_delta_step = options.gbt_max_delta_step; }
    if (options.gbt_intercept_zero) { params.gbt.intercept_zero = true; }
    if (options.gbt_clip_gain) { params.gbt.clip_gain = true; }
    if (options.gbt_warmup_rounds > 0) { params.gbt.warmup_rounds = options.gbt_warmup_rounds; }
    params.seed = static_cast<unsigned>(options.classifier_seed);
    params.gbt.fixed_bins = options.gbt_fixed_bins;
    // More SHALLOW trees rather than fewer deep ones. Depth 2 is the only
    // configuration measured to attribute an added column correctly, and it
    // costs about 11% of identifications against depth 4's centre. Boosting
    // recovers capacity additively -- many depth-2 trees approximate a richer
    // function without ever fitting a leaf of a handful of rows, which is the
    // surface gbt.h:441 names and the one that lets a noise column win a split
    // at depth 4.
    if (options.gbt_n_trees > 0) { params.gbt.n_trees = options.gbt_n_trees; }
    if (options.gbt_learning_rate > 0.0)
    { params.gbt.learning_rate = options.gbt_learning_rate; }
    // `-classifier_iterations` was a DEAD OPTION: registered, parsed, stored on
    // Options, and never read here, so the loop always ran LDAParams' own
    // n_iter = 3 whatever the flag said. Found because a frozen-trajectory arm
    // at `-classifier_iterations 1` came back BYTE-IDENTICAL to the default --
    // which is the only reason a dead knob ever gets noticed.
    //
    // It matters beyond tidiness. Re-selecting the positive set each iteration
    // is the path by which an added column changes the labels and therefore
    // everything downstream; n_iter = 1 is the cheapest way to cut it, and
    // until now that experiment could not be run at all.
    // >= 0, not > 0: zero is a real request (seed-only, no semi-supervised iterations) that the
    // old guard made unreachable; -1 is the "leave the LDAParams default" sentinel now.
    if (options.classifier_iterations >= 0)
    { params.n_iter = options.classifier_iterations; }
    // Mechanism 5 rides the same pattern: the stop rule and its threshold existed in LDAParams
    // from the start but had no CLI path, so `-classifier_iterations` was the cap AND the rule.
    // With the stop armed, the cap becomes the backstop and composition decides.
    params.stop_on_composition = options.classifier_stop_on_composition;
    if (options.classifier_stop_jaccard > 0.0)
    { params.stop_jaccard = options.classifier_stop_jaccard; }
    if (options.classifier_stop_shrink_floor > 0.0)
    { params.stop_shrink_floor = options.classifier_stop_shrink_floor; }
    if (options.classifier_stop_patience > 0)
    { params.stop_patience = options.classifier_stop_patience; }
    params.iteration_log = options.classifier_iteration_log;
    // Two of the sub-scores are lower-is-better by construction, so their
    // weights may never come out positive. XCORR_COELUTION is the mean |lag|
    // between fragment maxima -- a peak group IS a co-elution, so more lag is
    // less evidence -- and LIBRARY_RMSD is a deviation from the library
    // spectrum. Without the constraint the discriminant can fit, in-sample,
    // that being further from the expectation argues FOR a peptide.
    //
    // This is DIA-NN's check_weights (diann.cpp:6592) applied to the features
    // we actually have. Note what we do NOT have: it clips pdRT and pAcc, and
    // this scorer carries neither an rt_delta nor a mass_error sub-score --
    // both deliberately absent, see the SubScore comments. So the constraint
    // transfers in principle and covers different columns.
    // LOWER-IS-BETTER features, whose weights may never come out positive.
    //
    // IM_DELTA is stored as +|error|, and nothing stopped the semi-supervised
    // fit assigning it a POSITIVE weight -- i.e. rewarding candidates for
    // sitting FURTHER from their predicted mobility. That is not a theoretical
    // risk: the initial target class is heavily contaminated on a sparse
    // library (see the seeding note below), so a positive coefficient can be
    // learned whenever the contaminating set happens to carry larger errors
    // than the decoys.
    //
    // RT_DELTA was pinned here too, and the pin turned out to be the wrong
    // half of the problem: measured, true peaks sit FURTHER from the predicted
    // retention time than decoys do (AUC 0.215), so forcing "further is worse"
    // made a useless feature into a harmful one. The feature is now removed
    // rather than re-signed -- a group-level retention-time delta cannot
    // separate a real group from an interference group at all.
    params.nonpositive_features = {XCORR_COELUTION, LIBRARY_RMSD, IM_DELTA};
    params.match_decoy_candidate_counts = options.match_decoy_candidate_counts;
    params.fold_pool_rank = options.fold_pool_rank;

    // Seed the semi-supervised loop on CORR_SUM alone.
    //
    // The loop ignites by ranking with a single feature and taking whatever
    // clears q <= train_fdr_initial as its first positive set. lda.h says what
    // happens when that feature is too weak: "iteration 0 selects positives
    // with this w, finds none at q<=train_fdr_initial, skips the fit, so w is
    // unchanged and every later iteration skips too."
    //
    // That is not hypothetical here. On a REALISTIC library -- a stride sample
    // of the human proteome, where only ~1.5% of targets are present in the
    // sample -- ODIA identified NOTHING from 1,639,188 peak groups while DIA-NN
    // found 738 of the same 50,000 precursors. The automatic seed picks the
    // feature with the largest |t|, and against a positive class that is 98.5%
    // absent peptides, no feature's class means separate enough to ignite.
    //
    // CORR_SUM is the detector's own summed pairwise fragment correlation --
    // the quantity that moved rank-1 accuracy from 40.7% to 75.7%. DIA-NN
    // mandates the same choice rather than deriving it: reset_weights sets
    // w = (1, 0, 0, ...) so iteration 0 ranks by pTimeCorr, its co-elution sum,
    // and nothing else (diann.cpp:6584).
    //
    // Only when the picker actually computed it. The amplitude detector leaves
    // CORR_SUM at 0, and seeding on a constant would guarantee the failure this
    // is meant to prevent.
    //
    // The OpenSWATH picker leaves it at 0 too, and the guard did not say so --
    // `coelution_picking` stays true under `-picker openswath` because it is
    // driven by a different flag, so the pure-OpenSWATH arm seeded its
    // semi-supervised loop on an identically-zero column. That is not a
    // hypothetical: every recorded number for that arm was produced this way,
    // which means the measured verdict on OpenSWATH-style picking (+4.7% on
    // IH1, -38.8% on Astral) is confounded with this defect and cannot be read
    // as a property of the picker until it is re-measured.
    //
    // `union_openswath` is deliberately NOT excluded: it mixes co-elution
    // candidates that carry a real corr_sum with OpenSWATH ones that carry
    // zero, so the column is informative rather than constant.
    const bool pure_openswath = options.openswath_picking && !options.union_picking;
    if (options.coelution_picking && !pure_openswath)
    {
      params.seed_mask.assign(N_SUB_SCORES, 0);
      params.seed_mask[CORR_SUM] = 1;
    }
    // Threading is per-classifier, not on LDAParams: the LDA solve is a small
    // dense Cholesky and does not want threads, while the tree and network
    // fits do.
    if (options.classifier == "xgboost")
    {
      // Same learner as "gbt" -- the algorithm here IS XGBoost's, second-order
      // Newton leaf values and the same split gain -- configured as pyProphet
      // 3.0.15 configures XGBoost 3.2.0. Six of its nine parameters already
      // match; this changes max_depth 4 -> 6 and eta 0.1 -> 0.3.
      params.classifier = Scoring::LDAParams::Classifier::GBT;
      // WHOLESALE replacement, so it must not silently discard the per-parameter
      // overrides applied above -- `-classifier xgboost -gbt_max_depth 8` used to
      // ignore the depth entirely. Take the preset first, then re-apply anything
      // the caller set explicitly.
      params.gbt = Scoring::GBTParams::pyprophet();
      if (options.gbt_max_depth > 0) { params.gbt.max_depth = options.gbt_max_depth; }
      if (options.gbt_min_child_rows > 0) { params.gbt.min_child_rows = options.gbt_min_child_rows; }
      if (options.gbt_lambda > 0.0) { params.gbt.lambda = options.gbt_lambda; }
      if (options.gbt_max_delta_step > 0.0) { params.gbt.max_delta_step = options.gbt_max_delta_step; }
      if (options.gbt_intercept_zero) { params.gbt.intercept_zero = true; }
      if (options.gbt_clip_gain) { params.gbt.clip_gain = true; }
      if (options.gbt_warmup_rounds > 0) { params.gbt.warmup_rounds = options.gbt_warmup_rounds; }
      params.seed = static_cast<unsigned>(options.classifier_seed);
      if (options.gbt_n_trees > 0) { params.gbt.n_trees = options.gbt_n_trees; }
      if (options.gbt_learning_rate > 0.0) { params.gbt.learning_rate = options.gbt_learning_rate; }
      params.gbt.fixed_bins = options.gbt_fixed_bins;
      params.gbt.n_threads = static_cast<int>(options.threads);
    }
    else if (options.classifier == "gbt")
    {
      params.classifier = Scoring::LDAParams::Classifier::GBT;
      params.gbt.n_threads = static_cast<int>(options.threads);
    }
    else if (options.classifier == "nn")
    {
      params.classifier = Scoring::LDAParams::Classifier::NN;
      params.nn.n_threads = static_cast<int>(options.threads);
    }

    // THE SCORING-ENGINE SEAM.
    //
    // Every engine takes the same (features, labels, group) and returns the
    // same ScoredGroups, so swapping one is a branch here and nothing
    // downstream can tell them apart. A common function signature is the whole
    // interface; a class hierarchy would add ceremony and make A/B harder, not
    // easier.
    std::string engine_note;
    // The sub-score names travel with a saved gbt model so a frozen application
    // refuses a run whose sub-score set differs, whatever its width.
    const std::vector<std::string> frozen_names = subScoreNames();
    const auto scored =
      options.classifier == "percolator"
        ? Scoring::scorePercolator(features, labels, group, params,
                                   subScoreNames(), &engine_note,
                                   options.classifier_model_out,
                                   options.classifier_model_in)
        // The gbt engine takes the same two paths as percolator: train-and-save,
        // or load-and-apply with no training. Measured 2026-09-05: the native
        // retraining flips between a compact and a saturated score regime
        // under small changes to the binary or the feature distribution
        // (pick/wf_hist.txt), while the same evidence scored by ONE model held
        // fixed across arms does not (pick/wf_fixedmodel.txt). A frozen model
        // is how two arms are compared on evidence rather than on retraining.
        : Scoring::scoreSemiSupervisedLDA(features, labels, group, params,
                                          options.classifier_model_out,
                                          options.classifier_model_in,
                                          &frozen_names);
    // Reported through the same channel the picker census uses, so an engine
    // swap is visible in the run log rather than only in the numbers.
    // Same channel the picker census uses, so an engine swap is visible in the
    // run log rather than only in the resulting numbers.
    if (!engine_note.empty()) { std::fprintf(stderr, "%s\n", engine_note.c_str()); }
    result.iterations_trained = scored.n_iterations_trained;
    result.iterations_skipped = scored.n_iterations_skipped;
    // A fit that never trained calibrates its q-values against an
    // initialisation rather than a model, so those are not an FDR either.
    result.fdr_valid = scored.n_iterations_trained > 0;

    for (std::size_t i = 0; i < result.groups.size(); ++i)
    {
      result.groups[i].dscore = scored.dscore[i];
      result.groups[i].qvalue = scored.qvalue[i];
      result.groups[i].pep = scored.pep[i];
    }

    // Counted on the best group per precursor: FDR is a per-precursor quantity,
    // and counting rows would report a precursor several times over.
    std::vector<double> best(n_precursors, 1.0);
    std::vector<char> is_target(n_precursors, 0);
    for (const auto& g : result.groups)
    {
      if (g.qvalue < best[g.precursor]) { best[g.precursor] = g.qvalue; }
      if (!g.decoy) { is_target[g.precursor] = 1; }
    }
    if (result.fdr_valid)
    {
      for (std::size_t i = 0; i < n_precursors; ++i)
      {
        if (is_target[i] && best[i] <= 0.01) { ++result.identified_at_1pct; }
      }
    }  }

  PeakGroupScorer::Result PeakGroupScorer::Session::finish()
  {
    Result& result = result_;
    const Options& options = options_;
    const std::size_t n_precursors = library_->precursorCount();

    // PEAK_WIDTH_RATIO was stored as a raw cycle count per candidate, because a
    // ratio needs the run's median and no single precursor knows it. Normalise
    // now: a peptide elutes on the chromatography's timescale, so what carries
    // information is being wide or narrow RELATIVE to this run, not absolutely.
    {
      std::vector<double> widths;
      widths.reserve(result.groups.size());
      for (const auto& g : result.groups)
      {
        const double w = g.sub_scores[PEAK_WIDTH_RATIO];
        if (std::isfinite(w) && w > 0.0) { widths.push_back(w); }
      }
      if (!widths.empty())
      {
        std::nth_element(widths.begin(), widths.begin() + widths.size() / 2, widths.end());
        const double median = widths[widths.size() / 2];
        if (median > 0.0)
        {
          for (auto& g : result.groups) { g.sub_scores[PEAK_WIDTH_RATIO] /= median; }
        }
      }
    }

    // MASS_ACCURACY, for the same reason and by the same route: a deviation of
    // -3 ppm means nothing until you know where this run sits.
    //
    // Centred on the median over ALL candidates, not over the accepted ones.
    // That is deliberate and it is this morning's lesson: a statistic pooled
    // over q<=0.01 groups describes the precursors that were already easy, and
    // sizing anything from it cost half the run when I tried it for the m/z
    // window. Here the centre only has to be robust, and the median over every
    // candidate -- most of which are wrong -- is a fine estimator of the
    // INSTRUMENT's offset precisely because it does not care which are right.
    {
      std::vector<double> dev;
      dev.reserve(result.groups.size());
      for (const auto& g : result.groups)
      {
        const double d = g.sub_scores[MASS_ACCURACY];
        if (std::isfinite(d)) { dev.push_back(d); }
      }
      if (!dev.empty())
      {
        std::nth_element(dev.begin(), dev.begin() + dev.size() / 2, dev.end());
        // v1.16: a frozen centre (Options::mass_accuracy_centre, finite) replaces the median so that an
        // intervention on a few candidates (a transition mask) cannot move every row's MASS_ACCURACY
        // through this run-wide statistic. Off (NaN) = the native median, identical arithmetic.
        const double centre = std::isfinite(options.mass_accuracy_centre) ? options.mass_accuracy_centre
                                                                         : dev[dev.size() / 2];
        for (auto& g : result.groups)
        {
          double& s = g.sub_scores[MASS_ACCURACY];
          // Negated so larger is better, like every other column.
          if (std::isfinite(s)) { s = -std::abs(s - centre); }
        }
        std::fprintf(stderr,
                     "fragment mass accuracy as a sub-score: %zu of %zu candidates "
                     "carried a deviation, centred on %.3f ppm (exact %.17g%s)\n",
                     dev.size(), result.groups.size(), centre, centre,
                     std::isfinite(options.mass_accuracy_centre) ? ", FROZEN by -mass_accuracy_centre" : "");
      }
    }

    // Ablation, last, so it survives every normalisation above rather than
    // being recomputed by one of them.
    for (const int idx : options.disabled_sub_scores)
    {
      if (idx < 0 || idx >= static_cast<int>(N_SUB_SCORES)) { continue; }
      for (auto& g : result.groups) { g.sub_scores[idx] = 0.0; }
      std::fprintf(stderr, "ablated sub-score %s (flattened to a constant)\n",
                   subScoreNames()[static_cast<std::size_t>(idx)].c_str());
    }

    // Peak groups arrive in whatever order the extractor finished their
    // precursors, which is retention-time order rather than library order. Put
    // them back in library order before anything downstream sees them: the rows
    // of the feature matrix are these groups, and a classifier fitted on rows
    // ordered by elution time is a classifier that can depend on it. Stable, so
    // a precursor's candidates keep the order the search produced.
    // Canonical order: (precursor, apex RT, apex intensity). NOT "precursor,
    // then whatever order the search produced".
    //
    // PART 4 step 0 of the OpenSWATH handoff requires (group, apex RT, feature
    // id) before anything else, and calls skipping it the determinism bug of
    // section 3.3 -- which that project rates its most important non-result.
    // Sorting by precursor alone leaves a precursor's candidates in discovery
    // order, so a thread-count change reorders rows, which reassigns folds,
    // which changes the fit. Two runs of the same input would then disagree for
    // a reason invisible in any output.
    //
    // apex_intensity substitutes for "feature id": we have no stable per-feature
    // identifier, and it breaks ties that apex RT alone leaves.
    //
    // Sorted through an index PERMUTATION rather than in place, because
    // `MassAnchor::group` indexes into this vector. Sorting the groups directly
    // left every harvested anchor pointing at whatever group moved into its old
    // slot -- so an anchor from an accepted target could be read with a decoy's
    // flag and a rejected group's q-value. Silent, and it invalidated a whole
    // measurement before it was caught.
    std::vector<std::uint32_t> order(result.groups.size());
    std::iota(order.begin(), order.end(), 0u);
    std::stable_sort(order.begin(), order.end(),
                     [&](std::uint32_t ia, std::uint32_t ib) {
                       const PeakGroup& a = result.groups[ia];
                       const PeakGroup& b = result.groups[ib];
                       if (a.precursor != b.precursor) { return a.precursor < b.precursor; }
                       if (a.apex_rt != b.apex_rt) { return a.apex_rt < b.apex_rt; }
                       return a.apex_intensity < b.apex_intensity; });

    if (!result.mass_anchors.empty())
    {
      std::vector<std::uint32_t> moved_to(order.size());
      for (std::size_t n = 0; n < order.size(); ++n) { moved_to[order[n]] = static_cast<std::uint32_t>(n); }
      for (MassAnchor& a : result.mass_anchors)
      { if (a.group < moved_to.size()) { a.group = moved_to[a.group]; } }
    }

    {
      std::vector<PeakGroup> reordered;
      reordered.reserve(result.groups.size());
      for (const std::uint32_t o : order) { reordered.push_back(std::move(result.groups[o])); }
      result.groups.swap(reordered);
    }

    // The fragvec rows are parallel to `groups`, so they take the SAME
    // permutation. This is the whole reason the 78 floats are retained rather
    // than streamed as each candidate finishes: the row order -out publishes --
    // and therefore the ordinal the export is keyed on -- does not exist until
    // the sort above has run. Costs one transient copy of 312 B x groups
    // (6.8 GB at 21.8 M rows, against the arm's measured 716 GB peak), the same
    // shape of transient the group reorder just above pays.
    if (!result.fragvec.empty())
    {
      std::vector<float> fv_reordered;
      fv_reordered.reserve(result.fragvec.size());
      for (const std::uint32_t o : order)
      {
        const float* src = result.fragvec.data() + static_cast<std::size_t>(o) * N_FRAGVEC;
        fv_reordered.insert(fv_reordered.end(), src, src + N_FRAGVEC);
      }
      result.fragvec.swap(fv_reordered);
    }

    for (const auto& g : result.groups)
    {
      (g.decoy ? result.decoy_groups : result.target_groups) += 1;
    }
    if (result.groups.empty()) { return result; }

    // No decoys means no negative class. The scorer would still return numbers
    // -- q = 0 for everything -- and they would be read as an FDR. Refuse
    // instead, leaving q at 1 so nothing downstream mistakes silence for
    // confidence.
    if (result.decoy_groups == 0 || result.target_groups == 0)
    {
      return result;
    }

    // Which criterion actually rejected. Reported unconditionally: the
    // aggregate "N precursors yielded no candidate peak group" was true and
    // useless, and on a realistic library it is 19,150 of 100,000.
    {
      const PickerRejects r = totalRejects();
      std::ostringstream w;
      auto pc = [](std::size_t t, std::size_t d) {
        std::ostringstream o; o.setf(std::ios::fixed); o.precision(2);
        o << t << "/" << d;
        if (t > 0) { o << " (" << (double(d) / double(t)) << "x)"; }
        return o.str();
      };
      // target/decoy at EVERY stage, because the imbalance has to arise
      // somewhere and only a per-stage split says where.
      w << "picker rejections, target/decoy (decoy:target ratio):"
        << "\n  no points (<3)      " << pc(r.no_points[0], r.no_points[1])
        << "\n  empty trace         " << pc(r.empty_trace[0], r.empty_trace[1])
        << "\n    of which gate C   " << pc(r.gate_c[0], r.gate_c[1])
        << "\n    too few excursions" << pc(r.few_excursions[0], r.few_excursions[1])
        << "\n    truly zero trace  " << pc(r.zero_trace[0], r.zero_trace[1])
        << "\n  precursors reached  " << pc(r.reached[0], r.reached[1])
        << "\n  scan positions      " << pc(r.scans[0], r.scans[1])
        << "\n  <2 fragments        " << pc(r.too_few_present[0], r.too_few_present[1])
        << "\n  <2 transitions      " << pc(r.too_few_transitions[0], r.too_few_transitions[1])
        << "\n  below min_corr      " << pc(r.below_corr[0], r.below_corr[1])
        << "\n  reference zero      " << pc(r.reference_zero[0], r.reference_zero[1])
        << "\n  not a local max     " << pc(r.not_local_max[0], r.not_local_max[1])
        << "\n  below apex_evidence " << pc(r.below_apex_evidence[0], r.below_apex_evidence[1])
        // THE one to watch: max_corr_diff is a margin around each precursor's
        // OWN best correlation, so an absent precursor -- whose best IS noise --
        // keeps far more positions than a present one whose best is a real
        // peak. Decoys are absent by construction and targets are a mixture, so
        // this stage is where the classes should diverge if that is the cause.
        << "\n  outside max_corr_diff " << pc(r.outside_margin[0], r.outside_margin[1])
        // These print target/decoy like every line above. They used to print the
        // ARRAY, which after the split-by-class change is a pointer -- the log
        // read "0x7fff6c83b1b8 with <3 points". A diagnostic that silently
        // prints an address is worse than one that prints nothing, because it
        // still looks like a reading.
        // Not split by class, so it prints as one number (and on its own line --
        // it previously ran on from the line above with no separator).
        << "\n  min_fragments_at_apex " << r.too_few_at_apex
        << "\n  transition_mask candidates " << r.masked_candidates << " fragments removed " << r.masked_fragments
        << "\n  EXTRACTION losses (no usable chromatogram): "
        << pc(r.no_points[0], r.no_points[1]) << " with <3 points, "
        << pc(r.empty_trace[0], r.empty_trace[1]) << " with an all-zero trace"
        << "\n  precursors that never reached the correlation loop: "
        << pc(r.too_few_transitions[0], r.too_few_transitions[1]) << " with <2 transitions, "
        << r.too_few_cycles << " with too few cycles; "
        << r.no_hit_anywhere << " entered it and found no qualifying position";

      // Why MS1_COELUTION arrives constant. It is the only sub-score measured to
      // discriminate (13.7x top bin against 1.0x saturated for every presence
      // statistic), so it being dropped as uninformative needs a reason on the
      // record rather than a plausible story.
      {
        const std::size_t spans = ms1_spans.load(std::memory_order_relaxed);
        w << "\n  MS1_COELUTION: " << ms1_ok.load(std::memory_order_relaxed) << " finite, "
          << ms1_unavailable.load(std::memory_order_relaxed) << " no MS1 available, "
          << ms1_no_signal.load(std::memory_order_relaxed) << " all-zero in the window, "
          << ms1_too_short.load(std::memory_order_relaxed) << " under 5 cycles, "
          << ms1_flat.load(std::memory_order_relaxed) << " flat (zero variance)"
          << "\n    unavailable breakdown: "
          << ms1_null.load(std::memory_order_relaxed) << " null ptr, "
          << ms1_empty.load(std::memory_order_relaxed) << " empty, "
          << ms1_index_oob.load(std::memory_order_relaxed) << " index past matrix, "
          << ms1_degenerate_span.load(std::memory_order_relaxed) << " hi<=lo";
        if (spans > 0)
        {
          w << "; mean distinct MS1 bins per candidate "
            << (double(ms1_bins_spanned.load(std::memory_order_relaxed)) / double(spans));
        }
      }
      std::fprintf(stderr, "%s\n", w.str().c_str());
    }

    fitAndAssign_(*library_, result, options);

    return result;
  }

  PeakGroupScorer::Result PeakGroupScorer::score(const Library& library,
                                                 const Chromatograms& chromatograms,
                                                 const Options& options)
  {
    Session session(library, options);
    const auto& p = library.precursors();
    for (std::size_t i = 0; i < library.precursorCount(); ++i)
    {
      const std::uint32_t tb = p.transition_begin[i];
      const std::uint32_t tc = p.transition_count[i];
      // The streaming path records these inside add(); this one skips before
      // calling it, so the SAME precursor would come out NotReached here and
      // NoTransitions there. Record it at the skip so the two paths agree.
      if (tc == 0 || tb >= chromatograms.begin.size())
      {
        if (options.terminal_reason)
        {
          options.terminal_reason[i] =
            static_cast<std::uint8_t>(TerminalReason::NoTransitions);
        }
        continue;
      }

      // A window onto the flat array in the same shape the extractor hands a
      // streaming sink, so both paths run the same code on the same data.
      PrecursorChromatogram trace;
      trace.precursor = static_cast<std::uint32_t>(i);
      trace.transition_begin = tb;
      trace.transition_count = tc;
      trace.axis = chromatograms.precursor_axis[i];
      trace.axis_begin = chromatograms.precursor_axis_begin[i];
      trace.cycles = chromatograms.precursor_cycles[i];
      if (trace.axis < chromatograms.axes.size())
      {
        trace.rt = chromatograms.axes[trace.axis].data() + trace.axis_begin;
      }
      trace.points = chromatograms.intensity.data();
      trace.offset = chromatograms.begin.data() + tb;
      trace.count = chromatograms.count.data() + tb;
      session.add(trace);
    }
    return session.finish();
  }


  std::vector<MassResidual> PeakGroupScorer::acceptedMassResiduals(const Result& result,
                                                                   double q_threshold,
                                                                   bool include_decoys)
  {
    std::vector<MassResidual> out;
    out.reserve(result.mass_anchors.size() / 2 + 1);
    for (const MassAnchor& a : result.mass_anchors)
    {
      if (a.group >= result.groups.size()) { continue; }   // cannot happen; not worth trusting
      const PeakGroup& g = result.groups[a.group];
      if (g.decoy && !include_decoys) { continue; }
      if (!(g.qvalue <= q_threshold)) { continue; }
      MassResidual r = a.residual;
      r.decoy = g.decoy;          // the single source of truth
      out.push_back(r);
    }
    return out;
  }

  // ------------------------------------------------ Gate C hash calibration (D6)

  std::uint64_t PeakGroupScorer::gateIdentityHash(std::string_view modified_sequence, int charge,
                                                  bool decoy, std::uint64_t seed)
  {
    std::uint64_t h = 0xcbf29ce484222325ULL;               // FNV-1a 64 offset basis
    const auto mix = [&h](unsigned char b) { h ^= b; h *= 0x100000001b3ULL; };
    for (int k = 0; k < 8; ++k) { mix(static_cast<unsigned char>((seed >> (8 * k)) & 0xffU)); }
    for (const char c : modified_sequence) { mix(static_cast<unsigned char>(c)); }
    mix(0);
    mix(static_cast<unsigned char>(charge & 0xff));
    mix(decoy ? 1 : 0);
    // splitmix64 finaliser: FNV's low bits are weak for short keys, and the
    // selection orders by the whole word.
    h ^= h >> 30; h *= 0xbf58476d1ce4e5b9ULL;
    h ^= h >> 27; h *= 0x94d049bb133111ebULL;
    h ^= h >> 31;
    return h;
  }

  bool PeakGroupScorer::gateUsesHashCalibration(const Options& options)
  {
    return options.gate_alpha > 0.0 && options.gate_mode != "prominence" &&
           options.gate_calibration == "hash";
  }

  std::vector<std::uint32_t> PeakGroupScorer::selectGateCalibrationDecoys(
    const Library& library, const std::vector<std::uint32_t>& plan_points,
    const Options& options, std::size_t* eligible_out)
  {
    const auto& p = library.precursors();
    const std::size_t n = std::min(plan_points.size(), library.precursorCount());
    std::vector<std::pair<std::uint64_t, std::uint32_t>> keyed;
    for (std::size_t i = 0; i < n; ++i)
    {
      // Exactly the decoys Session::add brings to the gate in this pass:
      // transitions present, pointCount(0) >= 3, not bypassed by the oracle.
      if (p.decoy[i] == 0 || p.transition_count[i] == 0 || plan_points[i] < 3) { continue; }
      if (options.oracle_rt != nullptr && std::isfinite(options.oracle_rt[i])) { continue; }
      keyed.emplace_back(gateIdentityHash(library.strings().get(p.modified_sequence[i]),
                                          static_cast<int>(p.charge[i]), true,
                                          options.gate_calibration_seed),
                         static_cast<std::uint32_t>(i));
    }
    if (eligible_out != nullptr) { *eligible_out = keyed.size(); }
    const std::size_t want = options.gate_calibration_n;
    if (want == 0 || keyed.size() < want) { return {}; }
    // (hash, index) is a strict total order: the selection is a set, whatever
    // order the library or the extraction presents precursors in.
    std::nth_element(keyed.begin(), keyed.begin() + std::ptrdiff_t(want - 1), keyed.end());
    std::vector<std::uint32_t> out;
    out.reserve(want);
    const auto pivot = keyed[want - 1];
    for (const auto& k : keyed) { if (k <= pivot) { out.push_back(k.second); } }
    std::sort(out.begin(), out.end());
    return out;
  }

  namespace
  {
    /// Measures Gate C's statistic for the selected decoys, exactly as
    /// Session::add computes it (same points, same smoothing, same RT centre),
    /// and nothing else. Every other precursor arrives empty and is ignored.
    class GateCalibrationSink final : public ChromatogramSink
    {
    public:
      GateCalibrationSink(const PeakGroupScorer::Options& options, const std::vector<char>& selected)
        : options_(options), selected_(selected) {}

      void accept(const PrecursorChromatogram& c) override
      {
        const std::size_t i = c.precursor;
        if (i >= selected_.size() || selected_[i] == 0) { return; }
        if (c.transition_count == 0) { return; }
        const std::size_t points = c.pointCount(0);
        if (points < 3) { return; }
        if (options_.oracle_rt != nullptr && std::isfinite(options_.oracle_rt[i])) { return; }
        std::size_t contributing = 0;
        const double m = coelutionEvidence(c, points, options_.gate_smooth_half, &contributing);
        const double rt = c.rt != nullptr
          ? 0.5 * (double(c.rt[0]) + double(c.rt[points - 1]))
          : std::numeric_limits<double>::quiet_NaN();
        stats.push_back(m);
        rts.push_back(rt);
      }

      std::vector<double> stats;
      std::vector<double> rts;

    private:
      const PeakGroupScorer::Options& options_;
      const std::vector<char>& selected_;
    };
  }

  PeakGroupScorer::GateCalibrationReport PeakGroupScorer::calibrateGateByHash(
    const Library& library, SpectrumSource& source,
    const ChromatogramExtractor::Options& extract, Sink& sink)
  {
    const auto t0 = std::chrono::steady_clock::now();
    const Options& so = sink.options();
    GateCalibrationReport r;
    r.requested = so.gate_calibration_n;

    // 1. The pass's plan: which precursors it assigns, over how many cycles.
    std::vector<std::uint32_t> plan;
    {
      ChromatogramExtractor::Options po = extract;
      po.plan_points = &plan;
      po.plan_only = true;
      po.terminal_reason = nullptr;
      po.progress_every = 0;
      NullChromatogramSink none;
      ChromatogramExtractor::extract(library, source, po, none, nullptr);
    }

    // 2. The sample, chosen before a peak is read.
    const std::vector<std::uint32_t> chosen =
      selectGateCalibrationDecoys(library, plan, so, &r.eligible);
    r.selected = chosen.size();
    if (chosen.empty())
    {
      // Insufficient sample: as in arrival mode, a null that cannot reach
      // gate_calibration_n never arms and the gate admits everything.
      r.armed = false;
      r.tau = -std::numeric_limits<double>::infinity();
      sink.freezeGate(r.tau);
      r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      return r;
    }

    // 3. Extract ONLY them, with the production pass's options otherwise.
    std::vector<char> keep(library.precursorCount(), 0);
    for (const std::uint32_t i : chosen) { keep[i] = 1; }
    if (extract.precursor_keep != nullptr)
    {
      // The plan already honoured the caller's mask; AND it in anyway so the
      // pre-pass can never extract something the pass excludes.
      for (const std::uint32_t i : chosen)
      { if (extract.precursor_keep[i] == 0) { keep[i] = 0; } }
    }
    GateCalibrationSink cal(so, keep);
    ChromatogramExtractor::Stats cst;
    {
      ChromatogramExtractor::Options co = extract;
      co.precursor_keep = keep.data();
      co.terminal_reason = nullptr;
      co.plan_points = nullptr;
      co.plan_only = false;
      ChromatogramExtractor::extract(library, source, co, cal, &cst);
    }
    r.chunks = cst.chunks;
    r.spectra_decoded = cst.spectra_decoded;
    r.measured = cal.stats.size();
    if (r.measured != r.selected)
    {
      throw std::runtime_error(
        "Gate C hash calibration: " + std::to_string(r.selected) +
        " decoys were selected from the pass plan but " + std::to_string(r.measured) +
        " reached the gate statistic -- the plan and the hand-over disagree");
    }

    // 4. tau exactly as arrival mode computes it from its sample.
    std::vector<double> v = cal.stats;
    std::sort(v.begin(), v.end());
    const std::size_t k = std::min(v.size() - 1,
      std::size_t((1.0 - so.gate_alpha) * double(v.size())));
    r.tau = v[k];
    r.zeros = std::size_t(std::count(v.begin(), v.end(), 0.0));
    std::vector<double> rt;
    for (const double x : cal.rts) { if (std::isfinite(x)) { rt.push_back(x); } }
    std::sort(rt.begin(), rt.end());
    const double qs[7] = {0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0};
    for (int j = 0; j < 7; ++j)
    {
      r.rt_q[j] = rt.empty() ? std::numeric_limits<double>::quiet_NaN()
                             : rt[std::min(rt.size() - 1, std::size_t(qs[j] * double(rt.size() - 1) + 0.5))];
    }
    r.armed = true;
    sink.freezeGate(r.tau);
    r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return r;
  }

  std::string PeakGroupScorer::GateCalibrationReport::describe() const
  {
    char b[768];
    if (!armed)
    {
      std::snprintf(b, sizeof b,
                    "gate C NOT armed (hash): %zu eligible decoys < -gate_calibration_n %zu; "
                    "the gate admits every precursor this pass (arrival mode would not arm "
                    "either), %.1f s",
                    eligible, requested, seconds);
      return b;
    }
    std::snprintf(b, sizeof b,
                  "gate C null armed (hash): n=%zu of %zu eligible decoys, tau=%.6g, "
                  "zeros=%.4f (%zu), calibration RT p0/p5/p25/p50/p75/p95/p100 "
                  "%.1f/%.1f/%.1f/%.1f/%.1f/%.1f/%.1f s, pre-pass %zu chunk(s), "
                  "%zu spectra decoded, %.1f s; the frozen tau gates every precursor",
                  measured, eligible, tau, measured ? double(zeros) / double(measured) : 0.0,
                  zeros, rt_q[0], rt_q[1], rt_q[2], rt_q[3], rt_q[4], rt_q[5], rt_q[6],
                  chunks, spectra_decoded, seconds);
    return b;
  }

} // namespace ODIA
