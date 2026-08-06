// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/PeakGroupScorer.h>

#include <odia/scoring/gbt.h>
#include <odia/scoring/lda.h>
#include <odia/scoring/score.h>

#include <algorithm>
#include <cmath>
#include <numeric>

namespace ODIA
{

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
    std::vector<double> noiseNormalisedTrace(const PrecursorChromatogram& c,
                                             std::size_t points)
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
        if (!(mad > 0.0)) { continue; }
        const double scale = 1.0 / (1.4826 * mad);   // MAD -> sigma for normal noise
        for (std::uint32_t i = 0; i < n && i < points; ++i)
        {
          total[i] += (at[i] - median) * scale;
        }
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
    double localBackground(const std::vector<double>& trace,
                           std::size_t lo, std::size_t hi)
    {
      std::vector<double> outside;
      outside.reserve(trace.size());
      for (std::size_t i = 0; i < trace.size(); ++i)
      {
        if (i < lo || i > hi) { outside.push_back(trace[i]); }
      }
      if (outside.size() < 3)
      {
        outside.assign(trace.begin(), trace.end());
      }
      if (outside.empty()) { return 0.0; }
      std::sort(outside.begin(), outside.end());
      return outside[outside.size() / 2];
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

    struct Candidate
    {
      std::size_t apex = 0, left = 0, right = 0;
      double apex_value = 0.0;
    };

    /// Local maxima of the smoothed trace, strongest first, with boundaries
    /// descending to a fraction of the apex.
    std::vector<Candidate> findCandidates(const std::vector<double>& smoothed,
                                          std::size_t max_candidates,
                                          double boundary_fraction)
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
        const double floor_value = boundary_fraction * smoothed[i];

        // Scale-aware descent. A bare "stop on any uptick" ends the window on
        // the first noise wobble; "tolerate upticks" without a scale walks
        // across a valley into the neighbouring peak and integrates both. So:
        // stop at the floor, stop if the trace rebounds more than a fraction
        // of the apex above the running minimum, and never run further than a
        // bounded number of points.
        const double rebound_limit = 0.25 * smoothed[i];
        const std::size_t max_span = std::max<std::size_t>(4, n / 4);

        std::size_t l = i, guard = 0;
        double run_min = smoothed[i];
        while (l > 0 && guard++ < max_span)
        {
          const double v = smoothed[l - 1];
          if (v <= floor_value) { break; }
          if (v > run_min + rebound_limit) { break; }
          run_min = std::min(run_min, v);
          --l;
        }
        std::size_t r = plateau_end;
        guard = 0;
        run_min = smoothed[plateau_end];
        while (r + 1 < n && guard++ < max_span)
        {
          const double v = smoothed[r + 1];
          if (v <= floor_value) { break; }
          if (v > run_min + rebound_limit) { break; }
          run_min = std::min(run_min, v);
          ++r;
        }
        i = plateau_end;
        c.left = l;
        c.right = r;
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
    /// S08 against DIA-NN's confident set, our amplitude picker put the true
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
      const PrecursorChromatogram& c, std::size_t half_window,
      double min_corr_score, double max_corr_diff, double apex_evidence,
      std::size_t smooth_half_width, double boundary_fraction,
      std::size_t max_candidates)
    {
      std::vector<Candidate> found;
      const std::uint32_t tc = c.transition_count;
      const std::size_t n = c.cycles;
      const std::size_t S = std::max<std::size_t>(1, half_window);
      if (n < 2 * S + 4 || tc < 2) { return found; }

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
      for (std::size_t k = S + 1; k + S + 2 < n; ++k)
      {
        // Cheap rejects first: something must be here, and it must persist
        // across neighbouring cycles rather than being a single spike.
        std::size_t present = 0;
        for (std::uint32_t f = 0; f < tc; ++f)
        {
          if (tr[f][k - 1] > 0.0 || tr[f][k] > 0.0 || tr[f][k + 1] > 0.0) { ++present; }
        }
        if (present < 2) { continue; }

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
          if (score[ref] < min_corr_score) { break; }
          if (!(sm[ref][k] > 0.0)) { continue; }

          const std::size_t half = std::max<std::size_t>(S / 3, 1);
          bool is_max = true;
          for (std::size_t j = (k > half ? k - half : 0); j <= k + half && j < n; ++j)
          {
            if (sm[ref][j] > sm[ref][k]) { is_max = false; break; }
          }
          if (!is_max) { continue; }

          double best_near = 0.0;
          const std::size_t e = S > 1 ? S - 1 : 1;
          for (std::size_t j = (k > e ? k - e : 0); j <= k + e && j < n; ++j)
          {
            best_near = std::max(best_near, sm[ref][j]);
          }
          if (best_near > 0.0 && sm[ref][k] < apex_evidence * best_near) { continue; }

          hits.push_back({k, score[ref]});
          break;
        }
      }
      if (hits.empty()) { return found; }

      // Keep by MARGIN from the best, not by rank.
      double best = 0.0;
      for (const auto& h : hits) { best = std::max(best, h.corr_sum); }
      std::stable_sort(hits.begin(), hits.end(),
                       [](const Hit& x, const Hit& y) { return x.corr_sum > y.corr_sum; });

      // Boundaries from the summed trace, as before: the extent of a peak is
      // not what changed here, only which positions are peaks.
      std::vector<double> total(n, 0.0);
      for (std::uint32_t f = 0; f < tc; ++f)
      {
        for (std::size_t j = 0; j < n; ++j) { total[j] += tr[f][j]; }
      }
      for (const auto& h : hits)
      {
        if (h.corr_sum < best - max_corr_diff) { break; }
        if (found.size() >= max_candidates) { break; }
        Candidate cd;
        cd.apex = h.apex;
        cd.apex_value = total[h.apex];
        const double floor_value = boundary_fraction * total[h.apex];
        std::size_t l = h.apex, r = h.apex;
        while (l > 0 && total[l - 1] > floor_value) { --l; }
        while (r + 1 < n && total[r + 1] > floor_value) { ++r; }
        cd.left = l; cd.right = r;
        // One candidate per apex; DIA-NN accepts at most one per scan position
        // and we must not emit two peaks that share one.
        bool dup = false;
        for (const auto& g : found) { if (g.apex == cd.apex) { dup = true; break; } }
        if (!dup) { found.push_back(cd); }
      }
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

  const std::vector<std::string>& PeakGroupScorer::subScoreNames()
  {
    static const std::vector<std::string> names{
      "var_xcorr_shape", "var_xcorr_coelution", "var_library_corr",
      "var_library_dotprod", "var_intensity_score", "var_log_sn",
      "var_usable_fragments", "var_library_rmsd", "var_yseries_score",
      "var_fragment_coverage"};
    return names;
  }

  PeakGroupScorer::Session::Session(const Library& library, const Options& options)
    : library_(&library), options_(options)
  {
  }

  void PeakGroupScorer::Session::add(const PrecursorChromatogram& chromatogram)
  {
    Result& result = result_;
    const Options& options = options_;
    const Library& library = *library_;
    const auto& p = library.precursors();
    const auto& t = library.transitions();

    const std::size_t i = chromatogram.precursor;
    const std::uint32_t tb = chromatogram.transition_begin;
    const std::uint32_t tc = chromatogram.transition_count;
    if (tc == 0) { return; }

    const std::size_t points = chromatogram.pointCount(0);
    if (points < 3) { ++result.precursors_without_candidate; return; }

    // D8: each transition standardised against its own local noise before
    // summing, so no transition dominates by being loud and none is boosted
    // by what the library expects. See the option's comment for why library
    // weighting was rejected.
    const auto total = options.noise_normalised_picking
      ? noiseNormalisedTrace(chromatogram, points)
      : summedTrace(chromatogram, points);
    const double window_total = std::accumulate(total.begin(), total.end(), 0.0);
    if (window_total <= 0.0) { ++result.precursors_without_candidate; return; }

    const auto candidates = options.coelution_picking
      ? findCandidatesByCorrelation(chromatogram, options.corr_half_window,
                                    options.min_corr_score, options.max_corr_diff,
                                    options.apex_evidence, options.smooth_half_width,
                                    options.boundary_fraction, options.max_candidates)
      : findCandidates(smooth(total, options.smooth_half_width),
                       options.max_candidates, options.boundary_fraction);
    if (candidates.empty()) { ++result.precursors_without_candidate; return; }

    // Library intensities, in the transition order the chromatograms use.
    std::vector<double> library_intensity(tc, 0.0);
    for (std::uint32_t k = 0; k < tc; ++k)
    {
      library_intensity[k] = t.library_intensity[tb + k];
    }

    for (const auto& cand : candidates)
    {
      const std::size_t lo = cand.left, hi = cand.right;
      const std::size_t width = hi - lo + 1;

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
      std::size_t at_apex = 0;
      for (std::uint32_t k = 0; k < tc; ++k)
      {
        const std::uint32_t n = chromatogram.pointCount(k);
        const float* points_k = n ? chromatogram.trace(k) : nullptr;
        std::vector<double> whole(n, 0.0);
        for (std::uint32_t j = 0; j < n; ++j) { whole[j] = points_k[j]; }
        const double bg = localBackground(whole, lo, hi);
        corrected[k] = std::max(0.0, observed[k] - bg * static_cast<double>(width));
        if (cand.apex < n && points_k[cand.apex] > bg) { ++at_apex; }
      }

      // D6/D8 gate: a peak group is a co-elution. One transition above its
      // own background is a spike, and emitting it as a candidate is what
      // let single-fragment interference into the score matrix.
      if (at_apex < options.min_fragments_at_apex) { continue; }

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
      g.left_rt = chromatogram.retentionTime(static_cast<std::uint32_t>(lo));
      g.right_rt = chromatogram.retentionTime(static_cast<std::uint32_t>(hi));
      g.apex_intensity = static_cast<float>(cand.apex_value);

      g.sub_scores.assign(N_SUB_SCORES, 0.0);
      g.sub_scores[XCORR_SHAPE] = shape;
      // Negated so that, like every other column, larger is better. A
      // classifier would learn the sign either way, but a human reading a
      // weight vector should not have to remember which column is inverted.
      g.sub_scores[XCORR_COELUTION] = -coelution;
      g.sub_scores[LIBRARY_CORR] = pearson(corrected, library_intensity);
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
      const double floor_bg = std::max(background, 0.01 * cand.apex_value);
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
          if (t.type[tb + k] == FragmentType::Y) { y_area += corrected[k]; }
        }
        g.sub_scores[YSERIES_SCORE] = all_area > 0.0 ? y_area / all_area : 0.0;
        g.sub_scores[FRAGMENT_COVERAGE] =
          tc > 0 ? static_cast<double>(at_apex) / static_cast<double>(tc) : 0.0;
      }
      // A candidate whose spectrum does not resemble the library is not this
      // peptide, wherever it eluted.
      //
      // Measured on S08 against DIA-NN's confident set, splitting our own
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
      result.groups.push_back(std::move(g));
    }
  }

  PeakGroupScorer::Result PeakGroupScorer::Session::finish()
  {
    Result& result = result_;
    const Options& options = options_;
    const std::size_t n_precursors = library_->precursorCount();

    // Peak groups arrive in whatever order the extractor finished their
    // precursors, which is retention-time order rather than library order. Put
    // them back in library order before anything downstream sees them: the rows
    // of the feature matrix are these groups, and a classifier fitted on rows
    // ordered by elution time is a classifier that can depend on it. Stable, so
    // a precursor's candidates keep the order the search produced.
    std::stable_sort(result.groups.begin(), result.groups.end(),
                     [](const PeakGroup& a, const PeakGroup& b) {
                       return a.precursor < b.precursor; });

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
    // Threading is per-classifier, not on LDAParams: the LDA solve is a small
    // dense Cholesky and does not want threads, while the tree and network
    // fits do.
    if (options.classifier == "gbt")
    {
      params.classifier = Scoring::LDAParams::Classifier::GBT;
      params.gbt.n_threads = static_cast<int>(options.threads);
    }
    else if (options.classifier == "nn")
    {
      params.classifier = Scoring::LDAParams::Classifier::NN;
      params.nn.n_threads = static_cast<int>(options.threads);
    }

    const auto scored = Scoring::scoreSemiSupervisedLDA(features, labels, group, params);
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
    }
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
      if (tc == 0 || tb >= chromatograms.begin.size()) { continue; }

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

} // namespace ODIA
