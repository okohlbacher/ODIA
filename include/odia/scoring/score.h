// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Adopted from okohlbacher/OpenDIAlyzer, tag odia-v0.3.0 (0ed5fb2), file
// src/odia_score.h. Same author, same BSD-3 licence. Changed here: namespace
// `odia` -> `ODIA::Scoring`. The algorithms are untouched.
//
// Worth knowing why this is adopted rather than rewritten: every function is
// checked upstream against the exact reference values in OpenSWATH's
// Scoring_test.cpp, so these are OpenSWATH's score definitions rather than a
// plausible re-derivation of them. Rewriting would have thrown that away and
// left the cross-check against OpenSWATH still to do.
// odia-score -- the DIA scoring kernel, CPU reference implementation.
//
// This is the hot path a future CUDA port targets: for a candidate's F fragment
// traces over T retention-time points, reduce the F x T matrix to a score
// vector. The expensive term is the all-pairs cross-correlation (O(F^2 * D * T),
// D = 2*maxdelay+1), which is what OpenSWATH's MRMScoring computes and what
// carries enough arithmetic intensity (36-110 flop/byte FP32) to suit an H100.
//
// Correctness is defined by OpenSWATH's own algorithms: every function here is
// checked against the exact reference values in
// ext/OpenMS/src/tests/class_tests/openswathalgo/Scoring_test.cpp, so the CUDA
// port has a bit-comparable oracle. Idea and reference values from OpenMS
// (BSD-3); implemented independently in a layout suited to GPU porting (SoA,
// dense, no per-pair heap allocation).
//
// Header-only so the same source compiles for host and, later, device.

#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

namespace ODIA::Scoring
{

// Standardize in place: subtract mean, divide by *population* standard deviation
// (÷n, matching OpenSWATH's standardize_data, not the ÷(n-1) sample form).
inline void standardize(std::vector<double>& x)
{
  const size_t n = x.size();
  if (n == 0) return;
  double mean = 0.0;
  for (double v : x) mean += v;
  mean /= n;
  double var = 0.0;
  for (double v : x) var += (v - mean) * (v - mean);
  var /= n;                                  // population variance
  const double sd = std::sqrt(var);
  if (sd == 0.0) { for (double& v : x) v -= mean; return; }
  for (double& v : x) v = (v - mean) / sd;
}

// Cross-correlation at each integer delay in [-maxdelay, +maxdelay], normalised
// by 1/n. Inputs are assumed already standardized (the "...Post" variant in
// OpenSWATH). Returns one value per delay; the delay for slot k is k-maxdelay.
//
// GPU-shaped: output is a dense array indexed by delay, not a vector of
// (delay,value) pairs. One thread per delay is the natural device mapping.
inline std::vector<double> xcorr_post(const std::vector<double>& a,
                                      const std::vector<double>& b, int maxdelay)
{
  const int n = static_cast<int>(a.size());
  std::vector<double> out(2 * maxdelay + 1, 0.0);
  for (int delay = -maxdelay; delay <= maxdelay; ++delay)
  {
    // overlap of a[i] with b[i+delay]
    const int lo = delay < 0 ? -delay : 0;
    const int hi = delay < 0 ? n : n - delay;
    double s = 0.0;
    for (int i = lo; i < hi; ++i) s += a[i] * b[i + delay];
    out[delay + maxdelay] = s / n;
  }
  return out;
}

// Convenience: standardize copies, then cross-correlate. Mirrors
// normalizedCrossCorrelation.
inline std::vector<double> xcorr(std::vector<double> a, std::vector<double> b, int maxdelay)
{
  standardize(a);
  standardize(b);
  return xcorr_post(a, b, maxdelay);
}

// Largest |value| in a delay array, and its delay. This is the co-elution
// score OpenSWATH extracts from the xcorr array (xcorrArrayGetMaxPeak).
struct XCorrPeak { int delay; double value; };
inline XCorrPeak xcorr_max(const std::vector<double>& xc, int maxdelay)
{
  XCorrPeak best{0, 0.0};
  double bestabs = -1.0;
  for (size_t k = 0; k < xc.size(); ++k)
    if (std::abs(xc[k]) > bestabs) { bestabs = std::abs(xc[k]); best = {int(k) - maxdelay, xc[k]}; }
  return best;
}

// Normalised Manhattan distance (mQuest delta_ratio_sum). OpenSWATH's
// normalize_sum divides each vector by its SUM (not its mean), then takes the
// mean absolute difference. Dividing by the mean instead is an n-fold error.
inline double normalized_manhattan(std::vector<double> x, std::vector<double> y)
{
  const size_t n = x.size();
  if (n == 0) return 0.0;
  double sx = 0.0, sy = 0.0;
  for (size_t i = 0; i < n; ++i) { sx += x[i]; sy += y[i]; }
  if (sx != 0.0) for (double& v : x) v /= sx;
  if (sy != 0.0) for (double& v : y) v /= sy;
  double s = 0.0;
  for (size_t i = 0; i < n; ++i) s += std::abs(x[i] - y[i]);
  return s / n;
}

// Root-mean-square deviation. Matches RootMeanSquareDeviation.
inline double rmsd(const std::vector<double>& x, const std::vector<double>& y)
{
  const size_t n = x.size();
  if (n == 0) return 0.0;
  double s = 0.0;
  for (size_t i = 0; i < n; ++i) s += (x[i] - y[i]) * (x[i] - y[i]);
  return std::sqrt(s / n);
}

// Spectral angle = acos(normalised dot product). Matches SpectralAngle.
inline double spectral_angle(const std::vector<double>& x, const std::vector<double>& y)
{
  const size_t n = x.size();
  double dot = 0.0, nx = 0.0, ny = 0.0;
  for (size_t i = 0; i < n; ++i) { dot += x[i] * y[i]; nx += x[i] * x[i]; ny += y[i] * y[i]; }
  const double denom = std::sqrt(nx) * std::sqrt(ny);
  if (denom == 0.0) return 0.0;
  double c = dot / denom;
  if (c > 1.0) c = 1.0; else if (c < -1.0) c = -1.0;
  return std::acos(c);
}

// The all-pairs cross-correlation matrix over F fragment traces -- the term
// that dominates cost and that a CUDA port parallelises. Upper triangle only
// (the matrix is symmetric in |delay|), matching MRMScoring::initializeXCorrMatrix.
// traces: F vectors, each length T, standardized in place. Returns, for each
// unordered pair (i<=j), the max-peak co-elution value.
struct PairScore { int i, j, delay; double value; };
inline std::vector<PairScore> allpairs_xcorr(std::vector<std::vector<double>> traces, int maxdelay)
{
  const int F = static_cast<int>(traces.size());
  for (auto& t : traces) standardize(t);
  std::vector<PairScore> out;
  out.reserve(static_cast<size_t>(F) * (F + 1) / 2);
  for (int i = 0; i < F; ++i)
    for (int j = i; j < F; ++j)
    {
      const auto xc = xcorr_post(traces[i], traces[j], maxdelay);
      const auto pk = xcorr_max(xc, maxdelay);
      out.push_back({i, j, pk.delay, pk.value});
    }
  return out;
}


// ---------------------------------------------------------------------------
// ODIA additions. The primitives above are checked against OpenSWATH's own
// reference values and are deliberately NOT modified -- that verification is
// the reason this file was adopted rather than rewritten. What follows are
// separate entry points for the places where OpenSWATH's exact definition is
// wrong for our data, so both remain available and the divergence is explicit.

/// Like xcorr_max, but selects the largest SIGNED value.
///
/// xcorr_max selects on |value| and returns the signed one, so a strong
/// anti-correlation at some lag beats a moderate genuine co-elution at lag 0
/// -- and the negative number then lands in the shape score, which is supposed
/// to reward co-elution. That is right for "find the dominant feature of the
/// correlogram" and wrong for "how well do these two fragments co-elute".
inline XCorrPeak xcorr_max_signed(const std::vector<double>& xc, int maxdelay)
{
  XCorrPeak best{0, -2.0};
  bool any = false;
  for (size_t k = 0; k < xc.size(); ++k)
  {
    if (!any || xc[k] > best.value)
    {
      best = {static_cast<int>(k) - maxdelay, xc[k]};
      any = true;
    }
  }
  return any ? best : XCorrPeak{0, 0.0};
}

/// True when a trace carries no information: constant, or all zero.
///
/// standardize() maps such a trace to all zeros, after which every pair it
/// joins scores 0.0 at every lag -- and xcorr_max's `bestabs = -1.0` seed then
/// reports lag -maxdelay for it. Counted as evidence, these measure how many
/// transitions are dead rather than how well the live ones agree.
inline bool degenerate(const std::vector<double>& t)
{
  if (t.size() < 2) { return true; }
  const double first = t.front();
  for (const double v : t) { if (v != first) { return false; } }
  return true;
}

struct PairOptions
{
  /// Self-pairs are autocorrelations: exactly 1.0 at lag 0, always. With 12
  /// transitions they are 12 of 78 pairs, adding a constant 0.1538 to the mean
  /// shape of every group, good or bad, and compressing the measured
  /// target/decoy separation by 18%.
  bool exclude_self = true;

  /// A candidate window is 5-9 points; a lag of 10 evaluates delays whose
  /// overlap is empty or one point. Capped at (n-1)/2 so every reported lag
  /// rests on at least half the trace.
  bool cap_delay_to_trace = true;

  /// Drop degenerate traces rather than scoring them as agreement.
  bool skip_degenerate = true;

  /// Shape wants the best co-elution, not the biggest excursion.
  bool signed_selection = true;
};

/// Pairwise cross-correlation with the corrections above.
///
/// `usable` receives the number of traces that carried information, which is
/// worth having as a feature in its own right: a precursor scored from three
/// live fragments is not the same evidence as one scored from twelve, and
/// without it the two are indistinguishable in the score.
inline std::vector<PairScore> allpairs_xcorr_ex(std::vector<std::vector<double>> traces,
                                                int maxdelay, const PairOptions& opt,
                                                std::size_t* usable = nullptr)
{
  std::vector<std::vector<double>> live;
  live.reserve(traces.size());
  for (auto& t : traces)
  {
    if (opt.skip_degenerate && degenerate(t)) { continue; }
    live.push_back(std::move(t));
  }
  if (usable) { *usable = live.size(); }

  const int F = static_cast<int>(live.size());
  std::vector<PairScore> out;
  if (F < 2) { return out; }

  int md = maxdelay;
  if (opt.cap_delay_to_trace)
  {
    std::size_t shortest = live.front().size();
    for (const auto& t : live) { shortest = std::min(shortest, t.size()); }
    md = std::max(1, std::min(md, static_cast<int>((shortest - 1) / 2)));
  }

  for (auto& t : live) { standardize(t); }
  for (int i = 0; i < F; ++i)
  {
    for (int j = opt.exclude_self ? i + 1 : i; j < F; ++j)
    {
      const auto xc = xcorr_post(live[i], live[j], md);
      const auto pk = opt.signed_selection ? xcorr_max_signed(xc, md) : xcorr_max(xc, md);
      out.push_back({i, j, pk.delay, pk.value});
    }
  }
  return out;
}

} // namespace ODIA::Scoring
