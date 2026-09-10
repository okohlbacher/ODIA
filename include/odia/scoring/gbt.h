// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Adopted from okohlbacher/OpenDIAlyzer, branch odia-engine-and-openms-boundary,
// tag odia-v0.3.0 (0ed5fb2), file src/odia_gbt.h. Same author, same BSD-3 licence.
//
// Changed here: namespace `odia` -> `ODIA::Scoring` to match this project, and
// the sibling includes rewritten to this tree's layout. The algorithms are
// untouched.
//
// Provenance that travels with the code, per that repository's CLEAN-ROOM.md:
// this is an independent reimplementation of a published confidence model
// (Reiter et al. mProphet, Nat Methods 2011; Teleman et al. pyprophet,
// Bioinformatics 2015). No DIA-NN code is present in it.
// odia_gbt.h — histogram gradient-boosted trees, as a drop-in replacement for the LDA's linear
// discriminant inside the same semi-supervised loop.
//
// WHY NOT XGBoost/LightGBM. Neither is available to link here (checked: nothing in the OpenMS 3.6
// ML/ tree, nothing in the environment), and vendoring one would put a large build dependency into a
// tool whose scoring layer is otherwise three standalone headers with self-checks. A histogram GBT
// is ~400 lines, has no dependencies, and — unlike a vendored library with its own threading and
// RNG — is trivially deterministic, which matters because these scores feed an FDR estimate.
//
// WHY BOTHER. The LDA fits ONE hyperplane through the sub-score space. OpenSWATH's sub-scores
// interact non-linearly in ways a plane cannot express: a high xcorr_shape means something very
// different at high vs low signal-to-noise, and library correlation matters more when many
// transitions are present. mokapot (Fondrie & Noble 2021) is exactly this substitution — Percolator
// with boosted trees instead of the linear model — and reports +15% modified PSMs, +19% peptides,
// +11% proteins at 1% FDR. pyprophet ships XGBoost as a non-default classifier for the same reason.
// Expect less on vanilla tryptic bulk DIA than those headline numbers, which come from
// immunopeptidomics and PTM work; see docs/OpenDIAlyzer-scoring-fdr-backlog.md, which flags exactly
// this and says to read the gains skeptically.
//
// WHAT IT IS. Binary logistic gradient boosting, second-order (Newton) leaf values, feature values
// pre-binned so split finding is O(bins) per feature per node rather than O(rows log rows):
//
//   - bin edges from quantiles of the TRAINING rows only, computed once per fit
//   - level-wise growth to a fixed depth (not leaf-wise): bounded tree shape, no surprise depth
//   - gain = 1/2 [ GL^2/(HL+lambda) + GR^2/(HR+lambda) - (GL+GR)^2/(HL+HR+lambda) ] - gamma
//   - leaf value = -G/(H+lambda), shrunk by the learning rate
//   - no row or column subsampling, no early stopping on a random split -> bit-identical every run
//
// DETERMINISM IS A REQUIREMENT, NOT A PREFERENCE. These scores produce q-values. A classifier whose
// output moves between runs makes the FDR irreproducible, so every source of randomness that
// ordinary GBT implementations use for regularisation is deliberately absent here. Ties in split
// selection break on the lowest (feature, bin) index.

#ifndef ODIA_GBT_H
#define ODIA_GBT_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <ios>
#include <istream>
#include <limits>
#include <locale>
#include <ostream>
#include <string>
#include <utility>
#include <numeric>
#include <vector>

namespace ODIA::Scoring
{

struct GBTParams
{
  int n_trees = 120;              ///< boosting rounds
  int max_depth = 4;              ///< level-wise depth; 4 gives <=16 leaves, enough for pairwise
                                  ///< interactions between sub-scores without memorising rows
  int n_bins = 64;                ///< histogram resolution per feature
  double learning_rate = 0.1;     ///< shrinkage
  double lambda = 1.0;            ///< L2 on leaf values
  double gamma = 0.0;             ///< minimum gain to split
  double min_child_weight = 1.0;  ///< minimum summed hessian in a child
  int min_child_rows = 20;        ///< minimum rows in a child; guards tiny leaves on small folds
  double max_delta_step = 0.0;    ///< XGBoost's max_delta_step: a leaf's Newton step G/(H+lambda) is
                                  ///< clipped to +/-this BEFORE the learning rate. 0 = off (the native
                                  ///< unbounded step). Bounds the first rounds of a fit whose prior is
                                  ///< extreme, where a pure-positive leaf's hessian is tiny against its
                                  ///< gradient and the unbounded step is |G|/lambda. EXPERIMENT knob.
  bool intercept_zero = false;    ///< Intercept 0 (XGBoost base_score 0.5) instead of logit(training prior).
                                  ///< At a 1:155 prior a pure-positive leaf's Newton step is ~1/p (~160);
                                  ///< at intercept 0 it is ~2. EXPERIMENT knob (v1.11).
  int warmup_rounds = 0;          ///< Learning-rate ramp: round r (1-based) uses lr * min(1, r/warmup_rounds);
                                  ///< 0 = off (constant lr, identical arithmetic). EXPERIMENT knob (v1.12).
  bool clip_gain = false;         ///< With max_delta_step > 0, evaluate split gains with the CLIPPED leaf
                                  ///< weight (XGBoost CalcGain semantics) instead of the unclipped G^2/(H+lambda),
                                  ///< so splits are ranked by the improvement they can actually realise.
                                  ///< Identical arithmetic when off. EXPERIMENT knob (v1.11).
  /// Bin edges from ALL rows rather than from the current training set.
  ///
  /// The default recomputes quantile edges from `pos + neg`, which in a
  /// semi-supervised loop CHANGES every iteration. That is an amplification
  /// path: adding one column changes the model, which changes the selected
  /// positives, which re-discretises EVERY OTHER column -- so a perturbation
  /// that should be local becomes global. Measured consequence: a column of
  /// pure noise moves identifications by up to 11.9%, with Spearman rho 0.79
  /// between a run and the same run plus that column.
  ///
  /// Edges over all rows are label-independent and therefore identical across
  /// iterations and across feature sets, which closes the loop. Costs one extra
  /// pass over the unlabelled rows per fit.
  bool fixed_bins = false;
  /// Threads for the histogram pass (0 = whatever OpenMP gives). fit() is normally called from
  /// inside the fold loop's parallel region, where a nested team defaults to ONE thread -- so the
  /// caller must both set this and enable a second active level, or the parallelism does nothing.
  /// Results do NOT depend on this value; see the chunking note in growTree_.
  int n_threads = 0;

  /// pyProphet 3.0.15's XGBoost 3.2.0 configuration, verbatim.
  ///
  /// Six of the nine parameters it sets already match this struct's defaults --
  /// lambda 1, gamma 0, min_child_weight 1, and subsample / colsample_bytree /
  /// alpha at 1/1/0, which are no-ops and therefore equivalent to this
  /// implementation not having them. The algorithm is already XGBoost's: same
  /// second-order Newton leaf values, same split gain.
  ///
  /// The two that differ are the two that matter, and both make the model
  /// stronger and more prone to overfitting: depth 4 -> 6 takes the leaf count
  /// from <=16 to <=64, and eta 0.1 -> 0.3 triples the step. This project's
  /// defaults were chosen for "pairwise interactions between sub-scores without
  /// memorising rows" with 10 features; we now have 15.
  static GBTParams pyprophet()
  {
    GBTParams p;
    p.max_depth = 6;          // XGBoost "max_depth": 6
    p.learning_rate = 0.3;    // XGBoost "eta": 0.3
    p.lambda = 1.0;           // "lambda": 1   -- already the default here
    p.gamma = 0.0;            // "gamma": 0
    p.min_child_weight = 1.0; // "min_child_weight": 1
    return p;
  }
};

namespace gbt_detail
{

/// Gradient and hessian for one histogram bin, INTERLEAVED.
///
/// They were two separate arrays. Every accumulation touches the same index k in both, and k is the
/// data-dependent bin, so the access is effectively random: two arrays meant two cache lines pulled
/// in per (row, feature) in the loop that is 92.6% of a fit -- 13.9e9 potentially-missing accesses
/// over a fit. Interleaved, g and h are 16 adjacent bytes and land in ONE line, halving the misses
/// for identical arithmetic in identical order (so the result is bit-unchanged).
struct GH
{
  double g = 0.0;
  double h = 0.0;
  /// Rows in this bin. Needed to enforce min_child_rows PER CHILD: without it
  /// `nL` below has nothing to accumulate, which is exactly how the predecessor
  /// left the same guard dead (its handoff, open defect 2: "nL is declared,
  /// never incremented, discarded with (void)nL ... leaves of ~4 rows are
  /// possible -- an overfitting surface in exactly the score tail that sets the
  /// FDR threshold").
  std::size_t n = 0;
};

/// One level-wise regression tree over pre-binned features.
struct Tree
{
  /// Nodes are stored as a complete binary tree indexed 1..2^(d+1)-1 so children of i are 2i, 2i+1;
  /// this costs a little memory at depth 4 (31 slots) and removes all pointer chasing.
  std::vector<int> feature;      ///< split feature, or -1 for a leaf
  std::vector<uint8_t> bin;      ///< split threshold: go left when bin <= this
  std::vector<double> value;     ///< leaf value (valid where feature == -1)
  int max_depth = 0;

  double predict(const uint8_t* row) const
  {
    int node = 1;
    for (int d = 0; d < max_depth; ++d)
    {
      if (feature[node] < 0) { break; }
      node = (row[feature[node]] <= bin[node]) ? 2 * node : 2 * node + 1;
    }
    return value[node];
  }
};

/// Quantile bin edges per feature, from the training rows. Returns edges[f] with <= n_bins-1
/// thresholds; a value goes to bin b when it is <= edges[f][b] (last bin catches the rest).
inline std::vector<std::vector<double>> computeBinEdges(
  const std::vector<std::vector<double>>& X, const std::vector<std::size_t>& rows, int n_bins)
{
  const std::size_t m = X.empty() ? 0 : X[0].size();
  std::vector<std::vector<double>> edges(m);
  std::vector<double> col;
  col.reserve(rows.size());
  for (std::size_t f = 0; f < m; ++f)
  {
    col.clear();
    for (const std::size_t r : rows)
    {
      const double v = X[r][f];
      if (std::isfinite(v)) { col.push_back(v); }
    }
    if (col.empty()) { continue; }
    std::sort(col.begin(), col.end());
    // Distinct quantile cut points. Duplicates are dropped, so a feature with few distinct values
    // simply gets few bins instead of many identical ones (which would waste histogram slots and
    // create splits that separate nothing).
    for (int b = 1; b < n_bins; ++b)
    {
      const std::size_t idx =
        static_cast<std::size_t>(static_cast<double>(col.size()) * b / n_bins);
      const double e = col[std::min(idx, col.size() - 1)];
      if (edges[f].empty() || e > edges[f].back()) { edges[f].push_back(e); }
    }
  }
  return edges;
}

/// Map every row to its bin index per feature. Non-finite values go to the LAST bin, so "missing"
/// is a value the tree can split on rather than something silently imputed to a mean.
inline std::vector<std::vector<uint8_t>> binFeatures(
  const std::vector<std::vector<double>>& X, const std::vector<std::vector<double>>& edges)
{
  const std::size_t n = X.size(), m = edges.size();
  std::vector<std::vector<uint8_t>> B(n, std::vector<uint8_t>(m, 0));
  for (std::size_t i = 0; i < n; ++i)
  {
    for (std::size_t f = 0; f < m; ++f)
    {
      const double v = X[i][f];
      if (!std::isfinite(v))
      {
        B[i][f] = static_cast<uint8_t>(edges[f].size());
        continue;
      }
      const auto it = std::lower_bound(edges[f].begin(), edges[f].end(), v);
      B[i][f] = static_cast<uint8_t>(it - edges[f].begin());
    }
  }
  return B;
}

} // namespace gbt_detail

/// Gradient-boosted trees for binary classification, trained on explicit positive/negative row sets
/// so it plugs into the same semi-supervised loop as the LDA.
class GBT
{
public:
  /// Fit on rows `pos` (label 1) and `neg` (label 0) of X. Returns false if the data cannot support
  /// a fit, in which case the caller should keep its previous model rather than use a degenerate one.
  bool fit(const std::vector<std::vector<double>>& X, const std::vector<std::size_t>& pos,
           const std::vector<std::size_t>& neg, const GBTParams& p)
  {
    if (pos.size() < 2 || neg.size() < 2 || X.empty()) { return false; }
    params_ = p;
    n_features_ = X[0].size();
    if (n_features_ == 0) { return false; }

    std::vector<std::size_t> rows;
    rows.reserve(pos.size() + neg.size());
    rows.insert(rows.end(), pos.begin(), pos.end());
    rows.insert(rows.end(), neg.begin(), neg.end());
    std::vector<double> y(rows.size(), 0.0);
    for (std::size_t i = 0; i < pos.size(); ++i) { y[i] = 1.0; }

    if (p.fixed_bins)
    {
      std::vector<std::size_t> all_rows(X.size());
      std::iota(all_rows.begin(), all_rows.end(), std::size_t{0});
      edges_ = gbt_detail::computeBinEdges(X, all_rows, p.n_bins);
    }
    else
    {
      edges_ = gbt_detail::computeBinEdges(X, rows, p.n_bins);
    }
    // A feature with no edges cannot split; that is fine, it just never wins a gain comparison.
    //
    // FLAT, not vector<vector>: the histogram pass below is 93% of a fit and reads every feature of
    // every row. One heap block per row would put a pointer chase in front of each of those reads
    // and leave consecutive rows scattered; a flat buffer makes the whole pass a contiguous sweep.
    const std::size_t n_rows = rows.size();
    std::vector<uint8_t> B(n_rows * n_features_, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(p.n_threads > 0 ? p.n_threads : omp_get_max_threads())
#endif
    for (long long ii = 0; ii < static_cast<long long>(n_rows); ++ii)
    {
      const std::size_t i = static_cast<std::size_t>(ii);
      for (std::size_t f = 0; f < n_features_; ++f)
      {
        const double v = X[rows[i]][f];
        if (!std::isfinite(v)) { B[i * n_features_ + f] = static_cast<uint8_t>(edges_[f].size()); continue; }
        const auto it = std::lower_bound(edges_[f].begin(), edges_[f].end(), v);
        B[i * n_features_ + f] = static_cast<uint8_t>(it - edges_[f].begin());
      }
    }

    // Base score = log-odds of the training prior, so the first tree corrects a residual rather
    // than having to discover the class balance.
    const double frac = static_cast<double>(pos.size()) / static_cast<double>(rows.size());
    const double clamped = std::min(1.0 - 1e-6, std::max(1e-6, frac));
    base_ = p.intercept_zero ? 0.0 : std::log(clamped / (1.0 - clamped));

    std::vector<double> pred(n_rows, base_);
    std::vector<double> grad(n_rows), hess(n_rows);
    trees_.clear();
    trees_.reserve(static_cast<std::size_t>(std::max(0, p.n_trees)));

    std::size_t n_bins_max = 1;
    for (const auto& e : edges_) { n_bins_max = std::max(n_bins_max, e.size() + 1); }

    // One scratch buffer for the whole fit. growTree_ clears only the slice each level needs, so
    // the shape is set by the WIDEST level (2^(depth-1) nodes) and every level after the first
    // reuses the same allocation.
    std::vector<gbt_detail::GH> scratch;

    const int T = p.n_threads > 0 ? p.n_threads :
#ifdef _OPENMP
      omp_get_max_threads();
#else
      1;
#endif
    for (int t = 0; t < p.n_trees; ++t)
    {
      // Element-wise, so no reduction and no ordering question.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(T)
#endif
      for (long long ii = 0; ii < static_cast<long long>(n_rows); ++ii)
      {
        const std::size_t i = static_cast<std::size_t>(ii);
        const double pi = 1.0 / (1.0 + std::exp(-pred[i]));
        grad[i] = pi - y[i];
        hess[i] = std::max(1e-12, pi * (1.0 - pi));
      }
      round_lr_ = params_.learning_rate;
      if (p.warmup_rounds > 0)
      { round_lr_ = params_.learning_rate * std::min(1.0, static_cast<double>(t + 1) / static_cast<double>(p.warmup_rounds)); }
      gbt_detail::Tree tree;
      if (!growTree_(B, n_rows, grad, hess, n_bins_max, T, scratch, tree)) { break; }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(T)
#endif
      for (long long ii = 0; ii < static_cast<long long>(n_rows); ++ii)
      {
        const std::size_t i = static_cast<std::size_t>(ii);
        pred[i] += tree.predict(&B[i * n_features_]);
      }
      trees_.push_back(std::move(tree));
    }
    return !trees_.empty();
  }

  /// Raw margin (log-odds). Monotone in the probability, so it can be ranked directly.
  double score(const std::vector<double>& x) const
  {
    if (trees_.empty()) { return 0.0; }
    std::vector<uint8_t> b(n_features_, 0);
    for (std::size_t f = 0; f < n_features_; ++f)
    {
      const double v = x[f];
      if (!std::isfinite(v)) { b[f] = static_cast<uint8_t>(edges_[f].size()); continue; }
      const auto it = std::lower_bound(edges_[f].begin(), edges_[f].end(), v);
      b[f] = static_cast<uint8_t>(it - edges_[f].begin());
    }
    double s = base_;
    for (const auto& t : trees_) { s += t.predict(b.data()); }
    return s;
  }

  bool trained() const { return !trees_.empty(); }
  std::size_t featureCount() const { return n_features_; }

  /// Serialise everything `score()` consumes -- the bin edges, the base score and the trees --
  /// as plain text, one token per value, doubles at full precision. The training parameters are
  /// not needed to apply a model and are not written. Returns false on a write error.
  bool save(std::ostream& os) const
  {
    if (trees_.empty()) { return false; }
    os.imbue(std::locale::classic());   // "1.5", never "1,5", whatever the process locale
    os.unsetf(std::ios::floatfield);    // default format: a caller's std::fixed would zero 1e-20
    os.precision(17);
    os << "ODIA-GBT 1\n" << n_features_ << ' ' << base_ << ' ' << trees_.size() << '\n';
    for (const auto& e : edges_)
    {
      os << e.size();
      for (const double v : e) { os << ' ' << v; }
      os << '\n';
    }
    for (const auto& t : trees_)
    {
      os << t.max_depth << ' ' << t.feature.size() << '\n';
      for (std::size_t i = 0; i < t.feature.size(); ++i)
      {
        os << t.feature[i] << ' ' << static_cast<int>(t.bin[i]) << ' ' << t.value[i] << '\n';
      }
    }
    return static_cast<bool>(os);
  }

  /// Inverse of save(). A model loaded here scores exactly as the one that was saved, provided
  /// the caller feeds it features standardised the way the SAVING run standardised them -- that
  /// transform is the caller's to carry (see scoreSemiSupervisedLDA), not the tree's.
  bool load(std::istream& is)
  {
    // Every count and index is bounded before it sizes a vector or indexes one:
    // a file that parses but describes an impossible tree must fail here, not
    // in predict(). Bounds are generous against anything fit() can produce.
    is.imbue(std::locale::classic());
    std::string magic;
    int version = 0;
    if (!(is >> magic >> version) || magic != "ODIA-GBT" || version != 1) { return false; }
    std::size_t n_trees = 0;
    if (!(is >> n_features_ >> base_ >> n_trees) || n_features_ == 0 || n_features_ > 4096 ||
        n_trees == 0 || n_trees > 100000 || !std::isfinite(base_))
    { return false; }
    edges_.assign(n_features_, {});
    for (auto& e : edges_)
    {
      std::size_t k = 0;
      if (!(is >> k) || k > 255) { return false; }   // bins are uint8_t; <= 255 edges
      e.resize(k);
      for (double& v : e) { if (!(is >> v) || !std::isfinite(v)) { return false; } }
      // score() bins with lower_bound: the edges must be non-decreasing or a
      // value lands in an arbitrary bin.
      for (std::size_t i = 1; i < e.size(); ++i) { if (e[i] < e[i - 1]) { return false; } }
    }
    trees_.clear();
    trees_.reserve(n_trees);
    for (std::size_t t = 0; t < n_trees; ++t)
    {
      gbt_detail::Tree tree;
      std::size_t nodes = 0;
      if (!(is >> tree.max_depth >> nodes)) { return false; }
      // A complete binary tree indexed 1..2^(d+1)-1 is stored in 2^(d+1) slots.
      if (tree.max_depth < 0 || tree.max_depth > 16 ||
          nodes != (static_cast<std::size_t>(1) << (tree.max_depth + 1)))
      { return false; }
      tree.feature.resize(nodes);
      tree.bin.resize(nodes);
      tree.value.resize(nodes);
      for (std::size_t i = 0; i < nodes; ++i)
      {
        int b = 0;
        if (!(is >> tree.feature[i] >> b >> tree.value[i])) { return false; }
        if (tree.feature[i] < -1 || tree.feature[i] >= static_cast<int>(n_features_) ||
            b < 0 || b > 255 || !std::isfinite(tree.value[i]))
        { return false; }
        tree.bin[i] = static_cast<uint8_t>(b);
      }
      trees_.push_back(std::move(tree));
    }
    return !trees_.empty();
  }
  std::size_t nTrees() const { return trees_.size(); }

private:
  /// Level-wise growth. `node_of_row` tracks which node each row is in; a level is built by
  /// histogramming (grad, hess) per (node, feature, bin) and picking the best split per node.
  bool growTree_(const std::vector<uint8_t>& B, std::size_t n_rows, const std::vector<double>& grad,
                 const std::vector<double>& hess, std::size_t n_bins_max, int n_threads,
                 std::vector<gbt_detail::GH>& scratch, gbt_detail::Tree& tree)
  {
    const int D = std::max(1, params_.max_depth);
    const std::size_t n_nodes = (std::size_t(1) << (D + 1));
    tree.max_depth = D;
    tree.feature.assign(n_nodes, -1);
    tree.bin.assign(n_nodes, 0);
    tree.value.assign(n_nodes, 0.0);

    std::vector<int> node_of_row(n_rows, 1);
    bool any_split = false;

    for (int depth = 0; depth < D; ++depth)
    {
      const int first = 1 << depth, last = (1 << (depth + 1)) - 1;
      // histogram[node][feature][bin] -> (G, H)
      const std::size_t n_at_level = static_cast<std::size_t>(last - first + 1);
      const std::size_t hist_size = n_at_level * n_features_ * n_bins_max;
      std::vector<gbt_detail::GH> H(hist_size);
      std::vector<std::size_t> node_rows(n_at_level, 0);

      // This loop is ~93% of a fit. It is a REDUCTION over rows, so threads need private buffers.
      //
      // The partition is by a fixed number of CHUNKS, not by thread: a per-thread partition would
      // make the number of partial sums -- and therefore the floating-point reduction order --
      // depend on the thread count, and these scores become q-values. Chunk count depends only on
      // the row count, threads take chunks in any order, and the combination below runs in chunk
      // index order, so the result is bit-identical at 1, 8 or 224 threads.
      const int nchunk = static_cast<int>(std::min<std::size_t>(32, std::max<std::size_t>(1, n_rows / 1000)));
      // Scratch is owned by fit() and reused across every level of every tree: allocating it here
      // cost 960 vector constructions per fit (2 per level x 4 levels x 120 trees) for buffers of
      // identical shape. Only the portion this level uses is cleared.
      const std::size_t need = static_cast<std::size_t>(nchunk) * hist_size;
      if (scratch.size() < need) { scratch.resize(need); }
      std::fill(scratch.begin(), scratch.begin() + static_cast<std::ptrdiff_t>(need), gbt_detail::GH{});
      std::vector<std::size_t> node_rows_c(static_cast<std::size_t>(nchunk) * n_at_level, 0);

#ifdef _OPENMP
// ponytail: the loop has exactly nchunk iterations, so threads beyond nchunk get no work and only
// add barrier and scratch-page cost. Measured: 224 threads ran this fit 60% SLOWER than 32.
// Raise the cap by raising nchunk -- which changes the reduction order, so it is a deliberate act.
#pragma omp parallel for schedule(dynamic, 1) num_threads(std::min(n_threads > 0 ? n_threads : 1, nchunk))
#endif
      for (int c = 0; c < nchunk; ++c)
      {
        const std::size_t lo = n_rows * static_cast<std::size_t>(c) / static_cast<std::size_t>(nchunk);
        const std::size_t hi = n_rows * static_cast<std::size_t>(c + 1) / static_cast<std::size_t>(nchunk);
        gbt_detail::GH* __restrict h = scratch.data() + static_cast<std::size_t>(c) * hist_size;
        std::size_t* __restrict nr = node_rows_c.data() + static_cast<std::size_t>(c) * n_at_level;
        for (std::size_t i = lo; i < hi; ++i)
        {
          const int nd = node_of_row[i];
          if (nd < first || nd > last) { continue; }
          const std::size_t slot = static_cast<std::size_t>(nd - first);
          ++nr[slot];
          const std::size_t base = (slot * n_features_) * n_bins_max;
          const uint8_t* __restrict brow = B.data() + i * n_features_;
          const double gi = grad[i], hi_ = hess[i];
          for (std::size_t f = 0; f < n_features_; ++f)
          {
            const std::size_t k = base + f * n_bins_max + brow[f];
            h[k].g += gi;
            h[k].h += hi_;
            h[k].n += 1;
          }
        }
      }
      // Fixed-order combination -- this is what keeps the fit thread-count independent.
      for (int c = 0; c < nchunk; ++c)
      {
        const gbt_detail::GH* hc = scratch.data() + static_cast<std::size_t>(c) * hist_size;
        const std::size_t* nr = node_rows_c.data() + static_cast<std::size_t>(c) * n_at_level;
        for (std::size_t k = 0; k < hist_size; ++k)
        { H[k].g += hc[k].g; H[k].h += hc[k].h; H[k].n += hc[k].n; }
        for (std::size_t k = 0; k < n_at_level; ++k) { node_rows[k] += nr[k]; }
      }

      for (int nd = first; nd <= last; ++nd)
      {
        const std::size_t slot = static_cast<std::size_t>(nd - first);
        if (node_rows[slot] == 0) { continue; }
        double G = 0.0, Hs = 0.0;
        const std::size_t base = (slot * n_features_) * n_bins_max;
        for (std::size_t b = 0; b < n_bins_max; ++b)
        {
          G += H[base + b].g;
          Hs += H[base + b].h;
        }
        const double parent_obj = nodeObj_(G, Hs);

        int best_f = -1;
        std::size_t best_b = 0;
        double best_gain = params_.gamma;
        for (std::size_t f = 0; f < n_features_; ++f)
        {
          double GL = 0.0, HL = 0.0;
          std::size_t nL = 0;
          const std::size_t fb = base + f * n_bins_max;
          // Scan cut points left to right; the last bin is never a cut (it would put everything left).
          for (std::size_t b = 0; b + 1 < n_bins_max; ++b)
          {
            GL += H[fb + b].g;
            HL += H[fb + b].h;
            nL += H[fb + b].n;
            const double GR = G - GL, HR = Hs - HL;
            if (HL < params_.min_child_weight || HR < params_.min_child_weight) { continue; }
            // PER-CHILD row floor, now actually enforced. The parent check below
            // (node_rows >= 2*min_child_rows) permits a 39/1 split; this forbids
            // it. A leaf of a handful of rows is fitted noise, and it lands in
            // the score tail that sets the FDR threshold.
            const std::size_t nR = node_rows[slot] - nL;
            if (nL < static_cast<std::size_t>(params_.min_child_rows) ||
                nR < static_cast<std::size_t>(params_.min_child_rows)) { continue; }
            const double gain = 0.5 * (nodeObj_(GL, HL) + nodeObj_(GR, HR) - parent_obj);
            // Strictly greater: ties keep the lowest (feature, bin), which is what makes the fit
            // reproducible independent of iteration order.
            if (gain > best_gain)
            {
              best_gain = gain;
              best_f = static_cast<int>(f);
              best_b = b;
            }
          }
        }

        const bool leaf_level = (depth == D - 1);
        if (best_f < 0 || leaf_level || node_rows[slot] < static_cast<std::size_t>(2 * params_.min_child_rows))
        {
          tree.feature[nd] = -1;
          tree.value[nd] = -round_lr_ * leafStep_(G, Hs);
          continue;
        }
        tree.feature[nd] = best_f;
        tree.bin[nd] = static_cast<uint8_t>(best_b);
        any_split = true;
        // Children inherit a provisional leaf value in case the next level does not split them.
        for (int child : {2 * nd, 2 * nd + 1})
        {
          tree.feature[child] = -1;
          tree.value[child] = 0.0;
        }
      }

      // Route rows into the next level and give the leaves their values.
      for (std::size_t i = 0; i < n_rows; ++i)
      {
        const int nd = node_of_row[i];
        if (nd < first || nd > last || tree.feature[nd] < 0) { continue; }
        node_of_row[i] = (B[i * n_features_ + tree.feature[nd]] <= tree.bin[nd]) ? 2 * nd : 2 * nd + 1;
      }
    }

    // Final level: every row now sits in a node whose value must be set from its own G/H.
    const int first = 1 << D, last = (1 << (D + 1)) - 1;
    std::vector<double> G(static_cast<std::size_t>(last - first + 1), 0.0), H(G.size(), 0.0);
    for (std::size_t i = 0; i < n_rows; ++i)
    {
      const int nd = node_of_row[i];
      if (nd < first || nd > last) { continue; }
      G[static_cast<std::size_t>(nd - first)] += grad[i];
      H[static_cast<std::size_t>(nd - first)] += hess[i];
    }
    for (int nd = first; nd <= last; ++nd)
    {
      const std::size_t s = static_cast<std::size_t>(nd - first);
      if (H[s] <= 0.0) { continue; }
      tree.feature[nd] = -1;
      tree.value[nd] = -round_lr_ * leafStep_(G[s], H[s]);
    }
    return any_split || true;   // a depth-0 stump is still a valid (if useless) tree
  }

  GBTParams params_;
  double round_lr_ = 0.0;         ///< the learning rate in force for the tree being grown (== params_.learning_rate unless warmup_rounds > 0)
  std::vector<gbt_detail::Tree> trees_;
  std::vector<std::vector<double>> edges_;
  std::size_t n_features_ = 0;
  double base_ = 0.0;
  /// A node's objective reduction for the split search: G^2/(H+lambda) natively; with clip_gain and a cap,
  /// the reduction the CLIPPED weight realises, 2*G*step - (H+lambda)*step^2 (equal when unclipped).
  double nodeObj_(double G, double H) const
  {
    if (params_.clip_gain && params_.max_delta_step > 0.0)
    {
      const double step = leafStep_(G, H);
      return 2.0 * G * step - (H + params_.lambda) * step * step;
    }
    return G * G / (H + params_.lambda);
  }
  /// The leaf's Newton step before the learning rate, clipped to +/-max_delta_step when that is set.
  double leafStep_(double G, double H) const
  {
    double step = G / (H + params_.lambda);
    if (params_.max_delta_step > 0.0)
    { step = std::max(-params_.max_delta_step, std::min(params_.max_delta_step, step)); }
    return step;
  }
};

} // namespace ODIA::Scoring

#endif // ODIA_GBT_H
