// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Adopted from okohlbacher/OpenDIAlyzer, branch odia-engine-and-openms-boundary,
// tag odia-v0.3.0 (0ed5fb2), file src/odia_lda.h. Same author, same BSD-3 licence.
//
// Changed here: namespace `odia` -> `ODIA::Scoring` to match this project, and
// the sibling includes rewritten to this tree's layout. The algorithms are
// untouched.
//
// Provenance that travels with the code, per that repository's CLEAN-ROOM.md:
// this is an independent reimplementation of a published confidence model
// (Reiter et al. mProphet, Nat Methods 2011; Teleman et al. pyprophet,
// Bioinformatics 2015). No DIA-NN code is present in it.
// odia_lda.h — in-process semi-supervised LDA scorer for OpenDIAlyzer.
//
// A pyprophet/mProphet-style confidence model, reimplemented in-process so
// OpenDIAlyzer needs no external pyprophet. Pure C++ (no OpenMS, no Eigen), so it
// compiles and self-tests standalone. Used in two places:
//   (1) recalibration: pick RT anchors as targets at q < threshold (a real
//       semi-supervised d-score, not raw VAR_XCORR_SHAPE/VAR_LIBRARY_CORR cutoffs);
//   (2) final FDR: q-values on the last pass.
//
// Start with LDA (deterministic, dependency-light). Gradient-boosted trees are the
// backlog upgrade (see docs/OpenDIAlyzer-scoring-fdr-backlog.md).
//
// ALGORITHM (target for implementation):
//   Inputs: features (N rows x M sub-scores), labels (1=target, 0=decoy), group
//   (precursor id; a precursor has several candidate peak groups = rows). FDR is a
//   per-precursor quantity, so it is computed on the best-scoring row per group.
//   1. z-standardize each feature column (mean/sd over all rows).
//   2. k-fold cross-validation BY GROUP (all rows of a precursor fall in one fold),
//      so a row is never scored by a model trained on its own precursor.
//   3. For each fold, train on the other folds by semi-supervised iteration:
//        - initialise the discriminant with the single most target/decoy-separating
//          feature (max |two-sample t|);
//        - pick a confident-target training set: best row per target group whose
//          current score beats the decoy-derived cutoff at a lenient train FDR
//          (e.g. 0.05); all decoy rows are negatives;
//        - Fisher LDA: w = Sw^{-1} (mu_pos - mu_neg), Sw = pooled within-class
//          covariance (solve SPD system by Cholesky; ridge-regularise the diagonal
//          if needed); score = w . x; repeat n_iter times.
//      Apply the fold's final w to its held-out rows -> cross-validated d-score.
//   4. q-values: take the best d-score per precursor (targets and decoys), sort
//      descending; at each cut FDR = #decoys_above / #targets_above (target-decoy),
//      monotonise to q-values; broadcast each precursor's q to its rows.
//   Deterministic given `seed` (only the fold assignment is randomised).
//
#ifndef ODIA_LDA_H
#define ODIA_LDA_H

#include <odia/scoring/gbt.h>
#include <odia/scoring/anchor_training.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <functional>
#include <limits>
#include <locale>
#include <random>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ODIA::Scoring
{

struct ScoredGroups
{
  std::vector<double> dscore;   ///< per input row (peak group)
  std::vector<double> qvalue;   ///< per input row (its precursor's q-value)
  std::vector<double> pvalue;   ///< per input row: tail probability under the decoy null
  std::vector<double> pep;      ///< per input row: local FDR (posterior error probability)
  /// Diagnostics: how many semi-supervised iterations actually FITTED a discriminant, and how many
  /// bailed out for want of confident positives. If trained==0 the returned scores come from the
  /// single-feature initialisation, not from a learned model -- which looks like a working LDA and
  /// is not one.
  int n_iterations_trained = 0;
  int n_iterations_skipped = 0;
};

struct LDAParams
{
  int n_folds = 3;          ///< cross-validation folds (by group)
  int n_iter = 3;           ///< semi-supervised iterations
  double train_fdr_initial = 0.15;  ///< FDR for the FIRST training-set selection. Deliberately more
                            ///< lenient than the rest: the first pass is seeded by a single feature,
                            ///< so a strict cut can select too few positives to fit anything and the
                            ///< model never bootstraps. pyprophet uses 0.15 here for the same reason.
  double train_fdr = 0.05;  ///< FDR for subsequent iterations, once a real discriminant exists.
  double ridge = 1e-6;      ///< diagonal regularisation for the within-class covariance solve

  /// Feature columns whose weight may never be positive, because the caller
  /// knows the feature is lower-is-better. Empty means unconstrained.
  ///
  /// DIA-NN's check_weights does this for its RT-deviation and mass-accuracy
  /// features. Applies to the LDA weight vector only; a tree ensemble has no
  /// weights to clip.
  std::vector<std::size_t> nonpositive_features;
  unsigned seed = 42;       ///< RNG seed (fold assignment only) — determinism
  bool use_pi0 = false;     ///< Storey pi0 correction. false = HONEST/conservative (true 1% FDR,
                            ///< fewer IDs); true = PYPROPHET parity (more IDs, but a nominal 1% is
                            ///< ~2% actual). NOT DIA-NN: 1.7.x uses a plain target-decoy count
                            ///< ratio with NO pi0 term (verified by literature survey 2026-08-15),
                            ///< so citing DIA-NN here was wrong and is the kind of stale comment
                            ///< that leaks into documents.

  /// Draw each decoy's best score from as many candidates as a TARGET has,
  /// rather than from all of its own.
  ///
  /// Target-decoy competition assumes the two classes are exchangeable under
  /// the null. They are not here, and the gap is large: measured on Astral,
  /// targets carry 10.98 candidates per precursor (median 6) and decoys 20.99
  /// (median 24), because the picker's margin rule keeps candidates within
  /// `max_corr_diff` of a precursor's OWN best -- so a precursor with no real
  /// peak admits nearly every position and one with a strong peak admits few.
  ///
  /// A q-value then compares best-of-11 against best-of-21, over 9,453 decoy
  /// precursors. The maximum of that many best-of-21 draws reached 5.997
  /// against a target p99 of 6.517, which is what set the 1% threshold: ODIA
  /// ranked 9,563 Astral precursors that are ALL in DIA-NN's truth set and its
  /// own FDR refused 3,834 of them.
  ///
  /// The cap per decoy group is quantile-matched to the target count
  /// distribution, so the two classes end up drawing from the same
  /// distribution of N. Targets are untouched.
  ///
  /// Which candidates are kept is decided in CANONICAL ORDER (precursor, apex
  /// RT, apex intensity), never by score. Keeping a decoy's top-K by dscore
  /// would subsample exactly the rows that make the null hard and would
  /// manufacture identifications rather than measure them.
  bool match_decoy_candidate_counts = false;
  bool top_decoys_only = true;  ///< Train the negative class on each decoy precursor's BEST peak
                            ///< group only, not on all of its candidate peak groups. See the note at
                            ///< the negative-class construction: the positive class is top-peaks-only
                            ///< by definition, so using all decoy rows makes the two classes differ
                            ///< in peak RANK as well as in label, and the discriminant partly learns
                            ///< "is this the best peak of its group" instead of "is this a target".
                            ///< pyprophet trains on get_top_decoy_peaks() for exactly this reason.
  /// Which learner fills the semi-supervised loop. LDA fits one hyperplane; GBT fits an additive
  /// ensemble of depth-limited trees and can express interactions between sub-scores that no
  /// hyperplane can. Measured on synthetic data with the SAME loop around it: on linearly separable
  /// data LDA 0.981 vs GBT 0.986 AUC (no penalty for the extra capacity), on an interaction LDA
  /// 0.514 -- chance -- vs GBT 0.926. Everything else (folds, training-set selection, fold
  /// normalisation, q-values) is identical, so the two are directly comparable.
  ///
  /// NN adds a third: an ensemble of small tanh MLPs, the shape DIA-NN uses. It can express the
  /// same interactions the GBT can, plus smooth ones a tree approximates in steps.
  enum class Classifier { LDA, GBT, NN };
  Classifier classifier = Classifier::LDA;
  GBTParams gbt;            ///< only consulted when classifier == GBT
  NNParams nn;              ///< only consulted when classifier == NN

  // ---- the anti-circularity mechanisms of odia_anchor_training.h, each toggleable ALONE ---------
  // They are separable on purpose. Bundled, an improvement and a regression cancel and the net
  // result is "no effect"; separately, each one's contribution is attributable.

  /// Mechanism 1. Feature mask for the SEED fit only (empty = no exclusion). The seed selects the
  /// anchors; if it selects them using the very features they will be used to calibrate, the
  /// correction is biased toward the uncorrected state.
  std::vector<char> seed_mask;
  /// Mechanism 3. Share of precursor GROUPS each ensemble member trains on (1.0 = off, which is
  /// DIA-NN's arrangement: members differ only in weight init, so they share one seed bias).
  double bag_fraction = 1.0;
  /// Mechanism 5. Stop when the positive SET stops changing, or when it SHRINKS. Counting
  /// identifications cannot see a collapse -- they rise throughout one.
  bool stop_on_composition = false;
  double stop_jaccard = 0.98;
  /// Collapse/convergence policy forwarded to anchorIterationShouldContinue. The library defaults
  /// are the LEGACY single-shot rule (any shrink = collapse); the CLI registers the repaired
  /// high-water + patience policy (floor 0.01, patience 2) as ITS defaults, so only armed CLI
  /// runs get the new rule and non-CLI callers cannot drift. See AnchorTrainingParams.
  double stop_shrink_floor = 0.0;
  int stop_patience = 1;
  /// Stderr-only churn diagnostic: one line per (fold, iteration) with the positive-set size and
  /// its Jaccard overlap with the previous iteration's set. Never control flow, and it shares no
  /// state with mechanism 5 -- an iteration-sweep arm carries it precisely because it cannot
  /// change the output bytes. Exists because the offline replay (analysis77 R0) measured 11-36%
  /// positive-set churn per iteration with ID counts still rising; whether the REAL loop
  /// converges is what this makes visible.
  bool iteration_log = false;
  /// DIAGNOSTIC ONLY, and FDR-INVALID when true: every group trains the model that scores it.
  ///
  /// It exists because `n_folds = 1` cannot express this -- the fold count is clamped to >= 2 a few
  /// lines into the routine, so an ablation arm that set n_folds=1 silently ran 2-fold CV and
  /// "demonstrated" nothing while claiming to demonstrate that cross-validation is load-bearing.
  /// A flag that says what it does cannot be defeated by a clamp.
  bool disable_cv = false;
  bool normalize_folds = true;  ///< Rescale each fold's held-out scores to its own decoy null
  bool fold_pool_rank = false;  ///< v1.13: pool the folds by WITHIN-FOLD RANK FRACTION over all groups (targets: best
                                ///< over all rows; decoys: best over the drawn prefix, as assignQValues ranks them),
                                ///< DScore = -log10(rank/n_f), instead of (raw - mu_f)/sigma_f. The decoy null of an
                                ///< uncapped fold does not fix the scale of its top (w1_ctl70 fold 1: top compressed
                                ///< 0.75x, bulk sd 1.04x), so the pooled head is partly sorted by fold. mu/sigma are
                                ///< still computed and saved. Off = native, identical arithmetic.
                            ///< (mean 0, sd 1) before pooling. Each fold has its OWN weight vector,
                            ///< with its own arbitrary scale and offset, so the raw scores are not
                            ///< comparable across folds; pooling them into one ranking without this
                            ///< mixes incommensurable scales.
};

/// Semi-supervised LDA scoring with cross-validation and target-decoy q-values.
/// features[i] has the same length for all i (M sub-scores). labels[i] in {0,1}.
/// group[i] is the precursor id shared by that precursor's candidate peak groups.
/// Returns a d-score and q-value per input row. Deterministic given params.seed.
/// `model_out` / `model_in` (gbt engine only) save the trained fold ensemble, or apply a saved
/// one without training; see the definition.
ScoredGroups scoreSemiSupervisedLDA(const std::vector<std::vector<double>>& features,
                                    const std::vector<int>& labels,
                                    const std::vector<long long>& group,
                                    const LDAParams& params = LDAParams(),
                                    const std::string& model_out = "",
                                    const std::string& model_in = "",
                                    const std::vector<std::string>* feature_names = nullptr);

namespace lda_detail
{

struct RankedGroup
{
  std::size_t group_index;
  std::size_t best_row;
  int label;
  double score;
  double qvalue;
  double pvalue = 1.0;   ///< empirical p from the decoy null (TAIL probability at this score)
  double pep = 1.0;      ///< posterior error probability = LOCAL false-discovery rate at this score
};

inline double dot(const std::vector<double>& a, const std::vector<double>& b)
{
  double result = 0.0;
  for (std::size_t j = 0; j < a.size(); ++j) { result += a[j] * b[j]; }
  return result;
}

// Assign target-decoy q-values to a list containing one best score per group.
// Equal scores are treated as one threshold, avoiding order-dependent q-values.
inline void assignQValues(std::vector<RankedGroup>& ranked, bool use_pi0)
{
  std::sort(ranked.begin(), ranked.end(), [](const RankedGroup& a, const RankedGroup& b) {
    if (a.score != b.score) { return a.score > b.score; }
    return a.group_index < b.group_index;
  });

  std::size_t Ntar = 0, Ndec = 0;
  for (const auto& r : ranked) { if (r.label == 1) { ++Ntar; } else { ++Ndec; } }

  // Storey pi0 = estimated fraction of TRUE-NULL targets. Walking high->low score, a
  // target's empirical p-value from the decoy null is (decoys_seen_so_far / Ndec); null
  // targets have ~uniform p, so the mass with p>lambda estimates pi0. Without pi0 the
  // target-decoy estimator assumes pi0=1 (every target could be false) and is
  // over-conservative vs pyprophet/DIA-NN (measured: 22,959 vs 37,539 IDs on the same osw,
  // despite our discriminant separating MORE clean targets). pi0 recovers the honest count.
  double pi0 = 1.0;
  if (use_pi0 && Ntar > 0 && Ndec > 0)
  {
    const double lambda = 0.5;
    const std::size_t dthr = static_cast<std::size_t>(lambda * static_cast<double>(Ndec));
    std::size_t dec_seen = 0, tar_hi = 0;
    for (const auto& r : ranked)
    {
      if (r.label == 0) { ++dec_seen; }
      else if (dec_seen > dthr) { ++tar_hi; }
    }
    pi0 = static_cast<double>(tar_hi) / ((1.0 - lambda) * static_cast<double>(Ntar));
    if (!(pi0 > 0.0)) { pi0 = 1.0 / static_cast<double>(Ntar); }   // never 0
    if (pi0 > 1.0) { pi0 = 1.0; }                                   // never > 1
  }

  std::size_t targets = 0;
  std::size_t decoys = 0;
  for (std::size_t begin = 0; begin < ranked.size();)
  {
    std::size_t end = begin + 1;
    while (end < ranked.size() && ranked[end].score == ranked[begin].score) { ++end; }
    for (std::size_t i = begin; i < end; ++i)
    {
      if (ranked[i].label == 1) { ++targets; }
      else                      { ++decoys; }
    }
    // FDR(t) = pi0 * (decoys_above/Ndec) / (targets_above/Ntar)  -- p-value based, ratio-corrected.
    double fdr;
    if (targets == 0) { fdr = std::numeric_limits<double>::infinity(); }
    else
    {
      // +1 finite-sample correction on the decoy count (Käll): keeps the estimator
      // honest at the tail where single decoys otherwise make it anti-conservative.
      const double dr = (Ndec > 0) ? (static_cast<double>(decoys) + 1.0) / static_cast<double>(Ndec) : 0.0;
      const double tr = static_cast<double>(targets) / static_cast<double>(Ntar);
      fdr = std::min(1.0, pi0 * dr / tr);
    }
    for (std::size_t i = begin; i < end; ++i) { ranked[i].qvalue = fdr; }
    begin = end;
  }

  double running_min = 1.0;
  for (std::size_t end = ranked.size(); end > 0;)
  {
    std::size_t begin = end - 1;
    while (begin > 0 && ranked[begin - 1].score == ranked[end - 1].score) { --begin; }
    running_min = std::min(running_min, ranked[begin].qvalue);
    for (std::size_t i = begin; i < end; ++i) { ranked[i].qvalue = running_min; }
    end = begin;
  }

  // --- p-value and PEP -------------------------------------------------------------------
  // These were previously never computed: the caller wrote the q-value into the PVALUE, QVALUE and
  // PEP columns alike, so 100% of rows had all three identical and anything downstream reading PEP
  // (IPF, protein-level inference) was silently consuming a q-value. They are different statistics:
  //
  //   p-value : TAIL probability under the null -- P(a null score >= this one), from the decoys.
  //   q-value : minimum FDR of the SET selected at this threshold (already computed above).
  //   PEP     : LOCAL FDR -- the probability that THIS peak group specifically is null. A group at
  //             q = 0.01 sitting right at the threshold can easily have PEP ~ 0.3.
  //
  // p is the standard conservative empirical estimate (Käll's +1 on both counts).
  {
    std::size_t dec_at_or_above = 0;
    for (std::size_t begin = 0; begin < ranked.size();)
    {
      std::size_t end = begin + 1;
      while (end < ranked.size() && ranked[end].score == ranked[begin].score) { ++end; }
      for (std::size_t i = begin; i < end; ++i) { if (ranked[i].label == 0) { ++dec_at_or_above; } }
      const double p = (Ndec > 0)
                         ? (static_cast<double>(dec_at_or_above) + 1.0) / (static_cast<double>(Ndec) + 1.0)
                         : 1.0;
      for (std::size_t i = begin; i < end; ++i) { ranked[i].pvalue = std::min(1.0, p); }
      begin = end;
    }
  }

  // PEP by the LOCAL analogue of the global estimator used above: in a score neighbourhood holding
  // t targets and d decoys, the decoys estimate the null target density (scaled by Ntar/Ndec and
  // pi0), so the expected null-target count is pi0*d*(Ntar/Ndec) and PEP ~ that divided by t.
  // Window is a fixed fraction of the list so it adapts to size; a raw local ratio is very noisy,
  // hence the monotonicity pass afterwards.
  {
    const std::size_t n = ranked.size();
    const std::size_t half = std::max<std::size_t>(50, n / 200);   // ~0.5% of the list, >=100 wide
    const double scale = (Ndec > 0) ? (static_cast<double>(Ntar) / static_cast<double>(Ndec)) : 0.0;
    // SLIDING counts, not a recount per element. The window is n/200 wide, so recomputing it inside
    // the loop is O(n^2/100): at n = 2M peak groups that is ~4e10 operations and turns a seconds-long
    // step into a multi-minute one. Maintaining running counts makes the whole pass O(n) -- each
    // element enters and leaves the window exactly once.
    std::size_t t = 0, d = 0, lo = 0, hi = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
      const std::size_t new_lo = (i > half) ? i - half : 0;
      const std::size_t new_hi = std::min(n, i + half + 1);
      while (hi < new_hi) { if (ranked[hi].label == 1) { ++t; } else { ++d; } ++hi; }
      while (lo < new_lo) { if (ranked[lo].label == 1) { --t; } else { --d; } ++lo; }
      double pep = 1.0;
      if (t > 0) { pep = pi0 * static_cast<double>(d) * scale / static_cast<double>(t); }
      ranked[i].pep = std::min(1.0, std::max(0.0, pep));
    }
    // PEP must not increase with score. The list is sorted high->low by score, so PEP must be
    // NON-DECREASING IN INDEX, and the sweep therefore runs from the high-score end forwards.
    //
    // It ran backwards. Sweeping from the low-score end with a running max gave every entry the
    // maximum PEP of its SUFFIX -- i.e. of everything scoring below it -- which is monotone in
    // the wrong direction and, on a real list, catastrophic: [0.01, 0.02, 0.5, 0.9] came out
    // [0.9, 0.9, 0.9, 0.9], so the best-scoring precursor in the run was reported with the worst
    // PEP in the run. Found by external review.
    //
    // Confined to the reported PEP column: nothing thresholds on it (it reaches only the output
    // TSV via OpenDIAlyzer.cpp), so no q-value or identification count was affected.
    //
    // The old comment also called this "the isotonic projection". A cumulative maximum is the
    // greatest monotone minorant, not a least-squares isotonic regression; the claim is dropped
    // rather than repeated.
    double running_max = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
      running_max = std::max(running_max, ranked[i].pep);
      ranked[i].pep = running_max;
    }
  }
}

// Solve A x = b for symmetric positive-definite A. The input matrix is
// row-major and is replaced by its lower-triangular Cholesky factor.
inline bool choleskySolve(std::vector<double> a,
                          const std::vector<double>& b,
                          std::vector<double>& x)
{
  const std::size_t m = b.size();
  for (std::size_t i = 0; i < m; ++i)
  {
    for (std::size_t j = 0; j <= i; ++j)
    {
      double value = a[i * m + j];
      for (std::size_t k = 0; k < j; ++k) { value -= a[i * m + k] * a[j * m + k]; }
      if (i == j)
      {
        if (!(value > 0.0) || !std::isfinite(value)) { return false; }
        a[i * m + i] = std::sqrt(value);
      }
      else
      {
        a[i * m + j] = value / a[j * m + j];
      }
    }
  }

  std::vector<double> y(m, 0.0);
  for (std::size_t i = 0; i < m; ++i)
  {
    double value = b[i];
    for (std::size_t j = 0; j < i; ++j) { value -= a[i * m + j] * y[j]; }
    y[i] = value / a[i * m + i];
  }

  x.assign(m, 0.0);
  for (std::size_t ii = m; ii > 0; --ii)
  {
    const std::size_t i = ii - 1;
    double value = y[i];
    for (std::size_t j = i + 1; j < m; ++j) { value -= a[j * m + i] * x[j]; }
    x[i] = value / a[i * m + i];
    if (!std::isfinite(x[i])) { return false; }
  }
  return true;
}

} // namespace lda_detail

/// @param model_out  gbt engine only: after training, save the fold models, the
///                   standardisation and each fold's decoy normalisation here.
/// @param model_in   gbt engine only: do not train -- load that file and score every
///                   row with it. The two paths mirror `scorePercolator`'s.
inline ScoredGroups scoreSemiSupervisedLDA(
  const std::vector<std::vector<double>>& features,
  const std::vector<int>& labels,
  const std::vector<long long>& group,
  const LDAParams& params,
  const std::string& model_out,
  const std::string& model_in,
  const std::vector<std::string>* feature_names)
{
  const std::size_t n = features.size();
  ScoredGroups result;
  result.dscore.assign(n, 0.0);
  result.qvalue.assign(n, 1.0);
  result.pvalue.assign(n, 1.0);
  result.pep.assign(n, 1.0);
  if (n == 0 || labels.size() != n || group.size() != n) { return result; }

  const std::size_t m = features.front().size();
  if (m == 0) { return result; }
  for (const auto& row : features)
  {
    if (row.size() != m) { return result; }
  }

  // Global z-standardisation is unsupervised. Non-finite cells are treated as
  // missing and become zero (the column mean) after standardisation.
  std::vector<double> mean(m, 0.0);
  std::vector<double> count(m, 0.0);
  for (const auto& row : features)
  {
    for (std::size_t j = 0; j < m; ++j)
    {
      if (std::isfinite(row[j]))
      {
        mean[j] += row[j];
        count[j] += 1.0;
      }
    }
  }
  for (std::size_t j = 0; j < m; ++j)
  {
    if (count[j] > 0.0) { mean[j] /= count[j]; }
  }

  std::vector<double> sum_squared_deviation(m, 0.0);
  std::vector<double> sd(m, 1.0);
  std::vector<std::vector<double>> z(n, std::vector<double>(m, 0.0));
  for (const auto& row : features)
  {
    for (std::size_t j = 0; j < m; ++j)
    {
      if (std::isfinite(row[j]))
      {
        const double d = row[j] - mean[j];
        sum_squared_deviation[j] += d * d;
      }
    }
  }
  for (std::size_t j = 0; j < m; ++j)
  {
    sd[j] = count[j] > 1.0
              ? std::sqrt(sum_squared_deviation[j] / (count[j] - 1.0))
              : 0.0;
    if (!(sd[j] > std::numeric_limits<double>::epsilon()) || !std::isfinite(sd[j]))
    {
      sd[j] = 1.0;
    }
  }
  for (std::size_t i = 0; i < n; ++i)
  {
    for (std::size_t j = 0; j < m; ++j)
    {
      z[i][j] = std::isfinite(features[i][j]) ? (features[i][j] - mean[j]) / sd[j] : 0.0;
    }
  }

  // Build a stable, first-occurrence ordering of precursor groups.
  std::unordered_map<long long, std::size_t> group_lookup;
  group_lookup.reserve(n);
  std::vector<std::vector<std::size_t>> group_rows;
  std::vector<int> group_label;
  std::vector<long long> group_id;                 // the precursor id behind each group index
  for (std::size_t i = 0; i < n; ++i)
  {
    auto inserted = group_lookup.emplace(group[i], group_rows.size());
    if (inserted.second)
    {
      group_rows.emplace_back();
      group_label.push_back(labels[i] == 1 ? 1 : 0);
      group_id.push_back(group[i]);
    }
    const std::size_t g = inserted.first->second;
    group_rows[g].push_back(i);
  }

  const std::size_t group_count = group_rows.size();
  if (group_count < 2) { return result; } // leakage-free training is impossible

  int folds = params.n_folds;
  if (folds < 2) { folds = 2; }
  if (folds > static_cast<int>(group_count)) { folds = static_cast<int>(group_count); }

  std::vector<std::size_t> target_groups;
  std::vector<std::size_t> decoy_groups;
  for (std::size_t g = 0; g < group_count; ++g)
  {
    (group_label[g] == 1 ? target_groups : decoy_groups).push_back(g);
  }
  // Order these by the precursor's IDENTITY before shuffling. A group's index is its
  // first-occurrence position in row order, and row order is whatever the parallel extraction
  // happened to produce -- so a seeded shuffle of indices still put the same precursor in a
  // different fold on every run, training a different model and reporting a different ID count.
  // Measured spread on byte-identical input: 6487 / 6565 / 6433 (+/-1%), which is larger than
  // most of the effects being A/B tested. Sorting by id first costs one sort and makes the fold
  // assignment a function of the data alone; the shuffle still balances fold sizes exactly.
  const auto by_id = [&](std::size_t a, std::size_t b) { return group_id[a] < group_id[b]; };
  std::sort(target_groups.begin(), target_groups.end(), by_id);
  std::sort(decoy_groups.begin(), decoy_groups.end(), by_id);
  std::mt19937 rng(params.seed);
  std::shuffle(target_groups.begin(), target_groups.end(), rng);
  std::shuffle(decoy_groups.begin(), decoy_groups.end(), rng);
  std::vector<int> group_fold(group_count, 0);
  for (std::size_t i = 0; i < target_groups.size(); ++i)
  {
    group_fold[target_groups[i]] =
      static_cast<int>(i % static_cast<std::size_t>(folds));
  }
  for (std::size_t i = 0; i < decoy_groups.size(); ++i)
  {
    group_fold[decoy_groups[i]] =
      static_cast<int>(i % static_cast<std::size_t>(folds));
  }

  // Folds are independent BY CONSTRUCTION: fold f trains on the groups not assigned to f and writes
  // result.dscore only for the rows of groups that ARE assigned to f. So no two iterations read or
  // write the same dscore element, and none of them touch shared state except the two skip/train
  // counters, which are reduced. Everything else (`z`, `group_rows`, `group_fold`) is read-only here.
  // This loop was serial, and on the benchmark feature table the whole scoring step is 142 s of the
  // 1,284 s serial tail that caps the run's speedup at 1.90x -- see
  // docs/OpenDIAlyzer-parallel-efficiency.md. Cheap to fix, so fixed.
  //
  // Determinism is preserved: fold assignment comes from the seeded RNG above, each fold's model
  // depends only on its own training set, and the q-values are computed after the loop from the
  // pooled scores. The result does not depend on completion order.
  // The GBT parallelises internally over rows, but fit() is called from INSIDE the fold loop below.
  // A nested parallel region defaults to a team of one, so without raising the active-level limit
  // the inner pragmas would be dead code. Splitting the available threads as
  // folds x (threads/folds) uses the machine without oversubscribing it. The GBT's result does not
  // depend on this number (odia_gbt_test T8), so it is purely a speed knob.
  GBTParams gbt_params = params.gbt;
  NNParams nn_params = params.nn;
  // Threads for the NESTED regions inside the fold loop. Computed once and applied to every one of
  // them: the training loop already honoured it, but the per-group scans did not, so each of the 3
  // concurrent folds opened teams of the FULL thread count -- 540 threads on 224 cores at
  // OMP_NUM_THREADS=180. Oversubscription of that size costs more in scheduling than the
  // parallelism returns.
  int inner_threads = 1;
#ifdef _OPENMP
  // BOTH inner-parallel learners need the active-level raised, not just the GBT. With it set for
  // the GBT alone, the network's chunk loop ran as a team of ONE inside the fold loop: measured
  // 435% CPU on a 224-core node, i.e. the 3 folds and nothing else. Splitting as
  // folds x (threads/folds) uses the machine without oversubscribing it, and neither learner's
  // result depends on the number (odia_gbt_test T8; odia_nn_test thread invariance).
  if (params.classifier == LDAParams::Classifier::GBT || params.classifier == LDAParams::Classifier::NN)
  {
    omp_set_max_active_levels(2);
    inner_threads = std::max(1, omp_get_max_threads() / std::max(1, folds));
    gbt_params.n_threads = inner_threads;
    nn_params.n_threads = inner_threads;
  }
#endif

  int n_trained = 0, n_skipped = 0;

  // ---- FROZEN MODEL, gbt engine only ----
  //
  // Score every row with a fold ensemble trained elsewhere and saved by the
  // `model_out` path below: THAT run's standardisation, its fold models, and
  // each fold's decoy normalisation, applied here with no training at all.
  // Measured 2026-09-05 (pick/wf_hist.txt, wf_fixedmodel.txt): the native
  // retraining flips between a compact and a saturated score regime under
  // small changes to the binary or the feature distribution, and every
  // comparison across that flip is about the scorer, not the evidence; the
  // same evidence under ONE model held fixed across arms is stable. A frozen
  // model is how two runs are compared on evidence, and it is the documented
  // remedy for a run that cannot bootstrap its own positives.
  //
  // The per-row score is the mean over folds of each fold's DECOY-normalised
  // score, so it sits on the same "decoy sds above the null" scale the native
  // path pools folds on. Rows the saving run never saw are, by construction,
  // out of sample for every fold model.
  const bool gbt_engine = (params.classifier == LDAParams::Classifier::GBT);
  std::vector<GBT> fold_models;
  std::vector<double> fold_norm_mu, fold_norm_sigma;
  bool frozen_applied = false;
  if (!model_in.empty() && !gbt_engine)
  {
    std::fprintf(stderr, "[lda] -classifier_model_in is honoured by the gbt and percolator "
                         "engines only; this engine trains as usual\n");
  }
  if (!model_in.empty() && !model_out.empty() && gbt_engine)
  {
    std::fprintf(stderr, "[gbt] both -classifier_model_in and -classifier_model_out are set: "
                         "the frozen model is APPLIED and nothing is saved\n");
  }
  if (!model_in.empty() && gbt_engine)
  {
    // A frozen model that cannot be applied is a HARD failure, never a silent
    // fall-back to training: the retraining is the very thing this path exists
    // to hold fixed, and a run that quietly trained would enter a comparison
    // as if it had not. The failure surfaces as "fitted 0 iterations", exactly
    // as a failed percolator engine does.
    auto refuse = [&](const std::string& why) -> ScoredGroups {
      std::fprintf(stderr, "[gbt] FROZEN MODEL %s: REFUSED -- %s; no scores produced\n",
                   model_in.c_str(), why.c_str());
      result.n_iterations_skipped = 1;
      return result;
    };
    if (params.fold_pool_rank)
    { return refuse("-fold_pool_rank is not supported with -classifier_model_in (v1.13): the pooled scale is a within-fold rank, not the saved mu/sigma"); }
    std::ifstream is(model_in);
    is.imbue(std::locale::classic());
    if (!is) { return refuse("cannot open"); }
    std::string magic;
    int version = 0;
    std::size_t m_saved = 0, k_saved = 0;
    if (!(is >> magic >> version >> m_saved >> k_saved) || magic != "ODIA-GBT-FOLDS" || version != 1)
    { return refuse("not an ODIA-GBT-FOLDS version-1 file"); }
    if (m_saved != m)
    {
      return refuse("saved for " + std::to_string(m_saved) + " sub-scores, this run has " +
                    std::to_string(m));
    }
    if (k_saved == 0 || k_saved > 64) { return refuse("implausible fold count"); }
    // The sub-score NAMES travel with the model and must match exactly: a
    // different sub-score set of the same width would load and rank on
    // nonsense with nothing downstream the wiser.
    std::size_t n_names = 0;
    if (!(is >> n_names) || n_names != m) { return refuse("sub-score name list missing or wrong length"); }
    for (std::size_t j = 0; j < m; ++j)
    {
      std::string nm;
      if (!(is >> nm)) { return refuse("sub-score name list truncated"); }
      if (feature_names != nullptr && j < feature_names->size() && (*feature_names)[j] != nm)
      {
        return refuse("sub-score " + std::to_string(j) + " is '" + nm + "' in the model but '" +
                      (*feature_names)[j] + "' in this run");
      }
    }
    std::vector<double> mean_saved(m, 0.0), sd_saved(m, 1.0);
    std::vector<GBT> models(k_saved);
    std::vector<double> mus(k_saved, 0.0), sigmas(k_saved, 0.0);
    std::vector<char> present(k_saved, 0);
    for (std::size_t j = 0; j < m; ++j)
    { if (!(is >> mean_saved[j]) || !std::isfinite(mean_saved[j])) { return refuse("standardisation means truncated or non-finite"); } }
    for (std::size_t j = 0; j < m; ++j)
    {
      if (!(is >> sd_saved[j]) || !std::isfinite(sd_saved[j]) || !(sd_saved[j] > 0.0))
      { return refuse("standardisation sds truncated, non-finite or non-positive"); }
    }
    for (std::size_t f = 0; f < k_saved; ++f)
    {
      int has = 0;
      if (!(is >> has >> mus[f] >> sigmas[f])) { return refuse("fold header truncated"); }
      if (has)
      {
        if (!models[f].load(is)) { return refuse("fold " + std::to_string(f) + " model malformed"); }
        if (models[f].featureCount() != m) { return refuse("fold model feature count mismatch"); }
        present[f] = 1;
      }
    }
    std::size_t k_have = 0;
    for (const char p : present) { k_have += (p != 0); }
    if (k_have == 0) { return refuse("no trained fold in the file"); }
    // Optional (files written before v1.7 lack it): the saving run's precursor-
    // to-fold assignment. A precursor listed here was TRAINING data for every
    // fold but its own, and is scored below by that one excluded fold only.
    std::unordered_map<long long, int> fold_of;
    {
      std::string tag;
      std::size_t n_map = 0;
      if (is >> tag)
      {
        if (tag != "FOLDMAP" || !(is >> n_map) || n_map > 50000000)
        { return refuse("trailing content is not a valid FOLDMAP"); }
        fold_of.reserve(n_map * 2);
        for (std::size_t g = 0; g < n_map; ++g)
        {
          long long pid = 0;
          int fd = -1;
          if (!(is >> pid >> fd)) { return refuse("FOLDMAP truncated"); }
          if (fd < 0 || static_cast<std::size_t>(fd) >= k_saved) { return refuse("FOLDMAP fold index out of range"); }
          fold_of.emplace(pid, fd);
        }
      }
    }
    std::size_t scored_exact = 0;
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (std::size_t i = 0; i < n; ++i)
      {
        std::vector<double> zf(m, 0.0);
        for (std::size_t j = 0; j < m; ++j)
        {
          zf[j] = std::isfinite(features[i][j]) ? (features[i][j] - mean_saved[j]) / sd_saved[j] : 0.0;
        }
        auto fold_score = [&](std::size_t f) -> double {
          double v = models[f].score(zf);
          if (sigmas[f] > std::numeric_limits<double>::epsilon() && std::isfinite(sigmas[f]))
          { v = (v - mus[f]) / sigmas[f]; }
          return v;
        };
        const auto it = fold_of.find(group[i]);
        if (it != fold_of.end() && present[static_cast<std::size_t>(it->second)])
        {
          // Training data of every other fold: its own excluded fold only.
          result.dscore[i] = fold_score(static_cast<std::size_t>(it->second));
#ifdef _OPENMP
#pragma omp atomic
#endif
          ++scored_exact;
        }
        else
        {
          double s = 0.0;
          for (std::size_t f = 0; f < k_saved; ++f)
          {
            if (!present[f]) { continue; }
            s += fold_score(f);
          }
          result.dscore[i] = s / static_cast<double>(k_have);
        }
      }
      frozen_applied = true;
      // Reported as trained iterations, as scorePercolator reports its frozen
      // path: downstream derives fdr_valid from `> 0`, and a frozen model IS a
      // valid discriminant. It is no longer "this run bootstrapped its own
      // positives", which is what the name says -- the only consumer today is
      // the `> 0` test, so the conflation is documented here rather than fixed.
      n_trained = static_cast<int>(k_have);
      std::fprintf(stderr, "[gbt] FROZEN MODEL %s: %zu fold models over %zu features applied to "
                           "%zu rows (%zu rows of precursors the model trained on scored by their "
                           "own excluded fold, the rest by the fold mean), no training. NOTE: "
                           "q-values under a frozen model are an FDR "
                           "estimate only if the saving run shared no targets with this one; when "
                           "it did (same sample, same library) they are anti-conservative -- use "
                           "them to COMPARE runs, not to report an FDR\n",
                   model_in.c_str(), k_have, m, n, scored_exact);
    }
  }
  // The captures below are sized only when a save was requested, so a run with
  // neither path set allocates and copies nothing -- flag-off touches one bool.
  if (!frozen_applied && gbt_engine && !model_out.empty())
  {
    fold_models.resize(static_cast<std::size_t>(folds));
    fold_norm_mu.assign(static_cast<std::size_t>(folds), 0.0);
    fold_norm_sigma.assign(static_cast<std::size_t>(folds), 0.0);
  }

  // How many candidates each group may draw its best from. Targets always use
  // all of theirs; decoys are quantile-matched to the target distribution when
  // asked. See `match_decoy_candidate_counts`.
  std::vector<std::size_t> draw_from(group_count);
  for (std::size_t g = 0; g < group_count; ++g) { draw_from[g] = group_rows[g].size(); }
  if (params.match_decoy_candidate_counts)
  {
    std::vector<std::size_t> target_n;
    std::vector<std::size_t> decoys;
    for (std::size_t g = 0; g < group_count; ++g)
    {
      if (group_label[g] == 1) { target_n.push_back(group_rows[g].size()); }
      else { decoys.push_back(g); }
    }
    if (!target_n.empty() && !decoys.empty())
    {
      std::sort(target_n.begin(), target_n.end());
      // Both sides sorted by count, then matched by rank: the k-th smallest
      // decoy takes the k-th smallest target's count. That maps the whole
      // distribution rather than just its mean, and it is a function of the
      // data alone, so the run stays deterministic.
      std::stable_sort(decoys.begin(), decoys.end(),
                       [&](std::size_t a, std::size_t b)
                       { return group_rows[a].size() < group_rows[b].size(); });
      for (std::size_t i = 0; i < decoys.size(); ++i)
      {
        const std::size_t q = i * target_n.size() / decoys.size();
        draw_from[decoys[i]] = std::min(target_n[q], group_rows[decoys[i]].size());
        if (draw_from[decoys[i]] == 0) { draw_from[decoys[i]] = 1; }
      }
    }
  }


  if (!frozen_applied)
  {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : n_trained, n_skipped)
#endif
  for (int fold = 0; fold < folds; ++fold)
  {
    std::vector<std::size_t> train_rows;
    std::vector<std::size_t> train_groups;
    train_rows.reserve(n);
    train_groups.reserve(group_count);
    for (std::size_t g = 0; g < group_count; ++g)
    {
      if (!params.disable_cv && group_fold[g] == fold) { continue; }
      train_groups.push_back(g);
      train_rows.insert(train_rows.end(), group_rows[g].begin(), group_rows[g].end());
    }
    if (train_rows.empty()) { continue; }

    // w = S_W^-1 (mu_pos - mu_neg): Fisher's discriminant on two explicit row sets.
    auto fit_lda = [&](const std::vector<std::size_t>& positive_rows,
                       const std::vector<std::size_t>& negative_rows,
                       std::vector<double>& out) -> bool {
      if (positive_rows.size() < 2 || negative_rows.size() < 2) { return false; }
      std::vector<double> positive_mean(m, 0.0);
      std::vector<double> negative_mean(m, 0.0);
      for (const std::size_t row : positive_rows)
      {
        for (std::size_t j = 0; j < m; ++j) { positive_mean[j] += z[row][j]; }
      }
      for (const std::size_t row : negative_rows)
      {
        for (std::size_t j = 0; j < m; ++j) { negative_mean[j] += z[row][j]; }
      }
      for (std::size_t j = 0; j < m; ++j)
      {
        positive_mean[j] /= static_cast<double>(positive_rows.size());
        negative_mean[j] /= static_cast<double>(negative_rows.size());
      }

      std::vector<double> covariance(m * m, 0.0);
      auto add_within_class = [&](const std::vector<std::size_t>& rows,
                                  const std::vector<double>& class_mean) {
        for (const std::size_t row : rows)
        {
          for (std::size_t j = 0; j < m; ++j)
          {
            const double dj = z[row][j] - class_mean[j];
            for (std::size_t k = 0; k <= j; ++k)
            {
              covariance[j * m + k] += dj * (z[row][k] - class_mean[k]);
            }
          }
        }
      };
      add_within_class(positive_rows, positive_mean);
      add_within_class(negative_rows, negative_mean);
      const double dof =
        static_cast<double>(positive_rows.size() + negative_rows.size() - 2);
      for (std::size_t j = 0; j < m; ++j)
      {
        for (std::size_t k = 0; k <= j; ++k)
        {
          covariance[j * m + k] /= dof;
          covariance[k * m + j] = covariance[j * m + k];
        }
      }

      std::vector<double> difference(m);
      for (std::size_t j = 0; j < m; ++j)
      {
        difference[j] = positive_mean[j] - negative_mean[j];
      }

      // A positive user ridge is tried exactly; a zero/negative ridge starts
      // with a tiny numerical ridge. Increase it geometrically on failure.
      double ridge = params.ridge > 0.0 && std::isfinite(params.ridge)
                       ? params.ridge
                       : 1e-12;
      bool solved = false;
      for (int attempt = 0; attempt < 10 && !solved; ++attempt)
      {
        std::vector<double> regularised = covariance;
        for (std::size_t j = 0; j < m; ++j) { regularised[j * m + j] += ridge; }
        solved = lda_detail::choleskySolve(regularised, difference, out);
        ridge *= 10.0;
      }
      // A feature the caller declares lower-is-better may never earn a positive
      // weight, however the fold's data happens to fall.
      //
      // This is DIA-NN's check_weights (diann.cpp:6592), which hard-clips the
      // RT-deviation and mass-accuracy weights to <= 0 so that a larger
      // deviation can never raise a score. The principle is the same wherever a
      // feature has a known direction: without it the discriminant is free to
      // learn, in-sample, that being further from the prediction is evidence
      // FOR a peptide -- which fits the fold and generalises to nothing.
      //
      // Only the linear classifier is constrained here. A tree ensemble has no
      // weight vector to clip, and enforcing the same thing on GBT needs
      // monotone split constraints, which is a larger change. Since GBT is the
      // default, this currently protects the non-default path.
      for (const std::size_t j : params.nonpositive_features)
      {
        if (j < out.size() && out[j] > 0.0) { out[j] = 0.0; }
      }

      double norm_squared = 0.0;
      for (const double value : out) { norm_squared += value * value; }
      return solved && norm_squared > std::numeric_limits<double>::epsilon();
    };

    // The learner for this fold. LDA carries a weight vector; GBT carries a tree ensemble. Every
    // score in this fold goes through score_row(), so the choice is made in exactly ONE place and
    // the surrounding machinery -- training-set selection, fold normalisation, q-values -- is
    // literally the same code for both.
    GBT gbt;
    std::vector<NNEnsemble> nn_members;      // one per bag; size 1 when bagging is off
    const bool use_gbt = (params.classifier == LDAParams::Classifier::GBT);
    const bool use_nn = (params.classifier == LDAParams::Classifier::NN);


    // Initial direction: the signed feature with the largest absolute Welch
    // two-sample t statistic between all target and decoy training rows.
    std::vector<double> w(m, 0.0);
    auto score_row = [&](std::size_t row) -> double {
      if (use_nn && !nn_members.empty()) { return baggedScore(nn_members, z[row]); }
      return (use_gbt && gbt.trained()) ? gbt.score(z[row]) : lda_detail::dot(w, z[row]);
    };
    // Fit whichever learner this run selected, on the given rows. `mask` is honoured by the NN
    // only -- for a tree, a feature the fit never split on is already inert at scoring time, and
    // for LDA a zero weight is likewise inert, so neither needs one.
    auto fit_learner = [&](const std::vector<std::size_t>& pos, const std::vector<std::size_t>& neg,
                           const std::vector<char>& mask) -> bool {
      if (use_nn)
      {
        NNParams np = nn_params;
        np.mask = mask;
        if (params.bag_fraction >= 1.0)
        {
          NNEnsemble e;
          if (!e.fit(z, pos, neg, np)) { return false; }
          nn_members.assign(1, std::move(e));
          return true;
        }
        if (!(params.bag_fraction > 0.0)) { return false; }   // an empty bag is not a trained model
        AnchorTrainingParams ap;
        ap.nn = np;
        ap.bag_fraction = params.bag_fraction;
        ap.min_anchors = 1;                  // the outer loop already refuses to fit on nothing
        ap.seed = params.seed;
        auto members = trainBaggedOnAnchors(z, pos, neg, group, ap);
        if (members.empty()) { return false; }
        nn_members = std::move(members);
        return true;
      }
      if (use_gbt) { GBT g; if (!g.fit(z, pos, neg, gbt_params)) { return false; }
                     gbt = std::move(g); return true; }
      std::vector<double> next_w;
      if (!fit_lda(pos, neg, next_w)) { return false; }
      w.swap(next_w);
      return true;
    };

    // Best-scoring row of each training group, computed in parallel.
    //
    // All three per-group scans below (seed selection, ranking, negative selection) have this
    // shape, and for the NN each score_row() is a forward pass through TWELVE nets -- which made
    // them the serial tail that held the whole routine to 763% CPU on a 224-core node while
    // training itself was parallel.
    //
    // Output is PRE-SIZED and written by position, so nothing is appended and no ordering question
    // arises: out[i] is group train_groups[i]'s best row whatever order the iterations complete in.
    // The tie-break is `>` -- the FIRST maximum in the group's own row order wins -- which is what
    // the serial loops did, so the result is identical, not merely equivalent.
    auto best_rows_of = [&](const std::vector<std::size_t>& groups,
                            std::vector<std::size_t>& out, std::vector<double>& out_score) {
      out.assign(groups.size(), 0);
      out_score.assign(groups.size(), 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(inner_threads) if (groups.size() > 512)
#endif
      for (long long i = 0; i < static_cast<long long>(groups.size()); ++i)
      {
        const auto& rows_of_g = group_rows[groups[static_cast<std::size_t>(i)]];
        std::size_t best_row = rows_of_g.front();
        double best = score_row(best_row);
        for (const std::size_t row : rows_of_g)
        {
          const double sc = score_row(row);
          if (sc > best) { best = sc; best_row = row; }
        }
        out[static_cast<std::size_t>(i)] = best_row;
        out_score[static_cast<std::size_t>(i)] = best;
      }
    };
    std::vector<std::size_t> bg_row;
    std::vector<double> bg_score;

    std::size_t best_feature = 0;
    double best_abs_t = -1.0;
    double best_difference = 1.0;
    for (std::size_t j = 0; j < m; ++j)
    {
      // The seed mask applies HERE TOO, not only to the fit. This bootstrap picks one feature, and
      // that single-feature direction is what ranks each group to choose its seed row. Searching a
      // masked column would let the excluded feature choose the training examples, and zeroing its
      // value at fit time cannot undo which rows it selected -- the leak mechanism 1 exists to
      // close, reopened one step upstream. Found by adversarial review, not by a test.
      if (!params.seed_mask.empty() && j < params.seed_mask.size() && !params.seed_mask[j]) { continue; }
      double sum[2] = {0.0, 0.0};
      double sum_sq[2] = {0.0, 0.0};
      std::size_t class_n[2] = {0, 0};
      for (const std::size_t row : train_rows)
      {
        const int cls = labels[row] == 1 ? 1 : 0;
        sum[cls] += z[row][j];
        sum_sq[cls] += z[row][j] * z[row][j];
        ++class_n[cls];
      }
      if (class_n[0] == 0 || class_n[1] == 0) { continue; }
      const double class_mean0 = sum[0] / static_cast<double>(class_n[0]);
      const double class_mean1 = sum[1] / static_cast<double>(class_n[1]);
      const double var0 = class_n[0] > 1
                            ? std::max(0.0, (sum_sq[0] - sum[0] * class_mean0) /
                                                static_cast<double>(class_n[0] - 1))
                            : 0.0;
      const double var1 = class_n[1] > 1
                            ? std::max(0.0, (sum_sq[1] - sum[1] * class_mean1) /
                                                static_cast<double>(class_n[1] - 1))
                            : 0.0;
      const double difference = class_mean1 - class_mean0;
      const double standard_error =
        std::sqrt(var0 / static_cast<double>(class_n[0]) +
                  var1 / static_cast<double>(class_n[1]));
      const double abs_t = standard_error > 0.0
                             ? std::abs(difference) / standard_error
                             : (difference == 0.0 ? 0.0
                                                  : std::numeric_limits<double>::infinity());
      if (abs_t > best_abs_t)
      {
        best_abs_t = abs_t;
        best_feature = j;
        best_difference = difference;
      }
    }
    w[best_feature] = best_difference < 0.0 ? -1.0 : 1.0;

    // ...but a SINGLE feature is a far weaker seed than it looks, and when it is too weak the
    // semi-supervised loop never ignites: iteration 0 selects positives with this w, finds none at
    // q<=train_fdr_initial, skips the fit, so w is unchanged and every later iteration skips too.
    // Measured on testdata/lda_fixture.txt (44,568 real OpenSWATH rows / 2,500 precursors): the best
    // single feature reaches |t|=13.7, which sounds decisive but over 29,482 rows is a 0.16 sd
    // per-row effect -- and with ~18 candidate peak groups per precursor, the max-over-group that
    // actually does the ranking is dominated by extreme-value noise. Result: a FLAT q of 0.965 for
    // every target group, 0/9 iterations trained, and 0 IDs.
    //
    // The fix is to seed with a real multivariate direction. The target/decoy labels are known
    // outright -- no FDR estimate needed -- so an LDA of all top-target rows against all top-decoy
    // rows is available for free and is enormously stronger than one column. The target class is
    // contaminated (most target precursors are false), which shrinks the fitted direction toward
    // zero but does not rotate it systematically: the contaminant IS the decoy distribution, so it
    // biases mu_pos toward mu_neg and costs magnitude, not orientation. It only has to be good
    // enough to ignite the loop; iteration 0 then re-selects positives properly.
    // pyprophet has no equivalent problem because OpenSWATH hands it a composite `main_score` to
    // rank by (semi_supervised.py: train.rank_by("main_score")). ODIA has no such column.
    {
      std::vector<std::size_t> seed_pos, seed_neg;
      seed_pos.reserve(train_groups.size());
      seed_neg.reserve(train_groups.size());
      best_rows_of(train_groups, bg_row, bg_score);
      for (std::size_t i = 0; i < train_groups.size(); ++i)
      {
        (group_label[train_groups[i]] == 1 ? seed_pos : seed_neg).push_back(bg_row[i]);
      }
      // The seed must use the SAME learner as the loop. Seeding a GBT run with an LDA fit
      // reintroduces exactly the cold start this block exists to prevent, one level up: on data
      // whose signal is an interaction, the LDA seed is at chance by construction, so iteration 0
      // selects no positives, the GBT never fits, and the run silently falls back to ranking by a
      // hyperplane that cannot see the signal. Caught by odia_gbt_test T7, which reported 0 IDs
      // for BOTH classifiers before this.
      //
      // params.seed_mask applies HERE and only here (mechanism 1) -- the seed picks the anchors,
      // so it is the seed that must not see the calibrated features.
      if (!fit_learner(seed_pos, seed_neg, params.seed_mask))
      {
        // A failed seed is the cold start this block exists to prevent. Recording it means the
        // caller's "fitted 0 iterations" warning can fire instead of the run silently ranking by
        // the one-feature bootstrap while claiming a learned model.
        ++n_skipped;
      }
    }

    std::vector<std::size_t> prev_positives;   // mechanism 5's state, sorted
    bool prev_is_main_threshold = false;       // provenance of prev_positives: selected at
                                               // train_fdr (iteration >= 1), never at the 0.15
                                               // initial cut. The 2026-09-01 smoke showed a
                                               // loop-counter guard is not provenance: any skip
                                               // desynchronises them and the cross-threshold
                                               // comparison comes back.
    bool fit_ok_last_iter = false;             // whether the PREVIOUS iteration refit the model.
                                               // A failed or skipped fit leaves the model -- and
                                               // therefore the next selection -- unchanged, and an
                                               // unchanged selection reads as Jaccard 1.0: a fold
                                               // whose fits keep failing would report "converged".
    AnchorTrainingReport stop_rep;             // PERSISTS across iterations: the high-water mark
                                               // and patience strikes are cumulative state, and a
                                               // per-iteration report silently reduces patience
                                               // to single-shot.
    std::vector<std::size_t> log_prev;         // iteration_log's own state, sorted -- deliberately
                                               // NOT shared with mechanism 5, so the diagnostic
                                               // can never perturb the stop decision
    for (int iteration = 0; iteration < std::max(0, params.n_iter); ++iteration)
    {
      // Reduce training scores to the best candidate row per precursor, then
      // estimate group-level q-values for confident positive selection.
      std::vector<lda_detail::RankedGroup> ranked;
      ranked.reserve(train_groups.size());
      best_rows_of(train_groups, bg_row, bg_score);
      for (std::size_t i = 0; i < train_groups.size(); ++i)
      {
        const std::size_t g = train_groups[i];
        ranked.push_back({g, bg_row[i], group_label[g], bg_score[i], 1.0});
      }
      lda_detail::assignQValues(ranked, params.use_pi0);

      // First pass is seeded by a SINGLE feature, so a strict cut can select too few positives to
      // fit anything and the model never bootstraps -- the failure this split exists to prevent.
      const double raw_fdr = (iteration == 0) ? params.train_fdr_initial : params.train_fdr;
      const double train_fdr = std::max(0.0, std::min(1.0, raw_fdr));
      std::vector<std::size_t> positive_rows;
      for (const auto& candidate : ranked)
      {
        if (candidate.label == 1 && candidate.qvalue <= train_fdr)
        {
          positive_rows.push_back(candidate.best_row);
        }
      }
      // The churn diagnostic reads the selection BEFORE the too-few-positives skip below, so a
      // skipped iteration still shows its (tiny) positive set instead of vanishing from the log.
      // The iteration 0 -> 1 overlap spans the train_fdr_initial -> train_fdr threshold change and
      // is expected to be low on a HEALTHY run; it is printed rather than suppressed -- it is the
      // first comparison anyone asks about -- and marked so nobody reads it as a collapse.
      if (params.iteration_log)
      {
        std::vector<std::size_t> curr = positive_rows;
        std::sort(curr.begin(), curr.end());
        if (log_prev.empty())
        {
          std::fprintf(stderr, "[lda] fold %d iter %d: positives %zu jaccard --\n",
                       fold, iteration, curr.size());
        }
        else
        {
          std::fprintf(stderr, "[lda] fold %d iter %d: positives %zu jaccard %.4f%s\n",
                       fold, iteration, curr.size(), jaccardOverlap(log_prev, curr),
                       iteration == 1 ? " (crosses the initial->main train_fdr change)" : "");
        }
        log_prev.swap(curr);
      }
      // Mechanism 5, checked BEFORE this iteration's skips and fit. Three placement consequences,
      // each fixing a reviewed defect of the post-fit version (analysis77, 2026-09-01/02 reviews):
      //  * a collapse now breaks BEFORE a model is trained on the collapsed selection, so the
      //    previous model is what survives;
      //  * a catastrophic shrink below m+2 is seen by the comparison instead of being skipped
      //    around (the too-few-positives `continue` used to bypass the stop entirely, i.e. the
      //    rule went blind exactly at the largest collapse);
      //  * comparisons are gated on PROVENANCE (prev selected at train_fdr) and on the previous
      //    iteration having actually refit -- a failed/skipped fit repeats the same selection,
      //    and an unchanged selection is Jaccard 1.0, which would read as "converged".
      if (params.stop_on_composition)
      {
        std::vector<std::size_t> curr = positive_rows;
        std::sort(curr.begin(), curr.end());
        bool go = true;
        if (prev_is_main_threshold && iteration >= 1 && fit_ok_last_iter)
        {
          AnchorTrainingParams ap;
          ap.stop_jaccard = params.stop_jaccard;
          // n_iter + 1, NOT n_iter: the for-loop's own bound is the cap here, and the helper's
          // pre-fit cap check double-counted it -- an armed k12 run returned "cap" at iteration
          // 11 BEFORE fit 11 ran, silently delivering a k11 model under a k12 label (codex S1).
          // With the helper's cap unreachable, cap termination is the loop's natural exit, and
          // the verdict for it prints after the loop.
          ap.max_iterations = params.n_iter + 1;
          ap.shrink_floor = params.stop_shrink_floor;
          ap.stop_patience = params.stop_patience;
          stop_rep.iterations_run = iteration;
          go = anchorIterationShouldContinue(prev_positives, curr, ap, stop_rep);
        }
        else if (stop_rep.collapse_strikes >= 1 && positive_rows.size() < m + 2)
        {
          // A catastrophic collapse starves the patience of its second strike: the tiny set
          // skips the fit, the fit gate then blocks every later comparison (the frozen model
          // reproduces the same selection forever), and the loop would burn to the cap with no
          // verdict -- silence exactly at the tail event this mechanism exists to catch
          // (kimi F1 / codex S1, found before any armed run shipped). A starved selection
          // arriving with a floor breach already on record IS the second strike.
          stop_rep.collapsed = true;
          stop_rep.note = "positive set starved below m+2 (" + std::to_string(positive_rows.size()) +
                          ") with a floor breach already on record";
          go = false;
        }
        prev_positives.swap(curr);
        prev_is_main_threshold = (iteration >= 1);
        if (!go)
        {
          // The stop VERDICT is part of the run's output contract: without this line the b12-style
          // arm cannot say whether it converged or collapsed and the readout has to be
          // reverse-engineered from the churn log.
          std::fprintf(stderr, "[lda] fold %d iter %d: STOP %s -- %s\n", fold, iteration,
                       stop_rep.converged ? "converged" : "collapsed", stop_rep.note.c_str());
          break;
        }
      }
      // From here to the fit, every early exit means "the model did not change this iteration" --
      // recorded so the next iteration's stop comparison knows its selection is a repeat.
      fit_ok_last_iter = false;
      // Too few confident positives to fit an m-dimensional discriminant. Skipping is right, but it
      // used to be SILENT -- and silence here is dangerous: if every iteration skips, `w` stays at
      // its initialisation (a single feature, weight +/-1), so the "LDA" degenerates to ranking by
      // one sub-score and nothing says so. Record it; the caller reports it.
      if (positive_rows.size() < m + 2)
      {
        ++n_skipped;
        continue;
      }
      // n_trained is incremented AFTER the fit succeeds, further down -- not here. Counting it at
      // selection time reported a trained iteration for a fit that then failed (an empty bag from
      // bag_fraction <= 0, a singular covariance), and the scores in that case come from the
      // previous model or the one-feature fallback.

      // The positive class above is one row per precursor -- `candidate.best_row`, the highest
      // scoring peak group. Taking ALL decoy rows as negatives would therefore make the two classes
      // differ in two ways at once: target vs decoy (wanted) AND rank-1 vs runner-up (not wanted).
      // With ~4-6 candidate peak groups per precursor the negative class is then dominated by
      // runner-ups, so the fitted direction partly separates "best peak in its group" from "not the
      // best peak" -- a real, learnable axis that carries no target/decoy information. Both the
      // negative mean and the within-class covariance are pulled by it. Matching the positive side
      // (top peak per precursor) removes the confound; this is what pyprophet does
      // (semi_supervised.py: td_peaks = train.get_top_decoy_peaks()).
      std::vector<std::size_t> negative_rows;
      if (!params.top_decoys_only)
      {
        for (const std::size_t g : train_groups)
        {
          if (group_label[g] == 1) { continue; }
          negative_rows.insert(negative_rows.end(), group_rows[g].begin(), group_rows[g].end());
        }
      }
      else
      {
        // The ranking scan above already computed every training group's best row with THIS SAME
        // model -- nothing refits in between -- so re-scanning the decoys was pure duplicated work
        // on the most expensive operation in the loop. Reuse it.
        for (std::size_t i = 0; i < train_groups.size(); ++i)
        {
          if (group_label[train_groups[i]] == 1) { continue; }
          negative_rows.push_back(bg_row[i]);
        }
      }
      if (negative_rows.size() < 2) { continue; }

      // A failed fit leaves the previous model in place, exactly as a failed Cholesky leaves the
      // previous w -- the iteration is skipped, not replaced with something degenerate. And it is
      // COUNTED as skipped, so a run whose every iteration failed cannot report itself as trained.
      if (fit_learner(positive_rows, negative_rows, {}))
      { ++n_trained; fit_ok_last_iter = true; }
      else { ++n_skipped; }
    }
    // Every armed termination gets a verdict, not only the explicit breaks: cap exhaustion,
    // k <= 2 (no same-threshold pair ever forms), and a final comparison gated out by a failed
    // fit all used to end in silence, and "no verdict line" is indistinguishable from a broken
    // log (kimi F4 / codex "every armed termination: not printed").
    if (params.stop_on_composition && !stop_rep.converged && !stop_rep.collapsed)
    {
      std::fprintf(stderr,
                   "[lda] fold %d: STOP cap -- ran all %d iterations without a convergence or "
                   "collapse verdict (final positives %zu, last flow +%zu/-%zu)\n",
                   fold, std::max(0, params.n_iter), prev_positives.size(),
                   stop_rep.last_entries, stop_rep.last_exits);
    }

    // This model has seen no row from the groups scored.
    //
    // Parallel over groups: each iteration writes result.dscore at indices belonging to ITS OWN
    // group and reads only the (now fixed) model, so there is no reduction, no shared accumulator,
    // and no ordering question -- the output is bit-identical to the serial loop at any thread
    // count. Worth doing because for the NN this is one forward pass per row through 12 nets, and
    // it was the serial tail of an otherwise parallel routine.
    std::vector<std::size_t> score_groups;
    score_groups.reserve(group_count / static_cast<std::size_t>(folds) + 1);
    for (std::size_t g = 0; g < group_count; ++g)
    {
      if (group_fold[g] == fold) { score_groups.push_back(g); }
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(inner_threads) if (score_groups.size() > 1024)
#endif
    for (long long k = 0; k < static_cast<long long>(score_groups.size()); ++k)
    {
      for (const std::size_t row : group_rows[score_groups[static_cast<std::size_t>(k)]])
      {
        result.dscore[row] = score_row(row);
      }
    }

    // Every fold has its OWN weight vector, and an LDA direction is defined only up to scale (and
    // the mean-offset depends on that fold's training set), so fold A's score of 4.0 and fold B's
    // of 4.0 mean different things. The final q-values below pool all folds into ONE ranking, which
    // silently assumes they are commensurable. They are not, and the pooled ranking is then partly
    // sorted by which fold a precursor happened to land in. Rescaling each fold's held-out scores to
    // that fold's own DECOY null (mean 0, sd 1) puts every fold on the one scale target-decoy FDR
    // actually cares about -- "how many decoy sds above the null" -- and makes pooling valid.
    // pyprophet sidesteps the problem differently: it averages the fold weight vectors into a single
    // model and rescores everything with it (at the cost of the leakage-free property kept here).
    // Keep this fold's trained ensemble for `model_out`. Fold-indexed slots are
    // disjoint across the parallel loop, so no synchronisation is needed.
    if (use_gbt && gbt.trained() && static_cast<std::size_t>(fold) < fold_models.size())
    { fold_models[static_cast<std::size_t>(fold)] = gbt; }
    if (params.normalize_folds)
    {
      double sum = 0.0, sum_sq = 0.0;
      std::size_t decoy_n = 0;
      for (std::size_t g = 0; g < group_count; ++g)
      {
        if (group_fold[g] != fold || group_label[g] == 1) { continue; }
        std::size_t best_row = group_rows[g].front();
        for (const std::size_t row : group_rows[g])
        {
          if (result.dscore[row] > result.dscore[best_row]) { best_row = row; }
        }
        sum += result.dscore[best_row];
        sum_sq += result.dscore[best_row] * result.dscore[best_row];
        ++decoy_n;
      }
      double mu = 0.0, sigma = 0.0;
      bool affine_ok = false;
      if (decoy_n >= 2)
      {
        mu = sum / static_cast<double>(decoy_n);
        const double var = std::max(0.0, (sum_sq - sum * mu) / static_cast<double>(decoy_n - 1));
        sigma = std::sqrt(var);
        if (static_cast<std::size_t>(fold) < fold_norm_mu.size())
        {
          // Saved with the fold model so a frozen application normalises
          // exactly as this fold did (a degenerate sigma is saved as-is and
          // the loader applies the same validity rule below).
          fold_norm_mu[static_cast<std::size_t>(fold)] = mu;
          fold_norm_sigma[static_cast<std::size_t>(fold)] = sigma;
        }
        // A degenerate null (all decoys identical) carries no scale information; leaving those
        // scores unscaled is the only honest option, and shifting them alone would be worse.
        affine_ok = (sigma > std::numeric_limits<double>::epsilon() && std::isfinite(sigma));
      }
      if (params.fold_pool_rank)
      {
        // v1.13: within-fold rank pooling. Knots = every group's best in THIS fold (targets over all
        // rows, decoys over their drawn prefix -- exactly what assignQValues ranks), sorted descending.
        // A group's best row lands on its knot: frac = (number of knots >= v) / n_f, so tied groups
        // share the MAX rank (the conservative side; verified on F17: the equality branch below is
        // unreachable because upper_bound(greater<>) stops at the first knot < v, and the
        // interpolation branch then yields exactly k/n_f); other rows interpolate between
        // adjacent knots so the within-group order is kept.
        std::vector<double> knots;
        for (std::size_t g = 0; g < group_count; ++g)
        {
          if (group_fold[g] != fold) { continue; }
          const std::size_t take = (group_label[g] == 1) ? group_rows[g].size()
                                                         : std::min(draw_from[g], group_rows[g].size());
          double best = -std::numeric_limits<double>::infinity();
          for (std::size_t k = 0; k < take; ++k) { best = std::max(best, result.dscore[group_rows[g][k]]); }
          if (std::isfinite(best)) { knots.push_back(best); }
        }
        std::sort(knots.begin(), knots.end(), std::greater<double>());
        const std::size_t n_f = knots.size();
        if (n_f >= 2)
        {
          const double nf = static_cast<double>(n_f);
          for (std::size_t g = 0; g < group_count; ++g)
          {
            if (group_fold[g] != fold) { continue; }
            for (const std::size_t row : group_rows[g])
            {
              const double v = result.dscore[row];
              const std::size_t k = static_cast<std::size_t>(
                std::upper_bound(knots.begin(), knots.end(), v, std::greater<double>()) - knots.begin());
              double frac;
              if (k < n_f && v == knots[k]) { frac = static_cast<double>(k + 1) / nf; }
              else if (k == 0)
              { frac = (1.0 - (v - knots[0]) / (knots[0] - knots[1] + 1e-12)) / nf; }
              else if (k < n_f)
              { frac = (static_cast<double>(k) + (knots[k - 1] - v) / (knots[k - 1] - knots[k] + 1e-300)) / nf; }
              else
              { frac = 1.0 + (knots[n_f - 1] - v) / (1.0 + knots[n_f - 2] - knots[n_f - 1]); }
              frac = std::max(frac, 1e-12);
              result.dscore[row] = -std::log10(frac);
            }
          }
        }
      }
      else if (affine_ok)
      {
        for (std::size_t g = 0; g < group_count; ++g)
        {
          if (group_fold[g] != fold) { continue; }
          for (const std::size_t row : group_rows[g])
          {
            result.dscore[row] = (result.dscore[row] - mu) / sigma;
          }
        }
      }
    }
  }

  }   // if (!frozen_applied): the training path

  // ---- SAVE, gbt engine only ----
  // The fold ensemble plus everything a frozen application needs to score a
  // NEW run identically to how this one scored its own held-out folds: this
  // run's standardisation (mean, sd per feature) and each fold's decoy
  // normalisation (mu, sigma). Written after every fold has finished.
  if (!frozen_applied && !model_out.empty())
  {
    if (!gbt_engine)
    {
      std::fprintf(stderr, "[lda] -classifier_model_out is honoured by the gbt and percolator "
                           "engines only; nothing written\n");
    }
    else
    {
      // Written to a sibling temporary and renamed into place, so a failed or
      // interrupted save never leaves a partial file that a later
      // -classifier_model_in would refuse (or, worse, half-read).
      const std::string tmp = model_out + ".partial";
      std::size_t k_have = 0;
      for (const auto& g : fold_models) { k_have += g.trained() ? 1 : 0; }
      bool written = false;
      if (k_have > 0)
      {
        std::ofstream os(tmp);
        os.imbue(std::locale::classic());
        os.precision(17);
        os << "ODIA-GBT-FOLDS 1 " << m << ' ' << fold_models.size() << '\n';
        // The sub-score names travel with the model; the loader refuses a
        // run whose sub-scores differ, whatever their count.
        os << m;
        for (std::size_t j = 0; j < m; ++j)
        {
          os << ' ' << ((feature_names != nullptr && j < feature_names->size())
                          ? (*feature_names)[j] : std::string("var_") + std::to_string(j));
        }
        os << '\n';
        for (std::size_t j = 0; j < m; ++j) { os << (j ? " " : "") << mean[j]; }
        os << '\n';
        for (std::size_t j = 0; j < m; ++j) { os << (j ? " " : "") << sd[j]; }
        os << '\n';
        for (std::size_t f = 0; f < fold_models.size(); ++f)
        {
          const bool has = fold_models[f].trained();
          os << (has ? 1 : 0) << ' ' << fold_norm_mu[f] << ' ' << fold_norm_sigma[f] << '\n';
          if (has) { fold_models[f].save(os); }
        }
        // The precursor-to-fold assignment. A frozen application scores a
        // precursor this ensemble TRAINED ON with the one fold that excluded
        // it -- its native held-out score -- and only a precursor the ensemble
        // never saw with the fold mean. Without this, re-scoring the saving
        // run hands every row K-1 models that memorised it, and a control
        // scored that way is optimistic against any arm it is compared to.
        os << "FOLDMAP " << group_count << '\n';
        for (std::size_t g = 0; g < group_count; ++g)
        { os << group_id[g] << ' ' << group_fold[g] << '\n'; }
        os.flush();
        written = static_cast<bool>(os);
        os.close();
        written = written && !os.fail();
        if (written) { written = (std::rename(tmp.c_str(), model_out.c_str()) == 0); }
        if (!written) { std::remove(tmp.c_str()); }
      }
      if (written)
      {
        std::fprintf(stderr, "[gbt] TRAINED -> %s: %zu of %zu fold models saved, %zu features\n",
                     model_out.c_str(), k_have, fold_models.size(), m);
      }
      else
      {
        std::fprintf(stderr, "[gbt] TRAINED -> %s: NOT written (%zu trained folds%s)\n",
                     model_out.c_str(), k_have, k_have > 0 ? ", write or rename failed" : "");
      }
    }
  }

  result.n_iterations_trained = n_trained;
  result.n_iterations_skipped = n_skipped;

  std::vector<lda_detail::RankedGroup> final_ranked;
  final_ranked.reserve(group_count);
  for (std::size_t g = 0; g < group_count; ++g)
  {
    std::size_t best_row = group_rows[g].front();
    const std::size_t take = std::min(draw_from[g], group_rows[g].size());
    for (std::size_t k = 0; k < take; ++k)
    {
      const std::size_t row = group_rows[g][k];
      if (result.dscore[row] > result.dscore[best_row]) { best_row = row; }
    }
    final_ranked.push_back(
      {g, best_row, group_label[g], result.dscore[best_row], 1.0});
  }
  lda_detail::assignQValues(final_ranked, params.use_pi0);
  for (const auto& ranked_group : final_ranked)
  {
    for (const std::size_t row : group_rows[ranked_group.group_index])
    {
      result.qvalue[row] = ranked_group.qvalue;
      result.pvalue[row] = ranked_group.pvalue;
      result.pep[row]    = ranked_group.pep;
    }
  }
  return result;
}

} // namespace ODIA::Scoring

#endif // ODIA_LDA_H
