// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Does the adopted scoring stack work in this tree?
//
// Not a re-test of the algorithms -- they carry their own tests upstream. This
// checks the things adoption can break and upstream tests cannot see: that the
// namespace rewrite compiles, that the headers link into odia_core, and that
// the scorer separates a planted signal rather than returning the initialisation
// it falls back to when it cannot fit.
//
// That last one is the point. `ScoredGroups::n_iterations_trained == 0` means
// the scores came from a single-feature initialisation, which -- in the
// upstream author's own words -- "looks like a working LDA and is not one". A
// smoke test that only checked for finite q-values would pass on exactly that.
#include <odia/scoring/fdr.h>
#include <odia/scoring/lda.h>

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

int main()
{
  // Two populations with a real but noisy separation, and one pure-noise
  // feature so a scorer that latches onto anything would be caught.
  std::mt19937 rng(20260804);
  std::normal_distribution<double> noise(0.0, 1.0);

  const int n_precursors = 600;
  const int candidates = 3;
  std::vector<std::vector<double>> features;
  std::vector<int> labels;
  std::vector<long long> group;

  for (int p = 0; p < n_precursors; ++p)
  {
    const bool target = (p % 2) == 0;
    for (int c = 0; c < candidates; ++c)
    {
      // Only the first candidate of a target carries the signal: a precursor
      // has several candidate peak groups and at most one is right.
      // Strong on purpose. The point is to check that the stack fits and
      // calibrates, not to probe its sensitivity limit -- a marginal planted
      // effect makes the q <= 0.01 assertion a test of the effect size rather
      // than of the code. At 1.6 sigma the discriminant separated cleanly
      // (targets -0.26 against decoys -1.18) and still nothing reached
      // q <= 0.01, which would have read as a failure of the scorer.
      const double shift = (target && c == 0) ? 4.0 : 0.0;
      features.push_back({shift + noise(rng), 0.7 * shift + noise(rng), noise(rng)});
      labels.push_back(target ? 1 : 0);
      group.push_back(p);
    }
  }

  ODIA::Scoring::LDAParams params;
  const auto scored = ODIA::Scoring::scoreSemiSupervisedLDA(features, labels, group, params);

  if (scored.dscore.size() != features.size())
  {
    std::fprintf(stderr, "scorer returned %zu scores for %zu rows\n",
                 scored.dscore.size(), features.size());
    return 1;
  }
  std::printf("iterations trained %d, skipped %d\n",
              scored.n_iterations_trained, scored.n_iterations_skipped);
  if (scored.n_iterations_trained == 0)
  {
    std::fprintf(stderr, "no iteration ever fitted a discriminant: these scores are "
                         "the single-feature initialisation, not a model\n");
    return 1;
  }

  for (std::size_t i = 0; i < scored.dscore.size(); ++i)
  {
    if (!std::isfinite(scored.dscore[i]) || !std::isfinite(scored.qvalue[i]))
    {
      std::fprintf(stderr, "non-finite score at row %zu\n", i);
      return 1;
    }
  }

  // Targets must outrank decoys on average. A scorer that returned a constant
  // would still give finite q-values everywhere, so this is the check that
  // distinguishes working from merely well-formed.
  double t_mean = 0, d_mean = 0;
  std::size_t nt = 0, nd = 0;
  for (std::size_t i = 0; i < scored.dscore.size(); ++i)
  {
    (labels[i] ? t_mean : d_mean) += scored.dscore[i];
    (labels[i] ? nt : nd) += 1;
  }
  t_mean /= (nt ? nt : 1);
  d_mean /= (nd ? nd : 1);
  std::printf("mean d-score: targets %.4f, decoys %.4f\n", t_mean, d_mean);
  if (!(t_mean > d_mean))
  {
    std::fprintf(stderr, "targets did not outrank decoys; the discriminant is not working\n");
    return 1;
  }

  std::size_t at_one_pct = 0;
  for (std::size_t i = 0; i < scored.qvalue.size(); ++i)
  {
    if (labels[i] && scored.qvalue[i] <= 0.01) { ++at_one_pct; }
  }
  std::printf("target rows at q <= 0.01: %zu\n", at_one_pct);
  if (at_one_pct == 0)
  {
    std::fprintf(stderr, "nothing passed 1%% FDR on a planted separation\n");
    return 1;
  }

  // The entrapment estimator, which the backlog wanted before any q-value is
  // believed and which turns out to be already implemented.
  const auto e = ODIA::Scoring::fdr::entrapmentFdp(1000, 12, 20000, 20000);
  std::printf("entrapment: fdp %.4f, valid %d\n", e.fdp, static_cast<int>(e.valid));
  if (!e.valid || !(e.fdp >= 0.0 && e.fdp <= 1.0))
  {
    std::fprintf(stderr, "entrapment estimate is not a proportion\n");
    return 1;
  }

  std::printf("adopted scoring stack OK\n");
  return 0;
}
