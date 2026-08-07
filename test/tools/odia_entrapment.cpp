// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

/// Does a reported 1% FDR mean 1% of the identifications are wrong?
///
/// doc/07-scoring-plan.md lists this as ordering item 4 -- "Entrapment before
/// believing any q-value" -- and it has never been run. Every q-value this
/// project has produced is unvalidated against an independent null, which
/// matters because the target-decoy estimate has been measured failing in
/// three separate ways: it admitted off-RT targets that were statistically
/// indistinguishable from decoys on library_corr (median -0.036 vs -0.032), it
/// collapsed to zero identifications on three occasions while reporting healthy
/// diagnostics, and its output swings ~10% between refits on near-identical
/// input.
///
/// ENTRAPMENT is the check target-decoy cannot perform on itself. Decoys are
/// constructed (shuffled sequences), so the classifier can in principle learn
/// what construction looks like rather than what a wrong answer looks like.
/// Entrapment peptides are REAL peptides that are simply absent from the sample
/// -- they cannot be identified correctly, so every one reported is a genuine
/// false positive, and they carry no construction signature at all.
///
/// The measurement is the false discovery PROPORTION:
///
///     FDP = entrapment hits / (entrapment hits + target hits) * (1/r)
///
/// where r is the entrapment-to-target ratio. If the q-values are calibrated,
/// FDP at q <= 0.01 is about 0.01. If FDP is materially higher, the reported
/// FDR is optimistic by that factor and every recovery number scaled by it.
///
/// This fixture is synthetic and self-contained: it does not need a run, a
/// library on disk, or a search. It builds a score distribution with a KNOWN
/// true/false split, runs it through the same scorer the pipeline uses, and
/// checks the q-values against the truth it planted. That makes it a test of
/// the FDR machinery rather than of any particular dataset -- a real
/// entrapment search on S08 and Astral is separate and belongs in the backlog.

#include <odia/scoring/lda.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

namespace
{
  int failures = 0;

  void check(bool ok, const std::string& what)
  {
    std::printf("  [%s] %s\n", ok ? "ok" : "FAIL", what.c_str());
    if (!ok) { ++failures; }
  }

  /// Three populations, as an entrapment experiment has:
  ///
  ///  * TRUE targets   -- present in the sample, separable, label 1
  ///  * ENTRAPMENT     -- real peptides absent from the sample, label 1 (the
  ///                      scorer cannot know they are wrong), scored from the
  ///                      same distribution as decoys
  ///  * DECOYS         -- label 0
  ///
  /// The scorer sees only the labels. Entrapment rows are targets as far as it
  /// is concerned, so every entrapment row it reports at q <= 0.01 is a false
  /// positive it failed to catch.
  struct Planted
  {
    std::vector<std::vector<double>> features;
    std::vector<int> labels;
    std::vector<long long> group;
    std::vector<char> is_entrapment;
  };

  Planted plant(std::size_t n_true, std::size_t n_entrap, std::size_t n_decoy,
                double separation, unsigned seed)
  {
    // Fixed seed: an FDR test that varies run to run cannot fail reproducibly,
    // and this suite has already been misled once by a green result that came
    // from stale binaries.
    std::mt19937 rng(seed);
    std::normal_distribution<double> noise(0.0, 1.0);
    Planted p;
    long long g = 0;

    auto add = [&](std::size_t n, double shift, int label, char entrap) {
      for (std::size_t i = 0; i < n; ++i)
      {
        // Four correlated features, as sub-scores are: a shared signal plus
        // per-feature noise. A single feature would let the discriminant be
        // exact and would test nothing.
        const double common = shift + noise(rng);
        std::vector<double> row;
        row.reserve(4);
        for (int k = 0; k < 4; ++k) { row.push_back(common * 0.7 + noise(rng) * 0.6); }
        p.features.push_back(std::move(row));
        p.labels.push_back(label);
        p.group.push_back(g++);
        p.is_entrapment.push_back(entrap);
      }
    };

    add(n_true, separation, 1, 0);   // real, separable
    add(n_entrap, 0.0, 1, 1);        // real peptides, absent from the sample
    add(n_decoy, 0.0, 0, 0);         // constructed
    return p;
  }

  /// FDP at a q threshold, corrected for the entrapment-to-target ratio.
  double fdp(const Planted& p, const ODIA::Scoring::ScoredGroups& s,
             double q_threshold, double ratio)
  {
    std::size_t entrap = 0, target = 0;
    for (std::size_t i = 0; i < p.labels.size(); ++i)
    {
      if (p.labels[i] != 1 || s.qvalue[i] > q_threshold) { continue; }
      if (p.is_entrapment[i]) { ++entrap; } else { ++target; }
    }
    const std::size_t total = entrap + target;
    if (total == 0) { return 0.0; }
    return (static_cast<double>(entrap) / static_cast<double>(total)) / ratio;
  }
} // namespace

int main()
{
  // 1:1 entrapment, which makes the correction factor 1 and the arithmetic
  // legible. A real search would use a smaller ratio to avoid diluting the
  // library.
  const std::size_t n_true = 3000, n_entrap = 3000, n_decoy = 3000;
  const double ratio = static_cast<double>(n_entrap) / static_cast<double>(n_true);

  std::printf("entrapment fixture: %zu true, %zu entrapment, %zu decoy (ratio %.1f)\n",
              n_true, n_entrap, n_decoy, ratio);

  bool any_scored = false;
  for (const double separation : {3.0, 2.0, 1.2})
  {
    const auto p = plant(n_true, n_entrap, n_decoy, separation, 0xE47A9u);

    ODIA::Scoring::LDAParams params;
    params.classifier = ODIA::Scoring::LDAParams::Classifier::GBT;
    const auto scored =
      ODIA::Scoring::scoreSemiSupervisedLDA(p.features, p.labels, p.group, params);

    std::size_t reported = 0;
    for (std::size_t i = 0; i < p.labels.size(); ++i)
    {
      if (p.labels[i] == 1 && scored.qvalue[i] <= 0.01) { ++reported; }
    }
    if (reported == 0)
    {
      std::printf("separation %.1f: nothing at q<=0.01 (%d iterations trained)\n",
                  separation, scored.n_iterations_trained);
      continue;
    }
    any_scored = true;

    const double observed = fdp(p, scored, 0.01, ratio);
    std::printf("separation %.1f: %zu reported at q<=0.01, FDP %.4f (claimed 0.01)\n",
                separation, reported, observed);

    // The bar is deliberately loose. A factor of 3 is not "calibrated" by any
    // reasonable standard, but this fixture exists to catch the FDR being
    // wrong by an ORDER of magnitude -- which is the failure mode observed on
    // real data, where an independent positional null put the true rate at
    // 1.88% where the decoys claimed 1%. Tightening it is a follow-up, and
    // requires deciding what calibration this scorer actually promises.
    check(observed <= 0.03,
          "FDP at q<=0.01 is within 3x of the claim at separation " +
            std::to_string(separation));
  }

  check(any_scored, "at least one separation produced identifications at all -- a "
                    "scorer that reports nothing has a trivially perfect FDP and "
                    "must not be allowed to pass by silence");

  if (failures == 0) { std::printf("all entrapment cases passed\n"); }
  return failures == 0 ? 0 : 1;
}
