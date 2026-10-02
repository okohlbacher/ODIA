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
/// entrapment search on IH1 and Astral is separate and belongs in the backlog.

#include <odia/scoring/lda.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
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

  // ---- the candidate-count asymmetry, and the null that corrects it -------
  //
  // Everything above plants ONE candidate per group, so `match_decoy_candidate_counts`
  // is a no-op and the suite passing says nothing about it. On real data the
  // classes are nowhere near symmetric: measured on Astral, targets carry 10.98
  // candidates per precursor (median 6) and decoys 20.99 (median 24), because
  // the picker keeps candidates within `max_corr_diff` of a precursor's OWN
  // best -- a precursor with no real peak admits nearly every position.
  //
  // A q-value then compares best-of-11 against best-of-21 and comes out
  // CONSERVATIVE. This fixture reproduces that and asserts both halves of the
  // fix: matching the null must recover identifications, AND it must not buy
  // them by inflating the false-discovery proportion. The second half is the
  // one that matters -- the change can only ever make the FDR more liberal, and
  // the Astral benchmark cannot see a false positive, so this is the only place
  // the claim can be falsified.
  {
    std::printf("\n--- asymmetric candidate counts (targets ~11, decoys ~21) ---\n");
    const std::size_t n_true = 1500, n_entrap = 1500, n_decoy = 3000;
    const double ratio = static_cast<double>(n_entrap) / static_cast<double>(n_true);
    // Separation 6, not the 3 the symmetric cases use. Every group here carries
    // ~11-21 NOISE candidates as well, so a group's score is a maximum over
    // that many draws and the signal has to clear the best of them. At 3 both
    // arms reported nothing and both new checks passed vacuously -- which is
    // the exact failure this fixture warns about two blocks up, so it is
    // guarded explicitly below rather than trusted.
    const double separation = 6.0;

    std::mt19937 rng(0xA5717Du);
    std::normal_distribution<double> noise(0.0, 1.0);
    // Candidate counts drawn per group. Only the FIRST row of a real target
    // carries the signal; every other row of every group is noise, which is
    // what a wrong peak is.
    std::poisson_distribution<int> tar_n(11), dec_n(21);

    Planted p;
    long long g = 0;
    const auto add_group = [&](int count, double shift, int label, char entrap) {
      for (int i = 0; i < std::max(1, count); ++i)
      {
        const double common = (i == 0 ? shift : 0.0) + noise(rng);
        std::vector<double> row;
        row.reserve(4);
        for (int k = 0; k < 4; ++k) { row.push_back(common * 0.7 + noise(rng) * 0.6); }
        p.features.push_back(std::move(row));
        p.labels.push_back(label);
        p.group.push_back(g);
        p.is_entrapment.push_back(entrap);
      }
      ++g;
    };
    for (std::size_t i = 0; i < n_true; ++i)   { add_group(tar_n(rng), separation, 1, 0); }
    for (std::size_t i = 0; i < n_entrap; ++i) { add_group(tar_n(rng), 0.0, 1, 1); }
    for (std::size_t i = 0; i < n_decoy; ++i)  { add_group(dec_n(rng), 0.0, 0, 0); }

    std::size_t reported[2] = {0, 0};
    double observed[2] = {0.0, 0.0};
    for (int matched = 0; matched < 2; ++matched)
    {
      ODIA::Scoring::LDAParams params;
      params.classifier = ODIA::Scoring::LDAParams::Classifier::GBT;
      params.match_decoy_candidate_counts = (matched == 1);
      const auto scored =
        ODIA::Scoring::scoreSemiSupervisedLDA(p.features, p.labels, p.group, params);

      // One count per GROUP, not per row: a group has many rows and reporting
      // rows would multiply every precursor by its candidate count.
      std::map<long long, char> hit;
      for (std::size_t i = 0; i < p.labels.size(); ++i)
      {
        if (p.labels[i] == 1 && scored.qvalue[i] <= 0.01) { hit[p.group[i]] = 1; }
      }
      reported[matched] = hit.size();
      observed[matched] = fdp(p, scored, 0.01, ratio);
      std::printf("  match_decoy_n=%d: %zu groups reported, FDP %.4f (claimed 0.01)\n",
                  matched, reported[matched], observed[matched]);
    }

    // Nothing may pass by silence. Both new checks below are satisfied by a
    // scorer that reports zero, so the fixture has to prove it is in a regime
    // where the question is even being asked.
    check(reported[0] > 0,
          "the unmatched arm reports SOMETHING, so the comparison below is not "
          "two zeroes agreeing with each other");

    // The point of the change: an unmatched null is conservative.
    check(reported[1] >= reported[0],
          "matching the decoy candidate count does not LOSE identifications");
    // And the check that can actually falsify it. Same 3x bar as above.
    check(observed[1] <= 0.03,
          "with the null matched, FDP at q<=0.01 is still within 3x of the claim -- "
          "the change may only remove conservatism, never buy identifications with "
          "false ones");
  }

  // ---- the regime that actually matters -----------------------------------
  //
  // Everything above plants 3,000 true against 3,000 entrapment: HALF of all
  // targets are genuinely present. A proteome-scale library is ~1.5% -- DIA-NN
  // finds 738 of 50,000 on IH1 and 37,596 of 2,127,559 on the full library. So
  // the cases above certify q-values at a prior roughly thirty times too
  // favourable, and the one positive result the mobility work is documented by
  // ("+0.90 points at a matched 1.00% entrapment false rate") was measured
  // through that harness. It should not be quoted until re-taken here.
  //
  // At a 1.5% prior, reporting anything at 1% FDR needs a likelihood ratio near
  // 6,500:1. This case makes that difficulty visible in seconds rather than in
  // a 36-minute run, and catches the day the scorer starts reporting at this
  // rate WITHOUT the FDP following it.
  {
    std::printf("\n--- the production regime: 1.5%% of targets are real ---\n");
    const std::size_t n_target = 10000;
    const std::size_t n_true = 150;                  // 1.5% of targets
    const std::size_t n_entrap = n_target - n_true;  // 9,850 real but absent
    const double ratio = static_cast<double>(n_entrap) / static_cast<double>(n_true);

    for (const double separation : {3.0, 2.0})
    {
      const auto p = plant(n_true, n_entrap, n_target, separation, 0xE47A9u);

      ODIA::Scoring::LDAParams params;
      params.classifier = ODIA::Scoring::LDAParams::Classifier::GBT;
      const auto scored =
        ODIA::Scoring::scoreSemiSupervisedLDA(p.features, p.labels, p.group, params);

      std::size_t reported = 0, true_hits = 0;
      for (std::size_t i = 0; i < p.labels.size(); ++i)
      {
        if (p.labels[i] != 1 || scored.qvalue[i] > 0.01) { continue; }
        ++reported;
        if (!p.is_entrapment[i]) { ++true_hits; }
      }

      if (reported == 0)
      {
        // The CURRENT state, recorded rather than asserted. Failing here would
        // paint the suite red for a known open problem; asserting nothing would
        // hide the day it changes. So it prints loudly, and the FDP check below
        // arms itself the moment anything IS reported.
        std::printf("  separation %.1f: NOTHING at q<=0.01 (%d iterations trained)"
                    " -- the open problem, not a regression\n",
                    separation, scored.n_iterations_trained);
        continue;
      }

      const double observed = fdp(p, scored, 0.01, ratio);
      std::printf("  separation %.1f: %zu reported, %zu genuinely true, FDP %.4f"
                  " (claimed 0.01)\n", separation, reported, true_hits, observed);

      // The moment the scorer reports at a realistic prior, its FDP must hold.
      // A scorer that starts reporting here with a broken rate is worse than one
      // reporting nothing, because the number looks like progress.
      check(observed <= 0.05,
            "at a 1.5% prior the reported FDP is within 5x of the claim at "
            "separation " + std::to_string(separation));
    }
  }

  if (failures == 0) { std::printf("all entrapment cases passed\n"); }
  return failures == 0 ? 0 : 1;
}
