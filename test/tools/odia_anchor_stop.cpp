// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// anchorIterationShouldContinue: the stopping rule mechanism 5 delegates to.
//
// This helper spent its whole life described as "the tested helper" while no
// test in the tree called it -- and the first time a CLI path reached it
// (2026-09-01) it truncated a 12-iteration arm to 2 iterations in every fold,
// because the legacy rule reads ANY shrink of the positive set as collapse and
// a one-row downtick out of ~2,500 satisfies that. This suite pins BOTH
// policies: the legacy single-shot rule stays bit-compatible for the anchor
// -training call sites that rely on it, and the high-water + patience policy
// behaves as designed on exactly the trajectories that broke the legacy rule.

#include <odia/scoring/anchor_training.h>

#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <vector>

namespace
{
  int failures = 0;

  void expect(bool ok, const char* what)
  {
    if (!ok)
    {
      std::fprintf(stderr, "FAIL: %s\n", what);
      ++failures;
    }
  }

  // A sorted index set 0..n-1, with `swap_out` members replaced by fresh ids
  // above 10'000 -- same cardinality, controllable Jaccard.
  std::vector<std::size_t> members(std::size_t n, std::size_t swap_out = 0, std::size_t salt = 0)
  {
    std::vector<std::size_t> v(n);
    std::iota(v.begin(), v.end(), 0);
    for (std::size_t i = 0; i < swap_out && i < n; ++i) { v[i] = 10000 + salt * 100 + i; }
    std::sort(v.begin(), v.end());
    return v;
  }
}

int main()
{
  using namespace ODIA::Scoring;

  // T1. Legacy defaults: a one-row shrink IS collapse (the behavior the anchor
  // call sites were built against; changing it silently would move their
  // ground under them).
  {
    AnchorTrainingParams p;   // shrink_floor 0, patience 1 = legacy
    p.max_iterations = 12;
    AnchorTrainingReport rep;
    const bool go = anchorIterationShouldContinue(members(2507), members(2506, 30, 1), p, rep);
    expect(!go && rep.collapsed && !rep.converged, "T1 legacy: 1-row shrink stops as collapse");
  }

  // T2. Legacy defaults: one Jaccard >= threshold is convergence, single-shot.
  {
    AnchorTrainingParams p;
    p.max_iterations = 12;
    p.stop_jaccard = 0.98;
    AnchorTrainingReport rep;
    const bool go = anchorIterationShouldContinue(members(2500), members(2500, 10, 1), p, rep);
    // 10 of 2500 swapped: J = 2490/2510 = 0.992 >= 0.98
    expect(!go && rep.converged, "T2 legacy: single high-Jaccard comparison converges");
  }

  // T3. Repaired policy survives the measured mix10k trajectory that killed
  // the legacy rule: one-row net downticks under climbing Jaccard must NOT
  // stop the loop.
  {
    AnchorTrainingParams p;
    p.shrink_floor = 0.01;
    p.stop_patience = 2;
    p.stop_jaccard = 0.98;
    p.max_iterations = 12;
    AnchorTrainingReport rep;
    bool go = true;
    // 2507 -> 2506 -> 2505: two consecutive one-row downticks, high water 2507;
    // 2505 > 2507 * 0.99 = 2481.9, so no strike accrues at all.
    rep.iterations_run = 2;
    go = anchorIterationShouldContinue(members(2507), members(2506, 30, 1), p, rep);
    expect(go && rep.collapse_strikes == 0, "T3a repaired: 1-row shrink accrues no strike");
    rep.iterations_run = 3;
    go = anchorIterationShouldContinue(members(2506, 30, 1), members(2505, 30, 2), p, rep);
    expect(go && !rep.collapsed, "T3b repaired: second 1-row shrink still no collapse");
  }

  // T4. Repaired policy catches a real collapse: cumulative bleed below the
  // high-water floor on two consecutive comparisons.
  {
    AnchorTrainingParams p;
    p.shrink_floor = 0.01;
    p.stop_patience = 2;
    p.max_iterations = 12;
    AnchorTrainingReport rep;
    rep.iterations_run = 2;
    bool go = anchorIterationShouldContinue(members(2500), members(2400, 200, 1), p, rep);
    expect(go && rep.collapse_strikes == 1, "T4a repaired: first breach = strike, not stop");
    rep.iterations_run = 3;
    go = anchorIterationShouldContinue(members(2400, 200, 1), members(2350, 200, 2), p, rep);
    expect(!go && rep.collapsed, "T4b repaired: second consecutive breach stops");
    expect(rep.high_water == 2500, "T4c high-water mark held at the run maximum");
  }

  // T5. A slow bleed that GROWS between breaches resets the strike counter --
  // patience means consecutive, not cumulative-count.
  {
    AnchorTrainingParams p;
    p.shrink_floor = 0.01;
    p.stop_patience = 2;
    p.max_iterations = 12;
    AnchorTrainingReport rep;
    rep.iterations_run = 2;
    (void)anchorIterationShouldContinue(members(2500), members(2400, 200, 1), p, rep);   // strike 1
    rep.iterations_run = 3;
    const bool go = anchorIterationShouldContinue(members(2400, 200, 1), members(2490, 100, 2), p, rep);
    expect(go && rep.collapse_strikes == 0, "T5 repaired: recovery above floor resets strikes");
  }

  // T6. Repaired convergence needs patience consecutive high-Jaccard reads.
  {
    AnchorTrainingParams p;
    p.shrink_floor = 0.01;
    p.stop_patience = 2;
    p.stop_jaccard = 0.98;
    p.max_iterations = 12;
    AnchorTrainingReport rep;
    rep.iterations_run = 2;
    bool go = anchorIterationShouldContinue(members(2500), members(2500, 5, 1), p, rep);
    expect(go && rep.converge_strikes == 1, "T6a repaired: first high-Jaccard = strike, not stop");
    rep.iterations_run = 3;
    go = anchorIterationShouldContinue(members(2500, 5, 1), members(2500, 8, 1), p, rep);
    expect(!go && rep.converged, "T6b repaired: second consecutive high-Jaccard converges");
  }

  // T7. The iteration cap still fires under the repaired policy.
  {
    AnchorTrainingParams p;
    p.shrink_floor = 0.01;
    p.stop_patience = 2;
    p.max_iterations = 3;
    AnchorTrainingReport rep;
    rep.iterations_run = 2;   // helper increments to 3 == max -> stop by cap
    const bool go = anchorIterationShouldContinue(members(2500), members(2500, 200, 1), p, rep);
    expect(!go && !rep.collapsed && !rep.converged, "T7 repaired: cap stops without a verdict");
  }

  if (failures == 0) { std::printf("odia_anchor_stop: all checks passed\n"); return EXIT_SUCCESS; }
  std::fprintf(stderr, "odia_anchor_stop: %d failure(s)\n", failures);
  return EXIT_FAILURE;
}
