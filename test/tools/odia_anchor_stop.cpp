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
  // far above any tested n (a base of 10'000 collided with plain ids once n
  // exceeded it, silently turning a set into a multiset) -- same cardinality,
  // controllable Jaccard.
  std::vector<std::size_t> members(std::size_t n, std::size_t swap_out = 0, std::size_t salt = 0)
  {
    std::vector<std::size_t> v(n);
    std::iota(v.begin(), v.end(), 0);
    for (std::size_t i = 0; i < swap_out && i < n; ++i) { v[i] = 1000000 + salt * 10000 + i; }
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
  // SUSTAINED high-water mark (second-highest observed -- 2500 is seen once
  // and never sustained, so the mark settles at 2400) on two consecutive
  // comparisons.
  {
    AnchorTrainingParams p;
    p.shrink_floor = 0.01;
    p.stop_patience = 2;
    p.max_iterations = 12;
    AnchorTrainingReport rep;
    rep.iterations_run = 2;
    bool go = anchorIterationShouldContinue(members(2500), members(2400, 200, 1), p, rep);
    expect(go && rep.collapse_strikes == 0 && rep.high_water == 2400,
           "T4a repaired: one-off start is not the mark; no strike vs the sustained level");
    rep.iterations_run = 3;
    go = anchorIterationShouldContinue(members(2400, 200, 1), members(2350, 200, 2), p, rep);
    expect(go && rep.collapse_strikes == 1, "T4b repaired: breach of the sustained mark = strike");
    rep.iterations_run = 4;
    go = anchorIterationShouldContinue(members(2350, 200, 2), members(2300, 200, 3), p, rep);
    expect(!go && rep.collapsed, "T4c repaired: second consecutive breach stops");
    expect(rep.high_water == 2400, "T4d the sustained mark held while the top (2500) did not");
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
    (void)anchorIterationShouldContinue(members(2400), members(2400, 200, 1), p, rep);  // seed mark
    rep.iterations_run = 3;
    (void)anchorIterationShouldContinue(members(2400, 200, 1), members(2350, 200, 2), p, rep);  // strike
    rep.iterations_run = 4;
    const bool go = anchorIterationShouldContinue(members(2350, 200, 2), members(2390, 100, 3), p, rep);
    expect(go && rep.collapse_strikes == 0, "T5 repaired: recovery above floor resets strikes");
  }

  // T8. Spike amnesty (kimi F2): a single-iteration admission spike must not
  // pin the mark and amputate the healthy settle-back to equilibrium.
  {
    AnchorTrainingParams p;
    p.shrink_floor = 0.01;
    p.stop_patience = 2;
    p.max_iterations = 12;
    AnchorTrainingReport rep;
    rep.iterations_run = 2;
    (void)anchorIterationShouldContinue(members(10000), members(10200, 900, 1), p, rep);  // spike
    rep.iterations_run = 3;
    bool go = anchorIterationShouldContinue(members(10200, 900, 1), members(10050, 900, 2), p, rep);
    expect(go && rep.collapse_strikes == 0,
           "T8a spike amnesty: settle-back below a one-off top is no breach");
    rep.iterations_run = 4;
    go = anchorIterationShouldContinue(members(10050, 900, 2), members(10040, 900, 3), p, rep);
    expect(go && !rep.collapsed, "T8b spike amnesty: stable plateau above the sustained mark runs on");
  }

  // T9. A floor breach resets convergence strikes (codex S1): a smoothly
  // contracting set has high adjacent Jaccard by construction, and without
  // the reset a genuine collapse banks convergence strikes on the way down.
  {
    AnchorTrainingParams p;
    p.shrink_floor = 0.01;
    p.stop_patience = 2;
    p.stop_jaccard = 0.98;
    p.max_iterations = 12;
    AnchorTrainingReport rep;
    rep.iterations_run = 2;
    // 10000 -> 9950, 5 swapped: J = 9945/10005 = 0.994 >= 0.98 -> converge strike 1; mark 9950.
    bool go = anchorIterationShouldContinue(members(10000), members(9950, 5, 1), p, rep);
    expect(go && rep.converge_strikes == 1, "T9a smooth contraction banks a convergence strike");
    rep.iterations_run = 3;
    // 9950 -> 9840: below 9950 * 0.99 = 9850.5 -> collapse strike AND converge reset (J = 0.983).
    go = anchorIterationShouldContinue(members(9950, 5, 1), members(9840, 5, 1), p, rep);
    expect(go && rep.collapse_strikes == 1 && rep.converge_strikes == 0,
           "T9b floor breach resets convergence strikes; no converged verdict on the way down");
  }

  // T10. patience 0 from a library caller is clamped, not "always collapse"
  // (kimi F3: collapse_strikes >= 0 is true even after a reset).
  {
    AnchorTrainingParams p;
    p.shrink_floor = 0.01;
    p.stop_patience = 0;
    p.max_iterations = 12;
    AnchorTrainingReport rep;
    rep.iterations_run = 2;
    const bool go = anchorIterationShouldContinue(members(2400), members(2400, 30, 1), p, rep);
    expect(go && !rep.collapsed, "T10 patience 0 is clamped to 1, not an unconditional stop");
  }

  // T11. Flow accounting: entries/exits of the comparison are recorded for
  // the log (the observable that separates sharpening from seed-nesting).
  {
    AnchorTrainingParams p;
    p.shrink_floor = 0.01;
    p.stop_patience = 2;
    p.max_iterations = 12;
    AnchorTrainingReport rep;
    rep.iterations_run = 2;
    (void)anchorIterationShouldContinue(members(2500), members(2400, 200, 1), p, rep);
    // curr: 200 fresh ids (10100+) not in prev -> entries 200; prev loses 0..199 and keeps
    // nothing above 2399 -> exits 100 + 200 = 300.
    expect(rep.last_entries == 200 && rep.last_exits == 300,
           "T11 entry/exit flow recorded exactly");
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
