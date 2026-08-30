# 79 — The oracle counterfactual: localization buys +0.3 points

2026-08-31, 00:15. The paired experiment doc/77 section 6 called for, run
overnight per the pre-registered header of `shared/libv2/run_oracle.sh`.

**Design.** samelib3's exact configuration on the current binary (debf9f5a),
two arms differing ONLY in `-oracle_rt`: DIA-NN's 33,330 apex times (both id
naming variants; 32,684 = 98.06% matched onto library indices, the 646
remainder is the bare-cys-id naming gap), `-oracle_rt_tol 20`. The oracle
bypasses admission for listed precursors and synthesises a candidate at the
requested time when none is within tolerance. arm_assert: clean.

**Controls.** ora_ctl reproduced samelib3 EXACTLY: 13,735 at q<=0.01 — the
d7ab850d -> debf9f5a binary evolution is a no-op on this configuration, and
the depth read of ora_ctl matches the samelib3 rank walk to within the join
convention.

**Read** (max DScore per precursor; canonical strip-mods join, with the
alias join as a cross-check — both shown because the earlier chain used the
alias consistently on both sides):

| depth | ora_ctl | ora_on | delta |
|---|---|---|---|
| 13,735 | 12,321 (37.0%) | 12,345 (37.0%) | +24 |
| 33,330 | 21,269 (63.8%) | 21,367 (64.1%) | **+98 (+0.3 pts)** |
| 50,000 | 22,874 (68.6%) | 23,106 (69.3%) | +232 |

**Verdict (pre-registered: <3 pts confirms, >8 amends).** +0.3 points.
Perfect candidate placement at the true apexes — admission gates bypassed,
missing candidates synthesised — moves essentially nothing at the operating
depths. Localization and admission are CONFIRMED minor; the buried class
fails on the evidence content of the features computed at the (correct)
coordinates. doc/77's conclusion now rests on a counterfactual, not only on
the near-RT fork. The evidence-rebuild program (doc/78) is the answer; its
instrument runs were launched the same night.

Note on joins: the canonical (strip-parenthesized-mods) join counts ~490
cys-carrying precursors the alias join silently missed on this file pair
(report ids are bare); all doc/77 numbers used the alias join on BOTH sides,
so every delta and verdict there is unaffected; absolute recalls are ~1.5
points higher under the canonical accounting (62.3 -> 63.8 shipped).
