# iRT calibration: what "at least as good as DIA-NN" means

Written 2026-08-09 BEFORE the comparison numbers were in, so the bar cannot move
afterwards.

## The statistic

Post-calibration retention-time residual, in run seconds: the calibrated
prediction minus the observed apex, over the precursors the engine is confident
about. For ODIA that is `reportRtResiduals_`'s "mapped" line, over pass 1's
anchors. For DIA-NN it is `RT - Predicted.RT` at `Q.Value <= 0.01`, converted
from minutes.

Reported as median (does the map remove the OFFSET), robust sigma = 1.4826 x MAD
(the bulk), SD (the tail), and p95|e| (what sets the pass-2 window at
`-rt_window_p95_factor` x p95).

**Robust sigma is the headline** and SD is reported beside it. The residual
distribution is heavy-tailed by construction, and the ratio between the two is
itself the diagnostic: a large gap means anchor contamination rather than a bad
map.

## The targets

**Astral is the primary benchmark, because the comparison is exact.**
`astral_truth.parquet` is DIA-NN on the same raw file, and restricted to the
10,891 precursors of `astral_lib_own.tsv` -- the library ODIA searches -- it
gives:

    n 10,891   median -2.51 s   SD 29.29 s   robust sigma 26.09 s   p95|e| 59.35 s

So the bar on Astral is **robust sigma <= 26.09 s and p95|e| <= 59.35 s**.

**REVISED 2026-08-09 16:00, by the project owner: the bar is SD <= 30 s.**
That is a different statistic from the one above and a deliberately easier one
-- SD is the tail, and ours is inflated by anchor contamination rather than by a
bad map. It is also close to DIA-NN's own Astral SD of 29.29 s, so "at least as
good as DIA-NN" survives the change. Current state against it:

    S08 anchor_q 0.05   SD 71.23     S08 anchor_q 0.01   SD 41.36
    Astral anchor_q 0.05 SD 95.62    DIA-NN Astral       SD 29.29

Both files fail it today, and by a wide margin on Astral.

**S08 is secondary and the comparison is NOT exact.** The available DIA-NN
reference is its v6_50k run: 738 identifications, SD 17.36 s, robust sigma
13.62 s, p95|e| 34.58 s. But only 39 of `lib_targets`' 2,665 precursors appear
in that confident set, so the two populations are different libraries over the
same raw file and the numbers are indicative, not a like-for-like target. S08
must IMPROVE and must be reported; it does not gate on a threshold.

## Two traps this document exists to prevent

**Do not compare across q thresholds.** ODIA's anchors are taken at
`-anchor_q` (default 0.05) and DIA-NN's reference is q <= 0.01. A looser cut
admits worse peaks by construction, so any ODIA-vs-DIA-NN number must state the
anchor q it was measured at, and the headline comparison uses 0.01.

**Do not read the in-sample p95 as the map's quality.** `Calibration::fit`
returns the residual on the very anchors it was fitted to. The 80/20 held-out
probe beside it is what can see overfitting, and both must be quoted.

## When it is frozen

When Astral meets the bar above and S08 has improved, on a build that passes
ctest, with the three reviewers finding no further defect in the wiring. After
that the calibration is not to be touched without a measured reason.
