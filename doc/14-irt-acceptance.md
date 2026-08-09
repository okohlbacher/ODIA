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

## 2026-08-09 16:20: the SD <= 30 s bar is NOT reachable by calibration

`test/rt_calibration_residuals.py` already existed and answers the question the
whole thread was circling. It starts from the LIBRARY's iRT -- the only honest
source for "what did we predict before seeing the run" -- fits a calibration on
one set of stripped sequences and evaluates on another, so the number is not the
calibration's own flexibility.

Astral, 10,891 precursors at 1% FDR, 8,723 train / 2,168 held out:

    calibration              set          n        SD        p95    (minutes)
    linear                   in-sample  8,723    0.702      1.398
    linear                   held out   2,168    0.704      1.379
    monotone (LOESS-like)    in-sample  8,723    0.625      1.233
    monotone (LOESS-like)    held out   2,168    0.644      1.269

**0.644 min = 38.6 s is the CEILING for any monotone post-hoc calibration on
this library and run.** ODIA's map is a monotone post-hoc calibration. It cannot
beat 38.6 s held-out however it is fitted, so the SD <= 30 s bar is unreachable
by anything in the calibration path -- LOESS spans, anchor thresholds, bin
counts, isotonic tweaks. Those can close the gap from our current 95.6 s (in
sample, contaminated anchors) down towards 38.6 s and no further.

**And it corrects the DIA-NN comparison.** DIA-NN's reported `RT -
Predicted.RT` SD is 29.29 s, which is BELOW the ceiling a monotone map of the
library iRT can achieve. That is only possible because `Predicted.RT` is
DIA-NN's RUN-REFINED prediction, not a monotone function of the library value.
So the earlier statement "our map is 34.22 against DIA-NN's 26.09" compared a
calibration against a retrained model. The right reading is that DIA-NN is doing
something our architecture does not do at all.

**Therefore fine-tuning is not an optimisation, it is the only path to the bar.**
doc/06 measured exactly this: held-out 0.701 -> 0.449 min (42.1 -> 26.9 s) from
500 peptides in 31 s, 0.395 min (23.7 s) from 2,000. Those numbers straddle the
30 s bar where no calibration can approach it.

Caveat to check before quoting the ceiling as final: `astral_lib_own.tsv` was
derived from a DIA-NN run on this same file, so its iRT column may already carry
run information, which would make 38.6 s optimistic rather than pessimistic.

S08 cannot be measured this way yet: the v6 report yields 738 usable rows
against the script's 1,000-row floor.

## 2026-08-09 16:45: sequence information DOES beat the ceiling (measured locally)

`test/rt_sequence_gain.py`. Astral, 10,891 paired precursors, 8,723 train /
2,168 held out, split by stripped sequence:

    monotone map of library iRT only            held-out SD  38.55 s
    ridge(AA composition + calibrated iRT) l=1  held-out SD  33.19 s
                                          l=10  held-out SD  33.26 s
                                          l=100 held-out SD  35.62 s

The monotone figure reproduces `rt_calibration_residuals.py`'s 38.6 s from an
independent implementation, which is the cross-check that makes the rest
trustworthy.

**A 24-feature ridge with no neural network breaks the ceiling by 5.4 s.** So
the residual left after the best monotone map is not noise -- it carries
sequence-dependent structure, which is precisely the claim doc/06 makes for
retraining and the reason a calibration cannot get there.

Scaling from doc/06's own comparison on S08 (AA-composition ridge 0.616 min
against fine-tuned 0.449, a factor of 0.73), the same factor applied to the
0.553 min ridge here would put a peptdeep fine-tune near 0.40 min = 24 s --
under the 30 s bar. That is an extrapolation across datasets and must be
measured, not quoted.
