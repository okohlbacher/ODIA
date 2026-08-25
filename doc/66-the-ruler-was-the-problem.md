# The ruler was the problem

2026-08-26. A change that measured well was merged, then measured again on a
different ruler and turned off. This is what happened and why the second ruler
wins.

## What was merged

`corr_sum` -- summed pairwise Pearson of the raw fragment traces over +-4 cycles
-- is the only evidence this detector has ever had at candidate-selection time,
and it is a statement about SHAPE. Every threshold moved in either direction had
already failed (doc/65), which is what "one kind of evidence, exhausted" looks
like.

The library says something shape cannot: which fragments should be BRIGHT.
Pearson over a window is invariant to scale, so two co-eluting species with
equally good shape and contradictory relative intensities are indistinguishable
to `corr_sum` by construction. Ranking the margin survivors by normalised
`corr_sum` + library correlation, cap unchanged:

    statistic                 recall@3   selection accuracy
    corr_sum                    69.5%          59.9%
    library correlation         71.3%          63.6%
    both                        72.7%          63.3%

The first change to this detector to improve both numbers. Merged.

## Why it was turned back off

Both of those columns are agreement with DIA-NN, and this project has already
written down that such a metric is a ceiling and not a result -- recall of
another tool's ID set is positive-unlabelled by construction, and there is a
recorded case where that ceiling rose from 85.1% to 93.8% while realised
identifications **fell by 319 peptides**. The rule, in the vault's words: recall
rising alone is never sufficient evidence.

The prescribed replacement uses no external tool. Target fraction among the
top-N, decoys as the control, ranked by the score being proposed:

    discriminant     corr_sum   corr+lib     (top-500 / top-1000 / top-2000)
    library            96.8%      91.0%      94.7/87.6   86.7/79.3
    co-elution         93.8%      94.8%      92.3/92.8   86.2/86.5
    both               99.2%      99.2%      98.6/98.2   94.2/92.5

Read the three rows together and the pattern is not a defect in the ordering, it
is the **winner's curse**: an arm that selects the position maximising library
correlation and is then ranked by something containing library correlation lets
DECOYS shop for their best value exactly as freely as targets. The more the
discriminant overlaps the selection statistic, the worse the arm looks. On a
discriminant it does not touch, it is mildly ahead.

Which of those three rows the shipped pipeline resembles is the whole question,
and the answer is unfavourable: **`LIBRARY_CORR` is one of the nineteen features
the classifier sees.** Diluted one in nineteen, but on the wrong side.

So `-select_library_weight` ships at 0. Not because the idea is wrong -- the
evidence is genuinely orthogonal to shape, which is exactly what the detector
was short of -- but because what is unproven is that moving it into SELECTION
survives target-decoy competition. The measurement that settles it is a full run
gated on entrapment FDP, and the knob exists so that run needs no patch.

## The same treatment for basin separation

`-candidate_min_separation` was measured in the same session and shipped the same
way, at 1. Suppressing within a chromatographic basin so a cap of three means
three distinct PEAKS rather than three samples of one raises recall@3 from 72.7%
to 75.5% -- confirming that eviction is real, which the non-monotone recall under
a relaxed `apex_evidence` had already implied. It costs selection accuracy,
63.3% to 57.7%, because the freed slots go to genuinely different peaks that can
genuinely outscore the true one.

Three separate recall-increasing changes, three costs on the selection column.
That consistency is itself the finding: with this detector, recall and precision
at the cap trade against each other, and no proxy available offline can say
where the trade lands once nineteen features and target-decoy competition are
involved.

## Two corrections to doc/65 from the same review

**The cap number was measured on a detector that was briefly not the one
running.** While `select_library_weight` defaulted on, production ranked margin
survivors by the blended key while the replay ranked by `corr_sum`. Moot now
that the default is 0, but it was true when written.

**The scan loop never examined three positions per precursor.** The correlation
window is [k-S, k+S+1), which fits for every k from S to n-S-1; the loop ran
`k = S+1` to `n-S-3`. Those positions were not rejected -- they were never
considered, so they appear in no reject counter, which is why no amount of
reading the rejection table would have found them. 198 of 33,337 confident
positives, 0.59%, have their true apex on exactly those three cycles. Fixed.

## What holds

The scoring/quantification split, which now has a regression test asserting both
invariances: moving the quantification bounds leaves every sub-score
bit-identical, and moving the scoring window leaves the reported retention-time
range bit-identical. Each half fails for a different re-coupling, and each
carries its own control so that neither passes vacuously.

## And a number of my own that did not survive re-measurement

The reviewer's other flag was that the scoring/quantification split's AUC gain
might be the apex snap re-centring rather than the narrower width. Checked, and
it is not: both arms of that bootstrap already snapped, so the comparison was
controlled. Separating all three factors over the same paired candidates:

    arm                       median width      AUC
    old rule (no snap)                 128     0.528
    old rule + snap                    128     0.530
    walked bounds, no snap               7     0.784
    walked bounds + snap                 7     0.782
    fixed +-2, no snap                   5     0.790
    fixed +-2 + snap                     5     0.779

    snap alone, holding the old rule       +0.001
    boundaries alone (old -> walked)       +0.253
    width alone (walked -> fixed)          -0.003   (far band only)

The snap is worth nothing. **The boundary fix is worth +0.253 and is the whole
result.**

But re-running the fixed-vs-walked bootstrap at the SHIPPED defaults rather than
the ones in force when it was first run gives **+0.0047, CI [+0.0011, +0.0082]**,
not the +0.0136 quoted in doc/64 and in three comments. The difference is
`peak_min_cycles`: it was briefly 5 when the bootstrap first ran, and moving it
back to 7 -- for reasons that had nothing to do with this measurement -- widens
the walked bounds and lifts the walked arm from 0.751 to 0.760 against the fixed
window's 0.764.

Nobody did anything wrong; a parameter changed two commits after a number was
measured and the number was not re-derived. Corrected in the header, the
implementation comment and the CLI text, all of which quoted the stale figure.

The split stays, but the honest case for it is architectural rather than
numerical: the two intervals answer different questions, DIA-NN separates them
for the same reason, an invariance test can hold them apart, and it costs
nothing. A third of a point of AUC is not why.
