# The fixture cannot resolve a feature-sized change

2026-08-26. Six changes were measured on the s08_6x60 fixture at matched
entrapment FDP and every one came back slightly negative at the operating point
and strongly positive at 15%. A column of PURE NOISE does the same thing.

## The control

`var_null_control` is a deterministic splitmix64 hash of the precursor index and
apex cycle, scaled to [0,1). It varies per candidate, so it survives the
constant-column guard, and it is uniform with respect to label, so it carries
nothing a classifier may legitimately use. It is not a feature; it is an
instrument test.

Two arms, both against a baseline established in the same session:

    arm            IDs   entrap    FDP   sigma   DIA-NN
    terminals    3,419       51   8.75    1.23    3,073
    inertcol     3,419       51   8.75    1.23    3,073
    nullfeat     3,370       49   8.52    1.22    3,038

`inertcol` is HEAD with no flags, and it reproduces `terminals` on every result
column. Two things follow: the all-NaN columns really are dropped by the
constant-column guard, so the baseline used all day was valid; and determinism
holds across a rebuild, not just across a repeat.

`nullfeat` differs from `inertcol` by one column of noise.

## The profile, at matched entrapment FDP

    arm                     7.42%    10%      15%
    Gate C weighting        -2.0%   -1.2%   +8.2%
    PROFILE_FIT             -2.8%   -1.6%   +8.6%
    mass survival 4.5 ppm   -0.9%   -1.2%   +6.6%
    mass survival 2.0 ppm   -3.1%   -0.1%   +8.0%
    PURE NOISE              -1.5%   -3.9%   +5.7%

The noise column sits inside the range of the four real changes on every column.
Its profile is not distinguishable from theirs.

**Adding any column moves the semi-supervised classifier's trajectory further
than a feature's information content does.** This is not non-determinism -- a
repeat run is byte-identical, and `inertcol` shows that holds across rebuilds.
It is deterministic and arbitrary: a different feature set gives a different
solution path, and on a fixture holding 19% of the retention time and ~40
entrapment hits at the operating point, that path difference dominates.

## The band, measured with three seeds rather than inferred from one

One null arm shows the band exists. Its WIDTH needs more than one sample, and
different seeds give different arbitrary trajectories, so:

    arm                  IDs   entrap    FDP    DIA-NN
    inertcol (base)    3,419       51   8.75     3,073
    null seed 0        3,370       49   8.52     3,038
    null seed 7        3,458       54   9.17     3,092
    null seed 99       3,326       57  10.07     2,940

Pure noise spans **132 identifications (3.9%)** and **1.55 pp of FDP** before any
matched-FDP comparison is made. At matched FDP against the same baseline:

    null seed    5.72%    7.42%     10%     15%
       0        -13.7%    -1.5%   -3.9%   +5.7%
       7        -14.4%    -0.8%   -2.3%   +7.9%
      99        -39.6%   -11.9%   -9.9%   +3.1%

**The band at the operating point is -0.8% to -11.9%.** Every change measured
today sits inside it:

    change                  7.42%
    mass survival 4.5 ppm   -0.9%
    Gate C weighting        -2.0%
    PROFILE_FIT             -2.8%
    mass survival 2.0 ppm   -3.1%
    ---- noise band ----    -0.8% to -11.9%

They are not merely indistinguishable from noise; they are at the BETTER end of
it. A verdict of "does not convert" was never available from this instrument.

One more thing the three seeds show, which a single arm could not: all three are
NEGATIVE at the operating point. A column carrying no information is not free --
it is an overfitting surface for a semi-supervised loop with ~40 entrapment hits
to calibrate against. That is consistent with this project's recorded finding
that feature COUNT buys nothing, and sharpens it: count is not neutral, it is
mildly harmful.

## What this invalidates, and what it does not

**Invalidated as evidence about the features:** all six "does not convert"
verdicts. They are real numbers and they were correctly computed; they simply do
not measure what they were read as measuring. That includes the +1.2% earlier
attributed to the boundary work, which sits in the same band.

**NOT invalidated:** the two reverts. Gate C weighting was removed because its
label asymmetry was measured directly -- entrapment admission 81.07% -> 77.62%
against decoy admission 82.24% -> 76.44%, with the mechanism identified (decoys
copy their target's library intensities while their fragment m/z are
recomputed). PROFILE_FIT was removed because its central mechanism was measured
to contribute nothing (no-reference control within noise of the argmax) and
because it was 0.896-correlated with an existing column. Neither rests on an
arm.

**Also not invalidated:** effects far outside this band. The prominence gate's
-20.4%, the 2.6x gap to DIA-NN at its own operating point, and the terminal-
reason accounting that puts Gate C at ~98% of the loss and the picker at ~0.2%
are all an order of magnitude clear of it.

## What to do instead

* **Never judge a feature-sized change on one fixture arm again.** The measured
  band at the operating point is -0.8% to -11.9% over three seeds, so the
  fixture's feature-level resolution is worse than 10%. Nothing in the range any
  single sub-score plausibly delivers is measurable there.
* **Run a null arm alongside any arm whose expected effect is under ~5%**, and
  report the difference against the null rather than against the baseline.
  `-null_feature` with `-null_feature_seed` exists for this; different seeds give
  different arbitrary trajectories, which is how the band's width gets measured
  instead of guessed from the single sample above.
* **Prefer measurements that do not touch the feature vector at all** --
  terminal-reason accounting, admitted-volume counts, gate-log comparisons.
  Those were the ones that produced every durable result today.
* For a real feature verdict, the full run is not optional: ~40 entrapment hits
  at the operating point on the fixture against a much larger count, and no
  19%-of-gradient truncation.

## The uncomfortable part

This control cost one 23-minute arm. It should have been run before the second
feature-level arm of the day, not after the sixth. Every feature-level verdict
issued today was read off an instrument that had never been checked against a
known-null input.


## The cause is model capacity, and there is a dial

The churn is in the DISCRIMINANT, not the threshold. Measured from four arms'
rank files with no new run: Spearman rho between a run and the same run plus a
noise column is 0.79-0.82, mean rank displacement 12-14%, and agreement in the
top 500 is 0.38-0.42 -- worst exactly where it should be best -- recovering to
0.94 only by top-3000.

Two ways of constraining the model, each with a paired null arm at seed 99:

    setting                  rho(base, +noise)   identical
    depth 4 (default)              0.791            no
    min_child_rows 200             0.663            no
    depth 2                        1.000000         YES

**At depth 2 the output is bit-identical with and without the noise column.**
The instability is capacity: 120 trees of depth 4 over ~3,400 positives can find
a different-but-equally-good solution whenever the feature set moves.
Regularising the LEAVES made it worse, which rules out "too little
regularisation" -- min_child_rows 200 coarsens the fit without removing the
freedom to choose among equivalent split sets.

The cost is real. Depth 2 gives 3,385 IDs at 10.07% FDP, and at matched FDP it
is -14.0% at the operating point against one draw of depth 4. Depth 4's own
range across four draws there is 2,778 to 3,152; depth 2's 2,711 sits below all
of them.

So the trade is exactness against depth: **depth 4 cannot be measured, depth 2
can be measured but costs identifications.** Depth 3 is running.

This reframes the whole day. Every feature-level arm was run on a model whose
solution path moves further than the feature does. The right order was always:
fix the instrument, then measure the features.
