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

* **Never judge a feature-sized change on one fixture arm again.** The band is
  at least +-3% at the operating point and wider at 15%.
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
