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

Depth 3 does not split the difference. It improves rho to 0.893 but leaves the
HEAD of the list exactly as unstable -- top-500 agreement 0.410, against depth
4's 0.418 and depth 2's 1.000. There is no middle: rho rises smoothly while the
quantity that matters does not move at all until depth 2.

    depth   rho        top-500   identical
      4     0.790859   0.418     no
      3     0.892746   0.410     no
      2     1.000000   1.000     YES

So the trade is exactness against depth, with a threshold rather than a slope:
**depth 4 cannot be measured, depth 2 can be measured exactly.**

## The methodology that follows

Depth 4's headline number is one draw from a band. Estimating its central value
from the four draws available at the operating point -- 3,152 / 3,105 / 3,127 /
2,778 -- gives about 3,041, against depth 2's 2,711 and depth 3's 3,074. So
depth 2 costs roughly 11% against depth 4's centre, not the 14% a single draw
suggested.

That is a price worth paying for an instrument, not for production:

* **Screen features at depth 2.** The output is bit-identical to a noise
  perturbation, so a difference of any size is real. This is the first
  configuration in this project where a sub-score can be judged at all.
* **Ship at depth 4**, and confirm a screened winner there with a PAIRED NULL
  arm rather than against the baseline.
* A feature that helps at depth 2 and not at depth 4 is not thereby refuted --
  but it is a different claim, and it needs the paired-null comparison to make.

This reframes the whole day. Every feature-level arm was run on a model whose
solution path moves further than the feature does. The right order was always:
fix the instrument, then measure the features.


## The instrument works, and the first thing it measured said no

`var_mass_survival` at 4.5 ppm, screened at depth 2 against a control differing
only by that column:

    matched FDP    delta
       5.72%       -8.7%
       7.42%       -6.6%   <- operating point
      10.00%       +5.6%
      15.00%       +2.5%

**-6.6% at the operating point, and this one is attributable.** There is no band
to hide in: at depth 2 the output is bit-identical under a noise perturbation,
so the difference is the feature.

The same run validates the instrument, which matters more than the verdict. A
reviewer's objection to screening at depth 2 was that stability might be a
tautology -- a depth-2 model might simply never split on a twentieth column, in
which case "insensitive to noise" would mean "deaf to everything". It does not:

    depth 2 + pure noise        3,385 IDs   IDENTICAL to control
    depth 2 + mass survival     3,295 IDs   -90

It ignored the uninformative column and responded to the informative one. That
is the behaviour an instrument is supposed to have, and it was not designed in
-- it was tested for and could have failed.

One assumption corrected on the way: `-mass_features on` turned out to be a
no-op on this fixture (d2mf is identical to d2base), so `auto` already resolves
to ON here -- presumably the calibration gate fails on a 19% slice, where it
passes on the full file. The reasoning that put `var_mass_survival` on the
`!mass_on` withhold list still holds for the full run; it simply had no effect
on the arms.

## Where that leaves the day

Nothing built today earned its place. Gate C weighting was label-asymmetric,
PROFILE_FIT's mechanism was absent, and mass survival costs 6.6% measured
exactly. Four features, four rejections.

What was built instead is the ability to tell -- and the record of six earlier
verdicts that were never available. That is the more useful half.


## The duplicate-column test was inert, and the reason answers the objection

A reviewer's blocker: depth-2 bit-identity proves only that a column which never
wins a split changes nothing -- a fixed point, not stability -- because
var_null_control sits at index 20 of 21 and split ties break toward the LOWEST
feature index, so it is maximally disadvantaged. The prescribed test was a
column that DOES win splits while carrying nothing new: an epsilon-jittered copy
of var_corr_sum.

Built and run. It changed nothing anywhere:

    arm       IDs   entrap    FDP   vs its base
    d2base  3,385       58  10.07
    d2dup   3,385       58  10.07   IDENTICAL
    fbbase  3,338       50   8.79
    fbdup   3,338       50   8.79   IDENTICAL

Identical even under fixed bins at depth 4, which the HASH null does move
(3,338 / 3,366 / 3,263). So the duplicate is inert in a configuration that is
demonstrably not stable -- which means the test measured its own construction,
not the model. At 1e-6 relative jitter the copy falls in the same 64 quantile
bins as the column it copies, so its gain is exactly equal at every split, and
an exact tie is resolved against it by index. It can never win a split.

**A truly redundant column cannot win a split against its twin.** A tree reads
only order, so any monotone transform of a feature partitions identically; the
only way to make a copy win is to add enough noise that it is no longer
redundant. The reviewer's test, as specified, cannot be constructed.

What answers the objection instead is the mechanism, and the evidence is already
in hand. At depth 2 a tree has 3 internal nodes; at depth 4 it has 15. A uniform
column has near-zero gain, so it legitimately loses every gain comparison
against real features when the budget is 3 -- and wins some when the budget is
15 and deep nodes are fitting small subsets. That is not blindness, it is
selection working, and it is the textbook overfitting surface this codebase's
own comment names at gbt.h:441 ("a leaf of a handful of rows is fitted noise,
and it lands in the score tail that sets the FDR threshold").

The decisive counter-evidence is that depth 2 is NOT inert to a split-winning
column: var_mass_survival moved it by -90 identifications. So depth 2 ignored
the column with no gain and responded to the one with gain, which is the
behaviour required of an instrument.

Residual, and it should be stated rather than argued away: this is one feature
set at one column index. A feature whose gain sits near the depth-2 selection
threshold could be admitted or refused on a margin, and nothing here measures
how wide that margin is.


## What the working instrument found first: the semi-supervised loop is harmful

Every comparison below is EXACT -- depth 2, where base and null arms are
bit-identical -- so these are attributable in a way nothing measured on the
default configuration has been.

    vs d2base, at matched entrapment FDP     5.72%    7.42%     10%
    240 trees, lr 0.05                      -27.1%    +5.2%   +3.6%
    480 trees, lr 0.025                     -14.9%    +7.1%   +0.8%
    n_iter = 1                               +4.0%   +11.2%   +5.1%

**Three iterations of re-selecting the positive set cost 11.2% at the operating
point against a single fit.** The loop is not paying for itself; it is losing
identifications. That is a result the broken instrument could not have produced
-- at depth 4 the same comparison reads -2.3%, comfortably inside the -0.8% to
-11.9% band, i.e. indistinguishable from noise.

More shallow trees also help, +7.1% at 480 rounds, and -- the part that matters
structurally -- **they do not cost exactness**: d2t240 and d2t480 are both
bit-identical to their null arms. Capacity added additively does not recreate
the surface that lets a noise column win a split; capacity added by DEPTH does.
That is the cleanest statement of the mechanism yet.

Both effects point the same way at the operating point and should compose.
Depth 2 sits near 2,711 IDs there; +11.2% and +7.1% together would land near
3,240, against depth 4's 3,152. If that holds, the configuration that can be
MEASURED is also the better one, and the 11% "cost of exactness" was never a
cost -- it was two bad defaults.

Running: the stacked arm with its null, plus the same trees-and-iterations
change at depth 4 to separate "shallow is good" from "these two defaults were
bad at any depth".
