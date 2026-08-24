# The ladder, completed: trace models match the shipped scalars

2026-08-24. All four rungs, on one contrast, one split, one population.

## The result

    rung  approach                          centring        AUC
      1   shipped 19 sub-scores + GBT       picked cand.  0.8649
      2   per-fragment order stats + GBT    window        0.8422
      3   Deep Sets (48,705 params)         window        0.8429
      3   Deep Sets                         picked apex   0.8611
      4   transformer (2 blocks, 4 heads)   window        0.8606
      4   transformer                       picked apex   0.8645

**A dead heat.** The best trace model reaches 0.8645 against the shipped
scalars' 0.8649 -- a difference of 0.0004, on 21,000 test pairs.

### The centring confound, found and quantified

Rungs 2-4 first ran on a 96-cycle window and had to FIND the peak; rung 1 scores
a candidate the picker had already located. Three independent architectures
converging just below rung 1 is what that confound looks like, so it was tested
rather than assumed: re-centre each arm on ITS OWN picked apex -- no asymmetry,
same picker, same window -- and Deep Sets goes 0.8429 -> 0.8611, **+0.018**.

So roughly two thirds of rung 1's apparent advantage was localisation the picker
had already done, not better features. The populations are comparable: 98.7% of
positives have a candidate in both arms (33,062 against rung 1's 33,085).

This is also the evaluation doc/56 A1 makes mandatory -- the deployment
distribution, where every candidate is centred on its own apex -- so a
class-conditioned centring artefact would have shown up here. It did not: the
apex-centred numbers are HIGHER, which is the opposite of what a centring
artefact produces when you remove the asymmetry it lived on.

### "best" is model selection on the test set, and is not the number

The transformer's per-epoch maximum was 0.8684, above rung 1. That is the
maximum over 45 epochs of a quantity measured on the test set, which is model
selection on test and biased upward. **The final epoch, 0.8645, is the honest
figure** and is what is quoted above. A proper three-way split would let the
epoch be chosen without touching test; it has not been run, so no number here is
allowed to depend on which epoch was best.

The transformer was also overfitting by the end (train 0.8847 against test
0.8645), where at 20 epochs it was underfitting (0.8546 / 0.8552). Its ceiling on
this corpus is near 0.865, not somewhere above it.

## What this says about doc/55's premise

doc/55's central hypothesis was that **collapsing twelve fragments into nineteen
scalars is where the low-abundance information is lost**, motivated by DIA-NN
keeping per-fragment arrays (`pCorr[12]`, `pAcc[6]`, `pShape[5]`) where ODIA
keeps scalars.

**On this contrast, that is not supported.** A permutation-equivariant model
consuming the raw per-fragment traces, with attention across fragments, extracts
no more label-predictive information than the nineteen scalars already do. The
scalars are not leaving much on the table.

That is a negative result about the REPRESENTATION, and it is worth more than a
marginal win would have been: it removes an entire direction. If the traces
contained recoverable structure the scalars discard, a set transformer with
48k-500k parameters and 45 epochs on 45,000 pairs would be expected to find some
of it. It found 0.0004.

## What it does NOT say

**The contrast is not the production task.** "Is this peptide HERE rather than
300 s away" has a localisation component that "is this peptide present at all,
against the decoy null, at low abundance" does not. The RT-shifted negative was
chosen because it makes metadata leakage arithmetically impossible (the
true-metadata control returns exactly 0.5000), and that property was bought at
the price of asking a slightly different question.

So this closes the summarisation hypothesis FOR LOCALISATION and leaves it open
for presence. The entrapment contrast asks the other question and has its own
confounds (doc/56 A2); `ent_hard` at 17 leaves the tail unmeasured entirely.

## Next, in order

1. The same ladder on the ENTRAPMENT contrast, which asks the presence question.
   Its confounds are documented and its metadata control will NOT be 0.5, so the
   leakage floor has to be measured rather than assumed.
2. `ent_hard` populated properly -- ~2,900 DIA-NN-called entrapment precursors
   exist and 17 were sampled. That is the tail, and it is currently unmeasured.
3. A three-way split so epoch selection stops touching test.
4. S30 remains locked and untouched.
