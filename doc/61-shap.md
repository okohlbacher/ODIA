# What the models actually use: exact TreeSHAP

2026-08-24. xgboost's own `pred_contribs`, which computes exact Shapley values
for tree ensembles -- not the `shap` package's sampling approximation, and no
extra dependency.

## A. The 19 shipped sub-scores, shifted contrast (AUC 0.8646)

    var_corr_sum            1.0771   towards POSITIVE
    var_fragment_coverage   0.5180   towards POSITIVE
    var_rt_spread           0.3243   towards NEGATIVE
    var_ms1_coelution       0.2608   towards POSITIVE
    var_xcorr_shape         0.2469   towards POSITIVE
    var_candidate_margin    0.2213   towards POSITIVE
    ...
    top 3 = 55.7% of all |SHAP|

**The scorer is three features wearing nineteen.** `var_corr_sum` alone carries
more attribution than the next four combined, and the top three account for 56%.
That is consistent with every earlier measurement: corr_sum was the only
sub-score to hold up at low abundance (0.897 at Q1 against 0.532 for
library_corr), and here it dominates the attribution too.

**`var_rt_spread` is third, and it points towards NEGATIVE.** The feature added
three days ago is doing real work -- but on this contrast it fires as evidence
AGAINST, which is what it should do: a shifted window's picked candidate is the
most peak-like noise available, and its fragments' retention-time centroids
disagree. That is the feature behaving exactly as designed, and it is the first
independent confirmation that it measures what its comment claims.

**The three library-comparison features are weak here** (`library_corr` 0.1378,
`library_dotprod` 0.0973, `library_rmsd` 0.0953, together 15% of corr_sum's
attribution). That is worth holding next to doc/60: giving the NETWORK
per-fragment library intensities was worth +0.0204 AUC, while the GBT's three
scalar summaries of the same information contribute little. The information is
valuable; the scalar summary of it is not. That is the summarisation hypothesis
stated as an attribution rather than as an AUC.

## C. Descriptors only, entrapment contrast -- the 0.7060 taxonomy floor

    by block:   relint 41.6%   prodmz 28.3%   series 14.5%
                ordinal 13.2%  frcharge 2.4%  mask 0.0%

    relint[0]   0.5041   towards NEGATIVE   (single largest contributor)

**41.6% of the taxonomy leak is library relative intensity, and most of that is
one slot** -- `relint[0]`, the first fragment's expected share. doc/57 measured
the marginal difference (top-fragment share 0.2404 human against 0.2165
Arabidopsis) and this confirms it is the dominant route, not merely a detectable
one.

Product m/z adds 28.3% and series plus ordinal another 27.7%: exactly codex's
doc/56 A2 argument that product m/z + series + ordinal is a partial
representation of the peptide sequence. `mask` contributes 0.0%, so the
composition matching DID succeed on fragment count -- it just never touched the
quantities that matter.

**Consequence: the entrapment leak cannot be matched away.** It lives in the
expected fragment pattern itself, which is the thing the model must see to do
its job. Removing it means removing the feature. So on the entrapment contrast
the floor can be measured and reported but not eliminated -- which is the
argument for the RT-shifted contrast being primary, now with a mechanism rather
than a preference.

## Where the combined model's gain comes from -- and it is not where I predicted

    traces + descriptors                            0.9059 +- 0.0003
    traces + descriptors + the 5 SIDE-CHANNEL       0.9065 +- 0.0011
      scalars (MS1, mass accuracy, mass spread,
      im delta, im spread)
    traces + descriptors + ALL 19 scalars           0.9219 +- 0.0005

**The side channels add +0.0006 -- nothing, inside the seed spread.**

I predicted the opposite. The argument was that those five are computed from
data the tensor does not contain at all (MS1 survey trace, per-fragment mass
deviation, ion mobility) and reach AUC 0.681 by themselves, so a gain from
combining ought to live there. It does not.

So the +0.0160 comes from the fourteen scalars computed from the SAME intensity
matrix the traces already carry. That splits further, and the halves have very
different consequences:

* **Competition information the trace model structurally cannot see.**
  `var_candidate_margin` is the distance to the precursor's SECOND-best
  candidate. The trace model sees one window and knows nothing about rival
  candidates elsewhere in the run. Not a leak -- it is exactly DIA-NN's
  competition family (`pBestCorrDelta`, `pTotCorrSum`), which the vault's
  feature matrix lists among OpenSWATH's four structural blind spots and which
  ODIA carries as a single scalar. If the gain is here, the fix is to show the
  model the competing candidates, not to bolt a GBT onto the side.
* **Engineered statistics the network has not learned to compute.** With 16,568
  training precursors, a hand-written summed pairwise correlation may simply be
  a better estimator than one learned from scratch. That is a legitimate
  inductive-bias benefit and argues for keeping the scalars as inputs.

`allint` (traces + only those fourteen) is running to confirm the gain
reproduces there. Separating the two halves needs one more step: fit a small
tree on [trace-model score, the 19 scalars] and take TreeSHAP over it, which
says which scalars add OVER the trace model rather than which are useful alone.
The trainer does not save per-example scores yet, so that needs a small change.

**On my own prediction.** It was stated before the split was run, precisely so
it could fail, and it failed cleanly. The value is that "MS1, mass and mobility
are worth adding to the tensor" is now measured FALSE on this contrast -- which
removes an expensive corpus rebuild that doc/56 A5 had listed as owed work.
