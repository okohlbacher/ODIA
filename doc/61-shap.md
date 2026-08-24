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
