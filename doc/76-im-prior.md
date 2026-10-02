# The 1/K0 prior: why the library's mobility is corrected at run time

## The defect

ODIA's generated libraries carry a predicted 1/K0 that sits systematically LOW
on the vendor-calibrated axis. Measured on the 1,614-precursor IH1 consensus
(DIA-NN AND OpenSWATH agree a peak exists; membership is peak-evidence-based,
not mobility-based):

    observed - library:  median +0.0212, sd 0.0144   (+2.1% overall)
    by 1/K0:  +3.1% at 0.6-0.9  ->  +0.9% at 1.1-1.5  (affine, not a scale)

Three wrong attributions preceded the right one, all retracted in the record:
the Mason-Schamp constant (alphabase uses the same one to 0.0046%); a "94%
removed by rescaling" argument (any ~2% factor absorbs a ~2% bias -- it sized
the defect, never identified it); and "the peptdeep model is worse than
DIA-NN's".

## The split, and the calibrated-vs-raw trap

The last attribution fell to a provenance check (kimi's question): DIA-NN's
report column `Predicted.IM` differs from DIA-NN's own predicted library by
+0.0096 (0.2% identical) -- it is POST-RUN-CALIBRATION, not the raw model.
The honest table:

    peptdeep raw          +0.0212   sd 0.0144
    DIA-NN raw            +0.0082   sd 0.0176    <- also biased, same sign, WORSE scatter
    DIA-NN calibrated     -0.0028   sd 0.0141

Both independent predictors sit low on this instrument. DIA-NN reaches
"unbiased" by calibrating its predictions per run. So the bias splits into
**~+0.013 peptdeep-vs-DIA-NN model offset** (static) and **~+0.008
instrument-convention offset** (per-run calibration territory -- exactly what
DIA-NN does). peptdeep's scatter is the BETTER of the two; its centring is not.

Never compare a calibrated quantity to a raw one when arguing about model
quality; search engines recalibrate predictions silently.

## Why the runtime stage alone cannot absorb it

* Pass 1 extracts and harvests its calibration anchors AT the biased
  coordinates: the bias (0.021) eats 84% of the ±0.025 half-window, and a
  harvest truncated by the error it must fix recovers <40% of a scale error.
* The fragment-mass probe gates on library 1/K0 at ±0.010 -- the bias is 2.1x
  that, and this coupling once flipped the mass gate and cost 414 IDs before
  the auto-centre fix.
* Charges below the anchor floor are left uncorrected entirely.
* Measured: pre-correcting the library at a matched mass regime gained
  +14.9%/+16.6% over letting the runtime stage cope (two regimes, n=1 each).

## The fix (user decision: run time, not the generator)

`-im_prior` (default `timstof`) applies a per-charge affine to the library's
1/K0 at the head of the score workflow -- mirroring the external-iRT block --
BEFORE anything extracts on it:

    z2:  im' = 0.94060*im + 0.08319     (affine; the only identifiable slope)
    z3:  im' = im + 0.02445             (constant -- the slope is NOT
                                         identifiable: three anchor sets give
                                         0.981 / 0.941 / 0.910, and the
                                         affine's OOF gain fails a 0.0005
                                         materiality floor)
    z4:  im' = im + 0.02868             (constant; gain CI includes 0)
    z1, z>=5: identity                  (n=40; predictor slope 0.48 -- ship
                                         nothing where nothing was measured)

* Fitted on `shared_rt_coords_ih1.parquet` -- 39,115 q<=0.01 anchors, 24x the
  consensus set. y is the vendor axis at the reported apex, verified
  library-invariant (median |delta| 0.0009 over 28,463 cross-search
  precursors); x verified identical to ODIA's library IM to 0.000000. An
  independent re-fit reproduced every coefficient and every model-selection
  decision exactly. The consensus-set fit (z2 0.95482/0.06674) survives as a
  cross-check: it agrees with the shipped z2 to <=0.006 across the fitted
  range, and its z3 disagreement is precisely the non-identifiability that
  demoted z3 to a constant. Fitting is Theil-Sen, NOT OLS -- with noise on
  both axes OLS attenuation manufactures exactly this slope<1/intercept>0
  shape from nothing.
* Honest precision: the z2 slope carries a ~0.011 set-to-set systematic
  (~+/-0.006 applied shift at the range edges). The digits are
  reproducibility, not accuracy.
* Applied to targets AND decoys identically -- correcting one class alone
  breaks target-decoy exchangeability and with it the FDR (both reviewers'
  severity-1 finding).
* Outside the fitted range [0.742, 1.367] the boundary's correction is carried
  (constant extrapolation) and the count logged.
* On a run with no mobility axis every consumer of library 1/K0 is inert
  (Astral's own log drops var_im_delta as no-information), so the prior is a
  no-op there; the validation includes a byte-identity proof.
* One instrument (IH1 fit, IH2 direction check -- IH2's own runtime estimate is
  harvest-censored and UNDERSTATES its bias). The profile is named `timstof`
  and is provisional until a second instrument is measured.
* The runtime stage keeps its job: after the prior it fits the RUN residue.
  Success criterion is that it finds ~0.006-0.008 -- if it still finds ~0.02,
  the prior did not take. Two of its knobs became defaults on measurement:
  pass-1 harvest x2 (MSE removed 56.6% -> 67.1%) and anchor floor 100 (+5.4 MSE
  points on thin-charge-3 libraries, bit-identical elsewhere).
* Rollback: `-im_prior off` restores the previous behaviour exactly.
* NOT in the generator: libraries stay portable predictions; the correction is
  an instrument property and lives with the instrument-facing tool. CCS
  transfer learning stays deferred -- it risks baking one vendor's axis into a
  supposedly instrument-independent quantity.

## Validation (v2 profile, binary 3b6c4374, pre-registered rules)

**IH1 / mix10k (DIA-NN-confident-enriched, the design population): PASS on both
clauses.** Mechanism: mass-probe centre 0.0000 (on) vs 0.0242 (off); runtime
MSE-removed 2.3% vs 70.0%; z3 residual +0.0053 = the predicted run part.
Counts at matched entrapment FDP, e>=10 cells: +180 / +212 / +193 at
FDP<=5/7.5/10% against a 99-ID band.

**IH1 / 500k random (stress case): EQUAL within band.** +93 / +99 / -49 at the
same cells; probe centre improves 0.0430 -> 0.0142 but misses the <=0.012 rule.
Known weakness, recorded as open: this library is z4-heavy (213,864 rows under
a constant fitted on 1,194 anchors of a different population) and 27.5% of its
rows sit outside the profile's fitted [0.759,1.367], carrying the boundary
correction. The prior arm's PASS 1 halves (964 -> 482) while pass 2 and FDP are
neutral -- the runtime stage compensates. The profile's support does not cover
libraries of this range; a range-extended refit (or per-run fallback fit) is
the successor item.

**Astral (no mobility axis): byte-identical** score tables with the prior on
and off (cmp, an_IDENTICAL) -- the no-op proof, as required.

76/76 ctest green. Shipped at a2a5506.
