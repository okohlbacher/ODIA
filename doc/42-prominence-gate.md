# Replacing Gate C's global quantile with a per-precursor threshold

## The two changes

**1. A global 1/K0 offset from the CiRT seed, applied before pass 1** (`ca63a29`).
The extractor already centres its mobility window on
`lib_im + mobility_model->offsetFor(...)`, but that model is fitted from pass-1
peak groups, so pass 1 logged `DEFERRED -- this pass extracts on the library's
1/K0`. Pass 1 is where the gate first rejects 3.77 M targets, so the correction
arrived one pass after the decision it would change. Codex ranked this the #2
issue overall and stated the invariant:

> No irreversible precursor rejection until current-run mobility uncertainty is
> demonstrably small relative to the extraction window.

Measured worth, emulating pass 2's correction on a 1,000-precursor fixture:

| group | source | apex | baseline | prominence |
|---|---|---:|---:|---:|
| SHARED | DIA-NN | 16,783 | 1,144 | 0.937 |
| SHARED | ODIA uncalibrated | 5,031 | 518 | 0.908 |
| SHARED | **ODIA calibrated** | **9,020** | 693 | **0.937** |
| HARD | DIA-NN | 4,595 | 836 | 0.826 |
| HARD | ODIA uncalibrated | 2,008 | 434 | 0.793 |
| HARD | **ODIA calibrated** | **2,493** | 515 | **0.808** |

+79% apex on shared, +24% on hard, baseline only +34%, and prominence on shared
reaching exactly DIA-NN's. A single constant is deliberate -- 213 seed anchors
cannot support the per-charge m/z-shaped model, which is why the full calibration
refuses them, but they support one number, and at +0.0183 against a 0.025
half-window one number removes most of the clipping.

**2. `-gate_mode prominence`** (`ce0cf8b`). Gate C's threshold is the (1-alpha)
quantile of a run-wide decoy null, which makes one precursor's admission depend
on the library's composition and on ARRIVAL ORDER -- precursors reach it in
retention-time order, so tau is calibrated on the earliest eluters, and alpha
does not deliver its rate (0.05 nominal, 15.3% decoy admission).

`coelutionEvidence` already works in per-precursor units: it sums robust
z-scores across the contributing transitions and averages over
w = 2*gate_smooth_half+1 cycles, so under white noise its null is
N(0, contributing/w). A k-sigma cut therefore needs no null, no calibration
sample and no global state. **The statistic is unchanged; only the threshold
moves.**

k is a look-elsewhere threshold, not a significance level:
`k = Phi^-1((1-alpha_P)^(1/M_eff))`, giving 3.1 / 3.3 / 3.5 at M_eff 50 / 100 /
200 for alpha_P = 0.05. Default 3.3, and documented as the centre of a sweep.

## The test, at MATCHED empirical FDP

Nominal-cutoff comparisons are meaningless here -- the run reports 1% and the
entrapment population says 12-50% -- so both arms were walked down their own
ranked lists to equal measured FDP. 1,011,145-precursor library carrying 4,000
Arabidopsis proteins among 20,416 human (r = 0.1723), `-passes 1`, identical
otherwise.

| target FDP | Gate C | prominence k=3.3 | gain | Gate C recovers | prominence recovers |
|---|---:|---:|---:|---:|---:|
| 5.0% | 555 | **1,103** | **1.99x** | 84 | **480** |
| 10.0% | 659 | **1,242** | **1.88x** | 99 | **529** |
| 15.0% | 732 | **1,368** | **1.87x** | 102 | **554** |
| 23.7% | 1,693 | 1,734 | 1.02x | 253 | 614 |

"Recovers" counts the 11,145 precursors DIA-NN identifies that the production run
never formed a candidate for.

**~1.9x the identifications at the same empirical FDP through the useful range,
and 5.4-5.7x the recovery of the previously-lost precursors.** The arms converge
at 23.7% because by then both admit essentially everything they can.

At each arm's own nominal 1% the picture is the misleading one: 1,650 IDs at
23.7% FDP against 4,955 at 49.9%. Reported only to show why that comparison must
not be used.

## What this does NOT establish

- **One instrument, one fixture.** The library is 98% ballast and its FDP scale
  is inflated relative to a real search. doc/13's rule is two instruments;
  Astral is unrun.
- **It is a proposal rule, not a presence test.** Codex's strongest objection to
  the original design was semantic, and it stands: peak-shaped interference in a
  25 Th window passes at any k. Discrimination has to come downstream from
  co-elution, library agreement, mass accuracy, RT and mobility.
- **The full-library effect is unmeasured.** These arms are `-passes 1` on a 1 M
  library; the production run is two passes on 9.6 M.
- **k is not tuned.** 3.3 is the M_eff=100 look-elsewhere value, tested against
  one alternative that did not finish.

## Not adopted, recorded as the next step

Codex's preferred primitive is candidate times from the UNION of per-transition
maxima rather than from the summed trace: a summed trace can dilute one clean
fragment, shift the apex between incompatible transition peaks, or let a coherent
interferer dominate. That is OpenSWATH's actual rule (verified 3-0 against the
OpenMS 3.6 source: a detectable local maximum in at least ONE detecting
transition, no fragment-count or correlation gate at formation). It is a larger
change and is untested here.

---

# Sub-score audit, and a negative result on var_log_sn

## Per-feature discriminative power, on KNOWN-PRESENT precursors

> **CORRECTED.** The first version of this table used an AUC with ORDINAL ranks
> from a stable sort, which breaks ties in favour of the first array and biases
> heavily-tied features DOWNWARD. Codex caught it via an internal inconsistency
> -- var_usable_fragments was reported at AUC 0.353 while its distribution had
> present >= decoy, which cannot both be true. Recomputed with average ranks:
>
> | feature | corrected | as first reported |
> |---|---:|---:|
> | log_sn | **0.502** | 0.004 |
> | im_spread | **0.500** | 0.000 |
> | usable_fragments | **0.620** | 0.353 |
> | fragment_coverage | **0.823** | 0.781 |
> | peak_width_ratio | 0.337 | 0.241 |
>
> So log_sn and im_spread are not INVERTED, they are exactly uninformative --
> which is what a constant should score. usable_fragments is mildly INFORMATIVE
> and was never dead. fragment_coverage is the second-strongest feature, not the
> third. Untied features (corr_sum, xcorr_*, library_*, rt_delta) are unchanged
> to three decimals, and the headline -- classifier 0.901 against corr_sum 0.906
> -- is unaffected, DScore being effectively tie-free.
>
> Corrected ranking: corr_sum 0.906, fragment_coverage 0.823,
> xcorr_coelution 0.811, xcorr_shape 0.764, library_rmsd 0.747,
> library_dotprod 0.708, library_corr 0.688, intensity_score 0.630,
> usable_fragments 0.620. Dead: log_sn 0.502, im_spread 0.500, ms1_coelution
> all-NaN. Lower-is-better and fine: rt_delta 0.311, peak_width_ratio 0.337.

Measured on the 14,081-ID full run, positive class = DIA-NN's 21,055 confirmed
present precursors (using all targets gives ~0.5 for everything -- 99% of them
are absent, which is the prevalence trap this project keeps falling into).

| strong | AUC | weak / inverted | AUC |
|---|---:|---|---:|
| **corr_sum** | **0.906** | mass_spread | 0.553 |
| xcorr_coelution | 0.810 | candidate_margin | 0.548 |
| fragment_coverage | 0.781 | im_delta | 0.513 |
| xcorr_shape | 0.764 | mass_accuracy | 0.490 |
| library_rmsd | 0.747 | yseries_score | 0.465 |
| library_dotprod | 0.708 | usable_fragments | 0.353 |
| library_corr | 0.688 | rt_delta | 0.311 |
| intensity_score | 0.630 | peak_width_ratio | 0.241 |

`rt_delta` and `peak_width_ratio` are lower-is-better and therefore fine.

**The 19-feature classifier scores AUC 0.901 -- WORSE than corr_sum alone at
0.906.** Four features explain why:

| feature | state |
|---|---|
| `ms1_coelution` | all-NaN; the CiRT seed's 708-precursor MS1 traces persisted (fixed, 9a496fa) |
| `log_sn` | constant: TEN distinct values over 21,055 present precursors, p10 = p50 = p90 = log(100) |
| `im_spread` | literal constant: ONE distinct value across both classes |
| `usable_fragments` | pinned at 12, as its own comment admits |

## var_log_sn: un-saturating it changes nothing

`floor_bg = max(background, frac * apex)` caps signal-to-noise at 1/frac, and
every real peak exceeds 100:1, so the historical 0.01 pinned the feature. Swept
on the 1M fixture:

| frac | cap | distinct var_log_sn values | IDs |
|---|---:|---:|---:|
| 0.01 (default) | 100:1 | 2,406 | 1,650 |
| 0.001 | 1000:1 | **11,023** | **1,644** |

The feature genuinely un-saturated -- 4.6x the distinct values -- and
identifications did not move (-6, 0.4%). **Default left at 0.01.**

The lesson is about prioritisation: a feature can be visibly broken and still not
be a lever, because its information is already carried by correlated features
(intensity_score, corr_sum) or because the classifier cannot exploit it. Repair
work on the remaining dead features should be ordered by measured effect, not by
how broken they look. On that basis `ms1_coelution` is the one worth having --
doc/17 records it at 13.7x enrichment in its top bin, and unlike log_sn it was
returning NaN rather than a constant, so the classifier never saw it at all.

---

# Kimi's review: the sub-score comparison was circular, and corr_sum is load-bearing four times over

## The evaluation was selection-biased

The positive class was "precursors DIA-NN confirmed". DIA-NN's acceptance is
dominated by co-elution correlation, and `corr_sum` is our closest relative of
that statistic -- it IS the picker's own detection statistic. So corr_sum's 0.906
partly measures how well corr_sum predicts what a co-elution-based engine
accepts, not how well it detects presence. **Any feature that could recover
precursors DIA-NN MISSED is penalised by this metric**, and the 0.901-vs-0.906
gap is close to a tautology.

Having escaped the prevalence trap, the audit walked into a selection trap.

## corr_sum enters the pipeline FOUR times

Verified in code by kimi:

| # | as | where |
|---|---|---|
| 1 | the candidate gate | `min_corr_score`, `PeakGroupScorer.cpp:520` |
| 2 | the feature `var_corr_sum` | `:1127` |
| 3 | `candidate_margin` -- literally its first difference | `:1477-1495` |
| 4 | **the semi-supervised loop's seed mask** | `seed_mask` admits only CORR_SUM, `:1653-1657` |

Iteration-0 positives are targets whose corr_sum-ranked row clears q <= 0.15, so
**the positive class the model trains on was selected by corr_sum**. The loop
does not fail to ignite -- it ignites onto its own seed and re-expresses it. That
is a sharper diagnosis than "label contamination" and it is structural.

Kimi also puts effective dimensionality at ~5-6, not 19: library_corr,
library_dotprod, library_rmsd, intensity_score and yseries_score all read the
same background-corrected area vector.

And it exonerates the dead features: constant and all-NaN columns are zeroed in
place (`fitAndAssign_ :1551-1584`) and LDA z-scores with NaN -> 0 = mean
(`lda.h:412-461`). They are inert, not noise.

## var_log_sn: closed, negative

Kimi's disambiguator -- standalone AUC of the UNSATURATED feature against
confirmed-present precursors:

| config | distinct values among present | standalone AUC |
|---|---:|---:|
| 0.01 (saturated) | 7 | 0.500 |
| 0.001 (unsaturated) | 61 | **0.494** |

8.7x the dynamic range, AUC unmoved at ~0.5. By the stated criterion the feature
is genuinely weak and the null result was foreordained -- it is NOT independent
evidence of the distillation problem. **Default stays at 0.01 and log_sn needs no
further work.**

## The two experiments that would settle the rest

1. **Ablation.** Retrain with `corr_sum` AND `candidate_margin` excluded.
   Drops to ~0.65 -> the other features genuinely carry little against these
   labels. Stays >= 0.89 -> the loop was merely picking corr_sum first and the
   classifier was never the bottleneck. Then re-run the winner against
   ENTRAPMENT positives rather than DIA-NN's set, which removes the selection
   bias entirely.
2. **The swap.** Feed DIA-NN's `--xic` traces through OUR scorer, using the
   aligned 1,000-precursor dump that already exists. If ODIA-scored-on-DIA-NN-
   traces recovers the identification gap, extraction limits us; if it does not,
   scoring does.

## A concrete extraction mechanism kimi found

The mobility correction is clamped at +/-0.030 while the extraction half-window
is +/-0.025, and `ChromatogramExtractor.h:519-521` records that +/-0.025 rejects
**28.3% of target anchors even when calibration is correct**. The window is
centred on the PRECURSOR's corrected 1/K0, identically for every fragment
(`:986-993`). Fragments whose true 1/K0 sits off the precursor's are therefore
asymmetrically clipped -- which halves intensity while preserving prominence
shape, and bites the HARD group first. That is exactly the pattern measured:
apex 0.54x of DIA-NN's, prominence matching on SHARED (0.937) but not HARD
(0.808 vs 0.826).

Testable from data already on disk: `collect_im_residuals` is on by default, so
the per-fragment 1/K0 residual distribution can be read directly. Wider than
~0.01 confirms the mechanism and implies the fix -- widen the window, or centre
it per fragment rather than per precursor.
