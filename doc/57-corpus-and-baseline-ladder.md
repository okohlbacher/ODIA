# The corpus, and the first two rungs

2026-08-24. doc/55 planned, doc/56 revised it after codex, kimi and vibe.
This is what the corpus looks like and where the ladder stands.

## The corpus

    dump                215,342,936 rows -> (142,042, 12, 128) float32, 0.87 GB
    shifted arm         the same, row-aligned (asserted, not assumed)
    cycles              130 per trace at 1.385 s, padded to 128
    calibration         PINNED at +-10 ppm about -10.7372 ppm, production's own

    label          n        role
    pos       33,494        DIA-NN confident AND apex inside the window
    unlabelled 5,472        confident but apex OUTSIDE -- A11, excluded
    ent       20,086        entrapment DIA-NN does not call
    ent_hard      17        entrapment DIA-NN DOES call: verified false positives
    unid      11,952        human targets DIA-NN does not call (unverifiable)
    decoy     71,021        held out entirely

`ent_hard` at 17 is **too small to measure the tail**, which was its entire
purpose (A3). Composition-matched sampling drew 20k entrapment out of 733k and
caught almost none of the ~2,900 DIA-NN calls among them. Fixing it needs those
precursors added explicitly and re-extracted; until then the hard-negative
question is unanswered, not answered negatively.

## The paired design works, and the proof is a control that returned exactly 0.5

The primary contrast is a positive against **the same precursor's own window
shifted 300 s** -- kimi's proposal. Same sequence, charge, product m/z, series,
ordinals, library intensities; only the window moves.

    true-metadata control       AUC 0.5000     identical between arms by pairing
    label-shuffle control       AUC 0.4954
    intensity+coverage floor    AUC 0.6939     trace-derived, so a floor not a control

The first line is the point. With entrapment negatives, "no metadata leakage"
would have been an argument resting on composition matching, and the measured
library-intensity difference (top-fragment share 0.2404 human against 0.2165
Arabidopsis) shows that argument would have been wrong. With shifted negatives
it is arithmetic.

**A control caught my own error before any conclusion rested on it.** The first
pass called intensity and coverage "metadata-only" and predicted 0.5; it
returned 0.594. Those are TRACE-derived and legitimately differ between arms, so
the control was mislabelled, not the design. Separated, both now behave as
predicted -- which is the only reason the 0.5 above can be believed.

## The ladder

    RUNG 1  shipped 19 sub-scores + GBT     AUC 0.8649
    RUNG 2  per-fragment order stats + GBT  AUC 0.8422
    RUNG 3  Deep Sets                       not run (no torch on dax)
    RUNG 4  transformer                     not run

**Rung 2 does not clear rung 1.** Both reviewers named order statistics as the
real bar a transformer must beat; my implementation of them is worse than what
already ships, so the bar is 0.8649 and the ladder has not been climbed.

The comparison is not clean in rung 2's favour, and the asymmetry runs the other
way from what would flatter it: **rung 1 scores ODIA's PICKED candidate** -- a
localisation decision already made by the picker -- while rung 2 gets the raw
window and must find the peak itself. Rung 2 reaching 0.8422 from raw traces
against 0.8649 from a picked candidate plus nineteen engineered scores is
respectable. It is still not a win, and it is reported as a loss.

`var_rt_spread`, added two days ago, is the 6th most important of the nineteen
(0.045). On a contrast with a localisation flavour that is unsurprising and is
not evidence for it in production.

    abundance quintile   AUC (rung 2)
      Q1    1,009- 13,032   0.7627
      Q2   13,032- 22,884   0.8072
      Q3   22,884- 41,346   0.8386
      Q4   41,346- 94,440   0.8739
      Q5   94,440-13.3M     0.9239

## What these numbers are NOT

**They are not comparable to the production sub-score AUCs** (`var_library_corr`
0.532 at Q1 and so on). Those measure a target against the DECOY NULL; these
measure a precursor's own window against its own window 300 s away. The second
task has a localisation component the first does not. Quoting 0.7627 against
0.532 would be comparing two different questions, and the temptation to do it
is exactly why this section exists.

The shifted contrast answers "is this peptide HERE rather than there". The
production task is "is this peptide present at all". They overlap and are not
the same, which is why the entrapment contrast is kept despite its confounds --
and why `ent_hard` being 17 matters.

## Next

1. `ent_hard` properly populated -- the tail is currently unmeasured.
2. Rungs 3 and 4 need torch, so a GPU node (spock or data), not dax.
3. The entrapment contrast run alongside the shifted one, since they answer
   different questions and only their disagreement is informative.
4. IH2 stays locked. It has been converted and not touched.
