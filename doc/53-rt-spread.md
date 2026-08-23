# var_rt_spread, and the prediction stated before the measurement

2026-08-23. Added in 6b2b9bd. This section is written BEFORE the arm runs, so
the interpretation cannot be fitted to the result afterwards.

## What it is

The retention-time analogue of `var_mass_spread` and `var_im_spread`: each
informative fragment gets a background-subtracted, intensity-weighted
retention-time centroid over the candidate's own boundaries, and the feature is
the weighted scatter of those centroids about the group apex, negated.

The motivating argument (the user's): a group-level retention-time delta cannot
discriminate, because a group has ONE retention time and its traces co-elute by
construction -- both a real group and an interference group sit wherever the
signal they were built from sits. Retention time is diagnostic BETWEEN traces.
That is what killed `var_rt_delta`, and it is what this is built on.

## The reason it may be dead on arrival, stated in advance

**The candidate picker already selects for the property this measures.**
`findCandidatesByCorrelation` proposes a position as a peak only where the
fragments ALREADY correlate over the window -- that is its whole design. So
every group reaching the scorer, target or decoy, has been filtered for
fragments that agree about when they elute. A feature measuring that agreement
may therefore find both classes equally good and land at AUC ~0.5.

This is the same trap `var_rt_delta` fell into and the same one that produced
the retracted "gate selects on the score axis" claim: **when a feature measures
something the candidate FORMATION already required, the surviving population
carries no contrast.**

Two things could still save it, and they are the actual hypothesis under test:

* Correlation is shape similarity over a window; centroid agreement is stricter
  about POSITION. Two traces can correlate at 0.9 with apices two cycles apart.
* `corr_sum` is the best reference fragment's summed correlation to the others.
  A group passes on a subset that agrees, while the remaining fragments -- the
  interfering ones -- are free to sit anywhere. The scatter sees them; the
  picker's threshold did not.

## The acceptance test, fixed in advance

Not identifications. `d9_auc_by_abundance.py`, which reports AUC against the
decoy null per abundance quintile, with the `random` decoy-row mode as the
control for argmax selection.

* **Dead** if AUC <= 0.55 at every quintile. Then the picker-selection argument
  is right, the feature is removed, and that is a result about candidate
  formation worth having.
* **Interesting only if it holds at LOW abundance.** The deficit this project
  is chasing is discrimination at the faintest quintile, where
  `var_library_corr` measures 0.532 and `var_xcorr_shape` 0.598. A feature that
  is strong at Q5 and chance at Q1 adds nothing to the problem -- every existing
  shape feature is already strong at Q5.
* **Redundancy must be measured, not argued.** `var_xcorr_coelution` is a
  pairwise cross-correlation LAG, not an offset from the apex, and it collapses
  at low abundance (0.696 at Q1 against 0.961 at Q5). If `var_rt_spread` tracks
  it quintile for quintile it is the same statistic wearing a different name.

Identifications on the fixture are the SECOND test, and only if the first
passes. A feature that moves IDs while measuring nothing is a feature that has
found a leak.
