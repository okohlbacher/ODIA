# The morning full benchmark, 2026-08-22

Two full-scale arms launched 08:15, S08 diaPASEF, our 9.6M library.

## v9: var_rt_delta removed. A wash, confirmed at full scale.

    arm                       IDs     wall      peak
    full v5 (production)   13,268   3:00:00   124 GB
    full v9 (HEAD)         13,301   3:26:12   138 GB

**+33 identifications, +0.25%.** The fixture measured +20; the full run measures
+33. Both are inside noise, and they agree, which is the useful part: the
fixture predicted the full-run outcome for this change.

So the removal is neutral, not positive. It stands on the argument that
motivated it -- a peak group has ONE retention time and its traces co-elute by
construction, so a group-level delta cannot separate a real group from an
interference group -- and on the measurement that the feature was
uninformative (AUC 0.215). It was never worth identifications either way.

At MATCHED empirical FDP (`fdp_compare.py`), v9 against v5:

    FDP target      v5 IDs    v9 IDs    delta   entrapment   sigma
        2.00%        1,014     1,890   +86.4%        3 / 6   1.41 pp
        3.00%        1,782     2,628   +47.5%       9 / 13   1.30 pp
        5.00%       12,159    11,264    -7.4%     105 / 97   0.70 pp
        5.72%       13,140    12,342    -6.1%    130 / 122   0.72 pp
        7.42%       14,484    14,396    -0.6%    186 / 184   0.77 pp  <- DIA-NN
       10.00%       15,700    15,860    +1.0%    271 / 274   0.86 pp
       15.00%       17,916    17,782    -0.7%    465 / 461   0.99 pp

**Neutral at the operating point.** The +86% and +47% at 2-3% rest on 3-13
entrapment hits at sigma ~1.4 pp and are noise -- the tool prints the counts
beside every cell precisely so those cannot be quoted as a win. The 5-5.72%
band shows -6 to -7% on better statistics (sigma 0.70-0.72 pp) and is the one
cell that could be a small real deficit; it is not excluded and not confirmed.

Taken with the nominal +33 IDs, the removal is neutral at full scale.

### Correction: removing it did not "retire the refine loop"

I wrote that RT_DELTA's removal made the retention-time refinement loop vacuous
and called that a real simplification. `-refine_rounds` DEFAULTS TO 0, so the
loop was already not running in the production configuration, and the
`refitsChangeScores()` guard I added never fires in a default run -- the
`max_rounds == 0` early return precedes it. The guard is still correct and
still the right place to record that no sub-score reads the map. It is not a
simplification of anything that was running.

## prom9: refused at the CiRT seed, and the refusal is a result

    CiRT seed REFUSED: the decoy control fits as well as the targets
    (26.6 s vs 25.7 s, needs > 32.1 s)

Exit 13 after 5:33, zero identifications. It never reached the gate it was
meant to test.

This is the SAME failure v8 hit with `shuffle` decoys, and the same mechanism
the fixture measured for prominence: admit more, and the decoy population
acquires enough structure to be fitted to a retention-time line as well as the
targets can. Null inflation showing up at a third stage -- first the q-value
threshold, then candidate multiplicity, now the seed's own control gate.

The workaround is documented and is what v8 used: supply the map directly
(`-irt_slope 1086.50 -irt_intercept 473.77 -rt_window_pass1 75.4069`) instead
of `-rt_seed cirt`. The arm needs relaunching that way. **The doc/46
contradiction is therefore still unsettled.**

## The full-run terminal-reason table, with an unexplained discrepancy

    few_points              1,119,490
    gate_c                  6,574,288
    scored                  1,556,317
    no_candidate                1,129
    all_candidates_dropped        397
    (not_reached ~13; no_window_coverage and prefilter_excluded absent)

Gate C rejects 71% of the library and 16.8% is scored, against 77.3% and 3.9%
on the fixture -- the full run is the harsher regime, as expected.

**But `no_window_coverage` fired ZERO times while `few_points` equals
1,119,490, which is exactly the count the same log reports as "precursors
covered by no isolation window".** That coincidence needs explaining before any
of these numbers is used to attribute a stage. The extractor's `covering == 0`
branch `continue`s without emitting, so those precursors should never reach
`Session::add` and should never be marked `few_points` -- yet the two counts
match exactly.

RESOLVED, from the per-precursor file rather than from the coincidence.

Of the 559,745 TARGET precursors marked `few_points`, **100.0% lie outside every
isolation window's m/z range** and 0.0% lie inside it. Median m/z 1468.3 against
525.8 for scored precursors. The reason: **S08's 24 windows stop at 1400.62 Th**
while the library is generated out to `-precursor_mz_max 1800`.

So `few_points` was the uncovered population wearing the wrong label, and
`no_window_coverage` -- the label that exists for exactly this -- never fired.
The scorer's `mark()` was overwriting the extractor's more specific reason with
its own vaguer one, which is also true and much less useful. The extractor's
reasons now take precedence, and the fixture arm testing that is running; if
`few_points` does not move, the precedence fix was the wrong diagnosis and the
extractor is not marking at all.

### 12.8% of the target library cannot be searched on this file

    targets outside the acquired m/z range   641,481 of 4,991,901   12.8%
    of those, marked few_points              559,745
    scored targets outside the range         0 of 788,444

This is an acquisition/library mismatch, not an algorithm loss, and it is
symmetric across target and decoy (few_points totals 1,119,490 = 2 x 559,745).
It does not distort the FDR -- these precursors produce no candidates, so they
never enter the null -- but generating and extracting them is wasted work, and
any per-library-precursor rate quoted without excluding them is diluted by an
eighth. Restricting library generation to the run's acquired range is a cheap
correctness improvement rather than a scoring one.

The doc/51 fixture conclusions are unaffected: `few_points` accounted for 14 of
DIA-NN's 4,948 confident precursors there (0.3%), so relabelling it does not
move the funnel.
