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

Entrapment FDP for v9 is not yet computed; it needs the output on /scratch and
the Kerberos ticket expired before it could be read.

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

I am NOT diagnosing this from the coincidence. Twice this session a mechanism
inferred from matching counters turned out to be wrong, which is the whole
reason the per-precursor table exists. Resolving it needs
`reasons_full_v9.tsv`, which is on /scratch and currently unreachable.

Until then the fixture table (doc/51) stands as the measured attribution and
this one is quarantined.
