# Gate C is what rejects the precursors we never extract

> **CORRECTION 2026-08-18, same day.** This document first concluded that
> candidate formation "collapses with library size". That framing was WRONG and
> the experiment that produced it was confounded. `NullCalibration::admit`
> (`PeakGroupScorer.cpp:753`) admits EVERY precursor unconditionally until it has
> collected `gate_calibration_n = 20000` decoy statistics:
>
> ```cpp
> if (!ready) { if (is_decoy) decoy_stats.push_back(stat);
>               if (decoy_stats.size() >= n_needed) { ...fix tau...; ready = true; }
>               return true; }        // admit while the null is being built
> return stat >= tau;
> ```
>
> The isolated run had 11,144 decoys -- BELOW 20,000 -- so Gate C never engaged
> and admitted everything. The 94.8% was a run with no gate, not a run at small
> scale. Library size is a proxy for "did the calibration set fill", nothing
> more.
>
> The real finding is below: **Gate C rejects these precursors**, because their
> co-elution evidence sits below the (1-alpha) quantile of the run's decoy null.
> The sections after this correction are kept because their negative results
> still stand.

2026-08-18. **This overturns doc/37's conclusion.** The precursors we fail to
extract are not too dim to extract. They extract fine in isolation and stop
extracting as the library grows.

## The experiment

Take the 11,145 precursors DIA-NN identifies on IH1 that our full run never
formed a single candidate for. Run them through the SAME binary, the SAME raw
file, the SAME supplied map (`-irt_slope 1095.00 -irt_intercept 453.49`), the
SAME window (`-rt_window 92.4`), the SAME thresholds, `-passes 1`. Vary only how
many other precursors are in the library.

| library | total precursors | of the 11,145 forming a candidate |
|---|---:|---:|
| the missing set alone | 22,289 | **10,566 (94.8%)** |
| padded with 1 M others | 1,011,145 | **1,518 (13.6%)** |
| the production run | 9,617,705 | **0 (0.0%)** |

In isolation these precursors are not marginal at all:

```
split_miss:  no points (<3)  0/0        empty trace  0/0
             precursors reached  11144/11144 (1.00x)
             30,431 target groups from 11,145 precursors  -- 2.7 each
             1,182 (10.6%) yielded no candidate
```

**Zero empty traces. Zero short traces. Every one reached the picker.** The
signal is there, the mass calibration finds it (doc/37: 0.25 ppm), the RT window
contains it (doc/37: 75% inside +/-92.4 s).

## What this overturns

doc/37 concluded the loss was sensitivity, on the evidence that the missing
precursors are 2.0-3.3x less abundant than the found ones. That correlation is
real but it is a CONSEQUENCE, not the cause: when candidate formation degrades
with library size, the dimmest precursors are the ones that drop out first. The
abundance gap is what a scale-dependent failure looks like from the outside.

Ranked reasons from doc/37 SS1 -- abundance first, RT second, charge 1 third --
should be read as ranked SYMPTOMS. The cause is upstream of all of them.

## What it is not

Ruled out by the logs of the two runs:

- **Not the memory cap or chunking.** Both runs report
  `bound by retention-time overlap`, never `bound by precursor cap`, and neither
  chunked. `live(n_slots)` is sized to all assignments
  (`ChromatogramExtractor.cpp:698`).
- **Not isolation-window coverage.** 12.4% of the padded library is uncovered,
  matching the library-wide 12.5%, and every one of the 11,145 is inside the
  327.5-1400.6 Th span (doc/37).
- **Not the retention-time map.** Both runs were given the identical map, and the
  signed residual difference between found and missing precursors is 6.5 s.
- **Not the pass structure.** Both sub-runs used `-passes 1`, so the pass-1
  subsetting stride cannot be the differentiator between them.
- **Not abundance.** Same precursors, 94.8% -> 13.6% by changing only the count
  of OTHER precursors in the file.

The `empty trace 369,097/529,262` in the padded run cannot be attributed: that
count is dominated by the ~518 k ballast precursors, which are genuinely absent
from the sample. Localising the stage for the 11,145 specifically needs
per-precursor instrumentation, which does not exist yet.

## What the literature says our gates should be

Deep-research, verified 3-0 against the OpenMS 3.6 source we build against
(`MRMTransitionGroupPicker.h`) and the installed binary's `-write_ini` defaults:

> OpenSWATH picks peaks INDEPENDENTLY in each fragment chromatogram, then merges
> them into a consensus feature seeded by the single largest peak. **A detectable
> local maximum in at least ONE detecting-transition chromatogram is the only
> hard requirement.** There is NO minimum co-eluting-fragment count and NO
> minimum library-spectrum correlation at candidate formation -- library
> correlation, dotprod/manhattan and MS1 isotope scores are all computed
> downstream as SCORES.

Its optional candidate-stage quality gate (`compute_peak_quality`,
`minimal_quality`) is OFF by default in both the class and OpenSwathWorkflow.
MS1 evidence is extracted by default but `use_precursors` defaults to false, so
MS1 cannot seed a candidate -- MS1 rescue is a ranking effect, not a formation
effect.

ODIA by contrast applies `min_corr`, `<2 fragments`, `apex_evidence` (default
0.99) and `not a local max` as FORMATION gates. We reject at a stage where the
reference implementation only scores. That is a design divergence worth having
deliberately rather than by accident -- but note it does not by itself explain
the scale dependence, since those gates are per-precursor.

The research returned NO verified claims on DIA-NN, Spectronaut, EncyclopeDIA,
Skyline or DIA-Umpire candidate thresholds, on smoothing's effect on weak-peak
detection, or on whether deep-learning rescoring acts at formation or ranking.
That evidence base is empty, not settled -- do not cite it either way.

## Next: instrument, then fix

The single measurement that localises this is a per-precursor "furthest stage
reached" record over a library-size sweep (22 k / 100 k / 1 M / 10 M) on the same
11,145. Everything else is speculation until that exists.

Only then is it worth touching thresholds -- and any threshold change must be
paired with the entrapment FDP (doc/36: 11.76% at a nominal 1%), because
loosening candidate formation without fixing the decoy null would inflate a
number that is already 12x optimistic.
