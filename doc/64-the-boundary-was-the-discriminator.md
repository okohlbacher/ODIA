# The boundary was the discriminator

2026-08-25. Phases 1-3 of the boundary work, measured; and phase 4 measured
before being written.

## What was wrong

The boundary walk descended while the smoothed trace stayed above
`0.10 * apex`. On 3,427 confident positives the local baseline is a median of
**36% of the apex**, and exceeds a tenth of the apex for **87.6%** of them. The
floor was unreachable by construction, so the walk ran to the window edge: the
median candidate group covered **129 of 130 cycles** against a fitted FWHM of
about 3 cycles.

## Why it cost identifications

Every sub-score is computed *inside* the group. If the group is the whole
window, then a candidate at the true apex and a candidate 20+ cycles away
integrate the same interval and receive the same number. The scores were blind
to position.

Measured directly -- same sub-score at the true apex and at the tallest local
maximum at least 20 cycles away, 11,728 paired candidates, AUC correct vs
wrong:

|            | median width | AUC lib_corr | AUC coelution |
|------------|-------------:|-------------:|--------------:|
| old rule   |          128 |        0.545 |         0.555 |
| new rule   |            7 |        0.769 |         0.782 |

Co-elution, old rule: correct +0.074, wrong +0.060, **median paired gap
+0.000**. New rule: correct +0.176, wrong +0.007, gap **+0.159**.

This is the funnel's "discrimination at low abundance" deficit given a
mechanism. It is also why doc/51 found ranking healthy: ranking was being asked
to order candidates whose scores carried almost no positional information.

The old rule's 0.545 is close to a tautology -- both candidates integrate the
same interval, so of course they score alike. That is the point, not a
weakness of the probe: the tautology is what the code was doing.

## The three changes

1. A shared `peakBounds()` for both pickers -- the default picker had been
   running without the guards the amplitude picker already had. Apex snap to a
   local maximum within +-2 cycles; floor, rebound (0.25*apex above the running
   minimum) and span (`peak_max_half_cycles = 20`) guards; then symmetric
   widening to `peak_min_cycles = 7`.
2. A noise-relative floor: robust baseline + `boundary_sigmas`*sigma from
   points outside a +-6-cycle core, capped at `apex - eps` so it never becomes
   unreachable in turn.
3. `localBackground` from the group's +-20-cycle flanks rather than the whole
   trace, widening until it has >=8 points and degrading to the old estimate
   rather than failing. 48.7% of per-fragment backgrounds differed by >20%
   between the two estimates; median signed difference 0.0, so this is variance,
   not bias.

## Paired replay, old rule against new

    width p50        128 -> 7        p90  128 -> 11     mean 112.4 -> 8.2
    spans >80% of the window        83.0% -> 0.0%
    true apex inside the group      99.1% -> 63.2%
    terminating on the span bound rather than the floor:  0.2%

**The 99.1% is trivial** -- a 128-cycle group contains everything. The 63.2%
looked like a regression and is not one: it is a *picking* result. ODIA's
top-ranked apex is within 3 cycles of DIA-NN's for 62.9% of these precursors and
then jumps to p75 = 22, p90 = 50 cycles. Conditioned on the candidate apex being
within 3 cycles, the new bounds contain the true apex **99.4%** of the time -- at
every k in {0, 0.25, 0.5, 1.0, 2.0}. No boundary rule can rescue a candidate
apex 50 cycles from the peak, and none should be credited or blamed for it.

That bimodality is now the largest open item in the funnel, and it is a
different phase.

## Choosing k

    k     med width   p90 width   walk<=1   apex in   apex in | d<=3
    0.00       8          19        5.5%     63.7%       99.6%
    0.25       7          16        8.6%     63.6%       99.5%
    0.50       7          14       11.9%     63.5%       99.4%
    1.00       7          11       18.7%     63.2%       99.4%
    2.00       7           8       30.6%     63.1%       99.3%

k does not trade against apex coverage (63.7% -> 63.1% across the whole range);
it trades against the width of the tail. Empirical FWHM is p50 3, p90 7 cycles,
so k = 1.0's p90 width of 11 tracks the p90 peak; k = 2.0's p90 of 8 would start
truncating genuinely broad peaks. **k = 1.0 kept, now on evidence.**

## Phase 4, measured before writing it

OpenSWATH does not do what ODIA does. `MRMTransitionGroupPicker::recalculate-
PeakBorders_` determines boundaries **per transition** and takes a consensus:
if the best transition's border is a z-score outlier against the others, it is
replaced by the median. ODIA -- including ODIA's own OpenSWATH-backed picker,
which feeds `PeakPickerChromatogram` a single summed chromatogram -- determines
them once, on the sum.

The fragments of one eluting molecule agree about where the peak starts and
stops. An interference's do not. So the agreement is both a better estimator
and a discriminant. Both parts measured, same 11,728 pairs:

    lib_corr  @ summed bounds     0.769
    lib_corr  @ consensus         0.786
    coelution @ summed bounds     0.782
    coelution @ consensus         0.801
    border dispersion (new score) 0.693

    border dispersion: correct 2.97 cycles, wrong 4.45 cycles
    corr(dispersion, coelution) on correct candidates: +0.543

Partly new information, not the same signal renamed. Modest but real, and it
arrives with a new sub-score rather than only a better estimate.

## Status of these numbers

Every figure here is a **sign detector**, not acceptance. They come from a
Python restatement of the rules over dumped traces, not from the shipped C++,
and the "wrong" candidate is constructed rather than one the picker emitted.
Acceptance remains `fdp_compare.py` at matched entrapment FDP on a full run.
The fixture suite (51 tests) passes.
