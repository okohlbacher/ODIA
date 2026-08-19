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
