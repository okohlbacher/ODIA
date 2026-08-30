# The mass-gate fallback costs more than anything it protects against

> **RETRACTION OF THE RECOMMENDATION, 2026-08-28 09:05 -- with its own evidence
> corrected 2026-08-29.** The section below headed "Recommendation" proposed
> making `-mass_width_from_ids apply` the default. **Withdraw it.** But the
> counter-evidence cited below is STALE and must not be used again:
>
> **The 4,969 -> 2,422 figures are from the RETIRED Mag-Net Astral file, not from
> the current `astral_neat_e1`.** `shared/bench/PROVENANCE.md` and
> `shared/ref/diann/README.md` both state that every Astral number measured
> before 2026-08-17 was measured on Mag-Net and does not carry over; the
> benchmark changed identity that day. Worse, 4,969 is not a configuration at all
> -- it is a point in the fragment m/z window sweep (15/30/50/75/100 ppm ->
> 4,290/4,499/**4,969**/4,382/3,967), and 2,422 is that same arm with the window
> sized from its own pass-1 identifications, on a 10,891-precursor pre-selected
> library. **No `-mass_width_from_ids apply` arm has ever been run on
> `astral_neat_e1`.** The correct DIA-NN reference for the current Astral file is
> **7,462**, not 12,308.
>
> So the retraction's CONCLUSION stands -- a default must not change on one file,
> and this one had one -- but it stands on that principle alone, not on a
> cross-instrument reversal that was never measured on the current benchmark. The
> honest position is that the second instrument is UNTESTED for this flag, and
> `shared/libv2/run_astral.sh` now exists to test it.
>
> The original text of this notice follows, with its stale citations left visible
> rather than quietly edited:
>
> * `doc/26-PLAN.md:58` — "`-mass_width_from_ids` defaults to `measure`
>   **deliberately**; `apply` was measured to take 4,969 -> 2,422"
> * `doc/16-phase34-plan.md:316` — "centre only. Do NOT infer WIDTH from
>   accepted IDs (4,969 -> 2,422)"
> * `doc/18-joint-reference-map.md:46`
>
> On Astral the same flag **halves** identifications. The document quoted the
> option's own "opt-in until it has been measured on both instruments" caveat
> and then recommended overriding it, without searching for whether that
> measurement existed. It did, in three places, and it says the opposite.
>
> **The measurement in this document stands; the inference from it does not.**
> +59.1% on S08 is real and reproducible. It is a fact about one instrument,
> and the project already knew the other instrument goes the other way.
>
> **Why both are true, and what the actual fix is.** The wide-uncentred fallback
> inflicts TWO changes at once: the window widens (10 -> 50 ppm) and its centre
> moves (-10.02 -> 0). `apply` reverses both, and its WIDTH inference is
> circular -- the width is measured on the subpopulation the wrong window could
> already identify, which `OpenDIAlyzer.cpp:5156-5157` states is biased narrow by
> construction ("pass 1 identifies its cleanest precursors first, so their
> scatter is a lower bound on the run's, not an estimate of it"). On Astral,
> where the gate fails routinely and wide is genuinely right, that narrowing is
> what costs 2,547 identifications.
>
> `doc/16-phase34-plan.md:316` already named the resolution: **centre only, never
> infer the width.** Restoring the centre is not circular in the same way -- a
> mode is far more robust to a truncated sample than a scatter -- and it is the
> half of the fallback's damage that is instrument-independent.
>
> `-mass_calibration_offset_only` already exists for exactly this ("Apply the
> fitted mass offset but do NOT narrow the window from it. Separates the two
> things a passing gate does"). The arms `x8_base_wideC` and `x8_imfix_wideC`,
> running now, are that test: 50 ppm wide, centred on -10.021. If they recover
> most of what the narrow C1 path recovered, centre-only is the fix and needs no
> width inference on either instrument.
>
> The confirmation arm this document was conditioned on DID pass:
> `w7_apply_base.tsv` is byte-identical to `i5_base.tsv` (md5
> `d18b09a743cdc5a8903b0808e2866153`), so `apply` is a strict no-op where the
> gate passes. That was necessary and is nowhere near sufficient.


A run whose fragment mass calibration gate fails extracts at +/-50 ppm with no
centring offset, while holding, in the same run, a measurement of the fragment
mass error taken from its own identifications. Enabling that measurement
recovers **+461 identifications, +59.1%**, on a single-factor change.

The option exists (`-mass_width_from_ids apply`), the measurement is already
taken and logged by default, and the default configuration throws it away.

## Measured

500,000 random targets on the full `S08_diaPASEF.mzpeak`, pinned binary
`8a48387cca4537eb`, affine-corrected library, `-rt_window_pass1 400
-rt_window 300`. The two arms differ in one flag.

    arm                          pass 1        pass 2 window            pass 2 IDs
    i5_imfix   (default)     519 / 282,299   +/-50.000 ppm @  0.000          780
    w7_apply_imfix           519 / 282,299   +/-8.859 ppm @ -8.666         1,241

Pass 1 is **identical in both arms** -- same 519 identifications, same 282,299
peak groups -- because the flag only reaches pass 2. So the entire difference is
the pass-2 extraction window.

At matched decoy budget, which does not depend on how many decoys each arm
generated (`shared/libv2/decoy_budget.py`):

    decoy budget   i5_base   i5_imfix   w7_apply_imfix
               1       244         60              411
               2       364        257              518
               5       770        653            1,034
              10     1,182        824            1,234
              20     1,351        943            1,300
              50     1,460      1,061            1,480
             100     1,592      1,176            1,623
             200     1,736      1,339            1,804
             500     2,083      1,784            2,232
           1,000     2,737      2,418            2,840

The rescue arm is at or above the raw-library arm at nine of ten budgets, and is
far above it at the top of the ranking (411 against 244 at a single decoy).

## The path that was not taken

`OpenDIAlyzer.cpp:5555-5578`, when the gate rejects the fit, has three branches
in order:

    C1  pass 1's own identifications measured the error  -> centre and narrow
    C2  nothing measured                                 -> +/-50 ppm, offset 0

C1 is unreachable by default. `measuredWidthOr_` (`OpenDIAlyzer.cpp:5158-5165`)
returns the fallback for any mode but `apply`, and `-mass_width_from_ids`
defaults to `measure` (`OpenDIAlyzer.cpp:1168`) -- measure it, log it, discard
it. The option's own help says applying is "opt-in until it has been measured on
both instruments".

The discarded measurement was correct. i5_imfix logged:

    fragment window from identifications: 519 accepted precursors, centre -8.67
    ppm, per-fragment sigma 2.95 ppm (p95 6.62), measured through +/-50.00 ppm
    -> suggested half-width 8.86 ppm

Against the value the passing arm fitted independently, -10.021 ppm, that is
1.4 ppm out. The run measured the instrument correctly from its own results and
then extracted 5x wide and 8.7 ppm off-centre.

`ChromatogramExtractor.h:459` already warned what that costs on this instrument
class -- "a diaPASEF frame is a merged stack of TIMS scans, so a wide m/z window
admits vastly more interference" -- and the fallback was raised 15 -> 50 ppm on
an **Astral** sweep, where the gate fails routinely and the spectra are cleaner.
The two instruments want opposite widths, which the header says explicitly. What
was missing is that a failed gate on diaPASEF is not the same situation as a
failed gate on Astral, because on diaPASEF pass 1 usually still identifies
enough precursors to measure the error directly.

## What this does to doc/72

doc/72 reported that the affine 1/K0 correction "costs 414 identifications,
-34.7%". That total effect is real under the default configuration -- and it is
now almost entirely accounted for by this fallback rather than by mobility:
enabling C1 on the same library recovers 461.

The mechanism by which the corrected library reached the fallback at all is a
cross-phase coupling worth stating on its own. The fragment mass probe samples
through a **+/-0.010 window on the LIBRARY's 1/K0** (`MassCalibration.cpp:998-1001`,
`MassCalibration.h:439`). The affine displaces 1/K0 by `0.06916 - 0.04677*IM`,
i.e. +0.0224 at 1/K0 1.0 -- 2.24x that window -- putting **79.06% of the 500,000
target precursors outside the probe**. The mass fit is then made on a decimated
sample, its statistic degrades 0.186 -> 0.252, and the gate rejects it.

So editing a library's mobility column silently starves the *mass* calibration.
Nothing in the flag names or the documentation suggests those two stages are
coupled.

## What the gate itself is not guilty of

An earlier draft of this analysis called the 0.252-against-0.25 rejection a
"knife-edge accident". That was wrong twice. Across 213 gate evaluations in
`shared/libv2/*.log`, **208 passes all sit at or below 0.236 and all 6 failures
at or above 0.252** -- a clean 0.016 separation, 8.6x imfix's overshoot. And the
treatment moved the statistic by 0.066, 35x the overshoot; neither term alone
would have failed (2.275/11.510 = 0.198 passes, 2.899/12.227 = 0.237 passes).

The gate correctly reported that the fit it was handed was bad. The defect is
entirely in what happens next.

## Recommendation (RETRACTED -- see the notice at the top of this file)

~~Make `-mass_width_from_ids apply` the default~~, subject to the confirmation arm
`w7_apply_base` (running): where the gate PASSES, `mass_model_.fitted` is true
and control flow never reaches `measuredWidthOr_`, so `apply` must be a strict
no-op. If that holds, C1 costs nothing where the gate works and recovers 59%
where it does not.

The existing guards are the right ones and should stay: `mass_width_min_groups`
200, the `censored` test that refuses an estimate whose scatter reaches the
window it was measured through, and `mass_width_min_ppm` as a floor.

## Caveats

* One file, one instrument, one precursor draw. The Astral benchmark is the
  natural second test and is cheap: one arm, read the pass-2 window line.
* `w7_apply_imfix` vs `i5_base` (1,241 against 1,194) is a CROSS-FAMILY
  comparison -- different mass regimes -- and does not on its own establish that
  the affine correction is neutral. The 2x2 in `run_immass2x2.sh`, which pins
  the mass path identically in both arms, is what tests that.
* q <= 0.01 counts are low-resolution at this scale: the operating point rests
  on 6-10 decoys, so a single tail decoy moves the count 12-17%. +59.1% is far
  outside that band, and the matched-budget table above is the statistic that
  does not have the failure mode.
* w7's pass-2 mobility gate FAILED where i5_imfix's passed, so the rescue arm
  reached 1,241 with **no runtime 1/K0 correction at all**. That is a result,
  not a configuration difference, but it means the two arms' mobility handling
  also differs and the +461 is a total effect of the flag, not a pure window
  effect.
