# Extraction quality: working plan and running log

Started 2026-08-18 evening. Rule for this log: **no full reruns.** Every
experiment uses the 1,000-precursor apex-centred fixture or an existing dump.

## Where the loss is, as measured today

```
33,749 DIA-NN truth precursors present in our library
22,604 (67.0%)  ODIA scored at least one candidate
11,145 (33.0%)  ODIA formed NO candidate            <- the subject here
11,071 (32.8%)  identified at q <= 0.01
```

Rejected reasons, each measured and each NEGATIVE (doc/37, doc/38):
m/z coverage, RT window (75% of the missing are inside it), fragment mass
calibration (correct to 0.25 ppm), library scale, abundance (a correlate of the
real cause, not the cause).

The proximate cause is **Gate C**: co-elution evidence below the (1-alpha)
quantile of the decoy null rejects them before scoring. But Gate C is reading a
genuinely worse trace, which is the real problem.

## The trace-shape comparison -- the anchor for everything below

Same 1,000 precursors, same fragments, both aligned on DIA-NN's own apex
(`--xic 60` against ODIA `-out_chrom` with the library RT overwritten by
DIA-NN's apex and an identity map):

| config | baseline | peak | prominence | FWHM | abs peak |
|---|---:|---:|---:|---:|---:|
| **DIA-NN** | 0.16 | 1.00 | **6.41x** | **8 s** | 11,547 |
| ODIA 10 ppm, im 0.025 (production) | 0.30 | 0.64 | 2.12x | 24 s | 17,124 |
| ODIA 3.6 ppm, im 0.025 | 0.27 | 0.63 | 2.34x | 12 s | 7,827 |
| **ODIA 3.6 ppm, im 0.050** | 0.27 | 0.71 | **2.64x** | **8 s** | 13,818 |
| ODIA 3.6 ppm, im 0.0125 | 0.28 | 0.51 | 1.80x | 116 s | 3,574 |

**We find the peak; we cannot resolve it.** Two independent defects confirmed:

1. **The fragment window is too wide.** Pass 1 measures a 3.62 ppm half-width
   from 8,103 accepted precursors and then DISCARDS it; production extracts at
   +/-10 ppm. Applying the measured value halves FWHM, 24 s -> 12 s.
2. **The mobility window is too narrow and CLIPS the peak.**
   `precursor_im_window` defaults to 0.025 around the library's predicted 1/K0.
   Our CCS error is ~2.8% (~0.028 at 1/K0 1.0) and the production run measured a
   systematic centre offset of **+0.0194** -- together they put the true mobility
   at or past the window edge. Widening to 0.050 recovers **1.77x** the peak
   intensity and brings FWHM to 8 s, matching DIA-NN exactly. Narrowing to
   0.0125 collapses it (abs peak 3,574, FWHM 116 s), which is the control that
   proves the direction.

**Neither closes the gap.** Prominence 2.64x against 6.41x, because the BASELINE
does not move (0.27-0.30 throughout, against DIA-NN's 0.16). Widening the
mobility window did not raise it either, so the excess baseline is not simply
co-isolated interference admitted by a loose window.

## Ranked worklist

| # | item | status | why it is ranked here |
|---|---|---|---|
| 1 | Apply the measured fragment window instead of discarding it | **measured, not yet wired** | free; halves FWHM; the number already exists in the run |
| 2 | Widen / recentre `precursor_im_window` | **measured** | 1.77x peak intensity, FWHM to DIA-NN's 8 s |
| 2b | Recentre rather than widen: apply the +0.0194 offset | OPEN | widening admits more interference in principle; recentring should not. Needs the mobility calibration to run BEFORE extraction, i.e. a re-extraction after calibration |
| 3 | Explain the 0.27 vs 0.16 baseline | OPEN | the whole residual gap. Candidates: `aggregate=sum` over the mobility box vs max; no background subtraction; genuine interference |
| 4 | Per-fragment apex spread | OPEN | if fragments disagree on apex, the SUM broadens even when each is sharp. Measurable from dumps already on disk |
| 5 | Gate C: feature not hard cut | OPEN | doc/38; couple to FDP (11.76%) before touching |
| 6 | Decoy generation's role | OPEN | mutation decoys keep most of the b/y series, so a decoy re-finds its target's peak -- this is why the RT-seed decoy control failed on DIA-NN's library. Does the same similarity inflate the Gate C null? |
| 7 | Re-extraction after recalibration | OPEN | items 1, 2b and the mobility calibration all only take effect on a SECOND extraction. Quantify what a third pass would buy before building it |

## CORRECTION: the first shape table used a broken statistic

Kimi caught it. Under max-normalisation the apex must read 1.00 by
construction; the table above reads **0.64 for ODIA and 1.00 for DIA-NN**, which
means ODIA's global max is usually NOT at the aligned apex -- 55-61% of ODIA
traces have their global max more than 10 s away (DIA-NN: 37.9%). Every ODIA row
was divided by an off-apex spike. The 116 s "FWHM" was a pure artefact of that.

Re-measured with a robust statistic (apex = mean of the top 3 points within
+/-10 s of the aligned apex; baseline = q10 of the trace; ABSOLUTE units):

| config | apex | baseline | prominence |
|---|---:|---:|---:|
| **DIA-NN** | 6,164 | **914** | **0.868** |
| ODIA 10 ppm, im 0.025 (production) | 9,478 | 3,173 | 0.659 |
| **ODIA 3.6 ppm, im 0.025** | 3,871 | **1,061** | **0.736** |
| ODIA 3.6 ppm, im 0.050 | 7,302 | 2,204 | 0.698 |
| ODIA 3.6 ppm, im 0.050, aggregate=max | 2,164 | 1,074 | 0.524 |

**Three conclusions above are now reversed or corrected:**

1. **The baseline DOES scale with box volume.** Mass 10 -> 3.6 ppm cuts it
   3,173 -> 1,061 (3.0x against a 2.8x window ratio); mobility 0.025 -> 0.050
   raises it 1,061 -> 2,204 (2.1x against a 2x ratio). The "baseline insensitive
   to both filters" claim -- the entire premise of the review brief -- was the
   normalisation artefact.
2. **Widening the mobility window is WRONG.** It grows apex and baseline
   together and prominence FALLS (0.736 -> 0.698). The default 0.025 is right.
   Kimi also found why the widening looked free: there are TWO mobility gates,
   the per-precursor window AND the isolation window's own IM band
   (`ChromatogramExtractor.cpp:1139-1145`), so widening past the band admits
   nothing. And recentring already exists -- the IM centre is
   library 1/K0 + `mobility_model->offsetFor(...)` at index build (`:986-993`),
   so worklist item 2b was already implemented.
3. **The gap to DIA-NN is much smaller than reported.** Prominence 0.736 vs
   0.868; baseline 1,061 vs 914, i.e. **1.16x apart**, not 2.8x.

**`aggregate=max` is worse** (0.524): it collapses the apex 3.4x while only
halving the baseline, so the real signal genuinely is spread over several TIMS
scans and `sum` integrates it correctly. But the baseline SURVIVES max
(1,074 vs 1,061), which is kimi's diagnostic that what remains is a persistent
isobaric, IM-coincident interferent rather than a diffuse noise floor -- i.e.
real co-isolated signal that no box tightening can reject.

## Where this leaves the worklist

| # | item | status |
|---|---|---|
| 1 | **`-mass_width_from_ids apply`** | **measured on S08: prominence 0.659 -> 0.736, baseline 3x lower.** The option exists and is opt-in pending a second instrument (doc/13's rule). Astral is the blocker, not the code |
| 2 | mobility window | **CLOSED -- the default is correct**; widening hurts, recentring already implemented |
| 3 | baseline mechanism | **mostly explained**: box volume plus a persistent isobaric interferent. Residual gap to DIA-NN is 1.16x |
| 4 | per-fragment apex spread | OPEN, measurable from dumps on disk |
| 5 | Gate C as a feature | OPEN, couple to FDP |
| 6 | decoy generation's role | OPEN |
| 7 | re-extraction after recalibration | **now concrete**: the mass width is measured in pass 1 and can only be applied in pass 2, which IS the re-extraction. Item 1 is that question |

## Log

- **2026-08-18 21:xx** items 1-2 measured on the 1,000-precursor fixture.
  First conclusion (`im 0.050` best) was WRONG -- see the correction above.
- **2026-08-18 22:xx** statistic fixed after kimi's review; mobility conclusion
  reversed; `aggregate=max` tested and rejected. Best configuration is
  **`-fragment_ppm 3.6` at the DEFAULT mobility width**, prominence 0.736
  against DIA-NN's 0.868. Still not adopted as a default: needs Astral, and
  needs the entrapment FDP measured alongside, because everything that has
  raised sensitivity today has also raised FDP.

## Iteration 2: the selection-bias control, and where the deficit actually is

Codex attacked the comparison itself: the precursors were chosen BY DIA-NN and
aligned on DIA-NN's OWN apex, estimated from the same signal being plotted, so
its 1.00 at zero is close to built in. Real concern, and testable.

**Measured: alignment bias is empirically nil here.**

| aligned on | tool | apex | baseline | prominence |
|---|---|---:|---:|---:|
| DIA-NN's apex | DIA-NN | 16,783 | 1,144 | 0.937 |
| DIA-NN's apex | ODIA | 6,515 | 1,137 | 0.834 |
| ODIA's apex | DIA-NN | 16,806 | 1,144 | 0.937 |
| ODIA's apex | ODIA | 6,538 | 1,138 | 0.835 |

Swapping the anchor moves nothing, because on precursors both tools identify the
apexes agree: **median difference 0.0 s, p90 |diff| 2.8 s, 96% within 10 s.**
Our retention-time determination is not the problem, and the comparison is fair.

**The deficit is uniform, not concentrated in the hard cases:**

| group | tool | apex | baseline | prominence | gap |
|---|---|---:|---:|---:|---:|
| SHARED (297, both identify) | DIA-NN | 16,783 | 1,144 | 0.937 | |
| | ODIA | 6,515 | 1,137 | 0.834 | **0.103** |
| HARD (703, DIA-NN only) | DIA-NN | 4,595 | 836 | 0.826 | |
| | ODIA | 3,280 | 1,011 | 0.697 | **0.129** |

Three things follow.

1. **ODIA is ~0.10-0.13 prominence below DIA-NN everywhere.** It is not that hard
   precursors are specially bad for us; we are uniformly slightly worse, and on
   bright precursors that costs nothing because both tools clear their
   thresholds anyway.
2. **The hard set is genuinely dim for BOTH tools** -- DIA-NN's own apex falls
   16,783 -> 4,595 (3.7x). These are low-abundance precursors, not ones DIA-NN
   finds easy.
3. **The mechanism is now arithmetic.** DIA-NN identifies the hard set at
   prominence 0.826 -- almost exactly what ODIA achieves on the SHARED set
   (0.834). On those same hard precursors ODIA delivers 0.697. So our constant
   ~0.13 deficit is precisely what pushes a dim precursor below the bar that
   DIA-NN still clears. Closing the deficit, not lowering a threshold, is the
   fix.

On SHARED precursors our baseline EQUALS DIA-NN's (1,137 vs 1,144); on HARD it
is 1.21x worse (1,011 vs 836). So the interferent floor bites hardest exactly
where the signal is weakest, which is consistent with a persistent co-isolated
species rather than a detector floor.

Absolute apex values are NOT comparable between the tools (DIA-NN reads 2.6x
ours on SHARED, 1.4x on HARD); `--xic` scaling is unknown. Only the ratio
statistics are used above.

### Next

- The ~0.13 prominence deficit is now the single number to attack. Candidates in
  order: per-fragment apex spread (a sum over fragments whose apexes disagree is
  flatter than each fragment); the interferent floor on dim precursors; whether
  DIA-NN drops interfered fragments before summing (its log reports "Removing
  interfering precursors").
- `-mass_width_from_ids apply` still needs Astral before it can be a default.
