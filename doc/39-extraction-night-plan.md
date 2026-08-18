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

## Iteration 4: the trace comparison, redone target-only

`-out_chrom` had no Decoy column and a decoy reconstructs its target's
Precursor.Id, so every ODIA trace statistic in iterations 1-3 was computed over
target+decoy against DIA-NN's target-only XICs. **Exactly 50.0% of the rows
(1,015,290 of 2,030,581) were decoy.** Fixed in `2438328`; re-measured:

| group | tool | apex | baseline | prominence | frag spread p75 | frags on apex | n frag |
|---|---|---:|---:|---:|---:|---:|---:|
| SHARED (297) | DIA-NN | 16,783 | 1,144 | 0.937 | 1.4 s | 83.3% | 10 |
| SHARED | **ODIA** | 5,031 | **518** | **0.908** | 9.0 s | 66.7% | 12 |
| HARD (703) | DIA-NN | 4,595 | 836 | 0.826 | 13.2 s | 55.6% | 10 |
| HARD | **ODIA** | 2,008 | **434** | **0.793** | 17.3 s | 33.3% | 12 |

### What this retracts

- **"Our baseline is 1.7-2.8x DIA-NN's."** FALSE, and backwards. Target-only our
  baseline is **518 against DIA-NN's 1,144** on shared precursors and **434
  against 836** on hard ones -- we are 1.9-2.2x CLEANER. Every "interferent
  floor" conclusion built on that number is withdrawn, including the reading of
  the `aggregate=max` result.
- **"A constant ~0.13 prominence deficit."** The real gap is **0.029 on shared
  and 0.033 on hard** -- four times smaller, i.e. we are within ~3% of DIA-NN.
- **"We use 24 fragments to DIA-NN's 10."** We use 12 to its 10. The 24 was
  target+decoy.
- **"41.7% of our fragments peak at the apex against DIA-NN's 83.3%."** The true
  figure is 66.7%.

### What survives, and is now the whole remaining gap

**Per-fragment apex agreement.** 66.7% vs 83.3% on shared, 33.3% vs 55.6% on
hard -- a real 17-22 point deficit, with spread 9.0 s against 1.4 s. Our summed
trace is nearly as good as DIA-NN's (prominence within 3%) while our INDIVIDUAL
fragments scatter far more about the apex. That pattern is per-fragment noise,
not systematic misplacement: a sum over noisy-but-unbiased fragments still peaks
in the right place.

### The consequence for the whole investigation

**Extraction quality is not what loses the 11,145 precursors.** Our traces are
cleaner than DIA-NN's, our prominence is within 3%, our apexes agree with its to
a median of 0.0 s. The loss therefore sits where iteration 0 first put it --
Gate C and the scoring that follows -- and the extraction detour has ruled out
the alternative rather than found the cause.

Standing items 1 (`-mass_width_from_ids apply`, real: 3x lower baseline) and the
mobility default (confirmed correct) survive as improvements on their own terms,
but neither is the explanation for the missing third.

## Iteration 5: Gate C measured directly. It rejects 83.6% of TRUE positives.

With the target/decoy split now trustworthy, `coelutionEvidence` was ported
faithfully (sqrt, robust z per transition against its own MAD, sum, smooth over
2*half+1, take the max) and evaluated on the fixture. Every one of these 1,000
targets is a precursor DIA-NN identifies at q <= 0.01, and every decoy is ODIA's
own generated decoy for it.

```
  targets  n=1000  median statistic 4.94
  decoys   n=1000  median statistic 3.20
  tau = 95th percentile of decoys = 13.24     (alpha = 0.05)

  targets ADMITTED  16.4%          decoys admitted  5.0%
    SHARED (both tools find)  median 9.96   admitted 38.0%
    HARD   (DIA-NN only)      median 4.08   admitted  7.3%
```

**Cross-validated independently.** The alpha sweep on the 1,011,145-precursor
library formed candidates for 86,832 of 529,262 targets = **16.4%**. This
fixture, a different library and a different code path, gives targets >= tau =
**16.4%**. Two independent measurements, same number.

### What is wrong with Gate C

1. **It rejects 83.6% of known-true precursors.** Not marginal ones -- these are
   all DIA-NN identifications at 1% FDR.
2. **Its statistic barely separates the classes.** Target median 4.94 against
   decoy median 3.20: a factor of 1.54, with heavily overlapping distributions.
   A threshold at the decoy 95th percentile therefore sits at 13.24, which is
   2.7x ABOVE the target median. Rejecting most targets is the arithmetic
   consequence, not a tuning accident.
3. **`alpha` is documented as "the false-admit rate, by construction"** and it is
   -- for decoys. Nothing in the design bounds the false-REJECT rate on targets,
   and nobody measured it until now. It is 83.6%.

### What follows

Tuning alpha treats the symptom: at alpha=0.50 the sweep recovered 4,178 IDs but
at 43.7% entrapment FDP. The statistic is what fails, so the options are

- **use the statistic as a FEATURE, not a hard cut** -- let the classifier weigh
  it against the other 14 sub-scores, where a 1.54x separation is worth
  something rather than fatal (this is also what OpenSWATH does: candidate
  formation applies no correlation gate at all, and library correlation enters
  only as a score -- deep-research, verified 3-0 against the OpenMS 3.6 source);
- or **replace the statistic** with one that separates better before any
  thresholding is considered.

Both are scoring-side changes. Extraction has now been excluded twice: our
traces are cleaner than DIA-NN's (iteration 4) and our apexes agree with its to
a median of 0.0 s (iteration 2).

### Caveats on this measurement

- The statistic was re-implemented in Python from the C++; it is a faithful port
  but not the same code.
- The fixture's traces span +/-60 s around DIA-NN's apex, narrower than
  production's window, so absolute values are not production's.
- Production calibrates tau from the first 20,000 decoys across the whole
  library, a more diverse population than these 1,000 matched decoys.

None of these affect the target-versus-decoy comparison, which is the point, and
the independent 16.4% agreement argues the port is sound.

## Iteration 6: decoy generation is NOT the problem; the gate is, and here is the arithmetic

### Decoy generation: hypothesis REFUTED by measurement

The standing suspicion -- DIA-NN-style mutation changes only two residues, so a
decoy keeps most of its target's b/y series and re-extracts its signal, inflating
the null -- is **false for the transitions that are actually extracted**:

```
1,000 target/decoy pairs
  decoy fragments shared with its target: median 0.0%, p90 0.0%
  940 of 1,000 pairs share NO fragment m/z at all
  correlation(shared fraction, decoy statistic) = +0.085

  tau from ALL decoys                   13.24  -> targets admitted 16.4%
  tau from decoys sharing NO fragments   12.63 -> targets admitted 17.2%
```

Mutation does change two residues, but each library entry stores its own top-12
fragments by predicted intensity, and the target's twelve and the decoy's twelve
essentially never coincide. **Switching to shuffle or reverse decoys would move
tau by 5% and target admission by 0.8 points.** Worklist item 6 is closed
negative, and the same reasoning retracts my earlier explanation for why the
RT-seed decoy control failed on DIA-NN's library -- that was also premised on
shared fragments, and it is wrong.

What IS confirmed in code (kimi): the decoy inherits its target's precursor m/z,
iRT, IM, CCS and charge (`LibraryGenerator.cpp:1014-1024`), so it is extracted in
the same isolation window at the same retention time and therefore sees the same
CO-ISOLATED INTERFERENCE. That, not shared fragments, is why the decoy median
sits at 3.20 rather than at a noise floor. The null is elevated by shared
interference, which is arguably correct behaviour for a null.

### The gate's arithmetic, from kimi

Fitting log-normals to the measured medians and the 16.4% admission:

```
  sigma ~ 0.645 (one sigma is a 1.9x fold change)
  d' = 0.43/0.645 = 0.67      AUC = 0.68
  to admit 90% of targets, decoy admission must rise to ~73%
```

**5% decoy admission and 84% target rejection are the same fact.** d' = 0.67
describes a usable FEATURE and an indefensible GATE; there is no operating point
that both bounds compute and retains targets.

### The strongest evidence is tool-independent

Kimi's sharpest point: forget DIA-NN. **The gate admits only 38.0% of the
precursors ODIA ITSELF identifies.** The pipeline rejects 62% of its own output
before scoring it. No cross-tool comparison is needed to see that is wrong.

### Two further code facts

- **The gate was validated on SYNTHETIC data only** (`PeakGroupScorer.cpp:864-871`
  records the Astral measurement it replaced, and the "100% present-precursor
  retention" claim comes from inserted signal). It has never been checked
  against real targets on a real run until tonight.
- **My warm-up premise was wrong.** Precursors reach the gate in RETENTION-TIME
  order, not m/z order -- `ChromatogramExtractor.cpp:307` says so outright. So
  tau is calibrated from the ~20,000 earliest-eluting decoys and applied to late
  eluters, and every target emitted before the null is ready is admitted
  unconditionally, i.e. early-RT precursors are gate-free. That is an RT-dependent
  sensitivity artefact in the final report.
- **Kimi also claims tau is nondeterministic** because `admit` runs under a mutex
  with multi-threaded emission. That CONTRADICTS a measurement of my own: the
  emit loop calling `sink.accept` runs on the driver thread outside
  `pool_impl.run`, which is exactly why the sink was 52.6% of pass 1 (doc/34).
  Recorded as unresolved; I side with serial until someone measures tau twice.

### Decision

Demote the statistic to a sub-score. It is already computed per precursor, so it
is a free 20th feature. If compute must be bounded, rank and cap (top-K) rather
than reject outright -- a ranker degrades gracefully under distribution shift, a
threshold calibrated on the first 20,000 arrivals does not. The documented risk
is the semi-supervised loop failing to ignite at low true-positive rates
(`lda.h:807-813`), which must be measured rather than assumed away.

## Iteration 7: RETRACT "83.6% of true positives" as a production figure

Codex found what kimi and I both missed. Three problems, the third decisive.

**1. The "independent cross-validation" was a coincidence of unequal
populations.** My fixture measured *of DIA-NN-positive precursors, 16.4% exceed
tau*. The C++ sweep measured *of ALL 529,262 library targets -- mostly ballast
that is genuinely absent from the sample -- 16.4% formed candidates*. Different
denominators, different estimands. Equal percentages validate nothing.

**2. Warm-up inflates the C++ number.** With 481,883 decoys and 529,262 targets
interleaved, roughly 20,000 x 529,262/481,883 ~ 22,000 targets are admitted
UNCONDITIONALLY before the null is ready. Steady-state target admission is then
(86,832-21,970)/(529,262-21,970) ~ **12.8%**, not 16.4%.

**3. The fixture does not reproduce deployed decisions.** All 297 SHARED
precursors were identified by the production run, so every one of them PASSED the
real Gate C. Yet only 38.0% clear the fixture's tau. That is a direct
contradiction and it means the fixture's threshold is too high, most likely
because its decoys are paired with KNOWN-PRESENT targets and therefore inherit
parent-signal leakage, while production calibrates from library decoys whose
parents are mostly absent.

### What is retracted

- **"Gate C rejects 83.6% of true positives"** as a statement about production.
  Not established.
- **"The gate admits only 38.0% of the precursors ODIA itself identifies"** --
  kimi's "strongest, tool-independent" point. It is a fixture artefact; those
  precursors demonstrably passed the real gate.
- The claimed independent confirmation between the two 16.4% figures.

### What survives

- **836 of 1,000 reference-positive precursors fail this fixture-calibrated
  threshold** -- true as stated, 95% binomial interval 81.3-85.9%, and not
  explicable by DIA-NN false positives (10 false among the 1,000 would move it
  to ~83.4-84.4%).
- **The statistic separates weakly**: target median 4.94 vs decoy median 3.20,
  d' ~ 0.67, AUC ~ 0.68. This comes from the two distributions, not from tau, so
  the fixture's threshold error does not touch it. A d' of 0.67 is a usable
  feature and a poor gate at ANY operating point.
- **Decoy generation is not the cause** (iteration 6) -- that measurement is
  about fragment overlap and tau shifts within the fixture, and its conclusion
  (0.8 points) is unaffected.
- The 5.0% decoy admission was always tautological -- those same decoys defined
  their own 95th percentile.

### The measurement that would settle it

Per-precursor C++ versus Python concordance, not aggregate percentages:
instrument the real gate to log (precursor, statistic, tau, admitted) and join
against the fixture. That also yields the production false-negative rate
directly, which is the number I claimed and do not have.

**Implementing gate-as-feature is deferred until that exists.** Building on an
unvalidated 83.6% is how the last four retractions happened.
