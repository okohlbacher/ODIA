# Pass-2 extraction support: the semantics, verified from source

Prerequisite that review round 13 named before D0/D3 may be touched. Established
2026-08-15 by reading the code and the reference logs, not by inference.

## 1. `rt_window_seconds` is a HALF-width in the extractor

`src/extract/ChromatogramExtractor.cpp:509-510`

    lo = std::max(lo, axis[w].lowerBound(centre - options.rt_window_seconds));
    hi = std::min(hi, axis[w].lowerBound(centre + options.rt_window_seconds));

Extraction support is `centre ± rt_window_seconds`, i.e. FULL width
`2 x pass2_window`.

## 2. The prefilter treats the same value as a FULL width -- factor-2 disagreement

`src/OpenDIAlyzer.cpp:2045`

    po.rt_half_window = (pw > 0.0 ? pw : pass2_window) * 0.5;

So with the same `pass2_window`, the prefilter admits `± pass2_window/2` while
the extractor admits `± pass2_window`. One of the two is wrong by 2x.

**Latent, not active.** `prefilter` defaults to `"off"`
(`src/OpenDIAlyzer.cpp:473`) and neither reference run enabled it. This has NOT
cost identifications in any measurement on record. It is a bug to fix before the
prefilter is ever switched on, not an explanation for anything measured so far.

## 3. The cap was binding in BOTH reference runs -- the derived window never took effect

`ref_astral.log:96`

    pass 2 extraction window 60 s (2 x p95 47.9551 s, floor 20 s, cap 60 s)

`pass2_window = min(cap, max(floor, factor x p95)) = min(60, 95.9) = 60`.

Both reference runs were invoked with `-rt_window 60 -rt_window_pass1 60`
(verified in `ref_astral.log` and `49.log`). The residual-driven narrowing that
D0 and D3 are built around **did not run**: the cap I passed on the command line
overrode it, and the runs extracted at the flat 60 s.

This is the fifth flag-semantics trap of the session and the same shape as
`-pass1_precursors` being a STRIDE: a measurement attributed to the algorithm
that was in fact produced by a flag I set. See PART E.

## 4. What this does and does not say about "11.7% truth outside the support"

Actual support was `± 60 s`. The anchors' p95 residual was `47.96 s`. So the
window was WIDER than the anchor p95 -- the 11.7% cannot be explained by a
too-narrow window.

If the 11.7% figure was measured on a run with this configuration, the
implication is that **anchor residuals understate library residuals**: a window
covering >95% of anchors covers only ~88% of all precursors. That is selection
bias in the fitted subset -- exactly the failure
`memory: odia-accepted-groups-are-a-biased-sample` records, and exactly what
`doc/14:160-166` forbids sizing a window from.

**The provenance does not exist.** Searched every `.log`, `.py` and `.txt` in
`bench/` for the bucket figures co-occurring. Three files matched on substring
and all three are coincidental -- they are OpenSWATH progress percentages
(`73.76 %`, `11.76 %` in `osw_assay.log`). No saved artefact contains the table,
and no script reproduces it.

So the candidate-flow table

    0.2% no candidate / 11.7% truth outside support /
    12.7% reachable-not-picked / 1.7% ranking / 73.7% correct

exists ONLY in conversation. It is the most-cited result in the plan -- the
entire front-of-pipeline diagnosis rests on it -- and it cannot be reproduced,
checked, or attributed to a command line.

This is not a claim that it is wrong. It is a claim that it is **unverifiable as
it stands**, and therefore cannot carry the weight the plan puts on it. It must
be REGENERATED from a scripted, logged run before any step depends on it. That
regeneration is now a prerequisite of D3, alongside the held-out window
derivation.

## Consequence for the plan

D3 cannot derive the window from anchor percentiles at all, whether the
percentile is 95 or 99, because the anchors are the biased sample. The
derivation needs a HELD-OUT set under the DEPLOYED transform. Rebuilding D3 on
that basis is the next step, and it must precede D0.

---

# The regenerated flow table (2026-08-15)

Script: `bench/flow_table.py`. Artefacts: `bench/flow_table.log`, `.json`.

## The join bug, and that the ORIGINAL analysis had it too

DIA-NN writes cysteine bare (`AACLPLPGYR2`); ODIA writes
`C(Carbamidomethyl)`. Joining on `Precursor.Id` therefore drops EVERY
cysteine-containing precursor -- 2,269 of 12,308, 100% cysteine-containing,
97.8% rescued by rewriting the notation.

The non-cysteine subset numbers **exactly 10,040** -- the denominator the
original table quoted ("of DIA-NN's 10,040"). The original flow table was
computed on cysteine-free peptides only, and nobody knew.

## Result, denominator = 12,308 DIA-NN confident precursors

| bucket | quoted | corrected |
|---|---|---|
| no candidate | 0.2% | **0.41%** |
| truth outside support | 11.7% | **27.54%** |
| reachable, not picked | 12.7% | **13.67%** |
| ranking | 1.7% | **3.17%** |
| correct | 73.7% | **55.22%** |

## Cysteine peptides fail almost completely

| subset | n | no cand | outside | not picked | rank | correct |
|---|---|---|---|---|---|---|
| no cysteine | 10,040 | 0.01% | 20.46% | 10.29% | 2.59% | **66.65%** |
| cysteine | 2,268 | 2.16% | 58.91% | 28.62% | 5.73% | **4.59%** |

**4.59% vs 66.65%.** 18.4% of the truth set is being lost almost entirely, and
the analysis join hid exactly those peptides. Candidates DO get built for them
(only 2.16% have none), but at the wrong retention time -- the signature of
extracting at a wrong m/z. Carbamidomethyl is +57.021464 Da; a fixed
modification dropped or double-applied in the library reader or the precursor
m/z calculation would produce precisely this. **Not yet diagnosed -- this is the
next thing to test, and it is the single largest identified loss in the tool.**

## The global RT offset is NOT the problem -- a correction

Over all 12,258 precursors with candidates, truth-minus-deployed-centre:

    p5 -78.5   p25 +0.0   p50 +19.7   p75 +47.5   p95 +148.7   mean +26.4  sd 64.5

    within +/-60 s as-is                : 72.34%
    within +/-60 s after a -19.7 s shift:  73.18%

Removing the global bias buys **0.8 points**. The earlier reading of "median
+84.6 s, so the axis is shifted" was measured ON THE OUTSIDE-SUPPORT GROUP --
conditioning on the outcome. Selection, not bias. Same error as
`memory: odia-accepted-groups-are-a-biased-sample`.

The distribution is WIDE (sd 64.5 s) and asymmetric (p5 -78.5, p95 +148.7), not
displaced. So per-precursor centre accuracy or a wider support is the lever;
re-centring globally is not.

## Consequence

Reaching 95% coverage on this axis needs roughly +/-150 s, against the 60 s
used. That is the opposite direction from D0/D3, which were both written to
NARROW the window. Neither may proceed on the current evidence.

---

# The carbamidomethyl RT defect (2026-08-15)

## Instrument correction first

`RT_DELTA` is stored as `std::fabs(apex_rt - predicted)`
(`PeakGroupScorer.cpp:1136`, and `:1614` says so). So `centre = RT - var_rt_delta`
is only right when the apex is LATER than predicted. The first pass of this
analysis used that and its offsets were biased downward. Centres below are
solved exactly instead: candidates of one precursor share a centre, each gives
`RT +/- |delta|`, and with >=2 candidates only one value is consistent
(12,230 of 12,258 solved).

## The measurement

Offset = observed truth - deployed extraction centre:

| Cys | n | median | within +/-60 s |
|---|---|---|---|
| 0 | 10,024 | -9.1 s | 76.1% |
| 1 | 1,748 | +42.6 s | 57.0% |
| 2 | 399 | +91.1 s | 30.6% |
| 3 | 53 | +96.2 s | 18.9% |

~+50 s per cysteine, additive. Globally the axis is UNBIASED (median -2.1 s,
sd 63.9), so this is cysteine-specific, not a calibration offset.

## Where it is NOT

Regressing the deployed centre on library RT (fit on non-Cys):
`centre = 11.1658 * libRT + 913.00`. Residuals by Cys count: -0.7, +1.0, -0.2,
-8.8 s. **The runtime maps library RT to extraction centre faithfully.** The lag
is in the library's own predicted RT.

The library is ODIA's own (`doc/20:123`), so the RT came from
`LibraryGenerator::predictRetentionTimes` -> PeptDeep. But the encoder is NOT
miswired: `test/peptdeep_reference.py:47` encodes Carbamidomethyl as C2H3NO and
`test/check_peptdeep_encoding.py:84` tests `PEPC(Carbamidomethyl)TIDEK` against
it. Mod site indexing (`PeptDeepEncoder.cpp:144-163`) is the standard layout --
residue k at k+1, N-term 0, C-term n+1. So this is a MODEL-LEVEL systematic, not
a wiring bug: PeptDeep predicts carbamidomethylated peptides too early.

## Why the fine-tuning does not fix it

`RtRefiner::features()` is built from the STRIPPED sequence
(`RtRefiner.cpp:155`, stripper at `:28-34`), so the modification is invisible to
it. Cysteine COUNT is still available (`AAS` includes `C`) and carbamidomethyl is
fixed here, so count is a perfect proxy -- the refiner COULD learn it.

It does not, because the anchor pool has almost no cysteine peptides:

| set | cysteine-containing |
|---|---|
| all scored precursors | 24.66% |
| DIA-NN confident set | 18.43% |
| **ODIA anchors (q<=0.01)** | **2.95%** |

The refiner ran and converged (`ref_astral.log:88-95`, 4,615 then 4,735 anchors,
median residual improving 0.058 s) on a pool that is 97% cysteine-free.

**The loop:** predicted ~50 s early -> outside the +/-60 s window -> not
identified -> not an anchor -> no correction learned -> predicted early.
Fine-tuning cannot break this: the evidence it needs is excluded by the window
it would correct.

## Fixes, in order

1. **Break the loop** -- widen pass 1 enough to admit cysteine peptides
   (+96 s at 3 Cys, so +/-150 s), let them become anchors, then let the refiner
   learn the coefficient. This is the unblocker and needs no new code.
2. **Make the modification explicit in the refiner** -- stop stripping; add a
   modification-composition term. The fixed-mod proxy is luck; variable mods
   (oxidation) are invisible today and would have the same failure with no proxy
   to save them.
3. Seed a per-carbamidomethyl offset analytically as a prior.

This is further evidence AGAINST D0/D3, which both narrow the window.

---

# The refiner is trained on self-confirming false anchors (2026-08-15)

## Hypotheses tested and REFUTED on the way

- "Cysteine anchors are selected to have ~zero offset." REFUTED: cysteine
  anchors show median +35.4 / mean +47.4 s against truth. The signal IS present.
- "Too few cysteine anchors (2.95%)." Codex: ~136 anchors suffice to detect a
  50 s effect. Arithmetic agrees -- standardised, ridge 1.0 against a column
  sum-of-squares in the thousands, clamp ~88 s -- the refiner SHOULD learn ~+47 s.
- "MAD trimming removes them." REFUTED: the run logs `0 trimmed`.

## What is actually happening

| | cysteine anchors | non-cysteine anchors |
|---|---|---|
| n | 81 | 3,378 |
| picked RT inside DIA-NN peak bounds | **23.5%** | 97.5% |
| picked - truth | **-60.2 s** | -0.0 s |
| picked - centre (THE REFINER'S TARGET) | **-17.0 s** | -7.4 s |

**76.5% of cysteine identifications at q<=0.01 are FALSE** -- peaks picked ~60 s
early, at the wrong prediction. They pass FDR because the decoys share the same
wrong prediction. They become anchors. The refiner fits
`observed_rt - calibrated_irt` where `observed_rt` is ODIA'S OWN PICK, so those
rows contribute residuals near zero.

**The refiner is not starved of signal. It is taught the opposite of the truth.**

    wrong prediction -> false peak AT the wrong RT -> anchor says "prediction was
    right" -> no correction learned -> wrong prediction persists

A self-confirming loop. Every internal metric improves under it: anchor count
rises, anchor residual SD falls, the refiner's held-out gate passes. Only an
EXTERNAL truth or an entrapment measure can see it.

## Consequences

1. The reported 3,625 IDs contain a contaminated cysteine subset. Precursor-level
   FDR does not catch it.
2. Any refiner validated on anchor residuals is validating the contamination.
   **Validation must be identifications at 1% FDR, not anchors, not residual SD.**
3. This generalises past cysteine: ANY systematic library RT error self-confirms
   the same way. Cysteine is merely where it is large enough to see.

---

# Review round 14 (codex + kimi) -- what survived, what did not

## CONFIRMED, and strengthened

**The loop, measured at the coefficient level.** Kimi replicated `RtRefiner::fit`
on the real anchors (`q2_refiner_coef.py`):

    as ODIA runs it                  Cys coef +0.69 iRT u ~ +7.5 s
    same machinery, truth anchors             +4.62 u    ~ +50.2 s
    truth injected for the 98 Cys anchors     +3.09 u    ~ +33.6 s

The refiner learns ~15% of the needed correction. This is a measurement, not the
inference I had made from the affine residual.

**E4 reproduces digit-for-digit** (23.5% vs 97.5% inside bounds; -60.19 vs -0.00
s). M/W/P/H residue controls all ~0, so it is cysteine-specific.

## REFUTED -- my errors, do not cite these again

1. **Scale.** `libRT = 175.064*out - 41.089` matches NO model on disk. The true
   relation is `152.2356*out - 39.2322`, r=1.0000, resid_sd=0, from the **stock
   OpenMS model** (`05_odia_lib.sh` passes no `-rt_model`) -- not either
   rtfinetune model. My value came from regressing a sample.
2. **Control coefficients.** My OLS controls (length -2.32, hydrophobic +80.36)
   do not reproduce under any definition kimi tried (it gets +7.5 and ~0). The
   matched-pair result and the dose-response DO survive; the controls do not.
3. **The mechanism.** "The model predicts the unmodified peptide" is
   SIGN-CONTRADICTED: free Cys is MORE hydrophobic than CAM-Cys, so predicting
   the unmodified form would make the library LATE, not early. Kimi verified the
   stored library RTs are reproduced EXACTLY by feeding CAM-encoded sequences to
   the stock model -- **CAM was fed to the predictor**. This is "the model never
   learned a CAM contrast", not "the encoding was never fed".
   Cross-check: DIA-NN's own report errs -5.0/-9.1/-11.3 iRT per Cys -- same
   magnitude, OPPOSITE sign. **The per-Cys deficit is measured; its CAUSE is not
   established.** Act on the deficit, not on the chemistry story.
4. **Oxidation contrast** -0.122 u / -239 s was ~40% overstated; it is -0.088 u
   / ~-156 s.

## The directive's premise does not hold -- measured

Feeding every candidate cannot work until the window moves:

- pass-1 support is +/-60 s; Cys truth sits at **+67.5 s median (p90 121 s)**
  from centre -- outside the window for the majority.
- only **20.4%** of false-picked Cys anchors have ANY candidate within 10 s of
  truth.
- where a truth candidate DOES exist (20/98), **ODIA picked it 100% of the
  time**.

**The picker is innocent; the window censors truth.** "The true peak is still
among the candidates" holds ~1 time in 5, exactly for the precursors that need
it. Estimator work (hard-EM with per-precursor weights is the right one) is
premature until the axis moves, and starting it from a zero correction
reinforces the same loop one level down.

## Phase 0 gate -- PASSED

Offset `RT += 4.3114 iRT * nCys` (fitted fold 0 +48.14 s/Cys, held-out fold 1
+52.53 s/Cys):

    n_Cys coefficient  +49.43 s (t=+43.5)  ->  +2.33 s (t=+2.0)
    1 Cys residual     +49.4 s -> +2.3 s
    2 Cys residual     +96.8 s -> +2.6 s
    3 Cys residual     +97.2 s -> -44.1 s   OVER-CORRECTED

**The effect saturates**: 3 Cys costs +97 s, not 3x49. The additive model
over-shoots on 61 of 2,261 Cys precursors. Both reviewers warned linear
additivity past two cysteines was unsafe. A capped or per-count correction is
the follow-up.

## Ordering, settled by both reviewers

Library/axis first, refiner second. Refiner-first is futile because there is
nothing truthful in its training candidates to learn from. After the offset,
re-run pass 1 -> anchors -> REFIT the map (a stale axis bakes in), and the
refiner's Cys coefficient should now fit ~0. If it relearns +4-5 u, the offset
did not reach the axis the refiner sees.

---

# IH1: the second instrument (2026-08-15)

## There is NO DIA-NN truth for IH1

`diann_ih1_report.parquet` has **0 rows**. `diann_ih1_report.log.txt`:

    ERROR: cannot load the file, skipping
    0 MS1 and 0 MS2 scans in 0 (inferred) and 0 (encoded) cycles
    ERROR: DIA-NN tried but failed to load .../IH1_diaPASEF.mzML
    0 precursors saved

DIA-NN never ran on IH1. **Any IH1-vs-DIA-NN comparison in the project's
reference numbers must be audited** -- see `memory: odia-reference-numbers`.
The cysteine defect CANNOT be validated on IH1 against external truth.

## IH1 reference run

    identified 8514 precursors at 1% FDR
    q <= 0.010   8131 target + 104 entrapment   FDP 7.422%
    q <= 0.050  10449 target + 213 entrapment   FDP 11.829%
    pass 2 extraction window 60 s (2 x p95 42.5019 s, floor 20, cap 60)

More than twice Astral's identifications (8,514 vs 3,625) at SEVEN times the
nominal FDP (7.42% vs Astral's 1.285% measured here). The window was capped on
BOTH instruments -- doc/27 section 3 holds for IH1 too.

## Truth-free cross-instrument test: cysteine share of identifications

The library composition is known, so depletion needs no external truth:

| | IDs@1% | Cys share | available | depletion |
|---|---|---|---|---|
| Astral | 3,659 | 3.01% | 24.59% | **8.2x** |
| IH1 | 8,587 | 8.58% | 20.42% | **2.4x** |

**Present on both instruments** -- the defect is not Astral-specific. But 3.4x
weaker on IH1, which the iRT-domain mechanism PREDICTS: the deficit is a fixed
library property of ~4.6 iRT units, and the damage depends on the run's
seconds-per-iRT slope. A shallower slope turns the same iRT error into fewer
seconds and pushes fewer truths outside +/-60 s. IH1's ion mobility also
suppresses the interference that otherwise fills the wrong window.

Prediction: the iRT-unit correction helps both runs, and helps Astral more.
This is falsifiable and should be checked once the Astral rerun lands.

---

# The CAM rerun: correction validated, pass 2 collapses (2026-08-15, night)

## Result

| | pass 1 | pass 2 |
|---|---|---|
| baseline (`ref_astral.log:71,140`) | 554 | 3,625 |
| CAM-corrected (`52.log`) | **3,109** | **1,759** |

Verified from log context that these ARE pass 1 and pass 2 (13.3M and 13.5M peak
groups respectively), not a probe stage.

**Pass 1 improved 5.6x.** Pass 1 uses the SUPPLIED map directly on the library,
so this is a clean measurement of the correction on identifications -- not on
agreement with DIA-NN. The +50 s/Cys deficit is confirmed by IDs.

**Pass 2 lost more than pass 1 gained**, so the run FAILS kimi's Phase-1 stop
condition (IDs must not drop >5%). The CAM library is a validated DIAGNOSTIC,
NOT an accepted production change.

## The refiner is exonerated by measurement

`cam_astral_refiner.txt` (25 folded weights, `AAS="ACDEFGHIKLMNPQRSTVWY"`, so
index 1 is cysteine):

    Cys -3.960 s   length +0.714   charge +1.377   irt -0.01429   bias +17.300

**Cys coefficient -3.96 s, against +0.69 iRT (~+7.5 s) on the uncorrected
library.** This is kimi's falsifier PASSING: the offset reached the axis the
refiner sees, and it no longer tries to relearn the correction.

Correction magnitude over 40,000 precursors:

    p5 -36.0   p50 -4.3   p95 +17.5   mean -6.1   sd 16.9
    |corr|>60 s: 0.7%     at clamp: 0.0%
    by Cys: 0 -> -3.0 s, 1 -> -6.4 s, 2 -> -10.3 s

Two rounds give ~-12 s mean. **Far too small to explain a 50% ID loss.**

## Where pass 2's gain actually comes from -- and what that implicates

Pass 1 and pass 2 score almost the SAME number of peak groups (13.32M vs 13.55M)
but pass 2 gets 6.5x more IDs on the baseline. The only material difference is
that pass 2's `var_rt_delta` is computed against the FITTED map rather than the
supplied one. So pass 2's advantage IS the fitted map.

In the CAM run the fitted map underperformed the supplied map (3,109 -> 1,759).
With the refiner ruled out by magnitude, **the map fit is the prime suspect**.

## Structural hazard found while reading (not yet shown causal)

`RtRefiner::apply` (`RtRefiner.cpp:362-382`) mutates the library IN PLACE
(`p.irt[i] += corr`) and is called once per refine round
(`OpenDIAlyzer.cpp:1760`), while `features()` takes `p.irt[i]` AS AN INPUT
FEATURE. Corrections therefore accumulate across rounds and each round's input
is the previous round's output. Benign at the measured magnitudes here, but it
is an unguarded feedback path.

## Running

Arm A: CAM library, `-rt_refine off`  -- does the refiner cause the collapse?
Arm B: baseline library, `-rt_refine off` -- control for the same ablation.

---

# ABLATION: the refiner is the cause (2026-08-16)

| configuration | pass 1 | pass 2 |
|---|---|---|
| baseline library, refiner ON | 554 | 3,625 |
| CAM library, refiner ON | 3,109 | **1,759** |
| **CAM library, refiner OFF (Arm A)** | 3,109 | **3,764** |

**H2 REFUTED -- my own hypothesis.** I argued the refiner could not be causal
because its corrections were small (mean -6.1 s, sd 16.9, 0% at clamp). Wrong:
turning it off recovers 2,005 identifications. The refiner is harmless on the
UNCORRECTED library and destructive on the corrected one. Magnitude reasoning
from a reconstructed model was not a substitute for the ablation.

## Arm A full validation vs baseline

| metric | baseline | Arm A |
|---|---|---|
| IDs @1% | 3,659 | **3,799** (+3.8%) |
| entrapment FDP | 1.285% | 1.369% |
| precision vs truth | 95.75% | **96.31%** |
| recovery of DIA-NN | 28.11% | **29.09%** |
| cysteine IDs | 81 | **50** |
| cysteine precision | 23.46% | **28.0%** |
| cysteine recovery | 3.57% | **2.2%** |

**Strictly a REJECT under the stated rule** (FDP must not rise; 1.285 -> 1.369).
The rise is 47 -> 52 entrapment hits, inside Poisson noise at n~50, so it reads
as flat -- but the rule says reject and that is recorded rather than relaxed.

## The correction did NOT do what it was for

Cysteine recovery FELL (3.57% -> 2.2%) and cysteine IDs fell 81 -> 50 while
cysteine precision ROSE (23.46% -> 28.0%). So the fix removed FALSE cysteine
identifications without recovering TRUE ones. Kimi's Phase-1 expectation
("Cys share moves 3.0% toward 15-18%") FAILED: it went to 1.3%.

The +140 net identifications are almost entirely NON-cysteine (3,578 -> 3,749),
consistent with false cysteine anchors having poisoned the map fit and the
classifier for every other peptide.

## Why: the second, independent defect

After correction the cysteine RT residual is ~+2 s, so the window is no longer
the binding constraint. But kimi measured `var_library_corr` **0.149 vs 0.743**
at the true peak for in-window Cys1, with `xcorr_shape` near normal. **The
retention time is now right and the fragments still do not match.** The library
m/z are verified correct (100% CAM-inclusive, 3,001 precursors / 16,750
Cys-spanning fragments), and predicted intensities correlate as well as non-Cys
(r 0.778 vs 0.802). So this is an extraction-level spectral failure that is
NOT explained by RT, mass, or intensity prediction. It is the next target.

## Standing questions

- WHY does the refiner destroy a corrected axis but not an uncorrected one?
  (double correction / in-place cumulative application / gauge freedom between
  map and refiner -- all named by codex, none yet demonstrated)
- Is `-rt_refine off` acceptable as an interim default, given it now BEATS the
  refiner-on baseline on IDs, precision and recovery?

---

# The full 2x2 (2026-08-16) -- the refiner INTERACTS with the library fix

| library | refiner | pass 1 | pass 2 | IDs@1% | FDP | precision | recovery | cys prec | cys rec |
|---|---|---|---|---|---|---|---|---|---|
| baseline | ON  | 554   | 3,625 | 3,659 | 1.285% | 95.75% | 28.11% | 23.46% | 3.57% |
| baseline | OFF | 547   | 2,093 | **2,112** | 1.373% | 96.35% | 16.26% | 26.53% | 2.16% |
| CAM      | ON  | 3,109 | 1,759 | **1,775** | 1.183% | 96.28% | 13.76% | 23.08% | 1.15% |
| CAM      | OFF | 3,109 | 3,764 | **3,799** | 1.369% | 96.31% | 29.09% | 28.00% | 2.20% |

**The refiner is +1,547 IDs on the uncorrected library and -2,024 on the
corrected one.** An INTERACTION, not a main effect. Every single-arm reading
tonight -- including mine -- was misleading because of it.

This is the double-correction signature codex named: the refiner had been
learning to compensate for the library's systematic error, and once that error
is fixed at source its compensation becomes over-correction. It is NOT the
cysteine term (-3.96 s, small). The refiner is harmful precisely when pass 1
hands it GOOD anchors (3,109 vs 547), which points at the global bias (+17.3 s)
and length (+0.714/residue) terms overfitting a better-conditioned anchor set.

**Best configuration: CAM library + `-rt_refine off` -- 3,799 IDs (+3.8%),
precision 96.31% (+0.56), recovery 29.09% (+0.98), FDP 1.369% vs 1.285%.**

The FDP difference is 47 -> 52 entrapment hits, inside Poisson noise at n~50.

## What this does NOT fix

Cysteine recovery is still 2.20% (baseline 3.57%). The RT correction removed
FALSE cysteine IDs without recovering TRUE ones, because of the independent
spectral defect: at the true peak, cysteine `library_corr` is **0.014** vs 0.713
for non-cysteine, with `xcorr_shape` NORMAL (0.650 vs 0.703), and
`usable_fragments` 7.0 vs 10.0 against a library census of 12 fragments
(6 spanning a cysteine, 6 not). The CAM-bearing fragments contribute nothing.

**Retention time and spectra are independent defects, and the spectral one is
larger.** Root cause not established -- five mechanisms proposed tonight, four
already refuted by measurement.

---

# ROOT CAUSE: the two engines searched different chemistry (2026-08-16)

`diann_lib.parquet` (DIA-NN's own search library, 112,783 rows), cysteine
precursors:

| | CAM-inclusive | CAM-FREE |
|---|---|---|
| precursor m/z | 0 (0.00%) | **21,355 (100.00%)** |
| Cys-spanning fragment m/z | 0 (0.00%) | **12,399 (100.00%)** |

e.g. `AACAQLNDFLQEYGTQGCQV` y7 = 692.3032, not 749.3247.
DIA-NN's `Modified.Sequence` equals `Stripped.Sequence` for these peptides --
no modification annotated at all.

**DIA-NN searched this data WITHOUT carbamidomethylation.** ODIA's library is
100% CAM-INCLUSIVE (3,001 precursors / 16,750 Cys-spanning fragments, zero
exceptions), because `LibraryGenerator.h:52` defaults to
`fixed_modifications{"Carbamidomethyl (C)"}`.

## This explains everything at once

1. **The spectral failure.** We extract Cys-spanning fragments +57.021 Da per
   cysteine from where there is no signal. `library_corr` 0.014 vs 0.713;
   `usable_fragments` 7.0 vs 10.0, matching the ~6 NON-spanning fragments of a
   12-fragment assay.
2. **The RT deficit AND ITS SIGN.** Kimi correctly objected that free cysteine
   is MORE hydrophobic than CAM-cysteine. Exactly so: our library predicts the
   CAM form (less hydrophobic, EARLIER) while the real peptide is free (more
   hydrophobic, LATER). Library early by ~50 s. The sign that refuted the
   previous story CONFIRMS this one.
3. **Why candidates exist but never match.** DIA isolation windows are wide, so
   a precursor 57.021/z Da off still falls in a window: candidates are built,
   fragments find nothing, the picker takes interference.
4. **Why the RT offset helped without recovering cysteine IDs.** It corrected
   the retention time of peptides whose FRAGMENTS were still being extracted
   from empty m/z. Symptom, not cause.

## Consequences

- The +4.3114 iRT/Cys library offset is symptom treatment. The mass fix
  supersedes it; whether to keep it as well is an open question (the CAM/free
  retention difference is roughly constant, so it partly stands in for the real
  correction).
- **The DIA-NN benchmark has been apples-to-oranges for ~18-25% of the library
  from the beginning.** `memory: odia-reference-numbers` needs auditing on this
  point, alongside the IH1 finding that DIA-NN never ran there at all.
- The real fix is to generate ODIA's library to match the sample's actual
  alkylation state. That state is not yet independently established -- DIA-NN's
  configuration is strong evidence but is not the wet-lab record.

---

# Review round 15 (codex + kimi) on the pass-2 collapse

Both reviewers independently reached the same verdicts. Codex's runner failed a
third time (`codex-code-mode-host: No such file or directory`) so its review is
of the pasted sources only -- the paste DID work (it quoted the loop code); only
command execution failed.

## Confirmed

- **H1 (map fit) REFUTED.** The maps are equivalent: `52.log` (refiner on) and
  Arm A (off) both print `fitted the retention-time map from 4795 anchors;
  p95 52.4203 s`. `-rt_refine off` removes the LOOP, not the map.
- **H2 (refiner exonerated) REFUTED** -- my magnitude numbers were right
  (kimi recomputed over all 4.97M targets: mean -3.09 s, 0.035% at clamp) but
  the exoneration does not follow.
- **H3 (in-place accumulation) REFUTED.** The loop resets
  `library.precursors().irt = original_irt` before every `apply` (`:1759`),
  builds a fresh `RtRefiner` each round (`:1726`), and rolls back a worse round
  (`:1884-1890`). Nothing accumulates. My "structural hazard" was not real.

## The actual damage mechanism (kimi, measured)

For the 2,320 precursors identified in Arm A and lost with the refiner on:

    apex still covered      |var_rt_delta| p50 IMPROVED 23.6 -> 19.1 s
    not Cys-enriched        2.2% Cys vs 1.9% among Arm A IDs
    shift uncorrelated      r = 0.105 with the model's own correction
    q-ladder                1,970 of them on ONE rung at exactly q=0.01248
    top decoy DScore        16.4 (Arm A) -> 23.5 (refiner on)
    decoy FEATURES unmoved  decoy |var_rt_delta| p50 25.6 vs 25.2 s

**Small corrections flip the semi-supervised GBT into a different optimum; the
decoy tail gets heavier and the target cutoff moves 12.66 -> 17.53.** The axis
is fine; the ladder is not. The loop's acceptance guard (anchor p50 on a frozen,
self-referential set) is structurally blind to this -- the documented failure
mode firing in production.

## Numbers that do NOT reconcile -- do not quote effect sizes precisely

| | tool log | validate_run.py |
|---|---|---|
| CAM + refiner off | 3,764 | 3,799 |
| baseline + refiner off | 2,093 | 2,112 |

And baseline pass 1 differs across two runs that should be identical:
**554 vs 547**, though `-rt_refine` should not touch pass 1. So the pipeline is
NOT bit-reproducible and differences of tens of precursors are noise.
Reconcile count semantics (my validator canonicalises CAM and dedups by
best-DScore) before treating any small effect as real.

## Pass-1 counts are NOT 1% FDR

Pass-1 entrapment FDP is **6.2%-8.6%** at nominal q<=0.01 -- the ladder is 6-8x
anti-conservative there. **"Pass 1 improved 5.6x" must not be quoted as a
1%-FDR result.** Retracted as stated.

## The isolation-window rescue is impossible

mzpeak isolation offsets are +/-2.0 Th (4 Th full). A precursor 57.021/z Da off
is 28.5 Th away at z=2 -- it CANNOT be co-isolated. So the 475/623 cysteine
candidates at the "true peak" are co-isolated interferents plus selection
(their apex was conditioned to fall in DIA-NN's window): real traces of OTHER
peptides, hence normal `xcorr_shape`, ~7 usable (the non-spanning) fragments,
and `library_corr` ~ 0. **"Shape is near-normal, so the peak is real" does not
follow** -- my inference, withdrawn.

## Code defect to land regardless of any arm

The refine loop's acceptance is anchor-p50 convergence. It accepted a round that
HALVED pass-2 identifications while reporting improved residuals. Replace with
the `validate_run.py` rule: IDs at 1% FDR rise AND entrapment FDP does not rise
AND truth-precision does not fall. Anchor statistics may be logged, never used
to accept.

---

# Implemented: `-fixed_modifications` (2026-08-16)

The root cause was unreachable from the command line: `LibraryGenerator.h:52`
hard-coded `{"Carbamidomethyl (C)"}` and no option overrode it, so there was no
way to search this dataset correctly.

Changes in `src/OpenDIAlyzer.cpp`:

1. **`-fixed_modifications`** registered (comma-separated OpenMS UniMod names,
   default `Carbamidomethyl (C)`, empty string = none), parsed into
   `params.fixed_modifications`, and logged as
   `fixed modifications: ...` / `NONE`.
2. **Added to the library cache fingerprint** as `;fixmod=`. It was absent, so
   flipping the alkylation would have silently reused a library built with the
   other one -- kimi's stale-cache hazard, confirmed in the source.

## Two self-inflicted failures worth recording

**a) A bare `cmake --build` destroyed the binary.** Without
`source scripts/env.sh` the link fails on transitive OpenMS/mzpeak dependencies,
and the failed link DELETES the target first. A 2.5 h benchmark was running at
the time and survived only because Linux keeps a running process's inode alive;
recovered via `cp /proc/<pid>/exe`. Now in
`memory: odia-build-needs-env-sh`.

**b) The first version of the option broke every run.**
`registerStringOption_`'s `required` parameter defaults to TRUE, and TOPPBase
forbids a required option with a non-empty default -- so the tool threw
`InvalidValue` at startup before doing anything. 17 of 67 ctest cases failed.
Passing `false` fixes it; **67/67 now pass**. The lesson is the ordinary one:
run the test suite after touching argument registration, not just the compiler.

---

# ARM C: the mass fix works, and exposes the real defect (2026-08-16)

No-CAM library (`odia_lib_nocam.parquet`), `-rt_refine off`.

| metric | baseline | Arm A (CAM offset) | **Arm C (no-CAM)** |
|---|---|---|---|
| **cysteine precision** | 23.46% | 28.0% | **97.04%** |
| **cysteine recovery** | 3.57% | 2.2% | **8.95%** |
| cysteine IDs | 81 | 50 | **203** |
| precision vs truth | 95.75% | 96.31% | **98.07%** |
| entrapment FDP | 1.285% | 1.369% | **1.241%** |
| total IDs @1% | 3,659 | 3,799 | **1,289** |
| recovery of DIA-NN | 28.11% | 29.09% | 10.09% |

**The root cause is confirmed beyond argument: cysteine identifications went
from 23% correct to 97% correct.** Precision and FDP both improved globally.

## But total identifications collapsed, and it is NOT the strip

The strip touched only cysteine rows, yet non-cysteine IDs fell ~3,578 -> ~1,086
(-70%). Verified directly: **41,157 non-cysteine rows byte-identical between the
two libraries, 0 changed; every cysteine row correctly shifted.**

So this is the CLASSIFIER/Q-LADDER re-optimising -- the same amplifier kimi
identified behind the refiner collapse ("the amplifier is the classifier+q-ladder,
not the map"), where the retrained GBT moved the top decoy DScore 16.4 -> 23.5
and the target cutoff 12.66 -> 17.53.

## The real defect is now clear

**ODIA's semi-supervised classifier is unstable to changes in the input
distribution.** Three separate interventions -- a per-cysteine RT offset, the
refiner, and a mass correction -- each swung identifications by 2-3x in ways
uncorrelated with the quality of the change:

| intervention | effect on IDs | effect on precision |
|---|---|---|
| refiner ON (uncorrected lib) | +1,547 | -0.6 pt |
| refiner ON (corrected lib) | -2,024 | ~0 |
| mass fix (Arm C) | -2,510 | **+2.3 pt** |

Arm C is the CORRECT library and yields the FEWEST but CLEANEST
identifications. Estimated true positives: baseline ~3,503, Arm A ~3,659,
Arm C ~1,264. So Arm C is genuinely worse in yield while being right in
chemistry -- the ladder is throwing away real peptides.

**Next target: the classifier and q-value ladder, not the RT axis.** Everything
upstream is now measurably closer to correct than it has ever been.

---

# ARM E: the cysteine defect is FIXED (2026-08-16)

no-CAM masses **and** the +4.3114 iRT/Cys offset, refiner off.

| metric | baseline | ArmA (RT only) | ArmC (mass only) | **ArmE (both)** |
|---|---|---|---|---|
| IDs@1% | 3,625 | 3,764 | 1,278 | 2,676 |
| entrapment FDP | 1.159% | 1.249% | 1.017% | **0.859%** |
| precision vs truth | 95.75% | 96.31% | 98.07% | 97.17% |
| recovery of DIA-NN | 28.1% | 29.09% | 10.08% | 20.95% |
| cysteine IDs | 81 | 50 | 203 | **543** |
| cysteine precision | 23.46% | 28.0% | 97.04% | **95.21%** |
| **cysteine recovery** | 3.57% | 2.2% | 8.95% | **23.94%** |

**Cysteine recovery 3.57% -> 23.94%, a 6.7x improvement, and now ABOVE the
overall recovery rate (20.95%).** The stratum is no longer broken -- it is
recovered slightly better than average.

**It required BOTH halves.** Masses alone (ArmC) gave 8.95%; RT alone (ArmA)
gave 2.2%. The two defects were independent, which is why every single-variable
arm looked like a failure.

ArmE also has the lowest entrapment FDP of any arm (0.859%) -- the only one
comfortably below nominal.

## Per-stratum q-values were a WORKAROUND, not an improvement

| arm | pooled | per-stratum | gain |
|---|---|---|---|
| ArmC | 1,278 | 1,952 | +674 |
| **ArmE** | 2,676 | 2,640 | **-36** |
| baseline | 3,625 | 3,637 | +12 |

In ArmC the cysteine decoys were enriched at high scores (masses fixed, RT not),
so the strata had genuinely different nulls and separating them helped. In ArmE
the strata are exchangeable again and stratifying costs a little. **Do not adopt
per-stratum q; it was compensating for a defect ArmE removes at source.**

## Still open

Total identifications remain below baseline (2,676 vs 3,625). Estimated true
positives: baseline ~3,471, ArmA ~3,625, ArmE ~2,600. The correct-chemistry arm
still finds fewer real peptides, and the cause is now downstream of the library.
