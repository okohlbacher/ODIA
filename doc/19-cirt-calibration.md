# CiRT-seeded calibration and fine-tuning

Status: PLAN + first implementation, 2026-08-14. Supersedes `doc/18` §5's
bootstrap; `doc/18` §9 carries the measurements this rests on.

## 1. Why this and not the four things that failed

One afternoon, four unsupervised anchoring mechanisms each failed against its own
control:

| mechanism | result | control |
|---|---|---|
| unique library m/z match | r = 0.045 | shifted library gave **more** "unique" matches |
| mass-offset mode | prominence 1.0–1.1× | identical |
| whole-library Hough, MS1 | 2.3× | 2.3× |
| dscore-ranked Hough, MS2 | 4.6× | 4.0×, slope wandering 4.2–10.7 |

Then **17 CiRT peptides fitted `RT = 907.0 + 10.873·iRT` at r = 0.9913**, median
residual 32.3 s.

The reason the unsupervised routes cannot work is countable and is not going to
change: 5M library precursors over ~380–980 Th **saturate the m/z axis**, so
matching on precursor m/z carries no identifying information at any tolerance or
map density. CiRT works because it substitutes *external knowledge* — a short list
of peptides known a priori to be present, abundant and conserved — for
information the data cannot supply.

**And the enabling property is smallness.** 149 CiRT precursors can be searched
**blind over the entire gradient**, which is precisely what is unaffordable for
10M. The circularity that has blocked Phase 3 all along — calibration needs
identifications, identifications need calibration — dissolves because a
149-precursor blind search needs no calibration to be affordable.

## 2. Coverage, measured

| | |
|---|---|
| CiRT list (Parker et al., MCP 2015) | 73 sequences |
| present in our predicted library | **47 (64%)** |
| precursors (sequence × charge) | **149** |
| identified by DIA-NN at q ≤ 0.01 | **27** |

No Biognosys iRT spike-ins exist in this sample (0/11), so the endogenous set is
the only option — which is also the point: it needs no experimental change.

27 is the realistic ceiling on Astral and is ample; the r = 0.9913 fit used 17.

## 3. The pipeline

```
  1. CiRT BLIND SEARCH     149 precursors, full RT range, no calibration.
                           Cheap because it is 149 and not 10^7.
  2. ROBUST LINE FIT       RANSAC / Theil-Sen on (library iRT, observed apex RT).
                           Must tolerate ~50% wrong apexes -- a blind search on
                           an uncalibrated run will place some of them on
                           interference. This is the regime robust fitting is
                           FOR, unlike the 98%-contaminated whole-library case
                           where it demonstrably failed.
  3. PREDICT LIBRARY RT    Apply the line to all precursors.
  4. EXTRACT IN A WINDOW   +/-60 s rather than blind. The MS2 population is now
                           mostly correctly placed, which is the precondition
                           everything downstream needed.
  5. MS1 CONFIRMATION      Measured to lift r 0.164 -> 0.420 with its control at
                           0.211 -- real, and previously insufficient only
                           because it was consolidating a ~98%-wrong population.
                           Given step 4 it consolidates a mostly-right one.
  6. ROBUST REFIT          On the union of MS1-confirmed and MS2 anchors.
  7. FINE-TUNING           Per-run refinement in iRT space, the frozen path
                           (`-rt_refine auto`, akima, 0.1 s delta convergence).
```

## 4. Acceptance criteria, pre-registered

**Set before looking, because deciding afterwards is how four anchor statistics
pointed the wrong way during the original RT work.**

- **Step 2 passes** if the fitted slope is within **20%** of the DIA-NN CiRT line
  (10.873 on Astral) AND at least **12** anchors are inliers AND the inlier
  median residual is **≤ 60 s**. The DIA-NN line is an external reference here,
  not a target to optimise — it is used once, as a sanity bound, not iterated
  against.
- **Step 6 passes** if the slope is **stable within 10% across data halves**.
  Stability, not fit quality, per `doc/16` §4: it is immune to both measured
  failure modes (IDs oscillate, total scatter is blind).
- **The phase passes** if identifications at 1% FDR are **> 0** on both
  instruments with entrapment FDP calibrated at nominal 1%. IDs are the product
  criterion, never the optimisation target.

**Kill condition:** if the blind CiRT search cannot place ≥ 12 of 149 precursors
within 60 s of DIA-NN's RT for them, the seed does not exist and this route is as
dead as the previous four. Say so rather than loosening the bar.

## 5. Known risks

- **CiRT coverage is organism- and sample-dependent.** 64% here. A sample with
  poor coverage has no seed, and the fallback is DIA-NN's shape: wide windows →
  simple fixed-criterion seeds → refit → tighten.
- **Astral-only.** Everything in `doc/18` §9 was measured on one instrument. S08
  coverage is unmeasured, and the whole project has been caught by
  instrument-conditional results before (`MS1_COELUTION`: +17.6% on S08, +0.4% on
  Astral).
- **Linearity is assumed.** r = 0.9913 over 17 points supports a line on this
  gradient; a curved gradient would need the monotone/akima path that already
  exists.
- **The blind search may place CiRT peptides on interference.** Step 2's robust
  fit is the mitigation, and its inlier count is the diagnostic.
- **DIA-NN's CiRT line is itself an identification-derived quantity.** Using it as
  the step-2 bound imports DIA-NN's biases as a *bound*, which is acceptable; using
  it as a fitting target would not be.

---

## 6. First blind run — result, and an honest epistemic problem

Astral, 149 CiRT precursors, no calibration, whole gradient. 34:20 wall, 2.2 GB
peak RSS, 529 peak groups (258 target / 271 decoy).

**The seed works.** Of 149 CiRT targets, 92 got a peak group and 27 are
verifiable against DIA-NN:

    |ODIA apex - DIA-NN RT|:  median 1.1 s,  22/27 within 30 s

**The pre-registered gate nevertheless FAILED**, because it specified Theil-Sen
over *all* CiRT peak groups:

    RT = 1218.8 + 1.814*iRT   vs truth 907.0 + 10.873   -> 83.3% slope error

Theil-Sen breaks down above ~29% contamination, and 92 groups contain ~27 real
ones — it was asked to work at ~71%. **The estimator collapsed, not the seed.**

Selecting by ODIA's own dscore (not by truth, which would be circular):

| top-N | fitted line | slope error | inliers |
|---|---|---|---|
| **20** | **RT = 918.8 + 10.835·iRT** | **0.3%** | 13 |
| 30 | +8.433 | 22.4% | 8 |
| 40 | +6.201 | 43.0% | 8 |
| 60 | +3.467 | 68.1% | 5 |
| 92 | +1.814 | 83.3% | 10 |

Top-20 recovers the line to **0.3% on slope, 12 s on intercept**, and the
degradation with N is smooth and monotonic — a breakdown-point mechanism, not a
coincidence.

### The problem with that result

**The passing configuration was chosen AFTER the pre-registered one failed.**
That is post-hoc, and pre-registration existed precisely to prevent it. The
mechanism is demonstrated; the *test* of it is not valid.

**Therefore, pre-registered NOW, for S08 — untouched data, a genuine independent
test:**

- blind CiRT search, same settings;
- **Theil-Sen on the top 20 by dscore**, fixed in advance, no sweep;
- **PASS** = slope within 20% of S08's own DIA-NN CiRT line, ≥ 12 inliers within
  60 s, inlier median ≤ 60 s.

If S08 needs a different N, the Astral result was tuned and must be discarded.
Do not sweep N on S08 and then report the best one.

### Also worth recording

- 92 of 149 CiRT targets got a peak group but only ~27 are real. **CiRT presence
  in a predicted library does not mean presence in the sample** — the list is
  "commonly observed", not "always observed", and 57 of 149 produced no group at
  all.
- "identified nothing at 1% FDR" appears in this run's log and is **irrelevant
  here**: 149 targets cannot support an FDR estimate. Calibration needs correct
  apexes, not q-values — which is the whole reason this route works where the
  full-library bootstrap does not.

---

## 7. AMENDED GATE — frozen 2026-08-14 BEFORE the S08 output was opened

Review round 8 (codex effort max, kimi; vibe 0 bytes for the fifth time) reached
the same verdict independently: **sound seed, sound architecture, unsound gate.**
The S08 run had completed when these reviews landed and **its output was not read
until this section was written and committed.**

### What was wrong with §4

1. **The intercept was unconstrained.** Both reviewers led with this. A line of
   slope 11.0 and intercept 1,157 — 250 s from truth — passes every criterion in
   §4. Kimi: *"the plan HAS the right test, it just filed it under the wrong
   heading"* — the kill condition tests absolute placement; the pass criteria do
   not, so **step 2 could pass while the kill condition fires.**
2. **Two of three criteria were internal-consistency checks.** Inlier count and
   inlier median residual are both measured against the fitted line, so any
   coherent structure satisfies them by construction — and if the fitter defines
   inliers at 60 s, "inlier median ≤ 60 s" is tautological.
3. **Theil–Sen was disqualified by the plan's own numbers.** Breakdown point
   1−1/√2 ≈ 29.3% against a contamination rate that is not "~50%" but **82–89%**
   (27, or the reference fit's 17, of 149). Writing "RANSAC / Theil–Sen" also left
   an outcome-dependent choice inside something called a pre-registration.
4. **No mechanism-level control**, in a project where four mechanisms died to
   controls the same day.
5. **Key spaces mixed.** 149 precursors are 47 sequences; charge states of one
   peptide sit at the *same* RT and iRT and are not independent evidence. "12
   inliers" could be 5 peptides.

### The frozen criteria

**Fitter:** RANSAC only. Residual threshold 60 s, 2-point samples, 2,000
iterations, seed 0, inliers refitted by least squares. No alternative estimator.

**PASS requires all of:**

- **P1 absolute placement** — ≥ **12 distinct stripped sequences** (not
  precursors) whose predicted RT is within **60 s of that instrument's own DIA-NN
  RT**. This folds the old kill condition into the pass, which is where it
  belonged.
- **P2 intercept** — |intercept − that instrument's DIA-NN CiRT intercept| ≤ **90 s**.
- **P3 slope** — within **20%** of that instrument's own DIA-NN CiRT slope.
  *Astral's 10.873 is an Astral number and is not a bound for S08.*
- **P4 permutation null** — the identical fitter run on **999 iRT-label
  permutations** must not reach the observed inlier count. Passing means beating
  the empirical maximum-consensus null, not resembling a slope.
- **P5 shifted-CiRT control** — the same blind search on an m/z-shifted CiRT list
  must **fail** P1–P3.
- **P6 independent cross-check** — the whole-library MS2 Hough, an estimator that
  never saw the CiRT list, must agree on slope within **10%**. On Astral it gave
  10.53 / 890.6 against CiRT's 10.873 / 907.0 (~3%). Kimi is right that this was
  buried as a "retrospect" when it is the strongest evidence in the document that
  the line is real: **two mechanisms that failed certification can still confirm
  each other.**
- **P7 span** — anchors must cover ≥ 6 of 10 RT deciles. A mid-gradient fit
  extrapolates silently at the ends, where the ±60 s window then misses.

### Corrections to record

- **"27 is the realistic ceiling" was wrong** — 27 is *externally checkable*, not
  present. A blind search may legitimately place peptides DIA-NN filtered out.
- **"0/11 Biognosys, therefore absent" was wrong.** I searched the **library**,
  which is FASTA-predicted from the sample proteome — spike-ins would not appear
  there even if they were in the sample. Absence must be tested against the raw
  data by the same blind search, or established from prep records.
- **The r = 0.9913 fit used 17 of 27 and the selection rule was never stated.**
  If those 17 are the well-fitting subset, the number is a selected-sample
  statistic. This must be resolved before the figure is cited again.
- **Step 6's stability test is circular**: groups extracted within ±60 s of the
  seed line are conditioned to lie near it, so both halves reproduce the slope
  even from noise. Needs a held-out set searched over the full gradient.
- **The +600 s MS1-shift control goes vacuous after step 4** — a feature shifted
  600 s can never confirm inside a ±60 s window, so it collapses to zero by
  construction and "beats its control" stops meaning anything. Step 5 needs a
  within-window micro-shift or m/z permutation instead.
- **">0 IDs with FDP calibrated at 1%" is internally inconsistent.** A handful of
  identifications cannot validate a 1% FDP; pre-register a 95% upper bound
  (≤1.5% at nominal 1%) and a powered minimum count.
- **±60 s is asserted, not defended.** Window width has an optimum, not a maximum.
  Sweep ±30/±60/±90 judged on the frozen criteria, never on ID count.
- **If step 7 fails, run the oracle ladder** (`doc/16` §1) to discriminate "the
  CiRT line was wrong" from "the classifier never ignited". Pre-registered now as
  the failure-analysis branch so it is not negotiated after a red result.
