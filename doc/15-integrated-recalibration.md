> **SUPERSEDED BY doc/26-PLAN.md (2026-08-15).** Kept as the working record. Where this document and doc/26 disagree, doc/26 wins — it carries the corrections this one predates.

# Integrated recalibration: MS1 + MS2 mass, then RT

Status: DESIGN v5, not implemented. Written 2026-08-10, after RT recalibration
was frozen (`rt-calibration-frozen-2026-08-10`).

**Read §9, §12 and §13 first. They are the demotion logs, and they are the
content.** Each version of this document was mostly wrong, and what killed each
one is more useful than what it proposed.

- **v1** — reviewed by Kimi 0.34.0, Codex 0.147.0 and Vibe 2.24.0. All three
  rejected it; six proposals cut (§9).
- **v2** — answered them.
- **v3** — rewrote the architecture after a web survey (§11): the mass↔RT
  iteration became an additive covariate model.
- **v4** — round-2 review demoted every part of v3 that v3 was pleased with,
  including its own headline (§12).
- **v5** — the research vault was ingested, `diann.cpp` 1.7.x was read, and
  **three of §11's survey findings turned out to be wrong** (§13).

Written *before* measurement, in the discipline of `14-irt-acceptance.md`.

**A withdrawn headline, kept visible on purpose.** v3 opened by claiming the
ion-mobility term was "almost certainly the largest single miss", on the strength
of *52% of modelled mass-error variance on timsTOF* (Prianichnikov et al., MCP
2020, 19(6) 1058-1069). **That claim is withdrawn** — §12.1. It is a share of
MaxQuant's *modelled* part, measured on **DDA precursors**, and 1/K0 is not a
fragment property at all. Propagated through ODIA's own IH1 numbers it predicts
**~0.02 ppm** on total scatter. An IM term may still be worth having; it is now a
measurement with a permutation null, not a headline.

---

## 0. What is wrong with today

| | today |
|---|---|
| MS2 mass model | one global `intercept + a·log(m/z)`, fitted once |
| MS1 mass model | **none** — `Ms1Traces` borrows `fragment_ppm`, applies no offset |
| anchors | `MassCalibration::collect`'s own probe: brightest co-occurring cells, chosen without reference to any identification |
| gate | peakedness ≥ 3 against an m/z-shifted control; on failure the model is discarded and the window opens to 50 ppm |
| when | before pass 1; optionally re-measured once after the RT map |
| iteration | none |

Two real deficits, and one that turned out not to be a deficit at all.

1. **MS1 is uncalibrated.** `MassCalibration.h:90-93` records that the MS1 arm of
   the ported reference was deliberately dropped, on the reasoning that "ODIA
   extracts fragments; there is no MS1 window to size". That reasoning was true
   when written and is false now: `Ms1Traces` matches precursor peaks through
   `options.fragment_ppm` (`OpenDIAlyzer.cpp:774-776`) and MS1 co-elution is a
   scored feature. There *is* an MS1 window, it is just borrowed and uncentred.
2. **The anchor source is the probe, not the identifications.** The probe chooses
   cells by brightest co-occurrence; on a library that does not match the run it
   correctly refuses (12_80: peakedness 2.44, no offset applied), but on a
   partially-matching library it is describing a sample its own search shaped.
3. **RT-dependent mass drift is NOT a deficit — it has been measured to be
   absent.** See §2. This is what killed most of v1.

---

## 1. The measurement that constrains this design

`MassCalibration.h:104-136` records, for IH1:

```
  correction   LOG m/z: -9.35 ppm at 599 Th, +3.32 ppm per e-fold (t=7.4)
  systematic   per-m/z-bin modes sit 1.71 ppm from a constant and 0.54 ppm
               from this -- 68% of the systematic error removed. The linear
               basis leaves 0.77, so log wins on the number.
  rt drift     -0.11 ppm across 1,600 s at t=0.26, i.e. NONE
  window       3 sigma would be 12.6 ppm, WIDER than the 10 ppm in force,
               so it is rejected -- calibration may only narrow
```

and two independent confirmations: an unrelated measurement over **124 M
peak-transition hits at DIA-NN's retention times** gives −9.78 ppm and +2.79 ppm
per e-fold and **also finds no RT trend**; and sweeping the brightness cut from
0 to 0.9 moves the offset only −10.51 → −9.50 ppm.

Three consequences, all binding on this design:

- **The m/z shape is real and is logarithmic.** t = 7.4, reproduced
  independently, and it beats the linear basis on residual (0.54 vs 0.77 ppm).
  A single constant would mis-centre by ~3 ppm at both ends of the range, in
  *opposite* directions, against a ±10 ppm window.
- **The RT shape is not real, on this instrument.** t = 0.26 across the whole
  gradient, confirmed by a 124 M-hit independent measurement. Fitting RT blocks
  to it is fitting noise.
- **Total scatter is the wrong statistic and was already known to be.** The
  shape correction removed 68% of the *systematic* error while moving total
  per-hit scatter only 4.66 → 4.19 ppm. Any convergence criterion built on
  total scatter would have scored this correction as nearly worthless.

---

## 2. The mass model — additive and separable, not iterated

v1 and v2 both treated RT as something to *iterate against*. The literature
treats it as a **covariate inside the mass model**. That is a better answer and
it is the one with prior art (§11, findings 2 and 3).

**The model:**

```
    ppm = f1(m/z) + f2(RT) + f3(log I) + f4(1/K0)
```

each term a function of one variable, no interactions — MaxQuant's form, which
is the only published multi-covariate mass model and is stated verbatim as
`Δm/z = f1(m/z) + f2(RT) + f3(I) + f4(IM)`, with f1 and f2 piecewise linear and
f3 a low-degree polynomial in log intensity. **Fitted by backfitting**, term by
term, IM first (as a per-window median, then subtracted), not as one joint
least-squares.

Term by term, against what ODIA has already measured:

| term | status in ODIA | evidence |
|---|---|---|
| `f1(m/z)` | **fitted today**, log basis | real: t = 7.4, removes 68% of systematic error, log beats linear 0.54 vs 0.77 ppm residual |
| `f2(RT)` | not fitted | **measured absent on IH1**: t = 0.26 over 1,600 s, confirmed over 124 M hits |
| `f3(log I)` | **never fitted** | space-charge/saturation; MaxQuant fits it as a low-degree polynomial in log I |
| `f4(1/K0)` | **never fitted** | **52% of modelled variance on timsTOF** |

`MassResidual` already carries `intensity` (`MassCalibration.h:161-168`), so
`f3` costs one column of a fit and no new plumbing. `f4` needs the residual's
1/K0, which the extractor has at match time (the mobility gate already reads
it) but does not currently store — **add `im` to `MassResidual`**.

**Keep the log basis for f1**, against the field. No verified engine uses
log(m/z) — the implemented menu is constant shift, binned-median piecewise
curves, or linear/quadratic in raw m/z (§11 finding 1). ODIA's log basis is
unusual. It is kept because it was *measured* to beat the linear basis on this
data (0.54 vs 0.77 ppm residual to the per-bin modes), which outranks
convention. **But MaxQuant's `f1` is piecewise linear, which strictly subsumes
both**, so the honest test is log vs piecewise-linear, and that test has not been
run. Added to §8.

**`f2(RT)` is fitted but expected to be dropped by its own t-test** on IH1
(|t| > 3 to survive). Astral is a different machine and is untested. Note this
is exactly the mzRefinery design — it selects *one* of an RT-dependent or an
m/z-dependent shift per run, never both — except that an additive model does not
force the choice.

**Deleted from v1 and v2:** `MassRecalibration`'s RT-*blocked* piecewise-constant
model, and the akima interpolation across its 8 block centres. Blocks were
overfitting a signal measured at t = 0.26; the akima was additionally incoherent,
since it would have treated blocks that *inherit the global fit* as measured
control points (`MassRecalibration.h:76-83` flags them with `Block::fitted`;
v1's interpolation ignored the flag). If `f2(RT)` ever survives its t-test, it is
**piecewise linear** per MaxQuant, not blocks.

---

## 3. Where the residuals come from

### 3.1 MS2 — the hook already exists

`ChromatogramExtractor` keeps two `[row][cycle]` planes under
`collect_mass_residuals`: `ppm_num` (Σ intensity·ppm) and `ppm_den` (Σ
intensity). **`row` is the fragment**, so the per-fragment axis is already there.

`PeakGroupScorer.cpp:845-896` collapses it at the apex to `mass_ppm` (group mean)
and `mass_ppm_spread` (1.4826 × MAD). The collapse destroys the m/z axis, which
is the only axis the model is defined on.

**Change:** at that site, when a harvest sink is attached, also emit one
`MassResidual{mz, rt, ppm, intensity}` per contributing fragment. No new decode,
no new pass. ~12 fragments × ~15 k groups × 24 B (after adding `im`, with alignment) ≈ 4.3 MB.

### 3.2 MS1 — a new path

`Ms1Traces.cpp:82,87` computes `tol` and tests `|it->mz - m| > tol`; the
deviation is in scope and discarded, exactly as the MS2 deviation was before
`-collect_mass_residuals`.

**Change, with one correction forced by review:** the retained residual must be
the residual of **the peak that won the max**, not an intensity-weighted average
over all matching peaks. `Ms1Traces.cpp:99` stores `Max, not sum` because a
mobility-merged frame holds the same ion in several scans. A weighted mean over
matching peaks would report a ppm at which no ion exists: given a 1,000-count
ion at −2 ppm and an 800-count interferent at +4 ppm, the trace records the
−2 ppm ion and a weighted mean would record +0.67 ppm. **Carry the argmax peak's
ppm alongside the max intensity.**

**Threading:** `Ms1Traces::build` is not threaded — `Ms1Traces.cpp:72-107` is a
plain serial loop over 64-spectrum blocks, no OpenMP pragma, no `std::thread`.
So there is no race. The write index is already `slot * bins_ + (b + s)`, i.e.
per-(precursor, cycle), so it would stay race-free if parallelised **over
spectra**. Parallelising over *peaks within* a spectrum would race and must not
be done.

### 3.3 The null — the m/z-shifted control, NOT decoys

v1 proposed using decoy peak groups' fragments as the null. All three reviewers
rejected this and they are right.

Pseudo-reverse decoy fragments query a **different m/z distribution** than
target transitions (tryptic y-ions are not uniform over 200–2000 Th), so they
sample a different local peak density; and the component contaminating *target*
groups is structured interference — homologous peptides, isotope and
charge-state neighbours, co-eluting species with real ion-series overlap — which
decoys do not reproduce. The two distributions differ in width and in tail shape,
and the direction of the error is **not determined** (Kimi and Codex both note it
can go either way depending on where the decoy transitions land; Vibe asserted a
direction, and that assertion is the one claim among the three reviews I do not
accept).

Misspecifying `f_null` makes both `π` and the scale of `f_true` non-identifiable.

**Use the m/z-shifted control of the target transitions** — which
`MassCalibration` already builds (`MassCalibration.h:78-84`) and already
validates against (the IH1 gate is 4.34 target vs 2.98 control). It preserves the
target m/z distribution and the local peak environment, which is the whole point
of a null.

---

## 4. MS1 — a separate WINDOW for certain, a separate MODEL as an experiment

v1 and v2 asserted "two fully independent models". The survey does not support
that as the default (§11 finding 4). What engines actually do is **one shared
correction model, with a separately estimated MS1 window**:

- OpenSWATH fits its regression on **fragment anchors only** and applies it to
  all maps, but collects MS1 residuals separately (`delta_ppm_ms1`,
  `SwathMapMassCorrection.cpp:491,577-598`) and derives a **distinct precursor
  window** (`estimateWindow` at `:632`, separate from the fragment window at
  `:614`).
- mzRefinery applies its single chosen model to both MS1 spectra and MS2
  precursor m/z.
- MaxQuant's published model is precursor-only, with a global fragment tolerance.

So the claim "MS1 and MS2 need different error *functions*" is **not
established** — and note that a verifier panel specifically refuted the inference
from DIA-NN's differing `--mass-acc-ms1` / `--mass-acc` defaults (Astral 4 vs
10 ppm), because different *tolerances* do not demonstrate different *error
functions*, and DIA-NN publishes no MS1 mass model at all.

**Revised:**

1. **Ship the separate MS1 window unconditionally.** It is well-precedented, it
   is the actual deficit (`Ms1Traces` currently borrows `options.fragment_ppm`,
   `OpenDIAlyzer.cpp:774-776`), and it needs no new model.
2. **Fit a separate MS1 model as an experiment**, reported alongside the shared
   model, with the falsification test already in §10: if its coefficients agree
   with MS2's within their confidence intervals, delete it and share the model.

MS1 anchor count is ~1 per accepted precursor rather than ~12, so the MS1 model
gets `f1` only — no RT, intensity or IM terms, which it cannot support.

---

## 5. Sizing the window

v1 proposed "the 99% coverage quantile of the deconvolved true component". All
three reviewers rejected transplanting coverage-at-window from RT to m/z, on the
same argument, and it is a good argument:

> Widening an RT window admits a neighbouring peptide — a real signal the
> classifier can weigh. Widening an m/z window admits chemical noise for every
> fragment in every spectrum. Coverage is monotone in width and contains no term
> penalising noise, so coverage alone is optimised by "as wide as possible".

**Revised, and now with a precedent to copy rather than an invented scheme.**
OpenSWATH sizes its window as the **99th percentile of residuals × a 1.3 padding
factor** (`estimateWindow(residuals, quantile=0.99, full_width=true,
padding_factor=1.0)`, driven by `mz_estimation_percentile` default 99.0 and
`mz_estimation_padding_factor` default 1.3). Note the padding is a deliberate
*widening past* the quantile — the opposite of v1's instinct to take the quantile
as an upper bound.

- **Base:** p99 of held-out residuals × 1.3, matching OpenSWATH.
- **Constraint (noise):** false-match density in the **m/z-shifted control** at
  that width must not exceed its value at the current window by a stated factor.
  This is the term coverage lacks; no surveyed engine has it, so it is ODIA's
  own and must be justified by measurement, not assumed.
- **Rails kept:** may only narrow; ≤ 2× narrowing per round; 5 ppm floor.
- **Robust statistics throughout**: median for location, MAD for scale. SD
  appears nowhere as the primary residual statistic in any surveyed engine, and
  ODIA's own IH1 numbers show why — the shape correction removed 68% of
  systematic error while total scatter moved only 4.66 → 4.19 ppm.

The existing "3 sigma would be 12.6 ppm, WIDER than the 10 ppm in force, so it is
rejected" behaviour is preserved by construction.

**Widen the bootstrap.** Calibration windows in the field run **7×–75×** the
final tolerance: DIA-NN `--mass-acc-cal` defaults to **100 ppm** against a 4–15
ppm search window; MaxQuant used 70 ppm precursor / 40 ppm fragment against a
post-calibration precursor MAD of 0.94 ppm. ODIA's 50 ppm bootstrap against a
10 ppm window is **5×**, below the entire documented range. **Propose 100 ppm**,
matching DIA-NN. Caveat carried from the survey: DIA-NN's *stated* rationale for
the wide window is instrument-miscalibration tolerance, not selection-bias
avoidance — treat this as an existence proof of the mechanism, not an endorsement
of the statistical argument.

---

## 6. Iteration — and the evidence against it

**State this plainly, because it contradicts the brief this design was written
to.** The survey found **no verified source** that iterates mass calibration
against RT calibration to joint convergence — no iteration count, no convergence
criterion, in any engine or paper. The two documented architectures are
(a) RT as a covariate *inside* the mass model (MaxQuant, §2), and (b) RT- and
m/z-dependence as *mutually exclusive* alternatives, one selected per run
(mzRefinery). Everything else is single-pass: loose → fit → tighten.

And the one system that does auto-tune is cautionary. DIA-NN's README, on its own
automatic mass-accuracy determination:

> "This optimisation is inherently noisy: even replicate injections may not
> produce identical results, and therefore the analysis results will depend on
> which run is first in the list. It is preferable to fix these three parameters
> to values known to be optimal for a particular LC-MS setup."

That is an against-interest admission by the tool's own author, and it is direct
evidence for the runaway-narrowing and overfit-to-an-early-biased-set failure
modes §6.5 worries about.

**So: adopt the additive model of §2 as the primary design, which needs no
iteration at all.** The residual loop below is kept because it is cheap, because
the confirmatory pass is needed regardless, and because the covariate model and
the loop are not exclusive — but it is capped hard and it is *not* where the win
is expected to come from. If §2's terms do their job, §6 should converge in one
inner round, and §10 already names that outcome as the one that deletes the loop.

Note that v2's architecture — one wide extraction, fit, one confirmatory pass at
the tightened window — had already converged on the field's documented
loose→fit→tighten shape before the survey was read. The survey confirms it
rather than changing it.

### 6.1 Extract once, iterate on retained residuals

v1 said "re-extract each round", capped at 8 rounds. All three reviewers called
this the critical practical flaw, and they are right, though their cost figures
were guesses (Codex cited `MassCalibration.h:86-87` for "~10 minutes per pass";
those lines are a *what-is-not-ported* comment and carry no such number). The
defensible figure from the repo: the probe touched **11.9% of IH1's spectra in
84.6 s with only 3,000 sampled precursors**. A full pass over the whole library
is materially more than 84.6/0.119 ≈ 12 min. Eight of those is around
1.5 h per file, on top of a ~10.3 GiB decode floor paid each time.

**It is also unnecessary.** Fitting the mass model needs only
`(mz, rt, ppm, intensity)` per fragment. Once those are retained from **one wide
extraction**, evaluating any *narrower* candidate window is a **filter over
retained residuals**, not a re-extraction. So:

```
  ONE wide extraction  (bootstrap width per SS5, library iRT, collect residuals)
        |
        v
  score + FDR at 1%
        |
        v
  inner loop, NO decode:
        fit MS2 model  -> fit MS1 model  -> size both windows
        -> refit RT map + per-run refinement (frozen path)
        -> recompute residuals under the new models by
           re-centring the RETAINED residuals
        -> report, test convergence
        |
        v
  ONE confirmatory extraction at the converged windows
        |
        v
  if the accepted set moved materially, at most ONE more outer round
```

Two decodes, not eight. The cap on outer rounds is **2**.

What this does *not* capture: a narrower window changes chromatogram shape and
therefore co-elution scores and therefore which groups are accepted. That
feedback is real and is exactly what the confirmatory extraction measures. The
design's claim is that this feedback is second-order relative to the model fit —
and that claim is falsifiable at the confirmatory pass (§8).

### 6.2 Held-out validation is mandatory

Freezing identities does not create held-out evidence — a flexible model can
improve in-sample σ while worsening unseen fragments, and a stable Δσ can be a
stable *wrong* fixed point.

**Split at the stripped-sequence level** (the split the RT work already uses):
fit on 85%, validate on 15%. Every reported σ and every convergence test is on
the **held-out** part.

### 6.3 The frozen evaluation set, and why σ is now well-defined

v1 froze the reporting set after round 1 and counted lost members as "infinite
residual". Kimi and Codex both showed that is coherent for *coverage* (binary)
and incoherent for *σ* (a robust scale cannot absorb infinities; MAD ignores them
below 50% loss and becomes undefined above it). Codex's scenario is the clean
one: freeze 1,000 precursors, 800 at σ≈2 ppm and 200 difficult ones at σ≈8 ppm;
narrowing loses the 200, and σ "improves" 2.0 → 1.9 purely by truncation.

**The single-wide-extraction architecture of §6.1 dissolves this.** Residuals
come from the fixed wide extraction, so **every frozen member has a residual in
every round regardless of the current window**. Nothing is lost, nothing is
infinite, and σ is computed on a constant denominator over a constant
population. Both Codex and Kimi independently proposed a "constant wide shadow
extraction" as the fix; making the *primary* extraction that shadow is strictly
cheaper than adding one.

Coverage is still reported separately and still counts losses as failures.

Freezing does **not** make the set representative — it inherits round 1's
selection bias. It only prevents the denominator artefact (the one that faked a
"3.3× narrowing" during the RT work when the anchor population grew 1,183 →
1,408).

### 6.4 Convergence

On deltas, over the held-out frozen set:

```
    stop when  |Δσ_mz(MS2)| < 0.1 ppm  and  |Δσ_mz(MS1)| < 0.1 ppm
          and  |Δσ_RT| < 0.1 s
          or   inner round == 6, or outer round == 2
```

with a 1% relative floor as the alternative trigger, mirroring
`-rt_converge_rel`.

**"Best round" is pre-registered, lexicographically**, so oscillation recovery is
reproducible: (1) held-out m/z coverage within 1% of the best observed;
(2) then held-out RT coverage likewise; (3) then smallest MS2 window;
(4) then smallest MS1 window. σ is a *diagnostic*, not a tie-breaker — the IH1
measurement above is precisely a case where σ barely moved while the correction
did real work.

### 6.5 Round 1 is an initialisation, and its bias is bounded, not removed

The round-1 anchors are selected by the library's linear iRT and by a classifier
whose features include `MASS_SPREAD` — so mass-coherent groups are preferentially
accepted, and the sample is conditioned on the very quantity being fitted. A wide
m/z window removes *m/z truncation*; it does not remove *RT-conditioned
missingness*. Codex's case: true error −8 ppm early and +8 ppm late, library iRT
90 s wrong late, so the late population is under-sampled and a run whose true
trend crosses zero can be fitted flat.

This cannot be solved inside the loop. It can be bounded:

- **Round 1 fits at 1% FDR, not 5%.** v1 proposed 5% to widen the anchor base.
  Kimi and Codex both showed that is worse: it adds interference-rich groups
  rather than recovering true peptides the RT search never found, and 5%
  contamination at a +3 ppm systematic shifts the fitted median by ~0.15 ppm —
  larger than the 0.1 ppm convergence tolerance. It also would have frozen the
  evaluation set at 5% contamination while later rounds fit at 1%.
- **Report anchor counts stratified by RT and by m/z decile every round.** An
  RT region with no anchors is then visible rather than silently fitted.
- **Refuse to fit any m/z or RT stratum below a minimum count**, inheriting the
  global fit there — and mark it, so it is never used as a control point.
- **External validation via the entrapment library** (already on the backlog).
  That is the only evidence in this project that is not conditioned on ODIA's own
  acceptance.

---

## 7. The per-iteration report

One row per round, to the log and to a TSV under `-out_recalibration`:

```
round  n_fit  n_heldout  ms2_off  ms2_slope  ms2_sig_ho  ms2_win  ms2_cov_ho
       ms1_off  ms1_slope  ms1_sig_ho  ms1_win  ms1_cov_ho
       rt_sig_ho  rt_win  rt_cov_ho  ctrl_density  ids
       d_sig_mz2  d_sig_mz1  d_sig_rt
```

- `*_ho` are **held-out** (§6.2) on the **frozen** set (§6.3).
- `ctrl_density` is constraint 2 of §5, from the m/z-shifted control.
- `d_*` are the convergence quantities, so the stopping decision is auditable
  from the report alone.
- `ids` is reported and is **not** the objective and **not** the stopping
  criterion.

Plus, per round: anchor counts by RT decile and m/z decile (§6.5), and the
target-vs-control residual histograms for MS1 and MS2, so the mixture fit is
checked by eye rather than trusted.

---

## 8. Implementation order

Reordered in v3 so that the terms with the strongest evidence come first.

1. **MS2 per-fragment harvest** — emit `MassResidual` at
   `PeakGroupScorer.cpp:845-896`, and **add `im` to the struct**. Refit today's
   global log model from identification-derived anchors instead of the probe.
   *Measure: does it beat the probe's model on held-out fragments?*
2. **`f4(1/K0)` — the ion-mobility term.** Highest expected value of anything in
   this document: 52% of modelled variance on timsTOF, and IH1 is timsTOF.
   Fitted first and subtracted, per MaxQuant's backfitting order. *Measure: how
   much of the 4.19 ppm residual scatter does it remove?*
3. **`f3(log I)` — the intensity term.** Low-degree polynomial in log intensity.
   `MassResidual::intensity` already exists, so this is one column.
   *Measure: same.*
4. **`f1` basis test** — log(m/z) versus piecewise-linear in m/z. ODIA's log
   basis beat the linear basis (0.54 vs 0.77 ppm) but has never been tested
   against the piecewise-linear form that every surveyed engine uses and that
   subsumes both. *Measure: held-out residual to per-bin modes.*
5. **`f2(RT)` and its t-test** — fit the term, report |t|, drop it below 3. On
   IH1 this should drop (t = 0.26). *Measure: does it survive on Astral?*
6. **Widen the bootstrap to 100 ppm** (§5) and re-measure everything above.
   *Measure: do the fitted coefficients move? If they do, the 50 ppm bootstrap
   was censoring and every earlier number in this project is suspect.*
7. **Window sizing** — §5, p99 × 1.3 plus the control-density constraint.
   *Measure: does it recover the 4,969 that `-mass_width_from_ids apply`
   destroyed (Astral 4,969 → 2,422)?*
8. **MS1 window** — §4, separate window, shared model. Unconditional.
9. **m/z-shifted control as `f_null`** + the mixture fit. *Measure: judge against
   the entrapment population, NOT by "deconvolved σ tracks raw σ" — two wrong
   estimators can track each other.*
10. **MS1 model as an experiment** — §4, with the §10 falsification test.
11. **The inner loop** — §6.1, no decode. Then the confirmatory pass.
12. **The report** — §7, written before step 11 produces numbers.

Steps 1–10 ship without any loop and are where the evidence says the win is.
Step 11 is where the risk is, has no prior art, and is last.

---

## 9. What v1 got wrong, recorded so it is not re-proposed

| v1 proposal | why it is gone |
|---|---|
| RT-blocked mass model as the default | RT drift measured at t = 0.26 on IH1 and confirmed absent over 124 M hits. Fitting noise. Now opt-in behind a t-test. |
| Akima across 8 block centres | Overfits 8 noisy points; and would treat *inherited* blocks as measured control points, manufacturing end curvature from fallback policy. |
| Decoy fragments as `f_null` | Wrong m/z sampling and wrong contamination structure; makes `π` and `f_true` non-identifiable. Use the m/z-shifted control. |
| 99% coverage as the window objective | Coverage is monotone in width and has no noise term; its optimum is "as wide as possible". Now a constraint, paired with control-density. |
| Re-extract every round, cap 8 | ~12+ min per pass on IH1 and a 10.3 GiB decode floor each time. Now one wide extraction + one confirmatory pass. |
| Lenient 5% FDR in round 1 | Adds interference, not missing truths; 5% contamination moves the fit by more than the convergence tolerance. Now 1%. |
| σ on a frozen set with lost members counted as infinite | Coherent for coverage, incoherent for a robust scale. Dissolved by fitting from one fixed wide extraction. |
| Convergence measured in-sample | A flexible model improves in-sample σ while worsening held-out. Now held-out by stripped sequence. |
| "revert to the best round" | Undefined over five competing metrics. Now pre-registered lexicographic. |
| Intensity-weighted MS1 ppm | Conflicts with the `Max, not sum` trace semantics; reports a ppm at which no ion exists. Now the argmax peak's ppm. |

## 10. What would falsify v2

- If the identification-derived global model does not beat the probe's model on
  held-out fragments (step 1), the anchor-source change is not the win and the
  whole phase reduces to adding MS1.
- If |t| < 3 for the RT trend on every instrument tested (step 2),
  `MassRecalibration` should be **deleted**, not merely left unwired.
- If the MS1 fitted offset and slope agree with MS2's within their confidence
  intervals (step 5), one shared model is correct and the second is waste.
- If σ converges in the first inner round (step 7), the loop is unnecessary and
  this ships as a single pass with a better anchor source and an MS1 arm — which
  would still be a clear win over today.
- If the confirmatory extraction moves the accepted set materially (step 7), the
  "feedback is second-order" claim of §6.1 is false and the cheap inner loop is
  not a valid substitute for re-extraction. **This is the assumption v2 is most
  exposed on.**

Every one of these outcomes deletes code.

---

## 11. The survey this design was rewritten against

A 101-agent multi-source survey with 3-vote adversarial verification per claim,
run 2026-08-10. Findings that changed the design, with their vote records.
Claims that were voted down are recorded too, because they are the ones most
likely to be re-proposed.

**F1 (3-0).** No verified production implementation uses a log(m/z) term or a
spline/LOWESS in m/z. The implemented menu is: no-op or constant global ppm
shift; piecewise binned-median curves; linear or quadratic regression in *raw*
m/z. OpenSWATH's 8 `mz_correction_function` values all use raw `exp_mz` with
`Y = a + bX + cX²`. LOWESS and b-spline appear in OpenSwathWorkflow only under
`RTNormalization:alignmentMethod` — splines are used for RT, never for m/z.
mzRefinery is fully non-parametric (median shift over 0.5 ppm bins; scan-time
bins of 75 s; m/z bins of 25 units). MaxQuant's `f1(m/z)` is piecewise linear.
→ *ODIA's log basis is unusual; kept on measurement, tested in step 4.*

**F2 (3-0).** MaxQuant: `Δm/z = f1(m/z) + f2(RT) + f3(I) + f4(IM)`, additive and
separable, no interaction terms; f1 and f2 piecewise linear, f3 a low-degree
polynomial in log intensity; fitted by backfitting with the IM term determined
first as a window median and subtracted. `f4(IM)` explains **52% of modelled
variance** on timsTOF. Scope caveats: MS1/precursor only, and DDA (PASEF) — never
described for fragments. → *§2, and step 2.*

**F3 (3-0).** **Negative finding, directly against the original brief.** No
verified source iterates mass against RT to joint convergence. No iteration
count, no convergence criterion, anywhere. → *§6.*

**F4 (2-1/3-0).** One shared correction model with a *separately estimated MS1
window* is the pattern — not two independent models. **Refuted 0-3:** the
inference that DIA-NN's differing `--mass-acc-ms1`/`--mass-acc` defaults prove
different MS1/MS2 error functions. → *§4.*

**F5 (3-0).** Calibrant anchors are the sample's own peptides from a bootstrap
search, selected by a **simple fixed criterion, not an FDR-controlled set**
(MaxQuant: Andromeda score ≥ 70), with no spike-in standards. mzRefinery requires
≥ 500 confident points before attempting any dependent shift. **Refuted (0-3 and
1-2):** that mzRefinery uses a q < 0.01 cutoff — the paper states no FDR cutoff,
so do not cite one. → *tension with §6.5's 1% FDR choice, noted not resolved.*

**F6 (3-0).** The only documented remedy for tolerance-fitted-to-selected-IDs is
"use a much wider bootstrap window", 7×–75× the final tolerance. **No verified
engine implements an EM/mixture model, a decoy-based null, or any explicit
selection-bias correction.** → *§5 widens the bootstrap; and note this makes
ODIA's §3.3 mixture fit novel, hence step 9's demand that it be judged against
entrapment rather than against itself.*

**F7 (3-0/2-1).** Windows come from residual **quantiles**, not multiples of SD:
OpenSWATH p99 × 1.3 padding. MaxQuant abandons a global precursor tolerance
entirely for per-PSM individual tolerances derived from the variability of
multiple mass measurements within the 4D feature — precursor only; fragments keep
one global tolerance. OpenSWATH silently no-ops if fewer than 3 anchors are
available. → *§5.*

**F8 (3-0).** Robust statistics dominate: median for location, MAD for scale and
for model selection. SD appears nowhere as the primary residual statistic.
mzRefinery selects among its three forms by `percentImprovementSmoothedMAD > 3.0`.
**Caveat:** mzRefinery is a precedent for a robust estimator of *location* only —
the paper has no MAD, no scale estimation, no tolerance setting. Do not cite it
for MAD-over-SD tolerance setting. → *§5.*

**F9 (3-0).** DIA-NN's own README on its automatic mass-accuracy determination:
"inherently noisy", replicate injections may not reproduce, results depend on run
order, "it is preferable to fix these three parameters". Its auto-tightening
"may, in rare cases, misinterpret MS/MS data with atypical characteristics."
→ *§6; the strongest available evidence against auto-narrowing loops.*

**F10 (3-0/2-1).** Magnitude anchors: across 91 Orbitrap and QqTOF files, per-file
**median** mass error spanned −2.8 to +8.4 ppm before calibration and −0.59 to
+0.28 ppm after. MaxQuant on timsTOF moved precursor MAD 2.79 → 0.94 ppm.
**Four caveats bound their use:** those are medians (location), not widths, so
they license nothing about window width; they are in-sample and partly
tautological; both datasets are precursor-only and DDA, with no fragment or DIA
benchmark; and the MaxQuant 3× is dominated by the IM term, so on non-IM
instruments the same model yields materially less.

### What the survey did NOT settle

- No fragment/MS2 calibration benchmark exists in any surveyed source. Every
  quantitative anchor above is precursor-only and mostly DDA. **ODIA's MS2 arm
  is operating without prior art**, which is a reason for the entrapment
  validation in §6.5, not a reason for confidence.
- Whether an iterated mass↔RT loop beats a single pass is **unknown**, not
  refuted. There is no evidence either way, only the absence of anyone trying and
  DIA-NN's warning about a related mechanism.

---

## 12. v4 amendments — round-2 adversarial review

Round 2 (Kimi 0.34.0, Codex 0.147.0) reviewed v3 at `4b67160`. It did not reject
the architecture, but it removed the confidence from the part v3 was most
pleased with. **Every amendment below is a demotion.**

### 12.1 The IM claim was overstated. Measure before building.

v3's headline — "the IM term is almost certainly the largest single miss" — does
not follow from the cited number and is **withdrawn**.

- 52% is a share of MaxQuant's *modelled* variance, not of total mass error, and
  its denominator is that model on that dataset.
- It is a **DDA PASEF precursor** result. ODIA would apply it to **DIA
  fragments**, where 1/K0 is not a fragment property at all: fragments inherit
  the precursor's mobility (`ChromatogramExtractor.cpp:1111` tests a fragment
  peak against `x.precursor_im`). So `f4` is a *per-precursor* covariate wearing
  a per-fragment label.
- Propagating it through ODIA's own IH1 numbers gives an expected improvement in
  total scatter of **~0.02 ppm** (random variance ≈ 4.19² − 0.54² ≈ 17.3 ppm²;
  systematic residual after a 52% cut ≈ 0.37 ppm; new total ≈ 4.17 ppm). That is
  five times *below* v3's own 0.1 ppm convergence tolerance.

**Step 2 is therefore a measurement, not an implementation:**

1. Take existing fragment residuals, remove the fitted `f1(m/z)`.
2. Estimate the 1/K0 slope **within narrow m/z × isolation-window × charge
   strata** — the marginal trend is confounded, the conditional one is not.
3. **Permute 1/K0 within those same strata as the null.** If the conditional
   slope survives permutation, the effect is real.
4. Compare with `MobilityCalibration` on and off (see 12.2).
5. Validate by **held-out precursor**, never held-out fragment — fragments of one
   precursor share mobility, apex and interference, so a fragment-level split
   leaks.

Only a reproducible *conditional* trend justifies building `f4`.

### 12.2 A confounder in the existing call order

`applyMassCalibration_` runs at `OpenDIAlyzer.cpp:790`; `applyMobilityCalibration_`
runs at `:811`. **The mass model is fitted on residuals collected before the
mobility axis is centred.** `MassCalibration` gates its matches on an IM window
around the *library* 1/K0, so any IM-dependent mass trend measured today is
confounded with the uncorrected mobility gate — and `MobilityCalibration` may
already be absorbing it by fixing which peak is selected.

This is a pre-existing property of the pipeline, not a v3 invention, and it means
**`f4` cannot be honestly measured without reversing the order or re-measuring
after mobility calibration.** Added as a precondition to step 2.

### 12.3 Backfitting order is not innocent; fit m/z first

v3 copied MaxQuant's IM-first order. Both reviewers independently showed it can
alias: m/z and 1/K0 covary through CCS, so fitting IM first and subtracting can
assign several ppm of the *measured, strong* `f1` shape (t = 7.4) to `f4`, and a
subtracted-and-frozen term cannot give it back.

- Fit **`f1` first** — it is ODIA's established term, not a borrowed one.
- Use **true backfitting** (cycle to stability) with **identifiability
  constraints** — each term zero weighted mean — not one-pass subtraction.
- Test all three of: IM-first, m/z-first, and a joint penalised fit with centred
  terms. **If they disagree materially, the effect is confounding, not physics.**

### 12.4 `f3(log I)` is endogenous, not "one column"

Fragment intensity determines which peak wins the match, whether the fragment
contributes, the group score, and ultimately FDR acceptance. It is simultaneously
a proposed predictor **and part of the selection mechanism**.

Failure mode: low-intensity genuine fragments get displaced by higher-intensity
background peaks nearby; `f3` then learns an apparent intensity-dependent ppm
shift **from peak-selection error**, and applying it moves genuine low-intensity
fragments away from where they belong.

Also, MaxQuant's `f3` models space charge on the *intact* ion; fragment intensity
is governed by fragmentation propensity, not by the current that caused the
shift. Two fragments of one precursor would receive different `f3` corrections,
distorting their relative masses.

`f3` is demoted to a candidate covariate requiring the same conditional and
permutation tests as `f4`.

### 12.5 §6.1's "second order" claim is false; the inner loop is model-selection only

Both reviewers produced the same counter-example, and it is decisive. A fragment
has the true peak at −8 ppm and a stronger interferer at +18 ppm; the wide
extraction sees both and the retained residual is a weighted mixture near 0 ppm.
Once the model centres at −8 ppm and the window narrows, the interferer leaves —
and **re-extraction yields a clean residual near 0, while re-centring the
retained residual yields +10 ppm for a peak that no longer exists.** The apex
moves, the co-elution score changes, and acceptance can flip.

So retained residuals are **not a sufficient statistic** for a narrower
extraction. Amendments:

- The inner loop fits and *compares* models. It does **not** size the final
  window from re-centred residuals alone.
- The confirmatory pass must compare **fitted coefficients, residual
  distributions, apex movement and contributing-peak identity** — not merely
  "did the accepted set move".
- If those disagree, the outer cap of 2 is **not** sufficient, and the honest
  fallback is the field's plain single pass: bootstrap wide → fit once → search
  narrow. §11 finding 3 says that is what everyone does anyway.

### 12.6 The bootstrap widening is largely inert as specified

`search_ppm` defaults to 50 but **`gate_ppm` defaults to 30**
(`MassCalibration.h:184,201`), and the gate, mode and scale estimator all work
over `gate_ppm`. Widening `search_ppm` to 100 while leaving `gate_ppm` at 30
collects residuals the fit then ignores — and the badly-miscalibrated instruments
that motivate a 100 ppm search are exactly the ones whose peak lies outside a
30 ppm gate. **Widen both, or widening buys nothing.**

**A disagreement, resolved against Codex.** Codex argued the m/z-shifted control
band overlaps the target band at 100 ppm, on `S ≤ 2W`. That is wrong:
`decoy_shifts` are `{7.33, −7.19}` **Th**, not ppm
(`MassCalibration.h:365`). At 2000 Th a 100 ppm half-width is 0.20 Th, so
`2W = 0.40 Th` against a 7.33 Th shift — the control sits **18× clear** of the
target band and does not overlap at any width under consideration. Kimi's framing
(the real cost is match-loop work and mis-assignment opportunity) is the correct
one.

### 12.7 The 0.1 ppm convergence tolerance is not usable

It fails in **both** directions:

- **Too loose** to see the win it was written for: §12.1's propagation puts
  `f4`'s effect on total scatter at ~0.02 ppm.
- **Not resolvable**: `SE(σ̂_MAD) ≈ 1.17 σ/√n`. At σ ≈ 4.19 ppm, MS1 has
  ~15,000 precursors × 15% held out ≈ 2,250 residuals → **SE ≈ 0.10 ppm**, i.e.
  the whole tolerance. MS2's nominal ~27,000 residuals give SE ≈ 0.03 ppm, but
  fragments of one precursor share apex, mobility and interference, so with a
  design effect of 10–12 the effective SE approaches 0.1 ppm there too.
- And the window it feeds is noisier still: the p99 quantile at n ≈ 7,000,
  σ ≈ 4 ppm has SE ≈ 0.2 ppm, ~0.25 ppm after the 1.3 padding. Converging σ to
  0.1 ppm while sizing the window from a quantity with 0.25 ppm sampling noise is
  incoherent.

**Replace the criterion.** Converge on the **stability of the correction function
evaluated on a fixed m/z grid**, with uncertainty from a **precursor-cluster
bootstrap**, and stop when successive corrections agree within that bootstrap
uncertainty. This also matches §1's warning that total scatter is the wrong
statistic — the IH1 correction removed 68% of systematic error while moving
scatter 4.66 → 4.19.

### 12.8 Smaller corrections

- **`MassResidual` grows** to ≥20 B, 24 B aligned, once `im` is added; the
  worked estimate is now ~4.3 MB, and 15,000 groups is not an upper bound at a
  widened bootstrap. Use measured record counts.
- **The MS1 model is underpowered.** ~1 residual per precursor means the §10
  falsification test ("coefficients agree within CIs") will often fail to reject
  a shared model even when MS1 needs its own. The separate MS1 *window* (§4.1)
  is unaffected and still ships.
- **No interaction terms in a co-isolated world.** Fragments of co-eluting
  precursors share RT and IM and their intensities are not independent; an
  additive model can attribute cross-precursor interference structure to `f3` or
  `f4`. Another reason both terms need permutation nulls.

### 12.9 What survived round 2 intact

Keeping the measured log-m/z term pending a held-out comparison against
piecewise-linear; making `f2(RT)` earn inclusion per instrument; the separate MS1
window; using the m/z-shifted control rather than decoys as the null; and the
confirmatory extraction as the falsification point — provided its comparison is
broadened per §12.5.

---

## 13. v5 corrections — the reference implementation was readable all along

Added after the research vault was ingested (2026-08-10). `diann.cpp` 1.7.x @ `cff0408` was
present in the reference directory, under CC BY 4.0, and had been sitting unread while §11's
survey inferred the same facts from public documentation. **Three of §11's findings are wrong.**
Read as documentation only; cite Demichev et al., *Nat Methods* 17:41–44 (2020).

### 13.1 The actual function

`predicted_mz` (`diann.cpp:6685`) and `predicted_mz_ms1` (`:6697`):

```
    corrected_mz = mz
                 + t[0] * mz^2                        <- ONE global quadratic term
                 + interp_RT( t[1+2i] + t[2+2i]*mz )  <- per-RT-bin intercept and slope in m/z
```

RT enters by **linear interpolation between adjacent bin centres**, weighted by distance, clamped
outside the outermost centres. The correction is in **absolute m/z** (`return s + mz`), so in ppm
the shape is `t0*mz + b + a/mz` — a different family from both the log basis and a linear-in-ppm
basis.

Note what this does to §2 and §12.3: the reference implementation makes the **curvature global and
the offset/slope RT-local**, and interpolates them **linearly**. v1 proposed akima across block
centres — *more* flexible than the reference. v2/v3 deleted RT dependence altogether — *less*. The
reference sits between, and it defaults to a single bin (below), so out of the box it is exactly
what ODIA already does.

### 13.2 §11 finding 4 is REFUTED — MS1 and MS2 get separate models

`MassCorrection` and `MassCorrectionMs1` are **separate coefficient vectors**, with separate bin
counts (`MassCalBins` / `MassCalBinsMs1`), separate centres, and a separate application site
(`:6774`). Two independent models, unconditionally.

The verifier panel specifically refuted (0-3) the inference that DIA-NN fits separate MS1/MS2
models, on the correct grounds that differing tolerances do not prove differing error functions.
The reasoning was sound and the conclusion was wrong: **absence of evidence was scored as evidence
of absence.**

**Consequence for §4:** demoting the separate MS1 model to "an experiment" rested on a false
premise. Restore it to the default. The MS1 *window* was already unconditional and stays so.

### 13.3 §11 finding 6 is REFUTED — the bootstrap ratio is 5×, not 7–75×

`CalibrationMassAccuracy = 100 ppm` (`:217`) against `GlobalMassAccuracy = 20 ppm` (`:218`) and
`GlobalMassAccuracyMs1 = 20 ppm` (`:219`). **5×.** The 7–25× figures come from DIA-NN **2.x**
README defaults (4–15 ppm search windows) — a different program. MaxQuant's "74×" divides a
*window* by a *MAD*, which are not the same kind of quantity.

**Consequence for §5:** ODIA's existing 50 ppm bootstrap against a 10 ppm window is **also 5×**,
i.e. already matching the reference. **The proposal to widen to 100 ppm is withdrawn** — it loses
its justification, and §12.6 had already shown it was inert anyway while `gate_ppm` stays at 30.

### 13.4 §11 finding 3 is PARTLY REFUTED — DIA-NN does iterate

Calibration in 1.7.x is **interleaved with the search over 12 iterations**, each feature declaring
the earliest iteration at which it may be used and at which it may be fitted. So "no verified
source iterates mass calibration against RT" is false as stated.

**What survives is narrower and still supports §6:** no source documents a *convergence criterion*
for such a loop. DIA-NN runs a fixed 12 and stops. ODIA's problem — when to stop — remains
unsolved by prior art.

### 13.5 One finding the source strengthens rather than refutes

`MassCalBinsMax` defaults to **1** (`:162`) — **RT binning is off by default**. That is convergent
with ODIA's own IH1 measurement (`rt drift −0.11 ppm at t = 0.26, i.e. none`) and with §2's
decision to make `f2(RT)` earn inclusion per instrument. Two independent programs concluding the
same thing about the same axis is the strongest evidence in this document.

### 13.6 The divergence this project is choosing deliberately

After fitting, DIA-NN runs a **grid search over mass accuracy, stepping ×1.2, maximising
identifications at 10% FDR**. That is a fourth window-sizing option beyond k·σ, p99 and coverage —
and it optimises **exactly the quantity this project forbids as an objective**.

This is a deliberate divergence, not an oversight. The reason is Wen et al. 2025 (no DIA tool
consistently controls peptide-level FDR; DIA-NN's true precursor FDP > 2.3% at a nominal 1%), so
maximising nominal IDs is not evidence of a better window. **But it should be recorded that the
reference implementation disagrees with us here**, and that our position costs us the one
window-sizing rule that is known to work in practice.

### 13.7 The process failure, and the guard

Three wrong claims propagated through a design document and **two adversarial review rounds**
without either reviewer catching them — because both reviewers were handed the design and not the
evidence, while the evidence sat in the reference directory.

**Guard, now implemented:** `scripts/review_with_vault.sh` injects vault context into every
reviewer prompt — absolute paths for Kimi (which has a shell but runs in a worktree the vault is
absent from) and pasted note text for Codex (which has no working shell on this box). Verified
end-to-end: both reviewers independently recovered §13.1–13.3 from the vault with citations.

Still open: `MassCalSplit` versus `MassCalCenter`, how the RT bins are placed (equal-count vs
equal-width), and which anchors feed the fit, are all **unread** in the source.

---

## 14. The window sweep, and what it settled — 2026-08-10

Step 1 and its follow-ups produced one clear result, one refuted diagnosis, and
one open question that is larger than the phase that found it.

### 14.1 Shape gating is gone

The mass calibration used to be gated on the SHAPE of its residuals — an
absolute peakedness bar, and a comparison against the m/z-shifted control. Both
are proxies for "are these real fragments"; the residual of the CENTRED
distribution answers it directly. Acceptance is now `after/before <= 0.25` with
a 5 ppm backstop.

The proxy was measurably wrong. On Astral the control produced 98 residuals to
the target's 4,033, so its peakedness came out 12.00 from an edge count of one
to three — 12.00 ± 7 to ±12 against a target statistic of 6.90 ± 0.74. The gate
failed the run on that. Meanwhile the same probe measured a match RATE of 1.344
per target cell against 0.0170 per control cell: a **79× enrichment, 558 sigma**
under the null. Astral's evidence that its matches are real is nine times
stronger than IH1's, and Astral is the run that failed.

**Astral now calibrates**: offset −1.2705 ppm, confirmed independently by the
ID-anchor median (−1.202) and by the historical `min_control_residuals = 400`
experiment (−1.27). Three routes, same number.

### 14.2 The collapse was the WIDTH, and specifically the wrong sigma

Calibrated Astral, narrowed to the model's own 2.89 ppm: pass 1 fell from
**12,211 target groups to 665**. The offset was right and the run still died.

The mechanism is that 2.89 ppm is `3 x 0.96`, and 0.96 ppm is the PROBE's sigma
— measured through a 50 ppm search on a brightness-cut, mobility-gated
population. The run's real per-fragment sigma, from FDR-accepted fragments, is
**2.012 ppm**. On IH1 the same two estimators give **4.19** and **2.087**.

So the probe's sigma differs by a factor of FOUR between two runs whose true
precision agrees to 4%. **One multiple times that sigma cannot mean the same
thing on two instruments.** That is the defect, not narrowing as such.

### 14.3 The sweep

k = 1..10 multiples of each run's own ID-derived sigma, offset and shape applied
in every arm, only the width varying (`-mass_calibration_offset_only` — without
it every arm clamps to the model's width and the sweep measures nothing):

```
    k        1    2    3    4    5    6    7    8    9   10
    IH1    226  665  721  832  817  809  787  842  743  775
    Astral   0  318  617  866 1069 1125 1172 1148 1194 1222
```

IH1 plateaus from k=4 (743-842, ±6% noise). **Astral never turns over.** k=8
maximises the worst case across the two at 93.9%, and is now the default
(`-mass_sigma_multiple`).

Purity moves the other way throughout, measured three ways:
- per-fragment sigma inflates 0.96 → 1.95 ppm (IH1) as the window widens 10x;
- decoy groups outnumber target groups from k=3 on IH1;
- fragment-level purity 99.3% at ±2 ppm falls to 89.1% at ±50 (Astral).

But narrowing buys little of it: ±10 → ±2 on Astral gains 1.9 points of purity
and costs 28% of the evidence. Interference was never the dominant term.

### 14.4 THE OPEN QUESTION, which is bigger than this phase

**On Astral the uncalibrated 50 ppm fallback still beats every calibrated arm
in the sweep: 2,078 identifications against 1,222 at k=10.**

So k=8 is a REGRESSION on Astral against what ships today, and the default is
provisional on that basis. Astral's optimum lies out near 50 ppm (k≈25), which
is what `ChromatogramExtractor.h`'s own sweep already found (IDs 4,290 / 4,499 /
4,969 / 4,382 / 3,967 at 15 / 30 / 50 / 75 / 100 ppm).

Two readings, and they are not distinguishable from identification counts:

1. Astral genuinely needs a wide window — its true fragment component is
   heavy-tailed, and the deconvolution puts 99.9% coverage at ±30 ppm.
2. The extra identifications at 50 ppm are **false**, admitted by a wide window
   and passed by a classifier that cannot tell. Wen et al. 2025 put DIA-NN's
   true precursor FDP above 2.3% at a nominal 1%, and this project has
   explicitly refused identification count as a calibration objective.

**Entrapment settles this and nothing else does.** Until it is run, the window
default is being chosen by the metric this project disowned.
