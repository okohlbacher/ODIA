# Integrated recalibration: MS1 + MS2 mass, then RT

Status: DESIGN v3, not implemented. Written 2026-08-10, after RT recalibration
was frozen (`rt-calibration-frozen-2026-08-10`).

v1 was reviewed adversarially by Kimi 0.34.0, Codex 0.147.0 and Vibe 2.24.0; all
three rejected it. v2 answered them. v3 then folded in a literature and
source-code survey of how the field actually does this (§11), which overturned
the *architecture* rather than the details: the iteration is replaced by an
additive model with RT, intensity and ion mobility as covariates.

§9 records what was cut and why; the deletions are the useful part.

Written *before* measurement, in the discipline of `14-irt-acceptance.md`.

**The headline change in v3:** on timsTOF, the ion-mobility term explains **52%
of modelled mass-error variance** (Prianichnikov et al., MCP 2020, 19(6)
1058-1069). S08 is timsTOF. ODIA's mass model has **no IM term**. That is almost
certainly the largest single miss in this whole area, and it is not something any
amount of iteration would have found.

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

`MassCalibration.h:104-136` records, for S08:

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
| `f2(RT)` | not fitted | **measured absent on S08**: t = 0.26 over 1,600 s, confirmed over 124 M hits |
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

**`f2(RT)` is fitted but expected to be dropped by its own t-test** on S08
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
no new pass. ~12 fragments × ~15 k groups × 16 B ≈ 2.9 MB.

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
validates against (the S08 gate is 4.34 target vs 2.98 control). It preserves the
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
  ODIA's own S08 numbers show why — the shape correction removed 68% of
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
defensible figure from the repo: the probe touched **11.9% of S08's spectra in
84.6 s with only 3,000 sampled precursors**. A full pass over the whole library
is materially more than 84.6/0.119 ≈ 12 min. Eight of those is around
1.5 h per file, on top of a ~10.3 GiB decode floor paid each time.

**It is also unnecessary.** Fitting the mass model needs only
`(mz, rt, ppm, intensity)` per fragment. Once those are retained from **one wide
extraction**, evaluating any *narrower* candidate window is a **filter over
retained residuals**, not a re-extraction. So:

```
  ONE wide extraction  (50 ppm, library iRT, collect residuals)
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
(4) then smallest MS1 window. σ is a *diagnostic*, not a tie-breaker — the S08
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
   this document: 52% of modelled variance on timsTOF, and S08 is timsTOF.
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
   S08 this should drop (t = 0.26). *Measure: does it survive on Astral?*
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
| RT-blocked mass model as the default | RT drift measured at t = 0.26 on S08 and confirmed absent over 124 M hits. Fitting noise. Now opt-in behind a t-test. |
| Akima across 8 block centres | Overfits 8 noisy points; and would treat *inherited* blocks as measured control points, manufacturing end curvature from fallback policy. |
| Decoy fragments as `f_null` | Wrong m/z sampling and wrong contamination structure; makes `π` and `f_true` non-identifiable. Use the m/z-shifted control. |
| 99% coverage as the window objective | Coverage is monotone in width and has no noise term; its optimum is "as wide as possible". Now a constraint, paired with control-density. |
| Re-extract every round, cap 8 | ~12+ min per pass on S08 and a 10.3 GiB decode floor each time. Now one wide extraction + one confirmatory pass. |
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
