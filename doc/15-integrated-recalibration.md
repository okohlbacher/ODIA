# Integrated recalibration: MS1 + MS2 mass, then RT, then both to convergence

Status: DESIGN, not implemented. Written 2026-08-10, after RT recalibration was
frozen (`rt-calibration-frozen-2026-08-10`).

This document is the plan of record for replacing today's single-shot mass
calibration with a joint m/z + RT recalibration that iterates. It is written
*before* measurement, in the same discipline as `14-irt-acceptance.md`, because
every time in this project a calibration statistic was chosen after seeing it,
it was chosen wrong.

---

## 0. What is wrong with today

| | today |
|---|---|
| MS2 mass model | one global `intercept + a·log(m/z) + b·m/z`, fitted once |
| MS1 mass model | **none** — `Ms1Traces` borrows `fragment_ppm` and applies no offset |
| anchors | `MassCalibration::collect`'s own probe: brightest co-occurring cells over the whole run, chosen without reference to any identification |
| gate | peakedness ratio ≥ 3 against an m/z-shifted control; on failure the model is discarded and the window opens to 50 ppm |
| when | before pass 1, so no RT map exists; optionally re-measured once after the map under `-mass_calibration_remeasure` |
| iteration | none |

Three consequences, each observed:

1. **The probe fits a sample its own search ruined.** On S08 the gate fails and
   reports "a mostly-noise sample" (`MassRecalibration.h` records the
   measurement: an apparent −8.44 ppm offset was −4.98 ppm at a control 300 s
   away, so the real offset is ≈ −3.5 ppm and the rest was the ±50 ppm search).
2. **MS1 is uncalibrated and shares a tolerance derived from fragments.** There
   is no reason precursor and fragment error share a model; the literature and
   every other engine fit them separately.
3. **The one attempt to size the window from identifications was harmful**
   (`-mass_width_from_ids apply`: Astral 4,969 → 2,422), because FDR-accepted
   groups are a *selected* sample, not a truncated one. Their residual spread is
   biased narrow by exactly the window that admitted them.

`MassRecalibration` — RT-blocked, linear in log m/z, MAD-trimmed, bounded — is
already written and unit-tested (`test/tools/odia_mass_recalibration.cpp`) and
is **not wired to anything**, because nothing produces the per-fragment
residuals it consumes.

---

## 1. The shape of the new phase

```
  bootstrap pass          wide m/z (50 ppm), wide RT (library linear iRT)
        |                 -pass1_precursors sample, collect_mass_residuals ON
        v
  score + FDR             existing semi-supervised path
        |
        v
  ROUND r = 1..N:
        |
        |-- harvest        per-fragment MS2 residuals + per-precursor MS1
        |                  residuals, at the apex of accepted groups,
        |                  PLUS the same from decoy groups (the null)
        |
        |-- fit mass       MS2 model and MS1 model, INDEPENDENTLY
        |                  deconvolve true component against the decoy null
        |
        |-- size windows   from the deconvolved true component at a coverage
        |                  target, not from k*SD of the accepted sample
        |
        |-- fit RT         existing frozen path: isotonic/akima map, then
        |                  per-run refinement to 0.1 s convergence
        |
        |-- re-extract     narrower m/z, narrower RT, corrected masses
        |
        |-- REPORT         one row per round (section 5)
        |
        `-- converge?      on DELTAS of sigma_mz(MS1), sigma_mz(MS2), sigma_RT
                           over a FROZEN evaluation set
```

The user's ordering constraint — mass recalibration *before* RT recalibration —
is honoured: within each round, mass is fitted first and the RT map is fitted
against mass-corrected extraction. Round 1's mass fit necessarily runs on a
library-iRT bootstrap, because confident fragments require a search and a search
requires some RT estimate. That is a bootstrap, not a circularity, provided the
bootstrap window is wide enough not to censor (section 4).

---

## 2. Where the residuals come from

### 2.1 MS2 — the hook already exists

`ChromatogramExtractor` already keeps two `[row][cycle]` planes under
`collect_mass_residuals`: `ppm_num` (Σ intensity·ppm) and `ppm_den` (Σ
intensity). `row` **is the fragment**, so the per-fragment axis we need is
already there and is being thrown away.

`PeakGroupScorer.cpp:845-896` walks those planes at the apex and collapses them
to two scalars — `mass_ppm` (intensity-weighted group mean) and
`mass_ppm_spread` (1.4826 × MAD across fragments). The collapse destroys the m/z
axis, which is the axis the recalibration function is defined on.

**Change:** at the same site, when a harvest sink is attached, also emit one
`MassResidual{mz = row's theoretical m/z, rt = apex RT, ppm, intensity}` per
contributing fragment. No new decode, no new pass, no new memory in the hot
loop — the values are in registers at that point.

Cost estimate: ~12 fragments × ~15 k accepted groups × 16 B ≈ 2.9 MB. Negligible
against the 10.3 GiB decode floor.

### 2.2 MS1 — a new path, deliberately symmetric

`Ms1Traces::build` (`src/extract/Ms1Traces.cpp:82,87`) computes
`tol = m * fragment_ppm * 1e-6` and tests `|it->mz - m| > tol` — the deviation
`it->mz - m` is in scope and discarded, exactly as the MS2 deviation was before
`-collect_mass_residuals`.

**Change:** retain an intensity-weighted ppm per (precursor, cycle), mirroring
the MS2 planes, and emit `MassResidual` at the accepted group's apex with
`mz = precursor m/z`. One residual per accepted precursor per round, so ~15 k
records — smaller than MS2 by the fragment multiplicity.

This is the piece that does not exist today in any form. It is the single
largest structural gap against DIA-NN.

### 2.3 The null — decoys, not a second probe

`MassCalibration::collect` builds its own m/z-shifted control. That control is
the right *idea* and the wrong *source*: it samples cells the search never
scored.

We already run target-decoy. **Decoy peak groups' fragments are the null we
need**, at zero extra cost: they were extracted through the same window, scored
by the same classifier, and accepted by nothing. Their ppm distribution is the
random-match distribution *shaped by the window*, which is precisely the
contaminating component in the target distribution.

Set `MassResidual::decoy = true` for them. `MassRecalibration::fit` already
ignores flagged residuals when fitting; we additionally *use* them, for
deconvolution (section 3.3).

---

## 3. The recalibration functions

### 3.1 MS2

Reuse `MassRecalibration`'s form, with one upgrade:

```
    ppm(m/z, rt) = intercept(rt) + slope(rt) * log(m/z / m_ref)
```

- **linear in log m/z**, not m/z: the existing measurement reports the shape
  that way ("3.322 ppm per e-fold in m/z"), and a linear-in-m/z fit over
  200–1800 Th is dominated by the top of the range.
- **RT dependence**: today's `blocks = 8` equal-count blocks are *piecewise
  constant*, which puts a step discontinuity at each block edge. A mass
  spectrometer does not step. **Fit the blocks as now, then interpolate
  `intercept` and `slope` across block centres with the same akima
  interpolation the RT map uses** (`TransformationModelInterpolated`,
  `interpolation_type = akima`), clamped to the end blocks outside the range.
  This reuses a dependency already in the build and keeps one interpolation
  story across the tool.
- **Rails kept**: `max_correction_ppm = 20`, `max_slope_ppm = 15`,
  `trim_mads = 4`, `min_per_block = 150`, `min_for_slope = 400`. A block below
  `min_per_block` inherits the global fit rather than inventing a local one.

### 3.2 MS1 — an independent model of the same form

Fitted from precursor residuals only, with its own blocks, its own rails, and
its own reported σ. **It must not be constrained to agree with MS2.** If the
two models come out the same, that is a finding; forcing it is an assumption.

Practical difference to expect: far fewer anchors (one per precursor, not
twelve), so `min_per_block` must scale down — propose `blocks = 4`,
`min_per_block = 100`, `min_for_slope = 300`, and fall back to a single global
block when the run does not support four.

### 3.3 Sizing the window — the part that was wrong before

`-mass_width_from_ids apply` failed because it took the spread of accepted
groups at face value. The fix has two parts.

**(a) Deconvolve.** The target residual histogram is a mixture:

```
    f_target(ppm) = pi * f_true(ppm) + (1 - pi) * f_null(ppm)
```

with `f_null` measured directly from decoy fragments (same window, same
classifier) rather than assumed. Estimate `pi` and the scale of `f_true` by
fitting the mixture on the *uncensored* bootstrap window. A robust scale
(MAD-based, or Huber) on the deconvolved true component, not an SD on the raw
sample.

**(b) Size to a coverage target, not to k·σ.** The objective that survived the
RT work was *coverage at window*, and it should be the objective here too:
choose the window as the quantile of `f_true` that admits a target fraction
(propose 99%) of true fragments, floored at a minimum (propose 5 ppm) so a
pathologically confident round cannot collapse the window.

**(c) The "may only narrow" rail is retained but bounded per round.** Limit any
single round's narrowing to a factor (propose 2×) so an early biased round
cannot slam the window shut before the anchor set matures. This directly
addresses the observed failure mode of `apply`.

---

## 4. Iteration and convergence

### 4.1 What is refit each round

| refit each round | held fixed |
|---|---|
| MS1 mass model | the library |
| MS2 mass model | the classifier's feature set |
| m/z windows (MS1, MS2) | the FDR procedure |
| RT map + per-run refinement | the frozen evaluation set (4.3) |
| RT window | |
| classifier | |

### 4.2 The convergence criterion

By the precedent set in `14-irt-acceptance.md`, convergence is tested on the
**delta**, not the absolute:

```
    stop when   |Δσ_mz(MS2)| < 0.1 ppm
          and   |Δσ_mz(MS1)| < 0.1 ppm
          and   |Δσ_RT|      < 0.1 s
          or    round == cap (propose 8)
```

with a relative floor (1%) as the alternative trigger, mirroring
`-rt_converge_rel`. Rationale for a *delta* criterion is the one already
learned: an absolute target is a statement about the sample, and every absolute
statistic taken over our own anchors described the precursors we already found.

### 4.3 The frozen evaluation set — non-negotiable

The "3.3× narrowing" that turned out to be a denominator artefact happened
because the anchor population grew between rounds (1,183 → 1,408). The same trap
is available here and is *worse*, because narrowing the m/z window changes which
precursors are accepted.

**Freeze, at the end of round 1, the set of precursors over which σ is
reported.** The fitting set may grow; the reporting set may not. Any round that
loses members of the frozen set reports them as failures (infinite residual for
coverage purposes), not as absences. Without this, every number in section 5 is
uninterpretable.

### 4.4 Failure modes to instrument, not to hope about

- **Runaway narrowing** — bounded by 4.3's frozen set, the per-round 2× cap, and
  the 5 ppm floor.
- **Oscillation** — detect a σ that increases twice consecutively; revert to the
  best round's models and stop. Report it loudly.
- **Early bias lock-in** — round 1's accepted set is the most biased. Propose
  fitting round 1's mass model at a *deliberately lenient* FDR (5%, matching the
  existing `-train_fdr_initial` philosophy) to widen the anchor base, tightening
  to 1% from round 2.
- **A block fitting interference** — the rails, plus: report per-block anchor
  counts every round so a block that inherits the global fit is visible rather
  than silent.

---

## 5. The per-iteration report

One row per round, to the log and to a TSV under `-out_recalibration`:

```
round  n_acc  n_frag  ms2_off  ms2_slope  ms2_sig  ms2_win  ms1_off  ms1_sig  ms1_win  rt_sig  rt_win  cov_mz  cov_rt  ids  d_sig_mz2  d_sig_mz1  d_sig_rt
```

- `*_sig` are **deconvolved, robust** scales on the **frozen** set (4.3).
- `cov_*` are coverage-at-window on the frozen set — the objective, stated
  first-class rather than inferred from σ.
- `d_*` are the convergence quantities of 4.2, so the stopping decision is
  auditable from the report alone.
- `ids` is reported but is **not** the objective and must not be the stopping
  criterion. The stated goal is clean extraction and converged calibration;
  identification counts are a later concern (standing decision from the RT work).

Additionally, dump the round-1 and final residual histograms (target and decoy,
MS1 and MS2) so the mixture fit of 3.3(a) can be checked by eye rather than
trusted.

---

## 6. Implementation order

1. **MS2 per-fragment harvest** — emit `MassResidual` at `PeakGroupScorer.cpp:845-896`.
   Wire `MassRecalibration` to it. Single round, no iteration. *Measure: does the
   RT-blocked model beat today's global model on held-out fragments?*
2. **Decoy null + deconvolution** — 3.3(a). *Measure: is the deconvolved σ
   materially different from the raw σ? If not, the mixture model is not earning
   its complexity and should be dropped.*
3. **Coverage-based window sizing** — 3.3(b,c), replacing `MassWidth`'s
   `apply`. *Measure: does it recover the 4,969 that `apply` destroyed?*
4. **MS1 residual collection** — 2.2. *Measure: is the MS1 error the same
   function as MS2? This is a real question with a real answer.*
5. **MS1 model** — 3.2.
6. **The loop** — section 4, with the frozen evaluation set built first, not
   retrofitted.
7. **The report** — section 5. Written before step 6 produces numbers.

Steps 1–3 are useful independently and can ship without the loop. Step 6 is
where the risk is and is deliberately last.

---

## 7. What would falsify this design

Stated now, so it cannot be rationalised later:

- If the RT-blocked MS2 model does not beat the global model on held-out
  fragments (step 1), **the RT dependence is not real in this data** and blocks
  should be dropped to one.
- If the deconvolved σ tracks the raw σ within 10% (step 2), the decoy null is
  not needed and the mixture fit should be deleted.
- If the MS1 model's fitted offset and slope agree with MS2's within their
  confidence intervals (step 4), a single shared model is correct and the second
  model is waste.
- If σ_mz converges in one round (step 6), the loop is unnecessary and the phase
  is a single pass with a better anchor source — which would still be a win over
  today, and should be shipped as such.

The design is written to be *reducible*: each of these outcomes deletes code
rather than requiring more.
