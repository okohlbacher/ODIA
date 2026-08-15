# Linear anchor, then fine-tune the nonlinearity — measured

Status: MEASURED 2026-08-14, both instruments. Answers the question directly:
anchor linearly and fix nonlinearity in fine-tuning, or fit a monotone curve
first and fine-tune only its residual?

## Answer: LINEAR FIRST. The nonlinearity is real but a CiRT seed cannot support it.

### 1. The relation IS nonlinear

Linear fit on all DIA-NN precursors with a library iRT, residual median per iRT
decile:

| | decile medians (s) | curvature |
|---|---|---|
| Astral (n=10,040) | +13 +4 −3 −6 −3 −4 −10 −31 −57 **+78** | **135 s** |
| S08 (n=32,684) | +15 +10 +4 −7 −22 −38 −50 −55 +0 **+120** | **175 s** |

Flat through the middle, sagging late, then a sharp positive jump in the last
decile — gradient end-effects. A straight line mis-places the extremes by well
over a ±60 s window.

### 2. But fitting that curve from CiRT anchors makes things WORSE

From the blind CiRT search only (24 RANSAC inliers on Astral, 47 on S08),
evaluated externally on DIA-NN's confident set:

| | ±30 s | ±60 s | median |
|---|---|---|---|
| Astral **linear** | **37.6%** | **65.5%** | **42.0 s** |
| Astral monotone PCHIP | 32.8% | 58.9% | 48.5 s |
| S08 **linear** | **56.2%** | **81.2%** | **25.6 s** |
| S08 monotone PCHIP | 51.9% | 76.7% | 28.7 s |

**Monotone loses 5–7 pp on both instruments.** ~30 anchors cannot constrain a
curve; the fit chases anchor noise and generalises worse. This is the same
overfitting the vault records elsewhere (akima across 8 noisy block centres,
`doc/15` §9).

### 3. The nonlinearity is worth having — but only from thousands of IDs

Fitted on half the DIA-NN IDs, evaluated on the held-out half:

| | best LINEAR | best MONOTONE | nonlinearity worth |
|---|---|---|---|
| Astral (held-out 5,020) | 71.8% @±60 s | 76.8% | **+5.1 pp** |
| S08 (held-out 16,342) | 74.8% @±60 s | **87.4%** | **+12.6 pp** |

So the curve is worth 5–13 points at ±60 s — real, and on S08 substantial — and
it is only reachable once many identifications exist. That is precisely
fine-tuning's job.

**And the CiRT linear seed is already close to the linear ceiling**: 65.5% against
71.8% on Astral, 81.2% against 74.8% on S08 (the S08 comparison is approximate —
the seed was scored on all points, the ceiling on a held-out half). The seed is
near-optimal *for a line*; the residual gap is exactly the nonlinear part.

### 4. Computational effort is not a consideration

| | fit time |
|---|---|
| linear (lstsq) | **0.09 – 0.26 ms** |
| monotone PCHIP, ~30 anchors | 0.5 – 1,222 ms |
| monotone PCHIP, binned from thousands | **15 – 43 ms** |

Against a 34-minute blind extraction, every option is free. **Compute does not
decide this; generalisation does.** Do not let a cost argument enter here.

## The resulting architecture

```
  1. CiRT blind search        -> ~30-47 robust anchors
  2. RANSAC LINEAR            -> near the linear ceiling, no curve fitted
  3. extract in a window      -> identifications now exist
  4. MONOTONE fine-tune       -> on the RESIDUAL of the line, binned medians +
     from the IDs                PCHIP, monotone enforced. This is where the
                                +5 to +13 pp lives.
  5. iterate                  -> stop on parameter stability, never ID count
```

The user's second option — monotone first, fine-tune the residual — is refuted by
§2: there are not enough seed anchors to fit a monotone curve at all, and
attempting it costs 5–7 pp before fine-tuning even starts.

## Caveats

- Fine-tuning at step 4 anchors on the *detected* subset, which is the regime
  where four anchor statistics pointed the wrong way during the original RT work.
  It must be validated externally (held-out DIA-NN coverage), never on its own
  anchors.
- The last-decile jump (+78 s Astral, +120 s S08) is where a monotone fit will
  help most and where extrapolation is most dangerous. Anchor coverage in the
  final decile must be reported, not assumed.
- Both instruments agree on the qualitative shape, which is reassuring, but the
  magnitude differs 2.5× (5.1 vs 12.6 pp) — so the value of fine-tuning is
  instrument-dependent and must be measured per run, not assumed from Astral.

---

## 5. REVIEW ROUND 9 — the architecture survives, most of my reasoning does not

Codex (effort max) and kimi, both with vault access; kimi ran its own measurements
(`bench/rt_step4_sim.py`). Vibe: 0 bytes, sixth time.

**Both accept the architecture** — linear CiRT seed, monotone fine-tune from IDs —
and kimi notes DIA-NN's own pipeline is an existence proof of exactly that
division of labour. Almost everything else below is a correction.

### 5.1 My central experiment was structurally unfair (codex)

The PCHIP in §2 was fitted on **RANSAC inliers — anchors selected for agreeing
with a line.** The curve was trained on a line-conditioned subset and could not
have won. §2 therefore does **not** refute monotone-first as a model class; it
refutes *"an interpolating PCHIP fitted de novo on these particular
line-selected anchors"*. Untested: a 3–4 degree-of-freedom constrained curve, a
curve shrunk toward a line, or an affine alignment of a template learned on other
runs.

### 5.2 The nonlinearity is mostly OUR PREDICTOR, not the gradient (kimi, measured)

§1 measured curvature on iRT values from `odia_lib_parity.parquet` — ODIA's own
predictor. Redone on DIA-NN's independent iRT axis, same reports, same filter:

| axis | Astral curvature | S08 curvature |
|---|---|---|
| ODIA parity-lib iRT (§1) | 135 s | **175 s** |
| DIA-NN report `iRT` | 90 s | **22 s** |
| DIA-NN `Predicted.iRT` | 88 s | **12 s** |

**On S08 the "sharp +120 s last-decile jump" collapses to 12–22 s — inside a
±60 s window — on an independent axis.** It is largely ODIA-predictor error at
extreme iRT rank, not chromatography. Consequences:

- **"Both instruments agree on the shape" was not independent confirmation** —
  both were measured on the *same predictor's* axis. The agreement is evidence
  about our predictor, not about two chromatographies.
- **The fine-tuning prize is confounded with predictor quality.** Much of S08's
  +12.6 pp is fixing our own iRT prediction. Swap in peptdeep (already in the
  tree, and it beat everything in the frozen work) and the prize shrinks. The
  caveat should read *predictor-axis-dependent, re-measure when the library
  generator changes.*
- On Astral 88–90 s survives, so a real chromatographic component exists — but
  its shape is a **frown**, not my "flat middle, sag late, jump at the end". My
  physical story does not survive contact with an independent axis.

### 5.3 Step 4 truncation, simulated (kimi)

Seed = the actual frozen-gate lines; "detected" = DIA-NN points within ±60 s of
the seed, i.e. what a windowed extraction can see; fine-tune on a parity half;
evaluated externally.

| | seed | fine-tune on **detected** | unconditioned, same n | half of all |
|---|---|---|---|---|
| Astral ±60 s | 65.5% | **71.3%** | 77.4% | 77.8% |
| S08 ±60 s | 81.2% | **84.1%** | 87.5% | 87.5% |

Step 4 **helps** (+5.8 / +2.9 pp) but loses **~6 pp on Astral and ~3.4 pp on S08**
to truncation. Codex's formalisation: detection admits an apex only when
|residual| ≤ 60, so step 4 estimates `median(r | x, detected)` while it needs
`median(r | x)`. The tail needing the largest correction supplies the least
correct evidence.

### 5.4 Arithmetic errors in §3 and §4 (codex)

- **S08's "linear ceiling" (74.8%) is BELOW the CiRT seed it supposedly bounds
  (81.2%).** A ceiling cannot be lower than its candidate. I wrote this off as
  "different evaluation subsets"; that invalidates the comparison rather than
  rescuing it.
- **Astral's seed-to-ceiling gap is 6.3 pp — larger than the 5.1 pp nonlinear
  gain.** Calling the seed "near-optimal" while calling the smaller difference
  architecturally decisive is inconsistent.
- **"5–7 pp on both instruments" is wrong for S08**, which is 4.3 / 4.5 pp.
- **OLS is not a ±60-coverage ceiling** — it minimises squared residuals while
  coverage is thresholded consensus. A real ceiling needs a frozen objective and
  identical evaluation rows.
- **A 12-point discrepancy is unexplained**: the frozen vault note reports 88.89%
  at ±60 s for Astral library iRT; my "best monotone" is 76.8%. And my 71.8%
  coincidentally equals the frozen 71.88% *±30* figure. Both demand an artifact
  audit before any of these numbers are cited again.

### 5.5 The gate claims were overstated (codex)

P5 and P6 are unrun and P7 was not reported in the S08/Astral results. **Both
instruments have partial P1–P4 passes, not full PASSes.** I reported "PASS".
Also, P4's permutation should block by stripped sequence — permuting precursor
rows lets charge duplicates inflate both consensus and apparent n.

### 5.6 A mathematical fix (codex)

The measured residual **decreases then jumps** — it is not monotone. Applying
PAVA/monotone regression to the *residual* cannot represent that shape. Fit
`f(x) = a + bx + δ(x)` subject to `f'(x) > 0` on the **total** map, shrinking the
curvature of δ; do not constrain δ itself.

### 5.7 What stands

The architecture, the CiRT seed, and the direction of §2's result for the
*specific* estimator tested. Everything quantitative needs re-measuring on an
independent iRT axis with frozen artifacts and identical evaluation rows.

---

## 6. §1 RE-MEASURED on independent axes — kimi confirmed, quantified

Reproduced 2026-08-15 (`bench/axis_recheck.py`), same reports, same q≤0.01 filter,
identical precursors, only the iRT axis changed:

| axis | Astral | S08 |
|---|---|---|
| ODIA parity-lib iRT (**what §1 used**) | 135.2 s | **174.6 s** |
| DIA-NN report `iRT` | 90.8 s | **20.8 s** |
| DIA-NN `Predicted.iRT` | 88.9 s | **11.8 s** |

**S08's nonlinearity is ~90% our own predictor's error.** 174.6 s → 11.8 s on an
independent axis is inside a ±60 s window, i.e. not worth correcting per-run at
all. §1's "sharp +120 s last-decile jump" is our predictor failing at extreme iRT
rank.

**Astral retains a real ~89 s chromatographic nonlinearity** — but its shape is a
**frown** (−4 −43 −32 −20 +4 +18 +35 +48 +40 −35), not §1's "flat middle, sag
late, jump at the end". The physical story in §1 is wrong on both instruments.

**Consequence for priorities:** on S08 the fine-tuning prize is mostly recoverable
by **improving the iRT predictor**, not by per-run refinement — and peptdeep is
already in the tree and beat everything in the frozen work. That is a cheaper and
more portable lever than the per-run monotone fit, and it was invisible while the
measurement used our own predictor as its x-axis.

§5.2's caveat stands and is now measured: **predictor-axis-dependent, re-measure
whenever the library generator changes.**
