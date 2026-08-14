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
