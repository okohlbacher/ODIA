# CCS / ion-mobility calibration: adversarial review, and the measurement that settled it

Reviewers: Codex (codex-cli 0.147.0, effort max, two rounds) and Kimi (0.36.1).
Both worked from pasted sources; Kimi ran in a throwaway worktree, which came
back clean. HEAD f7315c7.

**Outcome: agreement on all six questions, one genuine dispute, and the dispute
was settled by a measurement on data already in hand rather than by argument.**

## 1. The physics is CORRECT. Do not touch it.

Both reviewers independently confirmed `include/odia/Library.h`:

- The relation `1/K0 = CCS * sqrt(mu*T) / (18509*z)` is right for a Bruker
  timsTOF's reported 1/K0. Kimi re-derived the constant from Mason-Schamp and
  got **18509.9**.
- **T = 305 K is not an STP double-count.** Reduction to N0 normalises number
  density at 273.15 K; the `T` under the square root is the collision
  temperature and survives the reduction, because `K0 = K*N/N0 ∝ T^-1/2`.
  Substituting 273.15 K would multiply mobility by `sqrt(273.15/305) = 0.946`,
  a 5.4% move in the WRONG DIRECTION -- the data want +1.4 to +2.8%.
- `m_ion = mz * charge` is **correct**. My brief proposed subtracting
  `z * 1.00728`; Codex: that recovers the neutral peptide mass, not the ion's,
  so the proposed "fix" is backwards. The code was already right.

## 2. Field heating is dead as an explanation for the charge dependence

I had proposed (H2) that `T_eff = T_gas + m_gas*v_d^2/(3k_B)` could make the
effective temperature legitimately rise with charge, matching the observed
implied 313.7 / 319.8 / 322.3 K for z = 2/3/4.

**Both reviewers killed it on the same mechanism, reached independently:**

> In a TIMS device the ion is held stationary against the gas flow, so the
> ion-neutral relative velocity IS the gas velocity -- common to every ion at a
> given position, regardless of charge. Unlike a drift tube, where
> `v_d = K*E`, charge does not enter.

Both note the MAGNITUDE is fine (v_g ~ 100-200 m/s gives dT ~ 11-45 K, and the
common offset from 305 K to ~317 K is plausibly real field heating). It is only
the z-DEPENDENCE that field heating cannot produce; the reduced-mass change from
z=2 to z=4 moves the heating term by ~1%, not the 8.6 K required.

So: **the shared offset may be physics; the charge dependence is not.** Codex
adds the correct caution that "proves AlphaPeptDeep bias" is too strong -- a
charge-dependent instrument calibration, anchor selection, or population
confounding could also produce it.

## 3. The dispute: is the fitted slope mostly a statistical artifact?

**Kimi's finding (which Codex round 1 had listed but not developed).** The
calibration regresses `delta = observed - library` ON `library`. Since `library`
appears on both axes, any library error induces a negative slope by
construction. With `L = T + eps`, `O = T + eta`:

```
    E[slope] = -Var(eps) / Var(L)
```

Kimi estimated this at about **-0.11** from a 2.8% library error, i.e. the same
order as every measured slope (-0.0245 z=2, -0.0423 pooled, -0.113 z=3), and
concluded the slopes are largely fictitious.

**Codex's rebuttal.** The formula needs classical measurement error
(`Cov(T,eps) = 0`), which "2.8% prediction error" does not establish -- for a
calibrated predictor `T = L + u` with `Cov(L,u) = 0` the induced slope is ZERO.
And genuine scale and attenuation are confounded: with `O = a + cT + eta`,

```
    1 + beta = c * lambda,    lambda = Var(T)/(Var(T)+Var(eps))
```

so the arithmetic was defensible but the empirical input was not, and Kimi's
"same order as every slope" was **unsupported**.

## 4. The measurement that settled it

Codex named the decisive test: estimate the artifact directly, per charge, by
regressing `(measured - library)` on `library` over the paired predicted/measured
set, with no classical-error assumption. Run on 33,749 precursors paired between
our library and DIA-NN's S08 empirical library:

| z | n | **a_q (measured artifact)** | mean relative bias | sd(library 1/K0) |
|---|---:|---:|---:|---:|
| 2 | 23,650 | **-0.0768** | +1.68% | 0.1245 |
| 3 | 8,329 | **-0.0777** | +2.59% | 0.0948 |
| 4 | 971 | **-0.0581** | +2.97% | 0.0890 |

**Kimi was right in substance and wrong in magnitude.** The artifact is real and
is the same order as the fitted slopes -- but it is **-0.077, not -0.11**, and
Codex was right that the 2.8% input did not establish it.

Applying Codex's separation `g_q = (1+d_q)/(1+a_q)` to the run's fitted slopes:

| z | fitted slope d_q | artifact a_q | **genuine scale g_q - 1** |
|---|---:|---:|---:|
| 2 | -0.0245 | -0.0768 | **+5.67%** |
| 3 | +0.0528 | -0.0777 | **+14.15%** |

**This resolves the sign paradox.** The apparent disagreement between charge 2
(-2.5%) and charge 3 (+5.3%) was the shrinkage term. Once removed, both charges
show a POSITIVE genuine scale error, consistent in sign with the per-charge mean
biases (+1.68/+2.59/+2.97%) and with the global coefficient fit. Nothing
contradictory remains.

Note the magnitudes are slopes, not mean ratios: a +14% slope on charge 3 is
compatible with a +2.6% mean bias, because the slope describes how the bias
varies across the mobility range. They are different quantities and must not be
compared directly -- Codex round 1 made exactly this point.

## 5. What is actually broken

**a) The "coefficient error" gloss in the log is wrong. Remove it.**
The line `1/K0-linear -0.0245 per 1/K0 about 0.991, i.e. a -2.5% error in the
CCS->1/K0 coefficient` is invalid. Codex's arithmetic: if -0.0245 were a pure
multiplier, the correction at the pivot would have to be `-0.0245*0.991 =
-0.0243`, not the logged **+0.017**. Same for charge 3: a +5.3% scale implies
+0.046 at 0.871, not the logged +0.010. The gloss is not a measurement.

**b) The stated physical justification for pooling the slope is FALSE.**
All three reviews agree, and the measurement proves it. The header claims:

> "Pooling the SLOPE is physically justified ... The slope is a relative scale
> error in the CCS->1/K0 conversion coefficient, and that coefficient is a
> property of the conversion, shared by every charge."

The fitted slope is `c*lambda_z`, and `lambda_z` is the charge's own reliability
ratio. Measured, `a_q` spans **-0.058 to -0.078, a 28% spread across charge** --
so the shrinkage component is demonstrably NOT shared. Pooling may still be
useful as statistical regularisation, but the physics argument for it is wrong
and must be struck.

This also explains the pooled tier's signature (out-of-fold MSE improves, IDs
fall): it applies one charge's shrinkage fraction to another.

**c) `min_anchors_per_charge = 120` is not defensible as a single gate.**
Both agree. For a constant-only correction at 38 charge-4 anchors, SE = 0.0013
against a ~0.018 correction. But 38 anchors cannot support a slope plus 8 m/z
bins, and both warn a naive `SE/|correction|` rule has a selection problem -- it
preferentially fires on accidentally large estimates. The agreed shape is a
TIERED, model-specific gate: constant-only from n ~ 30 with per-charge
out-of-fold improvement as the acceptance test; keep 120 for shaped fits.

**d) The 20-ID result is underpowered.** 20/1232 = 1.6% from one run. Kimi gives
the test: McNemar on discordant pairs, plus fold-seed replicates (the fold
assignment is a precursor hash -- perturbing the salt is free). If >= 12 of the
20 lost IDs are charge 3, the dilution mechanism is confirmed.

**e) Extrapolation beyond anchor support.** Codex: a slope or m/z shape fitted on
198 or 38 anchors is applied to a 5-million-precursor population. Needs explicit
clamping or shrinkage outside the anchors' range.

## 6. What NOT to do

- Do not replace 305 K with a fitted or charge-dependent temperature. A fitted
  "temperature" from predicted CCS is not a thermometer -- it is non-identifiable
  from a multiplicative CCS-model error.
- Do not change 18509 or the reduced-mass definition.
- Do not discard the slope correction. Even as partly shrinkage, it is the best
  linear predictor of observed from library on the same population; out-of-fold
  MSE rewards it legitimately. Only the INTERPRETATION and the POOLING are wrong.

## 7. Ranked actions

| | action | basis |
|---|---|---|
| 1 | Strike the "i.e. an s% coefficient error" gloss from the log and the header | both reviewers, plus the pivot arithmetic |
| 2 | Strike the physical justification for pooling; re-document as regularisation | both, plus a_q's 28% charge spread |
| 3 | Tiered gate: constant-only tier at n ~ 30, gated on per-charge out-of-fold gain; 120 kept for shaped fits | both |
| 4 | Clamp corrections outside anchor support | Codex |
| 5 | McNemar + fold-seed replicates before believing the 20-ID divergence | both |
| 6 | Report `a_q` alongside the fitted slope so the two are never conflated again | this measurement |

Nothing here changes `Library.h`. The conversion was right the whole time; what
was wrong was the story told about the correction fitted on top of it.
