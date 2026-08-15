# Integrating the previous project: what to port, and the stepped plan

Written 2026-08-15 from the predecessor's handoff (`odia-v0.4.0`, commit
`32c9dcd`, measured 2026-08-07), review round 11 (codex effort max, kimi), and
this project's own measurements. Vault:
`50-Benchmarks/Previous ODIA vs current ODIA.md`.

## 1. Where the two projects actually stand

| | previous (v0.4.0) | current |
|---|---|---|
| precursors @ nominal 1% | **7,008** | 3,625 |
| fraction of its DIA-NN | **75.7%** | 29.5% |
| **entrapment FDP** | **not measured** | **5.98%** |
| wall | **22:33** (224 thr) | 2:33:35 (48 thr) |
| peak RSS | **82.4 GB** | 212 GB |
| OpenMS | 20-file vendored fork | **read-only, unmodified** |

**The previous attempt was scientifically ahead.** Not a clean comparison —
different DIA-NN versions, libraries and core counts — but 75.7% against 29.5%
is too large to attribute to that.

**Both headline counts are nominal.** The predecessor's own closing caveat is the
strongest statement of the principle in either project: *no DIA search tool
consistently controls FDR at the peptide level (Wen et al., Nat Methods 22:1454,
2025), with DIA-NN's true precursor FDP above 2.3% at a nominal 1%. Gate any
change on entrapment FDP, not on identification counts.*

## 2. Corrections to this project's claims (round 11)

- **"3,625 at 1% FDR" was wrong wording.** It is 3,625 at nominal q ≤ 0.01 with a
  measured FDP of **5.98%** — six times nominal, and level with OpenSWATH's 5.93%.
- **The tail statistic mixed units.** `q` operates on the best row per precursor,
  so the same-unit progression is **27 → 490 precursors**, not 56 → 1,098 rows.
  Recomputed floors reproduce the observed values exactly (0.03654 / 0.00204
  against 0.0365 / 0.00199).
- **"Clean seed" is unearned.** Targets above every decoy are a *target-enriched
  tail*, not known true positives — non-exchangeable decoys, leakage and lucky
  absent targets all produce it.
- **The ~100 threshold does not certify anything.** It says when this estimator
  *can emit* q ≤ 0.01; it is not a calibration or a confidence bound. Zero
  failures in 100 trials still carries a ~3% upper 95% bound.
- **`targets-above-the-top-decoy` is tied to the floor algebraically**, so it is a
  cheap reachability diagnostic and **not** independent evidence that an
  intervention improved anything.

## 3. What to port, in order of value per unit of risk

**P1 — tcmalloc via `LD_PRELOAD`.** Measured there: **−44% peak RSS for +21% wall,
identical identifications.** Live data never exceeded ~33 GB against an 82.4 GB
peak — 88% allocator debris. We peak at **212 GB**. One environment variable, no
code, and it is verifiable by the invariant that IDs must not move.

**P2 — the phase-timing breakdown and its Amdahl analysis.** The predecessor
found 367.8 s single-threaded, *exceeding DIA-NN's entire 331 s runtime*, so no
amount of parallelism could win. We are 8.6× slower and **have never measured
where our time goes.** Same instrumentation, same question.

**P3 — audit our code for the three defects it names**, because two are bug
*classes* we have already hit:
- `min_child_rows` not enforced per child → leaves of ~4 rows, "an overfitting
  surface in exactly the score tail that sets the FDR threshold". **Check our
  GBT** — this bears directly on our 5.98% FDP.
- a composite score silently dropped by a **prefix filter** (`MAIN_VAR_` failing a
  `VAR_` test) — the identical class as our `ENTRAP_` starts-with failure that hid
  14.7% of the library. **Grep for every prefix test we have.**
- `-fdr_pi0` named as one of "the two most seductive fake gains", resting on a
  uniformity assumption the data violates. **Our `lda.h` uses pi0 by default.**

**P4 — the prefilter, with explicit caution.** It made the predecessor's search
tractable (7.15M → 423k). But this project measured that **growing the budget
lowers identifications monotonically** while the recall ceiling rises, and that
oracle force-admission bought only +215 precursors and **−45 proteins**. Port it
as a *budget* mechanism, not a recall mechanism, and only after P1–P3.

**Do NOT port:** the vendored OpenMS fork. Read-only OpenMS is a user constraint
and the cleanest thing this project has.

## 4. The stepped plan

```
  S1  FDR CALIBRATION            5.98% -> nominal. DIA-NN achieves 0.914% on this
      (the priority)             file, so it is attainable. Audit pi0, the GBT
                                 leaf floor, and decoy exchangeability. Gate on
                                 entrapment FDP, never on counts.
  S2  MEMORY + TIMING            tcmalloc (P1) and the phase breakdown (P2).
                                 Cheap, measurable, and independent of S1.
  S3  TWO INSTRUMENTS            S08 at the known-good configuration. Running now.
                                 Everything to date is Astral-only.
  S4  PASS-1 YIELD               pass 1 contributes 554 of 3,625; the refit does
                                 the rest. Find out why before adding features.
  S5  THE ID GAP                 only here does chasing 3,625 -> 7,008+ make
                                 sense, and only with S1 holding FDP at nominal.
```

**S1 first, and this is the argument:** 3,625 at 5.98% is worth less than 2,000 at
1%. The predecessor's 7,008 is a nominal number with no FDP at all, so "catch up
to 7,008" is not even a well-posed target until our own error rate is honest.
Every count in both projects is currently a nominal count.

## 5. What neither project has

An entrapment-validated number from **either** codebase on the **same** file with
the **same** library. The predecessor measured empirical FDR only on a frozen
2.07M-feature table (4,367 precursors at 0.96% with LDA; 6,798 with GBT at matched
empirical FDR). Producing one comparable pair is the only way to settle which
architecture is actually ahead.
