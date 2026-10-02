> **SUPERSEDED BY doc/26-PLAN.md (2026-08-15).** Kept as the working record. Where this document and doc/26 disagree, doc/26 wins — it carries the corrections this one predates.

# The known-good configuration (2026-08-15)

The first configuration that identifies at proteome scale. Recorded so it is not
lost and not re-derived.

```bash
OpenDIAlyzer \
  -in  astral.mzpeak \
  -tr  odia_lib_parity.parquet \
  -irt_slope 11.399 -irt_intercept 915.0   # CiRT-seeded RANSAC line, per-run
  -rt_window_pass1 60 -rt_window 60        # BOTH passes narrowed
  -gate_alpha 0.05                         # Gate C
  -live_memory_gb 400 -threads 48
  # NOTE: no -pass1_precursors. See below -- it is a STRIDE.
```

**Result:** 3,625 precursors at 1% FDR, 13.3 M peak groups, 212 GB peak RSS,
2:33:35. DIA-NN on the same file: 12,308 in 17:47. OpenSWATH: 3,920.

## The three flags that mattered, and why

**`-pass1_precursors` is a STRIDE, not a cap.** At 100000 against 9,983,789
precursors it extracts every ~99th, so only ~101 of the 10,040 findable peptides
were searched. Per-present-precursor efficiency was never the problem — 55%
strided against the oracle's 53% — the *search* was 1% of the library. Removing
it took the clean seed 56 → 1,098 and identifications 0 → 3,625.

**`-rt_window` does not govern pass 1.** Pass 1 defaults to `1.0e9` seconds
(`OpenDIAlyzer.cpp:1432`) so it can find anchors anywhere. `-rt_window_pass1` is
the flag. Before this was found, three consecutive runs extracted 119,088,506
transitions each, identical to the digit, and "calibration changes nothing" was
an artefact of never having restricted anything.

**A run-specific iRT line is required and is cheap.** CiRT-seeded RANSAC, ~35 min
blind search of 149 precursors, recovers the line to 3.9% slope / 13 s intercept
on Astral and 0.6% / 7.3 s on IH1 (`doc/19`). The line is per-run: Astral
915.0 + 11.399·iRT, IH1 741.5 + 7.617·iRT.

## The q-value floor, which explains every previous zero

`lda.h:209-265` uses the Käll +1 finite-sample correction, so at the extreme tail

```
    q_min = pi0 * (Ntar / Ndec) * (1 / targets_above_the_top_decoy)
```

**~100 clean targets is the hard minimum for 1% FDR to be reachable.** Every past
zero was this floor refusing to certify a claim the data could not support — 27
observations cannot demonstrate a 1% error rate, and an estimator that said
otherwise would be lying.

**Use `targets-above-the-top-decoy` as the diagnostic.** It is computable from any
scored TSV with no FDR run, it predicted every outcome measured so far, and it
converts "did this help?" into a number:

| run | above top decoy | q floor | IDs |
|---|---|---|---|
| oracle (10,040 present) | 5,357 | 0.00041 | 6,815 |
| full library, strided | 56 | 0.0365 | **0** |
| full library, unstrided | 1,098 | 0.00199 | **3,625** |

## What is NOT established

- **FDP is unmeasured.** The run reports `-entrapment_prefix 'ENTRAP_' matched NO
  library precursor`, so 3,625 has no demonstrated error rate and is **not** yet
  comparable to DIA-NN's 12,308 at 0.914% FDP. Rerun on the entrapment library.
- Pass 1 contributes 554 of the 3,625; the pass-2 refit does the rest.
- Astral only. 8.6× slower than DIA-NN.
