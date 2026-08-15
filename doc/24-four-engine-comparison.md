# Four targeted-DIA engines, phase by phase

Draft 1, 2026-08-15. Sources: `DIA-NN-workflow-handoff.md` (read from `diann.cpp`
1.7.x, CC BY — documentation, not a code source), `OpenSWATH-workflow-and-findings-handoff.md`,
`ODIA-agent-handoff.md` (V1, tag `odia-v0.4.0`), and the V2 source tree plus this
project's own measurements. Literature-gap section pending a web survey.

**Scope caveat, stated first.** DIA-NN's description is from reading 1.7.x source;
2.x (QuantUMS, tunable decoy models, proteoform confidence) is closed and not
represented. OpenSWATH is described as the OpenMS 3.6 implementation. V1 and V2
are the two OpenDIAlyzer attempts and are private.

---

## 1. The four systems in one line each

| | what it is |
|---|---|
| **DIA-NN 1.7.x** | Purpose-built C++ engine. In-silico library from FASTA, its own extraction/picking/scoring, a linear classifier then a neural ensemble, iterative tolerance optimisation. |
| **OpenSWATH** | OpenMS-based targeted workflow: `OpenSwathWorkflow::performExtraction` + `MRMFeatureFinderScoring`, sub-scores persisted to `.osw`, statistical validation delegated to PyProphet. |
| **ODIA V1** | A *harness* around OpenSWATH's kernel — own prefilter, input layer, two-pass recalibration, in-process classifier and FDR, parquet output. 20-file vendored OpenMS fork. |
| **ODIA V2** | Full reimplementation. Own extractor, picker, scorer, calibration and FDR. OpenMS 3.6 consumed **read-only**; mzPeak native. |

---

## 2. Phase-by-phase

### 2.1 Library construction

| | approach |
|---|---|
| **DIA-NN** | FASTA digest → in-silico prediction (deep model or the non-DL path), elution groups, `Entry::generate_decoy` (`diann.cpp:3590`), proteotypicity annotation. |
| **OpenSWATH** | Consumes a prepared assay library (TraML/PQP/tsv); decoys via `OpenSwathDecoyGenerator`. Library construction is out of scope for the engine. |
| **V1** | OSWPQ/PQP/TSV → `LightTargetedExperiment`. `library_load` is **108.6 s single-threaded**, dominated by ~471 M small allocations rather than bytes. |
| **V2** | Own generator (`LibraryGenerator`) from FASTA + peptdeep; **compact Parquet** at 174 B/precursor against DIA-NN's 238; decoys stored, not regenerated. |

### 2.2 Raw-data access

- **DIA-NN**: own reader; `build_index` (`:7101`), `scan_RT[]`/`ms1_RT[]` cached for binary search; MS1 "windows" = observed span shrunk by `MinMs1RangeOverlap` 0.25 Th.
- **OpenSWATH**: `SwathMap` vector via OpenMS I/O; whole-map residency is the memory story.
- **V1**: mzML/mzPeak/`.d` → `vector<SwathMap>`; on mzPeak an mmap archive + `SpectrumStore`. `dia_run_load` 105.4 s at 45.6 cores.
- **V2**: mzPeak native, streaming; MS1 kept in a **separate index** (`MzPeakSource`) because merging would break the isolation-window assumption.

### 2.3 Chromatogram extraction

| | structure |
|---|---|
| **DIA-NN** | Per precursor: ≤6 fragment traces, plus `ms2_min` (running 3-point minimum), **shadow** traces at `−1.00335 Th`, **heavy** at `+1.00335/z` from the neighbouring window, **9 MS1 channels** (mono + ¹³C₁ + ¹³C₂ × 3 tolerances), and `nf` (precursor read out of MS2). |
| **OpenSWATH** | `performExtraction` into `MSChromatogram` per transition; RT windows from the normalised RT map. |
| **V1** | Delegates to `performExtraction`. Two passes: wide → recalibrate → narrow. **43% of wall at 138–146 of 224 cores — the best-behaved part of the tool.** |
| **V2** | Own streaming extractor with a `LiveSlot`/`BlockPool` model: a precursor's block is activated at RT-window entry, scored at emit, and returned to the pool. Residency scales with RT overlap, not library size. |

### 2.4 Peak picking

- **DIA-NN** (`Searcher::peaks`, `:7423`): at every scan position, require ≥1 fragment above `MinPeakHeight` and ≥2 in `{k−1,k,k+1}`; compute the **full pairwise Pearson matrix** of ≤6 traces over `[k−S,k+S]`; rank fragments by summed correlation; optionally add `corr(fragment, MS1)` (**`MS1PeakSelection`, default on**); require total ≥ `MinCorrScore` 0.5; smooth (¼–½–¼) and require a local maximum within `±max(S/3,1)`.
- **OpenSWATH**: `MRMFeatureFinderScoring` — peak-group detection on the summed trace, then sub-scoring.
- **V1**: OpenSWATH's.
- **V2**: ported DIA-NN's shape — pairwise Pearson **at detection time**, `MinCorrScore` as a hard gate, candidates kept **by margin not rank**. Replacing an amplitude picker with this moved apex-within-30 s from 41.0% to 84.1% and rank-1 accuracy 40.7% → 75.7%.

### 2.5 Sub-scores

| | count and character |
|---|---|
| **DIA-NN** | ~26 in 1.7.x incl. `pTimeCorr`, `pCos/pCosCube`, `pMs1TimeCorr`, `pNFCorr`, `pdRT`, `pResCorr(Norm)`, `pTightCorrOne/Two`, `pShadow`, `pHeavy`, `pMs1TightOne/Two`, `pMs1Iso*`, `pMs1Ratio`, `pBestCorrDelta`, `pTotCorrSum`, plus peptide properties. |
| **OpenSWATH** | `FEATURE_MS2` 38 columns, `FEATURE_MS1` **17**, `FEATURE_TRANSITION` 43 — xcorr shape/coelution (+weighted, contrast, combined), library similarity (dotprod, manhattan, RMSD, sangle), isotope correlation/overlap, massdev, MI variants, IM terms, elution-model fit. |
| **V1** | OpenSWATH's ~8–9 evidence families; **its composite main score never reached the classifier** — `MAIN_VAR_` failed a `VAR_` prefix filter. |
| **V2** | **19**, of which **18 read MS2 fragment traces** and one is MS1 (`MS1_COELUTION`). Measured: with a shared latent, feature *count* buys nothing (0 IDs at 4/8/15/30/60/110 alike) — **orthogonality is the lever**. |

### 2.6 Calibration

- **DIA-NN**: RT map refitted inside the iteration schedule; RT window from residuals; peak width → scan window; **mass calibration in bins with MS1 handled identically but `MassAccuracyMs1 × 5`**, MS1 delta kept only where `pMs1TimeCorr ≥ 0.90`; Q1 calibration; `--ref` reference-run mode. Its own README warns auto-tolerance is *"inherently noisy"*.
- **OpenSWATH**: `SwathMapMassCorrection` fits on **fragment anchors only** but derives a **separate precursor window** (`estimateWindow`, p99 × 1.3 padding).
- **V1**: `recalibrate_` fits library-RT → observed-RT from confident pass-1 anchors; `calibrateMassFromPass_` narrows ppm — **may only narrow**. **CiRT calibration was a routine 20.2 s phase.**
- **V2**: CiRT-seeded RANSAC line (this project, 2026-08-15) — Astral `915.0 + 11.399·iRT`, S08 `741.5 + 7.617·iRT`; mass from `MassCalibration::collect`'s own probe with a peakedness gate; MS1 mass calibration **not implemented** (extraction ran before calibration until yesterday).

### 2.7 Classifier and FDR

| | approach |
|---|---|
| **DIA-NN** | Phase A bootstrap/calibration → Phase B batch growth + tolerance optimisation (`MinCal` 1000, `MinClassifier` 2000) → Phase C linear refinement → Phase D final search, interference removal, **neural ensemble**. |
| **OpenSWATH** | Sub-scores → `.osw` → **PyProphet** semi-supervised LDA, then context-level FDR. |
| **V1** | **In-process** classifier (LDA or GBT), 3-fold **group** CV, 3 iterations; own q-values and peptide/protein rollup with picked competition. On a frozen 2.07 M-feature table: **4,367 precursors at 0.96% empirical FDR (LDA), 6,798 (GBT)**. |
| **V2** | Semi-supervised LDA/GBT/NN, 9 iterations; q-values with **Storey π₀** and the **Käll +1** finite-sample correction. Derived here: `q_min = π₀·(Ntar/Ndec)/targets_above_top_decoy`, so **~100 clean targets is the floor for a 1% claim**. |

---

## 3. What is common to all four

1. **Targeted extraction against a library**, not spectrum-centric search.
2. **Co-elution of fragments is the primary evidence.** Every engine's strongest sub-score family is some correlation among a precursor's own fragment traces.
3. **Two-stage tolerance**: search wide, calibrate, search narrow. Nobody searches once.
4. **Target-decoy with a semi-supervised discriminant.** All four learn weights from the data rather than fixing them.
5. **Peak-group competition**: several candidates per precursor, best-of-N kept — with the winner's-curse cost paid in the decoy tail.

## 4. Where they genuinely differ

| axis | DIA-NN | OpenSWATH | V1 | V2 |
|---|---|---|---|---|
| MS1 channels | **9** | mono + 4 isotopes × 4 charges | OpenSWATH's | **1** |
| MS1 at *detection* | **yes** | no | no | no |
| interference handling | explicit removal phase | none explicit | none | none |
| FDR engine | in-process | PyProphet | in-process | in-process |
| tolerance adaptation | iterative, per-run | one correction | narrow-only | narrow-only |
| decoy m/z | recomputed | recomputed | recomputed | **inherits target's** |

That last row is load-bearing and is V2's alone: because decoys carry the target's
precursor m/z, **any MS1-envelope-only feature is numerically identical for a
decoy and its target** and has zero discriminative power in production.

## 5. Measured standing (Astral, same file)

| | precursors @ nominal 1% | entrapment FDP | wall |
|---|---|---|---|
| DIA-NN | 12,308 | **0.914%** | 17:47 |
| OpenSWATH | 3,920 | 5.927% | days |
| V1 | 7,008 (vs its DIA-NN 9,261 = 75.7%) | not measured | 22:33 |
| **V2** | **3,625** | **5.983%** | 2:33:35 |

**Every count here is nominal.** Wen et al. (*Nat Methods* 22:1454, 2025) report
that no DIA search tool consistently controls FDR at the peptide level, with
DIA-NN's true precursor FDP above 2.3% at a nominal 1%.

## 6. Open routes in the literature — from a 105-agent verified survey

Run 2026-08-15: 105 agents, 4.44 M tokens, 1,137 tool calls, 3-vote adversarial
verification per claim. 12 findings survived; 9 candidate "gaps" were **refuted**
as already implemented.

### 6.1 Genuinely open

**Deep scoring on raw signal, replacing the handcrafted sub-score vector.**
The clearest verified gap (3-0, merged from three independent claims). DreamDIA
(2021) feeds a 16-dim LSTM embedding — two layers of 128/64 neurons, input dropout
0.4, recurrent 0.3 — to an XGBoost discriminator; DIA-BERT (2025) uses an
encoder-only 8-block transformer. **No mainstream engine ingests raw XIC tensors.**
DIA-BERT's authors name the handcrafted-feature dependency as the design limit of
current tools — a direct critique of the OpenSWATH/DIA-NN lineage that both ODIA
attempts inherit. Reported gains over DIA-NN: **+22% precursors / +51% proteins**
library-based, **+56% / +73%** library-free, pretrained on 276 M precursors from
952 files.

**A much broader XIC set per precursor.** Mainstream engines score on **6–10**
library fragments. DreamDIA extracts **170 XICs across six types** — library, self,
qt3, ms1, iso, light — explicitly *not* implemented in DIA-NN, OpenSWATH, Skyline
or Spectronaut. V2 currently uses ~12 fragments plus one MS1 channel.

**Deep-learning peak DETECTION.** Open question, and the survey found no
implementation: DIA-BERT and DreamDIA are both **scorers over library-anchored
candidates**. Proposing candidate peak groups from raw signal appears unexploited.

**Per-precursor adaptive m/z or RT tolerance.** OpenSWATH's tolerances are
confirmed **global** (a single `mz_extraction_window`, default 50 ppm, uniform
across precursors). The one candidate found — dia-PASEF IM windows — was refuted
(see below). Whether truly-per-precursor adaptive tolerance exists anywhere is
still open.

### 6.2 Refuted — assumed gaps that are already CLOSED

**Do not propose these as novel.**

- **Per-fragment ion-mobility consistency scoring** — shipped in DIA-NN 1.8.1 for
  dia-PASEF, with outlier fragments penalised even at high elution correlation.
- **2D peak picking in 1/K0 × m/z space** — DIA-NN's dia-PASEF module finds local
  maxima in the combined space rather than reducing to m/z-only chromatograms.
- **Adaptive IM extraction windows** — set per run from the alignment between
  confident identifications and the library, i.e. data-driven adaptive tolerance
  in an extraction dimension already exists.
- **Per-peak interference deconvolution** — DIA-NN 1.7–1.8 selects the fragment
  least affected by interference (highest summed Pearson against the others, from
  the top six by library intensity), treats its profile as the true elution shape,
  and subtracts. **OpenSWATH has no analogue**, and neither does V2.

### 6.3 FDR — the finding that bears directly on our 5.98%

**DIA-NN's FDR estimate is a plain target-decoy count ratio with NO π₀
correction**, claimed conservative by construction. **V2 applies Storey π₀ by
default** (`lda.h:227`), and V1's handoff independently names `-fdr_pi0` one of
"the two most seductive fake gains… resting on a uniformity assumption this data
violates."

Two independent lines now point at π₀ as a prime suspect for our FDP being 6×
nominal while DIA-NN achieves 0.914% on the same file. **This is the first
concrete, testable hypothesis for S1** and it costs one rerun with π₀ disabled.

Structural evidence of DIA FDR failure, verified:

- **Rosenberger et al. 2017**: 1% control at the peptide-query level alone lets
  inferred proteins accumulate toward the *library's* protein count across 229
  runs; corrected by a **global analyte constraint** — best peak group per analyte
  across all runs, used as a master list.
- **The OpenSWATH paper itself** reports its estimated FDR underestimating the
  manually determined false-positive rate by **0.9 points at nominal 1%**.

Wen et al. (2025) remains an open question the survey did not fully resolve — what
exactly it measures, by how much each engine exceeds nominal, at which level, and
which corrections are proposed.

### 6.4 Lineage, for the paper's framing

mProphet (2011) established combining heterogeneous chromatographic sub-scores
into one semi-supervised discriminant trained against decoys → OpenSWATH (2014)
transplanted it to DIA → PyProphet made the LDA experiment-wide with per-level
q-values → DIA-NN 1.7–1.8 kept the same handcrafted vector but replaced the linear
combiner with a **12-network feed-forward ensemble** and added interference
subtraction. **All four systems in this paper sit on that single lineage**, which
is precisely what the deep-scoring literature proposes to leave.
