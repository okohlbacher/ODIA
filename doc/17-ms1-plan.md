# The MS1 arm: rebuild plan

Status: PLAN. Written 2026-08-14 after adversarial review by Codex `gpt-5.6-sol`
at `model_reasoning_effort=max` and Kimi 0.36.0, both **with vault access** (the
first round to have it). Vibe returned 0 bytes for the third time and did not
participate. Raw: `bench/review_ms1/`.

**This round found three real code defects by reading, before any measurement.**
That is the argument for `doc/review-brief-preamble.md`.

---

## 0. The defects, all verified in the tree

**1. The scorer never receives the MS1 traces.** `extractAndScore_`
(`OpenDIAlyzer.cpp:2792-2796`) calls `scoringOptions_()` and constructs
`PeakGroupScorer::Sink(library, options)` **before** calling `extractInto_`,
which builds the matrix at `:1048`. `scoringOptions_()` sets
`options.ms1 = ms1_traces_.empty() ? nullptr : &ms1_traces_` — and at that moment
it *is* empty. The Sink therefore holds a null pointer for the entire run.
Measured consequence: `MS1_COELUTION: 0 finite, 11,877 no MS1 available`.

Kimi predicted this from the source and predicted the instrumented run would
report `ms1_null = 11,877, everything else 0`.

**2. `-live_memory_gb` does not gate the MS1 build.** The skip check at `:1033`
reads `options.live_memory_budget_bytes`, which is not assigned until `:1112`.
`cap` is therefore always 0, `cap > 0.0` is false, and the **148 GiB matrix is
allocated unconditionally** whatever the user asked for. Kimi: "a matrix nobody
approved."

**3. MS1 is extracted before mass calibration.** Build at `:1048`,
`applyMassCalibration_` at `:1066`. So the traces are matched at uncalibrated
masses through a borrowed `fragment_ppm`. **The 62.1%-with-signal figure was
measured under a defective configuration and must not be used to size anything.**

## 0.1 Corrections to numbers this project has been quoting

From Kimi, checked against `BACKLOG.md`:

- **The 13.7× MS1/MS2 co-elution enrichment is a 14-precursor top bin, ±3.6×.**
  The bulk effect is **1.4–1.5×** (SMD 0.296 sd). It has been repeatedly cited as
  if robust. It is a good *discovery* statistic and not an acceptance criterion.
- "49,999 of 50,000" conflates **three** saturation counts (49,999 whole-frame,
  49,835 mobility-sliced, 49,916 isotope depth).
- The **616.7 GB** OOM figure is uncorroborated; the code says **588 GB** twice.
- Rules about search-budget growth and winner's curse come from **unmarked seeded
  vault notes**, i.e. the *previous* codebase. Transferable in mechanism, not as
  measurements of this one.

---

## 1. THE TRAP THAT RESHAPES THE ISOTOPE PLAN

**ODIA's production decoys carry the target's precursor m/z** — a decoy row keeps
its target's sequence with only the *fragment* m/z shifted
(`DIANNLibraryFile.cpp:757-765`). Therefore **any feature of the MS1 envelope
alone — `isotope_corr`, `isotope_overlap`, MS1 massdiff — is numerically
IDENTICAL for a decoy and its target, and has exactly zero discriminative
power.**

On the `v6_50k` search benchmark decoys are sequence-shuffled with recomputed
precursor m/z (`LibraryGenerator.cpp:775-776`), so the same features *do*
separate there. **A pure-MS1 feature can therefore look excellent on the search
benchmark and be inert in production, and nobody has measured the two side by
side.** Codex reached the same place from the other direction: non-exchangeable
decoys let an isotope feature "look spectacular while making q-values
anti-conservative."

**Consequence: only JOINT MS1×MS2 quantities can discriminate against production
decoys.** Every isotope feature below is defined as a joint shape statistic, not
an envelope statistic. This is the single most important finding of the round.

---

## 2. Storage: kill the matrix; do not replace it with CSR

Both reviewers rejected the sparse-CSR proposal, for different and both-correct
reasons.

- **Codex:** CSR beats dense only below ~67% cell density with 16-bit indices;
  ten million independently allocated sparse vectors bring metadata,
  fragmentation and cache costs; and zeros *inside* a trace are meaningful for
  shape. Also: 62.1% is the fraction of precursor ROWS with some signal, which
  says nothing about cell density — "calling that 62.1% occupancy is already too
  strong."
- **Kimi:** the codebase already rejected sparsity for this exact problem at
  51.3% occupancy (`ChromatogramExtractor.h:118-122`). Keeping a persistent,
  whole-run, spectra-keyed matrix in any encoding keeps the wrong *shape*.

**The MS2 side solved this already and MS1 was exempted when libraries were
small** (`OpenDIAlyzer.cpp:1026`: "At 5,330 precursors it is 83 MB and
invisible"). The streaming extractor activates a per-precursor block at
RT-window entry, scores at emit, and returns memory to the pool.

**Fold MS1 into the same live block, keyed to the same RT window.** A ±60 s
window at a ~1.8 s MS1 cycle is ~67 bins; 4 channels × 67 × 4 B ≈ **1 KB per
live precursor**, against ~40 KB for the MS2 block. Residency scales with RT
overlap, not library size. **148 GiB → GB-scale**, and it inherits the existing
budget machinery — which becomes honest once defect 2 is fixed.

Nothing in DIA-NN's four uses needs persistence: detection, scoring, calibration
deltas and `ms1_area` are all computable at emit time. Random access to another
precursor's MS1 is never required.

Engineering cost: MS1 spectra live on a separate index interleaved in acquisition
order, so the forward pass needs a second matching path — a second `MzIndex`, not
a redesign. Keep the frame-wise **max** aggregation; it is deliberate.

**Explicitly rejected:** memory-mapping the 148 GiB matrix ("converts RSS into
page-cache thrashing"), and persisting DIA-NN's three tolerance widths as
separate channels — query the widest interval once and derive nested-window
statistics in the same scan.

---

## 3. The MS1 feature map (FeatureFinderCentroided) — a DIFFERENT object

Measured on Astral MS1 (137 MiB extracted in 39 s; the whole raw file is 3.1 GiB
against the matrix's 144.6 GiB — a 1,080× inflation):

| | |
|---|---|
| features at default sensitivity | 18,446 (30 MiB, 1:42, 768 MB RSS) |
| DIA-NN q≤0.01 precursors matched (±20 ppm, ±30 s, charge) | **3,670 — 29.8%** |
| misses **with** non-zero DIA-NN `Ms1.Area` | **8,169** |
| features with a UNIQUE library match at ±5 ppm | **365 (2.0%)**, RT span 57–2314 s, **0 empty deciles** |
| unique at ±10 / ±20 ppm | 260 / 196 |

**Role, decided by these numbers:**

- **Prefilter — NO.** 29.8% recall and losses are unrecoverable. (Distinct from
  the budget rule: that concerns admitting noise, this is discarding signal.)
- **Anchors — YES.** The requirement inverts: anchors need precision and RT
  coverage, not recall. 365 unambiguous, unconditioned `(library iRT, observed
  RT)` pairs covering the whole gradient is a viable calibration set, and it is
  the **exogenous `a,b` transform** that `doc/16` §2.2 identified as missing —
  peptdeep emits iRT, a rank, not run seconds.
- **Global run overview, expanded by MS2 — YES, as the destination.** The 73%
  ambiguous features are not waste; they are shapes awaiting evidence that names
  them, and MS1 alone cannot provide it.

**Uniqueness improves as tolerance tightens** (365 at 5 ppm vs 196 at 20 ppm), so
mass calibration first buys RT anchors — independent support for `doc/16`'s
reordering.

**Tension to respect:** sensitivity helps the overview role and *hurts* the
anchor role (more features → more ambiguity). Run the finder twice at different
settings rather than seeking one map for both.

**Unmeasured and required before the anchor role is real:** precision. "Unique"
is not "correct". Test: do those 365 pairs' library iRT correlate with observed
RT? If not, the role collapses.

This is complementary to §2, not competing: §2 is the per-candidate scoring path,
§3 is a run-level summary. Neither reviewer was asked about §3 — their brief
predates it.

---

## 4. Features: what to build, what never to build

**Build, in order:**

1. **`corr(mono MS1, summed fragments)`** — the existing one. Fix, don't replace.
2. **Native-grid alignment (E′).** Interpolate onto the MS2 cycle grid **at score
   time** from raw per-scan values in the window — no persistent upsampling, and
   not a literal copy of DIA-NN's stored interpolation. On S08 the grid is 1,343
   MS1 vs 16,105 MS2 (~12:1), so several MS2 cycles snap into one MS1 bin and the
   Pearson leg goes flat.
3. **Three physical channels only: `M−1`, `M`, `M+1`**, at the library charge.
   **`M−1` before `M+2`** — Codex, explicitly against the expected preference for
   `M+2`: `M−1` asks a different *identity* question ("is the nominated
   monoisotope somebody else's isotope?"), while `M+2` is more of the same
   envelope at lower signal. `M+2` only if it adds held-out conditional
   information.
4. **Common chromatographic-shape fit** — fit the summed-fragment shape to `M`
   and `M+1` with one nuisance amplitude per channel; score explained shape or
   normalised residual. **This replaces my proposed "isotope-ratio time series
   correlated with the monoisotope", which Codex refuted:** a true isotope ratio
   is ~constant across elution, so correlating it with `M` measures denominator
   noise and signal strength — the abundance trap.
5. **Coherent pre-isotope overlap penalty** at `M − 1.00335/z`, requiring **both**
   substantial normalised amplitude **and** co-elution before penalising. A lone
   lower-mass peak is not enough. Use as a **gate**, not a standalone sub-score —
   §1 means it cannot separate target from decoy in production.
6. **One tight-tolerance mono variant** (DIA-NN's 0.45×). MS1 is
   interference-rich and tightening is the cheap selectivity lever.

**Never build:** isotope depth or count (measured 1.0×, three times — stop
re-testing this class); MS1 apex/area as classifier inputs; `pMs1Ratio`-style
level ratios (rule 6 killed `ms1_max` at 1.7× for exactly this); DIA-NN's nine
stored traces; OpenSWATH's 17 columns wholesale; four alternative charges (the
library charge is known); maximised cross-correlation over many RT lags (each lag
is another winning opportunity); MS1 quantification as part of this work.

---

## 5. D — MS1 at detection — is safe ONLY as fixed-budget preselection

Kimi: DIA-NN's `MinMs1Corr` defaults to **−INF**, so `MS1PeakSelection` rejects
nothing and merely re-ranks within the existing candidate set — N unchanged, null
width unchanged, symmetric across labels.

Codex sharpens this into three implementations that behave differently:

| implementation | effect on the null |
|---|---|
| MS1 proposes extra peaks / rescues precursors | **widens — unsafe** |
| add MS1 correlation to every candidate score, then max over the same N | adds variance to the decoy score before a maximum — **can widen** |
| MS1 as a fixed veto/ranker returning at most one of the existing candidates | **can shrink** |

Only the third is safe. The protocol: freeze the candidate set before observing
MS1; apply a frozen out-of-fold gate; never add a peak, widen an interval or
admit a precursor; identical missing-data fallback for targets and decoys;
target-decoy competition **after** selection; keep the preselection audit table.

**Absent MS1 must be NaN-neutral, never penalised** — ~38% of precursors have no
MS1 signal, and penalising absence imports the abundance bias into *detection*,
where it is far harder to audit than a dropped sub-score.

**D is downstream of §0 and the RT bootstrap, not parallel.** If candidate windows
do not overlap MS1 signal, detection-time MS1 is noise.

---

## 6. Acceptance gates — pre-registered, layered

Rejected: ID counts (629 → 928 → 799); total mass scatter (blind: 4.66 → 4.19 ppm
while 68% of systematic error was removed); DIA-NN's confident set for feature
selection (and that reference list was itself measured wrong — 1,230 omitted real
IDs, 205 non-IDs, Jaccard ~66%).

1. **Functional.** Every scored candidate lands in exactly one diagnostic
   category; no null-pointer, empty-store or row-mapping failure; a silently
   dropped feature column is a run **warning or error**. Today's `0 finite,
   11,877 unavailable` fails here before discrimination is even considered.
2. **Fixed-population incremental information.** Out-of-fold
   `ΔLL = LL(label | existing scores, MS1) − LL(label | existing scores)`, folds
   separated **by run and by peptide/protein**, not random rows. Paired CI must
   exclude zero. Repeat within bins of intensity, charge, m/z, candidate count and
   MS1 availability — this is what exposes an abundance proxy.
3. **Negative controls** (Codex, and the best single idea of the round):
   RT-circularly-shifted MS1 traces outside the peak width; wrong isotope spacing;
   mass/RT-matched wrong pairings. These preserve smoothness, occupancy and
   abundance, so incremental information must collapse to baseline. **Stronger
   evidence of real co-elution than any comparison with DIA-NN's IDs.**
4. **Candidate-null gate for D.** At fixed candidate set: candidates retained per
   label, decoy max-score high quantiles, top-minus-second, RT displacement,
   availability by label. Passes only if target excess improves while the decoy
   tail is non-inferior. Margin set from repeated baseline runs **before**
   inspecting D. Assert candidate count within ±0.5%.
5. **Entrapment FDP at nominal 1% AND 5%, with Poisson intervals.** The
   `var_im_spread` lesson: 3 vs 6 entrapment hits "halved FDP", reversed at
   q≤0.05, and cost 60 IDs. The 1% threshold is set by ~70 decoy events, ±23% at
   95% — a feature must move *that tail*, not the bulk.
6. **Two instruments, always.** `MS1_COELUTION` itself was +17.6% on S08 and
   **+0.4% on Astral**. A feature that wins on one instrument is flagged, not
   shipped.

---

## 7. Order of work

```
  0. FIX THE PLUMBING          null Sink ptr; dead budget guard; MS1 after
                               calibration; -ms1_ppm; split guard counters
                               -> gate: MS1_COELUTION produces finite values
  1. E': native-grid alignment interpolate at score time
  2. VALIDATE THE MONO ARM     gates 2 + 3 on the EXISTING single feature
                               -> if 13.7% bulk 1.4x does not survive negative
                                  controls, the whole MS1 premise is wrong
  3. B': streaming MS1         fold into live blocks, ~1 KB/precursor
  4. FEATURE MAP AS ANCHORS    precision check first, then the a,b transform
  5. D': fixed-budget          veto/ranker only, protocol in section 5
     preselection
  6. C': three channels        M-1, M, M+1; joint shape statistics only
```

Steps 0–2 are cheap, and step 2 is a genuine kill gate for everything after it.

## 8. What the plan still does not address

- **The MS1 arm's real job may be ignition, not discrimination.** If the
  classifier cannot bootstrap at a 1.5% prior, A–E all measure as zero again.
  That argues for *one* orthogonal feature done properly, not a broad programme.
- **"DIA-NN leans on MS1" is source-read, never ablated.** Cheap prior worth
  running first: disable DIA-NN's MS1 features on one file and measure the delta.
  We are rebuilding an arm whose value in the reference engine is assumed.
- **MS1 mass calibration is absent from this plan** and may be the largest payoff
  — DIA-NN recommended 4.6 ppm MS1 having calibrated at 22 ppm, and its third use
  of MS1 (deltas kept where `pMs1TimeCorr ≥ 0.90`) is impossible today because of
  defect 3.
- **Every MS1 feature's value is decoy-scheme-dependent** (§1) and that dependence
  has never been measured.
