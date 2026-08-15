> **SUPERSEDED BY doc/26-PLAN.md (2026-08-15).** Kept as the working record. Where this document and doc/26 disagree, doc/26 wins — it carries the corrections this one predates.

# The joint reference map: one anchor structure for calibration and fine-tuning

Status: **v1 REJECTED by review round 7 (codex at effort max, kimi). v2 below.**
Written 2026-08-14. Raw: `vault/70-Adversarial/round7-*`.

**The decision to build one joint structure stands — the user set it and both
reviewers endorsed it. What failed is v1's ANCHOR MECHANISM, and it failed a
measurement, not an argument.**

**This document is the spine.** `doc/15` (mass calibration), `doc/16` (phases 3+4)
and `doc/17` (the MS1 arm) each solved a piece and each assumed an anchor source
that did not exist. This one names the anchor source and re-hangs the other three
off it.

---

## 1. What it is

A **run-level map of detected features, built from MS1 and MS2 evidence together**,
independent of the library, and used as the anchor set for RT calibration, mass
calibration and fine-tuning.

One entry:

```
  m/z, RT apex, RT start/stop, charge
  elution profile (shape, FWHM)
  isotope envelope (M-1, M, M+1 amplitudes)      [MS1 side]
  fragment co-elution shape + transitions        [MS2 side]
  provenance: MS1 | MS2 | BOTH
  identity: unassigned | unique library match | ambiguous(n)
  1/K0 + mobility width                          [diaPASEF only]
```

Cost is trivial: 95,216 features × ~64–128 B is **6–12 MB**. This is not a memory
question. It is an architecture question.

## 2. Why this is the right spine — four problems it closes at once

**(a) Every calibration in ODIA is currently anchored on something conditioned.**
`MassCalibration::collect` uses its own brightest-co-occurrence probe; the RT map
is fitted from FDR-accepted peak groups. The measured cost of the latter class of
error is on the record twice: sizing the fragment m/z window from pass-1 IDs took
Astral from **4,969 to 2,422**, and four separate anchor statistics pointed the
wrong way during the RT work. A data-driven map is **unconditioned by
construction** — it does not know what the library contains or what scoring
accepted.

**(b) ~~It is the exogenous transform `doc/16` §2.2 said Option B required.~~**
**RETRACTED — see §4.** The claim rested on 1,794 unique `(library iRT, observed
RT)` pairs, and those pairs measured at **r = 0.045** against a shifted-library
control that produced **more** of them than the real library. The map does not
supply the `a, b` transform. That comes from the free metadata iRT→gradient
transform instead (88.89% at ±60 s), per §5.

The underlying problem `doc/16` §2.2 identified is unchanged and real: peptdeep
emits **iRT — a rank, not run seconds** — so `RT = a + b·iRT` needs its
coefficients from something available before any identification. This document
proposed the wrong source for them.

**(c) It makes MS1 mass calibration possible at all.** `doc/17` defect 3: MS1 is
extracted before `applyMassCalibration_`, so DIA-NN's third use of MS1 —
calibration deltas kept where `pMs1TimeCorr ≥ 0.90` — cannot be implemented today.
The vault's own measurement says this is where the payoff is: DIA-NN recommended
**4.6 ppm MS1** on this file having calibrated at 22 ppm, and *"m/z is the cheap
selectivity lever, RT the expensive one."*

**(d) It is the persistent structure `doc/16` §6.1 needs to iterate without
re-decoding.** Two decodes instead of eight was that design's central claim; the
map is what survives between them.

## 3. The one thing it is NOT

**Not a prefilter.** Measured, and this is settled:

| sensitivity | features | recall of DIA-NN q≤0.01 | misses **with** real `Ms1.Area` |
|---|---|---|---|
| default (10/1) | 18,446 | 29.8% | 8,169 |
| **mid (5/2)** | **95,216** | **51.0%** | 5,598 |
| hot2 (4/1) | 64,981 | 41.8% | 6,703 |

Even at best it loses half of what DIA-NN finds, and prefilter losses are
unrecoverable. **The map indexes; it never gates.**

**COST DECIDES THE ANCHOR CONFIGURATION, and it is the cheapest one.** Every
setting below reaches all ten RT deciles, so coverage does not discriminate:

| config | `min_spectra`/`max_missing` | wall | features | anchors (±5 ppm) | **anchors/min** |
|---|---|---|---|---|---|
| **default** | **10 / 1** | **1:42** | 18,446 | 365 | **215** |
| g10_4 | 10 / 4 | 5:47 | 40,979 | 871 | 151 |
| g8_3 | 8 / 3 | 8:11 | 50,185 | 1,049 | 128 |
| mid | 5 / 2 | 22:00 | 95,216 | 1,794 | 82 |
| hot2 | 4 / 1 | 28:14 | 64,981 | 1,154 | 41 |

Anchors scale **sub-linearly** with compute -- 13x the runtime buys 4.9x the
anchors -- and an RT calibration curve is monotone and smooth, so it needs
*coverage*, not thousands of points. **Use default settings for the anchor map:
1:42 wall, 749 MB, 365 anchors spanning 57-2314 s with no empty decile.** Plus a
one-off 39 s / 218 MB `FileFilter` to extract MS1.

`mid` remains the candidate for the run-level *overview* map (role 3) if the
extra features are ever shown to earn their place -- a separate question, decided
separately. `hot2` is strictly dominated: 28 minutes for fewer features than
`mid`'s 22.

Two measured lessons from the sweep: both trace parameters matter and
**`min_spectra` is the stronger lever** — at `max_missing`≈1, taking it 10→5
gives 18,446→95,216 features; at `min_spectra`=10, `max_missing` 1→4 gives
18,446→40,979. An earlier draft of this section called `max_missing` dominant on
the strength of mid-beating-hot2, but those two differ in **both** trace
parameters and cannot isolate either. What survives from that claim is narrower:
dropping `min_spectra` 5→4 did not compensate for losing one gap of tolerance.
And recall and
anchor count move **in lockstep** (0.30/365, 0.51/1,794, 0.42/1,154), so the
anchor and overview roles do **not** compete; a single sensitivity axis serves
both. (I predicted they would compete. They do not.)

## 4. WHAT KILLED v1 — measured, not argued

I proposed selecting anchors as MS1 features matching **exactly one** library
precursor within tolerance, and queued a precision check I had not yet run. Kimi
ran it.

| test | result |
|---|---|
| Pearson / Spearman of the 1,794 pairs, library iRT vs observed RT | **0.045 / 0.047** |
| **+30 ppm shifted (all-false) library, density matched** | **2,848 "unique" matches — 158% of the real 1,794** |
| same control on the default map | 605 vs 365 real — **166%** |
| overlap with DIA-NN's confident set | 3 of 1,794 |

**A deliberately wrong library yields MORE unique matches than the real one.** The
unique bucket sits *below* its own chance baseline, so it contains essentially no
true pairs. The mechanism: **real present peptides land in DENSE library
neighbourhoods and therefore fall into the AMBIGUOUS bucket. The unique bucket is
where sparse-neighbourhood coincidences live.**

Codex reached the same place structurally, without the measurement:

> The bootstrap is **residual-censored**. Step 2 selects matches within ±5 ppm
> *before* mass calibration, then step 4 estimates calibration from them. True
> anchors outside the window are excluded while chance matches near zero residual
> survive. The fitted model recovers **the selection window**, not the instrument
> error. Tightening 20 → 5 ppm mechanically increases "uniqueness"; it says
> nothing about correctness.

v1 cited that very trend (365 at ±5 ppm vs 196 at ±20) as *evidence*. It is an
artefact.

**And the bootstrap would have been worse than doing nothing:** it replaces a free
degenerate transform covering **88.89% at ±60 s** with a fitted one covering
~4.5%.

**The conceptual error, in kimi's words:** *the plan confused "unambiguous" with
"unconditioned". Unconditioned anchoring is possible, but via robust consensus
over ambiguous evidence, not via uniqueness filtering.*

Four further blocking findings from codex, all accepted:

- **"Unique library match" cannot mean identity at all.** Production decoys share
  the target's precursor m/z, so a mass match is *at minimum* target/decoy
  ambiguous. Only a precursor-mass equivalence class is defensible, and that is
  not peptide identity.
- **The map is not library-independent after step 1.** MS1 *detection* is;
  library mass *matching* is not, and MS2 candidate extraction is explicitly
  conditioned. Keep the claim narrowly, at the observation layer.
- **The `BOTH` merge is unsound.** MS1↔MS2 is many-to-many, and an MS2-only DIA
  observation has no observed precursor m/z or charge — filling those from the
  library turns a hypothesis into an observation. Use a versioned association
  *edge* (isolation-window compatibility, RT-profile overlap on native time
  grids, charge/isotope hypotheses, mobility overlap). Never average or overwrite
  apexes; disagreement is QC evidence.
- **"Indexes, never gates" is operationally false.** A map-derived calibration
  that narrows extraction windows excludes evidence — that *is* a gate, and it is
  precisely the 4,969 → 2,422 failure. Needs a narrow primary path plus a fixed
  **wide rescue/audit path**, with unsupported strata keeping wide windows.

Also caught: v1 fits RT before mass in §5 while §6 argues mass-before-RT. A real
internal contradiction.

---

## 4b. THE DETECTABILITY RISK (v1's §4, still live but no longer the sharpest)

**Fitting calibration only where the map has entries conditions the calibration on
what the feature finder detects.** That is the
[[accepted-groups-are-a-biased-sample]] trap wearing a new costume: the map is
unconditioned *by the library and by scoring*, but it is emphatically **conditioned
by detectability** — abundance, isotope-pattern regularity, elution-peak shape.
A calibration fitted from its entries describes well-behaved abundant precursors,
which is precisely the population that was never the problem.

The 4,969 → 2,422 collapse happened because a sound estimator was fitted on a
selected sample. Nothing about this map prevents a rerun of that.

**Mitigations to be reviewed, none of them proven:**
- Report anchor counts stratified by RT decile, m/z decile, charge and intensity
  quartile every round, and refuse to fit a stratum below a minimum count,
  inheriting the global fit there and marking it.
- Compare the map's intensity distribution against the run's whole peak
  distribution and state the shortfall rather than assuming it away.
- Validate the fitted map on precursors the feature finder did **not** detect —
  which is the only test that speaks for the population at risk.

## 4c. v2's STEP 2 IS ALSO REFUTED — measured 2026-08-14

The consensus/mode repair was tested before implementation, and it fails the same
way v1 did.

For each of 1,828,424 distinct library target precursor m/z, histogram the ppm
deviation to **every** map feature within the window; a true instrument offset
should appear as a peak above a flat background.

| map | window | mode | **prominence over background** |
|---|---|---|---|
| default (18,446) | ±50 ppm | +2.5 ppm | **1.1×** |
| default | ±15 ppm | +2.5 ppm | **1.0×** |
| default | ±8 ppm | +2.5 ppm | **1.0×** |
| mid (95,216) | ±15 ppm | +3.0 ppm | **1.0×** |

Shifted-library controls give the same 1.0–1.1× in every case. (Their modes track
the shift exactly, e.g. +30 ppm shift → −27.5 ppm mode, but that is arithmetic,
not evidence.) **There is no peak. The histogram is flat.**

**A first version of this test took the NEAREST feature per precursor and found a
7.6× "mode" at 0 ppm — with the shifted controls at 7.0×.** Nearest-neighbour
distance peaks at zero by construction whatever the truth. The control caught the
estimator before it reached C++. Keep the controls.

### Why, and the rule it re-teaches

1.8M distinct precursor m/z over ~380–980 Th means the m/z axis is **saturated**:
at ±15 ppm every position has library entries, so 4.1M candidate pairs contain at
most ~18k true ones — 0.4% signal. **On a proteome-scale library, precursor m/z
ALONE carries no identifying information**, and no window width or map density
changes that. Both v1 (uniqueness) and v2 (mode) were single-coordinate matches
against a saturated axis.

This is [[odia-shape-not-presence]] again. The existing MS2 probe works precisely
because it uses **fragment co-occurrence** — many coordinates per precursor.

**v3, therefore: break the m/z degeneracy with a second coordinate.** Match
features to library precursors on m/z **and** retention time, using the free
metadata iRT→gradient transform for the RT constraint. A ±60 s window on a
~1,770 s gradient is ~3.4% of the run, which cuts background ~30× while keeping
true matches. That is the next measurement, and it must carry the same
shifted-library control.

## 5. v2 BOOTSTRAP — consensus, not uniqueness (STEP 2 REFUTED, see 4c)

**Replace uniqueness with mode/consensus estimation on BOTH axes. Neither needs
identity.**

```
  1. MS1 features            FeatureFinderCentroided, DEFAULT settings
                             (1:42, 749 MB, 215 anchors/min -- section 3).
  2. MASS offset, IDENTITY-FREE   For each library precursor, collect the ppm
                             deviation to nearby map features. TRUE offsets
                             CLUSTER; false matches are FLAT in ppm. Take the
                             MODE over thousands of precursors. This is
                             MassCalibration::collect's existing probe logic with
                             map features as cleaner, denser input -- a genuine
                             improvement available today.
  3. RT transform            DEFAULT: the FREE metadata iRT->gradient transform.
                             88.89% at +/-60 s with no anchors at all, and the
                             bar any fitted alternative must beat.
                             Fallback: DIA-NN's shape -- wide windows -> simple
                             fixed-criterion seeds (survey F5) -> refit -> tighten.
                             UNPROVEN option: RANSAC / median-of-residuals
                             consensus line fit over AMBIGUOUS candidates. Kimi's
                             first attempt at this failed; it is not the default.
  4. EXTRACT + detect MS2    Narrow primary path PLUS a fixed wide rescue/audit
                             path (section 4). Unsupported strata keep wide windows.
  5. MERGE as EDGES          Versioned association edges, not fused rows. Never
                             overwrite an apex.
  6. REFIT from the graph    Convergence on parameter STABILITY over data halves.
  7. ONE confirmatory pass
```

**The map's role is DEMOTED**: from "the escape from the circle" to "the MS1
mass-calibration residual source, and the persistent index". Step 2 is the part
that works today and needs no identity.

Codex's shape for the structure itself: **one typed, versioned evidence graph
with task-specific anchor views**, not one fused feature row and not one
universal anchor set. Eligibility differs per task — precursor-mass calibration
needs only that all hypotheses imply the same theoretical mass; RT calibration
needs a sufficiently concentrated iRT; fine-tuning needs independently supported
sequence identity.

### v1's bootstrap, kept visible because it was wrong in an instructive way

```
  1. MS1 features            FeatureFinderCentroided. No library, no calibration.
  2. UNIQUE m/z matches      1,794 unambiguous (library iRT, observed RT) pairs.
  3. FIT the coarse RT map   From those. Unconditioned by scoring.
  4. FIT the mass model      Same anchors, log(m/z) basis (t = 7.4). Tightening
                             ppm RAISES uniqueness -- 365 at 5 ppm vs 196 at 20 --
                             so mass and RT reinforce each other here.
  5. EXTRACT + detect MS2    Now at approximately the right RT and m/z.
  6. MERGE MS2 evidence      Co-elution-detected candidates enter the same map;
                             ambiguous MS1 entries acquire identity; unmatched
                             MS2 entries are added with provenance MS2.
  7. REFIT from the merged   Both calibrations, now from joint evidence.
     map                     Convergence per doc/16 section 4: parameter
                             STABILITY over data halves, never ID count.
  8. ONE confirmatory pass
```

Steps 1–4 need no identifications, no classifier ignition, and no MS2. That is
what makes this a bootstrap rather than another circle.

## 6. How the other plans re-hang off it

| plan | what changes |
|---|---|
| **`doc/15`** mass calibration | The additive model `f1(m/z) + f2(RT) + f3(log I) + f4(1/K0)` stands. Its anchor source changes from `MassCalibration::collect`'s probe to the map. Its "one wide extraction, retain residuals" inner loop is **demoted** (`doc/16` §2.1) — codex showed the narrow-window feedback is first-order, not second. |
| **`doc/16`** phases 3+4 | Steps 0 (oracle ladder) and 0.5 (equalise opportunities) are unchanged and remain gates. Step 2's "predicted iRT" gets its `a,b` from step 2 above. Mass-before-RT ordering is reinforced: tighter ppm buys more unique anchors. |
| **`doc/17`** MS1 arm | **Split cleanly in two.** The *scoring* path is per-candidate and streaming (~1 KB per live precursor, no persistent matrix) — that is kimi's and codex's answer and it stands. The *map* is run-level and 6–12 MB. They are different objects with different lifetimes and must not be conflated. |

**And a load-bearing point about `doc/17`'s kill gate:** the MS1 sub-score must
still pass its negative-control test before any MS1 *scoring* feature is built.
**But the map's anchor role does not depend on that gate.** Anchoring needs
correct positions; discrimination needs target-versus-decoy separation. A feature
finder can be excellent at the first and useless at the second — and given
`doc/17` §1, an MS1-envelope feature is *provably* useless at the second against
production decoys, since ODIA's decoys carry the target's precursor m/z.

## 7. Open questions for review

1. **Is §4's detectability bias fatal, or manageable?** This is the question. If
   fatal, the map can still index but must never anchor, and phases 3+4 need a
   different source.
2. **Identity assignment.** Only ~1.9% of MS1 features match a library precursor
   uniquely. Is uniqueness-at-tight-ppm the right criterion, or is there a better
   one — and does the answer change once MS2 evidence merges in at step 6?
3. **Is the merge at step 6 sound?** MS1 and MS2 features have different position
   semantics (precursor m/z vs fragment co-elution apex). What is the join, and
   what happens when they disagree?
4. **Precision is unmeasured.** "Unique" is not "correct". The pending test is
   whether those 1,794 pairs' library iRT correlates with observed RT. What
   correlation would be sufficient, decided **before** seeing it?
5. **Two instruments.** Astral has no ion mobility; S08 is 4D and
   `FeatureFinderCentroided` has no IM concept, so on diaPASEF it sees
   mobility-merged frames. Does the map need `Biosaur2Algorithm` for S08, and does
   a schema serving both compromise either?
6. **What breaks if the MS1 negative-control gate fails?** Which parts of this
   plan survive an MS1 sub-score that carries no information?

---

## 9. CiRT peptides: the seed that works (measured 2026-08-14)

Four anchoring mechanisms failed their controls in one afternoon — uniqueness,
mass-mode, whole-library Hough on MS1, and dscore-ranked Hough on MS2. Then:

| fit | n | line | r | median residual |
|---|---|---|---|---|
| **CiRT × DIA-NN RT** | **17** | **RT = 907.0 + 10.873·iRT** | **+0.9913** | **32.3 s** |
| CiRT × nearest MS1 feature | 26 | RT = 1119.7 + 2.718·iRT | +0.320 | 258.7 s |

**Seventeen endogenous peptides define the calibration line that 5M library
precursors could not.** No Biognosys spike-ins exist in this sample (0/11), so
the CiRT set (Parker et al., MCP 2015) — conserved, abundant, endogenous — is the
usable one.

Retrospect worth keeping: the whole-library MS2 Hough returned slope **+10.53,
intercept +890.6** against the truth **+10.873, +907.0**. The fit was essentially
RIGHT; its prominence (4.6× against a 4.0× permuted control) simply could not
certify it. The signal was there and the statistic could not see it.

### Why "a better MS1→sequence mapping rule" is the wrong question

Selection criteria compared on CiRT peptides with known truth, all drawing from
the same ±20 ppm pool:

| criterion | median err | within 60 s |
|---|---|---|
| nearest-m/z | 53.8 s | 54% |
| most-intense | 79.2 s | 46% |
| nearest predicted-RT | 53.8 s | 54% |
| both | 53.8 s | 54% |

All identical, because the **median candidate pool is 2 features (max 4)**. The
correct feature is usually ABSENT, which is the default map's 29.8% recall showing
through. No selection rule recovers a feature that was never detected.

### MS1-confirmed MS2 — real, controlled, and insufficient

An MS2 peak group is confirmed when an MS1 feature sits at the same precursor m/z
within ±W s of its apex (mid map, 95,216 features):

| set | n | slope | r | within 60 s of truth |
|---|---|---|---|---|
| all MS2 groups | 5,964 | +1.59 | 0.164 | 8% |
| confirmed ±10 s | 310 | +4.39 | **0.420** | 12% |
| **CONTROL, MS1 RTs +600 s** | 343 | +1.98 | **0.211** | 7% |

Confirmation lifts r 0.164 → 0.420 and the control collapses it to 0.211 — **the
first intervention of the day to beat its own control by a clear margin.** But
slope 4.39 against a true 10.873, and 12% within 60 s, is not a usable
calibration. Confirmation can only re-weight the population it is given, and an
uncalibrated pass supplies one that is ~98% misplaced.

### v4: seed from CiRT, then confirm

```
  1. CiRT targeted search      ~17 anchors, r ~ 0.99, no calibration needed
  2. predict library RT        from that line
  3. MS2 extraction            in a +/-60 s window instead of blind
  4. MS1 confirmation          NOW powerful: the population is mostly right
  5. robust refit on the union stop on SLOPE STABILITY, not ID count
```

The union is the right structure. It was being asked to bootstrap from nothing;
given a CiRT seed it consolidates a mostly-correct population instead.

**Open:** CiRT coverage on S08 (this is Astral-only), whether a targeted CiRT
search finds them without calibration (it should — they are abundant and few),
and how many of the ~120-peptide CiRT set are present rather than the 20 probed.
