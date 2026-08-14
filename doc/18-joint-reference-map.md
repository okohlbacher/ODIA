# The joint reference map: one anchor structure for calibration and fine-tuning

Status: PLAN v1, for adversarial review. Written 2026-08-14.

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

**(b) It is the exogenous transform `doc/16` §2.2 said Option B required.** Both
codex and kimi independently found that peptdeep emits **iRT — a rank, not run
seconds** — so `RT = a + b·iRT` needs `a` and `b` from something available before
any identification, or Option B is as circular as the thing it replaces.
Measured: **1,794 unique `(library iRT, observed RT)` pairs** at ±5 ppm from the
mid-sensitivity map, spread evenly over all ten RT deciles — more than the frozen
RT work's ~1,183–1,408 anchors, and unconditioned where those were not.

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

## 4. THE SHARPEST RISK, stated first

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

## 5. Bootstrap order — the circularity, and the way out

MS2 features come from candidate detection, which needs RT windows, which needs
calibration, which needs anchors. The escape is that **the MS1 side is
library-independent and needs no calibration to exist**:

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
