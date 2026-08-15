# ODIA: the plan of record

**v2, 2026-08-15**, after review round 13 (codex effort max, kimi) said *do not
implement v1 as written*. Both found over-corrections; kimi found that v1's
motivating number for D0 was itself the contaminated statistic D3 disqualifies.
Supersedes doc/15–25, which remain as the working record. Where
this contradicts them, this wins. FDR/FDP is the LAST phase and is out of scope
until the front of the pipeline is aligned.

---

## PART A — ESTABLISHED (measured, and survived adversarial review)

| # | fact | evidence |
|---|---|---|
| A1 | **NARROWED:** 3,625 precursors is an exact-run fact. The **FDP is not settled**: it is a hand computation, it used protein-level r = 0.19592 where the tool's precursor-level r = 0.172 gives **6.80%** on the same 42 hits, it carries **no Poisson interval**, the run **predates the A13 GBT fix**, and "4,991,901 targets" counts 733,780 entrapment precursors as targets | needs a rerun |
| A2 | **`-pass1_precursors` is a STRIDE**, not a cap. At 100,000 it searched ~101 of 10,040 findable precursors. Removing it: **0 → 3,625 IDs** | `doc/22` |
| A3 | q-floor: `q_min = π₀·(Ntar/Ndec)/targets_above_top_decoy`, so **~100 clean targets** is the minimum for any 1% claim. Predicted 0.0361 vs 0.0365 observed | `lda.h:209-265` |
| A4 | **Candidate flow** (10,040 denominator, mutually exclusive): 0.2% no candidate · **11.7% truth outside support** · **12.7% reachable but not picked** · 1.7% ranking · 73.7% correct | codex-corrected |
| ~~A5~~ | **WITHDRAWN → PART C.** ±112 s is a quantile of *selected candidates*, not the extraction support, and the statistic appears in no log, script or note — "a number with a citation to itself" | — |
| A6 | The **CiRT line is WORSE than raw library iRT** (p50 42.0 vs 36.7 s; decile bias to −102 s) | kimi |
| A7 | **NARROWED:** PAVA removes the measured **decile-level LOCATION bias** (range −4…+6 s, median signed +0.1 s). It does *not* establish removal of bias across charge, modification, intensity or within-decile structure. p50 stays 32 s | held-out |
| A8 | Trim-percentile sweep (50/67/95/99/none) is **flat** — binned medians + PAVA are already robust | held-out |
| A9 | **CORRECTED:** level 3 **has** been invoked — `doc/14:251` records an end-to-end Astral run on 2026-08-09 with `-repredict_irt -rt_model <tuned>`: **6,382 → 6,836 IDs, window 213.2 → 134.1 s**. "Never invoked" is true only of the recent parity campaign. **And the 97.27% must NOT be used for sizing** — `doc/14:160-166` forbids it: measured through a non-production alignment, and the peptdeep row includes training peptides (held-out by stripped sequence: 81.15% vs 88.24% at ±30) | `doc/14` |
| A10 | The **oracle ladder works**: 6,815 IDs on DIA-NN's own confident set — extractor, picker, gate, classifier, refiner and FDR all function | `doc/21` |
| A11 | **VERSION-SCOPED:** DIA-NN **1.7.x** uses no π₀ (2.x is closed and not represented); ours defaults off. The stale `lda.h:101` comment claiming "DIA-NN parity" is now **fixed in code**, not only in docs | survey + code |
| A12 | Entrapment matching used starts-with against `sp\|ENTRAP_…\|` — **14.7% of the library was invisible**. Fixed | `5ccf26e` |
| A13 | GBT `min_child_rows` was never enforced per child (V1's defect, same variable). Fixed | `34930e9` |

## PART B — REFUTED, after v1's OVER-CORRECTIONS were reversed

**I asked both reviewers whether I had over-corrected. Both said yes, to four
entries. Restored below.**

| v1 said REFUTED | v2 verdict |
|---|---|
| "no monotone function of iRT can fix the 32 s" | **RESTORED, reworded.** `doc/14:85-90` measures 38.55 s as the *ceiling* for any monotone post-hoc map, from an independent implementation. Fine-tuning is not a counterexample: it changes peptide *ordering*, so it is not a monotone remapping of the original iRT. Correct form: **"a monotone remapping of the fixed stock predictor cannot remove its rank errors."** |
| "availability is window-limited" | **RESTORED.** If centre error puts truth outside support, availability *is* operationally window-limited; calibration and widening are alternative remedies. The split per axis remains unverified. |
| "the 3-candidate cap is not binding" | **MOVED — but to ESTABLISHED, not unverified.** A K sweep already exists: `BACKLOG.md:2834`, S08 K = 12/6/3/1 → **1,307 / 1,245 / 1,274 / 1,141**. K=3 ≈ K=12. My planned K replay was redundant. The 54% statistic still cannot carry the conclusion; the sweep can. |
| "envelope-only MS1 features are identical for target and decoy" | **RESTORED.** The refutation inverted codex's position — it said envelope-*alone* features are nondiscriminative; co-elution features are not envelope-only. Premise is code-verified (`DIANNLibraryFile.cpp:757-763`). |
| three untraceable numbers | stay dead, **wrong reason**: unsupported, not disproven. Reclassify as *withdrawn/untraceable* and flag `doc/24` §2.1/§2.4. |

**Still refuted:** "the blocker is the classifier prior" and "competition/null
width" — **but only as explanations of the A2 zero.** As live suspects for the
remaining FDP gap they are open, and the crossed design that would separate them
has never been run.

## PART B′ — genuinely refuted (all mine)

| claim | why it died |
|---|---|
| "the blocker is the classifier prior" | explained an artefact of A2 |
| "the blocker is competition / null width" | same artefact, second explanation |
| "V2 applies Storey π₀" | `use_pi0 = false`, and off in the run |
| "no monotone function of iRT can fix the 32 s" | a substantial share is stock-peptdeep error, which fine-tuning addresses |
| "availability is window-limited" | it is **half centre error** (A4); confounded by A6 |
| "the 3-candidate cap is not binding" | 54% of misses **at** K=3 is evidence **for** censoring |
| "m/z recalibration is dead wiring" | **`-mass_width_from_ids` defaults to `measure` deliberately**; `apply` was measured to take 4,969 → 2,422 |
| "envelope-only MS1 features are identical for target and decoy" | candidate features depend on fragment-derived apexes; refuted by codex |
| 174 B/precursor; 41.0→84.1; 40.7→75.7 | untraceable to any source; two collide with other measurements |

## PART C — UNVERIFIED (do not build on these)

- The **65.5% vs 88.89%** discrepancy for "library iRT at ±60 s" on the same instrument.
- Reachability on the **curve** and **fine-tuned** axes. The curve I have is the raw-axis curve and **cannot size a window**.
- Whether the 12.7% "reachable but not picked" is cap-censoring or picker failure.
- **Why the CiRT seed degrades the axis** (A6). A seed that makes things worse poisons everything hung off it.
- Everything on **S08**. All of the above is Astral.

---

## PART D — THE PLAN

### D0 · Fix the one genuine wiring defect
**v1 cited the wrong number.** `ref_astral.log:96` reads `pass 2 extraction window
60 s (2 x p95 47.9551 s, floor 20 s, cap 60 s)` — the discarded measurement was
**95.9 s, not 128 s.** The 128 s came from the refiner's *"window for 99%
coverage"* line, which is precisely the anchor-contaminated statistic D3
disqualifies (anchor p99 455 s vs true 108 s). **v1 motivated its first
implementation step with a quantity it elsewhere calls contamination.**

The defect is still real — `BACKLOG.md:2636` records 300/600/1200 s giving
bit-identical output — but it is "the cap discards a 95.9 s measurement", and the
replacement must come from D3's core-scale estimator, not from either number. **Do not touch the m/z arm** — its default is
deliberate and its estimator, not its wiring, is what needs work.

### D1 · Answer PART C before designing further
Four cheap measurements, no new machinery:
1. reconcile 65.5% vs 88.89% (grep the artefact, state the axis and the n);
2. re-measure reachability on the curve axis;
3. **K = 3/5/10/all replay on identical chromatograms** — the only test that splits A4's 12.7%;
4. audit why the CiRT seed degrades the axis.

### D2 · The recalibration order
```
  1  m/z          once, no loop. Synergistic: tighter ppm raises anchor
                  uniqueness 196 → 365, expanding and purifying the RT anchors.
                  Estimator must use the UNCENSORED pass-1 harvest, not accepted
                  IDs -- that is why `apply` failed.
  2  linear seed  robust (RANSAC/jackknife), NOT OLS on contaminated anchors
  3  monotone map binned medians + PAVA + akima
  4  FINE-TUNE    after the map: labels are observed RTs expressed through the
                  current map, so a stale axis would be baked in
  5  REFIT map    mandatory -- peptdeep emits iRT, a RANK not seconds
  6  window + ONE final extraction
```
**Never carried across a predictor change:** the window, and any residual statistic.

### D3 · Window derivation, without truth
Anchor p99 reads **455 s against a true 108 s** — high quantiles measure
misassignment. Therefore:
1. **Fit the core, not the tail** — MAD or p25–p75 IQR; DIA-NN uses the **80th
   percentile**, deliberately below the contamination tail.
2. **Use the uncensored pass-1 harvest** (1e9 s window), never pass-2 anchors.
3. **Edge pile-up diagnostic** — fraction of apexes within ε of the boundary;
   a runtime test for "the window is binding", requiring no truth.
4. **Per-bin, not global** — error varies systematically along iRT.
5. Validate on precursors the feature finder did **not** detect.

With the fine-tuned axis (p95 45.3 s), **±90 s buys ~99% core coverage** — which
is why D2 step 4 must precede any sizing.

### D4 · Convergence
**Convergence quantities must be functions of model PARAMETERS, or of predictions
on a set frozen before round 1 — never of the fitting sample.** All four
statistics that failed during the original RT work were functions of a changing,
self-selected population.

### D5 · Lossless final extraction
**Narrow primary path + fixed wide rescue path.** Unsupported strata keep wide
windows. A single narrow window is never right.

---

## PART E — HOW THIS GOES WRONG AGAIN

The failure mode of this project, measured across ~12 review rounds: **acting
before verifying the premise.** Four flag-semantics traps (`-pass1_precursors`,
`-rt_window` × 2, `-mass_width_from_ids`), two committed conclusions explaining an
artefact, and one review dispatched on a false premise that made both reviewers
reason from my error.

**Before any run: prove the thing you believe ran, ran.** Compare extracted-count
invariants between configurations; identical counts mean the knob did nothing.
Read the option's registration text before assuming its semantics.
