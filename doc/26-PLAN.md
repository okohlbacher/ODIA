# ODIA: the plan of record

2026-08-15. **Supersedes doc/15–25**, which remain as the working record. Where
this contradicts them, this wins. FDR/FDP is the LAST phase and is out of scope
until the front of the pipeline is aligned.

---

## PART A — ESTABLISHED (measured, and survived adversarial review)

| # | fact | evidence |
|---|---|---|
| A1 | ODIA identifies **3,625** precursors on the full 4,991,901-target Astral library, entrapment **FDP 5.98%** at nominal 1% | `fullpass1`, `ref_astral` |
| A2 | **`-pass1_precursors` is a STRIDE**, not a cap. At 100,000 it searched ~101 of 10,040 findable precursors. Removing it: **0 → 3,625 IDs** | `doc/22` |
| A3 | q-floor: `q_min = π₀·(Ntar/Ndec)/targets_above_top_decoy`, so **~100 clean targets** is the minimum for any 1% claim. Predicted 0.0361 vs 0.0365 observed | `lda.h:209-265` |
| A4 | **Candidate flow** (10,040 denominator, mutually exclusive): 0.2% no candidate · **11.7% truth outside support** · **12.7% reachable but not picked** · 1.7% ranking · 73.7% correct | codex-corrected |
| A5 | **Effective extraction support is ±112 s**, measured — not the ±60 assumed | p99.9 of \|cand−pred\| |
| A6 | The **CiRT line is WORSE than raw library iRT** (p50 42.0 vs 36.7 s; decile bias to −102 s) | kimi |
| A7 | **Monotone PAVA removes ALL systematic bias** (decile medians → −0…+6 s, median signed +0.1 s) but leaves p50 32 s | held-out |
| A8 | Trim-percentile sweep (50/67/95/99/none) is **flat** — binned medians + PAVA are already robust | held-out |
| A9 | **Three RT levels exist; level 3 was never invoked.** Frozen: library iRT 88.89% → ridge 93.62% → **fine-tuned 97.27%** at ±60 s | `doc/14`, run scripts |
| A10 | The **oracle ladder works**: 6,815 IDs on DIA-NN's own confident set — extractor, picker, gate, classifier, refiner and FDR all function | `doc/21` |
| A11 | **DIA-NN uses no π₀**; ours defaults **off** (`lda.h:100`) | survey + code |
| A12 | Entrapment matching used starts-with against `sp\|ENTRAP_…\|` — **14.7% of the library was invisible**. Fixed | `5ccf26e` |
| A13 | GBT `min_child_rows` was never enforced per child (V1's defect, same variable). Fixed | `34930e9` |

## PART B — REFUTED (do not re-propose; all of these were mine)

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
The RT window cap: `pass2_window = min(cap, max(floor, 2×p95))` silently discarded
a measured 128 s in favour of 60 s. **Do not touch the m/z arm** — its default is
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
