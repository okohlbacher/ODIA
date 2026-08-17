# The RT mapping framework: model space, iRT space, run seconds

Written 2026-08-16 after the deep-research finding that our fine-tuning
normalises into a different space than the pretrained model lives in.

## The defect this replaces

`rtfinetune/run_exp.py:113-114`

    lo, hi = rt[tr].min(), rt[tr].max()
    tdf['rt_norm'] = (rt[tr] - lo) / (hi - lo)

Two compounding errors:

1. **Wrong space.** AlphaPeptDeep's `AlphaRTModel` trains on `rt_norm` produced
   by `rt / max_rt` with `min_rt` FORCED TO 0 (`_min_max_rt_norm = False` in
   alphabase `_normalize_rt`). We used a true min-max with a non-zero minimum
   (`rt_norm_min_minutes: 7.502` in the provenance), so the pretrained weights
   are asked to relearn a shifted AND rescaled target. That burns exactly the
   transfer-learning advantage the paper's 500-peptide result depends on
   (R^2 0.927 -> 0.986, Zeng et al., Nat Commun 13:7238).
2. **The space is sample-dependent.** `lo`/`hi` come from whichever 500-600
   peptides were randomly drawn, so every fine-tune defines a DIFFERENT target
   space, and the inverse mapping at deployment does not match the one used at
   training.

Symptom on record: the fine-tune scores well inside its own experiment
(held-out sd 0.4149 min) while the exported ONNX is no better than stock on real
data (p95 121.4 s vs 115.1 s).

## The three domains

| domain | symbol | units | who owns it |
|---|---|---|---|
| model space | `y` | [0,1] | AlphaPeptDeep. `y = t / R`, min pinned at 0. |
| run seconds | `t` | s | this run's gradient |
| calibrated run seconds | `t'` | s | after the per-run map `C` |

**iRT is deliberately NOT a pipeline domain.** It is a portable *storage* unit
only. The 11-peptide Biognosys rescale currently in `LibraryGenerator`
(`152.236*raw - 39.2322`, worst standard off by 8.76 iRT) is a THIRD transform
fitted on 11 spiked peptides; it is removed from the live path and kept, if at
all, only for writing portable libraries.

## The two chains, and the invariant that binds them

**Forward** (library -> extraction centre):

    y  --(x R)-->  t  --(C)-->  t'

**Backward** (observed apex -> fine-tuning target):

    t_obs  --(C^-1)-->  t  --(/ R)-->  y_target

`R` and `C` are the SAME objects in both directions, for the same run. This is
the whole framework; everything else is bookkeeping.

## Separation of concerns -- what learns what

- **`C` learns the run.** Global, peptide-INDEPENDENT gradient shape: column
  ageing, flow, temperature. Linear first, then monotone nonlinear.
- **The model learns the peptide.** Composition, length, modifications.
  Peptide-DEPENDENT structure that no global monotone map can express -- e.g.
  the +49.6 s cysteine bias, which survives an ORACLE monotone map (kimi
  measured the oracle floor at p95 122.1 s, in-sample and held-out).

This is why fine-tuning must target `C^-1(t_obs)/R` and not `t_obs/R`: training
on raw observed RT forces the model to absorb this run's calibration
idiosyncrasy, destroying portability. Training on the inverse-calibrated value
leaves `C` to do the per-run work and the model to do the chemistry.

## Invariants, to be asserted in code

- **I1.** `R` is fixed ONCE per run from the spectrum axis (max observed MS2 RT),
  never from a peptide subset, and is identical forward and backward.
- **I2.** `rt_norm = t / R`, minimum pinned to 0. No subset min-max anywhere.
- **I3.** `C` is STRICTLY monotone. PAVA plateaus and non-monotone Akima
  segments are illegal: enforce strict increase or fail closed.
- **I4.** `C^-1` is defined only on `[C(0), C(R)]`. Outside it returns a STATUS
  (`outside_support`), never a clamped endpoint -- clamping manufactures
  high-leverage training rows at exactly the gradient ends where anchors are
  thinnest.
- **I5.** Round trip `C^-1(C(t)) = t` within 0.1 s on a dense grid, asserted at
  runtime, every round.
- **I6.** ONE code path produces `rt_norm` for both training and inference.
  Divergence between them is the defect being fixed.

## Order of operations, per run

1. Predict `y` for the library with the current model.
2. Seed `t = y * R` (or the supplied affine, which must be expressed as a `C`).
3. Pass 1 extraction -> anchors `(peptide, t_obs)`.
4. Fit `C`: robust linear, then monotone nonlinear on the residual. Assert I3-I5.
5. Measure the post-`C` residual. Test it for PEPTIDE-DEPENDENT structure
   (composition regression, cysteine stratum). If none, stop -- fine-tuning has
   nothing to learn.
6. Fine-tune on `y_target = C^-1(t_obs) / R`, with I2 normalisation, the paper's
   hyperparameters (epoch 30, warmup 10, lr 1e-4, batch 256, L1) and
   `psm_num_to_test_rt_ccs > 0` so overfitting is actually monitored -- the
   shipped default of 0 means it is NOT, and our runs used the shipped default.
7. Re-predict, recompute `t' = C(y' * R)`, refit `C` once.
8. Pass 2 extraction.

## Acceptance -- measured in the RT domain against DIA-NN

Report on the same 12,308 DIA-NN-confident precursors, p95 |residual| and sd:

| stage | current | target |
|---|---|---|
| library via supplied map | 153.3 / 75.0 | -- |
| after `C` (linear+monotone) | 130.1 / 62.7 | <= 122 (the oracle-monotone floor) |
| after fine-tuning | not achieved | **< 122**, ideally toward DIA-NN's 65.5 / 31.9 |

The 122 s figure is kimi's measured oracle-monotone floor for THIS library: no
calibration can beat it, so any improvement below 122 s is necessarily the
model, which is the point of the exercise.

**STOP conditions**

- If the post-`C` residual shows NO peptide-dependent structure (step 5), do not
  fine-tune -- the deficit is rank error the model cannot fix from this data.
- If a correctly-normalised fine-tune does not beat 122 s p95 on a
  PROTEIN-level-split held-out set, the RT-parity programme stops and we record
  "RT: ~2x DIA-NN, predictor-bound".
- Never accept on anchor residual SD. Accept on identifications at 1% FDR with
  entrapment FDP not rising (doc/27).

## Known prior that tempers all of this

`doc/06` records that pasting DIA-NN's PERFECT RT column into our library moved
identifications by **-108 net**. A better RT axis reshuffles more than it gains.
The justification for this work is Cys recovery and correctness, not a large ID
win -- and that should be stated before the runs, not after.

---

# REVISION after adversarial review (codex + kimi), 2026-08-16

**Codex verdict on the design above: DO NOT IMPLEMENT AS WRITTEN.** Its
objections are accepted. What survives, what changes:

## Accepted objections

**A. `R` is arbitrary and CANCELS.** For `R' = aR` with `C'(x) = C(x/a)`, both
`C'(R'y) = C(Ry)` and `C'^-1(t)/R' = C^-1(t)/R`. So invariants I1/I2's fixation
on the scale constant was misplaced. **The real defect was never the scale --
it was that `lo`/`hi` came from the sampled SUBSET, so the train-time and
deploy-time conventions differed.** Any convention works if it is used
consistently in both directions and is run-independent.
Revised I1/I2: the normalisation convention is FIXED, declared, and identical
at training and inference. Not "must be max MS2 RT".

**B. Gauge freedom breaks the separation of concerns.** For any increasing `h`,
`f' = h(f)` with `C'(x) = C(R h^-1(x/R))` yields identical calibrated RT, so
refitting `C` after fine-tuning can UNDO any globally monotone model change.
"The model learns the peptide, C learns the run" is NOT identifiable on one run.
Revised: cross-fit `C` (never derive a row's target from a map fitted using that
row) and RESIDUALISE the fine-tune target against the same basis `C` can
express. Only the component `C` cannot represent is attributable to the model.

**C. The 11-standard rescale is a COORDINATE ADAPTER, not a run warp.** It
predicts the standard sequences and maps those predictions onto known iRT; it
never reads their observed run RT. My "third transform fitted on 11 spiked
peptides" framing was wrong -- removing it costs interoperability and gains no
accuracy. **iRT stays as the canonical serialisation gauge.** `y` is
checkpoint- and preprocessing-specific and is not a portable unit.

**D. The map family cannot honour the invariants as written.** PAVA produces
plateaus (so strict monotonicity is unimplementable without changing the
family), Akima can overshoot between monotone knots, `invertAt()` clamps rather
than returning status, and its bisection assumes global monotonicity from two
endpoints. These must be fixed BEFORE any framework rests on them.

## Also revised by kimi's RT-model review

**E. This is orchestration, not new C++.** All pieces exist: the pass-1 anchors
sidecar (`OpenDIAlyzer.cpp:1536`), `finetune_rt.py`, `export_finetuned_rt.sh`,
`-repredict_irt`. Cost measured at 65-415 s for 500-600 peptides.

**F. The 2x gap is a missing per-run fine-tune, not predictor weakness.** Stock
+ best monotone map floors at sd ~41 s / p95 85-96 s; WITH per-run fine-tuning
held-out p95 is 41-52 s, matching DIA-NN's own refit (42.5 s, same protocol).
DIA-NN refits per run -- 1.7.x runs 12 interleaved calibration iterations, 2.x
fine-tunes the predictor. Benchmarking our stock model against their refit model
was not like-for-like.

**G. CAM's ~0 response is constant-mod unidentifiability**, not an encoder bug:
CAM sat on every cysteine in pretraining, so its effect is unidentifiable
between the C token and the mod pathway and was absorbed into the token.
Explicit CAM is a NO-OP (not double-counting); a CAM-free cysteine is
unrepresentable, wrong by ~+9 s/Cys. That is only ~9 s of our measured +49.6 s,
so it is NOT the whole cysteine bias -- the remainder is unexplained.

**H. The frozen affine degrades every search.** `-irt_slope 11.399
-irt_intercept 915.0` was passed in every benchmark; the "supplied affine
153.3 s" figure is contaminated by it.

## Unresolved, and it gates a headline

`rtfinetune/onnx/` has NO provenance sidecar (unlike `integrated/` and
`check_direct/`). Kimi asserts it is ft-23k, which would make my
"stock vs fine-tuned" comparison two fine-tuned models. Against that: a model
fine-tuned on 23k peptides from THIS run would not show p95 115 s on those same
peptides. **Unresolved; resolve before quoting any stock-vs-tuned number.**

## Revised order of work

1. Resolve the `onnx/` provenance.
2. Drop the frozen affine from all benchmarks; re-baseline honestly.
3. Fix the fine-tune normalisation to a FIXED, declared, run-independent
   convention identical at train and deploy (the actual bug).
4. Fix the map family: strict monotonicity, non-overshooting interpolation,
   status-returning inverse.
5. Only then wire the per-run fine-tune between passes, with >=2,000 anchors,
   cross-fitted `C`, and the target residualised against `C`'s basis.

Acceptance is unchanged: identifications at 1% FDR with entrapment FDP not
rising. Never anchor residual SD.

---

# DESIGN PRINCIPLE, stated by the project owner 2026-08-16

**Fine-tuning is intended to TRAIN ON THE SAME RUN it is applied to. Neither a
fine-tuned model nor a recalibration is transferable between runs.**

This is not a caveat, it is the architecture. Both objects are per-run:

| object | fitted on | applied to | lifetime |
|---|---|---|---|
| calibration `C` | this run's anchors | this run | discarded with the run |
| fine-tuned RT model | this run's confident IDs | this run | discarded with the run |

Consequences, and a correction to the previous section:

1. **"Train-on-test" is the WRONG objection** to the rtfinetune numbers. Training
   on the run's own IDs and predicting for that run is the intended use. My
   framing of the 115.1 s figure as "optimistic because it may be train-on-test"
   is withdrawn.
2. **The v4/v5 record CONFIRMS the principle rather than contradicting it.**
   `-rt_model`'s help text: applying a model fine-tuned on one run to ANOTHER
   cost 2,027 confident precursors (v5 35,556 vs v4 37,583). That is a
   cross-run transfer failure -- exactly what this principle forbids -- not
   evidence against per-run fine-tuning. The same text already prescribes the
   right design: "Per-run fine-tuning belongs inside that run's calibration
   loop, built from it and discarded with it."
3. **The real validation requirement is narrower and still binding:** held-out
   PEPTIDES within the run, split at PROTEIN level. Peptides of one protein
   share composition and co-elute, so a peptide-level random split leaks. My
   fold split hashed the peptide sequence and is therefore optimistic for THAT
   reason -- not because of the run.
4. **Nothing fine-tuned may ever be written into a shipped library.** The
   library is a cross-run artefact and must carry stock predictions; the per-run
   model exists only between pass 1 and pass 2.
5. It follows that `-rt_model`'s default (the OpenMS bundled model) is CORRECT
   for library generation, and the fix is not to change that default but to add
   the per-run loop.

---

# REVISION 2: FDR is POST-calibration and POST-final-extraction ONLY

Stated by the project owner 2026-08-16. **No q-value may be consulted anywhere
in the calibration or fine-tuning path.** Anchors are selected on RELATIVE RT
ACCURACY from a WIDE window, incrementally.

## What this removes, and why it matters

The current harvest is q-gated -- `OpenDIAlyzer.cpp:1500`,
`if (g.decoy || g.qvalue > anchor_q) continue;` at `-anchor_q 0.05`. That single
line is the root of four separate problems already measured in this project:

1. **Selection bias.** Anchors are an FDR-ACCEPTED sample, so every statistic
   computed on them is conditioned on the outcome
   (`memory: odia-accepted-groups-are-a-biased-sample`; sizing a window from
   them cost half a run).
2. **The self-confirming loop.** A wrong prediction produces a false peak AT the
   prediction, that peak passes FDR, becomes an anchor, and teaches the refiner
   that no correction is needed. Measured: 76.5% of cysteine anchors were false.
3. **A ladder that is not even valid there.** Pass-1 entrapment FDP is 3.4-8.6%
   at nominal 1% -- the q-values being used to gate anchors are wrong by 3-9x.
4. **Instability.** The q-ladder sits on O(10) decoys; any perturbation swings
   identifications by thousands and the anchor set with them.

None of these can occur if no q-value is read before the final extraction.

## The replacement: incremental selection on relative RT accuracy

**Seed (prediction-INDEPENDENT).** Cross-charge agreement: the same modified
sequence observed at two or more charge states must elute at the same time.
Charge states are separate library entries extracted independently, so their
agreement is evidence about the RUN, not about the prediction. Peptides whose
charge states agree within a tolerance are the seed set. This is what stops the
loop -- the seed cannot be manufactured by a wrong prediction.

**Grow (relative, not absolute).** Given the current map, compute residuals for
ALL precursors' best-scoring candidate (best by spectral dscore -- a per-
candidate quantity, NOT a q-value). Accept those within `k * MAD` of the current
residual median. Refit. Repeat until the residual spread stops shrinking.

**Wide window throughout.** The extraction that feeds anchor selection uses a
DELIBERATELY wide window; narrowing happens only for the final extraction, after
calibration and fine-tuning are complete. A narrow window during calibration is
what censors the truth (measured: cysteine truth at +67.5 s median against a
+/-60 s window, so only 20.4% of false-picked anchors had a truth candidate
at all).

**Order:** wide extract -> seed by cross-charge -> fit -> grow by relative
residual -> refit -> converge -> fine-tune -> re-predict -> refit map -> FINAL
narrow extraction -> scoring -> **FDR, once, here**.

## Consequences for the rest of the plan

- `-anchor_q` is not used in the new path. Keep it only for the legacy mode.
- The >=2,000-anchor gate stays, but counts RT-consistent anchors, not
  FDR-accepted ones.
- The acceptance rule for the WHOLE loop is still identifications at 1% FDR with
  entrapment FDP held -- that is a POST-hoc measurement of the finished run, and
  is not read by anything inside the loop.
- Every earlier statement of the form "pass 1 identified N precursors at 1% FDR"
  becomes irrelevant to calibration; those numbers were never load-bearing for
  the map and should stop being quoted as if they were.
