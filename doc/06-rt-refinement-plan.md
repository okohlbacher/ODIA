# Refining retention-time prediction: what to build, and what not to

**Status.** Plan. Sections 1-3 are measured on this project's data or verified
from primary sources; section 4 is the proposal; section 5 is what would
falsify it.

**One-line summary.** Fine-tuning the RT model on a run's own identifications
removes the defect that no calibration can touch, at 31 seconds per run, with
no change to ODIA's C++ and no runtime dependency on torch. It should be built.
It should *not* be built as an in-search feedback loop, because the risk of
doing so is unquantified in the literature and cheap for us to measure first.

---

## 1. What is actually wrong, measured

The RT model's raw output saturates. Its top two deciles are compressed 5x and
12x, and **12.8% of the human library lands in the single output bin 0.85-0.90**
(14.6% above 0.85) with essentially nothing above 0.90. Over the *identified*
subset the figure is 16.5%, which is identification-biased and must not be
quoted as the library number. Confirmed three ways: the decile slopes,
the library's output histogram, and synthetic ladders through the model
(`G14` 0.016 -> `G9L5` 0.876 -> `L14` 0.936, so five leucines spend 0.89 of the
range and the next nine spend 0.06).

Held-out residuals against observed RT on IH1, minutes, protein-level split:

| n_train | isotonic baseline | AA-composition ridge | **fine-tuned** | wall clock |
|---:|---:|---:|---:|---:|
| 500 | 0.701 | 0.616 | **0.449** | 31 s |
| 2 000 | 0.691 | 0.577 | **0.395** | 78 s |
| 23 179 | 0.678 | 0.555 | **0.329** | 17 min |

The middle column is a control that the first version of this plan did not run
and a reviewer did: ridge regression on amino-acid composition, same splits,
same seed. **The cheap run-specific correction buys 12-19%; the fine-tune buys
36-52%.** The torch dependency earns its keep by a factor of ~3 over the obvious
free alternative, which was the attack most likely to kill R-A.

n = 2000 is the sweet spot the first draft missed: 85% of the full-data gain for
1/13 of the wall time. Note the confound in reading this as a data-quantity
curve -- batch size is 1024, so these are 40, 80 and 920 optimiser steps.

Per decile, fine-tuned against the isotonic baseline: decile 9 0.918 -> 0.409,
decile 10 1.045 -> 0.518.

**The decisive observation is not the headline ratio.** It is that after
fine-tuning, a flexible calibration buys *nothing* over a straight line (0.4494
vs 0.4493 at 500 peptides), where before it was worth 33%. Isotonic regression
is the most flexible monotone map that exists; if the gain were recalibration,
isotonic on the raw model would have captured it, and it does not -- it barely
moves across 46x more training data (0.701 -> 0.678). What is left after the
best possible monotone map is elution-*order* error, which only retraining can
touch.

The curvature is **reduced, not removed**: the synthetic leucine ladder still
compresses at the top after fine-tuning (0.132 of range against 0.060 baseline,
2.2x better rather than gone), and the fine-tuned tail deciles remain 1.6x the
overall residual. That heteroscedasticity is what R-A' below exploits.

**And the thing that keeps this in proportion -- stated more carefully than in
the first draft:** an earlier controlled ablation pasted DIA-NN's RT column into
our library and changed identifications by **-108 precursors** net. That net
hides the actual effect. Decomposed: **3,703 precursors lost, 3,595 gained,
7,298 changed -- 20.4% of the answer** -- with lost and gained both concentrated
at 0.42-0.45x the shared median abundance, and roughly half of each
uncorroborated by DIA-NN's own run. A perfect RT column is worth approximately
nothing in the *count* and reshuffles a fifth of the *result*. Whatever is built here is
justified by extraction cost, not by sensitivity, and the plan says so
throughout.

---

## 2. What the three tools actually do

Verified from primary sources; where a source could not be found, that is
stated rather than filled in.

**AlphaPeptDeep** is the only one of the three doing true weight-updating
transfer learning with published hyperparameters. `ModelManager.train_rt_model(psm_df)`,
or `peptdeep transfer <settings.yaml>` from the CLI, fine-tunes the whole model
in place -- no layer freezing -- and saves `rt.pth` to
`{PEPTDEEP_HOME}/refined_models`, reused by pointing `model_mgr:external_rt_model`
at it. Shipped defaults have *diverged from the paper*: `default_settings.yaml`
sets 40 epochs / 10 warmup / batch 1024 / lr 1e-4, where Zeng et al. state
30 / 10 / 256. Input must be an AlphaBase `psm_df` with `sequence`, `mods`,
`mod_sites` and a normalised **`rt_norm`** column; raw retention time will not
work.

**DeepLC** separates the two operations explicitly. Calibration is post-hoc and
first-class -- `from deeplc import predict_and_calibrate`, default
`SplineTransformerCalibration` -- and is *a monotone mapping, so it cannot
change predicted elution order*. That is precisely the limitation we measured
independently, and it means DeepLC's headline feature cannot fix our defect.
Weight updating lives in a separate package, `DeepLCRetrainer.retrain()`,
defaults 75 epochs / batch 128 / all layers trainable, with opt-in layer
freezing keyed on the layer named `concatenate`.

**DIA-NN could not be verified.** Its core is closed -- `diann.cpp` is not in
the public repository -- so nothing about `--tune-rt`, `--tune-lib`, the
`rt.d{0,1,2}.pt` files, `--rt-model`, `--tokens`, the layers retrained, the
learning rate or the epoch count was confirmed by any source. The flags do
exist and run, which is direct evidence they work; what they do internally is
not knowable from here. **No claim about DIA-NN's fine-tuning mechanism should
enter this project's documents as fact.**

**Iterative refinement is standard practice.** DIA-NN, Spectronaut and MaxDIA
all refit RT calibration during the search on the precursors passing 1% FDR in
the current round. So a feedback loop is not novel and not disreputable.

**But its risk is unquantified.** Exactly one published statement of the hazard
surfaced -- that FDR estimated on a *subset* of targets is biased low, so the
anchor set used to refit carries more false identifications than 1% implies --
and it is an asserted mechanism, not a measurement. **No primary source
quantifies FDR inflation caused by refitting an RT model on a search's own
confident identifications.** That absence is the single most important input to
this plan.

---

## 3. What OpenMS gives us, and what it does not

- **No DeepLC interface exists** anywhere in OpenMS.
- The **PeptDeep binding is inference-only** -- no `train`, no optimiser, no
  gradients -- **is not installed** (compiled into `libOpenMS.so`, headers never
  exported), and **rejects modified peptides**. It cannot be used for any of
  this, which is why ODIA drives ONNX Runtime itself.
- OpenMS's own sequence-based RT predictors were **removed** (`RTModel`,
  `RTPredict`, `SVMWrapper`). What remains are alignment tools, which perform
  the monotone remapping we have measured to be insufficient.
- **`tools/scripts/export_peptdeep_models_to_onnx.py` is the bridge, and it
  already exists.** It loads `rt.pth` through `AlphaRTModel` and exports
  `peptdeep_rt_dynamic.onnx` with input names `input_sequences` / `mod_x` and
  output `rt_pred` -- exactly what ODIA consumes today.

**Consequence, and it is the reason this plan is cheap:** a fine-tuned model
reaches ODIA as an ONNX file it can already read. **No C++ changes are
required.** Torch is confined to an offline preparation step and never enters
ODIA's runtime, which was the main objection to this route.

---

## 4. The proposal

### R-A. Offline fine-tuning as a library-preparation step -- BUILD THIS

A script that takes identifications with observed retention times, fine-tunes
`rt.pth`, and exports ONNX through OpenMS's own exporter. The library is then
generated with `-rt_model <the new file>`, which ODIA already supports.

- **Measured gain:** held-out sd 0.678 -> 0.449 min from 500 peptides in 31 s
  of CPU, or 0.329 min from the full set in 17 min.
- **Cost:** a separate venv with peptdeep and torch, plus a script. No C++.
- **Risk: NOT low, and the first draft was wrong about why.** It claimed "the
  training identifications come from a search that did not use the library being
  built, so there is no circularity". That is true only of the historical
  accident that these anchors came from a DIA-NN search. R-A is specified as
  per-run -- it learns *this* chromatography -- and per-run identifications on
  the run about to be searched means a first-pass search of that run. **R-A is
  round one of R-B.** The FDR caveat below belongs here, not there.

  What lowers the risk, and is evidence rather than argument: the fine-tuned
  model generalises to held-out *sequences* as well as to held-out *proteins*
  (0.3287 peptide-split against 0.3312 protein-split). Decoys are permuted
  target sequences, i.e. held-out sequences. That is not proof of target-decoy
  symmetry, but it is the strongest available prior and it points away from the
  hazard.
- **What it buys:** see section 4a. The first draft said "roughly 1.5x", which
  was lifted from the perfect-RT ablation's window ratio and does not belong to
  the fine-tune at all.

**Provenance, scoped.** The first draft said it "ships with R-A or not at all",
which contradicted the same paragraph's "no C++" cost. The hazard is applying a
run-specific model to the wrong run, and provenance records that rather than
preventing it. The minimum that is genuinely load-bearing is a sidecar: model
file SHA-256, source path, and the identification set it was tuned from. A day,
not a project. Note also that the worst version of this hazard is *already*
closed: `LibraryGenerator::fitIrtCalibration` refits the iRT constants at run
time from the model itself, so swapping the checkpoint self-recalibrates.

**Why provenance matters at all:** Two libraries built
with different models are no longer comparable, and nothing currently records
which model produced a library. This is already a backlog item; R-A promotes it
from tidy to load-bearing, and it ships with R-A or not at all.

### 4a. What a narrower window is actually worth

The first draft asserted this. The arithmetic, from the run's own metadata
(24 windows per cycle, 1,342 spectra per window, cycle 0.0229 min) and the
measured Phase 2 slice (113 ns per transition probe):

| whole library x whole run | points | memory | match | decode after the reader fix |
|---|---:|---:|---:|---:|
| no RT window | 1.7e10 | 137 GB | ~32 min | 3.6 s |
| raw-model window (14% of run), targets + decoys | **4.9e9** | 39 GB | 552 s | 3.6 s |
| fine-tuned window (6.9%), targets + decoys | 2.4e9 | 19 GB | 268 s | 3.6 s |

Two consequences the first draft missed.

**`ChromatogramExtractor` throws above 2^32 points.** At the raw model's
accuracy a whole-run, target-and-decoy extraction is 4.9e9 and **the extractor
refuses to run**; with the fine-tuned window it is 2.4e9 and does. RT accuracy
is not a constant factor on Phase 2, it is a feasibility threshold. That is the
argument R-A should make.

**Once decode is fixed, window width becomes the only thing that matters.**
Match against decode goes from 1:1300 today to roughly 75:1 after the batched
reader. While decode stays slow, window width is worth nothing in time, because
the one-forward-pass design already amortises decode across every transition.

**And "no C++ changes are required" was true of the wrong thing.** It is true
that a fine-tuned ONNX drops into `-rt_model`. It is false that the benefit
lands without C++: `ChromatogramExtractor` has one global `rt_low`/`rt_high`,
sizes every transition for the whole range, and **never reads a retention time
from the library at all**. It is RT-blind today. The per-precursor window is the
deliverable; the model is an input to it.

### R-A'. Predict an interval, not a point -- BUILD THIS TOO

The residual is strongly heteroscedastic: fine-tuned deciles 9 and 10 are 0.409
and 0.518 against 0.329 overall. A single global window must be sized for the
worst decile, so the whole library pays the tail. A per-precursor window sized
from predicted uncertainty -- a quantile head, or just a lookup over the classes
already identified as bad (length 25-30, late RT, high GRAVY, charge 3) --
shrinks the mean window by roughly 0.64x *on top of* the accuracy gain.

Given that window width is the difference between 19 GB and a `uint32`
overflow, this is plausibly worth more than the accuracy improvement it rides
on, and it costs nothing at inference.

### R-B. In-search iterative refinement -- NOT DEFERRED, UNBUILDABLE

Search, fine-tune on the precursors passing 1% FDR, re-search. This is what
DIA-NN, Spectronaut and MaxDIA do, so it is defensible on precedent alone.

**The honest reason it is not being built is that it cannot be: "search,
fine-tune on the 1% FDR precursors, re-search" requires a search engine and an
FDR estimate, and ODIA has neither yet.** The first draft dressed a roadmap fact
as an FDR argument. The FDR argument is real, but it applies to R-A, where it
has been moved.

The measurement still needs making, and the entrapment machinery Phase 3 needs
anyway is the instrument: fine-tune on a first pass, re-search, and check
whether the entrapment-estimated FDR still matches the reported one. That is a
precondition for R-A fed by ODIA's own first pass -- not for R-B, which is
blocked on prior work regardless. Until then R-A is fed by an external search,
and a workflow whose first step is "run DIA-NN" is itself a problem this project
has already named.

### Closed questions -- recorded so they stay closed

- **Per-charge RT models.** Charge states of one peptide coelute to within a
  median of **0.0018 min** (p95 0.0355) across 3,800 multi-charge peptides.
  Against a 0.329 min residual the ceiling is ~10% and the realistic value is
  zero. Treat as an upper bound: DIA-NN may propagate RT within an elution
  group.
- **Calibrating against a spiked iRT standard kit.** Eleven anchor points can
  only fit a monotone map, and the best possible monotone map over 23,179 points
  still leaves 0.678 min. Eleven cannot beat 23,179.
- **MS1 for RT anchoring.** Circular for narrowing a window -- you need a window
  to extract the MS1 XIC in. Useful for refining an apex after a candidate hit,
  which is Phase 3, not here.

### R-C. DeepLC as an alternative predictor -- EVALUATE, CHEAPLY, FIRST

DeepLC ships four `.pt` files and depends on torch; it has no ONNX export, so
native integration would mean exporting the models ourselves and writing an
encoder for its atom-composition input scheme, with an independent oracle to
validate it. That is a substantial job.

**Before any of it: test whether DeepLC's raw predictions saturate the way
PeptDeep's do.** Its headline feature is a monotone calibration, which we have
measured cannot fix saturation, so DeepLC is only interesting if its untuned
output is better shaped. That is an afternoon in a scratch venv and it gates
everything else.

### R-D. Do nothing -- the honest baseline

Justified by the ablation: a perfect RT column is worth -108 precursors. If
Phase 2 profiling shows extraction width is not a leading cost, R-A's benefit
is close to zero and this is the correct choice. R-A is cheap enough that it is
worth doing anyway, but it must not be defended as a sensitivity improvement.

---

## 5. What would falsify this plan

- **R-A**: Phase 2 profiling shows extraction width is not a leading cost term.
  Then R-A buys nothing measurable and R-D is right.
- **R-A**: a model fine-tuned on one run degrades another run's library. Expected
  by construction -- it learns *this* chromatography -- which is why the model
  must be per-run and recorded, never shipped as a default.
- **R-B**: entrapment shows the feedback loop leaves FDR intact. Then it should
  be built, and the deferral cost us time.
- **R-C**: DeepLC's raw output turns out not to saturate. Then it is worth the
  integration and R-A is the weaker option.
- **The whole plan**: every number here is from **one run, one gradient, one
  organism, essentially one modification**. A second run could change the
  ranking, and four confident diagnoses on this project have now been overturned
  by a single measurement -- methionine excision, the RT gap, the multi-chunk
  Arrow claim, and this plan's own "there is no circularity".

- **IH1 is close to the least RT-sensitive configuration available**, and the
  plan's central "RT does not matter" rests on it. Ion mobility is present and
  is doing discrimination RT would otherwise do; the gradient is short; the
  windows are narrow; the sample is high-load human tryptic. **The counter-test
  is already on disk**: `data/astral.mzpeak` is an Orbitrap Astral run with no
  ion mobility and +/-2 Th windows -- 3D, not 4D -- which is precisely the
  regime where RT should matter more. `/scratch/kohlbach/rawd/` also holds IH3
  and IH2 from the same batch and method, which test run-to-run transfer of a
  fine-tuned model. These caveats are not limitations, they are unexecuted
  experiments.

- **RT prediction quality matters only in library-free mode**, because a
  DDA-derived library carries measured retention times. ODIA is library-free.
  That is a point in R-A's favour the first draft left on the table.

- **`prep.py` supports five modifications.** Its UniMod dictionary covers
  {4, 35, 1, 21, 7} and raises on anything else. That is the right failure mode,
  but it is a real limit on R-A as it stands, and the fix is to reuse ODIA's own
  verified 109-element encoder rather than a second hand-written parser.
