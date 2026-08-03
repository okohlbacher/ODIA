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
12x, and **16.5% of a human library lands in the single output bin 0.85-0.90**
with essentially nothing above 0.90. Confirmed three ways: the decile slopes,
the library's output histogram, and synthetic ladders through the model
(`G14` 0.016 -> `G9L5` 0.876 -> `L14` 0.936, so five leucines spend 0.89 of the
range and the next nine spend 0.06).

Held-out residuals against observed RT on S08, minutes, protein-level split:

| | linear calibration | flexible monotone | fine-tuned (linear) |
|---|---:|---:|---:|
| sd | 1.015 | 0.678 | **0.329** |
| p95 | 2.171 | 1.419 | 0.688 |
| decile 9 | | 0.918 | 0.409 |
| decile 10 | | 1.045 | 0.518 |

**The decisive observation is not the headline ratio.** It is that after
fine-tuning, a flexible calibration buys *nothing* over a straight line (0.449
vs 0.449 at 500 peptides), where before it was worth 33%. The curvature was in
the model, and retraining removed it rather than compensating for it.

**And the thing that keeps this in proportion:** an earlier controlled ablation
pasted DIA-NN's RT column into our library and changed identifications by
**-108 precursors and +21 proteins**. A *perfect* RT column is worth
approximately nothing in identifications on this run. Whatever is built here is
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
- **Risk:** low. The training identifications come from a search that did not
  use the library being built, so there is no circularity.
- **What it buys:** a narrower extraction window in Phase 2 -- roughly 1.5x on
  this run -- and nothing measurable in identifications.

**Provenance becomes mandatory the moment this exists.** Two libraries built
with different models are no longer comparable, and nothing currently records
which model produced a library. This is already a backlog item; R-A promotes it
from tidy to load-bearing, and it ships with R-A or not at all.

### R-B. In-search iterative refinement -- DO NOT BUILD YET

Search, fine-tune on the precursors passing 1% FDR, re-search. This is what
DIA-NN, Spectronaut and MaxDIA do, so it is defensible on precedent alone.

**It is deferred for one reason: the FDR risk is unquantified in the
literature, and we can measure it more cheaply than we can argue about it.**
The entrapment machinery this project needs for Phase 3 anyway is exactly the
instrument: fine-tune on a first pass's identifications, re-search, and see
whether the entrapment-estimated FDR still matches the reported one. If it
does, build the loop. If it does not, we will have the measurement the field
appears to lack.

Building it before that measurement would mean shipping a loop whose failure
mode -- more identifications, quietly wrong FDR -- looks exactly like success.

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
  ranking, and three confident diagnoses on this project have already been
  overturned by a single ablation -- methionine excision, the RT gap, and the
  claim that a multi-chunk Arrow column could not be reproduced at test scale.
