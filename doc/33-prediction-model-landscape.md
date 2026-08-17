# Better RT, CCS and MS2 models? A researched answer

Commissioned after doc/32 established that our fragment-selection bias is a
model training-regime gap rather than a plumbing bug. 101 agents, five search
angles, three-vote adversarial verification per claim. **10 claims survived, 12
were killed** — and several of the killed ones are the interesting part.

Confidence and vote counts below are the harness's, not mine.

## The one that matters most: a timsTOF MS2 model already exists

**HIGH, 3-0.** `Prosit_2023_intensity_timsTOF` was fine-tuned on **277,779
timsTOF spectra**, and it exists precisely because Orbitrap-trained HCD Prosit
2020 gives materially lower spectral contrast angle on timsTOF data (Adams et
al., Nat Commun 15, 2024). `ms2pip_timsTOF2024` is a second option. Rescoring
pipelines route by instrument as a matter of course.

That is a direct, named fix for the bias measured in doc/32 — and it confirms
the diagnosis was right about the *kind* of problem even though it was wrong
about the cause.

## And there is no Astral model, from anyone

**HIGH, 3-0.** MSBooster's Koina catalogue (updated 2026-01-14) lists MS2:
AlphaPept_ms2_generic, Prosit_2019/2020/2023/2024/2025 variants, four ms2pip
variants, UniSpec, PredFull. Ion mobility: exactly two — AlphaPept_ccs_generic
and IM2Deep. **Astral appears in no section, for RT, MS2 or IM.** The gap is
chronologically forced: the systematic benchmarks draw on PRIDE submissions that
predate the Astral's 2023 release.

For our new neat/Astral benchmark that means there is no instrument-matched
option, and QE/Lumos conditioning is the honest best available — defensible
because Astral MS2 is HCD in an Orbitrap-family instrument, and Prosit HCD
models trained on QE/Lumos are in routine published use for Astral DIA library
generation. Note the asymmetry the verifiers insisted on: the
"no evidence for this instrument" argument is **strong and demonstrated for
timsTOF, formally-true-but-weak for Astral**.

## Koina: capped at 1,000 peptides a request — but probably runnable locally

**HIGH, 3-0.** Verified against Koina's own Triton `config.pbtxt` files, not the
client docs: `max_batch_size: 1000` appears uniformly across
Prosit_2019_intensity, Prosit_2023_intensity_timsTOF_core, AlphaPept_ms2_generic,
AlphaPept_rt_generic, AlphaPept_ccs_generic, Chronologer_RT, UniSpec,
ms2pip_2021_HCD, Deeplc_hela_hf and IM2Deep. Triton hard-rejects anything larger.
**5,000,000 / 1,000 = 5,000 sequential round-trips per library build.** Remote
inference is not viable for us, and it would ship our sequences to a public
server.

**But the claim that Koina is remote-ONLY was killed 0-3.** Its model
definitions — `config.pbtxt` plus backends, including a `Chronologer_RT` entry —
live in a public GitHub repo, so a **local Triton deployment, and possibly
extraction of local model files, is plausible**. The verifiers flagged this as a
lead to check rather than an established finding. It is the most promising route
to the timsTOF Prosit weights, and it is the first thing to test.

## Retention time: accuracy and fine-tunability point at different models

**Chronologer** (MEDIUM, 3-0) is the strongest reported accuracy candidate,
Apache-2.0, with pretrained weights and the full 2.2M-peptide training DB shipped
in-repo. Caveat: the headline numbers are a 2023 unreviewed self-benchmark
against non-fine-tuned competitors.

But its adaptation story is awkward for us (**HIGH**, 3-0 on data efficiency,
2-1 on forgetting, with the forgetting half corroborated by shipped source):
transfer learning needs only ~100 peptides, yet training on a run's own IDs alone
causes **catastrophic forgetting**, so every refinement must replay the 72 MB /
2.2M-peptide Chronologer-DB. Its designed per-run mechanism is **alignment into
Hydrophobic Index space, not fine-tuning**.

**DeepLC** (HIGH, 3-0 across three claims) has the best-documented, peer-reviewed,
reusable per-run fine-tuning story of any RT model reviewed, with concrete data
volumes that our per-run ID counts clear by orders of magnitude.

Given doc/32's finding that ~35% of our RT error is monotone calibration shape
rather than rank error, note that Chronologer's *alignment* mechanism is
addressing the same thing we would get from a monotone recalibration.

## CCS is the weakest leg

**MEDIUM, 3-0 on existence; no accuracy claim survived.** Only two servable
models were found — AlphaPept_ccs_generic (what we already use) and IM2Deep — and
**no benchmark numbers, no fine-tuning API and no timsTOF-vs-generic comparison
survived verification** for either. There is no evidence base here to act on.

## Prosit_2025_intensity_MultiFrag is not the upgrade it sounds like

**HIGH, 3-0.** Nature Methods s41592-026-03042-9 (23 Mar 2026): a modified-GRU
Prosit predicting 815 fragment ions across length, charge and fragment type. It
is **MS2-intensity-only and Orbitrap-only**. It cannot replace AlphaPeptDeep for
RT or CCS and does not close the instrument gap; its contribution is alternative
fragmentation modes.

## The Prosit-over-AlphaPeptDeep result is weaker than it looks

**MEDIUM, 3-0.** PepSpecBench (arXiv:2605.01945, May 2026, third-party) retrained
six predictors in a shared ion space under leakage-controlled splits. Median
spectral angle on MassIVE-KB: Prosit Transformer 0.902, Prosit 0.901,
AlphaPeptDeep 0.871, PredFull 0.833, UniSpec 0.663. On PROSPECT: Prosit 0.862,
AlphaPeptDeep 0.722.

The caveat is decisive and the verifiers confirmed it: **all six were retrained
by the authors on standardised mini-corpora, so this measures architectures, not
the pretrained weights anyone would download.** The paper's own conclusion is
"no single universally optimal model", and it is an unreviewed v1 preprint.
Countervailing: AlphaPeptDeep predicts ~40× faster than Prosit-Transformer
(35 s vs 24 min on 1.4M peptides) — which matters directly at 5M precursors.

Also relevant: the main systematic MS2 benchmark (Ranjan/Lazar, JPR 2024) is
**Orbitrap-only** — zero hits for timsTOF, Astral, ion mobility, Bruker or Sciex
in a full-text scan — and its authors warn it may be biased toward Q Exactive.

## What the research did NOT establish, and it is what decides everything

**LOW.** The two constraints that actually govern an OpenMS integration — **a
working ONNX/C++ inference path**, and **any existing OpenMS-side integration to
build on** — remain **unestablished**. No surviving claim addresses either, and
**every ONNX-export claim put forward was voted down 0-3**, in both directions.

So "Chronologer ships Python-only with no ONNX export" and "AlphaPeptDeep has no
documented ONNX path" are *not* findings — they were killed as unsupported. Which
is consistent with what we know locally: OpenMS already ships AlphaPeptDeep as
ONNX, so an export path demonstrably exists whatever the upstream docs say.

## What I would actually do

1. **Test the local-Triton lead.** Koina's configs are public; if the
   `Prosit_2023_intensity_timsTOF` weights can be pulled and converted to ONNX,
   that is the single highest-value change available — a model trained on the
   instrument whose bias we measured.
2. **Do not chase Prosit on the strength of PepSpecBench.** It compares
   architectures, not shipped weights, and AlphaPeptDeep's 40× speed advantage is
   real at our scale.
3. **RT: take the monotone recalibration first** (doc/32: worth ~35% of our RT
   error, no new model needed), then evaluate DeepLC for per-run fine-tuning
   before Chronologer, whose adaptation mechanism needs a 2.2M-peptide replay.
4. **CCS: leave it.** No evidence base, and our derived 1/K0 already sits at 2.8%
   relative error against measured values.
5. **Astral: accept there is no matched model.** Use QE or Lumos conditioning and
   say so, rather than inventing an index — the mistake removed today.
