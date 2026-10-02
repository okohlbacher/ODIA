# Plan: a trace-level neural classifier for real vs absent precursors

2026-08-24. Written BEFORE any corpus is built, for review by vibe, kimi and
codex.

## Why this shape, and why now

Six measured findings define the problem, and they all point the same way.

1. **Extraction is matched to DIA-NN** (doc/44): 0.270 against 0.311 mean
   pairwise fragment coherence at production settings.
2. **Within-precursor ranking is healthy** (doc/51): 95.6% -- when the right
   peak is in the candidate list, the classifier puts it first.
3. **Admission is not the lever.** Relaxing the gate is -20.4% at matched FDP on
   the full run (doc/54); an ORACLE that bypasses it with DIA-NN's own answers,
   targets only, gains at most +20% (doc/51).
4. **The FDR bucket is not a threshold problem**: recovering it costs 51.3% FDP.
5. **What is left is discrimination at low abundance.** Sub-score AUC at the
   faintest abundance quintile against the brightest: `var_library_corr` 0.532
   -> 0.913, `var_xcorr_shape` 0.598 -> 0.962, `var_ms1_coelution` 0.618 ->
   0.961, `var_xcorr_coelution` 0.696 -> 0.961. Only `var_corr_sum` holds
   (0.897 -> 0.979).
6. **DIA-NN does not summarise where we do.** From the vault's feature matrix
   (`30-Algorithms/OpenSWATH vs DIA-NN feature matrix.md`), DIA-NN 1.7.x's 110
   network inputs are mostly PER-FRAGMENT ARRAYS -- `pCorr[12]`, `pRef[12]`,
   `pAcc[6]`, `pSig[6]`, `pShape[5]` (the raw profile, unsummarised), plus
   `pShadow`/`pShadowCorr[6]` (interference), `pBestCorrDelta`/`pTotCorrSum`
   (competition) and `pAAs[20]`/`pLength`/`pCharge`/`pMz`/`pFrNum` (sequence).
   ODIA collapses the same evidence into 19 scalars and has NO interference,
   competition or sequence family at all.

Collapsing twelve fragments into one number costs most when every fragment is
weak. That is exactly the regime where our features were measured to die. So the
hypothesis is not "a bigger model scores better" -- it is **the summarisation
step is where the low-abundance information is lost, and a model consuming raw
per-fragment traces need not lose it.**

## The label design, which is the crux

Everything else here is engineering. This is where the experiment can be
worthless without being obviously wrong.

**Rejected: target vs decoy as the primary contrast.** A decoy is a permuted
SEQUENCE, so its fragments sit at m/z no real peptide produces. doc/50 measured
the consequence: decoys pick up random coincidence while absent targets pick up
structured signal from real co-eluting species, which is why decoy-implied and
entrapment-implied false positives disagree ~100x in the tail. Training on
target-vs-decoy therefore teaches the wrong contrast, and invites the model to
learn sequence artefacts instead of trace quality.

**Chosen: positives = DIA-NN-confident (q <= 0.01); negatives = ENTRAPMENT
targets.** Entrapment precursors are Arabidopsis sequences in the searched
library (`arab_acc.txt`): real peptides, real fragment masses, and absent from a
human sample by construction. That is the honest "fake" -- nothing about
the trace of an absent real sequence is synthetic.

Decoys are held out entirely, as a calibration check: a model that learned
presence should also rank decoys low, and if it does not, it learned something
else.

### The biggest risk, stated plainly

**Entrapment negatives are a different ORGANISM.** Arabidopsis and human
peptides differ systematically in amino-acid composition, length distribution
and therefore precursor m/z. A model given sequence or composition inputs will
separate the classes perfectly by learning taxonomy, and report a spectacular
AUC that means nothing.

Three mitigations, all mandatory:

* **No sequence or composition inputs in the primary model.** This deliberately
  gives up DIA-NN's `pAAs[20]` family. It is the price of a negative class that
  is not synthetic.
* **Composition-matched sampling.** Negatives are drawn to match the positive
  distribution on precursor m/z, charge, peptide length and fragment count, in
  strata. If the classes cannot be separated on those, the model cannot cheat
  through them.
* **A metadata-only control model**, trained on the scalars alone with all
  intensities zeroed. Its AUC is the LEAKAGE FLOOR: any trace model must beat
  it by a margin, and if the control already scores well the corpus is broken
  rather than the model good.

## The comparative experiment

The point the project has never been able to settle: is our shortfall in the
traces or in what we do with them? doc/44 says the traces are matched on a
coherence statistic; this tests it end to end.

Train the SAME architecture, same split, same labels, on:

* **(a) ODIA-extracted traces** -- `-out_chrom`, which writes long-format
  per-fragment traces from the EXTRACTION path, so it emits them regardless of
  whether Gate C would admit the precursor. Essential: most negatives, and 41%
  of DIA-NN's confident set, are gate-rejected.
* **(b) DIA-NN's own XICs** for the same precursors -- `dn_xic.parquet`,
  45,163,185 rows of `(pr, feature, info, rt, value)`.

If (b) >> (a), extraction is the problem and doc/44's coherence match was
measuring the wrong thing. If they are equal, the traces are fine and the
scoring is where the gap lives. Either answer is worth the corpus on its own.

## Corpus

* **Train: IH1 diaPASEF** (`IH1_diaPASEF.mzpeak`), the run everything else on
  this project is measured on.
* **Validate: IH2** -- an in-house run (IH2), 9.2 GB, a DIFFERENT
  sample from the same acquisition batch. `mzpeak-convert` reads Bruker `.d`
  directly and DIA-NN reads `.d` natively, so this avoids the mzML conversion
  that presents no MS2 spectra to DIA-NN. Genuinely held out: different sample,
  different run, its own true-positive set.
* Sizes: all DIA-NN-confident precursors present in our library (~35k on IH1),
  a composition-matched entrapment sample of the same size, and an equal decoy
  sample for the held-out check. ~105k precursors x ~12 fragments x a 64-cycle
  window is ~80M points -- about 1 GB in float32.
* Extract a SUPERSET deliberately, so that a label redesign after review does
  not require re-extraction.

## Architecture

Fragments are a SET, not a sequence: their order in the library is arbitrary, so
a convolution over fragment index would learn a spurious ordering. The
architecture must be permutation-equivariant.

* Per-fragment token: a small 1-D convolutional encoder over that fragment's
  trace (fixed 64-cycle window centred on the candidate apex, intensity
  normalised per group so the model sees SHAPE not absolute abundance -- with
  the group's total intensity supplied separately as a scalar, since abundance
  is a real covariate and should be available but not smuggled in through every
  channel), concatenated with that fragment's scalars: product m/z, library
  relative intensity, measured ppm deviation, fragment charge, series, ordinal.
* Precursor-level tokens: the MS1 trace through the same encoder, plus
  precursor m/z, charge, and the group's width and total intensity.
* Two to four transformer blocks with self-attention across fragment tokens,
  then attention pooling, then an MLP head.
* This is the same evidence DIA-NN keeps per fragment, without the fixed-12
  truncation its arrays impose.

**Baseline to beat: the shipped 19 scalars through the current GBT, on the same
split.** A neural result that does not beat that is not a result.

## Acceptance tests, fixed in advance

1. **AUC by abundance quintile**, as in `d9_auc_by_abundance.py`. The deficit is
   at Q1; a model strong at Q5 and weak at Q1 adds nothing, because every
   existing shape feature is already strong at Q5.
2. **Held-out IH2**, never touched during training or model selection.
3. **Entrapment-calibrated FDP, never nominal q** -- but note the negatives ARE
   the entrapment set, so the evaluation entrapment population must be a
   DISJOINT held-out split of it, or the FDP is circular. This is a trap the
   corpus construction has to handle, not the evaluation.
4. **Leakage audits**: label shuffle must give AUC 0.5; the metadata-only
   control defines the floor; and the model must rank held-out DECOYS low
   despite never having seen one.
5. **Integration test, only if 1-4 pass**: the model's score as an additional
   sub-score, measured at matched FDP on the full run.

## What would make this fail, and be worth knowing anyway

* The metadata-only control scores as well as the trace model -> the corpus is
  separable without traces and the whole design is void.
* (a) and (b) are equal and both weak -> the information is not in the traces at
  all, and the ceiling is upstream of scoring.
* Q1 does not improve -> summarisation was not the loss, and the low-abundance
  deficit is in the DATA, not its representation. That would redirect the
  project harder than a success would.
