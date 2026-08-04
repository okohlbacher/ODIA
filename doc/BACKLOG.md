# Backlog

Open items, newest first within each group. Items that need a decision or an
external fix are marked **[you]**; the rest are mine to work through.

---

## Blocking someone else

- **[you] `kimi` and `codex` CLIs are unavailable, so the adversarial reviews
  are not the ones you asked for.** Neither binary is installed, and there is no
  `node`/`npm`/`npx` to install them, nor any `MOONSHOT`/`KIMI`/`OPENAI`
  credential in the environment — so this is a credentials problem, not just a
  missing binary. Overnight reviews are being done with independent subagents
  prompted to refute instead, which is a weaker substitute: same model family,
  so correlated blind spots. Install the CLIs and provide keys and I will switch.

- **[decided 2026-08-03, user] Phase 2 proceeds on the slow mzPeak reader.**
  The decode fix is being worked on elsewhere; ODIA does not wait for it. The
  extractor is written against `SpectrumSource` (D10), so when the batched
  decoder lands it is a substitution, not a rewrite.
  **Do not let this be forgotten:** every extraction timing measured before the
  replacement lands is a measurement of the reader, not of ODIA, and must be
  re-taken afterwards. `doc/05-mzpeak-batched-reader-handoff.md` is the spec.

- **[for the reader author, not us] mzPeak decode diagnosis.** At `f93f938` the batch API exists but does not amortise:
  284.4 ms/spectrum via `get_spectra_batch`, against 0.082 ms/spectrum if each
  Parquet row group were decoded once and its 873 spectra served from it.
  Measured, not estimated -- a row group decodes in 0.07 s with pyarrow, and a
  full pass over `12_80` would take 1.5 s, beating mzML's 3.5 s. The reader
  decodes a whole row group per spectrum and keeps ~0.1% of it. Everything that
  touches peaks inherits this, including `extract_ion_chromatogram`.

- **[superseded, kept for the measurement] mzPeak peak decoding is ~1000× too slow.** ~277 ms per
  spectrum against ~0.27 ms for OpenMS parsing mzML. Full requirements in
  `03-mzpeak-streaming-requirements.md` (R1, R2, R5 are the blocking set).
  Phase 2 cannot start until this moves.

- **[decided 2026-08-03, user] `.oswpq` writing is not needed.** Dropped from
  scope. Reading stays, and accepts either intensity width. The note below is
  kept only because the reader's behaviour depends on it.

- **[was: you] `library_intensity` `float64` vs `float32` in `.oswpq`.** The
  predecessor's `float32` change is a breaking on-disk change that was never
  upstreamed, and we cannot apply it without modifying OpenMS. Proposal on the
  table: read both, write upstream-compatible `float64` by default with
  `float32` behind an option. **The reader now does accept both**, verified
  against fixtures of each width, so only the writer is still blocked on this.

- **[resolved 2026-08-04] GPU access.** These nodes authenticate with
  **Kerberos**, not keys -- `ssh spock` works from an interactive login because
  PAM puts a TGT in the session keyring, which an agent session cannot reach.
  A file-based ticket (`KRB5CCNAME=FILE:<path> kinit`) hands one over without
  creating a standing credential, and it expires on its own. The CUDA path is
  now measured; see "GPU inference measured on an H100" below.

- **[you] GPU 0 is unusable on spock and data.** `cudaSetDevice(0)` returns
  error 46 with both devices idle, in Default compute mode, and `/dev/nvidia*`
  world-writable. On `data` device 1 works, so `CUDA_VISIBLE_DEVICES=1` is a
  workaround; on `spock` **both** devices fail, leaving it with no usable GPU.
  This is an infrastructure fault, not an ODIA one -- worth a ticket.

---

## Correctness, must be settled before results are trusted

- **Validate the PeptDeep modified-peptide path against AlphaPeptDeep reference
  values.** The encoding is reconstructed from upstream sources and documented in
  `04-peptdeep-encoding.md`, and the unmodified path can be validated against
  OpenMS exactly — but nothing here can confirm the modification path. A wrong
  `mod_x` yields plausible, quietly wrong intensities. **No library generated
  with modifications should be trusted until this is done.**

- **Isotope element mapping.** AlphaPeptDeep's element list contains `2H`,
  `13C`, `15N`, `18O`; OpenMS writes isotopes as `(13)C`. Without an explicit
  translation, isotope-labelled modifications land silently in the `?` bucket.
  Adversarial review confirmed these four renames are the *only* mapping needed:
  all 40 element symbols in the shipped `unimod.xml` are in the 109-element list.

- **Mod feature counts are signed.** 600 of 2859 rows in the UniMod-derived table
  have negative counts (`Deamidated@N` is `H(-1)N(-1)O(1)`), so an unsigned type
  or an `abs()` on the OpenMS side is silently wrong. Measured: RT 0.5231 signed
  against 0.5121 with `abs()`.

- **Batches must be length-homogeneous.** Trailing padding is not inert -- index 0
  is one-hot encoded and no model applies a padding mask. The same peptide padded
  to `seq_len` 13/20/30 predicts RT 0.8196/0.4737/0.2714. Group by encoded length
  and run one batch per length.

- **Determinism harness.** Permutation invariance and thread invariance at
  1/8/64 threads with OpenMP linked, per the cross-cutting invariants. Not yet
  built; it needs to exist before any A/B measurement is believed.

---

## From the Phase 1 adversarial review — still open

Four criticals from that review are fixed (decoy fragment masses, the
target/decoy round-trip merge, non-contiguous row detection, Parquet numeric and
string type coverage), plus the `toFixed` domain, the unchecked output stream and
the exit-code-6-on-success. These remain:

- **N-terminal modifications break the decoy tokeniser.** `(UniMod:1)PEPTIDE`
  has its modification name tokenised as residues. **Reachable from the supplied
  DIA-NN fixture**, which holds 441 such precursors — 428 targets (1.3%) get no
  decoy as a result. The count is now reported rather than silently dropped, but
  the tokeniser still needs fixing.
- **The Parquet reader still cannot open a library past 2 GB of characters.** The
  `large_string` handling casts down to `arrow::utf8()`, which has the 32-bit
  offsets `large_string` exists to avoid, and the plain-`string` path fails in
  `CombineChunks`. Both emit a message blaming the file for a missing column that
  is present. `LargeStringArray` needs to be a first-class case.
- **`sortByPrecursorMz` doubles peak memory** (+483 MiB, +39% on the human
  proteome) by materialising complete copies before moving. An in-place
  permutation or a block-streamed rebuild avoids it.
- **The contiguity check is order-dependent.** Two blocks of one precursor that
  are *adjacent* still merge silently, discarding the second block's RT, IM,
  precursor m/z and protein group; only non-adjacent repeats are detected.
- **Null cells are indistinguishable from zero** for `Decoy`, `Precursor.Charge`,
  `Fragment.Series.Number` and the string columns — a null `Decoy` column makes
  every decoy a target, the same effect the BOOL case was added to fix.
- **Unusual loss labels are not round-tripped.** `H2O+H2O`, `CH3SOH`,
  `H3PO4+H2O` all parse to `LossType::Other` and are written back as `other`, so
  the label is lost and the transition is dropped from decoys. Either carry the
  label as an interned string or refuse the input.
- **A decoy's `Product.Mz` can no longer be checked against its own row.** Since
  the decoy stores the target's sequence (DIA-NN convention), reproducing its
  fragment m/z requires applying the mutation table — which has become an
  unwritten part of the file format. Defensible, but it should be stated in D7.
- **The `_decoy` id suffix deviates from DIA-NN**, whose decoy `Precursor.Id` is
  exactly `sequence+charge`, identical in form to a target's. Ours has to differ
  because we emit both rows; worth documenting as an intentional divergence.
- **Contiguity-check hash collisions are a hard failure.** 4.3e-7 at 4 M
  precursors (fine), but 1.7e-4 at the 78.6 M design scale — one load in 6,000 —
  and the failure is a deterministic `throw` whose message sends the user to
  re-sort an already-sorted file.
- **`generate()` interns before deciding whether a precursor survives**, leaving
  5,532 permanently unreferenced arena entries on the human proteome and
  over-reporting the distinct-string count by that much.
- **D5's mandated warning is not implemented** — the reader should warn once when
  an input's transition names disagree with what ODIA would synthesise. The
  accepted round-trip risk is currently unmitigated.
- **`Decoy` given as `True`/`False` parses to 0**, and `toLong` cannot
  distinguish absent from unparseable.
- Smaller: `char buf[64]` truncates long numeric fields and `strtod` accepts
  partial parses (`500,1` becomes 500); no `static_assert` that the two mutation
  tables are the same length; `loadTSV` never reserves; a handful of
  `pseudo_reverse` decoys are identical to their target for palindromic
  prefixes.

---

## Library generation is ~5x slower than it needs to be (measured 2026-08-04)

Library generation is 53 min for the human proteome and dominates end-to-end
cost by ~40x. Two causes, both in how ONNX Runtime is driven, neither in the
models themselves:

**1. `-threads` never reaches inference.** `LibraryGenerator` constructs all
four predictors as `PeptDeepPredictor(model_path, prefer_gpu)` — two arguments,
so `intra_op_threads` takes its default of 0 and `SetIntraOpNumThreads` is
never called. ONNX Runtime then sizes its pool to the whole machine (128 cores
here) regardless of what the user asked for. The tool's `-threads` flag governs
digestion and decoy construction and nothing else, which is ~7% of the phase.
Sites: `LibraryGenerator.cpp:343` (RT), `:377` (MS2), `:561` (iRT), `:631` (CCS).

**2. Intra-op parallelism is the wrong axis for these models.** PeptDeep's MS2
network is recurrent; its per-op tensors are too small to spread across many
threads. Measured on 16 pinned cores, one session, MS2 throughput against
intra-op thread count:

| intra-op | peptides/s |     | processes x threads | peptides/s |
|---------:|-----------:|-----|--------------------:|-----------:|
|        1 |      168.7 |     |          16 x 1     |   **1760** |
|        4 |      338.7 |     |           8 x 2     |     1253   |
|        8 |  375.6 (peak) |  |           4 x 4     |      815   |
|       16 |      317.9 |     |           2 x 8     |      505   |
|       32 |      231.2 |     |                     |            |

Intra-op saturates at 8 threads for 2.2x and then *degrades*. Data parallelism
— N independent sessions, one thread each — is monotonically better and beats
the best intra-op configuration by **4.7x on identical cores**. It keeps
scaling: 8 cores 1063/s, 16 cores 1480/s, 32 cores 2716/s, 64 cores 4585/s
(per-core falls 133 -> 72 as memory bandwidth saturates, but aggregate climbs).

Against the 897 peptides/s the v5 run actually achieved, 64 data-parallel
sessions measure **5.1x**. Projected: MS2 2372 s -> ~465 s, and library
generation 53 min -> **~14 min** if RT and CCS gain proportionally (they share
the same defect; not separately measured).

**The fix**: give each worker thread its own `Ort::Session` with
`intra_op_threads = 1` and slice the block across them, replacing the single
shared session. Cost is memory — 129 MB RSS per standalone process, though
in-process the marginal cost is the arena plus a 16 MB weight copy, so a 64-way
split should land near 3-5 GiB against today's 2.1 GiB peak. Worth capping the
session count independently of `-threads` for that reason.

Both measurements were taken while another user held ~114 of the 128 cores, so
the ratios are back-to-back under equal contention and the absolute rates are
floors.

## Benchmark ODIA against other DIA tools via ProteoBench [user, 2026-08-04]

Explore <https://proteobench.cubimed.rub.de/>, download the DIA benchmark, and
design a workflow that compares ODIA to the other DIA tools on it.

Why this is worth doing properly rather than quickly: every comparison made so
far has been against DIA-NN on S08, with DIA-NN's own library as the reference
and DIA-NN as the search engine. That measures agreement with one tool on one
file, and it cannot distinguish "ODIA is good" from "ODIA resembles DIA-NN".
ProteoBench supplies a defined ground truth and a published leaderboard, which
is the first thing here that could falsify a claim rather than confirm it.

Points to settle when this is picked up:

* Which module — the DIA modules differ in organism mix and in whether the
  quantitative ratios or the identification counts are the scored quantity.
* ODIA does not yet do quantification, so the honest first submission may be
  identifications only. Check what the leaderboard requires before building to it.
* The entrapment measurement already on this list is a prerequisite, not a
  parallel task: a leaderboard position computed on an uncalibrated FDR is
  worse than no position.
* Submission means publishing a result under our name. Do not submit anything
  without asking first; running the benchmark locally needs no such permission.

## [you] The pyProphet/XGBoost scoring code is not reachable (2026-08-04)

The decision is to use the pyProphet-like XGBoost scoring already written at
`/Users/kohlbach/Claude/mzPeak/OpenDIAlyzer`. That is a **macOS path**; the
nodes this project builds on are Linux and have no `/Users` at all, so the code
cannot be read, let alone called. Copy it to `/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/`
(or anywhere on Ceph) and Phase 3's classifier is unblocked.

Upstream pyProphet 3.0.15 with XGBoost 3.2.0 *is* installed at
`/ceph/ibmi/abi/oliver/envs/pyprophet` -- a useful cross-check, but not the
code the decision named, so it is not a substitute without saying so.

Peak-group detection and the sub-scores do not depend on this and proceed
meanwhile; see `doc/07-scoring-plan.md`.

## Implementation, unblocked

*Still uncovered by any fixture, from the mutation-testing review: neutral loss
on 3+ charged fragments; a phosphorylated modification (the checker's MOD table
has no UniMod 21, so it cannot yet be pointed at a phospho library, which is the
library type `lossMass`'s own comment cites); `pseudo_reverse` decoys, which no
test invokes and which `check_invariants.py` cannot validate since its `mutate()`
implements only the substitution table; and a precursor split across two
**adjacent** blocks.*

*Fixture coverage is otherwise addressed: `scripts/make_adversarial_fixtures.py` plus
`test/check_invariants.py` give 20 tests over losses, modifications, decoy
collisions, unusual Parquet types, nulls and malformed input, with masses
verified against a residue table independent of OpenMS and pinned by a
self-test of literal values. Mutation testing: 21 injected defects, all of the
ones re-tested after repair are now caught.*

- **iRT prediction is wired into library generation.** Predicted once per
  distinct modified sequence rather than per precursor, since the RT model has
  no charge input: 1.25 M predictions instead of 4.0 M on the human proteome,
  206 s, and 232 s for the whole run at 1.9 GB peak. Unpredictable peptides are
  left NaN and counted, never given a made-up value.
- **iRT prediction is done and validated**: `PeptDeepPredictor` runs the RT model
  with CUDA attempted and CPU fallback, and its predictions match the
  independent Python reference to 1e-5 on mixed-length batches including
  residue, N-terminal and C-terminal modifications. Still to write: MS2 and CCS,
  which add `charges` (×0.1), `nce` (×0.01) and a rank-1 `instrument_indices`,
  and the wiring that replaces the placeholder intensities and iRT in
  `LibraryGenerator`.
- Earlier note, retained: the **encoder is validated**: `PeptDeepEncoder` agrees
  exactly with the independent Python reference on residue, N-terminal and
  C-terminal modifications, multiple modifications, signed counts and
  `aa_indices`.
- **Predicted iRT is the model's normalised output, and nothing says so in the
  file.** Measured range over 109,864 human sequences: -0.039 to 0.940, with
  14.6% above 0.85 (saturation) and 1.24% negative. Nothing converts it to iRT
  units or minutes, and the invariant checker treats a literal `0` as missing
  while accepting -0.039. Decide the unit and declare it.
- **Prediction costs 3 h 07 m of CPU time** (11,267 CPU-seconds at 4688% for
  206 s wall) for one proteome library. Fine here; worth knowing before it runs
  anywhere metered, and an argument for the GPU path.
- **ONNX Runtime's arena accounts for +610 MiB of the +688 MiB** prediction adds
  at proteome scale, and it never shrinks. It is bounded by `MAX_BATCH_ROWS` x
  encoded length, so `-max_peptide_length` scales it linearly with nothing
  capping or testing that.
- **Peptides longer than the model was validated on are unreachable only by
  accident**: `precursor_mz_max = 1200` at charges 2-3 caps peptides at ~32
  residues before a precursor exists. `PeptDeepEncoder` has no length guard, so
  widening the m/z window would feed the model lengths it has never seen.
- **The library records neither the model nor the execution provider used.**
  CUDA and CPU differ in the last bits, so a written library is currently
  unattributable.
- **Predictions depend on batch composition at the last bit.** Within one build,
  repeated calls are bit-identical, but batch size moves the result by ~1 ULP
  (1.5e-8) and the two ONNX Runtime builds compared in the tests differ by up to
  2.7e-7. `MAX_BATCH_ROWS` pins the batch size, which bounds it, but D1 asks for
  determinism as a precondition and this is not zero. Decide whether that
  matters at the iRT-window scale.
- **The CUDA path is still untested** — `spock`/`data` remain unreachable, so
  only the CPU branch has ever run. `Ort::GetAvailableProviders()` is queried
  before attempting, so a CPU-only build no longer logs a scary provider-load
  error, but that also means the CUDA branch is skipped rather than exercised
  here.
- **A doubly-modified residue is under-encoded.** OpenMS's `AASequence` keeps
  only the last modification on a residue, so it is lost before the encoder sees
  it; AlphaPeptDeep would accumulate. Decide whether to detect and refuse, or to
  carry modifications outside `AASequence`. The reference encoder and RT
  predictor (`test/peptdeep_reference.py`) are in place as the validation oracle,
  written from the spec rather than from the C++; the element list is generated
  from the authoritative yaml into `data/peptdeep_mod_elements.txt` and
  `include/odia/PeptDeepElements.h`. The C++ side is still to be written.
- **MS2 prediction works and matches the reference**; still to do is *using* it:
  `LibraryGenerator` caps fragments by descending m/z because it has no
  intensities yet. Wiring it in means mapping the `[positions, 8]` output onto
  b/y ions at charge 1-2, choosing a cap by predicted intensity, and deciding
  what NCE and instrument to assume when the caller does not say.
- **NCE and instrument are guesses.** They default to 30.0 and "QE", and nothing
  derives them from the data. They materially change the spectrum, so a library
  generated for one instrument is not right for another.
- CCS prediction, for the ion-mobility column.
- `.oswpq` **write** (still gated on the `float64`/`float32` decision above).
- Phase instrumentation: wall, CPU, RSS, `mallinfo2`, and node load per run,
  charging un-phased time to the preceding phase.
- mzPeak chromatogram writer, using the fork's writer.

---

## Decisions taken 2026-08-03 (user)

- **Development proceeds against a DIA-NN library**, not ours, so that Phase 2
  work is not confounded by library differences. Ours stays the deliverable;
  it is simply not the variable under test while the extractor is built.

- **Decoy construction is deferred, not solved.** Ours cost 14,082 precursors
  (45%) when DIA-NN searched with them instead of building its own. They keep
  the target's m/z, iRT and intensity pattern with only fragment masses moved,
  which is DIA-NN's published design, so the fault is more likely in how the
  intensity pattern is copied than in the mutation scheme. Until this is
  understood a library written for another engine should probably carry no
  decoys at all. **This is an FDR question, and it blocks Phase 3, not Phase 2.**

- **NCE and instrument stay guesses for now** (30.0 and "QE"). Nothing derives
  them from the data and they materially change every predicted spectrum, so a
  library generated for one instrument is not right for another. Revisit when
  the extractor can measure the mismatch rather than assume it.

---

## Phase 1 measured against DIA-NN on S08 (2026-08-03)

Both libraries predicted from the same FASTA with the same digest settings;
DIA-NN searched the same diaPASEF run with each, everything identical but the
library.

| library | precursors | proteins |
|---|---:|---:|
| DIA-NN's own | 37,170 | 5,255 |
| ODIA, first attempt | 17,478 | 3,263 |
| + targets only (DIA-NN makes its own decoys) | 31,560 | 4,776 |
| + ion mobility from predicted CCS | 33,732 | 5,021 |
| + fragment charge up to precursor charge | 34,609 | 5,053 |

**Yield per library precursor is identical**: DIA-NN 1.732%, ODIA 1.736%. The
assays are as good, one for one; every remaining difference was coverage.

Still open from this:

- **Our decoys cost 14,082 precursors, 45% of the total.** Shipping them made
  DIA-NN use ours instead of building its own, and ours are too target-like to
  separate: same m/z, same iRT, same intensity pattern, only fragment masses
  moved. That is DIA-NN's own published design, so the fault is more likely in
  how the pattern is copied than in the mutation scheme. **Until this is
  understood, a library written for another engine should probably not carry
  decoys at all** -- which is a decision, not a default to pick silently.

- **Retention time is the weakest remaining signal.** DIA-NN set an RT window of
  2.29 min for our library against 1.44 for its own, so our iRT is about half as
  predictive after calibration despite a rank correlation of 0.992 with theirs.
  Worth attacking next; it is a scoring dimension, not a coverage problem.

- **Fragment selection overlaps only 0.60** with DIA-NN's, and base peaks agree
  62% of the time, while intensities on the fragments both keep correlate at
  r = 0.861. So the models mostly agree and the caps diverge. Not obviously
  wrong, but not understood either.

- **CCS agrees with DIA-NN's ion mobility to 1.43% median** over 1,873,932
  precursors, with no charge or length dependence and a split-half held-out
  median of -0.004%. Both predictors are sound; see `test/compare_ccs_to_diann.py`.

---

## RT and m/z calibration are independent -- measured 2026-08-03

Asked whether the m/z tolerance has to be recalibrated whenever the retention
times are. **No.** Across three libraries differing only in their RT column,
DIA-NN's achieved mass accuracy is identical to six significant figures while
RT accuracy varies 2.2-fold:

| library | RTPredAcc | median MS1 | MS1 corrected | median MS2 | MS2 corrected |
|---|---:|---:|---:|---:|---:|
| v3 stock RT | 0.361 | 4.85479 | 1.16167 | 1.48522 | 1.06719 |
| v5 fine-tuned | 0.211 | 4.85479 | 1.16167 | 1.48522 | 1.06719 |
| v4 fine-tuned | 0.164 | 4.85479 | 1.16167 | 1.48522 | 1.06719 |

Mass error is a property of the instrument and of the library's m/z values;
recalibrating retention times changes neither. The only coupling observed is
second-order: DIA-NN's *chosen search tolerance* moved between 7 and 10 ppm
across the variants and not monotonically in RT accuracy, which is its coarse
parameter grid reacting to a different identification set rather than a change
in the underlying mass error.

**The finding that matters for ODIA is that it calibrates neither.**

- `ChromatogramExtractor::Options::fragment_ppm` is a fixed **20.0**. The
  measured median MS2 mass error on this run is **1.49 ppm**, and DIA-NN
  optimised its own search tolerance to 7-10 ppm. Ours is roughly 2-3x wider
  than the data supports, and every extra ppm admits proportionally more
  interfering peaks into every transition's chromatogram. It is a guess, and
  nothing currently measures it.
- `rt_low`/`rt_high` are caller-supplied, global, and the extractor never reads
  a retention time from the library at all (see section 4a of
  `doc/06-rt-refinement-plan.md`).

Both should be derived from the run. Because they are independent, they can be
derived independently -- a joint or iterative calibration is not needed, which
is worth knowing before someone builds one.

---

## Phase 2 first slice, measured 2026-08-03

Extraction runs end to end: library -> window assignment -> one forward pass ->
chromatograms. On `12_80.mzpeak`, 200 precursors / 2,400 transitions over a 60 s
retention-time slice (797 spectra):

| | |
|---|---|
| points | 176,268, 14.9% non-zero |
| **decode** | **262.03 s** (328.8 ms/spectrum) |
| **match** | **0.02 s** (0.021 ms/spectrum) |
| memory | 1.4 MiB |

**The reader costs 13,000x what the matching costs.** ODIA's own work is 0.006%
of the runtime. With the batched decoder specified in
`05-mzpeak-batched-reader-handoff.md` the same extraction would take about
0.09 s rather than 262 s, and matching would become the dominant term -- which
is the point at which optimising ODIA starts to be worth anything.

Consequences to keep in view:

- **No timing taken through the current reader means anything about ODIA.**
  Recorded here so it is not quoted later as an ODIA benchmark.
- **The retention-time range must be resolved to a spectrum index range before
  decoding**, not filtered afterwards. Filtering afterwards read the whole run
  and discarded most of it -- an hour instead of seconds. Fixed, and worth
  remembering as the shape of mistake this reader punishes.
- **Still missing before this is useful**: an iRT-to-RT calibration, so a
  precursor is extracted over a window around where it should elute rather than
  over the whole run; and an output format. Without the first, a whole library
  over a whole run is billions of points.

---

## RT fine-tuning: the controlled result (2026-08-03)

Five searches of S08, everything identical but the library's RT column.

| library | RT source | RTPredAcc | DIA-NN's window | precursors | proteins |
|---|---|---:|---:|---:|---:|
| DIA-NN's own | DIA-NN | 0.2133 | 1.4397 | 37,170 | 5,255 |
| ODIA v3 | stock PeptDeep | 0.3610 | 2.18905 | 35,296 | 5,109 |
| ODIA + DIA-NN's RT | DIA-NN, pasted | 0.2189 | 1.36581 | 35,523 | 5,130 |
| **ODIA v5** | **fine-tuned on 500** | **0.2107** | -- | **35,131** | **5,091** |
| ODIA v4 | fine-tuned on 23,179 | 0.1642 | 1.05543 | 37,466 | 5,295 |

**v4 is contaminated and must not be quoted.** 62.5% of the peptides it
identified had their observed retention time in its own training set. It is
train-on-test at the run level, and the +296 over DIA-NN is memorisation.

**v5 is the honest one** -- 500 training peptides, 1.5% of the run's
identifiable peptides -- and it settles the question:

**RT accuracy improved from 0.361 to 0.211, reaching DIA-NN's own 0.213, and
identifications went DOWN by 165.**

That reproduces the earlier RT-swap ablation (+227, also noise) with an
independent mechanism, and it means the original conclusion was right: **on this
run, retention-time accuracy does not buy identifications.** Two libraries with
equal RT accuracy -- v5 at 0.2107 and DIA-NN at 0.2133 -- differ by 2,039
precursors, so the remaining gap is somewhere else entirely.

**Where "somewhere else" is** remains the fragment-charge and intensity question
recorded below: we emit 78.7% singly-charged fragments where DIA-NN emits 72.4%,
and the precursors we miss are enriched in charge 3.

**What fine-tuning is still for:** the extraction window, which v4 measured at
1.055 min against v3's 2.189. That is a Phase 2 compute and feasibility argument
(`doc/06-rt-refinement-plan.md` section 4a), not a sensitivity one -- and it now
has an experiment behind it rather than an assertion.

---

## Retention-time fine-tuning is integrated (2026-08-03)

`scripts/finetune_rt.sh <report.parquet> <outdir> [n]` fine-tunes the RT model on
a run's own identifications and exports ONNX that ODIA reads with `-rt_model`.
No C++ changed, and torch stays out of the runtime -- OpenMS's own exporter
already emits the input names ODIA expects.

Two things the integration fixed over the prototype:

- **Modifications resolve through alphabase's own table**, 1,524 UniMod ids
  against the prototype's five, with the residue disambiguating candidates and
  a refusal rather than a guess when it cannot.
- **A provenance sidecar** records the model's SHA-256, the identifications it
  was tuned from, the peptide count, and the minute range defining the scale it
  predicts on -- which nothing downstream can recover otherwise.

Still open, and both are in `doc/06-rt-refinement-plan.md`:

- **Method choices settled by measurement**, not preference:
  - *Direct target beats residual learning ON S08, which is one file.* Both
    mechanisms are implemented and selectable (`--method direct|residual`)
    precisely because a single benchmark should not decide it. Measured here:
    direct 0.429 against residual 0.635 at n=500 / 500 epochs, and 0.415 against
    0.591 at n=600 / 60 epochs. Residual IS a real improvement over the stock
    model (42%), just a smaller one than direct (59%). The likely reason is that
    retargeting pretrained weights at a small centred residual fights the
    initialisation; freezing the trunk was not tried and might reverse it.
    `--evaluate` makes every run report which method won on its own data.
  - *Data beats optimisation steps.* 500 peptides at 500 epochs reaches 0.429;
    10,000 peptides at 40 epochs reaches 0.359. Twelve times the gradient steps
    does not close what more peptides close easily.
  - *40 epochs is about right.* At n=2000, going to 120 buys 4.7% while the
    train/held-out gap doubles from 0.097 to 0.198 -- overfitting onset without
    the held-out curve having turned up yet.

- **The library does not carry the model's identity.** The sidecar sits beside
  the model, not beside the library, so two libraries built with different
  models still look identical. That is the remaining half of provenance.
- **Feeding this from ODIA's own first pass closes a loop.** The script cannot
  tell an external search from its own and records the source rather than
  judging it. The entrapment measurement is the precondition.

---

## Retention time — accurate diagnosis, and it does NOT cost identifications

Researched 2026-08-03. **The framing that sent this investigation was wrong,
and the correction matters more than the investigation.** RT was called "the
most likely source of the remaining 5%". It is not.

**The controlled ablation.** ODIA's library with DIA-NN's RT column pasted in,
everything else identical:

| library | RT column | DIA-NN's window | precursors | proteins |
|---|---|---:|---:|---:|
| DIA-NN's own | DIA-NN | 1.4397 | 37,247 | 5,255 |
| ODIA v3 as shipped | PeptDeep raw | 2.18905 | 35,796 | 5,109 |
| **ODIA + DIA-NN's RT** | DIA-NN | **1.36581** | **35,688** | **5,130** |
| ODIA + linear iRT calibration | 152.2356x - 39.2322 | 2.18905 | 35,831 | 5,088 |

A *perfect* RT column buys **-108 precursors and +21 proteins**. The gap to
DIA-NN survives intact. The window is a faithful measure of RT accuracy -- the
peptide sets are identical, so library size cannot be the confound -- but RT
accuracy is not what is costing identifications on this run.

**Implemented 2026-08-03.** The calibration is applied, refitted at run time
from `data/irt_standards.tsv` rather than pinned, and verified against OpenMS's
own fixture (`peptdeep_irt_peptides_predicted.csv`), which carries both the raw
ONNX prediction and the calibrated iRT for the 11 standards. ODIA reproduces the
published line to six digits: 152.2356 x raw - 39.2322, worst standard off by
8.7570 iRT. The invariant checker no longer demands the library's RT equal the
model's output -- it recovers the affine map from the data and requires it to be
the one the standards give, which catches both a missing calibration and a
dropped intercept.

**The iRT calibration is a no-op for accuracy, bit-exactly.** Applying the AlphaPeptDeep
Biognosys linear calibration produced an *identical* window, 2.18905. This is
structural, not luck: DIA-NN fits its own monotone calibration, so any monotone
reparametrisation of our iRT axis cannot change anything. Worth applying anyway
for **unit hygiene** -- a field named `irt` holding a training-gradient
coordinate is the same silent-units hazard the CCS decision exists to avoid --
but it must never be described as an accuracy fix. If applied, recompute the
constants at runtime from the 11 Biognosys peptides; they are checkpoint-specific
and would rot silently on a model swap.

**The real RT defect is output saturation, and it is unfixable by calibration.**
Local calibration slope by decile of our predicted iRT, normalised to the median:

```
ODIA    1.05 0.97 0.97 0.96 0.90 0.95 0.95 1.22 5.11 11.88
DIA-NN  0.89 1.04 1.03 0.98 1.09 1.01 1.02 0.97 0.97 0.84
```

The top two deciles are compressed 5x and 12x; 12.8% of the library piles into
the single bin 0.85-0.90. Confirmed against the model directly with synthetic
ladders: `G14` 0.016 -> `G9L5` 0.876 -> `L14` 0.936, so five leucines spend 0.89
of the range and the next nine spend 0.06. This is the standard artefact of a
min-max-normalised training target. The information is gone, not misplaced,
which is exactly why the calibration experiment came back null.

Worst classes, all one axis seen four ways (long, hydrophobic, charge 3, late):
length 25-30 is 2.06x DIA-NN's spread, late RT 1.73x, GRAVY top decile 1.67x,
charge 3 1.55x. Nothing systematic in modifications or terminal residue.

**Where the effort should go instead (R1).** With RT equalised we still miss
4,937 of DIA-NN's precursors while finding 3,378 it does not. The missed ones
have **median quantity 0.40x** the shared ones and are **enriched in charge 3
(37.0% vs 24.5%)**. Libraries are structurally matched (12.00 vs 11.92
fragments/precursor, b/y 29/71 both) but ODIA emits **78.7% singly-charged
fragments against DIA-NN's 72.4%**. That is the thread to pull: fragment
selection for multiply-charged precursors, not retention time.

**Where RT still pays (R3).** Extraction cost. A 1.6x narrower window is 1.6x
less XIC per precursor -- ~128 cycles against ~80 on this run. That is a direct
constant factor on Phase 2's runtime and peak memory, which is the one place
this still matters.

**Deferred, with reasons.** Fine-tuning the RT model on a run's own
identifications is the only thing that can touch the saturation (Zeng et al.,
Nat Commun 13:7238, report R2 0.927 -> 0.986 from 500 peptides). But R1 says the
identification benefit is ~0, and the self-contained path costs a torch
dependency ODIA does not have. DIA-NN can do it in minutes with bundled libtorch
(`--tune-rt`), which would make ODIA depend on DIA-NN at library-build time --
probably unacceptable for a tool positioning itself as an alternative.
A run-learned residual correction is worth -21% sd, cross-validated by peptide,
and becomes useful only once ODIA scores its own data.

**Caveats that must travel with these numbers.** Single run, single gradient,
single library size -- every result is n=1 at the experiment level, and the
RT-swap null deserves a second run before being treated as settled. No verified
cross-tool RT accuracy comparison exists in a common unit on a common dataset;
do not let one into this document as fact. The tensor-level encoding is verified
correct: ODIA's outputs are bit-identical to OpenMS's reference values and to
AlphaPeptDeep's torch predictions.

---

## Fragment charge ranking — measured, and it is a MODEL difference (2026-08-04)

The suspicion was that raising the enumeration cap put doubly-charged fragments
into the candidate list without the ranking ever keeping them. Measured, per
precursor charge, on ODIA v3 against DIA-NN's library:

| precursor z | frag z | ODIA share | DIA-NN share | ODIA median intensity | DIA-NN median |
|---|---|---:|---:|---:|---:|
| 2 | 1 | 82.5% | 78.3% | 0.2547 | 0.2330 |
| 2 | **2** | **17.5%** | **21.7%** | **0.0983** | **0.1546** |
| 3 | 1 | 74.2% | 64.3% | 0.2617 | 0.2288 |
| 3 | **2** | **25.8%** | **35.7%** | 0.2386 | 0.2527 |

**Our selection is faithful to our own model.** For a 2+ precursor the PeptDeep
model puts a doubly-charged fragment at median 0.098 where DIA-NN's model puts
it at 0.155, so ranking by intensity demotes it exactly as it should given what
we believe. Both tools cap fragment charge at 2, so the candidate pools match
and the rule is not the difference -- the two intensity models are.

That reframes the fix. Forcing more z2 into the library **overrides the model
rather than correcting a defect**, which may still be worth doing if it finds
more peptides, and is worth nothing if it does not. So it is an option
(`-reserved_doubly_charged`, default 0) and an ablation rather than a change:
same build, same FASTA, same RT model, same search, only the quota differing.
Measured on a 200-protein library, a quota of 4 moves the z2 share from 21.6%
to 36.2%, which brackets DIA-NN's mix.

**Result: the quota costs identifications, so the default stays 0.** Full human
proteome, same build, same FASTA, same RT model, same search, only the quota
differing:

| reserved_doubly_charged | precursors | proteins |
|---|---:|---:|
| **0** (default) | **35,348** | **5,141** |
| 4 | 35,050 | 5,051 |

Forcing DIA-NN's charge mix loses 298 precursors (-0.84%) and 90 proteins
(-1.75%). So the model difference is real but overriding it is not an
improvement: our ranking is not merely faithful to our model, it is *better for
the search* than the mix it was being compared against. The option stays --
it is how this was measured and how it would be re-measured on other data --
but nothing should set it without repeating the ablation.

This closes the fragment-charge question.

---

## Fragment charge rules — the original suspicion, kept for the record

Raised by the RT research (2026-08-03), which ruled retention time out and left
this as the leading explanation for the remaining identification gap.

**The measurements.** With RT equalised between the two libraries, ODIA still
misses 4,937 of DIA-NN's precursors while finding 3,378 it does not. The missed
ones are **enriched in charge 3 -- 37.0% against 24.5% of the shared set** --
and have **median quantity 0.40x** the shared ones. The libraries are otherwise
structurally matched: 12.00 fragments per precursor against 11.92, b/y split
29/71 in both. But **ODIA emits 78.7% singly-charged fragments where DIA-NN
emits 72.4%**.

**Why this is not already fixed.** The enumeration rule was corrected earlier:
fragments may now carry up to the precursor's own charge, which was worth 44.7%
of the assays DIA-NN had and we lacked. What was never checked is whether the
*ranking* then selects them. Raising the cap only puts doubly-charged ions into
the candidate list; the top-12 cut is by predicted intensity, and if the model's
z2 channels are systematically weaker than its z1 channels, the cap change adds
candidates that never survive. The 78.7%/72.4% split says exactly that is
happening.

**Questions to answer, in order.**

1. Is the z2 deficit in the *model* or in our *use* of it? Compare the predicted
   intensity distribution of channels 1 and 3 (b_z2, y_z2) against channels 0
   and 2, for 3+ precursors specifically. The MS2 comparison already
   demonstrated that channel 1 is barely reached at charge 2 and needs a 4+
   precursor with two basic residues to appear at all.
2. Does DIA-NN cap fragments the same way? 11.92 per precursor is close enough
   to our 12 to suggest a similar cap, but the composition differs, so either
   the ranking or the candidate set does.
3. Would a charge-aware cap help -- for example reserving slots for z2 ions on
   3+ precursors rather than letting one global ranking decide?
4. Do the missed precursors' assays actually extract? Their 0.40x median
   quantity says they are low-abundance, so some of the gap may be sensitivity
   rather than assay choice. Phase 2 can answer this directly once it can
   extract at scale, and that is the cleanest test available.

**Do not assume the answer.** The last two confident diagnoses on this project
-- methionine excision for the coverage gap, retention time for the
identification gap -- were both wrong, and both were settled in minutes by an
ablation. Run the equivalent here before changing the ranking: rebuild the
library with a charge-aware cap and search it, rather than reasoning about it.

---

## MS2 prediction — what the review left open

- **`predictMS2` materialises every spectrum.** 100,000 peptides measured at
  94.8 s and +1717 MiB, the payload being 14.4 M floats spread over 100,000
  separate vectors. MS2 takes charge as an input, so the proteome case is ~4.0 M
  precursors rather than the 1.25 M distinct sequences RT needs: extrapolating,
  ~63 min and ~2.2 GiB of payload in ~4 M allocations, all live at once, on top
  of the generator's existing 1.23 GB. `predictRT` returns 4 bytes per peptide.
  There is no streaming or callback form, and this should be decided before the
  intensities are wired into `LibraryGenerator`.

- **`MAX_BATCH_ROWS` is not covered by any test.** Its purpose is to bound
  memory -- without it a 2 M-peptide run peaked at 26.7 GiB -- and removing the
  chunking entirely still passes, because values do not change. Correctness
  across the 2048 boundary was verified by hand (indices 2046-2050 bit-identical
  to solo runs) but nothing holds it.

- **CCS stays in square angstroms. Decided, not open** (user, 2026-08-03).
  `predictCCS` returns collision cross-section and must keep doing so. The
  conversion to 1/K0 is deliberately *not* ODIA's: it is instrument-specific and
  is applied downstream via the Mason-Schamp equation, where the drift gas and
  the instrument's calibration are known. ODIA must therefore never convert, and
  must not acquire a `-ion_mobility_unit` option that implies it could.
  Consequence to settle: the library's `im` column is 1/K0-shaped, so predicted
  CCS needs either its own column or an explicit unit tag -- writing CCS into a
  field consumers read as 1/K0 is precisely the silent unit error this decision
  avoids.

- **The model guard's reject path is untested.** All three shipped models are
  accepted, so nothing here exercises the case it exists for -- a model that is
  none of the three. Testing it needs a synthetic ONNX file and the `onnx`
  Python package is not in the environment. The per-call guards *are* tested,
  by feeding each `predict*` the other two models.

- **Two output-shape checks are equivalent mutants.** Removing the CCS rank
  check, or the batch-size check, changes nothing observable, because the real
  model always returns the right shape. They are cheap insurance against a
  re-exported model, not something any test here can hold.

- **A failed peptide is an empty spectrum**, which `Spectrum::at()` will index
  out of bounds. `predictRT` uses NaN for the same situation; there is no
  equivalent value here, so the contract is "check `positions` first" and it is
  only documented, not enforced.

---

## `.oswpq` reading — still open after the second review

- **`getInt64` on a `uint32`/`uint64` above the signed range is an equivalent
  mutant.** Narrowing the cast wraps the value, but ids are used only as join
  keys and the wrap is bijective, so precursors and transitions still meet.
  Only `UINT64_MAX` would alias the null sentinel. Nothing here can catch a
  narrowing cast, and nothing here needs to.

- **The two null-id guards mask each other.** A precursor with no id is kept out
  of the join map, and a transition with no id never looks one up; removing
  either alone changes nothing, and only removing both attaches an id-less
  fragment to an id-less precursor. Kept as belt and braces, with the pair
  covered by a test.

- **`ChunkedColumn` is not thread-safe** and now says so. The cursor is shared
  mutable state, so two threads scanning one column return each other's rows
  rather than colliding visibly. The header advertises the class as the hot
  path for a 64-thread scan, so this needs a per-thread cursor before any
  parallel reader is written.

- **A precursor whose m/z cannot be represented is kept at m/z 0** and can never
  match a window. `invalidMzCount()` reports it -- one on the real upstream
  bundle -- but `load()` does not surface it to the caller.

---

## `.oswpq` reading — what it does not yet do

- **The multi-chunk path is never exercised against a real file.** It is the
  reason `ChunkedColumn` exists, and it is unreachable below 2 GB of characters
  in one column: Arrow's Parquet reader concatenates row groups, so a fixture
  written with one row group per row still comes back as a single chunk. The
  cursor is therefore tested against chunk layouts built directly in
  `test/tools/odia_chunked_column_test.cpp`, which catches the defect but does
  not prove Arrow splits where we think it does. A proteome-scale bundle would
  settle it; none exists here.

- **Both tables are materialised whole before conversion.** For the benchmark
  library that is 78.6 M transitions of Arrow on top of the ODIA library being
  built. Column projection removes the largest contributor (`traml_id`), but the
  peak is still roughly double what it needs to be. Reading row group by row
  group would fix it and was not done, because nothing here can measure it: the
  largest bundle available is 18 transitions.

- **Fields read and then dropped**: `traml_id`, `unmodified_sequence`,
  `transition_id`, and the `detecting` / `identifying` / `quantifying` flags.
  The first is deliberate (D3: no per-transition strings) but it means ODIA
  cannot round-trip a bundle -- the human-facing precursor id is gone, and a
  written bundle would have to synthesise one. The three booleans matter for
  OpenSWATH scoring and will have to be carried before Phase 3.

- **The census check is advisory.** A disagreement between the row counts and
  `library/metadata.json` is recorded in `Stats` and nothing acts on it. It
  should probably be a hard failure by default, since the whole point is to
  catch a truncated read before an hour of extraction, but "probably" is not
  enough to make a load fail.

- **`schema_version` is checked; nothing else in the metadata is.** The
  `fragment_type_counts` and `charge_counts` blocks are a second, finer census
  that would catch a mis-parsed `type` column, and they are ignored.

- **The spec document disagrees with the file in three places** and should be
  corrected: it does not mention that `type`/`annotation` may be empty and
  `ordinal` `-1`; it does not mention that a precursor may have no transitions;
  and it says the transition `traml_id` is the precursor's, denormalised, where
  the sample bundle carries the transition's own id there instead.

---

## Smaller things noticed in passing

- Three decoys of 1,993,956 are skipped on the human proteome — peptides with no
  suitable unmodified position near both termini. Harmless, but the count should
  be reported rather than silently dropped.
- `LibraryGenerator` holds the whole peptide-to-protein map in memory before
  emitting. Fine at proteome scale (1.23 GB peak), worth revisiting if libraries
  get much larger.
- The mzpeak fork's test suite does not compile against Arrow 23
  (`parquet_writer_test.cpp` calls `FileReader::ReadTable()` with a signature
  this Arrow lacks), so `build_mzpeak.sh` builds only the library and tools.
- `ODIAInfo -peaks` is unusable on real files until the decode issue above is
  fixed; it is opt-in and documented as such.

---

## Upstream reports worth filing

- **OpenMS**: `WITH_ONNX=ON` cannot configure as shipped —
  `cmake/FindONNXRuntime.cmake` is not on the `CMAKE_MODULE_PATH` that OpenMS
  sets (`cmake/Modules` and `cmake/Windows` only).
- **OpenMS**: `WITH_ONNX=ON` compiles and exports the PeptDeep classes but never
  installs their headers. `libOpenMS.so` carries 18 PeptDeep symbols, and 26 ML
  headers are installed — every `ML/` subdirectory except `PEPTDEEP/` and
  `ONNX/`, the two the flag exists to enable. No external project can include
  them, so the feature is unusable outside the OpenMS tree.
- **OpenMS**: no mzPeak entry in `FileTypes`, so a TOPP tool cannot declare
  mzPeak as an input format.
- **OpenMS**: `ParquetFile::getColumn()` returns only the first chunk while
  callers loop to `num_rows` — a truncating read on any column past 2 GB.
- **mzpeak**: `EnumerableProxy::Iterator` defaults a constructor taking
  `const Iterator&&`, which is not a move constructor and cannot be defaulted;
  GCC rejects it. Carried as `patches/mzpeak-0001-defaulted-move-ctor.patch`.
- **mzpeak**: `install_headers(subdir: 'mzpeak')` flattens the header tree, so
  the installed headers do not compile.
- **mzpeak-convert**: `--layout chunked` silently produces a point-layout file
  while printing its own `BUG:` diagnostics about signal arrays spilling to
  `auxiliary_arrays`.

## GPU inference: built and verified, unrun for want of a credential (2026-08-04)

`opt/env-gpu` holds a CUDA 12.9 build of ONNX Runtime **1.26.0** -- the same
version the CPU path uses, so it is ABI-compatible with what OpenMS links --
plus the cudart/cublas/cufft/cudnn the provider dlopens, which the conda
package does *not* pull in on its own (it declares only the `cuda-version`
metapackage). It is a separate prefix: `opt/env` is what the read-only OpenMS
install resolves against and must not be swapped under.

`build-gpu/` is ODIA linked against it, on Ceph so the GPU nodes see it without
a rebuild. Verified on ibminode05: `GetAvailableProviders()` returns
`CUDAExecutionProvider, CPUExecutionProvider`, against `CPUExecutionProvider`
alone for the old build. So `prefer_gpu` finally has something to find.

`build-gpu/run_gpu_bench.sh` is staged and needs no arguments.

**What blocks the run is authentication, not routing.** On ibminode05:

- `ssh -v` reports every identity file as `type -1` -- none exist. `~/.ssh`
  symlinks to `/afs/wsi/home/oliver/.ssh` and holds `authorized_keys`,
  `config`, `known_hosts` and an unrelated GitHub deploy key.
- The two keys in `authorized_keys` have no private half on this filesystem,
  and the deploy key's public half does not match either of them.
- No Kerberos ticket (`klist`: no credentials cache), so `gssapi-*` is out too.

Every hop therefore refuses at the first one: `spock` and `data` directly,
`-J hive`, and `-J sshgw` all fail on publickey. Jump hosts do not help,
because the failure is not reachability.

Unblocking needs a private key on ibminode05 whose public half is in
`~/.ssh/authorized_keys` -- the setup step the `ibmi-hpc` skill documents.
Once that exists the benchmark is one command.

## GPU inference measured on an H100 (2026-08-04)

Full human proteome (20,416 proteins -> 2,127,559 target precursors), `data`,
one H100 PCIe, against the v5 CPU run on ibminode05:

| stage | CPU (v5) | GPU | |
|---|---:|---:|---:|
| retention time | 223.5 s | **17.1 s** | 13.1x |
| fragment intensities | 2,372.5 s | **67.3 s** | 35.2x |
| collision cross-sections | 368.1 s | **34.3 s** | 10.7x |
| **total wall** | **53:23** | **3:58** | **13.5x** |
| peak RSS | 2.04 GiB | 2.15 GiB | unchanged |
| CPU used | 5,527% | 105% | ~53x less |

**Read the 35x carefully.** It is against the *old* CPU path, which took ONNX
Runtime's intra-op default. Same node, same binary, same peptides, GPU against
the *fixed* CPU path: 54,500 peptides/s against 8,953 at 32 sessions -- **6.1x**.
That is the honest hardware comparison; the rest of the 35x was the threading
bug. Note `data` was at load 370 on 224 cores, so the CPU figure is depressed.

**The phase is no longer inference-dominated.** Inference is 119 s of a 238 s
run; digestion, decoy construction and writing 7.28 GB of TSV are the other
119 s. Further model speedups now buy at most 2x on this stage.

CPU and GPU are **not** bit-identical, and should not be expected to be --
different kernels, different reduction order. With the same MS2 model on both,
over 200,000 fragments: 98.05% agree within 1e-3, 99.99% within 1e-2, and 9
fragments (0.0045%) differ by more than 0.1 -- ranking ties near the top-12
boundary. CCS agrees to a median 0.019 A^2. The determinism test's bit-for-bit
guarantee is *within* a provider, which is what it claims.

RT differed by 3.2 min median between the two libraries, which is **not** a GPU
effect: v5 used the fine-tuned checkpoint
(`rtfinetune/integrated/peptdeep_rt_dynamic.onnx`) and the GPU run used the
stock OpenMS model. MS2 and CCS used the same stock models on both, which is
exactly why only RT moved.

### Two traps worth keeping

* **GPU 0 is unusable on both nodes.** `cudaSetDevice(0)` returns error 46
  ("devices busy or unavailable") on spock and data, with both GPUs in Default
  compute mode, zero memory used and no processes attached -- and world-writable
  `/dev/nvidia*`. On `data`, `cudaSetDevice(1)` succeeds. On `spock` **both**
  devices fail, so spock currently has no usable GPU at all. ONNX Runtime
  defaults to device 0, so **`CUDA_VISIBLE_DEVICES=1` is required** on data.
  Filed for the admin.
* **The conda `onnxruntime-cpp` cuda build declares only `cuda-version`**, a
  metapackage. cudart/cublas/cufft/cudnn must be installed explicitly or the
  provider fails to load and ONNX Runtime falls back to CPU silently.
