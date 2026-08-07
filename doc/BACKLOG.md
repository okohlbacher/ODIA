# Backlog

Open items, newest first within each group. Items that need a decision or an
external fix are marked **[you]**; the rest are mine to work through.

---

## The 1/K0 anchors are in; what is still open on that axis (2026-08-06)

The 1/K0 calibration now measures at anchors instead of guessing at them, and
the gate passes on its own margin. `collectAt()` visits only the cycle blocks
that hold a scored peak group's apex; the null is the library's own decoys at
the apexes THEIR groups claimed. On S08, with 1,534 confident targets and a
rank-matched null: peakedness **12.41 against 3.75, a 3.31x margin** where 1.25x
is required and where the blind probe managed 1.09x. 23.7% of the mean squared
1/K0 error removed out of fold. The fitted curve matches DIA-NN's observed 1/K0,
which the probe never sees, to within a couple of milli-1/K0 per m/z bin.

End to end on the frozen discriminant: **52.27% -> 53.17% (1393 -> 1417 of
2,665) at a matched, re-measured 1.00% entrapment false rate**; 54.15% at the
frozen threshold, where the entrapment null has itself moved to 1.23%. The gain
sits exactly where the lever was aimed -- the 388 outside-the-cell misses go
15.7% -> 23.2%, and every other bucket moves by less than a point.

What is still open:

- **The stage costs a second sequential decode: 22,440 of 32,210 spectra,
  ~450 s, roughly doubling a single-pass extraction.** Anchors are spread over
  the whole gradient, so "only the blocks with an apex" is 70% of the run. The
  saving available is that within a visited block only the ONE isolation window
  each anchor lives in is needed -- about 2 of 32 spectra per cycle. That is a
  `SpectrumSource::peaks` call pattern change, not an algorithm change, and it
  is worth ~10x if the reader serves sub-ranges without re-decoding a row group
  per spectrum. It does not today (see the mzPeak item below), so this is
  blocked behind the same decode fix everything else is.

- **The gain at 1% is the smallest point on its own sweep**, and the 1%
  entrapment quantile rests on 26 events. Each arm at its own entrapment
  quantile: +1.02 points at 5%, +1.51 at 2%, +0.90 at 1%, +1.16 at 0.5%, +1.62
  at 0.2%. Nothing here says 1% is special; it is where the operating point was
  fixed. Worth re-measuring on a second run before the +0.90 is quoted as the
  number.

- **66 precursors are LOST at the frozen threshold against 116 gained.** A
  recentred window that moves off a precursor whose library 1/K0 was already
  right is the obvious mechanism, and it is not yet measured. The out-of-fold
  scatter says the correction is right on average; it does not say it is right
  for those 66.

- ~~**The width lever is still untaken.**~~ **TAKEN AND CLOSED (2026-08-07):
  widening is worse, keep +/-0.025.** Measured on S08/lib_targets with the
  mobility slope fitted, identifications at 1% FDR: **0.025 -> 1232, 0.035 ->
  1215, 0.050 -> 1167.** Monotonic, and 0.050 gives back nearly the whole gain
  from the slope (1165 with no mobility-linear term at all).

  The residual argument that pointed the other way -- corrected SD 0.025, so
  +/-0.025 is ~1 sigma and rejects 28.3% of target anchors -- was right about
  the width and wrong about the remedy. Those 28.3% are the interference-prone
  tail; admitting them costs more than they bring, which is exactly the ~8.5x
  same-window mobility this window exists to exclude. The fix for a trend is to
  model the trend, not to widen the gate until the trend fits through it.
  Recorded at `ChromatogramExtractor.h`'s `precursor_im_window`.

- ~~**The production path has not been measured end to end on S08.**~~ **RUN
  (2026-08-07).** `-ion_mobility_calibration anchors` over the two-pass
  workflow, anchors from ODIA's own pass-1 scorer, S08/lib_targets: pass 1
  offers 729 anchors at q<=0.01 against a null of 729 decoys, the gate passes,
  and pass 2 reports **1232 at 1% FDR**. The stage works on its own anchors
  without an external scorer.

  Two caveats it exposed, both open below: charge 3 is left uncorrected, and
  the run bootstraps its own iRT map (no `-irt_slope`/`-irt_intercept`), so it
  warns that pass 1 extracts at approximately the wrong retention times. Both
  arms of the slope comparison shared that, so the differential is clean, but
  the absolute number is not a tuned-iRT number.

- **Charge 3 is never corrected on S08/lib_targets: 93 anchors against a
  minimum of 120.** The mobility slope is fitted per charge on purpose --
  another charge's offset is not this charge's answer -- but that means the
  whole correction, slope included, silently does nothing for charge 3. Charge
  3 is ~20% of the library and its measured residual is the WORST (offline:
  constant +0.0335 and slope -0.113, against +0.0019/-0.096 for charge 2), so
  the charge that most needs the correction is the one that cannot reach the
  anchor count for it. Options, in order of preference: pool charges for the
  SLOPE only while keeping the constant per charge (the slope is a property of
  the CCS->1/K0 conversion, which is shared, so this is physically justified in
  a way that pooling offsets is not); or lower the threshold with a widened
  confidence requirement. Do not simply lower the minimum.

- **The mobility slope is confirmed but its magnitude is not pinned down.**
  Fitted -0.0709 and -0.0748 per 1/K0 on the run's own 371-378 anchors, against
  -0.096 (charge 2) and -0.113 (charge 3) estimated offline on ~22k anchors,
  and -0.128/-0.109 from outlier-resistant bin medians. Same sign and order
  throughout, but a ~30% spread. Since this is a scale error in the
  `ccs_to_mobility.py` coefficient (1037.1902), the right fix is upstream and
  one-off: fit the coefficient properly against Mason-Schamp rather than
  re-deriving a slope per run. Worth doing once the Astral/second-run
  measurement says whether the slope is instrument-stable.

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

## Quantise chromatogram intensity to one byte

The shared retention-time axis landed; the intensity encoding did not. Storing
each point as a uint8 log step against a **per-transition** float32 scale is a
further 4x on the point arrays.

Recorded because the objection raised against it was wrong: the concern was
that quantisation would cost precision in the weak-signal regime currently
under investigation. It does not, because the scale is per transition rather
than global -- 255 log steps over one transition's own two-decade range is
~2.7% per step, well below the noise on any point that matters. The reference
implementation uses exactly this (`odia_chromstore.h`, `Quantised8Log`).

Reserve 0 for exact zero rather than for the smallest representable value: an
absent measurement and a very small one are different claims, and a
chromatogram is mostly the former.

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

  **Done, 2026-08-05** (`include/odia/MassCalibration.h`). The run's own
  fragment mass error is now measured before extraction and the window is
  centred on it, gated so that a run with no measurable error is a no-op. The
  finding that mattered was not the width but the CENTRE: on S08 the axis is
  about -10 ppm out, and a window narrowed about zero throws away more than it
  saves. `-fragment_ppm` therefore no longer has a fixed default at all -- it is
  10 ppm when the calibration centres the window and 15 when it cannot.
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

## A comprehensive score inventory, and a learned score for SELECTION (2026-08-07)

Two connected pieces of work. We have 15 sub-scores; DIA-NN has **110**, and
OpenSWATH/pyProphet's published set is ~20. The 15 were chosen conservatively --
"a sub-score computed from a placeholder is worse than an absent one" -- and that
was right while extraction was broken. It is now the binding constraint: the true
peak is available for **97.4%** of precursors and we rank it first for **75.4%**,
so ~22 points sit in discrimination we are not computing.

### Part 1: the score inventory

Test every candidate below on **both** benchmark files (S08 diaPASEF, Astral),
against each file's DIA-NN confident set, measuring **best-ranked-right** -- not
on-RT, which counts TSV rows and is inflated by ~1.5 candidates per precursor
sharing a broadcast q-value.

**From DIA-NN (`diann.cpp:780` enum, handoff section 6.1). Cheap and applicable:**

| DIA-NN name | what it is | why it may matter here | cost |
|---|---|---|---|
| `pTightCorrOne/Two` | Σ fragment correlations at 0.45x and 0.20x tolerance | interference survives a loose tolerance and dies at a tight one; this is the single cheapest orthogonal signal we lack | re-match at 2 extra tolerances |
| `pAcc+0..5` | per-fragment `(\|obs-exp m/z\| / tol) x correlation` | the mass-error feature we keep deferring, but PER FRAGMENT and correlation-weighted | needs the ppm retention already planned |
| `pBestCorrDelta` | `pTimeCorr - best_corr_sum` for this precursor | we compute CANDIDATE_MARGIN, which is close; theirs normalises against the run | free |
| `pTotCorrSum` | `log(pTimeCorr / (total_corr_sum + 1))` | normalises a candidate against how much correlation the whole run offered | free |
| `pResCorr`, `pResCorrNorm` | correlation of fragments BEYOND the top 6 | we use all 12 equally; splitting top-6 from the rest is a real distinction | free |
| `pShape+0..4` | elution profile in 5 symmetric bins | peak shape as a vector rather than a width scalar | free |
| `pSig+0..5` | per-fragment share of integrated signal | we have YSERIES_SCORE and INTENSITY_SCORE; this is finer | free |
| `pCorr+0..11` | RAW per-fragment correlations, not summed | a summary hides which fragment disagreed | free |
| `pMinCorr` | correlation against 3-point-minimum traces | a baseline-insensitive variant | cheap |
| `pShadow`, `pShadowCorr` | correlation with -1.00335 Th shadow traces | an isotope shadow that a real peptide has and interference does not | needs a second extraction offset |
| `pHeavy` | correlation with +1 isotope traces from the NEIGHBOURING window | same idea across windows | needs cross-window extraction |
| `pdRT`, `pRT` | `sqrt(\|dRT\|/span)`, and position in gradient | we have RT_DELTA raw; the sqrt and the position are different shapes | free |
| `pMz`, `pCharge`, `pLength`, `pFrNum`, `pMods`, `pAAs+0..19` | peptide properties and AA counts | 26 features that need no signal at all. DIA-NN marks these `p_none` -- NN-only, never used by its linear classifier. Test them, but expect them to need the NN | free |

**Requires MS1 extraction, which we do not do at all (a whole sub-project):**
`pMs1TimeCorr`, `pMs1TightOne/Two`, `pMs1Iso*`, `pMs1Ratio` -- MS1 correlation and
isotope agreement. DIA-NN devotes ~8 features to it. Note `MS1PeakSelection` is on
by default there, i.e. MS1 participates in DETECTION, not only scoring.

**Not applicable:** the Q1 block (`pQLeft/pQRight/pQPos/pQNFCorr/pQCorr` x3) is
Scanning SWATH only.

**From OpenSWATH / pyProphet (mProphet lineage), for the ones we lack:**
- `xx_swath_prelim_score`, `bseries_score` (we have y-series only)
- `massdev_score` and `massdev_score_weighted`
- `isotope_correlation_score`, `isotope_overlap_score` (needs MS1)
- `norm_rt_score` (we have the raw delta)
- `elution_model_fit_score` -- fit an EMG to the peak and score the residual;
  genuinely orthogonal to every correlation we compute
- `sn_ratio` variants beyond our single LOG_SN
- library dot-product and manhattan variants beyond our two

### Part 2: a learned score for SELECTION, not just for FDR

**This is the structural gap, and it is worth more than any individual score.**
Today the discriminant runs ONCE, after candidates are chosen: the picker emits
candidates by correlation, sub-scores are computed, and the classifier separates
target from decoy. Nothing learned ever feeds back into WHICH candidate wins.

DIA-NN's `cscore = Sum_i w_i * score_i` (`Precursor::seek`, `diann.cpp:8121`) is
used to SELECT the winning peak group, with the weights re-fitted every iteration
and features gated by `min_iter_seek`. Selection and classification are the same
learned function.

We now have `PeakGroupScorer::refit` and a convergence loop, so the machinery
exists. The missing piece is using the fitted discriminant to re-select the best
candidate per precursor between rounds, not only to score the one already chosen.

Watch for: this closes a loop that can inflate its own confidence -- the
classifier would be selecting the rows it is then trained on. DIA-NN's answer is
the feature schedule plus `check_weights` sign clipping. Ours must be an external
check (best-ranked-right against DIA-NN) rather than an internal one, because
every internal number moves together.

### How to test without fooling ourselves

Leave-one-out over the final set, on both files, reporting best-ranked-right.
Two of our existing 15 are known to be worthless: USABLE_FRAGMENTS was a constant
12.000 for a week, and IM_DELTA is currently all-NaN because
`Options::observed_im` was never wired to MobilityCalibration's output. Any new
score has to beat that bar, and the ablation is what proves it.

## Memory, parallel occupancy and runtime across all phases (2026-08-07)

Every performance number this project has is a spot measurement of whichever
phase was under suspicion at the time. There has never been one profile that
covers the whole pipeline on both benchmark files, and the gaps have already
cost us: a 3.6x memory win sat in a hardcoded `BLOCK = 1024` for weeks, and a
1.02 GB index was per-transition when every value in it was per-precursor.

What to measure, per phase and end to end, on **both** files:

**Phases:** library load -> decoy generation -> (library generation if from FASTA,
which is ONNX-dominated and separately GPU-capable) -> mass calibration ->
mobility calibration -> pass-1 extract+score -> RT fit -> pass-2 extract+score ->
refine rounds -> output writing.

**Per phase:**
- wall time, and its share of the total;
- peak RSS attributable to the phase (RSS is monotone, so use the delta plus a
  heap profile where the delta is ambiguous -- `footprintBytes()` was undercounting
  the string-arena map by 3x until it was checked against an RSS delta);
- **parallel occupancy**: threads requested vs mean threads actually running.
  This is the number we have never measured and the one most likely to be
  embarrassing. Known: mzPeak decode holds a mutex across the decode itself
  (`src/util/parquet.cpp`, deliberately -- concurrent readers would corrupt each
  other through one seek-and-read handle), and decode is ~98% of a scoring run on
  the tiny library. Measured indirectly: 16 threads 18:46 vs 96 threads 9:51, a
  1.9x for a 6x thread increase. Sublinear, unexplained in detail.

**Scaling curves, not single points:**
- threads: 1, 4, 16, 32, 64, 96, 224 (both nodes have 224 cores);
- library size: the existing lib_tiny / n5k / n100k / lib50k / 4.26M ladder;
- file: S08 (32,210 spectra, 12.75 GiB, diaPASEF) vs Astral (303,701 spectra,
  3.10 GiB, no IM, no co-packed frames). Astral is 9.4x the spectra and took
  27:55 against S08's ~10 -- 2.8x, which is either good news about per-spectrum
  cost or bad news about something else, and nobody has looked.

**Specific open questions this would answer:**
- Where does the ~2.9 GiB decode floor actually come from now? It was 10.3 GiB,
  `-decode_block 256` took it to 2.9, and the composition was never re-derived.
- Does the co-elution picker's per-position pairwise correlation scale badly with
  spectrum count? It does far more arithmetic per cycle than a local maximum of a
  sum, and has only ever been measured on the small file.
- What is the actual thread ceiling, and is it the decode mutex, memory
  bandwidth, or the serial scoring on the extractor's thread (noted as unmeasured
  when the sliding window landed)?
- At 4.26M precursors, pass 1 wants ~274 GiB. Sampling it is blocked on the
  classifier's viability floor (below ~5,000 precursors it fails, silently,
  measured: 635/1820/1211/0 identifications from four equal-sized subsets). What
  is that floor exactly?

**Deliverable:** one table per file, phases x (wall, share, peak RSS, mean
occupancy), plus the three scaling curves. Publish it as doc/14 and keep it
current -- the reason this backlog entry exists is that every previous
performance claim was made from a measurement taken for a different purpose.

### Part 1b: compute the EXPENSIVE scores too, then let SHAP decide

Addendum. The table above sorts candidates by cost and implicitly suggests doing
the cheap ones first. That ordering is a trap: it selects features by
implementation convenience rather than by information, and we have no evidence
the cheap ones carry the signal. Compute them ALL, including:

- **MS1 block** (~8 features): `pMs1TimeCorr`, `pMs1TightOne/Two`, `pMs1Iso*` at
  three tolerances, `pMs1Ratio`. Needs MS1 extraction, which we do not do at all.
  This is the single largest missing capability -- DIA-NN also uses MS1 in
  DETECTION (`MS1PeakSelection`, default on), so building it may pay twice.
- **Shadow and isotope traces**: `pShadow`, `pShadowCorr+0..5` (-1.00335 Th),
  `pHeavy` (+1 isotope from the NEIGHBOURING isolation window). Each needs an
  extra extraction offset; `pHeavy` needs cross-window extraction, which our
  one-precursor-one-window assignment currently forbids.
- **Tightened-tolerance re-matching**: `pTightCorrOne/Two`, and the `p_fit` class
  DIA-NN only enables with `--tight-mass-acc-aux-for-cal`.
- **`elution_model_fit_score`** (OpenSWATH): fit an EMG per peak group and score
  the residual. Expensive per candidate and orthogonal to every correlation.

### Score selection by ML, with SHAP -- and why the obvious way is wrong

Fit the GBT over the full candidate set and rank features by **mean |SHAP|**,
then keep the ones that carry the model. Straightforward, and there are two
failure modes to design around:

**1. SHAP measures what separates targets from DECOYS, not what finds the right
peak.** Those are different objectives and we have already measured them coming
apart: off-RT targets are statistically identical to decoys on `library_corr`
(median -0.036 vs -0.032), yet the classifier still admitted them, because it
learned signal-PRESENCE features that separate a real-but-misplaced target from a
shuffled decoy. A SHAP ranking would have rewarded exactly those features. So the
selection metric must be the **external** one -- best-ranked-right against each
file's DIA-NN confident set -- with SHAP used to generate the candidate ordering,
not to make the decision.

**2. Importance is dataset-specific, and our two files differ structurally.**
S08 is diaPASEF (ion mobility, co-packed two-windows-per-frame); Astral has
neither. `IM_DELTA` is definitionally worthless on Astral. Anything derived from
frame packing is S08-only. So the procedure is:
  - fit and rank on S08 and on Astral **independently**;
  - report both rankings side by side, plus the rank correlation between them;
  - keep the union of what is decisive on either, not the intersection -- a
    feature that only works on one instrument class is still worth having, gated;
  - treat a large ranking disagreement as a finding about the instruments, not as
    noise to average away.

A third file would make this much stronger, since two points cannot distinguish
"instrument-specific" from "this particular run". `12_80` (SCIEX SWATH, no IM,
different vendor) is on disk and unused.

**Guard against the loop closing on itself.** If SHAP selection feeds a
discriminant that also SELECTS candidates (Part 2), the model chooses the rows it
is trained on and every internal metric agrees with itself. Hold out precursors,
not rows, and judge on the external metric only.

## Fine-tuning belongs inside the run's loop, not in the library on disk (2026-08-07)

**Measured regression.** Library v4 scored **37,583** confident precursors in
DIA-NN; v5 scored **35,556**. v4 BEAT DIA-NN's own library (37,247). The
generation recipes are byte-identical except one path:

    v4:  -rt_model .../rtfinetune/onnx/peptdeep_rt_dynamic.onnx
    v5:  -rt_model .../rtfinetune/integrated/peptdeep_rt_dynamic.onnx

The `integrated/` model's own `rt_provenance.json`:

    "tuned_from": "search_diann/report.parquet",
    "peptides": 500, "epochs": 40,
    "warning": "This model is specific to the run it was tuned on and must not
                be reused across runs or gradients."

500 peptides at 40 epochs, and the file says not to reuse it. The script passed
it anyway, and it cost 2,027 precursors -- 5.4%.

A guard now warns when `-rt_model` has a provenance file carrying a warning.
That is a plaster; the design below is the fix.

### The intended architecture

Library generation moves INSIDE the calibration loop, and the fine-tuned model
is scoped to one run and never written anywhere another run can find it:

1. generate the first library with the **stock** model;
2. search, calibrate RT (and mass, and mobility);
3. fine-tune the RT model **on this run's own confident identifications**;
4. re-predict the library with the tuned model;
5. re-calibrate;
6. iterate 3-5 until identifications stop improving;
7. discard the tuned model with the run.

The distinction that matters: a model tuned on run A and used to build a library
for run B imports A's gradient into B's predictions. Tuned on run A and used
within run A it is legitimate refinement -- the model has only seen data from
the run it is predicting for.

### What this needs

- Library generation callable mid-run, not only as `-stop_after library`. Today
  it is a separate invocation writing a TSV.
- Fine-tuning in-process, or as a subprocess the loop drives. Today it is the
  standalone `rtfinetune/` pipeline (`prep.py` + training, ~413 s for 500
  peptides / 40 epochs on one run).
- 500 peptides / 40 epochs is almost certainly overfitting. Inside a loop with
  the run's own identifications there are thousands available -- 1,018 anchors
  at q<=0.05 on a 2,665-precursor subset, so a full library gives far more.
  Tune the count and epochs against held-out identifications from the SAME run,
  which is the check the standalone pipeline never had.
- A convergence criterion that is not the identification count, which swings
  ~10% between refits on fixed input.

### The cross-run trap this closes

`odia_v5.tsv` carries iRT from a model tuned on S08's DIA-NN results. Every
recovery number measured on S08 with that library has partly seen the answers.
It is not the reason the library underperforms -- it underperforms because the
tuning was bad, not because it leaked -- but both are reasons the tuned model
must not be a persisted artefact.

### Confirmed 2026-08-07: v6 reproduces v4 and beats DIA-NN

Regenerated with the stock RT model, via `odia_v4_ft.sh` with only the output
name changed (sourcing env.sh, applying ccs_to_mobility, calling search.sh):

  odia_v6  stock RT model, CURRENT code   37,596
  odia_v4  stock RT model                 37,583   (+13, 0.03%)
  DIA-NN own                              37,247   (v6 +349, +0.9%)
  odia_v5  run-tuned RT model             35,556   (v6 +2,040)

Settles three things. No code regression -- everything changed this week left
library generation untouched. The v5 loss was ENTIRELY the run-tuned model,
all 2,040 recovered by the path swap. And our generator beats DIA-NN's own
library, reproducibly, twice.

Generation is also ~5x faster than when v4 was made: RT 215.7 -> 44.1 s, CCS
364.8 -> 89.7 s, from the data-parallel ONNX sessions (2a4c73a, 43af3fc).
Both stages use `inferenceSessions(sessions)` with 1 intra-op thread.

**Process note, because it cost four wrong numbers.** Getting here took four
attempts, each failing because I rebuilt the pipeline from memory rather than
running the script sitting beside the data:

  1. no --fasta            -> 20,065, reported as a "53.9% library gap" that
                              does not exist; it drove a 2x2 cross-analysis
                              and a plan step before the correctly-configured
                              search turned up already on disk
  2. no ccs_to_mobility    -> 14,408; nearly read as refuting the diagnosis
  3. hand-written DIA-NN   -> search.sh exists so the settings are identical
  4. no env.sh             -> CCS all-zero, announced as "a real code
                              regression in HEAD" on the strength of a grep
                              matching "CCS model available" inside the string
                              "No CCS model available"

Every one produced a plausible number that was reasoned from before being
checked. **When comparing against a prior result, run the prior result's
script** -- not an equivalent command. `odia_v4_ft.sh` is four lines and would
have been right the first time.

## The scorer has never been run against a realistic library (2026-08-07)

**Measured, S08, `v6_50k.tsv` (a stride sample of our own v6, no run-tuned iRT,
neither tool privileged):**

    DIA-NN   738 of 50,000 at 1% FDR   1:43
    ODIA       0                      36:43, 9.09 GiB

The true positive rate is **~1.5%** -- a whole-proteome library is almost
entirely peptides absent from the sample. DIA-NN gets 738/50,000 here and
37,596/2,127,559 (1.8%) on the full library, so the sample is representative.

ODIA's guards fired correctly: 1,639,188 peak groups, **819,353 target and
819,835 decoy (1:1)**, then `pass 1 identified nothing at 1% FDR`. The positive
class is ~98.5% noise, so the semi-supervised loop has no separable seed at
iteration 0 and never ignites.

**Which picker criterion rejects most, instrumented on that same run**
(`PickerRejects`, over 107,634,884 scan positions -- note criteria are counted
independently, so they sum to more than the positions):

    not a local maximum        189,557,902
    below min_corr_score        50,689,548
    below apex_evidence         36,856,546
    fewer than 2 fragments      30,191,757
    reference trace zero        11,525,195
    outside max_corr_diff            25,997
    then min_fragments_at_apex        75,479 candidates dropped
    -> 19,150 precursors yielded NO candidate at all, 0 identified at 1% FDR

The local-maximum test dominates by ~4x over anything else, and
`max_corr_diff` is doing essentially nothing (26k, a 7000x smaller effect than
the leader) -- so `-max_corr_diff` is a tuning knob with no leverage here and
should not be spent effort on. **The correlation thresholds together
(min_corr_score + apex_evidence = 87.5M) are the second force**, and those are
the ones that encode "is this a co-eluting peptide", so they are where a
realistic-library fix has to be careful: loosening them re-admits the 98.5%
noise the FDR then cannot separate.

Note this does NOT implicate the ion-mobility window: the mobility work of
2026-08-07 (see the 1/K0 section) found +/-0.025 already optimal, so the
picker's loss is not a window-width problem.

**This invalidates the regime, not the work.** Every recovery figure this week
-- 75.4% best-ranked-right, 89.1% precision, 97.4% availability -- used
`lib_targets.tsv`, 2,665 precursors drawn from DIA-NN's own confident set, i.e.
~100% true positives. Those measure ranking quality among true positives, which
is real. They say nothing about a library anyone would search.

It also kills the "viability floor below ~5,000 precursors" diagnosis: this run
had 50,000. The governing quantity is the true positive RATE.

~~**The fix is DIA-NN's batching**~~ -- **NO. MEASURED AND REFUTED
(2026-08-07).** This was my inference from handoff 7.1, not a measurement, and
it is wrong. Batching cannot bootstrap from a subset that itself yields nothing,
and that is exactly what a subset yields.

Simulated at the real regime -- 10,000 targets, 10,000 decoys, a planted true
positive rate, the actual `scoreSemiSupervisedLDA` with GBT. A RANDOM subset of
2,000 (same rate, a tenth of the rows) reports **zero in all eight
rate x separation conditions**, including the one where the full 20,000 reports
676. Fewer rows means fewer decoys, coarser q resolution and a weaker fit.
Subsetting is uniformly worse, never better.

(The first version of this probe drew a PREFIX rather than a random subset, and
the true positives are planted at the front, so the "subset" was a pre-selected
library and every subset number was inflated. That is the ranking-vs-search
confusion below, reproduced inside the very experiment meant to study it.)

**The constraint is discrimination power, and the lever is FEATURE COUNT.**
Splitting ranking from certification at a 1.5% rate and 3 sigma: the top 150
targets by the learned dscore are **48.7% genuinely true against a 1.5% base
rate** -- a 32x enrichment, so the ranking works -- while q there is 0.38
against an actual FDP of 51%, so the estimator is roughly HONEST too. Neither
is broken. At a 1.5% prior, 1% FDR needs a likelihood ratio near 6500:1, and
that is simply more evidence than a handful of features carries.

Feature count crosses it. Same regime, 3 sigma, features as independent noisy
looks at the truth (the optimistic case -- real sub-scores are correlated, so
real counts must be HIGHER than this):

    features   top-150 purity   reported at q<=0.01
       4           52.0%              0
       8           72.0%              0
      15           84.7%              0        <- what ODIA has
      30           87.3%            104
      60           95.3%            116
     110           94.7%            137        <- what DIA-NN has

The cliff is between 15 and 30, and ODIA carries exactly 15 sub-scores. That is
a direct explanation for 0 against DIA-NN's 738 on the same file, and it makes
the score-inventory item below THE priority rather than a nice-to-have. At 2
sigma even 110 features report nothing, so per-feature strength matters as much
as count -- more good features, not merely more.

Probes: `scratchpad/batch/{probe,diag,feat}.cpp`, header-only against
`odia/scoring/lda.h`, seconds to run.

`test/tools/odia_entrapment.cpp` does NOT reproduce this regime: it plants a
33% true positive rate, where the real library is 1.5%. It is a valid FDR
calibration test and a poor model of the collapse. A low-rate case belongs in
it.

**Benchmark hygiene:** always state what fraction of a library's targets are
actually present. A library built from the comparator's confident set is a
RANKING benchmark; a proteome-scale library is a SEARCH benchmark. We have only
ever run the first.
