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

- **[you] mzPeak peak decoding is ~1000× too slow to extract with.** ~277 ms per
  spectrum against ~0.27 ms for OpenMS parsing mzML. Full requirements in
  `03-mzpeak-streaming-requirements.md` (R1, R2, R5 are the blocking set).
  Phase 2 cannot start until this moves.

- **[you] `library_intensity` `float64` vs `float32` in `.oswpq`.** The
  predecessor's `float32` change is a breaking on-disk change that was never
  upstreamed, and we cannot apply it without modifying OpenMS. Proposal on the
  table: read both, write upstream-compatible `float64` by default with
  `float32` behind an option. **The reader now does accept both**, verified
  against fixtures of each width, so only the writer is still blocked on this.

- **[you] GPU access.** `spock` and `data` both refuse: `Permission denied
  (publickey)`. The CUDA path in the predictor cannot be tested until then, and
  the environment currently carries the CPU `onnxruntime` build.

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
