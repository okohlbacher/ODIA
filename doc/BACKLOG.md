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
  `float32` behind an option.

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

- ONNX prediction: iRT, MS2 intensities, CCS — with CUDA attempted and CPU
  fallback, and ODIA's own `mod_x` encoder.
- Replace the placeholder fragment ranking. `LibraryGenerator` currently caps
  fragments by descending m/z because there are no predicted intensities yet;
  prediction should drive the selection.
- `.oswpq` read and write (gated on the `float64`/`float32` decision above for
  writing; reading can proceed and should accept both).
- Phase instrumentation: wall, CPU, RSS, `mallinfo2`, and node load per run,
  charging un-phased time to the preceding phase.
- mzPeak chromatogram writer, using the fork's writer.

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
