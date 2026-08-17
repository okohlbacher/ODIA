# DIALibraryGenerator: plan

A separate TOPP tool that turns a FASTA into an in-silico DIA library, with the
whole configuration in one JSON file that is embedded in the output.

## Why separate at all

Library generation is a DIFFERENT job from searching, on a different cadence: a
library is a cross-run artefact built once and reused, a search is per-run. They
are currently welded together in `OpenDIAlyzer` behind `-out_lib`, which is why
the library's own parameters are scattered across ~15 CLI options and why the
alkylation state -- the single defect that cost a quarter of the search space --
was hard-coded and unreachable (`LibraryGenerator.h:52`, doc/27).

## Scope

- `-in` FASTA, `-config` JSON, `-out` `.parquet` or `.tsv`.
- Everything content-affecting comes from the JSON. No content parameter is a
  bare CLI flag, so a library cannot be built with settings nobody wrote down.
- The JSON is EMBEDDED in the Parquet output, so a library always carries the
  recipe that produced it.

## What is reused, not rewritten

The tool is a thin driver. All of this already exists in `odia_core`:

| step | existing code |
|---|---|
| digest + precursor enumeration | `LibraryGenerator::generate` |
| iRT prediction | `LibraryGenerator::predictRetentionTimes` -> `PeptDeepPredictor` |
| fragment intensities | `LibraryGenerator::predictFragmentIntensities` |
| CCS | `LibraryGenerator::predictCollisionCrossSections` |
| Parquet write + KV metadata | `DIANNLibraryFile::storeParquetCompact` |
| fingerprint | `DIANNLibraryFile::Fingerprint` |

Target: the tool itself is a config parser plus a call sequence. If it grows a
second algorithm it is in the wrong file.

## The JSON config

```json
{
  "enzyme": "Trypsin",
  "missed_cleavages": 1,
  "peptide_length": [7, 30],
  "precursor_charges": [1, 2, 3, 4],
  "precursor_mz": [300.0, 1800.0],
  "fragment_mz": [200.0, 1800.0],
  "fragment_charges": [1, 2],
  "fragments": [4, 12],
  "fixed_modifications": ["Carbamidomethyl (C)"],
  "variable_modifications": [],
  "max_variable_modifications": 1,
  "n_terminal_methionine_excision": true,
  "decoys": "mutate",
  "rt_model": "",
  "ms2_model": "",
  "ccs_model": "",
  "instrument": "QE",
  "nce": 30.0
}
```

`fixed_modifications` defaults to `["Carbamidomethyl (C)"]` -- the common case --
but it is now WRITTEN DOWN and embedded, which is the point. An empty list means
no fixed modification, which is what the Astral benchmark actually needed.

## Embedding, and the compatibility constraint

`storeParquetCompact` already attaches `odia.fingerprint`,
`odia.target_fingerprint`, `odia.layout`, `odia.fasta_sha`, `odia.fasta_bytes`,
`odia.params` as Arrow schema key-value metadata. The plan adds ONE key,
`odia.config_json`, holding the canonical config verbatim.

**Two hard constraints:**

1. **ODIA must read the output unchanged.** Extra KV keys are ignored by Arrow
   readers, so this is additive -- but it must be VERIFIED, not assumed.
2. **The fingerprint must match what `OpenDIAlyzer` computes**, or the library
   cache misses and the tool silently regenerates. `OpenDIAlyzer.cpp` builds
   `;len=..;mc=..;z=..;pmz=..;fmz=..;fz=..;frag=..;varmod=..;nme=..;rdc=..;fixmod=..;rt=..`
   -- and `fixmod` was only added today (doc/27). The new tool must emit the
   IDENTICAL string, or the two tools disagree about what "the same library"
   means. `DIANNLibraryFile.cpp:927` already warns that metadata written in the
   wrong place is silently dropped and every cache lookup then misses.

## Predictors

**AlphaPeptDeep via ONNX** is the implemented path (`PeptDeepPredictor`, already
working, already tested against a Python reference in `test/`).

**DeepLC is NOT implemented in this pass.** The deep-research round found only a
nominal, unverified `tf2onnx` path for it and no ONNX export in wide use. Adding
a config key for a predictor that cannot be run would be a lie in the schema.
The config carries `"rt_model": ""` (a path), so a DeepLC ONNX can be dropped in
later WITHOUT a schema change if one materialises -- that is the extension point,
and it costs nothing now.

## Testing

- Round trip: write Parquet, read with `DIANNLibraryFile`, assert precursor and
  fragment counts and a sample of masses.
- **CAM on/off** on a small FASTA: assert that cysteine-containing precursor and
  fragment masses differ by exactly `57.021464 * nCys / z` and that non-cysteine
  rows are byte-identical. This is the regression test for doc/27's root cause.
- Config round trip: the embedded `odia.config_json` parses back to the input.
- Fingerprint parity: same FASTA + same config through both tools produces the
  same fingerprint string.
- `ctest` stays green (67/67 today).

## Deliverable

A library for the reviewed UniProt human FASTA at the default config, written as
ODIA-compatible Parquet with the config embedded.

## Open questions for review

1. Is one added KV key genuinely safe for every reader in the project
   (`DIANNLibraryFile`, `OSWPQLibraryFile`, the cache path)?
2. Should `odia.params` be DERIVED from the JSON so there is one source of
   truth, or kept separate for backward compatibility with existing libraries?
3. TSV flavour: DIA-NN-compatible column names, or ODIA's own? A TSV that cannot
   be read back is a dead end.
4. Should `OpenDIAlyzer -out_lib` be deprecated in favour of this tool, or kept?
   Two code paths producing libraries is how they drift apart.
