# Technical constraints discovered while bootstrapping

Facts established by building the stack, not design decisions. They bound what
ODIA can be built on, so they should be settled before the scope discussion.

## 1. mzpeak's high-level API is not yet sufficient for DIA extraction

`MzPeak::Spectrum` exposes exactly three things:

```c++
const std::vector<double>& mz() const;
const std::vector<float>&  intensity() const;
uint8_t                    ms_level() const;
```

`MzPeak::Metadata::Spectrum` adds only `ms_level()` and `delta_model()`.

Everything a DIA extraction actually keys on — **retention time, precursor
isolation window, ion mobility** — is absent from the high-level interface. It
is present in the archive: `mzpeak-convert` reports the member tables

```
spectra_metadata, spectra_metadata_scans,
spectra_metadata_precursors, spectra_metadata_selected_ions
```

and these are reachable through the low-level interface, `Index::parquet(const
Schema::File&)`, which hands back a `Util::Parquet` wrapping the Arrow/Parquet
reader. mzpeak's own `bin/mzp-inspect.cpp` is the worked example.

Consequences:

- ODIA cannot treat mzpeak as a drop-in spectrum source. It needs its own
  access layer over the Parquet metadata tables, or upstream must grow the
  high-level API. This is the single largest unknown in the build.
- That layer is also where any column-projection and row-group-pruning strategy
  belongs, which is precisely where mzPeak should beat mzML. `Util::Planner`
  and `Util::Projection` exist for this and are worth reading before designing
  anything.
- mzpeak's README states the API is **not stable** ("no stability is guaranteed
  at this point"; freeze targeted for summer 2026). ODIA pins a commit.

## 1. Reader: use okohlbacher/mzpeak-openms, not OpenMS/mzpeak

ODIA builds against the fork, pinned at `587a4fb` (branch `writer_test`). It
carries three things upstream does not have, all of which ODIA needs:

- **It reads the split-metadata (v2) layout**, which is what
  `mzpeak-convert >= 0.7.0` writes and therefore what every example archive
  uses. Upstream trunk throws on all of them (§1a).
- **It exposes per-spectrum retention time, precursors/isolation windows and
  ion mobility.** Upstream's `Spectrum` exposes only m/z, intensity and MS
  level, which cannot drive a DIA extraction (§1b).
- **It has a writer**, which is how ODIA will emit mzPeak chromatograms.

### Verified against OpenMS on `12_80`

`ODIAInfo` reports the same quantities through both readers, deliberately, so
the two can be cross-checked. They agree exactly:

| | mzPeak (fork) | mzML (OpenMS) |
|---|---|---|
| spectra | 13,009 | 13,009 |
| MS1 / MS2 | 1,083 / 11,926 | 1,083 / 11,926 |
| RT range | 30.1 – 929.995 s | 30.1 – 929.995 s |
| isolation windows | 11 | 11 |
| cycles | 1,085 | 1,085 |
| **wall / peak RSS** | **0.29 s / 95 MB** | 3.76 s / 412 MB |

So for the metadata pass — which is what run indexing and window/cycle
detection need — mzPeak is **13× faster and uses 4.3× less memory** than
parsing mzML. That part of the thesis holds.

### Two live gaps in the fork

1. **Peak decoding costs ~277 ms per spectrum.** Measured over 1,000 spectra of
   `12_80`. OpenMS parses the *entire* run — all 13,009 spectra and 21.2 M peaks
   — out of mzML in 3.55 s, i.e. ~0.27 ms/spectrum. The reader is therefore
   about **1,000× slower per spectrum than mzML parsing**, which puts a single
   pass over `12_80` at roughly an hour and over `astral` (307,590 spectra) at
   about a day. This is not a tuning matter; it makes extraction impossible
   until it is fixed, and it is the single most important number for Phase 2.
   `ODIAInfo` therefore decodes peaks only under `-peaks`.

2. **Not every encoding decodes.** Spectra 0–66 of `12_80` decode; spectrum 67
   (MS2) throws `chunked array decoding is not implemented (MS:1000515)` —
   the intensity array. 25 of the first 1,000 spectra are affected.

Both are WIP-branch issues rather than design faults, but Phase 2 cannot begin
extraction until at least the first is resolved.

## 1a. Upstream OpenMS/mzpeak cannot read our files at all

Worse than an API gap: `Index::spectra()` throws on every example archive.

```
MzPeak::ParquetError: metadata file missing or does not have the spectrum group
```

mzpeak's own `examples/read_spectra` fails identically, so this is not an ODIA
bug. `Metadata::Spectrum` requires a Parquet schema *group* named `spectrum`
inside `spectra_metadata.parquet`. Our archives have no such group — the
per-spectrum fields sit flat at the top level:

```
index, id, ms_level, time, scan_polarity, spectrum_representation,
spectrum_type, lowest_observed_mz, highest_observed_mz,
number_of_data_points, number_of_peaks, base_peak_mz,
base_peak_intensity, total_ion_current, data_processing_id,
+ groups: parameters, auxiliary_arrays, mz_delta_model
```

The writer confirms this is a deliberate format change. `mzpeak-convert 0.7.0`
contains the string:

> this mzPeak archive uses the **pre-0.7.0 packed metadata layout**, which this
> build cannot read; reconvert it from the source file with the current converter

So *packed* = the nested `spectrum` group; 0.7.0 flattened it. The example data
was written by 0.7.0 and is current; **the OpenMS C++ library (trunk,
2026-07-31) still implements the pre-0.7.0 layout and has not caught up.**

Consequences:

- The mzpeak high-level API is unusable for this data today, not merely
  insufficient. ODIA reads `spectra_metadata.parquet` through
  `Index::parquet()` → `Util::Parquet::reader()` → Arrow instead. `ODIAInfo`
  already does this and agrees with the mzML path on `12_80` (13 009 spectra,
  1 083 MS1 + 11 926 MS2).
- Arrow/Parquet is therefore a **direct** ODIA dependency, not a transitive one.
- Note the units trap: mzPeak's `time` column is in **minutes**, mzML's
  retention time in **seconds**.
- Worth raising with the mzpeak maintainers; it likely blocks anyone else
  trying to use the C++ library against current files.

## 1b. OpenMS has no mzPeak file type

`FileTypes` in OpenMS 3.6 has no mzPeak entry (it has `MZML`, `SQMASS`,
`CHROMPARQUET`, `PARQUET`, … but nothing for `.mzpeak`). A TOPP tool therefore
cannot declare mzPeak as a valid input format: calling `setValidFormats_` makes
TOPPBase sniff the file, see a multi-entry zip archive and reject it as a
malformed compressed mzML *before* `main_` runs. ODIA dispatches on the file
extension instead. Adding a `FileTypes` entry upstream is a small, obvious
contribution.

## 1c. Measured: m/z pruning does not work on the example archives

The project thesis is that mzPeak's row-group and page-index pruning makes
DIA-NN-style selective, repeated reads of a run affordable, removing the
constraint that pushed OpenSWATH to a single-tier design. That is testable
today, and the answer is currently **no for the m/z axis**.

All three example archives store their signal in `spectra_peaks.parquet` in the
flat **point** layout — `(spectrum_index u64, mz f64, intensity f32)`, one row
per peak. `spectra_data.parquet`, which would hold the *chunked* m/z layout
(`chunk.mz_chunk_start` / `chunk.mz_chunk_end`), is present in the schema but
has **zero rows** in all three.

Row-group statistics for `12_80` (21.17 M peaks, 21 row groups of 1,048,576):

| axis | clustered? | consequence |
|---|---|---|
| `point.spectrum_index` | **yes** — rg0 = spectra 0–897, rg1 = 897–2090, … | RT / spectrum-range pruning works |
| `point.mz` | **no** — every row group spans 100.0–1200.0 | m/z-band pruning is impossible |

Measured directly: a 10 ppm band touches **21 of 21** row groups at m/z 400 and
at m/z 700 (5 of 21 at m/z 1200, only because a few row groups top out just
below it). This is inherent to the layout — every spectrum covers the full m/z
range, so every row group does too.

**Re-converting with the chunked layout does not currently help.** Running
`mzpeak-convert --layout chunked` on `12_80.mzML` produced a file with the same
21,172,704 point rows and no chunk rows — the same layout, within 168 bytes of
the original — while emitting its own diagnostics:

```
[ERROR mzpeak_prototyping::chunk_series] BUG: signal array IntensityArray is
being spilled to auxiliary_arrays (metadata facet); signal arrays must live in
spectra_data/spectra_peaks
```

So the chunked path appears to fail and fall back to point. Cause unknown;
needs the converter's authors. This also explains the `S08_diaPASEF` size
anomaly — 13 GB of mzPeak against 1.3 GB of mzML is what the point layout costs
at 20 bytes per peak with no chunk-level encoding.

**What survives**, and it is not nothing: spectrum-index pruning is real and DIA
extraction does restrict to an RT window; column projection works; and Parquet
decode avoids XML parsing entirely (`ODIAInfo` reads the full spectrum inventory
of `12_80` in 0.07 s against 3.58 s for OpenMS to parse the mzML).

**What does not:** reading only the m/z bands a set of transitions needs, which
is the specific capability the two-tier extraction design would trade on.

This has to be settled before Phase 2 commits to an extraction strategy.

## 2. Everything must be built with one toolchain

- mzpeak's headers are C++23 and use `std::print`, `std::ranges::to` and
  `std::move_only_function`. The system compiler (GCC 13) cannot compile them;
  GCC 14 is the minimum.
- OpenMS 3.6 newly requires **Arrow/Parquet ≥ 23** (`cmake_findExternalLibs.cmake`
  hard-fails below that) — the same requirement mzpeak has. This is convenient:
  one Arrow satisfies both, and OpenMS 3.6 is the first release where linking
  ODIA against both libraries does not mean two Arrows in one process.
- Anything ODIA includes from mzpeak drags C++23 into that translation unit,
  so ODIA is a C++23 project throughout.

## 3. Dependency versions unavailable from the distribution

| Component | Needed by | Ubuntu 24.04 | Resolution |
|---|---|---|---|
| GCC ≥ 14 | mzpeak (C++23) | 13 | conda-forge 14.4.0 |
| Arrow/Parquet ≥ 23 | OpenMS 3.6, mzpeak | — | conda-forge 23.0.1 |
| Boost ≥ 1.89 | mzpeak | 1.83 | conda-forge 1.89.0 |
| libzip ≥ 1.11.4 | mzpeak | 1.10.1 | **built from source** (conda-forge stops at 1.11.2) |

No root on these machines, so nothing can come from `apt`.

## 4. Two upstream mzpeak bugs are in the way

Both are carried locally; both should go upstream.

1. **`EnumerableProxy::Iterator` cannot compile under GCC.** It declares
   `Iterator(const Iterator&&) = default` and `operator=(const Iterator&&) =
   default`. A `const` rvalue reference is not a move constructor/assignment, so
   defaulting it is ill-formed. mzpeak's `flake.nix` pins `clang_22`, which is
   why upstream has not hit it. Patched in
   `patches/mzpeak-0001-defaulted-move-ctor.patch`.

2. **The installed header tree does not compile.** `meson.build` calls
   `install_headers(..., subdir: 'mzpeak')`, which *flattens* every header into
   one directory, while the headers include each other by nested path
   (`mzpeak/io/archive.h`, `mzpeak/util/parquet.h`). The umbrella header also
   lands at `include/mzpeak/mzpeak.h` rather than `include/mzpeak.h`, so even
   mzpeak's own README example (`#include <mzpeak.h>`) does not resolve against
   an installed copy. `scripts/build_mzpeak.sh` installs a correct tree instead;
   the upstream fix is `preserve_path: true`.

## 5. OpenMS 3.6 build notes

- `WITH_THERMO_RAW` defaults ON and fetches `openms-thermo-bridge`, which
  requires a `dotnet` executable. Not present and not needed — ODIA reads mzPeak
  and mzML — so the build sets `WITH_THERMO_RAW=OFF`.
- `OPENMS_USE_VCPKG` defaults OFF, so pointing `CMAKE_PREFIX_PATH` at the
  conda-forge prefix is a supported configuration and avoids a multi-hour vcpkg
  build.
- conda-forge splits libxml2 headers into `libxml2-devel`; without it OpenMS
  fails at `find_package(LibXml2)`.

## 5a. Enabling ONNX in OpenMS: two build traps

`WITH_ONNX=ON` gives OpenMS's PeptDeep bindings *and* makes the build download
the three PeptDeep ONNX models from `archive.openms.de` against pinned SHA256s.
ODIA therefore vendors no model weights. Getting it to configure took two
workarounds, both in `build_openms.sh`:

1. **Upstream bug: the Find module is not on the module path.**
   `cmake/FindONNXRuntime.cmake` lives in `cmake/`, but OpenMS only appends
   `cmake/Modules` and `cmake/Windows` to `CMAKE_MODULE_PATH`. So
   `find_package(ONNXRuntime REQUIRED)` never sees it, falls through to config
   mode, and fails looking for `ONNXRuntimeConfig.cmake` — which conda-forge does
   not provide under that name. Passing `-DCMAKE_MODULE_PATH=<src>/cmake` fixes
   it, because OpenMS's own `list(APPEND ...)` preserves a command-line value.
   **`WITH_ONNX=ON` cannot configure as shipped; worth reporting upstream.**

2. **The build rediscovers its own previous install.** `find_package` searches
   `CMAKE_INSTALL_PREFIX` implicitly, so configuring picks up the opentims
   package exported by the last OpenMS install, whose exported target references
   a `zstd::libzstd_shared` target that does not exist in that scope. Restricting
   `CMAKE_PREFIX_PATH` is not sufficient; `-DCMAKE_FIND_USE_INSTALL_PREFIX=OFF`
   is. Non-destructive — the previous install stays usable if the build fails.

## 5b. The environment is pinned, and must stay that way

The conda prefix is not just a build environment; the built artefacts link
against its exact sonames. Installing `pyarrow` into it once resolved Arrow
23 → 25 and Boost 1.89 → 1.91, after which `libmzpeak.so` and the ODIA tools
could no longer find `libarrow.so.2300`, `libparquet.so.2300` and
`libboost_json.so.1.89.0`, and stopped running. Two rules follow:

- `bootstrap_env.sh` writes `conda-meta/pinned`. conda and micromamba refuse to
  move a pinned spec, so an incidental install can no longer drag the ABI
  forward. Verified: `pyarrow` now resolves to 23.0.1 rather than forcing 25.
- **Never mutate the environment while a build against it is running.** The same
  transaction unlinked Boost mid-compile and produced a burst of
  `boost/assert.hpp: No such file or directory` errors that look like a missing
  dependency rather than a moving one.

## 6. Benchmarking discipline

`data/` is on ceph, a shared network filesystem — reading from it measures the
network, not the code. Benchmarks must read from node-local NVMe
(`/scratch`), staged per node with `scripts/stage_data.sh`.

This node (`ibminode05`) has 128 cores, ~1 TB RAM and 12 TB free on `/scratch`.
Slurm is installed but not configured here (`sinfo` cannot reach a controller),
so benchmarks run directly on the node for now.

## 7. Example data

| Run | Instrument | mzML | mzPeak | Spectra |
|---|---|---|---|---|
| `12_80` | Orbitrap DIA | 244 MB | 130 MB | 13 009 |
| `astral` | Thermo Orbitrap Astral | 6.0 GB | 3.2 GB | 307 590 |
| `S08_diaPASEF` | Bruker timsTOF diaPASEF | 1.3 GB | **13 GB** | — |

The diaPASEF pair inverts the usual ratio: the mzPeak archive is ten times the
mzML. Worth explaining before any storage or speed claim is made — the mzML is
presumably reduced while the archive retains the full ion-mobility data, and the
archive additionally embeds the Bruker `.m` method directory (75 of its 85
members are vendor files). The other two runs compress to roughly half of mzML,
as expected.
