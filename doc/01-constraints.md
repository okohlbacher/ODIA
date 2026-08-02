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

## 1a. The mzpeak C++ reader cannot read our mzPeak files at all

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
