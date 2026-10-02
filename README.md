# ODIA — OpenDIAlyzer

DIA extraction tools built on [OpenMS](https://github.com/OpenMS/OpenMS) 3.6 and
the [mzPeak](https://github.com/OpenMS/mzpeak) file format.

ODIA is a **TOPP-compatible** tool set that lives **outside** the OpenMS source
tree. Its tools derive from `TOPPBase`, so they accept the usual TOPP command
line (`-ini`, `-write_ini`, `-threads`, …) and emit CTD, but OpenMS itself is
consumed as an installed, read-only dependency and is never modified.

> Status: infrastructure only. The scientific scope is still being defined —
> see `doc/`.

## Layout

| Path | Contents |
|---|---|
| `src/` | one `.cpp` per TOPP tool |
| `cmake/` | `FindMzPeak.cmake` |
| `scripts/` | dependency bootstrap, build, and data-staging scripts |
| `patches/` | local patches applied to dependencies while awaiting upstream fixes |
| `test/` | CTest definitions |
| `doc/` | design notes and project description |

## Dependencies

Both OpenMS 3.6 and mzpeak need software newer than Ubuntu 24.04 ships, so the
build provisions its own environment rather than using system packages:

- **Apache Arrow/Parquet ≥ 23** — required by OpenMS 3.6 *and* by mzpeak
- **Boost ≥ 1.89**, **libzip ≥ 1.11.4** — required by mzpeak
- **GCC ≥ 14** — mzpeak's headers use C++23 (`std::print`, `std::ranges::to`,
  `std::move_only_function`); the system GCC 13 cannot compile them

Everything comes from a single conda-forge prefix so that OpenMS, mzpeak and
ODIA share one C++ ABI. libzip is the exception: conda-forge stops at 1.11.2,
so it is built from source into the same prefix.

## Building

```bash
scripts/bootstrap_env.sh      # toolchain + libraries  (~2.5 GB, once)
scripts/bootstrap_libzip.sh   # libzip 1.11.4 from source
scripts/build_mzpeak.sh       # mzpeak, pinned commit
scripts/build_openms.sh       # OpenMS 3.6, pinned commit
scripts/build_odia.sh         # configure, build, test ODIA
```

Each script is idempotent and sources `scripts/env.sh`, which defines the
prefixes. Sources and build trees go to node-local `/scratch`; installed
dependencies go to `$ODIA_ROOT/opt` so they are shared across nodes.

To work in an already-provisioned checkout:

```bash
source scripts/env.sh
```

## Example data

The example runs are read-only at `$ODIA_ROOT/data`, each available as both
mzML and mzPeak:

| Run | Instrument | Spectra |
|---|---|---|
| `12_80` | Orbitrap DIA (small; used by the tests) | 13 009 |
| `astral` | Thermo Orbitrap Astral | 307 590 |
| `IH1_diaPASEF` | Bruker timsTOF diaPASEF (has ion mobility) | — |

Benchmarks must not read from `data/`: it is on ceph, a shared network
filesystem, so timings there measure the network. Stage to node-local NVMe
first, on every node used:

```bash
scripts/stage_data.sh          # smallest run
scripts/stage_data.sh --all    # everything (~24 GB)
```

## Local patches

`patches/` holds fixes carried against pinned dependency commits. Each should
be sent upstream and dropped from here once it lands.

- `mzpeak-0001-defaulted-move-ctor.patch` — `EnumerableProxy::Iterator` declares
  `Iterator(const Iterator&&) = default`. A const rvalue reference is not a move
  constructor, so it cannot be defaulted; GCC rejects it. mzpeak is developed
  against clang (its `flake.nix` pins `clang_22`), which does not.

`scripts/build_mzpeak.sh` additionally repairs the installed header layout:
mzpeak's `install_headers(subdir: 'mzpeak')` flattens the header tree, but the
headers include each other by nested path (`mzpeak/io/archive.h`), so the
installed tree does not compile as shipped.

## Pinned versions

| Component | Pin |
|---|---|
| OpenMS | `552802a` (develop, reports 3.6.0) |
| mzpeak | `80fff01` |
| GCC | 14.4.0 (conda-forge) |
| Arrow / Parquet | 23.0.1 |
| Boost | 1.89.0 |
| libzip | 1.11.4 |

mzpeak's README states its API is not yet stable, which is why it is pinned to
a commit rather than tracked.

## Licence

BSD-3-Clause, matching OpenMS.
