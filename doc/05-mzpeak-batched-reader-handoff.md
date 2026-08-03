# mzPeak C++ reader: building a batched, streaming peak decoder

**Audience.** Whoever implements peak decoding in the mzPeak C++ reader
(`okohlbacher/mzpeak-openms`, eventually `OpenMS/mzpeak`).

**The ask, in one sentence.** Decode each Parquet row group **once**, then hand
out the spectra it contains one at a time, instead of decoding a row group per
spectrum.

**Why it matters.** Every peak-touching call is currently ~2,450x slower than
the file layout allows. Nothing about the format or the files is wrong — the
data is already grouped the way a fast reader wants it. This is a traversal
problem, and it is the only thing standing between mzPeak and being usable as
the input to a DIA search engine.

**Status this was measured against:** `f93f938` (trunk), built with GCC 14.4,
Arrow/Parquet 23.0.1, Boost 1.89, meson release build, archives on node-local
NVMe, 128-core node.

---

## 1. The measurement

### What it costs today

| path | cost |
|---|---|
| `Spectra::get_spectra_batch`, batches of 512 | **284.4 ms/spectrum** (2,000 spectra of `12_80` in 568.7 s) |
| per-spectrum `Spectrum::mz()`, measured 2026-08-02 at `587a4fb` | 276.8 ms/spectrum |
| `Spectra::extract_ion_chromatogram` | 5 XICs over the full RT range did not finish in 10 minutes |
| OpenMS parsing the same run from mzML | ~0.27 ms/spectrum (3.55 s for 13,009 spectra) |

The batch entry point is therefore **not** a batched decode. It is a loop over
per-spectrum decodes: 512 spectra in one call cost 512 x the single-spectrum
price. The API shape asked for in `03-mzpeak-streaming-requirements.md` (R2)
landed; the throughput requirement (R1) did not.

### What it should cost

`12_80.mzpeak` → `spectra_peaks.parquet`: 21,172,704 rows in **21 row groups**
of 1,048,576 rows, **12.9 MB each**. Row group 0 covers **873 distinct
spectra** and decodes whole, with pyarrow, in **0.07 s**.

| | |
|---|---|
| one row group, decoded once | 0.07 s for 873 spectra ⇒ **0.082 ms/spectrum** |
| whole run (21 groups) | **1.5 s** for 13,009 spectra ⇒ 0.116 ms/spectrum |
| measured today | 284.4 ms/spectrum |
| **ratio** | **~2,450x** |

A correct implementation is not merely "as fast as mzML" — at 1.5 s against
mzML's 3.5 s for the same run, it is **more than twice as fast**, which is the
whole point of a columnar format.

### Reproducing both numbers

```bash
# what it costs now (C++, against the installed fork)
g++ -std=c++23 -O2 bench.cpp -o bench -I$ODIA_MZPEAK/include -I$ODIA_ENV/include \
    -L$ODIA_MZPEAK/lib -L$ODIA_ENV/lib -lmzpeak -larrow -lparquet
./bench 12_80.mzpeak 2000
```
```cpp
// bench.cpp -- batches of 512 through the batch entry point
auto spectra = MzPeak::open(argv[1]).spectra();
for (std::size_t b = 0; b < limit; b += 512) {
  std::vector<std::size_t> want(std::min<std::size_t>(512, limit - b));
  std::iota(want.begin(), want.end(), b);
  for (auto& s : spectra.get_spectra_batch(want)) peaks += s.mz().size();
}
```
```python
# what it should cost (python, same file, decode one row group whole)
import zipfile, io, time, pyarrow.parquet as pq
pf = pq.ParquetFile(io.BytesIO(zipfile.ZipFile(f).read('spectra_peaks.parquet')))
t = time.time(); tbl = pf.read_row_group(0); print(time.time() - t)   # 0.07 s
idx = tbl.column(0).combine_chunks().field('spectrum_index').to_numpy()
print(len(set(idx)))                                                  # 873
```

---

## 2. What the layout already gives you

These are properties of the files as written, verified on `12_80`. They are what
make the fix straightforward rather than a redesign.

**`point.spectrum_index` is monotonically non-decreasing within every row
group.** Checked on all 21 groups of `12_80`. A forward cursor is therefore
enough — no sorting, no hash grouping. The boundaries between spectra are found
with a single linear pass, or with `std::upper_bound` per spectrum if you prefer
random access within the decoded block.

**Every row group carries min/max statistics on `point.spectrum_index`.**

```
group  0: 1,048,576 rows  12.9 MB  spectrum_index 0..897
group  1: 1,048,576 rows  12.9 MB  spectrum_index 897..2090
group  2: 1,048,576 rows  12.9 MB  spectrum_index 2090..2877
group 20:   201,184 rows   2.7 MB  spectrum_index 12836..13008
```

So random access to spectrum *i* needs no scan: pick the row groups whose
[min, max] contains *i*. That is one footer read, already in memory.

**Spectra straddle row-group boundaries — 20 of the 21 groups in `12_80` begin
with a spectrum the previous group started.** A reader that assumes a spectrum
lives entirely in one group will silently truncate the first spectrum of nearly
every group. This is the single most likely way to get a batched reader subtly
wrong, and it will not look like a crash: it looks like slightly low intensities
on ~1 spectrum in 900.

**A decoded row group is small.** 1,048,576 rows x (8 + 8 + 4) bytes ≈ 21 MB
materialised, 12.9 MB on disk. Holding one, or two during a boundary carry-over,
is nothing. There is no reason to stream *within* a row group.

---

## 3. Three regimes, and why one implementation must handle all of them

The ratio of spectra to row groups is not a constant — it swings by a factor of
40 across the three runs in this project. A design tuned for one will fail on
another.

| run | peak rows | row groups | rows/group | spectra per group | note |
|---|---:|---:|---:|---:|---|
| `12_80` | 21,172,704 | 21 | 1,048,576 | **~873** | Orbitrap DIA, point layout |
| `astral` | 512,278,842 | 489 | 1,048,576 | **~1,143** | same layout, 24x larger |
| `S08_diaPASEF` | 3,688,828,005 | 3,518 | 1,048,576 | **~27** | Bruker TDF *ims-compact* |

`S08` is the case that breaks a naive "one group holds many spectra" loop. Its
row group 0 spans `spectrum_index` 0..26. A diaPASEF frame is one spectrum
carrying millions of points across the mobility dimension, so **a single
spectrum can span several consecutive row groups** — the mirror image of the
`12_80` case, where a single group holds hundreds of spectra.

The reader must therefore treat the two as the same problem: *a spectrum is a
contiguous run of rows that may begin and end anywhere, including across group
boundaries.* Both regimes fall out of one cursor that carries a partial spectrum
forward.

`S08` also uses the ims-compact schema, which has no `mz` column at all:

```
point: struct<spectrum_index: uint64, intensity: int32,
              mean_inverse_reduced_ion_mobility: double, tof: int32>
```

m/z is reconstructed from `tof`. The batching is orthogonal to this — it is the
same row-group traversal over different columns — but a batched decoder must
project only the columns the caller asked for, or it will pay for `tof`
reconstruction on a metadata-only pass.

---

## 4. How ODIA actually reads a run

This is the access pattern the interface should serve. It is not hypothetical:
it is what `ODIA::SpectrumSource` (`include/odia/SpectrumSource.h`) already
declares, with an mzML implementation behind it today and mzPeak intended to
substitute.

**One forward pass, in acquisition order, is the dominant mode.** For each MS2
spectrum, ODIA needs `(m/z, intensity)` — and `ion_mobility` when the run has
it — and matches the library transitions whose isolation window contains their
precursor. There is no seeking backwards, no revisiting. A run is read once.

```cpp
// The shape ODIA wants. Peaks for a RANGE, decoded once, handed out in order.
virtual void peaks(std::size_t begin, std::size_t end,
                   std::vector<SpectrumPeaks>& out) = 0;
```

**Metadata is needed first, and separately, for the whole run.** Retention
time, isolation window bounds, mobility limits and MS level for every spectrum,
without decoding a single peak. mzPeak already does this well — better than
mzML — and that must not regress. ODIA plans the whole extraction from metadata
before it touches a peak.

**Peaks are wanted sorted by m/z within a spectrum**, because the inner loop is
a binary search for each transition's m/z at a ppm tolerance. If the decoder can
guarantee that ordering it saves ODIA a sort per spectrum; if not, say so in the
API so the caller sorts once rather than defensively every time.

**Random access exists but is rare.** Re-reading a specific spectrum happens for
diagnostics and for verifying one library entry against the raw data, not in the
extraction loop. It should be correct and it need not be fast; the row-group
statistics above make it cheap anyway.

**Ion mobility is per peak, and diaPASEF needs it.** A frame carries several
isolation windows over disjoint mobility ranges, so ODIA filters peaks by
mobility as well as m/z. `Spectrum::ion_mobility_array()` already provides this;
it must survive batching.

**What ODIA does not need**, so the implementation need not carry the cost:
chromatograms, MS1 peaks in the main loop (metadata only), any writing, and any
per-spectrum object identity — a spectrum can be a view into the decoded block
as long as its lifetime is documented.

---

## 5. Suggested shape

Not prescriptive about internals, but these are the semantics ODIA depends on.

```cpp
/// Decode row groups once; hand out spectra in order.
class SpectrumStream {
public:
  /// Spectra [begin, end), ascending. Decodes each underlying row group once.
  /// Successive calls with contiguous ranges must not re-decode the group that
  /// spans the join.
  std::vector<Spectrum> next(std::size_t begin, std::size_t end);

  /// Columns to materialise. Omitting m/z on a mobility-only pass, or
  /// intensity on a planning pass, should cost nothing to decode.
  void project(PeakColumns columns);
};
```

Three properties that matter more than the exact signature:

1. **A row group is decoded at most once per pass.** This is the whole fix. A
   batch spanning groups 3 and 4 decodes each once; the next batch starting
   inside group 4 must reuse it, not decode it again.
2. **A spectrum straddling a boundary is returned whole**, with its rows from
   both groups concatenated, or the caller is told explicitly that it is
   partial. Silently short spectra are the failure mode to design against.
3. **The returned peaks' lifetime is stated.** A view into the decoded block is
   fine and is faster; it just has to be documented, because ODIA holds several
   spectra at once while matching a window.

---

## 6. How to know it worked

**Throughput.** A full forward pass over `12_80.mzpeak` (13,009 spectra,
21.2 M peaks) in **under 3.5 s** — mzML parity — with **1.5 s** as the figure
the layout supports. Same pass over `astral.mzpeak` (512 M peaks, 489 groups)
should scale linearly: roughly 35 s, not 24 hours.

**Correctness, which the speed work can easily break.** Two checks catch the
realistic mistakes:

- **Total peak count per spectrum, against the current reader**, for every
  spectrum in `12_80`. The straddling bug shows up here immediately and nowhere
  else: it makes the first spectrum of 20 of 21 row groups short.
- **Value-for-value against mzML** for a sample of spectra across the run. The
  same three runs exist in `data/` as both `.mzML` and `.mzpeak`, so this is a
  direct comparison rather than a self-consistency check. mzPeak's own e2e suite
  already compares against the Rust reference; this adds an independent third
  implementation.

**A note on what "same" means.** The fork's README records that decoded values
match the Rust reference exactly at stored points and to <= 8.7e-07 Da at
null-reconstructed ones. A batched decoder must not widen that: it changes
*when* rows are decoded, never *what* they decode to. A tolerance that has to be
loosened to make the batched path pass is a bug in the batching.

---

## 7. Where the numbers in this document come from

Everything above is measured on this machine, on 2026-08-03, against `f93f938`:

- `284.4 ms/spectrum`: 2,000 spectra of `12_80` via `get_spectra_batch` in
  batches of 512, 568.7 s wall.
- `0.07 s` / `873 spectra` / `0.082 ms/spectrum`: `pq.ParquetFile.read_row_group(0)`
  on `spectra_peaks.parquet` extracted from `12_80.mzpeak`.
- Row-group tables, monotonicity and straddle counts: all 21 groups of `12_80`
  enumerated; `astral` and `S08_diaPASEF` from their Parquet footers.
- `0.27 ms/spectrum` for mzML: OpenMS full parse of `12_80.mzML`, 3.55 s for
  13,009 spectra, measured 2026-08-02.

The prior document, `03-mzpeak-streaming-requirements.md`, states the
requirements this one proposes an implementation for; its R1 is the requirement
still open.
