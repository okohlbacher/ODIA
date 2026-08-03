# mzPeak C++ reader: streaming requirements for DIA extraction

**Audience.** Whoever works on the mzPeak C++ reader
(`okohlbacher/mzpeak-openms`, and eventually `OpenMS/mzpeak`). This states what
a targeted DIA extraction engine needs from the reader, what the reader does
today, and how to reproduce every number.

**Status update, 2026-08-03, against `f93f938` (trunk).** The fork now provides
the *entry points* this document asked for -- `get_spectra_batch`,
`indices_in_time_range`, `extract_ion_chromatogram`, per-peak ion mobility -- and
R2 is satisfied as an API. **R1 is not.** Peak decode is still
**284.4 ms/spectrum** through `get_spectra_batch` (2,000 spectra of `12_80` in
568.7 s, batches of 512), against 276.8 ms/spectrum for the per-spectrum path
measured a day earlier. The batch call is a loop over spectra, not a shared
decode.

The cost this *should* carry is now measured rather than asserted.
`spectra_peaks.parquet` holds 21,172,704 rows in 21 row groups of ~1,048,576,
and one row group covers **873 distinct spectra**. Decoding row group 0 whole,
with pyarrow, takes **0.07 s** -- so a reader that decodes each row group once
and serves every spectrum in it would cost **0.082 ms/spectrum**, and a full
pass over the run would take **1.5 s** (0.116 ms/spectrum). That is *better than
mzML parity* (0.27 ms/spectrum), and it is 2,450x faster than the reader
achieves today.

So the gap is entirely in how the reader traverses, not in the format or the
layout: each spectrum access decodes a whole row group and keeps ~0.1% of it.
`extract_ion_chromatogram` inherits the same cost -- five XICs over the full RT
range did not finish in 10 minutes.

**Original status.** Measured on 2026-08-02 against `587a4fb` (branch `writer_test`),
built with GCC 14.4, Arrow/Parquet 23.0.1, Boost 1.89, on a 128-core node with
the archives on node-local NVMe.

**Summary.** Metadata access is good — better than mzML by a wide margin. Peak
access is ~1000× slower per spectrum than parsing the same data out of mzML,
which is the only thing standing between ODIA and a working extractor.

---

## 1. What the reader is being asked to serve

DIA extraction is not a sequential scan. For each library precursor it needs the
intensity of ~6 predicted fragment m/z values, at a stated tolerance, in every
MS2 spectrum whose isolation window contains the precursor and whose retention
time lies in a window around the predicted elution time. The output is one
extracted-ion chromatogram per transition.

At the scale we care about — `astral`: 307,590 spectra, 150 isolation windows,
a 7.1 M-precursor library — that is on the order of 10⁸–10⁹ (m/z, tolerance,
spectrum) probes per run. Three access patterns follow.

### A1 — Metadata sweep (once per run)

Read MS level, retention time, isolation window and ion mobility for every
spectrum, to build the window/cycle index. No peak data.

**Status: good.** 13,009 spectra of `12_80` in 0.29 s and 95 MB, against 3.76 s
and 412 MB for OpenMS to parse the equivalent mzML. Nothing is needed here.

### A2 — Bulk peak access over a retention-time block

Given a contiguous span of spectra, obtain their peaks so that many transitions
can be extracted from one pass. This is the primitive ODIA's extractor is built
on, because retention time is the axis that prunes: row groups in
`spectra_peaks.parquet` are clustered by `spectrum_index` (row group 0 = spectra
0–897, row group 1 = 897–2090, …), so an RT window maps to a contiguous row-group
range.

**Status: unusable.** See §2.

### A3 — m/z-ranged access

Given an m/z interval and a spectrum range, obtain only the peaks inside both.
This is what would make a two-tier extractor (cheap first pass over all
candidates, expensive second pass only around surviving peaks) affordable, and
it is the specific capability that would let a tool built on mzPeak adopt the
iterative structure DIA-NN uses without converting to a private format first.

**Status: not available.** See §3.

---

## 2. R1 — Peak decode throughput (blocking)

### Measured

| | |
|---|---|
| mzPeak, per-spectrum `Spectrum::mz()` | **276.8 ms/spectrum** (1,000 spectra of `12_80`, 276.8 s total) |
| mzML, OpenMS full parse of the same run | 3.55 s for all 13,009 spectra ⇒ **~0.27 ms/spectrum** |
| ratio | **~1000×** |

Extrapolated: one pass over `12_80` ≈ 1 hour; over `astral` (307,590 spectra)
≈ 24 hours. For comparison, the whole DIA search these runs come from takes
minutes to tens of minutes.

The cost is **flat per spectrum**, not amortising. That is the diagnostic: a
lazy decode that amortised over a row group would show a large first call and
cheap successors. Flat cost is the signature of re-reading or re-scanning a
large slice per spectrum — 21.17 M point rows live in 21 row groups of 1,048,576,
so touching a whole row group per spectrum would produce roughly this.

### Requirement

**R1.** Sequential access to every spectrum's peaks must amortise to at most the
cost of decoding each Parquet row group once. Concretely: a full pass over
`12_80` (21.2 M peaks, 135 MB of Parquet) should complete in **seconds, not
hours** — mzML parity (3.5 s) is the floor, and columnar decode should beat it
comfortably.

**R2.** A **batch entry point** — decode peaks for a spectrum *range* in one
call, returning the whole block — so callers are not forced into a
per-spectrum API for what is inherently a columnar read. ODIA's extractor wants
exactly this: "give me the peaks for spectra `[i, j)`".

### Reproduce

```cpp
#include <mzpeak.h>
#include <chrono>
#include <cstdio>
int main(int argc, char** argv) {
  auto index = MzPeak::open(argv[1]);
  auto spectra = index.spectra();
  long n = 0, fail = 0;
  auto t0 = std::chrono::steady_clock::now();
  for (const auto& s : spectra) {
    try { (void)s.mz().size(); } catch (const std::exception&) { ++fail; }
    if (++n >= 1000) break;
  }
  auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
              std::chrono::steady_clock::now() - t0).count();
  std::printf("1000 spectra: %ld ms (%ld threw) -> %.2f ms/spectrum\n",
              ms, fail, ms / 1000.0);
}
```

Against `data/12_80.mzpeak`. The mzML baseline is
`ODIAInfo -in 12_80.mzML`, whose loader reports its own parse time.

---

## 3. R3 — m/z pruning requires the chunked layout to work

### Measured

All three of our archives store signal in `spectra_peaks.parquet` in the flat
**point** layout — `(spectrum_index u64, mz f64, intensity f32)`.
`spectra_data.parquet`, which carries the chunked layout
(`chunk.mz_chunk_start` / `chunk.mz_chunk_end`), has **zero rows** in all three.

Row-group statistics for `12_80`, 21 row groups:

| predicate | row groups touched |
|---|---|
| 10 ppm band at m/z 400 | **21 / 21** |
| 10 ppm band at m/z 700 | **21 / 21** |
| 10 ppm band at m/z 1200 | 5 / 21 |
| spectrum index range | prunes correctly |

Every spectrum spans the full m/z range, so every row group does too, and no m/z
predicate can prune. This is inherent to the point layout, not a statistics
problem.

Re-converting does not currently help. `mzpeak-convert 0.7.0` with
`--layout chunked` on `12_80.mzML` produced a file with the same 21,172,704
point rows and no chunk rows — the same layout within 168 bytes — while printing:

```
[ERROR mzpeak_prototyping::chunk_series] BUG: signal array IntensityArray is
being spilled to auxiliary_arrays (metadata facet); signal arrays must live in
spectra_data/spectra_peaks
```

This also explains a size anomaly worth flagging on its own: `S08_diaPASEF` is
**13 GB of mzPeak against 1.3 GB of mzML**. At 20 bytes per point with no
chunk-level encoding, the point layout is what costs that.

### Requirement

**R3.** The chunked layout must round-trip through the converter, and the reader
must decode it. Until then, m/z pruning is unavailable and any tool on mzPeak
must read all peaks in an RT range and filter in memory.

**R4.** Given the chunked layout, expose a **ranged query** —
`(spectrum range) × (m/z interval) → peaks` — that pushes the m/z predicate into
Parquet row-group and page-index selection. `Util::Query`, `Util::Planner`
(which already returns `{row_group, offset, length}` ranges) and
`Util::Projection` look like the right machinery; what is missing is data laid
out so the predicate can prune.

---

## 4. R5 — Not every encoding decodes

`12_80` spectra 0–66 decode. Spectrum 67 (MS2) throws:

```
chunked array decoding is not implemented (MS:1000515)
```

MS:1000515 is the intensity array. 25 of the first 1,000 spectra are affected.

**R5.** All encodings present in files produced by the current converter must
decode, or the reader must be able to say so *without throwing* — see R6.

---

## 5. R6 — No exceptions, and no `shared_ptr`, in the hot path

Two API-shape requirements that come from the extractor's inner loop rather than
from correctness.

**R6a — a non-throwing capability query.** ODIA must be able to ask "can this
spectrum's peaks be decoded?" and branch, without paying for a thrown exception
per spectrum. Exceptions are fine for genuine faults; they are not fine as
control flow at 10⁵–10⁶ calls per run.

**R6b — raw, borrowed access to decoded arrays.** `Slice::raw()` returns
`std::vector<std::shared_ptr<arrow::Array>>`. If an accessor of that shape is
called once per (row, column), the atomic refcounting on a control block shared
by every thread **scales negatively**: a comparable row loop in a sibling project
went **20.1 s serial → 34.5 s on 64 threads, burning 1,135 s of CPU**, and was
fixed only by returning `const arrow::Array*` valid for the owning table's
lifetime.

So the reader should offer a borrowed-pointer accessor alongside the owning one,
and document that the borrowed form is the one for loops. ODIA will hoist
resolution out of its inner loops regardless, but the API should not make the
negative-scaling form the path of least resistance.

**R7 — bounded memory.** Block-wise access must let a caller hold one RT block
at a time and release it, so peak memory is a function of block size rather than
run size. This matters at `astral` scale (3.3 GB of peaks) and much more at
`S08_diaPASEF` (9.4 GB).

---

## 6. Priority

| # | Requirement | Blocking? |
|---|---|---|
| R1 | Peak decode amortises to one row-group decode | **yes** — nothing works without it |
| R2 | Batch/range decode entry point | **yes** in practice; R1 is hard to reach without it |
| R5 | All present encodings decode | **yes** for correctness on current files |
| R6a | Non-throwing capability query | no, but cheap and avoids a nasty cost |
| R6b | Borrowed-pointer accessors | no, but shapes the API before callers depend on it |
| R3/R4 | Chunked layout + ranged m/z query | no — unlocks two-tier extraction, not basic function |
| R7 | Bounded memory over blocks | no at `12_80`; yes at `astral`/`S08` |

R1, R2 and R5 are the set that turns the reader from a metadata reader into
something an extractor can be built on.

---

## 7. What already works, and should not regress

- Split-metadata (v2) layout reading. Upstream `OpenMS/mzpeak` trunk throws
  `metadata file missing or does not have the spectrum group` on all three of
  our archives; this fork reads them.
- `retention_time()` (seconds), `precursors()` with isolation windows, and
  `ion_mobility()`. Verified against OpenMS on `12_80`: identical spectrum
  counts, MS-level split, RT range, isolation-window scheme and cycle count.
- Metadata sweep performance: 0.29 s and 95 MB for 13,009 spectra, against
  3.76 s and 412 MB for OpenMS on the mzML.
