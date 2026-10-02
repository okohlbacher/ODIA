# Where ODIA's memory actually goes, and the order to fix it

Measured 2026-08-05 overnight, IH1_diaPASEF.mzpeak (12.75 GiB, 32,210 spectra)
unless stated. Every number here came from a run, not a model.

## The decomposition

Peak RSS is two terms and nothing else:

| term | size | scales with |
|---|---:|---|
| decode floor | **~10.3 GiB** | nothing measured so far |
| chromatograms | **3.5 GiB / 95k precursors** | library size x RT window |

The 4.26M-precursor library (`odia_v5.tsv`) was observed **live at 154.2 GiB
anonymous**. That is the answer to "why can't we use our own library": not a
projection, an observation.

Chromatogram arithmetic, from the extractor's own log at n100k:
`983,561,430 points (3.66 GiB), 504,467,125 nonzero` = 9,836 points/precursor
at `-rt_window 600`. At 4.26M precursors that is ~156 GiB.

## Five hypotheses closed by measurement

Do not reopen these from priors. Each cost a run; each was plausible.

| hypothesis | test | verdict |
|---|---|---|
| per-thread parquet row-group buffers | 1 / 4 / 16 threads | **refuted** — 10.34 / 10.27 / 10.35 GiB |
| page cache / mmap counted in RSS | `RssAnon` vs `RssFile` | **refuted** — RssFile 31 MiB |
| glibc arena fragmentation | tcmalloc `LD_PRELOAD` | **refuted** — 10.49 GiB, no better |
| ditto, arena proliferation | `MALLOC_ARENA_MAX=2` | **refuted** — 10.34 GiB |
| chromatograms are sparse | nonzero fraction | **refuted** — 51.3% dense, RLE loses |

The first is the explanation from the mzPeak handoff (224 workers -> 105.8 GB,
capped to 26.5 GB). Its precondition does not hold here: `MzPeakSource` holds one
`Index` and one `Spectra`, and the floor is flat in threads. A
`-mzpeak_decode_threads` equivalent would buy nothing.

The third and fourth matter because the *prior* OpenDIAlyzer measured 88% of RSS
as arena fragmentation and got 189 -> 106 GB from tcmalloc alone. That result does
not transfer. See [[odia-memory-is-the-problem]].

## What the floor still could be

Surviving shape, from a 3-second RSS sample (`slide/shape2.rss`):

* rises monotonically 3.86 -> 10.26 GiB over the first ~120 s,
* then **flat for the remaining 60% of the run**,
* on a 1,200-precursor library, single-threaded, where `decode` is 324 s of 332 s.

A fixed-depth cache over large parquet row groups explains all of: the plateau,
thread-invariance, library-invariance, and the magnitude. `spectra_peaks.parquet`
is a single **8.80 GiB** member inside the archive, and the floor is 81% of the
whole file.

**Resolved.** Same tiny-library run across three files:

| file | type | size | spectra | floor | floor/file |
|---|---|---:|---:|---:|---:|
| 12_80 | SCIEX SWATH | 0.13 GiB | 11,926 | 0.27 GiB | 2.08 |
| astral | Astral | 3.10 GiB | 303,701 | 0.91 GiB | 0.29 |
| IH1 | **diaPASEF** | 12.75 GiB | 32,210 | 10.27 GiB | 0.81 |

The floor is **not** fixed — it varies 11x with the input, so "a ~10.3 GiB
constant" was invariant only *within* IH1. It does not track file size, spectrum
count, or bytes-per-spectrum cleanly. Subtracting the ~0.25 GiB process baseline
(which 12_80 is essentially all of), the two **non-mobility** files sit at 0.15
and 0.21 of file size and the one **ion-mobility** file sits at 0.78 — ~4x worse.

**The floor is specific to the ion-mobility path.** That is where to look, and it
is code with a history of silent defects: band derivation, the co-packed
two-windows-per-frame handling, the scan-order-vs-1/K0 swap. The 154.2 GiB on the
4.26M library is therefore a diaPASEF tax, not a generic scale limit — Astral and
SWATH runs would not pay it.

Also uncovered: **astral took 31 minutes** for 303,701 spectra against IH1's 6 for
32,210. IH1-only benchmarking has been hiding a Phase 2 throughput problem.

Blocked on tooling: no `pyarrow`, `duckdb`, or `parquet-tools` on this node, so the
row-group count and size are not yet known. Get them via mzPeak's own metadata API
or a manual footer read.

## The floor, located (tcmalloc heap profile, 2026-08-06)

Six hypotheses were refuted by measurement before profiling: per-thread buffers,
page cache, glibc fragmentation, arena proliferation, sparse chromatograms, and
mzPeak's row-group cache (bounded — `cache_.erase(cache_.begin())` at
`kCachedGroups`, ~21 MB/group). Reading code produced six wrong answers; the
profile produced the right one in one run.

`HEAPPROFILE` + `libtcmalloc.so.4`, tiny library, IH1. Peak **live** heap 9.45 GiB
(dump 162 of 197; the final dump is 4.7 MB — profile the peak, not the end).
Cumulative: **248 GiB allocated over 23.5M allocations**.

Top five sites hold 8.76 GiB of the 9.45 GiB peak. Resolved with `addr2line`:

| bytes | site |
|---:|---|
| 2.79 GiB | `MzPeakSource::peaks` — `dst.mz.assign` (double, 8 B) |
| 1.39 GiB | `MzPeakSource::peaks` — `dst.intensity.assign` (float, 4 B) |
| 1.39 GiB | `MzPeakSource::peaks` — `dst.ion_mobility.assign` (float, 4 B) |
| 1.95 GiB | Arrow parquet decode batch |
| 1.24 GiB | Arrow parquet decode batch |

The 2.79 : 1.39 : 1.39 split is exactly the 8:4:4 byte ratio of the three arrays,
which is how the three were identified before symbols confirmed it.

**So 5.57 GiB of the floor is the `out` block of `SpectrumPeaks`** at
`MzPeakSource.cpp:258-262`, and ~3.2 GiB is Arrow's decode buffers.

Two independent fixes, both targeting that 5.57 GiB:

1. **Stop duplicating co-packed frames.** `info_` holds one entry per
   (spectrum, isolation window); a diaPASEF frame co-packs several windows into
   one physical spectrum, so IH1's 32,210 entries cover 17,448 physical spectra —
   **1.85x**. Both entries of a co-packed frame get *byte-identical* peak arrays
   (the comment at line 253 says so outright) and are separated only by their
   mobility bands. Sharing one buffer: 5.57 -> **3.01 GiB**. Lossless. This is why
   the floor is diaPASEF-specific — Astral and SWATH have no co-packing and no
   `ion_mobility` array, so they copy two arrays once instead of three twice.
2. **Reduce `BLOCK`** (`ChromatogramExtractor.cpp:693`, `constexpr 1024`, not
   tunable). Linear in block size; 256 would give another 4x. Tune, not fix — but
   note per-request cost is flat, so the throughput cost should be small.

Deferred: `mz` as `double` is 8 of the 16 B/peak. float32 quantisation is
0.03-0.06 ppm against a 10 ppm window (169-328x margin), so it is *probably* safe
— but the prior project saw an unexplained 22-peptide delta on a float32 path, so
this must be measured alone, never bundled with a representation change.

## The compaction order — this is the load-bearing part

Applying these in any other sequence measures nothing, because the third is
*negative* until the first lands.

**1. Narrow the RT window.** Points per precursor scale directly with it.
Measured at n100k: w600 -> w60 takes RSS 14.00 -> 10.75 GiB, a **3.98 GiB** saving.
This is the dominant lever and it is a prerequisite for step 3. It depends on the
iRT calibration being trustworthy enough to place a narrow window.

**2. uint8 chromatograms with a per-transition f32 scale.** A clean 4x on what
remains, ~156 -> 39 GiB at 4.26M. Already in the backlog. The scale is per
transition, so two orders of dynamic range fit in 8 bits without a sensitivity
cost.

**3. The extract/score fusion.** Only pays *after* step 1:

| window | old | new (fused) | delta |
|---|---:|---:|---:|
| w600 | 14.00 GiB | 14.85 GiB | **+0.85** |
| w60 | 10.75 GiB | 10.56 GiB | **-0.19** |

Its overhead scales with *live* precursors, and at w600 every precursor is live
(the extractor reports `5330 of 5330 precursors live at once`), so it pays a cost
proportional to the whole library and frees nothing. It also underperforms its own
theory at w60 — ~8.6% of precursors should be live, so it should free ~91% of the
chromatogram term; it freed ~56%. Worth understanding before trusting at scale.

Steps 1+2 together take the chromatogram term from ~156 GiB to roughly 10 GiB.
**At that point the floor is the dominant term**, which is why it is the right
thing to chase first.

## Correction worth carrying

I called the fusion "a regression" on the strength of the w600 arms. It is one at
the sizes tested — but those sizes cannot show its benefit, since its overhead is
proportional to live precursors. The verdict needs a narrow-window run at a size
where chromatograms dominate. Measuring it at w600 was measuring the wrong thing.
