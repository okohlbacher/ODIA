# ODIA performance: where the time and the memory actually go

Measured 2026-08-05 on `ibminode05`. Everything below is measured on this node
unless a line says otherwise; nothing is extrapolated except where it says
"derived" and shows the arithmetic.

---

## 0. The node, and why these numbers are trustworthy anyway

`ibminode05`: 128 cores, 995 GiB RAM. **Load average 6,374–6,606 for the whole
session** — roughly 50x oversubscribed, from ~90 `phold` processes belonging to
one other user. `top` showed 99.7% user CPU and ~6,900 runnable tasks
throughout. The GPU nodes were unreachable (see §6), so there was no idle
machine to move to.

A wall-clock number taken naively here would measure that other user's job. It
does not, because `kernel.sched_autogroup_enabled = 1`: each login session is a
scheduling group, so my session competes as *one* entity, not as one task among
6,900. I calibrated exactly how much CPU that buys with a pure busy-loop:

| threads | wall (s) | CPU (s) | effective cores | efficiency |
|---:|---:|---:|---:|---:|
| 1 | 1.171 | 1.169 | 1.00 | 100% |
| 2 | 1.176 | 2.345 | 1.99 | 100% |
| 4 | 1.175 | 4.677 | 3.98 | 99% |
| 8 | 1.180 | 9.372 | 7.94 | 99% |
| 16 | 1.257 | 18.740 | 14.90 | 93% |
| 32 | 1.368 | 37.422 | 27.35 | 85% |
| 64 | 1.966 | 74.882 | 38.10 | 60% |
| 128 | 2.915 | 149.813 | 51.40 | 40% |

**Up to 8 threads this node is effectively idle for my purposes.** Phase 2 is
single-threaded in the part that matters (measured at 108–113% CPU), so its
numbers are clean. Phase 1 was run at `-threads 16`, i.e. 93% efficiency — its
wall time is therefore at most ~7% pessimistic. The 64/128-thread rows are why
no measurement below uses more than 16 threads.

Two further confounds, stated up front:

* For the first Phase-2 run only (the 500-precursor row), **two other
  `OpenDIAlyzer` processes belonging to the same account** were running from
  another session, reading the same file. Later runs were clean of those.
* The run file `IH1_diaPASEF.mzpeak` was staged from Ceph to node-local NVMe
  (`/scratch`) before the main series, per `scripts/stage_data.sh`. This turned
  out **not to matter**: decode took 602.2 s from Ceph and 593.8–599.8 s from
  NVMe. Decode is CPU-bound, not I/O-bound. That is itself a finding.

Inputs: `human-bench.fasta` (20,416 proteins); `IH1_diaPASEF.mzpeak`
(12.75 GiB, 17,448 physical spectra, 32,210 MS2 window-entries, 24 isolation
windows); libraries as noted per run.

---

## 1. Phase 1 — library generation

`OpenDIAlyzer -fasta human-bench.fasta -out_lib <tsv> -stop_after library
-threads 16` with all three PeptDeep ONNX models. Output library:
**2,127,559 target precursors → 4,255,113 with decoys, 51,060,552 transitions**,
828.6 MiB in memory, 7.28 GB on disk.

### Stage table

Wall times for RT/MS2/CCS are ODIA's own log values. `load time` (609.29 s) is
ODIA's timer covering FASTA load through decoy append. The TSV write is *not
instrumented at all* and is obtained as `total wall − load time`. Per-stage CPU%
is **derived**, not directly measured: total CPU was 7,648.8 s over 821.1 s wall
(931%); the non-inference stages are serial by construction (confirmed by
reading the code), so all but ~243 CPU-s of that belongs to the three inference
stages.

| # | Stage | Wall (s) | % of wall | CPU% (derived) | Threaded? |
|---|---|---:|---:|---:|---|
| 1 | FASTA load, digest, dedup, modifications, precursor + fragment enumeration | ~31.6 | 3.8% | ~100% | no |
| 2 | RT prediction | 44.4 | 5.4% | ~1280% | 16 ONNX sessions |
| 3 | MS2 prediction + fragment re-rank | 458.6 | **55.9%** | ~1280% | 16 ONNX sessions |
| 4 | iRT calibration fit | (inside row 1) | — | ~100% | no |
| 5 | CCS prediction | 74.7 | 9.1% | ~1280% | 16 ONNX sessions |
| 6 | Decoy generation | (inside row 1) | — | ~100% | no |
| — | *(subtotal = ODIA's `load time`)* | *609.3* | *74.2%* | | |
| 7 | **`-out_lib` TSV write (7.28 GB)** | **211.7** | **25.8%** | ~100% | no |
| | **Total** | **821.1** | 100% | **931%** | |

**Peak RSS: 13.28 GiB** (13,926,912 KB). The RSS trace (sampled every 2 s) rises
from 0.12 GiB to a peak of 13.3 GiB at t≈510 s — i.e. during MS2 prediction —
then drops to a flat **2.72 GiB** for the entire TSV-write phase.

Rows 1, 4 and 6 could not be separated from each other without adding timers to
the source; together they are only 31.6 s (3.8%), so I did not spend a rebuild
on splitting them. **The digest is not a bottleneck**, which contradicts the
prior suspicion that it might be.

### Threads: `-threads 16` does not limit threads

The process ran with **129–144 OS threads** despite `-threads 16`. ONNX Runtime
and/or OpenMP size their pools from `hardware_concurrency()` (128), which
`-threads` never reaches. This is harmless here only because autogroup capped
the damage; on a busy node with per-process limits it is a liability.

### CPU vs GPU, and against the earlier measurement

The earlier figures in the task (53 min at 5,527% CPU; MS2 2,372 s, CCS 368 s,
RT 224 s; peak RSS 2.04 GiB) are **superseded**:

| | earlier | now (measured) | change |
|---|---:|---:|---:|
| MS2 | 2,372 s | 458.6 s | **5.2x faster** |
| CCS | 368 s | 74.7 s | 4.9x faster |
| RT | 224 s | 44.4 s | 5.0x faster |
| total wall | 53:23 | **13:41** | 3.9x faster |
| peak RSS | 2.04 GiB | **13.28 GiB** | **6.5x more** |

The speedup is almost exactly uniform (4.9–5.2x) across all three models, and
the memory went up 6.5x. That signature matches the session-replica change
documented in `LibraryGenerator.cpp`'s own table (1 session: 409.7 s / 0.86 GiB;
16 sessions: 34.6 s / 8.82 GiB). The old numbers were the pre-fix, low-session
configuration. **The 5x speedup was bought with memory, and the report should
stop quoting 2.04 GiB.**

**GPU could not be re-measured this session** (§6). The recorded prior result is
53:23 → 3:58 on one H100. Note what today's numbers do to that claim: with the
CPU inference now 5x faster, and with the TSV write at 211.7 s *unaffected by
any accelerator*, a GPU run today would be bounded below by ~212 s of serial
formatting plus ~32 s of digest. The GPU's advantage on the whole stage is now
much smaller than 13x, and **the TSV write would be the largest single item in a
GPU library build**.

### Hypothesis tested: does raising `inferenceSessions` past 16 help?

`LibraryGenerator.cpp` caps sessions at 16 with a comment asserting 16 is the
knee and that 32 costs 16.11 GiB for 6% more speed. Measured with the existing
`odia_ms2_sessions` tool (MS2 model, 20,000 peptides, CPU):

| sessions | wall (s) | peptides/s | speedup vs 1 |
|---:|---:|---:|---:|
| 1 | 48.27 | 414.4 | 1.00x |
| 2 | 24.75 | 808.1 | 1.95x |
| 8 | 6.61 | 3,024.5 | 7.30x |
| 32 | 3.62 | 5,520.5 | **13.32x** |

Output was bit-identical at every session count (0 of 2,234,488 values differ).

**The cap looks too low, and the memory claim behind it looks stale.** 8 → 32
is still worth **1.83x**, and peak RSS for the whole sweep — whose largest
configuration is 32 sessions — was **6.57 GiB**, not the 16.11 GiB the comment
predicts. I could not measure 16 directly because the tool's sweep is hardcoded
to {1, 2, 8, 32}; so I cannot state the 16 → 32 gain exactly, only that the
curve has clearly not flattened by 8 and that 32 sessions cost far less memory
than the comment assumes.

### Phase 1 bottlenecks, ranked

1. **MS2 inference — 458.6 s (56%).** *ONNX Runtime*, with an ODIA-side policy
   knob (the 16-session cap). Evidence: ODIA's own timer, plus the session
   sweep above.
2. **`-out_lib` TSV write — 211.7 s (26%).** *ODIA's code, fully actionable.*
   Evidence: 7.28 GB in 211.7 s = **34.4 MB/s**. A micro-benchmark on this node
   (§3.5) gives `ofstream <<` at 25.4 MB/s and raw `fwrite` at 1,042 MB/s. This
   is formatting cost, not disk.
3. **CCS inference — 74.7 s (9%).** *ONNX Runtime.* Same lever as (1).

Digest, enumeration, iRT fit and decoy generation together are 31.6 s (3.8%) —
**not** a bottleneck.

---

## 2. Phase 2 — extraction

`OpenDIAlyzer -tr <lib> -in IH1_diaPASEF.mzpeak -stop_after extract -out_chrom
<tsv> -threads 8`, run file on local NVMe, `-mass_calibration off` except where
noted. No iRT calibration was supplied, so RT windows span essentially the whole
gradient (1,342 of 1,342 cycles) — this is the worst case for point count and it
is stated explicitly because every number below scales with it.

### Stage table and scaling with library size

| precursors | transitions | points (allocated) | nonzero | chrom MiB | decode s | index s | match s | extract s | wall | CPU% | peak RSS |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 500 † | 6,000 | 8,359,680 | 3,798,582 | 32.1 | 602.24 | 0.006 | 12.25 | 614.5 | 11:58 | 108% | 14.19 GiB |
| 2,500 | 30,000 | 41,508,708 | 18,486,918 | 158.9 | 599.76 | 0.134 | 14.21 | 614.1 | 11:37 | 110% | 14.32 GiB |
| 9,522 | 114,256 | 153,363,800 | 68,059,302 | 586.9 | 593.83 | 0.507 | 18.96 | 613.3 | 15:14 | 112% | **14.74 GiB** |
| 9,522 ‡ | 114,256 | 153,363,800 | 96,108,130 | 586.9 | 590.95 | 0.506 | 21.18 | 612.6 | 15:21 | 113% | 14.76 GiB |
| 100,000 | 1,199,962 | 1,610,662,790 | 695,012,675 | 6,162.6 | 596.30 | 5.219 | 61.72 | 663.2 | — § | — | **20.26 GiB** |

† mass calibration **on** (adds 84.45 s, see below) and read from Ceph, not NVMe.
‡ `-precursor_im_window 0`, i.e. today's mobility filter disabled.
§ stopped deliberately during the chromatogram write, which would have produced a
~76 GB TSV at ~24 MB/s; extraction itself completed and is reported in full.

**A 200x range of library size (500 → 100,000 precursors) moves decode by 1.1%**
(602.2 → 596.3 s). Points per transition is 1,342.3 at every size — exactly the
cycle count — confirming the CSR is allocated, not data-driven.

The 100,000-precursor row is also the sharpest test of the RSS attribution
below: 20.26 GiB total − 6.02 GiB of chromatogram CSR = **14.24 GiB**, against
**14.15 GiB** of non-CSR memory measured independently at 9,522 precursors.
**Agreement to 0.7%**, across a 10x change in library size.

Two stages are missing from ODIA's own instrumentation and had to be obtained by
subtraction from total wall:

* **Mass calibration** (`-mass_calibration auto`, the **default**): **84.45 s**,
  decoding 3,840 spectra. It is reported only inside the calibration report's
  `collect_seconds` and is *not* part of `decode_seconds`. It is single-threaded
  and `-threads` never reaches it.
* **Chromatogram TSV write**: at 9,522 precursors, wall 914 s − extract 613 s =
  **~301 s** to write 7.16 GB. Untimed and unlogged.

### What scales and what does not

* **decode is flat at 596–602 s across a 200x range of library size.** It is a
  function of the *file*, not of the library.
* **index (0.006 → 5.2 s) and match (12.2 → 61.7 s) grow with the library but
  stay ≤10% of extraction even at 100,000 precursors** — 3.5% at the realistic
  9,522. The only threaded stage in Phase 2 is not the problem; that is why the
  whole run sits at 108–113% CPU no matter what `-threads` says.
* **`points` scales linearly with transitions** at 1,342 points per transition
  (= every cycle in the run). It is *allocated*, not data-dependent.
* **Peak RSS = a constant ~14.2 GiB of reader + 4 bytes per allocated point.**
  That is the whole memory model of Phase 2, and both terms are measured.

### RSS attribution — measured, not guessed

Peak RSS is **14.2–14.8 GiB at every library size from 500 to 9,522
precursors**, while ODIA's chromatogram arrays over that range grow from 32 MiB
to 587 MiB. **At least 96% of resident memory is not ODIA's output.** Only at
100,000 precursors does ODIA's own CSR (6.02 GiB) become a comparable term.
Four independent probes locate the constant part:

1. **RSS tracks the run file, not the library.** Same library, different runs:

   | run file | file size | spectra | peak RSS |
   |---|---:|---:|---:|
   | `12_80.mzpeak` | 0.12 GiB | 11,926 | 0.28 GiB |
   | `IH1_diaPASEF.mzpeak` | 12.75 GiB | 17,448 | ~14.5 GiB |

2. **A standalone probe reproduces it with no ODIA code at all.** `mzsplit`
   calls `MzPeak::Spectra::get_spectra_batch` in ODIA's exact 1024-spectrum
   block pattern and discards each batch immediately:

   | spectra decoded | peaks | peak RSS | RSS after loop |
   |---:|---:|---:|---:|
   | 1,024 | 21.0 M | 0.58 GiB | 0.64 GiB |
   | 2,048 | 44.5 M | 0.66 GiB | 0.59 GiB |
   | 4,096 | 105.1 M | 0.93 GiB | 0.82 GiB |
   | 8,192 | 998.6 M | 6.00 GiB | — |
   | 16,384 | 3,670.0 M | 8.70 GiB | — |

   **It grows monotonically and never plateaus at a block-sized working set**,
   and RSS *after* every batch has been destroyed is essentially equal to the
   peak. The memory is held by the long-lived `MzPeak::Spectra` object, not by
   the live batch and not by ODIA.

3. **It is not allocator slack.** `MALLOC_MMAP_THRESHOLD_=131072`
   `MALLOC_TRIM_THRESHOLD_=131072` reduced peak RSS only 6.04 → 5.36 GiB (−11%)
   and cost 16% more time; `MALLOC_ARENA_MAX=1` changed nothing. **~89% is
   genuine retention.**

4. **ODIA's own share is small and matches theory.** Decode-only peak RSS at
   4,096 spectra is 0.93 GiB; adding ODIA's exact copy loop takes it to
   1.45 GiB. So ODIA's 1024-spectrum `SpectrumPeaks` block costs **~0.52 GiB**,
   against a predicted 1,024 × 25,659 peaks × 16 B = 0.42 GiB.

**Attribution at 9,522 precursors (14.74 GiB peak):**

| holder | GiB | share | whose |
|---|---:|---:|---|
| mzPeak `Spectra` retained decode buffers | ~13.5 | **92%** | mzPeak reader |
| ODIA's 1024-spectrum peak block | ~0.52 | 3.5% | ODIA |
| ODIA chromatogram CSR (`intensity` + offsets) | 0.587 | 4.0% | ODIA |
| ODIA library + m/z index | ~0.03 | 0.2% | ODIA |

### How much of decode is the reader, and how much is ours

**Plain-Arrow floor.** `spectra_peaks.parquet` inside the mzPeak ZIP is 8.80 GiB
holding **3,688,828,005 points** in 3,518 ZSTD row groups (columns:
`spectrum_index` u64, `intensity` i32, `mean_inverse_reduced_ion_mobility` f64,
`tof` i32). Read end-to-end with pyarrow 23.0.1, **single-threaded** for a fair
comparison against ODIA's single-threaded decode:

> **63.86 s wall, 57.32 s CPU (90% — genuinely one thread), peak RSS 0.19 GiB,
> 57.8 M points/s.**

**mzPeak through its own API**, same file, same block pattern: **2.4–3.0 M
points/s, 8.2–8.6 ms per spectrum.**

> **The mzPeak reader is ~19–24x slower and uses ~50–76x more memory than a
> plain Arrow scan of the identical bytes.**

**ODIA's copy loop is free.** Interleaved, repeated, on the same 4,096 spectra:

| mode | ms/spectrum | peak RSS |
|---|---:|---:|
| decode only (touch, no copy) | 8.546, 8.583 | 0.93 GiB |
| decode + ODIA's exact copy | 8.580, 8.493 | 1.45 GiB |

The difference is inside the noise. **Essentially 100% of `decode_seconds`
belongs to the mzPeak reader, not to ODIA.**

One correction to the brief: the quoted "measured plain-Arrow floor of 0.758 s
for a whole file" **does not reproduce**. A full single-threaded Arrow scan of
the peak table is 63.9 s. 0.758 s is plausible only for the metadata tables
(`spectra_metadata.parquet` is 449 KB against 8.80 GiB of peaks). The honest
floor for the peak data is ~64 s single-threaded, so decode is ~9.3x over
floor — still large, but not the ~780x the 0.758 s figure would imply.

### The half of decode that *is* ODIA's fault

The file contains **17,448 physical spectra**; ODIA reports **32,210 MS2
spectrum entries**. The arithmetic: IH1 packs two isolation windows per frame,
so 16,105 MS2 frames × 2 windows = 32,210 entries, plus 1,343 MS1 frames =
17,448. `MzPeakSource::peaks` builds its request list as `info_[i].index` for
each entry, so **every MS2 frame is requested twice**.

mzPeak does no deduplication and no inter-request caching — asking for the same
index twice costs twice:

| pattern | requests | distinct frames | wall (s) | ms/request |
|---|---:|---:|---:|---:|
| unique indices | 4,096 | 4,096 | 34.82 | 8.500 |
| each index twice | 4,096 | 2,048 | 33.78 | 8.248 |

Cost is per *request*, flat, regardless of whether the frame was just decoded.

> **~50% of ODIA's ~594 s decode — roughly 297 s, a third of Phase-2 wall — is
> decoding frames it has already decoded.** That is ODIA's code, and it is the
> single largest directly-actionable item in this report.

### Effect of today's changes (per-precursor mobility window)

Same library, same run, mobility window on (0.025) vs off:

| | points allocated | nonzero points | match (s) |
|---|---:|---:|---:|
| `-precursor_im_window 0.025` (default) | 153,363,800 | 68,059,302 | 18.96 |
| `-precursor_im_window 0` | 153,363,800 | 96,108,130 | 21.18 |

**The mobility window rejects 29.2% of populated points and makes matching 10.5%
faster — and saves no memory at all.** The expectation in the brief that it
would "reduce the point count" holds for *populated* points only. The CSR is
preallocated as `transitions × cycles_in_rt_window` before any peak is seen, so
neither the mobility window nor the narrower m/z window can shrink it. Memory is
set by the library and the RT window, full stop.

### The 32-bit CSR ceiling binds far below a real library

`ChromatogramExtractor.cpp` throws above 2^32 points. Measured points per
precursor here: 153,363,800 / 9,522 = **16,106** (12 transitions × 1,342
cycles). So the ceiling is **266,664 precursors**. Verified by bisection:

* `-max_precursors 260000` → proceeds to extract.
* `-max_precursors 270000` → fails in **5.1 s** with
  `more than 2^32 chromatogram points; ... The CSR index is 32-bit.`

**Phase 1's own output library is 4,255,113 precursors — 16.0x over the
ceiling.** ODIA cannot today extract against the library it just generated.
This is not a corner case; it is the main line.

### Phase 2 bottlenecks, ranked

Against 914 s wall at 9,522 precursors:

1. **mzPeak decode — 593.8 s (65%).** Splits into two very different halves:
   * **~297 s: duplicate frame requests.** *ODIA's code — actionable.*
   * **~297 s: the reader itself, ~19x slower than plain Arrow.** *mzPeak,
     upstream.*
2. **Chromatogram TSV write — ~301 s (33%).** *ODIA's code — actionable.*
   7.16 GB at ~24 MB/s; see §3.5 for proof this is formatting, not disk.
3. **Mass calibration — 84.45 s** whenever `-mass_calibration auto` (the
   default) is used. *ODIA's code — actionable.* Single-threaded, decodes 3,840
   spectra, and ignores `-threads` (`MassCalibration::Options::threads` exists
   and is never read).

`index` (0.51 s) and `match` (19.0 s) together are 2.1% of wall. **The one part
of Phase 2 that is parallelised is the one part that does not need to be.**

---

## 3. Ownership summary

| # | Bottleneck | Cost | Whose | Directly actionable? |
|---|---|---:|---|---|
| P1-1 | MS2 inference | 458.6 s | ONNX Runtime (+ ODIA's 16-session cap) | partly — the cap is ours |
| P1-2 | `-out_lib` TSV write | 211.7 s | **ODIA** | **yes** |
| P1-3 | CCS inference | 74.7 s | ONNX Runtime | partly |
| P2-1a | Duplicate frame decode | ~297 s | **ODIA** | **yes** |
| P2-1b | mzPeak reader speed | ~297 s | mzPeak (upstream) | no |
| P2-2 | Chromatogram TSV write | ~301 s | **ODIA** | **yes** |
| P2-3 | Mass calibration | 84.45 s | **ODIA** | **yes** |
| P2-M | 13.5 GiB retained decode buffers | 92% of RSS | mzPeak (upstream) | no — but ODIA can bound it |
| P2-C | 2^32 CSR ceiling | hard failure >266k precursors | **ODIA** | **yes** |
| — | Node contention | up to 60% loss at 128 threads | the machine | no |

### 3.5 The formatting evidence, once, since three findings rest on it

Micro-benchmark on this node, `/scratch` NVMe:

| method | throughput |
|---|---:|
| `std::ofstream operator<<` (ODIA's pattern) | **25.4 MB/s** |
| raw `fwrite` of pre-rendered bytes | **1,042 MB/s** |

**41x.** The observed rates — 34.4 MB/s for the library TSV and ~24 MB/s for the
chromatogram TSV — sit right on the `ofstream` line and nowhere near the disk
line. Both writes are CPU-bound number formatting.

---

## 4. Costed suggestions

Ordered by measured payoff per unit of risk.

**A. Deduplicate frame requests in `MzPeakSource::peaks`.**
*Buys ~297 s — a third of Phase-2 wall.* Build `want` as the unique physical
indices in the block, decode once, and fan each decoded frame out to both window
entries that reference it. ~20 lines in `src/io/MzPeakSource.cpp`, plus an
entry→physical map. Low risk; output must be bit-identical, which is directly
testable against the current `chrom.tsv`. **Best ratio in this report.**

**B. Replace `ofstream <<` with `std::to_chars` in both TSV writers.**
*Buys ~276 s in Phase 2 and ~192 s in Phase 1.* Measured headroom is 41x; a
realistic `to_chars` + manual buffer implementation lands 10–15x, taking the
chromatogram write ~301 s → ~25 s and the library write 211.7 s → ~20 s. One
helper plus two call sites (`DIANNLibraryFile::storeTSV`,
`TOPPOpenDIAlyzer::writeChromatograms_`). Low risk, mechanical.
**Note this is the *largest* Phase-1 item on a GPU, where inference is ~4 min.**

**C. Make the CSR offsets 64-bit and chunk extraction by precursor blocks.**
*Buys the ability to run at all above 266,664 precursors.* `begin`/`count` are
per *transition*, not per point, so widening them costs 8 B/transition — nothing.
The real work is chunking, because 4.3 G points is already 16 GiB of `intensity`;
a full human library needs the extraction streamed in precursor blocks with
chromatograms flushed per block. This is the difference between a demo and a
tool.

**D. Parallelise mass calibration.** *Buys ~74 s of the 84.45 s default cost.*
`MassCalibration::Options::threads` already exists and is never read; the tool
never sets it. Collect is an embarrassingly parallel loop over sampled spectra.
Moderate effort, contained.

**E. Raise the `inferenceSessions` cap and expose it.** *Buys an estimated
~100 s on Phase-1 MS2, but measure 16-vs-32 first.* Measured: 8 → 32 sessions is
1.83x, and 32 sessions peaked at 6.57 GiB, not the 16.11 GiB the code comment
assumes. The comment's memory table should be re-measured before the cap is
defended again. Cost: one line, plus a CLI flag so it is not a hidden constant.

**F. Quantise chromatogram intensity to `uint8`** (already in `BACKLOG.md`).
*Buys 4x on the CSR point array and ~4x on chromatogram TSV volume*, and pushes
the 2^32 ceiling out only if combined with (C). Independent of everything above.

**G. Bound the reader's retained memory (upstream, or work around).**
mzPeak retains ~13.5 GiB. Until that is fixed upstream, ODIA could re-open the
reader every N blocks to cap RSS — but **I did not measure re-open cost**, so
this is a hypothesis, not a recommendation. The upstream report should carry the
`mzsplit` table from §2: the retention is reproducible in ~40 lines with no ODIA
code in the picture.

**H. Do not spend on:** the digest (31.6 s, 3.8% of Phase 1), the m/z index
(0.51 s), or the match loop (19.0 s, and already threaded). The match loop is
the only part anyone has parallelised and it is 2% of the run.

---

## 5. Interaction worth noting

Fixes A and B are additive and together take Phase 2 at 9,522 precursors from
**914 s to ~341 s** (decode 594 → 297, write 301 → 25). At that point **the
mzPeak reader alone is ~87% of what remains**, and every further gain has to come
from upstream or from not decoding what is not needed. Adding D (mass
calibration) matters only when it is enabled, which by default it is.

Phase 1 with fixes B and E lands near **530 s** (write 212 → 20, MS2 459 → ~358),
at which point MS2 inference is ~68% of the run and only the GPU moves it — and
on the GPU, fix B is what determines the floor, not the model.

The through-line: **after the two mechanical fixes, both phases are bounded by
things outside ODIA** — ONNX Runtime in Phase 1, the mzPeak reader in Phase 2.
Before them, both phases are bounded by ODIA doing avoidable work.

---

## 6. What I could not measure, and why

* **GPU (Phase 1 and any GPU comparison).** The Kerberos ticket at
  `krb5cc` expired 2026-08-05 08:54; `spock` and `data` authenticate with
  Kerberos, not SSH keys, and `kinit` needs a password that is not available to a
  non-interactive session. Every node probe returned `UNREACHABLE`. The
  `build-gpu` tree and `opt/env-gpu` were therefore never exercised. **The
  "3:58 on one H100" figure in the project memory is unverified today** and, per
  §1, is now bounded below by the 211.7 s TSV write regardless.
* **The 16-session inference point.** `odia_ms2_sessions` hardcodes its sweep to
  {1, 2, 8, 32}. I measured those four; 16 would need a source change and a
  rebuild. So suggestion (E) is quantified as "8 → 32 is 1.83x", not as a 16 → 32
  number.
* **A per-stage CPU% for Phase 1.** ODIA has no per-stage CPU accounting and
  `perf` is unusable on this node (`kernel.perf_event_paranoid = 4` blocks
  unprivileged sampling). The per-stage CPU% column in §1 is derived from total
  CPU plus the knowledge that the non-inference stages are serial; it is
  arithmetic, not measurement.
* **Splitting digest from decoy generation from fragment enumeration.** They
  share one untimed 31.6 s block. Separating them needs timers in
  `LibraryGenerator.cpp`. At 3.8% of Phase 1 it was not worth a rebuild.
* **Wall time and CPU% for the 100,000-precursor row.** Its *extraction* is
  reported in full (decode/index/match/points/RSS all measured), but the run was
  stopped during the chromatogram write, which would have emitted ~76 GB at
  ~24 MB/s and taken ~50 min for no additional information. So that row has no
  end-to-end wall or CPU% figure.
* **A direct in-situ measurement of chromatogram write throughput.** I intended
  to sample the growing `c100k.tsv`, but had already unlinked it during cleanup.
  The ~24 MB/s figure is therefore derived (7.16 GB / ~301 s at 9,522
  precursors) and corroborated by the §3.5 micro-benchmark, not sampled live.
* **Re-open cost for the mzPeak reader**, which suggestion (G) depends on.
* **An idle machine.** Load was 6,374–6,606 throughout. The autogroup
  calibration in §0 bounds the error at ≤7% for everything reported (all runs
  used ≤16 threads), but no number here was taken on a quiet node.

---

## Appendix — reproduction

Artifacts under `/scratch/kohlbach/odia/prof/`: `runner.sh` (wraps
`/usr/bin/time -v` and samples RSS every 2 s), `p2series.sh`, per-run
`*.log` / `*.time` / `*.rss` / `*.meta`, `arrowfloor2.py` (plain-Arrow floor),
`lib50k.tsv`, `lib200k.tsv`. Probes under the session scratchpad: `calmt.c`
(contention calibration), `mzsplit.cpp` (decode vs copy, RSS growth),
`mzdup.cpp` (duplicate-request cost), `fmt.cpp` (formatting vs raw write).
