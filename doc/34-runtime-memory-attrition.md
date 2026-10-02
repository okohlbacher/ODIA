# Where the 6 hours and the 1.5 TiB went

IH1 diaPASEF, our 4,991,901-target library (9,617,705 with decoys),
`-threads 96 -live_memory_gb 900 -rt_seed cirt -rt_window 60`, dax, exit 0.

**Read the whole document against one fact: the RT seed refused.** The tool
printed `Treat the output as a smoke test, not a result`, and it meant it. Every
cost below is the cost of extracting at approximately the wrong retention time.

## 0. The root cause, because everything else descends from it

`rt seed REFUSED at every contiguity threshold` -> no iRT map -> pass 1 ran with
`restrict_rt=0`, i.e. **the entire gradient extracted for every one of 9.6 M
precursors**. That single fact produces the runtime, the memory and the result:

| | pass 1 (no RT restriction) | pass 2 (+/-60 s) | ratio |
|---|---:|---:|---:|
| points extracted | 111,523,257,174 | 7,142,920,370 | **15.6x** |
| bytes per live precursor | 265,011 | 17,037 | **15.6x** |
| live set | 930.97 GiB | 66.09 GiB | **14.1x** |
| chunks over the file | 3 | 1 | 3x decode |
| wall | 14,842.6 s | 2,010.5 s | **7.4x** |

Pass 2 is the same work with a calibrated axis. It is not 15x cheaper in wall
time only because its serial fraction does not shrink with the point count.

## 1. Runtime

Wall 6:01:03 (21,663 s). CPU 23:12:42. **386% of one core on a 96-thread run
= 4.0% of the machine.**

| phase | wall | share |
|---|---:|---:|
| **pass 1 extraction** | 14,842.6 s | **68.5%** |
| MS1 trace build | 3,952.0 s | 18.2% |
| pass 2 extraction | 2,010.5 s | 9.3% |
| RT map fit + 2 refine rounds | 1,273.2 s | 5.9% |
| IM calibration (pass 2) | 240.8 s | 1.1% |
| library load | 52.2 s | 0.2% |
| mass calibration (x2) | 56.8 s | 0.3% |

Inside extraction:

| | pass 1 | pass 2 | threaded? |
|---|---:|---:|---|
| decode | 497.0 s | 161.3 s | worker pool |
| index | 15.6 s | 18.1 s | serial |
| match | 5,457.8 s | 1,322.2 s | **worker pool** |
| assemble | 1,033.4 s | 37.8 s | serial |
| **sink (scoring)** | **7,802.5 s** | 450.8 s | **SERIAL** |

**The sink is single-threaded and is 52.6% of pass 1.**
`ChromatogramExtractor.cpp:1218-1252` runs the `emit()` loop *after*
`pool_impl.run(work, threads)` returns, on the driver thread; `emit` calls
`sink.accept(trace)`, which is the on-the-fly scorer. Match parallelises,
scoring does not.

Amdahl on the measured split: sink + assemble + index + RT fit + MS1 build is
~14,550 s of 21,663 s. **Even with infinite threads this run cannot go below
~3x faster.** Adding cores is not the fix; parallelising the sink is.

## 2. Memory

Peak RSS **1,536,632 MiB = 1.465 TiB**.

| what | size | share | source |
|---|---:|---:|---|
| live chromatogram blocks (pass 1) | 930.97 GiB | 63.5% | log |
| MS1 traces | 48.12 GiB | 3.3% | log |
| assay library | 1.41 GiB | 0.1% | log |
| **unaccounted** | **~485 GiB** | **33.1%** | RSS minus the above |

Two defects here.

**`-live_memory_gb` under-counts by 63%.** The user asked for a 900 GiB budget;
the live set alone came to 930.97 GiB (already 3.4% over its own budget) and the
process used 1,465 GiB. The knob governs 5 planes per live precursor -- points,
ppm_num, ppm_den, im_num, im_den -- and nothing else. Sizing a job to a machine
with this flag will OOM by a factor of 1.6.

**The 485 GiB is measured, not explained.** The strongest candidate is BlockPool
retention: `blocks.give()` returns blocks to a size-keyed pool rather than to the
OS, so each size class keeps its own high-water mark and the pool never shrinks.
3,263,300,398 minor page faults -- ~12.4 TiB of page touches against a 1.465 TiB
footprint, ~8.5 per page -- is consistent with continuous fresh-page allocation
rather than reuse. This is unproven and worth a heap profile before anyone acts
on it.

## 3. Where the candidates went

Per-precursor attrition, targets only, pass 1:

| stage | targets | of library |
|---|---:|---:|
| in the library | 4,991,901 | 100% |
| all-zero chromatogram | -4,217,475 | -84.5% |
| fewer than 3 points | -640,557 | -12.8% |
| **reached the picker** | **128,584** | **2.58%** |
| scored peak groups | 372,224 (2.89/precursor) | |
| **identified at 1% FDR** | **2,436** | **0.049%** |

Pass 2: 44,694 targets reached the picker (fewer), **6,660 identified** (2.7x
more). Narrowing to +/-60 s removed wrong-RT candidates *and* revived the RT
feature -- pass 1 logged `dropping 2 sub-score(s) carrying no information:
var_rt_delta var_im_spread`. **The single most discriminating feature in any DIA
scorer was dead in pass 1 because there was no RT map to compute it against.**

Separately: `1,201,226 precursors (12.5%) are covered by no isolation window` --
the library spans a wider precursor m/z range than the diaPASEF method. That is
free to filter before extraction and nobody is doing it.

### The target/decoy ratios say the survivors are mostly noise

The library is 92.7% decoys by count, so a filter that cannot tell targets from
decoys prints 0.927. Measured:

| criterion | decoy:target | target enrichment |
|---|---:|---:|
| empty trace | 0.94x | 0.99x |
| no points (<3) | 0.87x | 1.06x |
| **precursors reached** | **0.80x** | **1.16x** |
| below min_corr | 0.79x | 1.17x |
| not a local max | 0.82x | 1.13x |
| outside max_corr_diff | **0.44x** | **2.11x** |

Reaching the picker enriches targets by 1.16x. **Of the 128,584 targets that got
that far, the overwhelming majority are indistinguishable from decoys** -- which
is what extracting at the wrong RT produces. `max_corr_diff` is the only
criterion doing real work (2.11x), and it fires on 3,619 precursors.

Of targets that produced a scored peak group, **1.9% were identified**. Against
DIA-NN's 37,334 on the same file and the same library, ODIA returns 6,660
(17.8%) -- from a run with no retention-time calibration.

## 4. Defects found while reading the log

1. **The pass-2 picker rejection table is cumulative across both passes**, while
   the peak-group and identification counts are per-pass. `empty trace 8,502,599`
   exceeds the 4,991,901 targets in the library. Anyone reading pass 2 in
   isolation gets nonsense.
2. **The measured fragment window is computed and discarded.** Pass 1 concluded
   `suggested half-width 3.53 ppm`; pass 2 extracted at +/-10 ppm. The string
   appears only in `MassWidth.h:155` -- it is printed, never applied. A 2.8x
   tightening is sitting unused.
3. **Pass 2's mass calibration is byte-identical to pass 1** (same 5,354/2,092
   residuals, same 2.15/3.28 peakedness, same -10.74 ppm) -- it re-derives the
   same answer from the same sample rather than using pass 1's 1,963 accepted
   identifications.
4. **`var_im_spread` is dropped as uninformative in every pass of both passes.**
   It has never carried information; `Options::observed_im` is not wired to
   `MobilityCalibration`'s per-precursor output.
5. **The RT-seed contiguity histogram is misleading since the CiRT subsetting
   fix** (mine): 4,991,614 precursors are reported at contiguity 0 when they were
   never measured at all. Not measured and measured-empty print the same.
6. **The file is decoded 3x in pass 1** because the live budget forces 3 chunks.
   Decode is only 497 s, so this is cheap now, but it scales with the budget.

## 5. What actually follows

Ranked by measured cost, not by appeal:

| | fix | recovers |
|---|---|---|
| 1 | **Make the RT seed succeed** (or accept `-irt_slope/-irt_intercept`) | pass 1's 15.6x point count, ~931 GiB, ~4 h, and `var_rt_delta` in pass 1 |
| 2 | **Parallelise the sink** | 52.6% of pass 1; the 4%-of-machine utilisation |
| 3 | Drop the 1,201,226 precursors no window covers | 12.5% of everything, free |
| 4 | Apply the measured 3.53 ppm window | 2.8x narrower extraction |
| 5 | Account for the 485 GiB, then make `-live_memory_gb` mean total RSS | correctness of the only memory knob |

Nothing in this list is a scoring change. **The scorer has not yet been measured
on a calibrated run** -- 6,660 is not evidence about the scorer.
