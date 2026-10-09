# Codex gpt-5.6-sol (effort max) review of the v0.5.0 step-back (2026-10-01, verbatim; reviewed from the pasted brief only — codex tool host broken)

1. Q1 — memory diagnosis

- [HIGH] The arithmetic decomposition is consistent, but “~165 GiB measured and unattributed” overstates what was measured. `657,612 MiB = 642.20 GiB`; subtracting the logged 425.944 GiB live payload and 50,835.2 MiB = 49.644 GiB MS1 matrix leaves 166.61 GiB. Subtracting the claimed 1.95 GiB library gives 164.66 GiB, but `c3_cent_gate_s43.log:7` was not pasted, so that library term cannot be verified (`c3_cent_gate_s43.log:46,51-52,232`). More importantly, the log reports `BlockPool::peakPoints()`, not `reservedPoints()`: live payload is exported at `src/extract/ChromatogramExtractor.cpp:1330-1334`, while retained allocation is separately tracked but never printed at `src/extract/ChromatogramExtractor.cpp:230-248`. Thus ~165 GiB is an arithmetic residual using maxima, not a measured owner attribution.

- [HIGH] The BlockPool does retain blocks across chunks, but not as a C++ object across extraction calls, `finish()`, or pass 2. It is local to one `extract()` invocation at `src/extract/ChromatogramExtractor.cpp:737-740`, owns every allocation until destruction at function exit, and the chunk loop is inside that lifetime at `src/extract/ChromatogramExtractor.cpp:964-965,1298-1337`. `extractInto_` returns before `sink.finish()` at `src/OpenDIAlyzer.cpp:4724-4735`; therefore the pool has already been destroyed before that finish. Freed pages may remain resident in the allocator, but that is not established by the pasted source. Diagnosis D’s “freed BlockPool arena not returned and finish vectors allocated on top” is explicitly a hypothesis, not measured source fact.

- [HIGH] F03’s denominator identity is only conditional:

  - Under `Aggregate::Sum`, finite positive intensities are added to `base` and `ppm_den` by identical `float += intensity` operations in the same per-spectrum worker and loop order, so those cells are bit-identical under that data precondition (`src/extract/ChromatogramExtractor.cpp:1187-1219,1225-1231`; race ownership at `:759-761`).
  - The source does not guarantee non-negative, finite intensities. A negative or NaN intensity is added to `base` but rejected by `intensity > 0.0f`, so `ppm_den` then differs (`src/extract/ChromatogramExtractor.cpp:1217-1231`).
  - `im_den` additionally requires finite `peak_im`; `hasIonMobility()` only establishes that an array exists, not that each element is finite (`src/extract/ChromatogramExtractor.cpp:1160-1162,1178-1179,1233-1241`). “Every peak carries 1/K0” must therefore mean every contributing value is non-NaN, and should be asserted.
  - Under `Aggregate::Max`, `base` keeps a maximum while the denominators sum positive contributions, so aliasing is invalid (`src/extract/ChromatogramExtractor.cpp:1198,1215-1231`).
  - A naive pointer alias is wrong: if `ppm_den` points at `base` but line 1231 is left intact, every positive intensity is added twice; the same applies at line 1241. The aliased path must suppress those writes and suppress duplicate `give()` calls (`src/extract/ChromatogramExtractor.cpp:945-958,1218-1241`).
  - The `n == 0` case is already multiply aliased through `dummy_`, but no cell can be addressed when there are no valid rows/cycles; `give()` is also a no-op for it (`src/extract/ChromatogramExtractor.cpp:215-220,240-244,868-881`).

- [MED] The “400 GiB budget” is not a hard byte bound; it is converted into a precursor-count cap using mean cells over all assignments (`src/extract/ChromatogramExtractor.cpp:676-699`). Each activated precursor nevertheless receives its actual `valid × cycles` cells in each enabled plane (`src/extract/ChromatogramExtractor.cpp:868-881`). The logged live set averages about 6,878 cells per precursor per plane, versus 6,459 used by the cap, explaining the 6.5% payload overshoot (`c3_cent_gate_s43.log:51-52`). Calling this “precursor-cap saturation” is exact; calling it literal byte-budget saturation is not.

- [MED] “All peak-live precursors were full-width” is not established. Besides composition differing from the global mean, activation is rounded early and expiry late by `MATCH_BATCH`; the code explicitly says the measured live peak differs from the interval sweep for this reason (`src/extract/ChromatogramExtractor.cpp:836-848,1116-1134,1254-1268,1330-1334`). In this run the logged peak precursor count equals the derived cap, but the source does not isolate which mechanism caused the byte excess (`c3_cent_gate_s43.log:51-52`).

- [MED] Several non-live-block allocations are real but omitted from the simple three-term equation. They include the `LiveSlot` array, chunk slot arrays, two current per-window slot orderings, and the m/z index’s `double + uint32 + uint16 + float` arrays per indexed transition (`src/extract/ChromatogramExtractor.cpp:93-115,675-725,737-740,1058-1069`). The omitted `:964-1057` index-building body means its capacity and lifetime cannot be fully verified. The scorer also retains 21.77 million groups until `finish()` (`c3_cent_gate_s43.log:226`; append at `src/score/PeakGroupScorer.cpp:2651-2655`), but `PeakGroup`’s definition is not pasted, so the claimed 435–437 B/group is not verifiable.

- [LOW] The decoded `SpectrumPeaks block` is not library-sized: it is bounded by `decode_block`, default 256, and reused across the extraction (`src/extract/ChromatogramExtractor.cpp:764,829-848,1099-1106`). The pasted historical measurement puts the entire 256-block small-library peak at 2.89 GiB, not hundreds of GiB (`include/odia/ChromatogramExtractor.h:675-686`). It belongs in the residual, but cannot plausibly explain ~165 GiB merely through the 9.9M-precursor library.

- [FINE] The dense-MS1 attribution is exact: allocation is `rows × bins × sizeof(float)` at `src/extract/Ms1Traces.cpp:63-82`, it is built once and retained in `ms1_traces_` for scoring at `src/OpenDIAlyzer.cpp:1840-1842,1883-1922`, and the logged 50,835.2 MiB matches 9,922,669 × 1,343 × 4 bytes (`c3_cent_gate_s43.log:46`).

2. Q2 — output identity

- [HIGH] F02 changes emission order even when it should not change chromatogram values. Different caps create different consecutive chunks; each chunk is completed before the next, while within a chunk expiry is sorted only by `hi` and emitted window-by-window (`src/extract/ChromatogramExtractor.cpp:701-729,964-965,1058-1069,1259-1268`). A long-lived precursor in an early chunk can therefore be emitted before a short-lived precursor that would have expired earlier in a larger combined chunk.

- [HIGH] Gate C is definitively order-dependent: it uses the first `n_needed` arriving decoys, admits everything until the sample becomes ready, and then applies the resulting quantile (`src/score/PeakGroupScorer.cpp:1244-1281`). `rt_first/rt_last` are only diagnostic TSV-independent effects, but `tau` and the identity of “admitted uncalibrated” precursors affect scoring. `c3_cent_gate_s43.log:61` shows zero observed Gate-C rejections in that arm; it does not prove that a changed chunk order will form the same null or remain inert.

- [HIGH] The final stable sort does not prove earlier scoring is order-independent. Groups are appended in emission order at `src/score/PeakGroupScorer.cpp:2651-2655`; the shown sort occurs later in `finish()` and exists partly to repair already-created anchor indices (`src/score/PeakGroupScorer.cpp:3143-3157`). Fold construction, tie handling and model fitting between append and that sort are not pasted. F02 therefore remains “expected identical, cmp required,” not established identical.

- [MED] Chunk edges are designed to preserve traces. The chunk spectrum interval is the min/max required RT range (`src/extract/ChromatogramExtractor.cpp:1073-1093`), early activation cannot add out-of-window values because every match rechecks `c ∈ [lo,hi)` (`src/extract/ChromatogramExtractor.cpp:1204-1213`), and delayed expiry happens only after matching the batch (`src/extract/ChromatogramExtractor.cpp:1254-1268`). Construction of `slot_rt_lo/hi` is not pasted, so exact boundary coverage still requires byte comparison.

Fix-by-fix classification:

- [HIGH] **F01 — not output-identical.** Narrowing/sampling pass 1 changes which chromatograms reach the pass-1 scorer and therefore the fitted calibration; `-passes 1` removes the refinement/re-extraction altogether (`src/OpenDIAlyzer.cpp:644-657,1032-1035,2378-2395`).

- [HIGH] **F02 — conditionally identical, not proven.** Trace arithmetic should be invariant, but chunk/emission order changes Gate C and may change unshown pre-sort consumers (`src/extract/ChromatogramExtractor.cpp:701-729,1259-1268`; `src/score/PeakGroupScorer.cpp:1248-1281,2655,3149-3157`).

- [HIGH] **F03 — conditionally identical.** It requires Sum mode, positive finite intensity, finite mobility for the IM alias, suppressed denominator writes, no double return, and the same cap/chunks. At unchanged 400 GiB, reducing five planes to three raises the derived cap and changes chunk order (`src/extract/ChromatogramExtractor.cpp:679-699,868-880,1218-1241`). Pinning `-max_live_precursors 3324888`, or using exactly 240 GiB for three planes, preserves the cap.

- [HIGH] **F04 — not currently safe to classify as identical.** Ordered commit can preserve group order, but Gate C must not run in speculative workers (`src/score/PeakGroupScorer.cpp:1248-1281`). There is also a missed lifetime bug: worker calls to `Session::add` instantiate `thread_local rejects_` and register raw pointers globally, but no destructor unregisters them (`src/score/PeakGroupScorer.cpp:1310-1329`). The extractor pool threads terminate before `sink.finish()` (`src/extract/ChromatogramExtractor.cpp:818-828`; `src/OpenDIAlyzer.cpp:4724-4735`), leaving dangling registry pointers unless F04 redesigns that reducer.

- [MED] **F05 — output-identical if implemented as disjoint-bin work.** Different MS1 frames write different `b+s` columns, and each cell uses Max (`src/extract/Ms1Traces.cpp:129-170`). The capped residual prefix must retain original frame/peak/target order if the log is required to match (`src/extract/Ms1Traces.cpp:119-126,156-178`). A new bucket index must retain the exact per-target ppm test at lines 140–144.

- [MED] **F06 — instrumentation alone is output-identical.** Printing `reservedPoints`, RSS and rusage does not alter scored data (`src/extract/ChromatogramExtractor.cpp:247-248,1330-1334`). Its conditional pool redesign and `malloc_trim` are separate changes and still require raw-TSV `cmp`.

- [MED] **F07 — deliberately behavior-changing.** It changes the calibration sample and therefore `tau` (`src/score/PeakGroupScorer.cpp:1248-1281`). Merely fixing the sample membership is insufficient: because the present code admits all precursors until the last calibration member arrives, a deterministic implementation must precompute tau, buffer early admissions, or guarantee the complete sample is evaluated first (`src/score/PeakGroupScorer.cpp:1251-1279`).

- [MED] **F08 — thread-count measurements should be identical; changing `MATCH_BATCH` is only conditionally identical.** Match workers own spectra/cycles, so scheduling does not change a cell’s accumulation order (`src/extract/ChromatogramExtractor.cpp:759-761,1139-1250`). A larger batch changes cross-window emission order and therefore reintroduces the F02/Gate-C issue (`src/extract/ChromatogramExtractor.cpp:836-848,1259-1268`).

- [LOW] **F09 — intended to be output-identical.** NUMA binding and allocator environment do not intentionally change arithmetic; any TSV difference exposes latent undefined/order-dependent behavior.

- [HIGH] **F10 — behavior-changing.** Library RT/IM coordinates influence extraction, mobility calibration and the borrowed MS1 mass centre (`src/OpenDIAlyzer.cpp:1786-1792,1885-1906`).

- [HIGH] **F11 — stage 0 produces no engine TSV; wiring is behavior-changing.** Enabling FRAGVEC changes feature data appended beside groups and ultimately the learned ranking (`src/OpenDIAlyzer.cpp:4431-4434`; `src/score/PeakGroupScorer.cpp:2651-2655`).

- [HIGH] **F12 — only a correctly mapped no-window row drop can be identical; spans are behavior-changing unless they cover both passes exactly.** With a keep mask, `Ms1Traces::at()` uses compact row indices, not library indices (`include/odia/Ms1Traces.h:90-103`), while the scorer currently passes library index `i` directly (`src/score/PeakGroupScorer.cpp:2423-2437`). The proposed row-drop implementation therefore needs a library→row map. No-window rows are safely empty (`src/extract/ChromatogramExtractor.cpp:535-540,1310-1321`); the 24,576 initially out-of-run rows could move after the RT refit, and their assignment logic was not pasted.

- [HIGH] **F13 and F14 — behavior-changing.** F13 changes the candidate set/reselection beyond the current maximum-three retained candidates; F14 changes the population used by the iterative model (`src/OpenDIAlyzer.cpp:1131-1135`; `include/odia/scoring/lda.h:1142-1167`; `src/score/PeakGroupScorer.cpp:2721-2727`).

- [LOW] The verification must compare the unmodified scored TSV. Externally sorting it before `cmp` weakens the stated byte-identity criterion and can hide an ordering regression.

3. Q3 — best first experiments

- [HIGH] **Memory: F02, `-live_memory_gb 150`.** It directly exercises the proven count-cap mechanism without code changes (`src/extract/ChromatogramExtractor.cpp:679-729`; `c3_cent_gate_s43.log:51-53`). Primary pass threshold: `Peak Memory Usage ≤ 420,000 MiB`, versus 657,612 MiB baseline (`c3_cent_gate_s43.log:232`). Mandatory correctness gate: raw scored TSV must be byte-identical; any difference is a chunk/order defect, not an acceptable memory trade. F03 should not precede this because an unchanged byte budget absorbs most plane savings by increasing the cap.

- [HIGH] **Wall clock: F01, first try `-rt_window_pass1 200`, not F04.** Pass-1’s shown components consume about 20,650 s, dominated by its 14,637 s serial sink, while pass 2 at ±110 s completes extraction in 6,504 s with only 66.3 GiB live payload (`c3_cent_gate_s43.log:49,185-189`). Primary performance threshold: total wall `≤ 7.0 h`. Hard scientific gate: median matched-entrapment \(N(e)\) over FDP 2–10% must be at least 99% of the same-night control; retained-anchor p95 alone is insufficient because narrowing preferentially removes the anchors whose errors exceed the new window. F04 is not the first experiment because it is weeks of code, has the rejects-registry lifetime defect, and its 20–30× sink scaling is unmeasured (`src/score/PeakGroupScorer.cpp:1310-1329`).

- [HIGH] **IDs: F11 stage-0 oracle first, as a falsification experiment.** Use the existing final-pass export path and require held-out recall at DIA-NN depth 37,334 of `≥72%` for 42+72 features before wiring anything (`src/OpenDIAlyzer.cpp:4431-4434`; `src/score/PeakGroupScorer.cpp:2651-2655`). Below 72%, do not fund the native scorer change; proceed to coordinate/extraction experiments such as F10. The cited 64–65% baseline and FRAGVEC gains are in documents not pasted here, so they cannot be independently verified; this is precisely why the offline gate should precede an engine arm.

4. Q4 — missing or under-specified work

- [HIGH] **Missing: a calibration-only pass-1 pipeline using a deterministic RT/charge/m/z-stratified target-plus-decoy panel and a lightweight anchor scorer.** The CLI itself says pass 1 exists to harvest anchors, yet the current path constructs the full `PeakGroupScorer::Sink`, scores every admitted precursor and runs full `finish()` (`src/OpenDIAlyzer.cpp:644-650,4718-4735`). Pure target-only extraction is invalid because pass-1 positive selection uses target/decoy group labels and q-values (`include/odia/scoring/lda.h:1144-1166`). A weighted representative decoy panel is viable; “drop every decoy” is not. Rank: after the zero-code F01 experiment, before F04.

- [HIGH] **Missing: enforce the live budget in points/bytes rather than mean-sized precursor counts.** The code already knows every assignment’s `valid × (hi-lo)` cost but averages it, then partitions only on precursor count (`src/extract/ChromatogramExtractor.cpp:679-699,710-729`). A weighted interval sweep/chunker including the batch halo would make 400 GiB a real bound and remove the observed 25.9 GiB overshoot. Rank: after F02 for magnitude, before treating the cap as production-safe.

- [HIGH] **F12 needs redesign, not merely a keep mask.** Build a stable library→compact-MS1-row map and share it with the scorer; otherwise compaction silently reads another precursor’s row (`include/odia/Ms1Traces.h:90-103`; `src/score/PeakGroupScorer.cpp:2423-2437`). Computing assignment/static no-window eligibility before allocating MS1 would safely remove about 6.38 GiB; dropping initially out-of-run rows before the fitted second pass is not proven safe (`src/extract/ChromatogramExtractor.cpp:535-540`; `src/OpenDIAlyzer.cpp:1032-1035`).

- [MED] **A sparse or span-based MS1 builder is preferable to naïve on-demand extraction.** Current build performs random row-major writes into `slot × bins + frame` for every match (`src/extract/Ms1Traces.cpp:82,129-166`), while scoring reads only candidate-local bins (`src/score/PeakGroupScorer.cpp:2423-2437`). On-demand per precursor would repeatedly search the same MS1 peaks; an inverted sparse arena or assignment-span matrix preserves the one-pass peak inversion.

- [MED] **Decode/match double buffering is real but low priority on IH1.** Decode precedes and blocks each pool run (`src/extract/ChromatogramExtractor.cpp:1099-1106,1138-1250`), so a dedicated decoder and two owned buffers could overlap block \(N+1\) with matching \(N\). The maximum IH1 saving is only the logged 592 s of decode, about 1.7% of wall (`c3_cent_gate_s43.log:49,186,232`), so it ranks below F01/F04/F05. Astral-specific benefit is not verifiable because the cited Astral log was not pasted.

- [MED] **F05 must address false sharing, not merely races.** Parallel frames write distinct cells but adjacent frame cells of the same precursor occupy the same row-major cache line (`src/extract/Ms1Traces.cpp:129-166`). Scheduling separated frame blocks, or partitioning by m/z/target buckets, is more credible than assigning adjacent frames to 64 workers.

- [LOW] **Float16/uint8 chromatograms are not a first move.** The old proposal covers only intensity cells and is behavior-changing (`doc/11-memory-and-compaction-plan.md:131-144`); current residual numerators are signed weighted sums and cannot share the same quantization scheme (`src/extract/ChromatogramExtractor.cpp:1225-1241`). Exact denominator aliasing should precede lossy storage.

- [LOW] **Per-block mmap plus `MADV_DONTNEED` is a poor first allocator fix.** `BlockPool` intentionally reuses blocks and zeroes them on every take (`src/extract/ChromatogramExtractor.cpp:200-245`); discarding every returned page would cause the same pages to be faulted and zeroed again without reducing the simultaneous live peak. A chunk-sized arena with phase-boundary discard is more defensible, after F06 identifies the peak phase.

- [LOW] **The mzPeak reader is not the first IH1 lever.** It already deduplicates physical spectrum requests inside each call (`src/io/MzPeakSource.cpp:372-402`), and total decode is only 592 s (`c3_cent_gate_s43.log:49,186`). Cross-chunk caching may help, but cannot explain either 9.5 h or 657 GB.

- [LOW] **`-threads 0` is a real CLI inconsistency but irrelevant to this run.** The extractor interprets zero as hardware concurrency, while both driver paths clamp it to one (`src/extract/ChromatogramExtractor.cpp:762-763`; `src/OpenDIAlyzer.cpp:1793,4480`). IH1 explicitly used 64 threads.

5. Q5 — which fixes bite

Astral-specific measurements cited by the diagnoses, including `e2e_astral.log`, were not pasted and cannot be verified; Astral judgments below are source-conditional.

- [HIGH] **F01 — BITES on IH1.** Pass 1 has 425.9 GiB live payload and 14,637 s sink versus pass 2’s 66.3 GiB and 4,203 s sink (`c3_cent_gate_s43.log:49-53,185-189`). ID preservation remains undecidable.

- [HIGH] **F02 — BITES on IH1.** The run is explicitly precursor-cap bound, so lowering the cap must lower concurrent live payload (`c3_cent_gate_s43.log:51-53`; `src/extract/ChromatogramExtractor.cpp:679-729`). Total-RSS response remains unmeasured because `reservedPoints()` was not logged.

- [HIGH] **F03 — BITES only with a pinned cap or proportionally lowered budget.** At unchanged 400 GiB, three planes raise the cap by about 5/3 and refill almost the same byte budget; the claimed −170 GiB does not occur automatically (`src/extract/ChromatogramExtractor.cpp:679-699`). On a no-mobility run, current “auto” still allocates IM planes and never fills them, so the bug is source-real (`src/OpenDIAlyzer.cpp:1404-1413,1778-1782`; `src/extract/ChromatogramExtractor.cpp:1178-1179,1236-1241`).

- [HIGH] **F04 — BITES as a bottleneck; claimed 20–30× is THEORETICAL.** The sink is demonstrably serial and consumes 18,841 s (`src/extract/ChromatogramExtractor.cpp:1250-1268`; `c3_cent_gate_s43.log:49,186`). The same pool’s match phase uses only a small fraction of 64-core capacity, so a 20–30× sink speedup has no pasted scaling support.

- [HIGH] **F05 — BITES as a bottleneck; 10–20× is THEORETICAL.** The entire 5,590 s MS1 loop is serial and uses a per-peak binary search (`src/extract/Ms1Traces.cpp:129-170`; `c3_cent_gate_s43.log:46`). False sharing and only 1,343 frame columns constrain simple frame parallelism.

- [MED] **F06 — UNDECIDABLE-FROM-BRIEF.** `reserved_` exists but is not reported, so neither BlockPool retention nor heap non-return has been measured at this scale (`src/extract/ChromatogramExtractor.cpp:230-248,1330-1334`).

- [MED] **F07 — THEORETICAL for current IH1 output.** The reference arm reports no Gate-C rejects (`c3_cent_gate_s43.log:61`); making tau positive could reduce work, but that is a behavior-changing admission experiment, not merely determinism.

- [MED] **F08 — UNDECIDABLE-FROM-BRIEF.** The pool uses 128-spectrum barriers and 64 workers (`src/extract/ChromatogramExtractor.cpp:805-848,1139-1250`), but no controlled scaling curve for this binary is pasted. Process packing will bite only after memory is reduced.

- [MED] **F09 — UNDECIDABLE-FROM-BRIEF.** The 11,150 s system CPU is real (`c3_cent_gate_s43.log:232`), but no pasted evidence assigns it to the driver, allocator trimming, NUMA hint faults or first-touch.

- [HIGH] **F10 — UNDECIDABLE-FROM-BRIEF.** The source confirms strong RT/IM/mass coupling (`src/OpenDIAlyzer.cpp:1786-1792,1885-1906`), but all claimed CCS dose-response and fine-tuned-library gains cite documents not pasted.

- [HIGH] **F11 — UNDECIDABLE-FROM-BRIEF.** The export/wiring points exist (`src/OpenDIAlyzer.cpp:4431-4434`; `src/score/PeakGroupScorer.cpp:2651-2655`), but the cited AUC and +232/+298 measurements are not pasted.

- [MED] **F12 — BITES for static no-window rows; span storage remains THEORETICAL.** The dense matrix is 50.8 GiB and no-window precursors are explicitly skipped by assignment (`c3_cent_gate_s43.log:46,55`; `src/extract/ChromatogramExtractor.cpp:535-540`). The proposed compact-row mapping is currently incomplete.

- [MED] **F13 — THEORETICAL.** The structural limit of three retained candidates and no score-changing RT refit is source-real (`src/OpenDIAlyzer.cpp:1131-1135`; `src/score/PeakGroupScorer.cpp:2721-2727`), but the A1/A2 measurements and expected ID gain were not pasted.

- [HIGH] **F14 — THEORETICAL.** The iterative learner exists (`include/odia/scoring/lda.h:1142-1167`), but the pasted lines do not show that all 22M decoy rows enter training, nor is the proposed tail-truncation experiment pasted.

6. Q6 — outright errors and overclaims

- [HIGH] **F14’s stated imbalance arithmetic is wrong.** “10.9M decoy groups against ~17k positives” is roughly 641:1, not 155:1. The reference run actually logs 10,897,786 decoy groups and 18,995 q≤0.01 IDs, still about 574:1 (`c3_cent_gate_s43.log:227-228`). A 155:1 ratio would require about 70,000 positives; such a training count was not pasted.

- [HIGH] **F03 is not implementable by merely pointing denominator pointers at `base`.** Without suppressing lines 1231 and 1241, it doubles the intensity chromatogram (`src/extract/ChromatogramExtractor.cpp:1218-1241`). The canonical description omits this necessary write-path change.

- [HIGH] **F04 incorrectly treats the existing thread-local reject registry as parallel-sink-ready.** Worker TLS objects register raw pointers with no unregistering destructor (`src/score/PeakGroupScorer.cpp:1310-1329`), while the local pool destroys those worker threads before `finish()` (`src/extract/ChromatogramExtractor.cpp:818-828`; `src/OpenDIAlyzer.cpp:4724-4735`).

- [HIGH] **F12’s “existing keep mask makes row drop risk-free” is incomplete.** Compact MS1 rows are explicitly addressed by compact row number, whereas the scorer supplies library precursor indices (`include/odia/Ms1Traces.h:90-103`; `src/score/PeakGroupScorer.cpp:2423-2437`).

- [MED] **Diagnosis D’s unconditional “two planes are bit-identical copies” is too strong.** The guards make that false for negative/NaN intensity and missing/NaN mobility (`src/extract/ChromatogramExtractor.cpp:1218-1241`). It is an asserted-data invariant, not a source theorem.

- [MED] **The ~165 GiB remainder is not “MEASURED by subtraction” as an owner.** The logged live payload is `peakPoints()`, while retained allocation is separately available as `reservedPoints()` but unreported (`src/extract/ChromatogramExtractor.cpp:247-248,1330-1334`). Peak-phase coincidence was not measured.

- [MED] **Diagnosis A’s unique explanation of the 6.5% overshoot is unsupported.** Mean-cell inversion is proven, but batch-rounded liveness is another explicit source mechanism (`src/extract/ChromatogramExtractor.cpp:685-699,836-848,1330-1334`).

- [MED] **B/D label the 77–78% serial fraction “MEASURED,” but part is inferred.** Sink and MS1 alone are measured; decode/index/assemble are source-serial, but calibration, tail and `finish()` allocations are reconstructed from brackets not pasted (`c3_cent_gate_s43.log:46,49,186`; `src/OpenDIAlyzer.cpp:4724-4743`). The exact 77–78% and 1.28× ceiling should be labelled inferred.

- [MED] **F07 does not automatically make chunk/thread output deterministic.** A fixed calibration set still admits everything until its last member arrives; arrival order therefore remains observable unless tau is computed before admission or early calls are buffered (`src/score/PeakGroupScorer.cpp:1251-1279`).

- [LOW] **F02’s statement that match/assemble work is unchanged is too absolute.** Every precursor’s cell updates remain, but each additional chunk repeats per-spectrum dispatch/bucket overhead and can alter zeroing/reservation behavior (`src/extract/ChromatogramExtractor.cpp:964-1106,1116-1252`).

- [HIGH] **The full-scale ID rankings are not auditable from the pasted evidence.** Claims sourced to `doc/77`, `doc/82`, `FORENSICS_VERDICT.md`, `odia-fragvec-passes-six-columns.md`, and the other `odia-*.md` files are not pasted, so the 97.7%, 92%, 64.1%, AUC, oracle and seed-swing numbers cannot be verified. The pasted source supports only the structural facts of at most three candidates and no score-changing RT refit (`src/OpenDIAlyzer.cpp:1131-1135`; `src/score/PeakGroupScorer.cpp:2721-2727`).

- [FINE] The central source-backed diagnoses are correct: the cap is mean-cell-derived, the sink executes after the worker barrier, MS1 construction is serial, and Gate C consumes arrival-ordered decoys (`src/extract/ChromatogramExtractor.cpp:679-699,1250-1268`; `src/extract/Ms1Traces.cpp:129-170`; `src/score/PeakGroupScorer.cpp:1248-1281`).
