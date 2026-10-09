# Adversarial review brief: OpenDIAlyzer 2 (ODIA) v0.5.0 — why 657 GB, why 9.5 h, why 53.7% of DIA-NN's IDs

## READ THIS FIRST — your tool situation

You have NO working shell or file access in this session: every tool call fails with an IPC decode error. Do NOT attempt to run commands, open files, or search a repository. Review ONLY from the text pasted in this brief. Cite the pasted `file:line` numbers (the line numbers are the `cat -n` numbers shown at the left of every pasted source line, and they are the line numbers of commit bd73e65). If a claim in the diagnoses cites a line that is NOT pasted here, say "not pasted, cannot verify" rather than guessing. Do not invent citations.

You are reviewing, not editing. Your output is consumed by a program; be dense, rank by severity, give file:line for every claim, and when something is fine say so in one line.

## 0. CONTEXT (measured)

OpenDIAlyzer 2 (ODIA) v0.5.0 = commit bd73e65. C++ DIA proteomics search engine built on OpenMS, reads mzPeak (parquet-in-zip) files, extracts fragment chromatograms over a predicted library, scores peak groups with a semi-supervised GBT, target-decoy FDR.

Measured on IH1 diaPASEF (Bruker timsTOF; 32,210 (spectrum, window) entries over 17,448 physical frames, 1,343 MS1 frames, ~5,400 s gradient), library dn_pred_cam (9,922,669 precursors incl. 4.96M decoys; ~12 transitions each), same file / library / node as DIA-NN 2.0:

| | DIA-NN 2.0 | ODIA best arm (c3_cent_gate_s43) |
|---|---|---|
| wall | 14 m 22 s | 9 h 31 m 43 s (34,303 s) |
| CPU | 7 h 34 m | 36 h 32 m (131,524 s; 11,150 s of it system) |
| cores busy | 31.6 of 48 (66%) | 3.8 of 64 (6%) |
| peak RSS | 28.1 GB | 657,612 "MB" (OpenMS prints MiB; = 642.2 GiB). Within 407 MB of that on EVERY arm regardless of config |
| precursors at q<=0.01 | 37,334 (entrapment FDP 1.134%) | 18,995 at own q<=0.01; 20,040 = 53.7% at matched entrapment FDP |

Typical arm flags: `-threads 64 -live_memory_gb 400 -rt_window_pass1 400 -rt_window 110` (sometimes 300). Two passes: pass 1 extracts wide (±400 s) to harvest RT/IM anchors, pass 2 re-extracts narrow (±110 s) and scores.

Reference log c3_cent_gate_s43.log (shared/libv2/analysis77/), the lines the diagnoses cite, verbatim:

```
:46  MS1 traces: 1343 bins over 9922669 precursors, 8556564 with signal (86.2%), 50835.2 MiB, in 5589.599765 s
:49    decode 408.018 s, index 19.5254 s, match 5085.81 s, assemble 498.984 s, sink 14637.4 s
:50    points 55690608059 (33846277734 nonzero), all at once would be 207.464 GiB
:51    peak live 3324888 precursors, 114338442565 points (425.944 GiB), bound by precursor cap (3324888), 3 chunks; 5704152 would have been live at once
:52    live budget: budget 400 GiB / 129176 B per live precursor (5 planes) -> cap 3324888
:53    3 chunks, 59902 spectra decoded against 32210 in the run
:54    24576 precursors predicted to elute outside the run
:55    1275674 precursors covered by no isolation window
:59  fragment mass accuracy as a sub-score: 22596996 of 22596996 candidates carried a deviation, centred on -0.176 ppm
:61      of which gate C   0/0
:185 extracted 117152047 transitions of 9922669 precursors from 32210 spectra in 6.50402e+06 ms
:186   decode 183.591 s, index 23.7458 s, match 2004.71 s, assemble 76.4832 s, sink 4203.48 s
:187   points 15930173533 (9022684879 nonzero), all at once would be 59.3445 GiB
:188   peak live 1886108 precursors, 17804697520 points (66.3277 GiB), bound by retention-time overlap (1836986 of 8622717 precursors live at once)
:189   live budget: budget 400 GiB / 36949 B per live precursor (5 planes) -> cap 11623943
:226 scored on the fly 21765227 peak groups with gbt in 7.08407e+06 ms
:227   10867441 target / 10897786 decoy groups
:228   identified 18995 precursors at q <= 0.01
:231   2415473 precursors yielded no candidate peak group
:232 OpenDIAlyzer took 09:31:43 h (wall), 1d 12:32:04 h (CPU), 03:05:50 h (system), 1d 09:26:14 h (user); Peak Memory Usage: 657612 MB.
```

Rules the diagnoses were written under (and that you should hold them to): MEASURED (a run produced it; log/doc named) vs INFERRED must be separated; a hypothesis a doc closed by measurement may be re-opened only explicitly; ODIA is bitwise deterministic run-to-run so byte-identical outputs (`cmp` of the scored TSV) are the test for "output-identical" changes; comparisons of ID counts must be at matched entrapment budget N(e), never on q<=0.01 alone (q<=0.01 rests on 6-10 decoys and swings 12-17% on one decoy).

## 1. THE FOUR DIAGNOSES (written independently by four analysts; you are to attack them)

### Diagnosis A (memory-first)

UNIT. OpenMS prints peak RSS in MiB labelled "MB": 657,612 MiB = 642.2 GiB (c3_cent_gate_s43.log:232).

MEMORY (lens). Ranked by bytes, with status.
1. The 400 GiB budget is SATURATED, not a footprint — MEASURED. c3:51-53: "peak live 3324888 precursors ... 425.944 GiB, bound by precursor cap, 3 chunks; 5704152 would have been live"; cap = budget/(mean_cells×4 B×planes) (ChromatogramExtractor.cpp:679-692). 5.70M of 8.6M assigned precursors overlap at ±400 s (66%, as the header predicted, .h:619-622), so every full arm chunks and fills the cap. Peak is 6.5% over budget because the cap is inverted from MEAN cells (6,459) while the live set at the peak is all full-width (114,338,442,565/3,324,888/5 = 6,878 cells). Same line on all 40 dn_pred_cam arms (656,227-665,739 MiB). Dose-response exists: human_v2 at 350 vs 400 GiB → 456,506 vs 533,883 MiB (rep_human2.log:186 vs r5_p1.log:186), same 3 chunks, same 59,902 decodes, wall +2.5% (inside the 2% noise floor).
2. Of the 5 planes, 2 are byte-for-byte copies of the intensity plane — INFERRED from source, high confidence. Under the default Aggregate::Sum (.h:609; OpenDIAlyzer.cpp:1783-1785) the match loop does `at += intensity` (:1218), `ppm_den += intensity` for intensity>0 (:1231) and `im_den += intensity` for intensity>0 and finite 1/K0 (:1241) — identical sums, same thread, same order. On IH1 every peak carries 1/K0, so base ≡ ppm_den ≡ im_den; on Astral the two im planes are all-zero yet allocated ("5 planes" e2e_astral.log:49 beside "NO MOBILITY AXIS" :111). 40% of 425.9 GiB (170 GiB) stores nothing.
3. ~165 GiB (25.6%) unattributed — MEASURED by subtraction: 642.2 − 425.9 − 49.6 (MS1, c3:46) − 1.95 (library, c3:7). human_v2 leaves 59.6 GiB, r6_z234 79.5, e2e_astral ~290. Known pieces: retained PeakGroups 437 B each (MEASURED: p1_c6 vs p1_a1 +21.26M groups → +8,884 MiB; p1_k12c6 vs p1_k12 reproduces at 435 B) = 9.2 GiB for 22.6M groups — this refutes the "~5 KB/group" inference; decode working set 2.9 GiB (.h:678-682); bookkeeping ~2 GiB. ~150 GiB remains. Leading INFERRED candidate: BlockPool over-reservation. free_ is keyed by exact n = valid×cycles (:221), owned_ never shrinks (:230-232); every start-truncated precursor (lo=0) activates in the first batch (:1128-1132) with an odd size, expires, and its block never matches an interior full-width request, so reserved_ ≈ peak_live + truncated mass. reserved_ is tracked (:247) and printed nowhere. doc/34:80-86 named this, unproven. The dn-vs-hv contrast (165 vs 60 GiB at near-equal live sets) must be explained by whatever the instrument finds.
4. MS1 dense matrix 9,922,669×1,343×4 B = 49.6 GiB — MEASURED (c3:46; Ms1Traces.cpp:82); 1,275,674 uncovered + 24,576 out-of-run precursors (c3:54-55) own 6.6 GiB of rows that are never read (no blocks, ChromatogramExtractor.cpp:535-540). Checked against but not charged to the budget (OpenDIAlyzer.cpp:1860-1882).
5. Pass-1 volume is the only reason for >120 GiB: pass 2 at ±110 s is RT-overlap-bound at 66.3 GiB (c3:188). Pass 1 extracts all 9.9M incl. 4.96M decoys at ±400 s to harvest anchors; `-pass1_precursors` (target count→stride, OpenDIAlyzer.cpp:644-657, 2382-2395) was never used at scale.
Closed, not reopened: decode floor (now 2.9 GiB), page cache, per-thread buffers, IH1 sparsity (60.8% nonzero in pass 1, c3:50; Astral is 8.5% nonzero, e2e_astral.log:47 — instrument-specific, noted not proposed), tcmalloc/arena at the 10 GiB scale (doc/11). Allocator is re-opened ONLY as an instrument at ≥400 GiB.

Which single change: zero-code, `-live_memory_gb 150` removes the most bytes (−250..−340 GiB INFERRED) and on dn_pred_cam should be cmp-identical (Gate C tau=0, c3:61). Best code change: alias the two den planes (−170 GiB at identical cap, bit-identical by construction). Caveat MEASURED on human_v2: pass-1 group count moved by 1 (2,689,293 vs 2,689,292) between 350 and 400 GiB — Gate C's tau is calibrated from the first n_needed decoys in arrival order (bd73e65 PeakGroupScorer.cpp:1248-1279), and chunking changes arrival order; provenance of rep_human2 is unrecorded (run_r5.sh:15-19), so this is a pre-registered check, not a verdict.

CPU/WALL — MEASURED timers, serial status from source. Sink 18,841 s = 54.9% (c3:49,186), 100% serial: emit loop runs after pool_impl.run returns (:1250 → :1259-1268) and calls sink.accept (:941) on the driver thread. MS1 build 5,590 s = 16.3%, no threading (Ms1Traces.cpp:128-170). Match 7,090 s = 20.7% on the 64-thread pool at ~20% efficiency (INFERRED from CPU arithmetic). Serial ≈73% of wall → Amdahl ceiling 1.37x; 3.8 of 64 cores busy is this fraction, not a thread bug. Pass 1 alone is 60% of wall and exists only for calibration. Chunking is cheap: each extra chunk re-decodes 0.43 of a run (2 chunks 1.43x r6_z234:54; 3 chunks 1.86x c3:53) ≈ 95 s.

IDs — from the ids map, all MEASURED at full scale: admission 97.7%, correct peak in list 92%, but ranking inside depth 37,334 = 64.1% (doc/77, doc/82); the dominant loss (~33 pts) is evidence content at the correct apex — fragment co-elution (ref corr 0.85 vs 0.37, FORENSICS_VERDICT) — not localisation (+0.3), not admission (2.3), not the classifier (≤+3). Multi-candidate precursors mediate; library 1/K0 residual >0.030 is the one exogenous dose-response. The decoy null inflates the nominal count 1.56-2.24x; honest 1% reads below 18,995. None of the memory fixes below touch any of this; fixes 7-8 can and carry the frozen-RT gate.

### Diagnosis B (scheduling-first)

LENS: wall clock and parallel efficiency. Subject bd73e65; the four extractor files cited are md5-identical between HEAD and bd73e65 (verified), scorer files cited from `git show bd73e65:`. Reference arm c3_cent_gate_s43 (shared/libv2/analysis77/c3_cent_gate_s43.log:232): wall 34,303 s, CPU 131,524 s, system 11,150 s = 3.83 cores busy of 64 (6.0%).

MEASURED wall decomposition from the arm's own stage timers: library 44.1 s (log:8); MS1 trace build 5,589.6 s (log:46); pass-1 extraction 20,704.8 s = decode 408.0 + index 19.5 + match 5,085.8 + assemble 499.0 + sink 14,637.4 (log:48-49); pass-1 finish() ≈225 s (26,588.3 s bracket at log:89 minus MS1, extraction and the ~69 s mass probe; the bracket is `extractAndScore_`, OpenDIAlyzer.cpp bd73e65:4724-4743, so the "scored on the fly" line is NOT a scoring time); IM calibration 382.5 s (log:147); pass-2 extraction 6,504.0 s = decode 183.6 + index 23.7 + match 2,004.7 + assemble 76.5 + sink 4,203.5 (log:185-186); pass-2 finish ≈198 s (log:226 bracket 7,084.1); untimed tail (anchor selection, RT fit, mobility centring, 21.8M-row TSV) ≤586.5 s by subtraction. The "(602.1 s)" at log:96 is the MAD keep-width in seconds, not a time (OpenDIAlyzer.cpp bd73e65:2630-2634); the siblings map's "anchor selection 1.8%" is wrong.

Thread model (SOURCE). One worker pool (src/extract/ChromatogramExtractor.cpp:762-828, started at :828) is awake only inside `pool_impl.run(work, threads)` (:1250) for the per-spectrum match (:1139-1249). Everything else in extract() runs on the driver thread while 63 workers sleep on the condvar (:785-790): decode `source.peaks` (:1104), activate/zero-fill (:1131, `std::fill_n` :234), and the emit loop (:1259-1268 per batch, :1287-1292 flush) → `sink.accept(trace)` (:941) → `Sink::accept` → `Session::add` (include/odia/PeakGroupScorer.h bd73e65:1256; src/score/PeakGroupScorer.cpp bd73e65:1362) = picker plus every sub-score. The header documents the invariant: "`add` runs on one thread ... what already makes `Session::result_.groups` safe to push to" (PeakGroupScorer.h bd73e65:923-930). `Ms1Traces::build` has no thread/omp token (src/extract/Ms1Traces.cpp:129-170; per-peak `lower_bound` over a 9.9M-target index :140, random RMW into a slot-major 50.8 GiB matrix :162) and runs before the pool exists. OpenMP exists only in finish(): gbt.h bd73e65:286,505 (histogram loop capped at nchunk ≤32, :492), lda.h bd73e65:601-604 (inner_threads = omp_get_max_threads()/folds).

Does ODIA use the threads it is given? Yes: `-threads 64` → `options.threads` (OpenDIAlyzer.cpp bd73e65:1793, :4480) → pool; TOPPBase.cpp:411 → :96 `omp_set_num_threads`. The 6% is a serial-fraction fact, not a plumbing bug. Latent only: `max(1, threads)` at :1793 makes `-threads 0` a 1-worker pool while OMP gets all cores.

Ranked root causes, CPU/WALL:
1. MEASURED: 78% of wall is single-threaded by construction — sink 18,841 s (54.9%), MS1 build 5,590 s (16.3%), decode/assemble/index/IM+mass calibration/load/tail ≈2,300 s (6.7%). Amdahl ceiling on this structure with infinite threads: 34,303/26,723 = 1.28x. Replicated on ora_ctl/ora_on/samelib3 (sink 13,760-13,844 s, MS1 4,340-5,621 s; logs :45,:48) and on human_v2 hv_cap0 (sink 4,759+2,324 s, MS1 4,261 s of 17,973 s wall; hv_cap0.log:46,49,176,220). Fixture bench_d2t240n_ih1: sink 206 of 459 s pass-1 extraction, CPU 812% (bench_d2t240n_ih1.log:46,154) — the serial share GROWS with scale (45% → 71% of pass-1 extraction), so fixture speedups understate full-scale wins.
2. INFERRED (CPU arithmetic, no per-stage CPU timer exists): inside the only parallel region efficiency is 17-21%: 131,524 − ~26,700 serial − ≤27,000 finish() → match 78-95k CPU-s over 7,090 s wall = 11-13 effective cores of 64. Candidate causes: MATCH_BATCH 128 (:848) = 2 spectra/thread/barrier with heavy per-spectrum variance (co-packed diaPASEF frames); five RMW planes (:1218-1241) into a 426 GiB random-access live set; system time 8.5-17.7% of CPU. The only thread-scaling datum (doc/BACKLOG.md:1469, 16→96 threads 1.9x) is an old decode-dominated binary; the planned curve (BACKLOG.md:1473) was never run.
3. INFERRED: system time 11,150 s (c3) / 19,134 s (ora_ctl.log:191) is page-fault traffic, not futex: on the two `time -v` arms faults scale with RSS (117.6M at 40 GB, 255M at 195 GB) while voluntary context switches do not (bench_d2t240n_ih1.log:163-165; bench_astralw269_astral.log:168-170); 5-11 faults per resident page in every record (doc/34:84 3.26e9) → re-faulting (glibc heap-top trim thrash on the driver thread's transient vectors, or NUMA hint faults on dax), not first-touch. Stage attribution has never been measured; if it lands on the driver thread it is up to ~30% of wall. dax's NUMA/THP config is unverified (ssh refused; ibminode05 is single-socket, numa_balancing=0, but the arms ran on dax).
4. MEASURED: sink cost is linear in extracted points: 14,637 s/8.62M emitted = 1.70 ms per precursor at ±400 s vs 4,203 s/8.62M = 0.49 ms at ±110 s (log:188 gives 8,622,717 assigned). Pass 1 (60.4% of wall) exists only to harvest anchors and mobility centres; its width is a direct wall AND memory lever.
5. MEASURED: not levers — classifier/finish ≈1.2%; decode 1.7% (dedup landed, src/io/MzPeakSource.cpp:372-402); IM calibration 1.1%; load 0.1%.

MEMORY (other lens; coupling only): peak is the saturated budget — cap 3,324,888 live precursors × 129,176 B (5 planes) = 425.9 GiB (log:51-52), 4 of 5 planes residual accumulators (match loop :1218-1241), dense MS1 50.8 GiB (log:46), ~135-165 GiB unattributed. A parallel sink adds at most one batch of deferred `blocks.give` (:945-959) to the live set; a narrower pass-1 window moves wall and peak RSS together.

IDs (other lens): no scheduling fix buys an ID; the gap is evidence content at the correct apex (odia_ids map). The wall fix IS the ID programme's throughput: a 3-5x faster arm is 3-5x more seeds/baselines per week.

Amdahl ordering (INFERRED, 64 threads, ODIA bitwise deterministic so each step is `cmp`-gated): (0) instrument + thread-invariance harness; (1) parallel sink −18.2 ks → ~4.5 h (2.1x); (2) parallel MS1 −5.3 ks → ~3.0 h (3.2x); (3) pass-1 110 s, behaviour-changing, → ~1.8 h (5.2x) behind the RT-convergence/matched-entrapment gate; (4) match 20%→50% → ~1.2 h (8x). Floor at current algorithms ≈25 min (decode + IM calibration + tail + finish), i.e. ~1.8x DIA-NN's 14 m 22 s (dn_ih1_cam.log:77).

### Diagnosis C (evidence-first)

Subject: ODIA bd73e65, IH1/dn_pred_cam, best arm c3_cent_gate_s43 (analysis77/c3_cent_gate_s43.log:228 "identified 18995", :232 wall 09:31:43, CPU 1d12:32, "Peak Memory Usage: 657612 MB"). DIA-NN 37,334 at 1.134% entrapment FDP (doc/82:77-81). 53.7% = N(62) 20,040 (odia-vs-diann-measured.md:23).

IDS, ranked by share of the missing 46%:
1. Ranking with the correct peak already in hand (~33 pts, ~12,400 precursors; INFERRED decomposition from MEASURED stages). Admission 97.7% (doc/77:22), correct peak in the ≤3-candidate list 92% (doc/77:80), forced-apex oracle +0.3 pts (odia-oracle-localization-is-minor.md:90-95), recall at depth 37,334 only 64.1% (doc/82:89,103-107). Mechanism MEASURED: fragment co-elution at DIA-NN's own apex — per-fragment ref corr 0.85 (recovered) vs 0.37 (buried), AUC 0.932 inside DIA-NN's own quantity deciles; 12/12 fragments present in both classes (FORENSICS_VERDICT.md:16-25,59-70,524-525). Oracle plateau 65.2% on the 22 sub-scores, invariant to six recipe perturbations (doc/77:47-62) → feature content, not classifier. At the decision point the 42 columns reach matched-score AUC 0.811 and the co-elution family ≤0.61 (odia-matched-score-is-the-test.md:139,159-162). FRAGVEC (78 per-fragment columns vs a data-driven reference, DIA-NN's mechanic) lifts that to 0.893 and buys +298/+232 concordant ids at N(62) with 49/43 entrapments vs 62 — MEASURED offline at full scale, never wired into the scorer (odia-fragvec-passes-six-columns.md:103-116; bd73e65 msg 39-49).
2. Multi-candidate precursors: ≥2 candidates → missed 85% (0.915 adjusted), 41.6% vs 7.3% in the top abundance decile; LOBO "mediator, not cause" (FORENSICS:35-36, block table). Candidate SELECTION headroom: rt_dn oracle 93.9/87.4 vs deployable 87.4/73.9 cohort-projected (odia-tensor-verdict.md:62-75). Structural gap INFERRED from source: ODIA picks ≤3 candidates once per pass (OpenDIAlyzer.cpp:1131) and the semi-supervised loop only re-ranks those rows (lda.h:1142-1166; `refitsChangeScores()` false, PeakGroupScorer.cpp:2721-2727), while DIA-NN re-seeks with refitted weights across ~12 in-RAM iterations (DIA-NN-workflow-handoff.md:445-452,464-573). The deciding A1/A2 measurements were named and NOT run (FORENSICS:248-274).
3. Library 1/K0 residual >0.030: 2,967 precursors missed 80% at matched abundance AND matched DIA-NN evidence (FORENSICS:34-35; dose table) — the only clean exogenous handle. ODIA's own IM calibration recovers <40% of the scale error (odia-mobility-calibration-is-truncated.md:218-229); fine-tuned CCS reaches DIA-NN's run-refit (held-out sd 0.01475 vs 0.01480; odia-refined-library-diann-arms.md:136-149) and has never been run through ODIA. Every IM-column edit so far failed through the fragment-mass probe coupling (affine −34.7%, odia-im-affine-correction-wins.md:89-97; prior −5.81%, odia-width-default-regression.md:360).
4. Scorer regime instability (reliability, not a steady loss): seed 42→43 moves q≤0.01 16,965→4,105 with identical evidence (odia-scorer-regime-flips.md:278-285); cause a 1:155 positive:decoy imbalance and unbounded early leaves (odia-leaf-cap-b2.md:97-101,176). Cap +387 ranking-only but the human_v2 replicate FAILED through Gate C's arrival-order null (odia-leaf-cap-b2.md:151-158; odia-fixture-cannot-arm-gate-c.md:30-34). Intercept 0 MEASURED: +1.62% vs w1 but −1.58% vs b2_cap1, 333/333 negative (SESSION_NOTES:5687-5689).
5. FDR calibration: own q≤0.01 sits at 1.84-2.31% true FDP (odia-vs-diann-measured.md:24), over-count 1.56-2.24x and depth-dependent (odia-decoy-null-borrowed-evidence.md:129-134). It INFLATES the nominal count; it does not hide IDs — recalibration changes nothing in N(e) (ibid:61-63). Correction to the task premise: fixing the estimate alone changes the reported count, not the matched-FDP count.
Closed, not re-opened: m/z and IM windows both directions, admission/localisation, Gate C relaxation, candidate generation, training knobs, XGBoost, ownership mask (1.6-1.7% exposure), trace tensors, MS1 intensity.

MEMORY (MEASURED 657,612 MB, 400 MB spread across arms): pass 1 at ±400 s over 9.92M precursors (4.96M decoys) would have 5,704,152 live; cap 3,324,888 = 400 GiB / 129,176 B / 5 planes (c3 log:51-52; ChromatogramExtractor.cpp:680-699) → 425.9 GiB live blocks, 3 chunks (66% of RSS); MS1 dense matrix 50,835 MiB (log:46); library 2 GiB; ~165 GiB unattributed (INFERRED: BlockPool size-class retention, 22.6M PeakGroups with heap sub_scores, glibc non-return). 4 of 5 planes are residual accumulators; 2 are INFERRED byte-copies of the intensity plane (ChromatogramExtractor.cpp:1218-1241). The peak is the budget knob, not the data.

CPU/WALL (MEASURED): sink 14,637 + 4,204 s = 55% of wall, serial on the driver thread after `pool_impl.run` returns (ChromatogramExtractor.cpp:1250,1260-1267; log:49,186); MS1 build 5,590 s serial (log:46); match 7,090 s in the only pool at ~20% efficiency (INFERRED from CPU/wall); pass 1 = 60% of wall for 3.8M anchors. 78% single-threaded → 3.8/64 cores; Amdahl ceiling ~1.3x on threads. Decode (592 s) and the classifier (~420 s) are not levers.

Ordering by IDs per week: FRAGVEC wiring (days, measured channel) > fine-tuned library with pinned mass path (days, exogenous handle) > Gate C null correctness (days; unlocks replicates and 2/3 of the sink) > reselection (weeks; stage-0 replay in hours decides) > training-population stabilisation > honest q (zero IDs, mandatory). Cost fixes raise the arms-per-node-day that every ID experiment is rate-limited by.

### Diagnosis D (contrarian)

CONTRARIAN DIAGNOSIS of ODIA v0.5.0 (bd73e65) on IH1/dn_pred_cam. M = measured (log/doc named), I = inferred.

MEMORY (657,612 MB = 642 GiB), ranked
1. (M) The peak is the `-live_memory_gb 400` flag saturated, not a footprint: cap = 400 GiB / 129,176 B = 3,324,888 live precursors, 3 chunks, 425.9 GiB (c3_cent_gate_s43.log:51-53; derivation ChromatogramExtractor.cpp:679-695). A survey of all 49 full-library logs shows every completed scorer arm ran at budget 400; the one sub-400 arm (cohort_plain, 350 GiB -> 651,404 MB; analysis77/cohort_plain.log:52-53,189) carried -out_chrom/-out_ms1_iso exports (run_cohort_plain.sh:20-22), which add ~58 GB on own_ctl/x1_fragvec_r1 (715,750/715,365 MB). The budget->RSS relation at full library has ONE clean point; "within 407 MB on every arm" describes the flag.
2. (M, source) Two of the five planes are bit-identical copies of the chromatogram plane: `at += intensity` (:1218) vs `ppm_den += intensity` for intensity>0 (:1231) and `im_den += intensity` for non-NaN IM (:1241), same thread, same order. Under the default Sum, ppm_den == base cell for cell; im_den == base whenever every peak carries 1/K0 (diaPASEF). 170 GiB of the 426 GiB live set is redundant. But: at a fixed BUDGET the cap simply rises (5.54M < 5.70M overlap -> still cap-bound, 2 chunks), so the RSS win exists only with the cap pinned or the budget lowered.
3. (M by subtraction, I in split) 165 GiB remainder (642 - 426 - 49.6 - 2). It scales ~0.4-0.5x live (doc/34:62-71: 485 of 1,465 GiB) and moves 1:1 with retained scorer state (c6 arms +8 GB, export arms +58 GB) -> consistent with VmHWM being reached at/after pass-1 finish(), with the freed BlockPool arena not returned (glibc 2.39, nothing else linked) and finish()'s large vectors mmapped on top. Nobody has measured WHEN the peak occurs; `reserved_` is tracked (:232) and has had no caller (:247) since 2026-08-05 (doc/34:80-86 names this as unproven).
4. (M) MS1 dense matrix 50.8 GiB (c3:46) is 7.7%; not a lever, and -no_ms1 removes one of the two signal-carrying feature families.

CPU/WALL (9h31m wall, 36.5 h CPU, 3.1 h system), ranked
1. (M) 77% of wall is single-threaded by construction: sink 14,637+4,203 s (c3:49,186; emit() runs after pool_impl.run returns, :1250 -> :941 -> PeakGroupScorer::Session::add :1362), MS1 build 5,590 s (c3:46; Ms1Traces.cpp:130-170, no threading), assemble/index/decode/IM-probe/anchors ~2,200 s. Amdahl floor for ANY thread count is 7.4 h; "6% efficiency" is this serial fraction. Across 44 byte-identical arms the sink is 13,520-15,060 s (+-5%, thread- and load-invariant) while match spans 2,937-6,292 s (2.1x) -> the parallel region is memory-bound and contention-sensitive (I); a parallel sink inherits that risk.
2. (M) Pass 1 is 60% of wall and 100% of memory above 66 GiB (pass 2: c3:188), and its measured contribution in c3 is small: the library already carried the external linear map (c3:10; run_c3.sh:112), the fitted map moved the anchor median 45 s with p95 residual 363 s ~ the window (ora_ctl.log:90-91; c3:100), the IM calibration left robust scatter 0.0098 -> 0.0098 (c3:145), centring covered 1.2M of 7.7M entries (c3:98) for a per-charge-mixed +3.94%. The only narrow-pass-1 full arm (odia2_dnlib, -rt_window_pass1 120, run_dnlib_seeded.sh:23) died under 50x oversubscription (odia2_done.log "dnlib exit 13"; MS1 took 15,227 s, odia2_dnlib.log:36). OPEN, never refuted.
3. (M) System time is 8.5-18% of CPU (c3:232; ora_ctl.log:191); on small arms it scales with RSS (4.0% at 40 GB -> 13.9% at 195 GB) with flat context switches -> page-fault traffic, not futex (I). 3.26e9 faults on the 1.4 TiB arm (doc/34:84) is 8.5 faults/page, not first-touch; NUMA balancing on dax is unverified (ibminode05: off, THP madvise).
4. (I) Match runs at ~12-15 effective cores of 64: 5 scattered RMW per matched peak into a 426 GiB set. Only thread-scaling datum on record is 16->96 = 1.9x on an old binary (BACKLOG.md:1468-1470).

IDs (53.7% at matched count), ranked
1. (M) The headline is lenient: N(62) matches the false-discovery COUNT at ODIA's depth (FDP ~2.1% there); own q<=0.01 sits at 1.84-2.31% true FDP (odia-vs-diann-measured.md:24); at 1.134% ODIA reads below 18,995 and the cell rests on 1-12 entrapment events (odia-matched-fdp-is-depth-allowance.md:8-10). Real gap >= 2x; "fixing the null" lowers the number.
2. (M) 33 of 46 points are ranking with the correct peak in hand; a DIA-NN-labelled oracle on ODIA's columns plateaus at 65.2% (odia-feature-plateau-65pct.md:8-15). Every scorer-side lever at full scale is <= +10% at matched budget; FRAGVEC moves matched-score AUC 0.811 -> 0.893 but only +298/+232 at common depth (odia-fragvec-passes-six-columns.md:54-61) and no oracle-at-depth read for it exists.
3. (M) The only exogenous dose-response is library 1/K0: production slice +-0.025 vs median library error 0.0235 (odia-matched-coords-delete-aperture.md:14-16) puts the median true mobility at the slice edge; 41% of pick losses are >0.025 off (odia-pick-losses.md:18-21); calibration recovers <40% of the scale error; width refuted both ways -> centring, not width. Fine-tuned CCS (sd 0.01475 = DIA-NN's run-refit) was never run through ODIA.
4. (M) The extraction-parity verdict (doc/44:6-12: 0.270 vs 0.311) was 998 precursors non-degenerate in BOTH tools with ODIA's own +0.0170 seed -- a both-found sample. The buried class (ref corr 0.37 at DIA-NN's apex) has never been compared trace-vs-trace against DIA-NN; `--xic` export exists (diann_readme.md:1187). I re-open "extraction exonerated" only for class C and only at production coordinates.
5. (I) Structural: one picking pass, no re-seek with refitted weights, no interference-deletion phase -- untested in-engine; the P1 port's DIA-NN overlap stayed flat (63.8 -> 64.0), so "do what DIA-NN does" has so far bought complementary IDs, not DIA-NN's.

## 2. THE CANONICAL FIX LIST (merged from the four diagnoses; "votes" = how many of the four analysts proposed it)

[F01] Pass-1 volume: narrow -rt_window_pass1 (110/200), stride -pass1_precursors 2000000, and the -passes 1 bound — one variable per arm, judged on frozen-RT statistics then N(e) (target several, votes 4, effort hours)
  mechanism: Pass 1 extracts all 9.92M precursors at ±400 s solely to harvest RT anchors (3.87M used of 3.87M scored), IM anchors and pass-1 1/K0 centres; it is 60% of wall (20,705 s) and the only reason the live budget binds (5.70M overlapping > cap 3.32M → 3 chunks, 426 GiB live). Live bytes and sink scan positions are linear in window width × precursors (pass-1 129,176 B vs pass-2 36,949 B per precursor; 4.66G vs 1.28G positions; sink 14,637 vs 4,203 s). Arm a: ±110 s (the MEASURED pass-2 analogue: ~6,500 s, 66 GiB, 1 chunk, 32,210 decodes); arm a': ±200 s (library already in run seconds via the supplied iRT map, c3:12); arm b: stride 4 (~1.43M overlap < cap, 1 chunk, ~190 GiB at 5 planes); arm c: -passes 1 with the external linear map, which bounds what pass 1 buys at all. Pass 2 unchanged in a/a'/b.
  evidence: MEASURED: c3_cent_gate_s43.log:48-53 (pass-1 20,705 s, 425.9 GiB, 3 chunks, 59,902 decodes) vs :185-189 (pass-2 6,504 s, 66.3 GiB, 1 chunk); :96-100 (anchors 3,869,789/3,869,803; map p95 362.9 in-sample / 362.6 held-out); :98 (centring 1,211,024 of 7,739,024); :97 (1/K0 anchors 15,484); c3:10,12 (external map applied, library in run seconds); ora_ctl.log:85-91. Source: OpenDIAlyzer.cpp:644-657 (-pass1_precursors is a target count → stride), :2382-2395 (stride, indices not renumbered), :1032 (-passes 1), :5798-5860 (-im_center_from_pass1 reads each precursor's own pass-1 group); ChromatogramExtractor.h:614-622. Rules: odia-rt-calibration-frozen.md (acceptance = map convergence, not IDs); odia-verify-the-negative.md:17-20 (stride 100000 → 0 IDs was a flag artefact). The only narrow-pass-1 arm (120 s, run_dnlib_seeded.sh:23-24) died under 50x oversubscription (odia2_done.log 'dnlib exit 13', odia2_dnlib.log:36) — never refuted. doc/34:158-163 ranked recovering pass-1 volume #1 in Aug.
  expected win: INFERRED from the measured pass-2 analogue: arm a pass-1 extraction 20,705 → ~6,500 s (wall 9.5 h → ~5.5 h), live 426 → ~67 GiB, Peak RSS ~657 → ~250-300 GB, 1 chunk; arm b peak ~280k MiB (~200k with F03), wall ~5.5 h; arm c (if pass 1 is not earning its cost) wall ~3.5 h, RSS ~200 GB. Zero code for all arms. The largest cost lever available tonight; combinable with F02-F05 only after each is read alone.
  risk: Behaviour-changing. Anchors with |map error| in (110..400] s are lost — exactly the population anchor statistics cannot see (frozen note); stride starves -im_center_from_pass1 (75% uncentred), the per-charge 1/K0 anchor floor (120) and -mass_width_from_ids; -passes 1 loses centring (+3.94% vs own control) and the IM calibration anchors. Calibration changes re-land the fold partition → head ±5-9%; read at matched entrapment budget with two baselines (b2_cap1, w1_ctl70), never on q≤0.01 alone. Not a re-opening of the frozen RT decision (cost arm read on its statistics) nor of the REFUTED fragment m/z narrowing (different window). NOT compatible with F02's chunk-invariance read in the same arm (one variable).
  verify: Arms a (-rt_window_pass1 110), a' (200), b (-pass1_precursors 2000000), c (-passes 1, drop -rt_window_pass1/-im_center_from_pass1, keep -irt_slope 7.9012 -irt_intercept 759.3118 -rt_window 110 -im_prior off) vs a byte-identical c3 control and w1_ctl70, same binary (sha a568ddf0), idle node (ibmi-nodes.sh first), read 'pass 1 extracts every Nth' line before anything. Frozen-RT gate first: 'fitted the retention-time map ... p95' within +5% of 362.6 s (ideally 0.5 s), refine residuals equal to control, anchor count ≥50% of control, per-charge 1/K0 anchors ≥120 and IM gate PASSED, centring line ≥50% of 1,211,024 (a/a' only). Cost reads (all must hold or the arm is mis-specified): peak live ≤70 GiB (a) / ≤100 GiB (b), '1 chunk'/32,210 decodes, pass-1 sink ≤5,000 s, Peak Memory ≤300,000 MB, wall ≤6.5 h. Benefit: pass-2 region median N(e) over FDP 2-10% ≥ −1% (a/b) / ≥ −2% (c) with ≥90% ordinals ≥ −2%, q≤0.01 ≥0.9x control, DScore histogram by signature not regime-flipped; report three numbers. FAIL: p95 worse >5% or N(e) median < −1% (a/b) → pass-1 width is load-bearing, lever moves to the stride; arm c median ≤ −5% → pass 1 earns its cost and a/b must keep the IM anchors.
  already tried: Width 400 inherited from the uncalibrated-library design, never varied at scale; 120 s arm died (odia2_dnlib, no result); -pass1_precursors used only at 100000 (99x stride → 0-ID artefact) and absent from every analysis77 arm (run_c3.sh:111-114); no -passes 1 scoring arm in 49 logs (cohort_* single-pass arms were -out_chrom exports at ±600 s); -rt_window 110↔400 on PASS 2 is byte-inert (odia-pick-losses.md:43-46, different flag).

[F02] Lower -live_memory_gb (150, then 100) on the c3 configuration — zero code; byte-identity of the TSV is the gate (target memory, votes 3, effort hours)
  mechanism: The live set is a saturated budget: cap = budget/(mean_cells×4×planes) (ChromatogramExtractor.cpp:679-692, verified in bd73e65), 5.70M precursors overlap at ±400 s so pass 1 always fills the cap and chunks greedily (:699-725). Lowering the budget lowers the live set linearly; each extra chunk costs ~0.43 of a decode pass (~95 s; a chunk needs its band plus one 800 s window of spectra). Sink/match/assemble work are unchanged (every precursor emitted once). RSS/budget has been 1.53-1.63x on every full arm, so the arm also measures whether the ~165 GiB remainder is budget-proportional (feeds F06).
  evidence: MEASURED dose-response on human_v2: rep_human2.log:52-54,186 (350 GiB → 3 chunks, 59,902 decodes, 456,506 MiB, 6:23:34) vs r5_p1.log:52-54,186 (400 GiB → 533,883 MiB, 6:14:11); c3_cent_gate_s43.log:51-53 (425.9 GiB, 3 chunks); r6_z234.log:52-54 (2 chunks → 46,056 decodes = 1.43x); analysis77/cohort_heavy.log:22-26 (7 chunks, 135,768 decodes in 770 s; match per point 73 ns vs c3's 91 ns — chunking did not inflate match); default budget is 20 GiB (OpenDIAlyzer.cpp:623-636) while all 44 scorer arms ran 400; doc/34:75-79 (budget under-counts RSS 1.6x); odia-timing-noise-floor.md:19-24 (2% wall noise).
  expected win: INFERRED: 150 GiB → Peak 657,612 → 300-390k MiB (−40..−55%), ~7 chunks, pass-1 decode 408 → ~800 s (+1.1% wall); 100 GiB → 220-350k MiB, ~9-12 chunks, +1,200-2,600 s decode (≤ +8% wall). IDs unchanged on dn_pred_cam where Gate C is inert (c3:61 '0/0'). 3-4 arms per node instead of 1.
  risk: Chunk-count invariance of the output is UNTESTED (chunk edges move with the cap): a TSV byte difference is a defect report, not a result. On Gate-C-armed libraries (human_v2) emission order feeds tau (PeakGroupScorer.cpp:1248-1282 verified: tau from the first n_needed decoys by arrival) — rep_human2 vs r5_p1 differ by 1 pass-1 group and 3 Gate-C admissions (provenance of rep_human2 unrecorded, run_r5.sh:15-19); F07 is the precondition there. If the remainder is fixed rather than proportional RSS lands ~350 GB.
  verify: c3 command verbatim with -live_memory_gb 150 (then 100), same binary sha a568ddf0, solo on an idle node, same-night control. PASS: Peak Memory ≤420,000 MiB (150) / ≤350,000 MB (100) AND `cmp` of the sorted scored TSV vs c3_cent_gate_s43.tsv clean AND pass-1 decode ≤1,300 s (150) / ≤1,500 s (100) AND wall ≤1.10x. FAIL: any TSV byte differs (fix F07 first, report as defect); Peak >500,000 (150) / >400,000 (100) → remainder not budget-proportional → F06 decides; wall >1.25x → stop at 200 GiB. Secondary read: system time vs 11,150 s (feeds F09).
  already tried: No full-library (9,922,669-row) arm below budget 350 has completed (49 logs); rep_human2 at 350 (launcher lost); cohort_* at 250/350 confounded by -out_chrom tees (run_cohort_plain.sh:20-21; cohort_heavy 998,856 MiB); 1M-precursor arms were RT-overlap-bound at 85 GB (b_base.log) and say nothing about the cap regime; auto mode (60% MemAvailable) unused by every arm.

[F03] Alias ppm_den (and, where exact, im_den) to the intensity plane under Aggregate::Sum (5 → 3 planes); turn im planes off on no-mobility runs; pin cap or lower budget to realise RSS (target memory, votes 3, effort hours)
  mechanism: Under Sum the base cell receives `at += intensity` (ChromatogramExtractor.cpp:1218) and ppm_den receives `+= intensity` for intensity>0 (:1225-1231) in the same thread and loop iteration — adding 0 is a no-op, so ppm_den is bit-identical to base. im_den (:1236-1241) adds only when !isnan(peak_im), so it is identical to base only when no matched peak carries NaN mobility — VERIFIED in bd73e65 source during this merge: the im identity is conditional and must be asserted, not assumed. Scorer reads both only as num/den (PeakGroupScorer.cpp:2116-2145, 2191-2212, 2291-2301). Allocate base, ppm_num, im_num; point ppm_den (and im_den if the assert holds) at base; skip the second give-back. On runs with no mobility axis, collect_im_residuals=false (the documented 'auto' semantic, OpenDIAlyzer.cpp:1404-1413) instead of two zero planes. Cap uses the real plane count (:682-684), so at a fixed 400 GiB budget the cap rises to ~5.54M and the win is 3 → 2 chunks, NOT RSS; pin -max_live_precursors 3324888 or lower the budget to 240 to convert it into bytes.
  evidence: ChromatogramExtractor.cpp:1216-1242 (accumulation, verified), :682-684 (planes in cap), :868-883 (activate takes 5 blocks), :884-960 (emit hands pointers), :948-958 (five gives); ChromatogramExtractor.h:609 (Sum default), :655-657; OpenDIAlyzer.cpp:1775-1785 (both on by default); c3_cent_gate_s43.log:51-52 ('(5 planes)', 425.9 GiB → 85.2 GiB per plane); e2e_astral.log:49 '(5 planes)' with :111 'NO MOBILITY AXIS' (two zero planes of 416.9 GiB on Astral); history runs the other way (ChromatogramExtractor.cpp:170-186: 1 → 2 planes after review; the 'real but not binding' registration text predates the 9.9M library).
  expected win: INFERRED: −40% bytes per live precursor: 425.9 → 255.6 GiB at the same cap (−170 GiB; −85 GiB if only ppm_den aliases); pool slack shrinks proportionally if it scales with cap (−60 GiB). Peak 657,612 → 420-485k MiB; stacked with F02 at 100 GiB of 3-plane blocks = 1.67x more precursors per chunk. Astral: an additional −167 GiB with the 'auto' fix alone. Output bit-identical.
  risk: Implementation only: Max mode must keep separate planes; aliased pointers must not be given back twice; no-IM runs must keep im planes null (aliasing im_den→base would turn observed_im from NaN into 0.0); if any matched peak has NaN 1/K0 on diaPASEF the im_den alias is wrong — the fixture assert decides. No ID risk when the TSV cmp is clean.
  verify: Stage 1 (fixture, ~20 min, mix10k/IH1): assert in emit() `ppm_den[i]==base[i]` and count `im_den[i]!=base[i]` — PASS: 0 ppm violations (any violation falsifies the identity, stop); im violations 0 → alias both, >0 → alias ppm_den only. Stage 2 (19,995-precursor acq_* fixture, old vs new binary): `cmp` TSV clean, log prints '(3 planes)'. Stage 3 (one full arm, one variable): -live_memory_gb 240 → log 'cap 3324888', '3 chunks, 59902 spectra', TSV cmp clean vs c3_cent_gate_s43.tsv, Peak ≤490,000 MiB; or -max_live_precursors 3324888 pinned → 'peak live' ≤256 GiB. Astral: rep_astral 30-min arm, TSV cmp clean, '(1 planes)' or '(3 planes)' printed. FAIL: any cell mismatch beyond the assert's allowance or any TSV diff → revert.
  already tried: None found; the denominator/base duplication was never noticed. Residual planes shipped as defaults (collect_mass_residuals, im_features auto) and never measured for identity. mzPeak/decode path, sparse chromatograms and allocator tunables CLOSED as memory levers at the 10 GiB scale (not re-opened here).

[F04] Score emitted precursors on the idle worker pool (parallel sink by ordered speculation) — after F07 makes Gate C order-independent (target cpu_wall, votes 3, effort weeks)
  mechanism: The emit loop runs after pool_impl.run returns (ChromatogramExtractor.cpp:1250 → :1259-1268, verified) and calls sink.accept (:941) → Session::add (PeakGroupScorer.cpp:1362) serially on the driver while 63 workers sleep (:785-790): 55% of wall. Each emitted precursor is independent work over its own contiguous block plus read-only MS1 rows (PeakGroupScorer.cpp:2423-2436). Compute candidates and sub-scores per precursor on the pool (per-thread Result; finish() already re-sorts groups into (precursor, apex_rt, apex_intensity) order, :3145-3157); commit admissions (Gate C tau, result.groups.push_back :2655, the unguarded `++result.precursors_without_candidate` :1405,:1502,:1507,:1513,:1644) on the driver in the original emit order; defer blocks.give (:945-959) by one batch. `rejects_` is already thread_local with a registry mutex (:1310-1327); terminal_reason writes are per-precursor-index and disjoint (PeakGroupScorer.h:916-930). A -out_chrom/-out_fragvec tee must stay ordered. 78% of wall is single-threaded so -threads alone gives Amdahl ~1.3x; this is where the 40x wall gap lives.
  evidence: MEASURED timers: sink 14,637.4 + 4,203.5 s = 54.9% of 34,303 s (c3_cent_gate_s43.log:49,:186,:232); replicated 13,760-13,844 s on ora_ctl/ora_on/samelib3 (:48); 1.70 ms (pass 1) / 0.49 ms (pass 2) per emitted precursor over 8,622,717 (log:188); sink invariant 13,520-15,060 s across 44 arms. Serial by construction: PeakGroupScorer.h:923-930 states the one-thread invariant. doc/34-runtime-memory-attrition.md:42-60,160-167 (same split measured Aug 2026, 'parallelise the sink' ranked #2, never built); doc/39:378-382 (deferred over tau: 'I side with serial until someone measures tau twice'); odia-is-deterministic.md:34-46 (thread invariance holds BECAUSE the sink is serial); sibling lessons: DIAspeXtract PERF-ATTEMPT-1 lost 7.2% of peptides to order-dependent shared state in an 'output-identical' parallelisation; MS1 path 444.7 → 38.8 s at 80.5x is the precedent; FASTag DirecTag 192-cpu serial tail.
  expected win: INFERRED: sink 18,841 s → ~600-2,000 s at 20-30x on 64 threads (MS1 random reads cap it below 64x); wall 34,303 → ~16,000 s (9.5 h → ~4.5 h, 2.1x) with pass 1 as is, ~2 h combined with F01/F05; CPU unchanged; cores busy 6% → ~25%; 2.5x more arms per node-day for the ID programme. Memory: +≤1 batch of deferred blocks (bounded, tens of GiB at most).
  risk: Highest-risk engine change: Gate C tau is arrival-order dependent (F07 precondition) and inert on dn_pred_cam and un-armable on the 10k fixture, so the order dependence is INVISIBLE on both — determinism must be built by construction, not tested into existence. A single admission change re-lands the fold partition and makes IDs incomparable (compare on bytes). The match region's ~20% efficiency suggests the node's memory system may cap the sink too (F08 sizes the ceiling). Joint thread budget with the classifier's nested OMP (lda.h:589-604 sizes finish() from omp_get_max_threads(); OpenBLAS linked with no OPENBLAS_NUM_THREADS; `-threads 0` clamps the pool to ONE worker at OpenDIAlyzer.cpp:1793/:4480 — fix in the same change). Do after F01/F02/F03, which may shrink the sink 3.5x first.
  verify: GATE first: md5-identical sorted score TSV at -threads 1/8/64 on the 23-min fixture with the new binary AND bd73e65 at 64, AND on one full dn_pred_cam arm, AND one human_v2 arm (where 'gate C null armed ... tau=' must print the same tau at -threads 16 vs 64). Fixture: pass-1 sink line ~206 s → ≤25 s at 64 threads. Full IH1 with c3's exact flags: PASS = TSV md5-identical to c3_cent_gate_s43.tsv, pass-1 sink ≤2,000 s (target ≤1,000), pass-2 sink ≤600 s, wall ≤5 h; task count during pass-1 sink = pool + 1 + Arrow threads with no OMP team alive. FAIL: any byte difference (correctness bug, not a performance result; revert); pass-1 sink >4,000 s at 64 threads → memory-bound, stop.
  already tried: Streaming sink shipped (a0989fc) but serial; doc/34 ranked it #2 (Aug 2026), doc/39 parked it on tau; thread-invariance harness listed 'Not yet' (BACKLOG.md:456-457); only thread-scaling datum 16→96 = 1.9x on an old decode-dominated binary (BACKLOG.md:1469); no parallel-sink arm exists.

[F05] Frame-block-parallel Ms1Traces::build with a direct-address log-m/z bucket index (output-identical) (target cpu_wall, votes 3, effort days)
  mechanism: The build walks 64-frame blocks serially (Ms1Traces.cpp:129-170) and writes `values_[slot*bins_+(b+s)]` (:162/:165): the column index is the frame, so two blocks never write one cell — block-parallel is race-free by the index expression and Max aggregation is order-free per cell. Only `resid` (:121-126, capped at 2M) is shared: collect per block, concatenate in block order, truncate to RESID_CAP → the logged median is byte-identical. Decode the block on the driver (one Index + one Spectra, MzPeakSource.cpp:454-455, not thread-safe), match in parallel. Second, separately measured step: replace the per-peak `std::lower_bound` over 9.9M sorted targets (:140/:149-151) with the log-m/z bucket index the MS2 match already uses (ChromatogramExtractor.cpp:81-145, `x.bucket`).
  evidence: MEASURED: MS1 traces 5,589.6 s = 16.3% of wall (c3_cent_gate_s43.log:46); replicated 4,340-5,621 s (ora_ctl.log:45, ora_on.log:46, samelib3.log:45), 4,261 s on human_v2 (hv_cap0.log:46), 15,227 s under oversubscription (odia2_dnlib.log:36); serial: grep for thread|omp|parallel|pool in Ms1Traces.cpp returns nothing; built once per run before extract() so the pool does not exist yet. OpenDIAlyzer.cpp:1898 (median residual is log-only). Sibling precedent: DIAspeXtract RUNTIME-PLAN:194-204, 242-256 (MS1 path 444.7 → 38.8 s at 80.5x, byte-identical) and CHANGELOG:813-822 (lower_bound → bucket).
  expected win: INFERRED: 5,590 s → 300-600 s at 10-20x (DRAM-latency bound on random RMW into a 50.8 GiB slot-major matrix, so not 64x); −5,000 s wall (−14%, ~1.4 h); combined with F04: 9.5 h → ~3.0 h (3.2x). The serial pre-pass that delays pass 1 on every arm disappears.
  risk: Low: cells disjoint by index; resid sample must keep first-2M-in-frame-order semantics or the printed median changes (diagnostic only). Memory unchanged (dense matrix stays; F12 addresses it). Load-gate the node (14% per-stage noise).
  verify: Fixture (255 MS1 bins, 143 s MS1 line, bench_d2t240n_ih1.log:43): PASS = MS1 matrix checksum identical to the serial build, MS1 line ≤20 s at 64 threads (≥7x), score TSV md5-identical at -threads 1/8/64. Full IH1, one variable: 'MS1 traces: ... 8556564 with signal (86.2%)' count identical, median residual identical, sorted TSV md5-identical to c3, MS1 line ≤600 s. FAIL = checksum/count/md5 differs (race, revert) or speedup <4x (bandwidth-bound; the bucket index is then the next step).
  already tried: None found; no threading in the file, no BACKLOG item (doc/34 looked at extraction only); the lower_bound → bucket pattern is the sibling's win, not tried in ODIA's MS1 path.

[F06] Instrument before coding memory: per-stage rusage/VmRSS/VmHWM/reservedPoints() lines plus a phase-stamped /proc sampler — then the conditional one-liners (malloc_trim at the pass boundary; size-class-insensitive BlockPool) (target several, votes 3, effort hours)
  mechanism: Every number above the stage walls is reconstructed by subtraction: stage timers are wall-only (ChromatogramExtractor.h:778-791), the ~165 GiB remainder (RSS − live − MS1) is measured only by difference and has three owners with three different fixes. Add getrusage deltas (user, sys, minflt), VmRSS/VmHWM and task count at every existing timer boundary (ChromatogramExtractor.cpp:942, :1105, :1251, Ms1Traces.cpp:129, finish() at OpenDIAlyzer.cpp:4734), relabel the 'scored on the fly' bracket (:4743, which encloses MS1 build + extraction + finish), and print `blocks.reservedPoints()` (:247, tracked at :232, no caller since a0989fc) beside peakPoints() (:1334). Zero-rebuild alternative tonight: a 10-s sampler of /proc/<pid>/status + smaps_rollup joined to log timestamps on any running arm. Conditional follow-ups it decides: (b) VmHWM set after 'extracted ... pass 1' → the extract()-local BlockPool (:740) died as pinned brk-heap free chunks (27.5 KB blocks, below mmap threshold; MALLOC_MMAP_THRESHOLD_ blocked by vm.max_map_count 1,048,576 vs 16.6M blocks) → one malloc_trim(0) at the extract→finish boundary; (a) reserved−peak ≥15% → free_ is keyed by exact n=valid×cycles (:221) so start-truncated blocks (lo=0, ~5,800 sizes) never match interior requests (valid×584) and owned_ never shrinks (:226-233, :241-244) → allocate valid×full_window per precursor, zero-fill only n (stride unchanged, output bit-identical).
  evidence: ChromatogramExtractor.cpp:232,247,740,941,1334; doc/34:80-86 (retention named, unproven; 3.26e9 minor faults); BACKLOG.md:705,1444-1476 (phase instrumentation planned 2026-08-07, never run); doc/10-performance-report.md:66-69 ('Per-stage CPU% is derived'); c3_cent_gate_s43.log:232 (system 11,150 s with no stage owner), :51-52 (6,878 cells at peak vs 6,459 mean = the 6.5% overshoot from a MEAN-cells cap, :689); MEASURED per-group cost p1_c6 vs p1_a1 +8,884 MiB / +21,258,638 groups = 437 B; export arms +58 GB (own_ctl 715,750) — HWM moves with post-extraction state; small-library floor 2.89 GiB (.h:678-682); hv remainder 59.6 GiB vs dn 164.7 GiB (log subtraction); sibling: trim at phase boundaries −30% RSS / +8.5% wall (DIAspeXtract CHANGELOG:841-851).
  expected win: No bytes directly; converts ~150 GiB of unexplained RSS and every INFERRED CPU figure into named terms at the cost of one 23-min fixture run; decides between a 1-line fix worth 8-160 GiB (trim) and a pool redesign worth −100..−150 GiB (INFERRED). Prevents coding the wrong 165 GiB.
  risk: Print-only, none to output. Trap from the sibling record: nested brackets double-count — keep stages disjoint and check the sum against the footer CPU. Trim: +faults when pass 2 re-touches pages (sibling +26% minor faults). Pool: a wrong stride or zero-fill length corrupts traces — the fixture cmp catches it. Re-opens allocator behaviour ONLY at the ≥400 GiB scale that doc/11 never measured.
  verify: Fixture (bench_d2t240n_ih1 configuration) once: stage CPU-s sum within 5% of footer CPU; sink and MS1 lines show CPU/wall in [0.95,1.10]; TSV byte-identical to the uninstrumented binary; FAIL = any stage believed serial shows CPU/wall >1.1. Attach the sampler to the next already-planned full arm (clock-stamped). Pre-registered classification: (a) HWM inside pass-1 extraction AND (reserved−peak)/peak ≥0.15 (≥30 G points ≈ 110 GiB) → pool over-reservation → build the size-class fix: PASS = reserved−peak <2% of peak, Peak ≤540,000 MiB at the c3 configuration, TSV cmp clean; (b) HWM between 'extracted ... in' and 'scored on the fly' → heap non-return → one arm with malloc_trim(0): PASS = Peak lower by ≥8,000 MiB (target ≥40,000), wall within 2%, TSV cmp clean; FAIL = Peak within 1,000 MiB (peak is inside pass 1, close); (c) HWM inside extraction with reserved ≈ peak (<5 G points) → scorer/bookkeeping: size PeakGroup × 22.6M before touching anything. The hv/dn contrast (60 vs 165 GiB) must be reproduced by the printed term or the hypothesis fails.
  already tried: reservedPoints() never printed since a0989fc; tcmalloc heap profile only at 10 GiB (doc/11:84-110); doc/34 asked for a heap profile 'before anyone acts' — none taken; /usr/bin/time -v footers exist only on the two bench arms; allocator tunables (tcmalloc 10.49 vs 10.34 GiB, MALLOC_ARENA_MAX none, MALLOC_TRIM_THRESHOLD_ −11% at +16% time) all at the ≤10 GiB decode-floor scale (doc/10:254-257, doc/11:27-33); pool size-class fix named strongest candidate in doc/34:80-86, never attempted; BACKLOG.md:1791 reuses the shape only for MS1.

[F07] Make Gate C's calibration null deterministic and RT-stratified (RT-defined/strided decoy sample instead of the first n_needed by arrival) — precondition for F02 on armed libraries and for F04 (target several, votes 2, effort days)
  mechanism: NullCalibrationBody::admit (PeakGroupScorer.cpp:1248-1282, verified) pushes decoys in EMISSION order until n_needed, sorts, tau = quantile; emission order depends on chunk boundaries (budget) and would depend on thread interleaving in a parallel sink. On dn_pred_cam the first 20,000 decoy windows are 64-77% pre-gradient, statistic 0 → tau=0 → everything admitted (4.31M vs 0.5M precursors reach the picker, 3x sink time). On human_v2 the same sample sits on the 50%-zeros borderline, so a 0.7 s RT-refit change re-set tau and rejected 103k more targets — the whole leaf-cap 'replicate FAIL' and every human_v2 admission swing. Define the calibration set by RT (decoys whose window ends before the RT at which the n_needed-th decoy window ends, computed from assignments before extraction) or by stride, so tau is a function of library and run only.
  evidence: MEASURED: odia-fixture-cannot-arm-gate-c.md:21-34 (tau=0, 4.31M vs 0.5M, 'a strided/rolling or fixed-RT-band null is a correctness fix; v1.17's rt_min alone is not it'); odia-leaf-cap-b2.md:151-158 (cap replicate FAIL = 103k targets rejected through tau; common-set +1.14%); odia-width-default-regression.md:348-356 (human_v2 width gain is 100% Gate C admission); rep_human2.log:61,79 vs r5_p1.log:61,79 (pass-1 groups 2,689,293 vs 2,689,292; gate C 3600552/3607561 vs 3600551/3607558 at budgets 350 vs 400; provenance caveat run_r5.sh:15-19); c3:61,197 (inert '0/0' on dn_pred_cam); doc/39:374-382 (order dependence recorded); bd73e65 msg lines 30-32 (`-gate_calibration_rt_min` exists); BACKLOG.md:456-457 (thread-invariance harness 'Not yet').
  expected win: ids: removes a ±10% admission swing on gated libraries (MEASURED on human_v2), making cap/centring replicates decidable; cpu: dn_pred_cam sink −50..−66% of 18,841 s if an armed gate rejects as on human_v2 (INFERRED from '3x sink time'); makes every memory lever (F02/F03) output-identical across chunk counts on every library and unblocks F04. ID effect of arming on dn_pred_cam UNKNOWN (read at N(e)).
  risk: Arming the gate on dn_pred_cam removes candidates and could cost IDs (human_v2 gained +18.9% admission by loosening at width 1.0) — read at matched entrapment budget with two baselines; changes the training population and re-lands the fold partition; must be byte-identical when disabled. Gate C RELAXATION (prominence) is refuted on the fixture (−85%) — this re-opens only the null's SAMPLE, not its threshold, and is the opposite direction.
  verify: (a) dn_pred_cam: 'gate C null armed' prints tau>0, zeros<50%, calibration RT inside the gradient; ID clause PASS = region median N(e) over FDP 2-10% within the lineage band (≥ −1%, no ordinal < −2%) vs b2_cap1 and w1_ctl70 AND pass-1 sink ≤0.5 × 14,637 s; FAIL = region < −2%. (b) human_v2: two arms at -live_memory_gb 400 and 200 with the fix → identical 'gate C null armed' tau lines, pass-1 TSV cmp clean, 400-arm IDs within ±1% of r5_p1; re-run hv_cap1 vs hv_cap0 → scored targets differ <1% (was −9.4%) and the cap's common-set +1.14% equals its whole-arm read within 0.5%. (c) determinism: md5-identical TSV at -threads 8/64 and with the flag off vs bd73e65. FAIL = tau differs between budgets or thread counts.
  already tried: v1.17 `-gate_calibration_rt_min` built (excludes pre-gradient decoys), ID effect never read; `-gate_log_path` exists; kimi flagged the order dependence in the doc/39 review and it was parked; tau never measured twice; strided/fixed-RT-band null named as the correctness fix 2026-09-07, never built.

[F08] Thread-scaling curve + thread-invariance harness on the fixture for THIS binary; then pack 16-thread arms per node; match-loop granularity only if the curve says the pool scales (target cpu_wall, votes 2, effort hours)
  mechanism: 77% of wall is serial, so per-arm wall barely depends on thread count while idle cores are wasted; the match region's 17-21% efficiency is arithmetic, not a measurement, and its 2.1x spread across byte-identical arms (2,937-6,292 s) says it is bandwidth/contention-bound. Run the 23-min fixture at -threads 1/4/16/64 with bd73e65 (+/usr/bin/time -v): speedup(64)/speedup(16) tells whether the pool is barrier/imbalance-bound (flat above 16) or bandwidth-bound (flat above ~8). The same four runs are the md5 thread-invariance gate every later fix (F04, F05) must pass. If match(16) ≤1.5x match(64), pack 3-4 arms at -threads 16 per node (needs F02's footprint). Conditional on 'pool scales': MATCH_BATCH 128 → 256/512 (:848; 2 spectra per thread per barrier today, each batch costs max not mean of co-packed frames), then sub-spectrum work units or driver/worker overlap.
  evidence: ChromatogramExtractor.cpp:848 (MATCH_BATCH 128), :836-847 (window-granularity rationale), :1139-1250 (fetch_add per spectrum, barrier per batch), :764-768; ChromatogramExtractor.h:688-693 ('3.6x of the available 7.4x' on a 1,200-precursor library — the only measured efficiency figure); BACKLOG.md:1468-1470 (16→96 threads 1.9x, old decode-dominated binary), :1473 (curve planned, never run), :456-457 (1/8/64 md5 harness 'Not yet'), :2040-2042 (md5-identical 16 vs 64 on lib_targets, old binary); stage table from 44 arms: sink 13,520-15,060 s invariant, match 2,937-6,292 s; cohort_plain at -threads 48 matched 3,768 s (inside the 64-thread spread); odia-timing-noise-floor.md:41-46 (three concurrent arms ran pass 1 no slower); run_c3.sh:114 (-threads 64 always).
  expected win: No direct per-arm wall win; per-node throughput 3-4x with zero code if packing passes; decides whether match-loop work (INFERRED 7,090 → 3,500-4,700 s if efficiency 20% → 30-50%) is worth days and sizes F04's ceiling (the sink will run on the same pool and memory system).
  risk: The fixture understates the full-scale serial share (45% vs 71% of pass-1 extraction) and its match is 172 s, so per-stage noise (14%) requires an idle node and two replicates per point; shapes transport, magnitudes do not. Concurrent arms perturb each other (14%); memory-bandwidth saturation could make 4 arms slower than 2. Match-loop (b)/(c) touch the 'no allocation while a worker reads live' invariant (:1116-1118) and emission resets live[slot] (:960) that workers read (:1205-1207).
  verify: (a) Zero cost: `top -H -p <pid>` during pass 1 of any running arm → one task ~100%, 63 sleeping. (b) Fixture at -threads 1/4/16/64, two replicates, idle node: all four TSVs md5-identical is a precondition for believing any timing (a byte difference is a correctness finding that stops the programme); pre-registered read: match(16)/match(64) ≥1.5 → pool scales, full-scale 20% is load/NUMA (go to F09) and pursue MATCH_BATCH 256 (PASS = match −15% with md5-identical TSV and peak live +≤5%; FAIL <5% → bandwidth-bound, stop); ratio <1.2 → bandwidth/imbalance-bound, skip match-loop work. (c) Two concurrent arms at -threads 16 -live_memory_gb 100: PASS if each wall ≤1.25x the solo wall AND sorted TSVs identical; >1.5x → pack 2 not 4.
  already tried: Thread counts 48/64/96 used but never as a controlled variable (run_2x2.sh:11, run_cohort_plain.sh:20); 'cores busy 6%' computed, never decomposed into serial fraction vs match efficiency; no curve on the bd73e65 extractor; no batch-size or work-unit sweep at full scale.

[F09] Attribute the 3-5 h of system time by thread and stage (page faults vs NUMA balancing), zero code; then one less-eager-trim glibc arm and one NUMA-bound arm (target cpu_wall, votes 2, effort hours)
  mechanism: System time is 8.5-18% of all CPU (11,150-19,134 s) and scales with RSS, not context switches; 3.26e9 faults / 1.4 TiB ≈ 5-11 faults per resident page is not first-touch. Candidates: (i) glibc trimming the main-arena top after each of the sink's transient vectors on the driver thread and re-faulting millions of times → `MALLOC_TOP_PAD_=1073741824` / `MALLOC_TRIM_THRESHOLD_=4294967295` (LESS eager trimming — the opposite direction of the record's negative results); (ii) automatic NUMA balancing on dax's multi-socket board unmapping/re-faulting a 600 GiB RSS → `numactl --interleave=all` or a 16-thread arm bound to one socket (--cpunodebind --membind), combinable with F08's packing; (iii) THP splitting. If the driver thread owns the faults, system time is SERIAL wall (up to 32% of c3's wall); if first-touch zeroing (BlockPool zero-fill :234), only a smaller live set (F02/F03) helps.
  evidence: MEASURED: c3_cent_gate_s43.log:232 (system 3:05:50), ora_ctl.log:191 (5:18:54 of 29:55 CPU); bench_d2t240n_ih1.log:152-165 (sys 455 s, 117.6M faults, 3.58M voluntary switches at 40 GB) vs bench_astralw269_astral.log:157-170 (sys 1,000 s, 255M faults, 2.85M switches at 195 GB): faults scale with RSS, switches do not; doc/34:84-87 (3.26e9 faults 'consistent with continuous fresh-page allocation', unproven); ibminode05: single-socket, numa_balancing=0, THP=madvise (measured); dax unverified (ssh refused from ibminode05); ldd: glibc only, no alternative allocator (CMakeLists.txt). Sibling: DIAspeXtract found EAGER trims and MALLOC_ARENA_MAX harmful, tcmalloc_minimal −7% RSS / halved system time (different program, scale-dependent).
  expected win: Unknown until attributed; bounded above by the system time itself: up to 2-4 h CPU per run and some match wall if NUMA/trim-bound, 0 if spread over workers as first-touch. The measurement costs nothing; the env-var arm costs one fixture run.
  risk: Explicitly RE-OPENING allocator/page-fault behaviour only for the ≥400 GiB regime and only as a measurement: doc/11:27-33 closed tcmalloc/arena for the 10 GiB decode floor and doc/10:254-257 tested MALLOC_TRIM_THRESHOLD_=131072 (the glibc default, i.e. a null test for trimming) and MALLOC_ARENA_MAX=1. A top pad raises RSS by the pad; socket binding halves cores and bandwidth for that arm (only sensible with F08 packing and ≤500 GB footprint).
  verify: Step 1 (zero cost, next running arm, pass-1 sink and match): per-thread minflt/stime from /proc/<pid>/task/*/stat every 10 s; `perf stat -e minor-faults,major-faults,context-switches -p <pid> -- sleep 300`; `grep numa_faults /proc/<pid>/sched`; read /proc/sys/kernel/numa_balancing and THP state on dax. Pre-registered: driver thread ≥60% of minflt and stime → proceed; fault rate × ~1.5 µs within 2x of the system-time share → faults are the sink; numa_hint_faults growing >1e6/s → balancing. Step 2: fixture with/without MALLOC_TOP_PAD_/MALLOC_TRIM_THRESHOLD_, interleaved, two replicates: PASS = system time −50% (455 → ≤230 s), TSV md5-identical, peak RSS +≤2 GB; FAIL = within 14% noise. Step 3 (dax only, if hint faults ≥1e9): numactl-bound 16-thread arm vs unbound same-night control: PASS system time ≤50% of control AND match wall not worse AND identical TSV.
  already tried: doc/10:254-257 (MALLOC_TRIM_THRESHOLD_=131072 −11% RSS +16% time; MALLOC_ARENA_MAX=1 no effect) at 10 GiB; doc/11:27-33 closed allocator for the decode floor; doc/34:84-87 noted 3.26e9 faults, no follow-up; no fault or NUMA attribution at full scale.

[F10] Run the DIALibTune run-specific fine-tuned RT+CCS dn_pred_cam library through ODIA — one arm, single-factor mass-regime gates first (target ids, votes 2, effort hours)
  mechanism: |library 1/K0 residual| >0.030 is the one clean exogenous dose-response (2,967 precursors missed 80% at matched abundance AND matched DIA-NN evidence); the production slice is ±0.025 around a column with median error 0.0235; ODIA's runtime mobility calibration harvests anchors through the window it corrects and recovers <40% of the scale error; the aperture is settled at 0.025 (both directions refuted). A per-peptide fine-tuned CCS reaches DIA-NN's own run-refit (held-out sd 0.01475 vs stock 0.01796 vs DIA-NN 0.01480) plus fine-tuned RT (0.3035 vs DIA-NN refit 0.3524), i.e. centring moves from the slice edge to the centre at the source without touching aperture or calibration code. Remaining gap is z3 (0.02379 vs 0.02155); z2 already tied.
  evidence: MEASURED: FORENSICS_VERDICT.md:34-35,186-187,288-306 (dose table: |res|>0.030 miss 0.813/0.808/0.801, n=2,967); odia-mobility-calibration-is-truncated.md:12-27,218-233; odia-matched-coords-delete-aperture.md:14-16; odia-pick-losses.md:18-34,130-146 (oracle 1/K0 coordinates +1,100..+1,300 net at matched budget, 85% in the >0.025 mis-centred stratum; centring lifts 71% of pick losses); odia-refined-library-diann-arms.md:76-91,129-149,179-205 (DIA-NN arm E: fine-tuned library +7.8/+8.1% at matched e; fine-tuned CCS never run through ODIA); dialibtune-finetune-results.md:70-79 (H100 27.8 s); odia-im-affine-correction-wins.md:54-111 and odia-mass-calibration-reads-library-im.md:9-18 (why an IM edit is multi-factor); odia-charge3-is-a-free-control.md; odia-im-width-sweep-keep-025.md.
  expected win: INFERRED: bounded above by the oracle-centring prize (+1,100..+1,300 net at matched budget, ~+6%, which used DIA-NN's observed 1/K0 on half the ids); realistic +300..+1,000, concentrated in z3 and the mis-centred stratum, plus whatever the RT column adds to pass-2 placement. Not the 46-point gap; the first ID lever with an exogenous mechanism and the cheapest (library swap, zero code).
  risk: HIGH coupling: the fragment-mass probe gates on the LIBRARY 1/K0 at ±0.010 (MassCalibration.cpp:998-1001) — the affine edit moved 79% of precursors out of the probe and cost −34.7% (mechanism = mass gate, loss all z2); the runtime IM fit collapses when anchors carry zero residual (f6_cent); interacts with -im_center_from_pass1 and the anchor floor. Circularity: models tuned on this run's DIA-NN IDs bake run-specific RT/IM into the library searched on the same run; IH1→IH2 transfer is the deployable test and was never measured. Mason-Schamp constant is NOT the lever (do not touch).
  verify: c3 command verbatim except -tr = the fine-tuned library (one variable), vs same-binary stock-library control AND w1_ctl70/b2_cap1. Single-factor gates FIRST, pre-registered: fragment mass centre within max(1 SE, 0.1 ppm) (hard FAIL at >1 ppm from c3's −0.176 ppm) and the pass-1 'fragment mass accuracy' line intact; per-charge anchors on the same side of the 120 floor; IM gate PASSED; otherwise pin -fragment_ppm/-fragment_ppm_offset to the control's fitted values (-mass_calibration off) on BOTH arms and re-read. Report the ion-mobility calibration block (anchors/charge, offset, %MSE) beside it. PASS = region median N(e) over FDP 2-10% ≥ +2.0%, ≥90% ordinals positive, z3 stratum positive at matched depth (charge-3 free control), q≤0.01 ≥0.9x control, DScore histogram not regime-flipped; replicate on a second seed. FAIL = region <0 or mass centre outside tolerance (then it is a mass-regime arm). Second arm: models tuned on ODIA's own pass-1 anchors only (no DIA-NN labels) must reproduce ≥50% of the gain, else circular.
  already tried: Affine library pre-correction −34.7% (mass-gate mechanism); `-im_prior` per-charge affine −5.81% region; `-im_center_from_pass1` +3.94% vs own control but common-depth −77 vs b2_cap1 → NO; oracle centring prize measured on fixture and one full run; the fine-tuned library was run only in DIA-NN (+7.8-8.1% at matched e, per-run caveat); fine-tuned CCS never through ODIA.

[F11] FRAGVEC: at-depth oracle (offline, hours, on the existing x1_fragvec_r1 export) FIRST; wire the fragment-evidence block into sub_scores (six R1 scalars, then 72) only if the oracle clears 72% (target ids, votes 2, effort hours)
  mechanism: The 42-column oracle plateaus at 65.2% (invariant across six perturbations): the information is not in the scalar columns. The buried class fails on fragment co-elution at the correct apex while the co-elution family does no work at the decision point (matched-score AUC 0.811, family ≤0.61). Per-fragment evidence vectors against a data-driven best-fragment reference (DIA-NN's 73-score mechanic) lift matched-score AUC +0.083 but only +298/+232 concordant at common depth. A DIA-NN-labelled, family-grouped OOF GBT on 42+6 / 42+72 columns at depth 37,334 tells whether the extracted per-fragment traces CONTAIN the missing evidence (then scoring is the path) or not (then extraction geometry is, and scorer work is capped near 65%). The columns are already computed on the final pass (PeakGroupScorer.cpp:2488-2572, appended :2653-2654) and only exported (`options.fragvec = !calibrating_ && !out_fragvec_.empty()`, OpenDIAlyzer.cpp:4431-4434; writer :4864-4887); wiring = enum-append to sub_scores as the P1 port did (doc/81:20-25).
  evidence: MEASURED: odia-feature-plateau-65pct.md:8-15,34; odia-fragvec-passes-six-columns.md:48-67,103-122 (full run x1_fragvec_r1: dAUC +0.0825/+0.0842 vs +0.020 floor; +298/+232 concordant at N(62); 49/43 entrapments at comparator depth vs 62; six-column share 0.540 vs 0.60 bar UNDECIDED; no at-depth oracle reported); odia-matched-score-is-the-test.md:137-143,159-166; FORENSICS_VERDICT.md:19-25,59-70; odia-samelib3-ranking-is-the-gap.md (62-64% at matched depth); odia-p1-port-shipped.md:194-203,219-226 (precedent: 21 appended columns = +7.6/+8.8% region; null-append control pattern); x1_fragvec_r1 export exists (analysis77/x1_fragvec_r1.log:194; run_x1.sh:219); diann-id-alias-unimod4.md (canonical join).
  expected win: Oracle: decides the whole ID programme for a python job. Wiring (if funded): MEASURED offline +232..+298 concordant at common depth (ranking-only, the first channel not 97% depth allowance) plus depth allowance; INFERRED engine +1..+3% region N(e); raises the feature ceiling ~0.08 AUC. Does not close the 46% gap.
  risk: An oracle is an upper bound under reference labels; a PASS does not transfer to the semi-supervised scorer (regime flips, decoy null) — a FAIL is the stronger result and is what the plateau predicts. Wiring: native GBT regime flip on append (P1 de-saturated it; 6-72 new columns can re-saturate — read histograms by signature, not label); with DIA-NN-visible ids removed every arm sits at 0.505-0.520, so the gain may be DIA-NN-overlap only; cohort pools were 97.5% CAM-depleted (full run escapes it).
  verify: Stage 0 on x1_fragvec_r1's 21.85M rows: A = 42 columns (must reproduce 64-65% at depth 37,334, canonical UniMod:4 join, role-matched negatives as wf_v50), B = 42+6, C = 42+72; two fold partitions, MDE ~1 pt. Pre-registered: C ≥72% → wire (B ≥70% decides six columns); C ≤68% → evidence is not in the extracted traces, stop scorer work, go to extraction geometry (F10). Stage 1 (engine): arms F6 and F72 vs TWO baselines b2_cap1 and w1_ctl70 on one binary plus a 6-null-column append control. PASS = region median N(e) over FDP 2-10% ≥ +2.0% vs b2_cap1 AND ≥90% ordinals positive AND common-depth concordant ≥ +200 at N(62) with entrapments at comparator depth ≤55 AND q≤0.01 ≥0.95 × 18,072 AND null-append control within ±0.5%. FAIL = region <0, common-depth < +100, or q≤0.01 <0.75 × 18,072 (collapse fails regardless). Report three numbers plus a fixed-model re-score of the shared 42 columns. A PASS earns the human_v2 replicate, not a default.
  already tried: Export validated at full scale (x1_fragvec_r1 scored TSV sha-identical to control; engine block == sealed builder); offline acceptance PAIR passed (wf_v67/v69); at-depth oracle recall on FRAGVEC not in the record; 42-column oracle done twice (doc/77, wf_v50); never wired into the native scorer; MS1-iso (+2.4 buried, cohort) next in queue with no full-run read.

[F12] Shrink Ms1Traces: drop never-read rows (uncovered/outside precursors) risk-free, then store per-precursor [lo,hi) spans instead of the dense precursors×bins matrix (target memory, votes 2, effort days)
  mechanism: Ms1Traces is a dense 9.92M × 1,343 f32 matrix (values_.assign(rows*bins), Ms1Traces.cpp:82; dense write :162; Ms1Traces.h:107) built before pass 1 and alive through both passes: 50,835 MiB. Rows for the 1,275,674 precursors covered by no isolation window and the 24,576 predicted outside the run get no blocks → no candidate → no ms1_coelution call → 6.6 GiB of zeros never read; the `keep` mask path already exists for the isotope export (Ms1Traces.h:92-100). Each read precursor's trace is only read inside its pass-1 window (±400 s = 584 of 1,343 bins) → a per-row span with zero outside stores 43% (the DIAspeXtract 'entry frame + span arena' structure); pass-2 windows are re-centred by the fitted map and must stay inside the stored span. Becomes the largest single structure (~22%) once F02/F03/F06 land.
  evidence: MEASURED sizes: c3_cent_gate_s43.log:46 (50,835.2 MiB, 5,589.6 s), :54-55 (uncovered/outside counts); OpenDIAlyzer.cpp:1860-1882 (MS1 checked against the budget, not charged); ChromatogramExtractor.cpp:535-540 (uncovered precursors get no assignment); sibling: DIAspeXtract trace redesign 120 B + 12 B/pt → 20 B + 4 B/frame with identical digest; doc/11:131-161 (compaction order); odia-cohort-wave-results.md (the MS1 channel is exhausted as a FEATURE — this changes storage, not the feature).
  expected win: INFERRED: −6.6 GiB risk-free (row drop); −28..−45 GiB more with spans (49.6 → ~15-20 GiB). Secondary: the build loop touches fewer rows (small F05 synergy). Combined with F02/F03, sub-100 GB arms become plausible.
  risk: Row drop: none (never-read rows). Spans: a pass-2 candidate whose re-centred window leaves the span reads zeros → ms1_coelution changes → small but nonzero ID effect, must be read at matched entrapment budget. Zero-sentinel/presence semantics in consumers ('every consumer skips zeros' was false in the sibling project). Memory-only; no wall win.
  verify: Row drop: `cmp` TSV clean vs c3 (any diff is a bug), MS1 MiB line lower by ~6,600. Spans: PASS = ≥99.9% of ms1_coelution values identical to the dense arm AND pass-2 N(e) median over FDP 2-10% within the ±1% lineage band AND Peak −28 GiB ±3 AND MS1 MiB ≤20,000; FAIL = N(e) median < −1% or q≤0.01 collapse.
  already tried: None for the search path; the `keep` mask was built for -out_ms1_iso only; compaction order scoped in doc/11, only the fusion step shipped; MS1 sized and skipped when over budget (OpenDIAlyzer.cpp:1868-1882) but never charged or shrunk.

[F13] Offline stage 0 (A2 lane-join + A1 replay on the 16,349 multi-candidate class) to decide whether to build retained-group re-seeking at refined coordinates (DIA-NN's reselection) (target ids, votes 1, effort weeks)
  mechanism: DIA-NN re-seeks the winning peak group every iteration with refitted weights over several margin-kept candidates and refits RT/mass/mobility ~12 times against retained in-RAM data (DIA-NN-workflow-handoff.md:445-452,464-573). ODIA picks ≤3 candidates once per pass (`-max_candidates` 3, OpenDIAlyzer.cpp:1131; coelution picker PeakGroupScorer.cpp:1549-1585), the semi-supervised loop only re-ranks those rows (lda.h:1142-1166), `refitsChangeScores()` returns false (PeakGroupScorer.cpp:2721-2727), pass 2 re-extracts from the raw file. Multi-candidate precursors are missed 85%; the rt_dn oracle shows +13.5 buried points locked behind candidate selection. Retaining groups (~120 B each) makes re-seeking cheap; re-extract only candidates whose coordinates moved. This is the only lever whose promise is DIA-NN-overlap recall (P1 port and FRAGVEC are overlap-flat/complementary). Stage 0 is hours and offline; the build is weeks and competes with F10/F11 for the same weeks — fund it only on the stage-0 read.
  evidence: FORENSICS_VERDICT.md:35-36 (H2b unique +0.0239, 'mediator, not cause'), :248-274 (A1 replay and A2 lane-join named, NOT run); odia-tensor-verdict.md:62-75 (oracle 93.9/87.4 vs deployable 87.4/73.9); odia-tensor-parked-by-guards.md (P0 multi-candidate 88.1/75.4); doc/81:57-62; doc/45:60-64,334-341; BACKLOG.md:201-205 ('retain peak groups between passes ... never built'); odia-training-recipe-is-a-null.md:116-133 (`-max_candidates 6` −0.4..−0.6%, 'they rescore, they do not discover'); odia-why-we-miss.md (fragment co-elution failure; 1/K0 residual + multi-candidate are the handles).
  expected win: Upper bound cohort-projected +13.5 buried pts (oracle, non-deployable); realistic unknown — stage 0 decides at the cost of a python job on existing tables. INFERRED.
  risk: The lane-join may show ONE mechanism (weak signal → several plausible peaks → incoherent fragments → bad rank), in which case re-seeking buys nothing; more candidates without competition is a measured null; any change in the pass-2 candidate set re-lands the fold partition (compare on bytes or with a frozen model); large engine change.
  verify: Stage 0 (hours, offline, tables under /scratch/kohlbach/odia2x2/wf_v44/ on dax and analysis77/forensics): A2 lane-join (MS2 statistics × var_cand_count × C-near/far × 1/K0 residual) and A1 replay on the 16,349 multi-candidate class. Pre-registered: ≥40% of buried multi-candidate precursors had the correct candidate PROPOSED within 8.31 s but outranked → fund the build; <20% → do not build, record closed; between → build only retained-groups rescoring (no re-seek). Stage 1 (engine): vs b2_cap1 and w1_ctl70, region median N(e) ≥ +2.0% with ≥90% ordinals, common-depth ≥ +200 at N(62), q≤0.01 not collapsed, AND DIA-NN-overlap recall at depth 37,334 ≥65.1% (64.1% + 1.0) — overlap-flat = FAIL for this lever specifically.
  already tried: A1/A2 listed as next measurements 2026-09-01, not run; 'retain peak groups between passes' scoped, never built; `-max_candidates 6/10` null; k=1/12/24 iteration knobs null (re-rank, do not re-seek); multi-hypothesis mobility grid break-even, not built; single-anchor cohort experiments STOPPED (not re-proposed).

[F14] Stabilise the semi-supervised GBT's training population: DIA-NN-style decoy-tail truncation before training, read on ≥3 fold seeds (target ids, votes 1, effort days)
  mechanism: 10.9M decoy groups against ~17k positives (1:155) give pure-positive leaves a Newton step ~1/p; one fold's rounds-2-5 oscillation (leaves 214/127) collapses q≤0.01 16,965 → 4,105 at seed 43 with identical evidence. DIA-NN's `remove_rubbish` drops q>0.5 targets and truncates the decoy score tail at the minimum surviving target score before its NN sees 153k PSMs (handoff:588-592; dn_ih1_cam.log:39); ODIA trains on all 22M rows (lda.h:1142-1166). Truncating decoys below the lowest positive-candidate score leaves q untouched (q counts decoys ABOVE target scores) and changes only the training negatives — a clean single-factor arm, unlike the leaf cap (also moves pass 1) and intercept 0 (loses to the cap). Also the one place where IDs could move through the decoy POPULATION (the error-estimate recalibration itself changes nothing at matched budget — MEASURED, odia-decoy-null-borrowed-evidence.md:63).
  evidence: MEASURED: odia-scorer-regime-flips.md:278-285 (seed collapse); odia-leaf-cap-b2.md:97-101,151-158,176 (imbalance, runaway, cap replicate FAIL via gate, 'whether intercept 0 prevents the oscillation is UNTESTED'); SESSION_NOTES.md:5687-5689 (i1_int0: +1.62% vs w1, −1.58% vs b2_cap1, 333/333 negative); DIA-NN-workflow-handoff.md:588-592; lda.h:86-94 (train_fdr 0.15/0.05 selection); odia-gbt-is-8pct-behind-xgboost.md (framework fine, split finding).
  expected win: INFERRED: removes the partition collapse (a shipping defect worth up to −12,860 at q≤0.01 on a bad seed, 0 on a good one); region gain bounded by the cap's +3.2% vs w1; the value is reliability of every other arm's read, not depth.
  risk: Negatives must be role-matched (the 'decoy best rows only' set is a label leak — the learner reads var_cand_rank as the label, odia-feature-plateau-65pct.md:34); may re-saturate the scorer; changes the fold partition; DIA-NN's rule has preconditions (≥5,000 survivors and total >5x survivors). Must not be stacked with the cap or intercept 0 without a 2x2 (cap × intercept interfere, −4.5%). Entrapments must never enter training.
  verify: Seeds 42/43/44 × {b2_cap1 configuration, truncation-only (no cap)}. PASS = every seed's q≤0.01 ≥0.9 × the cross-seed median AND region median vs b2_cap1 ≥ −1.0% with no seed below −2% AND max |leaf| in any fold's first 5 trees <20 (read from the .model file, odia-reading-the-scorer.md) AND the human_v2 pair (after F07) ≥0 on the common candidate set. FAIL = any seed collapses (>25% below the cross-seed median) or region < −2%. Fixed-model re-score beside the native read on every arm.
  already tried: Leaf cap 1.0/2.0 (+3.2/+4.0% region vs w1; replicate FAIL on human_v2 through Gate C; not a default); intercept 0 (loses to the cap, 333/333); lr 0.05 FAIL; 240 rounds a seed draw; k=1/12/24 null; `-gbt_warmup_rounds` built, unread; decoy-tail truncation never tried; XGBoost import closed (NOT REPRODUCED).


## 3. PASTED SOURCE (commit bd73e65; `cat -n` line numbers; the extractor/reader/MS1 files are byte-identical between bd73e65 and HEAD)

Note: the four extractor/reader files (ChromatogramExtractor.cpp/.h, MzPeakSource.cpp, Ms1Traces.cpp/.h) are unchanged since bd73e65; OpenDIAlyzer.cpp, PeakGroupScorer.cpp/.h, lda.h, gbt.h are shown AT bd73e65 (`git show bd73e65:<path>`), so all line numbers below are bd73e65 line numbers and match the diagnoses.

### src/extract/ChromatogramExtractor.cpp

```text  (lines 81-145)
    81	    /// Transitions of one window, ordered by m/z and bucketed on log m/z.
    82	    ///
    83	    /// Bucketing on the LOG of m/z is what makes the lookup O(1): a ppm
    84	    /// tolerance is a constant width in log space, so one bucket covers one
    85	    /// tolerance at every mass, where a linear grid would be far too coarse at
    86	    /// 200 Th and far too fine at 1800.
    87	    ///
    88	    /// The point of this index is to invert the match. Probing the spectrum
    89	    /// once per transition is O(N log M) with N transitions and M peaks, and N
    90	    /// exceeds M by three orders of magnitude at proteome scale -- 4.6 M
    91	    /// against ~1,100. Asking instead which transitions each PEAK could belong
    92	    /// to is O(M + matches).
    93	    struct MzIndex
    94	    {
    95	      std::vector<double> mz;              ///< ascending
    96	      /// Which precursor's live block this transition writes into, and which
    97	      /// row of it.
    98	      ///
    99	      /// This replaces the absolute offset into one flat point array that the
   100	      /// index used to carry, and the replacement is the whole point: an
   101	      /// absolute offset only exists if the whole array exists. A slot is
   102	      /// resolved through `LiveSlot`, which holds a pointer that is null until
   103	      /// the pass reaches the precursor's window and null again once it has
   104	      /// left -- so the index survives the block being allocated and freed
   105	      /// underneath it, and does not have to be rebuilt when it is.
   106	      std::vector<std::uint32_t> slot;
   107	      /// A precursor has 12 transitions here and the library format allows far
   108	      /// fewer than 65,535, so the row is 16 bits. At 51 M transitions the two
   109	      /// bytes saved against a uint32 are 102 MB.
   110	      std::vector<std::uint16_t> row;
   111	      /// Expected 1/K0 of each transition's precursor. NaN when the library has
   112	      /// none, which disables the per-precursor mobility test for it rather
   113	      /// than rejecting it -- absent information is not evidence of mismatch.
   114	      std::vector<float> precursor_im;
   115	      std::vector<std::uint32_t> bucket;        ///< bucket -> first entry
   116	      double log_base = 0.0, log_lo = 0.0;
   117	
   118	      std::size_t bucketOf(double m) const
   119	      {
   120	        const auto b = static_cast<std::ptrdiff_t>((std::log(m) - log_lo) / log_base);
   121	        return std::size_t(std::clamp<std::ptrdiff_t>(b, 0, std::ptrdiff_t(bucket.size()) - 2));
   122	      }
   123	
   124	      void build(double tolerance_ppm)
   125	      {
   126	        if (mz.empty()) { return; }
   127	        log_base = std::log1p(tolerance_ppm * 1e-6);
   128	        if (!(log_base > 0.0)) { log_base = 1e-6; }
   129	        log_lo = std::log(mz.front());
   130	        const auto n = static_cast<std::size_t>(
   131	                         (std::log(mz.back()) - log_lo) / log_base) + 2;
   132	        bucket.assign(n + 1, static_cast<std::uint32_t>(mz.size()));
   133	        // Filled backwards so each bucket points at its FIRST entry and empty
   134	        // buckets inherit the next non-empty one, which keeps the lookup
   135	        // branchless.
   136	        for (std::size_t i = mz.size(); i-- > 0;)
   137	        {
   138	          bucket[bucketOf(mz[i])] = static_cast<std::uint32_t>(i);
   139	        }
   140	        for (std::size_t b = bucket.size() - 1; b-- > 0;)
   141	        {
   142	          bucket[b] = std::min(bucket[b], bucket[b + 1]);
   143	        }
   144	      }
   145	    };
```

```text  (lines 160-260)
   160	    /// Where a precursor's points are RIGHT NOW.
   161	    ///
   162	    /// `base` is null outside the precursor's retention-time window, which is
   163	    /// most of the run for most precursors and is the entire memory argument.
   164	    /// The cycle range is kept here too so the match loop resolves a peak in
   165	    /// one cache line rather than chasing back into the assignment table.
   166	    struct LiveSlot
   167	    {
   168	      float* base = nullptr;
   169	      /// Intensity-weighted m/z deviation per cell: sum(intensity * ppm) in
   170	      /// `ppm_num` and sum(intensity) in `ppm_den`, reduced where it is read.
   171	      ///
   172	      /// TWO planes, after a wrong shortcut. The first version stored a single
   173	      /// plane holding "the deviation of the peak that won the cell", justified
   174	      /// as exact under Aggregate::Max -- but **Max is not the default**;
   175	      /// Options::aggregate defaults to Sum and so does the CLI. Under Sum every
   176	      /// peak "wins", so that plane held the LAST matching peak in spectrum
   177	      /// order, which is arbitrary: peak arrays are not sorted, so it was
   178	      /// neither the brightest nor the nearest. Found by external review.
   179	      ///
   180	      /// The weighted form is correct under both modes and is the right
   181	      /// quantity regardless: a bright fragment's centroid is better determined
   182	      /// than a dim one's, and on DIA data the dim end is where interference
   183	      /// lives.
   184	      ///
   185	      /// Race-free for the same reason `base` is: a thread owns a spectrum, and
   186	      /// a spectrum owns a distinct cycle, so no two threads touch one cell.
   187	      ///
   188	      /// Null unless Options::collect_mass_residuals. Live blocks are the
   189	      /// extractor's dominant memory term, so this is opt-in.
   190	      float* ppm_num = nullptr;
   191	      float* ppm_den = nullptr;
   192	      /// The same construction for ION MOBILITY. Per-(row, cycle) for the same
   193	      /// reason: the match loop is threaded over spectra and is race-free only
   194	      /// because each spectrum owns a distinct cycle.
   195	      float* im_num = nullptr;
   196	      float* im_den = nullptr;
   197	      std::uint32_t lo = 0, hi = 0;
   198	    };
   199	
   200	    /// Fixed-size blocks, reused rather than returned to the allocator.
   201	    ///
   202	    /// Nearly every precursor asks for the same size -- 12 transitions by the
   203	    /// window's cycle count -- so a free list per size hits essentially always,
   204	    /// and the pass does one allocation per concurrently-live precursor for the
   205	    /// whole run instead of one per precursor. That matters at 4.26 M
   206	    /// precursors: 4.26 M allocate/free pairs of ~40 KB through malloc is where
   207	    /// a general allocator fragments, and fragmentation would put back exactly
   208	    /// the memory this change exists to remove.
   209	    ///
   210	    /// Blocks are handed out ZEROED, because the match accumulates into them
   211	    /// and a reused block still holds the previous precursor's peaks.
   212	    class BlockPool
   213	    {
   214	    public:
   215	      float* take(std::size_t n)
   216	      {
   217	        // A precursor all of whose transitions lack a product m/z has no rows.
   218	        // It is still emitted -- a consumer counting precursors that yielded
   219	        // nothing must see it -- so it needs a base that is not null.
   220	        if (n == 0) { return &dummy_; }
   221	        auto& bin = free_[n];
   222	        float* p = nullptr;
   223	        if (!bin.empty())
   224	        {
   225	          p = bin.back();
   226	          bin.pop_back();
   227	        }
   228	        else
   229	        {
   230	          owned_.push_back(std::make_unique_for_overwrite<float[]>(n));
   231	          p = owned_.back().get();
   232	          reserved_ += n;
   233	        }
   234	        std::fill_n(p, n, 0.0f);
   235	        live_ += n;
   236	        peak_ = std::max(peak_, live_);
   237	        return p;
   238	      }
   239	
   240	      void give(float* p, std::size_t n)
   241	      {
   242	        if (n == 0) { return; }
   243	        free_[n].push_back(p);
   244	        live_ -= n;
   245	      }
   246	
   247	      std::uint64_t reservedPoints() const { return reserved_; }
   248	      std::uint64_t peakPoints() const { return peak_; }
   249	
   250	    private:
   251	      std::unordered_map<std::size_t, std::vector<float*>> free_;
   252	      std::vector<std::unique_ptr<float[]>> owned_;
   253	      std::uint64_t live_ = 0, peak_ = 0, reserved_ = 0;
   254	      float dummy_ = 0.0f;
   255	    };
   256	  } // namespace
   257	
   258	  void ChromatogramCollector::begin(const ChromatogramLayout& layout)
   259	  {
   260	    // Refused with the number rather than left to a bad_alloc from an
```

```text  (lines 528-545)
   528	      for (std::size_t w = 0; w < windows.size(); ++w)
   529	      {
   530	        if (!windows[w].contains(mz) || axis[w].rt.empty()) { continue; }
   531	        ++covering;
   532	        const double offset_from_centre = std::abs(mz - windows[w].centre());
   533	        if (offset_from_centre < best_offset) { best_offset = offset_from_centre; best = w; }
   534	      }
   535	      if (covering == 0)
   536	      {
   537	        ++st.precursors_without_window;
   538	        if (options.terminal_reason)
   539	        { options.terminal_reason[i] = Options::kNoWindowCoverage; }
   540	        continue;
   541	      }
   542	      if (covering > 1) { ++st.precursors_in_several_windows; }
   543	
   544	      const std::size_t w = best;
   545	      std::size_t lo = global_lo[w], hi = global_hi[w];
```

```text  (lines 660-965)
   660	        else
   661	        {
   662	          --live;
   663	          live_points -= slot_points[by_end[j]];
   664	          ++j;
   665	        }
   666	      }
   667	    }
   668	
   669	    // The library, split into chunks whose live sets each fit under the cap.
   670	    //
   671	    // Greedy over start time with a running set of end times, so a chunk is a
   672	    // band of the gradient rather than an arbitrary slice of the library -- and
   673	    // a band needs only the spectra inside it, which is why chunking costs a
   674	    // fraction of a decode pass rather than a whole one per chunk.
   675	    std::vector<std::vector<std::uint32_t>> chunks;
   676	    // Derive the cap from the memory budget when one is given. Stated in
   677	    // bytes because that is the quantity a caller actually has, and inverted
   678	    // here because only the extractor knows the mean cells per precursor.
   679	    std::size_t cap = options.max_live_precursors;
   680	    if (options.live_memory_budget_bytes > 0)
   681	    {
   682	      std::size_t planes = 1;
   683	      if (options.collect_mass_residuals) { planes += 2; }
   684	      if (options.collect_im_residuals)   { planes += 2; }
   685	      double mean_cells = 0.0;
   686	      for (const Assignment& a : assignments)
   687	      { mean_cells += double(a.valid) * double(a.hi - a.lo); }
   688	      if (!assignments.empty()) { mean_cells /= double(assignments.size()); }
   689	      const double per = mean_cells * 4.0 * double(planes);
   690	      const std::size_t derived = per > 0.0
   691	        ? std::max<std::size_t>(1, std::size_t(double(options.live_memory_budget_bytes) / per))
   692	        : 0;
   693	      // The tighter of the two wins: an explicit -max_live_precursors is a
   694	      // caller's assertion and must not be loosened by a budget.
   695	      cap = (cap == 0) ? derived : std::min(cap, derived);
   696	      st.live_budget_note = "budget " +
   697	        std::to_string(options.live_memory_budget_bytes / (1024ull*1024*1024)) +
   698	        " GiB / " + std::to_string(std::size_t(per)) + " B per live precursor (" +
   699	        std::to_string(planes) + " planes) -> cap " + std::to_string(cap);
   700	    }
   701	    if (cap == 0 || overlap_precursors <= cap)
   702	    {
   703	      chunks.push_back(std::move(by_start));
   704	      st.memory_bound_by = "retention-time overlap (" +
   705	                           std::to_string(overlap_precursors) + " of " +
   706	                           std::to_string(n_slots) + " precursors live at once)";
   707	    }
   708	    else
   709	    {
   710	      using MinHeap = std::priority_queue<float, std::vector<float>, std::greater<float>>;
   711	      MinHeap ends;
   712	      std::vector<std::uint32_t> current;
   713	      for (const std::uint32_t s : by_start)
   714	      {
   715	        while (!ends.empty() && ends.top() < slot_rt_lo[s]) { ends.pop(); }
   716	        if (ends.size() + 1 > cap && !current.empty())
   717	        {
   718	          chunks.push_back(std::move(current));
   719	          current.clear();
   720	          ends = MinHeap{};
   721	        }
   722	        current.push_back(s);
   723	        ends.push(slot_rt_hi[s]);
   724	      }
   725	      if (!current.empty()) { chunks.push_back(std::move(current)); }
   726	      st.memory_bound_by = "precursor cap (" + std::to_string(cap) + "), " +
   727	                           std::to_string(chunks.size()) + " chunks; " +
   728	                           std::to_string(overlap_precursors) +
   729	                           " would have been live at once";
   730	    }
   731	    st.chunks = chunks.size();
   732	    if (chunks.empty()) { chunks.emplace_back(); st.chunks = 1; }
   733	
   734	    st.index_seconds = std::chrono::duration<double>(
   735	                         std::chrono::steady_clock::now() - t_index).count();
   736	
   737	    // Where every precursor's points are right now: null outside its window,
   738	    // which for most precursors is most of the run.
   739	    std::vector<LiveSlot> live(n_slots);
   740	    BlockPool blocks;
   741	    std::size_t live_now = 0, live_peak = 0;
   742	    std::atomic<std::size_t> nonzero{0};
   743	    // Points a peak fell on while the destination precursor was not live. It
   744	    // must be zero -- the cycle test below already excludes them -- and it is
   745	    // counted rather than assumed, because the alternative to catching it is a
   746	    // silently truncated chromatogram.
   747	    std::atomic<std::size_t> unhoused{0};
   748	
   749	    // One forward pass over the run in acquisition order, so the decode stays
   750	    // contiguous -- a window's spectra are strided through the file, and
   751	    // fetching them window by window would ask the reader for a stride it
   752	    // cannot serve cheaply.
   753	    //
   754	    // Each spectrum is matched by INVERTING the loop: for every peak, the m/z
   755	    // index says which transitions could contain it. The alternative -- probe
   756	    // the spectrum once per transition -- is O(N log M) with N transitions and
   757	    // M peaks, and at proteome scale N is 4.6 M against M of ~1,100.
   758	    //
   759	    // Matching runs in parallel with no synchronisation at all: a
   760	    // (transition, cycle) pair has exactly one destination index, computed
   761	    // rather than allocated, so no two threads ever write the same slot.
   762	    const unsigned threads = options.threads ? options.threads
   763	                                             : std::max(1u, std::thread::hardware_concurrency());
   764	    std::vector<SpectrumPeaks> block;
   765	    // Workers are created ONCE and parked between blocks. Spawning them per
   766	    // block cost more than the matching did: at 64 workers and a 128-spectrum
   767	    // block each got two spectra, and parallel ran slower than serial.
   768	    // Raising the block size hid that; creating them once removes it.
   769	    struct Pool
   770	    {
   771	      std::vector<std::thread> workers;
   772	      std::mutex m;
   773	      std::condition_variable cv, done_cv;
   774	      std::function<void()> job;
   775	      std::size_t generation = 0, finished = 0;
   776	      bool stop = false;
   777	
   778	      void start(unsigned n)
   779	      {
   780	        for (unsigned i = 0; i < n; ++i)
   781	        {
   782	          workers.emplace_back([this] {
   783	            std::size_t seen = 0;
   784	            for (;;)
   785	            {
   786	              std::function<void()> mine;
   787	              {
   788	                std::unique_lock<std::mutex> lock(m);
   789	                cv.wait(lock, [&] { return stop || generation != seen; });
   790	                if (stop) { return; }
   791	                seen = generation;
   792	                mine = job;
   793	              }
   794	              mine();
   795	              {
   796	                std::lock_guard<std::mutex> lock(m);
   797	                ++finished;
   798	              }
   799	              done_cv.notify_one();
   800	            }
   801	          });
   802	        }
   803	      }
   804	
   805	      void run(std::function<void()> f, unsigned n)
   806	      {
   807	        {
   808	          std::lock_guard<std::mutex> lock(m);
   809	          job = std::move(f);
   810	          finished = 0;
   811	          ++generation;
   812	        }
   813	        cv.notify_all();
   814	        std::unique_lock<std::mutex> lock(m);
   815	        done_cv.wait(lock, [&] { return finished >= n; });
   816	      }
   817	
   818	      ~Pool()
   819	      {
   820	        {
   821	          std::lock_guard<std::mutex> lock(m);
   822	          stop = true;
   823	        }
   824	        cv.notify_all();
   825	        for (auto& w : workers) { if (w.joinable()) { w.join(); } }
   826	      }
   827	    } pool_impl;
   828	    if (threads > 1) { pool_impl.start(threads); }
   829	    // Large enough that spawning workers is amortised. At 128 spectra and 64
   830	    // threads each worker got two spectra and thread creation cost more than
   831	    // the matching did -- measured 0.19 s against 0.14 s single-threaded.
   832	    // Was a constant; now `Options::decode_block`, because a heap profile put
   833	    // 5.57 GiB of a 9.45 GiB live peak in the block's `SpectrumPeaks` copies.
   834	    // Zero keeps the historical default rather than degenerating to no block.
   835	    const std::size_t BLOCK = options.decode_block ? options.decode_block : 256;
   836	    // Matching walks the decoded block in smaller batches, because the
   837	    // allocate/free cursors can only move between batches -- no chromatogram
   838	    // may appear or vanish while a worker is reading `live`. The batch is
   839	    // therefore the GRANULARITY OF THE SLIDING WINDOW: a precursor goes live up
   840	    // to one batch early and is freed up to one batch late.
   841	    //
   842	    // 128 spectra is ~5 cycles across IH1's 24 windows, about 7 s of gradient
   843	    // against retention-time windows of hundreds -- so the rounding is under 1%
   844	    // of what is resident. At the 1,024 the decode uses it would be 59 s, which
   845	    // is 5% of a 1,200 s window and 100% of a short one. The cost of the finer
   846	    // batch is one more pool wakeup per 128 spectra, microseconds against the
   847	    // ~1 s the decode of those spectra takes.
   848	    constexpr std::size_t MATCH_BATCH = 128;
   849	
   850	    std::size_t first_spectrum = 0, last_spectrum = info.size();
   851	    if (options.rt_high > options.rt_low)
   852	    {
   853	      first_spectrum = static_cast<std::size_t>(
   854	        std::lower_bound(info.begin(), info.end(), options.rt_low,
   855	                         [](const SpectrumInfo& s, double v) { return s.retention_time < v; })
   856	        - info.begin());
   857	      last_spectrum = static_cast<std::size_t>(
   858	        std::lower_bound(info.begin(), info.end(), options.rt_high,
   859	                         [](const SpectrumInfo& s, double v) { return s.retention_time < v; })
   860	        - info.begin());
   861	    }
   862	    st.spectra_read = last_spectrum - first_spectrum;
   863	
   864	    // A precursor leaving the pass: its chromatogram is final, so hand it over
   865	    // and give the memory back.
   866	    std::vector<std::uint64_t> off_scratch;
   867	    std::vector<std::uint32_t> count_scratch;
   868	    const auto activate = [&](std::uint32_t slot) {
   869	      const Assignment& a = assignments[slot];
   870	      const std::size_t cells = std::size_t(a.valid) * (a.hi - a.lo);
   871	      live[slot].base = blocks.take(cells);
   872	      if (options.collect_mass_residuals)
   873	      {
   874	        live[slot].ppm_num = blocks.take(cells);
   875	        live[slot].ppm_den = blocks.take(cells);
   876	      }
   877	      if (options.collect_im_residuals)
   878	      {
   879	        live[slot].im_num = blocks.take(cells);
   880	        live[slot].im_den = blocks.take(cells);
   881	      }
   882	      live[slot].lo = a.lo;
   883	      live[slot].hi = a.hi;
   884	      live_peak = std::max(live_peak, ++live_now);
   885	    };
   886	    const auto emit = [&](std::uint32_t slot) {
   887	      const Assignment& a = assignments[slot];
   888	      // A precursor whose window the pass never reached -- possible only when
   889	      // the caller restricted the run -- still owes the sink a trace of zeros,
   890	      // which is what the flat array used to hand it.
   891	      if (live[slot].base == nullptr) { activate(slot); }
   892	
   893	      const std::uint32_t tb = p.transition_begin[a.precursor];
   894	      const std::uint32_t tc = p.transition_count[a.precursor];
   895	      const std::uint32_t cycles = a.hi - a.lo;
   896	      off_scratch.assign(tc, 0);
   897	      count_scratch.assign(tc, 0);
   898	      std::uint32_t row = 0;
   899	      for (std::uint32_t k = 0; k < tc; ++k)
   900	      {
   901	        // Skipped on the same test as the counting loop above, and it has to
   902	        // BE the same test: a transition counted there and skipped here is a
   903	        // point run with no axis behind it.
   904	        if (t.product_mz[tb + k] == MZ_INVALID) { continue; }
   905	        off_scratch[k] = std::uint64_t(row) * cycles;
   906	        count_scratch[k] = cycles;
   907	        ++row;
   908	      }
   909	
   910	      // The representation's invariant, checked rather than assumed: every
   911	      // point of the precursor must land on the axis it names. Breaking it does
   912	      // not crash -- a consumer reads past the end of a vector and returns
   913	      // whatever is there -- so it has to be caught here or not at all.
   914	      if (a.window >= axes.size() ||
   915	          std::size_t(a.lo) + cycles > axes[a.window].size())
   916	      {
   917	        throw std::logic_error(
   918	          "chromatogram precursor " + std::to_string(a.precursor) + " has " +
   919	          std::to_string(cycles) + " points from cycle " + std::to_string(a.lo) +
   920	          " on axis " + std::to_string(a.window) + ", which holds " +
   921	          std::to_string(a.window < axes.size() ? axes[a.window].size() : 0) + " cycles");
   922	      }
   923	
   924	      PrecursorChromatogram trace;
   925	      trace.precursor = a.precursor;
   926	      trace.transition_begin = tb;
   927	      trace.transition_count = tc;
   928	      trace.axis = a.window;
   929	      trace.axis_begin = a.lo;
   930	      trace.cycles = cycles;
   931	      trace.rt = axes[a.window].data() + a.lo;
   932	      trace.points = live[slot].base;
   933	      trace.ppm_num = live[slot].ppm_num;
   934	      trace.ppm_den = live[slot].ppm_den;
   935	      trace.im_num = live[slot].im_num;
   936	      trace.im_den = live[slot].im_den;
   937	      trace.offset = off_scratch.data();
   938	      trace.count = count_scratch.data();
   939	
   940	      const auto t_sink = std::chrono::steady_clock::now();
   941	      sink.accept(trace);
   942	      st.sink_seconds += std::chrono::duration<double>(
   943	                           std::chrono::steady_clock::now() - t_sink).count();
   944	
   945	      blocks.give(live[slot].base, std::size_t(a.valid) * cycles);
   946	      // The residual planes are pool blocks too, and resetting the slot without
   947	      // giving them back leaks two per emitted precursor. Mild on a 2,665-
   948	      // precursor benchmark (3.55 -> 3.87 GiB) because the pool retains them;
   949	      // at the 4.26 M design scale it would be the whole live set again, twice.
   950	      if (live[slot].ppm_num != nullptr)
   951	      {
   952	        blocks.give(live[slot].ppm_num, std::size_t(a.valid) * cycles);
   953	        blocks.give(live[slot].ppm_den, std::size_t(a.valid) * cycles);
   954	      }
   955	      if (live[slot].im_num != nullptr)
   956	      {
   957	        blocks.give(live[slot].im_num, std::size_t(a.valid) * cycles);
   958	        blocks.give(live[slot].im_den, std::size_t(a.valid) * cycles);
   959	      }
   960	      live[slot] = LiveSlot{};
   961	      --live_now;
   962	    };
   963	
   964	    for (const auto& slots : chunks)
   965	    {
```

```text  (lines 1058-1337)
  1058	      // Activation and expiry cursors: this chunk's precursors ordered by the
  1059	      // cycle their window starts at, and by the cycle it ends at.
  1060	      std::vector<std::vector<std::uint32_t>> by_lo(windows.size()), by_hi(windows.size());
  1061	      for (const std::uint32_t slot : slots) { by_lo[assignments[slot].window].push_back(slot); }
  1062	      for (std::size_t w = 0; w < windows.size(); ++w)
  1063	      {
  1064	        std::sort(by_lo[w].begin(), by_lo[w].end(), [&](std::uint32_t a, std::uint32_t b) {
  1065	          return assignments[a].lo < assignments[b].lo; });
  1066	        by_hi[w] = by_lo[w];
  1067	        std::sort(by_hi[w].begin(), by_hi[w].end(), [&](std::uint32_t a, std::uint32_t b) {
  1068	          return assignments[a].hi < assignments[b].hi; });
  1069	      }
  1070	      std::vector<std::size_t> cur_lo(windows.size(), 0), cur_hi(windows.size(), 0);
  1071	      std::vector<std::uint32_t> seen(windows.size(), 0);
  1072	
  1073	      // One chunk needs only the spectra its own precursors elute in. With a
  1074	      // single chunk that is the caller's range unchanged, so nothing about a
  1075	      // run that fits today changes.
  1076	      std::size_t chunk_first = first_spectrum, chunk_last = last_spectrum;
  1077	      if (chunks.size() > 1 && !slots.empty())
  1078	      {
  1079	        double lo_rt = std::numeric_limits<double>::infinity();
  1080	        double hi_rt = -std::numeric_limits<double>::infinity();
  1081	        for (const std::uint32_t slot : slots)
  1082	        {
  1083	          lo_rt = std::min(lo_rt, double(slot_rt_lo[slot]));
  1084	          hi_rt = std::max(hi_rt, double(slot_rt_hi[slot]));
  1085	        }
  1086	        chunk_first = std::max(first_spectrum, static_cast<std::size_t>(
  1087	          std::lower_bound(info.begin(), info.end(), lo_rt,
  1088	                           [](const SpectrumInfo& s, double v) { return s.retention_time < v; })
  1089	          - info.begin()));
  1090	        chunk_last = std::min(last_spectrum, static_cast<std::size_t>(
  1091	          std::upper_bound(info.begin(), info.end(), hi_rt,
  1092	                           [](double v, const SpectrumInfo& s) { return v < s.retention_time; })
  1093	          - info.begin()));
  1094	      }
  1095	
  1096	      st.index_seconds += std::chrono::duration<double>(
  1097	                            std::chrono::steady_clock::now() - t_chunk_index).count();
  1098	
  1099	      for (std::size_t begin = chunk_first; begin < chunk_last; begin += BLOCK)
  1100	      {
  1101	        const std::size_t end = std::min(begin + BLOCK, chunk_last);
  1102	
  1103	        const auto t_decode = std::chrono::steady_clock::now();
  1104	        source.peaks(begin, end, block);
  1105	        st.decode_seconds += std::chrono::duration<double>(
  1106	                               std::chrono::steady_clock::now() - t_decode).count();
  1107	        st.spectra_decoded += end - begin;
  1108	
  1109	        // Decoding is per block, because that is what the reader wants;
  1110	        // matching walks it in smaller batches, because that is the
  1111	        // granularity at which the sliding window can move. See MATCH_BATCH.
  1112	        for (std::size_t batch = begin; batch < end; batch += MATCH_BATCH)
  1113	        {
  1114	          const std::size_t batch_end = std::min(batch + MATCH_BATCH, end);
  1115	
  1116	          // Everything whose window starts inside this batch goes live. Doing it
  1117	          // between batches rather than per spectrum is what keeps the pass free
  1118	          // of synchronisation: no allocation happens while a worker is running.
  1119	          const auto t_alloc = std::chrono::steady_clock::now();
  1120	          for (std::size_t si = batch; si < batch_end; ++si)
  1121	          {
  1122	            const std::uint32_t w = window_of[si];
  1123	            if (w == std::numeric_limits<std::uint32_t>::max()) { continue; }
  1124	            seen[w] = std::max(seen[w], cycle_of[si] + 1);
  1125	          }
  1126	          for (std::size_t w = 0; w < windows.size(); ++w)
  1127	          {
  1128	            while (cur_lo[w] < by_lo[w].size() &&
  1129	                   assignments[by_lo[w][cur_lo[w]]].lo < seen[w])
  1130	            {
  1131	              activate(by_lo[w][cur_lo[w]]);
  1132	              ++cur_lo[w];
  1133	            }
  1134	          }
  1135	          st.assemble_seconds += std::chrono::duration<double>(
  1136	                                   std::chrono::steady_clock::now() - t_alloc).count();
  1137	
  1138	          const auto t_match = std::chrono::steady_clock::now();
  1139	          std::atomic<std::size_t> next{batch};
  1140	          const auto work = [&]() {
  1141	            std::size_t local_nonzero = 0, local_unhoused = 0;
  1142	            for (;;)
  1143	            {
  1144	              const std::size_t si = next.fetch_add(1);
  1145	              if (si >= batch_end) { break; }
  1146	              const std::uint32_t w = window_of[si];
  1147	              if (w == std::numeric_limits<std::uint32_t>::max()) { continue; }
  1148	              const auto& x = index[w];
  1149	              if (x.mz.empty()) { continue; }
  1150	              const auto& peaks = block[si - begin];
  1151	              const std::uint32_t c = cycle_of[si];
  1152	              // The run's mobility and the FRAME BAND test are separate questions.
  1153	              // Reading the peak's mobility only inside the band's branch made
  1154	              // `use_ion_mobility` switch off the per-precursor test as well: the
  1155	              // mobility stayed NaN, and a NaN mobility skips that test under the
  1156	              // "absent information is not evidence of mismatch" rule meant for a
  1157	              // run that carries no mobility at all. -no_ion_mobility is the
  1158	              // control arm for measuring what mobility filtering buys, so it has
  1159	              // to turn off exactly the one filter it names.
  1160	              const bool has_im = peaks.hasIonMobility();
  1161	              const bool use_band = options.use_ion_mobility && has_im;
  1162	              const double im_low = info[si].window.im_low, im_high = info[si].window.im_high;
  1163	
  1164	              for (std::size_t k = 0; k < peaks.size(); ++k)
  1165	              {
  1166	                const double m = peaks.mz[k];
  1167	                // The bounds must allow for the tolerance. A peak just BELOW the
  1168	                // smallest transition can still be within ppm of it, and skipping
  1169	                // it silently returned zero for the lowest-m/z transition of every
  1170	                // window -- which a brute-force scan caught and nothing else would
  1171	                // have, because a chromatogram of zeros looks like an ion that is
  1172	                // simply not there.
  1173	                const double slack = m * options.fragment_ppm * 1e-6 * 1.01 + 1e-6;
  1174	                if (m + slack < x.mz.front() || m - slack > x.mz.back()) { continue; }
  1175	                // The frame's band separates the co-packed windows. It cannot
  1176	                // separate a precursor from its same-window neighbours, which is
  1177	                // what the per-transition test below does.
  1178	                const double peak_im = has_im ? double(peaks.ion_mobility[k])
  1179	                                              : std::numeric_limits<double>::quiet_NaN();
  1180	                // Half-open, [im_low, im_high). Co-packed windows are separated by
  1181	                // a DERIVED band -- the split is the midpoint between their two
  1182	                // mobility positions -- so adjacent bands share their boundary
  1183	                // exactly, and a peak sitting on it would enter both windows and be
  1184	                // integrated twice under Sum. Measure zero on real data, and the
  1185	                // reason it is stated as a convention rather than left to chance.
  1186	                if (use_band && (peak_im < im_low || peak_im >= im_high)) { continue; }
  1187	                // The tolerance belongs to the TRANSITION, not to the peak:
  1188	                // a match means |peak - transition| <= transition * ppm. Searching
  1189	                // by peak makes it tempting to size the window on the peak
  1190	                // instead, which differs at the boundary and silently includes or
  1191	                // drops edge matches -- 11 of 240 transitions disagreed with a
  1192	                // brute-force scan because of exactly that.
  1193	                //
  1194	                // So the SEARCH window is deliberately a little wide, and the
  1195	                // exact test is applied per candidate inside it.
  1196	                const double ppm = options.fragment_ppm * 1e-6;
  1197	                const double im_half = options.precursor_im_window;
  1198	                const bool sum_peaks = options.aggregate == Options::Aggregate::Sum;
  1199	                std::size_t i = x.bucket[x.bucketOf(std::max(m - slack, x.mz.front()))];
  1200	                const float intensity = peaks.intensity[k];
  1201	                for (; i < x.mz.size() && x.mz[i] <= m + slack; ++i)
  1202	                {
  1203	                  if (std::abs(m - x.mz[i]) > x.mz[i] * ppm) { continue; }
  1204	                  // Per-precursor mobility. Skipped when either side is unknown:
  1205	                  // absent information is not evidence of mismatch.
  1206	                  if (im_half > 0.0 && !std::isnan(peak_im))
  1207	                  {
  1208	                    const float want = x.precursor_im[i];
  1209	                    if (!std::isnan(want) && std::abs(peak_im - want) > im_half) { continue; }
  1210	                  }
  1211	                  const LiveSlot& s = live[x.slot[i]];
  1212	                  if (c < s.lo || c >= s.hi) { continue; }
  1213	                  if (s.base == nullptr) { ++local_unhoused; continue; }
  1214	                  float& at = s.base[std::size_t(x.row[i]) * (s.hi - s.lo) + (c - s.lo)];
  1215	                  // Maximum, not sum: two peaks inside one tolerance are the same
  1216	                  // ion split by centroiding far more often than they are two ions.
  1217	                  if (at == 0.0f && intensity > 0.0f) { ++local_nonzero; }
  1218	                  if (sum_peaks) { at += intensity; }
  1219	                  else if (intensity > at) { at = intensity; }
  1220	                  // The deviation was already computed to test the match above;
  1221	                  // it has been discarded here since the extractor was written.
  1222	                  // Accumulated for EVERY contributing peak, not just a winner:
  1223	                  // under Sum aggregation there is no winner, and picking one
  1224	                  // by arrival order is picking at random.
  1225	                  if (s.ppm_num != nullptr && intensity > 0.0f)
  1226	                  {
  1227	                    const std::size_t at_i =
  1228	                      std::size_t(x.row[i]) * (s.hi - s.lo) + (c - s.lo);
  1229	                    s.ppm_num[at_i] +=
  1230	                      intensity * static_cast<float>((m - x.mz[i]) / x.mz[i] * 1e6);
  1231	                    s.ppm_den[at_i] += intensity;
  1232	                  }
  1233	                  // The OBSERVED 1/K0 of whatever produced this peak. NaN
  1234	                  // mobility is skipped rather than accumulated as zero: a
  1235	                  // missing measurement is not a mobility of nothing.
  1236	                  if (s.im_num != nullptr && intensity > 0.0f && !std::isnan(peak_im))
  1237	                  {
  1238	                    const std::size_t at_i =
  1239	                      std::size_t(x.row[i]) * (s.hi - s.lo) + (c - s.lo);
  1240	                    s.im_num[at_i] += intensity * static_cast<float>(peak_im);
  1241	                    s.im_den[at_i] += intensity;
  1242	                  }
  1243	                }
  1244	              }
  1245	            }
  1246	            nonzero.fetch_add(local_nonzero);
  1247	            unhoused.fetch_add(local_unhoused);
  1248	          };
  1249	          if (threads <= 1) { work(); }
  1250	          else { pool_impl.run(work, threads); }
  1251	          st.match_seconds += std::chrono::duration<double>(
  1252	                                std::chrono::steady_clock::now() - t_match).count();
  1253	
  1254	          // Everything the pass has now passed the end of is FINAL. This is the
  1255	          // whole change: the chromatogram goes to the consumer and the memory
  1256	          // goes back, so what is resident is what is live at one retention time
  1257	          // rather than what the library contains.
  1258	          const auto t_free = std::chrono::steady_clock::now();
  1259	          const double sink_before = st.sink_seconds;
  1260	          for (std::size_t w = 0; w < windows.size(); ++w)
  1261	          {
  1262	            while (cur_hi[w] < by_hi[w].size() &&
  1263	                   assignments[by_hi[w][cur_hi[w]]].hi <= seen[w])
  1264	            {
  1265	              emit(by_hi[w][cur_hi[w]]);
  1266	              ++cur_hi[w];
  1267	            }
  1268	          }
  1269	          // The sink's own time is reported separately, so it is taken out here
  1270	          // rather than counted twice.
  1271	          st.assemble_seconds += std::chrono::duration<double>(
  1272	                                   std::chrono::steady_clock::now() - t_free).count()
  1273	                                 - (st.sink_seconds - sink_before);
  1274	        }
  1275	
  1276	        if (options.progress_every && (begin / BLOCK) % 16 == 0)
  1277	        {
  1278	          std::cerr << "\r  " << st.spectra_decoded << " / " << st.spectra_read
  1279	                    << " spectra"
  1280	                    << std::flush;
  1281	        }
  1282	      }
  1283	
  1284	      // Whatever the pass ended inside is final too.
  1285	      const auto t_flush = std::chrono::steady_clock::now();
  1286	      const double sink_before_flush = st.sink_seconds;
  1287	      for (std::size_t w = 0; w < windows.size(); ++w)
  1288	      {
  1289	        while (cur_hi[w] < by_hi[w].size())
  1290	        {
  1291	          emit(by_hi[w][cur_hi[w]]);
  1292	          ++cur_hi[w];
  1293	        }
  1294	      }
  1295	      st.assemble_seconds += std::chrono::duration<double>(
  1296	                               std::chrono::steady_clock::now() - t_flush).count()
  1297	                             - (st.sink_seconds - sink_before_flush);
  1298	    }
  1299	    if (options.progress_every) { std::cerr << "\r" << std::string(48, ' ') << "\r"; }
  1300	
  1301	    if (unhoused.load() != 0)
  1302	    {
  1303	      throw std::logic_error(
  1304	        std::to_string(unhoused.load()) +
  1305	        " matched points had no live chromatogram to go into. The sliding "
  1306	        "window released a precursor before the pass had left its retention-"
  1307	        "time range, which silently truncates its trace.");
  1308	    }
  1309	
  1310	    // Precursors nothing extracted still owe the sink a trace, empty, so that
  1311	    // a consumer counting them sees what it saw when every precursor had a row.
  1312	    const auto t_empty = std::chrono::steady_clock::now();
  1313	    for (std::size_t i = 0; i < n_prec; ++i)
  1314	    {
  1315	      if (p.transition_count[i] == 0 || assigned[i]) { continue; }
  1316	      PrecursorChromatogram trace;
  1317	      trace.precursor = static_cast<std::uint32_t>(i);
  1318	      trace.transition_begin = p.transition_begin[i];
  1319	      trace.transition_count = p.transition_count[i];
  1320	      sink.accept(trace);
  1321	    }
  1322	    st.sink_seconds += std::chrono::duration<double>(
  1323	                         std::chrono::steady_clock::now() - t_empty).count();
  1324	
  1325	    st.precursors = n_prec;
  1326	    st.transitions = n_trans;
  1327	    st.points = total_points;
  1328	    st.nonzero_points = nonzero.load();
  1329	    st.mean_live_transitions = n_trans ? double(live_sum) / double(n_trans) : 0.0;
  1330	    // MEASURED, from the pool, not from the interval sweep: activation moves
  1331	    // only between match batches, so a precursor goes live up to one batch
  1332	    // early, and a chunked run never reaches the sweep's number at all.
  1333	    st.peak_live_precursors = live_peak;
  1334	    st.peak_live_points = blocks.peakPoints();
  1335	  }
  1336	
  1337	} // namespace ODIA
```

### include/odia/ChromatogramExtractor.h

```text  (lines 600-712)
   600	      /// at every other step. Max over that returns the interference envelope,
   601	      /// and -- being an order statistic -- acts as a hard threshold rather
   602	      /// than a graded penalty: a peptide below the envelope contributes
   603	      /// nothing at all. That is the shape of the measured failure, where the
   604	      /// bottom four abundance deciles sat exactly at the decoy null.
   605	      ///
   606	      /// Sum is only correct once the box is small. Summing over the old
   607	      /// oversized box makes the trace worse, not better.
   608	      enum class Aggregate { Sum, Max };
   609	      Aggregate aggregate = Aggregate::Sum;
   610	
   611	      /// Extract only the first N precursors of the library, 0 for all.
   612	      std::size_t max_precursors = 0;
   613	
   614	      /// Cap on how many precursors may have their chromatograms live at one
   615	      /// instant. 0 lets the data decide.
   616	      ///
   617	      /// The sliding window bounds memory by RETENTION-TIME OVERLAP, which is
   618	      /// a property of the run and the window width rather than of the library
   619	      /// -- and on a wide window that bound is weak. With a 600 s half-width on
   620	      /// a 1,859 s gradient a precursor is live for ~1,200 s, so about 65% of
   621	      /// the library is live at once and the win over holding all of it is
   622	      /// 1.5x. At a calibrated 120 s it is ~13% and the win is 7x.
   623	      ///
   624	      /// This is the second mechanism, for when the first is not enough. Above
   625	      /// the cap the library is split into chunks whose live sets each fit, and
   626	      /// each chunk is a separate pass -- restricted to the spectra its own
   627	      /// precursors need, so the extra decode is the chunks' retention-time
   628	      /// spans and not a whole run per chunk. It is stated in precursors rather
   629	      /// than bytes because that is the number a caller can reason about
   630	      /// against a library size; `Stats::peak_live_points` reports what it
   631	      /// cost.
   632	      ///
   633	      /// A pass is expensive -- decode is ~600 s on IH1 and is 65% of Phase 2
   634	      /// -- so chunking is a fallback, not a default. Which mechanism actually
   635	      /// bound the memory is reported in `Stats::memory_bound_by`.
   636	      std::size_t max_live_precursors = 0;
   637	
   638	      /// Memory budget for the live blocks, BYTES. Non-zero overrides
   639	      /// `max_live_precursors`, which is derived from it.
   640	      ///
   641	      /// `max_live_precursors = 0` means "bounded only by retention-time
   642	      /// overlap", and that bound is a property of the RUN, not of the library:
   643	      /// it does not tighten as the library grows. On a 4,986,319-precursor
   644	      /// library with no iRT map -- so every precursor predicted across the
   645	      /// whole gradient -- essentially the entire library is live at once, and
   646	      /// the extractor was OOM-killed at 588 GB on a 995 GB node.
   647	      ///
   648	      /// A cap stated in precursors cannot be chosen without knowing the
   649	      /// transition count and window width, which is why nobody set one. A cap
   650	      /// stated in BYTES can: the extractor knows the mean cells per precursor
   651	      /// and how many planes each carries, so it inverts the budget into a
   652	      /// precursor count itself and reports what it chose.
   653	      ///
   654	      /// Per live precursor the cost is
   655	      ///     valid_transitions x cycles_in_window x 4 B x planes
   656	      /// where planes is 1, +2 with `collect_mass_residuals`, +2 with
   657	      /// `collect_im_residuals` -- so 5 with both, which is the default.
   658	      std::size_t live_memory_budget_bytes = 0;
   659	
   660	      /// How many spectra are decoded and held at once.
   661	      ///
   662	      /// This is the largest single term in the run's memory, and it was found
   663	      /// by profile rather than by reading: a tcmalloc heap profile of a
   664	      /// 1,200-precursor IH1 run put **5.57 GiB of a 9.45 GiB live peak** in the
   665	      /// three `assign` calls of `MzPeakSource::peaks` -- 2.79 / 1.39 / 1.39 GiB
   666	      /// across `mz` (double), `intensity` and `ion_mobility` (float), exactly
   667	      /// the 8:4:4 ratio of their element sizes.
   668	      ///
   669	      /// The cost is linear in this number and it is NOT the sliding window:
   670	      /// six other explanations were measured and refuted first (per-thread
   671	      /// parquet buffers, page cache, glibc fragmentation, arena count, sparse
   672	      /// chromatograms, and mzPeak's own row-group cache, which is bounded).
   673	      /// See doc/11-memory-and-compaction-plan.md.
   674	      ///
   675	      /// Measured on IH1 with a 1,200-precursor library at 4 threads, where the
   676	      /// only variable was this number:
   677	      ///
   678	      /// | block | peak RSS | wall |
   679	      /// |------:|---------:|-----:|
   680	      /// | 1024  | 10.34 GiB | 7:00.3 |
   681	      /// |  256  |  2.89 GiB | 6:53.7 |
   682	      /// |   64  |  1.40 GiB | 6:50.3 |
   683	      ///
   684	      /// There is no tradeoff to balance: smaller is both smaller and faster,
   685	      /// so 1024 was pure waste. The saving is better than linear because
   686	      /// Arrow's decode buffers scale with the batch as well as our peak arrays.
   687	      ///
   688	      /// 256 rather than 64 because of a DIFFERENT measurement, already in this
   689	      /// file: at 128 spectra and 64 threads each worker got two spectra and
   690	      /// thread creation cost more than the matching (0.19 s against 0.14 s
   691	      /// single-threaded). The arms above ran at 4 threads and so cannot see
   692	      /// that. 256 keeps four spectra per worker at 64 threads and still takes
   693	      /// 3.6x of the available 7.4x. Lower it when threads are few.
   694	      std::size_t decode_block = 256;
   695	
   696	      /// Extract only every Nth precursor. 1 is all of them.
   697	      ///
   698	      /// For pass 1 of the two-pass workflow, which extracts solely to harvest
   699	      /// RT anchors. Indices are NOT renumbered -- a strided pass still reports
   700	      /// full-library precursor indices, so the anchor harvest and the iRT it
   701	      /// is fitted against cannot drift apart.
   702	      std::size_t precursor_stride = 1;
   703	
   704	      /// Which residue class the stride keeps, so equal-sized subsets with
   705	      /// DIFFERENT members can be compared. Exists to separate "how many
   706	      /// anchors" from "which anchors" when measuring how stable the RT fit is.
   707	      std::size_t precursor_offset = 0;
   708	
   709	      /// Optional per-precursor keep-mask, in FULL-library indexing, or null
   710	      /// for all of them. Non-null means `PrecursorPrefilter` has already
   711	      /// decided this precursor cannot be supported by the run.
   712	      ///
```

```text  (lines 775-795)
   775	      std::size_t progress_every = 500;
   776	    };
   777	
   778	    struct Stats
   779	    {
   780	      std::size_t spectra_read = 0;
   781	      std::size_t precursors = 0;
   782	      std::size_t transitions = 0;
   783	      std::size_t points = 0;
   784	      std::size_t nonzero_points = 0;
   785	      double decode_seconds = 0.0;
   786	      double match_seconds = 0.0;
   787	      double index_seconds = 0.0;
   788	      /// Allocating, zeroing and releasing the live chromatograms.
   789	      double assemble_seconds = 0.0;
   790	      /// Time inside the sink -- copying, or scoring and discarding.
   791	      double sink_seconds = 0.0;
   792	
   793	      /// Precursors not extracted because their predicted elution fell outside
   794	      /// the run entirely, or could not be predicted at all (a NaN iRT).
   795	      std::size_t outside_rt_range = 0;
```

### src/OpenDIAlyzer.cpp (TOPP tool driver)

```text  (lines 620-660)
   620	    registerDoubleOption_("irt_intercept", "<b>", 0.0, "See -irt_slope.", false, true);
   621	    registerIntOption_("max_precursors", "<n>", 0,
   622	                       "Extract only the first N precursors, 0 for all.", false, true);
   623	    registerDoubleOption_("live_memory_gb", "<GB>", 20.0,
   624	                          "Memory budget for the extractor's live blocks, GiB. The live-"
   625	                          "precursor cap is DERIVED from this, because a cap in precursors "
   626	                          "cannot be chosen without knowing the transition count and window "
   627	                          "width -- which is why it was never set, and why a 4,986,319-"
   628	                          "precursor library OOM-killed the run at 588 GB.\n\n"
   629	                          "0 = auto: take 60% of MemAvailable at startup. A negative value "
   630	                          "restores the old behaviour of no cap at all, bounded only by "
   631	                          "retention-time overlap -- which does NOT tighten as the library "
   632	                          "grows and is therefore not a bound on a large library.\n\n"
   633	                          "-max_live_precursors still applies and the TIGHTER of the two "
   634	                          "wins: an explicit cap is a caller's assertion and a budget must "
   635	                          "not loosen it.",
   636	                          false);
   637	    registerIntOption_("max_live_precursors", "<n>", 0,
   638	                       "Cap how many precursors may have chromatograms in memory at "
   639	                       "once. 0 lets the retention-time overlap decide, which is the "
   640	                       "cheap bound; above the cap the library is split into chunks "
   641	                       "and each is a separate pass over the run, which costs a "
   642	                       "decode. The run reports which of the two bound it.",
   643	                       false, true);
   644	    registerIntOption_("pass1_precursors", "<n>", 0,
   645	                       "Sample about this many precursors for pass 1. It runs only to harvest "
   646	                       "retention-time anchors -- 692 came from 2,450 precursors and "
   647	                       "-min_anchors defaults to 20 -- so extracting the whole "
   648	                       "library costs memory proportional to it: ~274 GiB at 4.26 M "
   649	                       "precursors over the whole run, which is why the calibrated "
   650	                       "pass-2 window does not unblock a large library on its own. "
   651	                       "Indices are not renumbered, so anchors still refer to the "
   652	                       "full library. This is a TARGET COUNT rather than a stride "
   653	                       "because a stride is not scale-invariant: 4 starved a 2,450-"
   654	                       "precursor library (anchors 709 -> 149, identifications 1862 "
   655	                       "-> 635) while at 4.26 M it would still leave a million. 0 "
   656	                       "disables sampling.",
   657	                       false, true);
   658	    registerIntOption_("refine_rounds", "<n>", 0,
   659	                       "After pass 2, alternately refit the retention-time map from "
   660	                       "the current best identifications and refit the discriminant, "
```

```text  (lines 1031-1036)
  1031	    setValidStrings_("classifier", {"xgboost", "lda", "gbt", "nn", "percolator"});
  1032	    registerIntOption_("passes", "<n>", 2,
  1033	                       "1 extracts once with the given calibration. 2 extracts wide, "
  1034	                       "scores, fits the retention-time map from the confident "
  1035	                       "identifications, and re-extracts narrow.", false);
  1036	    registerDoubleOption_("rt_window_pass1", "<seconds>", 0.0,
```

```text  (lines 1131-1135)
  1131	    registerIntOption_("max_candidates", "<n>", 3,
  1132	                       "Candidate peak groups kept per precursor. More than one on "
  1133	                       "purpose: keeping only the best hides the true peak whenever "
  1134	                       "it ranks second, and leaves the decoys nothing to be wrong "
  1135	                       "about, which deflates the FDR.", false, true);
```

```text  (lines 1377-1413)
  1377	    registerFlag_("collect_mass_residuals",
  1378	                  "Keep the m/z deviation of every matched peak and report it per peak group "
  1379	                  "as Mass.Ppm. The deviation is computed anyway to test the match and has "
  1380	                  "always been discarded. This is the only way to measure the run's real "
  1381	                  "fragment mass error: a standalone probe asks whether SOME peak lies within "
  1382	                  "tolerance near a time, and on a mostly-absent library the answer is yes by "
  1383	                  "coincidence -- an RT-shifted control produced a LARGER apparent offset than "
  1384	                  "the true apex. Matches kept here are constrained by co-elution and, at "
  1385	                  "q<=0.01, by the whole discriminant. Costs TWO extra float planes per live "
  1386	                  "block -- an intensity-weighted sum and its denominator, because a single "
  1387	                  "plane could not tell 'no peak matched' from 'matched at exactly 0.000 ppm' "
  1388	                  "under the default Sum aggregation. The decode path, not the live blocks, is "
  1389	                  "what sets this tool's memory floor, so the cost is real but not binding.",
  1390	                  true);
  1391	    registerStringOption_("entrapment_prefix", "<text>", "",
  1392	                          "Protein-name prefix marking ENTRAPMENT precursors in the library: "
  1393	                          "real peptides known to be absent from the sample, so every one "
  1394	                          "reported is a genuine false positive. Given this, the run reports "
  1395	                          "the false discovery PROPORTION beside its own q-values.\n\n"
  1396	                          "This is the check target-decoy cannot perform on itself. Decoys "
  1397	                          "are CONSTRUCTED, so a classifier can learn what construction looks "
  1398	                          "like rather than what a wrong answer looks like; entrapment "
  1399	                          "peptides carry no construction signature. Wen et al. (Nat Methods "
  1400	                          "22:1454, 2025) measure no DIA tool controlling peptide-level FDR, "
  1401	                          "with DIA-NN's true precursor FDP above 2.3% at a nominal 1%, so a "
  1402	                          "reported q-value is not evidence until this has been run.",
  1403	                          false);
  1404	    registerStringOption_("im_features", "<mode>", "auto",
  1405	                          "Measure the run's OBSERVED 1/K0 per candidate from intensity-"
  1406	                          "weighted mobility planes, and feed var_im_delta from it. Until "
  1407	                          "this existed the feature was all-NaN on every run and was dropped "
  1408	                          "as carrying no information -- nothing in the pipeline measured "
  1409	                          "observed mobility per precursor. 'auto' collects wherever the run "
  1410	                          "has ion mobility and is inert where it does not; 'off' restores "
  1411	                          "the previous behaviour.",
  1412	                          false);
  1413	    setValidStrings_("im_features", {"auto", "spread", "off"});
```

```text  (lines 1770-1812)
  1770	                                    ? seed_im_window_
  1771	                                    : getDoubleOption_("precursor_im_window") * pass1_im_scale_;
  1772	    // The width measurement reads the ppm planes, so asking for it turns them
  1773	    // on. Making the user pass two flags that only work together is a way of
  1774	    // producing runs that silently measured nothing.
  1775	    options.collect_mass_residuals = getFlag_("collect_mass_residuals") ||
  1776	                                     getStringOption_("mass_width_from_ids") != "off" ||
  1777	                                     getStringOption_("mass_anchors") != "off";
  1778	    // The mobility planes are what make IM_DELTA computable at all -- it has
  1779	    // been all-NaN on every run since it was added, because nothing measured
  1780	    // observed 1/K0 per precursor. Inert on a run without mobility, so the
  1781	    // default is on: an absent feature costs more than two float planes.
  1782	    options.collect_im_residuals = getStringOption_("im_features") != "off";
  1783	    options.aggregate = getStringOption_("aggregate") == "max"
  1784	                          ? ODIA::ChromatogramExtractor::Options::Aggregate::Max
  1785	                          : ODIA::ChromatogramExtractor::Options::Aggregate::Sum;
  1786	    options.irt_slope = library_rt_is_run_seconds ? 1.0 : getDoubleOption_("irt_slope");
  1787	    options.irt_intercept = library_rt_is_run_seconds ? 0.0 : getDoubleOption_("irt_intercept");
  1788	    // After the mass calibration, because the mobility probe matches fragments
  1789	    // through the mass window and a mis-centred one fills its sample with the
  1790	    // interference that calibration exists to exclude; and after the iRT map,
  1791	    // because the probe uses it to look only where a precursor should elute.
  1792	    applyMobilityCalibration_(library, *source, options);
  1793	    options.threads = static_cast<unsigned>(std::max(1, getIntOption_("threads")));
  1794	    options.max_live_precursors = static_cast<std::size_t>(
  1795	      std::max(0, getIntOption_("max_live_precursors")));
  1796	    {
  1797	      const double gb = getDoubleOption_("live_memory_gb");
  1798	      if (gb < 0.0) { options.live_memory_budget_bytes = 0; }
  1799	      else if (gb > 0.0)
  1800	      { options.live_memory_budget_bytes = std::size_t(gb * 1024.0 * 1024.0 * 1024.0); }
  1801	      else
  1802	      {
  1803	        // auto: 60% of what the kernel says is available right now. Not of
  1804	        // MemTotal -- a node with other tenants has less than it owns, and the
  1805	        // vault records a run measuring 52 effective cores on a 128-core box
  1806	        // for the same reason.
  1807	        std::size_t avail_kb = 0;
  1808	        if (std::ifstream mi("/proc/meminfo"); mi)
  1809	        {
  1810	          std::string k; unsigned long long v; std::string unit;
  1811	          while (mi >> k >> v >> unit)
  1812	          { if (k == "MemAvailable:") { avail_kb = v; break; } }
```

```text  (lines 1838-1925)
  1838	    // its memory behaviour is the one part of this pipeline that is measured.
  1839	    // See Ms1Traces for why co-elution is the only MS1 quantity worth carrying.
  1840	    if (ms1_traces_.empty() && !getFlag_("no_ms1"))
  1841	    {
  1842	      const auto t0 = std::chrono::steady_clock::now();
  1843	      // Width RELATIVE to the MS2 window, because the right ratio is an
  1844	      // empirical question and both of my confident answers were wrong.
  1845	      //
  1846	      // It shipped at 2x on an unmeasured assumption. I then set it to 1x on a
  1847	      // consistency argument -- the two traces being correlated should sample
  1848	      // one ion population -- and that cost 76 identifications on
  1849	      // IH1/lib_targets (1306 -> 1230). The consistency argument is sound about
  1850	      // what the correlation MEANS and wrong about what it is worth; MS1 ions
  1851	      // are not mobility-selected by an isolation window, so the precursor's
  1852	      // MS1 mobility spread is genuinely wider than its fragments'.
  1853	      //
  1854	      // The NaN handling is NOT part of this knob and stays fixed: a precursor
  1855	      // with no library 1/K0 is ungated, matching MS2, rather than having every
  1856	      // peak rejected.
  1857	      // SIZE IT BEFORE ALLOCATING IT. Ms1Traces is a DENSE
  1858	      // precursors x MS1-spectra float matrix (Ms1Traces.cpp:52,
  1859	      // values_.assign(np * bins_, 0.0f)), so it grows linearly with the
  1860	      // library and is charged before extraction begins -- it is not covered by
  1861	      // the live-block budget.
  1862	      //
  1863	      // At 5,330 precursors it is 83 MB and invisible. At 4,986,319 it is
  1864	      // 77.6 GB, which is 13% of the 588 GB peak that OOM-killed a benchmark
  1865	      // run. Every memory figure this project published before that was
  1866	      // measured on a library ~1000x too small to show it.
  1867	      {
  1868	        const std::size_t ms1_bins = source->ms1Spectra().size();
  1869	        const double need = double(library.precursorCount()) * double(ms1_bins) * 4.0;
  1870	        const double cap = double(options.live_memory_budget_bytes);
  1871	        if (cap > 0.0 && need > cap)
  1872	        {
  1873	          std::ostringstream w;
  1874	          w.setf(std::ios::fixed); w.precision(1);
  1875	          w << "MS1 traces would need " << need / 1073741824.0 << " GiB ("
  1876	            << library.precursorCount() << " precursors x " << ms1_bins
  1877	            << " MS1 spectra x 4 B), above the " << cap / 1073741824.0
  1878	            << " GiB budget -- SKIPPING them. var_ms1_coelution will be absent "
  1879	               "rather than the run being killed. Raise -live_memory_gb to keep "
  1880	               "them, or narrow the library.";
  1881	          writeLogWarn_(w.str());
  1882	        }
  1883	        else
  1884	        {
  1885	          double ms1_resid = std::numeric_limits<double>::quiet_NaN();
  1886	          ms1_traces_ = ODIA::Ms1Traces::build(library, *source, options.fragment_ppm,
  1887	                                               options.precursor_im_window *
  1888	                                                 getDoubleOption_("ms1_im_scale"),
  1889	                                               ms1PpmCentre_(), &ms1_resid);
  1890	          {
  1891	            // Reported, not asserted. The MS1 axis borrows the FRAGMENT offset,
  1892	            // which is a hypothesis: the instrument need not err identically on
  1893	            // the two. If the residual does not sit near 0 the borrowed centre
  1894	            // is wrong and this line is how we find out.
  1895	            std::ostringstream m;
  1896	            m.setf(std::ios::fixed); m.precision(3);
  1897	            m << "MS1 mass axis centred on " << extracted_ppm_offset_
  1898	              << " ppm (borrowed from the fragment fit); median residual against "
  1899	                 "that centre " << ms1_resid << " ppm";
  1900	            if (std::isfinite(ms1_resid) && std::abs(ms1_resid) > 3.0)
  1901	            {
  1902	              m << " -- FAR FROM ZERO, so the fragment offset is not the MS1 offset "
  1903	                   "and var_ms1_coelution is measured through a mis-centred window";
  1904	              writeLogWarn_(m.str());
  1905	            }
  1906	            else { writeLogInfo_(m.str()); }
  1907	          }
  1908	        }
  1909	      }
  1910	      const double secs = std::chrono::duration<double>(
  1911	        std::chrono::steady_clock::now() - t0).count();
  1912	      if (ms1_traces_.empty())
  1913	      {
  1914	        writeLogInfo_("MS1: the run carries no MS1 spectra, so var_ms1_coelution "
  1915	                      "will be absent rather than zero");
  1916	      }
  1917	      else
  1918	      {
  1919	        writeLogInfo_(ms1_traces_.describe() + ", in " + std::to_string(secs) + " s");
  1920	      }
  1921	    }
  1922	    sink.ms1Available(ms1_traces_.empty() ? nullptr : &ms1_traces_);
  1923	
  1924	    // -out_ms1_iso: the cohort M/M+1/M+2 export. Once per run (the MS1 grid and
  1925	    // calibration are per-run, not per-pass), separate cohort-restricted builds
```

```text  (lines 2378-2400)
  2378	      // is not scale-invariant: 4 starves a 2,450-precursor library (measured:
  2379	      // anchors 709 -> 149, identifications 1862 -> 635) while at 4.26 M it
  2380	      // would still leave a million precursors, far more than any fit needs.
  2381	      pass_offset_ = static_cast<std::size_t>(std::max(0, getIntOption_("pass1_offset")));
  2382	      const int target = std::max(0, getIntOption_("pass1_precursors"));
  2383	      const std::size_t n_prec = library.precursorCount();
  2384	      pass_stride_ = (target > 0 && n_prec > static_cast<std::size_t>(target))
  2385	                       ? n_prec / static_cast<std::size_t>(target)
  2386	                       : 1;
  2387	      if (pass_stride_ > 1)
  2388	      {
  2389	        std::ostringstream st;
  2390	        st << "pass 1 extracts every " << pass_stride_ << "th precursor ("
  2391	           << (library.precursorCount() / pass_stride_) << " of "
  2392	           << library.precursorCount() << "); it exists to harvest anchors, and "
  2393	           << "extracting all of them costs memory proportional to the library";
  2394	        writeLogInfo_(st.str());
  2395	      }
  2396	      // external_irt_, NOT false. When -irt_slope/-irt_intercept were supplied
  2397	      // the library's irt was already converted to run seconds above, and
  2398	      // extractInto_ reads this flag to decide whether to apply that map --
  2399	      // so passing false applied it a SECOND time. On IH1 with the stored v6
  2400	      // map that is 7.77 * (run seconds) + 733, i.e. ~24,000 s for a 5,400 s
```

```text  (lines 4425-4440)
  4425	      std::max(1, getIntOption_("max_candidates")));
  4426	    options.match_decoy_candidate_counts = !getFlag_("no_match_decoy_n");
  4427	    // v1.14: rank pooling is a FINAL-scoring rule; pass 1 and the RT refinement keep the native
  4428	    // pooling so the calibration anchors (selected at pass-1 q) are identical to the flag-off run.
  4429	    options.fold_pool_rank = calibrating_ ? false : getFlag_("fold_pool_rank");
  4430	    options.transition_mask = calibrating_ ? nullptr : transition_mask_.get();
  4431	    // -out_fragvec, FINAL scoring only. Pass 1 and the RT refinement compute
  4432	    // the identical columns and nothing would ever read them: the export is
  4433	    // keyed on -out's row order, and -out is written from the final result.
  4434	    options.fragvec = !calibrating_ && !out_fragvec_.empty();
  4435	    {
  4436	      const double mc = getDoubleOption_("mass_accuracy_centre");
  4437	      options.mass_accuracy_centre = (calibrating_ || mc > 1e8) ? std::numeric_limits<double>::quiet_NaN() : mc;
  4438	    }
  4439	    options.train_fdr_initial = getDoubleOption_("train_fdr_initial");
  4440	    options.train_fdr = getDoubleOption_("train_fdr");
```

```text  (lines 4475-4490)
  4475	    options.gate_calibration_n =
  4476	      static_cast<std::size_t>(std::max(100, getIntOption_("gate_calibration_n")));
  4477	    options.gate_calibration_rt_min = getDoubleOption_("gate_calibration_rt_min");   // v1.17, both passes
  4478	    options.empty_trace_min_transitions =
  4479	      static_cast<std::size_t>(std::max(1, getIntOption_("empty_trace_min_transitions")));
  4480	    options.threads = static_cast<unsigned>(std::max(1, getIntOption_("threads")));
  4481	    // Built once per run and owned by the tool; null until then, and null
  4482	    // forever on a run with no MS1, in which case MS1_COELUTION is NaN for every
  4483	    // row and the constant-column guard drops it.
  4484	    options.ms1 = ms1_traces_.empty() ? nullptr : &ms1_traces_;
  4485	    // Sized and cleared by the caller, once per scoring pass. Left null when
  4486	    // the table was not asked for, which is the default.
  4487	    options.terminal_reason = terminal_reasons_.empty() ? nullptr
  4488	                                                        : terminal_reasons_.data();
  4489	    options.oracle_rt = oracle_rt_.empty() ? nullptr : oracle_rt_.data();
  4490	    options.oracle_rt_tol = getDoubleOption_("oracle_rt_tol");
```

```text  (lines 4705-4750)
  4705	  /// scorer consumes, from the SAME extraction -- that is `-out_chrom`. It is
  4706	  /// deliberately the only way a scoring run can obtain chromatograms, so the
  4707	  /// dump cannot fork the scoring: same sink, same sub-scores, same
  4708	  /// calibration surface, whether or not anything collects.
  4709	  ExitCodes extractAndScore_(const ODIA::Library& library, const std::string& run,
  4710	                             double rt_window_override, bool library_rt_is_run_seconds,
  4711	                             ODIA::PeakGroupScorer::Result& scored,
  4712	                             ODIA::ChromatogramCollector* collect = nullptr)
  4713	  {
  4714	    resetTerminalReasons_(library);
  4715	    loadOracleRt_(library);
  4716	    auto options = scoringOptions_();
  4717	    options.library_rt_is_run_seconds = scoring_rt_is_run_seconds_;
  4718	    ODIA::PeakGroupScorer::Sink sink(library, options);
  4719	    std::optional<TeeChromatogramSink> tee;
  4720	    if (collect != nullptr) { tee.emplace(sink, *collect); }
  4721	    ODIA::ChromatogramSink& into =
  4722	      tee ? static_cast<ODIA::ChromatogramSink&>(*tee)
  4723	          : static_cast<ODIA::ChromatogramSink&>(sink);
  4724	    const auto t = std::chrono::steady_clock::now();
  4725	    const auto rc = extractInto_(library, run, into, rt_window_override,
  4726	                                 library_rt_is_run_seconds);
  4727	    if (rc != EXECUTION_OK) { return rc; }
  4728	    // Only now does the fragment mass calibration's verdict exist -- the probe
  4729	    // runs during extraction setup, after the Sink was built. Deciding earlier
  4730	    // ran pass 1 with features pass 2 rejects, and pass 1 supplies the anchors.
  4731	    sink.disableSubScores(ablatedSubScores_());
  4732	    try
  4733	    {
  4734	      scored = sink.finish();
  4735	    }
  4736	    catch (const std::exception& e)
  4737	    {
  4738	      writeLogError_(std::string("Scoring failed: ") + e.what());
  4739	      return INTERNAL_ERROR;
  4740	    }
  4741	    const auto ms = std::chrono::duration<double, std::milli>(
  4742	                      std::chrono::steady_clock::now() - t).count();
  4743	    reportScoring_(scored, options.classifier, ms, "scored on the fly");
  4744	    return EXECUTION_OK;
  4745	  }
  4746	
  4747	  // runScoring_ -- score from a held `Chromatograms` -- is deliberately GONE.
  4748	  // It was the -out_chrom scoring fork: it lost the residual planes and the
  4749	  // MS1 hookup the streaming Sink gets (Mass.Ppm NaN on every row, 7
  4750	  // sub-scores dropped against 4), and its callers skipped the refinement
```

### src/io/MzPeakSource.cpp (reader)

```text  (lines 363-410)
   363	    void peaks(std::size_t begin, std::size_t end,
   364	               std::vector<SpectrumPeaks>& out) override
   365	    {
   366	      if (begin > end || end > info_.size())
   367	      {
   368	        throw std::out_of_range("peaks(): range outside the run");
   369	      }
   370	      out.assign(end - begin, SpectrumPeaks{});
   371	
   372	      // One request per DISTINCT physical spectrum, not one per entry.
   373	      //
   374	      // `info_` holds one entry per (spectrum, isolation window), because a
   375	      // diaPASEF frame co-packs several windows into one physical spectrum and
   376	      // each of them is a separate thing to extract. Those entries share a
   377	      // `index`. mzPeak does no deduplication and no inter-request caching, so
   378	      // asking for the same index twice costs twice: measured on IH1,
   379	      // 8.500 ms/request for 4,096 distinct indices against 8.248 ms/request
   380	      // for 2,048 indices asked for twice each -- flat per request, regardless
   381	      // of whether the frame was just decoded.
   382	      //
   383	      // IH1 has 17,448 physical spectra and 32,210 entries, so the un-deduped
   384	      // request list asked for every MS2 frame twice and spent ~297 s of its
   385	      // ~594 s decode re-decoding what it already had.
   386	      std::vector<std::size_t> want;
   387	      want.reserve(end - begin);
   388	      // Position in [begin, end) -> which entry of `want`, and therefore of the
   389	      // batch, carries its peaks. `out` is indexed by position in the range,
   390	      // NOT by file index, so this indirection is what keeps the two apart.
   391	      std::vector<std::size_t> slot(end - begin, 0);
   392	      {
   393	        std::unordered_map<std::size_t, std::size_t> first_request;
   394	        first_request.reserve((end - begin) * 2);
   395	        for (std::size_t i = begin; i < end; ++i)
   396	        {
   397	          const auto [it, fresh] = first_request.emplace(info_[i].index, want.size());
   398	          if (fresh) { want.push_back(info_[i].index); }
   399	          slot[i - begin] = it->second;
   400	        }
   401	      }
   402	      auto batch = spectra_.get_spectra_batch(want);
   403	
   404	      for (std::size_t k = 0; k < out.size(); ++k)
   405	      {
   406	        const std::size_t s = slot[k];
   407	        if (s >= batch.size()) { continue; }
   408	        auto& dst = out[k];
   409	        const auto& mz = batch[s].mz();
   410	        const auto& intensity = batch[s].intensity();
```

```text  (lines 451-457)
   451	  private:
   452	    std::string filename_;
   453	    std::vector<SpectrumInfo> ms1_info_;
   454	    MzPeak::Index index_;
   455	    MzPeak::Spectra spectra_;
   456	    std::vector<SpectrumInfo> info_;
   457	    std::vector<IsolationWindow> windows_;
```

### src/extract/Ms1Traces.cpp (dense MS1 matrix build; whole build loop)

```text  (lines 60-183)
    60	    // index over the kept precursors, and the caller gets that order back via
    61	    // kept_indices -- the writer reconstructs ids from it, so the mapping is
    62	    // never inferred twice.
    63	    std::vector<std::uint32_t> row;
    64	    std::size_t rows = np;
    65	    if (keep != nullptr)
    66	    {
    67	      row.assign(np, UINT32_MAX);
    68	      std::uint32_t r = 0;
    69	      for (std::size_t i = 0; i < np; ++i)
    70	      {
    71	        if (i < keep->size() && (*keep)[i]) { row[i] = r++; }
    72	      }
    73	      rows = r;
    74	      if (kept_indices != nullptr)
    75	      {
    76	        kept_indices->clear();
    77	        kept_indices->reserve(rows);
    78	        for (std::size_t i = 0; i < np; ++i)
    79	        { if (row[i] != UINT32_MAX) { kept_indices->push_back(static_cast<std::uint32_t>(i)); } }
    80	      }
    81	    }
    82	    out.values_.assign(rows * out.bins_, 0.0f);
    83	
    84	    // Search the sorted LIBRARY side and iterate the peaks: SpectrumSource
    85	    // documents that a peak array is not ascending in m/z (a mobility frame
    86	    // concatenates its TIMS scans), and a binary search over it "does not fail
    87	    // loudly -- it returns near-zero matches", which is indistinguishable from
    88	    // an ion that is not there.
    89	    struct Target { double mz; std::uint32_t slot; float im; };
    90	    std::vector<Target> idx;
    91	    idx.reserve(np);
    92	    for (std::size_t i = 0; i < np; ++i)
    93	    {
    94	      if (keep != nullptr && row[i] == UINT32_MAX) { continue; }
    95	      // CALIBRATED, like the fragment axis. This matched on the library's
    96	      // THEORETICAL m/z with a symmetric window and no offset, while the
    97	      // fragment extractor was centred on the fitted deviation -- on IH1 that
    98	      // is -10.0108 ppm against a +/-10 ppm half-width, so a precursor whose
    99	      // MS1 error resembles its MS2 error sat at the window EDGE and a weak one
   100	      // fell out entirely. That matters because ms1_coelution is the main
   101	      // evidence for calling a precursor ABSENT (median -0.093 for the 10,736
   102	      // DIA-NN precursors we reject, against -0.124 for the 801,458 bulk
   103	      // non-identifications and 0.473 for accepted ones), and absence cannot be
   104	      // concluded from an uncalibrated measurement.
   105	      //
   106	      // The isotope offset is applied BEFORE the calibration scaling, and per
   107	      // this precursor's own charge: the M+k target is `mz + k*dm/z`, and the
   108	      // instrument's relative (ppm) error then applies to that target as it
   109	      // does to any mass.
   110	      const int z = p.charge[i] > 0 ? static_cast<int>(p.charge[i]) : 1;
   111	      const double mz = (fromFixed(p.mz[i]) + isotope_offset_da / z) *
   112	                        (1.0 + ppm_offset * 1e-6);
   113	      const std::uint32_t slot = keep != nullptr ? row[i] : static_cast<std::uint32_t>(i);
   114	      if (mz > 0.0) { idx.push_back({mz, slot, p.im[i]}); }
   115	    }
   116	    std::sort(idx.begin(), idx.end(),
   117	              [](const Target& a, const Target& b) { return a.mz < b.mz; });
   118	
   119	    // CAPPED. The first version pushed one double per (peak, target) match over
   120	    // the whole run and reached 591 GB RSS against v3's 116 GB peak, blowing
   121	    // through -live_memory_gb on a shared node. The median of a bounded prefix
   122	    // is the same number to far more precision than it is worth: this is a
   123	    // diagnostic, not a fit.
   124	    static constexpr std::size_t RESID_CAP = 1u << 21;   // 2M samples, 16 MB
   125	    std::vector<double> resid;
   126	    if (observed_ppm_median != nullptr) { resid.reserve(RESID_CAP); }
   127	    std::vector<SpectrumPeaks> block;
   128	    const std::size_t STEP = 64;
   129	    for (std::size_t b = 0; b < ms1.size(); b += STEP)
   130	    {
   131	      const std::size_t e = std::min(b + STEP, ms1.size());
   132	      source.ms1Peaks(b, e, block);
   133	      for (std::size_t s = 0; s < block.size(); ++s)
   134	      {
   135	        const auto& sp = block[s];
   136	        const bool gated = im_window > 0.0 && sp.ion_mobility.size() == sp.mz.size();
   137	        for (std::size_t k = 0; k < sp.mz.size(); ++k)
   138	        {
   139	          const double m = sp.mz[k], tol = m * fragment_ppm * 1e-6;
   140	          auto it = std::lower_bound(idx.begin(), idx.end(), m - tol,
   141	                                     [](const Target& a, double v) { return a.mz < v; });
   142	          for (; it != idx.end() && it->mz <= m + tol; ++it)
   143	          {
   144	            if (std::abs(it->mz - m) > it->mz * fragment_ppm * 1e-6) { continue; }
   145	            // Same rule as the MS2 match loop, deliberately: skip the gate when
   146	            // EITHER side is unknown, because absent information is not evidence
   147	            // of mismatch. Testing `abs(NaN - x) <= w` is false, so a precursor
   148	            // with no library 1/K0 had every MS1 peak rejected and came out with
   149	            // a NaN MS1_COELUTION -- while its MS2 side was extracted ungated.
   150	            if (gated && !std::isnan(static_cast<double>(it->im)))
   151	            {
   152	              const double d = std::abs(static_cast<double>(sp.ion_mobility[k]) -
   153	                                        static_cast<double>(it->im));
   154	              if (!(d <= im_window)) { continue; }
   155	            }
   156	            if (observed_ppm_median != nullptr && resid.size() < RESID_CAP)
   157	            {
   158	              // Residual against the CALIBRATED target, so a correct offset
   159	              // centres this on 0 and a wrong one does not.
   160	              resid.push_back((m - it->mz) / it->mz * 1e6);
   161	            }
   162	            float& c = out.values_[it->slot * out.bins_ + (b + s)];
   163	            // Max, not sum: a mobility-merged frame holds the same ion in
   164	            // several scans, and summing would make the trace a function of how
   165	            // many scans it spans rather than of how much ion is present.
   166	            c = std::max(c, sp.intensity[k]);
   167	          }
   168	        }
   169	      }
   170	    }
   171	    if (observed_ppm_median != nullptr)
   172	    {
   173	      if (resid.empty()) { *observed_ppm_median = std::numeric_limits<double>::quiet_NaN(); }
   174	      else
   175	      {
   176	        const std::size_t h = resid.size() / 2;
   177	        std::nth_element(resid.begin(), resid.begin() + h, resid.end());
   178	        *observed_ppm_median = resid[h];
   179	      }
   180	    }
   181	    return out;
   182	  }
   183	} // namespace ODIA
```

### include/odia/Ms1Traces.h

```text  (lines 90-109)
    90	    /// @param keep  optional per-precursor mask (size precursorCount). When
    91	    ///        given, only masked-in precursors get a dense row -- the matrix is
    92	    ///        `kept x bins` instead of `np x bins`, which is what makes a 3x
    93	    ///        isotope build affordable (the full-library dense matrix is 77.6
    94	    ///        GB at 4.99M precursors and is charged before extraction begins).
    95	    ///        With @p keep, `at()` takes ROW indices, not library indices;
    96	    ///        @p kept_indices receives the library index of each row.
    97	    static Ms1Traces build(const Library& library, SpectrumSource& source,
    98	                           double fragment_ppm, double im_window,
    99	                           double ppm_offset = 0.0,
   100	                           double* observed_ppm_median = nullptr,
   101	                           double isotope_offset_da = 0.0,
   102	                           const std::vector<std::uint8_t>* keep = nullptr,
   103	                           std::vector<std::uint32_t>* kept_indices = nullptr);
   104	
   105	  private:
   106	    std::vector<float> times_;
   107	    std::vector<float> values_;      ///< precursor-major, `bins_` per precursor
   108	    std::size_t bins_ = 0;
   109	  };
```

### src/score/PeakGroupScorer.cpp (the sink: Gate C null, Session::add, plane reads, MS1 read, group push, sort)

```text  (lines 1240-1282)
  1240	      std::size_t zeros = 0;                     ///< sample statistics exactly 0
  1241	      double rt_first = std::numeric_limits<double>::quiet_NaN();
  1242	      double rt_last = std::numeric_limits<double>::quiet_NaN();
  1243	
  1244	      /// Returns true if the precursor should be admitted. @p rt_centre is the
  1245	      /// window's position in run seconds; with @p rt_min > 0 a decoy before it
  1246	      /// is admitted (as every precursor is while the null is built) but does
  1247	      /// not enter the calibration sample (v1.17). rt_min 0 = off.
  1248	      bool admit(double stat, bool is_decoy, std::size_t n_needed, double alpha,
  1249	                 double rt_centre, double rt_min)
  1250	      {
  1251	        std::lock_guard<std::mutex> g(mu);
  1252	        if (!ready)
  1253	        {
  1254	          if (is_decoy && !(rt_min > 0.0 && rt_centre < rt_min))
  1255	          {
  1256	            decoy_stats.push_back(stat);
  1257	            if (stat == 0.0) { ++zeros; }
  1258	            if (decoy_stats.size() == 1) { rt_first = rt_centre; }
  1259	            rt_last = rt_centre;
  1260	          }
  1261	          if (decoy_stats.size() >= n_needed)
  1262	          {
  1263	            std::sort(decoy_stats.begin(), decoy_stats.end());
  1264	            const std::size_t k = std::min(decoy_stats.size() - 1,
  1265	              std::size_t((1.0 - alpha) * double(decoy_stats.size())));
  1266	            tau = decoy_stats[k];
  1267	            ready = true;
  1268	            // v1.17: say so, once per Session, on the channel the picker census
  1269	            // uses (stderr -> the run log). Nothing was logged before, which is
  1270	            // how tau = 0 from pre-gradient windows went unnoticed.
  1271	            std::fprintf(stderr,
  1272	                         "gate C null armed: n=%zu decoys, tau=%.6g, zeros=%.4f (%zu), "
  1273	                         "first/last calibration RT %.1f-%.1f s, rt_min %.1f s\n",
  1274	                         decoy_stats.size(), tau,
  1275	                         double(zeros) / double(decoy_stats.size()), zeros,
  1276	                         rt_first, rt_last, rt_min);
  1277	          }
  1278	          ++admitted_uncalibrated;
  1279	          return true;                 // admit while the null is being built
  1280	        }
  1281	        return stat >= tau;
  1282	      }
```

```text  (lines 1310-1330)
  1310	    std::mutex rejects_registry_mutex_;
  1311	    std::vector<PickerRejects*> rejects_registry_;
  1312	
  1313	    struct RegisteredRejects : PickerRejects
  1314	    {
  1315	      RegisteredRejects()
  1316	      {
  1317	        std::lock_guard<std::mutex> g(rejects_registry_mutex_);
  1318	        rejects_registry_.push_back(this);
  1319	      }
  1320	    };
  1321	    thread_local RegisteredRejects rejects_;
  1322	
  1323	    /// Every thread's counters, summed. The only correct way to read them.
  1324	    PickerRejects totalRejects()
  1325	    {
  1326	      PickerRejects t;
  1327	      std::lock_guard<std::mutex> g(rejects_registry_mutex_);
  1328	      for (const PickerRejects* r : rejects_registry_)
  1329	      {
  1330	        for (int c = 0; c < 2; ++c)
```

```text  (lines 1355-1406)
  1355	  struct PeakGroupScorer::Session::GateNull : NullCalibrationBody {};
  1356	
  1357	  PeakGroupScorer::Session::Session(const Library& library, const Options& options)
  1358	    : gate_null_(std::make_shared<GateNull>()), library_(&library), options_(options)
  1359	  {
  1360	  }
  1361	
  1362	  void PeakGroupScorer::Session::add(const PrecursorChromatogram& chromatogram_in)
  1363	  {
  1364	    Result& result = result_;
  1365	    const Options& options = options_;
  1366	    const Library& library = *library_;
  1367	    const auto& p = library.precursors();
  1368	    const auto& t = library.transitions();
  1369	    const PrecursorChromatogram& chromatogram = chromatogram_in;   // rebound per candidate under a mask (v1.15)
  1370	
  1371	    const std::size_t i = chromatogram.precursor;
  1372	    const std::uint32_t tb = chromatogram.transition_begin;
  1373	    const std::uint32_t tc = chromatogram.transition_count;
  1374	    // Every return below records why. `mark` is a no-op unless
  1375	    // -out_terminal_reasons asked for the table.
  1376	    // The EXTRACTOR's reasons win. It knows things this function cannot -- that
  1377	    // no isolation window covers the precursor, or that the prefilter dropped
  1378	    // it -- and if the scorer overwrites them the specific reason is replaced
  1379	    // by a vaguer one that is also true.
  1380	    //
  1381	    // That is not hypothetical. On full_v9, `few_points` came out at 1,119,490
  1382	    // and `no_window_coverage` at ZERO, and 100.0% of the `few_points` targets
  1383	    // turned out to lie outside every window's m/z range (median m/z 1468.3
  1384	    // against 525.8 for scored precursors; the windows stop at 1400.62 Th).
  1385	    // The whole bucket was the uncovered population wearing the wrong label,
  1386	    // which is exactly the kind of misattribution this table exists to prevent.
  1387	    const auto mark = [&](TerminalReason r) {
  1388	      if (!options.terminal_reason) { return; }
  1389	      const std::uint8_t prior = options.terminal_reason[i];
  1390	      if (prior == static_cast<std::uint8_t>(TerminalReason::NoWindowCoverage) ||
  1391	          prior == static_cast<std::uint8_t>(TerminalReason::PrefilterExcluded))
  1392	      { return; }
  1393	      options.terminal_reason[i] = static_cast<std::uint8_t>(r);
  1394	    };
  1395	    if (tc == 0) { mark(TerminalReason::NoTransitions); return; }
  1396	
  1397	    const std::size_t points = chromatogram.pointCount(0);
  1398	    // The class is needed BEFORE these gates, not after: `reached` is
  1399	    // incremented past them, so an asymmetry here is invisible in every counter
  1400	    // downstream. One thread's numbers already put the entire target/decoy
  1401	    // imbalance in `reached` (137/254) while scans PER PRECURSOR were identical
  1402	    // (2014 both), so the divergence happens at exactly these two returns.
  1403	    const bool is_decoy = p.decoy[chromatogram.precursor] != 0;
  1404	    if (points < 3)
  1405	    { mark(TerminalReason::FewPoints); ++rejects_.no_points[is_decoy]; ++result.precursors_without_candidate; return; }
  1406	
```

```text  (lines 2115-2150)
  2115	      double im_obs = std::numeric_limits<double>::quiet_NaN();
  2116	      if (chromatogram.im_num != nullptr && chromatogram.im_den != nullptr && hi > lo)
  2117	      {
  2118	        double num = 0.0, den = 0.0;
  2119	        for (std::uint32_t k = 0; k < tc; ++k)
  2120	        {
  2121	          const std::uint32_t n = chromatogram.pointCount(k);
  2122	          const std::ptrdiff_t off = chromatogram.trace(k) - chromatogram.points;
  2123	          const float* inum = chromatogram.im_num + off;
  2124	          const float* iden = chromatogram.im_den + off;
  2125	          for (std::size_t j = lo; j <= hi && j < n; ++j)
  2126	          {
  2127	            if (iden[j] > 0.0f) { num += inum[j]; den += iden[j]; }
  2128	          }
  2129	        }
  2130	        if (den > 0.0) { im_obs = num / den; }
  2131	
  2132	        // Per-FRAGMENT observed mobility, and its scatter across the group.
  2133	        // Same construction as mass_ppm_spread: one number per fragment (the
  2134	        // intensity-weighted mean over the cycles it was seen in), then a MAD
  2135	        // across fragments. Averaging within a fragment is legitimate -- those
  2136	        // cells measure one ion -- whereas averaging across fragments is what
  2137	        // would hide the very disagreement being looked for.
  2138	        std::vector<double> per_fragment_im;
  2139	        per_fragment_im.reserve(tc);
  2140	        for (std::uint32_t k = 0; k < tc; ++k)
  2141	        {
  2142	          const std::uint32_t n = chromatogram.pointCount(k);
  2143	          const std::ptrdiff_t off = chromatogram.trace(k) - chromatogram.points;
  2144	          const float* inum = chromatogram.im_num + off;
  2145	          const float* iden = chromatogram.im_den + off;
  2146	          double fn = 0.0, fd = 0.0;
  2147	          for (std::size_t j = lo; j <= hi && j < n; ++j)
  2148	          { if (iden[j] > 0.0f) { fn += inum[j]; fd += iden[j]; } }
  2149	          if (fd > 0.0) { per_fragment_im.push_back(fn / fd); }
  2150	        }
```

```text  (lines 2421-2438)
  2421	      {
  2422	        double r = std::numeric_limits<double>::quiet_NaN();
  2423	        if (options.ms1 != nullptr && !options.ms1->empty() &&
  2424	            i < options.ms1->precursors() && hi > lo)
  2425	        {
  2426	          // `total` is already the per-cycle summed fragment intensity and
  2427	          // [lo, hi] the candidate's own cycle bounds -- reuse both rather than
  2428	          // recomputing a second, subtly different fragment sum.
  2429	          std::vector<double> f, m;
  2430	          f.reserve(hi - lo + 1);
  2431	          m.reserve(hi - lo + 1);
  2432	          for (std::size_t j = lo; j <= hi && j < total.size(); ++j)
  2433	          {
  2434	            const std::size_t b = options.ms1->binFor(chromatogram.retentionTime(
  2435	              static_cast<std::uint32_t>(j)));
  2436	            f.push_back(total[j]);
  2437	            m.push_back(static_cast<double>(options.ms1->at(i, b)));
  2438	          }
```

```text  (lines 2651-2656)
  2651	      // In lockstep with `groups`, and appended at the SAME statement so no
  2652	      // path can push one without the other.
  2653	      if (options.fragvec)
  2654	      { result.fragvec.insert(result.fragvec.end(), fv, fv + N_FRAGVEC); }
  2655	      result.groups.push_back(std::move(g));
  2656	    }
```

```text  (lines 2721-2727)
  2721	  /// Whether refitting after a new retention-time map can change any score.
  2722	  ///
  2723	  /// False since RT_DELTA was removed. Kept as a function rather than deleted
  2724	  /// at the call sites so that adding a map-dependent sub-score re-enables the
  2725	  /// loop by flipping one return, instead of by remembering that a loop was
  2726	  /// deleted somewhere.
  2727	  bool PeakGroupScorer::refitsChangeScores() { return false; }
```

```text  (lines 3143-3157)
  3143	    // Sorted through an index PERMUTATION rather than in place, because
  3144	    // `MassAnchor::group` indexes into this vector. Sorting the groups directly
  3145	    // left every harvested anchor pointing at whatever group moved into its old
  3146	    // slot -- so an anchor from an accepted target could be read with a decoy's
  3147	    // flag and a rejected group's q-value. Silent, and it invalidated a whole
  3148	    // measurement before it was caught.
  3149	    std::vector<std::uint32_t> order(result.groups.size());
  3150	    std::iota(order.begin(), order.end(), 0u);
  3151	    std::stable_sort(order.begin(), order.end(),
  3152	                     [&](std::uint32_t ia, std::uint32_t ib) {
  3153	                       const PeakGroup& a = result.groups[ia];
  3154	                       const PeakGroup& b = result.groups[ib];
  3155	                       if (a.precursor != b.precursor) { return a.precursor < b.precursor; }
  3156	                       if (a.apex_rt != b.apex_rt) { return a.apex_rt < b.apex_rt; }
  3157	                       return a.apex_intensity < b.apex_intensity; });
```

### include/odia/PeakGroupScorer.h

```text  (lines 920-930)
   920	      /// Where to record each precursor's `TerminalReason`, indexed by library
   921	      /// precursor index. Null disables the accounting entirely.
   922	      ///
   923	      /// Written without a lock, and the reason is NOT that sessions are
   924	      /// per-thread -- a Sink owns exactly one Session and hands it to
   925	      /// extraction. `ChromatogramExtractor` calls `accept` from its serial
   926	      /// emit loop in both paths (`ChromatogramExtractor.cpp:900`, `:1279`;
   927	      /// the threaded phase is matching, which finishes first), so `add` runs
   928	      /// on one thread. That same invariant is what already makes
   929	      /// `Session::result_.groups` safe to push to.
   930	      std::uint8_t* terminal_reason = nullptr;
```

```text  (lines 1250-1256)
  1250	      Sink(const Library& library, const Options& options) : session_(library, options) {}
  1251	
  1252	      /// See ChromatogramSink::ms1Available -- the traces do not exist yet when
  1253	      /// this sink is constructed.
  1254	      void ms1Available(const Ms1Traces* m) override { session_.setMs1Traces(m); }
  1255	      void accept(const PrecursorChromatogram& trace) override { session_.add(trace); }
  1256	      Result finish() { return session_.finish(); }
```

### include/odia/scoring/lda.h (semi-supervised loop; nested OMP sizing)

```text  (lines 587-605)
   587	  // Threads for the NESTED regions inside the fold loop. Computed once and applied to every one of
   588	  // them: the training loop already honoured it, but the per-group scans did not, so each of the 3
   589	  // concurrent folds opened teams of the FULL thread count -- 540 threads on 224 cores at
   590	  // OMP_NUM_THREADS=180. Oversubscription of that size costs more in scheduling than the
   591	  // parallelism returns.
   592	  int inner_threads = 1;
   593	#ifdef _OPENMP
   594	  // BOTH inner-parallel learners need the active-level raised, not just the GBT. With it set for
   595	  // the GBT alone, the network's chunk loop ran as a team of ONE inside the fold loop: measured
   596	  // 435% CPU on a 224-core node, i.e. the 3 folds and nothing else. Splitting as
   597	  // folds x (threads/folds) uses the machine without oversubscribing it, and neither learner's
   598	  // result depends on the number (odia_gbt_test T8; odia_nn_test thread invariance).
   599	  if (params.classifier == LDAParams::Classifier::GBT || params.classifier == LDAParams::Classifier::NN)
   600	  {
   601	    omp_set_max_active_levels(2);
   602	    inner_threads = std::max(1, omp_get_max_threads() / std::max(1, folds));
   603	    gbt_params.n_threads = inner_threads;
   604	    nn_params.n_threads = inner_threads;
   605	  }
```

```text  (lines 1142-1167)
  1142	    for (int iteration = 0; iteration < std::max(0, params.n_iter); ++iteration)
  1143	    {
  1144	      // Reduce training scores to the best candidate row per precursor, then
  1145	      // estimate group-level q-values for confident positive selection.
  1146	      std::vector<lda_detail::RankedGroup> ranked;
  1147	      ranked.reserve(train_groups.size());
  1148	      best_rows_of(train_groups, bg_row, bg_score);
  1149	      for (std::size_t i = 0; i < train_groups.size(); ++i)
  1150	      {
  1151	        const std::size_t g = train_groups[i];
  1152	        ranked.push_back({g, bg_row[i], group_label[g], bg_score[i], 1.0});
  1153	      }
  1154	      lda_detail::assignQValues(ranked, params.use_pi0);
  1155	
  1156	      // First pass is seeded by a SINGLE feature, so a strict cut can select too few positives to
  1157	      // fit anything and the model never bootstraps -- the failure this split exists to prevent.
  1158	      const double raw_fdr = (iteration == 0) ? params.train_fdr_initial : params.train_fdr;
  1159	      const double train_fdr = std::max(0.0, std::min(1.0, raw_fdr));
  1160	      std::vector<std::size_t> positive_rows;
  1161	      for (const auto& candidate : ranked)
  1162	      {
  1163	        if (candidate.label == 1 && candidate.qvalue <= train_fdr)
  1164	        {
  1165	          positive_rows.push_back(candidate.best_row);
  1166	        }
  1167	      }
```

### include/odia/scoring/gbt.h (the only OpenMP regions)

```text  (lines 283-290)
   283	    const std::size_t n_rows = rows.size();
   284	    std::vector<uint8_t> B(n_rows * n_features_, 0);
   285	#ifdef _OPENMP
   286	#pragma omp parallel for schedule(static) num_threads(p.n_threads > 0 ? p.n_threads : omp_get_max_threads())
   287	#endif
   288	    for (long long ii = 0; ii < static_cast<long long>(n_rows); ++ii)
   289	    {
   290	      const std::size_t i = static_cast<std::size_t>(ii);
```

```text  (lines 488-510)
   488	      // make the number of partial sums -- and therefore the floating-point reduction order --
   489	      // depend on the thread count, and these scores become q-values. Chunk count depends only on
   490	      // the row count, threads take chunks in any order, and the combination below runs in chunk
   491	      // index order, so the result is bit-identical at 1, 8 or 224 threads.
   492	      const int nchunk = static_cast<int>(std::min<std::size_t>(32, std::max<std::size_t>(1, n_rows / 1000)));
   493	      // Scratch is owned by fit() and reused across every level of every tree: allocating it here
   494	      // cost 960 vector constructions per fit (2 per level x 4 levels x 120 trees) for buffers of
   495	      // identical shape. Only the portion this level uses is cleared.
   496	      const std::size_t need = static_cast<std::size_t>(nchunk) * hist_size;
   497	      if (scratch.size() < need) { scratch.resize(need); }
   498	      std::fill(scratch.begin(), scratch.begin() + static_cast<std::ptrdiff_t>(need), gbt_detail::GH{});
   499	      std::vector<std::size_t> node_rows_c(static_cast<std::size_t>(nchunk) * n_at_level, 0);
   500	
   501	#ifdef _OPENMP
   502	// ponytail: the loop has exactly nchunk iterations, so threads beyond nchunk get no work and only
   503	// add barrier and scratch-page cost. Measured: 224 threads ran this fit 60% SLOWER than 32.
   504	// Raise the cap by raising nchunk -- which changes the reduction order, so it is a deliberate act.
   505	#pragma omp parallel for schedule(dynamic, 1) num_threads(std::min(n_threads > 0 ? n_threads : 1, nchunk))
   506	#endif
   507	      for (int c = 0; c < nchunk; ++c)
   508	      {
   509	        const std::size_t lo = n_rows * static_cast<std::size_t>(c) / static_cast<std::size_t>(nchunk);
   510	        const std::size_t hi = n_rows * static_cast<std::size_t>(c + 1) / static_cast<std::size_t>(nchunk);
```

## 4. DOC EXCERPTS (project documents, `cat -n` line numbers; both predate the 9.9M-precursor library and were measured at the ≤10 GiB scale on a ≤100k-precursor library)

### doc/10-performance-report.md §2 (Phase 2 — extraction), excerpts

```text  (lines 167-227)
   167	## 2. Phase 2 — extraction
   168	
   169	`OpenDIAlyzer -tr <lib> -in IH1_diaPASEF.mzpeak -stop_after extract -out_chrom
   170	<tsv> -threads 8`, run file on local NVMe, `-mass_calibration off` except where
   171	noted. No iRT calibration was supplied, so RT windows span essentially the whole
   172	gradient (1,342 of 1,342 cycles) — this is the worst case for point count and it
   173	is stated explicitly because every number below scales with it.
   174	
   175	### Stage table and scaling with library size
   176	
   177	| precursors | transitions | points (allocated) | nonzero | chrom MiB | decode s | index s | match s | extract s | wall | CPU% | peak RSS |
   178	|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
   179	| 500 † | 6,000 | 8,359,680 | 3,798,582 | 32.1 | 602.24 | 0.006 | 12.25 | 614.5 | 11:58 | 108% | 14.19 GiB |
   180	| 2,500 | 30,000 | 41,508,708 | 18,486,918 | 158.9 | 599.76 | 0.134 | 14.21 | 614.1 | 11:37 | 110% | 14.32 GiB |
   181	| 9,522 | 114,256 | 153,363,800 | 68,059,302 | 586.9 | 593.83 | 0.507 | 18.96 | 613.3 | 15:14 | 112% | **14.74 GiB** |
   182	| 9,522 ‡ | 114,256 | 153,363,800 | 96,108,130 | 586.9 | 590.95 | 0.506 | 21.18 | 612.6 | 15:21 | 113% | 14.76 GiB |
   183	| 100,000 | 1,199,962 | 1,610,662,790 | 695,012,675 | 6,162.6 | 596.30 | 5.219 | 61.72 | 663.2 | — § | — | **20.26 GiB** |
   184	
   185	† mass calibration **on** (adds 84.45 s, see below) and read from Ceph, not NVMe.
   186	‡ `-precursor_im_window 0`, i.e. today's mobility filter disabled.
   187	§ stopped deliberately during the chromatogram write, which would have produced a
   188	~76 GB TSV at ~24 MB/s; extraction itself completed and is reported in full.
   189	
   190	**A 200x range of library size (500 → 100,000 precursors) moves decode by 1.1%**
   191	(602.2 → 596.3 s). Points per transition is 1,342.3 at every size — exactly the
   192	cycle count — confirming the CSR is allocated, not data-driven.
   193	
   194	The 100,000-precursor row is also the sharpest test of the RSS attribution
   195	below: 20.26 GiB total − 6.02 GiB of chromatogram CSR = **14.24 GiB**, against
   196	**14.15 GiB** of non-CSR memory measured independently at 9,522 precursors.
   197	**Agreement to 0.7%**, across a 10x change in library size.
   198	
   199	Two stages are missing from ODIA's own instrumentation and had to be obtained by
   200	subtraction from total wall:
   201	
   202	* **Mass calibration** (`-mass_calibration auto`, the **default**): **84.45 s**,
   203	  decoding 3,840 spectra. It is reported only inside the calibration report's
   204	  `collect_seconds` and is *not* part of `decode_seconds`. It is single-threaded
   205	  and `-threads` never reaches it.
   206	* **Chromatogram TSV write**: at 9,522 precursors, wall 914 s − extract 613 s =
   207	  **~301 s** to write 7.16 GB. Untimed and unlogged.
   208	
   209	### What scales and what does not
   210	
   211	* **decode is flat at 596–602 s across a 200x range of library size.** It is a
   212	  function of the *file*, not of the library.
   213	* **index (0.006 → 5.2 s) and match (12.2 → 61.7 s) grow with the library but
   214	  stay ≤10% of extraction even at 100,000 precursors** — 3.5% at the realistic
   215	  9,522. The only threaded stage in Phase 2 is not the problem; that is why the
   216	  whole run sits at 108–113% CPU no matter what `-threads` says.
   217	* **`points` scales linearly with transitions** at 1,342 points per transition
   218	  (= every cycle in the run). It is *allocated*, not data-dependent.
   219	* **Peak RSS = a constant ~14.2 GiB of reader + 4 bytes per allocated point.**
   220	  That is the whole memory model of Phase 2, and both terms are measured.
   221	
   222	### RSS attribution — measured, not guessed
   223	
   224	Peak RSS is **14.2–14.8 GiB at every library size from 500 to 9,522
   225	precursors**, while ODIA's chromatogram arrays over that range grow from 32 MiB
   226	to 587 MiB. **At least 96% of resident memory is not ODIA's output.** Only at
   227	100,000 precursors does ODIA's own CSR (6.02 GiB) become a comparable term.
```

```text  (lines 254-272)
   254	3. **It is not allocator slack.** `MALLOC_MMAP_THRESHOLD_=131072`
   255	   `MALLOC_TRIM_THRESHOLD_=131072` reduced peak RSS only 6.04 → 5.36 GiB (−11%)
   256	   and cost 16% more time; `MALLOC_ARENA_MAX=1` changed nothing. **~89% is
   257	   genuine retention.**
   258	
   259	4. **ODIA's own share is small and matches theory.** Decode-only peak RSS at
   260	   4,096 spectra is 0.93 GiB; adding ODIA's exact copy loop takes it to
   261	   1.45 GiB. So ODIA's 1024-spectrum `SpectrumPeaks` block costs **~0.52 GiB**,
   262	   against a predicted 1,024 × 25,659 peaks × 16 B = 0.42 GiB.
   263	
   264	**Attribution at 9,522 precursors (14.74 GiB peak):**
   265	
   266	| holder | GiB | share | whose |
   267	|---|---:|---:|---|
   268	| mzPeak `Spectra` retained decode buffers | ~13.5 | **92%** | mzPeak reader |
   269	| ODIA's 1024-spectrum peak block | ~0.52 | 3.5% | ODIA |
   270	| ODIA chromatogram CSR (`intensity` + offsets) | 0.587 | 4.0% | ODIA |
   271	| ODIA library + m/z index | ~0.03 | 0.2% | ODIA |
   272	
```

```text  (lines 359-376)
   359	### Phase 2 bottlenecks, ranked
   360	
   361	Against 914 s wall at 9,522 precursors:
   362	
   363	1. **mzPeak decode — 593.8 s (65%).** Splits into two very different halves:
   364	   * **~297 s: duplicate frame requests.** *ODIA's code — actionable.*
   365	   * **~297 s: the reader itself, ~19x slower than plain Arrow.** *mzPeak,
   366	     upstream.*
   367	2. **Chromatogram TSV write — ~301 s (33%).** *ODIA's code — actionable.*
   368	   7.16 GB at ~24 MB/s; see §3.5 for proof this is formatting, not disk.
   369	3. **Mass calibration — 84.45 s** whenever `-mass_calibration auto` (the
   370	   default) is used. *ODIA's code — actionable.* Single-threaded, decodes 3,840
   371	   spectra, and ignores `-threads` (`MassCalibration::Options::threads` exists
   372	   and is never read).
   373	
   374	`index` (0.51 s) and `match` (19.0 s) together are 2.1% of wall. **The one part
   375	of Phase 2 that is parallelised is the one part that does not need to be.**
   376	
```

### doc/11-memory-and-compaction-plan.md, excerpts

```text  (lines 1-45)
     1	# Where ODIA's memory actually goes, and the order to fix it
     2	
     3	Measured 2026-08-05 overnight, IH1_diaPASEF.mzpeak (12.75 GiB, 32,210 spectra)
     4	unless stated. Every number here came from a run, not a model.
     5	
     6	## The decomposition
     7	
     8	Peak RSS is two terms and nothing else:
     9	
    10	| term | size | scales with |
    11	|---|---:|---|
    12	| decode floor | **~10.3 GiB** | nothing measured so far |
    13	| chromatograms | **3.5 GiB / 95k precursors** | library size x RT window |
    14	
    15	The 4.26M-precursor library (`odia_v5.tsv`) was observed **live at 154.2 GiB
    16	anonymous**. That is the answer to "why can't we use our own library": not a
    17	projection, an observation.
    18	
    19	Chromatogram arithmetic, from the extractor's own log at n100k:
    20	`983,561,430 points (3.66 GiB), 504,467,125 nonzero` = 9,836 points/precursor
    21	at `-rt_window 600`. At 4.26M precursors that is ~156 GiB.
    22	
    23	## Five hypotheses closed by measurement
    24	
    25	Do not reopen these from priors. Each cost a run; each was plausible.
    26	
    27	| hypothesis | test | verdict |
    28	|---|---|---|
    29	| per-thread parquet row-group buffers | 1 / 4 / 16 threads | **refuted** — 10.34 / 10.27 / 10.35 GiB |
    30	| page cache / mmap counted in RSS | `RssAnon` vs `RssFile` | **refuted** — RssFile 31 MiB |
    31	| glibc arena fragmentation | tcmalloc `LD_PRELOAD` | **refuted** — 10.49 GiB, no better |
    32	| ditto, arena proliferation | `MALLOC_ARENA_MAX=2` | **refuted** — 10.34 GiB |
    33	| chromatograms are sparse | nonzero fraction | **refuted** — 51.3% dense, RLE loses |
    34	
    35	The first is the explanation from the mzPeak handoff (224 workers -> 105.8 GB,
    36	capped to 26.5 GB). Its precondition does not hold here: `MzPeakSource` holds one
    37	`Index` and one `Spectra`, and the floor is flat in threads. A
    38	`-mzpeak_decode_threads` equivalent would buy nothing.
    39	
    40	The third and fourth matter because the *prior* OpenDIAlyzer measured 88% of RSS
    41	as arena fragmentation and got 189 -> 106 GB from tcmalloc alone. That result does
    42	not transfer. See [[odia-memory-is-the-problem]].
    43	
    44	## What the floor still could be
    45	
```

```text  (lines 84-130)
    84	## The floor, located (tcmalloc heap profile, 2026-08-06)
    85	
    86	Six hypotheses were refuted by measurement before profiling: per-thread buffers,
    87	page cache, glibc fragmentation, arena proliferation, sparse chromatograms, and
    88	mzPeak's row-group cache (bounded — `cache_.erase(cache_.begin())` at
    89	`kCachedGroups`, ~21 MB/group). Reading code produced six wrong answers; the
    90	profile produced the right one in one run.
    91	
    92	`HEAPPROFILE` + `libtcmalloc.so.4`, tiny library, IH1. Peak **live** heap 9.45 GiB
    93	(dump 162 of 197; the final dump is 4.7 MB — profile the peak, not the end).
    94	Cumulative: **248 GiB allocated over 23.5M allocations**.
    95	
    96	Top five sites hold 8.76 GiB of the 9.45 GiB peak. Resolved with `addr2line`:
    97	
    98	| bytes | site |
    99	|---:|---|
   100	| 2.79 GiB | `MzPeakSource::peaks` — `dst.mz.assign` (double, 8 B) |
   101	| 1.39 GiB | `MzPeakSource::peaks` — `dst.intensity.assign` (float, 4 B) |
   102	| 1.39 GiB | `MzPeakSource::peaks` — `dst.ion_mobility.assign` (float, 4 B) |
   103	| 1.95 GiB | Arrow parquet decode batch |
   104	| 1.24 GiB | Arrow parquet decode batch |
   105	
   106	The 2.79 : 1.39 : 1.39 split is exactly the 8:4:4 byte ratio of the three arrays,
   107	which is how the three were identified before symbols confirmed it.
   108	
   109	**So 5.57 GiB of the floor is the `out` block of `SpectrumPeaks`** at
   110	`MzPeakSource.cpp:258-262`, and ~3.2 GiB is Arrow's decode buffers.
   111	
   112	Two independent fixes, both targeting that 5.57 GiB:
   113	
   114	1. **Stop duplicating co-packed frames.** `info_` holds one entry per
   115	   (spectrum, isolation window); a diaPASEF frame co-packs several windows into
   116	   one physical spectrum, so IH1's 32,210 entries cover 17,448 physical spectra —
   117	   **1.85x**. Both entries of a co-packed frame get *byte-identical* peak arrays
   118	   (the comment at line 253 says so outright) and are separated only by their
   119	   mobility bands. Sharing one buffer: 5.57 -> **3.01 GiB**. Lossless. This is why
   120	   the floor is diaPASEF-specific — Astral and SWATH have no co-packing and no
   121	   `ion_mobility` array, so they copy two arrays once instead of three twice.
   122	2. **Reduce `BLOCK`** (`ChromatogramExtractor.cpp:693`, `constexpr 1024`, not
   123	   tunable). Linear in block size; 256 would give another 4x. Tune, not fix — but
   124	   note per-request cost is flat, so the throughput cost should be small.
   125	
   126	Deferred: `mz` as `double` is 8 of the 16 B/peak. float32 quantisation is
   127	0.03-0.06 ppm against a 10 ppm window (169-328x margin), so it is *probably* safe
   128	— but the prior project saw an unexplained 22-peptide delta on a float32 path, so
   129	this must be measured alone, never bundled with a representation change.
   130	
```

```text  (lines 131-161)
   131	## The compaction order — this is the load-bearing part
   132	
   133	Applying these in any other sequence measures nothing, because the third is
   134	*negative* until the first lands.
   135	
   136	**1. Narrow the RT window.** Points per precursor scale directly with it.
   137	Measured at n100k: w600 -> w60 takes RSS 14.00 -> 10.75 GiB, a **3.98 GiB** saving.
   138	This is the dominant lever and it is a prerequisite for step 3. It depends on the
   139	iRT calibration being trustworthy enough to place a narrow window.
   140	
   141	**2. uint8 chromatograms with a per-transition f32 scale.** A clean 4x on what
   142	remains, ~156 -> 39 GiB at 4.26M. Already in the backlog. The scale is per
   143	transition, so two orders of dynamic range fit in 8 bits without a sensitivity
   144	cost.
   145	
   146	**3. The extract/score fusion.** Only pays *after* step 1:
   147	
   148	| window | old | new (fused) | delta |
   149	|---|---:|---:|---:|
   150	| w600 | 14.00 GiB | 14.85 GiB | **+0.85** |
   151	| w60 | 10.75 GiB | 10.56 GiB | **-0.19** |
   152	
   153	Its overhead scales with *live* precursors, and at w600 every precursor is live
   154	(the extractor reports `5330 of 5330 precursors live at once`), so it pays a cost
   155	proportional to the whole library and frees nothing. It also underperforms its own
   156	theory at w60 — ~8.6% of precursors should be live, so it should free ~91% of the
   157	chromatogram term; it freed ~56%. Worth understanding before trusting at scale.
   158	
   159	Steps 1+2 together take the chromatogram term from ~156 GiB to roughly 10 GiB.
   160	**At that point the floor is the dominant term**, which is why it is the right
   161	thing to chase first.
```

## 5. YOUR QUESTIONS (answer each, numbered; file:line on every claim; severity-rank your findings; if something is fine, say so in one line)

Q1. Is the diagnosis of the 657 GB right? Specifically: (i) is "peak RSS = saturated -live_memory_gb budget (425.9 GiB live blocks) + MS1 dense matrix (49.6 GiB) + ~165 GiB unattributed" consistent with the pasted ChromatogramExtractor.cpp (cap derivation :679-700, chunking :701-732, BlockPool :212-255, activate :868-885, emit :886-962, match loop :1139-1250)? (ii) Does anything in the pasted code CONTRADICT it — e.g. is the cap really inverted from MEAN cells so that a full-width live set overshoots the budget, does the BlockPool actually retain blocks across chunks and across the two passes (it is a local at :740 — what is its lifetime relative to pass 1 / finish() / pass 2?), could the `block` vector of decoded SpectrumPeaks (:764, :1104) or the per-chunk index (:964-1057, not pasted in full) be large at this scale, is anything else allocated per-precursor that the diagnoses have not counted? (iii) The "two planes are byte-copies" claim (F03): check :1218 vs :1231 vs :1241 — is ppm_den bit-identical to base under Sum given the `intensity > 0.0f` guard at :1225 and the `at == 0.0f && intensity > 0.0f` counting at :1217? Is im_den identical to base on data where every peak carries 1/K0? What about float accumulation order — same loop, same thread, same order? Any case (Max mode, NaN intensities, negative intensities, `n == 0` dummy_ at :220) where aliasing would change output?

Q2. Which of the fourteen fixes would NOT be output-identical (byte-identical scored TSV), and why? The diagnoses claim F02 (lower budget → more chunks), F03 (alias planes), F04 (parallel sink), F05 (parallel MS1), F06 (instrumentation), F12-row-drop are output-identical and F01, F07, F10, F11, F13, F14 are behaviour-changing. Attack that classification from the pasted code: e.g. does chunking change emission ORDER (:1058-1069 by_lo/by_hi are per chunk; :1259-1268 emit per window in hi order) and does anything in the sink depend on emission order besides Gate C's null (:1248-1282) — e.g. `result.groups.push_back` order (:2655) is re-sorted at :3149-3157, but are there other order-dependent consumers (thread_local rejects_ :1321, counters, the rt_first/rt_last log line :1271-1276)? Does chunking change which spectra are decoded for a precursor at a chunk edge (:1076-1094) in a way that could change a trace? Does the `seen[w]` / activate-one-batch-early logic (:1120-1134) interact with chunk boundaries?

Q3. Which fix is the single best FIRST experiment for (a) memory, (b) wall clock, (c) IDs — and what is the pre-registered pass/fail you would set? Be concrete: name the flag or code change, the one number that decides, and the threshold. Challenge the diagnoses' own picks (memory: F02 then F03; wall: F04 after F06/F08; IDs: F11 stage-0 oracle then F10). If you think F01 (narrow pass 1) dominates everything, say so and say what the risk gate must be.

Q4. What is MISSING from the fourteen that an experienced DIA / HPC engineer would try FIRST? Examples to consider (reject the bad ones): decoys do not need pass-1 extraction at all if pass 1 only harvests anchors from targets (4.96M of 9.92M precursors — is that true from the pasted code, or does pass-1 scoring need decoys for its q-values?); uint8/float16 chromatogram cells with per-transition scale (doc/11:141-144 proposed, never built); building the MS1 matrix only for precursors that have an assignment (:535-540) rather than all 9.92M; overlapping decode of block N+1 with match of block N (:1104 is on the driver); a process-level fix (mmap'd blocks with MADV_DONTNEED on give) instead of BlockPool; running pass 1 on a strided TARGET-only subset; moving the sink onto the pool workers between batches (they are idle at :1259-1268 anyway); replacing the dense MS1 matrix with on-demand MS1 extraction inside the ±window; anything about the mzPeak reader (:363-410) at this scale; anything about `-threads 0` (:1793). Rank your additions against the fourteen.

Q5. Which fixes will actually BITE on real data (IH1 diaPASEF at 9.9M precursors; Astral at 303k spectra) vs are theoretical? For each of F01-F14 give one line: BITES / THEORETICAL / UNDECIDABLE-FROM-BRIEF, with the reason. Flag any fix whose expected win is inconsistent with the measured numbers in §0 (e.g. an F03 claim of −170 GiB when the budget, not the plane count, sets the live bytes; an F04 claim of 20-30x when the match region already runs at ~20% efficiency on the same memory system).

Q6 (short). Any outright ERRORS in the four diagnoses — wrong arithmetic, a cited line that does not say what is claimed, a MEASURED label on something that is actually inferred, or a "closed" hypothesis being re-opened without saying so? List them with file:line.

Format: numbered answers Q1-Q6; within each, findings ranked most-severe first with a severity tag [HIGH]/[MED]/[LOW]/[FINE]; file:line on every source claim using the pasted cat -n numbers; no preamble, no summary of the brief back to me.
