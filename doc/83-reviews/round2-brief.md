# Adversarial review brief #2: OpenDIAlyzer 2 (ODIA) — where we stand on 2026-10-02, and what to do next

## READ THIS FIRST
Review only — do not modify, create or delete any file in the repository. Cite file:line for every source claim
(for pasted text: the `cat -n` numbers shown). Severity-rank findings [HIGH]/[MED]/[LOW]/[FINE]; for each say whether it
would BITE on real data or is THEORETICAL. If something is fine, say so in one line rather than padding. Separate what
you VERIFIED (in pasted/opened source or numbers) from what you INFER. The subject is commit bd73e65 (tag v0.5.0, the
public release). Numbers marked MEASURED come from runs; do not re-derive them from priors.

## 0. THE STANDING (measured; IH1 diaPASEF, DIA-NN's predicted CAM library dn_pred_cam, 4.96M targets + 4.96M decoys)
| | DIA-NN 2.0 | ODIA v0.5.0 best arm (c3_cent_gate_s43) |
|---|---|---|
| wall / CPU | 14 m 25 s / 7 h 38 m | 9 h 31 m / 36 h 32 m |
| cores busy | 31.75 of 48 | 3.83 of 64 |
| peak RSS | 27.3 GiB | 642.2 GiB ("657,612 MB", MiB) |
| IDs at own q<=0.01 | 37,334 (entrapment FDP 1.13%) | 18,995 (true FDP 1.84-2.31%) |
| IDs at matched entrapment count N(62) | 37,334 | 20,040 = 53.7% |

## 1. WHAT THE LAST REVIEW (doc/83, 2026-10-01; you or a sibling model reviewed it) CONCLUDED — summary
- MEMORY: peak RSS is a SATURATED budget, not a footprint. `-live_memory_gb 400` is inverted into a precursor cap from
  mean cells x 4 B x 5 planes; pass 1 at +/-400 s has 66% of the 9.92M library in RT overlap, so 425.9 GiB live.
  Under Aggregate::Sum two of the five planes (ppm_den, im_den) are cell-for-cell copies of the intensity plane
  (170 GiB). MS1 dense matrix 49.6 GiB uncharged. ~165 GiB remainder of unknown owner (BlockPool::reservedPoints()
  never printed; 26 KB blocks on the brk heap; libarrow bundles jemalloc).
- WALL: 71% measured single-threaded by construction — the emit loop runs on the driver after the pool barrier and
  feeds the scorer serially (pass-1 sink 14,637 s + pass-2 4,203 s = 55%); the MS1 trace build 5,590 s (16%) is serial.
  Pass 1 (+/-400 s) is 60% of wall and exists only to harvest calibration anchors.
- IDS: admission 97.7%, correct peak in the <=3-candidate list 92%, ranked inside depth 37,334 64.1%: the loss is
  RANKING WITH THE RIGHT PEAK IN HAND — fragment co-elution failure at the correct apex (per-fragment reference corr
  0.85 recovered vs 0.37 buried; AUC 0.93 within DIA-NN's own quantity deciles). Closed/refuted (do not re-propose
  without new evidence): m/z window narrowing, IM half-window changes, affine IM pre-correction, Gate C relaxation,
  candidate generation, -max_candidates, training knobs k=1/12/24 and c6, XGBoost import, DIA-NN score injection,
  congestion, trace tensors, MS1-as-intensity, deterministic Gate C null (F07), fine-tuned library through ODIA (F10),
  retained-group re-seek (F13), decoy-tail truncation (F14).
- Codex found (verified): `thread_local RegisteredRejects rejects_` registers raw pointers in a global registry and
  never unregisters (PeakGroupScorer.cpp:1310-1332) — a dangling-pointer bug for any parallel sink. Codex proposed a
  calibration-only pass 1 (stratified target+decoy panel) and a byte-weighted chunker.
- The doc/83 PLAN (ordered): (1) -live_memory_gb 150 arm; (2) -rt_window_pass1 200 and 110 arms; (3) FRAGVEC at-depth
  oracle read; (4) instrumented binary; (5) frame-parallel MS1 build; (6) plane aliasing at a pinned cap + MS1
  uncovered-row drop + byte-weighted chunker; (7) parallel sink by ordered speculation after fixing the registry;
  (8) calibration-only pass 1.

## 2. WHAT HAS HAPPENED SINCE (2026-10-02)
A. Plan steps 1-2 are RUNNING on dax (results ~09:00-12:30 CEST), each pre-registered in its script header (pasted
   in section 4): sb1_budget150 (-live_memory_gb 150; raw-TSV cmp vs c3 is the identity gate; Peak <= 420,000 MB),
   sb2_p1w200 and sb3_p1w110 (-rt_window_pass1 200/110; gates: EXTERNAL RT coverage of DIA-NN's 37,334 apexes from the
   -out_lib calibrated library within 0.5 pt of sb1's; IM anchors >= 120 per charge; region median N(e) FDP 2-10%
   >= -1% vs c3). A /proc sampler rides every arm (VmHWM phase, driver vs worker faults, NUMA counters; dax is
   2-socket, numa_balancing=1).
B. Plan step 3 is DONE (pre-registered; definitions lifted from the record with sha256 checks). Recall of DIA-NN's
   37,334 inside ODIA's top 37,334 (identity recall, two-seed mean; reference-labelled role-matched GBT fits on the
   arm x1_fragvec_r1; MEASURED):
       native ODIA DScore                        64.086%
       A   oracle-fit on the arm's 42 columns    67.655%   (+3.57 over native)
       Ris A + 6 FRAGVEC scalars                 69.547%   (+1.891 over A)  -> DO NOT WIRE (bar +2.0)
       Ri  A + 73 fitted FRAGVEC columns         70.455%   (+2.799 over A)  -> CONDITIONAL (FUND bar +3.0)
       Aperm permuted-column control             67.658%   (+0.003)         MDE 0.072
   Ri reaches depth with ~224 fewer entrapments (1,188 vs 1,412). Ri − A grows with depth (+0.58 @20k, +2.30 @30k,
   +2.80 @37,334, +2.98 @45k); most of it is the per-rank columns (Riv +2.58). On 8-s concordant recall Ri reads
   +3.024 (FUND by 0.024, inside that measure's MDE 0.123); the user has DECIDED: identity recall, CONDITIONAL, WIRE IT.
   CONDITIONAL means: wire the block into the engine as sub-scores behind a flag (off = byte-identical), then the
   engine arms must pass BOTH baselines (b2_cap1, 18,072 at q<=0.01; w1_ctl70, 16,965) at matched entrapment budget.
C. Exploratory implementations are starting now, each on its own branch from v0.5.0, fixture-gated:
   (i) FRAGVEC sub-scores (-fragvec_scores, final pass); (ii) instrumented binary (rusage/VmHWM per stage,
   reservedPoints per chunk, plane-identity assert); (iii) frame-parallel Ms1Traces::build; (iv) ppm_den/im_den
   aliasing at a pinned cap + byte-weighted chunker + never-read MS1 row drop; (v) rejects_ registry fix + a
   parallel-sink prototype by ordered speculation.

## 3. YOUR QUESTIONS (numbered; answer each)
Q1 THE STANDING. Attack the conclusions in sections 1-2. What is wrong, overclaimed, or missing? In particular: is
   "the loss is ranking with the right peak in hand" still the right frame, given that an oracle refit on the SAME 42
   columns beats the shipped scorer by +3.57 pts at depth — i.e. the native semi-supervised training (lda.h: positive
   selection by q on its own decoy null, folds, iterations) leaves more on the table than the whole FRAGVEC block adds?
Q2 THE ORACLE-vs-NATIVE GAP (+3.57). What mechanisms in the pasted training code (lda.h, gbt.h, the sub-score
   construction) could produce it, ranked; which are testable cheaply (no new extraction) and how? The record says the
   training RECIPE knobs are a null (k=1/12/24, c6) and that the scorer flips between a compact and a saturated
   regime (max DScore ~20 vs 75-165) with fold seed alone (16,965 -> 4,105 at q<=0.01, no byte changed); the leaf cap
   (-gbt_max_delta_step) restores one seed but failed a second library. Is there a structural fix (e.g. a
   non-self-training target, monotone constraints, calibrated per-fold scores, a different negative set, ensembling
   folds/seeds, early stopping on the decoy null) that the record has not tried? Be concrete.
Q3 THE FRAGVEC WIRING. The 78-column R1 block (six 12-long vectors by LIBRARY-INTENSITY rank: per-fragment 5-cycle
   area, share, log-ratio vs library share, at-apex flag, missingness, + six scalars) is today an OUTPUT-ONLY export
   computed only on the final pass (OpenDIAlyzer.cpp:4431-4434; PeakGroupScorer.cpp ~2488-2654). Wiring it as
   sub-scores: what can go wrong? (a) decoy-construction artefacts: the GBT can learn to separate decoys from targets
   by how decoy library intensities/fragments were generated rather than by evidence (this project caught exactly
   that once: a trace model "learned the decoy artifact"); which of the 78 columns are at risk, and how to test it
   (entrapment FDP is the arbiter; what else?); (b) the semi-supervised loop re-landing in the saturated regime;
   (c) final-pass-only vs both passes; (d) NaN/missingness semantics and the constant-column guard; (e) memory
   (73 doubles x 21.7M groups); (f) any reason the oracle read (+2.80 at depth, reference-labelled) will NOT transfer
   to native self-training. Propose the acceptance arms and controls you would require.
Q4 THE NEXT ROUND FOR IDENTIFICATIONS. Beyond FRAGVEC, what are the 5-8 most promising next solutions to close the
   46-point gap, ranked by expected IDs per week of work, each with: mechanism, the measurement that would justify or
   kill it cheaply, and the risk. Consider what DIA-NN/AlphaDIA/Alpha-XIC/DreamDIA do structurally (per-fragment
   evidence vectors vs a data-driven best-fragment reference; tiny FFN ensembles trained per run; co-located
   candidate competition; learned co-elution scores; interference-robust recomputation; MS1 isotope channels; per-run
   library refinement) — and what this record has already refuted. No re-proposals of closed items without new
   evidence.
Q5 THE RUNNING ARMS. Attack the three registrations pasted in section 4: gates that cannot fail, confounds,
   thresholds that are wrong, what the result will and will not license. Anything you would add NOW (they are
   running; a read-side addition is still possible)?
Q6 THE EXPLORATORY IMPLEMENTATIONS (section 2C). Ordering, missing items, design traps, and the acceptance gate each
   needs. Which of them would you NOT do, and why? What is missing that should be started now?
Q7 WHAT TO STOP. What is this project doing that wastes node-days or produces numbers that cannot decide anything?
Format: numbered answers Q1-Q7, findings ranked within each, file:line on source claims, no preamble.

## 4. THE THREE RUNNING REGISTRATIONS (script headers, verbatim)
### shared/libv2/run_sb1.sh
```
#!/bin/bash
# ==============================================================================================================
# sb1_budget150 -- doc/83 PLAN STEP 1: the memory budget is the memory. ZERO CODE.
#
# WHAT IT IS.  run_c3.sh's command (c3_cent_gate_s43, the reference arm of doc/83) character for character, on
# the same binary (/scratch/kohlbach/odia/build_f1q/OpenDIAlyzer, sha256 a568ddf0, v1.15), same seed 43, with
# ONE behavioural change:  -live_memory_gb 400 -> 150.  Two OUTPUT-ONLY files added:
#   -out_rt_map $S/sb1_budget150.trafoXML  -- a pure store of the map the run used (OpenDIAlyzer.cpp:3077-3090 at
#      bd73e65: TransformationXMLFile().store + a log line; feeds nothing back);
#   -out_lib $S/sb1_budget150.lib.parquet  -- the CALIBRATED library (RT in run seconds, map + accepted refinement),
#      written after the scored TSV (:5619-5649). It is the 400-width calibration reference for the step-2 external
#      RT-coverage gate (analysis77/rt_cov_ext.py) -- valid as c3's stand-in when this arm's "fitted the
#      retention-time map" and "rt refine" lines equal c3's (they must, if pass 1 is unchanged).
#   The -out_lib write happens AFTER scoring; if the sampler dates the high-water mark to that write ("wrote scored"
#   stage), M2 is read from the last sample before it and the write's own cost is reported separately.
#
# WHY.  doc/83 2.2: peak RSS 657,612 MiB is the 400 GiB budget saturating -- the cap is DERIVED from the budget
# (ChromatogramExtractor.cpp:679-699), pass 1 at +/-400 s has 5,704,152 precursors in RT overlap against a cap of
# 3,324,888, so the pool fills to the cap in each of 3 chunks (c3 log: 425.944 GiB live). Lower the budget and the
# live set must fall; chromatogram VALUES are chunk-invariant by construction (:759-761); emission ORDER is not
# (cross-chunk drain + unstable sorts :1064-1068) -- and on dn_pred_cam every order-sensitive consumer is believed
# inert (Gate C 0/0, mass-anchor cap not reached). This arm tests that belief and the memory claim at once.
#
# ============================== PRE-REGISTERED READ (this header IS the registration; written 2026-10-02 before
# ============================== launch, never amended afterwards)
#  GATE A (else SUSPENDED): sha a568ddf0; both completion lines ("wrote scored peak groups", "OpenDIAlyzer took");
#    the budget line reads "budget 150 GiB / 129176 B per live precursor (5 planes) -> cap 1246833" (+/-1).
#  CLAUSE M1 DECISIVE -- IDENTITY: raw `cmp` of the UNMODIFIED scored TSV against c3_cent_gate_s43.tsv is clean
#    (codex, doc/83 5 row 19: sorting before cmp hides ordering regressions).
#  CLAUSE M2 DECISIVE -- MEMORY: footer "Peak Memory Usage" <= 420,000 MB (c3: 657,612).
#  CLAUSE M3 SUPPORTING -- COST: pass-1 chunks ~7 (6-9); pass-1 decode <= 1,300 s (c3 408); wall <= 1.10 x c3's
#    34,303 s = 37,733 s. Wall is read with the concurrency noted (sb2/sb3 run beside it; c3 also ran with 2 others).
#  REPORTED (no threshold): Peak - 1.065 x 153,600 MiB against c3's remainder 221,445 MiB (fixed remainder predicts
#    ~385k at 150; < 300k would mean the remainder is budget-proportional); the HWM phase from the sampler
#    (the stage column at the first sample whose VmHWM is within 1% of the final); driver vs worker minflt/stime;
#    numa_hint_faults delta over the run.
#  VERDICT. PASS = A, M1, M2 -> 150 GiB is a free 1/3 cut and every later arm runs 5-6 per node instead of 3.
#    M1 FAIL -> run analysis77/tsv_coldiff.py: ORDER or DSCORE<=1e-9 = order/reduction effect, reported, M2 still
#    read; ROWS or CONTENT = a chunk-edge defect -> STOP budget work, report; next arm is c3 at 400 with -threads 48
#    to clear the cohort_plain thread confound (doc/83 4 F02). M2 FAIL (> 500,000) -> the remainder is fixed, F06's
#    instrumentation is the only memory lever.
#  PREDICTION: Peak 360,000-420,000 MB; 7 chunks; decode ~1,000 s; wall 34,300-37,000 s; cmp clean.
# ==============================================================================================================
set -uo pipefail
```
### shared/libv2/run_sbw.sh
```
#!/bin/bash
# ==============================================================================================================
# run_sbw.sh <200|110> -- doc/83 PLAN STEP 2: is pass 1's +/-400 s window load-bearing?  ZERO CODE.
#   200 -> sb2_p1w200      110 -> sb3_p1w110
#
# WHAT IT IS.  run_c3.sh's command (c3_cent_gate_s43) character for character, same binary a568ddf0, same seed 43,
# with ONE behavioural change:  -rt_window_pass1 400 -> 200 (sb2) or 110 (sb3);  plus the two OUTPUT-ONLY files
# -out_rt_map and -out_lib (see run_sb1.sh for why each is output-only).  -live_memory_gb stays 400: at these widths
# pass 1 is predicted RT-overlap-bound (one chunk), so the budget is not a second factor.
#
# WHY.  doc/83 2.3: pass 1 is 20,705 s of c3's 34,303 s wall (60.4%) and the only reason the budget binds (426 GiB
# live at +/-400 s); bytes and sink positions are linear in width (129,176 vs 36,949 B/precursor, 1.70 vs 0.49 ms
# per emitted precursor at 400 vs 110 s). It exists to harvest anchors: 3,869,789 RT, 15,484 1/K0, 1,211,024
# centred entries. The question is whether a narrower pass 1 harvests an equally good calibration. Codex (doc/83
# 5 row 20) prefers 200 first, kimi 110; both run.
#
# ============================== PRE-REGISTERED READ (this header IS the registration; written 2026-10-02 before
# ============================== launch, never amended afterwards)
#  GATE A (else SUSPENDED): sha a568ddf0; both completion lines; the log reads "rt_window_seconds=<W>" for pass 1.
#  GATE A3 (else the arm is TWO-FACTOR, reported not decided): "pass 2 extraction window 110 s" unchanged -- the
#    window is min(cap 110, max(floor 20, 2 x p95)) (OpenDIAlyzer.cpp:3146-3152), so a narrower pass 1 that shrinks
#    the anchor p95 below 55 s would move pass 2 too.
#  CLAUSE C1 DECISIVE -- EXTERNAL RT COVERAGE (not own anchors: the anchor p95 is window-bounded and cannot fail;
#    refuters + codex, doc/83 5 row 10): analysis77/rt_cov_ext.py coverage of DIA-NN's CAM q<=0.01 precursors at
#    +/-60 s AND +/-110 s within 0.5 points of sb1_budget150's (the 400-width calibration; c3's stand-in) -- join
#    >= 95% required, else the clause is unread.
#  CLAUSE C2 DECISIVE -- THE OTHER ANCHORS SURVIVE: "ion-mobility calibration: GATE PASSED -- 4 of 4 charges
#    corrected"; every charge >= 120 1/K0 anchors (c3: z1 162, z2 7198, z3 2215, z4 297); the centring line present
#    with >= 605,512 centred entries (50% of c3's 1,211,024).
#  CLAUSE C3 DECISIVE -- IDS AT MATCHED ENTRAPMENT BUDGET: region median N(e) over FDP 2-10% vs c3_cent_gate_s43
#    >= -1.0% (stair_read.py med%); q<=0.01 >= 17,096 (0.9 x 18,995). Read beside w1_ctl70 (retained best; context,
#    no clause). Regime label reported (c3: SATURATED); a flip is reported, not decided.
#  CLAUSE C4 SUPPORTING -- COST:
#    sb2 (200): pass-1 "1 chunks" (RT-bound), pass-1 live <= 200 GiB, Peak <= 400,000 MB, wall <= 25,200 s (7.0 h).
#    sb3 (110): pass-1 "1 chunks", pass-1 live <= 70 GiB, pass-1 sink <= 5,000 s, Peak <= 300,000 MB,
#               wall <= 23,400 s (6.5 h).
#  VERDICT. PASS = A, A3, C1, C2, C3 -> that width is a DEFAULT CANDIDATE, not a default: a default is an experiment
#    (width-2.0 regression, memory) -- it needs a second seed and a human_v2 arm first. FAIL C1 or C3 -> the window
#    is load-bearing at that width; if sb2 passes and sb3 fails, 200 is the floor. C2 fails alone -> report which
#    anchor family starved (the step-8 calibration-only pass 1 inherits exactly that question).
#  PREDICTIONS.  sb2: C1 within 0.3 pt; region -0.5..+1.0%; Peak 250-350k MB; wall 22,000-25,000 s.
#                sb3: C1 drop 0-1.0 pt; region -2.0..+1.0%; Peak 170-230k MB; wall 18,000-21,000 s.
# ==============================================================================================================
set -uo pipefail
```

## 5. THE FRAGVEC READ — registration and result (verbatim)
### analysis77/pick/sb_fragvec_depth_registration.md
```
# sb_fragvec_depth — doc/83 plan step 3: does the FRAGVEC R1 block carry AT-DEPTH ranking information?

Registered 2026-10-02 ~03:05 CEST, BEFORE any at-depth number for these score files was computed.
Never amended afterwards; anything learned while computing goes into the result file, not here.

## The question
The FRAGVEC acceptance pair passed at full scale on x1_fragvec_r1 (wf_v69: matched-score dAUC +0.0825/+0.0842,
common-depth +298/+232 concordant at N(62), 49/43 entrapments vs 62). Every number on record is a
matched-budget or common-depth delta; nobody has read RECALL OF DIA-NN'S SET AT DEPTH 37,334 for the oracle-fitted
score sets. doc/83 §4 F11: this is a re-read of existing per-row scores, minutes of compute, and it decides whether
the engine wiring (3-5 full arms) is worth funding.

## Inputs (dax, node-local)
- `/scratch/kohlbach/odia2x2/wf_v63/x1_fragvec_r1_{A,Ris,Ri,Alog,Riv,Aperm}_N3_seed{1,2}.npz` — per-row scores of
  reference-labelled, role-matched (N3), family-fold GBT fits: A = the arm's 42 var_* columns; Ris = 42 + 6 FRAGVEC
  scalars; Ri = 42 + 72 FRAGVEC columns; Alog/Riv = variants; Aperm = permuted-column control.
- `/scratch/kohlbach/odia2x2/x1_fragvec_r1.tsv` (the arm's scored TSV; native DScore) and `.fragvec.tsv`.
- Reference: whichever DIA-NN set wf_v63/wf_v69 used — to be ESTABLISHED from their code, not assumed. If it is the
  superseded CAM-free 33,330 set, the read is done on BOTH that set at depth 33,330 (to reproduce the record) AND the
  CAM 37,334 set (`dn_ih1_cam.parquet`) at depth 37,334 (the decision), with the UniMod:4 alias join.

## Definitions — LIFTED, never re-implemented
Depth, the per-precursor ranking quantity (max score per precursor, never min q), the id join and SEL/conc are
lifted from the record's own instruments (wf_v61_refute_recompute.py `depth`/`read_ids`/`at`, as wf_v69 does, with
the same sha256 checks). The native comparator is x1_fragvec_r1's OWN DScore ranked the same way (in-arm), not the
64.1% of a different arm.

## Sanity gates (any failure → STOP, report, no verdict)
- S1 JOIN: the npz scores reproduce wf_v69's reported common-depth figures for this arm at N(62) (concordant gains
  and entrapment counts) to the digit for at least A and Ri, both seeds.
- S2 NATIVE: the in-arm native recall at depth is computed by the same code and reported.
- S3 A >= native − 0.5 pt (an oracle fit on the arm's own 42 columns that ranks below the arm's own scorer means the
  read is broken, not that the columns are weak).

## The decision (per seed, then the two-seed mean; MDE = the larger |seed1 − seed2| over A, Ris, Ri)
- FUND stage 1 (engine wiring):      Ri ≥ A + 3.0 pts  AND  (Ri − A) > MDE  AND  Aperm ≤ A + MDE
- CONDITIONAL (fund only if stage-1 arms pass both baselines): A + 2.0 < Ri < A + 3.0, same MDE and Aperm conditions
- DO NOT WIRE:                       Ri ≤ A + 2.0, or Aperm gains as much as Ri (the gain is capacity, not information)
  A DO-NOT-WIRE licenses ONLY "the R1 areas/missingness add no at-depth recall" — NOT "the evidence is not in the
  traces" (doc/83 §5 row 5; the R1 block has no trace coherence). The next offline step is then rung (ii), the
  coherence columns from the sealed builder (wf_v33_fragvec_contract.md:159-176).
- Ris (42 + 6 scalars) is reported beside Ri: it decides the SIX-COLUMN scoping question that is UNDECIDED in the
  record (share 0.54 vs a 0.60 bar), by the same rule with Ris in place of Ri.
- Alog, Riv: reported, no clause.

## Prediction (registered)
native ≈ 60-64%; A ≈ native + 0.5..+3; Ri − A ≈ +1.5..+3.5 pts (CONDITIONAL most likely); Aperm − A within ±0.5.
```
### analysis77/pick/sb_fragvec_depth_result.txt (sections 6-8)
```
6. GATE S2 (NATIVE) -- the arm's OWN DScore, ranked by the SAME code
   (a) wf_v69 l.1899 verbatim: Ranking('incumbent', g_selds, g_selord)  -- the group table's retained selection
   (b) -out's DScore column streamed and passed through the SAME group_winner the fitted sets use (max per group,
       first in file order) -> Ranking('native_tsv', ...).  (a) and (b) must agree.
======================================================================================================================================================
  PASS (a) incumbent N(62) / N(73) / N_class(60) / acc / conc at N(62) = 17,272 / 18,563 / 17,735 / 16,477 / 15,996 [M] == record 17,272 / 18,563 / 17,735 / 16,477 / 15,996 [Q wf_v69 s.8]
  PASS -out streamed: 21,849,969 data rows == the group table's 21,849,969 (16 s) [M]
  PASS ROW ALIGNMENT: -out's Decoy column == the group table's decoy flag on every one of the rows [M]
  -out rows with a non-finite DScore: 0 (set to -inf, as wf_v69 l.919 does for sel_DScore) [M]
  PASS (b) max DScore per target group from -out == the retained sel_DScore on all 3,762,850 target groups (mismatches: target 0, all groups 0) [M]
  PASS (b) the max-DScore row (first in file order) == the retained sel_ordinal on 3,762,850 target groups (mismatches 0) -- a mismatch would move only the 8-s localisation (conc), never acc [M]
  PASS (a) == (b) at D = 37,334: acc/conc/E (23926, 23144, 1691) vs (23926, 23144, 1691), and the two entrapment staircases are identical [M]
  PASS GATE S2: native recall at depth 37,334 = 23,926 / 37,334 = 64.086 % (conc 23,144 = 61.992 %), identical by both routes; incumbent reproduces the record [M]
    for scale only, never a comparator: doc/82:89 quotes 64.1 % @ depth 37,334 for "ODIA + P1 port, same library" [Q] -- the arm whose -out this one is byte-identical to (gate A1)

======================================================================================================================================================
7. RECALL OF THE 37,334-ID CAM REFERENCE AT DEPTH D = 37,334 (top D target precursors by max score; lifted at())
======================================================================================================================================================
  PASS native: the lifted wf_v61 at() == wf_v69's independent numpy path at D = 37,334 (23,926/23,144)
  PASS A s1: the lifted wf_v61 at() == wf_v69's independent numpy path at D = 37,334 (25,272/24,486)
  PASS A s2: the lifted wf_v61 at() == wf_v69's independent numpy path at D = 37,334 (25,245/24,456)
  PASS Alog s1: the lifted wf_v61 at() == wf_v69's independent numpy path at D = 37,334 (25,674/24,888)
  PASS Alog s2: the lifted wf_v61 at() == wf_v69's independent numpy path at D = 37,334 (25,669/24,881)
  PASS Ris s1: the lifted wf_v61 at() == wf_v69's independent numpy path at D = 37,334 (25,978/25,239)
  PASS Ris s2: the lifted wf_v61 at() == wf_v69's independent numpy path at D = 37,334 (25,951/25,193)
  PASS Riv s1: the lifted wf_v61 at() == wf_v69's independent numpy path at D = 37,334 (26,225/25,542)
  PASS Riv s2: the lifted wf_v61 at() == wf_v69's independent numpy path at D = 37,334 (26,219/25,523)
  PASS Ri s1: the lifted wf_v61 at() == wf_v69's independent numpy path at D = 37,334 (26,305/25,611)
  PASS Ri s2: the lifted wf_v61 at() == wf_v69's independent numpy path at D = 37,334 (26,302/25,589)
  PASS Aperm s1: the lifted wf_v61 at() == wf_v69's independent numpy path at D = 37,334 (25,281/24,493)
  PASS Aperm s2: the lifted wf_v61 at() == wf_v69's independent numpy path at D = 37,334 (25,238/24,449)

  ranking   seed    acc@D  RECALL %   conc@D conc-rec %   E(D)  FDP^(D) | rec-A pts  rec-nat |  crec-A crec-nat
  native      --   23,926    64.086   23,144     61.992   1691  30.932% |        --       -- |      --       --   [M]
  A            1   25,272    67.692   24,486     65.586 1415.0  25.883% |    +0.000   +3.605 |  +0.000   +3.595   [M]
  A            2   25,245    67.619   24,456     65.506 1408.0  25.755% |    +0.000   +3.533 |  +0.000   +3.514   [M]
  A         mean 25,258.5    67.655 24,471.0     65.546 1411.5  25.819% |    +0.000   +3.569 |  +0.000   +3.554   [M]

  Alog         1   25,674    68.768   24,888     66.663 1282.0  23.451% |    +1.077   +4.682 |  +1.077   +4.671   [M]
  Alog         2   25,669    68.755   24,881     66.644 1280.0  23.414% |    +1.136   +4.669 |  +1.138   +4.653   [M]
  Alog      mean 25,671.5    68.762 24,884.5     66.654 1281.0  23.432% |    +1.106   +4.675 |  +1.108   +4.662   [M]

  Ris          1   25,978    69.583   25,239     67.603 1236.0  22.609% |    +1.891   +5.496 |  +2.017   +5.612   [M]
  Ris          2   25,951    69.510   25,193     67.480 1246.0  22.792% |    +1.891   +5.424 |  +1.974   +5.488   [M]
  Ris       mean 25,964.5    69.547 25,216.0     67.542 1241.0  22.701% |    +1.891   +5.460 |  +1.996   +5.550   [M]

  Riv          1   26,225    70.244   25,542     68.415 1195.0  21.859% |    +2.553   +6.158 |  +2.829   +6.423   [M]
  Riv          2   26,219    70.228   25,523     68.364 1206.0  22.060% |    +2.609   +6.142 |  +2.858   +6.372   [M]
  Riv       mean 26,222.0    70.236 25,532.5     68.389 1200.5  21.960% |    +2.581   +6.150 |  +2.843   +6.398   [M]

  Ri           1   26,305    70.459   25,611     68.600 1202.0  21.987% |    +2.767   +6.372 |  +3.013   +6.608   [M]
  Ri           2   26,302    70.451   25,589     68.541 1173.0  21.457% |    +2.831   +6.364 |  +3.035   +6.549   [M]
  Ri        mean 26,303.5    70.455 25,600.0     68.570 1187.5  21.722% |    +2.799   +6.368 |  +3.024   +6.578   [M]

  Aperm        1   25,281    67.716   24,493     65.605 1413.0  25.847% |    +0.024   +3.629 |  +0.019   +3.613   [M]
  Aperm        2   25,238    67.601   24,449     65.487 1413.0  25.847% |    -0.019   +3.514 |  -0.019   +3.495   [M]
  Aperm     mean 25,259.5    67.658 24,471.0     65.546 1413.0  25.847% |    +0.003   +3.572 |  +0.000   +3.554   [M]

  informational, NOT registered, no clause: two-seed-mean recall delta vs A (pts, acc) at other depths [M]
         depth    20,000    30,000    33,330    37,334    45,000    60,000
        native    49.973    60.698    62.407    64.086    66.334    69.151   (native level, %)
             A    50.497    63.774    65.858    67.655    69.994    72.785   (A level, %)
          Alog    +0.181    +0.924    +1.071    +1.106    +1.266    +1.307
           Ris    +0.343    +1.453    +1.687    +1.891    +1.966    +2.066
           Riv    +0.589    +2.179    +2.405    +2.581    +2.772    +2.818
            Ri    +0.583    +2.302    +2.563    +2.799    +2.983    +2.949
         Aperm    -0.012    -0.003    +0.004    +0.003    +0.015    +0.047
  informational, NOT registered, no clause: (Ris - A) / (Ri - A) at D = 37,334, two-seed mean: acc 0.676, conc 0.660 [M].  This is NOT the record's SHARE statistic (wf_v69 prereg, 0.540 vs a 0.60 bar) and decides nothing; the Ris clause below does.

======================================================================================================================================================
8. GATE S3 -- A >= native - 0.5 pt (an oracle on the arm's own 42 columns ranking below the arm's own scorer means the
   read is broken, not that the columns are weak).  Applied per seed on the PRIMARY measure; conc shown beside.
======================================================================================================================================================
  A    1: recall 67.692 vs native 64.086 - 0.5 = 63.586  -> PASS   (conc: 65.586 vs 61.492 -> pass) [M]
  A    2: recall 67.619 vs native 64.086 - 0.5 = 63.586  -> PASS   (conc: 65.506 vs 61.492 -> pass) [M]
  A mean: recall 67.655 vs native 64.086 - 0.5 = 63.586  -> PASS   (conc: 65.546 vs 61.492 -> pass) [M]
  PASS GATE S3: A - native = +3.605 / +3.533 pts (seed 1 / 2), floor -0.5

======================================================================================================================================================
9. THE REGISTERED DECISION, APPLIED MECHANICALLY (per seed, then the two-seed mean; the mean decides)
   MDE = the larger |seed1 - seed2| over A, Ris, Ri.  gain = X - A,  pgain = Aperm - A  (pts of recall)
   FUND        X >= A + 3.0  AND  gain > MDE  AND  Aperm <= A + MDE
   CONDITIONAL A + 2.0 < X < A + 3.0, same MDE and Aperm conditions (fund only if stage-1 arms pass both baselines)
   DO NOT WIRE X <= A + 2.0, or Aperm gains as much as X (the gain is capacity, not information)
   anything else is printed as NO VERDICT BY THE LETTER, with the clause that failed (the rule assigns none)
   MEASURE -- A JUDGEMENT CALL, declared in this script's docstring BEFORE any at-depth number was computed (first
   execution, script sha256 c9bd7e9497399da8...): the registration says "recall of DIA-NN's set at depth 37,334" and does
   not choose between acc and conc.  PRIMARY = acc/37,334 ("in ODIA's top D"), because every bar the registration and
   doc/83 anchor to is that quantity: doc/82:89's 64.1 % @ depth 37,334 (reproduced above as 64.086 % by acc; conc
   gives 61.99 %), doc/83 F11's "A > 64.1 %", and the registered "native ~ 60-64 %".  SECONDARY = conc/37,334 (the
   record's canonical concordance, lifted SEL), read by the same rule and printed beside; it does not decide.
======================================================================================================================================================

  --- PRIMARY (recall = acc/37,334):  MDE = 0.072 pts  (|s1-s2|: A 0.072, Ris 0.072, Ri 0.008) [M]
      Ri      1: A 67.692  Ri 70.459  Aperm 67.716  ->  CONDITIONAL                (gain +2.767 > MDE 0.072; pgain +0.024 <= MDE) [M]
      Ri      2: A 67.619  Ri 70.451  Aperm 67.601  ->  CONDITIONAL                (gain +2.831 > MDE 0.072; pgain -0.019 <= MDE) [M]
      Ri   mean: A 67.655  Ri 70.455  Aperm 67.658  ->  CONDITIONAL                (gain +2.799 > MDE 0.072; pgain +0.003 <= MDE) [M]
      Ris     1: A 67.692  Ris 69.583  Aperm 67.716  ->  DO NOT WIRE                (gain +1.891 <= +2.0) [M]
      Ris     2: A 67.619  Ris 69.510  Aperm 67.601  ->  DO NOT WIRE                (gain +1.891 <= +2.0) [M]
      Ris  mean: A 67.655  Ris 69.547  Aperm 67.658  ->  DO NOT WIRE                (gain +1.891 <= +2.0) [M]
      Riv  mean: gain +2.581 pts (reported, no clause) [M]
      Alog mean: gain +1.106 pts (reported, no clause) [M]

  --- SECONDARY (conc-recall = conc/37,334):  MDE = 0.123 pts  (|s1-s2|: A 0.080, Ris 0.123, Ri 0.059) [M]
      Ri      1: A 65.586  Ri 68.600  Aperm 65.605  ->  FUND                       (gain +3.013 > MDE 0.123; pgain +0.019 <= MDE) [M]
      Ri      2: A 65.506  Ri 68.541  Aperm 65.487  ->  FUND                       (gain +3.035 > MDE 0.123; pgain -0.019 <= MDE) [M]
      Ri   mean: A 65.546  Ri 68.570  Aperm 65.546  ->  FUND                       (gain +3.024 > MDE 0.123; pgain +0.000 <= MDE) [M]
      Ris     1: A 65.586  Ris 67.603  Aperm 65.605  ->  CONDITIONAL                (gain +2.017 > MDE 0.123; pgain +0.019 <= MDE) [M]
      Ris     2: A 65.506  Ris 67.480  Aperm 65.487  ->  DO NOT WIRE                (gain +1.974 <= +2.0) [M]
      Ris  mean: A 65.546  Ris 67.542  Aperm 65.546  ->  DO NOT WIRE                (gain +1.996 <= +2.0) [M]
      Riv  mean: gain +2.843 pts (reported, no clause) [M]
      Alog mean: gain +1.108 pts (reported, no clause) [M]

  ------------------------------------------------------------------------------------------------------------------------
```

## 6. doc/83 §2.1 (IDs) and §6 (plan), verbatim
```
### 2.1 IDs — the 46.3-point gap, ranked by share

Decomposition basis (samelib3 chain, same library as DIA-NN; stages MEASURED, the split of the gap into points INFERRED by carrying stage fractions to 37,334 per `D/82:103-107`):

| stage | of DIA-NN's set | source |
|---|---|---|
| admission (reaches scorer) | 97.7% | `D/77:22-24` |
| correct peak in the ≤3-candidate list | 92.0% (±20 s) | `D/77:73-81`; `M/odia-pick-losses.md:13-16` |
| ranked inside depth 37,334 | 64.0-64.1% | `D/79:25`; `D/82:89,103-107` |
| at N(62) | 53.7% | `M/odia-vs-diann-measured.md:23` |

1. **Ranking with the correct peak in hand — ≈33 pts (~12,400 precursors). Mechanism MEASURED: fragment co-elution failure at the correct apex.** At DIA-NN's own apex, per-fragment reference correlation 0.85 (recovered) vs 0.37 (buried), AUC 0.957; MS1↔MS2 corr 0.78 vs 0.16; still AUC 0.932 inside DIA-NN's own quantity deciles, so not abundance (`L/analysis77/forensics/FORENSICS_VERDICT.md:19-25,59-70`). 12/12 fragments present with ≥4 points in both classes — presence refuted (`FORENSICS:17-18,524-525`; `D/45:458-499`). Interference is distributed (apex MAD 20.8 vs 2.8 s), "deconvolution, not selection" (`D/45:512-552`). Family ablation: only library-agreement (−7.42) and co-elution (−4.65) carry signal (`M/odia-three-tools-three-failures.md:20-23`). A DIA-NN-labelled oracle on ODIA's 22 sub-scores plateaus at 65.2% (22 columns, depth 33,330 against the superseded CAM-free reference, carried over per `D/82:103-107`), invariant across six recipe perturbations — a classifier-side ceiling on the shipped columns (`D/77:37-62`; `M/odia-feature-plateau-65pct.md:8-15`); the full-scale 42-column probe is a lower bound that does NOT attribute the residual to missing feature content (`M/odia-feature-plateau-65pct.md:34`). At the decision point the 42 columns reach matched-score AUC 0.654 pooled within DScore (`M/odia-matched-score-is-the-test.md:10`; 0.811/0.808 on the (N,2N] window under wf_v69, `M/odia-fragvec-passes-six-columns.md:54`) with the co-elution family 0.59-0.61 (`M/odia-matched-score-is-the-test.md:29-33`). Caveat the forensics records: A/C classes are defined by ODIA's DScore, so levels partly decompose ODIA's own decision; the direction rests on own-decoy nulls and §6 measurements (`FORENSICS:43-48`).
   - Mediator, share with #1 UNRESOLVED: multi-candidate precursors (≥2 candidates → 85% missed, 0.915 adjusted, n=16,349; LOBO "mediator, not cause", `FORENSICS:35-36,174,188-189`). The deciding lane-join is partly on record: on b2_cap1 at N(62) the "proposed within 8 s but outranked" state is 1,834 ids = 9.2% of the loss, the dominant cell is "near row SELECTED but ranked deep" 13,349 = 67.3% (`L/analysis77/LOSS_REVIEW_2026-09-08.md:107-119,246`) — a re-seek cannot act on 67%.
   - Exogenous dose-response (the only clean one): library 1/K0 residual >0.030 → 80% missed at matched abundance AND matched DIA-NN evidence, n=2,967 ≈ 8.9% of the set (`FORENSICS:34-35,186-187,288-306`). On dn_pred_cam the column is better than on APD libraries: median |lib−obs| 0.0085/0.0117 (`forensics/miss_hypotheses.md:189`), and c3's runtime IM calibration corrects 4/4 charges (`c3:163-193`). Every IM-column edit so far failed through coupling (affine −34.7%, `M/odia-im-affine-correction-wins.md:54-111`; oracle-centred half-library f7_w1cent region −1.84%, 317/317 ordinals negative, scorer saturated, `L/analysis77/SESSION_NOTES.md:2428-2434`).
2. **Between depth 37,334 and the operating point — ≈10.4 pts (~3,900).** DIA-NN precursors at ODIA ranks ~20k-37k interleaved with entrapments/decoys. Not a threshold problem: the band is a rate, a pure threshold move buys zero matched-budget recovery (`M/odia-matched-score-is-the-test.md:12`); admitting the correctly-ranked-but-rejected bucket costs 51.3% FDP (`D/51:132-144`).
3. **Admission — ≈2.3 pts (~860).** 646 absent (281 by `-min_library_fragments 3`), 131 no_candidate (`D/77:22-24`). Gate C is inert on dn_pred_cam: tau=0 because the first-20k-decoy null is 99.9% zeros (`L/analysis77/x1_fragvec_r1.log:50,155`; `c3:61,197` "gate C 0/0") — there is no gate to open.
4. **Localisation — ≈0.3 pts.** Forcing candidates at DIA-NN's apexes with admission bypassed: +98 ids at depth 33,330 (`D/79:22-35`; `M/odia-oracle-localization-is-minor.md:23-25`).
5. **Scorer regime instability — a reliability defect, not a steady loss.** Fold seed 42→43 moves q≤0.01 16,965→4,105 with no byte changed (`M/odia-scorer-regime-flips.md:102-108`); one fold's rounds 2-5 oscillate (leaves 214/127) at a 1:155 positive:negative prior per fold (~574:1 over the whole run — 10,897,786 decoy groups vs 18,995, `c3:227-228`; codex Q6 — the "155" is the per-fold figure, not an error, but quote it as per-fold) — negatives are already one best row per decoy (`top_decoys_only=true`, `include/odia/scoring/lda.h:137,1261-1289`), ~2.5M per fold vs ~16.5k positives (`c3:80-82`). The leaf cap restores seed 43 (17,727) but its human_v2 replicate FAILED through Gate C's tau re-set (`M/odia-leaf-cap-b2.md:77-84,93-102`).
6. **FDR calibration inflates, does not hide.** Own q≤0.01 sits at 1.84-2.31% true FDP (`M/odia-vs-diann-measured.md:24`), over-count 1.56-2.24x and depth-dependent (`M/odia-decoy-null-borrowed-evidence.md:90-91`); recalibration changes nothing in N(e) (ibid:22-24). An honest 1% would read LOWER than 18,995.

Where the record is silent: (a) class C (buried) traces have never been compared trace-vs-trace against DIA-NN's `--xic` export at production coordinates — the "extraction exonerated" parity (0.270 vs 0.311, `D/44:6-12`) was measured on a both-found sample (contrarian lens #4); (b) no at-depth-37,334 oracle exists for the 42 columns or for FRAGVEC (wf_v50 reports matched-budget deltas only, `L/analysis77/pick/wf_v50_oof_b2_cap1_report.txt:445`); (c) whether one mechanism or two sits behind the multi-candidate class (A2 lane-join unrun, `FORENSICS:263-274`); (d) the per-fragment apex-offset decomposition on the forensics' A/C classes at DIA-NN's apex (`M/odia-why-we-miss.md:52-55`, listed as UNRUN) — `D/45:512-552`'s MAD 20.8 vs 2.8 s was measured on v7's own picked apexes (4,000 v7-rejected precursors, pre-forensics), not on the C class, so "dispersed, not anchoring" is INFERRED for the buried set.

### 2.2 Memory — 642.2 GiB, ranked by bytes
## 6. Plan (ordered, ≤8)

Conventions on every item: `ibmi-nodes.sh` load gate first; `kinit` before any dax work (Kerberos expiry fakes "run died"); `setsid nohup`; binary sha recorded; `cmp` outputs before reading any timing; two baselines (same-night control + retained best w1_ctl70 / b2_cap1 lineage) and the engine's own q≤0.01 on every ID read; three numbers (region median N(e) over FDP 2-10%, registered cell ± 1-entrapment envelope, common-depth contrast) on every ID claim; output-identity gate = md5 of the sorted score table (layered as in F02 where floating-point reduction order is involved) on every performance change.

1. **F02 arm + riding instruments (zero code).** Goal: largest memory lever, chunk-invariance test, remainder attribution, system-time attribution. Change: c3 command verbatim (`run_c3.sh:110-115`) with `-live_memory_gb 150`, same binary a568ddf0, seed 43, solo-slot on dax; lower the launcher's MemAvailable wait (`run_c3.sh:107-108`, sized at 800 GB for 657-GB arms) to ~450 GB in the arm's copy of the script, or the 150-GiB arm will not start while two 657-GB arms are resident; attach the 10 s `/proc/<pid>/status`+`smaps_rollup` sampler stamping log growth, per-thread `/proc/<pid>/task/*/stat` minflt/stime bins, `/proc/vmstat` NUMA/compaction deltas, `top -H` once during pass-1 sink (expect 1 task ~100%, 64 sleeping); read dax `numa_balancing`/THP/sockets. Node/cost: dax, ~9.5-10 h, ≤420 GB. Pre-registered: PASS = Peak ≤420,000 MiB AND raw `cmp` of the UNMODIFIED scored TSV clean against c3_cent_gate_s43.tsv (codex; §5 row 19) AND ~7 chunks, decode ≤1,300 s, wall ≤1.10x. If raw `cmp` fails, run the diagnostic ladder before any verdict — identical sorted (Id,Decoy) q≤0.01 set and entrapment count; `cmp` with DScore rounded to 5 s.f.; max rel |ΔDScore| — to classify the difference as order-only, last-digit reduction (`lda.h:842-846` → `gbt.h:507-528`, not a defect) or content (a defect); report Peak − 1.065×budget vs 221,445 MiB. FAIL branches: any q/PEP/sub-score/ID-set change → chunk-edge defect, stop budget work, run the 48-vs-64 thread arm at 400 to clear the cohort_plain confound; Peak >500k → remainder fixed → F06 only. Unblocks: 5-6 arms per node for everything below; F03/F12 sizing; F06 classification (a)/(b)/(c) from the HWM phase; F09 step 2 arm choice.
2. **F01 arms a and a' (zero code, behaviour-changing).** Goal: halve pass 1's wall and memory, test whether ±400 s is load-bearing. Change: c3 verbatim with `-rt_window_pass1 110` (a) and 200 (a'), same binary/seed, same night as a c3 re-run control if the node allows three. Node/cost: dax, 2 × 5.5-9.5 h at ≤300 GB each; if only one slot is free, run 200 first (codex, §5 row 20; its cost gate: wall ≤7.0 h). Pre-registered: external RT coverage of DIA-NN's 37,334 CAM apexes (canonical UniMod:4 join) at ±60/±110 s within 0.5 pt of c3; "pass 2 extraction window 110 s" line unchanged; per-charge 1/K0 anchors ≥120 and IM gate PASSED; centring ≥50% of 1,211,024; live ≤70 GiB, 1 chunk, pass-1 sink ≤5,000 s, Peak ≤300,000 MiB, wall ≤6.5 h; pass-2 region median N(e) ≥ −1% vs c3 and w1_ctl70, ≥90% ordinals ≥ −2%, q≤0.01 ≥0.9x, histogram not regime-flipped. FAIL: coverage −>0.5 pt or N(e) < −1% → width load-bearing; then a' decides whether 200 is the floor. Unblocks: a 2x wall/memory default candidate (a default is still an experiment: replicate on human_v2 + sentinel before shipping); halves every later arm's cost.
3. **F11 Stage 0 (offline, on dax).** Goal: decide whether the R1 fragment block carries at-depth ranking information. Change: `kinit`; verify `/scratch/kohlbach/odia2x2/wf_v63/x1_fragvec_r1_{A,Ris,Ri,Alog,Riv,Aperm}_N3_seed{1,2}.npz` and `x1_fragvec_r1.fragvec.tsv` (sha b89d2c15078b4bf1) exist; compute recall at depth 37,334 vs the 37,334 CAM reference for A, Ris, Ri; refit (wf_v50 pipeline, ~1.5-3 h at 120 columns) only if gone. Node/cost: dax, 4 threads, minutes-hours, no engine. Pre-registered: A > 64.1% and A reproduces wf_v50 acc/conc at N(62) to the digit, else the join is broken (stop); MDE = two-partition spread; C ≥ A+3 → fund Stage 1 (F11 engine arms F6/F72 vs b2_cap1 and w1_ctl70 + 6-null-append control + x2_fragvec_r1_b2cap1); A+2 < C < A+3 → Stage 1 conditional; C ≤ A+2 → do not wire; next offline step is rung (ii) coherence columns on the own_ctl chrom export. Unblocks: the only scorer-side ID lever with measured ranking-only gain, or its closure.
4. **Instrumented binary (F06 + F05 timer + F03 asserts), fixture-gated.** Goal: replace every subtraction with a printed term before writing memory code. Change (new build under `env.sh`, never a bare `cmake --build`): getrusage/VmRSS/VmHWM/task-count lines at `ChromatogramExtractor.cpp:942,1105,1251`, `Ms1Traces.cpp:129`, `OpenDIAlyzer.cpp:4734`; relabel `:4743`; `reservedPoints()`/`peakPoints()` per chunk and per pass; MS1 decode-vs-match split; `ppm_den[i]==base[i]` assert and `im_den!=base` count in `emit()` behind a flag. Node/cost: ibminode05 fixture (bench_d2t240n_ih1 config, 23 min, 40 GB) ×2 (old/new binary); then ride on the next full arm. Pre-registered: stage CPU-s within 5% of footer; sink and MS1 CPU/wall in [0.95,1.10]; TSV byte-identical to bd73e65; 0 ppm-plane violations (any violation falsifies F03's identity — stop); im violations 0 or >0 recorded. On the full arm: classify (a)/(b)/(c) as in §4 F06. Unblocks: F03 (plane count), F05 (decode floor), the trim-vs-pool decision, honest CPU shares.
5. **F05 frame-parallel MS1 build.** Goal: −5,000 s wall, byte-identical. Change: parallel over contiguous ≥16-frame runs within each driver-decoded 64-frame block, per-unit `resid` capped at RESID_CAP and concatenated in frame order; bucket index as a separate second commit. Node/cost: fixture at `-threads 1/8/64` (ibminode05, 3 × 23 min); one dax full arm. Pre-registered: fixture checksum + TSV md5 identical across thread counts (timing on the fixture is NOT acceptance); full run "8556564 with signal" and median residual identical, sorted TSV md5-identical to same-binary control, MS1 line ≤600 s. FAIL: identity breach → race, revert; <4x → bucket index next, measured alone. Unblocks: the serial pre-pass on every arm; the thread-invariance harness F04 needs (and that step 6 reuses); ~3.2x total with step 7; lands before F03 because it needs no plane assert (Kimi, §5 row 17).
6. **F03 plane aliasing at pinned cap (+ F12 row drop as a rider, after a `-no_ms1` additivity arm).** Goal: −170 GiB live planes at identical output; −6.4 GiB never-read MS1 rows. Change: alias `ppm_den` (and `im_den` iff step 4 counted 0; step 4 is this step's precondition, not step 5's) to base, delete the aliased `+=`, single take/give, Max mode unchanged, pass-1 mobility probe for no-axis runs; `-max_live_precursors 3324888` pinned on the full arm; drop the 1,275,674 isolation-uncovered MS1 rows with a library→row map. Node/cost: fixture `cmp` (ibminode05, 23 min); one dax full arm 9.5 h at ~450 GB; one `-no_ms1` arm (timing/memory only, not ID-comparable) if not already read. Pre-registered: fixture TSV `cmp` clean, "(3 planes)" printed; full arm "cap 3324888, 3 chunks, 59902 spectra", TSV `cmp` clean vs same-binary control, Peak in 420-450k MB (not merely ≤490k), MS1 line −6,535 MiB; Astral rep arm (30 min) "(3 planes)"/"(1 planes)" and `cmp` clean. FAIL: any cell mismatch or TSV diff → revert. Unblocks: 1.67x more precursors per chunk; with step 1, sub-300 GB arms; with step 2, sub-150 GB.
7. **F04 parallel sink by ordered speculation.** Goal: −~17,000 s wall, byte-identical. Change: FIRST fix the `rejects_` registry (per-task `PickerRejects` reduced on commit, or a destructor that unregisters; `PeakGroupScorer.cpp:1310-1332`; codex, verified) — without it the worker thread_locals dangle in `totalRejects()` after the pool dies; then per §4 F04 (full commit set, chunk-end flush barrier, per-task scratch, `-threads 0` fix, joint OMP budget). Node/cost: fixture md5 at 1/8/64 with new binary and bd73e65 (ibminode05); one dn_pred_cam full arm vs same-binary serial control (dax, 9.5 h → target ≤5 h); one human_v2 pair at `-threads 16/64` with a ≥v1.17 binary printing tau. Pre-registered: TSV md5 and row count identical everywhere; fixture pass-1 sink ≤25 s; full pass-1 sink ≤2,000 s, pass-2 ≤600 s, wall ≤5 h; same tau line at 16 vs 64 on human_v2. FAIL: any byte/row difference → bug, revert; sink >4,000 s → measure (per-thread stat), do not pre-attribute. Unblocks: 2.1x alone, ~5x with steps 2+5 → 3-5x more seeds/baselines per week for the ID programme, which is rate-limited by arms per node-day.
8. **Calibration-only pass 1 (codex Q4; after step 2 has read).** Goal: remove the 60%-of-wall, budget-saturating pass instead of narrowing it. Change: pass 1 extracts and scores a deterministic RT/charge/m/z-stratified panel of targets PLUS a matched decoy panel (target-only is invalid: positive selection is `label == 1 && qvalue <= train_fdr`, `lda.h:1163-1166`), harvests RT/1/K0/centring anchors from it, and skips the full `Sink`/`finish()` (`OpenDIAlyzer.cpp:4718-4735`); pass 2 unchanged. Node/cost: days of code; fixture; one dax arm per panel size (e.g. 10% / 25% of the library). Pre-registered: the same EXTERNAL gates as step 2 — coverage of DIA-NN's 37,334 CAM apexes at ±60/±110 s within 0.5 pt of c3, per-charge 1/K0 anchors ≥120 and IM gate PASSED, centring survives (or is shown not to matter under the frozen model), "pass 2 extraction window" line unchanged; cost: pass-1 wall ≤2,000 s, live ≤50 GiB; benefit: pass-2 region median N(e) ≥ −1% vs c3 and w1_ctl70, q≤0.01 ≥0.9x. FAIL: any anchor floor missed at 25% → the panel cannot replace the census and step 2's width is the lever. Unblocks: a sub-200 GB, sub-5 h ODIA without touching the scorer; combined with 5 and 7 the 3-5x target.

Slot fillers — **two scorer-reliability arms (residues of F14 and F07), one variable each.** Goal: (i) settle whether the seed-43 collapse is the prior step; (ii) read the built-but-never-read `-gate_calibration_rt_min` on dn_pred_cam. Change: (i) r1_s43 flags + `-gbt_intercept_zero` (flag exists, `gbt.h:88`; `OpenDIAlyzer.cpp:1228-1247`), read vs r1_s43 (4,105) and i1_int0 (seed 42); (ii) c3/b2 flags + `-gate_calibration_rt_min 450-500` (just past DIA-NN's earliest ID at 408 s, `M/odia-fixture-cannot-arm-gate-c.md:23`; `L/analysis77/SESSION_NOTES.md:3375`; the first 20k decoys currently sample RT 264.9-320.1 s, `x1_fragvec_r1.log:50`). Node/cost: dax, 2 × 9.5 h at 657 GB (or cheaper after steps 1/2). Pre-registered: (i) PASS = q≤0.01 ≥0.9× the cross-seed median AND region vs b2_cap1 ≥ −1% AND max |leaf| in any fold's first 5 trees <20 (from the .model file); (ii) "gate C null armed" prints tau>0, zeros <50%, calibration RT inside the gradient; PASS = region median N(e) ≥ −1% vs b2_cap1 and w1_ctl70 with no ordinal < −2% AND pass-1 sink ≤0.5 × 14,637 s; FAIL = region < −2% (human_v2 analogue: admission slices carry real ids). Both with a fixed-model (`-classifier_model_in`) re-score beside the native read. Unblocks: whether scorer reads need ≥3 seeds forever; whether arming the gate buys 2/3 of the sink for free on dn_pred_cam.

Sequencing logic: 1-3 are zero-code and can run the same night (two node arms + one offline job; their gates tolerate the 14% cross-talk, but wall reads need a solo node); 4 is the precondition for 6-7 (5 needs only its own identity gate; step 4's MS1 decode/match split only says whether 4x is reachable); 8 follows step 2's read (it is the code form of the same question); the slot fillers take any idle slot. The byte-weighted chunker (§5 row 23) rides on step 6. Measurements before code wherever §2 said "silent": the HWM phase, the match efficiency, the sink scaling, the MS1 decode floor, the system-time owner all come out of steps 1 and 4 before a line of memory or scheduling code is written.

---

## 7. What this report could not establish
```

## 7. PASTED SOURCE at bd73e65 (cat -n numbers). You have NO working shell or file access if you are Codex: review from this text.
### include/odia/scoring/lda.h (whole file: semi-supervised loop, folds, positive selection, nested OMP sizing)
```
     1	// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
     2	// SPDX-License-Identifier: BSD-3-Clause
     3	//
     4	// Adopted from okohlbacher/OpenDIAlyzer, branch odia-engine-and-openms-boundary,
     5	// tag odia-v0.3.0 (0ed5fb2), file src/odia_lda.h. Same author, same BSD-3 licence.
     6	//
     7	// Changed here: namespace `odia` -> `ODIA::Scoring` to match this project, and
     8	// the sibling includes rewritten to this tree's layout. The algorithms are
     9	// untouched.
    10	//
    11	// Provenance that travels with the code, per that repository's CLEAN-ROOM.md:
    12	// this is an independent reimplementation of a published confidence model
    13	// (Reiter et al. mProphet, Nat Methods 2011; Teleman et al. pyprophet,
    14	// Bioinformatics 2015). No DIA-NN code is present in it.
    15	// odia_lda.h — in-process semi-supervised LDA scorer for OpenDIAlyzer.
    16	//
    17	// A pyprophet/mProphet-style confidence model, reimplemented in-process so
    18	// OpenDIAlyzer needs no external pyprophet. Pure C++ (no OpenMS, no Eigen), so it
    19	// compiles and self-tests standalone. Used in two places:
    20	//   (1) recalibration: pick RT anchors as targets at q < threshold (a real
    21	//       semi-supervised d-score, not raw VAR_XCORR_SHAPE/VAR_LIBRARY_CORR cutoffs);
    22	//   (2) final FDR: q-values on the last pass.
    23	//
    24	// Start with LDA (deterministic, dependency-light). Gradient-boosted trees are the
    25	// backlog upgrade (see docs/OpenDIAlyzer-scoring-fdr-backlog.md).
    26	//
    27	// ALGORITHM (target for implementation):
    28	//   Inputs: features (N rows x M sub-scores), labels (1=target, 0=decoy), group
    29	//   (precursor id; a precursor has several candidate peak groups = rows). FDR is a
    30	//   per-precursor quantity, so it is computed on the best-scoring row per group.
    31	//   1. z-standardize each feature column (mean/sd over all rows).
    32	//   2. k-fold cross-validation BY GROUP (all rows of a precursor fall in one fold),
    33	//      so a row is never scored by a model trained on its own precursor.
    34	//   3. For each fold, train on the other folds by semi-supervised iteration:
    35	//        - initialise the discriminant with the single most target/decoy-separating
    36	//          feature (max |two-sample t|);
    37	//        - pick a confident-target training set: best row per target group whose
    38	//          current score beats the decoy-derived cutoff at a lenient train FDR
    39	//          (e.g. 0.05); all decoy rows are negatives;
    40	//        - Fisher LDA: w = Sw^{-1} (mu_pos - mu_neg), Sw = pooled within-class
    41	//          covariance (solve SPD system by Cholesky; ridge-regularise the diagonal
    42	//          if needed); score = w . x; repeat n_iter times.
    43	//      Apply the fold's final w to its held-out rows -> cross-validated d-score.
    44	//   4. q-values: take the best d-score per precursor (targets and decoys), sort
    45	//      descending; at each cut FDR = #decoys_above / #targets_above (target-decoy),
    46	//      monotonise to q-values; broadcast each precursor's q to its rows.
    47	//   Deterministic given `seed` (only the fold assignment is randomised).
    48	//
    49	#ifndef ODIA_LDA_H
    50	#define ODIA_LDA_H
    51	
    52	#include <odia/scoring/gbt.h>
    53	#include <odia/scoring/anchor_training.h>
    54	
    55	#include <algorithm>
    56	#include <cmath>
    57	#include <cstddef>
    58	#include <cstdio>
    59	#include <fstream>
    60	#include <functional>
    61	#include <limits>
    62	#include <locale>
    63	#include <random>
    64	#include <string>
    65	#include <unordered_map>
    66	#include <utility>
    67	#include <vector>
    68	
    69	namespace ODIA::Scoring
    70	{
    71	
    72	struct ScoredGroups
    73	{
    74	  std::vector<double> dscore;   ///< per input row (peak group)
    75	  std::vector<double> qvalue;   ///< per input row (its precursor's q-value)
    76	  std::vector<double> pvalue;   ///< per input row: tail probability under the decoy null
    77	  std::vector<double> pep;      ///< per input row: local FDR (posterior error probability)
    78	  /// Diagnostics: how many semi-supervised iterations actually FITTED a discriminant, and how many
    79	  /// bailed out for want of confident positives. If trained==0 the returned scores come from the
    80	  /// single-feature initialisation, not from a learned model -- which looks like a working LDA and
    81	  /// is not one.
    82	  int n_iterations_trained = 0;
    83	  int n_iterations_skipped = 0;
    84	};
    85	
    86	struct LDAParams
    87	{
    88	  int n_folds = 3;          ///< cross-validation folds (by group)
    89	  int n_iter = 3;           ///< semi-supervised iterations
    90	  double train_fdr_initial = 0.15;  ///< FDR for the FIRST training-set selection. Deliberately more
    91	                            ///< lenient than the rest: the first pass is seeded by a single feature,
    92	                            ///< so a strict cut can select too few positives to fit anything and the
    93	                            ///< model never bootstraps. pyprophet uses 0.15 here for the same reason.
    94	  double train_fdr = 0.05;  ///< FDR for subsequent iterations, once a real discriminant exists.
    95	  double ridge = 1e-6;      ///< diagonal regularisation for the within-class covariance solve
    96	
    97	  /// Feature columns whose weight may never be positive, because the caller
    98	  /// knows the feature is lower-is-better. Empty means unconstrained.
    99	  ///
   100	  /// DIA-NN's check_weights does this for its RT-deviation and mass-accuracy
   101	  /// features. Applies to the LDA weight vector only; a tree ensemble has no
   102	  /// weights to clip.
   103	  std::vector<std::size_t> nonpositive_features;
   104	  unsigned seed = 42;       ///< RNG seed (fold assignment only) — determinism
   105	  bool use_pi0 = false;     ///< Storey pi0 correction. false = HONEST/conservative (true 1% FDR,
   106	                            ///< fewer IDs); true = PYPROPHET parity (more IDs, but a nominal 1% is
   107	                            ///< ~2% actual). NOT DIA-NN: 1.7.x uses a plain target-decoy count
   108	                            ///< ratio with NO pi0 term (verified by literature survey 2026-08-15),
   109	                            ///< so citing DIA-NN here was wrong and is the kind of stale comment
   110	                            ///< that leaks into documents.
   111	
   112	  /// Draw each decoy's best score from as many candidates as a TARGET has,
   113	  /// rather than from all of its own.
   114	  ///
   115	  /// Target-decoy competition assumes the two classes are exchangeable under
   116	  /// the null. They are not here, and the gap is large: measured on Astral,
   117	  /// targets carry 10.98 candidates per precursor (median 6) and decoys 20.99
   118	  /// (median 24), because the picker's margin rule keeps candidates within
   119	  /// `max_corr_diff` of a precursor's OWN best -- so a precursor with no real
   120	  /// peak admits nearly every position and one with a strong peak admits few.
   121	  ///
   122	  /// A q-value then compares best-of-11 against best-of-21, over 9,453 decoy
   123	  /// precursors. The maximum of that many best-of-21 draws reached 5.997
   124	  /// against a target p99 of 6.517, which is what set the 1% threshold: ODIA
   125	  /// ranked 9,563 Astral precursors that are ALL in DIA-NN's truth set and its
   126	  /// own FDR refused 3,834 of them.
   127	  ///
   128	  /// The cap per decoy group is quantile-matched to the target count
   129	  /// distribution, so the two classes end up drawing from the same
   130	  /// distribution of N. Targets are untouched.
   131	  ///
   132	  /// Which candidates are kept is decided in CANONICAL ORDER (precursor, apex
   133	  /// RT, apex intensity), never by score. Keeping a decoy's top-K by dscore
   134	  /// would subsample exactly the rows that make the null hard and would
   135	  /// manufacture identifications rather than measure them.
   136	  bool match_decoy_candidate_counts = false;
   137	  bool top_decoys_only = true;  ///< Train the negative class on each decoy precursor's BEST peak
   138	                            ///< group only, not on all of its candidate peak groups. See the note at
   139	                            ///< the negative-class construction: the positive class is top-peaks-only
   140	                            ///< by definition, so using all decoy rows makes the two classes differ
   141	                            ///< in peak RANK as well as in label, and the discriminant partly learns
   142	                            ///< "is this the best peak of its group" instead of "is this a target".
   143	                            ///< pyprophet trains on get_top_decoy_peaks() for exactly this reason.
   144	  /// Which learner fills the semi-supervised loop. LDA fits one hyperplane; GBT fits an additive
   145	  /// ensemble of depth-limited trees and can express interactions between sub-scores that no
   146	  /// hyperplane can. Measured on synthetic data with the SAME loop around it: on linearly separable
   147	  /// data LDA 0.981 vs GBT 0.986 AUC (no penalty for the extra capacity), on an interaction LDA
   148	  /// 0.514 -- chance -- vs GBT 0.926. Everything else (folds, training-set selection, fold
   149	  /// normalisation, q-values) is identical, so the two are directly comparable.
   150	  ///
   151	  /// NN adds a third: an ensemble of small tanh MLPs, the shape DIA-NN uses. It can express the
   152	  /// same interactions the GBT can, plus smooth ones a tree approximates in steps.
   153	  enum class Classifier { LDA, GBT, NN };
   154	  Classifier classifier = Classifier::LDA;
   155	  GBTParams gbt;            ///< only consulted when classifier == GBT
   156	  NNParams nn;              ///< only consulted when classifier == NN
   157	
   158	  // ---- the anti-circularity mechanisms of odia_anchor_training.h, each toggleable ALONE ---------
   159	  // They are separable on purpose. Bundled, an improvement and a regression cancel and the net
   160	  // result is "no effect"; separately, each one's contribution is attributable.
   161	
   162	  /// Mechanism 1. Feature mask for the SEED fit only (empty = no exclusion). The seed selects the
   163	  /// anchors; if it selects them using the very features they will be used to calibrate, the
   164	  /// correction is biased toward the uncorrected state.
   165	  std::vector<char> seed_mask;
   166	  /// Mechanism 3. Share of precursor GROUPS each ensemble member trains on (1.0 = off, which is
   167	  /// DIA-NN's arrangement: members differ only in weight init, so they share one seed bias).
   168	  double bag_fraction = 1.0;
   169	  /// Mechanism 5. Stop when the positive SET stops changing, or when it SHRINKS. Counting
   170	  /// identifications cannot see a collapse -- they rise throughout one.
   171	  bool stop_on_composition = false;
   172	  double stop_jaccard = 0.98;
   173	  /// Collapse/convergence policy forwarded to anchorIterationShouldContinue. The library defaults
   174	  /// are the LEGACY single-shot rule (any shrink = collapse); the CLI registers the repaired
   175	  /// high-water + patience policy (floor 0.01, patience 2) as ITS defaults, so only armed CLI
   176	  /// runs get the new rule and non-CLI callers cannot drift. See AnchorTrainingParams.
   177	  double stop_shrink_floor = 0.0;
   178	  int stop_patience = 1;
   179	  /// Stderr-only churn diagnostic: one line per (fold, iteration) with the positive-set size and
   180	  /// its Jaccard overlap with the previous iteration's set. Never control flow, and it shares no
   181	  /// state with mechanism 5 -- an iteration-sweep arm carries it precisely because it cannot
   182	  /// change the output bytes. Exists because the offline replay (analysis77 R0) measured 11-36%
   183	  /// positive-set churn per iteration with ID counts still rising; whether the REAL loop
   184	  /// converges is what this makes visible.
   185	  bool iteration_log = false;
   186	  /// DIAGNOSTIC ONLY, and FDR-INVALID when true: every group trains the model that scores it.
   187	  ///
   188	  /// It exists because `n_folds = 1` cannot express this -- the fold count is clamped to >= 2 a few
   189	  /// lines into the routine, so an ablation arm that set n_folds=1 silently ran 2-fold CV and
   190	  /// "demonstrated" nothing while claiming to demonstrate that cross-validation is load-bearing.
   191	  /// A flag that says what it does cannot be defeated by a clamp.
   192	  bool disable_cv = false;
   193	  bool normalize_folds = true;  ///< Rescale each fold's held-out scores to its own decoy null
   194	  bool fold_pool_rank = false;  ///< v1.13: pool the folds by WITHIN-FOLD RANK FRACTION over all groups (targets: best
   195	                                ///< over all rows; decoys: best over the drawn prefix, as assignQValues ranks them),
   196	                                ///< DScore = -log10(rank/n_f), instead of (raw - mu_f)/sigma_f. The decoy null of an
   197	                                ///< uncapped fold does not fix the scale of its top (w1_ctl70 fold 1: top compressed
   198	                                ///< 0.75x, bulk sd 1.04x), so the pooled head is partly sorted by fold. mu/sigma are
   199	                                ///< still computed and saved. Off = native, identical arithmetic.
   200	                            ///< (mean 0, sd 1) before pooling. Each fold has its OWN weight vector,
   201	                            ///< with its own arbitrary scale and offset, so the raw scores are not
   202	                            ///< comparable across folds; pooling them into one ranking without this
   203	                            ///< mixes incommensurable scales.
   204	};
   205	
   206	/// Semi-supervised LDA scoring with cross-validation and target-decoy q-values.
   207	/// features[i] has the same length for all i (M sub-scores). labels[i] in {0,1}.
   208	/// group[i] is the precursor id shared by that precursor's candidate peak groups.
   209	/// Returns a d-score and q-value per input row. Deterministic given params.seed.
   210	/// `model_out` / `model_in` (gbt engine only) save the trained fold ensemble, or apply a saved
   211	/// one without training; see the definition.
   212	ScoredGroups scoreSemiSupervisedLDA(const std::vector<std::vector<double>>& features,
   213	                                    const std::vector<int>& labels,
   214	                                    const std::vector<long long>& group,
   215	                                    const LDAParams& params = LDAParams(),
   216	                                    const std::string& model_out = "",
   217	                                    const std::string& model_in = "",
   218	                                    const std::vector<std::string>* feature_names = nullptr);
   219	
   220	namespace lda_detail
   221	{
   222	
   223	struct RankedGroup
   224	{
   225	  std::size_t group_index;
   226	  std::size_t best_row;
   227	  int label;
   228	  double score;
   229	  double qvalue;
   230	  double pvalue = 1.0;   ///< empirical p from the decoy null (TAIL probability at this score)
   231	  double pep = 1.0;      ///< posterior error probability = LOCAL false-discovery rate at this score
   232	};
   233	
   234	inline double dot(const std::vector<double>& a, const std::vector<double>& b)
   235	{
   236	  double result = 0.0;
   237	  for (std::size_t j = 0; j < a.size(); ++j) { result += a[j] * b[j]; }
   238	  return result;
   239	}
   240	
   241	// Assign target-decoy q-values to a list containing one best score per group.
   242	// Equal scores are treated as one threshold, avoiding order-dependent q-values.
   243	inline void assignQValues(std::vector<RankedGroup>& ranked, bool use_pi0)
   244	{
   245	  std::sort(ranked.begin(), ranked.end(), [](const RankedGroup& a, const RankedGroup& b) {
   246	    if (a.score != b.score) { return a.score > b.score; }
   247	    return a.group_index < b.group_index;
   248	  });
   249	
   250	  std::size_t Ntar = 0, Ndec = 0;
   251	  for (const auto& r : ranked) { if (r.label == 1) { ++Ntar; } else { ++Ndec; } }
   252	
   253	  // Storey pi0 = estimated fraction of TRUE-NULL targets. Walking high->low score, a
   254	  // target's empirical p-value from the decoy null is (decoys_seen_so_far / Ndec); null
   255	  // targets have ~uniform p, so the mass with p>lambda estimates pi0. Without pi0 the
   256	  // target-decoy estimator assumes pi0=1 (every target could be false) and is
   257	  // over-conservative vs pyprophet/DIA-NN (measured: 22,959 vs 37,539 IDs on the same osw,
   258	  // despite our discriminant separating MORE clean targets). pi0 recovers the honest count.
   259	  double pi0 = 1.0;
   260	  if (use_pi0 && Ntar > 0 && Ndec > 0)
   261	  {
   262	    const double lambda = 0.5;
   263	    const std::size_t dthr = static_cast<std::size_t>(lambda * static_cast<double>(Ndec));
   264	    std::size_t dec_seen = 0, tar_hi = 0;
   265	    for (const auto& r : ranked)
   266	    {
   267	      if (r.label == 0) { ++dec_seen; }
   268	      else if (dec_seen > dthr) { ++tar_hi; }
   269	    }
   270	    pi0 = static_cast<double>(tar_hi) / ((1.0 - lambda) * static_cast<double>(Ntar));
   271	    if (!(pi0 > 0.0)) { pi0 = 1.0 / static_cast<double>(Ntar); }   // never 0
   272	    if (pi0 > 1.0) { pi0 = 1.0; }                                   // never > 1
   273	  }
   274	
   275	  std::size_t targets = 0;
   276	  std::size_t decoys = 0;
   277	  for (std::size_t begin = 0; begin < ranked.size();)
   278	  {
   279	    std::size_t end = begin + 1;
   280	    while (end < ranked.size() && ranked[end].score == ranked[begin].score) { ++end; }
   281	    for (std::size_t i = begin; i < end; ++i)
   282	    {
   283	      if (ranked[i].label == 1) { ++targets; }
   284	      else                      { ++decoys; }
   285	    }
   286	    // FDR(t) = pi0 * (decoys_above/Ndec) / (targets_above/Ntar)  -- p-value based, ratio-corrected.
   287	    double fdr;
   288	    if (targets == 0) { fdr = std::numeric_limits<double>::infinity(); }
   289	    else
   290	    {
   291	      // +1 finite-sample correction on the decoy count (Käll): keeps the estimator
   292	      // honest at the tail where single decoys otherwise make it anti-conservative.
   293	      const double dr = (Ndec > 0) ? (static_cast<double>(decoys) + 1.0) / static_cast<double>(Ndec) : 0.0;
   294	      const double tr = static_cast<double>(targets) / static_cast<double>(Ntar);
   295	      fdr = std::min(1.0, pi0 * dr / tr);
   296	    }
   297	    for (std::size_t i = begin; i < end; ++i) { ranked[i].qvalue = fdr; }
   298	    begin = end;
   299	  }
   300	
   301	  double running_min = 1.0;
   302	  for (std::size_t end = ranked.size(); end > 0;)
   303	  {
   304	    std::size_t begin = end - 1;
   305	    while (begin > 0 && ranked[begin - 1].score == ranked[end - 1].score) { --begin; }
   306	    running_min = std::min(running_min, ranked[begin].qvalue);
   307	    for (std::size_t i = begin; i < end; ++i) { ranked[i].qvalue = running_min; }
   308	    end = begin;
   309	  }
   310	
   311	  // --- p-value and PEP -------------------------------------------------------------------
   312	  // These were previously never computed: the caller wrote the q-value into the PVALUE, QVALUE and
   313	  // PEP columns alike, so 100% of rows had all three identical and anything downstream reading PEP
   314	  // (IPF, protein-level inference) was silently consuming a q-value. They are different statistics:
   315	  //
   316	  //   p-value : TAIL probability under the null -- P(a null score >= this one), from the decoys.
   317	  //   q-value : minimum FDR of the SET selected at this threshold (already computed above).
   318	  //   PEP     : LOCAL FDR -- the probability that THIS peak group specifically is null. A group at
   319	  //             q = 0.01 sitting right at the threshold can easily have PEP ~ 0.3.
   320	  //
   321	  // p is the standard conservative empirical estimate (Käll's +1 on both counts).
   322	  {
   323	    std::size_t dec_at_or_above = 0;
   324	    for (std::size_t begin = 0; begin < ranked.size();)
   325	    {
   326	      std::size_t end = begin + 1;
   327	      while (end < ranked.size() && ranked[end].score == ranked[begin].score) { ++end; }
   328	      for (std::size_t i = begin; i < end; ++i) { if (ranked[i].label == 0) { ++dec_at_or_above; } }
   329	      const double p = (Ndec > 0)
   330	                         ? (static_cast<double>(dec_at_or_above) + 1.0) / (static_cast<double>(Ndec) + 1.0)
   331	                         : 1.0;
   332	      for (std::size_t i = begin; i < end; ++i) { ranked[i].pvalue = std::min(1.0, p); }
   333	      begin = end;
   334	    }
   335	  }
   336	
   337	  // PEP by the LOCAL analogue of the global estimator used above: in a score neighbourhood holding
   338	  // t targets and d decoys, the decoys estimate the null target density (scaled by Ntar/Ndec and
   339	  // pi0), so the expected null-target count is pi0*d*(Ntar/Ndec) and PEP ~ that divided by t.
   340	  // Window is a fixed fraction of the list so it adapts to size; a raw local ratio is very noisy,
   341	  // hence the monotonicity pass afterwards.
   342	  {
   343	    const std::size_t n = ranked.size();
   344	    const std::size_t half = std::max<std::size_t>(50, n / 200);   // ~0.5% of the list, >=100 wide
   345	    const double scale = (Ndec > 0) ? (static_cast<double>(Ntar) / static_cast<double>(Ndec)) : 0.0;
   346	    // SLIDING counts, not a recount per element. The window is n/200 wide, so recomputing it inside
   347	    // the loop is O(n^2/100): at n = 2M peak groups that is ~4e10 operations and turns a seconds-long
   348	    // step into a multi-minute one. Maintaining running counts makes the whole pass O(n) -- each
   349	    // element enters and leaves the window exactly once.
   350	    std::size_t t = 0, d = 0, lo = 0, hi = 0;
   351	    for (std::size_t i = 0; i < n; ++i)
   352	    {
   353	      const std::size_t new_lo = (i > half) ? i - half : 0;
   354	      const std::size_t new_hi = std::min(n, i + half + 1);
   355	      while (hi < new_hi) { if (ranked[hi].label == 1) { ++t; } else { ++d; } ++hi; }
   356	      while (lo < new_lo) { if (ranked[lo].label == 1) { --t; } else { --d; } ++lo; }
   357	      double pep = 1.0;
   358	      if (t > 0) { pep = pi0 * static_cast<double>(d) * scale / static_cast<double>(t); }
   359	      ranked[i].pep = std::min(1.0, std::max(0.0, pep));
   360	    }
   361	    // PEP must not increase with score. The list is sorted high->low by score, so PEP must be
   362	    // NON-DECREASING IN INDEX, and the sweep therefore runs from the high-score end forwards.
   363	    //
   364	    // It ran backwards. Sweeping from the low-score end with a running max gave every entry the
   365	    // maximum PEP of its SUFFIX -- i.e. of everything scoring below it -- which is monotone in
   366	    // the wrong direction and, on a real list, catastrophic: [0.01, 0.02, 0.5, 0.9] came out
   367	    // [0.9, 0.9, 0.9, 0.9], so the best-scoring precursor in the run was reported with the worst
   368	    // PEP in the run. Found by external review.
   369	    //
   370	    // Confined to the reported PEP column: nothing thresholds on it (it reaches only the output
   371	    // TSV via OpenDIAlyzer.cpp), so no q-value or identification count was affected.
   372	    //
   373	    // The old comment also called this "the isotonic projection". A cumulative maximum is the
   374	    // greatest monotone minorant, not a least-squares isotonic regression; the claim is dropped
   375	    // rather than repeated.
   376	    double running_max = 0.0;
   377	    for (std::size_t i = 0; i < n; ++i)
   378	    {
   379	      running_max = std::max(running_max, ranked[i].pep);
   380	      ranked[i].pep = running_max;
   381	    }
   382	  }
   383	}
   384	
   385	// Solve A x = b for symmetric positive-definite A. The input matrix is
   386	// row-major and is replaced by its lower-triangular Cholesky factor.
   387	inline bool choleskySolve(std::vector<double> a,
   388	                          const std::vector<double>& b,
   389	                          std::vector<double>& x)
   390	{
   391	  const std::size_t m = b.size();
   392	  for (std::size_t i = 0; i < m; ++i)
   393	  {
   394	    for (std::size_t j = 0; j <= i; ++j)
   395	    {
   396	      double value = a[i * m + j];
   397	      for (std::size_t k = 0; k < j; ++k) { value -= a[i * m + k] * a[j * m + k]; }
   398	      if (i == j)
   399	      {
   400	        if (!(value > 0.0) || !std::isfinite(value)) { return false; }
   401	        a[i * m + i] = std::sqrt(value);
   402	      }
   403	      else
   404	      {
   405	        a[i * m + j] = value / a[j * m + j];
   406	      }
   407	    }
   408	  }
   409	
   410	  std::vector<double> y(m, 0.0);
   411	  for (std::size_t i = 0; i < m; ++i)
   412	  {
   413	    double value = b[i];
   414	    for (std::size_t j = 0; j < i; ++j) { value -= a[i * m + j] * y[j]; }
   415	    y[i] = value / a[i * m + i];
   416	  }
   417	
   418	  x.assign(m, 0.0);
   419	  for (std::size_t ii = m; ii > 0; --ii)
   420	  {
   421	    const std::size_t i = ii - 1;
   422	    double value = y[i];
   423	    for (std::size_t j = i + 1; j < m; ++j) { value -= a[j * m + i] * x[j]; }
   424	    x[i] = value / a[i * m + i];
   425	    if (!std::isfinite(x[i])) { return false; }
   426	  }
   427	  return true;
   428	}
   429	
   430	} // namespace lda_detail
   431	
   432	/// @param model_out  gbt engine only: after training, save the fold models, the
   433	///                   standardisation and each fold's decoy normalisation here.
   434	/// @param model_in   gbt engine only: do not train -- load that file and score every
   435	///                   row with it. The two paths mirror `scorePercolator`'s.
   436	inline ScoredGroups scoreSemiSupervisedLDA(
   437	  const std::vector<std::vector<double>>& features,
   438	  const std::vector<int>& labels,
   439	  const std::vector<long long>& group,
   440	  const LDAParams& params,
   441	  const std::string& model_out,
   442	  const std::string& model_in,
   443	  const std::vector<std::string>* feature_names)
   444	{
   445	  const std::size_t n = features.size();
   446	  ScoredGroups result;
   447	  result.dscore.assign(n, 0.0);
   448	  result.qvalue.assign(n, 1.0);
   449	  result.pvalue.assign(n, 1.0);
   450	  result.pep.assign(n, 1.0);
   451	  if (n == 0 || labels.size() != n || group.size() != n) { return result; }
   452	
   453	  const std::size_t m = features.front().size();
   454	  if (m == 0) { return result; }
   455	  for (const auto& row : features)
   456	  {
   457	    if (row.size() != m) { return result; }
   458	  }
   459	
   460	  // Global z-standardisation is unsupervised. Non-finite cells are treated as
   461	  // missing and become zero (the column mean) after standardisation.
   462	  std::vector<double> mean(m, 0.0);
   463	  std::vector<double> count(m, 0.0);
   464	  for (const auto& row : features)
   465	  {
   466	    for (std::size_t j = 0; j < m; ++j)
   467	    {
   468	      if (std::isfinite(row[j]))
   469	      {
   470	        mean[j] += row[j];
   471	        count[j] += 1.0;
   472	      }
   473	    }
   474	  }
   475	  for (std::size_t j = 0; j < m; ++j)
   476	  {
   477	    if (count[j] > 0.0) { mean[j] /= count[j]; }
   478	  }
   479	
   480	  std::vector<double> sum_squared_deviation(m, 0.0);
   481	  std::vector<double> sd(m, 1.0);
   482	  std::vector<std::vector<double>> z(n, std::vector<double>(m, 0.0));
   483	  for (const auto& row : features)
   484	  {
   485	    for (std::size_t j = 0; j < m; ++j)
   486	    {
   487	      if (std::isfinite(row[j]))
   488	      {
   489	        const double d = row[j] - mean[j];
   490	        sum_squared_deviation[j] += d * d;
   491	      }
   492	    }
   493	  }
   494	  for (std::size_t j = 0; j < m; ++j)
   495	  {
   496	    sd[j] = count[j] > 1.0
   497	              ? std::sqrt(sum_squared_deviation[j] / (count[j] - 1.0))
   498	              : 0.0;
   499	    if (!(sd[j] > std::numeric_limits<double>::epsilon()) || !std::isfinite(sd[j]))
   500	    {
   501	      sd[j] = 1.0;
   502	    }
   503	  }
   504	  for (std::size_t i = 0; i < n; ++i)
   505	  {
   506	    for (std::size_t j = 0; j < m; ++j)
   507	    {
   508	      z[i][j] = std::isfinite(features[i][j]) ? (features[i][j] - mean[j]) / sd[j] : 0.0;
   509	    }
   510	  }
   511	
   512	  // Build a stable, first-occurrence ordering of precursor groups.
   513	  std::unordered_map<long long, std::size_t> group_lookup;
   514	  group_lookup.reserve(n);
   515	  std::vector<std::vector<std::size_t>> group_rows;
   516	  std::vector<int> group_label;
   517	  std::vector<long long> group_id;                 // the precursor id behind each group index
   518	  for (std::size_t i = 0; i < n; ++i)
   519	  {
   520	    auto inserted = group_lookup.emplace(group[i], group_rows.size());
   521	    if (inserted.second)
   522	    {
   523	      group_rows.emplace_back();
   524	      group_label.push_back(labels[i] == 1 ? 1 : 0);
   525	      group_id.push_back(group[i]);
   526	    }
   527	    const std::size_t g = inserted.first->second;
   528	    group_rows[g].push_back(i);
   529	  }
   530	
   531	  const std::size_t group_count = group_rows.size();
   532	  if (group_count < 2) { return result; } // leakage-free training is impossible
   533	
   534	  int folds = params.n_folds;
   535	  if (folds < 2) { folds = 2; }
   536	  if (folds > static_cast<int>(group_count)) { folds = static_cast<int>(group_count); }
   537	
   538	  std::vector<std::size_t> target_groups;
   539	  std::vector<std::size_t> decoy_groups;
   540	  for (std::size_t g = 0; g < group_count; ++g)
   541	  {
   542	    (group_label[g] == 1 ? target_groups : decoy_groups).push_back(g);
   543	  }
   544	  // Order these by the precursor's IDENTITY before shuffling. A group's index is its
   545	  // first-occurrence position in row order, and row order is whatever the parallel extraction
   546	  // happened to produce -- so a seeded shuffle of indices still put the same precursor in a
   547	  // different fold on every run, training a different model and reporting a different ID count.
   548	  // Measured spread on byte-identical input: 6487 / 6565 / 6433 (+/-1%), which is larger than
   549	  // most of the effects being A/B tested. Sorting by id first costs one sort and makes the fold
   550	  // assignment a function of the data alone; the shuffle still balances fold sizes exactly.
   551	  const auto by_id = [&](std::size_t a, std::size_t b) { return group_id[a] < group_id[b]; };
   552	  std::sort(target_groups.begin(), target_groups.end(), by_id);
   553	  std::sort(decoy_groups.begin(), decoy_groups.end(), by_id);
   554	  std::mt19937 rng(params.seed);
   555	  std::shuffle(target_groups.begin(), target_groups.end(), rng);
   556	  std::shuffle(decoy_groups.begin(), decoy_groups.end(), rng);
   557	  std::vector<int> group_fold(group_count, 0);
   558	  for (std::size_t i = 0; i < target_groups.size(); ++i)
   559	  {
   560	    group_fold[target_groups[i]] =
   561	      static_cast<int>(i % static_cast<std::size_t>(folds));
   562	  }
   563	  for (std::size_t i = 0; i < decoy_groups.size(); ++i)
   564	  {
   565	    group_fold[decoy_groups[i]] =
   566	      static_cast<int>(i % static_cast<std::size_t>(folds));
   567	  }
   568	
   569	  // Folds are independent BY CONSTRUCTION: fold f trains on the groups not assigned to f and writes
   570	  // result.dscore only for the rows of groups that ARE assigned to f. So no two iterations read or
   571	  // write the same dscore element, and none of them touch shared state except the two skip/train
   572	  // counters, which are reduced. Everything else (`z`, `group_rows`, `group_fold`) is read-only here.
   573	  // This loop was serial, and on the benchmark feature table the whole scoring step is 142 s of the
   574	  // 1,284 s serial tail that caps the run's speedup at 1.90x -- see
   575	  // docs/OpenDIAlyzer-parallel-efficiency.md. Cheap to fix, so fixed.
   576	  //
   577	  // Determinism is preserved: fold assignment comes from the seeded RNG above, each fold's model
   578	  // depends only on its own training set, and the q-values are computed after the loop from the
   579	  // pooled scores. The result does not depend on completion order.
   580	  // The GBT parallelises internally over rows, but fit() is called from INSIDE the fold loop below.
   581	  // A nested parallel region defaults to a team of one, so without raising the active-level limit
   582	  // the inner pragmas would be dead code. Splitting the available threads as
   583	  // folds x (threads/folds) uses the machine without oversubscribing it. The GBT's result does not
   584	  // depend on this number (odia_gbt_test T8), so it is purely a speed knob.
   585	  GBTParams gbt_params = params.gbt;
   586	  NNParams nn_params = params.nn;
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
   606	#endif
   607	
   608	  int n_trained = 0, n_skipped = 0;
   609	
   610	  // ---- FROZEN MODEL, gbt engine only ----
   611	  //
   612	  // Score every row with a fold ensemble trained elsewhere and saved by the
   613	  // `model_out` path below: THAT run's standardisation, its fold models, and
   614	  // each fold's decoy normalisation, applied here with no training at all.
   615	  // Measured 2026-09-05 (pick/wf_hist.txt, wf_fixedmodel.txt): the native
   616	  // retraining flips between a compact and a saturated score regime under
   617	  // small changes to the binary or the feature distribution, and every
   618	  // comparison across that flip is about the scorer, not the evidence; the
   619	  // same evidence under ONE model held fixed across arms is stable. A frozen
   620	  // model is how two runs are compared on evidence, and it is the documented
   621	  // remedy for a run that cannot bootstrap its own positives.
   622	  //
   623	  // The per-row score is the mean over folds of each fold's DECOY-normalised
   624	  // score, so it sits on the same "decoy sds above the null" scale the native
   625	  // path pools folds on. Rows the saving run never saw are, by construction,
   626	  // out of sample for every fold model.
   627	  const bool gbt_engine = (params.classifier == LDAParams::Classifier::GBT);
   628	  std::vector<GBT> fold_models;
   629	  std::vector<double> fold_norm_mu, fold_norm_sigma;
   630	  bool frozen_applied = false;
   631	  if (!model_in.empty() && !gbt_engine)
   632	  {
   633	    std::fprintf(stderr, "[lda] -classifier_model_in is honoured by the gbt and percolator "
   634	                         "engines only; this engine trains as usual\n");
   635	  }
   636	  if (!model_in.empty() && !model_out.empty() && gbt_engine)
   637	  {
   638	    std::fprintf(stderr, "[gbt] both -classifier_model_in and -classifier_model_out are set: "
   639	                         "the frozen model is APPLIED and nothing is saved\n");
   640	  }
   641	  if (!model_in.empty() && gbt_engine)
   642	  {
   643	    // A frozen model that cannot be applied is a HARD failure, never a silent
   644	    // fall-back to training: the retraining is the very thing this path exists
   645	    // to hold fixed, and a run that quietly trained would enter a comparison
   646	    // as if it had not. The failure surfaces as "fitted 0 iterations", exactly
   647	    // as a failed percolator engine does.
   648	    auto refuse = [&](const std::string& why) -> ScoredGroups {
   649	      std::fprintf(stderr, "[gbt] FROZEN MODEL %s: REFUSED -- %s; no scores produced\n",
   650	                   model_in.c_str(), why.c_str());
   651	      result.n_iterations_skipped = 1;
   652	      return result;
   653	    };
   654	    if (params.fold_pool_rank)
   655	    { return refuse("-fold_pool_rank is not supported with -classifier_model_in (v1.13): the pooled scale is a within-fold rank, not the saved mu/sigma"); }
   656	    std::ifstream is(model_in);
   657	    is.imbue(std::locale::classic());
   658	    if (!is) { return refuse("cannot open"); }
   659	    std::string magic;
   660	    int version = 0;
   661	    std::size_t m_saved = 0, k_saved = 0;
   662	    if (!(is >> magic >> version >> m_saved >> k_saved) || magic != "ODIA-GBT-FOLDS" || version != 1)
   663	    { return refuse("not an ODIA-GBT-FOLDS version-1 file"); }
   664	    if (m_saved != m)
   665	    {
   666	      return refuse("saved for " + std::to_string(m_saved) + " sub-scores, this run has " +
   667	                    std::to_string(m));
   668	    }
   669	    if (k_saved == 0 || k_saved > 64) { return refuse("implausible fold count"); }
   670	    // The sub-score NAMES travel with the model and must match exactly: a
   671	    // different sub-score set of the same width would load and rank on
   672	    // nonsense with nothing downstream the wiser.
   673	    std::size_t n_names = 0;
   674	    if (!(is >> n_names) || n_names != m) { return refuse("sub-score name list missing or wrong length"); }
   675	    for (std::size_t j = 0; j < m; ++j)
   676	    {
   677	      std::string nm;
   678	      if (!(is >> nm)) { return refuse("sub-score name list truncated"); }
   679	      if (feature_names != nullptr && j < feature_names->size() && (*feature_names)[j] != nm)
   680	      {
   681	        return refuse("sub-score " + std::to_string(j) + " is '" + nm + "' in the model but '" +
   682	                      (*feature_names)[j] + "' in this run");
   683	      }
   684	    }
   685	    std::vector<double> mean_saved(m, 0.0), sd_saved(m, 1.0);
   686	    std::vector<GBT> models(k_saved);
   687	    std::vector<double> mus(k_saved, 0.0), sigmas(k_saved, 0.0);
   688	    std::vector<char> present(k_saved, 0);
   689	    for (std::size_t j = 0; j < m; ++j)
   690	    { if (!(is >> mean_saved[j]) || !std::isfinite(mean_saved[j])) { return refuse("standardisation means truncated or non-finite"); } }
   691	    for (std::size_t j = 0; j < m; ++j)
   692	    {
   693	      if (!(is >> sd_saved[j]) || !std::isfinite(sd_saved[j]) || !(sd_saved[j] > 0.0))
   694	      { return refuse("standardisation sds truncated, non-finite or non-positive"); }
   695	    }
   696	    for (std::size_t f = 0; f < k_saved; ++f)
   697	    {
   698	      int has = 0;
   699	      if (!(is >> has >> mus[f] >> sigmas[f])) { return refuse("fold header truncated"); }
   700	      if (has)
   701	      {
   702	        if (!models[f].load(is)) { return refuse("fold " + std::to_string(f) + " model malformed"); }
   703	        if (models[f].featureCount() != m) { return refuse("fold model feature count mismatch"); }
   704	        present[f] = 1;
   705	      }
   706	    }
   707	    std::size_t k_have = 0;
   708	    for (const char p : present) { k_have += (p != 0); }
   709	    if (k_have == 0) { return refuse("no trained fold in the file"); }
   710	    // Optional (files written before v1.7 lack it): the saving run's precursor-
   711	    // to-fold assignment. A precursor listed here was TRAINING data for every
   712	    // fold but its own, and is scored below by that one excluded fold only.
   713	    std::unordered_map<long long, int> fold_of;
   714	    {
   715	      std::string tag;
   716	      std::size_t n_map = 0;
   717	      if (is >> tag)
   718	      {
   719	        if (tag != "FOLDMAP" || !(is >> n_map) || n_map > 50000000)
   720	        { return refuse("trailing content is not a valid FOLDMAP"); }
   721	        fold_of.reserve(n_map * 2);
   722	        for (std::size_t g = 0; g < n_map; ++g)
   723	        {
   724	          long long pid = 0;
   725	          int fd = -1;
   726	          if (!(is >> pid >> fd)) { return refuse("FOLDMAP truncated"); }
   727	          if (fd < 0 || static_cast<std::size_t>(fd) >= k_saved) { return refuse("FOLDMAP fold index out of range"); }
   728	          fold_of.emplace(pid, fd);
   729	        }
   730	      }
   731	    }
   732	    std::size_t scored_exact = 0;
   733	    {
   734	#ifdef _OPENMP
   735	#pragma omp parallel for schedule(static)
   736	#endif
   737	      for (std::size_t i = 0; i < n; ++i)
   738	      {
   739	        std::vector<double> zf(m, 0.0);
   740	        for (std::size_t j = 0; j < m; ++j)
   741	        {
   742	          zf[j] = std::isfinite(features[i][j]) ? (features[i][j] - mean_saved[j]) / sd_saved[j] : 0.0;
   743	        }
   744	        auto fold_score = [&](std::size_t f) -> double {
   745	          double v = models[f].score(zf);
   746	          if (sigmas[f] > std::numeric_limits<double>::epsilon() && std::isfinite(sigmas[f]))
   747	          { v = (v - mus[f]) / sigmas[f]; }
   748	          return v;
   749	        };
   750	        const auto it = fold_of.find(group[i]);
   751	        if (it != fold_of.end() && present[static_cast<std::size_t>(it->second)])
   752	        {
   753	          // Training data of every other fold: its own excluded fold only.
   754	          result.dscore[i] = fold_score(static_cast<std::size_t>(it->second));
   755	#ifdef _OPENMP
   756	#pragma omp atomic
   757	#endif
   758	          ++scored_exact;
   759	        }
   760	        else
   761	        {
   762	          double s = 0.0;
   763	          for (std::size_t f = 0; f < k_saved; ++f)
   764	          {
   765	            if (!present[f]) { continue; }
   766	            s += fold_score(f);
   767	          }
   768	          result.dscore[i] = s / static_cast<double>(k_have);
   769	        }
   770	      }
   771	      frozen_applied = true;
   772	      // Reported as trained iterations, as scorePercolator reports its frozen
   773	      // path: downstream derives fdr_valid from `> 0`, and a frozen model IS a
   774	      // valid discriminant. It is no longer "this run bootstrapped its own
   775	      // positives", which is what the name says -- the only consumer today is
   776	      // the `> 0` test, so the conflation is documented here rather than fixed.
   777	      n_trained = static_cast<int>(k_have);
   778	      std::fprintf(stderr, "[gbt] FROZEN MODEL %s: %zu fold models over %zu features applied to "
   779	                           "%zu rows (%zu rows of precursors the model trained on scored by their "
   780	                           "own excluded fold, the rest by the fold mean), no training. NOTE: "
   781	                           "q-values under a frozen model are an FDR "
   782	                           "estimate only if the saving run shared no targets with this one; when "
   783	                           "it did (same sample, same library) they are anti-conservative -- use "
   784	                           "them to COMPARE runs, not to report an FDR\n",
   785	                   model_in.c_str(), k_have, m, n, scored_exact);
   786	    }
   787	  }
   788	  // The captures below are sized only when a save was requested, so a run with
   789	  // neither path set allocates and copies nothing -- flag-off touches one bool.
   790	  if (!frozen_applied && gbt_engine && !model_out.empty())
   791	  {
   792	    fold_models.resize(static_cast<std::size_t>(folds));
   793	    fold_norm_mu.assign(static_cast<std::size_t>(folds), 0.0);
   794	    fold_norm_sigma.assign(static_cast<std::size_t>(folds), 0.0);
   795	  }
   796	
   797	  // How many candidates each group may draw its best from. Targets always use
   798	  // all of theirs; decoys are quantile-matched to the target distribution when
   799	  // asked. See `match_decoy_candidate_counts`.
   800	  std::vector<std::size_t> draw_from(group_count);
   801	  for (std::size_t g = 0; g < group_count; ++g) { draw_from[g] = group_rows[g].size(); }
   802	  if (params.match_decoy_candidate_counts)
   803	  {
   804	    std::vector<std::size_t> target_n;
   805	    std::vector<std::size_t> decoys;
   806	    for (std::size_t g = 0; g < group_count; ++g)
   807	    {
   808	      if (group_label[g] == 1) { target_n.push_back(group_rows[g].size()); }
   809	      else { decoys.push_back(g); }
   810	    }
   811	    if (!target_n.empty() && !decoys.empty())
   812	    {
   813	      std::sort(target_n.begin(), target_n.end());
   814	      // Both sides sorted by count, then matched by rank: the k-th smallest
   815	      // decoy takes the k-th smallest target's count. That maps the whole
   816	      // distribution rather than just its mean, and it is a function of the
   817	      // data alone, so the run stays deterministic.
   818	      std::stable_sort(decoys.begin(), decoys.end(),
   819	                       [&](std::size_t a, std::size_t b)
   820	                       { return group_rows[a].size() < group_rows[b].size(); });
   821	      for (std::size_t i = 0; i < decoys.size(); ++i)
   822	      {
   823	        const std::size_t q = i * target_n.size() / decoys.size();
   824	        draw_from[decoys[i]] = std::min(target_n[q], group_rows[decoys[i]].size());
   825	        if (draw_from[decoys[i]] == 0) { draw_from[decoys[i]] = 1; }
   826	      }
   827	    }
   828	  }
   829	
   830	
   831	  if (!frozen_applied)
   832	  {
   833	#ifdef _OPENMP
   834	#pragma omp parallel for schedule(dynamic, 1) reduction(+ : n_trained, n_skipped)
   835	#endif
   836	  for (int fold = 0; fold < folds; ++fold)
   837	  {
   838	    std::vector<std::size_t> train_rows;
   839	    std::vector<std::size_t> train_groups;
   840	    train_rows.reserve(n);
   841	    train_groups.reserve(group_count);
   842	    for (std::size_t g = 0; g < group_count; ++g)
   843	    {
   844	      if (!params.disable_cv && group_fold[g] == fold) { continue; }
   845	      train_groups.push_back(g);
   846	      train_rows.insert(train_rows.end(), group_rows[g].begin(), group_rows[g].end());
   847	    }
   848	    if (train_rows.empty()) { continue; }
   849	
   850	    // w = S_W^-1 (mu_pos - mu_neg): Fisher's discriminant on two explicit row sets.
   851	    auto fit_lda = [&](const std::vector<std::size_t>& positive_rows,
   852	                       const std::vector<std::size_t>& negative_rows,
   853	                       std::vector<double>& out) -> bool {
   854	      if (positive_rows.size() < 2 || negative_rows.size() < 2) { return false; }
   855	      std::vector<double> positive_mean(m, 0.0);
   856	      std::vector<double> negative_mean(m, 0.0);
   857	      for (const std::size_t row : positive_rows)
   858	      {
   859	        for (std::size_t j = 0; j < m; ++j) { positive_mean[j] += z[row][j]; }
   860	      }
   861	      for (const std::size_t row : negative_rows)
   862	      {
   863	        for (std::size_t j = 0; j < m; ++j) { negative_mean[j] += z[row][j]; }
   864	      }
   865	      for (std::size_t j = 0; j < m; ++j)
   866	      {
   867	        positive_mean[j] /= static_cast<double>(positive_rows.size());
   868	        negative_mean[j] /= static_cast<double>(negative_rows.size());
   869	      }
   870	
   871	      std::vector<double> covariance(m * m, 0.0);
   872	      auto add_within_class = [&](const std::vector<std::size_t>& rows,
   873	                                  const std::vector<double>& class_mean) {
   874	        for (const std::size_t row : rows)
   875	        {
   876	          for (std::size_t j = 0; j < m; ++j)
   877	          {
   878	            const double dj = z[row][j] - class_mean[j];
   879	            for (std::size_t k = 0; k <= j; ++k)
   880	            {
   881	              covariance[j * m + k] += dj * (z[row][k] - class_mean[k]);
   882	            }
   883	          }
   884	        }
   885	      };
   886	      add_within_class(positive_rows, positive_mean);
   887	      add_within_class(negative_rows, negative_mean);
   888	      const double dof =
   889	        static_cast<double>(positive_rows.size() + negative_rows.size() - 2);
   890	      for (std::size_t j = 0; j < m; ++j)
   891	      {
   892	        for (std::size_t k = 0; k <= j; ++k)
   893	        {
   894	          covariance[j * m + k] /= dof;
   895	          covariance[k * m + j] = covariance[j * m + k];
   896	        }
   897	      }
   898	
   899	      std::vector<double> difference(m);
   900	      for (std::size_t j = 0; j < m; ++j)
   901	      {
   902	        difference[j] = positive_mean[j] - negative_mean[j];
   903	      }
   904	
   905	      // A positive user ridge is tried exactly; a zero/negative ridge starts
   906	      // with a tiny numerical ridge. Increase it geometrically on failure.
   907	      double ridge = params.ridge > 0.0 && std::isfinite(params.ridge)
   908	                       ? params.ridge
   909	                       : 1e-12;
   910	      bool solved = false;
   911	      for (int attempt = 0; attempt < 10 && !solved; ++attempt)
   912	      {
   913	        std::vector<double> regularised = covariance;
   914	        for (std::size_t j = 0; j < m; ++j) { regularised[j * m + j] += ridge; }
   915	        solved = lda_detail::choleskySolve(regularised, difference, out);
   916	        ridge *= 10.0;
   917	      }
   918	      // A feature the caller declares lower-is-better may never earn a positive
   919	      // weight, however the fold's data happens to fall.
   920	      //
   921	      // This is DIA-NN's check_weights (diann.cpp:6592), which hard-clips the
   922	      // RT-deviation and mass-accuracy weights to <= 0 so that a larger
   923	      // deviation can never raise a score. The principle is the same wherever a
   924	      // feature has a known direction: without it the discriminant is free to
   925	      // learn, in-sample, that being further from the prediction is evidence
   926	      // FOR a peptide -- which fits the fold and generalises to nothing.
   927	      //
   928	      // Only the linear classifier is constrained here. A tree ensemble has no
   929	      // weight vector to clip, and enforcing the same thing on GBT needs
   930	      // monotone split constraints, which is a larger change. Since GBT is the
   931	      // default, this currently protects the non-default path.
   932	      for (const std::size_t j : params.nonpositive_features)
   933	      {
   934	        if (j < out.size() && out[j] > 0.0) { out[j] = 0.0; }
   935	      }
   936	
   937	      double norm_squared = 0.0;
   938	      for (const double value : out) { norm_squared += value * value; }
   939	      return solved && norm_squared > std::numeric_limits<double>::epsilon();
   940	    };
   941	
   942	    // The learner for this fold. LDA carries a weight vector; GBT carries a tree ensemble. Every
   943	    // score in this fold goes through score_row(), so the choice is made in exactly ONE place and
   944	    // the surrounding machinery -- training-set selection, fold normalisation, q-values -- is
   945	    // literally the same code for both.
   946	    GBT gbt;
   947	    std::vector<NNEnsemble> nn_members;      // one per bag; size 1 when bagging is off
   948	    const bool use_gbt = (params.classifier == LDAParams::Classifier::GBT);
   949	    const bool use_nn = (params.classifier == LDAParams::Classifier::NN);
   950	
   951	
   952	    // Initial direction: the signed feature with the largest absolute Welch
   953	    // two-sample t statistic between all target and decoy training rows.
   954	    std::vector<double> w(m, 0.0);
   955	    auto score_row = [&](std::size_t row) -> double {
   956	      if (use_nn && !nn_members.empty()) { return baggedScore(nn_members, z[row]); }
   957	      return (use_gbt && gbt.trained()) ? gbt.score(z[row]) : lda_detail::dot(w, z[row]);
   958	    };
   959	    // Fit whichever learner this run selected, on the given rows. `mask` is honoured by the NN
   960	    // only -- for a tree, a feature the fit never split on is already inert at scoring time, and
   961	    // for LDA a zero weight is likewise inert, so neither needs one.
   962	    auto fit_learner = [&](const std::vector<std::size_t>& pos, const std::vector<std::size_t>& neg,
   963	                           const std::vector<char>& mask) -> bool {
   964	      if (use_nn)
   965	      {
   966	        NNParams np = nn_params;
   967	        np.mask = mask;
   968	        if (params.bag_fraction >= 1.0)
   969	        {
   970	          NNEnsemble e;
   971	          if (!e.fit(z, pos, neg, np)) { return false; }
   972	          nn_members.assign(1, std::move(e));
   973	          return true;
   974	        }
   975	        if (!(params.bag_fraction > 0.0)) { return false; }   // an empty bag is not a trained model
   976	        AnchorTrainingParams ap;
   977	        ap.nn = np;
   978	        ap.bag_fraction = params.bag_fraction;
   979	        ap.min_anchors = 1;                  // the outer loop already refuses to fit on nothing
   980	        ap.seed = params.seed;
   981	        auto members = trainBaggedOnAnchors(z, pos, neg, group, ap);
   982	        if (members.empty()) { return false; }
   983	        nn_members = std::move(members);
   984	        return true;
   985	      }
   986	      if (use_gbt) { GBT g; if (!g.fit(z, pos, neg, gbt_params)) { return false; }
   987	                     gbt = std::move(g); return true; }
   988	      std::vector<double> next_w;
   989	      if (!fit_lda(pos, neg, next_w)) { return false; }
   990	      w.swap(next_w);
   991	      return true;
   992	    };
   993	
   994	    // Best-scoring row of each training group, computed in parallel.
   995	    //
   996	    // All three per-group scans below (seed selection, ranking, negative selection) have this
   997	    // shape, and for the NN each score_row() is a forward pass through TWELVE nets -- which made
   998	    // them the serial tail that held the whole routine to 763% CPU on a 224-core node while
   999	    // training itself was parallel.
  1000	    //
  1001	    // Output is PRE-SIZED and written by position, so nothing is appended and no ordering question
  1002	    // arises: out[i] is group train_groups[i]'s best row whatever order the iterations complete in.
  1003	    // The tie-break is `>` -- the FIRST maximum in the group's own row order wins -- which is what
  1004	    // the serial loops did, so the result is identical, not merely equivalent.
  1005	    auto best_rows_of = [&](const std::vector<std::size_t>& groups,
  1006	                            std::vector<std::size_t>& out, std::vector<double>& out_score) {
  1007	      out.assign(groups.size(), 0);
  1008	      out_score.assign(groups.size(), 0.0);
  1009	#ifdef _OPENMP
  1010	#pragma omp parallel for schedule(static) num_threads(inner_threads) if (groups.size() > 512)
  1011	#endif
  1012	      for (long long i = 0; i < static_cast<long long>(groups.size()); ++i)
  1013	      {
  1014	        const auto& rows_of_g = group_rows[groups[static_cast<std::size_t>(i)]];
  1015	        std::size_t best_row = rows_of_g.front();
  1016	        double best = score_row(best_row);
  1017	        for (const std::size_t row : rows_of_g)
  1018	        {
  1019	          const double sc = score_row(row);
  1020	          if (sc > best) { best = sc; best_row = row; }
  1021	        }
  1022	        out[static_cast<std::size_t>(i)] = best_row;
  1023	        out_score[static_cast<std::size_t>(i)] = best;
  1024	      }
  1025	    };
  1026	    std::vector<std::size_t> bg_row;
  1027	    std::vector<double> bg_score;
  1028	
  1029	    std::size_t best_feature = 0;
  1030	    double best_abs_t = -1.0;
  1031	    double best_difference = 1.0;
  1032	    for (std::size_t j = 0; j < m; ++j)
  1033	    {
  1034	      // The seed mask applies HERE TOO, not only to the fit. This bootstrap picks one feature, and
  1035	      // that single-feature direction is what ranks each group to choose its seed row. Searching a
  1036	      // masked column would let the excluded feature choose the training examples, and zeroing its
  1037	      // value at fit time cannot undo which rows it selected -- the leak mechanism 1 exists to
  1038	      // close, reopened one step upstream. Found by adversarial review, not by a test.
  1039	      if (!params.seed_mask.empty() && j < params.seed_mask.size() && !params.seed_mask[j]) { continue; }
  1040	      double sum[2] = {0.0, 0.0};
  1041	      double sum_sq[2] = {0.0, 0.0};
  1042	      std::size_t class_n[2] = {0, 0};
  1043	      for (const std::size_t row : train_rows)
  1044	      {
  1045	        const int cls = labels[row] == 1 ? 1 : 0;
  1046	        sum[cls] += z[row][j];
  1047	        sum_sq[cls] += z[row][j] * z[row][j];
  1048	        ++class_n[cls];
  1049	      }
  1050	      if (class_n[0] == 0 || class_n[1] == 0) { continue; }
  1051	      const double class_mean0 = sum[0] / static_cast<double>(class_n[0]);
  1052	      const double class_mean1 = sum[1] / static_cast<double>(class_n[1]);
  1053	      const double var0 = class_n[0] > 1
  1054	                            ? std::max(0.0, (sum_sq[0] - sum[0] * class_mean0) /
  1055	                                                static_cast<double>(class_n[0] - 1))
  1056	                            : 0.0;
  1057	      const double var1 = class_n[1] > 1
  1058	                            ? std::max(0.0, (sum_sq[1] - sum[1] * class_mean1) /
  1059	                                                static_cast<double>(class_n[1] - 1))
  1060	                            : 0.0;
  1061	      const double difference = class_mean1 - class_mean0;
  1062	      const double standard_error =
  1063	        std::sqrt(var0 / static_cast<double>(class_n[0]) +
  1064	                  var1 / static_cast<double>(class_n[1]));
  1065	      const double abs_t = standard_error > 0.0
  1066	                             ? std::abs(difference) / standard_error
  1067	                             : (difference == 0.0 ? 0.0
  1068	                                                  : std::numeric_limits<double>::infinity());
  1069	      if (abs_t > best_abs_t)
  1070	      {
  1071	        best_abs_t = abs_t;
  1072	        best_feature = j;
  1073	        best_difference = difference;
  1074	      }
  1075	    }
  1076	    w[best_feature] = best_difference < 0.0 ? -1.0 : 1.0;
  1077	
  1078	    // ...but a SINGLE feature is a far weaker seed than it looks, and when it is too weak the
  1079	    // semi-supervised loop never ignites: iteration 0 selects positives with this w, finds none at
  1080	    // q<=train_fdr_initial, skips the fit, so w is unchanged and every later iteration skips too.
  1081	    // Measured on testdata/lda_fixture.txt (44,568 real OpenSWATH rows / 2,500 precursors): the best
  1082	    // single feature reaches |t|=13.7, which sounds decisive but over 29,482 rows is a 0.16 sd
  1083	    // per-row effect -- and with ~18 candidate peak groups per precursor, the max-over-group that
  1084	    // actually does the ranking is dominated by extreme-value noise. Result: a FLAT q of 0.965 for
  1085	    // every target group, 0/9 iterations trained, and 0 IDs.
  1086	    //
  1087	    // The fix is to seed with a real multivariate direction. The target/decoy labels are known
  1088	    // outright -- no FDR estimate needed -- so an LDA of all top-target rows against all top-decoy
  1089	    // rows is available for free and is enormously stronger than one column. The target class is
  1090	    // contaminated (most target precursors are false), which shrinks the fitted direction toward
  1091	    // zero but does not rotate it systematically: the contaminant IS the decoy distribution, so it
  1092	    // biases mu_pos toward mu_neg and costs magnitude, not orientation. It only has to be good
  1093	    // enough to ignite the loop; iteration 0 then re-selects positives properly.
  1094	    // pyprophet has no equivalent problem because OpenSWATH hands it a composite `main_score` to
  1095	    // rank by (semi_supervised.py: train.rank_by("main_score")). ODIA has no such column.
  1096	    {
  1097	      std::vector<std::size_t> seed_pos, seed_neg;
  1098	      seed_pos.reserve(train_groups.size());
  1099	      seed_neg.reserve(train_groups.size());
  1100	      best_rows_of(train_groups, bg_row, bg_score);
  1101	      for (std::size_t i = 0; i < train_groups.size(); ++i)
  1102	      {
  1103	        (group_label[train_groups[i]] == 1 ? seed_pos : seed_neg).push_back(bg_row[i]);
  1104	      }
  1105	      // The seed must use the SAME learner as the loop. Seeding a GBT run with an LDA fit
  1106	      // reintroduces exactly the cold start this block exists to prevent, one level up: on data
  1107	      // whose signal is an interaction, the LDA seed is at chance by construction, so iteration 0
  1108	      // selects no positives, the GBT never fits, and the run silently falls back to ranking by a
  1109	      // hyperplane that cannot see the signal. Caught by odia_gbt_test T7, which reported 0 IDs
  1110	      // for BOTH classifiers before this.
  1111	      //
  1112	      // params.seed_mask applies HERE and only here (mechanism 1) -- the seed picks the anchors,
  1113	      // so it is the seed that must not see the calibrated features.
  1114	      if (!fit_learner(seed_pos, seed_neg, params.seed_mask))
  1115	      {
  1116	        // A failed seed is the cold start this block exists to prevent. Recording it means the
  1117	        // caller's "fitted 0 iterations" warning can fire instead of the run silently ranking by
  1118	        // the one-feature bootstrap while claiming a learned model.
  1119	        ++n_skipped;
  1120	      }
  1121	    }
  1122	
  1123	    std::vector<std::size_t> prev_positives;   // mechanism 5's state, sorted
  1124	    bool prev_is_main_threshold = false;       // provenance of prev_positives: selected at
  1125	                                               // train_fdr (iteration >= 1), never at the 0.15
  1126	                                               // initial cut. The 2026-09-01 smoke showed a
  1127	                                               // loop-counter guard is not provenance: any skip
  1128	                                               // desynchronises them and the cross-threshold
  1129	                                               // comparison comes back.
  1130	    bool fit_ok_last_iter = false;             // whether the PREVIOUS iteration refit the model.
  1131	                                               // A failed or skipped fit leaves the model -- and
  1132	                                               // therefore the next selection -- unchanged, and an
  1133	                                               // unchanged selection reads as Jaccard 1.0: a fold
  1134	                                               // whose fits keep failing would report "converged".
  1135	    AnchorTrainingReport stop_rep;             // PERSISTS across iterations: the high-water mark
  1136	                                               // and patience strikes are cumulative state, and a
  1137	                                               // per-iteration report silently reduces patience
  1138	                                               // to single-shot.
  1139	    std::vector<std::size_t> log_prev;         // iteration_log's own state, sorted -- deliberately
  1140	                                               // NOT shared with mechanism 5, so the diagnostic
  1141	                                               // can never perturb the stop decision
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
  1168	      // The churn diagnostic reads the selection BEFORE the too-few-positives skip below, so a
  1169	      // skipped iteration still shows its (tiny) positive set instead of vanishing from the log.
  1170	      // The iteration 0 -> 1 overlap spans the train_fdr_initial -> train_fdr threshold change and
  1171	      // is expected to be low on a HEALTHY run; it is printed rather than suppressed -- it is the
  1172	      // first comparison anyone asks about -- and marked so nobody reads it as a collapse.
  1173	      if (params.iteration_log)
  1174	      {
  1175	        std::vector<std::size_t> curr = positive_rows;
  1176	        std::sort(curr.begin(), curr.end());
  1177	        if (log_prev.empty())
  1178	        {
  1179	          std::fprintf(stderr, "[lda] fold %d iter %d: positives %zu jaccard --\n",
  1180	                       fold, iteration, curr.size());
  1181	        }
  1182	        else
  1183	        {
  1184	          std::fprintf(stderr, "[lda] fold %d iter %d: positives %zu jaccard %.4f%s\n",
  1185	                       fold, iteration, curr.size(), jaccardOverlap(log_prev, curr),
  1186	                       iteration == 1 ? " (crosses the initial->main train_fdr change)" : "");
  1187	        }
  1188	        log_prev.swap(curr);
  1189	      }
  1190	      // Mechanism 5, checked BEFORE this iteration's skips and fit. Three placement consequences,
  1191	      // each fixing a reviewed defect of the post-fit version (analysis77, 2026-09-01/02 reviews):
  1192	      //  * a collapse now breaks BEFORE a model is trained on the collapsed selection, so the
  1193	      //    previous model is what survives;
  1194	      //  * a catastrophic shrink below m+2 is seen by the comparison instead of being skipped
  1195	      //    around (the too-few-positives `continue` used to bypass the stop entirely, i.e. the
  1196	      //    rule went blind exactly at the largest collapse);
  1197	      //  * comparisons are gated on PROVENANCE (prev selected at train_fdr) and on the previous
  1198	      //    iteration having actually refit -- a failed/skipped fit repeats the same selection,
  1199	      //    and an unchanged selection is Jaccard 1.0, which would read as "converged".
  1200	      if (params.stop_on_composition)
  1201	      {
  1202	        std::vector<std::size_t> curr = positive_rows;
  1203	        std::sort(curr.begin(), curr.end());
  1204	        bool go = true;
  1205	        if (prev_is_main_threshold && iteration >= 1 && fit_ok_last_iter)
  1206	        {
  1207	          AnchorTrainingParams ap;
  1208	          ap.stop_jaccard = params.stop_jaccard;
  1209	          // n_iter + 1, NOT n_iter: the for-loop's own bound is the cap here, and the helper's
  1210	          // pre-fit cap check double-counted it -- an armed k12 run returned "cap" at iteration
  1211	          // 11 BEFORE fit 11 ran, silently delivering a k11 model under a k12 label (codex S1).
  1212	          // With the helper's cap unreachable, cap termination is the loop's natural exit, and
  1213	          // the verdict for it prints after the loop.
  1214	          ap.max_iterations = params.n_iter + 1;
  1215	          ap.shrink_floor = params.stop_shrink_floor;
  1216	          ap.stop_patience = params.stop_patience;
  1217	          stop_rep.iterations_run = iteration;
  1218	          go = anchorIterationShouldContinue(prev_positives, curr, ap, stop_rep);
  1219	        }
  1220	        else if (stop_rep.collapse_strikes >= 1 && positive_rows.size() < m + 2)
  1221	        {
  1222	          // A catastrophic collapse starves the patience of its second strike: the tiny set
  1223	          // skips the fit, the fit gate then blocks every later comparison (the frozen model
  1224	          // reproduces the same selection forever), and the loop would burn to the cap with no
  1225	          // verdict -- silence exactly at the tail event this mechanism exists to catch
  1226	          // (kimi F1 / codex S1, found before any armed run shipped). A starved selection
  1227	          // arriving with a floor breach already on record IS the second strike.
  1228	          stop_rep.collapsed = true;
  1229	          stop_rep.note = "positive set starved below m+2 (" + std::to_string(positive_rows.size()) +
  1230	                          ") with a floor breach already on record";
  1231	          go = false;
  1232	        }
  1233	        prev_positives.swap(curr);
  1234	        prev_is_main_threshold = (iteration >= 1);
  1235	        if (!go)
  1236	        {
  1237	          // The stop VERDICT is part of the run's output contract: without this line the b12-style
  1238	          // arm cannot say whether it converged or collapsed and the readout has to be
  1239	          // reverse-engineered from the churn log.
  1240	          std::fprintf(stderr, "[lda] fold %d iter %d: STOP %s -- %s\n", fold, iteration,
  1241	                       stop_rep.converged ? "converged" : "collapsed", stop_rep.note.c_str());
  1242	          break;
  1243	        }
  1244	      }
  1245	      // From here to the fit, every early exit means "the model did not change this iteration" --
  1246	      // recorded so the next iteration's stop comparison knows its selection is a repeat.
  1247	      fit_ok_last_iter = false;
  1248	      // Too few confident positives to fit an m-dimensional discriminant. Skipping is right, but it
  1249	      // used to be SILENT -- and silence here is dangerous: if every iteration skips, `w` stays at
  1250	      // its initialisation (a single feature, weight +/-1), so the "LDA" degenerates to ranking by
  1251	      // one sub-score and nothing says so. Record it; the caller reports it.
  1252	      if (positive_rows.size() < m + 2)
  1253	      {
  1254	        ++n_skipped;
  1255	        continue;
  1256	      }
  1257	      // n_trained is incremented AFTER the fit succeeds, further down -- not here. Counting it at
  1258	      // selection time reported a trained iteration for a fit that then failed (an empty bag from
  1259	      // bag_fraction <= 0, a singular covariance), and the scores in that case come from the
  1260	      // previous model or the one-feature fallback.
  1261	
  1262	      // The positive class above is one row per precursor -- `candidate.best_row`, the highest
  1263	      // scoring peak group. Taking ALL decoy rows as negatives would therefore make the two classes
  1264	      // differ in two ways at once: target vs decoy (wanted) AND rank-1 vs runner-up (not wanted).
  1265	      // With ~4-6 candidate peak groups per precursor the negative class is then dominated by
  1266	      // runner-ups, so the fitted direction partly separates "best peak in its group" from "not the
  1267	      // best peak" -- a real, learnable axis that carries no target/decoy information. Both the
  1268	      // negative mean and the within-class covariance are pulled by it. Matching the positive side
  1269	      // (top peak per precursor) removes the confound; this is what pyprophet does
  1270	      // (semi_supervised.py: td_peaks = train.get_top_decoy_peaks()).
  1271	      std::vector<std::size_t> negative_rows;
  1272	      if (!params.top_decoys_only)
  1273	      {
  1274	        for (const std::size_t g : train_groups)
  1275	        {
  1276	          if (group_label[g] == 1) { continue; }
  1277	          negative_rows.insert(negative_rows.end(), group_rows[g].begin(), group_rows[g].end());
  1278	        }
  1279	      }
  1280	      else
  1281	      {
  1282	        // The ranking scan above already computed every training group's best row with THIS SAME
  1283	        // model -- nothing refits in between -- so re-scanning the decoys was pure duplicated work
  1284	        // on the most expensive operation in the loop. Reuse it.
  1285	        for (std::size_t i = 0; i < train_groups.size(); ++i)
  1286	        {
  1287	          if (group_label[train_groups[i]] == 1) { continue; }
  1288	          negative_rows.push_back(bg_row[i]);
  1289	        }
  1290	      }
  1291	      if (negative_rows.size() < 2) { continue; }
  1292	
  1293	      // A failed fit leaves the previous model in place, exactly as a failed Cholesky leaves the
  1294	      // previous w -- the iteration is skipped, not replaced with something degenerate. And it is
  1295	      // COUNTED as skipped, so a run whose every iteration failed cannot report itself as trained.
  1296	      if (fit_learner(positive_rows, negative_rows, {}))
  1297	      { ++n_trained; fit_ok_last_iter = true; }
  1298	      else { ++n_skipped; }
  1299	    }
  1300	    // Every armed termination gets a verdict, not only the explicit breaks: cap exhaustion,
  1301	    // k <= 2 (no same-threshold pair ever forms), and a final comparison gated out by a failed
  1302	    // fit all used to end in silence, and "no verdict line" is indistinguishable from a broken
  1303	    // log (kimi F4 / codex "every armed termination: not printed").
  1304	    if (params.stop_on_composition && !stop_rep.converged && !stop_rep.collapsed)
  1305	    {
  1306	      std::fprintf(stderr,
  1307	                   "[lda] fold %d: STOP cap -- ran all %d iterations without a convergence or "
  1308	                   "collapse verdict (final positives %zu, last flow +%zu/-%zu)\n",
  1309	                   fold, std::max(0, params.n_iter), prev_positives.size(),
  1310	                   stop_rep.last_entries, stop_rep.last_exits);
  1311	    }
  1312	
  1313	    // This model has seen no row from the groups scored.
  1314	    //
  1315	    // Parallel over groups: each iteration writes result.dscore at indices belonging to ITS OWN
  1316	    // group and reads only the (now fixed) model, so there is no reduction, no shared accumulator,
  1317	    // and no ordering question -- the output is bit-identical to the serial loop at any thread
  1318	    // count. Worth doing because for the NN this is one forward pass per row through 12 nets, and
  1319	    // it was the serial tail of an otherwise parallel routine.
  1320	    std::vector<std::size_t> score_groups;
  1321	    score_groups.reserve(group_count / static_cast<std::size_t>(folds) + 1);
  1322	    for (std::size_t g = 0; g < group_count; ++g)
  1323	    {
  1324	      if (group_fold[g] == fold) { score_groups.push_back(g); }
  1325	    }
  1326	#ifdef _OPENMP
  1327	#pragma omp parallel for schedule(static) num_threads(inner_threads) if (score_groups.size() > 1024)
  1328	#endif
  1329	    for (long long k = 0; k < static_cast<long long>(score_groups.size()); ++k)
  1330	    {
  1331	      for (const std::size_t row : group_rows[score_groups[static_cast<std::size_t>(k)]])
  1332	      {
  1333	        result.dscore[row] = score_row(row);
  1334	      }
  1335	    }
  1336	
  1337	    // Every fold has its OWN weight vector, and an LDA direction is defined only up to scale (and
  1338	    // the mean-offset depends on that fold's training set), so fold A's score of 4.0 and fold B's
  1339	    // of 4.0 mean different things. The final q-values below pool all folds into ONE ranking, which
  1340	    // silently assumes they are commensurable. They are not, and the pooled ranking is then partly
  1341	    // sorted by which fold a precursor happened to land in. Rescaling each fold's held-out scores to
  1342	    // that fold's own DECOY null (mean 0, sd 1) puts every fold on the one scale target-decoy FDR
  1343	    // actually cares about -- "how many decoy sds above the null" -- and makes pooling valid.
  1344	    // pyprophet sidesteps the problem differently: it averages the fold weight vectors into a single
  1345	    // model and rescores everything with it (at the cost of the leakage-free property kept here).
  1346	    // Keep this fold's trained ensemble for `model_out`. Fold-indexed slots are
  1347	    // disjoint across the parallel loop, so no synchronisation is needed.
  1348	    if (use_gbt && gbt.trained() && static_cast<std::size_t>(fold) < fold_models.size())
  1349	    { fold_models[static_cast<std::size_t>(fold)] = gbt; }
  1350	    if (params.normalize_folds)
  1351	    {
  1352	      double sum = 0.0, sum_sq = 0.0;
  1353	      std::size_t decoy_n = 0;
  1354	      for (std::size_t g = 0; g < group_count; ++g)
  1355	      {
  1356	        if (group_fold[g] != fold || group_label[g] == 1) { continue; }
  1357	        std::size_t best_row = group_rows[g].front();
  1358	        for (const std::size_t row : group_rows[g])
  1359	        {
  1360	          if (result.dscore[row] > result.dscore[best_row]) { best_row = row; }
  1361	        }
  1362	        sum += result.dscore[best_row];
  1363	        sum_sq += result.dscore[best_row] * result.dscore[best_row];
  1364	        ++decoy_n;
  1365	      }
  1366	      double mu = 0.0, sigma = 0.0;
  1367	      bool affine_ok = false;
  1368	      if (decoy_n >= 2)
  1369	      {
  1370	        mu = sum / static_cast<double>(decoy_n);
  1371	        const double var = std::max(0.0, (sum_sq - sum * mu) / static_cast<double>(decoy_n - 1));
  1372	        sigma = std::sqrt(var);
  1373	        if (static_cast<std::size_t>(fold) < fold_norm_mu.size())
  1374	        {
  1375	          // Saved with the fold model so a frozen application normalises
  1376	          // exactly as this fold did (a degenerate sigma is saved as-is and
  1377	          // the loader applies the same validity rule below).
  1378	          fold_norm_mu[static_cast<std::size_t>(fold)] = mu;
  1379	          fold_norm_sigma[static_cast<std::size_t>(fold)] = sigma;
  1380	        }
  1381	        // A degenerate null (all decoys identical) carries no scale information; leaving those
  1382	        // scores unscaled is the only honest option, and shifting them alone would be worse.
  1383	        affine_ok = (sigma > std::numeric_limits<double>::epsilon() && std::isfinite(sigma));
  1384	      }
  1385	      if (params.fold_pool_rank)
  1386	      {
  1387	        // v1.13: within-fold rank pooling. Knots = every group's best in THIS fold (targets over all
  1388	        // rows, decoys over their drawn prefix -- exactly what assignQValues ranks), sorted descending.
  1389	        // A group's best row lands on its knot: frac = (number of knots >= v) / n_f, so tied groups
  1390	        // share the MAX rank (the conservative side; verified on F17: the equality branch below is
  1391	        // unreachable because upper_bound(greater<>) stops at the first knot < v, and the
  1392	        // interpolation branch then yields exactly k/n_f); other rows interpolate between
  1393	        // adjacent knots so the within-group order is kept.
  1394	        std::vector<double> knots;
  1395	        for (std::size_t g = 0; g < group_count; ++g)
  1396	        {
  1397	          if (group_fold[g] != fold) { continue; }
  1398	          const std::size_t take = (group_label[g] == 1) ? group_rows[g].size()
  1399	                                                         : std::min(draw_from[g], group_rows[g].size());
  1400	          double best = -std::numeric_limits<double>::infinity();
  1401	          for (std::size_t k = 0; k < take; ++k) { best = std::max(best, result.dscore[group_rows[g][k]]); }
  1402	          if (std::isfinite(best)) { knots.push_back(best); }
  1403	        }
  1404	        std::sort(knots.begin(), knots.end(), std::greater<double>());
  1405	        const std::size_t n_f = knots.size();
  1406	        if (n_f >= 2)
  1407	        {
  1408	          const double nf = static_cast<double>(n_f);
  1409	          for (std::size_t g = 0; g < group_count; ++g)
  1410	          {
  1411	            if (group_fold[g] != fold) { continue; }
  1412	            for (const std::size_t row : group_rows[g])
  1413	            {
  1414	              const double v = result.dscore[row];
  1415	              const std::size_t k = static_cast<std::size_t>(
  1416	                std::upper_bound(knots.begin(), knots.end(), v, std::greater<double>()) - knots.begin());
  1417	              double frac;
  1418	              if (k < n_f && v == knots[k]) { frac = static_cast<double>(k + 1) / nf; }
  1419	              else if (k == 0)
  1420	              { frac = (1.0 - (v - knots[0]) / (knots[0] - knots[1] + 1e-12)) / nf; }
  1421	              else if (k < n_f)
  1422	              { frac = (static_cast<double>(k) + (knots[k - 1] - v) / (knots[k - 1] - knots[k] + 1e-300)) / nf; }
  1423	              else
  1424	              { frac = 1.0 + (knots[n_f - 1] - v) / (1.0 + knots[n_f - 2] - knots[n_f - 1]); }
  1425	              frac = std::max(frac, 1e-12);
  1426	              result.dscore[row] = -std::log10(frac);
  1427	            }
  1428	          }
  1429	        }
  1430	      }
  1431	      else if (affine_ok)
  1432	      {
  1433	        for (std::size_t g = 0; g < group_count; ++g)
  1434	        {
  1435	          if (group_fold[g] != fold) { continue; }
  1436	          for (const std::size_t row : group_rows[g])
  1437	          {
  1438	            result.dscore[row] = (result.dscore[row] - mu) / sigma;
  1439	          }
  1440	        }
  1441	      }
  1442	    }
  1443	  }
  1444	
  1445	  }   // if (!frozen_applied): the training path
  1446	
  1447	  // ---- SAVE, gbt engine only ----
  1448	  // The fold ensemble plus everything a frozen application needs to score a
  1449	  // NEW run identically to how this one scored its own held-out folds: this
  1450	  // run's standardisation (mean, sd per feature) and each fold's decoy
  1451	  // normalisation (mu, sigma). Written after every fold has finished.
  1452	  if (!frozen_applied && !model_out.empty())
  1453	  {
  1454	    if (!gbt_engine)
  1455	    {
  1456	      std::fprintf(stderr, "[lda] -classifier_model_out is honoured by the gbt and percolator "
  1457	                           "engines only; nothing written\n");
  1458	    }
  1459	    else
  1460	    {
  1461	      // Written to a sibling temporary and renamed into place, so a failed or
  1462	      // interrupted save never leaves a partial file that a later
  1463	      // -classifier_model_in would refuse (or, worse, half-read).
  1464	      const std::string tmp = model_out + ".partial";
  1465	      std::size_t k_have = 0;
  1466	      for (const auto& g : fold_models) { k_have += g.trained() ? 1 : 0; }
  1467	      bool written = false;
  1468	      if (k_have > 0)
  1469	      {
  1470	        std::ofstream os(tmp);
  1471	        os.imbue(std::locale::classic());
  1472	        os.precision(17);
  1473	        os << "ODIA-GBT-FOLDS 1 " << m << ' ' << fold_models.size() << '\n';
  1474	        // The sub-score names travel with the model; the loader refuses a
  1475	        // run whose sub-scores differ, whatever their count.
  1476	        os << m;
  1477	        for (std::size_t j = 0; j < m; ++j)
  1478	        {
  1479	          os << ' ' << ((feature_names != nullptr && j < feature_names->size())
  1480	                          ? (*feature_names)[j] : std::string("var_") + std::to_string(j));
  1481	        }
  1482	        os << '\n';
  1483	        for (std::size_t j = 0; j < m; ++j) { os << (j ? " " : "") << mean[j]; }
  1484	        os << '\n';
  1485	        for (std::size_t j = 0; j < m; ++j) { os << (j ? " " : "") << sd[j]; }
  1486	        os << '\n';
  1487	        for (std::size_t f = 0; f < fold_models.size(); ++f)
  1488	        {
  1489	          const bool has = fold_models[f].trained();
  1490	          os << (has ? 1 : 0) << ' ' << fold_norm_mu[f] << ' ' << fold_norm_sigma[f] << '\n';
  1491	          if (has) { fold_models[f].save(os); }
  1492	        }
  1493	        // The precursor-to-fold assignment. A frozen application scores a
  1494	        // precursor this ensemble TRAINED ON with the one fold that excluded
  1495	        // it -- its native held-out score -- and only a precursor the ensemble
  1496	        // never saw with the fold mean. Without this, re-scoring the saving
  1497	        // run hands every row K-1 models that memorised it, and a control
  1498	        // scored that way is optimistic against any arm it is compared to.
  1499	        os << "FOLDMAP " << group_count << '\n';
  1500	        for (std::size_t g = 0; g < group_count; ++g)
  1501	        { os << group_id[g] << ' ' << group_fold[g] << '\n'; }
  1502	        os.flush();
  1503	        written = static_cast<bool>(os);
  1504	        os.close();
  1505	        written = written && !os.fail();
  1506	        if (written) { written = (std::rename(tmp.c_str(), model_out.c_str()) == 0); }
  1507	        if (!written) { std::remove(tmp.c_str()); }
  1508	      }
  1509	      if (written)
  1510	      {
  1511	        std::fprintf(stderr, "[gbt] TRAINED -> %s: %zu of %zu fold models saved, %zu features\n",
  1512	                     model_out.c_str(), k_have, fold_models.size(), m);
  1513	      }
  1514	      else
  1515	      {
  1516	        std::fprintf(stderr, "[gbt] TRAINED -> %s: NOT written (%zu trained folds%s)\n",
  1517	                     model_out.c_str(), k_have, k_have > 0 ? ", write or rename failed" : "");
  1518	      }
  1519	    }
  1520	  }
  1521	
  1522	  result.n_iterations_trained = n_trained;
  1523	  result.n_iterations_skipped = n_skipped;
  1524	
  1525	  std::vector<lda_detail::RankedGroup> final_ranked;
  1526	  final_ranked.reserve(group_count);
  1527	  for (std::size_t g = 0; g < group_count; ++g)
  1528	  {
  1529	    std::size_t best_row = group_rows[g].front();
  1530	    const std::size_t take = std::min(draw_from[g], group_rows[g].size());
  1531	    for (std::size_t k = 0; k < take; ++k)
  1532	    {
  1533	      const std::size_t row = group_rows[g][k];
  1534	      if (result.dscore[row] > result.dscore[best_row]) { best_row = row; }
  1535	    }
  1536	    final_ranked.push_back(
  1537	      {g, best_row, group_label[g], result.dscore[best_row], 1.0});
  1538	  }
  1539	  lda_detail::assignQValues(final_ranked, params.use_pi0);
  1540	  for (const auto& ranked_group : final_ranked)
  1541	  {
  1542	    for (const std::size_t row : group_rows[ranked_group.group_index])
  1543	    {
  1544	      result.qvalue[row] = ranked_group.qvalue;
  1545	      result.pvalue[row] = ranked_group.pvalue;
  1546	      result.pep[row]    = ranked_group.pep;
  1547	    }
  1548	  }
  1549	  return result;
  1550	}
  1551	
  1552	} // namespace ODIA::Scoring
  1553	
  1554	#endif // ODIA_LDA_H
```
### include/odia/scoring/gbt.h (whole file)
```
     1	// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
     2	// SPDX-License-Identifier: BSD-3-Clause
     3	//
     4	// Adopted from okohlbacher/OpenDIAlyzer, branch odia-engine-and-openms-boundary,
     5	// tag odia-v0.3.0 (0ed5fb2), file src/odia_gbt.h. Same author, same BSD-3 licence.
     6	//
     7	// Changed here: namespace `odia` -> `ODIA::Scoring` to match this project, and
     8	// the sibling includes rewritten to this tree's layout. The algorithms are
     9	// untouched.
    10	//
    11	// Provenance that travels with the code, per that repository's CLEAN-ROOM.md:
    12	// this is an independent reimplementation of a published confidence model
    13	// (Reiter et al. mProphet, Nat Methods 2011; Teleman et al. pyprophet,
    14	// Bioinformatics 2015). No DIA-NN code is present in it.
    15	// odia_gbt.h — histogram gradient-boosted trees, as a drop-in replacement for the LDA's linear
    16	// discriminant inside the same semi-supervised loop.
    17	//
    18	// WHY NOT XGBoost/LightGBM. Neither is available to link here (checked: nothing in the OpenMS 3.6
    19	// ML/ tree, nothing in the environment), and vendoring one would put a large build dependency into a
    20	// tool whose scoring layer is otherwise three standalone headers with self-checks. A histogram GBT
    21	// is ~400 lines, has no dependencies, and — unlike a vendored library with its own threading and
    22	// RNG — is trivially deterministic, which matters because these scores feed an FDR estimate.
    23	//
    24	// WHY BOTHER. The LDA fits ONE hyperplane through the sub-score space. OpenSWATH's sub-scores
    25	// interact non-linearly in ways a plane cannot express: a high xcorr_shape means something very
    26	// different at high vs low signal-to-noise, and library correlation matters more when many
    27	// transitions are present. mokapot (Fondrie & Noble 2021) is exactly this substitution — Percolator
    28	// with boosted trees instead of the linear model — and reports +15% modified PSMs, +19% peptides,
    29	// +11% proteins at 1% FDR. pyprophet ships XGBoost as a non-default classifier for the same reason.
    30	// Expect less on vanilla tryptic bulk DIA than those headline numbers, which come from
    31	// immunopeptidomics and PTM work; see docs/OpenDIAlyzer-scoring-fdr-backlog.md, which flags exactly
    32	// this and says to read the gains skeptically.
    33	//
    34	// WHAT IT IS. Binary logistic gradient boosting, second-order (Newton) leaf values, feature values
    35	// pre-binned so split finding is O(bins) per feature per node rather than O(rows log rows):
    36	//
    37	//   - bin edges from quantiles of the TRAINING rows only, computed once per fit
    38	//   - level-wise growth to a fixed depth (not leaf-wise): bounded tree shape, no surprise depth
    39	//   - gain = 1/2 [ GL^2/(HL+lambda) + GR^2/(HR+lambda) - (GL+GR)^2/(HL+HR+lambda) ] - gamma
    40	//   - leaf value = -G/(H+lambda), shrunk by the learning rate
    41	//   - no row or column subsampling, no early stopping on a random split -> bit-identical every run
    42	//
    43	// DETERMINISM IS A REQUIREMENT, NOT A PREFERENCE. These scores produce q-values. A classifier whose
    44	// output moves between runs makes the FDR irreproducible, so every source of randomness that
    45	// ordinary GBT implementations use for regularisation is deliberately absent here. Ties in split
    46	// selection break on the lowest (feature, bin) index.
    47	
    48	#ifndef ODIA_GBT_H
    49	#define ODIA_GBT_H
    50	
    51	#ifdef _OPENMP
    52	#include <omp.h>
    53	#endif
    54	
    55	#include <algorithm>
    56	#include <cmath>
    57	#include <cstddef>
    58	#include <cstdint>
    59	#include <ios>
    60	#include <istream>
    61	#include <limits>
    62	#include <locale>
    63	#include <ostream>
    64	#include <string>
    65	#include <utility>
    66	#include <numeric>
    67	#include <vector>
    68	
    69	namespace ODIA::Scoring
    70	{
    71	
    72	struct GBTParams
    73	{
    74	  int n_trees = 120;              ///< boosting rounds
    75	  int max_depth = 4;              ///< level-wise depth; 4 gives <=16 leaves, enough for pairwise
    76	                                  ///< interactions between sub-scores without memorising rows
    77	  int n_bins = 64;                ///< histogram resolution per feature
    78	  double learning_rate = 0.1;     ///< shrinkage
    79	  double lambda = 1.0;            ///< L2 on leaf values
    80	  double gamma = 0.0;             ///< minimum gain to split
    81	  double min_child_weight = 1.0;  ///< minimum summed hessian in a child
    82	  int min_child_rows = 20;        ///< minimum rows in a child; guards tiny leaves on small folds
    83	  double max_delta_step = 0.0;    ///< XGBoost's max_delta_step: a leaf's Newton step G/(H+lambda) is
    84	                                  ///< clipped to +/-this BEFORE the learning rate. 0 = off (the native
    85	                                  ///< unbounded step). Bounds the first rounds of a fit whose prior is
    86	                                  ///< extreme, where a pure-positive leaf's hessian is tiny against its
    87	                                  ///< gradient and the unbounded step is |G|/lambda. EXPERIMENT knob.
    88	  bool intercept_zero = false;    ///< Intercept 0 (XGBoost base_score 0.5) instead of logit(training prior).
    89	                                  ///< At a 1:155 prior a pure-positive leaf's Newton step is ~1/p (~160);
    90	                                  ///< at intercept 0 it is ~2. EXPERIMENT knob (v1.11).
    91	  int warmup_rounds = 0;          ///< Learning-rate ramp: round r (1-based) uses lr * min(1, r/warmup_rounds);
    92	                                  ///< 0 = off (constant lr, identical arithmetic). EXPERIMENT knob (v1.12).
    93	  bool clip_gain = false;         ///< With max_delta_step > 0, evaluate split gains with the CLIPPED leaf
    94	                                  ///< weight (XGBoost CalcGain semantics) instead of the unclipped G^2/(H+lambda),
    95	                                  ///< so splits are ranked by the improvement they can actually realise.
    96	                                  ///< Identical arithmetic when off. EXPERIMENT knob (v1.11).
    97	  /// Bin edges from ALL rows rather than from the current training set.
    98	  ///
    99	  /// The default recomputes quantile edges from `pos + neg`, which in a
   100	  /// semi-supervised loop CHANGES every iteration. That is an amplification
   101	  /// path: adding one column changes the model, which changes the selected
   102	  /// positives, which re-discretises EVERY OTHER column -- so a perturbation
   103	  /// that should be local becomes global. Measured consequence: a column of
   104	  /// pure noise moves identifications by up to 11.9%, with Spearman rho 0.79
   105	  /// between a run and the same run plus that column.
   106	  ///
   107	  /// Edges over all rows are label-independent and therefore identical across
   108	  /// iterations and across feature sets, which closes the loop. Costs one extra
   109	  /// pass over the unlabelled rows per fit.
   110	  bool fixed_bins = false;
   111	  /// Threads for the histogram pass (0 = whatever OpenMP gives). fit() is normally called from
   112	  /// inside the fold loop's parallel region, where a nested team defaults to ONE thread -- so the
   113	  /// caller must both set this and enable a second active level, or the parallelism does nothing.
   114	  /// Results do NOT depend on this value; see the chunking note in growTree_.
   115	  int n_threads = 0;
   116	
   117	  /// pyProphet 3.0.15's XGBoost 3.2.0 configuration, verbatim.
   118	  ///
   119	  /// Six of the nine parameters it sets already match this struct's defaults --
   120	  /// lambda 1, gamma 0, min_child_weight 1, and subsample / colsample_bytree /
   121	  /// alpha at 1/1/0, which are no-ops and therefore equivalent to this
   122	  /// implementation not having them. The algorithm is already XGBoost's: same
   123	  /// second-order Newton leaf values, same split gain.
   124	  ///
   125	  /// The two that differ are the two that matter, and both make the model
   126	  /// stronger and more prone to overfitting: depth 4 -> 6 takes the leaf count
   127	  /// from <=16 to <=64, and eta 0.1 -> 0.3 triples the step. This project's
   128	  /// defaults were chosen for "pairwise interactions between sub-scores without
   129	  /// memorising rows" with 10 features; we now have 15.
   130	  static GBTParams pyprophet()
   131	  {
   132	    GBTParams p;
   133	    p.max_depth = 6;          // XGBoost "max_depth": 6
   134	    p.learning_rate = 0.3;    // XGBoost "eta": 0.3
   135	    p.lambda = 1.0;           // "lambda": 1   -- already the default here
   136	    p.gamma = 0.0;            // "gamma": 0
   137	    p.min_child_weight = 1.0; // "min_child_weight": 1
   138	    return p;
   139	  }
   140	};
   141	
   142	namespace gbt_detail
   143	{
   144	
   145	/// Gradient and hessian for one histogram bin, INTERLEAVED.
   146	///
   147	/// They were two separate arrays. Every accumulation touches the same index k in both, and k is the
   148	/// data-dependent bin, so the access is effectively random: two arrays meant two cache lines pulled
   149	/// in per (row, feature) in the loop that is 92.6% of a fit -- 13.9e9 potentially-missing accesses
   150	/// over a fit. Interleaved, g and h are 16 adjacent bytes and land in ONE line, halving the misses
   151	/// for identical arithmetic in identical order (so the result is bit-unchanged).
   152	struct GH
   153	{
   154	  double g = 0.0;
   155	  double h = 0.0;
   156	  /// Rows in this bin. Needed to enforce min_child_rows PER CHILD: without it
   157	  /// `nL` below has nothing to accumulate, which is exactly how the predecessor
   158	  /// left the same guard dead (its handoff, open defect 2: "nL is declared,
   159	  /// never incremented, discarded with (void)nL ... leaves of ~4 rows are
   160	  /// possible -- an overfitting surface in exactly the score tail that sets the
   161	  /// FDR threshold").
   162	  std::size_t n = 0;
   163	};
   164	
   165	/// One level-wise regression tree over pre-binned features.
   166	struct Tree
   167	{
   168	  /// Nodes are stored as a complete binary tree indexed 1..2^(d+1)-1 so children of i are 2i, 2i+1;
   169	  /// this costs a little memory at depth 4 (31 slots) and removes all pointer chasing.
   170	  std::vector<int> feature;      ///< split feature, or -1 for a leaf
   171	  std::vector<uint8_t> bin;      ///< split threshold: go left when bin <= this
   172	  std::vector<double> value;     ///< leaf value (valid where feature == -1)
   173	  int max_depth = 0;
   174	
   175	  double predict(const uint8_t* row) const
   176	  {
   177	    int node = 1;
   178	    for (int d = 0; d < max_depth; ++d)
   179	    {
   180	      if (feature[node] < 0) { break; }
   181	      node = (row[feature[node]] <= bin[node]) ? 2 * node : 2 * node + 1;
   182	    }
   183	    return value[node];
   184	  }
   185	};
   186	
   187	/// Quantile bin edges per feature, from the training rows. Returns edges[f] with <= n_bins-1
   188	/// thresholds; a value goes to bin b when it is <= edges[f][b] (last bin catches the rest).
   189	inline std::vector<std::vector<double>> computeBinEdges(
   190	  const std::vector<std::vector<double>>& X, const std::vector<std::size_t>& rows, int n_bins)
   191	{
   192	  const std::size_t m = X.empty() ? 0 : X[0].size();
   193	  std::vector<std::vector<double>> edges(m);
   194	  std::vector<double> col;
   195	  col.reserve(rows.size());
   196	  for (std::size_t f = 0; f < m; ++f)
   197	  {
   198	    col.clear();
   199	    for (const std::size_t r : rows)
   200	    {
   201	      const double v = X[r][f];
   202	      if (std::isfinite(v)) { col.push_back(v); }
   203	    }
   204	    if (col.empty()) { continue; }
   205	    std::sort(col.begin(), col.end());
   206	    // Distinct quantile cut points. Duplicates are dropped, so a feature with few distinct values
   207	    // simply gets few bins instead of many identical ones (which would waste histogram slots and
   208	    // create splits that separate nothing).
   209	    for (int b = 1; b < n_bins; ++b)
   210	    {
   211	      const std::size_t idx =
   212	        static_cast<std::size_t>(static_cast<double>(col.size()) * b / n_bins);
   213	      const double e = col[std::min(idx, col.size() - 1)];
   214	      if (edges[f].empty() || e > edges[f].back()) { edges[f].push_back(e); }
   215	    }
   216	  }
   217	  return edges;
   218	}
   219	
   220	/// Map every row to its bin index per feature. Non-finite values go to the LAST bin, so "missing"
   221	/// is a value the tree can split on rather than something silently imputed to a mean.
   222	inline std::vector<std::vector<uint8_t>> binFeatures(
   223	  const std::vector<std::vector<double>>& X, const std::vector<std::vector<double>>& edges)
   224	{
   225	  const std::size_t n = X.size(), m = edges.size();
   226	  std::vector<std::vector<uint8_t>> B(n, std::vector<uint8_t>(m, 0));
   227	  for (std::size_t i = 0; i < n; ++i)
   228	  {
   229	    for (std::size_t f = 0; f < m; ++f)
   230	    {
   231	      const double v = X[i][f];
   232	      if (!std::isfinite(v))
   233	      {
   234	        B[i][f] = static_cast<uint8_t>(edges[f].size());
   235	        continue;
   236	      }
   237	      const auto it = std::lower_bound(edges[f].begin(), edges[f].end(), v);
   238	      B[i][f] = static_cast<uint8_t>(it - edges[f].begin());
   239	    }
   240	  }
   241	  return B;
   242	}
   243	
   244	} // namespace gbt_detail
   245	
   246	/// Gradient-boosted trees for binary classification, trained on explicit positive/negative row sets
   247	/// so it plugs into the same semi-supervised loop as the LDA.
   248	class GBT
   249	{
   250	public:
   251	  /// Fit on rows `pos` (label 1) and `neg` (label 0) of X. Returns false if the data cannot support
   252	  /// a fit, in which case the caller should keep its previous model rather than use a degenerate one.
   253	  bool fit(const std::vector<std::vector<double>>& X, const std::vector<std::size_t>& pos,
   254	           const std::vector<std::size_t>& neg, const GBTParams& p)
   255	  {
   256	    if (pos.size() < 2 || neg.size() < 2 || X.empty()) { return false; }
   257	    params_ = p;
   258	    n_features_ = X[0].size();
   259	    if (n_features_ == 0) { return false; }
   260	
   261	    std::vector<std::size_t> rows;
   262	    rows.reserve(pos.size() + neg.size());
   263	    rows.insert(rows.end(), pos.begin(), pos.end());
   264	    rows.insert(rows.end(), neg.begin(), neg.end());
   265	    std::vector<double> y(rows.size(), 0.0);
   266	    for (std::size_t i = 0; i < pos.size(); ++i) { y[i] = 1.0; }
   267	
   268	    if (p.fixed_bins)
   269	    {
   270	      std::vector<std::size_t> all_rows(X.size());
   271	      std::iota(all_rows.begin(), all_rows.end(), std::size_t{0});
   272	      edges_ = gbt_detail::computeBinEdges(X, all_rows, p.n_bins);
   273	    }
   274	    else
   275	    {
   276	      edges_ = gbt_detail::computeBinEdges(X, rows, p.n_bins);
   277	    }
   278	    // A feature with no edges cannot split; that is fine, it just never wins a gain comparison.
   279	    //
   280	    // FLAT, not vector<vector>: the histogram pass below is 93% of a fit and reads every feature of
   281	    // every row. One heap block per row would put a pointer chase in front of each of those reads
   282	    // and leave consecutive rows scattered; a flat buffer makes the whole pass a contiguous sweep.
   283	    const std::size_t n_rows = rows.size();
   284	    std::vector<uint8_t> B(n_rows * n_features_, 0);
   285	#ifdef _OPENMP
   286	#pragma omp parallel for schedule(static) num_threads(p.n_threads > 0 ? p.n_threads : omp_get_max_threads())
   287	#endif
   288	    for (long long ii = 0; ii < static_cast<long long>(n_rows); ++ii)
   289	    {
   290	      const std::size_t i = static_cast<std::size_t>(ii);
   291	      for (std::size_t f = 0; f < n_features_; ++f)
   292	      {
   293	        const double v = X[rows[i]][f];
   294	        if (!std::isfinite(v)) { B[i * n_features_ + f] = static_cast<uint8_t>(edges_[f].size()); continue; }
   295	        const auto it = std::lower_bound(edges_[f].begin(), edges_[f].end(), v);
   296	        B[i * n_features_ + f] = static_cast<uint8_t>(it - edges_[f].begin());
   297	      }
   298	    }
   299	
   300	    // Base score = log-odds of the training prior, so the first tree corrects a residual rather
   301	    // than having to discover the class balance.
   302	    const double frac = static_cast<double>(pos.size()) / static_cast<double>(rows.size());
   303	    const double clamped = std::min(1.0 - 1e-6, std::max(1e-6, frac));
   304	    base_ = p.intercept_zero ? 0.0 : std::log(clamped / (1.0 - clamped));
   305	
   306	    std::vector<double> pred(n_rows, base_);
   307	    std::vector<double> grad(n_rows), hess(n_rows);
   308	    trees_.clear();
   309	    trees_.reserve(static_cast<std::size_t>(std::max(0, p.n_trees)));
   310	
   311	    std::size_t n_bins_max = 1;
   312	    for (const auto& e : edges_) { n_bins_max = std::max(n_bins_max, e.size() + 1); }
   313	
   314	    // One scratch buffer for the whole fit. growTree_ clears only the slice each level needs, so
   315	    // the shape is set by the WIDEST level (2^(depth-1) nodes) and every level after the first
   316	    // reuses the same allocation.
   317	    std::vector<gbt_detail::GH> scratch;
   318	
   319	    const int T = p.n_threads > 0 ? p.n_threads :
   320	#ifdef _OPENMP
   321	      omp_get_max_threads();
   322	#else
   323	      1;
   324	#endif
   325	    for (int t = 0; t < p.n_trees; ++t)
   326	    {
   327	      // Element-wise, so no reduction and no ordering question.
   328	#ifdef _OPENMP
   329	#pragma omp parallel for schedule(static) num_threads(T)
   330	#endif
   331	      for (long long ii = 0; ii < static_cast<long long>(n_rows); ++ii)
   332	      {
   333	        const std::size_t i = static_cast<std::size_t>(ii);
   334	        const double pi = 1.0 / (1.0 + std::exp(-pred[i]));
   335	        grad[i] = pi - y[i];
   336	        hess[i] = std::max(1e-12, pi * (1.0 - pi));
   337	      }
   338	      round_lr_ = params_.learning_rate;
   339	      if (p.warmup_rounds > 0)
   340	      { round_lr_ = params_.learning_rate * std::min(1.0, static_cast<double>(t + 1) / static_cast<double>(p.warmup_rounds)); }
   341	      gbt_detail::Tree tree;
   342	      if (!growTree_(B, n_rows, grad, hess, n_bins_max, T, scratch, tree)) { break; }
   343	#ifdef _OPENMP
   344	#pragma omp parallel for schedule(static) num_threads(T)
   345	#endif
   346	      for (long long ii = 0; ii < static_cast<long long>(n_rows); ++ii)
   347	      {
   348	        const std::size_t i = static_cast<std::size_t>(ii);
   349	        pred[i] += tree.predict(&B[i * n_features_]);
   350	      }
   351	      trees_.push_back(std::move(tree));
   352	    }
   353	    return !trees_.empty();
   354	  }
   355	
   356	  /// Raw margin (log-odds). Monotone in the probability, so it can be ranked directly.
   357	  double score(const std::vector<double>& x) const
   358	  {
   359	    if (trees_.empty()) { return 0.0; }
   360	    std::vector<uint8_t> b(n_features_, 0);
   361	    for (std::size_t f = 0; f < n_features_; ++f)
   362	    {
   363	      const double v = x[f];
   364	      if (!std::isfinite(v)) { b[f] = static_cast<uint8_t>(edges_[f].size()); continue; }
   365	      const auto it = std::lower_bound(edges_[f].begin(), edges_[f].end(), v);
   366	      b[f] = static_cast<uint8_t>(it - edges_[f].begin());
   367	    }
   368	    double s = base_;
   369	    for (const auto& t : trees_) { s += t.predict(b.data()); }
   370	    return s;
   371	  }
   372	
   373	  bool trained() const { return !trees_.empty(); }
   374	  std::size_t featureCount() const { return n_features_; }
   375	
   376	  /// Serialise everything `score()` consumes -- the bin edges, the base score and the trees --
   377	  /// as plain text, one token per value, doubles at full precision. The training parameters are
   378	  /// not needed to apply a model and are not written. Returns false on a write error.
   379	  bool save(std::ostream& os) const
   380	  {
   381	    if (trees_.empty()) { return false; }
   382	    os.imbue(std::locale::classic());   // "1.5", never "1,5", whatever the process locale
   383	    os.unsetf(std::ios::floatfield);    // default format: a caller's std::fixed would zero 1e-20
   384	    os.precision(17);
   385	    os << "ODIA-GBT 1\n" << n_features_ << ' ' << base_ << ' ' << trees_.size() << '\n';
   386	    for (const auto& e : edges_)
   387	    {
   388	      os << e.size();
   389	      for (const double v : e) { os << ' ' << v; }
   390	      os << '\n';
   391	    }
   392	    for (const auto& t : trees_)
   393	    {
   394	      os << t.max_depth << ' ' << t.feature.size() << '\n';
   395	      for (std::size_t i = 0; i < t.feature.size(); ++i)
   396	      {
   397	        os << t.feature[i] << ' ' << static_cast<int>(t.bin[i]) << ' ' << t.value[i] << '\n';
   398	      }
   399	    }
   400	    return static_cast<bool>(os);
   401	  }
   402	
   403	  /// Inverse of save(). A model loaded here scores exactly as the one that was saved, provided
   404	  /// the caller feeds it features standardised the way the SAVING run standardised them -- that
   405	  /// transform is the caller's to carry (see scoreSemiSupervisedLDA), not the tree's.
   406	  bool load(std::istream& is)
   407	  {
   408	    // Every count and index is bounded before it sizes a vector or indexes one:
   409	    // a file that parses but describes an impossible tree must fail here, not
   410	    // in predict(). Bounds are generous against anything fit() can produce.
   411	    is.imbue(std::locale::classic());
   412	    std::string magic;
   413	    int version = 0;
   414	    if (!(is >> magic >> version) || magic != "ODIA-GBT" || version != 1) { return false; }
   415	    std::size_t n_trees = 0;
   416	    if (!(is >> n_features_ >> base_ >> n_trees) || n_features_ == 0 || n_features_ > 4096 ||
   417	        n_trees == 0 || n_trees > 100000 || !std::isfinite(base_))
   418	    { return false; }
   419	    edges_.assign(n_features_, {});
   420	    for (auto& e : edges_)
   421	    {
   422	      std::size_t k = 0;
   423	      if (!(is >> k) || k > 255) { return false; }   // bins are uint8_t; <= 255 edges
   424	      e.resize(k);
   425	      for (double& v : e) { if (!(is >> v) || !std::isfinite(v)) { return false; } }
   426	      // score() bins with lower_bound: the edges must be non-decreasing or a
   427	      // value lands in an arbitrary bin.
   428	      for (std::size_t i = 1; i < e.size(); ++i) { if (e[i] < e[i - 1]) { return false; } }
   429	    }
   430	    trees_.clear();
   431	    trees_.reserve(n_trees);
   432	    for (std::size_t t = 0; t < n_trees; ++t)
   433	    {
   434	      gbt_detail::Tree tree;
   435	      std::size_t nodes = 0;
   436	      if (!(is >> tree.max_depth >> nodes)) { return false; }
   437	      // A complete binary tree indexed 1..2^(d+1)-1 is stored in 2^(d+1) slots.
   438	      if (tree.max_depth < 0 || tree.max_depth > 16 ||
   439	          nodes != (static_cast<std::size_t>(1) << (tree.max_depth + 1)))
   440	      { return false; }
   441	      tree.feature.resize(nodes);
   442	      tree.bin.resize(nodes);
   443	      tree.value.resize(nodes);
   444	      for (std::size_t i = 0; i < nodes; ++i)
   445	      {
   446	        int b = 0;
   447	        if (!(is >> tree.feature[i] >> b >> tree.value[i])) { return false; }
   448	        if (tree.feature[i] < -1 || tree.feature[i] >= static_cast<int>(n_features_) ||
   449	            b < 0 || b > 255 || !std::isfinite(tree.value[i]))
   450	        { return false; }
   451	        tree.bin[i] = static_cast<uint8_t>(b);
   452	      }
   453	      trees_.push_back(std::move(tree));
   454	    }
   455	    return !trees_.empty();
   456	  }
   457	  std::size_t nTrees() const { return trees_.size(); }
   458	
   459	private:
   460	  /// Level-wise growth. `node_of_row` tracks which node each row is in; a level is built by
   461	  /// histogramming (grad, hess) per (node, feature, bin) and picking the best split per node.
   462	  bool growTree_(const std::vector<uint8_t>& B, std::size_t n_rows, const std::vector<double>& grad,
   463	                 const std::vector<double>& hess, std::size_t n_bins_max, int n_threads,
   464	                 std::vector<gbt_detail::GH>& scratch, gbt_detail::Tree& tree)
   465	  {
   466	    const int D = std::max(1, params_.max_depth);
   467	    const std::size_t n_nodes = (std::size_t(1) << (D + 1));
   468	    tree.max_depth = D;
   469	    tree.feature.assign(n_nodes, -1);
   470	    tree.bin.assign(n_nodes, 0);
   471	    tree.value.assign(n_nodes, 0.0);
   472	
   473	    std::vector<int> node_of_row(n_rows, 1);
   474	    bool any_split = false;
   475	
   476	    for (int depth = 0; depth < D; ++depth)
   477	    {
   478	      const int first = 1 << depth, last = (1 << (depth + 1)) - 1;
   479	      // histogram[node][feature][bin] -> (G, H)
   480	      const std::size_t n_at_level = static_cast<std::size_t>(last - first + 1);
   481	      const std::size_t hist_size = n_at_level * n_features_ * n_bins_max;
   482	      std::vector<gbt_detail::GH> H(hist_size);
   483	      std::vector<std::size_t> node_rows(n_at_level, 0);
   484	
   485	      // This loop is ~93% of a fit. It is a REDUCTION over rows, so threads need private buffers.
   486	      //
   487	      // The partition is by a fixed number of CHUNKS, not by thread: a per-thread partition would
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
   511	        gbt_detail::GH* __restrict h = scratch.data() + static_cast<std::size_t>(c) * hist_size;
   512	        std::size_t* __restrict nr = node_rows_c.data() + static_cast<std::size_t>(c) * n_at_level;
   513	        for (std::size_t i = lo; i < hi; ++i)
   514	        {
   515	          const int nd = node_of_row[i];
   516	          if (nd < first || nd > last) { continue; }
   517	          const std::size_t slot = static_cast<std::size_t>(nd - first);
   518	          ++nr[slot];
   519	          const std::size_t base = (slot * n_features_) * n_bins_max;
   520	          const uint8_t* __restrict brow = B.data() + i * n_features_;
   521	          const double gi = grad[i], hi_ = hess[i];
   522	          for (std::size_t f = 0; f < n_features_; ++f)
   523	          {
   524	            const std::size_t k = base + f * n_bins_max + brow[f];
   525	            h[k].g += gi;
   526	            h[k].h += hi_;
   527	            h[k].n += 1;
   528	          }
   529	        }
   530	      }
   531	      // Fixed-order combination -- this is what keeps the fit thread-count independent.
   532	      for (int c = 0; c < nchunk; ++c)
   533	      {
   534	        const gbt_detail::GH* hc = scratch.data() + static_cast<std::size_t>(c) * hist_size;
   535	        const std::size_t* nr = node_rows_c.data() + static_cast<std::size_t>(c) * n_at_level;
   536	        for (std::size_t k = 0; k < hist_size; ++k)
   537	        { H[k].g += hc[k].g; H[k].h += hc[k].h; H[k].n += hc[k].n; }
   538	        for (std::size_t k = 0; k < n_at_level; ++k) { node_rows[k] += nr[k]; }
   539	      }
   540	
   541	      for (int nd = first; nd <= last; ++nd)
   542	      {
   543	        const std::size_t slot = static_cast<std::size_t>(nd - first);
   544	        if (node_rows[slot] == 0) { continue; }
   545	        double G = 0.0, Hs = 0.0;
   546	        const std::size_t base = (slot * n_features_) * n_bins_max;
   547	        for (std::size_t b = 0; b < n_bins_max; ++b)
   548	        {
   549	          G += H[base + b].g;
   550	          Hs += H[base + b].h;
   551	        }
   552	        const double parent_obj = nodeObj_(G, Hs);
   553	
   554	        int best_f = -1;
   555	        std::size_t best_b = 0;
   556	        double best_gain = params_.gamma;
   557	        for (std::size_t f = 0; f < n_features_; ++f)
   558	        {
   559	          double GL = 0.0, HL = 0.0;
   560	          std::size_t nL = 0;
   561	          const std::size_t fb = base + f * n_bins_max;
   562	          // Scan cut points left to right; the last bin is never a cut (it would put everything left).
   563	          for (std::size_t b = 0; b + 1 < n_bins_max; ++b)
   564	          {
   565	            GL += H[fb + b].g;
   566	            HL += H[fb + b].h;
   567	            nL += H[fb + b].n;
   568	            const double GR = G - GL, HR = Hs - HL;
   569	            if (HL < params_.min_child_weight || HR < params_.min_child_weight) { continue; }
   570	            // PER-CHILD row floor, now actually enforced. The parent check below
   571	            // (node_rows >= 2*min_child_rows) permits a 39/1 split; this forbids
   572	            // it. A leaf of a handful of rows is fitted noise, and it lands in
   573	            // the score tail that sets the FDR threshold.
   574	            const std::size_t nR = node_rows[slot] - nL;
   575	            if (nL < static_cast<std::size_t>(params_.min_child_rows) ||
   576	                nR < static_cast<std::size_t>(params_.min_child_rows)) { continue; }
   577	            const double gain = 0.5 * (nodeObj_(GL, HL) + nodeObj_(GR, HR) - parent_obj);
   578	            // Strictly greater: ties keep the lowest (feature, bin), which is what makes the fit
   579	            // reproducible independent of iteration order.
   580	            if (gain > best_gain)
   581	            {
   582	              best_gain = gain;
   583	              best_f = static_cast<int>(f);
   584	              best_b = b;
   585	            }
   586	          }
   587	        }
   588	
   589	        const bool leaf_level = (depth == D - 1);
   590	        if (best_f < 0 || leaf_level || node_rows[slot] < static_cast<std::size_t>(2 * params_.min_child_rows))
   591	        {
   592	          tree.feature[nd] = -1;
   593	          tree.value[nd] = -round_lr_ * leafStep_(G, Hs);
   594	          continue;
   595	        }
   596	        tree.feature[nd] = best_f;
   597	        tree.bin[nd] = static_cast<uint8_t>(best_b);
   598	        any_split = true;
   599	        // Children inherit a provisional leaf value in case the next level does not split them.
   600	        for (int child : {2 * nd, 2 * nd + 1})
   601	        {
   602	          tree.feature[child] = -1;
   603	          tree.value[child] = 0.0;
   604	        }
   605	      }
   606	
   607	      // Route rows into the next level and give the leaves their values.
   608	      for (std::size_t i = 0; i < n_rows; ++i)
   609	      {
   610	        const int nd = node_of_row[i];
   611	        if (nd < first || nd > last || tree.feature[nd] < 0) { continue; }
   612	        node_of_row[i] = (B[i * n_features_ + tree.feature[nd]] <= tree.bin[nd]) ? 2 * nd : 2 * nd + 1;
   613	      }
   614	    }
   615	
   616	    // Final level: every row now sits in a node whose value must be set from its own G/H.
   617	    const int first = 1 << D, last = (1 << (D + 1)) - 1;
   618	    std::vector<double> G(static_cast<std::size_t>(last - first + 1), 0.0), H(G.size(), 0.0);
   619	    for (std::size_t i = 0; i < n_rows; ++i)
   620	    {
   621	      const int nd = node_of_row[i];
   622	      if (nd < first || nd > last) { continue; }
   623	      G[static_cast<std::size_t>(nd - first)] += grad[i];
   624	      H[static_cast<std::size_t>(nd - first)] += hess[i];
   625	    }
   626	    for (int nd = first; nd <= last; ++nd)
   627	    {
   628	      const std::size_t s = static_cast<std::size_t>(nd - first);
   629	      if (H[s] <= 0.0) { continue; }
   630	      tree.feature[nd] = -1;
   631	      tree.value[nd] = -round_lr_ * leafStep_(G[s], H[s]);
   632	    }
   633	    return any_split || true;   // a depth-0 stump is still a valid (if useless) tree
   634	  }
   635	
   636	  GBTParams params_;
   637	  double round_lr_ = 0.0;         ///< the learning rate in force for the tree being grown (== params_.learning_rate unless warmup_rounds > 0)
   638	  std::vector<gbt_detail::Tree> trees_;
   639	  std::vector<std::vector<double>> edges_;
   640	  std::size_t n_features_ = 0;
   641	  double base_ = 0.0;
   642	  /// A node's objective reduction for the split search: G^2/(H+lambda) natively; with clip_gain and a cap,
   643	  /// the reduction the CLIPPED weight realises, 2*G*step - (H+lambda)*step^2 (equal when unclipped).
   644	  double nodeObj_(double G, double H) const
   645	  {
   646	    if (params_.clip_gain && params_.max_delta_step > 0.0)
   647	    {
   648	      const double step = leafStep_(G, H);
   649	      return 2.0 * G * step - (H + params_.lambda) * step * step;
   650	    }
   651	    return G * G / (H + params_.lambda);
   652	  }
   653	  /// The leaf's Newton step before the learning rate, clipped to +/-max_delta_step when that is set.
   654	  double leafStep_(double G, double H) const
   655	  {
   656	    double step = G / (H + params_.lambda);
   657	    if (params_.max_delta_step > 0.0)
   658	    { step = std::max(-params_.max_delta_step, std::min(params_.max_delta_step, step)); }
   659	    return step;
   660	  }
   661	};
   662	
   663	} // namespace ODIA::Scoring
   664	
   665	#endif // ODIA_GBT_H
```
### include/odia/PeakGroupScorer.h lines 30-430 (sub-score enum, fragvec declaration) and 860-1120 (Options, groups, Result)
```
    30	  /// `mass_error` needs the observed m/z of each matched peak, which the
    31	  /// extractor discards, `im_delta` needs the CCS -> 1/K0 conversion that is
    32	  /// deliberately downstream, and `isotope_corr` needs MS1 extraction that
    33	  /// does not exist. A sub-score computed from a placeholder is worse than an
    34	  /// absent one: the classifier gives it weight and it is noise.
    35	  class PeakGroupScorer
    36	  {
    37	  public:
    38	    /// Order matters -- it is the column order of the feature matrix, and it is
    39	    /// what `subScoreNames()` reports. Adding one means retraining.
    40	    enum SubScore : std::size_t
    41	    {
    42	      XCORR_SHAPE = 0,     ///< mean all-pairs cross-correlation at the peak
    43	      XCORR_COELUTION,     ///< mean |lag| of those maxima; coelution, so lower is better
    44	      LIBRARY_CORR,        ///< Pearson, observed against library intensities
    45	      LIBRARY_DOTPROD,     ///< normalised dot product, same pair
    46	      /// Replaced the old group/window area ratio (D6). That ratio carried no
    47	      /// library or co-elution information at all: a narrow decoy spike in an
    48	      /// otherwise empty window approaches 1.0, while a real target peak on a
    49	      /// real baseline scores lower -- which is why targets measured WORSE on
    50	      /// it (0.1582 against 0.1613). This is the fraction of the group's area
    51	      /// contributed by the fragments the library says should be brightest,
    52	      /// which a single-transition spike cannot satisfy.
    53	      INTENSITY_SCORE,
    54	      LOG_SN,              ///< log apex over a data-derived background floor
    55	      /// How many fragments actually carried information. A precursor scored
    56	      /// from three live fragments is not the same evidence as one scored from
    57	      /// twelve, and without this the two are indistinguishable to the
    58	      /// classifier.
    59	      /// Fragments whose own maximum coincides with the group apex (within one
    60	      /// cycle). Co-elution depth, not signal presence: the previous definition
    61	      /// counted non-degenerate traces and returned a constant 12.000 for right
    62	      /// and wrong answers alike.
    63	      USABLE_FRAGMENTS,
    64	
    65	      // Added toward the 10-14 orthogonal scores every working implementation
    66	      // in this family uses. Rosenberger 2017 used 14, DIA-NN uses 73; seven
    67	      // extracting a 58/42 edge from data holding an 8x separation is a scorer
    68	      // roughly half the width it needs, not one missing dominant feature.
    69	      // Names follow pyProphet's, so the columns are comparable to a published
    70	      // weight vector rather than privately invented.
    71	
    72	      /// Root-mean-square deviation of normalised observed against normalised
    73	      /// library intensities. Carries what correlation discards: dot product
    74	      /// and Pearson are both scale-free, so a spectrum with the right SHAPE
    75	      /// but the wrong contrast scores identically to a correct one.
    76	      LIBRARY_RMSD,
    77	      /// Share of the group's corrected area in y-ions. Tryptic y-ions carry
    78	      /// the basic C-terminal residue and dominate a real spectrum; a decoy's
    79	      /// mutated sequence has no reason to preserve the ratio.
    80	      YSERIES_SCORE,
    81	      /// Fraction of the fragments the library lists that were seen at all.
    82	      /// Distinct from USABLE_FRAGMENTS, which counts traces carrying any
    83	      /// signal: this counts those carrying signal ABOVE their own background
    84	      /// inside the candidate, which is the co-elution the group claims.
    85	      FRAGMENT_COVERAGE,
    86	
    87	      // --- added 2026-08-07. The first four are free: every quantity was
    88	      // already computed and thrown away. See doc/13.
    89	
    90	      /// The summed pairwise correlation the co-elution detector used to accept
    91	      /// this position as a peak at all.
    92	      ///
    93	      /// This is the quantity that moved rank-1 accuracy from 40.7% to 75.7%,
    94	      /// and until now the picker computed it per candidate and discarded it.
    95	      /// DIA-NN keeps both `best_corr_sum` and `total_corr_sum` as normalising
    96	      /// context for exactly this reason. Zero for the amplitude picker, which
    97	      /// never computes it.
    98	      CORR_SUM,
    99	
   100	      /// How far ahead of its own runner-up this candidate was, on CORR_SUM.
   101	      ///
   102	      /// Distinguishes "this precursor had one obvious answer" from "three
   103	      /// equally plausible ones", which no per-candidate score can express.
   104	      /// DIA-NN encodes the same idea at detection time as MaxCorrDiff, keeping
   105	      /// candidates by margin rather than rank; as a feature it reaches the
   106	      /// classifier instead of only the picker.
   107	      CANDIDATE_MARGIN,
   108	
   109	      /// Width of the candidate in cycles, over the run's median candidate
   110	      /// width. A peptide elutes on the chromatography's timescale;
   111	      /// interference need not.
   112	      PEAK_WIDTH_RATIO,
   113	
   114	      // RT_DELTA -- |apex RT - predicted RT| -- was here and is REMOVED.
   115	      //
   116	      // A peak group has ONE retention time and its traces co-elute by
   117	      // construction, so a group-level retention-time delta cannot separate a
   118	      // real group from an interference group: both sit wherever the signal
   119	      // they were built from sits. Retention time is diagnostic BETWEEN traces,
   120	      // not for the group as a whole.
   121	      //
   122	      // Measured before removal (d9_auc_by_abundance.py, full IH1): AUC 0.215
   123	      // against the decoy null, flat across every abundance quintile, and it
   124	      // got WORSE under random decoy rows (0.293 -> 0.215), so it was not an
   125	      // artefact of picking decoys by argmax. Target and decoy row
   126	      // distributions were identical -- median 41.7 s against 42.6 s, in a
   127	      // window of half-width ~75 s -- i.e. a true peak sits no closer to the
   128	      // predicted time than a random candidate does.
   129	      //
   130	      // It was also pinned in `nonpositive_features` -- but that clip applies
   131	      // to the LDA weight vector only and the default classifier is GBT, so
   132	      // the pin never took effect on the shipped path. Removing the feature
   133	      // measured as a WASH on the fixture (+20 IDs, FDP +1.03 pp, both inside
   134	      // noise), which is what an uninformative feature a tree ensemble already
   135	      // ignores should do. The removal is right because the quantity cannot
   136	      // discriminate, not because it was costing identifications.
   137	
   138	      /// |library 1/K0 - observed 1/K0| for the precursor, or NaN where the run
   139	      /// or the library has no mobility. Orthogonal to everything above on
   140	      /// diaPASEF and simply absent elsewhere.
   141	      IM_DELTA,
   142	
   143	      /// Pearson of the MS1 precursor trace against the MS2 fragment consensus
   144	      /// over the candidate's own cycles.
   145	      ///
   146	      /// The first sub-score here that does not read MS2 fragment traces. That
   147	      /// is the whole point: the other fifteen share a failure mode, because a
   148	      /// co-eluting interferent corrupts all of them at once. This asks whether
   149	      /// the PRECURSOR rises and falls with its fragments, which no amount of
   150	      /// fragment-side interference can fake.
   151	      ///
   152	      /// Measured on IH1 + v6_50k before being added (doc/13): 13.7x enrichment
   153	      /// in the top bin, 1.4-1.5x in bulk, and it survives stratification by MS1
   154	      /// intensity so it is shape rather than brightness. NaN when the run
   155	      /// carries no MS1 or the precursor has no MS1 signal -- NaN, not zero,
   156	      /// because zero is a legitimate correlation and the two must not be
   157	      /// confused.
   158	      MS1_COELUTION,
   159	
   160	      // --- added 2026-08-09. Both reviewers, independently, named fragment
   161	      // mass accuracy as the highest-value score ODIA does not have, and both
   162	      // pointed at the same reason: an extracted chromatogram records intensity
   163	      // inside a mass window over time and throws away WHERE in that window the
   164	      // peak sat. An interferent can co-elute perfectly, correlate perfectly
   165	      // and be systematically displaced in m/z, and no chromatogram-shape
   166	      // statistic can see that. It is also the cheapest orthogonal channel
   167	      // available: both quantities are already computed and stored on the
   168	      // group, they were simply never offered to the classifier.
   169	
   170	      /// Negated |deviation - the run's own centre|, ppm. Larger is better.
   171	      ///
   172	      /// CENTRED ON THE RUN, which is what makes it safe. The objection that
   173	      /// kept this out of the feature vector was real: fed raw, the classifier
   174	      /// would learn "this run's fragments sit at -3 ppm" and penalise
   175	      /// correctly calibrated identifications. Centring on the run's own median
   176	      /// removes exactly that, and what is left is the per-group departure from
   177	      /// the instrument's systematic error -- a property of the GROUP.
   178	      MASS_ACCURACY,
   179	
   180	      /// Negated scatter of the group's per-fragment deviations, ppm. Larger is
   181	      /// better.
   182	      ///
   183	      /// The stronger of the two, and calibration-free by construction: it is a
   184	      /// spread WITHIN one group, so no centring can affect it and no run-level
   185	      /// offset can leak in. A real peptide's fragments are one molecule
   186	      /// measured through one calibration and sit together; interference is
   187	      /// unrelated species whose errors scatter.
   188	      MASS_SPREAD,
   189	
   190	      // --- added 2026-08-10. The mobility analogue of MASS_SPREAD, and the
   191	      // first sub-score that can separate interference m/z and RT cannot.
   192	
   193	      /// Scatter of this group's PER-FRAGMENT observed 1/K0, as a robust sigma.
   194	      ///
   195	      /// A precursor's fragments are produced from ONE ion packet and therefore
   196	      /// share its mobility. A co-eluting interferent is a different packet at a
   197	      /// different 1/K0, so a group whose fragments disagree about mobility is a
   198	      /// mixture -- and on diaPASEF that is detectable even when the species are
   199	      /// inseparable in m/z and in retention time, because a frame merges TIMS
   200	      /// scans that the raw data still resolves.
   201	      ///
   202	      /// Distinct from IM_DELTA, which asks whether the group sits where the
   203	      /// LIBRARY predicted. This asks whether the group is internally
   204	      /// consistent, which needs no prediction and so cannot be wrong for the
   205	      /// library's reasons.
   206	      ///
   207	      /// NaN on a run without ion mobility, and until at least
   208	      /// `min_im_fragments` fragments carried one -- NaN rather than zero,
   209	      /// because zero scatter is what a perfect group looks like.
   210	      IM_SPREAD,
   211	
   212	      // --- added 2026-08-23. The RETENTION-TIME analogue of MASS_SPREAD and
   213	      // IM_SPREAD, and the only form in which retention time can discriminate
   214	      // at all.
   215	      //
   216	      // RT_DELTA asked whether the GROUP sits where the library predicted, and
   217	      // was removed: a group has one retention time and its traces co-elute by
   218	      // construction, so both a real group and an interference group sit
   219	      // wherever the signal they were built from sits. There is nothing there
   220	      // to separate.
   221	      //
   222	      // BETWEEN fragments there is. A peptide's fragments are produced from one
   223	      // ion packet and share one elution profile, so their individual apices
   224	      // agree. An interfering fragment belongs to a different species and peaks
   225	      // somewhere else. That is a property of the group's internal consistency
   226	      // and needs no prediction, which is what makes it immune to the library
   227	      // being wrong -- the same argument that makes IM_SPREAD work.
   228	
   229	      /// Weighted scatter of the fragments' own retention-time centroids about
   230	      /// THEIR OWN WEIGHTED MEAN, in seconds, NEGATED so that higher is better
   231	      /// like every other sub-score.
   232	      ///
   233	      /// Dispersion BETWEEN fragments, not displacement from the group apex:
   234	      /// that apex comes from the summed trace, so one loud interferent would
   235	      /// drag the reference and the measurement together.
   236	      ///
   237	      /// NaN for candidates narrower than 3 cycles -- with one or two, every
   238	      /// centroid is forced to the same value and the scatter is 0, the BEST
   239	      /// possible score, so a three-fragment noise spike would look perfect.
   240	      ///
   241	      /// The centroid rather than the argmax: on a faint fragment the argmax
   242	      /// jumps between adjacent cycles on noise, and the whole point of the
   243	      /// feature is to work where the shape features stop working (measured:
   244	      /// at the faintest abundance quintile var_library_corr is 0.532 and
   245	      /// var_xcorr_shape 0.598, against 0.913 and 0.962 at the brightest).
   246	      ///
   247	      /// Weighted by each fragment's background-corrected area, so a fragment
   248	      /// that is mostly noise does not get an equal vote on where the peak is.
   249	      ///
   250	      /// NaN with fewer than `min_rt_spread_fragments` informative fragments --
   251	      /// NaN rather than zero, because zero scatter is what a perfect group
   252	      /// looks like and a group with one fragment would score perfectly.
   253	      ///
   254	      /// Related to XCORR_COELUTION but not the same statistic: that is a
   255	      /// cross-correlation LAG between fragment PAIRS, this is each fragment's
   256	      /// offset from the group apex. Whether it is redundant with it is a
   257	      /// question for measurement, not for this comment.
   258	      RT_SPREAD,
   259	
   260	      /// What FRACTION of a fragment's matched intensity survives tightening
   261	      /// the mass tolerance, averaged over fragments.
   262	      ///
   263	      /// The idea is DIA-NN's -- a real peak sits on its calibrated m/z and
   264	      /// survives tightening, while an interferent that merely fell inside the
   265	      /// window often does not -- but the provenance claimed here first was
   266	      /// wrong and is corrected. DIA-NN's 1x / 0.45x / 0.20x tolerances are its
   267	      /// NINE MS1 channels, feeding tight-tolerance CORRELATIONS
   268	      /// (pMs1TightOne/Two); its per-fragment mass feature is pAcc[6], a
   269	      /// deviation rather than a surviving-intensity fraction. So 4.5 and 2.0
   270	      /// ppm are an MS1 constant imported into MS2, and they are not what
   271	      /// DIA-NN does with fragments.
   272	      ///
   273	      /// Measured on real IH1 frames rather than assumed: a cell is a MIXTURE,
   274	      /// and 95.1% of bright real-peak cells contain more than one peak. The
   275	      /// weighted mean therefore pulls crowded cells toward the window centre,
   276	      /// which cuts targets from a nominal 1.0 to 0.925 at 4.5 ppm AND lifts
   277	      /// decoy-like cells from 0.45 to 0.564 -- a nominal gap of 0.55 measured
   278	      /// at 0.36. At 2.0 ppm, 43% of perfectly on-mass bright cells fail purely
   279	      /// from co-window neighbours, so that setting measures local crowding
   280	      /// rather than mass correctness. **4.5 is the value to test; 2.0 is not.**
   281	      ///
   282	      /// No re-extraction was needed. The extractor already stores
   283	      /// sum(intensity * ppm) and sum(intensity) per cell, so the tightened
   284	      /// query is a filter over cells that are already in memory rather than a
   285	      /// second pass over the raw data.
   286	      ///
   287	      /// LABEL-SYMMETRIC in the way the two changes reverted before it were
   288	      /// not. It consults no library intensity, which is what broke those: a
   289	      /// decoy copies its target's per-fragment intensities verbatim while its
   290	      /// fragment m/z ARE recomputed, so any intensity-weighted statistic
   291	      /// silently asks a different question of each class. This asks both the
   292	      /// same question -- how far is the matched signal from the m/z THIS
   293	      /// precursor's fragment should have -- and a decoy's recomputed m/z makes
   294	      /// that a genuine test rather than a scrambled one.
   295	      ///
   296	      /// Orthogonal by construction to everything already here, which the
   297	      /// project's own rule says is the lever rather than count: the existing
   298	      /// mass features are a deviation and a scatter, both summaries of WHERE
   299	      /// the matched peaks sit. This is how much intensity is still there when
   300	      /// the window closes, which a median deviation cannot express -- one
   301	      /// fragment can have a perfect median and lose most of its area.
   302	      ///
   303	      /// NaN when the ppm planes are absent (`-collect_mass_residuals` off) or
   304	      /// no fragment matched a peak. Note the imputation hazard that bit an
   305	      /// earlier feature: PercolatorEngine fills NaN with the column median, so
   306	      /// a candidate that could not be measured is scored as an average one.
   307	      /// Here the NaN case is a whole-run structural absence rather than a
   308	      /// per-candidate weakness, which is the case that imputation was designed
   309	      /// for.
   310	      MASS_SURVIVAL,
   311	
   312	      /// A DELIBERATELY UNINFORMATIVE COLUMN. Not a feature -- a control.
   313	      ///
   314	      /// Six unrelated changes have now been measured on the fixture at matched
   315	      /// entrapment FDP, and every one produced the same profile: slightly
   316	      /// negative at 5.72-10% and +6 to +8% at 15%. Six coincidences is not a
   317	      /// hypothesis. The alternative is that ADDING A COLUMN AT ALL perturbs the
   318	      /// semi-supervised classifier's trajectory, and that at the operating
   319	      /// point the perturbation is larger than anything a single feature's
   320	      /// information content contributes.
   321	      ///
   322	      /// This column decides between those. It is a deterministic hash of the
   323	      /// precursor index and apex cycle, scaled to [0,1): it varies per
   324	      /// candidate, so it is not constant and survives the constant-column
   325	      /// guard, and it is uniform with respect to label, so it carries no
   326	      /// information a classifier could legitimately use.
   327	      ///
   328	      /// If an arm carrying it reproduces that profile, then the fixture cannot
   329	      /// resolve feature-level changes at this effect size and six "failures to
   330	      /// convert" collapse into one measurement artefact. If it comes back flat,
   331	      /// the six results stand and the features really were not worth their
   332	      /// slots.
   333	      ///
   334	      /// Off unless `-null_feature` is given. It must never be on in a real run.
   335	      ///
   336	      /// IT REPRODUCED THE PROFILE. Pure noise gives -1.5% at 7.42%, -3.9% at
   337	      /// 10% and +5.7% at 15%, inside the range the four real changes gave
   338	      /// (-0.9 to -3.1%, +6.6 to +8.6%). So the fixture cannot resolve a
   339	      /// feature-sized change: adding any column moves the classifier's
   340	      /// trajectory further than a feature's information content does.
   341	      NULL_CONTROL,
   342	
   343	      // --- added 2026-08-31, the P1 engine port (doc/80 retrospective ->
   344	      // RETRO_VERDICT -> a79_multicand P0 factorial). The cohort-priced recipe
   345	      // -- per-fragment evidence against an interference-robust best-fragment
   346	      // reference, plus candidate-competition context -- measured +3.9 buried
   347	      // / +1.8 overall cohort-projected on top of the shipped 22 (88.1/75.4,
   348	      // best of program). DIA-NN 1.7.12 carries the same three mechanics
   349	      // (per-fragment vectors vs a data-driven reference; margin-kept
   350	      // candidates; competition context). All computed identically for
   351	      // targets and decoys; appended at the enum tail so every existing
   352	      // index is unchanged.
   353	
   354	      /// Sum over the six area-top fragments of the Pearson correlation of
   355	      /// each fragment's trace against the SMOOTHED (1-2-1) trace of the
   356	      /// reference fragment -- the area-top-6 member maximising its summed
   357	      /// pairwise correlation to the other five. Unlike CORR_SUM (all-pairs,
   358	      /// picker-derived), this is correlation TO ONE robust reference, so a
   359	      /// single contaminated fragment lowers its own term without dragging
   360	      /// the pairs of every clean fragment.
   361	      REF_CORR_SUM,
   362	      /// The per-fragment correlations to that reference, SORTED descending
   363	      /// and padded with zeros -- twelve separate columns, so the classifier
   364	      /// sees WHICH fragments agree instead of one aggregate a contaminant
   365	      /// can drag (DIA-NN's pCorr[] design).
   366	      REF_CORR_1, REF_CORR_2, REF_CORR_3, REF_CORR_4, REF_CORR_5, REF_CORR_6,
   367	      REF_CORR_7, REF_CORR_8, REF_CORR_9, REF_CORR_10, REF_CORR_11, REF_CORR_12,
   368	      /// Each fragment's share of the group's corrected area, sorted
   369	      /// descending, top six (DIA-NN's pSig[] design): the spectrum's
   370	      /// concentration profile, which a one-fragment spike cannot fake.
   371	      SIG_SHARE_1, SIG_SHARE_2, SIG_SHARE_3, SIG_SHARE_4, SIG_SHARE_5,
   372	      SIG_SHARE_6,
   373	      /// This candidate's CORR_SUM rank among its precursor's candidates
   374	      /// (1 = best) and how many candidates the precursor produced. With
   375	      /// CANDIDATE_MARGIN these complete the competition context: "the only
   376	      /// peak found" and "best of seven near-ties" are different evidence.
   377	      CAND_RANK,
   378	      CAND_COUNT,
   379	
   380	      N_SUB_SCORES
   381	    };
   382	
   383	    static const std::vector<std::string>& subScoreNames();
   384	
   385	    /// Rung (i) of the fragment-evidence contract: 78 per-candidate columns,
   386	    /// six 12-long vectors by LIBRARY-INTENSITY rank followed by six scalars.
   387	    ///
   388	    /// Written only by `-out_fragvec`, never fitted and never a sub-score. The
   389	    /// names and their order ARE the sealed definition in
   390	    /// analysis77/pick/wf_v33_fragvec_contract.md s.5.1, and the arithmetic is
   391	    /// the reference builder's (wf_v33_fragvec_features.py, `row_features`)
   392	    /// re-expressed on the locals the scorer already has. Kept out of
   393	    /// `subScoreNames()` on purpose: a column the discriminant can see would
   394	    /// change the scores, and this export exists to be provably score-neutral.
   395	    static const std::vector<std::string>& fragvecNames();
   396	
   397	    /// 6 x 12 + 6. A compile-time constant so the writer, the row stride and
   398	    /// the name table cannot drift apart.
   399	    static constexpr std::size_t N_FRAGVEC = 78;
   400	
   401	    /// Test seam for the boundary rule. The rule lives in an anonymous
   402	    /// namespace in the .cpp, which is right for it and leaves no way to assert
   403	    /// candidate GEOMETRY -- and geometry was where the defect was: 77.4% of
   404	    /// production peak groups covered more than 80% of the extraction window.
   405	    static std::pair<std::size_t, std::size_t> peakBoundsForTest(
   406	      const std::vector<double>& smoothed, std::size_t left_from,
   407	      std::size_t right_from, double boundary_fraction,
   408	      std::size_t min_cycles, std::size_t max_half, double sigmas = 1.0);
   409	
   410	    /// Why a library precursor produced no scored candidate.
   411	    ///
   412	    /// Every precursor gets EXACTLY ONE of these, which is the whole point: the
   413	    /// aggregate reject counters cannot be cross-tabulated against a list of
   414	    /// precursors, so "41% of DIA-NN's confident set yields no candidate" could
   415	    /// be attributed to a stage only by inference -- and that inference had
   416	    /// holes. `Session::add` returns silently on `tc == 0`, and a precursor
   417	    /// covered by no isolation window never reaches `add` at all, so neither
   418	    /// appears in any counter.
   419	    enum class TerminalReason : std::uint8_t
   420	    {
   421	      NotReached = 0,     ///< never handed to the scorer: no isolation window, or filtered upstream
   422	      NoTransitions = 1,  ///< reached it with zero transitions extracted
   423	      FewPoints = 2,      ///< fewer than 3 points in the chromatogram
   424	      GateC = 3,          ///< co-elution evidence below the decoy-null quantile
   425	      FewExcursions = 4,  ///< too few transitions rose above their own noise
   426	      ZeroTrace = 5,      ///< the summed trace really is zero
   427	      NoCandidate = 6,    ///< entered the picker and it returned nothing
   428	      Scored = 7,         ///< produced at least one scored candidate
   429	      /// The picker returned candidates and the scorer then discarded every
   430	      /// one of them -- `min_fragments_at_apex` is the only rule that does
   ...
   860	      /// where running either alternative alone left `var_corr_sum` and
   861	      /// `var_candidate_margin` constant and therefore dropped.
   862	      bool union_picking = false;
   863	
   864	      /// Signal-to-noise threshold for the OpenSWATH picker. Its default is 1.0;
   865	      /// OpenSwathWorkflow commonly runs it at 0.1 for DIA, where a "peak" sits
   866	      /// on far more background than in targeted MRM.
   867	      double openswath_sn = 1.0;
   868	
   869	      /// Gaussian smoothing rather than Savitzky-Golay in the OpenSWATH picker.
   870	      bool openswath_gauss = false;
   871	
   872	      /// Expected peak width in seconds for the OpenSWATH picker, or 0 to leave
   873	      /// its default. IH1's peaks are ~20-30 s.
   874	      double openswath_peak_width = 0.0;
   875	
   876	      /// Save the trained discriminant here, or load a frozen one from here.
   877	      ///
   878	      /// Static modelling: at a 1.5% true-positive rate the semi-supervised loop
   879	      /// has no confident positives to bootstrap from and certifies nothing,
   880	      /// while the same discriminant RANKS 232 of DIA-NN's 738 into its top 738.
   881	      /// Training where the loop ignites and applying the frozen model where it
   882	      /// does not is the documented remedy.
   883	      std::string classifier_model_out;
   884	      std::string classifier_model_in;
   885	
   886	      /// Compute and retain the 78 rung-(i) fragment columns per candidate.
   887	      ///
   888	      /// OUTPUT-ONLY: nothing reads `Result::fragvec` back into a score, a
   889	      /// sub-score, a calibration or a candidate decision, so a run with this
   890	      /// on must produce a byte-identical `-out`. Off on pass 1 and the RT
   891	      /// refinement -- they are the same arithmetic and would only pay for it.
   892	      bool fragvec = false;
   893	
   894	      /// MS1 traces for MS1_COELUTION, or null when the run has no MS1.
   895	      ///
   896	      /// A real pointer to real data, not a plumbing point. doc/07 is explicit
   897	      /// that a sub-score computed from a placeholder is worse than an absent
   898	      /// one, and IM_DELTA was all-NaN for a week because its hook was added and
   899	      /// never connected. When this is null MS1_COELUTION is NaN for every row
   900	      /// and the constant-column guard drops it, which is the honest behaviour.
   901	      const Ms1Traces* ms1 = nullptr;
   902	
   903	      /// How much say each fragment gets in RT_SPREAD's centroid scatter.
   904	      ///
   905	      /// Area weighting estimates the dominant ion packet but suppresses the
   906	      /// signal this feature exists to find: for two clusters the weighted
   907	      /// variance scales as W1*W2/(W1+W2)^2, so one weak interfering fragment
   908	      /// contributes almost nothing. None gives every informative fragment an
   909	      /// equal vote, including barely-detected noisy ones. Sqrt is between.
   910	      RtSpreadWeight rt_spread_weight = RtSpreadWeight::Area;
   911	
   912	      /// Informative fragments required before RT_SPREAD is computed at all.
   913	      ///
   914	      /// Three is the smallest number for which a scatter means anything: with
   915	      /// two, the weighted sigma is a rescaled |difference| and every group with
   916	      /// two agreeing fragments looks perfect. The same reasoning that put a
   917	      /// four-point floor under the library correlation.
   918	      std::size_t min_rt_spread_fragments = 3;
   919	
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
   931	
   932	      /// DIAGNOSTIC ORACLE. Per-precursor retention time, in run seconds, at
   933	      /// which a candidate is FORCED to exist: the admission gates are bypassed
   934	      /// and, if the picker finds nothing within `oracle_rt_tol` of it, a
   935	      /// candidate is synthesised there. NaN means no oracle for that precursor.
   936	      ///
   937	      /// This exists to answer one question and must never be a default: if
   938	      /// admission were perfect, how many precursors would we actually
   939	      /// identify? `d5_yield.py` estimates at most +40.9% by extrapolating
   940	      /// acceptance rates across abundance bins; this measures it instead.
   941	      ///
   942	      /// It uses the answer as input, so anything it produces is an UPPER
   943	      /// BOUND on a real method and its q-values are not meaningful: only
   944	      /// targets get injections, decoys cannot (they have no true retention
   945	      /// time), so the null is not comparable. Read the injected candidates'
   946	      /// DScore against a threshold from an UNORACLED run.
   947	      const float* oracle_rt = nullptr;
   948	      double oracle_rt_tol = 20.0;
   949	
   950	      /// Half-window, in cycles, for the pairwise correlation at each position.
   951	      std::size_t corr_half_window = 4;
   952	
   953	      /// The reference fragment's summed correlation to the others must reach
   954	      /// this for a position to be a peak at all. DIA-NN's MinCorrScore.
   955	      double min_corr_score = 0.5;
   956	
   957	      /// Candidates are kept by MARGIN from the best correlation sum rather
   958	      /// than by rank, so an unambiguous precursor yields one candidate and an
   959	      /// ambiguous one yields several. DIA-NN's MaxCorrDiff. This is why a
   960	      /// fixed top-N hurt: at 25 it manufactured 24 competitors regardless of
   961	      /// whether any was plausible, and precision fell to 21.4%.
   962	      double max_corr_diff = 2.0;
   963	
   964	      /// The apex must be this fraction of the maximum evidence nearby on the
   965	      /// reference fragment's own smoothed trace. DIA-NN's PeakApexEvidence.
   966	      double apex_evidence = 0.99;
   967	
   968	      /// Reject a candidate whose observed spectrum correlates with the
   969	      /// library below this. -1.0 disables it, which is the default.
   970	      ///
   971	      /// This is the one sub-score measured to separate real identifications
   972	      /// from misplaced ones: on IH1, median 0.582 for calls landing within
   973	      /// 30 s of the true apex against -0.036 for those that do not -- and
   974	      /// -0.032 for decoys, i.e. a misplaced target is indistinguishable from a
   975	      /// decoy here. See the comment at the gate in PeakGroupScorer.cpp.
   976	      double min_library_corr = -1.0;
   977	
   978	      /// Classifier for the semi-supervised loop.
   979	      /// "lda" is deterministic and dependency-light; "gbt" is what the
   980	      /// upstream benchmarks use.
   981	      std::string classifier = "gbt";
   982	
   983	      unsigned threads = 0;
   984	
   985	      /// True when `Library::precursors().irt` holds RUN SECONDS rather than
   986	      /// library iRT units -- i.e. after the retention-time map has been fitted
   987	      /// and applied. It centres pass 2's extraction window; no sub-score
   988	      /// reads it any more (RT_DELTA was removed). Before it,
   989	      /// the comparison would be between two different units.
   990	      bool library_rt_is_run_seconds = false;
   991	
   992	      /// Observed 1/K0 per precursor, indexed as the library is, or empty.
   993	      /// Supplied by the caller because the scorer never sees spectra.
   994	      const std::vector<float>* observed_im = nullptr;
   995	
   996	      /// Sub-score indices to withhold from the classifier, by name on the
   997	      /// command line. Ablation, which the project has needed for a while and
   998	      /// faked twice by comparing different binaries -- a comparison that also
   999	      /// changes whatever else moved between them.
  1000	      ///
  1001	      /// Implemented by flattening the column to a constant, so the existing
  1002	      /// constant-column guard drops it and SAYS SO in the log. A NaN column
  1003	      /// would be silently imputed somewhere and the arm would claim to have
  1004	      /// ablated something it did not.
  1005	      std::vector<int> disabled_sub_scores;
  1006	
  1007	      /// Draw each decoy's best score from as many candidates as a target has.
  1008	      /// See `Scoring::LDAParams::match_decoy_candidate_counts`.
  1009	      bool match_decoy_candidate_counts = false;
  1010	      bool fold_pool_rank = false;   ///< v1.13: within-fold rank pooling (LDAParams::fold_pool_rank)
  1011	
  1012	      /// Semi-supervised loop knobs, previously reachable only by recompiling.
  1013	      /// 0 / negative means "leave the LDAParams default alone".
  1014	      double train_fdr_initial = 0.0;
  1015	      double train_fdr = 0.0;
  1016	      /// -1 means "leave the LDAParams default alone". 0 is a REAL request -- seed-only, no
  1017	      /// semi-supervised iterations -- which the old `default 0, guard > 0` arrangement could
  1018	      /// not express (the one reviewer-flagged k value the loop supports but the CLI could not
  1019	      /// reach).
  1020	      int classifier_iterations = -1;
  1021	      /// Mechanism-5 composition stop (`Scoring::LDAParams::stop_on_composition`): implemented
  1022	      /// and tested in lda.h since it landed, but never reachable from the CLI -- every run to
  1023	      /// date ran the fixed-count loop. 0 / negative jaccard means "leave the LDAParams default".
  1024	      bool classifier_stop_on_composition = false;
  1025	      /// Negative means "leave the LDAParams default alone". 0 is a REAL boundary value (stop
  1026	      /// only on shrink/cap, never on Jaccard) -- with 0 as the sentinel an explicit
  1027	      /// `-classifier_stop_jaccard 0` passed validation and was then silently swallowed back
  1028	      /// to 0.98.
  1029	      double classifier_stop_jaccard = -1.0;
  1030	      /// Repaired collapse policy (high-water + patience); 0 / negative leave the legacy
  1031	      /// single-shot LDAParams defaults. See AnchorTrainingParams::shrink_floor.
  1032	      double classifier_stop_shrink_floor = 0.0;
  1033	      int classifier_stop_patience = 0;
  1034	      /// Stderr-only per-(fold, iteration) churn line; never changes output bytes.
  1035	      bool classifier_iteration_log = false;
  1036	      bool use_pi0 = false;
  1037	    };
  1038	
  1039	    struct PeakGroup
  1040	    {
  1041	      std::uint32_t precursor = 0;
  1042	      float apex_rt = 0.0f;
  1043	      float left_rt = 0.0f;
  1044	      float right_rt = 0.0f;
  1045	      float apex_intensity = 0.0f;
  1046	      /// Median m/z deviation, ppm, over this group's matched fragment peaks,
  1047	      /// and how many contributed. NaN/0 when the extractor did not collect it.
  1048	      /// Deliberately NOT a sub-score -- see where it is filled.
  1049	      float mass_ppm = std::numeric_limits<float>::quiet_NaN();
  1050	      std::uint16_t mass_ppm_n = 0;
  1051	
  1052	      /// SCATTER of this group's PER-FRAGMENT deviations, ppm, as a robust
  1053	      /// sigma (MAD x 1.4826) about the group's own median. NaN until at least
  1054	      /// `MassWidth::min_fragments` fragments contributed.
  1055	      ///
  1056	      /// This is the quantity a window width must be sized from, and it is NOT
  1057	      /// the spread of `mass_ppm` across groups. `mass_ppm` is a median over
  1058	      /// (fragment x cycle) cells, so its own precision is sigma/sqrt(N_eff) --
  1059	      /// sizing a window from how tightly group medians cluster would give a
  1060	      /// number several times too narrow and would look, wrongly, like a very
  1061	      /// well calibrated instrument. What has to fit inside the window is one
  1062	      /// FRAGMENT's deviation, so one fragment is the unit measured here.
  1063	      float mass_ppm_spread = std::numeric_limits<float>::quiet_NaN();
  1064	
  1065	      /// The run's OBSERVED 1/K0 for this group, intensity-weighted over the
  1066	      /// candidate's cycles. NaN on a run with no ion mobility.
  1067	      float observed_im = std::numeric_limits<float>::quiet_NaN();
  1068	      /// Robust sigma (MAD x 1.4826) of this group's per-fragment observed
  1069	      /// 1/K0, and how many fragments carried one.
  1070	      float im_spread = std::numeric_limits<float>::quiet_NaN();
  1071	      std::uint8_t im_frags = 0;
  1072	      /// How many fragments contributed to `mass_ppm_spread`.
  1073	      std::uint8_t mass_ppm_frags = 0;
  1074	      std::vector<double> sub_scores;
  1075	
  1076	      double dscore = 0.0;
  1077	      double qvalue = 1.0;
  1078	      double pep = 1.0;
  1079	      bool decoy = false;
  1080	    };
  1081	
  1082	    /// One accepted fragment's mass residual, tagged with the group it came
  1083	    /// from so the FDR can filter it afterwards.
  1084	    ///
  1085	    /// The tag is necessary because a residual is produced while scoring, and
  1086	    /// whether its group is worth fitting from is only known after `finish()`.
  1087	    /// Fitting a calibration from every candidate would fit it from the
  1088	    /// interference too.
  1089	    struct MassAnchor
  1090	    {
  1091	      /// `residual.decoy` is NOT populated here -- target/decoy status belongs
  1092	      /// to the group, and holding a second copy on the residual is what let
  1093	      /// the two disagree once. `acceptedMassResiduals` fills it from the group.
  1094	      MassResidual residual;
  1095	      std::uint32_t group = 0;   ///< index into Result::groups
  1096	    };
  1097	
  1098	    struct Result
  1099	    {
  1100	      std::vector<PeakGroup> groups;
  1101	
  1102	      /// `N_FRAGVEC` float32 per group, row-major, parallel to `groups` --
  1103	      /// empty unless `Options::fragvec` was on.
  1104	      ///
  1105	      /// A flat side array rather than a member of PeakGroup: an empty
  1106	      /// std::vector on the group would still cost 24 B x 21.8 M groups on
  1107	      /// every run that does NOT ask for this, and a fixed array would cost
  1108	      /// 312 B. Here the flag-off run pays one empty vector for the whole
  1109	      /// result. `finish()` permutes it with the groups; anything that
  1110	      /// reorders `groups` must reorder this too or the rows silently swap.
  1111	      std::vector<float> fragvec;
  1112	
  1113	      /// Per-fragment mass residuals, when `collect_mass_anchors` was on.
  1114	      ///
  1115	      /// One entry per (retained candidate x contributing fragment), NOT per
  1116	      /// accepted group -- acceptance is decided after these are produced. Use
  1117	      /// `acceptedMassResiduals` to reduce them.
  1118	      std::vector<MassAnchor> mass_anchors;
  1119	
  1120	      /// Residuals discarded because `max_mass_anchors` was reached. Non-zero
```
### src/score/PeakGroupScorer.cpp lines 1170-1200 (fragvecNames), 1240-1420 (Gate C, rejects registry, Session::add head), 1780-1830 (frag_at_apex), 2380-2730 (MS1 read, FRAGVEC rung (i), group push, refits), 3120-3200 (finish: sort, re-index, fragvec reorder)
```
  1170	      "var_rt_spread", "var_mass_survival", "var_null_control",
  1171	      "var_ref_corr_sum",
  1172	      "var_ref_corr_1", "var_ref_corr_2", "var_ref_corr_3", "var_ref_corr_4",
  1173	      "var_ref_corr_5", "var_ref_corr_6", "var_ref_corr_7", "var_ref_corr_8",
  1174	      "var_ref_corr_9", "var_ref_corr_10", "var_ref_corr_11", "var_ref_corr_12",
  1175	      "var_sig_share_1", "var_sig_share_2", "var_sig_share_3",
  1176	      "var_sig_share_4", "var_sig_share_5", "var_sig_share_6",
  1177	      "var_cand_rank", "var_cand_count"};
  1178	    return names;
  1179	  }
  1180	
  1181	  const std::vector<std::string>& PeakGroupScorer::fragvecNames()
  1182	  {
  1183	    // Built rather than spelled out: 72 of the 78 are `<base>_<k>` for k = 1..12
  1184	    // and writing them by hand is 72 chances to transpose a digit in a column
  1185	    // order that a comparison against the reference builder would then report
  1186	    // as an arithmetic disagreement.
  1187	    static const std::vector<std::string> names = [] {
  1188	      std::vector<std::string> n;
  1189	      n.reserve(N_FRAGVEC);
  1190	      for (const char* base : {"R1_LOGAREA", "R1_SHARE", "R1_LOGRATIO",
  1191	                               "R1_ATAPEX", "R1_MEAS", "R1_ABSENT"})
  1192	      {
  1193	        for (int k = 1; k <= 12; ++k)
  1194	        { n.push_back(std::string(base) + "_" + std::to_string(k)); }
  1195	      }
  1196	      for (const char* scalar : {"R1_N_MEAS", "R1_N_ABSENT", "R1_N_PRESENT",
  1197	                                 "R1_LOGTOT", "R1_LIB_CORR", "R1_LIB_CORR_LOO"})
  1198	      { n.emplace_back(scalar); }
  1199	      return n;
  1200	    }();
   ...
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
  1283	
  1284	      /// Record one decision. The path is passed in rather than stored because
  1285	      /// this struct is a file-scope singleton declared before Options is in
  1286	      /// scope here; the file is opened on first use and closed at exit.
  1287	      void note(const std::string& path, std::uint32_t precursor, bool is_decoy,
  1288	                double stat, bool admitted, bool was_ready)
  1289	      {
  1290	        if (path.empty()) { return; }
  1291	        std::lock_guard<std::mutex> g(log_mu);
  1292	        if (log == nullptr)
  1293	        {
  1294	          // APPEND, not truncate. Each Session owns its own calibration and so
  1295	          // opens this file independently; with "w" the second pass wiped the
  1296	          // first pass's decisions and the log appeared to contain a single tau.
  1297	          // The header is written only into an empty file.
  1298	          log = std::fopen(path.c_str(), "a");
  1299	          if (log == nullptr) { return; }
  1300	          if (std::ftell(log) == 0)
  1301	          { std::fprintf(log, "precursor\tdecoy\tstatistic\ttau\ttau_ready\tadmitted\n"); }
  1302	        }
  1303	        std::fprintf(log, "%u\t%d\t%.6g\t%.6g\t%d\t%d\n",
  1304	                     precursor, is_decoy ? 1 : 0, stat, tau,
  1305	                     was_ready ? 1 : 0, admitted ? 1 : 0);
  1306	        std::fflush(log);
  1307	      }
  1308	    };
  1309	
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
  1331	        {
  1332	          t.scans[c] += r->scans[c];
  1333	          t.too_few_present[c] += r->too_few_present[c];
  1334	          t.too_few_transitions[c] += r->too_few_transitions[c];
  1335	          t.below_corr[c] += r->below_corr[c];
  1336	          t.reference_zero[c] += r->reference_zero[c];
  1337	          t.not_local_max[c] += r->not_local_max[c];
  1338	          t.below_apex_evidence[c] += r->below_apex_evidence[c];
  1339	          t.outside_margin[c] += r->outside_margin[c];
  1340	          t.reached[c] += r->reached[c];
  1341	        }
  1342	        for (int c = 0; c < 2; ++c) { t.no_points[c] += r->no_points[c]; }
  1343	        for (int c = 0; c < 2; ++c) { t.empty_trace[c] += r->empty_trace[c]; }
  1344	        for (int c = 0; c < 2; ++c) { t.gate_c[c] += r->gate_c[c]; }
  1345	        for (int c = 0; c < 2; ++c) { t.few_excursions[c] += r->few_excursions[c]; }
  1346	        for (int c = 0; c < 2; ++c) { t.zero_trace[c] += r->zero_trace[c]; }
  1347	        t.too_few_at_apex += r->too_few_at_apex;
  1348	        t.masked_candidates += r->masked_candidates;
  1349	        t.masked_fragments += r->masked_fragments;
  1350	      }
  1351	      return t;
  1352	    }
  1353	  }
  1354	
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
  1407	    // D8: each transition standardised against its own local noise before
  1408	    // summing, so no transition dominates by being loud and none is boosted
  1409	    // by what the library expects. See the option's comment for why library
  1410	    // weighting was rejected.
  1411	    std::size_t excursions = 0;
  1412	    const auto total = options.noise_normalised_picking
  1413	      ? noiseNormalisedTrace(chromatogram, points,
  1414	                             options.empty_trace_sigma, &excursions)
  1415	      : summedTrace(chromatogram, points);
  1416	
  1417	    // Does anything rise above this precursor's own noise?
  1418	    //
  1419	    // The old test was `sum(median-subtracted trace) <= 0.0`. That sum is
  1420	    // ZERO-MEAN for a precursor with no real peak, so it was a coin flip on the
   ...
  1780	      }
  1781	
  1782	      // D4: a per-transition local background, subtracted before anything
  1783	      // compares observed intensities to the library. Without it `observed`
  1784	      // is a raw area over a 60 s window and is dominated by baseline and
  1785	      // interference, which is why it correlated with nothing. Clamped at 0
  1786	      // rather than allowed negative: a weak real fragment sitting below its
  1787	      // own local median is absent evidence, not negative evidence.
  1788	      std::vector<double> corrected(tc, 0.0);
  1789	      // Retained per fragment: RT_SPREAD needs each fragment's own baseline to
  1790	      // place its centroid, and recomputing localBackground for that would be
  1791	      // the same scan twice. Named for the fragment axis -- a plain
  1792	      // `background` collides with the scalar one the sub-scores below use.
  1793	      std::vector<double> frag_background(tc, 0.0);
  1794	      std::size_t at_apex = 0;
  1795	      // -out_fragvec only. R1_ATAPEX_k is the PER-FRAGMENT form of the counter
  1796	      // below, and this is the only place it exists: `at_apex` collapses it to
  1797	      // a scalar in the same statement that computes it, so the 12 ATAPEX
  1798	      // columns cannot be recovered at the per-candidate block downstream.
  1799	      // Empty, and never touched, with the flag off.
  1800	      std::vector<double> frag_at_apex;
  1801	      if (options.fragvec) { frag_at_apex.assign(tc, 0.0); }
  1802	      for (std::uint32_t k = 0; k < tc; ++k)
  1803	      {
  1804	        const std::uint32_t n = chromatogram.pointCount(k);
  1805	        const float* points_k = n ? chromatogram.trace(k) : nullptr;
  1806	        std::vector<double> whole(n, 0.0);
  1807	        for (std::uint32_t j = 0; j < n; ++j) { whole[j] = points_k[j]; }
  1808	        const double bg = localBackground(whole, lo, hi);
  1809	        frag_background[k] = bg;
  1810	        corrected[k] = std::max(0.0, observed[k] - bg * static_cast<double>(width));
  1811	        if (cand.apex < n && points_k[cand.apex] > bg)
  1812	        {
  1813	          ++at_apex;
  1814	          if (options.fragvec) { frag_at_apex[k] = 1.0; }
  1815	        }
  1816	      }
  1817	
  1818	      // RT_SPREAD: do this group's fragments agree about WHEN they elute?
  1819	      //
  1820	      // Each informative fragment gets a background-subtracted, intensity-
  1821	      // weighted retention-time centroid over the candidate's own boundaries;
  1822	      // the feature is the weighted scatter of those centroids about the
  1823	      // group's apex time. A peptide's fragments come from one ion packet and
  1824	      // agree; an interferent belongs to a different species and does not.
  1825	      //
  1826	      // The scatter is about the fragments' OWN weighted mean -- inter-fragment
  1827	      // dispersion, not displacement from `cand.apex`, which cancels out of the
  1828	      // arithmetic. Deliberate: the apex comes from the SUMMED trace, so one
  1829	      // loud interferent would move the reference and the measurement together,
  1830	      // and a feature about internal consistency must not be anchored to
   ...
  2380	          for (const double v : per_fragment) { abs_dev.push_back(std::abs(v - centre)); }
  2381	          const std::size_t m = abs_dev.size() / 2;
  2382	          std::nth_element(abs_dev.begin(), abs_dev.begin() + m, abs_dev.end());
  2383	          // 1.4826 makes the MAD an estimate of sigma for a Gaussian, so the
  2384	          // width the driver derives from it can be stated in sigmas.
  2385	          g.mass_ppm_spread = static_cast<float>(1.4826 * abs_dev[m]);
  2386	          g.mass_ppm_frags = static_cast<std::uint8_t>(
  2387	            std::min<std::size_t>(per_fragment.size(), 255));
  2388	        }
  2389	      }
  2390	
  2391	      // MASS_SPREAD is final here: a within-group scatter needs no context.
  2392	      // MASS_ACCURACY holds the RAW deviation for now and is re-centred against
  2393	      // the run's median in finish(), once every group has been seen.
  2394	      // Lower is better, so negate -- same convention as MASS_SPREAD.
  2395	      g.sub_scores[IM_SPREAD] = std::isfinite(g.im_spread)
  2396	        ? -static_cast<double>(g.im_spread)
  2397	        : std::numeric_limits<double>::quiet_NaN();
  2398	      g.sub_scores[MASS_SPREAD] = std::isfinite(g.mass_ppm_spread)
  2399	        ? -static_cast<double>(g.mass_ppm_spread)
  2400	        : std::numeric_limits<double>::quiet_NaN();
  2401	      g.sub_scores[MASS_ACCURACY] = std::isfinite(g.mass_ppm)
  2402	        ? static_cast<double>(g.mass_ppm)
  2403	        : std::numeric_limits<double>::quiet_NaN();
  2404	      // Already negated where it was computed; NaN passes through so the
  2405	      // constant-column guard can drop it on a run where it never fires.
  2406	      g.sub_scores[RT_SPREAD] = rt_spread;
  2407	
  2408	      // MS1_COELUTION: does the PRECURSOR rise and fall with its fragments?
  2409	      //
  2410	      // Over the candidate's own cycles, pair each cycle's summed fragment
  2411	      // intensity with the MS1 monoisotopic intensity at the nearest MS1 bin,
  2412	      // and correlate. The MS1 grid is ~1.8 s on IH1 against a ~0.4 s MS2
  2413	      // cycle, so several cycles map to one bin -- that is a real resolution
  2414	      // limit of the survey scan, not an approximation to be apologised for,
  2415	      // and a 20-30 s peak still spans ~15 bins.
  2416	      //
  2417	      // NaN, not zero, when there is no MS1 or no signal: zero is a legitimate
  2418	      // correlation (a precursor whose trace is flat where the fragments peak
  2419	      // is EVIDENCE AGAINST), and collapsing the two would feed the classifier
  2420	      // a placeholder dressed as a measurement.
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
  2439	          // Refuse a correlation that would be computed from a handful of
  2440	          // points: over 3 or 4 cycles almost anything correlates, and the
  2441	          // classifier cannot tell a well-supported 0.9 from a lucky one.
  2442	          bool any = false;
  2443	          for (const double v : m) { if (v > 0.0) { any = true; break; } }
  2444	          if (any && f.size() >= 5) { const double c = pearson(f, m);
  2445	                                      if (std::isfinite(c)) { r = c; } }
  2446	
  2447	          // WHICH guard kills it? On Astral this sub-score came out constant
  2448	          // across all 96,259 rows and was dropped as carrying no information --
  2449	          // and it is the one feature measured to discriminate (13.7x top bin),
  2450	          // so "it is NaN" is not a good enough answer. Atomics, not
  2451	          // thread_local: a thread_local counter read from one thread already
  2452	          // produced a retracted conclusion in this file once.
  2453	          ms1_no_signal.fetch_add(any ? 0 : 1, std::memory_order_relaxed);
  2454	          ms1_too_short.fetch_add(f.size() >= 5 ? 0 : 1, std::memory_order_relaxed);
  2455	          if (any && f.size() >= 5)
  2456	          {
  2457	            // Zero variance in the MS1 leg makes Pearson undefined however much
  2458	            // signal is present -- distinct from "no signal", and the case that
  2459	            // fires if several MS2 cycles share one MS1 bin.
  2460	            bool flat = true;
  2461	            for (std::size_t k = 1; k < m.size(); ++k)
  2462	            { if (m[k] != m[0]) { flat = false; break; } }
  2463	            ms1_flat.fetch_add(flat ? 1 : 0, std::memory_order_relaxed);
  2464	            ms1_ok.fetch_add(std::isfinite(pearson(f, m)) ? 1 : 0, std::memory_order_relaxed);
  2465	            ms1_bins_spanned.fetch_add(distinctBins_(chromatogram, lo, hi, *options.ms1),
  2466	                                       std::memory_order_relaxed);
  2467	            ms1_spans.fetch_add(1, std::memory_order_relaxed);
  2468	          }
  2469	        }
  2470	        else
  2471	        {
  2472	          // Split the guard. The first cut said "11,877 no MS1 available, 0
  2473	          // everything else" while the run had just built a 148 GiB MS1 matrix
  2474	          // over 9,983,789 precursors -- so the composite condition is useless
  2475	          // and each term has to be counted on its own.
  2476	          ms1_unavailable.fetch_add(1, std::memory_order_relaxed);
  2477	          if (options.ms1 == nullptr) { ms1_null.fetch_add(1, std::memory_order_relaxed); }
  2478	          else
  2479	          {
  2480	            if (options.ms1->empty()) { ms1_empty.fetch_add(1, std::memory_order_relaxed); }
  2481	            if (i >= options.ms1->precursors())
  2482	            { ms1_index_oob.fetch_add(1, std::memory_order_relaxed); }
  2483	            if (!(hi > lo)) { ms1_degenerate_span.fetch_add(1, std::memory_order_relaxed); }
  2484	          }
  2485	        }
  2486	        g.sub_scores[MS1_COELUTION] = r;
  2487	      }
  2488	      // ---- -out_fragvec: rung (i) of the sealed fragment-evidence contract ----
  2489	      //
  2490	      // 78 float32 per candidate: six 12-long vectors by LIBRARY-INTENSITY rank
  2491	      // (LOGAREA, SHARE, LOGRATIO, ATAPEX, MEAS, ABSENT) and six scalars
  2492	      // (N_MEAS, N_ABSENT, N_PRESENT, LOGTOT, LIB_CORR, LIB_CORR_LOO). Every
  2493	      // input is a local the scorer already computed and already discards; the
  2494	      // definition is analysis77/pick/wf_v33_fragvec_contract.md s.5.1 and the
  2495	      // arithmetic is wf_v33_fragvec_features.py's `row_features`, which is
  2496	      // what the measured result was produced with.
  2497	      //
  2498	      // Staged on the stack and appended at the push_back below, not here: the
  2499	      // `min_library_corr` gate a few lines down `continue`s AFTER this point,
  2500	      // and a row emitted for a candidate that never becomes a group would put
  2501	      // the export permanently out of step with -out.
  2502	      //
  2503	      // TWO PLACES THIS AND THE REFERENCE BUILDER CAN DIVERGE. Both are
  2504	      // MEASURED UNREACHABLE on the library the arm runs (dn_pred_cam.parquet,
  2505	      // 58,576,095 transitions), so they are recorded here to make a future
  2506	      // gate-A4 disagreement diagnosable in one step rather than fixed blind:
  2507	      //   (1) NaN. `corrected[k]` is built with std::max(0.0, x) above, which
  2508	      //       returns 0.0 for a NaN x, while the builder's np.maximum returns
  2509	      //       NaN. A single non-finite exported intensity would therefore give
  2510	      //       LOGAREA 0 here and NaN there (ABSENT agrees by accident: NaN > 0
  2511	      //       is false either way). That line is EXISTING scoring arithmetic --
  2512	      //       LIBRARY_CORR and the mass residuals read the same `corrected` --
  2513	      //       so it must not be changed for this export's convenience; the
  2514	      //       reference is what would have to move. Measured: 0 non-numeric and
  2515	      //       0 negative Intensity over 4,294,710 exported points.
  2516	      //   (2) A transition the extractor could not place (Product.Mz invalid)
  2517	      //       still occupies a rank HERE -- tc counts it, it enters tot, libsum
  2518	      //       and both Pearsons with corrected == 0 -- but contributes no row to
  2519	      //       -out_chrom, so a reference that rebuilds tc from the export would
  2520	      //       be shifted by one from that rank on. Measured: 0 such transitions
  2521	      //       (Product.Mz in [200.015427, 1799.998901], no NaN, none <= 0).
  2522	      float fv[N_FRAGVEC];
  2523	      if (options.fragvec)
  2524	      {
  2525	        // Its OWN permutation, and a STABLE one. The D6 block above also orders
  2526	        // by library intensity, but that order is `std::sort` (unstable) and is
  2527	        // destroyed with its braced block. Reusing it would assign ranks
  2528	        // differently from the reference builder's stable sort on every
  2529	        // Relative.Intensity tie, and the resulting permuted columns would read
  2530	        // as an arithmetic disagreement rather than as a different order.
  2531	        std::vector<std::uint32_t> lib_order(tc);
  2532	        std::iota(lib_order.begin(), lib_order.end(), 0u);
  2533	        std::stable_sort(lib_order.begin(), lib_order.end(),
  2534	                         [&](std::uint32_t a, std::uint32_t b_)
  2535	                         { return library_intensity[a] > library_intensity[b_]; });
  2536	        const std::size_t m12 = std::min<std::size_t>(tc, 12);
  2537	        double tot = 0.0, libsum = 0.0;
  2538	        for (std::uint32_t k = 0; k < tc; ++k)
  2539	        { tot += corrected[k]; libsum += library_intensity[k]; }
  2540	        const float nanf = std::numeric_limits<float>::quiet_NaN();
  2541	        double n_absent = 0.0;
  2542	        for (std::size_t r = 0; r < 12; ++r)
  2543	        {
  2544	          if (r >= m12)
  2545	          {
  2546	            // Rank not measured. NaN in the four value columns, ZERO in the two
  2547	            // indicators -- the one place the contract's two missing-data
  2548	            // conventions differ, and the reason MEAS/ABSENT are never NaN.
  2549	            fv[r] = nanf; fv[12 + r] = nanf; fv[24 + r] = nanf; fv[36 + r] = nanf;
  2550	            fv[48 + r] = 0.0f; fv[60 + r] = 0.0f;
  2551	            continue;
  2552	          }
  2553	          const std::uint32_t k = lib_order[r];
  2554	          const double share = tot > 0.0 ? corrected[k] / tot : 0.0;
  2555	          const double libshare = libsum > 0.0
  2556	            ? library_intensity[k] / libsum
  2557	            : 1.0 / static_cast<double>(tc);
  2558	          // A5: "measured absent" is a zero background-corrected 5-cycle AREA.
  2559	          // It says nothing about the raw trace, which is why ATAPEX above is
  2560	          // not forced to 0 here.
  2561	          const double absent = corrected[k] > 0.0 ? 0.0 : 1.0;
  2562	          n_absent += absent;
  2563	          fv[r] = static_cast<float>(std::log1p(corrected[k]));
  2564	          fv[12 + r] = static_cast<float>(share);
  2565	          fv[24 + r] = static_cast<float>(std::log((share + 1e-3) / (libshare + 1e-3)));
  2566	          fv[36 + r] = static_cast<float>(frag_at_apex[k]);
  2567	          fv[48 + r] = 1.0f;
  2568	          fv[60 + r] = static_cast<float>(absent);
  2569	        }
  2570	        fv[72] = static_cast<float>(m12);
  2571	        fv[73] = static_cast<float>(n_absent);
  2572	        fv[74] = static_cast<float>(static_cast<double>(m12) - n_absent);
  2573	        fv[75] = static_cast<float>(std::log1p(tot));
  2574	        // The same helper the LIBRARY_CORR sub-score uses: fewer than four
  2575	        // fragments with a positive corrected area report 0, not a Pearson over
  2576	        // three points.
  2577	        fv[76] = static_cast<float>(libraryCorrelation(corrected, library_intensity));
  2578	        // Leave-one-transition-out, FLOORED AT 0 by starting the max there --
  2579	        // the builder initialises at 0.0 and never lowers it, so a group whose
  2580	        // every LOO correlation is negative reports 0 rather than its maximum.
  2581	        // O(tc^2): tc Pearsons of tc-1 points, ~144 multiplies at tc = 12.
  2582	        double loo = 0.0;
  2583	        if (tc > 1)
  2584	        {
  2585	          std::vector<double> obs_loo(tc - 1), lib_loo(tc - 1);
  2586	          for (std::uint32_t drop = 0; drop < tc; ++drop)
  2587	          {
  2588	            std::size_t at = 0;
  2589	            for (std::uint32_t k = 0; k < tc; ++k)
  2590	            {
  2591	              if (k == drop) { continue; }
  2592	              obs_loo[at] = corrected[k];
  2593	              lib_loo[at] = library_intensity[k];
  2594	              ++at;
  2595	            }
  2596	            loo = std::max(loo, libraryCorrelation(obs_loo, lib_loo));
  2597	          }
  2598	        }
  2599	        fv[77] = static_cast<float>(loo);
  2600	      }
  2601	
  2602	      // A candidate whose spectrum does not resemble the library is not this
  2603	      // peptide, wherever it eluted.
  2604	      //
  2605	      // Measured on IH1 against DIA-NN's confident set, splitting our own
  2606	      // q<=0.01 calls by whether they land within 30 s of the true apex:
  2607	      //
  2608	      //   group             n     median library_corr   frac > 0.5
  2609	      //   on-RT targets   1360             0.582           54.1%
  2610	      //   off-RT targets  5724            -0.036           12.5%
  2611	      //   decoys          7007            -0.032           12.0%
  2612	      //
  2613	      // Off-RT targets and decoys are the SAME population on this feature. That
  2614	      // is why target-decoy FDR cannot see them: every target is trained as a
  2615	      // potential positive, so the classifier learns whatever separates off-RT
  2616	      // targets from decoys -- signal presence -- rather than library
  2617	      // agreement. A target at the wrong retention time still sits on real
  2618	      // co-eluting ions from a real peptide, so it looks alive on log_sn and
  2619	      // intensity_score in a way a shuffled decoy never does.
  2620	      //
  2621	      // This gate is safe for the FDR by the criterion doc/08 states for any
  2622	      // selection upstream of it: it must be label-symmetric. It is, and that
  2623	      // is measured rather than assumed -- 12.5% of off-RT targets and 12.0% of
  2624	      // decoys survive a 0.5 cut, while 54.1% of on-RT targets do. Targets and
  2625	      // decoys traverse identical code here.
  2626	      //
  2627	      // Off by default: it changes which candidates exist, so it must be turned
  2628	      // on deliberately and its effect measured, not inherited.
  2629	      if (options.min_library_corr > -1.0 &&
  2630	          g.sub_scores[LIBRARY_CORR] < options.min_library_corr)
  2631	      {
  2632	        ++result.candidates_below_library_corr;
  2633	        continue;
  2634	      }
  2635	      // Commit this candidate's staged mass anchors, now that it is certain to
  2636	      // become a group and its index is known. `decoy` is taken from the group
  2637	      // rather than the fragment: a residual's status is the status of the
  2638	      // sequence it was matched against.
  2639	      if (options.collect_mass_anchors && !staged_anchors.empty())
  2640	      {
  2641	        const std::uint32_t gi = static_cast<std::uint32_t>(result.groups.size());
  2642	        for (MassAnchor& a : staged_anchors)
  2643	        {
  2644	          if (result.mass_anchors.size() >= options.max_mass_anchors)
  2645	          { ++result.mass_anchors_dropped; continue; }
  2646	          a.group = gi;
  2647	          result.mass_anchors.push_back(a);
  2648	        }
  2649	        staged_anchors.clear();
  2650	      }
  2651	      // In lockstep with `groups`, and appended at the SAME statement so no
  2652	      // path can push one without the other.
  2653	      if (options.fragvec)
  2654	      { result.fragvec.insert(result.fragvec.end(), fv, fv + N_FRAGVEC); }
  2655	      result.groups.push_back(std::move(g));
  2656	    }
  2657	
  2658	    // Margin over this precursor's own runner-up, on the detector's correlation
  2659	    // sum. Computed here because it is the only place all of one precursor's
  2660	    // candidates are in hand; a per-candidate score cannot express "this
  2661	    // precursor had one obvious answer" versus "three equally plausible ones".
  2662	    mark(result.groups.size() > first_group ? TerminalReason::Scored
  2663	                                           : TerminalReason::AllCandidatesDropped);
  2664	    if (result.groups.size() > first_group)
  2665	    {
  2666	      double best = 0.0, second = 0.0;
  2667	      for (std::size_t g = first_group; g < result.groups.size(); ++g)
  2668	      {
  2669	        const double v = result.groups[g].sub_scores[CORR_SUM];
  2670	        if (v > best) { second = best; best = v; }
  2671	        else if (v > second) { second = v; }
  2672	      }
  2673	      const bool alone = (result.groups.size() - first_group) == 1;
  2674	      for (std::size_t g = first_group; g < result.groups.size(); ++g)
  2675	      {
  2676	        const double v = result.groups[g].sub_scores[CORR_SUM];
  2677	        // The sole candidate is compared against nothing, not against zero:
  2678	        // giving it its full corr_sum as a margin would make "only one peak
  2679	        // found" look like overwhelming evidence.
  2680	        result.groups[g].sub_scores[CANDIDATE_MARGIN] =
  2681	          alone ? 0.0 : (v >= best ? best - second : v - best);
  2682	        // P1: the rest of the competition context. Rank on the same CORR_SUM
  2683	        // this block already reads (1 = best), and the block size itself --
  2684	        // "the only peak found" and "best of seven near-ties" are different
  2685	        // evidence, and no per-candidate score can express either.
  2686	        std::size_t rank = 1;
  2687	        for (std::size_t h = first_group; h < result.groups.size(); ++h)
  2688	        { if (result.groups[h].sub_scores[CORR_SUM] > v) { ++rank; } }
  2689	        result.groups[g].sub_scores[CAND_RANK] = static_cast<double>(rank);
  2690	        result.groups[g].sub_scores[CAND_COUNT] =
  2691	          static_cast<double>(result.groups.size() - first_group);
  2692	      }
  2693	    }
  2694	  }
  2695	
  2696	  /// Fit the discriminant over the retained groups and assign q-values.
  2697	  ///
  2698	  /// Split out of `finish()` so it can be run more than once over the SAME
  2699	  /// groups. That is what makes iteration cheap: the candidate picker is
  2700	  /// Refit the discriminant over the groups already scored.
  2701	  ///
  2702	  /// This used to recompute RT_DELTA first -- the ONLY sub-score that depended
  2703	  /// on the fitted retention-time map, which is what made "refit the map, then
  2704	  /// rescore" cheap. RT_DELTA is gone (see the note where it was declared), so
  2705	  /// NO sub-score depends on the map any more and re-running this after a new
  2706	  /// map produces a bit-identical result. `refitsChangeScores()` says so, and
  2707	  /// the retention-time refinement loop asks before spending a round.
  2708	  ///
  2709	  /// The map still matters upstream: it centres pass 2's extraction window.
  2710	  /// That has already run by the time anything calls this.
  2711	  void PeakGroupScorer::refit(const Library& library, Result& result,
  2712	                              const Options& options)
  2713	  {
  2714	    if (result.groups.empty()) { return; }
  2715	    // Counters are recomputed by the fit; reset so they do not accumulate
  2716	    // across iterations.
  2717	    result.identified_at_1pct = 0;
  2718	    fitAndAssign_(library, result, options);
  2719	  }
  2720	
  2721	  /// Whether refitting after a new retention-time map can change any score.
  2722	  ///
  2723	  /// False since RT_DELTA was removed. Kept as a function rather than deleted
  2724	  /// at the call sites so that adding a map-dependent sub-score re-enables the
  2725	  /// loop by flipping one return, instead of by remembering that a loop was
  2726	  /// deleted somewhere.
  2727	  bool PeakGroupScorer::refitsChangeScores() { return false; }
  2728	
  2729	  void PeakGroupScorer::fitAndAssign_(const Library& library, Result& result,
  2730	                                      const Options& options)
   ...
  3120	                   subScoreNames()[static_cast<std::size_t>(idx)].c_str());
  3121	    }
  3122	
  3123	    // Peak groups arrive in whatever order the extractor finished their
  3124	    // precursors, which is retention-time order rather than library order. Put
  3125	    // them back in library order before anything downstream sees them: the rows
  3126	    // of the feature matrix are these groups, and a classifier fitted on rows
  3127	    // ordered by elution time is a classifier that can depend on it. Stable, so
  3128	    // a precursor's candidates keep the order the search produced.
  3129	    // Canonical order: (precursor, apex RT, apex intensity). NOT "precursor,
  3130	    // then whatever order the search produced".
  3131	    //
  3132	    // PART 4 step 0 of the OpenSWATH handoff requires (group, apex RT, feature
  3133	    // id) before anything else, and calls skipping it the determinism bug of
  3134	    // section 3.3 -- which that project rates its most important non-result.
  3135	    // Sorting by precursor alone leaves a precursor's candidates in discovery
  3136	    // order, so a thread-count change reorders rows, which reassigns folds,
  3137	    // which changes the fit. Two runs of the same input would then disagree for
  3138	    // a reason invisible in any output.
  3139	    //
  3140	    // apex_intensity substitutes for "feature id": we have no stable per-feature
  3141	    // identifier, and it breaks ties that apex RT alone leaves.
  3142	    //
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
  3158	
  3159	    if (!result.mass_anchors.empty())
  3160	    {
  3161	      std::vector<std::uint32_t> moved_to(order.size());
  3162	      for (std::size_t n = 0; n < order.size(); ++n) { moved_to[order[n]] = static_cast<std::uint32_t>(n); }
  3163	      for (MassAnchor& a : result.mass_anchors)
  3164	      { if (a.group < moved_to.size()) { a.group = moved_to[a.group]; } }
  3165	    }
  3166	
  3167	    {
  3168	      std::vector<PeakGroup> reordered;
  3169	      reordered.reserve(result.groups.size());
  3170	      for (const std::uint32_t o : order) { reordered.push_back(std::move(result.groups[o])); }
  3171	      result.groups.swap(reordered);
  3172	    }
  3173	
  3174	    // The fragvec rows are parallel to `groups`, so they take the SAME
  3175	    // permutation. This is the whole reason the 78 floats are retained rather
  3176	    // than streamed as each candidate finishes: the row order -out publishes --
  3177	    // and therefore the ordinal the export is keyed on -- does not exist until
  3178	    // the sort above has run. Costs one transient copy of 312 B x groups
  3179	    // (6.8 GB at 21.8 M rows, against the arm's measured 716 GB peak), the same
  3180	    // shape of transient the group reorder just above pays.
  3181	    if (!result.fragvec.empty())
  3182	    {
  3183	      std::vector<float> fv_reordered;
  3184	      fv_reordered.reserve(result.fragvec.size());
  3185	      for (const std::uint32_t o : order)
  3186	      {
  3187	        const float* src = result.fragvec.data() + static_cast<std::size_t>(o) * N_FRAGVEC;
  3188	        fv_reordered.insert(fv_reordered.end(), src, src + N_FRAGVEC);
  3189	      }
  3190	      result.fragvec.swap(fv_reordered);
  3191	    }
  3192	
  3193	    for (const auto& g : result.groups)
  3194	    {
  3195	      (g.decoy ? result.decoy_groups : result.target_groups) += 1;
  3196	    }
  3197	    if (result.groups.empty()) { return result; }
  3198	
  3199	    // No decoys means no negative class. The scorer would still return numbers
  3200	    // -- q = 0 for everything -- and they would be read as an FDR. Refuse
```
### src/OpenDIAlyzer.cpp lines 4400-4460 (scorer options incl. fragvec) and 4700-4760 (pass driver: Sink, extractInto_, finish)
```
  4400	      std::max(0, getIntOption_("max_mass_anchors")));
  4401	    options.min_library_corr = getDoubleOption_("min_library_corr");
  4402	    const std::string picker = getStringOption_("picker");
  4403	    options.union_picking = picker == "union" || picker == "union_openswath";
  4404	    options.openswath_picking = picker == "openswath" || picker == "union_openswath";
  4405	    // `-picker amplitude` used to be a NO-OP. It is in the valid-strings list,
  4406	    // so it validates and prints in the help, but the amplitude detector is
  4407	    // selected by the older `-amplitude_picking` flag and nothing here read
  4408	    // `picker == "amplitude"` -- so anyone choosing it from the documented
  4409	    // options silently got the co-elution picker instead, and any arm labelled
  4410	    // "amplitude" that was driven this way measured co-elution under a wrong
  4411	    // name. Both spellings now select it.
  4412	    options.coelution_picking =
  4413	      !getFlag_("amplitude_picking") && picker != "amplitude";
  4414	    options.openswath_sn = getDoubleOption_("openswath_sn");
  4415	    options.openswath_gauss = getFlag_("openswath_gauss");
  4416	    options.openswath_peak_width = getDoubleOption_("openswath_peak_width");
  4417	    options.min_corr_score = getDoubleOption_("min_corr_score");
  4418	    // Final scoring only: pass 1 always trains its own discriminant, so the
  4419	    // anchor harvest and the calibration stay native whatever model pass 2 is
  4420	    // scored with (see calibrating_).
  4421	    options.classifier_model_out = calibrating_ ? "" : getStringOption_("classifier_model_out");
  4422	    options.classifier_model_in = calibrating_ ? "" : getStringOption_("classifier_model_in");
  4423	    options.max_corr_diff = getDoubleOption_("max_corr_diff");
  4424	    options.max_candidates = static_cast<std::size_t>(
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
  4441	    options.classifier_iterations = getIntOption_("classifier_iterations");
  4442	    options.classifier_stop_on_composition = getFlag_("classifier_stop_on_composition");
  4443	    options.classifier_stop_jaccard = getDoubleOption_("classifier_stop_jaccard");
  4444	    options.classifier_stop_shrink_floor = getDoubleOption_("classifier_stop_shrink_floor");
  4445	    options.classifier_stop_patience = getIntOption_("classifier_stop_patience");
  4446	    if (options.classifier_stop_on_composition)
  4447	    {
  4448	      // Validated only when armed: a Jaccard outside [0,1] silently makes convergence impossible
  4449	      // (only shrink/cap could ever fire), which is a configuration error pretending to be a
  4450	      // scientific result. Same logic for the floor and patience.
  4451	      if (options.classifier_stop_jaccard < 0.0 || options.classifier_stop_jaccard > 1.0)
  4452	      { throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
  4453	          "classifier_stop_jaccard must be in [0,1]",
  4454	          std::to_string(options.classifier_stop_jaccard)); }
  4455	      if (options.classifier_stop_shrink_floor < 0.0 || options.classifier_stop_shrink_floor > 0.5)
  4456	      { throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
  4457	          "classifier_stop_shrink_floor must be in [0,0.5]",
  4458	          std::to_string(options.classifier_stop_shrink_floor)); }
  4459	      if (options.classifier_stop_patience < 1)
  4460	      { throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
   ...
  4700	      writeLogInfo_(m.str());
  4701	    }
  4702	  }
  4703	
  4704	  /// @p collect, when non-null, receives a copy of every chromatogram the
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
  4751	  // loop. Every scoring run now goes through extractAndScore_; a run that
  4752	  // wants the chromatogram dump passes a collector, see TeeChromatogramSink.
  4753	
  4754	  /// Everything the scoring stage has to say, whichever path produced it.
  4755	  void reportScoring_(const ODIA::PeakGroupScorer::Result& scored,
  4756	                      const std::string& classifier, double ms,
  4757	                      const std::string& how)
  4758	  {
  4759	    std::ostringstream msg;
  4760	    msg << how << " " << scored.groups.size() << " peak groups with "
```
