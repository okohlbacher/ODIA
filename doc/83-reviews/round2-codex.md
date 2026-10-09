# Codex gpt-5.6-sol (effort max) review #2 of ODIA standing, 2026-10-02 (verbatim; reviewed from the pasted brief only — codex tool host broken)

## Q1. The standing

1. **[HIGH] [BITE] “≈33 points with the correct peak in hand” is arithmetically overstated.**

   **VERIFIED:** Correct-peak availability is 92.0%, while native recall at depth is 64.086%; the attributable gap is therefore about 27.9 points, not 33. The ≈33-point figure is admission 97.7% minus ranking 64.1%, so it includes precursors whose correct peak is not in the candidate list at all (doc/83, pasted §2.1; `analysis77/pick/sb_fragvec_depth_result.txt`, pasted §7).

   **INFER:** Keep “ranking is the dominant bottleneck,” but rename it: “ranking/learning after admission.” Do not describe all ≈33 points as ranking with the right row available.

2. **[HIGH] [BITE] The +3.57 oracle gap changes the diagnosis, but it does not isolate a “training-recipe” defect.**

   **VERIFIED:** The native scorer selects positives as the current model’s best target rows passing a q-value generated from that same target/decoy ranking (`include/odia/scoring/lda.h:1144-1166`). The oracle instead receives DIA-NN-derived reference labels. The source itself explains that wrong-RT target rows can contain real peptide signal that separates them from synthetic decoys even though they are not the claimed peptide (`src/score/PeakGroupScorer.cpp:2613-2619`).

   **INFER:** The +3.57 points establish that the existing 42 columns contain more teacher-predictive information than native scoring extracts. They do **not** establish that ordinary hyperparameter tuning can recover it. The gap combines:

   - teacher labels versus positive–unlabelled self-training;
   - different negative construction and class balance;
   - possibly a different GBT implementation/objective;
   - fold pooling/calibration;
   - actual optimizer instability.

   The right frame is now: **candidate availability → label/objective gap → learner instability → missing evidence**.

3. **[HIGH] [BITE] The native seed assumption is internally contradicted by the scorer’s own forensic comment.**

   **VERIFIED:** The seed-fit rationale says false target rows are distributed like decoys and therefore only shrink the target/decoy mean difference without rotating it (`include/odia/scoring/lda.h:1087-1093`). Elsewhere, the scorer documents the opposite: wrong target rows lie on real co-eluting ions and separate from decoys through generic signal presence (`src/score/PeakGroupScorer.cpp:2613-2619`).

   **INFER:** This is probably the central structural failure. The seed can learn “real biological signal versus synthetic decoy” rather than “correct peptide versus interference,” and subsequent q-selected iterations reinforce that choice.

4. **[MED] [BITE] Fragment co-elution remains a strong mechanism, but it is overclaimed as the explanation for the whole residual.**

   **VERIFIED:** The shipped feature set already contains aggregate co-elution, robust-reference correlation vectors, signal-share vectors, and within-precursor competition fields (`include/odia/PeakGroupScorer.h:42-64`, `include/odia/PeakGroupScorer.h:343-380`). The oracle recovers +3.57 points without adding evidence, while FRAGVEC adds +2.80 beyond that (`analysis77/pick/sb_fragvec_depth_result.txt`, pasted §7).

   **INFER:** Co-elution failure is a measured phenotype of buried precursors, but the causal split between poor evidence, wrong labels, and poor optimization remains unresolved. The A/C classes being defined by native DScore makes the mechanistic attribution partly circular, as doc/83 already concedes.

5. **[MED] [BITE] The FRAGVEC decision is less statistically secure than “MDE” suggests.**

   **VERIFIED:** The registered “MDE” is merely the largest absolute difference between two fold seeds, not a sampling-error estimate, cross-run variance, or minimum detectable effect. Identity recall yields +2.799 and “CONDITIONAL”; concordant recall yields +3.024 and “FUND,” with the primary metric chosen only after the registration left it ambiguous (`analysis77/pick/sb_fragvec_depth_result.txt`, pasted §9).

   **INFER:** The user’s identity-recall decision resolves the procedural ambiguity, but the two-seed spread should be renamed a **seed-spread tolerance**. It provides no evidence about another sample, library, decoy generator, or native self-training.

6. **[MED] [BITE] “Peak RSS is a saturated budget” remains a hypothesis until sb1 reads.**

   **VERIFIED:** The current observation is one 400-GiB setting with a derived precursor cap; the unknown remainder and allocator ownership remain unmeasured. The running arm explicitly tests the prediction rather than having already established it (`shared/libv2/run_sb1.sh:16-21`, `shared/libv2/run_sb1.sh:29-43`).

   **INFER:** “Budget-dominated” is fair; “the budget is the memory” is premature. A budget-correlated pool, allocator-retained pages, and overlapping transients can all produce the observed peak. Also, sb1’s verdict calls the cut “free” while making its wall-time clause merely supporting; that conclusion would be false if memory falls but throughput regresses materially (`shared/libv2/run_sb1.sh:30-41`).

7. **[FINE]** The serial-wall diagnosis is credible as a prioritization, and the dangling `thread_local` registry finding is directly verified: pointers are registered without any unregister path and later dereferenced globally (`src/score/PeakGroupScorer.cpp:1310-1351`).

---

## Q2. The +3.57 oracle-versus-native gap

Ranked mechanisms:

1. **[HIGH] [BITE] Wrong learning target: self-confirming target/decoy training.**

   **VERIFIED:** Each iteration chooses the model’s current best row per group, computes q-values from those same scores, and treats qualifying target rows as positives (`include/odia/scoring/lda.h:1144-1166`). Negatives are the current best synthetic-decoy rows (`include/odia/scoring/lda.h:1262-1289`).

   **Why it can produce +3.57:** Native training never observes “correct DIA-NN precursor/row versus plausible wrong target row.” It sees “selected target versus synthetic decoy.” A wrong target peak containing real peptide signal can therefore become a very convincing positive.

   **Cheap test:** On the existing 42-column table, run a label ladder with identical folds and learner:

   - native pseudo-positives versus decoys;
   - DIA-NN-positive rows versus decoys;
   - DIA-NN-positive rows versus decoys plus entrapment targets;
   - correct target row versus the same precursor’s competing rows;
   - true targets versus empirical wrong-RT target rows, excluding synthetic decoys.

   This decomposes label target from learner implementation without extraction.

   **Structural fix:** Cross-fitted positive–unlabelled or listwise training using an independent multi-view anchor label, rather than q-values generated by the same model. Candidate-row labels should come from held-out RT/IM/mass/MS1 agreement, with the feature family used to create a label withheld from that fold’s learner.

2. **[HIGH] [BITE] Extreme-prior Newton boosting explains the compact/saturated flip.**

   **VERIFIED:** GBT starts at the logit of the selected positive/negative prior (`include/odia/scoring/gbt.h:300-305`). Gradients and Hessians are then computed at that extreme prior (`include/odia/scoring/gbt.h:325-337`). Leaf Newton steps are unbounded by default (`include/odia/scoring/gbt.h:83-90`, `include/odia/scoring/gbt.h:653-659`).

   **INFER:** At roughly 1:155 per fold, positive Hessians are tiny. Small changes in which bins isolate positives can create enormous early leaves, explaining fold-seed switches between max scores near 20 and 75–165. A leaf cap treats the symptom; failure on a second library is therefore unsurprising.

   **Cheap test:** Replay the exact table with:

   - equal total class weight, base margin zero;
   - deterministic negative subsampling matched by charge/RT/m/z/candidate count;
   - bounded steps with gain computed using the bounded step;
   - five or more fold assignments.

   The decisive read is collapse frequency, empirical-decoy-tail ranking, and entrapment performance—not maximum leaf size alone.

   **Structural fix:** A class-balanced objective, not only `max_delta_step`. Weight positives and negatives to equal total mass, set the intercept consistently, and scale `min_child_weight` for those weights.

3. **[HIGH] [BITE] Candidate-count matching does not govern training hard-negative selection.**

   **VERIFIED:** `match_decoy_candidate_counts` constructs `draw_from` (`include/odia/scoring/lda.h:797-828`) and is enabled by default unless explicitly disabled (`src/OpenDIAlyzer.cpp:4424-4426`). But training’s `best_rows_of()` scans all rows of every decoy group (`include/odia/scoring/lda.h:1005-1025`, `include/odia/scoring/lda.h:1148-1153`), and `top_decoys_only` reuses that unrestricted best row (`include/odia/scoring/lda.h:1280-1289`). The cap is applied only in final ranking/q-value assignment (`include/odia/scoring/lda.h:1525-1539`).

   **INFER:** The default claims to match best-of-N nulls but still trains against best-of-all decoys. That can materially harden and distort the negative class, particularly when candidate counts are asymmetric.

   **Cheap test/fix:** Apply the same deterministic `draw_from` prefix during every seed, ranking, and negative-selection scan. Replay only; no extraction required.

4. **[MED] [BITE] Fold scores remain badly calibrated in the tail.**

   **VERIFIED:** A held-out group is scored by only its one fold model (`include/odia/scoring/lda.h:1313-1334`). Native pooling normalizes each fold by the mean and standard deviation of its held-out best-decoy scores (`include/odia/scoring/lda.h:1337-1440`).

   **INFER:** Mean/SD normalization does not make saturated, heavy-tailed GBT margins comparable. The top of one fold can remain compressed or explosive while its bulk variance matches another fold.

   **Cheap test:** Preserve raw fold scores and transform each through its held-out empirical decoy survival function:

   \[
   s_f(x)=-\log_{10}\frac{1+\#\{d_f\ge x\}}{1+N_{d,f}}
   \]

   Pool those tail probabilities, then repeat across several independent fold partitions and average honest out-of-fold scores. The existing `fold_pool_rank` is not this: it constructs knots from both targets and decoys (`include/odia/scoring/lda.h:1385-1427`).

5. **[MED] [BITE] Iteration-dependent binning amplifies tiny feature changes globally.**

   **VERIFIED:** Default bin edges are rebuilt from the currently selected positive and negative rows (`include/odia/scoring/gbt.h:268-277`). The source explicitly notes that this makes every column’s discretization move after adding even a noise feature and records ID changes up to 11.9% (`include/odia/scoring/gbt.h:97-110`).

   **Cheap test:** Precompute edges once per training fold on all fold-training rows, reuse them through every semi-supervised iteration, and repeat the real/null append experiment.

   **Structural fix:** Stable bins should be a prerequisite for interpreting feature append arms. The current `fixed_bins` implementation recomputes full-data quantiles inside each fit; production code should cache them once rather than repeatedly sorting the entire table.

6. **[MED] [BITE] Missingness is silently converted to “exactly average” before GBT sees it.**

   **VERIFIED:** `lda.h` replaces every non-finite value with standardized zero (`include/odia/scoring/lda.h:460-509`). Consequently, GBT’s explicit missing-value binning (`include/odia/scoring/gbt.h:220-240`) is unreachable through this scorer.

   **INFER:** Existing IM/MS1 features and the new FRAGVEC values can lose informative missingness. This can also explain a difference from an oracle pipeline that uses native missing-value routing.

   **Cheap fix/test:** Preserve NaN into the binner or add explicit missing indicators before imputation. Compare exact replay and entrapment performance.

7. **[MED] [BITE] `max_depth` is off by one in effective tree capacity.**

   **VERIFIED:** The parameter says depth 4 permits at most 16 leaves (`include/odia/scoring/gbt.h:74-76`). But on `depth == D-1`, the implementation forces the current node to be a leaf instead of splitting it (`include/odia/scoring/gbt.h:541-595`); the nominal final level therefore receives no newly routed rows (`include/odia/scoring/gbt.h:616-633`). Depth 4 consequently permits at most eight populated leaves.

   **INFER:** This bites every fit, although the null `c6` result weakens it as the principal explanation for +3.57. It still invalidates comparisons to an oracle using conventional depth semantics.

8. **[LOW] [THEORETICAL unless armed] GBT seed masking is not implemented despite comments saying it is.**

   **VERIFIED:** The initial one-feature search honors `seed_mask` (`include/odia/scoring/lda.h:1034-1039`), but `fit_learner` applies the mask only to the NN; the GBT receives all columns (`include/odia/scoring/lda.h:959-987`). The later comment claims the mask applies to this seed fit (`include/odia/scoring/lda.h:1112-1114`).

   **INFER:** This is not the current gap if the mask is empty, but any GBT anti-circularity experiment using the mask is invalid.

9. **[LOW] [BITE only in degenerate fits] GBT reports a trained tree even when no split exists.**

   **VERIFIED:** `growTree_()` returns `any_split || true`, which is unconditionally true (`include/odia/scoring/gbt.h:633`). `fit()` consequently treats a constant stump as a successful learner (`include/odia/scoring/gbt.h:341-354`).

   **Fix:** Report whether a discriminative split or nontrivial margin was learned; otherwise the outer loop must count a skipped fit.

**Recommended structural sequence:** exact native replay → fix candidate-count matching → class-balanced objective and stable bins → empirical decoy-tail fold calibration → repeated out-of-fold ensemble → only then test alternative labels, monotone constraints, or FFNs. Early stopping on the decoy null alone is not sufficient: a one-class null cannot tell whether targets are correctly ranked. Use a disjoint inner group split with pseudo-positive/decoy partial AUC, entrapment holdout, and positive-set stability.

---

## Q3. FRAGVEC wiring

1. **[HIGH] [BITE] First resolve the exact feature contract.**

   **VERIFIED:** The engine computes 78 values and names them explicitly as six 12-vectors plus six scalars (`include/odia/PeakGroupScorer.h:385-399`, `src/score/PeakGroupScorer.cpp:1181-1199`). The registration’s input bullet describes `Ri` as “42 + 72 FRAGVEC columns,” while the standing describes 73 fitted columns.

   **INFER:** Before implementation, emit the exact retained-name list and a schema hash from the oracle NPZ. “Wire Ri” is otherwise ambiguous between 72 vectors, all 78 raw columns, and 73 columns surviving the oracle’s constant guard.

2. **[HIGH] [BITE] Decoy artefact risk is concentrated in library-ranked and library-agreement columns.**

   **VERIFIED:** The project documents that decoys copy their target’s per-fragment library intensities while fragment m/z values are recomputed (`include/odia/PeakGroupScorer.h:287-294`). FRAGVEC orders fragments by those library intensities (`src/score/PeakGroupScorer.cpp:2531-2536`), computes observed/library share ratios (`src/score/PeakGroupScorer.cpp:2553-2568`), and adds full and leave-one-out library correlations (`src/score/PeakGroupScorer.cpp:2573-2599`).

   Risk ranking:

   - **Highest:** `R1_LOGRATIO_*`, `R1_LIB_CORR`, `R1_LIB_CORR_LOO`.
   - **High:** rank-specific `LOGAREA`, `SHARE`, `ATAPEX`, and `ABSENT`; their rank is assigned using copied target intensities but their evidence comes from recomputed decoy fragments.
   - **Medium:** `N_ABSENT`, `N_PRESENT`, `LOGTOT`, which can learn generic real-signal-versus-synthetic-decoy differences.
   - **Lower but nonzero:** `MEAS` and `N_MEAS`; these mainly expose fragment count/padding, but can leak if decoy construction or invalid-product filtering changes counts. The present library reportedly has no invalid product transitions (`src/score/PeakGroupScorer.cpp:2516-2521`).

   Required tests:

   - Evaluate reference-positive versus **entrapment target** ranking without synthetic decoys in the metric.
   - Measure an AUC for synthetic decoys versus empirical target nulls, matched on charge, m/z, RT, fragment count, abundance, and candidate count. A high AUC is evidence of nonexchangeability.
   - Use wrong-RT or precursor-swapped target evidence as an empirical null.
   - Train on one decoy construction and evaluate on another.
   - Stratify entrapment FDP by each feature’s score quantile and by library intensity shape.
   - Require the benefit over a same-dimensional permuted control, not merely over the 42-column baseline.

3. **[HIGH] [BITE] A six-column null control is inadequate for a 78-column append.**

   **VERIFIED:** Iteration-dependent bins make even a pure-noise append change the rest of the model (`include/odia/scoring/gbt.h:97-110`).

   **INFER:** The engine control must append the same number of columns with the same marginal distributions, NaN patterns, rank covariance, and float32 rounding as the real block. Column-wise shuffling alone may destroy covariance and be too easy a control; use group-wise or within-stratum permutation as well.

4. **[HIGH] [BITE] Native self-training may recover none of the oracle gain.**

   **VERIFIED:** The oracle’s gain is reference-labelled; native positives are self-selected targets (`include/odia/scoring/lda.h:1144-1166`). Moreover, FRAGVEC’s incremental gain is only about +0.58 points at depth 20k and grows toward deeper ranks (`analysis77/pick/sb_fragvec_depth_result.txt`, pasted §7).

   **INFER:** The operational q≤0.01 head may gain much less than the +2.80 at depth, and the new columns may simply reinforce the current wrong pseudo-positives. This is exactly why matched-entrapment performance on both b2_cap1 and w1_ctl70 is the correct acceptance criterion.

5. **[MED] [BITE] Keep it final-pass-only initially.**

   **VERIFIED:** Existing model-in/model-out, fold pooling, transition masking, and FRAGVEC export are deliberately disabled during calibration and enabled only on the final pass (`src/OpenDIAlyzer.cpp:4418-4434`).

   **INFER:** The first integration should preserve pass-1 anchors byte-for-byte and change only final scoring. Computing FRAGVEC in both passes would confound feature gain with altered RT/IM/centring anchors and double the expensive work. A both-pass arm is justified only after final-only scoring passes.

6. **[HIGH] [BITE] Missingness semantics must be preserved explicitly.**

   **VERIFIED:** Padded ranks use NaN for value fields and zero for `MEAS`/`ABSENT`; measured zero-area fragments use finite zero-like values with `MEAS=1` and `ABSENT=1` (`src/score/PeakGroupScorer.cpp:2542-2568`). The scorer then mean-imputes every NaN to standardized zero (`include/odia/scoring/lda.h:460-509`).

   Required gates:

   - Constant-column elimination must occur before training and print every retained/dropped name.
   - Assert the reference arm retains exactly the oracle’s schema.
   - Preserve the float32 rounding used by the exported oracle block; silently switching to unrounded doubles can move bins and self-training trajectories.
   - Never drop a missing indicator merely because its paired value column varies.
   - Test libraries with fewer than 12 fragments and structural no-MS1/no-IM cases.

7. **[MED] [BITE] Memory cost is several copies, not “73 doubles per group.”**

   **VERIFIED:** `PeakGroup` stores its sub-scores in a per-group `std::vector<double>` (`include/odia/PeakGroupScorer.h:1074`). The learner allocates another row-wise double matrix for standardized values (`include/odia/scoring/lda.h:480-509`). The existing float32 sidecar costs 312 bytes per group and is transiently copied during reorder—reported as 6.8 GB at 21.8M rows (`src/score/PeakGroupScorer.cpp:3174-3190`).

   **INFER:** Appending the block directly to every group can add roughly 12–13 GiB per double copy, plus millions of vector reallocations, allocator overhead, the standardized copy, and possibly the float export. A 25–40 GiB peak increment is plausible. Keep a columnar float sidecar, drop constants before promotion to doubles, and do not imply `-out_fragvec` merely because scoring is enabled.

8. **[FINE]** The current sidecar/group lockstep and reorder logic are sound (`src/score/PeakGroupScorer.cpp:2651-2655`, `src/score/PeakGroupScorer.cpp:3174-3190`).

**Acceptance programme:**

- **Offline gate:** exact native-recipe replay on the existing table; 42, 42+real block, and 42+matched-null block over at least five fold seeds.
- **C0:** new binary, score flag off; scored TSV byte-identical to v0.5.0.
- **Cnull:** final-pass-only, 78 matched-null columns.
- **F78:** final-pass-only real block.
- **F78-repeat:** same arms on both b2_cap1 and w1_ctl70 lineages, plus human_v2 or another independent entrapment-bearing run.
- **Artifact control:** target-only reference/entrapment evaluation and, ideally, a different decoy construction.
- **Acceptance:** real block beats both C0 and Cnull at matched entrapment over the full registered FDP staircase; no material low-FDP ordinal regression; no increase in true FDP at own q≤0.01; no saturated-regime collapse in any registered seed; exact candidate rows/ordinals unchanged except through scoring; memory and scorer time reported.
- **Only afterward:** a both-pass arm to test whether better scoring improves calibration anchors.

---

## Q4. Next identification work, beyond FRAGVEC

1. **[HIGH] Native learner stabilization and repeated cross-fitting — highest expected IDs/week.**

   **Mechanism:** Class-balanced boosting, stable cached bins, matched hard negatives, empirical decoy-tail fold calibration, and an average of several honest out-of-fold partitions.

   **Cheap decision:** Recover a substantial fraction of the measured +3.57-point A-versus-native gap on the frozen 42-column table, while reducing the 16,965→4,105 seed collapse and entrapments. Kill if five-seed matched-entrapment performance is flat.

   **Risk:** It can stabilize the wrong target/decoy objective without learning correct peptide identity.

2. **[HIGH] Cross-fitted PU/listwise target construction.**

   **Mechanism:** Stop treating every model-selected target as a positive. Train candidate competition directly using orthogonal, fold-held-out RT/IM/mass/MS1 agreement and empirical target nulls, with targets otherwise treated as unlabelled.

   **Cheap decision:** On the existing feature/candidate table, measure purity and coverage of proposed positive labels against DIA-NN only in held-out families. Require improvement beyond repaired target/decoy GBT and beyond the same labels shuffled within strata.

   **Risk:** Hidden circularity. Every label-generating feature family must be excluded from the learner that consumes that label.

3. **[HIGH] Interference-robust recomputation from existing traces.**

   **Mechanism:** Fit a robust rank-one consensus or reference trace, identify fragment-specific residuals, downweight inconsistent fragments, and recompute library agreement/co-elution. This is deconvolution of the current candidate, not candidate re-seeking or a raw trace tensor.

   **Cheap decision:** Offline leave-one/two-fragment-out or robust-M-estimator features on existing chromatogram exports. Require incremental reference-labelled and entrapment gain beyond existing `REF_CORR_*`, which already uses a robust best-fragment reference (`include/odia/PeakGroupScorer.h:354-367`).

   **Risk:** Synthetic decoy artefacts and suppression of legitimately heterogeneous low-signal fragments.

4. **[HIGH] Rank-aligned per-fragment mass and mobility evidence.**

   **Mechanism:** Export per-fragment ppm residual, IM residual, tight-window survival, and availability indicators by library rank. Existing scoring exposes mainly group summaries such as mass accuracy/spread and IM spread (`include/odia/PeakGroupScorer.h:170-210`, `src/score/PeakGroupScorer.cpp:2391-2403`).

   **Cheap decision:** Output-only sidecar on retained groups, then the same oracle/null/entrapment ladder used for FRAGVEC. Kill if it adds less than a preregistered depth gain and no matched-entrapment improvement.

   **Risk:** Memory, missingness, and library-rank decoy artefacts.

5. **[MED] Buried-class ODIA-versus-DIA-NN XIC parity audit.**

   **Mechanism:** Compare trace by trace at the production coordinates for recovered and buried reference precursors. If ODIA specifically loses coherence in the buried class while DIA-NN retains it, fix extraction/resampling before inventing more scores.

   **Cheap decision:** A few thousand matched precursor-fragment traces, with identical RT support and normalization. Kill the entire extraction-fix lane if parity holds after abundance and interference matching.

   **Risk:** Export semantic or coordinate mismatches masquerading as extraction differences.

6. **[MED] Cross-precursor, co-located candidate competition.**

   **Mechanism:** Current competition fields compare candidates only within one precursor (`src/score/PeakGroupScorer.cpp:2658-2692`). Build a local graph of hypotheses sharing apex time, isolation window, and fragment evidence; penalize a weak hypothesis when the same evidence is better explained by a stronger neighbour.

   **Cheap decision:** On existing groups, measure how often buried/entrapment selections have a stronger conflicting neighbour after matching local density. Kill if the affected fraction is small or an oracle graph feature adds negligible held-out recall.

   **Risk:** Co-eluting true peptides share evidence too. This reopens “congestion” only if the explicit shared-evidence graph supplies new evidence; a scalar local-density feature remains closed.

7. **[MED] MS1 isotope and tight-mass channels.**

   **Mechanism:** Extract M, M+1, M+2 and tight/loose precursor traces, then compare isotope shape and co-elution with MS2. The current implementation uses one monoisotopic MS1 trace (`src/score/PeakGroupScorer.cpp:2408-2445`), while the header still notes isotope correlation as unavailable (`include/odia/PeakGroupScorer.h:30-33`).

   **Cheap decision:** Pilot a stratified candidate panel rather than a full dense matrix. Require conditional AUC inside abundance/charge strata and incremental oracle gain over current MS1 co-elution.

   **Risk:** Large memory multiplication, sparse survey scans, isotope interference.

8. **[LOW] Tiny FFN/DeepSets ensembles on engineered fragment vectors.**

   **Mechanism:** Smooth nonlinear interactions and ensemble several deterministic/bagged fits. The code already has an NN ensemble path and group bagging hooks (`include/odia/scoring/lda.h:151-168`, `include/odia/scoring/lda.h:946-984`).

   **Cheap decision:** Offline only, after repairing labels and folds. Compare against the repaired GBT over multiple seeds and samples.

   **Risk:** A stronger model amplifies the same decoy artefact and pseudo-label feedback. Do not fund full extraction arms unless it beats GBT on held-out entrapment.

I would not reopen m/z/IM windows, candidate generation, retained-group re-seeking, raw trace tensors, MS1-as-intensity, score injection, or per-run library fine-tuning without a new diagnostic specific to the previously failed mechanism.

---

## Q5. The three running arms

1. **[HIGH] [BITE] sb1’s verdict cannot call the cut “free” while excluding cost from PASS.**

   **VERIFIED:** M3 limits wall to 1.10× but is labelled supporting; PASS is A+M1+M2 only (`shared/libv2/run_sb1.sh:25-36`).

   **Fix to interpretation:** If M1/M2 pass but wall exceeds the cost gate, conclude “memory-saving and score-invariant,” not “free” or “new default.”

2. **[MED] [BITE] sb1 cannot attribute a fixed “remainder” from one peak subtraction.**

   **INFER:** `Peak − 1.065×budget` assumes allocations peak concurrently and scale linearly. Allocator retention, output-library writing, MS1, and extraction planes may peak at different times. Treat the subtraction as a descriptive residual until the stage sampler and instrumentation locate it.

3. **[FINE]** Raw `cmp` is a strong sb1 invariant because final groups are canonicalized by precursor, apex RT, and intensity before downstream scoring/output (`src/score/PeakGroupScorer.cpp:3123-3157`); the registered column-diff ladder is the right response to a failure.

4. **[HIGH] [BITE] C1’s `join >=95%` can hide nearly the entire allowed RT loss.**

   **VERIFIED:** C1 computes coverage only after requiring a 95% join (`shared/libv2/run_sbw.sh:22-25`).

   **INFER:** Missing precursors are likely nonrandom. Evaluate a fixed outer-joined 37,334-ID cohort and count missing entries as uncovered. Also report signed error and p50/p90/p95/p99 by charge, RT decile, m/z, abundance, and library-fragment count. The preregistered C1 should remain unchanged, but this read must accompany it.

5. **[HIGH] [BITE] C2’s uniform 120-anchor floor is far too weak for abundant charges.**

   **VERIFIED:** The same ≥120 threshold is applied to charges whose baseline counts range from 162 to 7,198 (`shared/libv2/run_sbw.sh:26-28`).

   **INFER:** z2 could lose more than 98% of its anchors and still pass. Add read-side anchor overlap/Jaccard, per-charge retained fraction, RT/m/z coverage, residual quantiles, and effective sample size. Likewise, “centring ≥50%” licenses a twofold collapse and cannot establish calibration equivalence.

6. **[HIGH] [BITE] C3 weakened the doc/83 gate and can hide localized damage.**

   **VERIFIED:** The running script uses region median ≥−1% and q≤0.01 ≥90% (`shared/libv2/run_sbw.sh:29-31`). The earlier plan additionally required at least 90% of registered ordinals to remain within −2% (doc/83, pasted §6 step 2).

   **INFER:** A median can pass while the low-FDP operating region or whole abundance/charge strata fail. Read every registered ordinal, the worst contiguous interval, acc/conc at depth, selected-apex changes, and entrapment counts. Do not retroactively change the verdict; mark these as secondary safety reads.

7. **[MED] [BITE] “sb2 passes and sb3 fails ⇒ 200 is the floor” is too strong.**

   **VERIFIED:** Only widths 110, 200, and baseline 400 are tested (`shared/libv2/run_sbw.sh:3-15`, `shared/libv2/run_sbw.sh:36-39`).

   **INFER:** This licenses “200 is the narrowest tested safe value,” not a floor. The actual boundary could be anywhere from 111 to 199 seconds.

8. **[MED] [BITE] Concurrent execution invalidates performance conclusions, though not deterministic score reads.**

   **VERIFIED:** The headers acknowledge concurrent sb1/sb2/sb3 execution and treat wall as supporting (`shared/libv2/run_sb1.sh:30-35`, `shared/libv2/run_sbw.sh:32-35`).

   **Add now:** Per-phase process CPU, `/proc/<pid>/io`, PSI, `memory.current/events`, `smaps_rollup` PSS/RssAnon/RssFile/Swap, NUMA placement, and achieved I/O bandwidth. A default speed claim still needs a solo repeat.

9. **[MED] [BITE] One seed cannot test whether narrowing the window changed regime susceptibility.**

   **INFER:** Report per-fold positive counts, iteration churn, score-tail quantiles, leaf maxima, and empirical decoy-tail calibration if available. “Regime label reported, not decided” is acceptable for this run, but a default must require the promised second seed and human_v2 replicate.

---

## Q6. Exploratory implementations

Recommended order:

1. **Instrumentation first.**
2. **Registry/counter ownership fix as a standalone correctness change.**
3. **Offline native-scorer decomposition and FRAGVEC artefact controls.**
4. **FRAGVEC final-pass-only integration.**
5. **MS1 parallel build, if instrumentation confirms the expected kernel.**
6. **Memory changes one at a time: ppm alias → optional IM alias → row compaction → byte chunker.**
7. **Parallel sink last.**

Findings and gates:

1. **[HIGH] [BITE] Instrumentation needs current RSS as well as monotone HWM.**

   `ru_maxrss` cannot identify which later stage still owns pages after an earlier peak. Print stage deltas for CPU, current RSS/PSS, HWM, task count, pool capacity, and allocator/mapping category. The plane-identity assertion should be behind a diagnostic flag because a full cell-by-cell scan can perturb timing.

   **Gate:** flag-off TSV byte identity; fixture overhead under a small preregistered bound; stage CPU totals reconcile with process CPU; external sampler and internal HWM agree.

2. **[HIGH] [BITE] Do not “fix” the rejects registry with only a destructor.**

   **VERIFIED:** It stores raw `thread_local` addresses globally without session ownership or unregistering (`src/score/PeakGroupScorer.cpp:1310-1351`).

   **INFER:** A destructor introduces global/thread-local destruction-order hazards and still leaves counters process-scoped. The pasted code also shows no reset in `totalRejects()`, so verify whether pass-2 reports include pass-1 counters.

   **Preferred design:** task/session-local `PickerRejects`, merged deterministically into the owning `Result`; no global registry.

   **Gate:** repeated worker create/destroy stress under ASan/TSan; exact counter totals across 1/8/64 threads and consecutive sessions.

3. **[HIGH] [BITE] FRAGVEC should not be represented as millions of enlarged per-group vectors.**

   Use a flat float32 sidecar, reorder it with groups, drop constants, then materialize only retained classifier columns. Separate `compute_fragvec_for_score` from `write_fragvec`.

   **Gate:** Q3’s C0/Cnull/F78 programme, exact retained schema, off-state byte identity, and measured HWM.

4. **[HIGH] [BITE] Frame-parallel MS1 can silently change residual selection and floating-point order.**

   A per-task `RESID_CAP` is not equivalent to a global serial cap unless merge truncates in original frame order. Multiple frames mapping to one output bin can also race or reorder sums.

   **Gate:** exact hash of the MS1 matrix, residual list/order, and downstream TSV at 1/8/64 threads; not merely equal residual medians. Require a meaningful full-stage speedup after accounting for serial decode.

5. **[HIGH] [BITE] Do not implement the aliasing/chunker/row-drop branch as one experimental unit.**

   - Alias `ppm_den` only after the diagnostic reports zero violations.
   - Alias `im_den` only after a separate zero-violation result.
   - Remove aliased writes exactly once; otherwise the shared plane will be double-incremented.
   - Test Sum/Max modes and no-mobility runs separately.
   - Row compaction needs a compact library→row mapping and an assertion that every removed row has zero reads.
   - The byte chunker needs to account for transitions, points, planes, pool-block rounding, and side arrays—not only mean cells.

   **Gates:** one same-cap control per individual change; cell hashes/asserts; raw TSV `cmp`; expected HWM delta; exact MS1 trace hash for row compaction; every precursor emitted exactly once for chunking.

6. **[HIGH] [BITE] Parallel sink is much more stateful than the current plan implies.**

   **VERIFIED:** `Session::add` and group push are currently justified by a serial emit invariant (`include/odia/PeakGroupScorer.h:920-930`). Gate C builds a stateful first-N-decoy null and changes behaviour once armed (`src/score/PeakGroupScorer.cpp:1244-1281`). Candidates, fragvec rows, and mass-anchor group indices are committed in lockstep (`src/score/PeakGroupScorer.cpp:2635-2655`).

   **INFER:** Private worker Sessions would each build a different Gate C null. Shared concurrent Sessions would race group vectors and make cap/order-dependent anchor selection nondeterministic. Ordered speculation must preserve:

   - Gate C observation and arming order;
   - mass-anchor cap order;
   - candidate order, including exact sort ties;
   - terminal reasons and reject counters;
   - group-index rebasing for anchors;
   - bounded speculative memory;
   - cancellation and exception semantics.

   **Gate:** equality across 1/8/64 threads for raw TSV, candidate ordinals, Gate C tau and sampled decoys, anchor identities/order, terminal reasons, reject counters, trained-model bytes, and all output sidecars. Only then apply the wall-time gate.

7. **[MED] [BITE] Start the scorer work now; it is missing from section 2C.**

   The same-column oracle gap is larger than the incremental FRAGVEC gain. Begin the exact replay, class-balance, matched-negative, cached-bin, empirical-fold-tail, and repeated-OOF experiments immediately. They need no extraction and can prevent several full arms from measuring classifier chaos.

8. **[MED] [BITE] Also start the buried-class XIC parity read now.**

   It is the cheapest way to decide whether the next identification week belongs to evidence engineering or extraction fidelity.

**What I would not do:**

- Do not run a bundled alias+chunker+row-drop acceptance arm.
- Do not alias `im_den` before the full diagnostic proves identity.
- Do not compute FRAGVEC in pass 1 initially.
- Do not begin parallel sink performance work before replacing global reject ownership and specifying Gate C/cap commit order.
- Do not fund a full MS1 parallel arm if instrumentation shows raw decode or memory bandwidth dominates and the build kernel cannot plausibly reach the target.

---

## Q7. What to stop

1. **[HIGH] [BITE] Stop spending full extraction node-days on scorer-only questions.** The candidates and existing feature rows are already available; replay class weighting, bins, folds, calibration, null appends, and most FRAGVEC decisions offline first.

2. **[HIGH] [BITE] Stop interpreting single-seed scorer arms.** A known 16,965→4,105 fold-seed flip is larger than almost every proposed feature gain. Require several offline seeds before one engine arm.

3. **[HIGH] [BITE] Stop using own q≤0.01 ID count as a decisive performance gate.** Native q is already known to correspond to materially higher true FDP. Matched-entrapment staircases must remain primary; q count is a reporting statistic.

4. **[HIGH] [BITE] Stop bundled performance changes.** Aliasing, row compaction, and byte chunking answer different questions and have different correctness surfaces.

5. **[MED] [BITE] Stop making speed claims from co-scheduled wall time.** Concurrent wall measurements are useful for operational sizing but cannot establish a kernel speedup. Use solo repeats or phase CPU/throughput.

6. **[MED] [BITE] Stop dimension-mismatched null controls.** Six null columns cannot control a 78-column feature append in a learner whose bins and split opportunities change globally when columns are added (`include/odia/scoring/gbt.h:97-110`).

7. **[MED] [BITE] Stop calling two-seed range an MDE or treating a DIA-NN-labelled oracle as native-transfer evidence.** Use it to screen information content; require native self-training, independent entrapment, and another sample for acceptance.

8. **[MED] [BITE] Stop full arms for `intercept_zero` or Gate C timing before offline simulation.** The former is completely replayable from the existing table; the latter can first be simulated from recorded Gate C statistics and candidate/anchor tables.

9. **[MED] [BITE] Stop using anchor-count floors detached from baseline retention.** A universal 120 permits catastrophic proportional losses in abundant charges. Use overlap, coverage, residuals, and retained fractions.

10. **[FINE]** Continue sb1 and the two width arms: they answer real zero-code questions. Just constrain their conclusions to score invariance, memory, and calibration noninferiority—not “free,” “floor,” or “default” until solo and cross-sample replication.
