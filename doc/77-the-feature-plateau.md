# 77 — The feature plateau: samelib3 localises the e2e gap to evidence content

2026-08-30. Adversarially reviewed (codex gpt-5.6-sol at ultra + kimi 0.39);
every review attack answered by measurement in the same afternoon. Analysis
scripts: `shared/libv2/sl3_{rank,ceiling,ceiling2,ceiling3,fork}.py` +
cached arrays in `dax:/scratch/kohlbach/odia2x2/sl3_cache/`.

## The experiment chain

**samelib3** (run overnight, binary d7ab850d): ODIA on S08 with DIA-NN's OWN
library (`dn_pred_cam.parquet`, 4,961,341 precursors — the library of the
33,330 reference run), RT window pinned at 110 s to equal e2e_s08 by
construction. Correctly pre-prior: the `-im_prior` is peptdeep-fitted and must
not touch DIA-NN's IM column. Result: **13,735 at q<=0.01** (pass 1: 14,709).
samelib2's earlier collapse to ~4k is confirmed as the window, nothing else.

With the library confound gone, the funnel decomposes cleanly
(UniMod:4 -> Carbamidomethyl aliased join; verbatim joins lose 10% silently):

| stage | of DIA-NN's 33,330 |
|---|---|
| reach ODIA's scorer ("scored") | 32,552 = **97.7%** |
| absent from ledger | 646 (281 trace to `-min_library_fragments 3`, a policy knob) |
| no_candidate / few_points | 131 / 1 |
| in ODIA's top 13,735 (its own depth) | 11,950 (35.9%, precision 87%) |
| in ODIA's top 33,330 (matched depth) | 20,781 = **62.3%** |
| p50 / p75 / p90 rank in ODIA's ordering | 19,531 / 112,900 / **1,148,221** |

Admission is worth 2.3 points, no more. The gap is the RANKING. The mix
fixture's "ranking is fine (95.6%)" does not transfer to real library scale
(pre-selected library; and the fixture cannot arm Gate C at all — the gate
calibrates on the first 20,000 decoys and mix10k carries 9,995, so its
`gate C 0/0` means NEVER ARMED, not "no loss").

## The oracle ceiling — and its robustness

Oracle test: HistGB on ODIA's 22 sub-score features, POSITIVES = rows of
DIA-NN's precursors (the oracle part — the shipped scorer never sees these
labels), negatives = decoy rows, 5-fold out-of-fold BY PRECURSOR, ranked by
max out-of-fold probability per precursor (never min q). Both sides of the
comparison are cross-fitted — shipped DScore comes from Percolator-style
precursor-keyed CV groups (`PercolatorEngine.cpp:66-127`).

Recall of DIA-NN's set at depth 33,330, across every perturbation the two
reviews demanded:

| arm | factor changed | recall @33,330 |
|---|---|---|
| shipped DScore | — | 62.3% |
| v1 oracle | — | 65.2% |
| A: labels RT-matched ±20 s | label noise (47% of rows were off-apex) | 65.2% |
| A40: labels ±40 s | label tolerance | 65.2% |
| B: 600 iters, 255 leaves | capacity ×2 | 65.1% |
| E: early_stopping=False | v1 had silently stopped at n_iter_=117/300 (codex caught it) | 65.2% |
| C: negatives = 2:1 decoys + 2:1 non-DN targets | the decoy-only loophole | 65.2% |
| D1: bag-size-corrected max (shipped) | codex's 2.91-vs-1.88 rows/bag multiplicity bias | 62.4% |
| D2: top-2-mean aggregation (shipped) | max-pooling | 61.7% |

Pre-registered (in `sl3_ceiling3.py`'s header before running): the feature-
deficit wording survives only if no arm beats 65.2 by >2.0 points and D lifts
shipped by <=1.0. Every arm landed within 0.1. This is the plateau codex's
review named as the requirement for a practical-ceiling claim.

Still unrun from the review list: an XGBoost/ranking-loss learner and a
peptide-grouped (cross-charge) outer CV. Both reviews agree cross-charge
leakage, if present, INFLATES the oracle — conservative for this conclusion.

## The fork: found-but-indistinguishable, not unfound

For every DIA-NN ID, does ODIA's candidate list contain a row within ±20 s of
DIA-NN's own apex? (`sl3_fork.py`)

| shipped rank stratum | DN IDs | near-RT candidate | none |
|---|---|---|---|
| top-13,735 | 11,950 | 11,922 | 28 |
| 13,736–33,330 | 8,831 | 8,743 | 88 |
| 33,331–200k | 4,817 | 4,595 | 222 |
| deeper | 6,954 | 4,684 | 2,270 |

**Picking failure is 8.0%** (2,608/32,552), concentrated in the deepest
stratum. For the other 92% the correct peak IS in the candidate list; its
features fail to separate it. Feature medians, near-RT rows, recovered vs
deep (delta/IQR): library_dotprod 0.92→0.64 (1.02), library_rmsd (1.00),
fragment_coverage 1.00→0.75 (1.00), library_corr 0.83→0.29 (0.87), corr_sum
7.4→3.8 (0.87), ms1_coelution 0.84→0.30 (0.83). Mass and IM contribute
almost nothing (0.17/0.05) — calibration is not the discriminator, agreement
and co-elution are. Six features have zero median separation; `var_log_sn`
sits at ln(100)=4.605 in both classes — the D5 floor cap, already documented
as a near-constant in `PeakGroupScorer.cpp:1847`.

## Committed conclusion (review-corrected wording)

On S08/samelib3, the shipped semi-supervised scorer is 2.9 points below a
supervised-oracle plateau of 65.2% that is invariant to label definition,
label tolerance, capacity, training depth, negative composition, and bag
aggregation. Multiple strong probes support a practical feature limitation:
the missing ~35% of DIA-NN's IDs are overwhelmingly precursors whose correct
peak was FOUND but whose extracted evidence — library agreement and
co-elution above all — does not separate them from the 3.7M-target bulk.
NOT claimed: an information-theoretic bound, or "no scorer can do better"
(both reviews); the FDR/threshold machinery is not validated by this metric
(codex finding 5 — the 13,735-vs-33,330 yield question also involves the
decoy null, measured separately in the funnel work).

## Next

Codex's paired counterfactual, (a) > (c) > (b): for deep-found DIA-NN IDs and
matched recovered controls, force re-extraction at DIA-NN's apex, then apply
interference-aware fragment masking, recompute the same features, rerank with
a frozen scorer. Forced-apex rescue → localisation; fragment-mask rescue →
interference; no rescue with coherent traces → representation; no coherent
trace → evidence genuinely absent. DIA-NN score injection is dropped as
circular (both reviews, independently).
