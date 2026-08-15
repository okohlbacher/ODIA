# The oracle ladder: the whole pipeline works. The blocker is the prior.

Measured 2026-08-15, Astral. This is the discrimination `doc/16` §1 was written
to provide and `doc/19` §7 pre-registered as the response to a red Phase 3.

## The result

Same binary, same calibration, same windows, same gate. Only the library changes.

| library | targets | present | peak groups | **IDs at 1% FDR** |
|---|---|---|---|---|
| full parity (`odia_lib_parity`) | 9,983,789 | ~1.5% | 84,950 | **0** |
| **oracle** (DIA-NN's confident set) | **10,040** | ~100% | 49,907 | **pass 1 5,665 → pass 2 6,815** |

**68% recovery**, and pass 2 adds 1,150 over pass 1, so the two-pass refinement
works as designed.

## What this rules in and out

The four outcomes `doc/16` §1 named, and which one fired:

| outcome | fired? |
|---|---|
| fails even on a positive-only library → extraction/scoring/FDR **bug** | **no** |
| scores separate but q-values do not → competition/q-value **bug** | **no** |
| **positive-only works, full library collapses → classifier-prior failure** | **YES** |
| correct RT fixes the full library → the classifier cause was an RT consequence | no |

**The extractor, the picker, Gate C, the semi-supervised classifier, the RT
refiner and the FDR estimator all work.** Nine semi-supervised iterations trained
with none skipped. The RT refiner converged properly: residual p50 20.67 → 17.82 s,
p95 52.82 → 48.89 s over two rounds on 6,357 anchors, and pass 2 sized its own
window at 60 s from 2 × p95.

**The blocker is, and has always been, the true-positive PRIOR.** At ~1.5% the
classifier cannot ignite; at ~100% the same code recovers 68%. This is exactly
what review rounds 5 and 6 predicted in advance — *"fixing RT does not fix it"* —
and it is now measured on this pipeline rather than inferred from simulation.

## What this retires

- **"ODIA identifies nothing" is dead as a framing.** It identifies 6,815 when the
  library contains the peptides. Every past zero was a statement about the
  library's true-positive rate, not about the engine.
- **Phase 3 is not the blocker.** The CiRT calibration is good enough — 3.9% slope,
  13 s intercept — and restricting the search to ±60 s changed peak groups
  96,259 → 84,950 without producing a single identification. Calibration is
  necessary, sufficient for its own job, and not the thing standing in the way.
- **The MS1 arm is alive at scale**: 42,446 finite `MS1_COELUTION` values here.

## What it promotes

The priority is now **classifier ignition at a low prior**, which is
`doc/17` §3's Option E and the thing all three reviewers named in round 5 as the
failure mode no plan addressed:

- seed iteration 0 by ranking on the most discriminating feature rather than by
  q-value, with the MS1/MS2 co-elution statistic the obvious candidate — subject
  to its own negative-control gate (task #24), which is still unrun;
- and the orthogonality work, since 18 of 19 sub-scores read the same MS2 traces.

## The ladder is not finished

Only the top rung was run. The dilution series — 100% → 50% → 10% → 5% → 1.5%
present — is what locates the prior at which ignition fails, and that number sizes
the problem. Run it next.

## Caveats

- Astral only. S08 unmeasured.
- The oracle library is DIA-NN's confident set, so it inherits DIA-NN's biases;
  it measures "can the engine find what DIA-NN found", not absolute capability.
  That is the right question for this diagnostic and the wrong one for a
  benchmark.
- Entrapment FDP was not reported in this run and must be checked before 6,815 is
  treated as a quality number rather than a diagnostic one.
- 1:09:48 wall, 2.2 GB peak RSS.

---

## REVIEW ROUND 10 — the conclusion does not survive. The liveness result does.

Codex (effort max) and kimi, both with vault access. **Verdict: red on
"the blocker is the prior". Green on "the pipeline is not dead".**

### What survives

10,040 positives in, 6,815 out at nominal 1% is not something a broken pipeline
produces. The pre-registered gate passed and outcome 1 and 2 (extraction/scoring
bug, q-value bug) are ruled out in the strong sense. That is real and it is the
only claim this document should have made.

### What does not

**1. The dilution series discriminates nothing (kimi, quantitatively).** The vault
already measured the confound in isolation: adding 7,995 entrapment precursors to
a 2,665-target library — **4× size at FIXED prior, no ignition question anywhere**
— took Astral identifications **2,078 → 754 (−64%)**. My first dilution step,
2× size, gave **6,815 → 2,306 (−66%)**. A pure size/competition effect sits almost
exactly on my curve. Add the mechanical effect that the FDR threshold is set by
the extreme tail of the decoy-argmax null, which grows with library size, and the
whole series is **consistent with zero classifier-prior effect.** It is also
consistent with a prior effect. That is the problem.

**2. My "1.5% present" was wrong by 7×.** Measured: the parity library has
**4,991,901 targets, 0.201% of them DIA-NN-confident.** The 1.5% figure came from
the vault's S08 `v6_50k` measurement — a different library — and I carried it
across without checking. The lowest rung (5%) is therefore **25× above** the
production regime; the ladder never reaches the condition it claims to explain.
(Codex flagged an inconsistency here from a different direction, reading library
*rows* as targets; each library is half decoys, so the 50/10/5% labels are
correct. The 1.5% was not.)

**3. The full library is not a rung.** It is a different library build, not the
same 10,040 positives plus absent assays. It should never have been in the same
table.

**4. Outcome 4 was never tested.** The full-library run identified nothing in
pass 1, so its pass 2 was uncalibrated — the converged RT refiner ran only inside
the oracle run. **Nobody applied the oracle-fitted map, or DIA-NN's empirical
RTs, to the full library and re-ran.** "Not an RT consequence" is overreach.

**5. "68% recovery" is a ceiling wearing a capability costume.** The vault:
*"an offline ROC measuring recall of another tool's IDs is a ceiling, not a
result"* — and that reference list has been measured wrong in both directions
(1,230 omitted real IDs, 205 non-IDs). The 32% missed is a mixture of "ODIA
failed" and "DIA-NN was wrong" and this design cannot split it. The honest
wording is **"re-detected comparator-selected precursors"**, not "identified".

**6. The oracle library is the maximally abundance-selected subset of the run**,
and at ~100% prior the classifier trains almost-supervised. That is exactly where
memorisation of target/decoy construction asymmetry is invisible to decoy FDR —
and one such asymmetry is already proven live here: ODIA's decoys carry the
target's precursor m/z.

**7. Entrapment FDP — the one measurement that cuts through all of it — was not
run.** So "the FDR estimator works" is asserted, not shown, and 6,815 is a
diagnostic count.

**8. A pre-registered reporting requirement was dropped.** `doc/16` §1 says
*"report target excess in the top score percentiles, not group counts."* This
document reported ID counts and group counts.

**9. An unexplained 8× asymmetry (kimi).** Oracle: 49,907 groups / 10,040 targets
≈ **5 per target**. Full library: 84,950 groups against an estimated ~150k present
precursors ≈ **0.57 per present target**. "Same binary, only the library changes"
should not do that. Either the picker is library-size-dependent — in which case
"only the library changes" is false in effect — or most present precursors in the
full library **never got a candidate at all**, making the zero partly an
extraction-opportunity result rather than a classifier result. This needs an
answer before any outcome table is believed.

**10. Prevalence is not the classifier prior** (codex). 10,040 targets produced
49,907 peak groups, so the true-event fraction at the classifier input is at most
~20%, not 100%.

### The experiment that actually separates them (kimi)

A **two-arm crossed** design, plus one number per rung that no new code requires:

- **Fixed size, varying prior:** pad to constant N (~200,000 entries) with absent
  entrapment assays; 10,040 / 2,008 / 402 positives → 5% / 1% / 0.2% at identical
  size and identical decoy-null width.
- **Fixed prior, varying size:** subsample the parity library to 100k / 1M / 10M,
  holding the present rate constant.
- **At every rung report target excess in the top score percentiles and the
  iteration-0 seed count** — that single number distinguishes outcome 2
  (scores separate, q does not → competition) from outcome 3 (no separable seed →
  prior) directly.

Plus: entrapment FDP at every rung, and state the key space of every count.

---

## THE PRE-REGISTERED NUMBER, finally measured (`bench/tail_excess.py`)

`doc/16` §1 required *"target excess in the top score percentiles, not group
counts"*, and this document dropped it. Kimi caught that. Measured now, on
existing outputs, no new runs:

| run | IDs | top 0.1% targets/decoys | T/D | **targets above the single best decoy** |
|---|---|---|---|---|
| oracle | 6,815 | 50 / 0 | ∞ | **5,357** |
| dilute_10 | 807 | 18 / 0 | ∞ | **1,268** |
| **full library** | **0** | **81 / 4** | **20.25** | **56** |

**The full-library run HAS a separable tail.** 95.3% target fraction at the top
0.1%, and 56 targets outscore all 43,028 of its decoys — while accepting nothing
at 1% FDR.

**So outcome 3 as I stated it — "no separable seed" — is wrong.** The seed is not
absent; it is small. Seed size tracks the outcome monotonically across all three
runs (5,357 → 6,815 IDs; 1,268 → 807; 56 → 0).

That leaves two live readings, and this measurement does not yet separate them:

- **the seed is small because the prior is low** (few present peptides → few
  strong candidates), or
- **the seed is small because competition destroyed it** — best-of-N over a
  larger library both raises the decoy tail and lets absent-but-lucky targets
  outrank present ones.

**And one anomaly demands an answer before either:** 56 peak-group rows clear the
entire decoy distribution, yet zero precursors are accepted at q ≤ 0.01. q is a
per-precursor quantity computed on the best row per group, so 56 rows may be
fewer distinct precursors — but not zero. Either the per-precursor collapse is
severe, or the q-value path behaves differently at this scale. **This is
outcome 2 territory — "scores separate but q-values do not follow" — which I
recorded as ruled out.** It is not ruled out.

Next: count distinct precursors in that 56, and trace the q-value path for them.
That is a debugging question with a definite answer, not another experiment.

---

## THE MECHANISM, exactly (`lda.h:209-265`)

The q-value estimator is p-value based with the **Käll +1 finite-sample
correction** on the decoy count:

```
    dr  = (decoys_above + 1) / Ndec
    tr  =  targets_above     / Ntar
    q   = pi0 * dr / tr
```

At the extreme tail `decoys_above = 0`, so

```
    q_min  =  pi0 * (Ntar / Ndec) * (1 / targets_above_the_top_decoy)
```

Full-library run: `1 × (41,922 / 43,028) × (1/27)` = **0.0361**, against **0.0365
measured.** The mechanism is exact.

| run | targets above the top decoy | q floor | IDs at 1% |
|---|---|---|---|
| oracle | 5,357 | 0.00019 | 6,815 |
| dilute_10 | 1,268 | 0.0008 | 807 |
| **full library** | **27** | **0.0365** | **0** |

### This is CORRECT behaviour, not a defect

**~100 clean targets is the hard minimum for 1% FDR to be reachable at all**
(q_min ≈ 1/N). With 27, the estimator cannot certify 1% and refuses — which is
what the +1 correction is for. Twenty-seven observations with zero decoys do not
support a 1% claim, and an estimator that returned one would be lying.

So the full-library zero is not an estimator bug, not a q-value path bug, and not
"the classifier never ignited". **The discriminant produced a real but 27-strong
clean seed, and 27 is below the statistical floor for the claim being made.**

### What this changes

- **Outcome 2 is ruled out properly this time** — the q-value path is correct and
  its behaviour is derived from first principles above, not asserted.
- **The target is now a NUMBER, not a direction.** Any intervention — Option E
  seeding, MS1 orthogonality, competition control — has to get
  *targets-above-the-top-decoy* from 27 to ≳100 on the full library. That is a
  cheap thing to measure on existing outputs after any change, with no FDR run
  required.
- **The prior-versus-competition question is unchanged and still unseparated.**
  This explains why 27 yields zero; it does not explain why there are only 27.
  Kimi's crossed design (fixed size / varying prior, fixed prior / varying size)
  remains the experiment that answers it.
