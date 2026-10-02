> **SUPERSEDED BY doc/26-PLAN.md (2026-08-15).** Kept as the working record. Where this document and doc/26 disagree, doc/26 wins — it carries the corrections this one predates.

# Phases 3+4: integrated RT recalibration, fine-tuning and m/z calibration

Status: PLAN, not implemented. Written 2026-08-14 after adversarial review by
Codex `gpt-5.6-sol` at `model_reasoning_effort=max`, Kimi 0.36.0 and Vibe 2.24.1.
Brief and raw reviews: `/scratch/kohlbach/odia/massanchor/bench/review_p34/`.

**This is the first review round run at codex effort `max`.** Every prior round
in this project ran at the config default `low` without my noticing. The
difference is visible in the output and the strongest single argument below
(§2.2) comes from that arm.

---

## 0. The finding that reorganised the plan

I framed the 0-identification failure as *one* problem (RT) and then as *two
independent* problems (RT + classifier ignition). **Codex refuted the
independence claim and is right:**

> The "second independent cause" is not yet independent unless the classifier
> has already been tested under an oracle or externally supplied RT map.

Nobody has ever run ODIA's scorer with known-correct RTs on a known-present
population. Until that is done, "the classifier cannot ignite at a 1.5% prior"
is an inference from simulation, not a measurement on this pipeline. It may also
be masking something duller and worse — a wrong transition-to-window
association, a score sign, NaN/tie handling, broken target-decoy competition, or
non-exchangeable decoys.

**Consequence: the first step is not a fix, it is a diagnostic ladder (§1).**

A related warning, also from Codex, that retires a number I have quoted
repeatedly: **a 1:1 target:decoy peak-group count is nearly uninformative.** It
is consistent with simply generating one group per target and one per decoy.
What matters is target excess *in the extreme score tail*, and exchangeability.
Stop citing 819,353 / 819,835 as evidence of anything.

---

## 1. Step 0 — the oracle ladder (do this first, it is cheap)

Run the scorer on a library of **DIA-NN's confident precursors only**, with
matched decoys, their **known RTs**, and generous fixed windows. Then dilute with
absent/entrapment assays — 100%, 50%, 10%, 5%, 1.5% present — down to the
production prior.

This separates four outcomes that the current evidence cannot distinguish:

| outcome | diagnosis |
|---|---|
| fails even on the positive-only library | extraction / scoring / FDR **bug** |
| scores separate but q-values do not | target-decoy competition or q-value bug |
| positive-only works, full library collapses | genuine classifier-prior failure |
| correct external RT fixes the full library | the "independent" classifier cause was an RT consequence |

**Do not proceed to any calibration work until this passes.** If known positives
at known RT and m/z do not dominate decoys, every option below is irrelevant.

Report target excess in the top score percentiles, not group counts.

---

## 2. The four reviewed options, after review

### 2.1 Consensus, 3-0

- **Only A and B break the RT circle.** C and D both require pass 1 to have
  succeeded. My reading was correct on this point.
- **C creates no information.** Batching re-chunks the same 98.5%-noise
  population; the iteration-0 seed is as non-separable in batch 1 as in the whole
  run. Already measured and refuted (`BACKLOG.md:1437`). Both Codex and Kimi add
  the same structural point independently: **DIA-NN's batching works because
  DIA-NN always has predicted RTs in hand.** C copies the scheduling without the
  precondition. C without B is cargo-cult.
- **D's central assumption is false.** I wrote that the feedback from a narrower
  window onto co-elution scores is second-order. Codex: at 100 → 10 ppm the
  accidental fragment-hit density changes by ~an order of magnitude, and at a
  1.5% prior changes in the extreme score tail are **first-order**. "No second
  disk decode" is defensible; "no re-extraction, re-peak-picking, re-scoring" is
  not. Retained residuals and apexes alone are insufficient.
- **All three named the same missing failure mode** — see §3.

### 2.2 B is the primary circle-breaker, and it has a hidden dependency

**Both Codex and Kimi caught something I had not seen.** peptdeep emits values in
**iRT space**. A predicted iRT is only a rank until `a` and `b` place it in run
seconds:

```
    RT_hat = a + b * iRT_hat
```

B breaks the circle **only if `a` and `b` come from something available before
any ODIA identification** — a frozen LC-method transform, gradient metadata,
spiked standards, a transferred run, or a predictor already emitting
method-specific seconds. If they were fitted from the anchors that do not exist,
B contains the same circularity it claims to remove.

Worse, the frozen coverage table (±30 s 88.24%, ±60 s 97.27%, p50 7.8 s) was
measured **against DIA-NN's set**, i.e. through an alignment that will not exist
in production. **Those numbers do not validate B as written.** They must be
re-measured through whatever exogenous transform is actually deployed.

**Kimi's mitigation is the reason B still survives, and it is checkable from the
frozen table itself:** plain library iRT, with no anchors at all, already gives
±30 s 60.75% and **±60 s 88.89%**. So B degrades gracefully — even a degenerate
alignment plus a fixed ±60 s window covers ~89% of the eventual set, and
peptdeep lifts that to 97.27%. Take the **±60 s** window, not ±30 s.

Codex's arithmetic on the bootstrap population, assuming the transform is
genuinely exogenous: `738 × 0.9727 ≈ 718` recoverable on Astral — a strong seed
population, and notably **still below `MinCal = 1000`**, which is a further
argument against C.

The "biased population" attack on B is **toothless**, and Kimi argues this
directly against the position Vibe and I both leaned toward: absent peptides
never become IDs or anchors regardless of prediction quality, so a wrong
predicted RT for them costs one wasted extraction and nothing else. What matters
is coverage on *present* precursors. Do not spend effort measuring peptdeep on
the absent 98.5%.

### 2.3 A — the saturation question splits the panel 1-2, and Codex wins the argument

This was the question I flagged as decisive, and it is the one place the
reviewers genuinely disagree.

**Kimi (pro-A):** the saturation finding indicts *max-over-frame depth*, not
gated temporal coincidence. Frame depth fails because it is `1-(1-p)^10000 ≈ 1`.
Single-spectrum joint presence is `~p`; contiguity-3 is `~p³` per starting cycle.
So contiguity converts a **spatial** extreme-value statistic into a **temporal
coherence** statistic — and temporal coherence is exactly what measured as
discriminating on this data (MS1/MS2 co-elution, 13.7×).

**Codex (anti-A), and this is the strongest technical point in the three
reviews:** the step `p₃ = p³` **assumes independence across cycles, and DIA
cycles are not independent.** Real interfering ions persist chromatographically,
the same isolation window recurs, and fragment matches caused by *another*
peptide persist for several cycles. So contiguity-3 rejects isolated electronic
noise while **accepting chromatographic interference — precisely the nuisance
that matters in DIA.** The algorithm also still maximises over all triplets.

Codex's numerical burden, which sets the bar: with 738 positives against ~49,262
effectively-null assays, 90% seed precision at 80% recall requires an
absent-assay seed rate below **~0.14%**; 99% purity needs **~0.012%**. Going from
49,999/50,000 to that is several orders of magnitude of null suppression.

**Resolution — run Kimi's experiment, judge it by Codex's criterion.** Kimi's
sweep is the right instrument and costs one afternoon and one flag: run
`-rt_seed_min_contiguity` = 1, 2, 3, 4 and plot seed count. But **plateau shape
alone is not the acceptance test** — the test is the absolute false-seed rate
against a null that preserves m/z density and isolation-window membership
(entrapment or fragment-shuffled assays), against the 0.14% / 0.012% bar.

Do **not** validate A using residuals from A's own selected seeds. That is the
anchor-statistics mistake that already cost this project the most.

**A's role, if it survives, is not B's competitor — it is B's alignment source.**
Search for prefilter evidence *inside B's predicted window*, and use it to
estimate a small run-level shift or slope. That removes most of A's whole-run
extreme-value burden, and it de-oracles the `a`/`b` transform that §2.2 shows B
needs. Both Codex and Kimi reached this same subordinate role independently.

---

## 3. Option E — the missing piece all three reviewers named

None of A–D addresses the **iteration-0 classifier seed**. The semi-supervised
LDA selects positives by q at iteration 0; if that initial ranking is not
enriched, the loop never ignites, and fixing RT does not fix it.

- Codex: "initial-score / FDR identifiability under a massive null library."
- Kimi: seed Percolator-style, by **ranking on the single most discriminating
  feature** — measured to be MS1/MS2 co-elution correlation, 13.7× top bin — with
  no q-values involved.
- Vibe: "the LDA ignition problem at ultra-low prior."

Codex's conversion: 13.7× against a 1.5% base is **~17% positive fraction**. Not
enough for 1% FDR by itself, but plausibly enough to *initialise a classifier*.

Requirements, from Codex: apply the same RT-window and search-opportunity rules
to paired targets and decoys, and recompute q-values **globally** after the
classifier is fixed. Do not retain successful batches with batch-local
q-values — that invalidates global FDR.

**E is orthogonal to A–D and may be the actual blocker.** Without it we may fix
the RT circle and still measure 0.

---

## 4. The objective function

Both rejected statistics fail for the same reason: they aggregate over the wrong
decomposition. IDs are a thresholded count with huge boundary variance
(629 → 928 → 799) and selection bias by construction. Total scatter mixes
systematic error with irreducible noise — hence 4.66 → 4.19 ppm while 68% of
systematic error vanished.

**Separate the fitting metric, the stopping metric and the success metric.**

| | metric |
|---|---|
| RT fit | robust Huber / quantile loss on fixed training calibrants |
| RT validation | p50/p90/p95 abs error, ±30 s / ±60 s coverage, **and worst-bin coverage across the gradient** |
| mass fit | robust conditional-mode regression, or narrow-signal-plus-uniform-background mixture |
| mass validation | residual conditional centre vs log m/z, **max absolute bin bias**, fixed-window transition coverage |
| **stopping** | **parameter stability** — refit on data halves / bootstrap resamples; converged when map coefficients and mass terms move by less than the window tolerance |
| success | independently checked FDP ≤ 1%, reproducibility, ID yield **after the configuration is frozen** |

Kimi's stability criterion is preferred as the stopping rule because it is immune
to *both* measured failure modes: it does not oscillate like IDs and it is not
blind like total scatter. Codex's mass-side form is equivalent and matches the
language already in use: **converged when no held-out covariate term retains
|t| > 3.** The existing t = 7.4 (shape real) and t = 0.26 (drift absent) already
*are* this statistic.

On IDs, the correct statement is narrow and worth quoting exactly:

> Do not optimise calibration iterations by IDs; do require the frozen calibrated
> pipeline to improve or preserve valid IDs.

**The 85/15 held-out split of `doc/15` §6.2 is necessary but insufficient**, and
Codex's objection is correct: if both folds consist of q-accepted groups, both
are conditioned on the same easy subset, reproducing the exact selection failure
that took Astral from 4,969 to 2,422. Define the split **before scoring**,
stratify across RT / m/z / charge / intensity, and keep a separate final test set
if the 15% fold is used repeatedly for stopping.

---

## 5. What the research vault adds — read AFTER §1-4, it changes two of them

`vault/` (79 notes, sibling of the repo, never `git add`ed) was under-used when
§1-4 were written. Its governing caveat: the "ODIA" in seeded notes is the
**previous attempt**. Claims about DIA-NN, OpenSWATH, the literature and
measurement discipline transfer; claims about that code do not. The two findings
below are properties of DIA data and of the FDR estimator, so they transfer.

### 5.1 Wider is not better — `[[Library prefiltering and its ceiling]]`

Measured end to end on Astral:

| configuration | search space | peptides |
|---|---:|---:|
| fixed rule | 423,129 | **6,024** |
| learned model @0.06 | 707,311 | 5,855 |
| learned model @0.20 | 1,960,309 | **5,705** |

Identifications fall **monotonically** as the budget grows, while the offline
recall ceiling *rises* 85.1% → 93.8%. Admitted candidates enlarge the
target-decoy null faster than they contribute true positives.

> **A better recall ceiling is not a better result.**

**This amends §2.2.** Kimi's "take ±60 s, not ±30 s" — which I accepted — treats
coverage as the objective. It is not: window width trades coverage against null
width and therefore has an **optimum**, not a maximum. Sweep it; do not set it
from the coverage table.

Related caution from the same note: at its usual operating point the depth rule
is barely better than a coin flip (depth ≥ 4: T/D **1.103**, 52.4% target), and
essentially all discrimination lives at full depth (depth ≥ 6: T/D **27.561**,
96.5% target). **This partially rescues Option A** against §2.3 — fragment depth
*does* discriminate — but only at a far deeper operating point than the current
default implies, and depth is a different axis from contiguity. Sweep both.

### 5.2 The winner's curse explains the ~40 groups per precursor

`[[Target-decoy competition and the winner's curse]]`: scoring several candidate
peak groups per precursor and keeping the **maximum** means the decoy null is the
null of "best of N", which is **wider** than the null of one draw. The estimator
stays unbiased *provided opportunity distributions are exchangeable between
targets and decoys* — but a wider null means a higher threshold and fewer
identifications for everyone.

That answers §6's open question, and it supplies a diagnostic that needs **no
calibration at all**:

> **Equalise opportunities** — one prespecified RT-nearest peak, or the same peak
> count for every label. Strong contraction of the decoy tail confirms winner's
> curse as a contributor.

Cheap, decisive, never run on this codebase. It is promoted into the order below.

### 5.3 m/z is the cheap selectivity lever; RT is the expensive one

From `[[M2 prefilter selectivity measured]]`: DIA-NN recommended **4.6 ppm** MS1
accuracy on a run it had calibrated at 22 ppm — *"extraction tolerance is worth
roughly another order of magnitude of selectivity where the calibration supports
it, and unlike RT it is cheap to obtain."*

**This inverts the §5 ordering as first drafted**, which put mass calibration
last. Mass moves early — which is also the genuine phase-3+4 integration.

The same note carries a calibration of expectations that must be settled before
±60 s becomes a default: **DIA-NN's RT window on an in-silico predicted library
is ±17.2 MINUTES (~±1030 s)**, ~11% of the gradient — 17× wider than §2.2
proposes. Either peptdeep is far better than DIA-NN's library RT, or ±60 s is
over-tight. Given §2.2's alignment caveat, do not assume the former.

And its conclusion, stated for the previous codebase but structurally relevant:
**the prefilter is only viable downstream of a calibration pass, never on a cold
predicted library**, because cell width follows library RT accuracy. Option A
proposes exactly the cold-library regime. Weigh that against §2.3.

---

## 6. The order of work

```
  0. ORACLE LADDER            positive-only -> diluted to 1.5%     [gate: passes]
  0.5 EQUALISE OPPORTUNITIES  one RT-nearest peak per label;       [gate: tail
      (winner's curse, 5.2)   does the decoy tail contract?         contracts]
  1. MASS FIRST (5.3)         loose -> fit log(m/z) -> tighten
      -> centre only. Do NOT infer WIDTH from accepted IDs (4,969 -> 2,422)
      -> no f2(RT) term on IH1: t = 0.26
  2. B: predicted iRT         exogenous a,b; window SWEPT not      [gate: IDs > 0]
                              maximised (5.1)
  3. E: co-elution seed       rank iteration-0 by MS1/MS2 corr     [gate: ignition]
  4. A: depth x contiguity    sweep BOTH; judge on absolute false-seed rate
                              against an m/z-preserving null       [gate: <0.14%]
     -> if it survives, A refines B's alignment; it does not replace it
  5. CONFIRMATORY EXTRACTION  exact, from a raw-centroid cache
```

Steps 0, 0.5 and 1 need no RT calibration and no classifier ignition, so they are
independent of the blocker and can start immediately. Steps 2-3 are the blocker
proper.

Steps 0–2 are the blocker. Steps 4–5 are `doc/15`, which stands, with two
amendments: its inner loop is demoted (§2.1) and its held-out split is
strengthened (§4).

**What is deliberately NOT in this plan:** joint RT×mass iteration to
convergence. The coupling is measured absent on IH1 (−0.11 ppm over 1,600 s,
t = 0.26, confirmed over 124 M hits), a 101-agent survey found no engine that
does it, and DIA-NN calls its own auto-tuning "inherently noisy" and recommends
fixing the parameters. `doc/15` §6 already said this; the review confirms it.

---

## 6. Open questions this plan does not answer

- **~40 peak groups per present precursor.** 1,639,188 groups against ~37,600
  present precursors. Nobody has explained that ratio; it inflates the
  multiple-testing burden the q-values must overcome, and no option here reduces
  it. Kimi's warning: a wide fixed RT window (B) makes it **worse** before it
  makes it better. Watch it explicitly in step 1.
- **What is the acceptance bar?** At a 1.5% prior, 1% FDR needs a likelihood
  ratio ~6500:1 at the boundary. DIA-NN's 738/50,000 proves it is reachable, but
  the margin is thin. Decide now whether "1% FDR on the full library" is the test
  or whether a sample-matched library is the realistic bar — before measuring,
  not after.
- **B is Astral-only.** The frozen RT table has no IH1 arm.
