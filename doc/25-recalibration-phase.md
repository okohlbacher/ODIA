# The full recalibration phase

Plan, 2026-08-15, from review round 12 (codex effort max, kimi, both with vault).
FDR/FDP is explicitly out of scope — it is the last phase.

## 0. What the review corrected in my diagnosis

**Three confounds in the candidate-flow analysis (kimi):**

1. **The axis in force is our worst predictor.** The CiRT line that ran is *worse
   than raw library iRT* — p50 42.0 s vs 36.7 s, decile bias to −102 s. So
   "availability is window-limited" is true only against the worst axis. The
   frozen table moves ±60 s coverage 88.89% → 93.62% → 97.27% across predictor
   levels **with no window change at all**. **The window deficit is largely a
   predictor deficit.**
2. **Truth noise.** DIA-NN's own FDP >2.3% at nominal 1% means ~230 of the 10,040
   "truths" have effectively arbitrary RTs, landing almost entirely in the
   outside-window bucket. Tail numbers to three significant figures are
   overprecise.
3. **Every candidate-flow number was measured with mass calibration half
   applied.** I diagnosed the m/z bug and then drew conclusions from data the
   buggy pipeline produced.

**One number of mine is uninterpretable as printed.** "Availability at ±30 s =
75.8%" — ±30 s of *what*? I measured ±30 s of the **true apex**; the reachability
sweep measured ±60 s of the **prediction**. Different axes, presented as
comparable. Availability cannot exceed reachability on the same axis, and mine
did — which is the tell. **Every residual number must name its axis.**

**And a discrepancy to reconcile before designing against either:** my raw library
iRT gives p50 36.7 s / 65.5% at ±60 s (n=10,040); the frozen note gives p50 23.4 s
/ 88.89% (n=10,891) for "ODIA library iRT" on the same instrument. Either the
frozen row already includes the monotone map, or one of the two is wrong.

**Codex's flow correction** (one denominator, mutually exclusive) — the loss is
two roughly equal halves, not one:

| bucket (of 10,040) | n | % |
|---|---|---|
| no candidate at all | 22 | 0.2% |
| truth outside the extraction support | 1,172 | 11.7% |
| reachable, but no candidate within 30 s of truth | 1,275 | 12.7% |
| candidate at truth, not ranked first | 171 | 1.7% |
| correct and ranked first | 7,400 | 73.7% |

Effective support measured from the data: **±112 s**, not the ±60 assumed.
And 54% of misses sitting exactly at K=3 is **evidence of censoring, not against
it** — I had that backwards. Only a K sweep on identical chromatograms settles it.

---

## 1. Step 0 — PLUMBING, before any design

Both reviewers: nothing else matters while measurements do not reach pass 2.

- pass 2 must consume the **measured mass model** (pass 1 measured 4.13 ppm;
  pass 2 re-extracted at 7.69 ppm, byte-identical)
- pass 2 must consume the **derived RT window** (128 s measured, 60 s used,
  silently capped)
- kimi: *"three measured quantities discarded by downstream code is one bug, not
  three"*

## 2. The order

```
  0  PLUMBING            measurements must reach pass 2
  1  m/z RECALIBRATION   first, once, NO loop
  2  ROBUST LINEAR SEED  RANSAC/jackknife, not OLS on contaminated anchors
  3  MONOTONE MAP        binned medians + PAVA + akima, on seed residuals
  4  FINE-TUNE the       AFTER the map: fine-tune labels are observed RTs
     sequence model      expressed through the current map, so fitting against a
                         stale axis bakes it in
  5  REFIT the map       mandatory and cheap -- peptdeep emits iRT, a RANK not
     on the new axis     seconds, so the model does not make the map redundant
  6  WINDOW + ONE final extraction
```

**Why m/z first, and it is synergistic not merely independent:** mass `f2(RT)` is
measured absent on S08 (drift −0.11 ppm at t = 0.26 over 124 M hits); DIA-NN 1.7.x
defaults to a **single RT bin** (`MassCalBinsMax = 1`), agreeing that mass needs no
RT term; and **tighter ppm raises anchor uniqueness** (365 unique pairs at 5 ppm vs
196 at 20), so mass calibration *expands and purifies the RT anchor set*.

**Never carried across a predictor change:** the window, and any residual
statistic. **May be carried:** the mass model — re-verify per instrument with the
f2 t-test rather than refitting by reflex.

## 3. Deriving the window without truth

Both current options rejected: `p95 × 2` with a cap is arbitrary *and* the cap
discarded the measurement; high anchor quantiles are disqualified because anchor
p99 ≈ 455 s against a true 108 s — they read **misassignment**, not calibration.

1. **Fit the core, not the tail.** Anchor residuals are a mixture: true
   assignments in a tight core, misassignments broad and roughly uniform over the
   searched region. Estimate the core scale robustly — MAD about the median, or
   the p25–p75 IQR (contamination-insensitive where p10–p90 is not). Window =
   k × core-scale, k chosen **offline** on the *final* predictor.
   Reference point: DIA-NN sizes its RT window from the **80th-percentile**
   residual — a mid quantile, deliberately below the contamination tail.
2. **Respect censoring.** Residuals are only observable inside the window in
   force, so a window derived from a ±60 s-truncated sample cannot see past 60 s.
   **The pass-1 harvest at the 1e9 s window is the right input — uncensored.**
   Pass-2 anchors are not.
3. **Edge pile-up diagnostic** — fraction of picked apexes within ε of the window
   boundary. Computable at runtime with no truth; an excess means the window is
   binding *on this run, right now*, whatever the anchors say.
4. **Heteroscedastic.** Error varies systematically along iRT (our own decile
   table). The map already bins — carry the per-bin residual scale out and make
   the window per-bin. A global window over-covers easy deciles, paying a
   winner's-curse tax, and under-covers hard ones.
5. **Validate on precursors the feature finder did NOT detect** — the only test
   that speaks for the population at risk.

**And with the fine-tuned axis (p95 = 45.3 s), ±90 s already buys ~99% core
coverage — the whole ±150/±200 discussion becomes moot.** That is the strongest
argument for doing step 4 before sizing anything.

## 4. Convergence

The four statistics that pointed the wrong way during the original RT work shared
one property: **each was a function of a changing, self-selected population.** SD,
p95 and p99 moved because the anchor set moved (1,183 → 1,408), and the anchor set
moved because the model changed.

**Rule: convergence quantities must be functions of model PARAMETERS, or of
predictions on a set frozen before round 1 — never of the fitting sample.**

- Freeze a disjoint, stratified evaluation set before any fitting. Never let it
  grow.
- Extend the RtRefiner's existing parameter-delta criterion (0.1 s, cap 10) to the
  map (max change in predicted RT over the frozen set) and the mass model
  (coefficient deltas **relative to their standard errors** — the 0.1 ppm
  tolerance was recorded as *below* the MAD's standard error, i.e. statistically
  invisible).

## 5. Lossless final extraction

A single narrow window is never right. **Narrow primary path + fixed wide rescue
path**, with unsupported strata keeping wide windows — codex's shape, and it
matches the measured heteroscedasticity.

## 6. Open, before implementing beyond step 0

- Reconcile the 65.5% vs 88.89% discrepancy.
- Re-measure the reachability curve on the **curve** and **fine-tuned** axes; the
  one I have is the raw-axis curve and cannot size anything.
- The K = 3/5/10/all replay on identical chromatograms — the only experiment that
  splits the 12.7% bucket into cap-censoring vs picker failure.
- Audit **why the current CiRT seed degrades the axis** relative to raw library
  iRT. A seed step that makes things worse poisons everything hung off it.
