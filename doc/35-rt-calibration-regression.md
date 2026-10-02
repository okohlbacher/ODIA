# Where the RT calibration went, and putting it back

Status: FIXED and measured on IH1, 2026-08-18. Companion to doc/34, which
measured what the loss cost.

## 1. The seed did not "fail". It was never the method doc/19 measured.

doc/19 §3 specifies a five-step seed, and step 1 is the load-bearing one:

```
1. CiRT BLIND SEARCH     149 precursors, full RT range, no calibration.
2. ROBUST LINE FIT       RANSAC / Theil-Sen on (library iRT, observed apex RT).
3. PREDICT LIBRARY RT    Apply the line to all precursors.
4. EXTRACT IN A WINDOW   +/-60 s rather than blind.
```

What shipped kept the *selection* of CiRT precursors and replaced the *search*
with `PrecursorPrefilter::measure`'s CONTIGUITY statistic -- how many
consecutive cycles a precursor holds a near-complete fragment match. No robust
line fit was ever written: `grep -ri "ransac\|theil"` over `src/` and `include/`
returned nothing. The seed instead used `Calibration::fit`, the binned-median
-> PAVA -> akima MONOTONE fit, which doc/20 measured as the wrong choice at
this anchor count.

So two independent regressions from the measured design, in one code path.

## 2. Why the contiguity statistic cannot work: it saturates

From the IH1 log:

```
rt seed contiguity histogram (target / decoy), 0..9,10+:
  10: 264 / 263
>=3 cycles: 286 / 286  1.0x  p95 602.5 s vs 601.3 s  refused
```

**286 targets and 286 decoys pass at >= 3 cycles -- 81% of each class.** Decoys
are mutated peptides that are not in the sample. When 81% of them show
contiguous near-complete fragment matches, the statistic is counting
coincidence, and doc/19 §1 had already recorded the reason: 5M precursors over
~380-980 Th saturate the m/z axis, so fragment PRESENCE carries no identifying
information at any tolerance.

Both fits therefore described noise (p95 602.5 s vs 601.3 s on a ~1,384 s run)
and the decoy control refused every threshold. **The gate was right.** The
anchors were the defect, and the gate is the only reason a garbage map did not
reach pass 1 wearing a confident face.

What the refusal then did was fall back to spreading the library evenly over the
run -- and doc/34 measured that bill: 111.5 G points, 931 GiB live, 1.465 TiB
RSS, 6 h, `var_rt_delta` dropped as uninformative, and a result the tool itself
labelled a smoke test.

## 3. The fix: search, don't count

`seedRtFromCirtSearch_` runs the REAL extractor, picker and scorer over a
708-precursor sub-library of the standards (`Library::subsetByIndex`), blind
over the whole gradient, and takes each precursor's best candidate by `dscore`.
Co-elution across fragments, not fragment presence.

doc/19's enabling property is smallness, and it holds:

| | full library | CiRT sub-library |
|---|---:|---:|
| precursors | 9,617,705 | 708 |
| wall | 6 h 01 m | **5 m 37 s** |
| peak RSS | 1.465 TiB | **11.5 GiB** |
| live set | 930.97 GiB | 0.15 GiB |

Of those 5m37, **152.8 s is decoding the file** and 2.9 s is matching. The seed
costs one pass over the raw data and essentially nothing else.

## 4. RANSAC, not Theil-Sen -- decided by measurement

doc/19 §3 writes "RANSAC / Theil-Sen" as if interchangeable. They are not here.

doc/20 records how many anchors survive a blind search: **24 inliers on Astral,
47 on IH1**, from ~149 precursors. That is 68-84% wrong -- far past Theil-Sen's
29.3% breakdown point. On a synthetic 50%-contaminated set Theil-Sen returned a
slope 3.7% off with a 221 s p95: it degrades QUIETLY, which is the worst
possible behaviour for a number the seed then gates on. Exhaustive-pairs RANSAC
on the same set returns the slope to 0.24% with a 17.1 s p95.

Two properties of the implementation are load-bearing:

- **Deterministic.** Pairs are enumerated exhaustively (strided above a cap),
  never sampled at random, so the seed does not move between runs of the same
  input -- everything downstream is compared against a previous run's number.
- **A consensus floor against chance.** RANSAC maximises over candidate lines,
  so pure noise still yields the luckiest few points: measured, ~10 of 120
  anchors with a 32 s p95, which would sail through a 10%-of-run residual gate.
  An uncorrelated anchor falls in a band of width `2*tol` within `y_span`, so
  chance supplies `n * 2*tol/y_span` inliers; the fit requires 3x that.
  **The residual alone is not a guard, and the first version of the test assumed
  it was.**

## 5. Measured on IH1

```
CiRT seed: blind search over 708 standards (358 target / 350 decoy control)
  targets: 295 anchors, 49 inliers (17% consensus),
           RT = 453.49 + 1095.00 x libRT, p95 |resid| 23.6 s, median 10.2 s
  decoy control: 291 anchors -- no line survives the consensus floor
CiRT seed ACCEPTED: p95 residual 23.6 s = 1.7% of the run
pass 1 window 70.8 s (3 x the seed's p95)
```

| | old seed | new seed |
|---|---:|---:|
| target p95 residual | 602.5 s | **23.6 s** |
| as % of the run | 43.5% | **1.7%** |
| decoy control | 601.3 s (1.0x) | no fit at all |
| pass 1 window | whole gradient | **70.8 s** |

**25x better, and 49 inliers reproduces doc/20's independently measured 47 for
IH1.** The decoy control now separates completely: 291 anchors, no trend.

## 6. The control gate had the bug the old design's comment warned about

The first version refused this seed. It inherited the prefilter path's rule that
an unfittable control is never a pass -- correct THERE, where a rising
contiguity threshold starved the decoys and "not fittable" meant "sample too
small". Here the decoys are searched blind exactly like the targets, the sample
size is fixed by the standards list, and "not fittable" means RANSAC examined a
full-sized null and found no trend. That is the strongest pass available.

Sample size is now checked first (>= 100 decoy anchors), and only then the fit.
Conflating "could not evaluate" with "evaluated and found nothing" rejected a
seed measured at 1.7% of the run.

## 7. Abort instead of burning hours

`-rt_seed_max_residual_frac` (default 0.10) aborts when the accepted seed's p95
exceeds that fraction of the run's span, and reaching pass 1 with no map at all
now aborts unless `-allow_uncalibrated` is given. The run span is taken from the
apex RTs of every candidate the blind search produced, targets and decoys --
interference is found across the whole gradient, so the extraction estimates its
own run length rather than being told.

This is not conservatism. A seed whose residual is a tenth of the gradient does
not restrict pass 1, so the run pays the unseeded price (doc/34: 15.6x points,
931 GiB, 4 h) for a seeded result. **Failing in 3 minutes strictly dominates
failing in 6 hours.**

## 8. Linear now, monotone later -- unchanged from doc/20

The seed applies a LINE and nothing else. doc/20 measured that a monotone PCHIP
fitted from this many anchors GENERALISES WORSE (Astral 65.5% -> 58.9%, IH1
81.2% -> 76.7% at +/-60 s): ~30-50 anchors cannot constrain a curve, so the fit
chases anchor noise.

The nonlinearity is real -- 135 s of curvature on Astral, 175 s on IH1, mostly
gradient end-effects -- and worth +5.1 pp (Astral) to +12.6 pp (IH1), but only
from thousands of identifications. That is `refineToConvergence_`'s job, and it
is already wired: it runs `Calibration::fit` (binned medians -> PAVA isotonic ->
akima), which is exactly the monotone map doc/20 asks for. It converged in 2
rounds on the 6-hour run; it had simply been refining a map fitted from noise.

**So the three stages the design calls for are now all present and in the right
order:** robust linear anchor from a blind CiRT search -> monotone calibration
from pass-1 identifications -> fine-tuning to convergence.

## 9. Tests

`odia_robust_line` (`robust_line_survives_half_wrong_anchors`) pins both halves:
recovery of a known line under 50% wrong apexes, and outright refusal of noise,
of 3 anchors, and of a degenerate abscissa. 74/74 tests pass.

## 10. Still open

- **The seed's own sub-search reports "identified NOTHING at 1% FDR".** Expected
  -- a few hundred precursors cannot support a target-decoy threshold -- and the
  seed reads only `dscore`, never a q-value. But it means the ranking behind the
  anchors is not itself validated. The decoy control is what stands in for it.
- **17% consensus is low in absolute terms.** It matches doc/20's measurement, so
  it is the method's normal regime rather than a fault, but a seed resting on 49
  anchors deserves the second instrument before anything is concluded from it.
- **Astral is unmeasured.** doc/13's rule is two files before a conclusion is
  recorded. Everything here is IH1 only.
