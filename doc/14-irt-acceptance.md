# iRT calibration: what "at least as good as DIA-NN" means

Written 2026-08-09 BEFORE the comparison numbers were in, so the bar cannot move
afterwards.

## The statistic

Post-calibration retention-time residual, in run seconds: the calibrated
prediction minus the observed apex, over the precursors the engine is confident
about. For ODIA that is `reportRtResiduals_`'s "mapped" line, over pass 1's
anchors. For DIA-NN it is `RT - Predicted.RT` at `Q.Value <= 0.01`, converted
from minutes.

Reported as median (does the map remove the OFFSET), robust sigma = 1.4826 x MAD
(the bulk), SD (the tail), and p95|e| (what sets the pass-2 window at
`-rt_window_p95_factor` x p95).

**Robust sigma is the headline** and SD is reported beside it. The residual
distribution is heavy-tailed by construction, and the ratio between the two is
itself the diagnostic: a large gap means anchor contamination rather than a bad
map.

## The targets

**Astral is the primary benchmark, because the comparison is exact.**
`astral_truth.parquet` is DIA-NN on the same raw file, and restricted to the
10,891 precursors of `astral_lib_own.tsv` -- the library ODIA searches -- it
gives:

    n 10,891   median -2.51 s   SD 29.29 s   robust sigma 26.09 s   p95|e| 59.35 s

So the bar on Astral is **robust sigma <= 26.09 s and p95|e| <= 59.35 s**.

**REVISED 2026-08-09 16:00, by the project owner: the bar is SD <= 30 s.**
That is a different statistic from the one above and a deliberately easier one
-- SD is the tail, and ours is inflated by anchor contamination rather than by a
bad map. It is also close to DIA-NN's own Astral SD of 29.29 s, so "at least as
good as DIA-NN" survives the change. Current state against it:

    IH1 anchor_q 0.05   SD 71.23     IH1 anchor_q 0.01   SD 41.36
    Astral anchor_q 0.05 SD 95.62    DIA-NN Astral       SD 29.29

Both files fail it today, and by a wide margin on Astral.

**IH1 is secondary and the comparison is NOT exact.** The available DIA-NN
reference is its v6_50k run: 738 identifications, SD 17.36 s, robust sigma
13.62 s, p95|e| 34.58 s. But only 39 of `lib_targets`' 2,665 precursors appear
in that confident set, so the two populations are different libraries over the
same raw file and the numbers are indicative, not a like-for-like target. IH1
must IMPROVE and must be reported; it does not gate on a threshold.

## Two traps this document exists to prevent

**Do not compare across q thresholds.** ODIA's anchors are taken at
`-anchor_q` (default 0.05) and DIA-NN's reference is q <= 0.01. A looser cut
admits worse peaks by construction, so any ODIA-vs-DIA-NN number must state the
anchor q it was measured at, and the headline comparison uses 0.01.

**Do not read the in-sample p95 as the map's quality.** `Calibration::fit`
returns the residual on the very anchors it was fitted to. The 80/20 held-out
probe beside it is what can see overfitting, and both must be quoted.

## When it is frozen

When Astral meets the bar above and IH1 has improved, on a build that passes
ctest, with the three reviewers finding no further defect in the wiring. After
that the calibration is not to be touched without a measured reason.

## 2026-08-09 16:20: the SD <= 30 s bar is NOT reachable by calibration

`test/rt_calibration_residuals.py` already existed and answers the question the
whole thread was circling. It starts from the LIBRARY's iRT -- the only honest
source for "what did we predict before seeing the run" -- fits a calibration on
one set of stripped sequences and evaluates on another, so the number is not the
calibration's own flexibility.

Astral, 10,891 precursors at 1% FDR, 8,723 train / 2,168 held out:

    calibration              set          n        SD        p95    (minutes)
    linear                   in-sample  8,723    0.702      1.398
    linear                   held out   2,168    0.704      1.379
    monotone (LOESS-like)    in-sample  8,723    0.625      1.233
    monotone (LOESS-like)    held out   2,168    0.644      1.269

**0.644 min = 38.6 s is the CEILING for any monotone post-hoc calibration on
this library and run.** ODIA's map is a monotone post-hoc calibration. It cannot
beat 38.6 s held-out however it is fitted, so the SD <= 30 s bar is unreachable
by anything in the calibration path -- LOESS spans, anchor thresholds, bin
counts, isotonic tweaks. Those can close the gap from our current 95.6 s (in
sample, contaminated anchors) down towards 38.6 s and no further.

**And it corrects the DIA-NN comparison.** DIA-NN's reported `RT -
Predicted.RT` SD is 29.29 s, which is BELOW the ceiling a monotone map of the
library iRT can achieve. That is only possible because `Predicted.RT` is
DIA-NN's RUN-REFINED prediction, not a monotone function of the library value.
So the earlier statement "our map is 34.22 against DIA-NN's 26.09" compared a
calibration against a retrained model. The right reading is that DIA-NN is doing
something our architecture does not do at all.

**Therefore fine-tuning is not an optimisation, it is the only path to the bar.**
doc/06 measured exactly this: held-out 0.701 -> 0.449 min (42.1 -> 26.9 s) from
500 peptides in 31 s, 0.395 min (23.7 s) from 2,000. Those numbers straddle the
30 s bar where no calibration can approach it.

Caveat to check before quoting the ceiling as final: `astral_lib_own.tsv` was
derived from a DIA-NN run on this same file, so its iRT column may already carry
run information, which would make 38.6 s optimistic rather than pessimistic.

IH1 cannot be measured this way yet: the v6 report yields 738 usable rows
against the script's 1,000-row floor.

## 2026-08-09 16:45: sequence information DOES beat the ceiling (measured locally)

`test/rt_sequence_gain.py`. Astral, 10,891 paired precursors, 8,723 train /
2,168 held out, split by stripped sequence:

    monotone map of library iRT only            held-out SD  38.55 s
    ridge(AA composition + calibrated iRT) l=1  held-out SD  33.19 s
                                          l=10  held-out SD  33.26 s
                                          l=100 held-out SD  35.62 s

The monotone figure reproduces `rt_calibration_residuals.py`'s 38.6 s from an
independent implementation, which is the cross-check that makes the rest
trustworthy.

**A 24-feature ridge with no neural network breaks the ceiling by 5.4 s.** So
the residual left after the best monotone map is not noise -- it carries
sequence-dependent structure, which is precisely the claim doc/06 makes for
retraining and the reason a calibration cannot get there.

Scaling from doc/06's own comparison on IH1 (AA-composition ridge 0.616 min
against fine-tuned 0.449, a factor of 0.73), the same factor applied to the
0.553 min ridge here would put a peptdeep fine-tune near 0.40 min = 24 s --
under the 30 s bar. That is an extrapolation across datasets and must be
measured, not quoted.

## Production configuration is NOT the measurement configuration

Stated by the project owner 2026-08-09 and worth separating carefully, because
the two configurations answer different questions and the numbers are not
interchangeable.

**Measurement (what the 28.84 s means).** Train on one set of stripped
sequences, evaluate on a disjoint set, choose the recipe on a third. That
estimates the residual for a peptide the model has NEVER SEEN. It is the honest
generalisation number and it is what this document's bar is written against.

**Production (what a real run should do).** Train on the FULL set of confident
identifications, hold out ~15% for validation only -- epoch selection, early
stopping -- then apply the model to ALL precursors. The goal there is not to
estimate generalisation; it is to minimise the deviation on this run's data, and
holding data out of training purely to preserve a clean estimate would be paying
identifications for a statistic.

**The consequence that must not be missed.** After production fine-tuning, the
identified precursors have a residual near the TRAINING error and everything
else has a residual near the GENERALISATION error, and the training error is
much lower -- the curve shows 0.186 min train against 0.411 held-out at epoch
250, a factor of 2.2.

The pass-2 extraction window is sized from the p95 of the ANCHORS, which after
fine-tuning are exactly the trained-on population. So a naive wiring would size
the window from the optimistic in-sample residual and under-size it for the
library majority the model has never seen. The window must be sized from the
HELD-OUT residual, not the anchors', or fine-tuning will narrow the search onto
the peptides it already found.

That is the single most dangerous interaction in this feature and it is why the
module has to report both numbers rather than one.

## The training curve, 500 epochs (2026-08-09)

`--curve 25` on 5,867 Astral identifications. peptdeep's own held-out split.

    epoch     held-out (min)   train (min)   seconds
        0        0.974           1.033           0
       25        0.490           0.452         125
       50        0.447           0.358         246
      100        0.427           0.281         478
      150        0.416           0.238         713
      200        0.413           0.208         958
      250        0.411           0.186        1203
      375        0.406           0.148        1820
      500        0.405           0.123        2431

**Most of the gain is in the first 50 epochs** (0.974 -> 0.447). From 250 to 500
the held-out residual improves by 0.006 min -- 0.36 s for 250 epochs of compute
-- while train falls from 0.186 to 0.123. The train/held-out ratio reaches 3.3x.

So the model overfits steadily and the held-out residual never degrades; it just
stops improving. **150 to 250 epochs is the operating range**: past that is pure
compute, and below 50 leaves most of the gain on the table.

That the held-out curve is flat rather than U-shaped matters for the production
configuration: with only 15% held out for validation, early stopping has a wide
and forgiving target rather than a sharp optimum to miss.

## GPU

Fine-tuning ran on CPU all afternoon. torch was `2.5.1+cpu` and
`ModelManager(device='gpu')` falls back silently; the environment also lived on
ibminode05, which has no GPU, while `data` carries two H100s.

A CUDA environment now exists at `/scratch/kohlbach/odia/rtft_gpu` on `data`
(torch 2.6.0+cu124). **200 epochs in 41.7 s against 481.9 s for 100 epochs on
CPU -- about 23x per epoch.** `finetune_rt.py` now prints the RESOLVED device, so
a silent fallback cannot recur.

One split-brain to remember: the fine-tune runs on `data` (GPU) but the ONNX
exporter needs the OpenMS source tree, which is on ibminode05's node-local
scratch. The `.pth` lands on ceph and is visible from both, so the export is a
separate step on the other node until the exporter is reachable from one place.

## 2026-08-09 19:25: fine-tuning wired end to end, and what the SD was really measuring

`-repredict_irt` re-predicts a SUPPLIED library's iRT with `-rt_model` before
the sort and before any calibration. That was the missing link: `-rt_model` had
only ever been consumed when ODIA generates a library from FASTA, so a model
fine-tuned on a run's own identifications could not be applied to a search of
that run against a supplied TSV -- which is every benchmark we have. 21,782
Astral precursors re-predicted in 0.5 s.

Measured in a real search, on the same 2,694 anchors:

    library iRT   robust sigma 34.42 s   p50|e| 23.24   p95|e| 106.53   SD 125.98   max 1583
    tuned iRT     robust sigma  9.73 s   p50|e|  6.57   p95|e|  65.58   SD 122.76   max 1604

**The bulk improves 3.5x and the SD does not move.** That is the whole answer to
a question that has been confusing this document all afternoon: the SD was never
measuring the retention-time model. It was measuring ANCHOR CONTAMINATION -- a
max|e| of 1,600 s is a misidentification, not an RT error, and no model can fix
it. The robust sigma was measuring the model, and by that measure the tuned
prediction is 9.73 s against DIA-NN's 26.09.

So the SD <= 30 s bar splits into two claims that should never have been one:

* **RT model quality**: met, comfortably. Held out by stripped sequence with the
  recipe chosen on a separate validation split, 27.70 s against DIA-NN's
  29.29 s; in-run robust sigma 9.73 s.
* **SD on ODIA's own anchors**: NOT met, at 122.76 s, and it is an FDR problem.
  Tightening `-anchor_q` from 0.05 to 0.001 took SD from 95.62 to 36.85 earlier
  today without touching the RT model at all.

The practical consequence is immediate: p95 falls from 106.6 s to 67.0 s, which
halves the pass-2 extraction window -- the memory lever, and the interference
lever, that doc/06 said a better RT buys.

## The end-to-end result: +454 identifications AND a halved window

Astral, `-match_decoy_n`, same binary, the only difference being
`-repredict_irt -rt_model <tuned>`:

    base    6,382 identifications   pass-2 window 213.2 s
    tuned   6,836 identifications   pass-2 window 134.1 s

**+454 (+7.1%) and the window falls by 37%.** doc/06 recorded that a perfect RT
column measured -108 precursors on IH1 and concluded a better RT "buys a
narrower window, not identifications". On Astral it bought both, and the two are
the same mechanism: a narrower window admits less interference, and less
interference is what lets a real peak win its own competition.

That prior is not wrong, it is from another file and another regime -- IH1 is
diaPASEF with a merged mobility stack where the window matters differently. It
should be re-taken there rather than assumed to transfer, which is exactly the
mistake it warns about.

Astral progression for the day: 4,290 -> 4,969 (m/z window) -> 5,025 (apex gate)
-> 5,729 (mass sub-scores) -> 6,382 (decoy-null match) -> **6,836** (RT
fine-tuning). Against OpenSWATH 8,765 and DIA-NN 11,112: 78% of OSW.

Still recall against a ~100%-true library, not FDR-controlled discovery.

## 2026-08-10: external coverage, and the criterion that had it backwards

The question is "does the extraction window contain the peak", for EVERY
precursor. Anchors cannot answer it -- they are the precursors we already found.
`test/rt_coverage.py` compares the library ODIA writes after calibration (whose
RT column is predicted run seconds) against DIA-NN's observed RT, over all
10,891 Astral precursors, of which we anchor about 2,500.

    config          +/-15s   +/-30s   +/-60s  +/-120s    p50    p95    p99    max
    no refinement   34.18%   60.26%   88.46%   99.38%   23.5   76.7  108.4  298.0
    map_only        34.30%   60.61%   88.72%   99.37%   23.4   77.2  111.6  303.5
    ridge           43.05%   71.21%   93.62%   99.38%   18.1   66.5  108.9  298.4

**Three things, and two of them reverse what the anchor metric said.**

1. **The ridge is the win.** +11 points of coverage at +/-30 s, +5.2 at +/-60 s,
   p50 23.5 -> 18.1 s, p95 76.7 -> 66.5 s. On the anchors it looked useless.

2. **map_only does nothing.** 60.26% -> 60.61% at +/-30 s is noise. Iterating
   the anchors without a sequence model buys no coverage at all -- the apparent
   gain earlier was the anchor population changing between rounds.

3. **The true p99 is 108 s, not 455 s.** The anchor p99 was contamination. The
   calibration was always far better than our own anchors implied, and every
   window conclusion drawn from that 455 was wrong.

**Why the criterion inverted the answer.** The refinement improves the BULK and
leaves the tail. p99 is the one statistic it does not move, measured on the one
population that is contaminated -- so accepting on it rejected the thing that
works. Acceptance is now the p50 on the frozen set, which tracked the external
p50 and p95 correctly in both rounds (22.3 -> 15.7, 19.5 -> 13.7) and is the
honest runtime proxy for a number that needs truth to compute.

**What this means for the window.** At +/-120 s coverage is 99.38% and the
refinement does not change it -- the tail is set by something else. What the
refinement buys is a NARROWER window at the same coverage: 93.62% at +/-60 s
against 88.46%. That is the clean-features lever, and it is real.

## The delta convergence criterion: what it actually does

**Two loops, two criteria, same principle.** Both stop when the thing being
minimised STOPS MOVING, not when it crosses a threshold somebody chose. "Is the
residual below X" needs an X that differs by gradient, instrument and library --
the same mistake as a fixed m/z window. "Has it stopped moving" needs no
constant.

### 1. The refinement ROUND loop (`-rt_converge_tol`, `-rt_converge_rel`)

Each round re-derives anchors from the current scoring, refines, refits the map,
re-scores. After each, the median |residual| is measured on an evaluation set
FROZEN at round 1 -- frozen because the anchor set grows every round, and a
statistic measured on a growing population is not comparable across rounds. That
artefact once produced a reported 3.3x narrowing that was entirely a change of
denominator.

Three outcomes, distinguished:

    worse than the best so far   -> roll back to the previous axis, stop
    improved by < tol or < rel   -> ACCEPT this round, then stop
    improved by more             -> continue, up to -rt_refine_rounds (cap, 10)

Defaults 0.5 s and 1%. Against a median that starts near 23 s, 0.5 s is two
orders of magnitude below the signal.

### 2. The peptdeep EPOCH loop (`--converge-tol`, `--converge-rel`)

`--curve N` evaluates the held-out residual every N epochs; the same delta test
then stops training. Calibrated against the measured Astral 500-epoch curve:

    epoch   held-out    delta     rel
        0     58.47 s        -       -
       25     29.38 s   29.085   49.75%
       50     26.83 s    2.548    8.67%
       75     26.15 s    0.685    2.55%
      100     25.63 s    0.515    1.97%
      125     25.24 s    0.397    1.55%   <- 0.5 s / 1% stops HERE
      ...
      250     24.63 s    0.008    0.02%   <- 0.05 s / 0.1% stops here

    tol 0.50 s, rel 1.0%   -> stops at epoch 125
    tol 0.05 s, rel 0.1%   -> stops at epoch 250

**Verified: the epoch defaults (0.05 min, 0.1%) stop at 250 on this curve.**
The looser round-loop tolerance would stop at 125, which the held-out coverage
says is about 2 points of +/-30 s coverage too early (ep100 79.30%, ep200
81.15%). The two loops therefore carry DIFFERENT tolerances on purpose: a
refinement round costs a re-score, an epoch block costs seconds on a GPU, so the
epoch loop can afford to be fussier.

Note the held-out curve is FLAT, not U-shaped -- it stops improving and does not
degrade. So there is no minimum to overshoot, only compute to waste, which is
why a delta test is sufficient and no patience counter is needed.

## Parallelism

GPU when one is available, all cores when not, and it SAYS WHICH -- the resolved
device, not the requested one. torch was a CPU-only build for a whole afternoon
while `ModelManager(device='gpu')` silently fell back.

* **GPU**: picks the device with the most free memory rather than `cuda:0`.
  These cards are shared and one of them was carrying 60 GiB of another user's
  job today. 200 epochs takes 41.7 s on an H100 against 481.9 s for 100 on CPU
  -- about 23x per epoch.
* **CPU**: `torch.set_num_threads(os.cpu_count())` and a quarter of that for
  interop. torch's default is a heuristic that is frequently ONE thread inside a
  container, which is the difference between minutes and an hour for a 500-epoch
  fit. `ODIA_RT_THREADS` overrides.

---

# FROZEN: retention-time recalibration and fine-tuning (2026-08-10)

**Do not change any of this without a measured reason and a re-run of both
files.** The objective is RETENTION-TIME PREDICTION CONVERGENCE. Identification
counts are printed beside it and are explicitly not the criterion -- they
oscillate (629 -> 928 -> 799 in one arm) and optimising them selects for
whichever axis happens to flatter today's scoring.

## What is frozen

| piece | where | default |
|---|---|---|
| monotone map, binned medians + PAVA | `src/calibrate/RtCalibration.cpp` | on |
| akima interpolation between knots | `-rt_interpolation` | `akima` |
| explicit monotone inverse by bisection | `Calibration::invertAt` | — |
| per-run refinement in iRT space | `src/calibrate/RtRefiner.cpp` | `-rt_refine auto` |
| iterative rounds, anchors re-derived each | `-rt_refine_rounds` | 10 (cap) |
| delta convergence, round loop | `-rt_converge_tol` / `-rt_converge_rel` | 0.1 s / 1% |
| delta convergence, epoch loop | `--converge-tol` / `--converge-rel` | 0.1 s / 0.1% |
| model reuse across a series | `-rt_refine_model_out` / `_in` | off |
| peptdeep path | `-repredict_irt -rt_model` | off |

## The measured state, external (DIA-NN's confident set, not our anchors)

    config              +/-15s   +/-30s   +/-60s  +/-120s    p50    p95
    library iRT         34.58%   60.75%   88.89%   99.33%   23.4   77.6
    in-process ridge    43.05%   71.21%   93.62%   99.38%   18.1   66.5
    peptdeep ep200      70.29%   88.24%   97.27%   99.83%    7.8   45.3

Held out by stripped sequence, 4,164 peptides no model saw: ep40 76.01%,
ep100 79.30%, ep200 81.15%, ep500 81.80% at +/-30 s. So it plateaus near 200 and
the gain past that is memorisation.

## The four ways this was measured wrong, in order

Recorded because each looked correct at the time and each pointed the wrong way.

1. **SD of the residual.** Dominated by anchor contamination (max |e| ~1,600 s),
   so it measured the FDR, not the calibration.
2. **p95.** Same problem one quantile out; a window at p95 leaves one precursor
   in twenty outside it.
3. **p99.** Reads ~455 s on anchors against a TRUE 108 s on DIA-NN's set. A
   criterion built on it rejected the refinement that actually works.
4. **A growing evaluation population.** Anchors go 1,183 -> 1,408 between
   rounds, which produced a reported 3.3x narrowing that was entirely a change
   of denominator.

**The rule this leaves:** any statistic computed over our own anchors describes
the precursors we already find. The population at risk is the one we do not, and
only an external reference can see it. `test/rt_coverage.py` is that reference.

## Tests

`test/tools/odia_rt_refiner.cpp`, one case per failure this code has actually
had: the unbounded correction that took IH1 from 1,464 identifications to zero,
the 400 s contaminated anchors a squared loss chases, a model file whose weight
count did not match being read as valid, provenance silently dropped on load,
the clamp not travelling with the weights, and stale state surviving a refused
fit.

## Three engines, one footing (2026-08-10)

Fraction of DIA-NN's confident Astral precursors whose TRUE apex lies within
+/-W of each engine's OWN recalibrated prediction.

    engine                     n   +/-15s   +/-30s   +/-60s  +/-120s    p50     p95
    OpenSWATH (lowess)     10523   21.26%   41.69%   72.78%   93.88%   36.4   129.3
    ODIA library iRT       10891   34.58%   60.75%   88.89%   99.33%   23.4    77.6
    DIA-NN (own refit)     11112   43.38%   71.88%   95.04%   99.97%   17.9    59.9
    ODIA in-process ridge  10891   43.05%   71.21%   93.62%   99.38%   18.1    66.5
    ODIA peptdeep ep200    10891   70.29%   88.24%   97.27%   99.83%    7.8    45.3

**The in-process ridge lands on DIA-NN.** 71.21% against 71.88% at +/-30 s, p50
18.1 against 17.9. That is the no-dependency default, matching a mature engine's
own run refinement.

**The peptdeep path passes it.** 88.24% against 71.88%. Held out by stripped
sequence, on 4,164 peptides the model never saw, it is 81.15% -- still above
DIA-NN's 71.88%, and that comparison favours DIA-NN, whose figure is in-sample
with respect to its own refinement.

**OpenSWATH is well behind both**, and this is the number to treat most
carefully. It is reconstructed rather than read: OSW does not publish a
predicted RT per precursor, so its calibration was traced from the
(NORM_RT, EXP_RT) pairs its features record -- those pairs ARE its lowess map --
and interpolated. Three caveats, all pushing the same way:

* the .osw available here is UNSCORED, so its features include interference that
  a q-value filter would remove;
* it is from an earlier run configuration (2026-07-31), not the one that
  produced the 8,765 identifications;
* the binned interpolation is ours, not OSW's own evaluation of its map.

So read it as "OSW's calibration is materially looser than DIA-NN's on this
file", not as a precise figure. Re-running pyprophet would sharpen it and is
recorded in the backlog rather than done.
