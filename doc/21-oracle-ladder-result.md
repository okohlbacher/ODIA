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
