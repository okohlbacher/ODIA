# The bottleneck is Gate C, and four attempts to move it failed

2026-08-26. A work programme -- replace the picker with OpenSWATH's, close the
gap, remove the structural defect, add interference removal -- run to its end.
Two of the four changes were built, measured and removed. The diagnosis that
came out of it is worth more than any of them.

## The picker is not the bottleneck. It is not close.

The terminal-reason instrument (`-out_terminal_reasons`) has existed since
doc/51 and had never been run. Against DIA-NN's confident set on the fixture:

    terminal reason      2026-08-23 dump      current build
    scored                  3,859  78.0%        3,822  77.2%
    gate_c                  1,069  21.6%        1,112  22.5%
    few_points                 14   0.3%           12   0.2%
    no_candidate                6   0.1%            2   0.0%

    share of the LOSS:   gate_c 98.2%           gate_c 98.8%
                         no_candidate 0.6%      no_candidate 0.2%

**The picker returns nothing for two precursors.** Gate C rejects and RETURNS
before any picker runs, so 97% of the library never reaches one. Replacing the
picker -- with OpenSWATH's arrangement or any other -- addresses at most 0.6% of
the loss. Both dumps are reported because they disagree slightly and neither is
wrong; the older one is what three code comments cite.

Worth recording about OpenSWATH's picker specifically, since it was the proposed
replacement: it does NOT form groups by co-elution. `MRMTransitionGroupPicker`
picks each transition independently, then greedily consumes the RT axis --
take the tallest remaining peak anywhere, adopt ITS borders as the group's, zero
every peak inside that interval so it cannot seed again, repeat. Co-elution
enters only as an optional post-hoc veto that is FALSE by default. Emission is a
union over transitions and uncapped. Run at ODIA's dimensions it emits 17 noise
groups on a three-peak window, and its Gaussian smoother annihilates edge peaks
outright -- a clean 1000-count peak at cycle 3 of 130 yields zero peaks.

## What was built, and why both were removed

### Gate C library weighting

`coelutionEvidence` gives every transition an equal vote. Weighting by library
intensity looked like the one lever available at the stage that owns the loss:
better admission rather than more of it.

Three blockers, each measured, in the review.

**The volume claim was false as implemented.** At the shipped k = 3.3 the
weighted arm admits 83.9% against 87.3%; matching volume needs k = 2.944, and at
truly matched volume the gain collapses from +1.6 pp recall to **+0.31 pp**.

**The argument fails under its own model.** With realised
`sum(w^2)/contributing` = 1.000 exactly, the equal arm admits 3.20% and the
weighted arm 4.23% -- 32% more. `1.4826 * MAD` is estimated per transition from
finite samples, so each z-score is t-like, not Gaussian. Equal weights average
twelve independent scale errors away; library weights drop the participation
ratio `(sum w^2)^2 / sum w^4` to 3.15 effective transitions, so the sum inherits
a few transitions' scale error and the tail fattens. **A maximum over positions
is a tail statistic, and its tail is governed by the participation ratio, not by
the variance.** Matching a second moment cannot preserve it. The sign is not
even fixed -- it adds volume under Gaussian noise and removes it on real traces.

**It is label-asymmetric, and that is what the arm measured.** Decoys copy their
target's per-fragment library intensities verbatim while their fragment m/z ARE
recomputed. So on a decoy the intensity-to-m/z pairing is scrambled and on an
entrapment target it is a genuine prediction:

    admission at k = 3.3      equal      weighted
    entrapment targets       81.07%       77.62%
    decoys                   82.24%       76.44%

Gate C returns before scoring, so a rejected decoy never enters the null the
q-value is estimated from. q deflates, nominal identifications rise, true FDP
rises -- exactly the +99 IDs at FDP 8.75% -> 9.51% observed. **The
target/decoy-among-admitted ruler that was the change's headline evidence cannot
detect this**, and the vault says so explicitly.

And the whole argument addressed the wrong mode: `gate_mode` defaults to
`quantile`, where `contributing` is never used and the renormalisation is a
per-precursor monotone rescale against a single global tau.

### PROFILE_FIT

DIA-NN's interference idea as a score: the least-interfered fragment carries the
true elution shape, so project the others onto it.

The mechanism does not exist. Plain mean pairwise squared cosine with **no
reference at all** matches the implementation within noise (AUC 0.7733 against
0.7740, 95% CI [-0.0040, +0.0023]), and a *fixed* brightest reference is better
(+1.29 pp at top-5000). A 2x2 over {cosine, Pearson} x {argmax, fixed, none}
puts the entire effect on the metric (+0.084 AUC) and none on the selection
(-0.0007).

So it was never DIA-NN's idea. It was XCORR_SHAPE with the centring removed and
the result squared -- and under this project's own rule that orthogonality
rather than count is the lever, it does not earn a column: R^2 0.761 against the
existing nineteen, residual AUC 0.5500, cross-validated AUC 0.8111 -> 0.8127.

Two defects worth remembering because they will recur: no width guard, so at
width 1 the score is identically 1.0 for every candidate including pure noise
(RT_SPREAD guards precisely this with `width >= 3`); and its NaN fires for the
weakest candidates and is median-imputed at 0.4821, the **69.9th percentile of
the scoreable decoys** -- a candidate that could not be scored outranks 70% of
decoys that could.

## The measurement question, settled: ODIA is deterministic

A reviewer observed that the same default configuration drifted **+117 IDs and
+0.23 pp FDP between 2026-08-23 and 2026-08-26**, larger than the +99 effect
under test, and concluded the arms were not attributable. That was the right
thing to worry about and the wrong conclusion, and the repeat run settles it:

    arm                IDs   entrap    FDP   sigma   DIA-NN    wall    peak
    profilefit_s08   3,540       58   9.62    1.26    3,147   22:41   40419
    replicate_s08    3,540       58   9.62    1.26    3,147   22:47   40385

Identical configuration, identical binary, **identical on every result column**.
Only wall time and peak RSS move, which are scheduling and allocator noise
rather than output. **Run-to-run variance is zero: the classifier's folds, the
extraction and the FDR are all deterministic.**

Two consequences, and they point opposite ways.

**Single runs are exact comparisons.** Two arms differing only in the change
under test differ by exactly that change. So the matched-FDP results above --
Gate C weighting -2.0%, PROFILE_FIT -2.8% at the operating point -- ARE
attributable, and the removals stand on them. There is no noise floor to hide
behind, in either direction.

**The three-day baseline drift was therefore real.** It was not variance; it was
the ~45 commits of boundary work between the two dates changing the output by
+117 IDs and +0.23 pp of FDP. At matched FDP that came to +1.2% at the operating
point, entrapment 40/40 -- small, but a real effect rather than a measurement
artefact.

What the reviewer's caution correctly identifies is a narrower hazard: an arm is
only clean if the binaries differ *by the change alone*. That held for both arms
here. It would not hold for an arm compared against a baseline from a different
build, which is exactly what the `build-gpu` trap made easy before it was fixed.

## What survives

* `bench.sh` executed a binary `build_odia.sh` never writes, and its stale-build
  guard compared against that same wrong file. Both now derive from the build's
  own path.
* `-picker amplitude` was in the valid-strings list and read by nothing, so it
  silently ran the co-elution picker; any arm labelled "amplitude" driven that
  way measured co-elution.
* `-picker openswath` seeded its semi-supervised loop on CORR_SUM, which is
  identically 0 for OpenSWATH candidates. Every recorded number for that arm --
  including the -38.8% on Astral that argues against OpenSWATH-style picking --
  is confounded with this and must be re-measured before it is quoted again.
* Decoy precursor m/z recomputation exists behind a flag, with the 1/K0
  re-derived from CCS at the new mass. It is OFF, and the framing that motivated
  it was wrong: inheriting the precursor m/z is the FIELD CONVENTION, not an
  ODIA defect -- DIA-NN moves fragment m/z only, OpenSWATH's MRMDecoy ships
  `precursor_mz_shift = 0`. doc/24's comparison row said otherwise and has been
  corrected. The consequence is real and shared by every library-based DIA
  engine, which is an argument for co-elution-based MS1 evidence rather than for
  recomputing the mass.

## The rule this round earns

Six claims made in this session's commit messages and comments were wrong, and
every one was wrong in the same direction -- each made a change look better or
better-supported than it was. None was caught by running the pipeline. All six
were caught by reading the code, or the vault, against the claim.
