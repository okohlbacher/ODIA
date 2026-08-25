# The candidate is never emitted

2026-08-25, continuing doc/64. The boundary work made the sub-scores able to
tell a correct candidate from a wrong one. This asks what that can be worth,
and the answer changes what the next phase should be.

## Ranking was never the problem, and neither is scoring

For each DIA-NN-confident precursor, take every candidate ODIA actually
emitted and ask whether ANY of them sits within +-3 cycles of DIA-NN's
retention time:

    38,483 DIA-NN-confident precursors that ODIA emitted candidates for
      candidates per precursor: mean 2.87, p50 3, p90 3, max 3

      a candidate exists within +-3 cycles (recall@all):   58.3%
      and it is ranked FIRST:                              51.8%
      recall@1 51.8%   recall@2 55.8%   recall@3 58.3%

Where a correct candidate exists it ranks first **88.9%** of the time. That is
consistent with doc/51's healthy ranking and it closes the question: ranking is
not the deficit, and neither is scoring, because a sub-score can only reorder a
list. **41.7 points are lost to the correct position never being emitted, against
6.5 to ranking.**

Not a calibration artefact. The signed distance of the closest candidate has
median -0.00 cycles, and widening the tolerance from 3 to 10 cycles adds only
6.1 points, to 20 cycles 14.3. The misses are far away, not near.

## Which gate throws it away

The detector's own reject counters cannot answer this -- they count over all
128 positions, and these precursors do receive candidates, just elsewhere. So
the gates were replayed at the true position, and at every other position too,
because the question is one of RANK and rank needs the competitors.

**CORRECTED.** The first version of this table charged a gate whenever the
EXACT DIA-NN apex failed it. That is the wrong question, and a reviewer caught
it: if ODIA's smoothed reference fragment peaks one cycle from DIA-NN's, the
exact position fails `apex_evidence` while a perfectly good candidate is emitted
at k+1. A gate should be charged only when NO position in the +-3-cycle
reference interval survives it. Re-attributed that way the picture changes
substantially, and the wrong version is left here because the size of the error
is the point:

    11,731 confident positives with the true apex inside the window
      positions that become hits: median 39 per precursor, p90 50

                                          exact apex   interval
      <2 fragments present                     0.1%       0.0%
      corr_sum below min_corr_score 0.5        1.7%       0.1%
      shape gates (local max, apex_evidence)  21.3%       2.1%
      outside max_corr_diff 2.0               11.1%       9.3%
      outside the top-3 cap                   15.6%      17.7%
      SURVIVES and is emitted                 49.0%      69.5%

      rank by corr_sum, interval attribution, cumulative:
        <=1 61.0%   <=2 66.2%   <=3 69.5%   <=5 74.2%   <=10 80.4%

The shape gates are barely implicated at all -- 2.1%, not 21.3%. **27 points of
the loss are CAPACITY**: the top-3 cap at 17.7% and the margin at 9.3%. That is
a different diagnosis and it points somewhere different.

Raising `max_candidates` from 3 to 10 would add **10.9 points** of recall on its
own.

## AMENDED THE SAME DAY: the reopening below did not survive its first test

The section that follows argued that "relaxing admission costs 20.4% at matched
FDP" should be reopened, because it was measured when the sub-scores could not
distinguish candidates at all. The argument is still correct about the premise.
It is wrong about the conclusion, and three probes say so.

Selection accuracy -- how often a real sub-score picks the true position out of
the top K -- contains both halves of the trade, the recall a wider gate buys and
the precision it costs:

    cap K      recall    selection accuracy
      3         69.5%           59.9%
      5         74.2%           58.8%
     10         80.4%           57.4%
     20         85.4%           56.5%
     none       87.2%           56.6%

Recall rises 17.7 points and selection accuracy falls monotonically. Relaxing
`apex_evidence` is worse still -- it loses on BOTH axes, because the extra
positions it admits displace the true one from the top three:

    apex_evidence   recall    selection accuracy
        0.99        69.5%          59.9%
        0.95        69.2%          59.9%
        0.90        69.0%          59.8%
        0.80        68.3%          59.2%
        0.00        66.4%          56.9%

And replacing the selection statistic does not help either: ranking the same hit
positions by the co-elution sub-score over the new scoring window gives recall@3
66.2% and accuracy 56.6%, an equal-weight blend 69.0%/59.3%, against corr_sum's
69.5%/59.9%.

So the detector sits at a local optimum. Every threshold moved in either
direction makes it worse, and `apex_evidence = 0.99` is not the crude filter it
looks like -- it is doing real work selecting which positions compete.

**The conclusion that survives is narrower and more useful than the one I
reached first: the 41.7% is not recoverable by tuning this detector.** It needs
evidence the detector does not currently have -- not a different threshold on
the evidence it has.

Narrower still, after review. "The detector is at a local optimum" is more than
the measurements support; what they support is that the CORRELATION-ONLY
proposal policy is near a local Pareto optimum under its present gates AND a
single-sub-score selector. The cap question specifically is NOT closed, because
the selector in those tables is one feature and the pipeline's is nineteen. The
arithmetic: to beat cap 3 in absolute terms the full model needs 80.7%
conditional accuracy at K=5, against the co-elution selector's 79.2% -- 1.5
points. At K=10 it needs 74.5% against 71.4%, 3.1 points. Both are entirely
plausible for a trained model, so a one-feature proxy cannot settle it.

What the corrected table adds: since the loss is capacity rather than gating,
the first thing to do with a fixed capacity is SPEND IT BETTER, which needs no
extra candidates and no FDP risk at all.

Caveat on all three tables: the selector is one sub-score, not the trained
19-feature classifier, which is stronger. These are lower bounds and mechanism
probes. They are consistent enough, and monotone enough, that I would want a
positive reason to expect the classifier to reverse all three.

## Why this reopens a closed question (SUPERSEDED -- see the amendment above)

The plan lists "relaxing admission at all" as closed by measurement: -20.4% at
matched entrapment FDP. That measurement is not wrong, but it was made under a
scoring regime that doc/64 has since shown was degenerate. With window-wide
boundaries a correct candidate and a wrong one received the SAME sub-scores --
median paired difference exactly 0, separation AUC 0.545. Admitting more
candidates into a classifier that cannot distinguish them can only add noise,
so -20.4% is what that experiment had to produce.

The premise has changed: separation is now 0.769/0.782, and 0.78 with the
scoring window split out. Whether relaxed admission still costs was therefore
once again an open question -- and the amendment above answers it: yes, it
still costs.

This is the "verify the negative is real" pattern for the third time in this
project: a committed negative that was measuring a defect elsewhere.

## What this does NOT say

- 49.0% here against 58.3% in production is a difference between a Python
  replay on dumped traces and the shipped C++. The replay is a mechanism probe;
  it locates the gate, it does not size the fix.
- Recall against DIA-NN is agreement with a comparator, not truth.
- No claim is made that raising the cap is free. It costs FDP, and the only
  measurement that settles it is `fdp_compare.py` at matched entrapment FDP.
  The claim was that the reason to believe it costs 20% is gone. Measured
  directly, the cost is still there, so that claim is withdrawn.


## The evidence the detector did not have

`corr_sum` is a statement about SHAPE -- do these fragments rise and fall
together -- and it is the same kind of evidence at every position. That is why
every threshold moved in either direction failed: shape was exhausted. Ranking
by the co-elution sub-score instead of corr_sum was worse for the same reason,
being more of the same.

The library says something shape cannot: which fragments should be BRIGHT. Two
co-eluting species have equally good shape and contradictory relative
intensities, and Pearson over a window is invariant to scale, so corr_sum is
blind to exactly that case. Ranking the margin survivors by normalised corr_sum
plus library correlation, cap unchanged at 3:

    statistic                 recall@3   selection accuracy
    corr_sum                    69.5%          59.9%
    library correlation         71.3%          63.6%
    both                        72.7%          63.3%
    library-weighted corr_sum   71.2%          62.9%

The first change to this detector that improves BOTH, and it admits nothing new
-- it only reorders what corr_sum already accepted, before the cap. Merged.

Label symmetry checked before merging, because the statistic reads library
intensities and a decoy without them would be reordered differently from a
target, which would corrupt FDR silently. `LibraryGenerator.h:243`: a decoy
copies its target's intensities. Symmetric.

## The RT prior is not the next lever

Both comparators use predicted RT at selection time, so it is the obvious next
candidate. Measured, it is too diffuse to help: over 33,333 confident positives
the true apex sits a median of -1.0 cycles from the window centre with a
standard deviation of **23.6 cycles (32.7 s)**.

    window     true apexes kept    positions removed
    +-16            50.5%               74%
    +-24            69.2%               62%
    +-32            82.3%               49%
    +-48            94.9%               24%

Keeping 95% costs +-48 cycles, which removes only a quarter of the search space.
The prior is real but weak, and the trade is poor while 41.7 points are already
being lost. This is also the approved plan's `-rt_window_p95_factor` item: the
measurement says it is worth little on S08.

## Still open

Codex's design point, which the corrected table supports and I have not built:
`apex_evidence` is doing two jobs, centring and non-maximum suppression, and the
non-monotone recall when it is relaxed proves eviction -- relaxing an admission
gate cannot reduce a pre-cap superset, so the extra positions must be consuming
cap slots and displacing the true one. The fix is to cluster proposals into
chromatographic basins and cap GROUPS rather than scan positions. Then a cap of
3 means three distinct peaks rather than three samples of possibly one.
