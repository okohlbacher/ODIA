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

    11,731 confident positives with the true apex inside the window
      positions that become hits: median 39 per precursor, p90 50

      <2 fragments present                      0.1%
      corr_sum below min_corr_score 0.5         1.7%
      shape gates (local max, apex_evidence)   21.3%
      outside max_corr_diff 2.0                11.1%
      outside the top-3 cap                    15.6%
      SURVIVES and is emitted                  49.0%

      rank by corr_sum, cumulative over all positives:
        <=1 19.8%   <=2 38.7%   <=3 49.0%   <=5 54.9%   <=10 59.6%

Raising `max_candidates` from 3 to 10 would add **10.5 points** of recall on its
own. The single largest loss is the shape pair: `apex_evidence = 0.99` requires
the reference fragment's smoothed value at the position to be within 1% of its
own maximum over +-3 cycles, on top of being the maximum over +-1.

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
the evidence it has. That is a different phase, and a bigger one.

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
