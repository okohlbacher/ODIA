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

## Why this reopens a closed question

The plan lists "relaxing admission at all" as closed by measurement: -20.4% at
matched entrapment FDP. That measurement is not wrong, but it was made under a
scoring regime that doc/64 has since shown was degenerate. With window-wide
boundaries a correct candidate and a wrong one received the SAME sub-scores --
median paired difference exactly 0, separation AUC 0.545. Admitting more
candidates into a classifier that cannot distinguish them can only add noise,
so -20.4% is what that experiment had to produce.

The premise has changed: separation is now 0.769/0.782, and 0.78 with the
scoring window split out. Whether relaxed admission still costs is once again an
open question, and it is the one worth asking, because it is where 41.7 points
of recall are.

This is the "verify the negative is real" pattern for the third time in this
project: a committed negative that was measuring a defect elsewhere.

## What this does NOT say

- 49.0% here against 58.3% in production is a difference between a Python
  replay on dumped traces and the shipped C++. The replay is a mechanism probe;
  it locates the gate, it does not size the fix.
- Recall against DIA-NN is agreement with a comparator, not truth.
- No claim is made that raising the cap is free. It costs FDP, and the only
  measurement that settles it is `fdp_compare.py` at matched entrapment FDP.
  The claim is narrower: the reason to believe it costs 20% is gone.
