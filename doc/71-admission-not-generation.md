# Admission, not generation: the correct candidate is present and refused

doc/70 established that ODIA's peak PICKING is sound -- where it accepts a
precursor its apex matches DIA-NN to a median 0.00 s -- and that the gap is
whatever it declines. This locates that gap one step further.

Measured 2026-08-28 on the 10,000-precursor stratified library, full
`S08_diaPASEF.mzpeak`, fixed linear RT map, `-precursor_im_window 0.025`,
`-rt_window 300`, `-passes 1 -gate_alpha 0`, one pinned binary
(`8a48387cca4537eb`), scored against DIA-NN's **reported** RT.

Of the 4,722 DIA-NN-confident precursors ODIA scores, it accepts 2,376 at
q <= 0.01 and declines 2,346. **For 39.3% of those it declines, a candidate
within 10 s of DIA-NN's answer is already in its own candidate list.**

Raising candidate generation confirms which half is binding:

    arm            cand/prec  accepted   delta  apex<=10s  PROPOSED  ADMITTED  conversion
    baseline            2.22     2,376      +0      95.0%     39.3%         -           -
    -max_candidates 10  5.32     2,489    +113      95.3%     58.7%      5.8%        9.9%
    -cand_min_sep 5     2.19     2,441     +65      95.2%     39.6%      4.7%       11.9%
    both                5.18     2,502    +126      95.2%     60.4%      6.6%       11.0%

PROPOSED and ADMITTED are both over the baseline's 2,346 declined precursors:
the fraction for which a correct candidate exists anywhere in the list, and the
fraction that end up accepted at q <= 0.01 with a correct apex. Conversion is
the ratio.

**Generation responds and admission does not.** `-max_candidates 10` lifts the
proposed rate by 19.4 points, from 39.3% to 58.7%. Admission moves from 0 to
5.8%. **Conversion sits at 9.9-11.9% in every arm**: when ODIA has the correct
peak in its hand, it takes it about one time in nine.

That is the deficiency, and it is not candidate generation, not peak selection
(apex accuracy is unchanged at 95.0-95.3% across every arm, so the extra
candidates cost nothing), and not the ranking quality of DScore, which doc/69's
review established beats its own best single input where it counts.

## What the extra candidates are worth anyway

`-max_candidates 10 -candidate_min_separation 5` is +126 accepted precursors on
this library, 2,376 -> 2,502 (+5.3%), with apex accuracy flat. It is a free
modest win and it should be taken, but it addresses a fifth of the gap at most.
The option text at `OpenDIAlyzer.cpp:1024` predicted the generation half of this
("recall of the correct position within the cap rises 72.7% to 75.5%" at
separation 5) and left the classifier half open; the classifier half is where
the loss is.

## Where this points

The remaining question is why a correct, proposed candidate fails to reach
q <= 0.01. Three things it is NOT, all measured: not the apex (0.00 s median
where accepted), not the trace (r = 0.995 against a 0.972 wrong-fragment null,
doc/68 lineage), not candidate availability (39-60% proposed). It is the
admission decision itself -- the score the candidate receives relative to the
decoy null, and the threshold applied to it.

That is the same place the project's earlier accounting landed independently:
Gate C carrying ~98.8% of the emission loss. Two different routes, one target.
