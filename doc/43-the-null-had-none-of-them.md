# The prominence gate failed, and it exposed a defect Gate C was hiding

2026-08-19. Full S08, our 9.6M library, both arms identical but for the gate.

## What the run said

Entrapment FDP is codex's corrected form, `(e/r)/t` over ordinary targets,
not the `(e/(e+t))/r` used in doc/41-42 -- see "A formula correction" below.

                  IDs   entrapment     FDP   est. true   DIA-NN truth recovered
    Gate C     14,081          300  12.63%      12,040                  11,000
    prominence 58,174        6,278  70.20%      17,332                  13,530
    DIA-NN     39,211            -   7.42%     ~36,301                       -

4h14 wall and 156 GB against 2h04 and 80 GB.

## The part that was not a threshold problem

FDP by rank depth, and the same list restricted to precursors the baseline
also reports:

      top N   Gate C   prominence   prominence, shared only
        500     2.3%        61.5%                    12.8%
      1,000     4.6%        53.4%                     7.5%
      5,000     4.8%        44.7%                    11.0%
     14,000    12.1%        59.5%                    14.8%

Prominence's FDP never falls below 20% at ANY depth. Its highest-scoring
precursors are its most contaminated. A gate that is merely too permissive
adds junk at the BOTTOM of a ranking; this added it at the top.

Two hypotheses were pre-registered and both are dead.

**Count / look-elsewhere is dead.** Candidates per precursor are 2.95 and 2.93 --
the extra 2.2M are new precursors getting candidates, not more candidates each.
At MATCHED candidate count the gap survives:

    top-2000 within stratum    n=1      n=2-3
    Gate C                    50.2%       3.8%
    prominence                73.7%      39.8%

Ten-fold worse with the same number of draws. Look-elsewhere depends only on N.

**"The gate selects on the axis the scorer rewards" is also dead** -- my own
claim, refuted by my own data. The newly admitted precursors at the top of the
prominence ranking are not score-alike; they are terrible on every axis:

    prominence top-500     N   entrap%   corr_sum   shape      apex
    new                  278    106.5%      0.629   0.464       951
    shared               222      5.2%      8.854   0.789    26,274

corr_sum 14x LOWER, 28x dimmer, entrapment at the null rate -- i.e. drawn at
chance -- and still ranked first. kimi and codex independently refuted the
claim for a different and better reason: `PeakGroupScorer.cpp:922` computes
`m = coelutionEvidence(...)` ONCE and both modes branch on the same `m`. Gate C
was never orthogonal to prominence. It is the same statistic with an
empirically calibrated heavy-tailed threshold instead of a white-noise one the
code's own option help already disavows ("real interference is heavy-tailed and
structured", "peak-shaped interference passes at any k").

## What was actually wrong

Nothing about the gate. Sorting the run by fragment count:

    usable fragments      2       3       4      12
    var_library_corr   1.0000  1.0000  1.0000  0.6667
    frac exactly 1.0    79.7%   63.8%   51.8%    5.6%
    median DScore      14.389   1.760   0.734  -0.273

The score is monotonically INVERTED in the evidence behind it. `pearson`
(`PeakGroupScorer.cpp:423`) guards only `n < 2`; at n = 2 it returns +-1
always, because two points are always collinear. `fragment_coverage` is
likewise trivially 1.0 there.

And the library is full of them -- asymmetrically:

    fragments        0        1        2        3
    targets      5,285  202,178  158,621  148,227
    decoys           0        0        0  148,227

`appendDecoys` applies the fragment-count bar; `predictFragmentIntensities`
does not re-check it after the MS2 model's intensity floor prunes fragments,
committing `ranked.size()` unconditionally (`LibraryGenerator.cpp:615`). So
366,084 targets (7.3%) sit in a regime containing ZERO decoys --
4,991,901 - 4,625,804 = 366,097 skipped decoys is that same population -- and
they are priced against a null drawn entirely from precursors with strictly
more evidence than they have. The header comment for that parameter names the
failure exactly: "Applying it to one class only is an anti-conservative FDR
(D7 rule 2)."

The two defects compound. A two-fragment target gets a mathematically
guaranteed |r| = 1 AND has no decoy anywhere near it to price that against.
The decoy-based q-value cannot see it, which is why nominal 1% delivered 12.6%
in the baseline and 70.2% under a gate that admitted four times more of them.

Gate C was not controlling this. It was hiding it, by keeping most of that
population out of the scorer.

## The fix and what it buys

`0ba7131`: shrink r^2 by its null expectation 1/(n-1) (exactly 0 at n=2, ~4.5%
at n=12); `-min_library_fragments` default 3, the decoy bar, applied to BOTH
classes; and the generator empties a precursor pruned below the bar.

Post-hoc on the baseline's own output:

                              IDs  entrap     FDP   est. true
          baseline as run  14,081     300  12.63%      12,040
    -min_library_fragments  12,835     132   6.03%      11,937

Half the FDP for 0.9% of the true identifications, and below DIA-NN's 7.42% on
this library. Bars of 4 and 6 measured identical to 3, so the bar belongs
exactly where the decoy builder already put it.

## A formula correction

doc/41 and doc/42 computed entrapment FDP as `(e/(e+t))/r`. The quantity wanted
is the false fraction among reported ordinary targets, `(e/r)/t`
(`OpenDIAlyzer.cpp:3403-3418` describes the right one). Baseline 12.36% ->
12.63%, prominence 62.62% -> 70.20%. Conclusions unchanged; the failure was
worse than reported.

## Still open

- Whether 6.03% holds when the decoy null and classifier are refit ON the
  filtered library rather than filtered after it (`run_full_v3.sh`, running).
- `var_log_sn` is a hard constant log(100) = 4.605 in every population measured
  here, targets and decoys alike. Dead feature, confirmed again.
- kimi: entrapment validity needs a homology check -- Arabidopsis precursors
  sharing a stripped sequence or fragment signature with a human library member
  are not unequivocally absent hypotheses.
- codex: `options.coelution_picking` is set from the legacy `-amplitude_picking`
  flag, not `-picker`, so picker experiments may silently keep the CORR_SUM seed.
- codex's frozen-model cross-over (candidate set x training set, 2x2) is still
  the only clean separation of admission effects from seed poisoning.
- `-gate_log` cannot be joined back to precursors, which blocks the k-sweep.

## VALIDATED, 2026-08-20: `run_full_v3`

The question this run existed to answer was whether 6.03% survives when the
decoy null and the classifier are REFIT on the filtered library rather than the
filtered set being taken from an unfiltered run's output. It does, and better:

                              IDs  entrap     FDP   est. true   truth rec.
          baseline (v2)    14,081     300  12.63%      12,040       11,000
    v3 (min_frag 3 + fix)  13,526     129   5.59%      12,648       11,645
    DIA-NN, same library   39,211       -   7.42%      36,301            -

Post-hoc filtering of v2's own output predicted 12,835 IDs at 6.03% with 11,937
estimated true. The refit beat that on every axis, so retraining on a symmetric
library is worth something beyond removing the rows: 4% fewer raw
identifications, MORE true ones (+5.0%), MORE of DIA-NN's confident set
recovered (+5.9%), and less than half the false discovery proportion. 5.59% is
below DIA-NN's 7.42% on this library.

By fragment count, the degenerate regime is simply gone:

    usable_frag   accepted     FDP
              3          3    0.0%
            4-6         28    0.0%
           7-11      1,302   ~14%
             12     12,193    4.7%

The two-fragment class that supplied 1,214 identifications at 78.4% FDP no
longer exists, and the twelve-fragment core grew 11,699 -> 12,193 while getting
cleaner, 5.0% -> 4.7%.

**Cost, and where it is NOT from.** 3h12 wall and 116 GB peak, against 2h04 and
80 GB. That is not the fragment floor -- the library it searches is SMALLER
(9,251,621 against 9,617,705). It is `9a496fa`, the ms1_coelution repair: this
is the first full run in which MS1 traces are actually built, 1,343 bins over
9.25M precursors, 8,293,585 with signal (89.6%), 47,397 MiB. The baseline's
`var_ms1_coelution` was dead, so it paid neither the memory nor the time.
Whether that feature earns 47 GB and an hour is a separate question, and it has
not been asked yet: its AUC in the v3 run has not been measured.
