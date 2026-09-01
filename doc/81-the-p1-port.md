# 81 — The P1 engine port: +9-11% at matched entrapment FDP

2026-09-01. Commit 6ad5f76. The first shipped identification gain of the
evidence-rebuild program, and the first that survived a full-run matched-FDP
read rather than a cohort projection.

## The chain that produced it

doc/77 (feature plateau) -> doc/79 (oracle: localization is minor) -> doc/80
(tensor parked) -> doc/80 amendment (RETRACTED "traces squeezed dry": every
cohort trace null was computed at ONE summed-trace argmax anchor, and 81.5% of
buried blocks were anchored >5 s off-apex) -> the P0 factorial
(a79_multicand): scoring MANY candidates with an RT prior, correlating each
fragment against an interference-robust reference, and emitting competition
context reached 88.1/75.4 cohort-projected, the best of the program -> this
port.

## What shipped

21 sub-scores appended to PeakGroupScorer: REF_CORR_SUM + REF_CORR_1..12
(per-fragment Pearson against the smoothed best-fragment reference, sorted),
SIG_SHARE_1..6 (sorted area shares), CAND_RANK, CAND_COUNT. No extraction
change, no window change, no picker change; the engine already scored multiple
candidates -- the features simply describe them now. Identical arithmetic for
targets and decoys.

## The measurement (pre-registered in shared/libv2/run_p1.sh before launch)

samelib3 configuration, control = ora_ctl (same lineage, 13,735), binary the
only difference (arm_assert clean), matched entrapment FDP, r = 0.1464:

| FDP | control | port | delta |
|---|---|---|---|
| <=2% | 14,894 (e=43) | 16,233 (e=47) | +1,339 (+9.0%) |
| <=3% | 17,152 (e=75) | 19,012 (e=83) | +1,860 (+10.8%) |
| <=5% | 20,626 | 21,668 | +1,042 |
| <=7.5% | 22,791 | 23,433 | +642 |
| <=10% | 24,572 | 24,880 | +308 |

(The <=1% cell is not a measurement: the control sits at e=2, below the
resolution floor.) Operating-point non-inflation: FDP at each arm's own
q<=0.01 set is 2.09% (control) vs 2.21% (port) -- the same ~2x miscalibration
both sides, no new inflation. Composition of the change at FDP<=2%: +2,462
precursors gained carrying 0.8% entrapment against a 14.6% base rate (real
identifications, not tail junk), -1,123 lost (16 entrapment).

## The honest caveat

DIA-NN-overlap recall is FLAT: 63.8% -> 64.0% at depth 33,330, 37.0% -> 37.2%
at 13,735. The port therefore does NOT close the overlap gap with DIA-NN; it
accepts ~1,300-1,900 additional precursors that DIA-NN does not report, and
entrapment says they are real. The gain is COMPLEMENTARY. Anyone quoting this
result must quote both facts.

## Next

Replication before any default flip: S08 with ODIA's own library, and Astral
neat (reference 7,462), both at matched entrapment FDP. Then the remaining
oracle headroom (the P0 rt_dn arm reached 93.9/87.4 cohort-projected against
the deployable 87.4/73.9) -- that gap lives in candidate SELECTION, i.e.
iterative reselection, the mechanic DIA-NN 1.7.12 runs 12 times and ODIA runs
once.
