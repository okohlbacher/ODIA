# The Astral window, and doc/46 settled

2026-08-24, overnight.

## Astral: the window was the cause, and narrowing it doubles identifications

Astral was 8.5x behind DIA-NN with no explanation. The funnel localised it: on
Astral the admission gate rejects **80.2%** of DIA-NN's confident precursors
while rejecting 92.1% of everything that reaches it -- almost no discrimination
-- against 21.6% vs 95.2% on IH1. Gate C's statistic is co-elution evidence,
and Astral's supplied window was 269 s against IH1's 75 s, so every trace
carried 3.6x more retention time to dilute it.

Three arms, one binary, one feature set, so the window is isolated:

    -rt_window_pass1     IDs   entrap     FDP    wall     peak
             269.0 s     180        1   3.23%   21:45   195 GB
             134.5 s     202        0   0.00%   13:23    50 GB
              90.0 s     391        5   7.48%   12:47    49 GB

**180 -> 391 identifications, +117%**, with peak memory down 4x (195 -> 49 GB)
and wall time down 41%. Against DIA-NN's 1,533 that moves Astral from 8.5x
behind to **3.9x** -- the single largest improvement measured on this project,
and it came from a supplied constant rather than from any algorithm.

**The caveat is not small.** Entrapment counts are 1, 0 and 5. FDP is not
measurable at that scale -- the printed 0.00% and 7.48% rest on zero and five
hits -- so this gain is NOT FDP-matched and could in principle be bought with
false positives. What makes that unlikely rather than unknown: the mechanism was
predicted from the gate-discrimination measurement BEFORE the arms ran, and the
memory and wall-time collapse are independent corroboration that the traces
really did get cleaner. Confirming it needs a full Astral run, where the
entrapment population is large enough to measure.

The obvious next question is why the window was 269 s: it is 2 x p95 of the
residual of the map fitted from DIA-NN's own confident set (134.5 s), against
IH1's 37.7 s. So Astral's iRT prediction is genuinely 3.6x worse, and 90 s is
BELOW the p95 -- meaning it now discards true peaks whose prediction is off,
and still wins. That says the interference cost of a wide window exceeds the
recall cost of a narrow one by a wide margin, which is a design rule, not a
tuning result.

## doc/46 is settled: prominence is worse, at matched FDP, on the full run

doc/46 measured `-gate_mode prominence` at **+11.7%** at matched FDP on the full
run. The fixture measured -85%. That disagreement is now resolved on the full
run, with the map supplied so the arm reaches the gate it exists to test
(prom9 was refused at the CiRT seed -- itself the null-inflation mechanism
showing up at a third stage).

    full run, nominal q <= 0.01:   quantile 13,301    prominence 9,280   -30%
    wall 3:26 -> 4:23, peak 138 -> 144 GB

At MATCHED empirical FDP:

    FDP target   quantile   prominence    delta   entrapment   sigma
        5.00%      11,264        4,111   -63.5%      97 / 35   0.99 pp
        5.72%      12,342        7,526   -39.0%     122 / 74   0.84 pp
        7.42%      14,396       11,455   -20.4%    184 / 147   0.82 pp  <- DIA-NN
       10.00%      15,860       12,703   -19.9%    274 / 219   0.91 pp
       15.00%      17,782       14,080   -20.8%    461 / 365   1.05 pp

**-20.4% at DIA-NN's operating point on well-populated cells.** doc/46's
+11.7% does not reproduce and should not be relied on.

### What this says about the fixture

The fixture got the SIGN right and overstated the MAGNITUDE by about 4x
(-85% against -20%). That is the opposite of the direction its guardrails warn
about -- they say FDP effects ATTENUATE ~5x there -- so the guardrail is about
FDP effects specifically and does not transfer to identification effects, which
can be exaggerated instead. Worth adding to the guardrail text: the fixture is
a reliable SIGN detector and an unreliable magnitude estimator, in both
directions depending on what is being measured.

## Where the gate stands now

Both directions are closed. Relaxing admission loses 20% at matched FDP because
it admits decoys as freely as targets. Bypassing admission with an ORACLE --
DIA-NN's own answers, targets only, no decoys -- gains at most 20%. The gate is
doing necessary work and is not the lever.

What the Astral result shows is that the same admission statistic can be made
much more discriminating WITHOUT touching the gate, by cleaning up what it
measures. That is the first thing on this project to move a large gap, and it
suggests the productive direction is the traces the statistic sees, not the
threshold applied to them.
