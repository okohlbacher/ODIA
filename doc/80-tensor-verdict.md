# 80 — The multi-channel tensor: parked by its own guards, which is the result

2026-08-31 evening. The full program of doc/78's L3 (design: analysis77/
TENSOR_PLAN.md + round-2 addendum; 30-run matrix on spock H100; frozen eval;
one-shot seal). Everything below was pre-registered before training.

## What happened

**Contrast (a) — decoy-supervised — is anti-conservative poison, and every
guard caught it.** Appending its score DROPS the strongest stack 87.3 -> 80.8
(-6.5, every seed negative). The seal's paired-bootstrap read shows why:
stacked candidates score decoys to near-zero (AUC decoy-vs-null 0.003-0.011
vs R0 0.4855) — a scorer that would silently detonate any decoy-based FDR.
The shadow tripwire fired (dsym -0.078); the gate battery REFUSED the
interference-recognition claim (partial Spearman ~0: the gate keys on
intensity, not interference); trace attribution FL5-FL4 = -7.6 (the deep
channels HURT under decoy supervision). The prior era's entrapment failure
(docs 48-63, parked 2026-08-25) would have been repeated exactly — this time
it was caught before a single full-run token was spent.

**Contrast (b) — DN vs non-DN targets, decoy-free — is real but under the
bar.** BASE22+A1+MS1ISO+TNS_b = 89.2/89.1 overall, 77.5/77.4 buried
(+1.8/+1.9 overall, +3.9/+4.0 buried vs 87.3/73.5). It is the ONLY candidate
that passed the seal's exchangeability band (AUC 0.5052, LB 0.4975,
LB[AUC-R0] +0.0101, PP 0.0239) and its buried condition (77.4 >= 71.0). It
was formally rejected FAIL-CLOSED on a missing S1 masked-vector artifact
(nan — the masking-ablation data was not produced for it), and the one-shot
seal is spent. Under the pre-registered rules: NO candidate advances, and
+1.9 is below the +2.0 bar regardless.

**The floors tell the mechanism:** FL0 (shallow/statics-adjacent + MS1)
appended under contrast (a) reaches 89.0/77.1 ~ (b)-full's 89.1/77.3. The
lift that exists does NOT live in the deep multi-channel trace encoder — it
lives in shallow structure the MS1-iso family and descriptors already carry.
The learned tensor added nothing the guards allow us to keep.

## Verdict (binding, per pre-registration)

PARKED — doc/63 style, with dignity: the run was decisive. If revisited, the
pre-registered path is a NEW round (fresh seal data), the S1 artifact fixed,
and the mandatory full-run entrapment read of a (b)-contrast candidate — the
+4.0 buried signal is the one number that argues for that round.

## What stands after the whole 48-hour arc

Scoring phase: closed (plateau; classifier exonerated). Localization: closed
(oracle +0.3). Windows: closed. The evidence rebuild's legitimate gains so
far: MS1-isotope family (+0.7 overall / +2.4 buried stacked, symmetric,
anchor-robust) and per-fragment vectors (+0.4) — cohort-projected, both
below the C++ bar individually; their full-run entrapment read together is
the natural next measurement. The guard apparatus (symmetry audits, sealed
decoy partitions, injection-supervised gates, entrapment-first verdicts) is
now the project's permanent instrument — and ahead of the published field.

## ADDENDUM (same night, ~22:00) — the anchor confound and what survives

The parallel analysis line retracted the unqualified "traces are squeezed
dry": every cohort trace experiment (wave-1, AXIC, and this tensor's derived
channels) scored ONE label-blind argmax anchor, and 81.5% of buried blocks
were anchored >5 s from the true apex (median 45.7 s) — the single-anchor
construction, not the trace evidence, produced the flat results. What
SURVIVES of this doc: the seal's asymmetry verdict (decoy-supervision poison
is anchor-independent), the park of THIS tensor as built, the guard
apparatus, and MS1-iso (+2.4 buried, measured under the anti-trace bias).
What is SOFTENED: "the lift lives in shallow structure" — the P0 factorial
(label-blind multi-candidate reference-correlated features, DIA-NN 1.7.12
mechanics) reaches 87.4/73.9 from BASE22 alone and 88.1/75.4 stacked (+3.9
buried, the program's best), with the rt_dn oracle at 93.9/87.4 showing the
trace information is real and was locked behind candidate selection. Next
step (pre-registered by that line): port the per-candidate features into
ODIA's own scorer and re-measure FULL-RUN at matched entrapment FDP.
