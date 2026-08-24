# The findings, reviewed: C2 is withdrawn

2026-08-24. doc/58's results reviewed by codex, kimi and vibe. All three
independently found the same primary defect, and it is fatal to the main claim.

## The one they all found

**The trace models never receive library intensities or fragment identity.**

Rung 1's `var_library_corr`, `var_library_dotprod` and `var_library_rmsd`
compare the observed fragment pattern against the EXPECTED library pattern. The
neural models get three channels -- log1p raw, background-subtracted, nonzero
mask -- plus `log1p(t.sum(2))`, which is the OBSERVED area, not a library-derived
scale. They are blind to which fragments should be strong, to product m/z, to
ion series and to ordinal.

So the comparison was never "19 scalars against raw traces". It was **scalars
WITH a library prior against anonymous traces WITHOUT one.**

Codex puts the mechanism precisely: an expected-intensity vector can have AUC
0.5 by itself and still be highly predictive through its RELATIONSHIP with the
observed traces. That is exactly the quantity `var_library_corr` computes, and
`var_library_corr` is the feature measured to collapse hardest at low abundance
(0.532 at Q1 against 0.913 at Q5) -- the precise signal the whole hypothesis is
about.

**This also invalidates my own ablation.** I split the 19 into "intensity-derived
14" and "side-channel 5" (mass, mobility, MS1) and concluded the comparison was
fair because the side channels add only +0.003. But three of my "intensity-
derived" fourteen are LIBRARY-COMPARISON features, which the tensor cannot
express either. The tensor-comparable set is smaller than 14 and I have not
measured it.

**C2 is withdrawn.** What the experiment licenses is far narrower:

> Anonymous per-fragment intensity traces, encoded through a temporally
> average-pooled CNN, on DIA-NN-confident and RT-well-predicted positives, do
> not beat 19 scalars that include library-comparison features.

That is not evidence that summarisation loses nothing.

## What else broke

**C3's arithmetic is wrong.** I wrote "two thirds of rung 1's advantage was
localisation". Deep Sets closes 0.0220 -> 0.0038, which is **83%**; the
transformer closes 0.0043 -> 0.0004, **91%**. Neither is two thirds. And codex
notes apex-centring changes both localisation AND the negative-selection
distribution, so even 83% cannot be attributed to localisation alone.

**C4 was overclaimed.** "Metadata leakage is arithmetically impossible, proven
by the 0.5000 control" -- the control tests fragment count alone, and a single
arm-invariant feature returning 0.5 is nearly tautological. Kimi: it proves the
pairing is not broken, not that leakage is absent. The correct statement is
codex's: *arm-invariant metadata alone cannot classify correctly paired
examples*. The real guarantee comes from the pairing argument, and the control
is a sanity check on it.

**C1 needs uncertainty it does not have.** One protein partition, one seed, one
run, test evaluated every epoch. All three reviewers say 0.0004 cannot support
"match". The defensible statement is **"no measured improvement in this run"**.
vibe adds that with ~21k test pairs the standard error is about 0.0016, so the
transformer's per-epoch maximum of 0.8684 is more than 2 SE above rung 1 --
which is not a licence to quote the maximum, but is a reason to think the
question is live rather than settled.

**The reviewers disagree on one thing, and the disagreement is the product.**
On whether apex-centring inflates or deflates the trace models: vibe says
INFLATES (positives centred on a real peak, negatives on noise, so a localisation
task becomes a detection task); codex and kimi say DEFLATES (negatives are
enriched for the most peak-like thing in a wrong window, which is a harder
negative). Both mechanisms are real and they act in opposite directions. Nobody
knows the net, and the measurement cannot separate them.

## What held

* **Arm alignment**: verified -- ids and decoy flags identical across the two
  arms. The invariant was real; codex was right that the neural script never
  checked it, and it now must.
* **Protein split**: measured rather than asserted. **1.90% of accessions appear
  in both train and test** (1,933 of 101,532), because `A;B` and `B;C` hash
  independently. Small, and unlikely to explain 0.0004, but codex was right that
  disjointness was an intention. A connected-component split over accessions
  fixes it.
* **The centring confound is real** whatever its sign: +0.018 for Deep Sets.
* **The side-channel ablation stands** on its own terms: mass, mobility and MS1
  add +0.003 over the other sixteen, reaching 0.681 alone.

## The false-negative mechanism to bet on

Codex bets representation (missing library prior). kimi bets representation too,
and adds a specific architectural fault: `AdaptiveAvgPool1d(1)` mean-pools each
fragment's 96 cycles to a single vector BEFORE any cross-fragment attention, so
peak shape, width and apex position are largely averaged away -- and the jitter
augmentation deliberately forbids using absolute position. The models are
structurally prevented from using the temporal dimension, then asked a question
about time. vibe bets label noise plus capacity, noting the excluded 14% are
0.83x abundance, depleting exactly the low-abundance tail the hypothesis lives in.

They are not competing explanations. All three are true simultaneously and all
three are fixable.
