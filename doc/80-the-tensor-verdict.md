# 80 — The multi-channel tensor verdict: traces parked, one statics survivor gated

2026-08-31 evening. The full arc — design (TENSOR_DESIGN/TENSOR_PLAN in
shared/libv2/analysis77/), adversarial review (kimi + codex, 12 CHANGED
amendments), build (all gates PASS), 23 training arms on spock GPU1 (50 min,
zero failures), pre-registered evaluation, one-shot D_seal read. Everything
below was pre-registered before the runs.

## Verdicts

1. **Contrast (a) — decoy-supervised, the like-for-like primary — FAILED
   exactly as the symmetry audit predicted.** D_select AUC(dec>null) = 0.0048:
   the model learned near-perfect decoy detection from the derived channels'
   construction asymmetry. Appended: −6.5 pts. Shadow tripwire FIRED
   (Δsym −0.078); even the no-shadow arm seals at 0.087 (tight/heavy asymmetry
   suffices). SL2: capacity retry run (still bad) → (a) is PARKED by its own
   rule.
2. **The interference-recognition claim is NOT made** (pre-registered battery:
   partial Spearman ≈ 0 — the learned gate does not track fragment coherence;
   only the counterfactual sub-read passed). The learned-gating mechanism for
   interference removal did not validate on this data.
3. **The multi-channel TRACES add nothing under the artifact-immune contrast:**
   FL0 (all traces zeroed, statics only) appends 89.0/77.1 — equal to the full
   (b) model (89.2/77.5). Trace attribution FL5−FL4 = −7.6 under (a). The
   tensor experiment's central hypothesis — that tight/shadow/heavy/MS1-iso
   TRACE structure carries recoverable evidence — is refuted at cohort scale.
   Park, doc/63 precedent.
4. **One survivor, gated:** the contrast-(b) score (DN-vs-nonDN targets,
   decoys never seen) is the largest buried lift measured in the whole
   program — stacked 89.2/89.1 overall, 77.5/77.4 buried (+1.9 overall,
   +4.0 buried vs 87.3/73.5) — and its STANDALONE score passed the D_seal
   symmetry band (0.505, LB 0.498). But: (i) +1.9 is UNDER the pre-registered
   +2.0 bar — formally not cleared, no C++ investment; (ii) the seal's S1
   masking panel was missing for it (instrumentation gap; the one-shot is
   consumed — no re-seal); (iii) FL0 shows the lift is statics-INTERACTION
   driven, kimi F7's exact warning (mass-defect-manifold risk), and the
   statics-shuffle control was run only on the (a) side. Sanctioned
   continuation per §5/SL4: the FULL-RUN entrapment read — the only instrument
   that can distinguish a real statics-interaction gain from label-manifold
   exploitation. Nothing ships from the cohort.

## What the program bought

The negative is clean and triple-locked (audit prediction → (a) collapse →
seal), the FDR discipline caught what AUC never would (prior era: +0.067 AUC,
~22% entrapment FDP — this round rejected the same failure mode BEFORE any
investment), and the (b)/statics direction is a genuinely new, symmetric,
buried-class-targeted lead with a pre-registered gate in front of it.
