# AI denoising / deconvolution of DIA mass traces: what is actually published

2026-08-21. Deep-research sweep, 101 agents, 19 sources surviving adversarial
verification (3-vote, 2/3 to refute). Motivated by doc/45: ODIA's rejected
DIA-NN-confident precursors carry all twelve fragments at 0.57x signal-to-noise
with apexes 20.8 s apart against 2.8 s -- distributed interference, which no
SELECTION-shaped design can touch.

## The headline

**Nothing in the literature is a "denoiser" in the image-processing sense.** No
published DIA tool takes a contaminated fragment trace and emits a cleaned one
for downstream use. The work splits into three families, and only one of them
survives the objection codex raised against our own low-rank plan.

## Family A -- blind low-rank factorisation. Cannot name its components

**CANDIA** (Buric, Zrimec, Zelezniak, *Patterns* 1(9):100137, 2020) partitions
all runs into (m/z x RT x sample) tensors, one per isolation window x RT window,
**with MS1 survey scans aligned into the same tensor as their MS2 fragments**,
and runs non-negativity-constrained PARAFAC on GPU. Each component is a triplet
(fragment spectrum, elution profile, sample abundance).

Two structural limits that rule it out for us:
  * **the third mode is SAMPLE.** Uniqueness comes from trilinearity ACROSS
    RUNS. This is not a single-run bilinear factorisation and cannot be applied
    to one run in isolation.
  * it decomposes all analytes in a slice jointly (F = 10-90 components), not
    the traces of one precursor.

And on identification it LOSES: "Crux identified a total of 2,014 proteins using
the DIA-Umpire pseudo-spectra... in contrast only 684 proteins (1,553 peptides)
were identified on the output from CANDIA." Its own identification pitch is
replicate reproducibility and TIC coverage, not yield.

**DIA-NMF** (Karaki et al., EUSIPCO 2024) is the closest published match to the
framing we independently arrived at: an explicit bilinear model of ONE
precursor's ion x RT matrix, `X ≈ WH`, rows being candidate MS1
parent/isotope/adduct traces AND MS2 fragment traces, rows of H the shared
elution profiles, solved by sparse non-negative BSS (nGMCA: L1 on W,
non-negativity on both factors, forward-backward proximal splitting).

But: the domain is untargeted **metabolomics**, not peptide DIA -- transfer is
our inference, not theirs. And **all three of its downstream claims were REFUTED
by verification** (0-3, 0-3, 1-2): its identification advantage, its advantage
in the low-intensity co-eluting regime, and its automatic rank selection. Only
the architecture survived.

## Family B -- known-dictionary unmixing. This is the one that survives

**Specter** (Peckner et al., *Nat Methods* 15:371-378, 2018) models each
multiplexed MS2 spectrum as

    S = Lc + N = c1*L1 + ... + cm*Lm + N

and solves per scan by non-negative least squares, using "only the library
precursors whose m/z ratios fall into the precursor isolation window". **The
interferers ARE columns of L.** Interference is an explicitly modelled component
of the fit, not an outlier to clip. The paper says so directly: "Ambiguous
shared features are typically deemed 'interferences' and excluded from
consideration in analysis of DIA data."

Solving scan by scan yields a coefficient matrix whose rows are deconvolved
elution profiles -- precursor i's contribution is the outer product
c_i(t) ⊗ L_i, i.e. **exactly the rank-1 collapse**, obtained without blind
factorisation.

**Siren** (Hu, Lu, Bilmes, Noble, *JPR* 18(1):86-94, 2019) does the same at MS1:
`argmin_{B>=0} ||Y - XB||²_F + λ||B||`, X columns being averagine theoretical
isotope distributions, solved by non-negative LARS-LASSO.

**Why this matters to us.** codex's objection to our low-rank plan was that no
trace-only factorisation can NAME which recovered component is the query peptide
and which the interferent -- fatal, because decoys co-elute with confident
targets at median dRT 0.21 s. Specter answers it by construction: the components
are not discovered, they are the library entries in the window. The naming
problem disappears because the dictionary is labelled. **This is the same
insight our own backlog reached from the other direction** -- "the weight comes
from OTHER precursors' evidence" -- expressed as linear unmixing instead of
counting.

Specter is also explicitly positioned AGAINST the fragment-correlation paradigm
that ODIA, OpenSWATH and Spectronaut all use, arguing correlation "cannot
rigorously account for precursor co-fragmentation and so cannot separate
precursors sharing fragment m/z".

## Family C -- learned scorers that never emit a clean trace

**Alpha-XIC** (Bi-GRU x2, hidden 64, dropout 0.5, self-attention; traces
intensity-normalised, interpolated to 32 points, Savitzky-Golay smoothed),
**DIA-BERT** (CNN + transformer over a 330x16 MS1+MS2 peak-group matrix), and
**Dear-DIA^XMBD** (VAE + triplet loss replacing point-to-point XIC correlation)
all learn co-elution implicitly and output a SCORE, not a trace.

**DIA-BERT contains no explicit denoising step at all** -- interference is
handled implicitly by self-attention, plus training data deliberately simulated
with high inter-peptide interference. And that interference-simulated corpus
trains only the QUANTIFICATION model, not the identification model.

**AutoMS** is the one true denoising autoencoder found (verified from source,
the article being paywalled): one trace at a time, resampled to 50 points around
the apex, MS1-only, and the denoised output is used as a reconstruction-distance
quality SCORE rather than being passed downstream.

## Identification evidence is thin, and mostly refuted

Claims that interference handling improves IDENTIFICATION rather than
quantification:

    Siren        1,711 -> 2,516 peptides at 1% FDR on raw MS2      VERIFIED 3-0
    Alpha-XIC    +9.4-16.2% precursors appended to DIA-NN          VERIFIED 2-1
    CANDIA       1,553 vs 1,111 peptides (unimodality rank)        VERIFIED, but
                 the control arm ALSO uses the shape prior, and CANDIA loses
                 to DIA-Umpire outright
    DIA-NMF      identification advantage                          REFUTED 0-3
    DIA-BERT     +22% precursors / +51% proteins vs DIA-NN         REFUTED 0-3
    Dear-DIA^XMBD 31,439 vs 15,784 peptides                        REFUTED 0-3

Every surviving case is single-group, internally controlled, and measured
against baselines several tool generations old.

## The gap that matters most to us

**No source in the verified evidence base tests whether a learned denoiser or
deconvolver inflates decoy scores as much as target scores.** Every tool places
deconvolution BEFORE scoring and trusts downstream target-decoy machinery.
Published practice has three tiers:

  * **Specter** deconvolves decoys JOINTLY with targets (two-pass hybrid
    target+decoy library), plus intrinsic checks -- 2% with an E. coli decoy
    library, 0.5% with synthetic decoys. But the authors' stated rationale is
    "to avoid distorting the quantifications of non-decoy library members" --
    protecting quant, not controlling the FDR hazard -- and pass 2 uses a
    target-only library, so the threshold comes from a different optimisation
    than the one it is applied to.
  * **Siren / DIA-BERT** add a separate FDR layer for the interference step, or
    an external two-species entrapment audit.
  * **Alpha-XIC** does nothing: the learned score is appended and DIA-NN's own
    validation is trusted. Its verifier noted the network "is trained on peak
    groups from the same run it is applied to... and its score is fed back into
    the same target-decoy validation".

This is precisely the hazard we measured on our own data -- pruning tripled the
target median and tripled the decoy median with it, leaving AUC unmoved. **The
literature does not solve it and largely does not test it.**

One REFUTED claim is worth recording because it points at our v8 experiment: a
verifier rejected the assertion that shuffled in-silico libraries containing no
real peptides still produced 2,691 (Skyline) / 642 (DIA-NN) "quantified"
proteins at 1% FDR, ~185x expected. Refuted, so not evidence -- but the shape of
the claim is exactly what our decoy diagnostic found independently (~100x
understatement in the confident tail).

## Also worth knowing: item (4) is unanswered

How DIA-NN, Spectronaut, OpenSWATH, DIA-Umpire, CsoDIAq and Avant-garde actually
implement interference correction and profile clipping is **essentially
UNADDRESSED** by the surviving evidence -- every statement about them in the
verified set is second-hand characterisation by a competitor. Our vault's
DIA-NN interference note and the 1.7.x source reading remain the better source.

## Sources

19 verified. Primary: CANDIA (Patterns 2020), Specter (Nat Methods 2018), Siren
(JPR 2019), DIA-NMF (EUSIPCO 2024), Alpha-XIC, DIA-BERT, Dear-DIA^XMBD, AutoMS.
Full list and per-claim evidence in the workflow output.
