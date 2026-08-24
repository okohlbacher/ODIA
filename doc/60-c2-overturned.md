# C2 overturned: the library prior was the whole thing

2026-08-24. doc/58 concluded that trace models do not beat the shipped scalars.
doc/59 withdrew that on review. This is the corrected measurement.

## The result

All on the SAME three-way protein split (folds 0-2 test, 3-4 validation, 5-9
train), the same apex-centred contrast, the same 10,688 test precursors.

    model                                          test AUC    vs scalars
    rung 1: 19 shipped sub-scores + GBT              0.8629            --
    trace transformer, ANONYMOUS traces      0.8662 +- 0.0004      +0.0033
    trace transformer, + LIBRARY INTENSITY   0.8818 +- 0.0013      +0.0189
    trace transformer, + FULL DESCRIPTORS    0.8855 +- 0.0014      +0.0226

**+0.0226 over the shipped scalars**, about 16x the seed spread and 10x the
sampling standard error (~0.0023 on 21,376 test pairs). The ladder is monotone
and every step is outside the spread.

**The descriptors cannot be leaking.** Library intensity, product m/z, ion
series, ordinal and fragment charge are IDENTICAL between a positive and its own
shifted negative -- same precursor, only the window moves -- so a model given
them alone scores exactly 0.5. All of their +0.019 therefore comes from
INTERACTION with the traces. That is the cleanest available demonstration that
the gain is structure rather than metadata, and it is a property of the paired
design rather than of an argument.

**These are undertrained.** Five of six runs picked the LAST epoch on
validation, so the schedule was too short and the numbers are a floor.

## What changed, and in what proportion

Two fixes, both from the review, and they are not equal:

**The encoder, worth +0.0017.** doc/58 ended each fragment's encoder in
`AdaptiveAvgPool1d(1)`, averaging the time axis away BEFORE any cross-fragment
attention -- kimi's objection. Replacing it with a conv front end plus a
bidirectional GRU takes anonymous traces from 0.8645 to 0.8662. Small, and it
was enough on its own to put the trace model nominally ahead of rung 1 while
training on 27% LESS data (a fold now goes to validation).

**The library prior, worth +0.0156, and the remaining descriptors +0.0037.** This is the one all three reviewers found
independently, and it is an order of magnitude larger than the encoder fix. The
models could not see which fragments SHOULD be strong, while rung 1's
`var_library_corr`, `var_library_dotprod` and `var_library_rmsd` all compare
observed against expected. Codex's framing was exactly right: an
expected-intensity vector is AUC 0.5 by itself and becomes predictive through
its RELATIONSHIP with the observed traces.

So doc/58's negative was an artefact of withholding one input, and the
withheld input was the one the hypothesis was about.

## What this says about the original hypothesis

doc/55 asked whether collapsing twelve fragments into nineteen scalars is where
the low-abundance information is lost. The answer, on this contrast, is **yes**:
a model that sees the same evidence per fragment, without the summarisation
step, extracts +0.0189 AUC that the scalars do not.

Note what the winning model has that `var_library_corr` does not. Both compare
observed against expected. The scalar reduces that comparison to ONE Pearson
correlation over the fragment vector; the model keeps each fragment's expected
intensity beside its own trace and lets attention weigh them against each other.
The information was always there -- the scalar threw away which fragments
disagreed, and that is precisely doc/45's finding that interference is
DISTRIBUTED across fragments rather than concentrated.

## Method notes that make the number believable

* **Test read once per seed**, at the epoch chosen on validation. doc/58
  evaluated test every epoch and reported the final, which was defensible but
  not clean.
* **Rung 1 rerun on the identical split.** It had been training on folds 3-9
  while the neural arm trained on 5-9; on the matched split rung 1 measures
  0.8629 rather than 0.8649, and the comparison is now about representation
  rather than about who got more data.
* **Three seeds, spread reported.** doc/58's 0.0004 was uninterpretable without
  one; this margin is 15x that spread.
* **Arm alignment asserted in the trainer**, which doc/58's script never did.

## Still open, and none of it is small

1. The contrast is still "here rather than 300 s away", not "present at all".
   The presence question needs the entrapment arm, whose leakage floor must be
   MEASURED because the true-metadata control will not be 0.5 there.
2. `ent_hard` is 17. The tail is unmeasured.
3. The protein split is 98.1% disjoint, not 100% -- `A;B` and `B;C` hash
   independently. 1.9% accession bleed, unlikely to move +0.0189, but it should
   be a connected-component split.
4. Positives are DIA-NN-confident AND RT-well-predicted; the 14% excluded are
   0.83x abundance, so the low-abundance tail is depleted in exactly the regime
   the hypothesis is about. The win may be understated, or may not survive there.
5. S30 remains locked and untouched.
6. No integration test. A sub-score that beats the scalars on a paired
   localisation contrast has not been shown to add identifications at matched
   FDP, and doc/56 A4 warns the entrapment instrument breaks if this model is
   integrated naively.

## The entrapment contrast has a leakage floor of 0.7060, measured

doc/56 A2 said the entrapment negatives' organism confound could not be argued
away and the floor had to be measured. It is now measured:
`scripts/analysis/ent_floor.py`, descriptors ONLY -- library relative intensity,
product m/z, ion series, ordinal, fragment charge, fragment mask -- and **no
trace information whatsoever**:

    LEAKAGE FLOOR on the entrapment contrast   AUC 0.7060
    the same floor on the shifted contrast     AUC 0.5000

**Composition matching did not remove the taxonomy signal.** Matching on charge,
precursor m/z, peptide length and fragment count left the joint structure --
which product masses, in which series, at which ordinals, with which expected
intensities -- fully able to separate human from Arabidopsis. That is exactly
codex's point in doc/56 A2: product m/z plus series plus ordinal is a partial
representation of the peptide sequence, and matching four marginal distributions
does not match it.

So any trace model on the entrapment contrast must beat **0.706**, not 0.5, and
the honest quantity is its increment over that floor rather than its raw AUC.

This retrospectively vindicates the RT-shifted design as the PRIMARY contrast.
Its floor is 0.5 by construction and measured at 0.5000, so the +0.0226 in the
table above is not competing with any taxonomy signal at all.

## The presence question, answered against its floor

The entrapment contrast asks what the project actually cares about -- present
against absent, not here against there.

    contrast   floor (descriptors, no traces)   trace model (full)
    shifted                    0.5000            0.8855 +- 0.0014
    entrapment                 0.7060            0.8712 +- 0.0022

The trace model is far above the taxonomy floor on the presence question, so the
traces carry real presence information and not merely organism. But **0.8712
against a 0.7060 floor is not "+0.165"** -- AUC does not decompose that way, and
quoting a difference of two AUCs as an increment would be exactly the kind of
arithmetic this project has had to retract before.

The clean isolation is a trace model with NO descriptors on the entrapment
contrast: whatever it scores is presence information that cannot have come from
the library pattern. That run is queued. Note it is not perfectly clean either
-- Arabidopsis precursors sit in different m/z windows and therefore different
interference environments, so a trace-only model retains an indirect route to
taxonomy -- but it is much tighter than the descriptor-bearing model.

Worth stating plainly: the shifted contrast gives the larger number (0.8855) AND
has the floor that is 0.5 by construction. The entrapment contrast gives the
question we want answered and a floor that has to be subtracted by argument
rather than by design. Neither alone is sufficient, which is why doc/56 kept
both.
