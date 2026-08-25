# Every contrast that teaches the right task leaks

2026-08-25. Measured before training anything on them, which is the only reason
the numbers below can be read as a floor rather than as a result.

## The three floors

Descriptor-only models -- library relative intensity, product m/z, fragment
charge, ion series, ordinal, precursor charge/m/z/length, fragment mask -- with
**no trace information at all**:

    contrast     floor    what the floor is made of
    shifted     0.5000    nothing: a positive and its own +300 s window are the
                          SAME precursor, so every descriptor is identical
    decoy       0.7086    mass-defect manifold + predicted intensities of an
                          unnatural sequence
    entrapment  0.7284    taxonomy: Arabidopsis against human

## The decoy floor was the surprise

I expected it near 0.5. A decoy is the same peptide's own permuted sequence, so
naively nothing about it should be visible without looking at the data. It is
0.7086.

The mechanism is doc/50's, arriving from a new direction: `mutate` maps every
residue onto `LLLVVLLLLTSSSSLLNDQE` -- eight residue masses, no aromatics, no
sulfur, no G/A/P. So a decoy's product m/z sit on a different mass-defect
manifold, and its library intensities are a predictor's output for a sequence
that cannot exist. Both are legible from the descriptors alone.

doc/50 said decoys "pick up random coincidence while absent targets pick up
structured signal" and closed the decoy branch on that. This is the same fact
measured at the input rather than at the output.

## The tension this creates, stated plainly

    contrast     floor     task it teaches
    shifted     0.5000     localisation -- measured 4x WORSE at matched FDP
    decoy       0.7086     presence, and it is what the shipped scorer uses
    entrapment  0.7284     presence, with an organism confound

**The only contrast with a floor of 0.5 is the one that teaches the wrong
task.** Every contrast that asks the presence question carries ~0.71 of
free discrimination that has nothing to do with whether a peptide is present.

That does not make the decoy contrast unusable -- the shipped discriminant
trains against exactly these decoys and is not thereby invalid, because the
floor is shared by every method and cancels in a comparison BETWEEN methods.
What it forbids is reading a raw AUC as evidence of presence discrimination, or
comparing an AUC across contrasts with different floors.

## Consequence for acceptance

Unchanged and reinforced: **the acceptance test is identifications at matched
entrapment FDP** (`score_corpus.py`), not AUC on any contrast. doc/62 measured
a model at 0.9302 AUC that was 4x worse in identifications. A floor of 0.71
means an AUC of 0.85 on the decoy contrast might be almost entirely floor.

## A defect this work fixed on the way

The first descriptor builder matched only targets, giving every decoy row
all-zero library intensities and no precursor charge -- because ODIA regenerates
decoys internally and they are absent from the input library, and because that
library names a decoy `<target_id>_decoy` while the trace dump reconstructs the
bare target id. Training the decoy contrast on that corpus would have learned
"descriptors present means target", scored near-perfectly, and passed every
leakage control, since those were built for the paired shifted contrast where
metadata is identical by construction.

Fixed by `-stop_after library -out_lib` plus suffix stripping
(`build_descriptors_v2.py`); coverage is now 100% on both classes. The same fix
corrected `Fragment.Type`, which had been read as the enum's numeric value
rather than y/b. The entrapment floor moved 0.7060 -> 0.7284 on the corrected
descriptors.

## Training on the right task: 3,377 -> 9,454, and why it cannot be believed

The transfer test, rerun with a model trained on the ENTRAPMENT (presence)
contrast instead of the shifted (localisation) one, GRU fixed:

    ranking                        1% FDP    2%       5%       10%
    ODIA shipped DScore            13,646   14,697   16,287   17,829
    shifted-trained (doc/62)        3,377    3,770    4,593    5,826
    entrapment-trained              9,454   11,209   13,912   16,505

**2.8x better at 1% FDP than the shifted model**, closing to 93% of the shipped
scorer by 10% FDP. doc/62's diagnosis was right: the task was the problem, not
the representation.

**And the number is inflated by construction.** doc/56 A4 wrote this down before
any of it was built: a model trained AGAINST entrapment suppresses entrapment in
whatever it ranks, so entrapment FDP -- which is (e/r)/t over ordinary targets --
reads optimistic. The model has learned the CONSTRUCT the ruler is made of. A
protein-disjoint split does not help, because the split separates individuals
and the leak is at the level of the construct.

So 9,454 is an upper bound of unknown tightness, and the honest comparison is
still missing.

**The decoy contrast is the one that can be measured on this ruler.** Decoys
appear in neither the numerator nor the denominator of entrapment FDP, so a
decoy-trained model has not seen the ruler at all. Its floor is 0.7086 rather
than 0.5, which affects how its AUC may be read, but not the validity of its
identification count. That run is what decides whether any of this transfers.

Note also what the GRU fix was worth on this contrast: 0.8712 -> 0.8827
(+0.0115), on a corpus where half the encoder had been contributing a state that
had seen one time point.

## The honest measurement: ~22% of the shipped scorer

The decoy-trained model, scored on entrapment FDP -- a ruler it has never seen,
because decoys enter neither side of it:

    model              trained against   1% FDP    2%       5%       10%
    ODIA shipped       decoys, full run  13,646   14,697   16,287   17,829
    entrapment-trained entrapment         9,454   11,209   13,912   16,505
    shifted-trained    own +300s window   3,377    3,770    4,593    5,826
    decoy-trained      decoys             2,989    3,802    5,257    7,121

**Both models that never saw entrapment land at ~3,000.** The one that looked
good was trained against the construct the ruler is made of.

**That difference measures the circularity doc/56 A4 warned about: ~3x.** 9,454
against 2,989 for models of the same architecture, same corpus, same acceptance
test, differing only in whether training touched the entrapment construct.
A4 predicted the direction before any of this was built; this is its size.

### What this establishes, and what it does not

Established: **the trace-model direction as constituted does not beat the
shipped scorer.** Two honest measurements agree at ~22% of it at 1% FDP. AUC
said 0.93 against 0.86 and was wrong about what mattered, three times now.

NOT established -- the comparison is not like-for-like, and in the shipped
scorer's favour:

* **Training volume.** The shipped GBT is a semi-supervised discriminant fitted
  on ~1.5M peak groups FROM THIS RUN. The trace model saw 16,568 precursors from
  a protein-disjoint subset. Two orders of magnitude, on the run being scored.
* **Population.** The trace model only ever sees apex-centred windows of
  precursors that HAVE a picked candidate. It has no way to rank the 35% that
  do not, and the shipped scorer implicitly does.
* **Corpus.** 86,809 precursors here against millions in a real run; the
  ranking problem is not the same shape.

So the honest summary is: on equal footing at the tail the trace model is far
behind, and the footing is not equal. Fixing the footing is a much larger piece
of work than anything attempted so far -- it means training on a full run's peak
groups, not on a curated 83k corpus.

### Where that leaves the direction

The plan put the neural work first and it has now failed its acceptance test
twice, honestly. The representation work is sound and reusable -- library prior,
precursor descriptors, Gaussian prior, the GRU fix, and a corpus with measured
floors for three contrasts. The task and the training volume are the problems,
and neither is cheap.

**Recommendation: park it and take the plan's step 2 and 3.** The v3 regression
is 1,302 identifications sitting in a flag that was flipped for an experiment
and never flipped back, and the Astral window is +117%. Both are measured, both
are cheap, and neither depends on any of this.
