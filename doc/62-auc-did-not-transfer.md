# The AUC did not transfer. The trace model is 4x worse at matched FDP.

2026-08-24. `scripts/analysis/score_corpus.py`, the first measurement of this
work in the project's actual currency.

## The result

Whole corpus, own-window arm, 86,809 precursors with a located candidate and a
DScore. Entrapment ratio r = 19,708/44,791 = 0.44000. FDP is (e/r)/t over
ordinary targets, doc/43's definition.

    ranking                    FDP<=1%    FDP<=2%    FDP<=5%   FDP<=10%
    ODIA shipped DScore         13,646     14,697     16,287     17,829
    trace model (AUC 0.9308)     3,377      3,770      4,593      5,826

**A 4x deficit at 1% FDP**, from the model that beats the shipped scorer by
+0.068 AUC on the contrast it was trained on.

## Why, and it was foreseeable

**The shifted contrast teaches the wrong task.** The model was trained to
separate a precursor's own window from THE SAME precursor's window 300 s away.
That is localisation: "is the peptide here rather than there". Ranking real
precursors against entrapment and decoys is presence: "is this peptide in the
sample at all". A model that has never been shown an absent peptide has no
reason to rank one low, and the transfer measures exactly that.

doc/56 A2 recorded the tension when the contrast was chosen -- "shifted
negatives teach localisation as well as presence... the two classes answer
different questions" -- and doc/58 repeated it. It was written down every time
and never treated as a risk that could be decisive. It is decisive.

**AUC is blind where identifications live.** At 1% FDP the shipped ranking is at
depth 13,794 of 86,809 -- the top 16%. AUC integrates over the whole ranking and
is dominated by the bulk, so 0.93 against 0.86 can coexist with a 4x deficit in
the top few percent. I flagged this in the review brief as a risk and then
reported seven AUC comparisons before measuring it.

## What survives

* The ABLATION is unaffected: library prior +0.0204, precursor descriptors
  +0.0073, Gaussian prior, the encoder fix. Those are internally valid
  comparisons on a fixed task, and they say what information a trace model can
  use. They do not say that task is the right one.
* The leakage controls are unaffected -- 0.5000 arm-invariant floor, 0.4954
  label shuffle, 0.7060 measured taxonomy floor on entrapment.
* The corpus, the tensors, the descriptors and the checkpointing are all reusable.

## What does not survive

Every sentence of the form "the trace model beats the shipped scorer" needs the
qualifier **on the shifted-window contrast**. It does not beat it at ranking
real precursors, and ranking real precursors is the job.

## Next, and the order matters

1. **Retrain on the ENTRAPMENT contrast and re-run this test.** That contrast
   asks the presence question. Its taxonomy floor is 0.7060 and must be
   subtracted by argument rather than design, which is worse epistemically --
   but a model that answers the right question badly is more useful than one
   that answers the wrong question well.
2. **Train against DECOYS**, which is what the shipped discriminant does. The
   decoys are already in the corpus (71,021, never used). That is the like-for-
   like comparison and it has been sitting unused the whole time.
3. Only then is an integration test meaningful.

The honest summary of the day: the representation work is sound and the task
definition was wrong.

## The attribution, finished on a common basis

    full   (traces + all descriptors, no scalars)   0.9169 +- 0.0011
    allside (+ the 5 MS1/mass/mobility scalars)     0.9196 +- 0.0003   +0.0027
    allint  (+ the 14 intensity-derived scalars)    0.9292 +- 0.0012   +0.0123
    all     (+ all 19)                              0.9302 +- 0.0003   +0.0133

Settled: the combined model's gain is in the INTENSITY-derived scalars, not the
side channels. My prediction was the opposite and is now refuted on a single
basis rather than across two. The corollary stands: extending the tensor with
MS1, mass deviations and mobility -- doc/56 A5's owed work -- is not worth it.

## Review of the transfer, and one real bug

**The bidirectional GRU discards most of the backward context** (codex).
`o[:, -1, :]` takes the output at the LAST timestep: for the forward direction
that summarises the whole sequence, but the backward direction's state there has
seen only the final time point. It should be `h_n`, or forward[-1] concatenated
with backward[0]. The 0.9302 stands as an empirical number; the explanation that
the GRU preserves the time axis in both directions is false of this code, and
the encoder has been running at roughly half its intended capacity.

**Codex reached the transfer diagnosis independently and sharpened it.** The
apex-centred paired task keeps only positives with a picked candidate in BOTH
arms, then crops each around its own apex -- so the model learns "which of two
ALREADY-PICKED candidates looks more peptide-like", among the subset where both
exist. That is candidate re-ranking. It says nothing about the 35% with no
candidate, and "presence discrimination" overstated it.

**Two protocol claims that the code does not honour**, both mine:
* "test is read once" -- test is evaluated at every validation improvement.
  Not algorithmic leakage while nothing acts on the intermediate values, but the
  comment asserts something the code does not do.
* the protein split falls back from protein group to precursor ID silently when
  a lookup misses; coverage was measured (71,021 misses, all regenerated decoys,
  which never train) but the assertion was never made explicit.

**And a warning for the integration that has not happened yet:** adding the
model as a 20th GBT feature needs OUT-OF-FOLD predictions. A network score
computed for rows whose proteins the network trained on makes the stacked
classifier's FDR optimistic. Cross-fitting is required, or the integration must
be a standalone re-ranker evaluated on untouched proteins.
