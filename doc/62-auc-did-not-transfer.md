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
