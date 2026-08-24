# The trace-classifier plan, revised after review

2026-08-24. doc/55 was reviewed by codex and kimi before any corpus was built.
Five findings change the design and two change what the result could ever mean.
Recorded as amendments so the original reasoning stays visible.

## A1. The apex-centring leak (both reviewers, independently, as the answer to
## "how does this produce an impressive number that means nothing?")

doc/55 said the model input is a 64-cycle window "centred on the candidate
apex". **Positives have an observed apex by construction; absent negatives do
not.** Centre positives on DIA-NN's peak and negatives on anything else -- the
predicted retention time, the strongest noise excursion, ODIA's best candidate
-- and the network learns the CENTRING PROCEDURE: coherent-peak-at-centre
against off-centre or diffuse. kimi's addition is the part that makes this
lethal: **every audit in doc/55 would have passed**, because the audit data
shares the construction. Great AUC, starkest at Q1, decoys ranked low, transfers
to S30 -- and useless in production, where every candidate is centred on its own
picked apex and the asymmetry does not exist.

**Fix, already applied to the corpus:** the dump centres the window on the
LIBRARY-PREDICTED retention time, by the same procedure for every class, and
there is no re-centring at training time. The model consumes the whole
predicted-RT-centred window.

**Fix to the evaluation:** the deployment distribution is ODIA's own candidates
at their own picked apex. That evaluation is now mandatory and is the only one
that can catch a residual centring artefact.

## A2. A taxonomy-free negative class (kimi)

Entrapment negatives are a different ORGANISM, and matching on charge, m/z,
length and fragment count does not fix it, because sequence survives in the
product-m/z set, the series and ordinal, the missing-fragment pattern and the
library relative intensities. Codex: product m/z + series + ordinal IS a partial
sequence representation, and a transformer can exploit fragment-spacing patterns
no scalar control would see.

**Measured, not argued:** the library's own relative intensities already differ
between the classes -- median top-fragment share of library intensity 0.2404 for
human against 0.2165 for Arabidopsis. The leak is real and present in the
corpus, whatever its cause.

**Fix: RT-SHIFTED SAME-PRECURSOR NEGATIVES**, kimi's proposal and the strongest
idea in either review. Extract the same human precursor's own fragment set with
the map intercept shifted +300 s. Same sequence, same charge, same product m/z
set, same series and ordinals, same library intensities -- so not one of the
organism confounds can operate. It is a real peptide's fragments measured where
that peptide is not.

Now dumped as a second arm. Entrapment is demoted to a secondary negative class,
and its honest description is codex's: "human-presence versus
Arabidopsis-entrapment under matched acquisition coordinates", which does NOT by
itself establish discrimination of absent HUMAN targets.

Caveat kept in view: shifted negatives teach localisation as well as presence,
so a model could learn "is there a coherent peak here" rather than "is this
peptide present". The two classes answer different questions and both are kept.

## A3. Hard negatives were being excluded (kimi)

Entrapment precursors DIA-NN calls at q <= 0.01 are VERIFIED FALSE POSITIVES:
real peptide masses, genuinely absent, and confidently misidentified. Drawing
negatives only from uncalled entrapment trains and tests on the EASY absent
population, which is why every metric would look good while the tail -- the only
region that matters -- goes unmeasured.

They become a held-out HARD-NEGATIVE set, never trained on. The corpus sampler
already routes entrapment before positives so they never contaminate the
positive class, but it pools called with uncalled; that split is now required.

## A4. The integration test would break the project's measuring instrument (kimi)

This is the finding with the longest reach. The model is trained to separate
entrapment from human. Integrate its score as a sub-score and it suppresses
entrapment precursors in the production run BY CONSTRUCTION -- so entrapment
stops functioning as an unmodelled null and **entrapment FDP reads optimistic**.
Disjoint precursor splits do not fix it: the model has learned the CONSTRUCT,
not the individuals.

Entrapment-calibrated FDP is the instrument every acceptance decision on this
project has rested on since doc/43. Consequences:

* An integrated model must be evaluated on DECOY-based FDR plus a held-out
  entrapment construct never seen in training -- not on the entrapment set it
  was trained against.
* The RT-shifted negative class does not have this problem at all, which is a
  second and independent reason to prefer it.
* A balanced-corpus precision is NOT an FDP and must never be reported as one
  (codex).

## A5. Feasibility: the inputs doc/55 assumed are not in the dump (codex, kimi)

`-out_chrom` writes `Precursor.Id, Decoy, Transition.Index, Product.Mz, RT,
Intensity` and nothing else. Absent: measured ppm deviation, MS1 traces,
mobility planes, candidate apex and bounds, and per-fragment library metadata
(joinable from the library parquet, unlike the first three, which live on side
planes and in `Ms1Traces` and are not serialised by this writer).

So v1 is traces plus library-joined metadata. MS1 -- one of DIA-NN's seven MS1
features and our only surviving orthogonal evidence -- is NOT in v1, and that
must be stated whenever a v1 number is quoted rather than discovered later.

Also fixed: the first dump attempt subsetted the library and the mass
calibration gate FAILED on the subset, leaving extraction at +-50 ppm against
production's +-10. Calibration is now PINNED (`-fragment_ppm 10
-fragment_ppm_offset -10.7372`) so the corpus is extracted under production
conditions, identically for every class.

## A6. The baseline ladder, which the transformer must climb (both)

doc/55 proposed beating the 19-scalar GBT. Both reviewers say that is too easy
and name the real baseline: **order statistics of the per-fragment quantities**
-- sorted per-fragment correlations to the consensus, per-fragment apex offsets,
per-fragment intensities -- through the existing GBT. Permutation-invariant by
construction, at roughly 1% of the cost, and it captures much of what attention
would. doc/45's finding that interference is DISTRIBUTED across fragments is
precisely why order statistics should work.

    1. 19 scalars + GBT                      (what ships)
    2. + per-fragment ORDER STATISTICS + GBT  (the real baseline)
    3. Deep Sets: shared encoder, masked pooling, MLP
    4. Transformer with attention across fragments

Each rung must beat the one below on held-out Q1, or the rung above is not
justified. Codex: "the transformer is not yet justified."

## A7. Representation: do not normalise away the variable under study (codex)

Per-group intensity normalisation removes exactly the low-abundance evidence
this exists to recover, and supplying the group total does not restore it --
signal-to-noise depends on per-fragment background, missingness and
heteroscedastic counting statistics, not the group total. Two normalised traces
with the same shape and total can have very different reliability.

Channels instead: `log1p` raw intensity at run-level scale; background-
subtracted intensity; a zero/missing mask; per-fragment local-noise estimate;
group total as an explicit covariate. Normalise for numerical stability only,
and report within fixed abundance bins so a gain cannot come from abundance.

## A8. Splitting and single-use validation (both)

* Split by PROTEIN / sequence-homology cluster before matching, not by
  precursor: paralogues and shared peptides leak across a precursor-level split.
* Fit matching and normalisation on the training split only.
* S30 is locked until ONE final evaluation. Repeated inspection makes it a
  second training set. It is converted and ready (10 GB mzpeak, same 24
  isolation windows as S08) and will not be touched until then.

## A9. What the comparative experiment can actually claim (codex)

"(b) >> (a) means our extraction is the problem" is not licensed. DIA-NN's XICs
come from its own mass windows, RT placement, calibration, resampling and
fragment selection, so the comparison confounds extraction QUALITY with
extraction SETTINGS. Licensed claim: *the exported DIA-NN representation carries
more label-predictive information under its own complete configuration.* To say
more requires matching precursor and fragment lists, RT centre and span,
resampling grid, and mass tolerances one family at a time.

## A10. Positives are a benchmark label, not ground truth (codex)

DIA-NN-confident means "what DIA-NN can confidently recognise", so the model
inherits its extraction, peak-shape and abundance preferences. This directly
undercuts the planned Q1 conclusion: the low-abundance positives are the subset
already easy enough for DIA-NN, while the hard real Q1 peptides are unlabelled
or mislabelled. Treat non-calls as UNLABELLED rather than negative
(positive-unlabelled learning), keep DIA-NN calls as one benchmark label, and
say "benchmark label" wherever "ground truth" is tempting.

## A11. Asymmetric label noise in the positives (vibe, MEASURED and large)

Both earlier reviewers missed this. A DIA-NN-confident positive is labelled from
DIA-NN's apex, but the window is centred on the LIBRARY-PREDICTED retention
time. When the prediction is off by more than the half-width, the window holds
only noise and still carries a POSITIVE label.

It is asymmetric by construction: RT-shifted negatives have their peak 300 s
away and entrapment negatives have none, so both are correctly "no peak here".
Only the positive class is diluted.

Measured on the corpus, |RT_DIA-NN - (1086.50 x libRT + 473.77)| over 38,983
positives:

    <=  30 s   54.8%
    <=  45 s   70.9%
    <=  60 s   79.3%
    <=  90 s   86.0%   <- the dump half-width
    <= 120 s   89.2%
    median 26.7 s   p90 128.0 s   p99 290.7 s

**14.0% -- 5,473 positives -- have no peak in their own window.** Training on
them teaches "no peak implies present", the exact inverse of the target concept,
in the class that has no other source of noise.

A second thing falls out that is bigger than the corpus. The 37.7 s p95 this
project quotes for the S08 map is a FIXTURE number; on the full library the same
constants give **p90 128 s and p99 291 s**. The map is far worse than the figure
in circulation, and doc/54 has just shown on Astral that window width costs more
than it buys. That deserves its own experiment independent of any classifier.

**Fix:** the 5,473 move OUT of the positive class and become UNLABELLED -- the
positive-unlabelled framing codex asked for in A10, not a new mechanism. They
are true positives our window does not contain, so they are uninformative for a
trace model rather than negative.

**The bias this induces, stated rather than hidden:** the surviving positives
are those whose iRT prediction is good, which correlates with sequence
properties. It is measurable -- composition and abundance before against after
the exclusion -- and that comparison is now part of corpus construction.

## A12. Training and evaluation would be centred differently (vibe)

A1 mandates evaluating on ODIA's own candidates at their own picked apex, to
catch a residual centring artefact. But training windows are centred on the
PREDICTED retention time. A peak at offset +X in training appears at offset 0 in
evaluation, so A1 removes a between-class leak and introduces a
train/evaluation distribution shift in its place.

**Fix: random jitter augmentation.** Shift the window centre by a random offset
each epoch, making the model invariant to where in the window the peak sits.
That simultaneously removes any residual ability to use centring as a class cue
and aligns training with an apex-centred evaluation -- cheaper and more robust
than matching the two centrings exactly, and it attacks A1's leak from the other
side: a model that cannot use position cannot be fooled by it.
