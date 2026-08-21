# Staging DIA-NN's confident set through ODIA

2026-08-21. Full S08 diaPASEF, our 9.6M library, DIA-NN 2.0 on the SAME library
(`dn_xic.parquet`, 39,149 precursors at q <= 0.01, with XIC output). ODIA arm is
`full_v5.tsv`, production configuration.

`scripts/analysis/d3_stages.py`.

## The measurement

Five nested stages, each a subset of the one above:

    A  in our library                                 39,149  100.0%
    B  ODIA proposes >= 1 candidate                   23,058   58.9% of A
    C  a candidate within +-20 s of DIA-NN's RT       20,315   51.9% of A
    D  that candidate is the top-DScore one           19,430   95.6% of C
    E  it passes q <= 0.01                            12,292   63.3% of D

The stage losses, largest first:

    A->B  no candidate at all      16,091
    D->E  rejected by the threshold 7,138
    B->C  wrong peaks only          2,743
    C->D  ranked below a competitor   885

## What this establishes

**Within-precursor ranking is not the gap.** 95.6% at D: when the right peak is
in the candidate list, the classifier puts it first. Weeks of feature work were
aimed here.

**The largest single loss is upstream of scoring**, and it is larger than every
other loss combined.

**The precursors we lose are faint.** Median DIA-NN quantity 16,749 against
45,697 for the ones we find -- 2.7x. Same fragment count in the library (12.0
vs 12.0), same |dIM| against DIA-NN's observed 1/K0 (0.022 vs 0.024), heavier
(median m/z 731.9 vs 669.7).

## What it does NOT establish, and this is the correction

I attributed A->B to Gate C. Both reviewers refused it and they were right.

The argument was: the picker reports `no_hit_anywhere = 0`, so every precursor
that reached the correlation loop produced a candidate; therefore anything
absent was rejected before the loop; and the only material pre-loop reject in
the log is Gate C at 7,024,139 targets.

Kimi found the holes by reading the code. `Session::add` returns silently on
`tc == 0` and increments nothing. A precursor covered by no isolation window
never reaches `add` at all. `min_fragments_at_apex` can discard every candidate
AFTER the picker returned a non-empty list. None of the three appears in any
counter, and the log's own `no points (<3) 1,119,490` is a second candidate
explanation that the heavier m/z of the lost set actively supports.

Codex made the general point: aggregate counters cannot be cross-tabulated
against a list of precursors at all, whatever their magnitudes. The question
needs one terminal reason per precursor key.

That instrument now exists (`-out_terminal_reasons`, `d7349ed`) and the
attribution stays open until it has run.

## MEASURED: Gate C owns 98.2% of the candidate loss

`-out_terminal_reasons` on the s08_6x60 fixture, full library, against DIA-NN's
own run on the same fixture (4,948 confident precursors). Every library
precursor carries exactly one reason, so this is a count, not an argument.

    terminal reason   DIA-NN confident   share    whole library   share
             scored              3,859   78.0%          180,227    3.9%
             gate_c              1,069   21.6%        3,578,052   77.3%
         few_points                 14    0.3%          865,215   18.7%
       no_candidate                  6    0.1%            2,260    0.0%

Of the 1,089 that never reach candidate formation, **gate_c accounts for 1,069
-- 98.2%**. The picker itself fails 6 times.

The alternative both reviewers raised is REFUTED by the same table.
`not_reached` is ~13 across a 9.25M-precursor library and zero among DIA-NN's
confident set: no confident precursor is lost to missing isolation-window
coverage. `few_points` is large in the library (18.7%) and negligible in the
confident set (0.3%) -- that is the fixture's retention-time gaps catching
precursors that do not elute in the slices, which is what it should do.

So the attribution withdrawn above is restored, on evidence this time. Note the
fixture is EASIER than the full run at this stage -- 78.0% reach candidate
formation here against 58.9% on the full run -- so the share is a fixture
number and the full-run confirmation is still owed.

The whole-library column is the other half of the story: Gate C rejects **77.3%
of the library** and 3.9% of it is scored. That is what the gate is for; the
question has only ever been what it costs.

## A second correction: "ranking is not the problem" was too broad

Stage D measures within-precursor argmax. Stage E is a GLOBAL question --
whether the score separates targets from decoys across the run -- and it loses
7,138, the second-largest bucket. Both reviewers made this point independently.
The defensible claim is that per-precursor ORDERING is healthy; cross-precursor
DISCRIMINATION is not, and it is the second-largest loss in the funnel.

## What a perfect fix at A->B would be worth

`scripts/analysis/d5_yield.py`. Bin the precursors we find by DIA-NN abundance,
read acceptance off each bin, apply those rates to the abundance distribution of
the ones we lose:

    DIA-NN quantity        found   accepted     lost   predicted
        656 ..   8,661       775       7.0%    3,140         219
      8,661 ..  14,798     1,978      15.9%    3,894         620
     14,798 ..  29,371     5,069      30.4%    4,718       1,433
     29,371 ..  71,149     6,911      56.8%    2,876       1,633
     71,149 .. 191,728     4,843      75.0%    1,030         773
        >  191,728         3,482      81.1%      433         351

**At most +5,029, +40.9%** -- 12,292 -> 17,321 against DIA-NN's 39,149, moving
31.4% to 44.2%. It is an UPPER bound: within a bin, the precursors we currently
admit are the ones that passed the gate, so they are the better-behaved members
and the lost ones would score below their bin-mates. Assuming no abundance
dependence at all gives 8,578, which is 1.7x too high and is the number an
unstratified estimate would have produced.

So the biggest single lever in the funnel is worth about 1.4x against a 3.2x
gap. **Nothing found so far closes it, and the funnel does not contain a
missing 3x.**

The monotone acceptance column is the more useful reading: 7.0% of the faintest
precursors and 81.1% of the brightest, among ones we ALREADY find. DIA-NN
identifies all of them. Admission and discrimination are not two problems -- one
sensitivity deficit shows up at both.

## D->E is a separation problem, and it is not recoverable by fixing the null

`scripts/analysis/d6_separation.py`. The 7,138 correctly-ranked true peaks that
fail q <= 0.01, against the 626,295 decoy precursors' best scores:

    accepted at q<=0.01   DScore median  9.647   decoy percentile ~100.0
    REJECTED              DScore median  3.222   decoy percentile   98.8
    decoy null            DScore median -0.189   p99 3.412   max 10.078

The rejected true peaks outrank 98.8% of decoys -- and that is not close to
enough. 1.2% of 626,295 is **7,516 decoys above them**, against 7,138 targets
to be gained: accepting the bucket costs **51.3% FDP**. This is not a
conservative threshold, and no repair of the null moves it.

I nearly attributed the heavy tail to a candidate-count asymmetry -- the emitted
rows show targets 2.735 against decoys 2.957 per precursor, 8.1% more draws in
the decoy argmax. The code refutes it: `match_decoy_candidate_counts`
(`lda.h:1052`) already rank-matches the decoy draw counts against the target
distribution when the null is built, and does so one-sidedly (capping down,
never up), so it errs toward deflating the null. The row counts are the OUTPUT,
not the draws.

So the second-largest bucket is the LEAST tractable of the four. It needs
features the classifier does not have, and doc/45 measured that classifier
within +1.8% of its ceiling on the present feature set. Priority stays with
A->B -- now for a measured reason rather than because it is biggest.

## A measurement bug this turned up

DIA-NN writes `C(UniMod:4)`; our library writes `C(Carbamidomethyl)`. It is the
only modification on either side. Every analysis matching `Precursor.Id`
verbatim therefore dropped EVERY cysteine peptide -- 3,902 of 39,149, 10.0% --
from numerator and denominator alike. Aliasing takes stage A from 90.0% to
100.0%. Fixed in `d3_stages.py` and in `scripts/fdp_compare.py`, whose DIA-NN
concordance column had it too. Codex asks for a collision audit of the
normalised key space, which is not yet done.

## Open

* `-out_terminal_reasons` on the fixture, then on the full run: which stage
  actually owns the 16,091.
* Codex's two oracles, which separate admission from discrimination without
  guessing: a SELECTION oracle (pick the candidate nearest DIA-NN's RT, score
  unchanged) and an ADMISSION oracle (inject a candidate at DIA-NN's RT for the
  missing precursors, score it normally). If the admission oracle produces
  candidates that then fail q, admission is not sufficient -- which is exactly
  what d5_yield.py predicts quantitatively.
* The faintness signature is an association measured with DIA-NN's own
  abundance estimate on peptides DIA-NN identified. An ODIA-independent
  summed-fragment measure at DIA-NN's RT would be the non-circular version.
* |dIM| was computed against the library's PREDICTED 1/K0, not the calibrated
  centre ODIA searches. Equality between the groups is suggestive, not decisive.
