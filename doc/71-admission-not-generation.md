# Admission, not generation: the correct candidate is present and refused

doc/70 established that ODIA's peak PICKING is sound -- where it accepts a
precursor its apex matches DIA-NN to a median 0.00 s -- and that the gap is
whatever it declines. This locates that gap one step further.

Measured 2026-08-28 on the 10,000-precursor stratified library, full
`S08_diaPASEF.mzpeak`, fixed linear RT map, `-precursor_im_window 0.025`,
`-rt_window 300`, `-passes 1 -gate_alpha 0`, one pinned binary
(`8a48387cca4537eb`), scored against DIA-NN's **reported** RT.

Of the 4,722 DIA-NN-confident precursors ODIA scores, it accepts 2,376 at
q <= 0.01 and declines 2,346. **For 39.3% of those it declines, a candidate
within 10 s of DIA-NN's answer is already in its own candidate list.**

Raising candidate generation confirms which half is binding:

    arm            cand/prec  accepted   delta  apex<=10s  PROPOSED  ADMITTED  conversion
    baseline            2.22     2,376      +0      95.0%     39.3%         -           -
    -max_candidates 10  5.32     2,489    +113      95.3%     58.7%      5.8%        9.9%
    -cand_min_sep 5     2.19     2,441     +65      95.2%     39.6%      4.7%       11.9%
    both                5.18     2,502    +126      95.2%     60.4%      6.6%       11.0%

PROPOSED and ADMITTED are both over the baseline's 2,346 declined precursors:
the fraction for which a correct candidate exists anywhere in the list, and the
fraction that end up accepted at q <= 0.01 with a correct apex. Conversion is
the ratio.

**Generation responds and admission does not.** `-max_candidates 10` lifts the
proposed rate by 19.4 points, from 39.3% to 58.7%. Admission moves from 0 to
5.8%. **Conversion sits at 9.9-11.9% in every arm**: when ODIA has the correct
peak in its hand, it takes it about one time in nine.

That is the deficiency, and it is not candidate generation, not peak selection
(apex accuracy is unchanged at 95.0-95.3% across every arm, so the extra
candidates cost nothing), and not the ranking quality of DScore, which doc/69's
review established beats its own best single input where it counts.

## What the extra candidates are worth anyway

`-max_candidates 10 -candidate_min_separation 5` is +126 accepted precursors on
this library, 2,376 -> 2,502 (+5.3%), with apex accuracy flat. It is a free
modest win and it should be taken, but it addresses a fifth of the gap at most.
The option text at `OpenDIAlyzer.cpp:1024` predicted the generation half of this
("recall of the correct position within the cap rises 72.7% to 75.5%" at
separation 5) and left the classifier half open; the classifier half is where
the loss is.

## Where this points

The remaining question is why a correct, proposed candidate fails to reach
q <= 0.01. Three things it is NOT, all measured: not the apex (0.00 s median
where accepted), not the trace (r = 0.995 against a 0.972 wrong-fragment null,
doc/68 lineage), not candidate availability (39-60% proposed). It is the
admission decision itself -- the score the candidate receives relative to the
decoy null, and the threshold applied to it.

That is the same place the project's earlier accounting landed independently:
Gate C carrying ~98.8% of the emission loss. Two different routes, one target.

---

# RETRACTED: the baseline was not a control

Adversarial review, 2026-08-28, confirmed at source. Everything above that
compares an arm to "baseline" is a TWO-factor comparison reported as one.

`shared/libv2/run_random_extract.sh` passes `-out_chrom`; `run_admit.sh` does
not. `src/OpenDIAlyzer.cpp:1864` forks on `if (out_chrom.empty())` into a
different scoring path, and the two paths do not score the same feature set:

    baseline (acq_i025, -out_chrom)      dropping 7 sub-score(s)
                                         Mass.Ppm NaN on 19,999 / 19,999 rows
    cap10 / sep5 / both (no -out_chrom)  dropping 4 sub-score(s)
                                         Mass.Ppm NaN on 0 / 19,999 rows

**The baseline ran a classifier with three fewer live features than every arm it
was compared against** (`var_im_delta`, `var_mass_accuracy`, `var_mass_spread`).
So the generation response (39.3% -> 58.7%), the threshold clearance column, and
the +126 accepted are all confounded with a feature-set change.

Three further corrections from the same review, independent of the fork:

* **The 82% refusal rate is conditioned on the outcome it indicts.** "Declined"
  is DEFINED as q > 0.01, so measuring refusal within that set is circular. Over
  all 4,731 DIA-NN-confident precursors, ODIA admits **76.2%** of its correct
  top-ranked candidates at q <= 0.01, in every arm (75.7-77.5%). The 15-18%
  figure is a flip rate between two deterministic arms and is withdrawn.

* **+/-10 s is not a correct-peak criterion.** It is 5.6x DIA-NN's median FWHM
  (3.57 s) and 2.4x ODIA's own median peak-group width (8.31 s). Against a
  decoy-candidate null the "has a correct candidate" rate falls from 39.3% to
  32.4% (baseline) and 60.4% to 43.5% (both), so the generation response is
  +11.1 real points rather than +19.4.

* **Generation loses more than admission, not less.** In absolute counts over
  the 2,353 declined, with the decoy null applied: generation 1,590 (67.6%),
  ranking 121 (5.1%), admission 642 (27.3%). A1 is refuted by its own numbers
  once they are expressed as losses rather than rates.

What survives: **ranking is not the constraint.** Winner-level accuracy is
30.0% (baseline) to 36.5% (both) against a decoy-winner null of 3.5-3.7% --
8-18x enrichment -- and it RISES with more candidates, so the 76.2% -> 62.9%
"ranking degradation" cited above is denominator inflation, not degradation.

`adm_scores_base.tsv` re-runs the baseline configuration without `-out_chrom`
so the comparison can be redone on one factor.

---

# Redone on one factor: the title of this document is wrong

The matched baseline (`adm_scores_base.tsv`, same configuration WITHOUT
`-out_chrom`, so all four arms drop the same 4 sub-scores and carry `Mass.Ppm`
populated) gives:

    arm     cand/prec  accept q<=.01  matched decoy budget  top-cand RT-ok  admitted of those
    base         2.22          2,405          2,391    (+0)          3,000             76.3%
    cap10        5.32          2,489          2,405   (+14)          3,130             75.8%
    sep5         2.19          2,441          2,417   (+26)          2,996             77.6%
    both         5.18          2,502          2,403   (+12)          3,122             76.3%

Matched decoy budget = the 26 decoy groups above the baseline's own acceptance
cut (DScore 4.3234), so every arm is read at the same realised false-positive
load rather than at its own q-value.

**The generation gain is mostly bought.** +126 against the broken baseline
becomes +97 against a matched one and **+12 to +26 at matched decoy budget**.
`sep5` is the best of the three at +26 and proposes FEWER candidates than
baseline (2.19 vs 2.22), which is the opposite of a generation effect.

**Admission is not the defect.** ODIA admits **76.3%** of the precursors whose
top candidate sits within 10 s of DIA-NN's reported RT, and that rate is flat at
75.8-77.6% across all four arms. The 15-18% clearance in the retracted section
was conditioned on the outcome it indicted.

**The gap is that the top candidate is often not the right peak.** Of 4,731
DIA-NN-confident precursors, ODIA's top-DScore candidate is RT-correct for
**3,000 (63.4%)**. It then admits 76.3% of those and declines most of the rest,
which is conservative and correct behaviour rather than a refusal of good
evidence. The loss is the 36.6% where ODIA's best candidate is somewhere else.

That is not admission, not the q-value cut, not candidate count, and not
ranking-among-what-was-proposed (winner-level enrichment is 8-18x over a decoy
null). It is which peak group ends up on top, on precursors where a correct one
may or may not have been generated at all -- and separating those two requires a
correctness criterion better than "RT within 10 s", which is 5.6x DIA-NN's
median FWHM.

---

# Resolved with a peak-identity criterion: it IS generation

The +/-10 s criterion used above is 5.6x DIA-NN's median FWHM and carries a
decoy null of 10-30%, which is why every decomposition built on it moved when
poked. Replacing it with **boundary IoU >= 0.5 against DIA-NN's own
`RT.Start`/`RT.Stop`** (median peak width 9.70 s, p95 15.24 s) drops the decoy
null to **3.4%** and the picture stops moving:

    arm     top-cand correct   ANY cand correct   decoy null   admitted of top
    base           52.5%             56.2%            3.4%          78.0%
    cap10          54.8%             62.4%            9.2%          77.4%
    sep5           52.2%             55.9%            3.4%          79.5%
    both           54.4%             62.3%           10.2%          78.3%

Over 4,731 DIA-NN-confident precursors:

* **Generation loses 43.8%.** No ODIA candidate overlaps DIA-NN's peak at all.
* **Ranking loses 3.7 points** (56.2% -> 52.5%), i.e. 6.6% of what is generated.
* **Admission keeps 78.0%** of what reaches the top.

So the ordering is generation >> admission > ranking, which is the review's
position and the opposite of this document's original title. The title is wrong
and is left standing only so the retraction is legible.

**Raising the candidate cap does not fix it.** Null-corrected -- (obs - null) /
(1 - null) -- `cap10` moves ANY-candidate 54.7% -> 58.6% but TOP-candidate
50.8% -> **50.2%**. More candidates are generated, none of them win. That is
exactly why the matched-decoy-budget gain was +12 to +26 rather than +126: the
extra candidates are mostly coincidences, which the rising decoy null (3.4% ->
9.2%) states directly.

## What this means for where to work

The 43.8% is not a threshold, not a ranking, and not a cap. For those
precursors ODIA never proposes a peak group overlapping the one DIA-NN reports,
inside a +/-300 s window that provably contains DIA-NN's answer for 99.8% of
them. The candidates it does propose are elsewhere.

That points at the peak DETECTOR -- what counts as a local maximum worth
emitting -- rather than at how many are kept or how they are scored. The
co-elution evidence used to seed candidates is the natural suspect, and it is
testable without a new run: for a sample of the 43.8%, the traces are already on
disk and DIA-NN's bounds are known, so the question "was there a detectable
peak there at all, by any criterion" can be answered directly from the
chromatograms.
