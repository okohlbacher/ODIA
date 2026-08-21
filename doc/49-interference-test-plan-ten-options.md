# Ten options for interference handling, as a test plan

2026-08-21. Drawn from doc/48 (research) after reviews 48 and 49. Three are
DELETED before implementation on the reviewers' arguments; the rest are ordered
so that each is only worth running once the one before it has landed.

## The gate that sits in front of all ten

**Nothing here runs until the decoy model is settled** (`full_v8`, shuffle vs
mutate, in flight). Both reviewers, independently:

  * the ruler is 109x wrong at depth 2,000 and 1.4x at depth 20,000, so any
    result measured now is measured against a null of unknown shape;
  * deconvolution's decoy-symmetry properties DEPEND on the decoy model --
    mutate decoys at unchanged precursor m/z collide with real fragments
    differently from composition-preserving shuffle decoys, so a method tuned
    under mutate may not transfer.

And the 42.5% no-candidate bucket is larger than the 18.5% interference bucket
and upstream of it: a precursor with no candidate cannot be rescued by any
amount of cleaning.

## Four controls every option must satisfy

Pre-registered, from review 49. An option that cannot satisfy these is not
implementable regardless of its merits.

  1. **Label-blind construction.** Any active set / dictionary / competitor list
     must be built by a rule that cannot see the target-decoy label.
     `PeakGroupScorer.cpp:1586` already enforces "targets and decoys traverse
     identical code"; a new layer must inherit that mechanically.
  2. **Entrapment-matched acceptance, pre-registered.** IDs up at matched
     empirical FDP. Nominal q is 5-7x off and must not appear in the evaluation.
     The decoy arm alone CANNOT validate this: the failure mode is signal
     assigned to an absent hypothesis, which is what entrapment represents and
     decoys represent poorly.
  3. **Replicated tail statistics.** Decoys above the operating threshold and
     decoy p99, before and after, across >= 3 independent decoy realisations.
     The threshold is set by ~70 tail events with ~12% counting error.
  4. **Fit quality is a diagnostic, never a feature.** Per-scan explained
     fraction as a report of where the dictionary is lying. Scoring on it is
     selection-on-the-response again.

---

## DELETED before implementation

**D1. Per-scan full-library NNLS (Specter as published).** ~400,000 columns per
isolation window against far fewer observed fragment bins per scan. rank(L) <=
bins << columns, so the solution is non-unique before it is expensive.
Infeasible as stated.

**D2. Learned denoiser or scorer fed back into target-decoy validation
(Alpha-XIC shape).** Trained on peak groups from the run it is applied to, score
returned to the same validation loop. On a ruler that is 109x off in the tail, a
transform that lifts targets and decoys together reads as a gain. This is the
same evaluation loop in which our prominence gate read -17% at nominal q and
+12% at matched FDP. Unvalidatable here.

**D3. Two-pass design whose competitor set comes from pass-1 IDENTIFICATIONS.**
This was the backlog's surviving direction and it is label-asymmetric by
construction: pass-1 IDs are overwhelmingly targets, so targets get cleaned and
decoys keep their contamination. It would widen the decoy/entrapment gap, not
close it.

---

## The census that must precede any unmixing work

**O1. Dictionary census and identifiability audit. LIBRARY HALF DONE.**

Measured on human_v2 (9,617,705 precursors, 93,010,926 fragments), S08's 24
windows over 300-1800 Th, RT gating via the fitted map at t = 1000 s:

    window      m/z        columns    +-60s RT    +-10s RT
         0  300-362        895,615      97,857      15,551
         6  675-738        551,079      61,929       9,866
        12 1050-1112       304,905      27,682       4,540
        18 1425-1488       155,916      17,985       2,938

    mean columns per window: 400,738    max: 925,529

codex predicted ~400,000 from the arithmetic; measured 400,738. Even an
aggressive +-10 s retention-time gate leaves **2,938 columns in the sparsest
window and 15,551 in the densest**, and the gate itself is a selection rule that
can delete the true explanation.

**Still needed: the observation count.** rank(L) is bounded by the number of
distinct fragment m/z positions carrying signal in one scan, and
`odia_run_info` does not report it -- it needs a small instrumented run. Until
that number exists the verdict is only half-formed, but the shape is already
visible: the system is plausibly determined only in the HIGH-m/z windows, and
the low-m/z windows would need gating so tight that selection bias becomes the
dominant concern. O2-O4 may turn out to be viable over part of the m/z range
only, which is a materially different proposition from "Specter for ODIA".

Original description: No code, one pass over the
library plus one over a fixture run. Per scan, measure: columns after
precursor-window filtering; columns after RT and fragment-evidence gating;
observed non-zero fragment bins; effective rank / mutual coherence of the
resulting L; and solution stability under small perturbations and library
subsampling.

Codex's position, which I accept: **without these numbers "Specter-like NNLS" is
not a specified algorithm for this repository.** If gating cannot bring columns
below observed bins with margin, options O2-O4 are all dead and we should know
that before building any of them.

Cost: hours. Fixture-safe (this is a mechanism measurement, not an FDP one).

---

## Options, in dependency order

**O2. Candidate-bounded NNLS with a label-blind active set.** Specter's
S = Lc + N, but the dictionary is restricted to precursors with a candidate peak
group in this RT neighbourhood, chosen by a label-blind evidence rule so decoys
enter on the same terms. Report the coefficient trace for the query precursor.
Depends on O1 showing the restricted system is determined.

**O3. Coefficient-as-feature, not coefficient-as-trace.** Rather than replacing
traces (which would let one self-fitted transform improve several nominally
separate score columns at once -- codex, review 48), add the NNLS coefficient
profile's shape statistics as NEW features alongside the existing ones. Cheaper,
strictly additive, and the winner's-curse interaction is confined to one column.
Depends on O2.

**O4. Explained-fraction diagnostic per window.** From the same fit: how much of
each scan the dictionary explains. A report of dictionary completeness, which
doc/48's correction identifies as the binding constraint on the whole family.
Never a feature (control 4). Depends on O2.

**O5. Shared-fragment arbitration. TESTED -- REFUTED by saturation, like the two before it.**

kimi's distinction was that a competitor should only count if it has its OWN
EVIDENCE in this run -- a candidate peak group nearby -- rather than merely
existing in the library. Implemented label-blind (candidacy is open to targets
and decoys alike) over v7's 775,183 evidenced precursors, 8,381,732 fragment
entries, competitor = same fragment m/z within 10 ppm whose own apex is within
10 s:

    mean evidence-clean fragments      0.4 of 11.6
    corr over ALL fragments            AUC 0.767
    corr over EVIDENCE-clean           AUC 0.841  -- on 100 of 16,000 precursors

The 0.841 is NOT a result. It is computed on the 100 precursors that retained
four or more clean fragments, which is the same selected-minority trap that
produced the +0.021 at a 60 s competitor tolerance and evaporated at a correct
one. **11.2 of 11.6 fragments are contested even when the competitor must have
its own evidence.**

That is the third refutation of the same shape (library competitors, evidenced
competitors, and the pruning that preceded both), and together they say
something structural rather than incidental: **at 9.6M library precursors with
775k of them evidenced, every BINARY "is this fragment contested" test returns
yes.** Saturation is the regime, not an artefact of one tolerance.

The consequence is the useful part. Arbitration cannot be binary -- keep or
drop. It has to be QUANTITATIVE: how much of this fragment's intensity is mine.
That is deconvolution, and it is gated on O1's census showing the restricted
system is determined. The cheap non-dictionary route out of the interference
bucket is now closed; what remains all imports outside information at a cost.

Original description: When two library
precursors in one window share a fragment m/z within tolerance and both have
candidates at overlapping RT, arbitrate that fragment's intensity between them
by their OTHER, unshared fragments. This is DIA-NN's documented mechanism, needs
no matrix factorisation, and attacks the measured mechanism directly (twelve
fragments each carrying a different interferent's timing). kimi keeps this one.
Independent of O1-O4.

**O6. Rank-1 residual as a feature. TESTED -- weak, do not build yet.**

Measured offline on the 8,000-precursor diagnostic dump (7,991 target / 7,991
decoy), each fragment scaled to unit norm so the leading component measures
SHAPE agreement rather than brightness:

    feature                        AUC   med target   med decoy
    rank-1 energy share          0.712        0.660       0.587
    mean pairwise corr (ours)    0.767        0.096       0.023
    Spearman(rank-1, pairwise) = 0.488

Partly independent (rho 0.49) but WEAKER than the feature it would join. That is
now the second structurally-motivated alternative to come in below the existing
pairwise correlation -- apex dispersion was 0.716 against xcorr_coelution's
0.810. **Two different functionals of the same fragment x RT matrix both lose to
the one we already compute**, which suggests the corr_sum / xcorr family is
close to what that matrix supports and that new features must bring information
from OUTSIDE it (a dictionary, MS1, or another pass) to help.

Not dead: the decisive test is INCREMENTAL -- does adding it raise the
supervised ceiling of +1.8%? -- and that needs a representative population,
which the diagnostic set (50% hard negatives by construction) is not. Run it
after O8, on a full-run feature matrix.

Original description:

Fit the best rank-1 approximation Fit the best rank-1 approximation to the
precursor's own fragment x RT matrix and score the RESIDUAL energy, not the
fitted profile. A real precursor is rank-1 plus noise; a contaminated one is
not. No dictionary, no naming problem, no cross-precursor cost. Must be checked
for collinearity with `xcorr_shape` and `corr_sum` before being believed --
apex-dispersion looked promising on the same reasoning and came in at AUC 0.716
against xcorr_coelution's 0.810.

**O7. MS1-seeded detection, Siren-shaped.** Siren's verified 1,711 -> 2,516 is an
MS1-level result on raw MS2 and belongs to our 42.5% DETECTION bucket, not the
18.5% interference bucket -- both reviewers say the research document filed it
wrongly. Use the calibrated MS1 traces (6b1b2fa) to propose candidates where the
fragment picker proposes none. Treat Siren's number as an existence proof, not
an effect size.

**O8. Composition-preserving decoys as the shipped default.** Not an
interference method, but the prerequisite for evaluating any of them, and `v8`
is already measuring it. If shuffle flattens the 109x tail disagreement, it
lands regardless of what happens to O1-O7.

**O9. Mild mass tightening as an ADDITIONAL feature.** The one surviving
fragment of the tightening idea: ~6 ppm (above our 4.04 ppm per-hit scatter),
kept ALONGSIDE the 10-12 ppm sum rather than replacing it, as DIA-NN does at
base / 0.45x / 0.2x. Cheap, fixture-testable, and the failure at 3.6 ppm is
explained by width rather than principle.

**O10. Entrapment-calibrated evaluation harness. DONE -- `scripts/fdp_compare.py`.**

    scripts/fdp_compare.py gateC=rank_a.tsv prominence=rank_b.tsv

Reports identifications at a grid of MATCHED empirical FDP targets, with the
per-cell entrapment count and the FDP sigma beside every number, plus DIA-NN
concordance at each depth. It computes r on the library each arm ACTUALLY
searched and prints which rule it used. Verified against doc/46: it returns
+11.7% at 7.42%, +5.6% at 5.72%, +17.6% at 10%, reproducing the Gate C result
that nominal q had reported as -17%.

It also surfaces what the ad-hoc analyses hid -- at 2-3% FDP the cells rest on
3/0 and 9/0 entrapment hits with sigma 2.3-3.2 pp, so the apparent -86% and -92%
there are noise, and the tool says so.

Original description:

Make control 2 mechanical:
one command that takes two run outputs and reports IDs at matched empirical FDP
with replicated decoy tails, so no option above can be accepted on nominal q.
Given that nominal q is 5-7x off and that the fixture attenuates FDP effects
~5x, this is the piece that decides whether any of the other nine can be
believed at all.

## Suggested order

    O8  (in flight)  ->  O10  ->  O1  ->  O5, O6, O9 in parallel  ->  O2 -> O3, O4
                                          O7 alongside, different bucket

O5, O6 and O9 are the cheap ones and none needs the dictionary machinery. O2-O4
are gated on O1 returning a determined system. O7 is the only one aimed at the
larger bucket.
