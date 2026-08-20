# Where the gap to DIA-NN actually is, after the extraction question closed

2026-08-20, S08 diaPASEF, our 9.6M library, DIA-NN run on the SAME library.

## Three questions, and only one of them is still open

**Do we extract as cleanly?** Essentially yes. Identity-matched fragments,
common pair set, at PRODUCTION settings (10 ppm, +/-0.025, mobility seed
applied): ODIA 0.270 against DIA-NN 0.311 mean pairwise fragment coherence,
median paired difference -0.010, and ODIA is better on 41.6% of precursors.
doc/44 has the full control set -- smoothing, apex placement, fragment identity,
aggregation, pointwise cell-by-cell agreement. The earlier "0.167 against 0.310"
was an artefact of dump flags production does not use.

**Is our FDR calibrated?** Now, roughly. Entrapment FDP 5.59% at nominal 1%
(doc/43), against DIA-NN's 7.42% on this library. We are STRICTER than DIA-NN,
not looser, which is the opposite of the position we were in twelve hours ago at
12.63%.

**Do we find as many?** No, and this is the whole remaining gap. Walking v3's
ranked list to matched empirical FDP, with the post-filter entrapment ratio
r = 0.17309:

      target FDP   ODIA IDs   DIA-NN confident recovered
           5.59%     13,454                       11,696
           7.42%     15,036                       12,755   <- DIA-NN's point
          10.00%     16,193                       13,398
          20.00%     19,855                       14,556
          30.00%     23,972                       15,255

    DIA-NN: 39,211 IDs at 7.42%

**2.6x at DIA-NN's own operating point**, and we do not close it even by
accepting a 30% false discovery proportion. That is not a threshold that can be
moved. The identifications are not in the ranked list to be admitted.

## So the deficit is upstream of scoring, and we already have a number for part of it

On a truth-only library, Gate C alone rejected 10,727 of 33,749 precursors
DIA-NN confidently identifies -- 31.8% -- and it was the ONLY material
pre-scoring loss (doc/43: 0 lost to window coverage, 0 to "no points"). That
caps recovery at 68.2%, i.e. ~23,000 of the truth set, before scoring is even
reached.

We recover 12,755 at DIA-NN's operating point. So Gate C explains roughly a
third of the shortfall and something else explains the remaining ~2x.

## What that leaves, in the order the evidence supports

1. **Gate C's 31.8% true-precursor rejection.** Measured, and the largest single
   named loss. The prominence replacement was tried and was catastrophically
   worse (doc/43: 70.2% FDP) -- but that failure is now known to have been the
   degenerate library correlation and the missing decoys amplifying it, both
   fixed since. **Gate C's replacement has NOT been re-tested since the library
   became label-symmetric.** That is the single highest-value open experiment.
2. **Candidate detection sensitivity.** DIA-NN's picker is co-elution-first over
   per-transition maxima; ours picks from a summed standardised trace, which one
   loud interferent dominates. codex's per-transition union primitive (doc/42)
   was designed for exactly this and never built.
3. **Iteration.** DIA-NN refits RT, mass and mobility across ~12 passes against
   retained chromatograms; we do two and re-extract from the raw file for each.
   The plan file in this repo already scopes retaining peak groups and iterating
   scoring, and notes the candidate picker is RT-agnostic so refitting the map
   costs one column, not a re-extraction.
4. **Scoring headroom.** Sub-score AUCs on v3 against DIA-NN-confident targets:
   corr_sum 0.909, fragment_coverage 0.825, xcorr_coelution 0.810, xcorr_shape
   0.764, ms1_coelution 0.757 (now live, AUC earns its 47 GB), library_rmsd
   0.750. Dead: log_sn 0.501 (third confirmation), candidate_margin 0.487.
   These are measured with DIA-NN-confident as positives and therefore carry
   that selection bias -- usable to condemn a dead feature, not to rank live ones.

## What is NOT the problem, so we stop paying attention to it

- Trace quality, per above.
- m/z tolerance. DIA-NN's own log records an optimised MS2 accuracy of 12 ppm
  against our 10; it extracts through a WIDER window than we do.
- Aggregation. `Max` loses to `Sum` at every width tested.
- Mobility window width. At a corrected centre, +/-0.025 beats +/-0.0304 and
  +/-0.050; widening only ever compensated for a centring error.

## The sensitivity gap splits in half: detection and ranking, ~equally

Partitioning DIA-NN's 33,749 confident precursors through v3's pipeline:

                                 stage    count   % of truth
              not in our library at all        0         0.0%
      dropped by -min_library_fragments      168         0.5%
      eligible, but NO candidate formed   11,655        34.5%
         candidate formed, not accepted   10,281        30.5%
                    ACCEPTED at q<=0.01   11,645        34.5%

Three things follow.

**Library coverage is not a factor.** Zero of DIA-NN's confident precursors are
absent from our library. Whatever else is wrong, we are searching for the right
things.

**The fragment floor costs 0.5%.** 168 precursors, against halving the FDP.
That is the whole price of doc/43's fix measured against the truth set.

**The remaining loss is two roughly equal halves.** 34.5% never form a candidate
at all -- a DETECTION loss, and the independently measured Gate C rejection of
31.8% on a truth-only library accounts for nearly all of it. 30.5% do form a
candidate and are then not accepted -- a RANKING loss, which no gate change can
touch.

That is the argument for doing both, and against expecting either alone to
close a 2.6x gap. If Gate C were perfect and nothing else changed, the ceiling
is 21,926 of 33,749 (65%) -- still well short of DIA-NN, because the 10,281
detected-but-unranked would remain unconverted. Conversely, a perfect scorer on
today's candidate set cannot exceed that same 21,926.

Ordering, on this evidence: the gate first, because it is one flag on an
experiment already queued (`run_full_v6.sh`) and because its 31.8% is the single
largest named loss anywhere in the pipeline; then the ranking half, where the
sub-score audit says corr_sum at AUC 0.909 is carrying the classifier almost
alone and every other live feature is between 0.62 and 0.83.

## The definitive partition, against DIA-NN's own 39,149 on OUR library

The earlier table used `truth_ids.txt` (33,749), which is a DIA-NN run on a
DIFFERENT library and overlaps ours on only 29,275. `dn_xic_report.parquet` is
the run on OUR library -- 39,149 precursors at q <= 1% -- and is the right
reference. Matching each to v3's candidate list, with a 15 s tolerance on
DIA-NN's reported apex:

    bucket                                        count       %   cause
    accepted by ODIA                             11,773   30.1%   --
    right peak, correctly picked, scored too low  7,229   18.5%   SCORING
    right peak present, wrong one picked            730    1.9%   picking
    candidate formed, none near DIA-NN's RT       2,777    7.1%   picking/detection
    no candidate formed at all                   16,640   42.5%   DETECTION

Of the 10,736 rejected precursors for which we DID form a candidate, 7,959
(74.1%) have one within 15 s of DIA-NN's apex, and in 7,229 of those (90.8%) it
is already the best-scoring candidate we hold. Sanity: 96.4% of ACCEPTED
precursors have their best candidate at DIA-NN's RT, so the measure behaves.

**7,229 precursors are located correctly, picked correctly, and not promoted.**
Nothing upstream is at fault for them. That is 18.5% of DIA-NN's set sitting in
our own output at the right retention time, and it is the cleanest scoring
target we have ever had -- no gate change, no picker change and no extraction
change can reach it.

What they look like (medians, against accepted truth and decoys):

                            accepted   rejected   decoy
    var_library_corr           0.734      0.040   0.000
    var_xcorr_shape            0.599      0.227   0.211
    var_corr_sum               8.482      4.279   2.805
    var_ms1_coelution          0.991      0.028  -0.107
    var_peak_width_ratio       0.355      1.000   1.000
    var_rt_delta              13.107     19.133  29.911
    var_im_delta               0.008      0.009   0.024

Decoy-like on every co-elution and shape feature; POSITION-wise as good as the
accepted set (`im_delta` 0.009 against 0.008) and clearly better than decoys
(`rt_delta` 19.1 against 29.9). We put the peptide in the right place and then
cannot tell its trace from noise. DIA-NN can.

Note also that `var_log_sn` is 4.605 for all three groups and
`var_usable_fragments` is 12.000 for all three -- two of nineteen features are
constants across the exact comparison they exist to make.

## Revised ordering

The detection bucket is still the largest (42.5%) and `run_full_v6` is testing
it. But the 18.5% scoring bucket is the better-defined problem: the inputs are
already correct, so any improvement is attributable, and it needs no 3 h run to
iterate -- the candidates and their sub-scores are already on disk in
`sub_v3.tsv`. A classifier experiment can be run offline against the exact rows
that fail.

## The features are the ceiling, not the classifier

The 7,229 correctly-located, correctly-picked, unpromoted precursors invite an
obvious response: train a better scorer. Measured, that response is wrong.

Trained on v3's own sub-score matrix, positives = DIA-NN's 22,509 confident
precursors on our library, negatives = 120,000 sampled decoys, 5-fold
cross-validated so nothing is fitted on its own test rows:

                                      AUC   promoted at decoy FDR 1%
    production DScore                0.9156                   13,609
    CV supervised GBT, same features 0.9388                   13,848
                                                     headroom +239 (+1.8%)

The supervised model is handed DIA-NN's answer key as its training labels --
an advantage production cannot have, since its own trainer is semi-supervised
on target/decoy alone. Even so it recovers 1.8%. That is an UPPER BOUND on
everything classifier-side: model family, hyperparameters, training procedure,
the semi-supervised seed, all of it. There is nothing there.

So the 18.5% bucket is not reachable by scoring as currently defined. What is
missing is INPUTS.

The sub-score profile says which kind. Those precursors are decoy-like on every
co-elution and shape feature and positionally correct on rt_delta and im_delta:
we put the peptide in the right place and its fragments do not co-elute. For a
weak precursor sharing an isolation window with an abundant one, that is what
contaminated fragments look like -- and every co-elution feature we compute uses
ALL of a precursor's fragments, contaminated ones included. One interfered
fragment out of twelve drags corr_sum and xcorr_shape toward the decoy
distribution regardless of how clean the other eleven are.

DIA-NN has explicit interference correction (vault: *DIA-NN interference
correction*). ODIA has none. That is the single named capability difference left
standing after tonight, now that extraction cleanliness, mass tolerance,
aggregation, mobility window width and FDR calibration have all been excluded
by measurement.

The direction, therefore, is per-fragment interference handling: identify which
fragments of a precursor are contaminated and exclude or down-weight them before
the co-elution features are computed, rather than after. Candidate primitives
already discussed and never built: codex's per-transition union picking
(doc/42), and `var_im_spread`, which is built, correct and shipped OFF because
it cost 60 identifications when the FDP was 12.63% and unmeasurable -- a
judgement made against a purity signal we did not have then and do have now.

## Naive interference pruning is refuted before implementation

The obvious form of the fix above -- drop a precursor's contaminated fragments,
then compute the co-elution features on what remains -- was prototyped offline
on the 1,000-precursor seeded dump (998 targets, 998 decoys, real traces, +/-30 s
around the apex, fragments dropped greedily by lowest mean correlation to the
rest until half remain):

    variant                                AUC   median target   median decoy
    corr over ALL fragments (current)   0.8288           0.150          0.021
    corr over consistent subset         0.8278           0.434          0.182
    fraction consistent                 0.7465           0.083          0.000

Pruning nearly triples the target median -- and triples the decoy median with
it. Discrimination does not move (0.8288 -> 0.8278). The reason is that the
pruning criterion IS the scored statistic: selecting a decoy's best-correlating
fragments flatters the decoy exactly as much as it flatters a target. This is
the same selection-on-the-response failure that invalidated the doc/42 sub-score
audit, in a new place.

So the correction cannot be per-precursor. It has to use something a decoy does
not share, and the vault says what that is: identification-level interference
correction requires SHARED-FRAGMENT logic -- knowing which OTHER library
precursor owns the contaminating signal, not merely that a fragment disagrees
with its siblings. A decoy's fragments are recomputed from a shuffled sequence,
so they do not systematically collide with a real co-eluting peptide's fragments
in the way a real precursor's do; that asymmetry is the only place the signal
can come from.

That is a substantial build -- a cross-precursor index over the library's
fragment m/z within each isolation window, evaluated at the candidate's RT --
and it should be designed before it is started. It is, however, the only
capability difference against DIA-NN still standing after tonight.

## What tonight excluded, by measurement

    trace cleanliness      0.270 against DIA-NN's 0.311, matched
    m/z tolerance          DIA-NN optimises to 12 ppm; we extract at 10
    aggregation            Max loses to Sum at every width tested
    mobility window width  +/-0.025 beats +/-0.0304 and +/-0.050 at a correct centre
    FDR calibration        5.56% against DIA-NN's 7.42%; we are stricter
    library coverage       0 of DIA-NN's confident precursors are absent
    the classifier         +1.8% headroom with DIA-NN's own labels
    per-precursor pruning  no gain; decoys benefit equally

What remains: the detection bucket (42.5%, `run_full_v6` testing it now) and
cross-precursor interference correction.

## The library-competitor design is refuted too: at 9.6M precursors, everything is contested

Second prototype, before building. For each fragment of each of 998 targets and
998 decoys, count OTHER library precursors in the same isolation window
(+/-12.5 Th) with a fragment within 10 ppm whose predicted RT maps within 60 s
of the candidate's apex and whose predicted 1/K0 is within 0.025. Then score
co-elution over the UNCONTESTED fragments only.

    contested fraction, median          1.000 target   1.000 decoy   AUC 0.603
    mean uncontested / total fragments  1.9 / 11.7
    precursors with >= 4 uncontested    312 of 2,000
    on those 312: corr ALL              AUC 0.8574
                  corr UNCONTESTED      AUC 0.8788

**CORRECTION (codex, review 46):** the +/-60 s above is the pass-2 EXTRACTION
window, which is a resource bound, not evidence that a competitor eluted at the
apex. The right test is peak-width scale. Re-run properly:

    RT tol   contested t/d   uncontested/total   corr ALL   corr UNCONTESTED     n
      60 s   1.000 / 1.000          1.9 / 11.7     0.8574             0.8788   312
      30 s             --           2.7 / 11.7     0.8507             0.8575   474
      10 s   0.667 / 0.750          4.6 / 11.7     0.8496             0.8287   986
       5 s   0.500 / 0.500          6.1 / 11.7     0.8402             0.8298 1,392

Tightening the tolerance does cut the contested fraction from 1.000 to 0.500 --
codex was right about the test. But as the tolerance becomes correct and the
sample becomes REPRESENTATIVE, dropping contested fragments makes discrimination
WORSE: 0.8287 against 0.8496 at 10 s, 0.8298 against 0.8402 at 5 s. The +0.021
"gain" recorded below existed only on the 312-precursor minority that survived a
wrong filter -- a small-n artifact of the same family this document spends its
length catching, produced here by me.

Note also that at 10 s DECOYS are more contested than targets (0.750 against
0.667), which is the class-construction hazard codex predicted: the feature
would partly encode how decoy fragments are generated rather than interference.

The original text, now superseded:

Every fragment is contested, in both classes. The m/z axis is saturated at this
library size -- doc/19 recorded exactly that ("5M precursors over ~380-980 Th
saturate the m/z axis and no unsupervised matching statistic survives its own
control"), and this design was drawn without connecting it. The +0.021 AUC on
the 312-precursor remainder is real, is the same size as the classifier
headroom, and is measured on a selected 16% minority.

So a BINARY contested/uncontested split cannot work here. Interference
correction needs to know which competitors are actually PRESENT in this run, not
which ones the library permits.

## Which leaves iteration, and makes the case for it concrete

"Present in this run" is available without circularity if it comes from a
PREVIOUS pass: pass 1's identifications say which competing precursors are real,
and a decoy cannot fake having been identified. Weight each fragment's
contamination by the competitors that were actually identified near that RT and
mobility, then recompute the co-elution features in pass 2 on the down-weighted
set.

That is the missing capability, stated in a form that survives both refutations
above:

  * it is not per-precursor pruning, because the weight comes from other
    precursors' evidence, not from this one's own correlations;
  * it is not library structure, because it counts identified competitors, not
    permitted ones.

It also explains the architectural difference the plan file already names:
DIA-NN refits across ~12 passes against retained chromatograms; we do 2 and
re-extract from the raw file each time. The interference correction is the
reason to iterate, not iteration for its own sake.

Prerequisite, already scoped in the plan file and never built: retain peak
groups between passes (~120 B each against tens of kB per trace) so the second
pass can rescore without re-extracting.

## The MS1 axis was mis-centred by 9 ppm, and the absence conclusion rested on it

kimi (review 46) challenged the interference direction on the grounds that the
7,229 look ABSENT rather than interfered: ms1_coelution -0.093 against 0.473 for
accepted and -0.124 for the 801,458 bulk non-identifications, i.e.
indistinguishable from precursors nobody claims are there. MS1 sits at the
precursor m/z, a far sparser space than the fragment bins, so a present but
fragment-contaminated peptide should still show its monoisotope. The argument is
sound and would retire this whole direction.

Except the measurement was broken. `Ms1Traces::build` matched on the library's
THEORETICAL precursor m/z with no calibration offset, while the fragment axis is
centred on the fitted deviation. Instrumented (6b1b2fa), the run now reports:

    MS1 mass axis centred on -8.949 ppm (borrowed from the fragment fit);
    median residual against that centre -0.378 ppm

So the MS1 error IS the fragment error, and the uncalibrated window was centred
0 ppm away from a signal sitting at -8.949 ppm with a +/-10 ppm half-width --
89% of the way to the edge. Bright precursors survive that. Weak ones do not,
and weak ones are precisely the population the absence conclusion is about.

`var_ms1_coelution` therefore cannot support "absent" in any run before
6b1b2fa, including the AUC 0.757 that justified its 47 GB. Both numbers need
re-measuring on v7.

kimi's other correction stands regardless and is recorded here because it
weakens a claim made in this document: `im_delta` is computed INSIDE the
+/-0.025 extraction acceptance window, so it is caged in [0, 0.025] by
construction, and the decoy median of 0.024 sitting at the window edge shows the
cage IS the dynamic range. "im_delta 0.009 against accepted 0.008" is close to
content-free as presence evidence. Of the two positional pillars in the
partition above, only rt_delta survives scrutiny, and it too is measured inside
a +/-60 s cage.

## Mass-tightening: refuted at 3.6 ppm, and the reason says where the limit is

kimi (review 46, item 3) proposed the cheapest remaining idea: DIA-NN computes
its co-elution sum at base, 0.45x and 0.2x tolerance, on the principle that "a
real peak survives tightening; an interferent often does not". Testable with no
new code, because the same 1,000 precursors are already dumped at several
tolerances, all centred on the same fitted -9.6479 ppm:

    tolerance          AUC (target vs decoy)   median target   median decoy
    3.6 ppm                           0.7537           0.058          0.010
    12 ppm                            0.7658           0.093          0.022
    10 ppm + IM +/-0.050              0.8103           0.121          0.021

    intensity survival ratio I(3.6)/I(12) around the apex: AUC 0.5720
      (against 0.6084 for the wide intensity alone and 0.6138 for the tight)

Tightening makes discrimination WORSE, and the survival ratio is a weaker
feature than either intensity by itself. Both windows were centred identically,
so this is not the mis-centring failure that has explained three other negatives
tonight.

**The reason matters more than the result.** The fragment mass calibration
reports total per-hit scatter of 4.04 ppm, "dominated by irreducible
per-fragment noise". A +/-3.6 ppm window is therefore NARROWER THAN 1 SIGMA of
the noise on a single fragment hit, so tightening to it starves real peaks
faster than interference -- which is what the zero fraction says directly:
19.7% at 12 ppm against 45.3% at 3.6 ppm.

So the principle is not refuted; the WIDTH is. DIA-NN tightens 12 -> 5.4 -> 2.4
ppm and keeps all three as features rather than replacing the base. For us
anything below ~4 ppm is inside the noise. The untested arm is a mild tightening
to ~6 ppm, kept ALONGSIDE the 10-12 ppm sum rather than instead of it, which is
one extraction and costs minutes.
