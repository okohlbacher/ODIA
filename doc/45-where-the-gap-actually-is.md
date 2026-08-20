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
