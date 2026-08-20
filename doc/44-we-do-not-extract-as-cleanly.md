# We do not extract as cleanly as DIA-NN, and it is not the obvious causes

2026-08-20. 1,000 target precursors, S08 diaPASEF. DIA-NN ran with OUR library
and `--xic 60`; ODIA re-extracted around DIA-NN's own apex (`apex_dump_36`), so
the retention-time map is out of the comparison by construction.

## The statistic

Per precursor, build the fragment x RT matrix, take +/-30 s, and compute the
MEAN PAIRWISE PEARSON correlation among fragment traces -- "coherence", which
is what `corr_sum` and `xcorr_shape` exist to reward.

                          mean   median   frac > 0.5
    ODIA                 0.167    0.063        12.2%
    DIA-NN               0.310    0.203        25.2%

That is the IDENTITY-MATCHED number: fragments matched by Product.Mz through
the library's (type, ordinal, charge) annotation, 9.5 common fragments per
precursor, and -- per kimi -- restricted to pairs NON-DEGENERATE IN BOTH tools,
because "fragment absent" and "fragment decorrelated" are the same event and
dropping constant traces silently conditions on a co-elution proxy. The
uncontrolled version gave 0.160 vs 0.317, so the NaN rule was not making the
gap; it is 0.167 vs 0.310 with the rule stated and symmetric.

## What it is NOT

**Not apex placement.** Split by how far our apex sits from DIA-NN's RT, and
additionally re-centre OUR traces on THEIR RT:

    apex agreement      N     ODIA   DIA-NN   ODIA re-centred on DIA-NN RT
        within 5 s    431    0.311    0.484                          0.311
            5-20 s    137    0.064    0.231                          0.063
            > 20 s    419    0.037    0.251                          0.067
               all    987    0.160    0.350                          0.173

The gap survives where the two tools agree within 5 s. Re-centring recovers
almost nothing. (kimi is right that the >20 s row is measured on an overlap
sliver and should not carry weight; the within-5 s row is the argument.)

**Not smoothing.** DIA-NN's dump is not raw counts by our own test -- 2.9% of
its values are near-integer against 100% of ours -- so the confound is real and
had to be excluded rather than assumed away. Smoothing OUR traces to well past
DIA-NN's own autocorrelation does not close it:

    ODIA smoothing   lag-1 ac   coherence
               raw      0.219       0.160
        3-point MA      0.724       0.207
        9-point MA      0.889       0.239
    DIA-NN as given      0.388       0.351

DIA-NN gets MORE coherence at LESS autocorrelation. If anything the raw
comparison understates the gap. kimi separately read `diann.cpp`: the `--xic`
path is `Run::ms2_XIC` -> `Scan::level<false>`, a plain binary search returning
the max-height centroid, with the profile-clipping interference machinery
downstream at quantification. That is verified for 1.7.x; we ran 2.x, which is
closed, so this stays a caveat rather than a settled fact.

**Not fragment count, window or sampling.** 11.7 vs 11.5 fragments per
precursor, 119.1 s span both, 87.0 vs 79.6 points per trace. DIA-NN reports
MORE all-zero traces than we do, 8.7% against 1.6%.

**Not aggregation, and not a too-tight mobility window.** codex's leading
candidate was `Options::aggregate`, which defaults to `Sum` while the comment at
the aggregation site argues for `Max` ("two peaks inside one tolerance are the
same ion split by centroiding far more often than they are two ions"). Measured:

    arm                        coherence   % zeros
    sum, +/-0.025 (baseline)       0.160     45.3%
    max, +/-0.025                  0.131     45.3%
    sum, +/-0.010                  0.089     67.4%
    max, +/-0.010                  0.078     67.4%

Max is WORSE. Narrowing the mobility window is much worse and drives the
zero fraction to 67.4%, so +/-0.025 is not too tight. Note the zero fraction is
identical across aggregation modes: sparsity is set by the MATCHING rule, not
by what is done with the matches.

**Not our m/z tolerance being too wide.** DIA-NN's log records "Optimised mass
accuracy: 12 ppm". The dump above was extracted at 3.6 ppm. DIA-NN extracts
through a window three times wider than ours -- admitting strictly more
interference -- and still doubles our coherence.

## CORRECTION: the headline overstated the gap for production

The dump above was extracted at **3.6 ppm**, chosen for that comparison script.
**Production runs at 10 ppm.** So 0.167 is not ODIA's production coherence, and
the sweep says so. Identity-matched, common-pair-set, against the same DIA-NN
number throughout:

    arm                          n    ODIA   DIA-NN   median diff   ODIA wins
    3.6 ppm, +/-0.025          987   0.167    0.310        -0.088       18.3%
    12 ppm,  +/-0.025          987   0.215    0.310        -0.030       34.4%
    10 ppm,  +/-0.050        1,000   0.253    0.310        -0.030       30.1%

The full sweep, with the zero fraction that explains it:

    arm                        coherence   % zeros
    3.6 ppm, +/-0.025 (dump)       0.160     45.3%
    10 ppm,  +/-0.025              0.206     22.7%
    12 ppm sum, +/-0.025           0.208     19.7%
    12 ppm MAX, +/-0.025           0.144     19.7%
    3.6 ppm, +/-0.050              0.200     30.3%
    10 ppm,  +/-0.050              0.242     12.8%
    10 ppm,  IM off                0.194      5.0%

So the acceptance was too TIGHT, not too loose, and the sparsity was the
symptom: widening m/z alone takes the zero fraction from 45.3% to ~20%, which
is DIA-NN's own 22.9%. Mobility at +/-0.025 is also too tight; +/-0.050 composes
with the wider m/z for the best arm measured. Turning the mobility filter OFF is
worse than +/-0.050, so the filter earns its place -- it is mis-sized, not
wrong. `Max` is worse than `Sum` at every width tested, so the aggregation
comment's argument does not survive measurement.

The honest statement is therefore: at production settings the gap is 0.215
against 0.310, and a one-line window change takes us to 0.253. Not the 0.167
this document led with.

**This is a coherence result, not an identification result.** Wider windows
admit more interference by construction, and doc/38 and doc/43 are both records
of a change that improved a trace statistic and cost identifications. Nothing
here should reach production until a full run measures identifications at
MATCHED empirical FDP.

Still uncontrolled, from the reviews:

- 87.0 vs 79.6 points per trace over an identical span is a 9% discrepancy on
  what should be the same cycle axis. One tool is seeing a different scan set.
  kimi calls this out and it is not explained.
- ODIA applies BOTH a frame-band filter and a per-precursor mobility filter
  (`ChromatogramExtractor.cpp:1111-1169`); only the latter was swept.
- DIA-NN 2.x may filter the dump in ways 1.7.x did not.
- The pointwise identity test -- can we reproduce DIA-NN's non-zero points
  exactly with matched ppm, max aggregation and no IM gate -- has not been run,
  and it is the test that would settle "raw vs processed" without source.

## Why this matters beyond the number

Median `var_library_corr` on the full run at TWELVE fragments is 0.0137. That
is not a scoring defect to be tuned away; it is the extraction telling the truth
about traces that do not co-elute. Every co-elution feature we have -- corr_sum
(AUC 0.906, and the picker's own statistic), xcorr_shape, the Gate C evidence
statistic -- is computed on these traces. A scorer cannot rank on structure that
extraction did not deliver.

## Pointwise: the dumps are on the same axis, and our windows are MIS-CENTRED

kimi's mandatory test -- can we reproduce DIA-NN's points ourselves -- run by
matching (precursor, fragment m/z, cycle RT) across 835,374 cells:

    our arm                        median ratio ODIA/DIA-NN   both non-zero
    12 ppm, MAX, +/-0.025                              0.53         588,164
    12 ppm, SUM, +/-0.025                              1.34         588,164
    10 ppm, SUM, IM off                                3.66         645,538

Two things fall out.

**The 9% points-per-trace discrepancy is not real.** Distinct RT values per
precursor are 87.0 for BOTH tools. The tools are on the same cycle axis; the
apparent difference came from how all-zero traces were counted. kimi flagged it
as one tool seeing a different scan set; it is not.

**DIA-NN's value sits BETWEEN our max and our sum**, nearest to summing inside a
mobility acceptance NARROWER than our +/-0.025. But our +/-0.010 arm measured
WORSE, with the zero fraction jumping 45.3% -> 67.4%. A narrower window cannot
be simultaneously better for DIA-NN and worse for us unless ours is pointed
somewhere slightly wrong: mis-centred, tightening clips real signal on one side
while interference still enters on the other.

That is consistent with what the run itself reports. The mobility seed applies
ONE GLOBAL offset -- 0.0170 1/K0, with a robust sigma of 0.0241, i.e. a spread
wider than the offset it is correcting -- and per-precursor ion-mobility
calibration is DEFERRED to pass 2 ("it is measured at the peak groups this run
scores, and none have been scored yet"). So every pass-1 candidate, which is
what the classifier trains on, is extracted through a window centred on the
library's predicted 1/K0 plus a constant.

Predicted, and not yet tested: pass-2 traces should be measurably more coherent
than pass-1 traces for the same precursors, because pass 2 is the first to
centre on observed mobility. If they are not, the mis-centring hypothesis is
wrong and the residual gap is elsewhere.

## The mis-centring is measured, and the run measures it too late

The prediction above did not need a new run. The baseline's OWN pass-2 report
quantifies exactly the error the pass-1 window carries:

    ion-mobility calibration: GATE PASSED -- 3 of 4 charges corrected, mz_shaped;
    71% of the mean squared 1/K0 error removed out of fold

    charge 2: 4045 anchors, +0.0150 at 330 Th -> +0.0146 at 1293 Th;
              1/K0-linear -0.0393 per 1/K0, i.e. a -3.9% error in the CCS->1/K0 coefficient
    charge 3:  980 anchors, +0.0092 at 334 Th -> +0.0124 at 1066 Th; +1.5% coefficient error
    charge 4:  129 anchors, constant +0.0192;                        +2.1% coefficient error

The library's predicted 1/K0 is wrong by **+0.0092 to +0.0192**, against a
window half-width of **+/-0.025**. That is an offset of up to 76% of the
half-width: for a large share of precursors the true mobility sits near or
outside the window edge, so the trace samples the shoulder of the mobility peak
and drops out whenever the ion drifts across the boundary. That is the 45.3%
zero fraction, and it is why NARROWING to +/-0.010 made things worse -- it
clipped what little was inside.

The correction is charge-dependent AND m/z-shaped. The mobility seed applies ONE
global constant (+0.0170), which roughly fits charges 2 and 4 while
over-correcting charge 3 by ~0.006. It cannot represent a -3.9% coefficient
error in the CCS->1/K0 conversion, which is a SLOPE error, not an offset.

And the ordering is the problem: the pass-1 attempt reports

    ion-mobility calibration: GATE FAILED -- the gate passed, but no charge has
    enough anchors to fit from (min_anchors_per_charge = 120)
      residuals: 213 target ... from 343 sampled precursors

so the machinery exists and is starved, while pass 2 gets 5,156 residuals and
succeeds. Pass 1 forms all 2,835,055 peak groups and trains the classifier
through the uncorrected window; pass 2 earns the correction afterwards.

The gate's own comment defends the safe direction -- "an uncentred window keeps
the library's own error, where a window recentred on a badly measured offset
moves off the precursor entirely" -- and that is right as a refusal rule. It is
not an argument for having only 213 anchors to decide on.

Candidate fixes, in order of cost:

1. Fit the seed offset PER CHARGE rather than globally. The CiRT seed already
   searches 298 standards and the mobility seed already reports its own robust
   sigma (0.0241, wider than the 0.0170 offset it applies) -- that sigma is the
   charge and m/z spread being averaged away.
2. Widen pass 1 only, to +/-0.05, so the uncorrected window still contains the
   peak; the measured worst-case offset is 0.0192, which +/-0.05 covers and
   +/-0.025 does not. The 10 ppm / +/-0.050 arm is the best coherence measured
   (0.253 against DIA-NN's 0.310), which is consistent with this being the
   mechanism rather than a coincidence.
3. Lower `min_anchors_per_charge` for the pass-1 fit, or pool charges for the
   offset while keeping the slope global. Cheapest to try, weakest justification.

Option 2 is the one to test first because it is one flag and it does not change
what is estimated -- only how much of the peak is inside the window while the
estimate is still unavailable.
