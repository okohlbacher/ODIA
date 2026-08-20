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

## What is still open

The first factorial never tested the cell that matches DIA-NN: WIDE window with
MAX aggregation. At 3.6 ppm there is little to choose between peaks inside the
tolerance, so `max` was tested where it cannot help; at 12 ppm the choice is the
whole point. `-fragment_ppm 12 -aggregate max` is running.

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
