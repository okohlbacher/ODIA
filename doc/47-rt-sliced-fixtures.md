# RT-sliced fixtures: 7.5x faster, and not a substitute for the full run

2026-08-21. Question: can a 10-20% retention-time slice of S08 and Astral give a
test cycle at least 5x faster while still producing meaningful numbers?

Short answer: **7.5x is achievable and both fixtures work in ODIA and DIA-NN,
but the numbers are NOT meaningful for anything FDR-shaped.** Details below,
including three separate attenuations that only appeared on measurement.

## The fixtures

    fixture        slices   window   RT span        MS2 spectra   of full
    s08_6x60       6 x 60s  487-1512s              6,120         19.0%
    astral_7x60    7 x 60s  224-2126s             54,533         17.9%

Built with OpenMS `FileFilter -rt A:B -rt_block_mode shrink_to_preserve_full_cycle`,
one pass per slice in parallel, then `FileMerger`, then `mzpeak-convert` for
ODIA. DIA-NN reads the merged mzML directly.

Verified on both: **all isolation windows survive** (24 for S08, 150 for
Astral), ion mobility survives for diaPASEF, and absolute retention time is
preserved so the full run's RT map still applies.

### Three things that had to be got right

**Slice where peptides actually elute.** S08's accepted identifications live
between ~400 and ~1600 s of an 1859 s span; the first layout put 2 of 6 slices
on empty gradient. Placed correctly, 19.0% of the FILE is 31% of the PRODUCTIVE
gradient, which roughly doubles the identifications per second retained.

**Slice width against real peak width.** Median FWHM on S08 is 5.5 s and p90 is
9.7 s -- not the 90 s the peak-group boundaries suggest, which is the extraction
window saturating. 60 s slices are ~6x p90, so edge truncation costs the outer
~10 s of each slice.

**`mzpeak-convert --rt` cannot be used.** It REBASES retention time: slices at
8.25-9.25 min and 20-21 min both come back reporting 1.0-60.8 s. Fine for one
slice, useless for a merge. FileFilter preserves absolute time.

**S08's mzML had to be repaired first.** It carries `<dataProcessingList
count="0">` and no `defaultDataProcessingRef` on EITHER `spectrumList` or
`chromatogramList`, so OpenMS refuses it. Patched by a single streaming pass
that fills the list, adds both refs, and drops the `indexedmzML` wrapper (the
trailing byte index cannot survive a header insertion).

## Speed: the target is met, and better than linear

    ODIA on S08          full (v5)      fixture      ratio
    wall                 3h00:00        24:08        7.5x
    peak memory          124 GB         40 GB        3.1x
    MS1 traces build     3,805 s        141 s        27x
    MS1 bins             1,343          255          5.3x

Memory and bins scale with the spectra fraction exactly, but the MS1 build is
27x faster rather than 5x -- the full run's build is memory-bandwidth-bound at
46 GB and 9 GB is not. That superlinearity is where 5.3x becomes 7.5x.

DIA-NN on the same fixtures: S08 9:58, Astral 0:07.

## But the numbers are not meaningful. Three attenuations

**1. Retention-time calibration does not survive slicing.**

    run       slope      intercept   p95 target   p95 decoy control
    full      1086.50    473.77      25.1 s       not fittable
    fixture      2.35   1084.05      18.2 s       15.1 s (fits BETTER)

The slope collapses to a horizontal line and the seed correctly refuses itself.
The mechanism is inherent, not a tuning failure: a CiRT standard whose true
elution falls in a GAP does not disappear, it acquires a spurious apex inside
some retained slice. At 31% coverage ~69% of anchors are garbage scattered at
random RT, and no RANSAC reaches consensus on that. Fixing it needs >50%
coverage, which is not a 10-20% fixture. **A sliced fixture must be given its RT
map explicitly** (`-irt_slope 1086.50 -irt_intercept 473.77`), and therefore
cannot test RT calibration at all.

**2. It cannot resolve an effect that halves the FDP.**

    arm                     IDs      entrap    FDP           effect
    FULL  floor off      14,081         300    12.63 +-0.73
    FULL  floor on       13,268         130     5.72 +-0.50   +6.91 pp
    FX    floor off       2,986          44     8.64 +-1.30
    FX    floor on        3,282          42     7.49 +-1.16   -1.15 pp

The fragment floor is worth **6.91 pp** on the full run. The fixture measures
**1.15 pp**, against a combined Poisson sigma of 1.74 pp -- 83% attenuated and
not distinguishable from zero.

And this is NOT a sample-size problem that more slices would fix. The floor
works by removing 366,084 weak 0-2-fragment targets that have no decoy
counterpart; on a fixture those precursors largely fail to produce candidates at
all, so removing them changes little. **The failure mode the floor repairs is
mostly absent from the fixture.**

**3. It compresses the gap to DIA-NN, in our favour.**

                        DIA-NN        ODIA          ratio
    full run            39,211        13,268        2.96x
    fixture              4,606         3,282        1.40x
    FDP: DIA-NN 7.42% -> 4.44%,  ODIA 5.72% -> 7.49%

The two tools' false-discovery proportions move in OPPOSITE directions on the
fixture. Any competitive number taken from it flatters us by roughly 2x.

## What the fixture is for

USE IT FOR -- 7.5x, 40 GB, 24 minutes:
  * crash, regression and does-it-still-run testing;
  * performance and memory work (the phase ratios are faithful, and the MS1
    superlinearity means it is if anything conservative about speedups);
  * extraction-level measurement that does not depend on the FDR -- trace
    coherence, mass and mobility calibration residuals, zero fractions,
    per-fragment statistics;
  * large-effect identification changes, where a 2x move is not in doubt.

DO NOT USE IT FOR:
  * any FDP or FDR acceptance decision -- effects attenuate ~5x and 1-2 pp
    differences are inside the noise;
  * DIA-NN competitive comparison -- the ratio compresses 2.96x -> 1.40x;
  * retention-time calibration work -- it cannot seed a map.

## Open

  * Whether a LARGER fixture (say 40%) removes the attenuations, and at what
    speedup. On the mechanism above the floor effect should partly return, since
    weak precursors would start producing candidates again. That is one build and
    two runs.
  * ODIA has not yet been run on `astral_7x60`; Astral's truth set is stale
    (doc/BACKLOG) and must be re-derived before its numbers mean anything.
  * The S08 fixture is 21 GB of mzML for 19% of the data because FileFilter
    writes uncompressed binary arrays. `mzpeak-convert` gets it to 4.4 GB. If
    the mzML is only ever fed to DIA-NN, converting with compression would make
    the fixture far more portable.
