# The mobility boundary is a midpoint, and the windows are not symmetric

IH1 co-packs two isolation windows into every MS2 frame, separated only in ion
mobility. The mzpeak conversion carries each window's mobility POSITION
(`ion_mobility_value`) but not its bounds -- `ion_mobility_lower_limit` and
`_upper_limit` do not exist in the schema at all, over all 32,210 entries. So
`MzPeakSource` derives the split as the **midpoint between the two positions**
(`MobilityBands.cpp`), which is what its comment says it does.

That rule is exact only when the two co-packed windows are equally wide in
mobility. On IH1 they are not, and the error is measurable.

## The measurement

The vendor states the geometry exactly, in `DiaFrameMsMsWindows` inside the .d:
12 window groups x 2 windows, each a (m/z x mobility-scan) box. The scan axis
converts to 1/K0 essentially linearly -- fitting the 24 (scan centre, mobility
position) pairs gives `im = -0.000881*scan + 1.43065` with a **maximum residual
of 0.0004** over all 24 windows, so the conversion is not the uncertainty here.

Evaluating that fit at each group's stated boundary scan, against the midpoint
ODIA derives:

    grp  split scan  scan widths   vendor 1/K0   ODIA midpoint     error
      1     602       568 / 342       0.9003        0.9498       +0.0494
      2     579       545 / 365       0.9206        0.9604       +0.0398
      3     568       534 / 376       0.9303        0.9648       +0.0345
      4     557       523 / 387       0.9400        0.9700       +0.0301
      5     545       511 / 399       0.9506        0.9753       +0.0248
      6     523       489 / 421       0.9699        0.9850       +0.0151
      7     511       477 / 433       0.9805        0.9903       +0.0098
      8     488       454 / 456       1.0008        1.0000       -0.0008
      9     477       443 / 467       1.0105        1.0053       -0.0052
     10     443       409 / 501       1.0404        1.0203       -0.0202
     11     409       375 / 535       1.0704        1.0352       -0.0351
     12     295       261 / 649       1.1708        1.0855       -0.0853

**The error tracks the asymmetry exactly.** Group 8 splits at scan 488 of the
34-944 range -- 454 against 456 scans, as near equal as the scheme gets -- and
the midpoint rule is off by **-0.0008**, i.e. right. Groups 1 and 12 are the
most lopsided (568/342 and 261/649) and are off by +0.0494 and -0.0853. The
sign flips where the asymmetry does. That is the mechanism confirming itself,
not a correlation.

## Why it matters

The full spread is 0.135 in 1/K0. DIA-NN, on this same run, auto-selected an IM
window of **0.059045**. So at the extremes the boundary is misplaced by more
than a whole mobility window -- 1.4x it at group 12.

The band decides which of the two co-packed windows a peak in the merged list
belongs to. A boundary that sits 0.085 too low hands a slab of one window's
mobility range to the other: signal extracted from the wrong band, and
interference that looks like the instrument's rather than like a bug. This is
precisely the failure mode the `MzPeakSource` comment anticipated for a refused
derivation -- except the derivation did not refuse, it succeeded and was wrong.

## Two things this cost before it was found

**Nothing wrote the windows down.** ODIA derived them, DIA-NN's report does not
carry them, and neither engine emitted the scheme, so no comparison was
possible in either direction. `-out_windows` exists now for that reason.

**The fixture cannot show it.** Measured first on `ih1_6x60.mzpeak`, every one
of the 24 windows came back with `im_low = -inf, im_high = +inf`, which reads
as "ODIA collapsed the 2-D geometry" -- a much more dramatic and completely
wrong conclusion. The fixture has exactly **one** selected ion and one
precursor per spectrum, all 6,120 of them, against the full file's two. The
fixture flattened the co-packing during construction, so `deriveMobilityBands`
correctly returned `NotNeeded` (fewer than two windows), which by design does
not warn. On the full file the bands are derived, complementary, and share
their boundary exactly: windows 0-11 get `(-inf, X]`, windows 12-23 get
`[X, +inf)`, and window k's `im_high` equals window k+12's `im_low` to the last
digit.

Add it to the list of things the fixture cannot answer, beside FDR, standing
against DIA-NN, and retention-time calibration: **it cannot answer anything
about the acquisition geometry, because it does not preserve it.**

## The fix: the true bands are already inside the .mzpeak

The claim above -- that the midpoint is "the best available guess from what the
mzpeak carries" -- is **wrong, and was corrected the same day**. The conversion
does not lose the bands. It keeps them, in the embedded vendor method:

    vendor/<method>.m/diaSettings.diasqlite  ->  DiaWindowsSpecification
    (Id, Type, CycleId, OneOverK0Start, OneOverK0End, IsolationMz,
     IsolationWidth, CollisionEnergy)          25 rows

    Id  Cycle  OneOverK0Start  OneOverK0End  IsolationMz  IsolationWidth
     2      1            0.90          1.40      718.845           23.05
     3      1            0.60          0.90      398.160          141.28
     4      2            0.92          1.40      742.620           26.50

Group 1's boundary is stated as **0.90 exactly**, against the 0.9498 ODIA
derives -- the +0.0498 this document measured, now confirmed against the
vendor's own number rather than against a fitted scan axis.

Two further things that table settles:

* **The true bands are bounded on BOTH sides** -- 0.60-0.90 and 0.90-1.40.
  ODIA's derived bands are one-sided, `(-inf, X]` and `[X, +inf)`. The outer
  edges are real and are currently unbounded.
* **No .d dependency.** The file is inside the .mzpeak, so reading it needs no
  vendor SDK, no re-conversion, and no format change -- only that the reader
  open one more entry of the archive it already has open.

So the fix is to read `OneOverK0Start`/`OneOverK0End` verbatim and stop
deriving anything. `deriveMobilityBands` stays as the fallback for files that
carry no vendor method.

Until that lands, the derived band is systematically wrong wherever the
co-packed windows are unequal, which on this instrument is ten groups of twelve.

## What this is NOT

An adversarial review claimed from these numbers that "ODIA's mobility aperture
excludes half the reference population before any scoring". Checked against the
run log, that is wrong, and the correction matters more than the claim:

    full_prom10.log:41   ion-mobility calibration: DEFERRED -- this pass
                         extracts on the library's 1/K0; the next one is where
                         the correction lands
    full_prom10.log:110  GATE PASSED -- 2 of 3 charges corrected, 70% of the
                         mean squared 1/K0 error removed OUT OF FOLD

Pass 1 does extract on the raw library axis, and that axis is genuinely bad:
`observed - library iIM` over DIA-NN's 39,115 confident precursors has a median
of **+0.0232**, from a scale error of **-4.66%** (`obs = 0.9534*lib + 0.0693`),
so a +/-0.025 window centred on the library value admits only **54.09%** of
them. ODIA measures the same defect itself and calls it "-3.7% error in the
CCS->1/K0 coefficient" for charge 2.

But pass 1 exists only to collect anchors, and it collected 5,622 at q <= 0.01.
**Pass 2 re-extracts everything on the calibrated axis** -- both passes walk all
32,210 spectra, both report the same 1,119,490 precursors covered by no
isolation window, and pass 2 ends with MORE candidates than pass 1 (15,333,736
against 13,338,052). The pass-1 mobility loss does not propagate to scoring.

What remains true is narrower: after ODIA's own correction the residual is p95
**0.0362**, still beyond the +/-0.025 half-window, so the aperture is tight
against what the calibration can deliver. That is a sizing question, not the
funnel.

## What the fix changed, measured on the library rather than end to end

The strategy this sits under judges a phase on its own job, so the question is
not "how many identifications" but **how many precursors change what the
instrument is asked for**. Over the 4,991,901 target precursors, comparing the
window set that admits each one under the derived bands against the stated ones:

    same window count                    4,659,405   93.34%
    admitted by MORE windows now            74,745    1.50%
    admitted by FEWER windows now          257,751    5.16%
    covered by NO window -- derived         768,923   15.40%
    covered by NO window -- stated          953,234   19.10%

The stated bands are **stricter**, and that is the correction, not a loss. Of
the 257,282 targets newly excluded:

    outside the acquisition range [0.6, 1.4]   241,636   93.92%
    inside it, outside their window's box       15,646    6.08%

The instrument acquired 1/K0 in [0.600, 1.400] (`GlobalMetadata`
`OneOverK0AcqRangeLower/Upper`). **13.45% of the library -- 671,394 targets,
666,908 of them above the ceiling -- predicts a mobility the instrument never
sampled.** The derived bands were one-sided, so every one of those was admitted,
placed, extracted, and given a full-length trace that could not contain its
signal. The stated bands reject them.

The remaining 15,646 are the real diaPASEF geometry rather than an edge case: a
window is a BOX, so m/z 400 at 1/K0 1.0 was never acquired even though both
coordinates are individually in range -- the scheme samples a diagonal, not a
rectangle. Splitting those by whether the mobility calibration recovers them:

    covered again once the IM calibration is applied   7,882
    genuinely outside the scheme's coverage            7,764

So the tight bands cost at most **7,764 targets (0.16%)** that the scheme does
not sample, plus 7,882 that only the *uncalibrated* library axis pushed out --
and pass 2 extracts on the calibrated axis, where they return. Against that,
~250,000 precursors stop being extracted from windows that never contained them.

**The library, not the extractor, is where the 13.45% belongs.** A CCS model
predicting 1/K0 up to 2.2962 for an instrument that stops at 1.400 is a
library-construction problem, and it is now visible because the geometry is.

## Correction: "precursors excluded" is an offline proxy, not what the tool does

The impact table above counts precursors whose (m/z, 1/K0) falls in no window.
An adversarial review checked the extractor and the count does not describe
pipeline behaviour:

* **Precursor-to-window assignment is m/z only** --
  `ChromatogramExtractor.cpp:492-501` tests `windows[w].contains(mz)`, and
  `contains()` is `mz >= mz_low && mz <= mz_high`. Mobility is not consulted.
* **The band filters PEAKS**, at `ChromatogramExtractor.cpp:1153`:
  `if (use_band && (peak_im < im_low || peak_im >= im_high)) continue;`

So no precursor is dropped by this change. Every one is still assigned and
still extracted; what changes is which peaks are allowed to enter its trace.
The right statement is that **the fix converts cross-talk into zeros**: under
the derived boundary a precursor assigned to window A, whose true mobility slice
was B's, integrated B's fragments and scored them as its own. Removing that is
the point, and because cross-talk that scores well is also how a false anchor
enters pass 1, anchor collection should improve rather than degrade.

The same correction disposes of a worry raised against the change -- that tight
bands in pass 1 would censor the precursors needed to fit the mobility
calibration. **The band test never sees the library 1/K0.** It compares a RAW
peak's mobility, in the instrument frame, against a band in the same frame. The
-4.66% library scale error cannot push a peak across a boundary; it enters only
the per-transition tolerance around `lib_im` (`:1164`), which this change does
not touch. The "7,882 return once calibration is applied" figure is a property
of the offline proxy, not of the extractor.

What the review did find, and what has been fixed:

* **The outer edges were a real regression.** Peaks run to 1.4007 (the TIMS scan
  range pads past the stated ceiling; see `SpectrumSource.h:82`) while the
  method states 1.400, and the band test is half-open -- so every peak in
  [1.400, 1.4007] was dropped from every outermost window, where the derived
  band had been unbounded. The internal splits are what the method is worth;
  its outer edges are where hardware padding lives. Those ends are now left
  open.
* Entry selection tested a SUBSTRING, so `diaSettings.diasqlite-wal` would have
  matched; and the backup rule was a substring on the whole path, so a live
  method named `backup-final.m` would have been skipped in silence. Now a
  suffix test and a path-component test.
* Two live methods in one archive, or one isolation centre stated with two
  different bands across cycles, now REFUSE and fall back to derivation with a
  warning, rather than taking whichever row the archive listed first.
* A band the file states per spectrum now outranks the method, rather than
  being overwritten by it.
* The empty path is no longer silent; a run with mobility that finds no usable
  method says so.
* `mkstemp` rather than a pid-derived name in a shared /tmp, and a size cap
  before allocating.

And one check that was worth running: each window's own mobility position must
lie inside the band the method states for it. Written naively this fired on
16,105 of 32,210 windows -- exactly half -- which is not a stale method but the
writer attaching every selected ion of a frame to precursor 0, a quirk
`MobilityBands.cpp` already documents. Matched by m/z instead, **all 32,210
pass**. That is meaningful: the m/z centres agreeing only proves the method's
m/z column is right, while the positions agreeing is evidence the 1/K0 columns
belong to this acquisition too.
