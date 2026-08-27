# The mobility boundary is a midpoint, and the windows are not symmetric

S08 co-packs two isolation windows into every MS2 frame, separated only in ion
mobility. The mzpeak conversion carries each window's mobility POSITION
(`ion_mobility_value`) but not its bounds -- `ion_mobility_lower_limit` and
`_upper_limit` do not exist in the schema at all, over all 32,210 entries. So
`MzPeakSource` derives the split as the **midpoint between the two positions**
(`MobilityBands.cpp`), which is what its comment says it does.

That rule is exact only when the two co-packed windows are equally wide in
mobility. On S08 they are not, and the error is measurable.

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

**The fixture cannot show it.** Measured first on `s08_6x60.mzpeak`, every one
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

## The fix, not yet made

The midpoint is the best available guess from what the mzpeak carries. Three
ways out, in increasing order of correctness:

1. Derive the boundary from the observed mobility distribution of the merged
   peak list -- the two windows' ion populations are separated by a real gap,
   and the gap's location does not assume symmetry.
2. Have the converter preserve `ion_mobility_lower_limit` / `_upper_limit`.
   They exist in the vendor table; nothing but the conversion loses them.
3. Read `DiaFrameMsMsWindows` directly when the run is a .d.

Until one of them lands, the derived band is systematically wrong wherever the
co-packed windows are unequal, which on this instrument is ten groups out of
twelve.
