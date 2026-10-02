# Phase 2: restore the ion-mobility dimension

Status: plan. Measured basis below.

## What is wrong, measured

Validation against DIA-NN's confident set on IH1 (2,665 precursors, stratified
over RT x m/z x charge, not top-by-confidence):

* **20.7% recovered** (552/2665) against a **5.5% decoy null**.
* Recovery by DIA-NN abundance decile: bottom four **5.2-6.8%** — statistically
  indistinguishable from decoys — against **89.5%** in the top decile.
* The RT axis is correct: median apex offset **0.00 s**, SD 10.5 s. This is not
  a calibration failure.
* For non-recovered precursors, only 2.6% have their maximum anywhere near the
  true RT. The peak is absent, not displaced.

Root cause, probed directly rather than inferred:

* **0 of 32,210 spectrum entries carry a finite ion-mobility band.** All were
  widened to (-inf, +inf) by the degenerate-band fallback in `MzPeakSource.cpp`
  — a fallback added earlier the same day to defend against an *inverted* band,
  which is a different file's problem.
* Extraction with and without `-no_ion_mobility` is **byte-identical**
  (md5 `7c9c6562ded89b25378eb3a4c53b4838`), so the filter is provably inert.
* IH1 packs **two isolation windows per frame, separated only in mobility**.
  `peaks()` hands both window entries the same merged peak list, so with no
  band each window's transitions are matched against the other window's ions
  across the whole 0.60-1.40 1/K0 range. **The mobility dimension is collapsed.**

That is what produces the 84% non-zero point fraction and the flat
~1500-2000-count background that buries everything below the top ~20%.

## Why the current test cannot work, even when the band is finite

`ChromatogramExtractor.cpp:443-444` tests each peak against **the frame's own
band**. Every peak in a frame came from that frame, so the test is a tautology
that can only ever reject nothing or — when the band is inverted — everything.
Both failure modes have now been observed in this project.

## The fix

Filter each peak against **the library's expected 1/K0 for the precursor whose
transition it is being matched to**, not against the frame.

The data exists on both sides and neither is used:

* peaks carry per-peak 1/K0 (probe: `hasIM=1`, range 0.6000-1.4007);
* `Library::precursors().im` carries per-precursor 1/K0 and **agrees with the
  observed value to SD 0.019** — read, stored, copied, written back out, and
  consumed by nothing during extraction.

Concretely, in `MzIndex`: alongside each transition's m/z, store the expected
1/K0 of its precursor. In the match inner loop, reject a peak whose 1/K0 falls
outside `expected +/- tolerance`. That keeps the existing log-m/z bucketing and
adds one comparison per candidate.

Tolerance: the measured library-vs-observed SD is 0.019, so a default of
**+/-0.05** is ~2.6 sigma. It must be an option, not a constant, and its
effect must be measured rather than assumed.

## What must not break

* **A run with no per-peak mobility must be unaffected.** `12_80` is SCIEX
  SWATH with no mobility at all; the filter has to disable itself and say so,
  not silently reject everything. This is the failure mode that already
  happened once.
* **A precursor with no library 1/K0 (NaN) must not be filtered out.** Absent
  information is not evidence of mismatch.
* **Targets and decoys must be treated identically.** A decoy inherits its
  target's precursor m/z and 1/K0, so this is symmetric by construction — but
  assert it rather than assume it.

## How to know whether it worked

The same validation, rerun unchanged:

1. **Recovery fraction** against DIA-NN's confident set, currently 20.7%.
2. **Recovery by abundance decile** — the bottom four deciles at 5.2-6.8% are
   the number that has to move. If overall recovery rises while the bottom
   deciles stay at the decoy null, the filter has only sharpened peaks that
   were already found.
3. **Decoy null** must stay near 5.5% or fall. A filter that raises both is
   selecting on something other than peptide evidence.
4. **Non-zero point fraction**, currently 84.1%, and the median background,
   currently ~1500-2000 counts. Both should fall substantially; if they do not,
   the band is not actually being applied.
5. **`12_80` output must be byte-identical to today's**, proving the no-mobility
   path is untouched.

A supplementary probe quantifying the expected gain from a +/-0.05 band was
still running when this was written; its number belongs here when it lands.
