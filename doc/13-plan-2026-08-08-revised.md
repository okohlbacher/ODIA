# Revised plan, 2026-08-08 04:00 — after a night of measurement

`doc/12` was written at 23:00 and its critical path (doc/08's prefilter) was
dead by 02:00. This replaces it. Everything below is measured, and the
measurements are in `doc/BACKLOG.md` under the 2026-08-07/08 headings.

## 1. The one structural finding

**Presence statistics saturate on this data; shape statistics discriminate.**

Seven quantities were measured on the SEARCH benchmark (IH1 + `v6_50k`, 670 true
of 50,000, base 13.4 per 1000):

    fragment depth, whole frame                1.0x   saturated (49,999/50,000)
    fragment depth, mobility-sliced            1.0x   saturated (49,835/50,000)
    qualifying_spectra                         0.5x   worse than random
    total_matches                              0.5x   worse than random
    MS1 isotope depth                          1.0x   saturated (49,916/50,000)
    MS1 monoisotopic intensity                 1.7x   weak, abundance-confounded
    MS1/MS2 co-elution correlation            13.7x   top bin; 1.4-1.5x in bulk

The first six ask "does this exist somewhere in the run". Each is a maximum or
a count over 10^3-10^4 spectra, so each is an EXTREME-VALUE statistic over
thousands of draws: with 15 ppm tolerance and mobility-merged frames (a diaPASEF
frame is ~600-810 TIMS scans concatenated), coincidence is near-certain and the
answer is yes for everything. Making each draw more selective does not fix it —
the draw count is the problem.

The seventh asks whether two things rise and fall TOGETHER. Coincidences would
have to arrive in the right ORDER, which they do not.

This also explains, retrospectively, the largest win this project has had:
replacing the amplitude peak-picker with co-elution detection fixed RT, FDR,
recovery and memory simultaneously. **Treat this as a design rule: prefer shape
over presence, and distrust any feature that is a max or a count over the run.**

## 2. What is dead, and why

* **doc/08's prefilter.** Its 1,800x separation was measured on the reference
  engine's non-mobility data and does not transfer. Every rescue is blocked or
  circular: a tighter tolerance needs the mass calibration that FAILS its gate
  on IH1; an RT neighbourhood is circular because supplying the RT seed was the
  filter's second purpose; a per-TIMS-scan unit is a reader change blocked on
  the same decode work as `ODIAInfo -peaks`.
* **Feature COUNT as a lever.** The 15-to-30 "cliff" was an artefact of drawing
  features independently. With a shared latent — ODIA's actual regime, since all
  15 sub-scores read MS2 fragment traces — the count buys nothing: 0
  identifications at 4, 8, 15, 30, 60 and 110 alike. **Adding more MS2-derived
  sub-scores is not worth doing.** Orthogonality is the lever, not count.
* **The mass calibration thread.** Settled (the window WIDTH costs 6.5x what the
  centring does) and closed. Both flags off by default.
* **`ms1_max` and `ms1_iso` as sub-scores.** Saturated, and an abundance proxy.

## 3. The critical path

### A. `MS1_COELUTION` as a sub-score

The only measured, orthogonal, non-saturating evidence available. Survives
stratification by intensity (1.4x within-stratum against 1.5x unstratified), so
it is shape and not brightness.

1. **Extract MS1 traces in `ChromatogramExtractor`**, on the SAME RT and 1/K0
   windows as the MS2 extraction, so the two are comparable cycle for cycle.
   `MzPeakSource::ms1Spectra()`/`ms1Peaks()` already exist (2026-08-08). The
   probe used a nearest-MS1-bin approximation; the extractor must not.
2. **Then** add the sub-score. NOT before — `doc/07` is explicit that a
   sub-score computed from a placeholder is worse than an absent one, and
   `IM_DELTA` sat all-NaN for a week because its plumbing point was added and
   never connected. Do not repeat it.
3. Correlate the MS1 precursor trace against the MS2 fragment consensus over
   bins where either is non-zero. Correlating over the whole gradient is
   dominated by jointly-empty bins and reproduces saturation.
4. Measure on `v6_50k`. **Expect a real but bounded gain** — the bulk effect is
   1.4-1.5x, and the 13.7x top bin is 14 precursors with a +/-3.6x error bar.
   One good feature does not close a 6,500:1 likelihood-ratio requirement.

### B. Entrapment, because the benchmark itself is biased

The target is DIA-NN's confident set, so **any feature that predicts "what
DIA-NN finds" scores well whether or not it predicts correctness.** Abundance
predicts exactly that. This is a hazard in every measurement in this document.

`test/tools/odia_entrapment.cpp` now plants the production 1.5% prior
(2026-08-08) and reproduces the collapse in seconds. It is the only check that
does not share the comparator's biases. Use it as the development loop, and get
a real entrapment search onto IH1 before any recovery number is published.

### C. The picker, which is where the losses actually are

Astral: 9,890 precursors yielded no candidate peak group over 12.9M
scan-position rejections; the chromatograms were extracted. IH1: 19,150. Three
criteria dominate on both — `not a local maximum`, `min_corr_score`,
`apex_evidence`. `max_corr_diff` is inert (26k of 189M) and is not worth tuning.

Not yet separated, and needing a counter rather than a grep: candidates emitted
and then dropped by `min_fragments_at_apex`, versus never emitted.

## 4. Still-unbuilt preconditions from the ORIGINAL plan

* `doc/07` step 2 — the **OpenSWATH cross-check** of the sub-scores. Never done.
  Two sub-scores were later found worthless (`USABLE_FRAGMENTS` constant 12.000;
  `IM_DELTA` all-NaN); this step exists to catch exactly that.
* **D10's mzML reference source.** Still no second `SpectrumSource`, so no
  reader bug can be told from an extractor bug.
* **D1's determinism harness** — partly discharged 2026-08-08: the
  extract-score pipeline is bit-identical across a 4x thread change. Library
  generation, ONNX, and permutation invariance remain unmeasured.
* **Expose the learned model** from `scoreSemiSupervisedLDA` so batching can be
  tested honestly (train on a batch, score the full set). The current refutation
  tested a straw man and is not established.

## 5. Order

C (counter, one change) → B (entrapment loop) → A1 → A2 → A4.

A before more scoring work of any kind, and nothing from the `doc/13` MS2 score
inventory until orthogonality is demonstrated to matter more than count —
which, measured, it is.
