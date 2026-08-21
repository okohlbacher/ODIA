# Backlog

Open items, newest first within each group. Items that need a decision or an
external fix are marked **[you]**; the rest are mine to work through.

---

## THE DEFAULT BENCHMARK: `scripts/bench.sh` (2026-08-21)

One command, ~25 minutes, for iterating on performance and trying ideas:

    scripts/bench.sh baseline
    scripts/bench.sh no-floor  -min_library_fragments 0
    scripts/bench.sh wide-im   -precursor_im_window 0.05

Runs the RT-sliced S08 fixture (`s08_6x60`: 6 slices x 60 s over 487-1512 s,
6,120 MS2 spectra = 19.0% of the file, all 24 isolation windows, ion mobility
intact), prints the arm against stored references, and appends to
`shared/libv2/bench_results.tsv`. Construction and evidence in doc/47.

    ODIA on S08        full (v5)   fixture    ratio
    wall               3h00        24:08      7.5x
    peak memory        124 GB      40 GB      3.1x

Reference arms, all on the same fixture:

                       IDs   entrap    FDP        est.true
    fx baseline      3,282       42   7.49 +-1.16    2,997
    fx no-floor      2,986       44   8.64 +-1.30    2,688
    fx DIA-NN        4,606       35   4.44 +-0.75    4,404

**Three things it cannot answer**, and the script prints them after every run
so a number cannot be quoted without them:

  * **FDR / FDP acceptance.** Effects attenuate ~5x -- the fragment floor is
    worth 6.91 pp on the full run and 1.15 pp here, against a Poisson sigma of
    1.74. Differences under ~2.5 pp are not differences.
  * **Where we stand against DIA-NN.** The fixture compresses the ratio from
    2.96x to 1.40x and moves the two tools' FDP in OPPOSITE directions
    (DIA-NN 7.42 -> 4.44, ODIA 5.72 -> 7.49).
  * **Retention-time calibration.** The fixture cannot seed its own map -- the
    fit collapses to a slope of 2.35 against 1086.50 -- so the map is supplied
    and that path is untested.

For any of those: `shared/libv2/run_full_v5.sh`, 3 h.

DIA-NN reads the same fixture as `s08_6x60.mzML` (9:58); the Astral fixture
`astral_7x60` (54,533 MS2, 17.9%, 150 windows) runs DIA-NN in 7 s. Rebuild
either with `shared/libv2/build_fixture.sh`.

**Open:** whether a 40% fixture removes the attenuations and at what speedup.
The mechanism predicts the floor effect partly returns, because weak
0-2-fragment precursors would start producing candidates again. One build, two
runs.

---

## Interference correction: the target population is real, three designs are dead, one survives (2026-08-21)

Written to be picked up cold. Everything below is measured on S08 with the 9.6M
`human_v2` library, DIA-NN run on the SAME library (39,149 precursors at
q <= 0.01). Full working in doc/45 and doc/46.

### Why this is worth doing at all

At matched entrapment FDP we report 15,036 against DIA-NN's 39,211 -- 2.4x
short. Everything upstream has been excluded BY MEASUREMENT, so the gap is not
where it was assumed to be:

  * trace cleanliness is MATCHED -- 0.270 against DIA-NN's 0.311 mean pairwise
    fragment coherence, identity-matched fragments, common pair set, median
    paired difference -0.010, and we are better on 41.6% of precursors;
  * m/z tolerance is not too wide -- DIA-NN's own log optimises to 12 ppm
    against our 10, i.e. it extracts through a WIDER window;
  * `Max` aggregation loses to `Sum` at every width tested;
  * mobility window +/-0.025 beats +/-0.0304 and +/-0.050 at a correct centre;
  * library coverage is complete -- 0 of DIA-NN's confident precursors absent;
  * the CLASSIFIER is at its ceiling -- a 5-fold CV supervised GBT trained with
    DIA-NN's own labels as positives, an advantage production cannot have,
    beats the shipped DScore by +1.8% (13,848 against 13,609 at decoy FDR 1%).
    No model, hyperparameter or training change can reach this. The FEATURES
    are the limit.

Partition of DIA-NN's 39,149 through our pipeline:

    accepted by ODIA                             11,773   30.1%
    right peak, correctly picked, scored too low  7,229   18.5%   <- SCORING
    right peak present, wrong one picked            730    1.9%
    candidate formed, none near DIA-NN's RT       2,777    7.1%
    no candidate formed at all                   16,640   42.5%   <- DETECTION

### The population is PRESENT, and this was nearly got wrong

The 18.5% bucket is decoy-like on every co-elution and shape feature
(library_corr 0.040 against decoys' 0.000; xcorr_shape 0.227 against 0.211)
while being positionally correct in retention time. The natural reading -- and
kimi's, in review 46 -- was ABSENCE: `var_ms1_coelution` put them at -0.093
against -0.124 for the 750k precursors nobody claims are present. MS1 sits at
the precursor m/z, a far sparser space than fragment bins, so a present-but-
contaminated peptide should still show its monoisotope. That argument would
retire this whole line of work.

It was wrong, and the reason is instructive. `Ms1Traces::build` matched on the
library's THEORETICAL precursor m/z with no calibration offset while the
fragment axis was centred on the fitted deviation. Fixed in `6b1b2fa`; the run
now reports the centre and the residual against it:

    MS1 mass axis centred on -8.949 ppm (borrowed from the fragment fit);
    median residual against that centre -0.426 ppm

So the MS1 error IS the fragment error, and the old window sat 89% of the way
to its edge -- fine for bright precursors, fatal for weak ones, which is
exactly the population the absence claim was about. Re-measured on the
calibrated axis (`full_v7`):

                              v3 (uncalibrated)   v7 (calibrated)
    DIA-NN + we accept                    0.473             0.982
    DIA-NN, we REJECT                    -0.093             0.151
    not DIA-NN, bulk reject              -0.124            -0.717
    gap (reject - bulk)                   0.031             0.868

They are not absent. They sit 0.868 above the bulk, carrying real MS1 presence
evidence, positioned correctly in RT, with MS2 fragments that do not co-elute.
That is the interference signature.

**Watch the statistic.** The AUC of `var_ms1_coelution` against decoys moved
+0.005 (0.757 -> 0.762) while the within-target separation moved +0.837. AUC
compares targets to decoys and is blind to structure AMONG targets, which is
where this question lives. Reporting the AUC alone -- the habit everywhere else
in the sub-score audit -- would have said the calibration achieved nothing.

### Three designs are dead. Do not re-propose them without new evidence

**1. Per-precursor pruning** -- drop a precursor's least mutually-consistent
fragments, rescore the rest. Prototyped on 998 real targets and 998 decoys:

    corr over ALL fragments (current)   AUC 0.8288   target 0.150  decoy 0.021
    corr over consistent subset         AUC 0.8278   target 0.434  decoy 0.182
    fraction consistent                 AUC 0.7465   target 0.083  decoy 0.000

Pruning triples the target median and triples the DECOY median with it.
Discrimination does not move, because the pruning criterion IS the scored
statistic: a decoy's best-correlating subset flatters the decoy exactly as
much. Selection on the response.

**2. Library-competitor counting** -- flag a fragment contested when another
library precursor in the same isolation window has a fragment within tolerance
whose predicted RT and 1/K0 could co-elute, then score the uncontested subset.
Measured, sweeping the co-elution tolerance (codex was right that the pass-2
extraction window is a resource bound, not a co-elution test):

    RT tol   contested t/d   uncontested/total   corr ALL   corr UNCONTESTED     n
      60 s   1.000 / 1.000          1.9 / 11.7     0.8574             0.8788   312
      10 s   0.667 / 0.750          4.6 / 11.7     0.8496             0.8287   986
       5 s   0.500 / 0.500          6.1 / 11.7     0.8402             0.8298 1,392

At a CORRECT tolerance and a representative sample, dropping contested
fragments makes discrimination WORSE. The apparent +0.021 at 60 s existed only
on the 312-precursor minority that survived a wrong filter. Note also that
decoys are MORE contested than targets at 10 s (0.750 against 0.667), so the
feature partly encodes decoy construction: `mutate` maps all 20 residues onto
`LLLVVLLLLTSSSSLLNDQE` (LibraryGenerator.cpp:37-38), eight residue masses, no
aromatics, no sulfur, no G/A/P, so decoy fragments sit on a different
mass-defect manifold from natural ones.

**3. Mass-tightening survival** -- DIA-NN computes its co-elution sum at base,
0.45x and 0.2x tolerance ("a real peak survives tightening; an interferent
often does not"). All arms centred identically on -9.6479 ppm:

    tolerance                AUC     median target   median decoy
    3.6 ppm               0.7537             0.058          0.010
    12 ppm                0.7658             0.093          0.022
    10 ppm + IM 0.050     0.8103             0.121          0.021

    intensity survival ratio I(3.6)/I(12): AUC 0.5720, against 0.6084 for the
    wide intensity alone and 0.6138 for the tight.

Tightening HURTS. The reason is the useful part: the fragment mass calibration
reports total per-hit scatter of 4.04 ppm, so +/-3.6 ppm is NARROWER THAN ONE
SIGMA of single-hit noise and starves real peaks faster than interference --
zero fraction 19.7% at 12 ppm against 45.3% at 3.6 ppm. The PRINCIPLE is not
refuted, the WIDTH is. **Untested and cheap: tighten to ~6 ppm (above the 4.04
ppm scatter) and keep it ALONGSIDE the 10-12 ppm sum rather than instead of
it.** One extraction on the existing 1,000-precursor subset, minutes.

### What survives: evidence-weighted, cross-precursor, from a previous pass

The two failures share one cause: **library structure says what COULD be
contested, not what WAS contaminated**, and per-precursor statistics cannot
tell the difference because a decoy's own fragments flatter it identically.
Only run evidence separates them, and the only non-circular source of run
evidence is a PREVIOUS PASS.

Sketch: weight each fragment's contamination by the competitors that were
actually IDENTIFIED near that RT and mobility in pass 1, then recompute the
co-elution features in pass 2 on the down-weighted set. It evades both failure
modes -- the weight comes from OTHER precursors' evidence rather than this
one's own correlations, and it counts identified competitors rather than merely
permitted ones. A decoy cannot fake having been identified.

This also reframes the iteration item elsewhere in this file: DIA-NN's ~12
passes are not valuable as iteration per se; interference correction is the
REASON to iterate. Prerequisite, scoped in the plan file and never built:
retain peak groups between passes (~120 B each against tens of kB per trace) so
pass 2 can rescore without re-extracting.

### Constraints any design must satisfy (all learned the hard way)

  * **Never select on the scored statistic.** Failure 1.
  * **Never use a candidate-derived property to define the selection** --
    observed fragment coverage, winning peak width, retained candidate count --
    without a separate selective-inference argument (codex, review 46). Charge,
    library fragment count and precursor m/z are precursor-fixed and safe.
  * **`corr_sum` is a SUM, not a mean.** Removing fragments changes its scale
    and its null; a clean score over 3 fragments is not comparable with one over
    10. Carry the clean-fragment count, the pair count, mean-per-pair and a
    shrinkage term, and an explicit missing value below minimum support --
    never a placeholder dressed as a measurement.
  * **Guard degenerate statistics.** `pearson` over n points is +-1 at n = 2;
    `libraryCorrelation` now reports 0 below four INFORMATIVE (non-zero) points.
    The same trap will exist in any new per-fragment statistic.
  * **Audit target/decoy exchangeability BEFORE looking at identification
    gains** (codex). If the classes differ on the new feature before run
    evidence is consulted, it is unsafe for target-decoy FDR whatever it adds.
  * **AUC is the wrong acceptance metric.** The pruning experiment moved medians
    3x while AUC did not move; the MS1 result moved within-target separation
    +0.837 while AUC moved +0.005. Acceptance is entrapment FDP at MATCHED
    threshold and IDs at DIA-NN's operating point.
  * **Never compare arms at nominal 1%.** Different arms sit at different
    empirical FDPs -- 5.08% against 5.72% for two gates -- so nominal comparison
    compares two different thresholds. This reversed the Gate C verdict once.
  * **Pass-1 counts do not predict final results.** Three instances: the
    mobility seed +43% -> +2%, the gate retest misread as a loss when it is a
    +11.7% win, MS1 calibration +2.4% -> -1.4%.

### Measurements to run first, cheapest first

  1. The ~6 ppm mild-tightening arm above. Minutes, existing subset.
  2. Per-fragment diagnostics on the 10,092 against accepted, matched on
     charge/m/z/RT/intensity: effective fragment count, top-1/2 area share,
     apex-RT MAD (vault: *Testing for interference cheaply*). Distinguishes
     "uniformly weak" from "one or two hijacked fragments with clean siblings".
     kimi calls this the table that decides more than any argument.
  3. Contest-flag precision/recall against ACTUAL decorrelation on true
     precursors -- the load-bearing number the competitor design never had.
  4. Target/decoy competitor-count distributions from the library alone, no run,
     stratified by fragment m/z, iRT and charge, to settle the mutate-alphabet
     hazard. If they differ, prefer composition-preserving decoys
     (pseudo-reverse/shuffle keep the residue multiset) over normalising the
     feature per class.

### Where the code is

  * `src/score/PeakGroupScorer.cpp` -- `coelutionEvidence` (the gate statistic,
    shared by both gate modes), `libraryCorrelation`, the CORR_SUM block
    (~:553-565), the sub-score assignment block (~:1147-1230).
  * `src/extract/ChromatogramExtractor.cpp:1150-1190` -- the match loop and
    aggregation; `LiveSlot` already carries per-cell `ppm_num`/`ppm_den` planes,
    which is what a tightened-tolerance feature would read.
  * `src/extract/Ms1Traces.cpp` -- now takes `ppm_offset` and reports the
    residual.
  * `include/odia/scoring/lda.h:200-280` -- q-values; `:1035-1082` -- the
    per-precursor argmax that a new feature would perturb (winner's curse).

### Related open items

Group-wise FDR was designed and KILLED in review 47 -- both reviewers
independently: it fixes the ALLOCATION of a miscalibrated error rate and cannot
fix its LEVEL, and no stratum is well enough calibrated to rescue (the dominant
12-fragment population sits at 4.88% against a claimed 1%). kimi's framing is
the one to keep: **decoy construction / null-shrinking is the LEVEL lever,
stratification is the ALLOCATION lever, and only the first can ever deliver a
true 1%.** The entrapment homology confound that would have voided the
motivating gradient was audited and is clean -- 0 of 733,780 entrapment
precursors share a stripped sequence with any of 1,325,594 human sequences, in
every stratum.

---

## The 1/K0 anchors are in; what is still open on that axis (2026-08-06)

The 1/K0 calibration now measures at anchors instead of guessing at them, and
the gate passes on its own margin. `collectAt()` visits only the cycle blocks
that hold a scored peak group's apex; the null is the library's own decoys at
the apexes THEIR groups claimed. On S08, with 1,534 confident targets and a
rank-matched null: peakedness **12.41 against 3.75, a 3.31x margin** where 1.25x
is required and where the blind probe managed 1.09x. 23.7% of the mean squared
1/K0 error removed out of fold. The fitted curve matches DIA-NN's observed 1/K0,
which the probe never sees, to within a couple of milli-1/K0 per m/z bin.

End to end on the frozen discriminant: **52.27% -> 53.17% (1393 -> 1417 of
2,665) at a matched, re-measured 1.00% entrapment false rate**; 54.15% at the
frozen threshold, where the entrapment null has itself moved to 1.23%. The gain
sits exactly where the lever was aimed -- the 388 outside-the-cell misses go
15.7% -> 23.2%, and every other bucket moves by less than a point.

What is still open:

- **The stage costs a second sequential decode: 22,440 of 32,210 spectra,
  ~450 s, roughly doubling a single-pass extraction.** Anchors are spread over
  the whole gradient, so "only the blocks with an apex" is 70% of the run. The
  saving available is that within a visited block only the ONE isolation window
  each anchor lives in is needed -- about 2 of 32 spectra per cycle. That is a
  `SpectrumSource::peaks` call pattern change, not an algorithm change, and it
  is worth ~10x if the reader serves sub-ranges without re-decoding a row group
  per spectrum. It does not today (see the mzPeak item below), so this is
  blocked behind the same decode fix everything else is.

- **The gain at 1% is the smallest point on its own sweep**, and the 1%
  entrapment quantile rests on 26 events. Each arm at its own entrapment
  quantile: +1.02 points at 5%, +1.51 at 2%, +0.90 at 1%, +1.16 at 0.5%, +1.62
  at 0.2%. Nothing here says 1% is special; it is where the operating point was
  fixed. Worth re-measuring on a second run before the +0.90 is quoted as the
  number.

- **66 precursors are LOST at the frozen threshold against 116 gained.** A
  recentred window that moves off a precursor whose library 1/K0 was already
  right is the obvious mechanism, and it is not yet measured. The out-of-fold
  scatter says the correction is right on average; it does not say it is right
  for those 66.

- ~~**The width lever is still untaken.**~~ **TAKEN AND CLOSED (2026-08-07):
  widening is worse, keep +/-0.025.** Measured on S08/lib_targets with the
  mobility slope fitted, identifications at 1% FDR: **0.025 -> 1232, 0.035 ->
  1215, 0.050 -> 1167.** Monotonic, and 0.050 gives back nearly the whole gain
  from the slope (1165 with no mobility-linear term at all).

  The residual argument that pointed the other way -- corrected SD 0.025, so
  +/-0.025 is ~1 sigma and rejects 28.3% of target anchors -- was right about
  the width and wrong about the remedy. Those 28.3% are the interference-prone
  tail; admitting them costs more than they bring, which is exactly the ~8.5x
  same-window mobility this window exists to exclude. The fix for a trend is to
  model the trend, not to widen the gate until the trend fits through it.
  Recorded at `ChromatogramExtractor.h`'s `precursor_im_window`.

- ~~**The production path has not been measured end to end on S08.**~~ **RUN
  (2026-08-07).** `-ion_mobility_calibration anchors` over the two-pass
  workflow, anchors from ODIA's own pass-1 scorer, S08/lib_targets: pass 1
  offers 729 anchors at q<=0.01 against a null of 729 decoys, the gate passes,
  and pass 2 reports **1232 at 1% FDR**. The stage works on its own anchors
  without an external scorer.

  Two caveats it exposed, both open below: charge 3 is left uncorrected, and
  the run bootstraps its own iRT map (no `-irt_slope`/`-irt_intercept`), so it
  warns that pass 1 extracts at approximately the wrong retention times. Both
  arms of the slope comparison shared that, so the differential is clean, but
  the absolute number is not a tuned-iRT number.

- **Charge 3 is never corrected on S08/lib_targets: 93 anchors against a
  minimum of 120.** The mobility slope is fitted per charge on purpose --
  another charge's offset is not this charge's answer -- but that means the
  whole correction, slope included, silently does nothing for charge 3. Charge
  3 is ~20% of the library and its measured residual is the WORST (offline:
  constant +0.0335 and slope -0.113, against +0.0019/-0.096 for charge 2), so
  the charge that most needs the correction is the one that cannot reach the
  anchor count for it. Options, in order of preference: pool charges for the
  SLOPE only while keeping the constant per charge (the slope is a property of
  the CCS->1/K0 conversion, which is shared, so this is physically justified in
  a way that pooling offsets is not); or lower the threshold with a widened
  confidence requirement. Do not simply lower the minimum.

- **The mobility slope is confirmed but its magnitude is not pinned down.**
  Fitted -0.0709 and -0.0748 per 1/K0 on the run's own 371-378 anchors, against
  -0.096 (charge 2) and -0.113 (charge 3) estimated offline on ~22k anchors,
  and -0.128/-0.109 from outlier-resistant bin medians. Same sign and order
  throughout, but a ~30% spread. Since this is a scale error in the
  `ccs_to_mobility.py` coefficient (1037.1902), the right fix is upstream and
  one-off: fit the coefficient properly against Mason-Schamp rather than
  re-deriving a slope per run. Worth doing once the Astral/second-run
  measurement says whether the slope is instrument-stable.

---

## Blocking someone else

- **[you] `kimi` and `codex` CLIs are unavailable, so the adversarial reviews
  are not the ones you asked for.** Neither binary is installed, and there is no
  `node`/`npm`/`npx` to install them, nor any `MOONSHOT`/`KIMI`/`OPENAI`
  credential in the environment — so this is a credentials problem, not just a
  missing binary. Overnight reviews are being done with independent subagents
  prompted to refute instead, which is a weaker substitute: same model family,
  so correlated blind spots. Install the CLIs and provide keys and I will switch.

- **[decided 2026-08-03, user] Phase 2 proceeds on the slow mzPeak reader.**
  The decode fix is being worked on elsewhere; ODIA does not wait for it. The
  extractor is written against `SpectrumSource` (D10), so when the batched
  decoder lands it is a substitution, not a rewrite.
  **Do not let this be forgotten:** every extraction timing measured before the
  replacement lands is a measurement of the reader, not of ODIA, and must be
  re-taken afterwards. `doc/05-mzpeak-batched-reader-handoff.md` is the spec.

- **[for the reader author, not us] mzPeak decode diagnosis.** At `f93f938` the batch API exists but does not amortise:
  284.4 ms/spectrum via `get_spectra_batch`, against 0.082 ms/spectrum if each
  Parquet row group were decoded once and its 873 spectra served from it.
  Measured, not estimated -- a row group decodes in 0.07 s with pyarrow, and a
  full pass over `12_80` would take 1.5 s, beating mzML's 3.5 s. The reader
  decodes a whole row group per spectrum and keeps ~0.1% of it. Everything that
  touches peaks inherits this, including `extract_ion_chromatogram`.

- **[superseded, kept for the measurement] mzPeak peak decoding is ~1000× too slow.** ~277 ms per
  spectrum against ~0.27 ms for OpenMS parsing mzML. Full requirements in
  `03-mzpeak-streaming-requirements.md` (R1, R2, R5 are the blocking set).
  Phase 2 cannot start until this moves.

- **[decided 2026-08-03, user] `.oswpq` writing is not needed.** Dropped from
  scope. Reading stays, and accepts either intensity width. The note below is
  kept only because the reader's behaviour depends on it.

- **[was: you] `library_intensity` `float64` vs `float32` in `.oswpq`.** The
  predecessor's `float32` change is a breaking on-disk change that was never
  upstreamed, and we cannot apply it without modifying OpenMS. Proposal on the
  table: read both, write upstream-compatible `float64` by default with
  `float32` behind an option. **The reader now does accept both**, verified
  against fixtures of each width, so only the writer is still blocked on this.

- **[resolved 2026-08-04] GPU access.** These nodes authenticate with
  **Kerberos**, not keys -- `ssh spock` works from an interactive login because
  PAM puts a TGT in the session keyring, which an agent session cannot reach.
  A file-based ticket (`KRB5CCNAME=FILE:<path> kinit`) hands one over without
  creating a standing credential, and it expires on its own. The CUDA path is
  now measured; see "GPU inference measured on an H100" below.

- **[you] GPU 0 is unusable on spock and data.** `cudaSetDevice(0)` returns
  error 46 with both devices idle, in Default compute mode, and `/dev/nvidia*`
  world-writable. On `data` device 1 works, so `CUDA_VISIBLE_DEVICES=1` is a
  workaround; on `spock` **both** devices fail, leaving it with no usable GPU.
  This is an infrastructure fault, not an ODIA one -- worth a ticket.

---

## Correctness, must be settled before results are trusted

- **Validate the PeptDeep modified-peptide path against AlphaPeptDeep reference
  values.** The encoding is reconstructed from upstream sources and documented in
  `04-peptdeep-encoding.md`, and the unmodified path can be validated against
  OpenMS exactly — but nothing here can confirm the modification path. A wrong
  `mod_x` yields plausible, quietly wrong intensities. **No library generated
  with modifications should be trusted until this is done.**

- **Isotope element mapping.** AlphaPeptDeep's element list contains `2H`,
  `13C`, `15N`, `18O`; OpenMS writes isotopes as `(13)C`. Without an explicit
  translation, isotope-labelled modifications land silently in the `?` bucket.
  Adversarial review confirmed these four renames are the *only* mapping needed:
  all 40 element symbols in the shipped `unimod.xml` are in the 109-element list.

- **Mod feature counts are signed.** 600 of 2859 rows in the UniMod-derived table
  have negative counts (`Deamidated@N` is `H(-1)N(-1)O(1)`), so an unsigned type
  or an `abs()` on the OpenMS side is silently wrong. Measured: RT 0.5231 signed
  against 0.5121 with `abs()`.

- **Batches must be length-homogeneous.** Trailing padding is not inert -- index 0
  is one-hot encoded and no model applies a padding mask. The same peptide padded
  to `seq_len` 13/20/30 predicts RT 0.8196/0.4737/0.2714. Group by encoded length
  and run one batch per length.

- **Determinism harness.** Permutation invariance and thread invariance at
  1/8/64 threads with OpenMP linked, per the cross-cutting invariants. Not yet
  built; it needs to exist before any A/B measurement is believed.

---

## From the Phase 1 adversarial review — still open

Four criticals from that review are fixed (decoy fragment masses, the
target/decoy round-trip merge, non-contiguous row detection, Parquet numeric and
string type coverage), plus the `toFixed` domain, the unchecked output stream and
the exit-code-6-on-success. These remain:

- **N-terminal modifications break the decoy tokeniser.** `(UniMod:1)PEPTIDE`
  has its modification name tokenised as residues. **Reachable from the supplied
  DIA-NN fixture**, which holds 441 such precursors — 428 targets (1.3%) get no
  decoy as a result. The count is now reported rather than silently dropped, but
  the tokeniser still needs fixing.
- **The Parquet reader still cannot open a library past 2 GB of characters.** The
  `large_string` handling casts down to `arrow::utf8()`, which has the 32-bit
  offsets `large_string` exists to avoid, and the plain-`string` path fails in
  `CombineChunks`. Both emit a message blaming the file for a missing column that
  is present. `LargeStringArray` needs to be a first-class case.
- **`sortByPrecursorMz` doubles peak memory** (+483 MiB, +39% on the human
  proteome) by materialising complete copies before moving. An in-place
  permutation or a block-streamed rebuild avoids it.
- **The contiguity check is order-dependent.** Two blocks of one precursor that
  are *adjacent* still merge silently, discarding the second block's RT, IM,
  precursor m/z and protein group; only non-adjacent repeats are detected.
- **Null cells are indistinguishable from zero** for `Decoy`, `Precursor.Charge`,
  `Fragment.Series.Number` and the string columns — a null `Decoy` column makes
  every decoy a target, the same effect the BOOL case was added to fix.
- **Unusual loss labels are not round-tripped.** `H2O+H2O`, `CH3SOH`,
  `H3PO4+H2O` all parse to `LossType::Other` and are written back as `other`, so
  the label is lost and the transition is dropped from decoys. Either carry the
  label as an interned string or refuse the input.
- **A decoy's `Product.Mz` can no longer be checked against its own row.** Since
  the decoy stores the target's sequence (DIA-NN convention), reproducing its
  fragment m/z requires applying the mutation table — which has become an
  unwritten part of the file format. Defensible, but it should be stated in D7.
- **The `_decoy` id suffix deviates from DIA-NN**, whose decoy `Precursor.Id` is
  exactly `sequence+charge`, identical in form to a target's. Ours has to differ
  because we emit both rows; worth documenting as an intentional divergence.
- **Contiguity-check hash collisions are a hard failure.** 4.3e-7 at 4 M
  precursors (fine), but 1.7e-4 at the 78.6 M design scale — one load in 6,000 —
  and the failure is a deterministic `throw` whose message sends the user to
  re-sort an already-sorted file.
- **`generate()` interns before deciding whether a precursor survives**, leaving
  5,532 permanently unreferenced arena entries on the human proteome and
  over-reporting the distinct-string count by that much.
- **D5's mandated warning is not implemented** — the reader should warn once when
  an input's transition names disagree with what ODIA would synthesise. The
  accepted round-trip risk is currently unmitigated.
- **`Decoy` given as `True`/`False` parses to 0**, and `toLong` cannot
  distinguish absent from unparseable.
- Smaller: `char buf[64]` truncates long numeric fields and `strtod` accepts
  partial parses (`500,1` becomes 500); no `static_assert` that the two mutation
  tables are the same length; `loadTSV` never reserves; a handful of
  `pseudo_reverse` decoys are identical to their target for palindromic
  prefixes.

---

## Library generation is ~5x slower than it needs to be (measured 2026-08-04)

Library generation is 53 min for the human proteome and dominates end-to-end
cost by ~40x. Two causes, both in how ONNX Runtime is driven, neither in the
models themselves:

**1. `-threads` never reaches inference.** `LibraryGenerator` constructs all
four predictors as `PeptDeepPredictor(model_path, prefer_gpu)` — two arguments,
so `intra_op_threads` takes its default of 0 and `SetIntraOpNumThreads` is
never called. ONNX Runtime then sizes its pool to the whole machine (128 cores
here) regardless of what the user asked for. The tool's `-threads` flag governs
digestion and decoy construction and nothing else, which is ~7% of the phase.
Sites: `LibraryGenerator.cpp:343` (RT), `:377` (MS2), `:561` (iRT), `:631` (CCS).

**2. Intra-op parallelism is the wrong axis for these models.** PeptDeep's MS2
network is recurrent; its per-op tensors are too small to spread across many
threads. Measured on 16 pinned cores, one session, MS2 throughput against
intra-op thread count:

| intra-op | peptides/s |     | processes x threads | peptides/s |
|---------:|-----------:|-----|--------------------:|-----------:|
|        1 |      168.7 |     |          16 x 1     |   **1760** |
|        4 |      338.7 |     |           8 x 2     |     1253   |
|        8 |  375.6 (peak) |  |           4 x 4     |      815   |
|       16 |      317.9 |     |           2 x 8     |      505   |
|       32 |      231.2 |     |                     |            |

Intra-op saturates at 8 threads for 2.2x and then *degrades*. Data parallelism
— N independent sessions, one thread each — is monotonically better and beats
the best intra-op configuration by **4.7x on identical cores**. It keeps
scaling: 8 cores 1063/s, 16 cores 1480/s, 32 cores 2716/s, 64 cores 4585/s
(per-core falls 133 -> 72 as memory bandwidth saturates, but aggregate climbs).

Against the 897 peptides/s the v5 run actually achieved, 64 data-parallel
sessions measure **5.1x**. Projected: MS2 2372 s -> ~465 s, and library
generation 53 min -> **~14 min** if RT and CCS gain proportionally (they share
the same defect; not separately measured).

**The fix**: give each worker thread its own `Ort::Session` with
`intra_op_threads = 1` and slice the block across them, replacing the single
shared session. Cost is memory — 129 MB RSS per standalone process, though
in-process the marginal cost is the arena plus a 16 MB weight copy, so a 64-way
split should land near 3-5 GiB against today's 2.1 GiB peak. Worth capping the
session count independently of `-threads` for that reason.

Both measurements were taken while another user held ~114 of the 128 cores, so
the ratios are back-to-back under equal contention and the absolute rates are
floors.

## Benchmark ODIA against other DIA tools via ProteoBench [user, 2026-08-04]

Explore <https://proteobench.cubimed.rub.de/>, download the DIA benchmark, and
design a workflow that compares ODIA to the other DIA tools on it.

Why this is worth doing properly rather than quickly: every comparison made so
far has been against DIA-NN on S08, with DIA-NN's own library as the reference
and DIA-NN as the search engine. That measures agreement with one tool on one
file, and it cannot distinguish "ODIA is good" from "ODIA resembles DIA-NN".
ProteoBench supplies a defined ground truth and a published leaderboard, which
is the first thing here that could falsify a claim rather than confirm it.

Points to settle when this is picked up:

* Which module — the DIA modules differ in organism mix and in whether the
  quantitative ratios or the identification counts are the scored quantity.
* ODIA does not yet do quantification, so the honest first submission may be
  identifications only. Check what the leaderboard requires before building to it.
* The entrapment measurement already on this list is a prerequisite, not a
  parallel task: a leaderboard position computed on an uncalibrated FDR is
  worse than no position.
* Submission means publishing a result under our name. Do not submit anything
  without asking first; running the benchmark locally needs no such permission.

## [you] The pyProphet/XGBoost scoring code is not reachable (2026-08-04)

The decision is to use the pyProphet-like XGBoost scoring already written at
`/Users/kohlbach/Claude/mzPeak/OpenDIAlyzer`. That is a **macOS path**; the
nodes this project builds on are Linux and have no `/Users` at all, so the code
cannot be read, let alone called. Copy it to `/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/`
(or anywhere on Ceph) and Phase 3's classifier is unblocked.

Upstream pyProphet 3.0.15 with XGBoost 3.2.0 *is* installed at
`/ceph/ibmi/abi/oliver/envs/pyprophet` -- a useful cross-check, but not the
code the decision named, so it is not a substitute without saying so.

Peak-group detection and the sub-scores do not depend on this and proceed
meanwhile; see `doc/07-scoring-plan.md`.

## Quantise chromatogram intensity to one byte

The shared retention-time axis landed; the intensity encoding did not. Storing
each point as a uint8 log step against a **per-transition** float32 scale is a
further 4x on the point arrays.

Recorded because the objection raised against it was wrong: the concern was
that quantisation would cost precision in the weak-signal regime currently
under investigation. It does not, because the scale is per transition rather
than global -- 255 log steps over one transition's own two-decade range is
~2.7% per step, well below the noise on any point that matters. The reference
implementation uses exactly this (`odia_chromstore.h`, `Quantised8Log`).

Reserve 0 for exact zero rather than for the smallest representable value: an
absent measurement and a very small one are different claims, and a
chromatogram is mostly the former.

## Implementation, unblocked

*Still uncovered by any fixture, from the mutation-testing review: neutral loss
on 3+ charged fragments; a phosphorylated modification (the checker's MOD table
has no UniMod 21, so it cannot yet be pointed at a phospho library, which is the
library type `lossMass`'s own comment cites); `pseudo_reverse` decoys, which no
test invokes and which `check_invariants.py` cannot validate since its `mutate()`
implements only the substitution table; and a precursor split across two
**adjacent** blocks.*

*Fixture coverage is otherwise addressed: `scripts/make_adversarial_fixtures.py` plus
`test/check_invariants.py` give 20 tests over losses, modifications, decoy
collisions, unusual Parquet types, nulls and malformed input, with masses
verified against a residue table independent of OpenMS and pinned by a
self-test of literal values. Mutation testing: 21 injected defects, all of the
ones re-tested after repair are now caught.*

- **iRT prediction is wired into library generation.** Predicted once per
  distinct modified sequence rather than per precursor, since the RT model has
  no charge input: 1.25 M predictions instead of 4.0 M on the human proteome,
  206 s, and 232 s for the whole run at 1.9 GB peak. Unpredictable peptides are
  left NaN and counted, never given a made-up value.
- **iRT prediction is done and validated**: `PeptDeepPredictor` runs the RT model
  with CUDA attempted and CPU fallback, and its predictions match the
  independent Python reference to 1e-5 on mixed-length batches including
  residue, N-terminal and C-terminal modifications. Still to write: MS2 and CCS,
  which add `charges` (×0.1), `nce` (×0.01) and a rank-1 `instrument_indices`,
  and the wiring that replaces the placeholder intensities and iRT in
  `LibraryGenerator`.
- Earlier note, retained: the **encoder is validated**: `PeptDeepEncoder` agrees
  exactly with the independent Python reference on residue, N-terminal and
  C-terminal modifications, multiple modifications, signed counts and
  `aa_indices`.
- **Predicted iRT is the model's normalised output, and nothing says so in the
  file.** Measured range over 109,864 human sequences: -0.039 to 0.940, with
  14.6% above 0.85 (saturation) and 1.24% negative. Nothing converts it to iRT
  units or minutes, and the invariant checker treats a literal `0` as missing
  while accepting -0.039. Decide the unit and declare it.
- **Prediction costs 3 h 07 m of CPU time** (11,267 CPU-seconds at 4688% for
  206 s wall) for one proteome library. Fine here; worth knowing before it runs
  anywhere metered, and an argument for the GPU path.
- **ONNX Runtime's arena accounts for +610 MiB of the +688 MiB** prediction adds
  at proteome scale, and it never shrinks. It is bounded by `MAX_BATCH_ROWS` x
  encoded length, so `-max_peptide_length` scales it linearly with nothing
  capping or testing that.
- **Peptides longer than the model was validated on are unreachable only by
  accident**: `precursor_mz_max = 1200` at charges 2-3 caps peptides at ~32
  residues before a precursor exists. `PeptDeepEncoder` has no length guard, so
  widening the m/z window would feed the model lengths it has never seen.
- **The library records neither the model nor the execution provider used.**
  CUDA and CPU differ in the last bits, so a written library is currently
  unattributable.
- **Predictions depend on batch composition at the last bit.** Within one build,
  repeated calls are bit-identical, but batch size moves the result by ~1 ULP
  (1.5e-8) and the two ONNX Runtime builds compared in the tests differ by up to
  2.7e-7. `MAX_BATCH_ROWS` pins the batch size, which bounds it, but D1 asks for
  determinism as a precondition and this is not zero. Decide whether that
  matters at the iRT-window scale.
- **The CUDA path is still untested** — `spock`/`data` remain unreachable, so
  only the CPU branch has ever run. `Ort::GetAvailableProviders()` is queried
  before attempting, so a CPU-only build no longer logs a scary provider-load
  error, but that also means the CUDA branch is skipped rather than exercised
  here.
- **A doubly-modified residue is under-encoded.** OpenMS's `AASequence` keeps
  only the last modification on a residue, so it is lost before the encoder sees
  it; AlphaPeptDeep would accumulate. Decide whether to detect and refuse, or to
  carry modifications outside `AASequence`. The reference encoder and RT
  predictor (`test/peptdeep_reference.py`) are in place as the validation oracle,
  written from the spec rather than from the C++; the element list is generated
  from the authoritative yaml into `data/peptdeep_mod_elements.txt` and
  `include/odia/PeptDeepElements.h`. The C++ side is still to be written.
- **MS2 prediction works and matches the reference**; still to do is *using* it:
  `LibraryGenerator` caps fragments by descending m/z because it has no
  intensities yet. Wiring it in means mapping the `[positions, 8]` output onto
  b/y ions at charge 1-2, choosing a cap by predicted intensity, and deciding
  what NCE and instrument to assume when the caller does not say.
- **NCE and instrument are guesses.** They default to 30.0 and "QE", and nothing
  derives them from the data. They materially change the spectrum, so a library
  generated for one instrument is not right for another.
- CCS prediction, for the ion-mobility column.
- `.oswpq` **write** (still gated on the `float64`/`float32` decision above).
- Phase instrumentation: wall, CPU, RSS, `mallinfo2`, and node load per run,
  charging un-phased time to the preceding phase.
- mzPeak chromatogram writer, using the fork's writer.

---

## Decisions taken 2026-08-03 (user)

- **Development proceeds against a DIA-NN library**, not ours, so that Phase 2
  work is not confounded by library differences. Ours stays the deliverable;
  it is simply not the variable under test while the extractor is built.

- **Decoy construction is deferred, not solved.** Ours cost 14,082 precursors
  (45%) when DIA-NN searched with them instead of building its own. They keep
  the target's m/z, iRT and intensity pattern with only fragment masses moved,
  which is DIA-NN's published design, so the fault is more likely in how the
  intensity pattern is copied than in the mutation scheme. Until this is
  understood a library written for another engine should probably carry no
  decoys at all. **This is an FDR question, and it blocks Phase 3, not Phase 2.**

- **NCE and instrument stay guesses for now** (30.0 and "QE"). Nothing derives
  them from the data and they materially change every predicted spectrum, so a
  library generated for one instrument is not right for another. Revisit when
  the extractor can measure the mismatch rather than assume it.

---

## Phase 1 measured against DIA-NN on S08 (2026-08-03)

Both libraries predicted from the same FASTA with the same digest settings;
DIA-NN searched the same diaPASEF run with each, everything identical but the
library.

| library | precursors | proteins |
|---|---:|---:|
| DIA-NN's own | 37,170 | 5,255 |
| ODIA, first attempt | 17,478 | 3,263 |
| + targets only (DIA-NN makes its own decoys) | 31,560 | 4,776 |
| + ion mobility from predicted CCS | 33,732 | 5,021 |
| + fragment charge up to precursor charge | 34,609 | 5,053 |

**Yield per library precursor is identical**: DIA-NN 1.732%, ODIA 1.736%. The
assays are as good, one for one; every remaining difference was coverage.

Still open from this:

- **Our decoys cost 14,082 precursors, 45% of the total.** Shipping them made
  DIA-NN use ours instead of building its own, and ours are too target-like to
  separate: same m/z, same iRT, same intensity pattern, only fragment masses
  moved. That is DIA-NN's own published design, so the fault is more likely in
  how the pattern is copied than in the mutation scheme. **Until this is
  understood, a library written for another engine should probably not carry
  decoys at all** -- which is a decision, not a default to pick silently.

- **Retention time is the weakest remaining signal.** DIA-NN set an RT window of
  2.29 min for our library against 1.44 for its own, so our iRT is about half as
  predictive after calibration despite a rank correlation of 0.992 with theirs.
  Worth attacking next; it is a scoring dimension, not a coverage problem.

- **Fragment selection overlaps only 0.60** with DIA-NN's, and base peaks agree
  62% of the time, while intensities on the fragments both keep correlate at
  r = 0.861. So the models mostly agree and the caps diverge. Not obviously
  wrong, but not understood either.

- **CCS agrees with DIA-NN's ion mobility to 1.43% median** over 1,873,932
  precursors, with no charge or length dependence and a split-half held-out
  median of -0.004%. Both predictors are sound; see `test/compare_ccs_to_diann.py`.

---

## RT and m/z calibration are independent -- measured 2026-08-03

Asked whether the m/z tolerance has to be recalibrated whenever the retention
times are. **No.** Across three libraries differing only in their RT column,
DIA-NN's achieved mass accuracy is identical to six significant figures while
RT accuracy varies 2.2-fold:

| library | RTPredAcc | median MS1 | MS1 corrected | median MS2 | MS2 corrected |
|---|---:|---:|---:|---:|---:|
| v3 stock RT | 0.361 | 4.85479 | 1.16167 | 1.48522 | 1.06719 |
| v5 fine-tuned | 0.211 | 4.85479 | 1.16167 | 1.48522 | 1.06719 |
| v4 fine-tuned | 0.164 | 4.85479 | 1.16167 | 1.48522 | 1.06719 |

Mass error is a property of the instrument and of the library's m/z values;
recalibrating retention times changes neither. The only coupling observed is
second-order: DIA-NN's *chosen search tolerance* moved between 7 and 10 ppm
across the variants and not monotonically in RT accuracy, which is its coarse
parameter grid reacting to a different identification set rather than a change
in the underlying mass error.

**The finding that matters for ODIA is that it calibrates neither.**

- `ChromatogramExtractor::Options::fragment_ppm` is a fixed **20.0**. The
  measured median MS2 mass error on this run is **1.49 ppm**, and DIA-NN
  optimised its own search tolerance to 7-10 ppm. Ours is roughly 2-3x wider
  than the data supports, and every extra ppm admits proportionally more
  interfering peaks into every transition's chromatogram. It is a guess, and
  nothing currently measures it.

  **Done, 2026-08-05** (`include/odia/MassCalibration.h`). The run's own
  fragment mass error is now measured before extraction and the window is
  centred on it, gated so that a run with no measurable error is a no-op. The
  finding that mattered was not the width but the CENTRE: on S08 the axis is
  about -10 ppm out, and a window narrowed about zero throws away more than it
  saves. `-fragment_ppm` therefore no longer has a fixed default at all -- it is
  10 ppm when the calibration centres the window and 15 when it cannot.
- `rt_low`/`rt_high` are caller-supplied, global, and the extractor never reads
  a retention time from the library at all (see section 4a of
  `doc/06-rt-refinement-plan.md`).

Both should be derived from the run. Because they are independent, they can be
derived independently -- a joint or iterative calibration is not needed, which
is worth knowing before someone builds one.

---

## Phase 2 first slice, measured 2026-08-03

Extraction runs end to end: library -> window assignment -> one forward pass ->
chromatograms. On `12_80.mzpeak`, 200 precursors / 2,400 transitions over a 60 s
retention-time slice (797 spectra):

| | |
|---|---|
| points | 176,268, 14.9% non-zero |
| **decode** | **262.03 s** (328.8 ms/spectrum) |
| **match** | **0.02 s** (0.021 ms/spectrum) |
| memory | 1.4 MiB |

**The reader costs 13,000x what the matching costs.** ODIA's own work is 0.006%
of the runtime. With the batched decoder specified in
`05-mzpeak-batched-reader-handoff.md` the same extraction would take about
0.09 s rather than 262 s, and matching would become the dominant term -- which
is the point at which optimising ODIA starts to be worth anything.

Consequences to keep in view:

- **No timing taken through the current reader means anything about ODIA.**
  Recorded here so it is not quoted later as an ODIA benchmark.
- **The retention-time range must be resolved to a spectrum index range before
  decoding**, not filtered afterwards. Filtering afterwards read the whole run
  and discarded most of it -- an hour instead of seconds. Fixed, and worth
  remembering as the shape of mistake this reader punishes.
- **Still missing before this is useful**: an iRT-to-RT calibration, so a
  precursor is extracted over a window around where it should elute rather than
  over the whole run; and an output format. Without the first, a whole library
  over a whole run is billions of points.

---

## RT fine-tuning: the controlled result (2026-08-03)

Five searches of S08, everything identical but the library's RT column.

| library | RT source | RTPredAcc | DIA-NN's window | precursors | proteins |
|---|---|---:|---:|---:|---:|
| DIA-NN's own | DIA-NN | 0.2133 | 1.4397 | 37,170 | 5,255 |
| ODIA v3 | stock PeptDeep | 0.3610 | 2.18905 | 35,296 | 5,109 |
| ODIA + DIA-NN's RT | DIA-NN, pasted | 0.2189 | 1.36581 | 35,523 | 5,130 |
| **ODIA v5** | **fine-tuned on 500** | **0.2107** | -- | **35,131** | **5,091** |
| ODIA v4 | fine-tuned on 23,179 | 0.1642 | 1.05543 | 37,466 | 5,295 |

**v4 is contaminated and must not be quoted.** 62.5% of the peptides it
identified had their observed retention time in its own training set. It is
train-on-test at the run level, and the +296 over DIA-NN is memorisation.

**v5 is the honest one** -- 500 training peptides, 1.5% of the run's
identifiable peptides -- and it settles the question:

**RT accuracy improved from 0.361 to 0.211, reaching DIA-NN's own 0.213, and
identifications went DOWN by 165.**

That reproduces the earlier RT-swap ablation (+227, also noise) with an
independent mechanism, and it means the original conclusion was right: **on this
run, retention-time accuracy does not buy identifications.** Two libraries with
equal RT accuracy -- v5 at 0.2107 and DIA-NN at 0.2133 -- differ by 2,039
precursors, so the remaining gap is somewhere else entirely.

**Where "somewhere else" is** remains the fragment-charge and intensity question
recorded below: we emit 78.7% singly-charged fragments where DIA-NN emits 72.4%,
and the precursors we miss are enriched in charge 3.

**What fine-tuning is still for:** the extraction window, which v4 measured at
1.055 min against v3's 2.189. That is a Phase 2 compute and feasibility argument
(`doc/06-rt-refinement-plan.md` section 4a), not a sensitivity one -- and it now
has an experiment behind it rather than an assertion.

---

## Retention-time fine-tuning is integrated (2026-08-03)

`scripts/finetune_rt.sh <report.parquet> <outdir> [n]` fine-tunes the RT model on
a run's own identifications and exports ONNX that ODIA reads with `-rt_model`.
No C++ changed, and torch stays out of the runtime -- OpenMS's own exporter
already emits the input names ODIA expects.

Two things the integration fixed over the prototype:

- **Modifications resolve through alphabase's own table**, 1,524 UniMod ids
  against the prototype's five, with the residue disambiguating candidates and
  a refusal rather than a guess when it cannot.
- **A provenance sidecar** records the model's SHA-256, the identifications it
  was tuned from, the peptide count, and the minute range defining the scale it
  predicts on -- which nothing downstream can recover otherwise.

Still open, and both are in `doc/06-rt-refinement-plan.md`:

- **Method choices settled by measurement**, not preference:
  - *Direct target beats residual learning ON S08, which is one file.* Both
    mechanisms are implemented and selectable (`--method direct|residual`)
    precisely because a single benchmark should not decide it. Measured here:
    direct 0.429 against residual 0.635 at n=500 / 500 epochs, and 0.415 against
    0.591 at n=600 / 60 epochs. Residual IS a real improvement over the stock
    model (42%), just a smaller one than direct (59%). The likely reason is that
    retargeting pretrained weights at a small centred residual fights the
    initialisation; freezing the trunk was not tried and might reverse it.
    `--evaluate` makes every run report which method won on its own data.
  - *Data beats optimisation steps.* 500 peptides at 500 epochs reaches 0.429;
    10,000 peptides at 40 epochs reaches 0.359. Twelve times the gradient steps
    does not close what more peptides close easily.
  - *40 epochs is about right.* At n=2000, going to 120 buys 4.7% while the
    train/held-out gap doubles from 0.097 to 0.198 -- overfitting onset without
    the held-out curve having turned up yet.

- **The library does not carry the model's identity.** The sidecar sits beside
  the model, not beside the library, so two libraries built with different
  models still look identical. That is the remaining half of provenance.
- **Feeding this from ODIA's own first pass closes a loop.** The script cannot
  tell an external search from its own and records the source rather than
  judging it. The entrapment measurement is the precondition.

---

## Retention time — accurate diagnosis, and it does NOT cost identifications

Researched 2026-08-03. **The framing that sent this investigation was wrong,
and the correction matters more than the investigation.** RT was called "the
most likely source of the remaining 5%". It is not.

**The controlled ablation.** ODIA's library with DIA-NN's RT column pasted in,
everything else identical:

| library | RT column | DIA-NN's window | precursors | proteins |
|---|---|---:|---:|---:|
| DIA-NN's own | DIA-NN | 1.4397 | 37,247 | 5,255 |
| ODIA v3 as shipped | PeptDeep raw | 2.18905 | 35,796 | 5,109 |
| **ODIA + DIA-NN's RT** | DIA-NN | **1.36581** | **35,688** | **5,130** |
| ODIA + linear iRT calibration | 152.2356x - 39.2322 | 2.18905 | 35,831 | 5,088 |

A *perfect* RT column buys **-108 precursors and +21 proteins**. The gap to
DIA-NN survives intact. The window is a faithful measure of RT accuracy -- the
peptide sets are identical, so library size cannot be the confound -- but RT
accuracy is not what is costing identifications on this run.

**Implemented 2026-08-03.** The calibration is applied, refitted at run time
from `data/irt_standards.tsv` rather than pinned, and verified against OpenMS's
own fixture (`peptdeep_irt_peptides_predicted.csv`), which carries both the raw
ONNX prediction and the calibrated iRT for the 11 standards. ODIA reproduces the
published line to six digits: 152.2356 x raw - 39.2322, worst standard off by
8.7570 iRT. The invariant checker no longer demands the library's RT equal the
model's output -- it recovers the affine map from the data and requires it to be
the one the standards give, which catches both a missing calibration and a
dropped intercept.

**The iRT calibration is a no-op for accuracy, bit-exactly.** Applying the AlphaPeptDeep
Biognosys linear calibration produced an *identical* window, 2.18905. This is
structural, not luck: DIA-NN fits its own monotone calibration, so any monotone
reparametrisation of our iRT axis cannot change anything. Worth applying anyway
for **unit hygiene** -- a field named `irt` holding a training-gradient
coordinate is the same silent-units hazard the CCS decision exists to avoid --
but it must never be described as an accuracy fix. If applied, recompute the
constants at runtime from the 11 Biognosys peptides; they are checkpoint-specific
and would rot silently on a model swap.

**The real RT defect is output saturation, and it is unfixable by calibration.**
Local calibration slope by decile of our predicted iRT, normalised to the median:

```
ODIA    1.05 0.97 0.97 0.96 0.90 0.95 0.95 1.22 5.11 11.88
DIA-NN  0.89 1.04 1.03 0.98 1.09 1.01 1.02 0.97 0.97 0.84
```

The top two deciles are compressed 5x and 12x; 12.8% of the library piles into
the single bin 0.85-0.90. Confirmed against the model directly with synthetic
ladders: `G14` 0.016 -> `G9L5` 0.876 -> `L14` 0.936, so five leucines spend 0.89
of the range and the next nine spend 0.06. This is the standard artefact of a
min-max-normalised training target. The information is gone, not misplaced,
which is exactly why the calibration experiment came back null.

Worst classes, all one axis seen four ways (long, hydrophobic, charge 3, late):
length 25-30 is 2.06x DIA-NN's spread, late RT 1.73x, GRAVY top decile 1.67x,
charge 3 1.55x. Nothing systematic in modifications or terminal residue.

**Where the effort should go instead (R1).** With RT equalised we still miss
4,937 of DIA-NN's precursors while finding 3,378 it does not. The missed ones
have **median quantity 0.40x** the shared ones and are **enriched in charge 3
(37.0% vs 24.5%)**. Libraries are structurally matched (12.00 vs 11.92
fragments/precursor, b/y 29/71 both) but ODIA emits **78.7% singly-charged
fragments against DIA-NN's 72.4%**. That is the thread to pull: fragment
selection for multiply-charged precursors, not retention time.

**Where RT still pays (R3).** Extraction cost. A 1.6x narrower window is 1.6x
less XIC per precursor -- ~128 cycles against ~80 on this run. That is a direct
constant factor on Phase 2's runtime and peak memory, which is the one place
this still matters.

**Deferred, with reasons.** Fine-tuning the RT model on a run's own
identifications is the only thing that can touch the saturation (Zeng et al.,
Nat Commun 13:7238, report R2 0.927 -> 0.986 from 500 peptides). But R1 says the
identification benefit is ~0, and the self-contained path costs a torch
dependency ODIA does not have. DIA-NN can do it in minutes with bundled libtorch
(`--tune-rt`), which would make ODIA depend on DIA-NN at library-build time --
probably unacceptable for a tool positioning itself as an alternative.
A run-learned residual correction is worth -21% sd, cross-validated by peptide,
and becomes useful only once ODIA scores its own data.

**Caveats that must travel with these numbers.** Single run, single gradient,
single library size -- every result is n=1 at the experiment level, and the
RT-swap null deserves a second run before being treated as settled. No verified
cross-tool RT accuracy comparison exists in a common unit on a common dataset;
do not let one into this document as fact. The tensor-level encoding is verified
correct: ODIA's outputs are bit-identical to OpenMS's reference values and to
AlphaPeptDeep's torch predictions.

---

## Fragment charge ranking — measured, and it is a MODEL difference (2026-08-04)

The suspicion was that raising the enumeration cap put doubly-charged fragments
into the candidate list without the ranking ever keeping them. Measured, per
precursor charge, on ODIA v3 against DIA-NN's library:

| precursor z | frag z | ODIA share | DIA-NN share | ODIA median intensity | DIA-NN median |
|---|---|---:|---:|---:|---:|
| 2 | 1 | 82.5% | 78.3% | 0.2547 | 0.2330 |
| 2 | **2** | **17.5%** | **21.7%** | **0.0983** | **0.1546** |
| 3 | 1 | 74.2% | 64.3% | 0.2617 | 0.2288 |
| 3 | **2** | **25.8%** | **35.7%** | 0.2386 | 0.2527 |

**Our selection is faithful to our own model.** For a 2+ precursor the PeptDeep
model puts a doubly-charged fragment at median 0.098 where DIA-NN's model puts
it at 0.155, so ranking by intensity demotes it exactly as it should given what
we believe. Both tools cap fragment charge at 2, so the candidate pools match
and the rule is not the difference -- the two intensity models are.

That reframes the fix. Forcing more z2 into the library **overrides the model
rather than correcting a defect**, which may still be worth doing if it finds
more peptides, and is worth nothing if it does not. So it is an option
(`-reserved_doubly_charged`, default 0) and an ablation rather than a change:
same build, same FASTA, same RT model, same search, only the quota differing.
Measured on a 200-protein library, a quota of 4 moves the z2 share from 21.6%
to 36.2%, which brackets DIA-NN's mix.

**Result: the quota costs identifications, so the default stays 0.** Full human
proteome, same build, same FASTA, same RT model, same search, only the quota
differing:

| reserved_doubly_charged | precursors | proteins |
|---|---:|---:|
| **0** (default) | **35,348** | **5,141** |
| 4 | 35,050 | 5,051 |

Forcing DIA-NN's charge mix loses 298 precursors (-0.84%) and 90 proteins
(-1.75%). So the model difference is real but overriding it is not an
improvement: our ranking is not merely faithful to our model, it is *better for
the search* than the mix it was being compared against. The option stays --
it is how this was measured and how it would be re-measured on other data --
but nothing should set it without repeating the ablation.

This closes the fragment-charge question.

---

## Fragment charge rules — the original suspicion, kept for the record

Raised by the RT research (2026-08-03), which ruled retention time out and left
this as the leading explanation for the remaining identification gap.

**The measurements.** With RT equalised between the two libraries, ODIA still
misses 4,937 of DIA-NN's precursors while finding 3,378 it does not. The missed
ones are **enriched in charge 3 -- 37.0% against 24.5% of the shared set** --
and have **median quantity 0.40x** the shared ones. The libraries are otherwise
structurally matched: 12.00 fragments per precursor against 11.92, b/y split
29/71 in both. But **ODIA emits 78.7% singly-charged fragments where DIA-NN
emits 72.4%**.

**Why this is not already fixed.** The enumeration rule was corrected earlier:
fragments may now carry up to the precursor's own charge, which was worth 44.7%
of the assays DIA-NN had and we lacked. What was never checked is whether the
*ranking* then selects them. Raising the cap only puts doubly-charged ions into
the candidate list; the top-12 cut is by predicted intensity, and if the model's
z2 channels are systematically weaker than its z1 channels, the cap change adds
candidates that never survive. The 78.7%/72.4% split says exactly that is
happening.

**Questions to answer, in order.**

1. Is the z2 deficit in the *model* or in our *use* of it? Compare the predicted
   intensity distribution of channels 1 and 3 (b_z2, y_z2) against channels 0
   and 2, for 3+ precursors specifically. The MS2 comparison already
   demonstrated that channel 1 is barely reached at charge 2 and needs a 4+
   precursor with two basic residues to appear at all.
2. Does DIA-NN cap fragments the same way? 11.92 per precursor is close enough
   to our 12 to suggest a similar cap, but the composition differs, so either
   the ranking or the candidate set does.
3. Would a charge-aware cap help -- for example reserving slots for z2 ions on
   3+ precursors rather than letting one global ranking decide?
4. Do the missed precursors' assays actually extract? Their 0.40x median
   quantity says they are low-abundance, so some of the gap may be sensitivity
   rather than assay choice. Phase 2 can answer this directly once it can
   extract at scale, and that is the cleanest test available.

**Do not assume the answer.** The last two confident diagnoses on this project
-- methionine excision for the coverage gap, retention time for the
identification gap -- were both wrong, and both were settled in minutes by an
ablation. Run the equivalent here before changing the ranking: rebuild the
library with a charge-aware cap and search it, rather than reasoning about it.

---

## MS2 prediction — what the review left open

- **`predictMS2` materialises every spectrum.** 100,000 peptides measured at
  94.8 s and +1717 MiB, the payload being 14.4 M floats spread over 100,000
  separate vectors. MS2 takes charge as an input, so the proteome case is ~4.0 M
  precursors rather than the 1.25 M distinct sequences RT needs: extrapolating,
  ~63 min and ~2.2 GiB of payload in ~4 M allocations, all live at once, on top
  of the generator's existing 1.23 GB. `predictRT` returns 4 bytes per peptide.
  There is no streaming or callback form, and this should be decided before the
  intensities are wired into `LibraryGenerator`.

- **`MAX_BATCH_ROWS` is not covered by any test.** Its purpose is to bound
  memory -- without it a 2 M-peptide run peaked at 26.7 GiB -- and removing the
  chunking entirely still passes, because values do not change. Correctness
  across the 2048 boundary was verified by hand (indices 2046-2050 bit-identical
  to solo runs) but nothing holds it.

- **CCS stays in square angstroms. Decided, not open** (user, 2026-08-03).
  `predictCCS` returns collision cross-section and must keep doing so. The
  conversion to 1/K0 is deliberately *not* ODIA's: it is instrument-specific and
  is applied downstream via the Mason-Schamp equation, where the drift gas and
  the instrument's calibration are known. ODIA must therefore never convert, and
  must not acquire a `-ion_mobility_unit` option that implies it could.
  Consequence to settle: the library's `im` column is 1/K0-shaped, so predicted
  CCS needs either its own column or an explicit unit tag -- writing CCS into a
  field consumers read as 1/K0 is precisely the silent unit error this decision
  avoids.

- **The model guard's reject path is untested.** All three shipped models are
  accepted, so nothing here exercises the case it exists for -- a model that is
  none of the three. Testing it needs a synthetic ONNX file and the `onnx`
  Python package is not in the environment. The per-call guards *are* tested,
  by feeding each `predict*` the other two models.

- **Two output-shape checks are equivalent mutants.** Removing the CCS rank
  check, or the batch-size check, changes nothing observable, because the real
  model always returns the right shape. They are cheap insurance against a
  re-exported model, not something any test here can hold.

- **A failed peptide is an empty spectrum**, which `Spectrum::at()` will index
  out of bounds. `predictRT` uses NaN for the same situation; there is no
  equivalent value here, so the contract is "check `positions` first" and it is
  only documented, not enforced.

---

## `.oswpq` reading — still open after the second review

- **`getInt64` on a `uint32`/`uint64` above the signed range is an equivalent
  mutant.** Narrowing the cast wraps the value, but ids are used only as join
  keys and the wrap is bijective, so precursors and transitions still meet.
  Only `UINT64_MAX` would alias the null sentinel. Nothing here can catch a
  narrowing cast, and nothing here needs to.

- **The two null-id guards mask each other.** A precursor with no id is kept out
  of the join map, and a transition with no id never looks one up; removing
  either alone changes nothing, and only removing both attaches an id-less
  fragment to an id-less precursor. Kept as belt and braces, with the pair
  covered by a test.

- **`ChunkedColumn` is not thread-safe** and now says so. The cursor is shared
  mutable state, so two threads scanning one column return each other's rows
  rather than colliding visibly. The header advertises the class as the hot
  path for a 64-thread scan, so this needs a per-thread cursor before any
  parallel reader is written.

- **A precursor whose m/z cannot be represented is kept at m/z 0** and can never
  match a window. `invalidMzCount()` reports it -- one on the real upstream
  bundle -- but `load()` does not surface it to the caller.

---

## `.oswpq` reading — what it does not yet do

- **The multi-chunk path is never exercised against a real file.** It is the
  reason `ChunkedColumn` exists, and it is unreachable below 2 GB of characters
  in one column: Arrow's Parquet reader concatenates row groups, so a fixture
  written with one row group per row still comes back as a single chunk. The
  cursor is therefore tested against chunk layouts built directly in
  `test/tools/odia_chunked_column_test.cpp`, which catches the defect but does
  not prove Arrow splits where we think it does. A proteome-scale bundle would
  settle it; none exists here.

- **Both tables are materialised whole before conversion.** For the benchmark
  library that is 78.6 M transitions of Arrow on top of the ODIA library being
  built. Column projection removes the largest contributor (`traml_id`), but the
  peak is still roughly double what it needs to be. Reading row group by row
  group would fix it and was not done, because nothing here can measure it: the
  largest bundle available is 18 transitions.

- **Fields read and then dropped**: `traml_id`, `unmodified_sequence`,
  `transition_id`, and the `detecting` / `identifying` / `quantifying` flags.
  The first is deliberate (D3: no per-transition strings) but it means ODIA
  cannot round-trip a bundle -- the human-facing precursor id is gone, and a
  written bundle would have to synthesise one. The three booleans matter for
  OpenSWATH scoring and will have to be carried before Phase 3.

- **The census check is advisory.** A disagreement between the row counts and
  `library/metadata.json` is recorded in `Stats` and nothing acts on it. It
  should probably be a hard failure by default, since the whole point is to
  catch a truncated read before an hour of extraction, but "probably" is not
  enough to make a load fail.

- **`schema_version` is checked; nothing else in the metadata is.** The
  `fragment_type_counts` and `charge_counts` blocks are a second, finer census
  that would catch a mis-parsed `type` column, and they are ignored.

- **The spec document disagrees with the file in three places** and should be
  corrected: it does not mention that `type`/`annotation` may be empty and
  `ordinal` `-1`; it does not mention that a precursor may have no transitions;
  and it says the transition `traml_id` is the precursor's, denormalised, where
  the sample bundle carries the transition's own id there instead.

---

## Smaller things noticed in passing

- Three decoys of 1,993,956 are skipped on the human proteome — peptides with no
  suitable unmodified position near both termini. Harmless, but the count should
  be reported rather than silently dropped.
- `LibraryGenerator` holds the whole peptide-to-protein map in memory before
  emitting. Fine at proteome scale (1.23 GB peak), worth revisiting if libraries
  get much larger.
- The mzpeak fork's test suite does not compile against Arrow 23
  (`parquet_writer_test.cpp` calls `FileReader::ReadTable()` with a signature
  this Arrow lacks), so `build_mzpeak.sh` builds only the library and tools.
- `ODIAInfo -peaks` is unusable on real files until the decode issue above is
  fixed; it is opt-in and documented as such.

---

## Upstream reports worth filing

- **OpenMS**: `WITH_ONNX=ON` cannot configure as shipped —
  `cmake/FindONNXRuntime.cmake` is not on the `CMAKE_MODULE_PATH` that OpenMS
  sets (`cmake/Modules` and `cmake/Windows` only).
- **OpenMS**: `WITH_ONNX=ON` compiles and exports the PeptDeep classes but never
  installs their headers. `libOpenMS.so` carries 18 PeptDeep symbols, and 26 ML
  headers are installed — every `ML/` subdirectory except `PEPTDEEP/` and
  `ONNX/`, the two the flag exists to enable. No external project can include
  them, so the feature is unusable outside the OpenMS tree.
- **OpenMS**: no mzPeak entry in `FileTypes`, so a TOPP tool cannot declare
  mzPeak as an input format.
- **OpenMS**: `ParquetFile::getColumn()` returns only the first chunk while
  callers loop to `num_rows` — a truncating read on any column past 2 GB.
- **mzpeak**: `EnumerableProxy::Iterator` defaults a constructor taking
  `const Iterator&&`, which is not a move constructor and cannot be defaulted;
  GCC rejects it. Carried as `patches/mzpeak-0001-defaulted-move-ctor.patch`.
- **mzpeak**: `install_headers(subdir: 'mzpeak')` flattens the header tree, so
  the installed headers do not compile.
- **mzpeak-convert**: `--layout chunked` silently produces a point-layout file
  while printing its own `BUG:` diagnostics about signal arrays spilling to
  `auxiliary_arrays`.

## GPU inference: built and verified, unrun for want of a credential (2026-08-04)

`opt/env-gpu` holds a CUDA 12.9 build of ONNX Runtime **1.26.0** -- the same
version the CPU path uses, so it is ABI-compatible with what OpenMS links --
plus the cudart/cublas/cufft/cudnn the provider dlopens, which the conda
package does *not* pull in on its own (it declares only the `cuda-version`
metapackage). It is a separate prefix: `opt/env` is what the read-only OpenMS
install resolves against and must not be swapped under.

`build-gpu/` is ODIA linked against it, on Ceph so the GPU nodes see it without
a rebuild. Verified on ibminode05: `GetAvailableProviders()` returns
`CUDAExecutionProvider, CPUExecutionProvider`, against `CPUExecutionProvider`
alone for the old build. So `prefer_gpu` finally has something to find.

`build-gpu/run_gpu_bench.sh` is staged and needs no arguments.

**What blocks the run is authentication, not routing.** On ibminode05:

- `ssh -v` reports every identity file as `type -1` -- none exist. `~/.ssh`
  symlinks to `/afs/wsi/home/oliver/.ssh` and holds `authorized_keys`,
  `config`, `known_hosts` and an unrelated GitHub deploy key.
- The two keys in `authorized_keys` have no private half on this filesystem,
  and the deploy key's public half does not match either of them.
- No Kerberos ticket (`klist`: no credentials cache), so `gssapi-*` is out too.

Every hop therefore refuses at the first one: `spock` and `data` directly,
`-J hive`, and `-J sshgw` all fail on publickey. Jump hosts do not help,
because the failure is not reachability.

Unblocking needs a private key on ibminode05 whose public half is in
`~/.ssh/authorized_keys` -- the setup step the `ibmi-hpc` skill documents.
Once that exists the benchmark is one command.

## GPU inference measured on an H100 (2026-08-04)

Full human proteome (20,416 proteins -> 2,127,559 target precursors), `data`,
one H100 PCIe, against the v5 CPU run on ibminode05:

| stage | CPU (v5) | GPU | |
|---|---:|---:|---:|
| retention time | 223.5 s | **17.1 s** | 13.1x |
| fragment intensities | 2,372.5 s | **67.3 s** | 35.2x |
| collision cross-sections | 368.1 s | **34.3 s** | 10.7x |
| **total wall** | **53:23** | **3:58** | **13.5x** |
| peak RSS | 2.04 GiB | 2.15 GiB | unchanged |
| CPU used | 5,527% | 105% | ~53x less |

**Read the 35x carefully.** It is against the *old* CPU path, which took ONNX
Runtime's intra-op default. Same node, same binary, same peptides, GPU against
the *fixed* CPU path: 54,500 peptides/s against 8,953 at 32 sessions -- **6.1x**.
That is the honest hardware comparison; the rest of the 35x was the threading
bug. Note `data` was at load 370 on 224 cores, so the CPU figure is depressed.

**The phase is no longer inference-dominated.** Inference is 119 s of a 238 s
run; digestion, decoy construction and writing 7.28 GB of TSV are the other
119 s. Further model speedups now buy at most 2x on this stage.

CPU and GPU are **not** bit-identical, and should not be expected to be --
different kernels, different reduction order. With the same MS2 model on both,
over 200,000 fragments: 98.05% agree within 1e-3, 99.99% within 1e-2, and 9
fragments (0.0045%) differ by more than 0.1 -- ranking ties near the top-12
boundary. CCS agrees to a median 0.019 A^2. The determinism test's bit-for-bit
guarantee is *within* a provider, which is what it claims.

RT differed by 3.2 min median between the two libraries, which is **not** a GPU
effect: v5 used the fine-tuned checkpoint
(`rtfinetune/integrated/peptdeep_rt_dynamic.onnx`) and the GPU run used the
stock OpenMS model. MS2 and CCS used the same stock models on both, which is
exactly why only RT moved.

### Two traps worth keeping

* **GPU 0 is unusable on both nodes.** `cudaSetDevice(0)` returns error 46
  ("devices busy or unavailable") on spock and data, with both GPUs in Default
  compute mode, zero memory used and no processes attached -- and world-writable
  `/dev/nvidia*`. On `data`, `cudaSetDevice(1)` succeeds. On `spock` **both**
  devices fail, so spock currently has no usable GPU at all. ONNX Runtime
  defaults to device 0, so **`CUDA_VISIBLE_DEVICES=1` is required** on data.
  Filed for the admin.
* **The conda `onnxruntime-cpp` cuda build declares only `cuda-version`**, a
  metapackage. cudart/cublas/cufft/cudnn must be installed explicitly or the
  provider fails to load and ONNX Runtime falls back to CPU silently.

## A comprehensive score inventory, and a learned score for SELECTION (2026-08-07)

Two connected pieces of work. We have 15 sub-scores; DIA-NN has **110**, and
OpenSWATH/pyProphet's published set is ~20. The 15 were chosen conservatively --
"a sub-score computed from a placeholder is worse than an absent one" -- and that
was right while extraction was broken. It is now the binding constraint: the true
peak is available for **97.4%** of precursors and we rank it first for **75.4%**,
so ~22 points sit in discrimination we are not computing.

### Part 1: the score inventory

Test every candidate below on **both** benchmark files (S08 diaPASEF, Astral),
against each file's DIA-NN confident set, measuring **best-ranked-right** -- not
on-RT, which counts TSV rows and is inflated by ~1.5 candidates per precursor
sharing a broadcast q-value.

**From DIA-NN (`diann.cpp:780` enum, handoff section 6.1). Cheap and applicable:**

| DIA-NN name | what it is | why it may matter here | cost |
|---|---|---|---|
| `pTightCorrOne/Two` | Σ fragment correlations at 0.45x and 0.20x tolerance | interference survives a loose tolerance and dies at a tight one; this is the single cheapest orthogonal signal we lack | re-match at 2 extra tolerances |
| `pAcc+0..5` | per-fragment `(\|obs-exp m/z\| / tol) x correlation` | the mass-error feature we keep deferring, but PER FRAGMENT and correlation-weighted | needs the ppm retention already planned |
| `pBestCorrDelta` | `pTimeCorr - best_corr_sum` for this precursor | we compute CANDIDATE_MARGIN, which is close; theirs normalises against the run | free |
| `pTotCorrSum` | `log(pTimeCorr / (total_corr_sum + 1))` | normalises a candidate against how much correlation the whole run offered | free |
| `pResCorr`, `pResCorrNorm` | correlation of fragments BEYOND the top 6 | we use all 12 equally; splitting top-6 from the rest is a real distinction | free |
| `pShape+0..4` | elution profile in 5 symmetric bins | peak shape as a vector rather than a width scalar | free |
| `pSig+0..5` | per-fragment share of integrated signal | we have YSERIES_SCORE and INTENSITY_SCORE; this is finer | free |
| `pCorr+0..11` | RAW per-fragment correlations, not summed | a summary hides which fragment disagreed | free |
| `pMinCorr` | correlation against 3-point-minimum traces | a baseline-insensitive variant | cheap |
| `pShadow`, `pShadowCorr` | correlation with -1.00335 Th shadow traces | an isotope shadow that a real peptide has and interference does not | needs a second extraction offset |
| `pHeavy` | correlation with +1 isotope traces from the NEIGHBOURING window | same idea across windows | needs cross-window extraction |
| `pdRT`, `pRT` | `sqrt(\|dRT\|/span)`, and position in gradient | we have RT_DELTA raw; the sqrt and the position are different shapes | free |
| `pMz`, `pCharge`, `pLength`, `pFrNum`, `pMods`, `pAAs+0..19` | peptide properties and AA counts | 26 features that need no signal at all. DIA-NN marks these `p_none` -- NN-only, never used by its linear classifier. Test them, but expect them to need the NN | free |

**Requires MS1 extraction, which we do not do at all (a whole sub-project):**
`pMs1TimeCorr`, `pMs1TightOne/Two`, `pMs1Iso*`, `pMs1Ratio` -- MS1 correlation and
isotope agreement. DIA-NN devotes ~8 features to it. Note `MS1PeakSelection` is on
by default there, i.e. MS1 participates in DETECTION, not only scoring.

**Not applicable:** the Q1 block (`pQLeft/pQRight/pQPos/pQNFCorr/pQCorr` x3) is
Scanning SWATH only.

**From OpenSWATH / pyProphet (mProphet lineage), for the ones we lack:**
- `xx_swath_prelim_score`, `bseries_score` (we have y-series only)
- `massdev_score` and `massdev_score_weighted`
- `isotope_correlation_score`, `isotope_overlap_score` (needs MS1)
- `norm_rt_score` (we have the raw delta)
- `elution_model_fit_score` -- fit an EMG to the peak and score the residual;
  genuinely orthogonal to every correlation we compute
- `sn_ratio` variants beyond our single LOG_SN
- library dot-product and manhattan variants beyond our two

### Part 2: a learned score for SELECTION, not just for FDR

**This is the structural gap, and it is worth more than any individual score.**
Today the discriminant runs ONCE, after candidates are chosen: the picker emits
candidates by correlation, sub-scores are computed, and the classifier separates
target from decoy. Nothing learned ever feeds back into WHICH candidate wins.

DIA-NN's `cscore = Sum_i w_i * score_i` (`Precursor::seek`, `diann.cpp:8121`) is
used to SELECT the winning peak group, with the weights re-fitted every iteration
and features gated by `min_iter_seek`. Selection and classification are the same
learned function.

We now have `PeakGroupScorer::refit` and a convergence loop, so the machinery
exists. The missing piece is using the fitted discriminant to re-select the best
candidate per precursor between rounds, not only to score the one already chosen.

Watch for: this closes a loop that can inflate its own confidence -- the
classifier would be selecting the rows it is then trained on. DIA-NN's answer is
the feature schedule plus `check_weights` sign clipping. Ours must be an external
check (best-ranked-right against DIA-NN) rather than an internal one, because
every internal number moves together.

### How to test without fooling ourselves

Leave-one-out over the final set, on both files, reporting best-ranked-right.
Two of our existing 15 are known to be worthless: USABLE_FRAGMENTS was a constant
12.000 for a week, and IM_DELTA is currently all-NaN because
`Options::observed_im` was never wired to MobilityCalibration's output. Any new
score has to beat that bar, and the ablation is what proves it.

## Memory, parallel occupancy and runtime across all phases (2026-08-07)

Every performance number this project has is a spot measurement of whichever
phase was under suspicion at the time. There has never been one profile that
covers the whole pipeline on both benchmark files, and the gaps have already
cost us: a 3.6x memory win sat in a hardcoded `BLOCK = 1024` for weeks, and a
1.02 GB index was per-transition when every value in it was per-precursor.

What to measure, per phase and end to end, on **both** files:

**Phases:** library load -> decoy generation -> (library generation if from FASTA,
which is ONNX-dominated and separately GPU-capable) -> mass calibration ->
mobility calibration -> pass-1 extract+score -> RT fit -> pass-2 extract+score ->
refine rounds -> output writing.

**Per phase:**
- wall time, and its share of the total;
- peak RSS attributable to the phase (RSS is monotone, so use the delta plus a
  heap profile where the delta is ambiguous -- `footprintBytes()` was undercounting
  the string-arena map by 3x until it was checked against an RSS delta);
- **parallel occupancy**: threads requested vs mean threads actually running.
  This is the number we have never measured and the one most likely to be
  embarrassing. Known: mzPeak decode holds a mutex across the decode itself
  (`src/util/parquet.cpp`, deliberately -- concurrent readers would corrupt each
  other through one seek-and-read handle), and decode is ~98% of a scoring run on
  the tiny library. Measured indirectly: 16 threads 18:46 vs 96 threads 9:51, a
  1.9x for a 6x thread increase. Sublinear, unexplained in detail.

**Scaling curves, not single points:**
- threads: 1, 4, 16, 32, 64, 96, 224 (both nodes have 224 cores);
- library size: the existing lib_tiny / n5k / n100k / lib50k / 4.26M ladder;
- file: S08 (32,210 spectra, 12.75 GiB, diaPASEF) vs Astral (303,701 spectra,
  3.10 GiB, no IM, no co-packed frames). Astral is 9.4x the spectra and took
  27:55 against S08's ~10 -- 2.8x, which is either good news about per-spectrum
  cost or bad news about something else, and nobody has looked.

**Specific open questions this would answer:**
- Where does the ~2.9 GiB decode floor actually come from now? It was 10.3 GiB,
  `-decode_block 256` took it to 2.9, and the composition was never re-derived.
- Does the co-elution picker's per-position pairwise correlation scale badly with
  spectrum count? It does far more arithmetic per cycle than a local maximum of a
  sum, and has only ever been measured on the small file.
- What is the actual thread ceiling, and is it the decode mutex, memory
  bandwidth, or the serial scoring on the extractor's thread (noted as unmeasured
  when the sliding window landed)?
- At 4.26M precursors, pass 1 wants ~274 GiB. Sampling it is blocked on the
  classifier's viability floor (below ~5,000 precursors it fails, silently,
  measured: 635/1820/1211/0 identifications from four equal-sized subsets). What
  is that floor exactly?

**Deliverable:** one table per file, phases x (wall, share, peak RSS, mean
occupancy), plus the three scaling curves. Publish it as doc/14 and keep it
current -- the reason this backlog entry exists is that every previous
performance claim was made from a measurement taken for a different purpose.

### Part 1b: compute the EXPENSIVE scores too, then let SHAP decide

Addendum. The table above sorts candidates by cost and implicitly suggests doing
the cheap ones first. That ordering is a trap: it selects features by
implementation convenience rather than by information, and we have no evidence
the cheap ones carry the signal. Compute them ALL, including:

- **MS1 block** (~8 features): `pMs1TimeCorr`, `pMs1TightOne/Two`, `pMs1Iso*` at
  three tolerances, `pMs1Ratio`. Needs MS1 extraction, which we do not do at all.
  This is the single largest missing capability -- DIA-NN also uses MS1 in
  DETECTION (`MS1PeakSelection`, default on), so building it may pay twice.
- **Shadow and isotope traces**: `pShadow`, `pShadowCorr+0..5` (-1.00335 Th),
  `pHeavy` (+1 isotope from the NEIGHBOURING isolation window). Each needs an
  extra extraction offset; `pHeavy` needs cross-window extraction, which our
  one-precursor-one-window assignment currently forbids.
- **Tightened-tolerance re-matching**: `pTightCorrOne/Two`, and the `p_fit` class
  DIA-NN only enables with `--tight-mass-acc-aux-for-cal`.
- **`elution_model_fit_score`** (OpenSWATH): fit an EMG per peak group and score
  the residual. Expensive per candidate and orthogonal to every correlation.

### Score selection by ML, with SHAP -- and why the obvious way is wrong

Fit the GBT over the full candidate set and rank features by **mean |SHAP|**,
then keep the ones that carry the model. Straightforward, and there are two
failure modes to design around:

**1. SHAP measures what separates targets from DECOYS, not what finds the right
peak.** Those are different objectives and we have already measured them coming
apart: off-RT targets are statistically identical to decoys on `library_corr`
(median -0.036 vs -0.032), yet the classifier still admitted them, because it
learned signal-PRESENCE features that separate a real-but-misplaced target from a
shuffled decoy. A SHAP ranking would have rewarded exactly those features. So the
selection metric must be the **external** one -- best-ranked-right against each
file's DIA-NN confident set -- with SHAP used to generate the candidate ordering,
not to make the decision.

**2. Importance is dataset-specific, and our two files differ structurally.**
S08 is diaPASEF (ion mobility, co-packed two-windows-per-frame); Astral has
neither. `IM_DELTA` is definitionally worthless on Astral. Anything derived from
frame packing is S08-only. So the procedure is:
  - fit and rank on S08 and on Astral **independently**;
  - report both rankings side by side, plus the rank correlation between them;
  - keep the union of what is decisive on either, not the intersection -- a
    feature that only works on one instrument class is still worth having, gated;
  - treat a large ranking disagreement as a finding about the instruments, not as
    noise to average away.

A third file would make this much stronger, since two points cannot distinguish
"instrument-specific" from "this particular run". `12_80` (SCIEX SWATH, no IM,
different vendor) is on disk and unused.

**Guard against the loop closing on itself.** If SHAP selection feeds a
discriminant that also SELECTS candidates (Part 2), the model chooses the rows it
is trained on and every internal metric agrees with itself. Hold out precursors,
not rows, and judge on the external metric only.

## Fine-tuning belongs inside the run's loop, not in the library on disk (2026-08-07)

**Measured regression.** Library v4 scored **37,583** confident precursors in
DIA-NN; v5 scored **35,556**. v4 BEAT DIA-NN's own library (37,247). The
generation recipes are byte-identical except one path:

    v4:  -rt_model .../rtfinetune/onnx/peptdeep_rt_dynamic.onnx
    v5:  -rt_model .../rtfinetune/integrated/peptdeep_rt_dynamic.onnx

The `integrated/` model's own `rt_provenance.json`:

    "tuned_from": "search_diann/report.parquet",
    "peptides": 500, "epochs": 40,
    "warning": "This model is specific to the run it was tuned on and must not
                be reused across runs or gradients."

500 peptides at 40 epochs, and the file says not to reuse it. The script passed
it anyway, and it cost 2,027 precursors -- 5.4%.

A guard now warns when `-rt_model` has a provenance file carrying a warning.
That is a plaster; the design below is the fix.

### The intended architecture

Library generation moves INSIDE the calibration loop, and the fine-tuned model
is scoped to one run and never written anywhere another run can find it:

1. generate the first library with the **stock** model;
2. search, calibrate RT (and mass, and mobility);
3. fine-tune the RT model **on this run's own confident identifications**;
4. re-predict the library with the tuned model;
5. re-calibrate;
6. iterate 3-5 until identifications stop improving;
7. discard the tuned model with the run.

The distinction that matters: a model tuned on run A and used to build a library
for run B imports A's gradient into B's predictions. Tuned on run A and used
within run A it is legitimate refinement -- the model has only seen data from
the run it is predicting for.

### What this needs

- Library generation callable mid-run, not only as `-stop_after library`. Today
  it is a separate invocation writing a TSV.
- Fine-tuning in-process, or as a subprocess the loop drives. Today it is the
  standalone `rtfinetune/` pipeline (`prep.py` + training, ~413 s for 500
  peptides / 40 epochs on one run).
- 500 peptides / 40 epochs is almost certainly overfitting. Inside a loop with
  the run's own identifications there are thousands available -- 1,018 anchors
  at q<=0.05 on a 2,665-precursor subset, so a full library gives far more.
  Tune the count and epochs against held-out identifications from the SAME run,
  which is the check the standalone pipeline never had.
- A convergence criterion that is not the identification count, which swings
  ~10% between refits on fixed input.

### The cross-run trap this closes

`odia_v5.tsv` carries iRT from a model tuned on S08's DIA-NN results. Every
recovery number measured on S08 with that library has partly seen the answers.
It is not the reason the library underperforms -- it underperforms because the
tuning was bad, not because it leaked -- but both are reasons the tuned model
must not be a persisted artefact.

### Confirmed 2026-08-07: v6 reproduces v4 and beats DIA-NN

Regenerated with the stock RT model, via `odia_v4_ft.sh` with only the output
name changed (sourcing env.sh, applying ccs_to_mobility, calling search.sh):

  odia_v6  stock RT model, CURRENT code   37,596
  odia_v4  stock RT model                 37,583   (+13, 0.03%)
  DIA-NN own                              37,247   (v6 +349, +0.9%)
  odia_v5  run-tuned RT model             35,556   (v6 +2,040)

Settles three things. No code regression -- everything changed this week left
library generation untouched. The v5 loss was ENTIRELY the run-tuned model,
all 2,040 recovered by the path swap. And our generator beats DIA-NN's own
library, reproducibly, twice.

Generation is also ~5x faster than when v4 was made: RT 215.7 -> 44.1 s, CCS
364.8 -> 89.7 s, from the data-parallel ONNX sessions (2a4c73a, 43af3fc).
Both stages use `inferenceSessions(sessions)` with 1 intra-op thread.

**Process note, because it cost four wrong numbers.** Getting here took four
attempts, each failing because I rebuilt the pipeline from memory rather than
running the script sitting beside the data:

  1. no --fasta            -> 20,065, reported as a "53.9% library gap" that
                              does not exist; it drove a 2x2 cross-analysis
                              and a plan step before the correctly-configured
                              search turned up already on disk
  2. no ccs_to_mobility    -> 14,408; nearly read as refuting the diagnosis
  3. hand-written DIA-NN   -> search.sh exists so the settings are identical
  4. no env.sh             -> CCS all-zero, announced as "a real code
                              regression in HEAD" on the strength of a grep
                              matching "CCS model available" inside the string
                              "No CCS model available"

Every one produced a plausible number that was reasoned from before being
checked. **When comparing against a prior result, run the prior result's
script** -- not an equivalent command. `odia_v4_ft.sh` is four lines and would
have been right the first time.

## The scorer has never been run against a realistic library (2026-08-07)

**Measured, S08, `v6_50k.tsv` (a stride sample of our own v6, no run-tuned iRT,
neither tool privileged):**

    DIA-NN   738 of 50,000 at 1% FDR   1:43
    ODIA       0                      36:43, 9.09 GiB

The true positive rate is **~1.5%** -- a whole-proteome library is almost
entirely peptides absent from the sample. DIA-NN gets 738/50,000 here and
37,596/2,127,559 (1.8%) on the full library, so the sample is representative.

ODIA's guards fired correctly: 1,639,188 peak groups, **819,353 target and
819,835 decoy (1:1)**, then `pass 1 identified nothing at 1% FDR`. The positive
class is ~98.5% noise, so the semi-supervised loop has no separable seed at
iteration 0 and never ignites.

**Which picker criterion rejects most, instrumented on that same run**
(`PickerRejects`, over 107,634,884 scan positions -- note criteria are counted
independently, so they sum to more than the positions):

    not a local maximum        189,557,902
    below min_corr_score        50,689,548
    below apex_evidence         36,856,546
    fewer than 2 fragments      30,191,757
    reference trace zero        11,525,195
    outside max_corr_diff            25,997
    then min_fragments_at_apex        75,479 candidates dropped
    -> 19,150 precursors yielded NO candidate at all, 0 identified at 1% FDR

The local-maximum test dominates by ~4x over anything else, and
`max_corr_diff` is doing essentially nothing (26k, a 7000x smaller effect than
the leader) -- so `-max_corr_diff` is a tuning knob with no leverage here and
should not be spent effort on. **The correlation thresholds together
(min_corr_score + apex_evidence = 87.5M) are the second force**, and those are
the ones that encode "is this a co-eluting peptide", so they are where a
realistic-library fix has to be careful: loosening them re-admits the 98.5%
noise the FDR then cannot separate.

Note this does NOT implicate the ion-mobility window: the mobility work of
2026-08-07 (see the 1/K0 section) found +/-0.025 already optimal, so the
picker's loss is not a window-width problem.

**This invalidates the regime, not the work.** Every recovery figure this week
-- 75.4% best-ranked-right, 89.1% precision, 97.4% availability -- used
`lib_targets.tsv`, 2,665 precursors drawn from DIA-NN's own confident set, i.e.
~100% true positives. Those measure ranking quality among true positives, which
is real. They say nothing about a library anyone would search.

It also kills the "viability floor below ~5,000 precursors" diagnosis: this run
had 50,000. The governing quantity is the true positive RATE.

~~**The fix is DIA-NN's batching**~~ -- **NO. MEASURED AND REFUTED
(2026-08-07).** This was my inference from handoff 7.1, not a measurement, and
it is wrong. Batching cannot bootstrap from a subset that itself yields nothing,
and that is exactly what a subset yields.

Simulated at the real regime -- 10,000 targets, 10,000 decoys, a planted true
positive rate, the actual `scoreSemiSupervisedLDA` with GBT. A RANDOM subset of
2,000 (same rate, a tenth of the rows) reports **zero in all eight
rate x separation conditions**, including the one where the full 20,000 reports
676. Fewer rows means fewer decoys, coarser q resolution and a weaker fit.
Subsetting is uniformly worse, never better.

(The first version of this probe drew a PREFIX rather than a random subset, and
the true positives are planted at the front, so the "subset" was a pre-selected
library and every subset number was inflated. That is the ranking-vs-search
confusion below, reproduced inside the very experiment meant to study it.)

**The constraint is discrimination power, and the lever is FEATURE COUNT.**
Splitting ranking from certification at a 1.5% rate and 3 sigma: the top 150
targets by the learned dscore are **48.7% genuinely true against a 1.5% base
rate** -- a 32x enrichment, so the ranking works -- while q there is 0.38
against an actual FDP of 51%, so the estimator is roughly HONEST too. Neither
is broken. At a 1.5% prior, 1% FDR needs a likelihood ratio near 6500:1, and
that is simply more evidence than a handful of features carries.

Feature count crosses it. Same regime, 3 sigma, features as independent noisy
looks at the truth (the optimistic case -- real sub-scores are correlated, so
real counts must be HIGHER than this):

    features   top-150 purity   reported at q<=0.01
       4           52.0%              0
       8           72.0%              0
      15           84.7%              0        <- what ODIA has
      30           87.3%            104
      60           95.3%            116
     110           94.7%            137        <- what DIA-NN has

The cliff is between 15 and 30, and ODIA carries exactly 15 sub-scores. That is
a direct explanation for 0 against DIA-NN's 738 on the same file, and it makes
the score-inventory item below THE priority rather than a nice-to-have. At 2
sigma even 110 features report nothing, so per-feature strength matters as much
as count -- more good features, not merely more.

Probes: `scratchpad/batch/{probe,diag,feat}.cpp`, header-only against
`odia/scoring/lda.h`, seconds to run.

`test/tools/odia_entrapment.cpp` does NOT reproduce this regime: it plants a
33% true positive rate, where the real library is 1.5%. It is a valid FDR
calibration test and a poor model of the collapse. A low-rate case belongs in
it.

**Benchmark hygiene:** always state what fraction of a library's targets are
actually present. A library built from the comparator's confident set is a
RANKING benchmark; a proteome-scale library is a SEARCH benchmark. We have only
ever run the first.

## MS1 is present and unused -- the largest untapped feature family (2026-08-07)

`ODIAInfo` on S08: **1,343 MS1 spectra against 16,105 MS2** (17,448 total).
Over the gradient that is a ~1.8 s duty cycle, so a 30 s peak is sampled ~15
times -- ample for a chromatographic trace, not merely a survey.

**ODIA reads none of it.** `src/io/MzPeakSource.cpp:58` drops `ms_level != 2`
at INDEX time, so MS1 is invisible everywhere downstream; nothing in ODIA's
output reveals the file even has it.

This is the concrete instance of the feature-count finding above. ODIA carries
15 sub-scores and the simulated cliff sits between 15 and 30; MS1 is the
biggest single family available, and it is orthogonal to everything present
today (every current sub-score is computed on MS2 fragment traces, so they
share their failure modes -- a co-eluting interferent corrupts all of them at
once). DIA-NN's 110 include a substantial MS1 block.

**Plan, in dependency order:**

1. `MzPeakSource`: keep MS1 spectra in a SEPARATE index rather than merging
   them into `info_`. Merging would break every consumer's assumption that an
   entry has an isolation window. The mobility handling is the same
   diaPASEF-derived-boundary problem already solved for MS2.
2. An MS1 extraction pass over the precursor monoisotopic m/z plus the first
   two isotopes, on the same RT and 1/K0 windows the MS2 extraction uses, so
   the traces are directly comparable cycle for cycle. Reuse the existing
   `LiveSlot`/`BlockPool` shape; the row count is ~3 per precursor against ~12
   fragments, so the memory term is small.
3. Features, all cheap once the traces exist:
   - MS1/MS2 co-elution: correlation of the precursor trace against the
     fragment consensus. A real peptide's precursor and fragments share one
     elution profile; an interferent's do not.
   - Isotope ratio agreement: observed vs theoretical (averagine) for M, M+1,
     M+2. This is strong orthogonal evidence and needs no new extraction.
   - MS1 apex RT delta against the MS2 apex.
   - MS1 log signal-to-noise.
4. Measure on `v6_50k` -- the SEARCH benchmark, where the current answer is 0.
   `lib_targets` cannot show this: it is ~100% true positives, so it measures
   ranking among true positives and the whole point here is discrimination at
   a 1.5% prior.

**Do not read a gain on `lib_targets` as progress on the real problem.** That
is the ranking-vs-search confusion this document keeps having to restate.

## The mass calibration gate is measured once, at the worst moment (2026-08-07)

Adversarial review prompted by BOTH benchmark runs failing the gate.

**1. It cannot change between rounds, by construction.** `mass_model_known_`
(`src/OpenDIAlyzer.cpp:1823,1841`) latches after the first call. Pass 2 skips
`calibrate()` and re-prints the cached `Model`. The peakedness numbers come off
`Model`, so they are byte-identical across rounds -- Astral reported
"6.90 vs 12.00" twice. The tell is the diagnostics parenthetical: round 1 prints
"(from 3000 target and 5769 control cells over 24000 spectra, 79.73 s)" and
round 2 prints nothing, because pass 2 hands `report()` a fresh empty
`Diagnostics`. So the answer to "how do the residuals change between rounds" is
**they do not, and cannot**.

**2. The one measurement happens before the retention-time map exists.**
`applyMassCalibration_` is called at `OpenDIAlyzer.cpp:432`, the map is fitted
at `:772`. The run's own log says it: "the library is being spread evenly over
the run, which will extract from approximately the wrong retention times." The
mass probe therefore looks in the wrong place and mostly matches noise -- which
is precisely what the S08 gate then reports as "residuals are FLAT ... that is
what a mostly-noise sample looks like".

`MassCalibration.h:95-97` documents the assumption: "`rt_trafo`. The reference
locates its anchors in time through a fitted retention-time map. **There is none
before the first pass**, which is what the apex selection above replaces." That
was true when written and is now stale -- there IS a second pass, and
`applyMobilityCalibration_` already defers to it ("DEFERRED -- it is measured at
the peak groups this run scores"). The mobility arm was upgraded; the mass arm
was not. `MassCalibration::Options` still has no RT-map field.

**3. On Astral the gate fails on a statistic that carries no information.**
`peakednessRatio` is count(central 0.2*window) / count(edge band), so its
relative error is ~sqrt(1/c + 1/e). Astral produced **4,033 target residuals and
98 control**:

    target   6.90  over 4033   -> ~ +/- 0.74
    control 12.00  over   98   -> ~ +/- 7 to 12  (e is 1-3 counts)

The gate rule `decoy_peakedness >= peakedness` then compares 6.90 against a
number whose error bar is larger than itself. That is a coin flip, not a test.
5,769 control cells yielded 98 residuals because a 7 Th-shifted query rarely
matches anything -- the null is undersampled 41x by construction.

S08 is NOT this failure: 7,452 target and 2,645 control, both well sampled, and
2.22 against 3.00 is a real flatness. On a 1.5% true-positive library that is
the correct verdict -- but it is confounded with (2), so we do not yet know how
much of the flatness is the library and how much is the wrong-RT probe.

**Fix, in order:**
1. Un-latch: re-measure in pass 2 once the map exists (mirror the mobility arm).
2. Give `MassCalibration::Options` the iRT map and gate the probe's cells on it,
   as `im_calib_rt_window` already does.
3. Require a minimum control count before the `decoy_peakedness` rule may fire
   at all; below it, fall back to the absolute `min_peakedness` test only.
4. Then re-measure on both files. Expect target peakedness to RISE on both.

## Both benchmarks, current HEAD (2026-08-07)

    file    library           TP rate  type      at 1% FDR  wall    peak RSS
    S08     v6_50k (ours)      1.5%    SEARCH        0      42:31   9,331 MB
    Astral  astral_lib_own    ~100%    RANKING   4,275      35:17   3,148 MB

**These are not comparable and must never be quoted side by side.**
`astral_lib_own` holds 10,891 precursors against a DIA-NN confident set of
11,112 -- 98% of the truth is in the library, so it is ~100% true positives and
measures RANKING. `v6_50k` is 738 of 50,000, i.e. 1.5%, and measures SEARCH.

**Astral: the bottleneck is EXTRACTION, not scoring.**
  * 10,891 of the truth are in our library (98.0%) -- library coverage is fine.
  * ODIA produced a candidate for **5,704 (52.4%)**. **5,187 precursors the
    library carries got no candidate at all** -- genuine extraction loss.
  * Of the 5,704 reachable, 4,275 (75%) are identified at 1% FDR, and **100% of
    ODIA's 4,275 are in DIA-NN's set** -- zero false agreement. The scorer is
    not the problem here; the extractor is.

**S08: the bottleneck is discrimination.**
  * Extraction reached 576 of 738 (78.0%) -- fine.
  * Their best q: median 0.783, **0 at <=0.01 but 74 at <=0.05**. Not a total
    collapse; the discriminant does rank some correctly and cannot certify them,
    exactly as the feature-count finding predicts.
  * Pass 1 identified nothing, so no RT map was fitted at all and pass 2 ran
    without one -- the same wrong-RT condition as the mass gate above.

## The mass calibration is fixed as a MEASUREMENT and broken as an ACTION (2026-08-07)

Following the gate review above, both defects it named were fixed and both fixes
LOST identifications. Astral, our own library, at 1% FDR:

    latched, gate fails both rounds (HEAD)      pass1 1248   pass2 4275
    + min_control_residuals = 400               pass1  119   pass2 1994
    + re-measure against the fitted RT map      pass1 1248   pass2 1895

**Every variant that makes the gate PASS loses badly.** That is the finding.

**The re-measurement works, by its own metrics.** Probing at the fitted
retention times took the control from 98 residuals to 149 and it stopped
out-peaking the data (the pass-1 gate's whole complaint), and the fitted offset
moved -1.27 -> -1.75 ppm. So the "12.00 vs 6.90" null that failed pass 1 was
itself partly an artefact of probing the whole gradient. The diagnosis was
right.

**And applying it costs 56% of the run.** A passing gate does two things: it
centres the window AND narrows `fragment_ppm` off the uncalibrated 15. On Astral
the narrowing is the more expensive error EVEN WHEN CORRECTLY CENTRED -- the
opposite of what `ChromatogramExtractor.h` assumes when it says "an uncentred
narrow window is the worse of the two errors". That assumption was measured on
S08 and does not transfer.

**The thin-control guard was wrong, and instructively so.** `peakednessRatio` on
98 residuals genuinely is 12.00 +/- ~9, so the statistic is uninformative --
but the SPARSITY is not. Almost nothing matching a 7 Th-shifted query is itself
evidence that the fragment matches are not clean enough to calibrate from, and
the gate was reading that correctly through a bad estimator. Ignoring a noisy
number because it is noisy threw away the signal carried by why it was noisy.
Defaulted to 0.

**Next, in order:**
1. **Separate the width from the offset.** Apply the fitted offset while KEEPING
   the wide window, and measure. This is one run and it decides whether the
   offset is good and only the width is wrong, or the whole correction is.
2. If the offset is good: make the calibrated width a function of the measured
   `sigma_after` with a floor, rather than a fixed narrow value.
3. Only then consider defaulting `-mass_calibration_remeasure` on.

Both fixes are IN the tree and OFF: `-mass_calibration_remeasure` (flag) and
`Options::min_control_residuals` (0). Neither is reverted, because the
measurement half is correct and will be wanted once the application is fixed.

## Three of my own claims, corrected by adversarial review (2026-08-07, late)

### 1. "Feature count is the lever" -- REFUTED. It was an artefact of my own generator.

`scratchpad/batch/feat.cpp:34` drew every feature as `shift + noise*2.0`:
**independent given the class**. The sufficient statistic is then the row mean
and d' grows as sqrt(n_features) BY CONSTRUCTION. The two earlier probes,
`probe.cpp:32` and `diag.cpp:31`, drew a SHARED latent
(`common*0.7 + noise*0.6`) -- so the 4-feature row of `diag.cpp` and the
4-feature row of `feat.cpp` were never the same experiment, and I did not
notice I had changed the generator.

Re-run with the shared latent, everything else identical (same 1.5% rate, same
3 sigma, same seed, same `scoreSemiSupervisedLDA` + GBT):

    features    4    8   15   30   60  110
    at q<=0.01  0    0    0    0    0    0

**Zero at every count.** The "cliff between 15 and 30" exists only under
independence. And ODIA is squarely in the correlated case -- by my own words
two entries above, all 15 sub-scores read MS2 fragment traces and "share their
failure modes", which is exactly the regime that reports nothing at 110.

The corrected lever is **ORTHOGONALITY, not count**. That still points at MS1,
but for a different and stronger reason: a second measurement with an
independent failure mode, not a bigger number of columns. It also says that
working through the `doc/13` inventory of additional MS2-derived sub-scores
buys approximately nothing, which the previous framing did not.

Two further defects in that simulation, both real: it plants ONE row per group,
so the best-of-N maximum over ~9 candidates per precursor -- which is what
`lda.h` actually computes FDR on -- cannot appear; and the reported count is a
threshold crossing rather than a dose-response (seed 0xE47A9 gives
30/60/110 -> 104/116/137, seed 1 gives 115/133/134, the 60-vs-110 order flips,
while `q at 150` moves only 0.166 -> 0.147 across the "cliff").

### 2. The mobility slope is REAL but INFLATED, and the cross-charge agreement was an artefact.

`delta = im_observed - im_library`, and the fit regresses `delta` on
`im_library` -- i.e. regresses (y - x) on x. With prediction noise in the
library value that is **negative by construction**, slope
= -Var(noise)/Var(im_library), no scale error required. Refitting on
`im_observed` instead separates them: a genuine scale error stays negative on
both axes, pure dilution flips positive.

    charge    n      on im_library    on im_observed    dilution predicts
      2    16498        -0.0974          -0.0703             -0.0382
      3     5986        -0.1105          -0.0243             -0.0931

**Both stay negative, so a real scale error exists** -- the claim is not
refuted. But the magnitude fitted on `im_library` is inflated by dilution, and
charge 3 is very nearly ALL dilution (-0.1105 observed against a -0.0931
dilution prediction, leaving -0.0243 on the observed axis).

So "charge 2 and charge 3 agree at -0.096 and -0.113, two independent fits" was
wrong twice over: they are not independent (one CCS model, one coefficient, one
run -- and if the shared-coefficient hypothesis is true, agreement is
guaranteed rather than evidence), and on the honest axis they do NOT agree
(-0.070 against -0.024).

**This predicts the anomalies already recorded.** Fitting on `im_library`
over-corrects, i.e. shrinks toward the population mean, which always improves
MSE while moving the window off precursors whose library value was already
right. That is exactly "18.8% of out-of-fold MSE removed but -20
identifications", and "66 precursors LOST against 116 gained". The
"MSE is not a proxy for identifications" finding has a mechanism now.

**Fix:** an attenuation-corrected estimator (Deming / geometric-mean
regression), since at apply time only `im_library` is available. Roughly, the
true slope is near the midpoint of the two regressions: about -0.084 for
charge 2 and -0.067 for charge 3.

### 3. Astral's bottleneck is the PICKER, not extraction -- and the number was already in the log.

I reported "5,187 precursors the library carries got no candidate at all --
genuine extraction loss". `assess.py` labels `ref_ids & all_target` as
"extraction produced a candidate", but that label is an assumption, and the
measurement that separates the stages was printed by the same run:

    picker rejections over 12,912,142 scan positions:
      3,031,555 <2 fragments present, 6,870,896 below min_corr_score,
      3,702,166 reference trace zero, 5,725,208 not a local maximum
    -> 9,890 precursors yielded no candidate peak group

The chromatograms were extracted; the PICKER returned nothing for them. Same
shape as S08's 19,150. Other candidates not yet excluded: precursors dropped by
`min_fragments_at_apex` AFTER a candidate was emitted (a scorer filter, not a
picker one), pass-2 extraction on an RT map fitted from only 1,248 anchors, and
`Precursor.Id` string mismatch -- `writeScores_` concatenates modified sequence
and charge with NO separator, and a systematically mismatched class would land
in neither `hit` nor `avail` while leaving the "100% of ODIA's are in DIA-NN's"
figure untouched.

**Split the 9,890 four ways before spending anything on the extractor.**

## The pipeline IS deterministic, including thread-invariant (2026-08-07, late)

D1 called determinism "a precondition, not a deliverable" and it had never been
demonstrated. It is now, and the answer is good.

Three runs, S08 + `lib_targets`, `-ion_mobility_calibration anchors`, same node,
differing only in thread count:

    det_a   64 threads   1046 at 1% FDR   md5 e0175e928f30a37f7bd9230011efc4a5
    det_b   64 threads   1046             md5 e0175e928f30a37f7bd9230011efc4a5
    det_c   16 threads   1046             md5 e0175e928f30a37f7bd9230011efc4a5

**Bit-identical output files, 70,836 rows each, across a 4x thread change.**

This matters for how every A/B in this document may be read:

* **n=1 per arm is sufficient at fixed configuration.** The adversarial review
  argued the +67 and -20 deltas sat inside a "~10% swing between refits" and
  were therefore unreadable. For fixed configuration that is refuted: the swing
  is exactly zero. The `635/1820/1211/0` figures elsewhere came from four
  DIFFERENT subsets of a library, i.e. different data, not repeats.
* **The `stable_sort` on `(precursor, apex_rt, apex_intensity)` is doing its
  job.** The concern that a thread-count change reorders rows, which reassigns
  folds, which changes the fit, is real in principle and does not occur here.
* **A difference between two runs is therefore always a real difference** --
  which also means every config difference must be controlled explicitly. The
  1046 here against 1232 recorded earlier for "the same" config is NOT
  nondeterminism: `run.sh` adds `-stop_after score -max_candidates 25` and the
  binary now clamps the composite mobility correction. Two changed inputs, two
  different answers, both reproducible.

**Still not covered:** permutation invariance (reordering the library), and
determinism of the library-generation and ONNX paths. This measures the
extract-score pipeline only.

## C1 done: Astral's 5,187 is real, and it is the PICKER (2026-08-07, late)

The four candidate explanations, settled against the run that already existed.

**ID join: clean.** `writeScores_` concatenates modified sequence and charge
with no separator, and the suspicion was that this mismatches DIA-NN. It does
not -- all three agree on `AAAEVAGQFVIK2`. **Zero** ODIA output rows carry an id
absent from our library, which is the direct test for a format bug.

**Library coverage: 10,891 of 11,112 (98.0%),** and the missing 221 are almost
exactly the 220 MODIFIED precursors in DIA-NN's set. **Our Astral library
contains zero modified peptides.** That is a library-generation gap, not a
search failure, and it is the whole of the coverage shortfall. Worth its own
item: the generator emits no modifications at all for this run.

**The 5,187 are all in-library and all unmodified**, so they are a real pipeline
loss, and the picker census from the same run locates it:

    picker rejections over 12,912,142 scan positions
    -> 9,890 precursors yielded no candidate peak group

5,187 is a subset of that 9,890 (the rest are precursors absent from DIA-NN's
set). **The chromatograms were extracted and the picker returned nothing.**

Still unseparated within that: precursors where a candidate WAS emitted and then
dropped by `min_fragments_at_apex` (a scorer filter), and the contribution of
pass-2 extraction running on an RT map fitted from only 1,248 anchors. Both need
a counter rather than a grep.

**So "Astral's bottleneck is extraction" is corrected to "the picker's
thresholds".** Same conclusion as S08's 19,150, and the same three criteria
dominate: `not a local maximum`, `min_corr_score`, `apex_evidence`.

## Mass calibration, settled with correct code, and CLOSED (2026-08-08)

Three matched arms on Astral + `astral_lib_own`, corrected code (the
out-of-bounds RT gate is fixed), all sharing pass 1 = 1,248 so everything
upstream is identical:

    baseline, gate never passes                  4,275
    offset + m/z shape, window NOT narrowed      3,895    -380   (-8.9%)
    offset + shape + narrowing                   1,434  -2,461  (-63%)

**The window width is the culprit, 6.5x the offset's cost.** The original
hypothesis was right; it just could not be supported at the time because the
measurement ran on undefined behaviour. The VOID figures (98->149 residuals,
4,275->1,895) are superseded by these.

Two things this does NOT say:

* `-mass_calibration_offset_only` overrides only `fragment_ppm`; the m/z slopes
  are already applied by then. So "-380" is the JOINT cost of the constant
  offset and the m/z shape, not the offset alone. Separating them needs a third
  mode that zeroes both slopes.
* It says nothing about S08 or about a search benchmark. Astral +
  `astral_lib_own` is ~100% true positives, so this is a ranking measurement.

**What to do with it:** `ChromatogramExtractor.h` says "an uncentred narrow
window is the worse of the two errors", measured on S08. On Astral the ordering
is reversed and by a wide margin. The rule is instrument-specific and must stop
being applied as if it were general. The honest fix is to size the width from
the fit's own `sigma_after` with a floor, rather than from a constant -- but
that is a tuning change on a subsystem that cannot move the headline number.

**CLOSED per doc/12 section D.** Five commits deep on a correction that remains
off by default. Both flags stay in the tree
(`-mass_calibration_remeasure`, `-mass_calibration_offset_only`); neither is
default. Revisit only after the prefilter gives pass 1 something to fit an RT
map from, which is the precondition the whole subsystem was found to lack.

## MS1 measured before building it, and the night's synthesis (2026-08-08)

`MzPeakSource` now keeps a SEPARATE MS1 index (`ms1Spectra()`, `ms1Peaks()`),
which is item 1 of the MS1 plan and the precondition for everything else. S08:
1,343 MS1 spectra reachable. 61/61 tests still pass.

`test/tools/odia_ms1_probe.cpp` then measured the premise before the subsystem,
S08 + `v6_50k` (SEARCH, 1.5% true positives, base 13.4/1000), mobility-gated at
+/-0.05:

    MS1 isotope depth (M, M+1, M+2 together)
      iso   targets   true/1000   enrichment
        3    49,916        13.4        1.0x     <- 99.8% saturated
       <3        84         0.0        0.0x

    MS1 monoisotopic intensity
      true median 4,844   absent median 3,056   ratio 1.6x
      top decile by intensity: 22.6/1000            1.7x

**Isotope presence saturates exactly as fragment depth did.** Intensity is the
first statistic all night to beat the base rate -- 1.7x -- and it is weak, and
it is probably confounded: DIA-NN identifies ABUNDANT peptides preferentially,
and MS1 intensity measures abundance, so this may re-measure abundance rather
than correctness. It is not evidence that MS1 discriminates a correct
identification from a wrong one.

### The synthesis: presence saturates, shape discriminates

Five statistics were measured tonight on the search benchmark:

    fragment depth (whole frame)              1.0x  saturated
    fragment depth (mobility-sliced)          1.0x  saturated
    qualifying_spectra                        0.5x  worse than random
    total_matches                             0.5x  worse than random
    MS1 isotope depth                         1.0x  saturated
    MS1 monoisotopic intensity                1.7x  weak, likely confounded

**Every "does this exist somewhere in the run" statistic is worthless here, and
the reason is structural.** Each is a maximum or a count over ~10^3-10^4
spectra, so it is an extreme-value statistic over thousands of draws: with 15
ppm tolerance and mobility-merged frames, coincidence is near-certain and the
answer is yes for everything. Increasing selectivity per draw does not fix it;
the draw count does.

This explains, retrospectively, the one large win this project has had. Replacing
the amplitude peak-picker with CO-ELUTION detection fixed RT, FDR, recovery and
memory at once. Co-elution is not a presence question -- it asks whether the
fragments rise and fall TOGETHER, which is a shape over time and cannot be
satisfied by an accumulation of unrelated coincidences.

**So the MS1 feature worth building is the one this probe did NOT test:
correlation of the MS1 precursor trace against the MS2 fragment consensus.**
That is a shape comparison, it is genuinely orthogonal (a different measurement
with an independent failure mode), and it is the only MS1 quantity with a reason
to work. `ms1_iso` and `ms1_max` should not be added as sub-scores -- one is
saturated and the other is an abundance proxy, and doc/07 is explicit that a
sub-score computed from a placeholder is worse than an absent one.

**Revised next step:** extract MS1 traces on the same RT/mobility windows as the
MS2 extraction, so the two are comparable cycle for cycle, and measure
MS1/MS2 co-elution correlation as a discriminator on `v6_50k` BEFORE wiring it
into the scorer. Same gate as everything else tonight: one measurement first.

## MS1/MS2 co-elution DISCRIMINATES: the first positive result (2026-08-08)

`test/tools/odia_ms1_coelution.cpp` builds the MS1 monoisotopic trace and the
MS2 top-6 fragment trace on one common time grid (the 1,343 MS1 acquisition
times; MS2 frames accumulate into their nearest bin) and correlates them.
Pearson over bins where EITHER trace is non-zero -- correlating over the whole
gradient would be dominated by jointly-empty bins, which agree perfectly and
mean nothing, and would reproduce the very saturation this probe exists to
escape.

S08 + `v6_50k` (SEARCH, 670 of 50,000 true, base 13.4/1000):

    corr >=     kept    true   true/1000   enrichment   recall
       0.8        76      14       184.2        13.7x     2.1%
       0.6       258      28       108.5         8.1x     4.2%
       0.4     2,042      62        30.4         2.3x     9.3%
       0.2    11,918     185        15.5         1.2x    27.6%
       0.0    42,925     605        14.1         1.1x    90.3%

**A monotone gradient with a top bin 13.7x the base rate** -- which is exactly
the success criterion `doc/08` set out ("a monotone gradient with a top bin far
above the base rate"), reached by a SHAPE statistic after every PRESENCE
statistic failed. Against the night's other six measurements (0.5x to 1.7x, all
saturated or worse than random), this is the first thing that discriminates.

### What it is NOT

**It is not a prefilter.** Recall at the discriminating thresholds is tiny: corr
>= 0.8 captures 14 of 670 true positives (2.1%), corr >= 0.6 captures 4.2%.
Used as a gate it would discard 90%+ of what we are trying to find. `doc/08`'s
filter role stays dead.

**The bulk separation is modest.** Medians are +0.115 true against +0.102
absent, and the standardised mean difference is 0.296 sd. The signal lives in
the upper TAIL, not in the centre of the distribution. A classifier will get
real but bounded value from it -- this is one informative feature, not a
solution to a 6,500:1 likelihood-ratio requirement.

### What it is

**A discriminating sub-score, and the first genuinely orthogonal one.** Every
existing sub-score reads MS2 fragment traces, so a co-eluting interferent
corrupts all fifteen together; this asks whether the PRECURSOR rises and falls
with them, which is a different measurement with an independent failure mode.
That is the orthogonality argument stated properly, and it is now measured
rather than asserted.

It also confirms the night's structural finding: presence saturates because it
is an extreme-value statistic over thousands of draws; shape does not, because
coincidences would have to arrive in the right ORDER.

### Next

1. Extract MS1 traces inside the extractor on the SAME RT/mobility windows as
   MS2, so the two are comparable cycle for cycle rather than through a
   nearest-bin approximation.
2. Add `MS1_COELUTION` as a sub-score. Do NOT add `ms1_iso` or `ms1_max` --
   saturated and an abundance proxy respectively.
3. Re-measure on `v6_50k` end to end. The honest expectation is a real but
   modest gain, not a fix for 0-versus-738.
4. Re-run this probe on Astral with a SEARCH library before generalising -- the
   instrument-conditionality trap has now caught this project twice.

### Co-elution is not an abundance artefact -- checked before building on it

Both traces are intensity-based, abundant precursors give cleaner traces, and
DIA-NN preferentially identifies abundant peptides. So the 13.7x could have been
abundance wearing a shape costume. Stratifying by MS1 monoisotopic intensity
into quintiles and re-measuring the top-decile-by-correlation enrichment WITHIN
each:

    quintile        n   base/1000   top decile by corr   enrichment
       1 (dim)  10,000        5.2                  3.0        0.6x
       2        10,000       10.7                 14.0        1.3x
       3        10,000       13.4                 22.0        1.6x
       4        10,000       18.2                 29.0        1.6x
       5 (bright) 10,000     19.5                 36.0        1.8x

    mean within-stratum 1.4x, against 1.5x unstratified

**It survives.** Controlling for intensity barely moves it, so the correlation
is carrying shape information and not brightness. The reverse check also holds:
intensity within a correlation stratum gives 1.8x against 1.7x unstratified, so
the two are largely independent.

**Error bar on the headline.** The 13.7x is 14 true of 76 precursors at
corr >= 0.8. Poisson error on 14 is about +/-3.7, so the enrichment is roughly
13.7x +/- 3.6x -- comfortably above 1, and imprecise. The bulk effect (top
decile, 1.4-1.5x) is the number to plan against; the far tail is real but thin.

**A caveat about the benchmark itself, which applies to every feature.** The
target is DIA-NN's confident set, so any feature that predicts "what DIA-NN
finds" scores well -- and abundance predicts exactly that, for reasons that have
nothing to do with whether an identification is correct. `ms1_max` is therefore
a hazardous feature to add even though it measures 1.7x: it would improve
agreement with the comparator without necessarily improving truth. Co-elution
does not have this problem, because a precursor co-eluting with its own
fragments is evidence about the identification itself.

**Decision: add MS1_COELUTION only.** Not `ms1_max` (abundance proxy, and the
benchmark rewards that spuriously), not `ms1_iso` (saturated at 99.8%).

## Picking or extraction? Both, and they are separable (2026-08-08)

`test/tools/odia_mz_residuals.cpp` probes at +/-50 ppm and records the
intensity-weighted signed m/z deviation per precursor, optionally restricted to
a KNOWN apex. Run on S08 + `v6_50k`, anchored to DIA-NN's own retention times
for its 670 reachable confident precursors (+/-30 s):

    group                  n     median ppm    sd      |ppm| > 15
    all true             670        -8.44    13.13        21.9%
    true & PICKED        576        -8.21    13.27        21.9%
    true & NOT PICKED     94        -9.38    11.94        22.3%

**Two separate answers.**

**1. The picking/not-picking split is NOT caused by m/z mis-calibration.** Picked
and unpicked true precursors have the same residual distribution -- medians
differ by 1.2 ppm against a 13 ppm spread, and the fraction outside the
extraction window is identical (21.9% vs 22.3%). Whatever decides that a
precursor yields no candidate, it is not that its fragments were extracted from
the wrong place in m/z. **The loss is genuinely in the picker.**

**2. Extraction IS mis-centred, by about -8 ppm, and it costs BOTH groups.** The
window is +/-15 ppm centred on ZERO because the mass calibration gate fails on
S08; the signal sits at -8.4. So the window effectively covers -6.6 to +23.4 ppm
around the truth, and **21.9% of true precursors fall outside it even at their
own apex**. That is a real, uniform loss, and it is exactly what the failing
gate was supposed to prevent. It does not explain the split, but it is a
first-order defect in its own right.

**Method note, and it is the same trap as everything else.** The first run of
this probe took the best co-occurrence over ALL 32,210 spectra with no RT
restriction, and produced medians of -5.8 (picked) against -9.6 (unpicked) with
sd 21.5 -- numbers that look like an answer and are not one. On a 1.5%-true
library the unrestricted maximum is a chance event, so its residuals are noise;
the tell was that DIA-NN's own true set showed sd 20.6 rather than the tight
cluster a real identification must have. Anchoring to the known apex is what
makes the measurement mean anything. Both panels are plotted in
`scratchpad/mz_residuals.png`.

**Caveat on the spread.** Median within-precursor fragment spread is 97 ppm at
the true apex, i.e. nearly the full +/-50 probe -- so a per-precursor mean still
mixes true fragments with wrong matches, and the sd of 13.1 overstates the
instrument's real scatter. The MEDIAN offset over 670 precursors is robust; the
width is not. A tighter probe, or a per-fragment rather than per-precursor
statistic, would sharpen it.

**Consequence:** fixing the mass calibration is worth doing and will help
uniformly -- but it will not recover the 9,275 unpicked precursors. For those,
the picker's own thresholds are the target: `min_corr_score` rejects 47% of scan
positions outright and the presence gate another 28%.

## MS1 co-elution is worth +25%, and the -8.4 ppm offset was mine (2026-08-08)

### MS1_COELUTION, wired and measured end to end

S08 + `lib_targets`, everything else identical, same binary:

    -no_ms1        1046 at 1% FDR
    default        1306 at 1% FDR      +260, +24.9%

**Far larger than predicted.** The probe measured a 1.4-1.5x bulk enrichment and
a 13.7x tail, and I wrote that the honest expectation was "a real but bounded
gain". +25% on a single feature is not bounded in that sense. The explanation is
that a classifier does not use a feature the way a threshold does: the tail
enrichment is what a FILTER would get, while the discriminant can exploit the
whole ordering, and this feature is the only one in the set whose errors are
independent of the other fifteen.

### The -8.4 ppm "instrument offset" was a probe artefact

Pinning it made things worse, which is what prompted the check:

    no offset                          1306
    -fragment_ppm_offset -8.4, wide    1219    -87
    -fragment_ppm_offset -8.4, narrow  1216    -90

The control that settles it: run the SAME anchored residual probe at a
retention time 300 s away, where the peptide cannot be.

    at DIA-NN's true apex     n=670   median -8.44 ppm   sd 13.13
    300 s away (control)      n=670   median -4.98 ppm   sd 14.96
    difference                                -3.46 ppm

**Most of the -8.44 is the probe, not the instrument.** Random matches inside a
+/-50 ppm window are not symmetric in ppm, so an unanchored median is biased
negative regardless of what the instrument does. The real offset is the
true-minus-control difference, about **-3.5 ppm** -- and pinning -8.4
over-corrected by ~5 ppm, hence the loss.

This is the third time this project has been caught by the same class of error:
a statistic computed over a large search space looks like a measurement and is
a property of the search. The fix each time is the same -- run the control that
cannot contain the signal. **Any m/z offset quoted from now on must come with
its RT-shifted control.**

Consequence: the mass calibration gate refusing to fire on S08 is closer to
correct than it looked. There IS an offset, it is ~-3.5 ppm against a +/-15 ppm
window, and it is worth roughly nothing compared to what a wrong one costs.

## OpenSWATH's picker BEATS ours, and every pinned offset loses (2026-08-08)

S08 + `lib_targets`, one binary, one config, only the picker changed:

    -picker openswath (defaults)        1368     <- best
    -picker coelution (default, ours)   1306
    -picker amplitude (ours)            1281
    -picker openswath -openswath_sn 0.1 -openswath_peak_width 25   1166

**The reference implementation wins by +62 (+4.7%)**, at its own default
settings, and my attempt to "tune" it for DIA (a lower signal-to-noise threshold
and an explicit 25 s peak width) cost 202. This is the doc/07 step 2 cross-check
that has been outstanding since the scoring plan was written, and the first time
it has been run.

That is uncomfortable, because the co-elution picker was this project's largest
single win -- it replaced the amplitude picker and fixed RT, FDR, recovery and
memory together. Both can be true: co-elution beat OUR amplitude picker (1306 vs
1281 here, and far more decisively when it was introduced), and OpenSWATH's
amplitude picker is better than ours because its boundaries come from
signal-to-noise rather than a fixed 10% of apex height. **The lesson is about
peak BOUNDARIES, not about detection.** Our fixed-fraction boundary is the
crudest part of the co-elution path and is the obvious thing to replace.

Caveat, stated because this document keeps having to: `lib_targets` is ~100%
true positives, so this is a RANKING benchmark. It measures which picker gives
the scorer better peaks among peptides that are present. It does not measure
search. The Astral arms and the `v6_50k` arm are running.

### Every pinned global offset loses

    no offset                              1306
    -fragment_ppm_offset -3.5              1164     -142
    -fragment_ppm_offset -8.4, wide        1219     -87
    -fragment_ppm_offset -8.4, narrow      1216     -90

Even the RT-shifted-control-corrected -3.5 ppm makes it worse. A single global
constant is simply the wrong model, and forcing one is worse than leaving the
window uncentred. **That is the argument for the RT-blocked recalibration**, not
against it: what these arms refute is the constant, which is the model the
existing MassCalibration fits.

### ...but on Astral the result INVERTS, and that is the real finding

    file / library            coelution (ours)   openswath (reference)
    S08 + lib_targets              1306               1368    (+62,  +4.7%)
    Astral + astral_lib_own        4290               2626   (-1664, -38.8%)

**The picker comparison is instrument-conditional, and a single-file conclusion
would have been exactly wrong.** Had only S08 been run, the honest reading would
have been "the reference picker is better, replace ours" -- and that would have
cost 39% on the other instrument.

The mechanism is consistent with what each picker is: OpenSWATH's sets
boundaries by signal-to-noise on the summed trace, which suits S08's
mobility-merged frames where the summed trace is comparatively clean; Astral's
higher-resolution, sparser spectra give co-elution far more to work with, and an
amplitude picker there fires on whatever is brightest in the window.

**Keep `-picker coelution` as the default.** OpenSWATH's is retained as a
selectable arm because it is the independent implementation doc/07 step 2
requires and because its BOUNDARY rule is measurably better on one of the two
files -- the fixed 10%-of-apex boundary in our co-elution path is the crudest
part of it and is the next thing to replace, taking the idea without the
amplitude detection.

This is the third instrument-conditionality trap in two days: the doc/08
prefilter's 1,800x came from non-mobility data and did not transfer; the
"uncentred narrow window is the worse error" rule was measured on S08 and did
not transfer to Astral. **Two files, always, before any picker or calibration
conclusion is recorded.**

### MS1 co-elution is worth far less on Astral

    S08 + lib_targets      1046 -> 1306   (+260, +24.9%)
    Astral + own library   4275 -> 4290   (+15,   +0.4%)

Not yet explained, and worth one measurement rather than a story. The obvious
candidates: Astral's bottleneck is availability (only 52.4% of covered
precursors yield any candidate at all), so a scoring feature cannot reach the
missing 47.6%; and the two runs differ in MS1 duty cycle. The S08 gain is the
one that has been reproduced.

## The m/z residual is not measurable with a standalone probe (2026-08-08)

Attempting "iterate recalibration strategies until the residuals are centred and
flat" stopped at step zero: **the residual cannot currently be measured well
enough to tell whether any strategy centres it.**

Per-fragment residuals, S08 + `v6_50k`, restricted to DIA-NN's TRUE precursors,
at their own apex, requiring >=4 co-occurring fragments, +/-20 ppm probe --
against the same thing 300 s away, where the peptide cannot be:

    estimator        at apex     RT-shifted control
    median            -6.94            -3.96
    MODE              -7.64           -13.90

**The control's mode is MORE extreme than the apex's.** Apex-minus-control by
cell is incoherent -- +0.6 to +8.4 ppm across RT octiles, -2.0 to +4.2 across
m/z octiles, no monotone structure in either. The probe yields 64,513 control
matches against 100,080 real ones, so roughly two thirds of what it measures is
chance, and both the median and the mode are dominated by that rather than by
fragments.

**So every m/z offset this project has quoted is unmeasured**, including the
-8.4 ppm and its "control-corrected" -3.5 ppm successor. That is consistent with
what the runs already said: pinning -3.5 cost 142 identifications and pinning
-8.4 cost 87-90, because both were fitted to the probe rather than to the
instrument.

**Why a standalone probe cannot work here.** It asks "is there a peak within
X ppm of this theoretical mass, somewhere near this time" -- and on a
mostly-absent library at DIA resolution the answer is yes by coincidence.
Tightening the window does not fix it; it biases the estimate toward zero by
clipping the true tail while still admitting chance matches near the centre.
This is the same extreme-value trap that made fragment depth, MS1 isotope depth,
`qualifying_spectra` and `total_matches` all saturate.

**The instrument that would work already exists and is thrown away.**
`ChromatogramExtractor.cpp:1068-1082` computes the deviation of every matched
peak and discards it. Those matches are constrained by co-elution and, for a
scored peak group at q<=0.01, by the whole discriminant -- so they are true
fragments in a way no standalone probe's matches are. Recording ppm per match,
and reading it only for confident peak groups, gives a residual whose
contamination is bounded by the FDR rather than by the search space.

**Next, in order:**
1. Instrument the extractor to carry per-match ppm (backlog item, already
   scoped: two parallel `[row][cycle]` blocks, sum(intensity*ppm) and
   sum(intensity), reduced at emit; a per-transition scalar would race).
2. Read residuals from peak groups at q<=0.01 and re-run this diagnostic. Only
   then is "centred and flat" observable.
3. THEN compare strategies -- constant, linear in log m/z (OpenSWATH's
   `weighted_regression`), quadratic (its `quadratic_regression_delta_ppm`),
   RT-blocked constant, RT-blocked linear (`MassRecalibration`, built and unit
   tested) -- and derive the extraction WINDOW from the residual quantiles,
   which both DIA-NN and OpenSWATH do and ODIA does not.

`SwathMapMassCorrection` is the reference for step 3: it extracts at the
calibrants' own elution time, regresses delta-ppm on m/z (optionally weighted,
optionally quadratic), and estimates the extraction window from the residuals.
Its anchoring is the part ODIA's `MassCalibration` lacks, and the window
estimate is the part nothing in ODIA has.

## Percolator wired in-process; three engines compared (2026-08-08)

**There is no mokapot in OpenMS -- checked the whole install.** What OpenMS 3.6
ships instead is better for a C++ tool: `Percolator` with a domain-agnostic
`RescoreInput` (features[n_rows][n_features], is_decoy, cv_group_keys) and a
`RescoreOutput` of scores/q-values/PEPs. In process, no Python, no file
interchange. mokapot is a Python reimplementation of the same semi-supervised
SVM; independent research confirms it is pure Python with a PIN-TSV file
interface and no C/C++ binding, so shelling out would be the only option there.

S08 + `lib_targets`, one binary, only `-classifier` changed:

    gbt (default, our hand-rolled histogram GBT)   1306
    percolator (OpenMS, cross-validated SVM)       1278    -28, -2.1%
    lda (ours)                                     1042   -264

**Percolator lands within 2% of our GBT and beats our LDA by 236.** Both are
linear, so that gap is implementation and cross-validation, not model class --
which is a useful independent check on our own LDA. Our GBT still edges
Percolator, consistent with a linear SVM being unable to represent the
interactions a tree ensemble can.

Percolator's own diagnostic on that run: 70,017 rows over 4,395
cross-validation groups, 16 features, pi0 0.335. `cv_group_keys` carries the
precursor grouping so a precursor's several candidate peak groups cannot be
split across folds -- they share a chromatogram, so splitting would leak.

### What the literature says about OUR failure, which is not the engine

Independently researched and adversarially verified. The finding that matters:
**mokapot's substantive difference from pyProphet is a PLUGGABLE LEARNER, not a
different FDR method** -- and ODIA already hand-rolled that substitution as its
GBT. Swapping engines was never going to fix the 1.5%-true-positive collapse,
and the measurement above says so: three engines, all within one regime.

The line that does speak to our regime is Noble/Keich/Kall:

* **Semi-supervised post-processing is documented to lose power and become
  unstable when confident positives are scarce.** That is our failure named in
  the literature rather than inferred from our own runs.
* **"Static modelling" is the stated remedy**: train once on a dataset where the
  loop does ignite, freeze the model, apply it. **OpenMS's Percolator already
  supports exactly this** -- `train(input)` returns a `PercolatorModel`,
  `score(input, model)` applies it, `saveModel`/`loadModel` persist it. So the
  remedy is reachable today with no new dependency.
* **RESET** -- a model-agnostic decoy-splitting wrapper addressing
  cross-validation label leakage -- can WRAP our existing LDA/GBT rather than
  replace either.
* Entrapment validation is the right check, and the same group sweeps ratios to
  0.04-0.16% true peptides, BELOW our 1.5%. Our fixture now plants 1.5%.

**Caveat carried from the research, not softened:** the surviving evidence base
is DDA/PSM-level. Nothing verified addresses peptide-centric DIA library search
at 1.5% library presence, so static modelling and RESET are well-motivated
transfers, not demonstrated results. crema, Triqler, PeptideProphet,
Oktoberfest/Prosit and transfer-learning rescorers produced no surviving claims
at all.

### The plug-and-play seam

`PeakGroupScorer.cpp` now branches on `options.classifier` between
`scoreSemiSupervisedLDA` and `scorePercolator`, which take the SAME
(features, labels, group) and return the SAME `ScoredGroups`. A common function
signature is the entire interface; a class hierarchy would add ceremony and make
A/B harder. NaNs -- legal in our sub-scores, meaning "not measured" -- are
replaced per column by that column's median before Percolator sees them, since
a NaN would propagate silently through an SVM dot product.

**Next: static modelling.** Train on `lib_targets` where the loop ignites, save
the model, apply it to `v6_50k` where it does not. That is the one strategy the
literature actually endorses for this regime, and it is now a few lines.

## The residuals are ALREADY centred and flat: -0.36 ppm, IQR 1.07 (2026-08-08)

`-collect_mass_residuals` keeps the m/z deviation the extractor computes to test
each match and has always discarded, records it per peak group as `Mass.Ppm`,
and is read only for groups the FDR has already accepted. S08 + `lib_targets`:

    group                       n        median ppm   IQR
    confident targets q<=0.01   4,743      -0.356     1.07
    decoys (built-in null)     47,617      -0.556     1.10

    flat in RT (octiles):   -0.19 -0.29 -0.32 -0.34 -0.43 -0.33 -0.54 -0.31
    flat in fragment count: -0.21 -0.43 -0.36 -0.45 -0.31 -0.47 -0.36 -0.27

**S08's fragment mass calibration is essentially perfect. There is nothing to
recalibrate**, and the loop's target -- "centred on zero and flat in RT and
m/z" -- is already met by the instrument. Total spread across RT octiles is
0.35 ppm, against an extraction window of 15,000 ppm-thousandths.

**This retroactively explains four earlier results.** The mass calibration gate
refusing to fire on S08 was correct: there is no offset to fit. Pinning -8.4 ppm
cost 87-90 identifications and pinning -3.5 cost 142 because the true offset is
-0.36, so both were 10x and 24x over-corrections. And the unanchored probe's
-8.44 ppm was never the instrument -- its RT-shifted control said so, and this
says so independently and with a hundredth of the error bar.

Identifications are unchanged at 1306 with the flag on, and peak RSS is 3.73
GiB, so the instrumentation is free in results and cheap in memory (one float
plane per live block, opt-in).

### The caveat that bounds this

**It measures matches that SUCCEEDED.** A fragment whose true m/z fell outside
the +/-15 ppm window records nothing, so this cannot by itself prove no fragment
is lost to miscalibration -- it is conditioned on being matched. What makes the
conclusion safe is the independent measurement from 2026-08-08: picked and
unpicked precursors have the SAME residual distribution (medians -8.21 vs -9.38
on the contaminated probe, identical fractions outside the window). Two
measurements with different selection effects agreeing that extraction is not
mis-centred is worth more than either alone.

Also note the decoys sit at -0.556 with the same IQR, which is at first sight
odd -- a shuffled sequence's fragments should match noise. They do not, because
a decoy peak group only survives the picker when real co-eluting peaks were
matched; both populations are therefore measuring the same instrument. That is a
consistency check, not a contamination.

### Consequence for MassRecalibration

`MassRecalibration` (RT-blocked, linear in log m/z, anchor-fitted, 10 unit
tests) is correct and **not needed on S08**. Keep it: it is the right model for
an instrument that DOES drift, it is unit-tested against a planted drift, and
the cost of having it is zero while the cost of rediscovering the need would be
another night. Wire it only if a run's `Mass.Ppm` shows structure this one does
not.

**Do not spend further effort on recalibration strategies for this data.** The
open problems are the picker (47% of scan positions rejected at
`min_corr_score`, 9,890 Astral precursors yielding no candidate) and the
1.5%-true-positive collapse -- neither of which is a mass problem.

## The extraction window is a SENSITIVITY knob, not a mass knob (2026-08-08)

Sweeping the fragment window at each instrument's own measured offset:

    window     S08 + lib_targets      Astral + own library
     4-5 ppm            0                     2,598
     6-8 ppm          853                     3,690
      10 ppm        1,054                       ---
      15 ppm        1,306  (baseline)         4,290  (baseline)

**Monotone on both instruments: narrowing always loses.** And mass accuracy does
not explain it. Per-fragment sigma is about 1.6 ppm (see the correction below),
so 10 ppm is ~6 sigma and should cost nothing. It costs 252.

The picker census says what actually happens, 15 ppm against 6 ppm on S08:

    scan positions evaluated      7,361,409  ->  6,246,070   (-15%)
    precursors yielding NO candidate     155  ->        410   (2.6x)

Narrowing does not reject worse peaks; it leaves fewer TRANSITIONS carrying any
signal, so precursors fall below the picker's `>=2 fragments present` and the
scorer's `min_fragments_at_apex >= 3` and never produce a candidate at all.

**So the window is not a mass-selectivity parameter in ODIA. It is feeding the
picker's fragment-count gates**, and a wide window is compensating for how
easily those gates starve. That is worth knowing before anyone tunes either:
the window and the count thresholds are one coupled system, and moving the
window alone moves sensitivity, not accuracy.

It also refutes the derived-window plan outright. DIA-NN and OpenSWATH size the
window from residual quantiles because for them it IS a mass parameter. Copying
that here would narrow the window to ~5 ppm on a correct reading of the
residuals and cost 40-70% of identifications.

### Two corrections to my own analysis, in one day

**(a) `Mass.Ppm` is a group median, not a fragment measurement.** It is the
median over a peak group's matched (fragment x cycle) cells, so its spread is
sigma/sqrt(N_eff) and NOT the per-fragment accuracy. Reporting "the instruments
are accurate to ~1 ppm" from it was reporting the precision of an average as the
precision of a measurement.

**(b) The first correction over-corrected.** `Mass.Ppm.N` has a median of 1,453,
and I read that as 1,453 independent fragments, giving sigma_fragment ~18 ppm.
But a precursor carries ~12 transitions; 1,453 counts CELLS, and cells of one
fragment across cycles share a calibration and are not independent. Effective
sqrt(N) is ~3.5, not 38, so sigma_fragment is about **1.6 ppm** -- which is
still small, and still does not explain the sweep.

Two wrong explanations in two rounds, both arithmetic on a quantity I had not
established the meaning of. **The per-fragment residual must be emitted
directly rather than inferred from a group statistic.** That is a small change
to code that already exists: keep the distribution, not its median.

## The picker's thresholds are NOT the lever -- measured, closed (2026-08-08)

"The picker rejects 47% of scan positions at min_corr_score" has been the top
open item for days. Both of its dominant thresholds are now swept, and neither
recovers anything. S08 + `lib_targets`, baseline 1306.

    -min_corr_score   0.0   0.3   0.4   0.5*  0.6    8.0
    identified       1306  1306  1306  1306  1306      0

    -apex_evidence   0.99*  0.90   0.70
    identified       1306   1230   1316
    apex rejections  2.74M  1.99M  0.90M

**min_corr_score does nothing between 0.0 and 0.6.** Not because the option is
dead -- 8.0 collapses the run to zero and drives below_corr from 3.65M to
4.11M, which is the falsification test that proves it live. It does nothing
because the best fragment's summed pairwise correlation is **bimodal**: either
negative (2.9M positions are rejected even at a 0.0 threshold, i.e. their best
correlation is negative) or comfortably above 0.6. Almost nothing lies between,
so the threshold sits in an empty region of the distribution.

**So the 50.7M rejections are genuine non-peaks, not lost signal.** That number
looked like a 47% loss waiting to be recovered and is nothing of the kind. The
filter is doing exactly its job.

**apex_evidence is nearly as inert and is non-monotone**: 0.99 -> 1306, 0.90 ->
1230, 0.70 -> 1316. Loosening it by a third of its range buys +10, and the
intermediate value LOSES 76, so it is interacting with candidate selection
(max_candidates and the margin rule) rather than acting as a simple gate.
Not worth changing a default on.

**Consequence:** the picker's thresholds are at reasonable settings and tuning
them is not a path to DIA-NN's 738. Combined with the window result -- the mass
window is a sensitivity knob feeding the picker's fragment COUNTS, not a mass
knob -- the whole extraction-side tuning surface is now measured and closed.

What remains is the 1.5%-true-positive collapse, where GBT, LDA and Percolator
all report zero, and where static modelling is the one literature-endorsed
remedy still untried.

### Method note

Three parameter values giving bit-identical output should have been read as
"the option is not reaching the code" and was, until a falsification test at an
extreme value distinguished "dead option" from "flat region of a bimodal
distribution". **Sweeping a plausible range cannot tell those apart; only a
value that MUST change the outcome can.** Worth doing first, next time, and it
costs one run.

## Window narrowing tested on the SEARCH benchmark: still worse (2026-08-08)

The previous window conclusion was drawn entirely from ~100%-present libraries,
where "fewer candidates" trivially means "lost true positives". The objection is
correct and important: **empty traces are EXPECTED on a realistic library**, and
narrowing is supposed to buy CLEAN traces rather than more of them. So the test
was re-run where it belongs, on S08 + `v6_50k` (1.5% present, DIA-NN finds 738),
scored by precision rather than by count:

    window   candidates   DIA-NN reachable   top-738 true   precision
    15 ppm      40,694        576 (78.0%)        232          31.4%
     6 ppm      31,134        448 (60.7%)        176          23.8%
     4 ppm      24,849        349 (47.3%)         94          12.7%

**Narrowing loses reachability AND precision.** The hypothesis that a narrow
window would trade quantity for cleanliness is not supported: it loses both. So
the earlier conclusion survives, but it now rests on the right benchmark and the
right metric instead of on a count from a library where every precursor is real.

### The gap to DIA-NN is not what "0 identifications" suggests

At 15 ppm, ranked by dscore, **ODIA's top 738 contains 232 of DIA-NN's 738 --
31.4% precision against a 1.5% base rate, a 21x enrichment.** The discriminant
is working. What fails is the q-value: it certifies none of them at 1% FDR.

That reframes the headline. "ODIA 0, DIA-NN 738" reads as no signal; the truth
is roughly one third of DIA-NN's precision at DIA-NN's own operating point, with
an FDR that cannot certify it. **The remaining work is calibration of the
q-values, not discovery of signal** -- which is exactly what static modelling
addresses and what three interchangeable engines (GBT, LDA, Percolator) could
not.

## OpenSWATH reference numbers, first run ever (2026-08-08)

`OpenSwathWorkflow` had never been run in this project. It now has.

* Astral + `astral_lib_own`: **54,289 features over all 10,891 precursors** --
  a feature for **100%** of the library, against ODIA's 52.4%.
* **Not FDR-controlled**: OpenSwathWorkflow emits features, pyprophet assigns
  q-values, and pyprophet is installed on no node. Do not quote it as
  "identifications at 1% FDR".
* Its TSV reader needs OpenSWATH column names; ours are DIA-NN's. Converter at
  `shared/to_osw_tsv.py`, runner at `shared/osw_run.sh`.

**The 100% versus 52.4% availability gap is the largest single deficit this
project has measured**, and it is on the extraction/picking side where the
tuning surface was just shown to be exhausted -- so it is structural, not a
threshold.

## Three-way reference at last: ODIA 4,290 | OpenSWATH 8,765 | DIA-NN 11,112

Astral, all three FDR-controlled at q <= 0.01. OpenSWATH had never been run in
this project; it now has, with pyprophet, and it is **twice ODIA**.

Four obstacles, all of which will recur and none of which is about ODIA:
1. OSW's TSV reader wants its own column names (`shared/to_osw_tsv.py`).
2. **Our libraries carry no decoys** -- ODIA generates them internally, so
   pyprophet refuses with "0 decoy and 10891 target groups". Fix:
   `OpenSwathDecoyGenerator -method shuffle`, 130,691 -> 261,370 transitions.
   Every earlier "OpenSWATH" figure quoted in this document was therefore
   unscoreable and should never have been offered as a comparison.
3. pyprophet 2.3.4 needs `pypdf<5`.
4. pyprophet 2.3.4 against modern numpy/pandas/scipy needs three patches:
   read-only array copies in `find_top_ranked`, `pandas<3`, and
   `rankdata(...).astype(int)`.

## Why S08 finds nothing: a self-reinforcing initialisation failure

Not the scorer, not the picker, not the window, not the mass model -- all of
which were tested at length and cleared. The ranking is FINE:

    top-N   targets  decoys  q=(D+1)/T   DIA-NN true   precision
       50        48       2     0.062          41        85.4%
      100        94       6     0.074          80        85.1%
      738       511     227     0.446         212        41.5%

**85% precision in the top 100.** What blocks certification is that a decoy sits
at rank 2 and six sit in the top 100, so q floors at 0.062 -- and Kall's +1
correction means even ZERO decoys needs >=100 clean targets before q can reach
0.01 at all.

**The mechanism is circular.** On `v6_50k` three of sixteen sub-scores are
IDENTICALLY ZERO -- `var_rt_delta`, `var_ms1_coelution`, `var_im_delta` -- while
on `lib_targets` the first two are 100% populated. Pass 1 has no retention-time
map, so it spreads the library evenly and extracts at approximately the wrong
times; the MS1 signal sits at the true elution time and the MS2 candidate does
not, so they never overlap and MS1_COELUTION is undefined. RT_DELTA likewise.
The scorer then runs on 13 features, missing the RT feature AND the only
orthogonal one, identifies nothing, and **no RT map is ever fitted**, so pass 2
never happens.

That also explains why swapping GBT for LDA for Percolator changed nothing: all
three were handed the same 13-feature input, and no learner recovers a feature
that is identically zero.

Two experiments are running: seeding pass 1 with an external iRT map fitted from
DIA-NN's 670 anchors (slope 7.774, intercept 732.7, p50 residual 9.1 s), and
applying a frozen Percolator model trained on `lib_targets`. Training on
`lib_targets` itself already beats every other engine there: **1,356** against
GBT's 1,306 and rescore's 1,278.

## The Astral deficit is EXTRACTION, and it is not the RT window (2026-08-08 night)

Counters added this round split `precursors_without_candidate` into its three
real causes for the first time. It had pooled extraction failures with picker
rejections, which is why the deficit read as a picking gap for two rounds.

Astral, 21,782 precursors (10,891 target + decoys):

    EXTRACTION losses:  0 with <3 points, 16,910 with an ALL-ZERO trace
    never reached the picker: 0 too-few-transitions, 0 too-few-cycles
    entered the picker and found nothing: 0

**Every precursor is assigned an isolation window and given cycles. The
chromatogram is allocated and no peak ever matches.** `no_points = 0` proves the
window assignment is not the problem; the fragments are simply not found in the
spectra searched.

### It is not the retention-time window

    baseline (p95-derived window)   4,290 IDs
    -rt_window_min 300              3,968 IDs   27,593 all-zero
    -passes 1                       1,248 IDs   16,910 all-zero

Widening the pass-2 window makes it WORSE. Two passes halve the all-zero count
against one (16,910 -> 9,890), so the RT map does help -- but the residue is not
recovered by looking in a wider window.

**Also: `-rt_window` is only the CAP.** `pass2_window = min(cap, max(floor,
p95_factor * p95))`, so 300/600/1200 s gave bit-identical output. The binding
knob is `-rt_window_min`. That is the second time this session a sweep returned
identical results because the option was not the one that binds -- the first was
`min_corr_score`, flat because its distribution is bimodal. **Identical output
across a sweep means "find out why" before it means "this does not matter".**

### What is left

The missing precursors have a window, have cycles, and have no matching peaks.
Since m/z coverage is fine (library 380-980 Th on an Astral run) and the mass
residual is ~1.6 ppm per fragment for the ones that DO match, the candidates are:

1. The isolation window a precursor is ASSIGNED to does not contain its true
   precursor m/z, so the right spectra are never searched.
2. The fragments genuinely are not there at the RT we look, i.e. the RT map is
   wrong for these specific precursors rather than globally.
3. OpenSWATH is finding a peak in noise where we correctly find nothing -- in
   which case its 8,765 at q<=0.01 would rest on peaks we would refuse. Its
   picker emits from any chromatogram; ours requires co-elution.

Candidate 3 is not idle: it would mean the availability gap is partly a
difference in what counts as a peak, not purely a deficit. Distinguishing it
needs a per-precursor join of OSW's features against our all-zero list, checking
whether OSW's q for those is good or marginal.

## The m/z extraction window was 3.3x too NARROW (2026-08-08 night)

**OpenSWATH's default `mz_extraction_window` is 50 ppm. ODIA's is 15
(`fragment_ppm_uncalibrated`), or 10 once calibrated.** That asymmetry had never
been compared, because every window sweep this project ran went NARROWER --
4, 6, 10, 15 ppm -- and read the monotone response as "narrowing is bad" rather
than as an arrow pointing the other way.

Astral, same library, only `-fragment_ppm` changed:

    window        IDs      all-zero traces
    15 ppm      4,290           ~26,800      <- our default
    30 ppm      4,499            23,203
    50 ppm      4,969            15,550      <- OSW's default

**+679 identifications (+15.8%) and 42% fewer all-zero traces.** The all-zero
count falling with width is the mechanism made visible: those precursors were
not absent, their fragment centroids simply fell outside a window we had set
too tight.

It also fits every other piece of evidence about the missing precursors: they
are systematically HIGHER m/z (missing fraction rises smoothly 40% at 400 Th to
64% at 950 Th, with no clustering at window boundaries, so not an assignment
bug) and 5.3x DIMMER by DIA-NN's own quantity. A dim, high-m/z fragment is
exactly the one whose centroid is least well determined.

### Why this was missed for so long

The `-fragment_ppm` sweeps were all bounded above by the default. Sweeping a
parameter's plausible range downward, seeing monotone degradation, and
concluding the current value is right is a mistake with a specific shape: it
never tests whether the current value is itself the constraint. **Sweep across
the reference implementation's value, not around your own.**

Two related lessons already recorded this session -- `min_corr_score` flat
because its distribution is bimodal, `-rt_window` inert because it is only the
CAP in `min(cap, max(floor, factor*p95))` -- have the same root: a sweep that
returns nothing informative usually means the experiment is wrong, not the
parameter.

### The curve turns over exactly at OpenSWATH's default

    ppm        15      30      50      75     100
    IDs     4,290   4,499   4,969   4,382   3,967
    all-0  26,800  23,203  15,550   7,341   2,436

**50 ppm is the optimum, which is precisely OSW's `mz_extraction_window`
default.** That two independent implementations land on the same number is
worth more than either measurement alone.

The all-zero count keeps falling monotonically all the way to 100 ppm, so wider
windows DO keep finding real signal -- 2,436 all-zero at 100 ppm against 26,800
at 15. Identifications nevertheless fall past 50 because the interference
admitted costs more than the signal recovered. The two curves separating like
that is the clean statement of the trade: availability and specificity are
genuinely opposed here, and 50 ppm is where they balance on this instrument.

## Ledger after the width fix: the bottleneck MOVED to scoring (2026-08-09)

Verified at the new default (50 ppm when the mass gate fails):
Astral 4,290 -> **4,969**, S08 unchanged at 1,306. Exactly as predicted.

    phase                    now      before      OSW
    library targets       10,891      10,891   10,891
    candidate emitted      9,479       5,704   10,891
    at q<=0.01             4,969       4,290    8,765

**Availability went 52.4% -> 87.0%.** The extraction deficit that dominated this
whole investigation is largely closed: 1,412 precursors still get no candidate,
against 5,187 before.

**But identifications rose only 679 on 3,775 new candidates.** The conversion
rate from candidate to identification is now the deficit:

    ODIA   4,969 / 9,479  = 52.4%
    OSW    8,765 / 10,891 = 80.5%

At OSW's conversion rate our current candidates would yield ~7,630. **That is
the headroom, and it is in scoring and FDR, not extraction.**

This is the ledger doing its job: fix the largest stage, re-measure, and the
constraint relocates. It also retro-justifies the earlier work that went
nowhere -- the scoring engines were interchangeable and the picker thresholds
inert *while availability was the binding constraint*. They are worth revisiting
now that it is not.

Next: MS1_COELUTION was worth +17.6% on S08 but only +0.4% on Astral -- measured
when Astral availability was 52%. With 87% it has far more to work with, so that
A/B is worth re-taking before anything else.

## Per-file parameters must be DETECTED, not defaulted (2026-08-09)

The two benchmark files want opposite extraction widths, and the difference is
not marginal:

    file     gate     50 ppm            10-15 ppm
    Astral   FAILS    4,969  (best)     4,290
    S08      PASSES   0 (collapses)     1,306  (best, calibrated ~10)

A single global default cannot serve both. Shipping 50 would destroy S08;
shipping 15 costs Astral 679 identifications. **The width is a property of the
run -- instrument, spectral density, whether frames are mobility-merged -- and
has to be measured per file.**

### What already exists, and why it is not enough

`MassCalibration` IS this mechanism: probe the run, measure the fragment mass
error, narrow the window to fit. When it works it is better than any constant --
S08 calibrated (~10 ppm) gives 1,306 against 922 for a hand-set 15 ppm.

It fails on Astral, and the fallback was the whole problem: a failed gate meant
"use 15 ppm", i.e. a moderately narrow window chosen for no reason. Fixed today
to 50 ppm on the principle that not knowing the error argues for a WIDE window,
not a middling one. But that is still a constant, just a better-chosen one.

### The design that would actually detect it

We now have the missing ingredient. `-collect_mass_residuals` keeps the m/z
deviation of every matched peak and reports it per peak group, read only for
groups the FDR has accepted -- so its contamination is bounded by the FDR rather
than by the search space, which is what every standalone probe failed at. That
gives a trustworthy per-run error distribution.

The loop writes itself:

    pass 1  extract WIDE (50 ppm), score, keep q<=0.01 groups
    then    take the per-fragment residual distribution from those groups
            width = a quantile of |ppm| (say 99th), floored and capped
    pass 2  extract at that width

This is what DIA-NN and OpenSWATH both do in spirit, and it removes the gate
from the critical path entirely: instead of a fragile peakedness test deciding
whether we are allowed to calibrate, the identifications themselves supply the
measurement. A run with no identifications stays wide, which is the correct
behaviour rather than a failure mode.

**Two things to be careful about, both already burnt once here.**

1. `Mass.Ppm` as currently reported is a MEDIAN over a group's (fragment x cycle)
   cells, so its spread is sigma/sqrt(N_eff) and NOT the per-fragment accuracy.
   Sizing a window from it would give something absurdly narrow. The width must
   come from the per-FRAGMENT distribution, which means emitting the
   distribution rather than its median.
2. Any quantile must be validated against an RT-shifted control, because a
   residual measured over a large search space is a property of the search. That
   discipline has already caught a -8.44 ppm "offset" that was -4.98 in the
   control.

### Also per-file, and not yet detected

* `precursor_im_window` -- meaningless on Astral (no mobility), load-bearing on
  S08.
* `min_fragments_at_apex` and `apex_evidence` -- swept as inert on S08 while
  availability was the constraint; unmeasured on Astral at 87% availability.
* The mass gate's own thresholds. It passes on S08 and fails on Astral, and
  which of those is "correct" was never established -- Astral's residual may be
  genuinely unmeasurable, or the gate may be mis-tuned for high-resolution data.

**Interim rule until detection exists: never quote a parameter as tuned without
naming the file it was tuned on.** Three instrument-conditionality traps have
already been recorded in this document.

## 2026-08-09: candidate generation is not the deficit, measured from every side

The loop's question was why OpenSWATH finds more candidates. It does not, and
no change to candidate generation recovers the gap. Astral unless noted; S08 is
`lib_targets`, baseline 1,306.

    baseline                                   4,969   1,306
    min_fragments_at_apex 1                    5,025   pending
    apex_evidence 0.50                         4,971
    max_candidates 50                          4,964
    max_candidates 12 / 6 / 3 / 1 (S08)                1,307 / 1,245 / 1,274 / 1,141
    OpenMS PeakPickerChromatogram              3,983   1,368
      ... at its DIA default sn 0.1            3,936
    amplitude picker                           3,848   1,281
    union (co-elution + amplitude)             3,002       0
    union + OpenMS picker at sn 0.1            3,817
    DIA-NN's published min_corr/corr_diff 1.0  pending  1,185

Read together: DEPTH is inert between 12 and 50 -- NOT below. Depth 1 costs
1,306 -> 1,141, a 12.6% loss, and depth 3 costs 32. Alternatives are needed;
what is inert is piling on more of them. DENSITY is harmful. A different
PICKER is neutral to harmful, and running OpenSWATH's at its own DIA S/N default
changes nothing. Tightening the margin to DIA-NN's published values costs 121 on
S08. The only gain all day came from REMOVING a gate that deleted rows over
`FRAGMENT_COVERAGE`, a feature the classifier already had.

Two hypotheses died here and should not be retried without new evidence:

* **Candidate depth.** Decoys average 18.44 candidates per precursor against
  targets' 8.95, because the margin rule is relative to each precursor's own
  best and a weak precursor admits nearly every position. That asymmetry is
  real, but it does not cost identifications: `max_candidates 1` equalises N and
  measures WORSE (1,141), and depth 12 to 50 is flat.
* **Picker fidelity.** `-openswath_sn 0.1`, kimi's lead from the OpenSwathWorkflow
  default, moved 3,983 to 3,936.

The bound that closes the question: 87.0% availability on 11,112 truth
precursors caps the candidate-generation gain at ~1,445 against a 4,969 -> 8,765
gap. Stated carefully, because the loose version of this went into a commit
message: **the ODIA-to-OpenSWATH deficit is 3,796**, and that is the number two
comparable pipelines differ by. The ~4,700 residual is against DIA-NN's 11,112,
which exceeds the 10,891-entry library, so the units do not obviously line up
and it should not be quoted as a count of recoverable precursors. Either way the
deficit is after candidates exist -- in scoring, ranking or FDR -- and that is
where the next phase goes.

## 2026-08-09: fragment mass accuracy, and the rule it establishes

Both reviewers independently named it as the highest-value score ODIA lacked,
and both gave the same mechanism: an extracted chromatogram records intensity
inside a mass window over time and discards WHERE in the window the peak sat, so
an interferent can co-elute perfectly, correlate perfectly, and be
systematically displaced in m/z with nothing chromatographic able to see it.

Two sub-scores, both already computed for the width measurement:
`var_mass_spread` (within-group scatter of per-fragment deviations,
calibration-free) and `var_mass_accuracy` (|deviation - the run's median|).

    Astral   5,025 -> 5,729   (+704, largest single gain measured on this file)
    S08      1,342 -> 1,168   (-174)

Leave-one-out on S08, one binary, via the new `-ablate`:

    both ablated      1,342     (reproduces the pre-feature number exactly)
    spread only       1,315
    accuracy only     1,285
    both active       1,168

**A FEATURE IS WORTH WHAT THE EXTRACTION HAS NOT ALREADY SPENT.** S08's gate
passes, so its window is +/-10 ppm centred on -9.35 with a per-fragment sigma of
1.10 -- the mass information has already been used as a filter and what remains
is noise the classifier overfits. Astral's gate fails, so its window is 50 ppm
and uncentred, none of the mass information has been spent, and the deviation is
the strongest thing available.

So `-mass_features auto` switches on the gate's own verdict, which is measured
from the data rather than set per file. This is the first working instance of
the per-file auto-detection in the 2026-08-08 backlog entry, and note what made
it work: the switch is an OBSERVABLE STATE OF THE PIPELINE, not a statistic
pooled over accepted identifications. The latter is what cost half the Astral
run when I tried it for the window width earlier the same day.

## 2026-08-09 13:30: the q-values are conservative, and the decoy null is why

Codex asked for one table: rank the targets by dscore, join the truth labels,
and print empirical FDR beside reported q. It separates three explanations that
had been indistinguishable. Scripts in `scripts/analysis/`.

### SEARCH (S08 + v6_50k, 738 true of 50,000) -- the 0 is HONEST

    rank   cum_true   emp_FDR   reported_q
      50         40     0.200       0.180
     100         62     0.380       0.322
     200         84     0.580       0.493
     738        113     0.847       0.701

Longest prefix at <=1% EMPIRICAL FDR: **2 precursors.** So reporting 0 at 1% is
not an FDR artefact -- there is genuinely nothing to certify. Reported q TRACKS
empirical FDR, slightly optimistic. The previous framing ("the discriminant
works and the q-values do not certify it") was half wrong: the discriminant
enriches 21x and the q-values are approximately right. **The discriminant is too
weak at a 1.5% prior**, and that is the whole of it.

Second finding, from the same table:

    decoy         p50 -0.802  p99 2.063  max 12.767
    FALSE TARGET  p50 -0.106  p99 3.468  max 15.044
    true target   p50  0.377  p99 14.452 max 15.405

**False targets outscore decoys.** The decoy null is easier than the real
negative class, which is exactly why reported q runs optimistic.

### Astral -- and here the SAME defect runs the other way

Every one of the 9,563 target precursors ODIA ranks is in DIA-NN's truth set.
There are NO false targets: `astral_lib_own.tsv` is pre-selected to what DIA-NN
found, so precision is 100% by construction and empirical FDR is 0.0000 at every
rank. **The Astral number is recall, not FDR-controlled discovery**, and any
statement of it must say so.

What that exposes: reported q <= 0.01 keeps 5,729 of those 9,563 -- so our own
FDR REFUSES 3,834 precursors that are all, in fact, true. The q-values are badly
CONSERVATIVE here. OpenSWATH reports 8,765 on this file; we rank 9,563 true
precursors and then decline to call them.

The mechanism is measured and it is the candidate asymmetry, now confirmed on
both files:

    Astral  target 10.98 candidates/precursor (median 6)   best-dscore p99 6.517
            decoy  20.99 candidates/precursor (median 24)  best-dscore p99 2.898, max 5.997

A decoy contributes the best of ~21 draws and a target the best of ~11, over
9,453 decoy precursors. The maximum of that many best-of-21 draws is what sets
the 1% threshold, and it reaches 5.997 against a target p99 of 6.517. Unequal N
breaks the exchangeability that target-decoy competition assumes.

Note this REVERSES the earlier dismissal. `max_candidates 1` equalises N and
measured worse (1,306 -> 1,141), and I concluded the asymmetry was harmless.
That experiment conflated two things: it equalised the null AND took the right
peak away from targets. The asymmetry does not hurt RANKING; it inflates the
DECOY NULL and makes q conservative.

**Next: build the null from decoys subsampled to the target candidate-count
distribution, leaving target scoring untouched.** Deterministic, non-circular
(subsample in canonical order, not by dscore), and it should move Astral toward
the 9,563 it has already ranked correctly.

## 2026-08-09 15:30: var_im_delta is a dead column, and the fix is now obvious

Every run today drops it: "dropping N sub-score(s) carrying no information:
... var_im_delta ..." appears in BOTH passes on BOTH files. `Options::observed_im`
(`include/odia/PeakGroupScorer.h:343`) was added as a plumbing point and has
never been filled, so `IM_DELTA` is NaN for every row and the constant-column
guard removes it. One of nineteen sub-scores is a slot pretending to carry
information.

**Why it is worth filling rather than deleting.** Ion mobility is a SEPARATE
PHYSICAL AXIS, not another statistic over the same twelve fragment traces. The
project's measured design rule is that orthogonality is the lever and count is
not -- 4 to 110 correlated features gave 0 identifications alike, while
`MS1_COELUTION`, the one feature reading a different channel, was worth +17.6%
on S08, and the fragment mass sub-scores were worth +704 on Astral for the same
reason. Mobility is the third such channel and it is already measured.

**The design, which is now known exactly because the mass work built it.**
`MobilityCalibration::Anchor` already carries `im_observed` per precursor, but
only for precursors identified in pass 1 -- far too few, and selected, which is
the trap that cost half the Astral run this morning. So do it the way
`collect_mass_residuals` does it: two extra float planes in `LiveSlot`,
Sum(intensity x 1/K0) and Sum(intensity), reduced at `emit()` to an
intensity-weighted observed 1/K0 PER CANDIDATE. Then `IM_DELTA` is
|observed - library| for every peak group, not just the identified ones.

Note what it cannot do: Astral has no ion mobility, so this is S08/diaPASEF
value only, and the larger gap is on Astral. Worth doing, not worth doing first.

Until it is filled, `var_im_delta` should be understood as absent rather than
uninformative -- the two look identical in the log and are not the same claim.

## 2026-08-09 15:50: iRT calibration, iteration 1 measured

Through the new `-stop_after calib`. Post-calibration residual over pass 1's
anchors; DIA-NN's comparable figure is `RT - Predicted.RT` at q <= 0.01.

    S08  anchor_q 0.05      n 1012  median -0.45  SD 71.23  robust 30.60  p95 103.53  max 1052.69
    S08  LOESS span 0.15    n 1012  median -1.57  SD 71.47  robust 30.94  p95 104.39  max 1062.38
    S08  LOESS span 0.30    n 1012  median -1.84  SD 74.34  robust 31.50  p95 119.71  max 1061.13
    S08  anchor_q 0.01      n  417  median -3.78  SD 41.36  robust 28.59  p95  91.63  max  243.72
    Ast  anchor_q 0.05      n 2469  median  0.35  SD 95.62  robust 34.22  p95  85.64  max 1586.51
    Ast  LOESS span 0.30    n 2469  median  0.67  SD 95.67  robust 34.45  p95  85.03  max 1585.69

    DIA-NN Astral, our library, q<=0.01
                            n 10891 median -2.51  SD 29.29  robust 26.09  p95  59.35  max  129.14

**LOESS is INERT.** +0.3 to +0.9 s on robust sigma on both files -- nothing, or
slightly worse. So codex's finding that it had never run was a real
documentation defect and NOT a performance defect, and the header's claim that
it is part of the fit was wrong in a way that cost nothing. `-rt_loess_span`
stays, defaulting to 0, because the negative is worth being able to re-take.

**Anchor quality is the live lever.** 0.05 -> 0.01 on S08 takes SD from 71.23 to
41.36 and max|e| from 1052 s to 244 s while robust sigma moves only 30.60 ->
28.59. That shape is diagnostic: the BULK of the map is fine and the tail is
misidentified anchors, not a bad fit.

**Which means the statistic is not yet comparable to DIA-NN's.** Ours is over
ANCHORS (peak groups pass 1 accepted at anchor_q); theirs is over CONFIDENT
IDENTIFICATIONS at q <= 0.01. Even at matched q the populations differ, because
an anchor is a peak group and an identification is a precursor. This is the
first thing to settle before claiming parity either way.

Astral stands at robust sigma 34.22 against the bar of 26.09
(`doc/14-irt-acceptance.md`). Not yet there.

## 2026-08-09 21:10: the exposed loop knobs do not ignite the SEARCH benchmark

v6_50k, all three arms, all zero:

    default                                                0
    train_fdr_initial 0.35, train_fdr 0.10, 6 iterations    0
    train_fdr_initial 0.35, use_pi0                         0

Kimi's Phase 1-2 finding #2 was right that LDAParams/NNParams/GBTParams were
reachable only by recompiling and worth exposing -- they are real, measurable
knobs now -- and wrong about the effect. Relaxing the semi-supervised loop does
not make it bootstrap on a 1.5%-prior library.

That should have been my prior, because the empirical-FDR table already said it:
the longest prefix at <=1% EMPIRICAL FDR on this benchmark is **2 precursors**.
There is nothing to certify, so no threshold or iteration count can certify it.
Only a better discriminant moves this.

Kimi's finding #4 (charge-3 fragments) is dead: every library on disk is 0.0%
charge >= 3, DIA-NN's included, against the claimed 35.6%.

Kimi's finding #1 (the decoy null) is partly wrong and partly open. Wrong: it
claims decoys keep the target's fragments; `LibraryGenerator.cpp:775-790`
recomputes every fragment from the decoy sequence, and records why -- shifting
b/y ions left 15.3% of decoy fragment m/z not matching the sequence stored
beside them, which is the asymmetric criterion an FDR must never have. Open: a
two-residue mutation still shares the target's retention time, mobility and
intensity pattern, and on the SEARCH benchmark false targets outscore decoys
(p99 3.468 against 2.063). `-decoys pseudoreverse` already exists, so the test
is a flag rather than new code, and is running.

## 2026-08-09 22:05: Phase 1-2 review triage — nothing survived

Seven findings across kimi and vibe (codex could not run: missing
codex-code-mode-host). Every one checked. None produced an actionable defect.

**kimi #4, charge-3 fragments — DEAD.** Claimed DIA-NN's library carries 35.6%
charge-3 fragments for 3+ precursors against our 25.8%. Every library on disk is
**0.0%** charge >= 3: `astral_lib_own`, `astral_lib_diann`, `v6_50k`,
`lib_targets`. Three minutes to check; days to have "fixed".

**kimi #2, hidden loop knobs — REAL but INERT.** LDAParams/NNParams/GBTParams
were reachable only by recompiling, which is worth fixing on principle and now
is (`-train_fdr_initial`, `-train_fdr`, `-classifier_iterations`, `-use_pi0`).
The predicted effect did not appear: v6_50k reports 0 at default, at
train_fdr_initial 0.35 / train_fdr 0.10 / 6 iterations, and with use_pi0.

**kimi #1, the decoy null — HALF WRONG, half under test.** Wrong that decoys
keep the target's fragments: `LibraryGenerator.cpp:775-790` recomputes every
fragment from the decoy sequence, and records that shifting b/y ions left 15.3%
of decoy fragment m/z not matching the sequence beside them. The empirical
asymmetry (false targets p99 3.468 against decoys 2.063) stands and
`-decoys pseudo_reverse` is measuring it.

**vibe F1, the double co-elution gate — ALREADY FIXED.** Its own confirming
experiment ("lower min_fragments_at_apex to 1") is what was done this morning:
Astral 4,969 -> 5,025, S08 1,306 -> 1,342. It also quotes 5,187 precursors
without a candidate, which is the pre-50-ppm figure.

**vibe F2, decoy capping by canonical order — MECHANISM WRONG.** It argues the
cap should take the top K by dscore instead, because canonical order may drop a
decoy's best candidates. That is backwards: taking the top K by score makes
best-of-K identical to best-of-all, which defeats the entire change. Canonical
order is a subset chosen INDEPENDENTLY of score, which is exactly what makes
best-of-K comparable between the classes. And the direction is measured the
other way in `odia_entrapment`: with matching on, FDP is 0.0020 against a
claimed 0.01 -- conservative, not anti-conservative.

**vibe F3, the censoring guard not enforced — WRONG.** It cites a
`fragmentPpm()` at OpenDIAlyzer.cpp:2102-2106 that does not exist. The guard is
at `OpenDIAlyzer.cpp:2124`:
`if (!mass_width_.valid || mass_width_.censored) { return fallback; }`.

**What to take from this.** A review round that finds nothing is a result, not a
failure, and manufacturing work from it would be worse than reporting it. But it
also says these reviewers are near the limit of what they can find from the
source alone: the real defects today came from MEASUREMENT -- the empirical-FDR
table, the leak guard, the resolved-device print -- not from reading code.

## 2026-08-09 22:30: pseudo-reverse decoys cost identifications, and we cannot yet say whether that is a correction

Astral, `-match_decoy_n`, same binary:

    decoys          base      + fine-tuned RT
    mutate         6,382          6,836
    pseudo_reverse 6,221          6,708

Pseudo-reverse costs 161 and 128 respectively -- about 2.5%. The direction is
what kimi's hypothesis predicts: a pseudo-reverse decoy is a HARDER null than a
two-residue mutation, so it scores higher, the 1% threshold moves up, and fewer
targets clear it.

**But fewer identifications is not evidence of a better FDR.** It is equally
consistent with pseudo-reverse being over-conservative. Astral cannot decide it:
its library is pre-selected to DIA-NN's hits, so it contains no false positives
to count. v6_50k reports 0 for both, so it cannot decide it either.

**What would decide it: entrapment on real data.** Spike a known-absent
population -- peptides from an organism not in the sample -- into the library at
a known ratio, run both decoy methods, and compare the false-discovery
proportion among the entrapment hits against the claimed q. The method whose FDP
matches its claim is the correct one. `test/tools/odia_entrapment.cpp` does this
synthetically and cannot see decoy GENERATION at all, which is exactly the part
in question.

Until then `mutate` stays the default, on the grounds that it is what DIA-NN
does and the alternative has no evidence behind it beyond scoring higher --
which is what an over-conservative null also does.

## 2026-08-10: the best configuration was not the default

Asked whether everything is wired in, the answer was no. Of the day's five
gains, three were defaults and two were opt-in flags, so a plain run gave 5,729
on Astral rather than 6,836:

    m/z window 15 -> 50 ppm            default   +679
    min_fragments_at_apex 3 -> 1       default    +56
    mass sub-scores (mass_features)    default   +704
    -match_decoy_n                     FLAG      +653
    -repredict_irt                     FLAG      +454

`-match_decoy_n` is now ON by default, disabled with `-no_match_decoy_n`. The
evidence supports it: +653 on Astral, +122 on S08, and the entrapment fixture
measures it CONSERVATIVE (FDP 0.0020 against a claimed 0.01), so leaving it
opt-in meant every default run paid for an asymmetry we know how to remove.

`-repredict_irt` cannot follow and should not. It needs a model fine-tuned on
THAT run, and a model reused across runs cost 2,027 precursors when it happened.
An opt-in flag is the correct shape for something that requires an artefact the
user has to produce first.

**The general lesson, which is worth more than the flag.** Every measured gain
should be checked against the default path the same day it is measured. A win
that lives behind a flag is a win nobody gets, and the gap is invisible in the
benchmark logs because the benchmark scripts pass the flags.

## The Astral decoy excess: located upstream of the picker (2026-08-13)

Decoys yield ~1.5x more peak groups than targets on Astral from a balanced
library, nothing is identified at 1% FDR, and on diaPASEF the same code is
balanced to a single group (142,321 / 142,320). FIVE explanations are now dead,
each by measurement rather than argument:

* **RT placement.** Target mean 30.00, decoy 30.47. Matched.
* **Fragment m/z in denser spectral regions.** -2.35 Th library-wide. Too small.
* **A library artefact.** Survives two independently generated libraries with
  different charges and m/z windows.
* **The `max_corr_diff` relative margin.** The best mechanistic proposal anyone
  made -- an absent precursor's "best" correlation IS noise, so a fixed margin
  around it admits nearly every position, and decoys are absent by construction.
  It predicts the right sign AND the diaPASEF contrast. It is still wrong:
  `outside_margin` fires TWICE in an entire run, and setting -max_corr_diff 0.0
  cuts candidates threefold while leaving the ratio at 1.536 against 1.539.
* **Skewed sampling.** Pass 1 samples every 499th precursor of an m/z-sorted
  library in which each decoy sits adjacent to the target whose precursor m/z it
  shares. Measured composition: 9,979 target / 10,029 decoy, ratio 1.005.

WHERE IT IS NOT, AND A RETRACTION. The per-class picker table below was read as
putting the imbalance at the FIRST stage, upstream of the picker. THAT WAS
WRONG. `PickerRejects` is `thread_local` and the reporter read one thread's
copy with no aggregation across 48 threads, so every ratio here is one thread's
slice of the work and nothing distributes chromatograms to threads in a
class-balanced way. Fixed by registering each thread's instance and summing.

The code, meanwhile, PROVES the classes must reach the scorer equally:
assignment is unconditional, the only precursor-level drop is `covering == 0`
which depends solely on precursor m/z, and every decoy shares its target's m/z
(verified 4,991,888 of 4,991,888). `scans` also cannot measure window width --
the `break` at the min_corr_score test exits the INNER fragment loop, not the
scan loop, so scans per precursor is fixed at n-2S-3 and `restrict_rt` is false
in pass 1 anyway, giving every precursor the identical whole-window range.

The peak-group excess itself is REAL and unaffected: it comes from the scored
result, not from these counters. Only its localisation was wrong. The numbers
below are retained as the record of a misleading measurement, not as evidence:

    scan positions        1,091,533 / 1,661,479   1.52x
    <2 fragments            428,178 /   642,481   1.50x
    below min_corr          509,556 /   786,456   1.54x
    reference zero          163,324 /   243,947   1.49x
    not a local max         221,241 /   347,826   1.57x
    below apex_evidence      75,535 /   115,838   1.53x
    outside max_corr_diff           2 / 0

Per precursor, over a balanced sample: 109.4 scan positions for a target, 165.7
for a decoy. `rej.scans` increments once per position in a precursor's
extraction window, so DECOYS ARE GETTING WIDER EXTRACTION WINDOWS. The picker is
symmetric; the asymmetry is already present when it starts.

The open question is therefore narrow and answerable: a decoy shares its
target's precursor m/z, hence its isolation window, and pass 1 spreads the
library evenly over the run -- so why is its window wider? Instrument
cycles-per-precursor by class next, and look at boundary clipping.

Note the counters could not have shown this before: they summed to 17,320,479
over 13,547,545 scan positions, 27.8% more than exist, mixing per-precursor and
per-scan-position granularity. They were never a partition.

### RESOLVED (2026-08-13): the excess is a zero-threshold on a zero-mean quantity

With the counters aggregated across threads and every stage split by class, the
whole library is accounted for and the divergence is at ONE gate:

    no points (<3)   4,990,048 / 4,990,026   1.00x   balanced
    empty trace          1,716 /     1,608   0.94x   balanced
    precursors reached     137 /       254   1.85x   <-- all of it

    of those WITH usable points:  1,853 target / 1,862 decoy   (1.005x)
    surviving to the picker:        7.4%     /   13.6%

Extraction is symmetric -- identical numbers of each class have points, exactly
as the code proves. The gate is:

    total[i] += (at[i] - median) * scale;    // median-subtracted, ZERO-CENTRED
    if (window_total <= 0.0) { empty_trace; return; }

Every transition is median-subtracted, so for a precursor with NO REAL PEAK the
summed trace is a zero-mean random variable and `<= 0.0` is a COIN FLIP. The
gate does not test for signal; it tests the sign of noise. Roughly half of all
absent precursors pass by chance, and because the threshold sits exactly on the
centre of the distribution it is maximally sensitive to any systematic
difference between the classes.

The systematic difference is the one dismissed earlier as too small: decoy
fragments average 2.35 Th LOWER in m/z, where Astral spectra are denser, so
marginally more of their points are non-zero and the sum tips positive slightly
more often. 2.35 Th is negligible against an 8 ppm TOLERANCE -- the comparison
made at the time -- and decisive against a threshold on a distribution's centre.

It also explains the diaPASEF contrast: mobility gating removes the chance
coincidences that make the sum jitter about zero, the coin stops flipping, and
the classes balance (142,321 / 142,320).

CONSEQUENCE. This is a defect in the gate, not a property of decoys. A test that
admits half of all absent precursors inflates the decoy null, and an inflated
null raises the 1% threshold above the few real targets -- which is exactly
"the scorer ran, the classifier trained, and the threshold rejected everything".
The fix is to test for actual signal (e.g. a positive excursion above local
noise in a minimum number of transitions) rather than the sign of a centred sum.
