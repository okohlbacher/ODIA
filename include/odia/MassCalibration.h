// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Per-run fragment mass calibration, derived from the run alone.
//
// Adopted from OpenDIAlyzer odia-v0.3.0 (`inferMassAccuracyPpm_`, and the
// robust estimators around it), same author and licence, with two deliberate
// changes recorded below.
//
// ---------------------------------------------------------------------------
// TWO SEPARATE PROBLEMS, AND THIS HEADER KEEPS THEM SEPARATE
// ---------------------------------------------------------------------------
//
// (A) WHERE the window sits. The observed m/z axis can be systematically wrong.
//     On S08_diaPASEF fragments sit about -10 ppm below theoretical: measured as
//     the centroid of the retention-time-specific excess over a local decoy-cell
//     null, -11.2 ppm on precursors ODIA recovers and -12.6 ppm on those it
//     misses, present at off-peak times too -- so an instrument term, not a
//     selection artefact. That is a RECALIBRATION of the mass axis, and its
//     natural shape on a TOF is a function of m/z, not a scalar.
//
// (B) HOW WIDE the window should be. That is the random error left AFTER the
//     recalibration, and it is what the reference implementation's
//     `inferMassAccuracyPpm_` sizes: window = max(k * robust sigma, floor).
//
// Conflating them is expensive in both directions. Narrowing a window that is
// centred on the wrong place throws away signal -- on S08, +/-10 ppm about zero
// keeps 0.51 of true fragments and the half it keeps is the tail, where +/-10
// about -10 keeps 0.81. And widening on a bad fit is worse still: in the
// reference, an ungated MS1 estimate of 72.6 ppm on an instrument measured at
// 1.66 ppm took identifications from 6,798 to 4,496.
//
// So `Model` carries a correction that is a function of m/z, and a window that
// is sized on the residual spread measured after that correction is applied.
// `sigma_before` and `sigma_after` are both reported so the recalibration's
// effect on the width is visible rather than asserted.
//
// ---------------------------------------------------------------------------
// WHY IT CAN RUN BEFORE ANY IDENTIFICATION EXISTS
// ---------------------------------------------------------------------------
//
// The reference fits from pass-1 confident peak groups. ODIA has none at
// extraction time -- that is the whole point of calibrating first -- so the
// anchors have to come from the data. The reference already records what
// happens if you simply widen the search and hope: handed the iRT anchor list
// over a 50 ppm bootstrap it measured a peakedness of 1.14-2.64 against a
// threshold of 3 and correctly refused to fit, because inverting the statistic
// puts that at ~3% real anchors. A mode fitted to 97% noise is meaningless.
//
// Purity therefore has to be bought some other way. Four things are tried here,
// and the honest accounting of what each one is worth on S08 is:
//
//   * ION MOBILITY. A diaPASEF frame's band is ~0.40 wide and a precursor
//     occupies ~0.05 of it, so testing each peak against the PRECURSOR's own
//     library 1/K0 removes most of the same-window interference. The reference
//     had no mobility dimension to use. Substantial.
//
//   * THE RUN'S OWN APEX. With no RT map, each precursor is probed at cycles
//     spread across the whole gradient and only its STRONGEST cell is kept.
//     For a real precursor that is its elution apex; for an absent one it is
//     the luckiest noise, which is exactly what the decoy control below
//     measures. Substantial.
//
//   * BRIGHTNESS (`min_intensity_quantile`). Decisive: it is the difference
//     between a gate that fails and one that passes, peakedness 1.96 on all
//     matches against 5.09 on the top quartile.
//
//   * CO-OCCURRENCE (`min_fragments_matched`). Nearly free, and measured to be
//     nearly worthless at this search width: at +/-50 ppm with 12 fragments per
//     precursor it passes 98.7% of target cells and 96.4% of control cells, so
//     it separates almost nothing. It is kept because it costs nothing and does
//     bite once the search narrows, but it is NOT where the purity comes from,
//     and saying otherwise would misattribute the result.
//
// The fifth lever is not a filter at all but where the gate is EVALUATED --
// see `gate_ppm`.
//
// The control is the second adopted-and-changed piece: every precursor is
// probed again with all of its fragments shifted by a few Da (`decoy_shifts`).
// Those cells cannot contain the real fragments but pass through the identical
// mobility test, co-occurrence rule and apex selection, so their residual
// distribution IS the null. This is the same construction that produced the
// -11.2 ppm figure above (probe3/present3.cpp), which is why it is trusted
// here.
//
// ---------------------------------------------------------------------------
// WHAT IS NOT PORTED, AND WHY
// ---------------------------------------------------------------------------
//
//   * The MS1 arm. ODIA extracts fragments; there is no MS1 window to size, so
//     the `ms1` parameter and its "fall back to the MS2 window" rule have no
//     referent here. `inferMassAccuracyPpm_`'s guard against counting the same
//     precursor observable six times over goes with it.
//
//   * `rt_trafo`. The reference locates its anchors in time through a fitted
//     retention-time map. There is none before the first pass, which is what
//     the apex selection above replaces.
//
//   * Adding |offset| back into the window. The reference deliberately does not
//     (the bias is corrected elsewhere, by SwathMapMassCorrection), and neither
//     does this -- here the bias is corrected by `Model` itself.
//
// ---------------------------------------------------------------------------
// WHAT THIS PRODUCES ON S08, 2026-08-05
// ---------------------------------------------------------------------------
//
// 3,000 sampled precursors, 160 stratified-random cycles = 3,840 of 32,210
// spectra (11.9% of the run), 84.6 s:
//
//   gate            peakedness 4.34 against a control at 2.98, threshold 3
//   residuals       7,066 target and 1,516 control, after the brightness cut
//   correction      LOG m/z: -9.35 ppm at 599 Th, +3.32 ppm per e-fold (t=7.4)
//   systematic      per-m/z-bin modes sit 1.71 ppm from a constant and 0.54 ppm
//                   from this -- 68% of the systematic error removed. The linear
//                   basis was fitted too and leaves 0.77, so log wins on the
//                   number rather than by assumption.
//   bin modes       -12.91 ppm at 238 Th rising to -7.43 at 1,102 Th
//   scatter         4.66 -> 4.19 ppm total per-hit, which is NOT what decided
//   rt drift        -0.11 ppm across 1,600 s at t=0.26, i.e. none
//   window          3 sigma would be 12.6 ppm, WIDER than the 10 ppm in force,
//                   so it is rejected -- calibration may only narrow
//
// The correction therefore runs from about -13.0 ppm at 200 Th to -6.4 at
// 1,473. A single constant near -9.9 would mis-centre by ~3 ppm at both ends in
// OPPOSITE directions, against a +/-10 ppm window, costing the lightest and
// heaviest fragments preferentially -- which is the case for modelling the
// shape, and it is not visible at all in the total scatter.
//
// Two independent checks that this is the instrument and not the method. An
// unrelated measurement over 124 M peak-transition hits at DIA-NN's retention
// times puts the offset at -9.78 ppm, the log-m/z fit at +2.79 ppm per e-fold,
// and finds no RT trend; all three agree with what is fitted here from the run
// alone. And sweeping the brightness cut from 0 to 0.9 walks the offset only
// from -10.51 to -9.50 ppm, so the selection is buying purity rather than
// manufacturing a number.
//
// On 12_80 with a library that does not match it, the same code refuses:
// peakedness 2.44 over 679 residuals against a control at 2.69, no offset
// applied, window left wide at 15 ppm, extraction proceeds normally.
//
#pragma once

#include <odia/Library.h>
#include <odia/SpectrumSource.h>

#include <cmath>
#include <functional>
#include <cstddef>
#include <limits>
#include <cstdint>
#include <string>
#include <vector>

namespace ODIA
{

  /// One signed mass residual: a candidate match found in a deliberately wide
  /// window, tagged with everything a shape fit needs.
  ///
  /// Deliberately NOT pre-aggregated. The whole question is whether the error
  /// depends on m/z or on retention time, and a histogram collapsed over either
  /// axis cannot answer it.
  struct MassResidual
  {
    float mz = 0.0f;          ///< theoretical fragment m/z, Th
    float rt = 0.0f;          ///< retention time of the spectrum, seconds
    float ppm = 0.0f;         ///< (observed - theoretical) / theoretical * 1e6
    float intensity = 0.0f;   ///< of the matched peak
    /// The PRECURSOR's 1/K0, not the fragment's -- a fragment has no mobility
    /// of its own, it inherits the packet it was produced in. NaN when the run
    /// has no mobility or the precursor has no library value.
    ///
    /// Carried because MaxQuant's mass model has an ion-mobility term worth
    /// 52% of its modelled variance on timsTOF (Prianichnikov et al., MCP 2020)
    /// -- but that is a DDA PRECURSOR result, and whether it transfers to DIA
    /// fragments is exactly what this field exists to let us measure. It is not
    /// fitted anywhere yet, deliberately. See doc/15 sections 12.1 and 13.
    float im = std::numeric_limits<float>::quiet_NaN();
    bool decoy = false;       ///< from an m/z-shifted control cell
  };

  /// Fits the run's fragment mass error, and sizes the extraction window from
  /// what is left after correcting it.
  class MassCalibration
  {
  public:
    struct Options
    {
      /// Half-width of the search used to COLLECT residuals, ppm.
      ///
      /// Deliberately much wider than any window that will be extracted with:
      /// the distribution's shape and its shoulders both have to be visible, and
      /// a search narrower than the error truncates the thing being measured.
      /// The GATE is evaluated over `gate_ppm` instead -- see below for why the
      /// two cannot be the same number.
      double search_ppm = 50.0;

      /// Half-width the LOCATION and the GATE are evaluated over. Distinct from
      /// `search_ppm`, which is what gets collected.
      ///
      /// They are separated because they want opposite things. The search must
      /// be wide or the distribution is truncated and its scale censored. The
      /// gate must be narrow, because a flat background contributes to its edge
      /// band in proportion to the width searched -- for a narrow peak on a
      /// uniform background the statistic works out at 1 + 5N/M with M growing
      /// linearly in the width -- so a wide window dilutes the very contrast the
      /// gate exists to measure. Measured on S08's own residuals, same sample:
      /// peakedness 3.18 evaluated over +/-50 against 5.09 over +/-30.
      ///
      /// 30 rather than tighter because the gate's edge band must still be
      /// background: at +/-30 it sits 18-24 ppm from the mode, more than 5
      /// sigma out on this run.
      double gate_ppm = 30.0;

      /// Keep only matches whose peak intensity is at or above this quantile of
      /// all matched intensities. 0 keeps everything.
      ///
      /// This is the reference's own argument about which peak to pick, applied
      /// to which matches to KEEP: the interferent population is numerous but
      /// individually weak, so brightness is evidence. It is the difference
      /// between a gate that fails and one that passes on S08 -- peakedness 1.96
      /// on all matches against 5.09 on the top quartile.
      ///
      /// The check that this SELECTS purity rather than MANUFACTURES an answer
      /// is that the answer does not move. Sweeping the quantile from 0 to 0.9
      /// on S08 walks the fitted offset from -10.51 to -9.50 ppm -- a 1 ppm
      /// range against a 10 ppm effect -- while peakedness goes 1.96 -> 5.88. A
      /// cut that was inventing the peak would move the estimate with it.
      double min_intensity_quantile = 0.75;

      /// Precursors sampled from the library, spread evenly through it.
      /// A scalar-plus-slope needs hundreds of anchors, not thousands.
      std::size_t max_precursors = 3000;

      /// Acquisition cycles probed. STRATIFIED RANDOM across the gradient: the
      /// run is cut into this many equal strata and one cycle is drawn at
      /// random from each.
      ///
      /// One cycle is a contiguous run of spectra covering every isolation
      /// window, so a cycle is one range request and every window gets probed at
      /// the same instant. It is also how the sample is spread in retention
      /// time, which is what makes an RT drift measurable at all.
      ///
      /// Stratified random rather than either alternative, and both alternatives
      /// are wrong in a way that matters here. A PREFIX (or any contiguous
      /// block) measures one stretch of the gradient and cannot see a drift at
      /// all. A fixed STRIDE covers the gradient but can alias against anything
      /// periodic in the acquisition -- and a DIA run is periodic by
      /// construction. Drawing at random inside each stratum keeps the guaranteed
      /// coverage and destroys the aliasing.
      ///
      /// On S08: 160 of 1,342 cycles = 3,840 of 32,210 spectra, 11.9% of the
      /// run, 84 s against a ~10 min full extraction pass.
      std::size_t cycles = 160;

      /// Seed for that draw. Fixed, not clock-derived: a calibration that
      /// returns a different number each time it is run cannot be checked
      /// against a previous run, and "the offset moved" would be
      /// indistinguishable from "the sample moved".
      std::uint64_t sample_seed = 0x0D1A0805u;

      /// Fragments used per precursor, most intense first in library order.
      std::size_t max_fragments = 12;

      /// A cell counts only when this many of one precursor's fragments match
      /// in the SAME spectrum. This is the co-occurrence purity filter; at 1 it
      /// is switched off and the sample reverts to the reference's behaviour.
      std::size_t min_fragments_matched = 3;

      /// Below this many residuals the fit is not attempted.
      int min_residuals = 200;

      /// Reject the fit unless residual density near the mode exceeds density
      /// at the edge of the search window by this factor. Uniform (all-noise)
      /// residuals give ~1. Ported unchanged from the reference, including the
      /// threshold, because the threshold is what caught the 72.6 ppm fit.
      /// Kept because the number is still worth reporting, but it no longer
      /// gates anything -- see the note in `fit`. Shape was a proxy for "are
      /// these real fragments"; `accept_ratio` tests that directly.
      double min_peakedness = 3.0;

      /// ACCEPTANCE: the FRACTION of the centred systematic error the
      /// correction must remove, measured about zero on the probe's residuals.
      ///
      /// Two wrong versions preceded this one and both are worth recording.
      ///
      /// A pure ratio at 0.50 accepted PURE NOISE: re-centring any distribution
      /// moves its bin modes toward zero, so a uniform sample "improved"
      /// 21.11 -> 8.83 ppm (ratio 0.42) and would have had a 29.83 ppm constant
      /// applied. The bar was simply too loose, not the statistic wrong.
      ///
      /// An ABSOLUTE bar at 1.0 ppm then rejected S08, whose model is known
      /// good. It was calibrated on residuals left over ID-ANCHOR populations
      /// (0.18-0.21 ppm) and over synthetic fixtures (0.03-0.26), but it is
      /// applied to the PROBE's residuals, which are noisier by construction --
      /// S08's probe leaves 2.03 ppm there while the same model leaves 0.55 ppm
      /// on held-out identifications. A threshold in ppm cannot be transplanted
      /// between populations; a fraction can.
      ///
      ///     Astral probe    1.394 -> 0.198   ratio 0.142
      ///     S08 probe      12.040 -> 2.027   ratio 0.168
      ///     fixtures                         ratio 0.003 - 0.025
      ///     pure noise     21.105 -> 8.833   ratio 0.419
      ///
      /// 0.25 sits 1.7x below the noise case and 1.5x above the worst real one.
      /// That is a narrower margin than the project usually accepts and it is
      /// set on two real runs, so it should be re-checked on a third instrument
      /// before it is trusted as general.
      double accept_ratio = 0.25;

      /// ...and an absolute backstop, so a model that removes 80% of an
      /// enormous error and leaves a still-enormous one cannot pass on the
      /// ratio alone.
      double accept_after_ppm = 5.0;

      /// Minimum control residuals before the "control is as peaked as the
      /// data" rule may fire at all.
      ///
      /// `peakednessRatio` is count(central band) / count(edge band), so its
      /// relative error is ~sqrt(1/c + 1/e) and it is dominated by the EDGE
      /// count. Measured on Astral: 5,769 control cells yielded **98** control
      /// residuals against 4,033 target, because a 7 Th-shifted query rarely
      /// matches anything at all. The statistic came out at 12.00 with an edge
      /// count of one to three, i.e. an error bar of roughly +/-7 to +/-12,
      /// and the gate failed the run by comparing the target's well-estimated
      /// 6.90 +/- 0.74 against it.
      ///
      /// That is not a test, it is a coin flip.
      ///
      /// DEFAULTED OFF (0) ANYWAY, because the coin flip was landing the right
      /// way up. Setting this to 400 let the gate PASS on Astral pass 1, which
      /// applied a -1.27 ppm offset and narrowed the window off the
      /// uncalibrated 15 ppm -- and pass 1 fell from **1,248 identifications to
      /// 119**, taking the run from 4,275 to 1,994.
      ///
      /// So the thin control was vetoing a calibration this data genuinely does
      /// not support. It is a bad ESTIMATOR and, here, a correct DECISION:
      /// almost nothing matching a 7 Th-shifted query is itself evidence that
      /// the fragment matches are not clean enough to calibrate from. The
      /// argument that a noisy statistic "has no opinion" ignored that its
      /// SPARSITY is the signal, independent of the ratio computed from it.
      ///
      /// A real fix measures the null better -- more control cells, or a
      /// shift chosen to match at a comparable rate -- rather than ignoring it.
      /// Do not raise this without re-running Astral end to end.
      std::size_t min_control_residuals = 0;

      /// The run's own iRT map, so the probe looks only where a precursor
      /// should elute: rt = irt_slope * iRT + irt_intercept, +/- rt_window.
      ///
      /// NOTE: the first measurements taken through this path (the Astral
      /// "control 98 -> 149 residuals" and "4,275 -> 1,895" figures) were
      /// produced with an out-of-bounds index in the gate and are VOID. See the
      /// fix commit; they must be re-taken.
      ///
      /// Zero means no map, which is the FIRST pass -- there is none yet, and
      /// the brightest-cell apex stand-in is all that is available. From the
      /// second pass there is one, and using it matters: without it the probe
      /// keeps the brightest cluster over every block it looks in, and an
      /// ABSENT precursor gets one draw from the interference per block. That
      /// is how a mostly-absent library produces flat residuals and fails the
      /// gate for a reason that is about the probe rather than the instrument.
      double irt_slope = 0.0;
      double irt_intercept = 0.0;
      double rt_window_seconds = 0.0;

      /// Window = k x robust sigma of the CORRECTED residuals.
      /// Window = this many robust sigmas of the CORRECTED residual.
      ///
      /// 3 -> 8 on 2026-08-10, from a k = 1..10 sweep on both benchmark files
      /// with the offset and shape held fixed and only the width varying
      /// (-mass_calibration_offset_only, otherwise every arm clamps to the
      /// model's own width and the sweep measures nothing):
      ///
      ///     k      1    2    3    4    5    6    7    8    9   10
      ///     S08  226  665  721  832  817  809  787  842  743  775
      ///     Ast    0  318  617  866 1069 1125 1172 1148 1194 1222
      ///
      /// S08 plateaus from k=4; Astral never turns over. k=8 maximises the
      /// worst case across the two (93.9% of each file's own maximum).
      ///
      /// TWO THINGS MAKE THIS PROVISIONAL, BOTH MEASURED:
      ///
      /// (a) It is a REGRESSION ON ASTRAL against what ships today. The gate
      ///     currently fails there, the window stays at the uncalibrated 50 ppm,
      ///     and that yields 2,078 identifications -- more than every arm in the
      ///     sweep. Astral's optimum is out near 50 ppm (k~25) and the sweep
      ///     simply did not reach it.
      ///
      /// (b) The sigma this multiplies is NOT the sigma the sweep was scaled on.
      ///     The sweep used the robust sigma of FDR-accepted fragments: 2.087
      ///     ppm on S08, 2.012 on Astral -- near-identical. `sigma_after` here is
      ///     the probe's, measured through a 50 ppm search on a brightness-cut,
      ///     mobility-gated population: 4.19 ppm on S08 and 0.96 on Astral, a
      ///     factor of FOUR apart on runs whose true precision agrees to 4%.
      ///     So one multiple times this sigma cannot mean the same thing on two
      ///     instruments, and k here is not the k of the sweep.
      ///
      /// The honest fix for (b) is to size the window from identifications
      /// (`MassWidth`, -mass_width_from_ids) where the sigma is comparable
      /// across runs, not from the probe. That path exists and its `apply` mode
      /// was measured harmful at 3 sigma -- which this sweep now explains.
      double sigma_multiple = 8.0;

      /// Never infer a window below this, whatever the fit says.
      double floor_ppm = 2.0;

      /// Half-width of the ion-mobility test against the precursor's library
      /// 1/K0. 0 disables it; it is also skipped per precursor when the library
      /// carries no 1/K0, because absent information is not evidence of
      /// mismatch.
      ///
      /// DELIBERATELY TIGHTER than the extractor's `precursor_im_window`
      /// (0.025), which is the window this offset is ultimately applied
      /// through. That looks wrong -- fitting through a narrower band than you
      /// correct through -- and it is not, because the two windows are sized
      /// for different jobs. The extractor wants COMPLETENESS: 0.025 is ~2.6
      /// sigma on the library-vs-observed agreement, so it keeps the peptide's
      /// real fragments. This wants PURITY: a mode in the residual histogram
      /// that an m/z-shifted control does not also have. The offset itself is a
      /// property of the mass axis, not of the mobility axis, so estimating it
      /// on the cleanest subset is unbiased -- what it costs is residuals, and
      /// there are thousands.
      ///
      /// Measured on S08, 3,000 precursors over the default 160 cycles:
      ///
      ///     0.010   GATE PASSED, peakedness 5.03 against a control at 2.94,
      ///             -9.02 ppm at 668.9 Th + 2.68 ppm per e-fold
      ///     0.025   GATE FAILED, peakedness 4.41 against a control at 4.60 --
      ///             the shifted control is as peaked as the data
      ///     0.050   GATE FAILED, residuals FLAT (2.89)
      ///
      /// So widening this to match the extractor does not loosen the fit, it
      /// removes it: at 0.025 the same-window interference the band cannot see
      /// makes the target sample indistinguishable from its own null, and the
      /// run gets no calibration at all. The two numbers are the same decision
      /// taken against different denominators, and reconciling them by making
      /// them equal would be reconciling away the measurement.
      double im_window = 0.010;
      bool use_ion_mobility = true;

      /// m/z shifts, in Th, that build the null. Each shifts ALL of a
      /// precursor's fragments together, so a decoy cell passes exactly the same
      /// mobility, co-occurrence and apex tests as a target cell and differs
      /// only in that it cannot contain the real ions.
      ///
      /// The values are the ones the S08 measurement used, and they are chosen
      /// to be neither a common neutral loss nor an isotope spacing.
      std::vector<double> decoy_shifts{7.33, -7.19};

      /// An m/z-dependent term is accepted only when it clears all three tests
      /// below. They ask one question between them: how much of the SYSTEMATIC
      /// error does modelling the shape actually remove?
      ///
      /// It is deliberately NOT judged on the total per-hit residual, and that
      /// correction is worth recording because the first version of this class
      /// got it wrong. Total scatter on S08 is ~5.1 ppm and is dominated by
      /// irreducible per-fragment noise that no calibration can touch. Judging a
      /// systematic correction against that denominator understates it
      /// structurally: the trend can be almost perfectly removed and still move
      /// the total by ~11%, because the total was mostly never going to move.
      /// That criterion refused a correction on S08 that removes 64% of the
      /// systematic term.
      ///
      /// So the comparison is between the per-m/z-bin MODES and the model: the
      /// weighted RMS of the bin modes about the fitted shape, against their
      /// weighted RMS about the best constant. On S08 that is 1.71 ppm about a
      /// constant and 0.54 ppm about the log fit -- 68% of the systematic error
      /// removed, where the total-scatter view saw 4.66 -> 4.19 and called it a
      /// wash.
      ///
      /// `min_slope_t`: the slope must be this many standard errors from zero.
      /// Necessary, not sufficient.
      ///
      /// `max_systematic_ratio`: the shape must leave at most this fraction of
      /// the systematic spread the constant leaves.
      ///
      /// `min_systematic_gain_ppm`: and it must do so by an absolute margin. A
      /// two-parameter fit always beats a one-parameter fit on the same bins, so
      /// on a run whose residual really is flat the ratio alone would eventually
      /// wave a slope through -- the bins would just be noise, and fitting noise
      /// better is not an improvement. The absolute floor is what makes that
      /// impossible: on a flat residual the constant already leaves only the
      /// bin-estimation error (a few tenths of a ppm), so no shape can gain half
      /// a ppm on it.
      double min_slope_t = 3.0;
      double max_systematic_ratio = 0.70;
      double min_systematic_gain_ppm = 0.50;

      /// Worker threads for the matching. 0 uses the hardware concurrency.
      unsigned threads = 0;

      bool verbose = false;
    };

    /// The fitted correction, as a function of m/z, plus the window it implies.
    ///
    /// `form` says which model was chosen, and the choice is recorded rather
    /// than inferred from the coefficients: a linear fit that happens to return
    /// a slope of zero and a constant fit are the same numbers and different
    /// decisions.
    struct Model
    {
      bool fitted = false;
      /// "none" | "constant" | "log_mz" | "linear_mz"
      std::string form = "none";
      std::string reason;               ///< why, in words, for the log

      double intercept_ppm = 0.0;       ///< correction at `reference_mz`

      /// The two candidate bases. At most one is ever non-zero; both zero is the
      /// constant form. They are separate named fields rather than one slope
      /// plus a basis flag because the UNITS differ, and a single field whose
      /// meaning depends on a neighbouring enum is exactly the silent unit error
      /// this codebase separates `im` from `ccs` to avoid.
      ///
      /// Which basis is used is decided by measurement, not assumed. On S08 the
      /// log fit leaves 0.54 ppm of bin-mode residual against the linear fit's
      /// 0.77, so log wins -- but it wins on the number, and a run that prefers
      /// the linear basis gets it (the synthetic suite checks both directions).
      double log_slope_ppm = 0.0;             ///< ppm per e-fold in m/z
      double linear_slope_ppm_per_1000 = 0.0; ///< ppm per 1000 Th
      double reference_mz = 700.0;

      std::size_t residuals = 0;        ///< target-cell residuals used
      std::size_t decoy_residuals = 0;

      double sigma_before = 0.0;        ///< robust sigma of raw residuals, ppm
      double sigma_after = 0.0;         ///< ... after the chosen correction
      /// Total per-hit scatter under each candidate. REPORTED ONLY -- these no
      /// longer decide anything. They are kept because the difference between
      /// them and `systematic_*` below is the whole point: on S08 they say 5.11
      /// against 4.56 and look like a wash, while the systematic numbers say
      /// 1.40 against 0.51 and say the opposite.
      double sigma_constant = 0.0;
      double sigma_shaped = 0.0;

      /// What DOES decide: weighted RMS of the per-m/z-bin modes about the best
      /// constant, and about the chosen shape. In ppm of systematic error.
      double systematic_before = 0.0;
      double systematic_after = 0.0;
      /// The acceptance quantities: weighted RMS of the per-m/z-bin modes about
      /// ZERO, before and after the model. Distinct from `systematic_*`, which
      /// are taken about a refitted constant and therefore cannot see a residual
      /// offset -- the very thing a window mis-centres on.
      double centred_before_ppm = 0.0;
      double centred_after_ppm = 0.0;
      /// The losing basis's systematic residual, so the basis choice is visible.
      double systematic_log = 0.0;
      double systematic_linear = 0.0;

      double slope_t = 0.0;             ///< standard errors the slope is from zero
      double shape_swing_ppm = 0.0;     ///< what the shape predicts across the range
      double peakedness = 0.0;          ///< target cells
      double decoy_peakedness = 0.0;    ///< the same statistic on the null

      /// The thresholds this model was judged against, carried so a report can
      /// state them beside the numbers instead of restating constants that then
      /// drift apart from the Options they came from.
      double gate_peakedness = 0.0;
      double gate_ppm = 0.0;
      double sigma_multiple = 0.0;
      double floor_ppm = 0.0;
      double search_ppm = 0.0;

      /// Window half-width the fit supports, ppm. <= 0 means "not inferred --
      /// keep whatever the caller configured".
      double window_ppm = -1.0;

      /// Diagnostic only: how far the corrected residual drifts, in ppm, across
      /// the retention-time span actually sampled, and how many standard errors
      /// that slope is from zero. Reported, never applied -- see the .cpp.
      double rt_drift_ppm = 0.0;
      double rt_drift_t = 0.0;

      /// m/z range the residuals actually covered, which is where the model is
      /// supported and outside which a slope is extrapolation.
      double mz_low = 0.0, mz_high = 0.0;

      /// The correction to add, in ppm, at this fragment m/z.
      double ppmAt(double mz) const
      {
        if (!fitted) { return 0.0; }
        double ppm = intercept_ppm;
        if (log_slope_ppm != 0.0 && mz > 0.0 && reference_mz > 0.0)
        {
          ppm += log_slope_ppm * std::log(mz / reference_mz);
        }
        if (linear_slope_ppm_per_1000 != 0.0)
        {
          ppm += linear_slope_ppm_per_1000 * (mz - reference_mz) / 1000.0;
        }
        return ppm;
      }
    };

    /// Per-bin summaries, so the shape can be reported as numbers rather than
    /// as a fitted line the caller has to take on trust.
    struct Bin
    {
      double centre = 0.0;      ///< mean m/z (or RT) of the bin
      double location = 0.0;    ///< robust centre of its residuals, ppm
      double stderr_ppm = 0.0;
      std::size_t n = 0;
    };

    struct Diagnostics
    {
      std::vector<Bin> by_mz;
      std::vector<Bin> by_rt;
      std::size_t cells_probed = 0;
      std::size_t target_cells = 0;
      std::size_t decoy_cells = 0;
      std::size_t spectra_decoded = 0;
      double collect_seconds = 0.0;
    };

    /// Probe the run and return the raw residuals, targets and controls.
    static std::vector<MassResidual> collect(const Library& library, SpectrumSource& source,
                                             const Options& options,
                                             Diagnostics* diagnostics = nullptr);

    /// Fit a model to residuals from anywhere -- a run, or a synthetic sample.
    /// Weighted RMS of the per-m/z-bin modes left behind by @p correction, ppm.
    ///
    /// This is the statistic the basis choice turns on -- on S08, 1.71 ppm
    /// about a constant against 0.54 about the log fit, while total per-hit
    /// scatter moved only 4.66 -> 4.19 and looked like a wash.
    ///
    /// Exposed because `fit` computes it against its OWN sample, which measures
    /// how well a model describes the data it was fitted from. Scoring a model
    /// on residuals it never saw is a different question and the one that
    /// matters here: a model fitted from identifications is fitted from a
    /// selected sample, and the only honest check is held-out. See doc/15
    /// section 6.2.
    ///
    /// The evaluation population is fixed ONCE from the UNCORRECTED residuals
    /// (mode +/- 4 sigma) and every candidate correction is scored on those same
    /// fragments. Re-cutting the core per candidate let a model that pushed
    /// awkward fragments outside its own core score better for having lost them.
    ///
    /// Returns NaN -- never 0.0 -- when the sample cannot support the statistic
    /// (under 40 residuals, a degenerate scale, or fewer than 4 populated bins).
    /// Zero is the best possible score, so using it as a failure value reports
    /// failure as perfection.
    static double systematicResidualPpm(const std::vector<MassResidual>& residuals,
                                        const std::function<double(double)>& correction);

    static Model fit(const std::vector<MassResidual>& residuals, const Options& options,
                     Diagnostics* diagnostics = nullptr);

    /// collect() then fit().
    static Model calibrate(const Library& library, SpectrumSource& source,
                           const Options& options, Diagnostics* diagnostics = nullptr);

    /// Multi-line human-readable report: residual count, fitted form and
    /// coefficients, sigma before and after, the window, and whether the gate
    /// passed. This is what the tool writes to the log.
    static std::string report(const Model& model, const Diagnostics* diagnostics = nullptr);

    // -- robust estimators, exposed because the synthetic test drives them
    //    directly and because they are the part most easily broken silently.

    /// Half-sample mode (Bickel & Fruhwirth 2006). @p sorted must be ascending.
    static double halfSampleMode(const std::vector<double>& sorted);

    /// Re-centre a starting location onto the peak it is near.
    ///
    /// `halfSampleMode` is exact on the case it was ported for -- a narrow peak
    /// on a FLAT background -- and drifts on the case ODIA actually has, a broad
    /// peak on a background that slopes. It keeps the shortest half-interval,
    /// and when the peak is wide relative to the window the shortest half can
    /// sit off-centre on the steeper flank. Measured on S08: -7.64 ppm where the
    /// peak is at -10.3, and that 2.6 ppm error was enough on its own to fail
    /// the gate, because the gate's centre band then straddles the background
    /// near zero instead of the peak.
    ///
    /// So the mode is used as a STARTING POINT and refined: take the median of
    /// everything within two robust sigmas, recompute, repeat. Where the mode
    /// was already right this changes nothing, which is why it is safe to apply
    /// unconditionally.
    static double refineLocation(const std::vector<double>& sorted, double start, double window);

    /// Robust scale about a KNOWN mode, estimated locally. Not the shorth: the
    /// shorth has 50% breakdown and blows up once the signal is a minority.
    static double localScaleAboutMode(const std::vector<double>& values, double mode,
                                      double initial_window);

    /// Density near zero against density at the edge of @p window. ~1 for a
    /// uniform sample, >> 1 for a genuine error distribution.
    static double peakednessRatio(const std::vector<double>& absolute_deviation, double window);

    /// Robust scale of the SIGNAL alone, with the flat background subtracted.
    ///
    /// `localScaleAboutMode` is the reference's estimator and is kept above
    /// unchanged, but it cannot be used at ODIA's purity. Its neighbourhood
    /// shrinks by 3 sigma per pass and stops when it stops shrinking, so when
    /// the first MAD is already dominated by background -- 12.5 ppm on a 2 ppm
    /// signal at 33% purity, giving 3 sigma = 55 ppm against a 50 ppm search
    /// window -- pass one cannot contract and it returns the MIXTURE's scale.
    /// Measured: 13.2 ppm for a signal whose true sigma is 2.0. Sizing a window
    /// from that is exactly the "wider answer means the fit failed" failure the
    /// reference warns about, arrived at from the other direction.
    ///
    /// So the background is subtracted instead of shrunk away from. Its density
    /// is read off the same edge band `peakednessRatio` uses, the excess is
    /// accumulated inwards, and the half-excess point is the signal's MAD.
    /// Exact on the synthetic case (2.0 ppm recovered from 2.0 ppm at 33%
    /// purity) and, unlike a contraction rule, it has no way to converge onto a
    /// dense sub-cluster of the background.
    ///
    /// @param sorted_absolute_deviation |residual - mode|, ascending.
    /// @param window the search half-width the sample was drawn from.
    static double backgroundCorrectedScale(const std::vector<double>& sorted_absolute_deviation,
                                           double window);
  };

} // namespace ODIA
