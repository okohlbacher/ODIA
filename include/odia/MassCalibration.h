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
// 3,000 sampled precursors, 160 cycles, 3,840 spectra decoded, 84 s:
//
//   gate            peakedness 5.10 against a control at 2.69, threshold 3
//   residuals       7,010 target and 1,602 control, after the brightness cut
//   correction      CONSTANT, -9.96 ppm
//   sigma           5.11 ppm, unchanged by the shape
//   m/z trend       real and refused: -12.90 ppm at 238 Th to -8.16 at 1,089,
//                   +2.62 ppm per e-fold at t=6.5, and it moves the residual
//                   only 5.11 -> 4.56 ppm (ratio 0.89, threshold 0.85)
//   rt drift        -0.89 ppm across 1,590 s at t=1.3, i.e. none
//   window          3 sigma would be 15.3 ppm, which is WIDER than the 10 ppm
//                   in force, so it is rejected -- calibration may only narrow
//
// Two independent checks that this is the instrument and not the method. An
// unrelated measurement over 124 M peak-transition hits at DIA-NN's retention
// times puts the offset at -9.78 ppm, the log-m/z fit at +2.79 ppm per e-fold,
// and finds no RT trend; all three agree. And sweeping the brightness cut from
// 0 to 0.9 walks the answer only from -10.51 to -9.50 ppm, so the selection is
// buying purity rather than manufacturing a number.
//
#pragma once

#include <odia/Library.h>
#include <odia/SpectrumSource.h>

#include <cmath>
#include <cstddef>
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

      /// Acquisition cycles probed, spread evenly over the gradient.
      ///
      /// One cycle is a contiguous run of spectra covering every isolation
      /// window, so this is also how the sample is spread in retention time --
      /// which is what makes an RT drift measurable at all.
      std::size_t cycles = 160;

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
      double min_peakedness = 3.0;

      /// Window = k x robust sigma of the CORRECTED residuals.
      double sigma_multiple = 3.0;

      /// Never infer a window below this, whatever the fit says.
      double floor_ppm = 2.0;

      /// Half-width of the ion-mobility test against the precursor's library
      /// 1/K0. 0 disables it; it is also skipped per precursor when the library
      /// carries no 1/K0, because absent information is not evidence of
      /// mismatch.
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

      /// An m/z-dependent term is accepted only when it clears BOTH tests.
      ///
      /// `min_slope_t`: the slope must be this many standard errors from zero.
      /// Necessary, and nowhere near sufficient -- with tens of thousands of
      /// residuals a slope can be overwhelmingly significant and still worth
      /// nothing.
      ///
      /// `max_sigma_ratio`: the residual spread AFTER the shaped correction must
      /// be at most this fraction of the spread after a plain constant. This is
      /// the test that actually decides, because it asks the only question that
      /// matters -- does modelling the shape leave a materially tighter
      /// distribution to size a window from?
      ///
      /// On S08 the m/z trend is real and measurable (-12.1 ppm at 200-288 Th
      /// rising to -5.8 ppm at 1067-1694 Th, ppm = -24.98 + 2.79 ln(m/z),
      /// measured over 124 M peak-transition hits) and STILL fails this test:
      /// it takes the residual MAD-SD from 7.82 to 7.62 ppm, a 2.6% improvement.
      /// A 2.6% gain does not justify a model that can contort at the ends of
      /// the m/z range, where a tryptic library has fewest fragments. So on this
      /// instrument the constant is the right answer and the gate says so.
      /// 0.85 rather than something closer to 1: a shape has to earn a clear
      /// reduction, not win a coin toss. S08 lands at 0.91 on this sampling and
      /// at 0.974 in the 124 M-hit measurement, so it is refused either way --
      /// but a threshold set where one of those two would have flipped it would
      /// be a threshold set by the noise in the estimate.
      double min_slope_t = 3.0;
      double max_sigma_ratio = 0.85;

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
      std::string form = "none";        ///< "none" | "constant" | "log_mz"
      std::string reason;               ///< why, in words, for the log

      double intercept_ppm = 0.0;       ///< correction at `reference_mz`

      /// ppm per e-fold in m/z, about `reference_mz`. 0 for the constant form.
      ///
      /// LOG, not linear in m/z, because that is the shape the error was
      /// measured to have: ppm = -24.98 + 2.79 ln(m/z) over 200-1694 Th on S08.
      /// A log basis is also the better-behaved of the two at the ends of the
      /// range, which is where the anchors run out and where an over-flexible
      /// model does its damage.
      double log_slope_ppm = 0.0;
      double reference_mz = 700.0;

      std::size_t residuals = 0;        ///< target-cell residuals used
      std::size_t decoy_residuals = 0;

      double sigma_before = 0.0;        ///< robust sigma of raw residuals, ppm
      double sigma_after = 0.0;         ///< ... after the chosen correction
      /// The two candidates, both reported, so the model choice is auditable
      /// rather than announced. `sigma_shaped` is 0 when no shape could be fitted.
      double sigma_constant = 0.0;
      double sigma_shaped = 0.0;
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
        if (log_slope_ppm == 0.0 || !(mz > 0.0) || !(reference_mz > 0.0)) { return intercept_ppm; }
        return intercept_ppm + log_slope_ppm * std::log(mz / reference_mz);
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
