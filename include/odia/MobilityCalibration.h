// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Per-run ion-mobility recalibration: the 1/K0 axis, calibrated from the run
// alone, the way MassCalibration already calibrates the fragment mass axis.
//
// ---------------------------------------------------------------------------
// WHY THE 1/K0 AXIS NEEDS THIS AT ALL
// ---------------------------------------------------------------------------
//
// ODIA recalibrates two of its three extraction axes per run -- fragment mass
// (MassCalibration) and retention time (RtCalibration) -- and left the third
// alone. The library's 1/K0 is not measured at all: it is a predicted CCS
// pushed through Mason-Schamp with ONE coefficient for every run, every charge
// and every m/z (1037.1902 on this project's timsTOF libraries). Nothing in
// that chain adapts it to the instrument that acquired the run, and the
// extractor then cuts a +/-0.025 window around it.
//
// Measured on S08_diaPASEF, iteration 3: of the precursors DIA-NN identifies
// and ODIA misses, 388 are missed because the signal lies OUTSIDE the
// extraction cell, and that bucket is a 1/K0 PREDICTION problem rather than a
// width problem. Two measurements say which:
//
//   * an oracle that replaces the library 1/K0 with DIA-NN's OBSERVED 1/K0,
//     keeping the window at +/-0.025, recovers 194 of those 388 against 61 at
//     baseline -- at ZERO cost in false positives, +7.99 recovery points;
//   * widening the window to +/-0.030 instead is worth only +0.56 points,
//     because a wider window admits interference about as fast as it admits
//     signal.
//
// So the centre is wrong, not the width -- the same shape of problem, and the
// same answer, as the mass axis. A prototype fitted per-charge, m/z-smooth
// offsets on the precursors ODIA already recovers on its OWN score and applied
// them to the ones it misses: 26% of the mean squared 1/K0 error removed,
// in-cell membership for the missed precursors 74.5% -> 81.0%, +2.21 recovery
// points end to end. This class is that, built the way the mass calibration is
// built.
//
// ---------------------------------------------------------------------------
// WHAT IS DIFFERENT FROM THE MASS AXIS, AND WHY THE CODE LOOKS DIFFERENT
// ---------------------------------------------------------------------------
//
// (1) THE RESIDUAL IS PER PRECURSOR, NOT PER FRAGMENT. Mass error is a property
//     of each fragment's own m/z, so MassCalibration collects one residual per
//     matched peak and gets thousands from a few hundred precursors. 1/K0 is a
//     property of the PRECURSOR -- all of its fragments travel through the TIMS
//     tunnel as the same ion -- so a precursor yields exactly one residual no
//     matter how many fragments confirm it. The sample is therefore an order of
//     magnitude smaller for the same probe, which is why `max_precursors`
//     defaults far higher here and why the model is deliberately coarse.
//
// (2) THE FRAGMENTS' AGREEMENT IS THE PURITY FILTER. MassCalibration buys purity
//     with, among other things, the precursor's library 1/K0 (`im_window`
//     0.010). This class cannot: that window is centred on the very number
//     being measured, and testing a peak against it would guarantee the answer
//     zero. What replaces it is that the fragments must agree with EACH OTHER:
//     a cell counts only when `min_fragments_matched` of one precursor's
//     fragments are found inside one `cluster_im`-wide slice of the mobility
//     axis in the same spectrum. That is a strong test -- a diaPASEF frame's
//     band is ~0.4 wide and the slice is 0.02 -- and it uses no prior about
//     where the precursor should be.
//
// (3) THE SEARCH IS SYMMETRIC ABOUT THE LIBRARY VALUE, AND MUST FIT INSIDE THE
//     BAND. The residual is (observed - predicted), so anything that truncates
//     it asymmetrically biases the location. The frame's mobility band is a
//     fixed interval that does NOT move with the precursor, so a precursor
//     sitting near a band edge can only produce residuals on one side. Those
//     precursors are excluded rather than corrected: `search_im` is a symmetric
//     half-width about the library 1/K0 and the whole of it has to lie inside
//     the band. Symmetric truncation costs residuals; asymmetric truncation
//     costs the answer.
//
// (4) THE MODEL IS PER CHARGE. Mobility depends on charge through Mason-Schamp
//     directly, and a single global coefficient absorbs that only if the
//     coefficient is right. Fitting one offset across charges would average two
//     populations that have no reason to share an error. A charge with too few
//     anchors gets NO correction rather than another charge's -- see
//     `min_anchors_per_charge`.
//
// (5) IT IS FITTED OUT OF SAMPLE. MassCalibration fits two parameters to
//     thousands of residuals and applies them to everything, and the
//     self-influence of one residual is negligible. Here the anchors ARE
//     precursors, the thing being corrected IS a precursor, and the selection
//     that makes something an anchor is "we found a bright coherent cluster for
//     it" -- so a precursor that entered its own correction has partly
//     memorised its own noise. Anchors are therefore split into `folds` by a
//     deterministic hash of their library index, a model is fitted for each fold
//     from the OTHER folds, and a precursor is corrected by the model that never
//     saw it. Precursors that were never anchors -- which is every precursor the
//     run misses, i.e. all of the ones this exists for -- are corrected by the
//     model fitted on all anchors, and were out of sample already.
//
// ---------------------------------------------------------------------------
// WHAT THE GATE DOES, AND ON WHAT
// ---------------------------------------------------------------------------
//
// The same three-part refusal as MassCalibration, plus one this axis needs:
//
//   NO MOBILITY AT ALL. A run whose spectra carry no 1/K0 -- 12_80 is SCIEX
//   SWATH and has none -- has no axis to calibrate. It is reported as exactly
//   that and the extraction is left untouched; it is NOT reported as a failed
//   fit, because "there is nothing here to measure" and "we measured and it was
//   noise" are different facts and only one of them is a warning. The same
//   holds for a library that carries no 1/K0.
//
//   TOO FEW RESIDUALS. Below `min_residuals` nothing is fitted.
//
//   FLAT RESIDUALS. The peakedness statistic ported from the reference, read on
//   the delta distribution over `gate_im`. Uniform (all-noise) deltas give ~1.
//
//   A CONTROL AS PEAKED AS THE DATA. Every precursor is probed again with all
//   of its fragments shifted in m/z (`decoy_shifts`). Those cells pass the
//   identical cluster, co-occurrence and apex tests and differ only in that
//   they cannot hold the real ions, so their delta distribution IS the null.
//
// And, per charge and per fold, the shape has to earn itself against a constant
// the same way the mass axis's does -- see `min_shape_chi2` and
// `min_shape_swing_im`.
//
// ---------------------------------------------------------------------------
// WHAT IT DOES NOT DO
// ---------------------------------------------------------------------------
//
//   * It does not resize the extraction window. `window_im` is fitted and
//     reported so the run says what width its corrected residual would support,
//     and it is NOT applied. The reason is measured, not stylistic: on S08 the
//     width lever is worth +0.56 points and the centring lever +2.21, and the
//     two are separate decisions. Applying an untested width change alongside a
//     tested centring change would make the measurement unattributable.
//
//   * It does not touch retention time or mass. Those have their own classes.
//
// ---------------------------------------------------------------------------
// WHAT IT COSTS, AND WHAT IT HAS SO FAR BOUGHT
// ---------------------------------------------------------------------------
//
// The BLIND probe on S08_diaPASEF, whole library, 200 cycles as 40 blocks of 5:
// 4,800 of 32,210 spectra decoded, 128 s against a 12.7 min extraction.
//
// The ANCHORED probe on the same run, 3,068 anchors (1,534 confident targets
// and a rank-matched null of the best-scoring decoys): 187 cycle blocks,
// 22,440 of 32,210 spectra, 446-453 s. That is the honest price of option (b)
// -- anchors are spread over the whole gradient, so "only the blocks that hold
// an apex" is 70% of the run, and the stage roughly doubles a single-pass
// extraction. It buys, out of fold:
//
//   * a gate that PASSES on its own margin: peakedness 12.41 against a null of
//     3.75, a margin of 3.31x where 1.25x is required. The blind probe managed
//     1.09x on the same run with the same gate.
//   * 23.7% of the mean squared 1/K0 error removed, scatter 0.0123 -> 0.0111;
//   * a curve that matches a truth it never sees. Against DIA-NN's OBSERVED
//     1/K0 on the 2,665 real targets: charge 3 constant +0.0065 against +0.0047
//     (the blind probe put it at +0.030), and the charge-2 m/z shape running
//     +0.0153 at 365 Th to -0.0144 at 1191 Th against a truth of +0.0163 to
//     -0.0199;
//   * end to end, on the frozen entrapment-calibrated discriminant: the 388
//     precursors iteration 3 measured as missed because the signal lies OUTSIDE
//     the extraction cell go from 15.7% to 23.2% recovered (+7.5 points, and
//     +8.3 on the 339 whose sub-reason is 1/K0), while every other miss bucket
//     moves by less than a point. Overall recovery 52.27% -> 53.17% at a
//     MATCHED, re-measured 1.00% entrapment false rate; 54.15% at the frozen
//     threshold, where the null itself has moved to 1.23% -- moving the
//     extraction cell moves the null, so the threshold has to be re-derived and
//     not inherited.
//
// 12_80, which has no mobility, stops after the first block in either mode:
// 82.3 s with 500 anchors read against 83.3 s with the stage switched off, i.e.
// it costs nothing to establish that there is nothing here.
//
// On S08 the BLIND probe bought a refusal. The gate read a peakedness of 3.99
// against a control at 3.65, a margin of 1.09x where 1.25x is required, and no
// offset was applied. That is the correct answer for that probe on that run:
// the lever is real (an oracle 1/K0 is worth +7.99 recovery points) but
// reaching it needs anchors the run has SCORED, not anchors a blind probe has
// guessed.
//
// ---------------------------------------------------------------------------
// THE TWO ANCHOR SOURCES, AND WHY THERE ARE TWO
// ---------------------------------------------------------------------------
//
// `collect()` is the BLIND probe: sample precursors, sample cycle blocks
// spread over the gradient, keep whatever cluster the fragments agree on. It
// needs nothing but the run, so it can be measured before the first extraction
// -- and that is also its defect. On the mass axis a cell that holds no
// precursor is harmless, because its residual is uniform in ppm and the mode
// steps over it. On the mobility axis it is not: a wrong cell's 1/K0 sits
// wherever the frame's peaks are DENSE, and peak density is a systematic.
// Three measurements on S08 say the blind probe was finding density and not
// precursors -- 7,995 Arabidopsis entrapment precursors, which cannot be in a
// human sample, gave the same residual distribution as the 2,665 real targets;
// the library's own decoys gave peakedness 3.65 against the targets' 3.99; and
// the fitted correction made the out-of-fold scatter WORSE, 0.0099 -> 0.0141.
//
// `collectAt()` is the ANCHORED probe, and it is the one the prototype's +2.21
// points came from. It is given a list of precursors that have ALREADY been
// identified confidently, each with the retention time of the peak group that
// identified it, and it looks only there -- one sequential walk of the run,
// visiting only the cycle blocks that hold an anchor's apex. Two things change
// and both matter:
//
//   * WHERE it looks is no longer a guess. A blind probe asked to find a
//     precursor somewhere in a 150 s window over 40 scattered blocks gets one
//     draw from the interference per block whether the precursor is there or
//     not. An anchored probe looks at the few cycles where a SCORED peak group
//     already put the peak, so an empty cell has nowhere to hide.
//
//   * WHAT the null is changes with it. The m/z-shifted control is too easy
//     here -- MassCalibration's +7.33 Th shift moves a fragment off the
//     amino-acid mass lattice into a part of the spectrum where peaks never
//     are, and it reported 99.7% purity on a sample that was almost entirely
//     noise. The anchored probe's null is the LIBRARY'S OWN DECOYS, probed at
//     the apexes their own peak groups claimed. A decoy's fragments are real
//     fragment masses of a real (shuffled) sequence and its claimed apex sits
//     on real signal, so the null is drawn from the same interference the
//     targets swim in.
//
// The chain is: pass 1 extracts wide and scores -> the confident targets and
// the best-scoring decoys become anchors -> `collectAt()` measures their 1/K0
// -> `fit()` gates and fits -> pass 2 extracts on the corrected axis. The fold
// machinery is unchanged and still load-bearing: an anchor never receives the
// correction its own fold produced.
//
// What the anchored probe COSTS is a second sequential pass over the run's
// spectra. Anchors are spread over the whole gradient by construction, so the
// blocks they select cover most of it; the alternative considered and rejected
// was to accumulate an intensity-weighted 1/K0 alongside every chromatogram
// point, which needs no extra decode but roughly doubles the chromatogram
// store (1.65 -> 3.3 GiB on the S08 combined library) and pays that cost for
// every point when only the apex is ever read.
//
#pragma once

#include <odia/Library.h>
#include <odia/SpectrumSource.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace ODIA
{

  /// One precursor's measured 1/K0 error: what the library predicted, what the
  /// run showed, and everything the fit needs to place it.
  ///
  /// Deliberately not pre-aggregated, for the same reason MassResidual is not:
  /// the question is whether the error depends on charge and on m/z, and a
  /// histogram collapsed over either cannot answer it.
  struct MobilityResidual
  {
    std::uint32_t precursor = 0;   ///< library index, which is also the fold key
    float mz = 0.0f;               ///< precursor m/z, Th
    std::uint8_t charge = 0;
    float im_library = 0.0f;       ///< the library's predicted 1/K0
    float im_observed = 0.0f;      ///< intensity-weighted 1/K0 of the winning cluster
    float delta = 0.0f;            ///< observed - library
    float intensity = 0.0f;        ///< summed intensity of the cluster
    float rt = 0.0f;               ///< retention time of the winning spectrum, s
    std::uint8_t fragments = 0;    ///< how many agreed
    std::uint8_t cycles = 0;       ///< consecutive cycles it was seen in
    bool decoy = false;            ///< a control cell -- shifted, or a library decoy
  };

  /// One place in the run where a peak group has already been scored, and the
  /// precursor it was scored for.
  ///
  /// This is what turns the blind probe into the anchored one. `rt` is the
  /// APEX of a peak group, measured, not predicted -- so the probe does not
  /// search a window, it visits a place. `decoy` marks the anchors that make up
  /// the null: the same procedure applied to the library's own decoys, whose
  /// best peak groups sit on real signal that is nevertheless not theirs.
  struct MobilityAnchor
  {
    std::uint32_t precursor = 0;   ///< library index
    float rt = 0.0f;               ///< the peak group's apex retention time, s
    bool decoy = false;            ///< true for the null arm
  };

  /// Fits the run's 1/K0 prediction error, per charge and as a function of m/z.
  class MobilityCalibration
  {
  public:
    /// Highest charge given its own curve. Above this a precursor is corrected
    /// by the top slot, which on a tryptic library holds nothing anyway.
    static constexpr std::size_t MAX_CHARGE = 8;

    struct Options
    {
      /// Half-width of the 1/K0 search about the library value, and the margin a
      /// precursor's library 1/K0 must keep from BOTH edges of the frame's
      /// mobility band before it is used at all.
      ///
      /// Symmetric on purpose -- see (3) in the header. 0.06 is ~6x the observed
      /// scatter and ~1.5x the p99 of the prediction error, so it is wide enough
      /// that the distribution's shoulders are visible, while still leaving room
      /// inside a ~0.4-wide diaPASEF band for most of the library.
      double search_im = 0.06;

      /// Half-width the LOCATION and the GATE are evaluated over.
      ///
      /// Distinct from `search_im` for the reason MassCalibration's `gate_ppm`
      /// is distinct from its `search_ppm`: the search must be wide or the
      /// distribution is truncated and its scale censored, while the gate must
      /// be narrow, because a flat background contributes to its edge band in
      /// proportion to the width searched and so dilutes the very contrast the
      /// gate exists to measure.
      double gate_im = 0.05;

      /// Half-width of the mobility slice a precursor's fragments must share
      /// before they count as one detection, and over which its 1/K0 is
      /// centroided.
      ///
      /// This is the purity filter that replaces the prior this class is not
      /// allowed to use (see (2)). It is a property of the TIMS peak, not of the
      /// prediction error: one ion's fragments arrive together to within the
      /// mobility peak width, which on this instrument class is ~0.01-0.02 in
      /// 1/K0. Deliberately NOT the extractor's 0.025, which is sized to survive
      /// a prediction error rather than to describe a peak.
      double cluster_im = 0.010;

      /// Keep only cells whose summed matched intensity is at or above this
      /// quantile of all target cells. 0 keeps everything.
      ///
      /// The same argument as MassCalibration's `min_intensity_quantile`: the
      /// interferent population is numerous but individually weak, so brightness
      /// is evidence. Applied to the control on the same threshold, so the null
      /// keeps measuring the same procedure.
      double min_intensity_quantile = 0.25;

      /// Precursors sampled from the library, spread evenly through it.
      ///
      /// 0 means the whole library, which is the default, because one precursor
      /// yields ONE residual here where it yields a dozen on the mass axis --
      /// see (1) -- and because the probe's cost is the spectra it decodes,
      /// which `cycles` sets, not the queries it evaluates against them.
      std::size_t max_precursors = 0;

      /// Acquisition cycles probed. This is what the measurement costs: one
      /// cycle is one decoded spectrum per isolation window.
      std::size_t cycles = 200;

      /// ...drawn as `cycles / cycle_block` CONTIGUOUS BLOCKS of this many
      /// cycles, and a cell counts only when the same precursor is found at the
      /// same mobility in `min_cycles_matched` of one block's cycles.
      ///
      /// THIS IS WHAT MAKES THE PROBE WORK AT ALL, and it is the one place the
      /// design departs from MassCalibration rather than following it. That
      /// class draws single scattered cycles and keeps each precursor's
      /// BRIGHTEST cell, which for an absent precursor is the luckiest noise --
      /// harmless there, because a wrong cell's MASS residual is uniform and the
      /// mode steps over it. It is not harmless here, because a wrong cell's
      /// MOBILITY residual is not uniform: it sits wherever the frame's peaks
      /// are dense, which is a systematic, and the mode walks straight into it.
      ///
      /// Measured on S08 with scattered cycles and no reproducibility rule, and
      /// measured against the entrapment library rather than against an
      /// m/z-shifted control -- Arabidopsis precursors that are absent from a
      /// human sample by construction:
      ///
      ///     human targets   charge 2 median +0.0067, charge 3 +0.0254
      ///     ENTRAPMENTS     charge 2 median +0.0035, charge 3 +0.0358
      ///
      /// Precursors that cannot be there gave the same answer as precursors
      /// that are, including the same large charge-3 offset. The sample was
      /// noise, and an m/z-shifted control had said it was 99.7% pure -- because
      /// shifting a fragment by 7.33 Th moves it off the peptide mass-defect
      /// line, where no peaks live at all, so that control is not a null for
      /// this axis. Both facts are why the rule below exists and why the header
      /// says what the control is worth.
      ///
      /// Consecutive cycles are ~1.4 s apart on this run and a chromatographic
      /// peak is tens of seconds wide, so a real precursor is in every cycle of
      /// a block that overlaps its elution and a coincidence is in one.
      std::size_t cycle_block = 5;
      std::size_t min_cycles_matched = 4;

      /// How far from an ANCHOR's apex the anchored probe may look, in seconds.
      ///
      /// Only `collectAt()` reads this; the blind probe uses
      /// `rt_window_seconds` around a PREDICTED retention time and needs it
      /// wide, because a predicted iRT is wrong by more than a peak is broad.
      /// An anchor's retention time is measured, so this only has to cover the
      /// blocks the peak actually spans. 0 means "derive it from the run": the
      /// median cycle duration times `cycle_block`, which is the width of the
      /// block the probe scores over.
      double anchor_rt_window_seconds = 0.0;

      /// Seed for that draw. Fixed, not clock-derived.
      std::uint64_t sample_seed = 0x0D1A1130u;

      /// Fragments used per precursor, most intense first in library order.
      std::size_t max_fragments = 12;

      /// Fragments of one precursor that must land in the SAME mobility slice of
      /// the SAME spectrum before the cell counts.
      ///
      /// Higher than MassCalibration's 3 because here it is doing the whole job
      /// of the purity filter rather than being the cheapest of four.
      std::size_t min_fragments_matched = 4;

      /// The fragment mass window the probe matches through, ppm, and the
      /// systematic correction to apply to it. Normally handed the values
      /// MassCalibration just fitted: this probe runs after it, and matching
      /// through a mis-centred mass window would fill the sample with the
      /// interference the mass calibration exists to exclude.
      double fragment_ppm = 15.0;
      double fragment_ppm_offset = 0.0;
      double fragment_ppm_log_slope = 0.0;
      double fragment_ppm_slope_per_1000 = 0.0;
      double fragment_ppm_ref_mz = 700.0;

      /// How the null is built.
      ///
      /// LIBRARY_DECOYS is the default, and the m/z-shifted control that
      /// MassCalibration uses is kept only as a fallback for a library that has
      /// none, because on this axis it was MEASURED to be invalid.
      ///
      /// The measurement: with an m/z-shifted control the S08 probe reported
      /// 99.7% purity, and the entrapment library -- 7,995 Arabidopsis
      /// precursors that cannot be in a human sample -- reported the same
      /// residual distribution as the real human targets, including the same
      /// charge-3 offset. Both cannot be true. The reason the shift lies is that
      /// a peptide fragment's m/z is not a free number: b and y ions lie on the
      /// amino-acid mass lattice, so shifting one by 7.33 Th moves it into a
      /// region of the spectrum where peaks essentially never are, and the
      /// control then measures the absence of peaks rather than the absence of
      /// THIS peptide.
      ///
      /// The library's own decoys do not have that defect. They are shuffled or
      /// mutated peptides, so their fragment masses sit on the same lattice as
      /// any real peptide's and are transmitted in the same windows, and they
      /// carry a predicted 1/K0 produced by the same predictor -- they differ
      /// from a target in exactly one way, which is that the peptide is not in
      /// the sample. That is what a null has to be.
      enum class Control { LibraryDecoys, MzShift, None };
      Control control = Control::LibraryDecoys;

      /// Restrict the probe to blocks near where the library says the precursor
      /// elutes: rt = irt_slope * iRT + irt_intercept, +/- rt_window_seconds.
      /// 0 for either means no restriction, which is the default because a run
      /// need not have an iRT map before its first pass.
      ///
      /// It is worth having when there IS one. The probe keeps the brightest
      /// cluster a precursor produces over every block it probes, so a precursor
      /// that is ABSENT gets one draw per block from the interference and the
      /// maximum of forty draws is a long way into that distribution, while one
      /// that is PRESENT has its real peak in one or two blocks only. Cutting
      /// the blocks a precursor may be found in cuts the absent precursor's
      /// advantage and not the present one's.
      double irt_slope = 0.0;
      double irt_intercept = 0.0;
      double rt_window_seconds = 0.0;

      /// Used only when `control` is MzShift, or when it is LibraryDecoys and
      /// the library has none. Each shifts ALL of a precursor's fragments
      /// together, so the cell passes the same cluster and co-occurrence tests.
      std::vector<double> decoy_shifts{7.33, -7.19};

      /// Below this many target residuals nothing is fitted.
      int min_residuals = 200;

      /// Reject the fit unless the density of deltas near the mode exceeds the
      /// density at the edge of `gate_im` by this factor. Uniform deltas give
      /// ~1. Threshold ported unchanged from MassCalibration.
      double min_peakedness = 3.0;

      /// And reject it unless the data is this many times more peaked than the
      /// NULL. MassCalibration only asks that the control be less peaked at all;
      /// that is too weak here, and S08 is why.
      ///
      /// With the library's own decoys as the control, the S08 probe measures a
      /// peakedness of 3.02 against a control at 3.00. Both clear the absolute
      /// threshold, and the bare "control must be lower" rule passes the run by
      /// 0.02 -- on a sample where an independent null (7,995 Arabidopsis
      /// entrapment precursors, which cannot be in a human sample) gives the
      /// same residual distribution as the real targets, and where the fitted
      /// correction makes the out-of-fold scatter WORSE, 0.0099 -> 0.0141.
      /// A ratio of 1.25 refuses that and would still have accepted the mass
      /// axis's own passing case, which measures 4.34 against 2.98.
      double min_control_margin = 1.25;

      /// Refuse the fitted correction when the fraction of mean squared 1/K0
      /// error it removes OUT OF FOLD is at or below this.
      ///
      /// `squared_error_removed` was computed out of fold over `folds` and then
      /// only PRINTED, several hundred lines after the gate had decided -- so
      /// the engine measured the one quantity that answers "is this correction
      /// worth applying" and did not consult it. Measured on S08, 500k
      /// precursors, every one of these reporting GATE PASSED on the peakedness
      /// heuristic alone:
      ///
      ///     arm          removed out of fold   robust scatter
      ///     base x2              67.1%         0.0118 -> 0.0113
      ///     base x1              56.5%         0.0085 -> 0.0090
      ///     affine x1             4.7%         0.0117 -> 0.0118
      ///     affine x2            -5.5%         0.0105 -> 0.0112
      ///
      /// 0.0 is a threshold and is named as one. It is defensible as a NULL --
      /// it compares the correction against doing nothing -- but it treats a
      /// point estimate of -0.1% and one of -50% identically, and both reviews
      /// of this change said the principled version is a fold-level interval:
      /// refuse when the upper bound is <= 0 if the policy is "apply only
      /// demonstrated benefit", or when the lower bound is > 0 is violated if it
      /// is "catch only demonstrated harm". That interval is computable from the
      /// `folds` already fitted and is the right next step; this is the
      /// provisional conservative version.
      ///
      /// Set very negative to disable -- which the fold-isolation test does,
      /// because it plants DIFFERENT offsets per fold on purpose and so
      /// constructs a model that is genuinely harmful out of fold. That test is
      /// asserting fold mechanics, not calibration quality.
      double min_squared_error_removed = 0.0;

      /// A charge gets its own curve only with at least this many anchors. Below
      /// it the charge is left UNCORRECTED -- not given the pooled offset, which
      /// would be another charge's answer applied to a population that has no
      /// reason to share it.
      std::size_t min_anchors_per_charge = 120;

      /// Anchors below which a charge gets a CONSTANT plus the POOLED slope,
      /// rather than nothing at all.
      ///
      /// Measured on S08/lib_targets: charge 3 brings 93 anchors against the
      /// 120 above, so it is left entirely uncorrected -- and charge 3's
      /// residual is the WORSE of the two (constant +0.0335, slope -0.113,
      /// against +0.0019 and -0.096 for charge 2). The charge that most needs
      /// the correction is the one that cannot reach the count for it.
      ///
      /// STRUCK 2026-08-28 (doc/36 section 5b, ordered removed 2026-08-18 and
      /// still live until today). The physical justification was FALSE:
      ///
      ///   ~~"Pooling the SLOPE is physically justified ... The slope is a
      ///   relative scale error in the CCS->1/K0 conversion coefficient, and
      ///   that coefficient is a property of the conversion, shared by every
      ///   charge."~~
      ///
      /// The fitted slope is `c*lambda_z`, where lambda_z is the CHARGE'S OWN
      /// reliability ratio. Measured, `a_q` spans -0.058 to -0.078, a 28%
      /// spread across charge, so the pooled quantity is demonstrably not
      /// shared. Pooling survives only as statistical regularisation of a small
      /// sample -- which is a weaker claim and does not by itself justify
      /// borrowing across charges. It also explains this tier's signature,
      /// out-of-fold MSE improving while identifications fall: it applies one
      /// charge's shrinkage fraction to another.
      ///
      /// The m/z shape is NOT borrowed either -- it needs bins, and a charge
      /// that cannot reach 120 anchors cannot fill 8 of them.
      ///
      /// DEFAULTED OFF, because the reasoning above is sound and the
      /// measurement still says no. S08/lib_targets, at 1% FDR:
      ///
      ///     off (0)    1232      MSE removed out of fold 14.1%, 1 of 2 charges
      ///     on  (40)   1212      MSE removed out of fold 18.8%, 2 of 2 charges
      ///
      /// It corrects charge 3, it improves the calibration's own out-of-fold
      /// error by a third, and it costs 20 identifications. That divergence is
      /// the point: **the gate metric and the outcome disagree**, for the second
      /// time in this file -- the m/z shape did the same, passing on out-of-fold
      /// MSE while the robust scatter got worse. Mean squared 1/K0 error is not
      /// a proxy for identifications and must not be used as one.
      ///
      /// The likely mechanism is dilution. Charge 2 brings 378 anchors against
      /// charge 3's 93, so the pooled slope is essentially charge 2's, and the
      /// run fitted it at -0.0423 where charge 3's own residual wants about
      /// -0.113. Applying a quarter of the needed correction, plus the noise of
      /// a shared estimate, is worse for charge 3 than leaving it alone -- so
      /// the shared-coefficient argument, which is physically right, is defeated
      /// by the estimator being too weak and too unequally weighted to realise
      /// it. An anchor-count-balanced or errors-in-variables fit might; a plain
      /// pooled least squares does not.
      ///
      /// Kept, off, because it is the right shape for the problem and needs a
      /// better estimator rather than a different idea. Set 40 to enable.
      std::size_t min_anchors_pooled_slope = 0;

      /// Equal-COUNT m/z bins per charge, and the floor on each one's
      /// occupancy. Equal count rather than equal width because a library's m/z
      /// distribution is very far from uniform.
      std::size_t max_mz_bins = 8;
      std::size_t min_per_bin = 40;

      /// Out-of-sample folds. A precursor that was an anchor is corrected by the
      /// model fitted from the other folds; one that was not is corrected by the
      /// model fitted from all of them.
      std::size_t folds = 4;

      /// The m/z shape has to earn itself against a per-charge constant, and
      /// these are the two conditions -- one relative, one absolute, exactly as
      /// MassCalibration pairs `max_systematic_ratio` with
      /// `min_systematic_gain_ppm`.
      ///
      /// `min_shape_chi2`: chi-square per degree of freedom of the bin locations
      /// about their own weighted mean. At 1 the bins are consistent with noise
      /// and there is no shape to model; the threshold is the assertion that the
      /// bins disagree by several times their own uncertainty.
      ///
      /// `min_shape_swing_im`: and they must disagree by an amount that MATTERS
      /// against the window the correction is applied through. Without an
      /// absolute floor, a large enough anchor set makes every bin's standard
      /// error small enough for chi-square alone to wave a shape through.
      double min_shape_chi2 = 3.0;
      double min_shape_swing_im = 0.004;

      /// No correction larger than this is ever applied, whatever the fit says.
      ///
      /// This is the mobility axis's version of "calibration may only narrow".
      /// A correction bigger than the window it is applied through does not
      /// refine where the extractor looks, it moves it somewhere else entirely
      /// on the strength of a fit -- and a fit that says that is far more likely
      /// to have found the interference than the instrument.
      double max_correction_im = 0.030;

      /// Bound on the mobility-linear term, in 1/K0 per 1/K0.
      ///
      /// A relative scale error in the CCS->1/K0 conversion: a slope of s is an
      /// s*100% error in the coefficient. Measured at ~0.10 on S08. This is
      /// deliberately NOT max_correction_im -- that bounds an offset, and reusing
      /// it here clamped the real trend on both charge states.
      double max_im_slope = 0.25;

      /// Window half-width the corrected residual would support: k x robust
      /// sigma, floored. FITTED AND REPORTED, NEVER APPLIED -- see the header.
      double sigma_multiple = 3.0;
      double floor_im = 0.005;

      bool verbose = false;
    };

    /// One charge's correction for one fold: either a constant or a piecewise
    /// linear function of m/z through the bin locations.
    ///
    /// Piecewise linear through measured bins rather than a fitted basis, which
    /// is the one place this deliberately diverges from MassCalibration. There
    /// the error is a smooth instrument property with a physical shape (a TOF's
    /// flight-time-to-mass relation), so a two-parameter basis both fits and
    /// says something. Here the residual is the error of a PREDICTOR --
    /// PeptDeep's CCS, pushed through one Mason-Schamp coefficient -- and a
    /// predictor's bias has no reason to be log-linear in m/z. So the shape is
    /// read off the data and interpolated, and is held flat outside the range
    /// where anchors were found so that it can never extrapolate.
    struct Curve
    {
      bool supported = false;      ///< false means "apply nothing to this charge"
      bool shaped = false;         ///< false means the constant won
      double constant = 0.0;
      std::vector<double> knot_mz;      ///< ascending
      std::vector<double> knot_offset;
      std::size_t anchors = 0;

      /// Linear term in 1/K0 itself: correction = constant + im_slope * (im - im_pivot).
      ///
      /// The constant and the m/z shape between them cannot express what the
      /// residuals actually do. Measured on S08 against DIA-NN's confident set,
      /// library 1/K0 against observed, 670 precursors:
      ///
      ///   library 1/K0   n    median residual   beyond +/-0.025
      ///     0.6-0.9    139       +0.0170             25.9%
      ///     0.9-1.0    189       +0.0061              8.5%
      ///     1.0-1.1    178       +0.0027              4.5%
      ///     1.1-1.2    114       -0.0079             12.3%
      ///     1.2-1.6     50       -0.0106             20.0%
      ///
      /// A monotone trend from +0.017 to -0.011, crossing zero near 1.05 -- a
      /// SLOPE, not an offset. It comes from upstream: ccs_to_mobility.py
      /// applies a single --coefficient, i.e. a pure proportionality, and the
      /// Mason-Schamp relation between CCS and 1/K0 is not proportional. So the
      /// conversion is right in the middle of the range and wrong at both ends.
      ///
      /// 12.5% of precursors land outside the +/-0.025 extraction window
      /// because of it, against a 14% picker loss -- close enough that this is
      /// a strong candidate for most of that loss.
      ///
      /// The m/z shape cannot substitute. CCS correlates with m/z, so an m/z
      /// spline absorbs SOME of this, which is exactly why it is dangerous: it
      /// fits a mobility error through a mass-shaped model and leaves the
      /// residual structured.
      double im_slope = 0.0;
      double im_pivot = 1.0;

      /// The correction including the mobility-linear term. `at(mz)` is kept
      /// for callers that have no 1/K0 to hand.
      double at(double mz, double im) const
      {
        if (!supported) { return 0.0; }
        const double base = at(mz);
        if (!std::isfinite(im) || im_slope == 0.0) { return base; }
        const double composite = base + im_slope * (im - im_pivot);
        // The COMPOSITE is what max_correction_im bounds, not each term.
        //
        // Every component is clamped on its own -- the constant, each knot, and
        // the slope against max_im_slope -- and before the slope existed that
        // was the same thing. It is not any more: a timsTOF library spans 1/K0
        // 0.6-1.6 about a pivot near 1.0, so a fitted -0.10 slope contributes
        // +/-0.05 at the ends, and a knot may add 0.030 on top. That is ~0.08
        // applied against a documented cap of 0.030 and an extraction half-
        // window of 0.025 -- the correction could move a precursor clean out of
        // the window it is meant to centre, and it would do it worst at the
        // extremes of the mobility range, which is the population this whole
        // class exists to rescue.
        return std::clamp(composite, -max_correction, max_correction);
      }

      /// The bound the composite is held to. Carried on the curve so `at()` can
      /// enforce it without reaching for Options, which it does not have.
      double max_correction = 0.030;

      double at(double mz) const
      {
        if (!supported) { return 0.0; }
        if (!shaped || knot_mz.size() < 2) { return constant; }
        if (!(mz > knot_mz.front())) { return knot_offset.front(); }
        if (!(mz < knot_mz.back())) { return knot_offset.back(); }
        const auto it = std::upper_bound(knot_mz.begin(), knot_mz.end(), mz);
        const std::size_t hi = static_cast<std::size_t>(it - knot_mz.begin());
        const std::size_t lo = hi - 1;
        const double span = knot_mz[hi] - knot_mz[lo];
        if (!(span > 0.0)) { return knot_offset[lo]; }
        const double t = (mz - knot_mz[lo]) / span;
        return knot_offset[lo] + t * (knot_offset[hi] - knot_offset[lo]);
      }
    };

    /// Per-bin summary, so the shape is reported as numbers rather than as a
    /// curve the caller has to take on trust.
    struct Bin
    {
      double centre = 0.0;      ///< mean m/z of the bin
      double location = 0.0;    ///< median delta, 1/K0
      double stderr_im = 0.0;
      std::size_t n = 0;
    };

    /// The fitted correction, plus everything needed to judge it.
    struct Model
    {
      bool fitted = false;
      /// True unless the RUN carries no 1/K0 at all, or the LIBRARY does. Those
      /// are reported separately from a failed fit because they are not failures:
      /// there is nothing on this axis to measure.
      bool run_has_mobility = true;
      bool library_has_mobility = true;

      /// "none" | "constant" | "mz_shaped"
      std::string form = "none";
      std::string reason;               ///< why, in words, for the log

      std::size_t residuals = 0;        ///< target cells that produced a delta
      std::size_t decoy_residuals = 0;
      double peakedness = 0.0;
      double decoy_peakedness = 0.0;
      double gate_peakedness = 0.0;     ///< the threshold it was judged against
      double gate_im = 0.0;
      double search_im = 0.0;

      double location = 0.0;            ///< robust centre of all target deltas
      double sigma_before = 0.0;        ///< robust scale of the raw deltas
      double sigma_after = 0.0;         ///< ... after the per-charge correction
      /// Fraction of the mean squared delta the correction removes, measured on
      /// the anchors OUT OF FOLD. The prototype's headline number.
      double squared_error_removed = 0.0;

      /// The same fraction computed within each fold. Its SPREAD is the only
      /// uncertainty available for `squared_error_removed`, and it matters: a
      /// ratio on a small denominator is unstable, and at a cell whose
      /// correction is near zero the pooled figure moved 6.4 points between two
      /// runs differing only by a 0.0004 ppm mass offset. A guard that reads the
      /// pooled point estimate alone flips on noise precisely where its decision
      /// is closest.
      std::vector<double> fold_error_removed;
      /// Window half-width the corrected residual supports. REPORTED, NEVER
      /// APPLIED -- see the header.
      double window_im = -1.0;

      /// Per charge, in charge order, for the report.
      struct ChargeReport
      {
        std::uint8_t charge = 0;
        std::size_t anchors = 0;
        bool supported = false;
        bool shaped = false;
        double constant = 0.0;
        double low = 0.0, high = 0.0;   ///< the correction at the m/z extremes
        double mz_low = 0.0, mz_high = 0.0;
        double chi2_per_dof = 0.0;
        double swing = 0.0;
        std::string reason;
        std::vector<Bin> bins;
      };
      std::vector<ChargeReport> by_charge;

      std::size_t folds = 1;
      /// The library indices that were anchors, ascending. Small -- bounded by
      /// `max_precursors` -- and carried so that `offsetFor` can tell an anchor
      /// (corrected out of fold) from everything else.
      std::vector<std::uint32_t> anchors;
      /// (MAX_CHARGE + 1) x (folds + 1) curves, charge-major. The last fold slot
      /// is the model fitted on ALL anchors, which is what a non-anchor gets.
      std::vector<Curve> curves;

      /// Which fold model corrects @p precursor: its own fold if it was an
      /// anchor, `folds` (the all-anchor model) if it was not.
      std::size_t foldOf(std::uint32_t precursor) const
      {
        if (folds == 0) { return 0; }
        if (!std::binary_search(anchors.begin(), anchors.end(), precursor)) { return folds; }
        return foldIndex(precursor, folds);
      }

      /// The 1/K0 to ADD to this precursor's library value. 0 when nothing was
      /// fitted, when this charge is unsupported, or when the precursor's
      /// mobility is unknown -- in every case, the uncorrected library value.
      /// @param library_im the precursor's LIBRARY 1/K0, which the linear term
      ///        is a function of. Defaults to NaN, which yields the constant
      ///        plus the m/z shape and nothing else -- so a caller that has no
      ///        mobility to hand gets exactly the pre-slope behaviour rather
      ///        than a silently wrong correction.
      double offsetFor(std::uint32_t precursor, double mz, int charge,
                       double library_im = std::numeric_limits<double>::quiet_NaN()) const
      {
        if (!fitted || curves.empty()) { return 0.0; }
        const std::size_t c = static_cast<std::size_t>(
          std::clamp<int>(charge, 0, static_cast<int>(MAX_CHARGE)));
        const std::size_t f = std::min(foldOf(precursor), folds);
        return curves[c * (folds + 1) + f].at(mz, library_im);
      }
    };

    /// Deterministic fold of a library index. Static and public because the
    /// out-of-sample test has to be able to ask which fold a precursor is in
    /// before it plants anything in it -- otherwise "the fit excluded it" is an
    /// assertion rather than a measurement.
    static std::size_t foldIndex(std::uint32_t precursor, std::size_t folds);

    struct Diagnostics
    {
      std::size_t blocks = 0;           ///< contiguous cycle blocks probed
      std::size_t precursors_sampled = 0;
      std::size_t precursors_with_library_im = 0;
      /// Excluded because `search_im` about their library 1/K0 did not fit
      /// inside their frame's mobility band -- see (3).
      std::size_t precursors_outside_band = 0;
      /// Cell-spectra skipped because the spectrum was outside the precursor's
      /// predicted retention-time window.
      std::size_t cells_off_rt = 0;
      std::size_t cells_probed = 0;
      std::size_t target_cells = 0;
      std::size_t decoy_cells = 0;
      std::size_t control_precursors = 0;   ///< library decoys probed as the null
      std::string control_kind = "none";

      /// Where the probed positions came from, for the log: the blind sample,
      /// or a named number of scored anchors.
      std::string anchor_kind = "a blind sample of the run";
      std::size_t anchors_given = 0;        ///< handed to collectAt(), before filtering
      std::size_t anchors_probed = 0;       ///< of those, actually turned into a cell
      std::size_t spectra_decoded = 0;
      std::size_t spectra_with_mobility = 0;
      double band_width = 0.0;          ///< median frame band width seen
      double collect_seconds = 0.0;
      bool run_has_mobility = true;
      bool library_has_mobility = true;
    };

    /// Probe the run and return the raw residuals, targets and controls.
    static std::vector<MobilityResidual> collect(const Library& library, SpectrumSource& source,
                                                 const Options& options,
                                                 Diagnostics* diagnostics = nullptr);

    /// The ANCHORED probe: measure 1/K0 only where @p anchors say a peak group
    /// was already scored.
    ///
    /// One sequential walk of the run, visiting the cycle blocks that hold an
    /// anchor's apex and no others. Anchors need not be sorted. Anchors marked
    /// `decoy`, and anchors that name a library decoy, produce control
    /// residuals; `Options::control` is not consulted, because the null on this
    /// axis is the library's own decoys and an m/z shift is measurably too easy.
    static std::vector<MobilityResidual> collectAt(const Library& library,
                                                   SpectrumSource& source,
                                                   const std::vector<MobilityAnchor>& anchors,
                                                   const Options& options,
                                                   Diagnostics* diagnostics = nullptr);

    /// Fit a model to residuals from anywhere -- a run, or a synthetic sample.
    static Model fit(const std::vector<MobilityResidual>& residuals, const Options& options,
                     Diagnostics* diagnostics = nullptr);

    /// collect() then fit().
    static Model calibrate(const Library& library, SpectrumSource& source,
                           const Options& options, Diagnostics* diagnostics = nullptr);

    /// collectAt() then fit(), with the same no-mobility reporting calibrate()
    /// has. This is the pass-2 entry point.
    static Model calibrateFrom(const Library& library, SpectrumSource& source,
                               const std::vector<MobilityAnchor>& anchors,
                               const Options& options, Diagnostics* diagnostics = nullptr);

    /// Multi-line human-readable report. This is what the tool writes to the log.
    static std::string report(const Model& model, const Diagnostics* diagnostics = nullptr);

  private:
    /// The one probe body behind collect() and collectAt(). @p anchors null is
    /// the blind sample; non-null is the anchored walk.
    static std::vector<MobilityResidual> probe(const Library& library, SpectrumSource& source,
                                               const Options& options,
                                               const std::vector<MobilityAnchor>* anchors,
                                               Diagnostics* diagnostics);
  };

} // namespace ODIA
