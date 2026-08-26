// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/ChromatogramExtractor.h>
#include <odia/MassCalibration.h>
#include <odia/Ms1Traces.h>
#include <odia/OpenSwathPicker.h>
#include <odia/Library.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace ODIA
{

  /// Finds candidate peak groups in extracted chromatograms, describes each
  /// with sub-scores, and hands the matrix to the semi-supervised scorer.
  ///
  /// The sub-scores are deliberately the subset that can be computed honestly
  /// from what extraction currently produces. Four of the standard OpenSWATH
  /// set are absent on purpose: `rt_delta` needs a real iRT calibration (the
  /// tool spreads the library evenly over the run when none is given),
  /// `mass_error` needs the observed m/z of each matched peak, which the
  /// extractor discards, `im_delta` needs the CCS -> 1/K0 conversion that is
  /// deliberately downstream, and `isotope_corr` needs MS1 extraction that
  /// does not exist. A sub-score computed from a placeholder is worse than an
  /// absent one: the classifier gives it weight and it is noise.
  class PeakGroupScorer
  {
  public:
    /// Order matters -- it is the column order of the feature matrix, and it is
    /// what `subScoreNames()` reports. Adding one means retraining.
    enum SubScore : std::size_t
    {
      XCORR_SHAPE = 0,     ///< mean all-pairs cross-correlation at the peak
      XCORR_COELUTION,     ///< mean |lag| of those maxima; coelution, so lower is better
      LIBRARY_CORR,        ///< Pearson, observed against library intensities
      LIBRARY_DOTPROD,     ///< normalised dot product, same pair
      /// Replaced the old group/window area ratio (D6). That ratio carried no
      /// library or co-elution information at all: a narrow decoy spike in an
      /// otherwise empty window approaches 1.0, while a real target peak on a
      /// real baseline scores lower -- which is why targets measured WORSE on
      /// it (0.1582 against 0.1613). This is the fraction of the group's area
      /// contributed by the fragments the library says should be brightest,
      /// which a single-transition spike cannot satisfy.
      INTENSITY_SCORE,
      LOG_SN,              ///< log apex over a data-derived background floor
      /// How many fragments actually carried information. A precursor scored
      /// from three live fragments is not the same evidence as one scored from
      /// twelve, and without this the two are indistinguishable to the
      /// classifier.
      /// Fragments whose own maximum coincides with the group apex (within one
      /// cycle). Co-elution depth, not signal presence: the previous definition
      /// counted non-degenerate traces and returned a constant 12.000 for right
      /// and wrong answers alike.
      USABLE_FRAGMENTS,

      // Added toward the 10-14 orthogonal scores every working implementation
      // in this family uses. Rosenberger 2017 used 14, DIA-NN uses 73; seven
      // extracting a 58/42 edge from data holding an 8x separation is a scorer
      // roughly half the width it needs, not one missing dominant feature.
      // Names follow pyProphet's, so the columns are comparable to a published
      // weight vector rather than privately invented.

      /// Root-mean-square deviation of normalised observed against normalised
      /// library intensities. Carries what correlation discards: dot product
      /// and Pearson are both scale-free, so a spectrum with the right SHAPE
      /// but the wrong contrast scores identically to a correct one.
      LIBRARY_RMSD,
      /// Share of the group's corrected area in y-ions. Tryptic y-ions carry
      /// the basic C-terminal residue and dominate a real spectrum; a decoy's
      /// mutated sequence has no reason to preserve the ratio.
      YSERIES_SCORE,
      /// Fraction of the fragments the library lists that were seen at all.
      /// Distinct from USABLE_FRAGMENTS, which counts traces carrying any
      /// signal: this counts those carrying signal ABOVE their own background
      /// inside the candidate, which is the co-elution the group claims.
      FRAGMENT_COVERAGE,

      // --- added 2026-08-07. The first four are free: every quantity was
      // already computed and thrown away. See doc/13.

      /// The summed pairwise correlation the co-elution detector used to accept
      /// this position as a peak at all.
      ///
      /// This is the quantity that moved rank-1 accuracy from 40.7% to 75.7%,
      /// and until now the picker computed it per candidate and discarded it.
      /// DIA-NN keeps both `best_corr_sum` and `total_corr_sum` as normalising
      /// context for exactly this reason. Zero for the amplitude picker, which
      /// never computes it.
      CORR_SUM,

      /// How far ahead of its own runner-up this candidate was, on CORR_SUM.
      ///
      /// Distinguishes "this precursor had one obvious answer" from "three
      /// equally plausible ones", which no per-candidate score can express.
      /// DIA-NN encodes the same idea at detection time as MaxCorrDiff, keeping
      /// candidates by margin rather than rank; as a feature it reaches the
      /// classifier instead of only the picker.
      CANDIDATE_MARGIN,

      /// Width of the candidate in cycles, over the run's median candidate
      /// width. A peptide elutes on the chromatography's timescale;
      /// interference need not.
      PEAK_WIDTH_RATIO,

      // RT_DELTA -- |apex RT - predicted RT| -- was here and is REMOVED.
      //
      // A peak group has ONE retention time and its traces co-elute by
      // construction, so a group-level retention-time delta cannot separate a
      // real group from an interference group: both sit wherever the signal
      // they were built from sits. Retention time is diagnostic BETWEEN traces,
      // not for the group as a whole.
      //
      // Measured before removal (d9_auc_by_abundance.py, full S08): AUC 0.215
      // against the decoy null, flat across every abundance quintile, and it
      // got WORSE under random decoy rows (0.293 -> 0.215), so it was not an
      // artefact of picking decoys by argmax. Target and decoy row
      // distributions were identical -- median 41.7 s against 42.6 s, in a
      // window of half-width ~75 s -- i.e. a true peak sits no closer to the
      // predicted time than a random candidate does.
      //
      // It was also pinned in `nonpositive_features` -- but that clip applies
      // to the LDA weight vector only and the default classifier is GBT, so
      // the pin never took effect on the shipped path. Removing the feature
      // measured as a WASH on the fixture (+20 IDs, FDP +1.03 pp, both inside
      // noise), which is what an uninformative feature a tree ensemble already
      // ignores should do. The removal is right because the quantity cannot
      // discriminate, not because it was costing identifications.

      /// |library 1/K0 - observed 1/K0| for the precursor, or NaN where the run
      /// or the library has no mobility. Orthogonal to everything above on
      /// diaPASEF and simply absent elsewhere.
      IM_DELTA,

      /// Pearson of the MS1 precursor trace against the MS2 fragment consensus
      /// over the candidate's own cycles.
      ///
      /// The first sub-score here that does not read MS2 fragment traces. That
      /// is the whole point: the other fifteen share a failure mode, because a
      /// co-eluting interferent corrupts all of them at once. This asks whether
      /// the PRECURSOR rises and falls with its fragments, which no amount of
      /// fragment-side interference can fake.
      ///
      /// Measured on S08 + v6_50k before being added (doc/13): 13.7x enrichment
      /// in the top bin, 1.4-1.5x in bulk, and it survives stratification by MS1
      /// intensity so it is shape rather than brightness. NaN when the run
      /// carries no MS1 or the precursor has no MS1 signal -- NaN, not zero,
      /// because zero is a legitimate correlation and the two must not be
      /// confused.
      MS1_COELUTION,

      // --- added 2026-08-09. Both reviewers, independently, named fragment
      // mass accuracy as the highest-value score ODIA does not have, and both
      // pointed at the same reason: an extracted chromatogram records intensity
      // inside a mass window over time and throws away WHERE in that window the
      // peak sat. An interferent can co-elute perfectly, correlate perfectly
      // and be systematically displaced in m/z, and no chromatogram-shape
      // statistic can see that. It is also the cheapest orthogonal channel
      // available: both quantities are already computed and stored on the
      // group, they were simply never offered to the classifier.

      /// Negated |deviation - the run's own centre|, ppm. Larger is better.
      ///
      /// CENTRED ON THE RUN, which is what makes it safe. The objection that
      /// kept this out of the feature vector was real: fed raw, the classifier
      /// would learn "this run's fragments sit at -3 ppm" and penalise
      /// correctly calibrated identifications. Centring on the run's own median
      /// removes exactly that, and what is left is the per-group departure from
      /// the instrument's systematic error -- a property of the GROUP.
      MASS_ACCURACY,

      /// Negated scatter of the group's per-fragment deviations, ppm. Larger is
      /// better.
      ///
      /// The stronger of the two, and calibration-free by construction: it is a
      /// spread WITHIN one group, so no centring can affect it and no run-level
      /// offset can leak in. A real peptide's fragments are one molecule
      /// measured through one calibration and sit together; interference is
      /// unrelated species whose errors scatter.
      MASS_SPREAD,

      // --- added 2026-08-10. The mobility analogue of MASS_SPREAD, and the
      // first sub-score that can separate interference m/z and RT cannot.

      /// Scatter of this group's PER-FRAGMENT observed 1/K0, as a robust sigma.
      ///
      /// A precursor's fragments are produced from ONE ion packet and therefore
      /// share its mobility. A co-eluting interferent is a different packet at a
      /// different 1/K0, so a group whose fragments disagree about mobility is a
      /// mixture -- and on diaPASEF that is detectable even when the species are
      /// inseparable in m/z and in retention time, because a frame merges TIMS
      /// scans that the raw data still resolves.
      ///
      /// Distinct from IM_DELTA, which asks whether the group sits where the
      /// LIBRARY predicted. This asks whether the group is internally
      /// consistent, which needs no prediction and so cannot be wrong for the
      /// library's reasons.
      ///
      /// NaN on a run without ion mobility, and until at least
      /// `min_im_fragments` fragments carried one -- NaN rather than zero,
      /// because zero scatter is what a perfect group looks like.
      IM_SPREAD,

      // --- added 2026-08-23. The RETENTION-TIME analogue of MASS_SPREAD and
      // IM_SPREAD, and the only form in which retention time can discriminate
      // at all.
      //
      // RT_DELTA asked whether the GROUP sits where the library predicted, and
      // was removed: a group has one retention time and its traces co-elute by
      // construction, so both a real group and an interference group sit
      // wherever the signal they were built from sits. There is nothing there
      // to separate.
      //
      // BETWEEN fragments there is. A peptide's fragments are produced from one
      // ion packet and share one elution profile, so their individual apices
      // agree. An interfering fragment belongs to a different species and peaks
      // somewhere else. That is a property of the group's internal consistency
      // and needs no prediction, which is what makes it immune to the library
      // being wrong -- the same argument that makes IM_SPREAD work.

      /// Weighted scatter of the fragments' own retention-time centroids about
      /// THEIR OWN WEIGHTED MEAN, in seconds, NEGATED so that higher is better
      /// like every other sub-score.
      ///
      /// Dispersion BETWEEN fragments, not displacement from the group apex:
      /// that apex comes from the summed trace, so one loud interferent would
      /// drag the reference and the measurement together.
      ///
      /// NaN for candidates narrower than 3 cycles -- with one or two, every
      /// centroid is forced to the same value and the scatter is 0, the BEST
      /// possible score, so a three-fragment noise spike would look perfect.
      ///
      /// The centroid rather than the argmax: on a faint fragment the argmax
      /// jumps between adjacent cycles on noise, and the whole point of the
      /// feature is to work where the shape features stop working (measured:
      /// at the faintest abundance quintile var_library_corr is 0.532 and
      /// var_xcorr_shape 0.598, against 0.913 and 0.962 at the brightest).
      ///
      /// Weighted by each fragment's background-corrected area, so a fragment
      /// that is mostly noise does not get an equal vote on where the peak is.
      ///
      /// NaN with fewer than `min_rt_spread_fragments` informative fragments --
      /// NaN rather than zero, because zero scatter is what a perfect group
      /// looks like and a group with one fragment would score perfectly.
      ///
      /// Related to XCORR_COELUTION but not the same statistic: that is a
      /// cross-correlation LAG between fragment PAIRS, this is each fragment's
      /// offset from the group apex. Whether it is redundant with it is a
      /// question for measurement, not for this comment.
      RT_SPREAD,

      /// What FRACTION of a fragment's matched intensity survives tightening
      /// the mass tolerance, averaged over fragments.
      ///
      /// The idea is DIA-NN's -- a real peak sits on its calibrated m/z and
      /// survives tightening, while an interferent that merely fell inside the
      /// window often does not -- but the provenance claimed here first was
      /// wrong and is corrected. DIA-NN's 1x / 0.45x / 0.20x tolerances are its
      /// NINE MS1 channels, feeding tight-tolerance CORRELATIONS
      /// (pMs1TightOne/Two); its per-fragment mass feature is pAcc[6], a
      /// deviation rather than a surviving-intensity fraction. So 4.5 and 2.0
      /// ppm are an MS1 constant imported into MS2, and they are not what
      /// DIA-NN does with fragments.
      ///
      /// Measured on real S08 frames rather than assumed: a cell is a MIXTURE,
      /// and 95.1% of bright real-peak cells contain more than one peak. The
      /// weighted mean therefore pulls crowded cells toward the window centre,
      /// which cuts targets from a nominal 1.0 to 0.925 at 4.5 ppm AND lifts
      /// decoy-like cells from 0.45 to 0.564 -- a nominal gap of 0.55 measured
      /// at 0.36. At 2.0 ppm, 43% of perfectly on-mass bright cells fail purely
      /// from co-window neighbours, so that setting measures local crowding
      /// rather than mass correctness. **4.5 is the value to test; 2.0 is not.**
      ///
      /// No re-extraction was needed. The extractor already stores
      /// sum(intensity * ppm) and sum(intensity) per cell, so the tightened
      /// query is a filter over cells that are already in memory rather than a
      /// second pass over the raw data.
      ///
      /// LABEL-SYMMETRIC in the way the two changes reverted before it were
      /// not. It consults no library intensity, which is what broke those: a
      /// decoy copies its target's per-fragment intensities verbatim while its
      /// fragment m/z ARE recomputed, so any intensity-weighted statistic
      /// silently asks a different question of each class. This asks both the
      /// same question -- how far is the matched signal from the m/z THIS
      /// precursor's fragment should have -- and a decoy's recomputed m/z makes
      /// that a genuine test rather than a scrambled one.
      ///
      /// Orthogonal by construction to everything already here, which the
      /// project's own rule says is the lever rather than count: the existing
      /// mass features are a deviation and a scatter, both summaries of WHERE
      /// the matched peaks sit. This is how much intensity is still there when
      /// the window closes, which a median deviation cannot express -- one
      /// fragment can have a perfect median and lose most of its area.
      ///
      /// NaN when the ppm planes are absent (`-collect_mass_residuals` off) or
      /// no fragment matched a peak. Note the imputation hazard that bit an
      /// earlier feature: PercolatorEngine fills NaN with the column median, so
      /// a candidate that could not be measured is scored as an average one.
      /// Here the NaN case is a whole-run structural absence rather than a
      /// per-candidate weakness, which is the case that imputation was designed
      /// for.
      MASS_SURVIVAL,

      /// A DELIBERATELY UNINFORMATIVE COLUMN. Not a feature -- a control.
      ///
      /// Six unrelated changes have now been measured on the fixture at matched
      /// entrapment FDP, and every one produced the same profile: slightly
      /// negative at 5.72-10% and +6 to +8% at 15%. Six coincidences is not a
      /// hypothesis. The alternative is that ADDING A COLUMN AT ALL perturbs the
      /// semi-supervised classifier's trajectory, and that at the operating
      /// point the perturbation is larger than anything a single feature's
      /// information content contributes.
      ///
      /// This column decides between those. It is a deterministic hash of the
      /// precursor index and apex cycle, scaled to [0,1): it varies per
      /// candidate, so it is not constant and survives the constant-column
      /// guard, and it is uniform with respect to label, so it carries no
      /// information a classifier could legitimately use.
      ///
      /// If an arm carrying it reproduces that profile, then the fixture cannot
      /// resolve feature-level changes at this effect size and six "failures to
      /// convert" collapse into one measurement artefact. If it comes back flat,
      /// the six results stand and the features really were not worth their
      /// slots.
      ///
      /// Off unless `-null_feature` is given. It must never be on in a real run.
      ///
      /// IT REPRODUCED THE PROFILE. Pure noise gives -1.5% at 7.42%, -3.9% at
      /// 10% and +5.7% at 15%, inside the range the four real changes gave
      /// (-0.9 to -3.1%, +6.6 to +8.6%). So the fixture cannot resolve a
      /// feature-sized change: adding any column moves the classifier's
      /// trajectory further than a feature's information content does.
      NULL_CONTROL,

      N_SUB_SCORES
    };

    static const std::vector<std::string>& subScoreNames();

    /// Test seam for the boundary rule. The rule lives in an anonymous
    /// namespace in the .cpp, which is right for it and leaves no way to assert
    /// candidate GEOMETRY -- and geometry was where the defect was: 77.4% of
    /// production peak groups covered more than 80% of the extraction window.
    static std::pair<std::size_t, std::size_t> peakBoundsForTest(
      const std::vector<double>& smoothed, std::size_t left_from,
      std::size_t right_from, double boundary_fraction,
      std::size_t min_cycles, std::size_t max_half, double sigmas = 1.0);

    /// Why a library precursor produced no scored candidate.
    ///
    /// Every precursor gets EXACTLY ONE of these, which is the whole point: the
    /// aggregate reject counters cannot be cross-tabulated against a list of
    /// precursors, so "41% of DIA-NN's confident set yields no candidate" could
    /// be attributed to a stage only by inference -- and that inference had
    /// holes. `Session::add` returns silently on `tc == 0`, and a precursor
    /// covered by no isolation window never reaches `add` at all, so neither
    /// appears in any counter.
    enum class TerminalReason : std::uint8_t
    {
      NotReached = 0,     ///< never handed to the scorer: no isolation window, or filtered upstream
      NoTransitions = 1,  ///< reached it with zero transitions extracted
      FewPoints = 2,      ///< fewer than 3 points in the chromatogram
      GateC = 3,          ///< co-elution evidence below the decoy-null quantile
      FewExcursions = 4,  ///< too few transitions rose above their own noise
      ZeroTrace = 5,      ///< the summed trace really is zero
      NoCandidate = 6,    ///< entered the picker and it returned nothing
      Scored = 7,         ///< produced at least one scored candidate
      /// The picker returned candidates and the scorer then discarded every
      /// one of them -- `min_fragments_at_apex` is the only rule that does
      /// this, and it runs AFTER the candidate list is non-empty. Marking
      /// `Scored` at the earlier point would have reported these as successes
      /// that simply produced no output rows.
      AllCandidatesDropped = 8,
      /// Split out of `NotReached`, because they need opposite fixes and the
      /// funnel question turns on which one it is: a precursor no isolation
      /// window covers is an acquisition-scheme fact, while one the prefilter
      /// dropped is our own choice. Written by the EXTRACTOR, which is the only
      /// stage that knows -- see `ChromatogramExtractor::Options`.
      NoWindowCoverage = 9,
      PrefilterExcluded = 10,
    };

    /// Weighting for RT_SPREAD's per-fragment centroids. See the option.
    enum class RtSpreadWeight { Area, Sqrt, None };

    struct Options
    {
      /// Candidate peak groups kept per precursor.
      ///
      /// Not one. Keeping only the best makes the classifier's job impossible
      /// whenever the true peak is second, and -- more quietly -- gives the
      /// decoy population nothing to be wrong about, which deflates the FDR.
      std::size_t max_candidates = 3;

      /// Half-width, in cycles, of the moving average applied before peak
      /// picking. Smoothing the summed trace only; the sub-scores see the
      /// unsmoothed data.
      std::size_t smooth_half_width = 2;

      /// A candidate's boundaries extend until the summed trace falls below
      /// this fraction of its apex.
      double boundary_fraction = 0.10;

      /// Smallest candidate width, in cycles. Boundaries are widened
      /// symmetrically to reach it.
      ///
      /// NOT cosmetic: several sub-scores stop existing below a width.
      /// MS1_COELUTION needs 5 cycles, and the mass and mobility blocks need
      /// hi > lo, so a 1-3 cycle candidate returns NaN for MASS_ACCURACY,
      /// MASS_SPREAD, IM_DELTA, IM_SPREAD and MS1_COELUTION at once. RT_SPREAD
      /// is worse than absent there -- over three time points the fragment
      /// centroids have almost no freedom, so the scatter is mechanically small
      /// and the feature reports agreement it has not measured.
      ///
      /// This was briefly 5, on a sweep that measured the minimum width while it
      /// still governed SCORING. `score_half_cycles` took that job, and the two
      /// jobs want opposite widths, so the sweep's answer moved with the job:
      /// its optimum -- a 5-cycle window -- is now `score_half_cycles = 2`.
      /// What remains here is quantification and the reported RT range, where
      /// the argument runs the other way, because a minimum that truncates a
      /// peak loses area that a narrower scoring window is free to ignore. Back
      /// to 7 accordingly: 9.7 s at S08's 1.385 s cycle, about 2.8x the
      /// measured 3.5 s FWHM, so roughly +-2 sigma either side of the apex.
      ///
      /// The sub-score list above is therefore no longer the reason for this
      /// number -- none of those features read this interval any more. It is
      /// kept because it still describes what a too-narrow group costs, and
      /// because `score_half_cycles = 0` restores the old behaviour and with it
      /// the old constraint.
      std::size_t peak_min_cycles = 7;

      /// Largest half-span a boundary walk may take, in cycles.
      ///
      /// Replaces a bound of n/4, which was a catastrophe guard rather than a
      /// peak-width constraint: it made the widest admissible peak depend on
      /// the EXTRACTION WINDOW, so identical chromatography admitted different
      /// peaks at different window settings, and at 130 cycles it still allowed
      /// ~90 s against a 3.5 s FWHM. 20 cycles is 27.7 s, about 8x the FWHM.
      std::size_t peak_max_half_cycles = 20;

      /// Half-width of the moving average used for BOUNDARY DETECTION only.
      ///
      /// Separate from `smooth_half_width` (2, a 5-point window) because that
      /// is twice the measured peak width: a moving average broader than the
      /// peak lowers its apex, broadens it, and merges it with its neighbours,
      /// which is precisely the failure a boundary rule must not have. 1 gives
      /// a 3-point window, the widest that does not exceed a 2.5-cycle peak.
      std::size_t boundary_smooth_half = 1;

      /// Boundary floor in robust sigmas of local noise above the local
      /// baseline; the larger of this and `boundary_fraction` is used.
      ///
      /// Measured on 3,427 confident positives: the baseline is a median of 36%
      /// of the apex and exceeds a tenth of it for 87.6% of them, so the
      /// fractional floor was unreachable for most real peptides.
      double boundary_sigmas = 1.0;

      /// Half-width, in cycles, of the window the SUB-SCORES are computed over.
      /// 0 restores the previous behaviour, where the sub-scores used the
      /// walked boundaries.
      ///
      /// Separate from the walked boundaries because the two intervals answer
      /// different questions. The walked pair measures how far the peak
      /// extends, which is what quantification and the reported RT range need.
      /// A sub-score is asking whether the evidence AT THIS POSITION fits the
      /// peptide, and for that a wider interval is strictly worse: it dilutes
      /// the correlation with neighbouring signal and makes the score depend
      /// less on where the candidate is -- in the limit, on a window-wide
      /// interval, a correct candidate and a wrong one 20 cycles away received
      /// the SAME sub-scores, median paired difference exactly 0.
      ///
      /// 2, giving a 5-cycle window. A fixed window beats the walked bounds by
      /// **0.0047** AUC, 95% CI [0.0011, 0.0082] over 2,000 precursor-level
      /// bootstrap resamples.
      ///
      /// That figure was first recorded as 0.0136 and is corrected here. The
      /// larger number was measured while `peak_min_cycles` was briefly 5;
      /// restoring it to 7 widens the walked bounds and moves the walked arm
      /// from 0.751 to 0.760 against the fixed window's 0.764. Same data, same
      /// code, different value of a parameter that was changed for unrelated
      /// reasons two commits later -- which is the ordinary way a quoted number
      /// goes stale.
      ///
      /// So the split is worth about a third of what was first claimed, and the
      /// case for it rests more on the architecture than on the margin: it is
      /// what DIA-NN does, it is what the invariance test can hold in place, and
      /// it costs nothing. The effect also reverses where the walk found a
      /// genuinely broad peak (walked width >= 9, 21% of candidates: 0.838
      /// walked against 0.832 fixed) and in the far-only band, so it is a small
      /// average gain over a mixed population, not a uniform one.
      ///
      /// For scale, the change this one sits on top of -- making the boundary
      /// floor reachable at all -- is worth **+0.253**. That is the result; this
      /// is a refinement of it.
      ///
      /// Cycles, not seconds, like every other width here -- which makes it a
      /// portability hazard across acquisitions with a different cycle time.
      /// See the note on `peak_min_cycles`.
      std::size_t score_half_cycles = 2;

      /// Smallest gap, in cycles, between two emitted candidates. 1 is one
      /// candidate per scan position, which is what DIA-NN does and what this
      /// has always done; above 1 it becomes non-maximum suppression, so a cap
      /// of `max_candidates` means that many distinct PEAKS rather than that
      /// many samples of possibly one.
      ///
      /// 1 by default even though the wider setting is measured to raise recall,
      /// because it is measured to cost something too. Suppressing within a
      /// basin at a separation of 5 raises recall@3 from 72.7% to 75.5% and
      /// drops selection accuracy from 63.3% to 57.7% -- the freed slots go to
      /// genuinely different peaks, which can genuinely outscore the true one.
      /// The selector in that measurement is ONE sub-score and the pipeline's is
      /// nineteen, so the trade may well reverse; but it would have to reverse
      /// by 7.4 points of conditional accuracy, and a one-feature proxy cannot
      /// establish that. This exists so the experiment that can is a flag and
      /// not a patch.
      std::size_t candidate_min_separation = 1;

      /// Weight of library-intensity agreement when the cap chooses which of
      /// the margin survivors to keep. 0 is pure `corr_sum`, which is what this
      /// detector has always used.
      ///
      /// Non-zero is measured BETTER on agreement with DIA-NN (recall of the
      /// correct position within a cap of 3 rises 69.5% to 72.7%) and NOT
      /// better reference-free, which is the ruler that governs: target
      /// fraction among the top-N, decoys as the control, is 91.0% against
      /// 96.8% when the discriminant is library correlation and 94.8% against
      /// 93.8% when it is co-elution. The split between those two is the
      /// winner's curse -- selecting on a statistic and then ranking by it lets
      /// decoys shop for their best value as freely as targets -- and the
      /// shipped classifier does see LIBRARY_CORR, so the pipeline sits on the
      /// unfavourable side of it.
      ///
      /// 0 by default for that reason, not because the idea is wrong. The
      /// evidence it adds is genuinely orthogonal to shape, which is what the
      /// detector was short of; what is unproven is that adding it to the
      /// SELECTION survives target-decoy competition.
      double select_library_weight = 0.0;

      /// Tightened mass tolerance, in ppm, for MASS_SURVIVAL. 0 leaves the
      /// sub-score NaN for every candidate, which the constant-column guard
      /// then drops -- so 0 is genuinely off rather than a column of zeros.
      ///
      /// DIA-NN tightens to 0.45x and 0.20x of the extraction window. On S08
      /// the window is 10 ppm, so 4.5 and 2.0 are the corresponding values.
      /// Expressed in absolute ppm rather than as a fraction because the
      /// scorer does not know the extraction window, and a fraction of an
      /// unknown is worse than a number someone had to choose.
      double mass_survival_ppm = 0.0;

      /// Where the run's fragment deviations actually sit, in ppm. MASS_SURVIVAL
      /// measures its window about this rather than about zero.
      ///
      /// 0 is correct whenever the mass calibration succeeded, because the
      /// offset is folded into the query m/z and the residual is centred by
      /// construction. It is NOT correct when the calibration gate fails: that
      /// path sets the offset to 0 with a 50 ppm window and leaves the
      /// distribution where the instrument put it -- about -9 ppm on S08 --
      /// where a window about zero scores real fragments near 0 and uniform
      /// noise higher, inverting the feature.
      double mass_survival_centre = 0.0;

      /// Emit the NULL_CONTROL column. A diagnostic, never a setting: see the
      /// enum. Off leaves it NaN, which the constant-column guard drops.
      /// GBT stability knobs. 0 keeps the model's own default (120 trees,
      /// depth 4, lambda 1.0, min_child_rows 20).
      ///
      /// Exposed because the discriminant is measurably unstable to an
      /// uninformative column -- Spearman rho 0.79-0.82 and top-500 agreement
      /// 0.38-0.42 between a run and the same run plus a noise feature. A
      /// shallower or more strongly regularised model has less freedom to pick
      /// a different equally-good solution when the feature set moves, and
      /// whether that buys stability is now a measurement rather than an
      /// argument.
      int gbt_max_depth = 0;
      int gbt_min_child_rows = 0;
      double gbt_lambda = 0.0;

      /// Compute the GBT's histogram bin edges from ALL rows rather than from
      /// the evolving semi-supervised training set. See GBTParams::fixed_bins:
      /// the default re-discretises every feature whenever the positive set
      /// moves, which turns a one-column perturbation into a global one.
      bool gbt_fixed_bins = false;

      /// Boosting rounds and shrinkage. 0 keeps the model's own 120 and 0.1.
      ///
      /// The pairing that matters is depth 2 with MORE trees: shallow trees
      /// added additively recover capacity without creating the small-leaf
      /// overfitting surface that lets an uninformative column win splits at
      /// depth 4. Halving the learning rate while doubling the rounds is the
      /// usual way to buy accuracy at fixed variance.
      int gbt_n_trees = 0;
      double gbt_learning_rate = 0.0;

      bool null_feature = false;

      /// Seed for NULL_CONTROL. Different seeds give different arbitrary
      /// trajectories, which is how the width of the band is measured rather
      /// than inferred from one sample.
      std::uint64_t null_feature_seed = 0;

      /// Make NULL_CONTROL an epsilon-jittered copy of CORR_SUM instead of a
      /// uniform hash. The hash loses every split tie by position and so proves
      /// only a fixed point; a near-duplicate of the strongest feature WILL win
      /// splits while carrying nothing new, which is the actual test of whether
      /// a configuration attributes an added column correctly.
      bool null_feature_dup = false;

      /// Maximum lag, in cycles, considered by the cross-correlations. Capped
      /// internally at (n-1)/2 of the shortest trace, so a 5-point candidate
      /// never reports a lag resting on one point of overlap.
      int max_delay = 10;

      /// Minimum fragments with real signal at the apex before a candidate is
      /// emitted. A peak group is a co-elution; one transition is a spike.
      std::size_t min_fragments_at_apex = 3;

      /// Candidate picking standardises each transition against its OWN local
      /// noise (median subtracted, divided by MAD) before summing.
      ///
      /// Three formulations were considered and two rejected:
      ///
      /// * raw sum -- one bright interference transition sets the apex, which
      ///   is what produced peak groups that were spikes rather than
      ///   co-elutions;
      /// * unit-max per trace -- rescales a pure-noise trace's largest
      ///   fluctuation to 1.0, so twelve noise traces sum to a taller spurious
      ///   peak than a weak real one, handing decoys free structure;
      /// * library-intensity weighting -- rejected on review. Decoys carry a
      ///   library too, so weighted picking preferentially finds interference
      ///   matching the DECOY's expected pattern, and the same weights then
      ///   feed library_corr and library_dotprod. That closes a feedback loop
      ///   which inflates decoy library scores, and confining the weighting to
      ///   candidate placement does not open it again -- where a candidate
      ///   sits determines what gets scored.
      ///
      /// MAD standardisation weights each transition relative to its own noise
      /// rather than to library expectation or to raw intensity, so it is
      /// independent of everything the sub-scores later measure.
      bool noise_normalised_picking = true;

      /// A transition counts as carrying signal when it rises this many sigma
      /// above its own local noise. 0 restores the old sum-against-zero test.
      double empty_trace_sigma = 3.0;

      /// Gate C: admit a precursor when its co-elution evidence reaches the
      /// (1-alpha) quantile of the run's OWN decoy null. 0 disables it and
      /// falls back to the excursion gate.
      ///
      /// alpha was documented as "the false-admit rate, by construction". It is
      /// NOT, in deployment: measured, 0.05 produced 16.5% decoy admission, 3.3x
      /// off, because tau is fixed from the first `gate_calibration_n` decoys
      /// and precursors arrive in retention-time order, so the calibration
      /// sample is the earliest eluters. The construction is what breaks it.
      /// Background floor for var_log_sn, as a fraction of the candidate apex.
      /// It caps the reported ratio at 1/frac, so 0.01 caps at 100:1 -- and
      /// measured on S08, ALL 21,055 known-present precursors sat at that cap
      /// (10 distinct values, p10 = p50 = p90 = log(100) = 4.605). The feature
      /// was therefore a constant for present precursors. 0.001 caps at 1000:1;
      /// `min(10.0, ...)` bounds the runaway case regardless.
      double log_sn_floor_frac = 0.01;

      double gate_alpha = 0.05;

      /// "quantile" (default) is Gate C: a run-wide decoy-null quantile.
      /// "prominence" admits per precursor from the statistic's OWN noise
      /// model, with no null, no calibration sample and no dependence on
      /// library composition or arrival order.
      std::string gate_mode = "quantile";

      /// Look-elsewhere threshold in sigma for gate_mode = "prominence".
      /// k = Phi^-1((1-alpha_P)^(1/M_eff)); at alpha_P = 0.05 that is 3.1 at
      /// M_eff 50, 3.3 at 100, 3.5 at 200. A starting point for a sweep, not a
      /// derived constant.
      double gate_k = 3.3;

      /// Write every Gate C decision here: precursor, decoy flag, the
      /// statistic, tau in force, whether tau was ready, and the verdict.
      /// Empty disables it.
      ///
      /// The gate's false-NEGATIVE rate on real targets has never been
      /// measured -- alpha bounds decoy admission by construction and nothing
      /// bounds target rejection. A Python re-implementation on a fixture put
      /// it near 84%, but that fixture could not reproduce deployed decisions:
      /// every precursor the production run identified had by definition passed
      /// the real gate, yet only 38% of them cleared the fixture's tau.
      /// Aggregate percentages cannot settle this; these per-precursor
      /// decisions from the deployed code can.
      std::string gate_log_path;

      /// Decoy statistics to collect before the threshold is fixed. Those
      /// precursors are admitted unconditionally and scored normally; 20,000
      /// against a ~10M library is 0.2%.
      std::size_t gate_calibration_n = 20000;

      /// Half-width of the smoothing window, in cycles. A real peak spans
      /// several; a single bright cycle in one transition must not carry it.
      std::size_t gate_smooth_half = 2;

      /// How many transitions must show that excursion. Two, matching the
      /// picker's own "at least 2 fragments present" bar -- one transition
      /// above noise is a spike, not a peak group.
      std::size_t empty_trace_min_transitions = 2;

      /// Detect candidates by co-elution among the precursor's own fragments,
      /// the way DIA-NN's Searcher::peaks does, instead of by the height of a
      /// summed trace.
      ///
      /// The ordering is the point. A standardised sum is dominated by whatever
      /// is loud, so an amplitude-detected list is ordered by how much signal
      /// is present rather than by whether it is THIS peptide. Measured on S08
      /// against DIA-NN's confident set, the amplitude picker put the true peak
      /// first for 24.4% of precursors, in the top 3 for 47.2%, in the top 25
      /// for 89.8% -- it finds the peak and cannot rank it.
      ///
      /// ON by default, because it is better on every axis measured and not
      /// marginally. Same run, same library, DIA-NN as truth, depth 25:
      ///
      ///                  cand/pr  on-RT  off-RT  precision  rank-1  avail
      ///   amplitude         7.88   1102    4038      21.4%   40.7%  69.1%
      ///   co-elution        8.81   1862     228      89.1%   75.7%  97.4%
      ///
      /// Set false to get the old amplitude detector back.
      bool coelution_picking = true;

      /// Use OpenSWATH's PeakPickerChromatogram instead of either ODIA picker.
      ///
      /// The independent implementation doc/07 step 2 has always required and
      /// nobody has run. Overrides `coelution_picking` when set.
      bool openswath_picking = false;

      /// Keep one `MassResidual` per contributing fragment of every retained
      /// candidate, so a mass model can be fitted from IDENTIFICATIONS rather
      /// than from `MassCalibration::collect`'s own probe.
      ///
      /// Costs nothing extra to compute: the per-fragment medians are already
      /// derived here to make `mass_ppm_spread`, and were being discarded. It
      /// costs memory, hence `max_mass_anchors`.
      bool collect_mass_anchors = false;

      /// Ceiling on retained anchors. ~24 B each, so the default is ~192 MB
      /// against a decode floor measured in gigabytes. Hitting it truncates the
      /// sample in RUN ORDER, which is a retention-time bias, so the count of
      /// what was dropped is reported.
      std::size_t max_mass_anchors = 8000000;

      /// The mass correction the EXTRACTOR already applied, so the harvest can
      /// take it back out. Same four coefficients as
      /// `ChromatogramExtractor::Options`, and they must be the same values.
      ///
      /// This matters more than it looks. The extractor shifts each
      /// transition's target m/z by the correction in force
      /// (`ChromatogramExtractor.cpp:905-916`) and then records the deviation
      /// against the SHIFTED target. So a harvested residual is what is left
      /// AFTER the current model, not the run's raw mass error: on S08 the
      /// harvested deviations centre on -0.26 ppm while the instrument's actual
      /// offset is about -9.4.
      ///
      /// Fitting a model from those without correcting for it produces an
      /// INCREMENT to the model in force, and scoring the model in force
      /// against them charges it twice for a correction it already made -- it
      /// would have looked catastrophically worse than doing nothing. Taking
      /// the correction back out here makes `MassResidual::ppm` mean the same
      /// thing regardless of what was in force, which is the only definition a
      /// calibration input can safely have.
      double applied_ppm_offset = 0.0;
      double applied_ppm_log_slope = 0.0;
      double applied_ppm_slope_per_1000 = 0.0;
      double applied_ppm_ref_mz = 0.0;

      /// Take the UNION of the co-elution candidates and the amplitude ones
      /// (or the OpenSWATH ones, when `openswath_picking` is also set), and
      /// give every candidate its co-elution sum as a number.
      ///
      /// The point is to stop using co-elution as a GATE. A gate cannot be
      /// undone by a classifier, so it has to be right every time; a feature
      /// only has to be informative. This also fixes the picker comparison,
      /// where running either alternative alone left `var_corr_sum` and
      /// `var_candidate_margin` constant and therefore dropped.
      bool union_picking = false;

      /// Signal-to-noise threshold for the OpenSWATH picker. Its default is 1.0;
      /// OpenSwathWorkflow commonly runs it at 0.1 for DIA, where a "peak" sits
      /// on far more background than in targeted MRM.
      double openswath_sn = 1.0;

      /// Gaussian smoothing rather than Savitzky-Golay in the OpenSWATH picker.
      bool openswath_gauss = false;

      /// Expected peak width in seconds for the OpenSWATH picker, or 0 to leave
      /// its default. S08's peaks are ~20-30 s.
      double openswath_peak_width = 0.0;

      /// Save the trained discriminant here, or load a frozen one from here.
      ///
      /// Static modelling: at a 1.5% true-positive rate the semi-supervised loop
      /// has no confident positives to bootstrap from and certifies nothing,
      /// while the same discriminant RANKS 232 of DIA-NN's 738 into its top 738.
      /// Training where the loop ignites and applying the frozen model where it
      /// does not is the documented remedy.
      std::string classifier_model_out;
      std::string classifier_model_in;

      /// MS1 traces for MS1_COELUTION, or null when the run has no MS1.
      ///
      /// A real pointer to real data, not a plumbing point. doc/07 is explicit
      /// that a sub-score computed from a placeholder is worse than an absent
      /// one, and IM_DELTA was all-NaN for a week because its hook was added and
      /// never connected. When this is null MS1_COELUTION is NaN for every row
      /// and the constant-column guard drops it, which is the honest behaviour.
      const Ms1Traces* ms1 = nullptr;

      /// How much say each fragment gets in RT_SPREAD's centroid scatter.
      ///
      /// Area weighting estimates the dominant ion packet but suppresses the
      /// signal this feature exists to find: for two clusters the weighted
      /// variance scales as W1*W2/(W1+W2)^2, so one weak interfering fragment
      /// contributes almost nothing. None gives every informative fragment an
      /// equal vote, including barely-detected noisy ones. Sqrt is between.
      RtSpreadWeight rt_spread_weight = RtSpreadWeight::Area;

      /// Informative fragments required before RT_SPREAD is computed at all.
      ///
      /// Three is the smallest number for which a scatter means anything: with
      /// two, the weighted sigma is a rescaled |difference| and every group with
      /// two agreeing fragments looks perfect. The same reasoning that put a
      /// four-point floor under the library correlation.
      std::size_t min_rt_spread_fragments = 3;

      /// Where to record each precursor's `TerminalReason`, indexed by library
      /// precursor index. Null disables the accounting entirely.
      ///
      /// Written without a lock, and the reason is NOT that sessions are
      /// per-thread -- a Sink owns exactly one Session and hands it to
      /// extraction. `ChromatogramExtractor` calls `accept` from its serial
      /// emit loop in both paths (`ChromatogramExtractor.cpp:900`, `:1279`;
      /// the threaded phase is matching, which finishes first), so `add` runs
      /// on one thread. That same invariant is what already makes
      /// `Session::result_.groups` safe to push to.
      std::uint8_t* terminal_reason = nullptr;

      /// DIAGNOSTIC ORACLE. Per-precursor retention time, in run seconds, at
      /// which a candidate is FORCED to exist: the admission gates are bypassed
      /// and, if the picker finds nothing within `oracle_rt_tol` of it, a
      /// candidate is synthesised there. NaN means no oracle for that precursor.
      ///
      /// This exists to answer one question and must never be a default: if
      /// admission were perfect, how many precursors would we actually
      /// identify? `d5_yield.py` estimates at most +40.9% by extrapolating
      /// acceptance rates across abundance bins; this measures it instead.
      ///
      /// It uses the answer as input, so anything it produces is an UPPER
      /// BOUND on a real method and its q-values are not meaningful: only
      /// targets get injections, decoys cannot (they have no true retention
      /// time), so the null is not comparable. Read the injected candidates'
      /// DScore against a threshold from an UNORACLED run.
      const float* oracle_rt = nullptr;
      double oracle_rt_tol = 20.0;

      /// Half-window, in cycles, for the pairwise correlation at each position.
      std::size_t corr_half_window = 4;

      /// The reference fragment's summed correlation to the others must reach
      /// this for a position to be a peak at all. DIA-NN's MinCorrScore.
      double min_corr_score = 0.5;

      /// Candidates are kept by MARGIN from the best correlation sum rather
      /// than by rank, so an unambiguous precursor yields one candidate and an
      /// ambiguous one yields several. DIA-NN's MaxCorrDiff. This is why a
      /// fixed top-N hurt: at 25 it manufactured 24 competitors regardless of
      /// whether any was plausible, and precision fell to 21.4%.
      double max_corr_diff = 2.0;

      /// The apex must be this fraction of the maximum evidence nearby on the
      /// reference fragment's own smoothed trace. DIA-NN's PeakApexEvidence.
      double apex_evidence = 0.99;

      /// Reject a candidate whose observed spectrum correlates with the
      /// library below this. -1.0 disables it, which is the default.
      ///
      /// This is the one sub-score measured to separate real identifications
      /// from misplaced ones: on S08, median 0.582 for calls landing within
      /// 30 s of the true apex against -0.036 for those that do not -- and
      /// -0.032 for decoys, i.e. a misplaced target is indistinguishable from a
      /// decoy here. See the comment at the gate in PeakGroupScorer.cpp.
      double min_library_corr = -1.0;

      /// Classifier for the semi-supervised loop.
      /// "lda" is deterministic and dependency-light; "gbt" is what the
      /// upstream benchmarks use.
      std::string classifier = "gbt";

      unsigned threads = 0;

      /// True when `Library::precursors().irt` holds RUN SECONDS rather than
      /// library iRT units -- i.e. after the retention-time map has been fitted
      /// and applied. It centres pass 2's extraction window; no sub-score
      /// reads it any more (RT_DELTA was removed). Before it,
      /// the comparison would be between two different units.
      bool library_rt_is_run_seconds = false;

      /// Observed 1/K0 per precursor, indexed as the library is, or empty.
      /// Supplied by the caller because the scorer never sees spectra.
      const std::vector<float>* observed_im = nullptr;

      /// Sub-score indices to withhold from the classifier, by name on the
      /// command line. Ablation, which the project has needed for a while and
      /// faked twice by comparing different binaries -- a comparison that also
      /// changes whatever else moved between them.
      ///
      /// Implemented by flattening the column to a constant, so the existing
      /// constant-column guard drops it and SAYS SO in the log. A NaN column
      /// would be silently imputed somewhere and the arm would claim to have
      /// ablated something it did not.
      std::vector<int> disabled_sub_scores;

      /// Draw each decoy's best score from as many candidates as a target has.
      /// See `Scoring::LDAParams::match_decoy_candidate_counts`.
      bool match_decoy_candidate_counts = false;

      /// Semi-supervised loop knobs, previously reachable only by recompiling.
      /// 0 / negative means "leave the LDAParams default alone".
      double train_fdr_initial = 0.0;
      double train_fdr = 0.0;
      int classifier_iterations = 0;
      bool use_pi0 = false;
    };

    struct PeakGroup
    {
      std::uint32_t precursor = 0;
      float apex_rt = 0.0f;
      float left_rt = 0.0f;
      float right_rt = 0.0f;
      float apex_intensity = 0.0f;
      /// Median m/z deviation, ppm, over this group's matched fragment peaks,
      /// and how many contributed. NaN/0 when the extractor did not collect it.
      /// Deliberately NOT a sub-score -- see where it is filled.
      float mass_ppm = std::numeric_limits<float>::quiet_NaN();
      std::uint16_t mass_ppm_n = 0;

      /// SCATTER of this group's PER-FRAGMENT deviations, ppm, as a robust
      /// sigma (MAD x 1.4826) about the group's own median. NaN until at least
      /// `MassWidth::min_fragments` fragments contributed.
      ///
      /// This is the quantity a window width must be sized from, and it is NOT
      /// the spread of `mass_ppm` across groups. `mass_ppm` is a median over
      /// (fragment x cycle) cells, so its own precision is sigma/sqrt(N_eff) --
      /// sizing a window from how tightly group medians cluster would give a
      /// number several times too narrow and would look, wrongly, like a very
      /// well calibrated instrument. What has to fit inside the window is one
      /// FRAGMENT's deviation, so one fragment is the unit measured here.
      float mass_ppm_spread = std::numeric_limits<float>::quiet_NaN();

      /// The run's OBSERVED 1/K0 for this group, intensity-weighted over the
      /// candidate's cycles. NaN on a run with no ion mobility.
      float observed_im = std::numeric_limits<float>::quiet_NaN();
      /// Robust sigma (MAD x 1.4826) of this group's per-fragment observed
      /// 1/K0, and how many fragments carried one.
      float im_spread = std::numeric_limits<float>::quiet_NaN();
      std::uint8_t im_frags = 0;
      /// How many fragments contributed to `mass_ppm_spread`.
      std::uint8_t mass_ppm_frags = 0;
      std::vector<double> sub_scores;

      double dscore = 0.0;
      double qvalue = 1.0;
      double pep = 1.0;
      bool decoy = false;
    };

    /// One accepted fragment's mass residual, tagged with the group it came
    /// from so the FDR can filter it afterwards.
    ///
    /// The tag is necessary because a residual is produced while scoring, and
    /// whether its group is worth fitting from is only known after `finish()`.
    /// Fitting a calibration from every candidate would fit it from the
    /// interference too.
    struct MassAnchor
    {
      /// `residual.decoy` is NOT populated here -- target/decoy status belongs
      /// to the group, and holding a second copy on the residual is what let
      /// the two disagree once. `acceptedMassResiduals` fills it from the group.
      MassResidual residual;
      std::uint32_t group = 0;   ///< index into Result::groups
    };

    struct Result
    {
      std::vector<PeakGroup> groups;

      /// Per-fragment mass residuals, when `collect_mass_anchors` was on.
      ///
      /// One entry per (retained candidate x contributing fragment), NOT per
      /// accepted group -- acceptance is decided after these are produced. Use
      /// `acceptedMassResiduals` to reduce them.
      std::vector<MassAnchor> mass_anchors;

      /// Residuals discarded because `max_mass_anchors` was reached. Non-zero
      /// here means the sample is truncated at an arbitrary point in the run,
      /// i.e. biased towards early retention times -- it is reported rather
      /// than silently absorbed for exactly that reason.
      std::size_t mass_anchors_dropped = 0;

      /// Precursors that produced no candidate at all -- no non-zero point in
      /// the window. Reported rather than silently absent.
      std::size_t precursors_without_candidate = 0;

      std::size_t target_groups = 0;
      std::size_t decoy_groups = 0;

      /// False when the q-values must not be read as an FDR at all.
      ///
      /// Target/decoy FDR is undefined without decoys, and the failure is
      /// silent rather than loud: with no negative class every row looks
      /// perfect, q comes out 0 everywhere, and the run reports that it
      /// identified everything it was given. Measured on a 400-precursor slice
      /// of a target-only library, that is exactly what happened -- 400 of 400
      /// "identified at 1% FDR" -- so this is a guard against a specific
      /// observed failure, not a hypothetical one.
      bool fdr_valid = false;

      /// Straight from the scorer. Zero trained iterations means the d-scores
      /// are a single-feature initialisation rather than a fitted model, which
      /// looks like a working classifier and is not one.
      int iterations_trained = 0;
      int iterations_skipped = 0;

      /// Targets at q <= 0.01, counted on the best group per precursor.
      std::size_t identified_at_1pct = 0;

      /// Candidates dropped by `min_library_corr`. Reported rather than
      /// silently absent: a filter that discards without saying so is
      /// indistinguishable from a search that found nothing.
      std::size_t candidates_below_library_corr = 0;
    };

    /// Reduce harvested anchors to those from groups the FDR accepted.
    ///
    /// Separate from the harvest because acceptance is not known when a
    /// residual is produced, and because the threshold is a caller's decision:
    /// a calibration fitted at 1% and one fitted at 5% are different
    /// experiments, and doc/15 section 12.5 records that the lenient choice was
    /// measured to be the worse one.
    ///
    /// Decoy groups are returned too, with `MassResidual::decoy` set, because
    /// they are the null a mixture fit needs -- but note doc/15 section 3.3:
    /// pseudo-reverse decoys are NOT the right null for the interference
    /// component. Use the m/z-shifted control for that.
    static std::vector<MassResidual> acceptedMassResiduals(const Result& result,
                                                           double q_threshold,
                                                           bool include_decoys = false);

    /// Scoring one precursor at a time.
    ///
    /// The candidate search and every sub-score are per precursor already --
    /// nothing in them looks at another precursor's trace. Only the last step,
    /// fitting the discriminant and calibrating the FDR, is global, and it
    /// needs the score MATRIX rather than the chromatograms. So a precursor can
    /// be scored the moment its chromatogram is complete and the points thrown
    /// away, which is what lets extraction hold only what is live.
    ///
    /// The whole-`Chromatograms` entry point below is this class driven over a
    /// flat array, so the two cannot drift apart.
    class Session
    {
    public:
      Session(const Library& library, const Options& options);

      /// Score one precursor. Nothing about @p trace is retained after this
      /// returns, which is the contract that lets the caller free it.
      void add(const PrecursorChromatogram& trace);

      /// Fit the discriminant over everything added, and calibrate.
      ///
      /// Peak groups come out ordered by precursor whatever order they were
      /// added in. Extraction hands them over in retention-time order, and a
      /// feature matrix whose row order depended on when a peptide eluted would
      /// make the fit depend on it too.
      Result finish();

      /// Sub-scores decided after extraction; see `Sink::disableSubScores`.
      void disableSubScores(std::vector<int> indices)
      { options_.disabled_sub_scores = std::move(indices); }
      void setAppliedMassCorrection(double offset, double log_slope,
                                    double slope_per_1000, double ref_mz)
      {
        options_.applied_ppm_offset = offset;
        options_.applied_ppm_log_slope = log_slope;
        options_.applied_ppm_slope_per_1000 = slope_per_1000;
        options_.applied_ppm_ref_mz = ref_mz;
      }

    private:
      /// Gate C's decoy null, OWNED BY THIS SESSION.
      ///
      /// It used to be a file-scope singleton, so tau was calibrated once per
      /// PROCESS and pass 2 inherited pass 1's threshold -- pass 1 extracts
      /// wide, pass 2 narrow, and the statistic is a max over the window, so
      /// pass 2's values are systematically smaller than the tau they were
      /// compared against. An epoch counter was tried first and is not
      /// sufficient: with two Sessions constructed up front, the first one
      /// calibrates under the second's epoch and the second then inherits its
      /// tau, which is the original bug. Ownership is the only version that
      /// cannot be defeated by construction order.
      struct GateNull;
      std::shared_ptr<GateNull> gate_null_;
      const Library* library_;
      Options options_;

    public:
      /// Install the run's MS1 traces after construction. See
      /// ChromatogramSink::ms1Available for why this cannot be a constructor
      /// argument.
      void setMs1Traces(const Ms1Traces* m) { options_.ms1 = m; }

    private:
      Result result_;
    };

    /// A `ChromatogramSink` that scores each precursor and drops it.
    ///
    /// This is the default path at scale: it has no term proportional to the
    /// library except the peak groups it produces, which are ~120 bytes each
    /// against a chromatogram's tens of kilobytes.
    class Sink final : public ChromatogramSink
    {
    public:
      Sink(const Library& library, const Options& options) : session_(library, options) {}

      /// See ChromatogramSink::ms1Available -- the traces do not exist yet when
      /// this sink is constructed.
      void ms1Available(const Ms1Traces* m) override { session_.setMs1Traces(m); }
      void accept(const PrecursorChromatogram& trace) override { session_.add(trace); }
      Result finish() { return session_.finish(); }

      /// Tell the harvest what mass correction the extractor is applying.
      ///
      /// Called after the Sink is built for the same reason `disableSubScores`
      /// is: the mass probe runs during extraction setup, so the coefficients do
      /// not exist at construction time. Must be called before the first
      /// `accept()`, which is where the harvest reads them.
      void setAppliedMassCorrection(double offset, double log_slope,
                                    double slope_per_1000, double ref_mz)
      { session_.setAppliedMassCorrection(offset, log_slope, slope_per_1000, ref_mz); }

      /// Withhold sub-scores decided AFTER extraction.
      ///
      /// `-mass_features auto` keys off the fragment mass calibration's
      /// verdict, and that verdict does not exist when the Sink is built -- the
      /// probe runs during extraction setup. Deciding at construction time
      /// therefore left pass 1 running with features that pass 2 would reject,
      /// and pass 1 is where the retention-time and mobility anchors come from,
      /// so it is not a harmless inconsistency: measured on S08 it cost 117
      /// identifications (1,225 against 1,342).
      void disableSubScores(std::vector<int> indices)
      { session_.disableSubScores(std::move(indices)); }

    private:
      Session session_;
    };

    static Result score(const Library& library, const Chromatograms& chromatograms,
                        const Options& options);

    /// Refit over groups that already carry their sub-scores.
    ///
    /// The candidate picker is retention-time agnostic -- neither
    /// `findCandidates` nor `findCandidatesByCorrelation` reads the library's
    /// iRT -- and RT_DELTA used to be the one sub-score that depended on the
    /// fitted map, which is what made "new map -> one recomputed column -> one
    /// refit" cheap enough to iterate.
    ///
    /// RT_DELTA is now REMOVED, so nothing in `sub_scores` depends on the map
    /// and refitting after a new one is bit-identical. `refitsChangeScores()`
    /// reports that, and the refinement loop consults it rather than spending
    /// rounds to discover it. The map still centres pass 2's extraction window,
    /// which is upstream of anything here.
    ///
    /// Note this retired nothing that was running: `-refine_rounds` defaults to
    /// 0, so the loop's `max_rounds == 0` early return fires first and the
    /// guard never executes in a default run. It is there for the day someone
    /// turns the loop on.
    static void refit(const Library& library, Result& result, const Options& options);

    /// Whether a new retention-time map can change any score through `refit`.
    /// False while no sub-score reads the map; flip the one return in the .cpp
    /// when a map-dependent sub-score is added.
    static bool refitsChangeScores();

  private:
    static void fitAndAssign_(const Library& library, Result& result,
                              const Options& options);

  public:
  };

} // namespace ODIA
