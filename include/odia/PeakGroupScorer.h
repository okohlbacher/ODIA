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

      /// |apex RT - predicted RT|, seconds, on the CALIBRATED axis.
      ///
      /// Absent until now because the header's own reason had not expired: with
      /// no iRT calibration the library was spread evenly over the run and this
      /// would have been noise. Pass 2 runs on a fitted map, so both numbers
      /// are real -- but only then, which is why it is gated on
      /// `library_rt_is_run_seconds` and left at NaN otherwise. That gating IS
      /// DIA-NN's min_iter_learn schedule, arrived at from the other direction.
      ///
      /// mProphet ranks its equivalent last of seven at AUC 0.850. Last of
      /// seven is not zero.
      RT_DELTA,

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

      N_SUB_SCORES
    };

    static const std::vector<std::string>& subScoreNames();

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
      /// alpha IS the false-admit rate, by construction -- that is the point of
      /// calibrating on the null rather than assuming a distribution.
      double gate_alpha = 0.05;

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
      /// and applied. RT_DELTA is only computed when this is set; before it,
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
      const Library* library_;
      Options options_;
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
    /// iRT -- and RT_DELTA is the ONLY sub-score that depends on the fitted
    /// map. So a new map costs one recomputed column and one refit, not a
    /// re-extraction of the run. That is what makes iterating to convergence
    /// affordable: DIA-NN runs twelve iterations for the same reason, against
    /// chromatograms it has already read.
    ///
    /// Recomputes RT_DELTA from the library's CURRENT irt (which must be in run
    /// seconds), then refits the discriminant and reassigns q-values.
    static void refit(const Library& library, Result& result, const Options& options);

  private:
    static void fitAndAssign_(const Library& library, Result& result,
                              const Options& options);

  public:
  };

} // namespace ODIA
