// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/ChromatogramExtractor.h>
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
    };

    struct PeakGroup
    {
      std::uint32_t precursor = 0;
      float apex_rt = 0.0f;
      float left_rt = 0.0f;
      float right_rt = 0.0f;
      float apex_intensity = 0.0f;
      std::vector<double> sub_scores;

      double dscore = 0.0;
      double qvalue = 1.0;
      double pep = 1.0;
      bool decoy = false;
    };

    struct Result
    {
      std::vector<PeakGroup> groups;

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

    private:
      Session session_;
    };

    static Result score(const Library& library, const Chromatograms& chromatograms,
                        const Options& options);
  };

} // namespace ODIA
