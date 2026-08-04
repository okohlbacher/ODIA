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
      INTENSITY_SCORE,     ///< this group's intensity over the window's total
      LOG_SN,              ///< log apex over local background
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

      /// Maximum lag, in cycles, considered by the cross-correlations.
      int max_delay = 10;

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
    };

    static Result score(const Library& library, const Chromatograms& chromatograms,
                        const Options& options);
  };

} // namespace ODIA
