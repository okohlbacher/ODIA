// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/Library.h>

#include <cstddef>
#include <string>
#include <vector>

namespace ODIA
{

  /// Per-run retention-time refinement, learned IN PROCESS from the run's own
  /// identifications.
  ///
  /// WHY THIS EXISTS. The retention-time map is a monotone function of the
  /// library's iRT, and a monotone function cannot beat the ordering it is
  /// given. Measured on Astral by fitting on one set of stripped sequences and
  /// evaluating on another, the best monotone calibration reaches SD 38.6 s and
  /// no amount of LOESS span, anchor threshold or bin count goes below it.
  /// DIA-NN's own residual on the same file is 29.29 s, which is BELOW that
  /// ceiling -- possible only because its prediction is run-refined rather than
  /// a monotone map of the library value.
  ///
  /// So the residual left after the best possible calibration is elution-ORDER
  /// error, and only something that reads the SEQUENCE can touch it. That is
  /// this class.
  ///
  /// WHY IT IS NOT THE PEPTDEEP FINE-TUNER. Fine-tuning peptdeep on the run
  /// reaches 27.70 s and is better -- but it needs torch, a separate python
  /// environment, a GPU node, and an ONNX export step on a third machine,
  /// none of which belong inside a search. `-repredict_irt -rt_model` remains
  /// for that path. This one has no dependency at all and runs by default,
  /// because a refinement nobody enables is a refinement nobody gets.
  ///
  /// WHAT IT LEARNS. Ridge regression on amino-acid composition, peptide
  /// length, precursor charge and the CALIBRATED iRT, predicting observed
  /// retention time. Composition is the cheap half of what a sequence model
  /// knows; measured on Astral it takes the held-out residual from 38.55 s to
  /// 33.19 s, against peptdeep's 27.70. A gradient-boosted regressor would sit
  /// between the two and is the obvious next step -- `Scoring::GBT` is a
  /// classifier today and would need a squared-loss path.
  class RtRefiner
  {
  public:
    struct Options
    {
      /// L2 penalty. 1.0 measured best of {0.3, 1, 3, 10, 30, 100} on Astral,
      /// and the curve is flat between 0.3 and 3.
      double ridge = 1.0;

      /// Fraction of TRAINING SEQUENCES held out to report a residual on.
      ///
      /// 0.15 as the project owner specified: production trains on everything
      /// it has because the goal is minimising deviation on THIS run, not
      /// preserving a clean estimate, and holding data out purely to protect a
      /// statistic would be paying identifications for it. The 15% exists to
      /// report an honest number and to refuse the model when it does not help.
      double holdout = 0.15;

      /// Anchors below which no refinement is attempted. A model fitted on a
      /// handful of peptides will predict their neighbours and nothing else.
      std::size_t min_anchors = 250;

      /// Anchors further than this many robust sigmas from the median residual
      /// are dropped before fitting. Pass 1's anchors reach 1,600 s of residual
      /// -- misidentifications, not chromatography -- and a squared loss chases
      /// them.
      double trim_mads = 4.0;

      /// Hard bound on the correction, as a FRACTION of the anchors' own iRT
      /// range. A linear model over composition extrapolates without limit
      /// outside its training range, and an unbounded correction put pass 2's
      /// windows off the gradient entirely: IH1 went 1,464 identifications -> 0
      /// before this existed.
      ///
      /// A fraction rather than an absolute, because this now works in iRT
      /// units, and an iRT axis may be normalised to [-50, 150] or already be
      /// run seconds when -irt_slope was given. 0.05 of the range is a large
      /// chromatographic shift and a small modelling error.
      double max_shift_fraction = 0.05;
    };

    struct Report
    {
      bool fitted = false;
      std::size_t anchors = 0;      ///< rows used
      std::size_t train = 0;
      std::size_t held_out = 0;
      std::size_t trimmed = 0;    ///< anchors dropped as outliers before fitting
      double sd_before = 0.0;       ///< held-out SD of (observed - calibrated iRT), seconds
      double sd_after = 0.0;        ///< held-out SD after refinement, seconds
      std::string note;
    };

    /// Fit from (stripped sequence, charge, LIBRARY iRT, TARGET iRT).
    ///
    /// IN iRT SPACE, not in retention-time space. The calibration is a monotone
    /// map iRT -> RT and it should be the ONLY transform on the RT axis; a
    /// correction applied after it is a second, unconstrained one, which is how
    /// an earlier version put pass 2's windows off the gradient and took IH1 to
    /// zero. So the caller inverts the map to get the iRT each anchor SHOULD
    /// have had -- `target = map^-1(observed RT)` -- and this learns
    /// `sequence -> (target - library)`. The map is then refitted on the
    /// refined axis, so its monotonicity is re-established rather than assumed.
    ///
    /// The split is by STRIPPED SEQUENCE, so a peptide's charge states never
    /// straddle it: 2+ and 3+ of one peptide elute together, and splitting
    /// between them leaks the answer into the held-out set.
    Report fit(const std::vector<std::string>& sequences,
               const std::vector<int>& charges,
               const std::vector<double>& calibrated_irt,
               const std::vector<double>& observed_rt,
               const Options& options);

    /// Same, with the defaults. A defaulted `Options{}` argument cannot be
    /// written here: the member initialisers are not complete until the
    /// enclosing class is, so an overload is the portable spelling.
    Report fit(const std::vector<std::string>& sequences,
               const std::vector<int>& charges,
               const std::vector<double>& calibrated_irt,
               const std::vector<double>& observed_rt);

    /// Rewrite the library's `irt` in place with the refined prediction.
    ///
    /// A no-op when the fit was refused or did not help, so a caller can apply
    /// unconditionally and get the unrefined axis when refinement is not
    /// warranted.
    std::size_t apply(Library& library) const;

    bool fitted() const { return fitted_; }

    /// Write the fitted model so it can be applied to OTHER runs.
    ///
    /// READ THE HAZARD FIRST. This model learns THIS run's chromatography. The
    /// project has already paid for reusing one across runs: library v5 was
    /// generated with a model tuned on a different run and cost 2,027 confident
    /// precursors against v4's stock model. The sidecar carried a warning and
    /// the generation script passed it anyway.
    ///
    /// It is still worth having, because a SERIES of runs on one gradient,
    /// instrument and method is the case this is actually for: fine-tune once
    /// on a representative run, apply to the rest, skip the per-run fit. What
    /// must not happen is a model crossing a gradient or an instrument.
    ///
    /// So the file records what it was fitted on -- run, anchor count, held-out
    /// SD before and after -- and `load` refuses nothing but says all of it
    /// loudly. A machine cannot tell whether two gradients are the same; the
    /// person running it can.
    bool save(const std::string& path, const std::string& provenance) const;

    /// Load a model written by `save`. Returns false if the file is unusable.
    bool load(const std::string& path, std::string* provenance_out = nullptr);

  private:
    bool fitted_ = false;
    std::vector<double> weights_;   ///< in the feature order of `features()`
    double max_shift_ = 120.0;
  };

} // namespace ODIA
