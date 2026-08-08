// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/MassCalibration.h>

#include <cstddef>
#include <string>
#include <vector>

namespace ODIA
{

  /// The run's fragment mass error as a function of BOTH m/z and retention time.
  ///
  /// `MassCalibration` fits one global correction: a constant plus an optional
  /// shape in m/z. That model cannot express drift, and drift is the normal
  /// behaviour of a mass spectrometer over an hour of gradient -- source
  /// contamination, temperature, and space charge all move the measured mass
  /// slowly and monotonically.
  ///
  /// This fits ppm(m/z, RT) as PIECEWISE CONSTANT IN RT and LINEAR IN log m/z:
  ///
  ///     ppm = intercept[block(rt)] + slope[block(rt)] * log(mz / ref_mz)
  ///
  /// log m/z rather than m/z because that is the axis the error is actually
  /// linear on -- `MassCalibration` already reports its shape that way
  /// ("log_mz: 3.322 ppm per e-fold in m/z"), and a linear-in-m/z fit over a
  /// 200-1800 Th range is dominated by the top of the range.
  ///
  /// FITTED FROM ANCHORS, which is the whole point. `MassCalibration`'s probe
  /// chooses its own cells by brightest-co-occurrence over the WHOLE run, and on
  /// a mostly-absent library that maximum is a chance event whose residuals are
  /// flat -- which is why its gate fails on S08 and reports "a mostly-noise
  /// sample". It is describing a sample its own search ruined. Anchors are
  /// precursors this run has already SCORED confidently, taken at their own
  /// apex, so the residuals are those of real fragments at the right moment.
  ///
  /// Measured 2026-08-08 with an RT-shifted control (the discipline this class
  /// must be held to): an unanchored median m/z residual of -8.44 ppm was almost
  /// entirely a property of the +/-50 ppm search, since the same probe 300 s away
  /// from the peptide gave -4.98 ppm. The real offset is the difference, ~-3.5
  /// ppm. **Any offset this class produces must be checked the same way.**
  class MassRecalibration
  {
  public:
    struct Options
    {
      /// Retention-time blocks over the gradient. Equal-COUNT, not equal-width:
      /// a gradient's identifications are not uniform in time, and equal-width
      /// blocks would fit the sparse ends from almost nothing.
      std::size_t blocks = 8;

      /// Minimum anchors in a block before it is fitted at all. Below this the
      /// block inherits the global fit rather than inventing a local one.
      std::size_t min_per_block = 150;

      /// Minimum anchors before the m/z SLOPE is fitted within a block. A
      /// slope needs a spread of masses, not just a count.
      std::size_t min_for_slope = 400;

      /// Largest correction, ppm. A fit that wants more than this is not a
      /// calibration; it is the fit chasing interference.
      double max_correction_ppm = 20.0;

      /// Largest slope magnitude, ppm per e-fold in m/z.
      double max_slope_ppm = 15.0;

      /// Trim residuals beyond this many MADs before fitting. The anchor set is
      /// confident but not pure, and one decade-wide outlier moves a least
      /// squares far more than it moves a median.
      double trim_mads = 4.0;
    };

    struct Block
    {
      double rt_low = 0.0, rt_high = 0.0;
      double intercept_ppm = 0.0;      ///< at ref_mz
      double slope_ppm = 0.0;          ///< per e-fold in m/z
      std::size_t anchors = 0;
      bool fitted = false;             ///< false: inherits the global fit
    };

    bool fitted() const { return fitted_; }
    double refMz() const { return ref_mz_; }
    const std::vector<Block>& blocks() const { return blocks_; }

    /// The ppm to ADD to a theoretical m/z at retention time @p rt.
    double ppmAt(double mz, double rt) const;

    /// The largest correction this model applies anywhere, ppm. The extractor
    /// needs it to widen its search slack, since the exact per-candidate test
    /// happens after the bucket lookup.
    double maxAbsCorrection() const;

    std::string describe() const;

    /// Fit from residuals collected AT SCORED APEXES. @p residuals with
    /// `decoy` set are ignored: a shifted control cannot inform a correction.
    static MassRecalibration fit(const std::vector<MassResidual>& residuals,
                                 const Options& options);

  private:
    std::vector<Block> blocks_;
    double ref_mz_ = 500.0;
    double global_intercept_ = 0.0, global_slope_ = 0.0;
    double max_correction_ppm_ = 20.0;   ///< bounds the COMPOSITE, see ppmAt
    bool fitted_ = false;
  };

} // namespace ODIA
