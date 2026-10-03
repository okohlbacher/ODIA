// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace ODIA
{

  /// Per-frame TOF -> m/z calibration for the Bruker "ims-compact" mzPeak layout.
  ///
  /// That layout stores the raw TOF index, not m/z, and the pinned reader
  /// (opt/mzpeak) rebuilds m/z with ONE run-wide square-root-linear chord,
  /// m/z = (a + b*tof)^2, from the archive index. The vendor model is the same
  /// form with per-frame coefficients: Bruker's ModelType 1 corrects C1 for the
  /// frame's temperature (Frames.T1), and the chord leaves that term out. The
  /// converter (mzpeak-convert 0.12.5) writes the exact pair per spectrum as
  /// `opt_MZP_1000003_tof_c0` / `opt_MZP_1000004_tof_c1` in spectra_metadata.
  ///
  /// Measured on PXD047793 009/010 (shared/pxd/R0B_MASS.md, 24.3 M + 22.0 M
  /// peaks against the Bruker SDK): chord minus SDK = +0.45..+1.4 ppm, smooth
  /// in m/z; the per-frame pair minus SDK = 0.000 ppm on every sampled peak.
  ///
  /// The remap below inverts the chord to recover the integer TOF index and
  /// re-evaluates it with the frame's own pair. Recovering the index is exact
  /// (the chord is strictly monotone and the index is an Int32), and the
  /// residual from the nearest integer is CHECKED: an m/z array that the chord
  /// did not produce would otherwise be "corrected" into nonsense silently.

  /// The largest distance from an integer an inverted TOF may have before the
  /// array is refused. Measured round-off is ~1e-10 index units; one unit is
  /// ~0.1 ppm at 1,000 Th, so 1e-3 cannot mask a wrong chord and cannot trip
  /// on double round-off.
  constexpr double kTofIndexTolerance = 1e-3;

  struct TofRemapStats
  {
    std::size_t peaks = 0;
    /// max |tof - round(tof)| over the remapped peaks.
    double max_index_residual = 0.0;
  };

  /// Rewrite m/z values decoded with the chord (a, b) to the per-frame pair
  /// (c0, c1). Throws std::runtime_error when a value does not invert to an
  /// integer TOF index within kTofIndexTolerance, or when a/b/c0/c1 are not
  /// usable (non-finite, b or c1 not positive).
  ///
  /// When (c0, c1) == (a, b) the output equals the input bit for bit only up
  /// to double round-off of the round trip; callers that must stay byte-
  /// identical do not call this at all (the option is off by default).
  void remapTofMz(std::vector<double>& mz, double a, double b, double c0, double c1,
                  TofRemapStats* stats = nullptr);

  /// The per-spectrum coefficient table, indexed by the file's spectrum index.
  ///
  /// A spectrum whose coefficients are NULL in the file stays on the run-wide
  /// chord: the archive documents that case (NULL Frames.T1, counted in its
  /// `per_spectrum_chord_frames`) and keeps a/b exactly for it.
  class TofCoefficientTable
  {
  public:
    TofCoefficientTable() = default;

    /// Build from parallel columns. @p index holds each row's spectrum index;
    /// @p c0_null / @p c1_null flag NULL cells. Refuses (std::runtime_error)
    /// an index outside [0, n_spectra), a repeated index with different
    /// coefficients, and a row where exactly one of the pair is NULL.
    static TofCoefficientTable fromColumns(const std::vector<std::uint64_t>& index,
                                           const std::vector<double>& c0,
                                           const std::vector<double>& c1,
                                           const std::vector<bool>& c0_null,
                                           const std::vector<bool>& c1_null,
                                           std::size_t n_spectra);

    /// Read `spectra_metadata.parquet` out of an mzPeak archive. Throws when the
    /// entry or either coefficient column is absent: a caller that asked for
    /// per-frame calibration and silently got the chord would measure nothing.
    static TofCoefficientTable fromArchive(const std::string& mzpeak_path, std::size_t n_spectra);

    /// True when spectrum @p i carries its own pair.
    bool has(std::size_t i) const { return i < has_.size() && has_[i]; }
    double c0(std::size_t i) const { return c0_[i]; }
    double c1(std::size_t i) const { return c1_[i]; }

    std::size_t size() const { return has_.size(); }
    /// Spectra with their own pair, and spectra left on the chord (NULL or no row).
    std::size_t withPair() const { return with_pair_; }
    std::size_t onChord() const { return has_.size() - with_pair_; }

  private:
    std::vector<double> c0_;
    std::vector<double> c1_;
    std::vector<bool> has_;
    std::size_t with_pair_ = 0;
  };

  /// ppm shift the per-frame pair applies relative to the chord at @p mz:
  /// ((c0 + c1*t)^2 / mz - 1) * 1e6 with t the chord's TOF for @p mz (not
  /// rounded). For reporting only.
  double tofShiftPpm(double mz, double a, double b, double c0, double c1);

} // namespace ODIA
