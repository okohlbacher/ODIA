// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/PeptDeepEncoder.h>

#include <memory>
#include <string>
#include <vector>

namespace OpenMS { class AASequence; }

namespace ODIA
{

  /// Runs AlphaPeptDeep's PeptDeep models through ONNX Runtime.
  ///
  /// ODIA drives the session itself rather than using OpenMS's PeptDeep
  /// binding, for two reasons that are not preference: the binding rejects
  /// modified peptides outright, and its headers are compiled into libOpenMS.so
  /// but never installed, so no external project can include them. The model
  /// weights still come from OpenMS -- WITH_ONNX=ON downloads them against
  /// pinned checksums -- so nothing is vendored here.
  class PeptDeepPredictor
  {
  public:
    /// Which execution provider the session ended up using.
    enum class Provider { CPU, CUDA };

    /// @param model_path a PeptDeep .onnx file.
    /// @param prefer_gpu attempt CUDA first, falling back to CPU (D4).
    /// @param intra_op_threads 0 leaves it to ONNX Runtime.
    explicit PeptDeepPredictor(const std::string& model_path,
                               bool prefer_gpu = true,
                               int intra_op_threads = 0);
    ~PeptDeepPredictor();

    PeptDeepPredictor(const PeptDeepPredictor&) = delete;
    PeptDeepPredictor& operator=(const PeptDeepPredictor&) = delete;

    /// The provider actually in use. Requesting CUDA on a machine without it is
    /// not an error -- the same binary has to work on both -- so callers that
    /// care must ask rather than assume.
    Provider provider() const;

    /// Largest number of peptides submitted to the model in one call.
    ///
    /// Not a tuning knob. Without a cap, a length group is submitted whole:
    /// 2 M peptides with a tryptic length distribution put 234,593 in the
    /// largest group and peaked at 26.7 GiB, and under a memory limit it did
    /// not degrade but aborted inside the BLAS allocator with no indication of
    /// which stage failed. Capping also fixes the batch size, which the result
    /// depends on at the last bit.
    static constexpr std::size_t MAX_BATCH_ROWS = 2048;

    /// Peptides that could not be encoded, by index into the input, with the
    /// reason. One unencodable peptide must not discard a whole run's work.
    struct Failure
    {
      std::size_t index;
      std::string reason;
    };

    /// Predict normalised iRT, one value per peptide, in the input order.
    ///
    /// Peptides are grouped by length internally: padding is not inert, because
    /// index 0 is one-hot encoded and no model applies a padding mask, so a
    /// mixed-length batch would give every peptide a different answer from the
    /// one it gets alone.
    /// @param failures if given, receives the peptides that could not be
    ///        encoded; their entries in the result are left NaN. If not given,
    ///        an unencodable peptide throws.
    std::vector<float> predictRT(const std::vector<OpenMS::AASequence>& peptides,
                                 std::vector<Failure>* failures = nullptr);

    /// One predicted fragment spectrum: (nAA - 1) positions x 8 channels.
    ///
    /// Channel order comes from alphabase's sort_charged_frag_types, which is
    /// sorted(no_loss) + sorted(loss):
    ///   0 b_z1  1 b_z2  2 y_z1  3 y_z2
    ///   4 b_modloss_z1  5 b_modloss_z2  6 y_modloss_z1  7 y_modloss_z2
    struct Spectrum
    {
      std::size_t positions = 0;
      static constexpr std::size_t CHANNELS = 8;
      std::vector<float> intensities;   ///< [positions, CHANNELS], row-major

      float at(std::size_t position, std::size_t channel) const
      {
        return intensities[position * CHANNELS + channel];
      }
    };

    /// Predict fragment intensities. @p charges is one per peptide.
    std::vector<Spectrum> predictMS2(const std::vector<OpenMS::AASequence>& peptides,
                                     const std::vector<int>& charges,
                                     float nce = 30.0f,
                                     const std::string& instrument = "QE",
                                     std::vector<Failure>* failures = nullptr);

  private:
    void runBatch_(const PeptDeepEncoder::Batch& batch,
                   const std::vector<std::size_t>& group, std::vector<float>& out);

    struct Impl;
    std::unique_ptr<Impl> impl_;
  };

} // namespace ODIA
