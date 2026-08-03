// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

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

    /// Predict normalised iRT, one value per peptide, in the input order.
    ///
    /// Peptides are grouped by length internally: padding is not inert, because
    /// index 0 is one-hot encoded and no model applies a padding mask, so a
    /// mixed-length batch would give every peptide a different answer from the
    /// one it gets alone.
    std::vector<float> predictRT(const std::vector<OpenMS::AASequence>& peptides);

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
  };

} // namespace ODIA
