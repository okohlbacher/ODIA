// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/PeptDeepPredictor.h>

#include <odia/PeptDeepEncoder.h>

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <onnxruntime_cxx_api.h>

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace ODIA
{

  struct PeptDeepPredictor::Impl
  {
    Ort::Env env{ORT_LOGGING_LEVEL_WARNING, "odia"};
    Ort::SessionOptions options;
    std::unique_ptr<Ort::Session> session;
    Ort::MemoryInfo memory{Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault)};
    Provider provider = Provider::CPU;
    std::vector<std::string> input_names;
    std::vector<std::string> output_names;
  };

  PeptDeepPredictor::PeptDeepPredictor(const std::string& model_path,
                                       bool prefer_gpu, int intra_op_threads)
    : impl_(std::make_unique<Impl>())
  {
    if (intra_op_threads > 0) { impl_->options.SetIntraOpNumThreads(intra_op_threads); }
    impl_->options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);

    if (prefer_gpu)
    {
      // Ask before attempting. Registering CUDA blindly on a CPU-only build
      // still falls back correctly, but ONNX Runtime logs a failure to load
      // libonnxruntime_providers_shared.so at error level first -- alarming
      // output for the expected case on a machine with no GPU, which is the
      // machine this is developed on.
      bool cuda_available = false;
      for (const auto& name : Ort::GetAvailableProviders())
      {
        if (name == "CUDAExecutionProvider") { cuda_available = true; break; }
      }

      if (cuda_available)
      {
        try
        {
          OrtCUDAProviderOptions cuda{};
          impl_->options.AppendExecutionProvider_CUDA(cuda);
          impl_->provider = Provider::CUDA;
        }
        catch (const Ort::Exception&)
        {
          impl_->provider = Provider::CPU;
        }
      }
    }

    try
    {
      impl_->session = std::make_unique<Ort::Session>(impl_->env, model_path.c_str(),
                                                      impl_->options);
    }
    catch (const Ort::Exception& e)
    {
      if (impl_->provider == Provider::CUDA)
      {
        // Registering the provider can succeed while creating the session with
        // it fails, so the fallback has to cover both.
        impl_->provider = Provider::CPU;
        impl_->options = Ort::SessionOptions{};
        impl_->options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
        if (intra_op_threads > 0) { impl_->options.SetIntraOpNumThreads(intra_op_threads); }
        impl_->session = std::make_unique<Ort::Session>(impl_->env, model_path.c_str(),
                                                        impl_->options);
      }
      else
      {
        throw std::runtime_error("cannot open model " + model_path + ": " + e.what());
      }
    }

    Ort::AllocatorWithDefaultOptions allocator;
    for (std::size_t i = 0; i < impl_->session->GetInputCount(); ++i)
    {
      impl_->input_names.emplace_back(
        impl_->session->GetInputNameAllocated(i, allocator).get());
    }
    for (std::size_t i = 0; i < impl_->session->GetOutputCount(); ++i)
    {
      impl_->output_names.emplace_back(
        impl_->session->GetOutputNameAllocated(i, allocator).get());
    }

    // Validate by name and count, not just count. Inputs are bound positionally,
    // and the CCS and MS2 models take three and five inputs respectively -- both
    // previously constructed cleanly here and then failed inside Run with a
    // red "Missing Input" from ONNX Runtime, which is exactly the alarming
    // output the provider probing was added to avoid.
    // Two inputs is the RT model, five is MS2. Anything else is neither.
    const bool looks_like_rt =
      impl_->input_names.size() == 2 &&
      (impl_->input_names[0] == "input_sequences" || impl_->input_names[0] == "aa_indices") &&
      impl_->input_names[1] == "mod_x";
    const bool looks_like_ms2 =
      impl_->input_names.size() == 5 &&
      (impl_->input_names[0] == "input_sequences" || impl_->input_names[0] == "aa_indices") &&
      impl_->input_names[1] == "mod_x" && impl_->input_names[2] == "charges" &&
      impl_->input_names[3] == "nce" && impl_->input_names[4] == "instrument_indices";
    const bool looks_like_ccs =
      impl_->input_names.size() == 3 &&
      (impl_->input_names[0] == "input_sequences" || impl_->input_names[0] == "aa_indices") &&
      impl_->input_names[1] == "mod_x" && impl_->input_names[2] == "charges";
    if (!looks_like_rt && !looks_like_ms2 && !looks_like_ccs)
    {
      std::string names;
      for (const auto& n : impl_->input_names) { names += (names.empty() ? "" : ", ") + n; }
      throw std::runtime_error(
        "model " + model_path + " is not a PeptDeep retention-time, CCS or MS2 model: "
        "found inputs (" + names + ")");
    }
  }

  PeptDeepPredictor::~PeptDeepPredictor() = default;

  PeptDeepPredictor::Provider PeptDeepPredictor::provider() const
  {
    return impl_->provider;
  }

  void PeptDeepPredictor::runBatch_(const PeptDeepEncoder::Batch& batch,
                                    const std::vector<std::size_t>& group,
                                    std::vector<float>& out)
  {
    const std::int64_t rows = static_cast<std::int64_t>(batch.rows);
    const std::int64_t length = static_cast<std::int64_t>(batch.sequence_length);
    const std::int64_t width = static_cast<std::int64_t>(PEPTDEEP_MOD_ELEMENTS.size());

    std::array<std::int64_t, 2> aa_shape{rows, length};
    std::array<std::int64_t, 3> mod_shape{rows, length, width};

    std::vector<Ort::Value> inputs;
    inputs.push_back(Ort::Value::CreateTensor<std::int64_t>(
      impl_->memory, const_cast<std::int64_t*>(batch.aa_indices.data()),
      batch.aa_indices.size(), aa_shape.data(), aa_shape.size()));
    inputs.push_back(Ort::Value::CreateTensor<float>(
      impl_->memory, const_cast<float*>(batch.mod_x.data()),
      batch.mod_x.size(), mod_shape.data(), mod_shape.size()));

    std::vector<const char*> in_names{impl_->input_names[0].c_str(),
                                      impl_->input_names[1].c_str()};
    std::vector<const char*> out_names{impl_->output_names[0].c_str()};

    auto results = impl_->session->Run(Ort::RunOptions{nullptr}, in_names.data(),
                                       inputs.data(), inputs.size(),
                                       out_names.data(), out_names.size());

    // Validate what came back rather than trusting the shape. The MS2 model
    // returns [batch, seq_len-3, 8]; read blindly, the same loop would hand back
    // the first few floats of its b_z1 channel as retention times.
    const auto info = results.front().GetTensorTypeAndShapeInfo();
    if (info.GetElementType() != ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT)
    {
      throw std::runtime_error("model output is not float32");
    }
    if (info.GetElementCount() != group.size())
    {
      throw std::runtime_error(
        "model returned " + std::to_string(info.GetElementCount()) + " values for " +
        std::to_string(group.size()) + " peptides; this does not look like a "
        "retention-time model");
    }

    const float* values = results.front().GetTensorData<float>();
    for (std::size_t i = 0; i < group.size(); ++i)
    {
      // Back to the caller's order: grouping reorders, and a caller that got
      // predictions silently permuted would have no way to notice.
      out[group[i]] = values[i];
    }
  }

  std::vector<PeptDeepPredictor::Spectrum>
  PeptDeepPredictor::predictMS2(const std::vector<OpenMS::AASequence>& peptides,
                                const std::vector<int>& charges, float nce,
                                const std::string& instrument,
                                std::vector<Failure>* failures)
  {
    if (charges.size() != peptides.size())
    {
      throw std::invalid_argument("one charge per peptide is required");
    }
    if (impl_->input_names.size() != 5)
    {
      throw std::runtime_error("this is not a PeptDeep MS2 model: it takes " +
                               std::to_string(impl_->input_names.size()) +
                               " inputs, not 5");
    }

    // NCE is one value for the whole call, so a bad one is not a per-peptide
    // failure and must not be reported as 2048 of them. It is also silent: the
    // model returns a plausible spectrum for any float it is given.
    if (!std::isfinite(nce) || nce <= 0.0f || nce > 1000.0f)
    {
      throw std::invalid_argument("collision energy " + std::to_string(nce) +
                                  " is not a usable NCE");
    }

    std::vector<Spectrum> out(peptides.size());

    for (const auto& whole_group : PeptDeepEncoder::groupByLength(peptides))
    {
     for (std::size_t offset = 0; offset < whole_group.size(); offset += MAX_BATCH_ROWS)
     {
      const std::vector<std::size_t> group(
        whole_group.begin() + static_cast<std::ptrdiff_t>(offset),
        whole_group.begin() + static_cast<std::ptrdiff_t>(
          std::min(offset + MAX_BATCH_ROWS, whole_group.size())));

      std::vector<OpenMS::AASequence> subset;
      std::vector<int> subset_charges;
      for (const auto index : group)
      {
        subset.push_back(peptides[index]);
        subset_charges.push_back(charges[index]);
      }

      PeptDeepEncoder::Batch batch;
      try
      {
        batch = PeptDeepEncoder::encode(subset, subset_charges, nce, instrument);
      }
      catch (const std::exception&)
      {
        // As in predictRT: one unencodable peptide, or one impossible charge,
        // must not discard its whole chunk and take 2047 good peptides with it.
        // Retry singly so only the offending rows are recorded -- reporting the
        // whole chunk also gave every peptide a reason naming another one's
        // modification.
        if (failures == nullptr) { throw; }
        for (std::size_t i = 0; i < group.size(); ++i)
        {
          try
          {
            const auto single = PeptDeepEncoder::encode(
              std::vector<OpenMS::AASequence>{subset[i]},
              std::vector<int>{subset_charges[i]}, nce, instrument);
            runMS2Batch_(single, {group[i]}, out);
          }
          catch (const std::exception& inner)
          {
            failures->push_back({group[i], inner.what()});
          }
        }
        continue;
      }

      runMS2Batch_(batch, group, out);
     }
    }
    return out;
  }

  void PeptDeepPredictor::runMS2Batch_(const PeptDeepEncoder::Batch& batch,
                                       const std::vector<std::size_t>& group,
                                       std::vector<Spectrum>& out)
  {
    {
      {
      const std::int64_t rows = static_cast<std::int64_t>(batch.rows);
      const std::int64_t length = static_cast<std::int64_t>(batch.sequence_length);
      const std::int64_t width = static_cast<std::int64_t>(PEPTDEEP_MOD_ELEMENTS.size());

      std::array<std::int64_t, 2> aa_shape{rows, length};
      std::array<std::int64_t, 3> mod_shape{rows, length, width};
      std::array<std::int64_t, 2> meta_shape{rows, 1};
      // Rank 1, unlike charges and nce. Feeding [rows, 1] is rejected outright.
      std::array<std::int64_t, 1> instrument_shape{rows};

      std::vector<Ort::Value> inputs;
      inputs.push_back(Ort::Value::CreateTensor<std::int64_t>(
        impl_->memory, const_cast<std::int64_t*>(batch.aa_indices.data()),
        batch.aa_indices.size(), aa_shape.data(), aa_shape.size()));
      inputs.push_back(Ort::Value::CreateTensor<float>(
        impl_->memory, const_cast<float*>(batch.mod_x.data()),
        batch.mod_x.size(), mod_shape.data(), mod_shape.size()));
      inputs.push_back(Ort::Value::CreateTensor<float>(
        impl_->memory, const_cast<float*>(batch.charges.data()),
        batch.charges.size(), meta_shape.data(), meta_shape.size()));
      inputs.push_back(Ort::Value::CreateTensor<float>(
        impl_->memory, const_cast<float*>(batch.nces.data()),
        batch.nces.size(), meta_shape.data(), meta_shape.size()));
      inputs.push_back(Ort::Value::CreateTensor<std::int64_t>(
        impl_->memory, const_cast<std::int64_t*>(batch.instrument_indices.data()),
        batch.instrument_indices.size(), instrument_shape.data(),
        instrument_shape.size()));

      std::vector<const char*> in_names;
      for (const auto& n : impl_->input_names) { in_names.push_back(n.c_str()); }
      std::vector<const char*> out_names{impl_->output_names[0].c_str()};

      auto results = impl_->session->Run(Ort::RunOptions{nullptr}, in_names.data(),
                                         inputs.data(), inputs.size(),
                                         out_names.data(), out_names.size());

      const auto info = results.front().GetTensorTypeAndShapeInfo();
      const auto shape = info.GetShape();
      // shape[1] is checked as well as shape[0] and shape[2]. The comparison
      // against the reference cannot detect a wrong positions count -- both
      // sides take it from the same model -- and upstream OpenMS does not
      // trust it either, slicing to size()-1 rather than using what came back.
      if (info.GetElementType() != ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT ||
          shape.size() != 3 || shape[0] != rows ||
          shape[1] != length - 3 ||
          shape[2] != static_cast<std::int64_t>(Spectrum::CHANNELS))
      {
        throw std::runtime_error(
          "MS2 output is not [batch, nAA - 1, 8]: the model returned a tensor of "
          "rank " + std::to_string(shape.size()) + " for " + std::to_string(rows) +
          " peptides of " + std::to_string(length - 2) + " residues");
      }

      const auto positions = static_cast<std::size_t>(shape[1]);
      const float* values = results.front().GetTensorData<float>();
      for (std::size_t i = 0; i < group.size(); ++i)
      {
        auto& spectrum = out[group[i]];
        spectrum.positions = positions;
        // Indexed by group[i], not i: the peptides were reordered into length
        // groups, and writing them back in batch order would silently permute
        // whole spectra between peptides of the same length.
        spectrum.intensities.assign(
          values + i * positions * Spectrum::CHANNELS,
          values + (i + 1) * positions * Spectrum::CHANNELS);
      }
     }
    }
  }

  std::vector<float>
  PeptDeepPredictor::predictRT(const std::vector<OpenMS::AASequence>& peptides,
                               std::vector<Failure>* failures)
  {
    // The constructor accepts both model kinds, so the per-call guard has to be
    // here: feeding an MS2 model two of its five inputs fails inside Run with a
    // red "Missing Input" instead of a sentence a user can act on.
    if (impl_->input_names.size() != 2)
    {
      throw std::runtime_error("predictRT was given a model with " +
                               std::to_string(impl_->input_names.size()) +
                               " inputs; the retention-time model takes 2");
    }

    // NaN, not 0.0: an unfilled entry must be recognisable as "no prediction"
    // rather than pass for a real retention time.
    std::vector<float> out(peptides.size(), std::numeric_limits<float>::quiet_NaN());

    for (const auto& whole_group : PeptDeepEncoder::groupByLength(peptides))
    {
     // Submitted in capped chunks; see MAX_BATCH_ROWS.
     for (std::size_t offset = 0; offset < whole_group.size(); offset += MAX_BATCH_ROWS)
     {
      const std::vector<std::size_t> group(
        whole_group.begin() + static_cast<std::ptrdiff_t>(offset),
        whole_group.begin() + static_cast<std::ptrdiff_t>(
          std::min(offset + MAX_BATCH_ROWS, whole_group.size())));

      std::vector<OpenMS::AASequence> subset;
      subset.reserve(group.size());
      for (const auto index : group) { subset.push_back(peptides[index]); }

      PeptDeepEncoder::Batch batch;
      try
      {
        batch = PeptDeepEncoder::encode(subset);
      }
      catch (const std::exception&)
      {
        // One unencodable peptide must not discard the whole run -- at proteome
        // scale that is minutes of inference lost to a single mass-shift
        // modification -- nor its whole chunk, which would take 2047 good
        // peptides with it. Retry the chunk one peptide at a time so only the
        // offending ones are recorded.
        if (failures == nullptr) { throw; }
        for (std::size_t i = 0; i < group.size(); ++i)
        {
          try
          {
            const auto single = PeptDeepEncoder::encode(
              std::vector<OpenMS::AASequence>{subset[i]});
            runBatch_(single, {group[i]}, out);
          }
          catch (const std::exception& inner)
          {
            failures->push_back({group[i], inner.what()});
          }
        }
        continue;
      }

      runBatch_(batch, group, out);
     }
    }
    return out;
  }

  std::vector<float>
  PeptDeepPredictor::predictCCS(const std::vector<OpenMS::AASequence>& peptides,
                                const std::vector<int>& charges,
                                std::vector<Failure>* failures)
  {
    if (charges.size() != peptides.size())
    {
      throw std::invalid_argument("one charge per peptide is required");
    }
    if (impl_->input_names.size() != 3)
    {
      throw std::runtime_error("predictCCS was given a model with " +
                               std::to_string(impl_->input_names.size()) +
                               " inputs; the CCS model takes 3");
    }

    // NaN, as in predictRT: an unfilled entry must be recognisable as "no
    // prediction" rather than pass for a collision cross-section of zero.
    std::vector<float> out(peptides.size(), std::numeric_limits<float>::quiet_NaN());

    for (const auto& whole_group : PeptDeepEncoder::groupByLength(peptides))
    {
     for (std::size_t offset = 0; offset < whole_group.size(); offset += MAX_BATCH_ROWS)
     {
      const std::vector<std::size_t> group(
        whole_group.begin() + static_cast<std::ptrdiff_t>(offset),
        whole_group.begin() + static_cast<std::ptrdiff_t>(
          std::min(offset + MAX_BATCH_ROWS, whole_group.size())));

      std::vector<OpenMS::AASequence> subset;
      std::vector<int> subset_charges;
      subset.reserve(group.size());
      subset_charges.reserve(group.size());
      for (const auto index : group)
      {
        subset.push_back(peptides[index]);
        subset_charges.push_back(charges[index]);
      }

      PeptDeepEncoder::Batch batch;
      try
      {
        // The CCS model takes no NCE and no instrument, but the encoder's meta
        // overload fills both; they are simply not bound below. Reusing it
        // keeps the charge validation and scaling in one place.
        batch = PeptDeepEncoder::encode(subset, subset_charges, 30.0f, "QE");
      }
      catch (const std::exception&)
      {
        // One unencodable peptide, or one impossible charge, must not discard
        // its chunk. As in predictRT and predictMS2.
        if (failures == nullptr) { throw; }
        for (std::size_t i = 0; i < group.size(); ++i)
        {
          try
          {
            const auto single = PeptDeepEncoder::encode(
              std::vector<OpenMS::AASequence>{subset[i]},
              std::vector<int>{subset_charges[i]}, 30.0f, "QE");
            runCCSBatch_(single, {group[i]}, out);
          }
          catch (const std::exception& inner)
          {
            failures->push_back({group[i], inner.what()});
          }
        }
        continue;
      }

      runCCSBatch_(batch, group, out);
     }
    }
    return out;
  }

  void PeptDeepPredictor::runCCSBatch_(const PeptDeepEncoder::Batch& batch,
                                       const std::vector<std::size_t>& group,
                                       std::vector<float>& out)
  {
    const std::int64_t rows = static_cast<std::int64_t>(batch.rows);
    const std::int64_t length = static_cast<std::int64_t>(batch.sequence_length);
    const std::int64_t width = static_cast<std::int64_t>(PEPTDEEP_MOD_ELEMENTS.size());

    std::array<std::int64_t, 2> aa_shape{rows, length};
    std::array<std::int64_t, 3> mod_shape{rows, length, width};
    std::array<std::int64_t, 2> charge_shape{rows, 1};

    std::vector<Ort::Value> inputs;
    inputs.push_back(Ort::Value::CreateTensor<std::int64_t>(
      impl_->memory, const_cast<std::int64_t*>(batch.aa_indices.data()),
      batch.aa_indices.size(), aa_shape.data(), aa_shape.size()));
    inputs.push_back(Ort::Value::CreateTensor<float>(
      impl_->memory, const_cast<float*>(batch.mod_x.data()),
      batch.mod_x.size(), mod_shape.data(), mod_shape.size()));
    inputs.push_back(Ort::Value::CreateTensor<float>(
      impl_->memory, const_cast<float*>(batch.charges.data()),
      batch.charges.size(), charge_shape.data(), charge_shape.size()));

    std::vector<const char*> in_names;
    for (const auto& n : impl_->input_names) { in_names.push_back(n.c_str()); }
    std::vector<const char*> out_names{impl_->output_names[0].c_str()};

    auto results = impl_->session->Run(Ort::RunOptions{nullptr}, in_names.data(),
                                       inputs.data(), inputs.size(),
                                       out_names.data(), out_names.size());

    const auto info = results.front().GetTensorTypeAndShapeInfo();
    const auto shape = info.GetShape();
    // Rank 1, like the RT model and unlike what the format notes claim. A
    // [batch, 1] output would still have the right element count, so checking
    // the count alone would not notice.
    if (info.GetElementType() != ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT ||
        shape.size() != 1 || shape[0] != rows)
    {
      throw std::runtime_error("CCS output is not [batch]");
    }

    const float* values = results.front().GetTensorData<float>();
    for (std::size_t i = 0; i < group.size(); ++i)
    {
      // Indexed by group[i]: the peptides were reordered into length groups.
      out[group[i]] = values[i];
    }
  }

} // namespace ODIA
