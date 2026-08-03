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
    if (impl_->input_names.size() != 2 ||
        (impl_->input_names[0] != "input_sequences" &&
         impl_->input_names[0] != "aa_indices") ||
        impl_->input_names[1] != "mod_x")
    {
      std::string names;
      for (const auto& n : impl_->input_names) { names += (names.empty() ? "" : ", ") + n; }
      throw std::runtime_error(
        "model " + model_path + " is not a PeptDeep retention-time model: expected "
        "inputs (input_sequences, mod_x) but found (" + names + ")");
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

  std::vector<float>
  PeptDeepPredictor::predictRT(const std::vector<OpenMS::AASequence>& peptides,
                               std::vector<Failure>* failures)
  {
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

} // namespace ODIA
