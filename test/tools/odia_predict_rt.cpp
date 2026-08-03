// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Predicts iRT for the given peptides so the values can be compared against
// test/peptdeep_reference.py. Tensor agreement is not prediction agreement if
// dtype, layout or ordering assumptions differ, so this compares the thing that
// actually matters.

#include <odia/PeptDeepPredictor.h>

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <iomanip>
#include <iostream>
#include <vector>

int main(int argc, char** argv)
{
  if (argc < 3)
  {
    std::cerr << "usage: odia_predict_rt <model.onnx> <sequence>...\n";
    return 2;
  }
  try
  {
    std::vector<OpenMS::AASequence> peptides;
    for (int i = 2; i < argc; ++i)
    {
      peptides.push_back(OpenMS::AASequence::fromString(argv[i]));
    }

    ODIA::PeptDeepPredictor predictor(argv[1]);
    const auto rt = predictor.predictRT(peptides);

    std::cerr << "provider: "
              << (predictor.provider() == ODIA::PeptDeepPredictor::Provider::CUDA
                    ? "CUDA" : "CPU")
              << "\n";
    for (std::size_t i = 0; i < rt.size(); ++i)
    {
      std::cout << argv[i + 2] << "\t" << std::fixed << std::setprecision(6)
                << rt[i] << "\n";
    }
  }
  catch (const std::exception& e)
  {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
