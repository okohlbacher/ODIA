// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Dumps a predicted fragment spectrum as JSON so it can be compared against
// test/peptdeep_reference.py.

#include <odia/PeptDeepPredictor.h>

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <iomanip>
#include <iostream>
#include <vector>

int main(int argc, char** argv)
{
  if (argc < 4)
  {
    std::cerr << "usage: odia_predict_ms2 <model.onnx> <sequence> <charge> [nce] [instrument]\n";
    return 2;
  }
  try
  {
    const auto peptide = OpenMS::AASequence::fromString(argv[2]);
    const int charge = std::atoi(argv[3]);
    const float nce = argc > 4 ? std::atof(argv[4]) : 30.0f;
    const std::string instrument = argc > 5 ? argv[5] : "QE";

    ODIA::PeptDeepPredictor predictor(argv[1]);
    const auto spectra = predictor.predictMS2({peptide}, {charge}, nce, instrument);
    const auto& s = spectra.front();

    std::cout << "{\"shape\": [" << s.positions << ", " << ODIA::PeptDeepPredictor::Spectrum::CHANNELS
              << "], \"values\": [";
    for (std::size_t p = 0; p < s.positions; ++p)
    {
      std::cout << (p ? ", " : "") << "[";
      for (std::size_t c = 0; c < ODIA::PeptDeepPredictor::Spectrum::CHANNELS; ++c)
      {
        std::cout << (c ? ", " : "") << std::fixed << std::setprecision(6) << s.at(p, c);
      }
      std::cout << "]";
    }
    std::cout << "]}\n";
  }
  catch (const std::exception& e)
  {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
