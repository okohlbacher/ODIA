// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Dumps ODIA's PeptDeep encoding as JSON so it can be diffed against
// test/peptdeep_reference.py, which implements the same spec independently.

#include <odia/PeptDeepEncoder.h>

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <iostream>

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::cerr << "usage: odia_encode_dump <modified-sequence>\n";
    return 2;
  }
  try
  {
    const auto peptide = OpenMS::AASequence::fromString(argv[1]);
    const auto batch = ODIA::PeptDeepEncoder::encode(peptide);
    const auto width = ODIA::PEPTDEEP_MOD_ELEMENTS.size();

    std::cout << "{\"aa_indices\": [";
    for (std::size_t i = 0; i < batch.aa_indices.size(); ++i)
    {
      std::cout << (i ? ", " : "") << batch.aa_indices[i];
    }
    std::cout << "], \"mod_x\": {";
    bool first_row = true;
    for (std::size_t r = 0; r < batch.sequence_length; ++r)
    {
      std::string entries;
      for (std::size_t c = 0; c < width; ++c)
      {
        const float v = batch.mod_x[r * width + c];
        if (v != 0.0f)
        {
          if (!entries.empty()) { entries += ", "; }
          entries += "\"" + std::to_string(c) + "\": " + std::to_string(v);
        }
      }
      if (!entries.empty())
      {
        std::cout << (first_row ? "" : ", ") << "\"" << r << "\": {" << entries << "}";
        first_row = false;
      }
    }
    std::cout << "}}\n";
  }
  catch (const std::exception& e)
  {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
