// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/ChromatogramTsv.h>
#include <odia/TextWriter.h>

#include <cstdint>
#include <string>

namespace ODIA
{

  void writeChromatogramTsv(const std::string& path, const Library& library,
                            const Chromatograms& chromatograms)
  {
    TextWriter out(path);
    out.put("Precursor.Id\tTransition.Index\tProduct.Mz\tRT\tIntensity\n");

    const auto& p = library.precursors();
    const auto& t = library.transitions();
    for (std::size_t i = 0; i < library.precursorCount(); ++i)
    {
      // The library does not store Precursor.Id; DIA-NN's convention is
      // sequence + charge, and the writer reconstructs it the same way the
      // library writer does so the two files join on it.
      const auto seq = library.strings().get(p.modified_sequence[i]);
      const std::string id = std::string(seq) + std::to_string(static_cast<int>(p.charge[i]));
      for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
      {
        const std::uint32_t tr = p.transition_begin[i] + k;
        if (tr >= chromatograms.begin.size()) { continue; }
        const std::uint64_t b = chromatograms.begin[tr];
        const std::uint32_t n = chromatograms.count[tr];
        // The stream this replaces never had its precision set, so every number
        // below is defaultfloat at the default precision of 6. Passing 6
        // explicitly is what keeps the file identical: to_chars' own default is
        // shortest-round-trip, which would widen every m/z here from 6 digits
        // to as many as 17.
        const double product_mz = fromFixed(t.product_mz[tr]);
        for (std::uint32_t j = 0; j < n; ++j)
        {
          out.put(id); out.put('\t');
          out.integer(tr); out.put('\t');
          out.number(product_mz, 6); out.put('\t');
          out.number(chromatograms.retentionTime(tr, j), 6); out.put('\t');
          out.number(chromatograms.intensity[b + j], 6); out.put('\n');
        }
      }
    }
    out.close();
  }

} // namespace ODIA
