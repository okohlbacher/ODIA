// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/ChromatogramTsv.h>
#include <odia/TextWriter.h>

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>

namespace ODIA
{

  void writeChromatogramTsv(const std::string& path, const Library& library,
                            const Chromatograms& chromatograms,
                            const std::unordered_set<std::string>* keep)
  {
    TextWriter out(path);
    // Decoy is NOT optional. A decoy precursor reconstructs the SAME
    // Precursor.Id as its target (DIA-NN's convention is sequence + charge, and
    // a decoy carries its target's sequence), so without this column a consumer
    // joining on Precursor.Id silently sums the target's transitions and its
    // decoy's into one trace. Measured 2026-08-18: a 1,000-precursor dump came
    // back with 24 fragments per precursor against 12 in the library, and every
    // trace statistic taken from it -- baseline, prominence, per-fragment apex
    // agreement -- was computed over target+decoy and was wrong.
    out.put("Precursor.Id\tDecoy\tTransition.Index\tProduct.Mz\tRT\tIntensity\n");

    const auto& p = library.precursors();
    const auto& t = library.transitions();
    for (std::size_t i = 0; i < library.precursorCount(); ++i)
    {
      // The library does not store Precursor.Id; DIA-NN's convention is
      // sequence + charge, and the writer reconstructs it the same way the
      // library writer does so the two files join on it.
      const auto seq = library.strings().get(p.modified_sequence[i]);
      const std::string id = std::string(seq) + std::to_string(static_cast<int>(p.charge[i]));
      // -out_chrom_ids. A listed id passes BOTH its target and its decoy block
      // through -- the Decoy column keeps them apart downstream, per the
      // warning at the top of this function.
      if (keep != nullptr && keep->count(id) == 0) { continue; }
      const char decoy_flag = p.decoy[i] ? '1' : '0';
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
          out.put(decoy_flag); out.put('\t');
          out.integer(tr); out.put('\t');
          out.number(product_mz, 6); out.put('\t');
          out.number(chromatograms.retentionTime(tr, j), 6); out.put('\t');
          out.number(chromatograms.intensity[b + j], 6); out.put('\n');
        }
      }
    }
    out.close();
  }

  std::unordered_set<std::string> readPrecursorIdList(const std::string& path)
  {
    std::unordered_set<std::string> ids;
    std::ifstream in(path);
    if (!in)
    { throw std::runtime_error("cannot read -out_chrom_ids file: " + path); }
    std::string line;
    while (std::getline(in, line))
    {
      if (!line.empty() && line.back() == '\r') { line.pop_back(); }
      const std::size_t tab = line.find('\t');
      std::string id = tab == std::string::npos ? line : line.substr(0, tab);
      if (id.empty() || id == "Precursor.Id") { continue; }
      ids.insert(std::move(id));
    }
    return ids;
  }

  void writeMs1TracesTsv(const std::string& path, const Library& library,
                         const Ms1Traces& iso0, const Ms1Traces& iso1,
                         const Ms1Traces& iso2,
                         const std::vector<std::uint32_t>& kept)
  {
    TextWriter out(path);
    out.put("Precursor.Id\tDecoy\tIsotope\tPrecursor.Mz\tRT\tIntensity\n");
    const auto& p = library.precursors();
    const Ms1Traces* iso[3] = {&iso0, &iso1, &iso2};
    // All three builds ran over the same MS1 grid and the same keep mask, so
    // row r means the same precursor in each. Verified here rather than
    // assumed: a silent mismatch would attribute one precursor's M+1 to
    // another's M.
    for (int k = 1; k < 3; ++k)
    {
      if (!iso[k]->empty() && iso[k]->bins() != iso0.bins())
      { throw std::runtime_error("MS1 isotope traces disagree on the time grid"); }
    }
    const auto& times = iso0.times();
    for (std::size_t r = 0; r < kept.size(); ++r)
    {
      const std::size_t i = kept[r];
      const auto seq = library.strings().get(p.modified_sequence[i]);
      const std::string id = std::string(seq) + std::to_string(static_cast<int>(p.charge[i]));
      const char decoy_flag = p.decoy[i] ? '1' : '0';
      const double mz = fromFixed(p.mz[i]);
      for (int k = 0; k < 3; ++k)
      {
        const Ms1Traces& t = *iso[k];
        if (t.empty()) { continue; }
        const std::size_t bins = t.bins();
        // Non-zero runs, padded with ONE flanking zero on each side: peak
        // boundaries survive, the empty remainder of the run does not. A dense
        // dump at cohort scale (138k rows x ~4k bins x 3 isotopes) would be
        // ~100 GB of zeros.
        std::size_t j = 0;
        std::size_t written_until = 0;   // exclusive; stops adjacent runs
                                         // double-writing their shared zero
        while (j < bins)
        {
          if (t.at(r, j) <= 0.0f) { ++j; continue; }
          std::size_t start = std::max(j > 0 ? j - 1 : 0, written_until);
          std::size_t end = j;
          while (end < bins && t.at(r, end) > 0.0f) { ++end; }
          const std::size_t stop = end < bins ? end + 1 : bins;  // one trailing zero
          for (std::size_t b = start; b < stop; ++b)
          {
            out.put(id); out.put('\t');
            out.put(decoy_flag); out.put('\t');
            out.integer(static_cast<std::uint32_t>(k)); out.put('\t');
            out.number(mz, 6); out.put('\t');
            out.number(times[b], 6); out.put('\t');
            out.number(t.at(r, b), 6); out.put('\n');
          }
          written_until = stop;
          j = stop;
        }
      }
    }
    out.close();
  }

} // namespace ODIA
