// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/FragvecTsv.h>
#include <odia/TextWriter.h>

#include <cstddef>
#include <stdexcept>
#include <string>

namespace ODIA
{

  void writeFragvecTsv(const std::string& path, const Library& library,
                       const PeakGroupScorer::Result& scored)
  {
    const std::size_t stride = PeakGroupScorer::N_FRAGVEC;
    if (scored.fragvec.size() != scored.groups.size() * stride)
    {
      throw std::runtime_error(
        "-out_fragvec: " + std::to_string(scored.fragvec.size()) +
        " values for " + std::to_string(scored.groups.size()) +
        " groups; expected " + std::to_string(stride) + " per group");
    }

    TextWriter out(path);
    out.put("Precursor.Id\tDecoy\tOrdinal");
    for (const auto& n : PeakGroupScorer::fragvecNames()) { out.put('\t'); out.put(n); }
    out.put('\n');

    const auto& p = library.precursors();
    // Precursor.Id is reconstructed exactly as writeScores_ and
    // writeChromatogramTsv reconstruct it -- modified sequence + charge, so a
    // decoy carries its target's id and the Decoy column separates them. Any
    // other spelling here would not join to -out.
    std::string prev_id;
    int prev_decoy = -1;
    long long ordinal = -1;
    for (std::size_t r = 0; r < scored.groups.size(); ++r)
    {
      const auto& g = scored.groups[r];
      const auto seq = library.strings().get(p.modified_sequence[g.precursor]);
      std::string id =
        std::string(seq) + std::to_string(static_cast<int>(p.charge[g.precursor]));
      const int decoy = static_cast<int>(g.decoy);
      // Reset on the (id, Decoy) pair rather than on the library precursor
      // index: two library entries can reconstruct the same id, and the
      // consumer's block is the id, not the index.
      if (decoy != prev_decoy || id != prev_id)
      {
        ordinal = 0; prev_decoy = decoy; prev_id = id;
      }
      else { ++ordinal; }
      out.put(id); out.put('\t');
      out.integer(decoy); out.put('\t');
      out.integer(ordinal);
      const float* row = scored.fragvec.data() + r * stride;
      for (std::size_t j = 0; j < stride; ++j)
      {
        out.put('\t');
        // NINE, not the 6 -out and -out_chrom use. The contract stores these as
        // float32 and the arm's gate A4 compares them cell for cell against the
        // reference builder's float32, so the file has to round-trip a float32
        // exactly -- and %.6g does not: MEASURED on the wf_v64 fixture, 5,164
        // rows x 78 columns came back EXACTLY the reference rounded through
        // %.6g, which cost up to 4.96e-05 absolute on R1_LOGAREA/R1_LOGTOT
        // (values near 10) and 5.01e-06 on R1_LOGRATIO -- a uniform ~5e-06
        // RELATIVE error, i.e. six significant digits and nothing else.
        // FLT_DECIMAL_DIG is 9, the shortest precision that round-trips every
        // float32, so at 9 the gate is an equality rather than a tolerance a
        // real defect could hide inside. It is a deliberate DEVIATION, not a
        // free choice: run_x1.sh's registered gate A4 asks only for |d| <= 1e-4
        // relative off the six exact indicator/count families, and %.6g's 5e-06
        // sits inside that by 20x -- precision 9 buys an equality gate at the
        // cost of file size. That cost is MEASURED on the wf_v64 fixture, not
        // estimated: 400.4 -> 497.2 B/row (+24%), i.e. 10.86 GB at 21,849,969
        // rows, inside the 10-16 GB wf_v63_check.txt s.4.8 budgeted.
        // `number` takes a double because that is what an ostream would have
        // formatted; rendering the float directly is a different set of digits.
        out.number(static_cast<double>(row[j]), 9);
      }
      out.put('\n');
    }
    out.close();
  }

} // namespace ODIA
