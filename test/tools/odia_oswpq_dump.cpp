// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Loads a .oswpq bundle and prints what ODIA made of it, so an independent
// reader (test/check_oswpq.py, which goes through pyarrow) can check it row by
// row rather than only comparing totals.

#include <odia/OSWPQLibraryFile.h>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::cerr << "usage: odia_oswpq_dump <bundle.oswpq>\n";
    return 2;
  }

  ODIA::Library library;
  ODIA::OSWPQLibraryFile::Stats stats;
  try
  {
    ODIA::OSWPQLibraryFile::load(argv[1], library, &stats);
  }
  catch (const std::exception& e)
  {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }

  std::printf("schema_version\t%d\n", stats.schema_version);
  std::printf("generator\t%s\n", stats.generator.c_str());
  std::printf("precursor_rows\t%zu\n", stats.precursor_rows);
  std::printf("transition_rows\t%zu\n", stats.transition_rows);
  std::printf("orphan_transitions\t%zu\n", stats.orphan_transitions);
  std::printf("childless_precursors\t%zu\n", stats.childless_precursors);
  std::printf("duplicate_precursor_ids\t%zu\n", stats.duplicate_precursor_ids);
  std::printf("decoy_mismatches\t%zu\n", stats.decoy_mismatches);
  std::printf("unusable_ordinals\t%zu\n", stats.unusable_ordinals);
  std::printf("unrecognised_losses\t%zu\n", stats.unrecognised_losses);
  std::printf("unusable_transition_charges\t%zu\n", stats.unusable_transition_charges);
  std::printf("unusable_precursor_charges\t%zu\n", stats.unusable_precursor_charges);
  std::printf("census_present\t%d\n", stats.census_present ? 1 : 0);
  std::printf("census_agrees\t%d\n", stats.census_agrees ? 1 : 0);
  std::printf("precursors\t%zu\n", library.precursorCount());
  std::printf("transitions\t%zu\n", library.transitionCount());
  std::printf("decoys\t%zu\n", library.decoyCount());
  std::printf("invalid_precursor_mz\t%zu\n", library.invalidMzCount());
  std::printf("invalid_transition_mz\t%zu\n", library.invalidMzTransitionCount());

  const auto& p = library.precursors();
  const auto& t = library.transitions();
  const auto& arena = library.strings();

  // Every precursor with every one of its fragments, in the order the CSR gives
  // them. Printing only totals would miss a scatter that puts the right number
  // of fragments under the wrong precursor.
  for (std::size_t i = 0; i < library.precursorCount(); ++i)
  {
    const auto seq = arena.get(p.modified_sequence[i]);
    const auto prot = arena.get(p.protein_group[i]);
    // NaN is printed as NaN rather than mapped to a number. Printing an absent
    // retention time as 0.0 made it indistinguishable from a real 0.0, so a
    // null read as zero -- which puts that precursor at the start of the RT
    // axis -- could not be detected by any test.
    char irt[32], im[32];
    if (std::isnan(p.irt[i])) { std::snprintf(irt, sizeof(irt), "NaN"); }
    else { std::snprintf(irt, sizeof(irt), "%.6g", double(p.irt[i])); }
    if (std::isnan(p.im[i])) { std::snprintf(im, sizeof(im), "NaN"); }
    else { std::snprintf(im, sizeof(im), "%.6g", double(p.im[i])); }
    std::printf("P\t%.*s\t%.*s\t%.5f\t%u\t%u\t%s\t%s\t%u\n",
                static_cast<int>(seq.size()), seq.data(),
                static_cast<int>(prot.size()), prot.data(),
                ODIA::fromFixed(p.mz[i]), unsigned(p.charge[i]), unsigned(p.decoy[i]),
                irt, im, p.transition_count[i]);
    for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
    {
      const std::size_t j = p.transition_begin[i] + k;
      const auto type = ODIA::toString(t.type[j]);
      const auto loss = ODIA::toString(t.loss[j]);
      std::printf("T\t%.5f\t%.6g\t%.*s\t%u\t%d\t%.*s\n",
                  ODIA::fromFixed(t.product_mz[j]), double(t.library_intensity[j]),
                  static_cast<int>(type.size()), type.data(),
                  unsigned(t.ordinal[j]), int(t.charge[j]),
                  static_cast<int>(loss.size()), loss.data());
    }
  }
  return 0;
}
