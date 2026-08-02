// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/Library.h>

#include <cstddef>
#include <string>
#include <vector>

namespace ODIA
{

  /// How decoy precursors are constructed.
  ///
  /// Kept selectable rather than fixed. DIA-NN's mutation scheme is the default
  /// (D7), but its own documentation identifies it as the weak point of the FDR
  /// model -- the decoy keeps the target's precursor m/z, iRT and library
  /// intensity pattern, so it is searched in the same isolation window over the
  /// same RT range with the same expected shape, and only fragment masses
  /// differ. Phase 3 needs to be able to measure that against the alternative
  /// under entrapment rather than inherit it.
  enum class DecoyMethod
  {
    None,
    Mutate,        ///< DIA-NN: substitute one residue near each terminus
    PseudoReverse  ///< reverse all but the C-terminal residue
  };

  DecoyMethod parseDecoyMethod(const std::string& s);

  struct DigestParams
  {
    std::string enzyme = "Trypsin";
    std::size_t missed_cleavages = 1;
    std::size_t min_length = 7;
    std::size_t max_length = 30;

    std::vector<int> charges{2, 3};

    /// UniMod-style names as OpenMS knows them, e.g. "Carbamidomethyl (C)".
    std::vector<std::string> fixed_modifications{"Carbamidomethyl (C)"};
    std::vector<std::string> variable_modifications{};
    std::size_t max_variable_modifications = 1;

    double precursor_mz_min = 350.0;
    double precursor_mz_max = 1200.0;

    double fragment_mz_min = 200.0;
    double fragment_mz_max = 1800.0;
    int max_fragment_charge = 2;

    /// A precursor with fewer usable fragments than this is not searchable.
    std::size_t min_fragments = 4;
    /// Cap per precursor, taken in descending fragment m/z.
    std::size_t max_fragments = 12;

    DecoyMethod decoy_method = DecoyMethod::Mutate;
  };

  /// Builds an assay library from protein sequences.
  ///
  /// Produces the peptide, precursor and fragment skeleton with placeholder
  /// intensities and retention times; those are filled in by prediction.
  class LibraryGenerator
  {
  public:
    struct Stats
    {
      std::size_t proteins = 0;
      std::size_t peptides = 0;
      std::size_t precursors = 0;
      std::size_t transitions = 0;
      std::size_t dropped_too_few_fragments = 0;
      std::size_t dropped_precursor_mz = 0;
    };

    static Stats generate(const std::string& fasta_file,
                          const DigestParams& params,
                          Library& library);

    /// Append a decoy for every target currently in @p library.
    ///
    /// Decoys keep the target's precursor m/z, iRT and intensity pattern; only
    /// fragment m/z values move. See DecoyMethod.
    /// @param skipped_out receives the number of targets for which no decoy
    ///        could be built. Reported rather than dropped silently: 1.3% of the
    ///        DIA-NN fixture is affected, all of them N-terminally modified.
    /// @param min_fragments the same bar the targets had to clear. Applying it
    ///        to one class only is an anti-conservative FDR (D7 rule 2).
    /// Idempotent: targets that already have a decoy are skipped.
    static std::size_t appendDecoys(Library& library, DecoyMethod method,
                                    std::size_t* skipped_out = nullptr,
                                    std::size_t min_fragments = 0);
  };

} // namespace ODIA
