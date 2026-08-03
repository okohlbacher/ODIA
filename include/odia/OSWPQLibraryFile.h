// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/Library.h>

#include <cstddef>
#include <string>

namespace ODIA
{

  /// Reads OpenMS's `.oswpq` assay-library bundle.
  ///
  /// `.oswpq` is an OpenMS format, not an ODIA one: it was introduced upstream
  /// in PR #9684 and is read by OpenSwathWorkflow and downstream OpenSWATH
  /// tooling. ODIA therefore reads what upstream writes rather than a private
  /// dialect, and does so without touching the OpenMS tree -- the bundle is
  /// opened here with libzip and Arrow directly, so no patched OpenMS is
  /// needed.
  class OSWPQLibraryFile
  {
  public:
    /// What the file gave up, beyond the library itself.
    ///
    /// These are reported rather than logged-and-forgotten because each one
    /// means peptides that will never be found, and a library that loads
    /// without complaint but half-empty is the failure mode this format has
    /// already produced once (a multi-chunk column read as chunk(0)).
    struct Stats
    {
      std::size_t precursor_rows = 0;
      std::size_t transition_rows = 0;

      /// Transitions whose precursor_id matches no precursor. They are dropped:
      /// a fragment with no precursor cannot be extracted.
      std::size_t orphan_transitions = 0;

      /// Precursors with no transitions at all. Kept, because the upstream
      /// writer emits them (the sample bundle has one) and dropping them would
      /// silently change the precursor count against the metadata census.
      std::size_t childless_precursors = 0;

      /// Precursor rows sharing an id with an earlier one. The id is the join
      /// key, so a repeat makes the join ambiguous: the first row wins and the
      /// later ones end up childless. Nothing else in the file would show this.
      std::size_t duplicate_precursor_ids = 0;

      /// Transitions whose own decoy flag contradicts their precursor's. The
      /// column is redundant, which makes it a free consistency check, and a
      /// target/decoy mix-up is exactly the corruption that produces a
      /// plausible-looking library and an unreachable FDR.
      std::size_t decoy_mismatches = 0;

      /// Rows whose charge did not fit the stored width, or whose ordinal was
      /// outside 1..255. Both are recorded as unknown rather than dropped: an
      /// OpenSWATH transition list is not obliged to carry an annotation, and
      /// the sample bundle carries none at all (type "", ordinal -1).
      std::size_t unusable_precursor_charges = 0;
      std::size_t unusable_transition_charges = 0;
      std::size_t unusable_ordinals = 0;

      /// Annotations naming a neutral loss this reader does not recognise.
      /// LossType::Other carries no mass, so such a fragment silently becomes
      /// one at the parent mass; without a count, a library written in a
      /// foreign annotation dialect loads with a clean stats block.
      std::size_t unrecognised_losses = 0;

      /// Whether library/metadata.json carried a census, and whether the rows
      /// actually read agree with it. A mismatch is the cheapest possible
      /// detector for a truncated read: it costs one comparison and it fires
      /// before an hour of extraction rather than after.
      bool census_present = false;
      bool census_agrees = false;
      std::size_t census_precursors = 0;
      std::size_t census_transitions = 0;

      int schema_version = 0;
      std::string generator;
      std::string openms_version;
    };

    /// Load the `library/` bundle of @p filename into @p library.
    ///
    /// Any pre-existing content of @p library is replaced.
    /// @param stats optional; receives what the file gave up.
    static void load(const std::string& filename, Library& library, Stats* stats = nullptr);

    /// True if @p filename is a ZIP archive carrying a `library/` bundle.
    static bool isLibraryBundle(const std::string& filename);
  };

} // namespace ODIA
