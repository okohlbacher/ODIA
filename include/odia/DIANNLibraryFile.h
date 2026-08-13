// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/Library.h>

#include <string>

namespace ODIA
{

  /// Reads and writes assay libraries in DIA-NN's dialect.
  ///
  /// DIA-NN's TSV and Parquet libraries share column names, so one column map
  /// serves both. Both are fully denormalised -- one row per transition, with
  /// every precursor field repeated -- which is exactly the redundancy the CSR
  /// layout in Library collapses on load.
  class DIANNLibraryFile
  {
  public:
    /// Column names, as DIA-NN writes them.
    struct Columns
    {
      static constexpr const char* PRECURSOR_ID = "Precursor.Id";
      static constexpr const char* MODIFIED_SEQUENCE = "Modified.Sequence";
      static constexpr const char* STRIPPED_SEQUENCE = "Stripped.Sequence";
      static constexpr const char* PRECURSOR_CHARGE = "Precursor.Charge";
      static constexpr const char* DECOY = "Decoy";
      static constexpr const char* RT = "RT";
      static constexpr const char* IM = "IM";
      /// Predicted collision cross-section, square angstroms.
      ///
      /// Its own column, not IM. DIA-NN's IM is 1/K0; converting to it needs
      /// the Mason-Schamp relation with the drift gas and the instrument's
      /// calibration, which is done downstream where the instrument is known.
      /// A consumer that wants ion mobility must convert this, deliberately.
      static constexpr const char* CCS = "CCS";
      static constexpr const char* PRECURSOR_MZ = "Precursor.Mz";
      static constexpr const char* PRODUCT_MZ = "Product.Mz";
      static constexpr const char* RELATIVE_INTENSITY = "Relative.Intensity";
      static constexpr const char* FRAGMENT_TYPE = "Fragment.Type";
      static constexpr const char* FRAGMENT_CHARGE = "Fragment.Charge";
      static constexpr const char* FRAGMENT_SERIES_NUMBER = "Fragment.Series.Number";
      static constexpr const char* FRAGMENT_LOSS_TYPE = "Fragment.Loss.Type";
      static constexpr const char* PROTEIN_GROUP = "Protein.Group";
    };

    /// Load a library. Dispatches on the extension: .parquet, else TSV.
    ///
    /// Rows are grouped into precursors by a change in Precursor.Id, so the
    /// input must keep a precursor's transitions together -- which DIA-NN does.
    static void load(const std::string& filename, Library& library);

    static void loadTSV(const std::string& filename, Library& library);
    static void loadParquet(const std::string& filename, Library& library);

    /// Write in DIA-NN's TSV dialect.
    ///
    /// Precursor.Id is synthesised as <Modified.Sequence><charge>, matching
    /// DIA-NN, since ODIA does not store per-transition identifiers.
    static void storeTSV(const std::string& filename, const Library& library);

    /// Write the library as Parquet. Same columns as `storeTSV`, same names,
    /// so `loadParquet` reads back what this writes.
    ///
    /// Worth having because the TSV is one row per TRANSITION with the sequence
    /// and protein group repeated 12 times: a proteome library is 12.35 GiB of
    /// text that takes 67-86 s to re-parse on every single search. Parquet
    /// dictionary-encodes exactly those repeated strings.
    static void storeParquet(const std::string& filename, const Library& library);

    /// Dispatches on the extension, mirroring `load`.
    static void store(const std::string& filename, const Library& library);
  };

} // namespace ODIA
