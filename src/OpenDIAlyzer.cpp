// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <odia/DIANNLibraryFile.h>
#include <odia/LibraryGenerator.h>
#include <odia/Library.h>

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace OpenMS;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page ODIA_OpenDIAlyzer OpenDIAlyzer

@brief Targeted DIA analysis.

The stages run in one process so that the feature table is built once and stays
flat. @p -stop_after ends the run early and writes that stage's output, which
makes each stage independently usable without a process boundary.

Currently implemented: the assay library stage.

<B>The command line parameters of this tool are:</B>
@verbinclude ODIA_OpenDIAlyzer.cli
*/

/// @cond TOPPCLASSES

class TOPPOpenDIAlyzer : public TOPPBase
{
public:
  TOPPOpenDIAlyzer() :
    TOPPBase("OpenDIAlyzer", "Targeted DIA analysis on OpenMS and mzPeak.", false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("tr", "<file>", "", "Assay library (DIA-NN TSV or Parquet).", false);
    registerInputFile_("fasta", "<file>", "", "Generate the library from these protein sequences instead.", false);
    setValidFormats_("fasta", {"fasta"}, false);

    registerStringOption_("decoys", "<method>", "mutate",
                          "Decoy construction when generating from FASTA.", false);
    setValidStrings_("decoys", {"mutate", "pseudo_reverse", "none"});

    registerIntOption_("missed_cleavages", "<n>", 1, "Maximum missed cleavages.", false, true);
    registerIntOption_("min_peptide_length", "<n>", 7, "Minimum peptide length.", false, true);
    registerIntOption_("max_peptide_length", "<n>", 30, "Maximum peptide length.", false, true);

    registerOutputFile_("out_lib", "<file>", "",
                        "Write the assay library here (DIA-NN TSV).", false);
    setValidFormats_("out_lib", {"tsv"}, false);

    registerStringOption_("stop_after", "<stage>", "",
                          "End the run after this stage and write its output.", false);
    setValidStrings_("stop_after", {"", "library"});

    registerFlag_("sort_library", "Sort precursors by m/z on load.", true);
  }

  ExitCodes main_(int, const char**) override
  {
    const std::string tr = getStringOption_("tr");
    const std::string fasta = getStringOption_("fasta");
    const std::string out_lib = getStringOption_("out_lib");
    const std::string stop_after = getStringOption_("stop_after");
    const bool sort_library = getFlag_("sort_library");

    if (tr.empty() == fasta.empty())
    {
      writeLogError_("Give exactly one of -tr <library> or -fasta <proteins>.");
      return ILLEGAL_PARAMETERS;
    }

    // Checked before any work is done. Doing it afterwards meant a run that
    // built and wrote a library still exited 6.
    if (stop_after != "library")
    {
      writeLogError_("Only the library stage is implemented; use -stop_after library.");
      return ILLEGAL_PARAMETERS;
    }

    ODIA::Library library;
    const auto t0 = std::chrono::steady_clock::now();
    try
    {
      if (!tr.empty())
      {
        ODIA::DIANNLibraryFile::load(tr, library);
      }
      else
      {
        ODIA::DigestParams params;
        params.missed_cleavages = static_cast<std::size_t>(getIntOption_("missed_cleavages"));
        params.min_length = static_cast<std::size_t>(getIntOption_("min_peptide_length"));
        params.max_length = static_cast<std::size_t>(getIntOption_("max_peptide_length"));
        params.decoy_method = ODIA::parseDecoyMethod(getStringOption_("decoys"));

        const auto stats = ODIA::LibraryGenerator::generate(fasta, params, library);
        std::size_t decoys_skipped = 0;
        const auto decoys = ODIA::LibraryGenerator::appendDecoys(
          library, params.decoy_method, &decoys_skipped);

        std::ostringstream gen;
        gen << "generated from " << stats.proteins << " proteins: "
            << stats.peptides << " peptides, " << stats.precursors << " target precursors, "
            << decoys << " decoys";
        if (decoys_skipped) { gen << " (" << decoys_skipped << " targets got none)"; }
        gen << "\n"
            << "  dropped: " << stats.dropped_precursor_mz << " outside the precursor m/z range, "
            << stats.dropped_too_few_fragments << " with too few fragments";
        writeLogInfo_(gen.str());
        writeLogWarn_("Fragment intensities and retention times are placeholders; "
                      "prediction is not wired up yet.");
      }
    }
    catch (const std::exception& e)
    {
      writeLogError_(std::string("Failed to build assay library: ") + e.what());
      return INPUT_FILE_CORRUPT;
    }
    const auto load_ms = std::chrono::duration<double, std::milli>(
                           std::chrono::steady_clock::now() - t0).count();

    if (sort_library) { library.sortByPrecursorMz(); }
    library.shrinkToFit();

    reportLibrary_(library, load_ms);

    if (!out_lib.empty())
    {
      try
      {
        ODIA::DIANNLibraryFile::storeTSV(out_lib, library);
      }
      catch (const std::exception& e)
      {
        writeLogError_(std::string("Failed to write assay library: ") + e.what());
        return CANNOT_WRITE_OUTPUT_FILE;
      }
    }

    return EXECUTION_OK;
  }

private:
  void reportLibrary_(const ODIA::Library& library, double load_ms)
  {
    const std::size_t decoys = library.decoyCount();
    const std::size_t targets = library.precursorCount() - decoys;

    std::ostringstream os;
    os << "assay library\n"
       << "  precursors:   " << library.precursorCount()
       << "  (" << targets << " target / " << decoys << " decoy)\n"
       << "  transitions:  " << library.transitionCount() << "\n";
    if (library.precursorCount())
    {
      os << std::fixed << std::setprecision(2)
         << "  per precursor: " << (static_cast<double>(library.transitionCount())
                                    / static_cast<double>(library.precursorCount()))
         << " transitions\n";
    }
    os << "  distinct strings: " << library.strings().size()
       << " (" << library.strings().bytes() << " bytes)\n"
       << std::fixed << std::setprecision(1)
       << "  in-memory:    " << (library.footprintBytes() / (1024.0 * 1024.0)) << " MiB";
    if (library.transitionCount())
    {
      os << "  (" << (static_cast<double>(library.footprintBytes())
                      / static_cast<double>(library.transitionCount()))
         << " bytes/transition)";
    }
    os << "\n  load time:    " << load_ms << " ms\n";

    // A target/decoy imbalance makes a 1% FDR unreachable by construction, and
    // seeing it here costs nothing compared with finding out after extraction.
    if (const auto bad = library.invalidMzCount(); bad != 0)
    {
      writeLogWarn_(std::to_string(bad) + " precursors have an unusable m/z "
                    "(empty, negative, NaN or out of range) and cannot match any "
                    "isolation window.");
    }

    if (decoys == 0)
    {
      writeLogWarn_("Library contains no decoys; target-decoy FDR will not be computable.");
    }
    else if (decoys > targets * 2 || targets > decoys * 2)
    {
      writeLogWarn_("Target and decoy counts are strongly imbalanced.");
    }

    writeLogInfo_(os.str());
  }
};

int main(int argc, const char** argv)
{
  TOPPOpenDIAlyzer tool;
  return tool.main(argc, argv);
}

/// @endcond
