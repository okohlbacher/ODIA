// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <odia/DIANNLibraryFile.h>
#include <odia/LibraryGenerator.h>
#include <odia/Library.h>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
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
                          "Decoy construction. Applied to a library read with -tr "
                          "as well, if it has none already.", false);
    setValidStrings_("decoys", {"mutate", "pseudo_reverse", "none"});

    registerInputFile_("ms2_model", "<file>", "",
                       "PeptDeep MS2 fragment-intensity model. Defaults to the one OpenMS "
                       "downloads when built with WITH_ONNX=ON.", false);
    registerDoubleOption_("nce", "<energy>", 30.0,
                          "Normalised collision energy assumed for fragment-intensity "
                          "prediction. It changes the spectrum materially and nothing here "
                          "derives it from the data.", false);
    registerStringOption_("instrument", "<name>", "QE",
                          "Instrument assumed for fragment-intensity prediction. An "
                          "unrecognised name uses the model's 'unknown' slot rather than "
                          "silently predicting for a different instrument.", false);

    registerInputFile_("rt_model", "<file>", "",
                       "PeptDeep retention-time model. Defaults to the one OpenMS "
                       "downloads when built with WITH_ONNX=ON.", false);

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

        // Decoys for a supplied library too, not only for a generated one. A
        // library without them cannot be scored, and appendDecoys is idempotent
        // so one that already has them is left alone.
        const auto method = ODIA::parseDecoyMethod(getStringOption_("decoys"));
        if (method != ODIA::DecoyMethod::None && library.decoyCount() == 0)
        {
          std::size_t skipped = 0;
          const auto made = ODIA::LibraryGenerator::appendDecoys(library, method, &skipped);
          std::ostringstream msg;
          msg << "added " << made << " decoys";
          if (skipped) { msg << " (" << skipped << " targets got none)"; }
          writeLogInfo_(msg.str());
        }

        // Retention times are not predicted for a supplied library: it is
        // expected to carry its own. Say so if it does not, rather than writing
        // an empty column silently.
        std::size_t without_rt = 0;
        for (const auto v : library.precursors().irt)
        {
          if (std::isnan(v)) { ++without_rt; }
        }
        if (without_rt)
        {
          writeLogWarn_(std::to_string(without_rt) + " precursors in the supplied "
                        "library have no retention time; prediction is only applied "
                        "to libraries generated with -fasta.");
        }
      }
      else
      {
        ODIA::DigestParams params;
        params.missed_cleavages = static_cast<std::size_t>(getIntOption_("missed_cleavages"));
        params.min_length = static_cast<std::size_t>(getIntOption_("min_peptide_length"));
        params.max_length = static_cast<std::size_t>(getIntOption_("max_peptide_length"));
        params.decoy_method = ODIA::parseDecoyMethod(getStringOption_("decoys"));

        const auto stats = ODIA::LibraryGenerator::generate(fasta, params, library);

        std::ostringstream gen;
        gen << "generated from " << stats.proteins << " proteins: "
            << stats.peptides << " peptides, " << stats.precursors << " target precursors";
        gen << "\n"
            << "  dropped: " << stats.dropped_precursor_mz << " outside the precursor m/z range, "
            << stats.dropped_too_few_fragments << " with too few fragments";
        writeLogInfo_(gen.str());

        // Predict retention times, if a model is available. Fragment
        // intensities still need the MS2 model.
        // No absolute path baked into the binary. OpenMS installs the model it
        // downloads under its own share directory, so derive it from the
        // environment or let the user say.
        std::string rt_model = getStringOption_("rt_model");
        if (rt_model.empty())
        {
          if (const char* prefix = std::getenv("ODIA_OPENMS"))
          {
            const auto candidate = std::filesystem::path(prefix) /
              "share/OpenMS/models/peptdeep_rt_dynamic.onnx";
            if (std::filesystem::exists(candidate)) { rt_model = candidate.string(); }
          }
        }

        if (rt_model.empty())
        {
          writeLogWarn_("No retention-time model available; iRT is left unset. "
                        "Give one with -rt_model.");
        }
        else
        {
          // A prediction failure must not discard a library that is already
          // built: the unpredicted path warns and writes anyway, so this one
          // should too, rather than throwing the work away and reporting it as
          // a corrupt input.
          try
          {
            const auto t_rt = std::chrono::steady_clock::now();
            const auto unpredicted =
              ODIA::LibraryGenerator::predictRetentionTimes(library, rt_model);
            const auto rt_ms = std::chrono::duration<double, std::milli>(
                                 std::chrono::steady_clock::now() - t_rt).count();
            std::ostringstream rt;
            rt << "predicted retention times in " << rt_ms << " ms";
            if (unpredicted)
            {
              rt << "; " << unpredicted << " precursors left unpredicted";
            }
            writeLogInfo_(rt.str());
          }
          catch (const std::exception& e)
          {
            writeLogWarn_(std::string("retention-time prediction failed, iRT left "
                                      "unset: ") + e.what());
          }
        }

        // Fragment intensities, and with them the choice of which fragments to
        // keep. This runs before decoys are appended, because a decoy copies
        // its target's intensity pattern and would otherwise copy the
        // placeholder.
        std::string ms2_model = getStringOption_("ms2_model");
        if (ms2_model.empty())
        {
          if (const char* prefix = std::getenv("ODIA_OPENMS"))
          {
            const auto candidate = std::filesystem::path(prefix) /
              "share/OpenMS/models/peptdeep_ms2_dynamic.onnx";
            if (std::filesystem::exists(candidate)) { ms2_model = candidate.string(); }
          }
        }

        if (ms2_model.empty())
        {
          writeLogWarn_("No MS2 model available; fragment intensities stay as "
                        "placeholders and fragments are chosen by descending m/z. "
                        "Give one with -ms2_model.");
        }
        else
        {
          try
          {
            const auto t_ms2 = std::chrono::steady_clock::now();
            const auto unpredicted = ODIA::LibraryGenerator::predictFragmentIntensities(
              library, ms2_model, params,
              static_cast<float>(getDoubleOption_("nce")), getStringOption_("instrument"));
            const auto ms2_ms = std::chrono::duration<double, std::milli>(
                                  std::chrono::steady_clock::now() - t_ms2).count();
            std::ostringstream ms2;
            ms2 << "predicted fragment intensities in " << ms2_ms << " ms at NCE "
                << getDoubleOption_("nce") << " for " << getStringOption_("instrument");
            if (unpredicted)
            {
              ms2 << "; " << unpredicted << " precursors kept m/z-ranked placeholders";
            }
            writeLogInfo_(ms2.str());
          }
          catch (const std::exception& e)
          {
            writeLogWarn_(std::string("fragment-intensity prediction failed, "
                                      "placeholders kept: ") + e.what());
          }
        }

        std::size_t decoys_skipped = 0;
        const auto decoys = ODIA::LibraryGenerator::appendDecoys(
          library, params.decoy_method, &decoys_skipped, params.min_fragments);
        std::ostringstream dec;
        dec << "appended " << decoys << " decoys";
        if (decoys_skipped) { dec << " (" << decoys_skipped << " targets got none)"; }
        writeLogInfo_(dec.str());
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
    if (const auto bad = library.invalidMzTransitionCount(); bad != 0)
    {
      writeLogWarn_(std::to_string(bad) + " transitions have an unusable m/z and "
                    "cannot be extracted.");
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
