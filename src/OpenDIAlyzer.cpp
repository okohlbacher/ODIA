// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <odia/DIANNLibraryFile.h>
#include <odia/LibraryGenerator.h>
#include <odia/Library.h>
#include <fstream>
#include <odia/SpectrumSource.h>
#include <odia/ChromatogramExtractor.h>

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

    registerIntOption_("reserved_doubly_charged", "<n>", 0,
                       "Reserve this many of the fragment cap for doubly-charged ions. "
                       "0 ranks purely by predicted intensity, which is faithful to the "
                       "model; above 0 overrides it, and is only justified by a search.",
                       false);

    registerInputFile_("irt_standards", "<file>", "",
                       "Biognosys iRT standard peptides, for rescaling the RT model's "
                       "raw 0..1 output onto the iRT scale. Defaults to data/irt_standards.tsv "
                       "beside the tool. This changes units, not accuracy: the mapping is "
                       "monotone, so a consumer that fits its own RT calibration is "
                       "unaffected.", false);

    registerInputFile_("ccs_model", "<file>", "",
                       "PeptDeep collision-cross-section model. Defaults to the one OpenMS "
                       "downloads when built with WITH_ONNX=ON. Predicts CCS in square "
                       "angstroms; converting that to the 1/K0 an instrument reports is "
                       "deliberately left to the consumer, which knows the drift gas and "
                       "the calibration.", false);

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
    setValidStrings_("stop_after", {"", "library", "extract"});

    // No setValidFormats_ here on purpose. OpenMS has no mzPeak entry in
    // FileTypes, so declaring one makes TOPPBase try to identify the file and
    // fail with "ZIP archive contains 11 file entries; expected exactly 1" --
    // a message about the mzPeak container that names neither mzPeak nor the
    // tool's own option. The format is dispatched by openRun instead.
    registerInputFile_("in", "<file>", "",
                       "Run to extract from (mzPeak).", false);
    registerOutputFile_("out_chrom", "<file>", "",
                        "Write extracted chromatograms here (TSV).", false);
    setValidFormats_("out_chrom", {"tsv"}, false);

    registerDoubleOption_("fragment_ppm", "<ppm>", 20.0,
                          "Fragment mass tolerance for extraction.", false, true);
    registerDoubleOption_("rt_window", "<seconds>", 60.0,
                          "Half-width of the retention-time window around the "
                          "predicted elution.", false, true);
    registerDoubleOption_("irt_slope", "<a>", 0.0,
                          "Maps library iRT onto this run: rt = a * iRT + b. "
                          "0 spreads the library evenly over the run, which is a "
                          "placeholder, not a calibration.", false, true);
    registerDoubleOption_("irt_intercept", "<b>", 0.0, "See -irt_slope.", false, true);
    registerIntOption_("max_precursors", "<n>", 0,
                       "Extract only the first N precursors, 0 for all.", false, true);
    registerFlag_("no_ion_mobility",
                  "Ignore the ion-mobility dimension when matching windows.", true);

    registerFlag_("sort_library", "Sort precursors by m/z on load.", true);
  }


  /// Extract every transition of @p library from @p run, and write the
  /// chromatograms if asked.
  ///
  /// The iRT calibration is the part to be careful with. A library carries iRT,
  /// a run carries seconds, and nothing in either says how they map. Given a
  /// slope the caller trusts, that is used; given none, the library is spread
  /// evenly across the run's time range so that something extracts -- but that
  /// is a placeholder for a calibration, not one, and it says so out loud
  /// rather than producing quietly meaningless chromatograms.
  ExitCodes runExtraction_(const ODIA::Library& library, const std::string& run,
                           const std::string& out_chrom)
  {
    std::unique_ptr<ODIA::SpectrumSource> source;
    try
    {
      source = ODIA::openRun(run);
    }
    catch (const std::exception& e)
    {
      writeLogError_(std::string("Cannot open run ") + run + ": " + e.what());
      return INPUT_FILE_NOT_FOUND;
    }

    ODIA::ChromatogramExtractor::Options options;
    options.fragment_ppm = getDoubleOption_("fragment_ppm");
    options.rt_window_seconds = getDoubleOption_("rt_window");
    options.max_precursors = static_cast<std::size_t>(
      std::max(0, getIntOption_("max_precursors")));
    options.use_ion_mobility = !getFlag_("no_ion_mobility");
    options.irt_slope = getDoubleOption_("irt_slope");
    options.irt_intercept = getDoubleOption_("irt_intercept");
    options.threads = static_cast<unsigned>(std::max(1, getIntOption_("threads")));

    if (options.irt_slope == 0.0)
    {
      writeLogWarn_("No iRT calibration given (-irt_slope/-irt_intercept). The "
                    "library is being spread evenly over the run, which will "
                    "extract from approximately the wrong retention times. "
                    "Treat the output as a smoke test, not a result.");
    }

    ODIA::ChromatogramExtractor::Stats stats;
    ODIA::Chromatograms chromatograms;
    const auto t = std::chrono::steady_clock::now();
    try
    {
      chromatograms = ODIA::ChromatogramExtractor::extract(library, *source, options, &stats);
    }
    catch (const std::exception& e)
    {
      writeLogError_(std::string("Extraction failed: ") + e.what());
      return INPUT_FILE_CORRUPT;
    }
    const auto ms = std::chrono::duration<double, std::milli>(
                      std::chrono::steady_clock::now() - t).count();

    std::ostringstream msg;
    msg << "extracted " << stats.transitions << " transitions of "
        << stats.precursors << " precursors from " << stats.spectra_read
        << " spectra in " << ms << " ms\n"
        << "  decode " << stats.decode_seconds << " s, index "
        << stats.index_seconds << " s, match " << stats.match_seconds << " s\n"
        << "  points " << stats.points << " (" << stats.nonzero_points
        << " nonzero), " << chromatograms.footprintBytes() / 1048576.0 << " MiB";
    if (stats.outside_rt_range)
    {
      msg << "\n  " << stats.outside_rt_range
          << " precursors predicted to elute outside the run";
    }
    if (chromatograms.precursors_without_window)
    {
      msg << "\n  " << chromatograms.precursors_without_window
          << " precursors covered by no isolation window";
    }
    writeLogInfo_(msg.str());

    if (!out_chrom.empty())
    {
      try
      {
        writeChromatograms_(out_chrom, library, chromatograms);
      }
      catch (const std::exception& e)
      {
        writeLogError_(std::string("Failed to write chromatograms: ") + e.what());
        return CANNOT_WRITE_OUTPUT_FILE;
      }
      writeLogInfo_("wrote chromatograms to " + out_chrom);
    }
    return EXECUTION_OK;
  }

  /// Long format, one row per point: transition, retention time, intensity.
  ///
  /// Deliberately not one row per transition with packed arrays. This file is
  /// what the scoring stage and any external check will read, and a long table
  /// is what every tool that might read it -- pandas, R, DuckDB -- takes
  /// without a parser of its own.
  static void writeChromatograms_(const std::string& path, const ODIA::Library& library,
                                  const ODIA::Chromatograms& chromatograms)
  {
    std::ofstream out(path);
    if (!out) { throw std::runtime_error("cannot open " + path); }
    out << "Precursor.Id\tTransition.Index\tProduct.Mz\tRT\tIntensity\n";

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
        const std::uint32_t b = chromatograms.begin[tr];
        const std::uint32_t n = chromatograms.count[tr];
        for (std::uint32_t j = 0; j < n; ++j)
        {
          out << id << '\t' << tr << '\t' << ODIA::fromFixed(t.product_mz[tr]) << '\t'
              << chromatograms.retention_time[b + j] << '\t'
              << chromatograms.intensity[b + j] << '\n';
        }
      }
    }
    if (!out) { throw std::runtime_error("write failed for " + path); }
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

    const std::string in_run = getStringOption_("in");
    const std::string out_chrom = getStringOption_("out_chrom");

    // Checked before any work is done. Doing it afterwards meant a run that
    // built and wrote a library still exited 6.
    if (stop_after != "library" && stop_after != "extract")
    {
      writeLogError_("Implemented stages are 'library' and 'extract'; "
                     "use -stop_after library or -stop_after extract.");
      return ILLEGAL_PARAMETERS;
    }
    if (stop_after == "extract" && in_run.empty())
    {
      writeLogError_("-stop_after extract needs a run to extract from: give -in <file>.");
      return ILLEGAL_PARAMETERS;
    }

    ODIA::Library library;

    // -threads is a TOPPBase option, and until now it reached digestion and
    // decoy construction only. Inference is ~93% of this stage, so leaving it
    // to ONNX Runtime's own heuristic meant the flag governed almost nothing.
    const auto inference_sessions =
      static_cast<unsigned>(std::max(1, getIntOption_("threads")));

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
        params.reserved_doubly_charged =
          static_cast<std::size_t>(getIntOption_("reserved_doubly_charged"));

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
              ODIA::LibraryGenerator::predictRetentionTimes(library, rt_model, true,
                                                           inference_sessions);
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
              static_cast<float>(getDoubleOption_("nce")), getStringOption_("instrument"),
              true, inference_sessions);
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

        // Rescale the raw retention times onto the iRT scale, before decoys so
        // a decoy inherits a calibrated value like everything else.
        //
        // Units, not accuracy. Measured on the human proteome: DIA-NN's search
        // window was 2.18905 min with and without this, identical to the last
        // digit, because DIA-NN fits its own monotone calibration. It is here
        // so that a column named iRT holds an iRT, which matters for any
        // consumer that applies a tolerance in those units without calibrating.
        if (!rt_model.empty())
        {
          std::string standards = getStringOption_("irt_standards");
          if (standards.empty())
          {
            for (const auto& candidate :
                 {std::filesystem::path("data/irt_standards.tsv"),
                  std::filesystem::path(ODIA_DATA_DIR) / "irt_standards.tsv"})
            {
              if (std::filesystem::exists(candidate)) { standards = candidate.string(); break; }
            }
          }
          if (standards.empty())
          {
            writeLogWarn_("No iRT standards available; retention times are left on the "
                          "model's raw 0..1 scale, which is NOT iRT. Give a file with "
                          "-irt_standards.");
          }
          else
          {
            try
            {
              const auto cal =
                ODIA::LibraryGenerator::fitIrtCalibration(rt_model, standards);
              ODIA::LibraryGenerator::applyIrtCalibration(library, cal);
              std::ostringstream msg;
              msg << "rescaled retention times to iRT: " << cal.slope << " * raw + "
                  << cal.intercept << " from " << cal.peptides << " standards"
                  << " (worst standard off by " << cal.max_abs_error << " iRT)";
              writeLogInfo_(msg.str());
            }
            catch (const std::exception& e)
            {
              writeLogWarn_(std::string("iRT calibration failed; retention times stay on "
                                        "the raw 0..1 scale: ") + e.what());
            }
          }
        }

        // Collision cross-sections, before decoys so a decoy inherits its
        // target's value the way it inherits iRT and the intensity pattern.
        std::string ccs_model = getStringOption_("ccs_model");
        if (ccs_model.empty())
        {
          if (const char* prefix = std::getenv("ODIA_OPENMS"))
          {
            const auto candidate = std::filesystem::path(prefix) /
              "share/OpenMS/models/peptdeep_ccs_dynamic.onnx";
            if (std::filesystem::exists(candidate)) { ccs_model = candidate.string(); }
          }
        }
        if (ccs_model.empty())
        {
          writeLogWarn_("No CCS model available; the cross-section column is left "
                        "empty. Give one with -ccs_model.");
        }
        else
        {
          try
          {
            const auto t_ccs = std::chrono::steady_clock::now();
            const auto unpredicted =
              ODIA::LibraryGenerator::predictCollisionCrossSections(library, ccs_model, true,
                                                                   inference_sessions);
            const auto ccs_ms = std::chrono::duration<double, std::milli>(
                                  std::chrono::steady_clock::now() - t_ccs).count();
            std::ostringstream ccs;
            ccs << "predicted collision cross-sections in " << ccs_ms << " ms";
            if (unpredicted) { ccs << "; " << unpredicted << " left unpredicted"; }
            writeLogInfo_(ccs.str());
          }
          catch (const std::exception& e)
          {
            writeLogWarn_(std::string("CCS prediction failed, the column is left "
                                      "empty: ") + e.what());
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

    if (stop_after == "extract")
    {
      const auto rc = runExtraction_(library, in_run, out_chrom);
      if (rc != EXECUTION_OK) { return rc; }
    }

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
