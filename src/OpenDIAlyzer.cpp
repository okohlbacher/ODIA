// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <odia/DIANNLibraryFile.h>
#include <odia/LibraryGenerator.h>
#include <odia/Library.h>
#include <fstream>
#include <odia/SpectrumSource.h>
#include <odia/ChromatogramExtractor.h>
#include <odia/MassCalibration.h>
#include <odia/PeakGroupScorer.h>
#include <odia/RtCalibration.h>

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
    setValidStrings_("stop_after", {"", "library", "extract", "score"});

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

    // NO DEFAULT. -1 means "not given", and the run's own calibration then
    // supplies it. Naming a number up front is a guess about the instrument, and
    // a wrong guess is expensive in both directions -- too wide admits
    // interference, too narrow discards real fragments before anything can score
    // them. Measured on S08: target-minus-decoy fragment presence falls
    // monotonically with tolerance (0.083 at 10 ppm, 0.051 at 20, 0.034 at 30,
    // 0.019 at 50), and +/-10 CENTRED on the fitted -9.8 measured x1.13 overall
    // and x1.24 in the weakest abundance decile -- while +/-10 about zero keeps
    // only 0.51 of true fragments against 0.78 for +/-15. So the right width
    // depends on whether the centring succeeded, which is not knowable here.
    //
    // The "calibration may only narrow" clamp therefore applies ONLY when a
    // value was given explicitly. Unset means "tell me what this run says", and
    // clamping that against a number nobody chose would silently reinstate a
    // default.
    registerDoubleOption_("fragment_ppm", "<ppm>", -1.0,
                          "Fragment mass tolerance, ppm half-width. UNSET by default: 10 when the "
                          "run's own mass calibration centres the window, 15 when it cannot. Give "
                          "a value to pin it, in which case calibration may only narrow from "
                          "there.", false, true);
    registerStringOption_("mass_calibration", "<mode>", "auto",
                          "auto: measure the run's systematic fragment mass error before "
                          "extracting and centre the window on it. off: extract uncalibrated, "
                          "which on this instrument class discards about half the fragment "
                          "evidence.", false);
    setValidStrings_("mass_calibration", {"auto", "off"});
    registerDoubleOption_("fragment_ppm_offset", "<ppm>", 0.0,
                          "Pin the systematic fragment mass offset instead of measuring it. "
                          "Non-zero also switches -mass_calibration off, because a measured "
                          "value and a pinned one cannot both be applied.", false, true);
    registerIntOption_("mz_calib_precursors", "<n>", 3000,
                       "Precursors sampled when measuring the mass error.", false, true);
    registerIntOption_("mz_calib_cycles", "<n>", 160,
                       "Acquisition cycles probed, spread over the gradient. This is what the "
                       "measurement costs: one cycle is one decoded spectrum per isolation "
                       "window.", false, true);
    registerDoubleOption_("mz_calib_search_ppm", "<ppm>", 50.0,
                          "Half-width searched while COLLECTING the residuals, before any window "
                          "is inferred. Deliberately far wider than anything extracted with: the "
                          "distribution's shoulders have to be visible. 0 disables inference.",
                          false, true);
    registerDoubleOption_("precursor_im_window", "<1/K0>", 0.025,
                          "Half-width of the ion-mobility window around the PRECURSOR's own "
                          "library 1/K0. 0 disables it, leaving only the isolation window's "
                          "band, which is ~8x wider than a precursor occupies.", false, true);
    registerStringOption_("aggregate", "<how>", "sum",
                          "How several peaks inside one transition's tolerance box become one "
                          "number. sum integrates; max takes the largest, which returns the "
                          "interference envelope once the box spans the mobility axis.",
                          false, true);
    setValidStrings_("aggregate", {"sum", "max"});
    registerDoubleOption_("rt_window", "<seconds>", 600.0,
                          "Half-width of the retention-time window around the "
                          "predicted elution. 600 s matches OpenSWATH's second "
                          "pass; 60 s was below our own measured calibration "
                          "residual (76.6 s SD), which put the true peak outside "
                          "the window for a third of precursors.", false, true);
    registerDoubleOption_("irt_slope", "<a>", 0.0,
                          "Maps library iRT onto this run: rt = a * iRT + b. "
                          "0 spreads the library evenly over the run, which is a "
                          "placeholder, not a calibration.", false, true);
    registerDoubleOption_("irt_intercept", "<b>", 0.0, "See -irt_slope.", false, true);
    registerIntOption_("max_precursors", "<n>", 0,
                       "Extract only the first N precursors, 0 for all.", false, true);
    registerOutputFile_("out", "<file>", "",
                        "Write scored peak groups here (TSV).", false);
    setValidFormats_("out", {"tsv"}, false);
    registerStringOption_("classifier", "<name>", "gbt",
                          "Discriminant for the semi-supervised scorer.", false, true);
    setValidStrings_("classifier", {"lda", "gbt", "nn"});
    registerIntOption_("passes", "<n>", 2,
                       "1 extracts once with the given calibration. 2 extracts wide, "
                       "scores, fits the retention-time map from the confident "
                       "identifications, and re-extracts narrow.", false);
    registerDoubleOption_("rt_window_pass1", "<seconds>", 0.0,
                          "Half-width for the anchor-collecting pass. Deliberately "
                          "wide: a predicted library's RT error is far larger than "
                          "the window pass 2 uses, so a narrow first pass finds "
                          "nothing to calibrate from. 0 means the whole run.",
                          false, true);
    registerDoubleOption_("anchor_q", "<q>", 0.05,
                          "q-value below which a pass-1 identification is used as a "
                          "calibration anchor. Lenient on purpose -- pass 1 is "
                          "uncalibrated, so demanding 1% there yields no anchors and "
                          "no second pass.", false, true);
    registerIntOption_("min_anchors", "<n>", 20,
                       "Below this many anchors the fit is not attempted and the run "
                       "says so rather than calibrating from noise.", false, true);
    registerIntOption_("max_candidates", "<n>", 3,
                       "Candidate peak groups kept per precursor. More than one on "
                       "purpose: keeping only the best hides the true peak whenever "
                       "it ranks second, and leaves the decoys nothing to be wrong "
                       "about, which deflates the FDR.", false, true);

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
                           const std::string& out_chrom,
                           ODIA::Chromatograms* keep = nullptr,
                           double rt_window_override = 0.0,
                           bool library_rt_is_run_seconds = false)
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
    applyMassCalibration_(library, *source, options);
    options.rt_window_seconds = rt_window_override != 0.0 ? rt_window_override
                                                          : getDoubleOption_("rt_window");
    options.max_precursors = static_cast<std::size_t>(
      std::max(0, getIntOption_("max_precursors")));
    options.use_ion_mobility = !getFlag_("no_ion_mobility");
    options.precursor_im_window = getDoubleOption_("precursor_im_window");
    options.aggregate = getStringOption_("aggregate") == "max"
                          ? ODIA::ChromatogramExtractor::Options::Aggregate::Max
                          : ODIA::ChromatogramExtractor::Options::Aggregate::Sum;
    options.irt_slope = library_rt_is_run_seconds ? 1.0 : getDoubleOption_("irt_slope");
    options.irt_intercept = library_rt_is_run_seconds ? 0.0 : getDoubleOption_("irt_intercept");
    options.threads = static_cast<unsigned>(std::max(1, getIntOption_("threads")));

    if (options.irt_slope == 0.0 && !library_rt_is_run_seconds)
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
    if (keep != nullptr) { *keep = std::move(chromatograms); }
    return EXECUTION_OK;
  }


  /// Extract, score, calibrate from what was confidently identified, and do it
  /// again on the corrected retention-time axis.
  ///
  /// Why two passes rather than one. A predicted library's retention times are
  /// wrong by far more than the window a second pass can afford: the residual
  /// is hundreds of seconds where the window is tens. Extracting narrow from an
  /// uncalibrated library therefore samples the wrong part of the run and finds
  /// nothing -- measured here as 0 identifications with the classifier fitting
  /// correctly and the FDR correctly refusing to call noise. So pass 1 extracts
  /// WIDE purely to find anchors, and pass 2 extracts narrow where they say.
  ExitCodes runScoreWorkflow_(ODIA::Library& library, const std::string& run,
                              const std::string& out_chrom, const std::string& out)
  {
    const int passes = std::max(1, getIntOption_("passes"));
    ODIA::Chromatograms chromatograms;

    if (passes == 1)
    {
      const auto rc = runExtraction_(library, run, out_chrom, &chromatograms);
      if (rc != EXECUTION_OK) { return rc; }
      return runScoring_(library, chromatograms, out);
    }

    // The library's own retention times, kept before anything is applied to
    // them. The fit maps library RT -> run RT, so it must always be evaluated
    // on the ORIGINAL values; applying it to already-transformed ones composes
    // the passes and puts pass 2's windows nowhere.
    const std::vector<float> original_irt = library.precursors().irt;

    writeLogInfo_("pass 1 of 2: wide extraction to collect calibration anchors");
    const double pass1_window = getDoubleOption_("rt_window_pass1");
    {
      // 0 means "the whole run", expressed as a window wider than any gradient
      // rather than as a sentinel the extractor would have to know about. A
      // negative value would simply extract nothing.
      const auto rc = runExtraction_(library, run, "", &chromatograms,
                                     pass1_window > 0.0 ? pass1_window : 1.0e9);
      if (rc != EXECUTION_OK) { return rc; }
    }

    ODIA::PeakGroupScorer::Options options;
    options.classifier = getStringOption_("classifier");
    options.max_candidates = static_cast<std::size_t>(std::max(1, getIntOption_("max_candidates")));
    options.threads = static_cast<unsigned>(std::max(1, getIntOption_("threads")));
    const auto pass1 = ODIA::PeakGroupScorer::score(library, chromatograms, options);

    std::ostringstream p1;
    p1 << "pass 1: " << pass1.groups.size() << " peak groups, "
       << pass1.target_groups << " target / " << pass1.decoy_groups << " decoy, "
       << pass1.iterations_trained << " iterations trained";
    writeLogInfo_(p1.str());

    // Anchors: the best group of each confidently identified target, paired
    // with the library RT it came from.
    const double anchor_q = getDoubleOption_("anchor_q");
    std::vector<std::pair<double, double>> anchors;
    if (pass1.fdr_valid)
    {
      std::vector<const ODIA::PeakGroupScorer::PeakGroup*> best(
        library.precursorCount(), nullptr);
      for (const auto& g : pass1.groups)
      {
        if (g.decoy || g.qvalue > anchor_q) { continue; }
        auto*& b = best[g.precursor];
        if (b == nullptr || g.qvalue < b->qvalue) { b = &g; }
      }
      for (std::size_t i = 0; i < best.size(); ++i)
      {
        if (best[i] != nullptr && std::isfinite(original_irt[i]))
        {
          anchors.emplace_back(static_cast<double>(original_irt[i]),
                               static_cast<double>(best[i]->apex_rt));
        }
      }
    }

    const int min_anchors = std::max(1, getIntOption_("min_anchors"));
    if (static_cast<int>(anchors.size()) < min_anchors)
    {
      std::ostringstream why;
      why << "pass 1 yielded " << anchors.size() << " calibration anchors, fewer than "
          << min_anchors << ". Not fitting a retention-time map from that -- a "
          << "calibration from a handful of uncertain anchors is worse than none, "
          << "because pass 2 would extract narrow around it and find nothing. ";
      why << (pass1.fdr_valid
                ? "Widen -rt_window_pass1, or relax -anchor_q."
                : "Pass 1 produced no usable FDR at all; see the warning above.");
      writeLogWarn_(why.str());
      return runScoring_(library, chromatograms, out);
    }

    double p95 = 0.0;
    const auto trafo = ODIA::Calibration::fit(anchors, &p95, 0.0);
    std::ostringstream fit;
    fit << "fitted the retention-time map from " << anchors.size()
        << " anchors; p95 residual " << p95 << " s";
    writeLogInfo_(fit.str());

    // Applied to the ORIGINAL values, for the reason above.
    auto& irt = library.precursors().irt;
    for (std::size_t i = 0; i < irt.size(); ++i)
    {
      if (std::isfinite(original_irt[i]))
      {
        irt[i] = static_cast<float>(trafo.apply(static_cast<double>(original_irt[i])));
      }
    }

    writeLogInfo_("pass 2 of 2: narrow extraction on the calibrated axis");
    chromatograms = ODIA::Chromatograms{};
    {
      // The library now carries run seconds, so the affine map is the identity.
      const auto rc = runExtraction_(library, run, out_chrom, &chromatograms, 0.0, true);
      if (rc != EXECUTION_OK) { return rc; }
    }
    return runScoring_(library, chromatograms, out);
  }

  /// Find peak groups, score them, and write them with their q-values.
  ExitCodes runScoring_(const ODIA::Library& library,
                        const ODIA::Chromatograms& chromatograms,
                        const std::string& out)
  {
    ODIA::PeakGroupScorer::Options options;
    options.classifier = getStringOption_("classifier");
    options.max_candidates = static_cast<std::size_t>(
      std::max(1, getIntOption_("max_candidates")));
    options.threads = static_cast<unsigned>(std::max(1, getIntOption_("threads")));

    const auto t = std::chrono::steady_clock::now();
    ODIA::PeakGroupScorer::Result scored;
    try
    {
      scored = ODIA::PeakGroupScorer::score(library, chromatograms, options);
    }
    catch (const std::exception& e)
    {
      writeLogError_(std::string("Scoring failed: ") + e.what());
      return INTERNAL_ERROR;
    }
    const auto ms = std::chrono::duration<double, std::milli>(
                      std::chrono::steady_clock::now() - t).count();

    std::ostringstream msg;
    msg << "scored " << scored.groups.size() << " peak groups with "
        << options.classifier << " in " << ms << " ms\n"
        << "  " << scored.target_groups << " target / " << scored.decoy_groups
        << " decoy groups\n";
    if (scored.fdr_valid)
    {
      msg << "  identified " << scored.identified_at_1pct << " precursors at 1% FDR\n";
    }
    else
    {
      msg << "  no FDR reported: see the warning below\n";
    }
    msg << "  semi-supervised iterations: " << scored.iterations_trained
        << " trained, " << scored.iterations_skipped << " skipped";
    if (scored.precursors_without_candidate)
    {
      msg << "\n  " << scored.precursors_without_candidate
          << " precursors yielded no candidate peak group";
    }
    writeLogInfo_(msg.str());

    // Said loudly because it is the failure that looks like success: with no
    // iteration fitted, the d-scores are a single-feature initialisation and
    // the q-values are calibrated against it rather than against a model.
    if (!scored.fdr_valid && !scored.groups.empty())
    {
      if (scored.decoy_groups == 0)
      {
        writeLogWarn_("No decoy peak groups, so target/decoy FDR is undefined and "
                      "no q-values were computed. This is what -max_precursors "
                      "does on a target-only library: it slices the first N "
                      "precursors, and decoys are appended after all the targets. "
                      "Score the whole library, or one with decoys interleaved.");
      }
      else if (scored.target_groups == 0)
      {
        writeLogWarn_("No target peak groups; nothing to score against the decoys.");
      }
      else
      {
        writeLogWarn_("No semi-supervised iteration fitted a discriminant. The "
                      "scores come from a single-feature initialisation, not a "
                      "trained model -- do not read these q-values as an FDR.");
      }
    }

    if (!out.empty())
    {
      try
      {
        writeScores_(out, library, scored);
      }
      catch (const std::exception& e)
      {
        writeLogError_(std::string("Failed to write scores: ") + e.what());
        return CANNOT_WRITE_OUTPUT_FILE;
      }
      writeLogInfo_("wrote scored peak groups to " + out);
    }
    return EXECUTION_OK;
  }

  static void writeScores_(const std::string& path, const ODIA::Library& library,
                           const ODIA::PeakGroupScorer::Result& scored)
  {
    std::ofstream out(path);
    if (!out) { throw std::runtime_error("cannot open " + path); }
    out << "Precursor.Id\tDecoy\tRT\tLeft.RT\tRight.RT\tApex.Intensity"
           "\tDScore\tQValue\tPEP";
    for (const auto& n : ODIA::PeakGroupScorer::subScoreNames()) { out << '\t' << n; }
    out << '\n';

    const auto& p = library.precursors();
    for (const auto& g : scored.groups)
    {
      const auto seq = library.strings().get(p.modified_sequence[g.precursor]);
      out << seq << static_cast<int>(p.charge[g.precursor]) << '\t'
          << static_cast<int>(g.decoy) << '\t' << g.apex_rt << '\t' << g.left_rt
          << '\t' << g.right_rt << '\t' << g.apex_intensity << '\t' << g.dscore
          << '\t' << g.qvalue << '\t' << g.pep;
      for (const auto v : g.sub_scores) { out << '\t' << v; }
      out << '\n';
    }
    if (!out) { throw std::runtime_error("write failed for " + path); }
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
              << chromatograms.retentionTime(tr, j) << '\t'
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
    if (stop_after != "library" && stop_after != "extract" && stop_after != "score")
    {
      writeLogError_("Implemented stages are 'library', 'extract' and 'score'.");
      return ILLEGAL_PARAMETERS;
    }
    if ((stop_after == "extract" || stop_after == "score") && in_run.empty())
    {
      writeLogError_("-stop_after " + stop_after +
                     " needs a run to work on: give -in <file>.");
      return ILLEGAL_PARAMETERS;
    }

    ODIA::Library library;
    ODIA::Chromatograms chromatograms;

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
      const auto rc = runExtraction_(library, in_run, out_chrom, nullptr);
      if (rc != EXECUTION_OK) { return rc; }
    }
    if (stop_after == "score")
    {
      const auto rc = runScoreWorkflow_(library, in_run, out_chrom,
                                        getStringOption_("out"));
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
  /// The run's fitted mass model, measured once and reused by every pass.
  ///
  /// Cached because it is a property of the RUN, not of the pass: measuring it
  /// again in pass 2 would decode the same spectra to reach the same answer, and
  /// -- worse -- a second measurement taken through a window the first one
  /// narrowed would be a feedback loop that can only shrink.
  ODIA::MassCalibration::Model mass_model_;
  bool mass_model_known_ = false;

  /// Decide the fragment window's CENTRE and its WIDTH, in that order.
  ///
  /// They are two questions and they are answered from different things. The
  /// centre is a recalibration of the mass axis and comes from the run's own
  /// residuals; the width is the scatter left after that correction, bounded by
  /// what the caller was willing to accept.
  void applyMassCalibration_(const ODIA::Library& library, ODIA::SpectrumSource& source,
                             ODIA::ChromatogramExtractor::Options& options)
  {
    const double configured = getDoubleOption_("fragment_ppm");
    const double pinned = getDoubleOption_("fragment_ppm_offset");
    const bool off = getStringOption_("mass_calibration") == "off";
    const double search = getDoubleOption_("mz_calib_search_ppm");

    if (pinned != 0.0 || off || !(search > 0.0))
    {
      options.fragment_ppm_offset = pinned;
      // Uncalibrated means the WIDE width, and that is not a hedge: +/-10 about
      // zero keeps 0.51 of true fragments on this instrument where +/-15 keeps
      // 0.78, because the window is centred on the wrong place. A pinned offset
      // is a centring the caller asserted, so it earns the narrow width.
      const double fallback = pinned != 0.0 ? options.fragment_ppm
                                            : options.fragment_ppm_uncalibrated;
      options.fragment_ppm = configured > 0.0 ? configured : fallback;
      std::ostringstream os;
      os << "fragment mass calibration: not measured ("
         << (pinned != 0.0 ? "offset pinned by -fragment_ppm_offset"
                           : (off ? "-mass_calibration off" : "-mz_calib_search_ppm 0"))
         << "); extracting at " << options.fragment_ppm << " ppm about "
         << options.fragment_ppm_offset << " ppm";
      writeLogInfo_(os.str());
      return;
    }

    ODIA::MassCalibration::Diagnostics diagnostics;
    if (!mass_model_known_)
    {
      ODIA::MassCalibration::Options mzc;
      mzc.search_ppm = search;
      mzc.max_precursors = static_cast<std::size_t>(
        std::max(1, getIntOption_("mz_calib_precursors")));
      mzc.cycles = static_cast<std::size_t>(std::max(1, getIntOption_("mz_calib_cycles")));
      mzc.use_ion_mobility = !getFlag_("no_ion_mobility");
      try
      {
        mass_model_ = ODIA::MassCalibration::calibrate(library, source, mzc, &diagnostics);
      }
      catch (const std::exception& e)
      {
        writeLogWarn_(std::string("Mass calibration failed (") + e.what() +
                      "); extracting uncalibrated.");
        mass_model_ = ODIA::MassCalibration::Model{};
      }
      mass_model_known_ = true;
    }
    writeLogInfo_(ODIA::MassCalibration::report(mass_model_, &diagnostics));

    if (!mass_model_.fitted)
    {
      options.fragment_ppm_offset = 0.0;
      options.fragment_ppm = configured > 0.0 ? configured : options.fragment_ppm_uncalibrated;
      writeLogWarn_("The mass calibration gate FAILED, so no offset is applied and the window "
                    "stays wide at " + std::to_string(options.fragment_ppm) + " ppm. An "
                    "uncentred narrow window is the worse of the two errors: it keeps the tail "
                    "of the true distribution rather than its peak.");
      return;
    }

    options.fragment_ppm_offset = mass_model_.intercept_ppm;
    options.fragment_ppm_log_slope = mass_model_.log_slope_ppm;
    options.fragment_ppm_slope_per_1000 = mass_model_.linear_slope_ppm_per_1000;
    options.fragment_ppm_ref_mz = mass_model_.reference_mz;

    // With the window centred, the narrow width is the right one -- see the
    // header for the presence-versus-tolerance numbers that say so.
    const double baseline = configured > 0.0 ? configured : options.fragment_ppm;

    // CALIBRATION MAY ONLY NARROW. A fit that says "actually, use a wider
    // window" is telling you the fit failed, not that the instrument is bad, and
    // acting on it is strictly worse than doing nothing because it admits
    // interference the caller excluded. In the reference this rail was added
    // after an ungated estimate widened a window to 72.6 ppm on an instrument
    // measured at 1.66 ppm and took identifications from 6,798 to 4,496.
    //
    // Divergence from the reference, deliberate: there, an UNSET window let the
    // estimate stand however wide. Here it does not, because on this instrument
    // wide is measurably the wrong direction -- target-minus-decoy fragment
    // presence falls monotonically with tolerance (0.083 at 10 ppm down to 0.019
    // at 50), so 3 sigma of a scatter that is dominated by interference rather
    // than by measurement error would size the window from the interference.
    if (mass_model_.window_ppm > 0.0 && mass_model_.window_ppm < baseline)
    {
      options.fragment_ppm = mass_model_.window_ppm;
    }
    else
    {
      options.fragment_ppm = baseline;
      if (mass_model_.window_ppm > 0.0)
      {
        std::ostringstream os;
        os << "The mass calibration's own width, " << mass_model_.window_ppm
           << " ppm, is WIDER than the " << baseline << " ppm in force -- rejecting it. "
           << "Calibration may only narrow.";
        writeLogInfo_(os.str());
      }
    }

    std::ostringstream os;
    os << "extracting at +/-" << options.fragment_ppm << " ppm centred on "
       << options.fragment_ppm_offset << " ppm";
    if (options.fragment_ppm_log_slope != 0.0)
    {
      os << " + " << options.fragment_ppm_log_slope << " ppm per e-fold in m/z about "
         << options.fragment_ppm_ref_mz << " Th";
    }
    if (options.fragment_ppm_slope_per_1000 != 0.0)
    {
      os << " + " << options.fragment_ppm_slope_per_1000 << " ppm per 1000 Th about "
         << options.fragment_ppm_ref_mz << " Th";
    }
    writeLogInfo_(os.str());
  }

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
