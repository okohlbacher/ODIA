// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <odia/DIANNLibraryFile.h>
#include <odia/LibraryGenerator.h>
#include <odia/Library.h>
#include <fstream>
#include <odia/SpectrumSource.h>
#include <odia/ChromatogramExtractor.h>
#include <odia/ChromatogramTsv.h>
#include <odia/MassCalibration.h>
#include <odia/MobilityCalibration.h>
#include <odia/PeakGroupScorer.h>
#include <odia/RtCalibration.h>

#include <fstream>
#include <unordered_map>

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
                          "End the run after this stage and write its output. The "
                          "default runs to the end: 'library' when no -in is given, "
                          "'score' when one is.", false);
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
    registerStringOption_("ion_mobility_calibration", "<mode>", "auto",
                          "How the run's own 1/K0 prediction error is measured before the "
                          "mobility window is recentred on it. anchors: from the peak groups "
                          "this run has already SCORED, which needs a second pass or "
                          "-im_calib_anchors. prepass: from a blind probe of the run before "
                          "any pass, which is cheaper and measurably finds peak DENSITY rather "
                          "than precursors. auto: anchors when a scored pass will supply them, "
                          "prepass otherwise. off: extract on the library's 1/K0 as supplied. A "
                          "run with no ion mobility is untouched in every mode and says so.",
                          false);
    setValidStrings_("ion_mobility_calibration", {"auto", "off", "anchors", "prepass"});
    registerInputFile_("im_calib_anchors", "<file>", "",
                       "Take the 1/K0 anchors from this TSV instead of from a pass of this "
                       "run: columns Precursor.Id, Decoy, Apex.RT. That is what an EXTERNAL "
                       "scorer's confident identifications look like, and it is how the "
                       "stage is measured against a frozen discriminant without also moving "
                       "the retention-time axis.", false, true);
    setValidFormats_("im_calib_anchors", {"tsv"}, false);
    registerDoubleOption_("im_anchor_q", "<q>", 0.01,
                          "q-value below which a pass-1 identification becomes a 1/K0 anchor. "
                          "Tighter than -anchor_q on purpose: a retention-time map is fitted "
                          "from hundreds of anchors and a wrong one is an outlier the fit "
                          "steps over, whereas a wrong 1/K0 anchor contributes the mobility of "
                          "whatever the frame is dense at, which is a SYSTEMATIC and does not "
                          "average away.", false, true);
    registerIntOption_("im_calib_precursors", "<n>", 0,
                       "Precursors probed when measuring the 1/K0 error. 0 is the whole "
                       "library, which is the default because one precursor yields one "
                       "residual on this axis and the cost is the spectra decoded, not the "
                       "queries.", false, true);
    registerDoubleOption_("im_calib_rt_window", "<seconds>", 150.0,
                          "Look for a precursor's 1/K0 only within this much of where "
                          "-irt_slope/-irt_intercept say it elutes. 0, or no iRT map, searches "
                          "the whole gradient -- which is measurably worse, because the probe "
                          "keeps the brightest cluster over every block it looks in and an "
                          "absent precursor gets one draw from the interference per block.",
                          false, true);
    registerIntOption_("im_calib_cycles", "<n>", 200,
                       "Acquisition cycles probed for the 1/K0 measurement, drawn as short "
                       "CONTIGUOUS blocks so that a precursor has to be at the same mobility "
                       "in consecutive cycles to count.", false, true);
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
    registerIntOption_("max_live_precursors", "<n>", 0,
                       "Cap how many precursors may have chromatograms in memory at "
                       "once. 0 lets the retention-time overlap decide, which is the "
                       "cheap bound; above the cap the library is split into chunks "
                       "and each is a separate pass over the run, which costs a "
                       "decode. The run reports which of the two bound it.",
                       false, true);
    registerIntOption_("pass1_precursors", "<n>", 0,
                       "Sample about this many precursors for pass 1. It runs only to harvest "
                       "retention-time anchors -- 692 came from 2,450 precursors and "
                       "-min_anchors defaults to 20 -- so extracting the whole "
                       "library costs memory proportional to it: ~274 GiB at 4.26 M "
                       "precursors over the whole run, which is why the calibrated "
                       "pass-2 window does not unblock a large library on its own. "
                       "Indices are not renumbered, so anchors still refer to the "
                       "full library. This is a TARGET COUNT rather than a stride "
                       "because a stride is not scale-invariant: 4 starved a 2,450-"
                       "precursor library (anchors 709 -> 149, identifications 1862 "
                       "-> 635) while at 4.26 M it would still leave a million. 0 "
                       "disables sampling.",
                       false, true);
    registerIntOption_("refine_rounds", "<n>", 0,
                       "After pass 2, alternately refit the retention-time map from "
                       "the current best identifications and refit the discriminant, "
                       "until the identification count stops improving (three "
                       "consecutive rounds, as DIA-NN does) or this many rounds. "
                       "Cheap because the candidate picker is RT-agnostic and "
                       "RT_DELTA is the only map-dependent sub-score, so a round "
                       "recomputes one column rather than re-extracting the run. "
                       "DEFAULT 0 -- MEASURED AND IT DOES NOT HELP: on S08 it took "
                       "the map's p95 from 104.6 to 87.4 s and left best-ranked-right "
                       "unchanged to the decimal (75.0%, 77.4%) while costing 1.8-3.0 "
                       "points of precision. The map's real consumer is the pass-2 "
                       "extraction window, which has already run by then; RT_DELTA is "
                       "one feature of fifteen and cannot carry the improvement. Kept "
                       "for experiments.",
                       false, true);
    registerIntOption_("pass1_offset", "<n>", 0,
                       "Which residue class -pass1_precursors keeps. Diagnostic: it "
                       "lets equal-sized pass-1 subsets with different members be "
                       "compared, separating how many anchors from which anchors.",
                       false, true);
    registerDoubleOption_("rt_window_p95_factor", "<x>", 2.0,
                          "Pass 2 extracts over this many times the fitted map's p95 "
                          "residual, instead of the flat -rt_window. Points per "
                          "precursor scale directly with the width, so this is the "
                          "dominant memory lever in phase 2. On p95 rather than an SD "
                          "because a window has to cover the tail it is meant to "
                          "catch, and bounded below by -rt_window_min and above by "
                          "-rt_window, so it can only narrow. 0 disables it and "
                          "restores the flat window.",
                          false, true);
    registerDoubleOption_("rt_window_min", "<s>", 20.0,
                          "Floor for the residual-driven pass-2 window, seconds. "
                          "Guards the case where few anchors happen to agree, giving "
                          "a p95 too small for the calibration to actually support.",
                          false, true);
    registerFlag_("amplitude_picking",
                  "Detect candidate peaks as local maxima of the summed trace, the "
                  "way this tool used to. The default detects them by pairwise "
                  "correlation among the precursor's own fragments (DIA-NN's "
                  "Searcher::peaks), which on S08 raised availability of the true "
                  "peak from 69.1% to 97.4%, rank-1 accuracy from 40.7% to 75.7% "
                  "and precision from 21.4% to 89.1%. This flag restores the old "
                  "behaviour for comparison.", true);
    registerDoubleOption_("min_corr_score", "<x>", 0.5,
                          "With -coelution_picking: the reference fragment's summed "
                          "correlation to the others required for a position to be a "
                          "peak. DIA-NN's MinCorrScore.", false, true);
    registerDoubleOption_("max_corr_diff", "<x>", 2.0,
                          "With -coelution_picking: keep candidates within this much "
                          "of the best correlation sum. A margin, not a rank, so an "
                          "unambiguous precursor yields one candidate. DIA-NN's "
                          "MaxCorrDiff.", false, true);
    registerDoubleOption_("min_library_corr", "<x>", -1.0,
                          "Reject peak groups whose observed spectrum correlates "
                          "with the library below this. -1 disables it. Measured on "
                          "S08 against DIA-NN: median 0.582 for calls within 30 s of "
                          "the true apex, -0.036 for those that miss, -0.032 for "
                          "decoys -- so a misplaced target is indistinguishable from "
                          "a decoy here, which is why target-decoy FDR does not catch "
                          "it. A 0.5 cut keeps 54.1% of on-RT targets, 12.5% of "
                          "off-RT ones and 12.0% of decoys; the last two matching is "
                          "the label symmetry that keeps the FDR valid.",
                          false, true);
    registerIntOption_("decode_block", "<n>", 256,
                       "How many spectra are decoded and held at once. This is "
                       "the largest single term in the run's memory: a heap "
                       "profile put 5.57 GiB of a 9.45 GiB live peak in the "
                       "block's peak arrays. Cost is linear in this number and "
                       "results do not change -- only how much is resident. 0 "
                       "means the default.",
                       false, true);
    registerOutputFile_("out", "<file>", "",
                        "Write scored peak groups here (TSV).", false);
    setValidFormats_("out", {"tsv"}, false);
    registerStringOption_("classifier", "<name>", "gbt",
                          "Discriminant for the semi-supervised scorer. 'gbt' and "
                          "'xgboost' are the same learner -- second-order Newton leaf "
                          "values and XGBoost's split gain -- differing only in "
                          "configuration: 'xgboost' uses pyProphet 3.0.15's settings "
                          "for XGBoost 3.2.0 (max_depth 6, eta 0.3), 'gbt' the more "
                          "conservative depth 4 / eta 0.1 chosen when there were ten "
                          "sub-scores rather than fifteen. 'lda' is the linear "
                          "fallback and the only one the lower-is-better weight "
                          "constraint applies to.", false, true);
    setValidStrings_("classifier", {"xgboost", "lda", "gbt", "nn"});
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
                  "Ignore the ISOLATION WINDOW's mobility band when matching. Leaves "
                  "-precursor_im_window alone: the two are separate filters and this "
                  "flag is the control arm for measuring what the band is worth.", true);

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
  ExitCodes extractInto_(const ODIA::Library& library, const std::string& run,
                         ODIA::ChromatogramSink& sink,
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
    options.precursor_stride = pass_stride_;
    options.precursor_offset = pass_offset_;
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
    // After the mass calibration, because the mobility probe matches fragments
    // through the mass window and a mis-centred one fills its sample with the
    // interference that calibration exists to exclude; and after the iRT map,
    // because the probe uses it to look only where a precursor should elute.
    applyMobilityCalibration_(library, *source, options);
    options.threads = static_cast<unsigned>(std::max(1, getIntOption_("threads")));
    options.max_live_precursors = static_cast<std::size_t>(
      std::max(0, getIntOption_("max_live_precursors")));
    options.decode_block = static_cast<std::size_t>(
      std::max(0, getIntOption_("decode_block")));

    if (options.irt_slope == 0.0 && !library_rt_is_run_seconds)
    {
      writeLogWarn_("No iRT calibration given (-irt_slope/-irt_intercept). The "
                    "library is being spread evenly over the run, which will "
                    "extract from approximately the wrong retention times. "
                    "Treat the output as a smoke test, not a result.");
    }

    ODIA::ChromatogramExtractor::Stats stats;
    const auto t = std::chrono::steady_clock::now();
    try
    {
      ODIA::ChromatogramExtractor::extract(library, *source, options, sink, &stats);
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
        << stats.index_seconds << " s, match " << stats.match_seconds
        << " s, assemble " << stats.assemble_seconds << " s, sink "
        << stats.sink_seconds << " s\n"
        << "  points " << stats.points << " (" << stats.nonzero_points
        << " nonzero), all at once would be "
        << double(stats.points) * sizeof(float) / 1073741824.0 << " GiB\n"
        // The measurement the sliding window exists to move. `points` is what
        // the old flat array allocated; this is what was ever resident.
        << "  peak live " << stats.peak_live_precursors << " precursors, "
        << stats.peak_live_points << " points ("
        << double(stats.peak_live_points) * sizeof(float) / 1073741824.0
        << " GiB), bound by " << stats.memory_bound_by;
    if (stats.chunks > 1)
    {
      msg << "\n  " << stats.chunks << " chunks, " << stats.spectra_decoded
          << " spectra decoded against " << stats.spectra_read << " in the run";
    }
    if (stats.outside_rt_range)
    {
      msg << "\n  " << stats.outside_rt_range
          << " precursors predicted to elute outside the run";
    }
    if (stats.precursors_without_window)
    {
      msg << "\n  " << stats.precursors_without_window
          << " precursors covered by no isolation window";
    }
    writeLogInfo_(msg.str());
    return EXECUTION_OK;
  }

  /// Extract into one flat `Chromatograms`, write it if asked, and hand it on.
  ///
  /// This is the memory-bounded path: it keeps every point, so it is for a
  /// precursor count that fits. `-out_chrom` and every diagnostic built on it
  /// need it; scoring does not, and takes `extractInto_` with a scoring sink.
  ExitCodes runExtraction_(const ODIA::Library& library, const std::string& run,
                           const std::string& out_chrom,
                           ODIA::Chromatograms* keep = nullptr,
                           double rt_window_override = 0.0,
                           bool library_rt_is_run_seconds = false)
  {
    ODIA::ChromatogramCollector collector;
    const auto rc = extractInto_(library, run, collector, rt_window_override,
                                 library_rt_is_run_seconds);
    if (rc != EXECUTION_OK) { return rc; }
    ODIA::Chromatograms chromatograms = collector.take();
    writeLogInfo_("held all of them: " +
                  std::to_string(chromatograms.footprintBytes() / 1048576) + " MiB");

    if (!out_chrom.empty())
    {
      // Timed and reported. It was neither, despite being 33% of Phase-2 wall
      // -- ~301 s against 613 s of extraction at 9,522 precursors -- which had
      // to be recovered by subtracting the extractor's own timers from the
      // total.
      long long chrom_write_ms = 0;
      try
      {
        const auto t_write = std::chrono::steady_clock::now();
        ODIA::writeChromatogramTsv(out_chrom, library, chromatograms);
        chrom_write_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                           std::chrono::steady_clock::now() - t_write).count();
      }
      catch (const std::exception& e)
      {
        writeLogError_(std::string("Failed to write chromatograms: ") + e.what());
        return CANNOT_WRITE_OUTPUT_FILE;
      }
      writeLogInfo_("wrote chromatograms to " + out_chrom + " in " +
                    std::to_string(chrom_write_ms) + " ms");
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
    // Pass 2 has scored peak groups to measure the 1/K0 axis at, so `auto`
    // waits for them instead of guessing from a blind probe in pass 1.
    mobility_anchors_expected_ = passes > 1;
    ODIA::Chromatograms chromatograms;

    if (passes == 1)
    {
      // Nobody asked for the chromatograms, so nobody has to hold them. This is
      // the difference between a memory bill proportional to the library and
      // one proportional to what elutes at once.
      if (out_chrom.empty())
      {
        ODIA::PeakGroupScorer::Result scored;
        const auto rc = extractAndScore_(library, run, 0.0, false, scored);
        if (rc != EXECUTION_OK) { return rc; }
        return writeScoreResult_(scored, out, library);
      }
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
    ODIA::PeakGroupScorer::Result pass1;
    {
      // 0 means "the whole run", expressed as a window wider than any gradient
      // rather than as a sentinel the extractor would have to know about. A
      // negative value would simply extract nothing.
      //
      // Pass 1 exists to produce anchors, and nothing else ever reads its
      // chromatograms -- so it scores them as they finish and keeps none. That
      // matters most here: this is the WIDE pass, where every precursor is live
      // over most of the gradient and holding them all is at its worst.
      // A TARGET COUNT, not a stride. A stride is the wrong control because it
      // is not scale-invariant: 4 starves a 2,450-precursor library (measured:
      // anchors 709 -> 149, identifications 1862 -> 635) while at 4.26 M it
      // would still leave a million precursors, far more than any fit needs.
      pass_offset_ = static_cast<std::size_t>(std::max(0, getIntOption_("pass1_offset")));
      const int target = std::max(0, getIntOption_("pass1_precursors"));
      const std::size_t n_prec = library.precursorCount();
      pass_stride_ = (target > 0 && n_prec > static_cast<std::size_t>(target))
                       ? n_prec / static_cast<std::size_t>(target)
                       : 1;
      if (pass_stride_ > 1)
      {
        std::ostringstream st;
        st << "pass 1 extracts every " << pass_stride_ << "th precursor ("
           << (library.precursorCount() / pass_stride_) << " of "
           << library.precursorCount() << "); it exists to harvest anchors, and "
           << "extracting all of them costs memory proportional to the library";
        writeLogInfo_(st.str());
      }
      const auto rc = extractAndScore_(library, run,
                                       pass1_window > 0.0 ? pass1_window : 1.0e9,
                                       false, pass1);
      if (rc != EXECUTION_OK) { return rc; }
    }

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
        // By dscore, NOT by qvalue.
        //
        // The q-value is a per-PRECURSOR quantity: it is computed on the best
        // row of each group and then broadcast to every candidate of that
        // precursor (scoring/lda.h, assignQValues over one row per group).
        // So `g.qvalue < b->qvalue` is never true between two candidates of the
        // same precursor, and this kept whichever the stable sort put first --
        // library order, not the best peak. The anchor for the retention-time
        // fit was therefore an arbitrary candidate.
        //
        // dscore is genuinely per-candidate, and harvestMobilityAnchors_ has
        // always used it. The two anchor harvests now agree.
        if (b == nullptr || g.dscore > b->dscore) { b = &g; }
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

    // The 1/K0 anchors, harvested from the same pass and the same peak groups
    // as the retention-time ones, and kept for pass 2 to measure at.
    //
    // Two arms. The TARGETS are the confidently identified precursors, at the
    // apex of the group that identified them. The NULL is the best-scoring
    // DECOYS, as many of them as there are targets -- the same procedure
    // applied to precursors that are not in the sample, staking their claim on
    // whatever the interference offered. That is what the gate has to be able
    // to tell the targets apart from, and it is a harder null than an
    // m/z-shifted control: a decoy's fragments are real fragment masses of a
    // real (shuffled) sequence, and its apex sits on real signal.
    harvestMobilityAnchors_(pass1, library.precursorCount());

    // A pass that identified nothing at 1% has no business supplying anchors.
    //
    // Anchors are harvested at -anchor_q (0.05 by default), so a pass can fail
    // completely at 1% and still hand over a few hundred marginal anchors. That
    // is exactly what produced the zero-identification run above: 176 anchors
    // from q<=0.05, none at q<=0.01, a map fitted from them, and pass 2 failed
    // too. Better to keep pass 1's own scores than to calibrate on a pass that
    // did not work.
    if (pass1.identified_at_1pct == 0)
    {
      writeLogWarn_("pass 1 identified nothing at 1% FDR, so its q<=" +
                    std::to_string(getDoubleOption_("anchor_q")) +
                    " anchors are not evidence of anything. Not fitting a "
                    "retention-time map from them; returning pass 1's scores.");
      return writeScoreResult_(pass1, out, library);
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
      // Pass 1's own scores, not a re-score of chromatograms that no longer
      // exist. Same options, same traces, so the same answer.
      return writeScoreResult_(pass1, out, library);
    }

    double p95 = 0.0;
    const auto trafo = ODIA::Calibration::fit(anchors, &p95, 0.0);

    // p95 above is IN-SAMPLE: it is the residual on the very anchors the map
    // was fitted to, so it cannot see overfitting. That is not hypothetical --
    // sampling pass 1 down to 149 anchors left p95 at 81.7 s against 83.4 s for
    // 709 anchors, and the window width (which derives from p95) barely moved,
    // while identifications fell 66%. The fit was locally wrong wherever
    // anchors were thin, so window PLACEMENT degraded and nothing being
    // measured looked at placement.
    //
    // So refit on 80% and measure on the held-out 20%. A gap between the two is
    // the signal that the anchors are too few for the map they are being asked
    // to support.
    if (anchors.size() >= 25)
    {
      std::vector<std::pair<double, double>> fit_set, held;
      for (std::size_t i = 0; i < anchors.size(); ++i)
      {
        (i % 5 == 4 ? held : fit_set).push_back(anchors[i]);
      }
      double dummy = 0.0;
      const auto probe = ODIA::Calibration::fit(fit_set, &dummy, 0.0);
      std::vector<double> resid;
      resid.reserve(held.size());
      for (const auto& a : held)
      {
        resid.push_back(std::fabs(probe.apply(a.first) - a.second));
      }
      std::sort(resid.begin(), resid.end());
      const double oos = resid.empty() ? 0.0 : resid[std::size_t(0.95 * (resid.size() - 1))];
      std::ostringstream v;
      v << "map generalisation: p95 " << dummy << " s in-sample on "
        << fit_set.size() << " anchors, " << oos << " s out-of-sample on "
        << held.size() << " held out";
      if (dummy > 0.0 && oos > 2.0 * dummy)
      {
        v << " -- OUT-OF-SAMPLE IS " << (oos / dummy)
          << "x WORSE, the map is fitted to too few anchors and its windows will "
             "be misplaced where they are sparse";
        writeLogWarn_(v.str());
      }
      else { writeLogInfo_(v.str()); }
    }
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

    // Pass 2's window comes from the fit's OWN residual, not from a constant.
    //
    // Until now the p95 was computed, logged, and thrown away: pass 2 re-used
    // the flat -rt_window, so a run that had just measured its calibration to
    // be good extracted exactly as wide as one that had not. Points per
    // precursor scale directly with this width, and it is the dominant memory
    // lever in phase 2 -- measured at n100k, w600 -> w60 takes peak RSS from
    // 14.00 to 10.75 GiB, against 0.19 GiB for the extract/score fusion.
    //
    // The factor is on p95 rather than on an SD because a window has to cover
    // the tail it is meant to catch. The floor exists for the opposite failure:
    // few anchors that happen to agree give a tiny p95 and would place a window
    // narrower than the calibration can actually support. -rt_window stops
    // being the value and becomes the CAP, so this can only ever narrow.
    //
    // The measurement this guards against is on record: a flat 60 s was below
    // the run's own 76.6 s SD residual and put the true peak outside the window
    // for a third of precursors. Deriving the width from p95 is what makes a
    // narrow window safe rather than a gamble.
    double pass2_window = 0.0;
    const double p95_factor = getDoubleOption_("rt_window_p95_factor");
    if (p95_factor > 0.0 && std::isfinite(p95) && p95 > 0.0)
    {
      const double cap = getDoubleOption_("rt_window");
      const double floor_s = getDoubleOption_("rt_window_min");
      pass2_window = std::min(cap, std::max(floor_s, p95_factor * p95));
      std::ostringstream w;
      w << "pass 2 extraction window " << pass2_window << " s (" << p95_factor
        << " x p95 " << p95 << " s, floor " << floor_s << " s, cap " << cap << " s)";
      writeLogInfo_(w.str());
    }
    else
    {
      writeLogInfo_("pass 2 extraction window: flat -rt_window "
                    "(-rt_window_p95_factor 0 disables residual-driven narrowing)");
    }

    pass_stride_ = 1;  // pass 2 is the real search and needs every precursor
    pass_offset_ = 0;
    writeLogInfo_("pass 2 of 2: narrow extraction on the calibrated axis");
    chromatograms = ODIA::Chromatograms{};
    // The library now carries run seconds, so the affine map is the identity.
    if (out_chrom.empty())
    {
      ODIA::PeakGroupScorer::Result scored;
      scoring_rt_is_run_seconds_ = true;
      const auto rc = extractAndScore_(library, run, pass2_window, true, scored);
      if (rc != EXECUTION_OK) { return rc; }
      refineToConvergence_(library, original_irt, scored);
      return writeScoreResult_(scored, out, library);
    }
    const auto rc = runExtraction_(library, run, out_chrom, &chromatograms, pass2_window, true);
    if (rc != EXECUTION_OK) { return rc; }
    return runScoring_(library, chromatograms, out);
  }

  ODIA::PeakGroupScorer::Options scoringOptions_()
  {
    ODIA::PeakGroupScorer::Options options;
    options.classifier = getStringOption_("classifier");
    options.min_library_corr = getDoubleOption_("min_library_corr");
    options.coelution_picking = !getFlag_("amplitude_picking");
    options.min_corr_score = getDoubleOption_("min_corr_score");
    options.max_corr_diff = getDoubleOption_("max_corr_diff");
    options.max_candidates = static_cast<std::size_t>(
      std::max(1, getIntOption_("max_candidates")));
    options.threads = static_cast<unsigned>(std::max(1, getIntOption_("threads")));
    return options;
  }

  /// Extract and score in one forward pass, holding only what is live.
  ///
  /// The chromatograms are never all in memory at once: each precursor is
  /// scored as the pass leaves its retention-time window and then freed. What
  /// survives is the peak-group table, which is ~120 bytes a group against tens
  /// of kilobytes a chromatogram.
  /// Precursor stride for the pass being run. Pass 1 sets it from
  /// -pass1_stride and pass 2 restores 1, because pass 1 extracts only to
  /// harvest RT anchors and needs a few hundred of them, not a library's worth.
  std::size_t pass_stride_ = 1;
  /// Set before pass 2, when the library's irt has been rewritten to run
  /// seconds. RT_DELTA is only meaningful then.
  bool scoring_rt_is_run_seconds_ = false;
  std::size_t pass_offset_ = 0;

  /// Refit the retention-time map and the discriminant, alternately, until the
  /// identification count stops moving.
  ///
  /// This is the iteration schedule DIA-NN runs twelve rounds of and we ran
  /// two. It is affordable for one reason: the candidate picker is
  /// retention-time agnostic, and RT_DELTA is the only sub-score that depends
  /// on the map. So a round costs one recomputed column and one classifier fit
  /// over groups already in memory -- no re-extraction, no second decode of a
  /// file that is ~98% of the run's time.
  ///
  /// Convergence rather than a fixed count: DIA-NN stops after three
  /// consecutive rounds without improvement (diann.cpp:10448-10476). The same
  /// rule here, with a hard cap so a pathological run cannot spin.
  ///
  /// The map is refit from the CURRENT best group per precursor, so each round
  /// draws anchors from a better-scored set than the last. That is the whole
  /// mechanism: better anchors -> better map -> better RT_DELTA -> better
  /// discriminant -> better anchors.
  void refineToConvergence_(ODIA::Library& library,
                            const std::vector<float>& original_irt,
                            ODIA::PeakGroupScorer::Result& scored)
  {
    const int max_rounds = std::max(0, getIntOption_("refine_rounds"));
    if (max_rounds == 0) { return; }
    const double anchor_q = getDoubleOption_("anchor_q");
    const int min_anchors = std::max(1, getIntOption_("min_anchors"));

    auto options = scoringOptions_();
    options.library_rt_is_run_seconds = true;

    std::size_t best_ids = scored.identified_at_1pct;
    int stagnant = 0;
    for (int round = 1; round <= max_rounds && stagnant < 3; ++round)
    {
      // Anchors from the current best group per precursor, by dscore -- the
      // q-value is broadcast across a precursor's candidates and cannot
      // discriminate between them.
      std::vector<const ODIA::PeakGroupScorer::PeakGroup*> best(
        library.precursorCount(), nullptr);
      for (const auto& g : scored.groups)
      {
        if (g.decoy || g.qvalue > anchor_q) { continue; }
        auto*& b = best[g.precursor];
        if (b == nullptr || g.dscore > b->dscore) { b = &g; }
      }
      std::vector<std::pair<double, double>> anchors;
      for (std::size_t i = 0; i < best.size(); ++i)
      {
        if (best[i] != nullptr && std::isfinite(original_irt[i]))
        {
          anchors.push_back({static_cast<double>(original_irt[i]),
                             static_cast<double>(best[i]->apex_rt)});
        }
      }
      if (static_cast<int>(anchors.size()) < min_anchors) { break; }

      double p95 = 0.0;
      const auto trafo = ODIA::Calibration::fit(anchors, &p95, 0.0);
      // Applied to the ORIGINAL iRT every round, never to the previous round's
      // output: composing maps would drift, and each fit is a map from library
      // units to run seconds, not a correction to the last one.
      auto& irt = library.precursors().irt;
      for (std::size_t i = 0; i < irt.size(); ++i)
      {
        if (std::isfinite(original_irt[i]))
        {
          irt[i] = static_cast<float>(trafo.apply(static_cast<double>(original_irt[i])));
        }
      }

      ODIA::PeakGroupScorer::refit(library, scored, options);

      std::ostringstream m;
      m << "refine round " << round << ": " << anchors.size() << " anchors, p95 "
        << p95 << " s, " << scored.identified_at_1pct << " identified at 1% FDR";
      if (scored.identified_at_1pct > best_ids)
      {
        best_ids = scored.identified_at_1pct;
        stagnant = 0;
      }
      else { ++stagnant; m << " (no gain, " << stagnant << " of 3)"; }
      writeLogInfo_(m.str());
    }
  }

  ExitCodes extractAndScore_(const ODIA::Library& library, const std::string& run,
                             double rt_window_override, bool library_rt_is_run_seconds,
                             ODIA::PeakGroupScorer::Result& scored)
  {
    auto options = scoringOptions_();
    options.library_rt_is_run_seconds = scoring_rt_is_run_seconds_;
    ODIA::PeakGroupScorer::Sink sink(library, options);
    const auto t = std::chrono::steady_clock::now();
    const auto rc = extractInto_(library, run, sink, rt_window_override,
                                 library_rt_is_run_seconds);
    if (rc != EXECUTION_OK) { return rc; }
    try
    {
      scored = sink.finish();
    }
    catch (const std::exception& e)
    {
      writeLogError_(std::string("Scoring failed: ") + e.what());
      return INTERNAL_ERROR;
    }
    const auto ms = std::chrono::duration<double, std::milli>(
                      std::chrono::steady_clock::now() - t).count();
    reportScoring_(scored, options.classifier, ms, "scored on the fly");
    return EXECUTION_OK;
  }

  /// Find peak groups, score them, and write them with their q-values.
  ExitCodes runScoring_(const ODIA::Library& library,
                        const ODIA::Chromatograms& chromatograms,
                        const std::string& out)
  {
    const auto options = scoringOptions_();

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
    reportScoring_(scored, options.classifier, ms, "scored");
    return writeScoreResult_(scored, out, library);
  }

  /// Everything the scoring stage has to say, whichever path produced it.
  void reportScoring_(const ODIA::PeakGroupScorer::Result& scored,
                      const std::string& classifier, double ms,
                      const std::string& how)
  {
    std::ostringstream msg;
    msg << how << " " << scored.groups.size() << " peak groups with "
        << classifier << " in " << ms << " ms\n"
        << "  " << scored.target_groups << " target / " << scored.decoy_groups
        << " decoy groups\n";
    if (scored.fdr_valid)
    {
      msg << "  identified " << scored.identified_at_1pct << " precursors at 1% FDR\n";
      // Zero identifications from a run that produced peak groups and trained a
      // classifier is a failure, and until now it was reported as a result.
      //
      // Measured: four runs differing ONLY in which 1,225 of 4,900 precursors
      // pass 1 saw gave 635, 1820, 1211 and 0 identifications. The zero arm had
      // the second-most anchors, a middling p95 and the WIDEST window -- every
      // number this tool prints looked healthy -- while emitting 65,733 peak
      // groups and identifying none of them. The same signature appeared twice
      // more tonight (a library-correlation gate at 0.5, and an early narrowing
      // arm), so it is a recurring mode rather than one bad configuration.
      if (scored.identified_at_1pct == 0 && !scored.groups.empty())
      {
        std::ostringstream z;
        z << "identified NOTHING at 1% FDR from " << scored.groups.size()
          << " peak groups (" << scored.target_groups << " target, "
          << scored.decoy_groups << " decoy). The scorer ran and the classifier "
          << "trained; the target-decoy threshold then rejected everything. Do not "
          << "read the other numbers as healthy -- p95, window width and anchor "
          << "counts are all reported normally in this state.";
        writeLogWarn_(z.str());
      }
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
  }

  ExitCodes writeScoreResult_(const ODIA::PeakGroupScorer::Result& scored,
                              const std::string& out, const ODIA::Library& library)
  {
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

  ExitCodes main_(int, const char**) override
  {
    const std::string tr = getStringOption_("tr");
    const std::string fasta = getStringOption_("fasta");
    const std::string out_lib = getStringOption_("out_lib");
    std::string stop_after = getStringOption_("stop_after");
    const bool sort_library = getFlag_("sort_library");

    if (tr.empty() == fasta.empty())
    {
      writeLogError_("Give exactly one of -tr <library> or -fasta <proteins>.");
      return ILLEGAL_PARAMETERS;
    }

    const std::string in_run = getStringOption_("in");
    const std::string out_chrom = getStringOption_("out_chrom");

    // Naming no stage means "run to the end", and which end that is depends on
    // whether there is a run to work on.
    //
    // This option is registered with a default of "" and `setValidStrings_`
    // accepts "", but the check below rejected it -- so the documented default
    // was rejected by the tool's own validation and NO invocation without an
    // explicit -stop_after could start. It is resolved here rather than by
    // changing the registered default because the right end differs: with only
    // -tr the run has nothing to extract from and stops at the library, and
    // erroring with "needs -in" would be wrong for a caller who only wanted one.
    if (stop_after.empty()) { stop_after = in_run.empty() ? "library" : "score"; }

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
      // With no -out_chrom nothing wants the points, so nothing holds them:
      // this is then a decode-and-match benchmark that runs at any library
      // size. With -out_chrom it is the old path, and it is bounded by the
      // library exactly as it always was.
      if (out_chrom.empty())
      {
        ODIA::NullChromatogramSink sink;
        const auto rc = extractInto_(library, in_run, sink);
        if (rc != EXECUTION_OK) { return rc; }
      }
      else
      {
        const auto rc = runExtraction_(library, in_run, out_chrom, nullptr);
        if (rc != EXECUTION_OK) { return rc; }
      }
    }
    if (stop_after == "score")
    {
      const auto rc = runScoreWorkflow_(library, in_run, out_chrom,
                                        getStringOption_("out"));
      if (rc != EXECUTION_OK) { return rc; }
    }

    if (!out_lib.empty())
    {
      // Timed and reported. It was neither, and it is 26% of Phase 1 -- 211.7 s
      // of 821.1 s on the human library -- so the stage table had to obtain it
      // by subtracting the generator's own `load time` from the total wall.
      // That also makes it the only Phase-1 item no accelerator touches: on a
      // GPU, where inference is ~4 min, this write is the largest single item.
      const auto t_write = std::chrono::steady_clock::now();
      try
      {
        ODIA::DIANNLibraryFile::storeTSV(out_lib, library);
      }
      catch (const std::exception& e)
      {
        writeLogError_(std::string("Failed to write assay library: ") + e.what());
        return CANNOT_WRITE_OUTPUT_FILE;
      }
      const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::steady_clock::now() - t_write).count();
      writeLogInfo_("wrote assay library to " + out_lib + " in " +
                    std::to_string(ms) + " ms");
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

  /// The run's fitted 1/K0 model, cached for the same reason the mass model is:
  /// it is a property of the RUN, and a second measurement taken through a
  /// window the first pass narrowed is a feedback loop that can only shrink.
  ODIA::MobilityCalibration::Model mobility_model_;
  bool mobility_model_known_ = false;

  /// Where the 1/K0 measurement is allowed to look, and whether any such place
  /// is coming. Empty with `mobility_anchors_expected_` set means "a scored
  /// pass will fill this, do not measure yet"; empty without it means the
  /// blind probe is the only source there is.
  std::vector<ODIA::MobilityAnchor> mobility_anchors_;
  bool mobility_anchors_expected_ = false;
  bool mobility_anchors_loaded_ = false;

  /// Turn pass 1's peak groups into 1/K0 anchors: the confident targets, and a
  /// rank-matched null of the best-scoring decoys.
  void harvestMobilityAnchors_(const ODIA::PeakGroupScorer::Result& pass1,
                               std::size_t precursors)
  {
    mobility_anchors_.clear();
    if (getStringOption_("ion_mobility_calibration") == "off" || !pass1.fdr_valid) { return; }
    const double q = getDoubleOption_("im_anchor_q");

    std::vector<const ODIA::PeakGroupScorer::PeakGroup*> best(precursors, nullptr);
    std::vector<const ODIA::PeakGroupScorer::PeakGroup*> best_decoy(precursors, nullptr);
    for (const auto& g : pass1.groups)
    {
      if (g.precursor >= precursors) { continue; }
      if (g.decoy)
      {
        auto*& b = best_decoy[g.precursor];
        if (b == nullptr || g.dscore > b->dscore) { b = &g; }
        continue;
      }
      if (g.qvalue > q) { continue; }
      auto*& b = best[g.precursor];
      if (b == nullptr || g.dscore > b->dscore) { b = &g; }
    }

    for (std::size_t i = 0; i < precursors; ++i)
    {
      if (best[i] != nullptr && std::isfinite(best[i]->apex_rt))
      {
        mobility_anchors_.push_back({static_cast<std::uint32_t>(i), best[i]->apex_rt, false});
      }
    }
    const std::size_t targets = mobility_anchors_.size();

    // As many decoys as there are targets, best first. Equal size on purpose:
    // the gate compares two peakedness statistics, and a null with a tenth of
    // the sample would be compared on its noise.
    std::vector<std::pair<double, std::uint32_t>> decoys;
    for (std::size_t i = 0; i < precursors; ++i)
    {
      if (best_decoy[i] != nullptr && std::isfinite(best_decoy[i]->apex_rt))
      {
        decoys.emplace_back(best_decoy[i]->dscore, static_cast<std::uint32_t>(i));
      }
    }
    std::sort(decoys.begin(), decoys.end(), std::greater<>());
    if (decoys.size() > targets) { decoys.resize(targets); }
    for (const auto& d : decoys)
    {
      mobility_anchors_.push_back({d.second, best_decoy[d.second]->apex_rt, true});
    }

    std::ostringstream os;
    os << "pass 1 offers " << targets << " 1/K0 anchors at q <= " << q << ", against a null of "
       << decoys.size() << " best-scoring decoys";
    writeLogInfo_(os.str());
  }

  /// Read anchors from -im_calib_anchors, once.
  ///
  /// The file names precursors the way every other TSV here does -- modified
  /// sequence followed by charge -- and a target and its decoy share that name,
  /// so the Decoy column is not optional.
  void loadMobilityAnchors_(const ODIA::Library& library)
  {
    if (mobility_anchors_loaded_) { return; }
    mobility_anchors_loaded_ = true;
    const std::string path = getStringOption_("im_calib_anchors");
    if (path.empty()) { return; }
    mobility_anchors_expected_ = true;

    std::unordered_map<std::string, std::uint32_t> index;
    const auto& p = library.precursors();
    index.reserve(library.precursorCount() * 2);
    for (std::size_t i = 0; i < library.precursorCount(); ++i)
    {
      std::string key(library.strings().get(p.modified_sequence[i]));
      key += std::to_string(static_cast<int>(p.charge[i]));
      key += p.decoy[i] ? '-' : '+';
      index.emplace(std::move(key), static_cast<std::uint32_t>(i));
    }

    std::ifstream in(path);
    if (!in)
    {
      writeLogWarn_("Cannot read -im_calib_anchors " + path +
                    "; the 1/K0 axis is left as the library supplies it.");
      return;
    }
    std::string line;
    if (!std::getline(in, line))
    {
      writeLogWarn_("-im_calib_anchors " + path + " is empty.");
      return;
    }
    const auto split = [](const std::string& row) {
      std::vector<std::string> out;
      std::size_t b = 0;
      for (std::size_t i = 0; i <= row.size(); ++i)
      {
        if (i == row.size() || row[i] == '\t') { out.push_back(row.substr(b, i - b)); b = i + 1; }
      }
      return out;
    };
    const auto header = split(line);
    const auto column = [&](const char* name) {
      for (std::size_t i = 0; i < header.size(); ++i)
      {
        if (header[i] == name) { return static_cast<int>(i); }
      }
      return -1;
    };
    const int c_id = column("Precursor.Id"), c_d = column("Decoy"), c_rt = column("Apex.RT");
    if (c_id < 0 || c_d < 0 || c_rt < 0)
    {
      writeLogWarn_("-im_calib_anchors " + path + " needs columns Precursor.Id, Decoy and "
                    "Apex.RT; the 1/K0 axis is left as the library supplies it.");
      return;
    }
    std::size_t unmatched = 0, targets = 0;
    while (std::getline(in, line))
    {
      if (line.empty()) { continue; }
      const auto f = split(line);
      if (static_cast<int>(f.size()) <= std::max(c_id, std::max(c_d, c_rt))) { continue; }
      const bool decoy = f[c_d] != "0" && !f[c_d].empty();
      const auto it = index.find(f[c_id] + (decoy ? '-' : '+'));
      if (it == index.end()) { ++unmatched; continue; }
      mobility_anchors_.push_back({it->second,
                                   static_cast<float>(std::strtod(f[c_rt].c_str(), nullptr)),
                                   decoy});
      if (!decoy) { ++targets; }
    }
    std::ostringstream os;
    os << "read " << mobility_anchors_.size() << " 1/K0 anchors from " << path << " ("
       << targets << " target, " << mobility_anchors_.size() - targets << " null)";
    if (unmatched) { os << "; " << unmatched << " named no precursor in this library"; }
    writeLogInfo_(os.str());
  }

  /// Measure the run's 1/K0 prediction error and, if the gate passes, hand the
  /// model to the extractor.
  ///
  /// Structured exactly like applyMassCalibration_, with one state it does not
  /// have: a run with no ion mobility, or a library with no 1/K0, is a NO-OP
  /// that is reported as such. That is not a failed calibration and logging it
  /// as one would train the reader to ignore a warning that on another run
  /// means something.
  void applyMobilityCalibration_(const ODIA::Library& library, ODIA::SpectrumSource& source,
                                 ODIA::ChromatogramExtractor::Options& options)
  {
    options.mobility_model = nullptr;
    const std::string mode = getStringOption_("ion_mobility_calibration");
    if (mode == "off")
    {
      writeLogInfo_("ion-mobility calibration: not measured (-ion_mobility_calibration off); "
                    "extracting on the library's 1/K0 as supplied");
      return;
    }
    loadMobilityAnchors_(library);
    // auto resolves to the anchored probe exactly when a scored pass will
    // supply anchors. It is not the default because it is better in principle
    // -- it is the default because the blind probe was measured on S08 and
    // found peak density rather than precursors; see MobilityCalibration.h.
    const bool anchored = mode == "anchors" ||
                          (mode == "auto" && (mobility_anchors_expected_ ||
                                              !mobility_anchors_.empty()));
    if (anchored && mobility_anchors_.empty() && !mobility_model_known_)
    {
      writeLogInfo_("ion-mobility calibration: DEFERRED -- it is measured at the peak groups "
                    "this run scores, and none have been scored yet. This pass extracts on the "
                    "library's 1/K0; the next one is where the correction can be earned.");
      return;
    }
    if (options.precursor_im_window <= 0.0)
    {
      writeLogInfo_("ion-mobility calibration: not measured (-precursor_im_window 0, so the "
                    "per-precursor mobility window is switched off and there is nothing for a "
                    "recentring to move)");
      return;
    }

    ODIA::MobilityCalibration::Diagnostics diagnostics;
    if (!mobility_model_known_)
    {
      ODIA::MobilityCalibration::Options imc;
      imc.max_precursors = static_cast<std::size_t>(
        std::max(0, getIntOption_("im_calib_precursors")));
      imc.cycles = static_cast<std::size_t>(std::max(1, getIntOption_("im_calib_cycles")));
      // Probe through the mass window that is about to be extracted with, so
      // the two calibrations cannot disagree about what a fragment match is.
      imc.fragment_ppm = options.fragment_ppm;
      imc.fragment_ppm_offset = options.fragment_ppm_offset;
      imc.fragment_ppm_log_slope = options.fragment_ppm_log_slope;
      imc.fragment_ppm_slope_per_1000 = options.fragment_ppm_slope_per_1000;
      imc.fragment_ppm_ref_mz = options.fragment_ppm_ref_mz;
      // The run's own iRT map, so the probe looks only where a precursor should
      // be. Measured on S08: it takes the data-against-control peakedness margin
      // from 1.01x to 1.08x and the charge-2 centre from +0.0060 to +0.0016,
      // where an independent check against DIA-NN's observed 1/K0 says +0.0017.
      // Still short of the 1.25x the gate wants, which is the run's answer, not
      // a reason to leave the information unused.
      imc.irt_slope = options.irt_slope;
      imc.irt_intercept = options.irt_intercept;
      imc.rt_window_seconds = getDoubleOption_("im_calib_rt_window");
      try
      {
        mobility_model_ = anchored
          ? ODIA::MobilityCalibration::calibrateFrom(library, source, mobility_anchors_, imc,
                                                     &diagnostics)
          : ODIA::MobilityCalibration::calibrate(library, source, imc, &diagnostics);
      }
      catch (const std::exception& e)
      {
        writeLogWarn_(std::string("Ion-mobility calibration failed (") + e.what() +
                      "); extracting on the library's 1/K0.");
        mobility_model_ = ODIA::MobilityCalibration::Model{};
      }
      mobility_model_known_ = true;
    }
    writeLogInfo_(ODIA::MobilityCalibration::report(mobility_model_, &diagnostics));

    if (!mobility_model_.fitted)
    {
      // Deliberately not a warning when there is no axis: nothing is wrong.
      if (mobility_model_.run_has_mobility && mobility_model_.library_has_mobility)
      {
        writeLogWarn_("The ion-mobility calibration gate FAILED, so the library's 1/K0 is used "
                      "uncorrected. That is the safe direction: an uncentred window keeps the "
                      "library's own error, where a window recentred on a badly measured "
                      "offset moves off the precursor entirely.");
      }
      return;
    }
    options.mobility_model = &mobility_model_;
  }

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
