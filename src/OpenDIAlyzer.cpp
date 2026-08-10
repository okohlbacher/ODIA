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
#include <odia/MassWidth.h>
#include <odia/Ms1Traces.h>
#include <odia/MobilityCalibration.h>
#include <odia/PeakGroupScorer.h>
#include <odia/RtCalibration.h>
#include <odia/RtRefiner.h>

#include <array>
#include <limits>

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
                       "downloads when built with WITH_ONNX=ON, which is the right "
                       "choice for building a library: a model fine-tuned on one run "
                       "must not predict for another. Doing so cost 2,027 confident "
                       "precursors the last time (library v5 at 35,556 against v4's "
                       "37,583, v4 having used the stock model). Per-run fine-tuning "
                       "belongs inside that run's calibration loop, built from it and "
                       "discarded with it.", false);

    registerIntOption_("missed_cleavages", "<n>", 1, "Maximum missed cleavages.", false, true);
    registerIntOption_("min_peptide_length", "<n>", 7, "Minimum peptide length.", false, true);
    registerIntOption_("max_peptide_length", "<n>", 30, "Maximum peptide length.", false, true);

    registerOutputFile_("out_lib", "<file>", "",
                        "Write the assay library here (DIA-NN TSV).", false);
    setValidFormats_("out_lib", {"tsv"}, false);

    registerStringOption_("stop_after", "<stage>", "",
                          "End the run after this stage and write its output. The "
                          "default runs to the end: 'library' when no -in is given, "
                          "'score' when one is. 'calib' runs pass 1, fits the "
                          "retention-time map, reports its residuals and stops -- the "
                          "calibration is the thing being measured, so pass 2 would only "
                          "cost time.", false);
    setValidStrings_("stop_after", {"", "library", "extract", "calib", "score"});

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
    registerFlag_("mass_calibration_offset_only",
                  "Apply the fitted mass offset but do NOT narrow the window from it. Separates "
                  "the two things a passing gate does, so that a loss can be attributed to the "
                  "centring or to the width. Diagnostic; not a production setting.",
                  true);
    registerFlag_("mass_calibration_remeasure",
                  "Re-measure the fragment mass calibration in pass 2, against the fitted "
                  "retention-time map, instead of reusing the pass-1 model taken before any map "
                  "existed. The re-measurement is strictly better as a MEASUREMENT -- on Astral "
                  "it takes the control from 98 residuals to 149 and stops it out-peaking the "
                  "data -- but applying it costs 4,275 -> 1,895 identifications there, because a "
                  "passing gate also narrows fragment_ppm off the uncalibrated 15. Off until the "
                  "width and the offset are separable.",
                  true);
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
    registerDoubleOption_("max_im_slope", "<1/K0 per 1/K0>", 0.25,
                          "Bound on the mobility-linear term of the 1/K0 recalibration. The "
                          "CCS->1/K0 conversion is a pure proportionality, so a scale error in "
                          "its coefficient shows up as a slope: measured -0.096 (charge 2) and "
                          "-0.113 (charge 3) on S08, i.e. a ~10% coefficient error, which moved "
                          "the residual from +0.017 at 1/K0 0.7 to -0.030 at 1.3. Set 0 to "
                          "disable the term and fit a constant offset only.",
                          false, true);
    registerIntOption_("im_calib_pooled_slope_min", "<n>", 0,
                       "Anchors below which a charge gets the REDUCED 1/K0 model -- its own "
                       "constant plus the slope pooled over all charges, no m/z shape -- rather "
                       "than no correction at all. The slope is a scale error in the CCS->1/K0 "
                       "coefficient and so is shared by every charge; a constant offset is not, "
                       "and is still fitted per charge. Measured need: charge 3 brings 93 anchors "
                       "against a 120 minimum on S08 and is the charge with the WORSE residual. "
                       "DEFAULTED OFF: it corrects charge 3 and improves the out-of-fold 1/K0 "
                       "error 14.1%% -> 18.8%%, and still COSTS 20 identifications (1232 -> 1212) "
                       "-- the pooled slope is dominated by charge 2's 378 anchors against charge "
                       "3's 93, so it undercorrects. Set 40 to enable.",
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
    registerOutputFile_("out_anchors", "<file>", "",
                        "Write pass 1's retention-time anchors here as a TSV: precursor, "
                        "modified sequence, charge, library iRT, observed apex RT, dscore and "
                        "q-value. These are the (sequence, observed RT) pairs a retention-time "
                        "model would be fine-tuned on, and nothing could emit them before -- so "
                        "any experiment on per-run RT refinement had to reconstruct them from a "
                        "different search. Written before the map is applied, so the iRT column "
                        "is the LIBRARY value.", false, true);
    setValidFormats_("out_anchors", {"tsv"}, false);
    registerStringOption_("rt_interpolation", "<type>", "akima",
                          "How the retention-time map joins its knots. akima: nonlinear and "
                          "outlier-resistant. cspline: nonlinear but rings around outliers, and "
                          "pass 1's anchors reach 1,600 s of residual. linear: straight segments, "
                          "which is what this did until 2026-08-10 -- the derivative jumps at "
                          "every knot and a retention-time map has no reason to be "
                          "piecewise-linear. OpenSWATH aligns iRT with LOWESS and DIA-NN fits a "
                          "nonlinear monotone regression; neither joins knots with line segments.",
                          false, true);
    setValidStrings_("rt_interpolation", {"akima", "cspline", "linear"});
    registerDoubleOption_("rt_loess_span", "<fraction>", 0.0,
                          "LOESS span for the retention-time map, as a fraction of the anchors. "
                          "0 uses binned medians alone. "
                          "\n\nNOTE: 0 has ALWAYS been what production passed, so the LOESS the "
                          "header describes has never run -- RtCalibration.h says '0 lets the fit "
                          "choose' and RtCalibration.cpp says 'loess_span <= 0 selects the "
                          "historical binned-median path'. The header was wrong and this option "
                          "is how the difference gets measured instead of assumed. LOESS also "
                          "needs at least 50 anchors; below that it is skipped whatever is set.",
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
                          "Discriminant for the semi-supervised scorer. 'percolator' is "
                          "OpenMS 3.6's in-process Percolator -- a cross-validated linear SVM "
                          "of the same mProphet lineage, and the algorithm mokapot "
                          "reimplements in Python; there is no mokapot in OpenMS, and "
                          "in-process Percolator is strictly better for a C++ tool than "
                          "shelling out to one. It cannot represent feature interactions the "
                          "way 'gbt' can, which is the point of comparing them. 'gbt' and "
                          "'xgboost' are the same learner -- second-order Newton leaf "
                          "values and XGBoost's split gain -- differing only in "
                          "configuration: 'xgboost' uses pyProphet 3.0.15's settings "
                          "for XGBoost 3.2.0 (max_depth 6, eta 0.3), 'gbt' the more "
                          "conservative depth 4 / eta 0.1 chosen when there were ten "
                          "sub-scores rather than fifteen. 'lda' is the linear "
                          "fallback and the only one the lower-is-better weight "
                          "constraint applies to.", false, true);
    setValidStrings_("classifier", {"xgboost", "lda", "gbt", "nn", "percolator"});
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
    registerStringOption_("classifier_model_out", "<file>", "",
                          "Train the discriminant on THIS run and write it here. Only with "
                          "-classifier percolator. Use on a run where the semi-supervised loop "
                          "ignites -- a library that is mostly present.", false, true);
    registerStringOption_("classifier_model_in", "<file>", "",
                          "Apply a frozen discriminant from this file instead of training. "
                          "Only with -classifier percolator. This is static modelling, the "
                          "documented remedy for a run whose true-positive rate is too low to "
                          "bootstrap: on v6_50k (1.5% present) every engine certifies nothing "
                          "at 1% FDR while the same discriminant ranks 232 of DIA-NN's 738 into "
                          "its top 738. A discriminant over these sub-scores is a property of "
                          "the instrument and the scoring code, not of which peptides happen to "
                          "be in one sample.", false, true);
    registerDoubleOption_("apex_evidence", "<fraction>", 0.99,
                          "A scan position is a peak only if the reference fragment's smoothed "
                          "intensity there is at least this fraction of its maximum nearby. "
                          "DIA-NN's PeakApexEvidence, and the second-largest picker filter "
                          "after -min_corr_score: 36.9M rejections against 50.7M on a 107.6M "
                          "scan-position S08 run. At 0.99 it demands the position be within 1% "
                          "of the local maximum, which on a noisy trace is nearly exact "
                          "equality.",
                          false, true);
    registerIntOption_("min_fragments_at_apex", "<n>", 1,
                       "A candidate peak group is dropped unless at least this many of its "
                       "fragments have signal at the apex. "
                       "\n\nDEFAULT LOWERED 3 -> 1 on 2026-08-09, measured on both files: "
                       "Astral 4,969 -> 5,025 and S08 1,306 -> 1,342 (at 2 it is 1,340, so the "
                       "curve is flat and 3 was simply too strict). The reason it costs nothing "
                       "to relax is that `var_fragment_coverage` IS at_apex/tc -- the classifier "
                       "already has the quantity this gate thresholds on, so the gate was "
                       "DELETING ROWS over a feature that could have ranked them. 1 keeps it as "
                       "a validity check (a group needs some signal at its own apex) and gives "
                       "the discrimination back to the discriminant. "
                       "\n\nCoupled to -fragment_ppm in a way "
                       "that is easy to miss: narrowing the mass window leaves fewer "
                       "transitions carrying any signal, so precursors fall below this gate "
                       "and yield no candidate at all. Measured on S08, 15 ppm against 6 ppm "
                       "took precursors-with-no-candidate from 155 to 410 while scan positions "
                       "evaluated fell 15%. Tune the two together or not at all.",
                       false, true);
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
    registerStringOption_("picker", "<mode>", "coelution",
                          "Which transition finder to use. coelution: ODIA's port of DIA-NN's "
                          "Searcher::peaks -- a position is a peak only where the fragments "
                          "already correlate. amplitude: local maxima on the summed trace. "
                          "openswath: OpenMS's PeakPickerChromatogram, the reference "
                          "implementation, which picks by amplitude and sets boundaries by "
                          "signal-to-noise, so co-elution enters only afterwards as a score. "
                          "union: both -- the co-elution candidates plus the amplitude ones the "
                          "gate refused, with the co-elution sum attached to every candidate as "
                          "a number rather than used as a threshold. union_openswath: the same "
                          "hybrid with OpenMS's picker supplying the second set.",
                          false);
    setValidStrings_("picker",
                     {"coelution", "amplitude", "openswath", "union", "union_openswath"});
    registerDoubleOption_("openswath_sn", "<ratio>", 1.0,
                          "Signal-to-noise threshold for -picker openswath. OpenMS defaults to "
                          "1.0; OpenSwathWorkflow commonly runs 0.1 on DIA, where a peak sits "
                          "on far more background than in targeted MRM.", false, true);
    registerDoubleOption_("openswath_peak_width", "<seconds>", 0.0,
                          "Expected peak width for -picker openswath, or 0 for its default. "
                          "S08's peaks are ~20-30 s.", false, true);
    registerFlag_("openswath_gauss",
                  "Gaussian rather than Savitzky-Golay smoothing in -picker openswath.", true);
    registerFlag_("collect_mass_residuals",
                  "Keep the m/z deviation of every matched peak and report it per peak group "
                  "as Mass.Ppm. The deviation is computed anyway to test the match and has "
                  "always been discarded. This is the only way to measure the run's real "
                  "fragment mass error: a standalone probe asks whether SOME peak lies within "
                  "tolerance near a time, and on a mostly-absent library the answer is yes by "
                  "coincidence -- an RT-shifted control produced a LARGER apparent offset than "
                  "the true apex. Matches kept here are constrained by co-elution and, at "
                  "q<=0.01, by the whole discriminant. Costs TWO extra float planes per live "
                  "block -- an intensity-weighted sum and its denominator, because a single "
                  "plane could not tell 'no peak matched' from 'matched at exactly 0.000 ppm' "
                  "under the default Sum aggregation. The decode path, not the live blocks, is "
                  "what sets this tool's memory floor, so the cost is real but not binding.",
                  true);
    registerStringOption_("im_features", "<mode>", "auto",
                          "Measure the run's OBSERVED 1/K0 per candidate from intensity-"
                          "weighted mobility planes, and feed var_im_delta from it. Until "
                          "this existed the feature was all-NaN on every run and was dropped "
                          "as carrying no information -- nothing in the pipeline measured "
                          "observed mobility per precursor. 'auto' collects wherever the run "
                          "has ion mobility and is inert where it does not; 'off' restores "
                          "the previous behaviour.",
                          false);
    setValidStrings_("im_features", {"auto", "off"});
    registerStringOption_("mass_anchors", "<mode>", "off",
                          "Harvest one fragment mass residual per contributing fragment of "
                          "every retained candidate, so a mass model can be fitted from "
                          "IDENTIFICATIONS rather than from the probe that runs before pass 1. "
                          "'measure' collects them and reports what a model fitted from them "
                          "looks like, without applying it -- which is the whole experiment: "
                          "the probe's model is what is in force today, and whether "
                          "ID-derived anchors beat it on HELD-OUT fragments is unmeasured. "
                          "Implies -collect_mass_residuals, since the per-fragment deviations "
                          "are what it reads. See doc/15.",
                          false);
    setValidStrings_("mass_anchors", {"off", "measure"});
    registerStringOption_("out_mass_anchors", "<file>", "",
                          "Write the harvested fragment mass residuals here as a TSV: mz, rt, "
                          "ppm (against the UNCORRECTED theoretical m/z), intensity, im, decoy, "
                          "group, qvalue. Fitting a calibration is a modelling question and "
                          "re-running a 10-minute search to try another model is the wrong loop.",
                          false);
    registerIntOption_("max_mass_anchors", "<n>", 8000000,
                       "Ceiling on harvested anchors, ~24 B each. Hitting it truncates the "
                       "sample in RUN ORDER, which is a retention-time bias, so the number "
                       "dropped is reported rather than absorbed.",
                       false);
    registerStringOption_("mass_width_from_ids", "<mode>", "measure",
                          "Size pass 2's fragment window from pass 1's own identifications, "
                          "instead of falling back to a constant when the mass calibration gate "
                          "fails. 'off' does not measure it; 'measure' measures and logs it but "
                          "keeps the constant; 'apply' uses it. Measuring is free -- pass 1 "
                          "already extracts wide and already scores -- so the default measures, "
                          "and applying is opt-in until it has been measured on both "
                          "instruments. Implies -collect_mass_residuals unless 'off', so the "
                          "default carries that option's two extra float planes; 'off' restores "
                          "the smaller footprint.",
                          false, true);
    setValidStrings_("mass_width_from_ids", {"off", "measure", "apply"});
    registerStringOption_("ablate", "<names>", "",
                          "Comma-separated sub-scores to withhold from the classifier, by the "
                          "name they carry in the output (e.g. var_mass_spread). Leave-one-out "
                          "ablation within ONE binary: comparing two builds also compares "
                          "everything else that changed between them, which is how a feature "
                          "gets credited with somebody else's gain.", false, true);
    registerDoubleOption_("train_fdr_initial", "<q>", 0.15,
                          "FDR for the semi-supervised loop's FIRST training-set selection. "
                          "Deliberately lenient: the first pass is seeded by one feature, so a "
                          "strict cut selects too few positives to fit anything and the model "
                          "never bootstraps. pyprophet uses 0.15 for the same reason. Raise it "
                          "on a low-prior library where the loop reports 0.", false, true);
    registerDoubleOption_("train_fdr", "<q>", 0.05,
                          "FDR for the loop's later iterations, once a real discriminant exists.",
                          false, true);
    registerIntOption_("classifier_iterations", "<n>", 3,
                       "Semi-supervised iterations.", false, true);
    registerFlag_("use_pi0",
                  "Storey pi0 correction in the q-value. OFF is the honest/conservative setting "
                  "-- a nominal 1% is a true 1%. ON matches pyprophet and DIA-NN, which report "
                  "more identifications at a nominal 1% that is nearer 2% in truth. It is "
                  "exposed because parity with their NUMBERS and parity with their CALIBRATION "
                  "are different goals and the choice should be visible.", true);
    registerStringOption_("rt_refine", "<mode>", "auto",
                          "Refine retention-time prediction from THIS RUN's own identifications, "
                          "in process. auto: fit after pass 1 and apply if it helps on held-out "
                          "sequences; off: do not. "
                          "\n\nThe map is a monotone function of the library's iRT and cannot beat "
                          "the ordering it is given -- measured on Astral, the best monotone "
                          "calibration reaches SD 38.6 s held out and no span, threshold or bin "
                          "count goes below it, while DIA-NN's own residual is 29.29 s. What is "
                          "left is elution-ORDER error, and only something that reads the SEQUENCE "
                          "can touch it. "
                          "\n\nRidge on amino-acid composition, length, charge and the calibrated "
                          "iRT: 38.55 -> 33.19 s held out by stripped sequence on Astral. The "
                          "peptdeep fine-tuner reaches 27.70 s and is better, but needs torch, a "
                          "separate environment and a GPU node, so it stays behind "
                          "-repredict_irt. This one has no dependency and runs by default, "
                          "because a refinement nobody enables is a refinement nobody gets. "
                          "\n\nmap_only: iterate the ANCHORS without the sequence model -- "
                          "re-derive anchors from the current scoring, refit the map, re-score, "
                          "repeat. That half is separable and on S08 it is the half that works: "
                          "the window for 99% coverage fell 913.7 -> 281.0 s between rounds purely "
                          "from better anchors, while the ridge improved the median and WIDENED "
                          "the p99 that sets the window.",
                          false, true);
    setValidStrings_("rt_refine", {"auto", "map_only", "off"});
    registerDoubleOption_("rt_converge_tol", "<seconds>", 0.1,
                          "Stop refining when a round improves the median |residual| by less "
                          "than this. "
                          "\n\nA DELTA criterion, not an absolute one, and that is the point: "
                          "'is the median below X' cannot be answered without knowing what X "
                          "should be, and X differs by gradient, instrument and library -- the "
                          "same mistake as a fixed m/z window. 'Has it stopped moving' needs no "
                          "such constant and is what convergence actually means. "
                          "\n\n0.1 s against a median that starts near 23 s: more than two "
                          "orders of magnitude below the signal, so the loop stops when the model "
                          "has nothing left to add rather than when it hits a number somebody "
                          "chose. Lowered from 0.5 s, which stopped the equivalent epoch curve at "
                          "125 -- about two points of +/-30 s coverage too early.",
                          false, true);
    registerDoubleOption_("rt_converge_rel", "<fraction>", 0.01,
                          "Stop refining when a round improves the median |residual| by less than "
                          "this FRACTION of it. Runs at both ends of the residual scale then get "
                          "the same criterion: 0.5 s is a real gain on a 5 s median and noise on "
                          "a 200 s one. Whichever of -rt_converge_tol and this triggers first "
                          "stops the loop.", false, true);
    registerIntOption_("rt_refine_rounds", "<n>", 10,
                       "Iterations of the per-run retention-time refinement. Each round derives "
                       "its anchors from the CURRENT scoring, refines, refits the map and "
                       "re-scores -- so a better axis finds better anchors, which is the whole "
                       "point and what a single pass cannot do. Stops early when a round does not "
                       "increase identifications. DIA-NN runs twelve for the same reason. "
                       "\n\nThis is the CAP, not the schedule: -rt_converge_tol and "
                       "-rt_converge_rel stop the loop when the median residual stops moving, "
                       "which is normally well before the cap.",
                       false, true);
    registerOutputFile_("rt_refine_model_out", "<file>", "",
                        "Write the fitted retention-time refinement so it can be reused on OTHER "
                        "runs, skipping the per-run fit. The file records what it was fitted on. "
                        "\n\nThe hazard is specific and this project has already paid it: library "
                        "v5 was built with a model tuned on a different run and cost 2,027 "
                        "confident precursors. A SERIES sharing gradient, instrument and method is "
                        "the case this is for; crossing any of the three is the case it is not.",
                        false, true);
    setValidFormats_("rt_refine_model_out", {"txt"}, false);
    registerInputFile_("rt_refine_model_in", "<file>", "",
                       "Apply a retention-time refinement written by -rt_refine_model_out instead "
                       "of fitting one from this run. Skips the fit entirely, so it works on runs "
                       "with too few identifications to fit their own -- which is the point.",
                       false, true);
    setValidFormats_("rt_refine_model_in", {"txt"}, false);
    registerFlag_("repredict_irt",
                  "Re-predict the SUPPLIED library's iRT with -rt_model, instead of using the "
                  "values the library file carries. "
                  "\n\nThis is the link that was missing between the RT fine-tuner and a real "
                  "search. -rt_model has only ever been consumed when ODIA GENERATES a library "
                  "from FASTA, so a model fine-tuned on a run's own identifications could not be "
                  "applied to a search of that run against a supplied TSV -- which is every "
                  "benchmark we have. "
                  "\n\nMeasured on Astral, held out by stripped sequence with the recipe chosen "
                  "on a separate validation split: the library's own iRT through our map gives "
                  "SD 40.03 s, and a model tuned on 5,867 of this run's identifications for 200 "
                  "epochs gives 27.70 s. DIA-NN's own residual on this file is 29.29 s. "
                  "\n\nThe model must be tuned on THIS run. Reusing one across runs cost 2,027 "
                  "precursors when it happened, and the sidecar says so; -rt_model already warns "
                  "when a provenance file carries that warning.", true);
    registerStringOption_("mass_features", "<mode>", "auto",
                          "Whether var_mass_accuracy and var_mass_spread reach the classifier. "
                          "auto: only when the fragment mass calibration gate FAILED. "
                          "\n\nMeasured 2026-08-09, and the two files disagree completely. "
                          "Astral, where the gate fails and the window is 50 ppm and uncentred: "
                          "5,025 -> 5,729, the largest single gain measured on that file. S08, "
                          "where the gate passes and the window is 10 ppm centred on -9.35 with "
                          "a per-fragment sigma of 1.10: 1,342 -> 1,168, and leave-one-out says "
                          "both are harmful there (spread alone 1,315, accuracy alone 1,285). "
                          "\n\nThat is not a contradiction, it is the rule: A FEATURE IS WORTH "
                          "WHAT THE EXTRACTION HAS NOT ALREADY SPENT. A centred 10 ppm window "
                          "has already used the mass information as a filter, so the residual is "
                          "noise for a classifier to overfit; a 50 ppm uncentred window has "
                          "spent none of it, and mass deviation is then the one thing separating "
                          "a real fragment from something that merely co-elutes. "
                          "\n\nThe gate's own verdict is therefore the right switch, and it is "
                          "measured from the data rather than set per file. In pass 1 the gate "
                          "has not run yet, so auto leaves them ON; pass 2 is where the decision "
                          "is real.", false, true);
    setValidStrings_("mass_features", {"auto", "on", "off"});
    registerFlag_("no_match_decoy_n",
                  "TURN OFF the decoy candidate-count matching described below, which is on by "
                  "default. Named for what the FLAG does, not for what the feature does: a flag "
                  "whose help opens by praising the thing it disables is how a user ends up "
                  "setting it to get the behaviour it removes. "
                  "\n\nThe feature: draw each decoy's best score from as many candidates as a TARGET has, "
                  "instead of from all of its own. Target-decoy competition assumes the two "
                  "classes are exchangeable and they are not: the picker keeps candidates within "
                  "-max_corr_diff of a precursor's OWN best, so a precursor with no real peak "
                  "admits nearly every position. Measured on Astral, targets carry 10.98 "
                  "candidates per precursor and decoys 20.99, so a q-value compares best-of-11 "
                  "against best-of-21 over 9,453 decoy precursors -- and ODIA's own FDR then "
                  "refused 3,834 Astral precursors that are all in DIA-NN's truth set. Caps are "
                  "quantile-matched to the target distribution and applied in canonical order, "
                  "never by score. "
                  "\n\nON BY DEFAULT since 2026-08-10; this flag turns it OFF. It was worth +653 "
                  "on Astral and +122 on S08, and the entrapment fixture measures it CONSERVATIVE "
                  "-- FDP 0.0020 against a claimed 0.01 -- so leaving it opt-in meant every "
                  "default run paid for an asymmetry we know how to remove.", true);
    registerDoubleOption_("mass_sigma_multiple", "<n>", 8.0,
                          "Fragment window as this many robust sigmas of the corrected "
                          "residual, when the mass calibration fits one. Swept 1-10 on both "
                          "benchmark files: S08 plateaus from 4, Astral never turns over, 8 "
                          "maximises the worst case. PROVISIONAL -- on Astral the "
                          "uncalibrated 50 ppm fallback still beats every swept width "
                          "(2,078 identifications against 1,222 at k=10), and the sigma this "
                          "multiplies is the probe's, which differs 4x between the two files "
                          "where the true per-fragment sigma differs by 4%. See "
                          "MassCalibration.h.",
                          false);
    registerDoubleOption_("mass_width_sigmas", "<n>", 3.0,
                          "Half-width, in robust sigmas of the per-fragment deviation, for "
                          "-mass_width_from_ids apply. 3 covers 99.7% of a Gaussian; the "
                          "distribution has heavier tails than that, which is why the result is "
                          "also floored by -mass_width_min_ppm.", false, true);
    registerDoubleOption_("mass_width_min_ppm", "<ppm>", 5.0,
                          "Floor for -mass_width_from_ids apply. A pass 1 that identified only "
                          "its cleanest precursors measures their scatter, not the run's.",
                          false, true);
    registerIntOption_("mass_width_min_groups", "<n>", 200,
                       "Accepted precursors below which the measured width is not trusted.",
                       false, true);
    registerDoubleOption_("ms1_im_scale", "<factor>", 2.0,
                          "MS1 mobility half-width, as a multiple of -precursor_im_window. "
                          "MS1 ions are not mobility-selected by an isolation window, so the "
                          "precursor's MS1 mobility spread is wider than its fragments'. "
                          "Measured on S08/lib_targets: 2.0 gives 1306, 1.0 gives 1230. Set 0 "
                          "to disable the MS1 mobility gate entirely.",
                          false, true);
    registerFlag_("no_ms1",
                  "Do not read MS1 or compute var_ms1_coelution. MS1/MS2 co-elution is the "
                  "only sub-score that does not read MS2 fragment traces, and so the only one "
                  "a co-eluting interferent cannot corrupt along with the rest; measured at "
                  "13.7x enrichment in its top bin and 1.4-1.5x in bulk on S08 + v6_50k. Use "
                  "this to A/B it or on a run whose MS1 is not worth the pass.",
                  true);
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
    // MS1 traces, once per run, before the first extraction that will score.
    //
    // Built here rather than inside the extractor: the streaming extractor is
    // organised around isolation windows and MS2 cycles, MS1 has neither, and
    // its memory behaviour is the one part of this pipeline that is measured.
    // See Ms1Traces for why co-elution is the only MS1 quantity worth carrying.
    if (ms1_traces_.empty() && !getFlag_("no_ms1"))
    {
      const auto t0 = std::chrono::steady_clock::now();
      // Width RELATIVE to the MS2 window, because the right ratio is an
      // empirical question and both of my confident answers were wrong.
      //
      // It shipped at 2x on an unmeasured assumption. I then set it to 1x on a
      // consistency argument -- the two traces being correlated should sample
      // one ion population -- and that cost 76 identifications on
      // S08/lib_targets (1306 -> 1230). The consistency argument is sound about
      // what the correlation MEANS and wrong about what it is worth; MS1 ions
      // are not mobility-selected by an isolation window, so the precursor's
      // MS1 mobility spread is genuinely wider than its fragments'.
      //
      // The NaN handling is NOT part of this knob and stays fixed: a precursor
      // with no library 1/K0 is ungated, matching MS2, rather than having every
      // peak rejected.
      ms1_traces_ = ODIA::Ms1Traces::build(library, *source, options.fragment_ppm,
                                           options.precursor_im_window *
                                             getDoubleOption_("ms1_im_scale"));
      const double secs = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
      if (ms1_traces_.empty())
      {
        writeLogInfo_("MS1: the run carries no MS1 spectra, so var_ms1_coelution "
                      "will be absent rather than zero");
      }
      else
      {
        writeLogInfo_(ms1_traces_.describe() + ", in " + std::to_string(secs) + " s");
      }
    }

    applyMassCalibration_(library, *source, options, rt_window_override);
    // The harvest records deviations against the target the extractor actually
    // searched, which is the corrected one. Hand it the correction so it can
    // report against the uncorrected theoretical instead -- otherwise a model
    // fitted from the anchors is an increment to this one, and this one gets
    // charged twice when the two are compared. Safe to do here and nowhere
    // earlier: the coefficients were unknown when the Sink was built.
    if (auto* pgs = dynamic_cast<ODIA::PeakGroupScorer::Sink*>(&sink))
    {
      pgs->setAppliedMassCorrection(options.fragment_ppm_offset,
                                    options.fragment_ppm_log_slope,
                                    options.fragment_ppm_slope_per_1000,
                                    options.fragment_ppm_ref_mz);
    }
    options.rt_window_seconds = rt_window_override != 0.0 ? rt_window_override
                                                          : getDoubleOption_("rt_window");
    options.max_precursors = static_cast<std::size_t>(
      std::max(0, getIntOption_("max_precursors")));
    options.use_ion_mobility = !getFlag_("no_ion_mobility");
    options.precursor_im_window = getDoubleOption_("precursor_im_window");
    // The width measurement reads the ppm planes, so asking for it turns them
    // on. Making the user pass two flags that only work together is a way of
    // producing runs that silently measured nothing.
    options.collect_mass_residuals = getFlag_("collect_mass_residuals") ||
                                     getStringOption_("mass_width_from_ids") != "off" ||
                                     getStringOption_("mass_anchors") != "off";
    // The mobility planes are what make IM_DELTA computable at all -- it has
    // been all-NaN on every run since it was added, because nothing measured
    // observed 1/K0 per precursor. Inert on a run without mobility, so the
    // default is on: an absent feature costs more than two float planes.
    options.collect_im_residuals = getStringOption_("im_features") != "off";
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

    // An EXTERNAL iRT map must be worth as much as a fitted one.
    //
    // -irt_slope/-irt_intercept centred the extraction and stopped there: the
    // library kept iRT units, `library_rt_is_run_seconds` stayed false, and
    // RT_DELTA -- gated on exactly that flag -- was NaN for every row and then
    // zeroed by the constant-column guard. Measured on v6_50k: supplying a map
    // fitted from DIA-NN's own anchors (p50 residual 9.1 s) still left
    // var_rt_delta at 0.0% non-zero, so the run discarded the one feature the
    // map exists to enable.
    //
    // Converting the library up front makes the external path identical to the
    // internal one, which rewrites irt in place and sets the flag after pass 1.
    {
      const double sl = getDoubleOption_("irt_slope");
      const double ic = getDoubleOption_("irt_intercept");
      if (sl != 0.0)
      {
        auto& irt = library.precursors().irt;
        std::size_t n = 0;
        for (std::size_t i = 0; i < irt.size(); ++i)
        {
          if (!std::isfinite(irt[i])) { continue; }
          irt[i] = static_cast<float>(sl * irt[i] + ic);
          ++n;
        }
        external_irt_ = true;
        scoring_rt_is_run_seconds_ = true;
        writeLogInfo_("applied the supplied iRT map to " + std::to_string(n) +
                      " precursors: the library now carries RUN SECONDS, so "
                      "var_rt_delta is available. An external map is now worth "
                      "what a fitted one is.");
      }
    }

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
        const auto rc = extractAndScore_(library, run, 0.0, external_irt_, scored);
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
      // external_irt_, NOT false. When -irt_slope/-irt_intercept were supplied
      // the library's irt was already converted to run seconds above, and
      // extractInto_ reads this flag to decide whether to apply that map --
      // so passing false applied it a SECOND time. On S08 with the stored v6
      // map that is 7.77 * (run seconds) + 733, i.e. ~24,000 s for a 5,400 s
      // gradient: pass 1 extracted from beyond the end of the run and every
      // anchor it could have found was unreachable.
      const auto rc = extractAndScore_(library, run,
                                       pass1_window > 0.0 ? pass1_window : 1.0e9,
                                       external_irt_, pass1);
      if (rc != EXECUTION_OK) { return rc; }
    }

    // The fragment window, measured rather than defaulted.
    //
    // Pass 1 has just extracted wide and scored, so the per-fragment deviations
    // of everything the FDR accepted are sitting in `pass1` at no extra cost.
    // That population is what `MassCalibration`'s standalone probe could never
    // get: its matches are constrained by co-elution and by the whole
    // discriminant, where a probe over a mostly-absent library measures
    // coincidence. Measured always, applied only under -mass_width_from_ids
    // apply -- the measurement and its application are different claims and the
    // last time they were coupled a passing gate cost 56% of the run.
    if (getStringOption_("mass_width_from_ids") != "off")
    {
      mass_width_ = ODIA::MassWidth::measure(
        pass1, 0.01, extracted_ppm_, getDoubleOption_("mass_width_sigmas"),
        static_cast<std::size_t>(std::max(1, getIntOption_("mass_width_min_groups"))));
      writeLogInfo_(ODIA::MassWidth::report(mass_width_));
    }


    // The step-1 experiment of doc/15: are anchors taken from IDENTIFICATIONS a
    // better basis for the mass model than the probe that runs before pass 1?
    // Measured, never applied -- the probe's model is what is in force, and
    // this reports what would change if it were replaced.
    if (getStringOption_("mass_anchors") != "off")
    { reportMassAnchors_(pass1, library); }

    std::ostringstream p1;
    p1 << "pass 1: " << pass1.groups.size() << " peak groups, "
       << pass1.target_groups << " target / " << pass1.decoy_groups << " decoy, "
       << pass1.iterations_trained << " iterations trained";
    writeLogInfo_(p1.str());

    // Anchors: the best group of each confidently identified target, paired
    // with the library RT it came from.
    const double anchor_q = getDoubleOption_("anchor_q");
    std::vector<std::pair<double, double>> anchors;
    // (precursor index, observed apex RT) for the retention-time refiner. The
    // anchors themselves cannot serve: they carry the LIBRARY iRT, and the
    // refiner needs the CALIBRATED one, which does not exist until the map
    // below has been applied -- by which point `best` is out of scope.
    std::vector<std::pair<std::size_t, double>> refine_rows;
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
          refine_rows.emplace_back(i, static_cast<double>(best[i]->apex_rt));
        }
      }

      // The anchors themselves, for anything that wants to LEARN from them
      // rather than fit a monotone map through them. Written here because
      // `best` -- which carries the sequence's precursor index -- lives only in
      // this scope, and a (library iRT, observed RT) pair without the sequence
      // is useless to a retention-time model.
      const std::string anchor_path = getStringOption_("out_anchors");
      if (!anchor_path.empty())
      {
        const auto& pr = library.precursors();
        std::ofstream os(anchor_path);
        os << "Precursor.Id\tModified.Sequence\tPrecursor.Charge\tLibrary.iRT"
              "\tObserved.RT\tDScore\tQValue\n";
        std::size_t written = 0;
        for (std::size_t i = 0; i < best.size(); ++i)
        {
          if (best[i] == nullptr || !std::isfinite(original_irt[i])) { continue; }
          const auto seq = library.strings().get(pr.modified_sequence[i]);
          os << seq << '_' << static_cast<int>(pr.charge[i]) << '\t'
             << seq << '\t' << static_cast<int>(pr.charge[i]) << '\t'
             << original_irt[i] << '\t' << best[i]->apex_rt << '\t'
             << best[i]->dscore << '\t' << best[i]->qvalue << '\n';
          ++written;
        }
        writeLogInfo_("wrote " + std::to_string(written) + " anchors to " + anchor_path);
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
    const double loess_span = getDoubleOption_("rt_loess_span");
    const std::string rt_interp = getStringOption_("rt_interpolation");
    auto trafo = ODIA::Calibration::fit(anchors, &p95, loess_span, rt_interp);

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
      // Same span as the real fit, or the probe measures a different model.
      const auto probe = ODIA::Calibration::fit(fit_set, &dummy, loess_span, rt_interp);
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

    // What the calibration actually achieved, before and after, on the same
    // anchors. Reported always -- the p95 above says how tight the fit is but
    // not how much of the error it REMOVED, and those are different questions.
    //
    // Both a standard deviation and a robust sigma, deliberately. The residual
    // distribution is heavy-tailed by construction (pass 1's anchors come from
    // windows narrower than the RT error, so gross outliers are guaranteed),
    // and an SD over that is dominated by the tail rather than describing the
    // bulk. Quoting only the SD would understate the calibration; quoting only
    // the robust figure would hide the tail that sets the pass-2 window. The
    // gap between them IS the diagnostic.
    reportRtResiduals_(anchors, trafo);

    // Per-run retention-time refinement, IN iRT SPACE and ITERATIVELY.
    //
    // One-shot refinement measured on the anchors it was fitted to is circular:
    // the whole benefit of a better iRT is that it finds BETTER ANCHORS, and a
    // single pass never lets it. So each round re-derives the anchors from the
    // CURRENT scoring, refines against those, refits the map, and re-scores --
    // which is what DIA-NN does for twelve rounds and for the same reason.
    //
    // It costs no extraction. Pass 1's peak groups are already in memory and
    // the candidate picker is retention-time agnostic, so a new axis changes
    // RT_DELTA and nothing else; `PeakGroupScorer::refit` recomputes that one
    // column and refits the discriminant.
    if (getStringOption_("rt_refine") != "off" && !refine_rows.empty())
    {
      const int rounds = std::max(1, getIntOption_("rt_refine_rounds"));
      const double tol_abs = std::max(0.0, getDoubleOption_("rt_converge_tol"));
      const double tol_rel = std::max(0.0, getDoubleOption_("rt_converge_rel"));
      const double anchor_q_r = getDoubleOption_("anchor_q");
      auto refit_options = scoringOptions_();
      refit_options.library_rt_is_run_seconds = true;
      refit_options.disabled_sub_scores = ablatedSubScores_();

      const std::string model_in = getStringOption_("rt_refine_model_in");
      std::size_t best_ids = pass1.identified_at_1pct;
      ODIA::RtRefiner kept;
      // A FIXED evaluation set, frozen at round 1.
      //
      // The anchor set grows every round -- 1,183 -> 1,408 on S08 -- so a p99
      // measured on the current anchors is measured on a different population
      // each time and cannot be compared across rounds. Round 2's "before" p50
      // came out WORSE than round 1's "after" for exactly this reason, which
      // reads as a regression and is a change of denominator.
      //
      // Freezing the evaluation set makes the round-to-round numbers mean
      // something. It is still an anchor population and still self-referenced;
      // that limitation is reported separately and is not fixable from inside
      // the loop.
      std::vector<std::pair<double, double>> eval_set;
      double best_p50 = std::numeric_limits<double>::infinity();

      for (int round = 1; round <= rounds; ++round)
      {
        // ANCHORS FROM THE CURRENT SCORING, not from round zero.
        std::vector<const ODIA::PeakGroupScorer::PeakGroup*> cur(
          library.precursorCount(), nullptr);
        for (const auto& g : pass1.groups)
        {
          if (g.decoy || g.qvalue > anchor_q_r) { continue; }
          auto*& bp = cur[g.precursor];
          if (bp == nullptr || g.dscore > bp->dscore) { bp = &g; }
        }
        std::vector<std::pair<double, double>> a_r;
        std::vector<std::string> seqs;
        std::vector<int> chg;
        std::vector<double> lib_irt, target_irt;
        const auto& pr = library.precursors();
        for (std::size_t i = 0; i < cur.size(); ++i)
        {
          if (cur[i] == nullptr || !std::isfinite(original_irt[i])) { continue; }
          a_r.emplace_back(static_cast<double>(original_irt[i]),
                           static_cast<double>(cur[i]->apex_rt));
          seqs.emplace_back(library.strings().get(pr.modified_sequence[i]));
          chg.push_back(static_cast<int>(pr.charge[i]));
          lib_irt.push_back(static_cast<double>(original_irt[i]));
        }
        if (a_r.size() < 250) { break; }
        if (eval_set.empty()) { eval_set = a_r; }

        double p95_r = 0.0;
        const auto map_r = ODIA::Calibration::fit(a_r, &p95_r, loess_span, rt_interp);
        double ilo = a_r.front().first, ihi = a_r.front().first;
        for (const auto& q : a_r)
        { ilo = std::min(ilo, q.first); ihi = std::max(ihi, q.first); }
        for (const auto& q : a_r)
        { target_irt.push_back(ODIA::Calibration::invertAt(map_r, q.second, ilo, ihi)); }

        ODIA::RtRefiner refiner;
        ODIA::RtRefiner::Report rep;
        if (!model_in.empty() && round == 1)
        {
          std::string prov;
          if (!refiner.load(model_in, &prov))
          {
            writeLogError_("-rt_refine_model_in " + model_in + " is not a refinement model.");
            return INTERNAL_ERROR;
          }
          rep.fitted = true;
          writeLogWarn_("retention-time refinement LOADED from " + model_in +
                        " and NOT fitted on this run" +
                        (prov.empty() ? "" : " [" + prov + "]") + ".");
        }
        else if (!model_in.empty()) { break; }
        else if (getStringOption_("rt_refine") == "map_only")
        {
          // No sequence model: the round is anchors -> map -> re-score. The
          // refined axis IS the original one, so the refit below is a refit of
          // the same map on a better anchor set, which is the point.
          rep.fitted = true;
        }
        else { rep = refiner.fit(seqs, chg, lib_irt, target_irt); }
        if (!rep.fitted)
        {
          writeLogInfo_("retention-time refinement round " + std::to_string(round) +
                        ": not applied -- " + rep.note);
          break;
        }

        // Refine the axis and REFIT the map on it, so monotonicity is
        // re-established rather than assumed to survive the correction.
        library.precursors().irt = original_irt;
        if (getStringOption_("rt_refine") != "map_only") { refiner.apply(library); }
        std::vector<float> refined = library.precursors().irt;
        std::vector<std::pair<double, double>> refit_pts;
        refit_pts.reserve(a_r.size());
        {
          std::size_t k = 0;
          for (std::size_t i = 0; i < cur.size(); ++i)
          {
            if (cur[i] == nullptr || !std::isfinite(original_irt[i])) { continue; }
            refit_pts.emplace_back(static_cast<double>(refined[i]), a_r[k].second);
            ++k;
          }
        }
        double p95_ref = 0.0;
        const auto map_ref = ODIA::Calibration::fit(refit_pts, &p95_ref, loess_span, rt_interp);

        // Apply and RE-SCORE, then judge on identifications -- the quantity we
        // actually want -- rather than on a residual over the anchors the model
        // was just fitted to.
        auto& irt_r = library.precursors().irt;
        for (std::size_t i = 0; i < irt_r.size(); ++i)
        {
          if (std::isfinite(refined[i]))
          { irt_r[i] = static_cast<float>(map_ref.apply(static_cast<double>(refined[i]))); }
        }
        ODIA::PeakGroupScorer::refit(library, pass1, refit_options);

        // WHAT THIS IS FOR, and it is not identifications.
        //
        // The calibration exists so extraction does not MISS a peak and the
        // window it opens is CLEAN. Those are two numbers, not one:
        //
        //   COVERAGE  -- the half-width that contains a given fraction of true
        //                apices. Miss this and the peak is not in the window at
        //                any score.
        //   WIDTH     -- what that half-width costs. Every extra second admits
        //                interference into every trace, which is what makes a
        //                feature dirty and what the m/z recalibration downstream
        //                then has to survive.
        //
        // So the objective is the WIDTH REQUIRED FOR A FIXED COVERAGE. A
        // refinement that narrows it wins; one that widens it loses, whatever it
        // does to the identification count today.
        const auto quantiles = [&](const std::vector<std::pair<double, double>>& pts,
                                   const OpenMS::TransformationDescription& m) {
          std::vector<double> a;
          a.reserve(pts.size());
          for (const auto& q : pts) { a.push_back(std::fabs(q.second - m.apply(q.first))); }
          std::sort(a.begin(), a.end());
          const auto at = [&](double f) {
            return a[std::min(a.size() - 1,
                              static_cast<std::size_t>(f * (a.size() - 1)))];
          };
          return std::array<double, 4>{at(0.50), at(0.95), at(0.99), a.back()};
        };
        // Both maps scored on the SAME frozen set, in the axis each expects.
        const auto q_before = quantiles(eval_set, map_r);
        std::vector<std::pair<double, double>> eval_ref;
        eval_ref.reserve(eval_set.size());
        if (getStringOption_("rt_refine") == "map_only") { eval_ref = eval_set; }
        else
        {
          // The refined axis for the evaluation rows: apply the same correction
          // the library got, so the two maps are compared on one population.
          ODIA::Library& lib_ref = library;
          for (std::size_t i = 0, k = 0; i < lib_ref.precursorCount() && k < eval_set.size(); ++i)
          {
            if (!std::isfinite(original_irt[i])) { continue; }
            if (std::fabs(static_cast<double>(original_irt[i]) - eval_set[k].first) < 1e-9)
            { eval_ref.emplace_back(static_cast<double>(refined[i]), eval_set[k].second); ++k; }
          }
          if (eval_ref.size() != eval_set.size()) { eval_ref = refit_pts; }
        }
        const auto q_after = quantiles(eval_ref, map_ref);

        std::ostringstream ro;
        ro.setf(std::ios::fixed);
        ro.precision(2);
        ro << "rt refine round " << round << ": " << a_r.size() << " anchors ("
           << rep.trimmed << " trimmed)"
           << "\n  |residual| p50 " << q_before[0] << " -> " << q_after[0]
           << " s, p95 " << q_before[1] << " -> " << q_after[1]
           << " s, p99 " << q_before[2] << " -> " << q_after[2]
           << " s, max " << q_before[3] << " -> " << q_after[3] << " s"
           << "\n  window for 99% coverage: " << (2.0 * q_before[2]) << " s -> "
           << (2.0 * q_after[2]) << " s"
           << "   [identifications " << best_ids << " -> " << pass1.identified_at_1pct
           << ", not the criterion]";
        writeLogInfo_(ro.str());

        // p50, NOT p99 -- measured, and the reverse of what this said first.
        //
        // The anchor p99 is contamination: on the anchors it reads ~455 s, and
        // on DIA-NN's 10,891 confident precursors the true p99 is 108 s. So a
        // criterion built on it was reading misassignment, not calibration, and
        // it rejected the refinement that actually works.
        //
        // External coverage over all 10,891 (the population that includes the
        // ~8,600 we never anchor):
        //
        //   no refinement   +/-30s 60.26%   +/-60s 88.46%   p50 23.5  p95 76.7
        //   map_only        +/-30s 60.61%   +/-60s 88.72%   p50 23.4  p95 77.2
        //   ridge           +/-30s 71.21%   +/-60s 93.62%   p50 18.1  p95 66.5
        //
        // The refinement improves the BULK and leaves the tail, which is why
        // p99 could not see it and why p50 can. The anchor p50 tracked the
        // external p50 and p95 correctly in both rounds (22.3 -> 15.7,
        // 19.5 -> 13.7), so it is the honest runtime proxy for a number we
        // cannot compute without truth.
        // ACROSS rounds, on the frozen set. Within a round, `map_only` changes
        // nothing by construction and a within-round test rejects it at round 1
        // -- which is exactly what happened, and why its iteration never ran.
        // Its gain is between rounds: better anchors, better map.
        // CONVERGENCE ON THE DELTA, not on an absolute threshold.
        //
        // Three outcomes, and they are different. A round can make things
        // WORSE, in which case roll back. It can improve by so little that the
        // next round will not repay its cost, in which case accept and stop.
        // Or it can still be moving, in which case continue.
        const double improvement = best_p50 - q_after[0];
        const bool worse = !(q_after[0] < best_p50);
        const bool converged = !worse &&
          (improvement < tol_abs || improvement < tol_rel * best_p50);

        if (worse)
        {
          library.precursors().irt = original_irt;
          if (!refined_irt_.empty()) { library.precursors().irt = refined_irt_; }
          writeLogInfo_("rt refine: round " + std::to_string(round) +
                        " made the median residual worse; keeping the previous axis");
          break;
        }
        best_p50 = q_after[0];
        best_ids = pass1.identified_at_1pct;
        refined_irt_ = refined;
        trafo = map_ref;
        p95 = p95_ref;
        kept = refiner;

        if (converged)
        {
          std::ostringstream co;
          co.setf(std::ios::fixed);
          co.precision(3);
          co << "rt refine: CONVERGED after round " << round << " -- the median residual "
             << "improved by " << improvement << " s, below -rt_converge_tol " << tol_abs
             << " s and " << (100.0 * tol_rel) << "% of " << (best_p50 + improvement) << " s";
          writeLogInfo_(co.str());
          break;
        }
      }

      const std::string model_out = getStringOption_("rt_refine_model_out");
      if (!model_out.empty() && model_in.empty() && kept.fitted())
      {
        if (kept.save(model_out, "fitted iteratively on " + run))
        { writeLogInfo_("wrote the retention-time refinement to " + model_out); }
      }
    }

    // Applied to the ORIGINAL values, for the reason above.
    //
    // After this the library's `irt` holds RUN SECONDS, not normalised iRT,
    // which is what lets the mass probe be re-measured against it below with
    // an identity map.
    rt_map_fitted_ = true;
    auto& irt = library.precursors().irt;
    // The refined axis when refinement was accepted, the original otherwise.
    const std::vector<float>& source = refined_irt_.empty() ? original_irt : refined_irt_;
    for (std::size_t i = 0; i < irt.size(); ++i)
    {
      if (std::isfinite(source[i]))
      {
        irt[i] = static_cast<float>(trafo.apply(static_cast<double>(source[i])));
      }
    }

    // Per-run retention-time refinement, in process.
    //
    // HERE, and not earlier: the library's irt now holds RUN SECONDS, so the
    // anchors' x-values and the precursors' predictions are on one axis and the
    // model learns a correction to the CALIBRATED value rather than to a
    // normalised iRT the run has never seen.
    //
    // The refiner holds itself out by stripped sequence and REFUSES to apply a
    // model that does not beat the calibration on those rows -- so this can run
    // unconditionally and leave the axis alone when refinement is not warranted.
    // -stop_after calib: the calibration IS the deliverable here, and pass 2
    // would cost a full extraction to tell us nothing more about it. Placed
    // after the library's irt has been rewritten, so a caller that also asked
    // for -out_lib gets the CALIBRATED library rather than the original.
    if (getStringOption_("stop_after") == "calib")
    {
      writeLogInfo_("-stop_after calib: the retention-time map is fitted and its "
                    "residuals are reported above; not running pass 2");
      return EXECUTION_OK;
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


  /// Residuals of the fitted retention-time map, before and against it.
  ///
  /// `anchors` are (library iRT, observed apex RT). "before" is the library
  /// value against the observation on whatever axis the library arrived in;
  /// it is only meaningful when the library is already in run seconds (an
  /// -irt_slope/-irt_intercept calibration, or a previous run's map), so it is
  /// labelled as raw rather than presented as a like-for-like improvement.
  void reportRtResiduals_(const std::vector<std::pair<double, double>>& anchors,
                          const OpenMS::TransformationDescription& trafo)
  {
    if (anchors.empty()) { return; }
    std::vector<double> before, after;
    before.reserve(anchors.size());
    after.reserve(anchors.size());
    for (const auto& a : anchors)
    {
      before.push_back(a.second - a.first);
      after.push_back(a.second - trafo.apply(a.first));
    }
    const auto describe = [](std::vector<double> v, const char* label) {
      const std::size_t n = v.size();
      double mean = 0.0;
      for (const double x : v) { mean += x; }
      mean /= static_cast<double>(n);
      double ss = 0.0;
      for (const double x : v) { ss += (x - mean) * (x - mean); }
      const double sd = n > 1 ? std::sqrt(ss / static_cast<double>(n - 1)) : 0.0;
      std::vector<double> s = v;
      std::sort(s.begin(), s.end());
      const double median = s[n / 2];
      std::vector<double> ad;
      ad.reserve(n);
      for (const double x : s) { ad.push_back(std::fabs(x - median)); }
      std::sort(ad.begin(), ad.end());
      // 1.4826 x MAD is the Gaussian-consistent sigma; on a heavy-tailed sample
      // it describes the bulk where the SD describes the tail.
      const double robust = 1.4826 * ad[n / 2];
      std::vector<double> abs_v;
      abs_v.reserve(n);
      for (const double x : s) { abs_v.push_back(std::fabs(x)); }
      std::sort(abs_v.begin(), abs_v.end());
      std::ostringstream os;
      os.setf(std::ios::fixed);
      os.precision(2);
      os << "  " << label << ": n " << n << ", median " << median
         << " s, SD " << sd << " s, robust sigma " << robust
         << " s, p50|e| " << abs_v[n / 2] << " s, p95|e| "
         << abs_v[static_cast<std::size_t>(0.95 * (n - 1))] << " s, max|e| "
         << abs_v.back() << " s";
      return os.str();
    };
    writeLogInfo_("retention-time residuals on the anchors:");
    writeLogInfo_(describe(before, "raw   (library value vs observed)"));
    writeLogInfo_(describe(after,  "mapped(calibrated  vs observed)"));
  }

  /// Every sub-score withheld from the classifier: `-ablate` plus whatever
  /// `-mass_features` decides.
  ///
  /// Evaluated late, because the mass-feature arm depends on the fragment mass
  /// calibration's verdict and that does not exist until extraction has set up.
  std::vector<int> ablatedSubScores_()
  {
    std::vector<int> out;
    const auto& names = ODIA::PeakGroupScorer::subScoreNames();
    const auto index_of = [&](const std::string& n) {
      const auto it = std::find(names.begin(), names.end(), n);
      return it == names.end() ? -1 : static_cast<int>(std::distance(names.begin(), it));
    };

    std::stringstream ss(getStringOption_("ablate"));
    std::string name;
    while (std::getline(ss, name, ','))
    {
      name.erase(0, name.find_first_not_of(" \t"));
      const auto end = name.find_last_not_of(" \t");
      if (end != std::string::npos) { name.erase(end + 1); }
      if (name.empty()) { continue; }
      const int i = index_of(name);
      // Fatal, not ignored. An ablation arm that silently ablated nothing would
      // report the baseline and be read as "the feature is worthless".
      if (i < 0) { throw std::invalid_argument("-ablate: no sub-score is called '" + name + "'"); }
      out.push_back(i);
    }

    const std::string mode = getStringOption_("mass_features");
    const bool mass_on = mode == "on" ||
                         (mode == "auto" && (!mass_model_known_ || !mass_model_.fitted));
    if (!mass_on)
    {
      for (const char* n : {"var_mass_accuracy", "var_mass_spread"})
      {
        const int i = index_of(n);
        if (i >= 0) { out.push_back(i); }
      }
    }
    return out;
  }


  /// Fit the mass model from FDR-accepted identifications and score it on
  /// fragments it never saw, against the probe that is in force today.
  ///
  /// Split by STRIPPED SEQUENCE, not by fragment and not by precursor. Twelve
  /// fragments of one precursor share an apex, a mobility and whatever
  /// interference sits under them, and the same peptide at two charges shares
  /// its chemistry -- a fragment-level split would put near-copies on both
  /// sides and report a held-out number that is really an in-sample one.
  void reportMassAnchors_(const ODIA::PeakGroupScorer::Result& pass1,
                          const ODIA::Library& library)
  {
    if (pass1.mass_anchors.empty())
    {
      writeLogInfo_("mass anchors: none harvested -- either nothing was accepted or "
                    "the extractor did not collect per-fragment deviations");
      return;
    }
    if (pass1.mass_anchors_dropped != 0)
    {
      writeLogWarn_("mass anchors: " + std::to_string(pass1.mass_anchors_dropped) +
                    " residuals dropped at the -max_mass_anchors ceiling. The kept "
                    "sample is truncated in RUN ORDER, so it is biased towards early "
                    "retention times; raise the ceiling before trusting the numbers.");
    }

    const auto& pr = library.precursors();
    const double q = 0.01;   // 1%, not 5% -- doc/15 section 12.5

    std::vector<ODIA::MassResidual> train, test;
    train.reserve(pass1.mass_anchors.size());
    for (const auto& a : pass1.mass_anchors)
    {
      if (a.group >= pass1.groups.size()) { continue; }
      const auto& g = pass1.groups[a.group];
      if (g.decoy || !(g.qvalue <= q)) { continue; }
      if (g.precursor >= pr.modified_sequence.size()) { continue; }

      // Strip modifications so a peptide cannot appear on both sides wearing a
      // different mass. FNV-1a over the stripped residues, same hash the RT
      // refiner splits on, so the two stages agree about what a held-out
      // peptide is.
      const auto seq = library.strings().get(pr.modified_sequence[g.precursor]);
      std::uint32_t h = 2166136261u;
      for (const char c : seq)
      {
        if (c < 'A' || c > 'Z') { continue; }
        h ^= static_cast<std::uint8_t>(c);
        h *= 16777619u;
      }
      ((h % 100u) < 15u ? test : train).push_back(a.residual);
    }

    const std::string dump = getStringOption_("out_mass_anchors");
    if (!dump.empty())
    {
      std::ofstream out(dump);
      if (!out)
      { writeLogWarn_("cannot write -out_mass_anchors " + dump); }
      else
      {
        out << "mz\trt\tppm\tintensity\tim\tdecoy\tgroup\tqvalue\tprecursor\n";
        for (const auto& a : pass1.mass_anchors)
        {
          if (a.group >= pass1.groups.size()) { continue; }
          const auto& g = pass1.groups[a.group];
          out << a.residual.mz << '\t' << a.residual.rt << '\t' << a.residual.ppm << '\t'
              << a.residual.intensity << '\t' << a.residual.im << '\t'
              << (g.decoy ? 1 : 0) << '\t' << a.group << '\t' << g.qvalue << '\t'
              << g.precursor << '\n';
        }
        writeLogInfo_("wrote " + std::to_string(pass1.mass_anchors.size()) +
                      " mass anchors to " + dump);
      }
    }

    std::ostringstream os;
    os << "mass anchors: " << pass1.mass_anchors.size() << " harvested, "
       << train.size() << " train / " << test.size() << " held out by stripped sequence "
       << "from groups at q<=" << q;
    writeLogInfo_(os.str());

    if (train.size() < 200 || test.size() < 200)
    {
      writeLogInfo_("mass anchors: too few to fit and score (need 200 each side) -- "
                    "reporting nothing rather than a number from a handful of peptides");
      return;
    }

    // Defaults deliberately: the question is whether a DIFFERENT ANCHOR SOURCE
    // beats the probe, so every other knob has to stay where the probe had it.
    const ODIA::MassCalibration::Options fo;
    const auto id_model = ODIA::MassCalibration::fit(train, fo, nullptr);

    // Four corrections scored on the SAME held-out fragments. `none` is the
    // denominator; `constant` is what a scalar offset achieves, and is the bar
    // any shape has to clear; `probe` is what is in force today.
    const double none_ppm = ODIA::MassCalibration::systematicResidualPpm(
      test, [](double) { return 0.0; });
    const double id_ppm = id_model.fitted
      ? ODIA::MassCalibration::systematicResidualPpm(
          test, [&](double mz) { return id_model.ppmAt(mz); })
      : std::numeric_limits<double>::quiet_NaN();
    const double id_const = id_model.fitted
      ? ODIA::MassCalibration::systematicResidualPpm(
          test, [&](double) { return id_model.intercept_ppm; })
      : std::numeric_limits<double>::quiet_NaN();
    const double probe_ppm = (mass_model_known_ && mass_model_.fitted)
      ? ODIA::MassCalibration::systematicResidualPpm(
          test, [&](double mz) { return mass_model_.ppmAt(mz); })
      : std::numeric_limits<double>::quiet_NaN();

    {
      // The two counts that say whether the fit saw what it was meant to see.
      // The first attempt refused with a control peakedness of 1e9, which is
      // only reachable when a control population exists at all -- so it is
      // reported rather than inferred.
      std::ostringstream d;
      d << "mass anchors: fit saw " << id_model.residuals << " target and "
        << id_model.decoy_residuals << " control residuals";
      writeLogInfo_(d.str());
    }

    std::ostringstream r;
    r.setf(std::ios::fixed); r.precision(3);
    r << "mass anchors, HELD-OUT systematic residual (weighted RMS of per-m/z-bin "
         "modes, ppm -- lower is better):"
      << "\n  no correction        " << none_ppm
      << "\n  ID-fitted constant   " << id_const
      << "\n  ID-fitted " << id_model.form << "  " << id_ppm
      << "\n  probe (in force)     " << probe_ppm
      << "\n  ID model: " << (id_model.fitted ? id_model.form : std::string("REFUSED"))
      << ", " << id_model.reason;
    writeLogInfo_(r.str());
  }

  ODIA::PeakGroupScorer::Options scoringOptions_()
  {
    ODIA::PeakGroupScorer::Options options;
    options.classifier = getStringOption_("classifier");
    options.collect_mass_anchors = getStringOption_("mass_anchors") != "off";
    options.max_mass_anchors = static_cast<std::size_t>(
      std::max(0, getIntOption_("max_mass_anchors")));
    options.min_library_corr = getDoubleOption_("min_library_corr");
    options.coelution_picking = !getFlag_("amplitude_picking");
    const std::string picker = getStringOption_("picker");
    options.union_picking = picker == "union" || picker == "union_openswath";
    options.openswath_picking = picker == "openswath" || picker == "union_openswath";
    options.openswath_sn = getDoubleOption_("openswath_sn");
    options.openswath_gauss = getFlag_("openswath_gauss");
    options.openswath_peak_width = getDoubleOption_("openswath_peak_width");
    options.min_corr_score = getDoubleOption_("min_corr_score");
    options.classifier_model_out = getStringOption_("classifier_model_out");
    options.classifier_model_in = getStringOption_("classifier_model_in");
    options.max_corr_diff = getDoubleOption_("max_corr_diff");
    options.max_candidates = static_cast<std::size_t>(
      std::max(1, getIntOption_("max_candidates")));
    options.match_decoy_candidate_counts = !getFlag_("no_match_decoy_n");
    options.train_fdr_initial = getDoubleOption_("train_fdr_initial");
    options.train_fdr = getDoubleOption_("train_fdr");
    options.classifier_iterations = getIntOption_("classifier_iterations");
    options.use_pi0 = getFlag_("use_pi0");
    options.min_fragments_at_apex = static_cast<std::size_t>(
      std::max(1, getIntOption_("min_fragments_at_apex")));
    options.apex_evidence = getDoubleOption_("apex_evidence");
    options.threads = static_cast<unsigned>(std::max(1, getIntOption_("threads")));
    // Built once per run and owned by the tool; null until then, and null
    // forever on a run with no MS1, in which case MS1_COELUTION is NaN for every
    // row and the constant-column guard drops it.
    options.ms1 = ms1_traces_.empty() ? nullptr : &ms1_traces_;
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
  /// Guards the one-shot conversion of an externally supplied iRT map.
  bool external_irt_ = false;
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
    // The refine loop refits the same matrix, so it must withhold the same
    // columns the search did or it would refit against a different feature set.
    options.disabled_sub_scores = ablatedSubScores_();

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
      const auto trafo = ODIA::Calibration::fit(anchors, &p95,
                                                getDoubleOption_("rt_loess_span"),
                                                getStringOption_("rt_interpolation"));
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
    // Only now does the fragment mass calibration's verdict exist -- the probe
    // runs during extraction setup, after the Sink was built. Deciding earlier
    // ran pass 1 with features pass 2 rejects, and pass 1 supplies the anchors.
    sink.disableSubScores(ablatedSubScores_());
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
    auto options = scoringOptions_();
    options.disabled_sub_scores = ablatedSubScores_();
    // Without this RT_DELTA is NaN and the constant-column guard drops
    // var_rt_delta, so every -out_chrom path silently discarded the one feature
    // the retention-time map exists to enable. extractAndScore_ has always set
    // it; this path never did.
    options.library_rt_is_run_seconds = scoring_rt_is_run_seconds_;

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
           "\tDScore\tQValue\tPEP\tMass.Ppm\tMass.Ppm.N";
    for (const auto& n : ODIA::PeakGroupScorer::subScoreNames()) { out << '\t' << n; }
    out << '\n';

    const auto& p = library.precursors();
    for (const auto& g : scored.groups)
    {
      const auto seq = library.strings().get(p.modified_sequence[g.precursor]);
      out << seq << static_cast<int>(p.charge[g.precursor]) << '\t'
          << static_cast<int>(g.decoy) << '\t' << g.apex_rt << '\t' << g.left_rt
          << '\t' << g.right_rt << '\t' << g.apex_intensity << '\t' << g.dscore
          << '\t' << g.qvalue << '\t' << g.pep
          << '\t' << g.mass_ppm << '\t' << g.mass_ppm_n;
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
    if (stop_after != "library" && stop_after != "extract" && stop_after != "calib" &&
        stop_after != "score")
    {
      writeLogError_("Implemented stages are 'library', 'extract', 'calib' and 'score'.");
      return ILLEGAL_PARAMETERS;
    }
    if ((stop_after == "extract" || stop_after == "calib" || stop_after == "score") &&
        in_run.empty())
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
        // A model fine-tuned on one run must not build a library for another.
        //
        // This is not hypothetical. Library v5 was generated with a model from
        // rtfinetune/integrated/, tuned on 500 peptides for 40 epochs against
        // one run's DIA-NN results, and it cost 2,027 confident precursors
        // against v4 (35,556 vs 37,583) -- v4's model being the stock one. The
        // tuned model shipped with its own warning in rt_provenance.json saying
        // exactly this, and the generation script passed it anyway.
        //
        // So: if a provenance file sits beside the model and carries a warning,
        // say it loudly. Fine-tuning belongs INSIDE a run's calibration loop,
        // where the model is built from that run and discarded with it -- never
        // persisted and reused, which is what turns a per-run refinement into a
        // cross-run bias.
        if (!rt_model.empty())
        {
          const auto prov = std::filesystem::path(rt_model).parent_path() /
                            "rt_provenance.json";
          if (std::filesystem::exists(prov))
          {
            std::ifstream in(prov);
            const std::string text((std::istreambuf_iterator<char>(in)),
                                   std::istreambuf_iterator<char>());
            if (text.find("\"warning\"") != std::string::npos)
            {
              writeLogWarn_(
                "-rt_model " + rt_model + " has an rt_provenance.json carrying a "
                "warning, which means it was fine-tuned on a specific run. Using it "
                "to build a library for a DIFFERENT run cost 2,027 confident "
                "precursors when it was last done (v5 35,556 against v4's 37,583). "
                "Fine-tuning belongs inside a run's own calibration loop. Read " +
                prov.string() + " before trusting this library.");
            }
          }
        }

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

    // Re-predict the supplied library's iRT with a (fine-tuned) model.
    //
    // BEFORE the sort and before any calibration, because everything
    // downstream -- the anchors, the map, the pass-2 windows -- reads irt, and
    // changing it underneath a fitted map would compose two unrelated axes.
    if (getFlag_("repredict_irt"))
    {
      const std::string rt_model = getStringOption_("rt_model");
      if (rt_model.empty())
      {
        writeLogError_("-repredict_irt needs -rt_model; there is nothing to predict with.");
        return ILLEGAL_PARAMETERS;
      }
      try
      {
        const auto t_rt = std::chrono::steady_clock::now();
        const std::size_t missing = ODIA::LibraryGenerator::predictRetentionTimes(
          library, rt_model, true,
          static_cast<unsigned>(std::max(0, getIntOption_("threads"))));
        std::ostringstream os;
        os << "re-predicted iRT for the supplied library with " << rt_model << " in "
           << std::chrono::duration<double>(std::chrono::steady_clock::now() - t_rt).count()
           << " s";
        if (missing)
        {
          // NaN, not a made-up number: predictRetentionTimes leaves them NaN and
          // the extractor skips a precursor with no retention time rather than
          // extracting from an invented one.
          os << "; " << missing << " precursors could not be predicted and keep a NaN iRT";
        }
        writeLogInfo_(os.str());
      }
      catch (const std::exception& e)
      {
        writeLogError_(std::string("-repredict_irt failed: ") + e.what());
        return INTERNAL_ERROR;
      }
    }

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
    // "calib" enters the same workflow and returns from inside it once the
    // retention-time map is fitted. Listing it here is not optional: an
    // unrecognised stage falls through every branch and the tool exits 0 having
    // built the library and nothing else -- which is what the first three
    // -stop_after calib runs did, in 0.27 s, reporting success.
    if (stop_after == "score" || stop_after == "calib")
    {
      const auto rc = runScoreWorkflow_(library, in_run, out_chrom,
                                        getStringOption_("out"));
      if (rc != EXECUTION_OK) { return rc; }
    }

    if (!out_lib.empty())
    {
      // SAY WHAT THE RT COLUMN MEANS. After a search, `irt` no longer holds the
      // library's normalised retention time: the map rewrote it to RUN SECONDS
      // for this run, and the refiner may have rewritten it again. A DIA-NN TSV's
      // RT column is normally an iRT, so exporting this without saying so hands
      // out a file that looks transferable and is not.
      //
      // It is exactly what you want for a SERIES on one gradient -- the times
      // are already where the peptides elute -- and exactly wrong anywhere else.
      if (rt_map_fitted_)
      {
        writeLogWarn_(
          "-out_lib is being written AFTER the retention-time map was fitted, so its RT "
          "column holds RUN SECONDS for " + in_run + " -- not a normalised iRT. That is "
          "what makes it useful for other runs on the SAME gradient and wrong for any "
          "other. Run without -in to export the library's own retention times.");
      }
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
  /// The run's MS1 precursor traces, built once before the first scoring pass.
  ODIA::Ms1Traces ms1_traces_;

  bool mass_model_known_ = false;

  /// The fragment window sized from pass 1's identifications, and the
  /// half-width pass 1 was itself extracted through. See `MassWidth`.
  ODIA::MassWidth::Estimate mass_width_;
  double extracted_ppm_ = 0.0;

  /// The refined iRT axis, when per-run refinement was accepted. Empty
  /// otherwise, and the map is then applied to the original values.
  std::vector<float> refined_irt_;

  /// The measured half-width if it may be used, else `fallback`.
  ///
  /// Everything that could make the measurement inadmissible is checked in one
  /// place: the mode has to be `apply`, the estimate has to exist, it must not
  /// be censored by the window it was measured through, and it is floored --
  /// pass 1 identifies its cleanest precursors first, so their scatter is a
  /// lower bound on the run's, not an estimate of it.
  double measuredWidthOr_(double fallback)
  {
    if (getStringOption_("mass_width_from_ids") != "apply") { return fallback; }
    if (!mass_width_.valid || mass_width_.censored) { return fallback; }
    if (!std::isfinite(mass_width_.width_ppm) || !(mass_width_.width_ppm > 0.0))
    { return fallback; }
    return std::max(mass_width_.width_ppm, getDoubleOption_("mass_width_min_ppm"));
  }

  /// Whether the run's retention-time map has been fitted and written into the
  /// library's `irt`. Until it has, the mass probe has no way to know where a
  /// precursor should elute.
  bool rt_map_fitted_ = false;

  /// Whether the cached mass model was measured WITH that map. A model fitted
  /// before the map is not merely older, it was measured under a different and
  /// much worse condition, so it has to be discarded rather than kept.
  bool mass_model_used_rt_map_ = false;

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
      imc.max_im_slope = std::max(0.0, getDoubleOption_("max_im_slope"));
      imc.min_anchors_pooled_slope = static_cast<std::size_t>(
        std::max(0, getIntOption_("im_calib_pooled_slope_min")));
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
                             ODIA::ChromatogramExtractor::Options& options,
                             double rt_window_seconds = 0.0)
  {
    applyMassCalibrationImpl_(library, source, options, rt_window_seconds);
    // Remembered so the width measurement knows the window its residuals were
    // observed through, which is the only thing that makes its censoring test
    // possible. Set here rather than in each branch of the implementation so a
    // new early return cannot quietly leave it stale.
    extracted_ppm_ = options.fragment_ppm;
  }

  void applyMassCalibrationImpl_(const ODIA::Library& library, ODIA::SpectrumSource& source,
                                 ODIA::ChromatogramExtractor::Options& options,
                                 double rt_window_seconds)
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
                                            : measuredWidthOr_(options.fragment_ppm_uncalibrated);
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

    // Re-measure once the retention-time map exists.
    //
    // The first measurement happens before pass 1, when there is no map and the
    // library is spread evenly over the run -- the probe's own log says it will
    // "extract from approximately the wrong retention times". A model measured
    // under that condition is not a stale model to be reused, it is a model
    // taken in the worst available conditions, and the second pass can do
    // strictly better. `applyMobilityCalibration_` already defers for exactly
    // this reason; this arm was left behind.
    //
    // OPT-IN, because measuring it better and applying it are different things.
    // Astral, our own library, at 1% FDR:
    //
    //     latched (gate fails both rounds)   4,275
    //     re-measured (gate then passes)     1,895
    //
    // VOID, PENDING RE-MEASUREMENT: the figures below were taken with an
    // out-of-bounds read in the probe's RT gate (predicted_rt indexed by a
    // library index instead of a slot), so the "RT gating" they describe was a
    // pseudo-random cell filter over uninitialised memory. Kept here only so
    // the claim is not silently dropped. Re-run before quoting.
    //
    // The re-measurement is genuinely better BY ITS OWN METRICS -- probing at
    // the fitted retention times took the control from 98 residuals to 149 and
    // it stopped out-peaking the data, so the null the pass-1 gate choked on
    // was itself an artefact of probing everywhere. And applying the result
    // still costs 56% of the identifications, because a passing gate also
    // NARROWS fragment_ppm off the uncalibrated 15. On this run the narrow
    // window is the more expensive error even when it is correctly centred,
    // which is the opposite of what the extractor's fallback assumes.
    //
    // So: the measurement is fixed and the application is not. Separate the
    // width from the offset before turning this on by default.
    const bool remeasure_with_map =
      getFlag_("mass_calibration_remeasure") &&
      rt_map_fitted_ && !mass_model_used_rt_map_ && rt_window_seconds > 0.0;

    ODIA::MassCalibration::Diagnostics diagnostics;
    if (!mass_model_known_ || remeasure_with_map)
    {
      ODIA::MassCalibration::Options mzc;
      mzc.search_ppm = search;
      mzc.max_precursors = static_cast<std::size_t>(
        std::max(1, getIntOption_("mz_calib_precursors")));
      mzc.cycles = static_cast<std::size_t>(std::max(1, getIntOption_("mz_calib_cycles")));
      mzc.use_ion_mobility = !getFlag_("no_ion_mobility");
      mzc.sigma_multiple = getDoubleOption_("mass_sigma_multiple");
      if (remeasure_with_map)
      {
        // The map has already been written INTO the library's irt, so those
        // values are run seconds and the map the probe needs is the identity.
        mzc.irt_slope = 1.0;
        mzc.irt_intercept = 0.0;
        mzc.rt_window_seconds = rt_window_seconds;
      }
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
      mass_model_used_rt_map_ = remeasure_with_map;
      if (remeasure_with_map)
      {
        writeLogInfo_("fragment mass calibration: RE-MEASURED against the fitted "
                      "retention-time map (+/- " + std::to_string(rt_window_seconds) +
                      " s); the pass-1 model was taken before any map existed");
      }
    }
    writeLogInfo_(ODIA::MassCalibration::report(mass_model_, &diagnostics));

    if (!mass_model_.fitted)
    {
      // The measured centre is applied with the measured width or not at all.
      // A width is the scatter ABOUT a centre, so narrowing to it while leaving
      // the axis at zero is the exact error the header's presence table warns
      // about -- it would keep the tail of the distribution rather than its peak.
      const double from_ids = measuredWidthOr_(-1.0);
      if (from_ids > 0.0)
      {
        options.fragment_ppm_offset = mass_width_.centre_ppm;
        options.fragment_ppm = configured > 0.0 ? std::min(configured, from_ids) : from_ids;
        writeLogInfo_("The mass calibration gate FAILED, but pass 1's own identifications "
                      "measured the fragment error directly; extracting at +/-" +
                      std::to_string(options.fragment_ppm) + " ppm centred on " +
                      std::to_string(options.fragment_ppm_offset) + " ppm.");
        return;
      }
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
    if (getFlag_("mass_calibration_offset_only"))
    {
      // Centre without narrowing, to separate the two things a passing gate
      // does. On Astral a passing gate costs 4,275 -> 1,895 identifications and
      // we do not know whether the offset is wrong or merely the width, because
      // the gate applies both at once and -fragment_ppm only sets the BASELINE
      // that the model is then allowed to narrow below.
      options.fragment_ppm = baseline;
    }
    else if (mass_model_.window_ppm > 0.0 && mass_model_.window_ppm < baseline)
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
