// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

/// Discard precursors the run cannot support, BEFORE the pass that costs
/// memory proportional to how many survive.
///
/// ## Why this exists again after being declared dead
///
/// `doc/08` designed this filter and `doc/13` listed it under "what is dead".
/// That verdict was correct and is not being overturned by argument -- it was
/// measured on S08 and the gradient was flat: 99.7% of precursors reached the
/// maximum depth, and the two non-maximum statistics ranked WORSE than random
/// (0.5x top-decile enrichment).
///
/// `doc/08` also named the cause, and it is not the tolerance:
///
///   "The binding problem is that `depth` is a MAXIMUM over ~32,210 spectra.
///    An extreme-value statistic over thousands of draws saturates whatever
///    the per-draw probability is."
///
/// The number of DRAWS is the problem. `doc/13` recorded three possible
/// rescues and blocked all three; the one it called the rescue that "would
/// actually work" was a retention-time neighbourhood, dismissed as circular
/// because supplying the RT seed was the filter's own second purpose.
///
/// **That circularity binds only if the filter runs before pass 1.** It does
/// not run there. It runs BETWEEN the passes, where pass 1 has already produced
/// an RT map and a calibrated fragment window, and it inherits both. A +/-N
/// second neighbourhood of a ~1800 s gradient cuts the draws by more than an
/// order of magnitude, which is aimed at the stated binding problem rather than
/// at the per-draw probability.
///
/// So: the refutation stands for the regime it was taken in, and this is a
/// different regime. It is still a CONJECTURE that the gradient returns, and
/// the class reports the enrichment table needed to decide -- on two files,
/// per `doc/13`'s rule -- rather than assuming it.
///
/// ## What makes it safe to select before a target-decoy FDR
///
/// `doc/08`'s three rules, and each is load-bearing:
///
/// 1. **Label symmetry by COUNT, not by threshold.** Retain the same number of
///    targets and decoys. A shared threshold is the trap: targets clear an
///    evidence bar more often, so the surviving decoys would be a biased,
///    weaker sample and would then score lower than a fair null, making every
///    downstream q-value optimistic. Equal counts make the retained decoys the
///    BEST decoys, which biases the FDR conservative -- the safe direction.
/// 2. **No learned model here.** `doc/08` requires a learned filter to beat the
///    fixed rule on held-out data before replacing it, and a model fitted on
///    the same fragment evidence the scorer later uses needs k-fold
///    cross-fitting to avoid inflating the scorer's apparent separation. The
///    fixed rule has neither problem, so it is what is implemented.
/// 3. **Report what was dropped.** A filter that silently discards is
///    indistinguishable from a search that found nothing.

#pragma once

#include <odia/Library.h>
#include <odia/SpectrumSource.h>

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace ODIA
{
  class PrecursorPrefilter
  {
  public:
    struct Options
    {
      /// How many of the highest-intensity library fragments define a
      /// precursor's signature. `doc/08`'s measurement used 6.
      std::size_t top_n = 6;

      /// Fragment tolerance, ppm, and the centre it is applied about --
      /// the pass-1 calibration, not a default. `doc/08` was measured at 15 ppm
      /// about zero because no calibrated window existed then; one does now.
      double ppm = 15.0;
      double ppm_centre = 0.0;

      /// Mobility half-window. On diaPASEF this is MANDATORY and not an
      /// optimisation: a frame is ~600-810 concatenated TIMS scans, so without
      /// it "one spectrum" is a whole frame and every precursor saturates.
      /// Measured: 50,000 of 50,000 at maximum depth, 0.9x enrichment.
      double im_window = 0.0;

      /// Retention-time half-window, seconds, about where the pass-1 map says
      /// the precursor elutes. THIS is the term the refutation turned on.
      /// <= 0 disables the restriction, which reproduces the regime `doc/08`
      /// was refuted in.
      double rt_half_window = 0.0;
      double irt_slope = 0.0;      ///< rt = irt_slope * irt + irt_intercept
      double irt_intercept = 0.0;

      /// Fraction of each label class to retain. 1.0 keeps everything and makes
      /// this a pure measurement pass. Applied as a COUNT per class -- see the
      /// label-symmetry rule above.
      double keep_fraction = 1.0;

      /// A spectrum counts towards contiguity when it reaches this depth.
      /// Defaults to top_n - 1 when left at 0.
      std::size_t contiguity_depth = 0;

      /// Never retain fewer than this many per class, so a filter that fails on
      /// an unusual run degrades to "kept too much" rather than to an empty
      /// search.
      std::size_t min_keep_per_class = 1000;
    };

    struct Stats
    {
      std::size_t spectra_swept = 0;
      std::size_t targets_in = 0, decoys_in = 0;
      std::size_t targets_kept = 0, decoys_kept = 0;
      std::uint8_t depth_threshold = 0;   ///< the depth the retained set reached
      /// Depth histogram over ALL precursors, index = depth. Reported for the
      /// discarded set too, because that is where a filter's damage shows.
      std::vector<std::size_t> depth_hist_target, depth_hist_decoy;
      /// Same, over contiguity, capped at the last bin. This is the histogram
      /// that says whether the orthogonal statistic separates where depth did
      /// not.
      std::vector<std::size_t> contig_hist_target, contig_hist_decoy;
      double seconds = 0.0;
      std::string note;
    };

    /// Per precursor, in full-library indexing.
    struct Evidence
    {
      std::uint8_t depth = 0;
      std::uint32_t total_matches = 0;
      float best_rt = -1.0f;

      /// The longest run of CONSECUTIVE cycles of this precursor's isolation
      /// window in which it reached `contiguity_depth`.
      ///
      /// doc/08 listed this and it was never built: "whether the qualifying
      /// spectra are consecutive cycles. A peptide elutes over a peak; random
      /// coincidences scatter. Cheap, and orthogonal to everything else here."
      ///
      /// Orthogonal is the operative word. `depth` asks how good the best
      /// single moment was, and a chance coincidence can win that outright --
      /// which is why it saturated. A chance coincidence cannot easily repeat
      /// in the NEXT cycle and the one after, because the interfering ions
      /// that produced it are not eluting on this precursor's peak.
      std::uint16_t contiguity = 0;

      /// The centre of that run: an elution-peak midpoint rather than the
      /// single spectrum that happened to match best, and therefore the
      /// retention time worth calibrating against.
      float contiguous_rt = -1.0f;
    };

    /// Sweep the run and score every precursor. Does not select.
    /// Restrict the work to these precursors (1 = consider), or null for all.
    ///
    /// NOT an optimisation detail -- it is the difference between seeding from
    /// the CiRT standards and doing a whole-library blind search and then
    /// throwing 99.99% of it away. With -rt_seed cirt only 708 of 9,617,705
    /// precursors are ever read out of the result, but every one of the
    /// 9,617,705 was measured: the per-window index is rebuilt by scanning all
    /// of them and sorting ~115M fragment entries, once per isolation window.
    /// Measured before this existed: 4h11 and still running, single-threaded.
    ///
    /// Indices are NOT renumbered -- the returned Evidence vector stays the
    /// size of the library so every downstream index remains valid; entries
    /// outside the subset are simply left at their default.
    static std::vector<Evidence> measure(const Library& library, SpectrumSource& source,
                                         const Options& options, Stats& stats,
                                         const std::vector<char>* consider);

    static std::vector<Evidence> measure(const Library& library, SpectrumSource& source,
                                         const Options& options, Stats& stats);

    /// Turn evidence into a keep-mask, applying label symmetry by count.
    /// Separate from `measure` so the enrichment can be reported for a run that
    /// does not actually discard anything.
    static std::vector<char> select(const Library& library,
                                    const std::vector<Evidence>& evidence,
                                    const Options& options, Stats& stats);
  };

} // namespace ODIA
