// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/Library.h>
#include <odia/MobilityCalibration.h>
#include <odia/SpectrumSource.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace ODIA
{

  /// Extracted-ion chromatograms for a library's transitions, ALL OF THEM AT
  /// ONCE.
  ///
  /// Stored as one flat point array with a CSR index, for the same reason the
  /// library is: a vector-per-transition would be millions of small
  /// allocations, which is exactly the arena fragmentation D3 exists to avoid.
  ///
  /// This is no longer what extraction produces. It is what a caller asks for
  /// when it needs every chromatogram simultaneously -- `-out_chrom` and the
  /// diagnostics built on it -- and it is bounded by the LIBRARY, which is the
  /// thing `ChromatogramSink` exists to escape. A full human library at IH1's
  /// geometry is 274 GiB here. Everything else takes a
  /// `PrecursorChromatogram` as the pass finishes it.
  ///
  /// Two space decisions, both taken from the reference implementation:
  ///
  /// **Retention time is an index, not a number.** A transition's points sit on
  /// its window's acquisition grid, so storing a float per point stored the
  /// same ~1,000 timestamps once per transition. The axis is held once per
  /// window and a transition records where on it its run starts. At 157 M
  /// points that is 630 MB of duplicated timestamps not stored. The deeper
  /// reason is correctness rather than size: recalibrating a run rewrites one
  /// array of ~1,000 doubles and every point follows, instead of rewriting
  /// every structure that cached a time and silently missing one.
  ///
  /// That representation says a transition's points are ONE contiguous run on
  /// ONE window's axis, and the extractor now guarantees it by extracting each
  /// precursor from a single window (see `precursors_in_several_windows`). The
  /// invariant it rests on is
  ///
  ///     axis_begin[t] + count[t] <= axes[axis_of[t]].size()
  ///
  /// and the extractor checks it rather than assuming it: violating it is a
  /// silent heap over-read through retentionTime(), which is what happened when
  /// two windows' cycles were concatenated into one transition's run.
  ///
  /// Intensity stays float32 for now; one-byte log quantisation is in the
  /// backlog. It is a further 4x and costs nothing real, because the scale is
  /// per transition rather than global -- 255 log steps over one transition's
  /// own two-decade range is ~2.7% per step.
  ///
  /// The axis is behind an accessor, so callers index by (transition, position)
  /// and never see the representation.
  struct Chromatograms
  {
    /// One retention-time axis per isolation window, seconds, ascending.
    std::vector<std::vector<float>> axes;

    /// Where a precursor's points sit: which axis, where on it point 0 is, and
    /// how many cycles. **Indexed by precursor, not by transition.**
    ///
    /// These were per transition until it was measured what that costs. The
    /// extractor already guarantees a precursor is extracted from ONE window
    /// over ONE cycle range -- that is what `precursors_in_several_windows`
    /// enforces -- so all twelve of a precursor's transitions carried three
    /// identical values. `PeakGroupScorer` gave the game away by only ever
    /// reading index `transition_begin` and never the other eleven.
    ///
    /// The redundancy is not small at the size that matters. A 4.26 M-precursor
    /// library at twelve transitions is 51.1 M transitions, and per transition
    /// `begin`(8 B) + `count`(4 B) + `axis_of`(4 B) + `axis_begin`(4 B) is
    /// **1.02 GB of index for 4.26 M distinct facts**. Per precursor it is
    /// ~85 MB. The axes themselves are ~129 KB and were never the point.
    std::vector<std::uint32_t> precursor_axis;
    std::vector<std::uint32_t> precursor_axis_begin;
    std::vector<std::uint32_t> precursor_cycles;

    /// Which transition each precursor's run starts at, so a transition can be
    /// resolved back to its precursor without a per-transition map.
    std::vector<std::uint32_t> precursor_transition_begin;

    /// `begin` and `count` are still PER TRANSITION, and that is the other
    /// 614 MB of the 1.02 GB.
    ///
    /// They are not collapsed here because the stride is not simply
    /// `(t - transition_begin) * cycles`: a transition with no representable
    /// product m/z gets count 0 and is **packed out** of the flat array, so a
    /// later sibling's offset depends on how many of its predecessors were
    /// real. Deriving it needs a validity bitmap (1 bit per transition, 6.4 MB
    /// at 51.1 M) and a rank, plus `ChromatogramLayout` carrying per-precursor
    /// shape instead of the per-transition `counts` it passes today -- which is
    /// itself 204 MB and has the same problem. That is a bigger change than
    /// this one and belongs on its own.
    ///
    /// Getting it wrong is silent: multiplying by position rather than by rank
    /// reads every transition after the first absent one from its neighbour's
    /// points, with all offsets in range and all totals matching.
    std::vector<std::uint64_t> begin;
    std::vector<std::uint32_t> count;

    /// Every point of every transition, one flat array.
    ///
    /// `precursor_begin` is 64-bit deliberately. A 32-bit offset ran out at
    /// 2^32 points, which is 266,664 precursors at IH1's 12 transitions x 1,342
    /// cycles -- and phase 1's own human library is 4,255,113 precursors, so
    /// ODIA could not extract against the library it had just generated.
    /// Measured by bisection at 260,000 passing and 270,000 throwing in 5.1 s.
    ///
    /// The widening removed the index limit and nothing else. The limit that
    /// remains is memory: this array is sized transitions x cycles before a
    /// peak is seen, so a full human library at IH1's geometry wants
    /// 4.26 M x 12 x 1,342 x 4 B = 274 GiB -- and 51.3% of those points are
    /// non-zero, so sparse storage is not a way out either.
    ///
    /// That limit belongs to THIS TYPE and not to extraction any more. The
    /// extractor hands a precursor over when the pass leaves its retention-time
    /// window and frees it; `ChromatogramCollector` is what puts them all back
    /// into one array, and a caller that does not need them all never asks for
    /// one. A run that does ask is refused with the number rather than left to
    /// a bad_alloc.
    std::vector<float> intensity;

    /// Which precursor owns a transition. Needed only where a caller has a
    /// transition and no precursor; the extractor and the scorer both walk
    /// precursors and never call this.
    ///
    /// Binary search over `precursor_transition_begin` rather than a stored
    /// 4-byte-per-transition map, because storing that map would hand back
    /// 204 MB of the 1.02 GB this change was made to remove.
    std::uint32_t precursorOf(std::uint32_t transition) const
    {
      const auto it = std::upper_bound(precursor_transition_begin.begin(),
                                       precursor_transition_begin.end(), transition);
      return static_cast<std::uint32_t>(it - precursor_transition_begin.begin() - 1);
    }

    /// Which axis @p transition sits on, and where its first point is.
    ///
    /// Resolved through the precursor, because that is where the fact lives
    /// now. Callers that already know the precursor should read
    /// `precursor_axis` / `precursor_axis_begin` directly and skip the search.
    std::uint32_t axisOf(std::uint32_t transition) const
    {
      return precursor_axis[precursorOf(transition)];
    }
    std::uint32_t axisBegin(std::uint32_t transition) const
    {
      return precursor_axis_begin[precursorOf(transition)];
    }

    /// Precursors that no isolation window covered. They cannot be extracted
    /// and are reported rather than silently absent from the output.
    std::size_t precursors_without_window = 0;

    /// Precursors that MORE than one window covered -- 1.1% of them on 12_80
    /// and on IH1, whose schemes overlap adjacent windows by 1.0 Th.
    ///
    /// Each such window is a separate measurement of the same ion at
    /// interleaved times, and only one of them is extracted: the window whose
    /// centre the precursor is nearest, which is where it is furthest from the
    /// edge and best transmitted. The other is discarded, so it is counted
    /// here rather than silently dropped.
    ///
    /// The alternative -- concatenating the windows' cycles into one point run
    /// -- is what this field replaces. It broke the axis (a run cannot start at
    /// two places on two axes at once, so the last window won both fields and
    /// every point of the first was read from the wrong axis, past its end) and
    /// it also handed every consumer a "chromatogram" that runs backwards in
    /// time at the seam and elutes the same peptide twice. PeakGroupScorer
    /// reads these points as one time series -- median/MAD background,
    /// cross-correlation with a lag, apex/left/right retention time -- so
    /// repairing only the axis lookup would have left the scorer reading a
    /// trace that no acquisition ever produced.
    std::size_t precursors_in_several_windows = 0;

    std::size_t points() const { return intensity.size(); }
    std::size_t footprintBytes() const;

    /// Retention time of point @p j of @p transition, in seconds.
    float retentionTime(std::uint32_t transition, std::uint32_t j) const
    {
      return retentionTimeOfPrecursor(precursorOf(transition), j);
    }

    /// Retention time of point @p j of @p precursor, in seconds. Prefer this:
    /// the RT is a property of the precursor's cycle range, not of which of its
    /// transitions is being read, and this form skips the `precursorOf` search.
    float retentionTimeOfPrecursor(std::uint32_t precursor, std::uint32_t j) const
    {
      return axes[precursor_axis[precursor]][precursor_axis_begin[precursor] + j];
    }

  };

  /// One precursor's finished chromatograms, as a VIEW onto whoever owns the
  /// points -- the streaming extractor's block pool, or a whole `Chromatograms`.
  ///
  /// This is the unit the extractor produces and every consumer consumes. It
  /// exists because a precursor's chromatogram is complete the moment the
  /// forward pass leaves its retention-time window, and holding it after that
  /// is what made memory scale with the library instead of with the run: see
  /// `ChromatogramSink`.
  ///
  /// A precursor's transitions all sit on ONE window's axis over ONE cycle
  /// range, so `rt` and `axis_begin` are per precursor and not per transition.
  /// The per-transition (offset, count) pair is kept rather than a single
  /// stride because a transition with no representable product m/z gets no
  /// points at all -- count 0 -- and the flat `Chromatograms` it must also be
  /// able to describe packs those out.
  struct PrecursorChromatogram
  {
    std::uint32_t precursor = 0;
    std::uint32_t transition_begin = 0;   ///< first transition, library index
    std::uint32_t transition_count = 0;
    std::uint32_t axis = 0;               ///< which isolation window's axis
    std::uint32_t axis_begin = 0;         ///< first cycle on it
    std::uint32_t cycles = 0;             ///< points per live transition

    const float* rt = nullptr;            ///< `cycles` retention times, seconds
    const float* points = nullptr;        ///< base of the intensity storage
    /// Intensity-weighted m/z deviation per point: sum(intensity*ppm) and
    /// sum(intensity). Divide to get the weighted mean; a zero denominator
    /// means no peak matched that cell. Null when collect_mass_residuals was
    /// off. Same indexing as `points`.
    const float* ppm_num = nullptr;
    const float* ppm_den = nullptr;
    /// Intensity-weighted observed 1/K0 per point: sum(intensity * im) and
    /// sum(intensity). Divide to get the weighted mean; a zero denominator
    /// means no peak matched that cell. Null when collect_im_residuals was off.
    ///
    /// This is what makes IM_DELTA computable at all -- nothing else in the
    /// pipeline measures the OBSERVED mobility of a precursor. `MobilityAnchor`
    /// carries only {precursor, rt, decoy}, so the feature has been all-NaN on
    /// every run since it was added.
    ///
    /// And it is per-FRAGMENT, which is the point: a fragment whose observed
    /// 1/K0 disagrees with its precursor's is not that precursor's fragment.
    /// On diaPASEF a frame merges TIMS scans, so co-eluting interference that
    /// is inseparable in m/z and RT is often cleanly separable here.
    const float* im_num = nullptr;
    const float* im_den = nullptr;
    const std::uint64_t* offset = nullptr;///< transition_count offsets into it
    const std::uint32_t* count = nullptr; ///< transition_count point counts

    /// Points of transition @p k of this precursor, `pointCount(k)` of them.
    const float* trace(std::uint32_t k) const { return points + offset[k]; }
    std::uint32_t pointCount(std::uint32_t k) const { return count ? count[k] : 0; }
    float retentionTime(std::uint32_t j) const { return rt[j]; }

    /// False for a precursor nothing extracted -- no isolation window covered
    /// it, or it was predicted to elute outside the run. Such a precursor is
    /// still handed to the sink, empty, so that a consumer counting
    /// "precursors that yielded nothing" sees the same number it saw when every
    /// precursor had a row in one flat array.
    bool extracted() const { return cycles != 0; }
  };

  /// The shape of what a pass will produce, known before a peak is read.
  struct ChromatogramLayout
  {
    std::size_t transitions = 0;               ///< entries the CSR would have
    std::size_t precursors = 0;                ///< precursors it covers
    std::uint64_t points = 0;                  ///< points it would hold
    const std::vector<std::vector<float>>* axes = nullptr;  ///< one per window

    /// Points per transition, in transition order -- present only for a sink
    /// that asked (`ChromatogramSink::needsLayoutCounts`).
    ///
    /// It is four bytes per transition, 204 MB at 51 M transitions, and the
    /// only reason it is optional is that a streaming sink must have NO array
    /// proportional to the library. A sink that lays out one flat CSR needs it,
    /// because the offsets are that CSR's prefix sum in transition order and a
    /// precursor arriving in retention-time order cannot know its own offset.
    const std::vector<std::uint32_t>* counts = nullptr;
  };

  /// Where finished chromatograms go.
  ///
  /// The extractor no longer materialises every chromatogram before anything
  /// reads one. It hands each precursor over as the forward pass leaves that
  /// precursor's retention-time window, and then frees it -- so memory scales
  /// with how many windows overlap at one instant, not with the library. Two
  /// implementations, and the choice between them IS the memory decision:
  ///
  ///  * `ChromatogramCollector` keeps everything, which is what `-out_chrom`
  ///    and every diagnostic built on it need. It is bounded by a precursor
  ///    count, not by the run.
  ///  * `PeakGroupScorer::Sink` scores each precursor and drops it. That is
  ///    the default at scale and has no term proportional to the library at
  ///    all.
  ///
  /// The storage a `PrecursorChromatogram` points at is valid ONLY for the
  /// duration of `accept`. A sink that wants it afterwards must copy it, which
  /// is exactly what the collector does and what the scorer does not.
  class Ms1Traces;

  class ChromatogramSink
  {
  public:
    virtual ~ChromatogramSink() = default;

    /// Called once per extraction, before any trace, with what the whole run
    /// would produce -- including the parts that will never exist at once.
    virtual void begin(const ChromatogramLayout&) {}

    /// The run's MS1 traces, once they exist, or nullptr when there are none.
    ///
    /// This hook exists because the traces are built INSIDE the extraction, by
    /// which point the caller has already constructed its sink -- so a sink that
    /// captured the pointer at construction captured a null one and kept it. That
    /// is exactly what made MS1_COELUTION report `0 finite, 11,877 null ptr` on
    /// every candidate while a 148 GiB matrix sat fully built beside it. Default
    /// is a no-op: most sinks do not care.
    virtual void ms1Available(const Ms1Traces*) {}

    virtual void accept(const PrecursorChromatogram&) = 0;

    /// Batched delivery, for a sink that can consume several released
    /// precursors at once (the -parallel_sink scorer). Zero, the default,
    /// means "call accept() once per precursor", which is the historical path
    /// and the only one the extractor takes for such a sink.
    ///
    /// Non-zero: the extractor may hand over up to this many precursors in
    /// one `acceptBatch` call, in EXACTLY the order it would have called
    /// `accept`, and keeps their storage alive until the call returns. Each
    /// batch is a barrier: nothing later is released until it has returned.
    virtual std::size_t batchCapacity() const { return 0; }
    virtual void acceptBatch(const PrecursorChromatogram* traces, std::size_t n)
    {
      for (std::size_t i = 0; i < n; ++i) { accept(traces[i]); }
    }

    /// Whether `ChromatogramLayout::counts` has to be filled in. See it.
    virtual bool needsLayoutCounts() const { return false; }
  };

  /// The sink that keeps every point: the behaviour this class had before
  /// there was a choice.
  ///
  /// Memory is the whole library's worth of points, so this is for a bounded
  /// precursor count -- `-out_chrom`, the oracle experiments, the bucket
  /// attribution. It reproduces the flat `Chromatograms` exactly, including
  /// which transitions get zero points, so a file written through it is byte
  /// for byte the file written before the extractor streamed.
  class ChromatogramCollector final : public ChromatogramSink
  {
  public:
    void begin(const ChromatogramLayout& layout) override;
    void accept(const PrecursorChromatogram& trace) override;
    bool needsLayoutCounts() const override { return true; }

    Chromatograms& chromatograms() { return c_; }
    /// Repairs the precursor index before handing the table over. Done HERE
    /// rather than left to the caller because there are two extraction paths --
    /// ChromatogramExtractor::extract() and OpenDIAlyzer's own -out_chrom path,
    /// which builds a collector itself -- and only one of them was calling it.
    /// A caller that forgets gets a silently wrong retention-time column.
    Chromatograms take() { repairPrecursorIndex(); return std::move(c_); }

    /// Make `precursor_transition_begin` monotonic before anyone binary-searches
    /// it. Idempotent; take() calls it.
    void repairPrecursorIndex();

  private:
    Chromatograms c_;
  };

  /// A sink that counts and discards. Useful on its own -- `-stop_after
  /// extract` with no `-out_chrom` is a decode-and-match benchmark -- and it
  /// is the control that says what the sink costs.
  class NullChromatogramSink final : public ChromatogramSink
  {
  public:
    void accept(const PrecursorChromatogram&) override {}
  };

  /// Reads a run once and extracts every requested transition from it.
  ///
  /// One forward pass, not one pass per precursor. The cost of decoding a
  /// spectrum is paid once and shared by every transition that lands in it,
  /// which is what makes the extraction scale with the run rather than with
  /// the library. With the current mzPeak reader a spectrum costs ~294 ms to
  /// decode, so this ordering is the difference between one pass and none.
  ///
  /// That same forward order is what bounds the memory. A precursor's
  /// chromatogram is COMPLETE the moment the pass passes the end of its
  /// retention-time window: nothing later in the run can add to it. So a
  /// chromatogram is allocated when the pass reaches the start of its window,
  /// filled as the pass proceeds, handed to the sink when the pass leaves it,
  /// and freed. See `ChromatogramSink`.
  class ChromatogramExtractor
  {
  public:
    struct Options
    {
      /// Fragment mass tolerance. Relative, because that is how a mass
      /// spectrometer's accuracy behaves; an absolute window would be far too
      /// wide at 200 Th and far too narrow at 1800.
      /// Fragment mass tolerance, as a HALF-width in ppm: a peak matches when
      /// |peak - transition| <= transition * fragment_ppm * 1e-6.
      ///
      /// 10, which is correct ONCE THE WINDOW IS CENTRED. 15 was the
      /// uncalibrated fallback and is still the right answer when it is not:
      /// see `MassCalibration`, which fits the centring per run and gates
      /// itself, and the tool, which asks for 15 when that gate fails.
      ///
      /// This instrument has a systematic fragment mass offset of about
      /// -10 ppm, measured on IH1 as the centroid of the retention-time-
      /// specific excess over a local decoy-cell null: -11.2 ppm on precursors
      /// we recover, -12.6 ppm on those we miss, and present at off-peak times
      /// too, so it is a genuine instrument term and not an artefact of peak
      /// selection. Fraction of true fragments captured:
      ///
      ///     half-width   centred on 0   centred on -10 ppm
      ///        +/- 5        0.22             0.52
      ///        +/-10        0.51             0.81
      ///        +/-15        0.78             0.92
      ///
      /// So narrowing a window that is centred on the wrong place THROWS AWAY
      /// signal: +/-10 about zero keeps half the evidence, and the half it
      /// keeps is the tail rather than the peak. Until `fragment_ppm_offset` is
      /// fitted per run, the safe width is therefore the wider one.
      ///
      /// With the offset fitted, narrow is measurably right and wide measurably
      /// wrong. Target-minus-decoy fragment presence falls monotonically with
      /// tolerance -- 0.083 at 10 ppm, 0.051 at 20, 0.034 at 30, 0.019 at 50 --
      /// and +/-10 centred on -9.8 measured x1.13 overall and x1.24 in the
      /// weakest abundance decile, with a further x1.12 inside a per-precursor
      /// mobility band. So 10 it is, and the calibration is what earns it.
      double fragment_ppm = 10.0;

      /// The width to fall back to when no calibration could be fitted.
      ///
      /// Not a second tolerance: it is the same decision as `fragment_ppm`,
      /// taken with less information. An uncentred +/-10 keeps 0.51 of true
      /// fragments where an uncentred +/-15 keeps 0.78, so a run whose gate
      /// fails is strictly better off wide -- and a run that silently narrowed
      /// anyway would look like a calibration working.
      /// The width used when the mass calibration gate FAILS.
      ///
      /// Raised 15 -> 50 because 15 was the wrong direction. A failed gate means
      /// the run's mass error could NOT be measured, and the safe response to not
      /// knowing is a wide window, not a moderately narrow one. 15 ppm was narrow
      /// enough to lose real fragments while too wide to be a real tolerance.
      ///
      /// Measured on Astral, where the gate fails, sweeping -fragment_ppm:
      ///
      ///     ppm        15      30      50      75     100
      ///     IDs     4,290   4,499   4,969   4,382   3,967
      ///     all-0  26,800  23,203  15,550   7,341   2,436
      ///
      /// The optimum is 50, which is exactly OpenSWATH's `mz_extraction_window`
      /// default. All-zero traces fall monotonically to 100 ppm, so wider windows
      /// keep finding real signal; identifications turn over at 50 because past
      /// that the interference admitted costs more than the signal recovered.
      ///
      /// THIS IS NOT A GLOBAL WIDENING. On IH1, where the gate PASSES, the
      /// calibrated width (~10 ppm) is used and is far better: 1,306 identifications
      /// against 922 with calibration off at 15 ppm, and 50 ppm with calibration
      /// off collapses the run to zero. A diaPASEF frame is a merged stack of TIMS
      /// scans, so a wide m/z window admits vastly more interference there than in
      /// Astral's cleaner spectra. The two instruments want opposite widths, which
      /// is precisely why this value must apply ONLY when the gate could not
      /// measure the run.
      double fragment_ppm_uncalibrated = 50.0;

      /// Systematic fragment mass offset, ppm, added to the theoretical m/z
      /// before matching. 0 means uncalibrated.
      ///
      /// Should be fitted from the run rather than set by hand -- that is what
      /// per-run mass calibration is for, and it is the single measured
      /// difference most likely to account for the recovery gap. -10 is the
      /// measured value for IH1 and is NOT a default, because it is a property
      /// of that instrument and that acquisition. `MassCalibration` fits it;
      /// this is where its answer is applied.
      ///
      /// Note this is a RECALIBRATION of the mass axis and is a different thing
      /// from `fragment_ppm`, which is the width around it. Narrowing a window
      /// centred on the wrong place discards signal -- see the table above.
      double fragment_ppm_offset = 0.0;

      /// m/z dependence of that offset, about `fragment_ppm_ref_mz`. Both zero
      /// means the correction is a constant.
      ///
      /// Two named fields rather than one slope plus a basis flag, because the
      /// UNITS differ -- ppm per e-fold against ppm per 1000 Th -- and a single
      /// number whose meaning depends on a neighbouring enum is the same silent
      /// unit error this header separates `im` from `ccs` to avoid. At most one
      /// is ever non-zero; both are applied additively, so there is no invalid
      /// combination to guard against.
      ///
      /// A TOF's calibration error is characteristically a function of m/z, and
      /// on IH1 it measurably is: the per-m/z-bin modes run from -12.91 ppm at
      /// 238 Th to -7.43 at 1,102 (and to -3.5 at 1,327 Th in a measurement
      /// reaching further up the range). A single constant therefore mis-centres
      /// by ~3 ppm at BOTH ends in opposite directions, against a +/-10 ppm
      /// window -- which costs the lightest and heaviest fragments
      /// preferentially. Correcting it removes 68% of that systematic error:
      /// the bin modes go from 1.71 ppm off a constant to 0.54 ppm off the fit.
      ///
      /// Note this is nearly invisible in the total per-hit scatter (4.66 ppm ->
      /// 4.19), because that is dominated by irreducible per-fragment noise that
      /// no calibration can touch. Judging the correction by its effect on total
      /// scatter understates it structurally, and `MassCalibration` deliberately
      /// does not -- an earlier version did, and refused this correction.
      double fragment_ppm_log_slope = 0.0;         ///< ppm per e-fold in m/z
      double fragment_ppm_slope_per_1000 = 0.0;    ///< ppm per 1000 Th
      double fragment_ppm_ref_mz = 700.0;

      /// Half-width of the ion-mobility window around the PRECURSOR's own
      /// library 1/K0, in 1/K0 units. 0 disables it.
      ///
      /// Distinct from the isolation window's band, which is what
      /// `use_ion_mobility` gates and which only separates co-packed windows.
      /// A frame's band is ~0.40 wide on IH1; a precursor occupies ~0.05 of it.
      /// Filtering only by the band therefore admits the entire same-window
      /// mobility axis -- measured as ~8.5x more mobility than the reference
      /// accepts, and the reason a band-only fix bought 1.12x while the
      /// same-window interference it left behind is what dominates.
      ///
      /// 0.025 half-width was justified as ~2.6 sigma on the measured
      /// library-vs-observed agreement (SD 0.019, and 0.0186 on precursors we
      /// currently miss, so not a selection effect).
      ///
      /// THAT JUSTIFICATION IS OVERSTATED. The 0.019 was computed on residuals
      /// that still contained the mobility TREND -- the CCS->1/K0 conversion is
      /// a pure proportionality and carries a ~10% scale error, so the residual
      /// ran +0.017 at 1/K0 0.7 to -0.030 at 1.3. A single SD over a population
      /// with a trend in it is not a width, it is a trend and a scatter added
      /// together, and the scatter is the only part a symmetric window can
      /// cover. With the trend removed by Curve::im_slope, the residual SD on
      /// held-out anchors is 0.025 -- so this window is ~1 sigma, not 2.6, and
      /// rejects 28.3% of target anchors on IH1 even when the calibration is
      /// correct. At +/-0.05 that falls to 5.6%.
      ///
      /// MEASURED, and the answer is to leave it at 0.025. IH1, lib_targets,
      /// with the mobility slope fitted, identifications at 1% FDR:
      ///
      ///     0.025   1232      <- default
      ///     0.035   1215
      ///     0.050   1167
      ///
      /// Widening is monotonically worse. The 28.3% of anchors outside the
      /// window are real, but they are the interference-prone tail, and
      /// admitting them costs more than the precursors they bring: at 0.050 the
      /// loss is 65 identifications, back to roughly the no-calibration figure.
      /// This window's purpose is excluding the ~8.5x same-window mobility
      /// above, and it is earning that.
      ///
      /// So the residual argument -- "SD 0.025, therefore this window is 1
      /// sigma and too tight" -- was right about the width and wrong about the
      /// remedy. The fix for a trend is to model the trend, which is what
      /// Curve::im_slope does; it is not to widen the gate until the trend
      /// fits through it.
      ///
      /// `MassCalibration::Options::im_window` is 0.010 -- tighter than this,
      /// through which the offset it fits is then applied. That difference is
      /// deliberate and measured: this window is sized to KEEP a peptide's real
      /// fragments, that one to isolate a residual mode a shifted control does
      /// not have, and at 0.025 the calibration's gate fails outright on IH1
      /// (peakedness 4.41 against a control at 4.60). See that field.
      double precursor_im_window = 0.025;

      /// Keep the m/z deviation of every matched peak, for mass recalibration.
      ///
      /// The deviation is computed anyway, to test the match, and has always
      /// been thrown away. Keeping it is the only way to measure the run's real
      /// fragment mass error: a standalone probe asks "is there a peak within X
      /// ppm somewhere near this time", and on a mostly-absent library the
      /// answer is yes by coincidence -- measured 2026-08-08, an RT-shifted
      /// control produced a LARGER apparent offset than the true apex did.
      /// Matches recorded here are constrained by co-elution and, once scored
      /// at q<=0.01, by the whole discriminant, so their contamination is
      /// bounded by the FDR rather than by the search space.
      ///
      /// Costs one extra float plane per live block, so it is off by default.
      bool collect_mass_residuals = false;

      /// Keep the intensity-weighted observed 1/K0 planes. Costs two more float
      /// planes per live chromatogram, the same as the mass residuals. Silently
      /// inert on a run with no ion mobility.
      bool collect_im_residuals = false;

      /// Per-run recalibration of the library's 1/K0, or null for none.
      ///
      /// A POINTER, where the mass calibration's correction is three scalars
      /// copied into this struct. The difference is the shape of the two
      /// models: the mass correction is two coefficients and a reference m/z,
      /// which flatten into options; the mobility correction is a per-charge,
      /// per-fold table plus the anchor list that decides which fold a
      /// precursor is corrected by, and copying that into every Options would be
      /// copying a table per pass to save a lifetime rule.
      ///
      /// The rule is that it must outlive the extraction, which it does by
      /// construction: it is measured once per RUN and cached by the caller
      /// across passes, exactly like the mass model, because re-measuring it in
      /// pass 2 through a window pass 1 moved is a feedback loop.
      const MobilityCalibration::Model* mobility_model = nullptr;

      /// How several peaks inside one transition's box become one number.
      ///
      /// `Sum` integrates; `Max` takes the largest. Max was the original and is
      /// wrong for this data: a frame's peak array is the concatenation of
      /// 600-810 TIMS mobility scans, so peaks inside one m/z tolerance are the
      /// same ion across many mobility steps PLUS every co-isolated interferent
      /// at every other step. Max over that returns the interference envelope,
      /// and -- being an order statistic -- acts as a hard threshold rather
      /// than a graded penalty: a peptide below the envelope contributes
      /// nothing at all. That is the shape of the measured failure, where the
      /// bottom four abundance deciles sat exactly at the decoy null.
      ///
      /// Sum is only correct once the box is small. Summing over the old
      /// oversized box makes the trace worse, not better.
      enum class Aggregate { Sum, Max };
      Aggregate aggregate = Aggregate::Sum;

      /// Extract only the first N precursors of the library, 0 for all.
      std::size_t max_precursors = 0;

      /// Cap on how many precursors may have their chromatograms live at one
      /// instant. 0 lets the data decide.
      ///
      /// The sliding window bounds memory by RETENTION-TIME OVERLAP, which is
      /// a property of the run and the window width rather than of the library
      /// -- and on a wide window that bound is weak. With a 600 s half-width on
      /// a 1,859 s gradient a precursor is live for ~1,200 s, so about 65% of
      /// the library is live at once and the win over holding all of it is
      /// 1.5x. At a calibrated 120 s it is ~13% and the win is 7x.
      ///
      /// This is the second mechanism, for when the first is not enough. Above
      /// the cap the library is split into chunks whose live sets each fit, and
      /// each chunk is a separate pass -- restricted to the spectra its own
      /// precursors need, so the extra decode is the chunks' retention-time
      /// spans and not a whole run per chunk. It is stated in precursors rather
      /// than bytes because that is the number a caller can reason about
      /// against a library size; `Stats::peak_live_points` reports what it
      /// cost.
      ///
      /// A pass is expensive -- decode is ~600 s on IH1 and is 65% of Phase 2
      /// -- so chunking is a fallback, not a default. Which mechanism actually
      /// bound the memory is reported in `Stats::memory_bound_by`.
      std::size_t max_live_precursors = 0;

      /// Memory budget for the live blocks, BYTES. Non-zero overrides
      /// `max_live_precursors`, which is derived from it.
      ///
      /// `max_live_precursors = 0` means "bounded only by retention-time
      /// overlap", and that bound is a property of the RUN, not of the library:
      /// it does not tighten as the library grows. On a 4,986,319-precursor
      /// library with no iRT map -- so every precursor predicted across the
      /// whole gradient -- essentially the entire library is live at once, and
      /// the extractor was OOM-killed at 588 GB on a 995 GB node.
      ///
      /// A cap stated in precursors cannot be chosen without knowing the
      /// transition count and window width, which is why nobody set one. A cap
      /// stated in BYTES can: the extractor knows the mean cells per precursor
      /// and how many planes each carries, so it inverts the budget into a
      /// precursor count itself and reports what it chose.
      ///
      /// Per live precursor the cost is
      ///     valid_transitions x cycles_in_window x 4 B x planes
      /// where planes is 1, +2 with `collect_mass_residuals`, +2 with
      /// `collect_im_residuals` -- so 5 with both, which is the default.
      std::size_t live_memory_budget_bytes = 0;

      /// How many spectra are decoded and held at once.
      ///
      /// This is the largest single term in the run's memory, and it was found
      /// by profile rather than by reading: a tcmalloc heap profile of a
      /// 1,200-precursor IH1 run put **5.57 GiB of a 9.45 GiB live peak** in the
      /// three `assign` calls of `MzPeakSource::peaks` -- 2.79 / 1.39 / 1.39 GiB
      /// across `mz` (double), `intensity` and `ion_mobility` (float), exactly
      /// the 8:4:4 ratio of their element sizes.
      ///
      /// The cost is linear in this number and it is NOT the sliding window:
      /// six other explanations were measured and refuted first (per-thread
      /// parquet buffers, page cache, glibc fragmentation, arena count, sparse
      /// chromatograms, and mzPeak's own row-group cache, which is bounded).
      /// See doc/11-memory-and-compaction-plan.md.
      ///
      /// Measured on IH1 with a 1,200-precursor library at 4 threads, where the
      /// only variable was this number:
      ///
      /// | block | peak RSS | wall |
      /// |------:|---------:|-----:|
      /// | 1024  | 10.34 GiB | 7:00.3 |
      /// |  256  |  2.89 GiB | 6:53.7 |
      /// |   64  |  1.40 GiB | 6:50.3 |
      ///
      /// There is no tradeoff to balance: smaller is both smaller and faster,
      /// so 1024 was pure waste. The saving is better than linear because
      /// Arrow's decode buffers scale with the batch as well as our peak arrays.
      ///
      /// 256 rather than 64 because of a DIFFERENT measurement, already in this
      /// file: at 128 spectra and 64 threads each worker got two spectra and
      /// thread creation cost more than the matching (0.19 s against 0.14 s
      /// single-threaded). The arms above ran at 4 threads and so cannot see
      /// that. 256 keeps four spectra per worker at 64 threads and still takes
      /// 3.6x of the available 7.4x. Lower it when threads are few.
      std::size_t decode_block = 256;

      /// Extract only every Nth precursor. 1 is all of them.
      ///
      /// For pass 1 of the two-pass workflow, which extracts solely to harvest
      /// RT anchors. Indices are NOT renumbered -- a strided pass still reports
      /// full-library precursor indices, so the anchor harvest and the iRT it
      /// is fitted against cannot drift apart.
      std::size_t precursor_stride = 1;

      /// Which residue class the stride keeps, so equal-sized subsets with
      /// DIFFERENT members can be compared. Exists to separate "how many
      /// anchors" from "which anchors" when measuring how stable the RT fit is.
      std::size_t precursor_offset = 0;

      /// Optional per-precursor keep-mask, in FULL-library indexing, or null
      /// for all of them. Non-null means `PrecursorPrefilter` has already
      /// decided this precursor cannot be supported by the run.
      ///
      /// Like the stride, this does NOT renumber: a filtered pass still reports
      /// full-library indices, so anchors and the iRT they are fitted against
      /// cannot drift apart. The caller owns the storage and it must outlive
      /// the extraction.
      const char* precursor_keep = nullptr;

      /// Optional per-precursor disposition table, shared with the scorer's
      /// `PeakGroupScorer::Options::terminal_reason` and indexed the same way.
      ///
      /// The extractor writes ONLY the two reasons it alone can know: a
      /// precursor no isolation window covers, and one `precursor_keep`
      /// dropped. Everything downstream of a chromatogram being emitted is the
      /// scorer's to record. Null disables it.
      std::uint8_t* terminal_reason = nullptr;
      /// Must equal PeakGroupScorer::TerminalReason::NoWindowCoverage.
      static constexpr std::uint8_t kNoWindowCoverage = 9;
      /// Must equal PeakGroupScorer::TerminalReason::PrefilterExcluded.
      static constexpr std::uint8_t kPrefilterExcluded = 10;

      /// Maps the library's iRT onto this run's retention time, in seconds:
      /// rt = irt_slope * irt + irt_intercept.
      ///
      /// This is what lets a precursor be extracted around where it should
      /// elute instead of across the whole run, and it is the difference
      /// between 1.7e10 points and something that fits in memory. It cannot be
      /// derived here -- it needs the run to have been searched once -- so the
      /// caller supplies it. A slope of 0 disables the restriction and every
      /// precursor is extracted over the full range, which is the old
      /// behaviour and is kept because it is the only option before a first
      /// pass exists.
      double irt_slope = 0.0;
      double irt_intercept = 0.0;

      /// Half-width of the extraction window, seconds. Sized from the
      /// retention-time residual at a stated quantile -- NOT from its
      /// standard deviation, because a window has to cover the tail it is
      /// meant to catch.
      double rt_window_seconds = 60.0;

      /// Worker threads for the matching. 0 uses the hardware concurrency.
      /// The decode stays serial: the reader's thread-safety is unverified,
      /// and decode is a shared cost per spectrum rather than per transition.
      unsigned threads = 0;

      /// Restrict to this retention-time range, seconds. Both zero means all.
      double rt_low = 0.0;
      double rt_high = 0.0;

      /// Match a peak's ion mobility against the WINDOW's limits when the run
      /// carries it. A diaPASEF frame holds several windows over disjoint
      /// mobility ranges, so ignoring this mixes them.
      ///
      /// This gates the frame band and nothing else. `precursor_im_window` is
      /// a separate test with a separate switch, and turning this off must not
      /// turn that off -- it used to, because the peak's mobility was only read
      /// inside this branch and a NaN mobility skips the per-precursor test.
      /// Switching off both filters while claiming to switch off one is worse
      /// than either: this flag exists to be the control arm that says what
      /// mobility filtering is worth.
      bool use_ion_mobility = true;

      /// Report progress every this many spectra, 0 to stay quiet.
      std::size_t progress_every = 500;
    };

    struct Stats
    {
      std::size_t spectra_read = 0;
      std::size_t precursors = 0;
      std::size_t transitions = 0;
      std::size_t points = 0;
      std::size_t nonzero_points = 0;
      double decode_seconds = 0.0;
      double match_seconds = 0.0;
      double index_seconds = 0.0;
      /// Allocating, zeroing and releasing the live chromatograms.
      double assemble_seconds = 0.0;
      /// Time inside the sink -- copying, or scoring and discarding.
      double sink_seconds = 0.0;

      /// Precursors not extracted because their predicted elution fell outside
      /// the run entirely, or could not be predicted at all (a NaN iRT).
      std::size_t outside_rt_range = 0;
      /// Mean transitions live at one cycle, which is what the inverted match
      /// actually costs per spectrum.
      double mean_live_transitions = 0.0;

      /// What the sliding window actually bought, measured rather than
      /// assumed: the largest number of precursors -- and of points -- whose
      /// retention-time windows overlap at one instant. `peak_live_points x 4`
      /// bytes is the chromatogram term of peak RSS.
      std::size_t peak_live_precursors = 0;
      std::uint64_t peak_live_points = 0;
      /// Precursors that were extracted at all, i.e. had a window and a
      /// non-empty cycle range.
      std::size_t precursors_extracted = 0;
      /// Passes over the run. More than one means the live set did not fit
      /// under `max_live_precursors` and the library was chunked.
      std::size_t chunks = 1;
      /// Spectra decoded across all chunks, against `spectra_read` for one
      /// pass. Equal when there is one chunk.
      std::size_t spectra_decoded = 0;

      /// Also reported on the returned `Chromatograms`; here so a streaming
      /// caller, which never sees one, still gets them.
      std::size_t precursors_without_window = 0;
      std::size_t precursors_in_several_windows = 0;

      /// How many `precursor_keep` excluded. Reported because a filter that
      /// silently discards is indistinguishable from a search that found
      /// nothing -- doc/08's third safety rule.
      std::size_t precursors_prefiltered = 0;

      /// Which of the two mechanisms set the peak: "retention-time overlap"
      /// when the live set fitted, "precursor cap (N), C chunks" when it did
      /// not. Said out loud because the two have completely different costs --
      /// the first is free, the second buys memory with decode passes.
      std::string memory_bound_by;
      /// How the byte budget was inverted into a precursor cap, when one was
      /// given. Empty when the cap came from -max_live_precursors or from
      /// retention-time overlap alone.
      std::string live_budget_note;
    };

    /// Extract into a sink, holding only what is live.
    static void extract(const Library& library, SpectrumSource& source,
                        const Options& options, ChromatogramSink& sink,
                        Stats* stats = nullptr);

    /// Extract into one flat `Chromatograms`, as this class always did.
    ///
    /// Kept because `-out_chrom` and every diagnostic built on it need the
    /// whole thing at once. It is the streaming form with a collecting sink,
    /// so it is bounded by the library and says so: at IH1's geometry a full
    /// human library is 274 GiB here and a few GiB through the sink.
    static Chromatograms extract(const Library& library, SpectrumSource& source,
                                 const Options& options, Stats* stats = nullptr);
  };

} // namespace ODIA
