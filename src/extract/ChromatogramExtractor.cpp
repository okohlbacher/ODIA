// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/ChromatogramExtractor.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <thread>

namespace ODIA
{

  std::size_t Chromatograms::footprintBytes() const
  {
    const auto vec = [](const auto& v) { return v.capacity() * sizeof(v[0]); };
    std::size_t axes_bytes = 0;
    for (const auto& a : axes) { axes_bytes += vec(a); }
    return axes_bytes + vec(precursor_axis) + vec(precursor_axis_begin) +
           vec(precursor_cycles) + vec(precursor_transition_begin) + vec(begin) +
           vec(count) + vec(intensity);
  }

  namespace
  {
    bool sameWindow(const IsolationWindow& a, const IsolationWindow& b)
    {
      return std::abs(a.mz_low - b.mz_low) < 1e-6 && std::abs(a.mz_high - b.mz_high) < 1e-6;
    }

    /// The spectra of one isolation window, in acquisition order.
    ///
    /// A DIA run cycles through its windows, so a window's spectra ARE its
    /// retention-time axis: one entry per cycle. Measured on 12_80 the cycle is
    /// regular to 0.24%, which is why the index below can be computed rather
    /// than searched -- but one window there also drops a spectrum, giving a
    /// 1.66 s gap where 0.83 s was expected, and a dropped cycle shifts a pure
    /// formula by one FOR THE REST OF THE RUN. So the estimate is corrected
    /// against the real times rather than trusted.
    struct CycleAxis
    {
      std::vector<std::size_t> spectrum;   ///< index into SpectrumSource::spectra()
      /// POSITION of each cycle's spectrum in SpectrumSource::spectra(), the
      /// coordinate the pass decodes and batches in.
      ///
      /// Everything that decides WHICH spectra a pass reads is computed from
      /// these, never from `rt`. `rt` is float, the run's times are double, and
      /// a float key searched in the double times lands one spectrum early or
      /// late whenever the rounding goes the wrong way -- which is how chunked
      /// extraction used to skip the first or last frame of a chunk and leave
      /// those cells at zero (2,992 of 6,020 fixture spectra round up, 3,028
      /// down; none is exact). A position cannot round.
      std::vector<std::size_t> position;
      std::vector<float> rt;               ///< ascending, seconds
      double rt0 = 0.0;
      double step = 1.0;

      void finish()
      {
        if (rt.size() < 2) { return; }
        rt0 = rt.front();
        step = (rt.back() - rt.front()) / double(rt.size() - 1);
        if (!(step > 0.0)) { step = 1.0; }
      }

      /// First cycle at or after @p t. O(1) estimate plus a bounded local
      /// walk, so a dropped cycle costs a step or two rather than correctness.
      std::size_t lowerBound(double t) const
      {
        if (rt.empty() || t <= rt.front()) { return 0; }
        if (t > rt.back()) { return rt.size(); }
        auto i = static_cast<std::ptrdiff_t>((t - rt0) / step);
        i = std::clamp<std::ptrdiff_t>(i, 0, std::ptrdiff_t(rt.size()) - 1);
        while (i > 0 && rt[std::size_t(i) - 1] >= t) { --i; }
        while (std::size_t(i) < rt.size() && rt[std::size_t(i)] < t) { ++i; }
        return std::size_t(i);
      }
    };

    /// Transitions of one window, ordered by m/z and bucketed on log m/z.
    ///
    /// Bucketing on the LOG of m/z is what makes the lookup O(1): a ppm
    /// tolerance is a constant width in log space, so one bucket covers one
    /// tolerance at every mass, where a linear grid would be far too coarse at
    /// 200 Th and far too fine at 1800.
    ///
    /// The point of this index is to invert the match. Probing the spectrum
    /// once per transition is O(N log M) with N transitions and M peaks, and N
    /// exceeds M by three orders of magnitude at proteome scale -- 4.6 M
    /// against ~1,100. Asking instead which transitions each PEAK could belong
    /// to is O(M + matches).
    struct MzIndex
    {
      std::vector<double> mz;              ///< ascending
      /// Which precursor's live block this transition writes into, and which
      /// row of it.
      ///
      /// This replaces the absolute offset into one flat point array that the
      /// index used to carry, and the replacement is the whole point: an
      /// absolute offset only exists if the whole array exists. A slot is
      /// resolved through `LiveSlot`, which holds a pointer that is null until
      /// the pass reaches the precursor's window and null again once it has
      /// left -- so the index survives the block being allocated and freed
      /// underneath it, and does not have to be rebuilt when it is.
      std::vector<std::uint32_t> slot;
      /// A precursor has 12 transitions here and the library format allows far
      /// fewer than 65,535, so the row is 16 bits. At 51 M transitions the two
      /// bytes saved against a uint32 are 102 MB.
      std::vector<std::uint16_t> row;
      /// Expected 1/K0 of each transition's precursor. NaN when the library has
      /// none, which disables the per-precursor mobility test for it rather
      /// than rejecting it -- absent information is not evidence of mismatch.
      std::vector<float> precursor_im;
      std::vector<std::uint32_t> bucket;        ///< bucket -> first entry
      double log_base = 0.0, log_lo = 0.0;

      std::size_t bucketOf(double m) const
      {
        const auto b = static_cast<std::ptrdiff_t>((std::log(m) - log_lo) / log_base);
        return std::size_t(std::clamp<std::ptrdiff_t>(b, 0, std::ptrdiff_t(bucket.size()) - 2));
      }

      void build(double tolerance_ppm)
      {
        if (mz.empty()) { return; }
        log_base = std::log1p(tolerance_ppm * 1e-6);
        if (!(log_base > 0.0)) { log_base = 1e-6; }
        log_lo = std::log(mz.front());
        const auto n = static_cast<std::size_t>(
                         (std::log(mz.back()) - log_lo) / log_base) + 2;
        bucket.assign(n + 1, static_cast<std::uint32_t>(mz.size()));
        // Filled backwards so each bucket points at its FIRST entry and empty
        // buckets inherit the next non-empty one, which keeps the lookup
        // branchless.
        for (std::size_t i = mz.size(); i-- > 0;)
        {
          bucket[bucketOf(mz[i])] = static_cast<std::uint32_t>(i);
        }
        for (std::size_t b = bucket.size() - 1; b-- > 0;)
        {
          bucket[b] = std::min(bucket[b], bucket[b + 1]);
        }
      }
    };

    /// One precursor's place in the run: which window carried it, and over
    /// which of that window's cycles it is expected to elute.
    struct Assignment
    {
      std::uint32_t precursor = 0;
      std::uint32_t window = 0;
      std::uint32_t lo = 0, hi = 0;    ///< half-open cycle range on that window
      /// Transitions with a representable product m/z, which is how many rows
      /// the live block has. A transition with none can never receive a point,
      /// so giving it a row would reserve memory nothing can write to.
      std::uint32_t valid = 0;
    };

    /// Where a precursor's points are RIGHT NOW.
    ///
    /// `base` is null outside the precursor's retention-time window, which is
    /// most of the run for most precursors and is the entire memory argument.
    /// The cycle range is kept here too so the match loop resolves a peak in
    /// one cache line rather than chasing back into the assignment table.
    struct LiveSlot
    {
      float* base = nullptr;
      /// Intensity-weighted m/z deviation per cell: sum(intensity * ppm) in
      /// `ppm_num` and sum(intensity) in `ppm_den`, reduced where it is read.
      ///
      /// TWO planes, after a wrong shortcut. The first version stored a single
      /// plane holding "the deviation of the peak that won the cell", justified
      /// as exact under Aggregate::Max -- but **Max is not the default**;
      /// Options::aggregate defaults to Sum and so does the CLI. Under Sum every
      /// peak "wins", so that plane held the LAST matching peak in spectrum
      /// order, which is arbitrary: peak arrays are not sorted, so it was
      /// neither the brightest nor the nearest. Found by external review.
      ///
      /// The weighted form is correct under both modes and is the right
      /// quantity regardless: a bright fragment's centroid is better determined
      /// than a dim one's, and on DIA data the dim end is where interference
      /// lives.
      ///
      /// Race-free for the same reason `base` is: a thread owns a spectrum, and
      /// a spectrum owns a distinct cycle, so no two threads touch one cell.
      ///
      /// Null unless Options::collect_mass_residuals. Live blocks are the
      /// extractor's dominant memory term, so this is opt-in.
      float* ppm_num = nullptr;
      float* ppm_den = nullptr;
      /// The same construction for ION MOBILITY. Per-(row, cycle) for the same
      /// reason: the match loop is threaded over spectra and is race-free only
      /// because each spectrum owns a distinct cycle.
      float* im_num = nullptr;
      float* im_den = nullptr;
      std::uint32_t lo = 0, hi = 0;
    };

    /// How many precursors are live in each match batch of the run, for the
    /// chunk planner: add a precursor's span [first, last], ask the largest
    /// count inside a span. A segment tree with the pending add kept on the
    /// node (never pushed down), so both are O(log batches) -- the planner asks
    /// once per precursor, millions of times, over ~10^3 batches.
    ///
    /// reset() is O(1): every node carries the generation it was last written
    /// in, and a node from an older generation reads as zero. Clearing the
    /// whole tree at every chunk cut made planning O(chunks x batches) on top
    /// of the O(N log batches) of the adds -- at a cap of 1 with everything
    /// overlapping that is one full clear per precursor.
    class LiveCount
    {
    public:
      explicit LiveCount(std::size_t batches)
        : n_(std::max<std::size_t>(batches, 1)), max_(4 * n_, 0), add_(4 * n_, 0),
          gen_of_(4 * n_, 0) {}

      void reset() { ++gen_; }
      void add(std::size_t first, std::size_t last) { add(1, 0, n_ - 1, first, last); }
      std::uint32_t max(std::size_t first, std::size_t last) const
      { return max(1, 0, n_ - 1, first, last); }

    private:
      bool fresh(std::size_t node) const { return gen_of_[node] == gen_; }
      std::uint32_t maxOf(std::size_t node) const { return fresh(node) ? max_[node] : 0u; }
      std::uint32_t addOf(std::size_t node) const { return fresh(node) ? add_[node] : 0u; }
      void touch(std::size_t node)
      {
        if (fresh(node)) { return; }
        max_[node] = 0;
        add_[node] = 0;
        gen_of_[node] = gen_;
      }

      void add(std::size_t node, std::size_t lo, std::size_t hi,
               std::size_t first, std::size_t last)
      {
        if (last < lo || hi < first) { return; }
        touch(node);
        if (first <= lo && hi <= last) { ++max_[node]; ++add_[node]; return; }
        const std::size_t mid = lo + (hi - lo) / 2;
        add(2 * node, lo, mid, first, last);
        add(2 * node + 1, mid + 1, hi, first, last);
        max_[node] = add_[node] + std::max(maxOf(2 * node), maxOf(2 * node + 1));
      }
      std::uint32_t max(std::size_t node, std::size_t lo, std::size_t hi,
                        std::size_t first, std::size_t last) const
      {
        if (last < lo || hi < first || !fresh(node)) { return 0; }
        if (first <= lo && hi <= last) { return max_[node]; }
        const std::size_t mid = lo + (hi - lo) / 2;
        return addOf(node) + std::max(max(2 * node, lo, mid, first, last),
                                      max(2 * node + 1, mid + 1, hi, first, last));
      }

      std::size_t n_;
      std::vector<std::uint32_t> max_, add_;
      std::vector<std::uint64_t> gen_of_;
      std::uint64_t gen_ = 0;
    };

    /// Fixed-size blocks, reused rather than returned to the allocator.
    ///
    /// Nearly every precursor asks for the same size -- 12 transitions by the
    /// window's cycle count -- so a free list per size hits essentially always,
    /// and the pass does one allocation per concurrently-live precursor for the
    /// whole run instead of one per precursor. That matters at 4.26 M
    /// precursors: 4.26 M allocate/free pairs of ~40 KB through malloc is where
    /// a general allocator fragments, and fragmentation would put back exactly
    /// the memory this change exists to remove.
    ///
    /// Blocks are handed out ZEROED, because the match accumulates into them
    /// and a reused block still holds the previous precursor's peaks.
    class BlockPool
    {
    public:
      float* take(std::size_t n)
      {
        // A precursor all of whose transitions lack a product m/z has no rows.
        // It is still emitted -- a consumer counting precursors that yielded
        // nothing must see it -- so it needs a base that is not null.
        if (n == 0) { return &dummy_; }
        auto& bin = free_[n];
        float* p = nullptr;
        if (!bin.empty())
        {
          p = bin.back();
          bin.pop_back();
        }
        else
        {
          owned_.push_back(std::make_unique_for_overwrite<float[]>(n));
          p = owned_.back().get();
          reserved_ += n;
        }
        std::fill_n(p, n, 0.0f);
        live_ += n;
        peak_ = std::max(peak_, live_);
        return p;
      }

      void give(float* p, std::size_t n)
      {
        if (n == 0) { return; }
        free_[n].push_back(p);
        live_ -= n;
      }

      std::uint64_t reservedPoints() const { return reserved_; }
      std::uint64_t peakPoints() const { return peak_; }

    private:
      std::unordered_map<std::size_t, std::vector<float*>> free_;
      std::vector<std::unique_ptr<float[]>> owned_;
      std::uint64_t live_ = 0, peak_ = 0, reserved_ = 0;
      float dummy_ = 0.0f;
    };
  } // namespace

  void ChromatogramCollector::begin(const ChromatogramLayout& layout)
  {
    // Refused with the number rather than left to a bad_alloc from an
    // allocation whose size nothing reported. This sink keeps every point, so
    // this is knowable exactly here, and it is the point at which a run that
    // cannot fit should say so -- and say what to do instead, which is now a
    // real alternative rather than advice to shrink the problem.
    if (layout.points > c_.intensity.max_size())
    {
      throw std::runtime_error(
        "chromatograms need " + std::to_string(layout.points) + " points (" +
        std::to_string(layout.points * sizeof(float) / (1024ull * 1024 * 1024)) +
        " GiB) to be held at once. Collecting them all is what -out_chrom "
        "needs; scoring does not. Drop -out_chrom to score on the fly, narrow "
        "the retention-time window (-rt_window), or reduce the precursor "
        "count.");
    }

    if (layout.counts == nullptr)
    {
      throw std::logic_error("the chromatogram collector was given no layout counts");
    }

    // Said BEFORE the allocation, not after, because after is an OOM kill with
    // no message. The extractor's own summary reports this too, but only once
    // the pass has finished -- which is too late to be told what went wrong.
    constexpr std::uint64_t LOUD_GIB = 8;
    if (layout.points * sizeof(float) > LOUD_GIB * 1024 * 1024 * 1024)
    {
      std::cerr << "  holding every chromatogram at once: "
                << double(layout.points) * sizeof(float) / 1073741824.0
                << " GiB. Only a consumer that needs them all -- -out_chrom --"
                   " requires this; scoring does not.\n";
    }

    c_ = Chromatograms{};
    if (layout.axes != nullptr) { c_.axes = *layout.axes; }
    c_.count = *layout.counts;
    // Per PRECURSOR, not per transition. The extractor guarantees a precursor
    // is extracted from one window over one cycle range, so these were three
    // identical values repeated across its dozen transitions -- 612 MB of the
    // index at 51.1 M transitions, for 4.26 M distinct facts.
    c_.precursor_axis.assign(layout.precursors, 0);
    c_.precursor_axis_begin.assign(layout.precursors, 0);
    c_.precursor_cycles.assign(layout.precursors, 0);
    c_.precursor_transition_begin.assign(layout.precursors, 0);

    // The offsets are the prefix sum in TRANSITION order, which is what every
    // reader of a `Chromatograms` assumes and what the flat form guarantees.
    // Precursors arrive here in retention-time order instead, so the layout has
    // to be fixed before the first one does -- laying it out as they arrive
    // would produce the same points at different offsets, and a file written
    // from it would still be identical while the structure quietly was not.
    c_.begin.assign(layout.transitions, 0);
    std::uint64_t running = 0;
    for (std::size_t j = 0; j < layout.transitions; ++j)
    {
      c_.begin[j] = running;
      running += c_.count[j];
    }
    c_.intensity.assign(running, 0.0f);
  }

  void ChromatogramCollector::repairPrecursorIndex()
  {
    // `precursor_transition_begin` is zero-initialised and written only by
    // accept(), which returns early for a precursor that extracted nothing --
    // no isolation window, no points, or all-zero. Those leave ZERO HOLES in an
    // array that `precursorOf()` binary-searches with std::upper_bound, and
    // upper_bound on a non-monotonic array returns an arbitrary index. The
    // symptom is a precursor whose retention times are read from axis position
    // 0, i.e. the start of the run: on one 10,000-precursor arm 35.16% of the
    // written traces carried the run's first cycles as their RT column while
    // their intensities were correct, because intensity goes through the
    // separately-correct `begin[]` prefix sum and the scorer indexes
    // `precursor_axis_begin` directly. A TSV bug, not a scoring bug -- but a
    // silent one, and anything that reads the RT column (an apex-centred crop,
    // an RT-window check, a cross-engine join on time) is wrong for those rows.
    //
    // Filled from the RIGHT, not the left. A hole must take the begin of the
    // NEXT precursor so it occupies zero width; forward-filling it with the
    // PREVIOUS value would widen the predecessor's span and hand that
    // predecessor's own transitions to the empty precursor instead. A precursor
    // that extracted nothing owns no transition with data, so a zero-width
    // entry is exactly right for every transition that can be queried.
    auto& b = c_.precursor_transition_begin;
    if (b.empty()) { return; }
    std::uint32_t next = static_cast<std::uint32_t>(c_.count.size());
    for (std::size_t i = b.size(); i-- > 0;)
    {
      if (b[i] == 0 && i != 0) { b[i] = next; }
      else { next = b[i]; }
    }
  }


  void ChromatogramCollector::accept(const PrecursorChromatogram& trace)
  {
    if (!trace.extracted()) { return; }
    // Once per precursor. These used to be written once per transition with
    // the same three values each time.
    if (trace.precursor < c_.precursor_axis.size())
    {
      c_.precursor_axis[trace.precursor] = trace.axis;
      c_.precursor_axis_begin[trace.precursor] = trace.axis_begin;
      c_.precursor_cycles[trace.precursor] = trace.cycles;
      c_.precursor_transition_begin[trace.precursor] = trace.transition_begin;
    }
    for (std::uint32_t k = 0; k < trace.transition_count; ++k)
    {
      const std::uint32_t n = trace.pointCount(k);
      if (n == 0) { continue; }
      const std::uint32_t tr = trace.transition_begin + k;
      const float* from = trace.trace(k);
      std::copy(from, from + n, c_.intensity.begin() +
                                 static_cast<std::ptrdiff_t>(c_.begin[tr]));
    }
  }


  Chromatograms ChromatogramExtractor::extract(const Library& library,
                                               SpectrumSource& source,
                                               const Options& options, Stats* stats)
  {
    Stats local;
    Stats& st = stats == nullptr ? local : *stats;
    ChromatogramCollector collector;
    extract(library, source, options, collector, &st);
    Chromatograms out = collector.take();   // take() repairs the precursor index
    // Counted by the extractor, reported on the object, because that is where
    // every caller has always read them.
    out.precursors_without_window = st.precursors_without_window;
    out.precursors_in_several_windows = st.precursors_in_several_windows;
    return out;
  }

  void ChromatogramExtractor::extract(const Library& library, SpectrumSource& source,
                                      const Options& options, ChromatogramSink& sink,
                                      Stats* stats)
  {
    Stats local;
    Stats& st = stats == nullptr ? local : *stats;
    st = Stats{};
    const auto t_index = std::chrono::steady_clock::now();

    const auto& p = library.precursors();
    const auto& t = library.transitions();
    const auto& windows = source.windows();
    const auto& info = source.spectra();

    const std::size_t n_prec = options.max_precursors == 0
                                 ? library.precursorCount()
                                 : std::min(options.max_precursors, library.precursorCount());

    // The retention-time axis of every window, which is also its cycle index.
    std::vector<CycleAxis> axis(windows.size());
    std::vector<std::uint32_t> window_of(info.size(),
                                         std::numeric_limits<std::uint32_t>::max());
    std::vector<std::uint32_t> cycle_of(info.size(), 0);
    for (std::size_t si = 0; si < info.size(); ++si)
    {
      for (std::size_t w = 0; w < windows.size(); ++w)
      {
        if (sameWindow(info[si].window, windows[w]))
        {
          window_of[si] = static_cast<std::uint32_t>(w);
          cycle_of[si] = static_cast<std::uint32_t>(axis[w].rt.size());
          axis[w].spectrum.push_back(info[si].index);
          axis[w].position.push_back(si);
          axis[w].rt.push_back(static_cast<float>(info[si].retention_time));
          break;
        }
      }
    }
    for (auto& a : axis) { a.finish(); }

    // The caller's retention-time range, resolved to a cycle range per window
    // and to a spectrum range for the pass. Resolving it BEFORE decoding is the
    // difference between reading the slice and reading the run: filtering
    // afterwards decodes everything and throws most of it away.
    std::vector<std::size_t> global_lo(windows.size(), 0), global_hi(windows.size(), 0);
    for (std::size_t w = 0; w < windows.size(); ++w)
    {
      global_hi[w] = axis[w].rt.size();
      if (options.rt_high > options.rt_low)
      {
        global_lo[w] = axis[w].lowerBound(options.rt_low);
        global_hi[w] = axis[w].lowerBound(options.rt_high);
      }
    }

    // One axis per window, shared by every transition extracted from it. The
    // per-point timestamps this replaces were the same ~1,000 values repeated
    // once per transition.
    std::vector<std::vector<float>> axes(axis.size());
    for (std::size_t w = 0; w < axis.size(); ++w) { axes[w] = axis[w].rt; }

    const std::size_t n_trans =
      n_prec == 0 ? 0 : p.transition_begin[n_prec - 1] + p.transition_count[n_prec - 1];

    const bool restrict_rt = options.irt_slope != 0.0;
    // Three runs -- uncalibrated, calibrated at 600 s, calibrated at 60 s --
    // extracted 119,088,506 transitions each, identical to the digit. The RT
    // window has never restricted anything, and reading the code cannot say why.
    std::fprintf(stderr,
      "RT-RESTRICT DIAGNOSTIC: restrict_rt=%d irt_slope=%.6g irt_intercept=%.6g "
      "rt_window_seconds=%.6g\n",
      int(restrict_rt), options.irt_slope, options.irt_intercept,
      options.rt_window_seconds);
    std::uint64_t total_points = 0;
    std::uint64_t live_sum = 0;

    // Which window carries each precursor, and over which cycles.
    //
    // ONE window, even though several may cover it: the schemes on 12_80 and
    // IH1 overlap adjacent windows by 1.0 Th, which puts 1.1% of precursors in
    // two. Each of those is a separate measurement of the same ion at
    // interleaved times -- and interleaved times is exactly why they cannot be
    // one chromatogram. Concatenating them (which this replaces) made a
    // transition's points span two axes, so the single (axis_of, axis_begin)
    // pair the representation stores could only name one of them: every point
    // of the first window then read the second window's axis, past its end.
    // Even with a per-run axis table it would still hand the scorer a trace
    // that runs backwards in time at the seam and elutes the peptide twice,
    // and the scorer reads these points as one time series.
    //
    // So: the window whose centre the precursor is nearest, which is where it
    // is furthest from the edge and best transmitted. The comparison is strict,
    // so an exact tie would fall to the lower index -- but m/z arrives here
    // through 10 nTh fixed point, which makes an exact tie a thing to state
    // rather than a thing to rely on. The discarded measurement is counted, not
    // dropped silently.
    std::vector<Assignment> assignments;
    assignments.reserve(n_prec);
    // Which precursors got one. A precursor that did not is still handed to
    // the sink, empty, at the end: a consumer counting "precursors that yielded
    // nothing" used to see them as rows of zero-length in one flat array, and
    // must still see them now that the array does not exist.
    std::vector<char> assigned(n_prec, 0);

    for (std::size_t i = 0; i < n_prec; ++i)
    {
      // Pass 1 exists only to harvest retention-time anchors, and it needs very
      // few: 692 came from 2,450 precursors, and -min_anchors defaults to 20.
      // Extracting the whole library to find them costs memory proportional to
      // the library -- at 4.26 M precursors over the whole run that is ~274 GiB,
      // which is why the calibrated pass-2 window did not unblock the large
      // library on its own.
      //
      // A STRIDE rather than a random sample or a sub-library, because the
      // index `i` must stay a full-library index: anchors are harvested as
      // (original_irt[i], best[i]->apex_rt), and remapping them is how the RT
      // map would silently get fitted to the wrong iRTs.
      if (options.precursor_stride > 1 &&
          (i % options.precursor_stride) != (options.precursor_offset % options.precursor_stride))
      {
        continue;
      }
      // The prefilter's verdict, for the same reason and with the same
      // indexing rule as the stride above.
      if (options.precursor_keep != nullptr && options.precursor_keep[i] == 0)
      {
        ++st.precursors_prefiltered;
        if (options.terminal_reason)
        { options.terminal_reason[i] = Options::kPrefilterExcluded; }
        continue;
      }
      const double mz = fromFixed(p.mz[i]);
      std::size_t best = windows.size();
      double best_offset = std::numeric_limits<double>::infinity();
      std::size_t covering = 0;
      for (std::size_t w = 0; w < windows.size(); ++w)
      {
        if (!windows[w].contains(mz) || axis[w].rt.empty()) { continue; }
        ++covering;
        const double offset_from_centre = std::abs(mz - windows[w].centre());
        if (offset_from_centre < best_offset) { best_offset = offset_from_centre; best = w; }
      }
      if (covering == 0)
      {
        ++st.precursors_without_window;
        if (options.terminal_reason)
        { options.terminal_reason[i] = Options::kNoWindowCoverage; }
        continue;
      }
      if (covering > 1) { ++st.precursors_in_several_windows; }

      const std::size_t w = best;
      std::size_t lo = global_lo[w], hi = global_hi[w];
      if (restrict_rt)
      {
        const double centre = options.irt_slope * double(p.irt[i]) + options.irt_intercept;
        if (std::isnan(centre)) { ++st.outside_rt_range; continue; }
        lo = std::max(lo, axis[w].lowerBound(centre - options.rt_window_seconds));
        hi = std::min(hi, axis[w].lowerBound(centre + options.rt_window_seconds));
      }
      if (lo >= hi) { ++st.outside_rt_range; continue; }

      Assignment a;
      a.precursor = static_cast<std::uint32_t>(i);
      a.window = static_cast<std::uint32_t>(w);
      a.lo = static_cast<std::uint32_t>(lo);
      a.hi = static_cast<std::uint32_t>(hi);
      for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
      {
        const std::size_t j = p.transition_begin[i] + k;
        // A transition whose product m/z could not be represented has no mass
        // to match against, so the index build below skips it and it can never
        // receive a point. Counting it here anyway gave it a full-length point
        // run that nothing ever filled, the default axis_of = 0, and therefore
        // a read off the end of window 0's axis whenever window 0 was the
        // shorter one. Library::invalidMzTransitionCount() exists because this
        // input is expected -- a neutral loss can drive a small fragment to or
        // below zero -- so it is skipped, not rejected.
        if (t.product_mz[j] == MZ_INVALID) { continue; }
        ++a.valid;
        total_points += hi - lo;
        live_sum += hi - lo;
      }
      // The m/z index addresses a transition by its row in the precursor's
      // block, in 16 bits. Twelve is what a library holds here and the format
      // allows far fewer than 65,535, but a silent wrap would send every peak
      // of the 65,536th transition into the first one's trace, so it is
      // refused with the number instead.
      if (a.valid > std::numeric_limits<std::uint16_t>::max())
      {
        throw std::runtime_error(
          "precursor " + std::to_string(i) + " has " + std::to_string(a.valid) +
          " extractable transitions; the extractor addresses at most " +
          std::to_string(std::numeric_limits<std::uint16_t>::max()) + " per precursor");
      }
      assigned[i] = 1;
      assignments.push_back(a);
    }

    const std::size_t n_slots = assignments.size();
    st.precursors_extracted = n_slots;

    // Four bytes a transition, and only for a sink that lays out one flat CSR.
    // A streaming sink never sees it, which is the difference between a term
    // proportional to the library and none.
    std::vector<std::uint32_t> layout_counts;
    if (sink.needsLayoutCounts())
    {
      layout_counts.assign(n_trans, 0);
      for (const Assignment& a : assignments)
      {
        const std::uint32_t tb = p.transition_begin[a.precursor];
        for (std::uint32_t k = 0; k < p.transition_count[a.precursor]; ++k)
        {
          if (t.product_mz[tb + k] == MZ_INVALID) { continue; }
          layout_counts[tb + k] = a.hi - a.lo;
        }
      }
    }

    ChromatogramLayout layout;
    layout.transitions = n_trans;
    layout.precursors = library.precursorCount();
    layout.points = total_points;
    layout.axes = &axes;
    layout.counts = sink.needsLayoutCounts() ? &layout_counts : nullptr;
    sink.begin(layout);

    // The pass's coordinates: which spectra it reads, in which blocks it
    // decodes them, and in which batches it matches them. Fixed here, BEFORE
    // the live set is measured and the library chunked, because both are
    // properties of these coordinates and nothing else.
    //
    // Large enough that spawning workers is amortised. At 128 spectra and 64
    // threads each worker got two spectra and thread creation cost more than
    // the matching did -- measured 0.19 s against 0.14 s single-threaded.
    // Was a constant; now `Options::decode_block`, because a heap profile put
    // 5.57 GiB of a 9.45 GiB live peak in the block's `SpectrumPeaks` copies.
    // Zero keeps the historical default rather than degenerating to no block.
    //
    // The block decides ONLY how many spectra are decoded and held at once. It
    // decides nothing a sink can observe: see MATCH_BATCH.
    const std::size_t BLOCK = options.decode_block ? options.decode_block : 256;
    // Matching walks the run in LOGICAL batches, because the allocate/free
    // cursors can only move between batches -- no chromatogram may appear or
    // vanish while a worker is reading `live`. The batch is therefore the
    // GRANULARITY OF THE SLIDING WINDOW: a precursor goes live up to one batch
    // early and is freed up to one batch late.
    //
    // The logical batches are a fixed grid of MATCH_BATCH spectra from the
    // pass's first spectrum, independent of the decode block. They used to be
    // the block's own sub-batches, so a block that was not a multiple of 128
    // moved the batch boundaries -- 200 gave [128, 200) and [200, 328) where
    // 256 gave [128, 256) -- and with them the moment each precursor was handed
    // over, i.e. the stream's ORDER, which Gate C calibrates on. Now a block is
    // decoded at whatever size was asked for, each logical batch is matched
    // over its intersections with the decoded blocks, and activation, expiry
    // and hand-over happen only at logical-batch boundaries. 256 is a multiple
    // of 128, so the default block gives exactly the batches it always gave.
    //
    // 128 spectra is ~5 cycles across IH1's 24 windows, about 7 s of gradient
    // against retention-time windows of hundreds -- so the rounding is under 1%
    // of what is resident. At the 1,024 the decode uses it would be 59 s, which
    // is 5% of a 1,200 s window and 100% of a short one. The cost of the finer
    // batch is one more pool wakeup per 128 spectra, microseconds against the
    // ~1 s the decode of those spectra takes.
    constexpr std::size_t MATCH_BATCH = 128;

    std::size_t first_spectrum = 0, last_spectrum = info.size();
    if (options.rt_high > options.rt_low)
    {
      first_spectrum = static_cast<std::size_t>(
        std::lower_bound(info.begin(), info.end(), options.rt_low,
                         [](const SpectrumInfo& s, double v) { return s.retention_time < v; })
        - info.begin());
      last_spectrum = static_cast<std::size_t>(
        std::lower_bound(info.begin(), info.end(), options.rt_high,
                         [](const SpectrumInfo& s, double v) { return s.retention_time < v; })
        - info.begin());
      // The cycle range above was resolved on the FLOAT axis, this spectrum
      // range on the run's DOUBLE times, and the two disagree by one spectrum
      // whenever a time rounds across the bound. Widened to every cycle the
      // axis admitted, so a precursor's first or last cycle can never fall
      // outside what the pass reads.
      for (std::size_t w = 0; w < windows.size(); ++w)
      {
        if (global_lo[w] >= global_hi[w]) { continue; }
        first_spectrum = std::min(first_spectrum, axis[w].position[global_lo[w]]);
        last_spectrum = std::max(last_spectrum, axis[w].position[global_hi[w] - 1] + 1);
      }
    }
    st.spectra_read = last_spectrum - first_spectrum;

    // The batch grid is ABSOLUTE: batch boundaries sit at first_spectrum +
    // k*MATCH_BATCH whichever chunk is reading and whatever the decode block,
    // so a spectrum is in the same batch in every chunk plan and at every
    // block size. A chunk that started its grid at its own first spectrum
    // would move every boundary, and with them the moment each precursor goes
    // live and is handed over.
    const auto batchOf = [&](std::size_t si) { return (si - first_spectrum) / MATCH_BATCH; };
    const std::size_t n_batches =
      last_spectrum > first_spectrum ? batchOf(last_spectrum - 1) + 1 : 0;

    // Each precursor's span on that grid. It goes live in the batch that reads
    // its first cycle and is handed over after the batch that reads its last,
    // so it is resident for exactly the batches [first_batch, last_batch].
    // Computed from spectrum POSITIONS, never from retention times: see
    // CycleAxis::position.
    std::vector<float> slot_rt_lo(n_slots, 0.0f);
    std::vector<std::uint32_t> first_batch(n_slots, 0), last_batch(n_slots, 0);
    for (std::size_t s = 0; s < n_slots; ++s)
    {
      const Assignment& a = assignments[s];
      slot_rt_lo[s] = axis[a.window].rt[a.lo];
      const std::size_t si_lo = axis[a.window].position[a.lo];
      const std::size_t si_hi = axis[a.window].position[a.hi - 1];
      if (si_lo < first_spectrum || si_hi >= last_spectrum)
      {
        throw std::logic_error(
          "precursor " + std::to_string(a.precursor) + " spans spectra " +
          std::to_string(si_lo) + ".." + std::to_string(si_hi) +
          ", outside the pass's range " + std::to_string(first_spectrum) + ".." +
          std::to_string(last_spectrum));
      }
      first_batch[s] = static_cast<std::uint32_t>(batchOf(si_lo));
      last_batch[s] = static_cast<std::uint32_t>(batchOf(si_hi));
    }
    std::vector<std::uint32_t> by_start(n_slots);
    std::iota(by_start.begin(), by_start.end(), 0u);
    std::sort(by_start.begin(), by_start.end(),
              [&](std::uint32_t a, std::uint32_t b) { return slot_rt_lo[a] < slot_rt_lo[b]; });

    // How much of the library is live at once.
    //
    // This is the number the whole design turns on, so it is measured rather
    // than assumed -- and measured on the grid the pass actually moves on. It
    // used to be the overlap of the retention-time windows, which is what would
    // be live if the window slid per spectrum. It slides per batch, so that
    // undercounted by up to a batch's activations at each end: 1,461,670 swept
    // against 2,349,537 resident on the PXD fixture, and a chunked run peaked
    // at 734,070 live under a cap of 400,000 (worst where a batch straddles a
    // gap in the run and holds the slices on both sides live at once). The
    // live set at batch B is the precursors whose span contains B, exactly.
    std::size_t overlap_precursors = 0;
    {
      std::vector<std::int64_t> delta(n_batches + 1, 0);
      for (std::size_t s = 0; s < n_slots; ++s)
      {
        ++delta[first_batch[s]];
        --delta[std::size_t(last_batch[s]) + 1];
      }
      std::int64_t live = 0;
      for (std::size_t b = 0; b < n_batches; ++b)
      {
        live += delta[b];
        overlap_precursors = std::max(overlap_precursors, std::size_t(live));
      }
    }

    // The order the pass hands precursors to the sink, fixed ONCE for the
    // whole run and independent of how it is chunked.
    //
    // Unchunked, the pass hands over after each batch, window by window, every
    // precursor whose last cycle that batch read, in the order of its window's
    // end cursor. So the stream is sorted by (last_batch, window, rank in that
    // cursor). The cursors are built here exactly as a single-chunk pass
    // always built them -- the same start-ordered input, sorted by lo and then
    // by hi -- which keeps a run that fits byte-identical to before.
    //
    // It has to be fixed, because the sink is not order-free: Gate C's
    // threshold is the quantile of the FIRST gate_calibration_n decoys to
    // arrive, and everything arriving before that is admitted unconditionally.
    // Chunks that each emitted in their own order handed the gate a different
    // calibration sample -- tau 6.61 unchunked against 6.12 at 13 chunks on the
    // PXD fixture, and 2,634,793 peak groups against 3,083,453.
    //
    // The per-window cursors are kept: they ARE the hand-over order within a
    // window (a precursor's rank is its place in its window's cursor), so the
    // planner merges them rather than sorting the library again, and a single
    // chunk hands over straight from them.
    std::vector<std::vector<std::uint32_t>> cursor(windows.size());
    for (const std::uint32_t s : by_start) { cursor[assignments[s].window].push_back(s); }
    for (auto& c : cursor)
    {
      std::sort(c.begin(), c.end(), [&](std::uint32_t a, std::uint32_t b) {
        return assignments[a].lo < assignments[b].lo; });
      std::sort(c.begin(), c.end(), [&](std::uint32_t a, std::uint32_t b) {
        return assignments[a].hi < assignments[b].hi; });
    }

    // The library, split into chunks whose live sets each fit under the cap.
    //
    // A chunk is a CONTIGUOUS RUN OF THE HAND-OVER ORDER above, cut where
    // adding the next precursor would put more than `cap` live in some batch.
    // That makes chunking invisible to the sink by construction: each chunk
    // reads every spectrum its precursors span (positions, not times), so each
    // trace holds what the single pass would have put in it; it hands its
    // precursors over in the single pass's order; and the chunks follow one
    // another in that order, so their concatenation IS the single pass's
    // stream. The cap holds exactly, because the live count it is tested
    // against is the one the pass will have.
    //
    // A chunk's precursors END within a run of consecutive batches, but they
    // can START anywhere before that: a long-window precursor pulls the
    // chunk's first spectrum back towards the start of the run. So a chunk
    // reads from its earliest start to its latest end, which is usually a
    // fraction of a pass but is not bounded by one -- the worst case decodes
    // O(chunks x spectra). Measured on the PXD fixture: 13,348 spectra decoded
    // for 8 chunks of a 6,020-spectrum run (cap 400,000), 35,953 for 31.
    std::vector<std::vector<std::uint32_t>> chunks;
    // Derive the cap from the memory budget when one is given. Stated in
    // bytes because that is the quantity a caller actually has, and inverted
    // here because only the extractor knows the mean cells per precursor.
    // Inverted with the MEAN, so the result is exact in precursors and only
    // approximate in bytes: see Options::live_memory_budget_bytes.
    std::size_t cap = options.max_live_precursors;
    if (options.live_memory_budget_bytes > 0)
    {
      std::size_t planes = 1;
      if (options.collect_mass_residuals) { planes += 2; }
      if (options.collect_im_residuals)   { planes += 2; }
      double mean_cells = 0.0;
      for (const Assignment& a : assignments)
      { mean_cells += double(a.valid) * double(a.hi - a.lo); }
      if (!assignments.empty()) { mean_cells /= double(assignments.size()); }
      const double per = mean_cells * 4.0 * double(planes);
      // 0 when a live precursor costs nothing -- every assigned precursor has
      // zero valid transitions, which is supported input -- and then the budget
      // constrains nothing at all.
      const std::size_t derived = per > 0.0
        ? std::max<std::size_t>(1, std::size_t(double(options.live_memory_budget_bytes) / per))
        : 0;
      // The tighter of the constraints that actually BIND wins. An explicit
      // -max_live_precursors is a caller's assertion and must not be loosened
      // by a budget -- including a budget that binds nothing: min(cap, 0) used
      // to read that 0 as "unlimited" and erase the explicit cap.
      if (derived != 0) { cap = (cap == 0) ? derived : std::min(cap, derived); }
      st.live_budget_note = "budget " +
        std::to_string(options.live_memory_budget_bytes / (1024ull*1024*1024)) +
        " GiB / " + std::to_string(std::size_t(per)) + " B per live precursor (" +
        std::to_string(planes) + " planes) -> cap " + std::to_string(cap);
    }
    // Whether the chunks hand over straight from the window cursors -- only
    // when the planner did not cut the library, so the one chunk is all of it.
    bool from_cursors = false;
    if (cap == 0 || overlap_precursors <= cap)
    {
      from_cursors = true;
      chunks.push_back(std::move(by_start));
      st.memory_bound_by = "retention-time overlap (" +
                           std::to_string(overlap_precursors) + " of " +
                           std::to_string(n_slots) + " precursors live at once)";
    }
    else
    {
      // The hand-over order, (last_batch, window, rank), as a MERGE of the
      // window cursors: within a window the cursor is sorted by `hi`, and
      // last_batch is the batch of the spectrum at hi-1, so it never decreases
      // along a cursor. Walking the batches and, in each, the windows in index
      // order yields exactly the order a sort on the full key would -- the key
      // is unique, since a rank is unique within its window -- in O(N + batches
      // x windows) rather than another O(N log N) over the library.
      std::vector<std::uint32_t> order;
      order.reserve(n_slots);
      {
        std::vector<std::size_t> head(windows.size(), 0);
        for (std::size_t b = 0; b < n_batches; ++b)
        {
          for (std::size_t w = 0; w < windows.size(); ++w)
          {
            const auto& c = cursor[w];
            while (head[w] < c.size() && last_batch[c[head[w]]] == b) { order.push_back(c[head[w]++]); }
            if (head[w] < c.size() && last_batch[c[head[w]]] < b)
            {
              throw std::logic_error("hand-over cursor of window " + std::to_string(w) +
                                     " is not ordered by its last batch");
            }
          }
        }
        if (order.size() != n_slots)
        {
          throw std::logic_error("the hand-over merge placed " + std::to_string(order.size()) +
                                 " of " + std::to_string(n_slots) + " precursors");
        }
      }
      LiveCount live_count(n_batches);
      std::vector<std::uint32_t> current;
      for (const std::uint32_t s : order)
      {
        if (!current.empty() &&
            live_count.max(first_batch[s], last_batch[s]) + 1 > cap)
        {
          chunks.push_back(std::move(current));
          current.clear();
          live_count.reset();
        }
        live_count.add(first_batch[s], last_batch[s]);
        current.push_back(s);
      }
      if (!current.empty()) { chunks.push_back(std::move(current)); }
      std::vector<std::vector<std::uint32_t>>().swap(cursor);   // each chunk filters its own
      st.memory_bound_by = "precursor cap (" + std::to_string(cap) + "), " +
                           std::to_string(chunks.size()) + " chunks; " +
                           std::to_string(overlap_precursors) +
                           " would have been live at once";
    }
    st.chunks = chunks.size();
    if (chunks.empty()) { chunks.emplace_back(); st.chunks = 1; }

    st.index_seconds = std::chrono::duration<double>(
                         std::chrono::steady_clock::now() - t_index).count();

    // Where every precursor's points are right now: null outside its window,
    // which for most precursors is most of the run.
    //
    // The cycle range is set for every slot up front and KEPT after the slot
    // is handed over, so a peak that lands in a precursor's range while it is
    // not live reaches the `unhoused` test below. Resetting it to 0..0 on
    // hand-over (and leaving it 0..0 before activation) made the cycle test
    // skip such a peak first, and the guard could never fire.
    std::vector<LiveSlot> live(n_slots);
    for (std::size_t s = 0; s < n_slots; ++s)
    {
      live[s].lo = assignments[s].lo;
      live[s].hi = assignments[s].hi;
    }
    BlockPool blocks;
    std::size_t live_now = 0, live_peak = 0;
    std::atomic<std::size_t> nonzero{0};
    // Points a peak fell on while the destination precursor was not live. It
    // must be zero -- the cycle test below already excludes them -- and it is
    // counted rather than assumed, because the alternative to catching it is a
    // silently truncated chromatogram. Checked after EVERY match batch, before
    // anything that batch finished is handed over: a trace that lost a point
    // must never reach the sink, and a defect found only at the end of a
    // 30-hour extraction has already delivered every trace it damaged.
    std::atomic<std::size_t> unhoused{0};
    std::mutex unhoused_m;
    std::uint32_t unhoused_slot = std::numeric_limits<std::uint32_t>::max();
    std::uint32_t unhoused_cycle = 0;
    // The other half of that guarantee, and the half `unhoused` cannot see: a
    // cycle the pass never READ leaves no peak to count. Per window, the cycles
    // the current chunk has read are [first_seen, seen) -- a chunk reads a
    // contiguous run of spectra, so a window's cycles in it are contiguous --
    // and a precursor handed over with its range not inside that has cells that
    // are zero because nothing looked. That is how chunked extraction lost the
    // first or last frame of a chunk for as long as the chunk's spectrum range
    // was looked up by a float time. It is checked BEFORE every hand-over and
    // throws there, so an incomplete trace is never handed to the sink.
    std::vector<std::uint32_t> seen(windows.size(), 0);
    std::vector<std::uint32_t> first_seen(windows.size(),
                                          std::numeric_limits<std::uint32_t>::max());
    // Where the pass is, for the two errors above.
    std::size_t chunk_no = 0, chunk_first = first_spectrum, chunk_last = last_spectrum;
    const auto describe = [&](std::uint32_t slot) {
      const Assignment& a = assignments[slot];
      const auto seq = library.strings().get(p.modified_sequence[a.precursor]);
      return "precursor " + std::to_string(a.precursor) +
             (seq.empty() ? std::string() : " (" + std::string(seq) + "/" +
                                              std::to_string(int(p.charge[a.precursor])) + ")") +
             ", window " + std::to_string(a.window) + ", cycles [" + std::to_string(a.lo) +
             ", " + std::to_string(a.hi) + "), chunk " + std::to_string(chunk_no + 1) + " of " +
             std::to_string(chunks.size()) + " reading spectra [" + std::to_string(chunk_first) +
             ", " + std::to_string(chunk_last) + ")";
    };

    // One forward pass over the run in acquisition order, so the decode stays
    // contiguous -- a window's spectra are strided through the file, and
    // fetching them window by window would ask the reader for a stride it
    // cannot serve cheaply.
    //
    // Each spectrum is matched by INVERTING the loop: for every peak, the m/z
    // index says which transitions could contain it. The alternative -- probe
    // the spectrum once per transition -- is O(N log M) with N transitions and
    // M peaks, and at proteome scale N is 4.6 M against M of ~1,100.
    //
    // Matching runs in parallel with no synchronisation at all: a
    // (transition, cycle) pair has exactly one destination index, computed
    // rather than allocated, so no two threads ever write the same slot.
    const unsigned threads = options.threads ? options.threads
                                             : std::max(1u, std::thread::hardware_concurrency());
    std::vector<SpectrumPeaks> block;
    // Workers are created ONCE and parked between blocks. Spawning them per
    // block cost more than the matching did: at 64 workers and a 128-spectrum
    // block each got two spectra, and parallel ran slower than serial.
    // Raising the block size hid that; creating them once removes it.
    struct Pool
    {
      std::vector<std::thread> workers;
      std::mutex m;
      std::condition_variable cv, done_cv;
      std::function<void()> job;
      std::size_t generation = 0, finished = 0;
      bool stop = false;

      void start(unsigned n)
      {
        for (unsigned i = 0; i < n; ++i)
        {
          workers.emplace_back([this] {
            std::size_t seen = 0;
            for (;;)
            {
              std::function<void()> mine;
              {
                std::unique_lock<std::mutex> lock(m);
                cv.wait(lock, [&] { return stop || generation != seen; });
                if (stop) { return; }
                seen = generation;
                mine = job;
              }
              mine();
              {
                std::lock_guard<std::mutex> lock(m);
                ++finished;
              }
              done_cv.notify_one();
            }
          });
        }
      }

      void run(std::function<void()> f, unsigned n)
      {
        {
          std::lock_guard<std::mutex> lock(m);
          job = std::move(f);
          finished = 0;
          ++generation;
        }
        cv.notify_all();
        std::unique_lock<std::mutex> lock(m);
        done_cv.wait(lock, [&] { return finished >= n; });
      }

      ~Pool()
      {
        {
          std::lock_guard<std::mutex> lock(m);
          stop = true;
        }
        cv.notify_all();
        for (auto& w : workers) { if (w.joinable()) { w.join(); } }
      }
    } pool_impl;
    if (threads > 1) { pool_impl.start(threads); }
    // A precursor leaving the pass: its chromatogram is final, so hand it over
    // and give the memory back.
    std::vector<std::uint64_t> off_scratch;
    std::vector<std::uint32_t> count_scratch;
    const auto activate = [&](std::uint32_t slot) {
      const Assignment& a = assignments[slot];
      const std::size_t cells = std::size_t(a.valid) * (a.hi - a.lo);
      live[slot].base = blocks.take(cells);
      if (options.collect_mass_residuals)
      {
        live[slot].ppm_num = blocks.take(cells);
        live[slot].ppm_den = blocks.take(cells);
      }
      if (options.collect_im_residuals)
      {
        live[slot].im_num = blocks.take(cells);
        live[slot].im_den = blocks.take(cells);
      }
      live[slot].lo = a.lo;
      live[slot].hi = a.hi;
      live_peak = std::max(live_peak, ++live_now);
    };
    // -parallel_sink: a sink that asks for batches gets the released
    // precursors staged here and handed over together, in exactly the order
    // the one-at-a-time path would have called accept(). Their blocks go back
    // to the pool only after the batch returns, so the views stay valid. A
    // sink that does not ask (batchCapacity() == 0, every sink by default)
    // never reaches any of this: emit() below takes its historical path.
    const std::size_t sink_batch = sink.batchCapacity();
    struct PendingGive
    {
      float* base; float* ppm_num; float* ppm_den; float* im_num; float* im_den;
      std::size_t cells;
    };
    std::vector<PrecursorChromatogram> pending;
    std::vector<std::size_t> pending_at;        ///< start of each trace's offsets/counts
    std::vector<std::uint64_t> pending_off;
    std::vector<std::uint32_t> pending_count;
    std::vector<PendingGive> pending_give;
    const auto flush = [&] {
      if (pending.empty()) { return; }
      // Pointers fixed up only now: the flat arrays may have reallocated while
      // the batch was being staged.
      for (std::size_t q = 0; q < pending.size(); ++q)
      {
        pending[q].offset = pending_off.data() + pending_at[q];
        pending[q].count = pending_count.data() + pending_at[q];
      }
      const auto t_sink = std::chrono::steady_clock::now();
      sink.acceptBatch(pending.data(), pending.size());
      st.sink_seconds += std::chrono::duration<double>(
                           std::chrono::steady_clock::now() - t_sink).count();
      for (const PendingGive& g : pending_give)
      {
        blocks.give(g.base, g.cells);
        if (g.ppm_num != nullptr) { blocks.give(g.ppm_num, g.cells); blocks.give(g.ppm_den, g.cells); }
        if (g.im_num != nullptr) { blocks.give(g.im_num, g.cells); blocks.give(g.im_den, g.cells); }
      }
      pending.clear();
      pending_at.clear();
      pending_off.clear();
      pending_count.clear();
      pending_give.clear();
    };
    const auto emit = [&](std::uint32_t slot) {
      const Assignment& a = assignments[slot];
      // Every cycle the precursor spans must have been read by this chunk, and
      // that is checked HERE, before the sink sees anything.
      if (a.lo < first_seen[a.window] || a.hi > seen[a.window])
      {
        const std::uint32_t read_lo = first_seen[a.window], read_hi = seen[a.window];
        std::string missing;
        const auto range = [&](std::uint32_t lo, std::uint32_t hi) {
          if (lo >= hi) { return; }
          if (!missing.empty()) { missing += " and "; }
          missing += "[" + std::to_string(lo) + ", " + std::to_string(hi) + ")";
        };
        if (read_lo >= read_hi) { range(a.lo, a.hi); }
        else
        {
          range(a.lo, std::min(a.hi, read_lo));
          range(std::max(a.lo, read_hi), a.hi);
        }
        throw std::logic_error(
          describe(slot) + " was about to be handed over with cycles " + missing +
          " never read (the chunk read cycles [" +
          (read_lo >= read_hi ? std::string("none") :
           std::to_string(read_lo) + ", " + std::to_string(read_hi)) +
          ") of that window). Its cells there would be zero because nothing looked, "
          "not because nothing was there -- the chunk's spectrum range did not cover "
          "its precursors.");
      }
      // Unreachable once the test above has passed: a precursor whose first
      // cycle the chunk read was activated by the batch that read it.
      if (live[slot].base == nullptr)
      {
        throw std::logic_error(describe(slot) +
                               " was about to be handed over without ever going live");
      }

      const std::uint32_t tb = p.transition_begin[a.precursor];
      const std::uint32_t tc = p.transition_count[a.precursor];
      const std::uint32_t cycles = a.hi - a.lo;
      off_scratch.assign(tc, 0);
      count_scratch.assign(tc, 0);
      std::uint32_t row = 0;
      for (std::uint32_t k = 0; k < tc; ++k)
      {
        // Skipped on the same test as the counting loop above, and it has to
        // BE the same test: a transition counted there and skipped here is a
        // point run with no axis behind it.
        if (t.product_mz[tb + k] == MZ_INVALID) { continue; }
        off_scratch[k] = std::uint64_t(row) * cycles;
        count_scratch[k] = cycles;
        ++row;
      }

      // The representation's invariant, checked rather than assumed: every
      // point of the precursor must land on the axis it names. Breaking it does
      // not crash -- a consumer reads past the end of a vector and returns
      // whatever is there -- so it has to be caught here or not at all.
      if (a.window >= axes.size() ||
          std::size_t(a.lo) + cycles > axes[a.window].size())
      {
        throw std::logic_error(
          "chromatogram precursor " + std::to_string(a.precursor) + " has " +
          std::to_string(cycles) + " points from cycle " + std::to_string(a.lo) +
          " on axis " + std::to_string(a.window) + ", which holds " +
          std::to_string(a.window < axes.size() ? axes[a.window].size() : 0) + " cycles");
      }

      PrecursorChromatogram trace;
      trace.precursor = a.precursor;
      trace.transition_begin = tb;
      trace.transition_count = tc;
      trace.axis = a.window;
      trace.axis_begin = a.lo;
      trace.cycles = cycles;
      trace.rt = axes[a.window].data() + a.lo;
      trace.points = live[slot].base;
      trace.ppm_num = live[slot].ppm_num;
      trace.ppm_den = live[slot].ppm_den;
      trace.im_num = live[slot].im_num;
      trace.im_den = live[slot].im_den;
      trace.offset = off_scratch.data();
      trace.count = count_scratch.data();

      if (sink_batch != 0)
      {
        pending_at.push_back(pending_off.size());
        pending_off.insert(pending_off.end(), off_scratch.begin(), off_scratch.end());
        pending_count.insert(pending_count.end(), count_scratch.begin(), count_scratch.end());
        pending.push_back(trace);
        pending_give.push_back({live[slot].base, live[slot].ppm_num, live[slot].ppm_den,
                                live[slot].im_num, live[slot].im_den,
                                std::size_t(a.valid) * cycles});
        // The slot's bookkeeping is released now, exactly as below -- including
        // the kept cycle range, so `unhoused` still sees a late peak; only the
        // memory waits for the batch.
        live[slot] = LiveSlot{};
        live[slot].lo = a.lo;
        live[slot].hi = a.hi;
        --live_now;
        if (pending.size() >= sink_batch) { flush(); }
        return;
      }

      const auto t_sink = std::chrono::steady_clock::now();
      sink.accept(trace);
      st.sink_seconds += std::chrono::duration<double>(
                           std::chrono::steady_clock::now() - t_sink).count();

      blocks.give(live[slot].base, std::size_t(a.valid) * cycles);
      // The residual planes are pool blocks too, and resetting the slot without
      // giving them back leaks two per emitted precursor. Mild on a 2,665-
      // precursor benchmark (3.55 -> 3.87 GiB) because the pool retains them;
      // at the 4.26 M design scale it would be the whole live set again, twice.
      if (live[slot].ppm_num != nullptr)
      {
        blocks.give(live[slot].ppm_num, std::size_t(a.valid) * cycles);
        blocks.give(live[slot].ppm_den, std::size_t(a.valid) * cycles);
      }
      if (live[slot].im_num != nullptr)
      {
        blocks.give(live[slot].im_num, std::size_t(a.valid) * cycles);
        blocks.give(live[slot].im_den, std::size_t(a.valid) * cycles);
      }
      live[slot] = LiveSlot{};
      live[slot].lo = a.lo;     // kept, so `unhoused` can see a late peak
      live[slot].hi = a.hi;
      --live_now;
    };

    for (chunk_no = 0; chunk_no < chunks.size(); ++chunk_no)
    {
      const auto& slots = chunks[chunk_no];
      const auto t_chunk_index = std::chrono::steady_clock::now();
      std::fill(seen.begin(), seen.end(), 0u);
      std::fill(first_seen.begin(), first_seen.end(),
                std::numeric_limits<std::uint32_t>::max());

      // The m/z index of this chunk's transitions, per window. Rebuilt per
      // chunk rather than held for the whole library, because the index is
      // 18 bytes a transition and a chunk exists precisely because the whole
      // library did not fit.
      std::vector<MzIndex> index(windows.size());
      {
        std::vector<std::size_t> entries(windows.size(), 0);
        for (const std::uint32_t slot : slots)
        {
          entries[assignments[slot].window] += assignments[slot].valid;
        }
        for (std::size_t w = 0; w < windows.size(); ++w)
        {
          index[w].mz.reserve(entries[w]);
          index[w].slot.reserve(entries[w]);
          index[w].row.reserve(entries[w]);
          index[w].precursor_im.reserve(entries[w]);
        }
      }
      for (const std::uint32_t slot : slots)
      {
        const Assignment& a = assignments[slot];
        auto& x = index[a.window];
        std::uint16_t row = 0;
        for (std::uint32_t k = 0; k < p.transition_count[a.precursor]; ++k)
        {
          const std::size_t j = p.transition_begin[a.precursor] + k;
          if (t.product_mz[j] == MZ_INVALID) { continue; }
          const double theoretical = fromFixed(t.product_mz[j]);
          // The offset is applied to the TRANSITION, once, here -- not to every
          // peak in the match loop. Shifting the target is equivalent and costs
          // nothing per peak. The m/z-dependent term rides along for free for
          // the same reason: it is a property of the transition's own mass, so
          // it is known here and never has to be evaluated per peak.
          double correction_ppm = options.fragment_ppm_offset;
          if (options.fragment_ppm_log_slope != 0.0 && options.fragment_ppm_ref_mz > 0.0)
          {
            correction_ppm +=
              options.fragment_ppm_log_slope * std::log(theoretical / options.fragment_ppm_ref_mz);
          }
          if (options.fragment_ppm_slope_per_1000 != 0.0)
          {
            correction_ppm += options.fragment_ppm_slope_per_1000 *
                              (theoretical - options.fragment_ppm_ref_mz) / 1000.0;
          }
          x.mz.push_back(theoretical * (1.0 + correction_ppm * 1e-6));
          x.slot.push_back(slot);
          x.row.push_back(row);
          // The precursor's own 1/K0, carried onto each of its transitions so
          // the match loop can test a peak against the peptide it claims to be
          // from rather than against the whole frame. NaN when the library has
          // none, which disables the test for that transition.
          //
          // The run's mobility recalibration is applied HERE, once per
          // transition, for the same reason the mass offset is: it is a property
          // of the precursor, so it is known at index-build time and never has
          // to be evaluated per peak. offsetFor() returns 0 unless a model was
          // fitted AND this charge was supported AND -- for a precursor that was
          // one of the fit's own anchors -- from the fold that excluded it.
          const float lib_im = p.im[a.precursor];
          const double im_offset =
            options.mobility_model
              ? options.mobility_model->offsetFor(a.precursor, fromFixed(p.mz[a.precursor]),
                                                  p.charge[a.precursor],
                                                  p.im[a.precursor])
              : 0.0;
          x.precursor_im.push_back(static_cast<float>(lib_im + im_offset));
          ++row;
        }
      }

      // Sort each window's entries by m/z and build the bucket table. The
      // transitions of a window are fixed for the whole chunk, so this is paid
      // once rather than once per spectrum.
      for (std::size_t w = 0; w < windows.size(); ++w)
      {
        auto& x = index[w];
        std::vector<std::uint32_t> order(x.mz.size());
        for (std::size_t i = 0; i < order.size(); ++i) { order[i] = static_cast<std::uint32_t>(i); }
        std::sort(order.begin(), order.end(),
                  [&](std::uint32_t a, std::uint32_t b) { return x.mz[a] < x.mz[b]; });
        const auto permute = [&](auto& v) {
          std::decay_t<decltype(v)> tmp(v.size());
          for (std::size_t i = 0; i < order.size(); ++i) { tmp[i] = v[order[i]]; }
          v.swap(tmp);
        };
        permute(x.mz); permute(x.slot); permute(x.row); permute(x.precursor_im);
        x.build(options.fragment_ppm);
      }

      // Activation and expiry cursors: this chunk's precursors ordered by the
      // cycle their window starts at, and by the cycle it ends at.
      std::vector<std::vector<std::uint32_t>> by_lo(windows.size()), by_hi(windows.size());
      for (const std::uint32_t slot : slots) { by_lo[assignments[slot].window].push_back(slot); }
      for (std::size_t w = 0; w < windows.size(); ++w)
      {
        std::sort(by_lo[w].begin(), by_lo[w].end(), [&](std::uint32_t a, std::uint32_t b) {
          return assignments[a].lo < assignments[b].lo; });
      }
      // The hand-over cursors run in the run-wide order fixed above, not in an
      // order sorted afresh per chunk. An uncut library takes the window
      // cursors themselves. A planned chunk is a contiguous run of their merge, so
      // its precursors of one window are a contiguous run of that window's
      // cursor, already in order: filtering the chunk by window is enough.
      if (from_cursors) { by_hi = std::move(cursor); }
      else
      {
        for (const std::uint32_t slot : slots) { by_hi[assignments[slot].window].push_back(slot); }
      }
      std::vector<std::size_t> cur_lo(windows.size(), 0), cur_hi(windows.size(), 0);

      // One chunk needs only the spectra its own precursors elute in: from the
      // first cycle of any of them to the last, as spectrum POSITIONS. This
      // was a search of the run's double times for the float times of those
      // cycles, which lands one spectrum inside the range whenever the float
      // rounds the wrong way -- and that spectrum, the chunk's first or last
      // frame, was then never read. With a single chunk it is the pass's own
      // range [first_spectrum, last_spectrum), as it always was.
      chunk_first = first_spectrum;
      chunk_last = last_spectrum;
      if (chunks.size() > 1 && !slots.empty())
      {
        chunk_first = last_spectrum;
        chunk_last = first_spectrum;
        for (const std::uint32_t slot : slots)
        {
          const Assignment& a = assignments[slot];
          chunk_first = std::min(chunk_first, axis[a.window].position[a.lo]);
          chunk_last = std::max(chunk_last, axis[a.window].position[a.hi - 1] + 1);
        }
      }

      st.index_seconds += std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - t_chunk_index).count();

      // Logical batches on the run's ABSOLUTE grid (see batchOf), clipped to
      // the chunk. Decoding is per block, because that is what the reader
      // wants: blocks sit on their own absolute grid, first_spectrum +
      // j*BLOCK, clipped to the chunk (a chunk starting mid-block reads the
      // rest of that block), and a block is decoded when the pass first needs
      // one of its spectra. A logical batch that straddles two blocks is
      // matched in two parts; nothing is activated or handed over between
      // them, so the block size is invisible to everything but memory.
      std::size_t block_begin = chunk_first, block_end = chunk_first;  // held in `block`
      for (std::size_t batch_at =
             first_spectrum + ((chunk_first - first_spectrum) / MATCH_BATCH) * MATCH_BATCH;
           batch_at < chunk_last; batch_at += MATCH_BATCH)
      {
        const std::size_t batch = std::max(batch_at, chunk_first);
        const std::size_t batch_end = std::min(batch_at + MATCH_BATCH, chunk_last);
        if (batch >= batch_end) { continue; }

        // Everything whose window starts inside this batch goes live. Doing it
        // between batches rather than per spectrum is what keeps the pass free
        // of synchronisation: no allocation happens while a worker is running.
        const auto t_alloc = std::chrono::steady_clock::now();
        for (std::size_t si = batch; si < batch_end; ++si)
        {
          const std::uint32_t w = window_of[si];
          if (w == std::numeric_limits<std::uint32_t>::max()) { continue; }
          seen[w] = std::max(seen[w], cycle_of[si] + 1);
          first_seen[w] = std::min(first_seen[w], cycle_of[si]);
        }
        for (std::size_t w = 0; w < windows.size(); ++w)
        {
          while (cur_lo[w] < by_lo[w].size() &&
                 assignments[by_lo[w][cur_lo[w]]].lo < seen[w])
          {
            activate(by_lo[w][cur_lo[w]]);
            ++cur_lo[w];
          }
        }
        st.assemble_seconds += std::chrono::duration<double>(
                                 std::chrono::steady_clock::now() - t_alloc).count();

        // The workers match the part [next, part_end) of this batch that the
        // decoded block holds.
        std::atomic<std::size_t> next{batch};
        std::size_t part_end = batch;
        const auto work = [&]() {
          std::size_t local_nonzero = 0, local_unhoused = 0;
          std::uint32_t first_unhoused_slot = 0, first_unhoused_cycle = 0;
          for (;;)
          {
            const std::size_t si = next.fetch_add(1);
            if (si >= part_end) { break; }
            const std::uint32_t w = window_of[si];
            if (w == std::numeric_limits<std::uint32_t>::max()) { continue; }
            const auto& x = index[w];
            if (x.mz.empty()) { continue; }
            const auto& peaks = block[si - block_begin];
            const std::uint32_t c = cycle_of[si];
            // The run's mobility and the FRAME BAND test are separate questions.
            // Reading the peak's mobility only inside the band's branch made
            // `use_ion_mobility` switch off the per-precursor test as well: the
            // mobility stayed NaN, and a NaN mobility skips that test under the
            // "absent information is not evidence of mismatch" rule meant for a
            // run that carries no mobility at all. -no_ion_mobility is the
            // control arm for measuring what mobility filtering buys, so it has
            // to turn off exactly the one filter it names.
            const bool has_im = peaks.hasIonMobility();
            const bool use_band = options.use_ion_mobility && has_im;
            const double im_low = info[si].window.im_low, im_high = info[si].window.im_high;

            for (std::size_t k = 0; k < peaks.size(); ++k)
            {
              const double m = peaks.mz[k];
              // The bounds must allow for the tolerance. A peak just BELOW the
              // smallest transition can still be within ppm of it, and skipping
              // it silently returned zero for the lowest-m/z transition of every
              // window -- which a brute-force scan caught and nothing else would
              // have, because a chromatogram of zeros looks like an ion that is
              // simply not there.
              const double slack = m * options.fragment_ppm * 1e-6 * 1.01 + 1e-6;
              if (m + slack < x.mz.front() || m - slack > x.mz.back()) { continue; }
              // The frame's band separates the co-packed windows. It cannot
              // separate a precursor from its same-window neighbours, which is
              // what the per-transition test below does.
              const double peak_im = has_im ? double(peaks.ion_mobility[k])
                                            : std::numeric_limits<double>::quiet_NaN();
              // Half-open, [im_low, im_high). Co-packed windows are separated by
              // a DERIVED band -- the split is the midpoint between their two
              // mobility positions -- so adjacent bands share their boundary
              // exactly, and a peak sitting on it would enter both windows and be
              // integrated twice under Sum. Measure zero on real data, and the
              // reason it is stated as a convention rather than left to chance.
              if (use_band && (peak_im < im_low || peak_im >= im_high)) { continue; }
              // The tolerance belongs to the TRANSITION, not to the peak:
              // a match means |peak - transition| <= transition * ppm. Searching
              // by peak makes it tempting to size the window on the peak
              // instead, which differs at the boundary and silently includes or
              // drops edge matches -- 11 of 240 transitions disagreed with a
              // brute-force scan because of exactly that.
              //
              // So the SEARCH window is deliberately a little wide, and the
              // exact test is applied per candidate inside it.
              const double ppm = options.fragment_ppm * 1e-6;
              const double im_half = options.precursor_im_window;
              const bool sum_peaks = options.aggregate == Options::Aggregate::Sum;
              std::size_t i = x.bucket[x.bucketOf(std::max(m - slack, x.mz.front()))];
              const float intensity = peaks.intensity[k];
              for (; i < x.mz.size() && x.mz[i] <= m + slack; ++i)
              {
                if (std::abs(m - x.mz[i]) > x.mz[i] * ppm) { continue; }
                // Per-precursor mobility. Skipped when either side is unknown:
                // absent information is not evidence of mismatch.
                if (im_half > 0.0 && !std::isnan(peak_im))
                {
                  const float want = x.precursor_im[i];
                  if (!std::isnan(want) && std::abs(peak_im - want) > im_half) { continue; }
                }
                const LiveSlot& s = live[x.slot[i]];
                if (c < s.lo || c >= s.hi) { continue; }
                if (s.base == nullptr)
                {
                  if (local_unhoused++ == 0)
                  {
                    first_unhoused_slot = x.slot[i];
                    first_unhoused_cycle = c;
                  }
                  continue;
                }
                float& at = s.base[std::size_t(x.row[i]) * (s.hi - s.lo) + (c - s.lo)];
                // Maximum, not sum: two peaks inside one tolerance are the same
                // ion split by centroiding far more often than they are two ions.
                if (at == 0.0f && intensity > 0.0f) { ++local_nonzero; }
                if (sum_peaks) { at += intensity; }
                else if (intensity > at) { at = intensity; }
                // The deviation was already computed to test the match above;
                // it has been discarded here since the extractor was written.
                // Accumulated for EVERY contributing peak, not just a winner:
                // under Sum aggregation there is no winner, and picking one
                // by arrival order is picking at random.
                if (s.ppm_num != nullptr && intensity > 0.0f)
                {
                  const std::size_t at_i =
                    std::size_t(x.row[i]) * (s.hi - s.lo) + (c - s.lo);
                  s.ppm_num[at_i] +=
                    intensity * static_cast<float>((m - x.mz[i]) / x.mz[i] * 1e6);
                  s.ppm_den[at_i] += intensity;
                }
                // The OBSERVED 1/K0 of whatever produced this peak. NaN
                // mobility is skipped rather than accumulated as zero: a
                // missing measurement is not a mobility of nothing.
                if (s.im_num != nullptr && intensity > 0.0f && !std::isnan(peak_im))
                {
                  const std::size_t at_i =
                    std::size_t(x.row[i]) * (s.hi - s.lo) + (c - s.lo);
                  s.im_num[at_i] += intensity * static_cast<float>(peak_im);
                  s.im_den[at_i] += intensity;
                }
              }
            }
          }
          nonzero.fetch_add(local_nonzero);
          if (local_unhoused != 0)
          {
            unhoused.fetch_add(local_unhoused);
            std::lock_guard<std::mutex> lock(unhoused_m);
            if (unhoused_slot == std::numeric_limits<std::uint32_t>::max())
            {
              unhoused_slot = first_unhoused_slot;
              unhoused_cycle = first_unhoused_cycle;
            }
          }
        };
        for (std::size_t part = batch; part < batch_end; part = part_end)
        {
          if (part >= block_end)
          {
            const std::size_t block_at =
              first_spectrum + ((part - first_spectrum) / BLOCK) * BLOCK;
            block_begin = std::max(block_at, chunk_first);
            block_end = std::min(block_at + BLOCK, chunk_last);
            const auto t_decode = std::chrono::steady_clock::now();
            source.peaks(block_begin, block_end, block);
            st.decode_seconds += std::chrono::duration<double>(
                                   std::chrono::steady_clock::now() - t_decode).count();
            st.spectra_decoded += block_end - block_begin;
            if (options.progress_every && ((block_at - first_spectrum) / BLOCK) % 16 == 0)
            {
              std::cerr << "\r  " << st.spectra_decoded << " / " << st.spectra_read
                        << " spectra"
                        << std::flush;
            }
          }
          part_end = std::min(batch_end, block_end);
          next.store(part);
          const auto t_match = std::chrono::steady_clock::now();
          if (threads <= 1) { work(); }
          else { pool_impl.run(work, threads); }
          st.match_seconds += std::chrono::duration<double>(
                                std::chrono::steady_clock::now() - t_match).count();
        }

        // Checked before this batch hands anything over: see `unhoused`.
        if (unhoused.load() != 0)
        {
          throw std::logic_error(
            std::to_string(unhoused.load()) + " matched points in spectra [" +
            std::to_string(batch) + ", " + std::to_string(batch_end) +
            ") had no live chromatogram to go into, the first on cycle " +
            std::to_string(unhoused_cycle) + " of " + describe(unhoused_slot) +
            ". The sliding window released or had not yet allocated a precursor "
            "inside its retention-time range, which silently truncates its trace.");
        }

        // Everything the pass has now passed the end of is FINAL. This is the
        // whole change: the chromatogram goes to the consumer and the memory
        // goes back, so what is resident is what is live at one retention time
        // rather than what the library contains.
        const auto t_free = std::chrono::steady_clock::now();
        const double sink_before = st.sink_seconds;
        for (std::size_t w = 0; w < windows.size(); ++w)
        {
          while (cur_hi[w] < by_hi[w].size() &&
                 assignments[by_hi[w][cur_hi[w]]].hi <= seen[w])
          {
            emit(by_hi[w][cur_hi[w]]);
            ++cur_hi[w];
          }
        }
        flush();   // no-op unless the sink batches; a barrier when it does
        // The sink's own time is reported separately, so it is taken out here
        // rather than counted twice.
        st.assemble_seconds += std::chrono::duration<double>(
                                 std::chrono::steady_clock::now() - t_free).count()
                               - (st.sink_seconds - sink_before);
      }

      // Whatever the pass ended inside is final too.
      const auto t_flush = std::chrono::steady_clock::now();
      const double sink_before_flush = st.sink_seconds;
      for (std::size_t w = 0; w < windows.size(); ++w)
      {
        while (cur_hi[w] < by_hi[w].size())
        {
          emit(by_hi[w][cur_hi[w]]);
          ++cur_hi[w];
        }
      }
      flush();
      st.assemble_seconds += std::chrono::duration<double>(
                               std::chrono::steady_clock::now() - t_flush).count()
                             - (st.sink_seconds - sink_before_flush);
    }
    if (options.progress_every) { std::cerr << "\r" << std::string(48, ' ') << "\r"; }

    // Precursors nothing extracted still owe the sink a trace, empty, so that
    // a consumer counting them sees what it saw when every precursor had a row.
    const auto t_empty = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < n_prec; ++i)
    {
      if (p.transition_count[i] == 0 || assigned[i]) { continue; }
      PrecursorChromatogram trace;
      trace.precursor = static_cast<std::uint32_t>(i);
      trace.transition_begin = p.transition_begin[i];
      trace.transition_count = p.transition_count[i];
      sink.accept(trace);
    }
    st.sink_seconds += std::chrono::duration<double>(
                         std::chrono::steady_clock::now() - t_empty).count();

    st.precursors = n_prec;
    st.transitions = n_trans;
    st.points = total_points;
    st.nonzero_points = nonzero.load();
    st.mean_live_transitions = n_trans ? double(live_sum) / double(n_trans) : 0.0;
    // MEASURED, from the pool, not from the interval sweep: activation moves
    // only between match batches, so a precursor goes live up to one batch
    // early, and a chunked run never reaches the sweep's number at all.
    st.peak_live_precursors = live_peak;
    st.peak_live_points = blocks.peakPoints();
  }

} // namespace ODIA
