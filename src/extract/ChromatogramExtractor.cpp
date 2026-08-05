// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/ChromatogramExtractor.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
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
    return axes_bytes + vec(axis_of) + vec(axis_begin) + vec(begin) + vec(count) +
           vec(intensity);
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
      std::vector<std::uint32_t> transition;
      std::vector<std::uint32_t> first_live_cycle, last_live_cycle;
      std::vector<std::uint32_t> point_begin;   ///< where this transition's points start
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
  } // namespace

  Chromatograms ChromatogramExtractor::extract(const Library& library,
                                               SpectrumSource& source,
                                               const Options& options, Stats* stats)
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

    Chromatograms out;
    const std::size_t n_trans =
      n_prec == 0 ? 0 : p.transition_begin[n_prec - 1] + p.transition_count[n_prec - 1];
    out.begin.assign(n_trans, 0);
    out.count.assign(n_trans, 0);

    // Which window each precursor belongs to, and over which cycles.
    std::vector<MzIndex> index(windows.size());
    const bool restrict_rt = options.irt_slope != 0.0;
    std::uint64_t total_points = 0;
    std::uint64_t live_sum = 0;

    // Which window carries each precursor, and over which cycles.
    //
    // ONE window, even though several may cover it: the schemes on 12_80 and
    // S08 overlap adjacent windows by 1.0 Th, which puts 1.1% of precursors in
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
    struct Assignment
    {
      std::uint32_t precursor, window, lo, hi;
    };
    std::vector<Assignment> assignments;
    assignments.reserve(n_prec);

    for (std::size_t i = 0; i < n_prec; ++i)
    {
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
      if (covering == 0) { ++out.precursors_without_window; continue; }
      if (covering > 1) { ++out.precursors_in_several_windows; }

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

      assignments.push_back({static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(w),
                             static_cast<std::uint32_t>(lo), static_cast<std::uint32_t>(hi)});
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
        out.count[j] += static_cast<std::uint32_t>(hi - lo);
        total_points += hi - lo;
        live_sum += hi - lo;
      }
    }

    if (total_points > std::numeric_limits<std::uint32_t>::max())
    {
      throw std::runtime_error(
        "more than 2^32 chromatogram points; narrow the retention-time window "
        "(-rt_window) or reduce the precursor count. The CSR index is 32-bit.");
    }

    std::uint64_t running = 0;
    for (std::size_t j = 0; j < n_trans; ++j)
    {
      out.begin[j] = static_cast<std::uint32_t>(running);
      running += out.count[j];
    }
    out.intensity.assign(running, 0.0f);
    out.axis_of.assign(n_trans, 0);
    out.axis_begin.assign(n_trans, 0);

    // One axis per window, shared by every transition extracted from it. The
    // per-point timestamps this replaces were the same ~1,000 values repeated
    // once per transition.
    out.axes.resize(axis.size());
    for (std::size_t w = 0; w < axis.size(); ++w) { out.axes[w] = axis[w].rt; }

    std::vector<std::uint32_t> offset(n_trans, 0);
    for (const auto& a : assignments)
    {
      auto& x = index[a.window];
      for (std::uint32_t k = 0; k < p.transition_count[a.precursor]; ++k)
      {
        const std::size_t j = p.transition_begin[a.precursor] + k;
        // Skipped on the same test as the counting loop above, and it has to
        // BE the same test: a transition counted there and skipped here is a
        // point run with no axis behind it.
        if (t.product_mz[j] == MZ_INVALID) { continue; }
        const double theoretical = fromFixed(t.product_mz[j]);
        // The offset is applied to the TRANSITION, once, here -- not to every
        // peak in the match loop. Shifting the target is equivalent and costs
        // nothing per peak. The m/z-dependent term rides along for free for the
        // same reason: it is a property of the transition's own mass, so it is
        // known here and never has to be evaluated per peak.
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
        const double product = theoretical * (1.0 + correction_ppm * 1e-6);
        x.mz.push_back(product);
        x.transition.push_back(static_cast<std::uint32_t>(j));
        x.first_live_cycle.push_back(a.lo);
        x.last_live_cycle.push_back(a.hi);
        x.point_begin.push_back(out.begin[j] + offset[j]);
        // The precursor's own 1/K0, carried onto each of its transitions so the
        // match loop can test a peak against the peptide it claims to be from
        // rather than against the whole frame. NaN when the library has none,
        // which disables the test for that transition.
        x.precursor_im.push_back(p.im[a.precursor]);
        offset[j] += a.hi - a.lo;
      }
    }

    // Sort each window's entries by m/z and build the bucket table. The
    // transitions of a window are fixed for the whole run, so this is paid once
    // rather than once per spectrum.
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
      permute(x.mz); permute(x.transition); permute(x.precursor_im);
      permute(x.first_live_cycle); permute(x.last_live_cycle);
      permute(x.point_begin);
      x.build(options.fragment_ppm);

      // A point's time is fixed by its cycle, so nothing is written per point.
      // Each transition records its window's axis and the cycle its run starts
      // at, and retentionTime() resolves the rest -- one pair of uint32 per
      // transition in place of a float per point. Writing a time during the
      // pass instead meant walking every transition of a window for every
      // spectrum, the exact O(N) per spectrum this index exists to remove.
      //
      // Assigning (rather than accumulating) is correct only because a
      // transition appears in exactly ONE window's index; the loop over
      // assignments guarantees that and the invariant below re-checks it.
      for (std::size_t i = 0; i < x.transition.size(); ++i)
      {
        const std::uint32_t tr = x.transition[i];
        out.axis_of[tr] = static_cast<std::uint32_t>(w);
        out.axis_begin[tr] = x.first_live_cycle[i];
      }
    }

    // The representation's invariant, checked rather than assumed: every point
    // of every transition must land on the axis the transition names. Breaking
    // it does not crash -- retentionTime() reads past the end of a vector and
    // returns whatever is there -- so it has to be caught here or not at all.
    for (std::size_t j = 0; j < n_trans; ++j)
    {
      if (out.count[j] == 0) { continue; }
      const std::size_t a = out.axis_of[j];
      if (a >= out.axes.size() ||
          std::size_t(out.axis_begin[j]) + out.count[j] > out.axes[a].size())
      {
        throw std::logic_error(
          "chromatogram transition " + std::to_string(j) + " has " +
          std::to_string(out.count[j]) + " points from cycle " +
          std::to_string(out.axis_begin[j]) + " on axis " + std::to_string(a) +
          ", which holds " + std::to_string(a < out.axes.size() ? out.axes[a].size() : 0) +
          " cycles");
      }
    }
    st.index_seconds = std::chrono::duration<double>(
                         std::chrono::steady_clock::now() - t_index).count();

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
    // Large enough that spawning workers is amortised. At 128 spectra and 64
    // threads each worker got two spectra and thread creation cost more than
    // the matching did -- measured 0.19 s against 0.14 s single-threaded.
    constexpr std::size_t BLOCK = 1024;
    std::atomic<std::size_t> nonzero{0};

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
    }

    for (std::size_t begin = first_spectrum; begin < last_spectrum; begin += BLOCK)
    {
      const std::size_t end = std::min(begin + BLOCK, last_spectrum);
      const auto t_decode = std::chrono::steady_clock::now();
      source.peaks(begin, end, block);
      st.decode_seconds += std::chrono::duration<double>(
                             std::chrono::steady_clock::now() - t_decode).count();
      st.spectra_read += end - begin;

      const auto t_match = std::chrono::steady_clock::now();
      std::atomic<std::size_t> next{begin};
      const auto work = [&]() {
        std::size_t local_nonzero = 0;
        for (;;)
        {
          const std::size_t si = next.fetch_add(1);
          if (si >= end) { break; }
          const std::uint32_t w = window_of[si];
          if (w == std::numeric_limits<std::uint32_t>::max()) { continue; }
          const auto& x = index[w];
          if (x.mz.empty()) { continue; }
          const auto& peaks = block[si - begin];
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
            if (use_band && (peak_im < im_low || peak_im > im_high)) { continue; }
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
              if (c < x.first_live_cycle[i] || c >= x.last_live_cycle[i]) { continue; }
              const std::size_t at = x.point_begin[i] + (c - x.first_live_cycle[i]);
              // Maximum, not sum: two peaks inside one tolerance are the same
              // ion split by centroiding far more often than they are two ions.
              if (out.intensity[at] == 0.0f && intensity > 0.0f) { ++local_nonzero; }
              if (sum_peaks) { out.intensity[at] += intensity; }
              else if (intensity > out.intensity[at]) { out.intensity[at] = intensity; }
            }
          }
        }
        nonzero.fetch_add(local_nonzero);
      };
      if (threads <= 1) { work(); }
      else { pool_impl.run(work, threads); }
      st.match_seconds += std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - t_match).count();

      if (options.progress_every && (begin / BLOCK) % 16 == 0)
      {
        std::cerr << "\r  " << st.spectra_read << " / " << (last_spectrum - first_spectrum)
                  << " spectra"
                  << std::flush;
      }
    }
    if (options.progress_every) { std::cerr << "\r" << std::string(48, ' ') << "\r"; }

    st.precursors = n_prec;
    st.transitions = n_trans;
    st.points = out.points();
    st.nonzero_points = nonzero.load();
    st.mean_live_transitions = n_trans ? double(live_sum) / double(n_trans) : 0.0;
    return out;
  }

} // namespace ODIA
