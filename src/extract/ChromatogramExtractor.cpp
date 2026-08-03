// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/ChromatogramExtractor.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace ODIA
{

  std::size_t Chromatograms::footprintBytes() const
  {
    const auto vec = [](const auto& v) { return v.capacity() * sizeof(v[0]); };
    return vec(retention_time) + vec(intensity) + vec(begin) + vec(count);
  }

  namespace
  {
    /// Largest intensity within tolerance of @p target, or 0.
    ///
    /// The maximum rather than the sum: two peaks inside a 20 ppm window are
    /// the same ion split by centroiding far more often than they are two
    /// ions, and summing them double-counts.
    float probe(const SpectrumPeaks& peaks, double target, double tolerance,
                double im_low, double im_high, bool use_im)
    {
      const double lo = target - tolerance;
      const double hi = target + tolerance;
      auto it = std::lower_bound(peaks.mz.begin(), peaks.mz.end(), lo);
      float best = 0.0f;
      for (; it != peaks.mz.end() && *it <= hi; ++it)
      {
        const std::size_t k = static_cast<std::size_t>(it - peaks.mz.begin());
        if (use_im && peaks.hasIonMobility())
        {
          const double im = peaks.ion_mobility[k];
          if (im < im_low || im > im_high) { continue; }
        }
        best = std::max(best, peaks.intensity[k]);
      }
      return best;
    }

    bool sameWindow(const IsolationWindow& a, const IsolationWindow& b)
    {
      return std::abs(a.mz_low - b.mz_low) < 1e-6 && std::abs(a.mz_high - b.mz_high) < 1e-6;
    }
  } // namespace

  Chromatograms ChromatogramExtractor::extract(const Library& library,
                                               SpectrumSource& source,
                                               const Options& options, Stats* stats)
  {
    Stats local;
    Stats& st = stats == nullptr ? local : *stats;
    st = Stats{};

    const auto& p = library.precursors();
    const auto& t = library.transitions();
    const auto& windows = source.windows();
    const auto& info = source.spectra();

    const std::size_t n_prec = options.max_precursors == 0
                                 ? library.precursorCount()
                                 : std::min(options.max_precursors, library.precursorCount());

    const auto in_range = [&](const SpectrumInfo& s) {
      return !(options.rt_high > options.rt_low &&
               (s.retention_time < options.rt_low || s.retention_time > options.rt_high));
    };

    // Which precursors each window carries. Windows overlap in most schemes,
    // so a precursor legitimately belongs to several and is extracted from
    // each -- those are different measurements of the same ion, not duplicates.
    Chromatograms out;
    std::vector<std::vector<std::uint32_t>> by_window(windows.size());
    for (std::size_t i = 0; i < n_prec; ++i)
    {
      const double mz = fromFixed(p.mz[i]);
      bool placed = false;
      for (std::size_t w = 0; w < windows.size(); ++w)
      {
        if (windows[w].contains(mz))
        {
          by_window[w].push_back(static_cast<std::uint32_t>(i));
          placed = true;
        }
      }
      if (!placed) { ++out.precursors_without_window; }
    }

    const std::size_t n_trans =
      n_prec == 0 ? 0 : p.transition_begin[n_prec - 1] + p.transition_count[n_prec - 1];
    out.begin.assign(n_trans, 0);
    out.count.assign(n_trans, 0);

    // Points per transition: one per spectrum of every window carrying its
    // precursor. Counted first so the flat arrays are allocated once rather
    // than grown, which is what keeps this to two allocations at proteome
    // scale instead of millions.
    std::vector<std::uint32_t> per_transition(n_trans, 0);
    std::vector<std::size_t> window_spectra(windows.size(), 0);
    for (const auto& s : info)
    {
      if (!in_range(s)) { continue; }
      for (std::size_t w = 0; w < windows.size(); ++w)
      {
        if (sameWindow(s.window, windows[w])) { ++window_spectra[w]; break; }
      }
    }
    for (std::size_t w = 0; w < windows.size(); ++w)
    {
      for (const auto i : by_window[w])
      {
        for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
        {
          per_transition[p.transition_begin[i] + k] += window_spectra[w];
        }
      }
    }
    std::uint64_t total = 0;
    for (std::size_t j = 0; j < n_trans; ++j)
    {
      out.begin[j] = static_cast<std::uint32_t>(total);
      total += per_transition[j];
      if (total > std::numeric_limits<std::uint32_t>::max())
      {
        throw std::runtime_error("more than 2^32 chromatogram points; restrict the "
                                 "retention-time range or the precursor count");
      }
    }
    out.retention_time.assign(total, 0.0f);
    out.intensity.assign(total, 0.0f);
    std::vector<std::uint32_t> cursor = out.begin;

    // One forward pass. Every transition whose precursor falls in this
    // spectrum's window is probed while the spectrum is in hand, so the decode
    // is paid once however large the library is.
    // Resolve the retention-time range to an index range BEFORE decoding
    // anything. Filtering after the decode reads the whole run and throws most
    // of it away, which at 294 ms/spectrum is the difference between seconds
    // and an hour. The spectra are sorted by retention time, so this is a
    // binary search.
    std::size_t first = 0, last = info.size();
    if (options.rt_high > options.rt_low)
    {
      first = static_cast<std::size_t>(
        std::lower_bound(info.begin(), info.end(), options.rt_low,
                         [](const SpectrumInfo& s, double v) {
                           return s.retention_time < v;
                         }) - info.begin());
      last = static_cast<std::size_t>(
        std::upper_bound(info.begin(), info.end(), options.rt_high,
                         [](double v, const SpectrumInfo& s) {
                           return v < s.retention_time;
                         }) - info.begin());
    }

    std::vector<SpectrumPeaks> block;
    constexpr std::size_t BLOCK = 64;
    for (std::size_t begin = first; begin < last; begin += BLOCK)
    {
      const std::size_t end = std::min(begin + BLOCK, last);
      const auto t_decode = std::chrono::steady_clock::now();
      source.peaks(begin, end, block);
      st.decode_seconds += std::chrono::duration<double>(
                             std::chrono::steady_clock::now() - t_decode).count();

      const auto t_match = std::chrono::steady_clock::now();
      for (std::size_t si = begin; si < end; ++si)
      {
        const auto& s = info[si];
        if (!in_range(s)) { continue; }
        const auto& peaks = block[si - begin];
        ++st.spectra_read;

        for (std::size_t w = 0; w < windows.size(); ++w)
        {
          if (!sameWindow(s.window, windows[w])) { continue; }
          for (const auto i : by_window[w])
          {
            for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
            {
              const std::size_t j = p.transition_begin[i] + k;
              const double mz = fromFixed(t.product_mz[j]);
              const double tol = mz * options.fragment_ppm * 1e-6;
              const float value = probe(peaks, mz, tol, s.window.im_low, s.window.im_high,
                                        options.use_ion_mobility);
              const std::size_t at = cursor[j]++;
              out.retention_time[at] = static_cast<float>(s.retention_time);
              out.intensity[at] = value;
              if (value > 0.0f) { ++st.nonzero_points; }
            }
          }
          break;
        }
      }
      st.match_seconds += std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - t_match).count();

      if (options.progress_every && (begin / BLOCK) % 8 == 0)
      {
        std::cerr << "\r  " << st.spectra_read << " / " << (last - first) << " spectra"
                  << std::flush;
      }
    }
    if (options.progress_every) { std::cerr << "\r" << std::string(48, ' ') << "\r"; }

    for (std::size_t j = 0; j < n_trans; ++j) { out.count[j] = cursor[j] - out.begin[j]; }
    st.precursors = n_prec;
    st.transitions = n_trans;
    st.points = out.points();
    return out;
  }

} // namespace ODIA
