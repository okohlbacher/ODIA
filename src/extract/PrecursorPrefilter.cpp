// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/PrecursorPrefilter.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <sstream>

namespace ODIA
{
  namespace
  {
    struct Frag { double mz; float intensity; };

    /// One searchable fragment. 16 bytes, and deliberately WITHOUT the
    /// retention-time bounds: those are per precursor, not per fragment, and
    /// storing them here would add 8 bytes to every one of ~30M entries to
    /// avoid a lookup into a 40 MB array.
    struct Target
    {
      double mz;
      std::uint32_t slot;      ///< FULL-library precursor index, never renumbered
      float im;
    };
  }

  std::vector<PrecursorPrefilter::Evidence>
  PrecursorPrefilter::measure(const Library& library, SpectrumSource& source,
                              const Options& options, Stats& stats)
  {
    const auto t0 = std::chrono::steady_clock::now();
    const auto& p = library.precursors();
    const auto& t = library.transitions();
    const std::size_t n_prec = library.precursorCount();

    std::vector<Evidence> ev(n_prec);

    std::vector<std::vector<Frag>> want(n_prec);
    for (std::size_t i = 0; i < n_prec; ++i)
    {
      const std::uint32_t b = p.transition_begin[i], n = p.transition_count[i];
      std::vector<Frag> f;
      f.reserve(n);
      for (std::uint32_t k = 0; k < n; ++k)
      {
        const auto mz = t.product_mz[b + k];
        if (mz == MZ_INVALID) { continue; }
        f.push_back({fromFixed(mz), t.library_intensity[b + k]});
      }
      std::sort(f.begin(), f.end(),
                [](const Frag& a, const Frag& b2) { return a.intensity > b2.intensity; });
      if (f.size() > options.top_n) { f.resize(options.top_n); }
      std::sort(f.begin(), f.end(),
                [](const Frag& a, const Frag& b2) { return a.mz < b2.mz; });
      want[i] = std::move(f);
    }

    // The retention-time neighbourhood, per precursor. This is what running
    // BETWEEN the passes buys: pass 1 has already fitted the map.
    const bool rt_gated = options.rt_half_window > 0.0 && options.irt_slope != 0.0;
    const float ninf = -std::numeric_limits<float>::infinity();
    const float pinf =  std::numeric_limits<float>::infinity();
    std::vector<float> rt_lo(n_prec, ninf), rt_hi(n_prec, pinf);
    if (rt_gated)
    {
      for (std::size_t i = 0; i < n_prec; ++i)
      {
        const double c = options.irt_slope * double(p.irt[i]) + options.irt_intercept;
        rt_lo[i] = static_cast<float>(c - options.rt_half_window);
        rt_hi[i] = static_cast<float>(c + options.rt_half_window);
      }
    }

    const auto& info = source.spectra();
    const auto& windows = source.windows();

    std::vector<std::vector<Target>> window_index(windows.size());
    std::vector<char> window_built(windows.size(), 0);

    // Which window each spectrum belongs to, resolved once. Doing it per
    // spectrum inside the sweep is a linear scan over the window list for every
    // one of ~32,000 spectra.
    std::vector<std::uint32_t> window_of(info.size(),
                                         std::numeric_limits<std::uint32_t>::max());
    for (std::size_t si = 0; si < info.size(); ++si)
    {
      if (info[si].ms_level != 2) { continue; }
      for (std::size_t k = 0; k < windows.size(); ++k)
      {
        if (windows[k].mz_low == info[si].window.mz_low &&
            windows[k].mz_high == info[si].window.mz_high)
        { window_of[si] = static_cast<std::uint32_t>(k); break; }
      }
    }

    // `hits` is allocated ONCE and cleared SPARSELY -- only the slots this
    // spectrum touched. Re-assigning it per spectrum is O(precursors) per
    // spectrum, which at 4,986,319 precursors and ~32,000 spectra is 1.6e11
    // writes and would dominate everything else here by orders of magnitude.
    std::vector<std::uint8_t> hits(n_prec, 0);
    std::vector<std::uint32_t> slots_touched;
    slots_touched.reserve(1u << 16);

    // Contiguity state. A run is consecutive in CYCLES OF ITS OWN WINDOW, not
    // in spectrum index: a run interleaves its isolation windows, so successive
    // spectra of one precursor's window are tens of indices apart.
    const std::uint32_t NO_CYCLE = std::numeric_limits<std::uint32_t>::max();
    const std::size_t qual_depth = options.contiguity_depth > 0
                                     ? options.contiguity_depth
                                     : (options.top_n > 1 ? options.top_n - 1 : 1);
    std::vector<std::uint32_t> last_cycle(n_prec, NO_CYCLE);
    std::vector<std::uint16_t> cur_run(n_prec, 0);
    std::vector<float> run_start_rt(n_prec, -1.0f);
    std::vector<std::uint32_t> window_cycle(windows.size(), 0);

    constexpr std::size_t BLOCK = 256;
    std::vector<SpectrumPeaks> block;

    for (std::size_t begin = 0; begin < info.size(); begin += BLOCK)
    {
      const std::size_t end = std::min(begin + BLOCK, info.size());
      source.peaks(begin, end, block);

      for (std::size_t si = begin; si < end; ++si)
      {
        const std::uint32_t w = window_of[si];
        if (w == std::numeric_limits<std::uint32_t>::max()) { continue; }
        const SpectrumPeaks& sp = block[si - begin];
        if (sp.mz.empty()) { continue; }

        if (!window_built[w])
        {
          auto& idx = window_index[w];
          for (std::size_t i = 0; i < n_prec; ++i)
          {
            if (!info[si].window.contains(fromFixed(p.mz[i]))) { continue; }
            for (const auto& fr : want[i])
            { idx.push_back({fr.mz, static_cast<std::uint32_t>(i), p.im[i]}); }
          }
          std::sort(idx.begin(), idx.end(),
                    [](const Target& a, const Target& b) { return a.mz < b.mz; });
          window_built[w] = 1;
        }
        const auto& idx = window_index[w];
        if (idx.empty()) { continue; }

        const bool im_gated = options.im_window > 0.0 &&
                              sp.ion_mobility.size() == sp.mz.size();
        const float rt = static_cast<float>(info[si].retention_time);
        const std::uint32_t cycle = window_cycle[w]++;

        for (std::size_t pk = 0; pk < sp.mz.size(); ++pk)
        {
          // Peaks are NOT ascending in m/z -- a merged mobility frame restarts
          // m/z at every TIMS scan -- so the LIBRARY side is the sorted side
          // and the peaks are iterated. Binary-searching the peak array does
          // not fail loudly: it returns near-zero matches, indistinguishable
          // from an ion that is simply not there.
          //
          // The calibration is applied to the OBSERVED m/z rather than to each
          // target: one multiply here against six per precursor there, and the
          // two are equivalent to first order.
          const double m = sp.mz[pk] * (1.0 - options.ppm_centre * 1e-6);
          const double tol = m * options.ppm * 1e-6;
          auto it = std::lower_bound(idx.begin(), idx.end(), m - tol,
                                     [](const Target& a, double v) { return a.mz < v; });
          for (; it != idx.end() && it->mz <= m + tol; ++it)
          {
            if (std::abs(it->mz - m) > it->mz * options.ppm * 1e-6) { continue; }
            // The retention-time neighbourhood -- the reason this runs between
            // the passes. It cuts the number of DRAWS the depth maximum is
            // taken over, which doc/08 named as why the statistic saturated.
            if (rt < rt_lo[it->slot] || rt > rt_hi[it->slot]) { continue; }
            // A precursor with NO library 1/K0 is UNGATED, matching what the
            // extractor does and for the same reason. it->im is NaN whenever
            // the library carries no mobility -- every comparison against NaN
            // is false, so the negation below rejects every peak and the
            // precursor scores zero on a run it may well be present in.
            //
            // Measured: on S08 the seed sweep returned 0 target and 0 decoy
            // anchors at every contiguity threshold, against 13,360 / 8,264 on
            // Astral, purely because our generated library predicts CCS and
            // leaves 1/K0 unset. Astral has no ion mobility so the gate never
            // ran there, which is why the failure looked instrument-specific
            // rather than like a missing NaN case.
            if (im_gated && !std::isnan(it->im))
            {
              const float pim = sp.ion_mobility[pk];
              if (!(std::abs(double(pim) - double(it->im)) <= options.im_window))
              { continue; }
            }
            if (hits[it->slot] == 0) { slots_touched.push_back(it->slot); }
            // Depth counts DISTINCT fragments. One fragment matches many peaks
            // in a mobility-merged frame, so a fragment already counted must
            // not be able to raise the depth again.
            if (hits[it->slot] < options.top_n)
            { hits[it->slot] = static_cast<std::uint8_t>(hits[it->slot] + 1); }
          }
        }

        for (const std::uint32_t slot : slots_touched)
        {
          const std::uint8_t h = std::min<std::uint8_t>(
            hits[slot], static_cast<std::uint8_t>(want[slot].size()));
          if (h > 0)
          {
            ev[slot].total_matches += h;
            if (h > ev[slot].depth) { ev[slot].depth = h; ev[slot].best_rt = rt; }

            // A chance coincidence wins ONE cycle. It does not easily win the
            // next one too, because the interfering ions that produced it are
            // not eluting on this precursor's peak.
            if (h >= qual_depth)
            {
              if (last_cycle[slot] != NO_CYCLE && last_cycle[slot] + 1 == cycle)
              { ++cur_run[slot]; }
              else
              { cur_run[slot] = 1; run_start_rt[slot] = rt; }
              last_cycle[slot] = cycle;
              if (cur_run[slot] > ev[slot].contiguity)
              {
                ev[slot].contiguity = cur_run[slot];
                // The midpoint of the elution peak, not the single spectrum
                // that happened to match best.
                ev[slot].contiguous_rt = 0.5f * (run_start_rt[slot] + rt);
              }
            }
          }
          hits[slot] = 0;                      // the sparse clear
        }
        slots_touched.clear();
        ++stats.spectra_swept;
      }
    }

    stats.seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - t0).count();
    return ev;
  }

  std::vector<char>
  PrecursorPrefilter::select(const Library& library,
                             const std::vector<Evidence>& evidence,
                             const Options& options, Stats& stats)
  {
    const auto& p = library.precursors();
    const std::size_t n_prec = library.precursorCount();
    std::vector<char> keep(n_prec, 1);

    std::vector<std::uint32_t> tgt, dec;
    for (std::size_t i = 0; i < n_prec; ++i)
    { (p.decoy[i] ? dec : tgt).push_back(static_cast<std::uint32_t>(i)); }

    stats.targets_in = tgt.size();
    stats.decoys_in  = dec.size();

    const std::size_t max_depth = options.top_n + 1;
    stats.depth_hist_target.assign(max_depth, 0);
    stats.depth_hist_decoy.assign(max_depth, 0);
    for (std::size_t i = 0; i < n_prec; ++i)
    {
      const std::size_t d = std::min<std::size_t>(evidence[i].depth, max_depth - 1);
      ++(p.decoy[i] ? stats.depth_hist_decoy : stats.depth_hist_target)[d];
    }

    constexpr std::size_t CONTIG_BINS = 11;   // 0..9 and "10 or more"
    stats.contig_hist_target.assign(CONTIG_BINS, 0);
    stats.contig_hist_decoy.assign(CONTIG_BINS, 0);
    for (std::size_t i = 0; i < n_prec; ++i)
    {
      const std::size_t c = std::min<std::size_t>(evidence[i].contiguity, CONTIG_BINS - 1);
      ++(p.decoy[i] ? stats.contig_hist_decoy : stats.contig_hist_target)[c];
    }

    if (options.keep_fraction >= 1.0)
    {
      stats.targets_kept = tgt.size();
      stats.decoys_kept  = dec.size();
      stats.note = "measurement only -- nothing discarded";
      return keep;
    }

    // Rank within each label class and cut each at the same COUNT -- not at the
    // same depth. Targets clear an evidence bar more often, so a shared
    // threshold would retain a biased, weaker decoy sample, and the decoys
    // would then score below a fair null and make every downstream q-value
    // optimistic. Equal counts retain the BEST decoys, biasing the FDR
    // conservative, which is the safe direction. See the header.
    const std::size_t per_class = std::min(
      std::min(tgt.size(), dec.size()),
      std::max(options.min_keep_per_class,
               static_cast<std::size_t>(options.keep_fraction *
                                        double(std::min(tgt.size(), dec.size())))));

    auto better = [&](std::uint32_t a, std::uint32_t b)
    {
      if (evidence[a].depth != evidence[b].depth)
      { return evidence[a].depth > evidence[b].depth; }
      return evidence[a].total_matches > evidence[b].total_matches;
    };

    auto cut = [&](std::vector<std::uint32_t>& v) -> std::uint8_t
    {
      if (v.size() <= per_class) { return 0; }
      std::nth_element(v.begin(), v.begin() + std::ptrdiff_t(per_class), v.end(), better);
      for (std::size_t k = per_class; k < v.size(); ++k) { keep[v[k]] = 0; }
      // The retained set is v[0, per_class) but is only PARTITIONED, not
      // sorted, so v[per_class - 1] is an arbitrary retained element and not
      // the weakest one. The cut depth is the minimum over what was retained.
      std::uint8_t lowest = 255;
      for (std::size_t k = 0; k < per_class; ++k)
      { lowest = std::min(lowest, evidence[v[k]].depth); }
      v.resize(per_class);
      return lowest;
    };
    const std::uint8_t cut_t = cut(tgt);
    cut(dec);

    stats.targets_kept = tgt.size();
    stats.decoys_kept  = dec.size();
    stats.depth_threshold = cut_t;
    return keep;
  }

} // namespace ODIA
