// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/Ms1Traces.h>
#include <limits>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstring>
#include <exception>
#include <iterator>
#include <sstream>
#include <thread>

namespace ODIA
{
  namespace
  {
    /// FNV-1a over 32-bit words of the values' BIT PATTERNS, for the identity
    /// gates: a quarter of the multiplies of a byte-wise hash over 13.3G cells,
    /// and it is only ever compared against itself. Bits, not values, so -0.0
    /// against 0.0 or two NaN payloads count as the difference they are.
    struct Fnv
    {
      std::uint64_t h = 14695981039346656037ull;
      void mix(std::uint32_t w) { h ^= w; h *= 1099511628211ull; }
      void add(float x) { std::uint32_t w; std::memcpy(&w, &x, sizeof w); mix(w); }
      void add(double x)
      {
        std::uint64_t w;
        std::memcpy(&w, &x, sizeof w);
        mix(static_cast<std::uint32_t>(w));
        mix(static_cast<std::uint32_t>(w >> 32));
      }
    };
  }

  std::size_t Ms1Traces::binFor(double rt) const
  {
    if (times_.empty()) { return 0; }
    const auto it = std::lower_bound(times_.begin(), times_.end(), static_cast<float>(rt));
    std::size_t b = static_cast<std::size_t>(it - times_.begin());
    if (b >= times_.size()) { return times_.size() - 1; }
    if (b > 0 && std::abs(times_[b - 1] - rt) < std::abs(times_[b] - rt)) { --b; }
    return b;
  }

  std::string Ms1Traces::describe() const
  {
    std::ostringstream o;
    if (empty()) { return "MS1 traces: none (the run carries no MS1)"; }
    std::size_t live = 0;
    for (std::size_t i = 0; i < precursors(); ++i)
    {
      for (std::size_t b = 0; b < bins_; ++b)
      { if (values_[i * bins_ + b] > 0.0f) { ++live; break; } }
    }
    o.precision(1);
    o << std::fixed << "MS1 traces: " << bins_ << " bins over " << precursors()
      << " precursors, " << live << " with signal ("
      << (precursors() ? 100.0 * live / precursors() : 0.0) << "%), "
      << footprintBytes() / 1048576.0 << " MiB";
    return o.str();
  }

  std::string Ms1Traces::describeBuild() const
  {
    std::ostringstream o;
    o.setf(std::ios::fixed);
    o.precision(1);
    const double mib = 1048576.0;
    o << "MS1 build: " << stats_.frames << " frames in " << stats_.blocks
      << " decodes, " << stats_.match_blocks << " match blocks (largest "
      << stats_.max_block_bytes / mib << " MiB, at most " << stats_.max_held_bytes / mib
      << " MiB decoded at once, residual buffers <= " << stats_.max_resid_bytes / mib
      << " MiB; cap " << stats_.block_cap_bytes / mib << " MiB beyond the first decode, "
      << stats_.carried << " decode(s) carried), "
      << stats_.peaks << " peaks in all; decode " << stats_.decode_s
      << " s on the driver, match " << stats_.match_s << " s over " << stats_.units
      << " frame runs on " << stats_.threads << " thread(s)";
    return o.str();
  }

  std::uint64_t Ms1Traces::checksum() const
  {
    // Shape first, so a transposed or truncated matrix cannot collide with the
    // real one.
    Fnv h;
    h.mix(static_cast<std::uint32_t>(bins_));
    h.mix(static_cast<std::uint32_t>(precursors()));
    for (const float x : times_) { h.add(x); }
    for (const float x : values_) { h.add(x); }
    return h.h;
  }

  Ms1Traces Ms1Traces::build(const Library& library, SpectrumSource& source,
                             double fragment_ppm, double im_window,
                             double ppm_offset, double* observed_ppm_median,
                             double isotope_offset_da,
                             const std::vector<std::uint8_t>* keep,
                             std::vector<std::uint32_t>* kept_indices,
                             unsigned threads, std::size_t block_budget_bytes)
  {
    Ms1Traces out;
    const auto& ms1 = source.ms1Spectra();
    if (ms1.empty()) { return out; }

    const auto& p = library.precursors();
    const std::size_t np = library.precursorCount();
    out.bins_ = ms1.size();
    out.times_.reserve(ms1.size());
    for (const auto& s : ms1) { out.times_.push_back(static_cast<float>(s.retention_time)); }

    // Row assignment. Without a mask each library precursor owns row i (the
    // scorer's contract). With a mask, rows are assigned in ascending library
    // index over the kept precursors, and the caller gets that order back via
    // kept_indices -- the writer reconstructs ids from it, so the mapping is
    // never inferred twice.
    std::vector<std::uint32_t> row;
    std::size_t rows = np;
    if (keep != nullptr)
    {
      row.assign(np, UINT32_MAX);
      std::uint32_t r = 0;
      for (std::size_t i = 0; i < np; ++i)
      {
        if (i < keep->size() && (*keep)[i]) { row[i] = r++; }
      }
      rows = r;
      if (kept_indices != nullptr)
      {
        kept_indices->clear();
        kept_indices->reserve(rows);
        for (std::size_t i = 0; i < np; ++i)
        { if (row[i] != UINT32_MAX) { kept_indices->push_back(static_cast<std::uint32_t>(i)); } }
      }
    }
    out.values_.assign(rows * out.bins_, 0.0f);

    // Search the sorted LIBRARY side and iterate the peaks: SpectrumSource
    // documents that a peak array is not ascending in m/z (a mobility frame
    // concatenates its TIMS scans), and a binary search over it "does not fail
    // loudly -- it returns near-zero matches", which is indistinguishable from
    // an ion that is not there.
    struct Target { double mz; std::uint32_t slot; float im; };
    std::vector<Target> idx;
    idx.reserve(np);
    for (std::size_t i = 0; i < np; ++i)
    {
      if (keep != nullptr && row[i] == UINT32_MAX) { continue; }
      // CALIBRATED, like the fragment axis. This matched on the library's
      // THEORETICAL m/z with a symmetric window and no offset, while the
      // fragment extractor was centred on the fitted deviation -- on IH1 that
      // is -10.0108 ppm against a +/-10 ppm half-width, so a precursor whose
      // MS1 error resembles its MS2 error sat at the window EDGE and a weak one
      // fell out entirely. That matters because ms1_coelution is the main
      // evidence for calling a precursor ABSENT (median -0.093 for the 10,736
      // DIA-NN precursors we reject, against -0.124 for the 801,458 bulk
      // non-identifications and 0.473 for accepted ones), and absence cannot be
      // concluded from an uncalibrated measurement.
      //
      // The isotope offset is applied BEFORE the calibration scaling, and per
      // this precursor's own charge: the M+k target is `mz + k*dm/z`, and the
      // instrument's relative (ppm) error then applies to that target as it
      // does to any mass.
      const int z = p.charge[i] > 0 ? static_cast<int>(p.charge[i]) : 1;
      const double mz = (fromFixed(p.mz[i]) + isotope_offset_da / z) *
                        (1.0 + ppm_offset * 1e-6);
      const std::uint32_t slot = keep != nullptr ? row[i] : static_cast<std::uint32_t>(i);
      if (mz > 0.0) { idx.push_back({mz, slot, p.im[i]}); }
    }
    std::sort(idx.begin(), idx.end(),
              [](const Target& a, const Target& b) { return a.mz < b.mz; });

    // CAPPED. The first version pushed one double per (peak, target) match over
    // the whole run and reached 591 GB RSS against v3's 116 GB peak, blowing
    // through -live_memory_gb on a shared node. The median of a bounded prefix
    // is the same number to far more precision than it is worth: this is a
    // diagnostic, not a fit.
    static constexpr std::size_t RESID_CAP = 1u << 21;   // 2M samples, 16 MB
    std::vector<double> resid;
    if (observed_ppm_median != nullptr) { resid.reserve(RESID_CAP); }

    // Match decoded frames [lo, hi) of @p block, whose first frame is MS1
    // spectrum @p b. Residuals go to @p res while it holds fewer than @p cap.
    //
    // Safe to run on disjoint frame ranges at once, which is the whole of the
    // threading argument (doc/15:187-192, doc/83 F05): a write lands at
    // `slot * bins_ + (b + s)`, so its COLUMN is the frame and two ranges never
    // share a cell; and within a cell the updates arrive in the same peak and
    // target order as the serial loop, so even the Max is not being asked to
    // be order-free. `idx` and the decoded block are read-only here. The one
    // ordered thing is the residual prefix, which is why it is a parameter.
    auto match = [&](const std::vector<SpectrumPeaks>& block, std::size_t b,
                     std::size_t lo, std::size_t hi,
                     std::vector<double>* res, std::size_t cap)
    {
      for (std::size_t s = lo; s < hi; ++s)
      {
        const auto& sp = block[s];
        const bool gated = im_window > 0.0 && sp.ion_mobility.size() == sp.mz.size();
        for (std::size_t k = 0; k < sp.mz.size(); ++k)
        {
          const double m = sp.mz[k], tol = m * fragment_ppm * 1e-6;
          auto it = std::lower_bound(idx.begin(), idx.end(), m - tol,
                                     [](const Target& a, double v) { return a.mz < v; });
          for (; it != idx.end() && it->mz <= m + tol; ++it)
          {
            if (std::abs(it->mz - m) > it->mz * fragment_ppm * 1e-6) { continue; }
            // Same rule as the MS2 match loop, deliberately: skip the gate when
            // EITHER side is unknown, because absent information is not evidence
            // of mismatch. Testing `abs(NaN - x) <= w` is false, so a precursor
            // with no library 1/K0 had every MS1 peak rejected and came out with
            // a NaN MS1_COELUTION -- while its MS2 side was extracted ungated.
            if (gated && !std::isnan(static_cast<double>(it->im)))
            {
              const double d = std::abs(static_cast<double>(sp.ion_mobility[k]) -
                                        static_cast<double>(it->im));
              if (!(d <= im_window)) { continue; }
            }
            if (res != nullptr && res->size() < cap)
            {
              // Residual against the CALIBRATED target, so a correct offset
              // centres this on 0 and a wrong one does not.
              res->push_back((m - it->mz) / it->mz * 1e6);
            }
            float& c = out.values_[it->slot * out.bins_ + (b + s)];
            // Max, not sum: a mobility-merged frame holds the same ion in
            // several scans, and summing would make the trace a function of how
            // many scans it spans rather than of how much ion is present.
            c = std::max(c, sp.intensity[k]);
          }
        }
      }
    };

    // FRAME-PARALLEL over contiguous runs of frames, decode on this thread.
    //
    // This loop was serial and on the reference arm it is 5,589.6 s, 16.3% of
    // the wall, with 63 of 64 cores idle (an IH1 full run; doc/83 F05).
    // The per-frame cost is 3.17-4.16 s at full scale against 0.56 s on the
    // 6x60-s fixture -- a memory-size effect: the lookup over 9.9M targets and
    // the scattered writes into a 50 GiB matrix. The decode share is INFERRED
    // at 150-300 s of the 5,590 (no log split it; `BuildStats` now does), so
    // the matching is what is worth parallelising, and the fixture's speedup
    // will not transport to the full run.
    //
    // DECODE STAYS HERE. The mzPeak reader holds one Index and one Spectra
    // handle (MzPeakSource.cpp:454-455) and is not thread-safe; only the match
    // over a decoded block runs on several threads.
    //
    // THE UNIT IS A CONTIGUOUS RUN OF >= MIN_RUN FRAMES, not one frame. A
    // precursor's row is 1,343 floats = 5,372 B on IH1, not a multiple of a
    // cache line, and a co-eluting ion lights ~5 ADJACENT bins: with one frame
    // per thread, neighbouring threads write the same line of the same row
    // and the line ping-pongs between cores (doc/83 F05, from codex). Sixteen
    // frames are 64 B, so a unit shares at most its two edge lines.
    //
    // ONE WRITER PER CELL, checked rather than assumed (codex, doc/83 Q6.4: if
    // a frame spread over neighbouring bins, adjacent runs would read-max-write
    // the same cell). It does not: `bins_` is the MS1 frame count and the
    // column is `b + s`, the frame's own index, with no smoothing or binning
    // in between -- so a cell (row, frame) is written only by the unit that
    // owns the frame. Adjacent units may share a cache LINE, which costs time,
    // never a race: distinct floats are distinct memory locations.
    //
    // THE MATCH BLOCK GROWS WITH THE THREADS, MIN_RUN per thread, UP TO A
    // BYTE CAP; THE DECODE CALLS DO NOT GROW. At the serial 64 frames a block
    // holds only four units, which caps the speedup at 4x -- exactly F05's own
    // FAIL bar. So a match block is several decode blocks: `ms1Peaks` is still
    // asked for [64k, 64k + 64), the very ranges the serial loop asks for, in
    // the same order, and the frames are moved into one vector. A decoded
    // frame depends only on its own spectrum (MzPeakSource::ms1Peaks copies
    // batch[k] into out[k]) but this way that is not even needed.
    //
    // THE PRICE IS DECODED PEAKS HELD AT ONCE, AND IT IS CAPPED AND CHARGED
    // (review round 2, verified on the PXD fixture: an uncapped block was
    // MIN_RUN * threads frames, 735 MiB serial, 1,195 MiB at -threads 8,
    // 2,188 MiB = the whole 215-frame fixture at -threads 48 -- and on a full
    // run the whole decoded run, 23-40 GiB, once -threads reaches ~frames/16,
    // none of it checked against -live_memory_gb). A further decode joins the
    // block only while what the block holds BEYOND its first decode -- the
    // extra decodes plus the worst case of the per-unit residual buffers --
    // stays within `cap` = min(block_budget_bytes, MAX_BLOCK_BYTES). The
    // caller passes what the memory budget leaves after the matrix, so the
    // parallel path never holds more than the serial one plus that remainder.
    // A decode's size is known only after it is decoded, so the next one is
    // made only if a decode as large as the largest so far would fit; one
    // that is decoded and then does not fit (it outgrew every earlier decode)
    // is CARRIED to the next block rather than re-decoded, and only then are
    // a block plus one decode alive at once (`max_held_bytes`, `carried`).
    // The first decode of a block is always taken (the serial floor, which no
    // budget has ever governed), so cap 0 means one 64-frame decode per block
    // and nothing else held -- the serial footprint -- still matched as four
    // units. Block
    // boundaries never reach the output -- a cell's column is its frame, and
    // the residual prefix is rebuilt in frame order whatever the blocks -- so
    // the cap changes memory and speed, never a byte of the result.
    //
    // Threads are spawned per match block, not pooled as in the extractor: a
    // block is >= 64 frames against 0.5-4 s of matching per frame.
    //
    // THE RESIDUAL PREFIX IS REBUILT EXACTLY. Serial semantics: `resid` holds
    // the first RESID_CAP residuals in (frame, peak, target) order. Each unit
    // collects its own first `budget` residuals, `budget` being what the
    // prefix still lacks when the block starts; the units are then appended in
    // frame order and the total cut at RESID_CAP. No unit can need more than
    // `budget` of its own, so the cut prefix is the serial one sample for
    // sample -- the retained LIST, not merely its median (codex Q6.4(b)) --
    // and `resid_hash` prints it so a gate can compare it. Each buffer is
    // capped as the serial one is: the uncapped original reached 591 GB.
    // The price is transient -- while the prefix is still short, every unit of
    // a block may fill up to `budget` (16 MB) -- and it is charged to the
    // block cap above (one unit per started MIN_RUN frames of each decode),
    // so it shrinks the first blocks rather than adding to them; nothing
    // afterwards (budget 0 allocates no buffers).
    //
    // threads 0 and 1 are both serial, and the serial path is the old loop
    // with the old 64-frame block, so default and -threads 1 output is the old
    // output by construction.
    static constexpr std::size_t DECODE = 64;
    static constexpr std::size_t MIN_RUN = 16;
    const unsigned T = threads > 1 ? threads : 1u;
    // A multiple of DECODE, so the decode ranges below are the serial ones.
    // The frame target: no more frames than T threads have units for.
    const std::size_t STEP = T > 1 ? (MIN_RUN * T + DECODE - 1) / DECODE * DECODE : DECODE;
    const std::size_t cap = std::min(block_budget_bytes, MAX_BLOCK_BYTES);
    out.stats_ = BuildStats{};
    out.stats_.frames = ms1.size();
    out.stats_.block_cap_bytes = cap;
    auto bytesOf = [](const std::vector<SpectrumPeaks>& v)
    {
      std::size_t bytes = 0;
      for (const auto& sp : v)
      {
        bytes += sp.mz.capacity() * sizeof(double) +
                 sp.intensity.capacity() * sizeof(float) +
                 sp.ion_mobility.capacity() * sizeof(float);
      }
      return bytes;
    };
    std::vector<SpectrumPeaks> block, part;
    std::size_t next_frame = 0;      // first frame not yet decoded
    std::size_t part_lo = 0, part_bytes = 0, max_decode_bytes = 0;
    bool carried = false;            // `part` holds a decode the last block refused
    while (carried || next_frame < ms1.size())
    {
      const std::size_t b = carried ? part_lo : next_frame;
      // What the prefix still lacks: the per-unit residual buffers' bound.
      const std::size_t budget = observed_ppm_median != nullptr && resid.size() < RESID_CAP
                                   ? RESID_CAP - resid.size() : 0;
      const std::size_t resid_unit_bytes = T > 1 ? budget * sizeof(double) : 0;
      block.clear();
      std::size_t block_bytes = 0, charged = 0;
      double dec_s = 0.0;
      for (;;)
      {
        if (!carried)
        {
          if (next_frame >= ms1.size() || block.size() >= STEP) { break; }
          const std::size_t e = std::min(next_frame + DECODE, ms1.size());
          // Predicted from the largest decode so far, so a decode that cannot
          // fit is not made at all -- at cap 0 never -- and a carry happens
          // only when a decode outgrows every earlier one.
          if (!block.empty() &&
              charged + max_decode_bytes + (e - next_frame + MIN_RUN - 1) / MIN_RUN * resid_unit_bytes > cap)
          { break; }
          const auto t_dec = std::chrono::steady_clock::now();
          source.ms1Peaks(next_frame, e, part);
          dec_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_dec).count();
          ++out.stats_.blocks;
          part_lo = next_frame;
          part_bytes = bytesOf(part);
          max_decode_bytes = std::max(max_decode_bytes, part_bytes);
          next_frame = e;
        }
        carried = false;
        const std::size_t resid_cost = (part.size() + MIN_RUN - 1) / MIN_RUN * resid_unit_bytes;
        const std::size_t cost = (block.empty() ? 0 : part_bytes) + resid_cost;
        if (!block.empty() && charged + cost > cap)
        {
          // Did not fit: keep it decoded for the next block. Both are alive
          // now, which is the most this loop ever holds.
          out.stats_.max_held_bytes = std::max(out.stats_.max_held_bytes, block_bytes + part_bytes);
          ++out.stats_.carried;
          carried = true;
          break;
        }
        charged += cost;
        block_bytes += part_bytes;
        if (block.empty()) { block.swap(part); }
        else
        {
          block.insert(block.end(), std::make_move_iterator(part.begin()),
                       std::make_move_iterator(part.end()));
        }
      }
      const auto t_match = std::chrono::steady_clock::now();
      out.stats_.decode_s += dec_s;
      ++out.stats_.match_blocks;
      for (const auto& sp : block) { out.stats_.peaks += sp.mz.size(); }
      out.stats_.max_block_bytes = std::max(out.stats_.max_block_bytes, block_bytes);
      out.stats_.max_held_bytes = std::max(out.stats_.max_held_bytes, block_bytes);

      const std::size_t n = block.size();
      const std::size_t runs = T > 1 ? std::max<std::size_t>(1, n / MIN_RUN) : 1;
      std::vector<double>* const serial_res = observed_ppm_median != nullptr ? &resid : nullptr;
      if (runs == 1)
      {
        match(block, b, 0, n, serial_res, RESID_CAP);
        ++out.stats_.units;
      }
      else
      {
        out.stats_.max_resid_bytes = std::max(out.stats_.max_resid_bytes, runs * resid_unit_bytes);
        std::vector<std::vector<double>> unit_res(budget > 0 ? runs : 0);
        // Pulled, not dealt: frames in the middle of the gradient cost more
        // than those at its ends, and an atomic counter lets a thread that
        // drew a cheap run take another. Which thread does a run never
        // reaches the output -- only the run boundaries do, and those are
        // fixed by `n` and `runs`.
        std::atomic<std::size_t> next{0};
        const unsigned workers = static_cast<unsigned>(std::min<std::size_t>(T, runs));
        std::vector<std::exception_ptr> failed(workers);
        auto worker = [&](unsigned w)
        {
          try
          {
            for (std::size_t u = next.fetch_add(1); u < runs; u = next.fetch_add(1))
            {
              // Balanced split: every run has floor or ceil of n/runs frames,
              // and n/runs >= MIN_RUN because runs = floor(n / MIN_RUN).
              match(block, b, u * n / runs, (u + 1) * n / runs,
                    budget > 0 ? &unit_res[u] : nullptr, budget);
            }
          }
          catch (...) { failed[w] = std::current_exception(); }
        };
        {
          std::vector<std::thread> pool;
          // Joins on every exit, including a thread that fails to start: the
          // ones already running drain the counter and are waited for, rather
          // than a joinable std::thread reaching its destructor and aborting.
          struct Join
          {
            std::vector<std::thread>& t;
            ~Join() { for (auto& x : t) { if (x.joinable()) { x.join(); } } }
          } join{pool};
          pool.reserve(workers - 1);
          for (unsigned w = 1; w < workers; ++w) { pool.emplace_back(worker, w); }
          worker(0);   // the driver is one of the T threads, not a spectator
        }
        for (const auto& f : failed) { if (f) { std::rethrow_exception(f); } }
        for (const auto& r : unit_res)
        {
          const std::size_t take = std::min(r.size(), RESID_CAP - resid.size());
          resid.insert(resid.end(), r.begin(), r.begin() + static_cast<std::ptrdiff_t>(take));
          if (resid.size() >= RESID_CAP) { break; }
        }
        out.stats_.units += runs;
        out.stats_.threads = std::max(out.stats_.threads, workers);
      }
      out.stats_.match_s += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_match).count();
    }
    // Before nth_element, which reorders it: the gate is on the list in the
    // order it was kept, not merely on the median it yields.
    {
      Fnv h;
      for (const double r : resid) { h.add(r); }
      out.stats_.resid_kept = resid.size();
      out.stats_.resid_hash = h.h;
    }
    if (observed_ppm_median != nullptr)
    {
      if (resid.empty()) { *observed_ppm_median = std::numeric_limits<double>::quiet_NaN(); }
      else
      {
        const std::size_t h = resid.size() / 2;
        std::nth_element(resid.begin(), resid.begin() + h, resid.end());
        *observed_ppm_median = resid[h];
      }
    }
    return out;
  }
} // namespace ODIA
