// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Ms1Traces::build at several thread counts against a scripted run: the matrix,
// the retained residual list and its median must be IDENTICAL to the serial
// build, cell for cell and bit for bit, and the decoder must be asked for the
// same frame ranges as the serial loop asks for. Not close -- ODIA is deterministic, and "close" is how a
// race looks while it is still rare.
//
// The run is built to break the two things a frame-parallel build can get
// wrong:
//
//   cells     every frame writes many rows, ions sit in ~5 adjacent frames,
//             and a mobility gate rejects some peaks, so a unit writing the
//             wrong column or racing a neighbour changes a cell
//   residual  ~3.2M matches against a 2M-sample cap, crossing it inside a
//             block and inside a unit, with a residual that drifts with the
//             frame -- so a prefix rebuilt out of order, or cut per unit
//             instead of in total, changes the list hash and moves the median
//   decode    the match block grows with the threads but the decode calls
//             must not: every run records the ranges it was asked for
//
// Frame counts are chosen so blocks are not multiples of the 16-frame unit and
// the last block is short. No data file, a second or two.
//
// Usage: odia_ms1_parallel

#include <odia/Ms1Traces.h>
#include <odia/SpectrumSource.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace
{
  const float NA = std::numeric_limits<float>::quiet_NaN();

  /// MS1 only; the MS2 side is never touched by Ms1Traces.
  class ScriptedMs1Run : public ODIA::SpectrumSource
  {
  public:
    std::vector<ODIA::SpectrumInfo> ms1_info;
    std::vector<ODIA::SpectrumPeaks> ms1_peaks;
    std::vector<std::pair<std::size_t, std::size_t>> asked;   ///< ms1Peaks ranges, in call order

    const std::vector<ODIA::SpectrumInfo>& spectra() const override { return none_; }
    const std::vector<ODIA::IsolationWindow>& windows() const override { return windows_; }
    const std::vector<ODIA::SpectrumInfo>& ms1Spectra() const override { return ms1_info; }
    std::string describe() const override { return "scripted MS1 run"; }

    void ms1Peaks(std::size_t begin, std::size_t end,
                  std::vector<ODIA::SpectrumPeaks>& out) override
    {
      asked.emplace_back(begin, end);
      out.assign(ms1_peaks.begin() + static_cast<std::ptrdiff_t>(begin),
                 ms1_peaks.begin() + static_cast<std::ptrdiff_t>(end));
    }
    void peaks(std::size_t, std::size_t, std::vector<ODIA::SpectrumPeaks>& out) override
    { out.clear(); }

  private:
    std::vector<ODIA::SpectrumInfo> none_;
    std::vector<ODIA::IsolationWindow> windows_;
  };

  /// Deterministic, platform-independent pseudo-random stream.
  struct Lcg
  {
    std::uint64_t s;
    std::uint32_t next() { s = s * 6364136223846793005ull + 1442695040888963407ull; return std::uint32_t(s >> 33); }
    double unit() { return next() / 2147483648.0; }   // [0, 1)
  };

  int failures = 0;
  void check(bool ok, const std::string& what)
  {
    std::printf("%s  %s\n", ok ? "ok  " : "FAIL", what.c_str());
    if (!ok) { ++failures; }
  }

  bool sameBits(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }
}

int main()
{
  constexpr std::size_t PRECURSORS = 20000;
  constexpr std::size_t FRAMES = 403;           // 6 x 64 + 19: no block is a multiple of 16
  constexpr std::size_t PEAKS_PER_FRAME = 12000; // 4.8M peaks, ~3.2M pass the gate
  constexpr double PPM = 10.0, IM_WINDOW = 0.05;

  ODIA::Library lib;
  {
    auto& p = lib.precursors();
    Lcg r{42};
    for (std::size_t i = 0; i < PRECURSORS; ++i)
    {
      // ~40 mTh apart: one target per +/-10 ppm window, a few where two crowd.
      p.mz.push_back(ODIA::toFixed(400.0 + 800.0 * r.unit()));
      p.irt.push_back(0.0f);
      // One in 20 has no library 1/K0: the ungated path must parallelise too.
      p.im.push_back(i % 20 == 0 ? NA : static_cast<float>(0.7 + 0.6 * r.unit()));
      p.ccs.push_back(NA);
      p.charge.push_back(static_cast<std::uint8_t>(2 + i % 3));
      p.decoy.push_back(0);
      p.modified_sequence.push_back(0);
      p.protein_group.push_back(0);
      p.transition_begin.push_back(0);
      p.transition_count.push_back(0);
    }
  }

  ScriptedMs1Run run;
  {
    const auto& p = lib.precursors();
    Lcg r{7};
    for (std::size_t f = 0; f < FRAMES; ++f)
    {
      ODIA::SpectrumInfo s;
      s.index = f;
      s.retention_time = 1.8 * static_cast<double>(f);
      s.ms_level = 1;
      run.ms1_info.push_back(s);

      ODIA::SpectrumPeaks sp;
      // A residual that drifts with the frame, so the median depends on WHICH
      // 2M samples are kept and in what order they were cut.
      const double drift_ppm = 6.0 * std::sin(0.05 * static_cast<double>(f));
      for (std::size_t k = 0; k < PEAKS_PER_FRAME; ++k)
      {
        // Ions elute over ~5 adjacent frames: precursor chosen from a band
        // that slides with the frame, so neighbouring frames hit the same rows.
        const std::size_t i = (f * 37 + r.next() % 4000) % PRECURSORS;
        const double target = ODIA::fromFixed(p.mz[i]);
        const double jitter = (r.unit() - 0.5) * 8.0;
        sp.mz.push_back(target * (1.0 + (drift_ppm + jitter) * 1e-6));
        sp.intensity.push_back(static_cast<float>(1.0 + 1e5 * r.unit()));
        // Half the peaks inside the gate, half up to 3 widths outside.
        const double lib_im = std::isnan(p.im[i]) ? 1.0 : p.im[i];
        const double dim = (r.unit() < 0.5 ? 0.9 : 3.0) * IM_WINDOW * (2.0 * r.unit() - 1.0);
        sp.ion_mobility.push_back(static_cast<float>(lib_im + dim));
      }
      run.ms1_peaks.push_back(std::move(sp));
    }
  }

  double serial_median = 0.0;
  run.asked.clear();
  const auto serial = ODIA::Ms1Traces::build(lib, run, PPM, IM_WINDOW, 0.0, &serial_median);
  const auto serial_asked = run.asked;
  const auto& ss = serial.buildStats();
  std::size_t live = 0;
  for (std::size_t i = 0; i < serial.precursors(); ++i)
    for (std::size_t b = 0; b < serial.bins(); ++b) { if (serial.at(i, b) > 0.0f) { ++live; } }
  // A vacuous identity (an empty matrix equals an empty matrix) is not a test.
  check(serial.bins() == FRAMES && serial.precursors() == PRECURSORS,
        "serial build has the scripted shape");
  check(live > 100000, "serial build has " + std::to_string(live) + " non-zero cells");
  check(std::isfinite(serial_median) && serial_median != 0.0,
        "serial residual median is finite and nonzero (" + std::to_string(serial_median) + ")");
  check(ss.resid_kept == (std::size_t(1) << 21),
        "serial build kept exactly the 2M-sample cap (" + std::to_string(ss.resid_kept) + ")");
  check(serial_asked.size() == (FRAMES + 63) / 64 && serial_asked.front().second == 64,
        "serial build decodes in 64-frame calls (" + std::to_string(serial_asked.size()) + ")");

  // Was the cap actually crossed? Re-derive the match count independently of
  // the build: if fewer than 2M matched, the residual half of this test is
  // vacuous and must say so.
  {
    std::size_t matches = 0;
    const auto& p = lib.precursors();
    std::vector<std::pair<double, float>> t;
    for (std::size_t i = 0; i < PRECURSORS; ++i) { t.push_back({ODIA::fromFixed(p.mz[i]), p.im[i]}); }
    std::sort(t.begin(), t.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    for (const auto& sp : run.ms1_peaks)
      for (std::size_t k = 0; k < sp.mz.size(); ++k)
      {
        const double m = sp.mz[k], tol = m * PPM * 1e-6;
        auto it = std::lower_bound(t.begin(), t.end(), m - tol,
                                   [](const auto& a, double v) { return a.first < v; });
        for (; it != t.end() && it->first <= m + tol; ++it)
        {
          if (std::abs(it->first - m) > it->first * PPM * 1e-6) { continue; }
          if (!std::isnan(it->second) && !(std::abs(sp.ion_mobility[k] - it->second) <= IM_WINDOW))
          { continue; }
          ++matches;
        }
      }
    check(matches > (std::size_t(1) << 21) + 500000,
          "residual cap is crossed inside the run (" + std::to_string(matches) + " matches)");
  }

  // One decode's decoded bytes, as the build counts them (capacity, as the
  // scripted source assigns exact sizes): the unit of the block cap.
  const std::size_t frame_bytes = PEAKS_PER_FRAME * (sizeof(double) + 2 * sizeof(float));
  const std::size_t decode_bytes = 64 * frame_bytes;

  auto against_serial = [&](unsigned threads, std::size_t block_budget, const std::string& tag)
  {
    double median = 0.0;
    run.asked.clear();
    const auto par = ODIA::Ms1Traces::build(lib, run, PPM, IM_WINDOW, 0.0, &median,
                                            0.0, nullptr, nullptr, threads, block_budget);
    const auto& ps = par.buildStats();
    std::size_t diff = 0;
    for (std::size_t i = 0; i < PRECURSORS; ++i)
      for (std::size_t b = 0; b < FRAMES; ++b)
      {
        const float x = serial.at(i, b), y = par.at(i, b);
        if (std::memcmp(&x, &y, sizeof x) != 0) { ++diff; }
      }
    check(par.bins() == FRAMES && par.precursors() == PRECURSORS, tag + "same shape");
    check(diff == 0, tag + std::to_string(diff) + " cells differ from serial");
    check(par.checksum() == serial.checksum(), tag + "checksum equals serial");
    check(ps.resid_kept == ss.resid_kept && ps.resid_hash == ss.resid_hash,
          tag + "residual list identical (" + std::to_string(ps.resid_kept) + " kept, same hash)");
    check(sameBits(median, serial_median),
          tag + "median residual bit-identical (" + std::to_string(median) + ")");
    check(run.asked == serial_asked, tag + "decoder asked for the serial ranges");
    if (threads >= 2)
    {
      check(ps.threads > 1, tag + "actually matched on " + std::to_string(ps.threads) + " threads");
    }
    // The memory contract (review round 2): beyond its first decode a match
    // block holds at most the cap, residual buffers included, whatever the
    // threads; decoded bytes alive at once are at most a block plus one decode.
    const std::size_t cap = std::min(block_budget, ODIA::Ms1Traces::MAX_BLOCK_BYTES);
    check(ps.block_cap_bytes == cap, tag + "cap in force is min(budget, MAX_BLOCK_BYTES)");
    check(ps.max_block_bytes <= decode_bytes + cap,
          tag + "largest block " + std::to_string(ps.max_block_bytes) + " B <= one decode + cap");
    check(ps.max_held_bytes <= ps.max_block_bytes + decode_bytes,
          tag + "held at once " + std::to_string(ps.max_held_bytes) + " B <= block + one decode");
    std::printf("      %s\n", par.describeBuild().c_str());
    return ps;
  };

  for (const unsigned threads : {0u, 2u, 3u, 8u, 13u, 48u})
  {
    against_serial(threads, SIZE_MAX, "-threads " + std::to_string(threads) + ": ");
  }

  // THE CAP. Uncapped, the match block was 16 x threads frames -- the whole run
  // once threads reach frames/16 -- and nothing charged it to a budget. Driven
  // here at budgets that force every path: 0 (one decode per block, the serial
  // floor), one that refuses decodes while the residual prefix is open (four
  // 16 MB unit buffers per decode) and admits a second decode after it closes,
  // and two decodes' worth. Output must not move; blocks must.
  {
    const std::size_t decodes = (FRAMES + 63) / 64;
    const auto z = against_serial(48, 0, "-threads 48, block budget 0: ");
    check(z.match_blocks == decodes && z.max_block_bytes <= decode_bytes,
          "budget 0: one decode per match block (" + std::to_string(z.match_blocks) + " blocks)");
    check(z.carried == 0 && z.max_held_bytes == z.max_block_bytes,
          "budget 0: no decode carried, so no more held than the serial path's one decode");
    const auto one = against_serial(48, decode_bytes + decode_bytes / 2,
                                    "-threads 48, block budget 1.5 decodes: ");
    check(one.match_blocks > 1 && one.match_blocks < decodes,
          "1.5 decodes: carried decodes split the run into " +
          std::to_string(one.match_blocks) + " blocks (uncapped: 1)");
    // While the prefix is open each decode brings four 16 MiB unit buffers,
    // more than the cap: only the first decode of a block (the floor) gets in.
    check(one.max_resid_bytes > 0 && one.max_resid_bytes <= 4 * (std::size_t(1) << 21) * sizeof(double),
          "1.5 decodes: residual buffers no more than one decode's four units (" +
          std::to_string(one.max_resid_bytes) + " B)");
    const auto two = against_serial(13, 2 * decode_bytes, "-threads 13, block budget 2 decodes: ");
    check(two.match_blocks >= 3, "2 decodes at -threads 13: " +
          std::to_string(two.match_blocks) + " blocks");
  }

  // THE CARRY. A decode is predicted from the largest so far, so it is carried
  // only when it outgrows every earlier one: frames 64-127 are four times as
  // dense as the rest. No residual requested (the prefix would stay open and
  // keep every block at one decode). Blocks [0,64) | carried [64,128) |
  // [128,192) | [192,200): the dense decode is refused after being made.
  {
    ScriptedMs1Run vary;
    Lcg r{11};
    const auto& p = lib.precursors();
    for (std::size_t f = 0; f < 200; ++f)
    {
      ODIA::SpectrumInfo si;
      si.index = f;
      si.retention_time = 1.8 * static_cast<double>(f);
      si.ms_level = 1;
      vary.ms1_info.push_back(si);
      ODIA::SpectrumPeaks sp;
      const std::size_t npk = (f >= 64 && f < 128) ? 8000 : 2000;
      for (std::size_t k = 0; k < npk; ++k)
      {
        const std::size_t i = (f * 37 + r.next() % 4000) % PRECURSORS;
        const double target = ODIA::fromFixed(p.mz[i]);
        sp.mz.push_back(target * (1.0 + (r.unit() - 0.5) * 8.0e-6));
        sp.intensity.push_back(static_cast<float>(1.0 + 1e5 * r.unit()));
        const double lib_im = std::isnan(p.im[i]) ? 1.0 : p.im[i];
        sp.ion_mobility.push_back(static_cast<float>(lib_im + IM_WINDOW * (2.0 * r.unit() - 1.0)));
      }
      vary.ms1_peaks.push_back(std::move(sp));
    }
    const std::size_t small = 64 * 2000 * (sizeof(double) + 2 * sizeof(float));
    const auto ref = ODIA::Ms1Traces::build(lib, vary, PPM, IM_WINDOW, 0.0, nullptr, 0.0,
                                            nullptr, nullptr, 1);
    const auto serial_vary_asked = vary.asked;
    vary.asked.clear();
    const auto cut = ODIA::Ms1Traces::build(lib, vary, PPM, IM_WINDOW, 0.0, nullptr, 0.0,
                                            nullptr, nullptr, 48, small + small / 2);
    const auto& cs = cut.buildStats();
    std::printf("      %s\n", cut.describeBuild().c_str());
    check(cut.checksum() == ref.checksum(), "carry: checksum equals serial");
    check(vary.asked == serial_vary_asked, "carry: decoder asked for the serial ranges");
    check(cs.carried == 1, "carry: exactly the dense decode was carried (" +
          std::to_string(cs.carried) + ")");
    check(cs.match_blocks == 4, "carry: 4 match blocks (" + std::to_string(cs.match_blocks) + ")");
    check(cs.max_held_bytes == small + 4 * small && cs.max_block_bytes == 4 * small,
          "carry: held at once = refused block + carried decode, never more");
  }

  // The -out_ms1_iso path: a keep mask (compact rows, slot != library index)
  // with no residual requested. Monoisotopic, because the scripted peaks sit on
  // the monoisotopic targets and an M+1 build would compare two empty matrices.
  {
    std::vector<std::uint8_t> keep(PRECURSORS, 0);
    for (std::size_t i = 0; i < PRECURSORS; i += 3) { keep[i] = 1; }
    std::vector<std::uint32_t> rows1, rows8;
    const auto a = ODIA::Ms1Traces::build(lib, run, PPM, IM_WINDOW, 0.0, nullptr,
                                          0.0, &keep, &rows1, 1);
    const auto b = ODIA::Ms1Traces::build(lib, run, PPM, IM_WINDOW, 0.0, nullptr,
                                          0.0, &keep, &rows8, 8);
    std::size_t masked_live = 0;
    for (std::size_t i = 0; i < a.precursors(); ++i)
      for (std::size_t f = 0; f < a.bins(); ++f) { if (a.at(i, f) > 0.0f) { ++masked_live; } }
    check(rows1 == rows8 && !rows1.empty(), "masked build: same kept rows");
    check(masked_live > 10000, "masked build has " + std::to_string(masked_live) + " non-zero cells");
    check(a.checksum() == b.checksum(), "masked build: -threads 8 checksum equals serial");
  }

  std::printf("%s\n", failures ? "FAILED" : "all passed");
  return failures ? 1 : 0;
}
