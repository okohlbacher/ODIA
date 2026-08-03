// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// What SpectrumSource makes of a run: the window scheme and the RT range from
// metadata alone, then a small timed peak read. Metadata is nearly free and
// peaks are not, so the two are reported separately.

#include <odia/SpectrumSource.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::fprintf(stderr, "usage: odia_run_info <run> [spectra-to-decode]\n");
    return 2;
  }
  try
  {
    const auto t0 = std::chrono::steady_clock::now();
    auto run = ODIA::openRun(argv[1]);
    const double meta_s = std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - t0).count();

    const auto& s = run->spectra();
    std::printf("%s\n", run->describe().c_str());
    std::printf("metadata read in %.2f s\n", meta_s);
    if (!s.empty())
    {
      std::printf("RT range: %.1f .. %.1f s\n",
                  s.front().retention_time, s.back().retention_time);
    }
    std::printf("isolation windows: %zu\n", run->windows().size());
    std::size_t shown = 0;
    for (const auto& w : run->windows())
    {
      if (shown++ >= 5) { break; }
      std::printf("  %.3f .. %.3f  (width %.3f)", w.mz_low, w.mz_high, w.width());
      if (w.im_low > -1e30) { std::printf("  1/K0 %.3f .. %.3f", w.im_low, w.im_high); }
      std::printf("\n");
    }

    const std::size_t n = argc > 2 ? std::strtoul(argv[2], nullptr, 10) : 0;
    if (n > 0 && !s.empty())
    {
      std::vector<ODIA::SpectrumPeaks> block;
      const auto t1 = std::chrono::steady_clock::now();
      run->peaks(0, std::min(n, s.size()), block);
      const double sec = std::chrono::duration<double>(
                           std::chrono::steady_clock::now() - t1).count();
      std::size_t peaks = 0, with_im = 0;
      for (const auto& p : block)
      {
        peaks += p.size();
        if (p.hasIonMobility()) { ++with_im; }
      }
      std::printf("decoded %zu spectra, %zu peaks in %.2f s (%.1f ms/spectrum); "
                  "%zu carry ion mobility\n",
                  block.size(), peaks, sec, 1000.0 * sec / double(block.size()), with_im);
    }
  }
  catch (const std::exception& e)
  {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  return 0;
}
