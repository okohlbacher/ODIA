// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Extract chromatograms for a library from a run, and report what it cost.
// Decode and match are timed separately: with the current mzPeak reader the
// first dominates by three orders of magnitude, and conflating them would hide
// which one any future change actually improved.

#include <odia/ChromatogramExtractor.h>
#include <odia/DIANNLibraryFile.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>

int main(int argc, char** argv)
{
  if (argc < 3)
  {
    std::fprintf(stderr,
                 "usage: odia_extract <library.tsv> <run> [max_precursors] [rt_low rt_high]\n");
    return 2;
  }
  try
  {
    ODIA::Library library;
    ODIA::DIANNLibraryFile::load(argv[1], library);
    std::printf("library: %zu precursors, %zu transitions\n",
                library.precursorCount(), library.transitionCount());

    auto run = ODIA::openRun(argv[2]);
    std::printf("%s\n", run->describe().c_str());

    ODIA::ChromatogramExtractor::Options opt;
    if (argc > 3) { opt.max_precursors = std::strtoul(argv[3], nullptr, 10); }
    if (argc > 5)
    {
      opt.rt_low = std::atof(argv[4]);
      opt.rt_high = std::atof(argv[5]);
    }

    ODIA::ChromatogramExtractor::Stats stats;
    const auto t0 = std::chrono::steady_clock::now();
    auto xics = ODIA::ChromatogramExtractor::extract(library, *run, opt, &stats);
    const double sec = std::chrono::duration<double>(
                         std::chrono::steady_clock::now() - t0).count();

    std::printf("extracted %zu transitions of %zu precursors\n",
                stats.transitions, stats.precursors);
    std::printf("  spectra read:   %zu\n", stats.spectra_read);
    std::printf("  points:         %zu (%zu non-zero, %.1f%%)\n", stats.points,
                stats.nonzero_points,
                stats.points ? 100.0 * double(stats.nonzero_points) / double(stats.points) : 0.0);
    std::printf("  without window: %zu precursors\n", xics.precursors_without_window);
    std::printf("  memory:         %.1f MiB\n", xics.footprintBytes() / 1048576.0);
    std::printf("  decode:         %.2f s   match: %.2f s   total: %.2f s\n",
                stats.decode_seconds, stats.match_seconds, sec);
    if (stats.spectra_read)
    {
      std::printf("  per spectrum:   %.2f ms decode, %.3f ms match\n",
                  1000.0 * stats.decode_seconds / double(stats.spectra_read),
                  1000.0 * stats.match_seconds / double(stats.spectra_read));
    }
  }
  catch (const std::exception& e)
  {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  return 0;
}
