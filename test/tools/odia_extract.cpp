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
#include <cmath>
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

    if (const char* th = std::getenv("ODIA_EXTRACT_THREADS")) { opt.threads = std::atoi(th); }

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

    // Verify the m/z index against a brute-force scan of the same spectra.
    //
    // The index inverts the match and buckets on log m/z; a bucket off by one,
    // or a live-cycle range off by one, produces a plausible chromatogram that
    // is quietly wrong. Recomputing a sample the slow, obvious way is the only
    // thing that catches that.
    if (std::getenv("ODIA_EXTRACT_VERIFY"))
    {
      const auto& p = library.precursors();
      const auto& t = library.transitions();
      const auto& info = run->spectra();
      std::size_t checked = 0, bad = 0;
      double worst = 0.0;

      for (std::size_t i = 0; i < stats.precursors && checked < 400; ++i)
      {
        const double pmz = ODIA::fromFixed(p.mz[i]);
        for (std::uint32_t k = 0; k < p.transition_count[i] && checked < 400; ++k)
        {
          const std::size_t j = p.transition_begin[i] + k;
          if (xics.count[j] == 0) { continue; }
          const double fmz = ODIA::fromFixed(t.product_mz[j]);
          const double tol = fmz * opt.fragment_ppm * 1e-6;

          // A precursor can sit in several overlapping windows, and each is a
          // separate measurement. The extractor concatenates them in window
          // order, so the brute force must too -- interleaving them by
          // retention time compares two different orderings and reports a
          // disagreement that is only a convention.
          std::vector<float> want;
          for (const auto& win : run->windows())
          {
            if (!win.contains(pmz)) { continue; }
            for (std::size_t si = 0; si < info.size(); ++si)
            {
              const auto& s = info[si];
              if (std::fabs(s.window.mz_low - win.mz_low) > 1e-6 ||
                  std::fabs(s.window.mz_high - win.mz_high) > 1e-6) { continue; }
              if (opt.rt_high > opt.rt_low &&
                  (s.retention_time < opt.rt_low || s.retention_time >= opt.rt_high)) { continue; }
              ODIA::SpectrumPeaks pk;
              run->peaks(si, pk);
              float best = 0.0f;
              for (std::size_t q = 0; q < pk.size(); ++q)
              {
                if (pk.mz[q] >= fmz - tol && pk.mz[q] <= fmz + tol)
                { best = std::max(best, pk.intensity[q]); }
              }
              want.push_back(best);
            }
          }
          if (want.size() != xics.count[j])
          {
            std::printf("  MISMATCH transition %zu: %zu points, brute force says %zu\n",
                        j, std::size_t(xics.count[j]), want.size());
            ++bad;
          }
          else
          {
            bool reported = false;
            for (std::size_t q = 0; q < want.size(); ++q)
            {
              const double got = double(xics.intensity[xics.begin[j] + q]);
              const double d = std::fabs(double(want[q]) - got);
              if (d > worst) { worst = d; }
              if (d > 1e-3)
              {
                if (!reported && bad < 4)
                {
                  std::printf("  DIFF transition %zu (product %.5f, %u points) at %zu: "
                              "brute force %.4g, extractor %.4g\n",
                              j, fmz, xics.count[j], q, double(want[q]), got);
                  reported = true;
                }
                ++bad;
                break;
              }
            }
          }
          ++checked;
        }
      }
      std::printf("verify: %zu transitions brute-forced, %zu disagree, worst |diff| %.6g\n",
                  checked, bad, worst);
      if (bad) { return 1; }
    }
  }
  catch (const std::exception& e)
  {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  return 0;
}
