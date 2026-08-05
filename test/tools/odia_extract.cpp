// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Extract chromatograms for a library from a run, and report what it cost.
// Decode and match are timed separately: with the current mzPeak reader the
// first dominates by three orders of magnitude, and conflating them would hide
// which one any future change actually improved.

#include <odia/ChromatogramExtractor.h>
#include <odia/DIANNLibraryFile.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

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

    // Every option the matching depends on is reachable from the environment,
    // because the brute force below is only worth anything if it can be run
    // against the settings the workflow actually uses. A verification that can
    // only check the defaults stops covering the code the moment a default is
    // not what a run uses -- which is how it came to model an extractor that
    // no longer existed.
    if (const char* th = std::getenv("ODIA_EXTRACT_THREADS")) { opt.threads = std::atoi(th); }
    if (const char* v = std::getenv("ODIA_EXTRACT_PPM")) { opt.fragment_ppm = std::atof(v); }
    if (const char* v = std::getenv("ODIA_EXTRACT_PPM_OFFSET"))
    { opt.fragment_ppm_offset = std::atof(v); }
    if (const char* v = std::getenv("ODIA_EXTRACT_PPM_LOG_SLOPE"))
    { opt.fragment_ppm_log_slope = std::atof(v); }
    if (const char* v = std::getenv("ODIA_EXTRACT_PPM_SLOPE_PER_1000"))
    { opt.fragment_ppm_slope_per_1000 = std::atof(v); }
    if (const char* v = std::getenv("ODIA_EXTRACT_PRECURSOR_IM"))
    { opt.precursor_im_window = std::atof(v); }
    if (std::getenv("ODIA_EXTRACT_NO_IM")) { opt.use_ion_mobility = false; }
    if (const char* v = std::getenv("ODIA_EXTRACT_AGGREGATE"))
    {
      opt.aggregate = std::string(v) == "max"
                        ? ODIA::ChromatogramExtractor::Options::Aggregate::Max
                        : ODIA::ChromatogramExtractor::Options::Aggregate::Sum;
    }
    std::printf("options: %s, %.1f ppm about %+.2f ppm (log slope %.2f, slope/1000 %.2f), "
                "band %s, precursor 1/K0 window %.3f\n",
                opt.aggregate == ODIA::ChromatogramExtractor::Options::Aggregate::Sum
                  ? "sum" : "max",
                opt.fragment_ppm, opt.fragment_ppm_offset, opt.fragment_ppm_log_slope,
                opt.fragment_ppm_slope_per_1000, opt.use_ion_mobility ? "on" : "off",
                opt.precursor_im_window);

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
    std::printf("  several windows: %zu precursors (%.2f%%), extracted from the nearest centre\n",
                xics.precursors_in_several_windows,
                stats.precursors ? 100.0 * double(xics.precursors_in_several_windows) /
                                     double(stats.precursors) : 0.0);
    std::printf("  memory:         %.1f MiB\n", xics.footprintBytes() / 1048576.0);
    std::printf("  decode:         %.2f s   match: %.2f s   total: %.2f s\n",
                stats.decode_seconds, stats.match_seconds, sec);
    if (stats.spectra_read)
    {
      std::printf("  per spectrum:   %.2f ms decode, %.3f ms match\n",
                  1000.0 * stats.decode_seconds / double(stats.spectra_read),
                  1000.0 * stats.match_seconds / double(stats.spectra_read));
    }

    // Verify the extraction against a brute-force scan of the same spectra.
    //
    // The index inverts the match and buckets on log m/z; a bucket off by one,
    // or a live-cycle range off by one, produces a plausible chromatogram that
    // is quietly wrong, and recomputing a sample the slow, obvious way is the
    // only thing that catches that. But only while the slow way models the
    // fast one: this loop went on taking a MAXIMUM after the default became
    // Sum, and applied neither mobility test and no mass calibration, so it
    // agreed with an extractor that had not existed for four commits. A
    // brute force that models something else is worse than none, because it
    // reports agreement.
    //
    // What it must therefore reproduce, in the same order the extractor
    // applies it: the one window a precursor is extracted from, the corrected
    // target m/z, the tolerance about THAT target, the frame's mobility band,
    // the per-precursor mobility window, and the aggregate.
    if (std::getenv("ODIA_EXTRACT_VERIFY"))
    {
      const auto& p = library.precursors();
      const auto& t = library.transitions();
      const auto& info = run->spectra();
      const auto& windows = run->windows();
      const bool sum_peaks =
        opt.aggregate == ODIA::ChromatogramExtractor::Options::Aggregate::Sum;
      std::size_t limit = 400;
      if (const char* v = std::getenv("ODIA_EXTRACT_VERIFY_MAX")) { limit = std::strtoul(v, nullptr, 10); }

      const auto inRange = [&](const ODIA::SpectrumInfo& s) {
        return !(opt.rt_high > opt.rt_low) ||
               (s.retention_time >= opt.rt_low && s.retention_time < opt.rt_high);
      };

      // The window a precursor is extracted from: the one whose centre it is
      // nearest, of those covering it. Stated here rather than asked of the
      // extractor, so that a change of policy shows up as a disagreement
      // instead of being copied into both sides of the comparison.
      const auto windowOf = [&](double mz) {
        std::size_t best = windows.size();
        double best_offset = 0.0;
        for (std::size_t w = 0; w < windows.size(); ++w)
        {
          if (!windows[w].contains(mz)) { continue; }
          const double d = std::fabs(mz - windows[w].centre());
          if (best == windows.size() || d < best_offset) { best_offset = d; best = w; }
        }
        return best;
      };

      // Which spectra the check will look at, decoded ONCE. The reader costs
      // ~6 ms a spectrum on 12_80, so decoding per transition turned a 4 s
      // check into a 2-minute one and was why this never ran anywhere.
      std::vector<char> wanted(info.size(), 0);
      std::size_t planned = 0;
      for (std::size_t i = 0; i < stats.precursors && planned < limit; ++i)
      {
        const std::size_t w = windowOf(ODIA::fromFixed(p.mz[i]));
        if (w == windows.size()) { continue; }
        bool any = false;
        for (std::uint32_t k = 0; k < p.transition_count[i] && planned < limit; ++k)
        {
          if (xics.count[p.transition_begin[i] + k] == 0) { continue; }
          ++planned; any = true;
        }
        if (!any) { continue; }
        for (std::size_t si = 0; si < info.size(); ++si)
        {
          if (!inRange(info[si])) { continue; }
          if (std::fabs(info[si].window.mz_low - windows[w].mz_low) > 1e-6 ||
              std::fabs(info[si].window.mz_high - windows[w].mz_high) > 1e-6) { continue; }
          wanted[si] = 1;
        }
      }

      std::vector<ODIA::SpectrumPeaks> cache(info.size());
      std::size_t decoded = 0;
      {
        constexpr std::size_t BLOCK = 512;
        std::vector<ODIA::SpectrumPeaks> block;
        for (std::size_t b = 0; b < info.size(); b += BLOCK)
        {
          const std::size_t e = std::min(b + BLOCK, info.size());
          bool any = false;
          for (std::size_t si = b; si < e; ++si) { any = any || wanted[si]; }
          if (!any) { continue; }
          run->peaks(b, e, block);
          for (std::size_t si = b; si < e; ++si)
          {
            if (wanted[si]) { cache[si] = std::move(block[si - b]); ++decoded; }
          }
        }
      }

      std::size_t checked = 0, bad = 0, nonzero = 0, compared = 0;
      double worst = 0.0;

      for (std::size_t i = 0; i < stats.precursors && checked < limit; ++i)
      {
        const double pmz = ODIA::fromFixed(p.mz[i]);
        const std::size_t w = windowOf(pmz);
        if (w == windows.size()) { continue; }
        const double want_im = double(p.im[i]);

        for (std::uint32_t k = 0; k < p.transition_count[i] && checked < limit; ++k)
        {
          const std::size_t j = p.transition_begin[i] + k;
          if (xics.count[j] == 0) { continue; }

          // The corrected target, exactly as the extractor computes it once
          // per transition: the offset is applied to the TARGET, and the
          // tolerance is a fraction of that target and not of the peak.
          const double theoretical = ODIA::fromFixed(t.product_mz[j]);
          double correction_ppm = opt.fragment_ppm_offset;
          if (opt.fragment_ppm_log_slope != 0.0 && opt.fragment_ppm_ref_mz > 0.0)
          {
            correction_ppm +=
              opt.fragment_ppm_log_slope * std::log(theoretical / opt.fragment_ppm_ref_mz);
          }
          if (opt.fragment_ppm_slope_per_1000 != 0.0)
          {
            correction_ppm += opt.fragment_ppm_slope_per_1000 *
                              (theoretical - opt.fragment_ppm_ref_mz) / 1000.0;
          }
          const double target = theoretical * (1.0 + correction_ppm * 1e-6);
          const double tol = target * opt.fragment_ppm * 1e-6;

          std::vector<float> want;
          for (std::size_t si = 0; si < info.size(); ++si)
          {
            const auto& s = info[si];
            if (!wanted[si]) { continue; }
            if (std::fabs(s.window.mz_low - windows[w].mz_low) > 1e-6 ||
                std::fabs(s.window.mz_high - windows[w].mz_high) > 1e-6) { continue; }

            const auto& pk = cache[si];
            const bool has_im = pk.hasIonMobility();
            const bool use_band = opt.use_ion_mobility && has_im;
            float acc = 0.0f;
            for (std::size_t q = 0; q < pk.size(); ++q)
            {
              if (pk.mz[q] < target - tol || pk.mz[q] > target + tol) { continue; }
              const double peak_im = has_im ? double(pk.ion_mobility[q])
                                            : std::numeric_limits<double>::quiet_NaN();
              // Half-open at the top, as the extractor's band is.
              if (use_band &&
                  (peak_im < s.window.im_low || peak_im >= s.window.im_high)) { continue; }
              if (opt.precursor_im_window > 0.0 && !std::isnan(peak_im) &&
                  !std::isnan(want_im) &&
                  std::fabs(peak_im - want_im) > opt.precursor_im_window) { continue; }
              if (sum_peaks) { acc += pk.intensity[q]; }
              else { acc = std::max(acc, pk.intensity[q]); }
            }
            want.push_back(acc);
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
              ++compared;
              if (want[q] > 0.0f) { ++nonzero; }
              const double got = double(xics.intensity[xics.begin[j] + q]);
              const double d = std::fabs(double(want[q]) - got);
              if (d > worst) { worst = d; }
              if (d > 1e-3)
              {
                if (!reported && bad < 4)
                {
                  std::printf("  DIFF transition %zu (product %.5f, %u points) at %zu: "
                              "brute force %.4g, extractor %.4g\n",
                              j, theoretical, xics.count[j], q, double(want[q]), got);
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
      std::printf("verify: %zu transitions brute-forced over %zu decoded spectra, "
                  "%zu of %zu points non-zero, %zu disagree, worst |diff| %.6g\n",
                  checked, decoded, nonzero, compared, bad, worst);
      // An all-zero comparison agrees with anything. Most of a slice IS zeros
      // -- a precursor elutes for seconds in a run of minutes -- so the check
      // is only worth running where at least some of it carries signal, and
      // saying so is the difference between a test and a formality.
      if (checked == 0 || nonzero == 0)
      {
        std::printf("verify: nothing with signal was checked; "
                    "widen the range or raise ODIA_EXTRACT_VERIFY_MAX\n");
        return 1;
      }
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
