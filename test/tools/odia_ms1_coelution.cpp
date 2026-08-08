// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

/// Does the MS1 precursor trace rise and fall with the MS2 fragments?
///
/// The last untested statistic, and the only MS1 quantity with a reason to work.
/// Six measurements on 2026-08-07/08 established that every PRESENCE statistic
/// saturates on this data -- fragment depth, mobility-sliced depth,
/// qualifying_spectra, total_matches, MS1 isotope depth all give 0.5x to 1.0x
/// enrichment -- because each is a maximum or count over 10^3-10^4 spectra and
/// so is an extreme-value statistic over thousands of draws.
///
/// SHAPE is the alternative, and the project has one large piece of evidence for
/// it: replacing the amplitude peak-picker with co-elution detection fixed RT,
/// FDR, recovery and memory at once. Co-elution cannot be satisfied by an
/// accumulation of unrelated coincidences, because the coincidences would have
/// to arrive in the right ORDER.
///
/// So: build the MS1 monoisotopic trace and the MS2 top-fragment trace on one
/// common time grid, correlate them, and ask whether that separates DIA-NN's
/// confident set from the rest. If it does not, MS1 is not the lever either.
///
/// The grid is the MS1 acquisition times (S08: 1,343 points, ~1.8 s apart). MS2
/// spectra are accumulated into their nearest MS1 bin, which is the honest
/// alignment: an MS2 frame acquired between two survey scans belongs to the
/// nearer one, and a peptide's peak is ~20-30 s wide, i.e. ~15 bins.

#include <odia/DIANNLibraryFile.h>
#include <odia/Library.h>
#include <odia/SpectrumSource.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace
{
  constexpr std::size_t TOP_N = 6;

  /// Pearson over the bins where EITHER trace is non-zero.
  ///
  /// Correlating over the whole gradient would be dominated by the thousands of
  /// jointly-empty bins, which agree perfectly and mean nothing -- that alone
  /// would drive r towards 1 for every precursor and reproduce the saturation
  /// this probe exists to escape.
  double corr(const float* a, const float* b, std::size_t n)
  {
    std::size_t m = 0;
    double sa = 0.0, sb = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
      if (a[i] == 0.0f && b[i] == 0.0f) { continue; }
      sa += a[i]; sb += b[i]; ++m;
    }
    if (m < 5) { return 0.0; }
    const double ma = sa / m, mb = sb / m;
    double num = 0.0, da = 0.0, db = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
      if (a[i] == 0.0f && b[i] == 0.0f) { continue; }
      const double x = a[i] - ma, y = b[i] - mb;
      num += x * y; da += x * x; db += y * y;
    }
    if (!(da > 0.0) || !(db > 0.0)) { return 0.0; }
    return num / std::sqrt(da * db);
  }
}

int main(int argc, char** argv)
{
  if (argc < 4)
  {
    std::fprintf(stderr, "usage: %s <library.tsv> <run.mzpeak> <out.tsv> [ppm=15]"
                         " [im_window=0.05]\n", argv[0]);
    return 2;
  }
  const double ppm = argc > 4 ? std::atof(argv[4]) : 15.0;
  const double imw = argc > 5 ? std::atof(argv[5]) : 0.05;

  try
  {
    ODIA::Library lib;
    ODIA::DIANNLibraryFile::load(argv[1], lib);
    lib.sortByPrecursorMz();
    const auto& p = lib.precursors();
    const auto& tr = lib.transitions();
    const std::size_t NP = lib.precursorCount();

    auto run = ODIA::openRun(argv[2]);
    const auto& ms1 = run->ms1Spectra();
    const auto& ms2 = run->spectra();
    if (ms1.empty()) { std::fprintf(stderr, "no MS1\n"); return 1; }
    const std::size_t NB = ms1.size();
    std::fprintf(stderr, "%zu precursors, %zu MS1 bins, %zu MS2 entries\n", NP, NB, ms2.size());

    std::vector<float> t1(NP * NB, 0.0f), t2(NP * NB, 0.0f);
    std::fprintf(stderr, "traces: %.1f GiB\n",
                 2.0 * NP * NB * sizeof(float) / 1073741824.0);

    // ---- MS1 trace: monoisotopic intensity per bin -------------------------
    struct T1 { double mz; std::uint32_t slot; float im; };
    std::vector<T1> i1;
    i1.reserve(NP);
    for (std::size_t i = 0; i < NP; ++i)
    {
      const double mz = ODIA::fromFixed(p.mz[i]);
      if (mz > 0.0) { i1.push_back({mz, static_cast<std::uint32_t>(i), p.im[i]}); }
    }
    std::sort(i1.begin(), i1.end(), [](const T1& a, const T1& b) { return a.mz < b.mz; });

    std::vector<ODIA::SpectrumPeaks> blk;
    for (std::size_t b = 0; b < NB; b += 64)
    {
      const std::size_t e = std::min(b + 64, NB);
      run->ms1Peaks(b, e, blk);
      for (std::size_t s = 0; s < blk.size(); ++s)
      {
        const auto& sp = blk[s];
        const bool g = imw > 0.0 && sp.ion_mobility.size() == sp.mz.size();
        for (std::size_t k = 0; k < sp.mz.size(); ++k)
        {
          const double m = sp.mz[k], tol = m * ppm * 1e-6;
          auto it = std::lower_bound(i1.begin(), i1.end(), m - tol,
                                     [](const T1& a, double v) { return a.mz < v; });
          for (; it != i1.end() && it->mz <= m + tol; ++it)
          {
            if (g && !(std::abs(double(sp.ion_mobility[k]) - double(it->im)) <= imw)) { continue; }
            float& c = t1[it->slot * NB + (b + s)];
            c = std::max(c, sp.intensity[k]);
          }
        }
      }
      std::fprintf(stderr, "  MS1 %zu/%zu\r", e, NB);
    }
    std::fprintf(stderr, "\n");

    // ---- MS2 trace: summed top-N fragment intensity, into the nearest bin ---
    struct T2 { double mz; std::uint32_t slot; float im; };
    std::vector<std::vector<T2>> per_window(run->windows().size());
    std::vector<char> built(run->windows().size(), 0);

    for (std::size_t si = 0; si < ms2.size(); ++si)
    {
      const auto& in = ms2[si];
      std::size_t w = run->windows().size();
      for (std::size_t k = 0; k < run->windows().size(); ++k)
      {
        if (run->windows()[k].mz_low == in.window.mz_low &&
            run->windows()[k].mz_high == in.window.mz_high) { w = k; break; }
      }
      if (w == run->windows().size()) { continue; }
      if (!built[w])
      {
        auto& v = per_window[w];
        for (std::size_t i = 0; i < NP; ++i)
        {
          if (!in.window.contains(ODIA::fromFixed(p.mz[i]))) { continue; }
          const std::uint32_t bgn = p.transition_begin[i], n = p.transition_count[i];
          std::vector<std::pair<float, double>> f;
          for (std::uint32_t k = 0; k < n; ++k)
          {
            if (tr.product_mz[bgn + k] == ODIA::MZ_INVALID) { continue; }
            f.push_back({tr.library_intensity[bgn + k], ODIA::fromFixed(tr.product_mz[bgn + k])});
          }
          std::sort(f.begin(), f.end(), [](auto& a, auto& b) { return a.first > b.first; });
          if (f.size() > TOP_N) { f.resize(TOP_N); }
          for (const auto& x : f)
          { v.push_back({x.second, static_cast<std::uint32_t>(i), p.im[i]}); }
        }
        std::sort(v.begin(), v.end(), [](const T2& a, const T2& b) { return a.mz < b.mz; });
        built[w] = 1;
      }
      const auto& v = per_window[w];
      if (v.empty()) { continue; }

      // nearest MS1 bin by retention time
      const auto bit = std::lower_bound(ms1.begin(), ms1.end(), in.retention_time,
                                        [](const ODIA::SpectrumInfo& a, double t)
                                        { return a.retention_time < t; });
      std::size_t bin = static_cast<std::size_t>(bit - ms1.begin());
      if (bin >= NB) { bin = NB - 1; }
      if (bin > 0 && std::abs(ms1[bin - 1].retention_time - in.retention_time) <
                     std::abs(ms1[bin].retention_time - in.retention_time)) { --bin; }

      ODIA::SpectrumPeaks sp;
      run->peaks(si, sp);
      const bool g = imw > 0.0 && sp.ion_mobility.size() == sp.mz.size();
      for (std::size_t k = 0; k < sp.mz.size(); ++k)
      {
        const double m = sp.mz[k], tol = m * ppm * 1e-6;
        auto it = std::lower_bound(v.begin(), v.end(), m - tol,
                                   [](const T2& a, double x) { return a.mz < x; });
        for (; it != v.end() && it->mz <= m + tol; ++it)
        {
          if (g && !(std::abs(double(sp.ion_mobility[k]) - double(it->im)) <= imw)) { continue; }
          t2[it->slot * NB + bin] += sp.intensity[k];
        }
      }
      if (si % 500 == 0) { std::fprintf(stderr, "  MS2 %zu/%zu\r", si, ms2.size()); }
    }
    std::fprintf(stderr, "\n");

    std::FILE* out = std::fopen(argv[3], "w");
    if (!out) { return 1; }
    std::fprintf(out, "Precursor.Id\tDecoy\tms1_ms2_corr\tms1_bins\tms2_bins\n");
    for (std::size_t i = 0; i < NP; ++i)
    {
      const double r = corr(&t1[i * NB], &t2[i * NB], NB);
      std::size_t n1 = 0, n2 = 0;
      for (std::size_t b = 0; b < NB; ++b)
      { if (t1[i * NB + b] > 0.0f) { ++n1; } if (t2[i * NB + b] > 0.0f) { ++n2; } }
      const auto sq = lib.strings().get(p.modified_sequence[i]);
      std::fprintf(out, "%.*s%u\t%u\t%.4f\t%zu\t%zu\n",
                   static_cast<int>(sq.size()), sq.data(),
                   static_cast<unsigned>(p.charge[i]),
                   static_cast<unsigned>(p.decoy[i]), r, n1, n2);
    }
    std::fclose(out);
    std::fprintf(stderr, "wrote %s\n", argv[3]);
    return 0;
  }
  catch (const std::exception& e) { std::fprintf(stderr, "error: %s\n", e.what()); return 1; }
}
