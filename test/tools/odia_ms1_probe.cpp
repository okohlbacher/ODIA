// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

/// Does MS1 evidence discriminate a present peptide from an absent one?
///
/// MS1 is the last untested hypothesis for the 0-versus-738 collapse, and the
/// argument for it is ORTHOGONALITY: every one of ODIA's 15 sub-scores reads MS2
/// fragment traces, so a co-eluting interferent corrupts all fifteen at once.
/// The feature-COUNT argument that originally motivated it was measured to be an
/// artefact of independent feature draws and does not survive.
///
/// So the premise gets tested before the subsystem gets built. Two rounds were
/// just spent discovering that doc/08's prefilter premise does not transfer to
/// diaPASEF; the lesson is to spend one pass over the file first.
///
/// Per precursor, over MS1 spectra, optionally restricted to the precursor's
/// mobility slice:
///   * `ms1_iso`  -- the most of {M, M+1, M+2} seen together in one spectrum
///   * `ms1_max`  -- the largest monoisotopic intensity over the run
///   * `ms1_rt`   -- when that maximum occurred (the RT seed, for free)
///
/// If true and absent precursors are indistinguishable on these, MS1 is not the
/// lever either and the build is not worth starting.

#include <odia/DIANNLibraryFile.h>
#include <odia/Library.h>
#include <odia/SpectrumSource.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

int main(int argc, char** argv)
{
  if (argc < 4)
  {
    std::fprintf(stderr, "usage: %s <library.tsv> <run.mzpeak> <out.tsv> [ppm=15]"
                         " [im_window=0]\n", argv[0]);
    return 2;
  }
  const double ppm = argc > 4 ? std::atof(argv[4]) : 15.0;
  const double im_window = argc > 5 ? std::atof(argv[5]) : 0.0;

  try
  {
    ODIA::Library lib;
    ODIA::DIANNLibraryFile::load(argv[1], lib);
    lib.sortByPrecursorMz();
    const auto& p = lib.precursors();

    auto run = ODIA::openRun(argv[2]);
    const auto& ms1 = run->ms1Spectra();
    std::fprintf(stderr, "library %zu precursors; run has %zu MS1 spectra\n",
                 lib.precursorCount(), ms1.size());
    if (ms1.empty()) { std::fprintf(stderr, "no MS1 in this run\n"); return 1; }

    // Searchable side: the monoisotopic m/z and the first two isotopes, sorted.
    // Isotope spacing is 1.00335 Da / charge -- the C13 mass defect, not 1.0.
    struct Target { double mz; std::uint32_t slot; std::uint8_t iso; float im; };
    std::vector<Target> idx;
    idx.reserve(lib.precursorCount() * 3);
    for (std::size_t i = 0; i < lib.precursorCount(); ++i)
    {
      const double mz = ODIA::fromFixed(p.mz[i]);
      if (mz <= 0.0) { continue; }
      const double z = std::max<double>(1, p.charge[i]);
      for (std::uint8_t k = 0; k < 3; ++k)
      {
        idx.push_back({mz + 1.00335 * k / z, static_cast<std::uint32_t>(i), k, p.im[i]});
      }
    }
    std::sort(idx.begin(), idx.end(),
              [](const Target& a, const Target& b) { return a.mz < b.mz; });

    std::vector<std::uint8_t> iso(lib.precursorCount(), 0);
    std::vector<float> mx(lib.precursorCount(), 0.0f);
    std::vector<float> rt(lib.precursorCount(), -1.0f);

    std::vector<ODIA::SpectrumPeaks> block;
    std::vector<std::uint8_t> seen;
    std::vector<float> mono;
    std::vector<std::uint32_t> touched;

    const std::size_t STEP = 64;
    for (std::size_t b = 0; b < ms1.size(); b += STEP)
    {
      const std::size_t e = std::min(b + STEP, ms1.size());
      run->ms1Peaks(b, e, block);
      for (std::size_t s = 0; s < block.size(); ++s)
      {
        const auto& sp = block[s];
        if (sp.mz.empty()) { continue; }
        const bool gated = im_window > 0.0 && sp.ion_mobility.size() == sp.mz.size();
        seen.assign(lib.precursorCount(), 0);
        mono.assign(lib.precursorCount(), 0.0f);
        touched.clear();
        for (std::size_t k = 0; k < sp.mz.size(); ++k)
        {
          const double m = sp.mz[k];
          const double tol = m * ppm * 1e-6;
          auto it = std::lower_bound(idx.begin(), idx.end(), m - tol,
                                     [](const Target& a, double v) { return a.mz < v; });
          for (; it != idx.end() && it->mz <= m + tol; ++it)
          {
            if (std::abs(it->mz - m) > it->mz * ppm * 1e-6) { continue; }
            if (gated && !(std::abs(static_cast<double>(sp.ion_mobility[k]) -
                                    static_cast<double>(it->im)) <= im_window))
            { continue; }
            if (seen[it->slot] == 0) { touched.push_back(it->slot); }
            seen[it->slot] |= static_cast<std::uint8_t>(1u << it->iso);
            if (it->iso == 0) { mono[it->slot] = std::max(mono[it->slot], sp.intensity[k]); }
          }
        }
        for (const std::uint32_t t : touched)
        {
          const std::uint8_t n = static_cast<std::uint8_t>(
            (seen[t] & 1) + ((seen[t] >> 1) & 1) + ((seen[t] >> 2) & 1));
          if (n > iso[t]) { iso[t] = n; }
          if (mono[t] > mx[t])
          {
            mx[t] = mono[t];
            rt[t] = static_cast<float>(ms1[b + s].retention_time);
          }
        }
      }
      std::fprintf(stderr, "  %zu/%zu MS1\r", e, ms1.size());
    }
    std::fprintf(stderr, "\n");

    std::FILE* out = std::fopen(argv[3], "w");
    if (!out) { return 1; }
    std::fprintf(out, "Precursor.Id\tDecoy\tms1_iso\tms1_max\tms1_rt\n");
    for (std::size_t i = 0; i < lib.precursorCount(); ++i)
    {
      const auto sq = lib.strings().get(p.modified_sequence[i]);
      std::fprintf(out, "%.*s%u\t%u\t%u\t%.1f\t%.2f\n",
                   static_cast<int>(sq.size()), sq.data(),
                   static_cast<unsigned>(p.charge[i]),
                   static_cast<unsigned>(p.decoy[i]),
                   static_cast<unsigned>(iso[i]),
                   static_cast<double>(mx[i]), static_cast<double>(rt[i]));
    }
    std::fclose(out);
    std::fprintf(stderr, "wrote %s\n", argv[3]);
    return 0;
  }
  catch (const std::exception& e) { std::fprintf(stderr, "error: %s\n", e.what()); return 1; }
}
