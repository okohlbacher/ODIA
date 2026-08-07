// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

/// The measurement doc/08 rests on, taken on our own data.
///
/// The design claims unconditioned fragment match DEPTH separates 729.7 against
/// 0.39 true identifications per thousand -- 1,800x from one quantity -- where
/// our own `fragment_coverage` is flat at 16-29 per thousand. That claim came
/// from the reference engine's data, not ours, and the whole prefilter plan
/// rests on it reproducing here.
///
/// So measure it before building the classifier. If the gradient is flat on our
/// run, the prefilter is the wrong idea and we have learned that for one pass
/// over the file instead of for a subsystem.
///
/// DEPTH, precisely: over the spectra of a precursor's isolation window, the
/// MAXIMUM number of that precursor's top-6 library fragments (by library
/// intensity) that fall within tolerance in ONE spectrum. Unconditioned -- it
/// asks whether any single moment in the run looks like this peptide, with no
/// peak group, no RT map, and no candidate already chosen.
///
/// Emits one row per precursor so the enrichment table can be computed against
/// any external confident set.

#include <odia/DIANNLibraryFile.h>
#include <odia/Library.h>
#include <odia/SpectrumSource.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace
{
  constexpr std::size_t TOP_N = 6;

  struct Frag { double mz; float intensity; };
}

int main(int argc, char** argv)
{
  if (argc < 4)
  {
    std::fprintf(stderr,
                 "usage: %s <library.tsv> <run.mzpeak> <out.tsv> [ppm=15] [max_spectra=0]\n",
                 argv[0]);
    return 2;
  }
  const double ppm = argc > 4 ? std::atof(argv[4]) : 15.0;
  const std::size_t max_spectra = argc > 5 ? static_cast<std::size_t>(std::atol(argv[5])) : 0;

  try
  {
    ODIA::Library lib;
    ODIA::DIANNLibraryFile::load(argv[1], lib);
    lib.sortByPrecursorMz();
    const auto& p = lib.precursors();
    const auto& t = lib.transitions();
    std::fprintf(stderr, "library: %zu precursors, %zu transitions\n",
                 lib.precursorCount(), lib.transitionCount());

    auto run = ODIA::openRun(argv[2]);
    const auto& info = run->spectra();
    std::fprintf(stderr, "run: %zu spectra\n", info.size());

    // Top-N fragments per precursor, by library intensity, sorted by m/z so a
    // spectrum can be swept once against them.
    std::vector<std::vector<Frag>> want(lib.precursorCount());
    for (std::size_t i = 0; i < lib.precursorCount(); ++i)
    {
      const std::uint32_t b = p.transition_begin[i], n = p.transition_count[i];
      std::vector<Frag> f;
      f.reserve(n);
      for (std::uint32_t k = 0; k < n; ++k)
      {
        const auto mz = t.product_mz[b + k];
        if (mz == ODIA::MZ_INVALID) { continue; }
        f.push_back({ODIA::fromFixed(mz), t.library_intensity[b + k]});
      }
      std::sort(f.begin(), f.end(),
                [](const Frag& a, const Frag& b2) { return a.intensity > b2.intensity; });
      if (f.size() > TOP_N) { f.resize(TOP_N); }
      std::sort(f.begin(), f.end(), [](const Frag& a, const Frag& b2) { return a.mz < b2.mz; });
      want[i] = std::move(f);
    }

    std::vector<std::uint8_t> depth(lib.precursorCount(), 0);
    std::vector<std::uint32_t> qualifying(lib.precursorCount(), 0);
    std::vector<std::uint32_t> total_matches(lib.precursorCount(), 0);
    std::vector<float> best_rt(lib.precursorCount(), -1.0f);

    // Peaks are NOT ascending in m/z. SpectrumSource documents this as a
    // contract: a merged ion-mobility frame is the concatenation of its TIMS
    // scans, so m/z restarts at every mobility step. A consumer that binary-
    // searches the peak array "does not fail loudly -- it returns near-zero
    // matches and a chromatogram of zeros, which looks exactly like an ion that
    // is simply not there." So iterate the PEAKS and search the sorted library
    // side, which is what ChromatogramExtractor does.
    //
    // The searchable side is a flat (mz, precursor-slot) index over the top-N
    // fragments of every precursor in this spectrum's window, built once per
    // window rather than per spectrum.
    struct Target { double mz; std::uint32_t slot; };

    ODIA::SpectrumPeaks sp;
    const std::size_t limit = max_spectra ? std::min(max_spectra, info.size()) : info.size();
    std::size_t done = 0;

    // Cache the per-window target index; runs cycle through a handful of windows.
    std::vector<std::vector<Target>> window_index(run->windows().size());
    std::vector<char> window_built(run->windows().size(), 0);
    std::vector<std::uint8_t> hits;                 // per precursor slot, this spectrum
    std::vector<std::uint32_t> slots_touched;

    for (std::size_t si = 0; si < limit; ++si)
    {
      const auto& in = info[si];
      if (in.ms_level != 2) { continue; }

      // Which window is this? Match on the stated bounds.
      std::size_t w = run->windows().size();
      for (std::size_t k = 0; k < run->windows().size(); ++k)
      {
        if (run->windows()[k].mz_low == in.window.mz_low &&
            run->windows()[k].mz_high == in.window.mz_high) { w = k; break; }
      }
      if (w == run->windows().size()) { continue; }

      if (!window_built[w])
      {
        auto& idx = window_index[w];
        for (std::size_t i = 0; i < lib.precursorCount(); ++i)
        {
          const double mz = ODIA::fromFixed(p.mz[i]);
          if (!in.window.contains(mz)) { continue; }
          for (const auto& fr : want[i])
          { idx.push_back({fr.mz, static_cast<std::uint32_t>(i)}); }
        }
        std::sort(idx.begin(), idx.end(),
                  [](const Target& a, const Target& b) { return a.mz < b.mz; });
        window_built[w] = 1;
      }
      const auto& idx = window_index[w];
      if (idx.empty()) { continue; }

      run->peaks(si, sp);
      if (sp.mz.empty()) { continue; }

      hits.assign(lib.precursorCount(), 0);
      slots_touched.clear();
      for (const double m : sp.mz)
      {
        // Widest tolerance at this m/z; the per-target check below is exact.
        const double tol = m * ppm * 1e-6;
        auto it = std::lower_bound(idx.begin(), idx.end(), m - tol,
                                   [](const Target& a, double v) { return a.mz < v; });
        for (; it != idx.end() && it->mz <= m + tol; ++it)
        {
          if (std::abs(it->mz - m) > it->mz * ppm * 1e-6) { continue; }
          if (hits[it->slot] == 0) { slots_touched.push_back(it->slot); }
          // One fragment may match several peaks in a mobility-merged frame;
          // depth counts DISTINCT fragments, so cap per fragment by tracking
          // the count of matched targets rather than of matched peaks. A
          // fragment already counted cannot raise the depth again.
          if (hits[it->slot] < TOP_N) { ++hits[it->slot]; }
        }
      }

      for (const std::uint32_t slot : slots_touched)
      {
        const std::uint8_t h = std::min<std::uint8_t>(
          hits[slot], static_cast<std::uint8_t>(want[slot].size()));
        if (h == 0) { continue; }
        total_matches[slot] += h;
        if (h > depth[slot])
        {
          depth[slot] = h;
          best_rt[slot] = static_cast<float>(in.retention_time);
          qualifying[slot] = 1;
        }
        else if (h + 1 >= depth[slot]) { ++qualifying[slot]; }
      }
      if (++done % 2000 == 0) { std::fprintf(stderr, "  %zu spectra\r", done); }
    }
    std::fprintf(stderr, "\nswept %zu MS2 spectra\n", done);

    std::FILE* out = std::fopen(argv[3], "w");
    if (!out) { std::fprintf(stderr, "cannot write %s\n", argv[3]); return 1; }
    std::fprintf(out, "Precursor.Id\tDecoy\tCharge\tPrecursor.Mz\tdepth\tqualifying_spectra"
                      "\ttotal_matches\tbest_spectrum_rt\tn_top\n");
    for (std::size_t i = 0; i < lib.precursorCount(); ++i)
    {
      const auto seq = lib.strings().get(p.modified_sequence[i]);
      std::fprintf(out, "%.*s%u\t%u\t%u\t%.5f\t%u\t%u\t%u\t%.2f\t%zu\n",
                   static_cast<int>(seq.size()), seq.data(),
                   static_cast<unsigned>(p.charge[i]),
                   static_cast<unsigned>(p.decoy[i]),
                   static_cast<unsigned>(p.charge[i]),
                   ODIA::fromFixed(p.mz[i]),
                   static_cast<unsigned>(depth[i]), qualifying[i], total_matches[i],
                   static_cast<double>(best_rt[i]), want[i].size());
    }
    std::fclose(out);
    std::fprintf(stderr, "wrote %s\n", argv[3]);
    return 0;
  }
  catch (const std::exception& e)
  {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
}
