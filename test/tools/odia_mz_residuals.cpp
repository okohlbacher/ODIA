// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

/// Is the loss in PICKING, or is extraction already mis-centred in m/z?
///
/// The mass calibration gate FAILS on IH1, so extraction runs on an
/// UNCALIBRATED window: +/-15 ppm centred on zero. If the run's true fragment
/// error has an offset, that window is off-centre, the extracted traces lose
/// intensity, and the failure would present as a picking failure while actually
/// being an extraction one. The picker's rejection census cannot distinguish
/// these, because a trace that is present but attenuated looks like a trace
/// that correlates badly.
///
/// So: probe with a DELIBERATELY WIDE window (+/-50 ppm by default) and record
/// where the fragment matches actually sit. Per precursor, the
/// intensity-weighted mean ppm deviation and its spread, at the brightest
/// co-occurrence found. Joining that against which precursors ODIA emitted a
/// candidate for answers the question directly:
///
///   * if picked and unpicked precursors have the SAME residual distribution,
///     extraction is centred correctly and the loss is genuinely in picking;
///   * if the unpicked ones sit systematically off-centre, or outside +/-15 ppm,
///     extraction never gave the picker a fair trace and the mass calibration
///     is the thing to fix.

#include <odia/DIANNLibraryFile.h>
#include <odia/Library.h>
#include <odia/SpectrumSource.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <limits>
#include <string>
#include <vector>

namespace { constexpr std::size_t TOP_N = 6; }

int main(int argc, char** argv)
{
  if (argc < 4)
  {
    std::fprintf(stderr, "usage: %s <library.tsv> <run.mzpeak> <out.tsv> [ppm=50]"
                         " [im_window=0.05] [anchors.tsv] [rt_window_s=30]\n"
                 "  anchors.tsv: Precursor.Id<TAB>RT_seconds. Restricts the search to a window\n"
                 "  around a KNOWN apex. Without it the 'best co-occurrence' is a maximum over\n"
                 "  every spectrum in the run, which on a mostly-absent library is a chance\n"
                 "  event whose residuals are noise -- the extreme-value trap that makes every\n"
                 "  presence statistic on this data saturate.\n", argv[0]);
    return 2;
  }
  const double ppm = argc > 4 ? std::atof(argv[4]) : 50.0;
  const double imw = argc > 5 ? std::atof(argv[5]) : 0.05;
  const char* anchor_path = argc > 6 ? argv[6] : nullptr;
  const double rt_win = argc > 7 ? std::atof(argv[7]) : 30.0;

  try
  {
    ODIA::Library lib;
    ODIA::DIANNLibraryFile::load(argv[1], lib);
    lib.sortByPrecursorMz();
    const auto& p = lib.precursors();
    const auto& t = lib.transitions();
    const std::size_t NP = lib.precursorCount();

    // Optional known-apex anchors, keyed by the same Precursor.Id the output uses.
    std::vector<double> anchor_rt(lib.precursorCount(),
                                  std::numeric_limits<double>::quiet_NaN());
    std::size_t n_anchor = 0;
    if (anchor_path)
    {
      std::map<std::string, double> want;
      std::ifstream in(anchor_path);
      std::string line;
      std::getline(in, line);                     // header
      while (std::getline(in, line))
      {
        const auto tab = line.find('\t');
        if (tab == std::string::npos) { continue; }
        want[line.substr(0, tab)] = std::atof(line.c_str() + tab + 1);
      }
      for (std::size_t i = 0; i < lib.precursorCount(); ++i)
      {
        const auto sq = lib.strings().get(p.modified_sequence[i]);
        const std::string id = std::string(sq) + std::to_string(unsigned(p.charge[i]));
        const auto it = want.find(id);
        if (it != want.end()) { anchor_rt[i] = it->second; ++n_anchor; }
      }
      std::fprintf(stderr, "anchors: %zu of %zu matched, +/-%.0f s\n",
                   n_anchor, want.size(), rt_win);
    }

    auto run = ODIA::openRun(argv[2]);
    const auto& info = run->spectra();
    std::fprintf(stderr, "%zu precursors, %zu MS2 entries, +/-%.0f ppm probe\n",
                 NP, info.size(), ppm);

    // Best co-occurrence per precursor: the spectrum where the most of its
    // top-6 fragments appear at once, tie-broken by summed intensity. The
    // residuals are taken THERE -- at a random spectrum they would be noise.
    struct Best { std::uint8_t n = 0; float total = 0.0f;
                  double wsum = 0.0, wtot = 0.0;      // intensity-weighted ppm
                  double lo = 0.0, hi = 0.0;
                  float rt = 0.0f;
                  // Per-FRAGMENT residuals at the winning cell. A per-precursor
                  // mean cannot show a surface in (RT, m/z), and it hides the
                  // contamination: at +/-50 ppm the within-precursor spread is
                  // ~97 ppm, i.e. most matches are not the fragment.
                  std::vector<float> f_mz, f_ppm, f_int; };
    std::vector<Best> best(NP);

    struct Target { double mz; std::uint32_t slot; float im; };
    std::vector<std::vector<Target>> widx(run->windows().size());
    std::vector<char> built(run->windows().size(), 0);

    std::vector<std::uint8_t> cnt;
    std::vector<float> tot;
    std::vector<double> ws, wt, mn, mx;
    std::vector<std::uint32_t> touched, touched_prev;
    std::vector<std::vector<float>> fmz(NP), fppm(NP), fint(NP);
    ODIA::SpectrumPeaks sp;

    for (std::size_t si = 0; si < info.size(); ++si)
    {
      const auto& in = info[si];
      std::size_t w = run->windows().size();
      for (std::size_t k = 0; k < run->windows().size(); ++k)
      {
        if (run->windows()[k].mz_low == in.window.mz_low &&
            run->windows()[k].mz_high == in.window.mz_high) { w = k; break; }
      }
      if (w == run->windows().size()) { continue; }
      if (!built[w])
      {
        auto& v = widx[w];
        for (std::size_t i = 0; i < NP; ++i)
        {
          if (!in.window.contains(ODIA::fromFixed(p.mz[i]))) { continue; }
          const std::uint32_t b = p.transition_begin[i], n = p.transition_count[i];
          std::vector<std::pair<float, double>> f;
          for (std::uint32_t k = 0; k < n; ++k)
          {
            if (t.product_mz[b + k] == ODIA::MZ_INVALID) { continue; }
            f.push_back({t.library_intensity[b + k], ODIA::fromFixed(t.product_mz[b + k])});
          }
          std::sort(f.begin(), f.end(), [](auto& a, auto& c) { return a.first > c.first; });
          if (f.size() > TOP_N) { f.resize(TOP_N); }
          for (const auto& x : f)
          { v.push_back({x.second, static_cast<std::uint32_t>(i), p.im[i]}); }
        }
        std::sort(v.begin(), v.end(),
                  [](const Target& a, const Target& b2) { return a.mz < b2.mz; });
        built[w] = 1;
      }
      const auto& v = widx[w];
      if (v.empty()) { continue; }

      run->peaks(si, sp);
      if (sp.mz.empty()) { continue; }
      const bool g = imw > 0.0 && sp.ion_mobility.size() == sp.mz.size();

      for (const std::uint32_t s : touched_prev) { fmz[s].clear(); fppm[s].clear(); fint[s].clear(); }
      touched_prev.clear();
      cnt.assign(NP, 0); tot.assign(NP, 0.0f);
      ws.assign(NP, 0.0); wt.assign(NP, 0.0);
      mn.assign(NP, 1e9); mx.assign(NP, -1e9);
      touched.clear();
      for (std::size_t k = 0; k < sp.mz.size(); ++k)
      {
        const double m = sp.mz[k], tol = m * ppm * 1e-6;
        auto it = std::lower_bound(v.begin(), v.end(), m - tol,
                                   [](const Target& a, double x) { return a.mz < x; });
        for (; it != v.end() && it->mz <= m + tol; ++it)
        {
          const double d = (m - it->mz) / it->mz * 1e6;     // SIGNED ppm
          if (std::abs(d) > ppm) { continue; }
          if (g && !(std::abs(double(sp.ion_mobility[k]) - double(it->im)) <= imw)) { continue; }
          const std::uint32_t s = it->slot;
          // Only where this precursor is known to elute, when that is known.
          if (anchor_path)
          {
            if (!std::isfinite(anchor_rt[s])) { continue; }
            if (std::abs(in.retention_time - anchor_rt[s]) > rt_win) { continue; }
          }
          if (cnt[s] == 0) { touched.push_back(s); touched_prev.push_back(s); }
          fmz[s].push_back(static_cast<float>(it->mz));
          fppm[s].push_back(static_cast<float>(d));
          fint[s].push_back(sp.intensity[k]);
          if (cnt[s] < TOP_N) { ++cnt[s]; }
          tot[s] += sp.intensity[k];
          ws[s] += d * sp.intensity[k]; wt[s] += sp.intensity[k];
          mn[s] = std::min(mn[s], d); mx[s] = std::max(mx[s], d);
        }
      }
      for (const std::uint32_t s : touched)
      {
        if (cnt[s] > best[s].n || (cnt[s] == best[s].n && tot[s] > best[s].total))
        {
          best[s] = {cnt[s], tot[s], ws[s], wt[s], mn[s], mx[s],
                     static_cast<float>(in.retention_time),
                     fmz[s], fppm[s], fint[s]};
        }
      }
      if (si % 1000 == 0) { std::fprintf(stderr, "  %zu/%zu\r", si, info.size()); }
    }
    std::fprintf(stderr, "\n");

    std::FILE* out = std::fopen(argv[3], "w");
    if (!out) { return 1; }
    std::fprintf(out, "Precursor.Id\tDecoy\tn_matched\trt\tfrag_mz\tppm\tintensity\n");
    for (std::size_t i = 0; i < NP; ++i)
    {
      const auto& b = best[i];
      if (b.f_ppm.empty()) { continue; }
      const auto sq = lib.strings().get(p.modified_sequence[i]);
      for (std::size_t k = 0; k < b.f_ppm.size(); ++k)
      {
        std::fprintf(out, "%.*s%u\t%u\t%u\t%.2f\t%.4f\t%.3f\t%.1f\n",
                     static_cast<int>(sq.size()), sq.data(),
                     static_cast<unsigned>(p.charge[i]),
                     static_cast<unsigned>(p.decoy[i]),
                     static_cast<unsigned>(b.n),
                     static_cast<double>(b.rt), static_cast<double>(b.f_mz[k]),
                     static_cast<double>(b.f_ppm[k]), static_cast<double>(b.f_int[k]));
      }
    }
    std::fclose(out);
    std::fprintf(stderr, "wrote %s\n", argv[3]);
    return 0;
  }
  catch (const std::exception& e) { std::fprintf(stderr, "error: %s\n", e.what()); return 1; }
}
