// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// The invariants that make selecting precursors BEFORE a target-decoy FDR
// legitimate. These are not properties of a good filter -- they are the
// conditions under which a filter is allowed to exist at all, and every one of
// them fails silently: the run completes, the q-values look normal, and they
// are wrong.

#include <odia/PrecursorPrefilter.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdio>
#include <random>
#include <vector>

namespace
{
  int failures = 0;

  void check(bool ok, const char* what)
  {
    std::printf("%s  %s\n", ok ? "ok  " : "FAIL", what);
    if (!ok) { ++failures; }
  }

  /// A library of `n` targets and `n` decoys. Evidence is supplied directly so
  /// the selection rules can be tested without a run.
  ODIA::Library makeLibrary(std::size_t n)
  {
    ODIA::Library lib;
    auto& p = lib.precursors();
    for (std::size_t i = 0; i < 2 * n; ++i)
    {
      p.mz.push_back(ODIA::toFixed(500.0 + double(i) * 0.001));
      p.irt.push_back(float(i % 100));
      p.im.push_back(0.9f);
      p.ccs.push_back(0.0f);
      p.charge.push_back(2);
      p.decoy.push_back(i >= n ? 1 : 0);      // second half are decoys
      p.modified_sequence.push_back(0);
      p.protein_group.push_back(0);
      p.transition_begin.push_back(0);
      p.transition_count.push_back(0);
    }
    return lib;
  }
}

int main()
{
  std::mt19937 rng(20260811);
  constexpr std::size_t N = 20000;

  // Targets carry MORE evidence than decoys, which is the realistic case and
  // the one that makes a shared threshold dangerous.
  ODIA::Library lib = makeLibrary(N);
  std::vector<ODIA::PrecursorPrefilter::Evidence> ev(2 * N);
  std::binomial_distribution<int> tgt_depth(6, 0.55);
  std::binomial_distribution<int> dec_depth(6, 0.35);
  for (std::size_t i = 0; i < 2 * N; ++i)
  {
    const bool decoy = i >= N;
    ev[i].depth = static_cast<std::uint8_t>(decoy ? dec_depth(rng) : tgt_depth(rng));
    ev[i].total_matches = ev[i].depth * 10u + (i % 7);
  }

  ODIA::PrecursorPrefilter::Options po;
  po.top_n = 6;
  po.keep_fraction = 0.25;
  po.min_keep_per_class = 10;

  ODIA::PrecursorPrefilter::Stats ps;
  const auto keep = ODIA::PrecursorPrefilter::select(lib, ev, po, ps);

  // --- rule 1: label symmetry, by COUNT ------------------------------------
  {
    std::size_t kt = 0, kd = 0;
    for (std::size_t i = 0; i < 2 * N; ++i)
    { if (keep[i]) { (i >= N ? kd : kt)++; } }

    check(kt == kd, "exactly as many decoys as targets are retained");
    check(kt == std::size_t(0.25 * double(N)), "and the retained count is the requested fraction");
    check(ps.targets_kept == kt && ps.decoys_kept == kd, "the reported counts match reality");

    // The failure this rule exists to prevent. A SHARED depth threshold, on
    // the same evidence, keeps far fewer decoys than targets -- and those
    // survivors are the strongest decoys, so the null they form is not a fair
    // sample of the null the scorer will face. Demonstrate the asymmetry
    // directly, because it is the thing that silently breaks the FDR.
    std::size_t st = 0, sd = 0;
    for (std::size_t i = 0; i < 2 * N; ++i)
    { if (ev[i].depth >= 4) { (i >= N ? sd : st)++; } }
    check(st > sd * 2,
          "a shared threshold would have kept >2x more targets than decoys "
          "(which is why the cut is by count)");
  }

  // --- the retained set really is the best-ranked ---------------------------
  {
    // nth_element only PARTITIONS, so this checks the partition is the right
    // one: nothing discarded may outrank anything retained.
    std::uint8_t min_kept = 255, max_dropped = 0;
    for (std::size_t i = 0; i < N; ++i)          // targets only
    {
      if (keep[i]) { min_kept = std::min(min_kept, ev[i].depth); }
      else         { max_dropped = std::max(max_dropped, ev[i].depth); }
    }
    check(min_kept >= max_dropped,
          "no discarded target outranks a retained one");
    check(ps.depth_threshold == min_kept,
          "the reported cut depth is the MINIMUM retained, not an arbitrary "
          "element of the partition");
  }

  // --- rule 3: the histogram is reported, and per class ---------------------
  {
    std::size_t th = 0, dh = 0;
    for (const auto v : ps.depth_hist_target) { th += v; }
    for (const auto v : ps.depth_hist_decoy)  { dh += v; }
    check(th == N && dh == N, "every precursor appears in its class histogram");
    check(ps.depth_hist_target.size() == po.top_n + 1,
          "the histogram spans depth 0..top_n inclusive");
  }

  // --- keep_fraction >= 1 must discard NOTHING ------------------------------
  {
    // "measure" mode has to leave pass 2 bit-for-bit identical to a run with no
    // filter at all, or the measurement it exists to produce cannot be
    // compared against one.
    ODIA::PrecursorPrefilter::Options all = po;
    all.keep_fraction = 1.0;
    ODIA::PrecursorPrefilter::Stats s2;
    const auto k2 = ODIA::PrecursorPrefilter::select(lib, ev, all, s2);
    check(std::count(k2.begin(), k2.end(), char(0)) == 0, "measurement mode keeps everything");
    check(s2.targets_kept == N && s2.decoys_kept == N, "and says so");
    check(!s2.note.empty(), "and marks itself as a measurement");
  }

  // --- the floor cannot be undercut ----------------------------------------
  {
    // A filter that fails on an unusual run must degrade to "kept too much",
    // never to an empty search.
    ODIA::PrecursorPrefilter::Options tiny = po;
    tiny.keep_fraction = 0.0;
    tiny.min_keep_per_class = 500;
    ODIA::PrecursorPrefilter::Stats s3;
    ODIA::PrecursorPrefilter::select(lib, ev, tiny, s3);
    check(s3.targets_kept == 500 && s3.decoys_kept == 500,
          "keep_fraction 0 still retains min_keep_per_class from each class");
  }

  // --- an undiscriminating statistic must not look successful --------------
  {
    // The 2026-08-08 result: targets and decoys drawn from the SAME
    // distribution. The filter cannot detect this itself -- it will happily
    // retain a quarter of each -- so what has to hold is that the reported
    // histograms are identical in shape, which is the signal a reader uses to
    // refuse the filter. This test exists so that "the histograms look the
    // same" stays a detectable condition rather than an eyeballed one.
    std::vector<ODIA::PrecursorPrefilter::Evidence> flat(2 * N);
    std::binomial_distribution<int> same(6, 0.5);
    for (auto& e : flat) { e.depth = static_cast<std::uint8_t>(same(rng)); }
    ODIA::PrecursorPrefilter::Stats s4;
    ODIA::PrecursorPrefilter::select(lib, flat, po, s4);

    const double td = double(s4.depth_hist_target.back());
    const double dd = double(s4.depth_hist_decoy.back());
    const double ratio = dd > 0 ? td / dd : 0.0;
    check(ratio > 0.8 && ratio < 1.25,
          "an undiscriminating statistic yields a full-depth ratio near 1.0");
    check(s4.targets_kept == s4.decoys_kept,
          "and label symmetry still holds, so the FDR stays fair even then");
  }

  // --- contiguity is reported per class and is what the seed keys on --------
  {
    // The property that makes contiguity worth having: it is ORTHOGONAL to
    // depth. Give targets and decoys the SAME depth distribution -- so depth
    // cannot separate them at all, reproducing the 2026-08-08 result -- and
    // let only the targets repeat across consecutive cycles. A statistic that
    // is merely a proxy for depth would show nothing here.
    std::vector<ODIA::PrecursorPrefilter::Evidence> e2(2 * N);
    std::binomial_distribution<int> same(6, 0.5);
    std::binomial_distribution<int> tgt_run(8, 0.45);
    for (std::size_t i = 0; i < 2 * N; ++i)
    {
      e2[i].depth = static_cast<std::uint8_t>(same(rng));
      e2[i].contiguity = static_cast<std::uint16_t>(i >= N ? (rng() % 2) : tgt_run(rng));
    }

    ODIA::PrecursorPrefilter::Stats s5;
    ODIA::PrecursorPrefilter::select(lib, e2, po, s5);

    std::size_t th = 0, dh = 0;
    for (const auto v : s5.contig_hist_target) { th += v; }
    for (const auto v : s5.contig_hist_decoy)  { dh += v; }
    check(th == N && dh == N, "every precursor appears in its contiguity histogram");

    // Depth says nothing here, by construction.
    const double dratio = double(s5.depth_hist_target.back()) /
                          std::max<double>(1.0, double(s5.depth_hist_decoy.back()));
    check(dratio > 0.75 && dratio < 1.35,
          "depth does NOT separate this population (as on 2026-08-08)");

    // Contiguity does. Count at the seed's default threshold of 3 cycles.
    std::size_t t3 = 0, d3 = 0;
    for (std::size_t c = 3; c < s5.contig_hist_target.size(); ++c)
    { t3 += s5.contig_hist_target[c]; d3 += s5.contig_hist_decoy[c]; }
    check(d3 == 0 || double(t3) / double(d3) > 5.0,
          "contiguity DOES separate it -- the statistic is orthogonal to depth, "
          "not a proxy for it");
  }

  // --- a run must be CONSECUTIVE, which is the whole claim ------------------
  {
    // Contiguity counts consecutive cycles of one precursor's OWN isolation
    // window. Spectra of a run interleave its windows, so consecutive spectrum
    // INDICES are a different and wrong thing -- successive spectra of one
    // window are tens of indices apart. This models the invariant directly.
    const std::uint32_t NO_CYCLE = 0xFFFFFFFFu;
    auto longest = [&](const std::vector<std::uint32_t>& cycles) {
      std::uint32_t last = NO_CYCLE; std::uint16_t cur = 0, best = 0;
      for (const std::uint32_t c : cycles)
      {
        if (last != NO_CYCLE && last + 1 == c) { ++cur; } else { cur = 1; }
        last = c; best = std::max(best, cur);
      }
      return best;
    };
    check(longest({4, 5, 6}) == 3, "three consecutive cycles is a run of 3");
    check(longest({4, 6, 8}) == 1, "every other cycle is NOT a run -- it is three runs of 1");
    check(longest({1, 2, 9, 10, 11}) == 3, "the LONGEST run is reported, not the first");
    check(longest({}) == 0, "no qualifying cycles is a run of 0");
  }

  // --- a missing library 1/K0 must UNGATE, never reject ---------------------
  {
    // The gate is `!(|observed - library| <= window)`. Every comparison against
    // NaN is false, so the negation makes a NaN library value reject EVERY
    // peak -- a precursor with no predicted mobility scores zero on a run it
    // may well be present in, and does so silently.
    //
    // This shipped. On IH1 the seed sweep returned 0 target and 0 decoy anchors
    // at every contiguity threshold, against 13,360 / 8,264 on Astral, because
    // our generated library predicts CCS and leaves 1/K0 unset. Astral has no
    // ion mobility so the gate never ran, which made a missing NaN case look
    // like an instrument-specific result.
    const double window = 0.05;
    const float observed = 1.30f;
    const float absent = std::numeric_limits<float>::quiet_NaN();

    auto rejected = [&](float lib) {
      return !(std::abs(double(observed) - double(lib)) <= window);
    };
    check(rejected(absent), "the bare gate DOES reject a NaN library value");

    auto guarded = [&](float lib) {
      if (std::isnan(lib)) { return false; }          // ungated
      return !(std::abs(double(observed) - double(lib)) <= window);
    };
    check(!guarded(absent), "the guarded gate lets it through instead");
    check(!guarded(1.32f), "and still accepts a peak inside the window");
    check(guarded(0.90f), "and still rejects one outside it");
    check(guarded(0.0f),
          "0.0 is a VALUE, not 'absent' -- it must still gate, which is why the "
          "library writes an empty field rather than a zero for what it lacks");
  }

  std::printf("%s\n", failures == 0 ? "ALL PASSED" : "FAILURES");
  return failures == 0 ? 0 : 1;
}
