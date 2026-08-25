// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

/// The scoring interval and the quantification interval must stay separate.
///
/// They were one interval, and that was the defect: every sub-score integrated
/// the walked boundaries, the walked boundaries covered 129 of 130 cycles, and
/// so a correct candidate and a wrong one 20 cycles away received the same
/// feature vector -- median paired difference exactly 0. Splitting them is what
/// made the sub-scores position-sensitive again.
///
/// A split like that is easy to undo by accident. Any future change that reads
/// `cand.left`/`cand.right` from inside a sub-score, or reports the scoring
/// window as the group's retention-time range, silently re-couples them and no
/// existing test would notice -- the numbers would still be plausible. Hence
/// two invariances, which is the pair a reviewer asked for:
///
///   1. move the QUANTIFICATION bounds, hold the apex: every sub-score must be
///      bit-identical, and the reported RT range must change.
///   2. move the SCORING window, hold the bounds: the reported RT range must be
///      bit-identical, and the sub-scores must change.
///
/// Each half fails for a different re-coupling, which is why both are needed.

#include <odia/Library.h>
#include <odia/PeakGroupScorer.h>

#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

namespace
{
  int failures = 0;

  void check(bool ok, const std::string& what)
  {
    std::printf("  [%s] %s\n", ok ? "ok" : "FAIL", what.c_str());
    if (!ok) { ++failures; }
  }

  constexpr std::uint32_t CYCLES = 130;
  constexpr std::uint32_t TRANSITIONS = 6;
  constexpr std::uint32_t APEX = 65;

  ODIA::Library onePrecursorLibrary()
  {
    ODIA::Library lib;
    auto& p = lib.precursors();
    p.mz.push_back(ODIA::toFixed(500.0));
    p.irt.push_back(0.0f);
    p.im.push_back(std::numeric_limits<float>::quiet_NaN());
    p.ccs.push_back(std::numeric_limits<float>::quiet_NaN());
    p.charge.push_back(2);
    p.decoy.push_back(0);
    p.modified_sequence.push_back(0);
    p.protein_group.push_back(0);
    p.transition_begin.push_back(0);
    p.transition_count.push_back(TRANSITIONS);
    auto& t = lib.transitions();
    for (std::uint32_t k = 0; k < TRANSITIONS; ++k)
    {
      t.product_mz.push_back(ODIA::toFixed(300.0 + 50.0 * k));
      // Unequal, so LIBRARY_CORR is a real number rather than a degenerate one:
      // a correlation against a constant vector is undefined, and a test that
      // compares undefined against undefined passes for the wrong reason.
      t.library_intensity.push_back(static_cast<float>(1.0 + 0.4 * k));
      t.type.push_back(ODIA::FragmentType::Y);
      t.ordinal.push_back(static_cast<std::uint8_t>(4 + k));
      t.charge.push_back(1);
      t.loss.push_back(ODIA::LossType::None);
    }
    return lib;
  }

  /// One co-eluting peak of the measured width on a baseline, with the
  /// fragments carrying the library's intensity ratios.
  ODIA::Chromatograms scriptedRun()
  {
    ODIA::Chromatograms c;
    c.axes.assign(1, std::vector<float>(CYCLES, 0.0f));
    for (std::uint32_t j = 0; j < CYCLES; ++j)
    { c.axes[0][j] = static_cast<float>(j) * 1.385f; }
    c.precursor_axis.assign(1, 0);
    c.precursor_axis_begin.assign(1, 0);
    c.precursor_cycles.assign(1, CYCLES);
    c.precursor_transition_begin.assign(1, 0);
    c.begin.assign(TRANSITIONS, 0);
    c.count.assign(TRANSITIONS, CYCLES);
    c.intensity.assign(std::size_t(TRANSITIONS) * CYCLES, 0.0f);
    for (std::uint32_t k = 0; k < TRANSITIONS; ++k)
    {
      c.begin[k] = std::uint64_t(k) * CYCLES;
      for (std::uint32_t j = 0; j < CYCLES; ++j)
      {
        const double d = double(j) - double(APEX);
        const double v = (1.0 + 0.4 * k) * 100.0 * std::exp(-(d * d) / (2.0 * 1.07 * 1.07))
                       + 6.0;
        c.intensity[c.begin[k] + j] = static_cast<float>(v);
      }
    }
    return c;
  }

  const ODIA::PeakGroupScorer::PeakGroup* best(const ODIA::PeakGroupScorer::Result& r)
  {
    const ODIA::PeakGroupScorer::PeakGroup* b = nullptr;
    for (const auto& g : r.groups)
    { if (!b || g.dscore > b->dscore) { b = &g; } }
    return b;
  }

  ODIA::PeakGroupScorer::Options base()
  {
    ODIA::PeakGroupScorer::Options o;
    o.coelution_picking = true;
    o.max_candidates = 5;
    o.min_fragments_at_apex = 1;
    return o;
  }
} // namespace

int main()
{
  const auto lib = onePrecursorLibrary();
  const auto run = scriptedRun();

  // 1. The quantification bounds move; the scores must not.
  {
    auto a = base(); a.peak_min_cycles = 7;
    auto b = base(); b.peak_min_cycles = 41;
    const auto ra = ODIA::PeakGroupScorer::score(lib, run, a);
    const auto rb = ODIA::PeakGroupScorer::score(lib, run, b);
    const auto* ga = best(ra); const auto* gb = best(rb);
    if (!ga || !gb)
    {
      std::printf("      groups: %zu and %zu\n", ra.groups.size(), rb.groups.size());
      check(false, "both configurations produced a group"); return 1;
    }

    const double wa = double(ga->right_rt) - double(ga->left_rt);
    const double wb = double(gb->right_rt) - double(gb->left_rt);
    std::printf("quant width %.1f s -> %.1f s with min_cycles 7 -> 41\n", wa, wb);
    check(wb > wa + 1.0,
          "the reported RT range DID widen, so the comparison below is not vacuous");

    bool same = ga->sub_scores.size() == gb->sub_scores.size();
    std::size_t differing = 0;
    for (std::size_t i = 0; same && i < ga->sub_scores.size(); ++i)
    {
      const double x = ga->sub_scores[i], y = gb->sub_scores[i];
      // NaN == NaN is false and several sub-scores are legitimately NaN here.
      if (std::isnan(x) && std::isnan(y)) { continue; }
      if (!(x == y)) { ++differing; }
    }
    check(same && differing == 0,
          "and EVERY sub-score is unchanged: widening the quantification "
          "interval must not reach the identification features");
    if (differing) { std::printf("      %zu sub-scores moved\n", differing); }
  }

  // 2. The scoring window moves; the reported range must not.
  {
    auto a = base(); a.score_half_cycles = 2;
    auto b = base(); b.score_half_cycles = 8;
    const auto ra = ODIA::PeakGroupScorer::score(lib, run, a);
    const auto rb = ODIA::PeakGroupScorer::score(lib, run, b);
    const auto* ga = best(ra); const auto* gb = best(rb);
    if (!ga || !gb) { check(false, "both configurations produced a group"); return 1; }

    check(ga->left_rt == gb->left_rt && ga->right_rt == gb->right_rt,
          "the reported RT range is identical: the scoring window must not "
          "reach quantification either -- the invariance runs both ways");

    std::size_t differing = 0;
    for (std::size_t i = 0; i < ga->sub_scores.size() && i < gb->sub_scores.size(); ++i)
    {
      const double x = ga->sub_scores[i], y = gb->sub_scores[i];
      if (std::isnan(x) && std::isnan(y)) { continue; }
      if (!(x == y)) { ++differing; }
    }
    std::printf("scoring window 5 -> 17 cycles moved %zu sub-scores\n", differing);
    check(differing > 0,
          "and the sub-scores DID move, so the first test's zero means the "
          "interval is separate rather than that nothing depends on width");
  }

  if (failures == 0) { std::printf("the two intervals are independent\n"); }
  return failures == 0 ? 0 : 1;
}
