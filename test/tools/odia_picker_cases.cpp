// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

/// Does the peak picker find the peptide, or the loudest thing in the window?
///
/// This exists because swapping the entire detection algorithm -- amplitude
/// local-maxima to DIA-NN's pairwise-correlation search -- left all 59 tests
/// passing. Nothing in the suite was sensitive to WHICH peaks were found, only
/// to whether the plumbing around them held together, so the change that moved
/// rank-1 accuracy from 40.7% to 75.7% on real data was invisible to CI.
///
/// The scenario is the one that made the old detector fail, reduced to its
/// smallest form: a genuine co-eluting peak of moderate height, and a single
/// loud fragment elsewhere in the window. Summing the traces makes the spike
/// the tallest feature, so an amplitude detector picks it; the fragments of a
/// real peptide come from one eluting molecule and correlate, which is what
/// the correlation detector selects on.

#include <odia/Library.h>
#include <odia/PeakGroupScorer.h>

#include <cmath>
#include <limits>
#include <cstdio>
#include <cstdlib>
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

  constexpr std::uint32_t CYCLES = 120;
  constexpr std::uint32_t TRANSITIONS = 6;
  constexpr std::uint32_t REAL_APEX = 40;   ///< where the peptide elutes
  constexpr std::uint32_t SPIKE_APEX = 80;  ///< where one loud fragment sits

  /// A Gaussian-ish bump, so the traces have a shape to correlate rather than
  /// a single non-zero sample.
  double bump(std::uint32_t j, std::uint32_t centre, double height, double width)
  {
    const double d = static_cast<double>(j) - static_cast<double>(centre);
    return height * std::exp(-(d * d) / (2.0 * width * width));
  }

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
      t.library_intensity.push_back(1.0f);
      t.type.push_back(ODIA::FragmentType::Y);
      t.ordinal.push_back(static_cast<std::uint8_t>(4 + k));
      t.charge.push_back(1);
      t.loss.push_back(ODIA::LossType::None);
    }
    return lib;
  }

  /// Six transitions co-eluting at REAL_APEX; transition 0 additionally carries
  /// a much louder isolated bump at SPIKE_APEX that nothing else shares.
  ODIA::Chromatograms scriptedRun()
  {
    ODIA::Chromatograms c;
    c.axes.assign(1, std::vector<float>(CYCLES, 0.0f));
    for (std::uint32_t j = 0; j < CYCLES; ++j) { c.axes[0][j] = static_cast<float>(j); }
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
        double v = bump(j, REAL_APEX, 100.0, 3.0) + 5.0;   // shared peak + baseline
        if (k == 0) { v += bump(j, SPIKE_APEX, 900.0, 2.0); }  // the loud liar
        c.intensity[c.begin[k] + j] = static_cast<float>(v);
      }
    }
    return c;
  }

  /// Apex of the highest-scoring group, in cycles (the axis is the identity).
  double bestApex(const ODIA::PeakGroupScorer::Result& r)
  {
    double best_d = 0.0, apex = -1.0;
    bool first = true;
    for (const auto& g : r.groups)
    {
      if (first || g.dscore > best_d) { best_d = g.dscore; apex = g.apex_rt; first = false; }
    }
    return apex;
  }
} // namespace

int main()
{
  const auto lib = onePrecursorLibrary();
  const auto run = scriptedRun();

  std::printf("a co-eluting peak at cycle %u, a %gx louder single-fragment spike at %u\n",
              REAL_APEX, 9.0, SPIKE_APEX);

  {
    ODIA::PeakGroupScorer::Options o;
    o.coelution_picking = false;
    o.max_candidates = 1;
    o.min_fragments_at_apex = 1;
    const auto r = ODIA::PeakGroupScorer::score(lib, run, o);
    const double apex = bestApex(r);
    std::printf("amplitude picker: best apex at cycle %.0f\n", apex);
    check(std::fabs(apex - double(SPIKE_APEX)) <= 2.0,
          "the amplitude picker takes the spike -- this is the defect, asserted so "
          "the comparison below means something");
  }

  {
    ODIA::PeakGroupScorer::Options o;
    o.coelution_picking = true;
    o.max_candidates = 5;
    o.min_fragments_at_apex = 1;
    const auto r = ODIA::PeakGroupScorer::score(lib, run, o);
    const double apex = bestApex(r);
    std::printf("co-elution picker: best apex at cycle %.0f\n", apex);
    check(std::fabs(apex - double(REAL_APEX)) <= 2.0,
          "the co-elution picker takes the peptide, not the loudest feature");

    bool spike_absent = true;
    for (const auto& g : r.groups)
    {
      if (std::fabs(double(g.apex_rt) - double(SPIKE_APEX)) <= 2.0) { spike_absent = false; }
    }
    check(spike_absent,
          "and does not emit the single-fragment spike as a candidate at all: one "
          "fragment is not a co-elution, whatever its height");
  }

  if (failures == 0) { std::printf("all picker cases passed\n"); }
  return failures == 0 ? 0 : 1;
}
