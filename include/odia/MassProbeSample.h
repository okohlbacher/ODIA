// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// The fragment mass probe's calibration sample, drawn from SCORED peak groups
// rather than from the library (-mass_probe_source scored).
//
// The library stride (RunProbe::samplePrecursors) samples precursors the run
// mostly does not contain: on a 6.8 M-target predicted library ~0.1-1.4% are
// present, so the mass gate reads interference (memory
// odia-mass-probe-needs-present-sample; R0C_FIXTURE dilution series). Pass 1
// has already found where the present ones are. This takes:
//
//   targets  the best row of every target precursor whose best q-value is at or
//            below `q` -- confidently scored, i.e. present -- with its apex RT
//            and observed 1/K0;
//   decoys   the matching control: the same NUMBER of decoy precursors, their
//            best rows by d-score, probed the same way. Never fed to the gate;
//            it is reported, as the check that the selection itself does not
//            manufacture a mass mode (a decoy best row was chosen by the same
//            scorer through the same window).
//
// Both lists are capped by the stride the library sample uses, applied to the
// precursor-index-sorted list, so a cap spreads over the library's m/z rather
// than keeping the highest scores (a score-ordered cut is a selection on the
// very evidence being calibrated -- memory odia-accepted-groups-are-a-biased-sample).
// OpenDIAlyzer passes cap 0 and lets MassCalibration::collect apply the same
// stride (Options::max_precursors) AFTER its cycle cap, so the precursor cap is
// spent on entries the probe's visited cycles can reach.
// Deterministic: every tie is broken by group index.

#pragma once

#include <odia/PeakGroupScorer.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>

namespace ODIA::MassProbeSample
{
  /// Parallel arrays: library index, apex RT (run seconds), observed 1/K0
  /// (NaN when the run had none), sorted by library index.
  struct Entries
  {
    std::vector<std::uint32_t> precursor;
    std::vector<float> rt;
    std::vector<float> im;
    std::size_t size() const { return precursor.size(); }
  };

  struct Selection
  {
    Entries targets;
    Entries decoys;
    /// Target precursors whose best row passed `q`, before the cap.
    std::size_t eligible_targets = 0;
    /// Decoy precursors with any row, before matching and the cap.
    std::size_t available_decoys = 0;
  };

  namespace detail
  {
    /// Keep `want` of `idx` (already sorted) by the library sample's stride.
    inline std::vector<std::size_t> stride(const std::vector<std::size_t>& idx, std::size_t want)
    {
      if (want == 0 || want >= idx.size()) { return idx; }
      std::vector<std::size_t> out;
      out.reserve(want);
      for (std::size_t k = 0; k < want; ++k) { out.push_back(idx[k * idx.size() / want]); }
      return out;
    }
  } // namespace detail

  /// @param q    best-row q-value bar for a target to count as present
  /// @param cap  at most this many entries per list; 0 means no cap
  inline Selection fromGroups(const std::vector<PeakGroupScorer::PeakGroup>& groups, double q,
                              std::size_t cap)
  {
    // Best row per precursor. Targets: lowest q, then highest d-score, then
    // lowest group index. Decoys: highest d-score, then lowest group index.
    std::unordered_map<std::uint32_t, std::size_t> best;
    best.reserve(groups.size() / 4 + 16);
    const auto better = [&](std::size_t a, std::size_t b) {
      const auto& ga = groups[a];
      const auto& gb = groups[b];
      if (!ga.decoy)
      {
        if (ga.qvalue != gb.qvalue) { return ga.qvalue < gb.qvalue; }
      }
      if (ga.dscore != gb.dscore) { return ga.dscore > gb.dscore; }
      return a < b;
    };
    for (std::size_t g = 0; g < groups.size(); ++g)
    {
      const auto& pg = groups[g];
      // A row without a finite score was not scored; it is nobody's best row.
      if (!std::isfinite(pg.dscore)) { continue; }
      auto it = best.find(pg.precursor);
      if (it == best.end()) { best.emplace(pg.precursor, g); }
      else if (better(g, it->second)) { it->second = g; }
    }

    std::vector<std::size_t> tgt, dec;
    for (const auto& kv : best)
    {
      const auto& pg = groups[kv.second];
      if (pg.decoy) { dec.push_back(kv.second); }
      else if (pg.qvalue <= q) { tgt.push_back(kv.second); }
    }

    Selection sel;
    sel.eligible_targets = tgt.size();
    sel.available_decoys = dec.size();

    // The control: as many decoys as there are eligible targets, the
    // best-scoring ones -- the decoys that look most like an identification.
    std::sort(dec.begin(), dec.end(), [&](std::size_t a, std::size_t b) {
      if (groups[a].dscore != groups[b].dscore) { return groups[a].dscore > groups[b].dscore; }
      return a < b;
    });
    if (dec.size() > tgt.size()) { dec.resize(tgt.size()); }

    const auto by_precursor = [&](std::size_t a, std::size_t b) {
      return groups[a].precursor < groups[b].precursor;
    };
    std::sort(tgt.begin(), tgt.end(), by_precursor);
    std::sort(dec.begin(), dec.end(), by_precursor);

    const auto fill = [&](const std::vector<std::size_t>& idx, Entries& e) {
      for (const std::size_t g : detail::stride(idx, cap))
      {
        e.precursor.push_back(groups[g].precursor);
        e.rt.push_back(groups[g].apex_rt);
        e.im.push_back(groups[g].observed_im);
      }
    };
    fill(tgt, sel.targets);
    fill(dec, sel.decoys);
    return sel;
  }
} // namespace ODIA::MassProbeSample
