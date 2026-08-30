// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/ChromatogramExtractor.h>
#include <odia/Library.h>
#include <odia/Ms1Traces.h>

#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

namespace ODIA
{

  /// The extracted chromatograms, long format: one row per point.
  ///
  /// Deliberately not one row per transition with packed arrays. This file is
  /// what the scoring stage and any external check will read, and a long table
  /// is what every tool that might read it -- pandas, R, DuckDB -- takes
  /// without a parser of its own.
  ///
  /// Lives here rather than inside the tool because it is 33% of Phase 2 and
  /// nothing could reach it to test: at 9,522 precursors it writes 7.16 GB in
  /// ~301 s, which is ~24 MB/s -- `ofstream <<`'s measured 25.4 MB/s and not
  /// the NVMe's 1,042 MB/s. `odia_tsv_writers` now diffs its output byte for
  /// byte against the ofstream form it replaces.
  /// @p keep, when given, restricts the output to precursors whose
  /// reconstructed Precursor.Id (modified sequence + charge, the file's own
  /// join key) is in the set. A listed id emits BOTH its target and its decoy
  /// block -- the Decoy column separates them downstream, exactly as the
  /// header comment above warns.
  void writeChromatogramTsv(const std::string& path, const Library& library,
                            const Chromatograms& chromatograms,
                            const std::unordered_set<std::string>* keep = nullptr);

  /// First column of a TSV as a set, for -out_chrom_ids. Skips blank lines and
  /// a header line whose first cell is exactly "Precursor.Id".
  std::unordered_set<std::string> readPrecursorIdList(const std::string& path);

  /// MS1 M/M+1/M+2 traces for a cohort, long format:
  ///   Precursor.Id  Decoy  Isotope  Precursor.Mz  RT  Intensity
  /// One trace per (kept precursor, isotope); only non-zero runs are written,
  /// padded with one flanking zero on each side so peak boundaries survive
  /// without the file growing with the empty remainder of the run.
  /// @p kept  library precursor indices in ROW ORDER of the traces (the
  ///          `kept_indices` the cohort-restricted Ms1Traces::build filled).
  void writeMs1TracesTsv(const std::string& path, const Library& library,
                         const Ms1Traces& iso0, const Ms1Traces& iso1,
                         const Ms1Traces& iso2,
                         const std::vector<std::uint32_t>& kept);

} // namespace ODIA
