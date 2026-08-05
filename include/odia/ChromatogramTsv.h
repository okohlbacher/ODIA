// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/ChromatogramExtractor.h>
#include <odia/Library.h>

#include <string>

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
  void writeChromatogramTsv(const std::string& path, const Library& library,
                            const Chromatograms& chromatograms);

} // namespace ODIA
