// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/Library.h>
#include <odia/PeakGroupScorer.h>

#include <string>

namespace ODIA
{

  /// `-out_fragvec`: the 78 rung-(i) fragment columns, one row per scored
  /// candidate, in `-out`'s exact row order.
  ///
  /// A SEPARATE FILE, never extra columns on `-out`. The point of this export
  /// is that it is provably score-neutral -- the arm's gate is that `-out`'s
  /// sha256 with the flag ON equals the flag-off run's -- and columns appended
  /// to `-out` would make that gate unstateable. (`-fragvec_scores` does put
  /// all 78 of them on `-out`, as var_fv_* or var_fvperm_*, because there they
  /// ARE scored; this export stays the raw contract either way, in the same
  /// row order, so the two can be compared cell for cell.)
  ///
  /// Keyed (Precursor.Id, Decoy, Ordinal) so the join back to `-out` is
  /// checkable rather than assumed. `Ordinal` is the 0-based position of the
  /// row within its own (Precursor.Id, Decoy) block, which is contiguous in
  /// both files because `finish()` sorts by precursor first.
  ///
  /// Written at precision NINE, deliberately NOT the 6 `-out` and `-out_chrom`
  /// use: the contract stores these as float32 and the arm's gate A4 compares
  /// them cell for cell, so the file has to round-trip a float32 exactly.
  /// `%.6g` does not -- measured on the wf_v64 fixture it costs a uniform
  /// ~5e-06 RELATIVE error -- and 9 is FLT_DECIMAL_DIG, the shortest precision
  /// that round-trips every float32. See the comment at the `number` call.
  ///
  /// Throws when `Result::fragvec` is not `N_FRAGVEC` floats per group: that
  /// means something dropped groups without carrying the rows with them, and a
  /// silently misaligned export is worse than none. Note what that guard does
  /// NOT catch: two rows PAIRED with the wrong groups keep the count right. It
  /// cannot happen today because the only append site is the statement before
  /// the only `groups.push_back` and `Session::add` runs on one thread
  /// (PeakGroupScorer.h, "add runs on one thread"); a future parallel emit loop
  /// would have to carry an explicit row index instead.
  void writeFragvecTsv(const std::string& path, const Library& library,
                       const PeakGroupScorer::Result& scored);

} // namespace ODIA
