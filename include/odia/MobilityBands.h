// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/SpectrumSource.h>

#include <cstdint>
#include <vector>

namespace ODIA
{

  /// A selected ion as a diaPASEF file mis-attaches it: an m/z that NAMES an
  /// isolation window (it equals the window's centre) and the mobility position
  /// at which that window was acquired.
  ///
  /// Not a band. The position is a single 1/K0; the band has to be derived from
  /// the positions of the neighbouring windows, which is what this header is
  /// for.
  struct MobilityPosition
  {
    double mz = 0.0;
    double im = 0.0;
  };

  /// What became of a derivation. Anything but Derived leaves the windows
  /// exactly as they were.
  enum class MobilityBandResult : std::uint8_t
  {
    Derived,          ///< every window got a band
    NotNeeded,        ///< fewer than two windows: nothing to separate
    TooFewIons,       ///< fewer positions than windows
    Unattributable,   ///< a window's nearest position is too far to name it
    Ambiguous         ///< two attributions are equally good; guessing is not allowed here
  };

  const char* toString(MobilityBandResult r);

  /// Give each of @p windows the mobility band implied by @p ions.
  ///
  /// Each window is attributed the position whose m/z names it, and the band is
  /// the midpoint between its position and its neighbours'. The outermost
  /// windows keep an open end, because the frame says nothing about where the
  /// mobility range stops. A window that already carries a stated band keeps it
  /// -- a derivation never overrides what the file says.
  ///
  /// The attribution is by GLOBALLY smallest distance, not first-come: pairing
  /// each window in turn with its own nearest unclaimed position depends on the
  /// order the windows happen to be listed in, and can hand window A the
  /// position that names window B when A had an alternative and B did not.
  /// Taking the smallest distance first makes each pair mutually nearest among
  /// what is left, so the result is a property of the numbers.
  ///
  /// This is a guess about the instrument, and a wrong guess rejects real
  /// signal, so it is refused rather than approximated: a position further from
  /// a window's centre than half its width plus 1 Th does not name it, and two
  /// equally good attributions are not resolved by picking one. Refusing leaves
  /// the windows with the full mobility range, which admits everything -- the
  /// failure mode that loses nothing.
  MobilityBandResult deriveMobilityBands(const std::vector<MobilityPosition>& ions,
                                         std::vector<IsolationWindow>& windows);

} // namespace ODIA
