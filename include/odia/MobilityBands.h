// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/SpectrumSource.h>

#include <cstddef>
#include <cstdint>
#include <string>
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

  // ------------------------------------------------------------------------
  // -im_bands_from_params: the band a STOCK mzpeak-convert 0.12.5 file states
  // ------------------------------------------------------------------------
  //
  // mzpeak-convert 0.12.5 writes each diaPASEF window's mobility band, but not
  // where ODIA's reader looks: it is a pair of CV parameters on the selected
  // ion (MZP:1000006 lower, MZP:1000007 upper, both 1/K0 in Vs/cm2, the
  // DiaFrameMsMsWindows scan bounds through the same scan->1/K0 function that
  // produces every peak's mobility), the plain ion_mobility_lower/upper_limit
  // columns are absent, and precursor_index is NULL on every row -- so the
  // reader attaches every ion of a frame to precursor 0 (NULL == NULL) and
  // precursor 1 gets none. Without both pieces below the band is unreachable
  // and the midpoint split is derived instead, up to 0.116 1/K0 off on
  // PXD047793 (shared/pxd/R0C_CONVERT.md).

  /// Selected-ion CV parameter accessions carrying the isolation window's
  /// inverse-reduced-mobility band (mzPeak's own namespace).
  inline constexpr const char* MZP_IM_LOWER_LIMIT = "MZP:1000006";
  inline constexpr const char* MZP_IM_UPPER_LIMIT = "MZP:1000007";

  /// One CV parameter as the reader hands it: accession and stringified value.
  struct CvValue
  {
    std::string accession;
    std::string value;
  };

  /// The band stated by @p params, if both limits are present and parse
  /// COMPLETELY as finite numbers; ordered so that @p lo <= @p hi (a writer may
  /// record them in scan order, where 1/K0 falls). Returns false -- and leaves
  /// @p lo / @p hi untouched -- when either limit is missing, repeated with a
  /// different value, or not a number, or when the two are equal: a degenerate
  /// band is not a band.
  bool mobilityBandFromParameters(const std::vector<CvValue>& params, double& lo, double& hi);

  /// How a spectrum's selected ions were paired with its precursors.
  enum class IonPairing : std::uint8_t
  {
    AsAttached,     ///< every precursor already carries exactly one ion: use it
    ByPosition,     ///< mis-attached; the k-th ion of the frame names the k-th precursor
    CountMismatch,  ///< ion count != precursor count: refused
    MzMismatch      ///< a paired ion's m/z does not name its precursor's window: refused
  };

  const char* toString(IonPairing p);

  /// Pair the selected ions of one spectrum with its precursors.
  ///
  /// @p ions_per_precursor[k] is how many ions the reader attached to precursor
  /// k, @p targets[k] its isolation-window target m/z (NaN when absent), and
  /// @p ion_mz the m/z of every ion of the spectrum IN FRAME ORDER -- precursor
  /// 0's ions first, each precursor's in the order the reader attached them,
  /// which is the file's row order (NaN when absent). On success
  /// @p ion_of[k] is the frame-order index of precursor k's ion.
  ///
  /// Pairing by position is what the converter's own reader does since 0.12.x;
  /// it is checked, not trusted: an ion whose m/z and its precursor's target are
  /// both known must agree within 0.1 Th (window spacing here is 25 Th, float
  /// noise between the two columns ~1e-4), or the spectrum is refused and
  /// @p ion_of is cleared. A refusal leaves the caller on the legacy reading.
  IonPairing pairIonsByPosition(const std::vector<std::size_t>& ions_per_precursor,
                                const std::vector<double>& targets,
                                const std::vector<double>& ion_mz,
                                std::vector<std::size_t>& ion_of);

} // namespace ODIA
