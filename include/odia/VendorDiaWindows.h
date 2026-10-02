#pragma once

#include <string>
#include <vector>

namespace ODIA
{
  /// One diaPASEF window exactly as the instrument method states it.
  ///
  /// The mobility bounds are the point of this struct. Everything else about
  /// the scheme can be read from the spectra, but the per-window 1/K0 band
  /// cannot: the converted file carries each window's mobility POSITION and no
  /// bounds at all -- `ion_mobility_lower_limit` and `_upper_limit` are absent
  /// from the schema, not merely null. Deriving the split as the midpoint
  /// between two positions is exact only for equally wide windows, and on IH1
  /// ten of twelve groups are not (doc/69: the derived boundary is off by
  /// +0.0498 at group 1 and -0.0845 at group 12, against a 0.059 mobility
  /// window).
  struct VendorDiaWindow
  {
    double one_over_k0_start = 0.0;  ///< inclusive lower 1/K0 bound
    double one_over_k0_end = 0.0;    ///< inclusive upper 1/K0 bound
    double isolation_mz = 0.0;       ///< window centre, matches the spectra's target
    double isolation_width = 0.0;
  };

  /// Read the vendor's stated window table out of an .mzpeak.
  ///
  /// The table lives in the instrument method the conversion embeds --
  /// `vendor/<method>.m/diaSettings.diasqlite`, table `DiaWindowsSpecification`
  /// -- so this needs no vendor SDK, no .d beside the converted file, and no
  /// re-conversion. It is the same archive the reader already opens.
  ///
  /// Returns EMPTY rather than throwing when the file carries no vendor method,
  /// carries no such table, or is not an archive at all. A run on data from
  /// another instrument is the normal case for that, not an error, and the
  /// caller falls back to deriving the bands.
  std::vector<VendorDiaWindow> readVendorDiaWindows(const std::string& mzpeak_path);
}
