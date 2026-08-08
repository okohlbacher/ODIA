// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/OpenSwathPicker.h>

#include <OpenMS/ANALYSIS/OPENSWATH/PeakPickerChromatogram.h>
#include <OpenMS/KERNEL/MSChromatogram.h>

#include <algorithm>
#include <cmath>

namespace ODIA
{
  std::vector<OpenSwathCandidate> pickOpenSwath(const std::vector<double>& total,
                                                const std::vector<float>& rt,
                                                double sn_threshold, bool use_gauss,
                                                double peak_width,
                                                std::size_t max_candidates)
  {
    std::vector<OpenSwathCandidate> out;
    const std::size_t n = std::min(total.size(), rt.size());
    // The picker smooths over a window and then looks for a maximum with
    // neighbours on both sides; below a handful of points there is nothing for
    // it to do and OpenMS is entitled to be unhappy about it.
    if (n < 7) { return out; }

    OpenMS::MSChromatogram chrom;
    chrom.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
    {
      OpenMS::ChromatogramPeak p;
      p.setRT(static_cast<double>(rt[i]));
      // OpenMS peak intensities are float and non-negative by contract.
      p.setIntensity(static_cast<float>(std::max(0.0, total[i])));
      chrom.push_back(p);
    }

    OpenMS::PeakPickerChromatogram picker;
    OpenMS::Param p = picker.getDefaults();
    p.setValue("signal_to_noise", sn_threshold);
    p.setValue("use_gauss", use_gauss ? "true" : "false");
    if (peak_width > 0.0)
    {
      p.setValue("peak_width", peak_width);
    }
    // Let it report rather than throw on a chromatogram it dislikes.
    if (p.exists("write_sn_log_messages")) { p.setValue("write_sn_log_messages", "false"); }
    picker.setParameters(p);

    OpenMS::MSChromatogram picked;
    try { picker.pickChromatogram(chrom, picked); }
    catch (...) { return out; }
    if (picked.empty()) { return out; }

    // Boundaries come back as float data arrays named leftWidth/rightWidth, in
    // RETENTION TIME. Mapped back to cycle indices here, because everything
    // downstream in ODIA indexes by cycle.
    const OpenMS::MSChromatogram::FloatDataArray* lw = nullptr;
    const OpenMS::MSChromatogram::FloatDataArray* rw = nullptr;
    for (const auto& fa : picked.getFloatDataArrays())
    {
      if (fa.getName() == "leftWidth") { lw = &fa; }
      else if (fa.getName() == "rightWidth") { rw = &fa; }
    }

    const auto cycle_of = [&](double t) {
      const auto it = std::lower_bound(rt.begin(), rt.begin() + n, static_cast<float>(t));
      std::size_t i = static_cast<std::size_t>(it - rt.begin());
      if (i >= n) { return n - 1; }
      if (i > 0 && std::abs(rt[i - 1] - t) < std::abs(rt[i] - t)) { --i; }
      return i;
    };

    for (std::size_t k = 0; k < picked.size(); ++k)
    {
      OpenSwathCandidate c;
      c.apex = cycle_of(picked[k].getRT());
      c.apex_value = static_cast<double>(picked[k].getIntensity());
      c.left = (lw && k < lw->size()) ? cycle_of((*lw)[k]) : c.apex;
      c.right = (rw && k < rw->size()) ? cycle_of((*rw)[k]) : c.apex;
      if (c.left > c.apex) { c.left = c.apex; }
      if (c.right < c.apex) { c.right = c.apex; }
      if (c.right >= n) { c.right = n - 1; }
      out.push_back(c);
    }

    // OpenSWATH returns peaks in retention-time order; ODIA's convention is
    // best-first, and max_candidates must cut the WEAKEST rather than the
    // latest-eluting.
    std::stable_sort(out.begin(), out.end(),
                     [](const OpenSwathCandidate& a, const OpenSwathCandidate& b)
                     { return a.apex_value > b.apex_value; });
    if (out.size() > max_candidates) { out.resize(max_candidates); }
    return out;
  }
} // namespace ODIA
