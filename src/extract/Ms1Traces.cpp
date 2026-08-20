// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/Ms1Traces.h>
#include <limits>

#include <algorithm>
#include <cmath>
#include <sstream>

namespace ODIA
{
  std::size_t Ms1Traces::binFor(double rt) const
  {
    if (times_.empty()) { return 0; }
    const auto it = std::lower_bound(times_.begin(), times_.end(), static_cast<float>(rt));
    std::size_t b = static_cast<std::size_t>(it - times_.begin());
    if (b >= times_.size()) { return times_.size() - 1; }
    if (b > 0 && std::abs(times_[b - 1] - rt) < std::abs(times_[b] - rt)) { --b; }
    return b;
  }

  std::string Ms1Traces::describe() const
  {
    std::ostringstream o;
    if (empty()) { return "MS1 traces: none (the run carries no MS1)"; }
    std::size_t live = 0;
    for (std::size_t i = 0; i < precursors(); ++i)
    {
      for (std::size_t b = 0; b < bins_; ++b)
      { if (values_[i * bins_ + b] > 0.0f) { ++live; break; } }
    }
    o.precision(1);
    o << std::fixed << "MS1 traces: " << bins_ << " bins over " << precursors()
      << " precursors, " << live << " with signal ("
      << (precursors() ? 100.0 * live / precursors() : 0.0) << "%), "
      << footprintBytes() / 1048576.0 << " MiB";
    return o.str();
  }

  Ms1Traces Ms1Traces::build(const Library& library, SpectrumSource& source,
                             double fragment_ppm, double im_window,
                             double ppm_offset, double* observed_ppm_median)
  {
    Ms1Traces out;
    const auto& ms1 = source.ms1Spectra();
    if (ms1.empty()) { return out; }

    const auto& p = library.precursors();
    const std::size_t np = library.precursorCount();
    out.bins_ = ms1.size();
    out.times_.reserve(ms1.size());
    for (const auto& s : ms1) { out.times_.push_back(static_cast<float>(s.retention_time)); }
    out.values_.assign(np * out.bins_, 0.0f);

    // Search the sorted LIBRARY side and iterate the peaks: SpectrumSource
    // documents that a peak array is not ascending in m/z (a mobility frame
    // concatenates its TIMS scans), and a binary search over it "does not fail
    // loudly -- it returns near-zero matches", which is indistinguishable from
    // an ion that is not there.
    struct Target { double mz; std::uint32_t slot; float im; };
    std::vector<Target> idx;
    idx.reserve(np);
    for (std::size_t i = 0; i < np; ++i)
    {
      // CALIBRATED, like the fragment axis. This matched on the library's
      // THEORETICAL m/z with a symmetric window and no offset, while the
      // fragment extractor was centred on the fitted deviation -- on S08 that
      // is -10.0108 ppm against a +/-10 ppm half-width, so a precursor whose
      // MS1 error resembles its MS2 error sat at the window EDGE and a weak one
      // fell out entirely. That matters because ms1_coelution is the main
      // evidence for calling a precursor ABSENT (median -0.093 for the 10,736
      // DIA-NN precursors we reject, against -0.124 for the 801,458 bulk
      // non-identifications and 0.473 for accepted ones), and absence cannot be
      // concluded from an uncalibrated measurement.
      const double mz = fromFixed(p.mz[i]) * (1.0 + ppm_offset * 1e-6);
      if (mz > 0.0) { idx.push_back({mz, static_cast<std::uint32_t>(i), p.im[i]}); }
    }
    std::sort(idx.begin(), idx.end(),
              [](const Target& a, const Target& b) { return a.mz < b.mz; });

    std::vector<double> resid;
    if (observed_ppm_median != nullptr) { resid.reserve(1u << 20); }
    std::vector<SpectrumPeaks> block;
    const std::size_t STEP = 64;
    for (std::size_t b = 0; b < ms1.size(); b += STEP)
    {
      const std::size_t e = std::min(b + STEP, ms1.size());
      source.ms1Peaks(b, e, block);
      for (std::size_t s = 0; s < block.size(); ++s)
      {
        const auto& sp = block[s];
        const bool gated = im_window > 0.0 && sp.ion_mobility.size() == sp.mz.size();
        for (std::size_t k = 0; k < sp.mz.size(); ++k)
        {
          const double m = sp.mz[k], tol = m * fragment_ppm * 1e-6;
          auto it = std::lower_bound(idx.begin(), idx.end(), m - tol,
                                     [](const Target& a, double v) { return a.mz < v; });
          for (; it != idx.end() && it->mz <= m + tol; ++it)
          {
            if (std::abs(it->mz - m) > it->mz * fragment_ppm * 1e-6) { continue; }
            // Same rule as the MS2 match loop, deliberately: skip the gate when
            // EITHER side is unknown, because absent information is not evidence
            // of mismatch. Testing `abs(NaN - x) <= w` is false, so a precursor
            // with no library 1/K0 had every MS1 peak rejected and came out with
            // a NaN MS1_COELUTION -- while its MS2 side was extracted ungated.
            if (gated && !std::isnan(static_cast<double>(it->im)))
            {
              const double d = std::abs(static_cast<double>(sp.ion_mobility[k]) -
                                        static_cast<double>(it->im));
              if (!(d <= im_window)) { continue; }
            }
            if (observed_ppm_median != nullptr)
            {
              // Residual against the CALIBRATED target, so a correct offset
              // centres this on 0 and a wrong one does not.
              resid.push_back((m - it->mz) / it->mz * 1e6);
            }
            float& c = out.values_[it->slot * out.bins_ + (b + s)];
            // Max, not sum: a mobility-merged frame holds the same ion in
            // several scans, and summing would make the trace a function of how
            // many scans it spans rather than of how much ion is present.
            c = std::max(c, sp.intensity[k]);
          }
        }
      }
    }
    if (observed_ppm_median != nullptr)
    {
      if (resid.empty()) { *observed_ppm_median = std::numeric_limits<double>::quiet_NaN(); }
      else
      {
        const std::size_t h = resid.size() / 2;
        std::nth_element(resid.begin(), resid.begin() + h, resid.end());
        *observed_ppm_median = resid[h];
      }
    }
    return out;
  }
} // namespace ODIA
