// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/Ms1Traces.h>
#include <limits>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

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
    for (std::size_t i = 0; i < rows(); ++i)
    {
      for (std::size_t b = 0; b < bins_; ++b)
      { if (values_[i * bins_ + b] > 0.0f) { ++live; break; } }
    }
    o.precision(1);
    // Rows, not precursors(): the loop above walks stored rows, and with a
    // library->row map those are fewer than the library.
    const std::size_t rows = this->rows();
    o << std::fixed << "MS1 traces: " << bins_ << " bins over " << rows
      << " precursors, " << live << " with signal ("
      << (rows ? 100.0 * live / rows : 0.0) << "%), "
      << footprintBytes() / 1048576.0 << " MiB";
    if (!row_of_.empty())
    {
      o << "; " << dropped_rows_ << " of " << row_of_.size()
        << " library precursors given no row (no isolation window covers them)";
    }
    return o.str();
  }

  Ms1Traces Ms1Traces::build(const Library& library, SpectrumSource& source,
                             double fragment_ppm, double im_window,
                             double ppm_offset, double* observed_ppm_median,
                             double isotope_offset_da,
                             const std::vector<std::uint8_t>* keep,
                             std::vector<std::uint32_t>* kept_indices,
                             const std::vector<std::uint8_t>* library_rows)
  {
    if (keep != nullptr && library_rows != nullptr)
    {
      throw std::invalid_argument(
        "Ms1Traces::build: `keep` (row-indexed export) and `library_rows` (library-indexed "
        "scorer map) are different contracts for at() and cannot be combined");
    }
    Ms1Traces out;
    const auto& ms1 = source.ms1Spectra();
    if (ms1.empty()) { return out; }

    const auto& p = library.precursors();
    const std::size_t np = library.precursorCount();
    out.bins_ = ms1.size();
    out.times_.reserve(ms1.size());
    for (const auto& s : ms1) { out.times_.push_back(static_cast<float>(s.retention_time)); }

    // Row assignment. Without a mask each library precursor owns row i (the
    // scorer's contract). With a mask, rows are assigned in ascending library
    // index over the kept precursors, and the caller gets that order back via
    // kept_indices -- the writer reconstructs ids from it, so the mapping is
    // never inferred twice.
    std::vector<std::uint32_t> row;
    std::size_t rows = np;
    if (keep != nullptr)
    {
      row.assign(np, UINT32_MAX);
      std::uint32_t r = 0;
      for (std::size_t i = 0; i < np; ++i)
      {
        if (i < keep->size() && (*keep)[i]) { row[i] = r++; }
      }
      rows = r;
      if (kept_indices != nullptr)
      {
        kept_indices->clear();
        kept_indices->reserve(rows);
        for (std::size_t i = 0; i < np; ++i)
        { if (row[i] != UINT32_MAX) { kept_indices->push_back(static_cast<std::uint32_t>(i)); } }
      }
    }
    // -ms1_drop_uncovered: the same row assignment as `keep`, but the map is
    // KEPT on the object so at() can go on taking library indices.
    if (library_rows != nullptr)
    {
      out.row_of_.assign(np, UINT32_MAX);
      std::uint32_t r = 0;
      for (std::size_t i = 0; i < np; ++i)
      {
        if (i < library_rows->size() && (*library_rows)[i]) { out.row_of_[i] = r++; }
      }
      rows = r;
      out.dropped_rows_ = np - r;
      out.dropped_reads_ = std::make_shared<std::atomic<std::size_t>>(0);
    }
    out.values_.assign(rows * out.bins_, 0.0f);

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
      if (keep != nullptr && row[i] == UINT32_MAX) { continue; }
      if (library_rows != nullptr && out.row_of_[i] == UINT32_MAX) { continue; }
      // CALIBRATED, like the fragment axis. This matched on the library's
      // THEORETICAL m/z with a symmetric window and no offset, while the
      // fragment extractor was centred on the fitted deviation -- on IH1 that
      // is -10.0108 ppm against a +/-10 ppm half-width, so a precursor whose
      // MS1 error resembles its MS2 error sat at the window EDGE and a weak one
      // fell out entirely. That matters because ms1_coelution is the main
      // evidence for calling a precursor ABSENT (median -0.093 for the 10,736
      // DIA-NN precursors we reject, against -0.124 for the 801,458 bulk
      // non-identifications and 0.473 for accepted ones), and absence cannot be
      // concluded from an uncalibrated measurement.
      //
      // The isotope offset is applied BEFORE the calibration scaling, and per
      // this precursor's own charge: the M+k target is `mz + k*dm/z`, and the
      // instrument's relative (ppm) error then applies to that target as it
      // does to any mass.
      const int z = p.charge[i] > 0 ? static_cast<int>(p.charge[i]) : 1;
      const double mz = (fromFixed(p.mz[i]) + isotope_offset_da / z) *
                        (1.0 + ppm_offset * 1e-6);
      const std::uint32_t slot = keep != nullptr           ? row[i]
                                 : library_rows != nullptr ? out.row_of_[i]
                                                           : static_cast<std::uint32_t>(i);
      if (mz > 0.0) { idx.push_back({mz, slot, p.im[i]}); }
    }
    std::sort(idx.begin(), idx.end(),
              [](const Target& a, const Target& b) { return a.mz < b.mz; });

    // CAPPED. The first version pushed one double per (peak, target) match over
    // the whole run and reached 591 GB RSS against v3's 116 GB peak, blowing
    // through -live_memory_gb on a shared node. The median of a bounded prefix
    // is the same number to far more precision than it is worth: this is a
    // diagnostic, not a fit.
    static constexpr std::size_t RESID_CAP = 1u << 21;   // 2M samples, 16 MB
    std::vector<double> resid;
    if (observed_ppm_median != nullptr) { resid.reserve(RESID_CAP); }
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
            if (observed_ppm_median != nullptr && resid.size() < RESID_CAP)
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
