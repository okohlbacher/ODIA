// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/MassRecalibration.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>

namespace ODIA
{
  namespace
  {
    double medianOf(std::vector<double>& v)
    {
      if (v.empty()) { return 0.0; }
      const std::size_t h = v.size() / 2;
      std::nth_element(v.begin(), v.begin() + h, v.end());
      const double a = v[h];
      if (v.size() % 2) { return a; }
      std::nth_element(v.begin(), v.begin() + h - 1, v.end());
      return 0.5 * (a + v[h - 1]);
    }

    /// Weighted least squares of ppm on log(mz/ref), returning (intercept, slope).
    ///
    /// Weighted by intensity: a bright fragment's centroid is better determined
    /// than a dim one's, and on DIA data the dim end is where interference
    /// lives. Falls back to an unweighted median when the design is degenerate
    /// -- all one mass, or a single point -- because a slope fitted through a
    /// single m/z is not a slope.
    std::pair<double, double> fitLine(const std::vector<double>& x,
                                      const std::vector<double>& y,
                                      const std::vector<double>& w, bool want_slope)
    {
      double sw = 0.0, sx = 0.0, sy = 0.0;
      for (std::size_t i = 0; i < x.size(); ++i)
      { sw += w[i]; sx += w[i] * x[i]; sy += w[i] * y[i]; }
      if (!(sw > 0.0)) { return {0.0, 0.0}; }
      const double mx = sx / sw, my = sy / sw;
      if (!want_slope) { return {my, 0.0}; }
      double sxx = 0.0, sxy = 0.0;
      for (std::size_t i = 0; i < x.size(); ++i)
      {
        const double dx = x[i] - mx;
        sxx += w[i] * dx * dx;
        sxy += w[i] * dx * (y[i] - my);
      }
      if (!(sxx > 0.0)) { return {my, 0.0}; }
      const double slope = sxy / sxx;
      return {my - slope * mx, slope};
    }
  } // namespace

  double MassRecalibration::ppmAt(double mz, double rt) const
  {
    if (!fitted_ || blocks_.empty() || !(mz > 0.0)) { return 0.0; }
    const double l = std::log(mz / ref_mz_);
    // Blocks are contiguous and ascending; outside the fitted span a precursor
    // takes the nearest block rather than nothing. Extrapolating the TREND
    // would be worse: the ends of a gradient are where the fit is thinnest.
    const Block* chosen = &blocks_.front();
    for (const auto& b : blocks_)
    {
      if (rt >= b.rt_low && rt <= b.rt_high) { chosen = &b; break; }
      if (rt > b.rt_high) { chosen = &b; }
    }
    const double raw = chosen->fitted ? chosen->intercept_ppm + chosen->slope_ppm * l
                                      : global_intercept_ + global_slope_ * l;
    // The COMPOSITE is bounded, not merely each term.
    //
    // Clamping intercept and slope separately leaves the sum unbounded: with a
    // reference near 1000 Th and a slope at its own limit, m/z 150 is 1.9
    // e-folds away and the slope term alone reaches ~28 ppm on top of the
    // intercept. MobilityCalibration had exactly this defect and it was fixed
    // hours before this class was written; the unit test below is what caught
    // it being repeated.
    return std::clamp(raw, -max_correction_ppm_, max_correction_ppm_);
  }

  double MassRecalibration::maxAbsCorrection() const
  {
    if (!fitted_) { return 0.0; }
    // Evaluated at the extremes of a generous m/z range rather than analytically:
    // the model is piecewise, and the caller needs a bound it can trust to size
    // a search window, not a tight one.
    // Evaluated THROUGH ppmAt, not from the coefficients: the composite is
    // clamped there, and a bound computed from the raw coefficients would
    // overstate what is applied -- which for a caller sizing a search window
    // means a needlessly wide search, and for this class's own test meant a
    // failure that looked like a missing clamp when the clamp was present.
    double worst = 0.0;
    for (const auto& b : blocks_)
    {
      const double rt = b.fitted ? 0.5 * (std::max(b.rt_low, -1e6) + std::min(b.rt_high, 1e6))
                                 : 0.0;
      for (const double mz : {150.0, 400.0, 1000.0, 2000.0})
      {
        worst = std::max(worst, std::abs(ppmAt(mz, rt)));
      }
    }
    return worst;
  }

  std::string MassRecalibration::describe() const
  {
    std::ostringstream o;
    o.precision(2);
    o << std::fixed;
    if (!fitted_) { return "fragment mass recalibration: NOT FITTED"; }
    std::size_t live = 0;
    for (const auto& b : blocks_) { if (b.fitted) { ++live; } }
    o << "fragment mass recalibration: " << live << " of " << blocks_.size()
      << " retention-time blocks fitted, linear in log m/z about " << ref_mz_ << " Th"
      << " (global " << std::showpos << global_intercept_ << " ppm, slope "
      << global_slope_ << std::noshowpos << " per e-fold)\n";
    for (const auto& b : blocks_)
    {
      const bool first = (&b == &blocks_.front()), last = (&b == &blocks_.back());
      o << "  ";
      if (first) { o << "start"; } else { o << b.rt_low; }
      o << "-";
      if (last) { o << "end"; } else { o << b.rt_high; }
      o << " s: ";
      if (!b.fitted) { o << "only " << b.anchors << " anchors, inherits the global fit\n"; }
      else
      {
        o << std::showpos << b.intercept_ppm << " ppm" << std::noshowpos
          << ", slope " << std::showpos << b.slope_ppm << std::noshowpos
          << " per e-fold, from " << b.anchors << " anchors\n";
      }
    }
    return o.str();
  }

  MassRecalibration MassRecalibration::fit(const std::vector<MassResidual>& residuals,
                                           const Options& opt)
  {
    MassRecalibration out;

    struct P { double rt, l, ppm, w; };
    std::vector<P> pts;
    pts.reserve(residuals.size());
    std::vector<double> mzs;
    for (const auto& r : residuals)
    {
      if (r.decoy) { continue; }                 // a control cannot inform a correction
      if (!(r.mz > 0.0) || !std::isfinite(r.ppm)) { continue; }
      mzs.push_back(r.mz);
    }
    if (mzs.size() < opt.min_per_block) { return out; }
    out.ref_mz_ = medianOf(mzs);
    if (!(out.ref_mz_ > 0.0)) { return out; }

    for (const auto& r : residuals)
    {
      if (r.decoy || !(r.mz > 0.0) || !std::isfinite(r.ppm)) { continue; }
      pts.push_back({r.rt, std::log(r.mz / out.ref_mz_), r.ppm,
                     std::max(1.0, static_cast<double>(r.intensity))});
    }
    if (pts.size() < opt.min_per_block) { return out; }

    // Trim by MAD before anything is fitted. The anchor set is confident but
    // not pure, and least squares has no defence against a single wild point.
    {
      std::vector<double> p;
      p.reserve(pts.size());
      for (const auto& q : pts) { p.push_back(q.ppm); }
      const double m = medianOf(p);
      std::vector<double> ad;
      ad.reserve(pts.size());
      for (const auto& q : pts) { ad.push_back(std::abs(q.ppm - m)); }
      const double mad = medianOf(ad) * 1.4826;
      if (mad > 0.0)
      {
        const double lim = opt.trim_mads * mad;
        pts.erase(std::remove_if(pts.begin(), pts.end(),
                                 [&](const P& q) { return std::abs(q.ppm - m) > lim; }),
                  pts.end());
      }
    }
    if (pts.size() < opt.min_per_block) { return out; }

    std::sort(pts.begin(), pts.end(), [](const P& a, const P& b) { return a.rt < b.rt; });

    // The global fit first: every block that cannot support its own falls back
    // to this rather than to zero, so a thin block is still corrected.
    {
      std::vector<double> x, y, w;
      x.reserve(pts.size()); y.reserve(pts.size()); w.reserve(pts.size());
      for (const auto& q : pts) { x.push_back(q.l); y.push_back(q.ppm); w.push_back(q.w); }
      out.max_correction_ppm_ = opt.max_correction_ppm;
      const auto g = fitLine(x, y, w, pts.size() >= opt.min_for_slope);
      out.global_intercept_ = std::clamp(g.first, -opt.max_correction_ppm, opt.max_correction_ppm);
      out.global_slope_ = std::clamp(g.second, -opt.max_slope_ppm, opt.max_slope_ppm);
    }

    // Equal-COUNT blocks. A gradient's identifications are not uniform in time,
    // and equal-width blocks would fit the sparse ends from almost nothing.
    const std::size_t nb = std::max<std::size_t>(1, opt.blocks);
    const std::size_t per = std::max<std::size_t>(1, pts.size() / nb);
    for (std::size_t b = 0; b < nb; ++b)
    {
      const std::size_t lo = b * per;
      std::size_t hi = (b + 1 == nb) ? pts.size() : std::min(pts.size(), lo + per);
      if (lo >= pts.size()) { break; }
      Block blk;
      blk.rt_low = pts[lo].rt;
      blk.rt_high = pts[hi - 1].rt;
      blk.anchors = hi - lo;
      if (blk.anchors >= opt.min_per_block)
      {
        std::vector<double> x, y, w;
        x.reserve(blk.anchors); y.reserve(blk.anchors); w.reserve(blk.anchors);
        for (std::size_t i = lo; i < hi; ++i)
        { x.push_back(pts[i].l); y.push_back(pts[i].ppm); w.push_back(pts[i].w); }
        const auto f = fitLine(x, y, w, blk.anchors >= opt.min_for_slope);
        blk.intercept_ppm = std::clamp(f.first, -opt.max_correction_ppm, opt.max_correction_ppm);
        blk.slope_ppm = std::clamp(f.second, -opt.max_slope_ppm, opt.max_slope_ppm);
        blk.fitted = true;
      }
      out.blocks_.push_back(blk);
    }
    // The first and last blocks own everything outside the sampled span, so a
    // precursor eluting before the first anchor is corrected rather than not.
    if (!out.blocks_.empty())
    {
      out.blocks_.front().rt_low = -1e30;
      out.blocks_.back().rt_high = 1e30;
    }
    out.fitted_ = !out.blocks_.empty();
    return out;
  }
} // namespace ODIA
