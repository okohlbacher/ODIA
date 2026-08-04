// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Adopted from okohlbacher/OpenDIAlyzer, tag odia-v0.3.0 (0ed5fb2),
// src/opendialyzer.cpp, where these were private statics of the analyzer
// class. Same author, same BSD-3 licence. Changed here: lifted to free
// functions in ODIA::Calibration behind the fit()/identity() entry points.
// The algorithms and their comments are untouched.

#include <odia/RtCalibration.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <numeric>

using OpenMS::Param;
using OpenMS::TransformationDescription;

namespace ODIA::Calibration
{
  namespace
  {
  // --- recalibration fit (shared by the engine and the selftest) --------------
  // Pool-Adjacent-Violators: the optimal L2 monotone-nondecreasing fit to y (weighted by w).
  // Replaces the old forward-max clamp, which pinned y to a running max and so created flat
  // plateaus wherever one noisy bin dipped (distorting the map on both sides). PAVA instead
  // *averages* adjacent violators, the isotonic regression DIA-NN/Calib-RT use. O(n).
  std::vector<double> pava_(const std::vector<double>& y, const std::vector<double>& w)
  {
    assert(y.size() == w.size());                             // parallel; positive finite weights
    const size_t n = y.size();
    std::vector<double> val(n), wt(n);
    std::vector<size_t> len(n);
    size_t m = 0;                                             // active pooled blocks
    for (size_t i = 0; i < n; ++i)
    {
      val[m] = y[i]; wt[m] = w[i]; len[m] = 1;
      while (m > 0 && val[m] < val[m - 1])                    // downward violation -> pool
      {
        const double nw = wt[m - 1] + wt[m];
        val[m - 1] = (val[m - 1] * wt[m - 1] + val[m] * wt[m]) / nw;
        wt[m - 1] = nw; len[m - 1] += len[m];
        --m;
      }
      ++m;
    }
    std::vector<double> out; out.reserve(n);
    for (size_t b = 0; b < m; ++b) { for (size_t k = 0; k < len[b]; ++k) { out.push_back(val[b]); } }
    return out;
  }

  // LOESS: locally weighted linear regression with tricube weights, evaluated at xout.
  //
  // Why this and not the binned median it replaces: binning quantises the anchor set into <=60
  // fixed groups and takes a median per group, so the fit resolution is capped by bin count and a
  // bin straddling a curvature change averages across it. LOESS instead fits a LOCAL line at every
  // evaluation point, weighting neighbours by distance, so it follows gradient curvature (which is
  // exactly what RT warping between a predicted library and a real run looks like) without
  // imposing a global polynomial. It is what DIA-NN and OpenSWATH's own `lowess` alignment use.
  //
  // Output is NOT guaranteed monotone -- local fits can cross -- so callers must still run PAVA.
  // Elution order is physics: a non-monotone RT map is always wrong, however well it fits.
  std::vector<double> loess_(const std::vector<double>& x, const std::vector<double>& y,
                                    const std::vector<double>& w, const std::vector<double>& xout,
                                    double span, int robust_iters = 2)
  {
    const size_t n = x.size();
    std::vector<double> out(xout.size(), 0.0);
    if (n == 0) { return out; }
    if (n < 3)
    {
      for (size_t k = 0; k < xout.size(); ++k) { out[k] = y[n / 2]; }
      return out;
    }
    // Span rule. A FIXED fraction is wrong across the range of anchor counts this sees (19..2000):
    // at small n, 0.3*n is a handful of points and the local line is mostly variance. Require an
    // absolute floor of ~15 neighbours as well, so the effective span GROWS as n shrinks.
    const double eff_span = std::max(span, 15.0 / static_cast<double>(n));
    size_t q = static_cast<size_t>(std::llround(eff_span * static_cast<double>(n)));
    q = std::max<size_t>(std::min<size_t>(15, n), std::min(q, n));
    q = std::max<size_t>(3, q);

    // Robustness weights (Cleveland's bisquare loop). WITHOUT this it is not LOESS, just one pass
    // of local regression -- and that matters here specifically: the anchors are pass-1 peak groups
    // hunted in windows narrower than the RT error, so gross outliers (residuals of thousands of
    // seconds) are guaranteed in the set. A single tricube pass has breakdown ~1/q, i.e. ~2%,
    // whereas the binned-MEDIAN path it replaces has 50%. Skipping the robust loop would swap a
    // high-breakdown estimator for a low-breakdown one in exactly the dirty-anchor regime.
    std::vector<double> rob(n, 1.0);
    std::vector<double> fit_at_x(n, 0.0);
    for (int it = 0; it <= robust_iters; ++it)
    {
      // Evaluate at xout on the final pass, at the anchor x on earlier passes (to get residuals).
      const std::vector<double>& grid = (it == robust_iters) ? xout : x;
      std::vector<double> res(grid.size(), 0.0);
      for (size_t k = 0; k < grid.size(); ++k)
      {
        const double x0 = grid[k];
        size_t lo = static_cast<size_t>(std::lower_bound(x.begin(), x.end(), x0) - x.begin());
        size_t hi = lo;
        if (lo > 0) { --lo; } else { hi = std::min(n, q); }
        while (hi - lo < q && (lo > 0 || hi < n))
        {
          const bool take_left = (hi >= n) || (lo > 0 && (x0 - x[lo - 1]) <= (x[hi] - x0));
          if (take_left && lo > 0) { --lo; } else if (hi < n) { ++hi; } else { break; }
        }
        const double dmax = std::max(std::abs(x0 - x[lo]), std::abs(x[hi - 1] - x0));
        double sw = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (size_t i = lo; i < hi; ++i)
        {
          double u = (dmax > 0.0) ? std::abs(x[i] - x0) / dmax : 0.0;
          if (u >= 1.0) { u = 1.0; }
          const double tri = std::pow(1.0 - u * u * u, 3.0);
          const double wi = tri * rob[i] * (w.empty() ? 1.0 : w[i]);
          if (wi <= 0.0) { continue; }
          const double dx = x[i] - x0;                     // centre at x0 -> intercept IS the fit
          sw += wi; sx += wi * dx; sy += wi * y[i];
          sxx += wi * dx * dx; sxy += wi * dx * y[i];
        }
        if (sw <= 0.0) { res[k] = y[std::min(lo, n - 1)]; continue; }
        const double det = sw * sxx - sx * sx;
        // RELATIVE degeneracy test: det has units of weight^2 * x^4, so an absolute epsilon is
        // scale-dependent and would be miscalibrated by ~1e8 if RT were normalised to [0,1].
        if (std::abs(det) <= 1e-12 * std::abs(sw * sxx)) { res[k] = sy / sw; continue; }
        const double b = (sw * sxy - sx * sy) / det;
        res[k] = (sy - b * sx) / sw;
      }
      if (it == robust_iters) { out = res; break; }
      fit_at_x = res;
      // Bisquare reweight from the median absolute residual.
      std::vector<double> ar(n);
      for (size_t i = 0; i < n; ++i) { ar[i] = std::abs(y[i] - fit_at_x[i]); }
      std::vector<double> tmp = ar;
      std::nth_element(tmp.begin(), tmp.begin() + tmp.size() / 2, tmp.end());
      const double mad = tmp[tmp.size() / 2];
      const double s = 6.0 * mad;
      for (size_t i = 0; i < n; ++i)
      {
        if (s <= 0.0) { rob[i] = 1.0; continue; }          // perfect fit -> no downweighting
        const double u = ar[i] / s;
        rob[i] = (u >= 1.0) ? 0.0 : std::pow(1.0 - u * u, 2.0);
      }
    }
    return out;
  }

  // Write the nearest-rank 95th-pct |anchor residual| of a fitted map (diagnostic).
  void setP95Resid_(const TransformationDescription& td,
                           const std::vector<std::pair<double, double>>& pts, double* p95_resid)
  {
    if (!p95_resid) { return; }
    if (pts.empty()) { *p95_resid = 0.0; return; }
    std::vector<double> r; r.reserve(pts.size());
    for (const auto& pr : pts) { r.push_back(std::abs(pr.second - td.apply(pr.first))); }
    std::sort(r.begin(), r.end());
    const size_t nr = r.size();
    size_t k = (nr * 95 + 99) / 100;                         // ceil(0.95*nr) ...
    k = (k == 0 ? 0 : k - 1);                                // ... nearest-rank index
    *p95_resid = r[std::min(k, nr - 1)];
  }

  TransformationDescription identityTrafo_()
  {
    TransformationDescription td; td.fitModel("identity"); return td;
  }

  // Fit a robust, monotone library_RT -> observed_RT map. Binned medians (outlier-robust) +
  // PAVA isotonic (monotone) + piecewise-linear interpolation. If p95_resid != nullptr, the
  // 95th-pct |anchor residual| after the fit is written there (diagnostic only; this is the
  // IN-SAMPLE residual of the ANCHORS, not the predictive residual on unseen peptides).
  //
  // Robustness (codex review): OpenMS's interpolated model throws below 3 UNIQUE x-values even
  // for linear interpolation (TransformationModelInterpolated::preprocessDataPoints_). Small or
  // degenerate anchor sets are real (few IDs on a small library), so we guarantee >=3 unique x
  // or fall back to identity -- never let fitModel throw.
  // loess_span <= 0 selects the historical binned-median path; >0 fits LOESS at the bin centres
  // and then isotonises. Both end in PAVA + an interpolated model, so the only thing that changes
  // is how the control-point y values are estimated.
  TransformationDescription fitTrafo_(std::vector<std::pair<double, double>> pts,
                                             double* p95_resid = nullptr,
                                             double loess_span = 0.0)
  {
    if (pts.empty()) { if (p95_resid) { *p95_resid = 0.0; } return identityTrafo_(); }
    std::sort(pts.begin(), pts.end());
    const size_t n = pts.size();
    const int NB = std::max(1, std::min<int>(60, static_cast<int>(n / 10)));  // >=~10 anchors/bin (C10)
    std::vector<double> bx, by, bw;                           // per-bin median x, median y, count
    for (int b = 0; b < NB; ++b)
    {
      const size_t lo = n * b / NB, hi = n * (b + 1) / NB;
      if (hi - lo < 3) { continue; }                          // need >=3 for a robust median
      std::vector<double> xs, ys;
      xs.reserve(hi - lo); ys.reserve(hi - lo);
      for (size_t i = lo; i < hi; ++i) { xs.push_back(pts[i].first); ys.push_back(pts[i].second); }
      std::nth_element(xs.begin(), xs.begin() + xs.size() / 2, xs.end());
      std::nth_element(ys.begin(), ys.begin() + ys.size() / 2, ys.end());
      bx.push_back(xs[xs.size() / 2]);
      by.push_back(ys[ys.size() / 2]);
      bw.push_back(static_cast<double>(hi - lo));
    }
    // Coalesce tied bin-x (weighted) BEFORE isotonic so no bin is silently dropped, and the
    // control-point x-values are strictly increasing for the interpolated model.
    std::vector<double> ux, uy, uw;
    for (size_t i = 0; i < bx.size(); ++i)
    {
      if (!ux.empty() && bx[i] <= ux.back() + 1e-6)
      {
        const double nw = uw.back() + bw[i];
        uy.back() = (uy.back() * uw.back() + by[i] * bw[i]) / nw;   // weighted-mean y for tied x
        uw.back() = nw;
      }
      else { ux.push_back(bx[i]); uy.push_back(by[i]); uw.push_back(bw[i]); }
    }
    // LOESS path. Two things had to change from the first attempt:
    //
    // 1. GATING. Binning gives NB = min(60, n/10) bins, so at n=19 -- the low end actually observed
    //    on this data -- NB is 1, ux has ONE element, and the old `ux.size() >= 2` guard skipped
    //    LOESS silently. The feature was inert at exactly the anchor count that motivated it. Below
    //    ~50 anchors a local linear fit is mostly variance anyway, so refuse it EXPLICITLY and say
    //    so, rather than appearing to be on while doing nothing.
    // 2. GRID. Evaluating only at the <=60 bin centres keeps the old shape resolution: LOESS then
    //    improves the y at each knot but never the knot DENSITY, and curvature between knots stays
    //    linearly interpolated. Use a uniform grid over the anchor x-range instead (never beyond
    //    it -- extrapolation is where local fits are worst). Cost is trivial.
    const size_t n_anchor = pts.size();
    if (loess_span > 0.0 && n_anchor >= 50)
    {
      std::vector<double> rx, ry, rw;
      rx.reserve(n_anchor); ry.reserve(n_anchor);
      for (const auto& pr : pts) { rx.push_back(pr.first); ry.push_back(pr.second); }
      const double x_lo = rx.front(), x_hi = rx.back();
      if (x_hi > x_lo)
      {
        const size_t NG = std::min<size_t>(128, std::max<size_t>(8, n_anchor / 2));
        std::vector<double> gx(NG);
        for (size_t i = 0; i < NG; ++i)
        {
          gx[i] = x_lo + (x_hi - x_lo) * static_cast<double>(i) / static_cast<double>(NG - 1);
        }
        const std::vector<double> gy = loess_(rx, ry, rw, gx, loess_span);
        // Replace the binned control points wholesale; weights become uniform because they are no
        // longer bin counts (carrying the stale counts into PAVA would weight a LOESS value by how
        // many anchors happened to fall in an unrelated bin).
        ux = gx; uy = gy; uw.assign(NG, 1.0);
      }
    }
    else if (loess_span > 0.0)
    {
      OPENMS_LOG_INFO << "OpenDIAlyzer: only " << n_anchor << " RT anchors -- using the binned-median "
                      << "fit rather than LOESS (local regression needs >=50 to beat a median)."
                      << std::endl;
    }
    const std::vector<double> iso = (uy.size() >= 2) ? pava_(uy, uw) : uy;   // monotone fit
    TransformationDescription::DataPoints cps;
    for (size_t i = 0; i < ux.size(); ++i) { cps.emplace_back(ux[i], iso[i]); }

    // If binning collapsed to <2 unique x, rebuild from RAW anchors grouped by x -- per-group
    // MEDIAN y (outlier-robust) + PAVA, consistent with the main path (NOT first-y + a forward
    // clamp, which would reintroduce the plateau artefact PAVA exists to avoid). pts is sorted
    // by (x,y), so equal-x anchors are contiguous.
    if (cps.size() < 2)
    {
      std::vector<double> gx, gy, gw;
      size_t i = 0;
      while (i < pts.size())
      {
        size_t j = i;
        while (j < pts.size() && pts[j].first <= pts[i].first + 1e-6) { ++j; }
        std::vector<double> ys; ys.reserve(j - i);
        for (size_t k = i; k < j; ++k) { ys.push_back(pts[k].second); }
        std::nth_element(ys.begin(), ys.begin() + ys.size() / 2, ys.end());
        gx.push_back(pts[i].first);
        gy.push_back(ys[ys.size() / 2]);
        gw.push_back(static_cast<double>(j - i));
        i = j;
      }
      const std::vector<double> giso = (gy.size() >= 2) ? pava_(gy, gw) : gy;
      cps.clear();
      for (size_t k = 0; k < gx.size(); ++k) { cps.emplace_back(gx[k], giso[k]); }
    }
    if (cps.size() < 2)                                       // one unique x -> no map possible
    {
      TransformationDescription td = identityTrafo_();
      setP95Resid_(td, pts, p95_resid);
      return td;
    }
    if (cps.size() == 2)                                      // interpolated needs 3: add collinear midpoint
    {
      const double mx = cps[0].first + 0.5 * (cps[1].first - cps[0].first);   // overflow-safe
      const double my = cps[0].second + 0.5 * (cps[1].second - cps[0].second);
      if (!(mx > cps[0].first && mx < cps[1].first))          // FP degenerate -> can't interpolate
      {
        TransformationDescription td = identityTrafo_();
        setP95Resid_(td, pts, p95_resid);
        return td;
      }
      cps.insert(cps.begin() + 1, std::make_pair(mx, my));
    }
    TransformationDescription td;
    td.setDataPoints(cps);
    Param p;
    p.setValue("interpolation_type", "linear");
    p.setValue("extrapolation_type", "four-point-linear");
    td.fitModel("interpolated", p);
    setP95Resid_(td, pts, p95_resid);
    return td;
  }

  } // namespace

  TransformationDescription fit(std::vector<std::pair<double, double>> anchors,
                                double* p95_resid, double loess_span)
  {
    return fitTrafo_(std::move(anchors), p95_resid, loess_span);
  }

  TransformationDescription identity() { return identityTrafo_(); }

} // namespace ODIA::Calibration
