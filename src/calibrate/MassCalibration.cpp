// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/MassCalibration.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace ODIA
{

  // ==========================================================================
  // Robust location and scale.
  //
  // Ported unchanged in substance from OpenDIAlyzer odia-v0.3.0. The comments
  // are kept because they are the argument for the choice, and the choice is
  // not obvious: a mean would be wrong here for a reason (a narrow true peak on
  // a broad near-uniform background), and so would the shorth.
  // ==========================================================================

  double MassCalibration::halfSampleMode(const std::vector<double>& v)
  {
    // Recursively keep the SHORTEST interval containing half the remaining
    // points; the limit is the densest region, i.e. the mode. Chosen over a KDE
    // mode-finder because it has no bandwidth to select -- and Silverman's
    // critical-bandwidth calibration is not correct for general modality work,
    // so a bandwidth rule would itself need defending -- and because the dense
    // core it converges to is also where the scale has to be measured.
    const std::size_t n = v.size();
    if (n == 0) { return 0.0; }
    if (n == 1) { return v[0]; }
    if (n == 2) { return 0.5 * (v[0] + v[1]); }
    if (n == 3)
    {
      const double d1 = v[1] - v[0], d2 = v[2] - v[1];
      if (d1 < d2) { return 0.5 * (v[0] + v[1]); }
      if (d2 < d1) { return 0.5 * (v[1] + v[2]); }
      return v[1];
    }
    const std::size_t k = (n + 1) / 2;                 // ceil(n/2)
    std::size_t best = 0;
    double best_w = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i + k - 1 < n; ++i)
    {
      const double w = v[i + k - 1] - v[i];
      if (w < best_w) { best_w = w; best = i; }
    }
    return halfSampleMode(std::vector<double>(v.begin() + best, v.begin() + best + k));
  }

  double MassCalibration::localScaleAboutMode(const std::vector<double>& v, double mode,
                                              double init_window)
  {
    // The obvious choice -- the shorth, /1.349 -- is wrong here, and the
    // reference's selftest caught it: the shorth has 50% breakdown, so once the
    // signal is a MINORITY the shortest half is forced to include background and
    // the estimate blows up (measured 11.6 ppm for a 1 ppm signal at 40%
    // purity). The mode survives that because halfSampleMode recurses into the
    // dense core; the scale must be estimated the same way -- locally. Inside a
    // tight neighbourhood the signal is the majority even when it is a global
    // minority, which is exactly the condition MAD needs.
    if (v.size() < 8 || !(init_window > 0.0)) { return 0.0; }
    double w = init_window;
    double sigma = 0.0;
    for (int pass = 0; pass < 3; ++pass)
    {
      std::vector<double> dev;
      dev.reserve(v.size());
      for (double x : v) { if (std::abs(x - mode) <= w) { dev.push_back(std::abs(x - mode)); } }
      if (dev.size() < 8) { break; }                    // neighbourhood too thin to trust
      std::nth_element(dev.begin(), dev.begin() + dev.size() / 2, dev.end());
      const double mad = dev[dev.size() / 2];
      if (!(mad > 0.0)) { break; }
      sigma = 1.4826 * mad;                             // MAD -> Gaussian-equivalent sigma
      const double next = 3.0 * sigma;
      if (next >= w) { break; }                         // converged / not shrinking
      w = next;
    }
    return sigma;
  }

  double MassCalibration::peakednessRatio(const std::vector<double>& abs_err, double window)
  {
    // The guard that stops the estimator fitting noise. Most library targets are
    // absent from any given cell, so a wide search returns an unrelated centroid
    // and those errors are ~UNIFORM over the window. A uniform sample still has
    // a perfectly well-defined mode and scale -- they are just meaningless, and
    // acting on them produces a window that grows until it swallows everything.
    // Real measurement error is centrally peaked; noise is not.
    if (window <= 0.0 || abs_err.empty()) { return 0.0; }
    std::size_t c = 0, e = 0;
    for (double a : abs_err)
    {
      if (a <= 0.2 * window) { ++c; }
      else if (a >= 0.6 * window && a <= 0.8 * window) { ++e; }
    }
    // Both bands are 0.2*window wide, so raw counts are already densities.
    if (e == 0) { return c > 0 ? 1e9 : 0.0; }          // no edge mass at all
    return static_cast<double>(c) / static_cast<double>(e);
  }

  double MassCalibration::backgroundCorrectedScale(const std::vector<double>& dev, double window)
  {
    if (window <= 0.0 || dev.size() < 8) { return 0.0; }
    // Background density, per ppm of |deviation|, from the same 0.6W..0.8W band
    // the peakedness gate reads. Anything that far out is background by
    // assumption -- which is the assumption the gate has just tested.
    std::size_t edge = 0, inside = 0;
    for (double a : dev)
    {
      if (a <= window) { ++inside; }
      if (a >= 0.6 * window && a <= 0.8 * window) { ++edge; }
    }
    const double density = static_cast<double>(edge) / (0.2 * window);
    // Only what is INSIDE the window counts. The sample can extend past it --
    // `search_ppm` is deliberately wider than `gate_ppm` -- and counting those
    // points as excess makes every one of them look like signal. It inflated a
    // planted 2.0 ppm sigma to 5.1 ppm, i.e. an extraction window 2.5x too wide,
    // which is the failure this whole estimator exists to avoid.
    const double excess_total = static_cast<double>(inside) - density * window;
    if (!(excess_total > 0.0)) { return 0.0; }        // nothing above background

    // Walk inwards-out and stop where half the excess has accumulated. That
    // point is the signal's median absolute deviation.
    const double half = 0.5 * excess_total;
    for (std::size_t i = 0; i < dev.size(); ++i)
    {
      const double h = dev[i];
      if (h > window) { break; }
      const double excess = static_cast<double>(i + 1) - density * h;
      if (excess >= half) { return 1.4826 * h; }
    }
    return 0.0;
  }

  double MassCalibration::refineLocation(const std::vector<double>& sorted, double start,
                                         double window)
  {
    double centre = start;
    for (int pass = 0; pass < 6; ++pass)
    {
      std::vector<double> dev;
      dev.reserve(sorted.size());
      for (double x : sorted) { dev.push_back(std::abs(x - centre)); }
      std::sort(dev.begin(), dev.end());
      const double sigma = backgroundCorrectedScale(dev, window);
      if (!(sigma > 0.0)) { break; }
      std::vector<double> near;
      for (double x : sorted) { if (std::abs(x - centre) <= 2.0 * sigma) { near.push_back(x); } }
      if (near.size() < 8) { break; }
      std::nth_element(near.begin(), near.begin() + near.size() / 2, near.end());
      const double next = near[near.size() / 2];
      const bool settled = std::abs(next - centre) < 0.01;
      centre = next;
      if (settled) { break; }
    }
    return centre;
  }

  namespace
  {
    /// Median of a copy, so callers keep their ordering.
    double medianOf(std::vector<double> v)
    {
      if (v.empty()) { return 0.0; }
      std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
      return v[v.size() / 2];
    }

    double madSigma(const std::vector<double>& v, double centre)
    {
      if (v.size() < 4) { return 0.0; }
      std::vector<double> dev;
      dev.reserve(v.size());
      for (double x : v) { dev.push_back(std::abs(x - centre)); }
      std::nth_element(dev.begin(), dev.begin() + dev.size() / 2, dev.end());
      return 1.4826 * dev[dev.size() / 2];
    }

    /// A weighted straight-line fit, returning the slope's standard error too.
    ///
    /// The standard error is the point. Without it there is no way to tell a
    /// slope the data supports from one it merely permits, and the whole model
    /// choice below turns on that distinction.
    struct Line
    {
      double intercept = 0.0, slope = 0.0, slope_stderr = 0.0;
      bool ok = false;
    };

    Line weightedLine(const std::vector<double>& x, const std::vector<double>& y,
                      const std::vector<double>& sigma)
    {
      Line out;
      if (x.size() < 3) { return out; }
      double sw = 0, swx = 0, swy = 0, swxx = 0, swxy = 0;
      for (std::size_t i = 0; i < x.size(); ++i)
      {
        const double s = sigma[i] > 0.0 ? sigma[i] : 1.0;
        const double w = 1.0 / (s * s);
        sw += w; swx += w * x[i]; swy += w * y[i];
        swxx += w * x[i] * x[i]; swxy += w * x[i] * y[i];
      }
      const double det = sw * swxx - swx * swx;
      if (!(std::abs(det) > 0.0)) { return out; }
      out.slope = (sw * swxy - swx * swy) / det;
      out.intercept = (swxx * swy - swx * swxy) / det;
      // Scaled by the fit's own chi-square per degree of freedom: the per-bin
      // errors are themselves estimates, so trusting them unscaled would call a
      // slope significant whenever the bins happened to be optimistic.
      double chi2 = 0.0;
      for (std::size_t i = 0; i < x.size(); ++i)
      {
        const double s = sigma[i] > 0.0 ? sigma[i] : 1.0;
        const double r = (y[i] - (out.intercept + out.slope * x[i])) / s;
        chi2 += r * r;
      }
      const double dof = static_cast<double>(x.size()) - 2.0;
      const double scale = dof > 0.0 ? std::max(1.0, chi2 / dof) : 1.0;
      out.slope_stderr = std::sqrt(scale * sw / det);
      out.ok = std::isfinite(out.slope) && std::isfinite(out.intercept);
      return out;
    }

    /// Split residuals into equal-count bins along @p key and summarise each.
    ///
    /// Equal COUNT, not equal width: the m/z distribution of a tryptic library
    /// is very far from uniform, and equal-width bins put a handful of points in
    /// the tails, which is precisely where a slope fitted to them would do the
    /// most damage.
    std::vector<MassCalibration::Bin> binBy(std::vector<std::pair<double, double>> kv,
                                            std::size_t max_bins, std::size_t min_per_bin)
    {
      std::vector<MassCalibration::Bin> bins;
      if (kv.size() < min_per_bin) { return bins; }
      std::sort(kv.begin(), kv.end(),
                [](const auto& a, const auto& b) { return a.first < b.first; });
      std::size_t nb = std::min(max_bins, kv.size() / min_per_bin);
      if (nb == 0) { return bins; }
      for (std::size_t b = 0; b < nb; ++b)
      {
        const std::size_t lo = b * kv.size() / nb;
        const std::size_t hi = (b + 1) * kv.size() / nb;
        if (hi <= lo) { continue; }
        std::vector<double> vals;
        vals.reserve(hi - lo);
        double key_sum = 0.0;
        for (std::size_t i = lo; i < hi; ++i) { vals.push_back(kv[i].second); key_sum += kv[i].first; }
        MassCalibration::Bin bin;
        bin.n = hi - lo;
        bin.centre = key_sum / static_cast<double>(bin.n);
        bin.location = medianOf(vals);
        const double sd = madSigma(vals, bin.location);
        // 1.2533 = sqrt(pi/2), the median's standard error against the mean's.
        bin.stderr_ppm = sd > 0.0 ? 1.2533 * sd / std::sqrt(static_cast<double>(bin.n)) : 0.0;
        bins.push_back(bin);
      }
      return bins;
    }
  } // namespace

  // ==========================================================================
  // The fit.
  // ==========================================================================

  MassCalibration::Model MassCalibration::fit(const std::vector<MassResidual>& residuals,
                                              const Options& opt, Diagnostics* diag)
  {
    Model m;
    m.gate_peakedness = opt.min_peakedness;
    m.sigma_multiple = opt.sigma_multiple;
    m.floor_ppm = opt.floor_ppm;
    m.search_ppm = opt.search_ppm;
    const double gate_w = opt.gate_ppm > 0.0 ? std::min(opt.gate_ppm, opt.search_ppm)
                                             : opt.search_ppm;
    m.gate_ppm = gate_w;

    // Brightness selection, applied to the CONTROL on the same threshold so the
    // null keeps measuring the same procedure.
    double threshold = 0.0;
    if (opt.min_intensity_quantile > 0.0)
    {
      std::vector<double> intensity;
      intensity.reserve(residuals.size());
      for (const auto& r : residuals)
      {
        if (!r.decoy && std::isfinite(r.ppm)) { intensity.push_back(r.intensity); }
      }
      if (!intensity.empty())
      {
        const std::size_t k = static_cast<std::size_t>(
          std::min(opt.min_intensity_quantile, 1.0) * static_cast<double>(intensity.size() - 1));
        std::nth_element(intensity.begin(), intensity.begin() + k, intensity.end());
        threshold = intensity[k];
      }
    }

    std::vector<MassResidual> kept;
    kept.reserve(residuals.size());
    std::vector<double> target, decoy;
    std::vector<std::pair<double, double>> by_mz;   // (m/z, ppm)
    m.mz_low = std::numeric_limits<double>::max();
    m.mz_high = 0.0;
    for (const auto& r : residuals)
    {
      if (!std::isfinite(r.ppm)) { continue; }
      if (r.intensity < threshold) { continue; }
      kept.push_back(r);
      if (r.decoy) { decoy.push_back(r.ppm); continue; }
      target.push_back(r.ppm);
      by_mz.emplace_back(r.mz, r.ppm);
      m.mz_low = std::min<double>(m.mz_low, r.mz);
      m.mz_high = std::max<double>(m.mz_high, r.mz);
    }
    if (target.empty()) { m.mz_low = 0.0; }
    m.residuals = target.size();
    m.decoy_residuals = decoy.size();

    if (static_cast<int>(target.size()) < opt.min_residuals)
    {
      m.reason = "only " + std::to_string(target.size()) + " residuals (need " +
                 std::to_string(opt.min_residuals) + ")";
      return m;
    }

    std::vector<double> sorted = target;
    std::sort(sorted.begin(), sorted.end());
    const double mode0 = refineLocation(sorted, halfSampleMode(sorted), gate_w);

    std::vector<double> dev;
    dev.reserve(sorted.size());
    for (double e : sorted) { dev.push_back(std::abs(e - mode0)); }
    m.peakedness = peakednessRatio(dev, gate_w);
    std::sort(dev.begin(), dev.end());
    m.sigma_before = backgroundCorrectedScale(dev, gate_w);

    if (!decoy.empty())
    {
      std::vector<double> ds = decoy;
      std::sort(ds.begin(), ds.end());
      const double dm = refineLocation(ds, halfSampleMode(ds), gate_w);
      std::vector<double> dd;
      dd.reserve(ds.size());
      for (double e : ds) { dd.push_back(std::abs(e - dm)); }
      m.decoy_peakedness = peakednessRatio(dd, gate_w);
    }

    // ---- the gate ---------------------------------------------------------
    // A flat residual distribution is not a calibration waiting to be found; it
    // is the absence of one, and applying a mode fitted to it is strictly worse
    // than doing nothing.
    if (m.peakedness < opt.min_peakedness)
    {
      char buf[256];
      std::snprintf(buf, sizeof buf,
                    "residuals are FLAT (peakedness %.2f < %.2f over %zu residuals) -- "
                    "that is what a mostly-noise sample looks like; a mode and a width can "
                    "still be computed from it and would be meaningless",
                    m.peakedness, opt.min_peakedness, target.size());
      m.reason = buf;
      return m;
    }
    if (!(m.sigma_before > 0.0) || !std::isfinite(m.sigma_before) || !std::isfinite(mode0))
    {
      m.reason = "degenerate scale";
      return m;
    }
    // The control cells pass the same mobility, co-occurrence and apex tests and
    // differ only in that they cannot hold the real ions. If they are as peaked
    // as the data, the peak is coming from the procedure, not the instrument.
    if (!decoy.empty() && m.decoy_peakedness >= m.peakedness)
    {
      char buf[256];
      std::snprintf(buf, sizeof buf,
                    "the m/z-shifted control is as peaked as the data (%.2f vs %.2f) -- "
                    "whatever structure is there is not the fragments",
                    m.decoy_peakedness, m.peakedness);
      m.reason = buf;
      return m;
    }

    // ---- shape: is the error a function of m/z? ---------------------------
    // Restricted to the signal core. Outside it the sample is background, whose
    // location is the middle of the search window at every m/z -- including it
    // would flatten any real trend towards zero.
    std::vector<std::pair<double, double>> core;
    core.reserve(by_mz.size());
    for (const auto& kv : by_mz)
    {
      if (std::abs(kv.second - mode0) <= 4.0 * m.sigma_before) { core.push_back(kv); }
    }
    auto bins = binBy(core, 8, 40);
    if (diag) { diag->by_mz = bins; }

    m.reference_mz = 700.0;
    if (!bins.empty())
    {
      double s = 0.0;
      for (const auto& b : bins) { s += b.centre; }
      m.reference_mz = s / static_cast<double>(bins.size());
    }

    Line line;
    if (bins.size() >= 4)
    {
      std::vector<double> x, y, sg;
      for (const auto& b : bins)
      {
        x.push_back(std::log(b.centre / m.reference_mz));
        y.push_back(b.location);
        sg.push_back(b.stderr_ppm);
      }
      line = weightedLine(x, y, sg);
    }
    m.slope_t = (line.ok && line.slope_stderr > 0.0)
                  ? std::abs(line.slope) / line.slope_stderr : 0.0;
    m.shape_swing_ppm = line.ok && m.mz_low > 0.0
                          ? std::abs(line.slope) * std::log(m.mz_high / m.mz_low) : 0.0;

    // Both candidates are evaluated, and the CHOICE is made on the residual each
    // one leaves. `applied` centres the model on the global mode of its own
    // corrected residuals, so the two are compared on equal terms rather than
    // one of them carrying an accidental offset into its spread.
    const auto sigmaAfter = [&](double intercept, double slope, double* centred_intercept) {
      std::vector<double> c;
      c.reserve(target.size());
      for (const auto& r : kept)
      {
        if (r.decoy) { continue; }
        const double corr = (slope != 0.0 && r.mz > 0.0)
                              ? intercept + slope * std::log(r.mz / m.reference_mz)
                              : intercept;
        c.push_back(r.ppm - corr);
      }
      std::sort(c.begin(), c.end());
      // refineLocation, not the bare mode. The mode's flank bias is not a
      // starting-value problem that goes away once a shape is subtracted -- it
      // is a property of this distribution -- so the final centre needs the same
      // treatment the initial one got. Skipping it here left the S08 fit 2.3 ppm
      // high (-7.64 instead of -9.9) with a mass axis still visibly off centre
      // in every retention-time bin, which is the whole defect this class exists
      // to remove.
      const double mode = refineLocation(c, halfSampleMode(c), gate_w);
      if (centred_intercept) { *centred_intercept = intercept + mode; }
      std::vector<double> d;
      d.reserve(c.size());
      for (double x : c) { d.push_back(std::abs(x - mode)); }
      std::sort(d.begin(), d.end());
      return backgroundCorrectedScale(d, gate_w);
    };

    double constant_intercept = mode0;
    m.sigma_constant = sigmaAfter(mode0, 0.0, &constant_intercept);
    double shaped_intercept = line.ok ? line.intercept : 0.0;
    if (line.ok) { m.sigma_shaped = sigmaAfter(line.intercept, line.slope, &shaped_intercept); }

    // MODEL CHOICE. The simpler model wins unless the data insists TWICE: the
    // slope must be several standard errors from zero, and modelling the shape
    // must leave a materially tighter residual than the constant does.
    //
    // Significance alone is not enough, and S08 is the case that proves it. Its
    // m/z trend is real and large -- -12.1 ppm at 200-288 Th rising to -5.8 ppm
    // at 1067-1694 Th over 124 M hits, so a swing of ~6 ppm and a t of many tens
    // -- and correcting it still moves the residual MAD-SD only from 7.82 to
    // 7.62 ppm. The systematic part is small against the per-fragment scatter,
    // so the shape buys 2.6%, and 2.6% does not justify a model that can contort
    // where the library has no fragments to hold it down. A criterion based on
    // the swing, or on significance, would have shipped that curve.
    const bool take_shape = line.ok && bins.size() >= 4 &&
                            m.slope_t >= opt.min_slope_t &&
                            m.sigma_shaped > 0.0 && m.sigma_constant > 0.0 &&
                            m.sigma_shaped <= opt.max_sigma_ratio * m.sigma_constant;

    m.fitted = true;
    char buf[420];
    if (take_shape)
    {
      m.form = "log_mz";
      m.log_slope_ppm = line.slope;
      m.intercept_ppm = shaped_intercept;
      m.sigma_after = m.sigma_shaped;
      std::snprintf(buf, sizeof buf,
                    "log m/z: %.3f ppm per e-fold (t=%.1f), swing %.2f ppm over the %.0f-%.0f Th "
                    "sampled, and it tightens the residual from %.2f to %.2f ppm (ratio %.2f "
                    "<= %.2f)", line.slope, m.slope_t, m.shape_swing_ppm, m.mz_low, m.mz_high,
                    m.sigma_constant, m.sigma_shaped, m.sigma_shaped / m.sigma_constant,
                    opt.max_sigma_ratio);
    }
    else
    {
      m.form = "constant";
      m.log_slope_ppm = 0.0;
      m.intercept_ppm = constant_intercept;
      m.sigma_after = m.sigma_constant;
      if (!line.ok)
      {
        std::snprintf(buf, sizeof buf,
                      "constant: %zu m/z bins is too few to fit a shape through", bins.size());
      }
      else if (m.slope_t < opt.min_slope_t)
      {
        std::snprintf(buf, sizeof buf,
                      "constant: the m/z slope is %.3f ppm per e-fold at only t=%.1f (need %.1f), "
                      "so the residual is flat in m/z and a curve would only contort at the ends",
                      line.slope, m.slope_t, opt.min_slope_t);
      }
      else
      {
        std::snprintf(buf, sizeof buf,
                      "constant: the m/z trend is real (%.3f ppm per e-fold, t=%.1f, swing %.2f "
                      "ppm over %.0f-%.0f Th) but does not pay -- it moves the residual only "
                      "%.2f -> %.2f ppm (ratio %.2f > %.2f), which does not justify a model that "
                      "can contort where the library has no fragments",
                      line.slope, m.slope_t, m.shape_swing_ppm, m.mz_low, m.mz_high,
                      m.sigma_constant, m.sigma_shaped, m.sigma_shaped / m.sigma_constant,
                      opt.max_sigma_ratio);
      }
    }
    m.reason = buf;

    // ---- retention time, as a diagnostic ----------------------------------
    // Reported and never applied. Two reasons, and the first is the decisive
    // one: on S08 there is nothing to apply -- the corrected residual sits
    // between -8.4 and -9.0 ppm across the whole 1,859 s gradient. The second is
    // that a time-dependent correction cannot be applied where the m/z-dependent
    // one is, on the transition at index-build time; it would have to move every
    // PEAK, per spectrum. Paying that for a term this run does not have would be
    // speculative work, so the number is surfaced and the decision is left to the
    // run that actually shows a drift.
    std::vector<std::pair<double, double>> by_rt;
    for (const auto& r : kept)
    {
      if (r.decoy) { continue; }
      const double c = r.ppm - m.ppmAt(r.mz);
      if (std::abs(c) <= 4.0 * std::max(m.sigma_after, 0.5)) { by_rt.emplace_back(r.rt, c); }
    }
    auto rt_bins = binBy(by_rt, 8, 40);
    if (diag) { diag->by_rt = rt_bins; }
    if (rt_bins.size() >= 4)
    {
      std::vector<double> x, y, sg;
      for (const auto& b : rt_bins)
      {
        x.push_back(b.centre / 1000.0);
        y.push_back(b.location);
        sg.push_back(b.stderr_ppm);
      }
      const Line rt_line = weightedLine(x, y, sg);
      if (rt_line.ok)
      {
        const double rt_span = (rt_bins.back().centre - rt_bins.front().centre) / 1000.0;
        m.rt_drift_ppm = rt_line.slope * rt_span;
        m.rt_drift_t = rt_line.slope_stderr > 0.0
                         ? std::abs(rt_line.slope) / rt_line.slope_stderr : 0.0;
      }
    }

    if (!(m.sigma_after > 0.0) || !std::isfinite(m.sigma_after))
    {
      m.window_ppm = -1.0;
      m.reason += "; corrected scale degenerate, so no window is inferred";
      return m;
    }

    m.window_ppm = std::max(opt.sigma_multiple * m.sigma_after, opt.floor_ppm);
    // An estimate can never legitimately exceed the width it searched -- that is
    // a fit to the edge of the search window, i.e. noise. This rail is the
    // reference's, and it is there because an ungated fit once returned 72.6 ppm
    // on an instrument measured at 1.66 ppm and cost 2,302 identifications.
    if (m.window_ppm > gate_w)
    {
      char buf[220];
      std::snprintf(buf, sizeof buf,
                    "; window %.1f ppm is at or beyond the %.1f ppm width the scale was measured "
                    "over -- that is a fit to the window edge, not the instrument, so no window "
                    "is inferred", m.window_ppm, gate_w);
      m.reason += buf;
      m.window_ppm = -1.0;
    }
    return m;
  }

  // ==========================================================================
  // Collection: probing the run for candidate matches.
  // ==========================================================================

  namespace
  {
    /// One query: a fragment m/z (possibly shifted, for a control cell) and the
    /// cell it belongs to.
    struct Query
    {
      double mz;
      std::uint32_t cell;      ///< precursor slot * variants + variant
      std::uint8_t fragment;
    };

    struct WindowIndex
    {
      std::vector<Query> query;                ///< ascending in m/z
      std::vector<std::uint32_t> bucket;       ///< log-m/z bucket -> first query
      /// cell -> its queries, as a dense (local cell, fragment) -> query slot
      /// table. Built after the sort, because the sort is what scrambles it.
      std::vector<std::uint32_t> slot;
      std::vector<std::uint32_t> cells;        ///< global cell ids, in local order
      std::unordered_map<std::uint32_t, std::uint32_t> local_of_cell;
      std::size_t fragments = 0;
      double log_base = 0.0, log_lo = 0.0;

      std::size_t bucketOf(double m) const
      {
        const auto b = static_cast<std::ptrdiff_t>((std::log(m) - log_lo) / log_base);
        return std::size_t(std::clamp<std::ptrdiff_t>(b, 0, std::ptrdiff_t(bucket.size()) - 2));
      }

      void build(double tolerance_ppm)
      {
        if (query.empty()) { return; }
        std::sort(query.begin(), query.end(),
                  [](const Query& a, const Query& b) { return a.mz < b.mz; });
        log_base = std::log1p(tolerance_ppm * 1e-6);
        if (!(log_base > 0.0)) { log_base = 1e-6; }
        log_lo = std::log(query.front().mz);
        const auto n = static_cast<std::size_t>(
                         (std::log(query.back().mz) - log_lo) / log_base) + 2;
        bucket.assign(n + 1, static_cast<std::uint32_t>(query.size()));
        for (std::size_t i = query.size(); i-- > 0;)
        {
          bucket[bucketOf(query[i].mz)] = static_cast<std::uint32_t>(i);
        }
        for (std::size_t b = bucket.size() - 1; b-- > 0;)
        {
          bucket[b] = std::min(bucket[b], bucket[b + 1]);
        }
        slot.assign(cells.size() * fragments, std::numeric_limits<std::uint32_t>::max());
        for (std::size_t i = 0; i < query.size(); ++i)
        {
          const auto it = local_of_cell.find(query[i].cell);
          if (it == local_of_cell.end()) { continue; }
          slot[it->second * fragments + query[i].fragment] = static_cast<std::uint32_t>(i);
        }
      }
    };

    /// The best cell found so far for one (precursor, variant).
    ///
    /// "Best" is the largest summed matched intensity, which with no retention
    /// time map is the only available stand-in for the elution apex -- and it is
    /// the right one: a real precursor's fragments are all brightest together.
    struct Cell
    {
      float rt = 0.0f;
      float total = 0.0f;
      std::uint8_t n = 0;
      std::vector<float> mz, ppm, intensity;
    };
  } // namespace

  std::vector<MassResidual> MassCalibration::collect(const Library& library,
                                                     SpectrumSource& source,
                                                     const Options& opt, Diagnostics* diag)
  {
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<MassResidual> out;

    const auto& p = library.precursors();
    const auto& tr = library.transitions();
    const auto& windows = source.windows();
    const auto& info = source.spectra();
    if (windows.empty() || info.empty() || library.precursorCount() == 0) { return out; }

    const std::size_t variants = 1 + opt.decoy_shifts.size();
    // The (variant, fragment) pair is packed into one byte on each query, so the
    // product is bounded here rather than silently wrapping into another
    // fragment's slot.
    const std::size_t max_frag = std::clamp<std::size_t>(opt.max_fragments, 1, 255 / variants);

    // ---- which window each spectrum belongs to, and where the cycles are ----
    std::vector<std::uint32_t> window_of(info.size(), std::numeric_limits<std::uint32_t>::max());
    for (std::size_t si = 0; si < info.size(); ++si)
    {
      for (std::size_t w = 0; w < windows.size(); ++w)
      {
        if (std::abs(info[si].window.mz_low - windows[w].mz_low) < 1e-6 &&
            std::abs(info[si].window.mz_high - windows[w].mz_high) < 1e-6)
        {
          window_of[si] = static_cast<std::uint32_t>(w);
          break;
        }
      }
    }

    // A DIA cycle is a contiguous run of spectra covering each window once, so a
    // cycle boundary is where a window repeats. Contiguity is what makes the
    // sample cheap: one cycle is ONE range request, and every window is probed
    // at the same retention time.
    std::vector<std::pair<std::size_t, std::size_t>> cycles;
    {
      std::vector<char> seen(windows.size(), 0);
      std::size_t begin = 0;
      for (std::size_t si = 0; si < info.size(); ++si)
      {
        const std::uint32_t w = window_of[si];
        if (w == std::numeric_limits<std::uint32_t>::max()) { continue; }
        if (seen[w])
        {
          cycles.emplace_back(begin, si);
          std::fill(seen.begin(), seen.end(), 0);
          begin = si;
        }
        seen[w] = 1;
      }
      if (begin < info.size()) { cycles.emplace_back(begin, info.size()); }
    }
    if (cycles.empty()) { return out; }

    // ---- sample precursors, spread through the library ---------------------
    std::vector<std::uint32_t> sampled;
    {
      std::vector<std::uint32_t> eligible;
      eligible.reserve(library.precursorCount());
      for (std::size_t i = 0; i < library.precursorCount(); ++i)
      {
        if (p.decoy[i]) { continue; }
        if (p.mz[i] == MZ_INVALID) { continue; }
        if (p.transition_count[i] < opt.min_fragments_matched) { continue; }
        eligible.push_back(static_cast<std::uint32_t>(i));
      }
      if (eligible.empty()) { return out; }
      // A stride rather than the first N: the library may be sorted by m/z, and
      // taking a prefix would then measure the bottom of the mass range only --
      // which is precisely the axis the shape fit has to see.
      const std::size_t want = std::min(opt.max_precursors, eligible.size());
      for (std::size_t k = 0; k < want; ++k)
      {
        sampled.push_back(eligible[k * eligible.size() / want]);
      }
    }

    // ---- assign each to one window, and build that window's query index ----
    std::vector<WindowIndex> index(windows.size());
    for (auto& x : index) { x.fragments = max_frag * variants; }
    std::vector<std::uint32_t> window_of_precursor(sampled.size(),
                                                   std::numeric_limits<std::uint32_t>::max());
    for (std::size_t s = 0; s < sampled.size(); ++s)
    {
      const double mz = fromFixed(p.mz[sampled[s]]);
      // The most CENTRAL window, when several contain it. An edge precursor is
      // transmitted with reduced efficiency, so the central window is the one
      // whose spectra actually hold its fragments.
      std::uint32_t best = std::numeric_limits<std::uint32_t>::max();
      double best_margin = -1.0;
      for (std::size_t w = 0; w < windows.size(); ++w)
      {
        if (!windows[w].contains(mz)) { continue; }
        const double margin = std::min(mz - windows[w].mz_low, windows[w].mz_high - mz);
        if (margin > best_margin) { best_margin = margin; best = static_cast<std::uint32_t>(w); }
      }
      window_of_precursor[s] = best;
    }

    // Cells are numbered globally so a cell id carries its precursor and its
    // variant without a side table.
    const auto cell_id = [variants](std::size_t slot, std::size_t variant) {
      return static_cast<std::uint32_t>(slot * variants + variant);
    };

    for (std::size_t s = 0; s < sampled.size(); ++s)
    {
      const std::uint32_t w = window_of_precursor[s];
      if (w == std::numeric_limits<std::uint32_t>::max()) { continue; }
      auto& x = index[w];
      const std::size_t i = sampled[s];
      const std::size_t nfrag = std::min<std::size_t>(p.transition_count[i], max_frag);
      for (std::size_t v = 0; v < variants; ++v)
      {
        const double shift = v == 0 ? 0.0 : opt.decoy_shifts[v - 1];
        const std::uint32_t cid = cell_id(s, v);
        bool any = false;
        for (std::size_t k = 0; k < nfrag; ++k)
        {
          const double mz = fromFixed(tr.product_mz[p.transition_begin[i] + k]) + shift;
          if (!(mz > 1.0)) { continue; }
          x.query.push_back({mz, cid, static_cast<std::uint8_t>(v * max_frag + k)});
          any = true;
        }
        if (any)
        {
          x.local_of_cell.emplace(cid, static_cast<std::uint32_t>(x.cells.size()));
          x.cells.push_back(cid);
        }
      }
    }
    for (auto& x : index) { x.build(opt.search_ppm); }

    // ---- probe -------------------------------------------------------------
    std::vector<Cell> cell(sampled.size() * variants);
    std::vector<SpectrumPeaks> block;
    std::vector<float> best_int;
    std::vector<double> best_mz;

    const std::size_t want_cycles = std::min(opt.cycles, cycles.size());
    std::size_t decoded = 0, probed = 0;
    for (std::size_t c = 0; c < want_cycles; ++c)
    {
      const auto& cyc = cycles[c * cycles.size() / want_cycles];
      source.peaks(cyc.first, cyc.second, block);
      decoded += cyc.second - cyc.first;

      for (std::size_t si = cyc.first; si < cyc.second; ++si)
      {
        const std::uint32_t w = window_of[si];
        if (w == std::numeric_limits<std::uint32_t>::max()) { continue; }
        auto& x = index[w];
        if (x.query.empty()) { continue; }
        const auto& peaks = block[si - cyc.first];
        if (peaks.empty()) { continue; }

        best_int.assign(x.query.size(), -1.0f);
        best_mz.assign(x.query.size(), 0.0);

        const bool use_im = opt.use_ion_mobility && peaks.hasIonMobility();
        const double im_low = info[si].window.im_low, im_high = info[si].window.im_high;
        const double ppm = opt.search_ppm * 1e-6;
        const double front = x.query.front().mz, back = x.query.back().mz;

        for (std::size_t k = 0; k < peaks.size(); ++k)
        {
          const double m = peaks.mz[k];
          const double slack = m * ppm * 1.01 + 1e-6;
          if (m + slack < front || m - slack > back) { continue; }
          double peak_im = std::numeric_limits<double>::quiet_NaN();
          if (use_im)
          {
            peak_im = peaks.ion_mobility[k];
            if (peak_im < im_low || peak_im > im_high) { continue; }
          }
          const float intensity = peaks.intensity[k];
          std::size_t i = x.bucket[x.bucketOf(std::max(m - slack, front))];
          for (; i < x.query.size() && x.query[i].mz <= m + slack; ++i)
          {
            const double q = x.query[i].mz;
            if (std::abs(m - q) > q * ppm) { continue; }
            if (opt.im_window > 0.0 && !std::isnan(peak_im))
            {
              const float want = p.im[sampled[x.query[i].cell / variants]];
              if (!std::isnan(want) && std::abs(peak_im - want) > opt.im_window) { continue; }
            }
            // MOST INTENSE within the window, not nearest. On this instrument
            // class the interferents are numerous but individually weak, so
            // "nearest" preferentially selects single-ion noise sitting close to
            // the query centre, while "most intense" selects the real peak
            // whenever it is above that floor.
            if (intensity > best_int[i]) { best_int[i] = intensity; best_mz[i] = m; }
          }
        }

        // Co-occurrence, then apex. Both are applied identically to target and
        // control cells, which is what makes the control a null for the whole
        // procedure rather than for the mass window alone.
        for (std::size_t lc = 0; lc < x.cells.size(); ++lc)
        {
          const std::uint32_t cid = x.cells[lc];
          const std::size_t variant = cid % variants;
          std::uint8_t matched = 0;
          float total = 0.0f;
          for (std::size_t f = 0; f < max_frag; ++f)
          {
            const std::uint32_t qi = x.slot[lc * x.fragments + variant * max_frag + f];
            if (qi == std::numeric_limits<std::uint32_t>::max()) { continue; }
            if (best_int[qi] < 0.0f) { continue; }
            ++matched;
            total += best_int[qi];
          }
          ++probed;
          if (matched < opt.min_fragments_matched) { continue; }
          Cell& cl = cell[cid];
          if (total <= cl.total) { continue; }
          cl.total = total;
          cl.rt = static_cast<float>(info[si].retention_time);
          cl.n = matched;
          cl.mz.clear(); cl.ppm.clear(); cl.intensity.clear();
          for (std::size_t f = 0; f < max_frag; ++f)
          {
            const std::uint32_t qi = x.slot[lc * x.fragments + variant * max_frag + f];
            if (qi == std::numeric_limits<std::uint32_t>::max()) { continue; }
            if (best_int[qi] < 0.0f) { continue; }
            const double q = x.query[qi].mz;
            cl.mz.push_back(static_cast<float>(q));
            cl.ppm.push_back(static_cast<float>((best_mz[qi] - q) / q * 1e6));
            cl.intensity.push_back(best_int[qi]);
          }
        }
      }
    }

    std::size_t target_cells = 0, decoy_cells = 0;
    for (std::size_t cid = 0; cid < cell.size(); ++cid)
    {
      const Cell& cl = cell[cid];
      if (cl.n == 0) { continue; }
      const bool is_decoy = (cid % variants) != 0;
      if (is_decoy) { ++decoy_cells; } else { ++target_cells; }
      for (std::size_t k = 0; k < cl.ppm.size(); ++k)
      {
        MassResidual r;
        r.mz = cl.mz[k];
        r.rt = cl.rt;
        r.ppm = cl.ppm[k];
        r.intensity = cl.intensity[k];
        r.decoy = is_decoy;
        out.push_back(r);
      }
    }

    if (diag)
    {
      diag->cells_probed = probed;
      diag->target_cells = target_cells;
      diag->decoy_cells = decoy_cells;
      diag->spectra_decoded = decoded;
      diag->collect_seconds = std::chrono::duration<double>(
                                std::chrono::steady_clock::now() - t0).count();
    }
    return out;
  }

  MassCalibration::Model MassCalibration::calibrate(const Library& library,
                                                    SpectrumSource& source,
                                                    const Options& opt, Diagnostics* diag)
  {
    const auto residuals = collect(library, source, opt, diag);
    return fit(residuals, opt, diag);
  }

  std::string MassCalibration::report(const Model& m, const Diagnostics* diag)
  {
    std::ostringstream o;
    o.setf(std::ios::fixed);
    o.precision(2);
    o << "fragment mass calibration: " << (m.fitted ? "GATE PASSED" : "GATE FAILED")
      << " -- " << m.reason << "\n";
    o << "  residuals:   " << m.residuals << " target, " << m.decoy_residuals << " control";
    if (diag && diag->spectra_decoded)
    {
      o << " (from " << diag->target_cells << " target and " << diag->decoy_cells
        << " control cells over " << diag->spectra_decoded << " spectra, "
        << diag->collect_seconds << " s)";
    }
    o << "\n";
    o << "  peakedness:  " << m.peakedness << " (control " << m.decoy_peakedness
      << "), gate at " << m.gate_peakedness << " evaluated over +/-" << m.gate_ppm
      << " ppm of a +/-" << m.search_ppm << " ppm search\n";
    if (m.fitted)
    {
      o << "  correction:  " << m.form << ", " << m.intercept_ppm << " ppm at "
        << m.reference_mz << " Th";
      if (m.log_slope_ppm != 0.0)
      {
        o << " + " << m.log_slope_ppm << " ppm per e-fold in m/z";
      }
      o << "\n";
      o << "  model choice: constant leaves " << m.sigma_constant << " ppm, shaped leaves "
        << m.sigma_shaped << " ppm (slope t=" << m.slope_t << ", swing " << m.shape_swing_ppm
        << " ppm); chose " << m.form << "\n";
      o << "  sigma:       " << m.sigma_before << " ppm before correction, "
        << m.sigma_after << " ppm after\n";
      o << "  rt drift:    " << m.rt_drift_ppm << " ppm across the sampled gradient (t="
        << m.rt_drift_t << "), reported only\n";
      if (m.window_ppm > 0.0)
      {
        o << "  window:      +/-" << m.window_ppm << " ppm (" << m.sigma_multiple
          << " sigma, floor " << m.floor_ppm << ")\n";
      }
      else
      {
        o << "  window:      not inferred\n";
      }
    }
    if (diag && !diag->by_mz.empty())
    {
      o << "  residual vs m/z:\n";
      for (const auto& b : diag->by_mz)
      {
        o << "    " << b.centre << " Th: " << b.location << " +/- " << b.stderr_ppm
          << " ppm (n=" << b.n << ")\n";
      }
    }
    if (diag && !diag->by_rt.empty())
    {
      o << "  residual vs RT (after correction):\n";
      for (const auto& b : diag->by_rt)
      {
        o << "    " << b.centre << " s: " << b.location << " +/- " << b.stderr_ppm
          << " ppm (n=" << b.n << ")\n";
      }
    }
    return o.str();
  }

} // namespace ODIA
