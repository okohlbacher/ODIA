// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/MobilityCalibration.h>

#include <odia/MassCalibration.h>

#include "RunProbe.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace ODIA
{
  namespace
  {
    // ------------------------------------------------------------------------
    // EVERY ROBUST ESTIMATE BELOW IS TAKEN IN MILLI-1/K0.
    //
    // The location, scale and peakedness estimators are MassCalibration's,
    // reused rather than re-derived -- they are statements about a narrow peak
    // on a near-uniform background and know nothing about ppm. But they are not
    // perfectly scale-free: refineLocation() calls a fit settled once successive
    // centres agree to 0.01 of whatever unit it was handed, which in ppm is a
    // hundredth of the effect and in 1/K0 would be the whole of it -- the loop
    // would stop after one pass and the flank bias it exists to remove would
    // stay.
    //
    // Multiplying by 1000 fixes that and does something better: it puts this
    // axis in the same NUMERIC regime the estimators were tuned on. A 1/K0
    // scatter of 0.010 becomes 10, a gate window of 0.05 becomes 50 -- the same
    // numbers as a 5 ppm scatter in a 50 ppm window. Nothing is being assumed
    // about the physics; the units are being made comparable so that a shared
    // estimator is shared honestly.
    // ------------------------------------------------------------------------
    constexpr double MILLI = 1000.0;

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

    /// One anchor as the fit sees it: an m/z, a delta in milli-1/K0, and the
    /// library index that decides which fold it belongs to.
    struct Anchor
    {
      double mz = 0.0;
      double delta = 0.0;
      std::uint32_t precursor = 0;
    };

    /// Equal-COUNT bins along m/z, each summarised by the median delta and that
    /// median's standard error.
    ///
    /// Equal count rather than equal width, for MassCalibration's reason: a
    /// tryptic library's m/z distribution is very far from uniform, and
    /// equal-width bins put a handful of points in the tails -- which is
    /// precisely where a shape fitted to them would do the most damage.
    std::vector<MobilityCalibration::Bin> binByMz(std::vector<Anchor> a, std::size_t max_bins,
                                                  std::size_t min_per_bin)
    {
      std::vector<MobilityCalibration::Bin> bins;
      if (a.size() < min_per_bin) { return bins; }
      std::sort(a.begin(), a.end(),
                [](const Anchor& x, const Anchor& y) { return x.mz < y.mz; });
      const std::size_t nb = std::min(max_bins, a.size() / min_per_bin);
      if (nb == 0) { return bins; }
      for (std::size_t b = 0; b < nb; ++b)
      {
        const std::size_t lo = b * a.size() / nb;
        const std::size_t hi = (b + 1) * a.size() / nb;
        if (hi <= lo) { continue; }
        std::vector<double> v;
        v.reserve(hi - lo);
        double key = 0.0;
        for (std::size_t i = lo; i < hi; ++i) { v.push_back(a[i].delta); key += a[i].mz; }
        MobilityCalibration::Bin bin;
        bin.n = hi - lo;
        bin.centre = key / static_cast<double>(bin.n);
        bin.location = medianOf(v);
        const double sd = madSigma(v, bin.location);
        // 1.2533 = sqrt(pi/2), the median's standard error against the mean's.
        bin.stderr_im = sd > 0.0 ? 1.2533 * sd / std::sqrt(static_cast<double>(bin.n)) : 0.0;
        bins.push_back(bin);
      }
      return bins;
    }

    /// How far the bin locations sit from their own weighted mean, in units of
    /// their own uncertainty, and how far apart the extremes are.
    ///
    /// This is the shape test. At chi2/dof ~ 1 the bins are consistent with
    /// noise about a constant and there is nothing to model; the swing is the
    /// absolute floor that stops a large anchor set from making every standard
    /// error small enough for chi-square alone to wave a shape through.
    void shapeEvidence(const std::vector<MobilityCalibration::Bin>& bins, double& chi2_per_dof,
                       double& swing, double& weighted_mean)
    {
      chi2_per_dof = 0.0;
      swing = 0.0;
      weighted_mean = 0.0;
      if (bins.size() < 2) { return; }
      double sw = 0.0, swy = 0.0;
      for (const auto& b : bins)
      {
        const double e = b.stderr_im > 0.0 ? b.stderr_im : 1.0;
        const double w = 1.0 / (e * e);
        sw += w; swy += w * b.location;
      }
      weighted_mean = sw > 0.0 ? swy / sw : 0.0;
      double chi2 = 0.0;
      double lo = bins.front().location, hi = bins.front().location;
      for (const auto& b : bins)
      {
        const double e = b.stderr_im > 0.0 ? b.stderr_im : 1.0;
        const double r = (b.location - weighted_mean) / e;
        chi2 += r * r;
        lo = std::min(lo, b.location);
        hi = std::max(hi, b.location);
      }
      chi2_per_dof = chi2 / static_cast<double>(bins.size() - 1);
      swing = hi - lo;
    }

    /// Everything one (charge, fold) subset produces.
    struct CurveFit
    {
      MobilityCalibration::Curve curve;
      std::vector<MobilityCalibration::Bin> bins;
      double chi2_per_dof = 0.0;
      double swing = 0.0;
      double location = 0.0;
      double sigma = 0.0;
      std::size_t core = 0;
    };

    /// Fit one charge's correction from one subset of its anchors.
    ///
    /// Everything here -- the robust centre, the scale, the core selection, the
    /// bins and the knots -- is computed from the subset ALONE. That is what
    /// makes the leave-one-fold-out correction out of sample in fact and not
    /// just in the part that looks like a fit: a core selected from all anchors
    /// and then fitted per fold would still have let every precursor vote on
    /// whether it was an outlier.
    CurveFit fitCurve(const std::vector<Anchor>& a, bool want_shape,
                      const MobilityCalibration::Options& opt, double gate_milli)
    {
      CurveFit out;
      out.curve.anchors = a.size();
      if (a.size() < 8) { return out; }

      std::vector<double> d;
      d.reserve(a.size());
      for (const auto& x : a) { d.push_back(x.delta); }
      std::sort(d.begin(), d.end());
      out.location = MassCalibration::refineLocation(
        d, MassCalibration::halfSampleMode(d), gate_milli);

      std::vector<double> dev;
      dev.reserve(d.size());
      for (double x : d) { dev.push_back(std::abs(x - out.location)); }
      std::sort(dev.begin(), dev.end());
      out.sigma = MassCalibration::backgroundCorrectedScale(dev, gate_milli);
      if (!(out.sigma > 0.0) || !std::isfinite(out.sigma))
      {
        out.sigma = madSigma(d, out.location);
      }
      if (!(out.sigma > 0.0)) { return out; }

      // The core. Outside it the sample is background, whose location is the
      // middle of the search window at every m/z -- including it would flatten
      // any real trend towards zero, which is the same reason MassCalibration
      // restricts its shape fit.
      std::vector<Anchor> core;
      core.reserve(a.size());
      for (const auto& x : a)
      {
        if (std::abs(x.delta - out.location) <= 4.0 * out.sigma) { core.push_back(x); }
      }
      out.core = core.size();
      if (core.size() < 8) { return out; }

      std::vector<double> cd;
      cd.reserve(core.size());
      for (const auto& x : core) { cd.push_back(x.delta); }
      const double clamp_milli = opt.max_correction_im * MILLI;
      out.curve.constant = std::clamp(medianOf(cd), -clamp_milli, clamp_milli) / MILLI;
      out.curve.supported = true;

      out.bins = binByMz(core, opt.max_mz_bins, opt.min_per_bin);
      double mean = 0.0;
      shapeEvidence(out.bins, out.chi2_per_dof, out.swing, mean);
      if (want_shape && out.bins.size() >= 2)
      {
        out.curve.shaped = true;
        for (const auto& b : out.bins)
        {
          out.curve.knot_mz.push_back(b.centre);
          out.curve.knot_offset.push_back(
            std::clamp(b.location, -clamp_milli, clamp_milli) / MILLI);
        }
      }
      return out;
    }
  } // namespace

  std::size_t MobilityCalibration::foldIndex(std::uint32_t precursor, std::size_t folds)
  {
    if (folds <= 1) { return 0; }
    // splitmix64's finaliser. A plain modulus of the index would put every
    // stride-sampled precursor of one library region in the same fold, which is
    // exactly the correlation the folds exist to break.
    std::uint64_t z = static_cast<std::uint64_t>(precursor) + 0x9E3779B97F4A7C15ull;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    z = z ^ (z >> 31);
    return static_cast<std::size_t>(z % folds);
  }

  // ==========================================================================
  // The fit.
  // ==========================================================================

  MobilityCalibration::Model MobilityCalibration::fit(
    const std::vector<MobilityResidual>& residuals, const Options& opt, Diagnostics*)
  {
    Model m;
    m.folds = std::max<std::size_t>(1, opt.folds);
    m.gate_peakedness = opt.min_peakedness;
    m.search_im = opt.search_im;
    const double gate_w = opt.gate_im > 0.0 ? std::min(opt.gate_im, opt.search_im)
                                            : opt.search_im;
    m.gate_im = gate_w;
    const double gate_milli = gate_w * MILLI;
    m.curves.assign((MAX_CHARGE + 1) * (m.folds + 1), Curve{});

    // Brightness selection, applied to the CONTROL on the same threshold so the
    // null keeps measuring the same procedure.
    double threshold = 0.0;
    if (opt.min_intensity_quantile > 0.0)
    {
      std::vector<double> intensity;
      intensity.reserve(residuals.size());
      for (const auto& r : residuals)
      {
        if (!r.decoy && std::isfinite(r.delta)) { intensity.push_back(r.intensity); }
      }
      if (!intensity.empty())
      {
        const std::size_t k = static_cast<std::size_t>(
          std::min(opt.min_intensity_quantile, 1.0) * static_cast<double>(intensity.size() - 1));
        std::nth_element(intensity.begin(), intensity.begin() + k, intensity.end());
        threshold = intensity[k];
      }
    }

    std::vector<MobilityResidual> kept;
    std::vector<double> target_milli, decoy_milli;
    kept.reserve(residuals.size());
    for (const auto& r : residuals)
    {
      if (!std::isfinite(r.delta)) { continue; }
      if (r.intensity < threshold) { continue; }
      if (std::abs(r.delta) > opt.search_im) { continue; }
      if (r.decoy) { decoy_milli.push_back(r.delta * MILLI); continue; }
      kept.push_back(r);
      target_milli.push_back(r.delta * MILLI);
    }
    m.residuals = target_milli.size();
    m.decoy_residuals = decoy_milli.size();

    if (static_cast<int>(target_milli.size()) < opt.min_residuals)
    {
      m.reason = "only " + std::to_string(target_milli.size()) + " residuals (need " +
                 std::to_string(opt.min_residuals) + ")";
      return m;
    }

    std::vector<double> sorted = target_milli;
    std::sort(sorted.begin(), sorted.end());
    const double location_milli = MassCalibration::refineLocation(
      sorted, MassCalibration::halfSampleMode(sorted), gate_milli);
    m.location = location_milli / MILLI;

    // ---- the gate is read on PER-CHARGE-CENTRED deltas ---------------------
    //
    // The gate asks one question: is there a 1/K0 peak in this run at all, or
    // is every cluster a coincidence? Where the peak SITS is a nuisance
    // parameter for that question -- and it is a nuisance parameter that
    // differs by charge, because charge is the one thing Mason-Schamp says the
    // mobility depends on directly. Pooling the charges and then looking for a
    // peak measures the wrong thing: two charges offset from each other put one
    // population into the statistic's own background band and the gate refuses
    // a run whose calibration is not merely present but large. Measured on the
    // synthetic suite: +0.020 at 2+ against -0.015 at 3+, which is as clean a
    // signal as this axis can have, reads a peakedness of 2.09 pooled and
    // refuses.
    //
    // Subtracting a centre cannot manufacture a peak -- a uniform sample
    // shifted is still uniform -- and the control gets exactly the same
    // treatment from its own deltas, so the null still measures this procedure
    // and not a simpler one.
    const auto centreByCharge = [&](const std::vector<const MobilityResidual*>& rs,
                                    std::vector<double>& centred) {
      std::vector<std::vector<double>> per(MAX_CHARGE + 1);
      for (const auto* r : rs)
      {
        per[std::min<std::size_t>(r->charge, MAX_CHARGE)].push_back(r->delta * MILLI);
      }
      std::vector<double> all;
      all.reserve(rs.size());
      for (const auto* r : rs) { all.push_back(r->delta * MILLI); }
      std::sort(all.begin(), all.end());
      const double pooled = all.empty() ? 0.0 : MassCalibration::refineLocation(
        all, MassCalibration::halfSampleMode(all), gate_milli);
      std::vector<double> centre(MAX_CHARGE + 1, pooled);
      for (std::size_t c = 0; c <= MAX_CHARGE; ++c)
      {
        if (per[c].size() < 32) { continue; }
        std::sort(per[c].begin(), per[c].end());
        centre[c] = MassCalibration::refineLocation(
          per[c], MassCalibration::halfSampleMode(per[c]), gate_milli);
      }
      centred.clear();
      centred.reserve(rs.size());
      for (const auto* r : rs)
      {
        centred.push_back(r->delta * MILLI - centre[std::min<std::size_t>(r->charge, MAX_CHARGE)]);
      }
      std::sort(centred.begin(), centred.end());
    };

    std::vector<const MobilityResidual*> target_ptr, decoy_ptr;
    for (const auto& r : residuals)
    {
      if (!std::isfinite(r.delta) || r.intensity < threshold) { continue; }
      if (std::abs(r.delta) > opt.search_im) { continue; }
      (r.decoy ? decoy_ptr : target_ptr).push_back(&r);
    }

    std::vector<double> gate_sample;
    centreByCharge(target_ptr, gate_sample);
    std::vector<double> dev;
    dev.reserve(gate_sample.size());
    for (double x : gate_sample) { dev.push_back(std::abs(x)); }
    m.peakedness = MassCalibration::peakednessRatio(dev, gate_milli);
    std::sort(dev.begin(), dev.end());
    const double gate_sigma = MassCalibration::backgroundCorrectedScale(dev, gate_milli);
    // Reported as "before" because it is what the m/z shape and the folds still
    // have to work on. It is NOT the raw scatter -- the per-charge centre has
    // already come out of it -- and the raw scatter is not a useful number to
    // print anyway, since on a run whose charges are offset from each other it
    // is a property of the mixture rather than of the axis. What the whole
    // correction is worth against the RAW residual is `squared_error_removed`.
    m.sigma_before = gate_sigma / MILLI;

    if (!decoy_ptr.empty())
    {
      std::vector<double> ds;
      centreByCharge(decoy_ptr, ds);
      std::vector<double> dd;
      dd.reserve(ds.size());
      for (double x : ds) { dd.push_back(std::abs(x)); }
      m.decoy_peakedness = MassCalibration::peakednessRatio(dd, gate_milli);
    }

    // ---- the gate ---------------------------------------------------------
    if (m.peakedness < opt.min_peakedness)
    {
      char buf[320];
      std::snprintf(buf, sizeof buf,
                    "the 1/K0 residuals are FLAT (peakedness %.2f < %.2f over %zu residuals) -- "
                    "that is what a sample of mostly-noise clusters looks like; a centre and a "
                    "width can still be computed from it and would be meaningless",
                    m.peakedness, opt.min_peakedness, target_milli.size());
      m.reason = buf;
      return m;
    }
    if (!(gate_sigma > 0.0) || !std::isfinite(gate_sigma) || !std::isfinite(location_milli))
    {
      m.reason = "degenerate scale";
      return m;
    }
    if (!decoy_milli.empty() && m.peakedness < opt.min_control_margin * m.decoy_peakedness)
    {
      char buf[360];
      std::snprintf(buf, sizeof buf,
                    "the control is as informative as the data: peakedness %.2f against %.2f, "
                    "a margin of %.2fx where %.2fx is required. Whatever agrees in mobility is "
                    "not specifically this precursor's fragments, so a centre fitted to it is a "
                    "property of the interference and not of the instrument",
                    m.peakedness, m.decoy_peakedness,
                    m.decoy_peakedness > 0.0 ? m.peakedness / m.decoy_peakedness : 0.0,
                    opt.min_control_margin);
      m.reason = buf;
      return m;
    }

    // ---- anchors and their folds ------------------------------------------
    for (const auto& r : kept) { m.anchors.push_back(r.precursor); }
    std::sort(m.anchors.begin(), m.anchors.end());
    m.anchors.erase(std::unique(m.anchors.begin(), m.anchors.end()), m.anchors.end());

    std::vector<std::vector<Anchor>> by_charge(MAX_CHARGE + 1);
    for (const auto& r : kept)
    {
      const std::size_t c = std::min<std::size_t>(r.charge, MAX_CHARGE);
      by_charge[c].push_back({r.mz, r.delta * MILLI, r.precursor});
    }

    // ---- the FORM, decided once per charge from all of that charge's anchors
    //
    // Whether the correction has an m/z shape is one bit, and it is decided from
    // the whole charge rather than per fold so that two precursors of the same
    // charge are not corrected by structurally different models for no reason.
    // The COEFFICIENTS -- which is what actually gets applied -- are still fitted
    // per fold, from anchors that exclude the precursor being corrected.
    bool any_supported = false, any_shaped = false;
    for (std::size_t c = 0; c <= MAX_CHARGE; ++c)
    {
      const auto& all = by_charge[c];
      if (all.empty()) { continue; }
      Model::ChargeReport rep;
      rep.charge = static_cast<std::uint8_t>(c);
      rep.anchors = all.size();
      rep.mz_low = std::numeric_limits<double>::max();
      for (const auto& x : all)
      {
        rep.mz_low = std::min(rep.mz_low, x.mz);
        rep.mz_high = std::max(rep.mz_high, x.mz);
      }

      if (all.size() < opt.min_anchors_per_charge)
      {
        char buf[220];
        std::snprintf(buf, sizeof buf,
                      "%zu anchors is below %zu, so charge %zu is left UNCORRECTED -- another "
                      "charge's offset is not this charge's answer",
                      all.size(), opt.min_anchors_per_charge, c);
        rep.reason = buf;
        m.by_charge.push_back(rep);
        continue;
      }

      const CurveFit whole = fitCurve(all, false, opt, gate_milli);
      rep.chi2_per_dof = whole.chi2_per_dof;
      rep.swing = whole.swing / MILLI;
      const bool shaped = whole.bins.size() >= 4 &&
                          whole.chi2_per_dof >= opt.min_shape_chi2 &&
                          whole.swing >= opt.min_shape_swing_im * MILLI;

      // One fit per fold, each from the anchors of the OTHER folds, plus one
      // from all of them for the precursors that were never anchors.
      std::size_t supported_folds = 0;
      for (std::size_t f = 0; f <= m.folds; ++f)
      {
        std::vector<Anchor> subset;
        subset.reserve(all.size());
        for (const auto& x : all)
        {
          if (f < m.folds && foldIndex(x.precursor, m.folds) == f) { continue; }
          subset.push_back(x);
        }
        const CurveFit cf = fitCurve(subset, shaped, opt, gate_milli);
        m.curves[c * (m.folds + 1) + f] = cf.curve;
        if (cf.curve.supported) { ++supported_folds; }
        if (f == m.folds) { rep.bins = cf.bins; }
      }

      const Curve& full = m.curves[c * (m.folds + 1) + m.folds];
      rep.supported = full.supported && supported_folds == m.folds + 1;
      if (!rep.supported)
      {
        // Partial support is refused outright: a run where some folds could be
        // fitted and others could not would correct precursors inconsistently
        // depending on a hash, which is not a calibration.
        for (std::size_t f = 0; f <= m.folds; ++f) { m.curves[c * (m.folds + 1) + f] = Curve{}; }
        rep.reason = "the fold models could not all be fitted, so charge " + std::to_string(c) +
                     " is left UNCORRECTED";
        m.by_charge.push_back(rep);
        continue;
      }
      rep.shaped = full.shaped;
      rep.constant = full.constant;
      rep.low = full.at(rep.mz_low);
      rep.high = full.at(rep.mz_high);
      any_supported = true;
      any_shaped = any_shaped || full.shaped;
      char buf[320];
      if (full.shaped)
      {
        std::snprintf(buf, sizeof buf,
                      "m/z-shaped over %zu bins: %+.4f at %.0f Th to %+.4f at %.0f Th "
                      "(chi2/dof %.1f >= %.1f, swing %.4f >= %.4f)",
                      rep.bins.size(), rep.low, rep.mz_low, rep.high, rep.mz_high,
                      rep.chi2_per_dof, opt.min_shape_chi2, rep.swing, opt.min_shape_swing_im);
      }
      else
      {
        std::snprintf(buf, sizeof buf,
                      "constant %+.4f: the m/z bins are flat enough that a shape would be "
                      "fitting bin noise (chi2/dof %.1f against %.1f, swing %.4f against %.4f)",
                      rep.constant, rep.chi2_per_dof, opt.min_shape_chi2, rep.swing,
                      opt.min_shape_swing_im);
      }
      rep.reason = buf;
      m.by_charge.push_back(rep);
    }

    if (!any_supported)
    {
      m.reason = "the gate passed, but no charge has enough anchors to fit from "
                 "(min_anchors_per_charge = " + std::to_string(opt.min_anchors_per_charge) +
                 "), so the library 1/K0 is left as it is";
      m.curves.clear();
      return m;
    }

    m.fitted = true;
    m.form = any_shaped ? "mz_shaped" : "constant";

    // ---- what it is worth, measured OUT OF FOLD ---------------------------
    // Every anchor is scored with the model that never saw it, which is also
    // exactly how it will be corrected at extraction time.
    double sum_before = 0.0, sum_after = 0.0;
    std::vector<double> corrected;
    corrected.reserve(kept.size());
    for (const auto& r : kept)
    {
      const double off = m.offsetFor(r.precursor, r.mz, r.charge);
      const double a = r.delta, b = r.delta - off;
      sum_before += a * a;
      sum_after += b * b;
      corrected.push_back(b * MILLI);
    }
    m.squared_error_removed = sum_before > 0.0 ? 1.0 - sum_after / sum_before : 0.0;

    std::sort(corrected.begin(), corrected.end());
    const double cmode = MassCalibration::refineLocation(
      corrected, MassCalibration::halfSampleMode(corrected), gate_milli);
    std::vector<double> cdev;
    cdev.reserve(corrected.size());
    for (double x : corrected) { cdev.push_back(std::abs(x - cmode)); }
    std::sort(cdev.begin(), cdev.end());
    m.sigma_after = MassCalibration::backgroundCorrectedScale(cdev, gate_milli) / MILLI;

    char buf[420];
    std::snprintf(buf, sizeof buf,
                  "%zu of %zu charges corrected, %s; %.0f%% of the mean squared 1/K0 error "
                  "removed out of fold, robust scatter %.4f -> %.4f",
                  static_cast<std::size_t>(
                    std::count_if(m.by_charge.begin(), m.by_charge.end(),
                                  [](const Model::ChargeReport& r) { return r.supported; })),
                  m.by_charge.size(), m.form.c_str(), 100.0 * m.squared_error_removed,
                  m.sigma_before, m.sigma_after);
    m.reason = buf;

    // ---- the window this residual would support: REPORTED, NEVER APPLIED ---
    if (m.sigma_after > 0.0 && std::isfinite(m.sigma_after))
    {
      m.window_im = std::max(opt.sigma_multiple * m.sigma_after, opt.floor_im);
      // An estimate can never legitimately exceed the width it was measured
      // over; that is a fit to the edge of the window, not to the instrument.
      if (m.window_im > gate_w) { m.window_im = -1.0; }
    }
    return m;
  }

  // ==========================================================================
  // Collection: probing the run for a 1/K0 each precursor actually shows.
  // ==========================================================================

  namespace
  {
    /// One query: a fragment m/z (possibly shifted, for a control cell) and the
    /// cell it belongs to. `cell` is LOCAL to the window, so everything the
    /// match loop touches is an array index rather than a hash lookup.
    struct Query
    {
      double mz;
      std::uint32_t cell;
      std::uint8_t fragment;
    };

    struct WindowIndex
    {
      std::vector<Query> query;                ///< ascending in m/z
      std::vector<std::uint32_t> bucket;       ///< log-m/z bucket -> first query
      std::vector<float> cell_im;              ///< local cell -> library 1/K0
      std::vector<std::uint32_t> cell_slot;    ///< local cell -> sampled slot
      std::vector<std::uint8_t> cell_variant;  ///< local cell -> 0 target, else control
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
      }
    };

    /// One peak that matched one query, kept with its mobility.
    struct Match
    {
      std::uint32_t cell;
      float im;
      float intensity;
      std::uint8_t fragment;
    };

    /// The best REPRODUCED mobility cluster found for one cell.
    struct Cell
    {
      float rt = 0.0f;
      float total = 0.0f;      ///< summed intensity over the cycles that agreed
      float im = 0.0f;         ///< intensity-weighted 1/K0 over those cycles
      std::uint8_t n = 0;      ///< fragments that agreed, in the best cycle
      std::uint8_t cycles = 0; ///< cycles of the block that agreed
    };
  } // namespace

  std::vector<MobilityResidual> MobilityCalibration::collect(const Library& library,
                                                             SpectrumSource& source,
                                                             const Options& opt,
                                                             Diagnostics* diag)
  {
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<MobilityResidual> out;
    if (diag) { *diag = Diagnostics{}; }

    const auto& p = library.precursors();
    const auto& tr = library.transitions();
    const auto& windows = source.windows();
    const auto& info = source.spectra();
    if (windows.empty() || info.empty() || library.precursorCount() == 0) { return out; }

    const std::size_t max_frag = std::clamp<std::size_t>(opt.max_fragments, 1, 255);

    const auto window_of = RunProbe::windowOfSpectrum(source);
    const auto cycles = RunProbe::cycles(source, window_of);
    if (cycles.empty()) { return out; }

    // require_im: a precursor with no library 1/K0 has no residual to give, and
    // a library with none at all is the second no-op this class has to name.
    auto sampled = RunProbe::samplePrecursors(library, opt.max_precursors,
                                              opt.min_fragments_matched, true, false);
    const std::size_t n_targets = sampled.size();
    if (diag)
    {
      diag->precursors_sampled = n_targets;
      diag->precursors_with_library_im = n_targets;
      diag->library_has_mobility = n_targets != 0;
    }
    if (sampled.empty()) { return out; }

    // ---- the null ----------------------------------------------------------
    // The library's own decoys where it has them, an m/z shift where it does
    // not. See Options::control for the measurement that made that the order.
    std::vector<double> shifts;
    if (opt.control == Options::Control::LibraryDecoys)
    {
      const auto controls = RunProbe::samplePrecursors(library, opt.max_precursors,
                                                       opt.min_fragments_matched, true, true);
      sampled.insert(sampled.end(), controls.begin(), controls.end());
      if (diag)
      {
        diag->control_precursors = controls.size();
        diag->control_kind = controls.empty() ? "none (the library has no decoys)"
                                              : "the library's own decoy precursors";
      }
      if (controls.empty() && !opt.decoy_shifts.empty())
      {
        shifts = opt.decoy_shifts;
        if (diag) { diag->control_kind = "m/z-shifted fragments (no library decoys)"; }
      }
    }
    else if (opt.control == Options::Control::MzShift)
    {
      shifts = opt.decoy_shifts;
      if (diag) { diag->control_kind = "m/z-shifted fragments"; }
    }
    const std::size_t variants = 1 + shifts.size();
    const auto& lib_decoy = p.decoy;

    // ---- one window each, and that window's query index --------------------
    std::vector<WindowIndex> index(windows.size());
    std::vector<std::uint32_t> window_of_precursor(sampled.size(), RunProbe::NO_WINDOW);
    for (std::size_t s = 0; s < sampled.size(); ++s)
    {
      window_of_precursor[s] = RunProbe::centralWindow(windows, fromFixed(p.mz[sampled[s]]));
    }

    const auto correctedMz = [&](double theoretical) {
      double ppm = opt.fragment_ppm_offset;
      if (opt.fragment_ppm_log_slope != 0.0 && opt.fragment_ppm_ref_mz > 0.0)
      {
        ppm += opt.fragment_ppm_log_slope * std::log(theoretical / opt.fragment_ppm_ref_mz);
      }
      if (opt.fragment_ppm_slope_per_1000 != 0.0)
      {
        ppm += opt.fragment_ppm_slope_per_1000 *
               (theoretical - opt.fragment_ppm_ref_mz) / 1000.0;
      }
      return theoretical * (1.0 + ppm * 1e-6);
    };

    // Cells are numbered per window and densely, so the per-spectrum work below
    // is all array indexing.
    std::vector<std::uint32_t> cell_base(windows.size(), 0);
    for (std::size_t s = 0; s < sampled.size(); ++s)
    {
      const std::uint32_t w = window_of_precursor[s];
      if (w == RunProbe::NO_WINDOW) { continue; }
      auto& x = index[w];
      const std::size_t i = sampled[s];
      const std::size_t nfrag = std::min<std::size_t>(p.transition_count[i], max_frag);
      for (std::size_t v = 0; v < variants; ++v)
      {
        const double shift = v == 0 ? 0.0 : shifts[v - 1];
        const auto local = static_cast<std::uint32_t>(x.cell_im.size());
        bool any = false;
        for (std::size_t k = 0; k < nfrag; ++k)
        {
          const MzFixed fixed = tr.product_mz[p.transition_begin[i] + k];
          if (fixed == MZ_INVALID) { continue; }
          const double mz = correctedMz(fromFixed(fixed)) + shift;
          if (!(mz > 1.0)) { continue; }
          x.query.push_back({mz, local, static_cast<std::uint8_t>(k)});
          any = true;
        }
        if (!any) { continue; }
        x.cell_im.push_back(p.im[i]);
        x.cell_slot.push_back(static_cast<std::uint32_t>(s));
        x.cell_variant.push_back(static_cast<std::uint8_t>(v));
      }
    }
    std::size_t total_cells = 0;
    for (std::size_t w = 0; w < windows.size(); ++w)
    {
      cell_base[w] = static_cast<std::uint32_t>(total_cells);
      total_cells += index[w].cell_im.size();
      index[w].build(opt.fragment_ppm);
    }
    if (total_cells == 0) { return out; }

    // CONTIGUOUS BLOCKS, not scattered cycles -- see Options::cycle_block for the
    // measurement that says the scattered version cannot work on this axis.
    const std::size_t block_len = std::max<std::size_t>(1, opt.cycle_block);
    const std::size_t want_blocks = std::max<std::size_t>(1, opt.cycles / block_len);
    const auto starts = RunProbe::stratifiedBlocks(cycles.size(), want_blocks, block_len,
                                                   opt.sample_seed);
    if (diag) { diag->blocks = starts.size(); }

    std::vector<Cell> cell(total_cells);
    std::vector<char> ever_usable(sampled.size(), 0);
    std::vector<SpectrumPeaks> block;
    std::vector<Match> match;
    std::vector<char> usable;
    std::vector<float> frag_im(max_frag), frag_int(max_frag);
    std::vector<double> band_widths;
    // Per cell, what each cycle of the CURRENT block saw. This is the whole
    // point of the block: an elution profile the noise has to reproduce.
    std::vector<float> seen_im(total_cells * block_len, 0.0f);
    std::vector<float> seen_int(total_cells * block_len, 0.0f);
    std::vector<std::uint8_t> seen_n(total_cells * block_len, 0);
    std::vector<std::uint32_t> touched;

    std::size_t decoded = 0, with_mobility = 0, probed = 0;
    for (std::size_t bi = 0; bi < starts.size(); ++bi)
    {
      const std::size_t first_cycle = starts[bi];
      const std::size_t last_cycle = std::min(first_cycle + block_len, cycles.size());
      const std::size_t lo_spec = cycles[first_cycle].first;
      const std::size_t hi_spec = cycles[last_cycle - 1].second;
      // One range request per block: the cycles are contiguous, which is what
      // makes reproducibility cheaper than scattering the same cycle budget.
      source.peaks(lo_spec, hi_spec, block);
      decoded += hi_spec - lo_spec;

      for (auto t : touched) { seen_im[t] = 0.0f; seen_int[t] = 0.0f; seen_n[t] = 0; }
      touched.clear();

      std::size_t j = 0;
      for (std::size_t si = lo_spec; si < hi_spec; ++si)
      {
        while (j + 1 < last_cycle - first_cycle && si >= cycles[first_cycle + j].second) { ++j; }
        const std::uint32_t w = window_of[si];
        if (w == RunProbe::NO_WINDOW) { continue; }
        auto& x = index[w];
        if (x.query.empty()) { continue; }
        const auto& peaks = block[si - lo_spec];
        if (peaks.empty()) { continue; }
        // A spectrum with no per-peak 1/K0 says nothing about the mobility
        // axis. On a run where NO spectrum has one there is no axis at all,
        // which is a different fact from a failed fit -- see run_has_mobility.
        if (!peaks.hasIonMobility()) { continue; }
        ++with_mobility;

        const double im_low = info[si].window.im_low, im_high = info[si].window.im_high;
        const bool bounded = std::isfinite(im_low) && std::isfinite(im_high);
        if (bounded) { band_widths.push_back(im_high - im_low); }

        // Which cells this spectrum can produce an UNBIASED residual for: the
        // whole of the symmetric search interval about the library 1/K0 has to
        // lie inside the frame's band, or the residual is truncated on one side
        // only and its centre is pulled inwards.
        usable.assign(x.cell_im.size(), 1);
        for (std::size_t lc = 0; lc < x.cell_im.size(); ++lc)
        {
          const float want = x.cell_im[lc];
          if (!std::isfinite(want)) { usable[lc] = 0; continue; }
          if (bounded && (want - opt.search_im < im_low || want + opt.search_im > im_high))
          {
            usable[lc] = 0;
            continue;
          }
          ever_usable[x.cell_slot[lc]] = 1;
        }

        match.clear();
        const double ppm = opt.fragment_ppm * 1e-6;
        const double front = x.query.front().mz, back = x.query.back().mz;
        for (std::size_t k = 0; k < peaks.size(); ++k)
        {
          const double mzp = peaks.mz[k];
          const double slack = mzp * ppm * 1.01 + 1e-6;
          if (mzp + slack < front || mzp - slack > back) { continue; }
          const float peak_im = peaks.ion_mobility[k];
          if (bounded && (peak_im < im_low || peak_im >= im_high)) { continue; }
          const float intensity = peaks.intensity[k];
          std::size_t i = x.bucket[x.bucketOf(std::max(mzp - slack, front))];
          for (; i < x.query.size() && x.query[i].mz <= mzp + slack; ++i)
          {
            const double q = x.query[i].mz;
            if (std::abs(mzp - q) > q * ppm) { continue; }
            const std::uint32_t lc = x.query[i].cell;
            if (!usable[lc]) { continue; }
            // The ONLY mobility prior used: a symmetric search about the
            // library value, which truncates the residual without moving it.
            if (std::abs(peak_im - x.cell_im[lc]) > opt.search_im) { continue; }
            match.push_back({lc, peak_im, intensity, x.query[i].fragment});
          }
        }
        if (match.empty()) { continue; }

        std::sort(match.begin(), match.end(), [](const Match& a, const Match& b) {
          return a.cell != b.cell ? a.cell < b.cell : a.im < b.im;
        });

        // ---- the cluster: fragments of one precursor agreeing in mobility ---
        for (std::size_t lo = 0; lo < match.size();)
        {
          std::size_t hi = lo;
          while (hi < match.size() && match[hi].cell == match[lo].cell) { ++hi; }
          const std::uint32_t lc = match[lo].cell;
          ++probed;

          float best_total = 0.0f, best_im = 0.0f;
          std::uint8_t best_n = 0;
          for (std::size_t a = lo; a < hi; ++a)
          {
            const float top = match[a].im + 2.0f * static_cast<float>(opt.cluster_im);
            std::fill(frag_int.begin(), frag_int.end(), 0.0f);
            std::fill(frag_im.begin(), frag_im.end(), 0.0f);
            std::uint8_t n = 0;
            float total = 0.0f;
            for (std::size_t b = a; b < hi && match[b].im <= top; ++b)
            {
              const std::uint8_t f = match[b].fragment;
              if (f >= max_frag) { continue; }
              // The most intense peak of this fragment inside the slice. A TIMS
              // ion is spread over several mobility scans, so a fragment appears
              // many times and only its apex should weigh.
              if (match[b].intensity > frag_int[f])
              {
                if (frag_int[f] == 0.0f) { ++n; }
                total += match[b].intensity - frag_int[f];
                frag_int[f] = match[b].intensity;
                frag_im[f] = match[b].im;
              }
            }
            if (n < opt.min_fragments_matched) { continue; }
            if (total <= best_total) { continue; }
            double num = 0.0, den = 0.0;
            for (std::size_t f = 0; f < max_frag; ++f)
            {
              if (frag_int[f] <= 0.0f) { continue; }
              num += double(frag_int[f]) * double(frag_im[f]);
              den += double(frag_int[f]);
            }
            if (!(den > 0.0)) { continue; }
            best_total = total;
            best_n = n;
            best_im = static_cast<float>(num / den);
          }

          if (best_n >= opt.min_fragments_matched)
          {
            const std::uint32_t at =
              static_cast<std::uint32_t>((cell_base[w] + lc) * block_len + j);
            if (seen_int[at] == 0.0f) { touched.push_back(at); }
            if (best_total > seen_int[at])
            {
              seen_int[at] = best_total;
              seen_im[at] = best_im;
              seen_n[at] = best_n;
            }
          }
          lo = hi;
        }
      }

      // ---- consolidate the block: what was REPRODUCED across its cycles -----
      //
      // A precursor that is really eluting is in every cycle of a block that
      // overlaps its peak, at the same mobility. A coincidence is in one. This
      // is the test that a pre-identification probe of the mobility axis has to
      // pass and that a brightest-cell-anywhere probe does not.
      std::vector<std::uint32_t> cells_touched;
      cells_touched.reserve(touched.size());
      for (auto t : touched) { cells_touched.push_back(t / static_cast<std::uint32_t>(block_len)); }
      std::sort(cells_touched.begin(), cells_touched.end());
      cells_touched.erase(std::unique(cells_touched.begin(), cells_touched.end()),
                          cells_touched.end());
      for (std::uint32_t gc : cells_touched)
      {
        const std::size_t base = std::size_t(gc) * block_len;
        std::uint8_t best_k = 0, best_frag = 0;
        double best_sum = 0.0, best_im = 0.0;
        for (std::size_t a = 0; a < block_len; ++a)
        {
          if (seen_int[base + a] <= 0.0f) { continue; }
          std::uint8_t k = 0;
          double sum = 0.0, num = 0.0;
          std::uint8_t frag = 0;
          for (std::size_t b = 0; b < block_len; ++b)
          {
            if (seen_int[base + b] <= 0.0f) { continue; }
            if (std::abs(seen_im[base + b] - seen_im[base + a]) > opt.cluster_im) { continue; }
            ++k;
            sum += seen_int[base + b];
            num += double(seen_int[base + b]) * double(seen_im[base + b]);
            frag = std::max(frag, seen_n[base + b]);
          }
          if (k < best_k || (k == best_k && sum <= best_sum)) { continue; }
          best_k = k; best_sum = sum; best_im = num / sum; best_frag = frag;
        }
        if (best_k < opt.min_cycles_matched) { continue; }
        Cell& cl = cell[gc];
        // Across blocks the brightest reproduced elution wins, which with no
        // retention-time map is the best available stand-in for the apex.
        if (static_cast<float>(best_sum) <= cl.total) { continue; }
        cl.total = static_cast<float>(best_sum);
        cl.im = static_cast<float>(best_im);
        cl.n = best_frag;
        cl.cycles = best_k;
        cl.rt = static_cast<float>(info[cycles[first_cycle].first].retention_time);
      }

      // A run with no mobility at all is not worth decoding the rest of.
      if (bi == 0 && with_mobility == 0) { break; }
    }

    for (std::size_t w = 0; w < windows.size(); ++w)
    {
      const auto& x = index[w];
      for (std::size_t lc = 0; lc < x.cell_im.size(); ++lc)
      {
        const Cell& cl = cell[cell_base[w] + lc];
        if (cl.n == 0) { continue; }
        const std::size_t i = sampled[x.cell_slot[lc]];
        MobilityResidual r;
        r.precursor = static_cast<std::uint32_t>(i);
        r.mz = static_cast<float>(fromFixed(p.mz[i]));
        r.charge = p.charge[i];
        r.im_library = x.cell_im[lc];
        r.im_observed = cl.im;
        r.delta = cl.im - x.cell_im[lc];
        r.intensity = cl.total;
        r.rt = cl.rt;
        r.fragments = cl.n;
        r.cycles = cl.cycles;
        r.decoy = x.cell_variant[lc] != 0 || lib_decoy[i] != 0;
        out.push_back(r);
      }
    }

    if (diag)
    {
      diag->cells_probed = probed;
      diag->spectra_decoded = decoded;
      diag->spectra_with_mobility = with_mobility;
      diag->run_has_mobility = with_mobility > 0;
      for (const auto& r : out) { (r.decoy ? diag->decoy_cells : diag->target_cells) += 1; }
      diag->precursors_outside_band = static_cast<std::size_t>(
        std::count(ever_usable.begin(), ever_usable.begin() + n_targets, 0));
      if (!band_widths.empty())
      {
        std::nth_element(band_widths.begin(), band_widths.begin() + band_widths.size() / 2,
                         band_widths.end());
        diag->band_width = band_widths[band_widths.size() / 2];
      }
      diag->collect_seconds = std::chrono::duration<double>(
                                std::chrono::steady_clock::now() - t0).count();
    }
    return out;
  }

  MobilityCalibration::Model MobilityCalibration::calibrate(const Library& library,
                                                            SpectrumSource& source,
                                                            const Options& opt, Diagnostics* diag)
  {
    Diagnostics local;
    Diagnostics& d = diag ? *diag : local;
    const auto residuals = collect(library, source, opt, &d);

    // The two no-ops come first and are named, because "there is no 1/K0 axis
    // here" is not a failed calibration and must not be logged as one.
    if (!d.run_has_mobility)
    {
      Model m;
      m.run_has_mobility = false;
      m.form = "none";
      m.reason = "the run carries no ion mobility (no spectrum in the probed cycles has a "
                 "per-peak 1/K0), so there is no mobility axis to calibrate and extraction "
                 "is left exactly as it was";
      return m;
    }
    if (!d.library_has_mobility)
    {
      Model m;
      m.library_has_mobility = false;
      m.form = "none";
      m.reason = "the library carries no 1/K0, so there is nothing to correct and extraction "
                 "is left exactly as it was";
      return m;
    }
    return fit(residuals, opt, &d);
  }

  std::string MobilityCalibration::report(const Model& m, const Diagnostics* diag)
  {
    std::ostringstream o;
    o.setf(std::ios::fixed);
    o.precision(4);
    const char* verdict = !m.run_has_mobility || !m.library_has_mobility
                            ? "NO MOBILITY AXIS"
                            : (m.fitted ? "GATE PASSED" : "GATE FAILED");
    o << "ion-mobility calibration: " << verdict << " -- " << m.reason << "\n";
    if (!m.run_has_mobility || !m.library_has_mobility) { return o.str(); }

    o << "  residuals:   " << m.residuals << " target, " << m.decoy_residuals << " control"
      << (diag ? " from " + diag->control_kind : std::string());
    if (diag && diag->spectra_decoded)
    {
      o << " (from " << diag->precursors_sampled << " sampled precursors over "
        << diag->spectra_decoded << " spectra, " << diag->spectra_with_mobility
        << " with 1/K0, ";
      o.precision(2);
      o << diag->collect_seconds << " s)";
      o.precision(4);
    }
    o << "\n";
    if (diag && diag->precursors_outside_band)
    {
      o << "  band margin: " << diag->precursors_outside_band << " of "
        << diag->precursors_sampled
        << " precursors never had +/-" << m.search_im
        << " of clear band around their library 1/K0 (median band "
        << diag->band_width << ") and were not used\n";
    }
    o.precision(2);
    o << "  peakedness:  " << m.peakedness << " (control " << m.decoy_peakedness
      << "), gate at " << m.gate_peakedness << " evaluated over +/-";
    o.precision(4);
    o << m.gate_im << " of a +/-" << m.search_im << " search\n";
    if (m.fitted)
    {
      o << "  centre:      " << m.location << " overall; scatter about the per-charge "
           "centre " << m.sigma_before << " before the m/z shape, " << m.sigma_after
        << " after\n";
      o.precision(1);
      o << "  removed:     " << 100.0 * m.squared_error_removed
        << "% of the mean squared 1/K0 error, measured OUT OF FOLD over " << m.folds
        << " folds\n";
      o.precision(4);
      if (m.window_im > 0.0)
      {
        o << "  window:      the corrected residual would support +/-" << m.window_im
          << " -- REPORTED ONLY, the extraction window is not changed here\n";
      }
    }
    for (const auto& c : m.by_charge)
    {
      o << "  charge " << int(c.charge) << ":    " << c.anchors << " anchors, "
        << (c.supported ? "" : "NOT APPLIED -- ") << c.reason << "\n";
      for (const auto& b : c.bins)
      {
        o.precision(0);
        o << "    " << b.centre << " Th: ";
        o.precision(4);
        o << b.location << " +/- " << b.stderr_im << " (n=" << b.n << ")\n";
      }
    }
    return o.str();
  }

} // namespace ODIA
