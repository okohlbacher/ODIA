#pragma once

#include <odia/PeakGroupScorer.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace ODIA
{
  /// Size the fragment m/z window from the run's own identifications.
  ///
  /// WHY THIS EXISTS. `MassCalibration` already detects the right width per
  /// file, and when its gate passes it beats any constant -- on IH1 the
  /// calibrated ~10 ppm gives 1,306 identifications against 922 for a hand-set
  /// 15. The defect was never the calibration, it was the FALLBACK: a failed
  /// gate meant "use 15 ppm", a middling width chosen for no reason. Raising
  /// that to 50 was the right direction (not knowing the error argues for wide)
  /// but it is still a constant, just a better-chosen one -- and the two
  /// instruments want opposite widths, so no constant can be right for both.
  ///
  /// This is the measurement that replaces the constant. It runs on peak groups
  /// the FDR has ALREADY ACCEPTED, which is what every standalone probe failed
  /// at: a probe over the whole library measures mostly interference, whereas
  /// the contamination here is bounded by the FDR threshold itself. Pass 1
  /// already extracts wide and already scores, so the population exists at no
  /// extra cost beyond the two ppm planes.
  ///
  /// TWO TRAPS, both guarded below.
  ///
  /// 1. THE UNIT. `PeakGroup::mass_ppm` is a median over (fragment x cycle)
  ///    cells; its precision is sigma/sqrt(N_eff), so the spread of group
  ///    medians across groups is several times narrower than the spread of the
  ///    fragments a window has to admit. We therefore pool `mass_ppm_spread`,
  ///    which is measured BETWEEN FRAGMENTS WITHIN a group, and never the
  ///    scatter of `mass_ppm` between groups.
  ///
  /// 2. CENSORING. A residual can only be observed where a peak matched, and a
  ///    peak can only match inside the window that was extracted. So the
  ///    observed distribution is truncated at the current half-width, every
  ///    estimate drawn from it is biased narrow, and narrowing on that estimate
  ///    and re-measuring would spiral inward run after run. `censored` is set
  ///    when the distribution is not comfortably inside the window it was
  ///    measured through, and a censored estimate must not be applied.
  struct MassWidth
  {
    /// Fragments per group below which the group's spread is not estimated.
    static constexpr std::size_t min_fragments = 3;

    struct Estimate
    {
      bool valid = false;
      /// Groups that contributed.
      std::size_t groups = 0;
      /// Median of the accepted groups' own median deviations, ppm. The CENTRE.
      double centre_ppm = std::numeric_limits<double>::quiet_NaN();
      /// Median per-fragment robust sigma, ppm. The SCALE.
      double sigma_ppm = std::numeric_limits<double>::quiet_NaN();
      /// 95th percentile of |per-fragment sigma|, ppm; the censoring probe.
      double p95_sigma_ppm = std::numeric_limits<double>::quiet_NaN();
      /// The half-width the run was extracted at when this was measured.
      double measured_through_ppm = 0.0;
      /// True when the observed scatter runs into the extraction window, so the
      /// estimate is truncated and must not be acted on.
      bool censored = false;
      /// The recommended half-width, `sigmas` x `sigma_ppm`, unclamped.
      double width_ppm = std::numeric_limits<double>::quiet_NaN();
    };

    /// Pool the accepted target groups of a scored pass.
    ///
    /// `q_max` selects the accepted set, `extracted_ppm` is the half-width they
    /// were extracted through, and `sigmas` converts the robust scale into a
    /// half-width. `min_groups` is a floor below which the pooled estimate is
    /// not worth trusting.
    static Estimate measure(const PeakGroupScorer::Result& result, double q_max,
                            double extracted_ppm, double sigmas,
                            std::size_t min_groups)
    {
      Estimate e;
      e.measured_through_ppm = extracted_ppm;
      if (!result.fdr_valid) { return e; }

      // One row per PRECURSOR, the best-scoring accepted candidate. Peak groups
      // are per candidate and a precursor can contribute several; pooling them
      // all would weight precursors by how ambiguous they were.
      std::vector<const PeakGroupScorer::PeakGroup*> best;
      for (const auto& g : result.groups)
      {
        if (g.decoy || g.qvalue > q_max) { continue; }
        if (!(g.mass_ppm_frags >= min_fragments)) { continue; }
        if (!std::isfinite(g.mass_ppm) || !std::isfinite(g.mass_ppm_spread)) { continue; }
        if (best.size() <= g.precursor) { best.resize(g.precursor + 1, nullptr); }
        auto*& b = best[g.precursor];
        if (b == nullptr || g.dscore > b->dscore) { b = &g; }
      }

      std::vector<double> centres, sigmas_v;
      for (const auto* g : best)
      {
        if (g == nullptr) { continue; }
        centres.push_back(static_cast<double>(g->mass_ppm));
        sigmas_v.push_back(static_cast<double>(g->mass_ppm_spread));
      }
      e.groups = centres.size();
      if (e.groups < min_groups) { return e; }

      const auto quantile = [](std::vector<double>& v, double p) {
        const std::size_t i = std::min(v.size() - 1,
          static_cast<std::size_t>(p * static_cast<double>(v.size())));
        std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(i), v.end());
        return v[i];
      };

      e.centre_ppm = quantile(centres, 0.5);
      e.sigma_ppm = quantile(sigmas_v, 0.5);
      e.p95_sigma_ppm = quantile(sigmas_v, 0.95);
      e.width_ppm = sigmas * e.sigma_ppm;

      // The censoring test. A window of half-width W truncates every deviation
      // at |W|, so a scatter whose upper tail approaches W was measured through
      // a wall. 0.5 rather than something tighter because the quantity compared
      // is a robust sigma, not a raw deviation: at 0.5 W the 95th-percentile
      // group's own fragments are already reaching the window edge.
      if (!(extracted_ppm > 0.0) ||
          !(e.p95_sigma_ppm < 0.5 * extracted_ppm) ||
          !std::isfinite(e.sigma_ppm) || !(e.sigma_ppm > 0.0))
      {
        e.censored = true;
      }
      e.valid = true;
      return e;
    }

    static std::string report(const Estimate& e)
    {
      std::ostringstream os;
      os.setf(std::ios::fixed);
      os.precision(2);
      if (!e.valid)
      {
        os << "fragment window from identifications: not measured ("
           << (e.groups == 0 ? "no accepted group carried per-fragment residuals"
                             : "too few accepted groups")
           << "; " << e.groups << " usable)";
        return os.str();
      }
      os << "fragment window from identifications: " << e.groups
         << " accepted precursors, centre " << e.centre_ppm << " ppm, per-fragment sigma "
         << e.sigma_ppm << " ppm (p95 " << e.p95_sigma_ppm << "), measured through +/-"
         << e.measured_through_ppm << " ppm -> suggested half-width " << e.width_ppm
         << " ppm";
      if (e.censored)
      {
        os << " [CENSORED: the scatter reaches the window it was measured "
              "through, so this estimate is biased narrow and is not applied]";
      }
      return os.str();
    }
  };
} // namespace ODIA
