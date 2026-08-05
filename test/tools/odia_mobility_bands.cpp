// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// deriveMobilityBands: the guess that separates two isolation windows packed
// into one diaPASEF frame.
//
// It is a guess about the instrument -- the file states no band, only a
// mobility POSITION per window, mis-attached to the wrong precursor -- and a
// wrong guess rejects real signal in a way that looks like the instrument
// rather than like a bug. So the interesting cases are the ones where it must
// refuse, and where the attribution must not depend on the order the windows
// happen to be listed in.
//
// The derivation used to live inside MzPeakSource's constructor, where nothing
// could reach it without a file.

#include <odia/MobilityBands.h>

#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

namespace
{
  int failures = 0;

  void check(bool ok, const std::string& what)
  {
    std::printf("%s %s\n", ok ? "  ok  " : "  FAIL", what.c_str());
    if (!ok) { ++failures; }
  }

  void checkResult(ODIA::MobilityBandResult got, ODIA::MobilityBandResult want,
                   const std::string& what)
  {
    const bool ok = got == want;
    std::printf("%s %s (want %s, got %s)\n", ok ? "  ok  " : "  FAIL", what.c_str(),
                ODIA::toString(want), ODIA::toString(got));
    if (!ok) { ++failures; }
  }

  ODIA::IsolationWindow window(double centre, double width)
  {
    ODIA::IsolationWindow w;
    w.mz_low = centre - 0.5 * width;
    w.mz_high = centre + 0.5 * width;
    return w;
  }

  // a == b first, so that an open end (+/-inf) compares equal to itself:
  // inf - inf is NaN, and NaN < 1e-9 is false.
  bool near(double a, double b) { return a == b || std::abs(a - b) < 1e-9; }
  const double INF = std::numeric_limits<double>::infinity();
}

int main()
{
  // ------------------------------------------------------------------ shape
  {
    // Two windows, two positions, each naming its own centre. The split is the
    // midpoint; the outer ends stay open, because the frame says nothing about
    // where the mobility range stops.
    std::vector<ODIA::IsolationWindow> w{window(500.0, 25.0), window(525.0, 25.0)};
    const std::vector<ODIA::MobilityPosition> ions{{500.0, 1.10}, {525.0, 0.90}};
    checkResult(ODIA::deriveMobilityBands(ions, w), ODIA::MobilityBandResult::Derived,
                "two windows named by two positions");
    // Sorted by MOBILITY, not by m/z: window 1 is the lower band here.
    check(near(w[1].im_low, -INF) && near(w[1].im_high, 1.00), "lower band is (-inf, 1.00)");
    check(near(w[0].im_low, 1.00) && near(w[0].im_high, INF), "upper band is [1.00, +inf)");
    check(near(w[0].im_low, w[1].im_high),
          "adjacent bands share the boundary exactly (the band is half-open)");
  }

  // ------------------------------------------------- order independence (S6)
  {
    // The defect: pairing each window in turn with its own nearest unclaimed
    // position. Window 0 (centre 500) takes the position at 502 because it is
    // nearer than 503 -- but 502 names window 1 exactly, and window 1 is then
    // left with the position that does not name it. The two windows end up
    // with each other's mobility, i.e. with each other's band, which rejects
    // every real peak in both.
    //
    // Taking the globally smallest distance first pairs (window 1, 502) and
    // leaves window 0 with 503, which is the only self-consistent answer.
    std::vector<ODIA::IsolationWindow> w{window(500.0, 6.0), window(502.0, 6.0)};
    const std::vector<ODIA::MobilityPosition> ions{{502.0, 0.80}, {503.0, 1.20}};
    checkResult(ODIA::deriveMobilityBands(ions, w), ODIA::MobilityBandResult::Derived,
                "overlapping centres, positions that name them");
    // Window 1 took 502 (1/K0 0.80) so it is the LOWER band; window 0 took 503.
    check(near(w[1].im_high, 1.00) && near(w[1].im_low, -INF),
          "the window a position names gets that position");
    check(near(w[0].im_low, 1.00) && near(w[0].im_high, INF),
          "and the other window gets the other one");

    // The same input with the windows listed the other way round must give the
    // same bands. First-come matching does not.
    std::vector<ODIA::IsolationWindow> r{window(502.0, 6.0), window(500.0, 6.0)};
    checkResult(ODIA::deriveMobilityBands(ions, r), ODIA::MobilityBandResult::Derived,
                "the same frame with the windows listed in the other order");
    check(near(r[0].im_low, w[1].im_low) && near(r[0].im_high, w[1].im_high) &&
          near(r[1].im_low, w[0].im_low) && near(r[1].im_high, w[0].im_high),
          "the attribution does not depend on the order the windows are listed in");
  }

  // ------------------------------------------------------- a spare position
  {
    // Positions are collected over EVERY precursor of the spectrum, because
    // the file mis-attaches them all to precursor 0 -- but a window is only
    // built for a precursor that states one. A precursor with no isolation
    // window therefore leaves more positions than windows, and demanding
    // exactly as many of one as of the other abandoned the whole derivation
    // for the frame, silently.
    std::vector<ODIA::IsolationWindow> w{window(500.0, 25.0), window(525.0, 25.0)};
    const std::vector<ODIA::MobilityPosition> ions{
      {500.0, 1.10}, {525.0, 0.90}, {900.0, 0.70}};
    checkResult(ODIA::deriveMobilityBands(ions, w), ODIA::MobilityBandResult::Derived,
                "a position naming no window does not abandon the derivation");
    check(near(w[1].im_high, 1.00) && near(w[0].im_low, 1.00),
          "and the bands are the ones the two naming positions imply");
  }

  // --------------------------------------------------------------- refusals
  {
    std::vector<ODIA::IsolationWindow> w{window(500.0, 25.0), window(525.0, 25.0)};
    const std::vector<ODIA::MobilityPosition> one{{500.0, 1.10}};
    checkResult(ODIA::deriveMobilityBands(one, w), ODIA::MobilityBandResult::TooFewIons,
                "fewer positions than windows");
    check(!std::isfinite(w[0].im_low) && !std::isfinite(w[0].im_high),
          "a refusal leaves the windows untouched");
  }
  {
    // A position 40 Th from a 25 Th window's centre is not naming it. Half the
    // width plus 1 Th is the limit, and beyond it the attribution is a guess
    // about a guess.
    std::vector<ODIA::IsolationWindow> w{window(500.0, 25.0), window(525.0, 25.0)};
    const std::vector<ODIA::MobilityPosition> far{{500.0, 1.10}, {565.0, 0.90}};
    checkResult(ODIA::deriveMobilityBands(far, w), ODIA::MobilityBandResult::Unattributable,
                "a position too far from any centre");
  }
  {
    // Two positions equidistant from the same centre: whichever is chosen, the
    // other window inherits the mobility of the wrong one. Refuse.
    std::vector<ODIA::IsolationWindow> w{window(500.0, 25.0), window(510.0, 25.0)};
    const std::vector<ODIA::MobilityPosition> tie{{495.0, 1.10}, {505.0, 0.90}};
    checkResult(ODIA::deriveMobilityBands(tie, w), ODIA::MobilityBandResult::Ambiguous,
                "two equally good attributions");
  }
  {
    std::vector<ODIA::IsolationWindow> w{window(500.0, 25.0)};
    const std::vector<ODIA::MobilityPosition> ions{{500.0, 1.10}};
    checkResult(ODIA::deriveMobilityBands(ions, w), ODIA::MobilityBandResult::NotNeeded,
                "one window needs no split");
  }

  // ------------------------------------------------------- a stated band wins
  {
    std::vector<ODIA::IsolationWindow> w{window(500.0, 25.0), window(525.0, 25.0)};
    w[0].im_low = 0.60; w[0].im_high = 0.75;
    const std::vector<ODIA::MobilityPosition> ions{{500.0, 1.10}, {525.0, 0.90}};
    checkResult(ODIA::deriveMobilityBands(ions, w), ODIA::MobilityBandResult::Derived,
                "a frame where one window states its own band");
    check(near(w[0].im_low, 0.60) && near(w[0].im_high, 0.75),
          "what the file states is never overridden by what we derive");
  }

  std::printf(failures ? "FAILED (%d)\n" : "ok (%d)\n", failures);
  return failures ? 1 : 0;
}
