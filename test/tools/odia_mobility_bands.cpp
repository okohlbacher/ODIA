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

  void checkPairing(ODIA::IonPairing got, ODIA::IonPairing want, const std::string& what)
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

  // ===================================================== -im_bands_from_params
  // A stock mzpeak-convert 0.12.5 diaPASEF frame, as ODIA's reader hands it
  // (shared/pxd/imbands, raw dump of PXD047793 run 009): two precursors 763 and
  // 413 with NULL precursor_index, so BOTH selected ions land on precursor 0 in
  // file order (763 then 413) and precursor 1 gets none; the band is only in the
  // ions' MZP:1000006/7 parameters, stringified to 17 significant digits.
  const double NAN_ = std::numeric_limits<double>::quiet_NaN();
  {
    std::vector<std::size_t> ion_of;
    checkPairing(ODIA::pairIonsByPosition({2, 0}, {763.0, 413.0}, {763.0, 413.0}, ion_of),
                 ODIA::IonPairing::ByPosition, "stock 0.12.5 frame: both ions on precursor 0");
    check(ion_of.size() == 2 && ion_of[0] == 0 && ion_of[1] == 1,
          "the k-th ion of the frame is the k-th precursor's (763 -> 763, 413 -> 413)");

    // Legacy reading of the same frame gives precursor 0 the LAST ion's band
    // (413's) and precursor 1 nothing -- the wrong band, which R0C measured as
    // worse than none. Positional pairing is what removes that.
    checkPairing(ODIA::pairIonsByPosition({1, 1}, {763.0, 413.0}, {763.0, 413.0}, ion_of),
                 ODIA::IonPairing::AsAttached, "patched converter: one ion per precursor");
    check(ion_of.size() == 2 && ion_of[0] == 0 && ion_of[1] == 1, "as attached keeps the order");

    checkPairing(ODIA::pairIonsByPosition({1}, {763.0}, {763.0}, ion_of),
                 ODIA::IonPairing::AsAttached, "a single-window spectrum");
  }
  {
    // Refusals: never guess.
    std::vector<std::size_t> ion_of{7};
    checkPairing(ODIA::pairIonsByPosition({3, 0}, {763.0, 413.0}, {763.0, 413.0, 500.0}, ion_of),
                 ODIA::IonPairing::CountMismatch, "three ions for two precursors");
    check(ion_of.empty(), "a refusal leaves no pairing behind");
    checkPairing(ODIA::pairIonsByPosition({1, 0}, {763.0, 413.0}, {763.0}, ion_of),
                 ODIA::IonPairing::CountMismatch, "one ion for two precursors");
    checkPairing(ODIA::pairIonsByPosition({2, 0}, {763.0, 413.0}, {763.0}, ion_of),
                 ODIA::IonPairing::CountMismatch, "attached count disagrees with the ion list");
    // The writer listed the ions in the OTHER order: position would hand each
    // window its neighbour's band. The m/z check is what catches it.
    checkPairing(ODIA::pairIonsByPosition({2, 0}, {763.0, 413.0}, {413.0, 763.0}, ion_of),
                 ODIA::IonPairing::MzMismatch, "ions in the reverse order are refused");
    check(ion_of.empty(), "an m/z refusal leaves no pairing behind");
    checkPairing(ODIA::pairIonsByPosition({2, 0}, {763.0, 413.0}, {763.0, 413.11}, ion_of),
                 ODIA::IonPairing::MzMismatch, "0.11 Th off is not the window");
    checkPairing(ODIA::pairIonsByPosition({2, 0}, {763.0, 413.0}, {763.0, 413.0001}, ion_of),
                 ODIA::IonPairing::ByPosition, "float noise between the columns is the window");
    // Unknown m/z on either side cannot contradict; position stands.
    checkPairing(ODIA::pairIonsByPosition({2, 0}, {763.0, NAN_}, {NAN_, 413.0}, ion_of),
                 ODIA::IonPairing::ByPosition, "an absent m/z does not refuse");
  }
  {
    // The parameters themselves.
    const std::vector<ODIA::CvValue> stock_763{
      {"MZP:1000006", "0.89155395131847848"}, {"MZP:1000007", "1.6375174050443457"}};
    double lo = -1.0, hi = -1.0;
    check(ODIA::mobilityBandFromParameters(stock_763, lo, hi), "the stock 763 band parses");
    // The patched converter's Float64 column holds 0.8915539513184785 (R0C);
    // 17 significant digits round-trip to the same double, bit for bit.
    check(lo == 0.8915539513184785 && hi == 1.6375174050443457,
          "parsed band equals the column's double exactly");

    const std::vector<ODIA::CvValue> reversed{{"MZP:1000006", "1.6375"}, {"MZP:1000007", "0.8915"}};
    check(ODIA::mobilityBandFromParameters(reversed, lo, hi) && lo == 0.8915 && hi == 1.6375,
          "limits in scan order (lower > upper) are swapped, not rejected");

    const std::vector<ODIA::CvValue> other_params{
      {"MS:1000045", "45.7"}, {"MZP:1000007", "0.9"}, {"", ""}, {"MZP:1000006", "0.6"}};
    check(ODIA::mobilityBandFromParameters(other_params, lo, hi) && lo == 0.6 && hi == 0.9,
          "unrelated parameters and order are ignored");

    lo = hi = -1.0;
    check(!ODIA::mobilityBandFromParameters({{"MZP:1000006", "0.6"}}, lo, hi) && lo == -1.0 &&
            hi == -1.0,
          "one limit alone is not a band, and nothing is written");
    check(!ODIA::mobilityBandFromParameters({}, lo, hi), "no parameters, no band");
    check(!ODIA::mobilityBandFromParameters({{"MZP:1000006", "0.6abc"}, {"MZP:1000007", "0.9"}}, lo,
                                            hi),
          "a numeric prefix is not a number");
    check(!ODIA::mobilityBandFromParameters({{"MZP:1000006", "nan"}, {"MZP:1000007", "0.9"}}, lo, hi),
          "NaN is not a bound");
    check(!ODIA::mobilityBandFromParameters({{"MZP:1000006", "0.6"}, {"MZP:1000007", "inf"}}, lo, hi),
          "an infinite bound is not stated");
    check(!ODIA::mobilityBandFromParameters({{"MZP:1000006", ""}, {"MZP:1000007", "0.9"}}, lo, hi),
          "an empty value is not stated");
    check(!ODIA::mobilityBandFromParameters({{"MZP:1000006", "0.9"}, {"MZP:1000007", "0.9"}}, lo, hi),
          "a degenerate band is not a band");
    check(!ODIA::mobilityBandFromParameters(
            {{"MZP:1000006", "0.6"}, {"MZP:1000007", "0.9"}, {"MZP:1000006", "0.7"}}, lo, hi),
          "two different lower limits contradict each other");
    check(ODIA::mobilityBandFromParameters(
            {{"MZP:1000006", "0.6"}, {"MZP:1000007", "0.9"}, {"MZP:1000006", "0.6"}}, lo, hi),
          "a repeated identical limit is harmless");
  }

  std::printf(failures ? "FAILED (%d)\n" : "ok (%d)\n", failures);
  return failures ? 1 : 0;
}
