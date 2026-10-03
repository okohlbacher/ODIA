// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// -tof_calibration frame: rebuild m/z from a spectrum's own TOF calibration
// pair instead of the archive's run-wide chord.
//
// The coefficients are the real ones of PXD047793 run 009 (the archive's chord
// and the first frame's tof_c0/tof_c1), so the size of the shift checked here
// is the size measured against the Bruker SDK in shared/pxd/R0B_MASS.md.

#include <odia/TofCalibration.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <stdexcept>
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

  template <class F>
  bool throws(F f)
  {
    try { f(); }
    catch (const std::runtime_error&) { return true; }
    return false;
  }

  // PXD047793 009: archive chord, and spectrum 0's per-frame pair.
  const double A = 9.746500807982319;
  const double B = 7.886936881362772e-05;
  const double C0 = 9.746506817194948;
  const double C1 = 7.886927114780291e-05;

  double chord(double t) { const double r = A + B * t; return r * r; }
  double pair(double t) { const double r = C0 + C1 * t; return r * r; }
}

int main()
{
  // ---------------------------------------------------- the remap itself
  {
    // TOF indices across the whole acquired range of this instrument
    // (~95 to ~1,705 Th), plus every index of a dense stretch.
    std::vector<double> tofs;
    for (double t = 1000; t <= 400000; t += 997) { tofs.push_back(t); }
    for (double t = 200000; t < 200500; t += 1) { tofs.push_back(t); }
    std::vector<double> mz;
    for (double t : tofs) { mz.push_back(chord(t)); }

    ODIA::TofRemapStats st;
    ODIA::remapTofMz(mz, A, B, C0, C1, &st);
    bool exact = true;
    for (std::size_t i = 0; i < mz.size(); ++i) { exact = exact && mz[i] == pair(tofs[i]); }
    check(exact, "every value equals (c0 + c1*tof)^2 at its integer TOF, bit for bit");
    check(st.peaks == tofs.size(), "stats count every peak");
    check(st.max_index_residual < 1e-6,
          "inverted TOF lands on the integer within 1e-6 (max " + std::to_string(st.max_index_residual) + ")");

    // The shift has the measured size and sign: chord - SDK was +0.47 ppm in
    // the lowest m/z decile and +1.38 at the top (R0B_MASS), so the pair moves
    // m/z DOWN by about that much.
    const double lo = ODIA::tofShiftPpm(400.0, A, B, C0, C1);
    const double hi = ODIA::tofShiftPpm(1200.0, A, B, C0, C1);
    check(lo < -0.3 && lo > -1.0, "shift at 400 Th is -0.3..-1.0 ppm (" + std::to_string(lo) + ")");
    check(hi < -1.0 && hi > -1.8, "shift at 1200 Th is -1.0..-1.8 ppm (" + std::to_string(hi) + ")");
    check(hi < lo, "the shift grows with m/z");
  }

  // ------------------------------------------- identity pair is (near) identity
  {
    std::vector<double> mz{chord(12345), chord(234567), chord(399999)};
    const auto before = mz;
    ODIA::remapTofMz(mz, A, B, A, B);
    bool same = true;
    for (std::size_t i = 0; i < mz.size(); ++i) { same = same && mz[i] == before[i]; }
    check(same, "the chord as its own pair reproduces the input exactly");
  }

  // -------------------------------------------------------- refusals
  {
    // Half an index off: an m/z the chord did not produce.
    std::vector<double> bad{chord(100000.5)};
    check(throws([&] { ODIA::remapTofMz(bad, A, B, C0, C1); }),
          "an m/z between two TOF indices is refused, not silently corrected");

    // A mis-stated chord: values come from (A, B) but the remap is told A*1.0001.
    std::vector<double> wrong{chord(150000), chord(250000)};
    check(throws([&] { ODIA::remapTofMz(wrong, A * 1.0001, B, C0, C1); }),
          "a wrong run-wide chord is refused");

    std::vector<double> ok{chord(150000)};
    check(throws([&] { ODIA::remapTofMz(ok, A, B, std::nan(""), C1); }), "a NaN coefficient is refused");
    check(throws([&] { ODIA::remapTofMz(ok, A, B, C0, 0.0); }), "a zero slope is refused");
    check(throws([&] { ODIA::remapTofMz(ok, A, 0.0, C0, C1); }), "a zero chord slope is refused");

    std::vector<double> empty;
    ODIA::remapTofMz(empty, A, B, C0, C1);
    check(empty.empty(), "an empty spectrum stays empty");
  }

  // ------------------------------------------------ the coefficient table
  {
    // Rows out of order; spectrum 2 NULL (stays on the chord); spectrum 3 has
    // no row at all; spectrum 1 repeated with the same pair.
    const std::vector<std::uint64_t> idx{1, 0, 2, 1};
    const std::vector<double> c0{C0, A, 0.0, C0};
    const std::vector<double> c1{C1, B, 0.0, C1};
    const std::vector<bool> n0{false, false, true, false};
    const std::vector<bool> n1{false, false, true, false};
    const auto t = ODIA::TofCoefficientTable::fromColumns(idx, c0, c1, n0, n1, 4);
    check(t.size() == 4, "table is sized by the run, not the rows");
    check(t.has(0) && t.c0(0) == A && t.c1(0) == B, "spectrum 0 keeps its own row");
    check(t.has(1) && t.c0(1) == C0 && t.c1(1) == C1, "a repeated identical row is accepted");
    check(!t.has(2), "a NULL pair leaves the spectrum on the chord");
    check(!t.has(3), "a spectrum without a row stays on the chord");
    check(!t.has(99), "an index past the run has no pair");
    check(t.withPair() == 2 && t.onChord() == 2, "counts: 2 with a pair, 2 on the chord");

    check(throws([&] {
            ODIA::TofCoefficientTable::fromColumns({5}, {C0}, {C1}, {false}, {false}, 4);
          }),
          "an index outside the run is refused");
    check(throws([&] {
            ODIA::TofCoefficientTable::fromColumns({0, 0}, {C0, A}, {C1, B}, {false, false},
                                                   {false, false}, 1);
          }),
          "two different pairs for one spectrum are refused");
    check(throws([&] {
            ODIA::TofCoefficientTable::fromColumns({0, 0}, {C0, 0.0}, {C1, 0.0}, {false, true},
                                                   {false, true}, 1);
          }),
          "a pair and a NULL for one spectrum are refused");
    check(throws([&] {
            ODIA::TofCoefficientTable::fromColumns({0}, {C0}, {0.0}, {false}, {true}, 1);
          }),
          "half a pair is refused");
    check(throws([&] {
            ODIA::TofCoefficientTable::fromColumns({0}, {C0}, {-C1}, {false}, {false}, 1);
          }),
          "a negative slope is refused");
  }

  // --------------------------------------------- archive without the columns
  {
    check(throws([] { ODIA::TofCoefficientTable::fromArchive("/nonexistent/run.mzpeak", 1); }),
          "a missing archive is refused");
  }

  std::printf("%d failure(s)\n", failures);
  return failures == 0 ? 0 : 1;
}
