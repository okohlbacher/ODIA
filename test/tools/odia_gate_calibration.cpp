// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Gate C hash calibration (ROUND1_PLAN D6) against a scripted run.
//
// Arrival mode calibrates tau on the first gate_calibration_n decoys to reach
// the gate and admits everything before them unconditionally, so its decisions
// depend on hand-over order -- which chunking and the batch grid change. Hash
// mode selects the sample from the pass PLAN by identity hash, measures it in a
// pre-pass and freezes tau, so:
//
//   plan       the plan's pointCount(0) per precursor equals what the
//              extraction actually hands over (the selection's eligibility test)
//   select     the selection is the n smallest (hash, index) eligible decoys,
//              by an independent brute force, and the same IDENTITIES when the
//              library is presented in reverse order
//   invariant  tau and every gate decision are identical across chunk caps,
//              decode blocks and thread counts; every decision is made against
//              the frozen tau (no warm-up row), and admitted == (stat >= tau)
//              and when the library (hence hand-over order) is reversed
//   warmup     arrival mode on the same data DOES have warm-up rows (the
//              defect the mode removes is present in the control)
//   short      fewer eligible decoys than requested: the gate does not arm and
//              admits every precursor
//   pin        the identity hash is pinned, so a silent change of the sample
//              definition fails here
//
// Usage: odia_gate_calibration

#include <odia/ChromatogramExtractor.h>
#include <odia/Library.h>
#include <odia/PeakGroupScorer.h>
#include <odia/SpectrumSource.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace
{
  const float NA = std::numeric_limits<float>::quiet_NaN();
  int failures = 0;

  void check(bool ok, const std::string& what)
  {
    std::printf("%s %s\n", ok ? "  ok  " : "  FAIL", what.c_str());
    if (!ok) { ++failures; }
  }

  class ScriptedRun : public ODIA::SpectrumSource
  {
  public:
    std::size_t addWindow(double lo, double hi)
    {
      ODIA::IsolationWindow w;
      w.mz_low = lo; w.mz_high = hi;
      w.im_low = -std::numeric_limits<double>::infinity();
      w.im_high = std::numeric_limits<double>::infinity();
      windows_.push_back(w);
      return windows_.size() - 1;
    }
    std::size_t addSpectrum(std::size_t window, double rt)
    {
      ODIA::SpectrumInfo s;
      s.index = info_.size();
      s.retention_time = rt;
      s.window = windows_[window];
      info_.push_back(s);
      peaks_.emplace_back();
      return info_.size() - 1;
    }
    void addPeak(std::size_t spectrum, double mz, float intensity)
    {
      peaks_[spectrum].mz.push_back(mz);
      peaks_[spectrum].intensity.push_back(intensity);
    }
    const std::vector<ODIA::SpectrumInfo>& spectra() const override { return info_; }
    const std::vector<ODIA::IsolationWindow>& windows() const override { return windows_; }
    std::string describe() const override { return "scripted gate run"; }
    void peaks(std::size_t begin, std::size_t end, std::vector<ODIA::SpectrumPeaks>& out) override
    { out.assign(peaks_.begin() + begin, peaks_.begin() + end); }

  private:
    std::vector<ODIA::SpectrumInfo> info_;
    std::vector<ODIA::IsolationWindow> windows_;
    std::vector<ODIA::SpectrumPeaks> peaks_;
  };

  constexpr std::size_t N_PREC = 360;     ///< half decoys
  constexpr std::size_t N_TRANS = 4;
  constexpr std::size_t N_SPEC = 240;
  constexpr double RT0 = 100.0, DRT = 1.37;   ///< not exact in float, on purpose
  constexpr double HALF = 25.0;               ///< extraction half-window, s

  double productMz(std::size_t i, std::size_t k) { return 150.0 + 3.0 * double(i) + 0.7 * double(k); }
  double centreRt(std::size_t i) { return RT0 + 5.0 + double(i % 120) * (DRT * N_SPEC - 10.0) / 120.0; }
  std::string sequenceOf(std::size_t i)
  {
    // Decoys reuse a target's letters reversed, as real decoy libraries do, so
    // the decoy flag is part of what keeps identities distinct.
    std::string s = "PEPTIDE" + std::to_string(i / 2) + "K";
    if (i % 2 == 1) { std::reverse(s.begin(), s.end() - 1); }
    return s;
  }

  /// Precursor order: identity j is library row order[j].
  ODIA::Library makeLibrary(bool reversed)
  {
    ODIA::Library lib;
    auto& p = lib.precursors();
    auto& t = lib.transitions();
    for (std::size_t r = 0; r < N_PREC; ++r)
    {
      const std::size_t i = reversed ? N_PREC - 1 - r : r;   // identity
      p.mz.push_back(ODIA::toFixed(500.0));
      p.irt.push_back(static_cast<float>(centreRt(i)));
      p.im.push_back(NA);
      p.ccs.push_back(NA);
      p.charge.push_back(static_cast<std::uint8_t>(2 + (i % 3 == 0)));
      p.decoy.push_back(static_cast<std::uint8_t>(i % 2));
      p.modified_sequence.push_back(lib.strings().intern(sequenceOf(i)));
      p.protein_group.push_back(0);
      p.transition_begin.push_back(static_cast<std::uint32_t>(t.product_mz.size()));
      p.transition_count.push_back(0);
      // A few precursors whose transition 0 is unrepresentable: transition 0
      // decides pointCount(0), so the plan must say 0 for them.
      for (std::size_t k = 0; k < N_TRANS; ++k)
      {
        const bool invalid0 = (k == 0 && i % 37 == 5);
        t.product_mz.push_back(invalid0 ? ODIA::toFixed(0.0) : ODIA::toFixed(productMz(i, k)));
        t.library_intensity.push_back(1.0f);
        t.type.push_back(ODIA::FragmentType::Y);
        t.ordinal.push_back(static_cast<std::uint8_t>(4 + k));
        t.charge.push_back(1);
        t.loss.push_back(ODIA::LossType::None);
        ++p.transition_count.back();
      }
    }
    return lib;
  }

  std::uint32_t lcg(std::uint32_t& s) { s = s * 1664525u + 1013904223u; return s >> 8; }

  ScriptedRun makeRun()
  {
    ScriptedRun run;
    const std::size_t w = run.addWindow(400.0, 600.0);
    std::uint32_t seed = 12345u;
    for (std::size_t s = 0; s < N_SPEC; ++s)
    {
      const double rt = RT0 + DRT * double(s);
      const std::size_t si = run.addSpectrum(w, rt);
      for (std::size_t i = 0; i < N_PREC; ++i)
      {
        // Co-eluting signal on every third identity (targets and decoys), of
        // graded height, and a decoy-like spike pattern elsewhere.
        const bool has_peak = (i % 3 == 0);
        const double d = (rt - centreRt(i)) / 3.0;
        const double peak = has_peak ? (5.0 + double(i % 17) * 3.0) * std::exp(-0.5 * d * d) : 0.0;
        for (std::size_t k = 0; k < N_TRANS; ++k)
        {
          const float noise = 1.0f + float(lcg(seed) % 1000) / 250.0f;
          run.addPeak(si, productMz(i, k), static_cast<float>(noise + peak * (1.0 + 0.2 * k)));
        }
      }
    }
    return run;
  }

  ODIA::ChromatogramExtractor::Options extractOptions(unsigned threads, std::size_t cap,
                                                      std::size_t block)
  {
    ODIA::ChromatogramExtractor::Options o;
    o.threads = threads;
    o.progress_every = 0;
    o.irt_slope = 1.0;
    o.irt_intercept = 0.0;
    o.rt_window_seconds = HALF;
    o.max_live_precursors = cap;
    o.decode_block = block;
    return o;
  }

  ODIA::PeakGroupScorer::Options scorerOptions(const std::string& mode, std::size_t n,
                                               const std::string& log_path, std::uint8_t* reasons)
  {
    ODIA::PeakGroupScorer::Options o;
    o.gate_alpha = 0.2;
    o.gate_calibration = mode;
    o.gate_calibration_n = n;
    o.gate_log_path = log_path;
    o.terminal_reason = reasons;
    o.threads = 1;
    return o;
  }

  struct Decision { int decoy; double stat; double tau; int ready; int admitted; };

  std::map<std::uint32_t, Decision> readGateLog(const std::string& path)
  {
    std::map<std::uint32_t, Decision> out;
    std::ifstream in(path);
    std::string line;
    std::getline(in, line);   // header
    while (std::getline(in, line))
    {
      // strtod, not operator>>: the log prints a disarmed gate's tau as -inf,
      // which an istream refuses, silently zeroing every field after it.
      std::istringstream ss(line);
      std::string f[6];
      for (auto& x : f) { std::getline(ss, x, '\t'); }
      Decision d{};
      d.decoy = std::atoi(f[1].c_str());
      d.stat = std::strtod(f[2].c_str(), nullptr);
      d.tau = std::strtod(f[3].c_str(), nullptr);
      d.ready = std::atoi(f[4].c_str());
      d.admitted = std::atoi(f[5].c_str());
      out[static_cast<std::uint32_t>(std::strtoul(f[0].c_str(), nullptr, 10))] = d;
    }
    return out;
  }

  struct Arm
  {
    ODIA::PeakGroupScorer::GateCalibrationReport rep;
    std::map<std::uint32_t, Decision> log;
    std::vector<std::uint8_t> reasons;
  };

  /// One scoring extraction exactly as the tool drives it: calibrate (hash
  /// mode only), then extract into the same sink. Session::finish is not
  /// needed -- the gate's decisions are what is under test.
  Arm runArm(const ODIA::Library& lib, ScriptedRun& run, const std::string& mode, std::size_t n,
             unsigned threads, std::size_t cap, std::size_t block, const std::string& tag)
  {
    Arm a;
    a.reasons.assign(lib.precursorCount(), 0);
    const std::string path = "odia_gate_calibration_" + tag + ".tsv";
    std::remove(path.c_str());
    const auto so = scorerOptions(mode, n, path, a.reasons.data());
    ODIA::PeakGroupScorer::Sink sink(lib, so);
    const auto eo = extractOptions(threads, cap, block);
    if (ODIA::PeakGroupScorer::gateUsesHashCalibration(sink.options()))
    { a.rep = ODIA::PeakGroupScorer::calibrateGateByHash(lib, run, eo, sink); }
    ODIA::ChromatogramExtractor::Stats st;
    ODIA::ChromatogramExtractor::extract(lib, run, eo, sink, &st);
    std::printf("    arm %-14s chunks %zu, decoded %zu, tau %.6g\n", tag.c_str(), st.chunks,
                st.spectra_decoded, a.rep.tau);
    // The log file stays open in the Session's gate null until the Sink dies.
    return a;
  }

  Arm finishArm(Arm a, const std::string& tag)
  {
    a.log = readGateLog("odia_gate_calibration_" + tag + ".tsv");
    return a;
  }

  class PointRecorder final : public ODIA::ChromatogramSink
  {
  public:
    explicit PointRecorder(std::size_t n) : points(n, 0) {}
    void accept(const ODIA::PrecursorChromatogram& c) override
    { points[c.precursor] = c.transition_count ? c.pointCount(0) : 0; }
    std::vector<std::uint32_t> points;
  };
}

int main()
{
  const ODIA::Library lib = makeLibrary(false);
  ScriptedRun run = makeRun();

  // ---------------------------------------------------------------- plan
  std::printf("plan\n");
  std::vector<std::uint32_t> plan;
  {
    auto eo = extractOptions(1, 0, 256);
    eo.plan_points = &plan;
    eo.plan_only = true;
    ODIA::NullChromatogramSink none;
    ODIA::ChromatogramExtractor::extract(lib, run, eo, none, nullptr);
  }
  {
    PointRecorder rec(lib.precursorCount());
    ODIA::ChromatogramExtractor::extract(lib, run, extractOptions(1, 3, 1), rec, nullptr);
    check(plan == rec.points, "plan pointCount(0) == handed-over pointCount(0) for every precursor");
    std::size_t zero = 0;
    for (std::size_t i = 0; i < plan.size(); ++i) { zero += plan[i] == 0; }
    check(zero > 0, "the plan has unassigned / invalid-transition-0 precursors (" + std::to_string(zero) + ")");
  }

  // ---------------------------------------------------------------- select
  std::printf("select\n");
  const std::size_t N = 40;
  {
    const auto so = scorerOptions("hash", N, "", nullptr);
    std::size_t eligible = 0;
    const auto sel = ODIA::PeakGroupScorer::selectGateCalibrationDecoys(lib, plan, so, &eligible);
    // Independent brute force over the same eligibility rule.
    const auto& p = lib.precursors();
    std::vector<std::pair<std::uint64_t, std::uint32_t>> all;
    for (std::size_t i = 0; i < lib.precursorCount(); ++i)
    {
      if (!p.decoy[i] || plan[i] < 3) { continue; }
      all.emplace_back(ODIA::PeakGroupScorer::gateIdentityHash(
                         sequenceOf(i), int(p.charge[i]), true, so.gate_calibration_seed),
                       std::uint32_t(i));
    }
    std::sort(all.begin(), all.end());
    std::vector<std::uint32_t> want;
    for (std::size_t j = 0; j < N && j < all.size(); ++j) { want.push_back(all[j].second); }
    std::sort(want.begin(), want.end());
    check(eligible == all.size() && eligible > N, "eligible count " + std::to_string(eligible));
    check(sel == want, "selection == the N smallest identity hashes (brute force)");

    // Same identities when the library is presented in reverse.
    const ODIA::Library rev = makeLibrary(true);
    std::vector<std::uint32_t> rplan;
    auto eo = extractOptions(1, 0, 256);
    eo.plan_points = &rplan;
    eo.plan_only = true;
    ODIA::NullChromatogramSink none;
    ODIA::ChromatogramExtractor::extract(rev, run, eo, none, nullptr);
    const auto rsel = ODIA::PeakGroupScorer::selectGateCalibrationDecoys(rev, rplan, so, nullptr);
    std::set<std::string> a, b;
    for (auto i : sel) { a.insert(std::string(lib.strings().get(lib.precursors().modified_sequence[i]))); }
    for (auto i : rsel) { b.insert(std::string(rev.strings().get(rev.precursors().modified_sequence[i]))); }
    check(a == b && a.size() == N, "the selected identities do not depend on library order");

    // A different seed draws a different sample.
    auto so2 = so; so2.gate_calibration_seed = 7;
    check(ODIA::PeakGroupScorer::selectGateCalibrationDecoys(lib, plan, so2, nullptr) != sel,
          "a different seed selects a different sample");
  }

  // ---------------------------------------------------------------- invariant
  std::printf("invariant\n");
  struct Cfg { unsigned threads; std::size_t cap, block; const char* tag; };
  const Cfg cfgs[] = {{1, 0, 256, "h_c0"}, {1, 7, 256, "h_c7"}, {2, 3, 1, "h_c3_b1"},
                      {2, 11, 129, "h_c11_b129"}, {1, 2, 200, "h_c2_b200"}};
  std::vector<Arm> arms;
  for (const auto& c : cfgs)
  {
    Arm a;
    { a = runArm(lib, run, "hash", N, c.threads, c.cap, c.block, c.tag); }
    arms.push_back(finishArm(std::move(a), c.tag));
  }
  const Arm& ref = arms.front();
  check(ref.rep.armed && ref.rep.selected == N && ref.rep.measured == N, "hash gate armed on N measured decoys");
  bool all_ready = true, rule = true, some_reject = false, some_admit = false;
  for (const auto& [pr, d] : ref.log)
  {
    all_ready = all_ready && d.ready == 1;
    // The log prints %.6g, so compare against tau with that rounding allowed.
    const double tol = 1e-5 * std::max(1.0, std::abs(ref.rep.tau));
    rule = rule && (d.admitted == 1 ? d.stat >= ref.rep.tau - tol : d.stat <= ref.rep.tau + tol);
    some_reject = some_reject || d.admitted == 0;
    some_admit = some_admit || d.admitted == 1;
  }
  check(!ref.log.empty() && all_ready, "every gate decision made against the frozen tau (no warm-up)");
  check(rule, "admitted == (statistic >= tau) for every precursor");
  check(some_reject && some_admit, "the gate both admits and rejects on this run");
  {
    // The pre-pass measures the SAME statistic the production pass gates on:
    // tau recomputed from the production log's values for the selected decoys.
    const auto sel = ODIA::PeakGroupScorer::selectGateCalibrationDecoys(
      lib, plan, scorerOptions("hash", N, "", nullptr), nullptr);
    std::vector<double> v;
    for (auto i : sel) { const auto it = ref.log.find(i); if (it != ref.log.end()) { v.push_back(it->second.stat); } }
    std::sort(v.begin(), v.end());
    const std::size_t k = std::min(v.size() - 1, std::size_t((1.0 - 0.2) * double(v.size())));
    check(v.size() == N && std::abs(v[k] - ref.rep.tau) <= 1e-5 * std::max(1.0, std::abs(ref.rep.tau)),
          "tau == quantile of the production-pass statistics of the selected decoys");
  }
  for (std::size_t j = 1; j < arms.size(); ++j)
  {
    const Arm& a = arms[j];
    bool same_log = a.log.size() == ref.log.size();
    if (same_log)
    {
      for (const auto& [pr, d] : ref.log)
      {
        const auto it = a.log.find(pr);
        if (it == a.log.end() || it->second.admitted != d.admitted || it->second.stat != d.stat)
        { same_log = false; break; }
      }
    }
    check(a.rep.tau == ref.rep.tau, std::string(cfgs[j].tag) + ": tau identical");
    check(same_log, std::string(cfgs[j].tag) + ": every gate decision identical");
    check(a.reasons == ref.reasons, std::string(cfgs[j].tag) + ": every terminal reason identical");
  }

  // Hand-over ORDER itself, not only chunking: the reversed library hands
  // co-ending precursors over in the opposite order. Decisions keyed by identity.
  {
    const ODIA::Library rev = makeLibrary(true);
    Arm a;
    { a = runArm(rev, run, "hash", N, 1, 0, 256, "h_rev"); }
    a = finishArm(std::move(a), "h_rev");
    bool same = a.log.size() == ref.log.size();
    for (const auto& [pr, d] : a.log)
    {
      const auto it = ref.log.find(static_cast<std::uint32_t>(N_PREC - 1 - pr));
      if (it == ref.log.end() || it->second.admitted != d.admitted || it->second.stat != d.stat)
      { same = false; break; }
    }
    check(a.rep.tau == ref.rep.tau, "reversed library: tau identical");
    check(same, "reversed library: every gate decision identical by identity");
  }

  // ---------------------------------------------------------------- warmup (control)
  std::printf("warmup\n");
  {
    Arm a;
    { a = runArm(lib, run, "arrival", N, 1, 0, 256, "arrival_c0"); }
    a = finishArm(std::move(a), "arrival_c0");
    std::size_t unready = 0;
    for (const auto& [pr, d] : a.log) { unready += d.ready == 0; }
    check(unready > 0, "arrival mode admits a warm-up tranche unconditionally (" +
                       std::to_string(unready) + " rows) -- the defect hash mode removes");
  }

  // ---------------------------------------------------------------- short sample
  std::printf("short\n");
  {
    std::size_t eligible = 0;
    ODIA::PeakGroupScorer::selectGateCalibrationDecoys(lib, plan, scorerOptions("hash", N, "", nullptr), &eligible);
    Arm a;
    { a = runArm(lib, run, "hash", eligible + 1, 1, 0, 256, "h_short"); }
    a = finishArm(std::move(a), "h_short");
    bool all_admitted = !a.log.empty();
    for (const auto& [pr, d] : a.log) { all_admitted = all_admitted && d.admitted == 1; }
    check(!a.rep.armed && a.rep.selected == 0 && std::isinf(a.rep.tau) && a.rep.tau < 0,
          "insufficient sample: the gate does not arm");
    check(all_admitted, "insufficient sample: every precursor admitted");
    std::printf("    %s\n", a.rep.describe().c_str());
  }
  std::printf("    %s\n", ref.rep.describe().c_str());

  // ---------------------------------------------------------------- pin
  std::printf("pin\n");
  {
    const std::uint64_t h = ODIA::PeakGroupScorer::gateIdentityHash("PEPTIDEK", 2, true, 0x0D1A5EEDULL);
    std::printf("    hash(PEPTIDEK, 2, decoy, default seed) = 0x%016llx\n", (unsigned long long)h);
    check(h == 0xa0a26262425e84d3ULL, "identity hash pinned (a change redefines every hash-mode sample)");
  }

  std::printf(failures ? "FAILED (%d)\n" : "ok (%d)\n", failures);
  return failures ? 1 : 0;
}
