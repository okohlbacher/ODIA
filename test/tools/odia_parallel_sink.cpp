// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

/// -parallel_sink must be byte-identical to the serial sink.
///
/// The parallel sink runs each precursor's work on a pool and commits it on
/// the driver in the serial order. Everything that depends on WHICH precursors
/// came before is committed there: Gate C's decoy null (sampling, the arming
/// boundary, the gate log), group indices, the mass-anchor cap, counters and
/// terminal reasons. Each of those is exercised here on purpose:
///
///  * Gate C in quantile mode with a small calibration sample, so the null
///    arms INSIDE the first parallel batch and the precursors after the
///    boundary were speculated as admitted and must be discarded at commit
///    (the test requires that this actually happened);
///  * a mass-anchor cap far below what the run produces, so the first-N cut
///    falls in the middle of a batch;
///  * terminal reasons, including precursors with no transitions and with an
///    empty trace;
///  * sub-batches smaller than the caller's batches.
///
/// The serial reference is `Session::add` one precursor at a time; the arms
/// are `addBatch` at 1 (flag off), 8 and 48 threads. Everything the Result
/// carries is serialised to bytes and compared, after `finish()`, together
/// with the terminal-reason table and the Gate C log file.

#include <odia/Library.h>
#include <odia/PeakGroupScorer.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace
{
  int failures = 0;

  void check(bool ok, const std::string& what)
  {
    std::printf("  [%s] %s\n", ok ? "ok" : "FAIL", what.c_str());
    if (!ok) { ++failures; }
  }

  constexpr std::uint32_t N = 1200;
  constexpr std::uint32_t TRANSITIONS = 6;
  constexpr std::uint32_t CYCLES = 60;

  std::uint64_t mix(std::uint64_t z)
  {
    z += 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
  }
  double uniform(std::uint64_t key) { return double(mix(key) >> 11) / double(1ULL << 53); }

  /// Precursor i: decoy when odd; every 97th has no transitions.
  std::uint32_t transitionsOf(std::uint32_t i) { return i % 97 == 5 ? 0 : TRANSITIONS; }

  ODIA::Library makeLibrary()
  {
    ODIA::Library lib;
    auto& p = lib.precursors();
    auto& t = lib.transitions();
    std::uint32_t tb = 0;
    for (std::uint32_t i = 0; i < N; ++i)
    {
      p.mz.push_back(ODIA::toFixed(400.0 + 0.5 * i));
      p.irt.push_back(0.0f);
      p.im.push_back(std::numeric_limits<float>::quiet_NaN());
      p.ccs.push_back(std::numeric_limits<float>::quiet_NaN());
      p.charge.push_back(2);
      p.decoy.push_back(static_cast<std::uint8_t>(i % 2));
      p.modified_sequence.push_back(i);
      p.protein_group.push_back(0);
      p.transition_begin.push_back(tb);
      const std::uint32_t tc = transitionsOf(i);
      p.transition_count.push_back(tc);
      for (std::uint32_t k = 0; k < tc; ++k)
      {
        t.product_mz.push_back(ODIA::toFixed(300.0 + 50.0 * k + 0.01 * i));
        t.library_intensity.push_back(static_cast<float>(1.0 + 0.4 * k + 0.3 * uniform(i * 31 + k)));
        t.type.push_back(k % 2 ? ODIA::FragmentType::B : ODIA::FragmentType::Y);
        t.ordinal.push_back(static_cast<std::uint8_t>(4 + k));
        t.charge.push_back(1);
        t.loss.push_back(ODIA::LossType::None);
      }
      tb += tc;
    }
    return lib;
  }

  /// Flat storage the traces point into, like the extractor's blocks.
  struct Run
  {
    std::vector<float> rt;
    std::vector<float> points, ppm_num, ppm_den;
    std::vector<std::uint64_t> offset;
    std::vector<std::uint32_t> count;
    std::vector<ODIA::PrecursorChromatogram> traces;
  };

  Run makeRun(const ODIA::Library& lib)
  {
    Run r;
    r.rt.resize(CYCLES);
    for (std::uint32_t j = 0; j < CYCLES; ++j) { r.rt[j] = 100.0f + 1.5f * float(j); }
    const auto& p = lib.precursors();
    // Offsets/counts per transition, laid out first so pointers stay valid.
    std::vector<std::size_t> first(N, 0);
    for (std::uint32_t i = 0; i < N; ++i)
    {
      first[i] = r.offset.size();
      const bool empty = i % 89 == 7;   // an extracted-nothing precursor
      for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
      {
        r.offset.push_back(r.points.size());
        const std::uint32_t c = empty ? 0 : CYCLES;
        r.count.push_back(c);
        for (std::uint32_t j = 0; j < c; ++j)
        {
          const std::uint64_t key = (std::uint64_t(i) << 20) ^ (std::uint64_t(k) << 10) ^ j;
          double v = 50.0 * uniform(key);
          // Targets: two in three carry a co-eluting peak; decoys: one in five
          // carries a peak in a subset of fragments (interference).
          const bool target = p.decoy[i] == 0;
          const double apex = 10.0 + 40.0 * uniform(i * 7919ULL);
          const bool peak = target ? (i % 3 != 0) : (i % 5 == 0 && k < 3);
          if (peak)
          {
            const double d = double(j) - apex;
            v += (300.0 + 200.0 * k) * (0.5 + uniform(i)) * std::exp(-(d * d) / (2.0 * 2.2 * 2.2));
          }
          if (uniform(key ^ 0xABCDEFULL) < 0.1) { v = 0.0; }
          r.points.push_back(static_cast<float>(v));
          const double ppm = -3.0 + 4.0 * (uniform(key ^ 0x5555ULL) - 0.5);
          r.ppm_num.push_back(v > 0.0 ? static_cast<float>(v * ppm) : 0.0f);
          r.ppm_den.push_back(static_cast<float>(v));
        }
      }
    }
    for (std::uint32_t i = 0; i < N; ++i)
    {
      ODIA::PrecursorChromatogram c;
      c.precursor = i;
      c.transition_begin = p.transition_begin[i];
      c.transition_count = p.transition_count[i];
      const bool empty = i % 89 == 7;
      c.cycles = empty ? 0 : CYCLES;
      c.rt = r.rt.data();
      c.points = r.points.data();
      c.ppm_num = r.ppm_num.data();
      c.ppm_den = r.ppm_den.data();
      c.offset = r.offset.data() + first[i];
      c.count = r.count.data() + first[i];
      r.traces.push_back(c);
    }
    return r;
  }

  template <class T> void put(std::string& out, const T& v)
  {
    out.append(reinterpret_cast<const char*>(&v), sizeof(T));
  }

  /// Every byte the Result carries that a consumer could read.
  std::string serialise(const ODIA::PeakGroupScorer::Result& r)
  {
    std::string out;
    put(out, r.groups.size());
    for (const auto& g : r.groups)
    {
      put(out, g.precursor); put(out, g.apex_rt); put(out, g.left_rt); put(out, g.right_rt);
      put(out, g.apex_intensity); put(out, g.mass_ppm); put(out, g.mass_ppm_n);
      put(out, g.mass_ppm_spread); put(out, g.observed_im); put(out, g.im_spread);
      put(out, g.im_frags); put(out, g.mass_ppm_frags); put(out, g.dscore);
      put(out, g.qvalue); put(out, g.pep); put(out, g.decoy);
      put(out, g.sub_scores.size());
      for (const double s : g.sub_scores) { put(out, s); }
    }
    put(out, r.mass_anchors.size());
    for (const auto& a : r.mass_anchors)
    {
      put(out, a.group); put(out, a.residual.mz); put(out, a.residual.rt);
      put(out, a.residual.ppm); put(out, a.residual.intensity); put(out, a.residual.im);
    }
    put(out, r.fragvec.size());
    put(out, r.mass_anchors_dropped);
    put(out, r.precursors_without_candidate);
    put(out, r.target_groups); put(out, r.decoy_groups);
    put(out, r.identified_at_1pct);
    put(out, r.candidates_below_library_corr);
    put(out, r.picker_rejects);   // all size_t, no padding (static_assert in the scorer)
    put(out, r.ms1_census);
    return out;
  }

  std::string slurp(const std::string& path)
  {
    std::ifstream in(path, std::ios::binary);
    std::ostringstream s;
    s << in.rdbuf();
    return s.str();
  }

  struct Arm
  {
    std::string bytes, reasons, gate_log;
    std::size_t discarded = 0;
    std::size_t gate_rejected = 0;
    std::size_t anchors_dropped = 0;
  };

  Arm runArm(const ODIA::Library& lib, const Run& run, const std::string& tmp,
             bool batched, bool parallel, unsigned threads)
  {
    const std::string log = tmp + "/gate_" + std::to_string(batched) + "_" +
                            std::to_string(parallel) + "_" + std::to_string(threads) + ".tsv";
    std::remove(log.c_str());
    std::vector<std::uint8_t> reasons(N, 0);
    reasons[11] = static_cast<std::uint8_t>(ODIA::PeakGroupScorer::TerminalReason::NoWindowCoverage);

    ODIA::PeakGroupScorer::Options o;
    o.classifier = "lda";
    o.coelution_picking = true;
    o.max_candidates = 5;
    o.min_fragments_at_apex = 1;
    o.gate_alpha = 0.3;
    o.gate_calibration_n = 100;
    o.gate_log_path = log;
    o.collect_mass_anchors = true;
    o.max_mass_anchors = 700;
    o.terminal_reason = reasons.data();
    o.threads = threads;
    o.parallel_sink = parallel;
    o.sink_batch = 64;

    ODIA::PeakGroupScorer::Session session(lib, o);
    if (!batched)
    {
      for (const auto& t : run.traces) { session.add(t); }
    }
    else
    {
      // Caller batches of 250 (as the extractor hands them over), split by the
      // session into sub-batches of 64.
      for (std::size_t b = 0; b < run.traces.size(); b += 250)
      {
        session.addBatch(run.traces.data() + b, std::min<std::size_t>(250, run.traces.size() - b));
      }
    }
    Arm a;
    a.discarded = session.speculationDiscarded();
    const auto r = session.finish();
    a.bytes = serialise(r);
    a.reasons.assign(reasons.begin(), reasons.end());
    a.gate_log = slurp(log);
    a.gate_rejected = r.picker_rejects.gate_c[0] + r.picker_rejects.gate_c[1];
    a.anchors_dropped = r.mass_anchors_dropped;
    std::remove(log.c_str());
    return a;
  }
} // namespace

int main(int argc, char** argv)
{
  const std::string tmp = argc > 1 ? argv[1] : ".";
  const auto lib = makeLibrary();
  const auto run = makeRun(lib);

  const Arm serial = runArm(lib, run, tmp, false, false, 1);
  std::printf("serial: %zu result bytes, %zu gate-log bytes, Gate C rejected %zu, "
              "%zu anchors dropped by the cap\n",
              serial.bytes.size(), serial.gate_log.size(), serial.gate_rejected,
              serial.anchors_dropped);
  check(serial.gate_rejected > 0, "the serial reference has Gate C rejections (the gate armed)");
  check(serial.anchors_dropped > 0, "and the mass-anchor cap was reached");
  check(serial.gate_log.find("\t1\t1\n") != std::string::npos &&
        serial.gate_log.find("\t0\t1\n") != std::string::npos,
        "and the gate log holds decisions from both sides of the arming boundary");

  struct Spec { bool parallel; unsigned threads; const char* name; };
  for (const Spec s : {Spec{false, 1, "addBatch, flag off"},
                       Spec{true, 8, "addBatch, -parallel_sink, 8 threads"},
                       Spec{true, 48, "addBatch, -parallel_sink, 48 threads"}})
  {
    const Arm a = runArm(lib, run, tmp, true, s.parallel, s.threads);
    std::printf("%s: %zu speculative bodies discarded at commit\n", s.name, a.discarded);
    check(a.bytes == serial.bytes, std::string(s.name) + ": Result byte-identical to serial add()");
    check(a.reasons == serial.reasons, std::string(s.name) + ": terminal reasons identical");
    check(a.gate_log == serial.gate_log, std::string(s.name) + ": Gate C log byte-identical");
    if (s.parallel)
    {
      check(a.discarded > 0, std::string(s.name) +
            ": the null armed INSIDE a batch, so commit-time rejection of speculated work was exercised");
    }
  }

  if (failures == 0) { std::printf("the parallel sink is byte-identical to the serial sink\n"); }
  return failures == 0 ? 0 : 1;
}
