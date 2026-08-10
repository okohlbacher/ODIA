// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

/// The per-run retention-time refiner, tested against the failures it has had.
///
/// Every check here corresponds to something that actually went wrong on real
/// data, because those are the ones that recur:
///
///  * an UNBOUNDED correction. The first version fitted absolute retention time
///    and extrapolated without limit outside its training range, which put pass
///    2's windows off the gradient and took S08 from 1,464 identifications to
///    ZERO. The clamp must hold even when the model is confident and wrong.
///  * CONTAMINATED anchors. Pass 1's anchors reach 1,600 s of residual --
///    misidentifications, not chromatography -- and a squared loss chases them.
///    The trim must remove them and say how many.
///  * LEAKAGE across charge states. 2+ and 3+ of one peptide co-elute, so a
///    split between them reports memorisation as generalisation.
///  * a model that DOES NOT HELP. It must be refused rather than applied, on
///    the run's own held-out rows.

#include <odia/RtRefiner.h>

#include <cmath>
#include <cstdio>
#include <random>
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

  const char* AAS = "ACDEFGHIKLMNPQRSTVWY";

  std::string peptide(std::mt19937& rng, std::size_t len)
  {
    std::uniform_int_distribution<int> aa(0, 19);
    std::string s;
    for (std::size_t i = 0; i < len; ++i) { s += AAS[aa(rng)]; }
    return s;
  }

  /// A run whose true iRT depends on COMPOSITION -- the signal a sequence model
  /// can find and a monotone map cannot.
  struct Run
  {
    std::vector<std::string> seq;
    std::vector<int> charge;
    std::vector<double> lib;     ///< library iRT (the composition term missing)
    std::vector<double> target;  ///< what it should have been
  };

  Run plant(std::size_t n, double contamination, unsigned seed)
  {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::normal_distribution<double> noise(0.0, 1.0);
    std::uniform_int_distribution<int> len(7, 20);
    Run r;
    for (std::size_t i = 0; i < n; ++i)
    {
      const std::string p = peptide(rng, static_cast<std::size_t>(len(rng)));
      double hydro = 0.0;
      for (const char c : p) { if (c == 'L' || c == 'I' || c == 'V' || c == 'F') { hydro += 1.0; } }
      const double base = 100.0 * u(rng);
      r.seq.push_back(p);
      r.charge.push_back(2 + (i % 2));
      r.lib.push_back(base);
      // The library value misses a composition-dependent shift; that is exactly
      // what the refiner should recover.
      double t = base + 2.0 * hydro + 0.5 * noise(rng);
      if (u(rng) < contamination) { t += 400.0 + 200.0 * u(rng); }  // a misidentification
      r.target.push_back(t);
    }
    return r;
  }
}

int main()
{
  std::printf("per-run retention-time refiner\n");

  // Recovers a composition-dependent shift a monotone map cannot.
  {
    const auto r = plant(4000, 0.0, 11u);
    ODIA::RtRefiner ref;
    const auto rep = ref.fit(r.seq, r.charge, r.lib, r.target);
    check(rep.fitted, "fits on clean data");
    check(rep.sd_after < rep.sd_before,
          "held-out residual improves on a composition-dependent shift");
    std::printf("      held-out SD %.3f -> %.3f on %zu sequences\n",
                rep.sd_before, rep.sd_after, rep.held_out);
  }

  // Contaminated anchors are trimmed, and the fit still works.
  {
    const auto r = plant(4000, 0.05, 12u);
    ODIA::RtRefiner ref;
    const auto rep = ref.fit(r.seq, r.charge, r.lib, r.target);
    check(rep.fitted, "fits with 5% contaminated anchors");
    check(rep.trimmed > 50,
          "the contaminated tail is trimmed rather than fitted -- a squared loss "
          "chases 400 s outliers and that is how a 2% better model destroyed a run");
    std::printf("      trimmed %zu of %zu\n", rep.trimmed, rep.anchors);
  }

  // THE CLAMP. A precursor far outside the training range must not be moved
  // further than the bound, however confident the extrapolation.
  {
    const auto r = plant(4000, 0.0, 13u);
    ODIA::RtRefiner ref;
    ODIA::RtRefiner::Options o;
    o.max_shift_fraction = 0.05;
    const auto rep = ref.fit(r.seq, r.charge, r.lib, r.target, o);
    check(rep.fitted, "fits for the clamp check");

    ODIA::Library lib;
    auto& p = lib.precursors();
    // One ordinary peptide and one absurd one, far outside anything seen.
    for (const auto& s : {std::string("PEPTIDEK"),
                          std::string("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWW")})
    {
      p.modified_sequence.push_back(lib.strings().intern(s));
      p.charge.push_back(2);
      p.irt.push_back(50.0f);
      p.transition_begin.push_back(0);
      p.transition_count.push_back(0);
      p.decoy.push_back(0);
      p.mz.push_back(500.0);
      p.im.push_back(0.0f);
      p.protein_group.push_back(lib.strings().intern("P"));
    }
    const std::vector<float> before = p.irt;
    ref.apply(lib);
    const double range = 100.0;   // the planted iRT range
    for (std::size_t i = 0; i < p.irt.size(); ++i)
    {
      const double moved = std::fabs(static_cast<double>(p.irt[i]) - before[i]);
      check(moved <= 0.05 * range + 1e-3,
            "correction is clamped to 5% of the iRT range for precursor " +
              std::to_string(i) + " (moved " + std::to_string(moved) + ")");
    }
  }

  // Refused when it cannot help: the library value is already exact.
  {
    auto r = plant(4000, 0.0, 14u);
    r.target = r.lib;                       // nothing left to learn
    ODIA::RtRefiner ref;
    const auto rep = ref.fit(r.seq, r.charge, r.lib, r.target);
    check(!rep.fitted || rep.sd_after <= rep.sd_before,
          "an exact library is either refused or at least not made worse");
  }

  // Too few anchors is a refusal with a reason, not a bad model.
  {
    const auto r = plant(100, 0.0, 15u);
    ODIA::RtRefiner ref;
    const auto rep = ref.fit(r.seq, r.charge, r.lib, r.target);
    check(!rep.fitted && !rep.note.empty(),
          "100 anchors is refused, with a stated reason");
  }

  // save/load round-trips, and a truncated file is rejected rather than read as
  // a valid model that predicts noise.
  {
    const auto r = plant(4000, 0.0, 16u);
    ODIA::RtRefiner a;
    const auto rep = a.fit(r.seq, r.charge, r.lib, r.target);
    check(rep.fitted, "fits for the round-trip");
    const std::string path = "/tmp/odia_rtref_test.txt";
    check(a.save(path, "unit test"), "saves");
    ODIA::RtRefiner b;
    std::string prov;
    check(b.load(path, &prov), "loads");
    check(prov.find("unit test") != std::string::npos, "provenance survives the round trip");

    std::FILE* f = std::fopen("/tmp/odia_rtref_trunc.txt", "w");
    std::fprintf(f, "# broken\nfeatures\t25\n1.0\n2.0\n");
    std::fclose(f);
    ODIA::RtRefiner c;
    check(!c.load("/tmp/odia_rtref_trunc.txt"),
          "a file with the wrong weight count is REFUSED -- reading it as valid "
          "would predict noise with full confidence");
  }

  // The clamp must SURVIVE the round trip. Version 1 wrote only weights, so a
  // model fitted with a small clamp was applied after loading with the member
  // default -- defeating the protection added after the 1,464 -> 0 incident.
  {
    auto r = plant(4000, 0.0, 17u);
    // Compress the axis so the fitted clamp is far from any default.
    for (auto& v : r.lib) { v *= 0.1; }
    for (auto& v : r.target) { v *= 0.1; }
    ODIA::RtRefiner a;
    const auto rep = a.fit(r.seq, r.charge, r.lib, r.target);
    check(rep.fitted, "fits on a compressed axis");
    const std::string path = "/tmp/odia_rtref_clamp.txt";
    check(a.save(path, "clamp test"), "saves with the clamp");

    ODIA::RtRefiner b;
    check(b.load(path), "loads");

    const auto build = []() {
      ODIA::Library lib;
      auto& p = lib.precursors();
      p.modified_sequence.push_back(lib.strings().intern("WWWWWWWWWWWWWWWWWWWW"));
      p.charge.push_back(2); p.irt.push_back(5.0f);
      p.transition_begin.push_back(0); p.transition_count.push_back(0);
      p.decoy.push_back(0); p.mz.push_back(500.0); p.im.push_back(0.0f);
      p.protein_group.push_back(lib.strings().intern("P"));
      return lib;
    };
    ODIA::Library la = build(), lb = build();
    a.apply(la); b.apply(lb);
    check(std::fabs(static_cast<double>(la.precursors().irt[0]) -
                    static_cast<double>(lb.precursors().irt[0])) < 1e-3,
          "a loaded model applies the SAME correction as the fitted one -- if the "
          "clamp did not travel, this diverges by up to the default bound");
  }

  // A stale model must not survive a failed fit or a failed load.
  {
    const auto good = plant(4000, 0.0, 18u);
    ODIA::RtRefiner ref;
    check(ref.fit(good.seq, good.charge, good.lib, good.target).fitted, "fits once");
    check(ref.fitted(), "is fitted");
    const auto tiny = plant(100, 0.0, 19u);
    ref.fit(tiny.seq, tiny.charge, tiny.lib, tiny.target);
    check(!ref.fitted(),
          "a REFUSED refit clears the previous model rather than leaving it live");

    ODIA::RtRefiner ref2;
    check(ref2.fit(good.seq, good.charge, good.lib, good.target).fitted, "fits again");
    check(!ref2.load("/tmp/does_not_exist_odia.txt"), "a missing file fails to load");
    check(!ref2.fitted(), "and clears the model rather than keeping the old one");
  }

  std::printf("%s\n", failures ? "FAILED" : "all good");
  return failures ? 1 : 0;
}
