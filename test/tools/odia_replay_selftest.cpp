// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// v1.18 learner fixes and the offline replay: what can be checked without a run.
//
//   odia_replay_selftest <table.tsv>
//
// 1. Scores a synthetic candidate table with the ENGINE's parameter construction
//    (classifierParamsFor at the engine CLI's defaults) and writes it in -out_scorer_input format
//    (17 digits, Group.Index) to <table.tsv>. The companion ctest then runs odia_scorer_replay
//    -expect_identical on it: the replay must reproduce every DScore bit-for-bit, which exercises
//    the 17-digit round trip, the group-index fold assignment and the shared params builder.
// 2. Runs every v1.18 flag through scoreSemiSupervisedLDA on the same table: each must train,
//    give finite scores and identify something (counted at q <= 0.05: with ~1,200 decoy groups the
//    1% operating point rests on a handful of decoys and swings to 0 on a near-tie, which is a
//    property of the threshold, not of the flag), and each must CHANGE the output (else it is not
//    wired). Tail calibration must additionally leave every group's within-group order unchanged
//    (it is a strictly monotone per-fold transform of the same models' scores).
// 3. GBT-level unit checks that pin the defects the flags fix: max_depth 1 natively grows NO split
//    (the off-by-one) and grows one with depth_fix; stop_on_stump turns that silent no-split fit
//    into a reported failure; missing_bin lets a tree route NaN on its own.
#include <odia/PeakGroupScorer.h>
#include <odia/scoring/lda.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <locale>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
int failures = 0;
void check(bool ok, const char* what)
{
  std::printf("%s: %s\n", ok ? "ok  " : "FAIL", what);
  if (!ok) { ++failures; }
}

std::size_t idsAt5pct(const ODIA::Scoring::ScoredGroups& s, const std::vector<int>& labels,
                      const std::vector<long long>& group)
{
  std::unordered_map<long long, double> best;
  for (std::size_t i = 0; i < labels.size(); ++i)
  {
    if (labels[i] != 1) { continue; }
    auto it = best.emplace(group[i], s.qvalue[i]).first;
    it->second = std::min(it->second, s.qvalue[i]);
  }
  std::size_t n = 0;
  for (const auto& kv : best) { n += kv.second <= 0.05; }
  return n;
}
} // namespace

int main(int argc, char** argv)
{
  if (argc != 2) { std::fprintf(stderr, "usage: odia_replay_selftest <table.tsv>\n"); return 2; }
  using ODIA::PeakGroupScorer;
  namespace S = ODIA::Scoring;
  const auto& names = PeakGroupScorer::subScoreNames();
  const std::size_t m = names.size();

  // ---- synthetic candidate table -----------------------------------------------------------------
  // Targets carry 1-6 candidates, decoys 2-11, so decoy candidate-count matching has work to do.
  // 40% of targets have a true candidate; it carries the signal on CORR_SUM, XCORR_SHAPE and
  // MS1_COELUTION, and MS1_COELUTION is missing (NaN) far more often on wrong rows.
  std::mt19937 rng(20261002);
  std::normal_distribution<double> nz(0.0, 1.0);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  std::vector<std::vector<double>> X;
  std::vector<int> labels;
  std::vector<long long> group;
  std::vector<std::string> pid;
  const int n_prec = 2400;
  for (int p = 0; p < n_prec; ++p)
  {
    const bool target = (p % 2) == 0;
    // library-index-like ids: NOT the first-occurrence order, so the fold assignment depends on
    // the Group.Index column actually being read
    const long long gid = (static_cast<long long>(p) * 7919LL) % 100003LL;
    const int k = target ? 1 + static_cast<int>(rng() % 6) : 2 + static_cast<int>(rng() % 10);
    const int truth = (target && u(rng) < 0.4) ? static_cast<int>(rng() % static_cast<unsigned>(k)) : -1;
    for (int c = 0; c < k; ++c)
    {
      const bool real = (c == truth);
      std::vector<double> row(m);
      for (std::size_t j = 0; j < m; ++j) { row[j] = nz(rng); }
      row[PeakGroupScorer::CORR_SUM] += real ? 3.5 : 0.0;
      row[PeakGroupScorer::XCORR_SHAPE] += real ? 2.0 : 0.0;
      row[PeakGroupScorer::MS1_COELUTION] += real ? 1.0 : 0.0;
      if (u(rng) < (real ? 0.05 : 0.35)) { row[PeakGroupScorer::MS1_COELUTION] = std::nan(""); }
      row[PeakGroupScorer::CAND_COUNT] = static_cast<double>(k);
      X.push_back(std::move(row));
      labels.push_back(target ? 1 : 0);
      group.push_back(gid);
      pid.push_back("PEPTIDE" + std::to_string(p) + "2");
    }
  }
  const std::size_t n = X.size();

  // ---- the engine's parameters at its CLI defaults (as odia_scorer_replay sets them) -------------
  PeakGroupScorer::Options o;
  o.classifier = "gbt";
  o.classifier_seed = 42;
  o.classifier_iterations = 3;
  o.train_fdr_initial = 0.15;
  o.train_fdr = 0.05;
  o.match_decoy_candidate_counts = true;
  o.coelution_picking = true;
  o.classifier_stop_jaccard = 0.98;
  o.classifier_stop_shrink_floor = 0.01;
  o.classifier_stop_patience = 2;
#ifdef _OPENMP
  o.threads = static_cast<unsigned>(std::max(1, omp_get_max_threads()));
#else
  o.threads = 1;
#endif
  const S::LDAParams base = ODIA::classifierParamsFor(o);
  check(!base.match_training_draw && !base.gbt_keep_missing && !base.fold_tail_calibration &&
        base.oof_repeats == 1 && !base.gbt.class_balance && !base.gbt.missing_bin &&
        !base.gbt.depth_fix && !base.gbt.stop_on_stump,
        "every v1.18 field is at its default when no v1.18 option is set");
  const auto ref = S::scoreSemiSupervisedLDA(X, labels, group, base, "", "", &names);
  const std::size_t ref_ids = idsAt5pct(ref, labels, group);
  std::printf("baseline: trained %d skipped %d, %zu ids at q<=0.05\n", ref.n_iterations_trained,
              ref.n_iterations_skipped, ref_ids);
  check(ref.n_iterations_trained > 0 && ref_ids > 200, "baseline trains and identifies");

  // Explicitly-set defaults are the default: K=1 is the native path.
  {
    S::LDAParams p = base;
    p.oof_repeats = 1;
    const auto r = S::scoreSemiSupervisedLDA(X, labels, group, p, "", "", &names);
    check(r.dscore == ref.dscore && r.qvalue == ref.qvalue, "oof_repeats = 1 reproduces the native scores bit-for-bit");
  }

  // ---- write the -out_scorer_input-format table -------------------------------------------------
  {
    std::ofstream out(argv[1]);
    out.imbue(std::locale::classic());
    out.precision(17);
    out << "Precursor.Id\tDecoy\tGroup.Index\tDScore\tQValue\tPEP";
    for (const auto& nm : names) { out << '\t' << nm; }
    out << '\n';
    for (std::size_t i = 0; i < n; ++i)
    {
      out << pid[i] << '\t' << (labels[i] == 1 ? 0 : 1) << '\t' << group[i] << '\t' << ref.dscore[i]
          << '\t' << ref.qvalue[i] << '\t' << ref.pep[i];
      for (const double v : X[i]) { out << '\t' << v; }
      out << '\n';
    }
    check(static_cast<bool>(out), "wrote the replay table");
  }

  // ---- every v1.18 flag: trains, finite, identifies, and changes the output ----------------------
  struct Arm { const char* name; void (*set)(S::LDAParams&); };
  const Arm arms[] = {
    {"class_balance", [](S::LDAParams& p) { p.gbt.class_balance = true; }},
    {"matched_train_draw", [](S::LDAParams& p) { p.match_training_draw = true; }},
    {"keep_missing", [](S::LDAParams& p) { p.gbt_keep_missing = true; }},
    {"depth_fix", [](S::LDAParams& p) { p.gbt.depth_fix = true; }},
    {"stop_on_stump", [](S::LDAParams& p) { p.gbt.stop_on_stump = true; }},
    {"tail_calibration", [](S::LDAParams& p) { p.fold_tail_calibration = true; }},
    {"oof_repeats_3", [](S::LDAParams& p) { p.oof_repeats = 3; }},
  };
  for (const auto& a : arms)
  {
    S::LDAParams p = base;
    a.set(p);
    const auto r = S::scoreSemiSupervisedLDA(X, labels, group, p, "", "", &names);
    bool finite = true;
    for (const double d : r.dscore) { finite = finite && std::isfinite(d); }
    const std::size_t ids = idsAt5pct(r, labels, group);
    std::printf("%-20s trained %d skipped %d, %zu ids at q<=0.05\n", a.name, r.n_iterations_trained,
                r.n_iterations_skipped, ids);
    const std::string what = std::string(a.name) + ": trains, finite, identifies";
    check(r.n_iterations_trained > 0 && finite && ids > 200, what.c_str());
    // stop_on_stump only acts on a stump, and nothing here grows one: it must be a no-op instead.
    if (std::string(a.name) == "stop_on_stump")
    { check(r.dscore == ref.dscore, "stop_on_stump: no-op when every tree splits"); }
    else
    {
      const std::string w2 = std::string(a.name) + ": changes the scores (wired)";
      check(r.dscore != ref.dscore, w2.c_str());
    }
    if (std::string(a.name) == "tail_calibration")
    {
      // same models, strictly monotone per-fold transform: every group's internal order survives
      std::unordered_map<long long, std::vector<std::size_t>> rows_of;
      for (std::size_t i = 0; i < n; ++i) { rows_of[group[i]].push_back(i); }
      bool order_kept = true;
      for (const auto& kv : rows_of)
      {
        const auto& rs = kv.second;
        for (std::size_t x = 0; x < rs.size(); ++x)
        {
          for (std::size_t y = x + 1; y < rs.size(); ++y)
          {
            const double d0 = ref.dscore[rs[x]] - ref.dscore[rs[y]];
            const double d1 = r.dscore[rs[x]] - r.dscore[rs[y]];
            if ((d0 > 0) != (d1 > 0) || (d0 < 0) != (d1 < 0)) { order_kept = false; }
          }
        }
      }
      check(order_kept, "tail_calibration: within-group order identical to the native pooling");
    }
  }

  // ---- model I/O with the loop-level flags ---------------------------------------------------------
  // A model SAVE is output-only: with a loop-level flag the file is skipped and the scores are those
  // of the same run without the save (it used to zero every row). APPLYING a model is refused.
  {
    const std::string path = std::string(argv[1]) + ".model";
    std::remove(path.c_str());
    for (int k = 0; k < 3; ++k)
    {
      S::LDAParams p = base;
      if (k == 0) { p.fold_tail_calibration = true; }
      if (k == 1) { p.oof_repeats = 2; }
      if (k == 2) { p.gbt_keep_missing = true; }
      const auto plain = S::scoreSemiSupervisedLDA(X, labels, group, p, "", "", &names);
      const auto saved = S::scoreSemiSupervisedLDA(X, labels, group, p, path, "", &names);
      const std::string w = std::string(k == 0 ? "tail_calibration" : k == 1 ? "oof_repeats_2"
                                                                       : "keep_missing") +
                            " + model_out: scores unchanged, model file not written";
      check(saved.n_iterations_trained > 0 && saved.dscore == plain.dscore &&
              saved.qvalue == plain.qvalue && !std::ifstream(path).good(),
            w.c_str());
    }
    S::LDAParams p = base;
    p.fold_tail_calibration = true;
    const auto r = S::scoreSemiSupervisedLDA(X, labels, group, p, "", "/nonexistent.model", &names);
    check(r.n_iterations_trained == 0 && r.n_iterations_skipped == 1,
          "tail calibration with a model APPLY is refused, not silently applied");
  }

  // ---- GBT unit checks ----------------------------------------------------------------------------
  {
    // planted single-feature signal
    std::vector<std::vector<double>> G;
    std::vector<std::size_t> pos, neg;
    for (int i = 0; i < 4000; ++i)
    {
      const bool y = (i % 4) == 0;
      G.push_back({(y ? 1.5 : 0.0) + nz(rng), nz(rng)});
      (y ? pos : neg).push_back(static_cast<std::size_t>(i));
    }
    auto spread = [&](const S::GBT& g) {
      double lo = 1e300, hi = -1e300;
      for (const auto& x : G) { const double s = g.score(x); lo = std::min(lo, s); hi = std::max(hi, s); }
      return hi - lo;
    };
    S::GBTParams gp;
    gp.max_depth = 1;
    gp.n_trees = 20;
    S::GBT native;
    const bool native_ok = native.fit(G, pos, neg, gp);
    check(native_ok && spread(native) == 0.0,
          "max_depth 1, native: fit 'succeeds' yet no tree splits (the off-by-one; every score equal)");
    S::GBTParams gs = gp;
    gs.stop_on_stump = true;
    S::GBT stump;
    check(!stump.fit(G, pos, neg, gs), "stop_on_stump: that no-split fit is reported as a failure");
    S::GBTParams gd = gp;
    gd.depth_fix = true;
    S::GBT fixed;
    check(fixed.fit(G, pos, neg, gd) && spread(fixed) > 0.0, "depth_fix: max_depth 1 grows one split level");
    S::GBTParams gb;
    gb.class_balance = true;
    S::GBT bal;
    check(bal.fit(G, pos, neg, gb) && spread(bal) > 0.0, "class_balance: fits the planted signal");

    // missingness IS the signal: positives NaN, negatives finite and spread over the whole range
    std::vector<std::vector<double>> M;
    std::vector<std::size_t> mp, mn;
    for (int i = 0; i < 4000; ++i)
    {
      const bool y = (i % 5) == 0;
      M.push_back({y ? std::nan("") : nz(rng)});
      (y ? mp : mn).push_back(static_cast<std::size_t>(i));
    }
    auto separated = [&](const S::GBT& g) {
      double min_nan = 1e300, max_fin = -1e300;
      for (const auto& x : M)
      {
        const double s = g.score(x);
        if (std::isnan(x[0])) { min_nan = std::min(min_nan, s); } else { max_fin = std::max(max_fin, s); }
      }
      return min_nan > max_fin;
    };
    S::GBTParams g0;
    S::GBT nat;
    nat.fit(M, mp, mn, g0);
    S::GBTParams g1;
    g1.missing_bin = true;
    S::GBT mb;
    mb.fit(M, mp, mn, g1);
    check(!separated(nat), "native: NaN shares the top value bin, so it cannot be isolated");
    check(separated(mb), "missing_bin: every NaN row outscores every finite row");
    std::ofstream sink("/dev/null");
    check(!mb.save(sink), "missing_bin: the version-1 model format refuses to save it");
  }

  std::printf("%d failure(s)\n", failures);
  return failures == 0 ? 0 : 1;
}
