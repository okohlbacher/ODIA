// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// odia_scorer_replay -- re-run ODIA's FINAL-pass classifier offline on a scored table.
//
// WHY. Every question about the learner (objective, negatives, missingness, depth, fold pooling,
// ensembles) used to cost a full engine run, and an engine run cannot separate a learner change
// from everything upstream of it. The classifier's input is a plain table -- one row per candidate
// peak group: its precursor, target/decoy, and the var_* sub-scores -- so the learner can be re-run
// on that table alone, with any options, in minutes, over as many fold seeds as wanted.
//
// WHAT IT RUNS. Exactly the engine's code: ODIA::classifierParamsFor() (the function the engine's
// PeakGroupScorer::fitAndAssign_ calls) turns the options into LDAParams, and
// Scoring::scoreSemiSupervisedLDA (include/odia/scoring/lda.h + gbt.h) trains and scores. Nothing
// about the learner is reimplemented here.
//
// INPUT. A TSV with columns Precursor.Id, Decoy and every var_* sub-score of this binary's
// subScoreNames(), in any column order (they are read by name and handed to the classifier in
// subScoreNames() order, which the seed mask and the sign constraints index into). Two sources:
//   * -out_scorer_input of an engine run (v1.18): 17 significant digits and a Group.Index column,
//     so the replay sees bit-for-bit the engine's classifier input and its fold assignment.
//   * a plain -out score table: 6 significant digits and no library index. The features are then
//     rounded, and the fold assignment -- which sorts groups by LIBRARY PRECURSOR INDEX before the
//     seeded shuffle -- cannot be reproduced unless -foldmap supplies the indices (the FOLDMAP block
//     of a -classifier_model_out file from the same run lists them in exactly this table's group
//     order). Without either, groups are numbered by first occurrence and the replay is a
//     different-but-equivalent partition; it says so.
//
// OUTPUT. -out: Row, Precursor.Id, Decoy, DScore, QValue, PEP per input row (17 digits). Stdout: a
// summary, and -- when the input carries the engine's DScore/QValue -- the agreement: bit-identical
// rows, rows equal at the plain TSV's 6-digit print precision, max |delta|, the q<=0.01 sets.
#include <odia/PeakGroupScorer.h>
#include <odia/scoring/lda.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <locale>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{

void usage()
{
  std::fprintf(stderr,
    "usage: odia_scorer_replay -in <scored.tsv> [-out <replay.tsv>] [options]\n"
    "  options take the ENGINE's names and the engine CLI's defaults:\n"
    "  -classifier gbt|xgboost|lda|nn (gbt)   -classifier_seed <n> (42)\n"
    "  -classifier_iterations <n> (3)  -train_fdr_initial <q> (0.15)  -train_fdr <q> (0.05)  -use_pi0\n"
    "  -no_match_decoy_n  -fold_pool_rank  -picker coelution|amplitude|openswath|union|union_openswath\n"
    "  -gbt_max_depth -gbt_n_trees -gbt_learning_rate -gbt_min_child_rows -gbt_lambda\n"
    "  -gbt_max_delta_step -gbt_warmup_rounds <v> (0 = model default)\n"
    "  -gbt_intercept_zero -gbt_clip_gain -gbt_fixed_bins\n"
    "  -classifier_stop_on_composition -classifier_stop_jaccard (0.98) -classifier_stop_shrink_floor (0.01)\n"
    "  -classifier_stop_patience (2) -classifier_iteration_log\n"
    "  v1.18: -classifier_class_balance -classifier_matched_train_draw -gbt_missing_bin -gbt_depth_fix\n"
    "         -gbt_stop_on_stump -fold_tail_calibration -classifier_oof_repeats <k> (1)\n"
    "  -foldmap <model file>   group ids from the FOLDMAP of a -classifier_model_out file\n"
    "  -threads <n>            OpenMP threads (results do not depend on it)\n"
    "  -expect_identical       exit 4 unless every DScore is bit-identical to the input's and the\n"
    "                          q<=0.01 sets agree (the replay-reproduces-the-engine gate)\n");
}

std::vector<std::string> splitTabs(const std::string& line)
{
  std::vector<std::string> out;
  std::size_t b = 0;
  for (;;)
  {
    const std::size_t e = line.find('\t', b);
    if (e == std::string::npos) { out.push_back(line.substr(b)); break; }
    out.push_back(line.substr(b, e - b));
    b = e + 1;
  }
  return out;
}

double parseDouble(const std::string& s)
{
  // strtod reads what an ostream prints, including "nan", "-nan" and "inf".
  const char* c = s.c_str();
  char* end = nullptr;
  const double v = std::strtod(c, &end);
  if (end == c) { return std::numeric_limits<double>::quiet_NaN(); }
  return v;
}

/// The default ostream rendering (precision 6) the plain -out table uses.
std::string printed6(double v)
{
  std::ostringstream o;
  o.imbue(std::locale::classic());
  o << v;
  return o.str();
}

} // namespace

int main(int argc, char** argv)
{
  std::map<std::string, std::string> opt;
  std::unordered_set<std::string> flags;
  const std::unordered_set<std::string> known_flags = {
    "-use_pi0", "-no_match_decoy_n", "-fold_pool_rank", "-gbt_intercept_zero", "-gbt_clip_gain",
    "-gbt_fixed_bins", "-classifier_stop_on_composition", "-classifier_iteration_log",
    "-classifier_class_balance", "-classifier_matched_train_draw", "-gbt_missing_bin",
    "-gbt_depth_fix", "-gbt_stop_on_stump", "-fold_tail_calibration", "-expect_identical"};
  const std::unordered_set<std::string> known_opts = {
    "-in", "-out", "-classifier", "-classifier_seed", "-classifier_iterations", "-train_fdr_initial",
    "-train_fdr", "-picker", "-gbt_max_depth", "-gbt_n_trees", "-gbt_learning_rate",
    "-gbt_min_child_rows", "-gbt_lambda", "-gbt_max_delta_step", "-gbt_warmup_rounds",
    "-classifier_stop_jaccard", "-classifier_stop_shrink_floor", "-classifier_stop_patience",
    "-classifier_oof_repeats", "-foldmap", "-threads"};
  for (int i = 1; i < argc; ++i)
  {
    const std::string a = argv[i];
    if (a == "-h" || a == "--help") { usage(); return 0; }
    if (known_flags.count(a)) { flags.insert(a); continue; }
    if (known_opts.count(a) && i + 1 < argc) { opt[a] = argv[++i]; continue; }
    std::fprintf(stderr, "unknown or incomplete argument: %s\n", a.c_str());
    usage();
    return 2;
  }
  if (!opt.count("-in")) { usage(); return 2; }
  auto get = [&](const char* k, const std::string& d) { return opt.count(k) ? opt[k] : d; };
  auto has = [&](const char* k) { return flags.count(k) > 0; };

#ifdef _OPENMP
  if (opt.count("-threads")) { omp_set_num_threads(std::max(1, std::atoi(opt["-threads"].c_str()))); }
#endif

  // ---- options -> PeakGroupScorer::Options, as the engine's scoringOptions_ fills them ----------
  ODIA::PeakGroupScorer::Options o;
  o.classifier = get("-classifier", "gbt");
  if (o.classifier == "percolator")
  { std::fprintf(stderr, "the percolator engine is not replayable here\n"); return 2; }
  if (o.classifier != "gbt" && o.classifier != "xgboost" && o.classifier != "lda" && o.classifier != "nn")
  { std::fprintf(stderr, "unknown -classifier %s\n", o.classifier.c_str()); return 2; }
  o.classifier_seed = std::atoi(get("-classifier_seed", "42").c_str());
  o.classifier_iterations = std::atoi(get("-classifier_iterations", "3").c_str());
  o.train_fdr_initial = parseDouble(get("-train_fdr_initial", "0.15"));
  o.train_fdr = parseDouble(get("-train_fdr", "0.05"));
  o.use_pi0 = has("-use_pi0");
  o.match_decoy_candidate_counts = !has("-no_match_decoy_n");
  o.fold_pool_rank = has("-fold_pool_rank");
  {
    const std::string picker = get("-picker", "coelution");
    o.union_picking = picker == "union" || picker == "union_openswath";
    o.openswath_picking = picker == "openswath" || picker == "union_openswath";
    o.coelution_picking = picker != "amplitude";
  }
  o.gbt_max_depth = std::atoi(get("-gbt_max_depth", "0").c_str());
  o.gbt_n_trees = std::atoi(get("-gbt_n_trees", "0").c_str());
  o.gbt_learning_rate = parseDouble(get("-gbt_learning_rate", "0"));
  o.gbt_min_child_rows = std::atoi(get("-gbt_min_child_rows", "0").c_str());
  o.gbt_lambda = parseDouble(get("-gbt_lambda", "0"));
  o.gbt_max_delta_step = parseDouble(get("-gbt_max_delta_step", "0"));
  o.gbt_warmup_rounds = std::atoi(get("-gbt_warmup_rounds", "0").c_str());
  o.gbt_intercept_zero = has("-gbt_intercept_zero");
  o.gbt_clip_gain = has("-gbt_clip_gain");
  o.gbt_fixed_bins = has("-gbt_fixed_bins");
  o.classifier_stop_on_composition = has("-classifier_stop_on_composition");
  o.classifier_stop_jaccard = parseDouble(get("-classifier_stop_jaccard", "0.98"));
  o.classifier_stop_shrink_floor = parseDouble(get("-classifier_stop_shrink_floor", "0.01"));
  o.classifier_stop_patience = std::atoi(get("-classifier_stop_patience", "2").c_str());
  o.classifier_iteration_log = has("-classifier_iteration_log");
  o.classifier_class_balance = has("-classifier_class_balance");
  o.classifier_matched_train_draw = has("-classifier_matched_train_draw");
  o.gbt_missing_bin = has("-gbt_missing_bin");
  o.gbt_depth_fix = has("-gbt_depth_fix");
  o.gbt_stop_on_stump = has("-gbt_stop_on_stump");
  o.fold_tail_calibration = has("-fold_tail_calibration");
  o.classifier_oof_repeats = std::max(1, std::atoi(get("-classifier_oof_repeats", "1").c_str()));
#ifdef _OPENMP
  o.threads = static_cast<unsigned>(std::max(1, omp_get_max_threads()));
#else
  o.threads = 1;
#endif
  const ODIA::Scoring::LDAParams params = ODIA::classifierParamsFor(o);

  // ---- read the table ---------------------------------------------------------------------------
  std::ifstream in(opt["-in"]);
  if (!in) { std::fprintf(stderr, "cannot open %s\n", opt["-in"].c_str()); return 1; }
  std::string line;
  if (!std::getline(in, line)) { std::fprintf(stderr, "empty input\n"); return 1; }
  const auto header = splitTabs(line);
  std::unordered_map<std::string, std::size_t> col;
  for (std::size_t i = 0; i < header.size(); ++i) { col.emplace(header[i], i); }
  auto need = [&](const std::string& c) -> long {
    const auto it = col.find(c);
    return it == col.end() ? -1L : static_cast<long>(it->second);
  };
  const long c_id = need("Precursor.Id"), c_dec = need("Decoy");
  const long c_grp = need("Group.Index"), c_ds = need("DScore"), c_q = need("QValue");
  if (c_id < 0 || c_dec < 0) { std::fprintf(stderr, "input lacks Precursor.Id/Decoy\n"); return 1; }
  const auto& names = ODIA::PeakGroupScorer::subScoreNames();
  std::vector<std::size_t> feat_col;
  for (const auto& n : names)
  {
    const long c = need(n);
    if (c < 0)
    {
      std::fprintf(stderr, "input lacks sub-score column %s: it was written by a binary with a "
                           "different sub-score set, and the seed mask / sign constraints index by "
                           "position -- refusing\n", n.c_str());
      return 1;
    }
    feat_col.push_back(static_cast<std::size_t>(c));
  }
  for (const auto& h : header)
  {
    if (h.rfind("var_", 0) == 0 && std::find(names.begin(), names.end(), h) == names.end())
    { std::fprintf(stderr, "input has a sub-score %s this binary does not know -- refusing\n", h.c_str()); return 1; }
  }
  const std::size_t m = names.size();

  std::vector<std::vector<double>> features;
  std::vector<int> labels;
  std::vector<long long> group;
  std::vector<std::string> ids;
  std::vector<double> eng_ds, eng_q;
  std::vector<std::string> eng_ds_text;
  std::unordered_map<std::string, long long> first_index;   // (id, decoy) -> first-occurrence index
  std::vector<long long> first_occ;                         // per row
  std::size_t lineno = 1;
  while (std::getline(in, line))
  {
    ++lineno;
    if (line.empty()) { continue; }
    const auto f = splitTabs(line);
    if (f.size() != header.size())
    { std::fprintf(stderr, "line %zu: %zu fields, header has %zu\n", lineno, f.size(), header.size()); return 1; }
    const int decoy = std::atoi(f[static_cast<std::size_t>(c_dec)].c_str());
    labels.push_back(decoy ? 0 : 1);
    ids.push_back(f[static_cast<std::size_t>(c_id)]);
    std::vector<double> row(m);
    for (std::size_t j = 0; j < m; ++j) { row[j] = parseDouble(f[feat_col[j]]); }
    features.push_back(std::move(row));
    const std::string key = ids.back() + (decoy ? "\t1" : "\t0");
    const auto ins = first_index.emplace(key, static_cast<long long>(first_index.size()));
    first_occ.push_back(ins.first->second);
    if (c_grp >= 0) { group.push_back(std::atoll(f[static_cast<std::size_t>(c_grp)].c_str())); }
    if (c_ds >= 0)
    {
      eng_ds.push_back(parseDouble(f[static_cast<std::size_t>(c_ds)]));
      eng_ds_text.push_back(f[static_cast<std::size_t>(c_ds)]);
    }
    if (c_q >= 0) { eng_q.push_back(parseDouble(f[static_cast<std::size_t>(c_q)])); }
  }
  const std::size_t n = features.size();
  if (n == 0) { std::fprintf(stderr, "no rows\n"); return 1; }

  std::string group_source = "Group.Index column";
  if (c_grp < 0)
  {
    if (opt.count("-foldmap"))
    {
      // The FOLDMAP lists group_id per group INDEX, and the engine indexes groups by first
      // occurrence in row order -- this table's row order.
      std::ifstream fm(opt["-foldmap"]);
      std::string tok;
      bool found = false;
      while (fm >> tok) { if (tok == "FOLDMAP") { found = true; break; } }
      std::size_t n_map = 0;
      if (!found || !(fm >> n_map) || n_map != first_index.size())
      {
        std::fprintf(stderr, "-foldmap: no FOLDMAP, or it lists %zu groups against this table's %zu\n",
                     n_map, first_index.size());
        return 1;
      }
      std::vector<long long> gid(n_map);
      for (std::size_t g = 0; g < n_map; ++g)
      {
        int fd = 0;
        if (!(fm >> gid[g] >> fd)) { std::fprintf(stderr, "-foldmap: truncated\n"); return 1; }
      }
      group.resize(n);
      for (std::size_t i = 0; i < n; ++i) { group[i] = gid[static_cast<std::size_t>(first_occ[i])]; }
      group_source = "-foldmap " + opt["-foldmap"];
    }
    else
    {
      group = first_occ;
      group_source = "FIRST-OCCURRENCE index (no Group.Index, no -foldmap): the fold partition is "
                     "NOT the engine's";
    }
  }

  std::fprintf(stderr, "[replay] %zu rows, %zu groups, %zu sub-scores; groups from %s\n", n,
               first_index.size(), m, group_source.c_str());

  const auto scored = ODIA::Scoring::scoreSemiSupervisedLDA(features, labels, group, params, "", "", &names);

  // ---- report -----------------------------------------------------------------------------------
  // identifications: target groups whose best q <= 0.01, as the engine counts them
  auto ids_at = [&](const std::vector<double>& q) {
    std::unordered_map<long long, double> best;
    for (std::size_t i = 0; i < n; ++i)
    {
      if (labels[i] != 1) { continue; }
      auto it = best.emplace(group[i], q[i]).first;
      it->second = std::min(it->second, q[i]);
    }
    std::unordered_set<long long> s;
    for (const auto& kv : best) { if (kv.second <= 0.01) { s.insert(kv.first); } }
    return s;
  };
  const bool fdr_valid = scored.n_iterations_trained > 0;
  const auto rep_ids = ids_at(scored.qvalue);
  std::printf("replay: rows %zu groups %zu  iterations trained %d skipped %d  fdr_valid %d\n", n,
              first_index.size(), scored.n_iterations_trained, scored.n_iterations_skipped,
              fdr_valid ? 1 : 0);
  std::printf("replay: identified at q<=0.01: %zu\n", fdr_valid ? rep_ids.size() : std::size_t{0});
  if (!eng_ds.empty())
  {
    std::size_t bit_equal = 0, print_equal = 0;
    double max_abs = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
      const double a = eng_ds[i], b = scored.dscore[i];
      if (a == b || (std::isnan(a) && std::isnan(b))) { ++bit_equal; }
      if (printed6(b) == eng_ds_text[i] || printed6(b) == printed6(a)) { ++print_equal; }
      if (std::isfinite(a) && std::isfinite(b)) { max_abs = std::max(max_abs, std::abs(a - b)); }
    }
    std::printf("agreement: DScore bit-identical %zu / %zu rows; equal at 6-digit print precision %zu / %zu; "
                "max |delta| %.17g\n", bit_equal, n, print_equal, n, max_abs);
  }
  bool identical = !eng_ds.empty() && !eng_q.empty();
  if (!eng_ds.empty())
  {
    for (std::size_t i = 0; i < n && identical; ++i)
    {
      const double a = eng_ds[i], b = scored.dscore[i];
      identical = (a == b) || (std::isnan(a) && std::isnan(b));
    }
  }
  if (!eng_q.empty())
  {
    const auto eng_ids = ids_at(eng_q);
    std::size_t both = 0;
    for (const auto g : eng_ids) { both += rep_ids.count(g); }
    std::printf("agreement: q<=0.01 engine %zu replay %zu both %zu engine-only %zu replay-only %zu\n",
                eng_ids.size(), rep_ids.size(), both, eng_ids.size() - both, rep_ids.size() - both);
    identical = identical && both == eng_ids.size() && both == rep_ids.size();
  }
  if (has("-expect_identical"))
  {
    std::printf("expect_identical: %s\n", identical ? "PASS" : "FAIL");
    if (!identical) { return 4; }
  }

  if (opt.count("-out"))
  {
    std::ofstream out(opt["-out"]);
    if (!out) { std::fprintf(stderr, "cannot open %s\n", opt["-out"].c_str()); return 1; }
    out.imbue(std::locale::classic());
    out.precision(17);
    out << "Row\tPrecursor.Id\tDecoy\tGroup\tDScore\tQValue\tPEP\n";
    for (std::size_t i = 0; i < n; ++i)
    {
      out << i << '\t' << ids[i] << '\t' << (labels[i] == 1 ? 0 : 1) << '\t' << group[i] << '\t'
          << scored.dscore[i] << '\t' << scored.qvalue[i] << '\t' << scored.pep[i] << '\n';
    }
    if (!out) { std::fprintf(stderr, "write failed: %s\n", opt["-out"].c_str()); return 1; }
  }
  return fdr_valid ? 0 : 3;
}
