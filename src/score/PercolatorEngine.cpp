// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/scoring/PercolatorEngine.h>

#include <OpenMS/ANALYSIS/ID/Percolator.h>
#include <OpenMS/ANALYSIS/ID/PercolatorTypes.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <unordered_map>

namespace ODIA
{
  namespace Scoring
  {
    ScoredGroups scorePercolator(const std::vector<std::vector<double>>& features,
                                 const std::vector<int>& labels,
                                 const std::vector<long long>& group,
                                 const LDAParams& params,
                                 const std::vector<std::string>& feature_names,
                                 std::string* diagnostic)
    {
      ScoredGroups out;
      const std::size_t n = features.size();
      out.dscore.assign(n, 0.0);
      out.qvalue.assign(n, 1.0);
      out.pvalue.assign(n, 1.0);
      out.pep.assign(n, 1.0);
      if (n == 0 || labels.size() != n || group.size() != n) { return out; }

      OpenMS::RescoreInput in;
      in.features.reserve(n);
      in.is_decoy.reserve(n);
      in.cv_group_keys.reserve(n);

      // NaN is legal in our sub-scores and means "not measured" -- IM_DELTA on a
      // run without mobility, MS1_COELUTION without MS1. Percolator has no such
      // convention, and a NaN would propagate silently through an SVM dot
      // product and poison every score. Replaced by the column's own median,
      // which is the least informative value it can take rather than a zero
      // that the learner would read as a measurement.
      const std::size_t m = features[0].size();
      std::vector<double> fill(m, 0.0);
      for (std::size_t j = 0; j < m; ++j)
      {
        std::vector<double> col;
        col.reserve(n);
        for (std::size_t i = 0; i < n; ++i)
        {
          if (j < features[i].size() && std::isfinite(features[i][j]))
          { col.push_back(features[i][j]); }
        }
        if (!col.empty())
        {
          const std::size_t h = col.size() / 2;
          std::nth_element(col.begin(), col.begin() + h, col.end());
          fill[j] = col[h];
        }
      }

      // Precursor -> a small dense fold key. Percolator keeps rows sharing a key
      // in one fold; the candidate peak groups of one precursor share a
      // chromatogram, so splitting them across folds would leak.
      std::unordered_map<long long, int> key_of;
      key_of.reserve(n * 2);
      for (std::size_t i = 0; i < n; ++i)
      {
        std::vector<double> row(m, 0.0);
        for (std::size_t j = 0; j < m; ++j)
        {
          const double v = (j < features[i].size()) ? features[i][j] : fill[j];
          row[j] = std::isfinite(v) ? v : fill[j];
        }
        in.features.push_back(std::move(row));
        in.is_decoy.push_back(labels[i] != 1);
        const auto it = key_of.emplace(group[i], static_cast<int>(key_of.size()));
        in.cv_group_keys.push_back(it.first->second);
      }
      in.feature_names = feature_names;
      if (in.feature_names.size() != m) { in.feature_names.assign(m, "var"); }

      try
      {
        OpenMS::Percolator perc;
        const OpenMS::RescoreOutput res = perc.rescore(in);
        if (res.scores.size() != n) { return out; }
        for (std::size_t i = 0; i < n; ++i)
        {
          out.dscore[i] = res.scores[i];
          out.qvalue[i] = i < res.q_values.size() ? res.q_values[i] : 1.0;
          out.pep[i] = i < res.peps.size() ? res.peps[i] : 1.0;
        }
        out.n_iterations_trained = 1;
        if (diagnostic)
        {
          std::size_t called = 0;
          for (std::size_t i = 0; i < n; ++i)
          { if (labels[i] == 1 && out.qvalue[i] <= 0.01) { ++called; } }
          std::ostringstream o;
          o << "percolator: " << n << " rows over " << key_of.size()
            << " cross-validation groups, " << m << " features, pi0 " << perc.getPi0()
            << ", " << called << " target rows at q<=0.01";
          *diagnostic = o.str();
        }
      }
      catch (const std::exception& e)
      {
        // A failed engine must not look like a scored run of zero. The caller's
        // "fitted 0 iterations" warning is what surfaces this.
        if (diagnostic) { *diagnostic = std::string("percolator FAILED: ") + e.what(); }
        out.n_iterations_skipped = 1;
      }
      (void)params;
      return out;
    }
  } // namespace Scoring
} // namespace ODIA
