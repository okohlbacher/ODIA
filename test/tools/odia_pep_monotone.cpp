// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

/// PEP must be low where the score is high.
///
/// `assignQValues` sorts high->low by score, so the posterior error probability
/// must be NON-DECREASING IN INDEX. The monotonicity sweep ran the other way:
/// from the low-score end with a running max, which gave every entry the
/// maximum PEP of everything scoring BELOW it. On a clean list that is not a
/// small distortion -- [0.01, 0.02, 0.5, 0.9] came out [0.9, 0.9, 0.9, 0.9],
/// reporting the run's best precursor with the run's worst PEP.
///
/// Two properties are pinned, because only checking monotonicity would have
/// passed the broken version too -- it is monotone, in the wrong direction:
///   * PEP is non-decreasing from the best score to the worst;
///   * the best-scoring group's PEP is no worse than the worst-scoring group's.

#include <odia/scoring/lda.h>

#include <cstdlib>
#include <iostream>
#include <vector>

int main()
{
  using ODIA::Scoring::lda_detail::RankedGroup;

  // A clean separation: targets on top, decoys below, so the local ratio at the
  // high-score end is genuinely near zero and a correct sweep must keep it there.
  std::vector<RankedGroup> ranked;
  for (int i = 0; i < 400; ++i)
  {
    RankedGroup g;
    g.group_index = static_cast<std::size_t>(i);
    g.best_row = static_cast<std::size_t>(i);
    g.label = i < 200 ? 1 : -1;
    g.score = 10.0 - 0.02 * i;      // strictly decreasing
    ranked.push_back(g);
  }
  // Interleave a few decoys into the confident region, so the raw local ratio is
  // noisy and the sweep has something to do.
  for (int i : {17, 43, 88}) { ranked[static_cast<std::size_t>(i)].label = -1; }

  ODIA::Scoring::lda_detail::assignQValues(ranked, /*use_pi0=*/false);

  double prev = -1.0;
  for (std::size_t i = 0; i < ranked.size(); ++i)
  {
    if (ranked[i].pep < prev - 1e-12)
    {
      std::cerr << "PEP decreased with falling score at index " << i << ": "
                << prev << " -> " << ranked[i].pep << "\n";
      return EXIT_FAILURE;
    }
    prev = ranked[i].pep;
  }

  const double best = ranked.front().pep;
  const double worst = ranked.back().pep;
  if (!(best <= worst))
  {
    std::cerr << "best-scoring group has PEP " << best
              << ", worse than the worst-scoring group's " << worst << "\n";
    return EXIT_FAILURE;
  }
  // The broken sweep collapsed the whole column onto one value. A correct one
  // must retain range on a list this cleanly separated.
  if (!(worst - best > 1e-6))
  {
    std::cerr << "PEP column is constant (" << best
              << "): the monotonicity sweep flattened it\n";
    return EXIT_FAILURE;
  }

  std::cout << "PEP rises " << best << " -> " << worst
            << " from best to worst score over " << ranked.size() << " groups\n";
  return EXIT_SUCCESS;
}
