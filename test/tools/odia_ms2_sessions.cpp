// Spreading batches across sessions must not change a single bit. The chunking
// is built in serial order precisely so it cannot, but "cannot" is a claim
// about code that has to be checked against the model's actual output.
#include <odia/PeptDeepPredictor.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
  if (argc < 3)
  {
    std::fprintf(stderr, "usage: odia_ms2_sessions <model.onnx> <library.tsv> [peptides]\n");
    return 2;
  }
  const std::size_t n = argc > 3 ? std::stoul(argv[3]) : 3000;
  std::vector<OpenMS::AASequence> peptides;
  std::vector<int> charges;
  std::ifstream in(argv[2]);
  if (!in)
  {
    std::fprintf(stderr, "cannot read %s\n", argv[2]);
    return 2;
  }
  std::string line, last;
  std::getline(in, line);
  while (std::getline(in, line) && peptides.size() < n)
  {
    const auto a = line.find('\t');
    const auto b = line.find('\t', a + 1);
    const auto c = line.find('\t', b + 1);
    auto s = line.substr(a + 1, b - a - 1);
    if (s == last) { continue; }
    last = s;
    try { peptides.push_back(OpenMS::AASequence::fromString(s)); } catch (...) { continue; }
    charges.push_back(std::atoi(line.substr(b + 1, c - b - 1).c_str()));
  }

  auto run = [&](int sessions)
  {
    ODIA::PeptDeepPredictor p(argv[1], false, 1, sessions);
    std::vector<ODIA::PeptDeepPredictor::Failure> f;
    const auto t0 = std::chrono::steady_clock::now();
    auto out = p.predictMS2(peptides, charges, 30.0f, "timsTOF", &f);
    const double dt = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - t0).count();
    std::printf("%2d session(s): %6.2f s  %7.1f peptides/s  (%zu sessions live)\n",
                sessions, dt, peptides.size() / dt, p.sessionCount());
    return out;
  };

  if (peptides.size() < 8)
  {
    std::fprintf(stderr, "only %zu peptides read from %s; need at least 8\n",
                 peptides.size(), argv[2]);
    return 2;
  }

  const auto reference = run(1);
  for (const int s : {2, 8, 32})
  {
    const auto other = run(s);
    if (other.size() != reference.size()) { std::printf("  SIZE MISMATCH\n"); return 1; }
    std::size_t differing = 0, compared = 0;
    for (std::size_t i = 0; i < reference.size(); ++i)
    {
      if (other[i].intensities.size() != reference[i].intensities.size())
      {
        std::printf("  length mismatch at %zu\n", i);
        return 1;
      }
      for (std::size_t k = 0; k < reference[i].intensities.size(); ++k)
      {
        ++compared;
        // Bit-for-bit, not a tolerance: any difference at all means the split
        // is leaking into the arithmetic.
        if (other[i].intensities[k] != reference[i].intensities[k]) { ++differing; }
      }
    }
    std::printf("  vs 1 session: %zu of %zu values differ%s\n", differing, compared,
                differing ? "  <-- NOT IDENTICAL" : "  (identical)");
    if (differing) { return 1; }
  }
  return 0;
}
