#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <random>
#include <vector>
// Gate C at REAL dimensions, against Gate A and Gate B on identical traces.
namespace {
constexpr int J = 12, T = 200, HALF = 2;

std::vector<std::vector<double>> traces(std::mt19937& rng, bool peak, double lvl)
{
  std::poisson_distribution<int> pois(lvl);          // counting noise
  std::vector<std::vector<double>> v(J, std::vector<double>(T, 0.0));
  for (int j = 0; j < J; ++j)
    for (int i = 0; i < T; ++i) v[j][i] = pois(rng);
  if (peak) {                                        // a co-eluting peak
    const int c = T / 2;
    for (int j = 0; j < J; ++j)
      for (int i = c - 3; i <= c + 3; ++i) {
        const double d = i - c;
        v[j][i] += 6.0 * lvl * std::exp(-d * d / 4.0);
      }
  }
  return v;
}
double madOf(std::vector<double> a)
{
  std::sort(a.begin(), a.end());
  const double m = a[a.size()/2];
  for (auto& x : a) x = std::abs(x - m);
  std::sort(a.begin(), a.end());
  return a[a.size()/2];
}
double medOf(std::vector<double> a)
{ std::sort(a.begin(), a.end()); return a[a.size()/2]; }

// Gate C
double gateC(const std::vector<std::vector<double>>& v)
{
  std::vector<double> s(T, 0.0);
  for (int j = 0; j < J; ++j) {
    std::vector<double> y(T);
    for (int i = 0; i < T; ++i) y[i] = std::sqrt(std::max(0.0, v[j][i]));
    const double med = medOf(y), mad = madOf(y);
    if (!(mad > 0.0)) continue;
    const double sc = 1.0 / (1.4826 * mad);
    for (int i = 0; i < T; ++i) s[i] += (y[i] - med) * sc;
  }
  double best = 0.0;
  for (int i = 0; i < T; ++i) {
    double acc = 0.0; int used = 0;
    for (int j = std::max(0, i-HALF); j <= std::min(T-1, i+HALF); ++j) { acc += s[j]; ++used; }
    best = std::max(best, acc / (2*HALF+1));
  }
  return best;
}
// Gate A: sum of z over everything, pass if > 0
bool gateA(const std::vector<std::vector<double>>& v)
{
  double sum = 0.0;
  for (int j = 0; j < J; ++j) {
    const double med = medOf(v[j]), mad = madOf(v[j]);
    if (!(mad > 0.0)) continue;
    const double sc = 1.0 / (1.4826 * mad);
    for (int i = 0; i < T; ++i) sum += (v[j][i] - med) * sc;
  }
  return sum > 0.0;
}
// Gate B: >=2 transitions with a 3-sigma excursion anywhere
bool gateB(const std::vector<std::vector<double>>& v)
{
  int fired = 0;
  for (int j = 0; j < J; ++j) {
    const double med = medOf(v[j]), mad = madOf(v[j]);
    if (!(mad > 0.0)) continue;
    const double sc = 1.0 / (1.4826 * mad);
    double mx = -1e9;
    for (int i = 0; i < T; ++i) mx = std::max(mx, (v[j][i] - med) * sc);
    if (mx >= 3.0) ++fired;
  }
  return fired >= 2;
}
}
int main()
{
  std::mt19937 rng(20260814);
  const int TR = 3000;
  // Calibrate Gate C's tau on the NULL, as Kimi specifies -- 95th percentile.
  std::vector<double> null_stats;
  for (int t = 0; t < TR; ++t) null_stats.push_back(gateC(traces(rng, false, 20.0)));
  std::sort(null_stats.begin(), null_stats.end());
  const double tau = null_stats[std::size_t(0.95 * null_stats.size())];
  std::printf("Gate C tau calibrated on the decoy null (alpha=0.05): %.3f\n\n", tau);

  int aN=0,bN=0,cN=0,aP=0,bP=0,cP=0;
  for (int t = 0; t < TR; ++t) {
    auto n = traces(rng, false, 20.0);
    if (gateA(n)) ++aN; if (gateB(n)) ++bN; if (gateC(n) >= tau) ++cN;
    auto p = traces(rng, true, 20.0);
    if (gateA(p)) ++aP; if (gateB(p)) ++bP; if (gateC(p) >= tau) ++cP;
  }
  std::printf("at REAL dimensions (J=12 traces, T=200 cycles), Poisson noise:\n");
  std::printf("  %-8s %12s %12s\n", "gate", "ABSENT", "PRESENT");
  std::printf("  %-8s %11.1f%% %11.1f%%\n", "A sum>0", 100.0*aN/TR, 100.0*aP/TR);
  std::printf("  %-8s %11.1f%% %11.1f%%\n", "B 3sig",  100.0*bN/TR, 100.0*bP/TR);
  std::printf("  %-8s %11.1f%% %11.1f%%\n", "C coel",  100.0*cN/TR, 100.0*cP/TR);
}
