// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/RtRefiner.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string_view>

namespace ODIA
{
  namespace
  {
    constexpr const char* AAS = "ACDEFGHIKLMNPQRSTVWY";
    constexpr std::size_t N_AA = 20;
    /// composition(20) + length + charge + calibrated iRT + intercept
    constexpr std::size_t N_FEAT = N_AA + 4;

    /// Strip modifications: everything inside brackets, and the brackets.
    std::string stripMods(std::string_view s)
    {
      std::string out;
      int depth = 0;
      for (const char c : s)
      {
        if (c == '(' || c == '[' || c == '{') { ++depth; }
        else if (c == ')' || c == ']' || c == '}') { depth = std::max(0, depth - 1); }
        else if (depth == 0 && std::isalpha(static_cast<unsigned char>(c))) { out += c; }
      }
      return out;
    }

    void features(const std::string& seq, int charge, double irt, double* v)
    {
      for (std::size_t i = 0; i < N_FEAT; ++i) { v[i] = 0.0; }
      for (const char c : seq)
      {
        const char* p = std::strchr(AAS, c);
        if (p) { v[static_cast<std::size_t>(p - AAS)] += 1.0; }
      }
      v[N_AA + 0] = static_cast<double>(seq.size());
      v[N_AA + 1] = static_cast<double>(charge);
      v[N_AA + 2] = irt;
      v[N_AA + 3] = 1.0;
    }

    double clampCorrection(double c, double limit)
    {
      if (!std::isfinite(c)) { return 0.0; }
      return std::clamp(c, -limit, limit);
    }

    double sd(const std::vector<double>& r)
    {
      if (r.size() < 2) { return 0.0; }
      double m = 0.0;
      for (const double x : r) { m += x; }
      m /= static_cast<double>(r.size());
      double s = 0.0;
      for (const double x : r) { s += (x - m) * (x - m); }
      return std::sqrt(s / static_cast<double>(r.size() - 1));
    }

    /// Solve (X'X + lambda I) w = X'y by Gauss-Jordan with partial pivoting.
    /// N_FEAT is 24, so this is microseconds and a library would be overkill.
    bool solve(std::vector<double> A, std::vector<double> b, std::vector<double>& w)
    {
      const std::size_t n = b.size();
      for (std::size_t c = 0; c < n; ++c)
      {
        std::size_t piv = c;
        for (std::size_t r = c + 1; r < n; ++r)
        {
          if (std::fabs(A[r * n + c]) > std::fabs(A[piv * n + c])) { piv = r; }
        }
        if (std::fabs(A[piv * n + c]) < 1e-12) { return false; }
        if (piv != c)
        {
          for (std::size_t k = 0; k < n; ++k) { std::swap(A[c * n + k], A[piv * n + k]); }
          std::swap(b[c], b[piv]);
        }
        const double d = A[c * n + c];
        for (std::size_t k = 0; k < n; ++k) { A[c * n + k] /= d; }
        b[c] /= d;
        for (std::size_t r = 0; r < n; ++r)
        {
          if (r == c) { continue; }
          const double f = A[r * n + c];
          if (f == 0.0) { continue; }
          for (std::size_t k = 0; k < n; ++k) { A[r * n + k] -= f * A[c * n + k]; }
          b[r] -= f * b[c];
        }
      }
      w = std::move(b);
      return true;
    }
  } // namespace

  RtRefiner::Report RtRefiner::fit(const std::vector<std::string>& sequences,
                                   const std::vector<int>& charges,
                                   const std::vector<double>& calibrated_irt,
                                   const std::vector<double>& observed_rt,
                                   const Options& options)
  {
    // Same reason as load(): every early return below leaves the object as it
    // was, so a failed refit after a successful one would keep applying the old
    // model while reporting that it was not fitted.
    fitted_ = false;
    weights_.clear();
    max_shift_ = 0.0;

    Report rep;
    const std::size_t n = sequences.size();
    if (n != charges.size() || n != calibrated_irt.size() || n != observed_rt.size())
    {
      rep.note = "input vectors differ in length";
      return rep;
    }
    rep.anchors = n;
    if (n < options.min_anchors)
    {
      std::ostringstream os;
      os << "only " << n << " anchors, below the " << options.min_anchors
         << " a per-run retention-time model needs";
      rep.note = os.str();
      return rep;
    }

    // Held out BY STRIPPED SEQUENCE. A peptide's charge states co-elute, so a
    // split that puts 2+ in train and 3+ in test reports memorisation.
    std::vector<std::string> stripped(n);
    for (std::size_t i = 0; i < n; ++i) { stripped[i] = stripMods(sequences[i]); }
    // FNV-1a rather than zlib's crc32: the split only has to be deterministic
    // and INDEPENDENT OF SCORE, and linking zlib into odia_core for a bucket
    // hash would be a dependency bought for nothing.
    const auto bucket = [](const std::string& s) {
      std::uint64_t h = 1469598103934665603ULL;
      for (const unsigned char c : s) { h = (h ^ c) * 1099511628211ULL; }
      return static_cast<unsigned>(h % 1000ULL);
    };
    const unsigned cut = static_cast<unsigned>(options.holdout * 1000.0);
    const double med_abs_limit = options.trim_mads;
    // The clamp, in the units of the axis actually in use.
    double lo = calibrated_irt[0], hi = calibrated_irt[0];
    for (const double v : calibrated_irt) { lo = std::min(lo, v); hi = std::max(hi, v); }
    const double max_shift = std::max(1e-6, options.max_shift_fraction * (hi - lo));

    std::vector<double> X(n * N_FEAT);
    for (std::size_t i = 0; i < n; ++i)
    {
      features(stripped[i], charges[i], calibrated_irt[i], &X[i * N_FEAT]);
    }

    // TRIM before fitting. Pass 1's anchors are peak groups accepted at
    // -anchor_q, and their residual tail reaches 1,600 s -- misidentifications,
    // not chromatography. A squared-loss fit chases them, which is how a model
    // that improved the held-out SD by 2% still destroyed the search.
    std::vector<double> resid(n);
    for (std::size_t i = 0; i < n; ++i) { resid[i] = observed_rt[i] - calibrated_irt[i]; }
    std::vector<double> sorted = resid;
    std::sort(sorted.begin(), sorted.end());
    const double med = sorted[n / 2];
    std::vector<double> ad(n);
    for (std::size_t i = 0; i < n; ++i) { ad[i] = std::fabs(resid[i] - med); }
    std::sort(ad.begin(), ad.end());
    const double mad = 1.4826 * ad[n / 2];
    const double keep = med_abs_limit * std::max(mad, 1.0);

    std::vector<std::size_t> train, test;
    std::size_t trimmed = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
      if (std::fabs(resid[i] - med) > keep) { ++trimmed; continue; }
      (bucket(stripped[i]) < cut ? test : train).push_back(i);
    }
    if (train.size() < options.min_anchors || test.size() < 30)
    {
      rep.note = "split left too few rows on one side";
      return rep;
    }
    rep.trimmed = trimmed;
    rep.train = train.size();
    rep.held_out = test.size();

    // Standardise on the TRAINING rows only. Composition counts and retention
    // times differ by three orders of magnitude, and an unstandardised ridge
    // penalty would fall almost entirely on the composition columns.
    std::vector<double> mu(N_FEAT, 0.0), sg(N_FEAT, 0.0);
    for (const std::size_t i : train)
    {
      for (std::size_t f = 0; f < N_FEAT; ++f) { mu[f] += X[i * N_FEAT + f]; }
    }
    for (std::size_t f = 0; f < N_FEAT; ++f) { mu[f] /= static_cast<double>(train.size()); }
    for (const std::size_t i : train)
    {
      for (std::size_t f = 0; f < N_FEAT; ++f)
      {
        const double d = X[i * N_FEAT + f] - mu[f];
        sg[f] += d * d;
      }
    }
    for (std::size_t f = 0; f < N_FEAT; ++f)
    {
      sg[f] = std::sqrt(sg[f] / static_cast<double>(train.size())) + 1e-9;
    }

    std::vector<double> A(N_FEAT * N_FEAT, 0.0), b(N_FEAT, 0.0);
    for (const std::size_t i : train)
    {
      double z[N_FEAT];
      for (std::size_t f = 0; f < N_FEAT; ++f) { z[f] = (X[i * N_FEAT + f] - mu[f]) / sg[f]; }
      for (std::size_t r = 0; r < N_FEAT; ++r)
      {
        // THE RESIDUAL, not the absolute retention time.
        //
        // Fitting observed_rt directly makes this an unbounded linear
        // extrapolator: a precursor whose calibrated iRT falls outside the
        // training range gets a prediction with nothing holding it near the
        // gradient, and pass 2 then extracts from a time the peptide cannot be
        // at. The monotone map it replaced was bounded by interpolation.
        // Measured: IH1 went from 1,464 identifications to ZERO.
        b[r] += z[r] * (observed_rt[i] - calibrated_irt[i]);
        for (std::size_t c = 0; c < N_FEAT; ++c) { A[r * N_FEAT + c] += z[r] * z[c]; }
      }
    }
    for (std::size_t f = 0; f < N_FEAT; ++f) { A[f * N_FEAT + f] += options.ridge; }

    std::vector<double> w;
    if (!solve(A, b, w)) { rep.note = "the normal equations are singular"; return rep; }

    std::vector<double> before, after;
    before.reserve(test.size());
    after.reserve(test.size());
    for (const std::size_t i : test)
    {
      double z[N_FEAT], pred = 0.0;
      for (std::size_t f = 0; f < N_FEAT; ++f) { z[f] = (X[i * N_FEAT + f] - mu[f]) / sg[f]; }
      for (std::size_t f = 0; f < N_FEAT; ++f) { pred += w[f] * z[f]; }
      pred = clampCorrection(pred, max_shift);
      before.push_back(observed_rt[i] - calibrated_irt[i]);
      after.push_back(observed_rt[i] - (calibrated_irt[i] + pred));
    }
    rep.sd_before = sd(before);
    rep.sd_after = sd(after);

    // REFUSED WHEN IT DOES NOT HELP, on the run's own held-out rows. A per-run
    // model that is worse than the calibration it sits on must not be applied,
    // and the only way to know is to measure it here rather than to assume it
    // from another file.
    if (!(rep.sd_after < rep.sd_before))
    {
      std::ostringstream os;
      os << "refinement did not improve the held-out residual (" << rep.sd_after
         << " s against " << rep.sd_before << " s); keeping the calibrated axis";
      rep.note = os.str();
      return rep;
    }

    // Fold the standardisation into the weights so `apply` is a dot product.
    weights_.assign(N_FEAT + 1, 0.0);
    double bias = 0.0;
    for (std::size_t f = 0; f < N_FEAT; ++f)
    {
      weights_[f] = w[f] / sg[f];
      bias -= w[f] * mu[f] / sg[f];
    }
    weights_[N_FEAT] = bias;
    max_shift_ = max_shift;
    fitted_ = true;
    rep.fitted = true;
    return rep;
  }

  RtRefiner::Report RtRefiner::fit(const std::vector<std::string>& sequences,
                                   const std::vector<int>& charges,
                                   const std::vector<double>& calibrated_irt,
                                   const std::vector<double>& observed_rt)
  {
    return fit(sequences, charges, calibrated_irt, observed_rt, Options{});
  }

  bool RtRefiner::save(const std::string& path, const std::string& provenance) const
  {
    if (!fitted_) { return false; }
    std::ofstream os(path);
    if (!os) { return false; }
    os << "# ODIA retention-time refinement model\n"
       << "# WARNING: fitted on ONE run's chromatography. Applying it to a run on a\n"
       << "# different gradient, instrument or method is the failure that cost 2,027\n"
       << "# precursors when it last happened. A series sharing all three is the case\n"
       << "# this exists for.\n"
       << "# provenance: " << provenance << "\n"
       << "version\t2\n"
       // THE CLAMP TRAVELS WITH THE WEIGHTS. Version 1 wrote only the weights,
       // so a model fitted with a 10-unit clamp was applied after loading with
       // the member default of 120 -- which on a normalised iRT axis is nearly
       // the whole gradient, and defeats the exact protection added after a
       // 1,464 -> 0 incident. A file without it is refused rather than guessed.
       << "max_shift\t" << max_shift_ << "\n"
       << "features\t" << (N_FEAT + 1) << "\n";
    os.precision(17);
    for (const double w : weights_) { os << w << "\n"; }
    return static_cast<bool>(os);
  }

  bool RtRefiner::load(const std::string& path, std::string* provenance_out)
  {
    // RESET FIRST. Neither this nor fit() cleared its own state, so a failed
    // load after a successful fit left the object `fitted_` with the old
    // weights, and two loads appended one file's weights to the other's.
    fitted_ = false;
    weights_.clear();
    max_shift_ = 0.0;

    std::ifstream is(path);
    if (!is) { return false; }
    std::string line, prov;
    std::size_t n = 0;
    int version = 0;
    double shift = 0.0;
    while (std::getline(is, line))
    {
      if (!line.empty() && line[0] == '#')
      {
        // An EXPLICIT marker. The first version matched lines containing
        // "fitted" or "anchors", so provenance written in any other words was
        // silently dropped -- and provenance is the whole safety mechanism for
        // reusing a model across runs. A guard that quietly discards its own
        // evidence is worse than no guard.
        static constexpr const char* kMark = "# provenance: ";
        if (line.rfind(kMark, 0) == 0) { prov = line.substr(std::strlen(kMark)); }
        continue;
      }
      if (line.rfind("features", 0) == 0)
      { std::istringstream(line.substr(8)) >> n; continue; }
      if (line.rfind("version", 0) == 0)
      { std::istringstream(line.substr(7)) >> version; continue; }
      if (line.rfind("max_shift", 0) == 0)
      { std::istringstream(line.substr(9)) >> shift; continue; }
      if (line.empty()) { continue; }
      weights_.push_back(std::stod(line));
    }
    // The weight count is the contract. A file written by a different feature
    // set would otherwise be read as a valid model and silently predict noise.
    // A version-1 file carries no clamp, and inventing one is how the
    // protection gets silently disabled. Refuse it.
    if (version < 2 || !(shift > 0.0) || n != N_FEAT + 1 || weights_.size() != N_FEAT + 1)
    {
      weights_.clear();
      return false;
    }
    max_shift_ = shift;
    if (provenance_out) { *provenance_out = prov; }
    fitted_ = true;
    return true;
  }

  std::size_t RtRefiner::apply(Library& library) const
  {
    if (!fitted_) { return 0; }
    auto& p = library.precursors();
    std::size_t changed = 0;
    for (std::size_t i = 0; i < library.precursorCount(); ++i)
    {
      if (!std::isfinite(p.irt[i])) { continue; }
      const std::string seq = stripMods(library.strings().get(p.modified_sequence[i]));
      double v[N_FEAT];
      features(seq, static_cast<int>(p.charge[i]), static_cast<double>(p.irt[i]), v);
      double corr = weights_[N_FEAT];
      for (std::size_t f = 0; f < N_FEAT; ++f) { corr += weights_[f] * v[f]; }
      // CLAMPED. The model may not move a precursor further than a real
      // chromatographic shift, however confident its extrapolation is.
      corr = clampCorrection(corr, max_shift_);
      p.irt[i] = static_cast<float>(static_cast<double>(p.irt[i]) + corr);
      ++changed;
    }
    return changed;
  }

} // namespace ODIA
