// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/TofCalibration.h>

#include <odia/ArrowColumn.h>
#include <odia/ZipArchive.h>

#include <arrow/api.h>
#include <parquet/arrow/reader.h>

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace ODIA
{

  namespace
  {
    bool usable(double a, double b) { return std::isfinite(a) && std::isfinite(b) && b > 0.0; }

    [[noreturn]] void fail(const std::string& what)
    {
      throw std::runtime_error("per-frame TOF calibration: " + what);
    }
  }

  void remapTofMz(std::vector<double>& mz, double a, double b, double c0, double c1,
                  TofRemapStats* stats)
  {
    if (!usable(a, b)) { fail("run-wide chord is not usable"); }
    if (!usable(c0, c1)) { fail("per-frame coefficients are not usable"); }
    double worst = 0.0;
    for (double& m : mz)
    {
      // The chord root a + b*tof is positive for every acquired index, so the
      // positive square root is the one the reader squared.
      const double tof = (std::sqrt(m) - a) / b;
      const double t = std::nearbyint(tof);
      const double r = std::abs(tof - t);
      if (!(r <= kTofIndexTolerance))
      {
        std::ostringstream os;
        os.precision(17);
        os << "m/z " << m << " does not invert to an integer TOF index on the chord a=" << a
           << " b=" << b << " (index " << tof << ", residual " << r
           << "); the array was not decoded with this chord";
        fail(os.str());
      }
      if (r > worst) { worst = r; }
      const double root = c0 + c1 * t;
      m = root * root;
    }
    if (stats)
    {
      stats->peaks += mz.size();
      if (worst > stats->max_index_residual) { stats->max_index_residual = worst; }
    }
  }

  double tofShiftPpm(double mz, double a, double b, double c0, double c1)
  {
    const double t = (std::sqrt(mz) - a) / b;
    const double root = c0 + c1 * t;
    return (root * root / mz - 1.0) * 1e6;
  }

  TofCoefficientTable TofCoefficientTable::fromColumns(const std::vector<std::uint64_t>& index,
                                                       const std::vector<double>& c0,
                                                       const std::vector<double>& c1,
                                                       const std::vector<bool>& c0_null,
                                                       const std::vector<bool>& c1_null,
                                                       std::size_t n_spectra)
  {
    const std::size_t rows = index.size();
    if (c0.size() != rows || c1.size() != rows || c0_null.size() != rows || c1_null.size() != rows)
    {
      fail("coefficient columns differ in length");
    }
    TofCoefficientTable t;
    t.c0_.assign(n_spectra, std::numeric_limits<double>::quiet_NaN());
    t.c1_.assign(n_spectra, std::numeric_limits<double>::quiet_NaN());
    t.has_.assign(n_spectra, false);
    std::vector<bool> seen(n_spectra, false);
    for (std::size_t r = 0; r < rows; ++r)
    {
      const std::uint64_t i = index[r];
      if (i >= n_spectra)
      {
        fail("row " + std::to_string(r) + " names spectrum " + std::to_string(i) + " of " +
             std::to_string(n_spectra));
      }
      if (c0_null[r] != c1_null[r])
      {
        fail("spectrum " + std::to_string(i) + " carries only one of tof_c0/tof_c1");
      }
      const bool pair = !c0_null[r];
      if (pair && !usable(c0[r], c1[r]))
      {
        fail("spectrum " + std::to_string(i) + " carries unusable coefficients");
      }
      if (seen[i])
      {
        // The same spectrum on two rows must say the same thing, or which one
        // wins would depend on row order.
        const bool same = (t.has_[i] == pair) && (!pair || (t.c0_[i] == c0[r] && t.c1_[i] == c1[r]));
        if (!same) { fail("spectrum " + std::to_string(i) + " has conflicting coefficient rows"); }
        continue;
      }
      seen[i] = true;
      if (pair)
      {
        t.c0_[i] = c0[r];
        t.c1_[i] = c1[r];
        t.has_[i] = true;
        ++t.with_pair_;
      }
    }
    return t;
  }

  TofCoefficientTable TofCoefficientTable::fromArchive(const std::string& mzpeak_path,
                                                       std::size_t n_spectra)
  {
    static const std::string entry = "spectra_metadata.parquet";
    static const char* kC0 = "opt_MZP_1000003_tof_c0";
    static const char* kC1 = "opt_MZP_1000004_tof_c1";

    ZipArchive zip(mzpeak_path);
    if (zip.count(entry) != 1)
    {
      fail(mzpeak_path + " holds " + std::to_string(zip.count(entry)) + " '" + entry + "' entries, need 1");
    }
    auto reader_result = parquet::arrow::OpenFile(zip.open(entry), arrow::default_memory_pool());
    if (!reader_result.ok()) { fail("cannot read " + entry + ": " + reader_result.status().ToString()); }
    std::unique_ptr<parquet::arrow::FileReader> reader = std::move(*reader_result);
    std::shared_ptr<arrow::Schema> schema;
    if (!reader->GetSchema(&schema).ok()) { fail("cannot read the schema of " + entry); }

    // ReadColumn takes a TOP-LEVEL field index, which is what GetFieldIndex
    // returns -- correct even though other columns of this table are nested.
    auto column = [&](const char* name) {
      const int f = schema->GetFieldIndex(name);
      if (f < 0)
      {
        fail(entry + " has no column '" + name +
             "' (not an ims-compact archive, or a converter that does not store the per-frame pair)");
      }
      std::shared_ptr<arrow::ChunkedArray> out;
      const auto st = reader->ReadColumn(f, &out);
      if (!st.ok()) { fail("cannot read column '" + std::string(name) + "': " + st.ToString()); }
      ChunkedColumn c(out);
      c.setName(name);
      return c;
    };
    ChunkedColumn idx = column("index");
    ChunkedColumn c0 = column(kC0);
    ChunkedColumn c1 = column(kC1);
    const std::int64_t rows = idx.length();
    if (c0.length() != rows || c1.length() != rows) { fail("coefficient columns differ in length"); }

    std::vector<std::uint64_t> vi(rows);
    std::vector<double> v0(rows), v1(rows);
    std::vector<bool> n0(rows), n1(rows);
    for (std::int64_t r = 0; r < rows; ++r)
    {
      if (idx.isNull(r)) { fail("row " + std::to_string(r) + " has a NULL spectrum index"); }
      const std::int64_t i = idx.getInt64(r, -1);
      if (i < 0) { fail("row " + std::to_string(r) + " has a negative spectrum index"); }
      vi[r] = static_cast<std::uint64_t>(i);
      n0[r] = c0.isNull(r);
      n1[r] = c1.isNull(r);
      v0[r] = c0.getDouble(r, std::numeric_limits<double>::quiet_NaN());
      v1[r] = c1.getDouble(r, std::numeric_limits<double>::quiet_NaN());
    }
    return fromColumns(vi, v0, v1, n0, n1, n_spectra);
  }

} // namespace ODIA
