// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/DIANNLibraryFile.h>

#include <arrow/api.h>
#include <arrow/compute/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>

#include <charconv>
#include <cstring>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace ODIA
{
  namespace
  {
    double toDouble(std::string_view s)
    {
      double v = 0.0;
      if (s.empty()) { return std::nan(""); }
      // from_chars for floating point is not universally available with the
      // required precision across libstdc++ versions in this project's range,
      // so parse through strtod on a NUL-terminated copy of the field.
      char buf[64];
      const std::size_t n = std::min(s.size(), sizeof(buf) - 1);
      std::memcpy(buf, s.data(), n);
      buf[n] = '\0';
      char* end = nullptr;
      v = std::strtod(buf, &end);
      return end == buf ? std::nan("") : v;
    }

    long toLong(std::string_view s)
    {
      long v = 0;
      const auto* first = s.data();
      const auto* last = s.data() + s.size();
      while (first != last && (*first == ' ' || *first == '+')) { ++first; }
      if (std::from_chars(first, last, v).ec != std::errc{}) { return 0; }
      return v;
    }

    /// One parsed row, in the order the builder needs it.
    struct Row
    {
      std::string_view precursor_id;
      std::string_view modified_sequence;
      std::string_view protein_group;
      double precursor_mz = 0, product_mz = 0;
      double rt = 0, im = 0, intensity = 0;
      long precursor_charge = 0, fragment_charge = 0, ordinal = 0, decoy = 0;
      std::string_view fragment_type, loss_type;
    };

    /// Appends rows to the library, starting a new precursor when the
    /// Precursor.Id changes.
    class Builder
    {
    public:
      explicit Builder(Library& lib) : lib_(lib) {}

      void append(const Row& r)
      {
        if (!have_current_ || r.precursor_id != current_id_ || (r.decoy != 0) != current_decoy_)
        {
          closeCurrent();
          // A precursor's rows must be contiguous, because grouping keys off a
          // change of Precursor.Id. If they are not, the same id reappears and
          // one precursor silently becomes several single-transition ones with
          // a plausible-looking count -- which is what any library that has been
          // sorted by m/z or concatenated will do.
          // Hash rather than store the id: a set<pair<string,bool>> measured
          // 128 B per precursor, ~0.9 GB at the scale D3 is written against.
          std::size_t key = std::hash<std::string_view>{}(r.precursor_id);
          key = key * 31 + (r.decoy != 0 ? 1u : 0u);
          if (!seen_.insert(key).second) { ++reopened_; }
          startPrecursor(r);
        }
        auto& t = lib_.transitions();
        t.product_mz.push_back(toFixed(r.product_mz));
        t.library_intensity.push_back(static_cast<float>(r.intensity));
        t.type.push_back(parseFragmentType(r.fragment_type));
        t.ordinal.push_back(static_cast<std::uint8_t>(r.ordinal));
        t.charge.push_back(static_cast<std::int8_t>(r.fragment_charge));
        t.loss.push_back(parseLossType(r.loss_type));
        ++count_;
      }

      void finish() { closeCurrent(); }

      /// Number of precursors whose rows were not contiguous.
      std::size_t reopened() const { return reopened_; }

    private:
      void startPrecursor(const Row& r)
      {
        auto& p = lib_.precursors();
        p.mz.push_back(toFixed(r.precursor_mz));
        p.irt.push_back(static_cast<float>(r.rt));
        p.im.push_back(r.im > 0.0 ? static_cast<float>(r.im) : std::nanf(""));
        p.charge.push_back(static_cast<std::uint8_t>(r.precursor_charge));
        p.decoy.push_back(r.decoy != 0 ? std::uint8_t{1} : std::uint8_t{0});
        p.modified_sequence.push_back(lib_.strings().intern(r.modified_sequence));
        p.protein_group.push_back(lib_.strings().intern(r.protein_group));
        p.transition_begin.push_back(static_cast<std::uint32_t>(lib_.transitionCount()));
        p.transition_count.push_back(0);
        current_id_owned_.assign(r.precursor_id);
        current_id_ = current_id_owned_;
        current_decoy_ = r.decoy != 0;
        have_current_ = true;
        count_ = 0;
      }

      void closeCurrent()
      {
        if (!have_current_) { return; }
        lib_.precursors().transition_count.back() = static_cast<std::uint32_t>(count_);
      }

      Library& lib_;
      std::unordered_set<std::size_t> seen_;
      std::size_t reopened_ = 0;
      std::string current_id_owned_;
      std::string_view current_id_;
      bool have_current_ = false;
      bool current_decoy_ = false;
      std::size_t count_ = 0;
    };
  } // namespace

  void DIANNLibraryFile::load(const std::string& filename, Library& library)
  {
    if (filename.ends_with(".parquet")) { loadParquet(filename, library); }
    else { loadTSV(filename, library); }
  }

  void DIANNLibraryFile::loadTSV(const std::string& filename, Library& library)
  {
    std::ifstream in(filename);
    if (!in) { throw std::runtime_error("cannot open library: " + filename); }

    std::string header;
    if (!std::getline(in, header)) { throw std::runtime_error("empty library: " + filename); }
    if (!header.empty() && header.back() == '\r') { header.pop_back(); }

    std::unordered_map<std::string, int> index;
    {
      int i = 0;
      std::size_t pos = 0;
      while (pos <= header.size())
      {
        const std::size_t tab = header.find('\t', pos);
        const std::size_t end = tab == std::string::npos ? header.size() : tab;
        index.emplace(header.substr(pos, end - pos), i++);
        if (tab == std::string::npos) { break; }
        pos = tab + 1;
      }
    }

    auto col = [&](const char* name) {
      const auto it = index.find(name);
      return it == index.end() ? -1 : it->second;
    };

    const int c_id = col(Columns::PRECURSOR_ID);
    const int c_seq = col(Columns::MODIFIED_SEQUENCE);
    const int c_pg = col(Columns::PROTEIN_GROUP);
    const int c_pmz = col(Columns::PRECURSOR_MZ);
    const int c_qmz = col(Columns::PRODUCT_MZ);
    const int c_rt = col(Columns::RT);
    const int c_im = col(Columns::IM);
    const int c_int = col(Columns::RELATIVE_INTENSITY);
    const int c_z = col(Columns::PRECURSOR_CHARGE);
    const int c_fz = col(Columns::FRAGMENT_CHARGE);
    const int c_ord = col(Columns::FRAGMENT_SERIES_NUMBER);
    const int c_dec = col(Columns::DECOY);
    const int c_ft = col(Columns::FRAGMENT_TYPE);
    const int c_lt = col(Columns::FRAGMENT_LOSS_TYPE);

    if (c_id < 0 || c_pmz < 0 || c_qmz < 0)
    {
      throw std::runtime_error("not a DIA-NN library (missing Precursor.Id / Precursor.Mz / "
                               "Product.Mz): " + filename);
    }

    Builder builder(library);
    std::string line;
    std::vector<std::string_view> f;
    while (std::getline(in, line))
    {
      if (line.empty()) { continue; }
      if (line.back() == '\r') { line.pop_back(); }

      f.clear();
      std::size_t pos = 0;
      while (pos <= line.size())
      {
        const std::size_t tab = line.find('\t', pos);
        const std::size_t end = tab == std::string::npos ? line.size() : tab;
        f.emplace_back(line.data() + pos, end - pos);
        if (tab == std::string::npos) { break; }
        pos = tab + 1;
      }

      auto get = [&](int i) -> std::string_view {
        return (i >= 0 && static_cast<std::size_t>(i) < f.size()) ? f[i] : std::string_view{};
      };

      Row r;
      r.precursor_id = get(c_id);
      r.modified_sequence = get(c_seq);
      r.protein_group = get(c_pg);
      r.precursor_mz = toDouble(get(c_pmz));
      r.product_mz = toDouble(get(c_qmz));
      r.rt = toDouble(get(c_rt));
      r.im = toDouble(get(c_im));
      r.intensity = toDouble(get(c_int));
      r.precursor_charge = toLong(get(c_z));
      r.fragment_charge = toLong(get(c_fz));
      r.ordinal = toLong(get(c_ord));
      r.decoy = toLong(get(c_dec));
      r.fragment_type = get(c_ft);
      r.loss_type = get(c_lt);
      builder.append(r);
    }
    builder.finish();
    if (builder.reopened())
    {
      throw std::runtime_error(
        "library rows are not grouped by precursor: " + std::to_string(builder.reopened())
        + " precursors reappear after another precursor's rows. Sort the file so each "
        "precursor's transitions are contiguous.");
    }
  }

  void DIANNLibraryFile::loadParquet(const std::string& filename, Library& library)
  {
    auto infile = arrow::io::ReadableFile::Open(filename);
    if (!infile.ok()) { throw std::runtime_error("cannot open library: " + filename); }

    auto reader_result = parquet::arrow::OpenFile(*infile, arrow::default_memory_pool());
    if (!reader_result.ok()) { throw std::runtime_error("not a Parquet file: " + filename); }
    std::unique_ptr<parquet::arrow::FileReader> reader = std::move(*reader_result);

    std::shared_ptr<arrow::Table> table;
    if (!reader->ReadTable(&table).ok())
    {
      throw std::runtime_error("cannot read Parquet table: " + filename);
    }

    // Combine chunks up front. A library at proteome scale exceeds the 2 GB
    // limit of a 32-bit-offset StringArray, so string columns arrive chunked;
    // resolving (global row) -> (chunk, index) per access is the alternative,
    // and is what the extractor will do for the far larger spectrum tables.
    // Here the fixture is small enough that one combine is simpler and safe.
    {
      auto combined = table->CombineChunks(arrow::default_memory_pool());
      if (!combined.ok()) { throw std::runtime_error("cannot combine Parquet chunks: " + filename); }
      table = *combined;
    }

    // Decode dictionary and large_string columns rather than rejecting them.
    // Parquet writes repeated strings dictionary-encoded by default -- which is
    // the very property that makes this path cheap -- and any library over 2 GB
    // of characters needs large_string. Casting straight to StringArray yields
    // nullptr for both, after which the reader blamed the file for "missing
    // Precursor.Id" when the column was present all along.
    std::vector<std::shared_ptr<arrow::Array>> keep_alive;
    auto str = [&](const char* n) -> std::shared_ptr<arrow::StringArray> {
      const int i = table->schema()->GetFieldIndex(n);
      if (i < 0) { return nullptr; }
      std::shared_ptr<arrow::Array> a = table->column(i)->chunk(0);
      if (a->type_id() != arrow::Type::STRING)
      {
        auto casted = arrow::compute::Cast(arrow::Datum(a), arrow::utf8());
        if (!casted.ok()) { return nullptr; }
        a = casted->make_array();
        keep_alive.push_back(a);
      }
      return std::dynamic_pointer_cast<arrow::StringArray>(a);
    };
    auto num = [&](const char* n) -> std::shared_ptr<arrow::Array> {
      const int i = table->schema()->GetFieldIndex(n);
      return i < 0 ? nullptr : table->column(i)->chunk(0);
    };
    auto at = [](const std::shared_ptr<arrow::Array>& a, int64_t i) -> double {
      if (!a || a->IsNull(i)) { return 0.0; }
      switch (a->type_id())
      {
        case arrow::Type::DOUBLE: return static_cast<const arrow::DoubleArray&>(*a).Value(i);
        case arrow::Type::FLOAT:  return static_cast<const arrow::FloatArray&>(*a).Value(i);
        case arrow::Type::INT64:  return static_cast<double>(static_cast<const arrow::Int64Array&>(*a).Value(i));
        case arrow::Type::INT32:  return static_cast<double>(static_cast<const arrow::Int32Array&>(*a).Value(i));
        case arrow::Type::INT16:  return static_cast<double>(static_cast<const arrow::Int16Array&>(*a).Value(i));
        case arrow::Type::INT8:   return static_cast<double>(static_cast<const arrow::Int8Array&>(*a).Value(i));
        case arrow::Type::UINT64: return static_cast<double>(static_cast<const arrow::UInt64Array&>(*a).Value(i));
        case arrow::Type::UINT32: return static_cast<double>(static_cast<const arrow::UInt32Array&>(*a).Value(i));
        case arrow::Type::UINT16: return static_cast<double>(static_cast<const arrow::UInt16Array&>(*a).Value(i));
        case arrow::Type::UINT8:  return static_cast<double>(static_cast<const arrow::UInt8Array&>(*a).Value(i));
        // Decoy is naturally a bool, and that is what pandas and polars write.
        // Falling through to 0.0 here silently made every decoy a target.
        case arrow::Type::BOOL:   return static_cast<const arrow::BooleanArray&>(*a).Value(i) ? 1.0 : 0.0;
        default: return std::nan("");
      }
    };
    auto sv = [](const std::shared_ptr<arrow::StringArray>& a, int64_t i) -> std::string_view {
      if (!a || a->IsNull(i)) { return {}; }
      return a->GetView(i);
    };

    const auto a_id = str(Columns::PRECURSOR_ID);
    const auto a_seq = str(Columns::MODIFIED_SEQUENCE);
    const auto a_pg = str(Columns::PROTEIN_GROUP);
    const auto a_ft = str(Columns::FRAGMENT_TYPE);
    const auto a_lt = str(Columns::FRAGMENT_LOSS_TYPE);
    const auto a_pmz = num(Columns::PRECURSOR_MZ);
    const auto a_qmz = num(Columns::PRODUCT_MZ);
    const auto a_rt = num(Columns::RT);
    const auto a_im = num(Columns::IM);
    const auto a_int = num(Columns::RELATIVE_INTENSITY);
    const auto a_z = num(Columns::PRECURSOR_CHARGE);
    const auto a_fz = num(Columns::FRAGMENT_CHARGE);
    const auto a_ord = num(Columns::FRAGMENT_SERIES_NUMBER);
    const auto a_dec = num(Columns::DECOY);

    if (!a_id || !a_pmz || !a_qmz)
    {
      throw std::runtime_error("not a DIA-NN library (missing Precursor.Id / Precursor.Mz / "
                               "Product.Mz): " + filename);
    }

    const int64_t rows = table->num_rows();
    library.reserve(static_cast<std::size_t>(rows) / 8, static_cast<std::size_t>(rows));
    library.strings().reserve(static_cast<std::size_t>(rows) / 8,
                              static_cast<std::size_t>(rows) * 4);

    Builder builder(library);
    for (int64_t i = 0; i < rows; ++i)
    {
      Row r;
      r.precursor_id = sv(a_id, i);
      r.modified_sequence = sv(a_seq, i);
      r.protein_group = sv(a_pg, i);
      r.fragment_type = sv(a_ft, i);
      r.loss_type = sv(a_lt, i);
      r.precursor_mz = at(a_pmz, i);
      r.product_mz = at(a_qmz, i);
      r.rt = at(a_rt, i);
      r.im = at(a_im, i);
      r.intensity = at(a_int, i);
      // at() yields NaN for a type it cannot decode or a null cell; casting
      // that to long is undefined behaviour, which is exactly what the toFixed
      // guard was added to eliminate two functions away.
      auto as_long = [](double v) -> long { return std::isnan(v) ? 0L : static_cast<long>(v); };
      r.precursor_charge = as_long(at(a_z, i));
      r.fragment_charge = as_long(at(a_fz, i));
      r.ordinal = as_long(at(a_ord, i));
      r.decoy = as_long(at(a_dec, i));
      builder.append(r);
    }
    builder.finish();
    if (builder.reopened())
    {
      throw std::runtime_error(
        "library rows are not grouped by precursor: " + std::to_string(builder.reopened())
        + " precursors reappear after another precursor's rows. Sort the file so each "
        "precursor's transitions are contiguous.");
    }
  }

  void DIANNLibraryFile::storeTSV(const std::string& filename, const Library& library)
  {
    std::ofstream out(filename);
    if (!out) { throw std::runtime_error("cannot write library: " + filename); }

    out << Columns::PRECURSOR_ID << '\t' << Columns::MODIFIED_SEQUENCE << '\t'
        << Columns::PRECURSOR_CHARGE << '\t' << Columns::DECOY << '\t'
        << Columns::RT << '\t' << Columns::IM << '\t'
        << Columns::PRECURSOR_MZ << '\t' << Columns::PRODUCT_MZ << '\t'
        << Columns::RELATIVE_INTENSITY << '\t' << Columns::FRAGMENT_TYPE << '\t'
        << Columns::FRAGMENT_CHARGE << '\t' << Columns::FRAGMENT_SERIES_NUMBER << '\t'
        << Columns::FRAGMENT_LOSS_TYPE << '\t' << Columns::PROTEIN_GROUP << '\n';

    out << std::defaultfloat;
    const auto& p = library.precursors();
    const auto& t = library.transitions();
    for (std::size_t i = 0; i < library.precursorCount(); ++i)
    {
      const auto seq = library.strings().get(p.modified_sequence[i]);
      const auto pg = library.strings().get(p.protein_group[i]);
      const int z = p.charge[i];
      const bool has_im = !std::isnan(p.im[i]);
      const std::uint32_t begin = p.transition_begin[i];
      for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
      {
        const std::uint32_t j = begin + k;
        // DIA-NN gives a decoy the same Modified.Sequence, charge and
        // Precursor.Mz as its target, so <sequence><charge> alone is not unique.
        // Without the suffix, reloading our own output merged each such pair and
        // re-labelled the decoy's transitions as target evidence -- 33,390
        // precursors became 33,386, silently, with no transitions lost.
        out << seq << z << (p.decoy[i] ? "_decoy" : "") << '\t'
            << seq << '\t' << z << '\t'
            << static_cast<int>(p.decoy[i]) << '\t';
        out.precision(9);
        out << p.irt[i] << '\t';
        if (has_im) { out << p.im[i]; } else { out << 0; }
        out << '\t';
        out.precision(10);
        out << fromFixed(p.mz[i]) << '\t' << fromFixed(t.product_mz[j]) << '\t';
        out.precision(9);
        out << t.library_intensity[j] << '\t' << toString(t.type[j]) << '\t'
            << static_cast<int>(t.charge[j]) << '\t'
            << static_cast<int>(t.ordinal[j]) << '\t'
            << toString(t.loss[j]) << '\t' << pg << '\n';
      }
    }

    // A stream checked only at open reports success on a full disk, an exceeded
    // quota or a broken mount, leaving a truncated library behind.
    out.flush();
    if (!out) { throw std::runtime_error("failed while writing library: " + filename); }
  }

} // namespace ODIA
