// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/DIANNLibraryFile.h>
#include <odia/TextWriter.h>

#include <arrow/api.h>
#include <arrow/compute/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <arrow/builder.h>
#include <arrow/table.h>

#include <charconv>
#include <cstring>
#include <cmath>
#include <fstream>
#include <iostream>
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

    /// Residues in a modified sequence, ignoring bracketed modification names.
    ///
    /// The name must not be scanned for residues: "UniMod" contains an M.
    std::size_t residueCount(std::string_view seq)
    {
      std::size_t n = 0;
      for (std::size_t i = 0; i < seq.size(); ++i)
      {
        if (seq[i] == '(' || seq[i] == '[')
        {
          const char close = seq[i] == '(' ? ')' : ']';
          while (i < seq.size() && seq[i] != close) { ++i; }
          continue;
        }
        if (seq[i] >= 'A' && seq[i] <= 'Z') { ++n; }
      }
      return n;
    }

    /// One parsed row, in the order the builder needs it.
    struct Row
    {
      std::string_view precursor_id;
      std::string_view modified_sequence;
      std::string_view protein_group;
      double precursor_mz = 0, product_mz = 0;
      double rt = 0, im = 0, ccs = 0, intensity = 0;
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
        // Reject fragments that cannot exist, rather than storing them.
        //
        // An ordinal past the stored width wraps into a *plausible* ordinal --
        // 300 becomes 44 -- so the row would survive every downstream check as
        // a different, wrong fragment. Ordinal 0, an ordinal at or beyond the
        // peptide length, and a non-positive charge are not fragments at all.
        if (r.ordinal <= 0 || r.ordinal > 255 ||
            r.fragment_charge <= 0 || r.fragment_charge > 127 ||
            static_cast<std::size_t>(r.ordinal) >= residueCount(r.modified_sequence))
        {
          ++out_of_range_;
          return;
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

      /// Transitions dropped because an ordinal or charge could not be stored.
      std::size_t outOfRange() const { return out_of_range_; }

    private:
      void startPrecursor(const Row& r)
      {
        auto& p = lib_.precursors();
        p.mz.push_back(toFixed(r.precursor_mz));
        p.irt.push_back(static_cast<float>(r.rt));
        p.im.push_back(r.im > 0.0 ? static_cast<float>(r.im) : std::nanf(""));
        // CCS was written by both writers and read by NEITHER reader, so a
        // library round-tripped through ODIA lost it silently: the generator
        // predicts it, storeTSV writes it, and every search that reloaded that
        // file saw an empty column. Measured on a round trip -- CCS 627.119385
        // went in and an empty field came out.
        //
        // It matters because CCS is the ONLY mobility information a generated
        // library carries: IM is left unset (NaN) because 1/K0 needs the drift
        // gas and the instrument calibration to derive. Dropping CCS therefore
        // discards the mobility dimension entirely on reload.
        p.ccs.push_back(r.ccs > 0.0 ? static_cast<float>(r.ccs) : std::nanf(""));
        p.charge.push_back(static_cast<std::uint8_t>(r.precursor_charge));
        p.decoy.push_back(r.decoy != 0 ? std::uint8_t{1} : std::uint8_t{0});
        p.modified_sequence.push_back(lib_.strings().intern(r.modified_sequence));
        p.protein_group.push_back(lib_.strings().intern(r.protein_group));
        p.transition_begin.push_back(static_cast<std::uint32_t>(lib_.transitionCount()));
        p.transition_count.push_back(0);
        // Appending breaks any m/z ordering; leaving the flag set would make
        // lowerBound binary-search an unsorted array and return an index
        // indistinguishable from a real answer.
        lib_.markUnsorted();
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
      std::size_t out_of_range_ = 0;
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
    // The interning index answered "have I seen this string" while reading and
    // is never consulted again -- handles resolve through the arena's entries.
    // Held for the library's lifetime it is one hash node per distinct string
    // plus a bucket array, for a question nobody asks twice.
    library.strings().releaseLookup();
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
    const int c_ccs = col(Columns::CCS);
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
      r.ccs = toDouble(get(c_ccs));
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
    if (builder.outOfRange())
    {
      std::cerr << "warning: dropped " << builder.outOfRange()
                << " transitions whose fragment ordinal or charge cannot describe "
                   "a fragment of the stated peptide\n";
    }
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
    const auto a_ccs = num(Columns::CCS);
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
      r.ccs = at(a_ccs, i);
      r.intensity = at(a_int, i);
      // at() yields NaN for a type it cannot decode or a null cell; casting
      // that to long is undefined behaviour, which is exactly what the toFixed
      // guard was added to eliminate two functions away.
      // NaN was guarded; +-inf and any |v| >= 2^63 still reached the cast and
      // UBSan flagged them. On x86-64 they yield INT64_MIN, which narrows to
      // charge 0 -- a wrong value rather than a crash.
      auto as_long = [](double v) -> long {
        if (!std::isfinite(v)) { return 0L; }
        if (v <= static_cast<double>(std::numeric_limits<long>::min())) { return 0L; }
        if (v >= static_cast<double>(std::numeric_limits<long>::max())) { return 0L; }
        return static_cast<long>(v);
      };
      r.precursor_charge = as_long(at(a_z, i));
      r.fragment_charge = as_long(at(a_fz, i));
      r.ordinal = as_long(at(a_ord, i));
      r.decoy = as_long(at(a_dec, i));
      builder.append(r);
    }
    builder.finish();
    if (builder.outOfRange())
    {
      std::cerr << "warning: dropped " << builder.outOfRange()
                << " transitions whose fragment ordinal or charge cannot describe "
                   "a fragment of the stated peptide\n";
    }
    if (builder.reopened())
    {
      throw std::runtime_error(
        "library rows are not grouped by precursor: " + std::to_string(builder.reopened())
        + " precursors reappear after another precursor's rows. Sort the file so each "
        "precursor's transitions are contiguous.");
    }
  }

  /// 7.28 GB for the human library, and it was the second-largest item in
  /// Phase 1 -- 211.7 s, 26% of the stage, and untimed until now. At 34.4 MB/s
  /// it sat on `ofstream <<`'s measured 25.4 MB/s line and nowhere near the
  /// NVMe's 1,042 MB/s, so the cost is number formatting, not disk. TextWriter
  /// renders the same digits through to_chars into a reusable buffer; see its
  /// header for why the precisions below are passed explicitly, and
  /// `odia_tsv_writers` for the byte-for-byte check against this function's
  /// previous ofstream form.
  void DIANNLibraryFile::store(const std::string& filename, const Library& library)
  {
    if (filename.ends_with(".parquet")) { storeParquet(filename, library); }
    else { storeTSV(filename, library); }
  }

  void DIANNLibraryFile::storeParquet(const std::string& filename, const Library& library)
  {
    const auto& p = library.precursors();
    const auto& t = library.transitions();

    // One row per TRANSITION, matching the TSV exactly, because loadParquet and
    // loadTSV populate the same structures from the same column names. A
    // narrower per-precursor schema would be smaller still but would not round
    // trip through the existing reader.
    arrow::StringBuilder b_id, b_seq, b_pg, b_ft, b_lt;
    arrow::Int32Builder b_z, b_dec, b_fz, b_ord;
    arrow::DoubleBuilder b_rt, b_im, b_ccs, b_pmz, b_qmz, b_int;

    auto ok = [](const arrow::Status& st) {
      if (!st.ok()) { throw std::runtime_error("parquet build: " + st.ToString()); }
    };

    for (std::size_t i = 0; i < library.precursorCount(); ++i)
    {
      const auto seq = library.strings().get(p.modified_sequence[i]);
      const auto pg = library.strings().get(p.protein_group[i]);
      const std::string seq_s(seq.data(), seq.size());
      const std::string pg_s(pg.data(), pg.size());
      const int z = static_cast<int>(p.charge[i]);
      // The suffix is what makes a decoy's id distinct from its target's: they
      // share sequence, charge and precursor m/z, so <sequence><charge> alone
      // collides and a reload silently merges the pair. See storeTSV.
      const std::string id = seq_s + std::to_string(z) + (p.decoy[i] ? "_decoy" : "");
      // NaN means ABSENT and is written as null, not as 0. A zero ion mobility
      // is a VALUE that gates every peak against 0 and rejects them all; the
      // TSV writer emits 0 here and only survives because the reader maps 0
      // back to NaN. Parquet has real nulls, so this path does not need that
      // round-trip coincidence.
      const bool has_im = i < p.im.size() && !std::isnan(p.im[i]);
      const bool has_ccs = i < p.ccs.size() && !std::isnan(p.ccs[i]);
      const bool has_rt = !std::isnan(p.irt[i]);

      const std::uint32_t begin = p.transition_begin[i];
      for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
      {
        const std::uint32_t j = begin + k;
        ok(b_id.Append(id));
        ok(b_seq.Append(seq_s));
        ok(b_pg.Append(pg_s));
        ok(b_z.Append(z));
        ok(b_dec.Append(static_cast<int>(p.decoy[i])));
        if (has_rt) { ok(b_rt.Append(double(p.irt[i]))); } else { ok(b_rt.AppendNull()); }
        if (has_im) { ok(b_im.Append(double(p.im[i]))); } else { ok(b_im.AppendNull()); }
        if (has_ccs) { ok(b_ccs.Append(double(p.ccs[i]))); } else { ok(b_ccs.AppendNull()); }
        ok(b_pmz.Append(fromFixed(p.mz[i])));
        ok(b_qmz.Append(fromFixed(t.product_mz[j])));
        ok(b_int.Append(double(t.library_intensity[j])));
        ok(b_ft.Append(std::string(toString(t.type[j]))));
        ok(b_fz.Append(static_cast<int>(t.charge[j])));
        ok(b_ord.Append(static_cast<int>(t.ordinal[j])));
        ok(b_lt.Append(std::string("noloss")));
      }
    }

    std::shared_ptr<arrow::Array> a_id, a_seq, a_pg, a_ft, a_lt, a_z, a_dec,
                                  a_fz, a_ord, a_rt, a_im, a_ccs, a_pmz, a_qmz, a_int;
    ok(b_id.Finish(&a_id));   ok(b_seq.Finish(&a_seq)); ok(b_pg.Finish(&a_pg));
    ok(b_ft.Finish(&a_ft));   ok(b_lt.Finish(&a_lt));   ok(b_z.Finish(&a_z));
    ok(b_dec.Finish(&a_dec)); ok(b_fz.Finish(&a_fz));   ok(b_ord.Finish(&a_ord));
    ok(b_rt.Finish(&a_rt));   ok(b_im.Finish(&a_im));   ok(b_ccs.Finish(&a_ccs));
    ok(b_pmz.Finish(&a_pmz)); ok(b_qmz.Finish(&a_qmz)); ok(b_int.Finish(&a_int));

    auto schema = arrow::schema({
      arrow::field(Columns::PRECURSOR_ID, arrow::utf8()),
      arrow::field(Columns::MODIFIED_SEQUENCE, arrow::utf8()),
      arrow::field(Columns::PRECURSOR_CHARGE, arrow::int32()),
      arrow::field(Columns::DECOY, arrow::int32()),
      arrow::field(Columns::RT, arrow::float64()),
      arrow::field(Columns::IM, arrow::float64()),
      arrow::field(Columns::CCS, arrow::float64()),
      arrow::field(Columns::PRECURSOR_MZ, arrow::float64()),
      arrow::field(Columns::PRODUCT_MZ, arrow::float64()),
      arrow::field(Columns::RELATIVE_INTENSITY, arrow::float64()),
      arrow::field(Columns::FRAGMENT_TYPE, arrow::utf8()),
      arrow::field(Columns::FRAGMENT_CHARGE, arrow::int32()),
      arrow::field(Columns::FRAGMENT_SERIES_NUMBER, arrow::int32()),
      arrow::field(Columns::FRAGMENT_LOSS_TYPE, arrow::utf8()),
      arrow::field(Columns::PROTEIN_GROUP, arrow::utf8()),
    });
    auto table = arrow::Table::Make(schema,
      {a_id, a_seq, a_z, a_dec, a_rt, a_im, a_ccs, a_pmz, a_qmz, a_int,
       a_ft, a_fz, a_ord, a_lt, a_pg});

    auto outfile = arrow::io::FileOutputStream::Open(filename);
    if (!outfile.ok()) { throw std::runtime_error("cannot write library: " + filename); }
    // Dictionary encoding is the whole point: the sequence and protein group
    // repeat once per transition, twelve times per precursor.
    auto props = parquet::WriterProperties::Builder()
                   .compression(parquet::Compression::ZSTD)
                   ->enable_dictionary()
                   ->build();
    const auto st = parquet::arrow::WriteTable(*table, arrow::default_memory_pool(),
                                               *outfile, 1 << 20, props);
    if (!st.ok()) { throw std::runtime_error("cannot write Parquet: " + st.ToString()); }
  }

  void DIANNLibraryFile::storeTSV(const std::string& filename, const Library& library)
  {
    TextWriter out(filename);

    for (const char* column : {Columns::PRECURSOR_ID, Columns::MODIFIED_SEQUENCE,
                               Columns::PRECURSOR_CHARGE, Columns::DECOY,
                               Columns::RT, Columns::IM,
                               Columns::PRECURSOR_MZ, Columns::PRODUCT_MZ,
                               Columns::RELATIVE_INTENSITY, Columns::FRAGMENT_TYPE,
                               Columns::FRAGMENT_CHARGE, Columns::FRAGMENT_SERIES_NUMBER,
                               Columns::FRAGMENT_LOSS_TYPE, Columns::PROTEIN_GROUP})
    {
      out.put(column);
      out.put('\t');
    }
    out.put(Columns::CCS);
    out.put('\n');

    const auto& p = library.precursors();
    const auto& t = library.transitions();
    for (std::size_t i = 0; i < library.precursorCount(); ++i)
    {
      const auto seq = library.strings().get(p.modified_sequence[i]);
      const auto pg = library.strings().get(p.protein_group[i]);
      const int z = p.charge[i];
      const bool has_im = !std::isnan(p.im[i]);
      const bool has_ccs = i < p.ccs.size() && !std::isnan(p.ccs[i]);
      const std::uint32_t begin = p.transition_begin[i];
      for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
      {
        const std::uint32_t j = begin + k;
        // DIA-NN gives a decoy the same Modified.Sequence, charge and
        // Precursor.Mz as its target, so <sequence><charge> alone is not unique.
        // Without the suffix, reloading our own output merged each such pair and
        // re-labelled the decoy's transitions as target evidence -- 33,390
        // precursors became 33,386, silently, with no transitions lost.
        out.put(seq); out.integer(z); out.put(p.decoy[i] ? "_decoy" : ""); out.put('\t');
        out.put(seq); out.put('\t'); out.integer(z); out.put('\t');
        out.integer(p.decoy[i]); out.put('\t');
        // An unpredicted retention time is written as an empty field, not as
        // "nan": no TSV consumer accepts the latter, and it propagates into
        // anything that reads the library back.
        if (!std::isnan(p.irt[i])) { out.number(p.irt[i], 9); }
        out.put('\t');
        if (has_im) { out.number(p.im[i], 9); } else { out.integer(0); }
        out.put('\t');
        out.number(fromFixed(p.mz[i]), 10); out.put('\t');
        out.number(fromFixed(t.product_mz[j]), 10); out.put('\t');
        out.number(t.library_intensity[j], 9); out.put('\t');
        out.put(toString(t.type[j])); out.put('\t');
        out.integer(t.charge[j]); out.put('\t');
        out.integer(t.ordinal[j]); out.put('\t');
        out.put(toString(t.loss[j])); out.put('\t');
        out.put(pg); out.put('\t');
        // Empty rather than "nan" when absent, as for RT: no TSV consumer
        // accepts the latter, and a 0 here would read as a real cross-section.
        if (has_ccs) { out.number(p.ccs[i], 9); }
        out.put('\n');
      }
    }

    // A file checked only at open reports success on a full disk, an exceeded
    // quota or a broken mount, leaving a truncated library behind.
    out.close();
  }

} // namespace ODIA
