// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/OSWPQLibraryFile.h>

#include <odia/ArrowColumn.h>
#include <odia/ZipArchive.h>

#include <arrow/api.h>
#include <parquet/arrow/reader.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace ODIA
{

  namespace
  {
    constexpr const char* PRECURSORS = "library/precursors.parquet";
    constexpr const char* TRANSITIONS = "library/transitions.parquet";
    constexpr const char* METADATA = "library/metadata.json";

    /// The schema version this reader was written against. Refusing a newer one
    /// is the point of the field: the format has already taken one breaking
    /// change on disk (library_intensity's width), and guessing at an unknown
    /// version would produce a library rather than an error.
    constexpr int SUPPORTED_SCHEMA_VERSION = 1;

    [[noreturn]] void fail(const std::string& what)
    {
      throw std::runtime_error("OSWPQLibraryFile: " + what);
    }

    /// Read one Parquet entry, materialising only the columns asked for.
    ///
    /// Projection is not a micro-optimisation here. traml_id is denormalised
    /// into the transition table -- 78.6 M copies of ~7.1 M distinct strings,
    /// the single largest contributor to bundle size -- and ODIA never needs
    /// it, because it addresses precursors by index. Reading it would dominate
    /// both the IO and the peak memory of a load that then discards it.
    std::shared_ptr<arrow::Table> readTable(const ZipArchive& zip, const std::string& entry,
                                            const std::vector<std::string>& wanted)
    {
      auto file = zip.open(entry);
      auto reader_result = parquet::arrow::OpenFile(file, arrow::default_memory_pool());
      if (!reader_result.ok())
      {
        fail("cannot read " + entry + ": " + reader_result.status().ToString());
      }
      std::unique_ptr<parquet::arrow::FileReader> reader = std::move(*reader_result);

      std::shared_ptr<arrow::Schema> schema;
      if (!reader->GetSchema(&schema).ok()) { fail("cannot read the schema of " + entry); }

      // Column indices are Arrow field indices only while the schema is flat.
      // Every column of this format is a primitive, but a nested one would make
      // the mapping wrong rather than absent, so fall back to reading whole.
      bool flat = true;
      for (const auto& f : schema->fields())
      {
        if (arrow::is_nested(f->type()->id())) { flat = false; break; }
      }

      std::vector<int> indices;
      if (flat)
      {
        for (const auto& name : wanted)
        {
          const int i = schema->GetFieldIndex(name);
          if (i >= 0) { indices.push_back(i); }
        }
      }

      std::shared_ptr<arrow::Table> table;
      const auto status = indices.empty() ? reader->ReadTable(&table)
                                          : reader->ReadTable(indices, &table);
      if (!status.ok()) { fail("cannot read " + entry + ": " + status.ToString()); }
      // Deliberately no CombineChunks: it is what makes a >2 GB text column look
      // readable as chunk(0) and then read short. ChunkedColumn handles chunks.
      return table;
    }

    /// Pull one field out of the metadata blob.
    ///
    /// A JSON dependency is not warranted for four scalars, and OpenMS writes
    /// this block itself with a fixed shape. The extraction is deliberately
    /// narrow: it finds `"key"` and takes the token after the colon.
    std::string jsonField(std::string_view json, std::string_view key)
    {
      const std::string needle = "\"" + std::string(key) + "\"";
      const auto k = json.find(needle);
      if (k == std::string_view::npos) { return {}; }
      auto p = json.find(':', k + needle.size());
      if (p == std::string_view::npos) { return {}; }
      ++p;
      while (p < json.size() && (json[p] == ' ' || json[p] == '\n' || json[p] == '\t' ||
                                 json[p] == '\r')) { ++p; }
      if (p >= json.size()) { return {}; }
      if (json[p] == '"')
      {
        const auto end = json.find('"', p + 1);
        if (end == std::string_view::npos) { return {}; }
        return std::string(json.substr(p + 1, end - p - 1));
      }
      const auto end = json.find_first_of(",}\n", p);
      return std::string(json.substr(p, (end == std::string_view::npos ? json.size() : end) - p));
    }

    /// The census for one table, e.g. counts.precursors.total.
    ///
    /// Scoped to the object that follows the table's name, so "total" is read
    /// from the right one of the five count blocks.
    std::size_t censusTotal(std::string_view json, std::string_view table)
    {
      const std::string needle = "\"" + std::string(table) + "\"";
      const auto k = json.find(needle);
      if (k == std::string_view::npos) { return 0; }
      const auto value = jsonField(json.substr(k), "total");
      if (value.empty()) { return 0; }
      try { return static_cast<std::size_t>(std::stoull(value)); }
      catch (const std::exception&) { return 0; }
    }

    /// The neutral loss, recovered from the annotation.
    ///
    /// The transition schema has no loss column: type, ordinal and charge are
    /// separate fields but the loss survives only inside the annotation string,
    /// as the suffix of e.g. "y7-H2O". Dropping it would silently turn every
    /// loss ion into its parent, which is a different mass.
    LossType lossFromAnnotation(std::string_view annotation)
    {
      const auto minus = annotation.find('-');
      if (minus == std::string_view::npos) { return LossType::None; }
      auto rest = annotation.substr(minus + 1);
      // A charge suffix ("y7-H2O^2") is not part of the loss.
      const auto caret = rest.find('^');
      if (caret != std::string_view::npos) { rest = rest.substr(0, caret); }
      return parseLossType(rest);
    }
  } // namespace

  bool OSWPQLibraryFile::isLibraryBundle(const std::string& filename)
  {
    try
    {
      ZipArchive zip(filename);
      return zip.has(PRECURSORS) && zip.has(TRANSITIONS);
    }
    catch (const std::exception&)
    {
      return false;
    }
  }

  void OSWPQLibraryFile::load(const std::string& filename, Library& library, Stats* stats)
  {
    Stats local;
    Stats& s = stats == nullptr ? local : *stats;
    s = Stats{};

    ZipArchive zip(filename);
    if (!zip.has(PRECURSORS) || !zip.has(TRANSITIONS))
    {
      std::string have;
      for (const auto& e : zip.entries()) { have += "\n  " + e; }
      fail(filename + " is a ZIP archive but carries no library bundle; its entries are:" + have);
    }

    // The version gate comes before any parsing, which is the whole reason the
    // field exists.
    std::string metadata;
    if (zip.has(METADATA))
    {
      metadata = zip.read(METADATA);
      const auto version = jsonField(metadata, "schema_version");
      if (!version.empty())
      {
        try { s.schema_version = std::stoi(version); }
        catch (const std::exception&) { s.schema_version = 0; }
      }
      if (s.schema_version > SUPPORTED_SCHEMA_VERSION)
      {
        fail(filename + " declares schema_version " + std::to_string(s.schema_version) +
             ", and this reader understands " + std::to_string(SUPPORTED_SCHEMA_VERSION) +
             ". Reading it anyway would produce a library rather than an error.");
      }
      s.generator = jsonField(metadata, "generator");
      s.openms_version = jsonField(metadata, "openms_version");
      s.census_precursors = censusTotal(metadata, "precursors");
      s.census_transitions = censusTotal(metadata, "transitions");
      s.census_present = s.census_precursors != 0 || s.census_transitions != 0;
    }

    library = Library{};

    // ---------------------------------------------------------------- precursors
    auto ptable = readTable(zip, PRECURSORS,
                            {"precursor_id", "precursor_mz", "charge", "library_rt",
                             "library_drift_time", "decoy", "modified_sequence",
                             "protein_accessions"});
    const std::int64_t nprec = ptable->num_rows();
    s.precursor_rows = static_cast<std::size_t>(nprec);

    const auto c_pid   = requiredColumn(ptable, "precursor_id", PRECURSORS);
    const auto c_pmz   = requiredColumn(ptable, "precursor_mz", PRECURSORS);
    const auto c_pchg  = requiredColumn(ptable, "charge", PRECURSORS);
    const auto c_prt   = column(ptable, "library_rt");
    const auto c_pim   = column(ptable, "library_drift_time");
    const auto c_pdec  = column(ptable, "decoy");
    const auto c_pseq  = requiredColumn(ptable, "modified_sequence", PRECURSORS);
    const auto c_pprot = column(ptable, "protein_accessions");

    auto& p = library.precursors();
    auto& arena = library.strings();
    library.reserve(static_cast<std::size_t>(nprec), 0);

    std::unordered_map<std::int64_t, std::uint32_t> id_to_index;
    id_to_index.reserve(static_cast<std::size_t>(nprec) * 2);

    for (std::int64_t r = 0; r < nprec; ++r)
    {
      const auto id = c_pid.getInt64(r, -1);
      // A duplicate id makes the join ambiguous, so the first row wins. Silence
      // here would be indistinguishable from a precursor that genuinely has no
      // fragments, so it is counted.
      if (!id_to_index.emplace(id, static_cast<std::uint32_t>(p.mz.size())).second)
      {
        ++s.duplicate_precursor_ids;
      }

      p.mz.push_back(toFixed(c_pmz.getDouble(r, 0.0)));

      const auto charge = c_pchg.getInt64(r, 0);
      if (charge < 0 || charge > 255) { ++s.unusable_precursor_charges; p.charge.push_back(0); }
      else { p.charge.push_back(static_cast<std::uint8_t>(charge)); }

      p.irt.push_back(static_cast<float>(c_prt.getDouble(r, std::numeric_limits<double>::quiet_NaN())));

      // -1 is the writer's "absent", and it is not a drift time. Storing it as
      // one would put every precursor of a bundle without ion mobility at a
      // negative 1/K0 that no frame can match.
      const double im = c_pim.getDouble(r, std::numeric_limits<double>::quiet_NaN());
      p.im.push_back(im < 0.0 ? std::numeric_limits<float>::quiet_NaN()
                              : static_cast<float>(im));

      p.decoy.push_back(c_pdec.getBool(r, false) ? 1 : 0);
      p.modified_sequence.push_back(arena.intern(c_pseq.getString(r)));
      p.protein_group.push_back(arena.intern(c_pprot.getString(r)));
    }
    p.transition_begin.assign(static_cast<std::size_t>(nprec), 0);
    p.transition_count.assign(static_cast<std::size_t>(nprec), 0);

    // --------------------------------------------------------------- transitions
    auto ttable = readTable(zip, TRANSITIONS,
                            {"precursor_id", "product_mz", "charge", "type", "annotation",
                             "ordinal", "library_intensity", "decoy"});
    const std::int64_t ntrans = ttable->num_rows();
    s.transition_rows = static_cast<std::size_t>(ntrans);

    const auto c_tpid  = requiredColumn(ttable, "precursor_id", TRANSITIONS);
    const auto c_tmz   = requiredColumn(ttable, "product_mz", TRANSITIONS);
    const auto c_tchg  = column(ttable, "charge");
    const auto c_ttype = column(ttable, "type");
    const auto c_tann  = column(ttable, "annotation");
    const auto c_tord  = column(ttable, "ordinal");
    // Upstream writes float64 here and the OpenDIAlyzer patch writes float32.
    // getDouble accepts either, so a bundle from either build reads; pinning
    // the width would reject half the files this reader exists to open.
    const auto c_tint  = column(ttable, "library_intensity");
    const auto c_tdec  = column(ttable, "decoy");

    // Owner index per transition, in file order. Resolving the join once costs
    // 4 bytes a row and saves a second pass over the hash map, which is the
    // expensive part -- ~78.6 M lookups on the benchmark library.
    std::vector<std::uint32_t> owner(static_cast<std::size_t>(ntrans),
                                     std::numeric_limits<std::uint32_t>::max());
    constexpr std::uint32_t NO_OWNER = std::numeric_limits<std::uint32_t>::max();

    for (std::int64_t r = 0; r < ntrans; ++r)
    {
      const auto it = id_to_index.find(c_tpid.getInt64(r, -1));
      if (it == id_to_index.end()) { ++s.orphan_transitions; continue; }
      owner[static_cast<std::size_t>(r)] = it->second;
      ++p.transition_count[it->second];
    }

    std::uint64_t running = 0;
    for (std::size_t i = 0; i < p.transition_count.size(); ++i)
    {
      if (p.transition_count[i] == 0) { ++s.childless_precursors; }
      if (running > std::numeric_limits<std::uint32_t>::max())
      {
        fail("more than 2^32 transitions; the CSR index is 32-bit");
      }
      p.transition_begin[i] = static_cast<std::uint32_t>(running);
      running += p.transition_count[i];
    }

    const std::size_t kept = static_cast<std::size_t>(running);
    auto& t = library.transitions();
    t.product_mz.assign(kept, MZ_INVALID);
    t.library_intensity.assign(kept, 0.0f);
    t.type.assign(kept, FragmentType::Unknown);
    t.ordinal.assign(kept, 0);
    t.charge.assign(kept, 0);
    t.loss.assign(kept, LossType::None);

    // Transitions are grouped by precursor in the file as written, but nothing
    // in the format requires it, and a library whose fragments belong to the
    // wrong precursor is worse than one that fails to load. Scattering by a
    // per-precursor cursor is correct whatever the row order.
    std::vector<std::uint32_t> cursor = p.transition_begin;

    for (std::int64_t r = 0; r < ntrans; ++r)
    {
      const std::uint32_t o = owner[static_cast<std::size_t>(r)];
      if (o == NO_OWNER) { continue; }
      const std::size_t at = cursor[o]++;

      t.product_mz[at] = toFixed(c_tmz.getDouble(r, 0.0));
      t.library_intensity[at] = static_cast<float>(c_tint.getDouble(r, 0.0));

      const auto charge = c_tchg.getInt64(r, 0);
      if (charge < std::numeric_limits<std::int8_t>::min() ||
          charge > std::numeric_limits<std::int8_t>::max())
      {
        ++s.unusable_transition_charges;
      }
      else { t.charge[at] = static_cast<std::int8_t>(charge); }

      // An OpenSWATH transition list need not be annotated at all: the upstream
      // sample bundle writes type "" and ordinal -1 for every one of its rows.
      // Dropping those, as the DIA-NN reader does for impossible ordinals,
      // would empty a legitimate library.
      const auto ordinal = c_tord.getInt64(r, 0);
      if (ordinal < 0 || ordinal > 255) { ++s.unusable_ordinals; }
      else { t.ordinal[at] = static_cast<std::uint8_t>(ordinal); }

      t.type[at] = parseFragmentType(c_ttype.getString(r));
      t.loss[at] = lossFromAnnotation(c_tann.getString(r));

      if (c_tdec.valid() && (c_tdec.getBool(r, false) ? 1 : 0) != p.decoy[o])
      {
        ++s.decoy_mismatches;
      }
    }

    s.census_agrees = s.census_present &&
                      s.census_precursors == s.precursor_rows &&
                      s.census_transitions == s.transition_rows;

    library.markUnsorted();
    library.shrinkToFit();
  }

} // namespace ODIA
