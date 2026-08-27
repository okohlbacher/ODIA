#include <odia/VendorDiaWindows.h>

#include <sqlite3.h>
#include <zip.h>

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

namespace ODIA
{
  namespace
  {
    /// The method directory is named after the acquisition method, so the entry
    /// cannot be named outright. Bruker also leaves dated backups beside it --
    /// `backup-2024-12-04.m/`, `backup-2024-12-04_1.m/` -- which hold the
    /// method as it was BEFORE the last edit. Taking the first match would
    /// therefore sometimes read a superseded scheme, so backups are skipped and
    /// only the live method is accepted.
    bool isLiveMethodEntry(const std::string& name)
    {
      if (name.find("diaSettings.diasqlite") == std::string::npos) { return false; }
      return name.find("backup-") == std::string::npos;
    }

    /// SQLite reads from a file, and the table sits inside an archive. Rather
    /// than link a VFS shim, spill the (24 kB) blob to a temp file and delete
    /// it again. Small enough that the copy is not worth avoiding.
    struct TempFile
    {
      std::filesystem::path path;
      ~TempFile()
      {
        std::error_code ec;
        if (!path.empty()) { std::filesystem::remove(path, ec); }
      }
    };
  }

  std::vector<VendorDiaWindow> readVendorDiaWindows(const std::string& mzpeak_path)
  {
    std::vector<VendorDiaWindow> out;

    int err = 0;
    zip_t* z = zip_open(mzpeak_path.c_str(), ZIP_RDONLY, &err);
    if (!z) { return out; }

    zip_int64_t found = -1;
    const zip_int64_t n = zip_get_num_entries(z, 0);
    for (zip_int64_t i = 0; i < n; ++i)
    {
      const char* name = zip_get_name(z, i, 0);
      if (name && isLiveMethodEntry(name)) { found = i; break; }
    }
    if (found < 0) { zip_close(z); return out; }

    zip_stat_t st;
    zip_stat_init(&st);
    if (zip_stat_index(z, found, 0, &st) != 0 || !(st.valid & ZIP_STAT_SIZE))
    {
      zip_close(z);
      return out;
    }

    std::vector<char> blob(static_cast<std::size_t>(st.size));
    zip_file_t* f = zip_fopen_index(z, found, 0);
    if (!f)
    {
      zip_close(z);
      return out;
    }
    const zip_int64_t got = zip_fread(f, blob.data(), blob.size());
    zip_fclose(f);
    zip_close(z);
    if (got < 0 || static_cast<std::size_t>(got) != blob.size()) { return out; }

    TempFile tmp;
    {
      std::error_code ec;
      auto dir = std::filesystem::temp_directory_path(ec);
      if (ec) { return out; }
      // Distinct per process AND per call, because two runs sharing a node
      // would otherwise race on one name and read each other's bytes.
      static unsigned counter = 0;
      tmp.path = dir / ("odia_diasettings_" + std::to_string(::getpid()) + "_" +
                        std::to_string(counter++) + ".sqlite");
      std::ofstream o(tmp.path, std::ios::binary);
      if (!o) { return out; }
      o.write(blob.data(), static_cast<std::streamsize>(blob.size()));
      if (!o) { return out; }
    }

    sqlite3* db = nullptr;
    if (sqlite3_open_v2(tmp.path.c_str(), &db, SQLITE_OPEN_READONLY, nullptr) != SQLITE_OK)
    {
      if (db) { sqlite3_close(db); }
      return out;
    }

    // Type 1 is a dia window; the table also carries a Type 0 header row whose
    // every column is null, and taking it would insert a window at 0 m/z with a
    // zero-width mobility band -- which the caller would then match against
    // nothing and silently keep deriving.
    const char* sql =
      "select OneOverK0Start, OneOverK0End, IsolationMz, IsolationWidth "
      "from DiaWindowsSpecification "
      "where OneOverK0Start is not null and OneOverK0End is not null "
      "  and IsolationMz is not null";
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) == SQLITE_OK)
    {
      while (sqlite3_step(stmt) == SQLITE_ROW)
      {
        VendorDiaWindow w;
        w.one_over_k0_start = sqlite3_column_double(stmt, 0);
        w.one_over_k0_end = sqlite3_column_double(stmt, 1);
        w.isolation_mz = sqlite3_column_double(stmt, 2);
        w.isolation_width = sqlite3_column_double(stmt, 3);
        // The method states them low-to-high, but a file that did not would
        // otherwise produce an empty band that rejects every peak in silence.
        if (w.one_over_k0_start > w.one_over_k0_end)
        {
          std::swap(w.one_over_k0_start, w.one_over_k0_end);
        }
        if (w.one_over_k0_start < w.one_over_k0_end) { out.push_back(w); }
      }
    }
    sqlite3_finalize(stmt);
    sqlite3_close(db);
    return out;
  }
}
