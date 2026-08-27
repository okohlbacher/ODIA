#include <odia/VendorDiaWindows.h>

#include <sqlite3.h>
#include <zip.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
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
      // ENDS with, not contains. SQLite leaves `-wal` and `-journal` beside a
      // database, and `foo.diasqlite-wal` contains the name -- reading it as the
      // database yields garbage, or worse, silence.
      static const std::string leaf = "diaSettings.diasqlite";
      if (name.size() < leaf.size()) { return false; }
      if (name.compare(name.size() - leaf.size(), leaf.size(), leaf) != 0) { return false; }

      // Bruker leaves the pre-edit method in `backup-<date>.m/`. Anchored to a
      // path COMPONENT: a substring test also skips a live method someone named
      // `backup-final.m`, and does it silently.
      std::size_t at = 0;
      while (at < name.size())
      {
        const std::size_t end = name.find('/', at);
        const std::string comp = name.substr(at, end == std::string::npos ? std::string::npos
                                                                          : end - at);
        if (comp.rfind("backup-", 0) == 0) { return false; }
        if (end == std::string::npos) { break; }
        at = end + 1;
      }
      return true;
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

    // Count rather than take the first: two non-backup method entries mean the
    // archive cannot say which one was acquired under, and picking by archive
    // order is the silent-failure mode this whole change exists to remove.
    zip_int64_t found = -1;
    int matches = 0;
    const zip_int64_t n = zip_get_num_entries(z, 0);
    for (zip_int64_t i = 0; i < n; ++i)
    {
      const char* name = zip_get_name(z, i, 0);
      if (name && isLiveMethodEntry(name))
      {
        ++matches;
        if (found < 0) { found = i; }
      }
    }
    if (matches > 1)
    {
      std::cerr << "warning: " << matches << " live instrument methods in " << mzpeak_path
                << "; cannot tell which was acquired under, deriving mobility bands instead\n";
      zip_close(z);
      return out;
    }
    if (found < 0) { zip_close(z); return out; }

    zip_stat_t st;
    zip_stat_init(&st);
    if (zip_stat_index(z, found, 0, &st) != 0 || !(st.valid & ZIP_STAT_SIZE))
    {
      zip_close(z);
      return out;
    }

    // The real table is ~24 kB. A corrupt entry claiming gigabytes should not
    // become a bad_alloc in a reader that is allowed to return empty.
    if (st.size > (16u << 20))
    {
      std::cerr << "warning: instrument method entry is " << st.size
                << " bytes, refusing to read it\n";
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
      // mkstemp, not a name built from the pid: a predictable path in a shared
      // /tmp is a symlink target, and two sources constructed concurrently
      // would race on one name and read each other's bytes.
      std::string tpl = (dir / "odia_diasettings_XXXXXX").string();
      std::vector<char> buf(tpl.begin(), tpl.end());
      buf.push_back('\0');
      const int fd = ::mkstemp(buf.data());
      if (fd < 0) { return out; }
      tmp.path = std::string(buf.data());
      const ssize_t wrote = ::write(fd, blob.data(), blob.size());
      ::close(fd);
      if (wrote < 0 || static_cast<std::size_t>(wrote) != blob.size()) { return out; }
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
      "where Type = 1 and OneOverK0Start is not null and OneOverK0End is not null "
      "  and IsolationMz is not null "
      "order by IsolationMz";
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

    // A method may state the same isolation centre in several cycles. If those
    // rows carry DIFFERENT mobility bands the scheme varies across the run, and
    // matching by centre alone would hand every cycle whichever row came first.
    // Collapsing two real bands into one is exactly the failure this change
    // exists to prevent, so refuse instead of picking.
    for (std::size_t i = 0; i + 1 < out.size(); ++i)
    {
      for (std::size_t j = i + 1; j < out.size(); ++j)
      {
        if (std::abs(out[i].isolation_mz - out[j].isolation_mz) > 1e-6) { continue; }
        if (std::abs(out[i].one_over_k0_start - out[j].one_over_k0_start) > 1e-9 ||
            std::abs(out[i].one_over_k0_end - out[j].one_over_k0_end) > 1e-9)
        {
          std::cerr << "warning: instrument method states isolation centre "
                    << out[i].isolation_mz << " with two different mobility bands ("
                    << out[i].one_over_k0_start << '-' << out[i].one_over_k0_end << " and "
                    << out[j].one_over_k0_start << '-' << out[j].one_over_k0_end
                    << "); the scheme varies across the run and matching by centre cannot "
                       "tell them apart, so mobility bands are being DERIVED instead\n";
          return {};
        }
      }
    }
    // Duplicates that agree are harmless; keep one of each so the match loop is
    // over distinct windows.
    std::vector<VendorDiaWindow> uniq;
    for (const auto& w : out)
    {
      bool seen = false;
      for (const auto& u : uniq)
      {
        if (std::abs(u.isolation_mz - w.isolation_mz) <= 1e-6) { seen = true; break; }
      }
      if (!seen) { uniq.push_back(w); }
    }
    return uniq;
  }
}
