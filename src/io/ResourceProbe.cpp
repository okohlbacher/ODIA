// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/ResourceProbe.h>

#include <arrow/memory_pool.h>

#include <chrono>
#include <cstdio>
#include <cstring>

#if defined(__GLIBC__)
#include <malloc.h>
#endif

namespace ODIA
{
  namespace
  {
    double msSince(std::chrono::steady_clock::time_point t0)
    {
      return std::chrono::duration<double, std::milli>(
               std::chrono::steady_clock::now() - t0).count();
    }
  } // namespace

  MemorySnapshot MemorySnapshot::now()
  {
    MemorySnapshot m;

    // smaps_rollup, not smaps: one record summed by the kernel instead of one
    // per mapping, which at this process's mapping count is the difference
    // between a page-table walk and a page-table walk plus a megabyte of text.
    {
      const auto t0 = std::chrono::steady_clock::now();
      if (std::FILE* f = std::fopen("/proc/self/smaps_rollup", "r"))
      {
        char line[256];
        const auto field = [&](const char* key, std::uint64_t& out) {
          const std::size_t n = std::strlen(key);
          unsigned long long v = 0;
          if (std::strncmp(line, key, n) == 0 && std::sscanf(line + n, "%llu", &v) == 1)
          {
            out = v;
            return true;
          }
          return false;
        };
        while (std::fgets(line, sizeof line, f) != nullptr)
        {
          if (field("Rss:", m.rss_kib) || field("Pss:", m.pss_kib) ||
              field("Pss_Anon:", m.pss_anon_kib) || field("Pss_File:", m.pss_file_kib) ||
              field("Anonymous:", m.anonymous_kib) || field("Swap:", m.swap_kib))
          { m.rollup_ok = true; }
        }
        std::fclose(f);
      }
      m.rollup_ms = msSince(t0);
    }

    // mallinfo2, not mallinfo: the int fields of the old call wrap at 2 GiB,
    // and every number of interest here is in the hundreds of GiB.
#if defined(__GLIBC__) && (__GLIBC__ > 2 || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 33))
    {
      const auto t0 = std::chrono::steady_clock::now();
      const struct mallinfo2 mi = ::mallinfo2();
      m.glibc_arena = mi.arena;
      m.glibc_mmap = mi.hblkhd;
      m.glibc_in_use = mi.uordblks;
      m.glibc_free = mi.fordblks;
      m.glibc_ok = true;
      m.glibc_ms = msSince(t0);
    }
#endif

    // Arrow's allocator. jemalloc_get_stat refreshes jemalloc's stats epoch
    // itself for the "stats.*" names; on a build without jemalloc it returns
    // NotImplemented and the fields stay -1.
    {
      const auto t0 = std::chrono::steady_clock::now();
      if (arrow::MemoryPool* pool = arrow::default_memory_pool())
      {
        m.arrow_backend = pool->backend_name();
        m.arrow_pool_bytes = pool->bytes_allocated();
      }
      const auto stat = [](const char* name) -> std::int64_t {
        const auto r = arrow::jemalloc_get_stat(name);
        return r.ok() ? *r : -1;
      };
      m.je_allocated = stat("stats.allocated");
      m.je_active = stat("stats.active");
      m.je_resident = stat("stats.resident");
      m.je_retained = stat("stats.retained");
      m.arrow_ms = msSince(t0);
    }
    return m;
  }

  std::string formatMemorySnapshot(const MemorySnapshot& m)
  {
    const auto kib = [](std::uint64_t v) { return double(v) / 1048576.0; };
    const auto b = [](std::int64_t v) { return v < 0 ? -1.0 : double(v) / 1073741824.0; };
    char rollup[320], glibc[320], arrow[384];
    if (m.rollup_ok)
    {
      std::snprintf(rollup, sizeof rollup,
                    "smaps_rollup rss %.3f GiB, pss %.3f GiB (anon %.3f, file %.3f), "
                    "anonymous %.3f GiB, swap %.3f GiB [%.1f ms]",
                    kib(m.rss_kib), kib(m.pss_kib), kib(m.pss_anon_kib), kib(m.pss_file_kib),
                    kib(m.anonymous_kib), kib(m.swap_kib), m.rollup_ms);
    }
    else { std::snprintf(rollup, sizeof rollup, "smaps_rollup unavailable"); }
    if (m.glibc_ok)
    {
      std::snprintf(glibc, sizeof glibc,
                    "glibc malloc (ODIA+OpenMS heap; Arrow's bundled jemalloc is NOT in it) "
                    "arena %.3f GiB, mmap %.3f GiB, in use %.3f GiB, free held %.3f GiB [%.1f ms]",
                    double(m.glibc_arena) / 1073741824.0, double(m.glibc_mmap) / 1073741824.0,
                    double(m.glibc_in_use) / 1073741824.0, double(m.glibc_free) / 1073741824.0,
                    m.glibc_ms);
    }
    else { std::snprintf(glibc, sizeof glibc, "glibc mallinfo2 unavailable"); }
    // -1 prints as -1.000: "not available", never confusable with an empty pool.
    std::snprintf(arrow, sizeof arrow,
                  "arrow pool '%s' %.3f GiB, jemalloc allocated %.3f GiB, active %.3f GiB, "
                  "resident %.3f GiB, retained %.3f GiB (-1 = unavailable) [%.1f ms]",
                  m.arrow_backend.empty() ? "?" : m.arrow_backend.c_str(),
                  b(m.arrow_pool_bytes), b(m.je_allocated), b(m.je_active),
                  b(m.je_resident), b(m.je_retained), m.arrow_ms);
    return std::string(rollup) + "; " + glibc + "; " + arrow;
  }

} // namespace ODIA
