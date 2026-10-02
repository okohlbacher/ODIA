// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/ResourceProbe.h>

#include <arrow/memory_pool.h>

#include <chrono>
#include <cstdio>
#include <algorithm>
#include <cstring>

#include <dlfcn.h>
#if defined(__GLIBC__)
#include <gnu/libc-version.h>
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

    /// glibc's `struct mallinfo2` (malloc.h, glibc >= 2.33), declared here
    /// because the sysroot's malloc.h is 2.28 and has no such struct. Ten
    /// size_t fields in this order are the ABI; it has not changed since 2.33.
    struct Mallinfo2
    {
      std::size_t arena, ordblks, smblks, hblks, hblkhd, usmblks, fsmblks,
                  uordblks, fordblks, keepcost;
    };
    using Mallinfo2Fn = Mallinfo2 (*)();

    /// Resolved once, from whatever libc the process actually runs on. Null
    /// when that libc has none -- then the line says so, with its version.
    Mallinfo2Fn mallinfo2Fn()
    {
      static const Mallinfo2Fn fn =
        reinterpret_cast<Mallinfo2Fn>(::dlsym(RTLD_DEFAULT, "mallinfo2"));
      return fn;
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
    // and every number of interest here is in the hundreds of GiB. Looked up
    // at run time, not tested at compile time -- see the header.
#if defined(__GLIBC__)
    m.libc_version = ::gnu_get_libc_version();
#endif
    if (const Mallinfo2Fn fn = mallinfo2Fn())
    {
      const auto t0 = std::chrono::steady_clock::now();
      const Mallinfo2 mi = fn();
      m.glibc_arena = mi.arena;
      m.glibc_mmap = mi.hblkhd;
      m.glibc_in_use = mi.uordblks;
      m.glibc_free = mi.fordblks;
      m.glibc_keepcost = mi.keepcost;
      m.glibc_ok = true;
      m.glibc_ms = msSince(t0);
    }

    // Arrow's default pool. jemalloc's stats only when jemalloc IS that pool:
    // on any other backend the first mallctl would initialise a dormant
    // jemalloc inside the process being measured (measured: +0.009 GiB
    // resident on a run whose pool is mimalloc) and report it as an owner.
    {
      const auto t0 = std::chrono::steady_clock::now();
      if (arrow::MemoryPool* pool = arrow::default_memory_pool())
      {
        m.arrow_backend = pool->backend_name();
        m.arrow_pool_bytes = pool->bytes_allocated();
        m.arrow_pool_peak = pool->max_memory();
      }
      if (m.arrow_backend == "jemalloc")
      {
        const auto stat = [](const char* name) -> std::int64_t {
          const auto r = arrow::jemalloc_get_stat(name);
          return r.ok() ? *r : -1;
        };
        m.je_allocated = stat("stats.allocated");
        m.je_active = stat("stats.active");
        m.je_resident = stat("stats.resident");
        m.je_retained = stat("stats.retained");
      }
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
                    "glibc malloc (libc %s; every new/malloc of ODIA, OpenMS, libstdc++; "
                    "NOT the arrow pool) arena %.3f GiB, mmap %.3f GiB, in use %.3f GiB, "
                    "free held %.3f GiB, keepcost %.3f GiB [%.1f ms]",
                    m.libc_version.empty() ? "?" : m.libc_version.c_str(),
                    double(m.glibc_arena) / 1073741824.0, double(m.glibc_mmap) / 1073741824.0,
                    double(m.glibc_in_use) / 1073741824.0, double(m.glibc_free) / 1073741824.0,
                    double(m.glibc_keepcost) / 1073741824.0, m.glibc_ms);
    }
    else
    {
      std::snprintf(glibc, sizeof glibc, "glibc mallinfo2 absent from the running libc (%s)",
                    m.libc_version.empty() ? "not glibc" : m.libc_version.c_str());
    }
    // -1 prints as -1.000: "not available", never confusable with an empty pool.
    if (m.je_allocated >= 0 || m.arrow_backend == "jemalloc")
    {
      std::snprintf(arrow, sizeof arrow,
                    "arrow pool '%s' %.3f GiB, pool peak %.3f GiB, jemalloc allocated %.3f GiB, "
                    "active %.3f GiB, resident %.3f GiB, retained %.3f GiB (-1 = unavailable) "
                    "[%.1f ms]",
                    m.arrow_backend.c_str(), b(m.arrow_pool_bytes), b(m.arrow_pool_peak),
                    b(m.je_allocated), b(m.je_active), b(m.je_resident), b(m.je_retained),
                    m.arrow_ms);
    }
    else
    {
      std::snprintf(arrow, sizeof arrow,
                    "arrow pool '%s' %.3f GiB, pool peak %.3f GiB (the backend's own "
                    "resident/retained totals are not exported; jemalloc not queried: it is "
                    "not the pool) [%.1f ms]",
                    m.arrow_backend.empty() ? "?" : m.arrow_backend.c_str(),
                    b(m.arrow_pool_bytes), b(m.arrow_pool_peak), m.arrow_ms);
    }
    // What neither allocator accounts for. glibc's arena is address space it
    // obtained, not all of it resident, so this is a LOWER bound on the
    // unattributed anonymous memory when positive and meaningless when the
    // glibc numbers are absent.
    char rest[160];
    if (m.rollup_ok && m.glibc_ok)
    {
      const double anon = double(m.anonymous_kib) * 1024.0;
      const double held = double(m.glibc_arena) + double(m.glibc_mmap) +
                          double(std::max<std::int64_t>(m.arrow_pool_bytes, 0));
      std::snprintf(rest, sizeof rest,
                    "anonymous minus (glibc arena+mmap + arrow pool) %.3f GiB",
                    (anon - held) / 1073741824.0);
    }
    else { std::snprintf(rest, sizeof rest, "anonymous not attributable (an input is missing)"); }
    return std::string(rollup) + "; " + glibc + "; " + arrow + "; " + rest;
  }

} // namespace ODIA
