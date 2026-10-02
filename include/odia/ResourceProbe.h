// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>

#include <sys/resource.h>

namespace ODIA
{
  /// What the PROCESS has consumed so far: wall clock plus getrusage(RUSAGE_SELF).
  ///
  /// WHY THIS EXISTS. Every stage timer in the pipeline was wall-only, so every
  /// CPU, memory and system-time share in doc/83 §2.2-2.3 was obtained by
  /// SUBTRACTION -- "match runs at 17-21% efficiency" is 131,524 CPU-s minus
  /// two inferred serial terms, divided by a wall time; "~165 GiB remainder" is
  /// a footer maximum minus two other maxima that need not coincide in time.
  /// A bracket that reads the same four counters the footer reads turns each of
  /// those into a printed term, and a sum of brackets that fails to reconcile
  /// with the footer says the decomposition is missing a stage.
  ///
  /// RUSAGE_SELF, not RUSAGE_THREAD, deliberately: a stage's cost is whatever
  /// the process spent while it ran, including the worker pool during match and
  /// OpenMP inside finish(). For a stage that is serial by construction (the
  /// sink, the MS1 build) the two are the same up to background threads -- and
  /// a CPU/wall ratio above 1 on such a stage is then itself the finding.
  ///
  /// Caveat that matters when reading the output: the kernel's user/sys SPLIT
  /// is tick-sampled and scaled to the precise total, so user+sys is accurate
  /// to the scheduler's runtime but the split of a sub-second bracket is not.
  /// Faults are exact counts. Every read is a syscall of a few microseconds;
  /// brackets are therefore placed per block, batch or chunk -- never per
  /// precursor, where 2 x 9.25 M reads a pass would be tens of seconds of the
  /// very sink time being measured.
  struct ResourceSample
  {
    double wall = 0.0;          ///< steady clock, seconds since its own epoch
    double user = 0.0;          ///< seconds, all threads of the process
    double sys = 0.0;
    std::uint64_t minflt = 0;   ///< minor (reclaim) page faults
    std::uint64_t majflt = 0;   ///< major (I/O) page faults
    std::uint64_t maxrss_kib = 0;  ///< ru_maxrss: what a `time -v` footer reports

    static ResourceSample now()
    {
      ResourceSample s;
      s.wall = std::chrono::duration<double>(
                 std::chrono::steady_clock::now().time_since_epoch()).count();
      rusage ru{};
      if (getrusage(RUSAGE_SELF, &ru) == 0)
      {
        s.user = double(ru.ru_utime.tv_sec) + double(ru.ru_utime.tv_usec) * 1e-6;
        s.sys = double(ru.ru_stime.tv_sec) + double(ru.ru_stime.tv_usec) * 1e-6;
        s.minflt = static_cast<std::uint64_t>(ru.ru_minflt);
        s.majflt = static_cast<std::uint64_t>(ru.ru_majflt);
        s.maxrss_kib = static_cast<std::uint64_t>(ru.ru_maxrss);
      }
      return s;
    }
  };

  /// Resident set, its split, its high-water mark and the thread count, from
  /// /proc/self/status. Zero where the file is absent (not Linux).
  ///
  /// VmHWM is monotone, so the first line on which it reaches the footer's
  /// peak dates the peak to that stage -- no arm at >=300,000 MB has ever had
  /// that. But HWM alone cannot say which LATER stage still owns the pages
  /// (codex review Q6.1): that is what the CURRENT VmRSS and its anon/file
  /// split are for. A stage that ends with rss far below hwm released what the
  /// peak was made of; one that ends at hwm still holds it.
  struct ProcessStatus
  {
    std::uint64_t rss_kib = 0;
    std::uint64_t rss_anon_kib = 0;   ///< RssAnon: heap, both allocators
    std::uint64_t rss_file_kib = 0;   ///< RssFile: mapped files, libraries
    std::uint64_t hwm_kib = 0;
    unsigned threads = 0;

    static ProcessStatus now()
    {
      ProcessStatus st;
      std::FILE* f = std::fopen("/proc/self/status", "r");
      if (f == nullptr) { return st; }
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
        std::uint64_t thr = 0;
        if (field("VmRSS:", st.rss_kib) || field("VmHWM:", st.hwm_kib) ||
            field("RssAnon:", st.rss_anon_kib) || field("RssFile:", st.rss_file_kib))
        { continue; }
        if (field("Threads:", thr)) { st.threads = static_cast<unsigned>(thr); }
      }
      std::fclose(f);
      return st;
    }
  };

  /// The summed cost of one stage over every bracket that was charged to it.
  struct StageCost
  {
    double wall = 0.0, user = 0.0, sys = 0.0;
    std::uint64_t minflt = 0, majflt = 0;
    std::uint64_t brackets = 0;

    double cpu() const { return user + sys; }

    void add(const ResourceSample& from, const ResourceSample& to)
    {
      wall += to.wall - from.wall;
      user += to.user - from.user;
      sys += to.sys - from.sys;
      minflt += to.minflt - from.minflt;
      majflt += to.majflt - from.majflt;
      ++brackets;
    }

    StageCost& operator+=(const StageCost& o)
    {
      wall += o.wall; user += o.user; sys += o.sys;
      minflt += o.minflt; majflt += o.majflt; brackets += o.brackets;
      return *this;
    }
  };

  /// The one line format every STAGE and PHASE line shares, so a single regex
  /// reads all of them:
  ///
  ///   wall W s, cpu C s (user U, sys S), cpu/wall R, minflt F, majflt M,
  ///   rss X GiB (anon A, file F), hwm Y GiB, threads T
  ///
  /// rss/hwm/threads are the state at the END of the stage (for a summed
  /// stage, at the moment it is printed), not averages over it.
  inline std::string formatStageCost(const StageCost& c, const ProcessStatus& s)
  {
    const auto gib = [](std::uint64_t kib) { return double(kib) / 1048576.0; };
    char buf[448];
    std::snprintf(buf, sizeof buf,
                  "wall %.3f s, cpu %.3f s (user %.3f, sys %.3f), cpu/wall %.3f, "
                  "minflt %llu, majflt %llu, rss %.3f GiB (anon %.3f, file %.3f), "
                  "hwm %.3f GiB, threads %u",
                  c.wall, c.cpu(), c.user, c.sys, c.wall > 0.0 ? c.cpu() / c.wall : 0.0,
                  static_cast<unsigned long long>(c.minflt),
                  static_cast<unsigned long long>(c.majflt),
                  gib(s.rss_kib), gib(s.rss_anon_kib), gib(s.rss_file_kib),
                  gib(s.hwm_kib), s.threads);
    return buf;
  }

  /// Who owns the resident memory at a stage BOUNDARY: the kernel's view
  /// (proportional set size), glibc's view and Arrow's jemalloc's view.
  ///
  /// The process runs TWO allocators. libarrow bundles jemalloc (doc/83 §2.2;
  /// `jemalloc_bg_thd` is visible inside the ODIA process), and the decode
  /// path's Arrow buffers go through it; everything ODIA and OpenMS allocate
  /// with `new` goes through glibc. So glibc's mallinfo2 covers only PART of
  /// the heap, and the line says so: `glibc in use + jemalloc allocated` is
  /// the heap both allocators hand out, `rss_anon - both` is what they hold
  /// without handing it out (free lists, retained arenas) plus anything else
  /// anonymous (thread stacks, OpenMP).
  ///
  /// Read at chunk ends and phase boundaries ONLY, never per bracket: the
  /// smaps_rollup read walks the page tables of the whole address space
  /// (~157 M PTEs at 600 GiB), and mallinfo2 walks every free chunk of every
  /// arena under that arena's lock. Each part times itself and prints it, so
  /// what the probe cost is on the line beside what it measured.
  struct MemorySnapshot
  {
    // /proc/self/smaps_rollup, KiB.
    bool rollup_ok = false;
    std::uint64_t rss_kib = 0, pss_kib = 0, pss_anon_kib = 0, pss_file_kib = 0;
    std::uint64_t anonymous_kib = 0, swap_kib = 0;
    double rollup_ms = 0.0;

    // glibc mallinfo2, bytes, summed over all arenas.
    bool glibc_ok = false;
    std::uint64_t glibc_arena = 0;     ///< obtained from the OS via brk/heaps
    std::uint64_t glibc_mmap = 0;      ///< in individually mmapped chunks
    std::uint64_t glibc_in_use = 0;    ///< handed out and not freed
    std::uint64_t glibc_free = 0;      ///< free but still held
    double glibc_ms = 0.0;

    // Arrow's default memory pool and, if it is jemalloc, jemalloc's own
    // totals (mallctl "stats.*", refreshed by Arrow). -1 where unavailable.
    std::string arrow_backend;
    std::int64_t arrow_pool_bytes = -1;
    std::int64_t je_allocated = -1, je_active = -1, je_resident = -1, je_retained = -1;
    double arrow_ms = 0.0;

    static MemorySnapshot now();
  };

  /// One line, all GiB:
  ///   smaps_rollup rss R, pss P (anon A, file F), anonymous N, swap S [t ms];
  ///   glibc (...) arena X, mmap Y, in use U, free held H [t ms];
  ///   arrow <backend> pool B, jemalloc allocated J, active .., resident ..,
  ///   retained .. [t ms]
  std::string formatMemorySnapshot(const MemorySnapshot& m);

} // namespace ODIA
