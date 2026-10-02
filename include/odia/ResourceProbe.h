// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <atomic>
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
  /// those into a printed term. NOTE what reconciles by construction: phases
  /// that telescope from one counter sum to that counter whatever they cover,
  /// so "phases == footer" checks nothing. The decomposition check is the
  /// nested one -- a phase against the sum of the STAGE brackets inside it.
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
  /// VmHWM is NOT monotone as read. The kernel prints max(mm->hiwater_rss,
  /// current RSS), and folds the current RSS into hiwater_rss only lazily (on
  /// munmap and a few other paths); pages released with madvise(DONTNEED) --
  /// how glibc trims non-main arenas and how mimalloc purges -- leave it
  /// behind. A read can therefore be LOWER than an earlier one, and the footer
  /// (ru_maxrss, the same hiwater_rss) can be lower than a VmRSS this process
  /// observed. StatusLedger keeps the maxima of every read so the closing line
  /// can say whether that happened. The first line on which the hwm reaches
  /// its final value dates the peak to that stage; whether the peak was
  /// RESIDENT at a boundary (rss within ~1% of it) or transient inside the
  /// stage is what the current VmRSS beside it says.
  struct ProcessStatus
  {
    std::uint64_t rss_kib = 0;
    std::uint64_t rss_anon_kib = 0;   ///< RssAnon: heap (every allocator), stacks
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
      record(st);
      return st;
    }

    /// Every read of /proc/self/status, summarised: how many, the largest
    /// VmRSS and VmHWM seen, and how often VmHWM read lower than the read
    /// before it. Process-wide; reads happen on the driver thread.
    struct Ledger
    {
      std::atomic<std::uint64_t> reads{0}, max_rss_kib{0}, max_hwm_kib{0};
      std::atomic<std::uint64_t> last_hwm_kib{0}, hwm_decreases{0}, max_hwm_drop_kib{0};
    };
    static Ledger& ledger()
    {
      static Ledger l;
      return l;
    }

  private:
    static void record(const ProcessStatus& st)
    {
      if (st.hwm_kib == 0) { return; }
      Ledger& l = ledger();
      l.reads.fetch_add(1);
      if (st.rss_kib > l.max_rss_kib.load()) { l.max_rss_kib.store(st.rss_kib); }
      if (st.hwm_kib > l.max_hwm_kib.load()) { l.max_hwm_kib.store(st.hwm_kib); }
      const std::uint64_t last = l.last_hwm_kib.exchange(st.hwm_kib);
      if (last != 0 && st.hwm_kib < last)
      {
        l.hwm_decreases.fetch_add(1);
        if (last - st.hwm_kib > l.max_hwm_drop_kib.load())
        { l.max_hwm_drop_kib.store(last - st.hwm_kib); }
      }
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
  /// (proportional set size), glibc malloc's view, and Arrow's memory pool's.
  ///
  /// Allocators. Everything ODIA, OpenMS and libstdc++ allocate with `new` or
  /// `malloc` goes through glibc. Arrow's buffers go through Arrow's DEFAULT
  /// POOL, whose backend is printed: on the conda libarrow 23 this build links
  /// it is 'mimalloc' (statically bundled; its own resident/retained totals
  /// are not exported by libarrow, so only the pool's bytes and peak are
  /// printed). libarrow also bundles jemalloc, but it is queried only when it
  /// IS the pool's backend -- the first mallctl initialises jemalloc and would
  /// add its own arenas to the process being measured. So the line separates
  /// `glibc arena+mmap` (what glibc took from the OS), `glibc in use` (what it
  /// handed out), the Arrow pool, and the rest of rss_anon (thread stacks,
  /// mimalloc's retained segments, anything mmapped directly).
  ///
  /// glibc's mallinfo2 is resolved at RUN time (dlsym): the conda sysroot the
  /// build compiles against is glibc 2.28, which predates mallinfo2 (2.33),
  /// so a compile-time test compiled it out on every env.sh build. The
  /// running libc (2.39 on the IBMI nodes) has it.
  ///
  /// Read at chunk ends and phase boundaries ONLY, never per bracket: the
  /// smaps_rollup read walks the page tables of the whole address space
  /// (~11-18 ms/GiB measured, i.e. seconds at 600 GiB), and mallinfo2 walks
  /// every free chunk of every arena under that arena's lock. Each part times
  /// itself, and the callers charge the whole probe to a PROBE ledger of its
  /// own rather than to the stage or phase that follows it.
  struct MemorySnapshot
  {
    // /proc/self/smaps_rollup, KiB.
    bool rollup_ok = false;
    std::uint64_t rss_kib = 0, pss_kib = 0, pss_anon_kib = 0, pss_file_kib = 0;
    std::uint64_t anonymous_kib = 0, swap_kib = 0;
    double rollup_ms = 0.0;

    // glibc mallinfo2, bytes, summed over all arenas. Absent (glibc_ok false)
    // only when the running libc has no mallinfo2.
    bool glibc_ok = false;
    std::uint64_t glibc_arena = 0;     ///< obtained from the OS via brk/heaps
    std::uint64_t glibc_mmap = 0;      ///< in individually mmapped chunks
    std::uint64_t glibc_in_use = 0;    ///< handed out and not freed
    std::uint64_t glibc_free = 0;      ///< free but still held
    std::uint64_t glibc_keepcost = 0;  ///< releasable from the top of the main arena
    double glibc_ms = 0.0;
    std::string libc_version;          ///< gnu_get_libc_version() of the RUNNING libc

    // Arrow's default memory pool: backend, bytes handed out now, the pool's
    // own peak. jemalloc's mallctl totals only when jemalloc IS the backend;
    // -1 otherwise ("not queried", never "zero").
    std::string arrow_backend;
    std::int64_t arrow_pool_bytes = -1, arrow_pool_peak = -1;
    std::int64_t je_allocated = -1, je_active = -1, je_resident = -1, je_retained = -1;
    double arrow_ms = 0.0;

    static MemorySnapshot now();
  };

  /// One line, all GiB:
  ///   smaps_rollup rss R, pss P (anon A, file F), anonymous N, swap S [t ms];
  ///   glibc malloc (...) arena X, mmap Y, in use U, free held H, keepcost K [t ms];
  ///   arrow pool '<backend>' B, pool peak P[, jemalloc ...] [t ms];
  ///   anonymous not held by glibc or the arrow pool R
  std::string formatMemorySnapshot(const MemorySnapshot& m);

} // namespace ODIA
