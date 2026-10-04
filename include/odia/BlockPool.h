// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// The extractor's chromatogram block pool. Internal to ChromatogramExtractor;
// it has a header of its own only so the test suite can drive take/give/release
// directly (test/tools/odia_extract_cases.cpp, case pool_capacity). Headers are
// not installed.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace ODIA
{

  /// Blocks, reused rather than returned to the allocator.
  ///
  /// Nearly every precursor asks for the same size -- 12 transitions by the
  /// window's cycle count -- so the free list hits essentially always,
  /// and the pass does one allocation per concurrently-live precursor for the
  /// whole run instead of one per precursor. That matters at 4.26 M
  /// precursors: 4.26 M allocate/free pairs of ~40 KB through malloc is where
  /// a general allocator fragments, and fragmentation would put back exactly
  /// the memory this change exists to remove.
  ///
  /// Blocks are handed out ZEROED, because the match accumulates into them
  /// and a reused block still holds the previous precursor's peaks.
  ///
  /// "The same size" holds within a stretch of the run, not across it. The
  /// cycles inside a fixed retention-time window fall as the cycle time
  /// lengthens, and on PXD047793 run 009 it drifts from 1.593 to 1.642 s
  /// (1,564 -> 1,531 cycles at +/-1255 s), so the chunks of a pass ask for
  /// sizes their predecessors never did. Kept for the whole pass, the pool
  /// then held every chunk's live set side by side -- eleven of them in pass
  /// 1, killed at 2.25 TB in the ninth. So the pass releases it at each chunk
  /// boundary, when nothing is live, and what it holds is one chunk's.
  ///
  /// Within a chunk the same drift still defeats an exact-size match, and a
  /// pass of few wide chunks -- pass 2, or any pass the cap does not cut --
  /// sees many sizes go by (pass 2 at -live_memory_gb 245 replayed up to
  /// ~567 GiB retained in one chunk). So a request takes the SMALLEST free
  /// block that holds it, and a block is filed under its capacity, not under
  /// the size it was last asked for. Only the first n cells are zeroed, and
  /// nothing reads past them: a precursor's cells are rows x its own cycles.
  ///
  /// Neither bounds what one chunk retains. A request larger than every free
  /// block allocates, so strictly growing sizes within a chunk (N, N+1, N+2,
  /// ... with each freed before the next is taken) retain their sum while
  /// only the largest is live.
  ///
  /// The backing arrays are counted as they are allocated and as they are
  /// DESTROYED -- by their deleter, i.e. when the memory goes back to the
  /// allocator, not by release()'s own bookkeeping -- so a release that reset
  /// its counters but kept the arrays shows up in blocksLeftAfterRelease().
  /// The deleter carries a pointer to the count: 8 B per backing array beside
  /// a block of ~40-74 KB.
  class BlockPool
  {
  public:
    BlockPool() = default;
    /// Not copyable, and so not movable: every owned array's deleter points
    /// at this pool's count.
    BlockPool(const BlockPool&) = delete;
    BlockPool& operator=(const BlockPool&) = delete;

    float* take(std::size_t n)
    {
      // A precursor all of whose transitions lack a product m/z has no rows.
      // It is still emitted -- a consumer counting precursors that yielded
      // nothing must see it -- so it needs a base that is not null.
      if (n == 0) { return &dummy_; }
      float* p = nullptr;
      const auto fit = free_.lower_bound(n);
      if (fit != free_.end())
      {
        p = fit->second;
        free_.erase(fit);
      }
      else
      {
        // Default-initialised, as make_unique_for_overwrite<float[]> was: the
        // fill below is the only write before use.
        Owned block(new float[n], Destroy{&freed_});
        ++allocated_;
        owned_.push_back(std::move(block));
        p = owned_.back().get();
        capacity_.emplace(p, n);
        held_ += n;
        reserved_ = std::max(reserved_, held_);
        blocks_peak_ = std::max(blocks_peak_, allocated_ - freed_);
      }
      std::fill_n(p, n, 0.0f);
      live_ += n;
      peak_ = std::max(peak_, live_);
      return p;
    }

    void give(float* p, std::size_t n)
    {
      if (n == 0) { return; }
      free_.emplace(capacity_.at(p), p);
      live_ -= n;
    }

    /// Frees every block back to the allocator and keeps the statistics.
    /// Only when none is live: a live block is a chromatogram still being
    /// filled, so releasing one is a defect and is refused with the number.
    void release()
    {
      if (live_ != 0)
      {
        throw std::logic_error("the block pool was released with " +
                               std::to_string(live_) + " points still live");
      }
      free_.clear();
      capacity_.clear();
      owned_.clear();
      held_ = 0;
      // Measured after the fact, from the deleters' count, not from what this
      // function believes it did.
      left_after_release_ = std::max(left_after_release_, allocated_ - freed_);
      between_releases_max_ = std::max(between_releases_max_, allocated_ - allocated_at_release_);
      allocated_at_release_ = allocated_;
    }

    /// Peak pool-owned float capacity: the most floats the pool's backing
    /// arrays held at once, live and free. Not RSS -- freed arrays go back to
    /// the allocator, which may keep their pages resident, and the pool's own
    /// map nodes and the allocation overhead are not in it. peakPoints() is
    /// the live part of it, in requested points.
    std::uint64_t reservedPoints() const { return reserved_; }
    std::uint64_t peakPoints() const { return peak_; }

    /// Backing arrays allocated so far.
    std::uint64_t blocksAllocated() const { return allocated_; }
    /// Backing arrays destroyed so far, counted by their deleter.
    std::uint64_t blocksFreed() const { return freed_; }
    /// The most backing arrays alive at once.
    std::uint64_t blocksPeak() const { return blocks_peak_; }
    /// The most backing arrays allocated between two releases (or since the
    /// last one): with a release per chunk, the largest chunk's.
    std::uint64_t blocksBetweenReleasesMax() const
    {
      return std::max(between_releases_max_, allocated_ - allocated_at_release_);
    }
    /// The most backing arrays still alive right after a release(): 0 when
    /// every release destroyed every array.
    std::uint64_t blocksLeftAfterRelease() const { return left_after_release_; }

  private:
    struct Destroy
    {
      std::uint64_t* freed = nullptr;
      void operator()(float* p) const noexcept
      {
        delete[] p;
        ++*freed;
      }
    };
    using Owned = std::unique_ptr<float[], Destroy>;

    // Declared before owned_, so they outlive it: the pool's destructor
    // destroys owned_ first, and its deleters still count into freed_.
    std::uint64_t allocated_ = 0, freed_ = 0, blocks_peak_ = 0;
    std::uint64_t left_after_release_ = 0, between_releases_max_ = 0, allocated_at_release_ = 0;
    /// Free blocks by CAPACITY, so lower_bound(n) is the best fit.
    std::multimap<std::size_t, float*> free_;
    std::unordered_map<const float*, std::size_t> capacity_;
    std::vector<Owned> owned_;
    std::uint64_t live_ = 0, peak_ = 0, held_ = 0, reserved_ = 0;
    float dummy_ = 0.0f;
  };

} // namespace ODIA
