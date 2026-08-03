// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// ChunkedColumn's whole reason to exist is that a column can arrive in several
// chunks, so that is what has to be tested -- and it cannot be tested through a
// Parquet fixture. Arrow's Parquet reader concatenates row groups: a file
// written with one row group per row still comes back as a single chunk, which
// is why a deliberate "chunk(0) is the whole column" defect survived a fixture
// built that way. In production the split is forced by a StringArray's 32-bit
// offsets overflowing at 2 GB of characters, which no test can afford to
// reproduce.
//
// So the chunking is constructed directly here. The values are a known
// function of the global row index, so reading the right row of the wrong
// chunk gives a wrong value rather than merely a different one.

#include <odia/ArrowColumn.h>

#include <arrow/api.h>

#include <iostream>
#include <string>
#include <vector>

namespace
{
  int failures = 0;

  void fail(const std::string& what)
  {
    ++failures;
    std::cout << "  FAIL " << what << "\n";
  }

  std::int64_t expectedInt(std::int64_t row) { return row * 7 + 3; }
  std::string expectedStr(std::int64_t row) { return "row-" + std::to_string(row); }

  /// A column split into @p sizes, carrying values by global row index. Empty
  /// chunks are included by the callers because Arrow produces them and a
  /// cursor that assumes non-empty chunks skips a row or loops.
  ODIA::ChunkedColumn makeInt(const std::vector<std::int64_t>& sizes)
  {
    std::vector<std::shared_ptr<arrow::Array>> chunks;
    std::int64_t row = 0;
    for (const auto n : sizes)
    {
      arrow::Int64Builder b;
      for (std::int64_t i = 0; i < n; ++i) { (void)b.Append(expectedInt(row++)); }
      std::shared_ptr<arrow::Array> a;
      (void)b.Finish(&a);
      chunks.push_back(a);
    }
    ODIA::ChunkedColumn c(std::make_shared<arrow::ChunkedArray>(chunks, arrow::int64()));
    c.setName("ints");
    return c;
  }

  ODIA::ChunkedColumn makeStr(const std::vector<std::int64_t>& sizes)
  {
    std::vector<std::shared_ptr<arrow::Array>> chunks;
    std::int64_t row = 0;
    for (const auto n : sizes)
    {
      arrow::StringBuilder b;
      for (std::int64_t i = 0; i < n; ++i) { (void)b.Append(expectedStr(row++)); }
      std::shared_ptr<arrow::Array> a;
      (void)b.Finish(&a);
      chunks.push_back(a);
    }
    ODIA::ChunkedColumn c(std::make_shared<arrow::ChunkedArray>(chunks, arrow::utf8()));
    c.setName("strings");
    return c;
  }

  std::int64_t total(const std::vector<std::int64_t>& sizes)
  {
    std::int64_t n = 0;
    for (const auto s : sizes) { n += s; }
    return n;
  }

  void checkLayout(const std::vector<std::int64_t>& sizes, const std::string& label)
  {
    const auto n = total(sizes);
    auto ints = makeInt(sizes);
    auto strs = makeStr(sizes);

    if (ints.length() != n) { fail(label + ": length " + std::to_string(ints.length())); }

    // Forward, the access pattern every caller has.
    for (std::int64_t r = 0; r < n; ++r)
    {
      if (ints.getInt64(r, -1) != expectedInt(r))
      {
        fail(label + ": forward row " + std::to_string(r) + " gave " +
             std::to_string(ints.getInt64(r, -1)));
        break;
      }
      if (strs.getString(r) != expectedStr(r))
      {
        fail(label + ": forward string row " + std::to_string(r) + " gave " +
             std::string(strs.getString(r)));
        break;
      }
    }

    // Backward, which forces the cursor to reset rather than walk.
    for (std::int64_t r = n - 1; r >= 0; --r)
    {
      if (ints.getInt64(r, -1) != expectedInt(r))
      {
        fail(label + ": backward row " + std::to_string(r) + " gave " +
             std::to_string(ints.getInt64(r, -1)));
        break;
      }
    }

    // Jumping, in case the cursor is only correct for adjacent rows.
    for (std::int64_t step : {1, 2, 3, 5, 7, 11})
    {
      for (std::int64_t r = 0; r < n; r += step)
      {
        if (ints.getInt64(r, -1) != expectedInt(r))
        {
          fail(label + ": step " + std::to_string(step) + " row " + std::to_string(r));
          break;
        }
      }
    }

    // Past the end must throw. Returning the last row instead would let a
    // reader that miscounts silently duplicate it across a whole library.
    bool threw = false;
    try { (void)ints.getInt64(n, -1); }
    catch (const std::out_of_range&) { threw = true; }
    if (!threw) { fail(label + ": row " + std::to_string(n) + " did not throw"); }

    // A negative row must throw too. The chunk walk cannot catch it -- it is
    // below the first chunk, not above the last -- so it would resolve to
    // chunk 0 at a negative index and read outside the buffer.
    threw = false;
    try { (void)ints.getInt64(-1, -1); }
    catch (const std::out_of_range&) { threw = true; }
    if (!threw) { fail(label + ": row -1 did not throw"); }
  }
} // namespace

int main()
{
  // One chunk: the shape Arrow actually returns for a Parquet file that fits.
  checkLayout({10}, "single");
  // Many equal chunks.
  checkLayout({4, 4, 4}, "equal");
  // Uneven, with the interesting boundaries at 1.
  checkLayout({1, 5, 1, 9, 1}, "uneven");
  // Empty chunks at the front, in the middle and at the end.
  checkLayout({0, 3, 0, 0, 4, 0}, "with-empties");
  // A single row split so that chunk 0 is empty -- the case where "chunk(0) is
  // the whole column" does not merely read the wrong row but reads nothing.
  checkLayout({0, 1}, "empty-first");

  // An absent column must answer with the caller's default rather than throw,
  // because optional columns are read through the same accessors.
  ODIA::ChunkedColumn absent;
  if (absent.valid()) { fail("a default-constructed column claims to be valid"); }
  if (absent.length() != 0) { fail("a default-constructed column has length"); }
  if (absent.getInt64(0, 42) != 42) { fail("an absent column ignored the default"); }
  if (!absent.isNull(0)) { fail("an absent column reported a value as present"); }

  std::cout << "chunked column: " << failures << " failures\n";
  return failures == 0 ? 0 : 1;
}
