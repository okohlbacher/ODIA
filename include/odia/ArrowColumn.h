// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <arrow/array.h>
#include <arrow/table.h>

#include <cstdint>
#include <memory>
#include <string>
#include <string_view>

namespace ODIA
{

  /// A cursor over one column of an arrow::Table that never assumes chunk(0).
  ///
  /// An Arrow StringArray addresses its characters with 32-bit offsets, so it
  /// cannot hold more than 2 GB of them. A proteome-scale transition table
  /// exceeds that, which means a large text column is *necessarily* written as
  /// several chunks -- the normal state of the data, not an edge case. Reading
  /// chunk(0) while looping to num_rows walks off the end of the column; in the
  /// predecessor project that died inside basic_string::_M_create, and the
  /// milder version is a library that silently loads short.
  ///
  /// The cursor caches the chunk it last resolved, so a sequential scan costs
  /// O(1) amortised per row rather than O(chunks).
  ///
  /// resolve() deliberately returns a raw pointer. Its shared_ptr equivalent is
  /// called once per (row, column) -- ~1e9 times on the benchmark library -- and
  /// every call is an atomic pair on a control block shared by all threads, so
  /// the cost *grows* with thread count: a measured 20.1 s serial became 34.5 s
  /// on 64 threads. The pointer is valid as long as the owning table is.
  /// Not thread-safe: the cursor is shared mutable state, so two threads
  /// scanning one column return each other's rows rather than colliding
  /// visibly. One cursor per thread.
  class ChunkedColumn
  {
  public:
    ChunkedColumn() = default;

    /// @param column may be null, which leaves the cursor invalid; callers that
    ///        require the column must check.
    explicit ChunkedColumn(std::shared_ptr<arrow::ChunkedArray> column);

    bool valid() const { return column_ != nullptr; }
    std::int64_t length() const { return column_ == nullptr ? 0 : column_->length(); }
    arrow::Type::type typeId() const;

    /// Resolves a global row index to its chunk and the index within it.
    /// Throws std::out_of_range for a row past the end.
    const arrow::Array* resolve(std::int64_t row, std::int64_t& index_in_chunk) const;

    bool isNull(std::int64_t row) const;

    /// Numeric value of a row, whatever the column's own arithmetic type.
    ///
    /// Widths and signedness vary between writers for the same logical column
    /// (charge has been seen as int32 and as int8), and refusing them would
    /// reject valid files. Values that are not numeric at all are refused
    /// rather than coerced.
    double getDouble(std::int64_t row, double if_null) const;
    std::int64_t getInt64(std::int64_t row, std::int64_t if_null) const;
    bool getBool(std::int64_t row, bool if_null) const;

    /// The characters of a row. Valid while the owning table is.
    std::string_view getString(std::int64_t row, std::string_view if_null = {}) const;

    const std::string& name() const { return name_; }
    void setName(std::string name) { name_ = std::move(name); }

  private:
    std::shared_ptr<arrow::ChunkedArray> column_;
    std::string name_;
    // Cursor state. Mutable so the accessors stay const: the cache is not part
    // of the value, and a non-const accessor would force every caller holding
    // the table by const reference to copy.
    mutable int chunk_ = 0;
    mutable std::int64_t chunk_begin_ = 0;
    mutable std::int64_t chunk_end_ = 0;
  };

  /// Look a column up by name. Returns an invalid cursor if it is absent, so a
  /// caller can distinguish an optional column from a required one.
  ChunkedColumn column(const std::shared_ptr<arrow::Table>& table, const std::string& name);

  /// Look a column up by name, throwing if it is absent.
  ChunkedColumn requiredColumn(const std::shared_ptr<arrow::Table>& table,
                               const std::string& name, const std::string& context);

} // namespace ODIA
