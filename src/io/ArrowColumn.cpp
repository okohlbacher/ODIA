// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/ArrowColumn.h>

#include <stdexcept>
#include <string>

namespace ODIA
{

  namespace
  {
    /// Resolve a dictionary-encoded cell to its underlying array and index.
    ///
    /// Dictionary encoding is what any Arrow-native writer reaches for on a
    /// column like `type` or `annotation`, where a handful of values repeat
    /// across millions of rows -- exactly the columns this format has. It also
    /// arrives as one chunk per Parquet row group rather than one chunk for
    /// the file, which makes it the only way to reach the multi-chunk path at
    /// test scale.
    const arrow::Array* undictionary(const arrow::Array* a, std::int64_t& i)
    {
      if (a->type_id() != arrow::Type::DICTIONARY) { return a; }
      const auto* d = static_cast<const arrow::DictionaryArray*>(a);
      const auto index = d->GetValueIndex(i);
      i = index;
      return d->dictionary().get();
    }

    [[noreturn]] void badType(const std::string& name, const arrow::Array& a,
                              const char* wanted)
    {
      throw std::runtime_error("column '" + name + "' has type " + a.type()->ToString() +
                               ", which is not " + wanted);
    }
  }

  ChunkedColumn::ChunkedColumn(std::shared_ptr<arrow::ChunkedArray> column)
    : column_(std::move(column))
  {
    if (column_ != nullptr && column_->num_chunks() > 0)
    {
      chunk_end_ = column_->chunk(0)->length();
    }
  }

  arrow::Type::type ChunkedColumn::typeId() const
  {
    if (column_ == nullptr) { return arrow::Type::NA; }
    return column_->type()->id();
  }

  const arrow::Array* ChunkedColumn::resolve(std::int64_t row,
                                             std::int64_t& index_in_chunk) const
  {
    if (column_ == nullptr) { throw std::out_of_range("resolve on an absent column"); }
    if (row < 0 || row >= column_->length())
    {
      throw std::out_of_range("row " + std::to_string(row) + " outside column '" + name_ +
                              "' of length " + std::to_string(column_->length()));
    }

    // Sequential scans, which is every caller, stay in the cached chunk. The
    // pair (chunk_begin_, chunk_end_) always describes chunk_ itself, so
    // advancing means stepping the index first and then the bounds.
    if (row < chunk_begin_ || row >= chunk_end_)
    {
      if (row < chunk_begin_)
      {
        chunk_ = 0;
        chunk_begin_ = 0;
        chunk_end_ = column_->chunk(0)->length();
      }
      // Walking forward from the cached chunk keeps a scan linear. Empty chunks
      // are legal and the loop steps over them without a special case.
      while (row >= chunk_end_)
      {
        ++chunk_;
        if (chunk_ >= column_->num_chunks())
        {
          throw std::out_of_range("row " + std::to_string(row) + " past the last chunk of '" +
                                  name_ + "'");
        }
        chunk_begin_ = chunk_end_;
        chunk_end_ = chunk_begin_ + column_->chunk(chunk_)->length();
      }
    }
    index_in_chunk = row - chunk_begin_;
    return column_->chunk(chunk_).get();
  }

  bool ChunkedColumn::isNull(std::int64_t row) const
  {
    if (column_ == nullptr) { return true; }
    std::int64_t i = 0;
    return resolve(row, i)->IsNull(i);
  }

  double ChunkedColumn::getDouble(std::int64_t row, double if_null) const
  {
    if (column_ == nullptr) { return if_null; }
    std::int64_t i = 0;
    const arrow::Array* a = resolve(row, i);
    if (a->IsNull(i)) { return if_null; }
    a = undictionary(a, i);
    switch (a->type_id())
    {
      case arrow::Type::DOUBLE: return static_cast<const arrow::DoubleArray*>(a)->Value(i);
      case arrow::Type::FLOAT:  return static_cast<const arrow::FloatArray*>(a)->Value(i);
      case arrow::Type::INT64:  return static_cast<double>(static_cast<const arrow::Int64Array*>(a)->Value(i));
      case arrow::Type::INT32:  return static_cast<double>(static_cast<const arrow::Int32Array*>(a)->Value(i));
      case arrow::Type::INT16:  return static_cast<double>(static_cast<const arrow::Int16Array*>(a)->Value(i));
      case arrow::Type::INT8:   return static_cast<double>(static_cast<const arrow::Int8Array*>(a)->Value(i));
      case arrow::Type::UINT64: return static_cast<double>(static_cast<const arrow::UInt64Array*>(a)->Value(i));
      case arrow::Type::UINT32: return static_cast<double>(static_cast<const arrow::UInt32Array*>(a)->Value(i));
      case arrow::Type::UINT16: return static_cast<double>(static_cast<const arrow::UInt16Array*>(a)->Value(i));
      case arrow::Type::UINT8:  return static_cast<double>(static_cast<const arrow::UInt8Array*>(a)->Value(i));
      case arrow::Type::BOOL:   return static_cast<const arrow::BooleanArray*>(a)->Value(i) ? 1.0 : 0.0;
      default: badType(name_, *a, "numeric");
    }
  }

  std::int64_t ChunkedColumn::getInt64(std::int64_t row, std::int64_t if_null) const
  {
    if (column_ == nullptr) { return if_null; }
    std::int64_t i = 0;
    const arrow::Array* a = resolve(row, i);
    if (a->IsNull(i)) { return if_null; }
    a = undictionary(a, i);
    switch (a->type_id())
    {
      case arrow::Type::INT64:  return static_cast<const arrow::Int64Array*>(a)->Value(i);
      case arrow::Type::INT32:  return static_cast<const arrow::Int32Array*>(a)->Value(i);
      case arrow::Type::INT16:  return static_cast<const arrow::Int16Array*>(a)->Value(i);
      case arrow::Type::INT8:   return static_cast<const arrow::Int8Array*>(a)->Value(i);
      case arrow::Type::UINT64: return static_cast<std::int64_t>(static_cast<const arrow::UInt64Array*>(a)->Value(i));
      case arrow::Type::UINT32: return static_cast<const arrow::UInt32Array*>(a)->Value(i);
      case arrow::Type::UINT16: return static_cast<const arrow::UInt16Array*>(a)->Value(i);
      case arrow::Type::UINT8:  return static_cast<const arrow::UInt8Array*>(a)->Value(i);
      case arrow::Type::BOOL:   return static_cast<const arrow::BooleanArray*>(a)->Value(i) ? 1 : 0;
      default: badType(name_, *a, "integral");
    }
  }

  bool ChunkedColumn::getBool(std::int64_t row, bool if_null) const
  {
    if (column_ == nullptr) { return if_null; }
    std::int64_t i = 0;
    const arrow::Array* a = resolve(row, i);
    if (a->IsNull(i)) { return if_null; }
    a = undictionary(a, i);
    if (a->type_id() == arrow::Type::BOOL)
    {
      return static_cast<const arrow::BooleanArray*>(a)->Value(i);
    }
    // Writers disagree about the decoy flag's type: DIA-NN spells it as an
    // integer, OpenSWATH as a bool. Anything non-zero is a decoy. Read from the
    // already-resolved array rather than calling getInt64, which would resolve
    // the row again and undo the dictionary hop above.
    switch (a->type_id())
    {
      case arrow::Type::INT64:  return static_cast<const arrow::Int64Array*>(a)->Value(i) != 0;
      case arrow::Type::INT32:  return static_cast<const arrow::Int32Array*>(a)->Value(i) != 0;
      case arrow::Type::INT16:  return static_cast<const arrow::Int16Array*>(a)->Value(i) != 0;
      case arrow::Type::INT8:   return static_cast<const arrow::Int8Array*>(a)->Value(i) != 0;
      case arrow::Type::UINT64: return static_cast<const arrow::UInt64Array*>(a)->Value(i) != 0;
      case arrow::Type::UINT32: return static_cast<const arrow::UInt32Array*>(a)->Value(i) != 0;
      case arrow::Type::UINT16: return static_cast<const arrow::UInt16Array*>(a)->Value(i) != 0;
      case arrow::Type::UINT8:  return static_cast<const arrow::UInt8Array*>(a)->Value(i) != 0;
      default: badType(name_, *a, "boolean or integral");
    }
  }

  std::string_view ChunkedColumn::getString(std::int64_t row, std::string_view if_null) const
  {
    if (column_ == nullptr) { return if_null; }
    std::int64_t i = 0;
    const arrow::Array* a = resolve(row, i);
    if (a->IsNull(i)) { return if_null; }
    a = undictionary(a, i);
    switch (a->type_id())
    {
      case arrow::Type::STRING:
        return static_cast<const arrow::StringArray*>(a)->GetView(i);
      case arrow::Type::LARGE_STRING:
        return static_cast<const arrow::LargeStringArray*>(a)->GetView(i);
      case arrow::Type::STRING_VIEW:
        return static_cast<const arrow::StringViewArray*>(a)->GetView(i);
      case arrow::Type::BINARY:
        return static_cast<const arrow::BinaryArray*>(a)->GetView(i);
      case arrow::Type::LARGE_BINARY:
        return static_cast<const arrow::LargeBinaryArray*>(a)->GetView(i);
      default: badType(name_, *a, "a string");
    }
  }

  ChunkedColumn column(const std::shared_ptr<arrow::Table>& table, const std::string& name)
  {
    const int i = table->schema()->GetFieldIndex(name);
    if (i < 0) { return ChunkedColumn{}; }
    ChunkedColumn c(table->column(i));
    c.setName(name);
    return c;
  }

  ChunkedColumn requiredColumn(const std::shared_ptr<arrow::Table>& table,
                               const std::string& name, const std::string& context)
  {
    ChunkedColumn c = column(table, name);
    if (!c.valid())
    {
      std::string have;
      for (const auto& f : table->schema()->fields())
      {
        have += (have.empty() ? "" : ", ") + f->name();
      }
      throw std::runtime_error(context + " has no column '" + name + "'; it has: " + have);
    }
    return c;
  }

} // namespace ODIA
