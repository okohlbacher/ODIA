// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <arrow/io/interfaces.h>

#include <memory>
#include <string>
#include <vector>

struct zip;

namespace ODIA
{

  /// Read-only access to a ZIP archive, with entries exposed to Arrow.
  ///
  /// `.oswpq` is an archive, not a directory. Treating it as one is the
  /// documented commonest mistake with the format and it fails late -- the
  /// predecessor project lost a 47-minute run to it, because the error only
  /// surfaces when something tries to open `<bundle>/library/x.parquet` as a
  /// filesystem path.
  class ZipArchive
  {
  public:
    explicit ZipArchive(const std::string& path);
    ~ZipArchive();

    ZipArchive(const ZipArchive&) = delete;
    ZipArchive& operator=(const ZipArchive&) = delete;

    bool has(const std::string& entry) const;

    /// How many entries carry this exact name. ZIP permits repeats and
    /// readers disagree about which wins, so a consumer that cares must ask.
    std::size_t count(const std::string& entry) const;

    /// Every entry name, in archive order, including directory entries.
    std::vector<std::string> entries() const;

    /// The whole entry as bytes. For metadata, not for Parquet.
    ///
    /// @param max_bytes refuse an entry declaring more than this. The size
    ///        comes from the archive's central directory, which is the file's
    ///        own unverified claim about itself: a 4 MB bundle declaring a
    ///        4 GB entry allocated 4.25 GB here and then loaded normally.
    std::string read(const std::string& entry,
                     std::size_t max_bytes = 64u << 20) const;

    /// A seekable view of one entry, for Parquet, which reads its footer first
    /// and then seeks to the column chunks it wants. Nothing is materialised:
    /// the entries in a bundle run to gigabytes and only some columns are read.
    std::shared_ptr<arrow::io::RandomAccessFile> open(const std::string& entry) const;

    const std::string& path() const { return path_; }

  private:
    std::string path_;
    ::zip* za_ = nullptr;
  };

} // namespace ODIA
