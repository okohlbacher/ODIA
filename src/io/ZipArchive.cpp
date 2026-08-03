// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/ZipArchive.h>

#include <arrow/buffer.h>
#include <arrow/result.h>
#include <arrow/status.h>

#include <zip.h>

#include <cstring>
#include <mutex>
#include <stdexcept>

namespace ODIA
{

  namespace
  {
    [[noreturn]] void fail(const std::string& what, const std::string& path)
    {
      throw std::runtime_error("ZipArchive: " + what + " (" + path + ")");
    }

    zip_int64_t locate(::zip* za, const std::string& entry)
    {
      return zip_name_locate(za, entry.c_str(), 0);
    }

    /// One entry of a ZIP archive, presented to Arrow as a file.
    ///
    /// It carries its own zip handle rather than borrowing the archive's. Two
    /// entries of a bundle are commonly open at once (precursors and
    /// transitions), libzip's per-archive state is not designed to be shared
    /// across concurrent readers, and Arrow may hand the file to its own IO
    /// threads. A second open costs one central-directory parse.
    class ZipEntryFile : public arrow::io::RandomAccessFile
    {
    public:
      ZipEntryFile(const std::string& path, const std::string& entry)
        : path_(path), entry_(entry)
      {
        int err = 0;
        za_ = zip_open(path_.c_str(), ZIP_RDONLY, &err);
        if (za_ == nullptr) { fail("cannot reopen archive", path_); }

        const auto index = locate(za_, entry_);
        if (index < 0) { zip_close(za_); fail("no entry " + entry_, path_); }

        zip_stat_t st;
        zip_stat_init(&st);
        if (zip_stat_index(za_, index, 0, &st) != 0 || (st.valid & ZIP_STAT_SIZE) == 0)
        {
          zip_close(za_);
          fail("entry " + entry_ + " has no recorded size", path_);
        }
        size_ = static_cast<std::int64_t>(st.size);

        zf_ = zip_fopen_index(za_, index, 0);
        if (zf_ == nullptr) { zip_close(za_); fail("cannot open entry " + entry_, path_); }

        // Parquet is not readable as a stream: the footer lives at the end and
        // must be read before anything else. Entries in a bundle are STORED
        // precisely so they stay seekable, but a deflated bundle from some other
        // writer would land here, and failing now beats a truncated read.
        if (zip_file_is_seekable(zf_) != 1)
        {
          zip_fclose(zf_);
          zip_close(za_);
          fail("entry " + entry_ + " is not seekable; Parquet needs random access "
               "(is it DEFLATEd? bundle entries must be STORED)", path_);
        }
      }

      ~ZipEntryFile() override
      {
        if (zf_ != nullptr) { zip_fclose(zf_); }
        if (za_ != nullptr) { zip_close(za_); }
      }

      arrow::Status Close() override
      {
        std::lock_guard<std::mutex> lock(mutex_);
        if (zf_ != nullptr) { zip_fclose(zf_); zf_ = nullptr; }
        if (za_ != nullptr) { zip_close(za_); za_ = nullptr; }
        closed_ = true;
        return arrow::Status::OK();
      }

      bool closed() const override { return closed_; }

      arrow::Result<std::int64_t> GetSize() override { return size_; }

      arrow::Result<std::int64_t> Tell() const override
      {
        std::lock_guard<std::mutex> lock(mutex_);
        return position_;
      }

      arrow::Status Seek(std::int64_t position) override
      {
        std::lock_guard<std::mutex> lock(mutex_);
        if (position < 0 || position > size_)
        {
          return arrow::Status::IOError("ZipArchive: seek out of range in ", entry_);
        }
        position_ = position;
        return arrow::Status::OK();
      }

      arrow::Result<std::int64_t> Read(std::int64_t nbytes, void* out) override
      {
        std::lock_guard<std::mutex> lock(mutex_);
        ARROW_ASSIGN_OR_RAISE(const auto read, readAtLocked(position_, nbytes, out));
        position_ += read;
        return read;
      }

      arrow::Result<std::shared_ptr<arrow::Buffer>> Read(std::int64_t nbytes) override
      {
        ARROW_ASSIGN_OR_RAISE(auto buffer, arrow::AllocateResizableBuffer(nbytes));
        ARROW_ASSIGN_OR_RAISE(const auto read, Read(nbytes, buffer->mutable_data()));
        ARROW_RETURN_NOT_OK(buffer->Resize(read, /*shrink_to_fit=*/true));
        return std::shared_ptr<arrow::Buffer>(std::move(buffer));
      }

      arrow::Result<std::int64_t> ReadAt(std::int64_t position, std::int64_t nbytes,
                                         void* out) override
      {
        std::lock_guard<std::mutex> lock(mutex_);
        return readAtLocked(position, nbytes, out);
      }

      arrow::Result<std::shared_ptr<arrow::Buffer>> ReadAt(std::int64_t position,
                                                           std::int64_t nbytes) override
      {
        ARROW_ASSIGN_OR_RAISE(auto buffer, arrow::AllocateResizableBuffer(nbytes));
        ARROW_ASSIGN_OR_RAISE(const auto read, ReadAt(position, nbytes, buffer->mutable_data()));
        ARROW_RETURN_NOT_OK(buffer->Resize(read, /*shrink_to_fit=*/true));
        return std::shared_ptr<arrow::Buffer>(std::move(buffer));
      }

    private:
      arrow::Result<std::int64_t> readAtLocked(std::int64_t position, std::int64_t nbytes,
                                               void* out)
      {
        if (closed_) { return arrow::Status::IOError("ZipArchive: read after close"); }
        if (position < 0 || nbytes < 0)
        {
          return arrow::Status::IOError("ZipArchive: negative read in ", entry_);
        }
        if (position >= size_) { return std::int64_t{0}; }
        const std::int64_t want = std::min(nbytes, size_ - position);

        if (zip_fseek(zf_, position, SEEK_SET) != 0)
        {
          return arrow::Status::IOError("ZipArchive: seek failed in ", entry_, ": ",
                                        zip_file_strerror(zf_));
        }

        // zip_fread is not obliged to return everything asked for in one call.
        // A short read treated as EOF would truncate a column chunk, and Parquet
        // would report a corrupt file rather than an IO problem.
        std::int64_t done = 0;
        auto* dst = static_cast<std::uint8_t*>(out);
        while (done < want)
        {
          const auto n = zip_fread(zf_, dst + done, static_cast<zip_uint64_t>(want - done));
          if (n < 0)
          {
            return arrow::Status::IOError("ZipArchive: read failed in ", entry_, ": ",
                                          zip_file_strerror(zf_));
          }
          if (n == 0) { break; }
          done += n;
        }
        if (done != want)
        {
          return arrow::Status::IOError("ZipArchive: entry ", entry_, " ended after ", done,
                                        " of ", want, " bytes at offset ", position,
                                        "; the archive is truncated");
        }
        return done;
      }

      std::string path_;
      std::string entry_;
      ::zip* za_ = nullptr;
      zip_file_t* zf_ = nullptr;
      std::int64_t size_ = 0;
      std::int64_t position_ = 0;
      bool closed_ = false;
      mutable std::mutex mutex_;
    };
  } // namespace

  ZipArchive::ZipArchive(const std::string& path) : path_(path)
  {
    int err = 0;
    za_ = zip_open(path_.c_str(), ZIP_RDONLY, &err);
    if (za_ == nullptr)
    {
      zip_error_t error;
      zip_error_init_with_code(&error, err);
      const std::string what = zip_error_strerror(&error);
      zip_error_fini(&error);
      fail("cannot open archive: " + what, path_);
    }
  }

  ZipArchive::~ZipArchive()
  {
    if (za_ != nullptr) { zip_close(za_); }
  }

  bool ZipArchive::has(const std::string& entry) const
  {
    return locate(za_, entry) >= 0;
  }

  std::vector<std::string> ZipArchive::entries() const
  {
    std::vector<std::string> names;
    const auto n = zip_get_num_entries(za_, 0);
    names.reserve(static_cast<std::size_t>(n < 0 ? 0 : n));
    for (zip_int64_t i = 0; i < n; ++i)
    {
      const char* name = zip_get_name(za_, i, 0);
      if (name != nullptr) { names.emplace_back(name); }
    }
    return names;
  }

  std::string ZipArchive::read(const std::string& entry) const
  {
    const auto index = locate(za_, entry);
    if (index < 0) { fail("no entry " + entry, path_); }

    zip_stat_t st;
    zip_stat_init(&st);
    if (zip_stat_index(za_, index, 0, &st) != 0 || (st.valid & ZIP_STAT_SIZE) == 0)
    {
      fail("entry " + entry + " has no recorded size", path_);
    }

    zip_file_t* zf = zip_fopen_index(za_, index, 0);
    if (zf == nullptr) { fail("cannot open entry " + entry, path_); }

    std::string out(static_cast<std::size_t>(st.size), '\0');
    std::size_t done = 0;
    while (done < out.size())
    {
      const auto n = zip_fread(zf, out.data() + done, out.size() - done);
      if (n <= 0) { break; }
      done += static_cast<std::size_t>(n);
    }
    zip_fclose(zf);
    if (done != out.size())
    {
      fail("entry " + entry + " is truncated", path_);
    }
    return out;
  }

  std::shared_ptr<arrow::io::RandomAccessFile> ZipArchive::open(const std::string& entry) const
  {
    if (!has(entry)) { fail("no entry " + entry, path_); }
    return std::make_shared<ZipEntryFile>(path_, entry);
  }

} // namespace ODIA
