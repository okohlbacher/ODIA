// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <mzpeak.h>
#include <mzpeak/util/parquet.h>

#include <arrow/api.h>
#include <arrow/compute/api.h>
#include <parquet/arrow/reader.h>

#include <algorithm>
#include <ranges>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace OpenMS;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page ODIA_ODIAInfo ODIAInfo

@brief Summarises the acquisition layout of a DIA run.

Reads a DIA run in either mzPeak or mzML format and reports the structure a
DIA extraction needs to know about: how many spectra there are per MS level,
the retention-time range, the isolation-window scheme, and whether an ion
mobility dimension is present.

mzPeak files are read through the OpenMS/mzpeak library; mzML files are read
through OpenMS. The two paths report different levels of detail, because the
mzpeak high-level API currently exposes only m/z, intensity and MS level per
spectrum -- retention time, precursor isolation windows and ion mobility live
in the archive's Parquet metadata tables and are reachable only through the
low-level interface.

<B>The command line parameters of this tool are:</B>
@verbinclude ODIA_ODIAInfo.cli
*/

/// @cond TOPPCLASSES

namespace
{
  /// One row of the isolation-window scheme, as counted over a run.
  struct WindowStat
  {
    double lower = 0.0;
    double upper = 0.0;
    Size count = 0;
  };

  /// Quantise a window edge so that near-identical windows collapse into one row.
  Int64 quantise(double mz)
  {
    return static_cast<Int64>(std::llround(mz * 1000.0));
  }
} // namespace

class TOPPODIAInfo : public TOPPBase
{
public:
  TOPPODIAInfo() :
    // official=false: ODIA lives outside the OpenMS source tree, so this tool
    // is deliberately not in OpenMS's ToolHandler list.
    TOPPBase("ODIAInfo", "Summarises the acquisition layout of a DIA run.", false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input DIA run (mzPeak or mzML).");
    // Deliberately no setValidFormats_ here. OpenMS 3.6 has no mzPeak entry in
    // FileTypes, so declaring formats makes TOPPBase sniff the file, see a
    // multi-entry zip archive and reject it as a malformed compressed mzML
    // before main_ ever runs. Extension dispatch below handles it instead.

    registerOutputFile_("out", "<file>", "", "Optional text report; written to stdout if omitted.", false);
    setValidFormats_("out", {"txt"}, false);

    registerFlag_("windows", "List the isolation-window scheme (mzML input only).");
  }

  ExitCodes main_(int, const char**) override
  {
    const std::string in = getStringOption_("in");
    const std::string out = getStringOption_("out");
    const bool list_windows = getFlag_("windows");

    std::ostringstream report;

    if (in.ends_with(".mzpeak"))
    {
      if (ExitCodes rc = describeMzPeak_(in, report); rc != EXECUTION_OK) { return rc; }
    }
    else
    {
      if (ExitCodes rc = describeMzML_(in, list_windows, report); rc != EXECUTION_OK) { return rc; }
    }

    if (out.empty())
    {
      std::cout << report.str();
    }
    else
    {
      std::ofstream os(out.c_str());
      if (!os) { writeLogError_("Cannot write to '" + out + "'."); return CANNOT_WRITE_OUTPUT_FILE; }
      os << report.str();
    }
    return EXECUTION_OK;
  }

private:
  /// Read one column of a table as doubles, whatever numeric type it is stored as.
  static std::vector<double> readNumericColumn_(const std::shared_ptr<arrow::Table>& table,
                                                const std::string& name)
  {
    std::vector<double> values;
    const int idx = table->schema()->GetFieldIndex(name);
    if (idx < 0) { return values; }

    auto casted = arrow::compute::Cast(arrow::Datum(table->column(idx)),
                                       arrow::float64());
    if (!casted.ok()) { return values; }

    const auto chunked = casted->chunked_array();
    values.reserve(static_cast<size_t>(chunked->length()));
    for (const auto& chunk : chunked->chunks())
    {
      const auto& a = static_cast<const arrow::DoubleArray&>(*chunk);
      for (int64_t i = 0; i < a.length(); ++i)
      {
        if (!a.IsNull(i)) { values.push_back(a.Value(i)); }
      }
    }
    return values;
  }

  /// Report the layout of an mzPeak archive.
  ///
  /// Deliberately goes through the low-level Parquet interface rather than
  /// Index::spectra(). Two reasons: the high-level Spectrum type exposes only
  /// m/z, intensity and MS level, so retention time would be unreachable; and
  /// on archives written by mzpeak-convert >= 0.7.0 the high-level path throws
  /// outright, because the C++ library still expects the pre-0.7.0 layout that
  /// nested the per-spectrum metadata under a "spectrum" group, while 0.7.0
  /// writes those fields flat. Reading the metadata table directly is both
  /// version-tolerant and the access pattern ODIA needs anyway.
  ExitCodes describeMzPeak_(const std::string& in, std::ostringstream& report)
  {
    try
    {
      MzPeak::Index index = MzPeak::open(in.c_str());

      report << "file:          " << in << "\n";
      report << "format:        mzPeak\n";
      report << "archive members:\n";
      for (const auto& f : index.files()) { report << "  - " << f.file_name << "\n"; }

      constexpr const char* kMetadata = "spectra_metadata.parquet";
      auto it = std::ranges::find(index.files(), kMetadata,
                                  &MzPeak::Schema::File::file_name);
      if (it == index.files().end())
      {
        writeLogError_(std::string("archive has no ") + kMetadata);
        return INPUT_FILE_CORRUPT;
      }

      auto parquet = index.parquet(*it);
      std::shared_ptr<arrow::Table> table;
      if (auto st = parquet->reader().ReadTable(&table); !st.ok())
      {
        writeLogError_("Cannot read " + std::string(kMetadata) + ": " + st.ToString());
        return INPUT_FILE_CORRUPT;
      }

      report << "spectra:       " << table->num_rows() << "\n";

      const auto levels = readNumericColumn_(table, "ms_level");
      const auto times = readNumericColumn_(table, "time");

      if (!levels.empty())
      {
        std::map<int, Size> level_counts;
        for (double l : levels) { ++level_counts[static_cast<int>(l)]; }
        report << "MS levels:\n";
        for (const auto& [level, count] : level_counts)
        {
          report << "  MS" << level << ": " << count << " spectra\n";
        }
      }

      if (!times.empty())
      {
        // mzPeak stores the scan start time in minutes, unlike mzML's seconds.
        const auto [lo, hi] = std::ranges::minmax(times);
        report << "RT range:      " << lo << " .. " << hi << " min\n";
      }

      report << "\nnote: read through the low-level Parquet interface; the\n"
                "mzpeak high-level API exposes only m/z, intensity and MS\n"
                "level per spectrum.\n";
    }
    catch (const std::exception& e)
    {
      writeLogError_(std::string("Failed to read mzPeak file: ") + e.what());
      return INPUT_FILE_CORRUPT;
    }
    return EXECUTION_OK;
  }

  /// Report the DIA layout of an mzML run using OpenMS.
  ExitCodes describeMzML_(const std::string& in, bool list_windows, std::ostringstream& report)
  {
    PeakMap exp;
    FileHandler().loadExperiment(in, exp, {FileTypes::MZML}, log_type_);

    std::map<int, Size> level_counts;
    std::map<Int64, WindowStat> windows;
    double rt_min = std::numeric_limits<double>::max();
    double rt_max = std::numeric_limits<double>::lowest();
    Size im_spectra = 0;
    UInt64 peaks = 0;

    for (const auto& s : exp)
    {
      ++level_counts[static_cast<int>(s.getMSLevel())];
      peaks += s.size();
      rt_min = std::min(rt_min, s.getRT());
      rt_max = std::max(rt_max, s.getRT());

      if (s.containsIMData() || s.getDriftTime() >= 0.0) { ++im_spectra; }

      for (const auto& prec : s.getPrecursors())
      {
        const double lower = prec.getMZ() - prec.getIsolationWindowLowerOffset();
        const double upper = prec.getMZ() + prec.getIsolationWindowUpperOffset();
        // Key on the lower edge; DIA schemes repeat the same windows each cycle.
        auto& w = windows[quantise(lower)];
        w.lower = lower;
        w.upper = upper;
        ++w.count;
      }
    }

    report << "file:          " << in << "\n";
    report << "format:        mzML\n";
    report << "spectra:       " << exp.size() << "\n";
    report << "peaks:         " << peaks << "\n";
    if (!exp.empty())
    {
      report << "RT range:      " << rt_min << " .. " << rt_max << " s\n";
    }
    report << "MS levels:\n";
    for (const auto& [level, count] : level_counts)
    {
      report << "  MS" << level << ": " << count << " spectra\n";
    }
    report << "IM spectra:    " << im_spectra
           << (im_spectra ? " (ion mobility present)" : " (no ion mobility)") << "\n";
    report << "distinct isolation windows: " << windows.size() << "\n";

    // A DIA run repeats its window scheme once per cycle, so the number of
    // cycles is the number of times each window was acquired.
    if (!windows.empty())
    {
      Size max_count = 0;
      for (const auto& [key, w] : windows) { max_count = std::max(max_count, w.count); }
      report << "cycles (max window repeats): " << max_count << "\n";
    }

    if (list_windows && !windows.empty())
    {
      report << "\nisolation windows:\n";
      report << "  lower_mz    upper_mz    width    count\n";
      for (const auto& [key, w] : windows)
      {
        report << "  " << std::fixed << std::setprecision(4)
               << std::setw(10) << w.lower << "  "
               << std::setw(10) << w.upper << "  "
               << std::setw(7) << (w.upper - w.lower) << "  "
               << std::setw(7) << w.count << "\n";
      }
    }
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPODIAInfo tool;
  return tool.main(argc, argv);
}

/// @endcond
