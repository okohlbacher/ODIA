// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <mzpeak.h>

#include <algorithm>
#include <cstdint>
#include <map>
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
    setValidFormats_("in", {"mzpeak", "mzML"}, false);

    registerOutputFile_("out", "<file>", "", "Optional text report; written to stdout if omitted.", false);
    setValidFormats_("out", {"txt"}, false);

    registerFlag_("windows", "List the isolation-window scheme (mzML input only).");
  }

  ExitCodes main_(int, const char**) override
  {
    const String in = getStringOption_("in");
    const String out = getStringOption_("out");
    const bool list_windows = getFlag_("windows");

    std::ostringstream report;

    if (in.hasSuffix(".mzpeak"))
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
  /// Report what the mzpeak library exposes about an mzPeak archive.
  ExitCodes describeMzPeak_(const String& in, std::ostringstream& report)
  {
    try
    {
      MzPeak::Index index = MzPeak::open(in.c_str());
      MzPeak::Spectra spectra = index.spectra();

      report << "file:          " << in << "\n";
      report << "format:        mzPeak\n";
      report << "archive members:\n";
      for (const auto& f : index.files()) { report << "  - " << f.file_name << "\n"; }

      // One pass over the spectra: MS-level histogram and peak counts. This is
      // also the honest way to time full sequential access to an archive.
      std::map<int, Size> level_counts;
      std::map<int, UInt64> level_peaks;
      startProgress(0, spectra.size(), "reading spectra");
      Size i = 0;
      for (const auto& s : spectra)
      {
        setProgress(i++);
        const int level = static_cast<int>(s.ms_level());
        ++level_counts[level];
        level_peaks[level] += s.mz().size();
      }
      endProgress();

      report << "spectra:       " << spectra.size() << "\n";
      report << "MS levels:\n";
      for (const auto& [level, count] : level_counts)
      {
        report << "  MS" << level << ": " << count << " spectra, "
               << level_peaks[level] << " peaks\n";
      }
      report << "\nnote: retention time, precursor isolation windows and ion\n"
                "mobility are not exposed by the mzpeak high-level API; they\n"
                "live in the archive's Parquet metadata tables.\n";
    }
    catch (const std::exception& e)
    {
      writeLogError_(String("Failed to read mzPeak file: ") + e.what());
      return INPUT_FILE_CORRUPT;
    }
    return EXECUTION_OK;
  }

  /// Report the DIA layout of an mzML run using OpenMS.
  ExitCodes describeMzML_(const String& in, bool list_windows, std::ostringstream& report)
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
