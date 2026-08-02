// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <mzpeak.h>


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

mzPeak files are read through okohlbacher/mzpeak-openms, mzML files through
OpenMS. Both paths report the same quantities on purpose: running the tool over
the two encodings of one run is the cheapest cross-check that the two readers
agree.

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

    registerFlag_("windows", "List the isolation-window scheme.");
    registerFlag_("peaks", "Also decode and count peaks (mzPeak: slow, see below).");
  }

  ExitCodes main_(int, const char**) override
  {
    const std::string in = getStringOption_("in");
    const std::string out = getStringOption_("out");
    const bool list_windows = getFlag_("windows");
    const bool count_peaks = getFlag_("peaks");

    std::ostringstream report;

    if (in.ends_with(".mzpeak"))
    {
      if (ExitCodes rc = describeMzPeak_(in, list_windows, count_peaks, report); rc != EXECUTION_OK) { return rc; }
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
  /// Report the DIA layout of an mzPeak archive.
  ///
  /// Reports the same quantities as the mzML path, deliberately: running both
  /// on the two encodings of the same run is the cheapest cross-check we have
  /// that the mzPeak reader agrees with OpenMS.
  ExitCodes describeMzPeak_(const std::string& in, bool list_windows,
                            bool count_peaks, std::ostringstream& report)
  {
    try
    {
      MzPeak::Index index = MzPeak::open(in.c_str());
      MzPeak::Spectra spectra = index.spectra();

      report << "file:          " << in << "\n";
      report << "format:        mzPeak\n";
      report << "archive members:\n";
      for (const auto& f : index.files()) { report << "  - " << f.file_name << "\n"; }

      std::map<int, Size> level_counts;
      std::map<Int64, WindowStat> windows;
      double rt_min = std::numeric_limits<double>::max();
      double rt_max = std::numeric_limits<double>::lowest();
      Size im_spectra = 0;
      UInt64 peaks = 0;
      Size undecodable = 0;

      ProgressLogger progress;
      progress.setLogType(log_type_);
      progress.startProgress(0, spectra.size(), "reading spectra");
      Size i = 0;
      for (const auto& s : spectra)
      {
        progress.setProgress(i++);
        ++level_counts[static_cast<int>(s.ms_level())];

        // Peak decoding is off by default, and deliberately so. It is not
        // needed for an acquisition summary, and in the current reader it costs
        // ~277 ms per spectrum -- about a thousand times slower than OpenMS
        // parses the same run out of mzML -- so counting peaks over 12_80 takes
        // an hour. It also cannot decode every encoding in these archives:
        // "chunked array decoding is not implemented (MS:1000515)" on the first
        // MS2 spectrum. Both are recorded in doc/01-constraints.md.
        if (count_peaks)
        {
          try { peaks += s.mz().size(); }
          catch (const std::exception&) { ++undecodable; }
        }

        // Retention time is in seconds here, matching mzML; the underlying
        // table stores minutes and the reader converts.
        if (const auto rt = s.retention_time())
        {
          rt_min = std::min(rt_min, *rt);
          rt_max = std::max(rt_max, *rt);
        }
        if (s.ion_mobility()) { ++im_spectra; }

        for (const auto& prec : s.precursors())
        {
          const auto& w = prec.isolation_window;
          if (!w.target_mz) { continue; }
          const double lower = *w.target_mz - (w.lower_offset ? *w.lower_offset : 0.0f);
          const double upper = *w.target_mz + (w.upper_offset ? *w.upper_offset : 0.0f);
          auto& stat = windows[quantise(lower)];
          stat.lower = lower;
          stat.upper = upper;
          ++stat.count;
        }
      }
      progress.endProgress();

      report << "spectra:       " << spectra.size() << "\n";
      report << "peaks:         " << (count_peaks ? std::to_string(peaks) : std::string("(not decoded; -peaks to enable)"));
      if (undecodable)
      {
        report << "  (" << undecodable << " spectra could not be decoded)";
      }
      report << "\n";
      if (rt_max >= rt_min)
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
      writeWindows_(windows, list_windows, report);
    }
    catch (const std::exception& e)
    {
      writeLogError_(std::string("Failed to read mzPeak file: ") + e.what());
      return INPUT_FILE_CORRUPT;
    }
    return EXECUTION_OK;
  }

  /// Shared tail of both reports: cycle count and, optionally, the window list.
  static void writeWindows_(const std::map<Int64, WindowStat>& windows,
                            bool list_windows, std::ostringstream& report)
  {
    if (windows.empty()) { return; }

    // A DIA run repeats its window scheme once per cycle, so the number of
    // cycles is the number of times each window was acquired.
    Size max_count = 0;
    for (const auto& [key, w] : windows) { max_count = std::max(max_count, w.count); }
    report << "cycles (max window repeats): " << max_count << "\n";

    if (!list_windows) { return; }
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
    writeWindows_(windows, list_windows, report);

    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPODIAInfo tool;
  return tool.main(argc, argv);
}

/// @endcond
