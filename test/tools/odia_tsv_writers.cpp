// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// The two TSV writers must be BYTE-IDENTICAL to the std::ofstream form they
// replace.
//
// Both were formatting-bound rather than disk-bound -- `ofstream <<` measured
// 25.4 MB/s against 1,042 MB/s for fwrite of pre-rendered bytes on the same
// NVMe -- so they now render through std::to_chars into a reusable buffer. That
// is a ~40x headroom bought against exactly one risk: to_chars' default is
// shortest-round-trip and the stream's is 6 significant digits, so getting the
// precision wrong does not round differently, it writes a different file. A
// library whose m/z quietly gained eleven digits, or whose intensities quietly
// lost three, would pass every other test in the suite.
//
// Two cases, because they fail differently:
//
//   files       the real writers against an ofstream reference, over a
//               deliberately awkward fixture, compared byte for byte. This is
//               the one that would catch a field reordered, a tab dropped, or
//               a precision changed at one call site.
//   formatting  TextWriter::number against the stream directly, over ~4.8 M
//               values including random bit patterns, denormals, infinities,
//               NaNs and realistic m/z / intensity / RT magnitudes, at every
//               precision the writers use. This is the one that would catch
//               to_chars and the stream disagreeing on a value the fixture
//               does not happen to contain.
//
// The reference below is a transcription of the ofstream code these writers
// replace. It is an oracle rather than a copy: it goes through num_put and the
// locale, which is the thing being replaced.
//
// Usage: odia_tsv_writers <files|formatting> [scratch-directory]

#include <odia/ChromatogramExtractor.h>
#include <odia/ChromatogramTsv.h>
#include <odia/DIANNLibraryFile.h>
#include <odia/FragvecTsv.h>
#include <odia/Library.h>
#include <odia/PeakGroupScorer.h>
#include <odia/TextWriter.h>

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace
{
  const float NA = std::numeric_limits<float>::quiet_NaN();

  int failures = 0;

  void check(bool ok, const std::string& what)
  {
    std::printf("%s %s\n", ok ? "  ok  " : "  FAIL", what.c_str());
    if (!ok) { ++failures; }
  }

  // ------------------------------------------------------- the ofstream oracle

  /// DIANNLibraryFile::storeTSV as it was written with std::ofstream.
  void referenceLibraryTsv(const std::string& filename, const ODIA::Library& library)
  {
    using Columns = ODIA::DIANNLibraryFile::Columns;
    std::ofstream out(filename);
    if (!out) { throw std::runtime_error("cannot write library: " + filename); }

    out << Columns::PRECURSOR_ID << '\t' << Columns::MODIFIED_SEQUENCE << '\t'
        << Columns::PRECURSOR_CHARGE << '\t' << Columns::DECOY << '\t'
        << Columns::RT << '\t' << Columns::IM << '\t'
        << Columns::PRECURSOR_MZ << '\t' << Columns::PRODUCT_MZ << '\t'
        << Columns::RELATIVE_INTENSITY << '\t' << Columns::FRAGMENT_TYPE << '\t'
        << Columns::FRAGMENT_CHARGE << '\t' << Columns::FRAGMENT_SERIES_NUMBER << '\t'
        << Columns::FRAGMENT_LOSS_TYPE << '\t' << Columns::PROTEIN_GROUP << '\t'
        << Columns::CCS << '\n';

    out << std::defaultfloat;
    const auto& p = library.precursors();
    const auto& t = library.transitions();
    for (std::size_t i = 0; i < library.precursorCount(); ++i)
    {
      const auto seq = library.strings().get(p.modified_sequence[i]);
      const auto pg = library.strings().get(p.protein_group[i]);
      const int z = p.charge[i];
      const bool has_im = !std::isnan(p.im[i]);
      const bool has_ccs = i < p.ccs.size() && !std::isnan(p.ccs[i]);
      const std::uint32_t begin = p.transition_begin[i];
      for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
      {
        const std::uint32_t j = begin + k;
        out << seq << z << (p.decoy[i] ? "_decoy" : "") << '\t'
            << seq << '\t' << z << '\t'
            << static_cast<int>(p.decoy[i]) << '\t';
        out.precision(9);
        if (std::isnan(p.irt[i])) { out << ""; } else { out << p.irt[i]; }
        out << '\t';
        if (has_im) { out << p.im[i]; } else { out << 0; }
        out << '\t';
        out.precision(10);
        out << ODIA::fromFixed(p.mz[i]) << '\t' << ODIA::fromFixed(t.product_mz[j]) << '\t';
        out.precision(9);
        out << t.library_intensity[j] << '\t' << toString(t.type[j]) << '\t'
            << static_cast<int>(t.charge[j]) << '\t'
            << static_cast<int>(t.ordinal[j]) << '\t'
            << toString(t.loss[j]) << '\t' << pg << '\t';
        if (has_ccs) { out << p.ccs[i]; }
        out << '\n';
      }
    }
    out.flush();
    if (!out) { throw std::runtime_error("failed while writing library: " + filename); }
  }

  /// writeChromatogramTsv as it was written with std::ofstream, inside the tool.
  void referenceChromatogramTsv(const std::string& path, const ODIA::Library& library,
                                const ODIA::Chromatograms& chromatograms)
  {
    std::ofstream out(path);
    if (!out) { throw std::runtime_error("cannot open " + path); }
    out << "Precursor.Id\tDecoy\tTransition.Index\tProduct.Mz\tRT\tIntensity\n";

    const auto& p = library.precursors();
    const auto& t = library.transitions();
    for (std::size_t i = 0; i < library.precursorCount(); ++i)
    {
      const auto seq = library.strings().get(p.modified_sequence[i]);
      const std::string id = std::string(seq) + std::to_string(static_cast<int>(p.charge[i]));
      for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
      {
        const std::uint32_t tr = p.transition_begin[i] + k;
        if (tr >= chromatograms.begin.size()) { continue; }
        const std::uint64_t b = chromatograms.begin[tr];
        const std::uint32_t n = chromatograms.count[tr];
        for (std::uint32_t j = 0; j < n; ++j)
        {
          out << id << '\t' << (p.decoy[i] ? '1' : '0') << '\t'
              << tr << '\t' << ODIA::fromFixed(t.product_mz[tr]) << '\t'
              << chromatograms.retentionTime(tr, j) << '\t'
              << chromatograms.intensity[b + j] << '\n';
        }
      }
    }
    if (!out) { throw std::runtime_error("write failed for " + path); }
  }

  // ------------------------------------------------------------------ fixture

  /// A library whose numbers are chosen to be awkward for `%g`, not realistic.
  ///
  /// Values that sit on a rounding boundary, that cross the fixed/exponent
  /// switch in both directions, that are absent (NaN -> an empty field), and a
  /// protein group long enough to exercise the buffered writer's overflow path.
  ODIA::Library awkwardLibrary()
  {
    ODIA::Library lib;
    auto& p = lib.precursors();
    auto& t = lib.transitions();

    const std::vector<double> mz = {
      0.0, 1.0, 100.0, 999.9995, 1000.00005, 0.1, 1.0 / 3.0, 2.0 / 3.0,
      1e-5, 1e-4, 123456.789012, 1999.9999999, 500.5, 0.5, 1234.5678901234
    };
    const std::vector<float> small = {
      0.0f, 1.0f, 0.5f, 1.0f / 3.0f, 1e-5f, 1e-4f, 9999995.0f, 0.000999999f,
      123456.789f, -0.0f, -1.5f, 1e20f, 1e-20f, 999999.5f, 100.0f
    };

    const std::string long_group(4096, 'P');
    for (std::size_t i = 0; i < mz.size(); ++i)
    {
      p.mz.push_back(ODIA::toFixed(mz[i]));
      // Every "absent" field exercised: no RT, no 1/K0, no CCS.
      p.irt.push_back(i % 4 == 0 ? NA : small[i]);
      p.im.push_back(i % 3 == 0 ? NA : small[(i + 1) % small.size()]);
      p.ccs.push_back(i % 5 == 0 ? NA : small[(i + 2) % small.size()]);
      p.charge.push_back(static_cast<std::uint8_t>(1 + i % 4));
      p.decoy.push_back(static_cast<std::uint8_t>(i % 2));
      p.modified_sequence.push_back(lib.strings().intern("PEPT(Phospho)IDEK"));
      // One group far longer than a line, so the writer's large-field path runs.
      p.protein_group.push_back(lib.strings().intern(i == 3 ? long_group : "sp|P00001|A_HUMAN"));
      p.transition_begin.push_back(static_cast<std::uint32_t>(t.product_mz.size()));
      p.transition_count.push_back(3);
      for (int k = 0; k < 3; ++k)
      {
        t.product_mz.push_back(ODIA::toFixed(mz[(i + std::size_t(k)) % mz.size()]));
        t.library_intensity.push_back(small[(i + std::size_t(k)) % small.size()]);
        t.type.push_back(k == 0 ? ODIA::FragmentType::B : ODIA::FragmentType::Y);
        t.ordinal.push_back(static_cast<std::uint8_t>(2 + k));
        t.charge.push_back(static_cast<std::uint8_t>(1 + k % 2));
        t.loss.push_back(k == 2 ? ODIA::LossType::Water : ODIA::LossType::None);
      }
    }
    return lib;
  }

  /// Chromatograms over that library, filled with the same awkward values.
  ODIA::Chromatograms awkwardChromatograms(const ODIA::Library& lib)
  {
    ODIA::Chromatograms c;
    const auto& t = lib.transitions();
    const std::size_t n_trans = t.product_mz.size();
    const std::uint32_t cycles = 7;

    c.axes.resize(1);
    for (std::uint32_t j = 0; j < cycles; ++j)
    {
      // Times that are not round: 6 significant digits is where a float's
      // decimal expansion stops being obvious.
      c.axes[0].push_back(100.0f + float(j) * 0.7333333f);
    }

    const std::vector<float> values = {
      0.0f, 1.0f, 1.0f / 3.0f, 1e-5f, 9999995.0f, 1e20f, 123456.789f
    };
    std::uint64_t running = 0;
    for (std::size_t j = 0; j < n_trans; ++j)
    {
      c.begin.push_back(running);
      c.count.push_back(cycles);
      running += cycles;
    }
    // One precursor owning every transition: the axis fields are per precursor
    // now, so a hand-built Chromatograms has to say which precursor that is.
    c.precursor_axis.assign(1, 0);
    c.precursor_axis_begin.assign(1, 0);
    c.precursor_cycles.assign(1, static_cast<std::uint32_t>(cycles));
    c.precursor_transition_begin.assign(1, 0);
    c.intensity.reserve(running);
    for (std::size_t k = 0; k < running; ++k)
    {
      c.intensity.push_back(values[k % values.size()]);
    }
    return c;
  }

  std::string slurp(const std::string& path)
  {
    std::ifstream in(path, std::ios::binary);
    if (!in) { throw std::runtime_error("cannot read " + path); }
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
  }

  /// Where the two differ, if they do -- a byte offset alone says nothing about
  /// which field moved.
  void checkIdentical(const std::string& a_path, const std::string& b_path,
                      const std::string& what)
  {
    const auto a = slurp(a_path), b = slurp(b_path);
    if (a == b)
    {
      check(true, what + " (" + std::to_string(a.size()) + " bytes identical)");
      return;
    }
    std::size_t at = 0;
    while (at < a.size() && at < b.size() && a[at] == b[at]) { ++at; }
    const std::size_t from = at > 60 ? at - 60 : 0;
    std::printf("  first difference at byte %zu of %zu / %zu\n", at, a.size(), b.size());
    std::printf("    to_chars: %s\n", a.substr(from, 160).c_str());
    std::printf("    ofstream: %s\n", b.substr(from, 160).c_str());
    check(false, what);
  }

  // -------------------------------------------------------------------- cases

  void caseFiles(const std::string& dir)
  {
    const auto lib = awkwardLibrary();
    const auto chrom = awkwardChromatograms(lib);

    const std::string lib_new = dir + "/tsv_lib_tochars.tsv";
    const std::string lib_ref = dir + "/tsv_lib_ofstream.tsv";
    ODIA::DIANNLibraryFile::storeTSV(lib_new, lib);
    referenceLibraryTsv(lib_ref, lib);
    checkIdentical(lib_new, lib_ref, "the library TSV is byte-identical to the ofstream form");

    const std::string chrom_new = dir + "/tsv_chrom_tochars.tsv";
    const std::string chrom_ref = dir + "/tsv_chrom_ofstream.tsv";
    ODIA::writeChromatogramTsv(chrom_new, lib, chrom);
    referenceChromatogramTsv(chrom_ref, lib, chrom);
    checkIdentical(chrom_new, chrom_ref,
                   "the chromatogram TSV is byte-identical to the ofstream form");

    // A fixture that is empty, or all-identical, agrees with anything. Say what
    // it actually covered.
    const auto text = slurp(lib_new);
    check(text.size() > 2000, "the library fixture is not trivially small");
    check(text.find("\t\t") != std::string::npos,
          "and it contains the empty fields NaN must produce");
  }

  /// TextWriter::number against the stream, value by value.
  ///
  /// The file comparison above can only cover the values the fixture holds.
  /// This covers the format: every double the writers could ever be handed, at
  /// the three precisions they use.
  void caseFormatting()
  {
    const auto viaStream = [](double v, int prec) {
      std::ostringstream os;
      os << std::defaultfloat;
      os.precision(prec);
      os << v;
      return os.str();
    };
    const auto viaChars = [](double v, int prec) {
      char b[256];
      const auto r = std::to_chars(b, b + sizeof b, v, std::chars_format::general, prec);
      if (r.ec != std::errc{}) { return std::string("<did not fit>"); }
      return std::string(b, r.ptr);
    };

    const std::vector<double> pinned = {
      0.0, -0.0, 1.0, -1.0, 0.1, 1.0 / 3.0, 2.0 / 3.0, 1e-5, 1e-4, 1e21, 1e22,
      123456789.0, 1234567890123.0, 0.5, 1.5, 2.5, 1e308, 5e-324,
      2.2250738585072014e-308, 999999.5, 9999999.5, 0.000999999949, 1e-300,
      -3.14159265358979,
      std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::quiet_NaN(),
      -std::numeric_limits<double>::quiet_NaN()
    };

    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> mz(100.0, 2000.0), intensity(0.0, 1e9),
                                           rt(0.0, 3600.0), irt(-50.0, 200.0);
    long long compared = 0, disagree = 0;
    std::string first;

    // 6 is the stream's default, which the chromatogram writer never changed;
    // 9 and 10 are what the library writer sets.
    for (const int prec : {6, 9, 10})
    {
      std::vector<double> values = pinned;
      // Random BIT PATTERNS, not random numbers: that is what reaches the
      // subnormals and the extreme exponents where the two implementations
      // would most plausibly part company.
      for (int i = 0; i < 400000; ++i)
      {
        double v;
        const std::uint64_t bits = rng();
        std::memcpy(&v, &bits, sizeof v);
        if (std::isfinite(v)) { values.push_back(v); }
      }
      // And the magnitudes the writers actually see, including the float->double
      // promotion the stream performed and TextWriter::number reproduces.
      for (int i = 0; i < 200000; ++i)
      {
        values.push_back(mz(rng));
        values.push_back(intensity(rng));
        values.push_back(rt(rng));
        values.push_back(irt(rng));
        values.push_back(double(float(intensity(rng))));
        values.push_back(double(float(rt(rng))));
      }

      for (const double v : values)
      {
        ++compared;
        const auto want = viaStream(v, prec), got = viaChars(v, prec);
        if (want != got)
        {
          ++disagree;
          if (first.empty())
          {
            first = "precision " + std::to_string(prec) + ", stream '" + want +
                    "' against to_chars '" + got + "'";
          }
        }
      }
    }

    std::printf("  compared %lld renderings at precision 6, 9 and 10\n", compared);
    if (!first.empty()) { std::printf("  first: %s\n", first.c_str()); }
    check(compared > 3000000, "the sweep is not vacuously small");
    check(disagree == 0, "to_chars general matches the stream's defaultfloat exactly");
  }

  // ------------------------------------------------------------ -out_fragvec

  /// The 78 names in the SEALED order, spelled out here rather than taken from
  /// `PeakGroupScorer::fragvecNames()`.
  ///
  /// The point of this list is that it is an independent transcription of
  /// analysis77/pick/wf_v33_fragvec_contract.md s.5.1: comparing the writer's
  /// header against the generator that produced it would agree with any
  /// consistent reordering, and the consumer joins these columns BY POSITION.
  std::vector<std::string> sealedFragvecNames()
  {
    return {
      "R1_LOGAREA_1", "R1_LOGAREA_2", "R1_LOGAREA_3", "R1_LOGAREA_4",
      "R1_LOGAREA_5", "R1_LOGAREA_6", "R1_LOGAREA_7", "R1_LOGAREA_8",
      "R1_LOGAREA_9", "R1_LOGAREA_10", "R1_LOGAREA_11", "R1_LOGAREA_12",
      "R1_SHARE_1", "R1_SHARE_2", "R1_SHARE_3", "R1_SHARE_4",
      "R1_SHARE_5", "R1_SHARE_6", "R1_SHARE_7", "R1_SHARE_8",
      "R1_SHARE_9", "R1_SHARE_10", "R1_SHARE_11", "R1_SHARE_12",
      "R1_LOGRATIO_1", "R1_LOGRATIO_2", "R1_LOGRATIO_3", "R1_LOGRATIO_4",
      "R1_LOGRATIO_5", "R1_LOGRATIO_6", "R1_LOGRATIO_7", "R1_LOGRATIO_8",
      "R1_LOGRATIO_9", "R1_LOGRATIO_10", "R1_LOGRATIO_11", "R1_LOGRATIO_12",
      "R1_ATAPEX_1", "R1_ATAPEX_2", "R1_ATAPEX_3", "R1_ATAPEX_4",
      "R1_ATAPEX_5", "R1_ATAPEX_6", "R1_ATAPEX_7", "R1_ATAPEX_8",
      "R1_ATAPEX_9", "R1_ATAPEX_10", "R1_ATAPEX_11", "R1_ATAPEX_12",
      "R1_MEAS_1", "R1_MEAS_2", "R1_MEAS_3", "R1_MEAS_4",
      "R1_MEAS_5", "R1_MEAS_6", "R1_MEAS_7", "R1_MEAS_8",
      "R1_MEAS_9", "R1_MEAS_10", "R1_MEAS_11", "R1_MEAS_12",
      "R1_ABSENT_1", "R1_ABSENT_2", "R1_ABSENT_3", "R1_ABSENT_4",
      "R1_ABSENT_5", "R1_ABSENT_6", "R1_ABSENT_7", "R1_ABSENT_8",
      "R1_ABSENT_9", "R1_ABSENT_10", "R1_ABSENT_11", "R1_ABSENT_12",
      "R1_N_MEAS", "R1_N_ABSENT", "R1_N_PRESENT", "R1_LOGTOT",
      "R1_LIB_CORR", "R1_LIB_CORR_LOO"
    };
  }

  std::vector<std::string> splitTabs(const std::string& line)
  {
    std::vector<std::string> f;
    std::size_t at = 0;
    for (;;)
    {
      const std::size_t t = line.find('\t', at);
      if (t == std::string::npos) { f.push_back(line.substr(at)); return f; }
      f.push_back(line.substr(at, t - at));
      at = t + 1;
    }
  }

  std::vector<std::string> readLines(const std::string& path)
  {
    std::vector<std::string> out;
    std::ifstream in(path);
    if (!in) { throw std::runtime_error("cannot read " + path); }
    std::string line;
    while (std::getline(in, line)) { out.push_back(line); }
    return out;
  }

  /// `-out_fragvec` is a JOIN, so the things that can break it are the header
  /// order, the key, and the precision -- not the arithmetic, which lives in
  /// the scorer. Each is checked against something other than the code that
  /// produces it.
  void caseFragvec(const std::string& dir)
  {
    const auto lib = awkwardLibrary();
    const auto& p = lib.precursors();

    // Precursor 0 and precursor 4 reconstruct the SAME Precursor.Id (the
    // fixture gives every precursor the same modified sequence and charge
    // 1 + i % 4), which is exactly the case decision 4 turns on: the ordinal
    // resets on the id STRING, not on the library index, so these two must
    // share one block and keep counting.
    const std::vector<std::uint32_t> from = {0, 0, 4, 0, 1, 2, 1};
    const std::vector<bool>          dec  = {false, false, false, true, true, false, true};
    const std::vector<std::string> want_id = {
      "PEPT(Phospho)IDEK1", "PEPT(Phospho)IDEK1", "PEPT(Phospho)IDEK1",
      "PEPT(Phospho)IDEK1", "PEPT(Phospho)IDEK2", "PEPT(Phospho)IDEK3",
      "PEPT(Phospho)IDEK2" };
    const std::vector<long long> want_ord = {0, 1, 2, 0, 0, 0, 0};

    ODIA::PeakGroupScorer::Result r;
    for (std::size_t i = 0; i < from.size(); ++i)
    {
      ODIA::PeakGroupScorer::PeakGroup g;
      g.precursor = from[i];
      g.decoy = dec[i];
      r.groups.push_back(std::move(g));
    }
    check(p.charge[0] == 1 && p.charge[4] == 1,
          "the fixture really does give two library rows the same Precursor.Id");

    // Values chosen so a precision that is too low is a WRONG FILE rather than
    // a rounded one: each needs more than six significant digits to come back
    // as the same float32.
    const std::vector<float> awkward = {
      0.0f, -0.0f, 1.0f, 1.0f / 3.0f, 2.0f / 3.0f, 0.1f, 1e-5f, 1e20f, 1e-20f,
      9.99999905f, 123456.789f, 0.000999999931f, 16777217.0f, 1.00000012f,
      -3.14159274f, 5.87747175e-39f, NA
    };
    for (std::size_t i = 0; i < r.groups.size() * ODIA::PeakGroupScorer::N_FRAGVEC; ++i)
    { r.fragvec.push_back(awkward[i % awkward.size()]); }

    const std::string path = dir + "/tsv_fragvec.tsv";
    ODIA::writeFragvecTsv(path, lib, r);
    const auto lines = readLines(path);

    check(lines.size() == r.groups.size() + 1,
          "one header plus one row per group (" + std::to_string(lines.size()) + ")");

    const auto head = splitTabs(lines.at(0));
    const auto sealed = sealedFragvecNames();
    check(sealed.size() == ODIA::PeakGroupScorer::N_FRAGVEC,
          "the sealed contract really is 78 columns");
    bool head_ok = head.size() == 3 + sealed.size() && head[0] == "Precursor.Id" &&
                   head[1] == "Decoy" && head[2] == "Ordinal";
    std::string head_first;
    for (std::size_t j = 0; head_ok && j < sealed.size(); ++j)
    {
      if (head[3 + j] != sealed[j])
      { head_ok = false; head_first = head[3 + j] + " != " + sealed[j]; }
    }
    if (!head_first.empty()) { std::printf("  first: %s\n", head_first.c_str()); }
    check(head_ok, "the header is the sealed s.5.1 order, independently transcribed");

    long long bad_key = 0, bad_width = 0, bad_roundtrip = 0, nan_cells = 0;
    std::string first;
    for (std::size_t i = 1; i < lines.size(); ++i)
    {
      const auto f = splitTabs(lines[i]);
      if (f.size() != 3 + ODIA::PeakGroupScorer::N_FRAGVEC) { ++bad_width; continue; }
      if (f[0] != want_id[i - 1] ||
          f[1] != std::to_string(dec[i - 1] ? 1 : 0) ||
          f[2] != std::to_string(want_ord[i - 1]))
      {
        if (first.empty())
        { first = "row " + std::to_string(i) + ": " + f[0] + "/" + f[1] + "/" + f[2]; }
        ++bad_key;
      }
      for (std::size_t j = 0; j < ODIA::PeakGroupScorer::N_FRAGVEC; ++j)
      {
        const float wrote = r.fragvec[(i - 1) * ODIA::PeakGroupScorer::N_FRAGVEC + j];
        const float back = std::strtof(f[3 + j].c_str(), nullptr);
        if (std::isnan(wrote))
        {
          ++nan_cells;
          if (f[3 + j] != "nan" || !std::isnan(back)) { ++bad_roundtrip; }
          continue;
        }
        std::uint32_t a = 0, b = 0;
        std::memcpy(&a, &wrote, 4);
        std::memcpy(&b, &back, 4);
        if (a != b)
        {
          if (first.empty())
          { first = "cell " + head[3 + j] + " wrote " + f[3 + j]; }
          ++bad_roundtrip;
        }
      }
    }
    if (!first.empty()) { std::printf("  first: %s\n", first.c_str()); }
    check(bad_width == 0, "every row carries exactly 81 fields");
    check(bad_key == 0,
          "the (Precursor.Id, Decoy, Ordinal) key is the one the join expects, "
          "and the ordinal resets on the id string rather than the library index");
    check(nan_cells > 0, "the fixture actually exercised the NaN rendering");
    check(bad_roundtrip == 0,
          "every value round-trips to the same float32 -- the precision is not 6");

    // The guard that says the export is still parallel to `groups`. It is the
    // only thing standing between a future reorder and a silently misaligned
    // 10.9 GB file, so it has to actually fire.
    bool threw = false;
    r.fragvec.pop_back();
    try { ODIA::writeFragvecTsv(dir + "/tsv_fragvec_short.tsv", lib, r); }
    catch (const std::exception&) { threw = true; }
    check(threw, "a fragvec that is not 78 floats per group is refused, not written");
  }

  // --------------------------------------------------------- -fragvec_scores

  /// The block `-fragvec_scores` appends to the classifier is defined by WHICH
  /// columns it takes, by the guard that is the only filter on them, by the
  /// permutation its control reads through, and by the one function both the
  /// fit and the -out writer read. None of that is arithmetic, so none of it
  /// is covered by gate A4; all of it is checked here against things other
  /// than the code that produces it.
  void caseFragvecScores(const std::string& /*dir*/)
  {
    using PGS = ODIA::PeakGroupScorer;
    using AB = PGS::AppendedBlock;
    const auto sealed = sealedFragvecNames();

    // 1. The schema: the whole sealed contract, 78 columns in contract order,
    //    against the independent transcription, not against fragvecNames().
    const auto& fv = PGS::appendedNames(AB::Fragvec);
    const auto& pm = PGS::appendedNames(AB::Permuted);
    bool names_ok = sealed.size() == 78 && fv.size() == sealed.size() && pm.size() == sealed.size();
    for (std::size_t k = 0; names_ok && k < sealed.size(); ++k)
    {
      names_ok = fv[k] == "var_fv_" + sealed[k] && pm[k] == "var_fvperm_" + sealed[k];
    }
    check(names_ok, "var_fv_* / var_fvperm_* are the 78 sealed contract names, in order");
    check(PGS::appendedNames(AB::None).empty(), "off appends nothing");

    // 1b. The guard: the only filter. Built on rows whose columns are each one
    //     of five shapes, chosen per column by index, so every shape lands on
    //     indicator AND value columns alike -- and the MEAS/ABSENT families get
    //     no special treatment either way.
    {
      const std::size_t n = 6;
      std::vector<float> block(n * PGS::N_FRAGVEC, 0.0f);
      std::vector<std::size_t> want_kept;
      for (std::size_t j = 0; j < PGS::N_FRAGVEC; ++j)
      {
        const int shape = static_cast<int>(j % 5);
        for (std::size_t r = 0; r < n; ++r)
        {
          float v = 0.0f;
          switch (shape)
          {
            case 0: v = 1.0f; break;                                     // constant
            case 1: v = NA; break;                                       // all NaN
            case 2: v = (r == 2) ? 0.5f : NA; break;                     // one finite value
            case 3: v = (r % 2) ? 1.0f : 0.0f; break;                    // 0/1 indicator
            case 4: v = (r == 0) ? NA : static_cast<float>(r) * 0.25f; break;  // NaN + varying
          }
          block[r * PGS::N_FRAGVEC + j] = v;
        }
        if (shape == 3 || shape == 4) { want_kept.push_back(j); }
      }
      const auto kept = PGS::fragvecKeptColumns(block, n);
      check(kept == want_kept, "the guard keeps exactly the columns with two different finite values (" +
                               std::to_string(kept.size()) + " of 78)");
      // Constant 1 on R1_MEAS_1 (every precursor has a rank-1 fragment) is
      // dropped; a varying R1_MEAS_6 is kept -- by value, not by name.
      std::vector<float> meas(n * PGS::N_FRAGVEC, 0.25f);
      std::size_t i_meas1 = 0, i_meas6 = 0;
      for (std::size_t j = 0; j < sealed.size(); ++j)
      {
        if (sealed[j] == "R1_MEAS_1") { i_meas1 = j; }
        if (sealed[j] == "R1_MEAS_6") { i_meas6 = j; }
      }
      for (std::size_t r = 0; r < n; ++r)
      {
        meas[r * PGS::N_FRAGVEC + i_meas1] = 1.0f;
        meas[r * PGS::N_FRAGVEC + i_meas6] = (r < 3) ? 1.0f : 0.0f;
      }
      const auto kept2 = PGS::fragvecKeptColumns(meas, n);
      check(kept2.size() == 1 && kept2[0] == i_meas6,
            "a constant R1_MEAS_1 is dropped and a varying R1_MEAS_6 kept, by value");
      // Row order cannot change the answer (so `permuted` keeps what `final` keeps).
      std::vector<float> rev(block.size());
      for (std::size_t r = 0; r < n; ++r)
      {
        std::copy_n(block.begin() + static_cast<std::ptrdiff_t>((n - 1 - r) * PGS::N_FRAGVEC),
                    PGS::N_FRAGVEC, rev.begin() + static_cast<std::ptrdiff_t>(r * PGS::N_FRAGVEC));
      }
      check(PGS::fragvecKeptColumns(rev, n) == kept, "the guard is invariant to row order");
      bool threw = false;
      try { (void)PGS::fragvecKeptColumns(block, n + 1); }
      catch (const std::exception&) { threw = true; }
      check(threw, "the guard refuses an array that is not 78 floats per row");
    }

    // 2. Fragvec: the value is the exported float32 widened, cell for cell,
    //    NaN included.
    PGS::Result r;
    const std::vector<std::uint32_t> prec = {3, 3, 1, 4, 0, 2, 2};
    for (std::size_t i = 0; i < prec.size(); ++i)
    {
      PGS::PeakGroup g;
      g.precursor = prec[i];
      g.apex_rt = 100.0f + 1.385f * static_cast<float>(i);
      g.decoy = (i % 2) == 1;
      r.groups.push_back(std::move(g));
    }
    for (std::size_t i = 0; i < r.groups.size() * PGS::N_FRAGVEC; ++i)
    {
      r.fragvec.push_back(i % 11 == 0 ? NA
                                      : 1.0f / static_cast<float>(i + 3) + 0.1f * static_cast<float>(i % 7));
    }
    r.appended = AB::Fragvec;
    const auto same = [](float src, double got) {
      return std::isnan(src) ? std::isnan(got) : got == static_cast<double>(src);
    };
    long long bad = 0;
    for (std::size_t row = 0; row < r.groups.size(); ++row)
    {
      for (std::size_t j = 0; j < PGS::N_FRAGVEC; ++j)
      {
        if (!same(r.fragvec[row * PGS::N_FRAGVEC + j], PGS::appendedRaw(r, row, j))) { ++bad; }
      }
    }
    check(bad == 0, "every var_fv_* cell is its own row's -out_fragvec cell widened to double");

    // 3. The permutation: a bijection, a function of (n, seed) only, moved by
    //    the seed, and not the identity at a realistic size.
    const auto p0 = PGS::fragvecPermutation(100000, 0);
    std::vector<std::uint32_t> sorted = p0;
    std::sort(sorted.begin(), sorted.end());
    bool bijection = sorted.size() == 100000;
    for (std::size_t k = 0; bijection && k < sorted.size(); ++k) { bijection = sorted[k] == k; }
    check(bijection, "fragvecPermutation(n, seed) is a permutation of 0..n-1");
    check(PGS::fragvecPermutation(100000, 0) == p0, "the permutation is deterministic");
    check(PGS::fragvecPermutation(100000, 7) != p0, "-null_feature_seed moves the permutation");
    std::size_t fixed_points = 0;
    for (std::size_t k = 0; k < p0.size(); ++k) { fixed_points += p0[k] == k; }
    check(fixed_points < 20, "the permutation is not near the identity (" +
                             std::to_string(fixed_points) + " fixed points of 100,000)");

    // 4. Permuted reads WHOLE ROWS through it: every cell of a row comes from
    //    the same source row (so NaN patterns and within-row covariance travel
    //    together), and the label plays no part in which row that is.
    PGS::Result q = r;
    q.appended = AB::Permuted;
    q.appended_perm = PGS::fragvecPermutation(q.groups.size(), 3);
    long long bad_perm = 0;
    for (std::size_t row = 0; row < q.groups.size(); ++row)
    {
      const std::size_t src = q.appended_perm[row];
      for (std::size_t j = 0; j < PGS::N_FRAGVEC; ++j)
      {
        if (!same(q.fragvec[src * PGS::N_FRAGVEC + j], PGS::appendedRaw(q, row, j))) { ++bad_perm; }
      }
    }
    check(bad_perm == 0, "every var_fvperm_* row is one whole source row of the export");
    PGS::Result flipped = q;
    for (auto& g : flipped.groups) { g.decoy = !g.decoy; }
    flipped.appended_perm = PGS::fragvecPermutation(flipped.groups.size(), 3);
    check(flipped.appended_perm == q.appended_perm,
          "flipping every Decoy flag leaves the permutation unchanged");
  }
}

int main(int argc, char** argv)
{
  const std::string which = argc > 1 ? argv[1] : "";
  const std::string dir = argc > 2 ? argv[2] : ".";
  std::printf("case: %s\n", which.c_str());

  try
  {
    if (which == "files") { caseFiles(dir); }
    else if (which == "formatting") { caseFormatting(); }
    else if (which == "fragvec") { caseFragvec(dir); }
    else if (which == "fragvec_scores") { caseFragvecScores(dir); }
    else
    {
      std::fprintf(stderr,
                   "usage: odia_tsv_writers <files|formatting|fragvec|fragvec_scores> "
                   "[directory]\n");
      return 2;
    }
  }
  catch (const std::exception& e)
  {
    std::fprintf(stderr, "  FAIL threw: %s\n", e.what());
    return 1;
  }

  std::printf(failures ? "FAILED (%d)\n" : "ok (%d)\n", failures);
  return failures ? 1 : 0;
}
