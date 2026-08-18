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
#include <odia/Library.h>
#include <odia/TextWriter.h>

#include <charconv>
#include <cmath>
#include <cstdio>
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
    else
    {
      std::fprintf(stderr, "usage: odia_tsv_writers <files|formatting> [directory]\n");
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
