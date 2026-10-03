// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// Dump the peaks of chosen spectra exactly as ODIA's SpectrumSource serves
// them, under -tof_calibration run or frame. The probe that compares them with
// the Bruker SDK lives outside the tree (it needs the vendor library and the
// .d); this only fixes WHAT ODIA reads, through the same openRun the tool uses.
//
// usage: odia_tof_dump <run.mzpeak> <run|frame> <out.bin> <file_index>...
// out.bin, per spectrum: u64 file index, u8 ms level, f64 rt, u64 n,
//                        n*f64 m/z, n*f32 intensity, n*f32 1/K0 (or n zeros)

#include <odia/SpectrumSource.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
  if (argc < 5)
  {
    std::fprintf(stderr, "usage: %s run.mzpeak run|frame out.bin file_index...\n", argv[0]);
    return 2;
  }
  const std::string mode = argv[2];
  if (mode != "run" && mode != "frame") { std::fprintf(stderr, "mode must be run or frame\n"); return 2; }
  ODIA::RunOptions opt;
  opt.per_frame_tof_calibration = mode == "frame";
  auto run = ODIA::openRun(argv[1], opt);

  std::FILE* f = std::fopen(argv[3], "wb");
  if (!f) { return 3; }
  const auto& ms1 = run->ms1Spectra();
  const auto& ms2 = run->spectra();
  int missing = 0;
  for (int a = 4; a < argc; ++a)
  {
    const std::uint64_t want = std::strtoull(argv[a], nullptr, 10);
    std::vector<ODIA::SpectrumPeaks> out;
    std::uint8_t lvl = 0;
    double rt = -1.0;
    for (std::size_t k = 0; k < ms1.size() && !lvl; ++k)
    {
      if (ms1[k].index == want) { run->ms1Peaks(k, k + 1, out); lvl = 1; rt = ms1[k].retention_time; }
    }
    // A co-packed MS2 frame has one entry per window and the same peaks in each.
    for (std::size_t k = 0; k < ms2.size() && !lvl; ++k)
    {
      if (ms2[k].index == want) { run->peaks(k, k + 1, out); lvl = 2; rt = ms2[k].retention_time; }
    }
    if (!lvl || out.empty()) { ++missing; std::fprintf(stderr, "spectrum %llu not found\n", (unsigned long long)want); continue; }
    const auto& p = out.front();
    const std::uint64_t n = p.mz.size();
    std::fwrite(&want, 8, 1, f); std::fwrite(&lvl, 1, 1, f); std::fwrite(&rt, 8, 1, f); std::fwrite(&n, 8, 1, f);
    std::fwrite(p.mz.data(), 8, n, f);
    std::fwrite(p.intensity.data(), 4, n, f);
    if (p.ion_mobility.size() == n) { std::fwrite(p.ion_mobility.data(), 4, n, f); }
    else { std::vector<float> z(n, 0.0f); std::fwrite(z.data(), 4, n, f); }
  }
  std::fclose(f);
  return missing ? 1 : 0;
}
