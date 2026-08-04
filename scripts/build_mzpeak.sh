#!/usr/bin/env bash
#
# Build and install OpenMS/mzpeak, the mzPeak format reader ODIA uses for
# spectrum access.
#
# Pinned to a specific commit: mzpeak's README states the API is not yet
# stable ("no stability is guaranteed at this point", API freeze targeted for
# summer 2026), so ODIA must build against a known revision rather than trunk.

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${here}/env.sh"

# okohlbacher/mzpeak-openms is a fork of OpenMS/mzpeak carrying the work ODIA
# needs and upstream does not yet have:
#
#   - reads the split-metadata (v2) layout, which is what mzpeak-convert >= 0.7.0
#     writes and therefore what all of our example archives use. Upstream trunk
#     throws "metadata file missing or does not have the spectrum group" on them.
#   - exposes per-spectrum retention time, precursors/isolation windows and ion
#     mobility. Upstream's Spectrum exposes only m/z, intensity and MS level,
#     which is not enough to drive a DIA extraction.
#   - has a writer, which is how ODIA emits mzPeak chromatograms.
#
# Pinned to a commit because this is an active WIP branch, not a release.
MZPEAK_REPO="${MZPEAK_REPO:-https://github.com/okohlbacher/mzpeak-openms.git}"
MZPEAK_COMMIT="${MZPEAK_COMMIT:-7c1b6e53470a1f674a2c579a201fd1070b74efd9}"
src="${ODIA_SCRATCH}/src/mzpeak-openms"
build="${ODIA_SCRATCH}/build/mzpeak"

mkdir -p "${ODIA_SCRATCH}/src"
if [[ ! -d "${src}/.git" ]]; then
  echo "==> cloning mzpeak"
  git clone "${MZPEAK_REPO}" "${src}"
fi
git -C "${src}" fetch --all --quiet
git -C "${src}" checkout --quiet --force "${MZPEAK_COMMIT}"
git -C "${src}" clean -qfd
echo "==> mzpeak at $(git -C "${src}" log -1 --format='%h %ci %s')"

# mzpeak is developed against clang (see its flake.nix, which pins clang_22),
# so a construct clang accepts and GCC rejects can appear at any time. Local
# patches go in patches/ and are dropped once upstream takes them -- the
# defaulted-move-constructor fix was carried here and is now upstream, so this
# loop is currently a no-op and the directory is empty.
for p in "${here}/../patches"/mzpeak-*.patch; do
  [[ -e "${p}" ]] || continue
  echo "==> applying $(basename "${p}")"
  git -C "${src}" apply "${p}"
done

rm -rf "${build}"
echo "==> configuring mzpeak"
meson setup "${build}" "${src}" \
  --prefix "${ODIA_MZPEAK}" \
  --libdir lib \
  --buildtype release \
  -Dwerror=false \
  -Dcpp_args="-I${ODIA_ENV}/include" \
  -Dcpp_link_args="-L${ODIA_ENV}/lib -Wl,-rpath,${ODIA_ENV}/lib"

# Build only the library and tools. The test suite does not currently compile
# against Arrow 23 (parquet_writer_test.cpp calls FileReader::ReadTable() with a
# signature this Arrow does not have) -- a WIP-branch issue, not ours.
echo "==> building mzpeak"
ninja -C "${build}" libmzpeak.so libmzpeak.a mzp-inspect mzp-bench read_spectra

echo "==> installing mzpeak to ${ODIA_MZPEAK}"
meson install -C "${build}" --no-rebuild

# Repair the installed header layout.
#
# mzpeak's meson.build calls install_headers(subdir: 'mzpeak'), which flattens
# every header into a single directory. The headers themselves use nested
# includes ("mzpeak/io/archive.h", "mzpeak/util/parquet.h", ...), and the
# umbrella header belongs at <mzpeak.h>, so the installed tree does not
# compile as shipped. Mirror the source layout instead, which does.
# Upstream fix would be install_headers(..., preserve_path: true) plus keeping
# include/mzpeak.h at the include root.
echo "==> correcting installed header layout"
rm -rf "${ODIA_MZPEAK}/include"
mkdir -p "${ODIA_MZPEAK}/include"
cp -a "${src}/include/." "${ODIA_MZPEAK}/include/"
test -f "${ODIA_MZPEAK}/include/mzpeak.h"
test -f "${ODIA_MZPEAK}/include/mzpeak/io/archive.h"
