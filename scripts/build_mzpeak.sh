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

MZPEAK_COMMIT="${MZPEAK_COMMIT:-80fff0152328dccc30a0e83a057e2e01d82d5f96}"
src="${ODIA_SCRATCH}/src/mzpeak"
build="${ODIA_SCRATCH}/build/mzpeak"

mkdir -p "${ODIA_SCRATCH}/src"
if [[ ! -d "${src}/.git" ]]; then
  echo "==> cloning mzpeak"
  git clone https://github.com/OpenMS/mzpeak.git "${src}"
fi
git -C "${src}" fetch --all --quiet
git -C "${src}" checkout --quiet --force "${MZPEAK_COMMIT}"
git -C "${src}" clean -qfd
echo "==> mzpeak at $(git -C "${src}" log -1 --format='%h %ci %s')"

# mzpeak is developed against clang (see its flake.nix, which pins clang_22).
# GCC rejects one construct that clang accepts, so we carry local patches
# rather than switching compilers -- ODIA, OpenMS and mzpeak must all be built
# with the same toolchain. These should be sent upstream; drop each patch here
# once it lands.
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

echo "==> building mzpeak"
meson compile -C "${build}"

echo "==> installing mzpeak to ${ODIA_MZPEAK}"
meson install -C "${build}"

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
