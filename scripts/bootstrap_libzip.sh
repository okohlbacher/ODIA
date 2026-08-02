#!/usr/bin/env bash
#
# Build libzip from source into the ODIA environment prefix.
#
# mzpeak's meson.build requires libzip >= 1.11.4; conda-forge only ships
# 1.11.2. Installing into ${ODIA_ENV} rather than a separate prefix keeps a
# single libzip on the loader path, so there is no chance of two versions
# being mixed at runtime.

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${here}/env.sh"

LIBZIP_VERSION="${LIBZIP_VERSION:-1.11.4}"
src_root="${ODIA_SCRATCH}/src"
build_root="${ODIA_SCRATCH}/build"
mkdir -p "${src_root}" "${build_root}"

tarball="libzip-${LIBZIP_VERSION}.tar.gz"
if [[ ! -d "${src_root}/libzip-${LIBZIP_VERSION}" ]]; then
  echo "==> fetching libzip ${LIBZIP_VERSION}"
  curl -sSL -o "${src_root}/${tarball}" \
    "https://libzip.org/download/${tarball}"
  tar -xzf "${src_root}/${tarball}" -C "${src_root}"
fi

echo "==> configuring libzip"
cmake -S "${src_root}/libzip-${LIBZIP_VERSION}" \
      -B "${build_root}/libzip" \
      -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="${ODIA_ENV}" \
      -DCMAKE_PREFIX_PATH="${ODIA_ENV}" \
      -DCMAKE_INSTALL_RPATH="${ODIA_ENV}/lib" \
      -DBUILD_SHARED_LIBS=ON \
      -DBUILD_TOOLS=OFF \
      -DBUILD_REGRESS=OFF \
      -DBUILD_EXAMPLES=OFF \
      -DBUILD_DOC=OFF \
      -DENABLE_GNUTLS=OFF \
      -DENABLE_MBEDTLS=OFF

echo "==> building libzip"
cmake --build "${build_root}/libzip" --parallel "$(nproc)"
cmake --install "${build_root}/libzip"

echo "==> installed:"
pkg-config --modversion libzip
