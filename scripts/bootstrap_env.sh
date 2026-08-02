#!/usr/bin/env bash
#
# Create the ODIA build/runtime environment.
#
# ODIA links against OpenMS 3.6 and OpenMS/mzpeak. Both need dependencies that
# are newer than Ubuntu 24.04 provides:
#
#   - OpenMS 3.6 requires Apache Arrow/Parquet >= 23 (cmake_findExternalLibs.cmake)
#   - mzpeak requires C++23 (std::print, std::ranges::to), Arrow >= 23,
#     Boost >= 1.89, libzip >= 1.11.4
#
# The system compiler here is GCC 13, which lacks <print> and std::ranges::to.
# So we provision a self-contained conda-forge environment: one toolchain and
# one set of shared libraries for OpenMS, mzpeak and ODIA, which keeps the C++
# ABI consistent across all three.
#
# libzip is the one exception: conda-forge tops out at 1.11.2 and mzpeak needs
# >= 1.11.4, so it is built from source into the same prefix (see
# bootstrap_libzip.sh).
#
# Idempotent: safe to re-run. Creates nothing outside $ODIA_OPT.

set -euo pipefail

ODIA_ROOT="${ODIA_ROOT:-/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer}"
ODIA_OPT="${ODIA_OPT:-${ODIA_ROOT}/opt}"
ENV_PREFIX="${ENV_PREFIX:-${ODIA_OPT}/env}"
MM="${ODIA_OPT}/tools/bin/micromamba"

GCC_VERSION="${GCC_VERSION:-14}"
ARROW_VERSION="${ARROW_VERSION:-23}"
BOOST_VERSION="${BOOST_VERSION:-1.89}"

mkdir -p "${ODIA_OPT}/tools"

# ---------------------------------------------------------------- micromamba
if [[ ! -x "${MM}" ]]; then
  echo "==> fetching micromamba"
  tmp="$(mktemp -d)"
  curl -sSL -o "${tmp}/micromamba.tar.bz2" \
    "https://micro.mamba.pm/api/micromamba/linux-64/latest"
  tar -xjf "${tmp}/micromamba.tar.bz2" -C "${ODIA_OPT}/tools" bin/micromamba
  rm -rf "${tmp}"
fi
"${MM}" --version

export MAMBA_ROOT_PREFIX="${ODIA_OPT}/mamba"

# ---------------------------------------------------------------- the env
# Package set, grouped by who needs it:
#
#   toolchain      GCC 14 (C++23) + build drivers. meson/ninja build mzpeak,
#                  cmake/ninja build OpenMS and ODIA.
#   arrow/parquet  required >= 23 by both OpenMS and mzpeak
#   boost          required >= 1.89 by mzpeak; OpenMS needs many components
#   OpenMS only    xerces-c, libsvm, eigen, coin-or (LP solver), hdf5, curl,
#                  libxml2, qt6 (Core/Network are used even with WITH_GUI=OFF)
#   shared         zlib, bzip2, zstd, xz, sqlite, openssl
echo "==> creating environment at ${ENV_PREFIX}"
"${MM}" create --yes --prefix "${ENV_PREFIX}" \
  -c conda-forge \
  "gcc_linux-64=${GCC_VERSION}.*" \
  "gxx_linux-64=${GCC_VERSION}.*" \
  "sysroot_linux-64=2.28" \
  cmake ninja meson pkg-config make patchelf \
  "libarrow=${ARROW_VERSION}.*" \
  "libparquet=${ARROW_VERSION}.*" \
  "libarrow-dataset=${ARROW_VERSION}.*" \
  "libarrow-acero=${ARROW_VERSION}.*" \
  "libboost-devel=${BOOST_VERSION}.*" \
  "libboost-headers=${BOOST_VERSION}.*" \
  xerces-c libsvm eigen hdf5 libcurl libxml2 libxml2-devel \
  coin-or-cbc coin-or-clp coin-or-cgl coin-or-osi coin-or-utils \
  qt6-main \
  zlib bzip2 zstd xz sqlite openssl

echo "==> environment created"

# Verify the environment actually provides every header OpenMS and mzpeak look
# for. This is not paranoia: an interrupted micromamba transaction leaves the
# prefix silently incomplete -- packages vanish while conda-meta still looks
# healthy -- and the failure then surfaces much later as a confusing CMake
# error. Re-running this script repairs such a prefix.
echo "==> verifying environment"
fail=0
check() {
  if [[ -e "$2" ]]; then
    printf '  ok    %s\n' "$1"
  else
    printf '  MISSING %-12s (expected %s)\n' "$1" "$2"
    fail=1
  fi
}
check xerces-c "${ENV_PREFIX}/include/xercesc/util/XercesVersion.hpp"
check libxml2  "${ENV_PREFIX}/include/libxml2/libxml/parser.h"
check boost    "${ENV_PREFIX}/include/boost/version.hpp"
check arrow    "${ENV_PREFIX}/include/arrow/api.h"
check parquet  "${ENV_PREFIX}/include/parquet/api/reader.h"
check eigen    "${ENV_PREFIX}/include/eigen3/Eigen/Core"
check libsvm   "${ENV_PREFIX}/include/svm.h"
check hdf5     "${ENV_PREFIX}/include/hdf5.h"
check curl     "${ENV_PREFIX}/include/curl/curl.h"
check zlib     "${ENV_PREFIX}/include/zlib.h"
check bzip2    "${ENV_PREFIX}/include/bzlib.h"
check coin-or  "${ENV_PREFIX}/include/coin/CbcModel.hpp"
check qt6      "${ENV_PREFIX}/lib/cmake/Qt6/Qt6Config.cmake"
check g++      "${ENV_PREFIX}/bin/x86_64-conda-linux-gnu-g++"

if [[ "${fail}" -ne 0 ]]; then
  echo "environment is incomplete -- re-run this script to repair it" >&2
  exit 1
fi
echo "==> environment OK; next: bootstrap_libzip.sh, build_mzpeak.sh, build_openms.sh"
