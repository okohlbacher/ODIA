#!/usr/bin/env bash
#
# Build and install OpenMS 3.6 as ODIA's base library.
#
# ODIA uses OpenMS but never modifies it: this script clones a pinned upstream
# commit and installs it to ${ODIA_OPENMS}, which is then treated as read-only.
# Re-run only to move to a new pinned commit.
#
# CMAKE_MODULE_PATH is passed explicitly to work around an upstream bug:
# cmake/FindONNXRuntime.cmake sits in cmake/, but OpenMS only appends
# cmake/Modules and cmake/Windows to CMAKE_MODULE_PATH, so with WITH_ONNX=ON
# find_package(ONNXRuntime REQUIRED) never sees the Find module, falls through
# to config mode and fails looking for ONNXRuntimeConfig.cmake. Prepending the
# directory on the command line survives OpenMS's later list(APPEND ...).
# Report upstream; drop this once fixed.
#
# WITH_ONNX=ON is required, not optional: it enables OpenMS's PeptDeep bindings
# (OpenMS/ML/PEPTDEEP) and makes the build download the three PeptDeep ONNX
# models from archive.openms.de against pinned SHA256s. Those models are the
# ones ODIA's library generator predicts with, and they are OpenMS artefacts --
# obtaining them through OpenMS's own build is what keeps ODIA from vendoring
# weights of its own.
#
# Built against the conda-forge environment from bootstrap_env.sh rather than
# vcpkg (OPENMS_USE_VCPKG defaults to OFF): the environment already supplies
# Arrow/Parquet 23, which is OpenMS 3.6's binding new requirement, and using
# one environment for OpenMS, mzpeak and ODIA keeps the C++ ABI consistent.

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${here}/env.sh"

# OpenMS develop as of 2026-07-31; reports itself as 3.6.0.
OPENMS_COMMIT="${OPENMS_COMMIT:-552802a781c15f909c33e044514c517ebc53df25}"
src="${ODIA_SCRATCH}/src/OpenMS"
build="${ODIA_SCRATCH}/build/OpenMS"

mkdir -p "${ODIA_SCRATCH}/src"
if [[ ! -d "${src}/.git" ]]; then
  echo "==> cloning OpenMS"
  git clone --filter=blob:none https://github.com/OpenMS/OpenMS.git "${src}"
fi
git -C "${src}" fetch --all --quiet
git -C "${src}" checkout --quiet --force "${OPENMS_COMMIT}"
echo "==> OpenMS at $(git -C "${src}" log -1 --format='%h %ci %s')"

# env.sh puts ${ODIA_OPENMS} on CMAKE_PREFIX_PATH so that ODIA can find the
# installed OpenMS. When building OpenMS itself that is actively harmful: the
# configure step then rediscovers packages exported by its own previous install
# (opentims, whose exported target references a zstd target that does not exist
# in this scope) and fails. Restricting CMAKE_PREFIX_PATH is not enough --
# find_package also searches CMAKE_INSTALL_PREFIX implicitly -- so that lookup
# is disabled too. Non-destructive: the previous install stays usable if this
# build fails.
export CMAKE_PREFIX_PATH="${ODIA_ENV}"

echo "==> configuring OpenMS"
cmake -S "${src}" -B "${build}" \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${ODIA_OPENMS}" \
  -DCMAKE_PREFIX_PATH="${ODIA_ENV}" \
  -DCMAKE_C_COMPILER="${CC}" \
  -DCMAKE_CXX_COMPILER="${CXX}" \
  -DCMAKE_INSTALL_RPATH="${ODIA_OPENMS}/lib;${ODIA_ENV}/lib" \
  -DCMAKE_BUILD_WITH_INSTALL_RPATH=OFF \
  -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON \
  -DOPENMS_USE_VCPKG=OFF \
  -DWITH_GUI=OFF \
  -DENABLE_TUTORIALS=OFF \
  -DENABLE_DOCS=OFF \
  -DENABLE_UPDATE_CHECK=OFF \
  -DPYOPENMS=OFF \
  -DENABLE_PREPARE_DOCS=OFF \
  -DARROW_USE_STATIC=OFF \
  -DBOOST_USE_STATIC=OFF \
  -DWITH_THERMO_RAW=OFF \
  -DWITH_ONNX=ON \
  -DCMAKE_MODULE_PATH="${src}/cmake" \
  -DCMAKE_FIND_USE_INSTALL_PREFIX=OFF

echo "==> building OpenMS"
cmake --build "${build}" --parallel "$(nproc)"

echo "==> installing OpenMS to ${ODIA_OPENMS}"
cmake --install "${build}"

echo "==> done"
"${ODIA_OPENMS}/bin/FileInfo" --help 2>&1 | head -5 || true
