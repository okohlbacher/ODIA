#!/usr/bin/env bash
#
# Configure, build and test ODIA itself.
#
# Assumes the dependencies are already in place:
#   scripts/bootstrap_env.sh
#   scripts/bootstrap_libzip.sh
#   scripts/build_mzpeak.sh
#   scripts/build_openms.sh

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "${here}/.." && pwd)"
# shellcheck source=/dev/null
source "${here}/env.sh"

build="${ODIA_BUILD:-${ODIA_SCRATCH}/build/odia}"

for dep in "${ODIA_OPENMS}/lib" "${ODIA_MZPEAK}/lib"; do
  if [[ ! -d "${dep}" ]]; then
    echo "missing dependency: ${dep} -- run the bootstrap/build scripts first" >&2
    exit 1
  fi
done

echo "==> configuring ODIA in ${build}"
cmake -S "${repo}" -B "${build}" \
  -G Ninja \
  -DCMAKE_BUILD_TYPE="${ODIA_BUILD_TYPE:-Release}" \
  -DCMAKE_C_COMPILER="${CC}" \
  -DCMAKE_CXX_COMPILER="${CXX}" \
  -DCMAKE_PREFIX_PATH="${ODIA_MZPEAK};${ODIA_OPENMS};${ODIA_ENV}" \
  -DCMAKE_INSTALL_PREFIX="${ODIA_INSTALL:-${repo}/install}" \
  -DCMAKE_BUILD_RPATH="${ODIA_MZPEAK}/lib;${ODIA_OPENMS}/lib;${ODIA_ENV}/lib" \
  -DCMAKE_INSTALL_RPATH="${ODIA_MZPEAK}/lib;${ODIA_OPENMS}/lib;${ODIA_ENV}/lib"

echo "==> building ODIA"
cmake --build "${build}" --parallel "$(nproc)"

echo "==> testing ODIA"
ctest --test-dir "${build}" --output-on-failure

echo "==> built tools:"
ls -1 "${build}"/ODIA* 2>/dev/null || true
