# Shell fragment: source this to get the ODIA build/runtime environment.
#
#   source scripts/env.sh
#
# Sets the compiler (GCC 14 from conda-forge, needed for mzpeak's C++23),
# and puts the ODIA prefixes on CMake's, pkg-config's and the loader's
# search paths.

ODIA_ROOT="${ODIA_ROOT:-/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer}"
ODIA_OPT="${ODIA_OPT:-${ODIA_ROOT}/opt}"

export ODIA_ROOT ODIA_OPT
export ODIA_ENV="${ODIA_OPT}/env"
export ODIA_OPENMS="${ODIA_OPT}/openms-3.6.0"
export ODIA_MZPEAK="${ODIA_OPT}/mzpeak"

# Scratch is node-local NVMe; sources and build trees live there, not on ceph.
export ODIA_SCRATCH="${ODIA_SCRATCH:-/scratch/$(id -un)/odia}"

export PATH="${ODIA_ENV}/bin:${ODIA_OPT}/tools/bin:${PATH}"

export CC="${ODIA_ENV}/bin/x86_64-conda-linux-gnu-gcc"
export CXX="${ODIA_ENV}/bin/x86_64-conda-linux-gnu-g++"

export CMAKE_PREFIX_PATH="${ODIA_MZPEAK}:${ODIA_OPENMS}:${ODIA_ENV}${CMAKE_PREFIX_PATH:+:${CMAKE_PREFIX_PATH}}"
export PKG_CONFIG_PATH="${ODIA_MZPEAK}/lib/pkgconfig:${ODIA_ENV}/lib/pkgconfig${PKG_CONFIG_PATH:+:${PKG_CONFIG_PATH}}"
export LD_LIBRARY_PATH="${ODIA_MZPEAK}/lib:${ODIA_OPENMS}/lib:${ODIA_ENV}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

# Boost has no pkg-config file; meson finds it through BOOST_ROOT.
export BOOST_ROOT="${ODIA_ENV}"
