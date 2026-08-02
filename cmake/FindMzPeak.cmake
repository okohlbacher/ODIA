# Locate the OpenMS/mzpeak C++ library.
#
# mzpeak is built with meson and ships neither a CMake package config nor a
# pkg-config file, so we look for the header and library directly.
#
# Set MZPEAK_ROOT (or put the install prefix on CMAKE_PREFIX_PATH) to point
# this at a specific build; scripts/env.sh does the latter.
#
# Defines the imported target MzPeak::MzPeak, plus MzPeak_FOUND,
# MzPeak_INCLUDE_DIR and MzPeak_LIBRARY.

find_path(MzPeak_INCLUDE_DIR
  NAMES mzpeak.h
  HINTS ${MZPEAK_ROOT} ENV MZPEAK_ROOT
  PATH_SUFFIXES include
  DOC "mzpeak include directory (the one containing mzpeak.h)"
)

find_library(MzPeak_LIBRARY
  NAMES mzpeak
  HINTS ${MZPEAK_ROOT} ENV MZPEAK_ROOT
  PATH_SUFFIXES lib lib64
  DOC "mzpeak library"
)

# The nested headers (mzpeak/io/archive.h and friends) must be present for the
# umbrella header to compile. mzpeak's meson install flattens them, so guard
# against being pointed at a broken install prefix -- the error is otherwise a
# wall of missing-include noise. scripts/build_mzpeak.sh installs a corrected
# layout.
if(MzPeak_INCLUDE_DIR AND NOT EXISTS "${MzPeak_INCLUDE_DIR}/mzpeak/io/archive.h")
  message(FATAL_ERROR
    "Found mzpeak.h in ${MzPeak_INCLUDE_DIR}, but the nested headers are "
    "missing (no mzpeak/io/archive.h). This is the flattened layout produced "
    "by mzpeak's own 'meson install'; rebuild it with scripts/build_mzpeak.sh, "
    "which installs a usable header tree.")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MzPeak
  REQUIRED_VARS MzPeak_LIBRARY MzPeak_INCLUDE_DIR
)

if(MzPeak_FOUND AND NOT TARGET MzPeak::MzPeak)
  add_library(MzPeak::MzPeak UNKNOWN IMPORTED)
  set_target_properties(MzPeak::MzPeak PROPERTIES
    IMPORTED_LOCATION "${MzPeak_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${MzPeak_INCLUDE_DIR}"
    # mzpeak's headers are C++23 (std::print, std::ranges::to,
    # std::move_only_function); anything including them needs the same.
    INTERFACE_COMPILE_FEATURES cxx_std_23
  )
endif()

mark_as_advanced(MzPeak_INCLUDE_DIR MzPeak_LIBRARY)
