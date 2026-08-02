cmake_minimum_required(VERSION 3.21 FATAL_ERROR)

project(ODIA
  VERSION 0.1.0
  DESCRIPTION "OpenDIAlyzer -- DIA extraction built on OpenMS and mzPeak"
  LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# mzpeak's headers require C++23 (std::print, std::ranges::to,
# std::move_only_function), which is the binding constraint for the whole
# project. GCC 14 or newer is needed; see scripts/bootstrap_env.sh.
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ODIA depends on OpenMS but never modifies it: OpenMS is consumed as an
# installed, read-only package, exactly like any other external library.
find_package(OpenMS REQUIRED)
find_package(MzPeak REQUIRED)

message(STATUS "OpenMS: ${OpenMS_DIR}")
message(STATUS "mzPeak: ${MzPeak_LIBRARY}")

# ---------------------------------------------------------------------------
# Shared code lives in a library so tools stay thin and are testable.
# ---------------------------------------------------------------------------
set(ODIA_LIB_SOURCES)

if(ODIA_LIB_SOURCES)
  add_library(odia_core STATIC ${ODIA_LIB_SOURCES})
  target_include_directories(odia_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
  target_link_libraries(odia_core PUBLIC OpenMS MzPeak::MzPeak)
  set_target_properties(odia_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

# ---------------------------------------------------------------------------
# TOPP-compatible tools. Each is a TOPPBase subclass, so it speaks the usual
# TOPP command line (-ini, -write_ini, -threads, ...) and emits CTD, but is
# built and shipped outside the OpenMS source tree.
# ---------------------------------------------------------------------------
set(ODIA_TOOLS
  ODIAInfo
)

foreach(tool ${ODIA_TOOLS})
  add_executable(${tool} "src/${tool}.cpp")
  target_link_libraries(${tool} PRIVATE OpenMS MzPeak::MzPeak)
  if(TARGET odia_core)
    target_link_libraries(${tool} PRIVATE odia_core)
  endif()
endforeach()

install(TARGETS ${ODIA_TOOLS} RUNTIME DESTINATION bin)

enable_testing()
add_subdirectory(test)
