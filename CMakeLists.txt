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

# Arrow/Parquet are a direct dependency, not just a transitive one: reaching
# the per-spectrum metadata (retention time, isolation windows) means going
# through mzpeak's low-level interface, which hands back Arrow objects.
# Version 23 is the floor OpenMS 3.6 and mzpeak both impose.
# The adopted scoring headers (include/odia/scoring) parallelise their tree and
# network fits with OpenMP. Optional: without it they still compile and run,
# single-threaded, because the pragmas are simply ignored -- so a missing
# OpenMP is a slower build, not a broken one.
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP: enabled for the scoring headers")
else()
  message(STATUS "OpenMP: not found; scoring will fit single-threaded")
endif()

# ONNX Runtime drives the PeptDeep models. OpenMS's own binding cannot be used:
# it rejects modified peptides and its headers are never installed.
find_path(ONNXRuntime_INCLUDE_DIR NAMES onnxruntime_cxx_api.h
          PATH_SUFFIXES include onnxruntime/core/session)
find_library(ONNXRuntime_LIBRARY NAMES onnxruntime PATH_SUFFIXES lib)
if(NOT ONNXRuntime_INCLUDE_DIR OR NOT ONNXRuntime_LIBRARY)
  message(FATAL_ERROR "ONNX Runtime not found; scripts/bootstrap_env.sh provides it")
endif()
message(STATUS "ONNX Runtime: ${ONNXRuntime_LIBRARY}")

# libzip reads the .oswpq container. It is already in the environment as a
# dependency of mzpeak; ODIA uses it directly because a bundle must be opened
# without a patched OpenMS.
find_package(libzip CONFIG QUIET)
if(NOT libzip_FOUND)
  find_path(LIBZIP_INCLUDE_DIR NAMES zip.h)
  find_library(LIBZIP_LIBRARY NAMES zip)
  if(NOT LIBZIP_INCLUDE_DIR OR NOT LIBZIP_LIBRARY)
    message(FATAL_ERROR "libzip not found; scripts/bootstrap_libzip.sh provides it")
  endif()
  add_library(libzip::zip UNKNOWN IMPORTED)
  set_target_properties(libzip::zip PROPERTIES
    IMPORTED_LOCATION "${LIBZIP_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${LIBZIP_INCLUDE_DIR}")
endif()

# SQLite reads the vendor instrument method embedded in an .mzpeak. The stated
# diaPASEF mobility bands live there and nowhere else in the converted file --
# see include/odia/VendorDiaWindows.h. Already in the environment; OpenMS uses
# it for .osw/.pqp, ODIA links it directly rather than through a patched OpenMS.
find_path(SQLite3_INCLUDE_DIR NAMES sqlite3.h)
find_library(SQLite3_LIBRARY NAMES sqlite3)
if(NOT SQLite3_INCLUDE_DIR OR NOT SQLite3_LIBRARY)
  message(FATAL_ERROR "sqlite3 not found; scripts/bootstrap_env.sh provides it")
endif()
message(STATUS "SQLite3: ${SQLite3_LIBRARY}")

find_package(Arrow CONFIG REQUIRED)
find_package(Parquet CONFIG REQUIRED)
if(Arrow_VERSION VERSION_LESS 23)
  message(FATAL_ERROR "Apache Arrow >= 23 required (found ${Arrow_VERSION}).")
endif()

if(TARGET Arrow::arrow_shared)
  set(ODIA_ARROW_TARGET Arrow::arrow_shared)
else()
  set(ODIA_ARROW_TARGET Arrow::arrow_static)
endif()
if(TARGET Parquet::parquet_shared)
  set(ODIA_PARQUET_TARGET Parquet::parquet_shared)
else()
  set(ODIA_PARQUET_TARGET Parquet::parquet_static)
endif()

message(STATUS "OpenMS: ${OpenMS_DIR}")
message(STATUS "mzPeak: ${MzPeak_LIBRARY}")
message(STATUS "Arrow:  ${Arrow_VERSION}")

# ---------------------------------------------------------------------------
# Shared code lives in a library so tools stay thin and are testable.
# ---------------------------------------------------------------------------
set(ODIA_LIB_SOURCES
  src/library/Library.cpp
  src/library/DIANNLibraryFile.cpp
  src/library/LibraryGenerator.cpp
  src/library/OSWPQLibraryFile.cpp
  src/io/ArrowColumn.cpp
  src/io/ResourceProbe.cpp
  src/io/MiniJson.cpp
  src/io/MobilityBands.cpp
  src/io/MzPeakSource.cpp
  src/io/VendorDiaWindows.cpp
  src/extract/ChromatogramExtractor.cpp
  src/extract/Ms1Traces.cpp
  src/extract/PrecursorPrefilter.cpp
  src/score/OpenSwathPicker.cpp
  src/score/PercolatorEngine.cpp
  src/calibrate/MassRecalibration.cpp
  src/score/PeakGroupScorer.cpp
  src/calibrate/RtCalibration.cpp
  src/calibrate/RtRefiner.cpp
  src/calibrate/MassCalibration.cpp
  src/calibrate/MobilityCalibration.cpp
  src/io/ChromatogramTsv.cpp
  src/io/FragvecTsv.cpp
  src/io/ZipArchive.cpp
  src/predict/PeptDeepEncoder.cpp
  src/predict/PeptDeepPredictor.cpp
)

if(ODIA_LIB_SOURCES)
  add_library(odia_core STATIC ${ODIA_LIB_SOURCES})
  target_include_directories(odia_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
  target_link_libraries(odia_core PUBLIC OpenMS MzPeak::MzPeak
                                         ${ODIA_ARROW_TARGET} ${ODIA_PARQUET_TARGET}
                                         ${ONNXRuntime_LIBRARY} libzip::zip
                                         ${SQLite3_LIBRARY})
  if(OpenMP_CXX_FOUND)
    target_link_libraries(odia_core PUBLIC OpenMP::OpenMP_CXX)
  endif()
  target_include_directories(odia_core PUBLIC ${ONNXRuntime_INCLUDE_DIR} ${SQLite3_INCLUDE_DIR})
  set_target_properties(odia_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

# ---------------------------------------------------------------------------
# TOPP-compatible tools. Each is a TOPPBase subclass, so it speaks the usual
# TOPP command line (-ini, -write_ini, -threads, ...) and emits CTD, but is
# built and shipped outside the OpenMS source tree.
# ---------------------------------------------------------------------------
set(ODIA_TOOLS
  DIALibraryGenerator
  ODIAInfo
  OpenDIAlyzer
)

# Where the shipped data tables live, so the tool finds them without a path
# baked in at an absolute location.
foreach(tool ${ODIA_TOOLS})
  add_executable(${tool} "src/${tool}.cpp")
  target_compile_definitions(${tool} PRIVATE
    ODIA_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
  target_link_libraries(${tool} PRIVATE OpenMS MzPeak::MzPeak
                                        ${ODIA_ARROW_TARGET} ${ODIA_PARQUET_TARGET})
  if(TARGET odia_core)
    target_link_libraries(${tool} PRIVATE odia_core)
  endif()
endforeach()

install(TARGETS ${ODIA_TOOLS} RUNTIME DESTINATION bin)

# Built only for the test suite: it dumps ODIA's encoding so it can be compared
# against the independent Python reference.
foreach(helper odia_encode_dump odia_group_dump odia_mod_symbols odia_predict_rt
               odia_predict_ms2 odia_predict_ccs odia_oswpq_dump odia_run_info odia_extract odia_irt_calibration
               odia_chunked_column_test odia_ms2_sessions odia_scoring_smoke odia_picker_cases
               odia_peak_bounds
               odia_interval_split
               odia_entrapment
               odia_rt_calibration odia_mass_calibration odia_mobility_calibration
               odia_extract_axis
               odia_extract_cases odia_mobility_bands odia_tsv_writers
               odia_prefilter_depth odia_ms1_probe odia_ms1_coelution odia_mz_residuals odia_mass_recalibration
               odia_mass_width odia_rt_refiner odia_mass_anchors odia_prefilter
               odia_empty_trace odia_robust_line odia_library_subset odia_pep_monotone
               odia_anchor_stop)
  add_executable(${helper} test/tools/${helper}.cpp)
  target_link_libraries(${helper} PRIVATE OpenMS odia_core
                        ${ODIA_ARROW_TARGET} ${ODIA_PARQUET_TARGET})
endforeach()

enable_testing()
add_subdirectory(test)
