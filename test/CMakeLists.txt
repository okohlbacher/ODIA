# Tests run against the example data when it is available. ODIA_DATA_DIR
# defaults to the read-only example directory on this cluster; point it
# elsewhere to run the same tests on another machine.
set(ODIA_DATA_DIR "/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/data"
    CACHE PATH "Directory holding the example DIA runs")

# Every TOPP tool must answer --help and round-trip its own INI file. This
# catches parameter-registration mistakes without needing any data.
foreach(tool ${ODIA_TOOLS})
  add_test(NAME ${tool}_help COMMAND ${tool} --help)

  add_test(NAME ${tool}_write_ini
           COMMAND ${tool} -write_ini "${CMAKE_CURRENT_BINARY_DIR}/${tool}.ini")

  # (An INI round-trip test belongs here; the previous placeholder ran
  # `cmake -E true` and asserted nothing at all.)
endforeach()

if(EXISTS "${ODIA_DATA_DIR}/12_80.mzpeak")
  # 12_80 is the smallest example run, so it is the one used for tests that
  # must stay fast.
  add_test(NAME ODIAInfo_mzpeak
           COMMAND ODIAInfo -in "${ODIA_DATA_DIR}/12_80.mzpeak")
  set_tests_properties(ODIAInfo_mzpeak PROPERTIES
    PASS_REGULAR_EXPRESSION "format:[ ]+mzPeak")
endif()

if(EXISTS "${ODIA_DATA_DIR}/12_80.mzML")
  add_test(NAME ODIAInfo_mzml
           COMMAND ODIAInfo -in "${ODIA_DATA_DIR}/12_80.mzML" -windows)
  set_tests_properties(ODIAInfo_mzml PROPERTIES
    PASS_REGULAR_EXPRESSION "distinct isolation windows")
endif()

# --------------------------------------------------------------------- library
# Fixtures are derived, not committed; scripts/make_fixtures.sh builds them onto
# node-local scratch. Skip these tests where they are absent rather than fail.
set(ODIA_FIXTURE_DIR "$ENV{ODIA_SCRATCH}/fixtures"
    CACHE PATH "Directory holding the assay-library fixtures")

if(EXISTS "${ODIA_FIXTURE_DIR}/diann_library_small.parquet")
  add_test(NAME library_parquet
           COMMAND OpenDIAlyzer -tr "${ODIA_FIXTURE_DIR}/diann_library_small.parquet"
                                -stop_after library)
  set_tests_properties(library_parquet PROPERTIES
    PASS_REGULAR_EXPRESSION "precursors:[ ]+500")
endif()

if(EXISTS "${ODIA_FIXTURE_DIR}/diann_library_small.tsv")
  add_test(NAME library_tsv
           COMMAND OpenDIAlyzer -tr "${ODIA_FIXTURE_DIR}/diann_library_small.tsv"
                                -stop_after library)
  set_tests_properties(library_tsv PROPERTIES
    PASS_REGULAR_EXPRESSION "precursors:[ ]+500")

  # The two readers must agree, and a TSV round-trip must be lossless in the
  # quantities we keep. library_roundtrip.cmake asserts both.
  add_test(NAME library_roundtrip
           COMMAND ${CMAKE_COMMAND}
                   -DTOOL=$<TARGET_FILE:OpenDIAlyzer>
                   -DFIXTURES=${ODIA_FIXTURE_DIR}
                   -DWORKDIR=${CMAKE_CURRENT_BINARY_DIR}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/library_roundtrip.cmake)
endif()

# ----------------------------------------------------- adversarial fixtures
# Every critical defect found so far lived in a path the DIA-NN fixture cannot
# reach. These are its complement: neutral losses, N-terminal and multiple
# modifications, decoy-sequence collisions, unusual Parquet types, nulls,
# malformed numbers and out-of-order rows.
#
# check_invariants.py recomputes fragment masses from its own residue table, not
# from OpenMS, so it cannot agree with ODIA by sharing its arithmetic.
find_program(ODIA_PYTHON python3 HINTS "$ENV{ODIA_ENV}/bin" NO_DEFAULT_PATH)
if(NOT ODIA_PYTHON)
  find_program(ODIA_PYTHON python3)
endif()

set(ODIA_ADVERSARIAL
  adv_losses adv_mods adv_collisions adv_fragments adv_target_decoy)

foreach(fx ${ODIA_ADVERSARIAL})
  if(EXISTS "${ODIA_FIXTURE_DIR}/${fx}.tsv" AND ODIA_PYTHON)
    add_test(NAME adversarial_${fx}
             COMMAND ${CMAKE_COMMAND}
                     -DTOOL=$<TARGET_FILE:OpenDIAlyzer>
                     -DPYTHON=${ODIA_PYTHON}
                     -DCHECKER=${CMAKE_CURRENT_SOURCE_DIR}/check_invariants.py
                     -DINPUT=${ODIA_FIXTURE_DIR}/${fx}.tsv
                     -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/${fx}.out.tsv
                     -DCHECK_DECOYS=--decoy-table
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/run_invariants.cmake)
  endif()
endforeach()

# The same invariants must hold on the real library -- targets only, since
# DIA-NN's own decoy rows were mutated by DIA-NN's position rule, not ours, and
# cannot be reproduced from the substitution table alone.
if(EXISTS "${ODIA_FIXTURE_DIR}/diann_library_small.tsv" AND ODIA_PYTHON)
  add_test(NAME adversarial_real_library
           COMMAND ${CMAKE_COMMAND}
                   -DTOOL=$<TARGET_FILE:OpenDIAlyzer>
                   -DPYTHON=${ODIA_PYTHON}
                   -DCHECKER=${CMAKE_CURRENT_SOURCE_DIR}/check_invariants.py
                   -DINPUT=${ODIA_FIXTURE_DIR}/diann_library_small.tsv
                   -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/real.out.tsv
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/run_invariants.cmake)
endif()

# Rows out of precursor order must be refused, not silently split.
if(EXISTS "${ODIA_FIXTURE_DIR}/adv_interleaved.tsv")
  # WILL_FAIL alone passes on any non-zero exit -- it passed even with the
  # contiguity check deleted outright -- and it cannot be combined with
  # PASS_REGULAR_EXPRESSION, which it inverts. Assert both from a script.
  add_test(NAME adversarial_interleaved_rejected
           COMMAND ${CMAKE_COMMAND}
                   -DTOOL=$<TARGET_FILE:OpenDIAlyzer>
                   -DINPUT=${ODIA_FIXTURE_DIR}/adv_interleaved.tsv
                   -DEXPECT=not grouped by precursor
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/expect_rejection.cmake)
endif()

# Unusual Parquet types must decode, including a bool decoy flag.
if(EXISTS "${ODIA_FIXTURE_DIR}/adv_types.parquet")
  add_test(NAME adversarial_parquet_types
           COMMAND OpenDIAlyzer -tr "${ODIA_FIXTURE_DIR}/adv_types.parquet"
                                -stop_after library)
  set_tests_properties(adversarial_parquet_types PROPERTIES
    PASS_REGULAR_EXPRESSION "1 target / 1 decoy")
endif()

# A library ODIA generated itself, decoys included. This is the only test where
# the decoy substitution table is ours end to end, so the mutation check is
# meaningful.
# adv_proteins.fasta comes from make_fixtures.sh, so this test travels. It
# previously pointed at a file in one user's scratch that no script produces,
# which made it the only test reaching -fasta *and* the only one that silently
# vanished on any other machine -- taking the sole guard for the substitution
# table, the fragment ordinals and the RT write path with it.
if(EXISTS "${ODIA_FIXTURE_DIR}/adv_proteins.fasta" AND ODIA_PYTHON)
  add_test(NAME adversarial_generated
           COMMAND ${CMAKE_COMMAND}
                   -DTOOL=$<TARGET_FILE:OpenDIAlyzer>
                   -DPYTHON=${ODIA_PYTHON}
                   -DCHECKER=${CMAKE_CURRENT_SOURCE_DIR}/check_invariants.py
                   -DFASTA=${ODIA_FIXTURE_DIR}/adv_proteins.fasta
                   -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/generated.out.tsv
                   -DCHECK_DECOYS=--decoy-table
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/run_invariants.cmake)
endif()

# Malformed numbers must be reported, not silently accepted.
if(EXISTS "${ODIA_FIXTURE_DIR}/adv_badnumbers.tsv")
  add_test(NAME adversarial_badnumbers
           COMMAND OpenDIAlyzer -tr "${ODIA_FIXTURE_DIR}/adv_badnumbers.tsv"
                                -stop_after library)
  set_tests_properties(adversarial_badnumbers PROPERTIES
    PASS_REGULAR_EXPRESSION "precursors have an unusable m/z")
endif()

# All-null columns must not silently become zeros without a word.
if(EXISTS "${ODIA_FIXTURE_DIR}/adv_nulls.parquet")
  add_test(NAME adversarial_nulls
           COMMAND OpenDIAlyzer -tr "${ODIA_FIXTURE_DIR}/adv_nulls.parquet"
                                -stop_after library)
  set_tests_properties(adversarial_nulls PROPERTIES
    PASS_REGULAR_EXPRESSION "precursors:")
endif()

# Sorting must be a property of the data: re-sorting our own output must be a
# fixed point. Nothing detected a tie-break keyed on input order before this.
if(EXISTS "${ODIA_FIXTURE_DIR}/diann_library_small.tsv")
  add_test(NAME sort_is_a_fixed_point
           COMMAND ${CMAKE_COMMAND}
                   -DTOOL=$<TARGET_FILE:OpenDIAlyzer>
                   -DINPUT=${ODIA_FIXTURE_DIR}/diann_library_small.tsv
                   -DWORKDIR=${CMAKE_CURRENT_BINARY_DIR}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/sort_fixed_point.cmake)
endif()

# m/z must survive load-store to the precision D6 claims. A mutation making
# toFixed truncate, or the writer drop two significant digits, passed every
# other test in the suite.
if(EXISTS "${ODIA_FIXTURE_DIR}/diann_library_small.tsv" AND ODIA_PYTHON)
  add_test(NAME mz_precision_round_trip
           COMMAND ${CMAKE_COMMAND}
                   -DTOOL=$<TARGET_FILE:OpenDIAlyzer>
                   -DPYTHON=${ODIA_PYTHON}
                   -DCHECKER=${CMAKE_CURRENT_SOURCE_DIR}/check_precision.py
                   -DINPUT=${ODIA_FIXTURE_DIR}/diann_library_small.tsv
                   -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/precision.out.tsv
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/run_precision.cmake)
endif()
