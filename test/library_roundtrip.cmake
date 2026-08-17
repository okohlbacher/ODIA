# Assert that the Parquet and TSV readers agree, and that a TSV round-trip
# preserves precursor count, transition count and the interned string set.
function(summarise input outvar)
  # -decoys none: this test measures READER/WRITER fidelity. With decoy
  # completion active the fixture's 5-decoy set is discarded and regenerated on
  # load, and dropDecoys does not release interned strings, so the arena count
  # is inflated on the first pass and not on the second -- a 9-string
  # "difference" that is arena residue, not lost content.
  execute_process(COMMAND ${TOOL} -tr ${input} -stop_after library -decoys none
                  OUTPUT_VARIABLE out ERROR_VARIABLE err RESULT_VARIABLE rc)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "OpenDIAlyzer failed on ${input}:\n${err}")
  endif()
  string(REGEX MATCH "precursors:[ ]+([0-9]+)" _p "${out}${err}")
  set(p ${CMAKE_MATCH_1})
  string(REGEX MATCH "transitions:[ ]+([0-9]+)" _t "${out}${err}")
  set(t ${CMAKE_MATCH_1})
  string(REGEX MATCH "distinct strings:[ ]+([0-9]+)" _s "${out}${err}")
  set(s ${CMAKE_MATCH_1})
  set(${outvar} "${p}/${t}/${s}" PARENT_SCOPE)
endfunction()

summarise("${FIXTURES}/diann_library_small.parquet" from_parquet)
summarise("${FIXTURES}/diann_library_small.tsv" from_tsv)
if(NOT from_parquet STREQUAL from_tsv)
  message(FATAL_ERROR "readers disagree: parquet=${from_parquet} tsv=${from_tsv}")
endif()

execute_process(COMMAND ${TOOL} -tr ${FIXTURES}/diann_library_small.parquet
                        -stop_after library -decoys none
                        -out_lib ${WORKDIR}/roundtrip.tsv
                RESULT_VARIABLE rc ERROR_VARIABLE err OUTPUT_QUIET)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "writing the library failed:\n${err}")
endif()

summarise("${WORKDIR}/roundtrip.tsv" after_roundtrip)
if(NOT after_roundtrip STREQUAL from_parquet)
  message(FATAL_ERROR "round-trip lost content: before=${from_parquet} after=${after_roundtrip}")
endif()

message(STATUS "readers agree and round-trip is lossless: ${from_parquet}")
