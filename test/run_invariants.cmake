# Load a library, write it back, and assert the invariants on what was written.
if(DEFINED FASTA)
  set(SOURCE_ARGS -fasta ${FASTA})
else()
  set(SOURCE_ARGS -tr ${INPUT})
endif()

execute_process(COMMAND ${TOOL} ${SOURCE_ARGS} -stop_after library -out_lib ${OUTPUT}
                RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "OpenDIAlyzer failed on ${INPUT}${FASTA}:\n${out}${err}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${PYTHON} ${CHECKER} ${OUTPUT} ${CHECK_DECOYS}
                RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
message(STATUS "${out}${err}")
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "invariant check failed for ${INPUT}")
endif()
