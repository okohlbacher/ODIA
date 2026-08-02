# The tool must reject INPUT, and for the stated reason.
execute_process(COMMAND ${TOOL} -tr ${INPUT} -stop_after library
                RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(rc EQUAL 0)
  message(FATAL_ERROR "expected ${INPUT} to be rejected, but the tool succeeded")
endif()
string(FIND "${out}${err}" "${EXPECT}" found)
if(found EQUAL -1)
  message(FATAL_ERROR "rejected, but not for the stated reason.\n"
                      "  expected to see: ${EXPECT}\n  got: ${out}${err}")
endif()
message(STATUS "rejected for the right reason: ${EXPECT}")
