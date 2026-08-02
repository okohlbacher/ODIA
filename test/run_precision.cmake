execute_process(COMMAND ${TOOL} -tr ${INPUT} -stop_after library -out_lib ${OUTPUT}
                RESULT_VARIABLE rc ERROR_VARIABLE err OUTPUT_QUIET)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "OpenDIAlyzer failed:\n${err}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${PYTHON} ${CHECKER} ${INPUT} ${OUTPUT}
                RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
message(STATUS "${out}${err}")
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "m/z precision check failed")
endif()
