# Sorting our own output must change nothing, and every numeric column must
# survive the round-trip. The previous round-trip test compared three counts
# scraped from stdout and called that "lossless".
execute_process(COMMAND ${TOOL} -tr ${INPUT} -sort_library -stop_after library
                        -out_lib ${WORKDIR}/sort1.tsv
                RESULT_VARIABLE rc ERROR_VARIABLE err OUTPUT_QUIET)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "first pass failed:\n${err}")
endif()

execute_process(COMMAND ${TOOL} -tr ${WORKDIR}/sort1.tsv -sort_library -stop_after library
                        -out_lib ${WORKDIR}/sort2.tsv
                RESULT_VARIABLE rc ERROR_VARIABLE err OUTPUT_QUIET)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "second pass failed:\n${err}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files
                        ${WORKDIR}/sort1.tsv ${WORKDIR}/sort2.tsv
                RESULT_VARIABLE differ)
if(NOT differ EQUAL 0)
  message(FATAL_ERROR "sorting is not a fixed point: re-sorting our own output "
                      "changed it, so the order depends on the input's row order "
                      "rather than on the data")
endif()
message(STATUS "sort is a fixed point and the round-trip is byte-identical")
