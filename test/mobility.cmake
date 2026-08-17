# Generates a small library, then round-trips it through the TSV writer and the
# reader, and asserts IM and CCS survive in the right units at every step.
# Expects: GEN TOOL PYTHON CHECKER FASTA WORKDIR RT_MODEL MS2_MODEL CCS_MODEL
set(cfg "${WORKDIR}/mobility.json")
file(WRITE "${cfg}"
  "{ \"rt_model\":\"${RT_MODEL}\", \"ms2_model\":\"${MS2_MODEL}\","
  " \"ccs_model\":\"${CCS_MODEL}\", \"decoys\":\"none\" }")
execute_process(COMMAND ${GEN} -in ${FASTA} -config ${cfg}
                        -out ${WORKDIR}/mobility.parquet
                RESULT_VARIABLE rc OUTPUT_VARIABLE o ERROR_VARIABLE e)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "DIALibraryGenerator failed: ${rc}\n${o}\n${e}")
endif()
# through the TSV interchange format and back -- this is the path DIA-NN uses
# -decoys none: the reader now completes a partial decoy set, so without this
# the round trip would legitimately gain decoys and the comparison would be
# against a different library.
execute_process(COMMAND ${TOOL} -tr ${WORKDIR}/mobility.parquet -stop_after library
                        -decoys none -out_lib ${WORKDIR}/mobility.tsv
                RESULT_VARIABLE rc OUTPUT_VARIABLE o ERROR_VARIABLE e)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "TSV export failed: ${rc}\n${o}\n${e}")
endif()
execute_process(COMMAND ${TOOL} -tr ${WORKDIR}/mobility.tsv -stop_after library
                        -decoys none -out_lib ${WORKDIR}/mobility.rt.parquet
                RESULT_VARIABLE rc OUTPUT_VARIABLE o ERROR_VARIABLE e)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "TSV re-read failed: ${rc}\n${o}\n${e}")
endif()
execute_process(COMMAND ${PYTHON} ${CHECKER}
                        ${WORKDIR}/mobility.parquet ${WORKDIR}/mobility.rt.parquet
                RESULT_VARIABLE rc OUTPUT_VARIABLE o ERROR_VARIABLE e)
message(STATUS "${o}${e}")
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "mobility check failed")
endif()
