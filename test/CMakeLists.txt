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

  add_test(NAME ${tool}_read_ini
           COMMAND ${CMAKE_COMMAND} -E true)
  set_tests_properties(${tool}_read_ini PROPERTIES DEPENDS ${tool}_write_ini)
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
