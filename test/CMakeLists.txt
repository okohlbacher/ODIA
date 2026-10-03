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

# Where the PeptDeep retention-time model lives. Declared once, used by every
# test that needs it, and every such test is guarded on its existence -- a test
# that silently vanishes elsewhere is bad, but one that fails elsewhere is worse.
set(ODIA_RT_MODEL_PATH
    "$ENV{ODIA_OPENMS}/share/OpenMS/models/peptdeep_rt_dynamic.onnx"
    CACHE FILEPATH "PeptDeep RT model")

# Derived from the RT path, NOT from $ENV{ODIA_OPENMS} again: the environment is
# read at configure time, so a `cmake .` in a shell that has not sourced
# scripts/env.sh silently bakes "/share/OpenMS/models/..." into every test that
# refers to it -- which is exactly what happened, and it fails at RUN time with
# a missing-file error that looks like a broken model rather than a broken path.
string(REPLACE "peptdeep_rt_dynamic.onnx" "peptdeep_ms2_dynamic.onnx"
       ODIA_MS2_MODEL_PATH "${ODIA_RT_MODEL_PATH}")
string(REPLACE "peptdeep_rt_dynamic.onnx" "peptdeep_ccs_dynamic.onnx"
       ODIA_CCS_MODEL_PATH "${ODIA_RT_MODEL_PATH}")

# --------------------------------------------------------------------- library
# Fixtures are derived, not committed; scripts/make_fixtures.sh builds them onto
# node-local scratch. Skip these tests where they are absent rather than fail.
set(ODIA_FIXTURE_DIR "$ENV{ODIA_SCRATCH}/fixtures"
    CACHE PATH "Directory holding the assay-library fixtures")

# Absent fixtures used to be SILENT. Every test below sits inside if(EXISTS ...),
# ODIA_SCRATCH is node-local, and the source library they are derived from lived
# on one node's /scratch -- so on the compute node where the work actually
# happens, ctest reported "100% tests passed, 0 tests failed out of 52" while
# TWENTY-FOUR tests were never configured, among them every library-reader,
# decoy-method and adversarial-input case. With the fixtures present the same
# tree reports 76, and all 76 pass: the code was fine, it simply was not checked.
#
# A skip that reports success is worse than a failure. Say so at configure time,
# loudly, and name the fix. scripts/make_fixtures.sh now defaults to a copy on
# shared storage so any node can derive them.
if(NOT EXISTS "${ODIA_FIXTURE_DIR}/diann_library_small.parquet")
  message(WARNING
    "Assay-library fixtures are ABSENT from ${ODIA_FIXTURE_DIR}.\n"
    "  ~24 tests -- every library-reader, decoy-method and adversarial-input "
    "case -- will be SILENTLY SKIPPED, and ctest will still report success.\n"
    "  Fix: bash scripts/make_fixtures.sh   (reads the shared source under "
    "shared/ref/fixture_src/, falls back to node-local /scratch)\n"
    "  Configure with -DODIA_REQUIRE_FIXTURES=ON to make this fatal instead.")
  if(ODIA_REQUIRE_FIXTURES)
    message(FATAL_ERROR "ODIA_REQUIRE_FIXTURES=ON and the fixtures are absent.")
  endif()
endif()

if(EXISTS "${ODIA_FIXTURE_DIR}/diann_library_small.parquet")
  # -decoys none keeps this a READER test. Decoy completion has its own test.
  add_test(NAME library_parquet
           COMMAND OpenDIAlyzer -tr "${ODIA_FIXTURE_DIR}/diann_library_small.parquet"
                                -stop_after library -decoys none)
  # 904, not 500: the fixture is a truncated real library, so its 5 decoys are
  # an artefact of where the cut fell. A decoy set that small cannot support an
  # FDR estimate, so the reader now discards and regenerates it (495 targets,
  # 409 decoys -- the rest collide with an existing target under mutation).
  set_tests_properties(library_parquet PROPERTIES
    PASS_REGULAR_EXPRESSION "precursors:[ ]+500")

  # Decoy completion, on its own. A library with a TOKEN decoy set -- which is
  # what DIA-NN's empirical library has, 365 for 37,193 targets -- must have it
  # replaced, not topped up: appendDecoys keys on the decoy it would generate,
  # so a foreign one is invisible and the target gets a second.
  add_test(NAME library_partial_decoys_replaced
           COMMAND OpenDIAlyzer -tr "${ODIA_FIXTURE_DIR}/diann_library_small.parquet"
                                -stop_after library)
  set_tests_properties(library_partial_decoys_replaced PROPERTIES
    PASS_REGULAR_EXPRESSION "495 target / 409 decoy")
endif()

if(EXISTS "${ODIA_FIXTURE_DIR}/diann_library_small.tsv")
  add_test(NAME library_tsv
           COMMAND OpenDIAlyzer -tr "${ODIA_FIXTURE_DIR}/diann_library_small.tsv"
                                -stop_after library -decoys none)
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
if(EXISTS "${ODIA_FIXTURE_DIR}/adv_proteins.fasta" AND ODIA_PYTHON
   AND EXISTS "${ODIA_RT_MODEL_PATH}")
  add_test(NAME adversarial_generated
           COMMAND ${CMAKE_COMMAND}
                   -DTOOL=$<TARGET_FILE:OpenDIAlyzer>
                   -DPYTHON=${ODIA_PYTHON}
                   -DCHECKER=${CMAKE_CURRENT_SOURCE_DIR}/check_invariants.py
                   -DFASTA=${ODIA_FIXTURE_DIR}/adv_proteins.fasta
                   -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/generated.out.tsv
                   -DCHECK_DECOYS=--decoy-table
                   -DREQUIRE_RT=--require-rt
                   -DRT_MODEL=${ODIA_RT_MODEL_PATH}
                   -DIRT_STANDARDS=${CMAKE_CURRENT_SOURCE_DIR}/../data/irt_standards.tsv
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/run_invariants.cmake)
endif()

# The free-cysteine RT correction is a table of measured constants (doc/30), and
# a table of constants is exactly the kind of thing that drifts unnoticed. This
# also guards the self-gating property that lets it default to on: an alkylated
# library must be bit-identical whether the flag is set or not.
if(EXISTS "${ODIA_FIXTURE_DIR}/adv_proteins.fasta" AND ODIA_PYTHON
   AND EXISTS "${ODIA_RT_MODEL_PATH}")
  add_test(NAME free_cysteine_correction
           COMMAND ${CMAKE_COMMAND}
                   -DTOOL=$<TARGET_FILE:DIALibraryGenerator>
                   -DPYTHON=${ODIA_PYTHON}
                   -DCHECKER=${CMAKE_CURRENT_SOURCE_DIR}/check_free_cysteine.py
                   -DFASTA=${ODIA_FIXTURE_DIR}/adv_proteins.fasta
                   -DWORKDIR=${CMAKE_CURRENT_BINARY_DIR}
                   -DRT_MODEL=${ODIA_RT_MODEL_PATH}
                   -DMS2_MODEL=${ODIA_MS2_MODEL_PATH}
                   -DCCS_MODEL=${ODIA_CCS_MODEL_PATH}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/free_cysteine.cmake)
endif()

# IM and CCS are different quantities in different units, and a library that
# swapped them would still look plausible. Guards presence, physical range,
# mutual consistency under Mason-Schamp, and survival of the TSV round trip --
# the interchange path DIA-NN reads, since it rejects our Parquet outright
# (DIA-NN 2.0 requires a 'Flags' column we do not write).
if(EXISTS "${ODIA_FIXTURE_DIR}/adv_proteins.fasta" AND ODIA_PYTHON
   AND EXISTS "${ODIA_RT_MODEL_PATH}")
  add_test(NAME mobility_units
           COMMAND ${CMAKE_COMMAND}
                   -DGEN=$<TARGET_FILE:DIALibraryGenerator>
                   -DTOOL=$<TARGET_FILE:OpenDIAlyzer>
                   -DPYTHON=${ODIA_PYTHON}
                   -DCHECKER=${CMAKE_CURRENT_SOURCE_DIR}/check_mobility.py
                   -DFASTA=${ODIA_FIXTURE_DIR}/adv_proteins.fasta
                   -DWORKDIR=${CMAKE_CURRENT_BINARY_DIR}
                   -DRT_MODEL=${ODIA_RT_MODEL_PATH}
                   -DMS2_MODEL=${ODIA_MS2_MODEL_PATH}
                   -DCCS_MODEL=${ODIA_CCS_MODEL_PATH}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/mobility.cmake)
endif()

# The decoy methods must actually be different methods. They were not: the
# generator hand-rolled its own name lookup and mapped every unknown name to
# Mutate, so three of the four produced identical libraries in silence.
if(EXISTS "${ODIA_FIXTURE_DIR}/adv_proteins.fasta" AND ODIA_PYTHON
   AND EXISTS "${ODIA_RT_MODEL_PATH}")
  add_test(NAME decoy_methods_distinct
           COMMAND ${CMAKE_COMMAND}
                   -DGEN=$<TARGET_FILE:DIALibraryGenerator>
                   -DPYTHON=${ODIA_PYTHON}
                   -DCHECKER=${CMAKE_CURRENT_SOURCE_DIR}/check_decoy_methods.py
                   -DFASTA=${ODIA_FIXTURE_DIR}/adv_proteins.fasta
                   -DWORKDIR=${CMAKE_CURRENT_BINARY_DIR}
                   -DRT_MODEL=${ODIA_RT_MODEL_PATH}
                   -DMS2_MODEL=${ODIA_MS2_MODEL_PATH}
                   -DCCS_MODEL=${ODIA_CCS_MODEL_PATH}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/decoy_methods.cmake)
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

# The PeptDeep reference encoder is the oracle ODIA's C++ encoder will be checked
# against, so it carries its own guard.
if(ODIA_PYTHON)
  add_test(NAME peptdeep_reference_encoding
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/check_peptdeep_encoding.py)
  # The C++ header and the Python data file are the same table for two
  # languages, and nothing compiled or tested the header before.
  add_test(NAME peptdeep_element_table
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/check_element_table.py)
endif()

# ODIA's C++ encoder against the independent Python reference. Written from the
# same spec but from different inputs -- OpenMS AASequence versus the raw string
# -- so agreement is evidence the spec was implemented, not that one copied the
# other.
if(ODIA_PYTHON)
  add_test(NAME peptdeep_encoder_matches_reference
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON}
                   ${CMAKE_CURRENT_SOURCE_DIR}/compare_encoders.py
                   $<TARGET_FILE:odia_encode_dump>)
endif()

# groupByLength and the homogeneity guard. Collapsing every peptide into one
# group, or removing the guard, previously passed the whole suite -- and mixing
# lengths changes a prediction by more than most modifications do.
add_test(NAME peptdeep_length_grouping
         COMMAND odia_group_dump PEPTIDEK PEPTIDER PEPTIDEKK PEPTIDEKKK PEPTIDEM)
set_tests_properties(peptdeep_length_grouping PROPERTIES
  PASS_REGULAR_EXPRESSION "groups: \\[0,1,4\\] \\[2\\] \\[3\\]\nmixed_batch_refused: yes")

# A modification with no elemental composition encodes to all zeros, which is
# indistinguishable from no modification. OpenMS accepts bare mass shifts.
add_test(NAME peptdeep_rejects_massshift_modification
         COMMAND odia_encode_dump "PEPTIDEC[999]")
set_tests_properties(peptdeep_rejects_massshift_modification PROPERTIES
  WILL_FAIL TRUE)

# Every element symbol in every modification OpenMS knows must either map into
# AlphaPeptDeep's 109-element list or be a deliberate '?' bucket entry.
#
# This is the only check that covers the isotope translation: the Python
# reference cannot parse a name like Label:13C(6)15N(2) -- its parser splits on
# the first ')' -- so deleting the rename table entirely passed the encoder
# comparison. The four symbols pinned here are halogen isotopes from PSI-MOD,
# which AlphaPeptDeep does not load at all, so there is no upstream behaviour to
# match; they are recorded so a new unmapped symbol fails instead of silently
# landing in '?'.
add_test(NAME peptdeep_modification_symbols COMMAND odia_mod_symbols)
set_tests_properties(peptdeep_modification_symbols PROPERTIES
  PASS_REGULAR_EXPRESSION "unmapped symbols: \\(35\\)Cl \\(37\\)Cl \\(79\\)Br \\(81\\)Br")

# Predictions, not just tensors. Encoded values can match while dtype, layout or
# ordering assumptions differ, and the sequence list mixes lengths so a broken
# grouping or a lost reordering shows up as wrong values rather than an error.
if(EXISTS "${ODIA_RT_MODEL_PATH}" AND ODIA_PYTHON)
  add_test(NAME peptdeep_predictions_match_reference
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON}
                   ${CMAKE_CURRENT_SOURCE_DIR}/compare_predictions.py
                   $<TARGET_FILE:odia_predict_rt> "${ODIA_RT_MODEL_PATH}")
endif()

# MS2 intensities, covering the inputs the RT model does not have: charge, NCE
# and the rank-1 instrument index. Each is varied across the cases, because a
# defect that ignores one of them passes any fixture that holds it constant.
set(ODIA_MS2_MODEL_PATH
    "$ENV{ODIA_OPENMS}/share/OpenMS/models/peptdeep_ms2_dynamic.onnx"
    CACHE FILEPATH "PeptDeep MS2 model")
if(EXISTS "${ODIA_MS2_MODEL_PATH}" AND ODIA_PYTHON)
  add_test(NAME peptdeep_ms2_matches_reference
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON}
                   ${CMAKE_CURRENT_SOURCE_DIR}/compare_ms2.py
                   $<TARGET_FILE:odia_predict_ms2> "${ODIA_MS2_MODEL_PATH}")
endif()

# ---------------------------------------------------------------------------
# .oswpq, OpenMS's Parquet-in-ZIP assay-library bundle.
#
# The fixtures are generated rather than committed because each one exists to
# encode a specific hazard, and a binary blob cannot say which. They are built
# with pyarrow and Python's zipfile -- independent implementations of both
# formats -- so agreement with ODIA is evidence, not a tautology.
# ---------------------------------------------------------------------------
# Multi-chunk columns cannot be produced through a Parquet fixture -- Arrow
# concatenates row groups -- so the cursor is exercised against chunk layouts
# built directly. Without this, a "chunk(0) is the whole column" defect, which
# is the exact failure this class exists to prevent, passes every other test.
add_test(NAME chunked_column COMMAND odia_chunked_column_test)

if(ODIA_PYTHON)
  set(ODIA_OSWPQ_FIXTURES "${CMAKE_CURRENT_BINARY_DIR}/oswpq")

  add_test(NAME oswpq_fixtures
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON}
                   ${CMAKE_CURRENT_SOURCE_DIR}/make_oswpq_fixtures.py
                   "${ODIA_OSWPQ_FIXTURES}")
  set_tests_properties(oswpq_fixtures PROPERTIES FIXTURES_SETUP oswpq)

  # The one bundle written by an OpenMS build nobody here controls. It is the
  # only evidence that this reads what upstream actually writes, as opposed to
  # what the fixtures and the specification agree it writes -- and it already
  # disagrees with the specification in three places (unannotated transitions,
  # a precursor with m/z 0, and a precursor with no fragments at all).
  set(ODIA_UPSTREAM_OSWPQ
      "$ENV{ODIA_OPENMS_SOURCE}/src/tests/class_tests/openms/data/OpenSwathWorkflow_tworuns_1_17.output.oswpq"
      CACHE FILEPATH "An .oswpq bundle written by upstream OpenMS")
  set(ODIA_UPSTREAM_ARG "")
  if(EXISTS "${ODIA_UPSTREAM_OSWPQ}")
    set(ODIA_UPSTREAM_ARG "${ODIA_UPSTREAM_OSWPQ}")
  else()
    message(STATUS "oswpq: no upstream bundle at ${ODIA_UPSTREAM_OSWPQ}; "
                   "testing against generated fixtures only")
  endif()

  add_test(NAME oswpq_matches_pyarrow
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON}
                   ${CMAKE_CURRENT_SOURCE_DIR}/check_oswpq.py
                   $<TARGET_FILE:odia_oswpq_dump> "${ODIA_OSWPQ_FIXTURES}"
                   ${ODIA_UPSTREAM_ARG})
  set_tests_properties(oswpq_matches_pyarrow PROPERTIES FIXTURES_REQUIRED oswpq)
endif()

# The MS2 meta-input constants, pinned. A comparison between ODIA and the
# oracle cannot see a constant that is wrong in both, and one of these -- the
# instrument order -- is documented incorrectly in OpenMS's own headers.
if(ODIA_PYTHON)
  add_test(NAME peptdeep_meta_inputs
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON}
                   ${CMAKE_CURRENT_SOURCE_DIR}/check_meta_inputs.py)
endif()

# CCS, the third model kind. Its three inputs match neither the RT model's two
# nor the MS2 model's five, so the model guard rejected it outright until it
# grew a third arm.
set(ODIA_CCS_MODEL_PATH
    "$ENV{ODIA_OPENMS}/share/OpenMS/models/peptdeep_ccs_dynamic.onnx"
    CACHE FILEPATH "PeptDeep CCS model")
if(EXISTS "${ODIA_CCS_MODEL_PATH}" AND ODIA_PYTHON)
  add_test(NAME peptdeep_ccs_matches_reference
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON}
                   ${CMAKE_CURRENT_SOURCE_DIR}/compare_ccs.py
                   $<TARGET_FILE:odia_predict_ccs> "${ODIA_CCS_MODEL_PATH}")
endif()

# The iRT calibration, against OpenMS's own fixture rather than against a
# constant this project wrote down. The fixture carries both the raw ONNX
# prediction and the calibrated iRT for the 11 Biognosys standards, so it tests
# our inference and our line separately.
set(ODIA_IRT_FIXTURE
    "$ENV{ODIA_OPENMS_SOURCE}/src/tests/class_tests/openms/data/peptdeep_irt_peptides_predicted.csv"
    CACHE FILEPATH "OpenMS's predicted iRT standards fixture")
if(EXISTS "${ODIA_IRT_FIXTURE}" AND EXISTS "${ODIA_RT_MODEL_PATH}" AND ODIA_PYTHON)
  add_test(NAME irt_calibration_matches_openms
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 ${ODIA_PYTHON}
                   ${CMAKE_CURRENT_SOURCE_DIR}/check_irt.py
                   $<TARGET_FILE:odia_irt_calibration> "${ODIA_RT_MODEL_PATH}"
                   "${CMAKE_CURRENT_SOURCE_DIR}/../data/irt_standards.tsv"
                   "${ODIA_IRT_FIXTURE}")
else()
  message(STATUS "iRT calibration test skipped; no fixture at ${ODIA_IRT_FIXTURE}")
endif()

# A standards file whose peptides all predict the same retention time has no
# line to fit. Rejecting it matters because the alternative is a slope of
# infinity silently applied to every precursor in the library.
if(EXISTS "${ODIA_RT_MODEL_PATH}")
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/degenerate_irt.tsv"
       "# every peptide identical, so every prediction is identical\nPEPTIDEK\t0.0\nPEPTIDEK\t50.0\nPEPTIDEK\t100.0\n")
  add_test(NAME irt_calibration_rejects_degenerate
           COMMAND odia_irt_calibration "${ODIA_RT_MODEL_PATH}"
                   "${CMAKE_CURRENT_BINARY_DIR}/degenerate_irt.tsv")
  set_tests_properties(irt_calibration_rejects_degenerate PROPERTIES WILL_FAIL TRUE)

  # Two standards is not a calibration, it is a line through two points with no
  # residual to check. Three is the documented floor.
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/too_few_irt.tsv"
       "LGGNEQVTR\t-24.92\nLFLQFGAQGSPFLK\t100.0\n")
  add_test(NAME irt_calibration_rejects_too_few
           COMMAND odia_irt_calibration "${ODIA_RT_MODEL_PATH}"
                   "${CMAKE_CURRENT_BINARY_DIR}/too_few_irt.tsv")
  set_tests_properties(irt_calibration_rejects_too_few PROPERTIES WILL_FAIL TRUE)
endif()

# Standards that parse but cannot be encoded. One bad standard must cost
# precision, not the calibration; three bad ones out of four must fail loudly
# rather than fit a line through the survivor.
if(EXISTS "${ODIA_RT_MODEL_PATH}")
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/unencodable_irt.tsv"
       "LGGNEQVTR\t-24.92\nPEPTIDEC[999]K\t20.0\nSAMPLEC[999]R\t60.0\nAAAAC[999]K\t90.0\n")
  add_test(NAME irt_calibration_rejects_unencodable
           COMMAND odia_irt_calibration "${ODIA_RT_MODEL_PATH}"
                   "${CMAKE_CURRENT_BINARY_DIR}/unencodable_irt.tsv")
  set_tests_properties(irt_calibration_rejects_unencodable PROPERTIES WILL_FAIL TRUE)
endif()

# The fine-tuning modification parser. It had four defects and no test; the
# worst passed a named modification through into the sequence column, which
# would have been fine-tuned on as though it were a peptide.
if(ODIA_FINETUNE_PYTHON OR EXISTS "$ENV{ODIA_FINETUNE_PYTHON}")
  set(ODIA_FT_PY "$ENV{ODIA_FINETUNE_PYTHON}")
else()
  set(ODIA_FT_PY "/scratch/kohlbach/odia/rtfinetune/env/bin/python")
endif()
if(EXISTS "${ODIA_FT_PY}")
  add_test(NAME finetune_parser
           COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 "${ODIA_FT_PY}"
                   ${CMAKE_CURRENT_SOURCE_DIR}/check_finetune_parser.py)
else()
  message(STATUS "finetune parser test skipped; no python with peptdeep at ${ODIA_FT_PY}")
endif()

# Spreading MS2 batches across sessions must be bit-for-bit identical to
# running them on one. The chunk list is built in serial order so that it is,
# but that is a claim about the code, and the model's arithmetic is what
# actually has to agree -- a 4.7x speedup that quietly perturbs intensities
# would change which fragments rank into the top twelve.
if(EXISTS "${ODIA_MS2_MODEL_PATH}" AND EXISTS "${ODIA_FIXTURE_DIR}/diann_library_small.tsv")
  add_test(NAME ms2_sessions_are_deterministic
           COMMAND odia_ms2_sessions "${ODIA_MS2_MODEL_PATH}"
                   "${ODIA_FIXTURE_DIR}/diann_library_small.tsv" 400)
endif()

# The adopted scoring stack (include/odia/scoring, from OpenDIAlyzer
# odia-v0.3.0). Guards the failure mode that a well-formed-but-unfitted scorer
# would otherwise pass silently -- see the tool's own comment.
add_test(NAME scoring_stack_adopted COMMAND odia_scoring_smoke)

# Does the picker find the peptide, or the loudest thing in the window?
#
# Added because swapping the entire detection algorithm -- amplitude
# local-maxima to DIA-NN's pairwise-correlation search -- left all 59 tests
# passing. Nothing was sensitive to WHICH peaks were found, so a change that
# moved rank-1 accuracy from 40.7% to 75.7% on real data was invisible here.
# The case asserts BOTH pickers' behaviour, so it fails if either the old
# detector stops taking the spike or the new one stops rejecting it.
add_test(NAME picker_prefers_coelution_over_amplitude COMMAND odia_picker_cases)

# Does a peak group cover the peak, or the whole extraction window?
#
# It covered the window: 77.4% of 400k production peak groups spanned more than
# 80% of the extraction window, against a fitted 3.5 s FWHM -- about 2.5 cycles.
# Nothing in this suite asserted candidate GEOMETRY, only that the plumbing
# around a candidate held, so the defect was invisible for as long as it existed.
add_test(NAME peak_bounds_cover_the_peak_not_the_window COMMAND odia_peak_bounds)
add_test(NAME scoring_and_quantification_intervals_are_independent COMMAND odia_interval_split)

# Does a reported 1% FDR mean 1% of the identifications are wrong?
#
# doc/07 lists entrapment as ordering item 4 and it had never been run. The
# target-decoy estimate has been measured failing three ways: admitting off-RT
# targets statistically identical to decoys, collapsing to zero identifications
# while reporting healthy diagnostics, and swinging ~10% between refits. Decoys
# are constructed, so the classifier can learn what construction looks like;
# entrapment peptides are real peptides absent from the sample, so every one
# reported is a genuine false positive with no construction signature.
add_test(NAME entrapment_fdp_matches_claimed_fdr COMMAND odia_entrapment)

# The adopted RT calibration (src/calibrate, from OpenDIAlyzer odia-v0.3.0).
# Checks it recovers a non-linear warp from anchors that include gross
# outliers, and that the result is monotone -- a map that reorders the
# gradient is wrong however well it fits.
add_test(NAME rt_calibration_recovers_warp COMMAND odia_rt_calibration)

# The CiRT seed's LINE fit, in the regime doc/19 §3 budgets for: ~50% of apexes
# placed on interference by a blind search over an uncalibrated run. Also pins
# the refusal half of the contract -- noise must produce a residual large enough
# to trip both the decoy control and the 10%-of-run abort. What shipped before
# accepted a statistic that passed 81% of decoys, so the fit was never the thing
# being checked.
add_test(NAME robust_line_survives_half_wrong_anchors COMMAND odia_robust_line)
add_test(NAME anchor_stop_policies_pinned COMMAND odia_anchor_stop)

# A subset must outlive the library it came from. `StringArena::Entry` is a raw
# pointer into the arena's blocks, so the old wholesale arena copy left the
# subset resolving into the PARENT's storage -- correct until the parent dies.
# `-min_library_fragments` is written `library = library.subsetByIndex(...)`,
# which kills the parent, and the run segfaulted 52 s in. The test destroys the
# parent before reading a single string.
add_test(NAME library_subset_outlives_parent COMMAND odia_library_subset)

# PEP must be low where the score is high. The monotonicity sweep ran from the
# low-score end with a running max, giving every entry the worst PEP below it --
# monotone, in the wrong direction, so a monotonicity-only check would pass it.
# This pins the direction AND that the column keeps its range.
add_test(NAME pep_rises_as_score_falls COMMAND odia_pep_monotone)

# v1.18 learner fixes and the offline scorer replay. The selftest scores a synthetic candidate
# table with the engine's own parameter construction, writes it in -out_scorer_input format, runs
# every new flag (each must train, identify and change the output) and pins the GBT defects the
# flags fix. The replay must then reproduce that table's DScore BIT-FOR-BIT -- at the default
# thread count and at 3 threads -- which is the replay-reproduces-the-engine gate in miniature.
add_test(NAME replay_selftest
         COMMAND odia_replay_selftest "${CMAKE_CURRENT_BINARY_DIR}/replay_table.tsv")
set_tests_properties(replay_selftest PROPERTIES FIXTURES_SETUP replay_table)
add_test(NAME replay_reproduces_engine
         COMMAND odia_scorer_replay -in "${CMAKE_CURRENT_BINARY_DIR}/replay_table.tsv" -expect_identical)
add_test(NAME replay_reproduces_engine_3_threads
         COMMAND odia_scorer_replay -in "${CMAKE_CURRENT_BINARY_DIR}/replay_table.tsv" -threads 3 -expect_identical)
set_tests_properties(replay_reproduces_engine replay_reproduces_engine_3_threads PROPERTIES
                     FIXTURES_REQUIRED replay_table PASS_REGULAR_EXPRESSION "expect_identical: PASS")

# The retention-time axis of a precursor that falls in more than one isolation
# window. Data-free: two overlapping windows and four cycles each, built in the
# test. Nothing in the suite constructed a Chromatograms before this.
add_test(NAME extract_rt_axis COMMAND odia_extract_axis)

# Extraction against a scripted run: every spectrum, peak, mobility and window
# written out in the test, so a case is exactly the hazard it names. Each is a
# defect that reached main because nothing here called extract() at all.
add_test(NAME extract_invalid_product_mz COMMAND odia_extract_cases invalid_mz)
add_test(NAME extract_aggregate COMMAND odia_extract_cases aggregate)
add_test(NAME extract_precursor_mobility COMMAND odia_extract_cases mobility)
add_test(NAME extract_im_gating_is_independent COMMAND odia_extract_cases im_gating)
add_test(NAME extract_mobility_band_edge COMMAND odia_extract_cases band_edge)

# Ms1Traces::build matches frames on several threads (doc/83 F05). The matrix,
# the retained residual list (hashed) and its median must be bit-identical to
# the serial build, and the decoder asked for the serial ranges, at every
# thread count, on a scripted run that crosses the residual cap inside a unit
# and has blocks that are not multiples of the 16-frame unit. Identity on
# the fixture TSV alone would not show a wrong cell outside every precursor's
# scoring window; this compares every cell.
add_test(NAME ms1_parallel_identity COMMAND odia_ms1_parallel)

# The sliding window, which is what stops memory scaling with the library: a
# precursor's chromatogram is allocated when the pass reaches the start of its
# retention-time window and freed when the pass leaves it. Twenty precursors
# tiled across the gradient with narrow windows, so "only the overlapping ones
# are resident" is a claim the geometry can actually falsify -- at the 600 s
# window a real run uses, 65% of a library is live at once and a test built at
# that ratio would pass with the whole mechanism removed.
#
# Also covers the two things that must NOT change with it: a streaming consumer
# is handed exactly what the flat array holds, and capping the live set into
# chunks changes the number of passes and not the answer.
add_test(NAME extract_sliding_window COMMAND odia_extract_cases sliding)

# Chunking must be invisible to whatever consumes the stream. Retention times
# that are NOT exact in float, because chunk spectrum ranges were searched with
# float keys in the double times and `sliding`'s whole seconds hid it: a chunk's
# first or last frame went unread and its cells stayed zero. Also asserts the
# hand-over ORDER is the unchunked one (Gate C calibrates on arrival order) and
# that the cap bounds the live set exactly. Fails on a48c223 on all three.
add_test(NAME extract_chunk_invariant COMMAND odia_extract_cases chunk_invariant)

# An explicit -max_live_precursors must survive a byte budget that binds
# nothing: two overlapping precursors with zero valid transitions cost zero
# bytes, the budget inverted to 0, and min(cap, 0) erased the explicit cap of 1
# (peak live 2). Fails on b7c602f.
add_test(NAME extract_cap_zero_rows COMMAND odia_extract_cases cap_zero_rows)

# The chunking guarantee against an INDEPENDENT brute force, not only against
# the unchunked path: decode_block 1/64/127/128/129/200/256 x caps 0/2/4/7 x
# threads 1/4, every trace compared bit for bit in all five planes (intensity
# and both residual pairs) and in hand-over order; restricted rt_low/rt_high
# between a spectrum's double time and its float rounding, on both sides; chunks
# that begin or end between the two windows of one diaPASEF frame; empty
# spectra, missing acquisitions, no-window / outside-run / NaN-iRT /
# no-transition / no-valid-transition precursors; explicit and byte caps
# together. Fails on a48c223 (cells, order, cap, restricted bounds) and on
# b7c602f (decode blocks off the multiples of 128).
add_test(NAME extract_chunk_matrix COMMAND odia_extract_cases chunk_matrix)
add_test(NAME extract_no_assignments COMMAND odia_extract_cases no_assignments)

# Gate C hash calibration (ROUND1_PLAN D6, -gate_calibration hash): the plan's
# pointCount(0) equals the hand-over's; the sample is the N smallest identity
# hashes (brute force) and the same identities in a reversed library; tau and
# every gate decision identical across chunk caps 0/2/3/7/11, decode blocks
# 1/129/200/256 and 1-2 threads, all made against the frozen tau (no warm-up);
# arrival mode on the same data does have a warm-up (control); too few eligible
# decoys leaves the gate unarmed and admitting; the identity hash is pinned.
add_test(NAME gate_calibration_hash_is_order_free COMMAND odia_gate_calibration
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# More than 2^32 chromatogram points, which the 32-bit CSR offsets refused
# outright at 266,664 precursors -- against the 4,255,113 the library generator
# emits. The point array is preallocated as transitions x cycles, so crossing
# that line costs 17.2 GiB and there is no cheaper way to cross it: the defect
# is an OFFSET that overflows, and an offset past 2^32 only exists once 2^32
# points do.
#
# Guarded on the machine having the memory, and SAID rather than skipped
# silently, because a test that quietly vanishes is how the brute-force check
# came to model an extractor four commits out of date.
cmake_host_system_information(RESULT ODIA_TOTAL_RAM_MB QUERY TOTAL_PHYSICAL_MEMORY)
if(ODIA_TOTAL_RAM_MB GREATER 32768)
  add_test(NAME extract_past_32bit_csr COMMAND odia_extract_cases wide_csr)
  set_tests_properties(extract_past_32bit_csr PROPERTIES TIMEOUT 900
                       RUN_SERIAL TRUE)
else()
  message(STATUS "extract_past_32bit_csr skipped; it needs ~17.2 GiB and this "
                 "machine reports ${ODIA_TOTAL_RAM_MB} MiB")
endif()

# The extraction against a brute-force scan of the same spectra, on real runs.
#
# odia_extract has carried this check since the extractor existed and it had NO
# add_test: ODIA_EXTRACT_VERIFY never ran anywhere, which is how the brute force
# came to model an extractor four commits out of date (it still took a maximum
# after the default became Sum, and applied no mobility test and no mass
# calibration). Registered here in three arms, because a single arm only ever
# covers the defaults:
#
#   .._brute_force            12_80, the shipped defaults
#   .._brute_force_calibrated 12_80 with a fitted mass offset, an m/z-dependent
#                             term and Max, over a deliberately wide 100 ppm
#                             box -- the paths a default run never reaches, and
#                             wide enough that Sum and Max differ at all (at 10
#                             ppm on a centroided run they almost never do, so
#                             a narrow arm cannot tell them apart: forcing the
#                             extractor to Sum passes the default arm and fails
#                             this one on 270 points)
#   .._brute_force_mobility   IH1, whose frames co-pack two windows and whose
#                             library carries 1/K0, so both mobility tests are
#                             exercised -- also at 100 ppm, so the mobility
#                             tests rather than the mass box are what decides.
#                             Dropping the per-precursor window disagrees on
#                             108 points here.
#
# Each arm asserts that some of the points it compared were NON-ZERO: most of a
# slice is zeros, an all-zero comparison agrees with anything, and a check that
# quietly became vacuous is what this whole exercise is about.
#
# Each is a slice of the run rather than the whole of it: the mzPeak reader
# costs ~6-7 ms a spectrum, so the slice is what keeps these in the seconds.
if(EXISTS "${ODIA_DATA_DIR}/12_80.mzpeak" AND
   EXISTS "${ODIA_FIXTURE_DIR}/diann_library_small.tsv")
  add_test(NAME extract_matches_brute_force
           COMMAND ${CMAKE_COMMAND} -E env ODIA_EXTRACT_VERIFY=1
                   $<TARGET_FILE:odia_extract>
                   "${ODIA_FIXTURE_DIR}/diann_library_small.tsv"
                   "${ODIA_DATA_DIR}/12_80.mzpeak" 500 300 340)
  set_tests_properties(extract_matches_brute_force PROPERTIES
    PASS_REGULAR_EXPRESSION "decoded spectra, [1-9][0-9]* of [0-9]+ points non-zero, 0 disagree")

  add_test(NAME extract_matches_brute_force_calibrated
           COMMAND ${CMAKE_COMMAND} -E env ODIA_EXTRACT_VERIFY=1
                   ODIA_EXTRACT_PPM_OFFSET=-9.8 ODIA_EXTRACT_PPM_LOG_SLOPE=3.5
                   ODIA_EXTRACT_AGGREGATE=max ODIA_EXTRACT_PPM=100
                   $<TARGET_FILE:odia_extract>
                   "${ODIA_FIXTURE_DIR}/diann_library_small.tsv"
                   "${ODIA_DATA_DIR}/12_80.mzpeak" 500 300 340)
  set_tests_properties(extract_matches_brute_force_calibrated PROPERTIES
    PASS_REGULAR_EXPRESSION "decoded spectra, [1-9][0-9]* of [0-9]+ points non-zero, 0 disagree")
endif()

if(EXISTS "${ODIA_DATA_DIR}/IH1_diaPASEF.mzpeak" AND
   EXISTS "${ODIA_FIXTURE_DIR}/diann_library_small.tsv")
  add_test(NAME extract_matches_brute_force_mobility
           COMMAND ${CMAKE_COMMAND} -E env ODIA_EXTRACT_VERIFY=1
                   ODIA_EXTRACT_VERIFY_MAX=200 ODIA_EXTRACT_PPM=100
                   $<TARGET_FILE:odia_extract>
                   "${ODIA_FIXTURE_DIR}/diann_library_small.tsv"
                   "${ODIA_DATA_DIR}/IH1_diaPASEF.mzpeak" 0 300 320)
  set_tests_properties(extract_matches_brute_force_mobility PROPERTIES
    PASS_REGULAR_EXPRESSION "decoded spectra, [1-9][0-9]* of [0-9]+ points non-zero, 0 disagree")
endif()

# The diaPASEF mobility band derived from the mis-attached selected ions. It
# used to live inside MzPeakSource's constructor, where nothing could reach it
# without a file; these cases are the refusals and the order independence.
add_test(NAME mobility_band_derivation COMMAND odia_mobility_bands)

# The per-run fragment mass calibration (src/calibrate/MassCalibration.cpp).
# Runs on synthesised residuals only, so it needs no data: it checks that a
# known offset is recovered from a narrow peak on a uniform background, that a
# purely uniform sample is REJECTED rather than fitted, that an m/z-dependent
# offset comes back as a trend rather than as its average, and that a residual
# which is flat in m/z yields the constant model instead of a spurious slope.
add_test(NAME mass_calibration_recovers_offset COMMAND odia_mass_calibration)

# The per-run ion-mobility calibration (src/calibrate/MobilityCalibration.cpp).
# Runs on a SCRIPTED run written out inside the test -- every spectrum, peak and
# mobility -- so it needs no data and exercises collect() as well as fit(): the
# half that can be silently wrong is the half that decides what an observation
# of 1/K0 is. It checks that a per-charge offset is recovered and not pooled,
# that an m/z-dependent one comes back as a shape, that fragments agreeing at a
# MEANINGLESS mobility are refused rather than fitted, that a control as peaked
# as the data is refused, that a run with no 1/K0 at all (the SCIEX SWATH case)
# is a named no-op rather than a failed fit, and that each fold's anchors are
# corrected by a fit that never saw them.
add_test(NAME mobility_calibration_recovers_offset COMMAND odia_mobility_calibration)

# The RT-blocked mass recalibration: does it recover a drift that is there, and
# refuse one that is not? It caught the composite-clamp defect on its first run
# -- the same defect MobilityCalibration had, repeated hours later in a sibling
# class, which is exactly what a test for a new model is for.
add_test(NAME mass_recalibration_recovers_drift COMMAND odia_mass_recalibration)

# Sizing the fragment window from the run's own identifications. Guards the two
# ways this estimator can lie: measuring the precision of a group median rather
# than the scatter of a fragment, and reading a distribution that was truncated
# by the window it was observed through.
add_test(NAME mass_width_from_identifications COMMAND odia_mass_width)

# Scoring a mass model on fragments it was NOT fitted from. `fit` reports its
# systematic residual against its own sample, which measures description rather
# than prediction; doc/15 step 1 asks whether identification-derived anchors
# beat the probe, and only a held-out number can answer that. The wrong-shape
# case is the one that matters: a model fitted to the opposite trend must come
# out worse than its own constant.
add_test(NAME mass_anchors_heldout COMMAND odia_mass_anchors)
add_test(NAME prefilter_selection COMMAND odia_prefilter)
add_test(NAME empty_trace_gate COMMAND odia_empty_trace)

# The per-run retention-time refiner. Every case is a failure it has actually
# had on real data: an unbounded correction that took IH1 to zero, contaminated
# anchors a squared loss chased, and a model file whose weight count did not
# match being read as valid.
add_test(NAME rt_refiner COMMAND odia_rt_refiner)

# The two TSV writers, byte for byte against the std::ofstream form they
# replace.
#
# Both were formatting-bound and not disk-bound -- `ofstream <<` measured at
# 25.4 MB/s against 1,042 MB/s for fwrite of pre-rendered bytes on the same
# NVMe -- so both now render through std::to_chars. That buys ~40x against
# exactly one risk: to_chars defaults to shortest round-trip where the stream
# defaults to 6 significant digits, so a wrong precision does not round
# differently, it writes a different file. Nothing else in the suite reads
# either file's numbers closely enough to notice.
#
# `files` is the end-to-end diff; `formatting` is the sweep that covers the
# values a fixture cannot hold -- subnormals, extreme exponents, NaN and inf.
add_test(NAME tsv_writers_are_byte_identical
         COMMAND odia_tsv_writers files "${CMAKE_CURRENT_BINARY_DIR}")
# -out_fragvec is a JOIN into -out, so its header order, its key and its
# precision are load-bearing in a way no other writer's are.
add_test(NAME fragvec_tsv_joins_to_out
         COMMAND odia_tsv_writers fragvec "${CMAKE_CURRENT_BINARY_DIR}")
# -fragvec_scores: the appended schema is the 78-column contract, the constant
# guard is its only filter (indicators included), the cell the fit and -out read
# is the exported float32, and the permuted control moves whole rows by a fixed,
# label-blind permutation.
add_test(NAME fragvec_scores_block
         COMMAND odia_tsv_writers fragvec_scores "${CMAKE_CURRENT_BINARY_DIR}")
add_test(NAME tsv_number_formatting_matches_ostream
         COMMAND odia_tsv_writers formatting)
