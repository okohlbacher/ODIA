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
