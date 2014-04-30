#---------------------------------------------------------------------------##
## Utils/cmake/UtilsTest.cmake
## Thomas M. Evans
## Wednesday July 11 14:35:42 2012
##---------------------------------------------------------------------------##
## Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## Setup Utils test expressions and services
##---------------------------------------------------------------------------##
## Setup PASS/FAIL expressions for Utils test harness

# Utils-test-harness pass/fail
SET(UtilsPass "Test: PASSED")
SET(UtilsFail "Test: FAILED")

# Utils google test harness pass/fail
SET(UtilsGtestPass "overall test result: PASSED")
SET(UtilsGtestFail "overall test result: FAILED")

# Python harness test/fail
SET(ExnihiloPyPass "Test: passed")
SET(ExnihiloPyFail "Test: failed")

# Turn off variadic-macros warnings when using gtest
IF (CMAKE_COMPILER_IS_GNUCXX AND NOT WIN32)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-variadic-macros")
ENDIF()

# Set number of parallel tests to run
IF (TPL_ENABLE_MPI)
  SET(UtilsNP 1 2 4)
ELSE()
  SET(UtilsNP "1")
ENDIF()

INCLUDE(CMakeParseArguments)
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsAddAdvancedTest)

##---------------------------------------------------------------------------##
## ADDING UNIT TEST
##---------------------------------------------------------------------------##
# ADD_UTILS_TEST(
#   SOURCE_FILE
#   [NP 1 [2 [...]]]
#   [DEPLIBS lib1 [lib2 ...]]
#   [ENVIRONMENT VAR=value [VAR2=value2 ...]]
#   [ISOLATE]
#   [LEGACY]
#   [DISABLE]
#   )
#
# Create and add a unit test from the source file SOURCE_FILE.
#
# NP specifies the number of processors to use for this unit test. The default
# is to use UtilsNP (1, 2, and 4) for MPI builds and 1 for serial builds.
#
# DEPLIBS specifies extra libraries to link to. By default, unit tests will link
# against the package's current library; Gtest unit tests will also link against
# Utils_gtest.
#
# ENVRIONMENT sets the given environmental variables when the test is run.
#
# If ISOLATE is specified, the test will be run in its own directory.
#
# If LEGACY is specified, we use the old Utils harness. Otherwise, use the
# Utils google test harness.
#
# If DISABLE is specified, we will build the test executable but omit it from
# the list of tests to run through CTest.
#
FUNCTION(ADD_UTILS_TEST SOURCE_FILE)
  cmake_parse_arguments(PARSE
    "LEGACY;ISOLATE;DISABLE"
    ""
    "DEPLIBS;NP;ENVIRONMENT" ${ARGN})

  # Set googletest/harness options
  IF (PARSE_LEGACY)
    SET(PASS_RE ${UtilsPass})
    SET(FAIL_RE ${UtilsFail})
    SET(DEPLIBS Utils_comm)
  ELSE()
    SET(PASS_RE ${UtilsGtestPass})
    SET(FAIL_RE ${UtilsGtestFail})
    SET(DEPLIBS Utils_gtest)
  ENDIF()

  # Add additional library dependencies
  LIST(APPEND DEPLIBS ${PARSE_DEPLIBS})

  # Set number of processors, defaulting to UtilsNP
  SET(NUM_PROCS ${PARSE_NP})
  IF (NOT NUM_PROCS)
    SET(NUM_PROCS ${UtilsNP})
  ENDIF()

  # Check to see if MPI-only unit test
  LIST(FIND NUM_PROCS 1 HAS_SERIAL)
  IF (HAS_SERIAL EQUAL -1)
    SET(COMM mpi)
    IF(NOT TPL_ENABLE_MPI)
      # return early to avoid potential set_property on nonexistent test
      RETURN()
    ENDIF()
  ELSE()
    SET(COMM serial mpi)
  ENDIF()

  # add the test executable
  GET_FILENAME_COMPONENT(EXE_NAME ${SOURCE_FILE} NAME_WE)
  TRIBITS_ADD_EXECUTABLE(
    ${EXE_NAME}
    SOURCES ${SOURCE_FILE}
    DEPLIBS ${DEPLIBS}
    COMM ${COMM}
    )

  # If the test is disabled, print a small message and omit it from the CTest
  IF (PARSE_DISABLE)
    MESSAGE("Disabling testing for ${SOURCE_FILE} in ${SUBPACKAGE_FULLNAME}")
    RETURN()
  ENDIF()

  # Loop over processors for parallel tests
  FOREACH(np ${NUM_PROCS})
    IF (PARSE_ISOLATE)
      # Add an "advanced" test
      TRIBITS_ADD_ADVANCED_TEST(
        ${EXE_NAME}_MPI_${np}
        TEST_0
          EXEC ${EXE_NAME}
          NUM_MPI_PROCS ${np}
          PASS_REGULAR_EXPRESSION "${PASS_RE}"
          FAIL_REGULAR_EXPRESSION "${FAIL_RE}"
        OVERALL_WORKING_DIRECTORY TEST_NAME
        )
    ELSE()
      # Add a normal test
      TRIBITS_ADD_TEST(
        ${EXE_NAME}
        NUM_MPI_PROCS ${np}
        PASS_REGULAR_EXPRESSION "${PASS_RE}"
        FAIL_REGULAR_EXPRESSION "${FAIL_RE}"
        )
    ENDIF()
  ENDFOREACH()

  # set environmental variables if necessary
  IF(PARSE_ENVIRONMENT)
    FOREACH(np ${NUM_PROCS})
      IF (TPL_ENABLE_MPI)
        SET(TEST_NAME "${SUBPACKAGE_FULLNAME}_${EXE_NAME}_MPI_${np}")
      ELSE()
        SET(TEST_NAME "${SUBPACKAGE_FULLNAME}_${EXE_NAME}")
      ENDIF()

      # Modify environment
      SET_PROPERTY(TEST "${TEST_NAME}"
        PROPERTY ENVIRONMENT
        ${PARSE_ENVIRONMENT}
        )
    ENDFOREACH()
  ENDIF()

ENDFUNCTION()

# SETUP_UTILS_PYTEST(
#   PKGDEPS [subpackage [subpackage2 ...]]
#   [ROOT dir]
#   )
#
# Set up dependencies for Python tests. Include the current subpackage as a
# default dependency.
#
# The ROOT directory option specifies where the test files are to be found. The
# default directory is 'python'.

MACRO(SETUP_UTILS_PYTEST)
  cmake_parse_arguments(PARSE "" "ROOT" "PKGDEPS" ${ARGN})

  # Construct python include path from PACKAGES argument.
  SET(UTILS_PYTHONPATH "${${SUBPACKAGE_FULLNAME}_BINARY_DIR}")
  FOREACH(package ${PARSE_PKGDEPS})
    SET(UTILS_PYTHONPATH "${UTILS_PYTHONPATH}:${${package}_BINARY_DIR}")
  ENDFOREACH()

  SET(UTILS_PYTHONTESTROOT ${PARSE_ROOT})
  IF(NOT UTILS_PYTHONTESTROOT)
    SET(UTILS_PYTHONTESTROOT "python/")
  ENDIF()
ENDMACRO()

# ADD_UTILS_PYTEST(
#   PYTHON_FILE
#   [NP 1 [2 [...]]]
#   [ENVIRONMENT VAR=value [VAR2=value2 ...]]
#   [LEGACY]
#   )
#
# Create and add a unit test from the source file at python/PYTHON_FILE. We also
# modify the test PYTHONPATH to include the subpackages specfied in
# SETUP_UTILS_PYTEST.
#
# NP specifies the number of processors to use for this unit test. The default
# is to use UtilsNP (1, 2, and 4) for MPI builds and 1 for serial builds.
#
# If LEGACY is specified, we use the old "tester.py" module. Otherwise we
# use the pythonic "denovo_unittest.py" harness.

FUNCTION(ADD_UTILS_PYTEST SOURCE_FILE)
  cmake_parse_arguments(PARSE
    "LEGACY"
    ""
    "NP;ENVIRONMENT" ${ARGN})

  # Set harness options
  IF (NOT PARSE_LEGACY)
    SET(PASS_RE ${UtilsGtestPass})
    SET(FAIL_RE ${UtilsGtestFail})
    SET(SYS_ARGV "-v")
  ELSE()
    SET(PASS_RE ${ExnihiloPyPass})
    SET(FAIL_RE ${ExnihiloPyFail})
    SET(SYS_ARGV)
  ENDIF()

  # Set number of processors, defaulting to UtilsNP
  SET(NUM_PROCS ${PARSE_NP})
  IF (NOT NUM_PROCS)
    SET(NUM_PROCS ${UtilsNP})
  ENDIF()

  # Check to see if MPI-only unit test
  LIST(FIND NUM_PROCS 1 HAS_SERIAL)
  IF (HAS_SERIAL EQUAL -1)
    SET(COMM mpi)
    IF(NOT TPL_ENABLE_MPI)
      # return early to avoid set_property on nonexistent test
      RETURN()
    ENDIF()
  ELSE()
    SET(COMM serial mpi)
  ENDIF()

  GET_FILENAME_COMPONENT(PYNAME ${SOURCE_FILE} NAME_WE)

  # Loop over processors for parallel tests
  FOREACH(np ${NUM_PROCS})
    # add the tests
    TRIBITS_ADD_TEST(
      "${PYTHON_EXECUTABLE}"
      NOEXEPREFIX NOEXESUFFIX
      NAME "${PYNAME}"
      ARGS "${CMAKE_CURRENT_SOURCE_DIR}/${UTILS_PYTHONTESTROOT}${PYNAME}.py ${SYS_ARGV}"
      COMM ${COMM}
      NUM_MPI_PROCS ${np}
      PASS_REGULAR_EXPRESSION "${PASS_RE}"
      FAIL_REGULAR_EXPRESSION "${FAIL_RE}"
      )
    IF (TPL_ENABLE_MPI)
      SET(TEST_NAME "${SUBPACKAGE_FULLNAME}_${PYNAME}_MPI_${np}")
    ELSE()
      SET(TEST_NAME "${SUBPACKAGE_FULLNAME}_${PYNAME}")
    ENDIF()

    # Change the python path and any other user environment settings
    SET_PROPERTY(TEST "${TEST_NAME}"
      PROPERTY ENVIRONMENT
      PYTHONPATH=${UTILS_PYTHONPATH}
      ${PARSE_ENVIRONMENT}
      )
  ENDFOREACH()

ENDFUNCTION()

##---------------------------------------------------------------------------##
##                   end of UtilsTest.cmake
##---------------------------------------------------------------------------##
