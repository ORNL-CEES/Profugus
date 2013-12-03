##---------------------------------------------------------------------------##
## gtest/CMakeLists.txt
## Thomas M. Evans
## Tuesday June 21 11:21:46 2011
##---------------------------------------------------------------------------##
## Copyright (C) 2012 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for gtest
##---------------------------------------------------------------------------##

INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsSubPackageMacros)

TRIBITS_SUBPACKAGE(Gtest)

##---------------------------------------------------------------------------##
## PACKAGE-SPECIFIC OPTIONS
##---------------------------------------------------------------------------##

# for now always turn off pthreads
SET (GTEST_HAS_PTHREAD OFF CACHE BOOL "Set GTEST pthreads")

# Build shared libraries?
SET(GTEST_CREATE_SHARED_LIBRARY OFF)
IF (BUILD_SHARED_LIBS)
  SET(GTEST_CREATE_SHARED_LIBRARY ON)
ENDIF()

# Turn off TUPLES when using INTEL
SET(GTEST_HAS_TR1_TUPLE ON)
IF (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  SET(GTEST_HAS_TR1_TUPLE OFF)
ENDIF()

# Enable c++11 support if appropriate
SET(GTEST_LANG_CXX11 OFF)
IF (${PROJECT_NAME}_ENABLE_CXX11)
  SET(GTEST_LANG_CXX11 ON)
ENDIF()

# Turn off variadic-macros warnings when using gtest
IF (CMAKE_COMPILER_IS_GNUCXX AND NOT WIN32)
  SET(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -Wno-variadic-macros")
ENDIF()

##---------------------------------------------------------------------------##
## PACKAGE CONFIGURE FILE
##---------------------------------------------------------------------------##

TRIBITS_CONFIGURE_FILE(config.h)

##---------------------------------------------------------------------------##
## HEADER AND SOURCE FILES
##---------------------------------------------------------------------------##

FILE(GLOB HEADERS *.hh)
LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/config.h)

SET(SOURCES
  Gtest_Functions.cc
  gtest-all.cc
  )

##---------------------------------------------------------------------------##
## PACKAGE TARGETS AND INSTALL
##---------------------------------------------------------------------------##

# Add main library that gets explicitly linked downstream

TRIBITS_ADD_LIBRARY(
  Nemesis_gtest
  SOURCES ${SOURCES}
  TESTONLY
  )

##---------------------------------------------------------------------------##
#
# Add test directory
#
TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_SUBPACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
##                   end of gtest/CMakeLists.txt
##---------------------------------------------------------------------------##
