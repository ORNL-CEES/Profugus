##---------------------------------------------------------------------------##
## harness/CMakeLists.txt
## Thomas M. Evans
## Tuesday June 21 11:21:46 2011
##---------------------------------------------------------------------------##
## Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for harness
##---------------------------------------------------------------------------##

INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsSubPackageMacros)
INCLUDE(CheckTypeSize)

TRIBITS_SUBPACKAGE(Harness)

##---------------------------------------------------------------------------##
## PACKAGE-SPECIFIC OPTIONS
##---------------------------------------------------------------------------##
## DBC

# set DBC to full on if enable-debug
IF (${PROJECT_NAME}_ENABLE_DEBUG)
  SET(UTILS_DBC 7 CACHE STRING "Value of Design-by-Contract [0-7]")
ELSE()
  SET(UTILS_DBC 0 CACHE STRING "Value of Design-by-Contract [0-7]")
ENDIF()

MESSAGE(STATUS "Value of DBC in Nemesis set to ${UTILS_DBC}")

##---------------------------------------------------------------------------##
## DIAGNOSTICS

SET(UTILS_DIAGNOSTICS 1
  CACHE STRING
  "Value of nemesis diagnostics, 1 by default")

##---------------------------------------------------------------------------##
## WARNINGS

TRIBITS_ADD_OPTION_AND_DEFINE(ENABLE_UTILS_WARNINGS
  UTILS_WARNINGS_ENABLED
  "Enable storage of warnings for later retreival"
  ON)

TRIBITS_ADD_OPTION_AND_DEFINE(ENABLE_UTILS_IMMEDIATE_WARN
  UTILS_WARNINGS_IMMEDIATE
  "Turn on immediate output of all warnings (default: no)"
  OFF)

IF(UTILS_IMMEDIATE_WARN AND NOT UTILS_WARNINGS)
  MESSAGE(STATUS "Disabling 'immediate warn' option because warnings are off")
  SET(UTILS_WARNINGS_IMMEDIATE OFF)
ENDIF()

##---------------------------------------------------------------------------##
## TYPE SIZES

CHECK_TYPE_SIZE(short            SIZEOF_SHORT)
CHECK_TYPE_SIZE(unsigned         SIZEOF_UNSIGNED)
CHECK_TYPE_SIZE("unsigned short" SIZEOF_UNSIGNED_SHORT)
CHECK_TYPE_SIZE(int              SIZEOF_INT)
CHECK_TYPE_SIZE(long             SIZEOF_LONG)
CHECK_TYPE_SIZE("long long"      SIZEOF_LONG_LONG)
CHECK_TYPE_SIZE(float            SIZEOF_FLOAT)
CHECK_TYPE_SIZE(double           SIZEOF_DOUBLE)

##---------------------------------------------------------------------------##
## PACKAGE CONFIGURE FILE
##---------------------------------------------------------------------------##

TRIBITS_CONFIGURE_FILE(config.h)

##---------------------------------------------------------------------------##
## HEADER AND SOURCE FILES
##---------------------------------------------------------------------------##

# C/C++ headers
FILE(GLOB HEADERS *.hh)
LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/config.h)

# C/C++ source
SET(SOURCES
  DBC.cc
  Diagnostics.cc
  Warnings.cc
  )

##---------------------------------------------------------------------------##
## PACKAGE TARGETS AND INSTALL
##---------------------------------------------------------------------------##

TRIBITS_ADD_LIBRARY(
  Utils_harness
  SOURCES ${SOURCES})

INSTALL(FILES
  ${HEADERS}
  DESTINATION
  include/Utils/harness)

TRIBITS_ADD_TEST_DIRECTORIES(test)

#
# Add documentation for this package
#
DENOVO_ADD_DOC()

TRIBITS_SUBPACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
##                   end of harness/CMakeLists.txt
##---------------------------------------------------------------------------##
