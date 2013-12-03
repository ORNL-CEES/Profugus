##---------------------------------------------------------------------------##
## Utils/CMakeLists.txt
## Thomas M. Evans
## Monday June 20 15:40:52 2011
##---------------------------------------------------------------------------##
## Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## TriBITS package file for Utils package
##---------------------------------------------------------------------------##

INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsLibraryMacros)
INCLUDE(CheckTypeSize)

##---------------------------------------------------------------------------##
## A) Define the package
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE(Utils)

# determine the MPI implementation if it hasn't been given by the user
IF ((NOT DEFINED Utils_MPI_IMP) AND (TPL_ENABLE_MPI))

  # check MPI version--see if its OPENMPI
  FIND_LIBRARY(OPENMPI_DEFINED NAMES open-rte open-palg
    PATHS ${MPI_BASE_DIR}/lib ${TPL_MPI_LIBRARY_DIRS}
    NO_DEFAULT_PATHS)
  IF (NOT (${OPENMPI_DEFINED} MATCHES NOTFOUND))
    GLOBAL_SET(Utils_MPI_IMP MPI_OPENMPI)
  ENDIF()

ENDIF()

#
# Resolve min/max conflict on Windows
#
IF( CMAKE_SYSTEM_NAME MATCHES "Windows" )
   ADD_DEFINITIONS(-DNOMINMAX)
ENDIF( CMAKE_SYSTEM_NAME MATCHES "Windows" )

IF (DEFINED Utils_MPI_IMP)
  MESSAGE(STATUS "Found MPI implementation ${Utils_MPI_IMP}")
ENDIF()

##---------------------------------------------------------------------------##
## B) Set up package-specific options
##---------------------------------------------------------------------------##

# Setup M4 for FORTRAN processing
APPEND_SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Setup debug option
TRIBITS_ADD_DEBUG_OPTION()

# to allow includes like #include "comm/Comm.h"
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# to allow includes like #include "Utils/config.h"
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/..)

##---------------------------------------------------------------------------##
## DBC

# set DBC to full on if enable-debug
IF (${PROJECT_NAME}_ENABLE_DEBUG)
  SET(UTILS_DBC 7 CACHE STRING "Value of Design-by-Contract [0-7]")
ELSE()
  SET(UTILS_DBC 0 CACHE STRING "Value of Design-by-Contract [0-7]")
ENDIF()

MESSAGE(STATUS "Value of DBC in Utils set to ${UTILS_DBC}")

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
# C) Process the subpackages
##---------------------------------------------------------------------------##

TRIBITS_CONFIGURE_FILE(config.h)

##---------------------------------------------------------------------------##
## SOURCES AND HEADERS

SET(SOURCES "")
SET(HEADERS "")

# harness sources

FILE(GLOB HARNESS_HEADERS harness/*.hh)
SET(HARNESS_SOURCES
  harness/DBC.cc
  harness/Diagnostics.cc
  harness/Warnings.cc)
LIST(APPEND HEADERS ${HARNESS_HEADERS})
LIST(APPEND SOURCES ${HARNESS_SOURCES})

##---------------------------------------------------------------------------##
## LIBRARY AND INSTALL TARGETS

TRIBITS_ADD_LIBRARY(
  Utils
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES})

##---------------------------------------------------------------------------##
## FINISH SETUP

# Exclude files for source package.
TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

##---------------------------------------------------------------------------##
## D) Do standard postprocessing
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
##                      end of Utils/CMakeLists.txt
##---------------------------------------------------------------------------##
