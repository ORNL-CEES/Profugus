##---------------------------------------------------------------------------##
## Utils/CMakeLists.txt
## Thomas M. Evans
## Monday June 20 15:40:52 2011
##---------------------------------------------------------------------------##
## Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## TriBITS package file for Utils package
##---------------------------------------------------------------------------##

INCLUDE(CheckTypeSize)
INCLUDE(ProfugusAddCFlags)
INCLUDE(TestBigEndian)

# Add Utils cmake services
APPEND_SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

##---------------------------------------------------------------------------##
## A) Define the package
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE(Utils)

##---------------------------------------------------------------------------##
## B) Set up package-specific options
##---------------------------------------------------------------------------##

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

IF(CMAKE_COMPILER_IS_GNUCXX AND NOT WIN32)
  PROFUGUS_ADD_CXX_FLAGS(
    -Wno-sign-compare
    -Wno-unused-but-set-variable
    -Wno-unused-variable)
ENDIF()

IF (DEFINED Utils_MPI_IMP)
  MESSAGE(STATUS "Found MPI implementation ${Utils_MPI_IMP}")
ENDIF()

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
## DIAGNOSTICS

SET(UTILS_DIAGNOSTICS 1
  CACHE STRING
  "Value of utils diagnostics, 1 by default")

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
## MPI

IF (TPL_ENABLE_MPI)
  SET(COMM_MPI ON)
ELSE()
  SET(COMM_SCALAR ON)
ENDIF()

## Set Comm MPI implementation
IF (DEFINED Utils_MPI_IMP)
  SET(MPI_IMP ${Utils_MPI_IMP})
ENDIF()

##---------------------------------------------------------------------------##
## TIMING

IF (NOT UTILS_TIMING)
  SET(UTILS_TIMING 1
          CACHE STRING
          "Value for Utils timing, 1 by default")
ENDIF()

##---------------------------------------------------------------------------##
## PARALLEL OUTPUT

IF (NOT UTILS_POUT)
  SET(UTILS_POUT ON CACHE BOOL "Turn on Utils Parallel Output.")
ENDIF()

##---------------------------------------------------------------------------##
## INCLUDES

CHECK_INCLUDE_FILE(sys/times.h HAVE_SYS_TIMES_H)
CHECK_INCLUDE_FILE(unistd.h    HAVE_UNISTD_H)
CHECK_INCLUDE_FILE(windows.h   HAVE_WINDOWS_H)

##---------------------------------------------------------------------------##
## GTEST options

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
## POSIX FLAGS

IF (DEFINED _POSIX_C_SOURCE)
  OPTION(_POSIX_SOURCE "Do not include non-POSIX code." ON)
ENDIF()

##---------------------------------------------------------------------------##
## RNG options

# NOTE: Currently, the lagged-fibonaci generator does not require any
# LittleEndian flag declarators; however, the congruential generator (lcg)
# does thus we may need to check for LittleEndian in any SPRNG implementation
# we bring over later.
SET(SPRNG_BigEndian OFF)
SET(SPRNG_LittleEndian ON)
TEST_BIG_ENDIAN(SPRNG_BigEndian)
IF(SPRNG_BigEndian)
  SET(SPRNG_LittleEndian OFF)
ENDIF()

##---------------------------------------------------------------------------##
## Optional TPLs

IF (TPL_ENABLE_HDF5)
  MESSAGE(STATUS "Building Utils with HDF5 support")
  SET(USE_HDF5 ON)
ENDIF()

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
INSTALL(FILES ${HARNESS_HEADERS} DESTINATION include/Utils/harness)

# comm sources

FILE(GLOB COMM_HEADERS comm/*.hh)
SET(COMM_SOURCES
  comm/Logger.cc
  comm/Logger_Statement.cc
  comm/MPI.cc
  comm/MPI_Blocking.pt.cc
  comm/MPI_Gathers.pt.cc
  comm/MPI_Non_Blocking.pt.cc
  comm/MPI_Reductions.pt.cc
  comm/P_Stream.cc
  comm/Parallel_Utils.cc
  comm/Request.cc
  comm/Serial.cc
  comm/SpinLock.cc
  comm/Timer.cc
  comm/Timing_Diagnostics.cc)
LIST(APPEND HEADERS ${COMM_HEADERS})
LIST(APPEND SOURCES ${COMM_SOURCES})
INSTALL(FILES ${COMM_HEADERS} DESTINATION include/Utils/comm)

# utils sources

FILE(GLOB UTILS_HEADERS utils/*.hh)
SET(UTILS_SOURCES
  utils/Hash_Functions.cc
  utils/String_Functions.cc)

IF(TPL_ENABLE_HDF5)
  APPEND_SET(UTILS_SOURCES
    utils/HDF5_IO.cc
    utils/HDF5_Reader.cc
    utils/Hash_Functions.cc
    utils/Member_Manager.cc
    utils/Parallel_HDF5_Writer.cc
    utils/Serial_HDF5_Writer.cc
    utils/String_Functions.cc
    utils/Vector_Functions.cc
    )
ENDIF()

LIST(APPEND HEADERS ${UTILS_HEADERS})
LIST(APPEND SOURCES ${UTILS_SOURCES})
INSTALL(FILES ${UTILS_HEADERS} DESTINATION include/Utils/utils)

# rng sources

FILE(GLOB RNG_HEADERS rng/*.hh)
SET(RNG_SOURCES
  rng/RNG.cc
  rng/RNG_Control.cc
  rng/sprng/checkid.c
  rng/sprng/communicate.c
  rng/sprng/lfg.c
  rng/sprng/makeseed.c
  rng/sprng/memory.c
  rng/sprng/simple_mpi.c
)

LIST(APPEND HEADERS ${RNG_HEADERS})
LIST(APPEND SOURCES ${RNG_SOURCES})
INSTALL(FILES ${RNG_HEADERS} DESTINATION include/Utils/rng)

# gtest sources

FILE(GLOB GTEST_HEADERS gtest/*.hh)
SET(GTEST_SOURCES
  gtest/Gtest_Functions.cc
  gtest/gtest-all.cc)

LIST(APPEND SOURCES ${GTEST_SOURCES})

LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/config.h)
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/config.h DESTINATION include/Utils)

##---------------------------------------------------------------------------##
## LIBRARY AND INSTALL TARGETS

TRIBITS_ADD_LIBRARY(
  Utils
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES})

##---------------------------------------------------------------------------##
# Add tests to this package

TRIBITS_ADD_TEST_DIRECTORIES(
  harness/test
  utils/test
  rng/test
  cxx11
  )

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
