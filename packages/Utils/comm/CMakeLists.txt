##---------------------------------------------------------------------------##
## comm/CMakeLists.txt
## Thomas M. Evans
## Wednesday June 22 11:14:40 2011
##---------------------------------------------------------------------------##
## Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for comm
##---------------------------------------------------------------------------##

INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsSubPackageMacros)
INCLUDE(CheckTypeSize)
INCLUDE(CheckIncludeFile)
INCLUDE(CheckLibraryExists)

TRIBITS_SUBPACKAGE(Comm)

##---------------------------------------------------------------------------##
## PACKAGE-SPECIFIC OPTIONS
##---------------------------------------------------------------------------##
## MPI

IF (TPL_ENABLE_MPI)
  SET(COMM_MPI ON)
ELSE()
  SET(COMM_SCALAR ON)
ENDIF()

## Set Comm MPI implementation
IF (DEFINED Nemesis_MPI_IMP)
  SET(MPI_IMP ${Nemesis_MPI_IMP})
ENDIF()

##---------------------------------------------------------------------------##
## TIMING

IF (NOT NEMESIS_TIMING)
  SET(NEMESIS_TIMING 1
          CACHE STRING
          "Value for Nemesis timing, 1 by default")
ENDIF()

##---------------------------------------------------------------------------##
## PARALLEL OUTPUT

IF (NOT NEMESIS_POUT)
  SET(NEMESIS_POUT ON CACHE BOOL "Turn on Nemesis Parallel Output.")
ENDIF()

##---------------------------------------------------------------------------##
## INCLUDES

CHECK_INCLUDE_FILE(sys/times.h HAVE_SYS_TIMES_H)
CHECK_INCLUDE_FILE(unistd.h    HAVE_UNISTD_H)
CHECK_INCLUDE_FILE(windows.h   HAVE_WINDOWS_H)

##---------------------------------------------------------------------------##
## POSIX FLAGS

IF (DEFINED _POSIX_C_SOURCE)
  OPTION(_POSIX_SOURCE "Do not include non-POSIX code." ON)
ENDIF()

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
  MPI.cc
  MPI_Blocking.pt.cc
  MPI_Gathers.pt.cc
  MPI_Non_Blocking.pt.cc
  MPI_Reductions.pt.cc
  P_Stream.cc
  Parallel_Unit_Test.cc
  Parallel_Utils.cc
  Request.cc
  Serial.cc
  SpinLock.cc
  Timer.cc
  Timing_Diagnostics.cc
  )

##---------------------------------------------------------------------------##
## PACKAGE TARGETS AND INSTALL
##---------------------------------------------------------------------------##

TRIBITS_ADD_LIBRARY(
  Nemesis_comm
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES})

INSTALL(FILES
  ${HEADERS}
  DESTINATION
  include/Nemesis/comm)

#
# Add test directory for this package
#
TRIBITS_ADD_TEST_DIRECTORIES(test)

#
# Add documentation for this package
#
DENOVO_ADD_DOC()

TRIBITS_SUBPACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
##                   end of comm/CMakeLists.txt
##---------------------------------------------------------------------------##
