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

##---------------------------------------------------------------------------##
## A) Define the package
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE_DECL(Utils)

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

# to allow includes like #include "comm/config.h"
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

##---------------------------------------------------------------------------##
# C) Process the subpackages for Thyra
##---------------------------------------------------------------------------##

TRIBITS_PROCESS_SUBPACKAGES()

# Exclude files for source package.

TRIBITS_PACKAGE_DEF()

# Exclude files for source package.
TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

##---------------------------------------------------------------------------##
## D) Do standard postprocessing
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
##                      end of nemesis/CMakeLists.txt
##---------------------------------------------------------------------------##
