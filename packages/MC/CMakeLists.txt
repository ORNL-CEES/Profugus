##---------------------------------------------------------------------------##
## MC/CMakeLists.txt
## Thomas M. Evans
## Monday June 20 15:40:52 2011
##---------------------------------------------------------------------------##
## Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## TriBITS package file for MC package
##---------------------------------------------------------------------------##

INCLUDE(AppendSet)
INCLUDE(ProfugusAddCFlags)

# Add Utils cmake to path
APPEND_SET(CMAKE_MODULE_PATH ${Utils_SOURCE_DIR}/cmake)

##---------------------------------------------------------------------------##
## A) Define the package
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE(MC)

##---------------------------------------------------------------------------##
## B) Set up package-specific options
##---------------------------------------------------------------------------##

## ML OPTIONS

## Turn on ML if that package is enabled
IF (${PROJECT_NAME}_ENABLE_ML)
  SET(USE_ML 1)
ENDIF()

# Setup M4 for FORTRAN processing
APPEND_SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

## Trilinos Timing options.
IF (${PROJECT_NAME}_ENABLE_TRILINOS_TIMING)
   MESSAGE(STATUS "Building Profugus with Trilinos timing")
   SET(USE_TRILINOS_TIMING 1)
ENDIF()
 
# Setup debug option
TRIBITS_ADD_DEBUG_OPTION()

IF(CMAKE_COMPILER_IS_GNUCXX AND NOT WIN32)
  PROFUGUS_ADD_CXX_FLAGS(
    -Wno-sign-compare
    -Wno-unused-but-set-variable
    -Wno-unused-variable)
ENDIF()

#
# Resolve min/max conflict on Windows
#
IF( CMAKE_SYSTEM_NAME MATCHES "Windows" )
   ADD_DEFINITIONS(-DNOMINMAX)
ENDIF( CMAKE_SYSTEM_NAME MATCHES "Windows" )

##---------------------------------------------------------------------------##
# C) Process the subpackages
##---------------------------------------------------------------------------##

TRIBITS_CONFIGURE_FILE(config.h)
TRIBITS_PROCESS_SUBPACKAGES()

##---------------------------------------------------------------------------##
## D) Do standard postprocessing
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
##                      end of Utils/CMakeLists.txt
##---------------------------------------------------------------------------##
