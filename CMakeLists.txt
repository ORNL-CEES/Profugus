##---------------------------------------------------------------------------##
## Profugus/CMakeLists.txt
##---------------------------------------------------------------------------##
## Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
# A) Define your project name and set up major project options
##---------------------------------------------------------------------------##

INCLUDE(${CMAKE_SOURCE_DIR}/ProjectName.cmake)

# CMake requires that you declare the CMake project in the top-level file and
# not in an include file :-(
PROJECT(${PROJECT_NAME} NONE)

# Add support for extra repositories; everything handled by native above
SET(${PROJECT_NAME}_SUPPORT_EXTRA_REPOS FALSE)

# Disable Tribits generating HTML dependencies webpage and xml files
SET(${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE OFF CACHE BOOL "" )

# Disable Tribits export system to save time configuring
SET(${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES OFF CACHE BOOL "")

# Disable Tribits export makefiles system to save time configuring
SET(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES OFF CACHE BOOL "")

##---------------------------------------------------------------------------##
# B) Pull in the TriBITS system and execute
##---------------------------------------------------------------------------##

# Just use the TriBITS system in Trilinos
SET(${PROJECT_NAME}_TRIBITS_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/TriBITS" CACHE PATH "" )
MESSAGE("${PROJECT_NAME}_TRIBITS_DIR = '${${PROJECT_NAME}_TRIBITS_DIR}'")
INCLUDE(${${PROJECT_NAME}_TRIBITS_DIR}/TriBITS.cmake)

# CMake requires this be in the top file and not in an include file :-(
CMAKE_MINIMUM_REQUIRED(VERSION ${TRIBITS_CMAKE_MINIMUM_REQUIRED})

# Default to building shared libraries if user hasn't specified.
#
# Note that if you define BUILD_SHARED_LIBS using a configure option file
# with -D PROFUGUS_PROJECT_CONFIGURE_FILE:PATH=config.cmake,
# then you'll have to use the FORCE keyword in that file because it gets read
# during TRIBITS_PROJECT().
#
#   SET(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
#
# (Alternatively, use 'cmake -C config.cmake' to read in BUILD_SHARED_LIBS
# *before* this point.)
IF(NOT DEFINED BUILD_SHARED_LIBS)
  SET(BUILD_SHARED_LIBS_DEFAULT ON)
  IF(DEFINED BUILD_SHARED)
    SET(BUILD_SHARED_LIBS_DEFAULT ${BUILD_SHARED})
  ENDIF()
  SET(BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS_DEFAULT}
    CACHE BOOL
    "Build shared libraries.")
  UNSET(BUILD_SHARED_LIBS_DEFAULT)
ENDIF()

# Do all of the processing for this Tribits project
TRIBITS_PROJECT()

# Add the final documentation
#add_subdirectory(doc)

##---------------------------------------------------------------------------##
## end of CMakeLists.txt
##---------------------------------------------------------------------------##
