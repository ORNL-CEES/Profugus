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

##---------------------------------------------------------------------------##
# B) Pull in the TriBITS system and execute
##---------------------------------------------------------------------------##

# Just use the TriBITS system in Trilinos
SET(${PROJECT_NAME}_TRIBITS_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/Trilinos/cmake/tribits" CACHE PATH "" )
MESSAGE("${PROJECT_NAME}_TRIBITS_DIR = '${${PROJECT_NAME}_TRIBITS_DIR}'")
INCLUDE(${${PROJECT_NAME}_TRIBITS_DIR}/package_arch/TribitsProject.cmake)

# CMake requires this be in the top file and not in an include file :-(
CMAKE_MINIMUM_REQUIRED(VERSION ${TRIBITS_CMAKE_MINIMUM_REQUIRED})

# Add install RPATH when building shared
IF(BUILD_SHARED_LIBS)
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    print_var(CMAKE_INSTALL_RPATH)
    IF (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
      SET(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib"
        CACHE STRING "On MacOSX the location of installed shared libs.")
      print_var(CMAKE_INSTALL_NAME_DIR)
    ENDIF()
ENDIF()

# Disable EpetraExt HDF5 if we are loading HDF5
SET(EpetraExt_ENABLE_HDF5 OFF CACHE BOOL "Turn off HDF5 in Trilinos.")

# Do all of the processing for this Tribits project
TRIBITS_PROJECT()

# Add the final documentation
add_subdirectory(doc)

##---------------------------------------------------------------------------##
## end of CMakeLists.txt
##---------------------------------------------------------------------------##
