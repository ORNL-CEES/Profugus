##---------------------------------------------------------------------------##
## MC/CMakeLists.txt
## Thomas M. Evans
## Monday June 20 15:40:52 2011
##---------------------------------------------------------------------------##
## Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## TriBITS package file for MC package
##---------------------------------------------------------------------------##

INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TribitsAddExecutable)
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

# Setup debug option
TRIBITS_ADD_DEBUG_OPTION()

# to allow includes like #include "comm/Comm.h"
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# to allow includes like #include "Utils/config.h"
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/..)

IF(CMAKE_COMPILER_IS_GNUCXX AND NOT WIN32)
  PROFUGUS_ADD_CXX_FLAGS(
    -Wno-sign-compare
    -Wno-unused-but-set-variable
    -Wno-unused-variable)
ENDIF()

IF(${PROJECT_NAME}_ENABLE_ACC)
  SET(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -acc -ta=tesla:kepler+,cuda5.5 -Minfo=accel")
  SET(CMAKE_EXE_LINKER_FLAGS
    "-acc -ta=tesla:kepler+,cuda5.5 ${CMAKE_EXE_LINKER_FLAGS}")
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

##---------------------------------------------------------------------------##
## SOURCES AND HEADERS

SET(SOURCES "")
SET(HEADERS "")

# geometry sources

FILE(GLOB GEOMETRY_HEADERS mesh/*.hh)
SET(GEOMETRY_SOURCES
  geometry/Cartesian_Mesh.cc
  geometry/Geometry.pt.cc
  geometry/Mesh_Geometry.cc
  geometry/Mesh_State.cc
  geometry/RTK_Array.cc
  geometry/RTK_Array.pt.cc
  geometry/RTK_Cell.cc
  geometry/RTK_Functions.cc
  geometry/RTK_State.cc
  )
LIST(APPEND HEADERS ${GEOMETRY_HEADERS})
LIST(APPEND SOURCES ${GEOMETRY_SOURCES})

# mc sources

FILE(GLOB MC_HEADERS mesh/*.hh)
SET(MC_SOURCES
  mc/Box_Shape.cc
  mc/Domain_Transporter.cc
  mc/Fission_Rebalance.cc
  mc/Fission_Source.cc
  mc/Fixed_Source_Solver.cc
  mc/Global_RNG.cc
  mc/Group_Bounds.cc
  mc/KCode_Solver.cc
  mc/Keff_Tally.cc
  mc/Particle.cc
  mc/Physics.cc
  mc/Sampler.cc
  mc/Solver.cc
  mc/Source.cc
  mc/Source_Transporter.cc
  mc/Tallier.cc
  mc/Tally.cc
  mc/Uniform_Source.cc
  mc/VR_Roulette.cc
  )
LIST(APPEND HEADERS ${MC_HEADERS})
LIST(APPEND SOURCES ${MC_SOURCES})

# driver sources

FILE(GLOB DRIVER_HEADERS driver/*.hh)
SET(DRIVER_SOURCES
  driver/Manager.cc
  driver/Problem_Builder.cc
  driver/mc.cc
  )
LIST(APPEND HEADERS ${DRIVER_HEADERS})
LIST(APPEND SOURCES ${DRIVER_SOURCES})

LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/config.h)

##---------------------------------------------------------------------------##
## LIBRARY AND INSTALL TARGETS

TRIBITS_ADD_LIBRARY(
  MC
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES})

TRIBITS_ADD_EXECUTABLE(
  xmc
  NOEXESUFFIX
  NOEXEPREFIX
  SOURCES driver/mc.cc
  INSTALLABLE
  )

##---------------------------------------------------------------------------##
# Add tests to this package

TRIBITS_ADD_TEST_DIRECTORIES(
  geometry/test
  mc/test)

##---------------------------------------------------------------------------##
## FINISH SETUP

TRIBITS_PROCESS_SUBPACKAGES()

# add the main sub-packages to the package
TRIBITS_PACKAGE_DEF()

# Exclude files for source package.
TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

##---------------------------------------------------------------------------##
## D) Do standard postprocessing
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
##                      end of Utils/CMakeLists.txt
##---------------------------------------------------------------------------##
