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

# Setup debug option
TRIBITS_ADD_DEBUG_OPTION()

# to allow includes like #include "comm/Comm.h"
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..)

# to allow includes like #include "Utils/config.h"
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/..)

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

##---------------------------------------------------------------------------##
## SOURCES AND HEADERS

SET(SOURCES "")
SET(HEADERS "")

# geometry sources

FILE(GLOB GEOMETRY_HEADERS geometry/*.hh)
SET(GEOMETRY_SOURCES
  geometry/Bounding_Box.cc
  geometry/Cartesian_Mesh.cc
  geometry/RTK_Geometry.pt.cc
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

# sim_ce sources

FILE(GLOB SIM_CE_HEADERS sim_ce/*.hh)
SET(SIM_CE_SOURCES
  )
LIST(APPEND HEADERS ${SIM_CE_HEADERS})
LIST(APPEND SOURCES ${SIM_CE_SOURCES})

# mc sources

FILE(GLOB MC_HEADERS mc/*.hh)
SET(MC_SOURCES
  mc/Anderson_Operator.pt.cc
  mc/Anderson_Solver.pt.cc
  mc/Box_Shape.cc
  mc/Cell_Tally.pt.cc
  mc/Domain_Transporter.pt.cc
  mc/Fission_Matrix_Acceleration.pt.cc
  mc/Fission_Matrix_Processor.cc
  mc/Fission_Matrix_Solver.pt.cc
  mc/Fission_Matrix_Tally.pt.cc
  mc/Fission_Rebalance.pt.cc
  mc/Fission_Source.pt.cc
  mc/Fixed_Source_Solver.pt.cc
  mc/Global_RNG.cc
  mc/Group_Bounds.cc
  mc/KCode_Solver.pt.cc
  mc/Keff_Tally.pt.cc
  mc/Particle.pt.cc
  mc/Particle_Metaclass.pt.cc
  mc/Physics.pt.cc
  mc/Sampler.cc
  mc/Solver.pt.cc
  mc/Source.pt.cc
  mc/Source_Diagnostic_Tally.pt.cc
  mc/Source_Transporter.pt.cc
  mc/Tallier.pt.cc
  mc/Tally.pt.cc
  mc/Uniform_Source.pt.cc
  mc/VR_Roulette.pt.cc
  )
LIST(APPEND HEADERS ${MC_HEADERS})
LIST(APPEND SOURCES ${MC_SOURCES})

# driver sources

FILE(GLOB DRIVER_HEADERS mc_driver/*.hh)
SET(DRIVER_SOURCES
  mc_driver/Manager.pt.cc
  mc_driver/Manager_Builder.cc
  mc_driver/Problem_Builder.pt.cc
  mc_driver/Geometry_Builder.pt.cc
  mc_driver/mc.cc
  )
LIST(APPEND HEADERS ${DRIVER_HEADERS})
LIST(APPEND SOURCES ${DRIVER_SOURCES})

LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/config.h)

##---------------------------------------------------------------------------##
## Check for CUDA

SET(USE_CUDA OFF)
IF(${PROJECT_NAME}_ENABLE_CUDA)
  MESSAGE(STATUS "Building ${PACKAGE_NAME} with CUDA support")
  SET(USE_CUDA ON)
ENDIF()

SET(CUDA_HEADERS "")
SET(CUDA_SOURCES "")
IF(USE_CUDA)
  FILE(GLOB CUDA_GEOMETRY_HEADERS cuda_geometry/*.hh)
  SET(CUDA_GEOMETRY_SOURCES
    cuda_geometry/Cartesian_Mesh.cu
    cuda_geometry/Mesh_Geometry.cu
    )
  LIST(APPEND CUDA_HEADERS ${CUDA_GEOMETRY_HEADERS})
  LIST(APPEND CUDA_SOURCES ${CUDA_GEOMETRY_SOURCES})

  FILE(GLOB CUDA_MC_HEADERS cuda_mc/*.hh)
  SET(CUDA_MC_SOURCES
    cuda_mc/Box_Shape.cu
    cuda_mc/Cell_Tally.pt.cu
    cuda_mc/Domain_Transporter.pt.cu
    cuda_mc/Physics.pt.cu
    cuda_mc/Source.pt.cu
    cuda_mc/Source_Transporter.pt.cu
    cuda_mc/Uniform_Source.pt.cu
    cuda_mc/VR_Roulette.pt.cu
    )
  LIST(APPEND CUDA_HEADERS ${CUDA_MC_HEADERS})
  LIST(APPEND CUDA_SOURCES ${CUDA_MC_SOURCES})
ENDIF()

##---------------------------------------------------------------------------##
## LIBRARY AND INSTALL TARGETS

TRIBITS_ADD_LIBRARY(
  MC
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES})

IF(USE_CUDA)
  PACKAGE_ADD_CUDA_LIBRARY(
    MC_cuda
    NOINSTALLHEADERS ${CUDA_HEADERS}
    SOURCES ${CUDA_SOURCES})
ENDIF()

TRIBITS_ADD_EXECUTABLE(
  xmc
  NOEXESUFFIX
  NOEXEPREFIX
  SOURCES mc_driver/mc.cc
  INSTALLABLE
  )

##---------------------------------------------------------------------------##
# Add tests to this package

TRIBITS_ADD_TEST_DIRECTORIES(
  geometry/test
  sim_ce/test
  mc/test)

IF(USE_CUDA)
  TRIBITS_ADD_TEST_DIRECTORIES(
    cuda_geometry/test
    cuda_mc/test)
ENDIF()

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
