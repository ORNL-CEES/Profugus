##---------------------------------------------------------------------------##
## SPn/CMakeLists.txt
## Thomas M. Evans
## Monday June 20 15:40:52 2011
##---------------------------------------------------------------------------##
## Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## TriBITS package file for SPn package
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

TRIBITS_PACKAGE(SPn)

##---------------------------------------------------------------------------##
## B) Set up package-specific options
##---------------------------------------------------------------------------##

## ML OPTIONS

## Turn on ML if that package is enabled
IF (${PROJECT_NAME}_ENABLE_ML)
  SET(USE_ML 1)
ENDIF()

## MCLS options
IF (${PROJECT_NAME}_ENABLE_MCLS)
   MESSAGE(STATUS "Building Profugus with MCLS support")
   SET(USE_MCLS 1)
ENDIF()

## MueLu
IF (${PROJECT_NAME}_ENABLE_MueLu)
    MESSAGE(STATUS "Building Profugus with MueLu support")
    SET(USE_MUELU 1)
ENDIF()

## Trilinos Timing options.
IF (${PROJECT_NAME}_ENABLE_TRILINOS_TIMING)
   MESSAGE(STATUS "Building Profugus with Trilinos timing")
   SET(USE_TRILINOS_TIMING 1)
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

# mesh sources

FILE(GLOB MESH_HEADERS mesh/*.hh)
SET(MESH_SOURCES
  mesh/LG_Indexer.cc
  mesh/Mesh.cc
  mesh/Partitioner.cc
  )
LIST(APPEND HEADERS ${MESH_HEADERS})
LIST(APPEND SOURCES ${MESH_SOURCES})

# solvers sources

FILE(GLOB SOLVERS_HEADERS solvers/*.hh)
SET(SOLVERS_SOURCES
  solvers/Arnoldi.pt.cc
  solvers/Davidson_Eigensolver.pt.cc
  solvers/Decomposition.cc
  solvers/EigenvalueSolverBuilder.pt.cc
  solvers/InverseOperator.pt.cc
  solvers/LinearSolverBuilder.pt.cc
  solvers/PowerIteration.pt.cc
  solvers/PreconditionerBuilder.pt.cc
  solvers/RayleighQuotient.pt.cc
  solvers/Richardson.pt.cc
  solvers/ShiftedInverseOperator.pt.cc
  solvers/ShiftedOperator.pt.cc
  solvers/StratimikosSolver.pt.cc
  )
LIST(APPEND HEADERS ${SOLVERS_HEADERS})
LIST(APPEND SOURCES ${SOLVERS_SOURCES})

IF( USE_MUELU )
  SET(MUELU_SOURCES
    solvers/MueLuPreconditioner.cc
    solvers/MueLuPreconditioner.pt.cc
  )
  LIST(APPEND SOURCES ${MUELU_SOURCES})
ENDIF()

# spn sources

FILE(GLOB SPN_HEADERS spn/*.hh)
SET(SPN_SOURCES
  spn/Dimensions.cc
  spn/Eigenvalue_Solver.pt.cc
  spn/Energy_Multigrid.pt.cc
  spn/Energy_Prolongation.pt.cc
  spn/Energy_Restriction.pt.cc
  spn/FV_Bnd_Indexer.cc
  spn/FV_Gather.cc
  spn/Fixed_Source_Solver.pt.cc
  spn/Isotropic_Source.cc
  spn/Linear_System.pt.cc
  spn/Linear_System_FV.pt.cc
  spn/Moment_Coefficients.cc
  spn/SDM_Face_Field.cc
  spn/Solver_Base.pt.cc
  spn/SpnSolverBuilder.cc
  spn/Time_Dependent_Solver.pt.cc
  )
LIST(APPEND HEADERS ${SPN_HEADERS})
LIST(APPEND SOURCES ${SPN_SOURCES})

# driver sources

FILE(GLOB DRIVER_HEADERS driver/*.hh)
SET(DRIVER_SOURCES
  driver/Manager.cc
  driver/Problem_Builder.cc
  driver/spn.cc
  )
LIST(APPEND HEADERS ${DRIVER_HEADERS})
LIST(APPEND SOURCES ${DRIVER_SOURCES})

LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/config.h)

##---------------------------------------------------------------------------##
## LIBRARY AND INSTALL TARGETS

TRIBITS_ADD_LIBRARY(
  SPn
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES})

TRIBITS_ADD_EXECUTABLE(
  xspn
  NOEXESUFFIX
  NOEXEPREFIX
  SOURCES driver/spn.cc
  DEPLIBS SPn
  INSTALLABLE
  )

##---------------------------------------------------------------------------##
# Add tests to this package

TRIBITS_ADD_TEST_DIRECTORIES(
  mesh/test
  solvers/test
  spn/test)

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
