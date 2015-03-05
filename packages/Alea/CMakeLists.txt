##---------------------------------------------------------------------------##
## Alea/CMakeLists.txt
## Thomas M. Evans
## Monday June 20 15:40:52 2011
##---------------------------------------------------------------------------##
## Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## TriBITS package file for SPn package
##---------------------------------------------------------------------------##

INCLUDE(AppendSet)
INCLUDE(ProfugusAddCFlags)

# Add Utils cmake to path
APPEND_SET(CMAKE_MODULE_PATH ${Utils_SOURCE_DIR}/cmake)

##---------------------------------------------------------------------------##
## A) Define the package
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE(Alea)

##---------------------------------------------------------------------------##
## B) Set up package-specific options
##---------------------------------------------------------------------------##

## MCLS options
IF (${PROJECT_NAME}_ENABLE_MCLS)
   MESSAGE(STATUS "Building Profugus with MCLS support")
   SET(USE_MCLS 1)
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

# mc_solvers sources

FILE(GLOB MC_SOLVERS_HEADERS mc_solvers/*.hh)
SET(MC_SOLVERS_SOURCES
  mc_solvers/AdditiveSchwarzWrapper.cc
  mc_solvers/AleaSolver.cc
  mc_solvers/BelosSolver.cc
  mc_solvers/ChebyshevIteration.cc
  mc_solvers/ChebyshevPolynomial.cc
  mc_solvers/GmresPolynomial.cc
  mc_solvers/LinearSolverFactory.cc
  mc_solvers/LinearSystem.cc
  mc_solvers/LinearSystemFactory.cc
  mc_solvers/MC_Data.cc
  mc_solvers/MonteCarloSolver.cc
  mc_solvers/NeumannPolynomial.cc
  mc_solvers/Polynomial.cc
  mc_solvers/PolynomialBasis.cc
  mc_solvers/PolynomialFactory.cc
  mc_solvers/PolynomialPreconditioner.cc
  mc_solvers/PolynomialUtils.cc
  mc_solvers/RichardsonIteration.cc
  mc_solvers/SyntheticAcceleration.cc
  )
LIST(APPEND HEADERS ${MC_SOLVERS_HEADERS})
LIST(APPEND SOURCES ${MC_SOLVERS_SOURCES})
INSTALL(FILES ${MC_SOLVERS_HEADERS} DESTINATION include/Alea/mc_solvers)

LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/config.h)
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/config.h DESTINATION include/Alea)

##---------------------------------------------------------------------------##
## LIBRARY AND INSTALL TARGETS

TRIBITS_ADD_LIBRARY(
  Alea
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES})

TRIBITS_ADD_EXECUTABLE(
  xalea
  NOEXESUFFIX
  NOEXEPREFIX
  SOURCES mc_solvers/alea.cc
  INSTALLABLE
  )

TRIBITS_ADD_EXECUTABLE(
  xaleae_multae
  NOEXESUFFIX
  NOEXEPREFIX
  SOURCES mc_solvers/aleae_multae.cc
  INSTALLABLE
  )

# MCLS executable
IF (${PROJECT_NAME}_ENABLE_MCLS)
  TRIBITS_ADD_EXECUTABLE(
    xalea_mcls
    NOEXESUFFIX
    NOEXEPREFIX
    SOURCES mc_solvers/alea_mcls.cc
    INSTALLABLE
    )
ENDIF()

##---------------------------------------------------------------------------##
# Add tests to this package

TRIBITS_ADD_TEST_DIRECTORIES(
  mc_solvers/test)

##---------------------------------------------------------------------------##
## FINISH SETUP

# Exclude files for source package.
TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

##---------------------------------------------------------------------------##
## D) Do standard postprocessing
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
##                      end of Alea/CMakeLists.txt
##---------------------------------------------------------------------------##
