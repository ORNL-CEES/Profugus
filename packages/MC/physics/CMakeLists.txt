##---------------------------------------------------------------------------##
## MC/phyiscs/CMakeLists.txt
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

TRIBITS_SUBPACKAGE(mc_physics)

##---------------------------------------------------------------------------##
## B) Set up package-specific options
##---------------------------------------------------------------------------##

# to allow includes like #include "comm/Comm.h"
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..)

# to allow includes like #include "Utils/config.h"
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/..)

##---------------------------------------------------------------------------##
# C) Process the subpackages
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## SOURCES AND HEADERS

SET(SOURCES "")
SET(HEADERS "")

# physics sources

FILE(GLOB PHYSICS_HEADERS src/*.hh)
SET(PHYSICS_SOURCES
  src/Anderson_Operator.pt.cc
  src/Anderson_Solver.pt.cc
  src/Box_Shape.cc
  src/Cell_Tally.cc
  src/Domain_Transporter.cc
  src/Fission_Matrix_Acceleration.pt.cc
  src/Fission_Matrix_Processor.cc
  src/Fission_Matrix_Solver.pt.cc
  src/Fission_Matrix_Tally.cc
  src/Fission_Rebalance.cc
  src/Fission_Source.cc
  src/Fixed_Source_Solver.cc
  src/Global_RNG.cc
  src/Group_Bounds.cc
  src/KCode_Solver.cc
  src/Keff_Tally.cc
  src/Particle.cc
  src/Particle_Metaclass.pt.cc
  src/Physics.cc
  src/Sampler.cc
  src/Solver.cc
  src/Source.cc
  src/Source_Diagnostic_Tally.cc
  src/Source_Transporter.cc
  src/Tallier.cc
  src/Tally.cc
  src/Uniform_Source.cc
  src/VR_Roulette.cc
  )
LIST(APPEND HEADERS ${PHYSICS_HEADERS})
LIST(APPEND SOURCES ${PHYSICS_SOURCES})

##---------------------------------------------------------------------------##
## LIBRARY AND INSTALL TARGETS

TRIBITS_ADD_LIBRARY(
  mc_physics
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES})

##---------------------------------------------------------------------------##
# Add tests to this package
TRIBITS_ADD_TEST_DIRECTORIES( test )

##---------------------------------------------------------------------------##
## D) Do standard postprocessing
##---------------------------------------------------------------------------##

TRIBITS_SUBPACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
##                      end of Utils/CMakeLists.txt
##---------------------------------------------------------------------------##
