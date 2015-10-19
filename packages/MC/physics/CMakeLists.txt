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

FILE(GLOB PHYSICS_HEADERS *.hh)
SET(PHYSICS_SOURCES
  Box_Shape.cc
  Cell_Tally.cc
  Domain_Transporter.cc
  Fixed_Source_Solver.cc
  Global_RNG.cc
  Group_Bounds.cc
  Particle.cc
  Particle_Metaclass.pt.cc
  Physics.cc
  Sampler.cc
  Solver.cc
  Source.cc
  Source_Transporter.cc
  Tallier.cc
  Tally.cc
  Uniform_Source.cc
  VR_Roulette.cc
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
