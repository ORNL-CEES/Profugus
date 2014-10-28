##---------------------------------------------------------------------------##
## MC/core/CMakeLists.txt
## Thomas M. Evans
## Monday June 20 15:40:52 2011
##---------------------------------------------------------------------------##
## Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## TriBITS package file for MC core package
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## Tribits Package name
##---------------------------------------------------------------------------##

TRIBITS_SUBPACKAGE(Core)

##---------------------------------------------------------------------------##
## TRIBITS_PACKAGE CONFIGURE FILE
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
  MC_core
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

TRIBITS_SUBPACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
##                      end of MC/core/CMakeLists.txt
##---------------------------------------------------------------------------##
