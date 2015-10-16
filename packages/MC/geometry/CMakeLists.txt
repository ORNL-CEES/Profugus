##---------------------------------------------------------------------------##
## MC/geometry/CMakeLists.txt
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

TRIBITS_SUBPACKAGE(mc_geometry)

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

# geometry sources

FILE(GLOB GEOMETRY_HEADERS *.hh)
SET(GEOMETRY_SOURCES
  Cartesian_Mesh.cc
  Geometry.pt.cc
  Mesh_Geometry.cc
  Mesh_State.cc
  RTK_Array.cc
  RTK_Array.pt.cc
  RTK_Cell.cc
  RTK_Functions.cc
  RTK_State.cc
  )
LIST(APPEND HEADERS ${GEOMETRY_HEADERS})
LIST(APPEND SOURCES ${GEOMETRY_SOURCES})

##---------------------------------------------------------------------------##
## LIBRARY AND INSTALL TARGETS

TRIBITS_ADD_LIBRARY(
  mc_geometry
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
