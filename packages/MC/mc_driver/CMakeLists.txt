##---------------------------------------------------------------------------##
## MC/mc_driver/CMakeLists.txt
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

TRIBITS_SUBPACKAGE(mc_driver)

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

# driver sources

FILE(GLOB DRIVER_HEADERS driver/*.hh)
SET(DRIVER_SOURCES
  src/Manager.cc
  src/Problem_Builder.cc
  src/mc.cc
  )
LIST(APPEND HEADERS ${DRIVER_HEADERS})
LIST(APPEND SOURCES ${DRIVER_SOURCES})

##---------------------------------------------------------------------------##
## LIBRARY AND INSTALL TARGETS

TRIBITS_ADD_LIBRARY(
  mc_driver
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES})

TRIBITS_ADD_EXECUTABLE(
  xmc
  NOEXESUFFIX
  NOEXEPREFIX
  SOURCES src/mc.cc
  INSTALLABLE
  )

##---------------------------------------------------------------------------##
## D) Do standard postprocessing
##---------------------------------------------------------------------------##

TRIBITS_SUBPACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
##                      end of Utils/CMakeLists.txt
##---------------------------------------------------------------------------##
