##---------------------------------------------------------------------------##
## Multiphysics/CMakeLists.txt
## Steven Hamilton
##---------------------------------------------------------------------------##
## Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## TriBITS package file for Multiphysics package
##---------------------------------------------------------------------------##

INCLUDE(AppendSet)
INCLUDE(ProfugusAddCFlags)

# Add Utils cmake to path
APPEND_SET(CMAKE_MODULE_PATH ${Utils_SOURCE_DIR}/cmake)

##---------------------------------------------------------------------------##
## A) Define the package
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE(Multiphysics)

##---------------------------------------------------------------------------##
## B) Set up package-specific options
##---------------------------------------------------------------------------##

# Setup debug option
TRIBITS_ADD_DEBUG_OPTION()

# to allow includes like #include "comm/Comm.h"
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

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


##---------------------------------------------------------------------------##
## SOURCES AND HEADERS

SET(SOURCES "")
SET(HEADERS "")

# surrogate sources

FILE(GLOB MC_HEADERS mesh/*.hh)
SET(SURROGATE_SOURCES
  surrogate/Assembly_Model.cc
  surrogate/Multi_Pin_Conduction.cc
  surrogate/Multi_Pin_Subchannel.cc
  surrogate/Multiphysics_Driver.cc
  surrogate/Single_Pin_Conduction.cc
  surrogate/Single_Pin_Subchannel.cc
  surrogate/Two_Group_Diffusion.cc 
  surrogate/Two_Group_Cross_Sections.cc
  surrogate/Water_Properties.cc
  )
LIST(APPEND HEADERS ${SURROGATE_HEADERS})
LIST(APPEND SOURCES ${SURROGATE_SOURCES})

##---------------------------------------------------------------------------##
## LIBRARY AND INSTALL TARGETS

TRIBITS_ADD_LIBRARY(
  Surrogate
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES})

##---------------------------------------------------------------------------##
# Add tests to this package

TRIBITS_ADD_TEST_DIRECTORIES(
  surrogate/test)

##---------------------------------------------------------------------------##
## FINISH SETUP

# Exclude files for source package.
TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

##---------------------------------------------------------------------------##
## D) Do standard postprocessing
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
##                      end of Multiphysics/CMakeLists.txt
##---------------------------------------------------------------------------##
