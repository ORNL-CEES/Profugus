##---------------------------------------------------------------------------##
## acc/CMakeLists.txt
## Thomas M. Evans
## Wednesday June 22 11:14:40 2011
##---------------------------------------------------------------------------##
## Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for acc
##---------------------------------------------------------------------------##

TRIBITS_SUBPACKAGE(ACC)

##---------------------------------------------------------------------------##
## TRIBITS_PACKAGE-SPECIFIC OPTIONS
##---------------------------------------------------------------------------##

SET(USE_ACC OFF)
IF(${PROJECT_NAME}_ENABLE_ACC)
  SET(USE_ACC ON)

  SET(ACC_FLAGS "-acc -ta=tesla:kepler+,cuda5.5,nodebug -Minfo=accel")

  # setup compiler flags
  SET(CMAKE_EXE_LINKER_FLAGS
    "${ACC_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")

  MESSAGE(STATUS "Building ${PACKAGE_NAME} with ACC support")
ELSE()

  # suppress pragma warnings
  IF(CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
  ENDIF()
ENDIF()

##---------------------------------------------------------------------------##
## TRIBITS_PACKAGE CONFIGURE FILE
##---------------------------------------------------------------------------##

TRIBITS_CONFIGURE_FILE(config.h)

##---------------------------------------------------------------------------##
## HEADER AND SOURCE FILES
##---------------------------------------------------------------------------##

FILE(GLOB ROOT_HEADERS *.hh)
LIST(APPEND ROOT_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/config.h)

# C/C++ sources
SET(SOURCES
  Particle_Functions.cc
  Physics.cc
  )

# C/C++ sources with ACC directives
SET(ACC_SOURCES
  Change_Direction.acc.cc
  Geometry.acc.cc
  RNG.acc.cc
  Physics.acc.cc
  )

##---------------------------------------------------------------------------##
## TRIBITS_PACKAGE TARGETS AND INSTALL
##---------------------------------------------------------------------------##

# Primary C++ library
TRIBITS_ADD_LIBRARY(
    MC_acc
    SOURCES ${SOURCES} ${ACC_SOURCES})

# Add OpenACC flags
IF(USE_ACC)
  SET_SOURCE_FILES_PROPERTIES(
    ${ACC_SOURCES}
    PROPERTIES
    COMPILE_FLAGS ${ACC_FLAGS}
    )
ENDIF()

#
# Add test directory
#
TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_SUBPACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
##                   end of acc/CMakeLists.txt
##---------------------------------------------------------------------------##
