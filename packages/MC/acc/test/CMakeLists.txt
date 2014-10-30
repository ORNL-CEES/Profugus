##---------------------------------------------------------------------------##
## acc/test/CMakeLists.txt
## Thomas M. Evans
## Wednesday July 11 14:35:42 2012
##---------------------------------------------------------------------------##
## Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for acc
##---------------------------------------------------------------------------##

INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(UtilsTest)

##---------------------------------------------------------------------------##
## LIBRARIES
##---------------------------------------------------------------------------##

SET(ACC_SOURCES
  ParticleTest.cc
  PhysicsTest.cc
  )

TRIBITS_ADD_LIBRARY(acc_test_lib
  SOURCES ${ACC_SOURCES}
  TESTONLY
  DEPLIBS MC_acc)

##---------------------------------------------------------------------------##
## TESTS
##---------------------------------------------------------------------------##

IF(USE_ACC)
  ADD_UTILS_TEST(tstACC.cc)
ENDIF()

ADD_UTILS_TEST(tstPhysics.cc  NP 1 DEPLIBS acc_test_lib)
ADD_UTILS_TEST(tstParticle.cc NP 1 DEPLIBS acc_test_lib)
ADD_UTILS_TEST(tstGeometry.cc      DEPLIBS acc_test_lib)

##---------------------------------------------------------------------------##
## FILES
##---------------------------------------------------------------------------##

IF(USE_ACC)
  SET_SOURCE_FILES_PROPERTIES(
    ${ACC_SOURCES}
    PROPERTIES
    COMPILE_FLAGS ${ACC_FLAGS}
    )
ENDIF()

##---------------------------------------------------------------------------##
##                    end of acc/test/CMakeLists.txt
##---------------------------------------------------------------------------##
