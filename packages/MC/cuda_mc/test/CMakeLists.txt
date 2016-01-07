##---------------------------------------------------------------------------##
## cuda_xs/test/CMakeLists.txt
## Stuart Slattery
##---------------------------------------------------------------------------##

INCLUDE(UtilsTest)

##---------------------------------------------------------------------------##
## TEST LIBRARIES
##---------------------------------------------------------------------------##

IF(USE_CUDA)
  SET(CUDA_SOURCES
    Particle_Vector_Tester.cu
    )

  PACKAGE_ADD_CUDA_LIBRARY(
    Profugus_cuda_mc_test_cuda
    SOURCES ${CUDA_SOURCES}
    DEPLIBS MC MC_cuda
    TESTONLY)
ENDIF()

##---------------------------------------------------------------------------##
## TESTS
##---------------------------------------------------------------------------##

IF(USE_CUDA)
  SET(DEPS DEPLIBS Profugus_cuda_mc_test_cuda)
ENDIF()

ADD_UTILS_TEST(tstParticle_Vector.cc              NP 1 ${DEPS} )

##---------------------------------------------------------------------------##
## end of cuda_xs/test/CMakeLists.txt
##---------------------------------------------------------------------------##
