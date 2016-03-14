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
    Tally_Tester.cu
    Particle_Vector_Tester.cu
    Physics_Tester.cu
    Uniform_Source_Tester.cu
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

ADD_UTILS_TEST(tstParticle_Vector_cuda.cc              NP 1 ${DEPS} )
ADD_UTILS_TEST(tstCollision_Tally_cuda.cc              NP 1 ${DEPS} )
ADD_UTILS_TEST(tstCell_Tally_cuda.cc                   NP 1 ${DEPS} )
ADD_UTILS_TEST(tstUniform_Source_cuda.cc               NP 1 ${DEPS} )
ADD_UTILS_TEST(tstBank_cuda.cc                         NP 1 ${DEPS} )
ADD_UTILS_TEST(tstPhysics_cuda.cc                      NP 1 ${DEPS} )
ADD_UTILS_TEST(tstKeff_Tally_cuda.cc                   NP 1 ${DEPS} )
ADD_UTILS_TEST(tstVR_Roulette_cuda.cc                  NP 1 ${DEPS} )
ADD_UTILS_TEST(tstFission_Rebalance.cc                 NP 1  )

##---------------------------------------------------------------------------##
## end of cuda_xs/test/CMakeLists.txt
##---------------------------------------------------------------------------##
