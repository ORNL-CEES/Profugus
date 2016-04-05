##---------------------------------------------------------------------------##
## cuda_mc/test/CMakeLists.txt
## Steven Hamilton
## Thu Apr 11 09:20:40 2013
##---------------------------------------------------------------------------##
## Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for cuda_mc/test
##---------------------------------------------------------------------------##

INCLUDE(UtilsTest)

##---------------------------------------------------------------------------##
## TEST LIBRARIES
##---------------------------------------------------------------------------##

# C++ code
SET(SOURCES
  )

SET(CUDA_SOURCES
  Box_Shape_Tester.cu
  Cell_Tally_Tester.cu
  Domain_Transporter_Tester.cu
  Fixed_Solver_Tester.cu
  Keff_Tally_Tester.cu
  Particle_Tester.cu
  Physics_Tester.cu
  Source_Transporter_Tester.cu
  Uniform_Source_Tester.cu
)

PACKAGE_ADD_CUDA_LIBRARY(
  mc_cuda_mc_test_cuda
  SOURCES ${CUDA_SOURCES}
  DEPLIBS MC MC_cuda
  TESTONLY)

##---------------------------------------------------------------------------##
## TESTS
##---------------------------------------------------------------------------##
SET(DEPS DEPLIBS mc_cuda_mc_test_cuda)

ADD_UTILS_TEST(tstBox_Shape_cuda.cc             NP 1 ${DEPS})
ADD_UTILS_TEST(tstCell_Tally_cuda.cc            NP 1 ${DEPS})
ADD_UTILS_TEST(tstDomain_Transporter_cuda.cc    NP 1 ${DEPS})
ADD_UTILS_TEST(tstFission_Rebalance_cuda.cc     NP 1  2 4 ${DEPS})
ADD_UTILS_TEST(tstFixed_Source_Solver_cuda.cc   NP 1 ${DEPS})
ADD_UTILS_TEST(tstKeff_Tally_cuda.cc            NP 1 ${DEPS})
ADD_UTILS_TEST(tstParticle_cuda.cc              NP 1 ${DEPS})
ADD_UTILS_TEST(tstPhysics_cuda.cc               NP 1 ${DEPS})
ADD_UTILS_TEST(tstSource_Transporter_cuda.cc    NP 1 ${DEPS})
ADD_UTILS_TEST(tstUniform_Source_cuda.cc        NP 1 ${DEPS})

##---------------------------------------------------------------------------##
##                        end of cuda_utils/test/CMakeLists.txt
##---------------------------------------------------------------------------##
