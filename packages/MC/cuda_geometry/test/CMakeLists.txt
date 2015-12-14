##---------------------------------------------------------------------------##
## cuda_geometry/test/CMakeLists.txt
## Seth R Johnson
## Thu Apr 11 09:20:40 2013
##---------------------------------------------------------------------------##
## Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for CUDA_utils
##---------------------------------------------------------------------------##

INCLUDE(UtilsTest)

##---------------------------------------------------------------------------##
## TEST LIBRARIES
##---------------------------------------------------------------------------##

# C++ code
SET(SOURCES
  )

#TRIBITS_ADD_LIBRARY(
#  mc_cuda_geometry_test
#  SOURCES ${SOURCES}
#  TESTONLY
#  DEPLIBS MC)

SET(CUDA_SOURCES
  Cartesian_Mesh_Tester.cu
)

PACKAGE_ADD_CUDA_LIBRARY(
  mc_cuda_geometry_test_cuda
  SOURCES ${CUDA_SOURCES}
  DEPLIBS MC MC_cuda
  TESTONLY)

##---------------------------------------------------------------------------##
## TESTS
##---------------------------------------------------------------------------##
SET(DEPS DEPLIBS mc_cuda_geometry_test_cuda)
#SET(DEPS ${DEPS} mc_cuda_geometry_test)

ADD_UTILS_TEST(tstCartesian_Mesh_cuda.cc        NP 1 ${DEPS})

##---------------------------------------------------------------------------##
##                        end of cuda_utils/test/CMakeLists.txt
##---------------------------------------------------------------------------##
