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

SET(CUDA_SOURCES
  RTK_State_Tester.cu
  RTK_Cell_Tester_1.cu
  RTK_Cell_Tester_2.cu
  RTK_Array_Tester_1.cu
  RTK_Array_Tester_2.cu
  RTK_Array_Tester_3.cu
)

PACKAGE_ADD_CUDA_LIBRARY(
  mc_cuda_rtk_test_cuda
  SOURCES ${CUDA_SOURCES}
  DEPLIBS MC MC_cuda
  TESTONLY)

##---------------------------------------------------------------------------##
## TESTS
##---------------------------------------------------------------------------##

SET(DEPS DEPLIBS mc_cuda_rtk_test_cuda)

ADD_UTILS_TEST(tstRTK_State_cuda.cc NP 1 ${DEPS})
ADD_UTILS_TEST(tstRTK_Cell_cuda.cc  NP 1 ${DEPS})
ADD_UTILS_TEST(tstRTK_Array_cuda.cc NP 1 ${DEPS})

##---------------------------------------------------------------------------##
##                        end of cuda_utils/test/CMakeLists.txt
##---------------------------------------------------------------------------##
