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
    XS_Device_Tester.cu
    )

  PACKAGE_ADD_CUDA_LIBRARY(
    Profugus_cuda_xs_test_cuda
    SOURCES ${CUDA_SOURCES}
    DEPLIBS Matprop Matprop_cuda
    TESTONLY)
ENDIF()

##---------------------------------------------------------------------------##
## TESTS
##---------------------------------------------------------------------------##

IF(USE_CUDA)
  SET(DEPS DEPLIBS Profugus_cuda_xs_test_cuda)
ENDIF()

ADD_UTILS_TEST(tstXS_Device.cc              NP 1 ${DEPS} )

##---------------------------------------------------------------------------##
## end of cuda_xs/test/CMakeLists.txt
##---------------------------------------------------------------------------##
