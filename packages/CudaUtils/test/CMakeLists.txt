##---------------------------------------------------------------------------##
## cuda_utils/test/CMakeLists.txt
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
  Atomic_Add_Kernel.cc
  Atomic_Lock_Test_Kernel.cc
  Run_Launch_Args.pt.cc
  Polyglot_Kernel.cc
  Profiler_Kernel.cc
  Stream_Test_Kernel.cc
  Texture_Vector_Test_Kernel.cc
  )

TRIBITS_ADD_LIBRARY(
  Profugus_cuda_utils_test
  SOURCES ${SOURCES}
  TESTONLY
  DEPLIBS Profugus_cuda_utils)

# CUDA-compiled code
IF(USE_CUDA)
  SET(CUDA_SOURCES
    Atomic_Add_Kernel.cu
    Atomic_Lock_Test_Kernel.cu
    Run_Launch_Args.pt.cu
    Polyglot_Kernel.cu
    Profiler_Kernel.cu
    Shared_Device_Ptr_Tester.cu
    Stream_Test_Kernel.cu
    Texture_Vector_Test_Kernel.cu
    )

  PACKAGE_ADD_CUDA_LIBRARY(
    Profugus_cuda_utils_test_cuda
    SOURCES ${CUDA_SOURCES}
    DEPLIBS Profugus_cuda_utils
    TESTONLY)
ENDIF()

##---------------------------------------------------------------------------##
## TESTS
##---------------------------------------------------------------------------##
SET(DEPS  DEPLIBS Profugus_cuda_utils_test)
IF(USE_CUDA)
  SET(DEPS ${DEPS} Profugus_cuda_utils_test_cuda)
ENDIF()

ADD_UTILS_TEST(tstAtomic_Add.cc           NP 1 ${DEPS})
ADD_UTILS_TEST(tstAtomic_Lock.cc          NP 1 ${DEPS})
ADD_UTILS_TEST(tstBLAS.cc                 NP 1        )
ADD_UTILS_TEST(tstDevice_Vector.cc        NP 1 ${DEPS})
ADD_UTILS_TEST(tstEvent.cc                NP 1 ${DEPS})
ADD_UTILS_TEST(tstHardware.cc             NP 1        )
ADD_UTILS_TEST(tstHost_Vector.cc          NP 1 ${DEPS})
ADD_UTILS_TEST(tstLaunch_Args.cc          NP 1 ${DEPS})
ADD_UTILS_TEST(tstMulti_Vector.cc         NP 1 ${DEPS})
ADD_UTILS_TEST(tstPseudo_Cuda.cc          NP 1 ${DEPS})
ADD_UTILS_TEST(tstPseudo_Cuda_Polyglot.cc NP 1 ${DEPS})
ADD_UTILS_TEST(tstShared_Device_Ptr.cc    NP 1 ${DEPS})
ADD_UTILS_TEST(tstStream.cc               NP 1 ${DEPS})
ADD_UTILS_TEST(tstTexture_Vector.cc       NP 1 ${DEPS})

IF(USE_CUDA)
  ADD_UTILS_TEST(tstBLAS_Handle.cc        NP 1        )
ENDIF()

##---------------------------------------------------------------------------##
##                        end of cuda_utils/test/CMakeLists.txt
##---------------------------------------------------------------------------##
