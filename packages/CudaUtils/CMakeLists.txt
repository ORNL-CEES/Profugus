##---------------------------------------------------------------------------##
## cuda_utils/CMakeLists.txt
## Seth R Johnson
## Thu Apr 11 09:09:00 2013
##---------------------------------------------------------------------------##
## Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for CUDA_utils
##---------------------------------------------------------------------------##

INCLUDE(AppendSet)
INCLUDE(ProfugusAddCFlags)

APPEND_SET(CMAKE_MODULE_PATH ${Utils_SOURCE_DIR}/cmake)

##---------------------------------------------------------------------------##
## Define the package
##---------------------------------------------------------------------------##
TRIBITS_PACKAGE(CudaUtils)

##---------------------------------------------------------------------------##
## PACKAGE-SPECIFIC OPTIONS
##---------------------------------------------------------------------------##

SET(USE_CUDA OFF)
IF(${PROJECT_NAME}_ENABLE_CUDA)
  MESSAGE(STATUS "Building ${PROJECT_NAME} with CUDA support")
  SET(USE_CUDA ON)
ENDIF()

##---------------------------------------------------------------------------##
## PACKAGE CONFIGURE FILE
##---------------------------------------------------------------------------##

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

TRIBITS_CONFIGURE_FILE(config.h)

##---------------------------------------------------------------------------##
## INCLUDES
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## HEADER AND SOURCE FILES
##---------------------------------------------------------------------------##

# C/C++ headers
FILE(GLOB HEADERS cuda_utils/*.hh)

SET(SOURCES
  cuda_utils/Atomic_Lock.cc
  cuda_utils/BLAS_host.cc
  cuda_utils/Device_Vector_host.pt.cc
  cuda_utils/Event.cc
  cuda_utils/Functions.cc
  cuda_utils/Hardware.cc
  cuda_utils/Host_Vector.pt.cc
  cuda_utils/Multi_Vector.pt.cc
  cuda_utils/Stream.cc
  cuda_utils/Texture_Vector_host.pt.cc
  )

IF(USE_CUDA)
  LIST(APPEND SOURCES
    cuda_utils/BLAS_device.cc
    cuda_utils/BLAS_Handle.cc
    cuda_utils/CudaDBC.cc
    cuda_utils/Device_Vector_device.pt.cc
    cuda_utils/Texture_Vector_device.pt.cc
    )
ENDIF()

LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/config.h)

##---------------------------------------------------------------------------##
## PACKAGE TARGETS AND INSTALL
##---------------------------------------------------------------------------##

TRIBITS_ADD_LIBRARY(
  Profugus_cuda_utils
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES})

IF(USE_CUDA)
  SET(CUDA_SOURCES
    cuda_utils/SerialDenseDeviceMatrix.cu
    cuda_utils/SerialDenseDeviceVector.cu
    )

  PACKAGE_ADD_CUDA_LIBRARY(
    Profugus_cuda_utils_cuda
    SOURCES ${CUDA_SOURCES}
    DEPLIBS Profugus_cuda_utils
    )
ENDIF()

#
# Add test directory for this package
#
TRIBITS_ADD_TEST_DIRECTORIES( test )

# postprocess
TRIBITS_PACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
##                   end of cuda_utils/CMakeLists.txt
##---------------------------------------------------------------------------##
