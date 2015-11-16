##---------------------------------------------------------------------------##
## cuda_utils/CMakeLists.txt
## Seth R Johnson
## Thu Apr 11 09:09:00 2013
##---------------------------------------------------------------------------##
## Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for CUDA_utils
##---------------------------------------------------------------------------##


TRIBITS_SUBPACKAGE(CUDA_utils)

##---------------------------------------------------------------------------##
## PACKAGE-SPECIFIC OPTIONS
##---------------------------------------------------------------------------##

SET(USE_CUDA OFF)
IF(${PROJECT_NAME}_ENABLE_CUDA)
  MESSAGE(STATUS "Building ${PACKAGE_NAME} with CUDA support")
  SET(USE_CUDA ON)
ENDIF()

##---------------------------------------------------------------------------##
## PACKAGE CONFIGURE FILE
##---------------------------------------------------------------------------##

TRIBITS_CONFIGURE_FILE(config.h)

##---------------------------------------------------------------------------##
## INCLUDES
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## HEADER AND SOURCE FILES
##---------------------------------------------------------------------------##

# C/C++ headers
FILE(GLOB HEADERS *.hh)

SET(SOURCES
  Atomic_Lock.cc
  BLAS_host.cc
  Device_Vector_host.pt.cc
  Event.cc
  Functions.cc
  Hardware.cc
  Host_Vector.pt.cc
  Multi_Vector.pt.cc
  Stream.cc
  Texture_Vector_host.pt.cc
  )

IF(USE_CUDA)
  LIST(APPEND SOURCES
    BLAS_device.cc
    BLAS_Handle.cc
    CudaDBC.cc
    Device_Vector_device.pt.cc
    Texture_Vector_device.pt.cc
    )
ENDIF()

##---------------------------------------------------------------------------##
## PACKAGE TARGETS AND INSTALL
##---------------------------------------------------------------------------##

TRIBITS_ADD_LIBRARY(
  Profugus_cuda_utils
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES})

#
# Add test directory for this package
#
TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_PACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
##                   end of cuda_utils/CMakeLists.txt
##---------------------------------------------------------------------------##
