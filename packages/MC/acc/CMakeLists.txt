##---------------------------------------------------------------------------##
## kba/CMakeLists.txt
## Thomas M. Evans
## Wednesday June 22 11:14:40 2011
##---------------------------------------------------------------------------##
## Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for kba
##---------------------------------------------------------------------------##


TRIBITS_SUBPACKAGE(KBA_CUDA)

##---------------------------------------------------------------------------##
## TRIBITS_PACKAGE-SPECIFIC OPTIONS
##---------------------------------------------------------------------------##

# Whether to use MPI with CUDA
SET(KBA_CUDA_USE_MPI ${TPL_ENABLE_MPI})

# We pass this variable to config.h because CUDA on-device assertions and
# some of the memory pointer functionality don't work on Mac systems.
IF(APPLE)
  SET(SYSTEM_IS_MACOSX ON)
ENDIF()

SET(USE_CUDA OFF)
IF(${PACKAGE_NAME}_ENABLE_CUDA)
  MESSAGE(STATUS "Building ${SUBPACKAGE_NAME} with CUDA support")
  SET(USE_CUDA ON)
ENDIF()

##---------------------------------------------------------------------------##
## TRIBITS_PACKAGE CONFIGURE FILE
##---------------------------------------------------------------------------##

TRIBITS_CONFIGURE_FILE(config.h)

##---------------------------------------------------------------------------##
## HEADER AND SOURCE FILES
##---------------------------------------------------------------------------##

FILE(GLOB ROOT_HEADERS *.hh)
LIST(APPEND ROOT_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/config.h)

# C/C++ sources for when CUDA is disabled
SET(NO_CUDA_SOURCES
  Cuda_Interface.cc
  Cuda_Sweeper_Imp.pt.cc
  LD_Equations.cc
  LD_Equations.pt.cc
  SC_Equations.cc
  SC_Equations.pt.cc
  Mpi_Interface.cc
  Parallel_Environment.cc
  Partition_Manager.cc
  Sweeper_Choices.cc
  Tuple.cc
  Utilities.cc
  )

# CUDA implementation
SET(CUDA_SOURCES
  Tuple.cu
  Utilities.cu
  Cuda_Interface.cu
  Mpi_Interface.cu
  Partition_Manager.cu
  Parallel_Environment.cu
  LD_Equations.cu
  LD_Equations.pt.cu
  SC_Equations.cu
  SC_Equations.pt.cu
  Sweeper_Choices.cu
  Cuda_Sweeper_Imp.pt.cu
  )

SET(SOURCES
  Cuda_Sweeper.pt.cc
  )

##---------------------------------------------------------------------------##
## TRIBITS_PACKAGE TARGETS AND INSTALL
##---------------------------------------------------------------------------##
IF(TPL_ENABLE_CUDA AND TPL_ENABLE_MPI)
  # Because the KBA CUDA code includes lots of files in Exnihilo but uses the
  # nvcc wrapper rather than the mpicxx wrapper, MPI headers aren't included by
  # default, causing stuff to fail. We manually add the include directories
  # here.
  INCLUDE_DIRECTORIES(${MPI_BASE_DIR}/include)
ENDIF()

# CUDA sweeper or C-only version
IF(TPL_ENABLE_CUDA)
  PACKAGE_ADD_CUDA_LIBRARY(
      Denovo_kba_cuda_sweeper
      SOURCES ${CUDA_SOURCES})
ELSE()
  TRIBITS_ADD_LIBRARY(
      Denovo_kba_cuda_sweeper
      SOURCES ${NO_CUDA_SOURCES})
ENDIF(TPL_ENABLE_CUDA)

# Primary C++ library
TRIBITS_ADD_LIBRARY(
    Denovo_kba_cuda
    SOURCES ${SOURCES}
    DEPLIBS Denovo_kba_cuda_sweeper)

INSTALL(FILES
  ${HEADERS}
  DESTINATION
  include/Denovo/kba_cuda)

#
# Add documentation for this package
#
DENOVO_ADD_DOC()

#
# Add test directory
#
TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_SUBPACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
##                   end of kba/CMakeLists.txt
##---------------------------------------------------------------------------##
