##---------------------------------------------------------------------------##
## Basic opt development build
##---------------------------------------------------------------------------##

SET(TPL_ENABLE_MPI ON CACHE BOOL "")
SET(MPI_BASE_DIR /opt/mpi/current CACHE PATH "")
SET(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "")
SET(MPI_EXEC_MAX_NUMPROCS 12 CACHE STRING "Max. num procs for testing")

INCLUDE(${CMAKE_CURRENT_LIST_DIR}/base.cmake)