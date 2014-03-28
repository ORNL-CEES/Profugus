##---------------------------------------------------------------------------##
## CMAKE REGRESSION TPLs FILE
##---------------------------------------------------------------------------##

# TPLS
SET(TPL_ENABLE_SILO ON CACHE BOOL "")
SET(SILO_INCLUDE_DIRS /vendors/gcc/silo/include CACHE PATH "")
SET(SILO_LIBRARY_DIRS /vendors/gcc/silo/lib     CACHE PATH "")

SET(TPL_ENABLE_HDF5 ON CACHE BOOL "")
SET(HDF5_INCLUDE_DIRS /vendors/hdf5/include CACHE PATH "")
SET(HDF5_LIBRARY_DIRS /vendors/hdf5/lib     CACHE PATH "")

SET(TPL_ENABLE_LAVA ON CACHE BOOL "")
SET(LAVA_INCLUDE_DIRS /vendors/lava/include CACHE PATH "")
SET(LAVA_LIBRARY_DIRS /vendors/lava/lib     CACHE PATH "")
SET(MCNP_EXECUTABLE /opt/advantg/mcnp/bin/linux/mcnp5_linux_x86_64_omp CACHE FILEPATH "")

SET(TPL_ENABLE_KGTLIB ON CACHE BOOL "")
SET(KGTLIB_INCLUDE_DIRS /vendors/kgtlib/debug/include CACHE PATH "")
SET(KGTLIB_LIBRARY_DIRS /vendors/kgtlib/debug/lib     CACHE PATH "")

SET(BLAS_LIBRARY_DIRS    /vendors/gcc/atlas/lib CACHE PATH "")
SET(LAPACK_LIBRARY_DIRS  /vendors/gcc/atlas/lib CACHE PATH "")
SET(BLAS_LIBRARY_NAMES   "f77blas;cblas;atlas"  CACHE STRING "")
SET(LAPACK_LIBRARY_NAMES "lapack"               CACHE STRING "")
