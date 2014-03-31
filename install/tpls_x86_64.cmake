##---------------------------------------------------------------------------##
## CMAKE REGRESSION TPLs FILE
##---------------------------------------------------------------------------##

# MPI BASE DIR
SET(MPI_BASE_DIR /opt/openmpi/gcc/current CACHE PATH "")

# TPLS
SET(TPL_ENABLE_HDF5 ON CACHE BOOL "")
SET(HDF5_INCLUDE_DIRS /vendors/hdf5_parallel/include CACHE PATH "")
SET(HDF5_LIBRARY_DIRS /vendors/hdf5_parallel/lib     CACHE PATH "")

SET(BLAS_LIBRARY_DIRS    /vendors/gcc/atlas/lib CACHE PATH "")
SET(LAPACK_LIBRARY_DIRS  /vendors/gcc/atlas/lib CACHE PATH "")
SET(BLAS_LIBRARY_NAMES   "f77blas;cblas;atlas"  CACHE STRING "")
SET(LAPACK_LIBRARY_NAMES "lapack"               CACHE STRING "")
