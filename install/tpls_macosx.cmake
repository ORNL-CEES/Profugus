##---------------------------------------------------------------------------##
## CMAKE TPLs FILE
##---------------------------------------------------------------------------##

# MPI BASE DIR
SET(MPI_BASE_DIR /opt/mpi/current CACHE PATH "")

# TPLS
SET(TPL_ENABLE_HDF5 ON CACHE BOOL "")
SET(HDF5_INCLUDE_DIRS /vendors/parallel/hdf5/include CACHE PATH "")
SET(HDF5_LIBRARY_DIRS /vendors/parallel/hdf5/lib     CACHE PATH "")
