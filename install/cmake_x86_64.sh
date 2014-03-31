#!/bin/sh
##---------------------------------------------------------------------------##
## CMAKE FOR X86_64
##---------------------------------------------------------------------------##

# SOURCE AND INSTALL
SOURCE=<SET_SOURCE_DIR>
INSTALL=<SET_INSTALL_DIR>

# BUILD OPTIONS
BUILD="DEBUG"
MPI="ON"

# TPL PATHS
HDF5_PATH="/vendors/hdf5_parallel"
MPI_PATH="/opt/openmpi/gcc/current"

##---------------------------------------------------------------------------##

cmake \
-DCMAKE_BUILD_TYPE:STRING="$BUILD" \
-DTPL_ENABLE_MPI:BOOL=$MPI \
-DCMAKE_INSTALL_PREFIX:PATH=$INSTALL \
\
-DMPI_BASE_DIR:PATH=$MPI_PATH \
\
-DTPL_ENABLE_HDF5:BOOL=ON \
-DHDF5_INCLUDE_DIRS:PATH=$HDF5_PATH/include \
-DHDF5_LIBRARY_DIRS:PATH=$HDF5_PATH/lib \
\
-DBLAS_LIBRARY_DIRS:PATH=/vendors/gcc/atlas/lib \
-DLAPACK_LIBRARY_DIRS:PATH=/vendors/gcc/atlas/lib \
-DBLAS_LIBRARY_NAMES:STRING="f77blas;cblas;atlas" \
-DLAPACK_LIBRARY_NAMES:STRING="lapack" \
\
-DProfugus_CONFIGURE_OPTIONS_FILE:FILEPATH="${SOURCE}/install/base.cmake" \
\
${SOURCE}

##---------------------------------------------------------------------------##
## end of cmake_x86_64.sh
##---------------------------------------------------------------------------##
