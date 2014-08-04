#!/bin/sh
##---------------------------------------------------------------------------##
## CMAKE FOR MACOSX
##---------------------------------------------------------------------------##

# CLEANUP
rm -rf CMakeCache.txt
rm -rf CMakeFiles

# SOURCE AND INSTALL
SOURCE=<SET_SOURCE_DIR>
INSTALL=<SET_INSTALL_DIR>

# BUILD OPTIONS
BUILD="DEBUG"
MPI="ON"

# TPL PATHS
HDF5_PATH="/vendors/parallel/hdf5"
MPI_PATH="/opt/mpi/current"

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
-DProfugus_CONFIGURE_OPTIONS_FILE:FILEPATH="${SOURCE}/install/base.cmake" \
-DProfugus_ASSERT_MISSING_PACKAGES:BOOL=OFF \
\
${SOURCE}

##---------------------------------------------------------------------------##
## end of cmake_macosx.sh
##---------------------------------------------------------------------------##
