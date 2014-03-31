#!/bin/sh
##---------------------------------------------------------------------------##
## CMAKE FOR MACOSX
##---------------------------------------------------------------------------##

# >>>>> EDIT PATHS HERE <<<<<

SOURCE=<SET_SOURCE_DIR>
INSTALL=<SET_INSTALL_DIR>
BUILD="DEBUG"
MPI="ON"

HDF5_PATH="/vendors/parallel/hdf5"

##---------------------------------------------------------------------------##

cmake \
-DCMAKE_BUILD_TYPE:STRING="$BUILD" \
-DTPL_ENABLE_MPI:BOOL=$MPI \
-DCMAKE_INSTALL_PREFIX:PATH=$INSTALL \
-DProfugus_CONFIGURE_OPTIONS_FILE:FILEPATH="${SOURCE}/install/base.cmake" \
\
-DMPI_BASE_DIR:PATH=/opt/mpi/current \
\
-DTPL_ENABLE_HDF5:BOOL=ON \
-DHDF5_INCLUDE_DIRS:PATH=$HDF5_PATH/include \
-DHDF5_LIBRARY_DIRS:PATH=$HDF5_PATH/lib \
${SOURCE}

##---------------------------------------------------------------------------##
## end of cmake_macosx.sh
##---------------------------------------------------------------------------##
