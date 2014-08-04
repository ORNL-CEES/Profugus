#!/bin/sh
##---------------------------------------------------------------------------##
## CMAKE FOR MACOSX WITH MCLS LIBRARY
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
-DProfugus_EXTRA_REPOSITORIES="ParaSails;MCLS" \
-DProfugus_ENABLE_MCLS:BOOL=ON \
-DProfugus_ENABLE_ML:BOOL=ON \
\
-DMCLS_ENABLE_DBC:BOOL=ON \
-DMCLS_ENABLE_TESTS:BOOL=OFF \
-DMCLS_ENABLE_EXAMPLES:BOOL=OFF \
-DMCLS_ENABLE_Epetra:BOOL=ON \
-DMCLS_ENABLE_EpetraExt:BOOL=ON \
-DMCLS_ENABLE_Ifpack:BOOL=ON \
-DMCLS_ENABLE_ParaSails:BOOL=ON \
-DMCLS_ENABLE_Thyra:BOOL=ON \
-DMCLS_ENABLE_Stratimikos:BOOL=ON \
\
${SOURCE}

##---------------------------------------------------------------------------##
## end of cmake_macosx_with_mcls.sh
##---------------------------------------------------------------------------##
