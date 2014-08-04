#!/bin/sh
##---------------------------------------------------------------------------##
## CMAKE FOR xk7
##---------------------------------------------------------------------------##

# CLEANUP
rm -rf CMakeCache.txt
rm -rf CMakeFiles

# SOURCE AND INSTALL
SOURCE=<SET_SOURCE_DIR>
INSTALL=<SET_INSTALL_DIR>

# BUILD OPTIONS
BUILD="RELEASE"

# TPL PATHS
COMP_PATH="/opt/cray/xt-asyncpe/5.24/bin"
PYTHON_EXE="/opt/sw/xk6/python/2.7.3/sles11.1_gnu4.3.4/bin/python"

##---------------------------------------------------------------------------##

cmake \
-DCMAKE_BUILD_TYPE:STRING="$BUILD" \
-DCMAKE_INSTALL_PREFIX:PATH=$INSTALL \
\
-DPYTHON_EXECUTABLE:FILEPATH=$PYTHON_EXE \
-DTPL_ENABLE_MPI:BOOL=ON \
-DMPI_C_COMPILER:STRING=${COMP_PATH}/cc
-DMPI_CXX_COMPILER:STRING=${COMP_PATH}/CC
-DMPI_Fortran_COMPILER:STRING=${COMP_PATH}/ftn
-DProgugus_ENABLE_TESTS:BOOL=OFF \
-DBUILD_SHARED_LIBS:BOOL=OFF \
-DBUILD_STATIC_LIBS:BOOL=ON
\
-DMPI_BASE_DIR:PATH=$MPI_PATH \
\
-DTPL_ENABLE_HDF5:BOOL=ON \
-DHDF5_INCLUDE_DIRS:PATH=$HDF5_DIR/include \
-DHDF5_LIBRARY_DIRS:PATH=$HDF5_DIR/lib \
\
-DTPL_BLAS_LIBRARIES:STRING="${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu.a" \
-DTPL_LAPACK_LIBRARIES:STRING="${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu.a" \
\
-DProfugus_CONFIGURE_OPTIONS_FILE:FILEPATH="${SOURCE}/install/base.cmake" \
-DProfugus_ASSERT_MISSING_PACKAGES:BOOL=OFF \
\
${SOURCE}

##---------------------------------------------------------------------------##
## end of cmake_xk7.sh
##---------------------------------------------------------------------------##
