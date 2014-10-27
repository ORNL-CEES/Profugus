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
BUILD="RELWITHDEBINFO"

# TPL PATHS
COMP_PATH="/opt/cray/xt-asyncpe/5.27/bin"
PYTHON_EXE="/opt/sw/xk6/python/2.7.3/sles11.1_gnu4.3.4/bin/python"

##---------------------------------------------------------------------------##

cmake \
-DCMAKE_BUILD_TYPE:STRING="$BUILD" \
-DCMAKE_INSTALL_PREFIX:PATH=$INSTALL \
\
-DPYTHON_EXECUTABLE:FILEPATH=$PYTHON_EXE \
-DTPL_ENABLE_MPI:BOOL=OFF \
-DCMAKE_C_COMPILER:STRING=pgcc \
-DCMAKE_CXX_COMPILER:STRING=pgc++ \
-DCMAKE_Fortran_COMPILER:STRING=pgf90 \
-DProfugus_ENABLE_TESTS:BOOL=ON \
-DBUILD_SHARED_LIBS:BOOL=OFF \
-DBUILD_STATIC_LIBS:BOOL=ON \
\
-DProfugus_ENABLE_INSTALL_CMAKE_CONFIG_FILES:BOOL=ON \
\
-DProfugus_ENABLE_SPn:BOOL=OFF \
-DProfugus_ENABLE_MueLu:BOOL=OFF \
\
-DCMAKE_CXX_FLAGS:STRING="-tp=bulldozer --c++11 -V14.7 -w -D_X86INTRIN_H_INCLUDED" \
-DCMAKE_C_FLAGS:STRING="-tp=bulldozer -V14.7 -w" \
-DCMAKE_Fortran_FLAGS:STRING="-tp=bulldozer -V14.7 -w" \
-DCMAKE_EXE_LINKER_FLAGS:STRING="-L/opt/gcc/4.8.2/snos/lib64 -latomic" \
\
-DMPI_EXEC:FILEPATH="/sw/xk6/site-aprun/bin/aprun" \
-DMPIEXEC_NUMPROC_FLAG:STRING="-n" \
\
-DTPL_ENABLE_HDF5:BOOL=ON \
-DHDF5_INCLUDE_DIRS:PATH=$HDF5_DIR/include \
-DHDF5_LIBRARY_DIRS:PATH=$HDF5_DIR/lib \
\
-DTPL_BLAS_LIBRARIES:STRING="${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_pgi.a" \
-DTPL_LAPACK_LIBRARIES:STRING="${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_pgi.a" \
\
-DProfugus_CONFIGURE_OPTIONS_FILE:FILEPATH="${SOURCE}/install/base.cmake" \
-DProfugus_ASSERT_MISSING_PACKAGES:BOOL=OFF \
\
${SOURCE}

##---------------------------------------------------------------------------##
## end of cmake_xk7.sh
##---------------------------------------------------------------------------##
