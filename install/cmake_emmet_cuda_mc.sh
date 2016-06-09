#!/bin/sh
##---------------------------------------------------------------------------##
## CMAKE FOR X86_64 WITH CUDA
##---------------------------------------------------------------------------##

# CLEANUP
rm -rf CMakeCache.txt
rm -rf CMakeFiles

# SOURCE AND INSTALL
SOURCE=<SOURCE_DIR>
INSTALL=<INSTALL_DIR>
MPI_HOME=<MPI_DIR>
HDF5_DIR=<HDF5_DIR>

# BUILD OPTIONS
BUILD="RELEASE"

##---------------------------------------------------------------------------##

cmake \
-DCMAKE_BUILD_TYPE:STRING="$BUILD" \
-DCMAKE_INSTALL_PREFIX:PATH=$INSTALL \
-DProfugus_ENABLE_CXX11:BOOL=ON \
-DCMAKE_C_COMPILER="${MPI_BIN}/mpicc" \
-DCMAKE_CXX_COMPILER="${MPI_BIN}/mpicxx" \
\
-DTPL_ENABLE_MPI:BOOL=ON \
-DMPI_BASE_DIR:PATH=${MPI_HOME} \
\
-DTPL_ENABLE_HDF5:BOOL=ON \
-DHDF5_INCLUDE_DIRS:PATH=${HDF5_DIR}/include \
-DHDF5_LIBRARY_DIRS:PATH=${HDF5_DIR}/lib \
\
-DProfugus_ASSERT_MISSING_PACKAGES:BOOL=OFF \
\
-DTPL_ENABLE_CUDA:BOOL=ON \
-DProfugus_ENABLE_CUDA:BOOL=ON \
-DCUDA_PROPAGATE_HOST_FLAGS:BOOL=OFF \
-DCUDA_NVCC_FLAGS:STRING="--std=c++11;-arch=sm_35;-lineinfo" \
\
-DProfugus_ENABLE_TESTS:BOOL=ON \
-DProfugus_ENABLE_MC:BOOL=ON \
-DProfugus_ENABLE_CUDA:BOOL=ON \
-DProfugus_ENABLE_CudaUtils:BOOL=ON \
\
-DKokkos_ENABLE_Cuda:BOOL=OFF \
\
${SOURCE}

##---------------------------------------------------------------------------##
## end of cmake_x86_64.sh
##---------------------------------------------------------------------------##
