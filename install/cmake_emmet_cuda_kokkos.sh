#!/bin/sh
##---------------------------------------------------------------------------##
## CMAKE FOR X86_64
##---------------------------------------------------------------------------##

# CLEANUP
rm -rf CMakeCache.txt
rm -rf CMakeFiles

# SOURCE AND INSTALL
SOURCE=<SOURCE_DIR>
INSTALL=<INSTALL_DIR>

# BUILD OPTIONS
BUILD="RELWITHDEBINFO"

HDF5_PATH=/usr/lib64
# TPL PATHS
MPI_PATH="/usr/lib64/openmpi"
HDF5_PATH="${MPI_PATH}"

# Cuda flags
CUDA_NVCC_FLAGS="-std=c++11 -Xcudafe --diag_suppress=code_is_unreachable"

##---------------------------------------------------------------------------##

cmake \
-DCMAKE_BUILD_TYPE:STRING="$BUILD" \
-DCMAKE_INSTALL_PREFIX:PATH=$INSTALL \
\
-DTPL_ENABLE_MPI:BOOL=ON \
-DMPI_BASE_DIR:PATH=${MPI_PATH} \
\
-DCMAKE_C_COMPILER="${MPI_PATH}/bin/mpicc" \
-DCMAKE_CXX_COMPILER="${PWD}/nvcc_wrapper" \
-DCMAKE_CXX_FLAGS="-Wno-unused-local-typedefs -lineinfo -Xcudafe --diag_suppress=code_is_unreachable" \
\
-DTPL_ENABLE_HDF5:BOOL=ON \
-DHDF5_INCLUDE_DIRS:PATH=${HDF5_PATH}/include \
-DHDF5_LIBRARY_DIRS:PATH=${HDF5_PATH}/lib \
\
-DProfugus_CONFIGURE_OPTIONS_FILE:FILEPATH="${SOURCE}/install/base.cmake" \
-DProfugus_ASSERT_MISSING_PACKAGES:BOOL=OFF \
\
-DTPL_ENABLE_CUDA:BOOL=ON \
-DCUDA_NVCC_FLAGS:STRING=${CUDA_NVCC_FLAGS} \
\
-DProfugus_ENABLE_CXX11:BOOL=ON \
-DProfugus_ENABLE_SPn:BOOL=ON \
-DProfugus_ENABLE_Alea:BOOL=ON \
-DProfugus_ENABLE_MC:BOOL=OFF \
-DProfugus_ENABLE_Utils:BOOL=ON \
-DProfugus_ENABLE_OpenMP:BOOL=ON \
-DProfugus_ENABLE_Teko:BOOL=OFF \
-DProfugus_ENABLE_Kokkos:BOOL=ON \
-DProfugus_ENABLE_Tpetra:BOOL=ON \
-DTpetra_ENABLE_TESTS:BOOL=OFF \
\
-DTPL_ENABLE_Pthread:BOOL=ON \
-DKokkos_ENABLE_Pthread:BOOL=ON \
-DKokkos_ENABLE_OpenMP:BOOL=ON \
-DKokkosClassic_DefaultNode:STRING="Kokkos::Compat::KokkosCudaWrapperNode" \
-DKokkos_ENABLE_CXX11:BOOL=ON \
-DKokkos_ENABLE_Cuda_UVM:BOOL=ON \
\
-DTpetra_INST_OPENMP:BOOL=ON \
-DTpetra_INST_PTHREAD:BOOL=ON \
\
${SOURCE}

##---------------------------------------------------------------------------##
## end of cmake_emmet_gpu.sh
##---------------------------------------------------------------------------##
