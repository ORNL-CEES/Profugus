#!/bin/sh
##---------------------------------------------------------------------------##
## CMAKE FOR emmet
##---------------------------------------------------------------------------##

# CLEANUP
rm -rf CMakeCache.txt
rm -rf CMakeFiles

# SOURCE AND INSTALL
SOURCE=/home/9te/profugus/Profugus
INSTALL=/home/9te/profugus/_openmp/debug

# BUILD OPTIONS
BUILD="DEBUG"

HDF5_PATH=/usr/lib64
# TPL PATHS
MPI_PATH="/usr/lib64/openmpi"
HDF5_LIB_PATH="${MPI_PATH}/lib"
HDF5_INC_PATH="/usr/include/openmpi-x86_64"

##---------------------------------------------------------------------------##

cmake \
-DCMAKE_BUILD_TYPE:STRING="$BUILD" \
-DCMAKE_INSTALL_PREFIX:PATH=$INSTALL \
\
-DTPL_ENABLE_MPI:BOOL=ON \
-DMPI_BASE_DIR:PATH=${MPI_PATH} \
\
-DCMAKE_C_COMPILER="${MPI_PATH}/bin/mpicc" \
-DCMAKE_CXX_COMPILER="${MPI_PATH}/bin/mpicxx" \
-DCMAKE_CXX_FLAGS="-Wno-unused-local-typedefs" \
\
-DTPL_ENABLE_HDF5:BOOL=ON \
-DHDF5_INCLUDE_DIRS:PATH=$HDF5_INC_PATH \
-DHDF5_LIBRARY_DIRS:PATH=$HDF5_LIB_PATH \
\
-DProfugus_CONFIGURE_OPTIONS_FILE:FILEPATH="${SOURCE}/install/base.cmake" \
-DProfugus_ASSERT_MISSING_PACKAGES:BOOL=OFF \
\
-DProfugus_ENABLE_CXX11:BOOL=ON \
-DProfugus_ENABLE_SPn:BOOL=ON \
-DProfugus_ENABLE_Alea:BOOL=OFF \
-DProfugus_ENABLE_MC:BOOL=ON \
-DProfugus_ENABLE_Utils:BOOL=ON \
\
-DProfugus_ENABLE_OpenMP:BOOL=ON \
-DUTILS_DIAGNOSTICS:STRING=0 \
\
-DProfugus_ENABLE_MueLu:BOOL=OFF \
-DKokkos_ENABLE_Pthread:BOOL=ON \
-DTpetra_INST_PTHREAD:BOOL=ON \
-DKokkosClassic_DefaultNode:STRING="Kokkos::Compat::KokkosThreadsWrapperNode" \
\
${SOURCE}

##---------------------------------------------------------------------------##
## end of cmake_emmet_cpu.sh
##---------------------------------------------------------------------------##
