#!/bin/sh
##---------------------------------------------------------------------------##
## CMAKE FOR Remus (DEBUG)
##---------------------------------------------------------------------------##

# CLEANUP
rm -rf CMakeCache.txt
rm -rf CMakeFiles

module load linux-x86_64/gcc@5.3.0%gcc@4.8.5+gold-5hy3c4b
module load mpi/openmpi/gcc5/1.10.2

# SOURCE AND INSTALL
SOURCE=/home/gqe/Codes/Profugus/source/
INSTALL=/home/gqe/Codes/Profugus/installs/release

# BUILD OPTIONS
BUILD="RELEASE"

# TPL PATHS
HDF5_PATH="/software/tools/apps/hdf5/gcc5/1.8.16/"
MPI_PATH="/software/tools/apps/openmpi/gcc5/1.10.2/"
HDF5_LIB_PATH="${HDF5_PATH}/lib/"
HDF5_INC_PATH="${HDF5_PATH}/include/"

##---------------------------------------------------------------------------##

cmake \
-D CMAKE_BUILD_TYPE:STRING="$BUILD" \
-D CMAKE_INSTALL_PREFIX:PATH=$INSTALL \
\
-D TPL_ENABLE_MPI:BOOL=ON \
-D MPI_BASE_DIR:PATH=${MPI_PATH} \
\
-D BLAS_LIBRARY_DIRS:PATH=/software/tools/apps/openblas/gcc5/nothreads/lib/ \
-D BLAS_LIBRARY_NAMES:STRING="openblas" \
-D LAPACK_LIBRARY_DIRS:PATH=/software/tools/apps/openblas/gcc5/nothreads/lib/ \
-D LAPACK_LIBRARY_NAMES:STRING="openblas" \
\
-D CMAKE_C_COMPILER="${MPI_PATH}/bin/mpicc" \
-D CMAKE_CXX_COMPILER="${MPI_PATH}/bin/mpicxx" \
-D CMAKE_C_FLAGS="-march=native -lgfortran" \
-D CMAKE_CXX_FLAGS="-march=native -Wno-unused-local-typedefs -lgfortran" \
-D CMAKE_EXE_LINKER_FLAGS="${CMAKE_EXE_LINKER_FLAGS} -lgfortran" \
\
-D TPL_ENABLE_HDF5:BOOL=ON \
-D HDF5_INCLUDE_DIRS:PATH=$HDF5_INC_PATH \
-D HDF5_LIBRARY_DIRS:PATH=$HDF5_LIB_PATH \
\
-D Profugus_CONFIGURE_OPTIONS_FILE:FILEPATH="${SOURCE}/install/base.cmake" \
-D Profugus_ASSERT_MISSING_PACKAGES:BOOL=OFF \
\
-D Profugus_ENABLE_CXX11:BOOL=ON \
-D Profugus_ENABLE_SPn:BOOL=ON \
-D Profugus_ENABLE_Alea:BOOL=ON \
-D Profugus_ENABLE_MC:BOOL=ON \
-D Profugus_ENABLE_Utils:BOOL=ON \
\
-D Kokkos_ENABLE_Pthread:BOOL=ON \
-D Tpetra_INST_PTHREAD:BOOL=ON \
-D KokkosClassic_DefaultNode:STRING="Kokkos::Compat::KokkosThreadsWrapperNode" \
\
${SOURCE}

##---------------------------------------------------------------------------##
## end of cmake for REMUS
##---------------------------------------------------------------------------##
