#!/bin/sh
##---------------------------------------------------------------------------##
## CMAKE FOR Remus (DEBUG)
##---------------------------------------------------------------------------##

# CLEANUP
rm -rf CMakeCache.txt
rm -rf CMakeFiles

# SOURCE AND INSTALL
SOURCE=/projects/Profugus/master/source
INSTALL=/projects/Profugus/master/installs/release

# BUILD OPTIONS
BUILD="RELEASE"

# TPL PATHS
HDF5_PATH="/opt/hdf5-1.8.13-gcc"
MPI_PATH="/opt/openmpi-1.6.5-gcc/"
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
-D BLAS_LIBRARY_DIRS:PATH=/opt/ATLAS-3.10-gcc/lib/ \
-D BLAS_LIBRARY_NAMES:STRING="cblas;f77blas" \
-D LAPACK_LIBRARY_DIRS:PATH=/opt/ATLAS-3.10-gcc/lib/ \
-D LAPACK_LIBRARY_NAMES:STRING="lapack;atlas" \
\
-D CMAKE_C_COMPILER="${MPI_PATH}/bin/mpicc" \
-D CMAKE_CXX_COMPILER="${MPI_PATH}/bin/mpicxx" \
-D CMAKE_CXX_FLAGS="-Wno-unused-local-typedefs" \
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
