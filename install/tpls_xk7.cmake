##---------------------------------------------------------------------------##
## CMAKE REGRESSION TPLs FILE
##---------------------------------------------------------------------------##

# TPLS
SET(TPL_ENABLE_SILO ON CACHE BOOL "")
SET(SILO_INCLUDE_DIRS ${PROJ_HOME}/tpls/silo/include CACHE PATH "")
SET(SILO_LIBRARY_DIRS ${PROJ_HOME}/tpls/silo/lib     CACHE PATH "")

SET(TPL_ENABLE_HDF5 ON CACHE BOOL "")
SET(HDF5_INCLUDE_DIRS $ENV{HDF5_DIR}/include CACHE PATH "")
SET(HDF5_LIBRARY_DIRS $ENV{HDF5_DIR}/lib     CACHE PATH "")

SET(TPL_ENABLE_LAVA ON CACHE BOOL "")
SET(LAVA_INCLUDE_DIRS ${PROJ_HOME}/tpls/lava/include CACHE PATH "")
SET(LAVA_LIBRARY_DIRS ${PROJ_HOME}/tpls/lava/lib     CACHE PATH "")
SET(MCNP_EXECUTABLE ${PROJ_HOME}/tpls/mcnp/Source/src/mcnp5
  CACHE FILEPATH "")
SET(MCNP_DATAPATH ${PROJ_PATH}/mcnp CACHE PATH "")

SET(TPL_BLAS_LIBRARIES "$ENV{CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu.a"
  CACHE STRING "")
SET(TPL_LAPACK_LIBRARIES "$ENV{CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu.a"
  CACHE STRING "")
