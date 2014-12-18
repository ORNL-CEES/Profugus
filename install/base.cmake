##---------------------------------------------------------------------------##
## CMAKE BASE FILE
##---------------------------------------------------------------------------##

# Default build all packages
SET(Profugus_ENABLE_Utils   ON CACHE BOOL "")
SET(Profugus_ENABLE_Matprop ON CACHE BOOL "")
SET(Profugus_ENABLE_SPn     ON CACHE BOOL "")
SET(Profugus_ENABLE_MC      ON CACHE BOOL "")

# Turn on tests
SET(Profugus_ENABLE_TESTS ON CACHE BOOL "")
SET(Profugus_TEST_CATEGORIES "BASIC" CACHE STRING "")

# Turn on SS code and optional packages by default
SET(Profugus_ENABLE_ALL_FORWARD_DEP_PACKAGES OFF CACHE BOOL "")
SET(Profugus_ENABLE_ALL_OPTIONAL_PACKAGES    ON  CACHE BOOL "")
SET(Profugus_ENABLE_SECONDARY_STABLE_CODE    ON  CACHE BOOL "")

# Set explicit instantiation options
SET(Profugus_ENABLE_EXPLICIT_INSTANTIATION ON CACHE BOOL "")
SET(Teuchos_ENABLE_FLOAT OFF  CACHE BOOL "")
SET(Teuchos_ENABLE_COMPLEX OFF CACHE BOOL "")
SET(Tpetra_INST_FLOAT OFF CACHE BOOL "")
SET(Tpetra_INST_COMPLEX_FLOAT OFF CACHE BOOL "")
SET(Tpetra_INST_COMPLEX_DOUBLE OFF CACHE BOOL "")
SET(Thyra_ENABLE_EXPLICIT_INSTANTIATION OFF CACHE BOOL "")
SET(Stratimikos_ENABLE_EXPLICIT_INSTANTIATION OFF CACHE BOOL "")

# Up the max num procs
SET(MPI_EXEC_MAX_NUMPROCS 8 CACHE STRING "")

# Turn off binutils
SET(Teuchos_ENABLE_BinUtils OFF CACHE BOOL "")

# Turn off Zoltan2
SET(Profugus_ENABLE_Zoltan2 OFF CACHE BOOL "")

# Compiler options
SET(BUILD_SHARED_LIBS ON CACHE BOOL "")
SET(CMAKE_CXX_FLAGS "-std=c++11 -Wno-deprecated-declarations" CACHE STRING "")
SET(Profugus_ENABLE_CXX11:BOOL=ON)

# TriBITS stuff
SET(Profugus_ENABLE_INSTALL_CMAKE_CONFIG_FILES OFF CACHE BOOL "")
SET(Profugus_DEPS_XML_OUTPUT_FILE "" CACHE FILEPATH "")
