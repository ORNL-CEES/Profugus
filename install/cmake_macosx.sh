#!/bin/sh
##---------------------------------------------------------------------------##
## CMAKE FOR MACOSX
##---------------------------------------------------------------------------##

# >>>>> EDIT PATHS HERE <<<<<

SOURCE=<SET_SOURCE_DIR>
INSTALL=<SET_INSTALL_DIR>
BUILD="DEBUG"
MPI="ON"

##---------------------------------------------------------------------------##
## No need to edit below here

cmake \
-DCMAKE_BUILD_TYPE:STRING="$BUILD" \
-DTPL_ENABLE_MPI:BOOL=$MPI \
-DCMAKE_INSTALL_PREFIX=$INSTALL \
-DProfugus_CONFIGURE_OPTIONS_FILE:FILEPATH="${SOURCE}/install/base.cmake;${SOURCE}/install/tpls_macosx.cmake" \
${SOURCE}

##---------------------------------------------------------------------------##
## end of cmake_macosx.sh
##---------------------------------------------------------------------------##
