##---------------------------------------------------------------------------##
## mesh/test/CMakeLists.txt
## Thomas M. Evans
## Wednesday July 11 14:35:42 2012
##---------------------------------------------------------------------------##
## Copyright (C) 2012 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for spn
##---------------------------------------------------------------------------##

INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(UtilsTest)

##---------------------------------------------------------------------------##
## TESTING
##---------------------------------------------------------------------------##

ADD_UTILS_TEST(tstDecomposition.cc      NP 1 2 4)
ADD_UTILS_TEST(tstStratimikosSolver.cc  NP 1   4)

##---------------------------------------------------------------------------##
## FILES
##---------------------------------------------------------------------------##

TRIBITS_COPY_FILES_TO_BINARY_DIR(XMLSOLVERINP
  SOURCE_FILES
  aztecoo.xml
  belos.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS tstStratimikosSolver
  )

##---------------------------------------------------------------------------##
##                    end of mesh/test/CMakeLists.txt
##---------------------------------------------------------------------------##
