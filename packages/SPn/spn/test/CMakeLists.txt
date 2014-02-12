##---------------------------------------------------------------------------##
## spn/test/CMakeLists.txt
## Thomas M. Evans
## Wednesday July 11 14:35:42 2012
##---------------------------------------------------------------------------##
## Copyright (C) 2012 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for spn
##---------------------------------------------------------------------------##

INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(UtilsTest)

##---------------------------------------------------------------------------##
## FILES
##---------------------------------------------------------------------------##

# Make testing library
TRIBITS_ADD_LIBRARY(spn_test_lib
  SOURCES Test_XS.cc
  TESTONLY
  DEPLIBS SPn)

##---------------------------------------------------------------------------##
## TESTING
##---------------------------------------------------------------------------##

ADD_UTILS_TEST(tstDimensions.cc          NP 1)
ADD_UTILS_TEST(tstMoment_Coefficients.cc NP 1     DEPLIBS spn_test_lib)
ADD_UTILS_TEST(tstSDM_Face_Field.cc      NP 1 2 4 DEPLIBS spn_test_lib)
ADD_UTILS_TEST(tstMesh.cc                NP 1 2 4)

##---------------------------------------------------------------------------##
##                    end of spn/test/CMakeLists.txt
##---------------------------------------------------------------------------##
