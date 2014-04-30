##---------------------------------------------------------------------------##
## geometry/test/CMakeLists.txt
## Thomas M. Evans
## Thursday July 12 10:52:17 2012
##---------------------------------------------------------------------------##
## Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for geometry/rtk/test
##---------------------------------------------------------------------------##

INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE(TribitsAddAdvancedTest)
INCLUDE(UtilsTest)

##---------------------------------------------------------------------------##
## FILES
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## TESTING
##---------------------------------------------------------------------------##

ADD_UTILS_TEST(tstRTK_Functions.cc NP 1)
ADD_UTILS_TEST(tstRTK_Cell.cc      NP 1)
ADD_UTILS_TEST(tstRTK_Array.cc     NP 1)
ADD_UTILS_TEST(tstLattice.cc       NP 1)
ADD_UTILS_TEST(tstCore.cc          NP 1)
ADD_UTILS_TEST(tstVessel.cc        NP 1)

##---------------------------------------------------------------------------##
##                      end of geometry/test/CMakeLists.txt
##---------------------------------------------------------------------------##
