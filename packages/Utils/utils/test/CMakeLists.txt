##---------------------------------------------------------------------------##
## utils/test/CMakeLists.txt
## Thomas M. Evans
## Thursday July 12 10:52:17 2012
##---------------------------------------------------------------------------##
## Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for utils//test
##---------------------------------------------------------------------------##

INCLUDE(UtilsTest)

##---------------------------------------------------------------------------##
## TESTS
##---------------------------------------------------------------------------##

ADD_UTILS_TEST(tstStatic_Map.cc          NP 1)
ADD_UTILS_TEST(tstVector_Lite.cc         NP 1)
ADD_UTILS_TEST(tstContainer_Functions.cc NP 1)
ADD_UTILS_TEST(tstContainer_Props.cc     NP 1)
ADD_UTILS_TEST(tstSoft_Comparators.cc    NP 1)
ADD_UTILS_TEST(tstPacking_Utils.cc       NP 1)
ADD_UTILS_TEST(tstVector_Functions.cc    NP 1)

IF(TPL_ENABLE_HDF5)
  ADD_UTILS_TEST(tstSerial_HDF5_Writer.cc   NP 1 2)
  ADD_UTILS_TEST(tstParallel_HDF5_Writer.cc NP 1 2 4)
  ADD_UTILS_TEST(tstHDF5_Reader.cc          NP 1 2 4)
  ADD_UTILS_TEST(tstHDF5.cc                 NP 1 2 4)
ENDIF()

##---------------------------------------------------------------------------##
##                      end of utils/test/CMakeLists.txt
##---------------------------------------------------------------------------##
