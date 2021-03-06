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

ADD_UTILS_TEST(tstContainer_Functions.cc NP 1)
ADD_UTILS_TEST(tstContainer_Props.cc     NP 1)
ADD_UTILS_TEST(tstHyperslab_Vector.cc    NP 1)
ADD_UTILS_TEST(tstHyperslab_View.cc      NP 1)
ADD_UTILS_TEST(tstMetaclass.cc           NP 1)
ADD_UTILS_TEST(tstPacking_Utils.cc       NP 1)
ADD_UTILS_TEST(tstRange.cc               NP 1)
ADD_UTILS_TEST(tstRegular_Indexer.cc     NP 1)
ADD_UTILS_TEST(tstSoft_Comparators.cc    NP 1)
ADD_UTILS_TEST(tstStatic_Map.cc          NP 1)
ADD_UTILS_TEST(tstVector_Functions.cc    NP 1)
ADD_UTILS_TEST(tstVector_Lite.cc         NP 1)
ADD_UTILS_TEST(tstView_Field.cc          NP 1)
ADD_UTILS_TEST(tstView_Field_Iterator.cc NP 1)

IF(TPL_ENABLE_HDF5)
  ADD_UTILS_TEST(tstSerial_HDF5_Writer.cc   NP 1 2)
  ADD_UTILS_TEST(tstParallel_HDF5_Writer.cc NP 1 2 4)
  ADD_UTILS_TEST(tstHDF5_Reader.cc          NP 1 2 4)
  ADD_UTILS_TEST(tstHDF5.cc                 NP 1 2 4)
ENDIF()

##---------------------------------------------------------------------------##
##                      end of utils/test/CMakeLists.txt
##---------------------------------------------------------------------------##
