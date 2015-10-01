##---------------------------------------------------------------------------##
## cxx11/CMakeLists.txt
## Thomas M. Evans
## Tuesday December 3 20:32:38 2013
##---------------------------------------------------------------------------##

INCLUDE(UtilsTest)

##---------------------------------------------------------------------------##
## TESTS
##---------------------------------------------------------------------------##

ADD_UTILS_TEST(tstTypes.cc             NP 1)
ADD_UTILS_TEST(tstPointers.cc          NP 1)
ADD_UTILS_TEST(tstLoops.cc             NP 1)
ADD_UTILS_TEST(tstClass.cc             NP 1)
ADD_UTILS_TEST(tstLambda.cc            NP 1)
ADD_UTILS_TEST(tstThread.cc            NP 1)
ADD_UTILS_TEST(tstMutex.cc             NP 1)
ADD_UTILS_TEST(tstConditionVariable.cc NP 1)

##---------------------------------------------------------------------------##
## end of cx11/CMakeLists.txt
##---------------------------------------------------------------------------##
