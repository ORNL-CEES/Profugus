##---------------------------------------------------------------------------##
## harness/test/CMakeLists.txt
## Thomas M. Evans
## Tuesday December 3 20:32:38 2013
##---------------------------------------------------------------------------##

INCLUDE(UtilsTest)

##---------------------------------------------------------------------------##
## TESTS
##---------------------------------------------------------------------------##

ADD_UTILS_TEST(tstDBC.cc NP 1)
ADD_UTILS_TEST(tstDiagnostics.cc NP 1)
#ADD_UTILS_TEST(tstTypes.cc NP 1)
#ADD_UTILS_TEST(tstSoft_Equiv.cc NP 1)
#ADD_UTILS_TEST(tstWarnings.cc NP 1)

##---------------------------------------------------------------------------##
## end of harness/test/CMakeLists.txt
##---------------------------------------------------------------------------##
