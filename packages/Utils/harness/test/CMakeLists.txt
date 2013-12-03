##---------------------------------------------------------------------------##
## harness/test/CMakeLists.txt
## Thomas M. Evans
## Wednesday July 11 14:35:42 2012
##---------------------------------------------------------------------------##
## Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for harness
##---------------------------------------------------------------------------##

INCLUDE(NemesisTest)

##---------------------------------------------------------------------------##
## TESTS
##---------------------------------------------------------------------------##

ADD_NEMESIS_TEST(tstCasts.cc            NP 1 LEGACY)
ADD_NEMESIS_TEST(tstDBC.cc              NP 1 LEGACY)
ADD_NEMESIS_TEST(tstDemangle_Type.cc    NP 1 LEGACY)
ADD_NEMESIS_TEST(tstDiagnostics.cc      NP 1 LEGACY)
ADD_NEMESIS_TEST(tstScalar_Unit_Test.cc NP 1 LEGACY)
ADD_NEMESIS_TEST(tstSoft_Equiv.cc       NP 1 LEGACY)
ADD_NEMESIS_TEST(tstTypes.cc            NP 1 LEGACY)
ADD_NEMESIS_TEST(tstWarnings.cc         NP 1 LEGACY)

##---------------------------------------------------------------------------##
##                   end of harness/test/CMakeLists.txt
##---------------------------------------------------------------------------##
