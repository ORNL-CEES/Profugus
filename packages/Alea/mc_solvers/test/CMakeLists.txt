##---------------------------------------------------------------------------##
## alea/test/CMakeLists.txt
## Thomas M. Evans
## Wednesday July 11 14:35:42 2012
##---------------------------------------------------------------------------##
## Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for spn
##---------------------------------------------------------------------------##

INCLUDE(UtilsTest)

##---------------------------------------------------------------------------##
## TESTING
##---------------------------------------------------------------------------##

ADD_UTILS_TEST(tstChebyshevIteration)
ADD_UTILS_TEST(tstLinearSolverFactory)
ADD_UTILS_TEST(tstLinearSystemFactory)
ADD_UTILS_TEST(tstPolynomial)
ADD_UTILS_TEST(tstPolynomialBasis)
ADD_UTILS_TEST(tstPolynomialUtils)
ADD_UTILS_TEST(tstRichardsonIteration)

##---------------------------------------------------------------------------##
##                    end of alea/test/CMakeLists.txt
##---------------------------------------------------------------------------##
