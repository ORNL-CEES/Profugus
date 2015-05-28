##---------------------------------------------------------------------------##
## mc_solvers/test/CMakeLists.txt
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

# General tests
ADD_UTILS_TEST(tstChebyshevIteration)
ADD_UTILS_TEST(tstLinearSolverFactory)
ADD_UTILS_TEST(tstLinearSystemFactory)
ADD_UTILS_TEST(tstPolynomial)
ADD_UTILS_TEST(tstPolynomialBasis)
ADD_UTILS_TEST(tstPolynomialUtils)
ADD_UTILS_TEST(tstRichardsonIteration)

# Single process only tests
ADD_UTILS_TEST(tstMCSA             NP 1)
ADD_UTILS_TEST(tstMC_Data          NP 1)
ADD_UTILS_TEST(tstMC_Components    NP 1)
ADD_UTILS_TEST(tstMonteCarloSolver NP 1)

##---------------------------------------------------------------------------##
##                    end of mc_solvers/test/CMakeLists.txt
##---------------------------------------------------------------------------##
